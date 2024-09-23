use basis::{BiVarPolyDistortions, DistortionBasis};
use pyo3::prelude::*;
use core::f64;
use std::fmt::Debug;
use glob;
use rayon::prelude::*;
use std::ops::{Add, AddAssign, Mul, Sub};
use rustfft::num_complex::ComplexFloat;

mod errors;
pub mod basis;
pub mod io;
pub use crate::io::{Image, Coordinate};
pub use crate::errors::{MavDACError, Result};

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

/// A Python module implemented in Rust.
#[pymodule]
fn mavdac(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(cli_main, m)?)?;
    m.add_function(wrap_pyfunction!(get_coordinates, m)?)?;
    m.add_function(wrap_pyfunction!(load_images, m)?)?;
    m.add_function(wrap_pyfunction!(measure_cogs, m)?)?;
    m.add_function(wrap_pyfunction!(cli_diff, m)?)?;
    Ok(())
}

// intended usage is to have the core functions exposed as much as possible,
// and to provide a CLI with the same API as the existing rust one.


#[pyfunction]
#[pyo3(signature = (pattern, coordinates=None))]
fn cli_main(pattern: &str, coordinates: Option<&str>) -> PyResult<()> {
    const FLUXTHRESH: f64 = 10000.0;
    
    // if coordinates exist, ingest them first since it is the quickest thing to do
    // and a likely source of user error
    let coords: Option<Vec<Coordinate>> = match coordinates {
        Some(filename) => Some(get_coordinates(filename)?),
        None => None,
    };

    println!("loading images");
    let images = load_images(pattern)?;
    let shape = images[0].shape.clone();

    let grid = Grid::Hex {
        pitch: 1.0/7.5e-3,
        rotation: 0.0,
        offset: Vec2D{x:0.0, y:0.0},
    };

    // compute all centroids, indexed by pinhole id and exposure id,
    eprintln!("computing centroids");
    const RAD: usize = 10;
    let cogs = measure_cogs(images, grid, RAD, FLUXTHRESH);

    eprintln!("solving for distortions");
    let distortions = solve_distortions(cogs, 3, shape)?;

    match coords {
        Some(mut coords) => {
            for coord in &mut coords {
                coord.dist = Some(distortions.eval(&coord.pos));
                println!("{:?}", coord);
            }
        },
        None => {
            for coeff in distortions.get_coeffs() {
                println!("{:13.8}, {:13.8}", coeff.x, coeff.y);
            }
        },
    }
    Ok(())
}

#[pyfunction]
#[pyo3(signature = (pattern, coordinates=None))]
fn cli_diff(pattern: &str, coordinates: Option<&str>) -> PyResult<()> {
    const FLUXTHRESH: f64 = 10000.0;
    
    // if coordinates exist, ingest them first since it is the quickest thing to do
    // and a likely source of user error
    let coords: Option<Vec<Coordinate>> = match coordinates {
        Some(filename) => Some(get_coordinates(filename)?),
        None => None,
    };

    println!("loading images");
    let images = load_images(pattern)?;
    let shape = images[0].shape.clone();
        
    let grid = Grid::Hex {
        pitch: 1.0/7.5e-3,
        rotation: 0.0,
        offset: Vec2D{x:0.0, y:0.0},
    };

    // do centroids for base set of pinholes over all exposures, then filter out
    // any that are lower than FLUXTHRESH in any exposure
    const RAD: usize = 50;
    let pinholes = grid.all_points(images[0].shape[1], images[0].shape[0]);
    // cogs should be a n_pinholes x n_images 
    let cogs: Vec<Option<Vec<Centroid>>> = pinholes.iter().map(|pinhole|
        images.iter().map(|image|
            image.cog(&(pinhole.clone()+image.shift), RAD)
        ).collect::<Vec<Centroid>>()
    ).map(|pinhole_cogs|
        if pinhole_cogs.iter().all(|cog|
            cog.flux > FLUXTHRESH
        ) {
            Some(pinhole_cogs)
        } else {
            None
        }
    ).collect();

    // calculate mean position of each pinhole over all exposures, and use this
    // to calculate residual distortion after pinhole position is corrected for.
    // this should give trivial distortion estimations if there is only one
    // exposure.
    let mut cleaned_cogs: Vec<Centroid> = vec![];
    for (i, _pinhole) in pinholes.iter().enumerate() {
        if let Some(pinhole_cogs) = &cogs[i] {
            let mut mean_pinhole_position = Vec2D {x: 0.0, y: 0.0};
            for (j, image) in images.iter().enumerate() {
                mean_pinhole_position += pinhole_cogs[j].cog - image.shift;
            }
            mean_pinhole_position.x /= images.len() as f64;
            mean_pinhole_position.y /= images.len() as f64;
            for (j, image) in images.iter().enumerate() {
                let mut cleaned_cog = pinhole_cogs[j].clone();
                cleaned_cog.pos = mean_pinhole_position.clone() + image.shift;
                cleaned_cogs.push(cleaned_cog);
            }
        }
    }
 
    eprintln!("solving for distortions");
    let distortions = solve_distortions(cleaned_cogs, 3, shape)?;

    match coords {
        Some(mut coords) => {
            for coord in &mut coords {
                coord.dist = Some(distortions.eval(&coord.pos));
                println!("{:?}", coord);
            }
        },
        None => {
            for coeff in distortions.get_coeffs() {
                println!("{:13.8}, {:13.8}", coeff.x, coeff.y);
            }
        },
    }
    
    Ok(())
}

#[pyfunction]
fn get_coordinates(filename: &str) -> PyResult<Vec<Coordinate>> {
    let contents = std::fs::read_to_string(filename)?;
    let mut coords: Vec<Coordinate> = vec![];
    for line in contents.split('\n') {
        if line.len() > 0 {
            coords.push(Coordinate::from_str(line)?);
        }
    }
    Ok(coords)
}

#[pyfunction]
fn load_images(pattern: &str) -> Result<Vec<Image>> {
    glob::glob(&pattern)?.into_iter()
    .map(|path| Image::from_fits(path?))
    .collect::<Result<Vec<Image>>>()
}

#[pyfunction]
fn measure_cogs(images: Vec<Image>, grid: Grid, rad: usize, flux_thresh: f64) -> Vec<Centroid> {
    images.par_iter().flat_map(|image| {
        let grid_shifted = grid + image.shift;
        image.cogs(&grid_shifted, rad).into_iter().filter_map(|cog|
            {
                //println!("{:?}", cog.flux);
                if cog.flux > flux_thresh {
                    Some(cog)
                } else {
                    None
                }
            }
        ).collect::<Vec<Centroid>>()
    }).collect()
}

#[pyfunction]
fn solve_distortions(cogs: Vec<Centroid>, degree: usize, shape: [usize;2]) -> PyResult<BiVarPolyDistortions> {
    let mut distortions = BiVarPolyDistortions::new(degree, shape);
    distortions.solve_for_coeffs(cogs)?;
    Ok(distortions)
}


#[derive(Clone,Debug,Copy)]
#[pyclass]
pub struct Vec2D {
    pub x: f64,
    pub y: f64,
}

impl AddAssign for Vec2D {
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
    }
}
impl Add for Vec2D {
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
    type Output = Self;
}
impl Sub for Vec2D {
    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
    type Output = Self;
}

impl Mul<f64> for Vec2D {
    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
    type Output = Self;
}


#[pyclass]
#[derive(Clone,Copy,Debug)]
pub enum Grid {
    Hex {
        pitch: f64,  // pixels
        rotation: f64,  // radians
        offset: Vec2D,  // pixels
    },
}

impl AddAssign<Vec2D> for Grid {
    fn add_assign(&mut self, rhs: Vec2D) {
        match self {
            Grid::Hex { offset , ..} =>  {
                *offset += rhs;
            },
        }
    }
}
impl Add<Vec2D> for Grid {
    fn add(self, rhs: Vec2D) -> Self::Output {
        match self {
            Grid::Hex { pitch, rotation, offset } => 
                Grid::Hex { pitch, rotation, offset: offset + rhs },
        }   
    }
    type Output = Self;
}

impl Grid {
    pub fn all_points(&self, width: usize, height: usize) -> Vec<Vec2D> {
        let max_rad = width.max(height)*2;
        match self {
            Grid::Hex { pitch, rotation, offset } => {
                // first make a square grid with way too many points
                (0..2*max_rad).map(|x| x as f64)
                .flat_map(|x| (0..2*max_rad)
                    // shift it to be centered at origin
                    .map(|y| y as f64).map(move |y|
                        (x - max_rad as f64,y - max_rad as f64)
                    )
                )
                // scale it to pixel units
                .map(|(x,y)| (x*pitch, y*pitch))
                // then map it to a hex grid with gaps
                .map(|(x,y)| (x+0.5*y,y*(3.0).sqrt()/2.0))
                // rotate
                .map(|(x,y)| (
                    x*rotation.cos()-y*rotation.sin(),
                    x*rotation.sin()+y*rotation.cos(),
                ))
                // now apply the offset:
                .map(|(x,y)| (x+offset.x, y+offset.y))
                // shift back to valid pixel range
                .map(|(x,y)| (x + (width/2) as f64 - 0.5, y + (height/2) as f64 - 0.5))
                .filter(|(x,y)| *x >= 0.0 && *x < width as f64 && *y >= 0.0 && *y < height as f64)
                .map(|(x,y)| Vec2D{x,y})
                .collect()
            },
        }
    }
}

#[derive(Debug, Clone)]
#[pyclass]
pub struct Centroid {
    pub cog: Vec2D,
    pub flux: f64,
    pub pos: Vec2D,
}
