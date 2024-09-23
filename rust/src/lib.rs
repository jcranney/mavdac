use basis::{BiVarPolyDistortions, DistortionBasis};
use pyo3::prelude::*;
use core::f64;
use std::fmt::Debug;
use glob;
use rayon::prelude::*;
use std::ops::{Add, AddAssign, Mul};
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
    Ok(())
}

// intended usage is to have the core functions exposed as much as possible,
// and to provide a CLI with the same API as the existing rust one.




#[pyfunction]
#[pyo3(signature = (pattern, coordinates=None))]
fn cli_main(pattern: &str, coordinates: Option<&str>) -> Result<()> {
    const FLUXTHRESH: f64 = 10000.0;
    
    // if coordinates exist, ingest them first since it is the quickest thing to do
    // and a likely source of user error
    let coords: Option<Vec<Coordinate>> = match coordinates {
        Some(filename) => {
            let contents = std::fs::read_to_string(filename)?;
            let mut coords: Vec<Coordinate> = vec![];
            for line in contents.split('\n') {
                if line.len() > 0 {
                    coords.push(Coordinate::from_str(line)?);
                }
            }
            Some(coords)
        },
        None => None,
    };

    println!("loading images");
    let images = glob::glob(&pattern)?.into_iter()
    .map(|path| Image::from_fits(path?))
    .collect::<Result<Vec<Image>>>()?;
        
    // define pinhole grid - needs to be changed to be user defined or fit from
    // data.
    /*
    let grid = Grid::Hex {
        pitch: 159.0,
        rotation: PI/2.0,
        offset: Vec2D{x:-35.0,y:20.0}, 
    };
    */
    let grid = Grid::Hex {
        pitch: 1.0/7.5e-3,
        rotation: 0.0,
        offset: Vec2D{x:0.0, y:0.0},
    };

    // compute all centroids, indexed by pinhole id and exposure id,
    eprintln!("computing centroids");
    const RAD: usize = 50;
    let cogs: Vec<Centroid> = images.par_iter().flat_map(|image| {
        let grid_shifted = grid + image.shift.unwrap();
        image.cogs(&grid_shifted, RAD).into_iter().filter_map(|cog|
            {
                //println!("{:?}", cog.flux);
                if cog.flux > FLUXTHRESH {
                    Some(cog)
                } else {
                    None
                }
            }
        ).collect::<Vec<Centroid>>()
    }).collect();

    
    let test_img = images[0].clone();
    let test_grid = grid + test_img.shift.unwrap();
    test_img.draw_on_circles(&test_grid, RAD as f64, 5_000.0)
    .to_fits("delete.fits")?;
    
    eprintln!("solving for distortions");
    let mut distortions = BiVarPolyDistortions::new(3);
    distortions.solve_for_coeffs(cogs)?;

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












#[derive(Clone,Debug,Copy)]
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

impl Mul<f64> for Vec2D {
    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
    type Output = Self;
}


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
    pub fn all_points(&self, image: &Image) -> Vec<Vec2D> {
        let width = image.shape[1];
        let height = image.shape[0];
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
pub struct Centroid {
    pub cogx: f64,
    pub cogy: f64,
    pub flux: f64,
    pub x: f64,
    pub y: f64,
}
