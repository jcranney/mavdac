use basis::{BiVarPolyDistortions, DistortionBasis};
use pyo3::prelude::*;
use core::f64;
use rayon::prelude::*;

mod errors;
pub mod basis;
pub mod io;
pub mod geom;
pub use crate::io::{Image, Coordinate};
pub use crate::errors::{MavDACError, Result};
pub use crate::geom::{Centroid,Vec2D,Grid};

/// A Python module implemented in Rust.
#[pymodule]
fn mavdac(m: &Bound<'_, PyModule>) -> PyResult<()> {
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

    eprintln!("loading images");
    let images = load_images(pattern)?;
    let shape = images[0].shape;

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

    eprintln!("loading images");
    let images = load_images(pattern)?;
    let shape = images[0].shape;

    
    // check if grid.yaml exists
    let grid = match Grid::from_yaml("grid.yaml") {
        Ok(grid) => grid,
        Err(MavDACError::YAMLError(e)) => {
            // file is badly formatted
            // exit with error:
            return Err(PyErr::from(MavDACError::YAMLError(e)));
        },
        Err(MavDACError::IOError(e)) => {
            match e.kind() {
                std::io::ErrorKind::NotFound => (),
                _ => {
                    eprintln!("{}", e);
                    return Err(MavDACError::IOError(e).into())
                },
            }
            // file doesn't exist, create one with defaults:
            let grid = Grid::Hex {
                pitch: 1.0/7.5e-3,
                rotation: 0.0,
                offset: Vec2D{x:0.0, y:0.0},
            };
            grid.to_yaml("grid.yaml")?;
            grid
        },
        Err(e) => {return Err(e.into());},
    };


    // do centroids for base set of pinholes over all exposures, then filter out
    // any that are lower than FLUXTHRESH in any exposure
    const RAD: usize = 90;

    // create sample image to see if fit is good
    images[0].clone()
    .draw_on_circles(&(grid+images[0].shift), RAD as f64, 30_000.0)
    .to_fits("sample.fits")?;
    eprintln!("? check sample.fits to verify grid alignment\n? modify grid.yaml if not");
    
    let pinholes = grid.all_points(images[0].shape[1], images[0].shape[0]);
    // cogs should be a n_pinholes x n_images 
    let cogs: Vec<Option<Vec<Centroid>>> = pinholes.iter().map(|pinhole|
        images.iter().map(|image|
            image.cog(&(*pinhole+image.shift), RAD)
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
                cleaned_cog.pos = mean_pinhole_position + image.shift;
                cleaned_cogs.push(cleaned_cog);
            }
        }
    }
 
    eprintln!("solving for distortions");
    let distortions = solve_distortions(cleaned_cogs, 7, shape)?;

    match coords {
        Some(mut coords) => {
            for coord in &mut coords {
                coord.dist = Some(distortions.eval(&coord.pos));
                println!("{}", coord);
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
        if !line.is_empty() {
            coords.push(Coordinate::try_from(line)?);
        }
    }
    Ok(coords)
}

#[pyfunction]
fn load_images(pattern: &str) -> Result<Vec<Image>> {
    glob::glob(pattern)?
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

