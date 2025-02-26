//! # mavdac
//! MAVIS Differential Astrometric Calibrator. This crate provides the low-level
//! functions used by the high-level Python API (see
//! [mavdac](https://github.com/jcranney/mavdac) on github). The intented usage
//! of this crate is through the Python API, but it may be possible to build
//! other calibration tools on top of this crate.
//! ## Examples
//! ```
//! use mavdac::{Image,Grid,Vec2D};
//! let imgs: Vec<Image> = mavdac::load_images("./imgs/img_*.fits").unwrap();
//! let grid = Grid::Hex {
//!     pitch: 100.0,  // pixels
//!     rotation: 0.0,  // radians
//!     offset: Vec2D{x:0.0,y:0.0}  // pixels
//! };
//! let cogs = mavdac::measure_cogs(
//!     imgs,  // vector of images
//!     grid,  // grid geometry
//!     10,  // radius for centroider
//!     10_000.0,  // flux threshold for "valid" cogs
//! );
//! // all of the remaining tasks are done in python, since numpy.linalg is more
//! // reliable than any rust linalg solution I tried. 
//! ```
use pyo3::prelude::*;
use core::f64;
use std::path::PathBuf;
use rayon::prelude::*;

mod errors;
mod basis;
mod io;
mod geom;
pub use crate::io::{Image, Coordinate};
pub use crate::errors::{MavDACError, Result};
pub use crate::geom::{Centroid,Vec2D,Grid};
pub use crate::basis::{BiVarPolyDistortions,BiVarFourierDistortions,DistortionBasis};

/// A Python module implemented in Rust.
#[pymodule]
fn mavdac(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_coordinates, m)?)?;
    m.add_function(wrap_pyfunction!(load_images, m)?)?;
    m.add_function(wrap_pyfunction!(measure_cogs, m)?)?;
    m.add_class::<Image>()?;
    m.add_class::<Grid>()?;
    m.add_class::<Centroid>()?;
    m.add_class::<BiVarPolyDistortions>()?;
    m.add_class::<BiVarFourierDistortions>()?;
    m.add_class::<Coordinate>()?;
    m.add_class::<Vec2D>()?;
    Ok(())
}

/// Load images from disk given a glob pattern
#[pyfunction]
pub fn load_images(pattern: &str) -> Result<Vec<Image>> {
    glob::glob(pattern)?
    .filter(|file| file.is_ok())
    .flatten()
    .collect::<Vec<PathBuf>>()
    .into_par_iter()
    .map(|path| Image::from_fits(path.to_str().unwrap()))
    .collect::<Result<Vec<Image>>>()
}

/// measured centroids from a set of images
#[pyfunction]
pub fn measure_cogs(
    images: Vec<Image>, grid: Grid, rad: usize, fluxthresh: f64
) -> Vec<Vec<Centroid>> {
    if images.is_empty() {
        return vec![];
    }
    let pinholes = grid.all_points(images[0].shape[1], images[0].shape[0]);
    
    // cogs should be a n_pinholes x n_images 
    let mut cogs: Vec<Option<Vec<Centroid>>> = pinholes.par_iter().map(|pinhole|
        images.iter().map(|image|
            image.cog(&(*pinhole+image.shift), rad)
        ).collect::<Vec<Centroid>>()
    ).map(|pinhole_cogs|
        if pinhole_cogs.iter().all(|cog|
            cog.flux > fluxthresh
        ) {
            Some(pinhole_cogs)
        } else {
            None
        }
    ).collect();
    for (i, _pinhole) in pinholes.iter().enumerate() {
        if let Some(pinhole_cogs) = &mut cogs[i] {
            let mut mean_pinhole_position = Vec2D {x: 0.0, y: 0.0};
            for (j, image) in images.iter().enumerate() {
                let a = pinhole_cogs[j].cog - image.shift;
                mean_pinhole_position += a;
            }
            mean_pinhole_position.x /= images.len() as f64;
            mean_pinhole_position.y /= images.len() as f64;
            for (j, image) in images.iter().enumerate() {
                pinhole_cogs[j].pos = mean_pinhole_position + image.shift;
            }
        }
    }
    cogs.into_iter().flatten().collect()
}

/// read coordinates from file
#[pyfunction]
pub fn get_coordinates(filename: &str) -> PyResult<Vec<Coordinate>> {
    let contents = std::fs::read_to_string(filename)?;
    let mut coords: Vec<Coordinate> = vec![];
    for line in contents.split('\n') {
        if !line.is_empty() {
            coords.push(Coordinate::try_from(line)?);
        }
    }
    Ok(coords)
}

