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
use core::f64;
use pyo3::prelude::*;

mod errors;
pub use crate::errors::MavDACError;
pub use crate::mavdac::{Vec2D, Grid, Image, measure_cogs, load_images};



#[derive(Debug, Clone)]
struct Pixel {
    x: usize,
    y: usize,
    val: f64,
}

/// A Python module implemented in Rust.
#[pymodule]
pub mod mavdac {
    use fitrs::{Fits, FitsData, Hdu, HeaderValue};
    use rayon::prelude::*;
    use serde::{Deserialize, Serialize};
    use std::{
        f64::consts::PI,
        fmt::Display,
        io::Write,
        ops::{Add, AddAssign, Mul, Sub},
        path::PathBuf,
    };

    use pyo3::prelude::*;

    use crate::{DistortionBasis, MavDACError, Pixel};

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

    /// Load images from disk given a glob pattern
    #[pyfunction]
    pub fn load_images(pattern: &str) -> PyResult<Vec<Image>> {
        glob::glob(pattern)
            .map_err(MavDACError::from)?
            .filter(|file| file.is_ok())
            .flatten()
            .collect::<Vec<PathBuf>>()
            .into_par_iter()
            .map(|path| Image::from_fits(path.to_str().unwrap()))
            .collect::<PyResult<Vec<Image>>>()
    }

    /// measured centroids from a set of images
    #[pyfunction]
    pub fn measure_cogs(
        images: Vec<Image>,
        grid: Grid,
        rad: usize,
        fluxthresh: f64,
    ) -> Vec<Vec<Centroid>> {
        if images.is_empty() {
            return vec![];
        }
        let pinholes = grid.all_points(images[0].shape[1], images[0].shape[0]);

        // cogs should be a n_pinholes x n_images
        let mut cogs: Vec<Option<Vec<Centroid>>> = pinholes
            .par_iter()
            .map(|pinhole| {
                images
                    .iter()
                    .map(|image| image.cog(&(*pinhole + image.shift), rad))
                    .collect::<Vec<Centroid>>()
            })
            .map(|pinhole_cogs| {
                if pinhole_cogs.iter().all(|cog| cog.flux > fluxthresh) {
                    Some(pinhole_cogs)
                } else {
                    None
                }
            })
            .collect();
        for (i, _pinhole) in pinholes.iter().enumerate() {
            if let Some(pinhole_cogs) = &mut cogs[i] {
                let mut mean_pinhole_position = Vec2D { x: 0.0, y: 0.0 };
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

    /// 2D vector, corresponding to float-valued pixel positions
    #[derive(Debug, Deserialize, PartialEq, Serialize, Clone, Copy)]
    #[pyclass(from_py_object)]
    pub struct Vec2D {
        pub x: f64,
        pub y: f64,
    }

    #[pymethods]
    impl Vec2D {
        #[new]
        fn new(x: f64, y: f64) -> Vec2D {
            Vec2D { x, y }
        }
        #[getter]
        fn x(&self) -> f64 {
            self.x
        }
        #[getter]
        fn y(&self) -> f64 {
            self.y
        }
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

    /// Grid type, defined from minimal parameters but able to determine all possible
    /// pinhole positions.
    #[pyclass(from_py_object)]
    #[derive(Debug, Serialize, Deserialize, PartialEq, Clone, Copy)]
    pub enum Grid {
        /// Hexagonal grid, with a pinhole at the centre of the image by default
        Hex {
            /// separation between adjacent pinholes (in pixels)
            pitch: f64, // pixels
            /// rotation of pinhole grid around centre of image (in radians)
            rotation: f64, // radians
            /// offset of pinhole grid (if {0,0} then there is a pinhole at the centre of the image)
            offset: Vec2D, // pixels
        },
    }

    impl Add<Vec2D> for Grid {
        fn add(self, rhs: Vec2D) -> Self::Output {
            match self {
                Grid::Hex {
                    pitch,
                    rotation,
                    offset,
                } => Grid::Hex {
                    pitch,
                    rotation,
                    offset: offset + rhs,
                },
            }
        }
        type Output = Self;
    }

    #[pymethods]
    impl Grid {
        /// load grid from yaml file
        #[new]
        pub fn from_yaml(filename: &str) -> PyResult<Grid> {
            let f = std::fs::File::open(filename)?;
            let grid: Grid = serde_yaml::from_reader(f).map_err(MavDACError::from)?;
            Ok(grid)
        }
        /// save grid to yaml file
        pub fn to_yaml(&self, filename: &str) -> PyResult<()> {
            let mut f = std::fs::File::create(filename)?;
            write!(
                f,
                "{}",
                serde_yaml::to_string(&self).map_err(MavDACError::from)?
            )?;
            Ok(())
        }
        /// determine all possible points of pinholes for a given grid
        pub fn all_points(&self, width: usize, height: usize) -> Vec<Vec2D> {
            let max_rad = width.max(height) * 4;
            match self {
                Grid::Hex {
                    pitch,
                    rotation,
                    offset,
                } => {
                    // first make a square grid with way too many points
                    (0..2 * max_rad)
                        .map(|x| x as f64)
                        .flat_map(|x| {
                            (0..2 * max_rad)
                                // shift it to be centered at origin
                                .map(|y| y as f64)
                                .map(move |y| (x - max_rad as f64, y - max_rad as f64))
                        })
                        // scale it to pixel units
                        .map(|(x, y)| (x * pitch, y * pitch))
                        // then map it to a hex grid with gaps
                        .map(|(x, y)| (x + 0.5 * y, y * (3.0f64).sqrt() / 2.0))
                        // rotate
                        .map(|(x, y)| {
                            (
                                x * rotation.cos() - y * rotation.sin(),
                                x * rotation.sin() + y * rotation.cos(),
                            )
                        })
                        // now apply the offset:
                        .map(|(x, y)| (x + offset.x, y + offset.y))
                        // shift back to valid pixel range
                        .map(|(x, y)| (x + (width / 2) as f64 - 0.5, y + (height / 2) as f64 - 0.5))
                        .filter(|(x, y)| {
                            *x >= 0.0 && *x < width as f64 && *y >= 0.0 && *y < height as f64
                        })
                        .map(|(x, y)| Vec2D { x, y })
                        .collect()
                }
            }
        }
    }

    /// centroid type, to be populated by centroider
    #[derive(Debug)]
    #[pyclass]
    pub struct Centroid {
        pub cog: Vec2D,
        pub flux: f64,
        pub pos: Vec2D,
    }

    #[pymethods]
    impl Centroid {
        /// measured centroid x coordinate
        #[getter]
        pub fn cogx(&self) -> f64 {
            self.cog.x
        }
        /// measured centroid y coordinate
        #[getter]
        pub fn cogy(&self) -> f64 {
            self.cog.y
        }
        /// nominal centroid x coordinate
        #[getter]
        pub fn posx(&self) -> f64 {
            self.pos.x
        }
        /// nominal centroid y coordinate
        #[getter]
        pub fn posy(&self) -> f64 {
            self.pos.y
        }
        /// flux of centroid (summed over all valid pixels)
        #[getter]
        pub fn flux(&self) -> f64 {
            self.flux
        }
    }

    /// Image struct, with metadata corresponding to calibration
    #[derive(Clone, Debug)]
    #[pyclass(from_py_object)]
    pub struct Image {
        pub data: Vec<f64>,
        pub shift: Vec2D,
        pub shape: [usize; 2], // numpy style
    }

    impl std::fmt::Display for Image {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(
                f,
                ", xshift {:0.4}, yshift {:0.4}",
                self.shift.x, self.shift.y
            )
        }
    }

    #[pymethods]
    impl Image {
        /// load image from fits file
        #[new]
        pub fn from_fits(filename: &str) -> PyResult<Image> {
            let fits = Fits::open(filename)?;
            match fits.get(0) {
                Some(hdu) => {
                    match hdu.value("NAXIS") {
                        Some(HeaderValue::IntegerNumber(2)) => (),
                        _ => {
                            return Err(
                                MavDACError::InvalidFITS("expected NAXIS==2".to_string()).into()
                            );
                        }
                    };
                    let mut shape: [usize; 2] = [0, 0];
                    match hdu.value("NAXIS2") {
                        Some(HeaderValue::IntegerNumber(x)) if *x > 0 => {
                            shape[0] = *x as usize;
                        }
                        _ => {
                            return Err(
                                MavDACError::InvalidFITS("invalid NAXIS2".to_string()).into()
                            );
                        }
                    }
                    match hdu.value("NAXIS1") {
                        Some(HeaderValue::IntegerNumber(x)) if *x > 0 => {
                            shape[1] = *x as usize;
                        }
                        _ => {
                            return Err(
                                MavDACError::InvalidFITS("invalid NAXIS1".to_string()).into()
                            );
                        }
                    };
                    let shift =
                        {
                            let xshift: f64 =
                                match hdu.value("XSHIFT").ok_or(MavDACError::InvalidFITS(
                                    format!("missing XSHIFT in fits header {}", filename),
                                ))? {
                                    HeaderValue::IntegerNumber(a) => *a as f64,
                                    HeaderValue::RealFloatingNumber(a) => *a,
                                    _ => {
                                        return Err(MavDACError::InvalidFITS(format!(
                                            "XSHIFT in fits header has invalid datatype, \
                                must be float or int {}",
                                            filename
                                        ))
                                        .into());
                                    }
                                };
                            let yshift: f64 =
                                match hdu.value("YSHIFT").ok_or(MavDACError::InvalidFITS(
                                    format!("missing YSHIFT in fits header {}", filename),
                                ))? {
                                    HeaderValue::IntegerNumber(a) => *a as f64,
                                    HeaderValue::RealFloatingNumber(a) => *a,
                                    _ => {
                                        return Err(MavDACError::InvalidFITS(format!(
                                            "YSHIFT in fits header has invalid datatype, \
                                must be float or int\n{}",
                                            filename
                                        ))
                                        .into());
                                    }
                                };
                            Vec2D {
                                x: xshift,
                                y: yshift,
                            }
                        };
                    let data: Vec<f64> = match hdu.read_data() {
                        FitsData::IntegersI32(array) => {
                            array.data.iter().flatten().map(|x| *x as f64).collect()
                        }
                        FitsData::IntegersU32(array) => {
                            array.data.iter().flatten().map(|x| *x as f64).collect()
                        }
                        FitsData::FloatingPoint32(array) => {
                            array.data.iter().map(|x| *x as f64).collect()
                        }
                        FitsData::FloatingPoint64(array) => array.data,
                        FitsData::Characters(array) => {
                            array.data.iter().map(|x| *x as u8 as f64).collect()
                        }
                    };
                    Ok(Image { data, shape, shift })
                }
                _ => {
                    Err(MavDACError::InvalidFITS(format!("no primary hdu in {}", filename)).into())
                }
            }
        }

        /// save image to fits file
        pub fn to_fits(&self, filename: &str) -> PyResult<()> {
            let primary_hdu = Hdu::new(&self.shape, self.data.clone());
            Fits::create(filename, primary_hdu)?;
            Ok(())
        }

        /// draw circles around COG region defined by grid (e.g., for alignment)
        pub fn draw_on_circles(&mut self, grid: &Grid, rad: f64, val: f64) {
            const NTHETA: usize = 1000;
            grid.all_points(self.shape[1], self.shape[0])
                .into_iter()
                .map(|v| v + self.shift)
                .flat_map(|v| {
                    (0..1000)
                        .map(|i| i as f64 / NTHETA as f64)
                        .map(|t| t * 2.0 * PI)
                        .map(move |theta| (v.x + theta.cos() * rad, v.y + theta.sin() * rad))
                        .map(|(x, y)| (x as usize, y as usize))
                        .filter(|(x, y)| *x < self.shape[1] && *y < self.shape[0])
                })
                .for_each(|(x, y)| {
                    self.data[y * self.shape[1] + x] += val;
                });
        }

        /// compute centroids for image given a grid and cog-radius
        pub fn cogs(&self, grid: &Grid, rad: usize) -> Vec<Centroid> {
            // get all nominal positions
            let points = grid.all_points(self.shape[1], self.shape[0]);

            // measure cog and intensity within radius at all points
            points
                .into_iter()
                .map(|v| v + self.shift)
                .map(|point| self.cog(&point, rad))
                .collect()
        }

        /// compute centroid for image given a point and cog-radius
        pub fn cog(&self, point: &Vec2D, rad: usize) -> Centroid {
            let (sumx, sumy, flux) = self
                .get_blob(&point.clone(), rad)
                .into_iter()
                .map(|pixel| {
                    (
                        pixel.x as f64 * pixel.val,
                        pixel.y as f64 * pixel.val,
                        pixel.val,
                    )
                })
                .fold((0.0, 0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2));
            Centroid {
                cog: Vec2D {
                    x: sumx / flux,
                    y: sumy / flux,
                },
                flux,
                pos: Vec2D {
                    x: point.x,
                    y: point.y,
                },
            }
        }

        /// get shape of image
        #[getter]
        fn shape(&self) -> (usize, usize) {
            self.shape.into()
        }
    }

    impl Image {
        fn get_blob(&self, pos: &Vec2D, rad: usize) -> Vec<Pixel> {
            let rad = rad as isize;
            let mut pixels: Vec<Pixel> = vec![];
            let xc = pos.x as isize;
            let yc = pos.y as isize;
            for x in -rad..rad + 1 {
                for y in -rad..rad + 1 {
                    if x.pow(2) + y.pow(2) > rad.pow(2) {
                        continue;
                    }
                    let px = (x + xc) as usize;
                    let py = (y + yc) as usize;
                    if px >= self.shape[1] || py >= self.shape[0] {
                        continue;
                    }
                    pixels.push(Pixel {
                        x: px,
                        y: py,
                        val: self.data[py * self.shape[1] + px],
                    });
                }
            }
            pixels
        }
    }

    /// coordinate struct for interfacing with coordinate files
    #[derive(Debug, Clone)]
    #[pyclass(from_py_object)]
    pub struct Coordinate {
        pub pos: Vec2D,
    }

    #[pymethods]
    impl Coordinate {
        /// get coordinate position
        #[getter]
        fn pos(&self) -> (f64, f64) {
            (self.pos.x, self.pos.y)
        }
    }

    impl TryFrom<&str> for Coordinate {
        /// parse string into coordinate
        fn try_from(s: &str) -> std::result::Result<Self, Self::Error> {
            let mut split = s.split(',');
            let x: f64;
            let y: f64;
            if let Some(a) = split.next() {
                x = match a.parse::<f64>() {
                    Ok(x) => x,
                    Err(_) => {
                        return Err(MavDACError::Coordinate(format!(
                            "failed to parse x-ordinate: {}",
                            a
                        )));
                    }
                };
            } else {
                return Err(MavDACError::Coordinate("missing x ordinate".to_string()));
            }
            if let Some(a) = split.next() {
                y = match a.parse::<f64>() {
                    Ok(y) => y,
                    Err(_) => {
                        return Err(MavDACError::Coordinate(format!(
                            "failed to parse y-ordinate: {}",
                            a
                        )));
                    }
                };
            } else {
                return Err(MavDACError::Coordinate("missing y ordinate".to_string()));
            }
            Ok(Self {
                pos: Vec2D { x, y },
            })
        }

        type Error = MavDACError;
    }

    impl Display for Coordinate {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "{},{},", self.pos.x, self.pos.y)?;
            Ok(())
        }
    }

    /// Bivariate Homogenous Polynomial to be used as distortion basis function
    ///
    /// See [wikipedia article](https://en.wikipedia.org/wiki/Homogeneous_polynomial) for
    /// more info.
    #[pyclass]
    pub struct BiVarPolyDistortions {
        /// degree of polynomial (maximum order)
        pub degree: usize,
        /// coefficients of distortions
        pub coeffs: Vec<Vec2D>,
        /// shape of image (numpy format)
        pub shape: [usize; 2],
        nk_lut: Vec<(usize, usize)>,
    }

    #[pymethods]
    impl BiVarPolyDistortions {
        /// construct a new set of bivariate homogenous polynomials
        #[new]
        pub fn new(degree: usize, shape: [usize; 2]) -> Self {
            let coeffs = vec![Vec2D { x: 0.0, y: 0.0 }; ((degree + 1) * (degree + 2)) / 2 - 1];
            let ncoeffs = coeffs.len();
            Self {
                degree,
                coeffs,
                shape,
                nk_lut: (0..ncoeffs).map(Self::index_to_nk).collect(),
            }
        }

        #[staticmethod]
        fn index_to_nk(index: usize) -> (usize, usize) {
            let n = ((-1.0 + (1.0 + 8.0 * ((index + 1) as f64)).sqrt()) / 2.0).floor() as usize;
            let k = index + 1 - n * (n + 1) / 2;
            (n, k)
        }

        /// sample basis function given index at x/y coordinates
        pub fn sample_xy(&self, x: f64, y: f64, ell: usize) -> f64 {
            self.sample(&Vec2D { x, y }, ell)
        }

        /// evaluate distortions (including coefficients) at x/y coordinates
        pub fn eval_xy(&self, x: f64, y: f64) -> (f64, f64) {
            let Vec2D { x, y } = self.eval(&Vec2D { x, y });
            (x, y)
        }

        #[getter]
        fn ncoeffs(&self) -> usize {
            self.coeffs.len()
        }

        #[getter]
        fn coeffs(&self) -> Vec<Vec<f64>> {
            self.coeffs
                .clone()
                .into_iter()
                .map(|v| {
                    let Vec2D { x, y } = v;
                    vec![x, y]
                })
                .collect::<Vec<Vec<f64>>>()
        }

        /// load coefficients (e.g.) from python
        pub fn load_coeffs(&mut self, coeffs: Vec<Vec<f64>>) {
            self.set_coeffs(
                coeffs
                    .into_iter()
                    .map(|p| Vec2D { x: p[0], y: p[1] })
                    .collect(),
            );
        }
    }

    impl DistortionBasis for BiVarPolyDistortions {
        fn sample(&self, pos: &Vec2D, index: usize) -> f64 {
            let (n, k) = self.nk_lut[index];
            let &Vec2D { mut x, mut y } = pos;
            x -= (self.shape[1] as f64) / 2.0;
            y -= (self.shape[0] as f64) / 2.0;
            x /= self.shape[1] as f64;
            y /= self.shape[0] as f64;
            x.powf(k as f64) * y.powf((n - k) as f64)
        }

        fn get_coeffs(&self) -> &Vec<Vec2D> {
            &self.coeffs
        }

        fn set_coeffs(&mut self, coeffs: Vec<Vec2D>) {
            self.coeffs = coeffs;
        }
    }

    /// Bivariate Fourier-Based functions to be used as distortion basis function
    #[pyclass]
    pub struct BiVarFourierDistortions {
        /// maximum frequency (0.5*periods per field)
        pub max_freq: usize,
        /// coefficients of distortions
        pub coeffs: Vec<Vec2D>,
        /// shape of image (numpy format)
        pub shape: [usize; 2],
    }

    #[pymethods]
    impl BiVarFourierDistortions {
        /// construct a new set of bivariate homogenous polynomials
        #[new]
        pub fn new(max_freq: usize, shape: [usize; 2]) -> Self {
            let coeffs = vec![Vec2D { x: 0.0, y: 0.0 }; 4 * (max_freq).pow(2)];
            Self {
                max_freq,
                coeffs,
                shape,
            }
        }

        /// sample basis function given index at x/y coordinates
        pub fn sample_xy(&self, x: f64, y: f64, ell: usize) -> f64 {
            self.sample(&Vec2D { x, y }, ell)
        }

        /// evaluate distortions (including coefficients) at x/y coordinates
        pub fn eval_xy(&self, x: f64, y: f64) -> (f64, f64) {
            let Vec2D { x, y } = self.eval(&Vec2D { x, y });
            (x, y)
        }

        #[getter]
        fn ncoeffs(&self) -> usize {
            self.coeffs.len()
        }

        #[getter]
        fn coeffs(&self) -> Vec<Vec<f64>> {
            self.coeffs
                .clone()
                .into_iter()
                .map(|v| {
                    let Vec2D { x, y } = v;
                    vec![x, y]
                })
                .collect::<Vec<Vec<f64>>>()
        }

        /// load coefficients (e.g.) from python
        pub fn load_coeffs(&mut self, coeffs: Vec<Vec<f64>>) {
            self.set_coeffs(
                coeffs
                    .into_iter()
                    .map(|p| Vec2D { x: p[0], y: p[1] })
                    .collect(),
            );
        }
    }

    impl DistortionBasis for BiVarFourierDistortions {
        fn sample(&self, pos: &Vec2D, index: usize) -> f64 {
            let &Vec2D { mut x, mut y } = pos;
            // x -= (self.shape[1] as f64)/2.0;
            // y -= (self.shape[0] as f64)/2.0;
            x /= self.shape[1] as f64;
            y /= self.shape[0] as f64;
            let freq_x: f64 = 1.0 * PI * (((index / 4) / self.max_freq) % self.max_freq) as f64;
            let freq_y: f64 = 1.0 * PI * ((index / 4) % self.max_freq) as f64;
            match index % 4 {
                0 => (freq_x * x).cos() * (freq_y * y).cos(),
                1 => (freq_x * x).cos() * (freq_y * y).sin(),
                2 => (freq_x * x).sin() * (freq_y * y).cos(),
                3 => (freq_x * x).sin() * (freq_y * y).sin(),
                _ => unreachable!(),
            }
        }

        fn get_coeffs(&self) -> &Vec<Vec2D> {
            &self.coeffs
        }

        fn set_coeffs(&mut self, coeffs: Vec<Vec2D>) {
            self.coeffs = coeffs;
        }
    }
}

/// Trait that allows standard evaluation of distortion functions
pub trait DistortionBasis {
    /// the sample function returns the value of the bivariate function with a given
    /// index, at coordinate (x,y)
    fn sample(&self, pos: &Vec2D, index: usize) -> f64;
    /// a method to get the coefficients vector
    fn get_coeffs(&self) -> &Vec<Vec2D>;
    /// a method to set the coefficients vector
    fn set_coeffs(&mut self, coeffs: Vec<Vec2D>);

    /// evaluate the basis function from the coefficients
    fn eval(&self, pos: &Vec2D) -> Vec2D {
        (0..self.get_coeffs().len())
            .map(|index| self.get_coeffs()[index] * self.sample(pos, index))
            .fold(Vec2D { x: 0.0, y: 0.0 }, |a, b| a + b)
    }
}
