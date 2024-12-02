use std::ops::{Add, AddAssign, Mul, Sub};
use std::io::Write;

use pyo3::{pyclass, pymethods};
use rustfft::num_traits::Float;
use serde::{Serialize,Deserialize};
use crate::Result;

/// 2D vector, corresponding to float-valued pixel positions
#[derive(Clone,Debug,Copy,Deserialize,PartialEq,Serialize)]
#[pyclass]
pub struct Vec2D {
    pub x: f64,
    pub y: f64,
}

#[pymethods]
impl Vec2D {
    #[new]
    fn new(x: f64, y: f64) -> Vec2D {
        Vec2D{x, y}
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
#[pyclass]
#[derive(Clone,Copy,Debug,Serialize,Deserialize,PartialEq)]
pub enum Grid {
    /// Hexagonal grid, with a pinhole at the centre of the image by default
    Hex {
        /// separation between adjacent pinholes (in pixels)
        pitch: f64,  // pixels
        /// rotation of pinhole grid around centre of image (in radians)
        rotation: f64,  // radians
        /// offset of pinhole grid (if {0,0} then there is a pinhole at the centre of the image)
        offset: Vec2D,  // pixels
    },
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


#[pymethods]
impl Grid {
    /// load grid from yaml file
    #[new]
    pub fn from_yaml(filename: &str) -> Result<Grid> {
        let f = std::fs::File::open(filename)?;
        let grid: Grid = serde_yaml::from_reader(f)?;
        Ok(grid)
    }
    /// save grid to yaml file
    pub fn to_yaml(&self, filename: &str) -> Result<()> {
        let mut f = std::fs::File::create(filename)?;
        write!(f,"{}",serde_yaml::to_string(&self)?)?;
        Ok(())
    }
    /// determine all possible points of pinholes for a given grid
    pub fn all_points(&self, width: usize, height: usize) -> Vec<Vec2D> {
        let max_rad = width.max(height)*4;
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

/// centroid type, to be populated by centroider
#[derive(Debug, Clone)]
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