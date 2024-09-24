use std::ops::{Add, AddAssign, Mul, Sub};

use pyo3::pyclass;
use rustfft::num_traits::Float;
use serde::{Serialize,Deserialize};

#[derive(Clone,Debug,Copy,Deserialize,PartialEq,Serialize)]
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
#[derive(Clone,Copy,Debug,Serialize,Deserialize,PartialEq)]
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