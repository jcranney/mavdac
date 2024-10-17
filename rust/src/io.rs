use std::{f64::consts::PI, fmt::Display};

use fitrs::{Fits, FitsData, Hdu, HeaderValue};
use pyo3::{pyclass, pymethods};

use crate::{Centroid, Grid, MavDACError, Result, Vec2D};


#[derive(Debug, Clone)]
struct Pixel {
    x: usize,
    y: usize,
    val: f64,
}

/// Image struct, with metadata corresponding to calibration
#[derive(Clone,Debug)]
#[pyclass]
pub struct Image {
    pub data: Vec<f64>,
    pub shift: Vec2D,
    pub shape: [usize;2], // numpy style
}

impl std::fmt::Display for Image {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, ", xshift {:0.4}, yshift {:0.4}", self.shift.x, self.shift.y)
    }
}


#[pymethods]
impl Image {
    /// load image from fits file
    #[new]
    pub fn from_fits(filename: &str) -> Result<Image> {
        let fits = Fits::open(filename)?;
        if let Some(hdu) = fits.get(0) {
            match hdu.value("NAXIS")  {
                Some(HeaderValue::IntegerNumber(2)) => (),
                _ => return Err(MavDACError::InvalidFITS("expected NAXIS==2".to_string())),
            };
            let mut shape: [usize;2] = [0,0];
            match hdu.value("NAXIS2")  {
                Some(HeaderValue::IntegerNumber(x)) if *x > 0 => {
                    shape[0] = *x as usize;
                },
                _ => return Err(MavDACError::InvalidFITS("invalid NAXIS2".to_string())),
            }
            match hdu.value("NAXIS1")  {
                Some(HeaderValue::IntegerNumber(x)) if *x > 0 => {
                    shape[1] = *x as usize;
                },
                _ => return Err(MavDACError::InvalidFITS("invalid NAXIS1".to_string())),
            };
            let shift = {
                let xshift: f64 = match hdu.value("XSHIFT").ok_or(MavDACError::InvalidFITS(
                    format!("missing XSHIFT in fits header {}", filename)
                ))? {
                    HeaderValue::IntegerNumber(a) => *a as f64,
                    HeaderValue::RealFloatingNumber(a) => *a,
                    _ => return Err(MavDACError::InvalidFITS(
                        format!("XSHIFT in fits header has invalid datatype, \
                                must be float or int {}", filename)
                    ))
                };
                let yshift: f64 = match hdu.value("YSHIFT").ok_or(MavDACError::InvalidFITS(
                    format!("missing YSHIFT in fits header {}", filename)
                ))? {
                    HeaderValue::IntegerNumber(a) => *a as f64,
                    HeaderValue::RealFloatingNumber(a) => *a,
                    _ => return Err(MavDACError::InvalidFITS(
                        format!("YSHIFT in fits header has invalid datatype, \
                                must be float or int\n{}", filename)
                    ))
                };
                Vec2D {
                    x: xshift,
                    y: yshift,
                }
            };
            let data: Vec<f64> = match hdu.read_data() {
                FitsData::IntegersI32(array) => {
                    array.data.iter().flatten().map(|x| *x as f64).collect()
                },
                FitsData::IntegersU32(array) => {
                    array.data.iter().flatten().map(|x| *x as f64).collect()
                },
                FitsData::FloatingPoint32(array) => {
                    array.data.iter().map(|x| *x as f64).collect()
                },
                FitsData::FloatingPoint64(array) => {
                    array.data
                },
                FitsData::Characters(array) => {
                    array.data.iter().map(|x| *x as u8 as f64).collect()
                },
            };
            Ok(Image { data, shape, shift })
        } else {
            Err (
                MavDACError::InvalidFITS(
                    format!("no primary hdu in {}", &filename)
                )
            )
        }
    }

    /// save image to fits file
    pub fn to_fits(&self, filename: &str) -> Result<()> 
    {
        let primary_hdu = Hdu::new(&self.shape, self.data.clone());
        Fits::create(filename, primary_hdu)?;
        Ok(())
    }
    
    /// draw circles around COG region defined by grid (e.g., for alignment)
    pub fn draw_on_circles(&mut self, grid: &Grid, rad: f64, val: f64) {
        const NTHETA: usize = 1000;
        grid.all_points(self.shape[1], self.shape[0])
        .into_iter()
        .map(|v| v+self.shift)
        .flat_map(|v| {
            (0..1000).map(|i| i as f64 / NTHETA as f64)
            .map(|t| t*2.0*PI)
            .map(move |theta| (v.x + theta.cos()*rad, v.y + theta.sin()*rad))
            .map(|(x,y)| (x as usize, y as usize))
            .filter(|(x,y)| *x < self.shape[1] && *y < self.shape[0])
        })
        .for_each(|(x,y)| {
            self.data[y*self.shape[1]+x] += val;
        });
    }

    /// compute centroids for image given a grid and cog-radius
    pub fn cogs(&self, grid: &Grid, rad: usize) -> Vec<Centroid> {
        // get all nominal positions
        let points = grid.all_points(self.shape[1], self.shape[0]);

        // measure cog and intensity within radius at all points
        points.into_iter().map(|v| v+self.shift).map(|point| self.cog(&point, rad)).collect()
    }

    /// compute centroid for image given a point and cog-radius
    pub fn cog(&self, point: &Vec2D, rad: usize) -> Centroid {
        let (sumx,sumy,flux) = self.get_blob(&point.clone(), rad).into_iter()
        .map(|pixel| (
            pixel.x as f64*pixel.val,
            pixel.y as f64*pixel.val,
            pixel.val
        )).fold((0.0,0.0,0.0), |a,b| (a.0+b.0, a.1+b.1, a.2+b.2));
        Centroid {
            cog: Vec2D{x: sumx / flux, y: sumy /flux},
            flux,
            pos: Vec2D{x: point.x, y: point.y},
        }
    }

    /// get shape of image
    #[getter]
    fn shape(&self) -> (usize,usize) {
        self.shape.into()
    }
}

impl Image {
    fn get_blob(&self, pos: &Vec2D, rad: usize) -> Vec<Pixel> {
        let rad = rad as isize;
        let mut pixels: Vec<Pixel> = vec![];
        let xc = pos.x as isize;
        let yc = pos.y as isize;
        for x in -rad..rad+1 {
            for y in -rad..rad+1 {
                if x.pow(2) + y.pow(2) > rad.pow(2) {
                    continue;
                }
                let px = (x + xc) as usize;
                let py = (y + yc) as usize;
                if px >= self.shape[1] || py >= self.shape[0] {
                    continue;
                }
                pixels.push(Pixel{
                    x: px,
                    y: py,
                    val: self.data[py*self.shape[1]+px]
                });
            }
        }
        pixels
    }
}
    

/// coordinate struct for interfacing with coordinate files
#[derive(Debug,Clone)]
#[pyclass]
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
                Err(_) => return Err(MavDACError::Coordinate(
                    format!("failed to parse x-ordinate: {}", a)
                ))
            };
        } else {
            return Err(MavDACError::Coordinate(
                "missing x ordinate".to_string(),
            ))
        }
        if let Some(a) = split.next() {
            y = match a.parse::<f64>() {
                Ok(y) => y,
                Err(_) => return Err(MavDACError::Coordinate(
                    format!("failed to parse y-ordinate: {}", a)
                ))
            };
        } else {
            return Err(MavDACError::Coordinate(
                "missing y ordinate".to_string(),
            ))
        }
        Ok(
            Self {
                pos: Vec2D{x,y},
            }
        )
    }
    
    type Error = MavDACError;
}

impl Display for Coordinate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{},{},", self.pos.x, self.pos.y)?;
        Ok(())
    }
}