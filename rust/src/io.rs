use std::{f64::consts::PI, path::{Path, PathBuf}};

use fitrs::{Fits, FitsData, Hdu, HeaderValue};

use crate::{Centroid, Grid, MavDACError, Result, Vec2D};


#[derive(Debug, Clone)]
pub struct Pixel {
    x: usize,
    y: usize,
    val: f64,
}

#[derive(Clone,Debug)]
pub struct Image {
    pub data: Vec<f64>,
    pub shift: Option<Vec2D>,
    pub shape: Vec<usize>, // numpy style
}

impl std::fmt::Display for Image {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "shape: {:?}", self.shape)?;
        if let Some(shift) = &self.shift {
            write!(f, ", xshift {:0.4}, yshift {:0.4}", shift.x, shift.y)?
        }
        Ok(())
    }
}

impl Image {
    pub fn from_fits(filename: PathBuf) -> Result<Image> {
        let fits = Fits::open(&filename)?;
        if let Some(hdu) = fits.get(0) {
            match hdu.value("NAXIS")  {
                Some(HeaderValue::IntegerNumber(2)) => (),
                _ => return Err(MavDACError::InvalidFITS("expected NAXIS==2".to_string())),
            };
            let mut shape: Vec<usize> = vec![];
            match hdu.value("NAXIS2")  {
                Some(HeaderValue::IntegerNumber(x)) if *x > 0 => {
                    shape.push(*x as usize);
                },
                _ => return Err(MavDACError::InvalidFITS("invalid NAXIS2".to_string())),
            }
            match hdu.value("NAXIS1")  {
                Some(HeaderValue::IntegerNumber(x)) if *x > 0 => {
                    shape.push(*x as usize);
                },
                _ => return Err(MavDACError::InvalidFITS("invalid NAXIS1".to_string())),
            };
            let shift = {
                let xshift: f64 = match hdu.value("XSHIFT").ok_or(MavDACError::InvalidFITS(
                    format!("missing XSHIFT in fits header {}", filename.display())
                ))? {
                    HeaderValue::IntegerNumber(a) => *a as f64,
                    HeaderValue::RealFloatingNumber(a) => *a,
                    _ => return Err(MavDACError::InvalidFITS(
                        format!("XSHIFT in fits header has invalid datatype, \
                                must be float or int {}", filename.display())
                    ))
                };
                let yshift: f64 = match hdu.value("YSHIFT").ok_or(MavDACError::InvalidFITS(
                    format!("missing YSHIFT in fits header {}", filename.display())
                ))? {
                    HeaderValue::IntegerNumber(a) => *a as f64,
                    HeaderValue::RealFloatingNumber(a) => *a,
                    _ => return Err(MavDACError::InvalidFITS(
                        format!("YSHIFT in fits header has invalid datatype, \
                                must be float or int\n{}", filename.display())
                    ))
                };
                Some(Vec2D {
                    x: xshift,
                    y: yshift,
                })
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
                    format!("no primary hdu in {}", &filename.display())
                )
            )
        }
    }

    pub fn to_fits<P>(self, filename: P) -> Result<()> 
    where P: AsRef<Path> {
        let primary_hdu = Hdu::new(&self.shape, self.data);
        Fits::create(filename, primary_hdu)?;
        Ok(())
    }
    
    pub fn draw_on_circles(mut self, grid: &Grid, rad: f64, val: f64) -> Self {
        const NTHETA: usize = 1000;
        grid.all_points(&self)
        .into_iter()
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
        self
    }

    pub fn cogs(&self, grid: &Grid, rad: usize) -> Vec<Centroid> {
        // get all nominal positions
        let points = grid.all_points(&self);

        // measure cog and intensity within radius at all points
        points.into_iter().map(|point| (point.clone(), self.get_blob(&point, rad)))
        .map(|(point,pixels)|
            (
                point,
                pixels.iter()
                .map(|pixel| (
                    pixel.x as f64*pixel.val,
                    pixel.y as f64*pixel.val,
                    pixel.val
                )).fold((0.0,0.0,0.0), |a,b| (a.0+b.0, a.1+b.1, a.2+b.2))
            )
        ).map(|(point, (sumx,sumy,flux))| Centroid {
            cogx: sumx / flux,
            cogy: sumy / flux,
            flux,
            x: point.x,
            y: point.y,
        }).collect()
    }

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

#[derive(Debug,Clone)]
pub struct Coordinate {
    pub pos: Vec2D,
    pub dist: Option<Vec2D>,
}

impl Coordinate {
    pub fn from_str(s: &str) -> Result<Self> {
        println!("{}", s);
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
                dist: None,
            }
        )
    }
}
