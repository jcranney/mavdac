use pyo3::pyclass;
use pyo3::pymethods;
use crate::Vec2D;

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
        (0..self.get_coeffs().len()).map(|index|
            {
                self.get_coeffs()[index] * self.sample(pos, index)
            }
        ).fold(Vec2D{x:0.0,y:0.0}, |a,b| a+b)
    }
}

/// Bivariate Homogenous Polynomial to be used as distortion basis function
/// 
/// See [wikipedia article](https://en.wikipedia.org/wiki/Homogeneous_polynomial) for
/// more info.
#[pyclass]
pub struct BiVarPolyDistortions{
    /// degree of polynomial (maximum order)
    pub degree: usize,
    /// coefficients of distortions
    pub coeffs: Vec<Vec2D>,
    /// shape of image (numpy format)
    pub shape: [usize; 2],
    nk_lut: Vec<(usize,usize)>,
}

#[pymethods]
impl BiVarPolyDistortions {
    /// construct a new set of bivariate homogenous polynomials
    #[new]
    pub fn new(degree: usize, shape: [usize; 2]) -> Self {
        let coeffs = vec![Vec2D{x:0.0,y:0.0}; ((degree+1)*(degree+2))/2-1];
        let ncoeffs = coeffs.len();
        Self {
            degree,
            coeffs,
            shape,
            nk_lut: (0..ncoeffs).map(Self::index_to_nk).collect()
        }
    }

    #[staticmethod]
    fn index_to_nk(index: usize) -> (usize, usize) {
        let n = ((-1.0 + (1.0+8.0*((index+1) as f64)).sqrt())/2.0).floor() as usize;
        let k = index + 1 - n*(n+1)/2;
        (n,k)
    }

    /// sample basis function given index at x/y coordinates
    pub fn sample_xy(&self, x: f64, y: f64, ell: usize) -> f64 {
        self.sample(&Vec2D{x,y}, ell)
    }

    /// evaluate distortions (including coefficients) at x/y coordinates
    pub fn eval_xy(&self, x: f64, y: f64) -> (f64,f64) {
        let Vec2D{x,y} = self.eval(&Vec2D{x,y});
        (x,y)
    }

    #[getter]
    fn ncoeffs(&self) -> usize {
        self.coeffs.len()
    }

    #[getter]
    fn coeffs(&self) -> Vec<Vec<f64>> {
        self.coeffs.clone().into_iter()
        .map(|v| {
            let Vec2D{x,y} = v;
            vec![x,y]
        }).collect::<Vec<Vec<f64>>>()
    }

    /// load coefficients (e.g.) from python
    pub fn load_coeffs(&mut self, coeffs: Vec<Vec<f64>>) {
        self.set_coeffs(
            coeffs.into_iter().map(|p|
                Vec2D{x: p[0], y: p[1]}
            ).collect()
        );
    }
}

impl DistortionBasis for BiVarPolyDistortions {
    fn sample(&self, pos: &Vec2D, index: usize) -> f64 {
        let (n,k) = self.nk_lut[index];
        let Vec2D{mut x, mut y} = pos;
        x -= (self.shape[1] as f64)/2.0;
        y -= (self.shape[0] as f64)/2.0;
        x /= self.shape[1] as f64;
        y /= self.shape[0] as f64;
        x.powf(k as f64)*y.powf((n-k) as f64)
    }
    
    fn get_coeffs(&self) -> &Vec<Vec2D> {
        &self.coeffs
    }
    
    fn set_coeffs(&mut self, coeffs: Vec<Vec2D>) {
        self.coeffs = coeffs;
    }
}
