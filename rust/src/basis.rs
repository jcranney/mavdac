use ndarray::Array2;
use ndarray_inverse::Inverse;

use crate::MavDACError;
use crate::Vec2D;
use crate::Centroid;

pub trait DistortionBasis {
    /// the sample function returns the value of the bivariate function with a given
    /// index, at coordinate (x,y)
    fn sample(&self, pos: &Vec2D, index: usize) -> f64;
    /// a method to get the coefficients vector
    fn get_coeffs(&self) -> &Vec<Vec2D>;
    /// a method to set the coefficients vector
    fn set_coeffs(&mut self, coeffs: Vec<Vec2D>);


    fn eval(&self, pos: &Vec2D) -> Vec2D {
        (0..self.get_coeffs().len()).map(|index|
            {
                self.get_coeffs()[index].clone() * self.sample(pos, index)
            }
        ).fold(Vec2D{x:0.0,y:0.0}, |a,b| a+b)
    }

    fn solve_for_coeffs(&mut self, cogs: Vec<Centroid>) -> Result<(), MavDACError> {
        let mut b: Vec<Vec<f64>> = vec![];
        for cog in &cogs {
            b.push(vec![
                cog.cogx - cog.x,
                cog.cogy - cog.y,
            ]);
        }
        let b_matrix: Array2<f64> = Array2::from_shape_vec(
            (cogs.len(), 2),
            b.into_iter().flatten().collect()
        ).unwrap().t().to_owned();
        
        let mut f: Vec<Vec<f64>> = vec![];
        for cog in &cogs {
            f.push(
                (0..self.get_coeffs().len()).map(|index|
                    self.sample(
                        &Vec2D {
                            x: cog.x,
                            y: cog.y,
                        },
                        index
                    )
                ).collect()
            );
        }
        let f_matrix: Array2<f64> = Array2::from_shape_vec(
            (cogs.len(), self.get_coeffs().len()),
            f.into_iter().flatten().collect()
        ).unwrap().t().to_owned();

        // 1. Compute F.T (transpose of F)
        let f_t = f_matrix.t();
        // 2. Compute F @ F.T
        let f_f_t: Array2<f64> = f_matrix.dot(&f_t);
        // 3. Invert (F @ F.T)
        let f_f_t_inv = match f_f_t.inv() {
            Some(inv) => inv,
            None => return Err(MavDACError::LinalgError("Gram matrix of sampled basis functions is not invertible".to_string())),
        };
        // 4. Compute result = B @ F.T @ inv(F @ F.T)
        let result = b_matrix.dot(&f_t).dot(&f_f_t_inv);
        
        self.set_coeffs(result.t().rows().into_iter().map(
            |row| Vec2D{x:row[0],y:row[1]}
        ).collect());        
        Ok(())
    }
}

pub struct BiVarPolyDistortions{
    pub degree: usize,
    pub coeffs: Vec<Vec2D>,
}

impl BiVarPolyDistortions {
    pub fn new(degree: usize) -> Self {
        Self {
            degree,
            coeffs: vec![Vec2D{x:0.0,y:0.0}; ((degree+1)*(degree+2))/2],
        }
    }
}

impl DistortionBasis for BiVarPolyDistortions {
    fn sample(&self, pos: &Vec2D, index: usize) -> f64 {
        let n = ((-1.0 + (1.0+8.0*(index as f64)).sqrt())/2.0).floor() as usize;
        let k = index - n*(n+1)/2;
        let mut pos = pos.clone();
        pos += Vec2D{x:-2048.0, y:-2048.0};
        pos.x /= 2048.0;
        pos.y /= 2048.0;
        pos.x.powf(k as f64)*pos.y.powf((n-k) as f64)
    }
    
    fn get_coeffs(&self) -> &Vec<Vec2D> {
        &self.coeffs
    }
    
    fn set_coeffs(&mut self, coeffs: Vec<Vec2D>) {
        self.coeffs = coeffs;
    }
}
