use num_complex::Complex;
use std::f64::consts::PI;

const DEG2RAD: f64 = PI / 180.;

// Define a type to represent geographic points
#[derive(Debug, Clone, Copy)]
struct GeographicPoint {
    latitude: f64,
    longitude: f64,
}

impl GeographicPoint {
    fn new(latitude: f64, longitude: f64) -> GeographicPoint {
        GeographicPoint { latitude,longitude }
    }
    fn n_stereographic( self: GeographicPoint ) -> Complex<f64> {
        let GeographicPoint {latitude, longitude} = self;
        let rcolat = DEG2RAD * (90.0 - latitude);
        let rlon = DEG2RAD * longitude;
        let x = rcolat.sin() * rlon.cos();
        let y = rcolat.sin() * rlon.sin();
        let z = rcolat.cos();

        let s_x = x / (1.0 - z);
        let s_y = y / (1.0 - z);

        Complex::<f64>::new(s_x,s_y)
    } 

    fn s_stereographic( self: GeographicPoint ) -> Complex<f64> {
        let GeographicPoint {latitude, longitude} = self;
        let rcolat = DEG2RAD * (90.0 + latitude);
        let rlon = DEG2RAD * longitude;
        let x = rcolat.sin() * rlon.cos();
        let y = rcolat.sin() * rlon.sin();
        let z = rcolat.cos();

        let s_x = x / (1.0 - z);
        let s_y = y / (1.0 - z);

        Complex::<f64>::new(s_x,s_y)
    } 

    fn from_s_stereographic( point: Complex<f64>  ) -> GeographicPoint {
        let s_x = point.re;
        let s_y = point.im;
        let x = (2. * s_x) / ( 1.0 + s_x.powi(2) + s_y.powi(2) );
        let y = (2. * s_y) / ( 1.0 + s_x.powi(2) + s_y.powi(2) );
        let z = (s_x.powi(2) + s_y.powi(2) - 1.0) / ( 1.0 + s_x.powi(2) + s_y.powi(2) );

        let rcolat = z.acos();
        let rlon = (y / rcolat.sin()).atan2(x / rcolat.sin());

        GeographicPoint::new( (90.0 - rcolat) / DEG2RAD , rlon / DEG2RAD )
    } 

    fn from_n_stereographic( point: Complex<f64>  ) -> GeographicPoint {
        let s_x = point.re;
        let s_y = point.im;
        let x = (2. * s_x) / ( 1.0 + s_x.powi(2) + s_y.powi(2) );
        let y = (2. * s_y) / ( 1.0 + s_x.powi(2) + s_y.powi(2) );
        let z = (s_x.powi(2) + s_y.powi(2) - 1.0) / ( 1.0 + s_x.powi(2) + s_y.powi(2) );

        let rcolat = z.acos();
        let rlon = (y / rcolat.sin()).atan2(x / rcolat.sin());

        GeographicPoint::new( (rcolat - 90.0) / DEG2RAD , rlon / DEG2RAD )
    } 
}

// Define a Möbius transformation struct
#[derive(Debug, Clone, Copy)]
struct MobiusTransformation {
    a: Complex<f64>,
    b: Complex<f64>,
    c: Complex<f64>,
    d: Complex<f64>,
}

impl MobiusTransformation {
    fn from_pole(latitude: f64, longitude: f64, a: f64) -> MobiusTransformation {

        let k = Complex::<f64>::from_polar(1.0, a);

        let point1 = GeographicPoint { latitude, longitude };
        let gamma1 = point1.n_stereographic();
        let gamma2 = compute_antipodal_point(point1).n_stereographic();

        MobiusTransformation::from_fixed_points(gamma1, gamma2, k)
    }

    fn from_fixed_points(gamma1: Complex<f64>, gamma2: Complex<f64>, k: Complex<f64>) -> MobiusTransformation {
        let a = gamma1 - k * gamma2;
        let b = (k - 1.0) * gamma1 * gamma2;
        let c = 1.0 - k;
        let d = k * gamma1 - gamma2;

        MobiusTransformation { a, b, c, d }
    }

    fn from_pairs(
        z1: Complex<f64>,
        z2: Complex<f64>,
        z3: Complex<f64>,
        w1: Complex<f64>,
        w2: Complex<f64>,
        w3: Complex<f64>,
    ) -> MobiusTransformation {
        let det_a = determinant(z1 * w1, w1, 1.0.into(), z2 * w2, w2, 1.0.into(), z3 * w3, w3, 1.0.into());
        let det_b = determinant(z1 * w1, z1, w1, z2 * w2, z2, w2, z3 * w3, z3, w3);
        let det_c = determinant(z1, w1, 1.0.into(), z2, w2, 1.0.into(), z3, w3, 1.0.into());
        let det_d = determinant(z1 * w1, z1, 1.0.into(), z2 * w2, z2, 1.0.into(), z3 * w3, z3, 1.0.into());

        let a = det_a / det_c;
        let b = det_b / det_c;
        let c = det_c / det_c;
        let d = det_d / det_c;

        MobiusTransformation { a, b, c, d }
    }

    fn apply(&self, z: Complex<f64>) -> Complex<f64> {
        let numerator = self.a * z + self.b;
        let denominator = self.c * z + self.d;

        numerator / denominator
    }
    
    fn add(self, other: MobiusTransformation) -> MobiusTransformation {
        let a = self.a + other.a;
        let b = self.b + other.b;
        let c = self.c + other.c;
        let d = self.d + other.d;
        MobiusTransformation { a, b, c, d }
    }
    
}


fn determinant(a: Complex<f64>, b: Complex<f64>, c: Complex<f64>, d: Complex<f64>, e: Complex<f64>, f: Complex<f64>, g: Complex<f64>, h: Complex<f64>, i: Complex<f64>) -> Complex<f64> {
    a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
}

// Implement Add trait for MobiusTransformation
use std::ops::Add;
impl Add for MobiusTransformation {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        self.add(other)
    }
}

// Function to compute the antipodal point
fn compute_antipodal_point(point: GeographicPoint) -> GeographicPoint {
    let GeographicPoint {latitude, longitude} = point;
    let rcolat = DEG2RAD * (90.0 - latitude);
    let rlon = DEG2RAD * longitude;
    let x = - rcolat.sin() * rlon.cos();
    let y = - rcolat.sin() * rlon.sin();
    let z = - rcolat.cos();

    let rcolat = z.acos();
    let rlon = (y / rcolat.sin()).atan2(x / rcolat.sin());
    
    GeographicPoint::new( 90.0 - rcolat / DEG2RAD , rlon / DEG2RAD )
}

#[cfg(test)]
mod tests {
    use super::*;
    // use rand_distr::{Normal};
    // use rand::{Rng,SeedableRng};
    // use rand_chacha::ChaCha8Rng;

    #[test]
    fn test_mobius_transformation() {
        let base_point = GeographicPoint {
            latitude: 45.0,
            longitude: 30.0,
        };

        let mobius_transformation = MobiusTransformation::from_pole(base_point.latitude,base_point.longitude, 1.0);

        // Assuming the perturbations are applied correctly, test a transformation's output
        let complex_input = Complex::new(3.0, 4.0);
        let transformed_output = mobius_transformation.apply(complex_input);

        // Example assertion: Check if the output is within an acceptable range or condition
        assert!((transformed_output.re).abs() < 1e6);
        assert!((transformed_output.im).abs() < 1e6);
    }

    #[test]
    fn test_antipodal_transformation() {
        let base_point = GeographicPoint {
            latitude: 45.0,
            longitude: 30.0,
        };

        let mobius_transformation = MobiusTransformation::from_pole(base_point.latitude,base_point.longitude, 1.0);

        let antipodal_point = compute_antipodal_point(base_point);

        let mobius_transformation_a = MobiusTransformation::from_pole(antipodal_point.latitude,antipodal_point.longitude, -1.0);

        // Assuming the perturbations are applied correctly, test a transformation's output
        //
        let complex_input = Complex::new(3.0, 4.0);

        let transformed_output = mobius_transformation.apply(complex_input);

        //
        let transformed_output_a = mobius_transformation_a.apply(complex_input);
        println!("{:?}",complex_input);
        println!("{:?}",base_point);
        println!("{:?}",antipodal_point);

        println!("{:?}",mobius_transformation);
        println!("{:?}",mobius_transformation_a);

        println!("{:?}",transformed_output);
        println!("{:?}",transformed_output_a);

        // Example assertion: Check if the output is within an acceptable range or condition
        assert!(( (transformed_output.re) - (transformed_output_a.re) ).abs() < 1e-6);
        assert!(( (transformed_output.im) - (transformed_output_a.im) ).abs() < 1e-6);
    }

}
