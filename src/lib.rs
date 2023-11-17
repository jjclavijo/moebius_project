use num_complex::Complex;
use rand_distr::{Normal};
use rand::{Rng,SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::f64::consts::PI;

const DEG2RAD: f64 = PI / 180.

// Define a type to represent geographic points
#[derive(Debug, Clone, Copy)]
struct GeographicPoint {
    latitude: f64,
    longitude: f64,
}

impl GeographicPoint {
    fn n_stereographic( self: GeographicPoint ) -> Complex<f64> {
        let (latitude, longitude) = self
        let rcolat = DEG2RAD * (90 - latitude)
        let rlon = DEG2RAD * longitude
        let x = rcolat.sin() * rlon.cos();
        let y = rcolat.sin() * rlon.sin();
        let z = rcolat.cos();

        let X = x / (1-z);
        let Y = y / (1-z);

        Complex<f64>::new(X,Y)
    } 

    fn s_stereographic( self: GeographicPoint ) -> Complex<f64> {
        let (latitude, longitude) = self
        let rcolat = DEG2RAD * (90 + latitude)
        let rlon = DEG2RAD * longitude
        let x = rcolat.sin() * rlon.cos();
        let y = rcolat.sin() * rlon.sin();
        let z = rcolat.cos();

        let X = x / (1-z);
        let Y = y / (1-z);

        Complex<f64>::new(X,Y)
    } 

    fn from_s_stereographic( point: Complex<f64>  ) -> GeographicPoint {
        let X = point.re
        let Y = point.im
        let x = 2X / ( 1 + X.pow(2) + Y.pow(2) )
        let y = 2Y / ( 1 + X.pow(2) + Y.pow(2) )
        let z = (X.pow(2) + Y.pow(2) - 1) / ( 1 + X.pow(2) + Y.pow(2) )

        let rcolat = z.acos();
        let rlon = (y / rcolat.sin()).atan2(x / rcolat.sin());

        GeographicPoint { (90 - rcolat) / DEG2RAD , rlon / DEG2RAD }
    } 

    fn from_n_stereographic( point: Complex<f64>  ) -> GeographicPoint {
        let X = point.re
        let Y = point.im
        let x = 2X / ( 1 + X.pow(2) + Y.pow(2) )
        let y = 2Y / ( 1 + X.pow(2) + Y.pow(2) )
        let z = (X.pow(2) + Y.pow(2) - 1) / ( 1 + X.pow(2) + Y.pow(2) )

        let rcolat = z.acos();
        let rlon = (y / rcolat.sin()).atan2(x / rcolat.sin());

        GeographicPoint { (rcolat - 90) / DEG2RAD , rlon / DEG2RAD }
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
    fn from_pole(latitude: f64, longitude: f64, k: f64) -> MobiusTransformation {
        let point1 = GeographicPoint { latitude, longitude };
        let gamma1 = point1::n_stereographic();
        let gamma2 = compute_antipodal_point(point1)::n_stereographic();

        MobiusTransformation::from_fixed_points(gamma1, gamma2, k)
    }

    fn from_fixed_points(gamma1: Complex<f64>, gamma2: Complex<f64>, k: f64) -> MobiusTransformation {
        let a = gamma1 - k * gamma2;
        let b = (k - 1.0) * gamma1 * gamma2;
        let c = Complex::new(1.0 - k, 0.0);
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
        let det_a = determinant(z1 * w1, w1, 1.0, z2 * w2, w2, 1.0, z3 * w3, w3, 1.0);
        let det_b = determinant(z1 * w1, z1, w1, z2 * w2, z2, w2, z3 * w3, z3, w3);
        let det_c = determinant(z1, w1, 1.0, z2, w2, 1.0, z3, w3, 1.0);
        let det_d = determinant(z1 * w1, z1, 1.0, z2 * w2, z2, 1.0, z3 * w3, z3, 1.0);

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

// Define a Möbius aligned transformation struct
#[derive(Debug, Clone, Copy)]
struct MobiusAlignedTransformation {
    base_transform: MobiusTransformation,
    noise_scale: MobiusTransformation,
}

impl MobiusAlignedTransformation {
    fn apply(&self, z: Complex<f64>) -> Complex<f64> {
        let perturbed_transform = self.base_transform.clone() + self.noise_scale;

        perturbed_transform.apply(z)
    }
}

// Function to compute the antipodal point
fn compute_antipodal_point(point: GeographicPoint) -> GeographicPoint {
    let (latitude, longitude) = point
    let rcolat = DEG2RAD * (90 - latitude)
    let rlon = DEG2RAD * longitude
    let x = - rcolat.sin() * rlon.cos();
    let y = - rcolat.sin() * rlon.sin();
    let z = - rcolat.cos();

    let rcolat = z.acos();
    let rlon = (y / rcolat.sin()).atan2(x / rcolat.sin());
    
    GeographicPoint { (90 - rcolat) / DEG2RAD , rlon / DEG2RAD }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mobius_transformation() {
        let base_point = GeographicPoint {
            latitude: 45.0,
            longitude: 30.0,
        };

        let mobius_transformation = MobiusTransformation::from_pole(base_point.latitude,base_point.longitude, 1.0);

        // Assuming the perturbations are applied correctly, test a transformation's output
        let complex_input = Complex::new(3.0, 4.0);
        let transformed_output = mobius_transformations[0].apply(complex_input);

        // Example assertion: Check if the output is within an acceptable range or condition
        assert!((transformed_output.re).abs() < 1e6);
        assert!((transformed_output.im).abs() < 1e6);
    }
}
