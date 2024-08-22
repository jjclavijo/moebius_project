use crate::geo::GeographicPoint;
use crate::geo::DEG2RAD;
use geo::Coord;
use geo::Point;
use num_complex::Complex;
use num_complex::ComplexFloat;
use num_traits::One;
use num_traits::Zero;

// Define a Möbius transformation struct
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MobiusTransformation {
    a: Complex<f64>,
    b: Complex<f64>,
    c: Complex<f64>,
    d: Complex<f64>,
}

impl MobiusTransformation {
    pub fn from_pole(latitude: f64, longitude: f64, a: f64) -> MobiusTransformation {

        let k = Complex::<f64>::from_polar(1.0, a);

        let point1 = GeographicPoint::new( latitude, longitude );
        let gamma1 = point1.n_stereographic();
        let gamma2 = compute_antipodal_point(point1).n_stereographic();

        MobiusTransformation::from_fixed_points(gamma1, gamma2, k)
    }

    pub fn from_fixed_points(gamma1: Complex<f64>, gamma2: Complex<f64>, k: Complex<f64>) -> MobiusTransformation {
        let a = gamma1 - k * gamma2;
        let b = (k - 1.0) * gamma1 * gamma2;
        let c = 1.0 - k;
        let d = k * gamma1 - gamma2;

        MobiusTransformation { a, b, c, d }
    }

    pub fn from_pairs(
        z1: Complex<f64>,
        z2: Complex<f64>,
        z3: Complex<f64>,
        w1: Complex<f64>,
        w2: Complex<f64>,
        w3: Complex<f64>,
    ) -> MobiusTransformation {

        // Esta cuenta es la misma que la del determinante.
        // let t = (w1 - w2) * (z1 - z3);
        // let k = (w1 - w3) * (z1 - z2);

        // let a = t * w3 - k * w2;
        // let b = k * w2 * z3 - t * w3 * z2;
        // let c = t - k;
        // let d = k * z3 - t * z2;

        let a = determinant(z1 * w1, w1, 1.0.into(),
                            z2 * w2, w2, 1.0.into(), 
                            z3 * w3, w3, 1.0.into());

        let b = determinant(z1 * w1, z1, w1, 
                            z2 * w2, z2, w2, 
                            z3 * w3, z3, w3);

        let c = determinant(z1, w1, 1.0.into(),
                            z2, w2, 1.0.into(),
                            z3, w3, 1.0.into());

        let d = determinant(z1 * w1, z1, 1.0.into(),
                            z2 * w2, z2, 1.0.into(),
                            z3 * w3, z3, 1.0.into());

        MobiusTransformation { a, b, c, d }
    }

    pub fn fixed_points(&self) -> (Complex<f64>,Complex<f64>) {
        if self.c == Complex::zero() {
            let p1 = Complex::one() * f64::INFINITY;
            if self.a == self.d {
                (p1.clone(), p1)
            } else {
                (p1, - self.b / (self.a - self.d))
            }
        } else {
            let det = self.a * self.d - self.b * self.c;
            let tr = self.a + self.d;
            let discr = tr * tr - 4. * det;

            let base_p = (self.a - self.d) / (2. * self.c);

            if discr == Complex::zero() {
                (base_p.clone(), base_p)
            } else {
                let addit = discr.sqrt() / (2. * self.c);
                (base_p + addit, base_p - addit)
            }
        }
    }

    pub fn apply(&self, z: Complex<f64>) -> Complex<f64> {
        let numerator = self.a * z + self.b;
        let denominator = self.c * z + self.d;

        numerator / denominator
    }
    
    pub fn add(self, other: MobiusTransformation) -> MobiusTransformation {
        let a = self.a + other.a;
        let b = self.b + other.b;
        let c = self.c + other.c;
        let d = self.d + other.d;
        MobiusTransformation { a, b, c, d }
    }

    pub fn normalizar(self) -> MobiusTransformation {
        let n = self.a;

        MobiusTransformation::new( self.a / n, self.b / n,
                                   self.c / n, self.d / n )

    }

    pub fn new(a: Complex<f64>, b: Complex<f64>, 
           c: Complex<f64>, d: Complex<f64>) -> MobiusTransformation {
        MobiusTransformation { a, b, c, d }
    }

}


fn determinant(a: Complex<f64>, b: Complex<f64>, c: Complex<f64>, d: Complex<f64>, e: Complex<f64>, f: Complex<f64>, g: Complex<f64>, h: Complex<f64>, i: Complex<f64>) -> Complex<f64> {
    // a b c
    // d e f
    // g h i
    
    //
    // Determinante con las matrices menores.
    //
   
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
    let Point (Coord {y:latitude, x:longitude}) = point.geometry;
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

    #[test]
    fn test_antipodal_transformation() {

        let base_point = GeographicPoint {
            geometry: Point ( Coord{
            y: 45.0,
            x: 30.0
            })
        };

        let mobius_transformation = MobiusTransformation::from_pole(base_point.latitude(),base_point.longitude(), 1.0);

        let antipodal_point = compute_antipodal_point(base_point);

        let mobius_transformation_a = MobiusTransformation::from_pole(antipodal_point.latitude(),antipodal_point.longitude(), -1.0);

        // Assuming the perturbations are applied correctly, test a transformation's output
        //
        let complex_input = Complex::new(3.0, 4.0);

        let transformed_output = mobius_transformation.apply(complex_input);

        //
        let transformed_output_a = mobius_transformation_a.apply(complex_input);


        println!("--- INICIO: Debug Antipoda y Contragiro ---");
        println!("|");

        println!("{:?}",complex_input);
        println!("{:?}",base_point);
        println!("{:?}",antipodal_point);

        println!("{:?}",mobius_transformation);
        println!("{:?}",mobius_transformation_a);

        println!("{:?}",transformed_output);
        println!("{:?}",transformed_output_a);

        println!("|");
        println!("--- FIN: Debug Antipoda y Contragiro ---");

        // Example assertion: Check if the output is within an acceptable range or condition
        assert!(( transformed_output - transformed_output_a ).norm() < 1e-6);
    }

    #[test]
    fn test_polo_vs_tres_puntos() {
        let base_point = GeographicPoint {
            geometry: Point ( Coord{
            y: 45.0,
            x: 30.0
            })
        };

        let mobius_transformation = MobiusTransformation::from_pole(base_point.latitude(),base_point.longitude(), 1.0);

        let complex_inputs : Vec<Complex<f64>> = vec![Complex::new(1.0, -4.0),
                                                    Complex::new(3.0, 4.0),
                                                    Complex::new(-3.0, 8.0)];

        let transformed_output : Vec<Complex<f64>> = complex_inputs.iter().map(|&i| mobius_transformation.apply(i)).collect();

        let a: Complex<f64> = Complex::new(f64::NAN,f64::NAN);
        let b: Complex<f64> = Complex::new(f64::NAN,f64::NAN);
        let c: Complex<f64> = Complex::new(f64::NAN,f64::NAN);
        let d: Complex<f64> = Complex::new(f64::NAN,f64::NAN);

        let mut other_mobius_transformation = MobiusTransformation { a, b, c, d };

        if let [z1,z2,z3] = complex_inputs[..] {
            if let [w1,w2,w3] = transformed_output[..] {
                other_mobius_transformation = MobiusTransformation::from_pairs(z1,z2,z3,w1,w2,w3);
            }
        }

        println!("--- INICIO: Debug Polos vs Tres Puntos ---");
        println!("|");
        println!("{:?}",other_mobius_transformation.normalizar());
        println!("{:?}",mobius_transformation.normalizar());
        println!("|");
        println!("--- FIN: Debug Polos vs Tres Puntos ---");

        assert!( (other_mobius_transformation.normalizar().a 
                      - mobius_transformation.normalizar().a).norm() < 1e-6);

        assert!( (other_mobius_transformation.normalizar().b 
                      - mobius_transformation.normalizar().b).norm() < 1e-6);

        assert!( (other_mobius_transformation.normalizar().c 
                      - mobius_transformation.normalizar().c).norm() < 1e-6);

        assert!( (other_mobius_transformation.normalizar().d 
                      - mobius_transformation.normalizar().d).norm() < 1e-6);

    }

    #[test]
    fn test_polo_y_vuelta() {
        let base_point = GeographicPoint {
            geometry: Point ( Coord{
            y: 45.0,
            x: 30.0
            })
        };

        let mobius_transformation = MobiusTransformation::from_pole(base_point.latitude(),base_point.longitude(), 1.0);

        let (f1,f2) = mobius_transformation.fixed_points();

        let f1g = GeographicPoint::from_n_stereographic(f1);
        let f2g = GeographicPoint::from_n_stereographic(f2);

        println!("{:?}, {:?}", f1g, f2g);

        assert!( (f1g.dist(base_point) < 1e12) | (f2g.dist(base_point) < 1e12));
    }
}
