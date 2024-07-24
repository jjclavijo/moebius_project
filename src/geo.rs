use num_complex::Complex;
use num_traits::Float;
use std::f64::consts::PI;

pub const DEG2RAD: f64 = PI / 180.;

// Define a type to represent geographic points
#[derive(Debug, Clone, Copy)]
pub struct GeographicPoint  {
    pub latitude: f64,
    pub longitude: f64,
}

impl GeographicPoint 
    {
    pub fn new(latitude: f64, longitude: f64) -> GeographicPoint {
        GeographicPoint { latitude,longitude }
    }
    pub fn n_stereographic( self: GeographicPoint ) -> Complex<f64> {
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

    pub fn s_stereographic( self: GeographicPoint ) -> Complex<f64> {
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

    pub fn from_s_stereographic( point: Complex<f64>  ) -> GeographicPoint {
        let s_x = point.re;
        let s_y = point.im;
        let x = (2. * s_x) / ( 1.0 + s_x.powi(2) + s_y.powi(2) );
        let y = (2. * s_y) / ( 1.0 + s_x.powi(2) + s_y.powi(2) );
        let z = (s_x.powi(2) + s_y.powi(2) - 1.0) / ( 1.0 + s_x.powi(2) + s_y.powi(2) );

        let rcolat = z.acos();
        let rlon = (y / rcolat.sin()).atan2(x / rcolat.sin());

        // GeographicPoint::new( (90.0 - rcolat) / DEG2RAD , rlon / DEG2RAD )
        GeographicPoint::new( rcolat / DEG2RAD - 90.0 , rlon / DEG2RAD )
    } 

    pub fn from_n_stereographic( point: Complex<f64>  ) -> GeographicPoint {
        let s_x = point.re;
        let s_y = point.im;
        let x = (2. * s_x) / ( 1.0 + s_x.powi(2) + s_y.powi(2) );
        let y = (2. * s_y) / ( 1.0 + s_x.powi(2) + s_y.powi(2) );
        let z = (s_x.powi(2) + s_y.powi(2) - 1.0) / ( 1.0 + s_x.powi(2) + s_y.powi(2) );

        let rcolat = z.acos();
        let rlon = (y / rcolat.sin()).atan2(x / rcolat.sin());

        GeographicPoint::new( 90.0 - rcolat / DEG2RAD , rlon / DEG2RAD )
    } 
}

mod test {
    // Parametrizar tests
    #[test]
    fn test_n_stereographic()
    {
       use super::GeographicPoint;

       let p = GeographicPoint::new(30.0,40.0);

       let pp = GeographicPoint::from_n_stereographic(p.n_stereographic());

       assert!((p.latitude - pp.latitude).abs() < 1.0e-6, "{} is not approximately equal to {}", p.latitude, pp.latitude)

    }

    #[test]
    fn test_s_stereographic()
    {
       use super::GeographicPoint;

       let p = GeographicPoint::new(-30.0,40.0);

       let pp = GeographicPoint::from_s_stereographic(p.s_stereographic());

       assert!((p.latitude - pp.latitude).abs() < 1.0e-6, "{} is not approximately equal to {}", p.latitude, pp.latitude)

    }
}
