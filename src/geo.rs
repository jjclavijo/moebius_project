use num_complex::Complex;
use num_traits::Float;
use std::f64::consts::PI;

pub const DEG2RAD: f64 = PI / 180.;

use geo::{Coord, Point};
use serde::Serialize;

// Define a type to represent geographic points
// pero deberia ser un Trait que se implementa para Point
#[derive(Debug, Clone, Copy, Serialize)]
pub struct GeographicPoint  {
    //pub latitude: f64,
    //pub longitude: f64,
    #[serde(serialize_with = "geojson::ser::serialize_geometry")]
    pub geometry: Point<f64>
}

impl GeographicPoint 
    {
    pub fn new(latitude: f64, longitude: f64) -> GeographicPoint {
        GeographicPoint { geometry: Point::new(longitude,latitude) }
    }
    pub fn n_stereographic( self: GeographicPoint ) -> Complex<f64> {
        let Point (Coord {x: longitude, y: latitude}) = self.geometry;
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
        let Point (Coord {x: longitude, y: latitude}) = self.geometry;
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
    
    pub fn latitude( self: GeographicPoint ) -> f64 {
        let Point (Coord {x: _, y: latitude}) = self.geometry;
        latitude
    }
    pub fn longitude( self: GeographicPoint ) -> f64 {
        let Point (Coord {x: longitude, y:_}) = self.geometry;
        longitude
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

       assert!((p.latitude() - pp.latitude()).abs() < 1.0e-6, "{} is not approximately equal to {}", p.latitude(), pp.latitude())

    }

    #[test]
    fn test_s_stereographic()
    {
       use super::GeographicPoint;

       let p = GeographicPoint::new(-30.0,40.0);

       let pp = GeographicPoint::from_s_stereographic(p.s_stereographic());

       assert!((p.latitude() - pp.latitude()).abs() < 1.0e-6, "{} is not approximately equal to {}", p.latitude(), pp.latitude())
    }
}

use rand::Rng;
use rand::distributions::{Uniform,Distribution};

// a < lat < b
// c < lon < d
// Sampleo con uniforme en el area al menos.
pub fn sample_rectangulo<R:Rng + ?Sized>(r: &mut R, a: f64, b: f64, c: f64, d:f64, s:usize) -> Vec<GeographicPoint>
    {
        let a_ = a * DEG2RAD;
        let b_ = b * DEG2RAD;
        let c_ = c * DEG2RAD;
        let d_ = d * DEG2RAD;

        let max_sen = f64::max(a_.sin(), b_.sin());
        let min_sen = f64::min(a_.sin(), b_.sin());

        let n = Uniform::new(a_,b_);
        let e = Uniform::new(c_ * min_sen,d_ * max_sen);
        
        //vec![GeographicPoint::new(x.sample(r),y.sample(r))]
        std::iter::repeat(()).map(|()| {
            let norte = n.sample(r);
            println!("{}",norte);

            let mut este = e.sample(r) / norte.sin();
            while este < c_ || este > d_
            {
                este = e.sample(r) / norte.sin();
                println!("{}",este);
            };
            GeographicPoint::new(este/DEG2RAD,norte/DEG2RAD)
        }).take(s).collect()
    }

