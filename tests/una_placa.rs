use mobius_project::mobius::MobiusTransformation;
use mobius_project::geo::GeographicPoint;


use rand::Rng;
use rand::distributions::Uniform;
use std::f64::consts::PI; 

fn sample_rectangulo<R:Rng + ?Sized>(r: &mut R,l: f64, a: f64, b: f64, c: f64, d:f64) -> Vec<GeographicPoint>
    {
        let x = Uniform::new(a,b);
        let y = Uniform::new(c,d);
        
        vec![GeographicPoint::new(x.sample(r),y.sample(r))]

    }

#[test]
fn test_placa() {

    // Simulación de una placa

    let transform = MobiusTransformation::from_pole(30.0,40.0,0.01);  


    // Se generan puntos
    let puntos: Vec<GeographicPoint> = vec![GeographicPoint::new(31.0,42.0)];

    // Se les aplica una transformación con polo en 30 40.
    let nuevos_puntos: Vec<GeographicPoint> = 
        puntos.iter().map(|p| p.s_stereographic())
                     .map(|z| transform.apply(z))
                     .map(|z| GeographicPoint::from_s_stereographic(z)).collect();

    // Se les agrega ruido aleatorio
    //

    // Se sortean ternas de puntos
    //

    // Se calculan las transformaciones correspondientes a cada.
    //

    // RANSAC
    //

    println!("{:?}",puntos);
    println!("{:?}",nuevos_puntos);

    assert!(true)
}
