use geo::{Coord, Point};
use mobius_project::mobius::MobiusTransformation;
use mobius_project::geo::{sample_rectangulo, GeographicPoint};

use num_complex::Complex;
use rand::thread_rng;
use rand::distributions::{Uniform,Distribution};

use itertools::Itertools;

#[test]
fn test_placa() {

    // Simulación de una placa

    let transform = MobiusTransformation::from_pole(30.0,40.0,0.1);  


    // Se generan puntos
    //let puntos: Vec<GeographicPoint> = vec![GeographicPoint::new(31.0,42.0)];
    let mut rngen = thread_rng();
    let puntos = sample_rectangulo(&mut rngen, 50., 80., 20., 30., 1000);

    // Se les aplica una transformación con polo en 30 40.
    let nuevos_puntos: Vec<GeographicPoint> = 
        puntos.iter().map(|p| p.s_stereographic())
                     .map(|z| transform.apply(z))
                     .map(|z| GeographicPoint::from_s_stereographic(z)).collect();

    // Se les agrega ruido aleatorio
    //

    let nuevos_puntos_r: Vec<GeographicPoint> = 
        nuevos_puntos.iter().map(
            |p| {
                let Point (Coord {x: longitude, y: latitude}) = p.geometry;
                let rdist = rand::distributions::Uniform::new(-1.,1.);
                let scale = 0.01;
                GeographicPoint::new(latitude+rdist.sample(&mut rngen)*scale,
                                     longitude+rdist.sample(&mut rngen)*scale)
            } ).collect();

    // Se sortean ternas de puntos
    //

    let uniforme_entera = Uniform::new(0,1000);
    let ternas: Vec<(usize,usize,usize)> = std::iter::repeat(()).map(|()| uniforme_entera.sample(&mut rngen)).tuples().filter_map(|(a,b,c)| if a != b && b != c {Some((a,b,c))} else {None}).take(3).collect();

    // Se calculan las transformaciones correspondientes a cada.
    //
    //
    let ts: Vec<MobiusTransformation> = ternas.clone().into_iter().map(|(a,b,c)|
    {
        let (z1,z2,z3,w1,w2,w3) = vec![puntos[a], puntos[b], puntos[c], 
                    nuevos_puntos[a], nuevos_puntos[b], nuevos_puntos[c]].into_iter()
            .map(|x| x.s_stereographic()).tuples().collect::<Vec<(Complex<f64>,Complex<f64>,Complex<f64>,Complex<f64>,Complex<f64>,Complex<f64>)>>()[0];
        MobiusTransformation::from_pairs(z1,z2,z3,w1,w2,w3).normalizar()
    }
    ).collect();

    // RANSAC
    //

    println!("{:?}",ts[0]);
    println!("{:?}",ts[1]);
    println!("{:?}",ts[2]);

    // Se calculan las transformaciones correspondientes a cada.
    //
    //
    let ts2: Vec<MobiusTransformation> = ternas.into_iter().map(|(a,b,c)|
    {
        let (z1,z2,z3,w1,w2,w3) = vec![puntos[a], puntos[b], puntos[c], 
                    nuevos_puntos_r[a], nuevos_puntos_r[b], nuevos_puntos_r[c]].into_iter()
            .map(|x| x.s_stereographic()).tuples().collect::<Vec<(Complex<f64>,Complex<f64>,Complex<f64>,Complex<f64>,Complex<f64>,Complex<f64>)>>()[0];
        MobiusTransformation::from_pairs(z1,z2,z3,w1,w2,w3).normalizar()
    }
    ).collect();

    // RANSAC
    //

    println!("{:?}",ts2[0]);
    println!("{:?}",ts2[1]);
    println!("{:?}",ts2[2]);

    // println!("{:?}",puntos);
    // println!("{:?}",nuevos_puntos);

    // println!("{}", geojson::ser::to_feature_collection_string(&puntos).unwrap());
    // println!("{}", geojson::ser::to_feature_collection_string(&nuevos_puntos).unwrap());

    assert!(true)
}
