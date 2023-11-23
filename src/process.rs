// Function to create a Möbius transformation from a geographic point and perturbations
fn create_mobius_transformation(
    base_point: GeographicPoint,
    perturbations: Vec<(f64, f64, f64, f64)>,
) -> Vec<MobiusTransformation> {
    perturbations
        .into_iter()
        .map(|perturbation| {
            let base_transform = MobiusTransformation::from_pole(
                base_point.latitude.into(),
                base_point.longitude.into(),
                1.0,
            );

            let a = base_transform.a + perturbation.0;
            let b = base_transform.b + perturbation.1;
            let c = base_transform.c + perturbation.2;
            let d = base_transform.d + perturbation.3;

            MobiusTransformation { a, b, c, d }
        })
        .collect()
}

// Function to create a Möbius aligned transformation from a geographic point and noise scale
fn create_mobius_aligned_transformation(
    base_point: GeographicPoint,
    noise_scale: f64,
    num_points: usize,
) -> Vec<MobiusAlignedTransformation> {
    let base_transform = MobiusTransformation::from_pole(
        base_point.latitude.into(),
        base_point.longitude.into(),
        1.0,
    );

    let normal = Normal::new(0.0, noise_scale).unwrap();
    let mut rng = ChaCha8Rng::seed_from_u64(42);

    let perturbations: Vec<MobiusTransformation> = (0..num_points)
        .map(|_| {
            let perturbation = (
                Complex::new(rng.sample(&normal), 0.0),
                Complex::new(rng.sample(&normal), 0.0),
                Complex::new(rng.sample(&normal), 0.0),
                Complex::new(rng.sample(&normal), 0.0),
            );

            MobiusTransformation {
                a: perturbation.0,
                b: perturbation.1,
                c: perturbation.2,
                d: perturbation.3,
            }
        })
        .collect();

    perturbations
        .into_iter()
        .map(|perturbation| MobiusAlignedTransformation {
            base_transform: base_transform.clone(),
            noise_scale: perturbation,
        })
        .collect()
}

// Function to generate a temporal series of points
fn generate_temporal_series(base_point: GeographicPoint, num_points: usize) -> Vec<GeographicPoint> {
    let mut rng = rand::thread_rng();
    (0..num_points)
        .map(|_| GeographicPoint {
            latitude: base_point.latitude + rng.gen_range(-5.0..5.0),
            longitude: base_point.longitude + rng.gen_range(-5.0..5.0),
        })
        .collect()
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
