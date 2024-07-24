use crate::geo::GeographicPoint;
// use crate::mobius::MobiusTransformation;
use num_complex::Complex;
use num_traits::Float;

// Definir el trait
// El trait de proceso tiene una función de aplicación que aplica la 
// transformación con la que el proceso se inicializó a una serie de 
// puntos geográficos.

trait Process {
    fn apply(&self, points: Vec<GeographicPoint>) -> Vec<Vec<GeographicPoint>>;
}

// Macro to generate the process structs with named inputs and a transformation function
macro_rules! define_process {
    // Nombre del struct, Transformación a aplicar, Nombres de las entradas y tipo.
    ($struct_name:ident,$transformation:ident, $($input:ident : $T:ty),+ ) => {
        // Crea un struct
        pub struct $struct_name {
            // con las inputs, que son los procesos de ruido 
            $( $input: Vec<$T>, )+
            // Y la transformación que aplicará
            $transformation: fn(&GeographicPoint, $($T),+) -> GeographicPoint,
        }

        // Implementa el trait process
        impl Process for $struct_name
        {
            // La aplicación de la transformaciòn
            fn apply(&self, points: Vec<GeographicPoint>) -> Vec<Vec<GeographicPoint>> {

                let mut instants = Vec::new();

                // toma la longitud de la entrada.
                let instant_count = {
                    let mut count = 0;
                    $(
                        count = self.$input.len();
                    )+
                    count
                };

                // para cada indice de la entrada
                for i in 0..instant_count {
                    
                    
                    let mut transformed_points = Vec::new();
                    
                    //toma las entradas.
                    $(
                        let $input = self.$input[i];
                    )+

                    // Para todos los puntos
                    for point in &points {
                        let transformed_point = (self.$transformation)(point, $($input),+);
                        // Agrega el punto transformado
                        transformed_points.push(transformed_point);
                    }

                    instants.push(transformed_points);
                }

                instants
            }
        }
    };
}

// Generate structs with 4 named inputs
define_process!(Process4, Transformation, Input1: f64, Input2: f64, Input3: f64, Input4: f64);

// Example transformation function that uses named inputs
fn example_transformation(
    point: &GeographicPoint,
    input1: f64,
    input2: f64,
    input3: f64,
    input4: f64,
) -> GeographicPoint
{
    // Your transformation logic here using named inputs
    // Example logic - summing up named inputs to the point coordinates
    let lat = point.latitude + input1;
    let lon = point.longitude + input2;

    GeographicPoint::new(lat,lon)
}

// Lo que sigue: Definir las transformaciones para los dos tipos de procesos
// que buscamos.


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn main() {
        // Example usage of Process4 struct
        let process = Process4 {
            Input1: vec![1.0, 2.0, 3.0],
            Input2: vec![4.0, 5.0, 6.0],
            Input3: vec![7.0, 8.0, 9.0],
            Input4: vec![10.0, 11.0, 12.0],
            Transformation: (example_transformation)
        };

        let points = vec![
            GeographicPoint::new(0.0, 0.0),
            GeographicPoint::new(1.0, 1.0),
            GeographicPoint::new(2.0, 2.0),
        ];

        let result = process.apply(points);
        println!("{:?}", result);
        assert!(true)
    }
}



