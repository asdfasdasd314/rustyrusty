use noise::{NoiseFn, Simplex};

pub fn create_simplex_noise() -> Vec<Vec<f64>> {
    let simplex = Simplex::new(0);
    let mut noise: Vec<Vec<f64>> = Vec::with_capacity(5);
    for y in 0..5 {
        let mut row: Vec<f64> = Vec::with_capacity(5);
        for x in 0..5 {
            row.push(simplex.get([x as f64 * 0.1, y as f64 * 0.1]));
        }
        noise.push(row);
    }

    noise
}
