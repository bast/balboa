#![allow(clippy::needless_return)]

use rand::{rngs::StdRng, Rng, SeedableRng};

use crate::point::Point;

pub fn get_rng(seed: [u8; 32]) -> StdRng {
    let rng = StdRng::from_seed(seed);
    return rng;
}

pub fn random_points<R: Rng>(rng: &mut R, num_points: usize, side_length: f64) -> Vec<Point> {
    let mut points = Vec::new();

    let xmax = 0.5 * side_length;
    let xmin = -xmax;
    let ymax = 0.5 * side_length;
    let ymin = -ymax;
    let zmax = 0.5 * side_length;
    let zmin = -zmax;

    for _ in 0..num_points {
        points.push(Point {
            x: rng.gen_range(xmin, xmax),
            y: rng.gen_range(ymin, ymax),
            z: rng.gen_range(zmin, zmax),
        });
    }

    return points;
}
