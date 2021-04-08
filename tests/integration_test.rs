use balboa;
use balboa::Point;

use rand::{rngs::StdRng, Rng, SeedableRng};
use std::fmt::Debug;
use std::fs;
use std::str::FromStr;
use std::time::Instant;

pub fn get_rng(seed: [u8; 32]) -> StdRng {
    StdRng::from_seed(seed)
}

pub fn random_points<R: Rng>(rng: &mut R, num_points: usize, side_length: f64) -> Vec<Point> {
    let xmax = 0.5 * side_length;
    let xmin = -xmax;
    let ymax = 0.5 * side_length;
    let ymin = -ymax;
    let zmax = 0.5 * side_length;
    let zmin = -zmax;

    (0..num_points)
        .map(|_| Point {
            x: rng.gen_range(xmin..xmax),
            y: rng.gen_range(ymin..ymax),
            z: rng.gen_range(zmin..zmax),
        })
        .collect()
}

fn read_vector<T: FromStr>(file_name: &str) -> Vec<T>
where
    <T as FromStr>::Err: Debug,
{
    let error_message = format!("something went wrong reading file {}", file_name);
    let contents = fs::read_to_string(file_name).expect(&error_message);
    contents
        .lines()
        .map(|s| s.trim().parse().unwrap())
        .collect()
}

fn read_square_matrix(file_name: &str, n: usize) -> Vec<f64> {
    let error_message = format!("something went wrong reading file {}", file_name);
    let contents = fs::read_to_string(file_name).expect(&error_message);

    let mut matrix = vec![0.0; n * n];

    for line in contents.lines() {
        let words: Vec<&str> = line.split_whitespace().collect();
        let i: usize = words[0].parse().unwrap();
        let f = words[1].parse().unwrap();
        matrix[i] = f;
    }

    matrix
}

fn floats_are_same(value: f64, reference: f64, threshold: f64) -> bool {
    let absolute_error = (value - reference).abs();
    if reference.abs() > threshold {
        let relative_error = (absolute_error / reference).abs();
        relative_error < threshold
    } else {
        absolute_error < threshold
    }
}

#[test]
fn density_noddy() {
    let basis = balboa::example_basis(false);

    let c_to_s_matrices = balboa::cartesian_to_spherical_matrices();

    let points_bohr = vec![
        Point {
            x: 1.7,
            y: 0.0,
            z: 0.0,
        },
        Point {
            x: 1.847135,
            y: 1.505297e-02,
            z: 1.505297e-02,
        },
        Point {
            x: 1.608568,
            y: -9.143180e-02,
            z: 1.437394e-01,
        },
        Point {
            x: 1.806738,
            y: -2.317019e-01,
            z: 3.402263e-01,
        },
    ];

    let aos = balboa::aos_noddy(0, &points_bohr, &basis, &c_to_s_matrices);

    let density_matrix = read_square_matrix("tests/dmat.txt", basis.num_ao_spherical);

    let densities_reference = vec![
        427.74880135855784,
        32.362597237925904,
        15.928938748468175,
        2.0615806308226245,
    ];

    let densities = balboa::densities_noddy(
        points_bohr.len(),
        &aos,
        &density_matrix,
        basis.num_ao_spherical,
    );
    for (&x, &x_reference) in densities.iter().zip(densities_reference.iter()) {
        assert!(floats_are_same(x, x_reference, 1.0e-12));
    }

    let densities = balboa::densities(
        points_bohr.len(),
        &aos,
        &density_matrix,
        basis.num_ao_spherical,
    );
    for (&x, &x_reference) in densities.iter().zip(densities_reference.iter()) {
        assert!(floats_are_same(x, x_reference, 1.0e-12));
    }
}

#[test]
fn ao_derivatives_noddy() {
    let basis = balboa::example_basis(true);

    let c_to_s_matrices = balboa::cartesian_to_spherical_matrices();

    let points_bohr = vec![
        Point {
            x: -1.46254302355,
            y: 1.38973494775,
            z: 1.05509847591,
        },
        Point {
            x: -0.979723897042,
            y: -0.0182596516322,
            z: -0.202035740845,
        },
        Point {
            x: 0.606371890891,
            y: 1.15489340454,
            z: -1.6245616529,
        },
        Point {
            x: -1.88661009391,
            y: 1.34306041568,
            z: -0.26893172838,
        },
    ];

    let aos = balboa::aos_noddy(2, &points_bohr, &basis, &c_to_s_matrices);

    // ao 1, point 1
    // ao 1, point 2
    // ao 1, point 3
    // ao 1, point 4
    // ao 2, point 1
    // ...
    let aos_reference: Vec<f64> = read_vector("tests/reference.txt");

    for (&ao, &ao_reference) in aos.iter().zip(aos_reference.iter()) {
        assert!(floats_are_same(ao, ao_reference, 1.0e-12));
    }
}

#[ignore]
#[test]
fn benchmark() {
    let mut rng = get_rng([0; 32]);
    let num_points = 50_000;
    let side_length = 2.0;
    let points_bohr = random_points(&mut rng, num_points, side_length);

    let basis = balboa::example_basis(true);

    let c_to_s_matrices = balboa::cartesian_to_spherical_matrices();

    let start = Instant::now();
    let _aos = balboa::aos_noddy(2, &points_bohr, &basis, &c_to_s_matrices);
    println!("time elapsed in aos_noddy: {:?}", start.elapsed());
}
