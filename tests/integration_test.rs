use balboa;

use rand::Rng;
use std::fmt::Debug;
use std::fs;
use std::str::FromStr;
use std::time::Instant;

fn random_points(num_points: usize, side_length: f64) -> Vec<(f64, f64, f64)> {
    let xmax = 0.5 * side_length;
    let xmin = -xmax;
    let ymax = 0.5 * side_length;
    let ymin = -ymax;
    let zmax = 0.5 * side_length;
    let zmin = -zmax;

    let mut rng = rand::thread_rng();

    (0..num_points)
        .map(|_| {
            (
                rng.gen_range(xmin..xmax),
                rng.gen_range(ymin..ymax),
                rng.gen_range(zmin..zmax),
            )
        })
        .collect()
}

fn random_matrix(n: usize, symmetric: bool) -> Vec<f64> {
    let mut rng = rand::thread_rng();

    let mut matrix = vec![0.0; n * n];

    if symmetric {
        for i in 0..n {
            for j in 0..=i {
                let r = rng.gen_range(-1.0..1.0);
                matrix[i * n + j] = r;
                matrix[j * n + i] = r;
            }
        }
    } else {
        for i in 0..n {
            for j in 0..n {
                matrix[i * n + j] = rng.gen_range(-1.0..1.0);
            }
        }
    }

    matrix
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

fn compare_vectors(v1: &[f64], v2: &[f64]) {
    for (&x1, &x2) in v1.iter().zip(v2.iter()) {
        assert!(floats_are_same(x1, x2, 1.0e-10));
    }
}

#[test]
fn density_noddy() {
    let basis = balboa::example_basis(false);
    let max_l_value = basis.shells.iter().fold(0, |m, s| m.max(s.l));
    let c_to_s_matrices = balboa::cartesian_to_spherical_matrices(max_l_value);

    let points_bohr = vec![
        (1.7, 0.0, 0.0),
        (1.847135, 1.505297e-02, 1.505297e-02),
        (1.608568, -9.143180e-02, 1.437394e-01),
        (1.806738, -2.317019e-01, 3.402263e-01),
    ];

    let aos = balboa::aos_noddy(1, &points_bohr, &basis, &c_to_s_matrices);

    let density_matrix = read_square_matrix("tests/dmat.txt", basis.num_ao_spherical);

    let densities_ref = vec![
        427.74880135855784,
        32.362597237925904,
        15.928938748468175,
        2.0615806308226245,
    ];

    let densities_x_ref = vec![
        -9.235703423787248,
        -543.7711358713716,
        118.94476088385002,
        -1.541565272853058,
    ];

    let densities_y_ref = vec![
        0.0,
        -55.08159215958256,
        117.73175332829445,
        3.012583667619674,
    ];

    let densities_z_ref = vec![
        0.0,
        -55.0815921595826,
        -185.08540337559853,
        -4.4236158386041105,
    ];

    let (densities, densities_x, densities_y, densities_z) = balboa::densities_noddy(
        points_bohr.len(),
        &aos,
        &density_matrix,
        basis.num_ao_spherical,
    );
    compare_vectors(&densities, &densities_ref);
    compare_vectors(&densities_x, &densities_x_ref);
    compare_vectors(&densities_y, &densities_y_ref);
    compare_vectors(&densities_z, &densities_z_ref);
}

#[test]
fn density() {
    let basis = balboa::example_basis_benchmark(10, 4);
    let max_l_value = basis.shells.iter().fold(0, |m, s| m.max(s.l));
    let c_to_s_matrices = balboa::cartesian_to_spherical_matrices(max_l_value);

    let num_points = 100;
    let side_length = 3.0;
    let points_bohr = random_points(num_points, side_length);

    let aos = balboa::aos_noddy(1, &points_bohr, &basis, &c_to_s_matrices);

    for symmetric in &[true, false] {
        let density_matrix = random_matrix(basis.num_ao_spherical, *symmetric);

        let start = Instant::now();
        let (densities_ref, densities_x_ref, densities_y_ref, densities_z_ref) =
            balboa::densities_noddy(
                points_bohr.len(),
                &aos,
                &density_matrix,
                basis.num_ao_spherical,
            );
        println!("time spent in densities_noddy: {:?}", start.elapsed());

        let start = Instant::now();
        let (densities, densities_x, densities_y, densities_z) = balboa::densities(
            points_bohr.len(),
            &aos,
            &density_matrix,
            *symmetric,
            basis.num_ao_spherical,
            1.0e-8,
        );
        compare_vectors(&densities, &densities_ref);
        compare_vectors(&densities_x, &densities_x_ref);
        compare_vectors(&densities_y, &densities_y_ref);
        compare_vectors(&densities_z, &densities_z_ref);
        println!(
            "time spent in blas version (symmetric={}): {:?}",
            symmetric,
            start.elapsed()
        );
    }
}

fn get_ijk_list(m: usize) -> Vec<(usize, usize, usize)> {
    let mut l = Vec::new();

    for a in 1..(m + 2) {
        for b in 1..=a {
            l.push((m + 1 - a, a - b, b - 1));
        }
    }

    l
}

#[test]
fn ao_derivatives_noddy() {
    let basis = balboa::example_basis(true);
    let max_l_value = basis.shells.iter().fold(0, |m, s| m.max(s.l));
    let c_to_s_matrices = balboa::cartesian_to_spherical_matrices(max_l_value);

    let points_bohr = vec![
        (-1.46254302355, 1.38973494775, 1.05509847591),
        (-0.979723897042, -0.0182596516322, -0.202035740845),
        (0.606371890891, 1.15489340454, -1.6245616529),
        (-1.88661009391, 1.34306041568, -0.26893172838),
    ];

    let aos = balboa::aos_noddy(2, &points_bohr, &basis, &c_to_s_matrices);

    let aos_ref: Vec<f64> = read_vector("tests/reference.txt");

    let mut n = 0;
    for geo_derv_order in 0..=2 {
        for (i, j, k) in get_ijk_list(geo_derv_order) {
            let batch = aos.get(&(i, j, k)).unwrap();
            let mut m = 0;
            for _ in 0..(basis.num_ao_spherical * points_bohr.len()) {
                assert!(floats_are_same(aos_ref[n], batch[m], 1.0e-10));
                n += 1;
                m += 1;
            }
        }
    }
}

#[ignore]
#[test]
fn ao_benchmark() {
    let basis = balboa::example_basis_benchmark(15, 4);
    let max_l_value = basis.shells.iter().fold(0, |m, s| m.max(s.l));
    let c_to_s_matrices = balboa::cartesian_to_spherical_matrices(max_l_value);

    let num_points = 50_000;
    let side_length = 2.0;
    let points_bohr = random_points(num_points, side_length);

    let start = Instant::now();
    let _aos = balboa::aos_noddy(2, &points_bohr, &basis, &c_to_s_matrices);
    println!("time spent in aos_noddy: {:?}", start.elapsed());
}
