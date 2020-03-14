use balboa;
use balboa::Point;

use std::fmt::Debug;
use std::fs;
use std::str::FromStr;

fn read_vector<T: FromStr>(file_name: &str) -> Vec<T>
where
    <T as FromStr>::Err: Debug,
{
    let error_message = format!("something went wrong reading file {}", file_name);
    let contents = fs::read_to_string(file_name).expect(&error_message);
    let v = contents.lines().map(|s| s.parse().unwrap()).collect();

    return v;
}

fn floats_are_same(f1: f64, f2: f64) -> bool {
    let d = f1 - f2;
    return d.abs() < 1.0e-7;
}

#[test]
fn undifferentiated() {
    let basis = balboa::example_basis();

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
        assert!(floats_are_same(ao, ao_reference));
    }
}
