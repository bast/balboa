use std::time::Instant;

fn main() {
    let mut rng = balboa::get_rng([0; 32]);
    let num_points = 50_000;
    let side_length = 2.0;
    let points_bohr = balboa::random_points(&mut rng, num_points, side_length);

    let basis = balboa::example_basis();

    let c_to_s_matrices = balboa::cartesian_to_spherical_matrices();

    let start = Instant::now();
    let _aos = balboa::aos_noddy(2, &points_bohr, &basis, &c_to_s_matrices);
    println!("time elapsed in aos_noddy: {:?}", start.elapsed());
}
