use balboa;

fn main() {
    let mut rng = balboa::get_rng([0; 32]);
    let num_points = 10;
    let side_length = 2.0;
    let _points = balboa::random_points(&mut rng, num_points, side_length);
}
