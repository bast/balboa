pub fn densities_noddy(
    num_points: usize,
    aos: &[f64],
    density_matrix: &[f64],
    num_ao: usize,
) -> Vec<f64> {
    let mut densities = vec![0.0; num_points];

    for k in 0..num_ao {
        for l in 0..num_ao {
            let d = 2.0 * density_matrix[k * num_ao + l];
            for p in 0..num_points {
                densities[p] += aos[k * num_points + p] * d * aos[l * num_points + p];
            }
        }
    }

    densities
}
