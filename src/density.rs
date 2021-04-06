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

pub fn densities_noddy2(
    num_points: usize,
    aos: &[f64],
    density_matrix: &[f64],
    num_ao: usize,
) -> Vec<f64> {
    let threshold = 0.003;

    // compress aos and keep track of mapping original -> compressed
    let mut aos_compressed: Vec<f64> = Vec::new();
    let mut compression_mapping: Vec<(usize, usize)> = Vec::new();
    let mut k_compressed = 0;
    for k in 0..num_ao {
        let mut keep_this_batch = false;
        for p in 0..num_points {
            if aos[k * num_points + p].abs() > threshold {
                keep_this_batch = true;
                break;
            }
        }
        if keep_this_batch {
            for p in 0..num_points {
                aos_compressed.push(aos[k * num_points + p]);
            }
            compression_mapping.push((k_compressed, k));
            k_compressed += 1;
        }
    }
    let num_ao_compressed = aos_compressed.len() / num_points;

    // compress density matrix
    let mut density_matrix_compressed = vec![0.0; num_ao_compressed * num_ao_compressed];
    for (kc, k) in &compression_mapping {
        let kc_offset = kc * num_ao_compressed;
        let k_offset = k * num_ao;
        for (lc, l) in &compression_mapping {
            density_matrix_compressed[kc_offset + lc] = density_matrix[k_offset + l];
        }
    }

    // this will be a matrix multiplication later
    let mut x_matrix = vec![0.0; num_ao_compressed * num_points];
    for k in 0..num_ao_compressed {
        for p in 0..num_points {
            for l in 0..num_ao_compressed {
                x_matrix[k * num_points + p] += density_matrix_compressed
                    [k * num_ao_compressed + l]
                    * aos_compressed[l * num_points + p];
            }
        }
    }

    // finally assemble densities in the second step
    let mut densities = vec![0.0; num_points];
    for k in 0..num_ao_compressed {
        let k_offset = k * num_points;
        for p in 0..num_points {
            densities[p] += 2.0 * aos[k_offset + p] * x_matrix[k_offset + p];
        }
    }

    densities
}
