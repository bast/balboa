use blas::*;

pub fn densities_noddy(
    num_points: usize,
    aos: &[f64],
    density_matrix: &[f64],
    num_ao: usize,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut densities = vec![vec![0.0; num_points]; 4];

    for k in 0..num_ao {
        for l in 0..num_ao {
            let d = 2.0 * density_matrix[k * num_ao + l];
            for p in 0..num_points {
                densities[0][p] += aos[k * num_points + p] * d * aos[l * num_points + p];
            }
            for ixyz in 1..=3 {
                let slice_offset = (num_ao * num_points) * ixyz;
                for p in 0..num_points {
                    densities[ixyz][p] += d
                        * (aos[slice_offset + k * num_points + p] * aos[l * num_points + p]
                            + aos[k * num_points + p] * aos[slice_offset + l * num_points + p]);
                }
            }
        }
    }

    (
        densities[0].to_vec(), // density
        densities[1].to_vec(), // gradient_x
        densities[2].to_vec(), // gradient_y
        densities[3].to_vec(), // gradient_z
    )
}

pub fn densities(
    num_points: usize,
    aos: &[f64],
    density_matrix: &[f64],
    density_matrix_is_symmetric: bool,
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

    let mut x_matrix = vec![0.0; num_ao_compressed * num_points];
    let block_length = num_points as i32;
    let k_aoc_num = num_ao_compressed as i32;
    let l_aoc_num = num_ao_compressed as i32;

    if density_matrix_is_symmetric {
        let si = b'R';
        let up = b'U';
        let m = block_length;
        let n = k_aoc_num;
        let lda = n;
        let ldb = m;
        let ldc = m;
        let alpha = 1.0;
        let beta = 0.0;
        unsafe {
            dsymm(
                si,
                up,
                m,
                n,
                alpha,
                &density_matrix_compressed,
                lda,
                &aos_compressed,
                ldb,
                beta,
                &mut x_matrix,
                ldc,
            );
        }
    } else {
        let ta = b'N';
        let tb = b'N';
        let m = block_length;
        let n = k_aoc_num;
        let k = l_aoc_num;
        let lda = m;
        let ldb = k;
        let ldc = m;
        let alpha = 1.0;
        let beta = 0.0;
        unsafe {
            dgemm(
                ta,
                tb,
                m,
                n,
                k,
                alpha,
                &aos_compressed,
                lda,
                &density_matrix_compressed,
                ldb,
                beta,
                &mut x_matrix,
                ldc,
            );
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
