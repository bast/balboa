use blas::*;
use std::collections::HashMap;

pub fn densities_noddy(
    num_points: usize,
    aos: &HashMap<(usize, usize, usize), Vec<f64>>,
    density_matrix: &[f64],
    num_ao: usize,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let ao_0 = aos.get(&(0, 0, 0)).unwrap();
    let ao_x = aos.get(&(1, 0, 0)).unwrap();
    let ao_y = aos.get(&(0, 1, 0)).unwrap();
    let ao_z = aos.get(&(0, 0, 1)).unwrap();

    let mut densities = vec![vec![0.0; num_points]; 4];
    for k in 0..num_ao {
        for l in 0..num_ao {
            let d = 2.0 * density_matrix[k * num_ao + l];
            for p in 0..num_points {
                densities[0][p] += ao_0[k * num_points + p] * d * ao_0[l * num_points + p];
                densities[1][p] += d
                    * (ao_x[k * num_points + p] * ao_0[l * num_points + p]
                        + ao_0[k * num_points + p] * ao_x[l * num_points + p]);
                densities[2][p] += d
                    * (ao_y[k * num_points + p] * ao_0[l * num_points + p]
                        + ao_0[k * num_points + p] * ao_y[l * num_points + p]);
                densities[3][p] += d
                    * (ao_z[k * num_points + p] * ao_0[l * num_points + p]
                        + ao_0[k * num_points + p] * ao_z[l * num_points + p]);
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
    aos: &HashMap<(usize, usize, usize), Vec<f64>>,
    density_matrix: &[f64],
    density_matrix_is_symmetric: bool,
    num_ao: usize,
    threshold: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let ao_0 = aos.get(&(0, 0, 0)).unwrap();
    let ao_x = aos.get(&(1, 0, 0)).unwrap();
    let ao_y = aos.get(&(0, 1, 0)).unwrap();
    let ao_z = aos.get(&(0, 0, 1)).unwrap();

    // compress aos and keep track of mapping original -> compressed
    let mut aos_compressed: Vec<Vec<f64>> = Vec::new();
    for _ in 0..4 {
        aos_compressed.push(Vec::new());
    }
    let mut compression_mapping: Vec<(usize, usize)> = Vec::new();
    let mut k_compressed = 0;
    for k in 0..num_ao {
        let mut keep_this_batch = false;
        for p in 0..num_points {
            if ao_0[k * num_points + p].abs() > threshold {
                keep_this_batch = true;
                break;
            }
        }
        if keep_this_batch {
            for p in 0..num_points {
                aos_compressed[0].push(ao_0[k * num_points + p]);
                aos_compressed[1].push(ao_x[k * num_points + p]);
                aos_compressed[2].push(ao_y[k * num_points + p]);
                aos_compressed[3].push(ao_z[k * num_points + p]);
            }
            compression_mapping.push((k_compressed, k));
            k_compressed += 1;
        }
    }
    let num_ao_compressed = aos_compressed[0].len() / num_points;

    // compress density matrix
    let mut density_matrix_compressed = vec![0.0; num_ao_compressed * num_ao_compressed];
    for (kc, k) in &compression_mapping {
        let kc_offset = kc * num_ao_compressed;
        let k_offset = k * num_ao;
        for (lc, l) in &compression_mapping {
            density_matrix_compressed[kc_offset + lc] = density_matrix[k_offset + l];
        }
    }

    let mut x_matrix = vec![vec![0.0; num_ao_compressed * num_points]; 4];
    let block_length = num_points as i32;
    let k_aoc_num = num_ao_compressed as i32;
    let l_aoc_num = num_ao_compressed as i32;

    let mut densities = vec![vec![0.0; num_points]; 4];

    if density_matrix_is_symmetric {
        let si = b'R';
        let up = b'U';
        let m = block_length;
        let n = k_aoc_num;
        let lda = n;
        let ldb = m;
        let ldc = m;
        let alpha = 2.0;
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
                &aos_compressed[0],
                ldb,
                beta,
                &mut x_matrix[0],
                ldc,
            );
        }

        let prefactors = vec![1.0, 2.0, 2.0, 2.0];
        for (ixyz, prefactor) in prefactors.iter().enumerate() {
            for k in 0..num_ao_compressed {
                let k_offset = k * num_points;
                for p in 0..num_points {
                    densities[ixyz][p] +=
                        prefactor * aos_compressed[ixyz][k_offset + p] * x_matrix[0][k_offset + p];
                }
            }
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
        let alpha = 2.0;
        let beta = 0.0;
        for ixyz in 0..4 {
            unsafe {
                dgemm(
                    ta,
                    tb,
                    m,
                    n,
                    k,
                    alpha,
                    &aos_compressed[ixyz],
                    lda,
                    &density_matrix_compressed,
                    ldb,
                    beta,
                    &mut x_matrix[ixyz],
                    ldc,
                );
            }
        }

        for k in 0..num_ao_compressed {
            let k_offset = k * num_points;
            for p in 0..num_points {
                densities[0][p] += aos_compressed[0][k_offset + p] * x_matrix[0][k_offset + p]
            }
        }
        for ixyz in 1..=3 {
            for k in 0..num_ao_compressed {
                let k_offset = k * num_points;
                for p in 0..num_points {
                    densities[ixyz][p] += aos_compressed[ixyz][k_offset + p]
                        * x_matrix[0][k_offset + p]
                        + aos_compressed[0][k_offset + p] * x_matrix[ixyz][k_offset + p];
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
