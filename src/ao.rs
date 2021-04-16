#![allow(clippy::too_many_arguments)]

use std::collections::HashMap;
use std::time::Instant;

use crate::basis::Basis;
use crate::diff;
use crate::generate;
use crate::transform;

fn coordinates(
    shell_center_coordinates: (f64, f64, f64),
    points_bohr: &[(f64, f64, f64)],
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let (cx, cy, cz) = shell_center_coordinates;

    let mut pxs = Vec::new();
    let mut pys = Vec::new();
    let mut pzs = Vec::new();
    let mut p2s = Vec::new();

    for &(x, y, z) in points_bohr {
        let px = x - cx;
        let py = y - cy;
        let pz = z - cz;
        let p2 = px * px + py * py + pz * pz;

        pxs.push(px);
        pys.push(py);
        pzs.push(pz);
        p2s.push(p2);
    }

    (pxs, pys, pzs, p2s)
}

pub fn aos_noddy(
    max_geo_derv_order: usize,
    points_bohr: &[(f64, f64, f64)],
    basis: &Basis,
    c_to_s_matrices: &HashMap<usize, Vec<(usize, usize, f64)>>,
) -> Vec<f64> {
    let max_l_value = basis.shell_l_quantum_numbers.iter().max().unwrap();
    assert!(
        c_to_s_matrices.contains_key(&max_l_value),
        "increase max l value in cartesian_to_spherical_matrices"
    );

    let num_points = points_bohr.len();
    let mut aos = Vec::new();

    let mut time_ms_prefactors: u128 = 0;
    let mut time_ms_multiply: u128 = 0;
    let mut time_ms_transform: u128 = 0;

    // FIXME loop over primitives should be outside the loop over geo derv orders
    // otherwise we recompute the s functions over and over
    for geo_derv_order in 0..=max_geo_derv_order {
        for &geo_derv_orders in generate::get_ijk_list(geo_derv_order).iter() {
            let mut offset = 0;
            for ishell in 0..basis.num_shells {
                let (pxs, pys, pzs, p2s) =
                    coordinates(basis.shell_centers_coordinates[ishell], &points_bohr);

                let num_primitives = basis.shell_num_primitives[ishell];
                let l = basis.shell_l_quantum_numbers[ishell];

                let timer = Instant::now();
                let mut gaussians = vec![vec![0.0; num_points]; max_geo_derv_order + 1];

                for ip in offset..(offset + num_primitives) {
                    let e = basis.primitive_exponents[ip];
                    let c = basis.contraction_coefficients[ip];

                    let ts: Vec<_> = p2s.iter().map(|p2| c * (-e * p2).exp()).collect();

                    for geo_derv_order in 0..=max_geo_derv_order {
                        let x = (-2.0 * e).powi(geo_derv_order as i32);
                        for ipoint in 0..num_points {
                            gaussians[geo_derv_order][ipoint] += x * ts[ipoint];
                        }
                    }
                }
                time_ms_prefactors += timer.elapsed().as_millis();

                let timer = Instant::now();
                let cartesian_deg = (l + 1) * (l + 2) / 2;
                let mut aos_c = vec![0.0; num_points * cartesian_deg];
                let mut n = 0;
                for &cartesian_orders in generate::get_ijk_list(l).iter() {
                    let map = diff::differentiate(cartesian_orders, geo_derv_orders);

                    for (key, value) in map {
                        let i = key[0] as i32;
                        let j = key[1] as i32;
                        let k = key[2] as i32;
                        let p = key[3];

                        let v = value as f64;

                        for ipoint in 0..num_points {
                            aos_c[n + ipoint] += v
                                * gaussians[p][ipoint]
                                * pxs[ipoint].powi(i)
                                * pys[ipoint].powi(j)
                                * pzs[ipoint].powi(k);
                        }
                    }

                    n += num_points;
                }
                time_ms_multiply += timer.elapsed().as_millis();

                let timer = Instant::now();
                if l < 2 {
                    aos.append(&mut aos_c);
                } else {
                    aos.append(&mut transform::transform_to_spherical(
                        num_points,
                        &aos_c,
                        l,
                        &c_to_s_matrices.get(&l).unwrap(),
                    ));
                }
                time_ms_transform += timer.elapsed().as_millis();

                offset += num_primitives;
            }
        }
    }

    println!("time spent in prefactors: {} ms", time_ms_prefactors);
    println!("time spent in multiply: {} ms", time_ms_multiply);
    println!("time spent in transform: {} ms", time_ms_transform);

    aos
}
