#![allow(clippy::needless_return)]

use crate::basis::Basis;
use crate::generate;
use crate::limits;
use crate::multiply;
use crate::point::Point;
use std::time::Instant;

fn g_batch(c: f64, e: f64, gaussians: &mut [Vec<f64>], p2s: &[f64], max_geo_derv_order: usize) {
    for ipoint in 0..limits::BATCH_LENGTH {
        let mut t = c * (-e * p2s[ipoint]).exp();
        for geo_derv_order in 0..=max_geo_derv_order {
            gaussians[ipoint][geo_derv_order] += t;
            t *= -2.0 * e;
        }
    }
}

fn compute_prefactors(
    max_geo_derv_order: usize,
    num_batches: usize,
    p2s: Vec<f64>,
    basis: &Basis,
    offset: usize,
    num_primitives: usize,
) -> Vec<Vec<f64>> {
    let num_points = p2s.len();
    let mut gaussians = vec![vec![0.0; max_geo_derv_order + 1]; num_points];

    for ip in offset..(offset + num_primitives) {
        let e = basis.primitive_exponents[ip];
        let c = basis.contraction_coefficients[ip];
        for ibatch in 0..num_batches {
            let offset = ibatch * limits::BATCH_LENGTH;
            let end = offset + limits::BATCH_LENGTH;
            g_batch(
                c,
                e,
                &mut gaussians[offset..end],
                &p2s[offset..end],
                max_geo_derv_order,
            );
        }
    }

    return gaussians;
}

fn compute_gaussians(
    geo_derv_orders: (usize, usize, usize),
    c_to_s_matrices: &[Vec<Vec<f64>>],
    num_batches: usize,
    l: usize,
    pxs: &[f64],
    pys: &[f64],
    pzs: &[f64],
    gaussians: &[Vec<f64>],
) -> Vec<f64> {
    let num_points = pxs.len();

    let cartesian_deg = (l + 1) * (l + 2) / 2;
    let mut aos_c = vec![0.0; num_points * cartesian_deg];

    // TODO create shortcut for s functions when multiplying
    let mut n = 0;
    for &cartesian_orders in generate::get_ijk_list(l).iter() {
        for ibatch in 0..num_batches {
            let ao_offset = n;
            let ao_end = n + limits::BATCH_LENGTH;
            let point_offset = ibatch * limits::BATCH_LENGTH;
            let point_end = point_offset + limits::BATCH_LENGTH;
            multiply::multiply_batch(
                cartesian_orders,
                geo_derv_orders,
                &gaussians[point_offset..point_end],
                &pxs[point_offset..point_end],
                &pys[point_offset..point_end],
                &pzs[point_offset..point_end],
                &mut aos_c[ao_offset..ao_end],
            );
            n += limits::BATCH_LENGTH;
        }
    }

    if l < 2 {
        return aos_c;
    } else {
        let aos_s = transform_to_spherical(num_points, &aos_c, l, &c_to_s_matrices[l]);
        return aos_s;
    }
}

fn transform_to_spherical(
    num_points: usize,
    aos_c: &[f64],
    l: usize,
    c_to_s_matrix: &[Vec<f64>],
) -> Vec<f64> {
    let spherical_deg = 2 * l + 1;
    let mut aos_s = vec![0.0; spherical_deg * num_points];

    for (i, row) in c_to_s_matrix.iter().enumerate() {
        let ioff = i * num_points;
        for (j, element) in row.iter().enumerate() {
            if element.abs() > std::f64::EPSILON {
                let joff = j * num_points;
                for ipoint in 0..num_points {
                    aos_s[joff + ipoint] += element * aos_c[ioff + ipoint];
                }
            }
        }
    }

    return aos_s;
}

fn coordinates(
    shell_centers_coordinates: (f64, f64, f64),
    points_bohr: &[Point],
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let (x, y, z) = shell_centers_coordinates;

    let mut pxs = Vec::new();
    let mut pys = Vec::new();
    let mut pzs = Vec::new();
    let mut p2s = Vec::new();

    for point in points_bohr.iter() {
        let px = point.x - x;
        let py = point.y - y;
        let pz = point.z - z;
        let p2 = px * px + py * py + pz * pz;

        pxs.push(px);
        pys.push(py);
        pzs.push(pz);
        p2s.push(p2);
    }

    return (pxs, pys, pzs, p2s);
}

pub fn aos_noddy(
    max_geo_derv_order: usize,
    points_bohr: &[Point],
    basis: &Basis,
    c_to_s_matrices: &[Vec<Vec<f64>>],
) -> Vec<f64> {
    let num_points = points_bohr.len();

    assert!(
        num_points % limits::BATCH_LENGTH == 0,
        "num_points must be multiple of BATCH_LENGTH"
    );
    let num_batches = num_points / limits::BATCH_LENGTH;

    let mut aos = Vec::new();

    let mut time_ms_gaussian: u128 = 0;
    let mut time_ms_multiply: u128 = 0;

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
                let gaussians = compute_prefactors(
                    max_geo_derv_order,
                    num_batches,
                    p2s,
                    &basis,
                    offset,
                    num_primitives,
                );
                time_ms_gaussian += timer.elapsed().as_millis();

                let timer = Instant::now();
                let mut aos_s = compute_gaussians(
                    geo_derv_orders,
                    c_to_s_matrices,
                    num_batches,
                    l,
                    &pxs,
                    &pys,
                    &pzs,
                    &gaussians,
                );
                time_ms_multiply += timer.elapsed().as_millis();

                aos.append(&mut aos_s);

                offset += num_primitives;
            }
        }
    }

    println!("time spent in exp: {} ms", time_ms_gaussian);
    println!("time spent in multiply: {} ms", time_ms_multiply);

    return aos;
}
