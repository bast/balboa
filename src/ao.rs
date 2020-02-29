use crate::basis::Basis;
use crate::generate;
use crate::point::Point;
use std::time::Instant;

fn compute_gaussian(
    p2s: Vec<f64>,
    basis: &Basis,
    offset: usize,
    num_primitives: usize,
) -> Vec<f64> {
    let num_points = p2s.len();
    let mut gaussians = vec![0.0; num_points];

    for ip in offset..(offset + num_primitives) {
        let e = basis.primitive_exponents[ip];
        let c = basis.contraction_coefficients[ip];
        for ipoint in 0..num_points {
            gaussians[ipoint] += c * (-e * p2s[ipoint]).exp();
        }
    }

    return gaussians;
}

fn transform_to_spherical(
    num_points: usize,
    aos_c: &Vec<f64>,
    spherical_deg: usize,
    c_to_s_matrix: &Vec<Vec<f64>>,
) -> Vec<f64> {
    let mut aos_s = vec![0.0; spherical_deg * num_points];

    for (i, row) in c_to_s_matrix.iter().enumerate() {
        let ioff = i * num_points;
        for (j, element) in row.iter().enumerate() {
            let joff = j * num_points;
            for ipoint in 0..num_points {
                aos_s[joff + ipoint] += element * aos_c[ioff + ipoint];
            }
        }
    }

    return aos_s;
}

fn coordinates(
    shell_centers_coordinates: (f64, f64, f64),
    points_bohr: &Vec<Point>,
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

fn aos_c_batch(
    l: usize,
    gaussians: Vec<f64>,
    pxs: Vec<f64>,
    pys: Vec<f64>,
    pzs: Vec<f64>,
) -> Vec<f64> {
    let mut aos_c = Vec::new();

    for (i, j, k) in generate::get_ijk_list(l).iter() {
        for ipoint in 0..pxs.len() {
            aos_c.push(
                gaussians[ipoint]
                    * pxs[ipoint].powi(*i as i32)
                    * pys[ipoint].powi(*j as i32)
                    * pzs[ipoint].powi(*k as i32),
            );
        }
    }

    return aos_c;
}

pub fn aos_noddy(
    points_bohr: Vec<Point>,
    basis: &Basis,
    c_to_s_matrices: &Vec<Vec<Vec<f64>>>,
) -> Vec<f64> {
    let num_points = points_bohr.len();

    let mut offset = 0;
    let mut aos = Vec::new();

    let mut time_ms_gaussian: u128 = 0;
    let mut time_ms_multiply: u128 = 0;
    let mut time_ms_transform: u128 = 0;

    for ishell in 0..basis.num_shells {
        let (pxs, pys, pzs, p2s) =
            coordinates(basis.shell_centers_coordinates[ishell], &points_bohr);

        let num_primitives = basis.shell_num_primitives[ishell];

        let timer = Instant::now();
        let gaussians = compute_gaussian(p2s, &basis, offset, num_primitives);
        time_ms_gaussian += timer.elapsed().as_millis();

        let l = basis.shell_l_quantum_numbers[ishell];

        let timer = Instant::now();
        let mut aos_c = aos_c_batch(l, gaussians, pxs, pys, pzs);
        time_ms_multiply += timer.elapsed().as_millis();

        let timer = Instant::now();
        if basis.cartesian_deg[ishell] == basis.spherical_deg[ishell] {
            aos.append(&mut aos_c);
        } else {
            let mut aos_s = transform_to_spherical(
                num_points,
                &aos_c,
                basis.spherical_deg[ishell],
                &c_to_s_matrices[l],
            );
            aos.append(&mut aos_s);
        }
        time_ms_transform += timer.elapsed().as_millis();

        offset += num_primitives;
    }

    println!("time spent in exp: {} ms", time_ms_gaussian);
    println!("time spent in multiply: {} ms", time_ms_multiply);
    println!("time spent in transform: {} ms", time_ms_transform);

    return aos;
}
