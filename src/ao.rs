use crate::basis::Basis;
use crate::generate;
use crate::point::Point;

fn get_s(p2: f64, basis: &Basis, offset: usize, num_primitives: usize) -> f64 {
    let mut s = 0.0;

    for ip in offset..(offset + num_primitives) {
        let e = basis.primitive_exponents[ip];
        let c = basis.contraction_coefficients[ip];
        s += c * (-e * p2).exp();
    }

    return s;
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

pub fn get_ao_noddy(
    points_bohr: Vec<Point>,
    basis: &Basis,
    c_to_s_matrices: &Vec<Vec<Vec<f64>>>,
) -> Vec<f64> {
    let num_points = points_bohr.len();

    let mut offset = 0;
    let mut aos = Vec::new();

    for ishell in 0..basis.num_shells {
        let (x, y, z) = basis.shell_centers_coordinates[ishell];

        let num_primitives = basis.shell_num_primitives[ishell];
        let l = basis.shell_l_quantum_numbers[ishell];

        let mut px = Vec::new();
        let mut py = Vec::new();
        let mut pz = Vec::new();
        let mut p2 = Vec::new();
        let mut s = Vec::new();

        for ipoint in 0..num_points {
            px.push(points_bohr[ipoint].x - x);
            py.push(points_bohr[ipoint].y - y);
            pz.push(points_bohr[ipoint].z - z);

            p2.push(px[ipoint] * px[ipoint] + py[ipoint] * py[ipoint] + pz[ipoint] * pz[ipoint]);

            s.push(get_s(p2[ipoint], &basis, offset, num_primitives));
        }

        let mut aos_c = Vec::new();
        for (i, j, k) in generate::get_ijk_list(l).iter() {
            for ipoint in 0..num_points {
                aos_c.push(
                    s[ipoint]
                        * px[ipoint].powi(*i as i32)
                        * py[ipoint].powi(*j as i32)
                        * pz[ipoint].powi(*k as i32),
                );
            }
        }

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

        offset += num_primitives;
    }

    return aos;
}
