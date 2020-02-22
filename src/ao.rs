use crate::basis::Basis;
use crate::generate;

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
    ao_c: &Vec<f64>,
    cartesian_deg: usize,
    spherical_deg: usize,
    c_to_s_matrix: &Vec<Vec<f64>>,
) -> Vec<f64> {
    let mut ao_s = vec![0.0; spherical_deg];

    for (i, row) in c_to_s_matrix.iter().enumerate() {
        for (j, element) in row.iter().enumerate() {
            ao_s[j] += element * ao_c[i];
        }
    }

    return ao_s;
}

pub fn get_ao_noddy(
    coordinates: (f64, f64, f64),
    basis: &Basis,
    c_to_s_matrices: &Vec<Vec<Vec<f64>>>,
) -> () {
    let mut offset = 0;

    for ishell in 0..basis.num_shells {
        let (x, y, z) = basis.shell_centers_coordinates[ishell];

        let px = coordinates.0 - x;
        let py = coordinates.1 - y;
        let pz = coordinates.2 - z;

        let p2 = px * px + py * py + pz * pz;

        let num_primitives = basis.shell_num_primitives[ishell];
        let s = get_s(p2, &basis, offset, num_primitives);
        offset += num_primitives;

        let l = basis.shell_l_quantum_numbers[ishell];

        let mut ao_c = Vec::new();
        for (i, j, k) in generate::get_ijk_list(l).iter() {
            ao_c.push(s * px.powf(*i as f64) * py.powf(*j as f64) * pz.powf(*k as f64));
        }

        if basis.cartesian_deg[ishell] == basis.spherical_deg[ishell] {
            dbg!(ao_c);
        } else {
            let ao_s = transform_to_spherical(
                &ao_c,
                basis.cartesian_deg[ishell],
                basis.spherical_deg[ishell],
                &c_to_s_matrices[l],
            );
            dbg!(ao_s);
        }
    }
}
