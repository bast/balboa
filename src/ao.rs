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

pub fn get_ao_noddy(basis: &Basis, coordinates: (f64, f64, f64)) -> () {
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

        let ijk_list = generate::get_ijk_list(l);

        for (i, j, k) in ijk_list.iter() {
            let ao = s * px.powf(*i as f64) * py.powf(*j as f64) * pz.powf(*k as f64);
            dbg!(ao);
        }
    }
}
