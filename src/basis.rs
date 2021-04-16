#[derive(Clone)]
pub struct Shell {
    pub coordinates: (f64, f64, f64),
    pub l: usize,
    pub primitive_exponents: Vec<f64>,
    pub contraction_coefficients: Vec<f64>,
}

pub struct Basis {
    pub num_ao_cartesian: usize,
    pub num_ao_spherical: usize,
    pub shells: Vec<Shell>,
}

impl Basis {
    pub fn new(shells: &[Shell]) -> Basis {
        let num_ao_cartesian = shells
            .iter()
            .fold(0, |sum, s| sum + ((s.l + 1) * (s.l + 2) / 2));
        let num_ao_spherical = shells.iter().fold(0, |sum, s| sum + (2 * s.l + 1));

        Basis {
            num_ao_cartesian,
            num_ao_spherical,
            shells: shells.to_vec(),
        }
    }
}

// // threshold and factors match Dalton implementation, see also pink book
// let f = vec![1.0, 1.3333, 1.6, 1.83, 2.03, 2.22, 2.39, 2.55, 2.70, 2.84];
// let shell_screening_threshold: f64 = 2.0e-12;

// let mut shell_extent_squared = vec![0.0; num_shells];
// let mut n = 0;
// for ishell in 0..num_shells {
//     let mut r: f64 = 0.0;
//     for _ in 0..shell_num_primitives[ishell] {
//         let e = primitive_exponents[n];
//         let c = contraction_coefficients[n];
//         n += 1;
//         r = r.max((c.abs().ln() - shell_screening_threshold.ln()) / e);
//     }

//     let l = shell_l_quantum_numbers[ishell];
//     if l < 10 {
//         r = r.sqrt() * f[l];
//     } else {
//         r = 1.0e10;
//     }

//     shell_extent_squared[ishell] = r * r;
// }
