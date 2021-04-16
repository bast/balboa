#![allow(clippy::redundant_field_names)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_range_loop)]

pub struct Basis {
    pub shell_centers_coordinates: Vec<(f64, f64, f64)>,
    pub shell_l_quantum_numbers: Vec<usize>,
    pub shell_num_primitives: Vec<usize>,
    pub primitive_exponents: Vec<f64>,
    pub contraction_coefficients: Vec<f64>,
    pub num_ao_cartesian: usize,
    pub num_ao_spherical: usize,
}

impl Basis {
    pub fn new(
        center_coordinates_bohr: Vec<(f64, f64, f64)>,
        shell_centers: Vec<usize>,
        shell_l_quantum_numbers: Vec<usize>,
        shell_num_primitives: Vec<usize>,
        primitive_exponents: Vec<f64>,
        contraction_coefficients: Vec<f64>,
    ) -> Basis {
        let shell_centers_coordinates = shell_centers
            .iter()
            .map(|&c| center_coordinates_bohr[c])
            .collect();

        // threshold and factors match Dalton implementation, see also pink book
        let f = vec![1.0, 1.3333, 1.6, 1.83, 2.03, 2.22, 2.39, 2.55, 2.70, 2.84];
        let shell_screening_threshold: f64 = 2.0e-12;

        let num_shells = shell_centers.len();

        let mut shell_extent_squared = vec![0.0; num_shells];
        let mut n = 0;
        for ishell in 0..num_shells {
            let mut r: f64 = 0.0;
            for _ in 0..shell_num_primitives[ishell] {
                let e = primitive_exponents[n];
                let c = contraction_coefficients[n];
                n += 1;
                r = r.max((c.abs().ln() - shell_screening_threshold.ln()) / e);
            }

            let l = shell_l_quantum_numbers[ishell];
            if l < 10 {
                r = r.sqrt() * f[l];
            } else {
                r = 1.0e10;
            }

            shell_extent_squared[ishell] = r * r;
        }

        let num_ao_cartesian = shell_l_quantum_numbers
            .iter()
            .fold(0, |sum, l| sum + ((l + 1) * (l + 2) / 2));
        let num_ao_spherical = shell_l_quantum_numbers
            .iter()
            .fold(0, |sum, l| sum + (2 * l + 1));

        Basis {
            shell_centers_coordinates: shell_centers_coordinates,
            shell_l_quantum_numbers: shell_l_quantum_numbers,
            shell_num_primitives: shell_num_primitives,
            primitive_exponents: primitive_exponents,
            contraction_coefficients: contraction_coefficients,
            num_ao_cartesian: num_ao_cartesian,
            num_ao_spherical: num_ao_spherical,
        }
    }
}
