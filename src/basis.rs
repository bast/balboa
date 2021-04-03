#![allow(clippy::redundant_field_names)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_range_loop)]

use crate::limits;

pub struct Basis {
    _is_spherical: bool,
    _num_centers: usize,
    _center_coordinates_bohr: Vec<(f64, f64, f64)>,
    pub num_shells: usize,
    _shell_centers: Vec<usize>,
    pub shell_centers_coordinates: Vec<(f64, f64, f64)>,
    pub shell_l_quantum_numbers: Vec<usize>,
    pub shell_num_primitives: Vec<usize>,
    pub primitive_exponents: Vec<f64>,
    pub contraction_coefficients: Vec<f64>,
    _shell_off: Vec<usize>,
    _num_ao_cartesian: usize,
    _num_ao_spherical: usize,
    _num_ao: usize,
    _ao_centers: Vec<usize>,
    _geo_offset: Vec<usize>,
    _geo_offset_size: usize,
}

impl Basis {
    pub fn new(
        is_spherical: bool,
        num_centers: usize,
        center_coordinates_bohr: Vec<(f64, f64, f64)>,
        num_shells: usize,
        shell_centers: Vec<usize>,
        shell_l_quantum_numbers: Vec<usize>,
        shell_num_primitives: Vec<usize>,
        primitive_exponents: Vec<f64>,
        contraction_coefficients: Vec<f64>,
    ) -> Basis {
        for &l in shell_l_quantum_numbers.iter() {
            assert!(l <= limits::MAX_L_VALUE, "increase MAX_L_VALUE");
        }

        let mut shell_centers_coordinates = Vec::new();
        for icenter in shell_centers.iter() {
            let (x, y, z) = center_coordinates_bohr[icenter - 1];
            shell_centers_coordinates.push((x, y, z));
        }

        // get approximate spacial shell extent
        const SHELL_SCREENING_THRESHOLD: f64 = 2.0e-12;

        // threshold and factors match Dalton implementation, see also pink book
        let f = vec![1.0, 1.3333, 1.6, 1.83, 2.03, 2.22, 2.39, 2.55, 2.70, 2.84];

        let mut shell_extent_squared = vec![0.0; num_shells];
        let mut n = 0;

        for ishell in 0..num_shells {
            let mut r: f64 = 0.0;
            for _ in 0..shell_num_primitives[ishell] {
                let e = primitive_exponents[n];
                let c = contraction_coefficients[n];
                n += 1;
                r = r.max((c.abs().ln() - SHELL_SCREENING_THRESHOLD.ln()) / e);
            }

            if shell_l_quantum_numbers[ishell] < 10 {
                r = r.sqrt() * f[shell_l_quantum_numbers[ishell]];
            } else {
                r = 1.0e10;
            }

            shell_extent_squared[ishell] = r * r;
        }

        let mut shell_off = Vec::new();

        let mut num_ao_cartesian = 0;
        let mut num_ao_spherical = 0;

        for l in shell_l_quantum_numbers.iter() {
            let kc = (l + 1) * (l + 2) / 2;
            let ks = 2 * l + 1;

            if is_spherical {
                shell_off.push(num_ao_spherical);
            } else {
                shell_off.push(num_ao_cartesian);
            }

            num_ao_cartesian += kc;
            num_ao_spherical += ks;
        }

        let num_ao = if is_spherical {
            num_ao_spherical
        } else {
            num_ao_cartesian
        };

        let mut ao_centers = vec![0; num_ao];
        let mut i = 0;

        for ishell in 0..num_shells {
            let l = shell_l_quantum_numbers[ishell];
            let deg = if is_spherical {
                2 * l + 1
            } else {
                (l + 1) * (l + 2) / 2
            };

            for j in i..(i + deg) {
                ao_centers[j] = shell_centers[ishell] - 1;
            }

            i += deg;
        }

        let g = limits::MAX_GEO_DIFF_ORDER + 1;
        let array_length = g * g * g;
        let mut geo_offset = vec![0; array_length];

        let mut id = 0;
        let mut m = 0;
        for l in 0..=limits::MAX_GEO_DIFF_ORDER {
            for a in 1..(l + 2) {
                for b in 1..=a {
                    let i = l + 1 - a;
                    let j = a - b;
                    let k = b - 1;

                    id = (limits::MAX_GEO_DIFF_ORDER + 1) * (limits::MAX_GEO_DIFF_ORDER + 1) * k;
                    id += (limits::MAX_GEO_DIFF_ORDER + 1) * j;
                    id += i;

                    geo_offset[id] = m * num_ao;

                    m += 1;
                }
            }
        }
        let geo_offset_size = id;

        Basis {
            _is_spherical: is_spherical,
            _num_centers: num_centers,
            _center_coordinates_bohr: center_coordinates_bohr,
            num_shells: num_shells,
            _shell_centers: shell_centers,
            shell_centers_coordinates: shell_centers_coordinates,
            shell_l_quantum_numbers: shell_l_quantum_numbers,
            shell_num_primitives: shell_num_primitives,
            primitive_exponents: primitive_exponents,
            contraction_coefficients: contraction_coefficients,
            _shell_off: shell_off,
            _num_ao_cartesian: num_ao_cartesian,
            _num_ao_spherical: num_ao_spherical,
            _num_ao: num_ao,
            _ao_centers: ao_centers,
            _geo_offset: geo_offset,
            _geo_offset_size: geo_offset_size,
        }
    }
}
