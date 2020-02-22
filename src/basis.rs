const MAX_L_VALUE: usize = 5;
const MAX_GEO_DIFF_ORDER: usize = 3;

struct Basis {
    is_spherical: bool,
    num_centers: usize,
    center_coordinates_bohr: Vec<(f64, f64, f64)>,
    num_shells: usize,
    shell_centers: Vec<usize>,
    shell_centers_coordinates: Vec<(f64, f64, f64)>,
    shell_l_quantum_numbers: Vec<usize>,
    shell_num_primitives: Vec<usize>,
    primitive_exponents: Vec<f64>,
    contraction_coefficients: Vec<f64>,
    cartesian_deg: Vec<usize>,
    shell_off: Vec<usize>,
    spherical_deg: Vec<usize>,
    num_ao_cartesian: usize,
    num_ao_spherical: usize,
    num_ao: usize,
    ao_centers: Vec<usize>,
    geo_offset: Vec<usize>,
    geo_offset_size: usize,
}

fn set_basis(
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
    for ishell in 0..num_shells {
        assert!(
            shell_l_quantum_numbers[ishell] <= MAX_L_VALUE,
            "increase MAX_L_VALUE"
        );
    }

    let mut shell_centers_coordinates = Vec::new();
    for ishell in 0..num_shells {
        let i = shell_centers[ishell];
        let (x, y, z) = center_coordinates_bohr[i - 1];
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
        for j in 0..shell_num_primitives[ishell] {
            let e = primitive_exponents[n];
            let c = contraction_coefficients[n];
            n += 1;
            r = r.max((c.abs().ln() - SHELL_SCREENING_THRESHOLD.ln()) / e);
        }

        if shell_l_quantum_numbers[ishell] < 10 {
            r = r.powf(0.5 as f64) * f[shell_l_quantum_numbers[ishell]];
        } else {
            r = 1.0e10;
        }

        shell_extent_squared[ishell] = r * r;
    }

    let mut cartesian_deg = Vec::new();
    let mut shell_off = Vec::new();
    let mut spherical_deg = Vec::new();

    let mut num_ao_cartesian = 0;
    let mut num_ao_spherical = 0;

    for ishell in 0..num_shells {
        let l = shell_l_quantum_numbers[ishell];

        let kc = (l + 1) * (l + 2) / 2;
        let ks = 2 * l + 1;

        cartesian_deg.push(kc);
        spherical_deg.push(ks);

        if is_spherical {
            shell_off.push(num_ao_spherical);
        } else {
            shell_off.push(num_ao_cartesian);
        }

        num_ao_cartesian += kc;
        num_ao_spherical += ks;
    }

    let num_ao = match is_spherical {
        true => num_ao_spherical,
        false => num_ao_cartesian,
    };

    let mut ao_centers = vec![0; num_ao];
    let mut i = 0;

    for ishell in 0..num_shells {
        let deg = match is_spherical {
            true => spherical_deg[ishell],
            false => cartesian_deg[ishell],
        };

        for j in i..(i + deg) {
            ao_centers[j] = shell_centers[ishell] - 1;
        }

        i += deg;
    }

    let g = MAX_GEO_DIFF_ORDER + 1;
    let array_length = g * g * g;
    let mut geo_offset = vec![0; array_length];

    let mut id = 0;
    let mut m = 0;
    for l in 0..(MAX_GEO_DIFF_ORDER + 1) {
        for a in 1..(l + 2) {
            for b in 1..(a + 1) {
                let i = l + 1 - a;
                let j = a - b;
                let k = b - 1;

                id = (MAX_GEO_DIFF_ORDER + 1) * (MAX_GEO_DIFF_ORDER + 1) * k;
                id += (MAX_GEO_DIFF_ORDER + 1) * j;
                id += i;

                geo_offset[id] = m * num_ao;

                m += 1;
            }
        }
    }
    let geo_offset_size = id;

    return Basis {
        is_spherical: is_spherical,
        num_centers: num_centers,
        center_coordinates_bohr: center_coordinates_bohr,
        num_shells: num_shells,
        shell_centers: shell_centers,
        shell_centers_coordinates: shell_centers_coordinates,
        shell_l_quantum_numbers: shell_l_quantum_numbers,
        shell_num_primitives: shell_num_primitives,
        primitive_exponents: primitive_exponents,
        contraction_coefficients: contraction_coefficients,
        cartesian_deg: cartesian_deg,
        shell_off: shell_off,
        spherical_deg: spherical_deg,
        num_ao_cartesian: num_ao_cartesian,
        num_ao_spherical: num_ao_spherical,
        num_ao: num_ao,
        ao_centers: ao_centers,
        geo_offset: geo_offset,
        geo_offset_size: geo_offset_size,
    };
}
