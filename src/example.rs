use rand::Rng;

use crate::basis::{Basis, Shell};

pub fn example_basis(extra_l_coverage: bool) -> Basis {
    let center_coordinates_bohr = vec![(1.7, 0.0, 0.0), (0.0, 0.0, 0.0)];

    let l_quantum_numbers = if extra_l_coverage {
        // these values have been modified for better code coverage for higher l
        vec![0, 1, 2, 3, 4, 5, 0, 1, 2]
    } else {
        vec![0, 0, 0, 1, 1, 2, 0, 0, 1]
    };

    let shells = vec![
        Shell {
            coordinates: center_coordinates_bohr[0],
            l: l_quantum_numbers[0],
            primitive_exponents: vec![
                1.471e+04, 2.207e+03, 5.028e+02, 1.426e+02, 4.647e+01, 1.670e+01, 6.356e+00,
                1.316e+00, 3.897e-01,
            ],
            contraction_coefficients: vec![
                6.86365e-01,
                1.27435e+00,
                2.13913e+00,
                3.13055e+00,
                3.63823e+00,
                2.64148e+00,
                7.55357e-01,
                1.34270e-02,
                -8.19760e-04,
            ],
        },
        Shell {
            coordinates: center_coordinates_bohr[0],
            l: l_quantum_numbers[1],
            primitive_exponents: vec![
                1.471e+04, 2.207e+03, 5.028e+02, 1.426e+02, 4.647e+01, 1.670e+01, 6.356e+00,
                1.316e+00, 3.897e-01,
            ],
            contraction_coefficients: vec![
                -1.57074e-01,
                -3.00172e-01,
                -4.91514e-01,
                -7.84991e-01,
                -9.34756e-01,
                -1.00548e+00,
                -3.20466e-01,
                4.92853e-01,
                1.99941e-01,
            ],
        },
        Shell {
            coordinates: center_coordinates_bohr[0],
            l: l_quantum_numbers[2],
            primitive_exponents: vec![3.897e-01],
            contraction_coefficients: vec![3.51526e-01],
        },
        Shell {
            coordinates: center_coordinates_bohr[0],
            l: l_quantum_numbers[3],
            primitive_exponents: vec![2.267e+01, 4.977e+00, 1.347e+00, 3.471e-01],
            contraction_coefficients: vec![3.16438e+00, 2.49771e+00, 1.05186e+00, 1.73975e-01],
        },
        Shell {
            coordinates: center_coordinates_bohr[0],
            l: l_quantum_numbers[4],
            primitive_exponents: vec![3.471e-01],
            contraction_coefficients: vec![3.79759e-01],
        },
        Shell {
            coordinates: center_coordinates_bohr[0],
            l: l_quantum_numbers[5],
            primitive_exponents: vec![1.640e+00],
            contraction_coefficients: vec![6.77559e+00],
        },
        Shell {
            coordinates: center_coordinates_bohr[1],
            l: l_quantum_numbers[6],
            primitive_exponents: vec![1.301e+01, 1.962e+00, 4.446e-01, 1.220e-01],
            contraction_coefficients: vec![9.61066e-02, 1.63020e-01, 1.85545e-01, 7.37438e-02],
        },
        Shell {
            coordinates: center_coordinates_bohr[1],
            l: l_quantum_numbers[7],
            primitive_exponents: vec![1.220e-01],
            contraction_coefficients: vec![1.47123e-01],
        },
        Shell {
            coordinates: center_coordinates_bohr[1],
            l: l_quantum_numbers[8],
            primitive_exponents: vec![7.270e-01],
            contraction_coefficients: vec![9.56881e-01],
        },
    ];

    Basis::new(&shells)
}

pub fn example_basis_benchmark(num_centers: usize, l: usize) -> Basis {
    let rmin = -2.0;
    let rmax = 2.0;

    let mut rng = rand::thread_rng();

    let mut shells = Vec::new();
    for _ in 0..num_centers {
        shells.push(Shell {
            coordinates: (
                rng.gen_range(rmin..rmax),
                rng.gen_range(rmin..rmax),
                rng.gen_range(rmin..rmax),
            ),
            l: l,
            primitive_exponents: vec![
                1.471e+04, 2.207e+03, 5.028e+02, 1.426e+02, 4.647e+01, 1.670e+01, 6.356e+00,
                1.316e+00, 3.897e-01,
            ],
            contraction_coefficients: vec![
                6.86365e-01,
                1.27435e+00,
                2.13913e+00,
                3.13055e+00,
                3.63823e+00,
                2.64148e+00,
                7.55357e-01,
                1.34270e-02,
                -8.19760e-04,
            ],
        });
    }

    Basis::new(&shells)
}
