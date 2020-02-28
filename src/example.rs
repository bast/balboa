use crate::basis;
use crate::basis::Basis;

pub fn example_basis() -> Basis {
    let is_spherical = false;
    let num_centers = 2;
    let center_coordinates_bohr = vec![(1.7, 0.0, 0.0), (0.0, 0.0, 0.0)];
    let num_shells = 9;
    let shell_centers = vec![1, 1, 1, 1, 1, 1, 2, 2, 2];

    // these values have been modified for better code coverage for higher l
    let shell_l_quantum_numbers = vec![0, 1, 2, 3, 4, 5, 0, 1, 2];

    let shell_num_primitives = vec![9, 9, 1, 4, 1, 1, 4, 1, 1];

    let primitive_exponents = vec![
        1.471e+04, 2.207e+03, 5.028e+02, 1.426e+02, 4.647e+01, 1.670e+01, 6.356e+00, 1.316e+00,
        3.897e-01, 1.471e+04, 2.207e+03, 5.028e+02, 1.426e+02, 4.647e+01, 1.670e+01, 6.356e+00,
        1.316e+00, 3.897e-01, 3.897e-01, 2.267e+01, 4.977e+00, 1.347e+00, 3.471e-01, 3.471e-01,
        1.640e+00, 1.301e+01, 1.962e+00, 4.446e-01, 1.220e-01, 1.220e-01, 7.270e-01,
    ];

    let contraction_coefficients = vec![
        6.86365e-01,
        1.27435e+00,
        2.13913e+00,
        3.13055e+00,
        3.63823e+00,
        2.64148e+00,
        7.55357e-01,
        1.34270e-02,
        -8.19760e-04,
        -1.57074e-01,
        -3.00172e-01,
        -4.91514e-01,
        -7.84991e-01,
        -9.34756e-01,
        -1.00548e+00,
        -3.20466e-01,
        4.92853e-01,
        1.99941e-01,
        3.51526e-01,
        3.16438e+00,
        2.49771e+00,
        1.05186e+00,
        1.73975e-01,
        3.79759e-01,
        6.77559e+00,
        9.61066e-02,
        1.63020e-01,
        1.85545e-01,
        7.37438e-02,
        1.47123e-01,
        9.56881e-01,
    ];

    let basis = basis::initialize_basis(
        is_spherical,
        num_centers,
        center_coordinates_bohr,
        num_shells,
        shell_centers,
        shell_l_quantum_numbers,
        shell_num_primitives,
        primitive_exponents,
        contraction_coefficients,
    );

    return basis;
}
