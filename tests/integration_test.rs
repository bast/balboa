use balboa;

fn floats_are_same(f1: f64, f2: f64) -> bool {
    let d = f1 - f2;
    return d.abs() < 1.0e-7;
}

#[test]
fn undifferentiated() {
    let is_spherical = false;
    let num_centers = 2;
    let center_coordinates_bohr = vec![(1.7, 0.0, 0.0), (0.0, 0.0, 0.0)];
    let num_shells = 9;
    let shell_centers = vec![1, 1, 1, 1, 1, 1, 2, 2, 2];
    let shell_l_quantum_numbers = vec![0, 0, 0, 1, 1, 2, 0, 0, 1];
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

    let basis = balboa::initialize_basis(
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

    let c_to_s_matrices = balboa::cartesian_to_spherical_matrices();

    let coordinates = (-1.46254302355, 1.38973494775, 1.05509847591);

    let aos = balboa::get_ao_noddy(coordinates, &basis, &c_to_s_matrices);

    let aos_reference = vec![
        -5.077332474866332e-06,
        0.001238502645354896,
        0.002177441460600934,
        -0.005941367297202692,
        0.0026108501003325557,
        0.00198217938330633,
        -0.012968864636764336,
        0.005698984736061902,
        0.0043267172053119765,
        -1.5198424629213896e-08,
        5.070550674917344e-09,
        -9.689594211058621e-09,
        -1.1538771971215681e-08,
        1.3953740808523413e-08,
        0.057703726467940794,
        0.07816846783547603,
        -0.032309396944527016,
        0.030700975869788635,
        0.023308439426962427,
    ];

    for (&ao, &ao_reference) in aos.iter().zip(aos_reference.iter()) {
        assert!(floats_are_same(ao, ao_reference));
    }
}
