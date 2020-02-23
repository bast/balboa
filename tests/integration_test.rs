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
        -0.0039168179007153464,
        0.0017211904091305233,
        0.001306742253524454,
        -0.009570065480906397,
        0.003192798146324918,
        -0.006101293610720681,
        -0.007265674306896439,
        0.008786319404741609,
        0.01496151378446815,
        -0.008711869795989362,
        -0.003087876706790655,
        -0.008591185010291626,
        0.00702691038492907,
        0.007998406154357254,
        -0.005102864786581584,
        -0.04198883331629265,
        0.034457525690197346,
        0.010331300718636494,
        -0.014541883483149953,
        0.005917105674812406,
        0.03309216065507327,
        -0.009485212840129847,
        -0.011752293050705426,
        -0.003592067180560306,
        3.407953753083849e-08,
        -3.7358672986733e-08,
        -6.947245123624579e-09,
        2.5942174680782685e-08,
        -5.37297588691292e-10,
        2.3952840752604505e-08,
        1.2226984314073409e-09,
        -2.3817625203685554e-08,
        2.369469628260252e-09,
        -3.19596551145082e-09,
        1.8590069026802844e-08,
        0.057703726467940794,
        -0.11432474729436803,
        0.10863345156303283,
        0.08247543127743061,
        -0.044901498074536264,
        0.03239255284916367,
        -0.011759096264421592,
        -0.03408959547374166,
        0.0022938320000254787,
    ];

    for (&ao, &ao_reference) in aos.iter().zip(aos_reference.iter()) {
        assert!(floats_are_same(ao, ao_reference));
    }
}
