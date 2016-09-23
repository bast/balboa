import sys
import os


def test_balboa():
    import balboa
    import numpy as np
    from cffi import FFI

    num_centers = 2
    center_coordinates = [
        1.7,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ]

    num_shells = 9
    shell_centers = [1, 1, 1, 1, 1, 1, 2, 2, 2]
    shell_l_quantum_numbers = [0, 0, 0, 1, 1, 2, 0, 0, 1]
    shell_num_primitives = [9, 9, 1, 4, 1, 1, 4, 1, 1]

    primitive_exponents = [
        1.471e+04,
        2.207e+03,
        5.028e+02,
        1.426e+02,
        4.647e+01,
        1.670e+01,
        6.356e+00,
        1.316e+00,
        3.897e-01,
        1.471e+04,
        2.207e+03,
        5.028e+02,
        1.426e+02,
        4.647e+01,
        1.670e+01,
        6.356e+00,
        1.316e+00,
        3.897e-01,
        3.897e-01,
        2.267e+01,
        4.977e+00,
        1.347e+00,
        3.471e-01,
        3.471e-01,
        1.640e+00,
        1.301e+01,
        1.962e+00,
        4.446e-01,
        1.220e-01,
        1.220e-01,
        7.270e-01,
    ]

    contraction_coefficients = [
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
    ]

    context = balboa.new_context()

    ierr = balboa.set_basis(context,
                            0,
                            num_centers,
                            center_coordinates,
                            num_shells,
                            shell_centers,
                            shell_l_quantum_numbers,
                            shell_num_primitives,
                            primitive_exponents,
                            contraction_coefficients)

    dir_path = os.path.dirname(os.path.realpath(__file__))

    with open(os.path.join(dir_path, 'grid.txt'), 'r') as f:
        num_points = 0
        p = []
        for line in f.readlines():
            num_points += 1
            for xyzw in line.split():
                p.append(float(xyzw))

    with open(os.path.join(dir_path, 'result.txt'), 'r') as f:
        ref_aos = []
        for line in f.readlines():
            ref_aos.append(float(line))

    max_geo_order = 1

    # buffer length is adjusted for number of cartesian aos
    # so possibly longer than the number of spherical aos
    l = balboa.get_buffer_len(context, max_geo_order, num_points)

    # allocate a numpy array of length l and zero it out
    aos = np.zeros(l, dtype=np.float64)

    # cast a pointer which points to the numpy array data
    ffi = FFI()
    aos_p = ffi.cast("double *", aos.ctypes.data)

    ierr = balboa.get_ao(context,
                         max_geo_order,
                         num_points,
                         p,
                         aos_p)

#   with open(os.path.join(dir_path, 'result.txt'), 'w') as f:
#       k = 0
#       for _diff in [0, 1, 2, 3]:
#           for _ao in range(20):
#               for _point in range(128):
#                   if _point > 0:
#                       f.write("{0}\n".format(ref_aos[k]))
#                   k += 1

    k = 0
    kr = 0
    for _diff in [0, 1, 2, 3]:
        for _ao in range(20):
            for _point in range(128):
                if _point > 0:
                    error = aos[k] - ref_aos[kr]
                    if abs(ref_aos[kr]) > 1.0e-20:
                        error /= ref_aos[kr]
                    assert abs(error) < 1.0e-14
                    kr += 1
                k += 1

    balboa.free_context(context)
