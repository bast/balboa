def get_ijk_list(m):
    l = []
    for a in range(1, m + 2):
        for b in range(1, a + 1):
            i = m + 1 - a
            j = a - b
            k = b - 1
            l.append([i, j, k])
    return l


def test_get_ijk_list():
    assert get_ijk_list(3) == [[3, 0, 0], [2, 1, 0], [2, 0, 1], [1, 2, 0], [1, 1, 1], [1, 0, 2], [0, 3, 0], [0, 2, 1], [0, 1, 2], [0, 0, 3]]


def print_line(exp, geo, m, r):

    _exp = exp[:]
    _exp[m] -= 1

    vec_r = 'buffer[OFFSET_%02d_%02d_%02d_%i%i%i]' % (exp[0], exp[1], exp[2], geo[0], geo[1], geo[2])
    vec_p = '%s' % r
    vec_a = 'buffer[OFFSET_%02d_%02d_%02d_%i%i%i]' % (_exp[0], _exp[1], _exp[2], geo[0], geo[1], geo[2])

    if (geo[m] > 0):
        geo_right = geo[:]
        geo_right[m] -= 1
        vec_b = 'buffer[OFFSET_%02d_%02d_%02d_%i%i%i]' % (_exp[0], _exp[1], _exp[2], geo_right[0], geo_right[1], geo_right[2])
        if geo[m] > 1:
            return 'get_pa_plus_sb_block(&{0}, {1}, {2}.0, &{3}, &{4});'.format(vec_a, vec_p, geo[m], vec_b, vec_r)
        else:
            return 'get_pa_plus_b_block(&{0}, {1}, &{2}, &{3});'.format(vec_a, vec_p, vec_b, vec_r)
    else:
        return 'get_pa_block(&{0}, {1}, &{2});'.format(vec_a, vec_p, vec_r)


def test_print_line():
    s = print_line([1, 1, 1], [1, 1, 1], 2, 'py')
    assert s == 'get_pa_plus_b_block(&buffer[OFFSET_01_01_00_111], py, &buffer[OFFSET_01_01_00_110], &buffer[OFFSET_01_01_01_111]);'


def get_offsets(max_l_value, ao_chunk_length, max_geo_diff_order):
    s = []
    s.append('#ifndef OFFSETS_H_INCLUDED')
    s.append('#define OFFSETS_H_INCLUDED')
    s.append('')
    offset = 0
    for l in range(0, max_l_value + 1):
        for exp in get_ijk_list(l):
            for g in range(0, max_geo_diff_order + 1):
                for geo in get_ijk_list(g):
                    s_geo = '%i%i%i' % (geo[0], geo[1], geo[2])
                    s.append('#define OFFSET_%02d_%02d_%02d_%s %i' % (exp[0], exp[1], exp[2], s_geo, offset))
                    offset += ao_chunk_length
    s.append('')
    s.append('#define BUFFER_LENGTH {0}'.format(offset))
    s.append('')
    s.append('#endif // OFFSETS_H_INCLUDED')
    return '\n'.join(s)


def test_get_offsets():
    s = '''#ifndef OFFSETS_H_INCLUDED
#define OFFSETS_H_INCLUDED

#define OFFSET_00_00_00_000 0
#define OFFSET_00_00_00_100 32
#define OFFSET_00_00_00_010 64
#define OFFSET_00_00_00_001 96
#define OFFSET_01_00_00_000 128
#define OFFSET_01_00_00_100 160
#define OFFSET_01_00_00_010 192
#define OFFSET_01_00_00_001 224
#define OFFSET_00_01_00_000 256
#define OFFSET_00_01_00_100 288
#define OFFSET_00_01_00_010 320
#define OFFSET_00_01_00_001 352
#define OFFSET_00_00_01_000 384
#define OFFSET_00_00_01_100 416
#define OFFSET_00_00_01_010 448
#define OFFSET_00_00_01_001 480

#define BUFFER_LENGTH 512

#endif // OFFSETS_H_INCLUDED'''
    assert get_offsets(1, 32, 1) == s


def write_routine(_maxg, file_name, max_l_value, ao_chunk_length, max_geo_diff_order, d_geo_ijk_slice, suffix):
    import cs_trans

    cs = cs_trans.get_cs_trans(max_l_value)

    sfoo = '''
//  this file is automatically generated by generate.py

#include <iostream>
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <cstring>

#include "autogenerated.h"
#include "offsets.h"
#include "ao_vector.h"
          \n'''

    s = get_signature(_maxg, suffix) + '\n    )'

    sfoo += '''
%s
{
         \n''' % (s)

    if _maxg > 0:
        for i in range(_maxg + 1):
            sfoo += '    double fx_%i;\n' % i
            sfoo += '    double fy_%i;\n' % i
            sfoo += '    double fz_%i;\n' % i
    sfoo += '    double a;\n'
    sfoo += '    double c;\n\n'

    if suffix == 'block':
        sfoo += '        get_p2_block(shell_centers_coordinates,\n'
    else:
        sfoo += '        get_p2(num_points_batch, shell_centers_coordinates,\n'
    sfoo += '                     pw,\n'
    sfoo += '                     px,\n'
    sfoo += '                     py,\n'
    sfoo += '                     pz,\n'
    sfoo += '                     p2);\n\n'

    sfoo += '        // screening\n'
    if suffix == 'block':
        sfoo += '        if (not calculate_chunk_block(extent_squared, p2)) return;\n\n'
    else:
        sfoo += '        if (not calculate_chunk(num_points_batch, extent_squared, p2)) return;\n\n'

    array = 'buffer[OFFSET_00_00_00_000]'
    if suffix == 'block':
        sfoo += '        memset(&%s, 0, %i*sizeof(double));\n' % (array, ao_chunk_length)
    else:
        sfoo += '        memset(&%s, 0, num_points_batch*sizeof(double));\n' % (array)
    for g in range(1, _maxg + 1):
        for geo in get_ijk_list(g):
            array = 'buffer[OFFSET_00_00_00_%i%i%i]' % (geo[0], geo[1], geo[2])
            if suffix == 'block':
                sfoo += '        memset(&%s, 0, %i*sizeof(double));\n' % (array, ao_chunk_length)
            else:
                sfoo += '        memset(&%s, 0, num_points_batch*sizeof(double));\n' % (array)

    if suffix == 'block':
        sfoo += '''
            for (int i = 0; i < num_primitives; i++)
            {
                a = -primitive_exponents[i];
                c = contraction_coefficients[i];
       
                get_exp_block(p2, c, a, s);
       
                #pragma ivdep
                #pragma vector aligned
                for (int k = 0; k < %i; k++)
                {
                    buffer[OFFSET_00_00_00_000 + k] += s[k];
                  \n''' % ao_chunk_length
    else:
        sfoo += '''
            for (int i = 0; i < num_primitives; i++)
            {
                a = -primitive_exponents[i];
                c = contraction_coefficients[i];
       
                get_exp(num_points_batch, p2, c, a, s);
       
                #pragma ivdep
                #pragma vector aligned
                for (int k = 0; k < num_points_batch; k++)
                {
                    buffer[OFFSET_00_00_00_000 + k] += s[k];
                  \n'''

    if (_maxg > 0):
        sfoo += '''                fx_0 = 1.0;
                fy_0 = 1.0;
                fz_0 = 1.0;
                fx_1 = 2.0*a*px[k];
                fy_1 = 2.0*a*py[k];
                fz_1 = 2.0*a*pz[k];
                buffer[OFFSET_00_00_00_100 + k] += fx_1*s[k];
                buffer[OFFSET_00_00_00_010 + k] += fy_1*s[k];
                buffer[OFFSET_00_00_00_001 + k] += fz_1*s[k];
                  \n'''

    for g in range(2, _maxg + 1):
        sfoo += '                fx_%i = fx_%i*fx_1 + %i.0*a*fx_%i;\n' % (int(g), int(g - 1), int(g - 1) * 2, int(g - 2))
        sfoo += '                fy_%i = fy_%i*fy_1 + %i.0*a*fy_%i;\n' % (int(g), int(g - 1), int(g - 1) * 2, int(g - 2))
        sfoo += '                fz_%i = fz_%i*fz_1 + %i.0*a*fz_%i;\n' % (int(g), int(g - 1), int(g - 1) * 2, int(g - 2))
        for geo in get_ijk_list(g):
            sfoo += '                buffer[OFFSET_00_00_00_%i%i%i + k] += fx_%i*fy_%i*fz_%i*s[k];\n' \
                    % (geo[0], geo[1], geo[2],
                       geo[0], geo[1], geo[2])
    sfoo += '            }\n'
    sfoo += '        }\n'

    for l in range(0, max_l_value + 1):
        sfoo += '\n        if (shell_l_quantum_numbers == ' + '%i)\n' % l
        sfoo += '        {\n'
        if l < 2:
            c = 0
            for exp in get_ijk_list(l):
                for s in range(len(cs[l][c])):
                    f = cs[l][c][s]
                    if abs(f) > 0.0:
                        for g in range(0, _maxg + 1):
                            for geo in get_ijk_list(g):
                                s_geo = '%i%i%i' % (geo[0], geo[1], geo[2])
                                if suffix == 'block':
                                    sfoo += '            memcpy(&ao_000[%i*xoff + %i*num_points], &buffer[OFFSET_%02d_%02d_%02d_%s], %i*sizeof(double));\n' % (d_geo_ijk_slice[tuple(geo)], s, exp[0], exp[1], exp[2], s_geo, ao_chunk_length)
                                else:
                                    sfoo += '            memcpy(&ao_000[%i*xoff + %i*num_points], &buffer[OFFSET_%02d_%02d_%02d_%s], num_points_batch*sizeof(double));\n' % (d_geo_ijk_slice[tuple(geo)], s, exp[0], exp[1], exp[2], s_geo)
                c += 1
        else:
            sfoo += '            if (is_spherical)\n'
            sfoo += '            {\n'
            c = 0
            for exp in get_ijk_list(l):
                for s in range(len(cs[l][c])):
                    f = cs[l][c][s]
                    if abs(f) > 0.0:
                        for g in range(0, _maxg + 1):
                            for geo in get_ijk_list(g):
                                s_geo = '%i%i%i' % (geo[0], geo[1], geo[2])
                                if suffix == 'block':
                                    sfoo += '                vec_daxpy_block(%20.16e, &buffer[OFFSET_%02d_%02d_%02d_%s], &ao_000[%i*xoff + %i*num_points]);\n' % (f, exp[0], exp[1], exp[2], s_geo, d_geo_ijk_slice[tuple(geo)], s)
                                else:
                                    sfoo += '                vec_daxpy(num_points_batch, %20.16e, &buffer[OFFSET_%02d_%02d_%02d_%s], &ao_000[%i*xoff + %i*num_points]);\n' % (f, exp[0], exp[1], exp[2], s_geo, d_geo_ijk_slice[tuple(geo)], s)
                c += 1
            sfoo += '            }\n'
            sfoo += '            else\n'
            sfoo += '            {\n'
            s = 0
            for exp in get_ijk_list(l):
                for g in range(0, _maxg + 1):
                    for geo in get_ijk_list(g):
                        s_geo = '%i%i%i' % (geo[0], geo[1], geo[2])
                        if suffix == 'block':
                            sfoo += '                memcpy(&ao_000[%i*xoff + %i*num_points], &buffer[OFFSET_%02d_%02d_%02d_%s], %i*sizeof(double));\n' % (d_geo_ijk_slice[tuple(geo)], s, exp[0], exp[1], exp[2], s_geo, ao_chunk_length)
                        else:
                            sfoo += '                memcpy(&ao_000[%i*xoff + %i*num_points], &buffer[OFFSET_%02d_%02d_%02d_%s], num_points_batch*sizeof(double));\n' % (d_geo_ijk_slice[tuple(geo)], s, exp[0], exp[1], exp[2], s_geo)
                s += 1
            sfoo += '            }\n'
        sfoo += '            return;\n'
        sfoo += '        }\n'
        sfoo += '        else\n'
        sfoo += '        {\n'
        if l + 1 < max_l_value + 1:
            for exp in get_ijk_list(l + 1):
                for g in range(0, _maxg + 1):
                    for geo in get_ijk_list(g):
                        if exp[0] > 0:
                            sfoo += print_line(exp, geo, 0, 'px') + '\n'
                        else:
                            if exp[1] > 0:
                                sfoo += print_line(exp, geo, 1, 'py') + '\n'
                            else:
                                if exp[2] > 0:
                                    sfoo += print_line(exp, geo, 2, 'pz') + '\n'
        else:
            sfoo += '             std::cout << "error: order too high";\n'
            sfoo += '             exit(1);\n'
        sfoo += '        }\n'
    sfoo += '\n}\n'

    with open(file_name, 'w') as f:
        f.write(sfoo)


def get_signature(g, suffix):

    s = []

    s.append('void get_ao_g{0}_{1}('.format(g, suffix))
    s.append('    const int    shell_l_quantum_numbers,')
    s.append('    const int    num_primitives,')
    s.append('    const bool   is_spherical,')
    s.append('    const double primitive_exponents[],')
    s.append('    const double contraction_coefficients[],')
    s.append('    const int    num_points,')
    if suffix == 'explicit':
        s.append('    const int    num_points_batch,')
    s.append('    const int    xoff,')
    s.append('          double s[],')
    s.append('          double buffer[],')
    s.append('    const double shell_centers_coordinates[],')
    s.append('    const double extent_squared,')
    s.append('    const double pw[],')
    s.append('          double px[],')
    s.append('          double py[],')
    s.append('          double pz[],')
    s.append('          double p2[],')
    s.append('          double ao_000[]')

    return '\n'.join(s)


def get_header(max_geo_diff_order):

    s = []

    s.append('#ifndef AUTOGENERATED_H_INCLUDED')
    s.append('#define AUTOGENERATED_H_INCLUDED')
    s.append('')

    for g in range(0, max_geo_diff_order + 1):
        s.append(get_signature(g, 'block'))
        s.append('    );')
        s.append('')
        s.append(get_signature(g, 'explicit'))
        s.append('    );')
        s.append('')

    s.append('')
    s.append('#endif // AUTOGENERATED_H_INCLUDED')

    return '\n'.join(s)


def generate_dispatcher(output_directory, max_geo_diff_order):
    import os
    with open(os.path.join(output_directory, 'ao_dispatch.h'), 'w') as f:
        f.write('''#ifndef AO_DISPATCH_H_INCLUDED
#define AO_DISPATCH_H_INCLUDED
void ao_dispatch(
    const int    max_geo_order,
    const int    shell_l_quantum_number,
    const int    num_primitives,
    const bool   is_spherical,
    const double primitive_exponents[],
    const double contraction_coefficients[],
    const int    num_points,
    const int    num_points_batch,
    const int    xoff,
          double s[],
          double buffer[],
    const double shell_centers_coordinates[],
    const double extent_squared,
    const double pw[],
          double px[],
          double py[],
          double pz[],
          double p2[],
          double ao_000[]
    );
#endif''')
    with open(os.path.join(output_directory, 'ao_dispatch.cpp'), 'w') as f:
        f.write('''#include <iostream>
#include <stdlib.h>     /* exit, EXIT_FAILURE */


#include "autogenerated.h"
#include "parameters.h"


void ao_dispatch(
    const int    max_geo_order,
    const int    shell_l_quantum_number,
    const int    num_primitives,
    const bool   is_spherical,
    const double primitive_exponents[],
    const double contraction_coefficients[],
    const int    num_points,
    const int    num_points_batch,
    const int    xoff,
          double s[],
          double buffer[],
    const double shell_centers_coordinates[],
    const double extent_squared,
    const double pw[],
          double px[],
          double py[],
          double pz[],
          double p2[],
          double ao_000[]
    )
{
    if (num_points_batch < AO_CHUNK_LENGTH)
    {
        switch (max_geo_order)
        {
''')

        for g in range(max_geo_diff_order + 1):
            f.write('''            case {order}:
                get_ao_g{order}_explicit(
                    shell_l_quantum_number,
                    num_primitives,
                    is_spherical,
                    primitive_exponents,
                    contraction_coefficients,
                    num_points,
                    num_points_batch,
                    xoff,
                    s,
                    buffer,
                    shell_centers_coordinates,
                    extent_squared,
                    pw,
                    px,
                    py,
                    pz,
                    p2,
                    ao_000
                    );
                break;
'''.format(order=g))

        f.write('''            default:
                std::cout << "ERROR: get_ao order too high";
                exit(1);
                break;
        }
    }
    else
    {
        switch (max_geo_order)
        {
''')

        for g in range(max_geo_diff_order + 1):
            f.write('''            case {order}:
                get_ao_g{order}_block(
                    shell_l_quantum_number,
                    num_primitives,
                    is_spherical,
                    primitive_exponents,
                    contraction_coefficients,
                    num_points,
                    xoff,
                    s,
                    buffer,
                    shell_centers_coordinates,
                    extent_squared,
                    pw,
                    px,
                    py,
                    pz,
                    p2,
                    ao_000
                    );
                break;
'''.format(order=g))

        f.write('''            default:
                std::cout << "ERROR: get_ao order too high";
                exit(1);
                break;
        }
    }
}
''')


def main(output_directory, max_l_value, ao_chunk_length, max_geo_diff_order):
    import os

    d_geo_ijk_slice = {}
    j = 0
    for _g in range(0, max_geo_diff_order + 1):
        for geo in get_ijk_list(_g):
            d_geo_ijk_slice[tuple(geo)] = j
            j += 1

    generate_dispatcher(output_directory, max_geo_diff_order)

    with open(os.path.join(output_directory, 'offsets.h'), 'w') as f:
        f.write(get_offsets(max_l_value, ao_chunk_length, max_geo_diff_order))

    for suffix in ['block', 'explicit']:
        for g in range(0, max_geo_diff_order + 1):
            write_routine(g,
                          os.path.join(output_directory, 'autogenerated_{0}_{1}.cpp'.format(g, suffix)),
                          max_l_value,
                          ao_chunk_length,
                          max_geo_diff_order,
                          d_geo_ijk_slice,
                          suffix)

    with open(os.path.join(output_directory, 'autogenerated.h'), 'w') as f:
        f.write(get_header(max_geo_diff_order))


if __name__ == '__main__':
    import sys
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
