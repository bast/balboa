#include <algorithm>
#include <assert.h>
#include <iostream>
#include <math.h>
#include <stdio.h>

#include "Main.h"
#include "autogenerated.h"
#include "balboa.h"
#include "offsets.h"
#include "parameters.h"

#define AS_TYPE(Type, Obj) reinterpret_cast<Type *>(Obj)
#define AS_CTYPE(Type, Obj) reinterpret_cast<const Type *>(Obj)


context_t *new_context()
{
    return AS_TYPE(context_t, new Main());
}
Main::Main()
{
    nullify();
}


void free_context(context_t *context)
{
    if (!context) return;
    delete AS_TYPE(Main, context);
}
Main::~Main()
{
    deallocate();
    nullify();
}


BALBOA_API int get_buffer_len(const context_t *context,
                              const int max_geo_order,
                              const int num_points)
{
    return AS_CTYPE(Main, context)->get_buffer_len(max_geo_order,
                                                   num_points);
}
int Main::get_buffer_len(const int max_geo_order,
                         const int num_points) const
{
    int m = 0;
    for (int l = 0; l <= max_geo_order; l++)
    {
        for (int a = 1; a <= (l + 1); a++)
        {
            for (int b = 1; b <= a; b++)
            {
                m++;
            }
        }
    }
    return m*num_points*num_ao_cartesian;
}


BALBOA_API int set_basis(context_t *context,
                         const int    basis_type,
                         const int    num_centers,
                         const double center_coordinates[],
                         const int    num_shells,
                         const int    shell_centers[],
                         const int    shell_l_quantum_numbers[],
                         const int    shell_num_primitives[],
                         const double primitive_exponents[],
                         const double contraction_coefficients[])
{
    return AS_TYPE(Main, context)->set_basis(basis_type,
                                             num_centers,
                                             center_coordinates,
                                             num_shells,
                                             shell_centers,
                                             shell_l_quantum_numbers,
                                             shell_num_primitives,
                                             primitive_exponents,
                                             contraction_coefficients);
}
int Main::set_basis(const int    in_basis_type,
                    const int    in_num_centers,
                    const double in_center_coordinates[],
                    const int    in_num_shells,
                    const int    in_shell_centers[],
                    const int    in_shell_l_quantum_numbers[],
                    const int    in_shell_num_primitives[],
                    const double in_primitive_exponents[],
                    const double in_contraction_coefficients[])
{
    int i, l, deg, kc, ks;

    // FIXME ugly hack to make basis set initialization idempotent
    if (is_initialized == 12345678) deallocate();
    nullify();

    num_centers = in_num_centers;

    switch (in_basis_type)
    {
        case 0:  // FIXME change to named integer parameters
            is_spherical = true;
            break;
        case 1:
            is_spherical = false;
            fprintf(stderr, "ERROR: XCINT_BASIS_CARTESIAN not tested.\n");
            exit(-1);
            break;
        default:
            fprintf(stderr, "ERROR: basis_type not recognized.\n");
            exit(-1);
            break;
    }

    num_shells = in_num_shells;

    center_coordinates = new double[3*num_centers];
    std::copy(&in_center_coordinates[0], &in_center_coordinates[3*num_centers], &center_coordinates[0]);

    shell_centers = new int[num_shells];
    std::copy(&in_shell_centers[0], &in_shell_centers[num_shells], &shell_centers[0]);

    shell_centers_coordinates = new double[3*num_shells];

    for (int ishell = 0; ishell < num_shells; ishell++)
    {
        for (int ixyz = 0; ixyz < 3; ixyz++)
        {
            shell_centers_coordinates[3*ishell + ixyz] = in_center_coordinates[3*(in_shell_centers[ishell]-1) + ixyz];
        }
    }

    shell_l_quantum_numbers = new int[num_shells];
    std::copy(&in_shell_l_quantum_numbers[0], &in_shell_l_quantum_numbers[num_shells], &shell_l_quantum_numbers[0]);

    for (int ishell = 0; ishell < num_shells; ishell++)
    {
        if (shell_l_quantum_numbers[ishell] > MAX_L_VALUE)
        {
            fprintf(stderr, "ERROR: increase MAX_L_VALUE.\n");
            exit(-1);
        }
    }

    shell_num_primitives = new int[num_shells];
    std::copy(&in_shell_num_primitives[0], &in_shell_num_primitives[num_shells], &shell_num_primitives[0]);

    int n = 0;
    for (int ishell = 0; ishell < num_shells; ishell++)
    {
        n += shell_num_primitives[ishell];
    }

    primitive_exponents = new double[n];
    contraction_coefficients = new double[n];
    std::copy(&in_primitive_exponents[0], &in_primitive_exponents[n], &primitive_exponents[0]);
    std::copy(&in_contraction_coefficients[0], &in_contraction_coefficients[n], &contraction_coefficients[0]);

    // get approximate spacial shell extent
    double SHELL_SCREENING_THRESHOLD = 2.0e-12;
    double e, c, r, r_temp;
    // threshold and factors match Dalton implementation, see also pink book
    double f[10] = {1.0, 1.3333, 1.6, 1.83, 2.03, 2.22, 2.39, 2.55, 2.70, 2.84};
    shell_extent_squared = new double[num_shells];
    n = 0;
    for (int ishell = 0; ishell < num_shells; ishell++)
    {
        r = 0.0;
        for (int j = 0; j < shell_num_primitives[ishell]; j++)
        {
            e = primitive_exponents[n];
            c = contraction_coefficients[n];
            n++;
            r_temp = (log(fabs(c)) - log(SHELL_SCREENING_THRESHOLD))/e;
            if (r_temp > r) r = r_temp;
        }
        if (shell_l_quantum_numbers[ishell] < 10)
        {
            r = pow(r, 0.5)*f[shell_l_quantum_numbers[ishell]];
        }
        else
        {
            r = 1.0e10;
        }
        shell_extent_squared[ishell] = r*r;
    }

    cartesian_deg = new int[num_shells];
    shell_off = new int[num_shells];
    spherical_deg = new int[num_shells];

    num_ao_cartesian = 0;
    num_ao_spherical = 0;
    for (int ishell = 0; ishell < num_shells; ishell++)
    {
        l = shell_l_quantum_numbers[ishell];
        kc = (l+1)*(l+2)/2;
        ks = 2*l + 1;
        cartesian_deg[ishell] = kc;
        spherical_deg[ishell] = ks;

        if (is_spherical)
        {
            shell_off[ishell] = num_ao_spherical;
        }
        else
        {
            shell_off[ishell] = num_ao_cartesian;
        }

        num_ao_cartesian += kc;
        num_ao_spherical += ks;
    }

    if (is_spherical)
    {
        num_ao = num_ao_spherical;
    }
    else
    {
        fprintf(stderr, "ERROR: XCint probably broken for cart basis, needs testing.\n");
        exit(-1);
        num_ao = num_ao_cartesian;
    }

    ao_center = new int[num_ao];
    i = 0;
    for (int ishell = 0; ishell < num_shells; ishell++)
    {
       if (is_spherical)
       {
           deg = spherical_deg[ishell];
       }
       else
       {
           deg = cartesian_deg[ishell];
       }
       for (int j = i; j < (i + deg); j++)
       {
           ao_center[j] = in_shell_centers[ishell] - 1;
       }

       i += deg;
    }

    is_initialized = 12345678;

    return 0;
}


BALBOA_API int get_ao(const context_t *context,
                      const int    max_geo_order,
                      const int    num_points,
                      const double p[],
                            double buffer[])
{
    return AS_CTYPE(Main, context)->get_ao(max_geo_order,
                                           num_points,
                                           p,
                                           buffer);
}
int Main::get_ao(const int    max_geo_order,
                 const int    num_points,
                 const double p[],
                       double ao_local[]) const
{
    double px[AO_CHUNK_LENGTH];
    double py[AO_CHUNK_LENGTH];
    double pz[AO_CHUNK_LENGTH];
    double p2[AO_CHUNK_LENGTH];
    double s[AO_CHUNK_LENGTH];
    double buffer[BUFFER_LENGTH];

    // FIXME can be optimized
    // we do this because p can be shorter than 4*AO_CHUNK_LENGTH
    // we pad it by very large numbers to let the code screen them away
    double p_block[4*AO_CHUNK_LENGTH];

    assert(max_geo_order <= MAX_GEO_DIFF_ORDER);

    int n = 0;
    for (int ishell = 0; ishell < num_shells; ishell++)
    {
        int num_points_left = num_points;

        for (int koff = 0; koff < num_points; koff += AO_CHUNK_LENGTH)
        {
            int num_points_batch = std::min(AO_CHUNK_LENGTH, num_points_left);

            num_points_left -= num_points_batch;

            std::fill(&p_block[0], &p_block[4*AO_CHUNK_LENGTH], 1.0e50);
            std::copy(&p[4*koff], &p[4*koff + 4*num_points_batch], &p_block[0]);

            int zoff = koff + shell_off[ishell]*num_points;
            int xoff = num_ao*num_points;

            switch (max_geo_order)
            {
                #include "aocalls.h"
                default:
                    std::cout << "ERROR: get_ao order too high\n";
                    exit(1);
                    break;
            }
        }
        n += shell_num_primitives[ishell];
    }

    return 0;
}


void Main::nullify()
{
    num_centers               = -1;
    num_shells                = -1;
    shell_l_quantum_numbers   = NULL;
    center_coordinates        = NULL;
    shell_centers             = NULL;
    shell_centers_coordinates = NULL;
    shell_extent_squared      = NULL;
    cartesian_deg             = NULL;
    shell_off                 = NULL;
    spherical_deg             = NULL;
    is_spherical              = false;
    num_ao                    = -1;
    num_ao_cartesian          = -1;
    num_ao_spherical          = -1;
    ao_center                 = NULL;
    shell_num_primitives      = NULL;
    primitive_exponents       = NULL;
    contraction_coefficients  = NULL;
    is_initialized            = 0;
}


void Main::deallocate()
{
    delete[] shell_l_quantum_numbers;
    delete[] center_coordinates;
    delete[] shell_centers;
    delete[] shell_centers_coordinates;
    delete[] shell_extent_squared;
    delete[] cartesian_deg;
    delete[] shell_off;
    delete[] spherical_deg;
    delete[] ao_center;
    delete[] shell_num_primitives;
    delete[] primitive_exponents;
    delete[] contraction_coefficients;
}
