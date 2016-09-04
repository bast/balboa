#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

#include <stdlib.h>
#include <vector>

#include "Basis.h"

class Main
{
    public:

        Main();
        ~Main();

        int set_basis(const int    basis_type,
                      const int    num_centers,
                      const double center_coordinates[],
                      const int    num_shells,
                      const int    shell_centers[],
                      const int    shell_l_quantum_numbers[],
                      const int    shell_num_primitives[],
                      const double primitive_exponents[],
                      const double contraction_coefficients[]);

        double *get_ao(const bool   use_gradient,
                   const int    max_ao_geo_order,
                   const int    block_length,
                   const double p[]);
               //  const double p[]) const;

        double *ao;

    private:

        Main(const Main &rhs);            // not implemented
        Main &operator=(const Main &rhs); // not implemented

        Basis basis;

        void get_ao(const Basis &basis,
                    const bool   use_gradient,
                    const int    max_ao_geo_order,
                    const int    block_length,
                    const double p[]);

        void get_ao_shell(const int        ishell,
                          const Basis &basis,
                                double     ao_local[],
                          const int        max_ao_geo_order,
                          const double     p[]);

        void get_ao_shell(const int        ishell,
                          const Basis &basis,
                          const int        max_ao_geo_order,
                          const double     p[]);

        bool is_same_center(const int c,
                            const std::vector<int> &carray) const;

        void transform_basis(const Basis &basis) const;

        void nullify();

        int ao_length;
};

#endif
