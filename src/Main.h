#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include "Basis.h"
#include "AOBatch.h"

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

//      int get_num_points() const;
//      double *get_grid() const;

    private:

        Main(const Main &rhs);            // not implemented
        Main &operator=(const Main &rhs); // not implemented

        Basis basis;
        AOBatch aobatch;
//      void nullify();

//      int get_closest_num_angular(int n) const;
//      int get_angular_order(int n) const;

//      int num_points;
//      double *xyzw;
};

#endif
