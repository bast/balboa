#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

#include <stdlib.h>
#include <vector>

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
//                 const double p[]) const;

        double *ao;

    private:

        Main(const Main &rhs);            // not implemented
        Main &operator=(const Main &rhs); // not implemented

        void get_ao_shell(const int        ishell,
                                double     ao_local[],
                          const int        max_ao_geo_order,
                          const double     p[]);

        void get_ao_shell(const int        ishell,
                          const int        max_ao_geo_order,
                          const double     p[]);

        bool is_same_center(const int c,
                            const std::vector<int> &carray) const;

        void transform_basis() const;

        void nullify();

        int ao_length;

        void init(const int    in_basis_type,
                  const int    in_num_centers,
                  const double in_center_coordinates[],
                  const int    in_num_shells,
                  const int    in_shell_centers[],
                  const int    in_shell_l_quantum_numbers[],
                  const int    in_shell_num_primitives[],
                  const double in_primitive_exponents[],
                  const double in_contraction_coefficients[]);
        int  get_num_centers() const;
        int  get_num_ao_slices() const;
        int  get_num_ao() const;
        int  get_num_ao_cartesian() const;
        int  get_ao_center(const int i) const;
        int  get_geo_off(const int i,
                         const int j,
                         const int k) const;
        void set_geo_off(const int geo_diff_order);

        int     num_centers; // FIXME
        int     num_shells; // FIXME
        int    *shell_l_quantum_numbers; // FIXME
        int    *shell_num_primitives; // FIXME
        double *primitive_exponents; // FIXME
        double *center_coordinates; // FIXME
        int    *shell_centers; // FIXME

        void    deallocate();

        double *shell_centers_coordinates;
        double *shell_extent_squared;
        int    *cartesian_deg;
        int    *shell_off;
        int    *spherical_deg;
        bool    is_spherical;
        int     num_ao;
        int     num_ao_cartesian;
        int     num_ao_spherical;
        int     num_ao_slices;
        int    *ao_center;
        int     geo_diff_order;
        int    *geo_off;
        double *contraction_coefficients;
        int     is_initialized;
};

#endif // MAIN_H_INCLUDED
