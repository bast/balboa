#ifndef AOBatch_h_
#define AOBatch_h_

#include <stdlib.h>
#include <vector>

#include "Basis.h"

class AOBatch
{
    public:

        AOBatch();
        ~AOBatch();

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

    private:

        AOBatch(const AOBatch &rhs);            // not implemented
        AOBatch &operator=(const AOBatch &rhs); // not implemented

        bool is_same_center(const int c,
                            const std::vector<int> &carray) const;

        void transform_basis(const Basis &basis) const;

        void nullify();

        int     ao_length;
        double *ao;

        double *ao_compressed;
        int     ao_compressed_num;
        int    *ao_compressed_index;

        double *k_ao_compressed;
        int     k_ao_compressed_num;
        int    *k_ao_compressed_index;

        double *l_ao_compressed;
        int     l_ao_compressed_num;
        int    *l_ao_compressed_index;
};

#endif // AOBatch_h_
