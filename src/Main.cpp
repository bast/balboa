#include "Main.h"
#include "aoeval.h"
//#include "AOBatch.h"

#define AS_TYPE(Type, Obj) reinterpret_cast<Type *>(Obj)
#define AS_CTYPE(Type, Obj) reinterpret_cast<const Type *>(Obj)


aoeval_context_t *aoeval_new()
{
    return AS_TYPE(aoeval_context_t, new Main());
}
Main::Main()
{
//  nullify();
}


void aoeval_free(aoeval_context_t *context)
{
    if (!context) return;
    delete AS_TYPE(Main, context);
}
Main::~Main()
{
//  MemAllocator::deallocate(xyzw);
//  nullify();
}


AOEVAL_API int aoeval_set_basis(aoeval_context_t *context,
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
int Main::set_basis(const int    basis_type,
                    const int    num_centers,
                    const double center_coordinates[],
                    const int    num_shells,
                    const int    shell_centers[],
                    const int    shell_l_quantum_numbers[],
                    const int    shell_num_primitives[],
                    const double primitive_exponents[],
                    const double contraction_coefficients[])
{
    basis.init(basis_type,
               num_centers,
               center_coordinates,
               num_shells,
               shell_centers,
               shell_l_quantum_numbers,
               shell_num_primitives,
               primitive_exponents,
               contraction_coefficients);

    return 0;
}


//AOEVAL_API double *aoeval_get_ao(const aoeval_context_t *context,
AOEVAL_API double *aoeval_get_ao(      aoeval_context_t *context,
                             const bool   use_gradient,
                             const int    max_ao_geo_order,
                             const int    block_length,
                             const double p[])
{
  //return AS_CTYPE(Main, context)->get_ao(use_gradient,
    return AS_TYPE(Main, context)->get_ao(use_gradient,
                                           max_ao_geo_order,
                                           block_length,
                                           p);
}
double *Main::get_ao(const bool   use_gradient,
                 const int    max_ao_geo_order,
                 const int    block_length,
                 const double p[])
              // const double p[]) const
{
    aobatch.get_ao(basis,
                   use_gradient,
                   max_ao_geo_order,
                   block_length,
                   p);

    return aobatch.ao;
}
