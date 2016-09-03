#include "Main.h"
#include "aoeval.h"

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
                                const int        basis_type,
                                const int        num_centers,
                                const double     center_coordinates[],
                                const int        num_shells,
                                const int        shell_centers[],
                                const int        shell_l_quantum_numbers[],
                                const int        shell_num_primitives[],
                                const double     primitive_exponents[],
                                const double     contraction_coefficients[])
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
