#ifndef BALBOA_H_INCLUDED
#define BALBOA_H_INCLUDED

#ifndef BALBOA_API
#  ifdef _WIN32
#    if defined(BALBOA_BUILD_SHARED) /* build dll */
#      define BALBOA_API __declspec(dllexport)
#    elif !defined(BALBOA_BUILD_STATIC) /* use dll */
#      define BALBOA_API __declspec(dllimport)
#    else /* static library */
#      define BALBOA_API
#    endif
#  else
#    if __GNUC__ >= 4
#      define BALBOA_API __attribute__((visibility("default")))
#    else
#      define BALBOA_API
#    endif
#  endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct balboa_context_s;
typedef struct balboa_context_s balboa_context_t;

BALBOA_API balboa_context_t *balboa_new();
BALBOA_API void balboa_free(balboa_context_t *context);

BALBOA_API int balboa_set_basis(balboa_context_t *context,
                                const int    basis_type,
                                const int    num_centers,
                                const double center_coordinates[],
                                const int    num_shells,
                                const int    shell_centers[],
                                const int    shell_l_quantum_numbers[],
                                const int    shell_num_primitives[],
                                const double primitive_exponents[],
                                const double contraction_coefficients[]);

BALBOA_API int balboa_get_buffer_len(const balboa_context_t *context,
                                     const int max_geo_order,
                                     const int num_points);

BALBOA_API int balboa_get_ao(const balboa_context_t *context,
                             const int    max_geo_order,
                             const int    num_points,
                             const double p[],
                                   double buf[]);

#ifdef __cplusplus
}
#endif

#endif /* BALBOA_H_INCLUDED */
