#ifndef AOEVAL_H_INCLUDED
#define AOEVAL_H_INCLUDED

#ifndef AOEVAL_API
#  ifdef _WIN32
#    if defined(AOEVAL_BUILD_SHARED) /* build dll */
#      define AOEVAL_API __declspec(dllexport)
#    elif !defined(AOEVAL_BUILD_STATIC) /* use dll */
#      define AOEVAL_API __declspec(dllimport)
#    else /* static library */
#      define AOEVAL_API
#    endif
#  else
#    if __GNUC__ >= 4
#      define AOEVAL_API __attribute__((visibility("default")))
#    else
#      define AOEVAL_API
#    endif
#  endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct aoeval_context_s;
typedef struct aoeval_context_s aoeval_context_t;

AOEVAL_API aoeval_context_t *aoeval_new();
AOEVAL_API void aoeval_free(aoeval_context_t *context);

#ifdef __cplusplus
}
#endif

#endif /* AOEVAL_H_INCLUDED */
