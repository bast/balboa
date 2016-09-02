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
