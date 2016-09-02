add_library(
    aoeval
    SHARED
    MemAllocator.cpp
    MemAllocator.h
    AOBatch.cpp
    AOBatch.h
    Basis.cpp
    Basis.h
    CMakeLists.txt
    ao_vector.cpp
    ao_vector.h
    )

install(TARGETS aoeval LIBRARY DESTINATION lib)
