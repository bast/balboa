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

if(ENABLE_FC_SUPPORT)
    add_library(
        numgrid_fortran
        ${CMAKE_CURRENT_LIST_DIR}/../api/numgrid.F90
        )

    target_link_libraries(
        numgrid_fortran
        numgrid
        )
endif()

install(TARGETS aoeval LIBRARY DESTINATION lib)
