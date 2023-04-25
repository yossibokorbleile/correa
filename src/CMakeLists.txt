# Boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
find_package(Threads)
#set (libraries ${libraries} ${CMAKE_THREAD_LIBS_INIT} correa)
find_package(BLAS)
find_package(LAPACK)
set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

foreach(exec 2DShape Comp2DShapes)
    add_executable(${exec} "${exec}.cpp")
    target_compile_options(${exec} PRIVATE -std=c++17 -Wall -v -lapackblas_libraries)
    target_link_libraries(${exec} ${libraries} ${lapackblas_libraries} hera)
endforeach()