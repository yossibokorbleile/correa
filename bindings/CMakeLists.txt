nanobind_add_module(_correa ../src/correa.cpp)
target_link_libraries(_correa PRIVATE ${libraries} ${lapackblas_libraries} hera)
