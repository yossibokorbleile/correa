/* ===============================================================================================
   Author:  Patrice Koehl & Yossi Bokor
   Date:    March 03 2023
   Version: 0.1
   =============================================================================================== */


#include <nanobind/nanobind.h>
#include <include/correa_bindings.cpp>

NB_MODULE(correa, m) {
    m.def("add", &add);
}