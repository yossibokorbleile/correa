/*
 	Correa.cpp

 	Authors: Patrice Koehl, Department of Computer Science, University of California, Davis
				Yossi Bokor Bleile, Department of Mathematical Sciences, University of Aalborg, Aalborg
 	Date: April 2023
	Version: 1
*/

#include <iostream>
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/array.h>
#include "correa_bindings.h"

namespace nb = nanobind;


NB_MODULE(correa, m) {
    nb::class_<correa::PyPolygon>(m, "PyPolygon")
        .def(nb::init<const std::string &>())
        //.def("extractVertices", &correa::PyPolygon::extractVertices)
        .def("vertices", &correa::PyPolygon::vertices)
        .def("ellipse_min", &correa::PyPolygon::ellipse_min)
        .def("ellipse_max", &correa::PyPolygon::ellipse_max)
        .def("ellipse_lsq", &correa::PyPolygon::ellipse_lsq)
        .def("willmore", &correa::PyPolygon::willmore)
        .def("size", &correa::PyPolygon::size);

    
    m.def("load_polygon", &correa::load_polygon);
    m.def("print_polygon", &correa::print_polygon);

}

