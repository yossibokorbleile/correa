/*
 	Correa.cpp

 	Authors: Patrice Koehl, Department of Computer Science, University of California, Davis
				Yossi Bokor Bleile, Department of Mathematical Sciences, University of Aalborg, Aalborg
 	Date: April 2023
	Version: 1
*/

/*!
* @file _correa.cpp
* @brief create bindings for Correa.
*/

#include <iostream>
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/array.h>
#include "correa_bindings.h"

namespace nb = nanobind;


NB_MODULE(_correa, m) {
	/*!
	* Bind PyPolygon
	*/
    nb::class_<correa::PyPolygon>(m, "PyPolygon")
        .def(nb::init<const std::string &>())
        //.def("extractVertices", &correa::PyPolygon::extractVertices)
        .def("vertices", &correa::PyPolygon::vertices)
        .def("ellipse_min", &correa::PyPolygon::ellipse_min)
        .def("ellipse_max", &correa::PyPolygon::ellipse_max)
        .def("ellipse_lsq", &correa::PyPolygon::ellipse_lsq)
        .def("willmore", &correa::PyPolygon::willmore)
        .def("size", &correa::PyPolygon::size)
        .def("persistence_diagram", &correa::PyPolygon::persistence_diagram);

    m.def("compare_polygons", &correa::compare_polygons);
    m.def("print_polygon", &correa::print_polygon);

}

