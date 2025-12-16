/*!
* @file _correa.cpp
* @brief expose bindings for Correa.
* @author Patrice Koehl
* @author Yossi Bokor Bleile
* @date December 2025
* @version 1.1
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
        .def(nb::init<const std::string , const bool&, const bool&, const double&>())
		.def(nb::init<const std::string &, const std::string&, const bool&, const bool&, const double&>())
		.def(nb::init<const std::string &, const std::vector<double>, const bool&, const bool&, const double&>())
		.def(nb::init<const bool&, const std::string &,  const std::vector<double>, const bool&, const double&>())
		.def(nb::init<const bool&, const std::string &,  const std::vector<double>, const bool&, const double&>())
        .def("vertices", &correa::PyPolygon::vertices)
        .def("ellipse_min_a", &correa::PyPolygon::ellipse_min_a)
		.def("ellipse_min_b", &correa::PyPolygon::ellipse_min_b)
		.def("ellipse_min_ratio", &correa::PyPolygon::ellipse_min_ratio)
        .def("ellipse_max", &correa::PyPolygon::ellipse_max)
        .def("ellipse_max_a", &correa::PyPolygon::ellipse_max_a)
		.def("ellipse_max_b", &correa::PyPolygon::ellipse_max_b)
		.def("ellipse_max_ratio", &correa::PyPolygon::ellipse_max_ratio)
        .def("ellipse_lsq", &correa::PyPolygon::ellipse_lsq)
		.def("ellipse_lsq_a", &correa::PyPolygon::ellipse_lsq_a)
		.def("ellipse_lsq_b", &correa::PyPolygon::ellipse_lsq_b)
		.def("ellipse_lsq_ratio", &correa::PyPolygon::ellipse_lsq_ratio)
        .def("willmore", &correa::PyPolygon::willmore)
        .def("size", &correa::PyPolygon::size)
		.def("area", &correa::PyPolygon::area)
        .def("persistence_diagram", &correa::PyPolygon::persistence_diagram)
        .def("convertPixelsToMicrometers", &correa::PyPolygon::convertPixelsToMicrometers)
		.def("originalArea", &correa::PyPolygon::originalArea);	

    m.def("compare_polygons", &correa::compare_polygons);
	m.def("curv_ot_distance", &correa::curv_ot_distance);
	m.def("frechet_distance", &correa::frechet_distance);
	m.def("max_ellipse_distance", &correa::max_ellipse_distance);
	m.def("min_ellipse_distance", &correa::min_ellipse_distance);
	m.def("lsq_ellipse_distance", &correa::lsq_ellipse_distance);
    m.def("print_polygon", &correa::print_polygon);
	m.def("wasserstein_distance", &correa::wasserstein_distance);
	m.def("willmore_distance", &correa::willmore_distance);
	m.def("hera_wasserstein_distance", &correa::hera_wasserstein_distance);

	// Verbosity controls
	m.def("set_verbose", &correa::set_verbose, "Enable or disable verbose logging from C++ bindings");
	m.def("get_verbose", &correa::get_verbose, "Return current verbosity flag of C++ bindings");
}

