/*!
* @file correa_bindgins.h
* @brief create bindings for Correa.
* @details Generate the Python bindings for Correa using nanobind. In particular the functions to comppare two polygons.
* @author Patrice Koehl
* @author Yossi Bokor Bleile
* @date April 2023
* @version 1
* @copyright BSD 3-Clause License.
*/



#ifndef _CORREABINDINGS_H_
#define _CORREABINDINGS_H_

/* ===============================================================================================
   System includes
   =============================================================================================== */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <ctime>
#include <unistd.h>
#include <cstdlib>
#include <limits>
#include <assert.h>
#include <algorithm>
#include <nanobind/stl/pair.h>
/*================================================================================================
 Definitions for multi-threading
================================================================================================== */

#define NUM_THREADS 32

int threadids[NUM_THREADS];
pthread_t threads[NUM_THREADS];

/* ===============================================================================================
   Local includes
   =============================================================================================== */

#include "InOut.h"
#include "PolygonBuilder.h"
#include "Polygon.h"
#include "Ellipse.h"
#include "OT1.h"
#include "PH0.h"
#include "Curvature.h"
#include "Component.h"
#include "hera/wasserstein.h"
#include "Frechet.h"
#include "Vector2D.h"

namespace correa{   

	/*!
	* Caculate the Wasserstein distance between two persistence diagrams using Hera.
	* @param pd1 first persistence diagram
	* @param pd2 second persistence diagram 
	* @param q	the qth power to use, default is 2
	* @return 	qth wasserstein distance between the two persistence diagrams
	*/
	auto hera_wasserstein_distance(const std::vector<std::pair<double,double>> pd1, const std::vector<std::pair<double,double>> pd2, int q = 2) {
		hera::AuctionParams<double> hera_params;
		hera_params.max_num_phases = 800;
		hera_params.wasserstein_power = q;
		hera_params.internal_p = q;
		const hera::AuctionParams<double> params = hera_params;
		double w_q_dist = hera::wasserstein_dist<std::vector<std::pair<double,double>>>(pd1, pd2, params);
		return w_q_dist;
	};
	// construct a polygon with focal point center of mass of the vertices supplied
	auto initialise_polygon(std::string path_to_vertices, bool clean_points, bool scale_by_area, double convert_to_microns_factor) {
		PolygonBuilder pbuilder;
		Polygon polygon;
		INOUT inout; 
		int ndim;
		int npoint;

		double *X;
		X = nullptr;

		inout.read(path_to_vertices, &ndim, &npoint, &X);
		std::cout << "Number of points before cleaning: " << npoint << std::endl;
		if (clean_points) {
			pbuilder.clean_points(&npoint, X);
		} else {
			std::cout << "No cleaning of points" << std::endl;
		
		}
		pbuilder.buildPolygon(npoint, X, polygon);
		// polygon.labelPolygon(polygon);
		polygon.originalArea_ = polygon.area();

		if (convert_to_microns_factor != 1.0) {
			pbuilder.convertPixelsToMicrometers(polygon, convert_to_microns_factor);
		}
		// Center polygon
		if (scale_by_area) {
			polygon.centerScaleArea();
		} else {
			double iscale = 0;
		double range = 100;
			polygon.centerScale(range,iscale);
		}
		return polygon;
	} 

	auto initialise_polygon(std::string path_to_vertices, std::string path_to_focal_point, bool clean_points, bool scale_by_area, double convert_to_microns_factor) {
		PolygonBuilder pbuilder;
		Polygon polygon;
		INOUT inout; 
		int ndim;
		int npoint;

		double *X;
		X = nullptr;

		inout.read(path_to_vertices , &ndim, &npoint, &X);
		std::string line;
   		double val;
		std::vector<double> focal_;
   		std::ifstream in(path_to_focal_point);

   		while (getline(in, line)){
				size_t start_pos = 0;
				while ((start_pos = line.find(",", start_pos)) != std::string::npos)
				{
					line.replace(start_pos, 1, " ");
					start_pos += 1; // Handles case where 'to' is a substring of 'from'
				}
				istringstream iss(line);
				while (iss >> val){
					focal_.push_back(val);
			};
		};
		Vector2D focal(focal_[0], focal_[1]);
		std::cout << "Number of points before cleaning: " << npoint << std::endl;
		if (clean_points) {
			pbuilder.clean_points(&npoint, X);
		} else {
			std::cout << "No cleaning of points" << std::endl;
		
		}
		pbuilder.buildPolygon(npoint, X, polygon);
		
		// Shift polygon to focal point
		polygon.shift(focal);
		polygon.originalArea_ = polygon.area();
		if (convert_to_microns_factor != 1.0) {
			pbuilder.convertPixelsToMicrometers(polygon, convert_to_microns_factor);
		}
		if (scale_by_area) {
			polygon.centerScaleArea();
		} else {
			double iscale = 0;
			double range = 100;
			polygon.centerScale(range,iscale);
		}
		return polygon;
	}

	auto initialise_polygon(std::string path_to_vertices, std::vector<double> focal_point, bool clean_points, bool scale_by_area, double convert_to_microns_factor) {
		PolygonBuilder pbuilder;
		Polygon polygon;
		INOUT inout; 
		int ndim;
		int npoint;

		double *X;
		X = nullptr;

		inout.read(path_to_vertices , &ndim, &npoint, &X);
		Vector2D focal (focal_point[0], focal_point[1]); 
		std::cout << "focal is (" << focal.x << ", " << focal.y << ")." << std::endl;
		std::cout << "Number of points before cleaning: " << npoint << std::endl;
		if (clean_points) {
			pbuilder.clean_points(&npoint, X);
		} else {
			std::cout << "No cleaning of points, so we still have " << npoint << " points" << std::endl;
		}
		
		pbuilder.buildPolygon(npoint, X, polygon);
		std::cout << "build the polygon with " << npoint << " points" << std::endl;
		// Shift polygon to focal point
		polygon.shift(focal);
		polygon.originalArea_ = polygon.area();
		if (scale_by_area) {
			polygon.centerScaleArea();
		} else {
			double iscale = 0;
			double range = 100;
			polygon.centerScale(range,iscale);
		}
		return polygon;
	}


	 auto load_polygon(std::string file_path, bool clean_points, bool scale_by_area, double convert_to_microns_factor = 1.0) {
		Polygon polygon = initialise_polygon(file_path, clean_points, scale_by_area, convert_to_microns_factor);
		return polygon;
	}

	auto load_polygon(std::string path_to_vertices, std::string path_to_focal, bool clean_points, bool scale_by_area, double convert_to_microns_factor) {
		Polygon polygon = initialise_polygon(path_to_vertices, path_to_focal, clean_points, scale_by_area, convert_to_microns_factor);
		return polygon;
	}

	auto load_polygon(std::string path_to_vertices, std::vector<double> focal, bool clean_points, bool scale_by_area, double convert_to_microns_factor = 1.0) {
		Polygon polygon = initialise_polygon(path_to_vertices, focal, clean_points, scale_by_area, convert_to_microns_factor);
		return polygon;
	}

	/*!
	* Expose polygons to python
	*/
	class PyPolygon{
		using PersistenceDiagram = std::vector<std::pair<double,double>>;
		private:
			std::tuple<double, double, double> ellipse_min_;
			std::tuple<double, double, double> ellipse_max_;
			std::tuple<double, double, double> ellipse_lsq_;
			double willmore_;
			PersistenceDiagram  persistence_diagram_;
			double area() {
				return polygon.area();
			}

		public:
			Polygon polygon;
			PyPolygon(bool test, std::string file_path, std::vector<double> focal_point, bool scale_by_area, double convert_to_microns_factor) {
				polygon = load_polygon(file_path, focal_point, test, scale_by_area, convert_to_microns_factor);
				std::cout << "polygon loaded with test " << test << std::endl;
				PH0 f(polygon.vertices);
				f.Persistence();
				std::cout << "persistence diagram calculated" << std::endl;
				persistence_diagram_ = f.persistence_diagram();
			}

			PyPolygon(std::string file_path, bool clean_points, bool scale_by_area, double convert_to_microns_factor) {
				polygon = load_polygon(file_path, clean_points, scale_by_area, convert_to_microns_factor);
				
				Ellipse ellipse;
				Curvature curv;
				double a, b; 
				ellipse.EllipseMin(polygon, &a, &b);
				std::get<0>(ellipse_min_)= a;
				std::get<1>(ellipse_min_) = b;
				std::get<2>(ellipse_min_)= a/b;
				ellipse.EllipseMax(polygon, &a, &b);
				std::get<0>(ellipse_max_)= a;
				std::get<1>(ellipse_max_) = b;
				std::get<2>(ellipse_max_)= a/b;
				ellipse.EllipseLSQ(polygon, &a, &b);
				std::get<0>(ellipse_lsq_)= a;
				std::get<1>(ellipse_lsq_) = b;
				std::get<2>(ellipse_lsq_)= a/b;

				willmore_ = curv.Willmore(polygon);
				PH0 f(polygon.vertices);
				f.Persistence();
				persistence_diagram_ = f.persistence_diagram();
			}

			PyPolygon(std::string file_path, std::string focal_path, bool clean_points, bool scale_by_area, double convert_to_microns_factor) {
				polygon = load_polygon(file_path, focal_path, clean_points, scale_by_area, convert_to_microns_factor);
				Ellipse ellipse;
				Curvature curv;
				double a, b; 
				ellipse.EllipseMin(polygon, &a, &b);
				std::get<0>(ellipse_min_)= a;
				std::get<1>(ellipse_min_) = b;
				std::get<2>(ellipse_min_)= a/b;
				ellipse.EllipseMax(polygon, &a, &b);
				std::get<0>(ellipse_max_)= a;
				std::get<1>(ellipse_max_) = b;
				std::get<2>(ellipse_max_)= a/b;
				ellipse.EllipseLSQ(polygon, &a, &b);
				std::get<0>(ellipse_lsq_)= a;
				std::get<1>(ellipse_lsq_) = b;
				std::get<2>(ellipse_lsq_)= a/b;

				willmore_ = curv.Willmore(polygon);
				PH0 f(polygon.vertices);
				f.Persistence();
				persistence_diagram_ = f.persistence_diagram();
			}

			PyPolygon(std::string file_path, std::vector<double> focal_point, bool clean_points, bool scale_by_area, double convert_to_microns_factor) {
				polygon = load_polygon(file_path, focal_point, clean_points, scale_by_area, convert_to_microns_factor);
				if (convert_to_microns_factor != 1.0) {
					PolygonBuilder pbuilder;
					pbuilder.convertPixelsToMicrometers(polygon, convert_to_microns_factor);
				}
				Ellipse ellipse;
				Curvature curv;
				double a, b; 
				ellipse.EllipseMin(polygon, &a, &b);
				std::get<0>(ellipse_min_)= a;
				std::get<1>(ellipse_min_) = b;
				std::get<2>(ellipse_min_)= a/b;
				ellipse.EllipseMax(polygon, &a, &b);
				std::get<0>(ellipse_max_)= a;
				std::get<1>(ellipse_max_) = b;
				std::get<2>(ellipse_max_)= a/b;
				ellipse.EllipseLSQ(polygon, &a, &b);
				std::get<0>(ellipse_lsq_)= a;
				std::get<1>(ellipse_lsq_) = b;
				std::get<2>(ellipse_lsq_)= a/b;

				willmore_ = curv.Willmore(polygon);
				PH0 f(polygon.vertices);
				f.Persistence();
				persistence_diagram_ = f.persistence_diagram();
			}

			//std::vector<std::vector<double>> vertices();

			//auto polygon() {
			//	return polygon;
			//}
			/*!
			* @return persistence diagram of the radial function from the center.
			*/
			PersistenceDiagram persistence_diagram() {
				return persistence_diagram_;
			}

			/*!
			* @return persistence diagram of the radial function from the center in a format for Python.
			*/
			//std::vector<std:: persistence_diagram() {
			//	return persistence_diagram_;
			//}
			/*!
			* @return number of vertices in the polygon
			*/
			int size() {
				return polygon.size();
			}

			/*!
			* @return length of the polygon
			*/
			double length() {
				return polygon.length();
			}

			/*!
			* @return length of the polygon
			*/
			double area() {
				return polygon.area();
			}

			/*!
			* @return std::vector<std::vector<double>> of vertices of the polygon
			*/
			auto vertices() {
				//std::cerr << "there are " << polygon.vertices.size() << " vertices in this polygon" << std::endl;
				std::vector<std::vector<double>> vertices;
				for (int i = 0; i < polygon.vertices.size(); i++) {
					std::vector<double> vert_i;
					vert_i.push_back(polygon.vertices[i].position.x);
					vert_i.push_back(polygon.vertices[i].position.y);
					vertices.push_back(vert_i);
				}
				return vertices;
			}

			/*!
			* @return parameters (major axis, minor axis, ratio) of the maximum inscribed ellipse
			*/
			auto ellipse_max() {
				//std::cerr << "ellipse_min_ is: (" << ellipse_min_[0] << ", " << ellipse_min_[1] << ", " << ellipse_min_[2] << ")." << std::endl;
				return ellipse_max_;
			}

			/*!
			* @return major axis of maximum inscribed ellipse
			*/
			auto ellipse_max_a() {
				//std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
				return std::get<0>(ellipse_max_);
			}

			/*!
			* @return minor axis of maximum inscribed ellipse
			*/
			auto ellipse_max_b() {
				//std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
				return std::get<1>(ellipse_max_);
			}

			/*!
			* @return ratio of maximum inscribed ellipse
			*/
			auto ellipse_max_ratio() {
				//std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
				return std::get<2>(ellipse_max_);
			}
		
			/*!
			* @return parameters (major axis, minor axis, ratio) of the minimum inscribing ellipse
			*/
		   auto ellipse_min() {
				//std::cerr << "ellipse_min_ is: (" << ellipse_min_[0] << ", " << ellipse_min_[1] << ", " << ellipse_min_[2] << ")." << std::endl;
				return ellipse_min_;
			}

			/*!
			* @return major axis of minimum inscribing ellipse
			*/
			auto ellipse_min_a() {
				//std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
				return std::get<0>(ellipse_min_);
			}

			/*!
			* @return minor axis of minimum inscribing ellipse
			*/
			auto ellipse_min_b() {
				//std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
				return std::get<1>(ellipse_min_);
			}

			/*!
			* @return ratio axis of minimum inscribing ellipse
			*/			
			auto ellipse_min_ratio() {
				//std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
				return std::get<2>(ellipse_min_);
			}

			/*!
			* @return parameters (major axis, minor axis, ratio) of the least squared ellipse
			*/
			auto ellipse_lsq() {
				//std::cerr << "ellipse_min_ is: (" << ellipse_min_[0] << ", " << ellipse_min_[1] << ", " << ellipse_min_[2] << ")." << std::endl;
				return ellipse_lsq_;
			}

			/*!
			* @return major axis of least square ellipse
			*/
			auto ellipse_lsq_a() {
				//std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
				return std::get<0>(ellipse_lsq_);
			};

			/*!
			* @return minor axis of least square ellipse
			*/
			auto ellipse_lsq_b() {
				//std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
				return std::get<1>(ellipse_lsq_);
			};

			/*!
			* @return ratio axis of least square ellipse
			*/
			auto ellipse_lsq_ratio() {
				//std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
				return std::get<2>(ellipse_lsq_);
			};

			/*!
			* @return Willmore energy of the polygon
			*/
			auto willmore() {
				//std::cerr << "willmore enegery is: " << willmore_ << "." << std::endl;
				return willmore_;
			};

			/*!
			* Convert polygon coordinates from pixels to micrometers using the given scale factor
			* @param scale the scale factor to convert pixels to micrometers
			*/
			void convertPixelsToMicrometers(double scale) {
				PolygonBuilder pbuilder;
				pbuilder.convertPixelsToMicrometers(polygon, scale);
				
				// Recompute persistence diagram after coordinate conversion
				PH0 f(polygon.vertices);
				f.Persistence();
				persistence_diagram_ = f.persistence_diagram();
			};

			/*!
			* @return original area of the polygon
			*/
			auto originalArea() {
				return polygon.originalArea();
			}

			/*!
			* print information about the polygon
			*/
			friend ostream &operator<<( ostream &out, PyPolygon &P ) { 
				out <<  "====================================\n";
				out << "Displaying Polygon Information:     \n";
				out << "Number of vertices:                 " << P.size() << "\n";
				out << "Length:                             " << P.length() << "\n";
				out << "Area:                               " << P.area() << "\n";
				out << "Sphericity (4*Pi*Area/L^2):         " <<  4*M_PI*P.area()/(P.length()*P.length()) << "\n";
				out << "Maximum volume inscribed ellipse:   a: " << P.ellipse_max_a() << " b: " << P.ellipse_max_b() << " ratio: " << P.ellipse_max_ratio()<< "\n";
				out << "Minimum volume inscribing ellipse:  a: " << P.ellipse_min_a() << " b: " << P.ellipse_min_b() << " ratio: " << P.ellipse_min_ratio() << "\n";
				out << "Least square ellipse:               a: " << P.ellipse_lsq_a() << " b: " << P.ellipse_lsq_b() << " ratio: " << P.ellipse_lsq_ratio() << "\n";
				out << "Wilmore energy:                     " << P.willmore() << "\n";
				out << "Persistence diagram: number of points:     " << P.persistence_diagram().size() << "\n";
				for (int i = 0; i < P.persistence_diagram().size(); i++) {
					out << "(" << P.persistence_diagram()[i].first << ", " << P.persistence_diagram()[i].second << ")" << std::endl;
				}
				return out;            
			};
			friend std::vector<double> compare_polygons(PyPolygon poly1, PyPolygon poly2);

			friend double wasserstein_distance(PyPolygon poly1, PyPolygon poly2, int q);

			friend double frechet_distance(PyPolygon poly1, PyPolygon poly2);

			friend double max_ellipse_distance(PyPolygon poly1, PyPolygon poly2);

			friend double min_ellipse_distance(PyPolygon poly1, PyPolygon poly2);

			friend double lsq_ellipse_distance(PyPolygon poly1, PyPolygon poly2);

			friend double willmore_distance(PyPolygon poly1, PyPolygon poly2);

			friend double curv_ot_distance(PyPolygon poly1, PyPolygon poly2);
	};  

 
	/*!
	* wrapper to compare pairs of polygons.
	* 
	* Has a variety of methods to calculate the distances between a pair of polygons. 
	*
	* Prints their comparisons
	*/
	std::vector<double> compare_polygons(PyPolygon poly1, PyPolygon poly2, int q, bool verbose) {
		Frechet frechet;
		Ellipse ellipse;
		Curvature curv;
		const std::vector<std::pair<double,double>> pd1 = poly1.persistence_diagram();
		const std::vector<std::pair<double,double>> pd2 = poly2.persistence_diagram();
		double dWasserstein = hera_wasserstein_distance(pd1, pd2, q=q);
		double dFrechet = frechet.dFD(poly1.polygon, poly2.polygon);
		double a1_M, b1_M, a2_M, b2_M;
		double dMax = ellipse.dEllipseMax(poly1.polygon, poly2.polygon, &a1_M, &b1_M, &a2_M, &b2_M);
		double a1_m, b1_m, a2_m, b2_m;
		double dMin = ellipse.dEllipseMin(poly1.polygon, poly2.polygon, &a1_m, &b1_m, &a2_m, &b2_m);
		double a1_lsq, b1_lsq, a2_lsq, b2_lsq;
		double dLSQ = ellipse.dEllipseLSQ(poly1.polygon, poly2.polygon, &a1_lsq, &b1_lsq, &a2_lsq, &b2_lsq);
		double dWillmore = std::abs(curv.Willmore(poly1.polygon) - curv.Willmore(poly2.polygon));
		double dCurvOT = curv.curvOT(poly1.polygon, poly2.polygon);
		if (verbose) {
			int digit = floor (log10 (q)) + 1;
			std::cerr << "The two polygons are:" << std::endl;
			std::cerr << "Polygon 1:" << std::endl;
			std::cerr << poly1 << std::endl;
			std::cerr << "Polygon 2:" << std::endl;
			std::cerr << poly2 << std::endl;
			std::cerr << "and the distances between them are:" << std::endl;
			std::cerr << "Wasserstein distance (q="<<q<<"):" << std::setw(10) << std::fixed  << dWasserstein << std::endl;
			std::cerr << "Frechet distance:"  << std::setw(19+digit) << std::fixed << dFrechet << std::endl;
			std::cerr << "Max Ellipse distance:"  << std::setw(14+digit) << std::fixed << dMax << std::endl;
			std::cerr << "Min Ellipse distance:"  << std::setw(14+digit) << std::fixed << dMin << std::endl;
			std::cerr << "LSQ Ellipse distance:"  << std::setw(14+digit) << std::fixed << dLSQ << std::endl;
			std::cerr << "Willmore distance:"  << std::setw(17+digit) << std::fixed << dWillmore << std::endl;
			std::cerr << "Curv OT distance:"  << std::setw(18+digit) << std::fixed << dCurvOT << std::endl;
		}
		return {dWasserstein, dFrechet, dMax, dMin, dLSQ, dWillmore, dCurvOT};
	};		

	double wasserstein_distance(PyPolygon poly1, PyPolygon poly2, int q) {
		const std::vector<std::pair<double,double>> pd1 = poly1.persistence_diagram();
		const std::vector<std::pair<double,double>> pd2 = poly2.persistence_diagram();
		return  hera_wasserstein_distance(pd1, pd2, q=q);
	};

	double frechet_distance(PyPolygon poly1, PyPolygon poly2) {
		Frechet frechet;
		return frechet.dFD(poly1.polygon, poly2.polygon);
	};

	double max_ellipse_distance(PyPolygon poly1, PyPolygon poly2) {
		Ellipse ellipse;
		double a1_M, b1_M, a2_M, b2_M;
		return ellipse.dEllipseMax(poly1.polygon, poly2.polygon, &a1_M, &b1_M, &a2_M, &b2_M);
	};

	double min_ellipse_distance(PyPolygon poly1, PyPolygon poly2) {
		Ellipse ellipse;
		double a1_m, b1_m, a2_m, b2_m;
		return ellipse.dEllipseMin(poly1.polygon, poly2.polygon, &a1_m, &b1_m, &a2_m, &b2_m);
	};

	double lsq_ellipse_distance(PyPolygon poly1, PyPolygon poly2) { 
		Ellipse ellipse;
		double a1_lsq, b1_lsq, a2_lsq, b2_lsq;
		return ellipse.dEllipseLSQ(poly1.polygon, poly2.polygon, &a1_lsq, &b1_lsq, &a2_lsq, &b2_lsq);
	};

	double willmore_distance(PyPolygon poly1, PyPolygon poly2) {
		Curvature curv;
		return  std::abs(curv.Willmore(poly1.polygon) - curv.Willmore(poly2.polygon));
	};
		
	double curv_ot_distance(PyPolygon poly1, PyPolygon poly2) {
		Curvature curv;
		return curv.curvOT(poly1.polygon, poly2.polygon);
	};

	/*!
	* print the information about a PyPolygon
	*/
	void print_polygon(PyPolygon P) {
		std::cerr << P << std::endl;
	};
}


#endif