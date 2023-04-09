/*
 	CorreaBindings.h

 	Authors: Patrice Koehl, Department of Computer Science, University of California, Davis
				Yossi Bokor Bleile, Department of Mathematical Sciences, University of Aalborg, Aalborg
 	Date: April 2023
	Version: 1
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
//#include "PersistenceDiagram.h"
#include "Curvature.h"
#include "Component.h"


namespace correa{   

    // construct a polygon with focal point center of mass of the vertices supplied
    auto initialise_polygon(std::string path_to_vertices) {

        PolygonBuilder pbuilder;
        Polygon polygon;
        INOUT inout; 
        int ndim;
        int npoint;

        double *X;
        X = nullptr;

        inout.read(path_to_vertices, &ndim, &npoint, &X);
        pbuilder.clean_points(&npoint, X);
        pbuilder.buildPolygon(npoint, X, polygon);

        // Center polygon
        int iscale = 0;
        double range = 100;
        polygon.centerScale(range,iscale);
        return polygon;
    } 

    auto initialise_polygon(std::string path_to_vertices, std::string path_to_focal_point) {

        PolygonBuilder pbuilder;
        Polygon polygon;
        INOUT inout; 
        int ndim;
        int npoint;

        double *X;
        X = nullptr;

        inout.read(path_to_vertices , &ndim, &npoint, &X);
        pbuilder.clean_points(&npoint, X);
        pbuilder.buildPolygon(npoint, X, polygon);

        // Center polygon
        int iscale = 0;
        double range = 100;
        polygon.centerScale(range,iscale);
        return polygon;
    }


     auto load_polygon(std::string file_path) {
        Polygon poly = initialise_polygon(file_path);
        return poly;
    }

    /*auto wasserstein_distance(std::vector<std::pair<double,double>> diagram_1, std::vector<std::pair<double,double>> diagram_2, double wasserstein_power, double delta) {
        hera::AuctionParams<double> params;
        params.max_num_phases = 800;
        params.wasserstein_power = wasserstein_power;
        params.delta = delta;
        params.internal_p = wasserstein_power;

        auto res = hera::wasserstein_cost_detailed(diagram_1, diagram_2, params);

        std::cout << "Relative error: " << res.final_relative_error << std::endl;

        return res.distance;
    }*/

    class PyPolygon{
        using PersistenceDiagram = std::vector<std::pair<double,double>>;
        private:
            Polygon polygon;
            std::tuple<double, double, double> ellipse_min_;
            std::tuple<double, double, double> ellipse_max_;
            std::tuple<double, double, double> ellipse_lsq_;
            double willmore_;
            PersistenceDiagram  persistence_diagram_;


        public:

            PyPolygon(std::string file_path) {
                polygon = load_polygon(file_path);
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

            auto persistence_diagram() {
                return persistence_diagram_;
            }
            int size() {
                return polygon.size();
            };

            double length() {
                return polygon.length();
            };

            double area() {
                return polygon.area();
            };

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
            };

            auto ellipse_max() {
                //std::cerr << "ellipse_min_ is: (" << ellipse_min_[0] << ", " << ellipse_min_[1] << ", " << ellipse_min_[2] << ")." << std::endl;
                return ellipse_max_;
            };

            auto ellipse_max_a() {
                //std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
                return std::get<0>(ellipse_max_);
            };
            auto ellipse_max_b() {
                //std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
                return std::get<1>(ellipse_max_);
            };
            auto ellipse_max_ratio() {
                //std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
                return std::get<2>(ellipse_max_);
            };
           
           auto ellipse_min() {
                //std::cerr << "ellipse_min_ is: (" << ellipse_min_[0] << ", " << ellipse_min_[1] << ", " << ellipse_min_[2] << ")." << std::endl;
                return ellipse_min_;
            };

            auto ellipse_min_a() {
                //std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
                return std::get<0>(ellipse_min_);
            };
            auto ellipse_min_b() {
                //std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
                return std::get<1>(ellipse_min_);
            };
            auto ellipse_min_ratio() {
                //std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
                return std::get<2>(ellipse_min_);
            };

            auto ellipse_lsq() {
                //std::cerr << "ellipse_min_ is: (" << ellipse_min_[0] << ", " << ellipse_min_[1] << ", " << ellipse_min_[2] << ")." << std::endl;
                return ellipse_lsq_;
            };

            auto ellipse_lsq_a() {
                //std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
                return std::get<0>(ellipse_lsq_);
            };
            auto ellipse_lsq_b() {
                //std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
                return std::get<1>(ellipse_lsq_);
            };
            auto ellipse_lsq_ratio() {
                //std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
                return std::get<2>(ellipse_lsq_);
            };

            auto willmore() {
                //std::cerr << "willmore enegery is: " << willmore_ << "." << std::endl;
                return willmore_;
            };

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
					out << "( " << P.persistence_diagram()[i].first << ", " << P.persistence_diagram()[i].second << ")" << std::endl;
				}
                return out;            
            };
    };  

 

    class ComparePolygons {
        private:
            Polygon poly1, poly2;

        public:

            ComparePolygons(PyPolygon p1, PyPolygon p2);
            
    };


    void print_polygon(PyPolygon P) {
        std::cout << P << std::endl;
    }

   
}


#endif