/* ===============================================================================================
   Author:  Patrice Koehl & Yossi Bokor
   Date:    March 03 2023
   Version: 0.1
   =============================================================================================== */

#ifndef _CORREA_BINDINGS_H_
#define _CORREA_BINDINGS_H_
#include <iostream>
#include <include/correa_bindings.h>
#include <Polygon.h>
#include <InOut.h>
#include <Ellipse.h>
#include <Curvature.h>
#include <PH0.h>
//#include <hera/wasserstein.h>

namespace correa{

    auto initialise_polygon(std::string path_to_file) {

        PolygonBuilder pbuilder;
        Polygon polygon;
        INOUT inout; 
        int ndim;
        int npoint;

        double *X;
        X = nullptr;

        inout.read(path_to_file, &ndim, &npoint, &X);
        pbuilder.clean_points(&npoint, X);
        pbuilder.buildPolygon(npoint, X, polygon);

        // Center polygon
        int iscale = 0;
        double range = 100;
        polygon.centerScale(range,iscale);
        return polygon;
    } 


     auto load_poly(std::string file_path) {
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
        private:
            Polygon polygon;
            std::vector<double> ellipse_min_;
            std::vector<double> ellipse_max_;
            std::vector<double> ellipse_lsq_;
            double willmore_;
            std::vector<std::vector<double>> persistence_diagram;


        public:

            PyPolygon(std::string file_path) {
                polygon = load_poly(file_path);
                Ellipse ellipse;
                Curvature curv;
                double a, b; 
                ellipse.EllipseMin(polygon, &a, &b);
                ellipse_min_.push_back(a);
                ellipse_min_.push_back(b);
                ellipse_min_.push_back(a/b);
                ellipse.EllipseMax(polygon, &a, &b);
                ellipse_max_.push_back(a);
                ellipse_max_.push_back(b);
                ellipse_max_.push_back(a/b);
                ellipse.EllipseLSQ(polygon, &a, &b);
                ellipse_lsq_.push_back(a);
                ellipse_lsq_.push_back(b);
                ellipse_lsq_.push_back(a/b);

                willmore_ = curv.Willmore(polygon);
                PH0 f(polygon.vertices);
                f.Persistence();
                persistence_diagram = f.pd.points;
                std::cerr << "there are " << persistence_diagram.size() << " points in the persistence diagram." << std::endl;
            }

            //std::vector<std::vector<double>> vertices();

            auto vertices() {
                std::cerr << "there are " << polygon.vertices.size() << " vertices in this polygon" << std::endl;
                std::vector<std::vector<double>> vertices;
                for (int i = 0; i < polygon.vertices.size(); i++) {
                    std::vector<double> vert_i;
                    vert_i.push_back(polygon.vertices[i].position.x);
                    vert_i.push_back(polygon.vertices[i].position.y);
                    vertices.push_back(vert_i);
                }
                return vertices;
            };

            auto ellipse_min() {
                std::cerr << "ellipse_min_ is: (" << ellipse_min_[0] << ", " << ellipse_min_[1] << ", " << ellipse_min_[2] << ")." << std::endl;
                return ellipse_min_;
            };

            auto ellipse_max() {
                std::cerr << "ellipse_max_ is: (" << ellipse_max_[0] << ", " << ellipse_max_[1] << ", " << ellipse_max_[2] << ")." << std::endl;
                return ellipse_max_;
            };

            auto ellipse_lsq() {
                std::cerr << "ellipse_lsq_ is: (" << ellipse_lsq_[0] << ", " << ellipse_lsq_[1] << ", " << ellipse_lsq_[2] << ")." << std::endl;
                return ellipse_lsq_;
            };

            auto willmore() {
                std::cerr << "willmore enegery is: " << willmore_ << "." << std::endl;
                return willmore_;
            };

    };

    class ComparePolygons {
        private:
            Polygon poly1, poly2;

        public:
            
    };




   
}
#endif