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

    class PyPolygon{
        private:
            Polygon polygon;
            std::vector<double> min_ellipse;
            std::vector<double> max_ellipse;
            std::vector<double> lsq_ellipse;
            double willmore_;
            std::vector<std::vector<double>> persistence_diagram;


        public:

            PyPolygon(std::string file_path) {
                polygon = load_poly(file_path);
                Ellipse ellipse;
                Curvature curv;
                double a, b; 
                ellipse.EllipseMin(polygon, &a, &b);
                min_ellipse.push_back(a);
                min_ellipse.push_back(b);
                min_ellipse.push_back(a/b);
                ellipse.EllipseMax(polygon, &a, &b);
                max_ellipse.push_back(a);
                max_ellipse.push_back(b);
                max_ellipse.push_back(a/b);
                ellipse.EllipseLSQ(polygon, &a, &b);
                lsq_ellipse.push_back(a);
                lsq_ellipse.push_back(b);
                lsq_ellipse.push_back(a/b);

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
                std::cerr << "min_ellipse is: (" << min_ellipse[0] << ", " << min_ellipse[1] << ", " << min_ellipse[2] << ")." << std::endl;
                return min_ellipse;
            };

            auto ellipse_max() {
                std::cerr << "max_ellipse is: (" << max_ellipse[0] << ", " << max_ellipse[1] << ", " << max_ellipse[2] << ")." << std::endl;
                return max_ellipse;
            };

            auto ellipse_lsq() {
                std::cerr << "lsq_ellipse is: (" << lsq_ellipse[0] << ", " << lsq_ellipse[1] << ", " << lsq_ellipse[2] << ")." << std::endl;
                return lsq_ellipse;
            };

            auto willmore() {
                std::cerr << "willmore enegery is: " << willmore_ << "." << std::endl;
                return willmore_;
            };

    };


   
}
#endif