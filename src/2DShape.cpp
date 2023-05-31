/*
 	2DShape.cpp

 	Authors: Patrice Koehl, Department of Computer Science, University of California, Davis
				Yossi Bokor Bleile, Department of Mathematical Sciences, University of Aalborg, Aalborg
 	Date: April 2023
	Version: 1
*/


#include "2DShape.h"
#include "InOut.h"
#include "PolygonBuilder.h"
#include "Polygon.h"
#include "Ellipse.h"
//#include "OT1.h"

//OT1 ot1;

#include "PH0.h"
#include "Curvature.h"
#include "Component.h"


/* ===============================================================================================
   Main program
   =============================================================================================== */

int main(int argc, char **argv)
{
	correa::INOUT inout;
	correa::PolygonBuilder pbuilder;
	correa::Ellipse ellipse;
	correa::Curvature curv;

/*	==========================================================================================
	Show usage if needed
	========================================================================================== */

	if( argc < 2 )
	{
		usage(argv);
		return -1;
	}

	std::string input = argv[1];
	if( input == "-h" || input == "-help" )
	{
		usage(argv);
		return -1;
	}

/*	==========================================================================================
	Read in all inputs (some values may have been preset)
	========================================================================================== */

	std::string INfile;

        if (!parse_args(argc, argv, &INfile)) return 1;

/*	==========================================================================================
	Read in the polygon from input file
	========================================================================================== */

	int ndim;
	int npoint;

	double *X;
	X = nullptr;
	inout.read(INfile, &ndim, &npoint, &X);

/*	==========================================================================================
	Store as a polygon
	========================================================================================== */


	correa::Polygon polygon;
	pbuilder.clean_points(&npoint, X);
	pbuilder.buildPolygon(npoint, X, polygon);

	// Center polygon
	int iscale = 0;
	double range = 100;
	polygon.centerScale(range,iscale);

/*	==========================================================================================
	Compute geometric properties
	========================================================================================== */

	double a_M, b_M, r_M;
	double a_m, b_m, r_m;
	double a_l, b_l, r_l;
	double willmore;

	ellipse.EllipseMin(polygon, &a_m, &b_m);
	r_m = a_m/b_m;
	ellipse.EllipseMax(polygon, &a_M, &b_M);
	r_M = a_M/b_M;
	ellipse.EllipseLSQ(polygon, &a_l, &b_l);
	r_l = a_l/b_l;
	willmore = curv.Willmore(polygon);
	correa::PH0 f(polygon.vertices);
	f.Persistence();
	
/*	==========================================================================================
	Print results
	========================================================================================== */

	// Info on polygon :

	double l = polygon.length();
	double A = polygon.area();

	std::cout << " " << std::endl;
	std::cout << "Polygon  : " << std::endl;
	std::cout << "===========" << std::endl;
	std::cout << "Number of points in polygon       			: " << npoint << std::endl;
	std::cout << "Length of polygon                 			: " << l << std::endl;
	std::cout << "Area of polygon                   			: " << A << std::endl;
	std::cout << "Sphericity (4*Pi*Area/L^2)        			: " << 4*M_PI*A/(l*l) << std::endl;
	std::cout << "Maximum volume inscribed ellipse  			: a : " << std::setw(7) << std::fixed << std::setprecision(3) << a_M << " b : " << b_M << " Aspect ratio: " << r_M << std::endl;
	std::cout << "Minimum volume inscribing ellipse 			: a : " << std::setw(7) << std::fixed << std::setprecision(3) << a_m << " b : " << b_m << " Aspect ratio: " << r_m << std::endl;
	std::cout << "Least square ellipse              			: a : " << std::setw(7) << std::fixed << std::setprecision(3) << a_l << " b : " << b_l << " Aspect ratio: " << r_l << std::endl;
	std::cout << "Willmore energy of polygon        			: " << willmore << std::endl;
	std::cout << "Number of points in the persistence diagram 	: " << f.persistence_diagram().size() ;
	f.printPD();
	std::cout << " " << std::endl;

	return 0;

}

/* ===============================================================================================
   Usage
   =============================================================================================== */

static void usage(char** argv)
{
    std::cout << "\n\n" <<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "     " << "=                                                                                              ="<<std::endl;
    std::cout << "     " << "=                                       2DShape                                                ="<<std::endl;
    std::cout << "     " << "=                                                                                              ="<<std::endl;
    std::cout << "     " << "=     Usage is:                                                                                ="<<std::endl;
    std::cout << "     " << "=          2DShape -i INFILE                                                                   ="<<std::endl;
    std::cout << "     " << "=                                                                                              ="<<std::endl;
    std::cout << "     " << "=     where:                                                                                   ="<<std::endl;
    std::cout << "     " << "=                 -i  INFILE    --> Input file (Curve; ascii or csv file with 1 point / line)  ="<<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "\n\n" <<std::endl;
}

/* ===============================================================================================
   Parse Argument from command line:

   =============================================================================================== */

bool parse_args(int argc, char **argv, std::string *infile)
{
//
// Make sure we have at least two parameters....
//
	std::string param;
	if (argc == 1)
	{
		return false;
	}
	else
	{
		for (int i = 1; i < argc - 1; i = i + 2)
		{
			param = argv[i];

			if (param == "-i") {
				*infile = argv[i + 1];
			}
		}
  	}
	return true;
}
