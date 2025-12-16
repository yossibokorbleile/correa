/*!
* @file	2DShape.cpp
* @brief extract infromation about a polygon
* @details Given a polygon stored in a file, calculate the following information about the polygon: length, area, sphericity, maximum volume inscribed ellipse, minimum volume inscribing ellipse, least square ellipse, Willmore energy, and the dimension 0 persistence diagram with essential pairing of the essential 0-cycle and essential 1-cycle.
* @author Patrice Koehl
* @author Yossi Bokor Bleile
* @date April 2023
* @version 1
* @copyright BSD 3-Clause License.
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

/*!
 * @brief Main function for analyzing a single 2D shape
 *
 * This program analyzes a single 2D polygon and computes various geometric and topological
 * properties. The polygon is centered and scaled to a normalized range before analysis.
 *
 * Computed properties include:
 * - Basic metrics: number of points, length, area, sphericity (4πA/L²)
 * - Ellipse fits: maximum inscribed, minimum inscribing, and least-squares ellipses with aspect ratios
 * - Willmore energy: measures the bending energy of the curve
 * - Persistence diagram: topological shape descriptor from persistent homology
 *
 * @param argc Number of command-line arguments
 * @param argv Array of command-line argument strings
 * @return 0 on success, 1 on failure
 *
 * @see parse_args() for command-line argument details
 */
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
	bool verbose = false; // default: no verbose output

        if (!parse_args(argc, argv, &INfile, &verbose)) return 1;

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

	if (verbose) {
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
	}

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
    std::cout << "     " << "=          2DShape -i INFILE [-v|--verbose]                                                    ="<<std::endl;
    std::cout << "     " << "=                                                                                              ="<<std::endl;
    std::cout << "     " << "=     where:                                                                                   ="<<std::endl;
    std::cout << "     " << "=                 -i  INFILE    --> Input file (Curve; ascii or csv file with 1 point / line)  ="<<std::endl;
    std::cout << "     " << "=                 -v, --verbose --> Enable verbose output (print detailed polygon info)        ="<<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "\n\n" <<std::endl;
}

/* ===============================================================================================
   Parse Argument from command line:

   =============================================================================================== */

/*!
 * @brief Parse command-line arguments for 2DShape
 *
 * Parses command-line arguments for the 2DShape program, extracting
 * the input file path and verbose flag. This is the simplest version
 * that analyzes a single polygon.
 *
 * @param argc Argument count from main()
 * @param argv Argument vector from main()
 * @param[out] infile Path to input polygon file (set by -i flag)
 * @param[out] verbose Enable verbose output (set by -v or --verbose flag, default: false)
 * @return true if arguments were parsed successfully, false otherwise
 */
bool parse_args(int argc, char **argv, std::string *infile, bool *verbose)
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
		for (int i = 1; i < argc; i++)
		{
			param = argv[i];

			// Handle flag-only arguments (no value)
			if (param == "-v" || param == "--verbose") {
				*verbose = true;
				continue;
			}

			// Skip if this is the last argument and it's not a flag
			if (i >= argc - 1) continue;

			if (param == "-i") {
				*infile = argv[i + 1];
				i++;
			}
		}
  	}
	return true;
}
