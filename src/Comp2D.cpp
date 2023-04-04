/* ===============================================================================================
   Author:  Patrice Koehl
   Date:    10/19/2021
   Version: 1
   =============================================================================================== */

/* ===============================================================================================
   Includes
   =============================================================================================== */

#include "Comp2D.h"


/* ===============================================================================================
   Main program
   =============================================================================================== */

int main(int argc, char **argv)
{

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

	std::string INfile1;
	std::string INfile2;
	std::string INfocal1;
	std::string INfocal2;
	int disttype = 0;

        if (!parse_args(argc, argv, &INfile1, &INfile2, &INfocal1, &INfocal2, &disttype)) return 1;

/*	==========================================================================================
	Read in the polygon1 from input file1
	========================================================================================== */

	int ndim;
	int npoint1;

	double *X1;
	X1 = nullptr;
	inout.read(INfile1, &ndim, &npoint1, &X1);

/*	==========================================================================================
	Store polygon1 as a polygon
	========================================================================================== */

	Polygon polygon1;
	pbuilder.clean_points(&npoint1, X1);
	pbuilder.buildPolygon(npoint1, X1, polygon1);

	// Center polygon
	int iscale = 0;
	double range = 100;
	polygon1.centerScale(range,iscale);

/*	==========================================================================================
	Read in the polygon2 from input file2
	========================================================================================== */

	int npoint2;

	double *X2;
	X2 = nullptr;
	inout.read(INfile2, &ndim, &npoint2, &X2);

/*	==========================================================================================
	Store polygon2 as a polygon
	========================================================================================== */

	Polygon polygon2;
	pbuilder.clean_points(&npoint2, X2);
	pbuilder.buildPolygon(npoint2, X2, polygon2);

	// Center polygon
	iscale = 0;
	range = 100;
	polygon2.centerScale(range,iscale);

/*	==========================================================================================
	Compute distances
	========================================================================================== */

	double dFrechet, dE_M, dE_m, dE_l, dW, dcOT1, pdWasserstein;
	double a1_M, b1_M, r1_M, a2_M, b2_M, r2_M;
	double a1_m, b1_m, r1_m, a2_m, b2_m, r2_m;
	double a1_l, b1_l, r1_l, a2_l, b2_l, r2_l;
	double willmore1, willmore2;

	if(disttype == 0) {
		dFrechet = frechet.dFD(polygon1, polygon2);
	} else if(disttype == 1) {
		dE_m = ellipse.dEllipseMin(polygon1, polygon2, &a1_m, &b1_m, &a2_m, &b2_m);
		r1_m = a1_m/b1_m; r2_m = a2_m/b2_m;
		dE_M = ellipse.dEllipseMax(polygon1, polygon2, &a1_M, &b1_M, &a2_M, &b2_M);
		r1_M = a1_M/b1_M; r2_M = a2_M/b2_M;
		dE_l = ellipse.dEllipseLSQ(polygon1, polygon2, &a1_l, &b1_l, &a2_l, &b2_l);
		r1_l = a1_l/b1_l; r2_l = a2_l/b2_l;
	} else if(disttype == 2) {
		willmore1 = curv.Willmore(polygon1);
		willmore2 = curv.Willmore(polygon2);
		dW = std::abs(willmore1-willmore2);
		dcOT1 = curv.curvOT(polygon1, polygon2);
	} else if(disttype == 3){
		PH0 f1(polygon1.vertices);
		PH0 f2(polygon2.vertices);
		f1.Persistence();
		f2.Persistence();
		pdWasserstein = WassersteinDistance(f1.pd, f2.pd);
	} else {
		dFrechet = frechet.dFD(polygon1, polygon2);
		dE_m = ellipse.dEllipseMin(polygon1, polygon2, &a1_m, &b1_m, &a2_m, &b2_m);
		r1_m = a1_m/b1_m; r2_m = a2_m/b2_m;
		dE_M = ellipse.dEllipseMax(polygon1, polygon2, &a1_M, &b1_M, &a2_M, &b2_M);
		r1_M = a1_M/b1_M; r2_M = a2_M/b2_M;
		dE_l = ellipse.dEllipseLSQ(polygon1, polygon2, &a1_l, &b1_l, &a2_l, &b2_l);
		r1_l = a1_l/b1_l; r2_l = a2_l/b2_l;
		willmore1 = curv.Willmore(polygon1);
		willmore2 = curv.Willmore(polygon2);
		dW = std::abs(willmore1-willmore2);
		dcOT1 = curv.curvOT(polygon1, polygon2);
		PH0 f1(polygon1.vertices);
		PH0 f2(polygon2.vertices);
		f1.Persistence();
		f2.Persistence();
		std::cout << "got to do the persistence" << std::endl;
		double pdWasserstein;
		pdWasserstein = WassersteinDistance(f1.pd, f2.pd);
	}

/*	==========================================================================================
	Print results
	========================================================================================== */

	// Info on polygon 1:

	double l1 = polygon1.length();
	double A1 = polygon1.area();

	std::cout << " " << std::endl;
	std::cout << "Polygon 1 : " << std::endl;
	std::cout << "===========" << std::endl;
	std::cout << "Number of points in polygon       : " << npoint1 << std::endl;
	std::cout << "Length of polygon                 : " << l1 << std::endl;
	std::cout << "Area of polygon                   : " << A1 << std::endl;
	std::cout << "Sphericity (4*Pi*Area/L^2)        : " << 4*M_PI*A1/(l1*l1) << std::endl;
	if(disttype==1 || disttype==4) {
		std::cout << "Maximum volume inscribed ellipse  : a : " << std::setw(7) << std::fixed << std::setprecision(3) << a1_M << " b : " << b1_M << " Aspect ratio: " << r1_M << std::endl;
		std::cout << "Minimum volume inscribing ellipse : a : " << std::setw(7) << std::fixed << std::setprecision(3) << a1_m << " b : " << b1_m << " Aspect ratio: " << r1_m << std::endl;
		std::cout << "Least square ellipse              : a : " << std::setw(7) << std::fixed << std::setprecision(3) << a1_l << " b : " << b1_l << " Aspect ratio: " << r1_l << std::endl;
	}
	if(disttype==2 || disttype==4) {
		std::cout << "Willmore energy of polygon        : " << willmore1 << std::endl;
	}
	if(disttype==3 || disttype==4) {
		std::cout << "Persistence diagram of polygon    : ";
		std::vector<std::vector<double>> pdpoints;
		int num_points;
		pdpoints =pd1.points;
		num_points = pdpoints.size();
		std::cout<< "There are " << num_points << " points in the persistence diagram." << std::endl;
		for (int i = 0; i < num_points; i++){
			std::vector<double> p = pdpoints[i];
			std::cout << "(" << p[0] << ", " << p[1] << ")" << std::endl;
		};
	}

	std::cout << " " << std::endl;

	// Info on polygon 2:

	double l2 = polygon2.length();
	double A2 = polygon2.area();

	std::cout << " " << std::endl;
	std::cout << "Polygon 2 : " << std::endl;
	std::cout << "===========" << std::endl;
	std::cout << "Number of points in polygon       : " << npoint2 << std::endl;
	std::cout << "Length of polygon                 : " << l2 << std::endl;
	std::cout << "Area of polygon                   : " << A2 << std::endl;
	std::cout << "Sphericity (4*Pi*Area/L^2)        : " << 4*M_PI*A2/(l2*l2) << std::endl;
	if(disttype==1 || disttype==4) {
		std::cout << "Maximum volume inscribed ellipse  : a : " << std::setw(7) << std::fixed << std::setprecision(3) << a2_M << " b : " << b2_M << " Aspect ratio: " << r2_M << std::endl;
		std::cout << "Minimum volume inscribing ellipse : a : " << std::setw(7) << std::fixed << std::setprecision(3) << a2_m << " b : " << b2_m << " Aspect ratio: " << r2_m << std::endl;
		std::cout << "Least square ellipse              : a : " << std::setw(7) << std::fixed << std::setprecision(3) << a2_l << " b : " << b2_l << " Aspect ratio: " << r2_l << std::endl;
	}
	if(disttype==2 || disttype==4) {
		std::cout << "Willmore energy of polygon        : " << willmore2 << std::endl;
	}
	if(disttype==3 || disttype==4) {
		std::cout << "Persistence diagram of polygon    : ";
		std::vector<std::vector<double>> pdpoints;
		pdpoints =pd2.points;
		int num_points = pdpoints.size();
		std::cout<< "There are " << num_points << " points in the persistence diagram." << std::endl;
		for (int i = 0; i < num_points; i++){
			std::vector<double> p = pdpoints[i];
			std::cout << "(" << p[0] << ", " << p[1] << ")" << std::endl;
		};
	}

	std::cout << " " << std::endl;

	std::cout << std::setw(10) << std::fixed << std::setprecision(5);

	// Distances between the two polygons
	std::cout << " " << std::endl;
	std::cout << "Distances between the 2 polygons : " << std::endl;
	std::cout << "===================================" << std::endl;
	std::cout << "Sphericity 4*Pi*(A1/l1^2 - A2/l2^2) : " << 4*M_PI*std::abs(A1/(l1*l1)-A2/(l2*l2)) << std::endl;
	if(disttype==0 || disttype == 4) {
		std::cout << "Frechet distance                                    : " << dFrechet << std::endl;
	}
	if(disttype==1 || disttype == 4) {
		std::cout << "Distance (inscribed ellipse)                        : " << dE_M << std::endl;
		std::cout << "Distance (inscribing ellipse)                       : " << dE_m << std::endl;
		std::cout << "Distance (least square ellipse)                     : " << dE_l << std::endl;
	}
	if(disttype==2 || disttype == 4) {
		std::cout << "Distance (Willmore)                                 : " << dW << std::endl;
	}
	if(disttype == 4) {
		std::cout << "Distance (Wasserstein-curvature)                    : " << std::scientific << dcOT1 << std::endl;
	}
	if(disttype==3 || disttype == 4) {
		std::cout << "Distance (2-Wasserstein between persitence diagrams): " << std::setprecision(10) << pdWasserstein << std::endl;
		cout << "The first persistence diagram is: " << endl;
	for (int i = 0; i < pd1.np; i++){
		vector<double> pt;
		pt = pd1.points[i];
		cout << "point " << i+1 << " out of " << pd1.np << " is (" << pt[0] << ", " << pt[1] << ")" << endl;
	};

	cout << "and the second persistence diagram is: " << endl;
	for (int i = 0; i < pd2.np; i++){
		vector<double> pt;
		pt = pd2.points[i];
		cout << "(" << pt[0] << ", " << pt[1] << ")" << endl;
	};
	cout << " " << endl;
	
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
    std::cout << "     " << "=                                  Comp2DShapes                                                ="<<std::endl;
    std::cout << "     " << "=                                                                                              ="<<std::endl;
    std::cout << "     " << "=     Usage is:                                                                                ="<<std::endl;
    std::cout << "     " << "=          Comp2DShapes -i1 FILE1 -i2 FILE2 -d disttype                                        ="<<std::endl;
    std::cout << "     " << "=                                                                                              ="<<std::endl;
    std::cout << "     " << "=     where:                                                                                   ="<<std::endl;
    std::cout << "     " << "=               -c1 FILE1     --> Input file (Curve; ascii or csv file with 1 point / line)    ="<<std::endl;
    std::cout << "     " << "=               -f1 FILE2     --> Input file (Point; ascii or csv file with 1 point, 1 line)   ="<<std::endl;
	std::cout << "     " << "=               -c2 FILE3     --> Input file (Curve; ascii or csv file with 1 point / line)    ="<<std::endl;
	std::cout << "     " << "=               -f2 FILE4     --> Input file (Point; ascii or csv file with 1 point, 1 line)   ="<<std::endl;
    std::cout << "     " << "=               -d disttype   --> Flag:                                                        ="<<std::endl;
    std::cout << "     " << "=                                   (0) Frechet distance                                       ="<<std::endl;
    std::cout << "     " << "=                                   (1) Aspect ratio distances (based on ellipses)             ="<<std::endl;
    std::cout << "     " << "=                                   (2) Curvature-based distances                              ="<<std::endl;
    std::cout << "     " << "=                                   (3) 2-Wasserstein distance between persistence diagrams    ="<<std::endl;
    std::cout << "     " << "=                                   (4) All                                                    ="<<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "     " << "================================================================================================"<<std::endl;
    std::cout << "\n\n" <<std::endl;
}

/* ===============================================================================================
   Parse Argument from command line:

   =============================================================================================== */

bool parse_args(int argc, char **argv, std::string *file1, std::string *file2, int *disttype)
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

			if (param == "-i1") {
				*file1 = argv[i + 1];
			}
			if (param == "-i2") {
				*file2 = argv[i + 1];
			}
			if (param == "-d") {
				*disttype = std::atoi(argv[i + 1]);
			}
			if (param == "-f1") {
				*focal1 = argv[i + 1];
			}
			if (param == "-f2") {
				*focal2 = argv[i + 1];
			}
		}
  	}
	return true;
}
