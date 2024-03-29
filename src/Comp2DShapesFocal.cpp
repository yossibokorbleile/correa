/*!
* @file	Comp2DShapesFocal.cpp
* @brief compare polygons with specific focal points
* @author Patrice Koehl
* @author Yossi Bokor Bleile
* @date April 2023
* @version 1
* @copyright BSD 3-Clause License.
* @copyright BSD 3-Clause License.
*/



/* ===============================================================================================
   Includes
   =============================================================================================== */

#include "Comp2DShapesFocal.h"


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

        if (!parse_args(argc, argv, &INfile1, &INfocal1, &INfile2, &INfocal2, &disttype)) return 1;

/*	==========================================================================================
	Read in the polygon1 from input file1
	========================================================================================== */

	int ndim;
	int npoint1;

	double *X1;
	X1 = nullptr;
	inout.read(INfile1, &ndim, &npoint1, &X1);

	std::vector<double> focal_1;
   	std::ifstream in1(INfocal1);
   	std::string line;
   	double val;

   	while (getline(in1, line)){
			size_t start_pos = 0;
			while ((start_pos = line.find(",", start_pos)) != std::string::npos)
			{
				line.replace(start_pos, 1, " ");
				start_pos += 1; // Handles case where 'to' is a substring of 'from'
			}
       		istringstream iss(line);
       		while (iss >> val){
				std::cout << "val is " << val << std::endl;
          		focal_1.push_back(val);
   		};
	};
	Vector2D focal1(focal_1[0], focal_1[1]);

/*	==========================================================================================
	Store polygon1 as a polygon
	========================================================================================== */

	Polygon polygon1;
	pbuilder.clean_points(&npoint1, X1);
	pbuilder.buildPolygon(npoint1, X1, polygon1);

	polygon1.shift(focal1);

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

	std::vector<double> focal_2;
   	std::ifstream in2(INfocal2);

   	while (getline(in2, line)){
			size_t start_pos = 0;
			while ((start_pos = line.find(",", start_pos)) != std::string::npos)
			{
				line.replace(start_pos, 1, " ");
				start_pos += 1; // Handles case where 'to' is a substring of 'from'
			}
       		istringstream iss(line);
       		while (iss >> val){
				std::cout << "val is " << val << std::endl;
          		focal_2.push_back(val);
   		};
	};
	Vector2D focal2(focal_2[0], focal_2[1]);

/*	==========================================================================================
	Store polygon2 as a polygon
	========================================================================================== */

	Polygon polygon2;
	pbuilder.clean_points(&npoint2, X2);
	pbuilder.buildPolygon(npoint2, X2, polygon2);

	polygon2.shift(focal2);
	// Center polygon
	iscale = 0;
	range = 100;
	polygon2.centerScale(range,iscale);

/*	==========================================================================================
	Compute distances
	========================================================================================== */

	double dFrechet, dE_M, dE_m, dE_l, dW, dcOT1;
	double a1_M, b1_M, r1_M, a2_M, b2_M, r2_M;
	double a1_m, b1_m, r1_m, a2_m, b2_m, r2_m;
	double a1_l, b1_l, r1_l, a2_l, b2_l, r2_l;
	double willmore1, willmore2;
	correa::PH0 f1(polygon1.vertices);
	correa::PH0 f2(polygon2.vertices);

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
		f1.Persistence();
		f2.Persistence();
		hera::AuctionParams<double> hera_params;
    	hera_params.max_num_phases = 800;
		hera_params.wasserstein_power = 2;
		hera_params.internal_p = 2;
		const std::vector<std::pair<double,double>> pd1 = f1.persistence_diagram();
		const std::vector<std::pair<double,double>> pd2 = f2.persistence_diagram();
		const hera::AuctionParams<double> params = hera_params;
		//double w_q_dist = hera::wasserstein_dist<std::vector<std::pair<double,double>>>(pd1, pd2, params);
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
		f1.Persistence();
		f2.Persistence();
		hera::AuctionParams<double> hera_params;
    	hera_params.max_num_phases = 800;
		hera_params.wasserstein_power = 2;
		hera_params.internal_p = 2;
		const std::vector<std::pair<double,double>> pd1 = f1.persistence_diagram();
		const std::vector<std::pair<double,double>> pd2 = f2.persistence_diagram();
		const hera::AuctionParams<double> params = hera_params;
		double w_q_dist = hera::wasserstein_dist<std::vector<std::pair<double,double>>>(pd1, pd2, params);
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
		std::cout << "Persistence diagram of polygon 1  : ";
		std::cout<< "There are " << f1.persistence_diagram().size() << " points in the persistence diagram." << std::endl;
		for (int i = 0; i < f1.persistence_diagram().size(); i++){
			std::cout << "(" << get<0>(f1.persistence_diagram()[i]) << ", " << get<1>(f1.persistence_diagram()[i]) << ")" << std::endl;
		};
	}

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
		std::cout << "Persistence diagram of polygon 2  : ";
		std::cout<< "There are " << f2.persistence_diagram().size() << " points in the persistence diagram." << std::endl;
		for (int i = 0; i < f2.persistence_diagram().size(); i++){
			std::cout << "(" << get<1>(f2.persistence_diagram()[i]) << ", " << get<1>(f2.persistence_diagram()[i]) << ")" << std::endl;
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
		std::cout << "Distance (2-Wasserstein between persitence diagrams): i need to figure out how hera works" <<  std::endl;
	
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
    std::cout << "     " << "=          Comp2DShapes -i1 FILE1 -f1 FOCAL1 -i2 FILE2 -f2 FILE2 -d disttype                   ="<<std::endl;
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

bool parse_args(int argc, char **argv, std::string *file1, std::string *focal1, std::string *file2, std::string *focal2, int *disttype)
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
