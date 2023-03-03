/*
File to test various bits of code.
*/

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

#define NUM_THREADS 32

int threadids[NUM_THREADS];
pthread_t threads[NUM_THREADS];

#include "InOut.h"
//#include "PDInOut.h"
#include "Polygon.h"
#include "PolygonBuilder.h"
#include "Vector2D.h"
#include "Vertex.h"
#include "Ellipse.h"
#include "Frechet.h"
#include "MatInOut.h"
#include "OT1.h"
OT1 ot1;

#include "PH0.h"
#include "PersistenceDiagram.h"

#include "Curvature.h"


#include <iostream>
#include <vector>

INOUT inout;
Polygon poly;
PolygonBuilder pbuilder;
Frechet frechet;
Ellipse ellipse;
Curvature curv;
MATINOUT mio;
using namespace std;

int main() //(int argc, char **argv)
{
	/*
	string INfile1;
	string INfile2;
	INfile1 = "/Users/yossibokor/USYD/ToMaCo/24hrs_hMSCs/part1/cell/csv_contour/001_contour.csv";
	INfile2 = "/Users/yossibokor/USYD/ToMaCo/24hrs_hMSCs/part1/cell/csv_contour/003_contour.csv";
	
	int ndim;
	int npoint1;

	double *X1;
	X1 = nullptr;
	inout.read(INfile1, &ndim, &npoint1, &X1);
	
	Polygon polygon1;
	pbuilder.clean_points(&npoint1, X1);
	pbuilder.buildPolygon(npoint1, X1, polygon1);

	// Center polygon
	int iscale = 0;
	double range = 100;
	polygon1.centerScale(range,iscale);
	
	int npoint2;

	double *X2;
	X2 = nullptr;
	inout.read(INfile2, &ndim, &npoint2, &X2);
	
	Polygon polygon2;
	pbuilder.clean_points(&npoint2, X2);
	pbuilder.buildPolygon(npoint2, X2, polygon2);

	// Center polygon
	iscale = 0;
	range = 100;
	polygon2.centerScale(range,iscale);
	
	PH0 f1(polygon1.vertices);
	PH0 f2(polygon2.vertices);
	cout << "Got the filtrations." << endl;
	cout << "Doing first persistence." << endl;
	f1.Persistence();
	cout << "Doing second persistence." << endl;
	f2.Persistence();
	cout << "Generated the persistence diagrams." << endl;
	
	PersistenceDiagram pd1(f1.pd);
	PersistenceDiagram pd2(f2.pd);
	cout << "First point in pd1 is (" << pd1.points[0][0] << ", " << pd1.points[0][1] <<") and last point in pd1 is (" << pd1.points.back()[0] << ", " << pd1.points.back()[1] << ")." << endl;
	cout << "Last point in pd2 is (" << pd2.points[0][0] << ", " << pd2.points[0][1] <<") and last point in pd2 is (" << pd2.points.back()[0] << ", " << pd2.points.back()[1] << ")." << endl;
	
	
	cout << "The first persistence diagram is: " << endl;
	for (int i = 0; i < pd1.np; i++){
		vector<double> pt;
		pt = pd1.points[i];
		cout << "point " << i+1 << " out of " << pd1.np << " is (" << pt[0] << ", " << pt[1] << ")" << endl;
	};

	cout << "\nand the second persistence diagram is: " << endl;;
	for (int i = 0; i < pd2.np; i++){
		vector<double> pt;
		pt = pd2.points[i];
		cout << "(" << pt[0] << ", " << pt[1] << ")" << endl;
	};
	cout << " " << endl;
	
	*/

	PersistenceDiagram pd1, pd2;
	vector<double> x0{1,2};
	vector<double> x1{3,4};
	vector<double> x2{4,5};
	vector<double> x3{6,7};
	vector<double> x4{8,9};
	pd1.addPoint(x0);
	pd1.addPoint(x1);
	pd1.addPoint(x2);
	pd2.addPoint(x0);
	pd2.addPoint(x1);
	//pd2.addPoint(x2);
	cout << "The persistence diagrams are: " << endl;
	cout << "pd1: " << endl;
	for (int i = 0; i < pd1.np; i++){
		cout << "(" << pd1.points[i][0] << ", " << pd1.points[i][1] << ")" << endl;
	};
	cout << "pd2: " << endl;
	for (int i = 0; i < pd2.np; i++){
		cout << "(" << pd2.points[i][0] << ", " << pd2.points[i][1] << ")" << endl;
	};
	cout << "pd1 has " << pd1.np << " points and pd2 has " << pd2.np << " points." << endl;

	//double pdWasserstein;
	//pdWasserstein = WassersteinDistance(f1.pd, f2.pd);
	cout << scientific << WassersteinDistance(pd2, pd2) << " is he 2-Wasserstein distance between pd2 and pd2 is " << endl;
	cout << scientific << WassersteinDistance(pd1, pd1) << " is he 2-Wasserstein distance between pd1 and pd1 is " << endl;
	cout << scientific << WassersteinDistance(pd1, pd2) << " is he 2-Wasserstein distance between pd1 and pd2 is " << endl;
	cout << scientific << WassersteinDistance(pd2, pd1) << " is he 2-Wasserstein distance between pd2 and pd1 is " << endl;
	
	
	return 0;

};
