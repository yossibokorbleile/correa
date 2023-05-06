/*
 	Comp2DShape.h

 	Authors: 	Patrice Koehl, Department of Computer Science, University of California, Davis
				Yossi Bokor Bleile, Department of Mathematical Sciences, University of Aalborg, Aalborg
 	Date: April 2023
	Version: 1
*/

#ifndef _COMP2DSHAPES_H_
#define _COMP2DSHAPES_H_

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

/*================================================================================================
 Definitions for multi-threading
================================================================================================== */

#define NUM_THREADS 32

int threadids[NUM_THREADS];
pthread_t threads[NUM_THREADS];

/* ===============================================================================================
   Local includes
   =============================================================================================== */

#include "hera/wasserstein.h"
#include "InOut.h"
#include "PolygonBuilder.h"
#include "Polygon.h"
#include "Frechet.h"
#include "Ellipse.h"
#include "Curvature.h"
#include "PH0.h"
#include "Curvature.h"

using namespace correa;
using PersistenceDiagram = std::vector<std::pair<double,double>>;

INOUT inout;
Polygon poly;
PolygonBuilder pbuilder;
Frechet frechet;
Ellipse ellipse;
Curvature curv;

/* ===============================================================================================
   Prototypes
   =============================================================================================== */

static void usage(char** argv);
bool parse_args(int argc, char **argv, std::string *file1, std::string *file2, int *disttype);

#endif
