/*!
* @file Comp2D.h
* @author Patrice Koehl
* @author Yossi Bokor Bleile
* @date April 2023
* @version 1
* @copyright BSD 3-Clause License.
*/

#ifndef _COMP2D_H_
#define _COMP2D_H_

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
#include <hera/wasserstein.h>

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
#include "PH0.h"
#include "Frechet.h"
#include "Ellipse.h"
//#include "OT1.h"
#include "Curvature.h"
#include "hera/wasserstein.h"

using namespace correa;

INOUT inout;
Polygon poly;
PolygonBuilder pbuilder;
Frechet frechet;
Ellipse ellipse;
Curvature curv;

//PDInOut PDIO;
/* ===============================================================================================
   Prototypes
   =============================================================================================== */

static void usage(char** argv);
bool parse_args(int argc, char **argv, std::string *file1, std::string *file2, std::string *focal1, std::string *focal2, int *disttype);

#endif


