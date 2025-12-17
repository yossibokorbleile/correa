/*!
* @file Comp2DShapesFocal.h
* @author Patrice Koehl
* @author Yossi Bokor Bleile
* @date December 2025
* @version 1.2
* @copyright BSD 3-Clause License.
*/
#ifndef _COMP2DSHAPESFOCAL_H_
#define _COMP2DSHAPESFOCAL_H_

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
bool parse_args(int argc, char **argv, std::string *file1, std::string *focal1, std::string *file2, std::string *focal2, int *disttype, double *microns_per_pixel1, double *microns_per_pixel2);

#endif
