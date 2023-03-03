/* ===============================================================================================
   Comp2DShapes.h

   Author:  Patrice Koehl
   Date:    3/18/2021
   Version: 1
   =============================================================================================== */

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

#include "InOut.h"
#include "PolygonBuilder.h"
#include "Polygon.h"
#include "Frechet.h"
#include "Ellipse.h"

#include "OT1.h"
OT1 ot1;

#include "Curvature.h"
#include "PH0.h"
#include "PersistenceDiagram.h"
#include "PDInOut.h"
#include "Curvature.h"



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
bool parse_args(int argc, char **argv, std::string *FILE1, std::string *FILE2, int *disttype);

#endif
