/* ===============================================================================================
   Comp2DShapes.h

   Author:  Patrice Koehl & Yossi Bokor Bleile
   Date:    April 4, 2023
   Version: 1
   =============================================================================================== */

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
#include "PersistenceDiagram.h"



INOUT inout;
Polygon poly;
PolygonBuilder pbuilder;

//PDInOut PDIO;
/* ===============================================================================================
   Prototypes
   =============================================================================================== */

static void usage(char** argv);
bool parse_args(int argc, char **argv, std::string *cell1, std::string *focal1, std::string *cell2, std::string *focal2);

#endif


