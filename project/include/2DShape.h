/* ===============================================================================================
   2DShape.h

   Author:  Patrice Koehl
   Date:    3/18/2021
   Version: 1
   =============================================================================================== */

#ifndef _2DSHAPE_H_
#define _2DSHAPE_H_

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
#include "Ellipse.h"
#include "OT1.h"

OT1 ot1;

#include "PH0.h"
#include "PersistenceDiagram.h"
#include "Curvature.h"
#include "Component.h"


INOUT inout;
Polygon poly;
PolygonBuilder pbuilder;
Ellipse ellipse;
Curvature curv;
PersistenceDiagram PD;


/* ===============================================================================================
   Prototypes
   =============================================================================================== */

static void usage(char** argv);
bool parse_args(int argc, char **argv, std::string *INFILE);

#endif
