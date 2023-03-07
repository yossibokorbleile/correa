/* ===============================================================================================
   correa_bindings.h

   Author:  Yossi Bokor Bleile
   Date:    March 03 2023
   Version: 0.1
   =============================================================================================== */

#ifndef _CORREABINDINGS_H_
#define _CORREABINDINGS_H_

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
#include <algorithm>


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
//#include "Ellipse.h"
//#include "OT1.h"
//#include "PH0.h"
#include "PersistenceDiagram.h"
//#include "Curvature.h"
//#include "Component.h"



#endif