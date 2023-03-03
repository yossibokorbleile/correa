/* ===============================================================================================
   Clustering.h

   Author:  Yossi Bokor
   Date:    February 2, 2022
   Version: 1
   =============================================================================================== */

#ifndef _CLUSTERING_H_
#define _CLUSTERING_H_

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
#include "InOut.h"
#include "kMeans.h"

INOUT InOut;
/*================================================================================================
 Definitions for multi-threading
================================================================================================== */

#define NUM_THREADS 32

int threadids[NUM_THREADS];
pthread_t threads[NUM_THREADS];

/* ===============================================================================================
   Local includes
   =============================================================================================== */

//#include "InOut.h"
#include "kMeans.h"

//INOUT inout;

/* ===============================================================================================
   Prototypes
   =============================================================================================== */

static void usage(char** argv);
bool parse_args(int argc, char **argv, std::string *infile, int *nobj, int *nclus);

#endif
