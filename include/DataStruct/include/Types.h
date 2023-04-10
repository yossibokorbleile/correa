/*
 	Types.h

 	Authors: 	Patrice Koehl, Department of Computer Science, University of California, Davis
				Yossi Bokor Bleile, Department of Mathematical Sciences, University of Aalborg, Aalborg
 	Date: April 2023
	Version: 1
*/


/*********************************************************************************
  Define types for iterators
 *********************************************************************************/

#ifndef TYPES_H
#define TYPES_H

#define _USE_MATH_DEFINES

#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include "math.h"
#include "Vector2D.h"

namespace correa {
class Vertex;
class Polygon;

typedef std::vector<Vertex>::iterator            VertexIter;
typedef std::vector<Vertex>::const_iterator      VertexCIter;
} //end namespace 
#endif
