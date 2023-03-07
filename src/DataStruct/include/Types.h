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

class Vertex;
class Polygon;

typedef std::vector<Vertex>::iterator            VertexIter;
typedef std::vector<Vertex>::const_iterator      VertexCIter;

#endif
