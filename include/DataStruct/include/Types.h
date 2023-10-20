/*!
* @file Types.h
* @author Patrice Koehl
* @author Yossi Bokor Bleile
* @date April 2023
* @version 1
* @copyright BSD 3-Clause License.
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
