/*!
* @file Vertex.h
* @author Patrice Koehl
* @author Yossi Bokor Bleile
* @date April 2023
* @version 1
* @copyright BSD 3-Clause License.
*/


/*********************************************************************************
 *	The Vertex class
 *********************************************************************************/

#ifndef VERTEX_H
#define VERTEX_H

#include "Types.h"
#include "Vector2D.h"

/*********************************************************************************
  Vertex class
 *********************************************************************************/
namespace correa {
  class Vertex {
	public:

		// Constructor
		Vertex();

		// position
		Vector2D position;

		// id between 0 and |V|-1
		int index;

		// previous vertex
		VertexIter prev;

		// next vertex
		VertexIter next;

		// return degree
		int degree();

		// returns exterior angle. Note: only valid for boundary vertices
		double exteriorAngle();
		
		//returns the distance form the origin
		double height();

		double x() {
			return position.x;
		};

		double y() {
			return position.y;
		};

		void shift(double c_x, double c_y) {
			position -= Vector2D(c_x, c_y);
		}
  }; 

/*********************************************************************************
  Constructor
 *********************************************************************************/

  Vertex::Vertex():
  index(-1)
  {

  }

/*********************************************************************************
  Degree of a vertex
 *********************************************************************************/

 int Vertex::degree()
 {
	int k = 2;

	return k;
  }

/*********************************************************************************
  Exterior angle at a vertex
 *********************************************************************************/

 double Vertex::exteriorAngle()
 {
	Vector2D a, b, c, u, v;
	double val1, val2, angle;

	a = prev->position;
	b = position;
	c = next->position;
	u = b - a;
	v = c - b;
	val1 = u[0]*v[1] - u[1]*v[0];
	val2 = u[0]*v[0] + u[1]*v[1];
	angle = std::atan2(val1, val2);

	return angle;
  };

/* Height of the vertex */
double Vertex::height(){
	std::cout << "position[0] is " << position[0] << std::endl;
	std::cout << "position[1] is " << position[1] << std::endl;
	double height = sqrt(pow(position[0],2) + pow(position[1], 2));
	std::cout << "height is " << height << std::endl;
	return height;
};

} //end namespace correa
#endif