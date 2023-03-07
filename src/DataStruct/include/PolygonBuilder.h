/* ===========================================================================================
 *
 * Class that contains tools to manipulate a 2D polygon
 *
 ==============================================================================================*/

#ifndef _POLYGONBUILDER_H
#define _POLYGONBUILDER_H

#include "Polygon.h"
#include <set>

/* ===========================================================================================
 The MeshBuilder class
 ==============================================================================================*/

  class PolygonBuilder {

	public:

		// clean polygon
		void clean_points(int* Np, double *coord);

		// Generate polygon structure
		bool buildPolygon(int npoint, double *coord, Polygon& polygon);

		// triangulate the polygon
		bool triangulatePolygon(Polygon& polygon, int *ntrig, int *triangles);

		// List boundary edges
		int boundaryEdges(Polygon& polygon);

		// List triangle edges
		int trigEdges(Polygon& polygon, int ntrig, int *indices);

	private:

		// index all elements in a polygon structure
		void indexElements(Polygon& polygon);

		 // compute the angle (in degree) defined by three points
		double angle_degree(Vector2D A, Vector2D B, Vector2D C);

		// true if C is between A and B
		bool between(Polygon& polygon, int ia, int ib, int ic);

		// returns a measure of collinearity for 3 points
		bool collinear(Polygon& polygon, int ia, int ib, int ic);

		// checks if the segment I to J is a proper internal diagonal
		bool diagonal(Polygon& polygon, int I, int J, int *prev_node,
		int *next_node);

		// checks if the segment I to J is a proper internal diagonal
		bool diagonalie(Polygon& polygon, int I, int J, int *next_node);

		bool in_cone(Polygon& polygon, int I, int J, int *prev_node,
                int *next_node);

		// true if lines A:B and C:D intersect
		bool intersect(Polygon& polygon, int A, int B, int C, int D);

		bool intersect_prop(Polygon& polygon, int A, int B, int C, int D);

		// signed area of a triangle
		double triangle_area(Polygon& polygon, int A, int B, int C);

  };

/* ===========================================================================================
 * Generate the polygon
 *
 ==============================================================================================*/

  bool PolygonBuilder::buildPolygon(int npoint, double *Coord, Polygon& polygon)
  {

	// Preallocate polygon elements
	polygon.vertices.clear(); 
	polygon.vertices.reserve(npoint); 

	// create and insert vertices
	std::vector<VertexIter> indexToVertex(npoint);

	double area = 0.0;

	int im1 = npoint-1;
	for(int i = 0; i < npoint; i++) {
		area += Coord[2*im1] * Coord[2*i+1] - Coord[2*i] * Coord[2*im1+1];
		im1 = i;
	}
	area *= 0.5;

	if(area >0) {
		for(int i = 0; i < npoint; i++) {
			VertexIter v = polygon.vertices.insert(polygon.vertices.end(), Vertex());
			v->position.x = Coord[2*i];
			v->position.y = Coord[2*i+1];
			indexToVertex[i] = v;
		}


	} else {
		for(int i = 0; i < npoint; i++) {
			VertexIter v = polygon.vertices.insert(polygon.vertices.end(), Vertex());
			v->position.x = Coord[2*(npoint-1-i)];
			v->position.y = Coord[2*(npoint-1-i)+1];
			indexToVertex[i] = v;
		}
	}



	VertexIter v;
	im1 = npoint-1;
	for(int i = 0; i < npoint-1; i++) {
		v = indexToVertex[i];
		v->prev = indexToVertex[im1];
		v->next = indexToVertex[i+1];
		im1 = i;
	}
	v = indexToVertex[npoint-1];
	v->prev = indexToVertex[npoint-2];
	v->next = indexToVertex[0];

	// index elements
	indexElements(polygon);

	return true;

  }

/* ===========================================================================================
 * Index all elements in a Mesh2D data structure
 ==============================================================================================*/

  void PolygonBuilder::indexElements(Polygon& polygon)
  {
	int index = 0;
	for (VertexIter v = polygon.vertices.begin(); v != polygon.vertices.end(); v++) {
		v->index = index++;
	}

  }

/* ===========================================================================================
 * Computes the angle (in degree) between three points:
 *
 *       IA
 *      /
 *     /
 *    /
 *   IB ------IC
 *
 ==============================================================================================*/

  double PolygonBuilder::angle_degree(Vector2D A, Vector2D B, Vector2D C)
  {
	Vector2D u,v;

	u = A-B; u.normalize();
	v = C-B; v.normalize();

	double angle = std::acos(std::max(-1.0, std::min(1.0, dot(u, v))));

	return angle*180./M_PI;

  }

/* ===========================================================================================
 * Computes the signed area of a triangle defined by three points:
 *
 *       IA
 *      /  \
 *     /    \
 *    /      \
 *   IB ------IC
 *
 ==============================================================================================*/

  double PolygonBuilder::triangle_area(Polygon& polygon, int ia, int ib, int ic)
  {
	Vector2D A = polygon.vertices[ia].position;
	Vector2D B = polygon.vertices[ib].position;
	Vector2D C = polygon.vertices[ic].position;

	double value = (B.x-A.x)*(C.y-A.y) - (C.x-A.x)*(B.y-A.y);

	return 0.5*value;
  }
	
	
/* ===========================================================================================
 * true if C is between A and B
 *
 * Note:
 *	The points must be (numerically) collinear.
 ==============================================================================================*/

  bool PolygonBuilder::between(Polygon& polygon, int ia, int ib, int ic)
  {

	Vector2D A = polygon.vertices[ia].position;
	Vector2D B = polygon.vertices[ib].position;
	Vector2D C = polygon.vertices[ic].position;

	bool value;
	double xmax, xmin, ymax, ymin;

	if(! collinear(polygon, ia, ib, ic)) {
		value = false;
	} else if (std::abs(A.y - B.y) < std::abs(A.x-B.x)) {
		xmax = std::max(A.x, B.x);
		xmin = std::min(A.x, B.x);
		value = (xmin <= C.x) && (C.x <= xmax);
	} else {
		ymax = std::max(A.y, B.y);
		ymin = std::min(A.y, B.y);
		value = (ymin <= C.y) && (C.y <= ymax);
	}

	return value;
  }

/* ===========================================================================================
 * returns a measure of collinearity for 3 points
 *
 * See John Burkardt's code for a discussion of the implementation
 ==============================================================================================*/

  bool PolygonBuilder::collinear(Polygon& polygon, int ia, int ib, int ic)
  {

	Vector2D A = polygon.vertices[ia].position;
	Vector2D B = polygon.vertices[ib].position;
	Vector2D C = polygon.vertices[ic].position;

	double eps = 2.220446049250313E-016;

	double area = triangle_area(polygon, ia, ib, ic);

	double side_ab_sq = (A.x-B.x)*(A.x-B.x) + (A.y-B.y)*(A.y-B.y);
	double side_bc_sq = (B.x-C.x)*(B.x-C.x) + (B.y-C.y)*(B.y-C.y);
	double side_ca_sq = (C.x-A.x)*(C.x-A.x) + (C.y-A.y)*(C.y-A.y);

	double size_max_sq = std::max(side_ab_sq, std::max(side_bc_sq, side_ca_sq));

	bool value;

	if(size_max_sq <= eps) {
		value = true;
	} else if(2.0*std::abs(area) <= eps*size_max_sq) {
		value = true;
	} else {
		value = false;
	}

	return value;
  }

/* ===========================================================================================
 * check if a segment is a proper diagonal
 *
 ==============================================================================================*/

  bool PolygonBuilder::diagonal(Polygon& polygon, int I, int J, int *prev_node,
		int *next_node)
  {

  	bool value1 = in_cone(polygon, I, J, prev_node, next_node);
  	bool value2 = in_cone(polygon, J, I, prev_node, next_node);
	bool value3 = diagonalie(polygon, I, J, next_node);

	bool value = value1 && value2 && value3;

	return value;
  }

/* ===========================================================================================
 * check if a segment is a proper diagonal
 *
 ==============================================================================================*/

  bool PolygonBuilder::diagonalie(Polygon& polygon, int I, int J, int *next_node)
  {

	int first = I;
	int j = first;
	int jp1 = next_node[first];

	bool value = true;

	while (1)
	{
		// Skip any edge that includes vertex I or J

		if(j== I || j == J || jp1 == I || jp1 == J) {
		} else {
			bool value2 = intersect(polygon, I, J, j, jp1);
			if(value2) {
				value = false;
				break;
			}
		}
		j = jp1;
		jp1 = next_node[j];

		if(j==first) break;
	}

	return value;
  }

/* ===========================================================================================
 * true if the diagonal I:J is strictly internal
 *
 ==============================================================================================*/

  bool PolygonBuilder::in_cone(Polygon& polygon, int I, int J, int *prev_node,
                int *next_node)
  {

	int IM1 = prev_node[I];
	int IP1 = next_node[I];

	double t1 = triangle_area(polygon, I, IP1, IM1);
	double t2, t3;
	bool value;

	if(0.0 <= t1) {
		t2 = triangle_area(polygon, I, J, IM1);
		t3 = triangle_area(polygon, J, I, IP1);
		value = (0.0 < t2) && (0.0<t3);
	} else {
		t2 = triangle_area(polygon, I, J, IP1);
		t3 = triangle_area(polygon, J, I, IM1);
		value = (0.0 <= t2) && (0.0<=t3);
		value = !value;
	}

	return value;
  }

/* ===========================================================================================
 * true if lines A:B and C:D intersect
 *
 ==============================================================================================*/

   bool PolygonBuilder::intersect(Polygon& polygon, int A, int B, int C, int D) 
   {

	bool value;

	if(intersect_prop(polygon, A, B, C, D)) {
		value = true;
	} else if (between(polygon, A, B, C) ) {
		value = true;
	} else if (between(polygon, A, B, D) ) {
		value = true;
	} else if (between(polygon, C, D, A) ) {
		value = true;
	} else if (between(polygon, C, D, B) ) {
		value = true;
	} else {
		value = false;
	}

	return value;
   }

/* ===========================================================================================
 * true if lines A:B and C:D intersect
 *
 ==============================================================================================*/

   bool PolygonBuilder::intersect_prop(Polygon& polygon, int A, int B, int C, int D) 
   {

	bool value;
	if(collinear(polygon, A, B, C)) {
		value = false;
	} else if (collinear(polygon, A, B, D) ) {
		value = false;
	} else if (collinear(polygon, C, D, A) ) {
		value = false;
	} else if (collinear(polygon, C, D, B) ) {
		value = false;
	} else {

		double t1 = triangle_area(polygon, A, B, C);
		double t2 = triangle_area(polygon, A, B, D);
		double t3 = triangle_area(polygon, C, D, A);
		double t4 = triangle_area(polygon, C, D, B);

		bool value1 = 0.0 < t1;
		bool value2 = 0.0 < t2;
		bool value3 = 0.0 < t3;
		bool value4 = 0.0 < t4;

		value = (value1^value2) && (value3^value4);
	}

	return value;
   }

/* ===========================================================================================
 * Generate a triangulation of a polygon
 *
 * Note: there are N-2 triangles in the triangulation, where N is the number of vertices defining
 *       the polygon
 *       the first edge listed is always an internal diagonal
 *
 ==============================================================================================*/

  bool PolygonBuilder::triangulatePolygon(Polygon& polygon, int *ntrig, int *triangle)
  {

	int Npoint = polygon.vertices.size();

	// We must have at least 3 vertices
	if(Npoint < 3) {
		return false;
	}

	int *prev_node = new int[Npoint];
	int *next_node = new int[Npoint];

	int ip = 0;
	prev_node[ip] = Npoint-1;
	next_node[ip] = ip + 1;

	for(int i = 1; i < Npoint-1; i++) {
		prev_node[i] = i-1; 
		next_node[i] = i+1; 
	}

	ip = Npoint-1;
	prev_node[ip] = ip - 1;
	next_node[ip] = 0;

	bool *ear = new bool[Npoint];

	for(int i = 0; i < Npoint; i++) {
		ear[i] = diagonal(polygon, prev_node[i], next_node[i], prev_node, next_node);
	}

	int Ntrig = 0;

	int i2 = 0;
	int i0, i1, i3, i4;

	while(Ntrig < Npoint-3)
	{
		if(ear[i2]) {
			i3 = next_node[i2];
			i4 = next_node[i3];
			i1 = prev_node[i2];
			i0 = prev_node[i1];

			next_node[i1] = i3;
			prev_node[i3] = i1;

			ear[i1] = diagonal(polygon, i0, i3, prev_node, next_node);
			ear[i3] = diagonal(polygon, i1, i4, prev_node, next_node);

			triangle[0+3*Ntrig] = i3;
			triangle[1+3*Ntrig] = i1;
			triangle[2+3*Ntrig] = i2;
			Ntrig++;
		}

		i2 = next_node[i2];
	}

	i3 = next_node[i2];
	i1 = prev_node[i2];

	triangle[0+3*Ntrig] = i3;
	triangle[1+3*Ntrig] = i1;
	triangle[2+3*Ntrig] = i2;
	Ntrig++;

	delete [] ear; delete [] next_node; delete [] prev_node;

	*ntrig = Ntrig;

	return true;
  }

/* ===========================================================================================
 * clean up polygon: 
 *	- if two consecutive vertices I and I+1 are equal, remove I+1
 *	- if a vertex I has angle 0, remove I
 *
 ==============================================================================================*/

  void PolygonBuilder::clean_points(int* Np, double *coord)
  {

	int Npoint = *Np;
	int Npoint2;

	std::vector<Vector2D> polygon;

	int nchanges = 0;
	double angle_tol = 1.e-4;
	double angle;
	int im1, ip1;

	std::cout << std::endl;
	std::cout << "Initial number of points in polygon      : " << Npoint << std::endl;

	do {

		nchanges = 0;
		polygon.clear();

		im1 = Npoint-1;
		for(int i = 0; i < Npoint; i++) {
			if(coord[2*i] != coord[2*im1] || coord[2*i+1] != coord[2*im1+1]) {
				Vector2D v;
				v.x = coord[2*i]; v.y = coord[2*i+1];
				polygon.push_back(v);
			} else {
				nchanges++;
			}
			im1 = i;
		}

		Npoint2 = polygon.size();

		im1 = Npoint2-1;
		Npoint = 0;
		for(int i = 0; i < Npoint2; i++) {
			ip1 = (i+1) %Npoint2;

			angle = angle_degree(polygon[im1], polygon[i], polygon[ip1]);

			if(std::abs(angle) <= angle_tol || std::abs(angle-180) <= angle_tol ) {
				nchanges++;
			} else {
				coord[2*Npoint] = polygon[i].x;
				coord[2*Npoint+1] = polygon[i].y;
				Npoint++;
			}
			im1 = i;
		}

	} while (nchanges !=0);

	*Np = Npoint;
	std::cout << "Number of points in polygon after cleanup: " << Npoint << std::endl;
	std::cout << std::endl;

  }

/* ===========================================================================================
 * Get list of boundary edges
 ==============================================================================================*/

  int PolygonBuilder::boundaryEdges(Polygon& polygon)
  {
	int idx1, idx2;
	for (VertexIter v = polygon.vertices.begin(); v != polygon.vertices.end(); v++)
	{
		idx1 = v->index;
		idx2 = v->next->index;
		polygon.edges.push_back(std::make_pair(idx1, idx2));
	}

	int nedges = polygon.edges.size();

	return nedges;

  }

/* ===========================================================================================
 * Get list of edges in triangulation
 ==============================================================================================*/

  int PolygonBuilder::trigEdges(Polygon& polygon, int ntrig, int *indices)
  {
	int npoints = polygon.vertices.size();

	std::vector<std::set<int> > data;

	data.resize(npoints);

	// build table
	for (int it = 0; it < ntrig; it++) {
		int i = indices[3*it];
		int j = indices[3*it+1];
		int k = indices[3*it+2];

		if(i<j) {
			data[i].insert(j);
		} else {
			data[j].insert(i);
		}
		if(i<k) {
			data[i].insert(k);
		} else {
			data[k].insert(i);
		}
		if(j<k) {
			data[j].insert(k);
		} else {
			data[k].insert(j);
		}
	}

	for(int i = 0; i < npoints; i++) {
		for (std::set<int>::iterator it = data[i].begin(); it != data[i].end(); it++) {
			polygon.edges.push_back(std::make_pair(i, *it));
		}
	}

	int nedges = polygon.edges.size();

	return nedges;

  }
#endif
