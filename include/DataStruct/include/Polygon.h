/*!
* @file Polygon.h
* @author Patrice Koehl
* @author Yossi Bokor Bleile
* @date December 2025
* @version 1.3
* @copyright BSD 3-Clause License.
*/

/* ===========================================================================================
 *
 * Generate polygon
 *
 ==============================================================================================*/

#ifndef _POLYGON_H
#define _POLYGON_H

#include "Types.h"
#include "Vertex.h"
#include "Vector2D.h"
/* ===========================================================================================
 * Class for polygon
 ==============================================================================================*/
namespace correa {

	/*!
	* Polygon class
	* Stores the polygon and calculates various information about it.
	*/
  	class Polygon {	

		public:
			// constructor
			Polygon() = default;

			// copy constructor
			Polygon(Polygon& polygon);

			// assignment operator
			Polygon& operator=(const Polygon& other);

			std::vector<Vertex> vertices;
			std::vector<double > bdLength0;
			std::vector<double > bdLength;
			double distorsion_; 
			double originalArea_;

			std::vector<std::pair<int, int> > edges;

			// area of the polygon
			double area();

			// unscaled area of the polygon
			double originalArea();

			// area of the polygon
			double length();

			// Center and scale the polygon
			void centerScale(double range, int iscale);

			// Center and scale the polygon by area
			void centerScaleArea();

			// Scale the polygon by area
			void scaleArea();


			void boundaryLength0();
			void boundaryLength();

			void distorsion();

			void shift(Vector2D center, bool verbose);

			int size() {
				return vertices.size();
			};

			void labelPolygon(const double EPSILON);

			// double signedDistanceToOutline(const Vector2D& point);

			Vector2D calculateCentroid();

			bool isPointInside(const Vector2D& point) const;

		private:
			Vector2D label;
	};

 /* ===========================================================================================
 * Computes the area inside a curve using the shoelace formula
 ==============================================================================================*/
	/*!
	* Calculates the area of polygon.
	*/
	double Polygon::area()
	{
		// for (const auto& v : vertices) {
		// 	std::cerr << "vertex: (" << v.position.x << ", " << v.position.y << ")" << std::endl;
		// }
		double surface = 0.0;
		Vector2D a,b;

		for(VertexIter v = vertices.begin(); v != vertices.end(); v++) {
			a = v->position;
			b = v->prev->position;
			surface += b.x * a.y - a.x * b.y;
		}

		surface *= 0.5;

		return std::abs(surface);
	}

	/*
	*
	* @return the unscaled area of the polygon
	*/
	double Polygon::originalArea()
	{
		return originalArea_;
	}

 /* ===========================================================================================
 * Computes the length of a polygon
 ==============================================================================================*/
	/*!
	* Calculates the length of polygon.
	*/
	double Polygon::length()
	{
		double l = 0.0;
		Vector2D a, b;

		for(VertexIter v = vertices.begin(); v != vertices.end(); v++) {
			a = v->position;
			b = v->next->position;
			l += (a-b).norm();
		}

		return l;
	}

 /* ===========================================================================================
 * Computes the length of boundary edges
 ==============================================================================================*/
	/*!
	* Calculates the length of boundary edges.
	*/
	void Polygon::boundaryLength0()
	{
		
		int n = vertices.size();
		bdLength0.clear();
		bdLength0.resize(n);
		Vector2D a, b;

		for(VertexIter v = vertices.begin(); v != vertices.end(); v++) {
			a = v->position;
			b = v->next->position;
			bdLength0[v->index] = (a-b).norm();
		}

	}

 /* ===========================================================================================
 * Computes the length of boundary edges
 ==============================================================================================*/
	/*!
	* Calculates the length of boundary edges.
	*/
	void Polygon::boundaryLength()
	{
		
		int n = vertices.size();
		bdLength.clear();
		bdLength.resize(n);
		Vector2D a, b;

		for(VertexIter v = vertices.begin(); v != vertices.end(); v++) {
			a = v->position;
			b = v->next->position;
			bdLength[v->index] = (a-b).norm();
		}

	}
 
 /* ===========================================================================================
 * Distorsion in edge lengths
 ==============================================================================================*/
	/*!
	* Calculates the distorsion in edge lengths.
	*/
	void Polygon::distorsion() {
		int n = vertices.size();
		double err=0;
		double val;
		for(int i = 0; i < n; i++) {
			val = ( (bdLength[i]-bdLength0[i]) /bdLength0[i]);
			err += 100*std::abs(val);
		}
		err = err/n;
		distorsion_ = err;
	}

 /* ===========================================================================================
 * Center and scale the polygon
 ==============================================================================================*/
	/*!
	* Center and rescale the polygon.
	* Automatically recenters to the center of mass of the vertices.
	* @param range
	* @param iscale
	*/
  	void Polygon::centerScale(double range, int iscale) {
		std::cout << "centerScale is called" << std::endl;
		Vector2D center, a;

		center[0] = 0; center[1] = 0;
		for(VertexIter v = vertices.begin(); v != vertices.end(); v++) {
			center += v->position;
		}
		center /= vertices.size();
		double r=0;

		for(VertexIter v = vertices.begin(); v != vertices.end(); v++) {
			v->position -= center;
	//		a = v->position;
	//		r = std::max(r, a.norm());
		}

		r = length();
		if(iscale !=0) {
			double scale = range/r;
			for(VertexIter v = vertices.begin(); v != vertices.end(); v++) {
				v->position *= scale;	
			}
		}
			void boundaryLength0();
			void boundaryLength();
			void distorsion();
	}

/* ===========================================================================================
 * Center and scale by area of the polygon
 ==============================================================================================*/
	/*!
	* Center and rescale the polygon.
	* Automatically recenters to the center of mass of the vertices.
	*/
  	void Polygon::centerScaleArea() {
		std::cout << "centerScaleArea is called" << std::endl;
		Vector2D center, a;

		center[0] = 0; center[1] = 0;
		for(VertexIter v = vertices.begin(); v != vertices.end(); v++) {
			center += v->position;
		}
		center /= vertices.size();
		double r=0;

		for(VertexIter v = vertices.begin(); v != vertices.end(); v++) {
			v->position -= center;
	//		a = v->position;
	//		r = std::max(r, a.norm());
		}

		double currentAreaBeforeScaling = area();
		std::cout << "area before scaling: " << currentAreaBeforeScaling << std::endl;
		double scale = sqrt(100/currentAreaBeforeScaling);
		std::cout << "scale is " << scale << std::endl;
		for(VertexIter v = vertices.begin(); v != vertices.end(); v++) {
			v->position *= scale;
		}

			void boundaryLength0();
			void boundaryLength();
			void distorsion();
		double currentArea = area();
		std::cout << "area after scaling: " << currentArea << std::endl;
		assert(std::abs(currentArea - 100.0) < 0.5 && "Area should be 100 after scaling");
	}


/*!
	* Rescale the polygon.
	* Automatically recenters to the center of mass of the vertices.
	*/
  	void Polygon::scaleArea() {
		std::cout << "scaleArea is called" << std::endl;
		double originalArea = std::abs(area());
		std::cerr << "originalArea is " << originalArea << std::endl;
		double scale = sqrt(100/originalArea);
		std::cerr << "scale is " << scale << std::endl;
		for(VertexIter v = vertices.begin(); v != vertices.end(); v++) {
			v->position *= scale;
		}

			void boundaryLength0();
			void boundaryLength();
			void distorsion();
		double currentArea = area();
		assert(std::abs(currentArea - 100.0) < 0.5 && "Area should be 100 after scaling");
	}

  /* ===========================================================================================
 * Shift the polygon
 ==============================================================================================*/

void Polygon::shift(Vector2D center, bool verbose = false) {
	for(VertexIter v = vertices.begin(); v != vertices.end(); v++) {
		v->position -= center;
		if (verbose) std::cout << "v was: (" << v->position.x << ", " << v->position.y << ")";
	};
	if (verbose) std::cout << "shifting has occured, center was selected as (" << center.x << ", " << center.y << ")." << std::endl;
	/*for(int v = 0; v < vertices.size(); v++) {
		if (verbose) std::cerr << "v was: (" << vertices[v].x() << ", " << vertices[v].y() << ")";
		vertices[v].shift(center.x, center.y);
		if (verbose) std::cerr << "v is now: (" << vertices[v].x() << ", " << vertices[v].y() << ")" << std::endl;

	};
	if (verbose) std::cerr << "shifting has occured, center was selected as (" << center.x << ", " << center.y << ")." << std::endl;*/
}
 /* ===========================================================================================
 * Constructor
 ==============================================================================================*/


 /*********************************************************************************
  Copy Constructor
 *********************************************************************************/
	/*!
	* Copy constructor.
	*/
  	Polygon::Polygon(Polygon& polygon) {

		// initialize vertices and copy from input polygon
		int nvertex = polygon.vertices.size();
		vertices.reserve(nvertex);
		std::vector<VertexIter> indexToVertex(nvertex);
		int i = 0;
		for (VertexIter v = polygon.vertices.begin(); v != polygon.vertices.end(); v++) {
			VertexIter vNew = vertices.insert(vertices.end(), Vertex());
			vNew->position = v->position;
			vNew->index = v->index;
			indexToVertex[i] = vNew;
			i++;
		}

		VertexIter v;
		int im1 = nvertex-1;
		for(int i = 0; i < nvertex-1; i++) {
			v = indexToVertex[i];
			v->prev = indexToVertex[im1];
			v->next = indexToVertex[i+1];
			im1 = i;
		}
		v = indexToVertex[nvertex-1];
		v->prev = indexToVertex[nvertex-2];
		v->next = indexToVertex[0];

		int nedges = polygon.edges.size();

		for(int i = 0; i < nedges; i++) {
			edges.push_back(polygon.edges[i]);
		}
	}

 /*********************************************************************************
  Assignment Operator
 *********************************************************************************/
	/*!
	* Assignment operator - properly handles prev/next iterator relinking
	*/
	Polygon& Polygon::operator=(const Polygon& other) {
		if (this == &other) {
			return *this;  // Handle self-assignment
		}

		// Clear existing data
		vertices.clear();
		edges.clear();
		bdLength0.clear();
		bdLength.clear();

		// Copy simple members
		distorsion_ = other.distorsion_;
		originalArea_ = other.originalArea_;

		// Copy vertices and relink prev/next iterators
		int nvertex = other.vertices.size();
		vertices.reserve(nvertex);
		std::vector<VertexIter> indexToVertex(nvertex);

		int i = 0;
		for (auto it = other.vertices.begin(); it != other.vertices.end(); ++it) {
			VertexIter vNew = vertices.insert(vertices.end(), Vertex());
			vNew->position = it->position;
			vNew->index = it->index;
			indexToVertex[i] = vNew;
			i++;
		}

		// Relink prev/next pointers
		int im1 = nvertex-1;
		for(int i = 0; i < nvertex-1; i++) {
			VertexIter v = indexToVertex[i];
			v->prev = indexToVertex[im1];
			v->next = indexToVertex[i+1];
			im1 = i;
		}
		VertexIter v = indexToVertex[nvertex-1];
		v->prev = indexToVertex[nvertex-2];
		v->next = indexToVertex[0];

		// Copy edges
		int nedges = other.edges.size();
		for(int i = 0; i < nedges; i++) {
			edges.push_back(other.edges[i]);
		}

		// Copy boundary lengths
		for(size_t i = 0; i < other.bdLength0.size(); i++) {
			bdLength0.push_back(other.bdLength0[i]);
		}
		for(size_t i = 0; i < other.bdLength.size(); i++) {
			bdLength.push_back(other.bdLength[i]);
		}

		return *this;
	}

	/*!
	* Computes the signed distance from a point to the polygon outline.
	* Positive distance indicates the point is outside the polygon,
	* negative distance indicates the point is inside the polygon.
	* @param point The point to compute the distance from.
	* @return The signed distance to the polygon outline.
	*/
	// double Polygon::signedDistanceToOutline(const Vector2D& point) const {
	// 	double minDistance = std::numeric_limits<double>::max();
	// 	bool isInside = false;
		
	// 	for (size_t i = 0; i < vertices.size(); ++i) {
	// 		const Vector2D& v1 = vertices[i].position;
	// 		const Vector2D& v2 = vertices[(i + 1) % vertices.size()].position;
			
	// 		// Compute the distance to the edge
	// 		Vector2D edge = v2 - v1;
	// 		Vector2D pointToV1 = point - v1;
	// 		double t = pointToV1.dot(edge) / edge.dot(edge);
	// 		t = std::max(0.0, std::min(1.0, t));
	// 		Vector2D closestPoint = v1 + edge * t;
	// 		double distance = (point - closestPoint).norm();
			
	// 		minDistance = std::min(minDistance, distance);
			
	// 		// Check if the point is inside the polygon (ray casting algorithm)
	// 		if (((v1.y > point.y) != (v2.y > point.y)) &&
	// 			(point.x < (v2.x - v1.x) * (point.y - v1.y) / (v2.y - v1.y) + v1.x)) {
	// 			isInside = !isInside;
	// 		}
	// 	}
		
	// 	return isInside ? -minDistance : minDistance;
	// }


	
	// void Polygon::labelPolygon(const double EPSILON = 1e-10) {
	// 	Vector2D centroid = calculateCentroid();
	// 	double minDist = signedDistanceToOutline(centroid);
	// 	Vector2D bestPoint = centroid;

	// 	std::vector<Vector2D> queue;
	// 	queue.push_back(centroid);

	// 	while (!queue.empty()) {
	// 		Vector2D point = queue.back();
	// 		queue.pop_back();

	// 		double dist = signedDistanceToOutline(point);
	// 		if (dist > minDist) {
	// 			minDist = dist;
	// 			bestPoint = point;
	// 		}

	// 		if (dist > minDist + EPSILON) {
	// 			double cellSize = dist / std::sqrt(2);
	// 			queue.push_back(Vector2D(point.x - cellSize, point.y - cellSize));
	// 			queue.push_back(Vector2D(point.x + cellSize, point.y - cellSize));
	// 			queue.push_back(Vector2D(point.x - cellSize, point.y + cellSize));
	// 			queue.push_back(Vector2D(point.x + cellSize, point.y + cellSize));
	// 		}
	// 	}

	// 	label = bestPoint;
	// }

 /* ===========================================================================================
 * Check if a point is inside the polygon using ray casting (Jordan curve theorem)
 ==============================================================================================*/
	/*!
	* Checks if a point is inside the polygon using the ray casting algorithm.
	* This implements the Jordan curve theorem: a point is inside a polygon if
	* a ray from the point crosses the polygon boundary an odd number of times.
	* @param point The point to check
	* @return true if the point is inside the polygon, false otherwise
	*/
	bool Polygon::isPointInside(const Vector2D& point) const {
		bool isInside = false;
		const double EPSILON = 1e-10;

		for (size_t i = 0; i < vertices.size(); ++i) {
			const Vector2D& v1 = vertices[i].position;
			const Vector2D& v2 = vertices[(i + 1) % vertices.size()].position;

			// Skip horizontal edges (parallel to the ray)
			if (std::abs(v1.y - v2.y) < EPSILON) {
				continue;
			}

			// Check if the ray intersects this edge
			// Only count vertices where the edge is going upward (v2.y > v1.y)
			// to avoid counting the same vertex twice
			if ((v1.y <= point.y && point.y < v2.y) ||
			    (v2.y <= point.y && point.y < v1.y)) {

				// Calculate the x-coordinate of the intersection point
				double x_intersect = v1.x + (point.y - v1.y) * (v2.x - v1.x) / (v2.y - v1.y);

				// Check for tangential crossing (ray passes exactly through a vertex)
				bool is_tangential = false;
				if (std::abs(point.y - v1.y) < EPSILON || std::abs(point.y - v2.y) < EPSILON) {
					// Ray passes through or very close to a vertex
					// Only count if both adjacent edges are on the same side of the horizontal line
					const Vector2D& v0 = vertices[(i + vertices.size() - 1) % vertices.size()].position;
					const Vector2D& v3 = vertices[(i + 2) % vertices.size()].position;

					if (std::abs(point.y - v1.y) < EPSILON) {
						// Near v1: check if v0 and v2 are on the same side
						is_tangential = ((v0.y > point.y) == (v2.y > point.y));
					} else {
						// Near v2: check if v1 and v3 are on the same side
						is_tangential = ((v1.y > point.y) == (v3.y > point.y));
					}
				}

				// Count the crossing if the ray intersects to the right and it's not tangential
				if (point.x < x_intersect && !is_tangential) {
					isInside = !isInside;
				}
			}
		}

		return isInside;
	}
} //end namespace correa
#endif
