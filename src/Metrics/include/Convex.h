/*
 	Convex.h

 	Authors: 	Patrice Koehl, Department of Computer Science, University of California, Davis
				Yossi Bokor Bleile, Department of Mathematical Sciences, University of Aalborg, Aalborg
 	Date: April 2023
	Version: 1
*/


#ifndef _CONVEX_H_
#define _CONVEX_H_

  #include <cmath>
  #include <algorithm>
  #include <functional>
namespace correa {
/* ===============================================================================================
   The ConvexHull class
   =============================================================================================== */

  class Convex {

  public:

	// Compute convex hull of a set of points
	void ConvexPoint(std::vector<Vector2D>& Points, std::vector<Vector2D>& CH);

	// Compute convex hull of a polygon
	void ConvexPoly(std::vector<Vector2D>& Points, std::vector<Vector2D>& CH);

  private:

	// Check orientation of 3 points
   	bool ccw(Vector2D& a, Vector2D& b, Vector2D& c);

 };

/* ==============================================================================================
   Checks if three points make a counter-clockwise turn
   =============================================================================================== */

   bool Convex::ccw(Vector2D& a, Vector2D& b, Vector2D& c)
   {
	double val1 = (b[0] - a[0]) * (c[1] - a[1]);
	double val2 = (b[1] - a[1]) * (c[0] - a[0]);

	return val1 > val2;
   }


/* ==============================================================================================
   Generate convex hull of a set of points using Graham scans
   =============================================================================================== */

   void Convex::ConvexPoint(std::vector<Vector2D>& Points, std::vector<Vector2D>& CH)
   {
	int n = Points.size();
	if (n == 0) return;

	std::vector<Vector2D> copy;

	for(int i = 0; i < n; i++) copy.push_back(Points[i]);

	std::sort(copy.begin(), copy.end(), [](Vector2D& a, Vector2D& b){
		if (a[0] < b[0] ) return true;
		return false;
	});
 
	CH.clear();
 
	// lower hull
	for (auto& pt : copy) {

		while (CH.size() >= 2 && !ccw(CH.at(CH.size() - 2), CH.at(CH.size() - 1), pt)) {
			CH.pop_back();
        	}
		CH.push_back(pt);
	}
 
	// upper hull
	auto t = CH.size() + 1;
	for (auto it = copy.crbegin(); it != copy.crend(); it = std::next(it)) {
		auto pt = *it;
		while (CH.size() >= t && !ccw(CH.at(CH.size() - 2), CH.at(CH.size() - 1), pt)) {
			CH.pop_back();
		}
		CH.push_back(pt);
	}
 
	CH.pop_back();
   }

/* ==============================================================================================
   Generate convex hull of a polyline using Melkman's algorithm
   =============================================================================================== */

   void Convex::ConvexPoly(std::vector<Vector2D>& Points, std::vector<Vector2D>& CH)
   {
	int n = Points.size();
	if (n == 0) return;

	std::vector<Vector2D> Dqueue;
	Dqueue.reserve(2*n+1);

	// initialize Dqueue from bottom to top so that the
	// 1st three vertices of Pointd[] are a ccw triangle

	int bot = n-2, top = bot+3;    		// initial bottom and top deque indices
	Dqueue[bot] = Dqueue[top] = Points[2];	// 3rd vertex is at both bot and top

	if (ccw(Points[0], Points[1], Points[2])) {
		Dqueue[bot+1] = Points[0];
		Dqueue[bot+2] = Points[1];           // ccw vertices are: 2,0,1,2
	}
	else {
		Dqueue[bot+1] = Points[1];
		Dqueue[bot+2] = Points[0];           // ccw vertices are: 2,1,0,2
	}

	// compute the hull on the deque Dqueue

	for (int i=3; i < n; i++) {   // process the rest of vertices
		// test if next vertex is inside the deque hull
		if ((ccw(Dqueue[bot], Dqueue[bot+1], Points[i])) &&
		   (ccw(Dqueue[top-1], Dqueue[top], Points[i])) )
			continue;         // skip an interior vertex

		// incrementally add an exterior vertex to the deque hull
		// get the rightmost tangent at the deque bot

		while (!ccw(Dqueue[bot], Dqueue[bot+1], Points[i])) {
			++bot;                  	// remove bot of deque
		}
		Dqueue[--bot] = Points[i];	        // insert P[i] at bot of deque

		// get the leftmost tangent at the deque top

		while (!ccw(Dqueue[top-1], Dqueue[top], Points[i]))
			--top;                 // pop top of deque
		Dqueue[++top] = Points[i];          // push P[i] onto top of deque
	}

	// transcribe deque Dqueue to the output hull CH

	for (int h=0; h <= (top-bot); h++)
        	CH.push_back(Dqueue[bot + h]);

	Dqueue.clear();
	return;

   }
} //end namespace correa
#endif


