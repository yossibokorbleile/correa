/* ====FRECHET_H =================================================================================
 *
 * Author: Patrice Koehl (in collaboration with Yossi Bokor, U. of Sydney), October 2021 
 * Department of Computer Science
 * University of California, Davis
 *
 * This file implements the discrete Frechet distance for comparing two polygonal curves
 *
 * Implementation based on:
 *      Eiter, Thomas; Mannila, Heikki (1994), Computing discrete Fréchet distance
 *      Tech. Report CD-TR 94/64, Christian Doppler Laboratory for Expert Systems, TU Vienna, Austria.
 *
 * LGPL licensing
 =============================================================================================== */

#ifndef _FRECHET_H_
#define _FRECHET_H_

  #include <cmath>
  #include <float.h>
  #include <algorithm>
  #include <functional>

/* ===============================================================================================
   The Frechet class
   =============================================================================================== */

  class Frechet {

  public:

	// Compute discrete Frechet distance between 2 curves
	double dFD(Polygon& curve1, Polygon& curve2);

  private:

	// Coupling search function
	double coupling(Polygon& curve1, Polygon& curve2, int i, int j, double *CA);

  };

/* ===============================================================================================
   Frechet

   Notes:
   =====
   curve1 and curve2 are two sets of points that define 2 polygonal curves
   The points along the curves are expected to be in the order as they appear.
  
   Returned is the discrete Frechet distance, i.e. the coupling
   measure, which is zero when curve1 equals curve2 and grows positively as the curves
   become more dissimilar.
  
   The L2 norm is used to compute the distances between points in curve1 and curve2.
  
   Explanation:
   ===========
   The Frechet distance is a measure of similarity between to curves.
   It is defined as the minimum cord-length sufficient to join a point
   traveling forward along curve1 and one traveling forward along curve2, although the
   rate of travel for either point may not necessarily be uniform.
  
   The Frechet distance, is not in general computable for any given
   continuous curves . However, the discrete Frechet Distance, also called
   the coupling measure, is a metric that acts on the endpoints of
   curves represented as polygonal chains. The magnitude of the coupling
   measure is bounded by the continuous Frechet distance plus the length of the 
   longest segment in either curve1 or curve2.
  
   This function implements the algorithm to calculate discrete Frechet
   distance outlined table 1 in:
   T. Eiter and H. Mannila. Computing discrete Frechet distance. Technical
   Report 94/64, Christian Doppler Laboratory, Vienna University of
   Technology, 1994.

   It loosely follows a matlab implementation:
   Zachary Danziger (2021). Discrete Frechet Distance 
   (https://www.mathworks.com/matlabcentral/fileexchange/31922-discrete-frechet-distance), 
   MATLAB Central File Exchange. Retrieved October 15, 2021.
   
   =============================================================================================== */

  double Frechet::dFD(Polygon& curve1, Polygon& curve2)
  {

	// Number of points in the two curves

	int n1 = curve1.vertices.size();
	int n2 = curve2.vertices.size();
	int n1n2 = n1*n2;

	// Initialize a matrix CA of size n1 x n2, with all -1

	double *CA = new double[n1*n2];
	for(int i = 0; i < n1n2; i++) CA[i] = -1;
	
	// The discrete Frechet distance is the coupling for the last points in curve1 and curve2
	double d = coupling(curve1, curve2, n1-1, n2-1, CA);	

	return d;
  }

/* ===============================================================================================
   Coupling search function
   Note: this function is recursive
   =============================================================================================== */

  double Frechet::coupling(Polygon& curve1, Polygon& curve2, int i, int j, double *CA)
  {
	int n1 = curve1.vertices.size();
	double val;

	if(CA[i+j*n1] > -1 ) return CA[i+j*n1];
	if (i==0 && j == 0) {
		CA[i+j*n1] = (curve1.vertices[0].position -curve2.vertices[0].position).norm();
	} else if (i>0 && j == 0) {
		CA[i+j*n1] = std::max( coupling(curve1, curve2, i-1, 0, CA), (curve1.vertices[i].position -curve2.vertices[0].position).norm());
	} else if (i==0 && j > 0) {
		CA[i+j*n1] = std::max( coupling(curve1, curve2, 0, j-1, CA), (curve1.vertices[0].position -curve2.vertices[j].position).norm());
	} else if (i>0 && j > 0) {
		val = std::min( coupling(curve1, curve2, i-1, j, CA), coupling(curve1, curve2, i-1, j-1, CA));
		val = std::min( val, coupling(curve1, curve2, i, j-1, CA));
		CA[i+j*n1] = std::max(val, (curve1.vertices[i].position -curve2.vertices[j].position).norm());
	} else {
		CA[i+j*n1] = DBL_MAX;
	}
	return CA[i+j*n1];
  }

#endif

