/*!
* @file	Curvature.h
* @brief  	This file looks at option to compare two curves based on the curvatures of their vertices. It is based on the paper: K. Crane, U. Pinkall, and P. Schroeder, "Robust fairing via conformal curvature flow", ACM Transactions on Graphics, Vol. 32, No. 4 (2013), and on the code SpinXform_OpenGL from Keenan Crane (available at: https://www.cs.cmu.edu/~kmcrane/index.html#code)
* @author Patrice Koehl
* @date October 2021
* @version 1
* @copyright BSD 3-Clause License.
*/


#ifndef _CURVATURE_H
#define _CURVATURE_H

  /* ===== INCLUDES AND TYPEDEFS =================================================================
   *
   =============================================================================================== */

#include <vector>
#include <cmath>
#include <OT1.h>

namespace correa {
/*================================================================================================
  Prototypes for BLAS and LAPACK
================================================================================================== */

  /* ===== The Curvature class    =================================================================
   *
   =============================================================================================== */

/*!
* Curvature class.
* Used to calculate various curvature related properties.
*/
class Curvature {

	public:
		// Computes Willmore energy for a polygon
		double Willmore(Polygon& polygon);

		// Computes Wasserstein distance between two curves, using the curvatures at
		// each vertex to compute a cost matrix
		double curvOT(Polygon& curve1, Polygon& curve2);

		// Build histogram of "curvature" (only external angle) for a polygon
		void curvHist(Polygon& polygon, int *nhist1, double *x1, double *hist1);

		// Computes Wasserstein distance between two histograms of curvatures
	
		double histOT(Polygon& curve1, Polygon& curve2);

	private:
		OT1 ot1;

	protected:

};


  /* =============================================================================================
   Computes Willmore energy for a polygon
   =============================================================================================== */
	/*!
	* Calculate the Willmore energy of a polygon.
	* @param polygon	Polygon object
	* @return E 		Willmore energy
	*/
  	double Curvature::Willmore(Polygon& polygon) {
		Vector2D a, b, c;
		double l, l_prev, curv;
		double E = 0;
		for (VertexIter v = polygon.vertices.begin(); v != polygon.vertices.end(); v++) {
			a = v->position;
			b = v->prev->position;
			c = v->next->position;
			l_prev = (a-b).norm();
			l      = (a-c).norm();
			curv = v->exteriorAngle()/(0.5*(l+l_prev));
			E += curv*curv*(0.5*(l+l_prev));
		}
		E /= polygon.length();

		return E;

	}

  /* =============================================================================================
   Computes Wasserstein distance between two curves, using the curvatures at
  each vertex to compute a cost matrix
   =============================================================================================== */
	/*!
	* Calculate the Wasserstein distance between two curves, using the curvatures at each vertex to compute a cost matrix
	* @param curve1 Polygon representing the first curve
	* @param curve2 Polygon representing the second curve
	* @return d_OT Wasserstein distance between two curves
	*/
   double Curvature::curvOT(Polygon& curve1, Polygon& curve2) {
		int n1 = curve1.vertices.size();
		int n2 = curve2.vertices.size();

		double *curv1 = new double[n1];
		double *curv2 = new double[n2];

		// Compute curvatures of all vertices in polygon 1
		Vector2D a, b, c;
		double l, l_prev ;
		for (VertexIter v = curve1.vertices.begin(); v != curve1.vertices.end(); v++) {
			a = v->position;
			b = v->prev->position;
			c = v->next->position;
			l_prev = (a-b).norm();
			l      = (a-c).norm();
			curv1[v->index] = v->exteriorAngle()/(0.5*(l+l_prev));
		}

		// Compute curvatures of all vertices in polygon 2
		for (VertexIter v = curve2.vertices.begin(); v != curve2.vertices.end(); v++) {
			a = v->position;
			b = v->prev->position;
			c = v->next->position;
			l_prev = (a-b).norm();
			l      = (a-c).norm();
			curv2[v->index] = v->exteriorAngle()/(0.5*(l+l_prev));
	}

	// Assign sizes of all arrays needed for OT1

	int nmax = std::max(n1,n2);
	int wsize = 2*nmax*nmax + 10*nmax;
	double *Cost = new double[n1*n2];
	double *G    = new double[n1*n2];
	double *Work = new double[wsize];
	double *m1   = new double[n1];
	double *m2   = new double[n2];
	double *ones = new double[nmax];

	// Initialize arrays
        for(int i = 0; i < n1; i++) m1[i] = 1./n1;
        for(int i = 0; i < n2; i++) m2[i] = 1./n2;
        for(int i = 0; i < nmax; i++) ones[i] = 1.;

	for(int i = 0; i < n1; i++) {
		for(int j = 0; j < n2; j++) {
			Cost[i+j*n1] = std::abs(curv1[i]-curv2[j]);
		}
	}

	// Perform OT
	double beta1 = 1.e3;
	double beta2 = 1.e10;
	int iprint = 0;
	int nthreads = sysconf( _SC_NPROCESSORS_ONLN );
	nthreads = std::min( nthreads, NUM_THREADS);
	if(nthreads==0) nthreads = 1;

	double d_OT = ot1.ot1(n1, m1, n2, m2, Cost, G, ones, Work,
                        beta1, beta2, iprint, nthreads);

	// Clean up
	delete [] curv1; delete [] curv2;
	delete [] Cost; delete [] G; delete [] m1; delete [] m2;
	delete [] Work; delete [] ones;

	return d_OT;

   }

  /* =============================================================================================
     Computes histogram of exterior angles for all vertices in a polygon
   =============================================================================================== */
	/*!
	* Build histogram of "curvature" for a polygon
	* describe what it does?
	* @param curve1 Polygon representing the first curve
	* @param curve2 Polygon representing the second curve
	* @param nhist
	* @param x
	* @param hist
	*/
   void Curvature::curvHist(Polygon& polygon, int *nhist, double *x, double *hist) {

		int nvert = polygon.vertices.size();

		// angles are in [0, 180], which we divide in 180 intervals of 1 degree...
		int nval = 180.0;
		*nhist = 180;
		double step = 1.0;
		for(int i = 0; i < nval; i++) x[i] = step/2 + i*step;

		memset(hist, 0, nval*sizeof(double));
		double ang;
		int n;
		for (VertexIter v = polygon.vertices.begin(); v != polygon.vertices.end(); v++) {
			ang = v->exteriorAngle() /M_PI;
			n = std::floor(ang);
			std::cout << "ang = " << ang*180 << " n = " << n << std::endl;
			hist[n] +=1;
		}
		exit(1);

		for(int i = 0; i < nval; i++) hist[i] /= nvert;

	}

  /* =============================================================================================
   Computes Wasserstein distance between the histograms of the exterior angles of two curves
   =============================================================================================== */
	/*! 
	* Calculate the Wasserstein distance between two historgrams of curvatures
	* @param curve1 Polygon representing the first curve
	* @param curve2 Polygon representing the second curve
	*/
   double Curvature::histOT(Polygon& curve1, Polygon& curve2) {
		int n1, n2;
		int nval = 180;
		double *x1 = new double[nval]; double *x2 = new double[nval];
		double *hist1 = new double[nval]; double *hist2 = new double[nval];

		curvHist(curve1, &n1, x1, hist1);
		curvHist(curve2, &n2, x2, hist2);

		// Assign sizes of all arrays needed for OT1

		int nmax = std::max(n1, n2);
		int wsize = 2*nmax*nmax + 10*nmax;
		double *Cost = new double[n1*n2];
		double *G    = new double[n1*n2];
		double *Work = new double[wsize];
		double *m1   = new double[n1];
		double *m2   = new double[n2];
		double *ones = new double[nmax];

		// Initialize arrays
			for(int i = 0; i < n1; i++) m1[i] = hist1[i];
			for(int i = 0; i < n2; i++) m2[i] = hist2[i];
			for(int i = 0; i < nmax; i++) ones[i] = 1.;

		double s1 = 0, s2=0;
			for(int i = 0; i < n1; i++) s1 +=m1[i];
			for(int i = 0; i < n2; i++) s2 +=m2[i];
		std::cout << "s1 = " << s1 << std::endl;
		std::cout << "s2 = " << s2 << std::endl;

		double Cmax = 0.0;
		for(int i = 0; i < n1; i++) {
			for(int j = 0; j < n2; j++) {
				Cost[i+j*n1] = std::abs(x1[i]-x2[j]);
				Cmax = std::max(Cmax, Cost[i+j*n1]);
			}
		}
		for(int i = 0; i < n1*n2; i++) Cost[i] = Cost[i]/Cmax;

		// Perform OT
		double beta1 = 0.01;
		double beta2 = 1.e6;
		int iprint = 1;
		int nthreads = sysconf( _SC_NPROCESSORS_ONLN );
		nthreads = std::min( nthreads, NUM_THREADS);
		if(nthreads==0) nthreads = 1;

		double d_OT = ot1.ot1(n1, m1, n2, m2, Cost, G, ones, Work,
							beta1, beta2, iprint, nthreads);

		// Clean up
		delete [] Cost; delete [] G; delete [] m1; delete [] m2;
		delete [] Work; delete [] ones;

		return d_OT;

	}
}//end namespace correa
#endif

