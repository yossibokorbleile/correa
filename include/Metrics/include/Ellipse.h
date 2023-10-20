/*!
* @file Ellipse.h
* @brief calculate various ellipse for a polygon
* @author Patrice Koehl
* @author Yossi Bokor Bleile
* @date April 2023
* @version 1
* @copyright BSD 3-Clause License.
*/

#ifndef _ELLIPSE_H_
#define _ELLIPSE_H_

#include <cmath>
#include <algorithm>
#include <functional>
#include "Convex.h"
namespace correa {
Convex convexhull;

/* ===============================================================================================
   prototypes for BLAS and LAPACK functions 
   =============================================================================================== */

  extern "C" {

	// BLAS1: copy one vector into another
	void dcopy_(int * n, double * X, int *incx, double * Y, int *incy);

	// BLAS1: dot product of two vectors    
	double ddot_(int * n, double * u, int * incu, double * v, int *incv);

	// BLAS1: norm of a vector
	double dnrm2_(int * n, double * X, int *incx);

	// BLAS1: Y <- alpha*X + Y
	void daxpy_( int * n, double *alpha, double *X, int *incx, double *Y, int *incy);

	// BLAS1: X <- alpha*X
	void dscal_(int * n, double * alpha, double * X, int *incx);

	// BLAS2: perform Y := alpha*op( A )* B  + beta*Y
	void dgemv_(char * trans, int * m, int * n, double * alpha, double *A,
	int *lda, double * X, int * incx, double *beta, double * Y, int * incy);

	// BLAS3: perform C := alpha*op( A )* op(B)  + beta*C
	void dgemm_(char * transa, char * transb, int * m, int * n, int * k,
		double * alpha, double * A, int * lda,
		double * B, int * ldb, double * beta, double * C, int * ldc);


	// LAPACK: diagonalize a symmetric matrix
	void dgeev_(char * JOBVL, char * JOBVR, int *N, double *A, int *LDA, double *WR,
	double *WI, double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO);

	// LAPACK: diagonalize a general matrix
	void dsyevd_(char * JOBZ, char * UPLO, int *N, double *A, int *LDA, double *W,
	double *WORK, int *LWORK, int *IWORK, int *LIWORK, int *INFO);

	// LAPACK: LU algorithm
	void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);

	// LAPACK: solve A*X=B or A^T * X = B based on LU decomposition
	void dgetrs_(char *Trans, int *N, int *Nrhs, double *A, int *LDA, 
	int *IPIV, double *B, int *LDB, int *INFO);

	void dgetri_(int *N, double *A, int *LDA, int *IPIV, double *WORK,
	int *LWORK, int *INFO);

}


/* ===============================================================================================
   The Ellipse class
   =============================================================================================== */

  class Ellipse {

  public:

	// Compute major/minor axes of the largest enclosed ellipse in a polygon
	void EllipseMax(Polygon& curve, double *a, double *b);

	// Compute major/minor axes of the smallest enclosing ellipse of a polygon
	void EllipseMin(Polygon& curve, double *a, double *b);

	// Compute major/minor axes of the least square ellipse of a polygon
	void EllipseLSQ(Polygon& curve, double *a, double *b);

	// Compute difference of aspect ratio between 2 curves based on the largest enclosed ellipse
	double dEllipseMax(Polygon& curve1, Polygon& curve2, double *a1, double *b1, double *a2, double *b2);

	// Compute difference of aspect ratio between 2 curves based on the smallest enclosing ellipse
	double dEllipseMin(Polygon& curve1, Polygon& curve2, double *a1, double *b1, double *a2, double *b2);

	// Compute difference of aspect ratio between 2 curves based on the least square fitting ellipse
	double dEllipseLSQ(Polygon& curve1, Polygon& curve2, double *a1, double *b1, double *a2, double *b2);

  private:

	void MaxVolEllipse(double *P, int N, int d, double *A, double *c, double tol);

	void MinVolEllipse(double *P, int N, int d, double *A, double *c, double tol);

	void lsqEllipse(double *P, int N, int d, double *a, double *b);

	void Ellipse2D(double *A, double *c, double *a, double *b, double *ellipse);

	void matInv2(double *A, double *Ainv);


  };

/* ===============================================================================================
   MinVolEllipse

   Finds the minimum volume enclosing ellipsoid of a set of data
   points stored in a matrix P. The following optimization problem is solved: 
  
   minimize       log(det(A))
   subject to     (P_i - c)' * A * (P_i - c) <= 1
                  
   in variables A and c, where P_i is the i-th column of the matrix P. 
   The solver is based on Khachiyan Algorithm, and the final solution 
   is different from the optimal value by the pre-specified amount of 'tol'.
  
   inputs:
  ---------
   P   : (d x N) dimensional matrix containing N points in R^d.
   tol : error in the solution with respect to the optimal value.
  
   outputs:
  ---------
   A            : matrix representing ellipse (in fact, inv(A))
   c            : center of ellipse

   Reference:
   ----------
   Code based on the Matlab package MinVolEllipse:

   	Nima Moshtagh Minimum Volume Enclosing Ellipsoid 
   	https://www.mathworks.com/matlabcentral/fileexchange/9542-minimum-volume-enclosing-ellipsoid, 
	MATLAB Central File Exchange. 

   and referencing paper:

	Nima Moshtagh, "Minimum volume enclosing ellipsoid"
	http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.116.7691

   =============================================================================================== */

  void Ellipse::MinVolEllipse(double *Points, int Np, int d, double *A, double *c, double tol)
  {
	// Build convex hull of the polygon defined by P

	std::vector<Vector2D> CH;
	std::vector<Vector2D> Poly;

	for(int i = 0; i < Np; i++) {
		Vector2D v(Points[2*i], Points[2*i+1]);
		Poly.push_back(v);
	}
	convexhull.ConvexPoly(Poly, CH);

	int N = CH.size()-1;
	double *P = new double[2*N];
	for(int i = 0; i < N; i++) {
		P[2*i + 0] = CH[i][0];
		P[2*i + 1] = CH[i][1];
	}

	double *Q = new double[(d+1)*N];
	memset(Q, 0.0, (d+1)*N*sizeof(double));

	for (int j = 0; j < N; j++) {
		for(int dim = 0; dim < d; dim++) {
			Q[dim+(d+1)*j] = P[dim+d*j];
		}
		Q[d+(d+1)*j] = 1.0;
	}

        /* =========================================================================================
	   Initialisation
        ========================================================================================= */

	double err = 2*tol;

	double *u = new double[N];
	double *new_u = new double[N];
	for(int i = 0; i < N; i++) u[i] = 1./N;

	double *X = new double[(d+1)*(d+1)];
	int *IPIV = new int[d+1];
	int lwork = 100*(d+1);
	double *WORK = new double[lwork];
	double *v = new double[d+1];
	double *w = new double[d+1];
	double *M = new double[N];

	int dp1 = d+1;
	int dd = d*d;
	int info;
	int inc = 1; int one = 1;
	char NoTrans = 'N'; char Trans = 'T';
	double alpha, beta;
	double dMax, step, sum;
	int idx;

        /* =========================================================================================
	Khachiyan Algorithm
        ========================================================================================= */

	while (err > tol) {

		// X = Q * diag(u) * Q'; is a (d+1)x(d+1) matrix

		for(int i = 0; i < d+1; i++) {
			for(int j = 0; j < d+1; j++) {
				sum = 0.0;
				for(int k = 0; k < N; k++) {
					sum += Q[i+k*(d+1)]*u[k]*Q[j+k*(d+1)];
				}
				X[i+(d+1)*j] = sum;
			}
		}

		// M = diag(Q' * inv(X) * Q);  % M the diagonal vector of an NxN matrix

		dgetrf_(&dp1, &dp1, X, &dp1, IPIV, &info);
		dgetri_(&dp1, X, &dp1, IPIV, WORK, &lwork, &info);

		alpha = 1.0; beta = 0;
		for(int i = 0; i < N; i++) {
			dcopy_(&dp1, &Q[i*dp1], &inc, v, &inc);
			dgemv_(&NoTrans, &dp1, &dp1, &alpha, X, &dp1, v, &inc, &beta, w, &inc);
			M[i] = ddot_(&dp1, v, &inc, w, &inc);
		}

		dMax = M[0]; idx = 0;
		for(int i = 1; i < N; i++) {
			if(M[i] > dMax) {
				dMax = M[i];
				idx = i;
			}
		}

		// Compute U_new
		
		step = (dMax -d -1)/(dp1*(dMax-1));
		alpha = 1.0 - step;
		memset(new_u, 0.0, N*sizeof(double));
		daxpy_(&N, &alpha, u, &inc, new_u, &inc);
		new_u[idx] += step;

		// Check error

		alpha = -1.0;
		daxpy_(&N, &alpha, new_u, &inc, u, &inc);
		err = dnrm2_(&N, u, &inc);

		// Update U
		
		dcopy_(&N, new_u, &inc, u, &inc);
	}

        /* =========================================================================================
	Find the ellipse equation
		(x-c)' * A * (x-c) = 1
	i.e. compute a dxd matrix 'A' and a d dimensional vector 'c' as the center
	of the ellipse:

	A = (1/d) * inv(P * diag(u) * P' - (P * u)*(P*u)' )
	c = P*u

	In fact, we only need inv(A) = d * ( P * diag(u) * P' - (P * u)*(P*u)' )
        ========================================================================================= */

	// Compute c = P*U
	alpha = 1.0; beta = 0.0;
	dgemv_(&NoTrans, &d, &N, &alpha, P, &d, u, &inc, &beta, c, &inc);

	// Compute inv(A) = d * ( P * diag(u) * P' - (P * u)*(P*u)' )
	for(int i = 0; i < d; i++) {
		for(int j = 0; j < d; j++) {
			sum = 0.0;
			for(int k = 0; k < N; k++) {
				sum += P[i+k*d]*u[k]*P[j+k*d];
			}
			A[i+j*d] = sum;
		}
	}

	alpha = -1.0; beta = 1.0;
	dgemm_(&NoTrans, &Trans, &d, &d, &one, &alpha, c, &d, c, &d, &beta, A, &d);

	alpha = d; 
	dscal_(&dd, &alpha, A, &inc);

	// Clean up
	delete [] P; delete [] Q; delete [] v; delete [] w; delete [] WORK; delete [] IPIV;

  }

/* ===============================================================================================
   Ellipse2D

   Given the matrix representation A of 2D ellipse and its center c
		(x-c)' * A * (x-c) = 1
   compute its major a and minor axis b, as well as a discrete representation
   of the ellipse over 100 points

   inputs:
  ---------
   A   : (d x d) dimensional matrix (in fact, A^{-1}
   c   : center of ellipse
   
   outputs:
  ---------
   a   : size along major axis of the ellipse
   b   : size along minor axis of the ellipse
   ellipse: discrete representation of the ellipse over 100 points

   =============================================================================================== */

  void Ellipse::Ellipse2D(double *A, double *c, double *a, double *b, double *ellipse)
  {

	// Diagonalize matrix inv(A) 
	int d = 2;
	double *w = new double[d];

	int lwork = 2*d*d+6*d+1; double *WORK = new double[lwork];
	int liwork = 5*d+3; int *IWORK = new int[liwork];

	char JOBZ='V'; char UPLO='U';
	int info;
	dsyevd_(&JOBZ, &UPLO, &d, A, &d, w, WORK, &lwork, IWORK, &liwork, &info);

	// Major and Minor axis are the square roots of the eigenvalues
	*a = std::sqrt(w[1]);
	*b = std::sqrt(w[0]);

	// Orientation of the ellipse defined from the eigenvector associated with major axis
	int inc = 1;
	dcopy_(&d, &A[d], &inc, w, &inc);

	double ang = std::atan2(w[1],w[0]);

	// Compute discrete represention of ellipse over 100 points
	int ntime = 100;
	double *t = new double[ntime];
	double step = 2*M_PI/(ntime-1);
	for(int i = 0; i < ntime; i++) t[i] = -M_PI + i*step;

	double x1, x2;
	for(int i = 0; i < ntime; i++) {
		x1 = (*a)*std::cos(t[i]);
		x2 = (*b)*std::sin(t[i]);
		ellipse[2*i] = std::cos(ang)*x1 - std::sin(ang)*x2 + c[0];
		ellipse[2*i+1] = std::sin(ang)*x1 + std::cos(ang)*x2 + c[1];
	}

	// Clean up
	delete [] WORK; delete [] IWORK; delete [] t;

  }

/* ===============================================================================================
   MaxVolEllipse

   Finds the maximum volume inclosing ellipsoid of a convex polygon.
   As the input is a general polygon, first computes the convex hull of that polygon
   The following optimization problem is solved: 

   Find the maximum volume ellipsoid
       {v:  v = x + Es, ||s|| <= 1}
   inscribing a full-dimensional polytope
          {v:  Av <= b}

   inputs:
  ---------
   P   : (d x N) dimensional matrix containing N points in R^d.
   tol : error in the solution with respect to the optimal value.
  
   outputs:
  ---------
   E            : matrix representing ellipse (in fact, E*E')
   c            : center of ellipse

   Reference:
   ----------
   Implementation based on:

     Yin Zhang and Liyan Gao. On Numerical Solution of the Maximum Volume Ellipsoid Problem. 
	SIAM Journal on Optimization, Vol.14, No.1, pp. 53-76, 2003.

    and corresponding MATLAB code available at:
	https://www.caam.rice.edu/~zhang/mve/index.html

   =============================================================================================== */

  void Ellipse::MaxVolEllipse(double *P, int N, int d, double *E2, double *x, double tol)
  {

	// parameters
	int maxiter = 50;
	double minmu = 1.e-8;
	double tau0 = 0.75;
	int iprint = 0;

	// Build convex hull of the polygon defined by P

	std::vector<Vector2D> CH;
	std::vector<Vector2D> Poly;

	for(int i = 0; i < N; i++) {
		Vector2D v(P[2*i], P[2*i+1]);
		Poly.push_back(v);
	}
	convexhull.ConvexPoly(Poly, CH);
	int nCH = CH.size();

	// Center of gravity of convex hull

	double x0[2];
	x0[0] = 0; x0[1] = 0;
	for(int i = 0; i < nCH-1; i++) {
		x0[0] += CH[i][0];
		x0[1] += CH[i][1];
	}
	x0[0] /= (nCH-1); x0[1] /= (nCH-1);

	// Generate affine constraints A and b

	int n = nCH-1; int nn=n*n;
	double *A = new double[n*d];
	double *b = new double[n];
	double *bmAx0 = new double[n];
	double r, s, t, u, val;
	for(int i = 0; i < n; i++) {
		r = CH[i+1][0]-CH[i][0]; s = CH[i+1][1]-CH[i][1];
		val = std::sqrt(r*r+s*s);
		t = r; r = -s/val; s = t/val; 
		t = CH[i+1][0]+CH[i][0];
		u = CH[i+1][1]+CH[i][1];
		A[i] = r; A[i+n] = s;
		b[i] = 0.5*(r*t+s*u);
		val = r*x0[0]+s*x0[1] -b[i];
		if(val > 0) {
			A[i] = -A[i]; A[i+n]=-A[i+n]; b[i] = -b[i];
		}
	}

	int inc = 1;
	double alpha = -1; double beta = 1;
	char NoTrans = 'N'; char Trans = 'T';

	double bnrm = dnrm2_(&n, b, &inc);
	dcopy_(&n, b, &inc, bmAx0, &inc);
	dgemv_(&NoTrans, &n, &d, &alpha, A, &n, x0, &inc, &beta, bmAx0, &inc);

	// Check that center of gravity, Xmean, is inside the convex polygon
	// Should never be a problem...

	for(int i = 0; i < n; i++) {
		if(bmAx0[i] <= 0) {
			std::cout << "Problem in MaxVolEllipse!" << std::endl;
			std::cout << "Center of gravity not in convex hull...." << std::endl;
			exit(1);
		}
	}

	for(int i = 0; i < n; i++) {
		for(int k = 0; k < d; k++) {
			A[i+k*n] /= bmAx0[i];
		}
		b[i] = 1.0;
	}

	double *y = new double[n]; double *z = new double[n];
	double *dx = new double[d]; double *dy = new double[n]; double *dz = new double[n];
	double *yz = new double[n]; double *yh = new double[n]; double *bmAx = new double[n];

	memset(x, 0, d*sizeof(double));
	for (int i = 0; i < n; i++) {
		y[i] = 1.0;
		bmAx[i] = 1.0;
	}

	double res = 1;
	double sum, gap, rmu;

	double *Ei = new double[d*d];
	double *W = new double[n*d];
	double *Q = new double[n*n]; double *G = new double[n*n];
	double *h = new double[n]; double *Adx = new double[n];
	double *R1 = new double[d]; double *R2 = new double[n]; double *R3 = new double[n];
	int *IPIV = new int[n];

	int info;
	double r1, r2, r3, objval;
	double ax, ay, az, tau, astep;

	if(iprint==1) {
		std::cout << "=========================================================" << std::endl;
		std::cout << "  Residuals:     Primal    Dual    Duality     log(detE) " << std::endl;
		std::cout << "=========================================================" << std::endl;
	}

	// Start iterations
	for(int iter = 0; iter < maxiter; iter++) {

		if(iter > 0) {
			alpha = -astep;
			daxpy_(&n, &alpha, Adx, &inc, bmAx, &inc);
		}

		for(int i = 0; i < d; i++) {
			for(int j = 0; j < d; j++) {
				sum = 0.;
				for(int k = 0; k < n; k++) {
					sum+= A[k+i*n]*y[k]*A[k+j*n];
				}
				Ei[i+d*j] = sum;
			}
		}
		matInv2(Ei, E2);
		
		alpha = 1.0; beta = 0.0;
		dgemm_(&NoTrans, &NoTrans, &n, &d, &d, &alpha, A, &n, E2, &d, &beta, W, &n);
		dgemm_(&NoTrans, &Trans, &n, &n, &d, &alpha, W, &n, A, &n, &beta, Q, &n);
		for(int i = 0; i < n; i++) {
			h[i] = std::sqrt(Q[i+n*i]);
		}

		if(iter == 0) {
			t = bmAx[0]/h[0]; 
			for(int i = 1; i < n; i++) t = std::min(t, bmAx[i]/h[i]);
			for(int i = 0; i < n; i++) {
				y[i] = y[i]/(t*t);
				h[i] = t*h[i];
				z[i] = std::max(0.1, bmAx[i]-h[i]);
			}
			alpha = t*t;
			dscal_(&nn, &alpha, Q, &inc);
		}

		sum = 0.0;
		for(int i = 0; i < n; i++) {
			yz[i] = y[i]*z[i];
			yh[i] = y[i]*h[i];
			sum += yz[i];
		}
		gap = sum/n;
		rmu = std::min(0.5, gap)*gap;
		rmu = std::max(rmu, minmu);

		alpha = -1.0; beta = 0.0;
		dgemv_(&Trans, &n, &d, &alpha, A, &n, yh, &inc, &beta, R1, &inc);
		r1 = std::max(std::abs(R1[0]), std::abs(R1[1]));
		r2 = 0; r3 = 0;
		for(int i = 0; i < n; i++) {
			R2[i] = bmAx[i] - h[i] - z[i];
			R3[i] = rmu - yz[i];
			r2 = std::max(r2, std::abs(R2[i]));
			r3 = std::max(r3, std::abs(R3[i]));
		}
		res = std::max(r1, std::max(r2, r3));
		val = E2[0]*E2[3]-E2[1]*E2[2];
		objval = std::log(val)/2;
		
		if(iprint==1) {
			std::cout << "iter = " << iter << "           " << r2 << "  " << r1 << "  " << r3 << "  " << objval << std::endl;
		}

		if( res < tol*(1+bnrm) && rmu <= minmu) {
			for(int i = 0; i < d; i++) x[i] += x0[i];
			break;
		}

		for(int i = 0; i < n; i++) {
			for(int j = 0; j < n; j++) {
				G[i+n*j] = y[i]*y[j]*Q[i+j*n]*Q[j+i*n];
			}
			G[i+i*n] += std::max(1.e-12, 2*yh[i]*z[i]);
			for(int j = 0; j < d; j++) {
				W[i+j*n] = y[i]*(h[i]+z[i])*A[i+j*n];
			}
		}

		dgetrf_(&n, &n, G, &n, IPIV, &info);
		dgetrs_(&NoTrans, &n, &d, G, &n, IPIV, W, &n, &info);

		for(int i = 0; i < n; i++) {
			for(int j = 0; j < d; j++) {
				W[i+j*n] = 2*yh[i]*W[i+j*n] - y[i]*A[i+j*n];
			}
		}

		for(int i = 0; i < n; i++) {
			R3[i] = R3[i]/y[i];
			R2[i] = R2[i] - R3[i];
		}

		alpha = 1.0; beta = 0.0;
		dgemm_(&Trans, &NoTrans, &d, &d, &n, &alpha, W, &n, A, &n, &beta, Ei, &d);
		dgemv_(&Trans, &n, &d, &alpha, W, &n, R2, &inc, &beta, h, &inc);
		daxpy_(&d, &alpha, R1, &inc, h, &inc);
		matInv2(Ei, E2);

		dgemv_(&NoTrans, &d, &d, &alpha, E2, &d, h, &inc, &beta, dx, &inc);

		dgemv_(&NoTrans, &n, &d, &alpha, A, &n, dx, &inc, &beta, Adx, &inc);
		for(int i = 0; i < n; i++) {
			h[i] = 2*yh[i]*(Adx[i]-R2[i]);
		}
		dgetrs_(&NoTrans, &n, &inc, G, &n, IPIV, h, &n, &info);
		for(int i = 0; i < n; i++) {
			dy[i] = y[i]*h[i];
			dz[i] = R3[i] - z[i]*h[i];
		}

		ax = Adx[0]/bmAx[0]; ay = h[0]; az = dz[0]/z[0];
		for(int i = 1; i < n; i++) {
			ax = std::min(ax, Adx[i]/bmAx[i]); ay = std::min(ay, h[i]); az = std::min(az, dz[i]/z[i]);
		}
		ax = std::min(-0.5, ax); ay = std::min(-0.5, ay); az = std::min(-0.5, az);
		ax = -1/ax; ay = -1/ay; az = -1/az;
		tau = std::max(tau0, 1-res);
		astep =  std::min(ax, std::min(ay, az)); astep = tau*std::min(1., astep);
		
		daxpy_(&d, &astep, dx, &inc, x, &inc);
		daxpy_(&n, &astep, dy, &inc, y, &inc);
		daxpy_(&n, &astep, dz, &inc, z, &inc);
	}


	if(iprint==1) {
		std::cout << "=========================================================" << std::endl;
		std::cout << std::endl;
	}

   }

/* ===============================================================================================
   Computes the inverse of a 2x2 matrix
   =============================================================================================== */

   void Ellipse::matInv2(double *A, double *Ainv)
   {
	double det = A[0]*A[3] - A[2]*A[1];
	det = 1.0/det;
	Ainv[0] = A[3]*det; Ainv[1] = -A[1]*det;
	Ainv[2] = -A[2]*det; Ainv[3] = A[0]*det;
   }

/* ===============================================================================================
   Computes the least square fitting helix

   inputs:
  ---------
   P   : (d x N) dimensional matrix containing N points in R^d.
  
   outputs:
  ---------
   a  : major axis of the ellipse
   b  : minor axis of the ellipse

   Implementation based on:

   	A. W. Fitzgibbon, M. Pilu and R. B. Fisher, "Direct least squares fitting of ellipses," 
	Proceedings of 13th International Conference on Pattern Recognition, 1996, pp. 253-257 vol.1, 
	doi: 10.1109/ICPR.1996.546029.

   =============================================================================================== */

   void Ellipse::lsqEllipse(double *P, int N, int d, double *a, double *b)
   {

	double *Dmat = new double[6*N];
	double *Smat = new double[6*6];
	double *Cmat = new double[6*6];
	double *VR   = new double[6*6];
	double *VL   = new double[6*6];
	int *IPIV = new int[6];
	int lwork = 200;
	double *WORK = new double[lwork];

	// Build design matrix (size 6xN)
	for(int i = 0; i < N; i++) {
		Dmat[i + 0*N] = P[2*i]*P[2*i];
		Dmat[i + 1*N] = P[2*i]*P[2*i+1];
		Dmat[i + 2*N] = P[2*i+1]*P[2*i+1];
		Dmat[i + 3*N] = P[2*i];
		Dmat[i + 4*N] = P[2*i+1];
		Dmat[i + 5*N] = 1.;
	}

	// Build scatter matrix Smat = Dmat^T Dmat (size 6x6)
	int info;
	int six = 6;
	char NoTrans = 'N'; char Trans = 'T';
	double alpha, beta;
	alpha = 1.0; beta = 0.;
	dgemm_(&Trans, &NoTrans, &six, &six, &N, &alpha, Dmat, &N, Dmat, &N, &beta, Smat, &six);

	// Build constraint matrix Cmat (size 6x6)
	memset(Cmat, 0, 6*6*sizeof(double));
	Cmat[12] = 2; Cmat[7] = -1; Cmat[2]=2;

	// Build eigen problem: matrix stored E = inv(Smat) * C stored in Dmat
	dgetrf_(&six, &six, Smat, &six, IPIV, &info);
	dgetri_(&six, Smat, &six, IPIV, WORK, &lwork, &info);

	alpha = 1.0; beta = 0.;
	dgemm_(&NoTrans, &NoTrans, &six, &six, &six, &alpha, Smat, &six, Cmat, &six, &beta, Dmat, &six);

	// Solve general eigen problem; we only need the right eigenvectors
	char JOBVL='N'; char JOBVR='V';
	dgeev_(&JOBVL, &JOBVR, &six, Dmat, &six, Cmat, Smat, VL, &six, VR, &six, WORK, &lwork, &info);

	// Find the positive eigenvalue
	int idx = -1;
	for(int i = 0; i < 6; i++) {
		if(Cmat[i] > 0) {
			idx = i;
			break;
		}
	}

	// Find coefficient of ellipse: A X^2 + B XY C Y^2 + D X + E Y + F
	// as the (right) eigenvector corresponding to the positive eigenvalue
	double A = VR[0+6*idx]; double B = VR[1+6*idx]; double C = VR[2+6*idx];
	double D = VR[3+6*idx]; double E = VR[4+6*idx]; double F = VR[5+6*idx];

	// Compute major axis (a) and minor axis (b) for the ellipse
	// see: https://en.wikipedia.org/wiki/Ellipse
	double den = B*B - 4*A*C;
	double val1 = (A-C)*(A-C) + B*B;
	double valp = A+C + std::sqrt(val1);
	double valn = A+C - std::sqrt(val1);
	double val2 = A*E*E + C*D*D - B*D*E + (B*B-4*A*C)*F;
	double nump = std::sqrt(2*val2*valp);
	double numn = std::sqrt(2*val2*valn);

	val1 = -nump/den;
	valp = -numn/den;

	*a = std::max(val1,valp);
	*b = std::min(val1,valp);

	// Clean up
	delete [] Dmat; delete [] Smat; delete [] Cmat; delete [] VL; delete [] VR; delete [] WORK;

   }

/* ===============================================================================================
   Compute largest enclosed ellipse for a polygon
   =============================================================================================== */

   void Ellipse::EllipseMax(Polygon& curve, double *a, double *b)
   {

	// Define point sets P

	int n = curve.vertices.size();

	double *P = new double[2*n];

	for(int i = 0; i < n; i++) {
		P[2*i+0] = curve.vertices[i].position[0];
		P[2*i+1] = curve.vertices[i].position[1];
	}

	int d = 2;
	double *A = new double[d*d];
	double *c = new double[d];
	double *ellipse = new double[1000];

	// Find maximum volume enclosed ellipse for curve 1
	double tol = 1.e-6;
	MaxVolEllipse(P, n, d, A, c, tol);
	Ellipse2D(A, c, a, b, ellipse);

	delete [] P; delete [] A; delete [] c; delete [] ellipse;
   }

/* ===============================================================================================
   Compute difference of aspect ratio between 2 curves based on the largest enclosed ellipse
   =============================================================================================== */

   double Ellipse::dEllipseMax(Polygon& curve1, Polygon& curve2, double *a1, double *b1, double *a2, double *b2)
   {

	EllipseMax(curve1, a1, b1);
	EllipseMax(curve2, a2, b2);

	// L1 norm of the difference of the aspect ratios
	double dist = std::abs( (*a1)/(*b1) - (*a2)/(*b2) );

	return dist;

  }
/* ===============================================================================================
   Compute smallest enclosing ellipse for a polygon
   =============================================================================================== */

   void Ellipse::EllipseMin(Polygon& curve, double *a, double *b)
   {

	// Define point sets P

	int n = curve.vertices.size();

	double *P = new double[2*n];

	for(int i = 0; i < n; i++) {
		P[2*i+0] = curve.vertices[i].position[0];
		P[2*i+1] = curve.vertices[i].position[1];
	}

	int d = 2;
	double *A = new double[d*d];
	double *c = new double[d];
	double *ellipse = new double[1000];

	// Find minimum volume enclosing ellipse for curve
	double tol = 1.e-6;
	MinVolEllipse(P, n, d, A, c, tol);
	Ellipse2D(A, c, a, b, ellipse);

	delete [] P; delete [] A; delete [] c; delete [] ellipse;
   }

/* ===============================================================================================
   Compute difference of aspect ratio between 2 curves based on the smallest enclosing ellipse
   =============================================================================================== */

   double Ellipse::dEllipseMin(Polygon& curve1, Polygon& curve2, double *a1, double *b1, double *a2, double *b2)
   {

	EllipseMin(curve1, a1, b1);
	EllipseMin(curve2, a2, b2);

	// L1 norm of the difference of the aspect ratios
	double dist = std::abs( (*a1)/(*b1) - (*a2)/(*b2) );

	return dist;

  }
/* ===============================================================================================
   Compute least square fitting ellipse for a polygon
   =============================================================================================== */

   void Ellipse::EllipseLSQ(Polygon& curve, double *a, double *b)
   {

	// Define point sets P

	int n = curve.vertices.size();

	double *P = new double[2*n];

	for(int i = 0; i < n; i++) {
		P[2*i+0] = curve.vertices[i].position[0];
		P[2*i+1] = curve.vertices[i].position[1];
	}

	int d = 2;

	// Find least square ellipse passing through Point set
	lsqEllipse(P, n, d, a, b);

	delete [] P;

   }

/* ===============================================================================================
   Compute difference of aspect ratio between 2 curves based on the lsq ellipse
   =============================================================================================== */

   double Ellipse::dEllipseLSQ(Polygon& curve1, Polygon& curve2, double *a1, double *b1, double *a2, double *b2)
   {

	EllipseLSQ(curve1, a1, b1);
	EllipseLSQ(curve2, a2, b2);

	// L1 norm of the difference of the aspect ratios
	double dist = std::abs( (*a1)/(*b1) - (*a2)/(*b2) );

	return dist;

  }
} //end namespace correa
#endif
