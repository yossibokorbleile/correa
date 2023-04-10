/* ====CONJGRADOT.H ==========================================================================
 *
 * Author: Patrice Koehl (in collaboration with Henri Orland), November 2018
 * Department of Computer Science
 * University of California, Davis
 *
 * This file implements a conjugate gradient method for soling a linear system of equation with
 * preconditioner
 *
 * It is based on cg_rc.cpp from John Burkardt 
 * (https://people.sc.fsu.edu/~jburkardt/cpp_src/cg_rc/cg_rc_prb.cpp)
 *
 =============================================================================================== */

#ifndef _CONJGRADOT_H_
#define _CONJGRADOT_H_

  #include <cmath>
  #include <algorithm>
  #include <functional>
  #include <cstdlib>

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
	// BLAS2: perform Y := alpha*op( A )* B  + beta*Y
	void dgemv_(char * trans, int * m, int * n, double * alpha, double *A,
		int *lda, double * X, int * incx, double *beta, double * Y, int * incy);
  }

/* ===============================================================================================
   Define the conjgrad class for OT
   =============================================================================================== */

  class ConjGradOT{

  public:

	// driver function for conjugate gradient minimizer
	int cgDriver(int n1, int n2, double *D1, double *D2, double *G, double *B, double *X,
		double *r, double *z, double *p, double *q, double *Temp, double *Temp2, double tol);

  private:

	// compute Y = A X
	void computeMatVect(int n1, int n2, double *D1, double *D2, double *G, double *X, double *Y,
			double *Temp);

	// apply preconditioner
	void preConditioner(int n1, int n2, double *D1, double *D2, double *G, double *B, double *X, double *wrk);

	// conjgrad routine with reverse communication
	int cgRC( int n, double *b, double *x, double *r, double *z, double *p, double *q, int job );

	// Internal variables for conjgrad
	int rlbl, iter;
	double rho, rho_old;

  };

/* ===============================================================================================
   Compute Y = A*X
   where A is of the form:

	A = [ D1    G ]
	    [ G^T   D2 ]

    with D1 a n1xn1 diagonal matrix, G a n1xn2 matrix, and D2 a n2xn2 diagonal matrix

    Temp is a work space array
   =============================================================================================== */

  void ConjGradOT::computeMatVect(int n1, int n2, double *D1, double *D2, double *G, double *X, double *Y,
			double *Temp)
  {

	int inc = 1;
	double a, b;

	char Trans   = 'T';
	char NoTrans = 'N';

	a = 1; b = 0;
	dgemv_(&NoTrans, &n1, &n2, &a, G, &n1, &X[n1], &inc, &b, Temp, &inc);

	for(int i = 0; i < n1; i++)
	{
		Y[i] = D1[i]*X[i] + Temp[i];
	}

	dgemv_(&Trans, &n1, &n2, &a, G, &n1, X, &inc, &b, Temp, &inc);

	for(int i = 0; i < n2; i++)
	{
		Y[i+n1] = D2[i]*X[i+n1] + Temp[i];
	}
  }

/* ===============================================================================================
   Precond

   Solve 
	M X = B
   where M is the chosen preconditionner

   Input:
   =====
	n1:	number of points for point set 1
	n2:	number of points for point set 2
	D1:	first diagonal matrix (see above)
	D2:	second diagonal matrix (see above)
	G:	upper matrix (see above)

   M = [ I             0 ] [ D1   0 ] [ I   D1^{-1} G ]
       [ G^T D1^{-1}   I ] [ 0   D2 ] [ 0           I ]
	
   Output:
   ======
	dX:		updates on lambda
	dm:		updates on mu

   =============================================================================================== */

  void ConjGradOT::preConditioner(int n1, int n2, double *D1, double *D2, double *G, double *B, double *X, double *wrk)
  {


	int inc = 1;
	double a, b;

	char Trans   = 'T';
	char NoTrans = 'N';

	for(int i = 0; i < n1; i++) {
		wrk[i] = B[i] / D1[i];
	}
	a = 1; b = 0;
        dgemv_(&Trans, &n1, &n2, &a, G, &n1, wrk, &inc, &b, &wrk[n1], &inc);
	for(int i = 0; i < n2; i++) {
		X[i+n1] = (B[i+n1]-wrk[i+n1])/D2[i];
	}
        dgemv_(&NoTrans, &n1, &n2, &a, G, &n1, &X[n1], &inc, &b, wrk, &inc);
	for(int i = 0; i < n1; i++) {
		X[i] = (B[i] - wrk[i]) / D1[i];
	}

/*
	for(int i = 0; i < n1; i++) {
		X[i] = B[i] / D1[i];
	}
	for(int i = 0; i < n2; i++) {
		X[i+n1] = B[i+n1] / D2[i];
	}

	for(int i = 0; i<n1+n2; i++) {
		X[i] = B[i];
	}
*/
  }
/* ===============================================================================================
   Conjugate gradient with reverse communication
   =============================================================================================== */

  int ConjGradOT::cgRC( int n, double *b, double *x, double *r, double *z, double *p, double *q, int job )

/* ===============================================================================================
* 
*   Purpose:
* 
*     CG_RC is a reverse communication conjugate gradient routine.
* 
*   Discussion:
* 
*     This routine seeks a solution of the linear system A*x=b
*     where b is a given right hand side vector, A is an n by n
*     symmetric positive definite matrix, and x is an unknown vector
*     to be determined.
* 
*     Under the assumptions that the matrix A is large and sparse,
*     the conjugate gradient method may provide a solution when
*     a direct approach would be impractical because of excessive
*     requirements of storage or even of time.
* 
*     The conjugate gradient method presented here does not require the 
*     user to store the matrix A in a particular way.  Instead, it only 
*     supposes that the user has a way of calculating
*       y = alpha * A * x + b * y
*     and of solving the preconditioned linear system
*       M * x = b
*     where M is some preconditioning matrix, which might be merely
*     the identity matrix, or a diagonal matrix containing the
*     diagonal entries of A.
* 
*     This routine was extracted from the "templates" package.
*     There, it was not intended for direct access by a user;
*     instead, a higher routine called "cg()" was called once by
*     the user.  The cg() routine then made repeated calls to 
*     cgrevcom() before returning the result to the user.
* 
*     The reverse communication feature of cgrevcom() makes it, by itself,
*     a very powerful function.  It allows the user to handle issues of
*     storage and implementation that would otherwise have to be
*     mediated in a fixed way by the function argument list.  Therefore,
*     this version of cgrecom() has been extracted from the templates
*     library and documented as a stand-alone procedure.
* 
*     The user sets the value of JOB to 1 before the first call,
*     indicating the beginning of the computation, and to the value of
*     2 thereafter, indicating a continuation call.  
*     The output value of JOB is set by cgrevcom(), which
*     will return with an output value of JOB that requests a particular
*     new action from the user.
* 
*   Licensing:
* 
*     This code is distributed under the GNU LGPL license.
* 
*   Modified:
* 
*     13 January 2013
* 
*   Author:
* 
*     John Burkardt
* 
*   Reference:
* 
*     Richard Barrett, Michael Berry, Tony Chan, James Demmel,
*     June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
*     Charles Romine, Henk van der Vorst,
*     Templates for the Solution of Linear Systems:
*     Building Blocks for Iterative Methods,
*     SIAM, 1994,
*     ISBN: 0898714710,
*     LC: QA297.8.T45.
* 
*   Parameters:
* 
*     Input, int N, the dimension of the matrix.
* 
*     Input, double B[N], the right hand side vector.
* 
*     Input/output, double X[N].  On first call, the user 
*     should store an initial guess for the solution in X.  On return with
*     JOB = 4, X contains the latest solution estimate.
* 
*     Input/output, double R[N], Z[N], P[N], Q[N],
*     information used by the program during the calculation.  The user
*     does not need to initialize these vectors.  However, specific
*     return values of JOB may require the user to carry out some computation
*     using data in some of these vectors.
* 
*     Input/output, int JOB, communicates the task to be done.
*     The user needs to set the input value of JOB to 1, before the first call,
*     and then to 2 for every subsequent call for the given problem.
*     The output value of JOB indicates the requested user action.  
*     * JOB = 1, compute Q = A * P;
*     * JOB = 2: solve M*Z=R, where M is the preconditioning matrix;
*     * JOB = 3: compute R = R - A * X;
*     * JOB = 4: check the residual R for convergence.  
*                If satisfactory, terminate the iteration.
*                If too many iterations were taken, terminate the iteration.
* 
 =============================================================================================== */

 {

	int job_next=0;
	int inc =1;

	double alpha, beta;
	double a;
	double pdotq;

/* ===============================================================================================
   Initialisation: ask user to compute the initial residual
   =============================================================================================== */

	if ( job == 1 ) {

		dcopy_(&n, b, &inc, r, &inc);
		job_next = 3;
		rlbl = 2;
	}

/* ===============================================================================================
   Begin first conjugate gradient loop. Ask the user for a preconditioner solve
   =============================================================================================== */

	else if ( rlbl == 2 ) {

		iter = 1;
		job_next = 2;
		rlbl = 3;
	}

/* ===============================================================================================
   Compute the direction; ask the user to compute ALPHA; save A*P to Q
   =============================================================================================== */

	else if ( rlbl == 3 ) {

		rho = ddot_(&n, r, &inc, z, &inc);

		if ( 1 < iter ) {

			beta = rho / rho_old;
			daxpy_(&n, &beta, p, &inc, z, &inc);
		}

		dcopy_(&n, z, &inc, p, &inc);

		job_next = 1;
		rlbl = 4;
	}

/* ===============================================================================================
   Compute current solution vector. Ask the user to check stopping criteria
   =============================================================================================== */

	else if ( rlbl == 4 ) {

		pdotq = ddot_(&n, p, &inc, q, &inc);
		alpha = rho / pdotq;
		a = alpha;
		daxpy_(&n, &a, p, &inc, x, &inc);
		a = -alpha;
		daxpy_(&n, &a, q, &inc, r, &inc);

		job_next = 4;
		rlbl = 5;
	}

/* ===============================================================================================
   Begin the next step; ask for a preconditioner solve
   =============================================================================================== */

	else if ( rlbl == 5 )
	{
		rho_old = rho;
		iter = iter + 1;

		job_next = 2;
		rlbl = 3;
	}

	return job_next;
}

/* ===============================================================================================
   Conjugate gradient driver
   =============================================================================================== */

   int ConjGradOT::cgDriver(int n1, int n2, double *D1, double *D2, double *G, double *B, double *X,
		double *r, double *z, double *p, double *q, double *Temp, double *Temp2, double tol)
   {

	int it, it_max;
	int job;
	int inc = 1;
	int n1pn2 = n1 + n2;

	double a;
	double bnrm2, rnrm2;

	double *X_old = new double[n1pn2];

/* ===============================================================================================
   Set initial guess for the solution
   =============================================================================================== */

	memset(X, 0, n1pn2*sizeof(double));
	memset(r, 0, n1pn2*sizeof(double));
	memset(z, 0, n1pn2*sizeof(double));
	memset(p, 0, n1pn2*sizeof(double));
	memset(q, 0, n1pn2*sizeof(double));

/* ===============================================================================================
   Parameters for the stopping test
   =============================================================================================== */

	it = 0;
	it_max = n1pn2;
	it_max = 100;
	bnrm2 = dnrm2_(&n1pn2, B, &inc);

/* ===============================================================================================
   Repeatedly call the cg_rc code, and, on return, do what JOB tells us
   =============================================================================================== */

	job = 1;

	int nm, npred;
	nm = 0;
	npred = 0;

	dcopy_(&n1pn2, X, &inc, X_old, &inc);

	for ( ; ; )
	{

		job = cgRC ( n1pn2, B, X, r, z, p, q, job );

/* 		==================================================================================
		Compute q = A *p
   		================================================================================== */

		if ( job == 1 ) {
			computeMatVect(n1, n2, D1, D2, G, p, q, Temp);
			nm++;
		}

/* 		==================================================================================
		Preconditioning: solve M * z = r
   		================================================================================== */

		else if ( job == 2 ) {
			preConditioner(n1, n2, D1, D2, G, r, z, Temp);
			npred++;
		}

/* 		==================================================================================
		Compute residual: r = r - A*x
   		================================================================================== */

		else if ( job == 3 ) {

			computeMatVect(n1, n2, D1, D2, G, X, Temp, Temp2);
			nm++;
			a = -1;
			daxpy_(&n1pn2, &a, Temp, &inc, r, &inc);

		}

/* 		==================================================================================
		Stopping test on r
   		================================================================================== */

		else if ( job == 4 ) {

			rnrm2 = dnrm2_(&n1pn2, r, &inc);
//			std::cout << "rnrm2 = " << rnrm2 << " bnrm2 = " << bnrm2 << " tol = " << tol << std::endl;


//			if ( bnrm2 == 0.0 ) {
			if ( bnrm2 <= 1.e-6 ) {
				if ( rnrm2 <= tol ) {
					break;
				}
			} else {
				if ( rnrm2 <= tol * bnrm2 ) {
					break;
				}
			}

			it = it + 1;
			if(it > it_max) {
				break;
			}

		}

		job = 2;

	}

	return it;
  
  }

#endif
