/* ====OT1_H ===================================================================================
 *
 * Author: Patrice Koehl (in collaboration with Henri Orland), November 2018
 * Department of Computer Science
 * University of California, Davis
 *
 * This file implements different methods needed to solve the optimal transport problem using
 * the minimization of a free energy
 *
 =============================================================================================== */

#ifndef _OT1_H_
#define _OT1_H_

  #include <cmath>
  #include <algorithm>
  #include <functional>
  #include "VectorOps.h"
  #include "ConjGradOT.h"

ConjGradOT conjgrad;

/* ===============================================================================================
   prototypes for BLAS and LAPACK functions 
   =============================================================================================== */

  extern "C" {

	// BLAS1: copy one vector into another
	void dcopy_(int * n, double * X, int *incx, double * Y, int *incy);

	// BLAS1: scale a vector
	void dscal_(int * n, double *scale, double *Y, int *incy);

	// BLAS1: dot product of two vectors	
	double ddot_(int * n, double * u, int * incu, double * v, int *incv);

	// BLAS1: norm of a vector
	double dnrm2_(int * n, double * X, int *incx);

	// BLAS2: perform Y := alpha*op( A )* B  + beta*Y
	void dgemv_(char * trans, int * m, int * n, double * alpha, double *A,
		int *lda, double * X, int * incx, double *beta, double * Y, int * incy);

	// BLAS3: perform C := alpha*op( A )* op(B)  + beta*C
	void dgemm_(char * transa, char * transb, int * m, int * n, int * k,
		double * alpha, double * A, int * lda,
		double * B, int * ldb, double * beta, double * C, int * ldc);

	// LAPACK: solve a real system of linear equations A * X = B, where A is a symmetric matrix 
	void dsysv_(char *uplo, int *n, int *nrhs, double *A, int *lda, int *ipiv, double *b,
		int *ldb, double *work, int *lwork, int *info);

	// LAPACK: diagonalize a symmetric matrix
        void dsyevd(char * JOBZ, char * UPLO, int *N, double *A, int *LDA, double *W,
        double *WORK, int *LWORK, int *IWORK, int *LIWORK, int *INFO);

        // LAPACK: solve a real system of linear equations A * X = B, where A is a symmetric matrix 
	void dsysvx_(char *fact, char *uplo, int *n, int *nrhs, double *A, int *lda, double *AF,
		int *ldaf, int *ipiv, double *b, int *ldb, double *x, int *ldx,
		double *rcond, double *ferr, double *berr,
		double *work, int *lwork, int *iwork, int *info);
  }

/* ===============================================================================================
   The OT1 class
   =============================================================================================== */

  class OT1{

  public:

	// Solve for G and d at infinite beta (iterative over beta)
	double ot1(int npoint1, double *m1, int npoint2, double *m2, double *C,
	double *G, double *ones, double *WorkSpace, double beta_init, double beta_final,
	int iprint, int nthreads);

	// Solve for G at a given beta value
	int solveG(int npoint1, double *m1, int npoint2, double *m2, double *C, double *G, double *lambda, 
	double *mu, double beta, double tol, double *ones, double *WorkSpace, int nthreads, int *info);

  private:

	// Check marginals: row sums and column sums of the transport plan G
	void checkMarginals(int n1, double *m1, int n2,  double *m2, double *G, double *err_row,
		double *err_col, double *ones, double *Space1D);

	// Compute the tranport plan G based on the auxiliary variables lambda and mu
	void computeG(int npoint1, int npoint2, double *C, double *G, 
	double *lambda, double *mu, double beta);

	// Check current Jacobian system
	double computeRC(int npoint1, double *m1, int npoint2, double *m2, double *C, double *G,
	double beta, double *lambda, double *mu, double *row, double *col, double *ones, int nthreads);

	// Solve Jacobian system for updates in Lambda and Mu, using direct solver
	void computedX_direct(int n1, int n2, double beta, double *C, double *G, double *lambda, double *mu, 
	double *ones, double *row, double *col, double *Space, double *Jac, double *B, double *Sol, int nthreads);

	// Solve Jacobian system for updates in Lambda and Mu, using iterative solver (CG)
	int computedX_iter(int n1, int n2, double beta, double *C, double *G, double *lambda, double *mu, 
	double *ones, double *row, double *col, double *B, double *Sol, double *Space, int nthreads);


  };

/* ===============================================================================================
   Common variables
   =============================================================================================== */

/* ===============================================================================================
   checkMarginals

   Notes:
   =====
   Check the marginals of a coupling matrix

   Input:
   =====
	n: 	number of points on space 1
	m1:	measure on space 1
	n2:	number of points on space 2
	m2:	measure on space 2
	G:	coupling matrix
  Output:
	err_row: error on row marginals
	err_col: error on col marginals
   =============================================================================================== */

  void OT1::checkMarginals(int n1, double *m1, int n2,  double *m2, double *G, double *err_row,
		double *err_col, double *ones, double *Space1D)
  {

	double alpha, beta;
	alpha = 1.0; beta = 0.0;
	int inc = 1;
	char Trans   = 'T';
	char NoTrans = 'N';

	dgemv_(&NoTrans, &n1, &n2, &alpha, G, &n1, ones, &inc, &beta, Space1D, &inc);
	for(int i = 0; i <n1; i++) { Space1D[i] = Space1D[i] - m1[i];};
	double val = ddot_(&n1, Space1D, &inc, Space1D, &inc);
	*err_row = std::sqrt(val);

	dgemv_(&Trans, &n1, &n2, &alpha, G, &n1, ones, &inc, &beta, Space1D, &inc);
	for(int i = 0; i <n2; i++) { Space1D[i] = Space1D[i] - m2[i];};
	val = ddot_(&n2, Space1D, &inc, Space1D, &inc);
	*err_col = std::sqrt(val);
  }
/* ===============================================================================================
   computeG

   Input:
   =====
	npoint1:	number of points for point set 1
	npoint2:	number of points for point set 2
	C:		cost matrix
	beta:		current beta
	lambda:		current values of Lagragians lambda
	mu:		current values of Lagragians mu
	nthreads:	number of threads for parallel computation
	
   Output:
   ======
	G:		coupling matrix
   =============================================================================================== */

  void OT1::computeG(int npoint1, int npoint2, double *C, double *G, 
	double *lambda, double *mu, double beta)
  {

	double val;
	double tol = 1.e-8;

	for(int j = 0; j < npoint2; j++) {
		for(int i = 0; i < npoint1; i++) {
			val = C[i+j*npoint1] + lambda[i] + mu[j];
			val = coupling(val, beta);
			if(val < tol) val = 0;
			if(val > 1.- tol) val = 1;
			G[i+npoint1*j] = val;
		}
	}

  }

/* ===============================================================================================
   computeRC

   Input:
   =====
	npoint1:	number of points for point set 1
	m1:		measure on points1
	npoint2:	number of points for point set 2
	m2:		measure on points2
	C:		cost matrix
	beta:		current beta
	lambda:		current values of Lagragians lambda
	mu:		current values of Lagragians mu
	nthreads:	number of threads for parallel computation
	
   Output:
   ======
	row:		errors on row sums
	col:		errors on col sums
	err_r:		total error on row sums
	err_m:		total error on col sums
   =============================================================================================== */

  double OT1::computeRC(int npoint1, double *m1, int npoint2, double *m2, double *C, double *G,
	double beta, double *lambda, double *mu, double *row, double *col, double *ones, int nthreads)
  {

	int n1n2 = npoint1*npoint2;
	char Trans   = 'T';
	char NoTrans = 'N';
	int inc = 1;
	int one = 1;
	double a, b;
	
	dcopy_(&n1n2, C, &inc, G, &inc);
	a = 1.0; b = 1;
	dgemm_(&NoTrans, &Trans, &npoint1, &npoint2, &one, &a, lambda, &npoint1, ones, &npoint2, &b, 
		G, &npoint1);
	dgemm_(&NoTrans, &Trans, &npoint1, &npoint2, &one, &a, ones, &npoint1, mu, &npoint2, &b, 
		G, &npoint1);

	vect_coupling(n1n2, G, beta, nthreads);

	dcopy_(&npoint1, m1, &inc, row, &inc);
	a = 1; b = -1;
	dgemv_(&NoTrans, &npoint1, &npoint2, &a, G, &npoint1, ones, &inc, &b, row, &inc);
	dcopy_(&npoint2, m2, &inc, col, &inc);
	dgemv_(&Trans, &npoint1, &npoint2, &a, G, &npoint1, ones, &inc, &b, col, &inc);

	double err_l = 0;
	double val = ddot_(&npoint1, row, &inc, row, &inc);
	err_l = std::sqrt(val);
//	for(int i = 0; i < npoint1; i++) err_l += std::abs(row[i]);

	double err_m = 0;
	int n2 = npoint2-1;
	val = ddot_(&n2, col, &inc, col, &inc);
	err_m = std::sqrt(val);
//	for(int j = 0; j < npoint2-1; j++) err_m += std::abs(col[j]);

	double err = err_l+err_m;

	return err;

  }

/* ===============================================================================================
   computedX_direct

   Input:
   =====
	npoint1:	number of points for point set 1
	npoint2:	number of points for point set 2
	beta:		current beta
	C:		Cost matrix
	lambda:		current values of Lagragians lambda
	mu:		current values of Lagragians mu
	row:		errors on row marginals
	col:		errors on col marginals
	nthreads:	number of threads for parallel computation
	
   Output:
   ======
	drow:		updates on lambda
	dcol:		updates on mu

   Method:
   ======
	We write the Jacobian as:
        Jac = [A     B]
              [B^T   D]
	where A and D are diagonal matrices
	To solve:
	Jac dF = B
	we write it as:
	A x + B y = a				(1)
	B^T x + D y = b				(2)
	We multiply (2) by BD^{-1}
	A x + B y = a
	BD^{-1}B^T x + B y = BD^{-1}b           (2')
	and (1) - (2') gives:
	(A - BD^{-1}B^T) x = a - BD^{-1}b
	which we solve using dsysv, and then:
	y = D^{-1}b - B^Tx

	Note that if dim(a) > dim(b), we reverse the system
	
   =============================================================================================== */

  void OT1::computedX_direct(int n1, int n2, double beta, double *C, double *G, double *lambda, double *mu, 
	double *ones, double *row, double *col, double *Space, double *Jac, double *B, double *Sol, int nthreads)
  {

	int npoint1 = n1;
	int npoint2 = n2+1;
	int n1n2 = npoint1*npoint2;
	int n1pn2 = n1+n2;
	char Trans   = 'T';
	char NoTrans = 'N';
	int inc = 1;
	int one = 1;
	double a, b;

	dcopy_(&n1, row, &inc, B, &inc);
	dcopy_(&n2, col, &inc, &B[npoint1], &inc);
	a = -1;
	dscal_(&n1pn2, &a, B, &inc);
	
	dcopy_(&n1n2, C, &inc, G, &inc);
	a = 1.0; b = 1;
	dgemm_(&NoTrans, &Trans, &npoint1, &npoint2, &one, &a, lambda, &npoint1, ones, &npoint2, &b, 
		G, &npoint1);
	dgemm_(&NoTrans, &Trans, &npoint1, &npoint2, &one, &a, ones, &npoint1, mu, &npoint2, &b, 
		G, &npoint1);

	vect_dcoupling(n1n2, G, beta, nthreads);

	a = 1; b = 0;
	dgemv_(&NoTrans, &npoint1, &npoint2, &a, G, &npoint1, ones, &inc, &b, row, &inc);
	dgemv_(&Trans, &npoint1, &npoint2, &a, G, &npoint1, ones, &inc, &b, col, &inc);

	int nrhs = 1; char U  = 'L'; int info;
	int lwork;
	int *IPIV = new int[std::max(npoint1, npoint2)];

	if(n1 <= n2+1) {

		for(int j = 0; j < n2; j++) col[j] = 1.0/col[j];

		for(int j = 0; j < n2; j++) {
			for(int i = 0; i < n1; i++) {
				Space[i+j*n1] = G[i+j*n1]*col[j];
			}
		}

		a = -1.0; b = 0;
		dgemm_(&NoTrans, &Trans, &n1, &n1, &n2, &a, Space, &n1, G, &n1, &b, Jac, &n1);
		for(int i = 0; i < n1; i++) Jac[i+n1*i] += row[i];

		a = -1.0; b = 1.0;
		dgemv_(&NoTrans, &n1, &n2, &a, Space, &n1, &B[npoint1], &inc, &b, B, &inc);

		lwork = 128*n1;
		dsysv_(&U, &n1, &nrhs, Jac, &n1, IPIV, B, &n1, Space, &lwork, &info);

		a = -1.0; b = 1.0;
		dgemv_(&Trans, &n1, &n2, &a, G, &n1, B, &inc, &b, &B[npoint1], &inc);

		for(int j = 0; j < n2; j++) B[npoint1+j] = col[j]*B[npoint1+j];
		dcopy_(&n1pn2, B, &inc, Sol, &inc);

	} else {

		for(int j = 0; j < n1; j++) row[j] = 1.0/row[j];

		for(int j = 0; j < n2; j++) {
			for(int i = 0; i < n1; i++) {
				Space[i+j*n1] = G[i+j*n1]*row[i];
			}
		}

		a = -1.0; b = 0;
		dgemm_(&Trans, &NoTrans, &n2, &n2, &n1, &a, G, &n1, Space, &n1, &b, Jac, &n2);
		for(int i = 0; i < n2; i++) Jac[i+n2*i] += col[i];

		a = 1.0; b = -1.0;
		dgemv_(&Trans, &n1, &n2, &a, Space, &n1, B, &inc, &b, &B[npoint1], &inc);

		lwork = 128*n2;
		dsysv_(&U, &n2, &nrhs, Jac, &n2, IPIV, &B[npoint1], &n2, Space, &lwork, &info);


		a = -1.0; b = -1.0;
		dgemv_(&NoTrans, &n1, &n2, &a, G, &n1, &B[npoint1], &inc, &b, B, &inc);

		for(int i = 0; i < n1; i++) B[i] = row[i]*B[i];
		dcopy_(&n1pn2, B, &inc, Sol, &inc);
	}

  }

/* ===============================================================================================
   computedX_iter

   Input:
   =====
	npoint1:	number of points for point set 1
	npoint2:	number of points for point set 2
	beta:		current beta
	C:		Cost matrix
	lambda:		current values of Lagragians lambda
	mu:		current values of Lagragians mu
	row:		errors on row marginals
	col:		errors on col marginals
	nthreads:	number of threads for parallel computation
	
   Output:
   ======
	drow:		updates on lambda
	dcol:		updates on mu

   Method:
   ======
	We write the Jacobian as:
        Jac = [A     B]
              [B^T   D]
	where A and D are diagonal matrices
	We solve:
	Jac dF = B
	using preconditioned conjugate gradient

   =============================================================================================== */

  int OT1::computedX_iter(int n1, int n2, double beta, double *C, double *G, double *lambda, double *mu, 
	double *ones, double *row, double *col, double *B, double *Sol, double *Space, int nthreads)
  {

	int npoint1 = n1;
	int npoint2 = n2+1;
	int n1n2 = npoint1*npoint2;
	int n1pn2 = n1+n2;
	char Trans   = 'T';
	char NoTrans = 'N';
	int inc = 1;
	int one = 1;
	double a, b;
	double tol = 1.e-4;

	dcopy_(&n1, row, &inc, B, &inc);
	dcopy_(&n2, col, &inc, &B[npoint1], &inc);
	a = -1;
	dscal_(&n1pn2, &a, B, &inc);
	
	dcopy_(&n1n2, C, &inc, G, &inc);
	a = 1.0; b = 1;
	dgemm_(&NoTrans, &Trans, &npoint1, &npoint2, &one, &a, lambda, &npoint1, ones, &npoint2, &b, 
		G, &npoint1);
	dgemm_(&NoTrans, &Trans, &npoint1, &npoint2, &one, &a, ones, &npoint1, mu, &npoint2, &b, 
		G, &npoint1);

	vect_dcoupling(n1n2, G, beta, nthreads);

	a = 1; b = 0;
	dgemv_(&NoTrans, &npoint1, &npoint2, &a, G, &npoint1, ones, &inc, &b, row, &inc);
	dgemv_(&Trans, &npoint1, &npoint2, &a, G, &npoint1, ones, &inc, &b, col, &inc);

	int niter = conjgrad.cgDriver(n1, n2, row, col, G, B, Sol, &Space[0], &Space[n1pn2],
	&Space[2*n1pn2], &Space[3*n1pn2], &Space[4*n1pn2], &Space[5*n1pn2], tol);

	return niter;

  }

/* ===============================================================================================
   solveG

   Input:
   =====
	npoint1:	number of points for point set 1
	m1:		measure on points1
	npoint2:	number of points for point set 2
	m2:		measure on points2
	C:		cost matrix
	beta:		parameter beta
	tol:		tolerance criteria
	nthreads:	number of threads for parallel computation
	
   Output:
   ======
	G:		coupling matrix
   =============================================================================================== */

  int OT1::solveG(int npoint1, double *m1, int npoint2, double *m2, double *C, double *G, double *lambda, 
	double *mu, double beta, double tol, double *ones, double *WorkSpace, int nthreads, int *info)
  {

	int nmax = std::max(npoint1, npoint2);
	double *row = &WorkSpace[0];
	double *col = &WorkSpace[npoint1];
	double *B   = &WorkSpace[npoint1+npoint2];
	double *Sol = &WorkSpace[2*(npoint1+npoint2)];
	double *Space = &WorkSpace[3*(npoint1+npoint2)];
	double *Jac   = &WorkSpace[3*(npoint1+npoint2) + nmax*nmax];

	int n1 = npoint1;
	int n2 = npoint2-1;

	double err, err_old;

	err=computeRC(npoint1, m1, npoint2, m2, C, G, beta, lambda, mu, row, col, ones, nthreads);
//	std::cout << "err0 = " << err << std::endl;

	*info = 0;
	int iter = 0;
	int nstep;
	double step;
	while (err > tol)
	{
		err_old = err;
//		int it = computedX_iter(n1, n2, beta, C, G, lambda, mu,
//		ones, row, col, B, Sol, Space, nthreads);
		computedX_direct(n1, n2, beta, C, G, lambda, mu,
		ones, row, col, Space, Jac, B, Sol, nthreads);

		step = 1.0;
		nstep = 0;
		for(int i = 0; i < 30; i++) {
			memset(Space, 0, (npoint1+npoint2)*sizeof(double));
			for(int i = 0; i < npoint1; i++) Space[i] = lambda[i] + step*Sol[i];
			for(int j = 0; j < npoint2-1; j++) Space[npoint1+j] = mu[j] + step*Sol[npoint1+j];
			err=computeRC(npoint1, m1, npoint2, m2, C, G, beta, Space, &Space[npoint1],
			row, col, ones, nthreads);
			if(err < err_old) break;
			step = step/2;
			nstep++;
		}
//		std::cout << "err_old = " << err_old << " err = " << err << " nstep = " << nstep << " it = " << iter << std::endl;
		if(iter==100 || nstep==30) {
			*info = 1;
			break;
		}

		for(int i = 0; i < npoint1; i++) lambda[i] = Space[i];
		for(int j = 0; j < npoint2-1; j++) mu[j] = Space[npoint1+j];

		iter++;
	}

	int niter = iter;

	computeG(npoint1, npoint2, C, G, lambda, mu, beta);

	return niter;

  }
		
/* ===============================================================================================
   OT1 (Wasserstein / Earth Mover)

   Input:
   =====
	npoint1:	number of points for point set 1
	m1:		measure on points1
	npoint2:	number of points for point set 2
	m2:		measure on points2
	C:		cost matrix
	beta1:		starting beta
	nthreads:	number of threads for parallel computation

	ones:		vector of ones, of size max(npoint1, npoint2)
	WorkSpace	working array of size 2*nmax*nmax + 10*nmax, where nmax=max(npoint1, npoint2)
	beta_init	starting value for beta (usually 1000)
	beta_init	final value for beta (usually 10^{10}
	iprint		flag: 
				0: no printing
				1: print info about annealing process
	
   Output:
   ======
	dist:		Optimal transport distance
	G:		coupling matrix
   =============================================================================================== */

  double OT1::ot1(int npoint1, double *m1, int npoint2, double *m2, double *C,
	double *G, double *ones, double *WorkSpace, double beta_init, double betaf,
	int iprint, int nthreads)
  {

	double *lambda = &WorkSpace[0];
	double *mu     = &WorkSpace[npoint1];

	memset(lambda, 0, npoint1*sizeof(double));
	memset(mu, 0, npoint2*sizeof(double));

	// Define all variables needed to compute Earh Mover's distance

	double beta_val, tol;
	double dist, err_row, err_col;
	int n1n2 = npoint1*npoint2;
	int inc = 1;

	dist = 1.0;

	int niter;
	int nitermax = 0;

	beta_val = beta_init;

	if(iprint == 1) {
		std::cout << " " << std::endl;
		std::cout << "        " << "=====================================================================================================" << std::endl;
		std::cout << "        " << "       Beta           Iter              U                  Ent             Err_row         Err_col   " << std::endl;
		std::cout << "        " << "=====================================================================================================" << std::endl;
	}

	int info;
	double coef = std::sqrt(std::sqrt(10));

	while(beta_val < betaf)
	{

		tol = std::max(1./beta_val, 1.e-5);
		tol = 1.e-5;
		niter = solveG(npoint1, m1, npoint2, m2, C, G, lambda, mu, beta_val, tol,
			ones, &WorkSpace[npoint1+npoint2], nthreads, &info);
		nitermax += niter;

		checkMarginals(npoint1, m1, npoint2,  m2, G, &err_row, &err_col, ones,
			&WorkSpace[npoint1+npoint2]);

		dist = ddot_(&n1n2, G, &inc, C, &inc); 
		double ent = 0;
		for(int i = 0; i < npoint1*npoint2; i++) {
			if(G[i]>0.0) ent -= G[i]*std::log(G[i]);
		}

		if(iprint==1) {
			std::cout << "        " << "   " << std::setw(10)<< beta_val << "    ";
			std::cout << std::setw(10) << niter << "        " << std::setw(10) << dist << "        ";
			std::cout << std::setw(10) << ent <<  "        ";
			std::cout << std::setw(10) << err_row <<  "        " << err_col << std::endl;
		}

		beta_val = beta_val*coef;
	}


	if(iprint==1) {
		std::cout << "        " << "=====================================================================================================" << std::endl;
		std::cout << " " << std::endl;
	}

	computeG(npoint1, npoint2, C, G, lambda, mu, betaf);

	return dist;

  }

#endif
