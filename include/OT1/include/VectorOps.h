/*================================================================================================
  VectorOps.h
  Version 1: 11/12/2018

  Purpose: Sets of routine for performing operations on vectors using pthreads

Copyright (c) Patrice Koehl.

>>> SOURCE LICENSE >>>

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

>>> END OF LICENSE >>>

================================================================================================== */

#ifndef _VECTOROPS_H_
#define _VECTOROPS_H_

/*================================================================================================
 Includes
================================================================================================== */

#include <math.h>
#include <cstdlib>
#include <pthread.h>

/*================================================================================================
 Definitions for multi-threading
================================================================================================== */

typedef struct thread_data {
	int Nval;
	double fact;
	double offset;
	double *vect;
} thread_data;

thread_data thread_datas[NUM_THREADS];

/*================================================================================================
 Function Definitions
================================================================================================== */

auto coupling = [](double xb, double beta)
{ 
	double vexp, val;
	double x = xb*beta;
	double tol = 1.e-10;
	if(std::abs(x) < tol) {
		val = 0.5;
	} else if(x < 0) {
		vexp = std::exp(x);
		val = 1./(1-vexp) + 1/x;
	} else {
		vexp = std::exp(-x);
		val = vexp/(vexp-1.0) + 1/x;
	}
	return val;
};

auto dcoupling = [](double xb, double beta)
{ 
	double vexp, val;
	double x = xb*beta;
	double tol = 1.e-10;
	if(std::abs(x) < tol) {
		val = -1.0/12;
	} else if(x < 0) {
		vexp = std::exp(x);
		val = vexp/((1-vexp)*(1-vexp)) - 1/(x*x);
	} else {
		vexp = std::exp(-x);
		val = vexp/((1-vexp)*(1-vexp)) - 1/(x*x);
	}
	return val;
};

/*================================================================================================
 Vect_exp_thread: computes the element-wise exponential of a vector
================================================================================================== */

void* vect_coupling_thread(void* data)
{
	int threadid = *((int *) data);

	int Nval = thread_datas[threadid].Nval;
	double beta = thread_datas[threadid].fact;

	double val;
	for (int i = 0; i < Nval; i++)
	{
		val = thread_datas[threadid].vect[i];
		val = coupling(val, beta);
		thread_datas[threadid].vect[i] = val;
	}
	return 0;
}

/*================================================================================================
 Vect_abs_thread: computes the element-wise absolute value of a vector
================================================================================================== */

void* vect_dcoupling_thread(void* data)
{
	int threadid = *((int *) data);

	int Nval = thread_datas[threadid].Nval;
	double beta = thread_datas[threadid].fact;

	double val;
	for (int i = 0; i < Nval; i++)
	{
		val = thread_datas[threadid].vect[i];
		val = beta*dcoupling(val, beta);
		thread_datas[threadid].vect[i] = val;
	}
	return 0;
}

/*================================================================================================
 vect_coupling: compute the element-wise exponential of an array using pthreads
================================================================================================== */

void vect_coupling(int N, double *X, double fact, int nthreads)
{
	int nval = N / nthreads;
	int nval2 = nval*(nthreads-1);

/*      ==========================================================================================
	Generate all threads
        ========================================================================================== */

	for(int i = 0; i < nthreads-1; i++) 
	{
		threadids[i]=i;
		thread_datas[i].Nval = nval;
		thread_datas[i].fact = fact;
		thread_datas[i].vect  = &X[nval*i];

		pthread_create(&threads[i], NULL, vect_coupling_thread, (void*) &threadids[i]);
	}
	
/*      ==========================================================================================
	Compute on main thread
        ========================================================================================== */

	double val;
	for(int i = nval2; i < N; i++) {
		val = coupling(X[i],fact);
		X[i] = val;
	}

/*      ==========================================================================================
	Join all the threads (to make sure they are all finished)
        ========================================================================================== */

	for (int i=0; i < nthreads-1; i++)
	{
		pthread_join(threads[i], NULL);
	}
}
/*================================================================================================
 vect_dcoupling: compute the element-wise absolute value of an array using pthreads
================================================================================================== */

void vect_dcoupling(int N, double *X, double fact, int nthreads)
{
	int nval = N / nthreads;
	int nval2 = nval*(nthreads-1);

/*      ==========================================================================================
	Generate all threads
        ========================================================================================== */

	for(int i = 0; i < nthreads-1; i++) 
	{
		threadids[i]=i;
		thread_datas[i].Nval = nval;
		thread_datas[i].fact = fact;
		thread_datas[i].vect  = &X[nval*i];

		pthread_create(&threads[i], NULL, vect_dcoupling_thread, (void*) &threadids[i]);
	}
	
/*      ==========================================================================================
	Compute on main thread
        ========================================================================================== */

	double val;
	for(int i = nval2; i < N; i++) {
		val = fact*dcoupling(X[i],fact);
		X[i] = val;
	}

/*      ==========================================================================================
	Join all the threads (to make sure they are all finished)
        ========================================================================================== */

	for (int i=0; i < nthreads-1; i++)
	{
		pthread_join(threads[i], NULL);
	}
}

#endif
