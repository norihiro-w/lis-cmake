/* Copyright (C) 2002-2012 The SSI Project. All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
   3. Neither the name of the project nor the names of its contributors 
      may be used to endorse or promote products derived from this software 
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE SCALABLE SOFTWARE INFRASTRUCTURE PROJECT
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE SCALABLE SOFTWARE INFRASTRUCTURE
   PROJECT BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
   OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef HAVE_CONFIG_H
	#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN32_H
	#include "lis_config_win32.h"
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
        #include <malloc.h>
#endif
#include <string.h>
#include <stdarg.h>
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

#define NWORK				4
/************************************************
 * lis_jacobi_check_params
 * lis_jacobi_malloc_work
 * lis_jacobi
 ************************************************/
#undef __FUNC__
#define __FUNC__ "lis_jacobi_check_params"
int lis_jacobi_check_params(LIS_SOLVER solver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_jacobi_malloc_work"
int lis_jacobi_malloc_work(LIS_SOLVER solver)
{
	LIS_VECTOR	*work;
	int			i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	worklen = NWORK;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_jacobi_malloc_work::work" );
	if( work==NULL )
	{
		LIS_SETERR_MEM(worklen*sizeof(LIS_VECTOR));
		return LIS_ERR_OUT_OF_MEMORY;
	}
	if( solver->precision==LIS_PRECISION_DEFAULT )
	{
		for(i=0;i<worklen;i++)
		{
			err = lis_vector_duplicate(solver->A,&work[i]);
			if( err ) break;
		}
	}
	else
	{
		for(i=0;i<worklen;i++)
		{
			err = lis_vector_duplicateex(LIS_PRECISION_QUAD,solver->A,&work[i]);
			if( err ) break;
		}
	}
	if( i<worklen )
	{
		for(j=0;j<i;j++) lis_vector_destroy(work[j]);
		lis_free(work);
		return err;
	}
	solver->worklen = worklen;
	solver->work    = work;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_jacobi"
int lis_jacobi(LIS_SOLVER solver)
{
	LIS_MATRIX A;
	LIS_PRECON M;
	LIS_VECTOR b,x;
	LIS_VECTOR d,r,t,s;
	LIS_REAL   bnrm2, nrm2, tol;
	int iter,maxiter,n,output,conv;
	double times,ptimes;

	LIS_DEBUG_FUNC_IN;

	A       = solver->A;
	M       = solver->precon;
	b       = solver->b;
	x       = solver->x;
	n       = A->n;
	maxiter = solver->options[LIS_OPTIONS_MAXITER];
	output  = solver->options[LIS_OPTIONS_OUTPUT];
	conv    = solver->options[LIS_OPTIONS_CONV_COND];
	tol     = solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN];
	ptimes  = 0.0;

	r       = solver->work[0];
	t       = solver->work[1];
	s       = solver->work[2];
	d       = solver->work[3];

	lis_vector_nrm2(b,&bnrm2);
	bnrm2   = 1.0 / bnrm2;

	lis_matrix_get_diagonal(A,d);
	lis_vector_reciprocal(d);

	for( iter=1; iter<=maxiter; iter++ )
	{
		/* x += D^{-1}(b - Ax) */
		times = lis_wtime();
		lis_psolve(solver,x,s);
		ptimes += lis_wtime() - times;
		LIS_MATVEC(A,s,t);
/*		LIS_MATVEC(A,x,t);*/
		lis_vector_axpyz(-1,t,b,r);
		lis_vector_nrm2(r,&nrm2);
		lis_vector_pmul(r,d,r);
		lis_vector_axpy(1,r,x);

		/* convergence check */
		nrm2 = nrm2 * bnrm2;

		if( output )
		{
			if( output & LIS_PRINT_MEM ) solver->residual[iter] = nrm2;
			if( output & LIS_PRINT_OUT && A->my_rank==0 ) printf("iter: %5d  residual = %e\n", iter, nrm2);
		}

		if( tol >= nrm2 )
		{
			times = lis_wtime();
			lis_psolve(solver,x,s);
			ptimes += lis_wtime() - times;
			lis_vector_copy(s,x);
			solver->retcode    = LIS_SUCCESS;
			solver->iter       = iter;
			solver->resid      = nrm2;
			solver->ptimes     = ptimes;
			LIS_DEBUG_FUNC_OUT;
			return LIS_SUCCESS;
		}
	}

	lis_psolve(solver,x,s);
	lis_vector_copy(s,x);
	solver->retcode   = LIS_MAXITER;
	solver->iter      = iter;
	solver->resid     = nrm2;
	LIS_DEBUG_FUNC_OUT;
	return LIS_MAXITER;
}

