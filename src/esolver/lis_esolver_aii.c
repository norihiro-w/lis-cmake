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
#include <math.h>
#include <string.h>
#include <stdarg.h>
#ifdef USE_SSE2
	#include <emmintrin.h>
#endif
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

#define NWORK 2
#undef __FUNC__
#define __FUNC__ "lis_eaii_check_params"
int lis_eaii_check_params(LIS_ESOLVER esolver)
{
	LIS_DEBUG_FUNC_IN;
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_eaii_malloc_work"
int lis_eaii_malloc_work(LIS_ESOLVER esolver)
{
	LIS_VECTOR	*work;
	int			i,j,worklen,err;

	LIS_DEBUG_FUNC_IN;

	worklen = NWORK;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_eaii_malloc_work::work" );
	if( work==NULL )
	{
		LIS_SETERR_MEM(worklen*sizeof(LIS_VECTOR));
		return LIS_ERR_OUT_OF_MEMORY;
	}
	if( esolver->eprecision==LIS_PRECISION_DEFAULT )
	{
		for(i=0;i<worklen;i++)
		{
			err = lis_vector_duplicate(esolver->A,&work[i]);
			if( err ) break;
		}
	}
	else
	{
		for(i=0;i<worklen;i++)
		{
			err = lis_vector_duplicateex(LIS_PRECISION_QUAD,esolver->A,&work[i]);
			if( err ) break;
		}
	}
	if( i<worklen )
	{
		for(j=0;j<i;j++) lis_vector_destroy(work[j]);
		lis_free(work);
		return err;
	}
	esolver->worklen = worklen;
	esolver->work    = work;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_eaii"
int lis_eaii(LIS_ESOLVER esolver)
{
  LIS_MATRIX        A;
  LIS_VECTOR        x;
  LIS_SCALAR        evalue, ievalue;
  LIS_SCALAR        lshift;
  int               emaxiter;
  LIS_REAL          tol;
  int               iter,iter2,i,n,gn,is,ie,output;
  int               nprocs,my_rank;
  LIS_REAL          nrm2,resid;
  LIS_VECTOR        z,q;
  LIS_MATRIX        A0;
  LIS_SOLVER        solver;
  LIS_PRECON        precon;
  double	    times,itimes,ptimes,p_c_times,p_i_times;

  LIS_DEBUG_FUNC_IN;

  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  lshift = esolver->lshift;
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];

  A = esolver->A;
  x = esolver->x;
  if (esolver->options[LIS_EOPTIONS_INITGUESS_ONES] ) 
    {
      lis_vector_set_all(1.0,x);
    }
  evalue = 1.0;
  z = esolver->work[0];
  q = esolver->work[1];

  iter=0;
  ievalue = 1/(evalue);
  if( A->my_rank==0 ) printf("\nlocal shift for preconditioned power iteration = %e\n", lshift);
  if (lshift != 0) lis_matrix_shift_diagonal(A, lshift);
  lis_solver_create(&solver);
  lis_solver_set_option("-i bicg -p ilu",solver);
  lis_solver_set_optionC(solver);

  lis_vector_set_all(1.0,q);
  lis_solve(A, q, x, solver);
  lis_precon_create(solver, &precon);
  solver->precon = precon;

  while (iter<emaxiter)
    {
      iter = iter+1;
      lis_vector_nrm2(x, &nrm2);
      lis_vector_scale(1/nrm2, x);
      lis_psolve(solver, x, z); 
      lis_vector_dot(x, z, &ievalue); 
      lis_vector_axpyz(-ievalue,x,z,q); 
      lis_vector_nrm2(q, &resid); 
      resid = fabs(resid/ievalue);
      if( output )
	{
	  if( output & LIS_EPRINT_MEM ) esolver->residual[iter] = resid;
	  if( output & LIS_EPRINT_OUT && A->my_rank==0 ) printf("iter: %5d  residual = %e\n", iter, resid);
	}
      lis_vector_copy(z,x);
      lis_solver_get_timeex(solver,&times,&itimes,&ptimes,&p_c_times,&p_i_times);
      esolver->ptimes += solver->ptimes;
      esolver->itimes += solver->itimes;
      esolver->p_c_times += solver->p_c_times;
      esolver->p_i_times += solver->p_i_times;
      if( tol >= resid ) 
	{
	  esolver->retcode    = LIS_SUCCESS;
	  esolver->iter       = iter;
	  esolver->resid      = resid;
	  esolver->evalue[0] = 1/ievalue;
	  if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
	  LIS_DEBUG_FUNC_OUT;
	  return LIS_SUCCESS;
	}
    }
  esolver->retcode    = LIS_MAXITER;
  esolver->iter       = iter;
  esolver->resid      = resid;
  esolver->evalue[0] = 1/ievalue;
  if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
  lis_precon_destroy(precon);
  lis_solver_destroy(solver);
  LIS_DEBUG_FUNC_OUT;
  return LIS_MAXITER;
}
