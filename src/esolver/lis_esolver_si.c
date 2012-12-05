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

#define NWORK 4
#undef __FUNC__
#define __FUNC__ "lis_esi_check_params"
int lis_esi_check_params(LIS_ESOLVER esolver)
{
        int ss;

	LIS_DEBUG_FUNC_IN;

	ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
	if( ss<0 )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_SUBSPACE(=%d) is less than 0\n",ss);
		return LIS_ERR_ILL_ARG;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esi_malloc_work"
int lis_esi_malloc_work(LIS_ESOLVER esolver)
{
	LIS_VECTOR	*work;
	int			i,j,worklen,err,ss;

	LIS_DEBUG_FUNC_IN;

	ss = esolver->options[LIS_EOPTIONS_SUBSPACE];

	worklen = NWORK + ss;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_esi_malloc_work::work" );
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
#define __FUNC__ "lis_esi"
int lis_esi(LIS_ESOLVER esolver)
{
  LIS_MATRIX        A;
  LIS_VECTOR        x, Ax;
  LIS_SCALAR        xAx, xx, mu, lshift;
  int               ss;
  int               emaxiter;
  LIS_REAL          tol;
  int               i,j,k;
  LIS_SCALAR        evalue,dotvr;
  int               iter,giter,output,niesolver;
  int               nprocs,my_rank;
  LIS_REAL          nrm2,dot,resid,resid0;
  LIS_VECTOR        *v,r,q;
  LIS_SOLVER        solver;
  LIS_PRECON        precon;
  double	    times,itimes,ptimes,p_c_times,p_i_times;
  int		    err;

  LIS_DEBUG_FUNC_IN;

  A = esolver->A;
  x = esolver->x;

  ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  lshift = esolver->lshift;
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];
  niesolver = esolver->options[LIS_EOPTIONS_INNER_ESOLVER];

  r = esolver->work[0];
  q = esolver->work[1];
  v = &esolver->work[2];
  Ax = esolver->work[3];
  lis_vector_set_all(1.0,r);
  lis_vector_nrm2(r, &nrm2);
  lis_vector_scale(1/nrm2,r);

  switch ( niesolver )
    {
    case LIS_ESOLVER_II:
      lis_solver_create(&solver);
      lis_solver_set_option("-i bicg -p ilu",solver);
      lis_solver_set_optionC(solver);
      if( A->my_rank==0 ) printf("\nlocal shift for inverse iteration = %e\n", lshift);
      if (lshift != 0) lis_matrix_shift_diagonal(A, lshift);
      break;
    case LIS_ESOLVER_AII:
      lis_solver_create(&solver);
      lis_solver_set_option("-i bicg -p ilu",solver);
      lis_solver_set_optionC(solver);
      if( A->my_rank==0 ) printf("\nlocal shift for approximate inverse iteration = %e\n", lshift);
      if (lshift != 0) lis_matrix_shift_diagonal(A, lshift);
      lis_vector_set_all(1.0,q);
      lis_solve(A, q, x, solver);
      lis_precon_create(solver, &precon);
      solver->precon = precon;
      break;
    case LIS_ESOLVER_RQI:
      lis_solver_create(&solver);
      lis_solver_set_option("-p ilu -maxiter 10",solver);
      lis_solver_set_optionC(solver);
      if( A->my_rank==0 ) printf("\nlocal shift for Rayleigh quotient iteration = %e\n", lshift);
      if (lshift != 0) lis_matrix_shift_diagonal(A, lshift);
      break;
    }

  giter=0;
  j=0;
  while (j<ss)
    {
      lis_vector_duplicate(A,&esolver->evector[j]); 
      j = j+1;
      lis_vector_copy(r, v[j]);

      if (niesolver==LIS_ESOLVER_II || niesolver==LIS_ESOLVER_RQI)
	{
	  /* create preconditioner */
	  solver->A = A;
	  err = lis_precon_create(solver, &precon);
	  if( err )
	    {
	      lis_solver_work_destroy(solver);
	      solver->retcode = err;
	      return err;
	    }
	}

      if (niesolver==LIS_ESOLVER_RQI)
	{
	  lis_vector_nrm2(x, &nrm2);
	  lis_vector_scale(1/nrm2, x);
	  lis_matvec(A, x, Ax);
	  lis_vector_dot(x, Ax, &xAx);
	  lis_vector_dot(x, x, &xx);
	  mu = xAx / xx;
	}

      iter = 0;
      while (iter<emaxiter)
	{
	  /* diagonalization */
	  iter = iter+1;
	  giter = giter+1;
	  for (k=1;k<j;k++)
	    { 
	      lis_vector_dot(v[j], v[k], &dot); 
	      lis_vector_axpy(-dot, v[k], v[j]);
	    }

	  switch( niesolver )
	    {
	    case LIS_ESOLVER_PI:
	      lis_matvec(A,v[j],r); 
	      break;
	    case LIS_ESOLVER_II:
	      lis_solve_kernel(A, v[j], r, solver, precon);
	      break;
	    case LIS_ESOLVER_AII:
	      lis_psolve(solver, v[j], r); 
	      break;
	    case LIS_ESOLVER_RQI:
	      lis_vector_nrm2(v[j], &nrm2);
	      lis_vector_scale(1/nrm2, v[j]);
	      lis_matrix_shift_diagonal(A, -mu);
	      lis_solve_kernel(A, v[j], r, solver, precon);
	      lis_matrix_shift_diagonal(A, mu);
	      break;
	    }

	  if ( j==1 && ( niesolver==LIS_ESOLVER_II || niesolver==LIS_ESOLVER_AII || niesolver==LIS_ESOLVER_RQI ))
	    {
	      lis_solver_get_timeex(solver,&times,&itimes,&ptimes,&p_c_times,&p_i_times);
	      esolver->ptimes += solver->ptimes;
	      esolver->itimes += solver->itimes;
	      esolver->p_c_times += solver->p_c_times;
	      esolver->p_i_times += solver->p_i_times;
	    }

	  lis_vector_nrm2(r, &nrm2);
	  lis_vector_dot(v[j],r,&dotvr);
	  mu = mu + 1/dotvr;
	  lis_vector_axpyz(-dotvr,v[j],r,q);
	  lis_vector_nrm2(q, &resid);
	  resid = fabs(resid / dotvr);
	  lis_vector_scale(1/nrm2,r);
	  lis_vector_copy(r, v[j]);
	  if ( j==1 ) 
	    {
	      if( output & LIS_PRINT_MEM ) esolver->residual[iter] = resid; 
	      if( output & LIS_PRINT_OUT ) printf("iter: %5d  residual = %e\n", iter, resid);
	      esolver->iter = iter;
	      esolver->resid = resid;
	    }
	  if (tol>resid) break;
	}

      if (niesolver==LIS_ESOLVER_II || niesolver==LIS_ESOLVER_RQI)
	{
	  lis_precon_destroy(precon);
	}

      switch ( niesolver )
	{
	case LIS_ESOLVER_PI:
	  esolver->evalue[j-1] = dotvr;
	  break;
	case LIS_ESOLVER_II:
	  esolver->evalue[j-1] = 1/dotvr;
	  break;
	case LIS_ESOLVER_AII:
	  esolver->evalue[j-1] = 1/dotvr;
	  break;
	case LIS_ESOLVER_RQI:
	  esolver->evalue[j-1] = mu;
	  break;
	}
      lis_vector_copy(v[j], esolver->evector[j-1]);  

      if (A->my_rank==0 && ss>1)
	{
	  printf("evalue[%d] = %e\n", j-1, esolver->evalue[j-1]);
	  printf("residual = %e\n", resid);
	  printf("iter = %d\n\n", iter);
	}
    }
  
  lis_vector_copy(esolver->evector[esolver->options[LIS_EOPTIONS_MODE]], esolver->x);

  switch ( niesolver )
    {
    case LIS_ESOLVER_II:
      if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
      lis_solver_destroy(solver);
      break;
    case LIS_ESOLVER_AII:
      if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
      lis_precon_destroy(precon);
      lis_solver_destroy(solver);
      break;
    case LIS_ESOLVER_RQI:
      if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
      lis_solver_destroy(solver);
      break;
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#ifdef USE_QUAD_PRECISION
#undef __FUNC__
#define __FUNC__ "lis_esi_quad"
int lis_esi_quad(LIS_ESOLVER esolver)
{
  LIS_MATRIX        A;
  LIS_VECTOR        x, Ax;
  LIS_SCALAR        xAx, xx, mu, lshift;
  int               ss;
  int               emaxiter;
  LIS_REAL          tol;
  int               i,j,k;
  LIS_SCALAR        evalue,dotvr;
  int               iter,giter,output,niesolver;
  int               nprocs,my_rank;
  LIS_REAL          nrm2,dot,resid,resid0;
  LIS_QUAD_PTR      qdot_vv, qdot_vr;
  LIS_VECTOR        *v,r,q;
  LIS_SOLVER        solver;
  LIS_PRECON        precon;
  double	    times,itimes,ptimes,p_c_times,p_i_times;
  int		    err;

  LIS_DEBUG_FUNC_IN;

  A = esolver->A;
  x = esolver->x;

  ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  lshift = esolver->lshift;
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];
  niesolver = esolver->options[LIS_EOPTIONS_INNER_ESOLVER];

  r = esolver->work[0];
  q = esolver->work[1];
  v = &esolver->work[2];
  Ax = esolver->work[3];

  LIS_QUAD_SCALAR_MALLOC(qdot_vv,0,1);
  LIS_QUAD_SCALAR_MALLOC(qdot_vr,1,1);

  lis_vector_set_all(1.0,r);
  lis_vector_nrm2(r, &nrm2);
  lis_vector_scale(1/nrm2,r);
	  
  switch ( niesolver )
    {
    case LIS_ESOLVER_II:
      lis_solver_create(&solver);
      lis_solver_set_option("-i bicg -p ilu -precision quad",solver);
      lis_solver_set_optionC(solver);
      if( A->my_rank==0 ) printf("\nlocal shift for inverse iteration = %e\n", lshift);
      if (lshift != 0) lis_matrix_shift_diagonal(A, lshift);
      break;
    case LIS_ESOLVER_AII:
      lis_solver_create(&solver);
      lis_solver_set_option("-i bicg -p ilu -precision quad",solver);
      lis_solver_set_optionC(solver);
      if( A->my_rank==0 ) printf("\nlocal shift for approximate inverse iteration = %e\n", lshift);
      if (lshift != 0) lis_matrix_shift_diagonal(A, lshift);
      lis_vector_set_all(1.0,q);
      lis_solve(A, q, x, solver);
      lis_precon_create(solver, &precon);
      solver->precon = precon;
      break;
    case LIS_ESOLVER_RQI:
      lis_solver_create(&solver);
      lis_solver_set_option("-p ilu -precision quad -maxiter 10",solver);
      lis_solver_set_optionC(solver);
      if( A->my_rank==0 ) printf("\nlocal shift for Rayleigh quotient iteration = %e\n", lshift);
      if (lshift != 0) lis_matrix_shift_diagonal(A, lshift);
      break;
    }

  giter=0;
  j=0;
  while (j<ss)
    {
      lis_vector_duplicate(A,&esolver->evector[j]); 
      j = j+1;
      lis_vector_copy(r, v[j]);

      if (niesolver==LIS_ESOLVER_II || niesolver==LIS_ESOLVER_RQI)
	{
	  /* create preconditioner */
	  solver->A = A;
	  err = lis_precon_create(solver, &precon);
	  if( err )
	    {
	      lis_solver_work_destroy(solver);
	      solver->retcode = err;
	      return err;
	    }
	}

      if (niesolver==LIS_ESOLVER_RQI)
	{
	  lis_vector_nrm2(x, &nrm2);
	  lis_vector_scale(1/nrm2, x);
	  lis_matvec(A, x, Ax);
	  lis_vector_dot(x, Ax, &xAx);
	  lis_vector_dot(x, x, &xx);
	  mu = xAx / xx;
	}

      iter = 0;
      while (iter<emaxiter)
	{
	  /* diagonalization */
	  iter = iter+1;
	  giter = giter+1;
	  for (k=1;k<j;k++)
	    { 
	      lis_vector_dotex_mmm(v[j], v[k], &qdot_vv);
	      lis_quad_minus((LIS_QUAD *)qdot_vv.hi);
	      lis_vector_axpyex_mmm(qdot_vv,v[k],v[j]);
	    }

	  switch( niesolver )
	    {
	    case LIS_ESOLVER_PI:
	      lis_matvec(A,v[j],r); 
	      break;
	    case LIS_ESOLVER_II:
	      lis_solve_kernel(A, v[j], r, solver, precon);
	      break;
	    case LIS_ESOLVER_AII:
	      lis_psolve(solver, v[j], r); 
	      break;
	    case LIS_ESOLVER_RQI:
	      lis_vector_nrm2(v[j], &nrm2);
	      lis_vector_scale(1/nrm2, v[j]);
	      lis_matrix_shift_diagonal(A, -mu);
	      lis_solve_kernel(A, v[j], r, solver, precon);
	      lis_matrix_shift_diagonal(A, mu);
	      break;
	    }

	  if ( j==1 && ( niesolver==LIS_ESOLVER_II || niesolver==LIS_ESOLVER_AII || niesolver==LIS_ESOLVER_RQI ))
	    {
	      lis_solver_get_timeex(solver,&times,&itimes,&ptimes,&p_c_times,&p_i_times);
	      esolver->ptimes += solver->ptimes;
	      esolver->itimes += solver->itimes;
	      esolver->p_c_times += solver->p_c_times;
	      esolver->p_i_times += solver->p_i_times;
	    }

	  lis_vector_nrm2(r, &nrm2);
	  lis_vector_dotex_mmm(v[j], r, &qdot_vr);
	  lis_quad_minus((LIS_QUAD *)qdot_vr.hi);
	  lis_vector_axpyzex_mmmm(qdot_vr,v[j],r,q);
	  lis_quad_minus((LIS_QUAD *)qdot_vr.hi);	  
	  dotvr = qdot_vr.hi[0];
	  mu = mu + 1/dotvr;

	  lis_vector_nrm2(q, &resid);
	  resid = fabs(resid / dotvr);
	  lis_vector_scale(1/nrm2,r);
	  lis_vector_copy(r, v[j]);
	  if ( j==1 ) 
	    {
	      if( output & LIS_PRINT_MEM ) esolver->residual[iter] = resid; 
	      if( output & LIS_PRINT_OUT ) printf("iter: %5d  residual = %e\n", iter, resid);
	      esolver->iter = iter;
	      esolver->resid = resid;
	    }
	  if (tol>resid) break;
	}

      if (niesolver==LIS_ESOLVER_II || niesolver==LIS_ESOLVER_RQI)
	{
	  lis_precon_destroy(precon);
	}

      switch ( niesolver )
	{
	case LIS_ESOLVER_PI:
	  esolver->evalue[j-1] = dotvr;
	  break;
	case LIS_ESOLVER_II:
	  esolver->evalue[j-1] = 1/dotvr;
	  break;
	case LIS_ESOLVER_AII:
	  esolver->evalue[j-1] = 1/dotvr;
	  break;
	case LIS_ESOLVER_RQI:
	  esolver->evalue[j-1] = mu;
	  break;
	}
      lis_vector_copy(v[j], esolver->evector[j-1]);  

      if (A->my_rank==0 && ss>1)
	{
	  printf("evalue[%d] = %e\n", j-1, esolver->evalue[j-1]);
	  printf("residual = %e\n", resid);
	  printf("iter = %d\n\n", iter);
	}
    }
  
  lis_vector_copy(esolver->evector[esolver->options[LIS_EOPTIONS_MODE]], esolver->x);

  switch ( niesolver )
    {
    case LIS_ESOLVER_II:
      if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
      lis_solver_destroy(solver);
      break;
    case LIS_ESOLVER_AII:
      if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
      lis_precon_destroy(precon);
      lis_solver_destroy(solver);
      break;
    case LIS_ESOLVER_RQI:
      if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
      lis_solver_destroy(solver);
      break;
    }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}
#endif

