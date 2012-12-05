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

#define NWORK 3
#undef __FUNC__
#define __FUNC__ "lis_eli_check_params"
int lis_eli_check_params(LIS_ESOLVER esolver)
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
#define __FUNC__ "lis_eli_malloc_work"
int lis_eli_malloc_work(LIS_ESOLVER esolver)
{
	LIS_VECTOR	*work;
	int			i,j,worklen,err,ss;

	LIS_DEBUG_FUNC_IN;

	ss = esolver->options[LIS_EOPTIONS_SUBSPACE];

	worklen = NWORK + ss;
	work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_eli_malloc_work::work" );
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
#define __FUNC__ "lis_eli"
int lis_eli(LIS_ESOLVER esolver)
{
  LIS_MATRIX        A;
  LIS_VECTOR        x;
  LIS_SCALAR        lshift;
  int               ss;
  int               emaxiter,iter0;
  LIS_REAL          tol;
  int               i,j,k;
  int               output, niesolver;
  LIS_REAL          nrm2,dot,resid,resid0,*residual0;
  LIS_VECTOR        *v,r;
  LIS_SCALAR        *t, *tx, tevalue, *tq, *tr, evalue, evalue0;
  LIS_SOLVER        solver;
  LIS_ESOLVER       esolver2;

  ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];
  niesolver = esolver->options[LIS_EOPTIONS_INNER_ESOLVER];

  t = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eli::t");
  tx = (LIS_SCALAR *)lis_malloc(ss*sizeof(LIS_SCALAR), "lis_eli::tx");
  tq = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eli::tq");
  tr = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eli::tr");
  
  A = esolver->A;
  x = esolver->work[0];
  r = esolver->work[1];
  v = &esolver->work[2];
  lis_vector_set_all(0.0,v[0]);
  lis_vector_set_all(1.0,r);
  lis_vector_nrm2(r, &nrm2);
  
  for (i=0;i<ss*ss;i++) t[i] = 0.0;

  j=0;
  while (j<ss-1)
    {
      j = j+1;
      lis_vector_copy(r, v[j]);
      if (j==1) {
	lis_vector_scale(1/nrm2, v[j]);
	lis_matvec(A, v[j], r);
      }
      else {
	lis_vector_scale(1/t[(j-2)*ss+j-1], v[j]);
	lis_matvec(A, v[j], r);
	lis_vector_axpy(-t[(j-2)*ss+j-1], v[j-1], r); 
      }
      lis_vector_dot(v[j], r, &t[(j-1)*ss+j-1]); 
      lis_vector_axpy(-t[(j-1)*ss+j-1], v[j], r); 
      /* diagonalization */
      for (k=1;k<j;k++)
	{ 
	  lis_vector_dot(v[j], v[k], &dot); 
	  lis_vector_axpy(-dot, v[k], v[j]);
	}
      lis_vector_nrm2(r, &t[(j-1)*ss+j]);
      if (t[(j-1)*ss+j]<tol) break;  
      t[j*ss+j-1] = t[(j-1)*ss+j];
    }

  lis_array_qr(ss, t, tq, tr); 

  for (i=0;i<ss;i++)
    {
      esolver->evalue[i] = t[i*ss+i];
    }
  lis_sort_d(0, ss-1, esolver->evalue);

  if( A->my_rank==0 ) 
    {
      printf("\n");
      printf("size of subspace = %d\n", ss);
      for (i=0;i<ss;i++)
	{
	  printf("evalue[%d] = %e\n", i, esolver->evalue[i]); 
	}
      printf("\n");
    }

  lis_esolver_create(&esolver2);
  esolver2->options[LIS_EOPTIONS_ESOLVER] = niesolver;
  esolver2->options[LIS_EOPTIONS_SUBSPACE] = 1;
  esolver2->options[LIS_EOPTIONS_MAXITER] = emaxiter;
  esolver2->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN] = tol; 

  for (i=0;i<ss;i++)
    {
      lis_vector_duplicate(A, &esolver->evector[i]); 
      esolver2->lshift = -(esolver->evalue[i]);
      lis_esolve(A, esolver->evector[i], &evalue, esolver2);
      esolver->evalue[i] = evalue - esolver2->lshift;

      if (i==0) 
	{
	  evalue0 = esolver->evalue[0];
	  iter0 = esolver2->iter;
	  resid0 = esolver2->resid;
	  residual0 = esolver2->residual;
	  esolver->ptimes = esolver2->ptimes;
	  esolver->itimes = esolver2->itimes;
	  esolver->p_c_times = esolver2->p_c_times;
	  esolver->p_i_times = esolver2->p_i_times;
	}

      if (A->my_rank==0) 
	{
	  printf("evalue[%d] = %e\n", i, esolver->evalue[i]);
	  printf("iter = %d\n", esolver2->iter);
	  printf("residual = %e\n\n", esolver2->resid);
	}
    }
  esolver->evalue[0] = evalue0; 
  esolver->iter = iter0;
  esolver->resid = resid0;
  esolver->residual = residual0;

  lis_vector_copy(esolver->evector[esolver->options[LIS_EOPTIONS_MODE]], esolver->x);

  lis_free(t); 
  lis_free(tx);
  lis_free(tq);
  lis_free(tr);

  return LIS_SUCCESS;
}

#ifdef USE_QUAD_PRECISION
#undef __FUNC__
#define __FUNC__ "lis_eli_quad"
int lis_eli_quad(LIS_ESOLVER esolver)
{
  LIS_MATRIX        A;
  LIS_VECTOR        x;
  LIS_SCALAR        lshift;
  int               ss;
  int               emaxiter,iter0;
  LIS_REAL          tol;
  int               i,j,k;
  int               output, niesolver;
  LIS_REAL          nrm2,dot,resid,resid0,*residual0;
  LIS_VECTOR        *v,r;
  LIS_SCALAR        *t, *tx, tevalue, *tq, *tr, evalue, evalue0;
  LIS_SOLVER        solver;
  LIS_ESOLVER       esolver2;

  ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];
  niesolver = esolver->options[LIS_EOPTIONS_INNER_ESOLVER];

  t = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eli::t");
  tx = (LIS_SCALAR *)lis_malloc(ss*sizeof(LIS_SCALAR), "lis_eli::tx");
  tq = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eli::tq");
  tr = (LIS_SCALAR *)lis_malloc(ss*ss*sizeof(LIS_SCALAR), "lis_eli::tr");
  
  A = esolver->A;
  x = esolver->work[0];
  r = esolver->work[1];
  v = &esolver->work[2];
  lis_vector_set_all(0.0,v[0]);
  lis_vector_set_all(1.0,r);
  lis_vector_nrm2(r, &nrm2);
  
  for (i=0;i<ss*ss;i++) t[i] = 0.0;

  j=0;
  while (j<ss-1)
    {
      j = j+1;
      lis_vector_copy(r, v[j]);
      if (j==1) {
	lis_vector_scale(1/nrm2, v[j]);
	lis_matvec(A, v[j], r);
      }
      else {
	lis_vector_scale(1/t[(j-2)*ss+j-1], v[j]);
	lis_matvec(A, v[j], r);
	lis_vector_axpy(-t[(j-2)*ss+j-1], v[j-1], r); 
      }
      lis_vector_dot(v[j], r, &t[(j-1)*ss+j-1]); 
      lis_vector_axpy(-t[(j-1)*ss+j-1], v[j], r); 
      /* diagonalization */
      for (k=1;k<j;k++)
	{ 
	  lis_vector_dot(v[j], v[k], &dot); 
	  lis_vector_axpy(-dot, v[k], v[j]);
	}
      lis_vector_nrm2(r, &t[(j-1)*ss+j]);
      if (t[(j-1)*ss+j]<tol) break;  
      t[j*ss+j-1] = t[(j-1)*ss+j];
    }

  lis_array_qr(ss, t, tq, tr); 

  for (i=0;i<ss;i++)
    {
      esolver->evalue[i] = t[i*ss+i];
    }
  lis_sort_d(0, ss-1, esolver->evalue);

  if( A->my_rank==0 ) 
    {
      printf("\n");
      printf("size of subspace = %d\n", ss);
      for (i=0;i<ss;i++)
	{
	  printf("evalue[%d] = %e\n", i, esolver->evalue[i]); 
	}
      printf("\n");
    }

  lis_esolver_create(&esolver2);
  esolver2->options[LIS_EOPTIONS_ESOLVER] = niesolver;
  esolver2->options[LIS_EOPTIONS_SUBSPACE] = 1;
  esolver2->options[LIS_EOPTIONS_MAXITER] = emaxiter;
  esolver2->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN] = tol; 
  esolver2->eprecision = LIS_PRECISION_QUAD;

  for (i=0;i<ss;i++)
    {
      lis_vector_duplicate(A, &esolver->evector[i]); 
      esolver2->lshift = -(esolver->evalue[i]);
      lis_esolve(A, esolver->evector[i], &evalue, esolver2);
      esolver->evalue[i] = evalue - esolver2->lshift;

      if (i==0) 
	{
	  evalue0 = esolver->evalue[0];
	  iter0 = esolver2->iter;
	  resid0 = esolver2->resid;
	  residual0 = esolver2->residual;
	  esolver->ptimes = esolver2->ptimes;
	  esolver->itimes = esolver2->itimes;
	  esolver->p_c_times = esolver2->p_c_times;
	  esolver->p_i_times = esolver2->p_i_times;
	}

      if (A->my_rank==0) 
	{
	  printf("evalue[%d] = %e\n", i, esolver->evalue[i]);
	  printf("iter = %d\n", esolver2->iter);
	  printf("residual = %e\n\n", esolver2->resid);
	}
    }
  esolver->evalue[0] = evalue0; 
  esolver->iter = iter0;
  esolver->resid = resid0;
  esolver->residual = residual0;

  lis_vector_copy(esolver->evector[esolver->options[LIS_EOPTIONS_MODE]], esolver->x);

  lis_free(t); 
  lis_free(tx);
  lis_free(tq);
  lis_free(tr);

  return LIS_SUCCESS;
}
#endif

