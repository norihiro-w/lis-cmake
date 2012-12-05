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


#undef __FUNC__
#define __FUNC__ "lis_precon_create_jacobi"
int lis_precon_create_jacobi(LIS_SOLVER solver, LIS_PRECON precon)
{
	int			err;

	LIS_DEBUG_FUNC_IN;

	if( solver->precision==LIS_PRECISION_DEFAULT )
	{
		err = lis_vector_duplicate(solver->A, &precon->D);
	}
	else
	{
		err = lis_vector_duplicateex(LIS_PRECISION_QUAD,solver->A, &precon->D);
	}
	if( err )
	{
		return err;
	}

	lis_matrix_get_diagonal(solver->A, precon->D);
	lis_vector_reciprocal(precon->D);

	LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_psolve_jacobi"
int lis_psolve_jacobi(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
	int i,n;
	LIS_SCALAR *b,*x,*d;
	LIS_PRECON	precon;
	LIS_QUAD_DECLAR;
	#ifdef USE_QUAD_PRECISION
		LIS_SCALAR *xl;
	#endif

	LIS_DEBUG_FUNC_IN;

	/*
	 *  Mx = b
	 *  M  = D
	 */

	precon = solver->precon;
	n = precon->D->n;
	d = precon->D->value;
	b = B->value;
	x = X->value;
	#ifdef USE_QUAD_PRECISION
		xl = X->value_lo;
	#endif

	#ifdef USE_QUAD_PRECISION
		if( B->precision==LIS_PRECISION_DEFAULT )
		{
	#endif
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0; i<n; i++)
			{
				x[i] = b[i] * d[i];
			}
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			#ifdef _OPENMP
			#ifndef USE_SSE2
				#pragma omp parallel for private(i,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel for private(i,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
			#endif
			#endif
			for(i=0; i<n; i++)
			{
				#ifndef USE_SSE2
					LIS_QUAD_MULD(x[i],xl[i],B->value[i],B->value_lo[i],d[i]);
				#else
					LIS_QUAD_MULD_SSE2(x[i],xl[i],B->value[i],B->value_lo[i],d[i]);
				#endif
				/* x[i] = b[i] * d[i]; */
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_psolvet_jacobi"
int lis_psolvet_jacobi(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
	int i,n;
	LIS_SCALAR *b,*x,*d;
	LIS_PRECON	precon;
	LIS_QUAD_DECLAR;
	#ifdef USE_QUAD_PRECISION
		LIS_SCALAR *xl;
	#endif

	LIS_DEBUG_FUNC_IN;

	/*
	 *  Mx = b
	 *  M  = D
	 */

	precon = solver->precon;
	n = precon->D->n;
	d = precon->D->value;
	b = B->value;
	x = X->value;
	#ifdef USE_QUAD_PRECISION
		xl = X->value_lo;
	#endif

	#ifdef USE_QUAD_PRECISION
		if( B->precision==LIS_PRECISION_DEFAULT )
		{
	#endif
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0; i<n; i++)
			{
				x[i] = b[i] * d[i];
			}
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			#ifdef _OPENMP
			#ifndef USE_SSE2
				#pragma omp parallel for private(i,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel for private(i,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
			#endif
			#endif
			for(i=0; i<n; i++)
			{
				#ifndef USE_SSE2
					LIS_QUAD_MULD(x[i],xl[i],B->value[i],B->value_lo[i],d[i]);
				#else
					LIS_QUAD_MULD_SSE2(x[i],xl[i],B->value[i],B->value_lo[i],d[i]);
				#endif
				/* x[i] = b[i] * d[i]; */
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_precon_create_bjacobi"
int lis_precon_create_bjacobi(LIS_SOLVER solver, LIS_PRECON precon)
{
	int			err;
	LIS_MATRIX	A;

	LIS_DEBUG_FUNC_IN;

	A = solver->A;

	err = lis_matrix_convert_self(solver);
	if( err ) return err;

	if( !A->is_block )
	{
		solver->options[LIS_OPTIONS_PRECON] = LIS_PRECON_TYPE_JACOBI;
		precon->precon_type = LIS_PRECON_TYPE_JACOBI;
		err = lis_precon_create_jacobi(solver,precon);
		return err;
	}

	err = lis_matrix_split(A);
	if( err ) return err;
	err = lis_matrix_diag_duplicate(A->D,&precon->WD);
	if( err ) return err;
	lis_matrix_diag_copy(A->D,precon->WD);
	lis_matrix_diag_inverse(precon->WD);


	LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_psolve_bjacobi"
int lis_psolve_bjacobi(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
	LIS_PRECON	precon;

	LIS_DEBUG_FUNC_IN;

	/*
	 *  Mx = b
	 *  M  = D
	 */

	precon = solver->precon;

	lis_matrix_diag_matvec(precon->WD,B,X);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_psolvet_bjacobi"
int lis_psolvet_bjacobi(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
	LIS_PRECON	precon;

	LIS_DEBUG_FUNC_IN;

	/*
	 *  Mx = b
	 *  M  = D
	 */

	precon = solver->precon;

	lis_matrix_diag_matvect(precon->WD,B,X);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
