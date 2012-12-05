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
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

LIS_MATVEC_FUNC LIS_MATVEC  = lis_matvec;
LIS_MATVEC_FUNC LIS_MATVECT = lis_matvect;

#undef __FUNC__
#define __FUNC__ "lis_matvec"
int lis_matvec(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_SCALAR *x,*y;

	LIS_DEBUG_FUNC_IN;

	if( X->precision==LIS_PRECISION_DEFAULT )
	{
		x = X->value;
		y = Y->value;

		switch( A->matrix_type )
		{
		case LIS_MATRIX_CRS:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_crs(A, x, y);
			break;
		case LIS_MATRIX_BSR:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			if( A->bnr<=4 && A->bnc<=4 )
			{
				lis_matvec_bsr_xxx[A->bnr-1][A->bnc-1](A, x, y);
			}
			else
			{
				lis_matvec_bsr(A, x, y);
			}
			break;
		case LIS_MATRIX_CCS:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_ccs(A, x, y);
			break;
		case LIS_MATRIX_BSC:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_bsc(A, x, y);
			break;
		case LIS_MATRIX_MSR:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_msr(A, x, y);
			break;
		case LIS_MATRIX_ELL:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_ell(A, x, y);
			break;
		case LIS_MATRIX_DIA:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_dia(A, x, y);
			break;
		case LIS_MATRIX_JDS:
			#ifdef USE_MPI
				#ifndef USE_OVERLAP
					LIS_MATVEC_SENDRECV;
				#else
					LIS_MATVEC_REALLOC;
				#endif
			#endif
			lis_matvec_jds(A, x, y);
			break;
		case LIS_MATRIX_VBR:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_vbr(A, x, y);
			break;
		case LIS_MATRIX_DNS:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_dns(A, x, y);
			break;
		case LIS_MATRIX_COO:
			#ifdef USE_MPI
				LIS_MATVEC_SENDRECV;
			#endif
			lis_matvec_coo(A, x, y);
			break;
		default:
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
		}
	}
#ifdef USE_QUAD_PRECISION
	else
	{

		switch( A->matrix_type )
		{
		case LIS_MATRIX_CRS:
			#ifdef USE_MPI
				lis_send_recv_mp(A->commtable,X);
			#endif
			#ifndef USE_FMA2_SSE2
				lis_matvec_crs_mp(A, X, Y);
			#else
				lis_matvec_crs_mp2(A, X, Y);
			#endif
			break;
		case LIS_MATRIX_CCS:
			#ifdef USE_MPI
				lis_send_recv_mp(A->commtable,X);
			#endif
			#ifndef USE_FMA2_SSE2
				lis_matvec_ccs_mp(A, X, Y);
			#else
				lis_matvec_ccs_mp2(A, X, Y);
			#endif
			break;
		default:
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
		}
	}
#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matvect"
int lis_matvect(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y)
{
	LIS_SCALAR *x,*y;

	LIS_DEBUG_FUNC_IN;

	x = X->value;
	y = Y->value;

	if( X->precision==LIS_PRECISION_DEFAULT )
	{
		switch( A->matrix_type )
		{
		case LIS_MATRIX_CRS:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvect_crs(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_BSR:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvect_bsr(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_CCS:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvect_ccs(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_BSC:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvect_bsc(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_MSR:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvect_msr(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_ELL:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvect_ell(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_JDS:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvect_jds(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_DIA:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvect_dia(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_VBR:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvect_vbr(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_DNS:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvect_dns(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		case LIS_MATRIX_COO:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			lis_matvect_coo(A, x, y);
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE;
			#endif
			break;
		default:
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
		}
	}
#ifdef USE_QUAD_PRECISION
	else
	{
		switch( A->matrix_type )
		{
		case LIS_MATRIX_CRS:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			#ifndef USE_FMA2_SSE2
				lis_matvect_crs_mp(A, X, Y);
			#else
				lis_matvect_crs_mp2(A, X, Y);
			#endif
			#ifdef USE_MPI
				lis_reduce_mp(A->commtable,Y);
			#endif
			break;
		case LIS_MATRIX_CCS:
			#ifdef USE_MPI
				LIS_MATVEC_REDUCE0;
			#endif
			#ifndef USE_FMA2_SSE2
				lis_matvect_ccs_mp(A, X, Y);
			#else
				lis_matvect_ccs_mp2(A, X, Y);
			#endif
			#ifdef USE_MPI
				lis_reduce_mp(A->commtable,Y);
			#endif
			break;
		default:
			LIS_SETERR_IMP;
			return LIS_ERR_NOT_IMPLEMENTED;
			break;
		}
	}
#endif
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
