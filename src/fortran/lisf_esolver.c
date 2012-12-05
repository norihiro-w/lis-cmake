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

/************************************************
 * lis_esolver_create
 * lis_esolver_destroy
 * lis_esolver_set_option
 * lis_esolver_get_option
 * lis_esolve
 ************************************************/
#ifdef USE_FORTRAN

#undef __FUNC__
#define __FUNC__ "lis_esolver_create_f"
void lis_esolver_create_f(LIS_ESOLVER_F *esolver, int *ierr)
{
	LIS_ESOLVER	s;

	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_create(&s);
	if( *ierr )	return;

	*esolver = LIS_P2V(s);
	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_destroy_f"
void lis_esolver_destroy_f(LIS_ESOLVER_F *esolver, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_destroy((LIS_ESOLVER)LIS_V2P(esolver));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolve_f"
void lis_esolve_f(LIS_MATRIX_F *A, LIS_VECTOR_F x, LIS_SCALAR *evalue0, LIS_ESOLVER_F *esolver, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolve((LIS_MATRIX)LIS_V2P(A), (LIS_VECTOR)LIS_V2P(x), evalue0, (LIS_ESOLVER)LIS_V2P(esolver));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_option_f"
void lis_esolver_set_option_f(char *text, LIS_ESOLVER_F *esolver, int *ierr, int len)
{
	char	buf[1024];
	LIS_DEBUG_FUNC_IN;

	strncpy(buf,text,len);
	buf[len] = '\0';
	*ierr = lis_esolver_set_option(buf,(LIS_ESOLVER)LIS_V2P(esolver));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_optionC_f"
void lis_esolver_set_optionC_f(LIS_ESOLVER esolver, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_set_optionC((LIS_ESOLVER)LIS_V2P(esolver));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_iters_f"
void lis_esolver_get_iters_f(LIS_ESOLVER_F *esolver, int *iters, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_iters((LIS_ESOLVER)LIS_V2P(esolver),iters);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_itersex_f"
void lis_esolver_get_itersex_f(LIS_ESOLVER_F *esolver, int *iters, int *iters_double, int *iters_quad, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_itersex((LIS_ESOLVER)LIS_V2P(esolver),iters,iters_double,iters_quad);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_time_f"
void lis_esolver_get_time_f(LIS_ESOLVER_F *esolver, double *times, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_time((LIS_ESOLVER)LIS_V2P(esolver),times);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_timeex_f"
void lis_esolver_get_timeex_f(LIS_ESOLVER_F *esolver, double *times, double *itimes, double *ptimes, double *p_c_times, double *p_i_times, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_timeex((LIS_ESOLVER)LIS_V2P(esolver),times,itimes,ptimes,p_c_times,p_i_times);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_evalues_f"
void lis_esolver_get_evalues_f(LIS_ESOLVER_F *esolver, LIS_VECTOR v, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_evalues((LIS_ESOLVER)LIS_V2P(esolver),(LIS_VECTOR)LIS_V2P(v));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_evectors_f"
void lis_esolver_get_evectors_f(LIS_ESOLVER_F *esolver, LIS_MATRIX_F A, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_evectors((LIS_ESOLVER)LIS_V2P(esolver),(LIS_MATRIX)LIS_V2P(A));
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_residualnorm_f"
void lis_esolver_get_residualnorm_f(LIS_ESOLVER_F *esolver, LIS_REAL *residual, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_residualnorm((LIS_ESOLVER)LIS_V2P(esolver),residual);
	if( *ierr )	return;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_esolver_f"
void lis_esolver_get_esolver_f(LIS_ESOLVER_F *esolver, int *nsol, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_esolver_get_esolver((LIS_ESOLVER)LIS_V2P(esolver),nsol);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_get_esolvername_f"
void lis_get_esolvername_f(int *esolver, char *name, int *ierr, int len)
{
	char	buf[1024];
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_get_esolvername(*esolver, buf);
	strncpy(name,buf,len);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#endif
