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
 * lis_vector_create_f
 * lis_vector_duplicate_f
 ************************************************/
#ifdef USE_FORTRAN

#undef __FUNC__
#define __FUNC__ "lis_vector_create_f"
void lis_vector_create_f(LIS_Comm_f *comm, LIS_VECTOR_F *vec, int *ierr)
{
	LIS_VECTOR	v;
	LIS_Comm	c_comm;

	LIS_DEBUG_FUNC_IN;

	#ifdef USE_MPI
		if( *comm==lis_comm_world_f )
		{
			c_comm = MPI_COMM_WORLD;
		}
		else
		{
			c_comm = MPI_Comm_f2c(*comm);
		}
	#else
		c_comm = *comm;
	#endif
	*ierr = lis_vector_create(c_comm,&v);
	if( *ierr )	return;

	v->origin = LIS_ORIGIN_1;

	*vec = LIS_P2V(v);
	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_size_f"
void lis_vector_set_size_f(LIS_VECTOR_F *vec, int *local_n, int *global_n, int *ierr)
{
	LIS_VECTOR	v;

	LIS_DEBUG_FUNC_IN;

	v = (LIS_VECTOR)LIS_V2P(vec); 
	v->origin = LIS_ORIGIN_1;

	*ierr = lis_vector_set_size(v,*local_n,*global_n);

	*vec = LIS_P2V(v);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_duplicate_f"
void lis_vector_duplicate_f(LIS_VECTOR_F *vin, LIS_VECTOR_F *vout, int *ierr)
{
	LIS_VECTOR	v;

	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_duplicate((LIS_VECTOR)LIS_V2P(vin),&v);
	*vout = LIS_P2V(v);

	LIS_DEBUG_FUNC_OUT;
	return;
}

/*
#undef __FUNC__
#define __FUNC__ "lis_vector_duplicateM_f"
void lis_vector_duplicateM_f(LIS_MATRIX_F *Amat, LIS_VECTOR_F *vout, int *ierr)
{
	LIS_VECTOR	v;

	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_duplicateM((LIS_MATRIX)LIS_V2P(Amat),&v);
	*vout = LIS_P2V(v);

	LIS_DEBUG_FUNC_OUT;
	return;
}
*/

#undef __FUNC__
#define __FUNC__ "lis_vector_destroy_f"
void lis_vector_destroy_f(LIS_VECTOR_F *vec, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_destroy((LIS_VECTOR)LIS_V2P(vec));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_get_range_f"
void lis_vector_get_range_f(LIS_VECTOR_F *vec, int *is, int *ie, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_get_range((LIS_VECTOR)LIS_V2P(vec),is,ie);
	(*is)++;
	(*ie)++;

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_get_size_f"
void lis_vector_get_size_f(LIS_VECTOR_F *vec, int *local_n, int *global_n, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_get_size((LIS_VECTOR)LIS_V2P(vec),local_n,global_n);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_value_f"
void lis_vector_set_value_f(int *flag, int *i, LIS_SCALAR *value, LIS_VECTOR_F *v, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_set_value(*flag,*i,*value,((LIS_VECTOR)LIS_V2P(v)));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_values_f"
void lis_vector_set_values_f(int *flag, int *count, int *index, LIS_SCALAR *values, LIS_VECTOR_F *v, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_set_values(*flag,*count,index,values,((LIS_VECTOR)LIS_V2P(v)));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_values2_f"
void lis_vector_set_values2_f(int *flag, int *start, int *count, LIS_SCALAR *values, LIS_VECTOR_F *v, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_set_values2(*flag,*start,*count,values,((LIS_VECTOR)LIS_V2P(v)));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_get_value_f"
void lis_vector_get_value_f(LIS_VECTOR_F *v, int *i, LIS_SCALAR *value, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_get_value((LIS_VECTOR)LIS_V2P(v),*i,value);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_get_values_f"
void lis_vector_get_values_f(LIS_VECTOR_F *v, int *start, int *count, LIS_SCALAR *values, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_get_values((LIS_VECTOR)LIS_V2P(v),*start,*count,values);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_scatter_f"
void lis_vector_scatter_f(LIS_SCALAR *values, LIS_VECTOR_F *v, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_scatter(values,(LIS_VECTOR)LIS_V2P(v));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_gather_f"
void lis_vector_gather_f(LIS_VECTOR_F *v, LIS_SCALAR *values, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_gather((LIS_VECTOR)LIS_V2P(v),values);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_print_f"
void lis_vector_print_f(LIS_VECTOR_F *v)
{
	LIS_DEBUG_FUNC_IN;

	lis_vector_print((LIS_VECTOR)LIS_V2P(v));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_all_f"
void lis_vector_set_all_f(LIS_SCALAR *alpha, LIS_VECTOR_F *v, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_set_all(*alpha,((LIS_VECTOR)LIS_V2P(v)));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_axpy_f"
void lis_vector_axpy_f(LIS_SCALAR *alpha, LIS_VECTOR_F *x, LIS_VECTOR_F *y, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_axpy(*alpha,(LIS_VECTOR)LIS_V2P(x),(LIS_VECTOR)LIS_V2P(y));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_xpay_f"
void lis_vector_xpay_f(LIS_VECTOR_F *x, LIS_SCALAR *alpha, LIS_VECTOR_F *y, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_xpay((LIS_VECTOR)LIS_V2P(x),*alpha,(LIS_VECTOR)LIS_V2P(y));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_axpyz_f"
void lis_vector_axpyz_f(LIS_SCALAR *alpha, LIS_VECTOR_F *x, LIS_VECTOR_F *y, LIS_VECTOR_F *z, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_axpyz(*alpha,(LIS_VECTOR)LIS_V2P(x),(LIS_VECTOR)LIS_V2P(y),(LIS_VECTOR)LIS_V2P(z));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_copy_f"
void lis_vector_copy_f(LIS_VECTOR_F *x, LIS_VECTOR_F *y, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_copy((LIS_VECTOR)LIS_V2P(x),(LIS_VECTOR)LIS_V2P(y));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_scale_f"
void lis_vector_scale_f(LIS_SCALAR *alpha, LIS_VECTOR_F *x, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_scale(*alpha,(LIS_VECTOR)LIS_V2P(x));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_pmul_f"
void lis_vector_pmul_f(LIS_VECTOR_F *x, LIS_VECTOR_F *y, LIS_VECTOR_F *z, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_pmul((LIS_VECTOR)LIS_V2P(x),(LIS_VECTOR)LIS_V2P(y),(LIS_VECTOR)LIS_V2P(z));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_pdiv_f"
void lis_vector_pdiv_f(LIS_VECTOR_F *x, LIS_VECTOR_F *y, LIS_VECTOR_F *z, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_pdiv((LIS_VECTOR)LIS_V2P(x),(LIS_VECTOR)LIS_V2P(y),(LIS_VECTOR)LIS_V2P(z));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_abs_f"
void lis_vector_abs_f(LIS_VECTOR_F *x, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_abs((LIS_VECTOR)LIS_V2P(x));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_reciprocal_f"
void lis_vector_reciprocal_f(LIS_VECTOR_F *x, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_reciprocal((LIS_VECTOR)LIS_V2P(x));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_shift_f"
void lis_vector_shift_f(LIS_SCALAR *alpha, LIS_VECTOR_F *x, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_shift(*alpha, (LIS_VECTOR)LIS_V2P(x));

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_dot_f"
void lis_vector_dot_f(LIS_VECTOR_F *x, LIS_VECTOR_F *y, LIS_SCALAR *val, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_dot((LIS_VECTOR)LIS_V2P(x),(LIS_VECTOR)LIS_V2P(y),val);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_nrm2_f"
void lis_vector_nrm2_f(LIS_VECTOR_F *x, LIS_SCALAR *val, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_nrm2((LIS_VECTOR)LIS_V2P(x),val);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_nrm1_f"
void lis_vector_nrm1_f(LIS_VECTOR_F *x, LIS_SCALAR *val, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_nrm1((LIS_VECTOR)LIS_V2P(x),val);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_nrmi_f"
void lis_vector_nrmi_f(LIS_VECTOR_F *x, LIS_SCALAR *val, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_nrmi((LIS_VECTOR)LIS_V2P(x),val);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_sum_f"
void lis_vector_sum_f(LIS_VECTOR_F *x, LIS_SCALAR *val, int *ierr)
{
	LIS_DEBUG_FUNC_IN;

	*ierr = lis_vector_sum((LIS_VECTOR)LIS_V2P(x),val);

	LIS_DEBUG_FUNC_OUT;
	return;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_is_null_f"
void lis_vector_is_null_f(LIS_VECTOR_F *vec, int *ierr)
{
	LIS_VECTOR v;

	LIS_DEBUG_FUNC_IN;

	v = (LIS_VECTOR)LIS_V2P(vec);
	if( !lis_is_malloc(v) )
	{
		*ierr = LIS_TRUE;
	}
	else
	{
		if( v->status==LIS_VECTOR_NULL )
		{
			*ierr = LIS_TRUE;
		}
		else
		{
			*ierr = LIS_FALSE;
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return;
}

#endif
