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
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_vector_init
 * lis_vector_create
 * lis_vector_duplicate
 * lis_vector_destroy
 * lis_vector_set_value
 * lis_vector_set_values
 * lis_vector_set_values2
 * lis_vector_get_value
 * lis_vector_get_values
 * lis_vector_get_range
 * lis_vector_get_size
 * lis_vector_scatter
 * lis_vector_gather
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_vector_init"
int lis_vector_init(LIS_VECTOR *vec)
{
	LIS_DEBUG_FUNC_IN;

	memset(*vec,0,sizeof(struct LIS_VECTOR_STRUCT));
	(*vec)->status = LIS_VECTOR_NULL;
	(*vec)->is_destroy = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_check"
int lis_vector_check(LIS_VECTOR v, int level)
{
	LIS_DEBUG_FUNC_IN;

	switch( level )
	{
	case LIS_VECTOR_CHECK_NULL:
		if( !lis_is_malloc(v) )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"vector v is undefined\n");
			return LIS_ERR_ILL_ARG;
		}
		break;
	default:
		if( !lis_is_malloc(v) )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"vector v is undefined\n");
			return LIS_ERR_ILL_ARG;
		}
		if( v->status<=LIS_VECTOR_ASSEMBLING )
		{
			LIS_SETERR(LIS_ERR_ILL_ARG,"vector v is assembling\n");
			return LIS_ERR_ILL_ARG;
		}
		break;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}



#undef __FUNC__
#define __FUNC__ "lis_vector_create"
int lis_vector_create(LIS_Comm comm, LIS_VECTOR *vec)
{
	int	err;

	LIS_DEBUG_FUNC_IN;

	err = lis_vector_createex(LIS_PRECISION_DEFAULT,comm,vec);
	if( err ) return err;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_createex"
int lis_vector_createex(int precision, LIS_Comm comm, LIS_VECTOR *vec)
{

	LIS_DEBUG_FUNC_IN;

	*vec = NULL;

	*vec = (LIS_VECTOR)lis_malloc( sizeof(struct LIS_VECTOR_STRUCT),"lis_vector_createex::vec" );
	if( NULL==*vec )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_VECTOR_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	lis_vector_init(vec);
	
	(*vec)->status      = LIS_VECTOR_NULL;
	(*vec)->precision   = precision;
	(*vec)->comm        = comm;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_vector_set_size"
int lis_vector_set_size(LIS_VECTOR vec, int local_n, int global_n)
{
	int nprocs,my_rank;
	int is,ie;
	int i,err,precision;
	int *ranges;

	LIS_DEBUG_FUNC_IN;

	if( global_n>0 && local_n>global_n )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"local n(=%d) is larger than global n(=%d)\n",local_n,global_n);
		return LIS_ERR_ILL_ARG;
	}
	if( local_n<0 || global_n<0 )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"local n(=%d) or global n(=%d) are less than 0\n",local_n,global_n);
		return LIS_ERR_ILL_ARG;
	}
	if( local_n==0 && global_n==0 )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"local n(=%d) and global n(=%d) are 0\n",local_n,global_n);
		return LIS_ERR_ILL_ARG;
	}


	err = lis_ranges_create(vec->comm,&local_n,&global_n,&ranges,&is,&ie,&nprocs,&my_rank);
	if( err )
	{
		return err;
	}
	vec->ranges      = ranges;

	precision = vec->precision;
	if( !precision )
	{
		vec->value = (LIS_SCALAR *)lis_malloc( local_n*sizeof(LIS_SCALAR),"lis_vector_set_size::vec->value" );
		if( NULL==vec->value )
		{
			LIS_SETERR_MEM(local_n*sizeof(LIS_SCALAR));
			return LIS_OUT_OF_MEMORY;
		}
		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0;i<local_n;i++)
		{
			vec->value[i] = 0.0;
		}
	}
	else
	{
		vec->value = (LIS_SCALAR *)lis_malloc( (2*local_n + local_n%2)*sizeof(LIS_SCALAR),"lis_vector_set_size::vec->value" );
		if( NULL==vec->value )
		{
			LIS_SETERR_MEM((2*local_n+local_n%2)*sizeof(LIS_SCALAR));
			return LIS_OUT_OF_MEMORY;
		}
		vec->value_lo = vec->value + local_n + local_n%2;
		vec->work = (LIS_SCALAR *)lis_malloc( 32*sizeof(LIS_SCALAR),"lis_vector_set_size::vec->work" );
		if( NULL==vec->work )
		{
			LIS_SETERR_MEM(32*sizeof(LIS_SCALAR));
			return LIS_OUT_OF_MEMORY;
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#endif
		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0;i<local_n;i++)
		{
			vec->value[i]    = 0.0;
			vec->value_lo[i] = 0.0;
		}
	}
	
	vec->is_copy     = LIS_TRUE;
	vec->status      = LIS_VECTOR_ASSEMBLED;
	vec->n           = local_n;
	vec->gn          = global_n;
	vec->np          = local_n;
	vec->my_rank     = my_rank;
	vec->nprocs      = nprocs;
	vec->is          = is;
	vec->ie          = ie;
	vec->origin      = LIS_ORIGIN_0;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_reuse"
int lis_vector_reuse(LIS_VECTOR *vec)
{
	int			err,np,precision;

	LIS_DEBUG_FUNC_IN;

	err = lis_vector_check(*vec,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	np = (*vec)->np;
	if( (*vec)->status==LIS_VECTOR_NULL )
	{
		precision = ((LIS_VECTOR)*vec)->precision;
		if( !precision )
		{
			(*vec)->value = (LIS_SCALAR *)lis_malloc( np*sizeof(LIS_SCALAR),"lis_vector_reuse::vec->value" );
			if( NULL==(*vec)->value )
			{
				LIS_SETERR_MEM(np*sizeof(LIS_SCALAR));
				return LIS_OUT_OF_MEMORY;
			}
			(*vec)->is_copy = LIS_TRUE;
		}
		else
		{
			(*vec)->value = (LIS_SCALAR *)lis_malloc( (2*np+np%2)*sizeof(LIS_SCALAR),"lis_vector_reuse::vec->value" );
			if( NULL==(*vec)->value )
			{
				LIS_SETERR_MEM((2*np+np%2)*sizeof(LIS_SCALAR));
				return LIS_OUT_OF_MEMORY;
			}
			(*vec)->is_copy = LIS_TRUE;
			(*vec)->value_lo = (*vec)->value + np + np%2;
			(*vec)->work = (LIS_SCALAR *)lis_malloc( 32*sizeof(LIS_SCALAR),"lis_vector_reuse::vec->work" );
			if( NULL==(*vec)->work )
			{
				LIS_SETERR_MEM(32*sizeof(LIS_SCALAR));
				lis_vector_destroy(*vec);
				*vec = NULL;
				return LIS_OUT_OF_MEMORY;
			}
		}
	}

	(*vec)->status = LIS_VECTOR_ASSEMBLED;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_unset"
int lis_vector_unset(LIS_VECTOR vec)
{
	int			err;

	LIS_DEBUG_FUNC_IN;

	err = lis_vector_check(vec,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	if( vec->is_copy ) lis_free(vec->value);
	vec->value  = NULL;
	vec->status = LIS_VECTOR_NULL;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set"
int lis_vector_set(LIS_VECTOR vec, LIS_SCALAR *value)
{
	int			err;
	int			n,np;

	LIS_DEBUG_FUNC_IN;

	err = lis_vector_check(vec,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	n  = vec->n;
	np = vec->np;
	if( vec->is_destroy ) lis_free(vec->value);
	vec->value   = value;
	vec->is_copy = LIS_FALSE;

	vec->status = LIS_VECTOR_ASSEMBLING;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_destroy"
int lis_vector_destroy(LIS_VECTOR vec)
{
	LIS_DEBUG_FUNC_IN;
	if( lis_is_malloc(vec) )
	{
		if( vec->value && vec->is_destroy ) lis_free( vec->value );
		if( vec->work ) lis_free( vec->work );
		if( vec->ranges ) lis_free( vec->ranges );
		if( vec ) lis_free(vec);
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_vector_duplicate"
int lis_vector_duplicate(void *vin, LIS_VECTOR *vout)
{
	int precision,err;

	LIS_DEBUG_FUNC_IN;

	precision = LIS_PRECISION_DEFAULT;
	if( ((LIS_VECTOR)vin)->label==LIS_LABEL_VECTOR)
	{
		precision = ((LIS_VECTOR)vin)->precision;
	}
	else if( ((LIS_VECTOR)vin)->label!=LIS_LABEL_MATRIX)
	{
		LIS_SETERR(LIS_ERR_ILL_ARG, "First argument is not LIS_VECTOR or LIS_MATRIX\n");
		return LIS_ERR_ILL_ARG;
	}
	err = lis_vector_duplicateex(precision,vin,vout);

	LIS_DEBUG_FUNC_OUT;
	return err;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_duplicateex"
int lis_vector_duplicateex(int precision, void *A, LIS_VECTOR *vout)
{
	int n,np,pad;
	int nprocs;
	int i;
	#ifdef USE_MPI
		int *ranges;
	#endif
	LIS_SCALAR *value;

	LIS_DEBUG_FUNC_IN;

	if( ((LIS_VECTOR)A)->label!=LIS_LABEL_VECTOR && ((LIS_VECTOR)A)->label!=LIS_LABEL_MATRIX)
	{
		LIS_SETERR(LIS_ERR_ILL_ARG, "Second argument is not LIS_VECTOR or LIS_MATRIX\n");
		return LIS_ERR_ILL_ARG;
	}
	nprocs = ((LIS_VECTOR)A)->nprocs;
	n      = ((LIS_VECTOR)A)->n;
	np     = ((LIS_VECTOR)A)->np;
	pad    = ((LIS_VECTOR)A)->pad;
	*vout  = NULL;
	*vout  = (LIS_VECTOR)lis_malloc( sizeof(struct LIS_VECTOR_STRUCT),"lis_vector_duplicateex::vout" );
	if( NULL==*vout )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_VECTOR_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	lis_vector_init(vout);


	if( !precision )
	{
		value = (LIS_SCALAR *)lis_malloc( (np+pad)*sizeof(LIS_SCALAR),"lis_vector_duplicateex::value" );
		if( NULL==value )
		{
			LIS_SETERR_MEM((np+pad)*sizeof(LIS_SCALAR));
			lis_vector_destroy(*vout);
			*vout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		(*vout)->value = value;
		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0;i<np+pad;i++)
		{
			(*vout)->value[i] = 0.0;
		}
	}
	else
	{
		value = (LIS_SCALAR *)lis_malloc( (2*(np+pad) + (np+pad)%2)*sizeof(LIS_SCALAR),"lis_vector_duplicateex::value" );
		if( NULL==value )
		{
			LIS_SETERR_MEM((2*(np+pad) + (np+pad)%2)*sizeof(LIS_SCALAR));
			lis_vector_destroy(*vout);
			*vout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		(*vout)->value = value;
		(*vout)->value_lo = value + np+pad + (np+pad)%2;
		(*vout)->work = (LIS_SCALAR *)lis_malloc( 32*sizeof(LIS_SCALAR),"lis_vector_duplicateex::vout->work" );
		if( NULL==(*vout)->work )
		{
			LIS_SETERR_MEM(32*sizeof(LIS_SCALAR));
			lis_vector_destroy(*vout);
			*vout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		#ifdef USE_VEC_COMP
		#pragma cdir nodep
		#endif
		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0;i<np+pad;i++)
		{
			(*vout)->value[i]    = 0.0;
			(*vout)->value_lo[i] = 0.0;
		}
	}

	#ifdef USE_MPI
		ranges = (int *)lis_malloc( (nprocs+1)*sizeof(int),"lis_vector_duplicateex::ranges" );
		if( ranges==NULL )
		{
			LIS_SETERR_MEM((nprocs+1)*sizeof(int));
			lis_vector_destroy(*vout);
			*vout = NULL;
			return LIS_OUT_OF_MEMORY;
		}
		for(i=0;i<nprocs+1;i++) ranges[i] = ((LIS_VECTOR)A)->ranges[i];
		(*vout)->ranges      = ranges;
	#else
		(*vout)->ranges      = NULL;
	#endif


	(*vout)->is_copy     = LIS_TRUE;
	(*vout)->status      = LIS_VECTOR_ASSEMBLED;
	(*vout)->precision   = precision;
	(*vout)->n           = ((LIS_VECTOR)A)->n;
	(*vout)->gn          = ((LIS_VECTOR)A)->gn;
	(*vout)->np          = ((LIS_VECTOR)A)->np;
	(*vout)->pad         = ((LIS_VECTOR)A)->pad;
	(*vout)->comm        = ((LIS_VECTOR)A)->comm;
	(*vout)->my_rank     = ((LIS_VECTOR)A)->my_rank;
	(*vout)->nprocs      = ((LIS_VECTOR)A)->nprocs;
	(*vout)->is          = ((LIS_VECTOR)A)->is;
	(*vout)->ie          = ((LIS_VECTOR)A)->ie;
	(*vout)->origin      = ((LIS_VECTOR)A)->origin;
	(*vout)->is_destroy  = ((LIS_VECTOR)A)->is_destroy;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_value0"
int lis_vector_set_value0(int flag, int i, LIS_SCALAR value, LIS_VECTOR v)
{
	int n,np,gn,is,ie;

	LIS_DEBUG_FUNC_IN;

	np  = v->np;
	n   = v->n;
	gn  = v->gn;
	is  = v->is;
	ie  = v->ie;
	if( v->origin ) i--;
	if( i<is || i>=ie )
	{
		if( v->origin )
		{
			is++;
			ie++;
			i++;
		}
		LIS_SETERR3(LIS_ERR_ILL_ARG, "i(=%d) is less than %d or larger than %d\n",i,is,ie);
		return LIS_ERR_ILL_ARG;
	}

	if(v->status==LIS_VECTOR_NULL)
	{
		v->value = (LIS_SCALAR *)lis_malloc( np*sizeof(LIS_SCALAR),"lis_vector_set_value::v->value" );
		if( NULL==v->value )
		{
			LIS_SETERR_MEM(np*sizeof(LIS_SCALAR));
			return LIS_OUT_OF_MEMORY;
		}
		v->is_copy = LIS_TRUE;
		v->status  = LIS_VECTOR_ASSEMBLING;
	}
	if(flag==LIS_INS_VALUE)
	{
		v->value[i-is] = value;
	}
	else
	{
		v->value[i-is] += value;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_value"
int lis_vector_set_value(int flag, int i, LIS_SCALAR value, LIS_VECTOR v)
{
	int n,np,gn,is,ie;

	LIS_DEBUG_FUNC_IN;

	np  = v->np;
	n   = v->n;
	gn  = v->gn;
	is  = v->is;
	ie  = v->ie;
	if( v->origin ) i--;
	if( i<is || i>=ie )
	{
		if( v->origin )
		{
			is++;
			ie++;
			i++;
		}
		LIS_SETERR3(LIS_ERR_ILL_ARG, "i(=%d) is less than %d or larger than %d\n",i,is,ie);
		return LIS_ERR_ILL_ARG;
	}

	if(v->status==LIS_VECTOR_NULL)
	{
		v->value = (LIS_SCALAR *)lis_malloc( np*sizeof(LIS_SCALAR),"lis_vector_set_value::v->value" );
		if( NULL==v->value )
		{

			LIS_SETERR_MEM(np*sizeof(LIS_SCALAR));
			return LIS_OUT_OF_MEMORY;
		}
		v->is_copy = LIS_TRUE;
		v->status  = LIS_VECTOR_ASSEMBLING;
	}
	if(flag==LIS_INS_VALUE)
	{
		v->value[i-is] = value;
	}
	else
	{
		v->value[i-is] += value;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_values"
int lis_vector_set_values(int flag, int count, int index[], LIS_SCALAR value[], LIS_VECTOR v)
{
	int n,np,gn,i,ii,is,ie;

	LIS_DEBUG_FUNC_IN;

	np  = v->np;
	n   = v->n;
	gn  = v->gn;
	is  = v->is;
	ie  = v->ie;
	if(v->status==LIS_VECTOR_NULL)
	{
		v->value = (LIS_SCALAR *)lis_malloc( np*sizeof(LIS_SCALAR),"lis_vector_set_values::v->value" );
		if( NULL==v->value )
		{
			LIS_SETERR_MEM(np*sizeof(LIS_SCALAR));
			return LIS_OUT_OF_MEMORY;
		}
		v->is_copy = LIS_TRUE;
		v->status  = LIS_VECTOR_ASSEMBLING;
	}
	if(flag==LIS_INS_VALUE)
	{
		for(i=0;i<count;i++)
		{
			ii = index[i];
			if( v->origin ) ii--;
			if( ii<is || ii>=ie )
			{
				if( v->origin )
				{
					is++;
					ie++;
					ii++;
					i++;
				}
				LIS_SETERR4(LIS_ERR_ILL_ARG, "index[%d](=%d) is less than %d or larger than %d\n",i,ii,is,ie);
				return LIS_ERR_ILL_ARG;
			}
			v->value[ii-is] = value[i];
		}
	}
	else
	{
		for(i=0;i<count;i++)
		{
			ii = index[i];
			if( v->origin ) ii++;
			if( ii<is || ii>=ie )
			{
				if( v->origin )
				{
					is++;
					ie++;
					ii++;
					i++;
				}
				LIS_SETERR4(LIS_ERR_ILL_ARG, "index[%d](=%d) is less than %d or larger than %d\n",i,ii,is,ie);
				return LIS_ERR_ILL_ARG;
			}
			v->value[ii-is] += value[i];
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_set_values"
int lis_vector_set_values2(int flag, int start, int count, LIS_SCALAR value[], LIS_VECTOR v)
{
	int n,np,gn,i,is,ie;

	LIS_DEBUG_FUNC_IN;

	np  = v->np;
	n   = v->n;
	gn  = v->gn;
	is  = v->is;
	ie  = v->ie;
	if(v->status==LIS_VECTOR_NULL)
	{
		v->value = (LIS_SCALAR *)lis_malloc( np*sizeof(LIS_SCALAR),"lis_vector_set_values::v->value" );
		if( NULL==v->value )
		{
			LIS_SETERR_MEM(np*sizeof(LIS_SCALAR));
			return LIS_OUT_OF_MEMORY;
		}
		v->is_copy = LIS_TRUE;
		v->status  = LIS_VECTOR_ASSEMBLING;
	}
	if(flag==LIS_INS_VALUE)
	{
		for(i=0;i<count;i++)
		{
			start = i;
			if( v->origin ) start--;
			if( start<is || start>=ie )
			{
				if( v->origin )
				{
					is++;
					ie++;
					start++;
					i++;
				}
				LIS_SETERR3(LIS_ERR_ILL_ARG, "%d is less than %d or larger than %d\n",start,is,ie);
				return LIS_ERR_ILL_ARG;
			}
			v->value[start-is] = value[i];
		}
	}
	else
	{
		for(i=0;i<count;i++)
		{
			start = i;
			if( v->origin ) start++;
			if( start<is || start>=ie )
			{
				if( v->origin )
				{
					is++;
					ie++;
					start++;
					i++;
				}
				LIS_SETERR3(LIS_ERR_ILL_ARG, "%d is less than %d or larger than %d\n",start,is,ie);
				return LIS_ERR_ILL_ARG;
			}
			v->value[start-is] += value[i];
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_get_size"
int lis_vector_get_size(LIS_VECTOR v, int *local_n, int *global_n)
{
	int		err;

	LIS_DEBUG_FUNC_IN;

	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	*local_n  = v->n;
	*global_n = v->gn;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_get_range"
int lis_vector_get_range(LIS_VECTOR v, int *is, int *ie)
{
	int		err;

	LIS_DEBUG_FUNC_IN;

	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	if( v->origin )
	{
		*is = v->is + 1;
		*ie = v->ie + 1;
	}
	else
	{
		*is = v->is;
		*ie = v->ie;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_get_value"
int lis_vector_get_value(LIS_VECTOR v, int i, LIS_SCALAR *value)
{
	int err,n,gn,is,ie;

	LIS_DEBUG_FUNC_IN;

	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	n   = v->n;
	gn  = v->gn;
	is  = v->is;
	ie  = v->ie;
	if( v->origin ) i--;
	if( i<is || i>=ie )
	{
		if( v->origin )
		{
			i++;
			is++;
			ie++;
		}
		LIS_SETERR3(LIS_ERR_ILL_ARG, "i(=%d) is less than %d or larger than %d\n",i,is,ie);
		return LIS_ERR_ILL_ARG;
	}

	*value = v->value[i-is];
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_get_values"
int lis_vector_get_values(LIS_VECTOR v, int start, int count, LIS_SCALAR value[])
{
	int err,n,gn,i,is,ie;

	LIS_DEBUG_FUNC_IN;

	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	n   = v->n;
	gn  = v->gn;
	is  = v->is;
	ie  = v->ie;
	if( v->origin ) start--;
	if( start<is || start>=ie )
	{
		if( v->origin )
		{
			start++;
			is++;
			ie++;
		}
		LIS_SETERR3(LIS_ERR_ILL_ARG, "start(=%d) is less than %d or larger than %d\n",start,is,ie);
		return LIS_ERR_ILL_ARG;
	}
	if( (start-is+count)>n )
	{
		LIS_SETERR3(LIS_ERR_ILL_ARG, "start(=%d) + count(=%d) exceeds the range of vector v(=%d).\n",start,count,ie);
		return LIS_ERR_ILL_ARG;
	}
	for(i=0;i<count;i++)
	{
		value[i] = v->value[start-is + i];
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

/*
#undef __FUNC__
#define __FUNC__ "lis_vector_set_destroyflag"
int lis_vector_set_destroyflag(LIS_VECTOR v, int flag)
{
	int			err;

	LIS_DEBUG_FUNC_IN;

	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	if( flag )
	{
		v->is_destroy = LIS_TRUE;
	}
	else
	{
		v->is_destroy = LIS_FALSE;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_vector_get_destroyflag"
int lis_vector_get_destroyflag(LIS_VECTOR v, int *flag)
{
	int			err;

	LIS_DEBUG_FUNC_IN;

	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	*flag = v->is_destroy;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
*/

int lis_vector_print(LIS_VECTOR x)
{
#ifdef USE_MPI
	int err,i,ii,is,n,k,nprocs,my_rank;

	err = lis_vector_check(x,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	nprocs  = x->nprocs;
	my_rank = x->my_rank;
	n       = x->n;
	is      = x->is;

	for(k=0;k<nprocs;k++)
	{
		if( k==my_rank )
		{
			for(i=0;i<n;i++)
			{
				ii = i+is;
				if( x->origin ) ii++;
				printf("%6d  %e\n",ii,x->value[i]);
			}
		}
		MPI_Barrier(x->comm);
	}
	return LIS_SUCCESS;
#else
	int err,i,ii,n;

	err = lis_vector_check(x,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	n = x->n;

	for(i=0; i<n; i++)
	{
		ii = i;
		if( x->origin ) ii++;
		if( x->precision==LIS_PRECISION_DEFAULT )
		{
			printf("%6d  %e\n",ii,x->value[i]);
		}
		else
		{
			printf("%6d  %e,%e\n",ii,x->value[i],x->value_lo[i]);
		}
	}

	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_vector_scatter"
int lis_vector_scatter(LIS_SCALAR value[], LIS_VECTOR v)
{
#ifdef USE_MPI
	int err,i,is,n,my_rank,nprocs,*sendcounts;

	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	my_rank = v->my_rank;
	nprocs  = v->nprocs;
	n       = v->n;
	is      = v->is;

	sendcounts = (int *)lis_malloc( (nprocs+1)*sizeof(int),"lis_vector_scatter::sendcounts" );
	for(i=0; i<nprocs; i++)
	{
	  sendcounts[i] = v->ranges[i+1] - v->ranges[i];
	}
	MPI_Scatterv(&value[0],sendcounts,v->ranges,MPI_DOUBLE,&value[is],n,MPI_DOUBLE,0,v->comm);
	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#endif
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0; i<n; i++)
	{
	  v->value[i] = value[i+is];
	}

	return LIS_SUCCESS;
#else
	int err,i,n;

	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	n = v->n;

	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#endif
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0; i<n; i++)
	{
	  v->value[i] = value[i];
	}

	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_vector_gather"
int lis_vector_gather(LIS_VECTOR v, LIS_SCALAR value[])
{
#ifdef USE_MPI
	int err,i,j,is,n,my_rank,nprocs,*recvcounts;

	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	my_rank = v->my_rank;
	nprocs  = v->nprocs;
	n       = v->n;
	is      = v->is;

	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#endif
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0; i<n; i++)
	{
	  value[i+is] = v->value[i];
	}
	recvcounts = (int *)lis_malloc( (nprocs+1)*sizeof(int),"lis_vector_gather::recvcounts" );
	for(i=0; i<nprocs; i++)
	{
	  recvcounts[i] = v->ranges[i+1] - v->ranges[i];
	}
	MPI_Allgatherv(&value[is],n,MPI_DOUBLE,&value[0],recvcounts,v->ranges,MPI_DOUBLE,v->comm);

	return LIS_SUCCESS;
#else
	int err,i,n;

	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	n = v->n;

	#ifdef USE_VEC_COMP
	#pragma cdir nodep
	#endif
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0; i<n; i++)
	{
	  value[i] = v->value[i];
	}

	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_vector_is_null"
int lis_vector_is_null(LIS_VECTOR v)
{
	LIS_DEBUG_FUNC_IN;

	if( v->status==LIS_VECTOR_NULL ) return LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_FALSE;
}
