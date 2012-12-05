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
 * function                    | SOM |
 *-----------------------------+-----+
 * lis_matrix_set              | o   |
 * lis_matrix_setDLU           | o   |
 * lis_matrix_malloc           | o   |
 * lis_matrix_elements_copy    | o   |
 * lis_matrix_transpose        | o   |
 * lis_matrix_split            | o   |
 * lis_matrix_merge            | o   |
 *-----------------------------+-----+-----+
 * function                    |merge|split|
 *-----------------------------+-----+-----|
 * lis_matrix_convert          | o   |     |
 * lis_matrix_copy             | o   | o   |
 * lis_matrix_get_diagonal     | o   | o   |
 * lis_matrix_scaling          | o   | o   |
 * lis_matrix_scaling_symm     | o   | o   |
 * lis_matrix_normf            | o   | o   |
 * lis_matrix_sort             | o   | o   |
 * lis_matrix_solve            | xxx | o   |
 * lis_matrix_solvet           | xxx | o   |
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_dns"
int lis_matrix_set_dns(LIS_SCALAR *value, LIS_MATRIX A)
{
	int err;

	LIS_DEBUG_FUNC_IN;

#if 0
	err = lis_matrix_check(A,LIS_MATRIX_CHECK_SET);
	if( err ) return err;
#else
	if(lis_matrix_is_assembled(A))  return LIS_SUCCESS;
	else {
	  err = lis_matrix_check(A,LIS_MATRIX_CHECK_SET);
	  if( err ) return err;
	}
#endif

	A->value       = value;
	A->is_copy     = LIS_FALSE;
	A->status      = -LIS_MATRIX_DNS;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_dns"
int lis_matrix_malloc_dns(int n, int np, LIS_SCALAR **value)
{
	LIS_DEBUG_FUNC_IN;

	*value   = NULL;

	*value = (LIS_SCALAR *)lis_malloc( n*np*sizeof(LIS_SCALAR),"lis_matrix_malloc_dns::value" );
	if( *value==NULL )
	{
		LIS_SETERR_MEM(n*np*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_elements_copy_dns"
int lis_matrix_elements_copy_dns(int n, int np, LIS_SCALAR *value, LIS_SCALAR *o_value)
{
	int			i,j,is,ie;
	int			nprocs,my_rank;

	LIS_DEBUG_FUNC_IN;

	#ifdef _OPENMP
		nprocs = omp_get_max_threads();
	#else
		nprocs = 1;
	#endif

	#ifdef _OPENMP
	#pragma omp parallel private(i,j,is,ie,my_rank)
	#endif
	{
		#ifdef _OPENMP
			my_rank = omp_get_thread_num();
		#else
			my_rank = 0;
		#endif
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

		for(j=0;j<np;j++)
		{
			for(i=is;i<ie;i++)
			{
				o_value[j*n + i]   = value[j*n + i];
			}
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_copy_dns"
int lis_matrix_copy_dns(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	int			err;
	int			i,n,np;
	LIS_SCALAR	*value;
	LIS_MATRIX_DIAG	D;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;
	np      = Ain->np;
	value   = NULL;
	
	err = lis_matrix_malloc_dns(n,np,&value);
	if( err )
	{
		return err;
	}

	lis_matrix_elements_copy_dns(n,np,Ain->value,value);

	if( Ain->is_splited )
	{
		err = lis_matrix_diag_duplicateM(Ain,&D);
		if( err )
		{
			lis_free(value);
			return err;
		}

		#ifdef _OPENMP
		#pragma omp parallel for private(i)
		#endif
		for(i=0;i<n;i++)
		{
			D->value[i] = Ain->value[i*n + i];
		}
		Aout->D = D;
	}

	err = lis_matrix_set_dns(value,Aout);
	if( err )
	{
		lis_free(value);
		return err;
	}

	err = lis_matrix_assemble(Aout);
	if( err )
	{
		lis_matrix_storage_destroy(Aout);
		return err;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_diagonal_dns"
int lis_matrix_get_diagonal_dns(LIS_MATRIX A, LIS_SCALAR d[])
{
	int i;
	int n;

	LIS_DEBUG_FUNC_IN;

	n    = A->n;
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0; i<n; i++)
	{
		d[i] = A->value[i*n + i];
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scaling_dns"
int lis_matrix_scaling_dns(LIS_MATRIX A, LIS_SCALAR d[])
{
	int i,j;
	int n,np;

	LIS_DEBUG_FUNC_IN;

	n  = A->n;
	np = A->np;
	for(j=0;j<np;j++)
	{
		for(i=0;i<n;i++)
		{
			A->value[j*n + i] *= d[i];
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scaling_symm_dns"
int lis_matrix_scaling_symm_dns(LIS_MATRIX A, LIS_SCALAR d[])
{
	int i,j;
	int n,np;

	LIS_DEBUG_FUNC_IN;

	n  = A->n;
	np = A->np;
	for(j=0;j<np;j++)
	{
		for(i=0;i<n;i++)
		{
			A->value[j*n + i] *= d[i]*d[j];
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_normf_dns"
int lis_matrix_normf_dns(LIS_MATRIX A, LIS_SCALAR *nrm)
{
	int i,j;
	int n;
	LIS_SCALAR sum;

	LIS_DEBUG_FUNC_IN;

	n    = A->n;
	sum  = (LIS_SCALAR)0;
	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for reduction(+:sum) private(i,j)
		#endif
		for(i=0; i<n; i++)
		{
			sum += A->D->value[i]*A->D->value[i];
			for(j=A->L->index[i];j<A->L->index[i+1];j++)
			{
				sum += A->L->value[j]*A->L->value[j];
			}
			for(j=A->U->index[i];j<A->U->index[i+1];j++)
			{
				sum += A->U->value[j]*A->U->value[j];
			}
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for reduction(+:sum) private(i,j)
		#endif
		for(i=0; i<n; i++)
		{
			sum += A->value[i]*A->value[i];
			for(j=A->index[i];j<A->index[i+1];j++)
			{
				sum += A->value[j]*A->value[j];
			}
		}
	}
	*nrm = sqrt(sum);
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_transpose_dns"
int lis_matrix_transpose_dns(LIS_MATRIX Ain, LIS_MATRIX *Aout)
{

	LIS_DEBUG_FUNC_IN;

/*	err = lis_matrix_convert_dns2ccs(Ain,Aout);*/
	(*Aout)->matrix_type = LIS_MATRIX_DNS;
	(*Aout)->status      = LIS_MATRIX_DNS;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_split_dns"
int lis_matrix_split_dns(LIS_MATRIX A)
{
	int i,n;
	int				err;
	LIS_MATRIX_DIAG	D;

	LIS_DEBUG_FUNC_IN;

	n = A->n;

	err = lis_matrix_diag_duplicateM(A,&D);
	if( err )
	{
		return err;
	}

	for(i=0;i<n;i++)
	{
		D->value[i] = A->value[i*n + i];
	}
	A->D          = D;
	A->is_splited = LIS_TRUE;
	A->is_save    = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_matrix_merge_dns"
int lis_matrix_merge_dns(LIS_MATRIX A)
{
	LIS_DEBUG_FUNC_IN;

	A->is_splited = LIS_FALSE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_sort_dns"
int lis_matrix_sort_dns(LIS_MATRIX A)
{
	LIS_DEBUG_FUNC_IN;

	A->is_sorted = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_solve_dns"
int lis_matrix_solve_dns(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag)
{
	int i,j,n,np;
	LIS_SCALAR	t;
	LIS_SCALAR	*b,*x;

	LIS_DEBUG_FUNC_IN;

	n  = A->n;
	np = A->np;
	b  = B->value;
	x  = X->value;

	switch(flag)
	{
	case LIS_MATRIX_LOWER:
		for(i=0;i<n;i++)
		{
			t = b[i];
			for(j=0;j<i;j++)
			{
				t -= A->value[j*n + i] * x[j];
			}
			x[i]   = t * A->WD->value[i];
		}
		break;
	case LIS_MATRIX_UPPER:
		for(i=n-1;i>=0;i--)
		{
			t = b[i];
			for(j=i+1;j<np;j++)
			{
				t -= A->value[j*n + i] * x[j];
			}
			x[i]   = t * A->WD->value[i];
		}
		break;
	case LIS_MATRIX_SSOR:
		for(i=0;i<n;i++)
		{
			t = b[i];
			for(j=0;j<i;j++)
			{

				t -= A->value[j*n + i] * x[j];
			}
			x[i]   = t * A->WD->value[i];
		}
		for(i=n-1;i>=0;i--)
		{
			t = 0.0;
			for(j=i+1;j<n;j++)
			{
				t += A->value[j*n + i] * x[j];
			}
			x[i]  -= t * A->WD->value[i];
		}
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_solvet_dns"
int lis_matrix_solvet_dns(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag)
{
	int i,j,n,np;
	LIS_SCALAR	t;
	LIS_SCALAR	*b,*x;

	LIS_DEBUG_FUNC_IN;

	n  = A->n;
	np = A->np;
	b  = B->value;
	x  = X->value;

	lis_vector_copy(B,X);
	switch(flag)
	{
	case LIS_MATRIX_LOWER:
		for(i=0;i<n;i++)
		{
			x[i]   = x[i] * A->WD->value[i];
			for(j=i+1;j<np;j++)
			{
				x[j] -= A->value[j*n + i] * x[i];
			}
		}
		break;
	case LIS_MATRIX_UPPER:
		for(i=n-1;i>=0;i--)
		{
			x[i]   = x[i] * A->WD->value[i];
			for(j=0;j<i;j++)
			{
				x[j] -= A->value[j*n + i] * x[i];
			}
		}
		break;
	case LIS_MATRIX_SSOR:
		for(i=0;i<n;i++)
		{
			t   = x[i] * A->WD->value[i];
			for(j=i+1;j<np;j++)
			{
				x[j] -= A->value[j*n + i] * t;
			}
		}
		for(i=n-1;i>=0;i--)
		{
			t    = x[i] * A->WD->value[i];
			x[i] = t;
			for(j=0;j<i;j++)
			{
				x[j] -= A->value[j*n + i] * t;
			}
		}
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_crs2dns"
int lis_matrix_convert_crs2dns(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	int			i,j;
	int			err;
	int			n,np,nprocs,my_rank;
	int			is,ie;
	LIS_SCALAR	*value;

	LIS_DEBUG_FUNC_IN;


	n       = Ain->n;
	np      = Ain->np;
	#ifdef _OPENMP
		nprocs  = omp_get_max_threads();
	#else
		nprocs  = 1;
	#endif

	value    = NULL;

	err = lis_matrix_malloc_dns(n,np,&value);
	if( err )
	{
		return err;
	}

	/* convert dns */
	#ifdef _OPENMP
	#pragma omp parallel private(i,j,is,ie,my_rank)
	#endif
	{
		#ifdef _OPENMP
			my_rank = omp_get_thread_num();
		#else
			my_rank = 0;
		#endif
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

		for(j=0;j<np;j++)
		{
			for(i=is;i<ie;i++)
			{
				value[j*n+i] = (LIS_SCALAR)0.0;
			}
		}
		for(i=is;i<ie;i++)
		{
			for(j=Ain->ptr[i];j<Ain->ptr[i+1];j++)
			{
				value[i + n*Ain->index[j]] = Ain->value[j];
			}
		}
	}

	err = lis_matrix_set_dns(value,Aout);
	if( err )
	{
		lis_free(value);
		return err;
	}
	err = lis_matrix_assemble(Aout);
	if( err )
	{
		lis_matrix_storage_destroy(Aout);
		return err;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_dns2crs"
int lis_matrix_convert_dns2crs(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	int			i,j,k;
	int			err;
	int			gn,n,np,nnz,is,ie;
	int			*ptr,*index;
	LIS_SCALAR	*value;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;
	np      = Ain->np;
	gn      = Ain->gn;
	is      = Ain->is;
	ie      = Ain->ie;

	ptr     = NULL;
	index   = NULL;
	value   = NULL;

	ptr = (int *)lis_malloc( (n+1)*sizeof(int),"lis_matrix_convert_dns2crs::ptr" );
	if( ptr==NULL )
	{
		LIS_SETERR_MEM((n+1)*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}

	#ifdef _OPENMP
	#pragma omp parallel for private(i,j)
	#endif
	for(i=0;i<n;i++)
	{
		ptr[i+1] = 0;
		for(j=0;j<np;j++)
		{
			if( Ain->value[j*n+i]!=(LIS_SCALAR)0.0 )
			{
				ptr[i+1]++;
			}
		}
	}
	ptr[0] = 0;
	for(i=0;i<n;i++)
	{
		ptr[i+1] += ptr[i];
	}
	nnz = ptr[n];

	index = (int *)lis_malloc( nnz*sizeof(int),"lis_matrix_convert_dns2crs::index" );
	if( index==NULL )
	{
		lis_free2(3,ptr,index,value);
		LIS_SETERR_MEM(nnz*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	value = (LIS_SCALAR *)lis_malloc( nnz*sizeof(LIS_SCALAR),"lis_matrix_convert_dns2crs::value" );
	if( value==NULL )
	{
		lis_free2(3,ptr,index,value);
		LIS_SETERR_MEM(nnz*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}

	/* convert crs */
	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,k)
	#endif
	for(i=0;i<n;i++)
	{
		k = ptr[i];
		for(j=0;j<np;j++)
		{
			if( Ain->value[j*n + i]!=(LIS_SCALAR)0.0 )
			{
				value[k] = Ain->value[j*n + i];
				index[k] = j;
				k++;
			}
		}
	}

	err = lis_matrix_set_crs(nnz,ptr,index,value,Aout);
	if( err )
	{
		lis_free2(3,ptr,index,value);
		return err;
	}
	err = lis_matrix_assemble(Aout);
	if( err )
	{
		lis_matrix_storage_destroy(Aout);
		return err;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

