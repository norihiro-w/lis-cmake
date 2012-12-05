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
#include <ctype.h>
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_input
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_input_mm"
int lis_input_mm(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, FILE *file)
{
	int			err;
	int			matrix_type;
	LIS_MATRIX	B;

	LIS_DEBUG_FUNC_IN;

	matrix_type = A->matrix_type;

	err = lis_input_mm_crs(A,b,x,file);
	if( err ) return err;

	if( matrix_type!=LIS_MATRIX_CRS )
	{
		err = lis_matrix_duplicate(A,&B);
		if( err ) return err;
		lis_matrix_set_type(B,matrix_type);
		err = lis_matrix_convert(A,B);
		if( err ) return err;
		lis_matrix_storage_destroy(A);
		lis_matrix_DLU_destroy(A);
		lis_matrix_diag_destroy(A->WD);
		if( A->l2g_map ) lis_free( A->l2g_map );
		if( A->commtable ) lis_commtable_destroy( A->commtable );
		if( A->ranges ) lis_free( A->ranges );
		err = lis_matrix_copy_struct(B,A);
		if( err ) return err;
		lis_free(B);
		if( A->matrix_type==LIS_MATRIX_JDS )
		{
			A->work = (LIS_SCALAR *)lis_malloc(A->n*sizeof(LIS_SCALAR),"lis_input_mm::A->work");
			if( A->work==NULL )
			{
				LIS_SETERR_MEM(A->n*sizeof(LIS_SCALAR));
				return LIS_OUT_OF_MEMORY;
			}
		}
	}


	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_input_mm_vec"
int lis_input_mm_vec(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, FILE *file, int isb, int isx, int isbin)
{
	char			buf[BUFSIZE];
	int				i;
	int				mode;
	int				gn,n,is,ie;
	int				idx;
	LIS_SCALAR		val;
	LIS_MM_VECFMT	vecfmt;

	LIS_DEBUG_FUNC_IN;

	if( isb==0 && isx==0 ) return LIS_SUCCESS;

	gn  = A->gn;
	n   = A->n;
	is  = A->is;
	ie  = A->ie;
	mode = 1;
	mode = *(char *)&mode;
	if( mode!=(isbin-1) )
	{
		mode = LIS_TRUE;			
	}
	else
	{
		mode = LIS_FALSE;
	}
	if( isb )
	{
		lis_vector_set_size(b,n,0);
		for(i=0;i<gn;i++)
		{
			if( isbin )
			{
				if( fread(&vecfmt, sizeof(vecfmt), 1, file)!=1 )
				{
					LIS_SETERR_FIO;
					return LIS_ERR_FILE_IO;
				}
				idx = vecfmt.i;
				val = vecfmt.value;
				if( mode )
				{
					lis_bswap_int(1,&idx);
					lis_bswap_double(1,&val);
				}
			}
			else
			{
				if( fgets(buf, BUFSIZE, file) == NULL )
				{
					LIS_SETERR_FIO;
					return LIS_ERR_FILE_IO;
				}
				if( sscanf(buf, "%d %lg", &idx, &val) != 2 )
				{
					LIS_SETERR_FIO;
					return LIS_ERR_FILE_IO;
				}
			}
			idx--;
			if( idx>=is && idx<ie )
			{
				b->value[idx-is] = val;
			}
		}
	}
	if( isx )
	{
		lis_vector_set_size(x,n,0);
		for(i=0;i<gn;i++)
		{
			if( isbin )
			{
				if( fread(&vecfmt, sizeof(vecfmt), 1, file)!=1 )
				{
					LIS_SETERR_FIO;
					return LIS_ERR_FILE_IO;
				}
				idx = vecfmt.i;
				val = vecfmt.value;
				if( mode )
				{
					lis_bswap_int(1,&idx);
					lis_bswap_double(1,&val);
				}
			}
			else
			{
				if( fgets(buf, BUFSIZE, file) == NULL )
				{
					LIS_SETERR_FIO;
					return LIS_ERR_FILE_IO;
				}
				if( sscanf(buf, "%d %lg", &idx, &val) != 2 )
				{
					LIS_SETERR_FIO;
					return LIS_ERR_FILE_IO;
				}
			}
			idx--;
			if( idx>=is && idx<ie )
			{
				x->value[idx-is] = val;
			}
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_input_mm_banner"
int lis_input_mm_banner(FILE *file, int *mmtype)
{
	char			buf[BUFSIZE];
	char			banner[64], mtx[64], fmt[64], dtype[64], dstruct[64];
	char			*p;

	LIS_DEBUG_FUNC_IN;

	/* check banner */
	if( fgets(buf, BUFSIZE, file) == NULL )
	{
		LIS_SETERR_FIO;
		return LIS_ERR_FILE_IO;
	}
	sscanf(buf, "%s %s %s %s %s", banner, mtx, fmt, dtype, dstruct);

	for(p=mtx;*p!='\0';p++)     *p = (char)tolower(*p);
	for(p=fmt;*p!='\0';p++)     *p = (char)tolower(*p);
	for(p=dtype;*p!='\0';p++)   *p = (char)tolower(*p);
	for(p=dstruct;*p!='\0';p++) *p = (char)tolower(*p);

	if( strncmp(banner, MM_BANNER, strlen(MM_BANNER))!=0 || strncmp(mtx, MM_MTX, strlen(MM_MTX))!=0 )
	{
		LIS_SETERR(LIS_ERR_FILE_IO,"Not Matrix Market banner\n");
		return LIS_ERR_FILE_IO;
	}
	if( strncmp(fmt, MM_FMT, strlen(MM_FMT))!=0 )
	{
		LIS_SETERR(LIS_ERR_FILE_IO,"Not Coodinate format\n");
		return LIS_ERR_FILE_IO;
	}
	if( strncmp(dtype, MM_TYPE_REAL, strlen(MM_TYPE_REAL))!=0 )
	{
		LIS_SETERR(LIS_ERR_FILE_IO,"Not real\n");
		return LIS_ERR_FILE_IO;
	}
	if( strncmp(dstruct, MM_TYPE_GENERAL, strlen(MM_TYPE_GENERAL))==0 )
	{
		*mmtype = MM_GENERAL;
	}
	else if( strncmp(dstruct, MM_TYPE_SYMM, strlen(MM_TYPE_SYMM))==0)
	{
		*mmtype = MM_SYMM;
	}
	else
	{
		LIS_SETERR(LIS_ERR_FILE_IO,"Not general or symmetric\n");
		return LIS_ERR_FILE_IO;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_input_mm_size"
int lis_input_mm_size(FILE *file, int *nr, int *nc, int *nnz, int *isb, int *isx, int *isbin)
{
	char			buf[BUFSIZE];
	int				err;

	LIS_DEBUG_FUNC_IN;
	/* check size */		
	do
	{
		if( fgets(buf, BUFSIZE, file) == NULL )
		{
			LIS_SETERR_FIO;
			return LIS_ERR_FILE_IO;
		}
	}while( buf[0]=='%' );
	err = sscanf(buf, "%d %d %d %d %d %d", nr, nc, nnz, isb, isx, isbin);
	if( err==3 )
	{
		*isb   = 0;
		*isx   = 0;
		*isbin = 0;
	}
	else if( err==5 )
	{
		*isbin = 0;
	}

	if( *nr!=*nc )
	{
		LIS_SETERR(LIS_ERR_FILE_IO,"matrix is not square\n");
		return LIS_ERR_FILE_IO;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

int lis_input_mm_crs(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, FILE *file)
{
	char			buf[BUFSIZE];
	int				nr,nc,nnz;
	int				i,j,my_rank;
	int				err;
	int				mmtype,mode;
	int				n,is,ie;
	int				ridx,cidx;
	int				*ptr, *index;
	int				*work;
	int				isb,isx,isbin;
	LIS_SCALAR		val;
	LIS_SCALAR		*value;
	LIS_MM_MATFMT	matfmt;

	LIS_DEBUG_FUNC_IN;

	#ifdef USE_MPI
		my_rank = A->my_rank;
	#else
		my_rank = 0;
	#endif
	
	/* check banner */
	err = lis_input_mm_banner(file,&mmtype);
	if( err ) return err;

	/* check size */		
	err = lis_input_mm_size(file,&nr,&nc,&nnz,&isb,&isx,&isbin);
	if( err ) return err;

	err = lis_matrix_set_size(A,0,nr);
	if( err ) return err;

	if( my_rank==0 ) printf("matrix size = %d x %d (%d nonzero entries)\n",nr,nc,nnz);

	n      = A->n;
	ptr    = NULL;
	index  = NULL;
	value  = NULL;
	work   = NULL;


	lis_matrix_get_range(A,&is,&ie);

	ptr   = (int *)lis_malloc( (n+1)*sizeof(int),"lis_input_mm_crs::ptr" );
	if( ptr==NULL )
	{
		LIS_SETERR_MEM((n+1)*sizeof(int));
		lis_free2(4,ptr,index,value,work);
		return LIS_OUT_OF_MEMORY;
	}
	work  = (int *)lis_malloc( (n+1)*sizeof(int),"lis_input_mm_crs::work" );
	if( work==NULL )
	{
		LIS_SETERR_MEM((n+1)*sizeof(int));
		lis_free2(4,ptr,index,value,work);
		return LIS_OUT_OF_MEMORY;
	}

	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for(i=0;i<n+1;i++)
	{
		ptr[i]  = 0;
		work[i]  = 0;
	}

	/* read data */
	mode = 1;
	mode = *(char *)&mode;
	if( mode!=(isbin-1) )
	{
		mode = LIS_TRUE;			
	}
	else
	{
		mode = LIS_FALSE;
	}
	for( i=0; i<nnz; i++ )
	{
		if( isbin )
		{
			if( fread(&matfmt, sizeof(matfmt), 1, file)!=1 )
			{
				LIS_SETERR_FIO;
				lis_free2(4,ptr,index,value,work);
				return LIS_ERR_FILE_IO;
			}
			ridx = matfmt.i;
			cidx = matfmt.j;
			if( mode )
			{
				lis_bswap_int(1,&ridx);
				lis_bswap_int(1,&cidx);
			}
		}
		else
		{
			if( fgets(buf, BUFSIZE, file)==NULL )
			{
				LIS_SETERR_FIO;
				lis_free2(4,ptr,index,value,work);
				return LIS_ERR_FILE_IO;
			}
			if( sscanf(buf, "%d %d %lg", &ridx, &cidx, &val) != 3 )
			{
				LIS_SETERR_FIO;
				lis_free2(4,ptr,index,value,work);
				return LIS_ERR_FILE_IO;
			}
		}
/*		if( val!=0.0 )*/
		{
			if( mmtype==MM_SYMM && ridx!=cidx )
			{
				if( cidx>is && cidx<=ie ) work[cidx-is-1]++;
			}
			if( ridx>is && ridx<=ie )
			{
				ptr[ridx-is]++;
			}
		}
	}


	ptr[0] = 0;
	for( i=0; i<n; i++ )
	{
		if( mmtype==MM_SYMM )
		{
			ptr[i+1] += ptr[i] + work[i];
		}
		else
		{
			ptr[i+1] += ptr[i];
		}
		work[i] = 0;
	}

	index   = (int *)lis_malloc( ptr[n]*sizeof(int),"lis_input_mm_crs::index" );
	if( index==NULL )
	{
		LIS_SETERR_MEM(ptr[n]*sizeof(int));
		lis_free2(4,ptr,index,value,work);
		return LIS_OUT_OF_MEMORY;
	}
	value   = (LIS_SCALAR *)lis_malloc( ptr[n]*sizeof(LIS_SCALAR),"lis_input_mm_crs::value" );
	if( value==NULL )
	{
		LIS_SETERR_MEM(ptr[n]*sizeof(LIS_SCALAR));
		lis_free2(4,ptr,index,value,work);
		return LIS_OUT_OF_MEMORY;
	}
	#ifdef _OPENMP
	#pragma omp parallel for private(i,j)
	#endif
	for(i=0;i<n;i++)
	{
		for(j=ptr[i];j<ptr[i+1];j++)
		{
			index[j] = 0;
			value[j] = 0.0;
		}
	}

	rewind(file);
	if( fgets(buf, BUFSIZE, file) == NULL )
	{
		LIS_SETERR_FIO;
		lis_free2(4,ptr,index,value,work);
		return LIS_ERR_FILE_IO;
	}
	do
	{
		if( fgets(buf, BUFSIZE, file) == NULL )
		{
			LIS_SETERR_FIO;
			lis_free2(4,ptr,index,value,work);
			return LIS_ERR_FILE_IO;
		}
	}while( buf[0]=='%' );

	for( i=0; i<nnz; i++ )
	{
		if( isbin )
		{
			if( fread(&matfmt, sizeof(matfmt), 1, file)!=1 )
			{
				LIS_SETERR_FIO;
				lis_free2(4,ptr,index,value,work);
				return LIS_ERR_FILE_IO;
			}
			ridx = matfmt.i;
			cidx = matfmt.j;
			val  = matfmt.value;
			if( mode )
			{
				lis_bswap_int(1,&ridx);
				lis_bswap_int(1,&cidx);
				lis_bswap_double(1,&val);
			}
		}
		else
		{
			if( fgets(buf, BUFSIZE, file) == NULL )
			{
				LIS_SETERR_FIO;
				lis_free2(4,ptr,index,value,work);
				return LIS_ERR_FILE_IO;
			}
			if( sscanf(buf, "%d %d %lg", &ridx, &cidx, &val) != 3 )
			{
				LIS_SETERR_FIO;
				lis_free2(4,ptr,index,value,work);
				return LIS_ERR_FILE_IO;
			}
		}
		ridx--;
		cidx--;
		if( ridx==cidx && val==0.0 )
		{
			printf("diagonal element is zero (i=%d)\n",ridx);
		}
/*		if( val!=0.0 )*/
		{
			if( mmtype==MM_SYMM && ridx!=cidx )
			{
				if( cidx>=is && cidx<ie )
				{
					value[ptr[cidx-is]+work[cidx-is]] = val;
					index[ptr[cidx-is]+work[cidx-is]] = ridx;
					work[cidx-is]++;
				}
			}
			if( ridx>=is && ridx<ie )
			{
				value[ptr[ridx-is]+work[ridx-is]] = val;
				index[ptr[ridx-is]+work[ridx-is]] = cidx;
				work[ridx-is]++;
			}
		}
	}
	#ifdef USE_MPI
		MPI_Barrier(A->comm);
	#endif

	err = lis_matrix_set_crs(ptr[n],ptr,index,value,A);
	if( err )
	{
		lis_free2(4,ptr,index,value,work);
		return err;
	}
	err = lis_matrix_assemble(A);
	if( err )
	{
		lis_matrix_storage_destroy(A);
		lis_free(work);
		return err;
	}

	if( b!=NULL && x!=NULL )
	{
		err = lis_input_mm_vec(A,b,x,file,isb,isx,isbin);
		if( err )
		{
			lis_matrix_storage_destroy(A);
			lis_free(work);
		}
	}
	lis_free(work);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
