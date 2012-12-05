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
 * lis_matrix_set
 * lis_matrix_malloc
 * lis_matrix_copy
 * lis_matrix_convert
 * lis_matrix_get_diagonal
 * lis_matrix_scaling
 * lis_matrix_scaling_symm
 * lis_matrix_normf
 * lis_matrix_transpose
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_vbr"
int lis_matrix_set_vbr(int nnz, int nr, int nc, int bnnz, int *row, int *col, int *ptr, int *bptr, int *bindex, LIS_SCALAR *value, LIS_MATRIX A)
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

	A->row         = row;
	A->col         = col;
	A->ptr         = ptr;
	A->bptr        = bptr;
	A->bindex      = bindex;
	A->value       = value;
	A->is_copy     = LIS_FALSE;
	A->status      = -LIS_MATRIX_VBR;
	A->is_block    = LIS_TRUE;
	A->nnz         = nnz;
	A->bnnz        = bnnz;
	A->nr          = nr;
	A->nc          = nc;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_malloc_vbr"
int lis_matrix_malloc_vbr(int n, int nnz, int nr, int nc, int bnnz, int **row, int **col, int **ptr, int **bptr, int **bindex, LIS_SCALAR **value)
{

	LIS_DEBUG_FUNC_IN;

	*row      = NULL;
	*col      = NULL;
	*ptr      = NULL;
	*bptr     = NULL;
	*bindex   = NULL;
	*value    = NULL;

	*row = (int *)lis_malloc( (nr+1)*sizeof(int),"lis_matrix_malloc_vbr::row" );
	if( *row==NULL )
	{
		LIS_SETERR_MEM((nr+1)*sizeof(int));
		lis_free2(6,*row,*col,*ptr,*bptr,*bindex,*value);
		return LIS_FAILS;
	}
	*col = (int *)lis_malloc( (nc+1)*sizeof(int),"lis_matrix_malloc_vbr::col" );
	if( *col==NULL )
	{
		LIS_SETERR_MEM((nc+1)*sizeof(int));
		lis_free2(6,*row,*col,*ptr,*bptr,*bindex,*value);
		return LIS_FAILS;
	}
	*ptr = (int *)lis_malloc( (bnnz+1)*sizeof(int),"lis_matrix_malloc_vbr::ptr" );
	if( *ptr==NULL )
	{
		LIS_SETERR_MEM((bnnz+1)*sizeof(int));
		lis_free2(6,*row,*col,*ptr,*bptr,*bindex,*value);
		return LIS_FAILS;
	}
	*bptr = (int *)lis_malloc( (nr+1)*sizeof(int),"lis_matrix_malloc_vbr::bptr" );
	if( *bptr==NULL )
	{
		LIS_SETERR_MEM((nr+1)*sizeof(int));
		lis_free2(6,*row,*col,*ptr,*bptr,*bindex,*value);
		return LIS_FAILS;
	}
	*bindex = (int *)lis_malloc( bnnz*sizeof(int),"lis_matrix_malloc_vbr::bindex" );
	if( *bindex==NULL )
	{
		LIS_SETERR_MEM(bnnz*sizeof(int));
		lis_free2(6,*row,*col,*ptr,*bptr,*bindex,*value);
		return LIS_OUT_OF_MEMORY;
	}
	*value = (LIS_SCALAR *)lis_malloc( nnz*sizeof(LIS_SCALAR),"lis_matrix_malloc_vbr::value" );
	if( *value==NULL )
	{
		LIS_SETERR_MEM(nnz*sizeof(LIS_SCALAR));
		lis_free2(6,*row,*col,*ptr,*bptr,*bindex,*value);
		return LIS_OUT_OF_MEMORY;
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_elements_copy_vbr"
int lis_matrix_elements_copy_vbr(int n, int nr, int nc, int bnnz, int *row, int *col, int *ptr, int *bptr, int *bindex, LIS_SCALAR *value,
								 int *o_row, int *o_col, int *o_ptr, int *o_bptr, int *o_bindex, LIS_SCALAR *o_value)
{
	int		bi,bj,i,j,k;

	LIS_DEBUG_FUNC_IN;

	#ifdef _OPENMP
	#pragma omp parallel private(bi,bj,i,j,k)
	#endif
	{
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<nr+1;i++)
		{
			o_row[i]  = row[i];
			o_bptr[i] = bptr[i];
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<nc+1;i++)
		{
			o_col[i] = col[i];
		}

		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(bi=0;bi<nr;bi++)
		{
			for(bj=bptr[bi];bj<bptr[bi+1];bj++)
			{
				k    = ptr[bj];
				for(j=col[bindex[bj]];j<col[bindex[bj]+1];j++)
				{
					for(i=row[bi];i<row[bi+1];i++)
					{
						o_value[k] = value[k];
						k++;
					}
				}
				o_bindex[bj]  = bindex[bj];
				o_ptr[bj+1]   = ptr[bj+1];
			}
		}
		o_ptr[0] = ptr[0];
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_copy_vbr"
int lis_matrix_copy_vbr(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	int			err;
	int			n,nnz,bnnz,nr,nc;
	int			*row,*col,*ptr,*bptr,*bindex;
	LIS_SCALAR	*value;

	LIS_DEBUG_FUNC_IN;

	n       = Ain->n;
	nnz     = Ain->nnz;
	bnnz    = Ain->bnnz;
	nr      = Ain->nr;
	nc      = Ain->nc;
    
	err = lis_matrix_malloc_vbr(n,nnz,nr,nc,bnnz,&row,&col,&ptr,&bptr,&bindex,&value);
	if( err )
	{
		return err;
	}

	lis_matrix_elements_copy_vbr(n,nr,nc,bnnz,Ain->row,Ain->col,Ain->ptr,Ain->bptr,Ain->bindex,Ain->value,row,col,ptr,bptr,bindex,value);

	err = lis_matrix_set_vbr(nnz,nr,nc,bnnz,row,col,ptr,bptr,bindex,value,Aout);
	if( err )
	{
		lis_free2(6,row,col,ptr,bptr,bindex,value);
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
#define __FUNC__ "lis_matrix_get_vbr_rowcol"
int lis_matrix_get_vbr_rowcol(LIS_MATRIX Ain, int *nr, int *nc, int **row, int **col)
{
	int			i,j,k,jj,kk,n;
	int			*iw;

	LIS_DEBUG_FUNC_IN;

	n        = Ain->n;
	iw       = NULL;

	iw = (int *)lis_malloc( (n+1)*sizeof(int),"lis_matrix_get_vbr_rowcol::iw" );
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}

	for(i=0;i<n+1;i++) iw[i] = 0;
	for(i=0;i<n;i++)
	{
		if( Ain->ptr[i]<Ain->ptr[i+1] )
		{
			jj = Ain->index[Ain->ptr[i]];
			iw[jj] = 1;
			for(j=Ain->ptr[i]+1;j<Ain->ptr[i+1];j++)
			{
				jj = Ain->index[j];
				kk = Ain->index[j-1];
				if( kk!=jj-1 )
				{

					iw[jj] = 1;
					iw[kk+1] = 1;
				}
			}
			iw[jj+1] = 1;
		}
	}

	k=0;
	iw[0] = 0;
	for(i=1;i<n+1;i++)
	{
		if( iw[i]!=0 )
		{
			k++;
			iw[k] = i;
		}
	}

	*nr = k;
	*nc = k;
	*row = (int *)lis_malloc((k+1)*sizeof(int),"lis_matrix_get_vbr_rowcol::row");
	if( (*row)==NULL )
	{
		LIS_SETERR_MEM((k+1)*sizeof(int));
		lis_free(iw);
		return LIS_OUT_OF_MEMORY;
	}
	*col = (int *)lis_malloc((k+1)*sizeof(int),"lis_matrix_get_vbr_rowcol::col");
	if( (*col)==NULL )
	{
		LIS_SETERR_MEM((k+1)*sizeof(int));
		lis_free2(2,iw,*row);
		return LIS_OUT_OF_MEMORY;
	}
	memcpy(*row,iw,(k+1)*sizeof(int));
	memcpy(*col,iw,(k+1)*sizeof(int));

	lis_free(iw);
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#if 0
#undef __FUNC__
#define __FUNC__ "lis_matrix_get_vbr_rowcol"
int lis_matrix_get_vbr_rowcol(LIS_MATRIX Ain, int *nr, int *nc, int **row, int **col)
{
	int			i,j,k,l,n;
	int			ii,jj,kk,ret;
	int			bnnz,bj,bnr,bnc,jpos,nnz,ij,kv,bi;
	int			err;
	int			gn,nprocs,my_rank;
	int			is,ie,pe;
	int			*iw;
	int			ac,oc,count;
	int			p[3][5],and[3],or[3];

	LIS_DEBUG_FUNC_IN;

	n        = Ain->n;
	gn       = Ain->gn;
	nprocs   = Ain->nprocs;
	my_rank  = Ain->my_rank;
	is       = Ain->is;
	ie       = Ain->ie;

	iw       = NULL;

	iw = (int *)lis_malloc( n*sizeof(int),"lis_matrix_get_vbr_rowcol::iw" );
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}

	memset(p[0],0,15*sizeof(int));
	count = 0;
	k = 0;
	for(i=0;i<n;i++)
	{
		kk = i - k;
		for(j=Ain->ptr[i];j<Ain->ptr[i+1];j++)
		{
			jj = Ain->index[j] - k;
			if( jj>=0 && jj<5 )
			{
				p[kk][jj] = 1;
			}
		}
		if( kk==1 )
		{
			and[0] = p[kk-1][0] & p[kk][0];
			and[1] = p[kk-1][1] & p[kk][1];
			and[2] = p[kk-1][2] & p[kk][2];
			ac = and[0] + and[1] + and[2];
			if( ac==0 )
			{
				memcpy(p[0],&p[1][1],3*sizeof(int));
				memset(p[1],0,5*sizeof(int));
				iw[count++] = 1;
				k++;

			}
			else if( ac==1 )
			{
				if( and[0]==1 || and[1]==1 )
				{
					memset(p[0],0,10*sizeof(int));
					iw[count++] = 2;
					k += 2;
				}
				else
				{
					memcpy(p[0],&p[1][1],3*sizeof(int));
					memset(p[1],0,5*sizeof(int));
					iw[count++] = 1;
					k++;
				}
			}
			else
			{
				or[0] = p[kk-1][0] | p[kk][0];
				or[1] = p[kk-1][1] | p[kk][1];
				or[2] = p[kk-1][2] | p[kk][2];
				oc = or[0] + or[1] + or[2];
				if( oc==2 )
				{
					memset(p[0],0,10*sizeof(int));
					iw[count++] = 2;
					k += 2;
				}
			}
		}
		else if( kk==2 )
		{
			oc = p[kk][0] + p[kk][1] + p[kk][2];
			if( ac==2 )
			{

				if( oc==3 )
				{
					memset(p[0],0,15*sizeof(int));
					iw[count++] = 3;
					k += 3;
				}
				else
				{
					memcpy(p[0],&p[2][2],3*sizeof(int));
					memset(p[1],0,10*sizeof(int));
					iw[count++] = 2;
					k += 2;
				}
			}
			else
			{
				if( oc==1 )
				{
					memcpy(p[0],&p[2][2],3*sizeof(int));
					memset(p[1],0,10*sizeof(int));
					iw[count++] = 2;
					k += 2;
				}
				else
				{
					memset(p[0],0,15*sizeof(int));
					iw[count++] = 3;
					k += 3;
				}
			}
		}
	}
	if( k<n )
	{
		iw[count++] = 1;
	}

	*nr = count;
	*nc = count;
	*row = (int *)lis_malloc((count+1)*sizeof(int),"lis_matrix_get_vbr_rowcol::row");
	if( (*row)==NULL )
	{
		LIS_SETERR_MEM((count+1)*sizeof(int));
		lis_free(iw);
		return LIS_OUT_OF_MEMORY;
	}
	*col = (int *)lis_malloc((count+1)*sizeof(int),"lis_matrix_get_vbr_rowcol::col");
	if( (*col)==NULL )
	{
		LIS_SETERR_MEM((count+1)*sizeof(int));
		lis_free2(2,iw,*row);
		return LIS_OUT_OF_MEMORY;
	}
	(*row)[0] = (*col)[0] = 0;
	for(i=0;i<count;i++)
	{
		(*row)[i+1] = (*row)[i] + iw[i];
		(*col)[i+1] = (*col)[i] + iw[i];
	}

	lis_free(iw);
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
#endif

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_crs2vbr"
int lis_matrix_convert_crs2vbr(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	int			i,j,k,n;
	int			ii,jj,kk,ret;
	int			bnnz,bj,bnr,jpos,nnz,ij,kv,bi;
	int			err;
	int			gn,nprocs,my_rank;
	int			nr,nc;
	int			is,ie;
	int			*iw,*iw2,*count,*p2bindex;
	int			*bptr,*bindex,*ptr;
	int			*row, *col;
	LIS_SCALAR	*value;
	LIS_Comm	comm;

	LIS_DEBUG_FUNC_IN;


	nr  = Aout->conv_bnr;
	nc  = Aout->conv_bnc;
	row = Aout->conv_row;
	col = Aout->conv_col;
	if( row==NULL || col==NULL )
	{
		lis_matrix_sort_crs(Ain);
		err = lis_matrix_get_vbr_rowcol(Ain,&nr,&nc,&row,&col);
		if( err ) return err;
	}

	n        = Ain->n;
	gn       = Ain->gn;
	nprocs   = Ain->nprocs;
	my_rank  = Ain->my_rank;
	comm     = Ain->comm;
	is       = Ain->is;
	ie       = Ain->ie;

	ptr      = NULL;
	value    = NULL;
	bptr     = NULL;
	bindex   = NULL;
	iw       = NULL;
	iw2      = NULL;
	count    = NULL;
	p2bindex = NULL;

	bptr = (int *)lis_malloc( (nr+1)*sizeof(int),"lis_matrix_convert_crs2vbr::bptr" );
	if( bptr==NULL )
	{
		lis_free2(6,ptr,value,bptr,bindex,count,p2bindex);
		LIS_SETERR_MEM((nr+1)*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	p2bindex   = (int *)lis_malloc( gn*sizeof(int),"lis_matrix_convert_crs2vbr::p2bindex" );
	if( p2bindex==NULL )
	{
		lis_free2(6,ptr,value,bptr,bindex,count,p2bindex);
		LIS_SETERR_MEM(gn*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	count   = (int *)lis_malloc( (nr+1)*sizeof(int),"lis_matrix_convert_crs2vbr::count" );
	if( count==NULL )
	{
		lis_free2(6,ptr,value,bptr,bindex,count,p2bindex);
		LIS_SETERR_MEM((nr+1)*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}

	#ifdef _OPENMP
	#pragma omp parallel for private(i,j)
	#endif
	for(i=0;i<nc;i++)
	{
		for(j=col[i];j<col[i+1];j++)
		{
			p2bindex[j] = i;
		}
	}
	#ifdef _OPENMP
	#pragma omp parallel private(i,bnr,k,ii,j,bj,kk,ij,jj,iw,iw2,kv,jpos)
	#endif
	{
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<nr+1;i++) count[i] = 0;
		iw    = (int *)lis_malloc( nc*sizeof(int),"lis_matrix_convert_crs2vbr::iw" );
		iw2   = (int *)lis_malloc( nc*sizeof(int),"lis_matrix_convert_crs2vbr::iw2" );
		memset(iw,0,nc*sizeof(int));
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<nr;i++)
		{
			k    = 0;
			kk   = row[i];
			bnr  = row[i+1] - row[i];
			jj   = 0;
			for(ii=0;ii+kk<n&&ii<bnr;ii++)
			{
				for(j=Ain->ptr[kk+ii];j<Ain->ptr[kk+ii+1];j++)
				{
					bj = p2bindex[Ain->index[j]];
					jpos = iw[bj];
					if( jpos==0 )
					{
						iw[bj] = 1;
						iw2[jj] = bj;
						jj++;
					}
				}
			}
			for(bj=0;bj<jj;bj++)
			{
				k++;
				ii = iw2[bj];
				iw[ii]=0;
				count[i+1] += bnr*(col[ii+1]-col[ii]);
			}
			bptr[i+1] = k;
		}
		lis_free(iw);
		lis_free(iw2);
	}

	bptr[0] = 0;
	for(i=0;i<nr;i++)
	{
		bptr[i+1] += bptr[i];
	}
	bnnz = bptr[nr];
	for(i=0;i<nr;i++)
	{
		count[i+1] += count[i];
	}
	nnz  = count[nr];

	ptr = (int *)lis_malloc( (bnnz+1)*sizeof(int),"lis_matrix_convert_crs2vbr::ptr" );
	if( ptr==NULL )
	{
		lis_free2(6,ptr,value,bptr,bindex,count,p2bindex);
		LIS_SETERR_MEM((bnnz+1)*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	bindex = (int *)lis_malloc( bnnz*sizeof(int),"lis_matrix_convert_crs2vbr::bindex" );
	if( bindex==NULL )
	{
		lis_free2(6,ptr,value,bptr,bindex,count,p2bindex);
		LIS_SETERR_MEM(bnnz*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	value = (LIS_SCALAR *)lis_malloc( nnz*sizeof(LIS_SCALAR),"lis_matrix_convert_crs2vbr::value" );
	if( value==NULL )
	{
		lis_free2(6,ptr,value,bptr,bindex,count,p2bindex);
		LIS_SETERR_MEM(nnz*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}

	/* convert vbr */
	#ifdef _OPENMP
	#pragma omp parallel private(bi,i,ii,k,j,bj,jpos,kv,kk,ij,jj,iw,bnr,ret)
	#endif
	{
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(i=0;i<nr;i++)
		{
			j      = bptr[i];
			ptr[j] = count[i];
		}
		iw = (int *)lis_malloc( nc*sizeof(int),"lis_matrix_convert_crs2vbr::iw" );
		memset(iw,0,nc*sizeof(int));

		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(bi=0;bi<nr;bi++)
		{
			i    = row[bi];
			ii   = 0;
			kk   = bptr[bi];
			kv   = ptr[kk];
			bnr  = row[bi+1] - row[bi];
			while( i+ii<n && ii<bnr )
			{
				for( k=Ain->ptr[i+ii];k<Ain->ptr[i+ii+1];k++)
				{
					bj = p2bindex[Ain->index[k]];
					j  = Ain->index[k] - col[bj];
					jpos = iw[bj];
					if( jpos==0 )
					{
						ret             = bnr * (col[bj+1]-col[bj]);
						ij              = j*bnr + ii;
						memset(&value[kv], 0, ret*sizeof(LIS_SCALAR));
						bindex[kk]      = bj;
						value[kv+ij]    = Ain->value[k];
						iw[bj]          = kv+1;
						kv             += ret;
						ptr[kk+1]       = kv;
						kk              = kk+1;
					}
					else
					{
						ij               = j*bnr + ii;
						value[jpos+ij-1] = Ain->value[k];
					}
				}
				ii = ii+1;
			}
			for(j=bptr[bi];j<bptr[bi+1];j++)
			{
				iw[bindex[j]] = 0;
			}
		}
		lis_free(iw);
	}

	err = lis_matrix_set_vbr(nnz,nr,nc,bnnz,row,col,ptr,bptr,bindex,value,Aout);
	if( err )
	{
		lis_free2(6,ptr,value,bptr,bindex,count,p2bindex);
		return err;
	}
	err = lis_matrix_assemble(Aout);
	if( err )
	{
		lis_free2(2,count,p2bindex);
		lis_matrix_storage_destroy(Aout);
		return err;
	}

	lis_free2(2,count,p2bindex);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_vbr2crs"
int lis_matrix_convert_vbr2crs(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
	int			i,j,k,l;
	int			nr,nc,bnr,bnc,bi,bj;
	int			err;
	int			n,nnz,is;
	int			*ptr,*index;
	LIS_SCALAR	*value;


	n       = Ain->n;
	nr      = Ain->nr;
	nc      = Ain->nc;
	is      = Ain->is;

	ptr     = NULL;
	index   = NULL;
	value   = NULL;

	ptr = (int *)lis_malloc( (n+1)*sizeof(int),"lis_matrix_convert_vbr2crs::ptr" );
	if( ptr==NULL )
	{
		LIS_SETERR_MEM((n+1)*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}

	#ifdef _OPENMP
	#pragma omp parallel private(i,j,k,bi,bj,bnr,bnc)
	#endif
	{
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(bi=0;bi<nr;bi++)
		{
			k   = Ain->row[bi];
			bnr = Ain->row[bi+1]-Ain->row[bi];
			for(i=0;i<bnr;i++)
			{
				ptr[k+i+1] = 0;
			}
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(bi=0;bi<nr;bi++)
		{
			k   = Ain->row[bi];
			bnr = Ain->row[bi+1]-Ain->row[bi];
			for(bj=Ain->bptr[bi];bj<Ain->bptr[bi+1];bj++)
			{
				bnc = Ain->col[Ain->bindex[bj]+1] - Ain->col[Ain->bindex[bj]];
				for(j=0;j<bnc;j++)
				{
					for(i=0;i<bnr;i++)
					{
						if( Ain->value[Ain->ptr[bj] + j*bnr + i] != (LIS_SCALAR)0.0 )
						{
							ptr[k+i+1]++;
						}
					}
				}
			}
		}
	}

	ptr[0] = 0;
	for(i=0;i<n;i++)
	{
		ptr[i+1] += ptr[i];
	}
	nnz = ptr[n];

	index = (int *)lis_malloc( nnz*sizeof(int),"lis_matrix_convert_vbr2crs::index" );
	if( index==NULL )
	{
		lis_free2(3,ptr,index,value);
		LIS_SETERR_MEM(nnz*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	value = (LIS_SCALAR *)lis_malloc( nnz*sizeof(LIS_SCALAR),"lis_matrix_convert_vbr2crs::value" );
	if( value==NULL )
	{
		lis_free2(3,ptr,index,value);
		LIS_SETERR_MEM(nnz*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}

	/* convert crs */
	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,k,l,bi,bj,bnr,bnc)
	#endif
	for(bi=0;bi<nr;bi++)
	{
		l   = Ain->row[bi];
		bnr = Ain->row[bi+1]-Ain->row[bi];
		for(i=0;i<bnr;i++)
		{
			k = ptr[l+i];
			for(bj=Ain->bptr[bi];bj<Ain->bptr[bi+1];bj++)
			{
				bnc = Ain->col[Ain->bindex[bj]+1] - Ain->col[Ain->bindex[bj]];
				for(j=0;j<bnc;j++)
				{
					if( Ain->value[Ain->ptr[bj] + j*bnr + i] != (LIS_SCALAR)0.0 )
					{
						value[k]   = Ain->value[Ain->ptr[bj] + j*bnr + i];
						index[k]   = Ain->col[Ain->bindex[bj]]+j;
						k++;
					}
				}
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
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_diagonal_vbr"
int lis_matrix_get_diagonal_vbr(LIS_MATRIX A, LIS_SCALAR d[])
{
	int i,j,k,bi,bj,bjj,nr,nc;
	int bnr,bnc;
	int n;

	LIS_DEBUG_FUNC_IN;


	n   = A->n;
	nr  = A->nr;
	nc  = A->nc;

	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j,bnr)
		#endif
		for(i=0;i<nr;i++)
		{
			bnr = A->D->bns[i];
			for(j=0;j<bnr;j++)
			{
				d[A->L->row[i]+j] = A->D->v_value[i][j*bnr+j];
			}

		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(bi,bj,bjj,bnr,bnc,i,j,k)
		#endif
		for(bi=0;bi<nr;bi++)
		{
			k = 0;
			i = A->row[bi];
			bnr = A->row[bi+1] - A->row[bi];
			for(bj=A->bptr[bi];bj<A->bptr[bi+1];bj++)
			{
				bjj = A->bindex[bj];
				bnc = A->col[bjj+1] - A->col[bjj];
				if( i>=bjj*bnc && i<(bjj+1)*bnc )
				{
					for(j=i%bnc;j<bnc&&k<bnr&&i<n;j++)
					{
						d[i] = A->value[A->ptr[bj] + j*bnr + k];
						i++;
						k++;
					}
				}
				if( k==bnr ) break;
			}
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scaling_vbr"
int lis_matrix_scaling_vbr(LIS_MATRIX A, LIS_SCALAR d[])
{
	int i,j,k;
	int bi,bj;
	int nr,nc;
	int n;

	LIS_DEBUG_FUNC_IN;

	n   = A->n;
	nr  = A->nr;
	nc  = A->nc;

	if( A->is_splited )
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(bi,bj,i,j,k)
		#endif
		for(bi=0;bi<nr;bi++)
		{
			k = A->L->ptr[A->L->bptr[bi]];
			for(bj=A->L->bptr[bi];bj<A->L->bptr[bi+1];bj++)
			{
				for(j=A->L->col[A->bindex[bj]];j<A->L->col[A->bindex[bj]+1];j++)
				{
					for(i=A->L->row[bi];i<A->L->row[bi+1];i++)
					{
						A->L->value[k] *= d[i];
						k++;
					}
				}
			}
			k = A->U->ptr[A->U->bptr[bi]];
			for(bj=A->U->bptr[bi];bj<A->U->bptr[bi+1];bj++)
			{
				for(j=A->U->col[A->U->bindex[bj]];j<A->U->col[A->U->bindex[bj]+1];j++)
				{
					for(i=A->U->row[bi];i<A->U->row[bi+1];i++)
					{
						A->U->value[k] *= d[i];
						k++;
					}
				}
			}
			k = 0;
			for(j=A->U->col[bi];j<A->U->col[bi+1];j++)
			{
				for(i=A->U->row[bi];i<A->U->row[bi+1];i++)
				{
					A->D->v_value[bi][k] *= d[i];
					k++;
				}
			}
		}
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(bi,bj,i,j,k)
		#endif
		for(bi=0;bi<nr;bi++)
		{
			k = A->ptr[A->bptr[bi]];
			for(bj=A->bptr[bi];bj<A->bptr[bi+1];bj++)
			{
				for(j=A->col[A->bindex[bj]];j<A->col[A->bindex[bj]+1];j++)
				{
					for(i=A->row[bi];i<A->row[bi+1];i++)
					{
						A->value[k] *= d[i];
						k++;
					}
				}
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scaling_symm_vbr"
int lis_matrix_scaling_symm_vbr(LIS_MATRIX A, LIS_SCALAR d[])
{
	int i,j,k;
	int bi,bj;
	int nr,nc;
	int n;

	LIS_DEBUG_FUNC_IN;

	n   = A->n;
	nr  = A->nr;
	nc  = A->nc;

	if( A->is_splited )
	{
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
	}
	else
	{
		#ifdef _OPENMP
		#pragma omp parallel for private(bi,bj,i,j,k)
		#endif
		for(bi=0;bi<nr;bi++)
		{
			k = A->ptr[A->bptr[bi]];
			for(bj=A->bptr[bi];bj<A->bptr[bi+1];bj++)
			{
				for(j=A->col[A->bindex[bj]];j<A->col[A->bindex[bj]+1];j++)
				{
					for(i=A->row[bi];i<A->row[bi+1];i++)
					{
						A->value[k] = A->value[k]*d[i]*d[j];
						k++;
					}
				}
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_matrix_split_vbr"
int lis_matrix_split_vbr(LIS_MATRIX A)
{
	int				i,j,jj,n;
	int				nr,nc,bs;
	int				nnzl,nnzu,bnnzl,bnnzu;
	int				err;
	int				*lrow,*lcol,*lptr,*lbptr,*lbindex;
	int				*urow,*ucol,*uptr,*ubptr,*ubindex;
	LIS_SCALAR		*lvalue,*uvalue;
	LIS_MATRIX_DIAG	D;
	#ifdef _OPENMP
		int			ku,kl,kbu,kbl;
		int			*liw,*uiw,*liw2,*uiw2;
	#endif

	LIS_DEBUG_FUNC_IN;

	n        = A->n;
	nr       = A->nr;
	nc       = A->nc;
	nnzl     = 0;
	nnzu     = 0;
	bnnzl    = 0;
	bnnzu    = 0;
	D        = NULL;
	lrow     = NULL;
	lcol     = NULL;
	lptr     = NULL;
	lbptr    = NULL;
	lbindex  = NULL;
	lvalue   = NULL;
	urow     = NULL;
	ucol     = NULL;
	uptr     = NULL;
	ubptr    = NULL;
	ubindex  = NULL;
	uvalue   = NULL;

	#ifdef _OPENMP
		liw = (int *)lis_malloc((nr+1)*sizeof(int),"lis_matrix_split_vbr::liw");
		if( liw==NULL )
		{
			LIS_SETERR_MEM((nr+1)*sizeof(int));
			return LIS_OUT_OF_MEMORY;
		}
		uiw = (int *)lis_malloc((nr+1)*sizeof(int),"lis_matrix_split_vbr::uiw");
		if( uiw==NULL )
		{
			LIS_SETERR_MEM((nr+1)*sizeof(int));
			return LIS_OUT_OF_MEMORY;
		}
		liw2 = (int *)lis_malloc((nr+1)*sizeof(int),"lis_matrix_split_vbr::liw2");
		if( liw2==NULL )
		{
			LIS_SETERR_MEM((nr+1)*sizeof(int));
			return LIS_OUT_OF_MEMORY;
		}
		uiw2 = (int *)lis_malloc((nr+1)*sizeof(int),"lis_matrix_split_vbr::uiw2");
		if( uiw2==NULL )
		{
			LIS_SETERR_MEM((nr+1)*sizeof(int));
			return LIS_OUT_OF_MEMORY;
		}

		#pragma omp parallel for private(i)
		for(i=0;i<nr+1;i++)
		{
			liw[i]  = 0;
			uiw[i]  = 0;
			liw2[i] = 0;
			uiw2[i] = 0;
		}
		#pragma omp parallel for private(i,j,jj)
		for(i=0;i<nr;i++)
		{
			for(j=A->bptr[i];j<A->bptr[i+1];j++)
			{
				jj = A->bindex[j];
				if( jj<i )
				{
					liw[i+1]++;
					liw2[i+1] += (A->row[i+1]-A->row[i]) * (A->col[jj+1]-A->col[jj]);
				}
				else if( jj>i )
				{
					uiw[i+1]++;
					uiw2[i+1] += (A->row[i+1]-A->row[i]) * (A->col[jj+1]-A->col[jj]);
				}
			}
		}
		for(i=0;i<nr;i++)
		{
			liw[i+1] += liw[i];
			uiw[i+1] += uiw[i];
			liw2[i+1] += liw2[i];
			uiw2[i+1] += uiw2[i];
		}
		bnnzl  = liw[nr];
		bnnzu  = uiw[nr];
		nnzl = liw2[nr];
		nnzu = uiw2[nr];
	#else
		for(i=0;i<nr;i++)
		{
			for(j=A->bptr[i];j<A->bptr[i+1];j++)
			{
				jj = A->bindex[j];
				if( jj<i )
				{
					nnzl++;
					bnnzl += (A->row[i+1]-A->row[i]) * (A->col[jj+1]-A->col[jj]);
				}
				else if( jj>i )
				{
					nnzu++;
					bnnzu += (A->row[i+1]-A->row[i]) * (A->col[jj+1]-A->col[jj]);
				}
			}
		}
	#endif

	err = lis_matrix_LU_create(A);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_vbr(n,nnzl,nr,nc,bnnzl,&lrow,&lcol,&lptr,&lbptr,&lbindex,&lvalue);
	if( err )
	{
		return err;
	}
	err = lis_matrix_malloc_vbr(n,nnzu,nr,nc,bnnzu,&urow,&ucol,&uptr,&ubptr,&ubindex,&uvalue);
	if( err )
	{
		lis_free2(6,lptr,lbindex,lvalue,uptr,ubindex,uvalue);
		return err;
	}
	err = lis_matrix_diag_duplicateM(A,&D);
	if( err )
	{
		lis_free2(6,lptr,lbindex,lvalue,uptr,ubindex,uvalue);
		return err;
	}

	#ifdef _OPENMP
		#pragma omp parallel for private(i)
		for(i=0;i<nr+1;i++)
		{
			lrow[i] = A->row[i];
			urow[i] = A->row[i];
		}
		#pragma omp parallel for private(i)
		for(i=0;i<nc+1;i++)
		{
			lcol[i] = A->col[i];
			ucol[i] = A->col[i];
		}
		#pragma omp parallel for private(i)
		for(i=0;i<nr+1;i++)
		{
			lbptr[i] = liw[i];
			ubptr[i] = uiw[i];
		}
		#pragma omp parallel for private(i)
		for(i=0;i<nr;i++)
		{
			lptr[lbptr[i]] = liw2[i];
			uptr[ubptr[i]] = uiw2[i];
		}
		#pragma omp parallel for private(i,j,kl,ku,kbl,kbu,jj,bs)
		for(i=0;i<nr;i++)
		{
			kbl = lbptr[i];
			kbu = ubptr[i];
			kl  = liw2[i];
			ku  = uiw2[i];
			for(j=A->bptr[i];j<A->bptr[i+1];j++)
			{
				jj = A->bindex[j];
				if( jj<i )
				{
					lbindex[kbl]  = jj;
					bs            = (A->row[i+1]-A->row[i]) * (A->col[jj+1]-A->col[jj]);
					lptr[kbl+1]   = lptr[kbl] + bs;
					memcpy(&lvalue[kl],&A->value[A->ptr[j]],bs*sizeof(LIS_SCALAR));;
					kbl++;
					kl           += bs;
				}
				else if( jj>i )
				{
					ubindex[kbu]  = jj;
					bs            = (A->row[i+1]-A->row[i]) * (A->col[jj+1]-A->col[jj]);
					uptr[kbu+1]   = uptr[kbu] + bs;
					memcpy(&uvalue[ku],&A->value[A->ptr[j]],bs*sizeof(LIS_SCALAR));;
					kbu++;
					ku           += bs;
				}
				else
				{
					bs            = (A->row[i+1]-A->row[i]) * (A->col[jj+1]-A->col[jj]);
					memcpy(D->v_value[i],&A->value[A->ptr[j]],bs*sizeof(LIS_SCALAR));
				}
			}
		}
		lis_free2(4,liw,uiw,liw2,uiw2);
	#else
		for(i=0;i<nr+1;i++)
		{
			lrow[i] = A->row[i];
			urow[i] = A->row[i];
		}
		for(i=0;i<nc+1;i++)
		{
			lcol[i] = A->col[i];
			ucol[i] = A->col[i];
		}
		nnzl  = 0;
		nnzu  = 0;
		bnnzl = 0;
		bnnzu = 0;
		lptr[0]  = 0;
		uptr[0]  = 0;
		lbptr[0] = 0;
		ubptr[0] = 0;
		for(i=0;i<nr;i++)
		{
			for(j=A->bptr[i];j<A->bptr[i+1];j++)
			{
				jj = A->bindex[j];
				if( jj<i )
				{
					lbindex[bnnzl] = jj;
					bs             = (A->row[i+1]-A->row[i]) * (A->col[jj+1]-A->col[jj]);
					memcpy(&lvalue[nnzl],&A->value[A->ptr[j]],bs*sizeof(LIS_SCALAR));;
					nnzl          += bs;
					bnnzl++;
					lptr[bnnzl]    = nnzl;
				}
				else if( jj>i )
				{
					ubindex[bnnzu] = jj;
					bs             = (A->row[i+1]-A->row[i]) * (A->col[jj+1]-A->col[jj]);
					memcpy(&uvalue[nnzu],&A->value[A->ptr[j]],bs*sizeof(LIS_SCALAR));;
					nnzu          += bs;
					bnnzu++;
					uptr[bnnzu]    = nnzu;
				}
				else
				{
					bs            = (A->row[i+1]-A->row[i]) * (A->col[jj+1]-A->col[jj]);
					memcpy(D->v_value[i],&A->value[A->ptr[j]],bs*sizeof(LIS_SCALAR));
				}
			}
			lbptr[i+1] = bnnzl;
			ubptr[i+1] = bnnzu;
		}
	#endif
	A->L->nr      = nr;
	A->L->nc      = nc;
	A->L->nnz     = nnzl;
	A->L->bnnz    = bnnzl;
	A->L->ptr     = lptr;
	A->L->row     = lrow;
	A->L->col     = lcol;
	A->L->bptr    = lbptr;
	A->L->bindex  = lbindex;
	A->L->value   = lvalue;

	A->U->nr      = nr;
	A->U->nc      = nc;
	A->U->nnz     = nnzu;
	A->U->bnnz    = bnnzu;
	A->U->ptr     = uptr;
	A->U->row     = urow;
	A->U->col     = ucol;
	A->U->bptr    = ubptr;
	A->U->bindex  = ubindex;
	A->U->value   = uvalue;

	A->D          = D;
	A->is_splited = LIS_TRUE;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_merge_vbr"
int lis_matrix_merge_vbr(LIS_MATRIX A)
{
	int				i,j,jj,n,nnz;
	int				bnnz,nr,nc,bs;
	int				err;
	int				*row,*col,*ptr,*bptr,*bindex;
	LIS_SCALAR		*value;

	LIS_DEBUG_FUNC_IN;


	n       = A->n;
	nr      = A->nr;
	nc      = A->nc;
	nnz     = A->nnz;
	row     = NULL;
	col     = NULL;
	ptr     = NULL;
	bptr    = NULL;
	bindex  = NULL;
	value   = NULL;
	bnnz    = A->L->bnnz + A->U->bnnz + nr;

	err = lis_matrix_malloc_vbr(n,nnz,nr,nc,bnnz,&row,&col,&ptr,&bptr,&bindex,&value);
	if( err )
	{
		return err;
	}

	bnnz    = 0;
	nnz     = 0;
	bptr[0] = 0;
	ptr[0]  = 0;
	for(i=0;i<nr+1;i++)
	{
		row[i] = A->L->row[i];
	}
	for(i=0;i<nc+1;i++)
	{
		col[i] = A->L->col[i];
	}
	for(i=0;i<nr;i++)
	{
		for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
		{
			jj             = A->L->bindex[j];
			bindex[bnnz]   = jj;
			bs             = (A->L->row[i+1]-A->L->row[i]) * (A->L->col[jj+1]-A->L->col[jj]);
			memcpy(&value[nnz],&A->L->value[A->L->ptr[j]],bs*sizeof(LIS_SCALAR));
			bnnz++;
			nnz           += bs;
			ptr[bnnz]      = nnz;
		}
		bindex[bnnz] = i;
		bs           = A->D->bns[i] * A->D->bns[i];
		memcpy(&value[nnz],A->D->v_value[i],bs*sizeof(LIS_SCALAR));
		bnnz++;
		nnz         += bs;
		ptr[bnnz]    = nnz;
		for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
		{
			jj             = A->U->bindex[j];
			bindex[bnnz]   = jj;
			bs             = (A->U->row[i+1]-A->U->row[i]) * (A->U->col[jj+1]-A->U->col[jj]);
			memcpy(&value[nnz],&A->U->value[A->U->ptr[j]],bs*sizeof(LIS_SCALAR));
			bnnz++;
			nnz           += bs;
			ptr[bnnz]      = nnz;
		}
		bptr[i+1] = bnnz;
	}

	A->bnnz       = bnnz;
	A->ptr        = ptr;
	A->row        = row;
	A->col        = col;
	A->bptr       = bptr;
	A->value      = value;
	A->bindex      = bindex;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_solve_vbr"
int lis_matrix_solve_vbr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag)
{
	int i,j,k,ii,jj,nr,bnr,bnc,bs,dim,sz;
	LIS_SCALAR	t0,t1,t2;
	LIS_SCALAR	*b,*x,w[1024];

	LIS_DEBUG_FUNC_IN;

	nr  = A->nr;
	b   = B->value;
	x   = X->value;

	
	switch(flag)
	{
	case LIS_MATRIX_LOWER:
		lis_vector_copy(B,X);
		for(i=0;i<nr;i++)
		{
			dim = A->L->row[i+1] - A->L->row[i];
			bnr = A->L->row[i];
			for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
			{
				jj     = A->L->bindex[j];
				sz     = A->L->col[jj+1] - A->L->col[jj];
				lis_array_matvec2(dim,sz,&A->L->value[A->L->ptr[j]],dim,&x[A->L->col[jj]],&x[bnr],LIS_SUB_VALUE);
			}
			lis_array_matvec2(dim,dim,A->WD->v_value[i],dim,&x[bnr],w,LIS_INS_VALUE);
			memcpy(&x[bnr],w,dim*sizeof(LIS_SCALAR));
		}
		break;
	case LIS_MATRIX_UPPER:
		lis_vector_copy(B,X);
		for(i=nr-1;i>=0;i--)
		{
			dim = A->U->row[i+1] - A->U->row[i];
			bnr = A->U->row[i];
			for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
			{
				jj     = A->U->bindex[j];
				sz     = A->U->col[jj+1] - A->U->col[jj];
				lis_array_matvec2(dim,sz,&A->U->value[A->U->ptr[j]],dim,&x[A->U->col[jj]],&x[bnr],LIS_SUB_VALUE);
			}
			lis_array_matvec2(dim,dim,A->WD->v_value[i],dim,&x[bnr],w,LIS_INS_VALUE);
			memcpy(&x[bnr],w,dim*sizeof(LIS_SCALAR));
		}
		break;
	case LIS_MATRIX_SSOR:
		lis_vector_copy(B,X);
		for(i=0;i<nr;i++)
		{
			dim = A->L->row[i+1] - A->L->row[i];
			bnr = A->L->row[i];
			for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
			{
				jj     = A->L->bindex[j];
				sz     = A->L->col[jj+1] - A->L->col[jj];
				lis_array_matvec2(dim,sz,&A->L->value[A->L->ptr[j]],dim,&x[A->L->col[jj]],&x[bnr],LIS_SUB_VALUE);
			}
			lis_array_matvec2(dim,dim,A->WD->v_value[i],dim,&x[bnr],w,LIS_INS_VALUE);
			memcpy(&x[bnr],w,dim*sizeof(LIS_SCALAR));
		}
		for(i=nr-1;i>=0;i--)
		{
			dim = A->U->row[i+1] - A->U->row[i];
			bnr = A->U->row[i];
			memset(w,0,dim*sizeof(LIS_SCALAR));
			for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
			{
				jj     = A->U->bindex[j];
				sz     = A->U->col[jj+1] - A->U->col[jj];
				lis_array_matvec2(dim,sz,&A->U->value[A->U->ptr[j]],dim,&x[A->U->col[jj]],w,LIS_ADD_VALUE);
			}
			lis_array_matvec2(dim,dim,A->WD->v_value[i],dim,w,&x[bnr],LIS_SUB_VALUE);
		}
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_solvet_vbr"
int lis_matrix_solvet_vbr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag)
{
	int i,j,k,ii,jj,nr,bnr,bnc,bs,dim,sz;
	LIS_SCALAR	t0,t1,t2;
	LIS_SCALAR	*b,*x,w[1024];

	LIS_DEBUG_FUNC_IN;

	nr  = A->nr;
	b   = B->value;
	x   = X->value;

	switch(flag)
	{
	case LIS_MATRIX_LOWER:
		lis_vector_copy(B,X);
		for(i=0;i<nr;i++)
		{
			dim = A->U->row[i+1] - A->U->row[i];
			bnr = A->U->row[i];
			lis_array_matvec2(dim,dim,A->WD->v_value[i],dim,&x[bnr],w,LIS_INS_VALUE);
			memcpy(&x[bnr],w,dim*sizeof(LIS_SCALAR));
			for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
			{
				jj     = A->U->bindex[j];
				sz     = A->U->col[jj+1] - A->U->col[jj];
				lis_array_matvec2(dim,sz,&A->U->value[A->U->ptr[j]],dim,&x[A->U->col[jj]],&x[bnr],LIS_SUB_VALUE);
			}
		}
		break;
	case LIS_MATRIX_UPPER:
		lis_vector_copy(B,X);
		for(i=nr-1;i>=0;i--)
		{
			dim = A->L->row[i+1] - A->L->row[i];
			bnr = A->L->row[i];
			lis_array_matvec2(dim,dim,A->WD->v_value[i],dim,&x[bnr],w,LIS_INS_VALUE);
			memcpy(&x[bnr],w,dim*sizeof(LIS_SCALAR));
			for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
			{
				jj     = A->L->bindex[j];
				sz     = A->L->col[jj+1] - A->L->col[jj];
				lis_array_matvec2(dim,sz,&A->L->value[A->L->ptr[j]],dim,&x[A->L->col[jj]],&x[bnr],LIS_SUB_VALUE);
			}
		}
		break;
	case LIS_MATRIX_SSOR:
		lis_vector_copy(B,X);
		for(i=0;i<nr;i++)
		{
			dim = A->U->row[i+1] - A->U->row[i];
			bnr = A->U->row[i];
			lis_array_matvec2(dim,dim,A->WD->v_value[i],dim,&x[bnr],w,LIS_INS_VALUE);
			for(j=A->U->bptr[i];j<A->U->bptr[i+1];j++)
			{
				jj     = A->U->bindex[j];
				sz     = A->U->col[jj+1] - A->U->col[jj];
				lis_array_matvec2(dim,sz,&A->U->value[A->U->ptr[j]],dim,w,&x[A->U->col[jj]],LIS_SUB_VALUE);
			}
		}
		for(i=nr-1;i>=0;i--)
		{
			dim = A->L->row[i+1] - A->L->row[i];
			bnr = A->L->row[i];
			lis_array_matvec2(dim,dim,A->WD->v_value[i],dim,&x[bnr],w,LIS_INS_VALUE);
			memcpy(&x[bnr],w,dim*sizeof(LIS_SCALAR));
			for(j=A->L->bptr[i];j<A->L->bptr[i+1];j++)
			{
				jj     = A->L->bindex[j];
				sz     = A->L->col[jj+1] - A->L->col[jj];
				lis_array_matvec2(dim,sz,&A->L->value[A->L->ptr[j]],dim,w,&x[A->L->col[jj]],LIS_SUB_VALUE);
			}
		}
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
