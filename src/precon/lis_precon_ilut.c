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

/*
 * This subroutine is made based on ITSOL.
 *
 * http://www-users.cs.umn.edu/~saad/software/ITSOL/
 *
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

#define EBABLE_BSR	0
#undef __FUNC__
#define __FUNC__ "lis_precon_create_ilut"
int lis_precon_create_ilut(LIS_SOLVER solver, LIS_PRECON precon)
{
	#ifdef ENABLE_BSR
		int				storage,block;
	#endif
	int				err;
	LIS_MATRIX		A,B;

	LIS_DEBUG_FUNC_IN;

#ifdef ENABLE_BSR
	storage     = solver->options[LIS_OPTIONS_STORAGE];
	block       = solver->options[LIS_OPTIONS_STORAGE_BLOCK];

	if( solver->A->matrix_type!=LIS_MATRIX_BSR && storage==LIS_MATRIX_BSR )
	{
		err = lis_matrix_convert_self(solver);
		if( err ) return err;
	}
#endif

	switch( solver->A->matrix_type )
	{
	case LIS_MATRIX_CRS:
		err = lis_precon_create_ilut_crs(solver,precon);
		lis_psolve_xxx[LIS_PRECON_TYPE_ILUT]  = lis_psolve_ilut_crs;
		lis_psolvet_xxx[LIS_PRECON_TYPE_ILUT]  = lis_psolvet_ilut_crs;
		break;
#ifdef ENABLE_BSR
	case LIS_MATRIX_BSR:
		err = lis_precon_create_ilut_bsr(solver,precon);
		lis_psolve_xxx[LIS_PRECON_TYPE_ILUT]  = lis_psolve_ilut_bsr;
		lis_psolvet_xxx[LIS_PRECON_TYPE_ILUT]  = lis_psolvet_ilut_bsr;
		break;
#endif
	default:
		A = solver->A;
		err = lis_matrix_duplicate(A,&B);
		if( err ) return err;
		lis_matrix_set_type(B,LIS_MATRIX_CRS);
		err = lis_matrix_convert(A,B);
		if( err ) return err;
		solver->A = B;
		err = lis_precon_create_ilut_crs(solver,precon);
		lis_psolve_xxx[LIS_PRECON_TYPE_ILUT]  = lis_psolve_ilut_crs;
		lis_psolvet_xxx[LIS_PRECON_TYPE_ILUT]  = lis_psolvet_ilut_crs;
		lis_matrix_destroy(B);
		solver->A = A;
		break;
	}

	LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
}



#undef __FUNC__
#define __FUNC__ "lis_precon_create_ilut_crs"
int lis_precon_create_ilut_crs(LIS_SOLVER solver, LIS_PRECON precon)
{
#ifdef _OPENMP
	int				err;
	int				i,j,k,ii,jj,kk;
	int				is,ie,my_rank,nprocs;
	int				n,nr,nnz,lfil,len;
	LIS_SCALAR		gamma,t,tol,toldd,m;
	LIS_MATRIX		A;
	LIS_MATRIX_ILU	L,U;
	LIS_VECTOR		D;

	LIS_SCALAR		tnorm, tolnorm;
	LIS_SCALAR		fact,lxu,*wn,*w;
	int				lenu,lenl,col,jpos,jrow,upos,para;
	int				*jbuf,*iw;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	tol    = solver->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN];
	m      = solver->params[LIS_PARAMS_RATE-LIS_OPTIONS_LEN];
	gamma  = solver->params[LIS_PARAMS_GAMMA-LIS_OPTIONS_LEN];
	lfil   = (int)((double)A->nnz/(2.0*n))*m;
	nprocs = omp_get_max_threads();

	L      = NULL;
	U      = NULL;


	err = lis_matrix_ilu_create(n,1,&L);
	if( err ) return err;
	err = lis_matrix_ilu_create(n,1,&U);
	if( err ) return err;
	err = lis_matrix_ilu_setCR(L);
	if( err ) return err;
	err = lis_matrix_ilu_setCR(U);
	if( err ) return err;
	err = lis_vector_duplicate(A,&D);
	if( err )
	{
		return err;
	}

	w   = (LIS_SCALAR *)lis_malloc(nprocs*(n+1)*sizeof(LIS_SCALAR),"lis_precon_create_ilut_crs::w");
	if( w==NULL )
	{
		LIS_SETERR_MEM(nprocs*(n+1)*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	wn = (LIS_SCALAR *)lis_malloc(nprocs*n*sizeof(LIS_SCALAR),"lis_precon_create_ilut_crs::w");
	if( wn==NULL )
	{
		LIS_SETERR_MEM(nprocs*n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}

	jbuf   = (int *)lis_malloc(nprocs*n*sizeof(int),"lis_precon_create_ilut_crs::iw");
	if( jbuf==NULL )
	{
		LIS_SETERR_MEM(nprocs*n*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	iw   = (int *)lis_malloc(nprocs*n*sizeof(int),"lis_precon_create_ilut_crs::iw");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(nprocs*n*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}


	#pragma omp parallel private(is,ie,my_rank,i,j,k,jj,tnorm,tolnorm,len,lenu,lenl,col,t,jpos,jrow,fact,lxu,upos)
	{
		my_rank  = omp_get_thread_num();
		LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

		for(i=is;i<ie;i++) iw[my_rank*n+i] = -1;

		for(i=is;i<ie;i++)
		{
			tnorm = 0;
			k = 0;
			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				jj = A->index[j];
				if( jj<is || jj>=ie ) continue;
				tnorm += fabs(A->value[j]);
				k++;
			}
			tnorm   = tnorm / (double)k;
			tolnorm = tol * tnorm;

			lenu = 0;
			lenl = 0;
			jbuf[my_rank*n+i] = i;
			w[my_rank*n+i] = 0;
			iw[my_rank*n+i] = i;

			for(j=A->ptr[i];j<A->ptr[i+1];j++)
			{
				col = A->index[j];
				if( col<is || col>=ie ) continue;
				t = A->value[j];
				if( col < i )
				{
					jbuf[my_rank*n+lenl] = col;
					iw[my_rank*n+col] = lenl;
					w[my_rank*n+lenl] = t;
					lenl++;
				}
				else if( col == i )
				{
					w[my_rank*n+i] = t;
				}
				else
				{
					lenu++;
					jpos = i + lenu;
					jbuf[my_rank*n+jpos] = col;
					iw[my_rank*n+col] = jpos;
					w[my_rank*n+jpos] = t;
				}
			}

			j = -1;
			len = 0;

			while( ++j < lenl )
			{
				jrow = jbuf[my_rank*n+j];
				jpos = j;
				for(k=j+1;k<lenl;k++)
				{
					if( jbuf[my_rank*n+k]<jrow )
					{
						jrow = jbuf[my_rank*n+k];
						jpos = k;
					}
				}
				if( jpos!=j )
				{
					col = jbuf[my_rank*n+j];
					jbuf[my_rank*n+j] = jbuf[my_rank*n+jpos];
					jbuf[my_rank*n+jpos] = col;
					iw[my_rank*n+jrow] = j;
					iw[my_rank*n+col] = jpos;
					t = w[my_rank*n+j];
					w[my_rank*n+j] = w[my_rank*n+jpos];
					w[my_rank*n+jpos] = t;
				}
				fact = w[my_rank*n+j] * D->value[jrow];
				w[my_rank*n+j] = fact;
				iw[my_rank*n+jrow] = -1;

				for(k=0;k<U->nnz[jrow];k++)
				{
					col = U->index[jrow][k];
					jpos = iw[my_rank*n+col];
					lxu = -fact * U->value[jrow][k];

					if( fabs(lxu) < tolnorm && jpos==-1 ) continue;
					if( col >= i )
					{
						if( jpos == -1 )
						{
							lenu++;
							upos = i + lenu;
							jbuf[my_rank*n+upos] = col;
							iw[my_rank*n+col] = upos;
							w[my_rank*n+upos] = lxu;
						}
						else
						{
							w[my_rank*n+jpos] += lxu;
						}
					}
					else
					{
						if( jpos == -1 )
						{
							jbuf[my_rank*n+lenl] = col;
							iw[my_rank*n+col] = lenl;
							w[my_rank*n+lenl] = lxu;
							lenl++;
						}
						else
						{
							w[my_rank*n+jpos] += lxu;
						}
					}
				}
			}

			iw[my_rank*n+i] = -1;
			for(j=0;j<lenu;j++)
			{
				iw[ my_rank*n+jbuf[my_rank*n+i+j+1] ] = -1;
			}

			D->value[i] = 1.0 / w[my_rank*n+i];


			len = _min(lfil,lenl);
			for(j=0;j<lenl;j++)
			{
				wn[my_rank*n+j] = fabs(w[my_rank*n+j]);
				iw[my_rank*n+j] = j;
			}
			lis_sort_di(0,lenl-1,&wn[my_rank*n],&iw[my_rank*n]);
			lis_sort_i(0,len-1,&iw[my_rank*n]);
			
			L->nnz[i] = len;
			if( len>0 )
			{
				L->index[i] = (int *)malloc(len*sizeof(int));
				L->value[i] = (LIS_SCALAR *)malloc(len*sizeof(LIS_SCALAR));
			}
			for(j=0;j<len;j++)
			{
				jpos = iw[my_rank*n+j];
				L->index[i][j] = jbuf[my_rank*n+jpos];
				L->value[i][j] = w[my_rank*n+jpos];
			}
			for(j=0;j<lenl;j++) iw[my_rank*n+j] = -1;

			len = _min(lfil,lenu);
			for(j=0;j<lenu;j++)
			{
				wn[my_rank*n+j] = fabs(w[my_rank*n+i+j+1]);
				iw[my_rank*n+j] = i+j+1;
			}
			lis_sort_di(0,lenu-1,&wn[my_rank*n],&iw[my_rank*n]);
			lis_sort_i(0,len-1,&iw[my_rank*n]);
			
			U->nnz[i] = len;
			if( len>0 )
			{
				U->index[i] = (int *)malloc(len*sizeof(int));
				U->value[i] = (LIS_SCALAR *)malloc(len*sizeof(LIS_SCALAR));
			}
			for(j=0;j<len;j++)
			{
				jpos = iw[my_rank*n+j];
				U->index[i][j] = jbuf[my_rank*n+jpos];
				U->value[i][j] = w[my_rank*n+jpos];
			}
			for(j=0;j<lenu;j++) iw[my_rank*n+j] = -1;
		}
	}

	precon->L  = L;
	precon->U  = U;
	precon->D  = D;

	lis_free2(4,w,iw,wn,jbuf);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	int				err;
	int				i,j,k;
	int				n,lfil,len;
	LIS_SCALAR		gamma,t,tol,m;
	LIS_MATRIX		A;
	LIS_MATRIX_ILU	L,U;
	LIS_VECTOR		D;

	LIS_SCALAR		tnorm, tolnorm;
	LIS_SCALAR		fact,lxu,*wn,*w;
	int				lenu,lenl,col,jpos,jrow,upos;
	int				*jbuf,*iw;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	tol    = solver->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN];
	m      = solver->params[LIS_PARAMS_RATE-LIS_OPTIONS_LEN];
	gamma  = solver->params[LIS_PARAMS_GAMMA-LIS_OPTIONS_LEN];
	lfil   = (int)(((double)A->nnz/(2.0*n))*m);

	L      = NULL;
	U      = NULL;


	err = lis_matrix_ilu_create(n,1,&L);
	if( err ) return err;
	err = lis_matrix_ilu_create(n,1,&U);
	if( err ) return err;
	err = lis_matrix_ilu_setCR(L);
	if( err ) return err;
	err = lis_matrix_ilu_setCR(U);
	if( err ) return err;
	err = lis_vector_duplicate(A,&D);
	if( err )
	{
		return err;
	}

	w   = (LIS_SCALAR *)lis_malloc((n+1)*sizeof(LIS_SCALAR),"lis_precon_create_ilut_crs::w");
	if( w==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	wn = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR),"lis_precon_create_ilut_crs::w");
	if( wn==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}

	jbuf   = (int *)lis_malloc(n*sizeof(int),"lis_precon_create_ilut_crs::iw");
	if( jbuf==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	iw   = (int *)lis_malloc(n*sizeof(int),"lis_precon_create_ilut_crs::iw");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}


	for(i=0;i<n;i++) iw[i] = -1;

	for(i=0;i<n;i++)
	{
		tnorm = 0;
		for(j=A->ptr[i];j<A->ptr[i+1];j++)
		{
			tnorm += fabs(A->value[j]);
		}
		tnorm   = tnorm / (double)(A->ptr[i+1]-A->ptr[i]);
		tolnorm = tol * tnorm;

		lenu = 0;
		lenl = 0;
		jbuf[i] = i;
		w[i] = 0;
		iw[i] = i;

		for(j=A->ptr[i];j<A->ptr[i+1];j++)
		{
			col = A->index[j];
			#ifdef USE_MPI
				if( col>n-1 ) continue;
			#endif
			t = A->value[j];
			if( col < i )
			{
				jbuf[lenl] = col;
				iw[col] = lenl;
				w[lenl] = t;
				lenl++;
			}
			else if( col == i )
			{
				w[i] = t;
			}
			else
			{
				lenu++;
				jpos = i + lenu;
				jbuf[jpos] = col;
				iw[col] = jpos;
				w[jpos] = t;
			}
		}

		j = -1;
		len = 0;

		while( ++j < lenl )
		{
			jrow = jbuf[j];
			jpos = j;
			for(k=j+1;k<lenl;k++)
			{
				if( jbuf[k]<jrow )
				{
					jrow = jbuf[k];
					jpos = k;
				}
			}
			if( jpos!=j )
			{
				col = jbuf[j];
				jbuf[j] = jbuf[jpos];
				jbuf[jpos] = col;
				iw[jrow] = j;
				iw[col] = jpos;
				t = w[j];
				w[j] = w[jpos];
				w[jpos] = t;
			}
			fact = w[j] * D->value[jrow];
			w[j] = fact;
			iw[jrow] = -1;

			for(k=0;k<U->nnz[jrow];k++)
			{
				col = U->index[jrow][k];
				jpos = iw[col];
				lxu = -fact * U->value[jrow][k];

				if( fabs(lxu) < tolnorm && jpos==-1 ) continue;
				if( col >= i )
				{
					if( jpos == -1 )
					{
						lenu++;
						upos = i + lenu;
						jbuf[upos] = col;
						iw[col] = upos;
						w[upos] = lxu;
					}
					else
					{
						w[jpos] += lxu;
					}
				}
				else
				{
					if( jpos == -1 )
					{
						jbuf[lenl] = col;
						iw[col] = lenl;
						w[lenl] = lxu;
						lenl++;
					}
					else
					{
						w[jpos] += lxu;
					}
				}
			}
/*			for(kk=0;kk<bs;kk++)
			{
				w[bs*len+kk] = -buf_fact[kk];
			}
			jbuf[len] = jrow;
			len++;*/
		}

		iw[i] = -1;
		for(j=0;j<lenu;j++)
		{
			iw[ jbuf[i+j+1] ] = -1;
		}

		D->value[i] = 1.0 / w[i];


		len = _min(lfil,lenl);
		for(j=0;j<lenl;j++)
		{
			wn[j] = fabs(w[j]);
			iw[j] = j;
		}
		lis_sort_di(0,lenl-1,wn,iw);
		lis_sort_i(0,len-1,iw);
		
		L->nnz[i] = len;
		if( len>0 )
		{
			L->index[i] = (int *)malloc(len*sizeof(int));
			L->value[i] = (LIS_SCALAR *)malloc(len*sizeof(LIS_SCALAR));
		}
		for(j=0;j<len;j++)
		{
			jpos = iw[j];
			L->index[i][j] = jbuf[jpos];
			L->value[i][j] = w[jpos];
		}
		for(j=0;j<lenl;j++) iw[j] = -1;

		len = _min(lfil,lenu);
		for(j=0;j<lenu;j++)
		{
			wn[j] = fabs(w[i+j+1]);
			iw[j] = i+j+1;
		}
		lis_sort_di(0,lenu-1,wn,iw);
		lis_sort_i(0,len-1,iw);
		
		U->nnz[i] = len;
		if( len>0 )
		{
			U->index[i] = (int *)malloc(len*sizeof(int));
			U->value[i] = (LIS_SCALAR *)malloc(len*sizeof(LIS_SCALAR));
		}
		for(j=0;j<len;j++)
		{
			jpos = iw[j];
			U->index[i][j] = jbuf[jpos];
			U->value[i][j] = w[jpos];
		}
		for(j=0;j<lenu;j++) iw[j] = -1;
	}

	precon->L  = L;
	precon->U  = U;
	precon->D  = D;

	lis_free2(4,w,iw,wn,jbuf);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_psolve_ilut_crs"
int lis_psolve_ilut_crs(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
#ifdef _OPENMP
	int i,j,jj,n;
	int nprocs,my_rank,is,ie;
	LIS_SCALAR t;
	LIS_SCALAR *b,*x;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;
	LIS_PRECON	precon;
	LIS_QUAD_DECLAR;
	#ifdef USE_QUAD_PRECISION
		LIS_SCALAR *xl;
	#endif

	/*
	 *  LUx = b
	 *  LU  = (D + L*A) * (I + D^-1 * U*A)
	 */

	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->D;
	b = B->value;
	x = X->value;
	#ifdef USE_QUAD_PRECISION
		xl = X->value_lo;
	#endif
	n = solver->A->n;

	#ifdef USE_QUAD_PRECISION
		if( B->precision==LIS_PRECISION_DEFAULT )
		{
	#endif
			lis_vector_copy(B,X);
			nprocs = omp_get_max_threads();
			#pragma omp parallel private(i,j,jj,is,ie,my_rank)
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
				for(i=is; i<ie; i++)
				{
					for(j=0;j<L->nnz[i];j++)
					{
						jj    = L->index[i][j];
						x[i] -= L->value[i][j] * x[jj];
					}
				}
				for(i=ie-1; i>=is; i--)
				{
					for(j=0;j<U->nnz[i];j++)
					{
						jj    = U->index[i][j];
						x[i] -= U->value[i][j] * x[jj];
					}
					x[i] = D->value[i] * x[i];
				}
			}
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_vector_copyex_mm(B,X);
			nprocs = omp_get_max_threads();
			#ifndef USE_SSE2
				#pragma omp parallel private(i,j,jj,is,ie,my_rank,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel private(i,j,jj,is,ie,my_rank,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
			#endif
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
				for(i=is; i<ie; i++)
				{
					for(j=0;j<L->nnz[i];j++)
					{
						jj    = L->index[i][j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(x[i],xl[i],x[i],xl[i],x[jj],xl[jj],-L->value[i][j]);
						#else
							LIS_QUAD_FMAD_SSE2(x[i],xl[i],x[i],xl[i],x[jj],xl[jj],-L->value[i][j]);
						#endif
/*						x[i] -= L->value[i][j] * x[jj];*/
					}
				}
				for(i=ie-1; i>=is; i--)
				{
					for(j=0;j<U->nnz[i];j++)
					{
						jj    = U->index[i][j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(x[i],xl[i],x[i],xl[i],x[jj],xl[jj],-U->value[i][j]);
						#else
							LIS_QUAD_FMAD_SSE2(x[i],xl[i],x[i],xl[i],x[jj],xl[jj],-U->value[i][j]);
						#endif
/*						x[i] -= U->value[i][j] * x[jj];*/
					}
					#ifndef USE_SSE2
						LIS_QUAD_MULD(x[i],xl[i],x[i],xl[i],D->value[i]);
					#else
						LIS_QUAD_MULD_SSE2(x[i],xl[i],x[i],xl[i],D->value[i]);
					#endif
/*					x[i] = D->value[i] * x[i];*/
				}
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	int i,j,jj,n;
	LIS_SCALAR *b,*x;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;
	LIS_PRECON	precon;
	LIS_QUAD_DECLAR;
	#ifdef USE_QUAD_PRECISION
		LIS_SCALAR *xl;
	#endif

	/*
	 *  LUx = b
	 *  LU  = (D + L*A) * (I + D^-1 * U*A)
	 */

	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->D;
	b = B->value;
	x = X->value;
	#ifdef USE_QUAD_PRECISION
		xl = X->value_lo;
	#endif
	n = solver->A->n;

	#ifdef USE_QUAD_PRECISION
		if( B->precision==LIS_PRECISION_DEFAULT )
		{
	#endif
			lis_vector_copy(B,X);
			for(i=0; i<n; i++)
			{
				for(j=0;j<L->nnz[i];j++)
				{
					jj    = L->index[i][j];
					x[i] -= L->value[i][j] * x[jj];
				}
			}
			for(i=n-1; i>=0; i--)
			{
				for(j=0;j<U->nnz[i];j++)
				{
					jj    = U->index[i][j];
					x[i] -= U->value[i][j] * x[jj];
				}
				x[i] = D->value[i] * x[i];
			}
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_vector_copy(B,X);
			for(i=0; i<n; i++)
			{
				for(j=0;j<L->nnz[i];j++)
				{
					jj     = L->index[i][j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(x[i],xl[i],x[i],xl[i],x[jj],xl[jj],-L->value[i][j]);
					#else
						LIS_QUAD_FMAD_SSE2(x[i],xl[i],x[i],xl[i],x[jj],xl[jj],-L->value[i][j]);
					#endif
/*					x[i] -= L->value[i][j] * x[jj];*/
				}
			}
			for(i=n-1; i>=0; i--)
			{
				for(j=0;j<U->nnz[i];j++)
				{
					jj = U->index[i][j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(x[i],xl[i],x[i],xl[i],x[jj],xl[jj],-U->value[i][j]);
					#else
						LIS_QUAD_FMAD_SSE2(x[i],xl[i],x[i],xl[i],x[jj],xl[jj],-U->value[i][j]);
					#endif
/*					x[i] -= U->value[i][j] * x[jj];*/
				}
				#ifndef USE_SSE2
					LIS_QUAD_MULD(x[i],xl[i],x[i],xl[i],D->value[i]);
				#else
					LIS_QUAD_MULD_SSE2(x[i],xl[i],x[i],xl[i],D->value[i]);
				#endif
/*				x[i] = D->value[i]*x[i];*/
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_psolvet_ilut_crs"
int lis_psolvet_ilut_crs(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
#ifdef _OPENMP
	int i,j,jj,n;
	int is,ie,my_rank,nprocs;
	LIS_SCALAR *b,*x;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;
	LIS_PRECON	precon;
	LIS_QUAD_DECLAR;
	#ifdef USE_QUAD_PRECISION
		LIS_SCALAR *xl;
	#endif

	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->D;
	b = B->value;
	x = X->value;
	#ifdef USE_QUAD_PRECISION
		xl = X->value_lo;
	#endif
	n = solver->A->n;
	nprocs = omp_get_max_threads();

	#ifdef USE_QUAD_PRECISION
		if( B->precision==LIS_PRECISION_DEFAULT )
		{
	#endif
			lis_vector_copy(B,X);
			#pragma omp parallel private(i,j,jj,is,ie,my_rank)
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

				for(i=is;i<ie;i++)
				{
					x[i] = D->value[i]*x[i];
					for(j=0;j<U->nnz[i];j++)
					{
						jj     = U->index[i][j];
						x[jj] -= U->value[i][j] * x[i];
					}
				}
				for(i=ie-1;i>=is;i--)
				{
					for(j=0;j<L->nnz[i];j++)
					{
						jj     = L->index[i][j];
						x[jj] -= L->value[i][j] * x[i];
					}
				}
			}
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_vector_copyex_mm(B,X);
			nprocs = omp_get_max_threads();
			#ifndef USE_SSE2
				#pragma omp parallel private(i,j,jj,is,ie,my_rank,p1,p2,tq,bhi,blo,chi,clo,sh,sl,th,tl,eh,el)
			#else
				#pragma omp parallel private(i,j,jj,is,ie,my_rank,bh,ch,sh,wh,th,bl,cl,sl,wl,tl,p1,p2,t0,t1,t2,eh)
			#endif
			{
				my_rank = omp_get_thread_num();
				LIS_GET_ISIE(my_rank,nprocs,n,is,ie);

				for(i=is;i<ie;i++)
				{
					#ifndef USE_SSE2
						LIS_QUAD_MULD(x[i],xl[i],x[i],xl[i],D->value[i]);
					#else
						LIS_QUAD_MULD_SSE2(x[i],xl[i],x[i],xl[i],D->value[i]);
					#endif
/*					x[i] = D->value[i]*x[i];*/
					for(j=0;j<U->nnz[i];j++)
					{
						jj     = U->index[i][j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-U->value[i][j]);
						#else
							LIS_QUAD_FMAD_SSE2(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-U->value[i][j]);
						#endif
/*						x[jj] -= U->value[i][j] * x[i];*/
					}
				}
				for(i=ie-1;i>=is;i--)
				{
					for(j=0;j<L->nnz[i];j++)
					{
						jj     = L->index[i][j];
						#ifndef USE_SSE2
							LIS_QUAD_FMAD(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-L->value[i][j]);
						#else
							LIS_QUAD_FMAD_SSE2(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-L->value[i][j]);
						#endif
/*						x[jj] -= L->value[i][j] * x[i];*/
					}
				}
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	int i,j,jj,n;
	LIS_SCALAR *b,*x;
	LIS_MATRIX_ILU L,U;
	LIS_VECTOR D;
	LIS_PRECON	precon;
	LIS_QUAD_DECLAR;
	#ifdef USE_QUAD_PRECISION
		LIS_SCALAR *xl;
	#endif


	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->D;
	b = B->value;
	x = X->value;
	#ifdef USE_QUAD_PRECISION
		xl = X->value_lo;
	#endif
	n = solver->A->n;

	#ifdef USE_QUAD_PRECISION
		if( B->precision==LIS_PRECISION_DEFAULT )
		{
	#endif
			lis_vector_copy(B,X);
			for(i=0; i<n; i++)
			{
				x[i] = D->value[i]*x[i];
				for(j=0;j<U->nnz[i];j++)
				{
					jj     = U->index[i][j];
					x[jj] -= U->value[i][j] * x[i];
				}
			}
			for(i=n-1; i>=0; i--)
			{
				for(j=0;j<L->nnz[i];j++)
				{
					jj     = L->index[i][j];
					x[jj] -= L->value[i][j] * x[i];
				}
			}
	#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_vector_copy(B,X);
			for(i=0; i<n; i++)
			{
				#ifndef USE_SSE2
					LIS_QUAD_MULD(x[i],xl[i],x[i],xl[i],D->value[i]);
				#else
					LIS_QUAD_MULD_SSE2(x[i],xl[i],x[i],xl[i],D->value[i]);
				#endif
/*				x[i] = D->value[i]*x[i];*/
				for(j=0;j<U->nnz[i];j++)
				{
					jj     = U->index[i][j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-U->value[i][j]);
					#else
						LIS_QUAD_FMAD_SSE2(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-U->value[i][j]);
					#endif
/*					x[jj] -= U->value[i][j] * x[i];*/
				}
			}
			for(i=n-1; i>=0; i--)
			{
				for(j=0;j<L->nnz[i];j++)
				{
					jj     = L->index[i][j];
					#ifndef USE_SSE2
						LIS_QUAD_FMAD(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-L->value[i][j]);
					#else
						LIS_QUAD_FMAD_SSE2(x[jj],xl[jj],x[jj],xl[jj],x[i],xl[i],-L->value[i][j]);
					#endif
/*					x[jj] -= L->value[i][j] * x[i];*/
				}
			}
		}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_precon_create_ilut_bsr"
int lis_precon_create_ilut_bsr(LIS_SOLVER solver, LIS_PRECON precon)
{
	int				err;
	int				i,j,k,kk,bnr,bs;
	int				n,nr,annz,lfil,len;
	LIS_SCALAR		gamma,t,tol,m;
	LIS_MATRIX		A;
	LIS_MATRIX_ILU	L,U;
	LIS_MATRIX_DIAG	D;

	LIS_SCALAR		tnorm, tolnorm;
	LIS_SCALAR		buf_ns[16],buf_fact[16],*xnrm,*wn,*w;
	int				lenu,lenl,col,jpos,jrow,upos,para;
	int				*jbuf,*iw;

	LIS_DEBUG_FUNC_IN;


	A      = solver->A;
	n      = A->n;
	nr     = A->nr;
	bnr    = A->bnr;
	bs     = bnr*bnr;
	tol    = solver->params[LIS_PARAMS_DROP-LIS_OPTIONS_LEN];
	m      = solver->params[LIS_PARAMS_RATE-LIS_OPTIONS_LEN];
	gamma  = solver->params[LIS_PARAMS_GAMMA-LIS_OPTIONS_LEN];
	annz   = 10+A->bnnz / A->nr;
	lfil   = (int)(((double)A->bnnz/(2.0*nr))*m);

	L      = NULL;
	U      = NULL;


	err = lis_matrix_ilu_create(nr,bnr,&L);
	if( err ) return err;
	err = lis_matrix_ilu_create(nr,bnr,&U);
	if( err ) return err;
	err = lis_matrix_ilu_setCR(L);
	if( err ) return err;
	err = lis_matrix_ilu_setCR(U);
	if( err ) return err;
	err = lis_matrix_diag_duplicateM(A,&D);
	if( err )
	{
		return err;
	}

	w   = (LIS_SCALAR *)lis_malloc(bs*(nr+1)*sizeof(LIS_SCALAR),"lis_precon_create_iluc_crs::w");
	if( w==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	xnrm = (LIS_SCALAR *)lis_malloc(nr*sizeof(LIS_SCALAR),"lis_precon_create_iluc_crs::w");
	if( xnrm==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}
	wn = (LIS_SCALAR *)lis_malloc(nr*sizeof(LIS_SCALAR),"lis_precon_create_iluc_crs::w");
	if( wn==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(LIS_SCALAR));
		return LIS_OUT_OF_MEMORY;
	}

	jbuf   = (int *)lis_malloc(n*sizeof(int),"lis_precon_create_iluc_crs::iw");
	if( jbuf==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	iw   = (int *)lis_malloc(nr*sizeof(int),"lis_precon_create_iluc_crs::iw");
	if( iw==NULL )
	{
		LIS_SETERR_MEM(n*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}


	for(i=0;i<nr;i++) iw[i] = -1;

	for(i=0;i<nr;i++)
	{
		tnorm = 0;
		for(j=A->bptr[i];j<A->bptr[i+1];j++)
		{
			lis_array_nrm2(bs,&A->value[bs*j],&t);
			tnorm = _max(t,tnorm);
		}
		tolnorm = tol * tnorm;

		lenu = 1;
		lenl = 0;
		jbuf[i] = i;
		memset(&w[bs*i],0,bs*sizeof(LIS_SCALAR));
		iw[i] = i;

		for(j=A->bptr[i];j<A->bptr[i+1];j++)
		{
			col = A->bindex[j];
			lis_array_nrm2(bs,&A->value[bs*j],&t);
			if( t<tolnorm && col!=i ) continue;
			if( col < i )
			{
				jbuf[lenl] = col;
				iw[col] = lenl;
				memcpy(&w[bs*lenl],&A->value[bs*j],bs*sizeof(LIS_SCALAR));
				lenl++;
			}
			else if( col == i )
			{
				memcpy(&w[bs*i],&A->value[bs*j],bs*sizeof(LIS_SCALAR));
			}
			else
			{
				jpos = i + lenu;
				jbuf[jpos] = col;
				iw[col] = jpos;
				memcpy(&w[bs*jpos],&A->value[bs*j],bs*sizeof(LIS_SCALAR));
				lenu++;
			}
		}

		j = -1;
		len = 0;

		while( ++j < lenl )
		{
			jrow = jbuf[j];
			jpos = j;
			for(k=j+1;k<lenl;k++)
			{
				if( jbuf[k]<jrow )
				{
					jrow = jbuf[k];
					jpos = k;
				}
			}
			if( jpos!=j )
			{
				col = jbuf[j];
				jbuf[j] = jbuf[jpos];
				jbuf[jpos] = col;
				iw[jrow] = j;
				iw[col] = jpos;
				memcpy(buf_ns,&w[bs*j],bs*sizeof(LIS_SCALAR));
				memcpy(&w[bs*j],&w[bs*jpos],bs*sizeof(LIS_SCALAR));
				memcpy(&w[bs*jpos],buf_ns,bs*sizeof(LIS_SCALAR));
			}
/*			lis_array_matmat(bnr,&D->value[bs*jrow],&w[bs*j],buf_fact,LIS_INS_VALUE);*/
			lis_array_matinv(bnr,&D->value[bs*jrow],&w[bs*j],buf_fact);
			iw[jrow] = -1;

			lis_array_nrm2(bs,buf_fact,&t);
			if( t * xnrm[jrow] <= tolnorm ) continue;

			for(k=0;k<U->nnz[jrow];k++)
			{
				col = U->index[jrow][k];
				lis_array_matmat(bnr,buf_fact,&U->value[jrow][bs*k],buf_ns,LIS_INS_VALUE);
				jpos = iw[col];

				lis_array_nrm2(bs,buf_ns,&t);
				if( t < tolnorm && jpos == -1 )
				{
					continue;
				}

				if( col >= i )
				{
					if( jpos == -1 )
					{
						upos = i + lenu;
						jbuf[upos] = col;
						iw[col] = upos;
						memcpy(&w[bs*upos],buf_ns,bs*sizeof(LIS_SCALAR));
						lenu++;
					}
					else
					{
						for(kk=0;kk<bs;kk++)
						{
							w[bs*jpos+kk] += buf_ns[kk];
						}
					}
				}
				else
				{
					if( jpos == -1 )
					{
						jbuf[lenl] = col;
						iw[col] = lenl;
						memcpy(&w[bs*lenl],buf_ns,bs*sizeof(LIS_SCALAR));
						lenl++;
					}
					else
					{
						for(kk=0;kk<bs;kk++)
						{
							w[bs*jpos+kk] += buf_ns[kk];
						}
					}
				}
			}
			for(kk=0;kk<bs;kk++)
			{
				w[bs*len+kk] = -buf_fact[kk];
			}
			jbuf[len] = jrow;
			len++;
		}
		lenl = len;
		len = _min(lfil,lenl);
		for(j=0;j<lenl;j++)
		{
			lis_array_nrm2(bs,&w[bs*j],&wn[j]);
			iw[j] = j;
		}
		lis_sort_di(0,lenl-1,wn,iw);
		lis_sort_i(0,len-1,iw);
		
		L->nnz[i] = len;
		if( len>0 )
		{
			L->index[i] = (int *)malloc(len*sizeof(int));
			L->value[i] = (LIS_SCALAR *)malloc(bs*len*sizeof(LIS_SCALAR));
		}
		for(j=0;j<len;j++)
		{
			jpos = iw[j];
			L->index[i][j] = jbuf[jpos];
			memcpy(&L->value[i][bs*j],&w[bs*jpos],bs*sizeof(LIS_SCALAR));
		}
		for(j=0;j<lenl;j++) iw[j] = -1;

		len = _min(lfil,lenu);
		for(j=1;j<lenu;j++)
		{
			jpos = i+j;
			lis_array_nrm2(bs,&w[bs*jpos],&wn[j-1]);
			iw[j-1] = jpos;
		}
		para = lenu - 1;
		lis_sort_di(0,para-1,wn,iw);
		lis_sort_i(0,len-2,iw);
		
		U->nnz[i] = len-1;
		if( len>1 )
		{
			U->index[i] = (int *)malloc((len-1)*sizeof(int));
			U->value[i] = (LIS_SCALAR *)malloc(bs*(len-1)*sizeof(LIS_SCALAR));
		}
		lis_array_nrm2(bs,&w[bs*i],&t);
		for(j=0;j<len-1;j++)
		{
			jpos = iw[j];
			U->index[i][j] = jbuf[jpos];
			memcpy(&U->value[i][bs*j],&w[bs*jpos],bs*sizeof(LIS_SCALAR));
			t = _max(t,wn[j]);
		}
		for(j=0;j<lenu-1;j++) iw[j] = -1;

		xnrm[i] = t;

		memcpy(&D->value[bs*i],&w[bs*i],bs*sizeof(LIS_SCALAR));

		if( i==nr-1 )
		{
			switch(bnr)
			{
			case 2:
				if( n%2!=0 )
				{
					D->value[4*(nr-1)+3] = 1.0;
				}
				break;
			case 3:
				if( n%3==1 )
				{
					D->value[9*(nr-1)+4] = 1.0;
					D->value[9*(nr-1)+8] = 1.0;
				}
				else if( n%3==2 )
				{
					D->value[9*(nr-1)+8] = 1.0;
				}
				break;
			}
		}
/*		lis_array_invGauss(bnr,&D->value[bs*i]);*/
		lis_array_LUdecomp(bnr,&D->value[bs*i]);

		for(j=0;j<lenu;j++)
		{
			iw[ jbuf[i+j] ] = -1;
		}
	}

	precon->L  = L;
	precon->U  = U;
	precon->WD  = D;

	lis_free2(5,w,iw,xnrm,wn,jbuf);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_psolve_ilut_bsr"
int lis_psolve_ilut_bsr(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
	int i,j,jj,nr,bnr,bs;
	LIS_SCALAR w[9];
	LIS_SCALAR *b,*x;
	LIS_MATRIX_ILU L,U;
	LIS_MATRIX_DIAG D;
	LIS_PRECON	precon;

	/*
	 *  LUx = b
	 *  LU  = (D + L*A) * (I + D^-1 * U*A)
	 */

	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->WD;
	b = B->value;
	x = X->value;
	nr = solver->A->nr;
	bnr = solver->A->bnr;
	bs  = bnr*bnr;

	lis_vector_copy(B,X);
	for(i=0; i<nr; i++)
	{
		for(j=0;j<L->nnz[i];j++)
		{
			jj     = L->index[i][j];
			lis_array_matvec(bnr,&L->value[i][bs*j],&x[bnr*jj],&x[bnr*i],LIS_SUB_VALUE);
		}
	}
	for(i=nr-1; i>=0; i--)
	{
		for(j=0;j<U->nnz[i];j++)
		{
			jj = U->index[i][j];
			lis_array_matvec(bnr,&U->value[i][bs*j],&x[bnr*jj],&x[bnr*i],LIS_SUB_VALUE);
		}
/*		lis_array_matvec(bnr,&D->value[bs*i],&x[bnr*i],w,LIS_INS_VALUE);*/
		lis_array_invvec(bnr,&D->value[bs*i],&x[bnr*i],w);
		memcpy(&x[bnr*i],w,bnr*sizeof(LIS_SCALAR));
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_psolvet_ilut_bsr"
int lis_psolvet_ilut_bsr(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
	int i,j,jj,nr,bnr,bs;
	LIS_SCALAR w[9];
	LIS_SCALAR *b,*x;
	LIS_MATRIX_ILU L,U;
	LIS_MATRIX_DIAG D;
	LIS_PRECON	precon;

	/*
	 *  LUx = b
	 *  LU  = (D + L*A) * (I + D^-1 * U*A)
	 */

	LIS_DEBUG_FUNC_IN;

	precon = solver->precon;
	L = precon->L;
	U = precon->U;
	D = precon->WD;
	b = B->value;
	x = X->value;
	nr = solver->A->nr;
	bnr = solver->A->bnr;
	bs  = bnr*bnr;

	lis_vector_copy(B,X);
	for(i=0; i<nr; i++)
	{
		lis_array_invvect(bnr,&D->value[bs*i],&x[bnr*i],w);
		memcpy(&x[bnr*i],w,bnr*sizeof(LIS_SCALAR));
		for(j=0;j<U->nnz[i];j++)
		{
			jj = U->index[i][j];
			lis_array_matvect(bnr,&U->value[i][bs*j],&x[bnr*i],&x[bnr*jj],LIS_SUB_VALUE);
		}
	}
	for(i=nr-1; i>=0; i--)
	{
		for(j=0;j<L->nnz[i];j++)
		{
			jj     = L->index[i][j];
			lis_array_matvect(bnr,&L->value[i][bs*j],&x[bnr*i],&x[bnr*jj],LIS_SUB_VALUE);
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
