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

/************************************************
 * lis_matrix_g2l
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_matrix_g2l"
int lis_matrix_g2l(LIS_MATRIX A)
{
	int err;

	LIS_DEBUG_FUNC_IN;

	switch( A->matrix_type )
	{
	case LIS_MATRIX_CRS:
		err = lis_matrix_g2l_crs(A);
		break;
	case LIS_MATRIX_RCO:
		err = lis_matrix_g2l_rco(A);
		break;
/*
	case LIS_MATRIX_CCS:
		err = lis_matrix_g2l_crs(A);
		break;
	case LIS_MATRIX_MSR:
		err = lis_matrix_g2l_msr(A);
		break;
	case LIS_MATRIX_DIA:
		err = lis_matrix_g2l_dia(A);
		break;
	case LIS_MATRIX_ELL:
		err = lis_matrix_g2l_ell(A);
		break;
	case LIS_MATRIX_JDS:
		err = lis_matrix_g2l_jds(A);
		break;
	case LIS_MATRIX_BJD:
		err = lis_matrix_g2l_bjd(A);
		break;
	case LIS_MATRIX_BSR:
		err = lis_matrix_g2l_bsr(A);
		break;
	case LIS_MATRIX_BSC:
		err = lis_matrix_g2l_bsc(A);
		break;
	case LIS_MATRIX_VBR:
		err = lis_matrix_g2l_vbr(A);
		break;
	case LIS_MATRIX_DNS:
		err = lis_matrix_g2l_dns(A);
		break;
	case LIS_MATRIX_COO:
		err = lis_matrix_g2l_coo(A);
		break;
*/
	default:
		LIS_SETERR_IMP;
		return LIS_ERR_NOT_IMPLEMENTED;
	}

	LIS_DEBUG_FUNC_OUT;
	return err;
}

#if 0
#undef __FUNC__
#define __FUNC__ "lis_matrix_g2l_crs"
int lis_matrix_g2l_crs(LIS_MATRIX A)
{
	int i,j,jj,k;
	int n,nnz,gn,np,ns;
	int is,ie;
	int *l2g_map;
	int *t;
	int err;
	LIS_HASHTABLE g2l_map;

	LIS_DEBUG_FUNC_IN;

	n       = A->n;
	gn      = A->gn;
	nnz     = A->nnz;
	is      = A->is;
	ie      = A->ie;
	np      = n;
	g2l_map = NULL;
	l2g_map = NULL;

	err = lis_hashtable_create(&g2l_map);
	if( err )
	{
		lis_free(l_index);
		return err;
	}

	/* check np */
	for(i=0;i<n;i++)
	{
		for(j=A->ptr[i];j<A->ptr[i+1];j++)
		{
			jj = A->index[j];
			if( jj<is || jj>=ie )
			{
				if( lis_hashtable_search(g2l_map,jj)==NULL )
				{
					np++;
					lis_hashtable_set_value(g2l_map,jj,1);
				}
			}
		}
	}

	l2g_map = (int *)lis_malloc( (np-n)*sizeof(int) );
	if( g2l_map==NULL )
	{
		lis_hashtable_destroy(g2l_map);
		LIS_SETERR_MEM((np-n)*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}

	/* make l2g_map */
	ns = 0;
	for(i=0;i<gn;i++)
	{
		if( i<is || i>=ie )
		{
			if( lis_hashtable_search(g2l_map,i)!=NULL )
			{
				l2g_map[ns++] = i;
			}
		}
	}
/*	lis_sort_i(0,ns-1,l2g_map);*/

	/* global index => local index */
	k = n;
	lis_hashtable_clear(g2l_map);
	for(i=0;i<ns;i++)
	{
		jj          = l2g_map[i];
		lis_hashtable_set_value(g2l_map,jj,k++);
	}
	for(i=0;i<n;i++)
	{
		for(j=A->ptr[i];j<A->ptr[i+1];j++)
		{
			jj = A->index[j];
			if( jj>=is && jj<ie )
			{
				A->index[j] = jj - is;
			}
			else
			{
				A->index[j] = lis_hashtable_get_value(g2l_map,jj);
			}
		}
		#if 0
			j = A->ptr[i];
			lis_sort_iid(j,A->ptr[i+1]-1,l_index,A->index,A->value);
		#endif
	}
	A->np       = np;
	A->l2g_map  = l2g_map;

	lis_hashtable_destroy(g2l_map);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
#else
#undef __FUNC__
#define __FUNC__ "lis_matrix_g2l_crs"
int lis_matrix_g2l_crs(LIS_MATRIX A)
{
	int i,j,jj,k;
	int n,nnz,gn,np,ns;
	int is,ie;
	int *g2l_map;
	int *l2g_map;

	LIS_DEBUG_FUNC_IN;

	n       = A->n;
	gn      = A->gn;
	nnz     = A->nnz;
	is      = A->is;
	ie      = A->ie;
	np      = n;
	g2l_map = NULL;
	l2g_map = NULL;

	g2l_map = (int *)lis_malloc( gn*sizeof(int),"lis_matrix_g2l_crs::g2l_map" );
	if( g2l_map==NULL )
	{
		LIS_SETERR_MEM(gn*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}

	/* check np */
	for(i=0;i<gn;i++) g2l_map[i] = 0;
	for(i=0;i<n;i++)
	{
		for(j=A->ptr[i];j<A->ptr[i+1];j++)
		{
			jj = A->index[j];
			if( jj<is || jj>=ie )
			{
				if( g2l_map[jj]==0 )
				{
					np++;
					g2l_map[jj]  = 1;
				}
			}
		}
	}

	l2g_map = (int *)lis_malloc( (np-n)*sizeof(int),"lis_matrix_g2l_crs::l2g_map" );
	if( g2l_map==NULL )
	{
		lis_free(g2l_map);
		LIS_SETERR_MEM((np-n)*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}

	/* make l2g_map */
	ns = 0;
	for(i=0;i<gn;i++)
	{
		if( g2l_map[i]==1 )
		{
			if( i<is || i>=ie )
			{
				l2g_map[ns++] = i;
			}
		}
	}
/*	lis_sort_i(0,ns-1,l2g_map);*/

	/* global index => local index */
	k = n;
	for(i=0;i<ns;i++)
	{
		jj          = l2g_map[i];
		g2l_map[jj] = k++;
	}
	for(i=0;i<n;i++)
	{
		for(j=A->ptr[i];j<A->ptr[i+1];j++)
		{
			jj = A->index[j];
			if( jj>=is && jj<ie )
			{
				A->index[j] = jj - is;
			}
			else
			{
				A->index[j] = g2l_map[jj];
			}
		}
		#if 0
			j = A->ptr[i];
			lis_sort_iid(j,A->ptr[i+1]-1,l_index,A->index,A->value);
		#endif
	}
	A->np       = np;
	A->l2g_map  = l2g_map;

	lis_free(g2l_map);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
#endif

#undef __FUNC__
#define __FUNC__ "lis_matrix_g2l_rco"
int lis_matrix_g2l_rco(LIS_MATRIX A)
{
	int i,j,jj,k;
	int n,nnz,gn,np,ns;
	int is,ie;
	int *g2l_map;
	int *l2g_map;

	LIS_DEBUG_FUNC_IN;

	n       = A->n;
	gn      = A->gn;
	nnz     = A->nnz;
	is      = A->is;
	ie      = A->ie;
	np      = n;
	g2l_map = NULL;
	l2g_map = NULL;

	g2l_map = (int *)lis_malloc( gn*sizeof(int),"lis_matrix_g2l_rco::g2l_map" );
	if( g2l_map==NULL )
	{
		LIS_SETERR_MEM(gn*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}

	/* check np */
	for(i=0;i<gn;i++) g2l_map[i] = 0;
	for(i=0;i<n;i++)
	{
		for(j=0;j<A->w_row[i];j++)
		{
			jj = A->w_index[i][j];
			if( jj<is || jj>=ie )
			{
				if( g2l_map[jj]==0 )
				{
					np++;
					g2l_map[jj]  = 1;
				}
			}
		}
	}

	l2g_map = (int *)lis_malloc( (np-n)*sizeof(int),"lis_matrix_g2l_rco::l2g_map" );
	if( g2l_map==NULL )
	{
		lis_free(g2l_map);
		LIS_SETERR_MEM((np-n)*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}

	/* make l2g_map */
	ns = 0;
	for(i=0;i<gn;i++)
	{
		if( g2l_map[i]==1 )
		{
			if( i<is || i>=ie )
			{
				l2g_map[ns++] = i;
			}
		}
	}
/*	lis_sort_i(0,ns-1,l2g_map);*/

	/* global index => local index */
	k = n;
	for(i=0;i<ns;i++)
	{
		jj          = l2g_map[i];
		g2l_map[jj] = k++;
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<A->w_row[i];j++)
		{
			jj = A->w_index[i][j];
			if( jj>=is && jj<ie )
			{
				A->w_index[i][j] = jj - is;
			}
			else
			{
				A->w_index[i][j] = g2l_map[jj];
			}
		}
	}
	A->np       = np;
	A->l2g_map  = l2g_map;

	lis_free(g2l_map);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_commtable_destroy"
void lis_commtable_destroy(LIS_COMMTABLE table)
{
	if( table )
	{
		#ifdef USE_MPI
			if( table->export_index ) lis_free( table->export_index );
			if( table->import_index ) lis_free( table->import_index );
			if( table->export_ptr ) lis_free( table->export_ptr );
			if( table->import_ptr ) lis_free( table->import_ptr );
			if( table->req1 ) lis_free( table->req1 );
			if( table->req2 ) lis_free( table->req2 );
			if( table->sta1 ) lis_free( table->sta1 );
			if( table->sta2 ) lis_free( table->sta2 );
			if( table->neibpe ) lis_free( table->neibpe );
			#ifndef USE_OVERLAP
				if( table->ws ) lis_free( table->ws );
				if( table->wr ) lis_free( table->wr );
			#else
				if( table->ws ) MPI_Free_mem( table->ws );
				if( table->wr ) MPI_Free_mem( table->wr );
			#endif
		#endif
		lis_free(table);
	}
}

#ifdef USE_MPI
#undef __FUNC__
#define __FUNC__ "lis_commtable_duplicate"
int lis_commtable_duplicate(LIS_COMMTABLE tin, LIS_COMMTABLE *tout)
{
	int i;
	int n,np,nprocs,neibpetot,imnnz,exnnz,wssize,wrsize;
	int *export_index,*export_ptr;
	int *import_index,*import_ptr;
	int *neibpe;
	int err;
	LIS_SCALAR *ws,*wr;
	MPI_Request *req1,*req2;
	MPI_Status	*sta1,*sta2;
	MPI_Comm	comm;

	LIS_DEBUG_FUNC_IN;

	neibpetot = tin->neibpetot;
	imnnz     = tin->imnnz;
	exnnz     = tin->exnnz;
	wssize    = tin->wssize;
	wrsize    = tin->wrsize;
	comm      = tin->comm;

	export_index = NULL;
	export_ptr   = NULL;
	import_index = NULL;
	import_ptr   = NULL;
	neibpe       = NULL;
	req1         = NULL;
	sta1         = NULL;
	req2         = NULL;
	sta2         = NULL;
	err          = LIS_TRUE;
	do
	{
		*tout = (LIS_COMMTABLE)lis_malloc( sizeof(struct LIS_COMMTABLE_STRUCT),"lis_commtable_duplicate::tout" );
		if( *tout==NULL ) break;
		neibpe       = (int *)lis_malloc(neibpetot*sizeof(int),"lis_commtable_duplicate::neibpe");
		if( neibpe==NULL ) break;
		export_ptr = (int *)lis_malloc((neibpetot+1)*sizeof(int),"lis_commtable_duplicate::export_ptr");
		if( export_ptr==NULL ) break;
		import_ptr = (int *)lis_malloc((neibpetot+1)*sizeof(int),"lis_commtable_duplicate::import_ptr");
		if( import_ptr==NULL ) break;
		import_index  = (int *)lis_malloc(imnnz*sizeof(int),"lis_commtable_duplicate::import_index");
		if( import_index==NULL ) break;
		export_index  = (int *)lis_malloc(exnnz*sizeof(int),"lis_commtable_duplicate::export_index");
		if( export_index==NULL ) break;
		#ifndef USE_OVERLAP
			ws = (LIS_SCALAR *)lis_malloc( wssize*sizeof(LIS_SCALAR),"lis_commtable_duplicate::ws" );
		#else
			ws = NULL;
			MPI_Alloc_mem(wssize*sizeof(LIS_SCALAR),MPI_INFO_NULL,&ws);
		#endif
		if( ws==NULL ) break;
		#ifndef USE_OVERLAP
			wr = (LIS_SCALAR *)lis_malloc( wrsize*sizeof(LIS_SCALAR),"lis_commtable_duplicate::wr" );
		#else
			wr = NULL;
			MPI_Alloc_mem(wrsize*sizeof(LIS_SCALAR),MPI_INFO_NULL,&wr);
		#endif
		if( wr==NULL ) break;
		req1 = (MPI_Request *)lis_malloc( neibpetot*sizeof(MPI_Request),"lis_commtable_duplicate::req1" );
		if( req1==NULL ) break;
		sta1 = (MPI_Status  *)lis_malloc( neibpetot*sizeof(MPI_Status) ,"lis_commtable_duplicate::sta1");
		if( sta1==NULL ) break;
		req2 = (MPI_Request *)lis_malloc( neibpetot*sizeof(MPI_Request),"lis_commtable_duplicate::req2" );
		if( req2==NULL ) break;
		sta2 = (MPI_Status  *)lis_malloc( neibpetot*sizeof(MPI_Status),"lis_commtable_duplicate::sta2" );
		if( sta2==NULL ) break;
		err = LIS_FALSE;
	}while(LIS_FALSE);
	if( err )
	{
		if( export_index ) lis_free(export_index);
		if( export_ptr  ) lis_free(export_ptr);
		if( import_index ) lis_free(import_index);
		if( import_ptr  ) lis_free(import_ptr);
		if( neibpe ) lis_free(neibpe);
		if( req1 ) lis_free(req1);
		if( sta1 ) lis_free(sta1);
		if( req2 ) lis_free(req2);
		if( sta2 ) lis_free(sta2);
		if( tout ) lis_free(tout);
		#ifndef USE_OVERLAP
			if( ws ) lis_free(ws);
			if( wr ) lis_free(wr);
		#else
			if( ws ) MPI_Free_mem(ws);
			if( wr ) MPI_Free_mem(wr);
		#endif
		printf("out of memory\n");
		return LIS_OUT_OF_MEMORY;
	}

	for(i=0;i<neibpetot;i++)   neibpe[i]       = tin->neibpe[i];
	for(i=0;i<neibpetot+1;i++) export_ptr[i]   = tin->export_ptr[i];
	for(i=0;i<neibpetot+1;i++) import_ptr[i]   = tin->import_ptr[i];
	for(i=0;i<exnnz;i++)       export_index[i] = tin->export_index[i];
	for(i=0;i<imnnz;i++)       import_index[i] = tin->import_index[i];

	(*tout)->comm         = comm;
	(*tout)->pad          = tin->pad;
	(*tout)->neibpetot    = neibpetot;
	(*tout)->imnnz        = imnnz;
	(*tout)->exnnz        = exnnz;
	(*tout)->wssize       = wssize;
	(*tout)->wrsize       = wrsize;
	(*tout)->wr           = wr;
	(*tout)->ws           = ws;
	(*tout)->export_index = export_index;
	(*tout)->import_index = import_index;
	(*tout)->export_ptr   = export_ptr;
	(*tout)->import_ptr   = import_ptr;
	(*tout)->neibpe       = neibpe;
	(*tout)->req1         = req1;
	(*tout)->req2         = req2;
	(*tout)->sta1         = sta1;
	(*tout)->sta2         = sta2;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_commtable_duplicateM"
int lis_commtable_duplicateM(LIS_MATRIX Ain, LIS_MATRIX *Aout)
{
	int err;
	LIS_COMMTABLE	tout;

	LIS_DEBUG_FUNC_IN;

	err = lis_commtable_duplicate(Ain->commtable,&tout);
	if( err ) return err;

	(*Aout)->commtable = tout;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_commtable_create"
int lis_commtable_create(LIS_MATRIX A)
{
	int i,j,jj,k,imcount,excount;
	int n,nnz,gn,np,ns,nprocs;
	int is,ie;
	int my_rank;
	int neibpetot,n1,n2,wssize,wrsize;
	int *imneibpe,*exneibpe;
	int *export_index,*export_ptr;
	int *import_index,*import_ptr;
	int *neibpe;
	int *pe_map;
	int *l_index;
	int *g2l_map;
	int *t;
	int err;
	LIS_SCALAR *ws,*wr;
	MPI_Request *req1,*req2;
	MPI_Status	*sta1,*sta2;
	MPI_Comm	comm;
	LIS_COMMTABLE commtable;

	LIS_DEBUG_FUNC_IN;

	n       = A->n;
	gn      = A->gn;
	nnz     = A->nnz;
	is      = A->is;
	ie      = A->ie;
	nprocs  = A->nprocs;
	my_rank = A->my_rank;
	comm    = A->comm;
	np      = A->np;
	ns		= np - n;

	
	import_index  = (int *)lis_malloc(ns*sizeof(int),"lis_commtable_create::import_index");							
	if( import_index==NULL )
	{
		LIS_SETERR_MEM(ns*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}

	for(i=0;i<ns;i++)
	{
		import_index[i] = A->l2g_map[i];
	}

	imneibpe = (int *)lis_malloc(nprocs*sizeof(int),"lis_commtable_create::imneibpe");			
	if( imneibpe==NULL )
	{
		LIS_SETERR_MEM(nprocs*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	exneibpe = (int *)lis_malloc(nprocs*sizeof(int),"lis_commtable_create::exneibpe");			
	if( exneibpe==NULL )
	{
		LIS_SETERR_MEM(nprocs*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}

	/* check import neibpes */
	memset(imneibpe,0,nprocs*sizeof(int));
	k = 0;
	for(i=0;i<ns;i++)
	{
		while( import_index[i]>=A->ranges[k+1] ) k++;
		imneibpe[k]++;
	}
	MPI_Alltoall(imneibpe,1,MPI_INT,exneibpe,1,MPI_INT,comm);
																	
	/* check total neibpes */
	neibpetot = 0;
	excount   = 0;
	for(i=0;i<nprocs;i++)
	{
		if( imneibpe[i]!=0 || exneibpe[i]!=0 )
		{
			neibpetot++;
			excount += exneibpe[i];
		}
	}

	commtable = (LIS_COMMTABLE)lis_malloc( sizeof(struct LIS_COMMTABLE_STRUCT),"lis_commtable_create::commtable" );	
	if( commtable==NULL )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_COMMTABLE_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	neibpe       = (int *)lis_malloc(neibpetot*sizeof(int),"lis_commtable_create::neibpe");						
	if( neibpe==NULL )
	{
		LIS_SETERR_MEM(neibpetot*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	export_ptr = (int *)lis_malloc((neibpetot+1)*sizeof(int),"lis_commtable_create::export_ptr");					
	if( export_ptr==NULL )
	{
		LIS_SETERR_MEM((neibpetot+1)*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	import_ptr = (int *)lis_malloc((neibpetot+1)*sizeof(int),"lis_commtable_create::import_ptr");					
	if( import_ptr==NULL )
	{
		LIS_SETERR_MEM((neibpetot+1)*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	export_index  = (int *)lis_malloc(excount*sizeof(int),"lis_commtable_create::export_index");							
	if( export_index==NULL )
	{
		LIS_SETERR_MEM(excount*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	req1 = (MPI_Request *)lis_malloc( neibpetot*sizeof(MPI_Request),"lis_commtable_create::req1" );
	if( req1==NULL )
	{
		LIS_SETERR_MEM(neibpetot*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	sta1 = (MPI_Status  *)lis_malloc( neibpetot*sizeof(MPI_Status),"lis_commtable_create::sta1" );
	if( sta1==NULL )
	{
		LIS_SETERR_MEM(neibpetot*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	req2 = (MPI_Request *)lis_malloc( neibpetot*sizeof(MPI_Request),"lis_commtable_create::req2" );
	if( req2==NULL )
	{
		LIS_SETERR_MEM(neibpetot*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	sta2 = (MPI_Status  *)lis_malloc( neibpetot*sizeof(MPI_Status),"lis_commtable_create::sta2" );
	if( sta2==NULL )
	{
		LIS_SETERR_MEM(neibpetot*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}

	/* set import_ptr and export_ptr */
	n1            = 0;
	import_ptr[0] = 0;
	export_ptr[0] = 0;
	for(i=0;i<nprocs && n1<neibpetot;i++)
	{
		if( imneibpe[i]!=0 || exneibpe[i]!=0 )
		{
			neibpe[n1]       = i;			
			import_ptr[n1+1] = import_ptr[n1] + imneibpe[i];
			export_ptr[n1+1] = export_ptr[n1] + exneibpe[i];
			n1++;
		}																
	}

	/**/
	k = 0;
	for(i=0;i<ns;i++)
	{
		while( import_index[i]>=A->ranges[k+1] ) k++;
		import_index[i] = import_index[i] - A->ranges[k];
	}

	for(i=0;i<neibpetot;i++)										
	{																
		k  = neibpe[i];												
		n1 = import_ptr[i];										
		n2 = import_ptr[i+1] - import_ptr[i];					
		MPI_Isend(&import_index[n1],n2,MPI_INT,k,0,comm,&req1[i]);			
	}																
	for(i=0;i<neibpetot;i++)										
	{																
		k  = neibpe[i];												
		n1 = export_ptr[i];										
		n2 = export_ptr[i+1] - export_ptr[i];					
		MPI_Irecv(&export_index[n1],n2,MPI_INT,k,0,comm,&req2[i]);	
	}																
	MPI_Waitall(neibpetot, req1, sta1);								
	MPI_Waitall(neibpetot, req2, sta2);								

	k = n;
	for(i=0;i<ns;i++)
	{
		import_index[i] = n++;
	}

	#ifndef USE_QUAD_PRECISION
		wssize = export_ptr[neibpetot];
		wrsize = import_ptr[neibpetot];
	#else
		wssize = export_ptr[neibpetot]*2;
		wrsize = import_ptr[neibpetot]*2;
	#endif
	#ifndef USE_OVERLAP
		ws = (LIS_SCALAR *)lis_malloc( wssize*sizeof(LIS_SCALAR),"lis_commtable_create::ws" );
	#else
		ws = NULL;
		MPI_Alloc_mem(wssize*sizeof(LIS_SCALAR),MPI_INFO_NULL,&ws);
	#endif
	if( ws==NULL )
	{
		LIS_SETERR_MEM(wssize*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	#ifndef USE_OVERLAP
		wr = (LIS_SCALAR *)lis_malloc( wrsize*sizeof(LIS_SCALAR),"lis_commtable_create::wr" );
	#else
		wr = NULL;
		MPI_Alloc_mem(wrsize*sizeof(LIS_SCALAR),MPI_INFO_NULL,&wr);
	#endif
	if( wr==NULL )
	{
		LIS_SETERR_MEM(wrsize*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}

/*	
	printf("neibpetot = %d\n",neibpetot);
	printf("excount = %d\n",excount);
	printf("wssize = %d\n",wssize);
	printf("wrsize = %d\n",wrsize);
*/

	commtable->comm         = comm;
	commtable->pad          = A->pad_comm;
	commtable->neibpetot    = neibpetot;	
	commtable->imnnz        = ns;		
	commtable->exnnz        = excount;		
	commtable->wssize       = wssize;		
	commtable->wrsize       = wrsize;		
	commtable->neibpe       = neibpe;		
	commtable->export_index = export_index;	
	commtable->import_index = import_index;	
	commtable->export_ptr   = export_ptr;	
	commtable->import_ptr   = import_ptr;	
	commtable->req1         = req1;			
	commtable->req2         = req2;			
	commtable->sta1         = sta1;			
	commtable->sta2         = sta2;			
	commtable->ws           = ws;			
	commtable->wr           = wr;			
											
	A->commtable = commtable;

	lis_free2(2,imneibpe,exneibpe);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
#endif

#ifdef USE_MPI
#undef __FUNC__
#define __FUNC__ "lis_send_recv"
int lis_send_recv(LIS_COMMTABLE commtable, LIS_SCALAR x[])
{
	int neib,i,is,inum,neibpetot,k,pad;
	LIS_SCALAR *ws,*wr;
	int	*iw,err;

	LIS_DEBUG_FUNC_IN;

	neibpetot = commtable->neibpetot;
	ws        = commtable->ws;
	wr        = commtable->wr;
	pad       = commtable->pad;

#if 0
	iw = (int *)malloc( neibpetot*sizeof(int),"lis_send_recv::iw" );
	for(neib=0;neib<neibpetot;neib++)
	{
		is = commtable->index_export[neib];
		inum = commtable->index_export[neib+1] - is;
		MPI_Isend(&inum,1,MPI_INT,commtable->neibpe[neib],0,commtable->comm,&commtable->req1[neib]);
	}
	for(neib=0;neib<neibpetot;neib++)
	{
		is = commtable->index_import[neib];
		inum = commtable->index_import[neib+1] - is;
		MPI_Irecv(&iw[neib],1,MPI_INT,commtable->neibpe[neib],0,commtable->comm,&commtable->req2[neib]);
	}
	MPI_Waitall(neibpetot, commtable->req2, commtable->sta2);
	MPI_Waitall(neibpetot, commtable->req1, commtable->sta1);
	for(i=0;i<neibpetot;i++)
	{
		is = commtable->index_import[i];
		inum = commtable->index_import[i+1] - is;
		err = 0;
		if(iw[i]!=inum )
		{
			printf("i=%d im_count=%d ex_count=%d\n",i,iw[i],inum);
			err = 1;
		}
		MPI_Allreduce(&err,&is,1,MPI_INT,MPI_SUM,commtable->comm);
		if( is ) CHKERR(1);
	}
	for(i=0;i<neibpetot;i++)
	{
		is = commtable->index_import[i];
		inum = commtable->index_import[i+1] - is;
		err = 0;
		if(is+inum>commtable->wrsize )
		{
			printf("i=%d is=%d inum=%d wrsize=%d\n",i,is,inum,commtable->wrsize);
			err = 1;
		}
		MPI_Allreduce(&err,&is,1,MPI_INT,MPI_SUM,commtable->comm);
		if( is ) CHKERR(1);
	}
	for(i=0;i<neibpetot;i++)
	{
		is = commtable->index_export[i];
		inum = commtable->index_export[i+1] - is;
		err = 0;
		if(is+inum>commtable->wssize )
		{
			printf("i=%d is=%d inum=%d wssize=%d\n",i,is,inum,commtable->wssize);
			err = 1;
		}
		MPI_Allreduce(&err,&is,1,MPI_INT,MPI_SUM,commtable->comm);
		if( is ) CHKERR(1);
	}
	free(iw);
#endif
	for(neib=0;neib<neibpetot;neib++)
	{
		is = commtable->export_ptr[neib];
		inum = commtable->export_ptr[neib+1] - is;
		for(i=is;i<is+inum;i++)
		{
			ws[i] = x[commtable->export_index[i]];
		}
		MPI_Isend(&ws[is],inum,MPI_DOUBLE,commtable->neibpe[neib],0,commtable->comm,&commtable->req1[neib]);
	}
	for(neib=0;neib<neibpetot;neib++)
	{
		is = commtable->import_ptr[neib];
		inum = commtable->import_ptr[neib+1] - is;
		MPI_Irecv(&wr[is],inum,MPI_DOUBLE,commtable->neibpe[neib],0,commtable->comm,&commtable->req2[neib]);
	}
	MPI_Waitall(neibpetot, commtable->req2, commtable->sta2);

#if 0
	for(neib=0;neib<neibpetot;neib++)
	{
		is = commtable->index_import[neib];
		inum = commtable->index_import[neib+1] - is;
		MPI_Get_count(&commtable->sta2[neib],MPI_DOUBLE,&is);
		printf("sndre rank=%d pe=%d rec_count=%d act_count=%d\n",0,commtable->neibpe[neib],inum,is);
	}
		printf("x[%d]=%f\n",commtable->node_import[0],wr[0]);
#endif

#if 1
	k = commtable->import_index[0] + pad;
	for(i=commtable->import_ptr[0];i<commtable->import_ptr[neibpetot];i++)
	{
		x[k++] = wr[i];
	}
#else
	for(neib=0;neib<neibpetot;neib++)
	{
		is = commtable->import_ptr[neib];
		inum = commtable->import_ptr[neib+1] - is;
		for(i=is;i<is+inum;i++)
		{
			x[commtable->import_index[i]+pad] = wr[i];
		}
	}
#endif

	MPI_Waitall(neibpetot, commtable->req1, commtable->sta1);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_reduce"
int lis_reduce(LIS_COMMTABLE commtable, LIS_SCALAR x[])
{
	int neib,i,is,inum,neibpetot,pad;
	LIS_SCALAR *ws,*wr;

	LIS_DEBUG_FUNC_IN;

	neibpetot = commtable->neibpetot;
	ws        = commtable->ws;
	wr        = commtable->wr;
	pad       = commtable->pad;

	for(neib=0;neib<neibpetot;neib++)
	{
		is = commtable->import_ptr[neib];
		inum = commtable->import_ptr[neib+1] - is;
		for(i=is;i<is+inum;i++)
		{
			wr[i] = x[commtable->import_index[i]+pad];
		}
		MPI_Isend(&wr[is],inum,MPI_DOUBLE,commtable->neibpe[neib],0,commtable->comm,&commtable->req1[neib]);
	}
	for(neib=0;neib<neibpetot;neib++)
	{
		is = commtable->export_ptr[neib];
		inum = commtable->export_ptr[neib+1] - is;
		MPI_Irecv(&ws[is],inum,MPI_DOUBLE,commtable->neibpe[neib],0,commtable->comm,&commtable->req2[neib]);
	}
	MPI_Waitall(neibpetot, commtable->req2, commtable->sta2);
	for(neib=0;neib<neibpetot;neib++)
	{
		is = commtable->export_ptr[neib];
		inum = commtable->export_ptr[neib+1] - is;
		for(i=is;i<is+inum;i++)
		{
			x[commtable->export_index[i]] += ws[i];
		}
	}
	MPI_Waitall(neibpetot, commtable->req1, commtable->sta1);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#if 0
#undef __FUNC__
#define __FUNC__ "lis_matrix_redistribute_crs"
int lis_matrix_redistribute_crs(LIS_MATRIX Ain, LIS_MATRIX *Aout, int re_n)
{
	int			i,j,l,k,n,nnz,count;
	int			is,ie,nprocs,my_rank;
	int			rec_n,snd_n,up_n,iwsize;
	int			ris,rie;
	int			err,err2;
	int			*ptr,*index;
	int			*rec_pe,*rec_count,*iw;
	int			*snd_pe,*snd_count,*snd_is;
	LIS_SCALAR	*value;
	MPI_Request *req1,*req2;
	MPI_Status	*sta1,*sta2;
	MPI_Comm	comm;

	LIS_DEBUG_FUNC_IN;

	nprocs   = Ain->nprocs;
	my_rank  = Ain->my_rank;
	is       = Ain->is;
	ie       = Ain->ie;
	comm     = Ain->comm;

	req1      = NULL;
	req2      = NULL;
	sta1      = NULL;
	sta2      = NULL;
	snd_pe    = NULL;
	snd_count = NULL;
	snd_is    = NULL;
	rec_pe    = NULL;
	rec_count = NULL;
	iw        = NULL;
	ptr       = NULL;
	index     = NULL;
	value     = NULL;

	rec_pe = (int *)lis_malloc( nprocs*sizeof(int) );
	if( rec_pe==NULL )
	{
		lis_free2(13,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count,iw,ptr,index,value);
		SETERR("out of memory\n");
		return LIS_OUT_OF_MEMORY;
	}
	rec_count = (int *)lis_malloc( nprocs*sizeof(int) );
	if( rec_count==NULL )
	{
		lis_free2(13,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count,iw,ptr,index,value);
		SETERR("out of memory\n");
		return LIS_OUT_OF_MEMORY;
	}
	snd_pe = (int *)lis_malloc( nprocs*sizeof(int) );
	if( snd_pe==NULL )
	{
		lis_free2(13,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count,iw,ptr,index,value);
		SETERR("out of memory\n");
		return LIS_OUT_OF_MEMORY;
	}
	snd_count = (int *)lis_malloc( nprocs*sizeof(int) );
	if( snd_count==NULL )
	{
		lis_free2(13,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count,iw,ptr,index,value);
		SETERR("out of memory\n");
		return LIS_OUT_OF_MEMORY;
	}
	snd_is = (int *)lis_malloc( nprocs*sizeof(int) );
	if( snd_is==NULL )
	{
		lis_free2(13,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count,iw,ptr,index,value);
		SETERR("out of memory\n");
		return LIS_OUT_OF_MEMORY;
	}

	err = lis_matrix_create(re_n,0,Ain->comm,Aout);
	if( err )
	{
		return err;
	}
	n     = (*Aout)->n;

	rec_n = 0;
	ris   = (*Aout)->ranges[my_rank];
	rie   = (*Aout)->ranges[my_rank+1];
	/* recv */
	if( rie <= is )
	{
		for(i=my_rank-1;i>=0;i--)
		{
			if( Ain->ranges[i] < rie )
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = rie - Ain->ranges[i];
				rec_n++;
				break;
			}
		}
		for(i-=1;i>=0;i--)
		{
			if( Ain->ranges[i] <= ris )
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = Ain->ranges[i+1] - ris;
				rec_n++;
				break;
			}
			else
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = Ain->ranges[i+1] - Ain->ranges[i];
				rec_n++;
			}
		}
	}
	else if( ris < is )
	{
		for(i=my_rank-1;i>=0;i--)
		{
			if( Ain->ranges[i] <= ris )
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = Ain->ranges[i+1] - ris;
				rec_n++;
				break;
			}
			else
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = Ain->ranges[i+1] - Ain->ranges[i];
				rec_n++;
			}
		}
	}
	up_n = rec_n;
	lis_sort_ii(0,rec_n-1,rec_pe,rec_count);
	if( ris >= ie )
	{
		for(i=my_rank+1;i<nprocs;i++)
		{
			if( Ain->ranges[i+1] > ris )
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = Ain->ranges[i+1] - ris;
				rec_n++;
				break;
			}
		}
		for(i+=1;i<nprocs;i++)
		{
			if( Ain->ranges[i+1] >= rie )
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = rie - Ain->ranges[i];
				rec_n++;
				break;
			}
			else
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = Ain->ranges[i+1] - Ain->ranges[i];
				rec_n++;
			}
		}
	}
	else if( rie > ie )
	{
		for(i=my_rank+1;i<nprocs;i++)
		{
			if( Ain->ranges[i+1] >= rie )
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = rie - Ain->ranges[i];
				rec_n++;
				break;
			}
			else
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = Ain->ranges[i+1] - Ain->ranges[i];
				rec_n++;
			}
		}
	}

	/* send */
	snd_n = 0;
	nnz   = 0;
	for(i=0;i<nprocs;i++)
	{
		if( i!=my_rank )
		{
			if( (*Aout)->ranges[i] < is && (*Aout)->ranges[i+1] > is && (*Aout)->ranges[i+1] <= ie )
			{
				snd_pe[snd_n]    = i;
				snd_is[snd_n]    = 0;
				snd_count[snd_n] = (*Aout)->ranges[i+1] - is;
				snd_n++;
			}
			else if( (*Aout)->ranges[i] >=is && (*Aout)->ranges[i+1] <= ie )
			{
				snd_pe[snd_n]    = i;
				snd_is[snd_n]    = (*Aout)->ranges[i] - is;
				snd_count[snd_n] = (*Aout)->ranges[i+1] - (*Aout)->ranges[i];
				snd_n++;
			}
			else if( (*Aout)->ranges[i] < ie && (*Aout)->ranges[i+1] > ie && (*Aout)->ranges[i] >= is )
			{
				snd_pe[snd_n]    = i;
				snd_is[snd_n]    = (*Aout)->ranges[i] - is;
				snd_count[snd_n] = ie - (*Aout)->ranges[i];
				snd_n++;
			}
			else if( (*Aout)->ranges[i] < is && (*Aout)->ranges[i+1] > ie )
			{
				snd_pe[snd_n]    = i;
				snd_is[snd_n]    = 0;
				snd_count[snd_n] = ie - is;
				snd_n++;
			}
		}
	}
	iwsize = 0;
	for(i=0;i<snd_n;i++)
	{
		iwsize += snd_count[i];
	}
	req1 = (MPI_Request *)lis_malloc( snd_n*sizeof(MPI_Request) );
	if( req1==NULL )
	{
		lis_free2(13,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count,iw,ptr,index,value);
		SETERR("out of memory\n");
		return LIS_OUT_OF_MEMORY;
	}
	req2 = (MPI_Request *)lis_malloc( rec_n*sizeof(MPI_Request) );
	if( req2==NULL )
	{
		lis_free2(13,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count,iw,ptr,index,value);
		SETERR("out of memory\n");
		return LIS_OUT_OF_MEMORY;
	}
	sta1 = (MPI_Status *)lis_malloc( snd_n*sizeof(MPI_Status) );
	if( sta1==NULL )
	{
		lis_free2(13,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count,iw,ptr,index,value);
		SETERR("out of memory\n");
		return LIS_OUT_OF_MEMORY;
	}
	sta2 = (MPI_Status *)lis_malloc( rec_n*sizeof(MPI_Status) );
	if( sta2==NULL )
	{
		lis_free2(13,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count,iw,ptr,index,value);
		SETERR("out of memory\n");
		return LIS_OUT_OF_MEMORY;
	}
	ptr = (int *)lis_malloc( (n+1)*sizeof(int) );
	if( ptr==NULL )
	{
		lis_free2(13,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count,iw,ptr,index,value);
		SETERR("out of memory\n");
		return LIS_FAILS;
	}
	iw = (int *)lis_malloc( iwsize*sizeof(int) );
	if( iw==NULL )
	{
		lis_free2(13,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count,iw,ptr,index,value);
		SETERR("out of memory\n");
		return LIS_OUT_OF_MEMORY;
	}

#if 0
	MPI_Barrier(Ain->comm);
	for(k=0;k<snd_n;k++)
	{
		MPI_Isend(&snd_count[k],1,MPI_INT,snd_pe[k],0,Ain->comm,&req1[k]);
	}
	for(k=0;k<rec_n;k++)
	{
		MPI_Irecv(&iw[k],1,MPI_INT,rec_pe[k],0,Ain->comm,&req2[k]);
	}
	if( rec_n>0 ) MPI_Waitall(rec_n, req2, sta2);
	if( snd_n>0 ) MPI_Waitall(snd_n, req1, sta1);
	for(k=0;k<rec_n;k++)
	{
		MPI_Get_count(&sta2[k],MPI_INT,&count);
		err = 0;
		if(iw[k]!=rec_count[k])
		{
			printf("rank=%d snd_count=%d rec_count=%d pe=%d\n",my_rank,iw[k],rec_count[k],rec_pe[k]);
			err = 1;
		}
	}
	if(Ain->bn==9) CHKERR(1);
	MPI_Barrier(Ain->comm);
#endif
	/* exchange ptr */
	ptr[0] = 0;
	i      = 0;
	if( (*Aout)->ranges[my_rank] < is && (*Aout)->ranges[my_rank+1] > is && (*Aout)->ranges[my_rank+1] <= ie )
	{
		k = is - (*Aout)->ranges[my_rank];
		for(i=0;i<(*Aout)->ranges[my_rank+1]-is;i++)
		{
			ptr[k+i+1] = Ain->ptr[i+1] - Ain->ptr[i];
		}
	}
	else if( (*Aout)->ranges[my_rank] < ie && (*Aout)->ranges[my_rank+1] > ie && (*Aout)->ranges[my_rank] >= is )
	{
		k = (*Aout)->ranges[my_rank] - is;
		for(i=0;i<ie-(*Aout)->ranges[my_rank];i++)
		{
			ptr[i+1] = Ain->ptr[k+i+1] - Ain->ptr[k+i];
		}
	}
	else if( (*Aout)->ranges[my_rank+1] <= ie && (*Aout)->ranges[my_rank] >=is )
	{
		k = (*Aout)->ranges[my_rank]-is;
		for(i=0;i<(*Aout)->ranges[my_rank+1]-(*Aout)->ranges[my_rank];i++)
		{
			ptr[i+1] = Ain->ptr[k+i+1] - Ain->ptr[k+i];
		}
	}
	else if( (*Aout)->ranges[my_rank] < is && (*Aout)->ranges[my_rank+1] > ie )
	{
		k = is - (*Aout)->ranges[my_rank];
		for(i=0;i<ie-is;i++)
		{
			ptr[k+i+1] = Ain->ptr[i+1] - Ain->ptr[i];
		}
	}
	count = i+1;
	j = 0;
	for(k=0;k<snd_n;k++)
	{
		for(i=0;i<snd_count[k];i++)
		{
			iw[j+i] = Ain->ptr[snd_is[k]+i+1] - Ain->ptr[snd_is[k]+i];
		}
		MPI_Isend(&iw[j],snd_count[k],MPI_INT,snd_pe[k],0,Ain->comm,&req1[k]);
		j += snd_count[k];
	}
	j = 1;
	for(k=0;k<up_n;k++)
	{
		MPI_Irecv(&ptr[j],rec_count[k],MPI_INT,rec_pe[k],0,Ain->comm,&req2[k]);
		j += rec_count[k];
	}
	for(;k<rec_n;k++)
	{
		MPI_Irecv(&ptr[count],rec_count[k],MPI_INT,rec_pe[k],0,Ain->comm,&req2[k]);
		count += rec_count[k];
	}
	if( rec_n>0 ) MPI_Waitall(rec_n, req2, sta2);
	if( snd_n>0 ) MPI_Waitall(snd_n, req1, sta1);
#if 0
	for(k=0;k<rec_n;k++)
	{
		MPI_Get_count(&sta2[k],MPI_INT,&count);
		printf("ptr   rank=%d pe=%d rec_count=%d act_count=%d\n",my_rank,rec_pe[k],rec_count[k],count);
	}
#endif

	for(i=0;i<n;i++)
	{
		ptr[i+1] += ptr[i];
	}
	nnz = ptr[n];
#if 0
	MPI_Allreduce(&nnz,&k,1,MPI_INT,MPI_SUM,Ain->comm);
#endif

	index = (int *)lis_malloc( nnz*sizeof(int) );
	if( index==NULL )
	{
		lis_free2(13,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count,iw,ptr,index,value);
		SETERR("out of memory\n");
		return LIS_FAILS;
	}
	value = (LIS_SCALAR *)lis_malloc( nnz*sizeof(LIS_SCALAR) );
	if( value==NULL )
	{
		lis_free2(13,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count,iw,ptr,index,value);
		SETERR("out of memory\n");
		return LIS_FAILS;
	}

	/* exchange index and value */
	i      = 0;
	if( (*Aout)->ranges[my_rank] < is && (*Aout)->ranges[my_rank+1] > is && (*Aout)->ranges[my_rank+1] <= ie )
	{
		k = is - (*Aout)->ranges[my_rank];
		for(i=0;i<(*Aout)->ranges[my_rank+1]-is;i++)
		{
			l = ptr[k+i];
			for(j=Ain->ptr[i];j<Ain->ptr[i+1];j++)
			{
				index[l] = Ain->g_index[j];
				value[l] = Ain->value[j];
				l++;
			}
		}
	}
	else if( (*Aout)->ranges[my_rank] < ie && (*Aout)->ranges[my_rank+1] > ie && (*Aout)->ranges[my_rank] >= is )
	{
		k = (*Aout)->ranges[my_rank] - is;
		for(i=0;i<ie-(*Aout)->ranges[my_rank];i++)
		{
			l = ptr[i];
			for(j=Ain->ptr[k+i];j<Ain->ptr[k+i+1];j++)
			{
				index[l] = Ain->g_index[j];
				value[l] = Ain->value[j];
				l++;
			}
		}
	}
	else if( (*Aout)->ranges[my_rank+1] <= ie && (*Aout)->ranges[my_rank] >=is )
	{
		k = (*Aout)->ranges[my_rank]-is;
		for(i=0;i<(*Aout)->ranges[my_rank+1]-(*Aout)->ranges[my_rank];i++)
		{
			l = ptr[i];
			for(j=Ain->ptr[k+i];j<Ain->ptr[k+i+1];j++)
			{
				index[l] = Ain->g_index[j];
				value[l] = Ain->value[j];
				l++;
			}
		}
	}
	else if( (*Aout)->ranges[my_rank] < is && (*Aout)->ranges[my_rank+1] > ie )
	{
		k = is - (*Aout)->ranges[my_rank];
		for(i=0;i<ie-is;i++)
		{
			l = ptr[k+i];
			for(j=Ain->ptr[i];j<Ain->ptr[i+1];j++)
			{
				index[l] = Ain->g_index[j];
				value[l] = Ain->value[j];
				l++;
			}
		}
	}
	/* index */
	for(k=0;k<snd_n;k++)
	{
		snd_count[k] = Ain->ptr[snd_is[k]+snd_count[k]] - Ain->ptr[snd_is[k]];
		MPI_Isend(&snd_count[k],1,MPI_INT,snd_pe[k],0,Ain->comm,&req1[k]);
	}
	for(k=0;k<rec_n;k++)
	{
		MPI_Irecv(&rec_count[k],1,MPI_INT,rec_pe[k],0,Ain->comm,&req2[k]);
	}
	if( rec_n>0 ) MPI_Waitall(rec_n, req2, sta2);
	if( snd_n>0 ) MPI_Waitall(snd_n, req1, sta1);
	for(k=0;k<snd_n;k++)
	{
		j     = Ain->ptr[snd_is[k]];
		MPI_Isend(&Ain->g_index[j],snd_count[k],MPI_INT,snd_pe[k],0,Ain->comm,&req1[k]);
	}
	j = 0;
	for(k=0;k<up_n;k++)
	{
		MPI_Irecv(&index[j],rec_count[k],MPI_INT,rec_pe[k],0,Ain->comm,&req2[k]);
		j += rec_count[k];
	}
	j = l;
	for(;k<rec_n;k++)
	{
		MPI_Irecv(&index[j],rec_count[k],MPI_INT,rec_pe[k],0,Ain->comm,&req2[k]);
		j += rec_count[k];
	}
	if( rec_n>0 ) MPI_Waitall(rec_n, req2, sta2);
	if( snd_n>0 ) MPI_Waitall(snd_n, req1, sta1);
#if 0
	for(k=0;k<rec_n;k++)
	{
		MPI_Get_count(&sta2[k],MPI_INT,&count);
		printf("index rank=%d pe=%d rec_count=%d act_count=%d\n",my_rank,rec_pe[k],rec_count[k],count);
	}
#endif
	/* value */
	for(k=0;k<snd_n;k++)
	{
		j     = Ain->ptr[snd_is[k]];
		MPI_Isend(&Ain->value[j],snd_count[k],MPI_DOUBLE,snd_pe[k],0,Ain->comm,&req1[k]);
	}
	j = 0;
	for(k=0;k<up_n;k++)
	{
		MPI_Irecv(&value[j],rec_count[k],MPI_DOUBLE,rec_pe[k],0,Ain->comm,&req2[k]);
		j += rec_count[k];
	}
	j = l;
	for(;k<rec_n;k++)
	{
		MPI_Irecv(&value[j],rec_count[k],MPI_DOUBLE,rec_pe[k],0,Ain->comm,&req2[k]);
		j += rec_count[k];
	}
	if( rec_n>0 ) MPI_Waitall(rec_n, req2, sta2);
	if( snd_n>0 ) MPI_Waitall(snd_n, req1, sta1);
#if 0
	for(k=0;k<rec_n;k++)
	{
		MPI_Get_count(&sta2[k],MPI_DOUBLE,&count);
		printf("value rank=%d pe=%d rec_count=%d act_count=%d\n",my_rank,rec_pe[k],rec_count[k],count);
	}
#endif
#if 0
	is = (*Aout)->ranges[my_rank];
	ie = (*Aout)->ranges[my_rank+1];
	printf("rank=%d is=%d ie=%d ris=%d rie=%d ptr[n]=%d\n",my_rank,is,ie,ris,rie,ptr[ie-is]);
	printf("rank=%d is=%d ie=%d ris=%d rie=%d index[0]=%d index[n]=%d\n",my_rank,is,ie,ris,rie,index[0],index[ptr[ie-is]-1]);
	printf("rank=%d is=%d ie=%d ris=%d rie=%d value[0]=%f value[n]=%f\n",my_rank,is,ie,ris,rie,value[0],value[ptr[ie-is]-1]);
#endif

	lis_free2(10,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count,iw);
	err = lis_matrix_set_crs(nnz,ptr,index,value,*Aout);
	if( err )
	{
		lis_free2(3,ptr,index,value);
		lis_matrix_destroy(*Aout);
		*Aout = NULL;
		return err;
	}
	err = lis_matrix_assemble(Aout);
	if( err )
	{
		lis_free2(3,ptr,index,value);
		lis_matrix_destroy(*Aout);
		*Aout = NULL;
		return err;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}
#endif

#undef __FUNC__
#define __FUNC__ "lis_vector_redistribute"
int lis_vector_redistribute(LIS_MATRIX Ain, LIS_VECTOR vin, LIS_VECTOR *vout)
{
	int			i,j,l,k,n,nnz,count;
	int			is,ie,nprocs,my_rank;
	int			rec_n,snd_n,up_n;
	int			ris,rie;
	int			err;
	int			*rec_pe,*rec_count;
	int			*snd_pe,*snd_count,*snd_is;
	MPI_Request *req1,*req2;
	MPI_Status	*sta1,*sta2;

	LIS_DEBUG_FUNC_IN;

	nprocs   = vin->nprocs;
	my_rank  = vin->my_rank;
	is       = vin->is;
	ie       = vin->ie;

	req1      = NULL;
	req2      = NULL;
	sta1      = NULL;
	sta2      = NULL;
	snd_pe    = NULL;
	snd_count = NULL;
	snd_is    = NULL;
	rec_pe    = NULL;
	rec_count = NULL;

	rec_pe = (int *)malloc( nprocs*sizeof(int) );
	if( rec_pe==NULL )
	{
		lis_free2(9,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count);
		LIS_SETERR_MEM(nprocs*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	rec_count = (int *)malloc( nprocs*sizeof(int) );
	if( rec_count==NULL )
	{
		lis_free2(9,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count);
		LIS_SETERR_MEM(nprocs*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	snd_pe = (int *)malloc( nprocs*sizeof(int) );
	if( snd_pe==NULL )
	{
		lis_free2(9,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count);
		LIS_SETERR_MEM(nprocs*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	snd_count = (int *)malloc( nprocs*sizeof(int) );
	if( snd_count==NULL )
	{
		lis_free2(9,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count);
		LIS_SETERR_MEM(nprocs*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}
	snd_is = (int *)malloc( nprocs*sizeof(int) );
	if( snd_is==NULL )
	{
		lis_free2(9,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count);
		LIS_SETERR_MEM(nprocs*sizeof(int));
		return LIS_OUT_OF_MEMORY;
	}

	err = lis_vector_duplicate(Ain,vout);
	if( err )
	{
		return err;
	}
	n     = (*vout)->n;

	rec_n = 0;
	ris   = (*vout)->ranges[my_rank];
	rie   = (*vout)->ranges[my_rank+1];
	/* recv */
	if( rie <= is )
	{
		for(i=my_rank-1;i>=0;i--)
		{
			if( vin->ranges[i] < rie )
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = rie - vin->ranges[i];
				rec_n++;
				break;
			}
		}
		for(i-=1;i>=0;i--)
		{
			if( vin->ranges[i] <= ris )
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = vin->ranges[i+1] - ris;
				rec_n++;
				break;
			}
			else
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = vin->ranges[i+1] - vin->ranges[i];
				rec_n++;
			}
		}
	}
	else if( ris < is )
	{
		for(i=my_rank-1;i>=0;i--)
		{
			if( vin->ranges[i] <= ris )
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = vin->ranges[i+1] - ris;
				rec_n++;
				break;
			}
			else
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = vin->ranges[i+1] - vin->ranges[i];
				rec_n++;
			}
		}
	}
	up_n = rec_n;
	lis_sort_ii(0,rec_n-1,rec_pe,rec_count);
	if( ris >= ie )
	{
		for(i=my_rank+1;i<nprocs;i++)
		{
			if( vin->ranges[i+1] > ris )
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = vin->ranges[i+1] - ris;
				rec_n++;
				break;
			}
		}
		for(i+=1;i<nprocs;i++)
		{
			if( vin->ranges[i+1] >= rie )
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = rie - vin->ranges[i];
				rec_n++;
				break;
			}
			else
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = vin->ranges[i+1] - vin->ranges[i];
				rec_n++;
			}
		}
	}
	else if( rie > ie )
	{
		for(i=my_rank+1;i<nprocs;i++)
		{
			if( vin->ranges[i+1] >= rie )
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = rie - vin->ranges[i];
				rec_n++;
				break;
			}
			else
			{
				rec_pe[rec_n]    = i;
				rec_count[rec_n] = vin->ranges[i+1] - vin->ranges[i];
				rec_n++;
			}
		}
	}

	/* send */
	snd_n = 0;
	nnz   = 0;
	for(i=0;i<nprocs;i++)
	{
		if( i!=my_rank )
		{
			if( (*vout)->ranges[i] < is && (*vout)->ranges[i+1] > is && (*vout)->ranges[i+1] <= ie )
			{
				snd_pe[snd_n]    = i;
				snd_is[snd_n]    = 0;
				snd_count[snd_n] = (*vout)->ranges[i+1] - is;
				snd_n++;
			}
			else if( (*vout)->ranges[i] >=is && (*vout)->ranges[i+1] <= ie )
			{
				snd_pe[snd_n]    = i;
				snd_is[snd_n]    = (*vout)->ranges[i] - is;
				snd_count[snd_n] = (*vout)->ranges[i+1] - (*vout)->ranges[i];
				snd_n++;
			}
			else if( (*vout)->ranges[i] < ie && (*vout)->ranges[i+1] > ie && (*vout)->ranges[i] >= is )
			{
				snd_pe[snd_n]    = i;
				snd_is[snd_n]    = (*vout)->ranges[i] - is;
				snd_count[snd_n] = ie - (*vout)->ranges[i];
				snd_n++;
			}
			else if( (*vout)->ranges[i] < is && (*vout)->ranges[i+1] > ie )
			{
				snd_pe[snd_n]    = i;
				snd_is[snd_n]    = 0;
				snd_count[snd_n] = ie - is;
				snd_n++;
			}
		}
	}
	req1 = (MPI_Request *)malloc( snd_n*sizeof(MPI_Request) );
	if( req1==NULL )
	{
		lis_free2(9,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count);
		LIS_SETERR_MEM(snd_n*sizeof(MPI_Request));
		return LIS_OUT_OF_MEMORY;
	}
	req2 = (MPI_Request *)malloc( rec_n*sizeof(MPI_Request) );
	if( req2==NULL )
	{
		lis_free2(9,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count);
		LIS_SETERR_MEM(rec_n*sizeof(MPI_Request));
		return LIS_OUT_OF_MEMORY;
	}
	sta1 = (MPI_Status *)malloc( snd_n*sizeof(MPI_Status) );
	if( sta1==NULL )
	{
		lis_free2(9,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count);
		LIS_SETERR_MEM(snd_n*sizeof(MPI_Request));
		return LIS_OUT_OF_MEMORY;
	}
	sta2 = (MPI_Status *)malloc( rec_n*sizeof(MPI_Status) );
	if( sta2==NULL )
	{
		lis_free2(9,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count);
		LIS_SETERR_MEM(rec_n*sizeof(MPI_Request));
		return LIS_OUT_OF_MEMORY;
	}

	/* exchange vector */
	i      = 0;
	if( (*vout)->ranges[my_rank] < is && (*vout)->ranges[my_rank+1] > is && (*vout)->ranges[my_rank+1] <= ie )
	{
		k = is - (*vout)->ranges[my_rank];
		for(i=0;i<(*vout)->ranges[my_rank+1]-is;i++)
		{
			(*vout)->value[k+i] = vin->value[i];
		}
	}
	else if( (*vout)->ranges[my_rank] < ie && (*vout)->ranges[my_rank+1] > ie && (*vout)->ranges[my_rank] >= is )
	{
		k = (*vout)->ranges[my_rank] - is;
		for(i=0;i<ie-(*vout)->ranges[my_rank];i++)
		{
			(*vout)->value[i] = vin->value[k+i];
		}
	}
	else if( (*vout)->ranges[my_rank+1] <= ie && (*vout)->ranges[my_rank] >=is )
	{
		k = (*vout)->ranges[my_rank]-is;
		for(i=0;i<(*vout)->ranges[my_rank+1]-(*vout)->ranges[my_rank];i++)
		{
			(*vout)->value[i] = vin->value[k+i];
		}
	}
	else if( (*vout)->ranges[my_rank] < is && (*vout)->ranges[my_rank+1] > ie )
	{
		k = is - (*vout)->ranges[my_rank];
		for(i=0;i<ie-is;i++)
		{
			(*vout)->value[k+i] = vin->value[i];
		}
	}
	count = i;
	for(k=0;k<snd_n;k++)
	{
		MPI_Isend(&vin->value[snd_is[k]],snd_count[k],MPI_DOUBLE,snd_pe[k],0,vin->comm,&req1[k]);
	}
	for(k=0;k<up_n;k++)
	{
		MPI_Irecv((*vout)->value,rec_count[k],MPI_DOUBLE,rec_pe[k],0,vin->comm,&req2[k]);
	}
	for(;k<rec_n;k++)
	{
		MPI_Irecv(&(*vout)->value[count],rec_count[k],MPI_DOUBLE,rec_pe[k],0,vin->comm,&req2[k]);
	}
	MPI_Waitall(rec_n, req2, sta2);
	MPI_Waitall(snd_n, req1, sta1);

#if 0
	is = (*vout)->ranges[my_rank];
	ie = (*vout)->ranges[my_rank+1];
	printf("rank=%d n=%d gn=%d is=%d ie=%d v[0]=%f v[n]=%f\n",my_rank,n,(*vout)->gn,is,ie,(*vout)->value[0],(*vout)->value[ie-is-1]);
#endif

	lis_free2(9,req1,req2,sta1,sta2,snd_pe,snd_count,snd_is,rec_pe,rec_count);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#endif
