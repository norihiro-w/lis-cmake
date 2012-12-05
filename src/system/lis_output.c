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
 * lis_output
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_output"
int lis_output(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, int format, char *path)
{
	int				err;
	LIS_MATRIX		B;

	LIS_DEBUG_FUNC_IN;

	#ifdef USE_MPI
		MPI_Barrier(A->comm);
	#endif
	err = lis_matrix_check(A,LIS_MATRIX_CHECK_ALL);
	if( err ) return err;

	if( format==LIS_FMT_MM || format==LIS_FMT_MMB )
	{
		switch( A->matrix_type )
		{
		case LIS_MATRIX_CRS:
			err = lis_output_mm_crs(A,b,x,format,path);
			break;
		default:
			err = lis_matrix_duplicate(A,&B);
			if( err ) return err;
			lis_matrix_set_type(B,LIS_MATRIX_CRS);
			err = lis_matrix_convert(A,B);
			if( err ) return err;
			err = lis_output_mm_crs(B,b,x,format,path);
			lis_matrix_destroy(B);
			break;
		}
	}

	LIS_DEBUG_FUNC_OUT;
	return err;
}

#undef __FUNC__
#define __FUNC__ "lis_output"
int lis_output_matrix(LIS_MATRIX A, int format, char *path)
{
	int			err;
	LIS_MATRIX		B;
	LIS_VECTOR		b,x;

	LIS_DEBUG_FUNC_IN;

	#ifdef USE_MPI
		MPI_Barrier(A->comm);
	#endif
	err = lis_matrix_check(A,LIS_MATRIX_CHECK_ALL);
	if( err ) return err;

	err = lis_vector_create(LIS_COMM_WORLD,&b);
	err = lis_vector_create(LIS_COMM_WORLD,&x);

	if( format==LIS_FMT_MM || format==LIS_FMT_MMB )
	{
		switch( A->matrix_type )
		{
		case LIS_MATRIX_CRS:
			err = lis_output_mm_crs(A,b,x,format,path);
			break;
		default:
			err = lis_matrix_duplicate(A,&B);
			if( err ) return err;
			lis_matrix_set_type(B,LIS_MATRIX_CRS);
			err = lis_matrix_convert(A,B);
			if( err ) return err;
			err = lis_output_mm_crs(B,b,x,format,path);
			lis_matrix_destroy(B);
			break;
		}
	}

	err = lis_vector_destroy(b);
	err = lis_vector_destroy(x);

	LIS_DEBUG_FUNC_OUT;
	return err;
}

#undef __FUNC__
#define __FUNC__ "lis_output_vector"
int lis_output_vector(LIS_VECTOR v, int format, char *filename)
{
	int				err;

	LIS_DEBUG_FUNC_IN;

	#ifdef USE_MPI
		MPI_Barrier(v->comm);
	#endif
	err = lis_vector_check(v,LIS_VECTOR_CHECK_NULL);
	if( err ) return err;

	switch( format )
	{
	case LIS_FMT_PLAIN:
		err = lis_output_vector_plain(v,filename);
		break;
	case LIS_FMT_MM:
		err = lis_output_vector_mm(v,filename);
		break;
	case LIS_FMT_LIS:
		err = lis_output_vector_lis_ascii(v,filename);
		break;
	default:
		LIS_SETERR(LIS_ERR_ILL_ARG,"ill format option\n");
		return LIS_ERR_ILL_ARG;
		break;
	}

	LIS_DEBUG_FUNC_OUT;
	return err;
}

#undef __FUNC__
#define __FUNC__ "lis_output_vector_plain"
int lis_output_vector_plain(LIS_VECTOR v, char *path)
{
	int				n,i;
	#ifdef USE_MPI
		int				pe,nprocs,my_rank;
		int				err,ret;
	#endif
	FILE			*file;

	LIS_DEBUG_FUNC_IN;

#ifdef USE_MPI
	nprocs    = v->nprocs;
	my_rank   = v->my_rank;
	n         = v->n;

	err = 0;
	for(pe=0;pe<nprocs;pe++)
	{
		if( my_rank==pe )
		{
			if( my_rank==0 )
			{
				file = fopen(path, "w");
			}
			else
			{
				file = fopen(path, "a");
			}
			if( file==NULL )
			{
				LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n", path);
				err = 1;
			}
			for(i=0;i<n;i++)
			{
				fprintf(file, "%28.20e\n", v->value[i]);
			}
			fclose(file);
		}
		MPI_Allreduce(&err,&ret,1,MPI_INT,MPI_SUM,v->comm);
		if( ret )
		{
			return LIS_FAILS;
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	n         = v->n;

	file = fopen(path, "w");
	if( file==NULL )
	{
		LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n", path);
		return LIS_ERR_FILE_IO;
	}
	for(i=0;i<n;i++)
	{
		fprintf(file, "%28.20e\n", v->value[i]);
	}
	fclose(file);
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_output_vector_mm"
int lis_output_vector_mm(LIS_VECTOR v, char *path)
{
	int				i,n,is;
	#ifdef USE_MPI
		int				pe,nprocs,my_rank;
		int				err,ret;
	#endif
	FILE			*file;

	LIS_DEBUG_FUNC_IN;
#ifdef USE_MPI
	nprocs    = v->nprocs;
	my_rank   = v->my_rank;
	n         = v->n;
	is        = v->is; 

	err = 0;
	if( my_rank==0 )
	{
		file = fopen(path, "w");
		if( file==NULL )
		{
			LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n", path);
			err = 1;
		}

		fprintf(file, "%%%%MatrixMarket vector coordinate real general\n");
		fprintf(file, "%d\n", v->gn);
		fclose(file);
	}
	MPI_Allreduce(&err,&ret,1,MPI_INT,MPI_SUM,v->comm);
	if( ret )
	{
		return LIS_FAILS;
	}

	err = 0;
	for(pe=0;pe<nprocs;pe++)
	{
		if( my_rank==pe )
		{
			file = fopen(path, "a");
			if( file==NULL )
			{
				LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n", path);
				err = 1;
			}
			else
			{
				for(i=0;i<n;i++)
				{
					fprintf(file, "%d %28.20e\n", i+is+1, v->value[i]); 
				}
				fclose(file);
			}
		}
		MPI_Allreduce(&err,&ret,1,MPI_INT,MPI_SUM,v->comm);
		if( ret )
		{
			return LIS_FAILS;
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	n         = v->n;
	is        = v->is; 

	file = fopen(path, "w");
	if( file==NULL )
	{
		LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n", path);
		return LIS_ERR_FILE_IO;
	}

	fprintf(file, "%%%%MatrixMarket vector coordinate real general\n");
	fprintf(file, "%d\n", v->gn);

	for(i=0;i<n;i++)
	{
	  fprintf(file, "%d %28.20e\n", i+is+1, v->value[i]);
	  /* fprintf(file, "%d %28.20e\n", i+1, v->value[i]); */
	}
	fclose(file);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}


#undef __FUNC__
#define __FUNC__ "lis_output_vector_lis_ascii"
int lis_output_vector_lis_ascii(LIS_VECTOR v, char *path)
{
	int				i,n;
	#ifdef USE_MPI
		int				pe,nprocs,my_rank;
		int				err,ret;
	#endif
	FILE			*file;

	LIS_DEBUG_FUNC_IN;
#ifdef USE_MPI
	nprocs    = v->nprocs;
	my_rank   = v->my_rank;
	n         = v->n;

	err = 0;
	if( my_rank==0 )
	{
		file = fopen(path, "w");
		if( file==NULL )
		{
			LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n", path);
			err = 1;
		}

		fprintf(file, "#LIS A vec\n");
		fprintf(file, "%d\n", nprocs);
		fclose(file);
	}
	MPI_Allreduce(&err,&ret,1,MPI_INT,MPI_SUM,v->comm);
	if( ret )
	{
		return LIS_FAILS;
	}

	err = 0;
	for(pe=0;pe<nprocs;pe++)
	{
		if( my_rank==pe )
		{
			file = fopen(path, "a");
			if( file==NULL )
			{
				LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n", path);
				err = 1;
			}
			else
			{
				fprintf(file, "# %d %d\n", pe, v->n);
				for(i=0;i<n;i++)
				{
					fprintf(file, "%28.20e ",v->value[i]);
					if( (i+1)%3==0 ) fprintf(file, "\n");
				}
				if( n%3!=0 ) fprintf(file, "\n");
				fclose(file);
			}
		}
		MPI_Allreduce(&err,&ret,1,MPI_INT,MPI_SUM,v->comm);
		if( ret )
		{
			return LIS_FAILS;
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	n         = v->n;

	file = fopen(path, "w");
	if( file==NULL )
	{
		LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n", path);
		return LIS_ERR_FILE_IO;
	}

	fprintf(file, "#LIS A vec\n");
	fprintf(file, "1\n");

	fprintf(file, "# 0 %d\n", v->n);
	for(i=0;i<n;i++)
	{
		fprintf(file, "%28.20e ",v->value[i]);
		if( (i+1)%3==0 ) fprintf(file, "\n");
	}
	if( n%3!=0 ) fprintf(file, "\n");
	fclose(file);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_solver_output_rhistory"
int lis_solver_output_rhistory(LIS_SOLVER solver, char *filename)
{
	int		i,maxiter;
	#ifdef USE_MPI
		int		my_rank,err;
	#endif
	FILE	*file;

	LIS_DEBUG_FUNC_IN;

	maxiter = solver->iter+1;
	if( solver->retcode!=LIS_SUCCESS )
	{
		maxiter--;
	}
#ifdef USE_MPI
	if( solver->residual==NULL )
	{
		LIS_SETERR(LIS_FAILS,"residual history is empty\n");
		return LIS_FAILS;
	}
	if( solver->A==NULL )
	{
		LIS_SETERR(LIS_FAILS,"matrix A is NULL\n");
		return LIS_FAILS;
	}

	MPI_Barrier(solver->A->comm);
	my_rank = solver->A->my_rank;
	err = 0;
	if( my_rank==0 )
	{
		file = fopen(filename, "w");
		if( file==NULL )
		{
			LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n", filename);
			err = 1;
		}
		else
		{
			for(i=0;i<maxiter;i++)
			{
				fprintf(file, "%e\n", solver->residual[i]);
			}
			fclose(file);
		}
	}
	MPI_Barrier(solver->A->comm);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	if( solver->residual==NULL )
	{
		LIS_SETERR(LIS_FAILS,"residual history is empty\n");
		return LIS_FAILS;
	}
	file = fopen(filename, "w");
	if( file==NULL )
	{
		LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n", filename);
		return LIS_ERR_FILE_IO;
	}
	for(i=0;i<maxiter;i++)
	{
		fprintf(file, "%e\n", solver->residual[i]);
	}
	fclose(file);
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_output_rhistory"
int lis_esolver_output_rhistory(LIS_ESOLVER esolver, char *filename)
{
	int		i,maxiter;
	#ifdef USE_MPI
		int		my_rank,err;
	#endif
	FILE	*file;

	LIS_DEBUG_FUNC_IN;

	maxiter = esolver->iter+1;
	if( esolver->retcode!=LIS_SUCCESS )
	{
		maxiter--;
	}
#ifdef USE_MPI
	if( esolver->residual==NULL )
	{
		LIS_SETERR(LIS_FAILS,"eigensolver's residual history is empty\n");
		return LIS_FAILS;
	}
	if( esolver->A==NULL )
	{
		LIS_SETERR(LIS_FAILS,"matrix A is NULL\n");
		return LIS_FAILS;
	}

	MPI_Barrier(esolver->A->comm);
	my_rank = esolver->A->my_rank;
	err = 0;
	if( my_rank==0 )
	{
		file = fopen(filename, "w");
		if( file==NULL )
		{
			LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n", filename);
			err = 1;
		}
		else
		{
			for(i=0;i<maxiter;i++)
			{
				fprintf(file, "%e\n", esolver->residual[i]);
			}
			fclose(file);
		}
	}
	MPI_Barrier(esolver->A->comm);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#else
	if( esolver->residual==NULL )
	{
		LIS_SETERR(LIS_FAILS,"eigensolver's residual history is empty\n");
		return LIS_FAILS;
	}
	file = fopen(filename, "w");
	if( file==NULL )
	{
		LIS_SETERR1(LIS_ERR_FILE_IO,"cannot open file %s\n", filename);
		return LIS_ERR_FILE_IO;
	}
	for(i=0;i<maxiter;i++)
	{
		fprintf(file, "%e\n", esolver->residual[i]);
	}
	fclose(file);
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
#endif
}

