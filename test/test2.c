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
#include <math.h>
#include "lis.h"

#undef __FUNC__
#define __FUNC__ "main"
int main(int argc, char* argv[])
{
	LIS_MATRIX		A0,A;
	LIS_VECTOR		x,b,u;
	LIS_SOLVER		solver;
	int				i,k,n,nnz,m,nn,gn,ii,j,jj;
	int				is,ie;
	int				nprocs,my_rank;
	int				np,nsol;
	int				err,iter,mtype,iter_double,iter_quad;
	double			times,itimes,ptimes,p_c_times,p_i_times;
	LIS_REAL		resid;
	char			solvername[128];
	int				*ptr,*index;
	LIS_SCALAR		*value;


	LIS_DEBUG_FUNC_IN;


	lis_initialize(&argc, &argv);

	#ifdef USE_MPI
		MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
		MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	#else
		nprocs  = 1;
		my_rank = 0;
	#endif

	if( argc < 6 )
	{
		if( my_rank==0 ) printf("Usage: test2 m n matrix_type solution_filename residual_filename [options]\n");
		lis_finalize();
		exit(0);
	}

	n  = atoi(argv[1]);
	m  = atoi(argv[2]);
	mtype  = atoi(argv[3]);
	if( n<=0 || m<=0 )
	{
		if( my_rank==0 ) printf("n=%d <=0 or m=%d <=0\n",n,m);
		lis_finalize();
		exit(0);
	}

		if( my_rank==0 )
		  {
		    printf("\n");
		    printf("number of processes = %d\n",nprocs);
		  }

#ifdef _OPENMP
		if( my_rank==0 )
		  {
		    printf("max number of threads = %d\n",omp_get_num_procs());
		    printf("number of threads = %d\n",omp_get_max_threads());
		  }
#endif
		
	/* create matrix and vectors */
	nn = n*m;
	err = lis_matrix_create(LIS_COMM_WORLD,&A);
	err = lis_matrix_set_size(A,0,nn);
	CHKERR(err);

	ptr   = (int *)malloc((A->n+1)*sizeof(int));
	if( ptr==NULL ) CHKERR(1);
	index = (int *)malloc(5*A->n*sizeof(int));
	if( index==NULL ) CHKERR(1);
	value = (LIS_SCALAR *)malloc(5*A->n*sizeof(LIS_SCALAR));
	if( value==NULL ) CHKERR(1);

	lis_matrix_get_range(A,&is,&ie);
	k = 0;
	for(ii=is;ii<ie;ii++)
	{
		i = ii/n;
		j = ii - i*n;
		if( i>0 )   { jj = ii - n; index[k] = jj; value[k++] = -1.0;}
		if( i<m-1 ) { jj = ii + n; index[k] = jj; value[k++] = -1.0;}
		if( j>0 )   { jj = ii - 1; index[k] = jj; value[k++] = -1.0;}
		if( j<n-1 ) { jj = ii + 1; index[k] = jj; value[k++] = -1.0;}
		index[k] = ii; value[k++] = 4.0;
		ptr[ii-is+1] = k;
	}
	ptr[0] = 0;
	err = lis_matrix_set_crs(ptr[ie-is],ptr,index,value,A);
	CHKERR(err);
	err = lis_matrix_assemble(A);
	CHKERR(err);

	nnz = A->nnz;
#ifdef USE_MPI
	MPI_Allreduce(&nnz,&i,1,MPI_INT,MPI_SUM,A->comm);
	nnz   = i;
#endif
	if( my_rank==0 ) printf("matrix size = %d x %d (%d nonzero entries)\n",nn,nn,nnz);

	err = lis_matrix_duplicate(A,&A0);
	CHKERR(err);
	lis_matrix_set_type(A0,mtype);
	err = lis_matrix_convert(A,A0);
	CHKERR(err);
	lis_matrix_destroy(A);
	A = A0;

	err = lis_vector_duplicate(A,&u);
	CHKERR(err);
	err = lis_vector_duplicate(A,&b);
	CHKERR(err);
	err = lis_vector_duplicate(A,&x);
	CHKERR(err);

	err = lis_vector_set_all(1.0,u);
	lis_matvec(A,u,b);

	err = lis_solver_create(&solver); CHKERR(err);
	lis_solver_set_option("-print mem",solver);
	lis_solver_set_optionC(solver);

	err = lis_solve(A,b,x,solver);
	CHKERR(err);
	lis_solver_get_itersex(solver,&iter,&iter_double,&iter_quad);
	lis_solver_get_timeex(solver,&times,&itimes,&ptimes,&p_c_times,&p_i_times);
	lis_solver_get_residualnorm(solver,&resid);
	lis_solver_get_solver(solver,&nsol);
	lis_get_solvername(nsol,solvername);
	if( my_rank==0 )
	{
		printf("%s: number of iterations     = %d (double = %d, quad = %d)\n",solvername,iter, iter_double, iter_quad);
		printf("%s: elapsed time             = %e sec.\n",solvername,times);
		printf("%s:   preconditioner         = %e sec.\n",solvername, ptimes);
		printf("%s:     matrix creation      = %e sec.\n",solvername, p_c_times);
		printf("%s:   linear solver          = %e sec.\n",solvername, itimes);
		printf("%s: relative residual 2-norm = %e\n\n",solvername,resid);
	}

	/* write solution */
	lis_output_vector(x,LIS_FMT_MM,argv[4]); 

	/* write residual */
	lis_solver_output_rhistory(solver, argv[5]); 

	lis_solver_destroy(solver);
	lis_matrix_destroy(A);
	lis_vector_destroy(b);
	lis_vector_destroy(x);
	lis_vector_destroy(u);

	lis_finalize();

	LIS_DEBUG_FUNC_OUT;
	return 0;
}


