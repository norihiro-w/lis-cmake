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
#include <string.h>
#include <math.h>
#include "lis.h"


#undef __FUNC__
#define __FUNC__ "main"
int main(int argc, char* argv[])
{
	LIS_MATRIX		A0,A,B;
	LIS_VECTOR		x,b,u;
	LIS_SOLVER		solver;
	int				nprocs,my_rank;
	int				nsol,rhs,len,i,j,k,jj,kk,p,nrow_p,l,n;
	int				err,iter,iter_double,iter_quad;
	int				*iw,*nrow,*index,*ptr;
	LIS_SCALAR		*s,t,*value;
	double			times,itimes,ptimes,p_c_times,p_i_times;
	LIS_REAL		resid;
	char			solvername[128];

	LIS_DEBUG_FUNC_IN;

	lis_initialize(&argc, &argv);

	#ifdef USE_MPI
		MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
		MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	#else
		nprocs  = 1;
		my_rank = 0;
	#endif

	if( argc < 5 )
	{
		if( my_rank==0 ) printf("Usage: test1 matrix_filename rhs_setting solution_filename residual_filename [options]\n");
		lis_finalize();
		exit(0);
	}

	len = (int)strlen(argv[2]);
	if( len==1 )
	{
		if( argv[2][0]=='0' || argv[2][0]=='1' || argv[2][0]=='2' )
		{
			rhs = atoi(argv[2]);
		}
		else
		{
			rhs = -1;
		}
	}
	else
	{
		rhs = -1;
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

	/* read matrix and vectors from file */
	err = lis_matrix_create(LIS_COMM_WORLD,&A); CHKERR(err);
	err = lis_vector_create(LIS_COMM_WORLD,&b); CHKERR(err);
	err = lis_vector_create(LIS_COMM_WORLD,&x); CHKERR(err);
	err = lis_input(A,b,x,argv[1]);
	CHKERR(err);

	err = lis_matrix_duplicate(A,&A0);
	CHKERR(err);
	lis_matrix_set_type(A0,LIS_MATRIX_CRS);
	err = lis_matrix_convert(A,A0);
	CHKERR(err);
	lis_matrix_destroy(A);
	A = A0;

	err = lis_vector_duplicate(A,&u);
	CHKERR(err);
	if( lis_vector_is_null(b) )
	{
		lis_vector_destroy(b);
		lis_vector_duplicate(A,&b);
		CHKERR(err);
		if( rhs==0 )
		{
			lis_finalize();
			exit(0);
		}
		else if( rhs==1 )
		{
			err = lis_vector_set_all(1.0,b);
		}
		else
		{
			err = lis_vector_set_all(1.0,u);
			lis_matvec(A,u,b);
		}
	}
	if( rhs==-1 )
	{
		lis_input_vector(b,argv[2]);
	}
	if( lis_vector_is_null(x) )
	{
		lis_vector_destroy(x);
		err = lis_vector_duplicate(A,&x);
		CHKERR(err);
	}

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
	

	/* write results */
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
	lis_output_vector(x,LIS_FMT_MM,argv[3]);

	/* write residual */
	lis_solver_output_rhistory(solver, argv[4]);

	lis_solver_destroy(solver);
	lis_vector_destroy(x);
	lis_vector_destroy(u);
	lis_vector_destroy(b);
	lis_matrix_destroy(A);

	lis_finalize();

	LIS_DEBUG_FUNC_OUT;

	return 0;
}


