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
#include <string.h>
#include "lis.h"

#undef __FUNC__
#define __FUNC__ "main"
int main(int argc, char* argv[])
{
    int               err,i,n,nnz,is,ie,maxiter,gn,ii,jj;
    int               j,k,m;
    int               nprocs,mtype,my_rank;
    int               nesol;
    LIS_MATRIX        A,B;
    LIS_VECTOR        x,y;
    LIS_REAL          evalue0;
    LIS_SCALAR        *tmpa;
    LIS_ESOLVER       esolver;
    LIS_REAL          residual;
    int               iters;
    double            times;
    double	      itimes,ptimes,p_c_times,p_i_times;
    int		      *ptr,*index;
    LIS_SCALAR	      *value;
    int               nesolver;
    char	      esolvername[128], solution_filename[128];
    
    LIS_DEBUG_FUNC_IN;
    
    n  = 12;
    lis_initialize(&argc, &argv);

#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
#else
    nprocs  = 1;
    my_rank = 0;
#endif

    lis_matrix_create(LIS_COMM_WORLD,&A);
    lis_matrix_set_size(A,0,n);
    lis_matrix_get_size(A,&n,&gn);
    lis_matrix_get_range(A,&is,&ie);
    for(i=is;i<ie;i++)
    {
        if( i>0   )  lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A);
        if( i<gn-1 ) lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0,A);
        lis_matrix_set_value(LIS_INS_VALUE,i,i,2.0,A);
    }
    lis_matrix_set_type(A,LIS_MATRIX_CRS);
    lis_matrix_assemble(A);

    if( argc < 3 )
      {
	if( my_rank==0 ) printf("Usage: etest5 evalue_filename evector_filename\n");
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
		
    lis_vector_duplicate(A,&x);

    err = lis_esolver_create(&esolver);
    CHKERR(err);

    lis_esolver_set_option("-e si -ss 2 -eprint mem",esolver);
    lis_esolve(A, x, &evalue0, esolver);
    lis_esolver_get_esolver(esolver,&nesol);
    lis_get_esolvername(nesol,esolvername);
    lis_esolver_get_residualnorm(esolver, &residual);
    lis_esolver_get_iters(esolver, &iters);
    lis_esolver_get_timeex(esolver,&times,&itimes,&ptimes,&p_c_times,&p_i_times);
    if( my_rank==0 ) {
      printf("%s: mode number              = %d\n", esolvername, esolver->options[LIS_EOPTIONS_MODE]);
      printf("%s: eigenvalue               = %e\n", esolvername, evalue0);
      printf("%s: number of iterations     = %d\n",esolvername, iters);
      printf("%s: elapsed time             = %e sec.\n", esolvername, times);
      printf("%s:   preconditioner         = %e sec.\n", esolvername, ptimes);
      printf("%s:     matrix creation      = %e sec.\n", esolvername, p_c_times);
      printf("%s:   linear solver          = %e sec.\n", esolvername, itimes);
      printf("%s: relative residual 2-norm = %e\n\n",esolvername, residual);
  }

    lis_vector_create(LIS_COMM_WORLD,&y);
    lis_matrix_create(LIS_COMM_WORLD,&B);
    lis_esolver_get_evalues(esolver,y);
    lis_esolver_get_evectors(esolver,B);

    /* output eigenvalues */
    lis_output_vector(y,LIS_FMT_MM,argv[1]);

    /* output eigenvectors */
    lis_output_matrix(B,LIS_FMT_MM,argv[2]);

    CHKERR(err);

    lis_esolver_destroy(esolver);
    lis_vector_destroy(x);
    lis_matrix_destroy(A);
    lis_vector_destroy(y);
    lis_matrix_destroy(B);

    lis_finalize();

    LIS_DEBUG_FUNC_OUT;

    return 0;
}

 
