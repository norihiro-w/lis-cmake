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
    int               err,i,n,nnz,is,ie,maxiter,nn,ii,jj;
    int               j,k,m;
    int               nprocs,mtype,my_rank;
    int               nesol;
    LIS_MATRIX        A,A0;
    LIS_VECTOR        x;
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
	if( my_rank==0 ) printf("Usage: etest2 n m matrix_type solution_filename residual_filename [options]\n");
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
		
    /* generate matrix */
    n  = atoi(argv[1]);
    m  = atoi(argv[2]);
    mtype  = atoi(argv[3]);
    if( n<=0 || m<=0 )
      {
	if( my_rank==0 ) printf("n=%d <=0 or m=%d <=0\n",n,m);
	lis_finalize();
	exit(0);
      }
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
    lis_vector_duplicate(A,&x);

    err = lis_esolver_create(&esolver);
    CHKERR(err);

    lis_esolver_set_option("-eprint mem",esolver);
    lis_esolver_set_optionC(esolver);
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

    /* write solution */
    lis_output_vector(x,LIS_FMT_MM,argv[4]);

    /* write residual */
    lis_esolver_output_rhistory(esolver, argv[5]);

    CHKERR(err);

    lis_esolver_destroy(esolver);
    lis_matrix_destroy(A);
    lis_matrix_destroy(A0);
    lis_vector_destroy(x);

    lis_finalize();

    LIS_DEBUG_FUNC_OUT;

    return 0;
}

 
