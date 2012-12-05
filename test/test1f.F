      implicit none
      
#include "lisf.h"

      LIS_MATRIX    A
      LIS_VECTOR    b,x,u
      LIS_SOLVER    solver
      integer      ierr
      integer      my_rank,nprocs
      integer      omp_get_num_procs,omp_get_max_threads
      integer      n,gn,i
      integer      nsol,iter,len,rhs,iter_double,iter_quad
      real*8       times,itimes,ptimes,p_c_times,p_i_times
      LIS_REAL      resid
      character*256 fname,solname,resname,argc
      character*20 solvername
      integer*4 iargc
      
      call lis_initialize(ierr) 

#ifdef USE_MPI
      call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
#else
      nprocs  = 1
      my_rank = 0
#endif

      i = iargc()
      if( i.lt.4 ) then
         if( my_rank.eq.0 ) then
          write(*,*) 'Usage: test1f matrix_filename rhs_setting ',
     .              'solution_filename residual_filename [options]'
          call lis_finalize(ierr)
        endif
        stop
      endif
      call getarg(2,argc)
      if( argc.eq.'0' ) then
        rhs = 0
      elseif( argc.eq.'1' ) then
        rhs = 1
      elseif( argc.eq.'2' ) then
        rhs = 2
      else
        rhs = -1
      endif
      
      if (my_rank .eq. 0) then
         write(*,*) ''
         write(*,*) 'number of processes = ',nprocs
      endif

#ifdef _OPENMP
      write(*,*) 'max number of threads = ',omp_get_num_procs()
      write(*,*) 'number of threads = ', omp_get_max_threads()
#endif

! read matrix and vectors from file 
      call getarg(1,fname)
      call lis_matrix_create(LIS_COMM_WORLD,A,ierr)
      call CHKERR(ierr)
      call lis_vector_create(LIS_COMM_WORLD,b,ierr)
      call CHKERR(ierr)
      call lis_vector_create(LIS_COMM_WORLD,x,ierr)
      call CHKERR(ierr)
      call lis_matrix_set_type(A,LIS_MATRIX_CRS,ierr)
      call lis_input(A,b,x,fname,ierr)
      call CHKERR(ierr);

      call lis_vector_duplicate(b,u,ierr)

      call lis_matrix_get_size(A,n,gn,ierr)

      call CHKERR(ierr)
      call lis_vector_is_null(b,ierr)
      if( ierr.eq.LIS_TRUE ) then
        call lis_vector_destroy(b,ierr)
        call lis_vector_duplicate(A,b,ierr)
        call CHKERR(ierr)
        if( rhs.eq.0 ) then
          call lis_finalize(ierr)
          stop
        elseif( rhs.eq.1 ) then
          call lis_vector_set_all(1.0d0,b,ierr)
        else
          call lis_vector_set_all(1.0d0,u,ierr)
          call lis_matvec(A,u,b,ierr)
        endif
      endif
      if( rhs.eq.-1 ) then
        call getarg(2,fname)
        call lis_input_vector(b,fname,ierr)
        call CHKERR(ierr)
      endif

      call lis_vector_is_null(x,ierr)
      if( ierr.eq.LIS_TRUE ) then
        call lis_vector_destroy(x,ierr)
        call lis_vector_duplicate(u,x,ierr)
        call CHKERR(ierr)
      endif

      call lis_solver_create(solver,ierr)
      call CHKERR(ierr)
      call lis_solver_set_option('-print mem',solver,ierr)
      call lis_solver_set_optionC(solver,ierr)
      
      call lis_solve(A,b,x,solver,ierr)

      call CHKERR(ierr)

	call lis_solver_get_itersex(solver,iter,iter_double,iter_quad,
     .                           ierr)
	call lis_solver_get_timeex(solver,times,itimes,ptimes,
     .                           p_c_times,p_i_times,ierr)
	call lis_solver_get_residualnorm(solver,resid,ierr)
      call lis_solver_get_solver(solver,nsol,ierr)
      call lis_get_solvername(nsol,solvername,ierr)

      if( my_rank.eq.0 ) then
        write(*,*) solvername,': iter     = ',iter
        write(*,*) solvername,':   double = ',iter_double
        write(*,*) solvername,':   quad   = ',iter_quad
        write(*,*) solvername,': times    = ',times
        write(*,*) solvername,': p_times  = ',ptimes
        write(*,*) solvername,':   p_c    = ',p_c_times
        write(*,*) solvername,':   p_i    = ',p_i_times
        write(*,*) solvername,': i_times  = ',itimes
        write(*,*) solvername,': Residual = ',resid
      endif

! write solution 
      call getarg(3,solname)
	call lis_output_vector(x,LIS_FMT_MM,solname,ierr);

! write residual 
      call getarg(4,resname)
      call lis_solver_output_rhistory(solver, resname, ierr)

      call lis_solver_destroy(solver,ierr)
      call lis_vector_destroy(u,ierr)
      call lis_vector_destroy(x,ierr)
      call lis_vector_destroy(b,ierr)
      call lis_matrix_destroy(A,ierr)

      call lis_finalize(ierr)

      stop
      end
      
