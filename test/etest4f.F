      implicit none
      
#include "lisf.h"
 
      integer           i,n,gn,is,ie,iter,ierr
      LIS_MATRIX        A,B
      LIS_VECTOR        x,y
      LIS_REAL          evalue0
      LIS_ESOLVER       esolver
      integer           my_rank,nprocs
      integer           nsol,len,rhs,iter_double,iter_quad
      real*8            tol
      real*8            times,itimes,ptimes,p_c_times,p_i_times
      LIS_REAL          resid
      character*256     fname,solname,resname,argc
      character*30      nchar,esolvername
      integer*4         iargc
      
      call lis_initialize(ierr) 

#ifdef USE_MPI
      call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierr)
#else
      nprocs  = 1
      my_rank = 0
#endif

      i = iargc()

      if( i.lt.1 ) then
        if( my_rank.eq.0 ) then
          write(*,*) 'etest4f n [options]'
          call lis_finalize(ierr)
        endif
        stop
      endif

      call getarg(1,nchar)
      read (nchar,*) n

      call lis_matrix_create(LIS_COMM_WORLD,A,ierr)
      call lis_matrix_set_size(A,0,n,ierr)
      call lis_matrix_get_size(A,n,gn,ierr)
      call lis_matrix_get_range(A,is,ie,ierr)
      do i=is,ie-1
        if( i>1  ) call lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0d0,
     .                                        A,ierr)
        if( i<gn ) call lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0d0,
     .                                        A,ierr)
        call lis_matrix_set_value(LIS_INS_VALUE,i,i,2.0d0,A,ierr)
      enddo
      call lis_matrix_set_type(A,LIS_MATRIX_CRS,ierr)
      call lis_matrix_assemble(A,ierr)
      
      call lis_vector_duplicate(A,x,ierr)
      call lis_vector_set_all(1.0d0,x,ierr)
      
      call lis_esolver_create(esolver,ierr)
      call lis_esolver_set_option('-eprint mem',esolver,ierr)
      call lis_esolver_set_optionC(esolver,ierr)
      call lis_esolve(A,x,evalue0,esolver,ierr)
      call CHKERR(ierr)
      call lis_esolver_get_itersex(esolver,iter,iter_double,iter_quad,
     .     ierr)
      call lis_esolver_get_timeex(esolver,times,itimes,ptimes,
     .     p_c_times,p_i_times,ierr)
      call lis_esolver_get_residualnorm(esolver,resid,ierr)
      call lis_esolver_get_esolver(esolver,nsol,ierr)
      call lis_get_esolvername(nsol,esolvername,ierr)

      if( my_rank.eq.0 ) then
        write(*,*) esolvername,': extreme eigenvalue = ',evalue0
        write(*,*) esolvername,': iter               = ',iter
        write(*,*) esolvername,': times              = ',times
        write(*,*) esolvername,': p_times            = ',ptimes
        write(*,*) esolvername,':   p_c              = ',p_c_times
        write(*,*) esolvername,':   p_i              = ',p_i_times
        write(*,*) esolvername,': i_times            = ',itimes
        write(*,*) esolvername,': Residual           = ',resid
        write(*,*) 
      endif

!      call lis_vector_create(LIS_COMM_WORLD,y,ierr);
!      call lis_matrix_create(LIS_COMM_WORLD,B,ierr);
!      call lis_esolver_get_evalues(esolver,y,ierr);
!      call lis_esolver_get_evectors(esolver,B,ierr);
!      call lis_output_vector(y,LIS_FMT_MM,'evalues.out',ierr);
!      call lis_output_matrix(B,LIS_FMT_MM,'evectors.out',ierr);
!      call lis_vector_destroy(y,ierr);
!      call lis_matrix_destroy(B,ierr);

      call lis_matrix_destroy(A,ierr)
      call lis_vector_destroy(x,ierr)
      call lis_esolver_destroy(esolver,ierr)
      call lis_finalize(ierr)

      stop
      end
      
