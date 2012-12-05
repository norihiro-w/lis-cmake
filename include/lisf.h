#ifndef __LISF_H__
#define __LISF_H__

#ifndef WIN32
#define SIZEOF_VOID_P 4
#else
#define SIZEOF_VOID_P 4
#endif

#if SIZEOF_VOID_P==8
#define LIS_MATRIX		integer*8
#define LIS_VECTOR		integer*8
#define LIS_SOLVER		integer*8
#define LIS_ESOLVER		integer*8
#else
#define LIS_MATRIX		integer
#define LIS_VECTOR		integer
#define LIS_SOLVER		integer
#define LIS_ESOLVER		integer
#endif

#define LIS_SCALAR		real*8
#define LIS_REAL		real*8
#define LIS_Comm		integer


#ifdef USE_MPI
      include 'mpif.h'
#define LIS_COMM_WORLD	MPI_COMM_WORLD
#else
#define LIS_COMM_WORLD	1
#endif

#define LIS_TRUE		1
#define LIS_FALSE		0

#define LIS_INS_VALUE	0
#define LIS_ADD_VALUE	1

#define LIS_MATRIX_ASSEMBLING	0
#define LIS_MATRIX_CRS			1
#define LIS_MATRIX_CSR			1
#define LIS_MATRIX_CCS			2
#define LIS_MATRIX_CSC			2
#define LIS_MATRIX_MSR			3
#define LIS_MATRIX_DIA			4
#define LIS_MATRIX_CDS			4
#define LIS_MATRIX_ELL			5
#define LIS_MATRIX_JDS			6
#define LIS_MATRIX_JAD			6
#define LIS_MATRIX_BSR			7
#define LIS_MATRIX_BSC			8
#define LIS_MATRIX_VBR			9
#define LIS_MATRIX_COO			10
#define LIS_MATRIX_DENSE		11
#define LIS_MATRIX_DNS			11
#define LIS_MATRIX_RCO			255

#define LIS_MATRIX_TJDS			12
#define LIS_MATRIX_BJDS			13
#define LIS_MATRIX_BCR			14
#define LIS_MATRIX_CJDS			15
#define LIS_MATRIX_PCRS			16
#define LIS_MATRIX_LCRS			17
#define LIS_MATRIX_LJDS			18
#define LIS_MATRIX_LBSR			19
#define LIS_MATRIX_CDIA			20
#define LIS_MATRIX_MSC			21
#define LIS_MATRIX_DECIDING_SIZE 	-(LIS_MATRIX_RCO+1)
#define LIS_MATRIX_NULL				-(LIS_MATRIX_RCO+2)

#define LIS_MATRIX_DEFAULT		LIS_MATRIX_CRS
#define LIS_MATRIX_POINT		LIS_MATRIX_CRS
#define LIS_MATRIX_BLOCK		LIS_MATRIX_BSR

#define LIS_SCALE_NONE				0
#define LIS_SCALE_JACOBI			1
#define LIS_SCALE_SYMM_DIAG			2
	
#define LIS_FMT_AUTO				0
#define LIS_FMT_PLAIN				1
#define LIS_FMT_MM					2
#define LIS_FMT_LIS					3
#define LIS_FMT_LIS_ASCII			3
#define LIS_FMT_LIS_BINARY			4
#define LIS_FMT_FREE				5
#define LIS_FMT_ITBL				6
#define LIS_FMT_HB					7
#define LIS_FMT_MMB					8

#define LIS_BINARY_BIG				0
#define LIS_BINARY_LITTLE			1

#define lis_matrix_create				lis_matrix_create_f
#define LIS_MATRIX_CREATE				lis_matrix_create_f

#define lis_matrix_destroy				lis_matrix_destroy_f
#define LIS_MATRIX_DESTROY				lis_matrix_destroy_f

#define lis_matrix_duplicate			lis_matrix_duplicate_f
#define LIS_MATRIX_DUPLICATE			lis_matrix_duplicate_f

#define lis_matrix_get_size				lis_matrix_get_size_f
#define LIS_MATRIX_GET_SIZE				lis_matrix_get_size_f

#define lis_matrix_get_range			lis_matrix_get_range_f
#define LIS_MATRIX_GET_RANGE			lis_matrix_get_range_f

#define lis_matrix_set_value			lis_matrix_set_value_f
#define LIS_MATRIX_SET_VALUE			lis_matrix_set_value_f

#define lis_matrix_get_type				lis_matrix_get_type_f
#define LIS_MATRIX_GET_TYPE				lis_matrix_get_type_f

#define lis_matrix_set_type				lis_matrix_set_type_f
#define LIS_MATRIX_SET_TYPE				lis_matrix_set_type_f

#define lis_matrix_get_size				lis_matrix_get_size_f
#define LIS_MATRIX_GET_SIZE				lis_matrix_get_size_f

#define lis_matrix_set_size				lis_matrix_set_size_f
#define LIS_MATRIX_SET_SIZE				lis_matrix_set_size_f
	
#define lis_matrix_set_crs				lis_matrix_set_crs_f
#define LIS_MATRIX_SET_CRS				lis_matrix_set_crs_f

#define lis_matrix_set_ccs				lis_matrix_set_ccs_f
#define LIS_MATRIX_SET_CCS				lis_matrix_set_ccs_f

#define lis_matrix_set_msr				lis_matrix_set_msr_f
#define LIS_MATRIX_SET_MSR				lis_matrix_set_msr_f

#define lis_matrix_set_dia				lis_matrix_set_dia_f
#define LIS_MATRIX_SET_DIA				lis_matrix_set_dia_f

#define lis_matrix_set_ell				lis_matrix_set_ell_f
#define LIS_MATRIX_SET_ELL				lis_matrix_set_ell_f

#define lis_matrix_set_jds				lis_matrix_set_jds_f
#define LIS_MATRIX_SET_JDS				lis_matrix_set_jds_f

#define lis_matrix_set_bsr				lis_matrix_set_bsr_f
#define LIS_MATRIX_SET_BSR				lis_matrix_set_bsr_f

#define lis_matrix_set_bsc				lis_matrix_set_bsc_f
#define LIS_MATRIX_SET_BSC				lis_matrix_set_bsc_f

#define lis_matrix_set_coo				lis_matrix_set_coo_f
#define LIS_MATRIX_SET_COO				lis_matrix_set_coo_f

#define lis_matrix_set_dns				lis_matrix_set_dns_f
#define LIS_MATRIX_SET_DNS				lis_matrix_set_dns_f

#define lis_matrix_set_vbr				lis_matrix_set_vbr_f
#define LIS_MATRIX_SET_VBR				lis_matrix_set_vbr_f

#define lis_matrix_assemble				lis_matrix_assemble_f
#define LIS_MATRIX_ASSEMBLE				lis_matrix_assemble_f

#define lis_matrix_is_assembled				lis_matrix_is_assembled_f
#define LIS_MATRIX_IS_ASSEMBLED				lis_matrix_is_assembled_f

#define lis_matrix_malloc				lis_matrix_malloc_f
#define LIS_MATRIX_MALLOC				lis_matrix_malloc_f

#define lis_matrix_malloc_crs				lis_matrix_malloc_crs_f
#define LIS_MATRIX_MALLOC_CRS				lis_matrix_malloc_crs_f

#define lis_matrix_malloc_ccs				lis_matrix_malloc_ccs_f
#define LIS_MATRIX_MALLOC_CCS				lis_matrix_malloc_ccs_f

#define lis_matrix_malloc_bsr				lis_matrix_malloc_bsr_f
#define LIS_MATRIX_MALLOC_BSR				lis_matrix_malloc_bsr_f

#define lis_matrix_malloc_msr				lis_matrix_malloc_msr_f
#define LIS_MATRIX_MALLOC_MSR				lis_matrix_malloc_msr_f

#define lis_matrix_malloc_ell				lis_matrix_malloc_ell_f
#define LIS_MATRIX_MALLOC_ELL				lis_matrix_malloc_ell_f

#define lis_matrix_malloc_jds				lis_matrix_malloc_jds_f
#define LIS_MATRIX_MALLOC_JDS				lis_matrix_malloc_jds_f

#define lis_matrix_malloc_dia				lis_matrix_malloc_dia_f
#define LIS_MATRIX_MALLOC_DIA				lis_matrix_malloc_dia_f

#define lis_matrix_malloc_bsc				lis_matrix_malloc_bsc_f
#define LIS_MATRIX_MALLOC_BSC				lis_matrix_malloc_bsc_f

#define lis_matrix_malloc_vbr				lis_matrix_malloc_vbr_f
#define LIS_MATRIX_MALLOC_VBR				lis_matrix_malloc_vbr_f

#define lis_matrix_malloc_coo				lis_matrix_malloc_coo_f
#define LIS_MATRIX_MALLOC_COO				lis_matrix_malloc_coo_f

#define lis_matrix_malloc_dns				lis_matrix_malloc_dns_f
#define LIS_MATRIX_MALLOC_DNS				lis_matrix_malloc_dns_f

#define lis_matrix_convert				lis_matrix_convert_f
#define LIS_MATRIX_CONVERT				lis_matrix_convert_f

#define lis_matrix_copy					lis_matrix_copy_f
#define LIS_MATRIX_COPY					lis_matrix_copy_f

#define lis_matrix_scaling				lis_matrix_scaling_f
#define LIS_MATRIX_SCALING				lis_matrix_scaling_f

#define lis_matrix_get_diagonal			lis_matrix_get_diagonal_f
#define LIS_MATRIX_GET_DIAGONAL			lis_matrix_get_diagonal_f

#define lis_matrix_set_blocksize		lis_matrix_set_blocksize_f
#define LIS_MATRIX_SET_BLOCKSIZE		lis_matrix_set_blocksize_f

#define lis_vector_create				lis_vector_create_f
#define LIS_VECTOR_CREATE				lis_vector_create_f

#define lis_vector_destroy				lis_vector_destroy_f
#define LIS_VECTOR_DESTROY				lis_vector_destroy_f

#define lis_vector_duplicate			lis_vector_duplicate_f
#define LIS_VECTOR_DUPLICATE			lis_vector_duplicate_f

#define lis_vector_get_size				lis_vector_get_size_f
#define LIS_VECTOR_GET_SIZE				lis_vector_get_size_f

#define lis_vector_set_size				lis_vector_set_size_f
#define LIS_VECTOR_SET_SIZE				lis_vector_set_size_f

#define lis_vector_get_range			lis_vector_get_range_f
#define LIS_VECTOR_GET_RANGE			lis_vector_get_range_f

#define lis_vector_set_value			lis_vector_set_value_f
#define LIS_VECTOR_SET_VALUE			lis_vector_set_value_f

#define lis_vector_set_values			lis_vector_set_values_f
#define LIS_VECTOR_SET_VALUES			lis_vector_set_values_f

#define lis_vector_set_values2			lis_vector_set_values2_f
#define LIS_VECTOR_SET_VALUES2			lis_vector_set_values2_f

#define lis_vector_get_value			lis_vector_get_value_f
#define LIS_VECTOR_GET_VALUE			lis_vector_get_value_f

#define lis_vector_get_values			lis_vector_get_values_f
#define LIS_VECTOR_GET_VALUES			lis_vector_get_values_f

#define lis_vector_scatter			lis_vector_scatter_f
#define LIS_VECTOR_SCATTER			lis_vector_scatter_f

#define lis_vector_gather			lis_vector_gather_f
#define LIS_VECTOR_GATHER			lis_vector_gather_f

#define lis_vector_print				lis_vector_print_f
#define LIS_VECTOR_PRINT				lis_vector_print_f

#define lis_vector_set_all				lis_vector_set_all_f
#define LIS_VECTOR_SET_ALL				lis_vector_set_all_f

#define lis_vector_axpy					lis_vector_axpy_f
#define LIS_VECTOR_AXPY					lis_vector_axpy_f

#define lis_vector_xpay					lis_vector_xpay_f
#define LIS_VECTOR_XPAY					lis_vector_xpay_f

#define lis_vector_axpyz				lis_vector_axpyz_f
#define LIS_VECTOR_AXPYZ				lis_vector_axpyz_f

#define lis_vector_copy					lis_vector_copy_f
#define LIS_VECTOR_COPY					lis_vector_copy_f

#define lis_vector_scale				lis_vector_scale_f
#define LIS_VECTOR_SCALE				lis_vector_scale_f

#define lis_vector_pmul					lis_vector_pmul_f
#define LIS_VECTOR_PMUL					lis_vector_pmul_f

#define lis_vector_pdiv					lis_vector_pdiv_f
#define LIS_VECTOR_PDIV					lis_vector_pdiv_f

#define lis_vector_abs					lis_vector_abs_f
#define LIS_VECTOR_ABS					lis_vector_abs_f

#define lis_vector_reciprocal			lis_vector_reciprocal_f
#define LIS_VECTOR_RECIPROCAL			lis_vector_reciprocal_f

#define lis_vector_shift				lis_vector_shift_f
#define LIS_VECTOR_SHIFT				lis_vector_shift_f

#define lis_vector_dot					lis_vector_dot_f
#define LIS_VECTOR_DOT					lis_vector_dot_f

#define lis_vector_nrm1					lis_vector_nrm1_f
#define LIS_VECTOR_NRM1					lis_vector_nrm1_f

#define lis_vector_nrm2					lis_vector_nrm2_f
#define LIS_VECTOR_NRM2					lis_vector_nrm2_f

#define lis_vector_nrmi					lis_vector_nrmi_f
#define LIS_VECTOR_NRMI					lis_vector_nrmi_f

#define lis_vector_sum					lis_vector_sum_f
#define LIS_VECTOR_SUM					lis_vector_sum_f

#define lis_vector_is_null				lis_vector_is_null_f
#define LIS_VECTOR_IS_NULL				lis_vector_is_null_f

#define lis_solver_create				lis_solver_create_f
#define LIS_SOLVER_CREATE				lis_solver_create_f

#define lis_solver_destroy				lis_solver_destroy_f
#define LIS_SOLVER_DESTROY				lis_solver_destroy_f

#define lis_solve						lis_solve_f
#define LIS_SOLVE						lis_solve_f

#define lis_solve_kernel					lis_solve_kernel_f
#define LIS_SOLVE_KERNEL					lis_solve_kernel_f

#define lis_solver_set_option			lis_solver_set_option_f
#define LIS_SOLVER_SET_OPTION			lis_solver_set_option_f

#define lis_solver_set_optionC			lis_solver_set_optionc_f
#define LIS_SOLVER_SET_OPTIONC			lis_solver_set_optionc_f

#define lis_solver_get_iters			lis_solver_get_iters_f
#define LIS_SOLVER_GET_ITERS			lis_solver_get_iters_f

#define lis_solver_get_itersex			lis_solver_get_itersex_f
#define LIS_SOLVER_GET_ITERSEX			lis_solver_get_itersex_f

#define lis_solver_get_time				lis_solver_get_time_f
#define LIS_SOLVER_GET_TIME				lis_solver_get_time_f

#define lis_solver_get_timeex			lis_solver_get_timeex_f
#define LIS_SOLVER_GET_TIMEEX			lis_solver_get_timeex_f

#define lis_solver_get_residualnorm		lis_solver_get_residualnorm_f
#define LIS_SOLVER_GET_RESIDUALNORM		lis_solver_get_residualnorm_f

#define lis_solver_get_solver			lis_solver_get_solver_f
#define LIS_SOLVER_GET_SOLVER			lis_solver_get_solver_f

#define lis_get_solvername				lis_get_solvername_f
#define LIS_GET_SOLVERNAME				lis_get_solvername_f

#define lis_esolver_create				lis_esolver_create_f
#define LIS_ESOLVER_CREATE				lis_esolver_create_f

#define lis_esolver_destroy				lis_esolver_destroy_f
#define LIS_ESOLVER_DESTROY				lis_esolver_destroy_f

#define lis_esolve						lis_esolve_f
#define LIS_ESOLVE						lis_esolve_f

#define lis_esolver_set_option			lis_esolver_set_option_f
#define LIS_ESOLVER_SET_OPTION			lis_esolver_set_option_f

#define lis_esolver_set_optionC			lis_esolver_set_optionc_f
#define LIS_ESOLVER_SET_OPTIONC			lis_esolver_set_optionc_f

#define lis_esolver_get_iters			lis_esolver_get_iters_f
#define LIS_ESOLVER_GET_ITERS			lis_esolver_get_iters_f

#define lis_esolver_get_itersex			lis_esolver_get_itersex_f
#define LIS_ESOLVER_GET_ITERSEX			lis_esolver_get_itersex_f

#define lis_esolver_get_time				lis_esolver_get_time_f
#define LIS_ESOLVER_GET_TIME				lis_esolver_get_time_f

#define lis_esolver_get_timeex			lis_esolver_get_timeex_f
#define LIS_ESOLVER_GET_TIMEEX			lis_esolver_get_timeex_f

#define lis_esolver_get_evalues			lis_esolver_get_evalues_f
#define LIS_ESOLVER_GET_EVALUES			lis_esolver_get_evalues_f

#define lis_esolver_get_evectors		lis_esolver_get_evectors_f
#define LIS_ESOLVER_GET_EVECTORS		lis_esolver_get_evectors_f

#define lis_esolver_get_residualnorm		lis_esolver_get_residualnorm_f
#define LIS_ESOLVER_GET_RESIDUALNORM		lis_esolver_get_residualnorm_f

#define lis_esolver_get_esolver			lis_esolver_get_esolver_f
#define LIS_ESOLVER_GET_ESOLVER			lis_esolver_get_esolver_f

#define lis_get_esolvername				lis_get_esolvername_f
#define LIS_GET_ESOLVERNAME				lis_get_esolvername_f

#define lis_matvec						lis_matvec_f
#define LIS_MATVEC						lis_matvec_f

#define lis_matvect						lis_matvect_f
#define LIS_MATVECT						lis_matvect_f

#define lis_initialize					lis_init_f
#define LIS_INITIALIZE					lis_init_f

#define lis_initializef					lis_initializef_f
#define LIS_INITIALIZEF					lis_initializef_f

#define lis_finalize					lis_finalize_f
#define LIS_FINALIZE					lis_finalize_f

#define lis_wtime						lis_wtime_f
#define LIS_WTIME						lis_wtime_f

#define chkerr							chkerr_f
#define CHKERR							chkerr_f

#define lis_set_argv_begin				lis_set_argv_begin_f
#define LIS_SET_ARGV_BEGIN				lis_set_argv_begin_f

#define lis_set_argv					lis_set_argv_f
#define LIS_SET_ARGV					lis_set_argv_f

#define lis_set_argv_end				lis_set_argv_end_f
#define LIS_SET_ARGV_END				lis_set_argv_end_f

#define lis_arg2args					lis_arg2args_f
#define LIS_ARG2ARGS					lis_arg2args_f

#define lis_input						lis_input_f
#define LIS_INPUT						lis_input_f

#define lis_output						lis_output_f
#define LIS_OUTPUT						lis_output_f

#define lis_output_matrix				lis_output_matrix_f
#define LIS_OUTPUT_MATRIX				lis_output_matrix_f

#define lis_output_vector				lis_output_vector_f
#define LIS_OUTPUT_VECTOR				lis_output_vector_f

#define lis_input_matrix				lis_input_matrix_f
#define LIS_INPUT_MATRIX				lis_input_matrix_f

#define lis_input_vector				lis_input_vector_f
#define LIS_INPUT_VECTOR				lis_input_vector_f

#define lis_solver_output_rhistory		lis_solver_output_rhistory_f
#define LIS_SOLVER_OUTPUT_RHISTORY		lis_solver_output_rhistory_f

#define lis_esolver_output_rhistory		lis_esolver_output_rhistory_f
#define LIS_ESOLVER_OUTPUT_RHISTORY		lis_esolver_output_rhistory_f

#endif
