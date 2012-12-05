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
#include <ctype.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_solver_init
 * lis_solver_create
 * lis_solver_destroy
 * lis_solver_work_destroy
 * lis_solver_set_option
 * lis_solver_get_option
 * lis_solve
 ************************************************/

#define LIS_SOLVERS_LEN			23
#define LIS_PRECON_TYPE_LEN		12


LIS_SOLVER_CHECK_PARAMS lis_solver_check_params[] = {
	NULL,
	lis_cg_check_params       , lis_bicg_check_params      , lis_cgs_check_params, 
	lis_bicgstab_check_params , lis_bicgstabl_check_params , lis_gpbicg_check_params, 
	lis_tfqmr_check_params    , lis_orthomin_check_params  , lis_gmres_check_params,
	lis_jacobi_check_params   , lis_gs_check_params        , lis_sor_check_params,
	lis_bicgsafe_check_params , lis_cr_check_params        , lis_bicr_check_params,
	lis_crs_check_params      , lis_bicrstab_check_params  , lis_gpbicr_check_params,
	lis_bicrsafe_check_params , lis_fgmres_check_params    , lis_idrs_check_params,
	lis_minres_check_params   , lis_idr1_check_params
};

LIS_SOLVER_MALLOC_WORK lis_solver_malloc_work[] = {
	NULL,
	lis_cg_malloc_work       , lis_bicg_malloc_work      , lis_cgs_malloc_work, 
	lis_bicgstab_malloc_work , lis_bicgstabl_malloc_work , lis_gpbicg_malloc_work, 
	lis_tfqmr_malloc_work    , lis_orthomin_malloc_work  , lis_gmres_malloc_work,
	lis_jacobi_malloc_work   , lis_gs_malloc_work        , lis_sor_malloc_work,
	lis_bicgsafe_malloc_work , lis_cr_malloc_work        , lis_bicr_malloc_work,
	lis_crs_malloc_work      , lis_bicrstab_malloc_work  , lis_gpbicr_malloc_work,
	lis_bicrsafe_malloc_work , lis_fgmres_malloc_work    , lis_idrs_malloc_work,
	lis_minres_malloc_work   , lis_idr1_malloc_work
};

LIS_SOLVER_EXECUTE lis_solver_execute[] = {
	NULL,
	lis_cg       , lis_bicg      , lis_cgs, 
	lis_bicgstab , lis_bicgstabl , lis_gpbicg, 
	lis_tfqmr    , lis_orthomin  , lis_gmres,
	lis_jacobi   , lis_gs        , lis_sor,
	lis_bicgsafe , lis_cr        , lis_bicr,
	lis_crs      , lis_bicrstab  , lis_gpbicr,
	lis_bicrsafe , lis_fgmres    , lis_idrs, 
	lis_minres   , lis_idr1
};

#ifdef USE_QUAD_PRECISION
	LIS_SOLVER_EXECUTE lis_solver_execute_quad[] = {
		NULL,
		lis_cg_quad       , lis_bicg_quad      , lis_cgs_quad, 
		lis_bicgstab_quad , lis_bicgstabl_quad , lis_gpbicg_quad,
		lis_tfqmr_quad    , lis_orthomin_quad  , lis_gmres_quad,
		NULL              , NULL               , NULL,
		lis_bicgsafe_quad , lis_cr_quad        , lis_bicr_quad,
		lis_crs_quad      , lis_bicrstab_quad  , lis_gpbicr_quad,
		lis_bicrsafe_quad , lis_fgmres_quad    , NULL
	};
	LIS_SOLVER_EXECUTE lis_solver_execute_switch[] = {
		NULL,
		lis_cg_switch       , lis_bicg_switch , lis_cgs_switch, 
		lis_bicgstab_switch , NULL            , lis_gpbicg_switch,
		NULL                , NULL            , lis_gmres_switch,
		NULL                , NULL            , NULL,
		NULL                , NULL            , NULL,
		NULL                , NULL            , NULL,
		NULL                , NULL            , NULL
	};
	/*
	LIS_SOLVER_EXECUTE lis_solver_execute_periodic[] = {
		NULL,
		lis_cg_periodic       , lis_bicg_periodic , lis_cgs_periodic, 
		lis_bicgstab_periodic , NULL              , lis_gpbicg_periodic,
		NULL                  , NULL              , NULL,
		NULL                  , NULL              , NULL,
		NULL
	};
	*/
#endif

LIS_SOLVER_GET_RESIDUAL lis_solver_get_residual[] = {
	lis_solver_get_residual_nrm2_r, 
	lis_solver_get_residual_nrm2_r, 
	lis_solver_get_residual_nrm1_b
};

int LIS_USE_AT_TYPE[] = {
	0,
	LIS_MATRIX_CCS,LIS_MATRIX_CRS
	};
#define LIS_SOLVER_OPTION_LEN	46
#define LIS_PRINT_LEN			4
#define LIS_SCALE_LEN			3
#define LIS_TRUEFALSE_LEN		2
#define LIS_PRECISION_LEN		3
#define LIS_STORAGE_LEN			11
#define LIS_CONV_COND_LEN		3

char *LIS_SOLVER_OPTNAME[] = {
	"-maxiter",           "-tol",           "-print",          "-scale",         "-ssor_w",
	"-ilu_fill",          "-ilu_relax",     "-is_alpha",       "-is_level",      "-is_m",
	"-hybrid_maxiter",    "-hybrid_ell",    "-hybrid_restart", "-hybrid_tol",    "-hybrid_w",
	"-hybrid_i",          "-sainv_drop",    "-ric2s_tau",      "-ric2s_sigma",   "-ric2s_gamma",
	"-restart",           "-ell",           "-omega",          "-i",             "-p",
	"-f",                 "-h",             "-ver",            "-hybrid_p",      "-initx_zeros",
	"-adds",              "-adds_iter",     "-f",              "-use_at",        "-switch_tol",
	"-switch_maxiter",    "-saamg_unsym",   "-iluc_drop",      "-iluc_gamma",    "-iluc_rate",
	"-storage",           "-storage_block", "-conv_cond",      "-tol_w",         "-saamg_theta",	"-irestart"
};

int LIS_SOLVER_OPTACT[] = {
	LIS_OPTIONS_MAXITER          , LIS_PARAMS_RESID          , LIS_OPTIONS_OUTPUT        , LIS_OPTIONS_SCALE        , LIS_PARAMS_SSOR_W ,
	LIS_OPTIONS_FILL             , LIS_PARAMS_RELAX          , LIS_PARAMS_ALPHA          , LIS_OPTIONS_ISLEVEL      , LIS_OPTIONS_M     ,
	LIS_OPTIONS_PMAXITER         , LIS_OPTIONS_PELL          , LIS_OPTIONS_PRESTART      , LIS_PARAMS_PRESID        , LIS_PARAMS_PW     ,
	LIS_OPTIONS_PSOLVER          , LIS_PARAMS_DROP           , LIS_PARAMS_TAU            , LIS_PARAMS_SIGMA         , LIS_PARAMS_GAMMA  ,
	LIS_OPTIONS_RESTART          , LIS_OPTIONS_ELL           , LIS_PARAMS_W              , LIS_OPTIONS_SOLVER       , LIS_OPTIONS_PRECON,
	LIS_OPTIONS_FILE             , LIS_OPTIONS_HELP          , LIS_OPTIONS_VER           , LIS_OPTIONS_PPRECON      , LIS_OPTIONS_INITGUESS_ZEROS,
	LIS_OPTIONS_ADDS             , LIS_OPTIONS_ADDS_ITER     , LIS_OPTIONS_PRECISION     , LIS_OPTIONS_USE_AT       , LIS_PARAMS_SWITCH_RESID,
	LIS_OPTIONS_SWITCH_MAXITER   , LIS_OPTIONS_SAAMG_UNSYM   , LIS_PARAMS_DROP           , LIS_PARAMS_GAMMA         , LIS_PARAMS_RATE, 
	LIS_OPTIONS_STORAGE          , LIS_OPTIONS_STORAGE_BLOCK , LIS_OPTIONS_CONV_COND     , LIS_PARAMS_RESID_WEIGHT  , LIS_PARAMS_SAAMG_THETA, LIS_OPTIONS_IDRS_RESTART
};

char *lis_solver_atoi[]    = {"cg", "bicg", "cgs", "bicgstab", "bicgstabl", "gpbicg", "tfqmr","orthomin", "gmres", "jacobi", "gs", "sor", "bicgsafe", "cr", "bicr", "crs", "bicrstab", "gpbicr", "bicrsafe", "fgmres", "idrs", "minres", "idr1"};
char *lis_precon_atoi[]    = {"none", "jacobi", "ilu", "ssor", "hybrid", "is", "sainv", "saamg", "iluc", "ilut", "bjacobi", ""};
char *lis_storage_atoi[]   = {"crs", "ccs", "msr", "dia", "ell", "jds", "bsr", "bsc", "vbr", "coo", "dns"};
char *lis_print_atoi[]     = {"none", "mem", "out", "all"};
char *lis_scale_atoi[]     = {"none", "jacobi", "symm_diag"};
char *lis_truefalse_atoi[] = {"false", "true"};
char *lis_precision_atoi[] = {"double", "quad", "switch"};
char *lis_conv_cond_atoi[] = {"nrm2_r", "nrm2_b", "nrm1_b"};

char *lis_solvername[] = {"", "CG", "BiCG", "CGS", "BiCGSTAB", "BiCGSTAB(l)", "GPBiCG", "TFQMR", "Orthomin", "GMRES", "Jacobi",	"Gauss-Seidel", "SOR", "BiCGSafe", "CR", "BiCR", "CRS", "BiCRSTAB", "GPBiCR", "BiCRSafe", "FGMRES", "IDR(s)", "MINRES", "IDR(1)"};
char *lis_preconname[] = {"none", "Jacobi", "ILU", "SSOR", "Hybrid", "I+S", "SAINV", "SAAMG", "Crout ILU", "ILUT", "Block Jacobi"};
char *lis_returncode[] = {"LIS_SUCCESS", "LIS_ILL_OPTION", "LIS_BREAKDOWN", "LIS_OUT_OF_MEMORY", "LIS_MAXITER", "LIS_NOT_IMPLEMENTED", "LIS_ERR_FILE_IO"};
char *lis_precisionname[] = {"double", "quad", "switch"};
char *lis_storagename[]   = {"CRS", "CCS", "MSR", "DIA", "ELL", "JDS", "BSR", "BSC", "VBR", "COO", "DNS"};

LIS_VECTOR lis_solver_residual_history = NULL;

#undef __FUNC__
#define __FUNC__ "lis_solver_init"
int lis_solver_init(LIS_SOLVER solver)
{

	LIS_DEBUG_FUNC_IN;

	solver->A        = NULL;
	solver->At       = NULL;
	solver->b        = NULL;
	solver->x        = NULL;
	solver->d        = NULL;
	solver->work     = NULL;
	solver->residual = NULL;
	solver->precon   = NULL;

	solver->worklen   = 0;
	solver->iter      = 0;
	solver->iter2     = 0;
	solver->resid     = 0;
	solver->times     = 0;
	solver->itimes    = 0;
	solver->ptimes    = 0;
	solver->precision = LIS_PRECISION_DOUBLE;

	solver->options[LIS_OPTIONS_SOLVER]               = LIS_SOLVER_BICG;
	solver->options[LIS_OPTIONS_PRECON]               = LIS_PRECON_TYPE_NONE;
	solver->options[LIS_OPTIONS_OUTPUT]               = LIS_FALSE;
	solver->options[LIS_OPTIONS_MAXITER]              = 1000;
	solver->options[LIS_OPTIONS_RESTART]              = 40;
	solver->options[LIS_OPTIONS_ELL]                  = 2;
	solver->options[LIS_OPTIONS_SCALE]                = LIS_SCALE_NONE;
	solver->options[LIS_OPTIONS_FILL]                 = 0;
	solver->options[LIS_OPTIONS_M]                    = 3;
	solver->options[LIS_OPTIONS_PSOLVER]              = LIS_SOLVER_SOR;
	solver->options[LIS_OPTIONS_PMAXITER]             = 25;
	solver->options[LIS_OPTIONS_PRESTART]             = 40;
	solver->options[LIS_OPTIONS_PELL]                 = 2;
	solver->options[LIS_OPTIONS_PPRECON]              = LIS_PRECON_TYPE_NONE;
	solver->options[LIS_OPTIONS_ISLEVEL]              = 1;
	solver->options[LIS_OPTIONS_INITGUESS_ZEROS]      = LIS_TRUE;
	solver->options[LIS_OPTIONS_ADDS]                 = LIS_FALSE;
	solver->options[LIS_OPTIONS_ADDS_ITER]            = 1;
	solver->options[LIS_OPTIONS_PRECISION]            = LIS_PRECISION_DOUBLE;
	solver->options[LIS_OPTIONS_USE_AT]               = LIS_FALSE;
	solver->options[LIS_OPTIONS_SWITCH_MAXITER]       = -1;
	solver->options[LIS_OPTIONS_SAAMG_UNSYM]          = LIS_FALSE;
	solver->options[LIS_OPTIONS_STORAGE]              = 0;
	solver->options[LIS_OPTIONS_STORAGE_BLOCK]        = 2;
	solver->options[LIS_OPTIONS_CONV_COND]            = 0;
	solver->options[LIS_OPTIONS_INIT_SHADOW_RESID]    = LIS_RESID;
	solver->options[LIS_OPTIONS_IDRS_RESTART]         = 2;

	solver->params[LIS_PARAMS_RESID        -LIS_OPTIONS_LEN] = 1.0e-12;
	solver->params[LIS_PARAMS_RESID_WEIGHT -LIS_OPTIONS_LEN] = 1.0;
	solver->params[LIS_PARAMS_W            -LIS_OPTIONS_LEN] = 1.9;
	solver->params[LIS_PARAMS_SSOR_W       -LIS_OPTIONS_LEN] = 1.0;
	solver->params[LIS_PARAMS_RELAX        -LIS_OPTIONS_LEN] = 1.0;
	solver->params[LIS_PARAMS_DROP         -LIS_OPTIONS_LEN] = 0.05;
	solver->params[LIS_PARAMS_ALPHA        -LIS_OPTIONS_LEN] = 1.0;
	solver->params[LIS_PARAMS_TAU          -LIS_OPTIONS_LEN] = 0.05;
	solver->params[LIS_PARAMS_SIGMA        -LIS_OPTIONS_LEN] = 2.0;
	solver->params[LIS_PARAMS_GAMMA        -LIS_OPTIONS_LEN] = 1.0;
	solver->params[LIS_PARAMS_PRESID       -LIS_OPTIONS_LEN] = 1.0e-3;
	solver->params[LIS_PARAMS_PW           -LIS_OPTIONS_LEN] = 1.5;
	solver->params[LIS_PARAMS_SWITCH_RESID -LIS_OPTIONS_LEN] = 1.0e-12;
	solver->params[LIS_PARAMS_RATE         -LIS_OPTIONS_LEN] = 5.0;
	solver->params[LIS_PARAMS_SAAMG_THETA  -LIS_OPTIONS_LEN] = 0.05;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_create"
int lis_solver_create(LIS_SOLVER *solver)
{
	LIS_DEBUG_FUNC_IN;

	*solver = NULL;

	*solver = (LIS_SOLVER)lis_malloc( sizeof(struct LIS_SOLVER_STRUCT),"lis_solver_create::solver" );
	if( NULL==*solver )
	{
		LIS_SETERR_MEM(sizeof(struct LIS_SOLVER_STRUCT));
		return LIS_OUT_OF_MEMORY;
	}
	lis_solver_init(*solver);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_work_destroy"
int lis_solver_work_destroy(LIS_SOLVER solver)
{
	int i;

	LIS_DEBUG_FUNC_IN;

	if( solver && solver->work )
	{
		for(i=0;i<solver->worklen;i++) lis_vector_destroy(solver->work[i]);
		lis_free(solver->work);
		solver->work    = NULL;
		solver->worklen = 0;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_destroy"
int lis_solver_destroy(LIS_SOLVER solver)
{
	LIS_DEBUG_FUNC_IN;

	if( solver )
	{
		lis_solver_work_destroy(solver);
		lis_vector_destroy(solver->d);
		if( solver->At ) lis_matrix_destroy(solver->At);
		if( solver->residual ) lis_free(solver->residual);
		lis_free(solver);
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}


#undef __FUNC__
#define __FUNC__ "lis_solve"
int lis_solve(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_SOLVER solver)
{
        int		err;
	LIS_PRECON      precon;

	LIS_DEBUG_FUNC_IN;

	solver->A = A;

	/* create preconditioner */

	if( solver->options[LIS_OPTIONS_PRECON] < 0 || solver->options[LIS_OPTIONS_PRECON] > LIS_PRECONNAME_MAX )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_PRECON is %d (Set between 0 to %d)\n",solver->options[LIS_OPTIONS_PRECON], LIS_PRECONNAME_MAX);
		return LIS_ERR_ILL_ARG;
	}

	err = lis_precon_create(solver, &precon);
	if( err )
	{
		lis_solver_work_destroy(solver);
		solver->retcode = err;
		return err;
	}

	/* Core Kernel of lis_solve() */
	lis_solve_kernel(A, b, x, solver, precon);

	lis_precon_destroy(precon);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solve_kernel"
int lis_solve_kernel(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_SOLVER solver, LIS_PRECON precon)
{
	int			nsolver, precon_type, maxiter;
	int			err;
	LIS_SCALAR	*residual;
	LIS_VECTOR	xx;

	int output;
	int scale;
	int conv_cond;
	int precision,is_use_at,storage,block;
	int i,n,np;
	double p_c_times, p_i_times,itimes;
	LIS_SCALAR nrm2,tol,tol_w;
	LIS_VECTOR t;
	LIS_VECTOR bb;
	LIS_MATRIX AA,B;
	LIS_MATRIX At;
	char buf[64];

	LIS_DEBUG_FUNC_IN;

	nsolver     = solver->options[LIS_OPTIONS_SOLVER];
	precon_type = solver->options[LIS_OPTIONS_PRECON];
	maxiter     = solver->options[LIS_OPTIONS_MAXITER];
	output      = solver->options[LIS_OPTIONS_OUTPUT];
	scale       = solver->options[LIS_OPTIONS_SCALE];
	precision   = solver->options[LIS_OPTIONS_PRECISION];
	is_use_at   = solver->options[LIS_OPTIONS_USE_AT];
	storage     = solver->options[LIS_OPTIONS_STORAGE];
	block       = solver->options[LIS_OPTIONS_STORAGE_BLOCK];
	conv_cond   = solver->options[LIS_OPTIONS_CONV_COND];
	tol         = solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN];
	tol_w       = solver->params[LIS_PARAMS_RESID_WEIGHT-LIS_OPTIONS_LEN];
	solver->precision = precision;

	if( nsolver < 1 || nsolver > LIS_SOLVERS_LEN )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_SOLVER is %d (Set between 1 to %d)\n",nsolver, LIS_SOLVERS_LEN);
		return LIS_ERR_ILL_ARG;
	}
	if( precon_type < 0 || precon_type > precon_register_type )
	{
		LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_PRECON is %d (Set between 0 to %d)\n",precon_type, precon_register_type-1);
		return LIS_ERR_ILL_ARG;
	}
	if( maxiter<0 )
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_OPTIONS_MAXITER(=%d) is less than 0\n",maxiter);
		return LIS_ERR_ILL_ARG;
	}
	#ifdef USE_MPI
	if( precon_type == LIS_PRECON_TYPE_SAAMG  && solver->A->nprocs < 2)
	{
		LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter A->nprocs (=%d) is less than 2 (Set more than 1 when using parallel version of SAAMG)\n",solver->A->nprocs);
		return LIS_ERR_ILL_ARG;
	}
	#endif
	#ifdef USE_QUAD_PRECISION
		if( precision==LIS_PRECISION_QUAD && lis_solver_execute_quad[nsolver]==NULL )
		{
			LIS_SETERR1(LIS_ERR_NOT_IMPLEMENTED,"Quad precision solver %s is not implemented\n",lis_solvername[nsolver]);
			return LIS_ERR_NOT_IMPLEMENTED;
		}
		else if( precision==LIS_PRECISION_SWITCH && lis_solver_execute_switch[nsolver]==NULL )
		{
			LIS_SETERR1(LIS_ERR_NOT_IMPLEMENTED,"Switch solver %s is not implemented\n",lis_solvername[nsolver]);
			return LIS_ERR_NOT_IMPLEMENTED;
		}
		if( solver->options[LIS_OPTIONS_SWITCH_MAXITER]==-1 )
		{
			solver->options[LIS_OPTIONS_SWITCH_MAXITER] = maxiter;
		}
	#endif

	err = lis_solver_check_params[nsolver](solver);
	if( err )
	{
		solver->retcode = err;
		return err;
	}
	/* end parameter check */

	solver->A        = A;
	solver->b        = b;

	/* create initial vector */
	#ifndef USE_QUAD_PRECISION
		err = lis_vector_duplicate(A,&xx);
	#else
		if( precision==LIS_PRECISION_DOUBLE )
		{
			err = lis_vector_duplicate(A,&xx);
		}
		else
		{
			err = lis_vector_duplicateex(LIS_PRECISION_QUAD,A,&xx);
		}
	#endif
	if( err )
	{
		solver->retcode = err;
		return err;
	}
	if( solver->options[LIS_OPTIONS_INITGUESS_ZEROS] )
	{
	  if( output ) lis_printf(A->comm,"initial vector x = 0\n");
		#ifndef USE_QUAD_PRECISION
			lis_vector_set_all(0.0,xx);
		#else
			if( precision==LIS_PRECISION_DOUBLE )
			{
				lis_vector_set_all(0.0,xx);
			}
			else
			{
				lis_vector_set_allex_nm(0.0,xx);
			}
		#endif
	}
	else
	{
	  if( output ) lis_printf(A->comm,"initial vector x = user defined\n"); 
		#ifndef USE_QUAD_PRECISION
			lis_vector_copy(x,xx);
		#else
			if( precision==LIS_PRECISION_DOUBLE )
			{
				lis_vector_copy(x,xx);
			}
			else
			{
				lis_vector_copyex_nm(x,xx);
			}
		#endif
	}

	/* create residual history vector */
	if( solver->residual ) lis_free(solver->residual);
	residual = (LIS_SCALAR *)lis_malloc((maxiter+2)*sizeof(LIS_SCALAR),"lis_solve::residual");
	if( residual==NULL )
	{
		LIS_SETERR_MEM((maxiter+2)*sizeof(LIS_SCALAR));
		lis_vector_destroy(xx);
		solver->retcode = err;
		return err;
	}
	residual[0] = 1.0;


	n       = A->n;
	np      = A->np;
	t       = NULL;
	At      = NULL;


	p_c_times = lis_wtime();
	if( precon_type==LIS_PRECON_TYPE_IS )
	{
		if( solver->d==NULL )
		{
			err = lis_vector_duplicate(A,&solver->d);
			if( err )
			{
				return err;
			}
		}
		if( !A->is_scaled )
		{
			lis_matrix_scaling(A,b,solver->d,LIS_SCALE_JACOBI);
		}
		else if( !b->is_scaled )
		{
			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for(i=0;i<n;i++)
			{
				b->value[i] = b->value[i]*solver->d->value[i];
			}
		}
		if( nsolver >= LIS_SOLVER_JACOBI && nsolver <= LIS_SOLVER_SOR )
		{
			solver->options[LIS_OPTIONS_ISLEVEL] = 0;
		}
	}
	else if( nsolver >= LIS_SOLVER_JACOBI && nsolver <= LIS_SOLVER_SOR && precon_type!=LIS_PRECON_TYPE_NONE )
	{
		if( solver->d==NULL )
		{
			err = lis_vector_duplicate(A,&solver->d);
			if( err )
			{
				return err;
			}
		}
		if( !A->is_scaled )
		{
			lis_matrix_scaling(A,b,solver->d,LIS_SCALE_JACOBI);
		}
	}
	else if( scale )
	{
		if( storage==LIS_MATRIX_BSR && scale==LIS_SCALE_JACOBI )
		{
			if( A->matrix_type!=LIS_MATRIX_BSR )
			{
				err = lis_matrix_duplicate(A,&B);
				if( err ) return err;
				lis_matrix_set_blocksize(B,block,block,NULL,NULL);
				lis_matrix_set_type(B,storage);
				err = lis_matrix_convert(A,B);
				if( err ) return err;
				lis_matrix_storage_destroy(A);
				lis_matrix_DLU_destroy(A);
				lis_matrix_diag_destroy(A->WD);
				if( A->l2g_map ) lis_free( A->l2g_map );
				if( A->commtable ) lis_commtable_destroy( A->commtable );
				if( A->ranges ) lis_free( A->ranges );
				err = lis_matrix_copy_struct(B,A);
				if( err ) return err;
				lis_free(B);
			}
			err = lis_matrix_split(A);
			if( err ) return err;
			err = lis_matrix_diag_duplicate(A->D,&solver->WD);
			if( err ) return err;
			lis_matrix_diag_copy(A->D,solver->WD);
			lis_matrix_diag_inverse(solver->WD);
			lis_matrix_bscaling_bsr(A,solver->WD);
			lis_vector_duplicate(A,&t);
			lis_matrix_diag_matvec(solver->WD,b,t);
			lis_vector_copy(t,b);
			lis_vector_destroy(t);
			t = NULL;
		}
		else
		{
			if( solver->d==NULL )
			{
				err = lis_vector_duplicate(A,&solver->d);
				if( err )
				{
					return err;
				}
			}
			if( scale==LIS_SCALE_JACOBI && nsolver==LIS_SOLVER_CG )
			{
				scale = LIS_SCALE_SYMM_DIAG;
			}
			if( !A->is_scaled )
			{
				lis_matrix_scaling(A,b,solver->d,scale);
			}
			else if( !b->is_scaled )
			{
				#ifdef _OPENMP
				#pragma omp parallel for
				#endif
				for(i=0;i<n;i++)
				{
					b->value[i] = b->value[i]*solver->d->value[i];
				}
			}
		}
	}

/*	precon_type = precon->precon_type;*/
	if( precon_type==LIS_PRECON_TYPE_IS )
	{
		if( nsolver < LIS_SOLVER_JACOBI || nsolver > LIS_SOLVER_SOR )
		{
			AA = solver->A;
			bb = solver->b;
		}
		else
		{
			AA = precon->A;
			bb = precon->Pb;
		}
	}
	else
	{
		AA = A;
		bb = b;
	}

	p_c_times = lis_wtime() - p_c_times;
	itimes = lis_wtime();

	/* Matrix Convert */
	solver->A  = AA;
	solver->b  = bb;
	err = lis_matrix_convert_self(solver);
	if( err )
	{
		lis_vector_destroy(xx);
		lis_solver_work_destroy(solver);
		lis_free(residual);
		solver->retcode = err;
		return err;
	}
	block = solver->A->bnr;

	if( A->my_rank==0 )
	{
	  if( output ) printf("precision : %s\n", lis_precisionname[precision]); 
	  if( output ) printf("solver    : %s %d\n", lis_solvername[nsolver],nsolver); 
		switch( precon_type )
		{
		case LIS_PRECON_TYPE_ILU:
			i = solver->options[LIS_OPTIONS_FILL];
			if( A->matrix_type==LIS_MATRIX_BSR || A->matrix_type==LIS_MATRIX_VBR )
			{
			  if( output ) sprintf(buf,"Block %s(%d)",lis_preconname[precon_type],i); 
			}
			else
			{
			  if( output ) sprintf(buf,"%s(%d)",lis_preconname[precon_type],i); 
			}
			break;
		default:
		  if( output ) sprintf(buf,"%s",lis_preconname[precon_type]); 
			break;
		}
		if( solver->options[LIS_OPTIONS_ADDS] && precon_type )
		{
		  if( output ) printf("precon    : %s + additive schwarz\n", buf); 
		}
		else
		{
		  if( output ) printf("precon    : %s\n", buf); 
		}
	}
	switch(conv_cond)
	{
	case LIS_CONV_COND_NRM2_R:
	case LIS_CONV_COND_NRM2_B:
		if( A->my_rank==0 )
		{
		  if( output ) ("CONV_COND : ||r||_2 <= %6.1e*||r_0||_2\n", tol); 
		}
		break;
	case LIS_CONV_COND_NRM1_B:
		lis_vector_nrm1(b,&nrm2);
		nrm2 = nrm2*tol_w + tol;
		if( A->my_rank==0 )
		{
		  if( output ) printf("conv_cond : ||r||_1 <= %6.1e*||b||_1 + %6.1e = %6.1e\n", tol_w,tol,nrm2);
		}
		break;
	}
	if( A->my_rank==0 )
	{
		if( AA->matrix_type==LIS_MATRIX_BSR || AA->matrix_type==LIS_MATRIX_BSC )
		{
		  if( output ) printf("storage   : %s(%d x %d)\n", lis_storagename[AA->matrix_type-1],block,block); 
		}
		else
		{
		  if( output ) printf("storage   : %s\n", lis_storagename[AA->matrix_type-1]); 
		}
	}


	/* create work vector */
	err = lis_solver_malloc_work[nsolver](solver); 
	if( err )
	{
		lis_vector_destroy(xx);
		lis_precon_destroy(precon);
		solver->retcode = err;
		return err;
	}
	if( nsolver==LIS_SOLVER_BICG && is_use_at )
	{
	  if( output ) lis_printf(A->comm,"Use At\n"); 
		lis_matrix_duplicate(AA,&At);
		lis_matrix_set_type(At,LIS_USE_AT_TYPE[AA->matrix_type]);
		lis_matrix_convert(AA,At);
		solver->At = At;
	}

	solver->x        = xx;
	solver->xx       = x;
	solver->precon   = precon;
	solver->residual = residual;

	/* execute solver */
	#ifndef USE_QUAD_PRECISION
		err = lis_solver_execute[nsolver](solver);
	#else
		if( precision==LIS_PRECISION_DOUBLE )
		{
			err = lis_solver_execute[nsolver](solver);
		}
		else if( precision==LIS_PRECISION_QUAD )
		{
			err = lis_solver_execute_quad[nsolver](solver);
		}
		else if( precision==LIS_PRECISION_SWITCH )
		{
			err = lis_solver_execute_switch[nsolver](solver);
		}
	#endif
	solver->retcode = err;

	if( scale==LIS_SCALE_SYMM_DIAG && precon_type!=LIS_PRECON_TYPE_IS)
	{
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(i=0;i<n;i++)
		{
			x->value[i] = xx->value[i]*solver->d->value[i];
		}
	}
	else
	{
		#ifndef USE_QUAD_PRECISION
			lis_vector_copy(xx,x);
		#else
			if( precision==LIS_PRECISION_DOUBLE )
			{
				lis_vector_copy(xx,x);
			}
			else
			{
				lis_vector_copyex_mn(xx,x);
			}
		#endif
	}
	itimes = lis_wtime() - itimes - solver->ptimes;
	p_i_times = solver->ptimes;
	solver->ptimes = p_c_times + p_i_times;
	solver->p_c_times = p_c_times;
	solver->p_i_times = p_i_times;
	solver->times  = solver->ptimes + itimes;
	solver->itimes = itimes;
	lis_solver_work_destroy(solver);
	lis_vector_duplicate(A,&t);
	xx->precision = LIS_PRECISION_DEFAULT;
	lis_matvec(A,xx,t);
	lis_vector_xpay(b,-1.0,t);
	if( scale==LIS_SCALE_SYMM_DIAG && precon_type!=LIS_PRECON_TYPE_IS)
	{
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(i=0;i<n;i++)
		{
			t->value[i] = t->value[i]/solver->d->value[i];
		}
	}
	lis_vector_nrm2(t,&nrm2);
	solver->resid = nrm2;
	if( A->my_rank==0 )
	{
		if( err )
		{
		  if( output ) printf("lis_solve : %s(code=%d)\n\n",lis_returncode[err],err); 

		}
		else
		{
		  if( output ) printf("lis_solve : normal end\n\n"); 
		}
	}
	if( precision==LIS_PRECISION_DOUBLE )
	{
		solver->iter2 = solver->iter;
	}
	else if( precision==LIS_PRECISION_QUAD )
	{
		solver->iter2 = 0;
	}


	lis_vector_destroy(t);
/*	lis_vector_destroy(d);*/
	lis_vector_destroy(xx);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_initial_residual"
int lis_solver_get_initial_residual(LIS_SOLVER solver, LIS_PRECON M, LIS_VECTOR t, LIS_VECTOR r, LIS_SCALAR *bnrm2)
{
	int			output,conv;
	#ifdef USE_QUAD_PRECISION
		int	i;
	#endif
	LIS_MATRIX	A;
	LIS_VECTOR	x,b,p,xx;
	LIS_SCALAR	nrm2;
	LIS_REAL	tol,tol_w,tol_switch;

	LIS_DEBUG_FUNC_IN;

	A  = solver->A;
	b  = solver->b;
	x  = solver->x;
	xx = solver->xx;
	output     = solver->options[LIS_OPTIONS_OUTPUT];
	conv       = solver->options[LIS_OPTIONS_CONV_COND];
	tol        = solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN];
	tol_w      = solver->params[LIS_PARAMS_RESID_WEIGHT-LIS_OPTIONS_LEN];
	tol_switch = solver->params[LIS_PARAMS_SWITCH_RESID-LIS_OPTIONS_LEN];


	/* Initial Residual */
	if( M==NULL )
	{
		p = r;
	}
	else
	{
		p = t;
	}

	if( !solver->options[LIS_OPTIONS_INITGUESS_ZEROS] )
	{
		#ifndef USE_QUAD_PRECISION
			lis_matvec(A,x,p);           /* p = Ax    */
			lis_vector_xpay(b,-1,p);     /* p = b - p */
		#else
			if( solver->precision==LIS_PRECISION_DOUBLE )
			{
				lis_matvec(A,x,p);           /* p = Ax    */
				lis_vector_xpay(b,-1,p);     /* p = b - p */
			}
			else
			{
				lis_matvec(A,xx,p);           /* p = Ax    */
				lis_vector_xpay(b,-1,p);     /* p = b - p */
				
				#ifdef _OPENMP
				#pragma omp parallel for private(i)
				#endif
				for(i=0;i<A->n;i++)
				{
					p->value_lo[i] = 0.0;
				}
				
			}
		#endif
	}
	else
	{
		#ifndef USE_QUAD_PRECISION
			lis_vector_copy(b,p);
		#else
			if( solver->precision==LIS_PRECISION_DOUBLE )
			{
				lis_vector_copy(b,p);
			}
			else
			{
				lis_vector_copyex_nm(b,p);
			}
		#endif
	}

	switch(conv)
	{
	case LIS_CONV_COND_NRM2_R:
		lis_vector_nrm2(p,&nrm2);
		*bnrm2 = nrm2;
		solver->tol = tol;
		solver->tol_switch = tol_switch;
		break;
	case LIS_CONV_COND_NRM2_B:
		lis_vector_nrm2(p,&nrm2);
		lis_vector_nrm2(b,bnrm2);
		solver->tol = tol;
		solver->tol_switch = tol_switch;
		break;
	case LIS_CONV_COND_NRM1_B:
		lis_vector_nrm1(p,&nrm2);
		lis_vector_nrm1(b,bnrm2);
		solver->tol = *bnrm2*tol_w + tol;
		solver->tol_switch = *bnrm2*tol_w + tol_switch;
		break;
	}
	if( *bnrm2 == 0.0 )
	{
		*bnrm2 = 1.0;
	}
	else
	{
		*bnrm2 = 1.0 / *bnrm2;
	}
	solver->bnrm = *bnrm2;
	nrm2 = nrm2 * *bnrm2;

	if( output && (r->precision==LIS_PRECISION_QUAD && solver->precision!=LIS_PRECISION_SWITCH) )
	{
		if( output & LIS_PRINT_MEM ) solver->residual[0] = nrm2;
		if( output & LIS_PRINT_OUT && A->my_rank==0 ) printf("iter: %5d  residual = %e\n", 0, nrm2); 
	}
	if( nrm2 <= solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN] )
	{
		solver->retcode = LIS_SUCCESS;
		solver->iter    = 1;
		solver->resid   = nrm2; 
		LIS_DEBUG_FUNC_OUT;
		return LIS_FAILS;
	}

	if( M!=NULL )
	{
		/* r = M^-1 * p */
		lis_psolve(solver, p, r);
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_optionC"
int lis_solver_set_optionC(LIS_SOLVER solver)
{
	LIS_ARGS	p;

	LIS_DEBUG_FUNC_IN;

	p = cmd_args->next;
	while( p!=cmd_args )
	{
		lis_solver_set_option2(p->arg1,p->arg2,solver);
		p = p->next;
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option"
int lis_solver_set_option(char *text, LIS_SOLVER solver)
{
	LIS_ARGS	args,p;

	LIS_DEBUG_FUNC_IN;

	lis_text2args(text,&args);
	p = args->next;
	while( p!=args )
	{
		lis_solver_set_option2(p->arg1,p->arg2,solver);
		p = p->next;
	}
	lis_args_free(args);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option2"
int lis_solver_set_option2(char* arg1, char *arg2, LIS_SOLVER solver)
{
	int i;

	LIS_DEBUG_FUNC_IN;

	for(i=0;i<LIS_SOLVER_OPTION_LEN;i++)
	{
		if( strcmp(arg1, LIS_SOLVER_OPTNAME[i])==0 )
		{
			switch( LIS_SOLVER_OPTACT[i] )
			{
			case LIS_OPTIONS_FILE:
				break;
			case LIS_OPTIONS_HELP:
				break;
			case LIS_OPTIONS_VER:
				break;
			case LIS_OPTIONS_SOLVER:
				lis_solver_set_option_solver(arg2,solver);
				break;
			case LIS_OPTIONS_PRECON:
				lis_solver_set_option_precon(arg2,solver);
				break;
			case LIS_OPTIONS_SCALE:
				lis_solver_set_option_scale(arg2,solver);
				break;
			case LIS_OPTIONS_OUTPUT:
				lis_solver_set_option_print(arg2,solver);
				break;
			case LIS_OPTIONS_PSOLVER:
				lis_solver_set_option_psolver(arg2,solver);
				break;
			case LIS_OPTIONS_PPRECON:
				lis_solver_set_option_pprecon(arg2,solver);
				break;
			case LIS_OPTIONS_INITGUESS_ZEROS:
				lis_solver_set_option_truefalse(arg2,LIS_OPTIONS_INITGUESS_ZEROS,solver);
				break;
			case LIS_OPTIONS_ADDS:
				lis_solver_set_option_truefalse(arg2,LIS_OPTIONS_ADDS,solver);
				break;
			case LIS_OPTIONS_PRECISION:
				lis_solver_set_option_precision(arg2,LIS_OPTIONS_PRECISION,solver);
				break;
			case LIS_OPTIONS_USE_AT:
				lis_solver_set_option_truefalse(arg2,LIS_OPTIONS_USE_AT,solver);
				break;
			case LIS_OPTIONS_SAAMG_UNSYM:
				lis_solver_set_option_truefalse(arg2,LIS_OPTIONS_SAAMG_UNSYM,solver);
				if (solver->options[LIS_OPTIONS_SAAMG_UNSYM])
				  {
				    solver->params[LIS_PARAMS_SAAMG_THETA  -LIS_OPTIONS_LEN] = 0.12;
				  }
				break;
			case LIS_OPTIONS_STORAGE:
				lis_solver_set_option_storage(arg2,solver);
				break;
			case LIS_OPTIONS_CONV_COND:
				lis_solver_set_option_conv_cond(arg2,solver);
				break;
			default:
				if( LIS_SOLVER_OPTACT[i] < LIS_OPTIONS_LEN )
				{
					sscanf(arg2, "%d", &solver->options[LIS_SOLVER_OPTACT[i]]);
				}
				else
				{
					sscanf(arg2, "%lg", &solver->params[LIS_SOLVER_OPTACT[i]-LIS_OPTIONS_LEN]);
				}
				break;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_solver"
int lis_solver_set_option_solver(char *argv, LIS_SOLVER solver)
{
	int  i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='9' )
	{
		sscanf(argv, "%d", &solver->options[LIS_OPTIONS_SOLVER]);
	}
	else
	{
		for(i=0;i<LIS_SOLVER_LEN;i++)
		{
			if( strcmp(argv,lis_solver_atoi[i])==0 )
			{
				solver->options[LIS_OPTIONS_SOLVER] = i+1;
				break;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_psolver"
int lis_solver_set_option_psolver(char *argv, LIS_SOLVER solver)
{
	int  i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='9' )
	{
		sscanf(argv, "%d", &solver->options[LIS_OPTIONS_PSOLVER]);
	}
	else
	{
		for(i=0;i<LIS_SOLVER_LEN;i++)
		{
			if( strcmp(argv,lis_solver_atoi[i])==0 )
			{
				solver->options[LIS_OPTIONS_PSOLVER] = i+1;
				break;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_precon"
int lis_solver_set_option_precon(char *argv, LIS_SOLVER solver)
{
	int  i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='9' )
	{
		sscanf(argv, "%d", &solver->options[LIS_OPTIONS_PRECON]);
	}
	else
	{
		for(i=0;i<LIS_PRECON_TYPE_LEN;i++)
		{
			if( strcmp(argv,lis_precon_atoi[i])==0 )
			{
				solver->options[LIS_OPTIONS_PRECON] = i;
				LIS_DEBUG_FUNC_OUT;
				return LIS_SUCCESS;
			}
		}
		for(i=0;i<precon_register_type-LIS_PRECON_TYPE_USERDEF;i++)
		{
			if( strcmp(argv,precon_register_top[i].name)==0 )
			{
				solver->options[LIS_OPTIONS_PRECON] = i+LIS_PRECON_TYPE_USERDEF;
				LIS_DEBUG_FUNC_OUT;
				return LIS_SUCCESS;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_pprecon"
int lis_solver_set_option_pprecon(char *argv, LIS_SOLVER solver)
{
	int  i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='9' )
	{
		sscanf(argv, "%d", &solver->options[LIS_OPTIONS_PPRECON]);
	}
	else
	{
		for(i=0;i<LIS_PRECON_TYPE_LEN;i++)
		{
			if( strcmp(argv,lis_precon_atoi[i])==0 )
			{
				solver->options[LIS_OPTIONS_PPRECON] = i;
				LIS_DEBUG_FUNC_OUT;
				return LIS_SUCCESS;
			}
		}
		for(i=0;i<precon_register_type-LIS_PRECON_TYPE_USERDEF;i++)
		{
			if( strcmp(argv,precon_register_top[i].name)==0 )
			{
				solver->options[LIS_OPTIONS_PPRECON] = i+LIS_PRECON_TYPE_USERDEF;
				LIS_DEBUG_FUNC_OUT;
				return LIS_SUCCESS;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_print"
int lis_solver_set_option_print(char *argv, LIS_SOLVER solver)
{
	int  i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='3' )
	{
		sscanf(argv, "%d", &solver->options[LIS_OPTIONS_OUTPUT]);
	}
	else
	{
		for(i=0;i<LIS_PRINT_LEN;i++)
		{
			if( strcmp(argv,lis_print_atoi[i])==0 )
			{
				solver->options[LIS_OPTIONS_OUTPUT] = i;
				break;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_scale"
int lis_solver_set_option_scale(char *argv, LIS_SOLVER solver)
{
	int  i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='2' )
	{
		sscanf(argv, "%d", &solver->options[LIS_OPTIONS_SCALE]);
	}
	else
	{
		for(i=0;i<LIS_SCALE_LEN;i++)
		{
			if( strcmp(argv,lis_scale_atoi[i])==0 )
			{
				solver->options[LIS_OPTIONS_SCALE] = i;
				break;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_truefalse"
int lis_solver_set_option_truefalse(char *argv, int opt, LIS_SOLVER solver)
{
	int  i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='1' )
	{
		sscanf(argv, "%d", &solver->options[opt]);
	}
	else
	{
		for(i=0;i<LIS_TRUEFALSE_LEN;i++)
		{
			if( strcmp(argv,lis_truefalse_atoi[i])==0 )
			{
				solver->options[opt] = i;
				break;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_precision"
int lis_solver_set_option_precision(char *argv, int opt, LIS_SOLVER solver)
{
	int  i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='1' )
	{
		sscanf(argv, "%d", &solver->options[opt]);
	}
	else
	{
		for(i=0;i<LIS_PRECISION_LEN;i++)
		{
			if( strcmp(argv,lis_precision_atoi[i])==0 )
			{
				solver->options[opt] = i;
				break;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_storage"
int lis_solver_set_option_storage(char *argv, LIS_SOLVER solver)
{
	int  i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='9' )
	{
		sscanf(argv, "%d", &solver->options[LIS_OPTIONS_STORAGE]);
	}
	else
	{
		for(i=0;i<LIS_STORAGE_LEN;i++)
		{
			if( strcmp(argv,lis_storage_atoi[i])==0 )
			{
				solver->options[LIS_OPTIONS_STORAGE] = i+1;
				break;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_conv_cond"
int lis_solver_set_option_conv_cond(char *argv, LIS_SOLVER solver)
{
	int  i;

	LIS_DEBUG_FUNC_IN;

	if( argv[0]>='0' && argv[0]<='3' )
	{
		sscanf(argv, "%d", &solver->options[LIS_OPTIONS_CONV_COND]);
	}
	else
	{
		for(i=0;i<LIS_CONV_COND_LEN;i++)
		{
			if( strcmp(argv,lis_conv_cond_atoi[i])==0 )
			{
				solver->options[LIS_OPTIONS_CONV_COND] = i;
				break;
			}
		}
	}
	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_iters"
int lis_solver_get_iters(LIS_SOLVER solver, int *iters)
{
	LIS_DEBUG_FUNC_IN;

	*iters = solver->iter;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_itersex"
int lis_solver_get_itersex(LIS_SOLVER solver, int *iters, int *iters_double, int *iters_quad)
{
	LIS_DEBUG_FUNC_IN;

	*iters = solver->iter;
	*iters_double = solver->iter2;
	*iters_quad = solver->iter - solver->iter2;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_time"
int lis_solver_get_time(LIS_SOLVER solver, double *times)
{
	LIS_DEBUG_FUNC_IN;

	*times  = solver->times;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_timeex"
int lis_solver_get_timeex(LIS_SOLVER solver, double *times, double *itimes, double *ptimes, double *p_c_times, double *p_i_times)
{
	LIS_DEBUG_FUNC_IN;

	*times  = solver->times;
	if( itimes ) *itimes = solver->itimes;
	if( ptimes ) *ptimes = solver->ptimes;
	if( p_c_times ) *p_c_times = solver->p_c_times;
	if( p_i_times ) *p_i_times = solver->p_i_times;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_residualnorm"
int lis_solver_get_residualnorm(LIS_SOLVER solver, LIS_REAL *residual)
{
	LIS_DEBUG_FUNC_IN;

	*residual  = solver->resid;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_solver"
int lis_solver_get_solver(LIS_SOLVER solver, int *nsol)
{
	LIS_DEBUG_FUNC_IN;

	*nsol = solver->options[LIS_OPTIONS_SOLVER];

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_status"
int lis_solver_get_status(LIS_SOLVER solver, int *status)
{
	LIS_DEBUG_FUNC_IN;

	*status = solver->retcode;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_rhistory"
int lis_solver_get_rhistory(LIS_SOLVER solver, LIS_VECTOR v)
{
	int		i,n,maxiter;

	maxiter = solver->iter+1;
        if( solver->retcode!=LIS_SUCCESS )
	  {
	    maxiter--;
	  }
	n = _min(v->n,maxiter);
	for(i=0;i<n;i++)
	{
		v->value[i] = solver->residual[i];
	}
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_get_solvername"
int lis_get_solvername(int solver, char *solvername)
{
	LIS_DEBUG_FUNC_IN;

	if( solver < 1 || solver > LIS_SOLVERS_LEN )
	{
		solvername = NULL;
		return LIS_FAILS;
	}
	strcpy(solvername,lis_solvername[solver]);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_residual_nrm2_r"
int lis_solver_get_residual_nrm2_r(LIS_VECTOR r, LIS_SOLVER solver, LIS_REAL *res)
{
	LIS_DEBUG_FUNC_IN;

	lis_vector_nrm2(r,res);
	*res = *res * solver->bnrm;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_get_residual_nrm1_b"
int lis_solver_get_residual_nrm1_b(LIS_VECTOR r, LIS_SOLVER solver, LIS_REAL *res)
{
	LIS_DEBUG_FUNC_IN;

	lis_vector_nrm1(r,res);

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_shadowresidual"
int lis_solver_set_shadowresidual(LIS_SOLVER solver, LIS_VECTOR r0, LIS_VECTOR rs0)
{
    unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
	int i,n,resid;

	LIS_DEBUG_FUNC_IN;

	resid = solver->options[LIS_OPTIONS_INIT_SHADOW_RESID];
	if( resid==LIS_RANDOM )
	{
		n     = solver->A->n;
		init_by_array(init, length);

		#ifdef USE_QUAD_PRECISION
		if( solver->precision==LIS_PRECISION_DEFAULT )
		#endif
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0;i<n;i++)
			{
				rs0->value[i] = genrand_real1();
			}
		}
		#ifdef USE_QUAD_PRECISION
		else
		{
			#ifdef _OPENMP
			#pragma omp parallel for private(i)
			#endif
			for(i=0;i<n;i++)
			{
				rs0->value[i]    = genrand_real1();
				rs0->value_lo[i] = 0.0;
			}
		}
		#endif
	}
	else
	{
		#ifdef USE_QUAD_PRECISION
		if( solver->precision==LIS_PRECISION_DEFAULT )
		#endif
		{
			lis_vector_copy(r0,rs0);
		}
		#ifdef USE_QUAD_PRECISION
		else
		{
			lis_vector_copyex_mm(r0,rs0);
		}
		#endif
	}

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

