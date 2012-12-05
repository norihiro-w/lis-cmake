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


extern LIS_SOLVER_MALLOC_WORK lis_solver_malloc_work[];
extern LIS_SOLVER_EXECUTE lis_solver_execute[];

#undef __FUNC__
#define __FUNC__ "lis_precon_create_hybrid"
int lis_precon_create_hybrid(LIS_SOLVER solver, LIS_PRECON precon)
{
	int			nsolver, precon_type, maxiter,precision;
	int			err;
	LIS_SCALAR	*residual;
	LIS_VECTOR	xx;
	LIS_SOLVER	psolver;
	LIS_MATRIX	A;
	LIS_PRECON	pprecon;

	LIS_DEBUG_FUNC_IN;


	A           = solver->A;

	err = lis_solver_create(&psolver);
	if( err )
	{
		return err;
	}

	psolver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN] = solver->params[LIS_PARAMS_PRESID-LIS_OPTIONS_LEN];
	psolver->params[LIS_PARAMS_W-LIS_OPTIONS_LEN]     = solver->params[LIS_PARAMS_PW-LIS_OPTIONS_LEN];
	psolver->options[LIS_OPTIONS_MAXITER]             = solver->options[LIS_OPTIONS_PMAXITER];
	psolver->options[LIS_OPTIONS_ELL]                 = solver->options[LIS_OPTIONS_PELL];
	psolver->options[LIS_OPTIONS_RESTART]             = solver->options[LIS_OPTIONS_PRESTART];
	psolver->options[LIS_OPTIONS_OUTPUT]              = 0;
	psolver->options[LIS_OPTIONS_SOLVER]              = solver->options[LIS_OPTIONS_PSOLVER];
	psolver->options[LIS_OPTIONS_PRECON]              = solver->options[LIS_OPTIONS_PPRECON];
	psolver->options[LIS_OPTIONS_INITGUESS_ZEROS]     = solver->options[LIS_OPTIONS_INITGUESS_ZEROS];
	psolver->options[LIS_OPTIONS_PRECISION]           = solver->options[LIS_OPTIONS_PRECISION];
	psolver->A                                        = solver->A;
	psolver->At                                       = solver->At;
	psolver->precision                                = solver->precision;

	nsolver     = psolver->options[LIS_OPTIONS_SOLVER];
	precon_type = psolver->options[LIS_OPTIONS_PRECON];
	maxiter     = psolver->options[LIS_OPTIONS_MAXITER];
	precision   = psolver->options[LIS_OPTIONS_PRECISION];
	A           = psolver->A;



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

	/* create residual history vector */
	residual = (LIS_SCALAR *)lis_malloc((maxiter+2)*sizeof(LIS_SCALAR),"lis_precon_create_hybrid::residual");
	if( residual==NULL )
	{
		LIS_SETERR_MEM((maxiter+2)*sizeof(LIS_SCALAR));
		lis_vector_destroy(xx);
		solver->retcode = err;
		return err;
	}

	/* create preconditioner */
	err = lis_precon_create(psolver, &pprecon);
	if( err )
	{
		lis_vector_destroy(xx);
		lis_solver_work_destroy(psolver);
		lis_free(residual);
		solver->retcode = err;
		return err;
	}

	/* create work vector */
	err = lis_solver_malloc_work[nsolver](psolver);
	if( err )
	{
		lis_vector_destroy(xx);
		lis_precon_destroy(pprecon);
		solver->retcode = err;
		return err;
	}


	psolver->x        = xx;
	psolver->precon   = pprecon;
	psolver->residual = residual;
	precon->solver    = psolver;

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_psolve_hybrid"
int lis_psolve_hybrid(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
	LIS_VECTOR		xx;
	LIS_SOLVER		solver2;
	int				nsolver;
	LIS_PRECON	precon;

	/*
	 *  Mx = b
	 *  M  = A
	 */

	LIS_DEBUG_FUNC_IN;

	precon      = solver->precon;
	solver2     = precon->solver;
	xx          = precon->solver->x;
	nsolver     = solver2->options[LIS_OPTIONS_SOLVER];
	solver2->b  = B;

	if( solver2->options[LIS_OPTIONS_INITGUESS_ZEROS] )
	{
		#ifdef USE_QUAD_PRECISION
		if( solver->precision==LIS_PRECISION_DEFAULT )
		{
		#endif
			lis_vector_set_all(0,xx);
		#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_vector_set_allex_nm(0,xx);
		}
		#endif
	}
	else
	{
		#ifdef USE_QUAD_PRECISION
		if( solver->precision==LIS_PRECISION_DEFAULT )
		{
		#endif
			lis_vector_copy(B,xx);
		#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_vector_copyex_mm(B,xx);
		}
		#endif
	}

	/* execute solver */
	lis_solver_execute[nsolver](solver2);
	#ifdef USE_QUAD_PRECISION
	if( solver->precision==LIS_PRECISION_DEFAULT )
	{
	#endif
		lis_vector_copy(solver2->x,X);
	#ifdef USE_QUAD_PRECISION
	}
	else
	{
		lis_vector_copyex_mm(solver2->x,X);
	}
	#endif

	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_psolvet_hybrid"
int lis_psolvet_hybrid(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X)
{
	LIS_VECTOR		xx;
	LIS_SOLVER		solver2;
	int				nsolver;
	LIS_PRECON		precon;

	/*
	 *  Mx = b
	 *  M  = A
	 */

	LIS_DEBUG_FUNC_IN;

	precon      = solver->precon;
	solver2     = precon->solver;
	xx          = precon->solver->x;
	nsolver     = solver2->options[LIS_OPTIONS_SOLVER];
	solver2->b  = B;
	LIS_MATVEC  = lis_matvect;
	LIS_MATVECT = lis_matvec;

	if( solver2->options[LIS_OPTIONS_INITGUESS_ZEROS] )
	{
		#ifdef USE_QUAD_PRECISION
		if( solver->precision==LIS_PRECISION_DEFAULT )
		{
		#endif
			lis_vector_set_all(0,xx);
		#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_vector_set_allex_nm(0,xx);
		}
		#endif
	}
	else
	{
		#ifdef USE_QUAD_PRECISION
		if( solver->precision==LIS_PRECISION_DEFAULT )
		{
		#endif
			lis_vector_copy(B,xx);
		#ifdef USE_QUAD_PRECISION
		}
		else
		{
			lis_vector_copyex_mm(B,xx);
		}
		#endif
	}

	/* execute solver */
	lis_solver_execute[nsolver](solver2);
	#ifdef USE_QUAD_PRECISION
	if( solver->precision==LIS_PRECISION_DEFAULT )
	{
	#endif
		lis_vector_copy(solver2->x,X);
	#ifdef USE_QUAD_PRECISION
	}
	else
	{
		lis_vector_copyex_mm(solver2->x,X);
	}
	#endif
	LIS_MATVEC  = lis_matvec;
	LIS_MATVECT = lis_matvect;


	LIS_DEBUG_FUNC_OUT;
	return LIS_SUCCESS;
}

