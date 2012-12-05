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


#ifndef __LIS_PRECON_H__
#define __LIS_PRECON_H__


#define lis_psolve(solver,b,x)		lis_psolve_xxx[solver->precon->precon_type](solver,b,x)
#define lis_psolvet(solver,b,x)		lis_psolvet_xxx[solver->precon->precon_type](solver,b,x)




#ifdef __cplusplus
extern "C"
{
#endif
	extern int lis_precon_init(LIS_PRECON precon);
	extern int lis_precon_create(LIS_SOLVER solver, LIS_PRECON *precon);
	extern int lis_precon_destroy(LIS_PRECON precon);

	extern LIS_PRECON_CREATE_XXX lis_precon_create_xxx[];
	extern LIS_PSOLVE_XXX lis_psolve_xxx[];
	extern LIS_PSOLVET_XXX lis_psolvet_xxx[];
	/*******************/
	/* NONE            */
	/*******************/
	extern int lis_precon_create_none(LIS_SOLVER solver, LIS_PRECON precon);
	extern int lis_psolve_none(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x);
	extern int lis_psolvet_none(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x);
	/*******************/
	/* Jacobi          */
	/*******************/
	extern int lis_precon_create_jacobi(LIS_SOLVER solver, LIS_PRECON precon);
	extern int lis_psolve_jacobi(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x);
	extern int lis_psolvet_jacobi(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x);
	extern int lis_precon_create_bjacobi(LIS_SOLVER solver, LIS_PRECON precon);
	extern int lis_psolve_bjacobi(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x);
	extern int lis_psolvet_bjacobi(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x);
	/*******************/
	/* ILU             */
	/*******************/
	extern int lis_precon_create_iluk(LIS_SOLVER solver, LIS_PRECON precon);
	extern int lis_psolve_iluk_crs(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	extern int lis_psolve_iluk_bsr(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	extern int lis_psolve_iluk_vbr(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	extern int lis_psolvet_iluk_crs(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	extern int lis_psolvet_iluk_bsr(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	/*******************/
	/* SSOR            */
	/*******************/
	extern int lis_precon_create_ssor(LIS_SOLVER solver, LIS_PRECON precon);
	extern int lis_psolve_ssor(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	extern int lis_psolvet_ssor(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	/*******************/
	/* Hybrid          */
	/*******************/
	extern int lis_precon_create_hybrid(LIS_SOLVER solver, LIS_PRECON precon);
	extern int lis_psolve_hybrid(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	extern int lis_psolvet_hybrid(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	/*******************/
	/* I+S             */
	/*******************/
	extern int lis_precon_create_is(LIS_SOLVER solver, LIS_PRECON precon);
	extern int lis_precon_create_is_crs(LIS_SOLVER solver, LIS_PRECON precon);
	extern int lis_psolve_is(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	extern int lis_psolvet_is(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	/*******************/
	/* SAINV           */
	/*******************/
	extern int lis_precon_create_sainv(LIS_SOLVER solver, LIS_PRECON precon);
	extern int lis_precon_create_sainv_crs(LIS_SOLVER solver, LIS_PRECON precon);
	extern int lis_psolve_sainv(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	extern int lis_psolvet_sainv(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	/*******************/
	/* SAAMG           */
	/*******************/
	extern int lis_precon_create_saamg(LIS_SOLVER solver, LIS_PRECON precon);
	extern int lis_psolve_saamg(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	extern int lis_psolvet_saamg(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	/*******************/
	/* Crout ILU       */
	/*******************/
	extern int lis_precon_create_iluc(LIS_SOLVER solver, LIS_PRECON precon);
	extern int lis_precon_create_iluc_crs(LIS_SOLVER solver, LIS_PRECON precon);
	extern int lis_precon_create_iluc_bsr(LIS_SOLVER solver, LIS_PRECON precon);
	extern int lis_psolve_iluc(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	extern int lis_psolvet_iluc(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	extern int lis_psolve_iluc_bsr(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	/*******************/
	/* ILUT            */
	/*******************/
	extern int lis_precon_create_ilut(LIS_SOLVER solver, LIS_PRECON precon);
	extern int lis_precon_create_ilut_crs(LIS_SOLVER solver, LIS_PRECON precon);
	extern int lis_precon_create_ilut_bsr(LIS_SOLVER solver, LIS_PRECON precon);
	extern int lis_psolve_ilut_crs(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	extern int lis_psolve_ilut_bsr(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	extern int lis_psolvet_ilut_crs(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	extern int lis_psolvet_ilut_bsr(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	/********************/
	/* Additive Schwarz */
	/********************/
	extern int lis_precon_create_adds(LIS_SOLVER solver, LIS_PRECON precon);
	extern int lis_psolve_adds(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);
	extern int lis_psolvet_adds(LIS_SOLVER solver, LIS_VECTOR B, LIS_VECTOR X);

#ifdef __cplusplus
}
#endif
#endif
