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


#ifndef __LIS_ESOLVER_H__
#define __LIS_ESOLVER_H__

typedef int (*LIS_ESOLVER_CHECK_PARAMS)(LIS_ESOLVER esolver);
typedef int (*LIS_ESOLVER_MALLOC_WORK)(LIS_ESOLVER esolver);
typedef int (*LIS_ESOLVER_EXECUTE)(LIS_ESOLVER esolver);

#ifdef __cplusplus
extern "C"
{
#endif
        extern LIS_ESOLVER_EXECUTE lis_esolver_execute[];
	extern int lis_esolver_work_destroy(LIS_ESOLVER esolver);
	extern int lis_esolver_set_option2(char* arg1, char *arg2, LIS_ESOLVER esolver);
	extern int lis_esolver_set_option_esolver(char *argv, LIS_ESOLVER esolver);
	extern int lis_esolver_set_option_iesolver(char *argv, LIS_ESOLVER esolver);
	extern int lis_esolver_set_option_print(char *argv, LIS_ESOLVER esolver);
	extern int lis_esolver_set_option_truefalse(char *argv, int opt, LIS_ESOLVER esolver);
	extern int lis_esolver_set_option_eprecision(char *argv, int opt, LIS_ESOLVER esolver);
	extern int lis_esolver_set_option_storage(char *argv, LIS_ESOLVER esolver);
	extern int lis_esolver_get_residual(LIS_VECTOR r, LIS_ESOLVER esolver, LIS_REAL *res);

        /*******************/
	/* Power Iteration */
        /*******************/
	extern int lis_epi(LIS_ESOLVER esolver);
	extern int lis_epi_quad(LIS_ESOLVER esolver);
	extern int lis_epi_check_params(LIS_ESOLVER esolver);
	extern int lis_epi_malloc_work(LIS_ESOLVER esolver);
        /*********************/
	/* Inverse Iteration */
	/*********************/
	extern int lis_eii(LIS_ESOLVER esolver);
	extern int lis_eii_quad(LIS_ESOLVER esolver);
	extern int lis_eii_check_params(LIS_ESOLVER esolver);
	extern int lis_eii_malloc_work(LIS_ESOLVER esolver);
        /*********************************/
	/* Approximate Inverse Iteration */
        /*********************************/
	extern int lis_eaii(LIS_ESOLVER esolver);
	extern int lis_eaii_check_params(LIS_ESOLVER esolver);
	extern int lis_eaii_malloc_work(LIS_ESOLVER esolver);
        /*******************************/
	/* Rayleigh Quotient Iteration */
        /*******************************/
	extern int lis_erqi(LIS_ESOLVER esolver);
	extern int lis_erqi_quad(LIS_ESOLVER esolver);
	extern int lis_erqi_check_params(LIS_ESOLVER esolver);
	extern int lis_erqi_malloc_work(LIS_ESOLVER esolver);
	/**********************/
	/* Subspace Iteration */
	/**********************/
	extern int lis_esi(LIS_ESOLVER esolver);
	extern int lis_esi_quad(LIS_ESOLVER esolver);
	extern int lis_esi_check_params(LIS_ESOLVER esolver);
	extern int lis_esi_malloc_work(LIS_ESOLVER esolver);
	/*********************/
	/* Lanczos Iteration */
	/*********************/
	extern int lis_eli(LIS_ESOLVER esolver);
	extern int lis_eli_quad(LIS_ESOLVER esolver);
	extern int lis_eli_check_params(LIS_ESOLVER esolver);
	extern int lis_eli_malloc_work(LIS_ESOLVER esolver);
        /**********************/
	/* Conjugate Gradient */
	/**********************/
	extern int lis_ecg(LIS_ESOLVER esolver);
	extern int lis_ecg_check_params(LIS_ESOLVER esolver);
	extern int lis_ecg_malloc_work(LIS_ESOLVER esolver);
        /**********************/
	/* Conjugate Residual */
	/**********************/
	extern int lis_ecr(LIS_ESOLVER esolver);
	extern int lis_ecr_check_params(LIS_ESOLVER esolver);
	extern int lis_ecr_malloc_work(LIS_ESOLVER esolver);

#ifdef __cplusplus
}
#endif
#endif
