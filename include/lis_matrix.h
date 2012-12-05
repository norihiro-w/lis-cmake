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


#ifndef __LIS_MATRIX_H__
#define __LIS_MATRIX_H__


#define LIS_MATRIX_CRS_STR		"crs"
#define LIS_MATRIX_CCS_STR		"ccs"
#define LIS_MATRIX_MSR_STR		"msr"
#define LIS_MATRIX_DIA_STR		"dia"
#define LIS_MATRIX_ELL_STR		"ell"
#define LIS_MATRIX_JDS_STR		"jds"
#define LIS_MATRIX_BSR_STR		"bsr"
#define LIS_MATRIX_BSC_STR		"bsc"
#define LIS_MATRIX_VBR_STR		"vbr"
#define LIS_MATRIX_DNS_STR		"dns"
#define LIS_MATRIX_COO_STR		"coo"
#define LIS_MATRIX_TJD_STR		"tjd"

#define LIS_MATRIX_CHECK_ALL			0
#define LIS_MATRIX_CHECK_SIZE			1
#define LIS_MATRIX_CHECK_NULL			2
#define LIS_MATRIX_CHECK_TYPE			3
#define LIS_MATRIX_CHECK_NOT_ASSEMBLED	4
#define LIS_MATRIX_CHECK_SET			5

#define LIS_MATRIX_OPTION_CALL_BY		0

#define LIS_CALL_BY_REFERENCE			0
#define LIS_CALL_BY_VALUE				1

#define LIS_MATRIX_W_ANNZ				10

#ifdef __cplusplus
extern "C"
{
#endif
	extern int lis_matrix_init(LIS_MATRIX *Amat);
	extern int lis_matrix_check(LIS_MATRIX A, int level);
	extern int lis_matrix_storage_destroy(LIS_MATRIX Amat);
	extern int lis_matrix_set_destroyflag(LIS_MATRIX A, int flag);
	extern int lis_matrix_get_destroyflag(LIS_MATRIX A, int *flag);
	extern int lis_matrix_LU_create(LIS_MATRIX A);
	extern int lis_matrix_DLU_destroy(LIS_MATRIX Amat);
	extern int lis_matrix_unset(LIS_MATRIX A);
	extern int lis_matrix_copy_struct(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_convert_self(LIS_SOLVER solver);
	/*******************/
	/* Operations      */
	/*******************/
	extern int lis_matrix_split(LIS_MATRIX A);
	extern int lis_matrix_merge(LIS_MATRIX A);
	extern int lis_matrix_solve(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, int flag);
	extern int lis_matrix_solvet(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, int flag);
	extern int lis_matrix_copyDLU(LIS_MATRIX Ain, LIS_MATRIX_DIAG *D, LIS_MATRIX *L, LIS_MATRIX *U);
        extern int lis_matrix_shift_diagonal(LIS_MATRIX A, LIS_SCALAR shift);
        extern int lis_matrix_shift_diagonal_crs(LIS_MATRIX A, LIS_SCALAR shift);
	extern int lis_matrix_shift_diagonal_ccs(LIS_MATRIX A, LIS_SCALAR shift);
	extern int lis_matrix_shift_diagonal_msr(LIS_MATRIX A, LIS_SCALAR shift);
	extern int lis_matrix_shift_diagonal_dia(LIS_MATRIX A, LIS_SCALAR shift);
	extern int lis_matrix_shift_diagonal_ell(LIS_MATRIX A, LIS_SCALAR shift);
	extern int lis_matrix_shift_diagonal_jds(LIS_MATRIX A, LIS_SCALAR shift);
	extern int lis_matrix_shift_diagonal_bsr(LIS_MATRIX A, LIS_SCALAR shift);
	extern int lis_matrix_shift_diagonal_bsc(LIS_MATRIX A, LIS_SCALAR shift);
	extern int lis_matrix_shift_diagonal_dns(LIS_MATRIX A, LIS_SCALAR shift);
	extern int lis_matrix_shift_diagonal_coo(LIS_MATRIX A, LIS_SCALAR shift);
	extern int lis_matrix_shift_diagonal_vbr(LIS_MATRIX A, LIS_SCALAR shift);

	/*******************/
	/* Array           */
	/*******************/
	extern void lis_array_LUdecomp(int n, LIS_SCALAR *a);
	extern void lis_array_invGauss(int n, LIS_SCALAR *a);
	extern void lis_array_matmat(int n, LIS_SCALAR *a, LIS_SCALAR *b, LIS_SCALAR *c, int op);
	extern void lis_array_matvec(int n, LIS_SCALAR *a, LIS_SCALAR *b, LIS_SCALAR *c, int op);
	extern void lis_array_matvect(int n, LIS_SCALAR *a, LIS_SCALAR *b, LIS_SCALAR *c, int op);
	extern void lis_array_nrm2(int n, LIS_SCALAR *v, LIS_SCALAR *nrm2);
	extern void lis_array_nrm1(int n, LIS_SCALAR *v, LIS_SCALAR *nrm1);
	extern void lis_array_dot(int n, LIS_SCALAR *v, LIS_SCALAR *dot);
	extern void lis_array_matmat2(int m, int n, int k, LIS_SCALAR *a, int lda, LIS_SCALAR *b, int ldb, LIS_SCALAR *c, int ldc, int op);
	extern void lis_array_matvec2(int m, int n, LIS_SCALAR *a, int lda, LIS_SCALAR *b, LIS_SCALAR *c, int op);
	extern void lis_array_matinv(int n, LIS_SCALAR *a, LIS_SCALAR *b, LIS_SCALAR *c);
	extern void lis_array_invvec(int n, LIS_SCALAR *a, LIS_SCALAR *x, LIS_SCALAR *y);
	extern void lis_array_invvect(int n, LIS_SCALAR *a, LIS_SCALAR *x, LIS_SCALAR *y);
	extern void lis_array_solve(int n, LIS_SCALAR *aa, LIS_SCALAR *b, LIS_SCALAR *x, LIS_SCALAR *a);
        extern void lis_array_set_all(int n, LIS_SCALAR alpha, LIS_SCALAR *v);
        extern void lis_array_scale(int n, LIS_SCALAR alpha, LIS_SCALAR *v);
        extern void lis_array_dot2(int n, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *alpha);
        extern void lis_array_copy(int n, LIS_SCALAR *x, LIS_SCALAR *y);
        extern void lis_array_axpyz(int n, LIS_SCALAR alpha, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *z);
        extern void lis_array_power(int n, LIS_SCALAR *a, LIS_SCALAR *x, LIS_SCALAR *mu, int maxiter, LIS_REAL tol, LIS_REAL *err);
        extern int lis_array_qr(int n, LIS_SCALAR *x, LIS_SCALAR *q, LIS_SCALAR *r);
        extern int lis_array_cgs(int n, LIS_SCALAR *x, LIS_SCALAR *q, LIS_SCALAR *r);
        extern int lis_array_mgs(int n, LIS_SCALAR *x, LIS_SCALAR *q, LIS_SCALAR *r);
	/*******************/
	/* Diagonal Matrix */
	/*******************/
	extern int lis_matrix_diag_init(LIS_MATRIX_DIAG *D);
	extern int lis_matrix_diag_check(LIS_MATRIX_DIAG D, int level);
	extern int lis_matrix_diag_create(int local_n, int global_n, LIS_Comm comm, LIS_MATRIX_DIAG *D);
	extern int lis_matrix_diag_destroy(LIS_MATRIX_DIAG D);
	extern int lis_matrix_diag_duplicate(LIS_MATRIX_DIAG Din, LIS_MATRIX_DIAG *Dout);
	extern int lis_matrix_diag_duplicateM(LIS_MATRIX Ain, LIS_MATRIX_DIAG *Dout);
	extern int lis_matrix_diag_get_range(LIS_MATRIX_DIAG D, int *is, int *ie);
	extern int lis_matrix_diag_get_size(LIS_MATRIX_DIAG D, int *local_n, int *global_n);
	extern int lis_matrix_diag_set_blocksize(LIS_MATRIX_DIAG D, int bn, int *bns);
	extern int lis_matrix_diag_copy(LIS_MATRIX_DIAG X, LIS_MATRIX_DIAG Y);
	extern int lis_matrix_diag_scale(LIS_SCALAR alpha, LIS_MATRIX_DIAG D);
	extern int lis_matrix_diag_inverse(LIS_MATRIX_DIAG D);
	extern int lis_matrix_diag_print(LIS_MATRIX_DIAG D);
	extern int lis_matrix_diag_mallocM(LIS_MATRIX A, LIS_SCALAR **diag);
	extern int lis_matrix_diag_matvec(LIS_MATRIX_DIAG D, LIS_VECTOR X, LIS_VECTOR Y);
	extern int lis_matrix_diag_matvect(LIS_MATRIX_DIAG D, LIS_VECTOR X, LIS_VECTOR Y);

	/*******************/
	/* CRS             */
	/*******************/
	extern int lis_matrix_setDLU_crs(int nnzl, int nnzu, LIS_SCALAR *diag, int *lptr, int *lindex, LIS_SCALAR *lvalue, int *uptr, int *uindex, LIS_SCALAR *uvalue, LIS_MATRIX A);
	extern int lis_matrix_elements_copy_crs(int n, int *ptr, int *index, LIS_SCALAR *value, int *o_ptr, int *o_index, LIS_SCALAR *o_value);
	extern int lis_matrix_copy_crs(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_copyDLU_crs(LIS_MATRIX Ain, LIS_MATRIX_DIAG *D, LIS_MATRIX *L, LIS_MATRIX *U);
	extern int lis_matrix_get_diagonal_crs(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_crs(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_symm_crs(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_sort_crs(LIS_MATRIX A);
	extern int lis_matrix_normf_crs(LIS_MATRIX A, LIS_SCALAR *nrm);
	extern int lis_matrix_transpose_crs(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_split_crs(LIS_MATRIX A);
	extern int lis_matrix_merge_crs(LIS_MATRIX A);
	extern int lis_matrix_solve_crs(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	extern int lis_matrix_solvet_crs(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	extern int lis_matrix_split2_crs(LIS_MATRIX A);
	/*******************/
	/* CCS             */
	/*******************/
	extern int lis_matrix_setDLU_ccs(int nnzl, int nnzu, LIS_SCALAR *diag, int *lptr, int *lindex, LIS_SCALAR *lvalue, int *uptr, int *uindex, LIS_SCALAR *uvalue, LIS_MATRIX A);
	extern int lis_matrix_elements_copy_ccs(int n, int *ptr, int *index, LIS_SCALAR *value, int *o_ptr, int *o_index, LIS_SCALAR *o_value);
	extern int lis_matrix_copy_ccs(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_get_diagonal_ccs(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_ccs(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_symm_ccs(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_sort_ccs(LIS_MATRIX A);
	extern int lis_matrix_normf_ccs(LIS_MATRIX A, LIS_SCALAR *nrm);
	extern int lis_matrix_transpose_ccs(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_split_ccs(LIS_MATRIX A);
	extern int lis_matrix_merge_ccs(LIS_MATRIX A);
	extern int lis_matrix_solve_ccs(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	extern int lis_matrix_solvet_ccs(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	extern int lis_matrix_convert_crs2ccs(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_convert_ccs2crs(LIS_MATRIX Ain, LIS_MATRIX Aout);
	/*******************/
	/* BSR             */
	/*******************/
	extern int lis_matrix_setDLU_bsr(int bnr, int bnc, int lbnnz, int ubnnz, LIS_MATRIX_DIAG D, int *lbptr, int *lbindex, LIS_SCALAR *lvalue, int *ubptr, int *ubindex, LIS_SCALAR *uvalue, LIS_MATRIX A);
	extern int lis_matrix_elements_copy_bsr(int n, int bnr, int bnc, int bnnz, int *ptr, int *index, LIS_SCALAR *value, int *o_ptr, int *o_index, LIS_SCALAR *o_value);
	extern int lis_matrix_copy_bsr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_convert_crs2bsr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_convert_bsr2crs(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_get_diagonal_bsr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_bsr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_symm_bsr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_bscaling_bsr(LIS_MATRIX A, LIS_MATRIX_DIAG D);
	extern int lis_matrix_split_bsr(LIS_MATRIX A);
	extern int lis_matrix_merge_bsr(LIS_MATRIX A);
	extern int lis_matrix_solve_bsr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	extern int lis_matrix_solvet_bsr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	extern int lis_matrix_sort_bsr(LIS_MATRIX A);
	/*******************/
	/* MSR             */
	/*******************/
	extern int lis_matrix_setDLU_msr(int lnnz, int unnz, int lndz, int undz, LIS_SCALAR *diag, int *lindex, LIS_SCALAR *lvalue, int *uindex, LIS_SCALAR *uvalue, LIS_MATRIX A);
	extern int lis_matrix_elements_copy_msr(int n, int *index, LIS_SCALAR *value, int *o_index, LIS_SCALAR *o_value);
	extern int lis_matrix_copy_msr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_merge_msr(LIS_MATRIX A);
	extern int lis_matrix_split_msr(LIS_MATRIX A);
	extern int lis_matrix_get_diagonal_msr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_msr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_symm_msr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_solve_msr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	extern int lis_matrix_solvet_msr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	extern int lis_matrix_convert_crs2msr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_convert_msr2crs(LIS_MATRIX Ain, LIS_MATRIX Aout);
	/*******************/
	/* ELL             */
	/*******************/
	extern int lis_matrix_setDLU_ell(int lmaxnzr, int umaxnzr, LIS_SCALAR *diag, int *lindex, LIS_SCALAR *lvalue, int *uindex, LIS_SCALAR *uvalue, LIS_MATRIX A);
	extern int lis_matrix_elements_copy_ell(int n, int maxnzr, int *index, LIS_SCALAR *value, int *o_index, LIS_SCALAR *o_value);
	extern int lis_matrix_copy_ell(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_merge_ell(LIS_MATRIX A);
	extern int lis_matrix_split_ell(LIS_MATRIX A);
	extern int lis_matrix_get_diagonal_ell(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_ell(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_symm_ell(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_solve_ell(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	extern int lis_matrix_solvet_ell(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	extern int lis_matrix_convert_crs2ell(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_convert_ell2crs(LIS_MATRIX Ain, LIS_MATRIX Aout);
	/*******************/
	/* JDS             */
	/*******************/
	extern int lis_matrix_setDLU_jds(int lnnz, int unnz, int lmaxnzr, int umaxnzr, LIS_SCALAR *diag, int *lperm, int *lptr, int *lindex, LIS_SCALAR *lvalue, int *uperm, int *uptr, int *uindex, LIS_SCALAR *uvalue, LIS_MATRIX A);
	extern int lis_matrix_elements_copy_jds(int n, int maxnzr, int *perm, int *ptr, int *index, LIS_SCALAR *value, int *o_perm, int *o_ptr, int *o_index, LIS_SCALAR *o_value);
	extern int lis_matrix_copy_jds(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_merge_jds(LIS_MATRIX A);
	extern int lis_matrix_split_jds(LIS_MATRIX A);
	extern int lis_matrix_get_diagonal_jds(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_jds(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_symm_jds(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_solve_jds(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	extern int lis_matrix_solvet_jds(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	extern int lis_matrix_convert_crs2jds(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_convert_jds2crs(LIS_MATRIX Ain, LIS_MATRIX Aout);
	/*******************/
	/* DIA             */
	/*******************/
	extern int lis_matrix_setDLU_dia(int lnnd, int unnd, LIS_SCALAR *diag, int *lindex, LIS_SCALAR *lvalue, int *uindex, LIS_SCALAR *uvalue, LIS_MATRIX A);
	extern int lis_matrix_elements_copy_dia(int n, int nnd, int *index, LIS_SCALAR *value, int *o_index, LIS_SCALAR *o_value);
	extern int lis_matrix_copy_dia(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_merge_dia(LIS_MATRIX A);
	extern int lis_matrix_split_dia(LIS_MATRIX A);
	extern int lis_matrix_get_diagonal_dia(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_dia(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_symm_dia(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_solve_dia(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	extern int lis_matrix_solvet_dia(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	extern int lis_matrix_convert_crs2dia(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_convert_dia2crs(LIS_MATRIX Ain, LIS_MATRIX Aout);
	/*******************/
	/* BSC             */
	/*******************/
	extern int lis_matrix_setDLU_bsc(int bnr, int bnc, int lbnnz, int ubnnz, LIS_MATRIX_DIAG D, int *lbptr, int *lbindex, LIS_SCALAR *lvalue, int *ubptr, int *ubindex, LIS_SCALAR *uvalue, LIS_MATRIX A);
	extern int lis_matrix_elements_copy_bsc(int n, int bnr, int bnc, int bnnz, int *ptr, int *index, LIS_SCALAR *value, int *o_ptr, int *o_index, LIS_SCALAR *o_value);
	extern int lis_matrix_copy_bsc(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_convert_ccs2bsc(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_convert_bsc2crs(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_get_diagonal_bsc(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_bsc(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_symm_bsc(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_split_bsc(LIS_MATRIX A);
	extern int lis_matrix_merge_bsc(LIS_MATRIX A);
	extern int lis_matrix_solve_bsc(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	extern int lis_matrix_solvet_bsc(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	/*******************/
	/* VBR             */
	/*******************/
	extern int lis_matrix_elements_copy_vbr(int n, int nr, int nc, int bnnz, int *row, int *col, int *ptr, int *bptr, int *bindex, LIS_SCALAR *value, int *o_row, int *o_col, int *o_ptr, int *o_bptr, int *o_bindex, LIS_SCALAR *o_value);
	extern int lis_matrix_copy_vbr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_convert_crs2vbr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_convert_vbr2crs(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_get_diagonal_vbr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_vbr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_symm_vbr(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_split_vbr(LIS_MATRIX A);
	extern int lis_matrix_merge_vbr(LIS_MATRIX A);
	extern int lis_matrix_solve_vbr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	extern int lis_matrix_solvet_vbr(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	/*******************/
	/* DNS             */
	/*******************/
	extern int lis_matrix_setDLU_dns(LIS_SCALAR *diag, LIS_SCALAR *lvalue, LIS_SCALAR *uvalue, LIS_MATRIX A);
	extern int lis_matrix_elements_copy_dns(int n, int gn, LIS_SCALAR *value, LIS_SCALAR *o_value);
	extern int lis_matrix_copy_dns(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_merge_dns(LIS_MATRIX A);
	extern int lis_matrix_split_dns(LIS_MATRIX A);
	extern int lis_matrix_get_diagonal_dns(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_dns(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_symm_dns(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_solve_dns(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	extern int lis_matrix_solvet_dns(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR X, int flag);
	extern int lis_matrix_convert_crs2dns(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_convert_dns2crs(LIS_MATRIX Ain, LIS_MATRIX Aout);
	/*******************/
	/* COO             */
	/*******************/
	extern int lis_matrix_copy_coo(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_get_diagonal_coo(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_coo(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_scaling_symm_coo(LIS_MATRIX A, LIS_SCALAR d[]);
	extern int lis_matrix_sort_coo(LIS_MATRIX A);
	extern int lis_matrix_normf_coo(LIS_MATRIX A, LIS_SCALAR *nrm);
	extern int lis_matrix_transpose_coo(LIS_MATRIX Ain, LIS_MATRIX *Aout);
	extern int lis_matrix_split_coo(LIS_MATRIX A);
	extern int lis_matrix_merge_coo(LIS_MATRIX A);
	extern int lis_matrix_convert_crs2coo(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_convert_coo2crs(LIS_MATRIX Ain, LIS_MATRIX Aout);
	/*******************/
	/* RCO             */
	/*******************/
	extern int lis_matrix_create_rco(int local_n, int global_n, LIS_Comm comm, int annz, int *nnz, LIS_MATRIX *Amat);
	extern int lis_matrix_malloc_rco(int n, int nnz[], int **row, int ***index, LIS_SCALAR ***value);
	extern int lis_matrix_realloc_rco(int row, int nnz, int ***index, LIS_SCALAR ***value);
	extern int lis_matrix_convert_rco2crs(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_convert_rco2bsr(LIS_MATRIX Ain, LIS_MATRIX Aout);
	extern int lis_matrix_convert_rco2ccs(LIS_MATRIX Ain, LIS_MATRIX Aout);
	/*******************/
	/* ILU             */
	/*******************/
	extern int lis_matrix_ilu_create(int n, int bs, LIS_MATRIX_ILU *A);
	extern int lis_matrix_ilu_setCR(LIS_MATRIX_ILU A);
	extern int lis_matrix_ilu_setVR(LIS_MATRIX_ILU A);
	extern int lis_matrix_ilu_destroy(LIS_MATRIX_ILU A);
	extern int lis_matrix_ilu_premalloc(int nnzrow, LIS_MATRIX_ILU A);
	extern int lis_matrix_ilu_realloc(int row, int nnz, LIS_MATRIX_ILU A);
	extern int lis_matvec_ilu(LIS_MATRIX A, LIS_MATRIX_ILU LU, LIS_VECTOR X, LIS_VECTOR Y);
	extern int lis_matvect_ilu(LIS_MATRIX A, LIS_MATRIX_ILU LU, LIS_VECTOR X, LIS_VECTOR Y);


#ifdef __cplusplus
}
#endif

#endif
