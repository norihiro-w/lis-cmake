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


#ifndef __LIS_MATRIX_IO_H__
#define __LIS_MATRIX_IO_H__

#define BUFSIZE				1024
#define MM_BANNER			"%%MatrixMarket"
#define MM_MTX				"matrix"
#define MM_VEC				"vector"
#define MM_FMT				"coordinate"
#define MM_TYPE_REAL		"real"
#define MM_TYPE_GENERAL		"general"
#define MM_TYPE_SYMM		"symmetric"
#define MM_REAL				0
#define MM_GENERAL			0
#define MM_SYMM				1

#define LISBanner			"#LIS"
#define ITBLBanner			"#ITBL"

#define LIS_INPUT_N				0
#define LIS_INPUT_NNZ			1
#define LIS_INPUT_NDZ			2
#define LIS_INPUT_NND			3
#define LIS_INPUT_MAXNZR			4
#define LIS_INPUT_NR				5
#define LIS_INPUT_NC				6
#define LIS_INPUT_BNNZ			7
#define LIS_INPUT_BNR			8
#define LIS_INPUT_BNC			9
#define LIS_INPUT_NNZL			10
#define LIS_INPUT_NNZU			11

#define LIS_INPUT_ROW			0
#define LIS_INPUT_COL			1
#define LIS_INPUT_PTR			2
#define LIS_INPUT_INDEX			3
#define LIS_INPUT_BPTR			4
#define LIS_INPUT_BINDEX			5
#define LIS_INPUT_VALUE			6
#define LIS_INPUT_SOLUTION		7
#define LIS_INPUT_RHS			8
#define LIS_INPUT_DIAG			9
#define LIS_INPUT_LPTR			10
#define LIS_INPUT_LCOL			11
#define LIS_INPUT_LINDEX			12
#define LIS_INPUT_LVALUE			13
#define LIS_INPUT_UPTR			14
#define LIS_INPUT_UCOL			15
#define LIS_INPUT_UINDEX			16
#define LIS_INPUT_UVALUE			17

typedef struct
{
	char		filename[1024];
	int			fileformat;
	int			matrix_type_in;
	int			matrix_type_out;
	int			is_matrix_split;
	int			is_array_zero_cut;
	int			is_solution;
	int			is_rhs;
	int			array_base_ptr;
	int			array_base_index;
	int			precision;
	int			param_len;
	int			array_len;
	int			param_order[9];
	int			array_order[9];
} LIS_MATRIX_INPUT_OPTION;

typedef struct
{
	int			i;
	LIS_SCALAR	value;
} LIS_MM_VECFMT;

typedef struct
{
	int			i;
	int			j;
	LIS_SCALAR	value;
} LIS_MM_MATFMT;

#ifdef __cplusplus
extern "C"
{
#endif

	/****************************/
	/* Matrix INPUT             */
	/****************************/
	extern int lis_input_mm(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, FILE *file);
		extern int lis_input_mm_crs(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, FILE *file);
/*		extern int lis_input_mm_ccs(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, FILE *file, LIS_Comm comm);
		extern int lis_input_mm_msr(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, FILE *file, LIS_Comm comm);
		extern int lis_input_mm_dia(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, FILE *file, LIS_Comm comm);
		extern int lis_input_mm_ell(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, FILE *file, LIS_Comm comm);
		extern int lis_input_mm_jds(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, FILE *file, LIS_Comm comm);
		extern int lis_input_mm_dns(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, FILE *file, LIS_Comm comm);
		extern int lis_input_mm_coo(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, FILE *file, LIS_Comm comm);*/
	extern int lis_input_hb(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, FILE *file);
		extern int lis_input_hb_crs(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, FILE *file);
/*	extern int lis_input_mmm(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, FILE *file, LIS_Comm comm, int matrix_type, ...);
	extern int lis_input_free(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, LIS_MATRIX_INPUT_OPTION option, FILE *file, LIS_Comm comm, int matrix_type, ...);
	extern int lis_input_lis(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, char *filename, FILE *file, LIS_Comm comm, int matrix_type, ...);
		extern int lis_input_lis_crs(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, int mode, FILE *file, LIS_Comm comm);
		extern int lis_input_lis_ccs(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, int mode, FILE *file, LIS_Comm comm);
		extern int lis_input_lis_msr(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, int mode, FILE *file, LIS_Comm comm);
		extern int lis_input_lis_dia(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, int mode, FILE *file, LIS_Comm comm);
		extern int lis_input_lis_ell(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, int mode, FILE *file, LIS_Comm comm);
		extern int lis_input_lis_jds(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, int mode, FILE *file, LIS_Comm comm);
		extern int lis_input_lis_bsr(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, int mode, FILE *file, LIS_Comm comm);
		extern int lis_input_lis_bsc(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, int mode, FILE *file, LIS_Comm comm);
		extern int lis_input_lis_vbr(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, int mode, FILE *file, LIS_Comm comm);
		extern int lis_input_lis_dns(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, int mode, FILE *file, LIS_Comm comm);
		extern int lis_input_lis_coo(LIS_MATRIX *A, LIS_VECTOR *b, LIS_VECTOR *x, int mode, FILE *file, LIS_Comm comm);*/
	/****************************/
	/* Matrix OUTPUT            */
	/****************************/
	extern int lis_output_mm_crs(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, int format, char *path);
	extern int lis_output_mm_ccs(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, int format, char *path);
/*	extern int lis_output_mm_msr(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x,  char *path);
	extern int lis_output_mm_dia(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x,  char *path);
	extern int lis_output_mm_ell(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x,  char *path);
	extern int lis_output_mm_jds(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x,  char *path);
	extern int lis_output_mm_bsr(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x,  char *path);
	extern int lis_output_mm_bsc(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x,  char *path);
	extern int lis_output_mm_vbr(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x,  char *path);
	extern int lis_output_mm_coo(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x,  char *path);
	extern int lis_output_mm_dns(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x,  char *path);
	extern int lis_output_lis_crs(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, int mode, char *path);
	extern int lis_output_lis_ccs(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, int mode, char *path);
	extern int lis_output_lis_msr(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, int mode, char *path);
	extern int lis_output_lis_dia(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, int mode, char *path);
	extern int lis_output_lis_ell(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, int mode, char *path);
	extern int lis_output_lis_jds(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, int mode, char *path);
	extern int lis_output_lis_bsr(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, int mode, char *path);
	extern int lis_output_lis_bsc(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, int mode, char *path);
	extern int lis_output_lis_vbr(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, int mode, char *path);
	extern int lis_output_lis_dns(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, int mode, char *path);
	extern int lis_output_lis_coo(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, int mode, char *path);*/
	/****************************/
	/* Vector OUTPUT            */
	/****************************/
	extern int lis_output_vector_plain(LIS_VECTOR v, char *path);
	extern int lis_output_vector_mm(LIS_VECTOR v, char *path);
	extern int lis_output_vector_lis_ascii(LIS_VECTOR v, char *path);
	/****************************/
	/* Vector INPUT             */
	/****************************/
	extern int lis_input_vector_plain(LIS_VECTOR v, FILE *file);
	extern int lis_input_vector_mm(LIS_VECTOR v, FILE *file);
	extern int lis_input_vector_lis(LIS_VECTOR v, char *filename, FILE *file);
	extern int lis_input_vector_lis_ascii(LIS_VECTOR v, FILE *file);

	extern int lis_fscan_int(int n, FILE *file, int val[], int origin);
	extern int lis_fscan_double(int n, FILE *file, double val[]);


#ifdef __cplusplus
}
#endif
#endif
