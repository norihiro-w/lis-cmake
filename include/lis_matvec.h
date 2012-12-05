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


#ifndef __LIS_MATVEC_H__
#define __LIS_MATVEC_H__

#define LIS_MATVEC_SENDRECV \
		if( (A->np+A->pad)>(X->np+X->pad) ) \
		{ \
			X->value = (LIS_SCALAR *)lis_realloc(X->value,(A->np+A->pad)*sizeof(LIS_SCALAR)); \
			if( X->value==NULL ) \
			{ \
				LIS_SETERR_MEM((A->np+A->pad)*sizeof(LIS_SCALAR)); \
				return LIS_ERR_OUT_OF_MEMORY; \
			} \
			X->np  = A->np; \
			X->pad = A->pad; \
			x = X->value; \
		} \
		lis_send_recv(A->commtable,x)

#define LIS_MATVEC_REALLOC \
		if( (A->np+A->pad)>(X->np+X->pad) ) \
		{ \
			X->value = (LIS_SCALAR *)lis_realloc(X->value,(A->np+A->pad)*sizeof(LIS_SCALAR)); \
			if( X->value==NULL ) \
			{ \
				LIS_SETERR_MEM((A->np+A->pad)*sizeof(LIS_SCALAR)); \
				return LIS_ERR_OUT_OF_MEMORY; \
			} \
			X->np  = A->np; \
			X->pad = A->pad; \
			x = X->value; \
		}

#define LIS_MATVEC_REDUCE0 \
		if( (A->np+A->pad)>(Y->np+Y->pad) ) \
		{ \
			Y->value = (LIS_SCALAR *)lis_realloc(Y->value,(A->np+A->pad)*sizeof(LIS_SCALAR)); \
			if( Y->value==NULL ) \
			{ \
				LIS_SETERR_MEM((A->np+A->pad)*sizeof(LIS_SCALAR)); \
				return LIS_ERR_OUT_OF_MEMORY; \
			} \
			Y->np  = A->np; \
			Y->pad = A->pad; \
			y      = Y->value; \
		}

#define LIS_MATVEC_REDUCE lis_reduce(A->commtable,y)

typedef void (*LIS_MATVEC_XXX)(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
typedef void (*LIS_MATVEC_XXXP)(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y);
typedef int (*LIS_MATVEC_FUNC)(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y);

#ifdef __cplusplus
extern "C"
{
#endif
	extern LIS_MATVEC_FUNC LIS_MATVEC;
	extern LIS_MATVEC_FUNC LIS_MATVECT;
	/*******************/
	/* CRS             */
	/*******************/
	extern void lis_matvec_crs(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvect_crs(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_crs_mp(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y);
	extern void lis_matvec_crs_mp2(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y);
	extern void lis_matvect_crs_mp(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y);
	extern void lis_matvect_crs_mp2(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y);
	/*******************/
	/* CCS             */
	/*******************/
	extern void lis_matvec_ccs(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvect_ccs(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvect_ccs_mp(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y);
	extern void lis_matvect_ccs_mp2(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y);
	extern void lis_matvec_ccs_mp(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y);
	extern void lis_matvec_ccs_mp2(LIS_MATRIX A, LIS_VECTOR X, LIS_VECTOR Y);
	/*******************/
	/* MSR             */
	/*******************/
	extern void lis_matvec_msr(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvect_msr(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	/*******************/
	/* ELL             */
	/*******************/
	extern void lis_matvec_ell(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvect_ell(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	/*******************/
	/* DIA             */
	/*******************/
	extern void lis_matvec_dia(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvect_dia(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	/*******************/
	/* JDS             */
	/*******************/
	extern void lis_matvec_jds(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
		extern void lis_matvec_jds_u7_1(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvect_jds(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	/*******************/
	/* DNS             */
	/*******************/
	extern void lis_matvec_dns(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvect_dns(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	/*******************/
	/* COO             */
	/*******************/
	extern void lis_matvec_coo(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvect_coo(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	/*******************/
	/* VBR             */
	/*******************/
	extern void lis_matvec_vbr(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvect_vbr(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	/*******************/
	/* BSR             */
	/*******************/
	extern LIS_MATVEC_XXX lis_matvec_bsr_xxx[4][4];
	extern void lis_matvec_bsr(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsr_1x1(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsr_1x2(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsr_1x3(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsr_1x4(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsr_2x1(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsr_2x2(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsr_2x3(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsr_2x4(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsr_3x1(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsr_3x2(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsr_3x3(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsr_3x4(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsr_4x1(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsr_4x2(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsr_4x3(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsr_4x4(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvect_bsr(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	/*******************/
	/* BSC             */
	/*******************/
	extern LIS_MATVEC_XXX lis_matvec_bsc_xxx[4][4];
	extern void lis_matvec_bsc(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsc_1x1(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsc_1x2(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsc_1x3(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsc_1x4(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsc_2x1(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsc_2x2(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsc_2x3(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsc_2x4(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsc_3x1(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsc_3x2(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsc_3x3(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsc_3x4(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsc_4x1(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsc_4x2(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsc_4x3(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvec_bsc_4x4(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);
	extern void lis_matvect_bsc(LIS_MATRIX A, LIS_SCALAR x[], LIS_SCALAR y[]);

#ifdef __cplusplus
}
#endif

#endif
