/*-
 * Copyright (c) 2010 Nathan Lay
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR(S) ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef SVD3_H
#define SVD3_H

#include <stdio.h>
#include <math.h>

/* ensure that gcc will recognise restrict even when not compining with -std=c99 */
#ifdef __GNUC__
#define restrict __restrict__
#endif

/* 
 * NOTE: All computations are done COLUMN MAJOR!
 */

#if defined(MATLAB_FMT)
#define printmat3_fmt(A)	#A " = [ %g %g %g; %g %g %g; %g %g %g ]\n"
#define readmat3_fmt(A)		"%lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define printcol3_fmt(A,j)	#A " = [ %g %g %g ]'\n"
#define printcolp3_fmt(A,P,j)	#P #A " = [ %g %g %g ]'\n"
#define printrow3_fmt(A,i)	#A " = [ %g %g %g ]\n"
#define printL3_fmt(LDU,P)	"L = [ 1 0 0; %g 1 0; %g %g 1 ]\n"
#define printD3_fmt(LDU,P)	"D = [ %g 0 0; 0 %g 0; 0 0 %g ]\n"
#define printU3_fmt(LDU,P)	"U = [ 1 %g %g; 0 1 %g; 0 0 1 ]\n"
#define printP3_fmt(P)		#P " = [ %d %d %d; %d %d %d; %d %d %d ]\n"
#define printmatp3_fmt(A,P)	#A #P " = [ %g %g %g; %g %g %g; %g %g %g ]\n"
#else /* MATLAB_FMT */
#define printmat3_fmt(A)	#A " =\n%g %g %g\n%g %g %g\n%g %g %g\n"
#define readmat3_fmt(A)		"%lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define printcol3_fmt(A,j)	#A " =\n%g\n%g\n%g\n"
#define printcolp3_fmt(A,P,j)	#P #A " =\n%g\n%g\n%g\n"
#define printrow3_fmt(A,i)	#A " =\n%g %g %g\n"
#define printL3_fmt(LDU,P)	"L =\n1 0 0\n%g 1 0\n%g %g 1\n"
#define printD3_fmt(LDU,P)	"D =\n%g 0 0\n0 %g 0\n0 0 %g\n"
#define printU3_fmt(LDU,P)	"U =\n1 %g %g\n0 1 %g\n0 0 1\n"
#define printP3_fmt(P)		#P " =\n%d %d %d\n%d %d %d\n%d %d %d\n"
#define printmatp3_fmt(A,P)	#A #P " =\n%g %g %g\n%g %g %g\n%g %g %g\n"
#endif /* !MATLAB_FMT */

#define printmat3(A)	( printf(printmat3_fmt(A), 					\
	((double *)(A))[3*0+0], ((double *)(A))[3*1+0], ((double *)(A))[3*2+0],		\
	((double *)(A))[3*0+1], ((double *)(A))[3*1+1], ((double *)(A))[3*2+1], 	\
	((double *)(A))[3*0+2], ((double *)(A))[3*1+2], ((double *)(A))[3*2+2]) )

#define readmat3(A)	( scanf(readmat3_fmt(A),				\
	(double *)(A)+(3*0+0), (double *)(A)+(3*1+0), (double *)(A)+(3*2+0),	\
	(double *)(A)+(3*0+1), (double *)(A)+(3*1+1), (double *)(A)+(3*2+1),	\
	(double *)(A)+(3*0+2), (double *)(A)+(3*1+2), (double *)(A)+(3*2+2)) )

#define printcol3(A,j)	( printf(printcol3_fmt(A,j),	\
	((double *)(A))[3*(j)+0], 			\
	((double *)(A))[3*(j)+1], 			\
	((double *)(A))[3*(j)+2]) )

#define printcolp3(A,P,j)	( printf(printcolp3_fmt(A,P,j),	\
	((double *)(A))[3*(j)+(P)[0]], 				\
	((double *)(A))[3*(j)+(P)[1]], 				\
	((double *)(A))[3*(j)+(P)[2]]) )

#define printrow3(A,i)	( printf(printrow3_fmt(A,i)					\
	((double *)(A))[3*0+(i)], ((double *)(A))[3*1+(i)], ((double *)(A))[3*2+(i)]) )

#define printL3(LDU,P)	( printf(printL3_fmt(LDU,P), 			\
	((double *)(LDU))[3*(P)[0]+1], 					\
	((double *)(LDU))[3*(P)[0]+2], ((double *)(LDU))[3*(P)[1]+2]) )

#define printD3(LDU,P)	( printf(printD3_fmt(LDU,P),	\
	((double *)(LDU))[3*(P)[0]+0], 			\
	((double *)(LDU))[3*(P)[1]+1], 			\
	((double *)(LDU))[3*(P)[2]+2]) )

#define printU3(LDU,P)	( printf(printU3_fmt(LDU,P),			\
	((double *)(LDU))[3*(P)[1]+0], ((double *)(LDU))[3*(P)[2]+0], 	\
	((double *)(LDU))[3*(P)[2]+1]) )

#define printP3(P)	( printf(printP3_fmt(P),	\
	(P)[0]==0, (P)[1]==0, (P)[2]==0,		\
	(P)[0]==1, (P)[1]==1, (P)[2]==1,		\
	(P)[0]==2, (P)[1]==2, (P)[2]==2) )

#define printmatp3(A,P)	( printf(printmatp3_fmt(A,P),						\
	((double *)(A))[3*(P)[0]+0], ((double *)(A))[3*(P)[1]+0], ((double *)(A))[3*(P)[2]+0],	\
	((double *)(A))[3*(P)[0]+1], ((double *)(A))[3*(P)[1]+1], ((double *)(A))[3*(P)[2]+1],	\
	((double *)(A))[3*(P)[0]+2], ((double *)(A))[3*(P)[1]+2], ((double *)(A))[3*(P)[2]+2]) )


/* Computes cross product of 3D vectors x, y and stores the result in z */
static inline void cross(double * restrict z, const double * restrict x, 
	const double * restrict y);

/* Sorts 3 elements */
static inline void sort3(double * restrict x);

/* Normalizes a 3D vector (with respect to L2) */
static inline void unit3(double * restrict x);

/*
 * Solves for the roots of a monic cubic polynomial with 3 coefficients 
 * ordered by degree that is assumed to have 3 real roots (D <= 0) 
 */
void solvecubic(double * restrict c);

/* Computes the LDUP decomposition in-place */
void ldu3(double * restrict A, int * restrict P);

/* Does the backward-solve step, or U*x = y */
static inline void ldubsolve3(double * restrict x, const double * restrict y, 
	const double * restrict LDU, const int * restrict P);

/* Explicitly computes the SVD of a 3x3 matrix */
void svd3(double * restrict U, double * restrict S, double * restrict V, 
	const double * restrict A);

/* Computes the matrix multiplication C = A*B */
static inline void matmul3(double * restrict C, const double * restrict A, 
	const double * restrict B);

/* Computes the matrix multiplication y = A*x */
static inline void matvec3(double * restrict y, const double * restrict A,
	const double * restrict x);

/* Computes the matrix multiplication AA = A^T*A */
static inline void ata3(double * restrict AA, const double * restrict A);

/* Computes the matrix multiplication AA = A*A^T */
static inline void aat3(double * restrict AA, const double * restrict A);

/* Computes the matrix transpose of A */
static inline void trans3(double * restrict A);

static inline void cross(double * restrict z, const double * restrict x, 
	const double * restrict y) {
	z[0] = x[1]*y[2]-x[2]*y[1];
	z[1] = -(x[0]*y[2]-x[2]*y[0]);
	z[2] = x[0]*y[1]-x[1]*y[0];
}

static inline void sort3(double * restrict x) {
	double tmp;

	if (x[0] < x[1]) {
		tmp = x[0];
		x[0] = x[1];
		x[1] = tmp;
	}
	if (x[1] < x[2]) {
		if (x[0] < x[2]) {
			tmp = x[2];
			x[2] = x[1];
			x[1] = x[0];
			x[0] = tmp;
		}
		else {
			tmp = x[1];
			x[1] = x[2];
			x[2] = tmp;
		}
	}
}

static inline void unit3(double * restrict x) {
	double tmp = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	x[0] /= tmp;
	x[1] /= tmp;
	x[2] /= tmp;
}

static inline void ldubsolve3(double * restrict x, const double * restrict y, 
	const double * restrict LDU, const int * restrict P) {
	x[P[2]] = y[2];
	x[P[1]] = y[1] - LDU[3*P[2]+1]*x[P[2]];
	x[P[0]] = y[0] - LDU[3*P[2]+0]*x[P[2]] - LDU[3*P[1]+0]*x[P[1]];

#ifdef DEBUG
	puts("\nbsolve");
	printcol3(y,0);
	putchar('\n');
	printcolp3(x,P,0);
	putchar('\n');
#endif
}

static inline void matmul3(double * restrict C, const double * restrict A, 
	const double * restrict B) {
	C[3*0+0] = A[3*0+0]*B[3*0+0] + A[3*1+0]*B[3*0+1] + A[3*2+0]*B[3*0+2];
	C[3*1+0] = A[3*0+0]*B[3*1+0] + A[3*1+0]*B[3*1+1] + A[3*2+0]*B[3*1+2];
	C[3*2+0] = A[3*0+0]*B[3*2+0] + A[3*1+0]*B[3*2+1] + A[3*2+0]*B[3*2+2];

	C[3*0+1] = A[3*0+1]*B[3*0+0] + A[3*1+1]*B[3*0+1] + A[3*2+1]*B[3*0+2];
	C[3*1+1] = A[3*0+1]*B[3*1+0] + A[3*1+1]*B[3*1+1] + A[3*2+1]*B[3*1+2];
	C[3*2+1] = A[3*0+1]*B[3*2+0] + A[3*1+1]*B[3*2+1] + A[3*2+1]*B[3*2+2];

	C[3*0+2] = A[3*0+2]*B[3*0+0] + A[3*1+2]*B[3*0+1] + A[3*2+2]*B[3*0+2];
	C[3*1+2] = A[3*0+2]*B[3*1+0] + A[3*1+2]*B[3*1+1] + A[3*2+2]*B[3*1+2];
	C[3*2+2] = A[3*0+2]*B[3*2+0] + A[3*1+2]*B[3*2+1] + A[3*2+2]*B[3*2+2];
}

static inline void matvec3(double * restrict y, const double * restrict A,
	const double * restrict x) {
	y[0] = A[3*0+0]*x[0] + A[3*1+0]*x[1] + A[3*2+0]*x[2];
	y[1] = A[3*0+1]*x[0] + A[3*1+1]*x[1] + A[3*2+1]*x[2];
	y[2] = A[3*0+2]*x[0] + A[3*1+2]*x[1] + A[3*2+2]*x[2];
}

static inline void ata3(double * restrict AA, const double * restrict A) {
	AA[3*0+0] = A[3*0+0]*A[3*0+0] + A[3*0+1]*A[3*0+1] + A[3*0+2]*A[3*0+2];
	AA[3*1+0] = A[3*0+0]*A[3*1+0] + A[3*0+1]*A[3*1+1] + A[3*0+2]*A[3*1+2];
	AA[3*2+0] = A[3*0+0]*A[3*2+0] + A[3*0+1]*A[3*2+1] + A[3*0+2]*A[3*2+2];

	AA[3*0+1] = AA[3*1+0];
	AA[3*1+1] = A[3*1+0]*A[3*1+0] + A[3*1+1]*A[3*1+1] + A[3*1+2]*A[3*1+2];
	AA[3*2+1] = A[3*1+0]*A[3*2+0] + A[3*1+1]*A[3*2+1] + A[3*1+2]*A[3*2+2];

	AA[3*0+2] = AA[3*2+0];
	AA[3*1+2] = AA[3*2+1];
	AA[3*2+2] = A[3*2+0]*A[3*2+0] + A[3*2+1]*A[3*2+1] + A[3*2+2]*A[3*2+2];
}

static inline void aat3(double * restrict AA, const double * restrict A) {
	AA[3*0+0] = A[3*0+0]*A[3*0+0] + A[3*1+0]*A[3*1+0] + A[3*2+0]*A[3*2+0];
	AA[3*1+0] = A[3*0+0]*A[3*0+1] + A[3*1+0]*A[3*1+1] + A[3*2+0]*A[3*2+1];
	AA[3*2+0] = A[3*0+0]*A[3*0+2] + A[3*1+0]*A[3*1+2] + A[3*2+0]*A[3*2+2];

	AA[3*0+1] = AA[3*1+0];
	AA[3*1+1] = A[3*0+1]*A[3*0+1] + A[3*1+1]*A[3*1+1] + A[3*2+1]*A[3*2+1];
	AA[3*2+1] = A[3*0+1]*A[3*0+2] + A[3*1+1]*A[3*1+2] + A[3*2+1]*A[3*2+2];

	AA[3*0+2] = AA[3*2+0];
	AA[3*1+2] = AA[3*2+1];
	AA[3*2+2] = A[3*0+2]*A[3*0+2] + A[3*1+2]*A[3*1+2] + A[3*2+2]*A[3*2+2];
}

static inline void trans3(double * restrict A) {
	double tmp;

	tmp = A[3*1+0];
	A[3*1+0] = A[3*0+1];
	A[3*0+1] = tmp;

	tmp = A[3*2+0];
	A[3*2+0] = A[3*0+2];
	A[3*0+2] = tmp;

	tmp = A[3*2+1];
	A[3*2+1] = A[3*1+2];
	A[3*1+2] = tmp;
}

#endif /* !SVD3_H */

