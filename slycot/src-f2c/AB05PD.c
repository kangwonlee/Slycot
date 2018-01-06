#line 1 "AB05PD.f"
/* AB05PD.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

#line 1 "AB05PD.f"
/* Table of constant values */

static doublereal c_b9 = 0.;
static integer c__0 = 0;
static doublereal c_b26 = 1.;
static integer c__1 = 1;

/* Subroutine */ int ab05pd_(char *over, integer *n1, integer *m, integer *p, 
	integer *n2, doublereal *alpha, doublereal *a1, integer *lda1, 
	doublereal *b1, integer *ldb1, doublereal *c1, integer *ldc1, 
	doublereal *d1, integer *ldd1, doublereal *a2, integer *lda2, 
	doublereal *b2, integer *ldb2, doublereal *c2, integer *ldc2, 
	doublereal *d2, integer *ldd2, integer *n, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, integer *info, ftnlen over_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, a1_dim1, a1_offset, a2_dim1, a2_offset, b_dim1, 
	    b_offset, b1_dim1, b1_offset, b2_dim1, b2_offset, c_dim1, 
	    c_offset, c1_dim1, c1_offset, c2_dim1, c2_offset, d_dim1, 
	    d_offset, d1_dim1, d1_offset, d2_dim1, d2_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, n1p1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical lover;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dlacpy_(char *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);


/*     SLICOT RELEASE 5.0. */

/*     Copyright (c) 2002-2009 NICONET e.V. */

/*     This program is free software: you can redistribute it and/or */
/*     modify it under the terms of the GNU General Public License as */
/*     published by the Free Software Foundation, either version 2 of */
/*     the License, or (at your option) any later version. */

/*     This program is distributed in the hope that it will be useful, */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/*     GNU General Public License for more details. */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program.  If not, see */
/*     <http://www.gnu.org/licenses/>. */

/*     PURPOSE */

/*     To compute the state-space model G = (A,B,C,D) corresponding to */
/*     the sum G = G1 + alpha*G2, where G1 = (A1,B1,C1,D1) and */
/*     G2 = (A2,B2,C2,D2).  G, G1, and G2 are the transfer-function */
/*     matrices of the corresponding state-space models. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     OVER    CHARACTER*1 */
/*             Indicates whether the user wishes to overlap pairs of */
/*             arrays, as follows: */
/*             = 'N':  Do not overlap; */
/*             = 'O':  Overlap pairs of arrays: A1 and A, B1 and B, */
/*                     C1 and C, and D1 and D, i.e. the same name is */
/*                     effectively used for each pair (for all pairs) */
/*                     in the routine call.  In this case, setting */
/*                     LDA1 = LDA, LDB1 = LDB, LDC1 = LDC, and LDD1 = LDD */
/*                     will give maximum efficiency. */

/*     Input/Output Parameters */

/*     N1      (input) INTEGER */
/*             The number of state variables in the first system, i.e. */
/*             the order of the matrix A1, the number of rows of B1 and */
/*             the number of columns of C1.  N1 >= 0. */

/*     M       (input) INTEGER */
/*             The number of input variables of the two systems, i.e. the */
/*             number of columns of matrices B1, D1, B2 and D2.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of output variables of the two systems, i.e. */
/*             the number of rows of matrices C1, D1, C2 and D2.  P >= 0. */

/*     N2      (input) INTEGER */
/*             The number of state variables in the second system, i.e. */
/*             the order of the matrix A2, the number of rows of B2 and */
/*             the number of columns of C2.  N2 >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The coefficient multiplying G2. */

/*     A1      (input) DOUBLE PRECISION array, dimension (LDA1,N1) */
/*             The leading N1-by-N1 part of this array must contain the */
/*             state transition matrix A1 for the first system. */

/*     LDA1    INTEGER */
/*             The leading dimension of array A1.  LDA1 >= MAX(1,N1). */

/*     B1      (input) DOUBLE PRECISION array, dimension (LDB1,M) */
/*             The leading N1-by-M part of this array must contain the */
/*             input/state matrix B1 for the first system. */

/*     LDB1    INTEGER */
/*             The leading dimension of array B1.  LDB1 >= MAX(1,N1). */

/*     C1      (input) DOUBLE PRECISION array, dimension (LDC1,N1) */
/*             The leading P-by-N1 part of this array must contain the */
/*             state/output matrix C1 for the first system. */

/*     LDC1    INTEGER */
/*             The leading dimension of array C1. */
/*             LDC1 >= MAX(1,P) if N1 > 0. */
/*             LDC1 >= 1 if N1 = 0. */

/*     D1      (input) DOUBLE PRECISION array, dimension (LDD1,M) */
/*             The leading P-by-M part of this array must contain the */
/*             input/output matrix D1 for the first system. */

/*     LDD1    INTEGER */
/*             The leading dimension of array D1.  LDD1 >= MAX(1,P). */

/*     A2      (input) DOUBLE PRECISION array, dimension (LDA2,N2) */
/*             The leading N2-by-N2 part of this array must contain the */
/*             state transition matrix A2 for the second system. */

/*     LDA2    INTEGER */
/*             The leading dimension of array A2.  LDA2 >= MAX(1,N2). */

/*     B2      (input) DOUBLE PRECISION array, dimension (LDB2,M) */
/*             The leading N2-by-M part of this array must contain the */
/*             input/state matrix B2 for the second system. */

/*     LDB2    INTEGER */
/*             The leading dimension of array B2.  LDB2 >= MAX(1,N2). */

/*     C2      (input) DOUBLE PRECISION array, dimension (LDC2,N2) */
/*             The leading P-by-N2 part of this array must contain the */
/*             state/output matrix C2 for the second system. */

/*     LDC2    INTEGER */
/*             The leading dimension of array C2. */
/*             LDC2 >= MAX(1,P) if N2 > 0. */
/*             LDC2 >= 1 if N2 = 0. */

/*     D2      (input) DOUBLE PRECISION array, dimension (LDD2,M) */
/*             The leading P-by-M part of this array must contain the */
/*             input/output matrix D2 for the second system. */

/*     LDD2    INTEGER */
/*             The leading dimension of array D2.  LDD2 >= MAX(1,P). */

/*     N       (output) INTEGER */
/*             The number of state variables (N1 + N2) in the resulting */
/*             system, i.e. the order of the matrix A, the number of rows */
/*             of B and the number of columns of C. */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N1+N2) */
/*             The leading N-by-N part of this array contains the state */
/*             transition matrix A for the resulting system. */
/*             The array A can overlap A1 if OVER = 'O'. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N1+N2). */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array contains the */
/*             input/state matrix B for the resulting system. */
/*             The array B can overlap B1 if OVER = 'O'. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N1+N2). */

/*     C       (output) DOUBLE PRECISION array, dimension (LDC,N1+N2) */
/*             The leading P-by-N part of this array contains the */
/*             state/output matrix C for the resulting system. */
/*             The array C can overlap C1 if OVER = 'O'. */

/*     LDC     INTEGER */
/*             The leading dimension of array C. */
/*             LDC >= MAX(1,P) if N1+N2 > 0. */
/*             LDC >= 1 if N1+N2 = 0. */

/*     D       (output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading P-by-M part of this array contains the */
/*             input/output matrix D for the resulting system. */
/*             The array D can overlap D1 if OVER = 'O'. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The matrices of the resulting systems are determined as: */

/*           ( A1   0  )             ( B1 ) */
/*       A = (         ) ,       B = (    ) , */
/*           ( 0    A2 )             ( B2 ) */

/*       C = ( C1  alpha*C2 ) ,  D = D1 + alpha*D2 . */

/*     REFERENCES */

/*     None */

/*     NUMERICAL ASPECTS */

/*     None */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Research Establishment, */
/*     Oberpfaffenhofen, Germany, and V. Sima, Katholieke Univ. Leuven, */
/*     Belgium, Nov. 1996. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, July 2003, */
/*     Feb. 2004. */

/*     KEYWORDS */

/*     Multivariable system, state-space model, state-space */
/*     representation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 235 "AB05PD.f"
    /* Parameter adjustments */
#line 235 "AB05PD.f"
    a1_dim1 = *lda1;
#line 235 "AB05PD.f"
    a1_offset = 1 + a1_dim1;
#line 235 "AB05PD.f"
    a1 -= a1_offset;
#line 235 "AB05PD.f"
    b1_dim1 = *ldb1;
#line 235 "AB05PD.f"
    b1_offset = 1 + b1_dim1;
#line 235 "AB05PD.f"
    b1 -= b1_offset;
#line 235 "AB05PD.f"
    c1_dim1 = *ldc1;
#line 235 "AB05PD.f"
    c1_offset = 1 + c1_dim1;
#line 235 "AB05PD.f"
    c1 -= c1_offset;
#line 235 "AB05PD.f"
    d1_dim1 = *ldd1;
#line 235 "AB05PD.f"
    d1_offset = 1 + d1_dim1;
#line 235 "AB05PD.f"
    d1 -= d1_offset;
#line 235 "AB05PD.f"
    a2_dim1 = *lda2;
#line 235 "AB05PD.f"
    a2_offset = 1 + a2_dim1;
#line 235 "AB05PD.f"
    a2 -= a2_offset;
#line 235 "AB05PD.f"
    b2_dim1 = *ldb2;
#line 235 "AB05PD.f"
    b2_offset = 1 + b2_dim1;
#line 235 "AB05PD.f"
    b2 -= b2_offset;
#line 235 "AB05PD.f"
    c2_dim1 = *ldc2;
#line 235 "AB05PD.f"
    c2_offset = 1 + c2_dim1;
#line 235 "AB05PD.f"
    c2 -= c2_offset;
#line 235 "AB05PD.f"
    d2_dim1 = *ldd2;
#line 235 "AB05PD.f"
    d2_offset = 1 + d2_dim1;
#line 235 "AB05PD.f"
    d2 -= d2_offset;
#line 235 "AB05PD.f"
    a_dim1 = *lda;
#line 235 "AB05PD.f"
    a_offset = 1 + a_dim1;
#line 235 "AB05PD.f"
    a -= a_offset;
#line 235 "AB05PD.f"
    b_dim1 = *ldb;
#line 235 "AB05PD.f"
    b_offset = 1 + b_dim1;
#line 235 "AB05PD.f"
    b -= b_offset;
#line 235 "AB05PD.f"
    c_dim1 = *ldc;
#line 235 "AB05PD.f"
    c_offset = 1 + c_dim1;
#line 235 "AB05PD.f"
    c__ -= c_offset;
#line 235 "AB05PD.f"
    d_dim1 = *ldd;
#line 235 "AB05PD.f"
    d_offset = 1 + d_dim1;
#line 235 "AB05PD.f"
    d__ -= d_offset;
#line 235 "AB05PD.f"

#line 235 "AB05PD.f"
    /* Function Body */
#line 235 "AB05PD.f"
    lover = lsame_(over, "O", (ftnlen)1, (ftnlen)1);
#line 236 "AB05PD.f"
    *n = *n1 + *n2;
#line 237 "AB05PD.f"
    *info = 0;

/*     Test the input scalar arguments. */

#line 241 "AB05PD.f"
    if (! lover && ! lsame_(over, "N", (ftnlen)1, (ftnlen)1)) {
#line 242 "AB05PD.f"
	*info = -1;
#line 243 "AB05PD.f"
    } else if (*n1 < 0) {
#line 244 "AB05PD.f"
	*info = -2;
#line 245 "AB05PD.f"
    } else if (*m < 0) {
#line 246 "AB05PD.f"
	*info = -3;
#line 247 "AB05PD.f"
    } else if (*p < 0) {
#line 248 "AB05PD.f"
	*info = -4;
#line 249 "AB05PD.f"
    } else if (*n2 < 0) {
#line 250 "AB05PD.f"
	*info = -5;
#line 251 "AB05PD.f"
    } else if (*lda1 < max(1,*n1)) {
#line 252 "AB05PD.f"
	*info = -8;
#line 253 "AB05PD.f"
    } else if (*ldb1 < max(1,*n1)) {
#line 254 "AB05PD.f"
	*info = -10;
#line 255 "AB05PD.f"
    } else if (*n1 > 0 && *ldc1 < max(1,*p) || *n1 == 0 && *ldc1 < 1) {
#line 257 "AB05PD.f"
	*info = -12;
#line 258 "AB05PD.f"
    } else if (*ldd1 < max(1,*p)) {
#line 259 "AB05PD.f"
	*info = -14;
#line 260 "AB05PD.f"
    } else if (*lda2 < max(1,*n2)) {
#line 261 "AB05PD.f"
	*info = -16;
#line 262 "AB05PD.f"
    } else if (*ldb2 < max(1,*n2)) {
#line 263 "AB05PD.f"
	*info = -18;
#line 264 "AB05PD.f"
    } else if (*n2 > 0 && *ldc2 < max(1,*p) || *n2 == 0 && *ldc2 < 1) {
#line 266 "AB05PD.f"
	*info = -20;
#line 267 "AB05PD.f"
    } else if (*ldd2 < max(1,*p)) {
#line 268 "AB05PD.f"
	*info = -22;
#line 269 "AB05PD.f"
    } else if (*lda < max(1,*n)) {
#line 270 "AB05PD.f"
	*info = -25;
#line 271 "AB05PD.f"
    } else if (*ldb < max(1,*n)) {
#line 272 "AB05PD.f"
	*info = -27;
#line 273 "AB05PD.f"
    } else if (*n > 0 && *ldc < max(1,*p) || *n == 0 && *ldc < 1) {
#line 275 "AB05PD.f"
	*info = -29;
#line 276 "AB05PD.f"
    } else if (*ldd < max(1,*p)) {
#line 277 "AB05PD.f"
	*info = -31;
#line 278 "AB05PD.f"
    }

#line 280 "AB05PD.f"
    if (*info != 0) {

/*        Error return. */

#line 284 "AB05PD.f"
	i__1 = -(*info);
#line 284 "AB05PD.f"
	xerbla_("AB05PD", &i__1, (ftnlen)6);
#line 285 "AB05PD.f"
	return 0;
#line 286 "AB05PD.f"
    }

/*     Quick return if possible. */

/* Computing MAX */
#line 290 "AB05PD.f"
    i__1 = *n, i__2 = min(*m,*p);
#line 290 "AB05PD.f"
    if (max(i__1,i__2) == 0) {
#line 290 "AB05PD.f"
	return 0;
#line 290 "AB05PD.f"
    }

#line 293 "AB05PD.f"
    n1p1 = *n1 + 1;

/*                       ( A1   0  ) */
/*     Construct     A = (         ) . */
/*                       ( 0    A2 ) */

#line 299 "AB05PD.f"
    if (lover && *lda1 <= *lda) {
#line 300 "AB05PD.f"
	if (*lda1 < *lda) {

#line 302 "AB05PD.f"
	    for (j = *n1; j >= 1; --j) {
#line 303 "AB05PD.f"
		for (i__ = *n1; i__ >= 1; --i__) {
#line 304 "AB05PD.f"
		    a[i__ + j * a_dim1] = a1[i__ + j * a1_dim1];
#line 305 "AB05PD.f"
/* L10: */
#line 305 "AB05PD.f"
		}
#line 306 "AB05PD.f"
/* L20: */
#line 306 "AB05PD.f"
	    }

#line 308 "AB05PD.f"
	}
#line 309 "AB05PD.f"
    } else {
#line 310 "AB05PD.f"
	dlacpy_("F", n1, n1, &a1[a1_offset], lda1, &a[a_offset], lda, (ftnlen)
		1);
#line 311 "AB05PD.f"
    }

#line 313 "AB05PD.f"
    if (*n2 > 0) {
#line 314 "AB05PD.f"
	dlaset_("F", n1, n2, &c_b9, &c_b9, &a[n1p1 * a_dim1 + 1], lda, (
		ftnlen)1);
#line 315 "AB05PD.f"
	dlaset_("F", n2, n1, &c_b9, &c_b9, &a[n1p1 + a_dim1], lda, (ftnlen)1);
#line 316 "AB05PD.f"
	dlacpy_("F", n2, n2, &a2[a2_offset], lda2, &a[n1p1 + n1p1 * a_dim1], 
		lda, (ftnlen)1);
#line 317 "AB05PD.f"
    }

/*                        ( B1 ) */
/*     Construct      B = (    ) . */
/*                        ( B2 ) */

#line 323 "AB05PD.f"
    if (lover && *ldb1 <= *ldb) {
#line 324 "AB05PD.f"
	if (*ldb1 < *ldb) {

#line 326 "AB05PD.f"
	    for (j = *m; j >= 1; --j) {
#line 327 "AB05PD.f"
		for (i__ = *n1; i__ >= 1; --i__) {
#line 328 "AB05PD.f"
		    b[i__ + j * b_dim1] = b1[i__ + j * b1_dim1];
#line 329 "AB05PD.f"
/* L30: */
#line 329 "AB05PD.f"
		}
#line 330 "AB05PD.f"
/* L40: */
#line 330 "AB05PD.f"
	    }

#line 332 "AB05PD.f"
	}
#line 333 "AB05PD.f"
    } else {
#line 334 "AB05PD.f"
	dlacpy_("F", n1, m, &b1[b1_offset], ldb1, &b[b_offset], ldb, (ftnlen)
		1);
#line 335 "AB05PD.f"
    }

#line 337 "AB05PD.f"
    if (*n2 > 0) {
#line 337 "AB05PD.f"
	dlacpy_("F", n2, m, &b2[b2_offset], ldb2, &b[n1p1 + b_dim1], ldb, (
		ftnlen)1);
#line 337 "AB05PD.f"
    }

/*     Construct      C = ( C1 alpha*C2 ) . */

#line 342 "AB05PD.f"
    if (lover && *ldc1 <= *ldc) {
#line 343 "AB05PD.f"
	if (*ldc1 < *ldc) {

#line 345 "AB05PD.f"
	    for (j = *n1; j >= 1; --j) {
#line 346 "AB05PD.f"
		for (i__ = *p; i__ >= 1; --i__) {
#line 347 "AB05PD.f"
		    c__[i__ + j * c_dim1] = c1[i__ + j * c1_dim1];
#line 348 "AB05PD.f"
/* L50: */
#line 348 "AB05PD.f"
		}
#line 349 "AB05PD.f"
/* L60: */
#line 349 "AB05PD.f"
	    }

#line 351 "AB05PD.f"
	}
#line 352 "AB05PD.f"
    } else {
#line 353 "AB05PD.f"
	dlacpy_("F", p, n1, &c1[c1_offset], ldc1, &c__[c_offset], ldc, (
		ftnlen)1);
#line 354 "AB05PD.f"
    }

#line 356 "AB05PD.f"
    if (*n2 > 0) {
#line 357 "AB05PD.f"
	dlacpy_("F", p, n2, &c2[c2_offset], ldc2, &c__[n1p1 * c_dim1 + 1], 
		ldc, (ftnlen)1);
#line 358 "AB05PD.f"
	if (*alpha != 1.) {
#line 358 "AB05PD.f"
	    dlascl_("G", &c__0, &c__0, &c_b26, alpha, p, n2, &c__[n1p1 * 
		    c_dim1 + 1], ldc, info, (ftnlen)1);
#line 358 "AB05PD.f"
	}
#line 361 "AB05PD.f"
    }

/*     Construct       D = D1 + alpha*D2 . */

#line 365 "AB05PD.f"
    if (lover && *ldd1 <= *ldd) {
#line 366 "AB05PD.f"
	if (*ldd1 < *ldd) {

#line 368 "AB05PD.f"
	    for (j = *m; j >= 1; --j) {
#line 369 "AB05PD.f"
		for (i__ = *p; i__ >= 1; --i__) {
#line 370 "AB05PD.f"
		    d__[i__ + j * d_dim1] = d1[i__ + j * d1_dim1];
#line 371 "AB05PD.f"
/* L70: */
#line 371 "AB05PD.f"
		}
#line 372 "AB05PD.f"
/* L80: */
#line 372 "AB05PD.f"
	    }

#line 374 "AB05PD.f"
	}
#line 375 "AB05PD.f"
    } else {
#line 376 "AB05PD.f"
	dlacpy_("F", p, m, &d1[d1_offset], ldd1, &d__[d_offset], ldd, (ftnlen)
		1);
#line 377 "AB05PD.f"
    }

#line 379 "AB05PD.f"
    i__1 = *m;
#line 379 "AB05PD.f"
    for (j = 1; j <= i__1; ++j) {
#line 380 "AB05PD.f"
	daxpy_(p, alpha, &d2[j * d2_dim1 + 1], &c__1, &d__[j * d_dim1 + 1], &
		c__1);
#line 381 "AB05PD.f"
/* L90: */
#line 381 "AB05PD.f"
    }

#line 383 "AB05PD.f"
    return 0;
/* *** Last line of AB05PD *** */
} /* ab05pd_ */

