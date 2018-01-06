#line 1 "MB02CV.f"
/* MB02CV.f -- translated by f2c (version 20100827).
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

#line 1 "MB02CV.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mb02cv_(char *typeg, char *strucg, integer *k, integer *
	n, integer *p, integer *q, integer *nb, integer *rnk, doublereal *a1, 
	integer *lda1, doublereal *a2, integer *lda2, doublereal *b, integer *
	ldb, doublereal *f1, integer *ldf1, doublereal *f2, integer *ldf2, 
	doublereal *g, integer *ldg, doublereal *cs, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen typeg_len, ftnlen strucg_len)
{
    /* System generated locals */
    integer a1_dim1, a1_offset, a2_dim1, a2_offset, b_dim1, b_offset, f1_dim1,
	     f1_offset, f2_dim1, f2_offset, g_dim1, g_offset, i__1, i__2, 
	    i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Local variables */
    static doublereal c__;
    static integer i__, j;
    static doublereal s;
    static integer ib, jj, nbl, len;
    static doublereal tau;
    static integer pos, col2, pst2;
    static doublereal beta;
    static logical lcol;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static logical ltri;
    static doublereal alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dlarf_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    ftnlen);
    static logical lrdef;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dlarfb_(char *, char *, char 
	    *, char *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, ftnlen, ftnlen, ftnlen, ftnlen), dlarft_(char *, char 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen), xerbla_(char *, integer 
	    *, ftnlen);
    static integer wrkmin;


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

/*     To apply the transformations created by the SLICOT Library routine */
/*     MB02CU on other columns / rows of the generator, contained in the */
/*     arrays F1, F2 and G. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYPEG   CHARACTER*1 */
/*             Specifies the type of the generator, as follows: */
/*             = 'D':  generator is column oriented and rank */
/*                     deficient; */
/*             = 'C':  generator is column oriented and not rank */
/*                     deficient; */
/*             = 'R':  generator is row oriented and not rank */
/*                     deficient. */
/*             Note that this parameter must be equivalent with the */
/*             used TYPEG in the call of MB02CU. */

/*     STRUCG  CHARACTER*1 */
/*             Information about the structure of the generators, */
/*             as follows: */
/*             = 'T':  the trailing block of the positive generator */
/*                     is upper / lower triangular, and the trailing */
/*                     block of the negative generator is zero; */
/*             = 'N':  no special structure to mention. */

/*     Input/Output Parameters */

/*     K       (input)  INTEGER */
/*             The number of rows in A1 to be processed.  K >= 0. */

/*     N       (input)  INTEGER */
/*             If TYPEG = 'D'  or  TYPEG = 'C', the number of rows in F1; */
/*             if TYPEG = 'R', the number of columns in F1.  N >= 0. */

/*     P       (input)  INTEGER */
/*             The number of columns of the positive generator.  P >= K. */

/*     Q       (input)  INTEGER */
/*             The number of columns in B. */
/*             If TYPEG = 'D',        Q >= K; */
/*             If TYPEG = 'C' or 'R', Q >= 0. */

/*     NB      (input)  INTEGER */
/*             On entry, if TYPEG = 'C'  or  TYPEG = 'R', NB specifies */
/*             the block size to be used in the blocked parts of the */
/*             algorithm. NB must be equivalent with the used block size */
/*             in the routine MB02CU. */

/*     RNK     (input)  INTEGER */
/*             If TYPEG = 'D', the number of linearly independent columns */
/*             in the generator as returned by MB02CU.  0 <= RNK <= K. */
/*             If TYPEG = 'C' or 'R', the value of this parameter is */
/*             irrelevant. */

/*     A1      (input)  DOUBLE PRECISION array, dimension */
/*             (LDA1, K) */
/*             On entry, if TYPEG = 'D', the leading K-by-K part of this */
/*             array must contain the matrix A1 as returned by MB02CU. */
/*             If TYPEG = 'C' or 'R', this array is not referenced. */

/*     LDA1    INTEGER */
/*             The leading dimension of the array A1. */
/*             If TYPEG = 'D',                   LDA1 >= MAX(1,K); */
/*             if TYPEG = 'C'  or  TYPEG = 'R',  LDA1 >= 1. */

/*     A2      (input)  DOUBLE PRECISION array, */
/*             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDA2, P-K); */
/*             if TYPEG = 'R',                   dimension (LDA2, K). */
/*             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             K-by-(P-K) part of this array must contain the matrix */
/*             A2 as returned by MB02CU. */
/*             On entry, if TYPEG = 'R', the leading (P-K)-by-K part of */
/*             this array must contain the matrix A2 as returned by */
/*             MB02CU. */

/*     LDA2    INTEGER */
/*             The leading dimension of the array A2. */
/*             If P = K,                  LDA2 >= 1; */
/*             If P > K and (TYPEG = 'D' or TYPEG = 'C'), */
/*                                        LDA2 >= MAX(1,K); */
/*             if P > K and TYPEG = 'R',  LDA2 >= P-K. */

/*     B       (input)  DOUBLE PRECISION array, */
/*             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDB, Q); */
/*             if TYPEG = 'R',                   dimension (LDB, K). */
/*             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             K-by-Q part of this array must contain the matrix B as */
/*             returned by MB02CU. */
/*             On entry, if TYPEG = 'R', the leading Q-by-K part of this */
/*             array must contain the matrix B as returned by MB02CU. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             If Q = 0,                  LDB >= 1; */
/*             If Q > 0 and (TYPEG = 'D' or TYPEG = 'C'), */
/*                                        LDB >= MAX(1,K); */
/*             if Q > 0 and TYPEG = 'R',  LDB >= Q. */

/*     F1      (input/output)  DOUBLE PRECISION array, */
/*             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDF1, K); */
/*             if TYPEG = 'R',                   dimension (LDF1, N). */
/*             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             N-by-K part of this array must contain the first part */
/*             of the positive generator to be processed. */
/*             On entry, if TYPEG = 'R', the leading K-by-N part of this */
/*             array must contain the first part of the positive */
/*             generator to be processed. */
/*             On exit, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             N-by-K part of this array contains the first part of the */
/*             transformed positive generator. */
/*             On exit, if TYPEG = 'R', the leading K-by-N part of this */
/*             array contains the first part of the transformed positive */
/*             generator. */

/*     LDF1    INTEGER */
/*             The leading dimension of the array F1. */
/*             If TYPEG = 'D'  or  TYPEG = 'C',   LDF1 >= MAX(1,N); */
/*             if TYPEG = 'R',                    LDF1 >= MAX(1,K). */

/*     F2      (input/output)  DOUBLE PRECISION array, */
/*             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDF2, P-K); */
/*             if TYPEG = 'R',                   dimension (LDF2, N). */
/*             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             N-by-(P-K) part of this array must contain the second part */
/*             of the positive generator to be processed. */
/*             On entry, if TYPEG = 'R', the leading (P-K)-by-N part of */
/*             this array must contain the second part of the positive */
/*             generator to be processed. */
/*             On exit, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             N-by-(P-K) part of this array contains the second part of */
/*             the transformed positive generator. */
/*             On exit, if TYPEG = 'R', the leading (P-K)-by-N part of */
/*             this array contains the second part of the transformed */
/*             positive generator. */

/*     LDF2    INTEGER */
/*             The leading dimension of the array F2. */
/*             If P = K,                  LDF2 >= 1; */
/*             If P > K and (TYPEG = 'D' or TYPEG = 'C'), */
/*                                        LDF2 >= MAX(1,N); */
/*             if P > K and TYPEG = 'R',  LDF2 >= P-K. */

/*     G       (input/output)  DOUBLE PRECISION array, */
/*             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDG, Q); */
/*             if TYPEG = 'R',                   dimension (LDG, N). */
/*             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             N-by-Q part of this array must contain the negative part */
/*             of the generator to be processed. */
/*             On entry, if TYPEG = 'R', the leading Q-by-N part of this */
/*             array must contain the negative part of the generator to */
/*             be processed. */
/*             On exit, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             N-by-Q part of this array contains the transformed */
/*             negative generator. */
/*             On exit, if TYPEG = 'R', the leading Q-by-N part of this */
/*             array contains the transformed negative generator. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G. */
/*             If Q = 0,                  LDG >= 1; */
/*             If Q > 0 and (TYPEG = 'D' or TYPEG = 'C'), */
/*                                        LDG >= MAX(1,N); */
/*             if Q > 0 and TYPEG = 'R',  LDG >= Q. */

/*     CS      (input)  DOUBLE PRECISION array, dimension (x) */
/*             If TYPEG = 'D' and P = K,                   x = 3*K; */
/*             If TYPEG = 'D' and P > K,                   x = 5*K; */
/*             If (TYPEG = 'C' or TYPEG = 'R') and P = K,  x = 4*K; */
/*             If (TYPEG = 'C' or TYPEG = 'R') and P > K,  x = 6*K. */
/*             On entry, the first x elements of this array must contain */
/*             Givens and modified hyperbolic rotation parameters, and */
/*             scalar factors of the Householder transformations as */
/*             returned by MB02CU. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = -23,  DWORK(1) returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             TYPEG = 'D':               LDWORK >= MAX(1,N); */
/*             (TYPEG = 'C' or TYPEG = 'R')  and  NB <= 0: */
/*                                        LDWORK >= MAX(1,N); */
/*             (TYPEG = 'C' or TYPEG = 'R')  and  NB >= 1: */
/*                                        LDWORK >= MAX(1,( N + K )*NB). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0(N*K*( P + Q )) floating point operations. */

/*     METHOD */

/*     The Householder transformations and modified hyperbolic rotations */
/*     computed by SLICOT Library routine MB02CU are applied to the */
/*     corresponding parts of the matrices F1, F2 and G. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, May 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, June 2001, */
/*     March 2004, March 2007. */

/*     KEYWORDS */

/*     Elementary matrix operations, Householder transformation, matrix */
/*     operations, Toeplitz matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

#line 276 "MB02CV.f"
    /* Parameter adjustments */
#line 276 "MB02CV.f"
    a1_dim1 = *lda1;
#line 276 "MB02CV.f"
    a1_offset = 1 + a1_dim1;
#line 276 "MB02CV.f"
    a1 -= a1_offset;
#line 276 "MB02CV.f"
    a2_dim1 = *lda2;
#line 276 "MB02CV.f"
    a2_offset = 1 + a2_dim1;
#line 276 "MB02CV.f"
    a2 -= a2_offset;
#line 276 "MB02CV.f"
    b_dim1 = *ldb;
#line 276 "MB02CV.f"
    b_offset = 1 + b_dim1;
#line 276 "MB02CV.f"
    b -= b_offset;
#line 276 "MB02CV.f"
    f1_dim1 = *ldf1;
#line 276 "MB02CV.f"
    f1_offset = 1 + f1_dim1;
#line 276 "MB02CV.f"
    f1 -= f1_offset;
#line 276 "MB02CV.f"
    f2_dim1 = *ldf2;
#line 276 "MB02CV.f"
    f2_offset = 1 + f2_dim1;
#line 276 "MB02CV.f"
    f2 -= f2_offset;
#line 276 "MB02CV.f"
    g_dim1 = *ldg;
#line 276 "MB02CV.f"
    g_offset = 1 + g_dim1;
#line 276 "MB02CV.f"
    g -= g_offset;
#line 276 "MB02CV.f"
    --cs;
#line 276 "MB02CV.f"
    --dwork;
#line 276 "MB02CV.f"

#line 276 "MB02CV.f"
    /* Function Body */
#line 276 "MB02CV.f"
    *info = 0;
/* Computing MAX */
#line 277 "MB02CV.f"
    i__1 = 0, i__2 = *p - *k;
#line 277 "MB02CV.f"
    col2 = max(i__1,i__2);
#line 278 "MB02CV.f"
    lrdef = lsame_(typeg, "D", (ftnlen)1, (ftnlen)1);
#line 279 "MB02CV.f"
    lcol = lsame_(typeg, "C", (ftnlen)1, (ftnlen)1);
#line 280 "MB02CV.f"
    ltri = lsame_(strucg, "T", (ftnlen)1, (ftnlen)1);
#line 281 "MB02CV.f"
    if (lrdef) {
#line 282 "MB02CV.f"
	wrkmin = max(1,*n);
#line 283 "MB02CV.f"
    } else {
#line 284 "MB02CV.f"
	if (*nb >= 1) {
/* Computing MAX */
#line 285 "MB02CV.f"
	    i__1 = 1, i__2 = (*n + *k) * *nb;
#line 285 "MB02CV.f"
	    wrkmin = max(i__1,i__2);
#line 286 "MB02CV.f"
	} else {
#line 287 "MB02CV.f"
	    wrkmin = max(1,*n);
#line 288 "MB02CV.f"
	}
#line 289 "MB02CV.f"
    }

/*     Check the scalar input parameters. */

#line 293 "MB02CV.f"
    if (! (lcol || lrdef || lsame_(typeg, "R", (ftnlen)1, (ftnlen)1))) {
#line 294 "MB02CV.f"
	*info = -1;
#line 295 "MB02CV.f"
    } else if (! (ltri || lsame_(strucg, "N", (ftnlen)1, (ftnlen)1))) {
#line 296 "MB02CV.f"
	*info = -2;
#line 297 "MB02CV.f"
    } else if (*k < 0) {
#line 298 "MB02CV.f"
	*info = -3;
#line 299 "MB02CV.f"
    } else if (*n < 0) {
#line 300 "MB02CV.f"
	*info = -4;
#line 301 "MB02CV.f"
    } else if (*p < *k) {
#line 302 "MB02CV.f"
	*info = -5;
#line 303 "MB02CV.f"
    } else if (*q < 0 || lrdef && *q < *k) {
#line 304 "MB02CV.f"
	*info = -6;
#line 305 "MB02CV.f"
    } else if (lrdef && (*rnk < 0 || *rnk > *k)) {
#line 306 "MB02CV.f"
	*info = -8;
#line 307 "MB02CV.f"
    } else if (*lda1 < 1 || lrdef && *lda1 < *k) {
#line 308 "MB02CV.f"
	*info = -10;
#line 309 "MB02CV.f"
    } else if (*p == *k && *lda2 < 1 || *p > *k && (lrdef || lcol) && *lda2 < 
	    max(1,*k) || *p > *k && ! (lrdef || lcol) && *lda2 < *p - *k) {
#line 314 "MB02CV.f"
	*info = -12;
#line 315 "MB02CV.f"
    } else if (*q == 0 && *ldb < 1 || *q > 0 && (lrdef || lcol) && *ldb < max(
	    1,*k) || *q > 0 && ! (lrdef || lcol) && *ldb < *q) {
#line 320 "MB02CV.f"
	*info = -14;
#line 321 "MB02CV.f"
    } else if ((lrdef || lcol) && *ldf1 < max(1,*n)) {
#line 322 "MB02CV.f"
	*info = -16;
#line 323 "MB02CV.f"
    } else if (! (lrdef || lcol) && *ldf1 < max(1,*k)) {
#line 325 "MB02CV.f"
	*info = -16;
#line 326 "MB02CV.f"
    } else if (*p == *k && *ldf2 < 1 || *p > *k && (lrdef || lcol) && *ldf2 < 
	    max(1,*n) || *p > *k && ! (lrdef || lcol) && *ldf2 < *p - *k) {
#line 331 "MB02CV.f"
	*info = -18;
#line 332 "MB02CV.f"
    } else if (*q == 0 && *ldg < 1 || *q > 0 && (lrdef || lcol) && *ldg < max(
	    1,*n) || *q > 0 && ! (lrdef || lcol) && *ldg < *q) {
#line 337 "MB02CV.f"
	*info = -20;
#line 338 "MB02CV.f"
    } else if (*ldwork < wrkmin) {
#line 339 "MB02CV.f"
	dwork[1] = (doublereal) wrkmin;
#line 340 "MB02CV.f"
	*info = -23;
#line 341 "MB02CV.f"
    }

/*     Return if there were illegal values. */

#line 345 "MB02CV.f"
    if (*info != 0) {
#line 346 "MB02CV.f"
	i__1 = -(*info);
#line 346 "MB02CV.f"
	xerbla_("MB02CV", &i__1, (ftnlen)6);
#line 347 "MB02CV.f"
	return 0;
#line 348 "MB02CV.f"
    }

/*     Quick return if possible. */

#line 352 "MB02CV.f"
    if (min(*k,*n) == 0 || ! lrdef && *q == 0 && *p == *k) {
#line 354 "MB02CV.f"
	return 0;
#line 355 "MB02CV.f"
    }

#line 357 "MB02CV.f"
    if (lrdef) {

/*        Deficient generator. */

#line 361 "MB02CV.f"
	if (col2 == 0) {
#line 362 "MB02CV.f"
	    pst2 = *k << 1;
#line 363 "MB02CV.f"
	} else {
#line 364 "MB02CV.f"
	    pst2 = *k << 2;
#line 365 "MB02CV.f"
	}

#line 367 "MB02CV.f"
	i__1 = *rnk;
#line 367 "MB02CV.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Apply elementary reflectors. */

#line 371 "MB02CV.f"
	    if (col2 > 1) {
#line 372 "MB02CV.f"
		tau = a2[i__ + a2_dim1];
#line 373 "MB02CV.f"
		a2[i__ + a2_dim1] = 1.;
#line 374 "MB02CV.f"
		dlarf_("Right", n, &col2, &a2[i__ + a2_dim1], lda2, &tau, &f2[
			f2_offset], ldf2, &dwork[1], (ftnlen)5);
#line 376 "MB02CV.f"
		a2[i__ + a2_dim1] = tau;
#line 377 "MB02CV.f"
	    }

#line 379 "MB02CV.f"
	    if (*k > i__) {
#line 380 "MB02CV.f"
		alpha = a1[i__ + i__ * a1_dim1];
#line 381 "MB02CV.f"
		a1[i__ + i__ * a1_dim1] = 1.;
#line 382 "MB02CV.f"
		i__2 = *k - i__ + 1;
#line 382 "MB02CV.f"
		dlarf_("Right", n, &i__2, &a1[i__ + i__ * a1_dim1], lda1, &cs[
			pst2 + i__], &f1[i__ * f1_dim1 + 1], ldf1, &dwork[1], 
			(ftnlen)5);
#line 384 "MB02CV.f"
		a1[i__ + i__ * a1_dim1] = alpha;
#line 385 "MB02CV.f"
	    }

#line 387 "MB02CV.f"
	    if (col2 > 0) {
#line 388 "MB02CV.f"
		c__ = cs[(*k << 1) + (i__ << 1) - 1];
#line 389 "MB02CV.f"
		s = cs[(*k << 1) + (i__ << 1)];
#line 390 "MB02CV.f"
		drot_(n, &f1[i__ * f1_dim1 + 1], &c__1, &f2[f2_offset], &c__1,
			 &c__, &s);
#line 391 "MB02CV.f"
	    }

#line 393 "MB02CV.f"
	    if (*q > 1) {
#line 394 "MB02CV.f"
		tau = b[i__ + b_dim1];
#line 395 "MB02CV.f"
		b[i__ + b_dim1] = 1.;
#line 396 "MB02CV.f"
		dlarf_("Right", n, q, &b[i__ + b_dim1], ldb, &tau, &g[
			g_offset], ldg, &dwork[1], (ftnlen)5);
#line 398 "MB02CV.f"
		b[i__ + b_dim1] = tau;
#line 399 "MB02CV.f"
	    }

/*           Apply hyperbolic rotation. */

#line 403 "MB02CV.f"
	    c__ = cs[(i__ << 1) - 1];
#line 404 "MB02CV.f"
	    s = cs[i__ * 2];
#line 405 "MB02CV.f"
	    d__1 = 1. / c__;
#line 405 "MB02CV.f"
	    dscal_(n, &d__1, &f1[i__ * f1_dim1 + 1], &c__1);
#line 406 "MB02CV.f"
	    d__1 = -s / c__;
#line 406 "MB02CV.f"
	    daxpy_(n, &d__1, &g[g_dim1 + 1], &c__1, &f1[i__ * f1_dim1 + 1], &
		    c__1);
#line 407 "MB02CV.f"
	    dscal_(n, &c__, &g[g_dim1 + 1], &c__1);
#line 408 "MB02CV.f"
	    d__1 = -s;
#line 408 "MB02CV.f"
	    daxpy_(n, &d__1, &f1[i__ * f1_dim1 + 1], &c__1, &g[g_dim1 + 1], &
		    c__1);
#line 409 "MB02CV.f"
/* L10: */
#line 409 "MB02CV.f"
	}

#line 411 "MB02CV.f"
	len = *q;
#line 412 "MB02CV.f"
	pos = 1;

#line 414 "MB02CV.f"
	i__1 = *k;
#line 414 "MB02CV.f"
	for (j = *rnk + 1; j <= i__1; ++j) {

/*           Apply the reductions working on singular rows. */

#line 418 "MB02CV.f"
	    if (col2 > 1) {
#line 419 "MB02CV.f"
		tau = a2[j + a2_dim1];
#line 420 "MB02CV.f"
		a2[j + a2_dim1] = 1.;
#line 421 "MB02CV.f"
		dlarf_("Right", n, &col2, &a2[j + a2_dim1], lda2, &tau, &f2[
			f2_offset], ldf2, &dwork[1], (ftnlen)5);
#line 423 "MB02CV.f"
		a2[j + a2_dim1] = tau;
#line 424 "MB02CV.f"
	    }
#line 425 "MB02CV.f"
	    if (*k > j) {
#line 426 "MB02CV.f"
		alpha = a1[j + j * a1_dim1];
#line 427 "MB02CV.f"
		a1[j + j * a1_dim1] = 1.;
#line 428 "MB02CV.f"
		i__2 = *k - j + 1;
#line 428 "MB02CV.f"
		dlarf_("Right", n, &i__2, &a1[j + j * a1_dim1], lda1, &cs[
			pst2 + j], &f1[j * f1_dim1 + 1], ldf1, &dwork[1], (
			ftnlen)5);
#line 430 "MB02CV.f"
		a1[j + j * a1_dim1] = alpha;
#line 431 "MB02CV.f"
	    }
#line 432 "MB02CV.f"
	    if (col2 > 0) {
#line 433 "MB02CV.f"
		c__ = cs[(*k << 1) + (j << 1) - 1];
#line 434 "MB02CV.f"
		s = cs[(*k << 1) + (j << 1)];
#line 435 "MB02CV.f"
		drot_(n, &f1[j * f1_dim1 + 1], &c__1, &f2[f2_offset], &c__1, &
			c__, &s);
#line 436 "MB02CV.f"
	    }
#line 437 "MB02CV.f"
	    if (len > 1) {
#line 438 "MB02CV.f"
		beta = b[j + pos * b_dim1];
#line 439 "MB02CV.f"
		b[j + pos * b_dim1] = 1.;
#line 440 "MB02CV.f"
		dlarf_("Right", n, &len, &b[j + pos * b_dim1], ldb, &cs[(j << 
			1) - 1], &g[pos * g_dim1 + 1], ldg, &dwork[1], (
			ftnlen)5);
#line 442 "MB02CV.f"
		b[j + pos * b_dim1] = beta;
#line 443 "MB02CV.f"
	    }
#line 444 "MB02CV.f"
	    --len;
#line 445 "MB02CV.f"
	    ++pos;
#line 446 "MB02CV.f"
/* L20: */
#line 446 "MB02CV.f"
	}

#line 448 "MB02CV.f"
    } else if (lcol) {

/*        Column oriented and not deficient generator. */

/*        Apply an LQ like hyperbolic/orthogonal blocked decomposition. */

#line 454 "MB02CV.f"
	if (ltri) {
/* Computing MAX */
#line 455 "MB02CV.f"
	    i__1 = *n - *k;
#line 455 "MB02CV.f"
	    len = max(i__1,0);
#line 456 "MB02CV.f"
	} else {
#line 457 "MB02CV.f"
	    len = *n;
#line 458 "MB02CV.f"
	}
#line 459 "MB02CV.f"
	if (col2 > 0) {

#line 461 "MB02CV.f"
	    nbl = min(col2,*nb);
#line 462 "MB02CV.f"
	    if (nbl > 0) {

/*              Blocked version. */

#line 466 "MB02CV.f"
		i__1 = *k - nbl + 1;
#line 466 "MB02CV.f"
		i__2 = nbl;
#line 466 "MB02CV.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 467 "MB02CV.f"
		    i__3 = *k - i__ + 1;
#line 467 "MB02CV.f"
		    ib = min(i__3,nbl);
#line 468 "MB02CV.f"
		    i__3 = *n + *k;
#line 468 "MB02CV.f"
		    dlarft_("Forward", "Rowwise", &col2, &ib, &a2[i__ + 
			    a2_dim1], lda2, &cs[(*k << 2) + i__], &dwork[1], &
			    i__3, (ftnlen)7, (ftnlen)7);
#line 470 "MB02CV.f"
		    i__3 = *n + *k;
#line 470 "MB02CV.f"
		    i__4 = *n + *k;
#line 470 "MB02CV.f"
		    dlarfb_("Right", "No Transpose", "Forward", "Rowwise", &
			    len, &col2, &ib, &a2[i__ + a2_dim1], lda2, &dwork[
			    1], &i__3, &f2[f2_offset], ldf2, &dwork[ib + 1], &
			    i__4, (ftnlen)5, (ftnlen)12, (ftnlen)7, (ftnlen)7)
			    ;

#line 475 "MB02CV.f"
		    i__3 = i__ + ib - 1;
#line 475 "MB02CV.f"
		    for (j = i__; j <= i__3; ++j) {
#line 476 "MB02CV.f"
			tau = a2[j + a2_dim1];
#line 477 "MB02CV.f"
			a2[j + a2_dim1] = 1.;
/* Computing MIN */
#line 478 "MB02CV.f"
			i__5 = col2, i__6 = j - i__ + 1;
#line 478 "MB02CV.f"
			i__4 = min(i__5,i__6);
#line 478 "MB02CV.f"
			dlarf_("Right", &len, &i__4, &a2[j + a2_dim1], lda2, &
				tau, &f2[f2_offset], ldf2, &dwork[1], (ftnlen)
				5);
#line 480 "MB02CV.f"
			a2[j + a2_dim1] = tau;
#line 481 "MB02CV.f"
			c__ = cs[(*k << 1) + (j << 1) - 1];
#line 482 "MB02CV.f"
			s = cs[(*k << 1) + (j << 1)];
#line 483 "MB02CV.f"
			drot_(&len, &f1[j * f1_dim1 + 1], &c__1, &f2[
				f2_offset], &c__1, &c__, &s);
#line 484 "MB02CV.f"
			if (ltri) {
#line 485 "MB02CV.f"
			    ++len;
#line 486 "MB02CV.f"
			    temp = f1[len + j * f1_dim1];
#line 487 "MB02CV.f"
			    f1[len + j * f1_dim1] = c__ * temp;
#line 488 "MB02CV.f"
			    f2[len + f2_dim1] = -s * temp;

#line 490 "MB02CV.f"
			    i__4 = col2;
#line 490 "MB02CV.f"
			    for (jj = 2; jj <= i__4; ++jj) {
#line 491 "MB02CV.f"
				f2[len + jj * f2_dim1] = 0.;
#line 492 "MB02CV.f"
/* L30: */
#line 492 "MB02CV.f"
			    }

#line 494 "MB02CV.f"
			}
#line 495 "MB02CV.f"
/* L40: */
#line 495 "MB02CV.f"
		    }

#line 497 "MB02CV.f"
/* L50: */
#line 497 "MB02CV.f"
		}

#line 499 "MB02CV.f"
	    } else {
#line 500 "MB02CV.f"
		i__ = 1;
#line 501 "MB02CV.f"
	    }

/*           Unblocked version for the last or only block. */

#line 505 "MB02CV.f"
	    i__2 = *k;
#line 505 "MB02CV.f"
	    for (j = i__; j <= i__2; ++j) {
#line 506 "MB02CV.f"
		if (col2 > 1) {
#line 507 "MB02CV.f"
		    tau = a2[j + a2_dim1];
#line 508 "MB02CV.f"
		    a2[j + a2_dim1] = 1.;
#line 509 "MB02CV.f"
		    dlarf_("Right", &len, &col2, &a2[j + a2_dim1], lda2, &tau,
			     &f2[f2_offset], ldf2, &dwork[1], (ftnlen)5);
#line 511 "MB02CV.f"
		    a2[j + a2_dim1] = tau;
#line 512 "MB02CV.f"
		}

#line 514 "MB02CV.f"
		c__ = cs[(*k << 1) + (j << 1) - 1];
#line 515 "MB02CV.f"
		s = cs[(*k << 1) + (j << 1)];
#line 516 "MB02CV.f"
		drot_(&len, &f1[j * f1_dim1 + 1], &c__1, &f2[f2_offset], &
			c__1, &c__, &s);
#line 517 "MB02CV.f"
		if (ltri) {
#line 518 "MB02CV.f"
		    ++len;
#line 519 "MB02CV.f"
		    temp = f1[len + j * f1_dim1];
#line 520 "MB02CV.f"
		    f1[len + j * f1_dim1] = c__ * temp;
#line 521 "MB02CV.f"
		    f2[len + f2_dim1] = -s * temp;

#line 523 "MB02CV.f"
		    i__1 = col2;
#line 523 "MB02CV.f"
		    for (jj = 2; jj <= i__1; ++jj) {
#line 524 "MB02CV.f"
			f2[len + jj * f2_dim1] = 0.;
#line 525 "MB02CV.f"
/* L60: */
#line 525 "MB02CV.f"
		    }

#line 527 "MB02CV.f"
		}
#line 528 "MB02CV.f"
/* L70: */
#line 528 "MB02CV.f"
	    }

#line 530 "MB02CV.f"
	    pst2 = *k * 5;
#line 531 "MB02CV.f"
	} else {
#line 532 "MB02CV.f"
	    pst2 = *k << 1;
#line 533 "MB02CV.f"
	}

#line 535 "MB02CV.f"
	if (ltri) {
#line 536 "MB02CV.f"
	    len = *n - *k;
#line 537 "MB02CV.f"
	} else {
#line 538 "MB02CV.f"
	    len = *n;
#line 539 "MB02CV.f"
	}

#line 541 "MB02CV.f"
	nbl = min(*q,*nb);
#line 542 "MB02CV.f"
	if (nbl > 0) {

/*           Blocked version. */

#line 546 "MB02CV.f"
	    i__2 = *k - nbl + 1;
#line 546 "MB02CV.f"
	    i__1 = nbl;
#line 546 "MB02CV.f"
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 547 "MB02CV.f"
		i__3 = *k - i__ + 1;
#line 547 "MB02CV.f"
		ib = min(i__3,nbl);
#line 548 "MB02CV.f"
		i__3 = *n + *k;
#line 548 "MB02CV.f"
		dlarft_("Forward", "Rowwise", q, &ib, &b[i__ + b_dim1], ldb, &
			cs[pst2 + i__], &dwork[1], &i__3, (ftnlen)7, (ftnlen)
			7);
#line 550 "MB02CV.f"
		i__3 = *n + *k;
#line 550 "MB02CV.f"
		i__4 = *n + *k;
#line 550 "MB02CV.f"
		dlarfb_("Right", "NonTranspose", "Forward", "Rowwise", &len, 
			q, &ib, &b[i__ + b_dim1], ldb, &dwork[1], &i__3, &g[
			g_offset], ldg, &dwork[ib + 1], &i__4, (ftnlen)5, (
			ftnlen)12, (ftnlen)7, (ftnlen)7);

#line 555 "MB02CV.f"
		i__3 = i__ + ib - 1;
#line 555 "MB02CV.f"
		for (j = i__; j <= i__3; ++j) {
#line 556 "MB02CV.f"
		    tau = b[j + b_dim1];
#line 557 "MB02CV.f"
		    b[j + b_dim1] = 1.;
#line 558 "MB02CV.f"
		    i__4 = j - i__ + 1;
#line 558 "MB02CV.f"
		    dlarf_("Right", &len, &i__4, &b[j + b_dim1], ldb, &tau, &
			    g[g_offset], ldg, &dwork[1], (ftnlen)5);
#line 560 "MB02CV.f"
		    b[j + b_dim1] = tau;

/*                 Apply hyperbolic rotation. */

#line 564 "MB02CV.f"
		    c__ = cs[(j << 1) - 1];
#line 565 "MB02CV.f"
		    s = cs[j * 2];
#line 566 "MB02CV.f"
		    d__1 = 1. / c__;
#line 566 "MB02CV.f"
		    dscal_(&len, &d__1, &f1[j * f1_dim1 + 1], &c__1);
#line 567 "MB02CV.f"
		    d__1 = -s / c__;
#line 567 "MB02CV.f"
		    daxpy_(&len, &d__1, &g[g_offset], &c__1, &f1[j * f1_dim1 
			    + 1], &c__1);
#line 568 "MB02CV.f"
		    dscal_(&len, &c__, &g[g_offset], &c__1);
#line 569 "MB02CV.f"
		    d__1 = -s;
#line 569 "MB02CV.f"
		    daxpy_(&len, &d__1, &f1[j * f1_dim1 + 1], &c__1, &g[
			    g_offset], &c__1);
#line 570 "MB02CV.f"
		    if (ltri) {
#line 571 "MB02CV.f"
			++len;
#line 572 "MB02CV.f"
			g[len + g_dim1] = -s / c__ * f1[len + j * f1_dim1];
#line 573 "MB02CV.f"
			f1[len + j * f1_dim1] /= c__;

#line 575 "MB02CV.f"
			i__4 = *q;
#line 575 "MB02CV.f"
			for (jj = 2; jj <= i__4; ++jj) {
#line 576 "MB02CV.f"
			    g[len + jj * g_dim1] = 0.;
#line 577 "MB02CV.f"
/* L80: */
#line 577 "MB02CV.f"
			}

#line 579 "MB02CV.f"
		    }
#line 580 "MB02CV.f"
/* L90: */
#line 580 "MB02CV.f"
		}

#line 582 "MB02CV.f"
/* L100: */
#line 582 "MB02CV.f"
	    }

#line 584 "MB02CV.f"
	} else {
#line 585 "MB02CV.f"
	    i__ = 1;
#line 586 "MB02CV.f"
	}

/*        Unblocked version for the last or only block. */

#line 590 "MB02CV.f"
	i__1 = *k;
#line 590 "MB02CV.f"
	for (j = i__; j <= i__1; ++j) {
#line 591 "MB02CV.f"
	    if (*q > 1) {
#line 592 "MB02CV.f"
		tau = b[j + b_dim1];
#line 593 "MB02CV.f"
		b[j + b_dim1] = 1.;
#line 594 "MB02CV.f"
		dlarf_("Right", &len, q, &b[j + b_dim1], ldb, &tau, &g[
			g_offset], ldg, &dwork[1], (ftnlen)5);
#line 596 "MB02CV.f"
		b[j + b_dim1] = tau;
#line 597 "MB02CV.f"
	    }
#line 598 "MB02CV.f"
	    if (*q > 0) {

/*              Apply hyperbolic rotation. */

#line 602 "MB02CV.f"
		c__ = cs[(j << 1) - 1];
#line 603 "MB02CV.f"
		s = cs[j * 2];
#line 604 "MB02CV.f"
		d__1 = 1. / c__;
#line 604 "MB02CV.f"
		dscal_(&len, &d__1, &f1[j * f1_dim1 + 1], &c__1);
#line 605 "MB02CV.f"
		d__1 = -s / c__;
#line 605 "MB02CV.f"
		daxpy_(&len, &d__1, &g[g_offset], &c__1, &f1[j * f1_dim1 + 1],
			 &c__1);
#line 606 "MB02CV.f"
		dscal_(&len, &c__, &g[g_offset], &c__1);
#line 607 "MB02CV.f"
		d__1 = -s;
#line 607 "MB02CV.f"
		daxpy_(&len, &d__1, &f1[j * f1_dim1 + 1], &c__1, &g[g_offset],
			 &c__1);
#line 608 "MB02CV.f"
		if (ltri) {
#line 609 "MB02CV.f"
		    ++len;
#line 610 "MB02CV.f"
		    g[len + g_dim1] = -s / c__ * f1[len + j * f1_dim1];
#line 611 "MB02CV.f"
		    f1[len + j * f1_dim1] /= c__;

#line 613 "MB02CV.f"
		    i__2 = *q;
#line 613 "MB02CV.f"
		    for (jj = 2; jj <= i__2; ++jj) {
#line 614 "MB02CV.f"
			g[len + jj * g_dim1] = 0.;
#line 615 "MB02CV.f"
/* L110: */
#line 615 "MB02CV.f"
		    }

#line 617 "MB02CV.f"
		}
#line 618 "MB02CV.f"
	    }
#line 619 "MB02CV.f"
/* L120: */
#line 619 "MB02CV.f"
	}

#line 621 "MB02CV.f"
    } else {

/*        Row oriented and not deficient generator. */

#line 625 "MB02CV.f"
	if (ltri) {
/* Computing MAX */
#line 626 "MB02CV.f"
	    i__1 = *n - *k;
#line 626 "MB02CV.f"
	    len = max(i__1,0);
#line 627 "MB02CV.f"
	} else {
#line 628 "MB02CV.f"
	    len = *n;
#line 629 "MB02CV.f"
	}

#line 631 "MB02CV.f"
	if (col2 > 0) {
#line 632 "MB02CV.f"
	    nbl = min(*nb,col2);
#line 633 "MB02CV.f"
	    if (nbl > 0) {

/*              Blocked version. */

#line 637 "MB02CV.f"
		i__1 = *k - nbl + 1;
#line 637 "MB02CV.f"
		i__2 = nbl;
#line 637 "MB02CV.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 638 "MB02CV.f"
		    i__3 = *k - i__ + 1;
#line 638 "MB02CV.f"
		    ib = min(i__3,nbl);
#line 639 "MB02CV.f"
		    i__3 = *n + *k;
#line 639 "MB02CV.f"
		    dlarft_("Forward", "Columnwise", &col2, &ib, &a2[i__ * 
			    a2_dim1 + 1], lda2, &cs[(*k << 2) + i__], &dwork[
			    1], &i__3, (ftnlen)7, (ftnlen)10);
#line 641 "MB02CV.f"
		    i__3 = *n + *k;
#line 641 "MB02CV.f"
		    i__4 = *n + *k;
#line 641 "MB02CV.f"
		    dlarfb_("Left", "Transpose", "Forward", "Columnwise", &
			    col2, &len, &ib, &a2[i__ * a2_dim1 + 1], lda2, &
			    dwork[1], &i__3, &f2[f2_offset], ldf2, &dwork[ib 
			    + 1], &i__4, (ftnlen)4, (ftnlen)9, (ftnlen)7, (
			    ftnlen)10);

#line 646 "MB02CV.f"
		    i__3 = i__ + ib - 1;
#line 646 "MB02CV.f"
		    for (j = i__; j <= i__3; ++j) {
#line 647 "MB02CV.f"
			tau = a2[j * a2_dim1 + 1];
#line 648 "MB02CV.f"
			a2[j * a2_dim1 + 1] = 1.;
/* Computing MIN */
#line 649 "MB02CV.f"
			i__5 = col2, i__6 = j - i__ + 1;
#line 649 "MB02CV.f"
			i__4 = min(i__5,i__6);
#line 649 "MB02CV.f"
			dlarf_("Left", &i__4, &len, &a2[j * a2_dim1 + 1], &
				c__1, &tau, &f2[f2_offset], ldf2, &dwork[1], (
				ftnlen)4);
#line 651 "MB02CV.f"
			a2[j * a2_dim1 + 1] = tau;
#line 652 "MB02CV.f"
			c__ = cs[(*k << 1) + (j << 1) - 1];
#line 653 "MB02CV.f"
			s = cs[(*k << 1) + (j << 1)];
#line 654 "MB02CV.f"
			drot_(&len, &f1[j + f1_dim1], ldf1, &f2[f2_offset], 
				ldf2, &c__, &s);
#line 655 "MB02CV.f"
			if (ltri) {
#line 656 "MB02CV.f"
			    ++len;
#line 657 "MB02CV.f"
			    temp = f1[j + len * f1_dim1];
#line 658 "MB02CV.f"
			    f1[j + len * f1_dim1] = c__ * temp;
#line 659 "MB02CV.f"
			    f2[len * f2_dim1 + 1] = -s * temp;

#line 661 "MB02CV.f"
			    i__4 = col2;
#line 661 "MB02CV.f"
			    for (jj = 2; jj <= i__4; ++jj) {
#line 662 "MB02CV.f"
				f2[jj + len * f2_dim1] = 0.;
#line 663 "MB02CV.f"
/* L130: */
#line 663 "MB02CV.f"
			    }

#line 665 "MB02CV.f"
			}
#line 666 "MB02CV.f"
/* L140: */
#line 666 "MB02CV.f"
		    }

#line 668 "MB02CV.f"
/* L150: */
#line 668 "MB02CV.f"
		}

#line 670 "MB02CV.f"
	    } else {
#line 671 "MB02CV.f"
		i__ = 1;
#line 672 "MB02CV.f"
	    }

/*           Unblocked version for the last or only block. */

#line 676 "MB02CV.f"
	    i__2 = *k;
#line 676 "MB02CV.f"
	    for (j = i__; j <= i__2; ++j) {
#line 677 "MB02CV.f"
		if (col2 > 1) {
#line 678 "MB02CV.f"
		    tau = a2[j * a2_dim1 + 1];
#line 679 "MB02CV.f"
		    a2[j * a2_dim1 + 1] = 1.;
#line 680 "MB02CV.f"
		    dlarf_("Left", &col2, &len, &a2[j * a2_dim1 + 1], &c__1, &
			    tau, &f2[f2_offset], ldf2, &dwork[1], (ftnlen)4);
#line 682 "MB02CV.f"
		    a2[j * a2_dim1 + 1] = tau;
#line 683 "MB02CV.f"
		}

#line 685 "MB02CV.f"
		c__ = cs[(*k << 1) + (j << 1) - 1];
#line 686 "MB02CV.f"
		s = cs[(*k << 1) + (j << 1)];
#line 687 "MB02CV.f"
		drot_(&len, &f1[j + f1_dim1], ldf1, &f2[f2_offset], ldf2, &
			c__, &s);
#line 688 "MB02CV.f"
		if (ltri) {
#line 689 "MB02CV.f"
		    ++len;
#line 690 "MB02CV.f"
		    temp = f1[j + len * f1_dim1];
#line 691 "MB02CV.f"
		    f1[j + len * f1_dim1] = c__ * temp;
#line 692 "MB02CV.f"
		    f2[len * f2_dim1 + 1] = -s * temp;

#line 694 "MB02CV.f"
		    i__1 = col2;
#line 694 "MB02CV.f"
		    for (jj = 2; jj <= i__1; ++jj) {
#line 695 "MB02CV.f"
			f2[jj + len * f2_dim1] = 0.;
#line 696 "MB02CV.f"
/* L160: */
#line 696 "MB02CV.f"
		    }

#line 698 "MB02CV.f"
		}
#line 699 "MB02CV.f"
/* L170: */
#line 699 "MB02CV.f"
	    }

#line 701 "MB02CV.f"
	    pst2 = *k * 5;
#line 702 "MB02CV.f"
	} else {
#line 703 "MB02CV.f"
	    pst2 = *k << 1;
#line 704 "MB02CV.f"
	}

#line 706 "MB02CV.f"
	if (ltri) {
#line 707 "MB02CV.f"
	    len = *n - *k;
#line 708 "MB02CV.f"
	} else {
#line 709 "MB02CV.f"
	    len = *n;
#line 710 "MB02CV.f"
	}

#line 712 "MB02CV.f"
	nbl = min(*q,*nb);
#line 713 "MB02CV.f"
	if (nbl > 0) {

/*           Blocked version. */

#line 717 "MB02CV.f"
	    i__2 = *k - nbl + 1;
#line 717 "MB02CV.f"
	    i__1 = nbl;
#line 717 "MB02CV.f"
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 718 "MB02CV.f"
		i__3 = *k - i__ + 1;
#line 718 "MB02CV.f"
		ib = min(i__3,nbl);
#line 719 "MB02CV.f"
		i__3 = *n + *k;
#line 719 "MB02CV.f"
		dlarft_("Forward", "Columnwise", q, &ib, &b[i__ * b_dim1 + 1],
			 ldb, &cs[pst2 + i__], &dwork[1], &i__3, (ftnlen)7, (
			ftnlen)10);
#line 721 "MB02CV.f"
		i__3 = *n + *k;
#line 721 "MB02CV.f"
		i__4 = *n + *k;
#line 721 "MB02CV.f"
		dlarfb_("Left", "Transpose", "Forward", "Columnwise", q, &len,
			 &ib, &b[i__ * b_dim1 + 1], ldb, &dwork[1], &i__3, &g[
			g_offset], ldg, &dwork[ib + 1], &i__4, (ftnlen)4, (
			ftnlen)9, (ftnlen)7, (ftnlen)10);

#line 726 "MB02CV.f"
		i__3 = i__ + ib - 1;
#line 726 "MB02CV.f"
		for (j = i__; j <= i__3; ++j) {
#line 727 "MB02CV.f"
		    tau = b[j * b_dim1 + 1];
#line 728 "MB02CV.f"
		    b[j * b_dim1 + 1] = 1.;
#line 729 "MB02CV.f"
		    i__4 = j - i__ + 1;
#line 729 "MB02CV.f"
		    dlarf_("Left", &i__4, &len, &b[j * b_dim1 + 1], &c__1, &
			    tau, &g[g_offset], ldg, &dwork[1], (ftnlen)4);
#line 731 "MB02CV.f"
		    b[j * b_dim1 + 1] = tau;

/*                 Apply hyperbolic rotation. */

#line 735 "MB02CV.f"
		    c__ = cs[(j << 1) - 1];
#line 736 "MB02CV.f"
		    s = cs[j * 2];
#line 737 "MB02CV.f"
		    d__1 = 1. / c__;
#line 737 "MB02CV.f"
		    dscal_(&len, &d__1, &f1[j + f1_dim1], ldf1);
#line 738 "MB02CV.f"
		    d__1 = -s / c__;
#line 738 "MB02CV.f"
		    daxpy_(&len, &d__1, &g[g_offset], ldg, &f1[j + f1_dim1], 
			    ldf1);
#line 739 "MB02CV.f"
		    dscal_(&len, &c__, &g[g_offset], ldg);
#line 740 "MB02CV.f"
		    d__1 = -s;
#line 740 "MB02CV.f"
		    daxpy_(&len, &d__1, &f1[j + f1_dim1], ldf1, &g[g_offset], 
			    ldg);
#line 741 "MB02CV.f"
		    if (ltri) {
#line 742 "MB02CV.f"
			++len;
#line 743 "MB02CV.f"
			g[len * g_dim1 + 1] = -s / c__ * f1[j + len * f1_dim1]
				;
#line 744 "MB02CV.f"
			f1[j + len * f1_dim1] /= c__;

#line 746 "MB02CV.f"
			i__4 = *q;
#line 746 "MB02CV.f"
			for (jj = 2; jj <= i__4; ++jj) {
#line 747 "MB02CV.f"
			    g[jj + len * g_dim1] = 0.;
#line 748 "MB02CV.f"
/* L180: */
#line 748 "MB02CV.f"
			}

#line 750 "MB02CV.f"
		    }
#line 751 "MB02CV.f"
/* L190: */
#line 751 "MB02CV.f"
		}

#line 753 "MB02CV.f"
/* L200: */
#line 753 "MB02CV.f"
	    }

#line 755 "MB02CV.f"
	} else {
#line 756 "MB02CV.f"
	    i__ = 1;
#line 757 "MB02CV.f"
	}

/*        Unblocked version for the last or only block. */

#line 761 "MB02CV.f"
	i__1 = *k;
#line 761 "MB02CV.f"
	for (j = i__; j <= i__1; ++j) {
#line 762 "MB02CV.f"
	    if (*q > 1) {
#line 763 "MB02CV.f"
		tau = b[j * b_dim1 + 1];
#line 764 "MB02CV.f"
		b[j * b_dim1 + 1] = 1.;
#line 765 "MB02CV.f"
		dlarf_("Left", q, &len, &b[j * b_dim1 + 1], &c__1, &tau, &g[
			g_offset], ldg, &dwork[1], (ftnlen)4);
#line 767 "MB02CV.f"
		b[j * b_dim1 + 1] = tau;
#line 768 "MB02CV.f"
	    }
#line 769 "MB02CV.f"
	    if (*q > 0) {

/*              Apply hyperbolic rotation. */

#line 773 "MB02CV.f"
		c__ = cs[(j << 1) - 1];
#line 774 "MB02CV.f"
		s = cs[j * 2];
#line 775 "MB02CV.f"
		d__1 = 1. / c__;
#line 775 "MB02CV.f"
		dscal_(&len, &d__1, &f1[j + f1_dim1], ldf1);
#line 776 "MB02CV.f"
		d__1 = -s / c__;
#line 776 "MB02CV.f"
		daxpy_(&len, &d__1, &g[g_offset], ldg, &f1[j + f1_dim1], ldf1)
			;
#line 777 "MB02CV.f"
		dscal_(&len, &c__, &g[g_offset], ldg);
#line 778 "MB02CV.f"
		d__1 = -s;
#line 778 "MB02CV.f"
		daxpy_(&len, &d__1, &f1[j + f1_dim1], ldf1, &g[g_offset], ldg)
			;
#line 779 "MB02CV.f"
		if (ltri) {
#line 780 "MB02CV.f"
		    ++len;
#line 781 "MB02CV.f"
		    g[len * g_dim1 + 1] = -s / c__ * f1[j + len * f1_dim1];
#line 782 "MB02CV.f"
		    f1[j + len * f1_dim1] /= c__;

#line 784 "MB02CV.f"
		    i__2 = *q;
#line 784 "MB02CV.f"
		    for (jj = 2; jj <= i__2; ++jj) {
#line 785 "MB02CV.f"
			g[jj + len * g_dim1] = 0.;
#line 786 "MB02CV.f"
/* L210: */
#line 786 "MB02CV.f"
		    }

#line 788 "MB02CV.f"
		}
#line 789 "MB02CV.f"
	    }
#line 790 "MB02CV.f"
/* L220: */
#line 790 "MB02CV.f"
	}

#line 792 "MB02CV.f"
    }

/* *** Last line of MB02CV *** */
#line 795 "MB02CV.f"
    return 0;
} /* mb02cv_ */

