#line 1 "MB02CU.f"
/* MB02CU.f -- translated by f2c (version 20100827).
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

#line 1 "MB02CU.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mb02cu_(char *typeg, integer *k, integer *p, integer *q, 
	integer *nb, doublereal *a1, integer *lda1, doublereal *a2, integer *
	lda2, doublereal *b, integer *ldb, integer *rnk, integer *ipvt, 
	doublereal *cs, doublereal *tol, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen typeg_len)
{
    /* System generated locals */
    integer a1_dim1, a1_offset, a2_dim1, a2_offset, b_dim1, b_offset, i__1, 
	    i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal c__;
    static integer i__, j;
    static doublereal s;
    static integer ib, jj, nbl, len, pdw, phv, pos, pvt, col2;
    static doublereal tau1, tau2;
    static integer pst2;
    static doublereal beta;
    static logical lcol;
    static doublereal dmax__;
    static integer imax, ierr;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal temp2;
    extern /* Subroutine */ int ma02fd_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    static doublereal alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dlarf_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    ftnlen);
    static logical lrdef;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer itemp;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal alpha2;
    extern /* Subroutine */ int dgelq2_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *), dgeqr2_(
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *);
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int dlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen), dlarfg_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlarft_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen), dlartg_(doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *), xerbla_(char *, integer *, ftnlen);
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

/*     To bring the first blocks of a generator to proper form. */
/*     The positive part of the generator is contained in the arrays A1 */
/*     and A2. The negative part of the generator is contained in B. */
/*     Transformation information will be stored and can be applied */
/*     via SLICOT Library routine MB02CV. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYPEG   CHARACTER*1 */
/*             Specifies the type of the generator, as follows: */
/*             = 'D':  generator is column oriented and rank */
/*                     deficiencies are expected; */
/*             = 'C':  generator is column oriented and rank */
/*                     deficiencies are not expected; */
/*             = 'R':  generator is row oriented and rank */
/*                     deficiencies are not expected. */

/*     Input/Output Parameters */

/*     K       (input)  INTEGER */
/*             The number of rows in A1 to be processed.  K >= 0. */

/*     P       (input)  INTEGER */
/*             The number of columns of the positive generator.  P >= K. */

/*     Q       (input)  INTEGER */
/*             The number of columns in B containing the negative */
/*             generators. */
/*             If TYPEG = 'D',        Q >= K; */
/*             If TYPEG = 'C' or 'R', Q >= 0. */

/*     NB      (input)  INTEGER */
/*             On entry, if TYPEG = 'C'  or  TYPEG = 'R', NB specifies */
/*             the block size to be used in the blocked parts of the */
/*             algorithm. If NB <= 0, an unblocked algorithm is used. */

/*     A1      (input/output)  DOUBLE PRECISION array, dimension */
/*             (LDA1, K) */
/*             On entry, the leading K-by-K part of this array must */
/*             contain the leading submatrix of the positive part of the */
/*             generator. If TYPEG = 'C', A1 is assumed to be lower */
/*             triangular and the strictly upper triangular part is not */
/*             referenced. If TYPEG = 'R', A1 is assumed to be upper */
/*             triangular and the strictly lower triangular part is not */
/*             referenced. */
/*             On exit, if TYPEG = 'D', the leading K-by-RNK part of this */
/*             array contains the lower trapezoidal part of the proper */
/*             generator and information for the Householder */
/*             transformations applied during the reduction process. */
/*             On exit, if TYPEG = 'C', the leading K-by-K part of this */
/*             array contains the leading lower triangular part of the */
/*             proper generator. */
/*             On exit, if TYPEG = 'R', the leading K-by-K part of this */
/*             array contains the leading upper triangular part of the */
/*             proper generator. */

/*     LDA1    INTEGER */
/*             The leading dimension of the array A1.  LDA1 >= MAX(1,K). */

/*     A2      (input/output)  DOUBLE PRECISION array, */
/*             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDA2, P-K); */
/*             if TYPEG = 'R',                   dimension (LDA2, K). */
/*             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             K-by-(P-K) part of this array must contain the (K+1)-st */
/*             to P-th columns of the positive part of the generator. */
/*             On entry, if TYPEG = 'R', the leading (P-K)-by-K part of */
/*             this array must contain the (K+1)-st to P-th rows of the */
/*             positive part of the generator. */
/*             On exit, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             K-by-(P-K) part of this array contains information for */
/*             Householder transformations. */
/*             On exit, if TYPEG = 'R', the leading (P-K)-by-K part of */
/*             this array contains information for Householder */
/*             transformations. */

/*     LDA2    INTEGER */
/*             The leading dimension of the array A2. */
/*             If P = K,                   LDA2 >= 1; */
/*             If P > K and (TYPEG = 'D' or TYPEG = 'C'), */
/*                                         LDA2 >= MAX(1,K); */
/*             if P > K and TYPEG = 'R',   LDA2 >= P-K. */

/*     B       (input/output)  DOUBLE PRECISION array, */
/*             if TYPEG = 'D'  or  TYPEG = 'C',  dimension (LDB, Q); */
/*             if TYPEG = 'R',                   dimension (LDB, K). */
/*             On entry, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             K-by-Q part of this array must contain the negative part */
/*             of the generator. */
/*             On entry, if TYPEG = 'R', the leading Q-by-K part of this */
/*             array must contain the negative part of the generator. */
/*             On exit, if TYPEG = 'D'  or  TYPEG = 'C', the leading */
/*             K-by-Q part of this array contains information for */
/*             Householder transformations. */
/*             On exit, if TYPEG = 'R', the leading Q-by-K part of this */
/*             array contains information for Householder transformations. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             If Q = 0,                  LDB >= 1; */
/*             if Q > 0 and (TYPEG = 'D' or TYPEG = 'C'), */
/*                                        LDB >= MAX(1,K); */
/*             if Q > 0 and TYPEG = 'R',  LDB >= Q. */

/*     RNK     (output)  INTEGER */
/*             If TYPEG = 'D', the number of columns in the reduced */
/*             generator which are found to be linearly independent. */
/*             If TYPEG = 'C' or TYPEG = 'R', then RNK is not set. */

/*     IPVT    (output)  INTEGER array, dimension (K) */
/*             If TYPEG = 'D', then if IPVT(i) = k, the k-th row of the */
/*             proper generator is the reduced i-th row of the input */
/*             generator. */
/*             If TYPEG = 'C' or TYPEG = 'R', this array is not */
/*             referenced. */

/*     CS      (output)  DOUBLE PRECISION array, dimension (x) */
/*             If TYPEG = 'D' and P = K,                   x = 3*K; */
/*             if TYPEG = 'D' and P > K,                   x = 5*K; */
/*             if (TYPEG = 'C' or TYPEG = 'R') and P = K,  x = 4*K; */
/*             if (TYPEG = 'C' or TYPEG = 'R') and P > K,  x = 6*K. */
/*             On exit, the first x elements of this array contain */
/*             necessary information for the SLICOT library routine */
/*             MB02CV (Givens and modified hyperbolic rotation */
/*             parameters, scalar factors of the Householder */
/*             transformations). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If TYPEG = 'D', this number specifies the used tolerance */
/*             for handling deficiencies. If the hyperbolic norm */
/*             of two diagonal elements in the positive and negative */
/*             generators appears to be less than or equal to TOL, then */
/*             the corresponding columns are not reduced. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = -17,  DWORK(1) returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,4*K),         if TYPEG = 'D'; */
/*             LDWORK >= MAX(1,MAX(NB,1)*K), if TYPEG = 'C' or 'R'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if TYPEG = 'D', the generator represents a */
/*                   (numerically) indefinite matrix; and if TYPEG = 'C' */
/*                   or TYPEG = 'R', the generator represents a */
/*                   (numerically) semidefinite matrix. */

/*     METHOD */

/*     If TYPEG = 'C' or TYPEG = 'R', blocked Householder transformations */
/*     and modified hyperbolic rotations are used to downdate the */
/*     matrix [ A1  A2  sqrt(-1)*B ], cf. [1], [2]. */
/*     If TYPEG = 'D', then an algorithm with row pivoting is used. In */
/*     the first stage it maximizes the hyperbolic norm of the active */
/*     row. As soon as the hyperbolic norm is below the threshold TOL, */
/*     the strategy is changed. Now, in the second stage, the algorithm */
/*     applies an LQ decomposition with row pivoting on B such that */
/*     the Euclidean norm of the active row is maximized. */

/*     REFERENCES */

/*     [1] Kailath, T. and Sayed, A. */
/*         Fast Reliable Algorithms for Matrices with Structure. */
/*         SIAM Publications, Philadelphia, 1999. */

/*     [2] Kressner, D. and Van Dooren, P. */
/*         Factorizations and linear system solvers for matrices with */
/*         Toeplitz structure. */
/*         SLICOT Working Note 2000-2, 2000. */

/*     NUMERICAL ASPECTS */
/*                               2 */
/*     The algorithm requires 0(K *( P + Q )) floating point operations. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, May 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, June 2001. */
/*     D. Kressner, Technical Univ. Berlin, Germany, July 2002. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2004. */

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

#line 260 "MB02CU.f"
    /* Parameter adjustments */
#line 260 "MB02CU.f"
    a1_dim1 = *lda1;
#line 260 "MB02CU.f"
    a1_offset = 1 + a1_dim1;
#line 260 "MB02CU.f"
    a1 -= a1_offset;
#line 260 "MB02CU.f"
    a2_dim1 = *lda2;
#line 260 "MB02CU.f"
    a2_offset = 1 + a2_dim1;
#line 260 "MB02CU.f"
    a2 -= a2_offset;
#line 260 "MB02CU.f"
    b_dim1 = *ldb;
#line 260 "MB02CU.f"
    b_offset = 1 + b_dim1;
#line 260 "MB02CU.f"
    b -= b_offset;
#line 260 "MB02CU.f"
    --ipvt;
#line 260 "MB02CU.f"
    --cs;
#line 260 "MB02CU.f"
    --dwork;
#line 260 "MB02CU.f"

#line 260 "MB02CU.f"
    /* Function Body */
#line 260 "MB02CU.f"
    *info = 0;
#line 261 "MB02CU.f"
    col2 = *p - *k;
#line 262 "MB02CU.f"
    lrdef = lsame_(typeg, "D", (ftnlen)1, (ftnlen)1);
#line 263 "MB02CU.f"
    lcol = lsame_(typeg, "C", (ftnlen)1, (ftnlen)1);
#line 264 "MB02CU.f"
    if (lrdef) {
/* Computing MAX */
#line 265 "MB02CU.f"
	i__1 = 1, i__2 = *k << 2;
#line 265 "MB02CU.f"
	wrkmin = max(i__1,i__2);
#line 266 "MB02CU.f"
    } else {
/* Computing MAX */
#line 267 "MB02CU.f"
	i__1 = 1, i__2 = *nb * *k, i__1 = max(i__1,i__2);
#line 267 "MB02CU.f"
	wrkmin = max(i__1,*k);
#line 268 "MB02CU.f"
    }

/*     Check the scalar input parameters. */

#line 272 "MB02CU.f"
    if (! (lcol || lrdef || lsame_(typeg, "R", (ftnlen)1, (ftnlen)1))) {
#line 273 "MB02CU.f"
	*info = -1;
#line 274 "MB02CU.f"
    } else if (*k < 0) {
#line 275 "MB02CU.f"
	*info = -2;
#line 276 "MB02CU.f"
    } else if (*p < *k) {
#line 277 "MB02CU.f"
	*info = -3;
#line 278 "MB02CU.f"
    } else if (*q < 0 || lrdef && *q < *k) {
#line 279 "MB02CU.f"
	*info = -4;
#line 280 "MB02CU.f"
    } else if (*lda1 < max(1,*k)) {
#line 281 "MB02CU.f"
	*info = -7;
#line 282 "MB02CU.f"
    } else if (*p == *k && *lda2 < 1 || *p > *k && (lrdef || lcol) && *lda2 < 
	    max(1,*k) || *p > *k && ! (lrdef || lcol) && *lda2 < *p - *k) {
#line 287 "MB02CU.f"
	*info = -9;
#line 288 "MB02CU.f"
    } else if (*q == 0 && *ldb < 1 || *q > 0 && (lrdef || lcol) && *ldb < max(
	    1,*k) || *q > 0 && ! (lrdef || lcol) && *ldb < *q) {
#line 293 "MB02CU.f"
	*info = -11;
#line 294 "MB02CU.f"
    } else if (*ldwork < wrkmin) {
#line 295 "MB02CU.f"
	dwork[1] = (doublereal) wrkmin;
#line 296 "MB02CU.f"
	*info = -17;
#line 297 "MB02CU.f"
    }

/*     Return if there were illegal values. */

#line 301 "MB02CU.f"
    if (*info != 0) {
#line 302 "MB02CU.f"
	i__1 = -(*info);
#line 302 "MB02CU.f"
	xerbla_("MB02CU", &i__1, (ftnlen)6);
#line 303 "MB02CU.f"
	return 0;
#line 304 "MB02CU.f"
    }

/*     Quick return if possible. */

#line 308 "MB02CU.f"
    if (*k == 0 || ! lrdef && *q == 0 && *p == *k) {
#line 309 "MB02CU.f"
	if (lrdef) {
#line 309 "MB02CU.f"
	    *rnk = 0;
#line 309 "MB02CU.f"
	}
#line 311 "MB02CU.f"
	return 0;
#line 312 "MB02CU.f"
    }

#line 314 "MB02CU.f"
    if (lrdef) {

/*        Deficient generator. */

#line 318 "MB02CU.f"
	if (col2 == 0) {
#line 319 "MB02CU.f"
	    pst2 = *k << 1;
#line 320 "MB02CU.f"
	} else {
#line 321 "MB02CU.f"
	    pst2 = *k << 2;
#line 322 "MB02CU.f"
	}

/*        Initialize partial hyperbolic row norms. */

#line 326 "MB02CU.f"
	*rnk = 0;
#line 327 "MB02CU.f"
	phv = *k * 3;

#line 329 "MB02CU.f"
	i__1 = *k;
#line 329 "MB02CU.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 330 "MB02CU.f"
	    ipvt[i__] = i__;
#line 331 "MB02CU.f"
	    dwork[i__] = dnrm2_(k, &a1[i__ + a1_dim1], lda1);
#line 332 "MB02CU.f"
/* L10: */
#line 332 "MB02CU.f"
	}

#line 334 "MB02CU.f"
	i__1 = *k;
#line 334 "MB02CU.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 335 "MB02CU.f"
	    d__1 = dnrm2_(&col2, &a2[i__ + a2_dim1], lda2);
#line 335 "MB02CU.f"
	    dwork[i__] = dlapy2_(&dwork[i__], &d__1);
#line 337 "MB02CU.f"
	    dwork[i__ + *k] = dwork[i__];
#line 338 "MB02CU.f"
/* L20: */
#line 338 "MB02CU.f"
	}

#line 340 "MB02CU.f"
	pdw = *k << 1;

#line 342 "MB02CU.f"
	i__1 = *k;
#line 342 "MB02CU.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 343 "MB02CU.f"
	    ++pdw;
#line 344 "MB02CU.f"
	    dwork[pdw] = dnrm2_(q, &b[i__ + b_dim1], ldb);
#line 345 "MB02CU.f"
/* L30: */
#line 345 "MB02CU.f"
	}

/*        Compute factorization. */

#line 349 "MB02CU.f"
	i__1 = *k;
#line 349 "MB02CU.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Determine pivot row and swap if necessary. */

#line 353 "MB02CU.f"
	    pdw = i__;
#line 354 "MB02CU.f"
	    alpha = (d__1 = dwork[pdw], abs(d__1));
#line 355 "MB02CU.f"
	    beta = (d__1 = dwork[pdw + (*k << 1)], abs(d__1));
#line 356 "MB02CU.f"
	    d__2 = sqrt((d__1 = alpha - beta, abs(d__1))) * sqrt(alpha + beta)
		    ;
#line 356 "MB02CU.f"
	    d__3 = alpha - beta;
#line 356 "MB02CU.f"
	    dmax__ = d_sign(&d__2, &d__3);
#line 358 "MB02CU.f"
	    imax = i__;

#line 360 "MB02CU.f"
	    i__2 = *k - i__;
#line 360 "MB02CU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 361 "MB02CU.f"
		++pdw;
#line 362 "MB02CU.f"
		alpha = (d__1 = dwork[pdw], abs(d__1));
#line 363 "MB02CU.f"
		beta = (d__1 = dwork[pdw + (*k << 1)], abs(d__1));
#line 364 "MB02CU.f"
		d__2 = sqrt((d__1 = alpha - beta, abs(d__1))) * sqrt(alpha + 
			beta);
#line 364 "MB02CU.f"
		d__3 = alpha - beta;
#line 364 "MB02CU.f"
		temp = d_sign(&d__2, &d__3);
#line 366 "MB02CU.f"
		if (temp > dmax__) {
#line 367 "MB02CU.f"
		    imax = i__ + j;
#line 368 "MB02CU.f"
		    dmax__ = temp;
#line 369 "MB02CU.f"
		}
#line 370 "MB02CU.f"
/* L40: */
#line 370 "MB02CU.f"
	    }

/*           Proceed with the reduction if the hyperbolic norm is */
/*           beyond the threshold. */

#line 375 "MB02CU.f"
	    if (dmax__ > *tol) {

#line 377 "MB02CU.f"
		pvt = imax;
#line 378 "MB02CU.f"
		if (pvt != i__) {
#line 379 "MB02CU.f"
		    dswap_(k, &a1[pvt + a1_dim1], lda1, &a1[i__ + a1_dim1], 
			    lda1);
#line 380 "MB02CU.f"
		    dswap_(&col2, &a2[pvt + a2_dim1], lda2, &a2[i__ + a2_dim1]
			    , lda2);
#line 381 "MB02CU.f"
		    dswap_(q, &b[pvt + b_dim1], ldb, &b[i__ + b_dim1], ldb);
#line 382 "MB02CU.f"
		    itemp = ipvt[pvt];
#line 383 "MB02CU.f"
		    ipvt[pvt] = ipvt[i__];
#line 384 "MB02CU.f"
		    ipvt[i__] = itemp;
#line 385 "MB02CU.f"
		    dwork[pvt] = dwork[i__];
#line 386 "MB02CU.f"
		    dwork[*k + pvt] = dwork[*k + i__];
#line 387 "MB02CU.f"
		    dwork[(*k << 1) + pvt] = dwork[(*k << 1) + i__];
#line 388 "MB02CU.f"
		}

/*              Generate and apply elementary reflectors. */

#line 392 "MB02CU.f"
		if (col2 > 1) {
#line 393 "MB02CU.f"
		    dlarfg_(&col2, &a2[i__ + a2_dim1], &a2[i__ + (a2_dim1 << 
			    1)], lda2, &tau2);
#line 394 "MB02CU.f"
		    alpha2 = a2[i__ + a2_dim1];
#line 395 "MB02CU.f"
		    if (*k > i__) {
#line 396 "MB02CU.f"
			a2[i__ + a2_dim1] = 1.;
#line 397 "MB02CU.f"
			i__2 = *k - i__;
#line 397 "MB02CU.f"
			dlarf_("Right", &i__2, &col2, &a2[i__ + a2_dim1], 
				lda2, &tau2, &a2[i__ + 1 + a2_dim1], lda2, &
				dwork[phv + 1], (ftnlen)5);
#line 399 "MB02CU.f"
		    }
#line 400 "MB02CU.f"
		    a2[i__ + a2_dim1] = tau2;
#line 401 "MB02CU.f"
		} else if (col2 > 0) {
#line 402 "MB02CU.f"
		    alpha2 = a2[i__ + a2_dim1];
#line 403 "MB02CU.f"
		    a2[i__ + a2_dim1] = 0.;
#line 404 "MB02CU.f"
		}

#line 406 "MB02CU.f"
		if (*k > i__) {
#line 407 "MB02CU.f"
		    i__2 = *k - i__ + 1;
#line 407 "MB02CU.f"
		    dlarfg_(&i__2, &a1[i__ + i__ * a1_dim1], &a1[i__ + (i__ + 
			    1) * a1_dim1], lda1, &tau1);
#line 408 "MB02CU.f"
		    alpha = a1[i__ + i__ * a1_dim1];
#line 409 "MB02CU.f"
		    a1[i__ + i__ * a1_dim1] = 1.;
#line 410 "MB02CU.f"
		    i__2 = *k - i__;
#line 410 "MB02CU.f"
		    i__3 = *k - i__ + 1;
#line 410 "MB02CU.f"
		    dlarf_("Right", &i__2, &i__3, &a1[i__ + i__ * a1_dim1], 
			    lda1, &tau1, &a1[i__ + 1 + i__ * a1_dim1], lda1, &
			    dwork[phv + 1], (ftnlen)5);
#line 412 "MB02CU.f"
		    cs[pst2 + i__] = tau1;
#line 413 "MB02CU.f"
		} else {
#line 414 "MB02CU.f"
		    alpha = a1[i__ + i__ * a1_dim1];
#line 415 "MB02CU.f"
		}

#line 417 "MB02CU.f"
		if (col2 > 0) {
#line 418 "MB02CU.f"
		    temp = alpha;
#line 419 "MB02CU.f"
		    dlartg_(&temp, &alpha2, &c__, &s, &alpha);
#line 420 "MB02CU.f"
		    if (*k > i__) {
#line 420 "MB02CU.f"
			i__2 = *k - i__;
#line 420 "MB02CU.f"
			drot_(&i__2, &a1[i__ + 1 + i__ * a1_dim1], &c__1, &a2[
				i__ + 1 + a2_dim1], &c__1, &c__, &s);
#line 420 "MB02CU.f"
		    }
#line 422 "MB02CU.f"
		    cs[(*k << 1) + (i__ << 1) - 1] = c__;
#line 423 "MB02CU.f"
		    cs[(*k << 1) + (i__ << 1)] = s;
#line 424 "MB02CU.f"
		}
#line 425 "MB02CU.f"
		a1[i__ + i__ * a1_dim1] = alpha;

#line 427 "MB02CU.f"
		if (*q > 1) {
#line 428 "MB02CU.f"
		    dlarfg_(q, &b[i__ + b_dim1], &b[i__ + (b_dim1 << 1)], ldb,
			     &tau2);
#line 429 "MB02CU.f"
		    beta = b[i__ + b_dim1];
#line 430 "MB02CU.f"
		    if (*k > i__) {
#line 431 "MB02CU.f"
			b[i__ + b_dim1] = 1.;
#line 432 "MB02CU.f"
			i__2 = *k - i__;
#line 432 "MB02CU.f"
			dlarf_("Right", &i__2, q, &b[i__ + b_dim1], ldb, &
				tau2, &b[i__ + 1 + b_dim1], ldb, &dwork[phv + 
				1], (ftnlen)5);
#line 434 "MB02CU.f"
		    }
#line 435 "MB02CU.f"
		    b[i__ + b_dim1] = tau2;
#line 436 "MB02CU.f"
		} else if (*q > 0) {
#line 437 "MB02CU.f"
		    beta = b[i__ + b_dim1];
#line 438 "MB02CU.f"
		    b[i__ + b_dim1] = 0.;
#line 439 "MB02CU.f"
		} else {
#line 440 "MB02CU.f"
		    beta = 0.;
#line 441 "MB02CU.f"
		}

/*              Create hyperbolic Givens rotation. */

#line 445 "MB02CU.f"
		ma02fd_(&a1[i__ + i__ * a1_dim1], &beta, &c__, &s, &ierr);
#line 446 "MB02CU.f"
		if (ierr != 0) {

/*                 Error return:  This should not happen. */

#line 450 "MB02CU.f"
		    *info = 1;
#line 451 "MB02CU.f"
		    return 0;
#line 452 "MB02CU.f"
		}

/*              Apply hyperbolic rotation. */

#line 456 "MB02CU.f"
		if (*k > i__) {
#line 457 "MB02CU.f"
		    i__2 = *k - i__;
#line 457 "MB02CU.f"
		    d__1 = 1. / c__;
#line 457 "MB02CU.f"
		    dscal_(&i__2, &d__1, &a1[i__ + 1 + i__ * a1_dim1], &c__1);
#line 458 "MB02CU.f"
		    i__2 = *k - i__;
#line 458 "MB02CU.f"
		    d__1 = -s / c__;
#line 458 "MB02CU.f"
		    daxpy_(&i__2, &d__1, &b[i__ + 1 + b_dim1], &c__1, &a1[i__ 
			    + 1 + i__ * a1_dim1], &c__1);
#line 459 "MB02CU.f"
		    i__2 = *k - i__;
#line 459 "MB02CU.f"
		    dscal_(&i__2, &c__, &b[i__ + 1 + b_dim1], &c__1);
#line 460 "MB02CU.f"
		    i__2 = *k - i__;
#line 460 "MB02CU.f"
		    d__1 = -s;
#line 460 "MB02CU.f"
		    daxpy_(&i__2, &d__1, &a1[i__ + 1 + i__ * a1_dim1], &c__1, 
			    &b[i__ + 1 + b_dim1], &c__1);
#line 461 "MB02CU.f"
		}
#line 462 "MB02CU.f"
		cs[(i__ << 1) - 1] = c__;
#line 463 "MB02CU.f"
		cs[i__ * 2] = s;

/*              Downdate the norms in A1. */

#line 467 "MB02CU.f"
		i__2 = *k;
#line 467 "MB02CU.f"
		for (j = i__ + 1; j <= i__2; ++j) {
/* Computing 2nd power */
#line 468 "MB02CU.f"
		    d__2 = (d__1 = a1[j + i__ * a1_dim1], abs(d__1)) / dwork[
			    j];
#line 468 "MB02CU.f"
		    temp = 1. - d__2 * d__2;
/* Computing 2nd power */
#line 469 "MB02CU.f"
		    d__1 = dwork[j] / dwork[*k + j];
#line 469 "MB02CU.f"
		    temp2 = temp * .05 * (d__1 * d__1) + 1.;
#line 471 "MB02CU.f"
		    if (temp2 == 1.) {
#line 472 "MB02CU.f"
			i__3 = *k - i__;
#line 472 "MB02CU.f"
			d__1 = dnrm2_(&i__3, &a1[j + (i__ + 1) * a1_dim1], 
				lda1);
#line 472 "MB02CU.f"
			d__2 = dnrm2_(&col2, &a2[j + a2_dim1], lda2);
#line 472 "MB02CU.f"
			dwork[j] = dlapy2_(&d__1, &d__2);
#line 474 "MB02CU.f"
			dwork[*k + j] = dwork[j];
#line 475 "MB02CU.f"
			dwork[(*k << 1) + j] = dnrm2_(q, &b[j + b_dim1], ldb);
#line 476 "MB02CU.f"
		    } else {
#line 477 "MB02CU.f"
			if (temp >= 0.) {
#line 478 "MB02CU.f"
			    dwork[j] *= sqrt(temp);
#line 479 "MB02CU.f"
			} else {
#line 480 "MB02CU.f"
			    dwork[j] = -dwork[j] * sqrt(-temp);
#line 481 "MB02CU.f"
			}
#line 482 "MB02CU.f"
		    }
#line 483 "MB02CU.f"
/* L50: */
#line 483 "MB02CU.f"
		}

#line 485 "MB02CU.f"
		++(*rnk);
#line 486 "MB02CU.f"
	    } else if (abs(dmax__) < *tol) {

/*              Displacement is positive semidefinite. */
/*              Do an LQ decomposition with pivoting of the leftover */
/*              negative part to find diagonal elements with almost zero */
/*              norm. These columns cannot be removed from the */
/*              generator. */

/*              Initialize norms. */

#line 496 "MB02CU.f"
		i__2 = *k;
#line 496 "MB02CU.f"
		for (j = i__; j <= i__2; ++j) {
#line 497 "MB02CU.f"
		    dwork[j] = dnrm2_(q, &b[j + b_dim1], ldb);
#line 498 "MB02CU.f"
		    dwork[j + *k] = dwork[j];
#line 499 "MB02CU.f"
/* L60: */
#line 499 "MB02CU.f"
		}

#line 501 "MB02CU.f"
		len = *q;
#line 502 "MB02CU.f"
		pos = 1;

#line 504 "MB02CU.f"
		i__2 = *k;
#line 504 "MB02CU.f"
		for (j = i__; j <= i__2; ++j) {

/*                 Generate and apply elementary reflectors. */

#line 508 "MB02CU.f"
		    i__3 = *k - j + 1;
#line 508 "MB02CU.f"
		    pvt = j - 1 + idamax_(&i__3, &dwork[j], &c__1);

/*                 Swap rows if necessary. */

#line 512 "MB02CU.f"
		    if (pvt != j) {
#line 513 "MB02CU.f"
			dswap_(k, &a1[pvt + a1_dim1], lda1, &a1[j + a1_dim1], 
				lda1);
#line 514 "MB02CU.f"
			dswap_(&col2, &a2[pvt + a2_dim1], lda2, &a2[j + 
				a2_dim1], lda2);
#line 515 "MB02CU.f"
			dswap_(q, &b[pvt + b_dim1], ldb, &b[j + b_dim1], ldb);
#line 516 "MB02CU.f"
			itemp = ipvt[pvt];
#line 517 "MB02CU.f"
			ipvt[pvt] = ipvt[j];
#line 518 "MB02CU.f"
			ipvt[j] = itemp;
#line 519 "MB02CU.f"
			dwork[pvt] = dwork[j];
#line 520 "MB02CU.f"
			dwork[*k + pvt] = dwork[*k + j];
#line 521 "MB02CU.f"
		    }

/*                 Annihilate second part of the positive generators. */

#line 525 "MB02CU.f"
		    if (col2 > 1) {
#line 526 "MB02CU.f"
			dlarfg_(&col2, &a2[j + a2_dim1], &a2[j + (a2_dim1 << 
				1)], lda2, &tau2);
#line 527 "MB02CU.f"
			alpha2 = a2[j + a2_dim1];
#line 528 "MB02CU.f"
			if (*k > j) {
#line 529 "MB02CU.f"
			    a2[j + a2_dim1] = 1.;
#line 530 "MB02CU.f"
			    i__3 = *k - j;
#line 530 "MB02CU.f"
			    dlarf_("Right", &i__3, &col2, &a2[j + a2_dim1], 
				    lda2, &tau2, &a2[j + 1 + a2_dim1], lda2, &
				    dwork[phv + 1], (ftnlen)5);
#line 532 "MB02CU.f"
			}
#line 533 "MB02CU.f"
			a2[j + a2_dim1] = tau2;
#line 534 "MB02CU.f"
		    } else if (col2 > 0) {
#line 535 "MB02CU.f"
			alpha2 = a2[j + a2_dim1];
#line 536 "MB02CU.f"
			a2[j + a2_dim1] = 0.;
#line 537 "MB02CU.f"
		    }

/*                 Transform first part of the positive generators to */
/*                 lower triangular form. */

#line 542 "MB02CU.f"
		    if (*k > j) {
#line 543 "MB02CU.f"
			i__3 = *k - j + 1;
#line 543 "MB02CU.f"
			dlarfg_(&i__3, &a1[j + j * a1_dim1], &a1[j + (j + 1) *
				 a1_dim1], lda1, &tau1);
#line 545 "MB02CU.f"
			alpha = a1[j + j * a1_dim1];
#line 546 "MB02CU.f"
			a1[j + j * a1_dim1] = 1.;
#line 547 "MB02CU.f"
			i__3 = *k - j;
#line 547 "MB02CU.f"
			i__4 = *k - j + 1;
#line 547 "MB02CU.f"
			dlarf_("Right", &i__3, &i__4, &a1[j + j * a1_dim1], 
				lda1, &tau1, &a1[j + 1 + j * a1_dim1], lda1, &
				dwork[phv + 1], (ftnlen)5);
#line 549 "MB02CU.f"
			cs[pst2 + j] = tau1;
#line 550 "MB02CU.f"
		    } else {
#line 551 "MB02CU.f"
			alpha = a1[j + j * a1_dim1];
#line 552 "MB02CU.f"
		    }

#line 554 "MB02CU.f"
		    if (col2 > 0) {
#line 555 "MB02CU.f"
			temp = alpha;
#line 556 "MB02CU.f"
			dlartg_(&temp, &alpha2, &c__, &s, &alpha);
#line 557 "MB02CU.f"
			if (*k > j) {
#line 557 "MB02CU.f"
			    i__3 = *k - j;
#line 557 "MB02CU.f"
			    drot_(&i__3, &a1[j + 1 + j * a1_dim1], &c__1, &a2[
				    j + 1 + a2_dim1], &c__1, &c__, &s);
#line 557 "MB02CU.f"
			}
#line 560 "MB02CU.f"
			cs[(*k << 1) + (j << 1) - 1] = c__;
#line 561 "MB02CU.f"
			cs[(*k << 1) + (j << 1)] = s;
#line 562 "MB02CU.f"
		    }
#line 563 "MB02CU.f"
		    a1[j + j * a1_dim1] = alpha;

/*                 Transform negative part to lower triangular form. */

#line 567 "MB02CU.f"
		    if (len > 1) {
#line 568 "MB02CU.f"
			dlarfg_(&len, &b[j + pos * b_dim1], &b[j + (pos + 1) *
				 b_dim1], ldb, &tau2);
#line 569 "MB02CU.f"
			beta = b[j + pos * b_dim1];
#line 570 "MB02CU.f"
			if (*k > j) {
#line 571 "MB02CU.f"
			    b[j + pos * b_dim1] = 1.;
#line 572 "MB02CU.f"
			    i__3 = *k - j;
#line 572 "MB02CU.f"
			    dlarf_("Right", &i__3, &len, &b[j + pos * b_dim1],
				     ldb, &tau2, &b[j + 1 + pos * b_dim1], 
				    ldb, &dwork[phv + 1], (ftnlen)5);
#line 574 "MB02CU.f"
			}
#line 575 "MB02CU.f"
			b[j + pos * b_dim1] = beta;
#line 576 "MB02CU.f"
			cs[(j << 1) - 1] = tau2;
#line 577 "MB02CU.f"
		    }

/*                 Downdate the norms of the rows in the negative part. */

#line 581 "MB02CU.f"
		    i__3 = *k;
#line 581 "MB02CU.f"
		    for (jj = j + 1; jj <= i__3; ++jj) {
#line 582 "MB02CU.f"
			if (dwork[jj] != 0.) {
/* Computing 2nd power */
#line 583 "MB02CU.f"
			    d__2 = (d__1 = b[jj + pos * b_dim1], abs(d__1)) / 
				    dwork[jj];
#line 583 "MB02CU.f"
			    temp = 1. - d__2 * d__2;
#line 585 "MB02CU.f"
			    temp = max(temp,0.);
/* Computing 2nd power */
#line 586 "MB02CU.f"
			    d__1 = dwork[jj] / dwork[*k + jj];
#line 586 "MB02CU.f"
			    temp2 = temp * .05 * (d__1 * d__1) + 1.;
#line 588 "MB02CU.f"
			    if (temp2 == 1.) {
#line 589 "MB02CU.f"
				i__4 = len - 1;
#line 589 "MB02CU.f"
				dwork[jj] = dnrm2_(&i__4, &b[jj + (pos + 1) * 
					b_dim1], ldb);
#line 590 "MB02CU.f"
				dwork[*k + jj] = dwork[jj];
#line 591 "MB02CU.f"
			    } else {
#line 592 "MB02CU.f"
				if (temp >= 0.) {
#line 593 "MB02CU.f"
				    dwork[jj] *= sqrt(temp);
#line 594 "MB02CU.f"
				} else {
#line 595 "MB02CU.f"
				    dwork[jj] = -dwork[jj] * sqrt(-temp);
#line 596 "MB02CU.f"
				}
#line 597 "MB02CU.f"
			    }
#line 598 "MB02CU.f"
			}
#line 599 "MB02CU.f"
/* L70: */
#line 599 "MB02CU.f"
		    }

#line 601 "MB02CU.f"
		    --len;
#line 602 "MB02CU.f"
		    ++pos;
#line 603 "MB02CU.f"
/* L80: */
#line 603 "MB02CU.f"
		}

#line 605 "MB02CU.f"
		return 0;
#line 606 "MB02CU.f"
	    } else {

/*              Error return: */

/*              Displacement is indefinite. */
/*              Due to roundoff error, positive semidefiniteness is */
/*              violated. This is a rather bad situation. There is no */
/*              meaningful way to continue the computations from this */
/*              point. */

#line 616 "MB02CU.f"
		*info = 1;
#line 617 "MB02CU.f"
		return 0;
#line 618 "MB02CU.f"
	    }
#line 619 "MB02CU.f"
/* L90: */
#line 619 "MB02CU.f"
	}

#line 621 "MB02CU.f"
    } else if (lcol) {

/*        Column oriented and not deficient generator. */

/*        Apply an LQ like hyperbolic/orthogonal blocked decomposition. */

#line 627 "MB02CU.f"
	if (col2 > 0) {
#line 628 "MB02CU.f"
	    nbl = min(col2,*nb);
#line 629 "MB02CU.f"
	    if (nbl > 0) {

/*              Blocked version. */

#line 633 "MB02CU.f"
		i__1 = *k - nbl + 1;
#line 633 "MB02CU.f"
		i__2 = nbl;
#line 633 "MB02CU.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 634 "MB02CU.f"
		    i__3 = *k - i__ + 1;
#line 634 "MB02CU.f"
		    ib = min(i__3,nbl);
#line 635 "MB02CU.f"
		    dgelq2_(&ib, &col2, &a2[i__ + a2_dim1], lda2, &cs[(*k << 
			    2) + i__], &dwork[1], &ierr);
#line 637 "MB02CU.f"
		    if (i__ + ib <= *k) {
#line 638 "MB02CU.f"
			dlarft_("Forward", "Rowwise", &col2, &ib, &a2[i__ + 
				a2_dim1], lda2, &cs[(*k << 2) + i__], &dwork[
				1], k, (ftnlen)7, (ftnlen)7);
#line 640 "MB02CU.f"
			i__3 = *k - i__ - ib + 1;
#line 640 "MB02CU.f"
			dlarfb_("Right", "No Transpose", "Forward", "Rowwise",
				 &i__3, &col2, &ib, &a2[i__ + a2_dim1], lda2, 
				&dwork[1], k, &a2[i__ + ib + a2_dim1], lda2, &
				dwork[ib + 1], k, (ftnlen)5, (ftnlen)12, (
				ftnlen)7, (ftnlen)7);
#line 644 "MB02CU.f"
		    }

/*                 Annihilate the remaining parts of A2. */

#line 648 "MB02CU.f"
		    i__3 = i__ + ib - 1;
#line 648 "MB02CU.f"
		    for (j = i__; j <= i__3; ++j) {
#line 649 "MB02CU.f"
			if (col2 > 1) {
/* Computing MIN */
#line 650 "MB02CU.f"
			    i__4 = col2, i__5 = j - i__ + 1;
#line 650 "MB02CU.f"
			    len = min(i__4,i__5);
#line 651 "MB02CU.f"
			    dlarfg_(&len, &a2[j + a2_dim1], &a2[j + (a2_dim1 
				    << 1)], lda2, &tau2);
#line 652 "MB02CU.f"
			    alpha2 = a2[j + a2_dim1];
#line 653 "MB02CU.f"
			    if (*k > j) {
#line 654 "MB02CU.f"
				a2[j + a2_dim1] = 1.;
#line 655 "MB02CU.f"
				i__4 = *k - j;
#line 655 "MB02CU.f"
				dlarf_("Right", &i__4, &len, &a2[j + a2_dim1],
					 lda2, &tau2, &a2[j + 1 + a2_dim1], 
					lda2, &dwork[1], (ftnlen)5);
#line 657 "MB02CU.f"
			    }
#line 658 "MB02CU.f"
			    a2[j + a2_dim1] = tau2;
#line 659 "MB02CU.f"
			} else {
#line 660 "MB02CU.f"
			    alpha2 = a2[j + a2_dim1];
#line 661 "MB02CU.f"
			    a2[j + a2_dim1] = 0.;
#line 662 "MB02CU.f"
			}
#line 663 "MB02CU.f"
			alpha = a1[j + j * a1_dim1];
#line 664 "MB02CU.f"
			dlartg_(&alpha, &alpha2, &c__, &s, &a1[j + j * 
				a1_dim1]);
#line 665 "MB02CU.f"
			if (*k > j) {
#line 665 "MB02CU.f"
			    i__4 = *k - j;
#line 665 "MB02CU.f"
			    drot_(&i__4, &a1[j + 1 + j * a1_dim1], &c__1, &a2[
				    j + 1 + a2_dim1], &c__1, &c__, &s);
#line 665 "MB02CU.f"
			}
#line 668 "MB02CU.f"
			cs[(*k << 1) + (j << 1) - 1] = c__;
#line 669 "MB02CU.f"
			cs[(*k << 1) + (j << 1)] = s;
#line 670 "MB02CU.f"
/* L100: */
#line 670 "MB02CU.f"
		    }

#line 672 "MB02CU.f"
/* L110: */
#line 672 "MB02CU.f"
		}

#line 674 "MB02CU.f"
	    } else {
#line 675 "MB02CU.f"
		i__ = 1;
#line 676 "MB02CU.f"
	    }

/*           Unblocked version for the last or only block. */

#line 680 "MB02CU.f"
	    i__2 = *k;
#line 680 "MB02CU.f"
	    for (j = i__; j <= i__2; ++j) {
#line 681 "MB02CU.f"
		if (col2 > 1) {
#line 682 "MB02CU.f"
		    dlarfg_(&col2, &a2[j + a2_dim1], &a2[j + (a2_dim1 << 1)], 
			    lda2, &tau2);
#line 683 "MB02CU.f"
		    alpha2 = a2[j + a2_dim1];
#line 684 "MB02CU.f"
		    if (*k > j) {
#line 685 "MB02CU.f"
			a2[j + a2_dim1] = 1.;
#line 686 "MB02CU.f"
			i__1 = *k - j;
#line 686 "MB02CU.f"
			dlarf_("Right", &i__1, &col2, &a2[j + a2_dim1], lda2, 
				&tau2, &a2[j + 1 + a2_dim1], lda2, &dwork[1], 
				(ftnlen)5);
#line 688 "MB02CU.f"
		    }
#line 689 "MB02CU.f"
		    a2[j + a2_dim1] = tau2;
#line 690 "MB02CU.f"
		} else {
#line 691 "MB02CU.f"
		    alpha2 = a2[j + a2_dim1];
#line 692 "MB02CU.f"
		    a2[j + a2_dim1] = 0.;
#line 693 "MB02CU.f"
		}
#line 694 "MB02CU.f"
		alpha = a1[j + j * a1_dim1];
#line 695 "MB02CU.f"
		dlartg_(&alpha, &alpha2, &c__, &s, &a1[j + j * a1_dim1]);
#line 696 "MB02CU.f"
		if (*k > j) {
#line 696 "MB02CU.f"
		    i__1 = *k - j;
#line 696 "MB02CU.f"
		    drot_(&i__1, &a1[j + 1 + j * a1_dim1], &c__1, &a2[j + 1 + 
			    a2_dim1], &c__1, &c__, &s);
#line 696 "MB02CU.f"
		}
#line 698 "MB02CU.f"
		cs[(*k << 1) + (j << 1) - 1] = c__;
#line 699 "MB02CU.f"
		cs[(*k << 1) + (j << 1)] = s;
#line 700 "MB02CU.f"
/* L120: */
#line 700 "MB02CU.f"
	    }

#line 702 "MB02CU.f"
	    pst2 = *k * 5;
#line 703 "MB02CU.f"
	} else {
#line 704 "MB02CU.f"
	    pst2 = *k << 1;
#line 705 "MB02CU.f"
	}

/*        Annihilate B with hyperbolic transformations. */

#line 709 "MB02CU.f"
	nbl = min(*nb,*q);
#line 710 "MB02CU.f"
	if (nbl > 0) {

/*           Blocked version. */

#line 714 "MB02CU.f"
	    i__2 = *k - nbl + 1;
#line 714 "MB02CU.f"
	    i__1 = nbl;
#line 714 "MB02CU.f"
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 715 "MB02CU.f"
		i__3 = *k - i__ + 1;
#line 715 "MB02CU.f"
		ib = min(i__3,nbl);
#line 716 "MB02CU.f"
		dgelq2_(&ib, q, &b[i__ + b_dim1], ldb, &cs[pst2 + i__], &
			dwork[1], &ierr);
#line 718 "MB02CU.f"
		if (i__ + ib <= *k) {
#line 719 "MB02CU.f"
		    dlarft_("Forward", "Rowwise", q, &ib, &b[i__ + b_dim1], 
			    ldb, &cs[pst2 + i__], &dwork[1], k, (ftnlen)7, (
			    ftnlen)7);
#line 721 "MB02CU.f"
		    i__3 = *k - i__ - ib + 1;
#line 721 "MB02CU.f"
		    dlarfb_("Right", "No Transpose", "Forward", "Rowwise", &
			    i__3, q, &ib, &b[i__ + b_dim1], ldb, &dwork[1], k,
			     &b[i__ + ib + b_dim1], ldb, &dwork[ib + 1], k, (
			    ftnlen)5, (ftnlen)12, (ftnlen)7, (ftnlen)7);
#line 725 "MB02CU.f"
		}

/*              Annihilate the remaining parts of B. */

#line 729 "MB02CU.f"
		i__3 = i__ + ib - 1;
#line 729 "MB02CU.f"
		for (j = i__; j <= i__3; ++j) {
#line 730 "MB02CU.f"
		    if (*q > 1) {
#line 731 "MB02CU.f"
			i__4 = j - i__ + 1;
#line 731 "MB02CU.f"
			dlarfg_(&i__4, &b[j + b_dim1], &b[j + (b_dim1 << 1)], 
				ldb, &tau2);
#line 732 "MB02CU.f"
			alpha2 = b[j + b_dim1];
#line 733 "MB02CU.f"
			if (*k > j) {
#line 734 "MB02CU.f"
			    b[j + b_dim1] = 1.;
#line 735 "MB02CU.f"
			    i__4 = *k - j;
#line 735 "MB02CU.f"
			    i__5 = j - i__ + 1;
#line 735 "MB02CU.f"
			    dlarf_("Right", &i__4, &i__5, &b[j + b_dim1], ldb,
				     &tau2, &b[j + 1 + b_dim1], ldb, &dwork[1]
				    , (ftnlen)5);
#line 737 "MB02CU.f"
			}
#line 738 "MB02CU.f"
			b[j + b_dim1] = tau2;
#line 739 "MB02CU.f"
		    } else {
#line 740 "MB02CU.f"
			alpha2 = b[j + b_dim1];
#line 741 "MB02CU.f"
			b[j + b_dim1] = 0.;
#line 742 "MB02CU.f"
		    }

/*                 Create hyperbolic rotation. */

#line 746 "MB02CU.f"
		    ma02fd_(&a1[j + j * a1_dim1], &alpha2, &c__, &s, &ierr);
#line 747 "MB02CU.f"
		    if (ierr != 0) {

/*                    Error return:  The matrix is not positive definite. */

#line 751 "MB02CU.f"
			*info = 1;
#line 752 "MB02CU.f"
			return 0;
#line 753 "MB02CU.f"
		    }

/*                 Apply hyperbolic rotation. */

#line 757 "MB02CU.f"
		    if (*k > j) {
#line 758 "MB02CU.f"
			i__4 = *k - j;
#line 758 "MB02CU.f"
			d__1 = 1. / c__;
#line 758 "MB02CU.f"
			dscal_(&i__4, &d__1, &a1[j + 1 + j * a1_dim1], &c__1);
#line 759 "MB02CU.f"
			i__4 = *k - j;
#line 759 "MB02CU.f"
			d__1 = -s / c__;
#line 759 "MB02CU.f"
			daxpy_(&i__4, &d__1, &b[j + 1 + b_dim1], &c__1, &a1[j 
				+ 1 + j * a1_dim1], &c__1);
#line 760 "MB02CU.f"
			i__4 = *k - j;
#line 760 "MB02CU.f"
			dscal_(&i__4, &c__, &b[j + 1 + b_dim1], &c__1);
#line 761 "MB02CU.f"
			i__4 = *k - j;
#line 761 "MB02CU.f"
			d__1 = -s;
#line 761 "MB02CU.f"
			daxpy_(&i__4, &d__1, &a1[j + 1 + j * a1_dim1], &c__1, 
				&b[j + 1 + b_dim1], &c__1);
#line 762 "MB02CU.f"
		    }
#line 763 "MB02CU.f"
		    cs[(j << 1) - 1] = c__;
#line 764 "MB02CU.f"
		    cs[j * 2] = s;
#line 765 "MB02CU.f"
/* L130: */
#line 765 "MB02CU.f"
		}

#line 767 "MB02CU.f"
/* L140: */
#line 767 "MB02CU.f"
	    }

#line 769 "MB02CU.f"
	} else {
#line 770 "MB02CU.f"
	    i__ = 1;
#line 771 "MB02CU.f"
	}

/*        Unblocked version for the last or only block. */

#line 775 "MB02CU.f"
	i__1 = *k;
#line 775 "MB02CU.f"
	for (j = i__; j <= i__1; ++j) {
#line 776 "MB02CU.f"
	    if (*q > 1) {
#line 777 "MB02CU.f"
		dlarfg_(q, &b[j + b_dim1], &b[j + (b_dim1 << 1)], ldb, &tau2);
#line 778 "MB02CU.f"
		alpha2 = b[j + b_dim1];
#line 779 "MB02CU.f"
		if (*k > j) {
#line 780 "MB02CU.f"
		    b[j + b_dim1] = 1.;
#line 781 "MB02CU.f"
		    i__2 = *k - j;
#line 781 "MB02CU.f"
		    dlarf_("Right", &i__2, q, &b[j + b_dim1], ldb, &tau2, &b[
			    j + 1 + b_dim1], ldb, &dwork[1], (ftnlen)5);
#line 783 "MB02CU.f"
		}
#line 784 "MB02CU.f"
		b[j + b_dim1] = tau2;
#line 785 "MB02CU.f"
	    } else if (*q > 0) {
#line 786 "MB02CU.f"
		alpha2 = b[j + b_dim1];
#line 787 "MB02CU.f"
		b[j + b_dim1] = 0.;
#line 788 "MB02CU.f"
	    }
#line 789 "MB02CU.f"
	    if (*q > 0) {

/*              Create hyperbolic rotation. */

#line 793 "MB02CU.f"
		ma02fd_(&a1[j + j * a1_dim1], &alpha2, &c__, &s, &ierr);
#line 794 "MB02CU.f"
		if (ierr != 0) {

/*                 Error return:  The matrix is not positive definite. */

#line 798 "MB02CU.f"
		    *info = 1;
#line 799 "MB02CU.f"
		    return 0;
#line 800 "MB02CU.f"
		}

/*              Apply hyperbolic rotation. */

#line 804 "MB02CU.f"
		if (*k > j) {
#line 805 "MB02CU.f"
		    i__2 = *k - j;
#line 805 "MB02CU.f"
		    d__1 = 1. / c__;
#line 805 "MB02CU.f"
		    dscal_(&i__2, &d__1, &a1[j + 1 + j * a1_dim1], &c__1);
#line 806 "MB02CU.f"
		    i__2 = *k - j;
#line 806 "MB02CU.f"
		    d__1 = -s / c__;
#line 806 "MB02CU.f"
		    daxpy_(&i__2, &d__1, &b[j + 1 + b_dim1], &c__1, &a1[j + 1 
			    + j * a1_dim1], &c__1);
#line 807 "MB02CU.f"
		    i__2 = *k - j;
#line 807 "MB02CU.f"
		    dscal_(&i__2, &c__, &b[j + 1 + b_dim1], &c__1);
#line 808 "MB02CU.f"
		    i__2 = *k - j;
#line 808 "MB02CU.f"
		    d__1 = -s;
#line 808 "MB02CU.f"
		    daxpy_(&i__2, &d__1, &a1[j + 1 + j * a1_dim1], &c__1, &b[
			    j + 1 + b_dim1], &c__1);
#line 809 "MB02CU.f"
		}
#line 810 "MB02CU.f"
		cs[(j << 1) - 1] = c__;
#line 811 "MB02CU.f"
		cs[j * 2] = s;
#line 812 "MB02CU.f"
	    }
#line 813 "MB02CU.f"
/* L150: */
#line 813 "MB02CU.f"
	}

#line 815 "MB02CU.f"
    } else {

/*        Row oriented and not deficient generator. */

#line 819 "MB02CU.f"
	if (col2 > 0) {
#line 820 "MB02CU.f"
	    nbl = min(*nb,col2);
#line 821 "MB02CU.f"
	    if (nbl > 0) {

/*              Blocked version. */

#line 825 "MB02CU.f"
		i__1 = *k - nbl + 1;
#line 825 "MB02CU.f"
		i__2 = nbl;
#line 825 "MB02CU.f"
		for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += 
			i__2) {
/* Computing MIN */
#line 826 "MB02CU.f"
		    i__3 = *k - i__ + 1;
#line 826 "MB02CU.f"
		    ib = min(i__3,nbl);
#line 827 "MB02CU.f"
		    dgeqr2_(&col2, &ib, &a2[i__ * a2_dim1 + 1], lda2, &cs[(*k 
			    << 2) + i__], &dwork[1], &ierr);
#line 829 "MB02CU.f"
		    if (i__ + ib <= *k) {
#line 830 "MB02CU.f"
			dlarft_("Forward", "Columnwise", &col2, &ib, &a2[i__ *
				 a2_dim1 + 1], lda2, &cs[(*k << 2) + i__], &
				dwork[1], k, (ftnlen)7, (ftnlen)10);
#line 832 "MB02CU.f"
			i__3 = *k - i__ - ib + 1;
#line 832 "MB02CU.f"
			dlarfb_("Left", "Transpose", "Forward", "Columnwise", 
				&col2, &i__3, &ib, &a2[i__ * a2_dim1 + 1], 
				lda2, &dwork[1], k, &a2[(i__ + ib) * a2_dim1 
				+ 1], lda2, &dwork[ib + 1], k, (ftnlen)4, (
				ftnlen)9, (ftnlen)7, (ftnlen)10);
#line 836 "MB02CU.f"
		    }

/*                 Annihilate the remaining parts of A2. */

#line 840 "MB02CU.f"
		    i__3 = i__ + ib - 1;
#line 840 "MB02CU.f"
		    for (j = i__; j <= i__3; ++j) {
#line 841 "MB02CU.f"
			if (col2 > 1) {
/* Computing MIN */
#line 842 "MB02CU.f"
			    i__4 = col2, i__5 = j - i__ + 1;
#line 842 "MB02CU.f"
			    len = min(i__4,i__5);
#line 843 "MB02CU.f"
			    dlarfg_(&len, &a2[j * a2_dim1 + 1], &a2[j * 
				    a2_dim1 + 2], &c__1, &tau2);
#line 844 "MB02CU.f"
			    alpha2 = a2[j * a2_dim1 + 1];
#line 845 "MB02CU.f"
			    if (*k > j) {
#line 846 "MB02CU.f"
				a2[j * a2_dim1 + 1] = 1.;
#line 847 "MB02CU.f"
				i__4 = *k - j;
#line 847 "MB02CU.f"
				dlarf_("Left", &len, &i__4, &a2[j * a2_dim1 + 
					1], &c__1, &tau2, &a2[(j + 1) * 
					a2_dim1 + 1], lda2, &dwork[1], (
					ftnlen)4);
#line 849 "MB02CU.f"
			    }
#line 850 "MB02CU.f"
			    a2[j * a2_dim1 + 1] = tau2;
#line 851 "MB02CU.f"
			} else {
#line 852 "MB02CU.f"
			    alpha2 = a2[j * a2_dim1 + 1];
#line 853 "MB02CU.f"
			    a2[j * a2_dim1 + 1] = 0.;
#line 854 "MB02CU.f"
			}
#line 855 "MB02CU.f"
			alpha = a1[j + j * a1_dim1];
#line 856 "MB02CU.f"
			dlartg_(&alpha, &alpha2, &c__, &s, &a1[j + j * 
				a1_dim1]);
#line 857 "MB02CU.f"
			if (*k > j) {
#line 857 "MB02CU.f"
			    i__4 = *k - j;
#line 857 "MB02CU.f"
			    drot_(&i__4, &a1[j + (j + 1) * a1_dim1], lda1, &
				    a2[(j + 1) * a2_dim1 + 1], lda2, &c__, &s)
				    ;
#line 857 "MB02CU.f"
			}
#line 860 "MB02CU.f"
			cs[(*k << 1) + (j << 1) - 1] = c__;
#line 861 "MB02CU.f"
			cs[(*k << 1) + (j << 1)] = s;
#line 862 "MB02CU.f"
/* L160: */
#line 862 "MB02CU.f"
		    }

#line 864 "MB02CU.f"
/* L170: */
#line 864 "MB02CU.f"
		}

#line 866 "MB02CU.f"
	    } else {
#line 867 "MB02CU.f"
		i__ = 1;
#line 868 "MB02CU.f"
	    }

/*           Unblocked version for the last or only block. */

#line 872 "MB02CU.f"
	    i__2 = *k;
#line 872 "MB02CU.f"
	    for (j = i__; j <= i__2; ++j) {
#line 873 "MB02CU.f"
		if (col2 > 1) {
#line 874 "MB02CU.f"
		    dlarfg_(&col2, &a2[j * a2_dim1 + 1], &a2[j * a2_dim1 + 2],
			     &c__1, &tau2);
#line 875 "MB02CU.f"
		    alpha2 = a2[j * a2_dim1 + 1];
#line 876 "MB02CU.f"
		    if (*k > j) {
#line 877 "MB02CU.f"
			a2[j * a2_dim1 + 1] = 1.;
#line 878 "MB02CU.f"
			i__1 = *k - j;
#line 878 "MB02CU.f"
			dlarf_("Left", &col2, &i__1, &a2[j * a2_dim1 + 1], &
				c__1, &tau2, &a2[(j + 1) * a2_dim1 + 1], lda2,
				 &dwork[1], (ftnlen)4);
#line 880 "MB02CU.f"
		    }
#line 881 "MB02CU.f"
		    a2[j * a2_dim1 + 1] = tau2;
#line 882 "MB02CU.f"
		} else {
#line 883 "MB02CU.f"
		    alpha2 = a2[j * a2_dim1 + 1];
#line 884 "MB02CU.f"
		    a2[j * a2_dim1 + 1] = 0.;
#line 885 "MB02CU.f"
		}
#line 886 "MB02CU.f"
		alpha = a1[j + j * a1_dim1];
#line 887 "MB02CU.f"
		dlartg_(&alpha, &alpha2, &c__, &s, &a1[j + j * a1_dim1]);
#line 888 "MB02CU.f"
		if (*k > j) {
#line 888 "MB02CU.f"
		    i__1 = *k - j;
#line 888 "MB02CU.f"
		    drot_(&i__1, &a1[j + (j + 1) * a1_dim1], lda1, &a2[(j + 1)
			     * a2_dim1 + 1], lda2, &c__, &s);
#line 888 "MB02CU.f"
		}
#line 891 "MB02CU.f"
		cs[(*k << 1) + (j << 1) - 1] = c__;
#line 892 "MB02CU.f"
		cs[(*k << 1) + (j << 1)] = s;
#line 893 "MB02CU.f"
/* L180: */
#line 893 "MB02CU.f"
	    }

#line 895 "MB02CU.f"
	    pst2 = *k * 5;
#line 896 "MB02CU.f"
	} else {
#line 897 "MB02CU.f"
	    pst2 = *k << 1;
#line 898 "MB02CU.f"
	}

/*        Annihilate B with hyperbolic transformations. */

#line 902 "MB02CU.f"
	nbl = min(*nb,*q);
#line 903 "MB02CU.f"
	if (nbl > 0) {

/*           Blocked version. */

#line 907 "MB02CU.f"
	    i__2 = *k - nbl + 1;
#line 907 "MB02CU.f"
	    i__1 = nbl;
#line 907 "MB02CU.f"
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing MIN */
#line 908 "MB02CU.f"
		i__3 = *k - i__ + 1;
#line 908 "MB02CU.f"
		ib = min(i__3,nbl);
#line 909 "MB02CU.f"
		dgeqr2_(q, &ib, &b[i__ * b_dim1 + 1], ldb, &cs[pst2 + i__], &
			dwork[1], &ierr);
#line 911 "MB02CU.f"
		if (i__ + ib <= *k) {
#line 912 "MB02CU.f"
		    dlarft_("Forward", "Columnwise", q, &ib, &b[i__ * b_dim1 
			    + 1], ldb, &cs[pst2 + i__], &dwork[1], k, (ftnlen)
			    7, (ftnlen)10);
#line 914 "MB02CU.f"
		    i__3 = *k - i__ - ib + 1;
#line 914 "MB02CU.f"
		    dlarfb_("Left", "Transpose", "Forward", "Columnwise", q, &
			    i__3, &ib, &b[i__ * b_dim1 + 1], ldb, &dwork[1], 
			    k, &b[(i__ + ib) * b_dim1 + 1], ldb, &dwork[ib + 
			    1], k, (ftnlen)4, (ftnlen)9, (ftnlen)7, (ftnlen)
			    10);
#line 918 "MB02CU.f"
		}

/*              Annihilate the remaining parts of B. */

#line 922 "MB02CU.f"
		i__3 = i__ + ib - 1;
#line 922 "MB02CU.f"
		for (j = i__; j <= i__3; ++j) {
#line 923 "MB02CU.f"
		    if (*q > 1) {
#line 924 "MB02CU.f"
			i__4 = j - i__ + 1;
#line 924 "MB02CU.f"
			dlarfg_(&i__4, &b[j * b_dim1 + 1], &b[j * b_dim1 + 2],
				 &c__1, &tau2);
#line 925 "MB02CU.f"
			alpha2 = b[j * b_dim1 + 1];
#line 926 "MB02CU.f"
			if (*k > j) {
#line 927 "MB02CU.f"
			    b[j * b_dim1 + 1] = 1.;
#line 928 "MB02CU.f"
			    i__4 = j - i__ + 1;
#line 928 "MB02CU.f"
			    i__5 = *k - j;
#line 928 "MB02CU.f"
			    dlarf_("Left", &i__4, &i__5, &b[j * b_dim1 + 1], &
				    c__1, &tau2, &b[(j + 1) * b_dim1 + 1], 
				    ldb, &dwork[1], (ftnlen)4);
#line 930 "MB02CU.f"
			}
#line 931 "MB02CU.f"
			b[j * b_dim1 + 1] = tau2;
#line 932 "MB02CU.f"
		    } else {
#line 933 "MB02CU.f"
			alpha2 = b[j * b_dim1 + 1];
#line 934 "MB02CU.f"
			b[j * b_dim1 + 1] = 0.;
#line 935 "MB02CU.f"
		    }

/*                 Create hyperbolic rotation. */

#line 939 "MB02CU.f"
		    ma02fd_(&a1[j + j * a1_dim1], &alpha2, &c__, &s, &ierr);
#line 940 "MB02CU.f"
		    if (ierr != 0) {

/*                    Error return:  The matrix is not positive definite. */

#line 944 "MB02CU.f"
			*info = 1;
#line 945 "MB02CU.f"
			return 0;
#line 946 "MB02CU.f"
		    }

/*                 Apply hyperbolic rotation. */

#line 950 "MB02CU.f"
		    if (*k > j) {
#line 951 "MB02CU.f"
			i__4 = *k - j;
#line 951 "MB02CU.f"
			d__1 = 1. / c__;
#line 951 "MB02CU.f"
			dscal_(&i__4, &d__1, &a1[j + (j + 1) * a1_dim1], lda1)
				;
#line 952 "MB02CU.f"
			i__4 = *k - j;
#line 952 "MB02CU.f"
			d__1 = -s / c__;
#line 952 "MB02CU.f"
			daxpy_(&i__4, &d__1, &b[(j + 1) * b_dim1 + 1], ldb, &
				a1[j + (j + 1) * a1_dim1], lda1);
#line 954 "MB02CU.f"
			i__4 = *k - j;
#line 954 "MB02CU.f"
			dscal_(&i__4, &c__, &b[(j + 1) * b_dim1 + 1], ldb);
#line 955 "MB02CU.f"
			i__4 = *k - j;
#line 955 "MB02CU.f"
			d__1 = -s;
#line 955 "MB02CU.f"
			daxpy_(&i__4, &d__1, &a1[j + (j + 1) * a1_dim1], lda1,
				 &b[(j + 1) * b_dim1 + 1], ldb);
#line 957 "MB02CU.f"
		    }
#line 958 "MB02CU.f"
		    cs[(j << 1) - 1] = c__;
#line 959 "MB02CU.f"
		    cs[j * 2] = s;
#line 960 "MB02CU.f"
/* L190: */
#line 960 "MB02CU.f"
		}

#line 962 "MB02CU.f"
/* L200: */
#line 962 "MB02CU.f"
	    }

#line 964 "MB02CU.f"
	} else {
#line 965 "MB02CU.f"
	    i__ = 1;
#line 966 "MB02CU.f"
	}

/*        Unblocked version for the last or only block. */

#line 970 "MB02CU.f"
	i__1 = *k;
#line 970 "MB02CU.f"
	for (j = i__; j <= i__1; ++j) {
#line 971 "MB02CU.f"
	    if (*q > 1) {
#line 972 "MB02CU.f"
		dlarfg_(q, &b[j * b_dim1 + 1], &b[j * b_dim1 + 2], &c__1, &
			tau2);
#line 973 "MB02CU.f"
		alpha2 = b[j * b_dim1 + 1];
#line 974 "MB02CU.f"
		if (*k > j) {
#line 975 "MB02CU.f"
		    b[j * b_dim1 + 1] = 1.;
#line 976 "MB02CU.f"
		    i__2 = *k - j;
#line 976 "MB02CU.f"
		    dlarf_("Left", q, &i__2, &b[j * b_dim1 + 1], &c__1, &tau2,
			     &b[(j + 1) * b_dim1 + 1], ldb, &dwork[1], (
			    ftnlen)4);
#line 978 "MB02CU.f"
		}
#line 979 "MB02CU.f"
		b[j * b_dim1 + 1] = tau2;
#line 980 "MB02CU.f"
	    } else if (*q > 0) {
#line 981 "MB02CU.f"
		alpha2 = b[j * b_dim1 + 1];
#line 982 "MB02CU.f"
		b[j * b_dim1 + 1] = 0.;
#line 983 "MB02CU.f"
	    }
#line 984 "MB02CU.f"
	    if (*q > 0) {

/*              Create hyperbolic rotation. */

#line 988 "MB02CU.f"
		ma02fd_(&a1[j + j * a1_dim1], &alpha2, &c__, &s, &ierr);
#line 989 "MB02CU.f"
		if (ierr != 0) {

/*                 Error return:  The matrix is not positive definite. */

#line 993 "MB02CU.f"
		    *info = 1;
#line 994 "MB02CU.f"
		    return 0;
#line 995 "MB02CU.f"
		}

/*              Apply hyperbolic rotation. */

#line 999 "MB02CU.f"
		if (*k > j) {
#line 1000 "MB02CU.f"
		    i__2 = *k - j;
#line 1000 "MB02CU.f"
		    d__1 = 1. / c__;
#line 1000 "MB02CU.f"
		    dscal_(&i__2, &d__1, &a1[j + (j + 1) * a1_dim1], lda1);
#line 1001 "MB02CU.f"
		    i__2 = *k - j;
#line 1001 "MB02CU.f"
		    d__1 = -s / c__;
#line 1001 "MB02CU.f"
		    daxpy_(&i__2, &d__1, &b[(j + 1) * b_dim1 + 1], ldb, &a1[j 
			    + (j + 1) * a1_dim1], lda1);
#line 1003 "MB02CU.f"
		    i__2 = *k - j;
#line 1003 "MB02CU.f"
		    dscal_(&i__2, &c__, &b[(j + 1) * b_dim1 + 1], ldb);
#line 1004 "MB02CU.f"
		    i__2 = *k - j;
#line 1004 "MB02CU.f"
		    d__1 = -s;
#line 1004 "MB02CU.f"
		    daxpy_(&i__2, &d__1, &a1[j + (j + 1) * a1_dim1], lda1, &b[
			    (j + 1) * b_dim1 + 1], ldb);
#line 1006 "MB02CU.f"
		}
#line 1007 "MB02CU.f"
		cs[(j << 1) - 1] = c__;
#line 1008 "MB02CU.f"
		cs[j * 2] = s;
#line 1009 "MB02CU.f"
	    }
#line 1010 "MB02CU.f"
/* L210: */
#line 1010 "MB02CU.f"
	}

#line 1012 "MB02CU.f"
    }

/* *** Last line of MB02CU *** */
#line 1015 "MB02CU.f"
    return 0;
} /* mb02cu_ */

