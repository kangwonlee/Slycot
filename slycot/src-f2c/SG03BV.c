#line 1 "SG03BV.f"
/* SG03BV.f -- translated by f2c (version 20100827).
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

#line 1 "SG03BV.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static doublereal c_b19 = -1.;
static doublereal c_b21 = 0.;
static doublereal c_b26 = 1.;
static integer c__4 = 4;

/* Subroutine */ int sg03bv_(char *trans, integer *n, doublereal *a, integer *
	lda, doublereal *e, integer *lde, doublereal *b, integer *ldb, 
	doublereal *scale, doublereal *dwork, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, e_dim1, e_offset, i__1, i__2, 
	    i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c__;
    static integer i__, j;
    static doublereal s, x, z__, m1[4]	/* was [2][2] */, m2[4]	/* was [2][2] 
	    */;
    static integer kb, kh, kl;
    static doublereal ui[4]	/* was [2][2] */, tm[4]	/* was [2][2] */, eps;
    static integer wpt, ypt;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer ldws, info1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sg03bw_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), sg03bx_(char *, char *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    dcopy_(integer *, doublereal *, integer *, doublereal *, integer *
	    ), drotg_(doublereal *, doublereal *, doublereal *, doublereal *),
	     dtrmm_(char *, char *, char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen, ftnlen);
    static integer uiipt;
    static doublereal scale1, delta1;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static doublereal bignum, smlnum;
    static logical notrns;


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

/*     To compute the Cholesky factor U of the matrix X, X = U**T * U or */
/*     X = U * U**T, which is the solution of the generalized c-stable */
/*     continuous-time Lyapunov equation */

/*         T            T                  2    T */
/*        A  * X * E + E  * X * A = - SCALE  * B  * B,                (1) */

/*     or the transposed equation */

/*                 T            T          2        T */
/*        A * X * E  + E * X * A  = - SCALE  * B * B ,                (2) */

/*     respectively, where A, E, B, and U are real N-by-N matrices. The */
/*     Cholesky factor U of the solution is computed without first */
/*     finding X. The pencil A - lambda * E must be in generalized Schur */
/*     form ( A upper quasitriangular, E upper triangular ). Moreover, it */
/*     must be c-stable, i.e. its eigenvalues must have negative real */
/*     parts. B must be an upper triangular matrix with non-negative */
/*     entries on its main diagonal. */

/*     The resulting matrix U is upper triangular. The entries on its */
/*     main diagonal are non-negative. SCALE is an output scale factor */
/*     set to avoid overflow in U. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANS   CHARACTER*1 */
/*             Specifies whether equation (1) or equation (2) is to be */
/*             solved: */
/*             = 'N':  Solve equation (1); */
/*             = 'T':  Solve equation (2). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N upper Hessenberg part of this array */
/*             must contain the quasitriangular matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     E       (input) DOUBLE PRECISION array, dimension (LDE,N) */
/*             The leading N-by-N upper triangular part of this array */
/*             must contain the matrix E. */

/*     LDE     INTEGER */
/*             The leading dimension of the array E.  LDE >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading N-by-N upper triangular part of this */
/*             array must contain the matrix B. */
/*             On exit, the leading N-by-N upper triangular part of this */
/*             array contains the solution matrix U. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor set to avoid overflow in U. */
/*             0 < SCALE <= 1. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (6*N-6) */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the generalized Sylvester equation to be solved in */
/*                   step II (see METHOD) is (nearly) singular to working */
/*                   precision;  perturbed values were used to solve the */
/*                   equation (but the matrices A and E are unchanged); */
/*             = 2:  the generalized Schur form of the pencil */
/*                   A - lambda * E contains a 2-by-2 main diagonal block */
/*                   whose eigenvalues are not a pair of conjugate */
/*                   complex numbers; */
/*             = 3:  the pencil A - lambda * E is not stable, i.e. there */
/*                   is an eigenvalue without a negative real part. */

/*     METHOD */

/*     The method [2] used by the routine is an extension of Hammarling's */
/*     algorithm [1] to generalized Lyapunov equations. */

/*     We present the method for solving equation (1). Equation (2) can */
/*     be treated in a similar fashion. For simplicity, assume SCALE = 1. */

/*     The matrix A is an upper quasitriangular matrix, i.e. it is a */
/*     block triangular matrix with square blocks on the main diagonal */
/*     and the block order at most 2. We use the following partitioning */
/*     for the matrices A, E, B and the solution matrix U */

/*               ( A11   A12 )        ( E11   E12 ) */
/*           A = (           ),   E = (           ), */
/*               (   0   A22 )        (   0   E22 ) */

/*               ( B11   B12 )        ( U11   U12 ) */
/*           B = (           ),   U = (           ).                  (3) */
/*               (   0   B22 )        (   0   U22 ) */

/*     The size of the (1,1)-blocks is 1-by-1 (iff A(2,1) = 0.0) or */
/*     2-by-2. */

/*     We compute U11 and U12**T in three steps. */

/*     Step I: */

/*        From (1) and (3) we get the 1-by-1 or 2-by-2 equation */

/*                T      T                  T      T */
/*             A11  * U11  * U11 * E11 + E11  * U11  * U11 * A11 */

/*                    T */
/*             = - B11  * B11. */

/*        For brevity, details are omitted here. The technique for */
/*        computing U11 is similar to those applied to standard Lyapunov */
/*        equations in Hammarling's algorithm ([1], section 6). */

/*        Furthermore, the auxiliary matrices M1 and M2 defined as */
/*        follows */

/*                               -1      -1 */
/*           M1 = U11 * A11 * E11   * U11 */

/*                         -1      -1 */
/*           M2 = B11 * E11   * U11 */

/*        are computed in a numerically reliable way. */

/*     Step II: */

/*        We solve for U12**T the generalized Sylvester equation */

/*              T      T      T      T */
/*           A22  * U12  + E22  * U12  * M1 */

/*                  T           T      T      T      T */
/*           = - B12  * M2 - A12  * U11  - E12  * U11  * M1. */

/*     Step III: */

/*        One can show that */

/*              T      T                  T      T */
/*           A22  * U22  * U22 * E22 + E22  * U22  * U22 * A22  = */

/*                T              T */
/*           - B22  * B22 - y * y                                     (4) */

/*        holds, where y is defined as follows */

/*                  T      T      T      T */
/*           w = E12  * U11  + E22  * U12 */
/*                  T         T */
/*           y = B12  - w * M2 . */

/*        If B22_tilde is the square triangular matrix arising from the */
/*        QR-factorization */

/*               ( B22_tilde )     ( B22  ) */
/*           Q * (           )  =  (      ), */
/*               (     0     )     ( y**T ) */

/*        then */

/*                T              T                T */
/*           - B22  * B22 - y * y   =  - B22_tilde  * B22_tilde. */

/*        Replacing the right hand side in (4) by the term */
/*        - B22_tilde**T * B22_tilde leads to a generalized Lyapunov */
/*        equation of lower dimension compared to (1). */

/*     The solution U of the equation (1) can be obtained by recursive */
/*     application of the steps I to III. */

/*     REFERENCES */

/*     [1] Hammarling, S.J. */
/*         Numerical solution of the stable, non-negative definite */
/*         Lyapunov equation. */
/*         IMA J. Num. Anal., 2, pp. 303-323, 1982. */

/*     [2] Penzl, T. */
/*         Numerical solution of generalized Lyapunov equations. */
/*         Advances in Comp. Math., vol. 8, pp. 33-48, 1998. */

/*     NUMERICAL ASPECTS */

/*     The routine requires 2*N**3 flops. Note that we count a single */
/*     floating point arithmetic operation as one flop. */

/*     FURTHER COMMENTS */

/*     The Lyapunov equation may be very ill-conditioned. In particular, */
/*     if the pencil A - lambda * E has a pair of almost degenerate */
/*     eigenvalues, then the Lyapunov equation will be ill-conditioned. */
/*     Perturbed values were used to solve the equation. */
/*     A condition estimate can be obtained from the routine SG03AD. */
/*     When setting the error indicator INFO, the routine does not test */
/*     for near instability in the equation but only for exact */
/*     instability. */

/*     CONTRIBUTOR */

/*     T. Penzl, Technical University Chemnitz, Germany, Aug. 1998. */

/*     REVISIONS */

/*     Sep. 1998 (V. Sima). */

/*     KEYWORDS */

/*     Lyapunov equation */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Decode input parameter. */

#line 278 "SG03BV.f"
    /* Parameter adjustments */
#line 278 "SG03BV.f"
    a_dim1 = *lda;
#line 278 "SG03BV.f"
    a_offset = 1 + a_dim1;
#line 278 "SG03BV.f"
    a -= a_offset;
#line 278 "SG03BV.f"
    e_dim1 = *lde;
#line 278 "SG03BV.f"
    e_offset = 1 + e_dim1;
#line 278 "SG03BV.f"
    e -= e_offset;
#line 278 "SG03BV.f"
    b_dim1 = *ldb;
#line 278 "SG03BV.f"
    b_offset = 1 + b_dim1;
#line 278 "SG03BV.f"
    b -= b_offset;
#line 278 "SG03BV.f"
    --dwork;
#line 278 "SG03BV.f"

#line 278 "SG03BV.f"
    /* Function Body */
#line 278 "SG03BV.f"
    notrns = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

#line 282 "SG03BV.f"
    if (! (notrns || lsame_(trans, "T", (ftnlen)1, (ftnlen)1))) {
#line 283 "SG03BV.f"
	*info = -1;
#line 284 "SG03BV.f"
    } else if (*n < 0) {
#line 285 "SG03BV.f"
	*info = -2;
#line 286 "SG03BV.f"
    } else if (*lda < max(1,*n)) {
#line 287 "SG03BV.f"
	*info = -4;
#line 288 "SG03BV.f"
    } else if (*lde < max(1,*n)) {
#line 289 "SG03BV.f"
	*info = -6;
#line 290 "SG03BV.f"
    } else if (*ldb < max(1,*n)) {
#line 291 "SG03BV.f"
	*info = -8;
#line 292 "SG03BV.f"
    } else {
#line 293 "SG03BV.f"
	*info = 0;
#line 294 "SG03BV.f"
    }
#line 295 "SG03BV.f"
    if (*info != 0) {
#line 296 "SG03BV.f"
	i__1 = -(*info);
#line 296 "SG03BV.f"
	xerbla_("SG03BV", &i__1, (ftnlen)6);
#line 297 "SG03BV.f"
	return 0;
#line 298 "SG03BV.f"
    }

#line 300 "SG03BV.f"
    *scale = 1.;

/*     Quick return if possible. */

#line 304 "SG03BV.f"
    if (*n == 0) {
#line 304 "SG03BV.f"
	return 0;
#line 304 "SG03BV.f"
    }

/*     Set constants to control overflow. */

#line 309 "SG03BV.f"
    eps = dlamch_("P", (ftnlen)1);
#line 310 "SG03BV.f"
    smlnum = dlamch_("S", (ftnlen)1) / eps;
#line 311 "SG03BV.f"
    bignum = 1. / smlnum;
#line 312 "SG03BV.f"
    dlabad_(&smlnum, &bignum);

/*     Set work space pointers and leading dimension of matrices in */
/*     work space. */

#line 317 "SG03BV.f"
    uiipt = 1;
#line 318 "SG03BV.f"
    wpt = (*n << 1) - 1;
#line 319 "SG03BV.f"
    ypt = (*n << 2) - 3;
#line 320 "SG03BV.f"
    ldws = *n - 1;

#line 322 "SG03BV.f"
    if (notrns) {

/*        Solve equation (1). */

/*        Main Loop. Compute block row U(KL:KH,KL:N). KB denotes the */
/*        number of rows in this block row. */

#line 329 "SG03BV.f"
	kh = 0;
/*        WHILE ( KH .LT. N ) DO */
#line 331 "SG03BV.f"
L20:
#line 331 "SG03BV.f"
	if (kh < *n) {
#line 332 "SG03BV.f"
	    kl = kh + 1;
#line 333 "SG03BV.f"
	    if (kl == *n) {
#line 334 "SG03BV.f"
		kh = *n;
#line 335 "SG03BV.f"
		kb = 1;
#line 336 "SG03BV.f"
	    } else {
#line 337 "SG03BV.f"
		if (a[kl + 1 + kl * a_dim1] == 0.) {
#line 338 "SG03BV.f"
		    kh = kl;
#line 339 "SG03BV.f"
		    kb = 1;
#line 340 "SG03BV.f"
		} else {
#line 341 "SG03BV.f"
		    kh = kl + 1;
#line 342 "SG03BV.f"
		    kb = 2;
#line 343 "SG03BV.f"
		}
#line 344 "SG03BV.f"
	    }

/*           STEP I: Compute block U(KL:KH,KL:KH) and the auxiliary */
/*                   matrices M1 and M2. (For the moment the result */
/*                   U(KL:KH,KL:KH) is stored in UI). */

#line 350 "SG03BV.f"
	    if (kb == 1) {
#line 351 "SG03BV.f"
		delta1 = a[kl + kl * a_dim1] * -2. * e[kl + kl * e_dim1];
#line 352 "SG03BV.f"
		if (delta1 <= 0.) {
#line 353 "SG03BV.f"
		    *info = 3;
#line 354 "SG03BV.f"
		    return 0;
#line 355 "SG03BV.f"
		}
#line 356 "SG03BV.f"
		delta1 = sqrt(delta1);
#line 357 "SG03BV.f"
		z__ = (d__1 = b[kl + kl * b_dim1], abs(d__1)) * 2. * smlnum;
#line 358 "SG03BV.f"
		if (z__ > delta1) {
#line 359 "SG03BV.f"
		    scale1 = delta1 / z__;
#line 360 "SG03BV.f"
		    *scale = scale1 * *scale;
#line 361 "SG03BV.f"
		    i__1 = *n;
#line 361 "SG03BV.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 362 "SG03BV.f"
			dscal_(&i__, &scale1, &b[i__ * b_dim1 + 1], &c__1);
#line 363 "SG03BV.f"
/* L40: */
#line 363 "SG03BV.f"
		    }
#line 364 "SG03BV.f"
		}
#line 365 "SG03BV.f"
		ui[0] = b[kl + kl * b_dim1] / delta1;
#line 366 "SG03BV.f"
		m1[0] = a[kl + kl * a_dim1] / e[kl + kl * e_dim1];
#line 367 "SG03BV.f"
		m2[0] = delta1 / e[kl + kl * e_dim1];
#line 368 "SG03BV.f"
	    } else {

/*              If a pair of complex conjugate eigenvalues occurs, apply */
/*              (complex) Hammarling algorithm for the 2-by-2 problem. */

#line 373 "SG03BV.f"
		sg03bx_("C", "N", &a[kl + kl * a_dim1], lda, &e[kl + kl * 
			e_dim1], lde, &b[kl + kl * b_dim1], ldb, ui, &c__2, &
			scale1, m1, &c__2, m2, &c__2, &info1, (ftnlen)1, (
			ftnlen)1);
#line 376 "SG03BV.f"
		if (info1 != 0) {
#line 377 "SG03BV.f"
		    *info = info1;
#line 378 "SG03BV.f"
		    return 0;
#line 379 "SG03BV.f"
		}
#line 380 "SG03BV.f"
		if (scale1 != 1.) {
#line 381 "SG03BV.f"
		    *scale = scale1 * *scale;
#line 382 "SG03BV.f"
		    i__1 = *n;
#line 382 "SG03BV.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 383 "SG03BV.f"
			dscal_(&i__, &scale1, &b[i__ * b_dim1 + 1], &c__1);
#line 384 "SG03BV.f"
/* L60: */
#line 384 "SG03BV.f"
		    }
#line 385 "SG03BV.f"
		}
#line 386 "SG03BV.f"
	    }

#line 388 "SG03BV.f"
	    if (kh < *n) {

/*              STEP II: Compute U(KL:KH,KH+1:N) by solving a generalized */
/*                       Sylvester equation. (For the moment the result */
/*                       U(KL:KH,KH+1:N) is stored in the workspace.) */

/*              Form right hand side of the Sylvester equation. */

#line 396 "SG03BV.f"
		i__1 = *n - kh;
#line 396 "SG03BV.f"
		dgemm_("T", "N", &i__1, &kb, &kb, &c_b19, &b[kl + (kh + 1) * 
			b_dim1], ldb, m2, &c__2, &c_b21, &dwork[uiipt], &ldws,
			 (ftnlen)1, (ftnlen)1);
#line 398 "SG03BV.f"
		i__1 = *n - kh;
#line 398 "SG03BV.f"
		dgemm_("T", "T", &i__1, &kb, &kb, &c_b19, &a[kl + (kh + 1) * 
			a_dim1], lda, ui, &c__2, &c_b26, &dwork[uiipt], &ldws,
			 (ftnlen)1, (ftnlen)1);
#line 400 "SG03BV.f"
		dgemm_("T", "N", &kb, &kb, &kb, &c_b26, ui, &c__2, m1, &c__2, 
			&c_b21, tm, &c__2, (ftnlen)1, (ftnlen)1);
#line 402 "SG03BV.f"
		i__1 = *n - kh;
#line 402 "SG03BV.f"
		dgemm_("T", "N", &i__1, &kb, &kb, &c_b19, &e[kl + (kh + 1) * 
			e_dim1], lde, tm, &c__2, &c_b26, &dwork[uiipt], &ldws,
			 (ftnlen)1, (ftnlen)1);

/*              Solve generalized Sylvester equation. */

#line 407 "SG03BV.f"
		dlaset_("A", &kb, &kb, &c_b21, &c_b26, tm, &c__2, (ftnlen)1);
#line 408 "SG03BV.f"
		i__1 = *n - kh;
#line 408 "SG03BV.f"
		sg03bw_("N", &i__1, &kb, &a[kh + 1 + (kh + 1) * a_dim1], lda, 
			tm, &c__2, &e[kh + 1 + (kh + 1) * e_dim1], lde, m1, &
			c__2, &dwork[uiipt], &ldws, &scale1, &info1, (ftnlen)
			1);
#line 411 "SG03BV.f"
		if (info1 != 0) {
#line 411 "SG03BV.f"
		    *info = 1;
#line 411 "SG03BV.f"
		}
#line 413 "SG03BV.f"
		if (scale1 != 1.) {
#line 414 "SG03BV.f"
		    *scale = scale1 * *scale;
#line 415 "SG03BV.f"
		    i__1 = *n;
#line 415 "SG03BV.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 416 "SG03BV.f"
			dscal_(&i__, &scale1, &b[i__ * b_dim1 + 1], &c__1);
#line 417 "SG03BV.f"
/* L80: */
#line 417 "SG03BV.f"
		    }
#line 418 "SG03BV.f"
		    dscal_(&c__4, &scale1, ui, &c__1);
#line 419 "SG03BV.f"
		}

/*              STEP III: Form the right hand side matrix */
/*                        B(KH+1:N,KH+1:N) of the (smaller) Lyapunov */
/*                        equation to be solved during the next pass of */
/*                        the main loop. */

/*              Compute auxiliary vectors (or matrices) W and Y. */

#line 428 "SG03BV.f"
		i__1 = *n - kh;
#line 428 "SG03BV.f"
		dlacpy_("A", &i__1, &kb, &dwork[uiipt], &ldws, &dwork[wpt], &
			ldws, (ftnlen)1);
#line 430 "SG03BV.f"
		i__1 = *n - kh;
#line 430 "SG03BV.f"
		dtrmm_("L", "U", "T", "N", &i__1, &kb, &c_b26, &e[kh + 1 + (
			kh + 1) * e_dim1], lde, &dwork[wpt], &ldws, (ftnlen)1,
			 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 432 "SG03BV.f"
		i__1 = *n - kh;
#line 432 "SG03BV.f"
		dgemm_("T", "T", &i__1, &kb, &kb, &c_b26, &e[kl + (kh + 1) * 
			e_dim1], lde, ui, &c__2, &c_b26, &dwork[wpt], &ldws, (
			ftnlen)1, (ftnlen)1);
#line 434 "SG03BV.f"
		i__1 = kh;
#line 434 "SG03BV.f"
		for (i__ = kl; i__ <= i__1; ++i__) {
#line 435 "SG03BV.f"
		    i__2 = *n - kh;
#line 435 "SG03BV.f"
		    dcopy_(&i__2, &b[i__ + (kh + 1) * b_dim1], ldb, &dwork[
			    ypt + ldws * (i__ - kl)], &c__1);
#line 437 "SG03BV.f"
/* L100: */
#line 437 "SG03BV.f"
		}
#line 438 "SG03BV.f"
		i__1 = *n - kh;
#line 438 "SG03BV.f"
		dgemm_("N", "T", &i__1, &kb, &kb, &c_b19, &dwork[wpt], &ldws, 
			m2, &c__2, &c_b26, &dwork[ypt], &ldws, (ftnlen)1, (
			ftnlen)1);

/*              Overwrite B(KH+1:N,KH+1:N) with the triangular matrix */
/*              from the QR-factorization of the (N-KH+KB)-by-(N-KH) */
/*              matrix */

/*                          (  B(KH+1:N,KH+1:N)  ) */
/*                          (                    ) */
/*                          (       Y**T         ) . */

#line 449 "SG03BV.f"
		i__1 = kb;
#line 449 "SG03BV.f"
		for (j = 1; j <= i__1; ++j) {
#line 450 "SG03BV.f"
		    i__2 = *n - kh;
#line 450 "SG03BV.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 451 "SG03BV.f"
			x = b[kh + i__ + (kh + i__) * b_dim1];
#line 452 "SG03BV.f"
			z__ = dwork[ypt + i__ - 1 + (j - 1) * ldws];
#line 453 "SG03BV.f"
			drotg_(&x, &z__, &c__, &s);
#line 454 "SG03BV.f"
			i__3 = *n - kh - i__ + 1;
#line 454 "SG03BV.f"
			drot_(&i__3, &b[kh + i__ + (kh + i__) * b_dim1], ldb, 
				&dwork[ypt + i__ - 1 + (j - 1) * ldws], &c__1,
				 &c__, &s);
#line 456 "SG03BV.f"
/* L120: */
#line 456 "SG03BV.f"
		    }
#line 457 "SG03BV.f"
/* L140: */
#line 457 "SG03BV.f"
		}

/*              Make main diagonal elements of B(KH+1:N,KH+1:N) positive. */

#line 461 "SG03BV.f"
		i__1 = *n;
#line 461 "SG03BV.f"
		for (i__ = kh + 1; i__ <= i__1; ++i__) {
#line 462 "SG03BV.f"
		    if (b[i__ + i__ * b_dim1] < 0.) {
#line 462 "SG03BV.f"
			i__2 = *n - i__ + 1;
#line 462 "SG03BV.f"
			dscal_(&i__2, &c_b19, &b[i__ + i__ * b_dim1], ldb);
#line 462 "SG03BV.f"
		    }
#line 464 "SG03BV.f"
/* L160: */
#line 464 "SG03BV.f"
		}

/*              Overwrite right hand side with the part of the solution */
/*              computed in step II. */

#line 469 "SG03BV.f"
		i__1 = kh;
#line 469 "SG03BV.f"
		for (j = kl; j <= i__1; ++j) {
#line 470 "SG03BV.f"
		    i__2 = *n - kh;
#line 470 "SG03BV.f"
		    dcopy_(&i__2, &dwork[uiipt + (j - kl) * ldws], &c__1, &b[
			    j + (kh + 1) * b_dim1], ldb);
#line 472 "SG03BV.f"
/* L180: */
#line 472 "SG03BV.f"
		}
#line 473 "SG03BV.f"
	    }

/*           Overwrite right hand side with the part of the solution */
/*           computed in step I. */

#line 478 "SG03BV.f"
	    dlacpy_("U", &kb, &kb, ui, &c__2, &b[kl + kl * b_dim1], ldb, (
		    ftnlen)1);

#line 480 "SG03BV.f"
	    goto L20;
#line 481 "SG03BV.f"
	}
/*        END WHILE 20 */

#line 484 "SG03BV.f"
    } else {

/*        Solve equation (2). */

/*        Main Loop. Compute block column U(1:KH,KL:KH). KB denotes the */
/*        number of columns in this block column. */

#line 491 "SG03BV.f"
	kl = *n + 1;
/*        WHILE ( KL .GT. 1 ) DO */
#line 493 "SG03BV.f"
L200:
#line 493 "SG03BV.f"
	if (kl > 1) {
#line 494 "SG03BV.f"
	    kh = kl - 1;
#line 495 "SG03BV.f"
	    if (kh == 1) {
#line 496 "SG03BV.f"
		kl = 1;
#line 497 "SG03BV.f"
		kb = 1;
#line 498 "SG03BV.f"
	    } else {
#line 499 "SG03BV.f"
		if (a[kh + (kh - 1) * a_dim1] == 0.) {
#line 500 "SG03BV.f"
		    kl = kh;
#line 501 "SG03BV.f"
		    kb = 1;
#line 502 "SG03BV.f"
		} else {
#line 503 "SG03BV.f"
		    kl = kh - 1;
#line 504 "SG03BV.f"
		    kb = 2;
#line 505 "SG03BV.f"
		}
#line 506 "SG03BV.f"
	    }

/*           STEP I: Compute block U(KL:KH,KL:KH) and the auxiliary */
/*                   matrices M1 and M2. (For the moment the result */
/*                   U(KL:KH,KL:KH) is stored in UI). */

#line 512 "SG03BV.f"
	    if (kb == 1) {
#line 513 "SG03BV.f"
		delta1 = a[kl + kl * a_dim1] * -2. * e[kl + kl * e_dim1];
#line 514 "SG03BV.f"
		if (delta1 <= 0.) {
#line 515 "SG03BV.f"
		    *info = 3;
#line 516 "SG03BV.f"
		    return 0;
#line 517 "SG03BV.f"
		}
#line 518 "SG03BV.f"
		delta1 = sqrt(delta1);
#line 519 "SG03BV.f"
		z__ = (d__1 = b[kl + kl * b_dim1], abs(d__1)) * 2. * smlnum;
#line 520 "SG03BV.f"
		if (z__ > delta1) {
#line 521 "SG03BV.f"
		    scale1 = delta1 / z__;
#line 522 "SG03BV.f"
		    *scale = scale1 * *scale;
#line 523 "SG03BV.f"
		    i__1 = *n;
#line 523 "SG03BV.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 524 "SG03BV.f"
			dscal_(&i__, &scale1, &b[i__ * b_dim1 + 1], &c__1);
#line 525 "SG03BV.f"
/* L220: */
#line 525 "SG03BV.f"
		    }
#line 526 "SG03BV.f"
		}
#line 527 "SG03BV.f"
		ui[0] = b[kl + kl * b_dim1] / delta1;
#line 528 "SG03BV.f"
		m1[0] = a[kl + kl * a_dim1] / e[kl + kl * e_dim1];
#line 529 "SG03BV.f"
		m2[0] = delta1 / e[kl + kl * e_dim1];
#line 530 "SG03BV.f"
	    } else {

/*              If a pair of complex conjugate eigenvalues occurs, apply */
/*              (complex) Hammarling algorithm for the 2-by-2 problem. */

#line 535 "SG03BV.f"
		sg03bx_("C", "T", &a[kl + kl * a_dim1], lda, &e[kl + kl * 
			e_dim1], lde, &b[kl + kl * b_dim1], ldb, ui, &c__2, &
			scale1, m1, &c__2, m2, &c__2, &info1, (ftnlen)1, (
			ftnlen)1);
#line 538 "SG03BV.f"
		if (info1 != 0) {
#line 539 "SG03BV.f"
		    *info = info1;
#line 540 "SG03BV.f"
		    return 0;
#line 541 "SG03BV.f"
		}
#line 542 "SG03BV.f"
		if (scale1 != 1.) {
#line 543 "SG03BV.f"
		    *scale = scale1 * *scale;
#line 544 "SG03BV.f"
		    i__1 = *n;
#line 544 "SG03BV.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 545 "SG03BV.f"
			dscal_(&i__, &scale1, &b[i__ * b_dim1 + 1], &c__1);
#line 546 "SG03BV.f"
/* L240: */
#line 546 "SG03BV.f"
		    }
#line 547 "SG03BV.f"
		}
#line 548 "SG03BV.f"
	    }

#line 550 "SG03BV.f"
	    if (kl > 1) {

/*              STEP II: Compute U(1:KL-1,KL:KH) by solving a generalized */
/*                       Sylvester equation. (For the moment the result */
/*                       U(1:KL-1,KL:KH) is stored in the workspace.) */

/*              Form right hand side of the Sylvester equation. */

#line 558 "SG03BV.f"
		i__1 = kl - 1;
#line 558 "SG03BV.f"
		dgemm_("N", "T", &i__1, &kb, &kb, &c_b19, &b[kl * b_dim1 + 1],
			 ldb, m2, &c__2, &c_b21, &dwork[uiipt], &ldws, (
			ftnlen)1, (ftnlen)1);
#line 560 "SG03BV.f"
		i__1 = kl - 1;
#line 560 "SG03BV.f"
		dgemm_("N", "N", &i__1, &kb, &kb, &c_b19, &a[kl * a_dim1 + 1],
			 lda, ui, &c__2, &c_b26, &dwork[uiipt], &ldws, (
			ftnlen)1, (ftnlen)1);
#line 562 "SG03BV.f"
		dgemm_("N", "T", &kb, &kb, &kb, &c_b26, ui, &c__2, m1, &c__2, 
			&c_b21, tm, &c__2, (ftnlen)1, (ftnlen)1);
#line 564 "SG03BV.f"
		i__1 = kl - 1;
#line 564 "SG03BV.f"
		dgemm_("N", "N", &i__1, &kb, &kb, &c_b19, &e[kl * e_dim1 + 1],
			 lde, tm, &c__2, &c_b26, &dwork[uiipt], &ldws, (
			ftnlen)1, (ftnlen)1);

/*              Solve generalized Sylvester equation. */

#line 569 "SG03BV.f"
		dlaset_("A", &kb, &kb, &c_b21, &c_b26, tm, &c__2, (ftnlen)1);
#line 570 "SG03BV.f"
		i__1 = kl - 1;
#line 570 "SG03BV.f"
		sg03bw_("T", &i__1, &kb, &a[a_offset], lda, tm, &c__2, &e[
			e_offset], lde, m1, &c__2, &dwork[uiipt], &ldws, &
			scale1, &info1, (ftnlen)1);
#line 572 "SG03BV.f"
		if (info1 != 0) {
#line 572 "SG03BV.f"
		    *info = 1;
#line 572 "SG03BV.f"
		}
#line 574 "SG03BV.f"
		if (scale1 != 1.) {
#line 575 "SG03BV.f"
		    *scale = scale1 * *scale;
#line 576 "SG03BV.f"
		    i__1 = *n;
#line 576 "SG03BV.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 577 "SG03BV.f"
			dscal_(&i__, &scale1, &b[i__ * b_dim1 + 1], &c__1);
#line 578 "SG03BV.f"
/* L260: */
#line 578 "SG03BV.f"
		    }
#line 579 "SG03BV.f"
		    dscal_(&c__4, &scale1, ui, &c__1);
#line 580 "SG03BV.f"
		}

/*              STEP III: Form the right hand side matrix */
/*                        B(1:KL-1,1:KL-1) of the (smaller) Lyapunov */
/*                        equation to be solved during the next pass of */
/*                        the main loop. */

/*              Compute auxiliary vectors (or matrices) W and Y. */

#line 589 "SG03BV.f"
		i__1 = kl - 1;
#line 589 "SG03BV.f"
		dlacpy_("A", &i__1, &kb, &dwork[uiipt], &ldws, &dwork[wpt], &
			ldws, (ftnlen)1);
#line 591 "SG03BV.f"
		i__1 = kl - 1;
#line 591 "SG03BV.f"
		dtrmm_("L", "U", "N", "N", &i__1, &kb, &c_b26, &e[e_dim1 + 1],
			 lde, &dwork[wpt], &ldws, (ftnlen)1, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 593 "SG03BV.f"
		i__1 = kl - 1;
#line 593 "SG03BV.f"
		dgemm_("N", "N", &i__1, &kb, &kb, &c_b26, &e[kl * e_dim1 + 1],
			 lde, ui, &c__2, &c_b26, &dwork[wpt], &ldws, (ftnlen)
			1, (ftnlen)1);
#line 595 "SG03BV.f"
		i__1 = kl - 1;
#line 595 "SG03BV.f"
		dlacpy_("A", &i__1, &kb, &b[kl * b_dim1 + 1], ldb, &dwork[ypt]
			, &ldws, (ftnlen)1);
#line 597 "SG03BV.f"
		i__1 = kl - 1;
#line 597 "SG03BV.f"
		dgemm_("N", "N", &i__1, &kb, &kb, &c_b19, &dwork[wpt], &ldws, 
			m2, &c__2, &c_b26, &dwork[ypt], &ldws, (ftnlen)1, (
			ftnlen)1);

/*              Overwrite B(1:KL-1,1:KL-1) with the triangular matrix */
/*              from the RQ-factorization of the (KL-1)-by-KH matrix */

/*                          (                        ) */
/*                          (  B(1:KL-1,1:KL-1)   Y  ) */
/*                          (                        ). */

#line 607 "SG03BV.f"
		i__1 = kb;
#line 607 "SG03BV.f"
		for (j = 1; j <= i__1; ++j) {
#line 608 "SG03BV.f"
		    for (i__ = kl - 1; i__ >= 1; --i__) {
#line 609 "SG03BV.f"
			x = b[i__ + i__ * b_dim1];
#line 610 "SG03BV.f"
			z__ = dwork[ypt + i__ - 1 + (j - 1) * ldws];
#line 611 "SG03BV.f"
			drotg_(&x, &z__, &c__, &s);
#line 612 "SG03BV.f"
			drot_(&i__, &b[i__ * b_dim1 + 1], &c__1, &dwork[ypt + 
				(j - 1) * ldws], &c__1, &c__, &s);
#line 614 "SG03BV.f"
/* L280: */
#line 614 "SG03BV.f"
		    }
#line 615 "SG03BV.f"
/* L300: */
#line 615 "SG03BV.f"
		}

/*              Make main diagonal elements of B(1:KL-1,1:KL-1) positive. */

#line 619 "SG03BV.f"
		i__1 = kl - 1;
#line 619 "SG03BV.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 620 "SG03BV.f"
		    if (b[i__ + i__ * b_dim1] < 0.) {
#line 620 "SG03BV.f"
			dscal_(&i__, &c_b19, &b[i__ * b_dim1 + 1], &c__1);
#line 620 "SG03BV.f"
		    }
#line 622 "SG03BV.f"
/* L320: */
#line 622 "SG03BV.f"
		}

/*              Overwrite right hand side with the part of the solution */
/*              computed in step II. */

#line 627 "SG03BV.f"
		i__1 = kl - 1;
#line 627 "SG03BV.f"
		dlacpy_("A", &i__1, &kb, &dwork[uiipt], &ldws, &b[kl * b_dim1 
			+ 1], ldb, (ftnlen)1);

#line 630 "SG03BV.f"
	    }

/*           Overwrite right hand side with the part of the solution */
/*           computed in step I. */

#line 635 "SG03BV.f"
	    dlacpy_("U", &kb, &kb, ui, &c__2, &b[kl + kl * b_dim1], ldb, (
		    ftnlen)1);

#line 637 "SG03BV.f"
	    goto L200;
#line 638 "SG03BV.f"
	}
/*        END WHILE 200 */

#line 641 "SG03BV.f"
    }

#line 643 "SG03BV.f"
    return 0;
/* *** Last line of SG03BV *** */
} /* sg03bv_ */

