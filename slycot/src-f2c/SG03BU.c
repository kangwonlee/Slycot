#line 1 "SG03BU.f"
/* SG03BU.f -- translated by f2c (version 20100827).
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

#line 1 "SG03BU.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static doublereal c_b19 = -1.;
static doublereal c_b21 = 0.;
static doublereal c_b24 = 1.;
static integer c__4 = 4;
static doublereal c_b77 = .5;
static doublereal c_b78 = 2.;
static integer c__32 = 32;

/* Subroutine */ int sg03bu_(char *trans, integer *n, doublereal *a, integer *
	lda, doublereal *e, integer *lde, doublereal *b, integer *ldb, 
	doublereal *scale, doublereal *dwork, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, e_dim1, e_offset, i__1, i__2, 
	    i__3, i__4, i__5, i__6, i__7;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c__;
    static integer i__, j, m;
    static doublereal s, x, z__, m1[4]	/* was [2][2] */, m2[4]	/* was [2][2] 
	    */, m3[16]	/* was [4][4] */;
    static integer kb, kh, kl;
    static doublereal ui[4]	/* was [2][2] */;
    static integer iw[24];
    static doublereal tm[4]	/* was [2][2] */, rw[32], m3c[16]	/* 
	    was [4][4] */, eps;
    static integer wpt, ypt;
    static doublereal m3ew[4];
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer ldws;
    static doublereal uflt;
    static integer info1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), sg03bw_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    sg03bx_(char *, char *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), dcopy_(integer *, doublereal *, 
	    integer *, doublereal *, integer *), drotg_(doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static integer uiipt;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen);
    static doublereal scale1, delta1;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static doublereal bignum, smlnum;
    extern /* Subroutine */ int dsyevx_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, ftnlen, ftnlen, ftnlen);
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
/*     X = U * U**T, which is the solution of the generalized d-stable */
/*     discrete-time Lyapunov equation */

/*         T            T                  2    T */
/*        A  * X * A - E  * X * E = - SCALE  * B  * B,                (1) */

/*     or the transposed equation */

/*                 T            T          2        T */
/*        A * X * A  - E * X * E  = - SCALE  * B * B ,                (2) */

/*     respectively, where A, E, B, and U are real N-by-N matrices. The */
/*     Cholesky factor U of the solution is computed without first */
/*     finding X. The pencil A - lambda * E must be in generalized Schur */
/*     form ( A upper quasitriangular, E upper triangular ). Moreover, it */
/*     must be d-stable, i.e. the moduli of its eigenvalues must be less */
/*     than one. B must be an upper triangular matrix with non-negative */
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
/*             = 3:  the pencil A - lambda * E is not d-stable, i.e. */
/*                   there are eigenvalues outside the open unit circle; */
/*             = 4:  the LAPACK routine DSYEVX utilized to factorize M3 */
/*                   failed to converge. This error is unlikely to occur. */

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

/*                T      T                   T      T */
/*             A11  * U11  * U11 * A11  - E11  * U11  * U11 * E11 */

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

/*              T      T           T      T */
/*           A22  * U12  * M1 - E22  * U12 */

/*                  T           T      T      T      T */
/*           = - B12  * M2 + E12  * U11  - A12  * U11  * M1. */

/*     Step III: */

/*        One can show that */

/*              T      T                  T      T */
/*           A22  * U22  * U22 * A22 - E22  * U22  * U22 * E22  = */

/*                T              T */
/*           - B22  * B22 - y * y                                     (4) */

/*        holds, where y is defined as follows */

/*                  T      T      T      T */
/*           w = A12  * U11  + A22  * U12 */

/*                    T */
/*           y = ( B12   w ) * M3EV, */

/*        where M3EV is a matrix which fulfils */

/*                ( I-M2*M2**T   -M2*M1**T )              T */
/*           M3 = (                        ) = M3EV * M3EV . */
/*                (  -M1*M2**T  I-M1*M1**T ) */

/*        M3 is positive semidefinite and its rank is equal to the size */
/*        of U11. Therefore, a matrix M3EV can be found by solving the */
/*        symmetric eigenvalue problem for M3 such that y consists of */
/*        either 1 or 2 rows. */

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
/*     if the pencil A - lambda * E has a pair of almost reciprocal */
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

#line 297 "SG03BU.f"
    /* Parameter adjustments */
#line 297 "SG03BU.f"
    a_dim1 = *lda;
#line 297 "SG03BU.f"
    a_offset = 1 + a_dim1;
#line 297 "SG03BU.f"
    a -= a_offset;
#line 297 "SG03BU.f"
    e_dim1 = *lde;
#line 297 "SG03BU.f"
    e_offset = 1 + e_dim1;
#line 297 "SG03BU.f"
    e -= e_offset;
#line 297 "SG03BU.f"
    b_dim1 = *ldb;
#line 297 "SG03BU.f"
    b_offset = 1 + b_dim1;
#line 297 "SG03BU.f"
    b -= b_offset;
#line 297 "SG03BU.f"
    --dwork;
#line 297 "SG03BU.f"

#line 297 "SG03BU.f"
    /* Function Body */
#line 297 "SG03BU.f"
    notrns = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

#line 301 "SG03BU.f"
    if (! (notrns || lsame_(trans, "T", (ftnlen)1, (ftnlen)1))) {
#line 302 "SG03BU.f"
	*info = -1;
#line 303 "SG03BU.f"
    } else if (*n < 0) {
#line 304 "SG03BU.f"
	*info = -2;
#line 305 "SG03BU.f"
    } else if (*lda < max(1,*n)) {
#line 306 "SG03BU.f"
	*info = -4;
#line 307 "SG03BU.f"
    } else if (*lde < max(1,*n)) {
#line 308 "SG03BU.f"
	*info = -6;
#line 309 "SG03BU.f"
    } else if (*ldb < max(1,*n)) {
#line 310 "SG03BU.f"
	*info = -8;
#line 311 "SG03BU.f"
    } else {
#line 312 "SG03BU.f"
	*info = 0;
#line 313 "SG03BU.f"
    }
#line 314 "SG03BU.f"
    if (*info != 0) {
#line 315 "SG03BU.f"
	i__1 = -(*info);
#line 315 "SG03BU.f"
	xerbla_("SG03BU", &i__1, (ftnlen)6);
#line 316 "SG03BU.f"
	return 0;
#line 317 "SG03BU.f"
    }

#line 319 "SG03BU.f"
    *scale = 1.;

/*     Quick return if possible. */

#line 323 "SG03BU.f"
    if (*n == 0) {
#line 323 "SG03BU.f"
	return 0;
#line 323 "SG03BU.f"
    }

/*     Set constants to control overflow. */

#line 328 "SG03BU.f"
    eps = dlamch_("P", (ftnlen)1);
#line 329 "SG03BU.f"
    uflt = dlamch_("S", (ftnlen)1);
#line 330 "SG03BU.f"
    smlnum = uflt / eps;
#line 331 "SG03BU.f"
    bignum = 1. / smlnum;
#line 332 "SG03BU.f"
    dlabad_(&smlnum, &bignum);

/*     Set work space pointers and leading dimension of matrices in */
/*     work space. */

#line 337 "SG03BU.f"
    uiipt = 1;
#line 338 "SG03BU.f"
    wpt = (*n << 1) - 1;
#line 339 "SG03BU.f"
    ypt = (*n << 2) - 3;
#line 340 "SG03BU.f"
    ldws = *n - 1;

#line 342 "SG03BU.f"
    if (notrns) {

/*        Solve equation (1). */

/*        Main Loop. Compute block row U(KL:KH,KL:N). KB denotes the */
/*        number of rows in this block row. */

#line 349 "SG03BU.f"
	kh = 0;
/*        WHILE ( KH .LT. N ) DO */
#line 351 "SG03BU.f"
L20:
#line 351 "SG03BU.f"
	if (kh < *n) {
#line 352 "SG03BU.f"
	    kl = kh + 1;
#line 353 "SG03BU.f"
	    if (kl == *n) {
#line 354 "SG03BU.f"
		kh = *n;
#line 355 "SG03BU.f"
		kb = 1;
#line 356 "SG03BU.f"
	    } else {
#line 357 "SG03BU.f"
		if (a[kl + 1 + kl * a_dim1] == 0.) {
#line 358 "SG03BU.f"
		    kh = kl;
#line 359 "SG03BU.f"
		    kb = 1;
#line 360 "SG03BU.f"
		} else {
#line 361 "SG03BU.f"
		    kh = kl + 1;
#line 362 "SG03BU.f"
		    kb = 2;
#line 363 "SG03BU.f"
		}
#line 364 "SG03BU.f"
	    }

/*           STEP I: Compute block U(KL:KH,KL:KH) and the auxiliary */
/*                   matrices M1 and M2. (For the moment the result */
/*                   U(KL:KH,KL:KH) is stored in UI). */

#line 370 "SG03BU.f"
	    if (kb == 1) {
/* Computing 2nd power */
#line 371 "SG03BU.f"
		d__1 = e[kl + kl * e_dim1];
/* Computing 2nd power */
#line 371 "SG03BU.f"
		d__2 = a[kl + kl * a_dim1];
#line 371 "SG03BU.f"
		delta1 = d__1 * d__1 - d__2 * d__2;
#line 372 "SG03BU.f"
		if (delta1 <= 0.) {
#line 373 "SG03BU.f"
		    *info = 3;
#line 374 "SG03BU.f"
		    return 0;
#line 375 "SG03BU.f"
		}
#line 376 "SG03BU.f"
		delta1 = sqrt(delta1);
#line 377 "SG03BU.f"
		z__ = (d__1 = b[kl + kl * b_dim1], abs(d__1)) * 2. * smlnum;
#line 378 "SG03BU.f"
		if (z__ > delta1) {
#line 379 "SG03BU.f"
		    scale1 = delta1 / z__;
#line 380 "SG03BU.f"
		    *scale = scale1 * *scale;
#line 381 "SG03BU.f"
		    i__1 = *n;
#line 381 "SG03BU.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 382 "SG03BU.f"
			dscal_(&i__, &scale1, &b[i__ * b_dim1 + 1], &c__1);
#line 383 "SG03BU.f"
/* L40: */
#line 383 "SG03BU.f"
		    }
#line 384 "SG03BU.f"
		}
#line 385 "SG03BU.f"
		ui[0] = b[kl + kl * b_dim1] / delta1;
#line 386 "SG03BU.f"
		m1[0] = a[kl + kl * a_dim1] / e[kl + kl * e_dim1];
#line 387 "SG03BU.f"
		m2[0] = delta1 / e[kl + kl * e_dim1];
#line 388 "SG03BU.f"
	    } else {

/*              If a pair of complex conjugate eigenvalues occurs, apply */
/*              (complex) Hammarling algorithm for the 2-by-2 problem. */

#line 393 "SG03BU.f"
		sg03bx_("D", "N", &a[kl + kl * a_dim1], lda, &e[kl + kl * 
			e_dim1], lde, &b[kl + kl * b_dim1], ldb, ui, &c__2, &
			scale1, m1, &c__2, m2, &c__2, &info1, (ftnlen)1, (
			ftnlen)1);
#line 396 "SG03BU.f"
		if (info1 != 0) {
#line 397 "SG03BU.f"
		    *info = info1;
#line 398 "SG03BU.f"
		    return 0;
#line 399 "SG03BU.f"
		}
#line 400 "SG03BU.f"
		if (scale1 != 1.) {
#line 401 "SG03BU.f"
		    *scale = scale1 * *scale;
#line 402 "SG03BU.f"
		    i__1 = *n;
#line 402 "SG03BU.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 403 "SG03BU.f"
			dscal_(&i__, &scale1, &b[i__ * b_dim1 + 1], &c__1);
#line 404 "SG03BU.f"
/* L60: */
#line 404 "SG03BU.f"
		    }
#line 405 "SG03BU.f"
		}
#line 406 "SG03BU.f"
	    }

#line 408 "SG03BU.f"
	    if (kh < *n) {

/*              STEP II: Compute U(KL:KH,KH+1:N) by solving a generalized */
/*                       Sylvester equation. (For the moment the result */
/*                       U(KL:KH,KH+1:N) is stored in the workspace.) */

/*              Form right hand side of the Sylvester equation. */

#line 416 "SG03BU.f"
		i__1 = *n - kh;
#line 416 "SG03BU.f"
		dgemm_("T", "N", &i__1, &kb, &kb, &c_b19, &b[kl + (kh + 1) * 
			b_dim1], ldb, m2, &c__2, &c_b21, &dwork[uiipt], &ldws,
			 (ftnlen)1, (ftnlen)1);
#line 418 "SG03BU.f"
		i__1 = *n - kh;
#line 418 "SG03BU.f"
		dgemm_("T", "T", &i__1, &kb, &kb, &c_b24, &e[kl + (kh + 1) * 
			e_dim1], lde, ui, &c__2, &c_b24, &dwork[uiipt], &ldws,
			 (ftnlen)1, (ftnlen)1);
#line 420 "SG03BU.f"
		dgemm_("T", "N", &kb, &kb, &kb, &c_b24, ui, &c__2, m1, &c__2, 
			&c_b21, tm, &c__2, (ftnlen)1, (ftnlen)1);
#line 422 "SG03BU.f"
		i__1 = *n - kh;
#line 422 "SG03BU.f"
		dgemm_("T", "N", &i__1, &kb, &kb, &c_b19, &a[kl + (kh + 1) * 
			a_dim1], lda, tm, &c__2, &c_b24, &dwork[uiipt], &ldws,
			 (ftnlen)1, (ftnlen)1);

/*              Solve generalized Sylvester equation. */

#line 427 "SG03BU.f"
		dlaset_("A", &kb, &kb, &c_b21, &c_b19, tm, &c__2, (ftnlen)1);
#line 428 "SG03BU.f"
		i__1 = *n - kh;
#line 428 "SG03BU.f"
		sg03bw_("N", &i__1, &kb, &a[kh + 1 + (kh + 1) * a_dim1], lda, 
			m1, &c__2, &e[kh + 1 + (kh + 1) * e_dim1], lde, tm, &
			c__2, &dwork[uiipt], &ldws, &scale1, &info1, (ftnlen)
			1);
#line 431 "SG03BU.f"
		if (info1 != 0) {
#line 431 "SG03BU.f"
		    *info = 1;
#line 431 "SG03BU.f"
		}
#line 433 "SG03BU.f"
		if (scale1 != 1.) {
#line 434 "SG03BU.f"
		    *scale = scale1 * *scale;
#line 435 "SG03BU.f"
		    i__1 = *n;
#line 435 "SG03BU.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 436 "SG03BU.f"
			dscal_(&i__, &scale1, &b[i__ * b_dim1 + 1], &c__1);
#line 437 "SG03BU.f"
/* L80: */
#line 437 "SG03BU.f"
		    }
#line 438 "SG03BU.f"
		    dscal_(&c__4, &scale1, ui, &c__1);
#line 439 "SG03BU.f"
		}

/*              STEP III: Form the right hand side matrix */
/*                        B(KH+1:N,KH+1:N) of the (smaller) Lyapunov */
/*                        equation to be solved during the next pass of */
/*                        the main loop. */

/*              Compute auxiliary matrices M3 and Y. The factorization */
/*              M3 = M3C * M3C**T is found by solving the symmetric */
/*              eigenvalue problem. */

#line 450 "SG03BU.f"
		i__1 = kb << 1;
#line 450 "SG03BU.f"
		i__2 = kb << 1;
#line 450 "SG03BU.f"
		dlaset_("U", &i__1, &i__2, &c_b21, &c_b24, m3, &c__4, (ftnlen)
			1);
#line 451 "SG03BU.f"
		dsyrk_("U", "N", &kb, &kb, &c_b19, m2, &c__2, &c_b24, m3, &
			c__4, (ftnlen)1, (ftnlen)1);
#line 452 "SG03BU.f"
		dgemm_("N", "T", &kb, &kb, &kb, &c_b19, m2, &c__2, m1, &c__2, 
			&c_b21, &m3[(kb + 1 << 2) - 4], &c__4, (ftnlen)1, (
			ftnlen)1);
#line 454 "SG03BU.f"
		dsyrk_("U", "N", &kb, &kb, &c_b19, m1, &c__2, &c_b24, &m3[kb 
			+ 1 + (kb + 1 << 2) - 5], &c__4, (ftnlen)1, (ftnlen)1)
			;
#line 456 "SG03BU.f"
		i__1 = kb << 1;
#line 456 "SG03BU.f"
		d__1 = uflt * 2.;
#line 456 "SG03BU.f"
		dsyevx_("V", "V", "U", &i__1, m3, &c__4, &c_b77, &c_b78, &
			c__1, &c__4, &d__1, &m, m3ew, m3c, &c__4, rw, &c__32, 
			&iw[4], iw, &info1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 459 "SG03BU.f"
		if (info1 != 0) {
#line 460 "SG03BU.f"
		    *info = 4;
#line 461 "SG03BU.f"
		    return 0;
#line 462 "SG03BU.f"
		}
#line 463 "SG03BU.f"
		i__1 = *n - kh;
#line 463 "SG03BU.f"
		dgemm_("T", "N", &i__1, &kb, &kb, &c_b24, &b[kl + (kh + 1) * 
			b_dim1], ldb, m3c, &c__4, &c_b21, &dwork[ypt], &ldws, 
			(ftnlen)1, (ftnlen)1);
#line 465 "SG03BU.f"
		i__1 = *n - kh;
#line 465 "SG03BU.f"
		dgemm_("T", "T", &i__1, &kb, &kb, &c_b24, &a[kl + (kh + 1) * 
			a_dim1], lda, ui, &c__2, &c_b21, &dwork[wpt], &ldws, (
			ftnlen)1, (ftnlen)1);
#line 467 "SG03BU.f"
		i__1 = *n - kh;
#line 467 "SG03BU.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
#line 468 "SG03BU.f"
		    i__3 = i__ + 1, i__4 = *n - kh;
#line 468 "SG03BU.f"
		    i__2 = min(i__3,i__4);
#line 468 "SG03BU.f"
		    dgemv_("T", &i__2, &kb, &c_b24, &dwork[uiipt], &ldws, &a[
			    kh + 1 + (kh + i__) * a_dim1], &c__1, &c_b24, &
			    dwork[wpt + i__ - 1], &ldws, (ftnlen)1);
#line 471 "SG03BU.f"
/* L100: */
#line 471 "SG03BU.f"
		}
#line 472 "SG03BU.f"
		i__1 = *n - kh;
#line 472 "SG03BU.f"
		dgemm_("N", "N", &i__1, &kb, &kb, &c_b24, &dwork[wpt], &ldws, 
			&m3c[kb], &c__4, &c_b24, &dwork[ypt], &ldws, (ftnlen)
			1, (ftnlen)1);

/*              Overwrite B(KH+1:N,KH+1:N) with the triangular matrix */
/*              from the QR-factorization of the (N-KH+KB)-by-(N-KH) */
/*              matrix */

/*                          (  B(KH+1:N,KH+1:N)  ) */
/*                          (                    ) */
/*                          (       Y**T         ) . */

#line 483 "SG03BU.f"
		i__1 = kb;
#line 483 "SG03BU.f"
		for (j = 1; j <= i__1; ++j) {
#line 484 "SG03BU.f"
		    i__2 = *n - kh;
#line 484 "SG03BU.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 485 "SG03BU.f"
			x = b[kh + i__ + (kh + i__) * b_dim1];
#line 486 "SG03BU.f"
			z__ = dwork[ypt + i__ - 1 + (j - 1) * ldws];
#line 487 "SG03BU.f"
			drotg_(&x, &z__, &c__, &s);
#line 488 "SG03BU.f"
			i__3 = *n - kh - i__ + 1;
#line 488 "SG03BU.f"
			drot_(&i__3, &b[kh + i__ + (kh + i__) * b_dim1], ldb, 
				&dwork[ypt + i__ - 1 + (j - 1) * ldws], &c__1,
				 &c__, &s);
#line 490 "SG03BU.f"
/* L120: */
#line 490 "SG03BU.f"
		    }
#line 491 "SG03BU.f"
/* L140: */
#line 491 "SG03BU.f"
		}

/*              Make main diagonal elements of B(KH+1:N,KH+1:N) positive. */

#line 495 "SG03BU.f"
		i__1 = *n;
#line 495 "SG03BU.f"
		for (i__ = kh + 1; i__ <= i__1; ++i__) {
#line 496 "SG03BU.f"
		    if (b[i__ + i__ * b_dim1] < 0.) {
#line 496 "SG03BU.f"
			i__2 = *n - i__ + 1;
#line 496 "SG03BU.f"
			dscal_(&i__2, &c_b19, &b[i__ + i__ * b_dim1], ldb);
#line 496 "SG03BU.f"
		    }
#line 498 "SG03BU.f"
/* L160: */
#line 498 "SG03BU.f"
		}

/*              Overwrite right hand side with the part of the solution */
/*              computed in step II. */

#line 503 "SG03BU.f"
		i__1 = kh;
#line 503 "SG03BU.f"
		for (j = kl; j <= i__1; ++j) {
#line 504 "SG03BU.f"
		    i__2 = *n - kh;
#line 504 "SG03BU.f"
		    dcopy_(&i__2, &dwork[uiipt + (j - kl) * ldws], &c__1, &b[
			    j + (kh + 1) * b_dim1], ldb);
#line 506 "SG03BU.f"
/* L180: */
#line 506 "SG03BU.f"
		}
#line 507 "SG03BU.f"
	    }

/*           Overwrite right hand side with the part of the solution */
/*           computed in step I. */

#line 512 "SG03BU.f"
	    dlacpy_("U", &kb, &kb, ui, &c__2, &b[kl + kl * b_dim1], ldb, (
		    ftnlen)1);

#line 514 "SG03BU.f"
	    goto L20;
#line 515 "SG03BU.f"
	}
/*        END WHILE 20 */

#line 518 "SG03BU.f"
    } else {

/*        Solve equation (2). */

/*        Main Loop. Compute block column U(1:KH,KL:KH). KB denotes the */
/*        number of columns in this block column. */

#line 525 "SG03BU.f"
	kl = *n + 1;
/*        WHILE ( KL .GT. 1 ) DO */
#line 527 "SG03BU.f"
L200:
#line 527 "SG03BU.f"
	if (kl > 1) {
#line 528 "SG03BU.f"
	    kh = kl - 1;
#line 529 "SG03BU.f"
	    if (kh == 1) {
#line 530 "SG03BU.f"
		kl = 1;
#line 531 "SG03BU.f"
		kb = 1;
#line 532 "SG03BU.f"
	    } else {
#line 533 "SG03BU.f"
		if (a[kh + (kh - 1) * a_dim1] == 0.) {
#line 534 "SG03BU.f"
		    kl = kh;
#line 535 "SG03BU.f"
		    kb = 1;
#line 536 "SG03BU.f"
		} else {
#line 537 "SG03BU.f"
		    kl = kh - 1;
#line 538 "SG03BU.f"
		    kb = 2;
#line 539 "SG03BU.f"
		}
#line 540 "SG03BU.f"
	    }

/*           STEP I: Compute block U(KL:KH,KL:KH) and the auxiliary */
/*                   matrices M1 and M2. (For the moment the result */
/*                   U(KL:KH,KL:KH) is stored in UI). */

#line 546 "SG03BU.f"
	    if (kb == 1) {
/* Computing 2nd power */
#line 547 "SG03BU.f"
		d__1 = e[kl + kl * e_dim1];
/* Computing 2nd power */
#line 547 "SG03BU.f"
		d__2 = a[kl + kl * a_dim1];
#line 547 "SG03BU.f"
		delta1 = d__1 * d__1 - d__2 * d__2;
#line 548 "SG03BU.f"
		if (delta1 <= 0.) {
#line 549 "SG03BU.f"
		    *info = 3;
#line 550 "SG03BU.f"
		    return 0;
#line 551 "SG03BU.f"
		}
#line 552 "SG03BU.f"
		delta1 = sqrt(delta1);
#line 553 "SG03BU.f"
		z__ = (d__1 = b[kl + kl * b_dim1], abs(d__1)) * 2. * smlnum;
#line 554 "SG03BU.f"
		if (z__ > delta1) {
#line 555 "SG03BU.f"
		    scale1 = delta1 / z__;
#line 556 "SG03BU.f"
		    *scale = scale1 * *scale;
#line 557 "SG03BU.f"
		    i__1 = *n;
#line 557 "SG03BU.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 558 "SG03BU.f"
			dscal_(&i__, &scale1, &b[i__ * b_dim1 + 1], &c__1);
#line 559 "SG03BU.f"
/* L220: */
#line 559 "SG03BU.f"
		    }
#line 560 "SG03BU.f"
		}
#line 561 "SG03BU.f"
		ui[0] = b[kl + kl * b_dim1] / delta1;
#line 562 "SG03BU.f"
		m1[0] = a[kl + kl * a_dim1] / e[kl + kl * e_dim1];
#line 563 "SG03BU.f"
		m2[0] = delta1 / e[kl + kl * e_dim1];
#line 564 "SG03BU.f"
	    } else {

/*              If a pair of complex conjugate eigenvalues occurs, apply */
/*              (complex) Hammarling algorithm for the 2-by-2 problem. */

#line 569 "SG03BU.f"
		sg03bx_("D", "T", &a[kl + kl * a_dim1], lda, &e[kl + kl * 
			e_dim1], lde, &b[kl + kl * b_dim1], ldb, ui, &c__2, &
			scale1, m1, &c__2, m2, &c__2, &info1, (ftnlen)1, (
			ftnlen)1);
#line 572 "SG03BU.f"
		if (info1 != 0) {
#line 573 "SG03BU.f"
		    *info = info1;
#line 574 "SG03BU.f"
		    return 0;
#line 575 "SG03BU.f"
		}
#line 576 "SG03BU.f"
		if (scale1 != 1.) {
#line 577 "SG03BU.f"
		    *scale = scale1 * *scale;
#line 578 "SG03BU.f"
		    i__1 = *n;
#line 578 "SG03BU.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 579 "SG03BU.f"
			dscal_(&i__, &scale1, &b[i__ * b_dim1 + 1], &c__1);
#line 580 "SG03BU.f"
/* L240: */
#line 580 "SG03BU.f"
		    }
#line 581 "SG03BU.f"
		}
#line 582 "SG03BU.f"
	    }

#line 584 "SG03BU.f"
	    if (kl > 1) {

/*              STEP II: Compute U(1:KL-1,KL:KH) by solving a generalized */
/*                       Sylvester equation. (For the moment the result */
/*                       U(1:KL-1,KL:KH) is stored in the workspace.) */

/*              Form right hand side of the Sylvester equation. */

#line 592 "SG03BU.f"
		i__1 = kl - 1;
#line 592 "SG03BU.f"
		dgemm_("N", "T", &i__1, &kb, &kb, &c_b19, &b[kl * b_dim1 + 1],
			 ldb, m2, &c__2, &c_b21, &dwork[uiipt], &ldws, (
			ftnlen)1, (ftnlen)1);
#line 594 "SG03BU.f"
		i__1 = kl - 1;
#line 594 "SG03BU.f"
		dgemm_("N", "N", &i__1, &kb, &kb, &c_b24, &e[kl * e_dim1 + 1],
			 lde, ui, &c__2, &c_b24, &dwork[uiipt], &ldws, (
			ftnlen)1, (ftnlen)1);
#line 596 "SG03BU.f"
		dgemm_("N", "T", &kb, &kb, &kb, &c_b24, ui, &c__2, m1, &c__2, 
			&c_b21, tm, &c__2, (ftnlen)1, (ftnlen)1);
#line 598 "SG03BU.f"
		i__1 = kl - 1;
#line 598 "SG03BU.f"
		dgemm_("N", "N", &i__1, &kb, &kb, &c_b19, &a[kl * a_dim1 + 1],
			 lda, tm, &c__2, &c_b24, &dwork[uiipt], &ldws, (
			ftnlen)1, (ftnlen)1);

/*              Solve generalized Sylvester equation. */

#line 603 "SG03BU.f"
		dlaset_("A", &kb, &kb, &c_b21, &c_b19, tm, &c__2, (ftnlen)1);
#line 604 "SG03BU.f"
		i__1 = kl - 1;
#line 604 "SG03BU.f"
		sg03bw_("T", &i__1, &kb, &a[a_offset], lda, m1, &c__2, &e[
			e_offset], lde, tm, &c__2, &dwork[uiipt], &ldws, &
			scale1, &info1, (ftnlen)1);
#line 606 "SG03BU.f"
		if (info1 != 0) {
#line 606 "SG03BU.f"
		    *info = 1;
#line 606 "SG03BU.f"
		}
#line 608 "SG03BU.f"
		if (scale1 != 1.) {
#line 609 "SG03BU.f"
		    *scale = scale1 * *scale;
#line 610 "SG03BU.f"
		    i__1 = *n;
#line 610 "SG03BU.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 611 "SG03BU.f"
			dscal_(&i__, &scale1, &b[i__ * b_dim1 + 1], &c__1);
#line 612 "SG03BU.f"
/* L260: */
#line 612 "SG03BU.f"
		    }
#line 613 "SG03BU.f"
		    dscal_(&c__4, &scale1, ui, &c__1);
#line 614 "SG03BU.f"
		}

/*              STEP III: Form the right hand side matrix */
/*                        B(1:KL-1,1:KL-1) of the (smaller) Lyapunov */
/*                        equation to be solved during the next pass of */
/*                        the main loop. */

/*              Compute auxiliary matrices M3 and Y. The factorization */
/*              M3 = M3C * M3C**T is found by solving the symmetric */
/*              eigenvalue problem. */

#line 625 "SG03BU.f"
		i__1 = kb << 1;
#line 625 "SG03BU.f"
		i__2 = kb << 1;
#line 625 "SG03BU.f"
		dlaset_("U", &i__1, &i__2, &c_b21, &c_b24, m3, &c__4, (ftnlen)
			1);
#line 626 "SG03BU.f"
		dsyrk_("U", "T", &kb, &kb, &c_b19, m2, &c__2, &c_b24, m3, &
			c__4, (ftnlen)1, (ftnlen)1);
#line 627 "SG03BU.f"
		dgemm_("T", "N", &kb, &kb, &kb, &c_b19, m2, &c__2, m1, &c__2, 
			&c_b21, &m3[(kb + 1 << 2) - 4], &c__4, (ftnlen)1, (
			ftnlen)1);
#line 629 "SG03BU.f"
		dsyrk_("U", "T", &kb, &kb, &c_b19, m1, &c__2, &c_b24, &m3[kb 
			+ 1 + (kb + 1 << 2) - 5], &c__4, (ftnlen)1, (ftnlen)1)
			;
#line 631 "SG03BU.f"
		i__1 = kb << 1;
#line 631 "SG03BU.f"
		d__1 = uflt * 2.;
#line 631 "SG03BU.f"
		dsyevx_("V", "V", "U", &i__1, m3, &c__4, &c_b77, &c_b78, &
			c__1, &c__4, &d__1, &m, m3ew, m3c, &c__4, rw, &c__32, 
			&iw[4], iw, &info1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 634 "SG03BU.f"
		if (info1 != 0) {
#line 635 "SG03BU.f"
		    *info = 4;
#line 636 "SG03BU.f"
		    return 0;
#line 637 "SG03BU.f"
		}
#line 638 "SG03BU.f"
		i__1 = kl - 1;
#line 638 "SG03BU.f"
		dgemm_("N", "N", &i__1, &kb, &kb, &c_b24, &b[kl * b_dim1 + 1],
			 ldb, m3c, &c__4, &c_b21, &dwork[ypt], &ldws, (ftnlen)
			1, (ftnlen)1);
#line 640 "SG03BU.f"
		i__1 = kl - 1;
#line 640 "SG03BU.f"
		dgemm_("N", "N", &i__1, &kb, &kb, &c_b24, &a[kl * a_dim1 + 1],
			 lda, ui, &c__2, &c_b21, &dwork[wpt], &ldws, (ftnlen)
			1, (ftnlen)1);
#line 642 "SG03BU.f"
		i__1 = kl - 1;
#line 642 "SG03BU.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
#line 643 "SG03BU.f"
		    i__3 = kl - i__ + 1, i__4 = kl - 1;
#line 643 "SG03BU.f"
		    i__2 = min(i__3,i__4);
/* Computing MAX */
#line 643 "SG03BU.f"
		    i__5 = uiipt, i__6 = uiipt + i__ - 2;
/* Computing MAX */
#line 643 "SG03BU.f"
		    i__7 = i__ - 1;
#line 643 "SG03BU.f"
		    dgemv_("T", &i__2, &kb, &c_b24, &dwork[max(i__5,i__6)], &
			    ldws, &a[i__ + max(i__7,1) * a_dim1], lda, &c_b24,
			     &dwork[wpt + i__ - 1], &ldws, (ftnlen)1);
#line 647 "SG03BU.f"
/* L280: */
#line 647 "SG03BU.f"
		}
#line 648 "SG03BU.f"
		i__1 = kl - 1;
#line 648 "SG03BU.f"
		dgemm_("N", "N", &i__1, &kb, &kb, &c_b24, &dwork[wpt], &ldws, 
			&m3c[kb], &c__4, &c_b24, &dwork[ypt], &ldws, (ftnlen)
			1, (ftnlen)1);

/*              Overwrite B(1:KL-1,1:KL-1) with the triangular matrix */
/*              from the RQ-factorization of the (KL-1)-by-KH matrix */

/*                          (                        ) */
/*                          (  B(1:KL-1,1:KL-1)   Y  ) */
/*                          (                        ). */

#line 658 "SG03BU.f"
		i__1 = kb;
#line 658 "SG03BU.f"
		for (j = 1; j <= i__1; ++j) {
#line 659 "SG03BU.f"
		    for (i__ = kl - 1; i__ >= 1; --i__) {
#line 660 "SG03BU.f"
			x = b[i__ + i__ * b_dim1];
#line 661 "SG03BU.f"
			z__ = dwork[ypt + i__ - 1 + (j - 1) * ldws];
#line 662 "SG03BU.f"
			drotg_(&x, &z__, &c__, &s);
#line 663 "SG03BU.f"
			drot_(&i__, &b[i__ * b_dim1 + 1], &c__1, &dwork[ypt + 
				(j - 1) * ldws], &c__1, &c__, &s);
#line 665 "SG03BU.f"
/* L300: */
#line 665 "SG03BU.f"
		    }
#line 666 "SG03BU.f"
/* L320: */
#line 666 "SG03BU.f"
		}

/*              Make main diagonal elements of B(1:KL-1,1:KL-1) positive. */

#line 670 "SG03BU.f"
		i__1 = kl - 1;
#line 670 "SG03BU.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 671 "SG03BU.f"
		    if (b[i__ + i__ * b_dim1] < 0.) {
#line 671 "SG03BU.f"
			dscal_(&i__, &c_b19, &b[i__ * b_dim1 + 1], &c__1);
#line 671 "SG03BU.f"
		    }
#line 673 "SG03BU.f"
/* L340: */
#line 673 "SG03BU.f"
		}

/*              Overwrite right hand side with the part of the solution */
/*              computed in step II. */

#line 678 "SG03BU.f"
		i__1 = kl - 1;
#line 678 "SG03BU.f"
		dlacpy_("A", &i__1, &kb, &dwork[uiipt], &ldws, &b[kl * b_dim1 
			+ 1], ldb, (ftnlen)1);

#line 681 "SG03BU.f"
	    }

/*           Overwrite right hand side with the part of the solution */
/*           computed in step I. */

#line 686 "SG03BU.f"
	    dlacpy_("U", &kb, &kb, ui, &c__2, &b[kl + kl * b_dim1], ldb, (
		    ftnlen)1);

#line 688 "SG03BU.f"
	    goto L200;
#line 689 "SG03BU.f"
	}
/*        END WHILE 200 */

#line 692 "SG03BU.f"
    }

#line 694 "SG03BU.f"
    return 0;
/* *** Last line of SG03BU *** */
} /* sg03bu_ */

