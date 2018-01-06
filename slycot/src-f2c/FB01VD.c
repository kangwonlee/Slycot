#line 1 "FB01VD.f"
/* FB01VD.f -- translated by f2c (version 20100827).
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

#line 1 "FB01VD.f"
/* Table of constant values */

static doublereal c_b5 = 1.;
static integer c__1 = 1;
static doublereal c_b15 = 2.;
static doublereal c_b26 = 0.;
static doublereal c_b40 = -1.;

/* Subroutine */ int fb01vd_(integer *n, integer *m, integer *l, doublereal *
	p, integer *ldp, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *q, integer *ldq, 
	doublereal *r__, integer *ldr, doublereal *k, integer *ldk, 
	doublereal *tol, integer *iwork, doublereal *dwork, integer *ldwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, k_dim1, 
	    k_offset, p_dim1, p_offset, q_dim1, q_offset, r_dim1, r_offset, 
	    i__1, i__2;

    /* Local variables */
    static integer j, n1, ldw;
    extern /* Subroutine */ int mb01rd_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), dscal_(integer *, doublereal *, 
	    doublereal *, integer *), dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal rcond;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *), dtrsm_(char *, char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer jwork;
    static doublereal rnorm;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dpocon_(char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, ftnlen);
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int dpotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);


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

/*     To compute one recursion of the conventional Kalman filter */
/*     equations. This is one update of the Riccati difference equation */
/*     and the Kalman filter gain. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The actual state dimension, i.e., the order of the */
/*             matrices P      and A .  N >= 0. */
/*                       i|i-1      i */

/*     M       (input) INTEGER */
/*             The actual input dimension, i.e., the order of the matrix */
/*             Q .  M >= 0. */
/*              i */

/*     L       (input) INTEGER */
/*             The actual output dimension, i.e., the order of the matrix */
/*             R .  L >= 0. */
/*              i */

/*     P       (input/output) DOUBLE PRECISION array, dimension (LDP,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain P     , the state covariance matrix at instant */
/*                      i|i-1 */
/*             (i-1). The upper triangular part only is needed. */
/*             On exit, if INFO = 0, the leading N-by-N part of this */
/*             array contains P     , the state covariance matrix at */
/*                             i+1|i */
/*             instant i. The strictly lower triangular part is not set. */
/*             Otherwise, the leading N-by-N part of this array contains */
/*             P     , its input value. */
/*              i|i-1 */

/*     LDP     INTEGER */
/*             The leading dimension of array P.  LDP >= MAX(1,N). */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain A , */
/*                                                                 i */
/*             the state transition matrix of the discrete system at */
/*             instant i. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain B , */
/*                                                                 i */
/*             the input weight matrix of the discrete system at */
/*             instant i. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading L-by-N part of this array must contain C , */
/*                                                                 i */
/*             the output weight matrix of the discrete system at */
/*             instant i. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,L). */

/*     Q       (input) DOUBLE PRECISION array, dimension (LDQ,M) */
/*             The leading M-by-M part of this array must contain Q , */
/*                                                                 i */
/*             the input (process) noise covariance matrix at instant i. */
/*             The diagonal elements of this array are modified by the */
/*             routine, but are restored on exit. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q.  LDQ >= MAX(1,M). */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR,L) */
/*             On entry, the leading L-by-L part of this array must */
/*             contain R , the output (measurement) noise covariance */
/*                      i */
/*             matrix at instant i. */
/*             On exit, if INFO = 0, or INFO = L+1, the leading L-by-L */
/*                                                                  1/2 */
/*             upper triangular part of this array contains (RINOV )   , */
/*                                                                i */
/*             the square root (left Cholesky factor) of the covariance */
/*             matrix of the innovations at instant i. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,L). */

/*     K       (output) DOUBLE PRECISION array, dimension (LDK,L) */
/*             If INFO = 0, the leading N-by-L part of this array */
/*             contains K , the Kalman filter gain matrix at instant i. */
/*                       i */
/*             If INFO > 0, the leading N-by-L part of this array */
/*             contains the matrix product P     C'. */
/*                                          i|i-1 i */

/*     LDK     INTEGER */
/*             The leading dimension of array K.  LDK >= MAX(1,N). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used to test for near singularity of */
/*             the matrix RINOV . If the user sets TOL > 0, then the */
/*                             i */
/*             given value of TOL is used as a lower bound for the */
/*             reciprocal condition number of that matrix; a matrix whose */
/*             estimated condition number is less than 1/TOL is */
/*             considered to be nonsingular. If the user sets TOL <= 0, */
/*             then an implicitly computed, default tolerance, defined by */
/*             TOLDEF = L*L*EPS, is used instead, where EPS is the */
/*             machine precision (see LAPACK Library routine DLAMCH). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (L) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, or INFO = L+1, DWORK(1) returns an */
/*             estimate of the reciprocal of the condition number (in the */
/*             1-norm) of the matrix RINOV . */
/*                                        i */

/*     LDWORK  The length of the array DWORK. */
/*             LDWORK >= MAX(1,L*N+3*L,N*N,N*M). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -k, the k-th argument had an illegal */
/*                   value; */
/*             = k:  if INFO = k, 1 <= k <= L, the leading minor of order */
/*                   k of the matrix RINOV  is not positive-definite, and */
/*                                        i */
/*                   its Cholesky factorization could not be completed; */
/*             = L+1: the matrix RINOV  is singular, i.e., the condition */
/*                                    i */
/*                   number estimate of RINOV  (in the 1-norm) exceeds */
/*                                           i */
/*                   1/TOL. */

/*     METHOD */

/*     The conventional Kalman filter gain used at the i-th recursion */
/*     step is of the form */

/*                            -1 */
/*        K  = P     C'  RINOV  , */
/*         i    i|i-1 i       i */

/*     where RINOV  = C P     C' + R , and the state covariance matrix */
/*                i    i i|i-1 i    i */

/*     P      is updated by the discrete-time difference Riccati equation */
/*      i|i-1 */

/*        P      = A  (P      - K C P     ) A'  + B Q B'. */
/*         i+1|i    i   i|i-1    i i i|i-1   i     i i i */

/*     Using these two updates, the combined time and measurement update */
/*     of the state X      is given by */
/*                   i|i-1 */

/*        X      = A X      + A K (Y  - C X     ), */
/*         i+1|i    i i|i-1    i i  i    i i|i-1 */

/*     where Y  is the new observation at step i. */
/*            i */

/*     REFERENCES */

/*     [1] Anderson, B.D.O. and Moore, J.B. */
/*         Optimal Filtering, */
/*         Prentice Hall, Englewood Cliffs, New Jersey, 1979. */

/*     [2] Verhaegen, M.H.G. and Van Dooren, P. */
/*         Numerical Aspects of Different Kalman Filter Implementations. */
/*         IEEE Trans. Auto. Contr., AC-31, pp. 907-917, 1986. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires approximately */

/*             3   2 */
/*      3/2 x N + N  x (3 x L + M/2) */

/*     operations. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */
/*     Supersedes Release 2.0 routine FB01JD by M.H.G. Verhaegen, */
/*     M. Vanbegin, and P. Van Dooren. */

/*     REVISIONS */

/*     February 20, 1998, November 20, 2003, April 20, 2004. */

/*     KEYWORDS */

/*     Kalman filtering, optimal filtering, recursive estimation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

#line 259 "FB01VD.f"
    /* Parameter adjustments */
#line 259 "FB01VD.f"
    p_dim1 = *ldp;
#line 259 "FB01VD.f"
    p_offset = 1 + p_dim1;
#line 259 "FB01VD.f"
    p -= p_offset;
#line 259 "FB01VD.f"
    a_dim1 = *lda;
#line 259 "FB01VD.f"
    a_offset = 1 + a_dim1;
#line 259 "FB01VD.f"
    a -= a_offset;
#line 259 "FB01VD.f"
    b_dim1 = *ldb;
#line 259 "FB01VD.f"
    b_offset = 1 + b_dim1;
#line 259 "FB01VD.f"
    b -= b_offset;
#line 259 "FB01VD.f"
    c_dim1 = *ldc;
#line 259 "FB01VD.f"
    c_offset = 1 + c_dim1;
#line 259 "FB01VD.f"
    c__ -= c_offset;
#line 259 "FB01VD.f"
    q_dim1 = *ldq;
#line 259 "FB01VD.f"
    q_offset = 1 + q_dim1;
#line 259 "FB01VD.f"
    q -= q_offset;
#line 259 "FB01VD.f"
    r_dim1 = *ldr;
#line 259 "FB01VD.f"
    r_offset = 1 + r_dim1;
#line 259 "FB01VD.f"
    r__ -= r_offset;
#line 259 "FB01VD.f"
    k_dim1 = *ldk;
#line 259 "FB01VD.f"
    k_offset = 1 + k_dim1;
#line 259 "FB01VD.f"
    k -= k_offset;
#line 259 "FB01VD.f"
    --iwork;
#line 259 "FB01VD.f"
    --dwork;
#line 259 "FB01VD.f"

#line 259 "FB01VD.f"
    /* Function Body */
#line 259 "FB01VD.f"
    *info = 0;
#line 260 "FB01VD.f"
    n1 = max(1,*n);
#line 261 "FB01VD.f"
    if (*n < 0) {
#line 262 "FB01VD.f"
	*info = -1;
#line 263 "FB01VD.f"
    } else if (*m < 0) {
#line 264 "FB01VD.f"
	*info = -2;
#line 265 "FB01VD.f"
    } else if (*l < 0) {
#line 266 "FB01VD.f"
	*info = -3;
#line 267 "FB01VD.f"
    } else if (*ldp < n1) {
#line 268 "FB01VD.f"
	*info = -5;
#line 269 "FB01VD.f"
    } else if (*lda < n1) {
#line 270 "FB01VD.f"
	*info = -7;
#line 271 "FB01VD.f"
    } else if (*ldb < n1) {
#line 272 "FB01VD.f"
	*info = -9;
#line 273 "FB01VD.f"
    } else if (*ldc < max(1,*l)) {
#line 274 "FB01VD.f"
	*info = -11;
#line 275 "FB01VD.f"
    } else if (*ldq < max(1,*m)) {
#line 276 "FB01VD.f"
	*info = -13;
#line 277 "FB01VD.f"
    } else if (*ldr < max(1,*l)) {
#line 278 "FB01VD.f"
	*info = -15;
#line 279 "FB01VD.f"
    } else if (*ldk < n1) {
#line 280 "FB01VD.f"
	*info = -17;
#line 281 "FB01VD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 281 "FB01VD.f"
	i__1 = 1, i__2 = *l * *n + *l * 3, i__1 = max(i__1,i__2), i__2 = *n * 
		*n, i__1 = max(i__1,i__2), i__2 = *n * *m;
#line 281 "FB01VD.f"
	if (*ldwork < max(i__1,i__2)) {
#line 282 "FB01VD.f"
	    *info = -21;
#line 283 "FB01VD.f"
	}
#line 283 "FB01VD.f"
    }

#line 285 "FB01VD.f"
    if (*info != 0) {

/*        Error return. */

#line 289 "FB01VD.f"
	i__1 = -(*info);
#line 289 "FB01VD.f"
	xerbla_("FB01VD", &i__1, (ftnlen)6);
#line 290 "FB01VD.f"
	return 0;
#line 291 "FB01VD.f"
    }

/*     Quick return if possible. */

#line 295 "FB01VD.f"
    if (max(*n,*l) == 0) {
#line 296 "FB01VD.f"
	dwork[1] = 1.;
#line 297 "FB01VD.f"
	return 0;
#line 298 "FB01VD.f"
    }

/*     Efficiently compute RINOV = CPC' + R in R and put CP in DWORK and */
/*     PC' in K. (The content of DWORK on exit from MB01RD is used.) */
/*     Workspace: need L*N. */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code.) */

#line 308 "FB01VD.f"
    mb01rd_("Upper", "No transpose", l, n, &c_b5, &c_b5, &r__[r_offset], ldr, 
	    &c__[c_offset], ldc, &p[p_offset], ldp, &dwork[1], ldwork, info, (
	    ftnlen)5, (ftnlen)12);
#line 310 "FB01VD.f"
    ldw = max(1,*l);

#line 312 "FB01VD.f"
    i__1 = *l;
#line 312 "FB01VD.f"
    for (j = 1; j <= i__1; ++j) {
#line 313 "FB01VD.f"
	dcopy_(n, &dwork[j], &ldw, &k[j * k_dim1 + 1], &c__1);
#line 314 "FB01VD.f"
/* L10: */
#line 314 "FB01VD.f"
    }

#line 316 "FB01VD.f"
    dlacpy_("Full", l, n, &c__[c_offset], ldc, &dwork[1], &ldw, (ftnlen)4);
#line 317 "FB01VD.f"
    dtrmm_("Right", "Upper", "Transpose", "Non-unit", l, n, &c_b5, &p[
	    p_offset], ldp, &dwork[1], &ldw, (ftnlen)5, (ftnlen)5, (ftnlen)9, 
	    (ftnlen)8);
#line 319 "FB01VD.f"
    i__1 = *ldp + 1;
#line 319 "FB01VD.f"
    dscal_(n, &c_b15, &p[p_offset], &i__1);

#line 321 "FB01VD.f"
    i__1 = *l;
#line 321 "FB01VD.f"
    for (j = 1; j <= i__1; ++j) {
#line 322 "FB01VD.f"
	daxpy_(n, &c_b5, &k[j * k_dim1 + 1], &c__1, &dwork[j], &ldw);
#line 323 "FB01VD.f"
	dcopy_(n, &dwork[j], &ldw, &k[j * k_dim1 + 1], &c__1);
#line 324 "FB01VD.f"
/* L20: */
#line 324 "FB01VD.f"
    }

/*     Calculate the Cholesky decomposition U'U of the innovation */
/*     covariance matrix RINOV, and its reciprocal condition number. */
/*     Workspace: need L*N + 3*L. */

#line 330 "FB01VD.f"
    jwork = *l * *n + 1;
#line 331 "FB01VD.f"
    rnorm = dlansy_("1-norm", "Upper", l, &r__[r_offset], ldr, &dwork[jwork], 
	    (ftnlen)6, (ftnlen)5);

#line 333 "FB01VD.f"
    toldef = *tol;
#line 334 "FB01VD.f"
    if (toldef <= 0.) {
#line 334 "FB01VD.f"
	toldef = (doublereal) (*l * *l) * dlamch_("Epsilon", (ftnlen)7);
#line 334 "FB01VD.f"
    }
#line 336 "FB01VD.f"
    dpotrf_("Upper", l, &r__[r_offset], ldr, info, (ftnlen)5);
#line 337 "FB01VD.f"
    if (*info != 0) {
#line 337 "FB01VD.f"
	return 0;
#line 337 "FB01VD.f"
    }

#line 340 "FB01VD.f"
    dpocon_("Upper", l, &r__[r_offset], ldr, &rnorm, &rcond, &dwork[jwork], &
	    iwork[1], info, (ftnlen)5);

#line 343 "FB01VD.f"
    if (rcond < toldef) {

/*        Error return: RINOV is numerically singular. */

#line 347 "FB01VD.f"
	*info = *l + 1;
#line 348 "FB01VD.f"
	dwork[1] = rcond;
#line 349 "FB01VD.f"
	return 0;
#line 350 "FB01VD.f"
    }

#line 352 "FB01VD.f"
    if (*l > 1) {
#line 352 "FB01VD.f"
	i__1 = *l - 1;
#line 352 "FB01VD.f"
	i__2 = *l - 1;
#line 352 "FB01VD.f"
	dlaset_("Lower", &i__1, &i__2, &c_b26, &c_b26, &r__[r_dim1 + 2], ldr, 
		(ftnlen)5);
#line 352 "FB01VD.f"
    }
/*                                                          -1 */
/*     Calculate the Kalman filter gain matrix  K = PC'RINOV . */
/*     Workspace: need L*N. */

#line 358 "FB01VD.f"
    dtrsm_("Right", "Upper", "No transpose", "Non-unit", n, l, &c_b5, &r__[
	    r_offset], ldr, &k[k_offset], ldk, (ftnlen)5, (ftnlen)5, (ftnlen)
	    12, (ftnlen)8);
#line 360 "FB01VD.f"
    dtrsm_("Right", "Upper", "Transpose", "Non-unit", n, l, &c_b5, &r__[
	    r_offset], ldr, &k[k_offset], ldk, (ftnlen)5, (ftnlen)5, (ftnlen)
	    9, (ftnlen)8);

/*     First part of the Riccati equation update: compute A(P-KCP)A'. */
/*     The upper triangular part of the symmetric matrix P-KCP is formed. */
/*     Workspace: need max(L*N,N*N). */

#line 367 "FB01VD.f"
    jwork = 1;

#line 369 "FB01VD.f"
    i__1 = *n;
#line 369 "FB01VD.f"
    for (j = 1; j <= i__1; ++j) {
#line 370 "FB01VD.f"
	dgemv_("No transpose", &j, l, &c_b40, &k[k_offset], ldk, &dwork[jwork]
		, &c__1, &c_b5, &p[j * p_dim1 + 1], &c__1, (ftnlen)12);
#line 372 "FB01VD.f"
	jwork += *l;
#line 373 "FB01VD.f"
/* L30: */
#line 373 "FB01VD.f"
    }

#line 375 "FB01VD.f"
    mb01rd_("Upper", "No transpose", n, n, &c_b26, &c_b5, &p[p_offset], ldp, &
	    a[a_offset], lda, &p[p_offset], ldp, &dwork[1], ldwork, info, (
	    ftnlen)5, (ftnlen)12);

/*     Second part of the Riccati equation update: add BQB'. */
/*     Workspace: need N*M. */

#line 381 "FB01VD.f"
    mb01rd_("Upper", "No transpose", n, m, &c_b5, &c_b5, &p[p_offset], ldp, &
	    b[b_offset], ldb, &q[q_offset], ldq, &dwork[1], ldwork, info, (
	    ftnlen)5, (ftnlen)12);
#line 383 "FB01VD.f"
    i__1 = *ldq + 1;
#line 383 "FB01VD.f"
    dscal_(m, &c_b15, &q[q_offset], &i__1);

/*     Set the reciprocal of the condition number estimate. */

#line 387 "FB01VD.f"
    dwork[1] = rcond;

#line 389 "FB01VD.f"
    return 0;
/* *** Last line of FB01VD *** */
} /* fb01vd_ */

