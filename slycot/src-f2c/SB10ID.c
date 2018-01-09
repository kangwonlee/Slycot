#line 1 "SB10ID.f"
/* SB10ID.f -- translated by f2c (version 20100827).
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

#line 1 "SB10ID.f"
/* Table of constant values */

static doublereal c_b5 = 1.;
static doublereal c_b6 = 0.;
static doublereal c_b39 = -1.;

/* Subroutine */ int sb10id_(integer *n, integer *m, integer *np, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, doublereal *d__, integer *ldd, doublereal *factor, 
	integer *nk, doublereal *ak, integer *ldak, doublereal *bk, integer *
	ldbk, doublereal *ck, integer *ldck, doublereal *dk, integer *lddk, 
	doublereal *rcond, integer *iwork, doublereal *dwork, integer *ldwork,
	 logical *bwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ak_dim1, ak_offset, b_dim1, b_offset, bk_dim1, 
	    bk_offset, c_dim1, c_offset, ck_dim1, ck_offset, d_dim1, d_offset,
	     dk_dim1, dk_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, i1, i2, i3, i4, i5, i6, i7, i8, i9, n2, i10, i11, 
	    i12, i13, ns, lwa;
    static doublereal sep;
    static integer sdim;
    static doublereal ferr;
    static char hinv[1];
    static integer iwrk, info2;
    static doublereal gamma;
    extern /* Subroutine */ int sb10jd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *), dgees_(char *, 
	    char *, L_fp, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, logical *, integer *, ftnlen, ftnlen), dgemm_(char *, 
	    char *, integer *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen), mb02vd_(char *, integer *, integer *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen), sb02rd_(char *, char *, char *, char *, char *
	    , char *, char *, char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, logical *, integer *, ftnlen, ftnlen, ftnlen, ftnlen, 
	    ftnlen, ftnlen, ftnlen, ftnlen, ftnlen), dtrsm_(char *, char *, 
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), dsyrk_(char *, char *, integer *, integer *, doublereal *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen), dlacpy_(char *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, ftnlen), dlaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen);
    extern logical select_();
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lwamax;
    extern /* Subroutine */ int dpotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer minwrk;
    extern /* Subroutine */ int dpotrs_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
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

/*     To compute the matrices of the positive feedback controller */

/*              | Ak | Bk | */
/*          K = |----|----| */
/*              | Ck | Dk | */

/*     for the shaped plant */

/*              | A | B | */
/*          G = |---|---| */
/*              | C | D | */

/*     in the McFarlane/Glover Loop Shaping Design Procedure. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the plant.  N >= 0. */

/*     M       (input) INTEGER */
/*             The column size of the matrix B.  M >= 0. */

/*     NP      (input) INTEGER */
/*             The row size of the matrix C.  NP >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             system state matrix A of the shaped plant. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             system input matrix B of the shaped plant. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading NP-by-N part of this array must contain the */
/*             system output matrix C of the shaped plant. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= max(1,NP). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading NP-by-M part of this array must contain the */
/*             system matrix D of the shaped plant. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= max(1,NP). */

/*     FACTOR  (input) DOUBLE PRECISION */
/*             = 1 implies that an optimal controller is required; */
/*             > 1 implies that a suboptimal controller is required, */
/*                 achieving a performance FACTOR less than optimal. */
/*             FACTOR >= 1. */

/*     NK      (output) INTEGER */
/*             The order of the positive feedback controller.  NK <= N. */

/*     AK      (output) DOUBLE PRECISION array, dimension (LDAK,N) */
/*             The leading NK-by-NK part of this array contains the */
/*             controller state matrix Ak. */

/*     LDAK    INTEGER */
/*             The leading dimension of the array AK.  LDAK >= max(1,N). */

/*     BK      (output) DOUBLE PRECISION array, dimension (LDBK,NP) */
/*             The leading NK-by-NP part of this array contains the */
/*             controller input matrix Bk. */

/*     LDBK    INTEGER */
/*             The leading dimension of the array BK.  LDBK >= max(1,N). */

/*     CK      (output) DOUBLE PRECISION array, dimension (LDCK,N) */
/*             The leading M-by-NK part of this array contains the */
/*             controller output matrix Ck. */

/*     LDCK    INTEGER */
/*             The leading dimension of the array CK.  LDCK >= max(1,M). */

/*     DK      (output) DOUBLE PRECISION array, dimension (LDDK,NP) */
/*             The leading M-by-NP part of this array contains the */
/*             controller matrix Dk. */

/*     LDDK    INTEGER */
/*             The leading dimension of the array DK.  LDDK >= max(1,M). */

/*     RCOND   (output) DOUBLE PRECISION array, dimension (2) */
/*             RCOND(1) contains an estimate of the reciprocal condition */
/*                      number of the X-Riccati equation; */
/*             RCOND(2) contains an estimate of the reciprocal condition */
/*                      number of the Z-Riccati equation. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension max(2*N,N*N,M,NP) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= 4*N*N + M*M + NP*NP + 2*M*N + N*NP + 4*N + */
/*                       max( 6*N*N + 5 + max(1,4*N*N+8*N), N*NP + 2*N ). */
/*             For good performance, LDWORK must generally be larger. */
/*             An upper bound of LDWORK in the above formula is */
/*             LDWORK >= 10*N*N + M*M + NP*NP + 2*M*N + 2*N*NP + 4*N + */
/*                       5 + max(1,4*N*N+8*N). */

/*     BWORK   LOGICAL array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the X-Riccati equation is not solved successfully; */
/*             = 2:  the Z-Riccati equation is not solved successfully; */
/*             = 3:  the iteration to compute eigenvalues or singular */
/*                   values failed to converge; */
/*             = 4:  the matrix Ip - D*Dk is singular; */
/*             = 5:  the matrix Im - Dk*D is singular; */
/*             = 6:  the closed-loop system is unstable. */

/*     METHOD */

/*     The routine implements the formulas given in [1]. */

/*     REFERENCES */

/*     [1] McFarlane, D. and Glover, K. */
/*         A loop shaping design procedure using H_infinity synthesis. */
/*         IEEE Trans. Automat. Control, vol. AC-37, no. 6, pp. 759-769, */
/*         1992. */

/*     NUMERICAL ASPECTS */

/*     The accuracy of the results depends on the conditioning of the */
/*     two Riccati equations solved in the controller design (see the */
/*     output parameter RCOND). */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 2000. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2000, */
/*     Feb. 2001. */

/*     KEYWORDS */

/*     H_infinity control, Loop-shaping design, Robust control. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test input parameters. */

#line 228 "SB10ID.f"
    /* Parameter adjustments */
#line 228 "SB10ID.f"
    a_dim1 = *lda;
#line 228 "SB10ID.f"
    a_offset = 1 + a_dim1;
#line 228 "SB10ID.f"
    a -= a_offset;
#line 228 "SB10ID.f"
    b_dim1 = *ldb;
#line 228 "SB10ID.f"
    b_offset = 1 + b_dim1;
#line 228 "SB10ID.f"
    b -= b_offset;
#line 228 "SB10ID.f"
    c_dim1 = *ldc;
#line 228 "SB10ID.f"
    c_offset = 1 + c_dim1;
#line 228 "SB10ID.f"
    c__ -= c_offset;
#line 228 "SB10ID.f"
    d_dim1 = *ldd;
#line 228 "SB10ID.f"
    d_offset = 1 + d_dim1;
#line 228 "SB10ID.f"
    d__ -= d_offset;
#line 228 "SB10ID.f"
    ak_dim1 = *ldak;
#line 228 "SB10ID.f"
    ak_offset = 1 + ak_dim1;
#line 228 "SB10ID.f"
    ak -= ak_offset;
#line 228 "SB10ID.f"
    bk_dim1 = *ldbk;
#line 228 "SB10ID.f"
    bk_offset = 1 + bk_dim1;
#line 228 "SB10ID.f"
    bk -= bk_offset;
#line 228 "SB10ID.f"
    ck_dim1 = *ldck;
#line 228 "SB10ID.f"
    ck_offset = 1 + ck_dim1;
#line 228 "SB10ID.f"
    ck -= ck_offset;
#line 228 "SB10ID.f"
    dk_dim1 = *lddk;
#line 228 "SB10ID.f"
    dk_offset = 1 + dk_dim1;
#line 228 "SB10ID.f"
    dk -= dk_offset;
#line 228 "SB10ID.f"
    --rcond;
#line 228 "SB10ID.f"
    --iwork;
#line 228 "SB10ID.f"
    --dwork;
#line 228 "SB10ID.f"
    --bwork;
#line 228 "SB10ID.f"

#line 228 "SB10ID.f"
    /* Function Body */
#line 228 "SB10ID.f"
    *info = 0;
#line 229 "SB10ID.f"
    if (*n < 0) {
#line 230 "SB10ID.f"
	*info = -1;
#line 231 "SB10ID.f"
    } else if (*m < 0) {
#line 232 "SB10ID.f"
	*info = -2;
#line 233 "SB10ID.f"
    } else if (*np < 0) {
#line 234 "SB10ID.f"
	*info = -3;
#line 235 "SB10ID.f"
    } else if (*lda < max(1,*n)) {
#line 236 "SB10ID.f"
	*info = -5;
#line 237 "SB10ID.f"
    } else if (*ldb < max(1,*n)) {
#line 238 "SB10ID.f"
	*info = -7;
#line 239 "SB10ID.f"
    } else if (*ldc < max(1,*np)) {
#line 240 "SB10ID.f"
	*info = -9;
#line 241 "SB10ID.f"
    } else if (*ldd < max(1,*np)) {
#line 242 "SB10ID.f"
	*info = -11;
#line 243 "SB10ID.f"
    } else if (*factor < 1.) {
#line 244 "SB10ID.f"
	*info = -12;
#line 245 "SB10ID.f"
    } else if (*ldak < max(1,*n)) {
#line 246 "SB10ID.f"
	*info = -15;
#line 247 "SB10ID.f"
    } else if (*ldbk < max(1,*n)) {
#line 248 "SB10ID.f"
	*info = -17;
#line 249 "SB10ID.f"
    } else if (*ldck < max(1,*m)) {
#line 250 "SB10ID.f"
	*info = -19;
#line 251 "SB10ID.f"
    } else if (*lddk < max(1,*m)) {
#line 252 "SB10ID.f"
	*info = -21;
#line 253 "SB10ID.f"
    }

/*     Compute workspace. */

/* Computing MAX */
/* Computing MAX */
#line 257 "SB10ID.f"
    i__3 = 1, i__4 = (*n << 2) * *n + (*n << 3);
#line 257 "SB10ID.f"
    i__1 = *n * 6 * *n + 5 + max(i__3,i__4), i__2 = *n * *np + (*n << 1);
#line 257 "SB10ID.f"
    minwrk = (*n << 2) * *n + *m * *m + *np * *np + (*m << 1) * *n + *n * *np 
	    + (*n << 2) + max(i__1,i__2);
#line 259 "SB10ID.f"
    if (*ldwork < minwrk) {
#line 260 "SB10ID.f"
	*info = -25;
#line 261 "SB10ID.f"
    }
#line 262 "SB10ID.f"
    if (*info != 0) {
#line 263 "SB10ID.f"
	i__1 = -(*info);
#line 263 "SB10ID.f"
	xerbla_("SB10ID", &i__1, (ftnlen)6);
#line 264 "SB10ID.f"
	return 0;
#line 265 "SB10ID.f"
    }

/*     Quick return if possible. */

#line 269 "SB10ID.f"
    if (*n == 0 || *m == 0 || *np == 0) {
#line 270 "SB10ID.f"
	rcond[1] = 1.;
#line 271 "SB10ID.f"
	rcond[2] = 1.;
#line 272 "SB10ID.f"
	dwork[1] = 1.;
#line 273 "SB10ID.f"
	return 0;
#line 274 "SB10ID.f"
    }

/*     Workspace usage. */

#line 278 "SB10ID.f"
    i1 = *n * *n;
#line 279 "SB10ID.f"
    i2 = i1 + *n * *n;
#line 280 "SB10ID.f"
    i3 = i2 + *m * *n;
#line 281 "SB10ID.f"
    i4 = i3 + *m * *n;
#line 282 "SB10ID.f"
    i5 = i4 + *m * *m;
#line 283 "SB10ID.f"
    i6 = i5 + *np * *np;
#line 284 "SB10ID.f"
    i7 = i6 + *np * *n;
#line 285 "SB10ID.f"
    i8 = i7 + *n * *n;
#line 286 "SB10ID.f"
    i9 = i8 + *n * *n;
#line 287 "SB10ID.f"
    i10 = i9 + *n * *n;
#line 288 "SB10ID.f"
    i11 = i10 + *n * *n;
#line 289 "SB10ID.f"
    i12 = i11 + (*n << 1);
#line 290 "SB10ID.f"
    i13 = i12 + (*n << 1);

#line 292 "SB10ID.f"
    iwrk = i13 + (*n << 2) * *n;

/*     Compute D'*C . */

#line 296 "SB10ID.f"
    dgemm_("T", "N", m, n, np, &c_b5, &d__[d_offset], ldd, &c__[c_offset], 
	    ldc, &c_b6, &dwork[i2 + 1], m, (ftnlen)1, (ftnlen)1);

/*     Compute S = Im + D'*D . */

#line 301 "SB10ID.f"
    dlaset_("U", m, m, &c_b6, &c_b5, &dwork[i4 + 1], m, (ftnlen)1);
#line 302 "SB10ID.f"
    dsyrk_("U", "T", m, np, &c_b5, &d__[d_offset], ldd, &c_b5, &dwork[i4 + 1],
	     m, (ftnlen)1, (ftnlen)1);

/*     Factorize S, S = T'*T, with T upper triangular. */

#line 306 "SB10ID.f"
    dpotrf_("U", m, &dwork[i4 + 1], m, &info2, (ftnlen)1);

/*              -1 */
/*     Compute S  D'*C . */

#line 311 "SB10ID.f"
    dpotrs_("U", m, n, &dwork[i4 + 1], m, &dwork[i2 + 1], m, &info2, (ftnlen)
	    1);

/*                -1 */
/*     Compute B*T  . */

#line 317 "SB10ID.f"
    dlacpy_("F", n, m, &b[b_offset], ldb, &dwork[i3 + 1], n, (ftnlen)1);
#line 318 "SB10ID.f"
    dtrsm_("R", "U", "N", "N", n, m, &c_b5, &dwork[i4 + 1], m, &dwork[i3 + 1],
	     n, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Compute R = Ip + D*D' . */

#line 323 "SB10ID.f"
    dlaset_("U", np, np, &c_b6, &c_b5, &dwork[i5 + 1], np, (ftnlen)1);
#line 324 "SB10ID.f"
    dsyrk_("U", "N", np, m, &c_b5, &d__[d_offset], ldd, &c_b5, &dwork[i5 + 1],
	     np, (ftnlen)1, (ftnlen)1);

/*     Factorize R, R = U'*U, with U upper triangular. */

#line 328 "SB10ID.f"
    dpotrf_("U", np, &dwork[i5 + 1], np, &info2, (ftnlen)1);

/*              -T */
/*     Compute U  C . */

#line 333 "SB10ID.f"
    dlacpy_("F", np, n, &c__[c_offset], ldc, &dwork[i6 + 1], np, (ftnlen)1);
#line 334 "SB10ID.f"
    dtrsm_("L", "U", "T", "N", np, n, &c_b5, &dwork[i5 + 1], np, &dwork[i6 + 
	    1], np, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                         -1 */
/*     Compute Ar = A - B*S  D'*C . */

#line 340 "SB10ID.f"
    dlacpy_("F", n, n, &a[a_offset], lda, &dwork[i7 + 1], n, (ftnlen)1);
#line 341 "SB10ID.f"
    dgemm_("N", "N", n, n, m, &c_b39, &b[b_offset], ldb, &dwork[i2 + 1], m, &
	    c_b5, &dwork[i7 + 1], n, (ftnlen)1, (ftnlen)1);

/*                                            -1 */
/*     Compute the upper triangle of Cr = C'*R  *C . */

#line 347 "SB10ID.f"
    dsyrk_("U", "T", n, np, &c_b5, &dwork[i6 + 1], np, &c_b6, &dwork[i8 + 1], 
	    n, (ftnlen)1, (ftnlen)1);

/*                                           -1 */
/*     Compute the upper triangle of Dr = B*S  B' . */

#line 353 "SB10ID.f"
    dsyrk_("U", "N", n, m, &c_b5, &dwork[i3 + 1], n, &c_b6, &dwork[i9 + 1], n,
	     (ftnlen)1, (ftnlen)1);

/*     Solution of the Riccati equation Ar'*X + X*Ar + Cr - X*Dr*X = 0 . */
/*     Workspace:    need   10*N*N + M*M + NP*NP + 2*M*N + N*NP + 4*N + */
/*                                   5 + max(1,4*N*N+8*N). */
/*                   prefer larger. */
/*                   AK is used as workspace. */

#line 362 "SB10ID.f"
    n2 = *n << 1;
#line 363 "SB10ID.f"
    i__1 = *ldwork - iwrk;
#line 363 "SB10ID.f"
    sb02rd_("A", "C", hinv, "N", "U", "G", "S", "N", "O", n, &dwork[i7 + 1], 
	    n, &dwork[i10 + 1], n, &ak[ak_offset], ldak, &dwork[i9 + 1], n, &
	    dwork[i8 + 1], n, &dwork[1], n, &sep, &rcond[1], &ferr, &dwork[
	    i11 + 1], &dwork[i12 + 1], &dwork[i13 + 1], &n2, &iwork[1], &
	    dwork[iwrk + 1], &i__1, &bwork[1], &info2, (ftnlen)1, (ftnlen)1, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
	    ftnlen)1);
#line 369 "SB10ID.f"
    if (info2 != 0) {
#line 370 "SB10ID.f"
	*info = 1;
#line 371 "SB10ID.f"
	return 0;
#line 372 "SB10ID.f"
    }
#line 373 "SB10ID.f"
    lwa = (integer) dwork[iwrk + 1] + iwrk;
#line 374 "SB10ID.f"
    lwamax = max(minwrk,lwa);

/*     Solution of the Riccati equation Ar*Z + Z*Ar' + Dr - Z*Cr*Z = 0 . */

#line 378 "SB10ID.f"
    i__1 = *ldwork - iwrk;
#line 378 "SB10ID.f"
    sb02rd_("A", "C", hinv, "T", "U", "G", "S", "N", "O", n, &dwork[i7 + 1], 
	    n, &dwork[i10 + 1], n, &ak[ak_offset], ldak, &dwork[i8 + 1], n, &
	    dwork[i9 + 1], n, &dwork[i1 + 1], n, &sep, &rcond[2], &ferr, &
	    dwork[i11 + 1], &dwork[i12 + 1], &dwork[i13 + 1], &n2, &iwork[1], 
	    &dwork[iwrk + 1], &i__1, &bwork[1], &info2, (ftnlen)1, (ftnlen)1, 
	    (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
	    (ftnlen)1);
#line 384 "SB10ID.f"
    if (info2 != 0) {
#line 385 "SB10ID.f"
	*info = 2;
#line 386 "SB10ID.f"
	return 0;
#line 387 "SB10ID.f"
    }
#line 388 "SB10ID.f"
    lwa = (integer) dwork[iwrk + 1] + iwrk;
#line 389 "SB10ID.f"
    lwamax = max(lwa,lwamax);

/*                      -1        -1 */
/*     Compute F1 = -( S  D'*C + S  B'*X ) . */

#line 394 "SB10ID.f"
    dtrsm_("R", "U", "T", "N", n, m, &c_b5, &dwork[i4 + 1], m, &dwork[i3 + 1],
	     n, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 396 "SB10ID.f"
    dgemm_("T", "N", m, n, n, &c_b39, &dwork[i3 + 1], n, &dwork[1], n, &c_b39,
	     &dwork[i2 + 1], m, (ftnlen)1, (ftnlen)1);

/*     Compute gamma . */

#line 401 "SB10ID.f"
    dgemm_("N", "N", n, n, n, &c_b5, &dwork[1], n, &dwork[i1 + 1], n, &c_b6, &
	    dwork[i7 + 1], n, (ftnlen)1, (ftnlen)1);
#line 403 "SB10ID.f"
    i__1 = *ldwork - iwrk;
#line 403 "SB10ID.f"
    dgees_("N", "N", (L_fp)select_, n, &dwork[i7 + 1], n, &sdim, &dwork[i11 + 
	    1], &dwork[i12 + 1], &dwork[iwrk + 1], n, &dwork[iwrk + 1], &i__1,
	     &bwork[1], &info2, (ftnlen)1, (ftnlen)1);
#line 406 "SB10ID.f"
    if (info2 != 0) {
#line 407 "SB10ID.f"
	*info = 3;
#line 408 "SB10ID.f"
	return 0;
#line 409 "SB10ID.f"
    }
#line 410 "SB10ID.f"
    lwa = (integer) dwork[iwrk + 1] + iwrk;
#line 411 "SB10ID.f"
    lwamax = max(lwa,lwamax);
#line 412 "SB10ID.f"
    gamma = 0.;
#line 413 "SB10ID.f"
    i__1 = *n;
#line 413 "SB10ID.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 414 "SB10ID.f"
	d__1 = gamma, d__2 = dwork[i11 + i__];
#line 414 "SB10ID.f"
	gamma = max(d__1,d__2);
#line 415 "SB10ID.f"
/* L10: */
#line 415 "SB10ID.f"
    }
#line 416 "SB10ID.f"
    gamma = *factor * sqrt(gamma + 1.);

/*     Workspace usage. */
/*     Workspace:    need   4*N*N + M*N + N*NP. */

#line 421 "SB10ID.f"
    i4 = i3 + *n * *n;
#line 422 "SB10ID.f"
    i5 = i4 + *n * *n;

/*     Compute Ac = A + B*F1 . */

#line 426 "SB10ID.f"
    dlacpy_("F", n, n, &a[a_offset], lda, &dwork[i4 + 1], n, (ftnlen)1);
#line 427 "SB10ID.f"
    dgemm_("N", "N", n, n, m, &c_b5, &b[b_offset], ldb, &dwork[i2 + 1], m, &
	    c_b5, &dwork[i4 + 1], n, (ftnlen)1, (ftnlen)1);

/*     Compute W1' = (1-gamma^2)*In + Z*X . */

#line 432 "SB10ID.f"
    d__1 = 1. - gamma * gamma;
#line 432 "SB10ID.f"
    dlaset_("F", n, n, &c_b6, &d__1, &dwork[i3 + 1], n, (ftnlen)1);
#line 433 "SB10ID.f"
    dgemm_("N", "N", n, n, n, &c_b5, &dwork[i1 + 1], n, &dwork[1], n, &c_b5, &
	    dwork[i3 + 1], n, (ftnlen)1, (ftnlen)1);

/*     Compute Bcp = gamma^2*Z*C' . */

#line 438 "SB10ID.f"
    d__1 = gamma * gamma;
#line 438 "SB10ID.f"
    dgemm_("N", "T", n, np, n, &d__1, &dwork[i1 + 1], n, &c__[c_offset], ldc, 
	    &c_b6, &bk[bk_offset], ldbk, (ftnlen)1, (ftnlen)1);

/*     Compute C + D*F1 . */

#line 443 "SB10ID.f"
    dlacpy_("F", np, n, &c__[c_offset], ldc, &dwork[i5 + 1], np, (ftnlen)1);
#line 444 "SB10ID.f"
    dgemm_("N", "N", np, n, m, &c_b5, &d__[d_offset], ldd, &dwork[i2 + 1], m, 
	    &c_b5, &dwork[i5 + 1], np, (ftnlen)1, (ftnlen)1);

/*     Compute Acp = W1'*Ac + gamma^2*Z*C'*(C+D*F1) . */

#line 449 "SB10ID.f"
    dgemm_("N", "N", n, n, n, &c_b5, &dwork[i3 + 1], n, &dwork[i4 + 1], n, &
	    c_b6, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);
#line 451 "SB10ID.f"
    dgemm_("N", "N", n, n, np, &c_b5, &bk[bk_offset], ldbk, &dwork[i5 + 1], 
	    np, &c_b5, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);

/*     Compute Ccp = B'*X . */

#line 456 "SB10ID.f"
    dgemm_("T", "N", m, n, n, &c_b5, &b[b_offset], ldb, &dwork[1], n, &c_b6, &
	    ck[ck_offset], ldck, (ftnlen)1, (ftnlen)1);

/*     Set Dcp = -D' . */

#line 461 "SB10ID.f"
    i__1 = *m;
#line 461 "SB10ID.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 462 "SB10ID.f"
	i__2 = *np;
#line 462 "SB10ID.f"
	for (j = 1; j <= i__2; ++j) {
#line 463 "SB10ID.f"
	    dk[i__ + j * dk_dim1] = -d__[j + i__ * d_dim1];
#line 464 "SB10ID.f"
/* L20: */
#line 464 "SB10ID.f"
	}
#line 465 "SB10ID.f"
/* L30: */
#line 465 "SB10ID.f"
    }

#line 467 "SB10ID.f"
    iwrk = i4;

/*     Reduce the generalized state-space description to a regular one. */
/*     Workspace:             need   3*N*N + M*N. */
/*     Additional workspace:  need   2*N*N + 2*N + N*MAX(5,N+M+NP). */
/*                            prefer larger. */

#line 474 "SB10ID.f"
    i__1 = *ldwork - iwrk;
#line 474 "SB10ID.f"
    sb10jd_(n, np, m, &ak[ak_offset], ldak, &bk[bk_offset], ldbk, &ck[
	    ck_offset], ldck, &dk[dk_offset], lddk, &dwork[i3 + 1], n, nk, &
	    dwork[iwrk + 1], &i__1, &info2);
#line 477 "SB10ID.f"
    if (info2 != 0) {
#line 478 "SB10ID.f"
	*info = 3;
#line 479 "SB10ID.f"
	return 0;
#line 480 "SB10ID.f"
    }
#line 481 "SB10ID.f"
    lwa = (integer) dwork[iwrk + 1] + iwrk;
#line 482 "SB10ID.f"
    lwamax = max(lwa,lwamax);

/*     Workspace usage. */
/*     Workspace:    need   4*N*N + M*M + NP*NP + 2*M*N + 2*N*NP. */
/*                          (NK <= N.) */

#line 488 "SB10ID.f"
    i2 = *np * *np;
#line 489 "SB10ID.f"
    i3 = i2 + *nk * *np;
#line 490 "SB10ID.f"
    i4 = i3 + *m * *m;
#line 491 "SB10ID.f"
    i5 = i4 + *n * *m;
#line 492 "SB10ID.f"
    i6 = i5 + *np * *nk;
#line 493 "SB10ID.f"
    i7 = i6 + *m * *n;

#line 495 "SB10ID.f"
    iwrk = i7 + (*n + *nk) * (*n + *nk);

/*     Compute Ip - D*Dk . */

#line 499 "SB10ID.f"
    dlaset_("Full", np, np, &c_b6, &c_b5, &dwork[1], np, (ftnlen)4);
#line 500 "SB10ID.f"
    dgemm_("N", "N", np, np, m, &c_b39, &d__[d_offset], ldd, &dk[dk_offset], 
	    lddk, &c_b5, &dwork[1], np, (ftnlen)1, (ftnlen)1);

/*                         -1 */
/*     Compute Bk*(Ip-D*Dk)  . */

#line 506 "SB10ID.f"
    dlacpy_("F", nk, np, &bk[bk_offset], ldbk, &dwork[i2 + 1], nk, (ftnlen)1);
#line 507 "SB10ID.f"
    mb02vd_("N", nk, np, &dwork[1], np, &iwork[1], &dwork[i2 + 1], nk, &info2,
	     (ftnlen)1);
#line 509 "SB10ID.f"
    if (info2 != 0) {
#line 510 "SB10ID.f"
	*info = 4;
#line 511 "SB10ID.f"
	return 0;
#line 512 "SB10ID.f"
    }

/*     Compute Im - Dk*D . */

#line 516 "SB10ID.f"
    dlaset_("Full", m, m, &c_b6, &c_b5, &dwork[i3 + 1], m, (ftnlen)4);
#line 517 "SB10ID.f"
    dgemm_("N", "N", m, m, np, &c_b39, &dk[dk_offset], lddk, &d__[d_offset], 
	    ldd, &c_b5, &dwork[i3 + 1], m, (ftnlen)1, (ftnlen)1);

/*                        -1 */
/*     Compute B*(Im-Dk*D)  . */

#line 523 "SB10ID.f"
    dlacpy_("F", n, m, &b[b_offset], ldb, &dwork[i4 + 1], n, (ftnlen)1);
#line 524 "SB10ID.f"
    mb02vd_("N", n, m, &dwork[i3 + 1], m, &iwork[1], &dwork[i4 + 1], n, &
	    info2, (ftnlen)1);
#line 526 "SB10ID.f"
    if (info2 != 0) {
#line 527 "SB10ID.f"
	*info = 5;
#line 528 "SB10ID.f"
	return 0;
#line 529 "SB10ID.f"
    }

/*     Compute D*Ck . */

#line 533 "SB10ID.f"
    dgemm_("N", "N", np, nk, m, &c_b5, &d__[d_offset], ldd, &ck[ck_offset], 
	    ldck, &c_b6, &dwork[i5 + 1], np, (ftnlen)1, (ftnlen)1);

/*     Compute Dk*C . */

#line 538 "SB10ID.f"
    dgemm_("N", "N", m, n, np, &c_b5, &dk[dk_offset], lddk, &c__[c_offset], 
	    ldc, &c_b6, &dwork[i6 + 1], m, (ftnlen)1, (ftnlen)1);

/*     Compute the closed-loop state matrix. */

#line 543 "SB10ID.f"
    i__1 = *n + *nk;
#line 543 "SB10ID.f"
    dlacpy_("F", n, n, &a[a_offset], lda, &dwork[i7 + 1], &i__1, (ftnlen)1);
#line 544 "SB10ID.f"
    i__1 = *n + *nk;
#line 544 "SB10ID.f"
    dgemm_("N", "N", n, n, m, &c_b5, &dwork[i4 + 1], n, &dwork[i6 + 1], m, &
	    c_b5, &dwork[i7 + 1], &i__1, (ftnlen)1, (ftnlen)1);
#line 546 "SB10ID.f"
    i__1 = *n + *nk;
#line 546 "SB10ID.f"
    dgemm_("N", "N", nk, n, np, &c_b5, &dwork[i2 + 1], nk, &c__[c_offset], 
	    ldc, &c_b6, &dwork[i7 + *n + 1], &i__1, (ftnlen)1, (ftnlen)1);
#line 548 "SB10ID.f"
    i__1 = *n + *nk;
#line 548 "SB10ID.f"
    dgemm_("N", "N", n, nk, m, &c_b5, &dwork[i4 + 1], n, &ck[ck_offset], ldck,
	     &c_b6, &dwork[i7 + (*n + *nk) * *n + 1], &i__1, (ftnlen)1, (
	    ftnlen)1);
#line 550 "SB10ID.f"
    i__1 = *n + *nk;
#line 550 "SB10ID.f"
    dlacpy_("F", nk, nk, &ak[ak_offset], ldak, &dwork[i7 + (*n + *nk) * *n + *
	    n + 1], &i__1, (ftnlen)1);
#line 552 "SB10ID.f"
    i__1 = *n + *nk;
#line 552 "SB10ID.f"
    dgemm_("N", "N", nk, nk, np, &c_b5, &dwork[i2 + 1], nk, &dwork[i5 + 1], 
	    np, &c_b5, &dwork[i7 + (*n + *nk) * *n + *n + 1], &i__1, (ftnlen)
	    1, (ftnlen)1);

/*     Compute the closed-loop poles. */
/*     Additional workspace:  need 3*(N+NK);  prefer larger. */
/*     The fact that M > 0, NP > 0, and NK <= N is used here. */

#line 560 "SB10ID.f"
    i__1 = *n + *nk;
#line 560 "SB10ID.f"
    i__2 = *n + *nk;
#line 560 "SB10ID.f"
    i__3 = *ldwork - iwrk;
#line 560 "SB10ID.f"
    dgees_("N", "N", (L_fp)select_, &i__1, &dwork[i7 + 1], &i__2, &sdim, &
	    dwork[1], &dwork[*n + *nk + 1], &dwork[iwrk + 1], n, &dwork[iwrk 
	    + 1], &i__3, &bwork[1], &info2, (ftnlen)1, (ftnlen)1);
#line 563 "SB10ID.f"
    if (info2 != 0) {
#line 564 "SB10ID.f"
	*info = 3;
#line 565 "SB10ID.f"
	return 0;
#line 566 "SB10ID.f"
    }
#line 567 "SB10ID.f"
    lwa = (integer) dwork[iwrk + 1] + iwrk;
#line 568 "SB10ID.f"
    lwamax = max(lwa,lwamax);

/*     Check the stability of the closed-loop system. */

#line 572 "SB10ID.f"
    ns = 0;
#line 573 "SB10ID.f"
    i__1 = *n + *nk;
#line 573 "SB10ID.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 574 "SB10ID.f"
	if (dwork[i__] >= 0.) {
#line 574 "SB10ID.f"
	    ++ns;
#line 574 "SB10ID.f"
	}
#line 575 "SB10ID.f"
/* L40: */
#line 575 "SB10ID.f"
    }
#line 576 "SB10ID.f"
    if (ns > 0) {
#line 577 "SB10ID.f"
	*info = 6;
#line 578 "SB10ID.f"
	return 0;
#line 579 "SB10ID.f"
    }

#line 581 "SB10ID.f"
    dwork[1] = (doublereal) lwamax;
#line 582 "SB10ID.f"
    return 0;
/* *** Last line of SB10ID *** */
} /* sb10id_ */

