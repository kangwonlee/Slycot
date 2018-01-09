#line 1 "SB10ZD.f"
/* SB10ZD.f -- translated by f2c (version 20100827).
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

#line 1 "SB10ZD.f"
/* Table of constant values */

static doublereal c_b5 = 0.;
static doublereal c_b6 = 1.;
static doublereal c_b42 = -1.;
static integer c__1 = 1;
static integer c__0 = 0;

/* Subroutine */ int sb10zd_(integer *n, integer *m, integer *np, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, doublereal *d__, integer *ldd, doublereal *factor, 
	doublereal *ak, integer *ldak, doublereal *bk, integer *ldbk, 
	doublereal *ck, integer *ldck, doublereal *dk, integer *lddk, 
	doublereal *rcond, doublereal *tol, integer *iwork, doublereal *dwork,
	 integer *ldwork, logical *bwork, integer *info)
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
	    i12, i13, i14, i15, i16, i17, i18, i19, i20, i21, i22, i23, i24, 
	    i25, i26, ns, sdim;
    static doublereal toll;
    static integer iwrk, info2;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal gamma;
    extern /* Subroutine */ int dgees_(char *, char *, L_fp, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen), dgemm_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen), mb02vd_(char *, integer *, integer *, doublereal 
	    *, integer *, integer *, doublereal *, integer *, integer *, 
	    ftnlen), sb02od_(char *, char *, char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen), 
	    mb01rx_(char *, char *, char *, integer *, integer *, doublereal *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), dtrsm_(char *, char *, char *, char *
	    , integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dsyev_(
	    char *, char *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen), dsyrk_(char *
	    , char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen), dlange_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgecon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dgetrf_(integer *, integer *, 
	    doublereal *, integer *, integer *, integer *), dlacpy_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    extern logical select_();
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dgetrs_(
	    char *, integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);
    static integer lwamax;
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int dpotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dsycon_(char *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen);
    static integer minwrk;
    extern /* Subroutine */ int dpotrs_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen), dsytrf_(char *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen), dtrtrs_(
	    char *, char *, char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen), dsytrs_(char *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen);


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

/*     in the Discrete-Time Loop Shaping Design Procedure. */

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
/*             system input/output matrix D of the shaped plant. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= max(1,NP). */

/*     FACTOR  (input) DOUBLE PRECISION */
/*             = 1  implies that an optimal controller is required */
/*                  (not recommended); */
/*             > 1  implies that a suboptimal controller is required */
/*                  achieving a performance FACTOR less than optimal. */
/*             FACTOR >= 1. */

/*     AK      (output) DOUBLE PRECISION array, dimension (LDAK,N) */
/*             The leading N-by-N part of this array contains the */
/*             controller state matrix Ak. */

/*     LDAK    INTEGER */
/*             The leading dimension of the array AK.  LDAK >= max(1,N). */

/*     BK      (output) DOUBLE PRECISION array, dimension (LDBK,NP) */
/*             The leading N-by-NP part of this array contains the */
/*             controller input matrix Bk. */

/*     LDBK    INTEGER */
/*             The leading dimension of the array BK.  LDBK >= max(1,N). */

/*     CK      (output) DOUBLE PRECISION array, dimension (LDCK,N) */
/*             The leading M-by-N part of this array contains the */
/*             controller output matrix Ck. */

/*     LDCK    INTEGER */
/*             The leading dimension of the array CK.  LDCK >= max(1,M). */

/*     DK      (output) DOUBLE PRECISION array, dimension (LDDK,NP) */
/*             The leading M-by-NP part of this array contains the */
/*             controller matrix Dk. */

/*     LDDK    INTEGER */
/*             The leading dimension of the array DK.  LDDK >= max(1,M). */

/*     RCOND   (output) DOUBLE PRECISION array, dimension (6) */
/*             RCOND(1) contains an estimate of the reciprocal condition */
/*                      number of the linear system of equations from */
/*                      which the solution of the P-Riccati equation is */
/*                      obtained; */
/*             RCOND(2) contains an estimate of the reciprocal condition */
/*                      number of the linear system of equations from */
/*                      which the solution of the Q-Riccati equation is */
/*                      obtained; */
/*             RCOND(3) contains an estimate of the reciprocal condition */
/*                      number of the matrix (gamma^2-1)*In - P*Q; */
/*             RCOND(4) contains an estimate of the reciprocal condition */
/*                      number of the matrix Rx + Bx'*X*Bx; */
/*             RCOND(5) contains an estimate of the reciprocal condition */
/*                                                  ^ */
/*                      number of the matrix Ip + D*Dk; */
/*             RCOND(6) contains an estimate of the reciprocal condition */
/*                                                ^ */
/*                      number of the matrix Im + Dk*D. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             Tolerance used for checking the nonsingularity of the */
/*             matrices to be inverted. If TOL <= 0, then a default value */
/*             equal to sqrt(EPS) is used, where EPS is the relative */
/*             machine precision.  TOL < 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension 2*max(N,M+NP) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= 16*N*N + 5*M*M + 7*NP*NP + 6*M*N + 7*M*NP + */
/*                        7*N*NP + 6*N + 2*(M + NP) + */
/*                        max(14*N+23,16*N,2*M-1,2*NP-1). */
/*             For good performance, LDWORK must generally be larger. */

/*     BWORK   LOGICAL array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    (output) INTEGER */
/*             =  0:  successful exit; */
/*             <  0:  if INFO = -i, the i-th argument had an illegal */
/*                    value; */
/*             =  1:  the P-Riccati equation is not solved successfully; */
/*             =  2:  the Q-Riccati equation is not solved successfully; */
/*             =  3:  the iteration to compute eigenvalues or singular */
/*                    values failed to converge; */
/*             =  4:  the matrix (gamma^2-1)*In - P*Q is singular; */
/*             =  5:  the matrix Rx + Bx'*X*Bx is singular; */
/*                                      ^ */
/*             =  6:  the matrix Ip + D*Dk is singular; */
/*                                    ^ */
/*             =  7:  the matrix Im + Dk*D is singular; */
/*             =  8:  the matrix Ip - D*Dk is singular; */
/*             =  9:  the matrix Im - Dk*D is singular; */
/*             = 10:  the closed-loop system is unstable. */

/*     METHOD */

/*     The routine implements the formulas given in [1]. */

/*     REFERENCES */

/*     [1] Gu, D.-W., Petkov, P.H., and Konstantinov, M.M. */
/*         On discrete H-infinity loop shaping design procedure routines. */
/*         Technical Report 00-6, Dept. of Engineering, Univ. of */
/*         Leicester, UK, 2000. */

/*     NUMERICAL ASPECTS */

/*     The accuracy of the results depends on the conditioning of the */
/*     two Riccati equations solved in the controller design. For */
/*     better conditioning it is advised to take FACTOR > 1. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, May 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, July 2001. */

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

#line 254 "SB10ZD.f"
    /* Parameter adjustments */
#line 254 "SB10ZD.f"
    a_dim1 = *lda;
#line 254 "SB10ZD.f"
    a_offset = 1 + a_dim1;
#line 254 "SB10ZD.f"
    a -= a_offset;
#line 254 "SB10ZD.f"
    b_dim1 = *ldb;
#line 254 "SB10ZD.f"
    b_offset = 1 + b_dim1;
#line 254 "SB10ZD.f"
    b -= b_offset;
#line 254 "SB10ZD.f"
    c_dim1 = *ldc;
#line 254 "SB10ZD.f"
    c_offset = 1 + c_dim1;
#line 254 "SB10ZD.f"
    c__ -= c_offset;
#line 254 "SB10ZD.f"
    d_dim1 = *ldd;
#line 254 "SB10ZD.f"
    d_offset = 1 + d_dim1;
#line 254 "SB10ZD.f"
    d__ -= d_offset;
#line 254 "SB10ZD.f"
    ak_dim1 = *ldak;
#line 254 "SB10ZD.f"
    ak_offset = 1 + ak_dim1;
#line 254 "SB10ZD.f"
    ak -= ak_offset;
#line 254 "SB10ZD.f"
    bk_dim1 = *ldbk;
#line 254 "SB10ZD.f"
    bk_offset = 1 + bk_dim1;
#line 254 "SB10ZD.f"
    bk -= bk_offset;
#line 254 "SB10ZD.f"
    ck_dim1 = *ldck;
#line 254 "SB10ZD.f"
    ck_offset = 1 + ck_dim1;
#line 254 "SB10ZD.f"
    ck -= ck_offset;
#line 254 "SB10ZD.f"
    dk_dim1 = *lddk;
#line 254 "SB10ZD.f"
    dk_offset = 1 + dk_dim1;
#line 254 "SB10ZD.f"
    dk -= dk_offset;
#line 254 "SB10ZD.f"
    --rcond;
#line 254 "SB10ZD.f"
    --iwork;
#line 254 "SB10ZD.f"
    --dwork;
#line 254 "SB10ZD.f"
    --bwork;
#line 254 "SB10ZD.f"

#line 254 "SB10ZD.f"
    /* Function Body */
#line 254 "SB10ZD.f"
    *info = 0;
#line 255 "SB10ZD.f"
    if (*n < 0) {
#line 256 "SB10ZD.f"
	*info = -1;
#line 257 "SB10ZD.f"
    } else if (*m < 0) {
#line 258 "SB10ZD.f"
	*info = -2;
#line 259 "SB10ZD.f"
    } else if (*np < 0) {
#line 260 "SB10ZD.f"
	*info = -3;
#line 261 "SB10ZD.f"
    } else if (*lda < max(1,*n)) {
#line 262 "SB10ZD.f"
	*info = -5;
#line 263 "SB10ZD.f"
    } else if (*ldb < max(1,*n)) {
#line 264 "SB10ZD.f"
	*info = -7;
#line 265 "SB10ZD.f"
    } else if (*ldc < max(1,*np)) {
#line 266 "SB10ZD.f"
	*info = -9;
#line 267 "SB10ZD.f"
    } else if (*ldd < max(1,*np)) {
#line 268 "SB10ZD.f"
	*info = -11;
#line 269 "SB10ZD.f"
    } else if (*factor < 1.) {
#line 270 "SB10ZD.f"
	*info = -12;
#line 271 "SB10ZD.f"
    } else if (*ldak < max(1,*n)) {
#line 272 "SB10ZD.f"
	*info = -14;
#line 273 "SB10ZD.f"
    } else if (*ldbk < max(1,*n)) {
#line 274 "SB10ZD.f"
	*info = -16;
#line 275 "SB10ZD.f"
    } else if (*ldck < max(1,*m)) {
#line 276 "SB10ZD.f"
	*info = -18;
#line 277 "SB10ZD.f"
    } else if (*lddk < max(1,*m)) {
#line 278 "SB10ZD.f"
	*info = -20;
#line 279 "SB10ZD.f"
    } else if (*tol >= 1.) {
#line 280 "SB10ZD.f"
	*info = -22;
#line 281 "SB10ZD.f"
    }

/*     Compute workspace. */

/* Computing MAX */
#line 285 "SB10ZD.f"
    i__1 = *n * 14 + 23, i__2 = *n << 4, i__1 = max(i__1,i__2), i__2 = (*m << 
	    1) - 1, i__1 = max(i__1,i__2), i__2 = (*np << 1) - 1;
#line 285 "SB10ZD.f"
    minwrk = (*n << 4) * *n + *m * 5 * *m + *np * 7 * *np + *m * 6 * *n + *m *
	     7 * *np + *n * 7 * *np + *n * 6 + (*m + *np << 1) + max(i__1,
	    i__2);
#line 287 "SB10ZD.f"
    if (*ldwork < minwrk) {
#line 288 "SB10ZD.f"
	*info = -25;
#line 289 "SB10ZD.f"
    }
#line 290 "SB10ZD.f"
    if (*info != 0) {
#line 291 "SB10ZD.f"
	i__1 = -(*info);
#line 291 "SB10ZD.f"
	xerbla_("SB10ZD", &i__1, (ftnlen)6);
#line 292 "SB10ZD.f"
	return 0;
#line 293 "SB10ZD.f"
    }

/*     Quick return if possible. */
/*     Note that some computation could be made if one or two of the */
/*     dimension parameters N, M, and P are zero, but the results are */
/*     not so meaningful. */

#line 300 "SB10ZD.f"
    if (*n == 0 || *m == 0 || *np == 0) {
#line 301 "SB10ZD.f"
	rcond[1] = 1.;
#line 302 "SB10ZD.f"
	rcond[2] = 1.;
#line 303 "SB10ZD.f"
	rcond[3] = 1.;
#line 304 "SB10ZD.f"
	rcond[4] = 1.;
#line 305 "SB10ZD.f"
	rcond[5] = 1.;
#line 306 "SB10ZD.f"
	rcond[6] = 1.;
#line 307 "SB10ZD.f"
	dwork[1] = 1.;
#line 308 "SB10ZD.f"
	return 0;
#line 309 "SB10ZD.f"
    }

/*     Set the default tolerance, if needed. */

#line 313 "SB10ZD.f"
    if (*tol <= 0.) {
#line 314 "SB10ZD.f"
	toll = sqrt(dlamch_("Epsilon", (ftnlen)7));
#line 315 "SB10ZD.f"
    } else {
#line 316 "SB10ZD.f"
	toll = *tol;
#line 317 "SB10ZD.f"
    }

/*     Workspace usage. */

#line 321 "SB10ZD.f"
    n2 = *n << 1;
#line 322 "SB10ZD.f"
    i1 = *n * *n + 1;
#line 323 "SB10ZD.f"
    i2 = i1 + *n * *n;
#line 324 "SB10ZD.f"
    i3 = i2 + *np * *np;
#line 325 "SB10ZD.f"
    i4 = i3 + *m * *m;
#line 326 "SB10ZD.f"
    i5 = i4 + *np * *np;
#line 327 "SB10ZD.f"
    i6 = i5 + *m * *m;
#line 328 "SB10ZD.f"
    i7 = i6 + *m * *n;
#line 329 "SB10ZD.f"
    i8 = i7 + *m * *n;
#line 330 "SB10ZD.f"
    i9 = i8 + *n * *n;
#line 331 "SB10ZD.f"
    i10 = i9 + *n * *n;
#line 332 "SB10ZD.f"
    i11 = i10 + n2;
#line 333 "SB10ZD.f"
    i12 = i11 + n2;
#line 334 "SB10ZD.f"
    i13 = i12 + n2;
#line 335 "SB10ZD.f"
    i14 = i13 + n2 * n2;
#line 336 "SB10ZD.f"
    i15 = i14 + n2 * n2;

#line 338 "SB10ZD.f"
    iwrk = i15 + n2 * n2;
#line 339 "SB10ZD.f"
    lwamax = 0;

/*     Compute R1 = Ip + D*D' . */

#line 343 "SB10ZD.f"
    dlaset_("U", np, np, &c_b5, &c_b6, &dwork[i2], np, (ftnlen)1);
#line 344 "SB10ZD.f"
    dsyrk_("U", "N", np, m, &c_b6, &d__[d_offset], ldd, &c_b6, &dwork[i2], np,
	     (ftnlen)1, (ftnlen)1);
#line 345 "SB10ZD.f"
    dlacpy_("U", np, np, &dwork[i2], np, &dwork[i4], np, (ftnlen)1);

/*     Factorize R1 = R'*R . */

#line 349 "SB10ZD.f"
    dpotrf_("U", np, &dwork[i4], np, &info2, (ftnlen)1);
/*                 -1 */
/*     Compute C'*R   in BK . */

#line 353 "SB10ZD.f"
    ma02ad_("F", np, n, &c__[c_offset], ldc, &bk[bk_offset], ldbk, (ftnlen)1);
#line 354 "SB10ZD.f"
    dtrsm_("R", "U", "N", "N", n, np, &c_b6, &dwork[i4], np, &bk[bk_offset], 
	    ldbk, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Compute R2 = Im + D'*D . */

#line 359 "SB10ZD.f"
    dlaset_("U", m, m, &c_b5, &c_b6, &dwork[i3], m, (ftnlen)1);
#line 360 "SB10ZD.f"
    dsyrk_("U", "T", m, np, &c_b6, &d__[d_offset], ldd, &c_b6, &dwork[i3], m, 
	    (ftnlen)1, (ftnlen)1);
#line 361 "SB10ZD.f"
    dlacpy_("U", m, m, &dwork[i3], m, &dwork[i5], m, (ftnlen)1);

/*     Factorize R2 = U'*U . */

#line 365 "SB10ZD.f"
    dpotrf_("U", m, &dwork[i5], m, &info2, (ftnlen)1);
/*               -1 */
/*     Compute (U  )'*B' . */

#line 369 "SB10ZD.f"
    ma02ad_("F", n, m, &b[b_offset], ldb, &dwork[i6], m, (ftnlen)1);
#line 370 "SB10ZD.f"
    dtrtrs_("U", "T", "N", m, n, &dwork[i5], m, &dwork[i6], m, &info2, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Compute D'*C . */

#line 375 "SB10ZD.f"
    dgemm_("T", "N", m, n, np, &c_b6, &d__[d_offset], ldd, &c__[c_offset], 
	    ldc, &c_b5, &dwork[i7], m, (ftnlen)1, (ftnlen)1);
/*               -1 */
/*     Compute (U  )'*D'*C . */

#line 380 "SB10ZD.f"
    dtrtrs_("U", "T", "N", m, n, &dwork[i5], m, &dwork[i7], m, &info2, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);
/*                          -1 */
/*     Compute Ar = A - B*R2  D'*C . */

#line 385 "SB10ZD.f"
    dlacpy_("F", n, n, &a[a_offset], lda, &dwork[i8], n, (ftnlen)1);
#line 386 "SB10ZD.f"
    dgemm_("T", "N", n, n, m, &c_b42, &dwork[i6], m, &dwork[i7], m, &c_b6, &
	    dwork[i8], n, (ftnlen)1, (ftnlen)1);
/*                       -1 */
/*     Compute Cr = C'*R1  *C . */

#line 391 "SB10ZD.f"
    dsyrk_("U", "N", n, np, &c_b6, &bk[bk_offset], ldbk, &c_b5, &dwork[i9], n,
	     (ftnlen)1, (ftnlen)1);
/*                      -1 */
/*     Compute Dr = B*R2  B' in AK . */

#line 395 "SB10ZD.f"
    dsyrk_("U", "T", n, m, &c_b6, &dwork[i6], m, &c_b5, &ak[ak_offset], ldak, 
	    (ftnlen)1, (ftnlen)1);
/*                                                       -1 */
/*     Solution of the Riccati equation Ar'*P*(In + Dr*P)  Ar - P + */
/*                                              Cr = 0 . */
#line 399 "SB10ZD.f"
    i__1 = *ldwork - iwrk + 1;
#line 399 "SB10ZD.f"
    sb02od_("D", "G", "N", "U", "Z", "S", n, m, np, &dwork[i8], n, &ak[
	    ak_offset], ldak, &dwork[i9], n, &dwork[1], m, &dwork[1], n, &
	    rcond[1], &dwork[1], n, &dwork[i10], &dwork[i11], &dwork[i12], &
	    dwork[i13], &n2, &dwork[i14], &n2, &dwork[i15], &n2, &c_b42, &
	    iwork[1], &dwork[iwrk], &i__1, &bwork[1], &info2, (ftnlen)1, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 405 "SB10ZD.f"
    if (info2 != 0) {
#line 406 "SB10ZD.f"
	*info = 1;
#line 407 "SB10ZD.f"
	return 0;
#line 408 "SB10ZD.f"
    }
/* Computing MAX */
#line 409 "SB10ZD.f"
    i__1 = lwamax, i__2 = (integer) dwork[iwrk] + iwrk - 1;
#line 409 "SB10ZD.f"
    lwamax = max(i__1,i__2);

/*     Transpose Ar . */

#line 413 "SB10ZD.f"
    i__1 = *n - 1;
#line 413 "SB10ZD.f"
    for (j = 1; j <= i__1; ++j) {
#line 414 "SB10ZD.f"
	dswap_(&j, &dwork[i8 + j], n, &dwork[i8 + j * *n], &c__1);
#line 415 "SB10ZD.f"
/* L10: */
#line 415 "SB10ZD.f"
    }
/*                                                      -1 */
/*     Solution of the Riccati equation Ar*Q*(In + Cr*Q)  *Ar' - Q + */
/*                                             Dr = 0 . */
#line 419 "SB10ZD.f"
    i__1 = *ldwork - iwrk + 1;
#line 419 "SB10ZD.f"
    sb02od_("D", "G", "N", "U", "Z", "S", n, m, np, &dwork[i8], n, &dwork[i9],
	     n, &ak[ak_offset], ldak, &dwork[1], m, &dwork[1], n, &rcond[2], &
	    dwork[i1], n, &dwork[i10], &dwork[i11], &dwork[i12], &dwork[i13], 
	    &n2, &dwork[i14], &n2, &dwork[i15], &n2, &c_b42, &iwork[1], &
	    dwork[iwrk], &i__1, &bwork[1], &info2, (ftnlen)1, (ftnlen)1, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 425 "SB10ZD.f"
    if (info2 != 0) {
#line 426 "SB10ZD.f"
	*info = 2;
#line 427 "SB10ZD.f"
	return 0;
#line 428 "SB10ZD.f"
    }
/* Computing MAX */
#line 429 "SB10ZD.f"
    i__1 = lwamax, i__2 = (integer) dwork[iwrk] + iwrk - 1;
#line 429 "SB10ZD.f"
    lwamax = max(i__1,i__2);

/*     Compute gamma. */

#line 433 "SB10ZD.f"
    dgemm_("N", "N", n, n, n, &c_b6, &dwork[i1], n, &dwork[1], n, &c_b5, &
	    dwork[i8], n, (ftnlen)1, (ftnlen)1);
#line 435 "SB10ZD.f"
    i__1 = *ldwork - iwrk + 1;
#line 435 "SB10ZD.f"
    dgees_("N", "N", (L_fp)select_, n, &dwork[i8], n, &sdim, &dwork[i10], &
	    dwork[i11], &dwork[iwrk], n, &dwork[iwrk], &i__1, &bwork[1], &
	    info2, (ftnlen)1, (ftnlen)1);
#line 438 "SB10ZD.f"
    if (info2 != 0) {
#line 439 "SB10ZD.f"
	*info = 3;
#line 440 "SB10ZD.f"
	return 0;
#line 441 "SB10ZD.f"
    }
/* Computing MAX */
#line 442 "SB10ZD.f"
    i__1 = lwamax, i__2 = (integer) dwork[iwrk] + iwrk - 1;
#line 442 "SB10ZD.f"
    lwamax = max(i__1,i__2);
#line 443 "SB10ZD.f"
    gamma = 0.;

#line 445 "SB10ZD.f"
    i__1 = *n - 1;
#line 445 "SB10ZD.f"
    for (i__ = 0; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 446 "SB10ZD.f"
	d__1 = gamma, d__2 = dwork[i10 + i__];
#line 446 "SB10ZD.f"
	gamma = max(d__1,d__2);
#line 447 "SB10ZD.f"
/* L20: */
#line 447 "SB10ZD.f"
    }

#line 449 "SB10ZD.f"
    gamma = *factor * sqrt(gamma + 1.);

/*     Workspace usage. */

#line 453 "SB10ZD.f"
    i5 = i4 + *np * *np;
#line 454 "SB10ZD.f"
    i6 = i5 + *m * *m;
#line 455 "SB10ZD.f"
    i7 = i6 + *np * *np;
#line 456 "SB10ZD.f"
    i8 = i7 + *np * *np;
#line 457 "SB10ZD.f"
    i9 = i8 + *np * *np;
#line 458 "SB10ZD.f"
    i10 = i9 + *np;
#line 459 "SB10ZD.f"
    i11 = i10 + *np * *np;
#line 460 "SB10ZD.f"
    i12 = i11 + *m * *m;
#line 461 "SB10ZD.f"
    i13 = i12 + *m;

#line 463 "SB10ZD.f"
    iwrk = i13 + *m * *m;

/*     Compute the eigenvalues and eigenvectors of R1 . */

#line 467 "SB10ZD.f"
    dlacpy_("U", np, np, &dwork[i2], np, &dwork[i8], np, (ftnlen)1);
#line 468 "SB10ZD.f"
    i__1 = *ldwork - iwrk + 1;
#line 468 "SB10ZD.f"
    dsyev_("V", "U", np, &dwork[i8], np, &dwork[i9], &dwork[iwrk], &i__1, &
	    info2, (ftnlen)1, (ftnlen)1);
#line 470 "SB10ZD.f"
    if (info2 != 0) {
#line 471 "SB10ZD.f"
	*info = 3;
#line 472 "SB10ZD.f"
	return 0;
#line 473 "SB10ZD.f"
    }
/* Computing MAX */
#line 474 "SB10ZD.f"
    i__1 = lwamax, i__2 = (integer) dwork[iwrk] + iwrk - 1;
#line 474 "SB10ZD.f"
    lwamax = max(i__1,i__2);
/*               -1/2 */
/*     Compute R1     . */

#line 478 "SB10ZD.f"
    i__1 = *np;
#line 478 "SB10ZD.f"
    for (j = 1; j <= i__1; ++j) {
#line 479 "SB10ZD.f"
	i__2 = *np;
#line 479 "SB10ZD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 480 "SB10ZD.f"
	    dwork[i10 - 1 + i__ + (j - 1) * *np] = dwork[i8 - 1 + j + (i__ - 
		    1) * *np] / sqrt(dwork[i9 + i__ - 1]);
#line 482 "SB10ZD.f"
/* L30: */
#line 482 "SB10ZD.f"
	}
#line 483 "SB10ZD.f"
/* L40: */
#line 483 "SB10ZD.f"
    }

#line 485 "SB10ZD.f"
    dgemm_("N", "N", np, np, np, &c_b6, &dwork[i8], np, &dwork[i10], np, &
	    c_b5, &dwork[i4], np, (ftnlen)1, (ftnlen)1);

/*     Compute the eigenvalues and eigenvectors of R2 . */

#line 490 "SB10ZD.f"
    dlacpy_("U", m, m, &dwork[i3], m, &dwork[i11], m, (ftnlen)1);
#line 491 "SB10ZD.f"
    i__1 = *ldwork - iwrk + 1;
#line 491 "SB10ZD.f"
    dsyev_("V", "U", m, &dwork[i11], m, &dwork[i12], &dwork[iwrk], &i__1, &
	    info2, (ftnlen)1, (ftnlen)1);
#line 493 "SB10ZD.f"
    if (info2 != 0) {
#line 494 "SB10ZD.f"
	*info = 3;
#line 495 "SB10ZD.f"
	return 0;
#line 496 "SB10ZD.f"
    }
/* Computing MAX */
#line 497 "SB10ZD.f"
    i__1 = lwamax, i__2 = (integer) dwork[iwrk] + iwrk - 1;
#line 497 "SB10ZD.f"
    lwamax = max(i__1,i__2);
/*               -1/2 */
/*     Compute R2     . */

#line 501 "SB10ZD.f"
    i__1 = *m;
#line 501 "SB10ZD.f"
    for (j = 1; j <= i__1; ++j) {
#line 502 "SB10ZD.f"
	i__2 = *m;
#line 502 "SB10ZD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 503 "SB10ZD.f"
	    dwork[i13 - 1 + i__ + (j - 1) * *m] = dwork[i11 - 1 + j + (i__ - 
		    1) * *m] / sqrt(dwork[i12 + i__ - 1]);
#line 505 "SB10ZD.f"
/* L50: */
#line 505 "SB10ZD.f"
	}
#line 506 "SB10ZD.f"
/* L60: */
#line 506 "SB10ZD.f"
    }

#line 508 "SB10ZD.f"
    dgemm_("N", "N", m, m, m, &c_b6, &dwork[i11], m, &dwork[i13], m, &c_b5, &
	    dwork[i5], m, (ftnlen)1, (ftnlen)1);

/*     Compute R1 + C*Q*C' . */

#line 513 "SB10ZD.f"
    dgemm_("N", "T", n, np, n, &c_b6, &dwork[i1], n, &c__[c_offset], ldc, &
	    c_b5, &bk[bk_offset], ldbk, (ftnlen)1, (ftnlen)1);
#line 515 "SB10ZD.f"
    mb01rx_("L", "U", "N", np, n, &c_b6, &c_b6, &dwork[i2], np, &c__[c_offset]
	    , ldc, &bk[bk_offset], ldbk, &info2, (ftnlen)1, (ftnlen)1, (
	    ftnlen)1);
#line 517 "SB10ZD.f"
    dlacpy_("U", np, np, &dwork[i2], np, &dwork[i8], np, (ftnlen)1);

/*     Compute the eigenvalues and eigenvectors of R1 + C*Q*C' . */

#line 521 "SB10ZD.f"
    i__1 = *ldwork - iwrk + 1;
#line 521 "SB10ZD.f"
    dsyev_("V", "U", np, &dwork[i8], np, &dwork[i9], &dwork[iwrk], &i__1, &
	    info2, (ftnlen)1, (ftnlen)1);
#line 523 "SB10ZD.f"
    if (info2 != 0) {
#line 524 "SB10ZD.f"
	*info = 3;
#line 525 "SB10ZD.f"
	return 0;
#line 526 "SB10ZD.f"
    }
/* Computing MAX */
#line 527 "SB10ZD.f"
    i__1 = lwamax, i__2 = (integer) dwork[iwrk] + iwrk - 1;
#line 527 "SB10ZD.f"
    lwamax = max(i__1,i__2);
/*                            -1 */
/*     Compute ( R1 + C*Q*C' )   . */

#line 531 "SB10ZD.f"
    i__1 = *np;
#line 531 "SB10ZD.f"
    for (j = 1; j <= i__1; ++j) {
#line 532 "SB10ZD.f"
	i__2 = *np;
#line 532 "SB10ZD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 533 "SB10ZD.f"
	    dwork[i10 - 1 + i__ + (j - 1) * *np] = dwork[i8 - 1 + j + (i__ - 
		    1) * *np] / dwork[i9 + i__ - 1];
#line 535 "SB10ZD.f"
/* L70: */
#line 535 "SB10ZD.f"
	}
#line 536 "SB10ZD.f"
/* L80: */
#line 536 "SB10ZD.f"
    }

#line 538 "SB10ZD.f"
    dgemm_("N", "N", np, np, np, &c_b6, &dwork[i8], np, &dwork[i10], np, &
	    c_b5, &dwork[i6], np, (ftnlen)1, (ftnlen)1);
/*               -1 */
/*     Compute Z2  . */

#line 543 "SB10ZD.f"
    i__1 = *np;
#line 543 "SB10ZD.f"
    for (j = 1; j <= i__1; ++j) {
#line 544 "SB10ZD.f"
	i__2 = *np;
#line 544 "SB10ZD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 545 "SB10ZD.f"
	    dwork[i10 - 1 + i__ + (j - 1) * *np] = dwork[i8 - 1 + j + (i__ - 
		    1) * *np] * sqrt(dwork[i9 + i__ - 1]);
#line 547 "SB10ZD.f"
/* L90: */
#line 547 "SB10ZD.f"
	}
#line 548 "SB10ZD.f"
/* L100: */
#line 548 "SB10ZD.f"
    }

#line 550 "SB10ZD.f"
    dgemm_("N", "N", np, np, np, &c_b6, &dwork[i8], np, &dwork[i10], np, &
	    c_b5, &dwork[i7], np, (ftnlen)1, (ftnlen)1);

/*     Workspace usage. */

#line 555 "SB10ZD.f"
    i9 = i8 + *n * *np;
#line 556 "SB10ZD.f"
    i10 = i9 + *n * *np;
#line 557 "SB10ZD.f"
    i11 = i10 + *np * *m;
#line 558 "SB10ZD.f"
    i12 = i11 + (*np + *m) * (*np + *m);
#line 559 "SB10ZD.f"
    i13 = i12 + *n * (*np + *m);
#line 560 "SB10ZD.f"
    i14 = i13 + *n * (*np + *m);
#line 561 "SB10ZD.f"
    i15 = i14 + *n * *n;
#line 562 "SB10ZD.f"
    i16 = i15 + *n * *n;
#line 563 "SB10ZD.f"
    i17 = i16 + (*np + *m) * *n;
#line 564 "SB10ZD.f"
    i18 = i17 + (*np + *m) * (*np + *m);
#line 565 "SB10ZD.f"
    i19 = i18 + (*np + *m) * *n;
#line 566 "SB10ZD.f"
    i20 = i19 + *m * *n;
#line 567 "SB10ZD.f"
    i21 = i20 + *m * *np;
#line 568 "SB10ZD.f"
    i22 = i21 + *np * *n;
#line 569 "SB10ZD.f"
    i23 = i22 + *n * *n;
#line 570 "SB10ZD.f"
    i24 = i23 + *n * *np;
#line 571 "SB10ZD.f"
    i25 = i24 + *np * *np;
#line 572 "SB10ZD.f"
    i26 = i25 + *m * *m;

#line 574 "SB10ZD.f"
    iwrk = i26 + *n * *m;

/*     Compute A*Q*C' + B*D' . */

#line 578 "SB10ZD.f"
    dgemm_("N", "T", n, np, m, &c_b6, &b[b_offset], ldb, &d__[d_offset], ldd, 
	    &c_b5, &dwork[i8], n, (ftnlen)1, (ftnlen)1);
#line 580 "SB10ZD.f"
    dgemm_("N", "N", n, np, n, &c_b6, &a[a_offset], lda, &bk[bk_offset], ldbk,
	     &c_b6, &dwork[i8], n, (ftnlen)1, (ftnlen)1);
/*                                                 -1 */
/*     Compute H = -( A*Q*C'+B*D' )*( R1 + C*Q*C' )   . */

#line 585 "SB10ZD.f"
    dgemm_("N", "N", n, np, np, &c_b42, &dwork[i8], n, &dwork[i6], np, &c_b5, 
	    &dwork[i9], n, (ftnlen)1, (ftnlen)1);
/*               -1/2 */
/*     Compute R1    D . */

#line 590 "SB10ZD.f"
    dgemm_("N", "N", np, m, np, &c_b6, &dwork[i4], np, &d__[d_offset], ldd, &
	    c_b5, &dwork[i10], np, (ftnlen)1, (ftnlen)1);

/*     Compute Rx . */

#line 595 "SB10ZD.f"
    i__1 = *np;
#line 595 "SB10ZD.f"
    for (j = 1; j <= i__1; ++j) {
#line 596 "SB10ZD.f"
	dcopy_(&j, &dwork[i2 + (j - 1) * *np], &c__1, &dwork[i11 + (j - 1) * (
		*np + *m)], &c__1);
#line 598 "SB10ZD.f"
	dwork[i11 - 1 + j + (j - 1) * (*np + *m)] = dwork[i2 - 1 + j + (j - 1)
		 * *np] - gamma * gamma;
#line 600 "SB10ZD.f"
/* L110: */
#line 600 "SB10ZD.f"
    }

#line 602 "SB10ZD.f"
    i__1 = *np + *m;
#line 602 "SB10ZD.f"
    dgemm_("N", "N", np, m, np, &c_b6, &dwork[i7], np, &dwork[i10], np, &c_b5,
	     &dwork[i11 + (*np + *m) * *np], &i__1, (ftnlen)1, (ftnlen)1);
#line 605 "SB10ZD.f"
    i__1 = *np + *m;
#line 605 "SB10ZD.f"
    dlaset_("U", m, m, &c_b5, &c_b6, &dwork[i11 + (*np + *m) * *np + *np], &
	    i__1, (ftnlen)1);

/*     Compute Bx . */

#line 610 "SB10ZD.f"
    dgemm_("N", "N", n, np, np, &c_b42, &dwork[i9], n, &dwork[i7], np, &c_b5, 
	    &dwork[i12], n, (ftnlen)1, (ftnlen)1);
#line 612 "SB10ZD.f"
    dgemm_("N", "N", n, m, m, &c_b6, &b[b_offset], ldb, &dwork[i5], m, &c_b5, 
	    &dwork[i12 + *n * *np], n, (ftnlen)1, (ftnlen)1);

/*     Compute Sx . */

#line 617 "SB10ZD.f"
    dgemm_("T", "N", n, np, np, &c_b6, &c__[c_offset], ldc, &dwork[i7], np, &
	    c_b5, &dwork[i13], n, (ftnlen)1, (ftnlen)1);
#line 619 "SB10ZD.f"
    dgemm_("T", "N", n, m, np, &c_b6, &c__[c_offset], ldc, &dwork[i10], np, &
	    c_b5, &dwork[i13 + *n * *np], n, (ftnlen)1, (ftnlen)1);

/*     Compute  (gamma^2 - 1)*In - P*Q . */

#line 624 "SB10ZD.f"
    d__1 = gamma * gamma - 1.;
#line 624 "SB10ZD.f"
    dlaset_("F", n, n, &c_b5, &d__1, &dwork[i14], n, (ftnlen)1);
#line 625 "SB10ZD.f"
    dgemm_("N", "N", n, n, n, &c_b42, &dwork[1], n, &dwork[i1], n, &c_b6, &
	    dwork[i14], n, (ftnlen)1, (ftnlen)1);
/*                                          -1 */
/*     Compute X =  ((gamma^2 - 1)*In - P*Q)  *gamma^2*P . */

#line 630 "SB10ZD.f"
    dlacpy_("F", n, n, &dwork[1], n, &dwork[i15], n, (ftnlen)1);
#line 631 "SB10ZD.f"
    d__1 = gamma * gamma;
#line 631 "SB10ZD.f"
    dlascl_("G", &c__0, &c__0, &c_b6, &d__1, n, n, &dwork[i15], n, info, (
	    ftnlen)1);
#line 633 "SB10ZD.f"
    anorm = dlange_("1", n, n, &dwork[i14], n, &dwork[iwrk], (ftnlen)1);
#line 634 "SB10ZD.f"
    dgetrf_(n, n, &dwork[i14], n, &iwork[1], &info2);
#line 635 "SB10ZD.f"
    if (info2 > 0) {
#line 636 "SB10ZD.f"
	*info = 4;
#line 637 "SB10ZD.f"
	return 0;
#line 638 "SB10ZD.f"
    }
#line 639 "SB10ZD.f"
    dgecon_("1", n, &dwork[i14], n, &anorm, &rcond[3], &dwork[iwrk], &iwork[*
	    n + 1], &info2, (ftnlen)1);

/*     Return if the matrix is singular to working precision. */

#line 644 "SB10ZD.f"
    if (rcond[3] < toll) {
#line 645 "SB10ZD.f"
	*info = 4;
#line 646 "SB10ZD.f"
	return 0;
#line 647 "SB10ZD.f"
    }
#line 648 "SB10ZD.f"
    dgetrs_("N", n, n, &dwork[i14], n, &iwork[1], &dwork[i15], n, &info2, (
	    ftnlen)1);

/*     Compute Bx'*X . */

#line 653 "SB10ZD.f"
    i__1 = *np + *m;
#line 653 "SB10ZD.f"
    i__2 = *np + *m;
#line 653 "SB10ZD.f"
    dgemm_("T", "N", &i__1, n, n, &c_b6, &dwork[i12], n, &dwork[i15], n, &
	    c_b5, &dwork[i16], &i__2, (ftnlen)1, (ftnlen)1);

/*     Compute Rx + Bx'*X*Bx . */

#line 658 "SB10ZD.f"
    i__1 = *np + *m;
#line 658 "SB10ZD.f"
    i__2 = *np + *m;
#line 658 "SB10ZD.f"
    i__3 = *np + *m;
#line 658 "SB10ZD.f"
    i__4 = *np + *m;
#line 658 "SB10ZD.f"
    dlacpy_("U", &i__1, &i__2, &dwork[i11], &i__3, &dwork[i17], &i__4, (
	    ftnlen)1);
#line 660 "SB10ZD.f"
    i__1 = *np + *m;
#line 660 "SB10ZD.f"
    i__2 = *np + *m;
#line 660 "SB10ZD.f"
    i__3 = *np + *m;
#line 660 "SB10ZD.f"
    mb01rx_("L", "U", "N", &i__1, n, &c_b6, &c_b6, &dwork[i17], &i__2, &dwork[
	    i16], &i__3, &dwork[i12], n, &info2, (ftnlen)1, (ftnlen)1, (
	    ftnlen)1);

/*     Compute  -( Sx' + Bx'*X*A ) . */

#line 665 "SB10ZD.f"
    i__1 = *np + *m;
#line 665 "SB10ZD.f"
    i__2 = *np + *m;
#line 665 "SB10ZD.f"
    ma02ad_("F", n, &i__1, &dwork[i13], n, &dwork[i18], &i__2, (ftnlen)1);
#line 666 "SB10ZD.f"
    i__1 = *np + *m;
#line 666 "SB10ZD.f"
    i__2 = *np + *m;
#line 666 "SB10ZD.f"
    i__3 = *np + *m;
#line 666 "SB10ZD.f"
    dgemm_("N", "N", &i__1, n, n, &c_b42, &dwork[i16], &i__2, &a[a_offset], 
	    lda, &c_b42, &dwork[i18], &i__3, (ftnlen)1, (ftnlen)1);

/*     Factorize Rx + Bx'*X*Bx . */

#line 671 "SB10ZD.f"
    i__1 = *np + *m;
#line 671 "SB10ZD.f"
    i__2 = *np + *m;
#line 671 "SB10ZD.f"
    anorm = dlansy_("1", "U", &i__1, &dwork[i17], &i__2, &dwork[iwrk], (
	    ftnlen)1, (ftnlen)1);
#line 673 "SB10ZD.f"
    i__1 = *np + *m;
#line 673 "SB10ZD.f"
    i__2 = *np + *m;
#line 673 "SB10ZD.f"
    i__3 = *ldwork - iwrk;
#line 673 "SB10ZD.f"
    dsytrf_("U", &i__1, &dwork[i17], &i__2, &iwork[1], &dwork[iwrk], &i__3, &
	    info2, (ftnlen)1);
#line 675 "SB10ZD.f"
    if (info2 != 0) {
#line 676 "SB10ZD.f"
	*info = 5;
#line 677 "SB10ZD.f"
	return 0;
#line 678 "SB10ZD.f"
    }
#line 679 "SB10ZD.f"
    i__1 = *np + *m;
#line 679 "SB10ZD.f"
    i__2 = *np + *m;
#line 679 "SB10ZD.f"
    dsycon_("U", &i__1, &dwork[i17], &i__2, &iwork[1], &anorm, &rcond[4], &
	    dwork[iwrk], &iwork[*np + *m + 1], &info2, (ftnlen)1);

/*     Return if the matrix is singular to working precision. */

#line 684 "SB10ZD.f"
    if (rcond[4] < toll) {
#line 685 "SB10ZD.f"
	*info = 5;
#line 686 "SB10ZD.f"
	return 0;
#line 687 "SB10ZD.f"
    }
/*                                   -1 */
/*     Compute F = -( Rx + Bx'*X*Bx )  ( Sx' + Bx'*X*A ) . */

#line 691 "SB10ZD.f"
    i__1 = *np + *m;
#line 691 "SB10ZD.f"
    i__2 = *np + *m;
#line 691 "SB10ZD.f"
    i__3 = *np + *m;
#line 691 "SB10ZD.f"
    dsytrs_("U", &i__1, n, &dwork[i17], &i__2, &iwork[1], &dwork[i18], &i__3, 
	    &info2, (ftnlen)1);

/*     Compute B'*X . */

#line 696 "SB10ZD.f"
    dgemm_("T", "N", m, n, n, &c_b6, &b[b_offset], ldb, &dwork[i15], n, &c_b5,
	     &dwork[i19], m, (ftnlen)1, (ftnlen)1);

/*     Compute  -( D' - B'*X*H ) . */

#line 701 "SB10ZD.f"
    i__1 = *np;
#line 701 "SB10ZD.f"
    for (j = 1; j <= i__1; ++j) {
#line 702 "SB10ZD.f"
	i__2 = *m;
#line 702 "SB10ZD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 703 "SB10ZD.f"
	    dwork[i20 - 1 + i__ + (j - 1) * *m] = -d__[j + i__ * d_dim1];
#line 704 "SB10ZD.f"
/* L120: */
#line 704 "SB10ZD.f"
	}
#line 705 "SB10ZD.f"
/* L130: */
#line 705 "SB10ZD.f"
    }

#line 707 "SB10ZD.f"
    dgemm_("N", "N", m, np, n, &c_b6, &dwork[i19], m, &dwork[i9], n, &c_b6, &
	    dwork[i20], m, (ftnlen)1, (ftnlen)1);
/*                   -1 */
/*     Compute C + Z2  *F1 . */

#line 712 "SB10ZD.f"
    dlacpy_("F", np, n, &c__[c_offset], ldc, &dwork[i21], np, (ftnlen)1);
#line 713 "SB10ZD.f"
    i__1 = *np + *m;
#line 713 "SB10ZD.f"
    dgemm_("N", "N", np, n, np, &c_b6, &dwork[i7], np, &dwork[i18], &i__1, &
	    c_b6, &dwork[i21], np, (ftnlen)1, (ftnlen)1);

/*     Compute R2 + B'*X*B . */

#line 718 "SB10ZD.f"
    mb01rx_("L", "U", "N", m, n, &c_b6, &c_b6, &dwork[i3], m, &dwork[i19], m, 
	    &b[b_offset], ldb, &info2, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Factorize R2 + B'*X*B . */

#line 723 "SB10ZD.f"
    dpotrf_("U", m, &dwork[i3], m, &info2, (ftnlen)1);
/*             ^                    -1 */
/*     Compute Dk = -( R2 + B'*X*B )  (D' - B'*X*H) . */

#line 727 "SB10ZD.f"
    dlacpy_("F", m, np, &dwork[i20], m, &dk[dk_offset], lddk, (ftnlen)1);
#line 728 "SB10ZD.f"
    dpotrs_("U", m, np, &dwork[i3], m, &dk[dk_offset], lddk, &info2, (ftnlen)
	    1);
/*             ^           ^ */
/*     Compute Bk = -H + B*Dk . */

#line 732 "SB10ZD.f"
    dlacpy_("F", n, np, &dwork[i9], n, &dwork[i23], n, (ftnlen)1);
#line 733 "SB10ZD.f"
    dgemm_("N", "N", n, np, m, &c_b6, &b[b_offset], ldb, &dk[dk_offset], lddk,
	     &c_b42, &dwork[i23], n, (ftnlen)1, (ftnlen)1);
/*               -1/2 */
/*     Compute R2    *F2  . */

#line 738 "SB10ZD.f"
    i__1 = *np + *m;
#line 738 "SB10ZD.f"
    dgemm_("N", "N", m, n, m, &c_b6, &dwork[i5], m, &dwork[i18 + *np], &i__1, 
	    &c_b5, &ck[ck_offset], ldck, (ftnlen)1, (ftnlen)1);
/*             ^      -1/2      ^          -1 */
/*     Compute Ck = R2    *F2 - Dk*( C + Z2  *F1 ) . */

#line 743 "SB10ZD.f"
    dgemm_("N", "N", m, n, np, &c_b42, &dk[dk_offset], lddk, &dwork[i21], np, 
	    &c_b6, &ck[ck_offset], ldck, (ftnlen)1, (ftnlen)1);
/*             ^                ^ */
/*     Compute Ak = A + H*C + B*Ck . */

#line 748 "SB10ZD.f"
    dlacpy_("F", n, n, &a[a_offset], lda, &ak[ak_offset], ldak, (ftnlen)1);
#line 749 "SB10ZD.f"
    dgemm_("N", "N", n, n, np, &c_b6, &dwork[i9], n, &c__[c_offset], ldc, &
	    c_b6, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);
#line 751 "SB10ZD.f"
    dgemm_("N", "N", n, n, m, &c_b6, &b[b_offset], ldb, &ck[ck_offset], ldck, 
	    &c_b6, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);
/*                    ^ */
/*     Compute Ip + D*Dk . */

#line 756 "SB10ZD.f"
    dlaset_("Full", np, np, &c_b5, &c_b6, &dwork[i24], np, (ftnlen)4);
#line 757 "SB10ZD.f"
    dgemm_("N", "N", np, np, m, &c_b6, &d__[d_offset], ldd, &dk[dk_offset], 
	    lddk, &c_b6, &dwork[i24], np, (ftnlen)1, (ftnlen)1);
/*                   ^ */
/*     Compute  Im + Dk*D . */

#line 762 "SB10ZD.f"
    dlaset_("Full", m, m, &c_b5, &c_b6, &dwork[i25], m, (ftnlen)4);
#line 763 "SB10ZD.f"
    dgemm_("N", "N", m, m, np, &c_b6, &dk[dk_offset], lddk, &d__[d_offset], 
	    ldd, &c_b6, &dwork[i25], m, (ftnlen)1, (ftnlen)1);
/*                  ^ ^    ^         ^    -1 */
/*     Compute Ck = M*Ck,  M = (Im + Dk*D)   . */

#line 768 "SB10ZD.f"
    anorm = dlange_("1", m, m, &dwork[i25], m, &dwork[iwrk], (ftnlen)1);
#line 769 "SB10ZD.f"
    dgetrf_(m, m, &dwork[i25], m, &iwork[1], &info2);
#line 770 "SB10ZD.f"
    if (info2 != 0) {
#line 771 "SB10ZD.f"
	*info = 7;
#line 772 "SB10ZD.f"
	return 0;
#line 773 "SB10ZD.f"
    }
#line 774 "SB10ZD.f"
    dgecon_("1", m, &dwork[i25], m, &anorm, &rcond[6], &dwork[iwrk], &iwork[*
	    m + 1], &info2, (ftnlen)1);

/*     Return if the matrix is singular to working precision. */

#line 779 "SB10ZD.f"
    if (rcond[6] < toll) {
#line 780 "SB10ZD.f"
	*info = 7;
#line 781 "SB10ZD.f"
	return 0;
#line 782 "SB10ZD.f"
    }
#line 783 "SB10ZD.f"
    dgetrs_("N", m, n, &dwork[i25], m, &iwork[1], &ck[ck_offset], ldck, &
	    info2, (ftnlen)1);
/*                  ^ ^ */
/*     Compute Dk = M*Dk . */

#line 787 "SB10ZD.f"
    dgetrs_("N", m, np, &dwork[i25], m, &iwork[1], &dk[dk_offset], lddk, &
	    info2, (ftnlen)1);
/*             ^ */
/*     Compute Bk*D . */

#line 791 "SB10ZD.f"
    dgemm_("N", "N", n, m, np, &c_b6, &dwork[i23], n, &d__[d_offset], ldd, &
	    c_b5, &dwork[i26], n, (ftnlen)1, (ftnlen)1);
/*                  ^    ^ */
/*     Compute Ak = Ak - Bk*D*Ck. */

#line 796 "SB10ZD.f"
    dgemm_("N", "N", n, n, m, &c_b42, &dwork[i26], n, &ck[ck_offset], ldck, &
	    c_b6, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);
/*                  ^          ^  -1 */
/*     Compute Bk = Bk*(Ip + D*Dk)   . */

#line 801 "SB10ZD.f"
    anorm = dlange_("1", np, np, &dwork[i24], np, &dwork[iwrk], (ftnlen)1);
#line 802 "SB10ZD.f"
    dlacpy_("Full", n, np, &dwork[i23], n, &bk[bk_offset], ldbk, (ftnlen)4);
#line 803 "SB10ZD.f"
    mb02vd_("N", n, np, &dwork[i24], np, &iwork[1], &bk[bk_offset], ldbk, &
	    info2, (ftnlen)1);
#line 805 "SB10ZD.f"
    if (info2 != 0) {
#line 806 "SB10ZD.f"
	*info = 6;
#line 807 "SB10ZD.f"
	return 0;
#line 808 "SB10ZD.f"
    }
#line 809 "SB10ZD.f"
    dgecon_("1", np, &dwork[i24], np, &anorm, &rcond[5], &dwork[iwrk], &iwork[
	    *np + 1], &info2, (ftnlen)1);

/*     Return if the matrix is singular to working precision. */

#line 814 "SB10ZD.f"
    if (rcond[5] < toll) {
#line 815 "SB10ZD.f"
	*info = 6;
#line 816 "SB10ZD.f"
	return 0;
#line 817 "SB10ZD.f"
    }

/*     Workspace usage. */

#line 821 "SB10ZD.f"
    i2 = *np * *np + 1;
#line 822 "SB10ZD.f"
    i3 = i2 + *n * *np;
#line 823 "SB10ZD.f"
    i4 = i3 + *m * *m;
#line 824 "SB10ZD.f"
    i5 = i4 + *n * *m;
#line 825 "SB10ZD.f"
    i6 = i5 + *np * *n;
#line 826 "SB10ZD.f"
    i7 = i6 + *m * *n;
#line 827 "SB10ZD.f"
    i8 = i7 + n2 * n2;
#line 828 "SB10ZD.f"
    i9 = i8 + n2;

#line 830 "SB10ZD.f"
    iwrk = i9 + n2;

/*     Compute Ip - D*Dk . */

#line 834 "SB10ZD.f"
    dlaset_("Full", np, np, &c_b5, &c_b6, &dwork[1], np, (ftnlen)4);
#line 835 "SB10ZD.f"
    dgemm_("N", "N", np, np, m, &c_b42, &d__[d_offset], ldd, &dk[dk_offset], 
	    lddk, &c_b6, &dwork[1], np, (ftnlen)1, (ftnlen)1);
/*                         -1 */
/*     Compute Bk*(Ip-D*Dk)   . */

#line 840 "SB10ZD.f"
    dlacpy_("Full", n, np, &bk[bk_offset], ldbk, &dwork[i2], n, (ftnlen)4);
#line 841 "SB10ZD.f"
    mb02vd_("N", n, np, &dwork[1], np, &iwork[1], &dwork[i2], n, &info2, (
	    ftnlen)1);
#line 842 "SB10ZD.f"
    if (info2 != 0) {
#line 843 "SB10ZD.f"
	*info = 8;
#line 844 "SB10ZD.f"
	return 0;
#line 845 "SB10ZD.f"
    }

/*     Compute Im - Dk*D . */

#line 849 "SB10ZD.f"
    dlaset_("Full", m, m, &c_b5, &c_b6, &dwork[i3], m, (ftnlen)4);
#line 850 "SB10ZD.f"
    dgemm_("N", "N", m, m, np, &c_b42, &dk[dk_offset], lddk, &d__[d_offset], 
	    ldd, &c_b6, &dwork[i3], m, (ftnlen)1, (ftnlen)1);
/*                        -1 */
/*     Compute B*(Im-Dk*D)    . */

#line 855 "SB10ZD.f"
    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[i4], n, (ftnlen)4);
#line 856 "SB10ZD.f"
    mb02vd_("N", n, m, &dwork[i3], m, &iwork[1], &dwork[i4], n, &info2, (
	    ftnlen)1);
#line 858 "SB10ZD.f"
    if (info2 != 0) {
#line 859 "SB10ZD.f"
	*info = 9;
#line 860 "SB10ZD.f"
	return 0;
#line 861 "SB10ZD.f"
    }

/*     Compute D*Ck . */

#line 865 "SB10ZD.f"
    dgemm_("N", "N", np, n, m, &c_b6, &d__[d_offset], ldd, &ck[ck_offset], 
	    ldck, &c_b5, &dwork[i5], np, (ftnlen)1, (ftnlen)1);

/*     Compute Dk*C . */

#line 870 "SB10ZD.f"
    dgemm_("N", "N", m, n, np, &c_b6, &dk[dk_offset], lddk, &c__[c_offset], 
	    ldc, &c_b5, &dwork[i6], m, (ftnlen)1, (ftnlen)1);

/*     Compute the closed-loop state matrix. */

#line 875 "SB10ZD.f"
    dlacpy_("F", n, n, &a[a_offset], lda, &dwork[i7], &n2, (ftnlen)1);
#line 876 "SB10ZD.f"
    dgemm_("N", "N", n, n, m, &c_b6, &dwork[i4], n, &dwork[i6], m, &c_b6, &
	    dwork[i7], &n2, (ftnlen)1, (ftnlen)1);
#line 878 "SB10ZD.f"
    dgemm_("N", "N", n, n, m, &c_b6, &dwork[i4], n, &ck[ck_offset], ldck, &
	    c_b5, &dwork[i7 + n2 * *n], &n2, (ftnlen)1, (ftnlen)1);
#line 880 "SB10ZD.f"
    dgemm_("N", "N", n, n, np, &c_b6, &dwork[i2], n, &c__[c_offset], ldc, &
	    c_b5, &dwork[i7 + *n], &n2, (ftnlen)1, (ftnlen)1);
#line 882 "SB10ZD.f"
    dlacpy_("F", n, n, &ak[ak_offset], ldak, &dwork[i7 + n2 * *n + *n], &n2, (
	    ftnlen)1);
#line 883 "SB10ZD.f"
    dgemm_("N", "N", n, n, np, &c_b6, &dwork[i2], n, &dwork[i5], np, &c_b6, &
	    dwork[i7 + n2 * *n + *n], &n2, (ftnlen)1, (ftnlen)1);

/*     Compute the closed-loop poles. */

#line 888 "SB10ZD.f"
    i__1 = *ldwork - iwrk + 1;
#line 888 "SB10ZD.f"
    dgees_("N", "N", (L_fp)select_, &n2, &dwork[i7], &n2, &sdim, &dwork[i8], &
	    dwork[i9], &dwork[iwrk], n, &dwork[iwrk], &i__1, &bwork[1], &
	    info2, (ftnlen)1, (ftnlen)1);
#line 891 "SB10ZD.f"
    if (info2 != 0) {
#line 892 "SB10ZD.f"
	*info = 3;
#line 893 "SB10ZD.f"
	return 0;
#line 894 "SB10ZD.f"
    }
/* Computing MAX */
#line 895 "SB10ZD.f"
    i__1 = lwamax, i__2 = (integer) dwork[iwrk] + iwrk - 1;
#line 895 "SB10ZD.f"
    lwamax = max(i__1,i__2);

/*     Check the stability of the closed-loop system. */

#line 899 "SB10ZD.f"
    ns = 0;

#line 901 "SB10ZD.f"
    i__1 = n2 - 1;
#line 901 "SB10ZD.f"
    for (i__ = 0; i__ <= i__1; ++i__) {
#line 902 "SB10ZD.f"
	if (dlapy2_(&dwork[i8 + i__], &dwork[i9 + i__]) > 1.) {
#line 902 "SB10ZD.f"
	    ++ns;
#line 902 "SB10ZD.f"
	}
#line 904 "SB10ZD.f"
/* L140: */
#line 904 "SB10ZD.f"
    }

#line 906 "SB10ZD.f"
    if (ns > 0) {
#line 907 "SB10ZD.f"
	*info = 10;
#line 908 "SB10ZD.f"
	return 0;
#line 909 "SB10ZD.f"
    }

#line 911 "SB10ZD.f"
    dwork[1] = (doublereal) lwamax;
#line 912 "SB10ZD.f"
    return 0;
/* *** Last line of SB10ZD *** */
} /* sb10zd_ */

