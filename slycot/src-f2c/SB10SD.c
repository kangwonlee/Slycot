#line 1 "SB10SD.f"
/* SB10SD.f -- translated by f2c (version 20100827).
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

#line 1 "SB10SD.f"
/* Table of constant values */

static doublereal c_b7 = -1.;
static doublereal c_b8 = 1.;
static doublereal c_b12 = 0.;
static integer c__1 = 1;

/* Subroutine */ int sb10sd_(integer *n, integer *m, integer *np, integer *
	ncon, integer *nmeas, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *d__, integer 
	*ldd, doublereal *ak, integer *ldak, doublereal *bk, integer *ldbk, 
	doublereal *ck, integer *ldck, doublereal *dk, integer *lddk, 
	doublereal *x, integer *ldx, doublereal *y, integer *ldy, doublereal *
	rcond, doublereal *tol, integer *iwork, doublereal *dwork, integer *
	ldwork, logical *bwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ak_dim1, ak_offset, b_dim1, b_offset, bk_dim1, 
	    bk_offset, c_dim1, c_offset, ck_dim1, ck_offset, d_dim1, d_offset,
	     dk_dim1, dk_offset, x_dim1, x_offset, y_dim1, y_offset, i__1, 
	    i__2, i__3, i__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, m1, m2, nd1, nd2, np1, np2, iw2, iwb, iwc, iwg, iwi, 
	    iwq, iwr, iws, iwt, iwu, iwv;
    static doublereal sepd, ferr, toll;
    static integer iwrk, info2;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     sb02od_(char *, char *, char *, char *, char *, char *, integer *
	    , integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, logical *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen, ftnlen, ftnlen), sb02sd_(char *, char *, 
	    char *, char *, char *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen), 
	    mb01rx_(char *, char *, char *, integer *, integer *, doublereal *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dsyrk_(
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal rcond2;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), dpocon_(char *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen);
    static integer lwamax;
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
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

/*     To compute the matrices of the H2 optimal controller */

/*              | AK | BK | */
/*          K = |----|----|, */
/*              | CK | DK | */

/*     for the normalized discrete-time system */

/*                   | A  | B1  B2  |   | A | B | */
/*               P = |----|---------| = |---|---| */
/*                   | C1 | D11 D12 |   | C | D | */
/*                   | C2 | D21  0  | */

/*     where B2 has as column size the number of control inputs (NCON) */
/*     and C2 has as row size the number of measurements (NMEAS) being */
/*     provided to the controller. */

/*     It is assumed that */

/*     (A1) (A,B2) is stabilizable and (C2,A) is detectable, */

/*     (A2) D12 is full column rank with D12 = | 0 | and D21 is */
/*                                             | I | */
/*          full row rank with D21 = | 0 I | as obtained by the */
/*          SLICOT Library routine SB10PD, */

/*               j*Theta */
/*     (A3) | A-e       *I  B2  | has full column rank for all */
/*          |    C1         D12 | */

/*          0 <= Theta < 2*Pi , */


/*               j*Theta */
/*     (A4) | A-e       *I  B1  | has full row rank for all */
/*          |    C2         D21 | */

/*          0 <= Theta < 2*Pi . */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The column size of the matrix B.  M >= 0. */

/*     NP      (input) INTEGER */
/*             The row size of the matrix C.  NP >= 0. */

/*     NCON    (input) INTEGER */
/*             The number of control inputs (M2).  M >= NCON >= 0, */
/*             NP-NMEAS >= NCON. */

/*     NMEAS   (input) INTEGER */
/*             The number of measurements (NP2).  NP >= NMEAS >= 0, */
/*             M-NCON >= NMEAS. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             system state matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             system input matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading NP-by-N part of this array must contain the */
/*             system output matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= max(1,NP). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading NP-by-M part of this array must contain the */
/*             system input/output matrix D. Only the leading */
/*             (NP-NP2)-by-(M-M2) submatrix D11 is used. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= max(1,NP). */

/*     AK      (output) DOUBLE PRECISION array, dimension (LDAK,N) */
/*             The leading N-by-N part of this array contains the */
/*             controller state matrix AK. */

/*     LDAK    INTEGER */
/*             The leading dimension of the array AK.  LDAK >= max(1,N). */

/*     BK      (output) DOUBLE PRECISION array, dimension (LDBK,NMEAS) */
/*             The leading N-by-NMEAS part of this array contains the */
/*             controller input matrix BK. */

/*     LDBK    INTEGER */
/*             The leading dimension of the array BK.  LDBK >= max(1,N). */

/*     CK      (output) DOUBLE PRECISION array, dimension (LDCK,N) */
/*             The leading NCON-by-N part of this array contains the */
/*             controller output matrix CK. */

/*     LDCK    INTEGER */
/*             The leading dimension of the array CK. */
/*             LDCK >= max(1,NCON). */

/*     DK      (output) DOUBLE PRECISION array, dimension (LDDK,NMEAS) */
/*             The leading NCON-by-NMEAS part of this array contains the */
/*             controller input/output matrix DK. */

/*     LDDK    INTEGER */
/*             The leading dimension of the array DK. */
/*             LDDK >= max(1,NCON). */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             The leading N-by-N part of this array contains the matrix */
/*             X, solution of the X-Riccati equation. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X.  LDX >= max(1,N). */

/*     Y       (output) DOUBLE PRECISION array, dimension (LDY,N) */
/*             The leading N-by-N part of this array contains the matrix */
/*             Y, solution of the Y-Riccati equation. */

/*     LDY     INTEGER */
/*             The leading dimension of the array Y.  LDY >= max(1,N). */

/*     RCOND   (output) DOUBLE PRECISION array, dimension (4) */
/*             RCOND contains estimates of the reciprocal condition */
/*             numbers of the matrices which are to be inverted and the */
/*             reciprocal condition numbers of the Riccati equations */
/*             which have to be solved during the computation of the */
/*             controller. (See the description of the algorithm in [2].) */
/*             RCOND(1) contains the reciprocal condition number of the */
/*                      matrix Im2 + B2'*X2*B2; */
/*             RCOND(2) contains the reciprocal condition number of the */
/*                      matrix Ip2 + C2*Y2*C2'; */
/*             RCOND(3) contains the reciprocal condition number of the */
/*                      X-Riccati equation; */
/*             RCOND(4) contains the reciprocal condition number of the */
/*                      Y-Riccati equation. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             Tolerance used in determining the nonsingularity of the */
/*             matrices which must be inverted. If TOL <= 0, then a */
/*             default value equal to sqrt(EPS) is used, where EPS is the */
/*             relative machine precision. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension max(M2,2*N,N*N,NP2) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal */
/*             LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= max(1, 14*N*N+6*N+max(14*N+23,16*N), */
/*                              M2*(N+M2+max(3,M1)), NP2*(N+NP2+3)), */
/*             where M1 = M - M2. */
/*             For good performance, LDWORK must generally be larger. */

/*     BWORK   LOGICAL array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the X-Riccati equation was not solved */
/*                   successfully; */
/*             = 2:  if the matrix Im2 + B2'*X2*B2 is not positive */
/*                   definite, or it is numerically singular (with */
/*                   respect to the tolerance TOL); */
/*             = 3:  if the Y-Riccati equation was not solved */
/*                   successfully; */
/*             = 4:  if the matrix Ip2 + C2*Y2*C2' is not positive */
/*                   definite, or it is numerically singular (with */
/*                   respect to the tolerance TOL). */

/*     METHOD */

/*     The routine implements the formulas given in [1]. The X- and */
/*     Y-Riccati equations are solved with condition estimates. */

/*     REFERENCES */

/*     [1] Zhou, K., Doyle, J.C., and Glover, K. */
/*         Robust and Optimal Control. */
/*         Prentice-Hall, Upper Saddle River, NJ, 1996. */

/*     [2] Petkov, P.Hr., Gu, D.W., and Konstantinov, M.M. */
/*         Fortran 77 routines for Hinf and H2 design of linear */
/*         discrete-time control systems. */
/*         Report 99-8, Department of Engineering, Leicester University, */
/*         April 1999. */

/*     NUMERICAL ASPECTS */

/*     The accuracy of the result depends on the condition numbers of the */
/*     matrices which are to be inverted and on the condition numbers of */
/*     the matrix Riccati equations which are to be solved in the */
/*     computation of the controller. (The corresponding reciprocal */
/*     condition numbers are given in the output array RCOND.) */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, April 1999. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 1999, */
/*     January 2003. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, H2 optimal control, LQG, LQR, optimal */
/*     regulator, robust control. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test input parameters. */

#line 296 "SB10SD.f"
    /* Parameter adjustments */
#line 296 "SB10SD.f"
    a_dim1 = *lda;
#line 296 "SB10SD.f"
    a_offset = 1 + a_dim1;
#line 296 "SB10SD.f"
    a -= a_offset;
#line 296 "SB10SD.f"
    b_dim1 = *ldb;
#line 296 "SB10SD.f"
    b_offset = 1 + b_dim1;
#line 296 "SB10SD.f"
    b -= b_offset;
#line 296 "SB10SD.f"
    c_dim1 = *ldc;
#line 296 "SB10SD.f"
    c_offset = 1 + c_dim1;
#line 296 "SB10SD.f"
    c__ -= c_offset;
#line 296 "SB10SD.f"
    d_dim1 = *ldd;
#line 296 "SB10SD.f"
    d_offset = 1 + d_dim1;
#line 296 "SB10SD.f"
    d__ -= d_offset;
#line 296 "SB10SD.f"
    ak_dim1 = *ldak;
#line 296 "SB10SD.f"
    ak_offset = 1 + ak_dim1;
#line 296 "SB10SD.f"
    ak -= ak_offset;
#line 296 "SB10SD.f"
    bk_dim1 = *ldbk;
#line 296 "SB10SD.f"
    bk_offset = 1 + bk_dim1;
#line 296 "SB10SD.f"
    bk -= bk_offset;
#line 296 "SB10SD.f"
    ck_dim1 = *ldck;
#line 296 "SB10SD.f"
    ck_offset = 1 + ck_dim1;
#line 296 "SB10SD.f"
    ck -= ck_offset;
#line 296 "SB10SD.f"
    dk_dim1 = *lddk;
#line 296 "SB10SD.f"
    dk_offset = 1 + dk_dim1;
#line 296 "SB10SD.f"
    dk -= dk_offset;
#line 296 "SB10SD.f"
    x_dim1 = *ldx;
#line 296 "SB10SD.f"
    x_offset = 1 + x_dim1;
#line 296 "SB10SD.f"
    x -= x_offset;
#line 296 "SB10SD.f"
    y_dim1 = *ldy;
#line 296 "SB10SD.f"
    y_offset = 1 + y_dim1;
#line 296 "SB10SD.f"
    y -= y_offset;
#line 296 "SB10SD.f"
    --rcond;
#line 296 "SB10SD.f"
    --iwork;
#line 296 "SB10SD.f"
    --dwork;
#line 296 "SB10SD.f"
    --bwork;
#line 296 "SB10SD.f"

#line 296 "SB10SD.f"
    /* Function Body */
#line 296 "SB10SD.f"
    m1 = *m - *ncon;
#line 297 "SB10SD.f"
    m2 = *ncon;
#line 298 "SB10SD.f"
    np1 = *np - *nmeas;
#line 299 "SB10SD.f"
    np2 = *nmeas;

#line 301 "SB10SD.f"
    *info = 0;
#line 302 "SB10SD.f"
    if (*n < 0) {
#line 303 "SB10SD.f"
	*info = -1;
#line 304 "SB10SD.f"
    } else if (*m < 0) {
#line 305 "SB10SD.f"
	*info = -2;
#line 306 "SB10SD.f"
    } else if (*np < 0) {
#line 307 "SB10SD.f"
	*info = -3;
#line 308 "SB10SD.f"
    } else if (*ncon < 0 || m1 < 0 || m2 > np1) {
#line 309 "SB10SD.f"
	*info = -4;
#line 310 "SB10SD.f"
    } else if (*nmeas < 0 || np1 < 0 || np2 > m1) {
#line 311 "SB10SD.f"
	*info = -5;
#line 312 "SB10SD.f"
    } else if (*lda < max(1,*n)) {
#line 313 "SB10SD.f"
	*info = -7;
#line 314 "SB10SD.f"
    } else if (*ldb < max(1,*n)) {
#line 315 "SB10SD.f"
	*info = -9;
#line 316 "SB10SD.f"
    } else if (*ldc < max(1,*np)) {
#line 317 "SB10SD.f"
	*info = -11;
#line 318 "SB10SD.f"
    } else if (*ldd < max(1,*np)) {
#line 319 "SB10SD.f"
	*info = -13;
#line 320 "SB10SD.f"
    } else if (*ldak < max(1,*n)) {
#line 321 "SB10SD.f"
	*info = -15;
#line 322 "SB10SD.f"
    } else if (*ldbk < max(1,*n)) {
#line 323 "SB10SD.f"
	*info = -17;
#line 324 "SB10SD.f"
    } else if (*ldck < max(1,m2)) {
#line 325 "SB10SD.f"
	*info = -19;
#line 326 "SB10SD.f"
    } else if (*lddk < max(1,m2)) {
#line 327 "SB10SD.f"
	*info = -21;
#line 328 "SB10SD.f"
    } else if (*ldx < max(1,*n)) {
#line 329 "SB10SD.f"
	*info = -23;
#line 330 "SB10SD.f"
    } else if (*ldy < max(1,*n)) {
#line 331 "SB10SD.f"
	*info = -25;
#line 332 "SB10SD.f"
    } else {

/*        Compute workspace. */

/* Computing MAX */
/* Computing MAX */
#line 336 "SB10SD.f"
	i__3 = *n * 14 + 23, i__4 = *n << 4;
#line 336 "SB10SD.f"
	i__1 = 1, i__2 = *n * 14 * *n + *n * 6 + max(i__3,i__4), i__1 = max(
		i__1,i__2), i__2 = m2 * (*n + m2 + max(3,m1)), i__1 = max(
		i__1,i__2), i__2 = np2 * (*n + np2 + 3);
#line 336 "SB10SD.f"
	minwrk = max(i__1,i__2);
#line 338 "SB10SD.f"
	if (*ldwork < minwrk) {
#line 338 "SB10SD.f"
	    *info = -30;
#line 338 "SB10SD.f"
	}
#line 340 "SB10SD.f"
    }
#line 341 "SB10SD.f"
    if (*info != 0) {
#line 342 "SB10SD.f"
	i__1 = -(*info);
#line 342 "SB10SD.f"
	xerbla_("SB10SD", &i__1, (ftnlen)6);
#line 343 "SB10SD.f"
	return 0;
#line 344 "SB10SD.f"
    }

/*     Quick return if possible. */

#line 348 "SB10SD.f"
    if (*n == 0 || *m == 0 || *np == 0 || m1 == 0 || m2 == 0 || np1 == 0 || 
	    np2 == 0) {
#line 350 "SB10SD.f"
	rcond[1] = 1.;
#line 351 "SB10SD.f"
	rcond[2] = 1.;
#line 352 "SB10SD.f"
	rcond[3] = 1.;
#line 353 "SB10SD.f"
	rcond[4] = 1.;
#line 354 "SB10SD.f"
	dwork[1] = 1.;
#line 355 "SB10SD.f"
	return 0;
#line 356 "SB10SD.f"
    }

#line 358 "SB10SD.f"
    nd1 = np1 - m2;
#line 359 "SB10SD.f"
    nd2 = m1 - np2;
#line 360 "SB10SD.f"
    toll = *tol;
#line 361 "SB10SD.f"
    if (toll <= 0.) {

/*        Set the default value of the tolerance for nonsingularity test. */

#line 365 "SB10SD.f"
	toll = sqrt(dlamch_("Epsilon", (ftnlen)7));
#line 366 "SB10SD.f"
    }

/*     Workspace usage. */

#line 370 "SB10SD.f"
    iwq = 1;
#line 371 "SB10SD.f"
    iwg = iwq + *n * *n;
#line 372 "SB10SD.f"
    iwr = iwg + *n * *n;
#line 373 "SB10SD.f"
    iwi = iwr + (*n << 1);
#line 374 "SB10SD.f"
    iwb = iwi + (*n << 1);
#line 375 "SB10SD.f"
    iws = iwb + (*n << 1);
#line 376 "SB10SD.f"
    iwt = iws + (*n << 2) * *n;
#line 377 "SB10SD.f"
    iwu = iwt + (*n << 2) * *n;
#line 378 "SB10SD.f"
    iwrk = iwu + (*n << 2) * *n;
#line 379 "SB10SD.f"
    iwc = iwr;
#line 380 "SB10SD.f"
    iwv = iwc + *n * *n;

/*     Compute Ax = A - B2*D12'*C1 in AK . */

#line 384 "SB10SD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &ak[ak_offset], ldak, (ftnlen)4);
#line 385 "SB10SD.f"
    dgemm_("N", "N", n, n, &m2, &c_b7, &b[(m1 + 1) * b_dim1 + 1], ldb, &c__[
	    nd1 + 1 + c_dim1], ldc, &c_b8, &ak[ak_offset], ldak, (ftnlen)1, (
	    ftnlen)1);

/*     Compute Cx = C1'*C1 - C1'*D12*D12'*C1 . */

#line 390 "SB10SD.f"
    if (nd1 > 0) {
#line 391 "SB10SD.f"
	dsyrk_("L", "T", n, &nd1, &c_b8, &c__[c_offset], ldc, &c_b12, &dwork[
		iwq], n, (ftnlen)1, (ftnlen)1);
#line 393 "SB10SD.f"
    } else {
#line 394 "SB10SD.f"
	dlaset_("L", n, n, &c_b12, &c_b12, &dwork[iwq], n, (ftnlen)1);
#line 395 "SB10SD.f"
    }

/*     Compute Dx = B2*B2' . */

#line 399 "SB10SD.f"
    dsyrk_("L", "N", n, &m2, &c_b8, &b[(m1 + 1) * b_dim1 + 1], ldb, &c_b12, &
	    dwork[iwg], n, (ftnlen)1, (ftnlen)1);

/*     Solution of the discrete-time Riccati equation */
/*        Ax'*inv(In + X2*Dx)*X2*Ax - X2 + Cx  = 0 . */
/*     Workspace:  need   14*N*N + 6*N + max(14*N+23,16*N); */
/*                 prefer larger. */

#line 407 "SB10SD.f"
    i__1 = *n << 1;
#line 407 "SB10SD.f"
    i__2 = *n << 1;
#line 407 "SB10SD.f"
    i__3 = *n << 1;
#line 407 "SB10SD.f"
    i__4 = *ldwork - iwrk + 1;
#line 407 "SB10SD.f"
    sb02od_("D", "G", "N", "L", "Z", "S", n, &m2, &np1, &ak[ak_offset], ldak, 
	    &dwork[iwg], n, &dwork[iwq], n, &dwork[iwrk], m, &dwork[iwrk], n, 
	    &rcond2, &x[x_offset], ldx, &dwork[iwr], &dwork[iwi], &dwork[iwb],
	     &dwork[iws], &i__1, &dwork[iwt], &i__2, &dwork[iwu], &i__3, &
	    toll, &iwork[1], &dwork[iwrk], &i__4, &bwork[1], &info2, (ftnlen)
	    1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 413 "SB10SD.f"
    if (info2 > 0) {
#line 414 "SB10SD.f"
	*info = 1;
#line 415 "SB10SD.f"
	return 0;
#line 416 "SB10SD.f"
    }
#line 417 "SB10SD.f"
    lwamax = (integer) dwork[iwrk] + iwrk - 1;

/*     Condition estimation. */
/*     Workspace:  need   4*N*N + max(N*N+5*N,max(3,2*N*N)+N*N); */
/*                 prefer larger. */

#line 423 "SB10SD.f"
    iwrk = iwv + *n * *n;
#line 424 "SB10SD.f"
    i__1 = *ldwork - iwrk + 1;
#line 424 "SB10SD.f"
    sb02sd_("C", "N", "N", "L", "O", n, &ak[ak_offset], ldak, &dwork[iwc], n, 
	    &dwork[iwv], n, &dwork[iwg], n, &dwork[iwq], n, &x[x_offset], ldx,
	     &sepd, &rcond[3], &ferr, &iwork[1], &dwork[iwrk], &i__1, &info2, 
	    (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 428 "SB10SD.f"
    if (info2 > 0) {
#line 428 "SB10SD.f"
	rcond[3] = 0.;
#line 428 "SB10SD.f"
    }
/* Computing MAX */
#line 429 "SB10SD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 429 "SB10SD.f"
    lwamax = max(i__1,lwamax);

/*     Workspace usage. */

#line 433 "SB10SD.f"
    iw2 = m2 * *n + 1;
#line 434 "SB10SD.f"
    iwrk = iw2 + m2 * m2;

/*     Compute B2'*X2 . */

#line 438 "SB10SD.f"
    dgemm_("T", "N", &m2, n, n, &c_b8, &b[(m1 + 1) * b_dim1 + 1], ldb, &x[
	    x_offset], ldx, &c_b12, &dwork[1], &m2, (ftnlen)1, (ftnlen)1);

/*     Compute Im2 + B2'*X2*B2 . */

#line 443 "SB10SD.f"
    dlaset_("L", &m2, &m2, &c_b12, &c_b8, &dwork[iw2], &m2, (ftnlen)1);
#line 444 "SB10SD.f"
    mb01rx_("Left", "Lower", "N", &m2, n, &c_b8, &c_b8, &dwork[iw2], &m2, &
	    dwork[1], &m2, &b[(m1 + 1) * b_dim1 + 1], ldb, &info2, (ftnlen)4, 
	    (ftnlen)5, (ftnlen)1);

/*     Compute the Cholesky factorization of Im2 + B2'*X2*B2 . */
/*     Workspace:  need   M2*N + M2*M2 + max(3*M2,M2*M1); */
/*                 prefer larger. */

#line 451 "SB10SD.f"
    anorm = dlansy_("I", "L", &m2, &dwork[iw2], &m2, &dwork[iwrk], (ftnlen)1, 
	    (ftnlen)1);
#line 452 "SB10SD.f"
    dpotrf_("L", &m2, &dwork[iw2], &m2, &info2, (ftnlen)1);
#line 453 "SB10SD.f"
    if (info2 > 0) {
#line 454 "SB10SD.f"
	*info = 2;
#line 455 "SB10SD.f"
	return 0;
#line 456 "SB10SD.f"
    }
#line 457 "SB10SD.f"
    dpocon_("L", &m2, &dwork[iw2], &m2, &anorm, &rcond[1], &dwork[iwrk], &
	    iwork[1], &info2, (ftnlen)1);

/*     Return if the matrix is singular to working precision. */

#line 462 "SB10SD.f"
    if (rcond[1] < toll) {
#line 463 "SB10SD.f"
	*info = 2;
#line 464 "SB10SD.f"
	return 0;
#line 465 "SB10SD.f"
    }

/*     Compute -( B2'*X2*A + D12'*C1 ) in CK . */

#line 469 "SB10SD.f"
    dlacpy_("Full", &m2, n, &c__[nd1 + 1 + c_dim1], ldc, &ck[ck_offset], ldck,
	     (ftnlen)4);
#line 470 "SB10SD.f"
    dgemm_("N", "N", &m2, n, n, &c_b7, &dwork[1], &m2, &a[a_offset], lda, &
	    c_b7, &ck[ck_offset], ldck, (ftnlen)1, (ftnlen)1);

/*     Compute F2 = -inv( Im2 + B2'*X2*B2 )*( B2'*X2*A + D12'*C1 ) . */

#line 475 "SB10SD.f"
    dpotrs_("L", &m2, n, &dwork[iw2], &m2, &ck[ck_offset], ldck, &info2, (
	    ftnlen)1);

/*     Compute -( B2'*X2*B1 + D12'*D11 ) . */

#line 479 "SB10SD.f"
    dlacpy_("Full", &m2, &m1, &d__[nd1 + 1 + d_dim1], ldd, &dwork[iwrk], &m2, 
	    (ftnlen)4);
#line 481 "SB10SD.f"
    dgemm_("N", "N", &m2, &m1, n, &c_b7, &dwork[1], &m2, &b[b_offset], ldb, &
	    c_b7, &dwork[iwrk], &m2, (ftnlen)1, (ftnlen)1);

/*     Compute F0 = -inv( Im2 + B2'*X2*B2 )*( B2'*X2*B1 + D12'*D11 ) . */

#line 486 "SB10SD.f"
    dpotrs_("L", &m2, &m1, &dwork[iw2], &m2, &dwork[iwrk], &m2, &info2, (
	    ftnlen)1);

/*     Save F0*D21' in DK . */

#line 491 "SB10SD.f"
    dlacpy_("Full", &m2, &np2, &dwork[iwrk + nd2 * m2], &m2, &dk[dk_offset], 
	    lddk, (ftnlen)4);

/*     Workspace usage. */

#line 496 "SB10SD.f"
    iwrk = iwu + (*n << 2) * *n;

/*     Compute Ay = A - B1*D21'*C2 in AK . */

#line 500 "SB10SD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &ak[ak_offset], ldak, (ftnlen)4);
#line 501 "SB10SD.f"
    dgemm_("N", "N", n, n, &np2, &c_b7, &b[(nd2 + 1) * b_dim1 + 1], ldb, &c__[
	    np1 + 1 + c_dim1], ldc, &c_b8, &ak[ak_offset], ldak, (ftnlen)1, (
	    ftnlen)1);

/*     Transpose Ay in-situ. */

#line 506 "SB10SD.f"
    i__1 = *n - 1;
#line 506 "SB10SD.f"
    for (j = 1; j <= i__1; ++j) {
#line 507 "SB10SD.f"
	dswap_(&j, &ak[j + 1 + ak_dim1], ldak, &ak[(j + 1) * ak_dim1 + 1], &
		c__1);
#line 508 "SB10SD.f"
/* L20: */
#line 508 "SB10SD.f"
    }

/*     Compute Cy = B1*B1' - B1*D21'*D21*B1' . */

#line 512 "SB10SD.f"
    if (nd2 > 0) {
#line 513 "SB10SD.f"
	dsyrk_("U", "N", n, &nd2, &c_b8, &b[b_offset], ldb, &c_b12, &dwork[
		iwq], n, (ftnlen)1, (ftnlen)1);
#line 515 "SB10SD.f"
    } else {
#line 516 "SB10SD.f"
	dlaset_("U", n, n, &c_b12, &c_b12, &dwork[iwq], n, (ftnlen)1);
#line 517 "SB10SD.f"
    }

/*     Compute Dy = C2'*C2 . */

#line 521 "SB10SD.f"
    dsyrk_("U", "T", n, &np2, &c_b8, &c__[np1 + 1 + c_dim1], ldc, &c_b12, &
	    dwork[iwg], n, (ftnlen)1, (ftnlen)1);

/*     Solution of the discrete-time Riccati equation */
/*        Ay*inv( In + Y2*Dy )*Y2*Ay' - Y2 + Cy = 0 . */

#line 527 "SB10SD.f"
    i__1 = *n << 1;
#line 527 "SB10SD.f"
    i__2 = *n << 1;
#line 527 "SB10SD.f"
    i__3 = *n << 1;
#line 527 "SB10SD.f"
    i__4 = *ldwork - iwrk + 1;
#line 527 "SB10SD.f"
    sb02od_("D", "G", "N", "U", "Z", "S", n, &np2, &m1, &ak[ak_offset], ldak, 
	    &dwork[iwg], n, &dwork[iwq], n, &dwork[iwrk], m, &dwork[iwrk], n, 
	    &rcond2, &y[y_offset], ldy, &dwork[iwr], &dwork[iwi], &dwork[iwb],
	     &dwork[iws], &i__1, &dwork[iwt], &i__2, &dwork[iwu], &i__3, &
	    toll, &iwork[1], &dwork[iwrk], &i__4, &bwork[1], &info2, (ftnlen)
	    1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 533 "SB10SD.f"
    if (info2 > 0) {
#line 534 "SB10SD.f"
	*info = 3;
#line 535 "SB10SD.f"
	return 0;
#line 536 "SB10SD.f"
    }
/* Computing MAX */
#line 537 "SB10SD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 537 "SB10SD.f"
    lwamax = max(i__1,lwamax);

/*     Condition estimation. */

#line 541 "SB10SD.f"
    iwrk = iwv + *n * *n;
#line 542 "SB10SD.f"
    i__1 = *ldwork - iwrk + 1;
#line 542 "SB10SD.f"
    sb02sd_("C", "N", "N", "U", "O", n, &ak[ak_offset], ldak, &dwork[iwc], n, 
	    &dwork[iwv], n, &dwork[iwg], n, &dwork[iwq], n, &y[y_offset], ldy,
	     &sepd, &rcond[4], &ferr, &iwork[1], &dwork[iwrk], &i__1, &info2, 
	    (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 546 "SB10SD.f"
    if (info2 > 0) {
#line 546 "SB10SD.f"
	rcond[4] = 0.;
#line 546 "SB10SD.f"
    }
/* Computing MAX */
#line 547 "SB10SD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 547 "SB10SD.f"
    lwamax = max(i__1,lwamax);

/*     Workspace usage. */

#line 551 "SB10SD.f"
    iw2 = *n * np2 + 1;
#line 552 "SB10SD.f"
    iwrk = iw2 + np2 * np2;

/*     Compute Y2*C2' . */

#line 556 "SB10SD.f"
    dgemm_("N", "T", n, &np2, n, &c_b8, &y[y_offset], ldy, &c__[np1 + 1 + 
	    c_dim1], ldc, &c_b12, &dwork[1], n, (ftnlen)1, (ftnlen)1);

/*     Compute Ip2 + C2*Y2*C2' . */

#line 561 "SB10SD.f"
    dlaset_("U", &np2, &np2, &c_b12, &c_b8, &dwork[iw2], &np2, (ftnlen)1);
#line 562 "SB10SD.f"
    mb01rx_("Left", "Upper", "N", &np2, n, &c_b8, &c_b8, &dwork[iw2], &np2, &
	    c__[np1 + 1 + c_dim1], ldc, &dwork[1], n, &info2, (ftnlen)4, (
	    ftnlen)5, (ftnlen)1);

/*     Compute the Cholesky factorization of Ip2 + C2*Y2*C2' . */

#line 567 "SB10SD.f"
    anorm = dlansy_("I", "U", &np2, &dwork[iw2], &np2, &dwork[iwrk], (ftnlen)
	    1, (ftnlen)1);
#line 568 "SB10SD.f"
    dpotrf_("U", &np2, &dwork[iw2], &np2, &info2, (ftnlen)1);
#line 569 "SB10SD.f"
    if (info2 > 0) {
#line 570 "SB10SD.f"
	*info = 4;
#line 571 "SB10SD.f"
	return 0;
#line 572 "SB10SD.f"
    }
#line 573 "SB10SD.f"
    dpocon_("U", &np2, &dwork[iw2], &np2, &anorm, &rcond[2], &dwork[iwrk], &
	    iwork[1], &info2, (ftnlen)1);

/*     Return if the matrix is singular to working precision. */

#line 578 "SB10SD.f"
    if (rcond[2] < toll) {
#line 579 "SB10SD.f"
	*info = 4;
#line 580 "SB10SD.f"
	return 0;
#line 581 "SB10SD.f"
    }

/*     Compute A*Y2*C2' + B1*D21' in BK . */

#line 585 "SB10SD.f"
    dlacpy_("Full", n, &np2, &b[(nd2 + 1) * b_dim1 + 1], ldb, &bk[bk_offset], 
	    ldbk, (ftnlen)4);
#line 586 "SB10SD.f"
    dgemm_("N", "N", n, &np2, n, &c_b8, &a[a_offset], lda, &dwork[1], n, &
	    c_b8, &bk[bk_offset], ldbk, (ftnlen)1, (ftnlen)1);

/*     Compute L2 = -( A*Y2*C2' + B1*D21' )*inv( Ip2 + C2*Y2*C2' ) . */

#line 591 "SB10SD.f"
    dtrsm_("R", "U", "N", "N", n, &np2, &c_b7, &dwork[iw2], &np2, &bk[
	    bk_offset], ldbk, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 593 "SB10SD.f"
    dtrsm_("R", "U", "T", "N", n, &np2, &c_b8, &dwork[iw2], &np2, &bk[
	    bk_offset], ldbk, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Compute F2*Y2*C2' + F0*D21' . */

#line 598 "SB10SD.f"
    dgemm_("N", "N", &m2, &np2, n, &c_b8, &ck[ck_offset], ldck, &dwork[1], n, 
	    &c_b8, &dk[dk_offset], lddk, (ftnlen)1, (ftnlen)1);

/*     Compute DK = L0 = ( F2*Y2*C2' + F0*D21' )*inv( Ip2 + C2*Y2*C2' ) . */

#line 603 "SB10SD.f"
    dtrsm_("R", "U", "N", "N", &m2, &np2, &c_b8, &dwork[iw2], &np2, &dk[
	    dk_offset], lddk, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 605 "SB10SD.f"
    dtrsm_("R", "U", "T", "N", &m2, &np2, &c_b8, &dwork[iw2], &np2, &dk[
	    dk_offset], lddk, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Compute CK = F2 - L0*C2 . */

#line 610 "SB10SD.f"
    dgemm_("N", "N", &m2, n, &np2, &c_b7, &dk[dk_offset], lddk, &c__[np1 + 1 
	    + c_dim1], ldc, &c_b8, &ck[ck_offset], ldck, (ftnlen)1, (ftnlen)1)
	    ;

/*     Find AK = A + B2*( F2 - L0*C2 ) + L2*C2 . */

#line 615 "SB10SD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &ak[ak_offset], ldak, (ftnlen)4);
#line 616 "SB10SD.f"
    dgemm_("N", "N", n, n, &m2, &c_b8, &b[(m1 + 1) * b_dim1 + 1], ldb, &ck[
	    ck_offset], ldck, &c_b8, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)
	    1);
#line 618 "SB10SD.f"
    dgemm_("N", "N", n, n, &np2, &c_b8, &bk[bk_offset], ldbk, &c__[np1 + 1 + 
	    c_dim1], ldc, &c_b8, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);

/*     Find BK = -L2 + B2*L0 . */

#line 623 "SB10SD.f"
    dgemm_("N", "N", n, &np2, &m2, &c_b8, &b[(m1 + 1) * b_dim1 + 1], ldb, &dk[
	    dk_offset], lddk, &c_b7, &bk[bk_offset], ldbk, (ftnlen)1, (ftnlen)
	    1);

#line 626 "SB10SD.f"
    dwork[1] = (doublereal) lwamax;
#line 627 "SB10SD.f"
    return 0;
/* *** Last line of SB10SD *** */
} /* sb10sd_ */

