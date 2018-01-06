#line 1 "SB10DD.f"
/* SB10DD.f -- translated by f2c (version 20100827).
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

#line 1 "SB10DD.f"
/* Table of constant values */

static doublereal c_b24 = 1.;
static doublereal c_b25 = 0.;
static doublereal c_b74 = -1.;
static integer c__1 = 1;

/* Subroutine */ int sb10dd_(integer *n, integer *m, integer *np, integer *
	ncon, integer *nmeas, doublereal *gamma, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *ak, integer *ldak, 
	doublereal *bk, integer *ldbk, doublereal *ck, integer *ldck, 
	doublereal *dk, integer *lddk, doublereal *x, integer *ldx, 
	doublereal *z__, integer *ldz, doublereal *rcond, doublereal *tol, 
	integer *iwork, doublereal *dwork, integer *ldwork, logical *bwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ak_dim1, ak_offset, b_dim1, b_offset, bk_dim1, 
	    bk_offset, c_dim1, c_offset, ck_dim1, ck_offset, d_dim1, d_offset,
	     dk_dim1, dk_offset, x_dim1, x_offset, z_dim1, z_offset, i__1, 
	    i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    static integer j, m1, m2, ir2, ir3, is2, is3, np1, np2, iwb, iwc, iwd, 
	    iwg, iwh, iwi, iwl, iwq, iwr, iws, iwt, iwu, iwv, iww;
    static doublereal sepd, ferr, toll;
    static integer iwrk, info2;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), dscal_(
	    integer *, doublereal *, doublereal *, integer *), dgemm_(char *, 
	    char *, integer *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen), sb02od_(char *, char *, char *, char 
	    *, char *, char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    logical *, integer *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, 
	    ftnlen), sb02sd_(char *, char *, char *, char *, char *, integer *
	    , doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen, ftnlen), mb01ru_(char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen), mb01rx_(char *, char *, 
	    char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dsyrk_(
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal rcond2;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgecon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dgetrf_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), dgesvd_(char *, char *, integer 
	    *, integer *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), dpocon_(char *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen), dgetrs_(char *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen);
    static integer lwamax;
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtrcon_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), dpotrf_(char *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);
    static integer minwrk;
    extern /* Subroutine */ int dsytrf_(char *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     dsytrs_(char *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen);


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

/*     To compute the matrices of an H-infinity (sub)optimal n-state */
/*     controller */

/*                           | AK | BK | */
/*                       K = |----|----|, */
/*                           | CK | DK | */

/*     for the discrete-time system */

/*                   | A  | B1  B2  |   | A | B | */
/*               P = |----|---------| = |---|---| */
/*                   | C1 | D11 D12 |   | C | D | */
/*                   | C2 | D21 D22 | */

/*     and for a given value of gamma, where B2 has as column size the */
/*     number of control inputs (NCON) and C2 has as row size the number */
/*     of measurements (NMEAS) being provided to the controller. */

/*     It is assumed that */

/*     (A1) (A,B2) is stabilizable and (C2,A) is detectable, */

/*     (A2) D12 is full column rank and D21 is full row rank, */

/*               j*Theta */
/*     (A3) | A-e       *I  B2  | has full column rank for all */
/*          |    C1         D12 | */

/*          0 <= Theta < 2*Pi , */

/*               j*Theta */
/*     (A4) | A-e       *I  B1  |  has full row rank for all */
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

/*     GAMMA   (input) DOUBLE PRECISION */
/*             The value of gamma. It is assumed that gamma is */
/*             sufficiently large so that the controller is admissible. */
/*             GAMMA > 0. */

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
/*             system input/output matrix D. */

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

/*     Z       (output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             The leading N-by-N part of this array contains the matrix */
/*             Z, solution of the Z-Riccati equation. */

/*     LDZ     INTEGER */
/*             The leading dimension of the array Z.  LDZ >= max(1,N). */

/*     RCOND   (output) DOUBLE PRECISION array, dimension (8) */
/*             RCOND contains estimates of the reciprocal condition */
/*             numbers of the matrices which are to be inverted and */
/*             estimates of the reciprocal condition numbers of the */
/*             Riccati equations which have to be solved during the */
/*             computation of the controller. (See the description of */
/*             the algorithm in [2].) */
/*             RCOND(1) contains the reciprocal condition number of the */
/*                      matrix R3; */
/*             RCOND(2) contains the reciprocal condition number of the */
/*                      matrix R1 - R2'*inv(R3)*R2; */
/*             RCOND(3) contains the reciprocal condition number of the */
/*                      matrix V21; */
/*             RCOND(4) contains the reciprocal condition number of the */
/*                      matrix St3; */
/*             RCOND(5) contains the reciprocal condition number of the */
/*                      matrix V12; */
/*             RCOND(6) contains the reciprocal condition number of the */
/*                      matrix Im2 + DKHAT*D22 */
/*             RCOND(7) contains the reciprocal condition number of the */
/*                      X-Riccati equation; */
/*             RCOND(8) contains the reciprocal condition number of the */
/*                      Z-Riccati equation. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             Tolerance used in neglecting the small singular values */
/*             in rank determination. If TOL <= 0, then a default value */
/*             equal to 1000*EPS is used, where EPS is the relative */
/*             machine precision. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension max(2*max(M2,N),M,M2+NP2,N*N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal */
/*             LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= max(LW1,LW2,LW3,LW4), where */
/*             LW1 = (N+NP1+1)*(N+M2) + max(3*(N+M2)+N+NP1,5*(N+M2)); */
/*             LW2 = (N+NP2)*(N+M1+1) + max(3*(N+NP2)+N+M1,5*(N+NP2)); */
/*             LW3 = 13*N*N + 2*M*M + N*(8*M+NP2) + M1*(M2+NP2) + 6*N + */
/*                   max(14*N+23,16*N,2*N+M,3*M); */
/*             LW4 = 13*N*N + M*M + (8*N+M+M2+2*NP2)*(M2+NP2) + 6*N + */
/*                   N*(M+NP2) + max(14*N+23,16*N,2*N+M2+NP2,3*(M2+NP2)); */
/*             For good performance, LDWORK must generally be larger. */
/*             Denoting Q = max(M1,M2,NP1,NP2), an upper bound is */
/*             max((N+Q)*(N+Q+6),13*N*N + M*M + 2*Q*Q + N*(M+Q) + */
/*                 max(M*(M+7*N),2*Q*(8*N+M+2*Q)) + 6*N + */
/*                 max(14*N+23,16*N,2*N+max(M,2*Q),3*max(M,2*Q)). */

/*     BWORK   LOGICAL array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*                                      j*Theta */
/*             = 1:  if the matrix | A-e       *I  B2  | had not full */
/*                                 |      C1       D12 | */
/*                   column rank; */
/*                                      j*Theta */
/*             = 2:  if the matrix | A-e       *I  B1  | had not full */
/*                                 |      C2       D21 | */
/*                   row rank; */
/*             = 3:  if the matrix D12 had not full column rank; */
/*             = 4:  if the matrix D21 had not full row rank; */
/*             = 5:  if the controller is not admissible (too small value */
/*                   of gamma); */
/*             = 6:  if the X-Riccati equation was not solved */
/*                   successfully (the controller is not admissible or */
/*                   there are numerical difficulties); */
/*             = 7:  if the Z-Riccati equation was not solved */
/*                   successfully (the controller is not admissible or */
/*                   there are numerical difficulties); */
/*             = 8:  if the matrix Im2 + DKHAT*D22 is singular. */
/*             = 9:  if the singular value decomposition (SVD) algorithm */
/*                   did not converge (when computing the SVD of one of */
/*                   the matrices |A   B2 |, |A   B1 |, D12 or D21). */
/*                                |C1  D12|  |C2  D21| */

/*     METHOD */

/*     The routine implements the method presented in [1]. */

/*     REFERENCES */

/*     [1] Green, M. and Limebeer, D.J.N. */
/*         Linear Robust Control. */
/*         Prentice-Hall, Englewood Cliffs, NJ, 1995. */

/*     [2] Petkov, P.Hr., Gu, D.W., and Konstantinov, M.M. */
/*         Fortran 77 routines for Hinf and H2 design of linear */
/*         discrete-time control systems. */
/*         Report 99-8, Department of Engineering, Leicester University, */
/*         April 1999. */

/*     NUMERICAL ASPECTS */

/*     With approaching the minimum value of gamma some of the matrices */
/*     which are to be inverted tend to become ill-conditioned and */
/*     the X- or Z-Riccati equation may also become ill-conditioned */
/*     which may deteriorate the accuracy of the result. (The */
/*     corresponding reciprocal condition numbers are given in */
/*     the output array RCOND.) */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, April 1999. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Sep. 1999. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2000. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, discrete-time H-infinity optimal */
/*     control, robust control. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */

/*     .. External Functions */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test input parameters. */

#line 329 "SB10DD.f"
    /* Parameter adjustments */
#line 329 "SB10DD.f"
    a_dim1 = *lda;
#line 329 "SB10DD.f"
    a_offset = 1 + a_dim1;
#line 329 "SB10DD.f"
    a -= a_offset;
#line 329 "SB10DD.f"
    b_dim1 = *ldb;
#line 329 "SB10DD.f"
    b_offset = 1 + b_dim1;
#line 329 "SB10DD.f"
    b -= b_offset;
#line 329 "SB10DD.f"
    c_dim1 = *ldc;
#line 329 "SB10DD.f"
    c_offset = 1 + c_dim1;
#line 329 "SB10DD.f"
    c__ -= c_offset;
#line 329 "SB10DD.f"
    d_dim1 = *ldd;
#line 329 "SB10DD.f"
    d_offset = 1 + d_dim1;
#line 329 "SB10DD.f"
    d__ -= d_offset;
#line 329 "SB10DD.f"
    ak_dim1 = *ldak;
#line 329 "SB10DD.f"
    ak_offset = 1 + ak_dim1;
#line 329 "SB10DD.f"
    ak -= ak_offset;
#line 329 "SB10DD.f"
    bk_dim1 = *ldbk;
#line 329 "SB10DD.f"
    bk_offset = 1 + bk_dim1;
#line 329 "SB10DD.f"
    bk -= bk_offset;
#line 329 "SB10DD.f"
    ck_dim1 = *ldck;
#line 329 "SB10DD.f"
    ck_offset = 1 + ck_dim1;
#line 329 "SB10DD.f"
    ck -= ck_offset;
#line 329 "SB10DD.f"
    dk_dim1 = *lddk;
#line 329 "SB10DD.f"
    dk_offset = 1 + dk_dim1;
#line 329 "SB10DD.f"
    dk -= dk_offset;
#line 329 "SB10DD.f"
    x_dim1 = *ldx;
#line 329 "SB10DD.f"
    x_offset = 1 + x_dim1;
#line 329 "SB10DD.f"
    x -= x_offset;
#line 329 "SB10DD.f"
    z_dim1 = *ldz;
#line 329 "SB10DD.f"
    z_offset = 1 + z_dim1;
#line 329 "SB10DD.f"
    z__ -= z_offset;
#line 329 "SB10DD.f"
    --rcond;
#line 329 "SB10DD.f"
    --iwork;
#line 329 "SB10DD.f"
    --dwork;
#line 329 "SB10DD.f"
    --bwork;
#line 329 "SB10DD.f"

#line 329 "SB10DD.f"
    /* Function Body */
#line 329 "SB10DD.f"
    m1 = *m - *ncon;
#line 330 "SB10DD.f"
    m2 = *ncon;
#line 331 "SB10DD.f"
    np1 = *np - *nmeas;
#line 332 "SB10DD.f"
    np2 = *nmeas;

#line 334 "SB10DD.f"
    *info = 0;
#line 335 "SB10DD.f"
    if (*n < 0) {
#line 336 "SB10DD.f"
	*info = -1;
#line 337 "SB10DD.f"
    } else if (*m < 0) {
#line 338 "SB10DD.f"
	*info = -2;
#line 339 "SB10DD.f"
    } else if (*np < 0) {
#line 340 "SB10DD.f"
	*info = -3;
#line 341 "SB10DD.f"
    } else if (*ncon < 0 || m1 < 0 || m2 > np1) {
#line 342 "SB10DD.f"
	*info = -4;
#line 343 "SB10DD.f"
    } else if (*nmeas < 0 || np1 < 0 || np2 > m1) {
#line 344 "SB10DD.f"
	*info = -5;
#line 345 "SB10DD.f"
    } else if (*gamma <= 0.) {
#line 346 "SB10DD.f"
	*info = -6;
#line 347 "SB10DD.f"
    } else if (*lda < max(1,*n)) {
#line 348 "SB10DD.f"
	*info = -8;
#line 349 "SB10DD.f"
    } else if (*ldb < max(1,*n)) {
#line 350 "SB10DD.f"
	*info = -10;
#line 351 "SB10DD.f"
    } else if (*ldc < max(1,*np)) {
#line 352 "SB10DD.f"
	*info = -12;
#line 353 "SB10DD.f"
    } else if (*ldd < max(1,*np)) {
#line 354 "SB10DD.f"
	*info = -14;
#line 355 "SB10DD.f"
    } else if (*ldak < max(1,*n)) {
#line 356 "SB10DD.f"
	*info = -16;
#line 357 "SB10DD.f"
    } else if (*ldbk < max(1,*n)) {
#line 358 "SB10DD.f"
	*info = -18;
#line 359 "SB10DD.f"
    } else if (*ldck < max(1,m2)) {
#line 360 "SB10DD.f"
	*info = -20;
#line 361 "SB10DD.f"
    } else if (*lddk < max(1,m2)) {
#line 362 "SB10DD.f"
	*info = -22;
#line 363 "SB10DD.f"
    } else if (*ldx < max(1,*n)) {
#line 364 "SB10DD.f"
	*info = -24;
#line 365 "SB10DD.f"
    } else if (*ldz < max(1,*n)) {
#line 366 "SB10DD.f"
	*info = -26;
#line 367 "SB10DD.f"
    } else {

/*        Compute workspace. */

/* Computing MAX */
#line 371 "SB10DD.f"
	i__1 = (*n + m2) * 3 + *n + np1, i__2 = (*n + m2) * 5;
#line 371 "SB10DD.f"
	iwb = (*n + np1 + 1) * (*n + m2) + max(i__1,i__2);
/* Computing MAX */
#line 373 "SB10DD.f"
	i__1 = (*n + np2) * 3 + *n + m1, i__2 = (*n + np2) * 5;
#line 373 "SB10DD.f"
	iwc = (*n + np2) * (*n + m1 + 1) + max(i__1,i__2);
/* Computing MAX */
#line 375 "SB10DD.f"
	i__1 = *n * 14 + 23, i__2 = *n << 4, i__1 = max(i__1,i__2), i__2 = (*
		n << 1) + *m, i__1 = max(i__1,i__2), i__2 = *m * 3;
#line 375 "SB10DD.f"
	iwd = *n * 13 * *n + (*m << 1) * *m + *n * ((*m << 3) + np2) + m1 * (
		m2 + np2) + *n * 6 + max(i__1,i__2);
/* Computing MAX */
#line 377 "SB10DD.f"
	i__1 = *n * 14 + 23, i__2 = *n << 4, i__1 = max(i__1,i__2), i__2 = (*
		n << 1) + m2 + np2, i__1 = max(i__1,i__2), i__2 = (m2 + np2) *
		 3;
#line 377 "SB10DD.f"
	iwg = *n * 13 * *n + *m * *m + ((*n << 3) + *m + m2 + (np2 << 1)) * (
		m2 + np2) + *n * 6 + *n * (*m + np2) + max(i__1,i__2);
/* Computing MAX */
#line 380 "SB10DD.f"
	i__1 = max(iwb,iwc), i__1 = max(i__1,iwd);
#line 380 "SB10DD.f"
	minwrk = max(i__1,iwg);
#line 381 "SB10DD.f"
	if (*ldwork < minwrk) {
#line 381 "SB10DD.f"
	    *info = -31;
#line 381 "SB10DD.f"
	}
#line 383 "SB10DD.f"
    }
#line 384 "SB10DD.f"
    if (*info != 0) {
#line 385 "SB10DD.f"
	i__1 = -(*info);
#line 385 "SB10DD.f"
	xerbla_("SB10DD", &i__1, (ftnlen)6);
#line 386 "SB10DD.f"
	return 0;
#line 387 "SB10DD.f"
    }

/*     Quick return if possible. */

#line 391 "SB10DD.f"
    if (*n == 0 || *m == 0 || *np == 0 || m1 == 0 || m2 == 0 || np1 == 0 || 
	    np2 == 0) {
#line 393 "SB10DD.f"
	rcond[1] = 1.;
#line 394 "SB10DD.f"
	rcond[2] = 1.;
#line 395 "SB10DD.f"
	rcond[3] = 1.;
#line 396 "SB10DD.f"
	rcond[4] = 1.;
#line 397 "SB10DD.f"
	rcond[5] = 1.;
#line 398 "SB10DD.f"
	rcond[6] = 1.;
#line 399 "SB10DD.f"
	rcond[7] = 1.;
#line 400 "SB10DD.f"
	rcond[8] = 1.;
#line 401 "SB10DD.f"
	dwork[1] = 1.;
#line 402 "SB10DD.f"
	return 0;
#line 403 "SB10DD.f"
    }

#line 405 "SB10DD.f"
    toll = *tol;
#line 406 "SB10DD.f"
    if (toll <= 0.) {

/*        Set the default value of the tolerance in rank determination. */

#line 410 "SB10DD.f"
	toll = dlamch_("Epsilon", (ftnlen)7) * 1e3;
#line 411 "SB10DD.f"
    }

/*     Workspace usage. */

#line 415 "SB10DD.f"
    iws = (*n + np1) * (*n + m2) + 1;
#line 416 "SB10DD.f"
    iwrk = iws + (*n + m2);

/*                      jTheta */
/*     Determine if |A-e      I  B2 | has full column rank at */
/*                  |     C1     D12| */
/*     Theta = Pi/2 . */
/*     Workspace: need (N+NP1+1)*(N+M2) + MAX(3*(N+M2)+N+NP1,5*(N+M2)); */
/*                prefer larger. */

#line 425 "SB10DD.f"
    i__1 = *n + np1;
#line 425 "SB10DD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[1], &i__1, (ftnlen)4);
#line 426 "SB10DD.f"
    i__1 = *n + np1;
#line 426 "SB10DD.f"
    dlacpy_("Full", &np1, n, &c__[c_offset], ldc, &dwork[*n + 1], &i__1, (
	    ftnlen)4);
#line 427 "SB10DD.f"
    i__1 = *n + np1;
#line 427 "SB10DD.f"
    dlacpy_("Full", n, &m2, &b[(m1 + 1) * b_dim1 + 1], ldb, &dwork[(*n + np1) 
	    * *n + 1], &i__1, (ftnlen)4);
#line 429 "SB10DD.f"
    i__1 = *n + np1;
#line 429 "SB10DD.f"
    dlacpy_("Full", &np1, &m2, &d__[(m1 + 1) * d_dim1 + 1], ldd, &dwork[(*n + 
	    np1) * *n + *n + 1], &i__1, (ftnlen)4);
#line 431 "SB10DD.f"
    i__1 = *n + np1;
#line 431 "SB10DD.f"
    i__2 = *n + m2;
#line 431 "SB10DD.f"
    i__3 = *n + np1;
#line 431 "SB10DD.f"
    i__4 = *n + np1;
#line 431 "SB10DD.f"
    i__5 = *n + m2;
#line 431 "SB10DD.f"
    i__6 = *ldwork - iwrk + 1;
#line 431 "SB10DD.f"
    dgesvd_("N", "N", &i__1, &i__2, &dwork[1], &i__3, &dwork[iws], &dwork[1], 
	    &i__4, &dwork[1], &i__5, &dwork[iwrk], &i__6, &info2, (ftnlen)1, (
	    ftnlen)1);
#line 434 "SB10DD.f"
    if (info2 > 0) {
#line 435 "SB10DD.f"
	*info = 9;
#line 436 "SB10DD.f"
	return 0;
#line 437 "SB10DD.f"
    }
#line 438 "SB10DD.f"
    if (dwork[iws + *n + m2] / dwork[iws] <= toll) {
#line 439 "SB10DD.f"
	*info = 1;
#line 440 "SB10DD.f"
	return 0;
#line 441 "SB10DD.f"
    }
#line 442 "SB10DD.f"
    lwamax = (integer) dwork[iwrk] + iwrk - 1;

/*     Workspace usage. */

#line 446 "SB10DD.f"
    iws = (*n + np2) * (*n + m1) + 1;
#line 447 "SB10DD.f"
    iwrk = iws + (*n + np2);

/*                      jTheta */
/*     Determine if |A-e      I   B1 | has full row rank at */
/*                  |     C2      D21| */
/*     Theta = Pi/2 . */
/*     Workspace: need   (N+NP2)*(N+M1+1) + */
/*                       MAX(3*(N+NP2)+N+M1,5*(N+NP2)); */
/*                prefer larger. */

#line 457 "SB10DD.f"
    i__1 = *n + np2;
#line 457 "SB10DD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[1], &i__1, (ftnlen)4);
#line 458 "SB10DD.f"
    i__1 = *n + np2;
#line 458 "SB10DD.f"
    dlacpy_("Full", &np2, n, &c__[np1 + 1 + c_dim1], ldc, &dwork[*n + 1], &
	    i__1, (ftnlen)4);
#line 460 "SB10DD.f"
    i__1 = *n + np2;
#line 460 "SB10DD.f"
    dlacpy_("Full", n, &m1, &b[b_offset], ldb, &dwork[(*n + np2) * *n + 1], &
	    i__1, (ftnlen)4);
#line 462 "SB10DD.f"
    i__1 = *n + np2;
#line 462 "SB10DD.f"
    dlacpy_("Full", &np2, &m1, &d__[np1 + 1 + d_dim1], ldd, &dwork[(*n + np2) 
	    * *n + *n + 1], &i__1, (ftnlen)4);
#line 464 "SB10DD.f"
    i__1 = *n + np2;
#line 464 "SB10DD.f"
    i__2 = *n + m1;
#line 464 "SB10DD.f"
    i__3 = *n + np2;
#line 464 "SB10DD.f"
    i__4 = *n + np2;
#line 464 "SB10DD.f"
    i__5 = *n + m1;
#line 464 "SB10DD.f"
    i__6 = *ldwork - iwrk + 1;
#line 464 "SB10DD.f"
    dgesvd_("N", "N", &i__1, &i__2, &dwork[1], &i__3, &dwork[iws], &dwork[1], 
	    &i__4, &dwork[1], &i__5, &dwork[iwrk], &i__6, &info2, (ftnlen)1, (
	    ftnlen)1);
#line 467 "SB10DD.f"
    if (info2 > 0) {
#line 468 "SB10DD.f"
	*info = 9;
#line 469 "SB10DD.f"
	return 0;
#line 470 "SB10DD.f"
    }
#line 471 "SB10DD.f"
    if (dwork[iws + *n + np2] / dwork[iws] <= toll) {
#line 472 "SB10DD.f"
	*info = 2;
#line 473 "SB10DD.f"
	return 0;
#line 474 "SB10DD.f"
    }
/* Computing MAX */
#line 475 "SB10DD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 475 "SB10DD.f"
    lwamax = max(i__1,lwamax);

/*     Workspace usage. */

#line 479 "SB10DD.f"
    iws = np1 * m2 + 1;
#line 480 "SB10DD.f"
    iwrk = iws + m2;

/*     Determine if D12 has full column rank. */
/*     Workspace: need   (NP1+1)*M2 + MAX(3*M2+NP1,5*M2); */
/*                prefer larger. */

#line 486 "SB10DD.f"
    dlacpy_("Full", &np1, &m2, &d__[(m1 + 1) * d_dim1 + 1], ldd, &dwork[1], &
	    np1, (ftnlen)4);
#line 487 "SB10DD.f"
    i__1 = *ldwork - iwrk + 1;
#line 487 "SB10DD.f"
    dgesvd_("N", "N", &np1, &m2, &dwork[1], &np1, &dwork[iws], &dwork[1], &
	    np1, &dwork[1], &m2, &dwork[iwrk], &i__1, &info2, (ftnlen)1, (
	    ftnlen)1);
#line 489 "SB10DD.f"
    if (info2 > 0) {
#line 490 "SB10DD.f"
	*info = 9;
#line 491 "SB10DD.f"
	return 0;
#line 492 "SB10DD.f"
    }
#line 493 "SB10DD.f"
    if (dwork[iws + m2] / dwork[iws] <= toll) {
#line 494 "SB10DD.f"
	*info = 3;
#line 495 "SB10DD.f"
	return 0;
#line 496 "SB10DD.f"
    }
/* Computing MAX */
#line 497 "SB10DD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 497 "SB10DD.f"
    lwamax = max(i__1,lwamax);

/*     Workspace usage. */

#line 501 "SB10DD.f"
    iws = np2 * m1 + 1;
#line 502 "SB10DD.f"
    iwrk = iws + np2;

/*     Determine if D21 has full row rank. */
/*     Workspace: need   NP2*(M1+1) + MAX(3*NP2+M1,5*NP2); */
/*                prefer larger. */

#line 508 "SB10DD.f"
    dlacpy_("Full", &np2, &m1, &d__[np1 + 1 + d_dim1], ldd, &dwork[1], &np2, (
	    ftnlen)4);
#line 509 "SB10DD.f"
    i__1 = *ldwork - iwrk + 1;
#line 509 "SB10DD.f"
    dgesvd_("N", "N", &np2, &m1, &dwork[1], &np2, &dwork[iws], &dwork[1], &
	    np2, &dwork[1], &m1, &dwork[iwrk], &i__1, &info2, (ftnlen)1, (
	    ftnlen)1);
#line 511 "SB10DD.f"
    if (info2 > 0) {
#line 512 "SB10DD.f"
	*info = 9;
#line 513 "SB10DD.f"
	return 0;
#line 514 "SB10DD.f"
    }
#line 515 "SB10DD.f"
    if (dwork[iws + np2] / dwork[iws] <= toll) {
#line 516 "SB10DD.f"
	*info = 4;
#line 517 "SB10DD.f"
	return 0;
#line 518 "SB10DD.f"
    }
/* Computing MAX */
#line 519 "SB10DD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 519 "SB10DD.f"
    lwamax = max(i__1,lwamax);

/*     Workspace usage. */

#line 523 "SB10DD.f"
    iwv = 1;
#line 524 "SB10DD.f"
    iwb = iwv + *m * *m;
#line 525 "SB10DD.f"
    iwc = iwb + *n * m1;
#line 526 "SB10DD.f"
    iwd = iwc + (m2 + np2) * *n;
#line 527 "SB10DD.f"
    iwq = iwd + (m2 + np2) * m1;
#line 528 "SB10DD.f"
    iwl = iwq + *n * *n;
#line 529 "SB10DD.f"
    iwr = iwl + *n * *m;
#line 530 "SB10DD.f"
    iwi = iwr + (*n << 1);
#line 531 "SB10DD.f"
    iwh = iwi + (*n << 1);
#line 532 "SB10DD.f"
    iws = iwh + (*n << 1);
#line 533 "SB10DD.f"
    iwt = iws + ((*n << 1) + *m) * ((*n << 1) + *m);
#line 534 "SB10DD.f"
    iwu = iwt + ((*n << 1) + *m << 1) * *n;
#line 535 "SB10DD.f"
    iwrk = iwu + (*n << 2) * *n;
#line 536 "SB10DD.f"
    ir2 = iwv + m1;
#line 537 "SB10DD.f"
    ir3 = ir2 + *m * m1;

/*     Compute R0 = |D11'||D11 D12| -|gamma^2*Im1 0| . */
/*                  |D12'|           |   0        0| */

#line 542 "SB10DD.f"
    dsyrk_("Lower", "Transpose", m, &np1, &c_b24, &d__[d_offset], ldd, &c_b25,
	     &dwork[1], m, (ftnlen)5, (ftnlen)9);
#line 544 "SB10DD.f"
    i__1 = *m * m1;
#line 544 "SB10DD.f"
    i__2 = *m + 1;
#line 544 "SB10DD.f"
    for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
#line 545 "SB10DD.f"
	dwork[j] -= *gamma * *gamma;
#line 546 "SB10DD.f"
/* L10: */
#line 546 "SB10DD.f"
    }

/*     Compute C1'*C1 . */

#line 550 "SB10DD.f"
    dsyrk_("Lower", "Transpose", n, &np1, &c_b24, &c__[c_offset], ldc, &c_b25,
	     &dwork[iwq], n, (ftnlen)5, (ftnlen)9);

/*     Compute C1'*|D11 D12| . */

#line 555 "SB10DD.f"
    dgemm_("Transpose", "NoTranspose", n, m, &np1, &c_b24, &c__[c_offset], 
	    ldc, &d__[d_offset], ldd, &c_b25, &dwork[iwl], n, (ftnlen)9, (
	    ftnlen)11);

/*     Solution of the X-Riccati equation. */
/*     Workspace: need   13*N*N + 2*M*M + N*(8*M+NP2) + M1*(M2+NP2) + */
/*                       6*N + max(14*N+23,16*N,2*N+M,3*M); */
/*                prefer larger. */

#line 563 "SB10DD.f"
    i__2 = (*n << 1) + *m;
#line 563 "SB10DD.f"
    i__1 = (*n << 1) + *m;
#line 563 "SB10DD.f"
    i__3 = *n << 1;
#line 563 "SB10DD.f"
    i__4 = *ldwork - iwrk + 1;
#line 563 "SB10DD.f"
    sb02od_("D", "B", "N", "L", "N", "S", n, m, np, &a[a_offset], lda, &b[
	    b_offset], ldb, &dwork[iwq], n, &dwork[1], m, &dwork[iwl], n, &
	    rcond2, &x[x_offset], ldx, &dwork[iwr], &dwork[iwi], &dwork[iwh], 
	    &dwork[iws], &i__2, &dwork[iwt], &i__1, &dwork[iwu], &i__3, &toll,
	     &iwork[1], &dwork[iwrk], &i__4, &bwork[1], &info2, (ftnlen)1, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 569 "SB10DD.f"
    if (info2 > 0) {
#line 570 "SB10DD.f"
	*info = 6;
#line 571 "SB10DD.f"
	return 0;
#line 572 "SB10DD.f"
    }
/* Computing MAX */
#line 573 "SB10DD.f"
    i__2 = (integer) dwork[iwrk] + iwrk - 1;
#line 573 "SB10DD.f"
    lwamax = max(i__2,lwamax);

/*     Condition estimation. */
/*     Workspace: need   4*N*N + 2*M*M + N*(3*M+NP2) + M1*(M2+NP2) + */
/*                       max(5*N,max(3,2*N*N)+N*N); */
/*                prefer larger. */

#line 580 "SB10DD.f"
    iws = iwr;
#line 581 "SB10DD.f"
    iwh = iws + *m * *m;
#line 582 "SB10DD.f"
    iwt = iwh + *n * *m;
#line 583 "SB10DD.f"
    iwu = iwt + *n * *n;
#line 584 "SB10DD.f"
    iwg = iwu + *n * *n;
#line 585 "SB10DD.f"
    iwrk = iwg + *n * *n;
#line 586 "SB10DD.f"
    dlacpy_("Lower", m, m, &dwork[1], m, &dwork[iws], m, (ftnlen)5);
#line 587 "SB10DD.f"
    i__2 = *ldwork - iwrk + 1;
#line 587 "SB10DD.f"
    dsytrf_("Lower", m, &dwork[iws], m, &iwork[1], &dwork[iwrk], &i__2, &
	    info2, (ftnlen)5);
#line 589 "SB10DD.f"
    if (info2 > 0) {
#line 590 "SB10DD.f"
	*info = 5;
#line 591 "SB10DD.f"
	return 0;
#line 592 "SB10DD.f"
    }
/* Computing MAX */
#line 593 "SB10DD.f"
    i__2 = (integer) dwork[iwrk] + iwrk - 1;
#line 593 "SB10DD.f"
    lwamax = max(i__2,lwamax);

#line 595 "SB10DD.f"
    ma02ad_("Full", n, m, &b[b_offset], ldb, &dwork[iwh], m, (ftnlen)4);
#line 596 "SB10DD.f"
    dsytrs_("Lower", m, n, &dwork[iws], m, &iwork[1], &dwork[iwh], m, &info2, 
	    (ftnlen)5);
#line 598 "SB10DD.f"
    mb01rx_("Left", "Lower", "NoTranspose", n, m, &c_b25, &c_b24, &dwork[iwg],
	     n, &b[b_offset], ldb, &dwork[iwh], m, &info2, (ftnlen)4, (ftnlen)
	    5, (ftnlen)11);
#line 600 "SB10DD.f"
    i__2 = *ldwork - iwrk + 1;
#line 600 "SB10DD.f"
    sb02sd_("C", "N", "N", "L", "O", n, &a[a_offset], lda, &dwork[iwt], n, &
	    dwork[iwu], n, &dwork[iwg], n, &dwork[iwq], n, &x[x_offset], ldx, 
	    &sepd, &rcond[7], &ferr, &iwork[1], &dwork[iwrk], &i__2, &info2, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 604 "SB10DD.f"
    if (info2 > 0) {
#line 604 "SB10DD.f"
	rcond[7] = 0.;
#line 604 "SB10DD.f"
    }
/* Computing MAX */
#line 605 "SB10DD.f"
    i__2 = (integer) dwork[iwrk] + iwrk - 1;
#line 605 "SB10DD.f"
    lwamax = max(i__2,lwamax);

/*     Workspace usage. */

#line 609 "SB10DD.f"
    iwrk = iwr;

/*     Compute the lower triangle of |R1  R2'| = R0 + B'*X*B . */
/*                                   |R2  R3 | */

#line 614 "SB10DD.f"
    i__2 = *m * *n;
#line 614 "SB10DD.f"
    mb01ru_("Lower", "Transpose", m, n, &c_b24, &c_b24, &dwork[1], m, &b[
	    b_offset], ldb, &x[x_offset], ldx, &dwork[iwrk], &i__2, &info2, (
	    ftnlen)5, (ftnlen)9);

/*     Compute the Cholesky factorization of R3, R3 = V12'*V12 . */
/*     Note that V12' is stored. */

#line 620 "SB10DD.f"
    anorm = dlansy_("1", "Lower", &m2, &dwork[ir3], m, &dwork[iwrk], (ftnlen)
	    1, (ftnlen)5);
#line 621 "SB10DD.f"
    dpotrf_("Lower", &m2, &dwork[ir3], m, &info2, (ftnlen)5);
#line 622 "SB10DD.f"
    if (info2 > 0) {
#line 623 "SB10DD.f"
	*info = 5;
#line 624 "SB10DD.f"
	return 0;
#line 625 "SB10DD.f"
    }
#line 626 "SB10DD.f"
    dpocon_("Lower", &m2, &dwork[ir3], m, &anorm, &rcond[1], &dwork[iwrk], &
	    iwork[1], &info2, (ftnlen)5);

/*     Return if the matrix is singular to working precision. */

#line 631 "SB10DD.f"
    if (rcond[1] < toll) {
#line 632 "SB10DD.f"
	*info = 5;
#line 633 "SB10DD.f"
	return 0;
#line 634 "SB10DD.f"
    }

#line 636 "SB10DD.f"
    dtrcon_("1", "Lower", "NonUnit", &m2, &dwork[ir3], m, &rcond[5], &dwork[
	    iwrk], &iwork[1], &info2, (ftnlen)1, (ftnlen)5, (ftnlen)7);

/*     Return if the matrix is singular to working precision. */

#line 641 "SB10DD.f"
    if (rcond[5] < toll) {
#line 642 "SB10DD.f"
	*info = 5;
#line 643 "SB10DD.f"
	return 0;
#line 644 "SB10DD.f"
    }

/*     Compute R2 <- inv(V12')*R2 . */

#line 648 "SB10DD.f"
    dtrsm_("Left", "Lower", "NoTranspose", "NonUnit", &m2, &m1, &c_b24, &
	    dwork[ir3], m, &dwork[ir2], m, (ftnlen)4, (ftnlen)5, (ftnlen)11, (
	    ftnlen)7);

/*     Compute -Nabla = R2'*inv(R3)*R2 - R1 . */

#line 653 "SB10DD.f"
    dsyrk_("Lower", "Transpose", &m1, &m2, &c_b24, &dwork[ir2], m, &c_b74, &
	    dwork[1], m, (ftnlen)5, (ftnlen)9);

/*     Compute the Cholesky factorization of -Nabla, -Nabla = V21t'*V21t. */
/*     Note that V21t' is stored. */

#line 659 "SB10DD.f"
    anorm = dlansy_("1", "Lower", &m1, &dwork[1], m, &dwork[iwrk], (ftnlen)1, 
	    (ftnlen)5);
#line 660 "SB10DD.f"
    dpotrf_("Lower", &m1, &dwork[1], m, &info2, (ftnlen)5);
#line 661 "SB10DD.f"
    if (info2 > 0) {
#line 662 "SB10DD.f"
	*info = 5;
#line 663 "SB10DD.f"
	return 0;
#line 664 "SB10DD.f"
    }
#line 665 "SB10DD.f"
    dpocon_("Lower", &m1, &dwork[1], m, &anorm, &rcond[2], &dwork[iwrk], &
	    iwork[1], &info2, (ftnlen)5);

/*     Return if the matrix is singular to working precision. */

#line 670 "SB10DD.f"
    if (rcond[2] < toll) {
#line 671 "SB10DD.f"
	*info = 5;
#line 672 "SB10DD.f"
	return 0;
#line 673 "SB10DD.f"
    }

#line 675 "SB10DD.f"
    dtrcon_("1", "Lower", "NonUnit", &m1, &dwork[1], m, &rcond[3], &dwork[
	    iwrk], &iwork[1], &info2, (ftnlen)1, (ftnlen)5, (ftnlen)7);

/*     Return if the matrix is singular to working precision. */

#line 680 "SB10DD.f"
    if (rcond[3] < toll) {
#line 681 "SB10DD.f"
	*info = 5;
#line 682 "SB10DD.f"
	return 0;
#line 683 "SB10DD.f"
    }

/*     Compute X*A . */

#line 687 "SB10DD.f"
    dgemm_("NoTranspose", "NoTranspose", n, n, n, &c_b24, &x[x_offset], ldx, &
	    a[a_offset], lda, &c_b25, &dwork[iwq], n, (ftnlen)11, (ftnlen)11);

/*     Compute |L1| = |D11'|*C1 + B'*X*A . */
/*             |L2| = |D12'| */

#line 693 "SB10DD.f"
    ma02ad_("Full", n, m, &dwork[iwl], n, &dwork[iwrk], m, (ftnlen)4);
#line 694 "SB10DD.f"
    dlacpy_("Full", m, n, &dwork[iwrk], m, &dwork[iwl], m, (ftnlen)4);
#line 695 "SB10DD.f"
    dgemm_("Transpose", "NoTranspose", m, n, n, &c_b24, &b[b_offset], ldb, &
	    dwork[iwq], n, &c_b24, &dwork[iwl], m, (ftnlen)9, (ftnlen)11);

/*     Compute L2 <- inv(V12')*L2 . */

#line 700 "SB10DD.f"
    dtrsm_("Left", "Lower", "NoTranspose", "NonUnit", &m2, n, &c_b24, &dwork[
	    ir3], m, &dwork[iwl + m1], m, (ftnlen)4, (ftnlen)5, (ftnlen)11, (
	    ftnlen)7);

/*     Compute L_Nabla = L1 - R2'*inv(R3)*L2 . */

#line 705 "SB10DD.f"
    dgemm_("Transpose", "NoTranspose", &m1, n, &m2, &c_b74, &dwork[ir2], m, &
	    dwork[iwl + m1], m, &c_b24, &dwork[iwl], m, (ftnlen)9, (ftnlen)11)
	    ;

/*     Compute L_Nabla <- inv(V21t')*L_Nabla . */

#line 711 "SB10DD.f"
    dtrsm_("Left", "Lower", "NoTranspose", "NonUnit", &m1, n, &c_b24, &dwork[
	    1], m, &dwork[iwl], m, (ftnlen)4, (ftnlen)5, (ftnlen)11, (ftnlen)
	    7);

/*     Compute Bt1 = B1*inv(V21t) . */

#line 716 "SB10DD.f"
    dlacpy_("Full", n, &m1, &b[b_offset], ldb, &dwork[iwb], n, (ftnlen)4);
#line 717 "SB10DD.f"
    dtrsm_("Right", "Lower", "Transpose", "NonUnit", n, &m1, &c_b24, &dwork[1]
	    , m, &dwork[iwb], n, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)7);

/*     Compute At . */

#line 722 "SB10DD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &ak[ak_offset], ldak, (ftnlen)4);
#line 723 "SB10DD.f"
    dgemm_("NoTranspose", "NoTranspose", n, n, &m1, &c_b24, &dwork[iwb], n, &
	    dwork[iwl], m, &c_b24, &ak[ak_offset], ldak, (ftnlen)11, (ftnlen)
	    11);

/*     Scale Bt1 . */

#line 728 "SB10DD.f"
    i__2 = *n * m1;
#line 728 "SB10DD.f"
    dscal_(&i__2, gamma, &dwork[iwb], &c__1);

/*     Compute |Dt11| = |R2 |*inv(V21t) . */
/*             |Dt21|   |D21| */

#line 733 "SB10DD.f"
    i__2 = m2 + np2;
#line 733 "SB10DD.f"
    dlacpy_("Full", &m2, &m1, &dwork[ir2], m, &dwork[iwd], &i__2, (ftnlen)4);
#line 735 "SB10DD.f"
    i__2 = m2 + np2;
#line 735 "SB10DD.f"
    dlacpy_("Full", &np2, &m1, &d__[np1 + 1 + d_dim1], ldd, &dwork[iwd + m2], 
	    &i__2, (ftnlen)4);
#line 737 "SB10DD.f"
    i__2 = m2 + np2;
#line 737 "SB10DD.f"
    i__1 = m2 + np2;
#line 737 "SB10DD.f"
    dtrsm_("Right", "Lower", "Transpose", "NonUnit", &i__2, &m1, &c_b24, &
	    dwork[1], m, &dwork[iwd], &i__1, (ftnlen)5, (ftnlen)5, (ftnlen)9, 
	    (ftnlen)7);

/*     Compute Ct = |Ct1| = |L2| + |Dt11|*inv(V21t')*L_Nabla . */
/*                  |Ct2| = |C2| + |Dt21| */

#line 743 "SB10DD.f"
    i__2 = m2 + np2;
#line 743 "SB10DD.f"
    dlacpy_("Full", &m2, n, &dwork[iwl + m1], m, &dwork[iwc], &i__2, (ftnlen)
	    4);
#line 745 "SB10DD.f"
    i__2 = m2 + np2;
#line 745 "SB10DD.f"
    dlacpy_("Full", &np2, n, &c__[np1 + 1 + c_dim1], ldc, &dwork[iwc + m2], &
	    i__2, (ftnlen)4);
#line 747 "SB10DD.f"
    i__2 = m2 + np2;
#line 747 "SB10DD.f"
    i__1 = m2 + np2;
#line 747 "SB10DD.f"
    i__3 = m2 + np2;
#line 747 "SB10DD.f"
    dgemm_("NoTranspose", "NoTranspose", &i__2, n, &m1, &c_b24, &dwork[iwd], &
	    i__1, &dwork[iwl], m, &c_b24, &dwork[iwc], &i__3, (ftnlen)11, (
	    ftnlen)11);

/*     Scale |Dt11| . */
/*           |Dt21| */

#line 754 "SB10DD.f"
    i__2 = (m2 + np2) * m1;
#line 754 "SB10DD.f"
    dscal_(&i__2, gamma, &dwork[iwd], &c__1);

/*     Workspace usage. */

#line 758 "SB10DD.f"
    iww = iwd + (m2 + np2) * m1;
#line 759 "SB10DD.f"
    iwq = iww + (m2 + np2) * (m2 + np2);
#line 760 "SB10DD.f"
    iwl = iwq + *n * *n;
#line 761 "SB10DD.f"
    iwr = iwl + *n * (m2 + np2);
#line 762 "SB10DD.f"
    iwi = iwr + (*n << 1);
#line 763 "SB10DD.f"
    iwh = iwi + (*n << 1);
#line 764 "SB10DD.f"
    iws = iwh + (*n << 1);
#line 765 "SB10DD.f"
    iwt = iws + ((*n << 1) + m2 + np2) * ((*n << 1) + m2 + np2);
#line 766 "SB10DD.f"
    iwu = iwt + ((*n << 1) + m2 + np2 << 1) * *n;
#line 767 "SB10DD.f"
    iwg = iwu + (*n << 2) * *n;
#line 768 "SB10DD.f"
    iwrk = iwg + (m2 + np2) * *n;
#line 769 "SB10DD.f"
    is2 = iww + (m2 + np2) * m2;
#line 770 "SB10DD.f"
    is3 = is2 + m2;

/*     Compute S0 = |Dt11||Dt11' Dt21'| -|gamma^2*Im2 0| . */
/*                  |Dt21|               |   0        0| */

#line 775 "SB10DD.f"
    i__2 = m2 + np2;
#line 775 "SB10DD.f"
    i__1 = m2 + np2;
#line 775 "SB10DD.f"
    i__3 = m2 + np2;
#line 775 "SB10DD.f"
    dsyrk_("Upper", "NoTranspose", &i__2, &m1, &c_b24, &dwork[iwd], &i__1, &
	    c_b25, &dwork[iww], &i__3, (ftnlen)5, (ftnlen)11);
#line 777 "SB10DD.f"
    i__2 = iww - 1 + (m2 + np2) * m2;
#line 777 "SB10DD.f"
    i__1 = m2 + np2 + 1;
#line 777 "SB10DD.f"
    for (j = iww; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
#line 778 "SB10DD.f"
	dwork[j] -= *gamma * *gamma;
#line 779 "SB10DD.f"
/* L20: */
#line 779 "SB10DD.f"
    }

/*     Compute Bt1*Bt1' . */

#line 783 "SB10DD.f"
    dsyrk_("Upper", "NoTranspose", n, &m1, &c_b24, &dwork[iwb], n, &c_b25, &
	    dwork[iwq], n, (ftnlen)5, (ftnlen)11);

/*     Compute Bt1*|Dt11' Dt21'| . */

#line 788 "SB10DD.f"
    i__1 = m2 + np2;
#line 788 "SB10DD.f"
    i__2 = m2 + np2;
#line 788 "SB10DD.f"
    dgemm_("NoTranspose", "Transpose", n, &i__1, &m1, &c_b24, &dwork[iwb], n, 
	    &dwork[iwd], &i__2, &c_b25, &dwork[iwl], n, (ftnlen)11, (ftnlen)9)
	    ;

/*     Transpose At in situ (in AK) . */

#line 794 "SB10DD.f"
    i__1 = *n;
#line 794 "SB10DD.f"
    for (j = 2; j <= i__1; ++j) {
#line 795 "SB10DD.f"
	i__2 = j - 1;
#line 795 "SB10DD.f"
	dswap_(&i__2, &ak[j + ak_dim1], ldak, &ak[j * ak_dim1 + 1], &c__1);
#line 796 "SB10DD.f"
/* L30: */
#line 796 "SB10DD.f"
    }

/*     Transpose Ct . */

#line 800 "SB10DD.f"
    i__1 = m2 + np2;
#line 800 "SB10DD.f"
    i__2 = m2 + np2;
#line 800 "SB10DD.f"
    ma02ad_("Full", &i__1, n, &dwork[iwc], &i__2, &dwork[iwg], n, (ftnlen)4);

/*     Solution of the Z-Riccati equation. */
/*     Workspace: need   13*N*N + M*M + (8*N+M+M2+2*NP2)*(M2+NP2) + */
/*                       N*(M+NP2) + 6*N + */
/*                       max(14*N+23,16*N,2*N+M2+NP2,3*(M2+NP2)); */
/*                prefer larger. */

#line 809 "SB10DD.f"
    i__1 = m2 + np2;
#line 809 "SB10DD.f"
    i__2 = m2 + np2;
#line 809 "SB10DD.f"
    i__3 = (*n << 1) + m2 + np2;
#line 809 "SB10DD.f"
    i__4 = (*n << 1) + m2 + np2;
#line 809 "SB10DD.f"
    i__5 = *n << 1;
#line 809 "SB10DD.f"
    i__6 = *ldwork - iwrk + 1;
#line 809 "SB10DD.f"
    sb02od_("D", "B", "N", "U", "N", "S", n, &i__1, np, &ak[ak_offset], ldak, 
	    &dwork[iwg], n, &dwork[iwq], n, &dwork[iww], &i__2, &dwork[iwl], 
	    n, &rcond2, &z__[z_offset], ldz, &dwork[iwr], &dwork[iwi], &dwork[
	    iwh], &dwork[iws], &i__3, &dwork[iwt], &i__4, &dwork[iwu], &i__5, 
	    &toll, &iwork[1], &dwork[iwrk], &i__6, &bwork[1], &info2, (ftnlen)
	    1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 815 "SB10DD.f"
    if (info2 > 0) {
#line 816 "SB10DD.f"
	*info = 7;
#line 817 "SB10DD.f"
	return 0;
#line 818 "SB10DD.f"
    }
/* Computing MAX */
#line 819 "SB10DD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 819 "SB10DD.f"
    lwamax = max(i__1,lwamax);

/*     Condition estimation. */
/*     Workspace: need   4*N*N + M*M + 2*(M2+NP2)*(M2+NP2)+ */
/*                       N*(M+2*M2+3*NP2) + (M2+NP2)*M1 + */
/*                       max(5*N,max(3,2*N*N)+N*N); */
/*                prefer larger. */

#line 827 "SB10DD.f"
    iws = iwr;
#line 828 "SB10DD.f"
    iwh = iws + (m2 + np2) * (m2 + np2);
#line 829 "SB10DD.f"
    iwt = iwh + *n * (m2 + np2);
#line 830 "SB10DD.f"
    iwu = iwt + *n * *n;
#line 831 "SB10DD.f"
    iwg = iwu + *n * *n;
#line 832 "SB10DD.f"
    iwrk = iwg + *n * *n;
#line 833 "SB10DD.f"
    i__1 = m2 + np2;
#line 833 "SB10DD.f"
    i__2 = m2 + np2;
#line 833 "SB10DD.f"
    i__3 = m2 + np2;
#line 833 "SB10DD.f"
    i__4 = m2 + np2;
#line 833 "SB10DD.f"
    dlacpy_("Upper", &i__1, &i__2, &dwork[iww], &i__3, &dwork[iws], &i__4, (
	    ftnlen)5);
#line 835 "SB10DD.f"
    i__1 = m2 + np2;
#line 835 "SB10DD.f"
    i__2 = m2 + np2;
#line 835 "SB10DD.f"
    i__3 = *ldwork - iwrk + 1;
#line 835 "SB10DD.f"
    dsytrf_("Upper", &i__1, &dwork[iws], &i__2, &iwork[1], &dwork[iwrk], &
	    i__3, &info2, (ftnlen)5);
#line 837 "SB10DD.f"
    if (info2 > 0) {
#line 838 "SB10DD.f"
	*info = 5;
#line 839 "SB10DD.f"
	return 0;
#line 840 "SB10DD.f"
    }
/* Computing MAX */
#line 841 "SB10DD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 841 "SB10DD.f"
    lwamax = max(i__1,lwamax);

#line 843 "SB10DD.f"
    i__1 = m2 + np2;
#line 843 "SB10DD.f"
    i__2 = m2 + np2;
#line 843 "SB10DD.f"
    i__3 = m2 + np2;
#line 843 "SB10DD.f"
    dlacpy_("Full", &i__1, n, &dwork[iwc], &i__2, &dwork[iwh], &i__3, (ftnlen)
	    4);
#line 845 "SB10DD.f"
    i__1 = m2 + np2;
#line 845 "SB10DD.f"
    i__2 = m2 + np2;
#line 845 "SB10DD.f"
    i__3 = m2 + np2;
#line 845 "SB10DD.f"
    dsytrs_("Upper", &i__1, n, &dwork[iws], &i__2, &iwork[1], &dwork[iwh], &
	    i__3, &info2, (ftnlen)5);
#line 847 "SB10DD.f"
    i__1 = m2 + np2;
#line 847 "SB10DD.f"
    i__2 = m2 + np2;
#line 847 "SB10DD.f"
    i__3 = m2 + np2;
#line 847 "SB10DD.f"
    mb01rx_("Left", "Upper", "Transpose", n, &i__1, &c_b25, &c_b24, &dwork[
	    iwg], n, &dwork[iwc], &i__2, &dwork[iwh], &i__3, &info2, (ftnlen)
	    4, (ftnlen)5, (ftnlen)9);
#line 850 "SB10DD.f"
    i__1 = *ldwork - iwrk + 1;
#line 850 "SB10DD.f"
    sb02sd_("C", "N", "N", "U", "O", n, &ak[ak_offset], ldak, &dwork[iwt], n, 
	    &dwork[iwu], n, &dwork[iwg], n, &dwork[iwq], n, &z__[z_offset], 
	    ldz, &sepd, &rcond[8], &ferr, &iwork[1], &dwork[iwrk], &i__1, &
	    info2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 854 "SB10DD.f"
    if (info2 > 0) {
#line 854 "SB10DD.f"
	rcond[8] = 0.;
#line 854 "SB10DD.f"
    }
/* Computing MAX */
#line 855 "SB10DD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 855 "SB10DD.f"
    lwamax = max(i__1,lwamax);

/*     Workspace usage. */

#line 859 "SB10DD.f"
    iwrk = iwr;

/*     Compute the upper triangle of */
/*             |St1   St2| = S0 + |Ct1|*Z*|Ct1' Ct2'| . */
/*             |St2'  St3|        |Ct2| */

#line 865 "SB10DD.f"
    i__1 = m2 + np2;
#line 865 "SB10DD.f"
    i__2 = m2 + np2;
#line 865 "SB10DD.f"
    i__3 = m2 + np2;
#line 865 "SB10DD.f"
    i__4 = (m2 + np2) * *n;
#line 865 "SB10DD.f"
    mb01ru_("Upper", "NoTranspose", &i__1, n, &c_b24, &c_b24, &dwork[iww], &
	    i__2, &dwork[iwc], &i__3, &z__[z_offset], ldz, &dwork[iwrk], &
	    i__4, &info2, (ftnlen)5, (ftnlen)11);

/*     Compute the Cholesky factorization of St3, St3 = U12'*U12 . */

#line 871 "SB10DD.f"
    i__1 = m2 + np2;
#line 871 "SB10DD.f"
    anorm = dlansy_("1", "Upper", &np2, &dwork[is3], &i__1, &dwork[iwrk], (
	    ftnlen)1, (ftnlen)5);
#line 873 "SB10DD.f"
    i__1 = m2 + np2;
#line 873 "SB10DD.f"
    dpotrf_("Upper", &np2, &dwork[is3], &i__1, &info2, (ftnlen)5);
#line 874 "SB10DD.f"
    if (info2 > 0) {
#line 875 "SB10DD.f"
	*info = 5;
#line 876 "SB10DD.f"
	return 0;
#line 877 "SB10DD.f"
    }
#line 878 "SB10DD.f"
    i__1 = m2 + np2;
#line 878 "SB10DD.f"
    dpocon_("Upper", &np2, &dwork[is3], &i__1, &anorm, &rcond[4], &dwork[iwrk]
	    , &iwork[1], &info2, (ftnlen)5);

/*     Return if the matrix is singular to working precision. */

#line 883 "SB10DD.f"
    if (rcond[4] < toll) {
#line 884 "SB10DD.f"
	*info = 5;
#line 885 "SB10DD.f"
	return 0;
#line 886 "SB10DD.f"
    }

/*     Compute St2 <- St2*inv(U12) . */

#line 890 "SB10DD.f"
    i__1 = m2 + np2;
#line 890 "SB10DD.f"
    i__2 = m2 + np2;
#line 890 "SB10DD.f"
    dtrsm_("Right", "Upper", "NoTranspose", "NonUnit", &m2, &np2, &c_b24, &
	    dwork[is3], &i__1, &dwork[is2], &i__2, (ftnlen)5, (ftnlen)5, (
	    ftnlen)11, (ftnlen)7);

/*     Check the negative definiteness of St1 - St2*inv(St3)*St2' . */

#line 895 "SB10DD.f"
    i__1 = m2 + np2;
#line 895 "SB10DD.f"
    i__2 = m2 + np2;
#line 895 "SB10DD.f"
    dsyrk_("Upper", "NoTranspose", &m2, &np2, &c_b24, &dwork[is2], &i__1, &
	    c_b74, &dwork[iww], &i__2, (ftnlen)5, (ftnlen)11);
#line 897 "SB10DD.f"
    i__1 = m2 + np2;
#line 897 "SB10DD.f"
    dpotrf_("Upper", &m2, &dwork[iww], &i__1, &info2, (ftnlen)5);
#line 898 "SB10DD.f"
    if (info2 > 0) {
#line 899 "SB10DD.f"
	*info = 5;
#line 900 "SB10DD.f"
	return 0;
#line 901 "SB10DD.f"
    }

/*     Restore At in situ . */

#line 905 "SB10DD.f"
    i__1 = *n;
#line 905 "SB10DD.f"
    for (j = 2; j <= i__1; ++j) {
#line 906 "SB10DD.f"
	i__2 = j - 1;
#line 906 "SB10DD.f"
	dswap_(&i__2, &ak[j + ak_dim1], ldak, &ak[j * ak_dim1 + 1], &c__1);
#line 907 "SB10DD.f"
/* L40: */
#line 907 "SB10DD.f"
    }

/*     Compute At*Z . */

#line 911 "SB10DD.f"
    dgemm_("NoTranspose", "NoTranspose", n, n, n, &c_b24, &ak[ak_offset], 
	    ldak, &z__[z_offset], ldz, &c_b25, &dwork[iwrk], n, (ftnlen)11, (
	    ftnlen)11);

/*     Compute Mt2 = Bt1*Dt21' + At*Z*Ct2' in BK . */

#line 916 "SB10DD.f"
    dlacpy_("Full", n, &np2, &dwork[iwl + *n * m2], n, &bk[bk_offset], ldbk, (
	    ftnlen)4);
#line 917 "SB10DD.f"
    i__1 = m2 + np2;
#line 917 "SB10DD.f"
    dgemm_("NoTranspose", "Transpose", n, &np2, n, &c_b24, &dwork[iwrk], n, &
	    dwork[iwc + m2], &i__1, &c_b24, &bk[bk_offset], ldbk, (ftnlen)11, 
	    (ftnlen)9);

/*     Compute St2 <- St2*inv(U12') . */

#line 923 "SB10DD.f"
    i__1 = m2 + np2;
#line 923 "SB10DD.f"
    i__2 = m2 + np2;
#line 923 "SB10DD.f"
    dtrsm_("Right", "Upper", "Transpose", "NonUnit", &m2, &np2, &c_b24, &
	    dwork[is3], &i__1, &dwork[is2], &i__2, (ftnlen)5, (ftnlen)5, (
	    ftnlen)9, (ftnlen)7);

/*     Compute DKHAT = -inv(V12)*St2 in DK . */

#line 928 "SB10DD.f"
    i__1 = m2 + np2;
#line 928 "SB10DD.f"
    dlacpy_("Full", &m2, &np2, &dwork[is2], &i__1, &dk[dk_offset], lddk, (
	    ftnlen)4);
#line 929 "SB10DD.f"
    dtrsm_("Left", "Lower", "Transpose", "NonUnit", &m2, &np2, &c_b74, &dwork[
	    ir3], m, &dk[dk_offset], lddk, (ftnlen)4, (ftnlen)5, (ftnlen)9, (
	    ftnlen)7);

/*     Compute CKHAT = -inv(V12)*(Ct1 - St2*inv(St3)*Ct2) in CK . */

#line 934 "SB10DD.f"
    i__1 = m2 + np2;
#line 934 "SB10DD.f"
    dlacpy_("Full", &m2, n, &dwork[iwc], &i__1, &ck[ck_offset], ldck, (ftnlen)
	    4);
#line 935 "SB10DD.f"
    i__1 = m2 + np2;
#line 935 "SB10DD.f"
    i__2 = m2 + np2;
#line 935 "SB10DD.f"
    dgemm_("NoTranspose", "NoTranspose", &m2, n, &np2, &c_b74, &dwork[is2], &
	    i__1, &dwork[iwc + m2], &i__2, &c_b24, &ck[ck_offset], ldck, (
	    ftnlen)11, (ftnlen)11);
#line 938 "SB10DD.f"
    dtrsm_("Left", "Lower", "Transpose", "NonUnit", &m2, n, &c_b74, &dwork[
	    ir3], m, &ck[ck_offset], ldck, (ftnlen)4, (ftnlen)5, (ftnlen)9, (
	    ftnlen)7);

/*     Compute Mt2*inv(St3) in BK . */

#line 943 "SB10DD.f"
    i__1 = m2 + np2;
#line 943 "SB10DD.f"
    dtrsm_("Right", "Upper", "NoTranspose", "NonUnit", n, &np2, &c_b24, &
	    dwork[is3], &i__1, &bk[bk_offset], ldbk, (ftnlen)5, (ftnlen)5, (
	    ftnlen)11, (ftnlen)7);
#line 945 "SB10DD.f"
    i__1 = m2 + np2;
#line 945 "SB10DD.f"
    dtrsm_("Right", "Upper", "Transpose", "NonUnit", n, &np2, &c_b24, &dwork[
	    is3], &i__1, &bk[bk_offset], ldbk, (ftnlen)5, (ftnlen)5, (ftnlen)
	    9, (ftnlen)7);

/*     Compute AKHAT in AK . */

#line 950 "SB10DD.f"
    dgemm_("NoTranspose", "NoTranspose", n, n, &m2, &c_b24, &b[(m1 + 1) * 
	    b_dim1 + 1], ldb, &ck[ck_offset], ldck, &c_b24, &ak[ak_offset], 
	    ldak, (ftnlen)11, (ftnlen)11);
#line 952 "SB10DD.f"
    i__1 = m2 + np2;
#line 952 "SB10DD.f"
    dgemm_("NoTranspose", "NoTranspose", n, n, &np2, &c_b74, &bk[bk_offset], 
	    ldbk, &dwork[iwc + m2], &i__1, &c_b24, &ak[ak_offset], ldak, (
	    ftnlen)11, (ftnlen)11);

/*     Compute BKHAT in BK . */

#line 957 "SB10DD.f"
    dgemm_("NoTranspose", "NoTranspose", n, &np2, &m2, &c_b24, &b[(m1 + 1) * 
	    b_dim1 + 1], ldb, &dk[dk_offset], lddk, &c_b24, &bk[bk_offset], 
	    ldbk, (ftnlen)11, (ftnlen)11);

/*     Compute Im2 + DKHAT*D22 . */

#line 962 "SB10DD.f"
    iwrk = m2 * m2 + 1;
#line 963 "SB10DD.f"
    dlaset_("Full", &m2, &m2, &c_b25, &c_b24, &dwork[1], &m2, (ftnlen)4);
#line 964 "SB10DD.f"
    dgemm_("NoTranspose", "NoTranspose", &m2, &m2, &np2, &c_b24, &dk[
	    dk_offset], lddk, &d__[np1 + 1 + (m1 + 1) * d_dim1], ldd, &c_b24, 
	    &dwork[1], &m2, (ftnlen)11, (ftnlen)11);
#line 966 "SB10DD.f"
    anorm = dlange_("1", &m2, &m2, &dwork[1], &m2, &dwork[iwrk], (ftnlen)1);
#line 967 "SB10DD.f"
    dgetrf_(&m2, &m2, &dwork[1], &m2, &iwork[1], &info2);
#line 968 "SB10DD.f"
    if (info2 > 0) {
#line 969 "SB10DD.f"
	*info = 8;
#line 970 "SB10DD.f"
	return 0;
#line 971 "SB10DD.f"
    }
#line 972 "SB10DD.f"
    dgecon_("1", &m2, &dwork[1], &m2, &anorm, &rcond[6], &dwork[iwrk], &iwork[
	    m2 + 1], &info2, (ftnlen)1);

/*     Return if the matrix is singular to working precision. */

#line 977 "SB10DD.f"
    if (rcond[6] < toll) {
#line 978 "SB10DD.f"
	*info = 8;
#line 979 "SB10DD.f"
	return 0;
#line 980 "SB10DD.f"
    }

/*     Compute CK . */

#line 984 "SB10DD.f"
    dgetrs_("NoTranspose", &m2, n, &dwork[1], &m2, &iwork[1], &ck[ck_offset], 
	    ldck, &info2, (ftnlen)11);

/*     Compute DK . */

#line 989 "SB10DD.f"
    dgetrs_("NoTranspose", &m2, &np2, &dwork[1], &m2, &iwork[1], &dk[
	    dk_offset], lddk, &info2, (ftnlen)11);

/*     Compute AK . */

#line 994 "SB10DD.f"
    dgemm_("NoTranspose", "NoTranspose", n, &m2, &np2, &c_b24, &bk[bk_offset],
	     ldbk, &d__[np1 + 1 + (m1 + 1) * d_dim1], ldd, &c_b25, &dwork[1], 
	    n, (ftnlen)11, (ftnlen)11);
#line 996 "SB10DD.f"
    dgemm_("NoTranspose", "NoTranspose", n, n, &m2, &c_b74, &dwork[1], n, &ck[
	    ck_offset], ldck, &c_b24, &ak[ak_offset], ldak, (ftnlen)11, (
	    ftnlen)11);

/*     Compute BK . */

#line 1001 "SB10DD.f"
    dgemm_("NoTranspose", "NoTranspose", n, &np2, &m2, &c_b74, &dwork[1], n, &
	    dk[dk_offset], lddk, &c_b24, &bk[bk_offset], ldbk, (ftnlen)11, (
	    ftnlen)11);

#line 1004 "SB10DD.f"
    dwork[1] = (doublereal) lwamax;
#line 1005 "SB10DD.f"
    return 0;
/* *** Last line of SB10DD *** */
} /* sb10dd_ */

