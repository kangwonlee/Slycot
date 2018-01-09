#line 1 "SB10FD.f"
/* SB10FD.f -- translated by f2c (version 20100827).
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

#line 1 "SB10FD.f"
/* Subroutine */ int sb10fd_(integer *n, integer *m, integer *np, integer *
	ncon, integer *nmeas, doublereal *gamma, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *ak, integer *ldak, 
	doublereal *bk, integer *ldbk, doublereal *ck, integer *ldck, 
	doublereal *dk, integer *lddk, doublereal *rcond, doublereal *tol, 
	integer *iwork, doublereal *dwork, integer *ldwork, logical *bwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ak_dim1, ak_offset, b_dim1, b_offset, bk_dim1, 
	    bk_offset, c_dim1, c_offset, ck_dim1, ck_offset, d_dim1, d_offset,
	     dk_dim1, dk_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, 
	    i__8, i__9, i__10, i__11, i__12;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer m1, m2, nd1, nd2, np1, np2, lw1, lw2, lw3, lw4, lw5, lw6, 
	    iwc, iwd, iwf, iwh, iwx, iwy;
    static doublereal toll;
    static integer iwrk, iwtu, iwty, info2;
    extern /* Subroutine */ int sb10pd_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *), sb10qd_(
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    logical *, integer *), sb10rd_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static integer lwamax, minwrk;


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

/*              | AK | BK | */
/*          K = |----|----|, */
/*              | CK | DK | */

/*     using modified Glover's and Doyle's 1988 formulas, for the system */

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

/*     (A3) | A-j*omega*I  B2  | has full column rank for all omega, */
/*          |    C1        D12 | */

/*     (A4) | A-j*omega*I  B1  |  has full row rank for all omega. */
/*          |    C2        D21 | */

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
/*             GAMMA >= 0. */

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

/*     RCOND   (output) DOUBLE PRECISION array, dimension (4) */
/*             RCOND(1) contains the reciprocal condition number of the */
/*                      control transformation matrix; */
/*             RCOND(2) contains the reciprocal condition number of the */
/*                      measurement transformation matrix; */
/*             RCOND(3) contains an estimate of the reciprocal condition */
/*                      number of the X-Riccati equation; */
/*             RCOND(4) contains an estimate of the reciprocal condition */
/*                      number of the Y-Riccati equation. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             Tolerance used for controlling the accuracy of the applied */
/*             transformations for computing the normalized form in */
/*             SLICOT Library routine SB10PD. Transformation matrices */
/*             whose reciprocal condition numbers are less than TOL are */
/*             not allowed. If TOL <= 0, then a default value equal to */
/*             sqrt(EPS) is used, where EPS is the relative machine */
/*             precision. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK), where */
/*             LIWORK = max(2*max(N,M-NCON,NP-NMEAS,NCON),N*N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal */
/*             LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= N*M + NP*(N+M) + M2*M2 + NP2*NP2 + */
/*                       max(1,LW1,LW2,LW3,LW4,LW5,LW6), where */
/*             LW1 = (N+NP1+1)*(N+M2) + max(3*(N+M2)+N+NP1,5*(N+M2)), */
/*             LW2 = (N+NP2)*(N+M1+1) + max(3*(N+NP2)+N+M1,5*(N+NP2)), */
/*             LW3 = M2 + NP1*NP1 + max(NP1*max(N,M1),3*M2+NP1,5*M2), */
/*             LW4 = NP2 + M1*M1 + max(max(N,NP1)*M1,3*NP2+M1,5*NP2), */
/*             LW5 = 2*N*N + N*(M+NP) + */
/*                   max(1,M*M + max(2*M1,3*N*N+max(N*M,10*N*N+12*N+5)), */
/*                       NP*NP + max(2*NP1,3*N*N + */
/*                                   max(N*NP,10*N*N+12*N+5))), */
/*             LW6 = 2*N*N + N*(M+NP) + */
/*                   max(1, M2*NP2 + NP2*NP2 + M2*M2 + */
/*                       max(D1*D1 + max(2*D1, (D1+D2)*NP2), */
/*                           D2*D2 + max(2*D2, D2*M2), 3*N, */
/*                           N*(2*NP2 + M2) + */
/*                           max(2*N*M2, M2*NP2 + */
/*                                       max(M2*M2+3*M2, NP2*(2*NP2+ */
/*                                              M2+max(NP2,N)))))), */
/*             with D1 = NP1 - M2, D2 = M1 - NP2, */
/*                 NP1 = NP - NP2, M1 = M - M2. */
/*             For good performance, LDWORK must generally be larger. */
/*             Denoting Q = max(M1,M2,NP1,NP2), an upper bound is */
/*             2*Q*(3*Q+2*N)+max(1,(N+Q)*(N+Q+6),Q*(Q+max(N,Q,5)+1), */
/*               2*N*(N+2*Q)+max(1,4*Q*Q+ */
/*                               max(2*Q,3*N*N+max(2*N*Q,10*N*N+12*N+5)), */
/*                                 Q*(3*N+3*Q+max(2*N,4*Q+max(N,Q))))). */

/*     BWORK   LOGICAL array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the matrix | A-j*omega*I  B2  | had not full */
/*                                 |    C1        D12 | */
/*                   column rank in respect to the tolerance EPS; */
/*             = 2:  if the matrix | A-j*omega*I  B1  |  had not full row */
/*                                 |    C2        D21 | */
/*                   rank in respect to the tolerance EPS; */
/*             = 3:  if the matrix D12 had not full column rank in */
/*                   respect to the tolerance TOL; */
/*             = 4:  if the matrix D21 had not full row rank in respect */
/*                   to the tolerance TOL; */
/*             = 5:  if the singular value decomposition (SVD) algorithm */
/*                   did not converge (when computing the SVD of one of */
/*                   the matrices |A   B2 |, |A   B1 |, D12 or D21). */
/*                                |C1  D12|  |C2  D21| */
/*             = 6:  if the controller is not admissible (too small value */
/*                   of gamma); */
/*             = 7:  if the X-Riccati equation was not solved */
/*                   successfully (the controller is not admissible or */
/*                   there are numerical difficulties); */
/*             = 8:  if the Y-Riccati equation was not solved */
/*                   successfully (the controller is not admissible or */
/*                   there are numerical difficulties); */
/*             = 9:  if the determinant of Im2 + Tu*D11HAT*Ty*D22 is */
/*                   zero [3]. */

/*     METHOD */

/*     The routine implements the Glover's and Doyle's 1988 formulas [1], */
/*     [2] modified to improve the efficiency as described in [3]. */

/*     REFERENCES */

/*     [1] Glover, K. and Doyle, J.C. */
/*         State-space formulae for all stabilizing controllers that */
/*         satisfy an Hinf norm bound and relations to risk sensitivity. */
/*         Systems and Control Letters, vol. 11, pp. 167-172, 1988. */

/*     [2] Balas, G.J., Doyle, J.C., Glover, K., Packard, A., and */
/*         Smith, R. */
/*         mu-Analysis and Synthesis Toolbox. */
/*         The MathWorks Inc., Natick, Mass., 1995. */

/*     [3] Petkov, P.Hr., Gu, D.W., and Konstantinov, M.M. */
/*         Fortran 77 routines for Hinf and H2 design of continuous-time */
/*         linear control systems. */
/*         Rep. 98-14, Department of Engineering, Leicester University, */
/*         Leicester, U.K., 1998. */

/*     NUMERICAL ASPECTS */

/*     The accuracy of the result depends on the condition numbers of the */
/*     input and output transformations and on the condition numbers of */
/*     the two Riccati equations, as given by the values of RCOND(1), */
/*     RCOND(2), RCOND(3) and RCOND(4), respectively. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 1998. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 1999, */
/*     Sept. 1999, Feb. 2000. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, H-infinity optimal control, robust */
/*     control. */

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

#line 315 "SB10FD.f"
    /* Parameter adjustments */
#line 315 "SB10FD.f"
    a_dim1 = *lda;
#line 315 "SB10FD.f"
    a_offset = 1 + a_dim1;
#line 315 "SB10FD.f"
    a -= a_offset;
#line 315 "SB10FD.f"
    b_dim1 = *ldb;
#line 315 "SB10FD.f"
    b_offset = 1 + b_dim1;
#line 315 "SB10FD.f"
    b -= b_offset;
#line 315 "SB10FD.f"
    c_dim1 = *ldc;
#line 315 "SB10FD.f"
    c_offset = 1 + c_dim1;
#line 315 "SB10FD.f"
    c__ -= c_offset;
#line 315 "SB10FD.f"
    d_dim1 = *ldd;
#line 315 "SB10FD.f"
    d_offset = 1 + d_dim1;
#line 315 "SB10FD.f"
    d__ -= d_offset;
#line 315 "SB10FD.f"
    ak_dim1 = *ldak;
#line 315 "SB10FD.f"
    ak_offset = 1 + ak_dim1;
#line 315 "SB10FD.f"
    ak -= ak_offset;
#line 315 "SB10FD.f"
    bk_dim1 = *ldbk;
#line 315 "SB10FD.f"
    bk_offset = 1 + bk_dim1;
#line 315 "SB10FD.f"
    bk -= bk_offset;
#line 315 "SB10FD.f"
    ck_dim1 = *ldck;
#line 315 "SB10FD.f"
    ck_offset = 1 + ck_dim1;
#line 315 "SB10FD.f"
    ck -= ck_offset;
#line 315 "SB10FD.f"
    dk_dim1 = *lddk;
#line 315 "SB10FD.f"
    dk_offset = 1 + dk_dim1;
#line 315 "SB10FD.f"
    dk -= dk_offset;
#line 315 "SB10FD.f"
    --rcond;
#line 315 "SB10FD.f"
    --iwork;
#line 315 "SB10FD.f"
    --dwork;
#line 315 "SB10FD.f"
    --bwork;
#line 315 "SB10FD.f"

#line 315 "SB10FD.f"
    /* Function Body */
#line 315 "SB10FD.f"
    m1 = *m - *ncon;
#line 316 "SB10FD.f"
    m2 = *ncon;
#line 317 "SB10FD.f"
    np1 = *np - *nmeas;
#line 318 "SB10FD.f"
    np2 = *nmeas;

#line 320 "SB10FD.f"
    *info = 0;
#line 321 "SB10FD.f"
    if (*n < 0) {
#line 322 "SB10FD.f"
	*info = -1;
#line 323 "SB10FD.f"
    } else if (*m < 0) {
#line 324 "SB10FD.f"
	*info = -2;
#line 325 "SB10FD.f"
    } else if (*np < 0) {
#line 326 "SB10FD.f"
	*info = -3;
#line 327 "SB10FD.f"
    } else if (*ncon < 0 || m1 < 0 || m2 > np1) {
#line 328 "SB10FD.f"
	*info = -4;
#line 329 "SB10FD.f"
    } else if (*nmeas < 0 || np1 < 0 || np2 > m1) {
#line 330 "SB10FD.f"
	*info = -5;
#line 331 "SB10FD.f"
    } else if (*gamma < 0.) {
#line 332 "SB10FD.f"
	*info = -6;
#line 333 "SB10FD.f"
    } else if (*lda < max(1,*n)) {
#line 334 "SB10FD.f"
	*info = -8;
#line 335 "SB10FD.f"
    } else if (*ldb < max(1,*n)) {
#line 336 "SB10FD.f"
	*info = -10;
#line 337 "SB10FD.f"
    } else if (*ldc < max(1,*np)) {
#line 338 "SB10FD.f"
	*info = -12;
#line 339 "SB10FD.f"
    } else if (*ldd < max(1,*np)) {
#line 340 "SB10FD.f"
	*info = -14;
#line 341 "SB10FD.f"
    } else if (*ldak < max(1,*n)) {
#line 342 "SB10FD.f"
	*info = -16;
#line 343 "SB10FD.f"
    } else if (*ldbk < max(1,*n)) {
#line 344 "SB10FD.f"
	*info = -18;
#line 345 "SB10FD.f"
    } else if (*ldck < max(1,m2)) {
#line 346 "SB10FD.f"
	*info = -20;
#line 347 "SB10FD.f"
    } else if (*lddk < max(1,m2)) {
#line 348 "SB10FD.f"
	*info = -22;
#line 349 "SB10FD.f"
    } else {

/*        Compute workspace. */

#line 353 "SB10FD.f"
	nd1 = np1 - m2;
#line 354 "SB10FD.f"
	nd2 = m1 - np2;
/* Computing MAX */
#line 355 "SB10FD.f"
	i__1 = (*n + m2) * 3 + *n + np1, i__2 = (*n + m2) * 5;
#line 355 "SB10FD.f"
	lw1 = (*n + np1 + 1) * (*n + m2) + max(i__1,i__2);
/* Computing MAX */
#line 357 "SB10FD.f"
	i__1 = (*n + np2) * 3 + *n + m1, i__2 = (*n + np2) * 5;
#line 357 "SB10FD.f"
	lw2 = (*n + np2) * (*n + m1 + 1) + max(i__1,i__2);
/* Computing MAX */
#line 359 "SB10FD.f"
	i__1 = np1 * max(*n,m1), i__2 = m2 * 3 + np1, i__1 = max(i__1,i__2), 
		i__2 = m2 * 5;
#line 359 "SB10FD.f"
	lw3 = m2 + np1 * np1 + max(i__1,i__2);
/* Computing MAX */
#line 360 "SB10FD.f"
	i__1 = max(*n,np1) * m1, i__2 = np2 * 3 + m1, i__1 = max(i__1,i__2), 
		i__2 = np2 * 5;
#line 360 "SB10FD.f"
	lw4 = np2 + m1 * m1 + max(i__1,i__2);
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
#line 361 "SB10FD.f"
	i__5 = *n * *m, i__6 = *n * 10 * *n + *n * 12 + 5;
#line 361 "SB10FD.f"
	i__3 = m1 << 1, i__4 = *n * 3 * *n + max(i__5,i__6);
/* Computing MAX */
/* Computing MAX */
#line 361 "SB10FD.f"
	i__9 = *n * *np, i__10 = *n * 10 * *n + *n * 12 + 5;
#line 361 "SB10FD.f"
	i__7 = np1 << 1, i__8 = *n * 3 * *n + max(i__9,i__10);
#line 361 "SB10FD.f"
	i__1 = 1, i__2 = *m * *m + max(i__3,i__4), i__1 = max(i__1,i__2), 
		i__2 = *np * *np + max(i__7,i__8);
#line 361 "SB10FD.f"
	lw5 = (*n << 1) * *n + *n * (*m + *np) + max(i__1,i__2);
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
#line 366 "SB10FD.f"
	i__5 = nd1 << 1, i__6 = (nd1 + nd2) * np2;
/* Computing MAX */
#line 366 "SB10FD.f"
	i__7 = nd2 << 1, i__8 = nd2 * m2;
/* Computing MAX */
/* Computing MAX */
#line 366 "SB10FD.f"
	i__11 = m2 * m2 + m2 * 3, i__12 = np2 * ((np2 << 1) + m2 + max(np2,*n)
		);
#line 366 "SB10FD.f"
	i__9 = (*n << 1) * m2, i__10 = m2 * np2 + max(i__11,i__12);
#line 366 "SB10FD.f"
	i__3 = nd1 * nd1 + max(i__5,i__6), i__4 = nd2 * nd2 + max(i__7,i__8), 
		i__3 = max(i__3,i__4), i__4 = *n * 3, i__3 = max(i__3,i__4), 
		i__4 = *n * ((np2 << 1) + m2) + max(i__9,i__10);
#line 366 "SB10FD.f"
	i__1 = 1, i__2 = m2 * np2 + np2 * np2 + m2 * m2 + max(i__3,i__4);
#line 366 "SB10FD.f"
	lw6 = (*n << 1) * *n + *n * (*m + *np) + max(i__1,i__2);
/* Computing MAX */
#line 374 "SB10FD.f"
	i__1 = max(1,lw1), i__1 = max(i__1,lw2), i__1 = max(i__1,lw3), i__1 = 
		max(i__1,lw4), i__1 = max(i__1,lw5);
#line 374 "SB10FD.f"
	minwrk = *n * *m + *np * (*n + *m) + m2 * m2 + np2 * np2 + max(i__1,
		lw6);
#line 376 "SB10FD.f"
	if (*ldwork < minwrk) {
#line 376 "SB10FD.f"
	    *info = -27;
#line 376 "SB10FD.f"
	}
#line 378 "SB10FD.f"
    }
#line 379 "SB10FD.f"
    if (*info != 0) {
#line 380 "SB10FD.f"
	i__1 = -(*info);
#line 380 "SB10FD.f"
	xerbla_("SB10FD", &i__1, (ftnlen)6);
#line 381 "SB10FD.f"
	return 0;
#line 382 "SB10FD.f"
    }

/*     Quick return if possible. */

#line 386 "SB10FD.f"
    if (*n == 0 || *m == 0 || *np == 0 || m1 == 0 || m2 == 0 || np1 == 0 || 
	    np2 == 0) {
#line 388 "SB10FD.f"
	rcond[1] = 1.;
#line 389 "SB10FD.f"
	rcond[2] = 1.;
#line 390 "SB10FD.f"
	rcond[3] = 1.;
#line 391 "SB10FD.f"
	rcond[4] = 1.;
#line 392 "SB10FD.f"
	dwork[1] = 1.;
#line 393 "SB10FD.f"
	return 0;
#line 394 "SB10FD.f"
    }

#line 396 "SB10FD.f"
    toll = *tol;
#line 397 "SB10FD.f"
    if (toll <= 0.) {

/*        Set the default value of the tolerance. */

#line 401 "SB10FD.f"
	toll = sqrt(dlamch_("Epsilon", (ftnlen)7));
#line 402 "SB10FD.f"
    }

/*     Workspace usage. */

#line 406 "SB10FD.f"
    iwc = *n * *m + 1;
#line 407 "SB10FD.f"
    iwd = iwc + *np * *n;
#line 408 "SB10FD.f"
    iwtu = iwd + *np * *m;
#line 409 "SB10FD.f"
    iwty = iwtu + m2 * m2;
#line 410 "SB10FD.f"
    iwrk = iwty + np2 * np2;

#line 412 "SB10FD.f"
    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[1], n, (ftnlen)4);
#line 413 "SB10FD.f"
    dlacpy_("Full", np, n, &c__[c_offset], ldc, &dwork[iwc], np, (ftnlen)4);
#line 414 "SB10FD.f"
    dlacpy_("Full", np, m, &d__[d_offset], ldd, &dwork[iwd], np, (ftnlen)4);

/*     Transform the system so that D12 and D21 satisfy the formulas */
/*     in the computation of the Hinf (sub)optimal controller. */

#line 419 "SB10FD.f"
    i__1 = *ldwork - iwrk + 1;
#line 419 "SB10FD.f"
    sb10pd_(n, m, np, ncon, nmeas, &a[a_offset], lda, &dwork[1], n, &dwork[
	    iwc], np, &dwork[iwd], np, &dwork[iwtu], &m2, &dwork[iwty], &np2, 
	    &rcond[1], &toll, &dwork[iwrk], &i__1, &info2);
#line 423 "SB10FD.f"
    if (info2 > 0) {
#line 424 "SB10FD.f"
	*info = info2;
#line 425 "SB10FD.f"
	return 0;
#line 426 "SB10FD.f"
    }
#line 427 "SB10FD.f"
    lwamax = (integer) dwork[iwrk] + iwrk - 1;

#line 429 "SB10FD.f"
    iwx = iwrk;
#line 430 "SB10FD.f"
    iwy = iwx + *n * *n;
#line 431 "SB10FD.f"
    iwf = iwy + *n * *n;
#line 432 "SB10FD.f"
    iwh = iwf + *m * *n;
#line 433 "SB10FD.f"
    iwrk = iwh + *n * *np;

/*     Compute the (sub)optimal state feedback and output injection */
/*     matrices. */

#line 438 "SB10FD.f"
    i__1 = *ldwork - iwrk + 1;
#line 438 "SB10FD.f"
    sb10qd_(n, m, np, ncon, nmeas, gamma, &a[a_offset], lda, &dwork[1], n, &
	    dwork[iwc], np, &dwork[iwd], np, &dwork[iwf], m, &dwork[iwh], n, &
	    dwork[iwx], n, &dwork[iwy], n, &rcond[3], &iwork[1], &dwork[iwrk],
	     &i__1, &bwork[1], &info2);
#line 443 "SB10FD.f"
    if (info2 > 0) {
#line 444 "SB10FD.f"
	*info = info2 + 5;
#line 445 "SB10FD.f"
	return 0;
#line 446 "SB10FD.f"
    }
/* Computing MAX */
#line 447 "SB10FD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 447 "SB10FD.f"
    lwamax = max(i__1,lwamax);

/*     Compute the Hinf (sub)optimal controller. */

#line 451 "SB10FD.f"
    i__1 = *ldwork - iwrk + 1;
#line 451 "SB10FD.f"
    sb10rd_(n, m, np, ncon, nmeas, gamma, &a[a_offset], lda, &dwork[1], n, &
	    dwork[iwc], np, &dwork[iwd], np, &dwork[iwf], m, &dwork[iwh], n, &
	    dwork[iwtu], &m2, &dwork[iwty], &np2, &dwork[iwx], n, &dwork[iwy],
	     n, &ak[ak_offset], ldak, &bk[bk_offset], ldbk, &ck[ck_offset], 
	    ldck, &dk[dk_offset], lddk, &iwork[1], &dwork[iwrk], &i__1, &
	    info2);
#line 457 "SB10FD.f"
    if (info2 == 1) {
#line 458 "SB10FD.f"
	*info = 6;
#line 459 "SB10FD.f"
	return 0;
#line 460 "SB10FD.f"
    } else if (info2 == 2) {
#line 461 "SB10FD.f"
	*info = 9;
#line 462 "SB10FD.f"
	return 0;
#line 463 "SB10FD.f"
    }
/* Computing MAX */
#line 464 "SB10FD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 464 "SB10FD.f"
    lwamax = max(i__1,lwamax);

#line 466 "SB10FD.f"
    dwork[1] = (doublereal) lwamax;
#line 467 "SB10FD.f"
    return 0;
/* *** Last line of SB10FD *** */
} /* sb10fd_ */
