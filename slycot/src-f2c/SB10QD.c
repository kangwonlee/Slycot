#line 1 "SB10QD.f"
/* SB10QD.f -- translated by f2c (version 20100827).
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

#line 1 "SB10QD.f"
/* Table of constant values */

static doublereal c_b5 = 0.;
static doublereal c_b8 = 1.;
static doublereal c_b17 = -1.;

/* Subroutine */ int sb10qd_(integer *n, integer *m, integer *np, integer *
	ncon, integer *nmeas, doublereal *gamma, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *f, integer *ldf, 
	doublereal *h__, integer *ldh, doublereal *x, integer *ldx, 
	doublereal *y, integer *ldy, doublereal *xycond, integer *iwork, 
	doublereal *dwork, integer *ldwork, logical *bwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, f_dim1, f_offset, h_dim1, h_offset, x_dim1, x_offset, 
	    y_dim1, y_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, 
	    i__9, i__10;
    doublereal d__1;

    /* Local variables */
    static integer m1, m2, n2, nn, nd1, nd2, np1, np2, iw2, iwa, iwg, iwi;
    static doublereal eps, sep;
    static integer iwq, iwr, iws, iwt, iwv;
    static doublereal ferr;
    static integer iwrk, info2;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     sb02rd_(char *, char *, char *, char *, char *, char *, char *, 
	    char *, char *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    logical *, integer *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, 
	    ftnlen, ftnlen, ftnlen, ftnlen);
    static doublereal rcond;
    extern /* Subroutine */ int mb01ru_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), mb01rx_(char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int dsymm_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen), dsyrk_(
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static integer lwamax;
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int dsycon_(char *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer minwrk;
    extern /* Subroutine */ int dsytrf_(char *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     dsytri_(char *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, ftnlen);


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

/*     To compute the state feedback and the output injection */
/*     matrices for an H-infinity (sub)optimal n-state controller, */
/*     using Glover's and Doyle's 1988 formulas, for the system */

/*                   | A  | B1  B2  |   | A | B | */
/*               P = |----|---------| = |---|---| */
/*                   | C1 | D11 D12 |   | C | D | */
/*                   | C2 | D21 D22 | */

/*     and for a given value of gamma, where B2 has as column size the */
/*     number of control inputs (NCON) and C2 has as row size the number */
/*     of measurements (NMEAS) being provided to the controller. */

/*     It is assumed that */

/*     (A1) (A,B2) is stabilizable and (C2,A) is detectable, */

/*     (A2) D12 is full column rank with D12 = | 0 | and D21 is */
/*                                             | I | */
/*          full row rank with D21 = | 0 I | as obtained by the */
/*          subroutine SB10PD, */

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

/*     F       (output) DOUBLE PRECISION array, dimension (LDF,N) */
/*             The leading M-by-N part of this array contains the state */
/*             feedback matrix F. */

/*     LDF     INTEGER */
/*             The leading dimension of the array F.  LDF >= max(1,M). */

/*     H       (output) DOUBLE PRECISION array, dimension (LDH,NP) */
/*             The leading N-by-NP part of this array contains the output */
/*             injection matrix H. */

/*     LDH     INTEGER */
/*             The leading dimension of the array H.  LDH >= max(1,N). */

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

/*     XYCOND  (output) DOUBLE PRECISION array, dimension (2) */
/*             XYCOND(1) contains an estimate of the reciprocal condition */
/*                       number of the X-Riccati equation; */
/*             XYCOND(2) contains an estimate of the reciprocal condition */
/*                       number of the Y-Riccati equation. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension max(2*max(N,M-NCON,NP-NMEAS),N*N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal */
/*             LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= max(1,M*M + max(2*M1,3*N*N + */
/*                                       max(N*M,10*N*N+12*N+5)), */
/*                           NP*NP + max(2*NP1,3*N*N + */
/*                                       max(N*NP,10*N*N+12*N+5))), */
/*             where M1 = M - M2 and NP1 = NP - NP2. */
/*             For good performance, LDWORK must generally be larger. */
/*             Denoting Q = MAX(M1,M2,NP1,NP2), an upper bound is */
/*             max(1,4*Q*Q+max(2*Q,3*N*N + max(2*N*Q,10*N*N+12*N+5))). */

/*     BWORK   LOGICAL array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the controller is not admissible (too small value */
/*                   of gamma); */
/*             = 2:  if the X-Riccati equation was not solved */
/*                   successfully (the controller is not admissible or */
/*                   there are numerical difficulties); */
/*             = 3:  if the Y-Riccati equation was not solved */
/*                   successfully (the controller is not admissible or */
/*                   there are numerical difficulties). */

/*     METHOD */

/*     The routine implements the Glover's and Doyle's formulas [1],[2] */
/*     modified as described in [3]. The X- and Y-Riccati equations */
/*     are solved with condition and accuracy estimates [4]. */

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

/*     [4] Petkov, P.Hr., Konstantinov, M.M., and Mehrmann, V. */
/*         DGRSVX and DMSRIC: Fortan 77 subroutines for solving */
/*         continuous-time matrix algebraic Riccati equations with */
/*         condition and accuracy estimates. */
/*         Preprint SFB393/98-16, Fak. f. Mathematik, Tech. Univ. */
/*         Chemnitz, May 1998. */

/*     NUMERICAL ASPECTS */

/*     The precision of the solution of the matrix Riccati equations */
/*     can be controlled by the values of the condition numbers */
/*     XYCOND(1) and XYCOND(2) of these equations. */

/*     FURTHER COMMENTS */

/*     The Riccati equations are solved by the Schur approach */
/*     implementing condition and accuracy estimates. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 1998. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 1999, */
/*     Sept. 1999. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, H-infinity optimal control, robust */
/*     control. */

/*     ****************************************************************** */

/*     .. Parameters .. */

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

#line 277 "SB10QD.f"
    /* Parameter adjustments */
#line 277 "SB10QD.f"
    a_dim1 = *lda;
#line 277 "SB10QD.f"
    a_offset = 1 + a_dim1;
#line 277 "SB10QD.f"
    a -= a_offset;
#line 277 "SB10QD.f"
    b_dim1 = *ldb;
#line 277 "SB10QD.f"
    b_offset = 1 + b_dim1;
#line 277 "SB10QD.f"
    b -= b_offset;
#line 277 "SB10QD.f"
    c_dim1 = *ldc;
#line 277 "SB10QD.f"
    c_offset = 1 + c_dim1;
#line 277 "SB10QD.f"
    c__ -= c_offset;
#line 277 "SB10QD.f"
    d_dim1 = *ldd;
#line 277 "SB10QD.f"
    d_offset = 1 + d_dim1;
#line 277 "SB10QD.f"
    d__ -= d_offset;
#line 277 "SB10QD.f"
    f_dim1 = *ldf;
#line 277 "SB10QD.f"
    f_offset = 1 + f_dim1;
#line 277 "SB10QD.f"
    f -= f_offset;
#line 277 "SB10QD.f"
    h_dim1 = *ldh;
#line 277 "SB10QD.f"
    h_offset = 1 + h_dim1;
#line 277 "SB10QD.f"
    h__ -= h_offset;
#line 277 "SB10QD.f"
    x_dim1 = *ldx;
#line 277 "SB10QD.f"
    x_offset = 1 + x_dim1;
#line 277 "SB10QD.f"
    x -= x_offset;
#line 277 "SB10QD.f"
    y_dim1 = *ldy;
#line 277 "SB10QD.f"
    y_offset = 1 + y_dim1;
#line 277 "SB10QD.f"
    y -= y_offset;
#line 277 "SB10QD.f"
    --xycond;
#line 277 "SB10QD.f"
    --iwork;
#line 277 "SB10QD.f"
    --dwork;
#line 277 "SB10QD.f"
    --bwork;
#line 277 "SB10QD.f"

#line 277 "SB10QD.f"
    /* Function Body */
#line 277 "SB10QD.f"
    m1 = *m - *ncon;
#line 278 "SB10QD.f"
    m2 = *ncon;
#line 279 "SB10QD.f"
    np1 = *np - *nmeas;
#line 280 "SB10QD.f"
    np2 = *nmeas;
#line 281 "SB10QD.f"
    nn = *n * *n;

#line 283 "SB10QD.f"
    *info = 0;
#line 284 "SB10QD.f"
    if (*n < 0) {
#line 285 "SB10QD.f"
	*info = -1;
#line 286 "SB10QD.f"
    } else if (*m < 0) {
#line 287 "SB10QD.f"
	*info = -2;
#line 288 "SB10QD.f"
    } else if (*np < 0) {
#line 289 "SB10QD.f"
	*info = -3;
#line 290 "SB10QD.f"
    } else if (*ncon < 0 || m1 < 0 || m2 > np1) {
#line 291 "SB10QD.f"
	*info = -4;
#line 292 "SB10QD.f"
    } else if (*nmeas < 0 || np1 < 0 || np2 > m1) {
#line 293 "SB10QD.f"
	*info = -5;
#line 294 "SB10QD.f"
    } else if (*gamma < 0.) {
#line 295 "SB10QD.f"
	*info = -6;
#line 296 "SB10QD.f"
    } else if (*lda < max(1,*n)) {
#line 297 "SB10QD.f"
	*info = -8;
#line 298 "SB10QD.f"
    } else if (*ldb < max(1,*n)) {
#line 299 "SB10QD.f"
	*info = -10;
#line 300 "SB10QD.f"
    } else if (*ldc < max(1,*np)) {
#line 301 "SB10QD.f"
	*info = -12;
#line 302 "SB10QD.f"
    } else if (*ldd < max(1,*np)) {
#line 303 "SB10QD.f"
	*info = -14;
#line 304 "SB10QD.f"
    } else if (*ldf < max(1,*m)) {
#line 305 "SB10QD.f"
	*info = -16;
#line 306 "SB10QD.f"
    } else if (*ldh < max(1,*n)) {
#line 307 "SB10QD.f"
	*info = -18;
#line 308 "SB10QD.f"
    } else if (*ldx < max(1,*n)) {
#line 309 "SB10QD.f"
	*info = -20;
#line 310 "SB10QD.f"
    } else if (*ldy < max(1,*n)) {
#line 311 "SB10QD.f"
	*info = -22;
#line 312 "SB10QD.f"
    } else {

/*        Compute workspace. */

/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
#line 316 "SB10QD.f"
	i__5 = *n * *m, i__6 = nn * 10 + *n * 12 + 5;
#line 316 "SB10QD.f"
	i__3 = m1 << 1, i__4 = nn * 3 + max(i__5,i__6);
/* Computing MAX */
/* Computing MAX */
#line 316 "SB10QD.f"
	i__9 = *n * *np, i__10 = nn * 10 + *n * 12 + 5;
#line 316 "SB10QD.f"
	i__7 = np1 << 1, i__8 = nn * 3 + max(i__9,i__10);
#line 316 "SB10QD.f"
	i__1 = 1, i__2 = *m * *m + max(i__3,i__4), i__1 = max(i__1,i__2), 
		i__2 = *np * *np + max(i__7,i__8);
#line 316 "SB10QD.f"
	minwrk = max(i__1,i__2);
#line 320 "SB10QD.f"
	if (*ldwork < minwrk) {
#line 320 "SB10QD.f"
	    *info = -26;
#line 320 "SB10QD.f"
	}
#line 322 "SB10QD.f"
    }
#line 323 "SB10QD.f"
    if (*info != 0) {
#line 324 "SB10QD.f"
	i__1 = -(*info);
#line 324 "SB10QD.f"
	xerbla_("SB10QD", &i__1, (ftnlen)6);
#line 325 "SB10QD.f"
	return 0;
#line 326 "SB10QD.f"
    }

/*     Quick return if possible. */

#line 330 "SB10QD.f"
    if (*n == 0 || *m == 0 || *np == 0 || m1 == 0 || m2 == 0 || np1 == 0 || 
	    np2 == 0) {
#line 332 "SB10QD.f"
	xycond[1] = 1.;
#line 333 "SB10QD.f"
	xycond[2] = 1.;
#line 334 "SB10QD.f"
	dwork[1] = 1.;
#line 335 "SB10QD.f"
	return 0;
#line 336 "SB10QD.f"
    }
#line 337 "SB10QD.f"
    nd1 = np1 - m2;
#line 338 "SB10QD.f"
    nd2 = m1 - np2;
#line 339 "SB10QD.f"
    n2 = *n << 1;

/*     Get the machine precision. */

#line 343 "SB10QD.f"
    eps = dlamch_("Epsilon", (ftnlen)7);

/*     Workspace usage. */

#line 347 "SB10QD.f"
    iwa = *m * *m + 1;
#line 348 "SB10QD.f"
    iwq = iwa + nn;
#line 349 "SB10QD.f"
    iwg = iwq + nn;
#line 350 "SB10QD.f"
    iw2 = iwg + nn;

/*     Compute |D1111'||D1111 D1112| - gamma^2*Im1 . */
/*             |D1112'| */

#line 355 "SB10QD.f"
    d__1 = -(*gamma) * *gamma;
#line 355 "SB10QD.f"
    dlaset_("L", &m1, &m1, &c_b5, &d__1, &dwork[1], m, (ftnlen)1);
#line 356 "SB10QD.f"
    if (nd1 > 0) {
#line 356 "SB10QD.f"
	dsyrk_("L", "T", &m1, &nd1, &c_b8, &d__[d_offset], ldd, &c_b8, &dwork[
		1], m, (ftnlen)1, (ftnlen)1);
#line 356 "SB10QD.f"
    }

/*     Compute inv(|D1111'|*|D1111 D1112| - gamma^2*Im1) . */
/*                 |D1112'| */

#line 362 "SB10QD.f"
    iwrk = iwa;
#line 363 "SB10QD.f"
    anorm = dlansy_("I", "L", &m1, &dwork[1], m, &dwork[iwrk], (ftnlen)1, (
	    ftnlen)1);
#line 364 "SB10QD.f"
    i__1 = *ldwork - iwrk + 1;
#line 364 "SB10QD.f"
    dsytrf_("L", &m1, &dwork[1], m, &iwork[1], &dwork[iwrk], &i__1, &info2, (
	    ftnlen)1);
#line 366 "SB10QD.f"
    if (info2 > 0) {
#line 367 "SB10QD.f"
	*info = 1;
#line 368 "SB10QD.f"
	return 0;
#line 369 "SB10QD.f"
    }

#line 371 "SB10QD.f"
    lwamax = (integer) dwork[iwrk] + iwrk - 1;
#line 372 "SB10QD.f"
    dsycon_("L", &m1, &dwork[1], m, &iwork[1], &anorm, &rcond, &dwork[iwrk], &
	    iwork[m1 + 1], &info2, (ftnlen)1);
#line 374 "SB10QD.f"
    if (rcond < eps) {
#line 375 "SB10QD.f"
	*info = 1;
#line 376 "SB10QD.f"
	return 0;
#line 377 "SB10QD.f"
    }

/*     Compute inv(R) block by block. */

#line 381 "SB10QD.f"
    dsytri_("L", &m1, &dwork[1], m, &iwork[1], &dwork[iwrk], &info2, (ftnlen)
	    1);

/*     Compute -|D1121 D1122|*inv(|D1111'|*|D1111 D1112| - gamma^2*Im1) . */
/*                                |D1112'| */

#line 386 "SB10QD.f"
    dsymm_("R", "L", &m2, &m1, &c_b17, &dwork[1], m, &d__[nd1 + 1 + d_dim1], 
	    ldd, &c_b5, &dwork[m1 + 1], m, (ftnlen)1, (ftnlen)1);

/*     Compute |D1121 D1122|*inv(|D1111'|*|D1111 D1112| - */
/*                               |D1112'| */

/*                  gamma^2*Im1)*|D1121'| + Im2 . */
/*                               |D1122'| */

#line 395 "SB10QD.f"
    dlaset_("Lower", &m2, &m2, &c_b5, &c_b8, &dwork[m1 * (*m + 1) + 1], m, (
	    ftnlen)5);
#line 396 "SB10QD.f"
    mb01rx_("Right", "Lower", "Transpose", &m2, &m1, &c_b8, &c_b17, &dwork[m1 
	    * (*m + 1) + 1], m, &d__[nd1 + 1 + d_dim1], ldd, &dwork[m1 + 1], 
	    m, &info2, (ftnlen)5, (ftnlen)5, (ftnlen)9);

/*     Compute D11'*C1 . */

#line 402 "SB10QD.f"
    dgemm_("T", "N", &m1, n, &np1, &c_b8, &d__[d_offset], ldd, &c__[c_offset],
	     ldc, &c_b5, &dwork[iw2], m, (ftnlen)1, (ftnlen)1);

/*     Compute D1D'*C1 . */

#line 407 "SB10QD.f"
    dlacpy_("Full", &m2, n, &c__[nd1 + 1 + c_dim1], ldc, &dwork[iw2 + m1], m, 
	    (ftnlen)4);

/*     Compute inv(R)*D1D'*C1 in F . */

#line 412 "SB10QD.f"
    dsymm_("L", "L", m, n, &c_b8, &dwork[1], m, &dwork[iw2], m, &c_b5, &f[
	    f_offset], ldf, (ftnlen)1, (ftnlen)1);

/*     Compute Ax = A - B*inv(R)*D1D'*C1 . */

#line 417 "SB10QD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[iwa], n, (ftnlen)4);
#line 418 "SB10QD.f"
    dgemm_("N", "N", n, n, m, &c_b17, &b[b_offset], ldb, &f[f_offset], ldf, &
	    c_b8, &dwork[iwa], n, (ftnlen)1, (ftnlen)1);

/*     Compute Cx = C1'*C1 - C1'*D1D*inv(R)*D1D'*C1 . */

#line 423 "SB10QD.f"
    if (nd1 == 0) {
#line 424 "SB10QD.f"
	dlaset_("L", n, n, &c_b5, &c_b5, &dwork[iwq], n, (ftnlen)1);
#line 425 "SB10QD.f"
    } else {
#line 426 "SB10QD.f"
	dsyrk_("L", "T", n, &np1, &c_b8, &c__[c_offset], ldc, &c_b5, &dwork[
		iwq], n, (ftnlen)1, (ftnlen)1);
#line 428 "SB10QD.f"
	mb01rx_("Left", "Lower", "Transpose", n, m, &c_b8, &c_b17, &dwork[iwq]
		, n, &dwork[iw2], m, &f[f_offset], ldf, &info2, (ftnlen)4, (
		ftnlen)5, (ftnlen)9);
#line 430 "SB10QD.f"
    }

/*     Compute Dx = B*inv(R)*B' . */

#line 434 "SB10QD.f"
    iwrk = iw2;
#line 435 "SB10QD.f"
    i__1 = *m * *n;
#line 435 "SB10QD.f"
    mb01ru_("Lower", "NoTranspose", n, m, &c_b5, &c_b8, &dwork[iwg], n, &b[
	    b_offset], ldb, &dwork[1], m, &dwork[iwrk], &i__1, &info2, (
	    ftnlen)5, (ftnlen)11);

/*     Solution of the Riccati equation Ax'*X + X*Ax + Cx - X*Dx*X = 0 . */
/*     Workspace:  need   M*M + 13*N*N + 12*N + 5; */
/*                 prefer larger. */

#line 443 "SB10QD.f"
    iwt = iw2;
#line 444 "SB10QD.f"
    iwv = iwt + nn;
#line 445 "SB10QD.f"
    iwr = iwv + nn;
#line 446 "SB10QD.f"
    iwi = iwr + n2;
#line 447 "SB10QD.f"
    iws = iwi + n2;
#line 448 "SB10QD.f"
    iwrk = iws + (nn << 2);

#line 450 "SB10QD.f"
    i__1 = *ldwork - iwrk + 1;
#line 450 "SB10QD.f"
    sb02rd_("All", "Continuous", "NotUsed", "NoTranspose", "Lower", "General"\
	    "Scaling", "Stable", "NotFactored", "Original", n, &dwork[iwa], n, 
	    &dwork[iwt], n, &dwork[iwv], n, &dwork[iwg], n, &dwork[iwq], n, &
	    x[x_offset], ldx, &sep, &xycond[1], &ferr, &dwork[iwr], &dwork[
	    iwi], &dwork[iws], &n2, &iwork[1], &dwork[iwrk], &i__1, &bwork[1],
	     &info2, (ftnlen)3, (ftnlen)10, (ftnlen)7, (ftnlen)11, (ftnlen)5, 
	    (ftnlen)14, (ftnlen)6, (ftnlen)11, (ftnlen)8);
#line 457 "SB10QD.f"
    if (info2 > 0) {
#line 458 "SB10QD.f"
	*info = 2;
#line 459 "SB10QD.f"
	return 0;
#line 460 "SB10QD.f"
    }

/* Computing MAX */
#line 462 "SB10QD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 462 "SB10QD.f"
    lwamax = max(i__1,lwamax);

/*     Compute F = -inv(R)*|D1D'*C1 + B'*X| . */

#line 466 "SB10QD.f"
    iwrk = iw2;
#line 467 "SB10QD.f"
    dgemm_("T", "N", m, n, n, &c_b8, &b[b_offset], ldb, &x[x_offset], ldx, &
	    c_b5, &dwork[iwrk], m, (ftnlen)1, (ftnlen)1);
#line 469 "SB10QD.f"
    dsymm_("L", "L", m, n, &c_b17, &dwork[1], m, &dwork[iwrk], m, &c_b17, &f[
	    f_offset], ldf, (ftnlen)1, (ftnlen)1);

/*     Workspace usage. */

#line 474 "SB10QD.f"
    iwa = *np * *np + 1;
#line 475 "SB10QD.f"
    iwq = iwa + nn;
#line 476 "SB10QD.f"
    iwg = iwq + nn;
#line 477 "SB10QD.f"
    iw2 = iwg + nn;

/*     Compute |D1111|*|D1111' D1121'| - gamma^2*Inp1 . */
/*             |D1121| */

#line 482 "SB10QD.f"
    d__1 = -(*gamma) * *gamma;
#line 482 "SB10QD.f"
    dlaset_("U", &np1, &np1, &c_b5, &d__1, &dwork[1], np, (ftnlen)1);
#line 483 "SB10QD.f"
    if (nd2 > 0) {
#line 483 "SB10QD.f"
	dsyrk_("U", "N", &np1, &nd2, &c_b8, &d__[d_offset], ldd, &c_b8, &
		dwork[1], np, (ftnlen)1, (ftnlen)1);
#line 483 "SB10QD.f"
    }

/*     Compute inv(|D1111|*|D1111' D1121'| - gamma^2*Inp1) . */
/*                 |D1121| */

#line 489 "SB10QD.f"
    iwrk = iwa;
#line 490 "SB10QD.f"
    anorm = dlansy_("I", "U", &np1, &dwork[1], np, &dwork[iwrk], (ftnlen)1, (
	    ftnlen)1);
#line 491 "SB10QD.f"
    i__1 = *ldwork - iwrk + 1;
#line 491 "SB10QD.f"
    dsytrf_("U", &np1, &dwork[1], np, &iwork[1], &dwork[iwrk], &i__1, &info2, 
	    (ftnlen)1);
#line 493 "SB10QD.f"
    if (info2 > 0) {
#line 494 "SB10QD.f"
	*info = 1;
#line 495 "SB10QD.f"
	return 0;
#line 496 "SB10QD.f"
    }

/* Computing MAX */
#line 498 "SB10QD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 498 "SB10QD.f"
    lwamax = max(i__1,lwamax);
#line 499 "SB10QD.f"
    dsycon_("U", &np1, &dwork[1], np, &iwork[1], &anorm, &rcond, &dwork[iwrk],
	     &iwork[np1 + 1], &info2, (ftnlen)1);
#line 501 "SB10QD.f"
    if (rcond < eps) {
#line 502 "SB10QD.f"
	*info = 1;
#line 503 "SB10QD.f"
	return 0;
#line 504 "SB10QD.f"
    }

/*     Compute inv(RT) . */

#line 508 "SB10QD.f"
    dsytri_("U", &np1, &dwork[1], np, &iwork[1], &dwork[iwrk], &info2, (
	    ftnlen)1);

/*     Compute -inv(|D1111||D1111' D1121'| - gamma^2*Inp1)*|D1112| . */
/*                  |D1121|                                |D1122| */

#line 513 "SB10QD.f"
    dsymm_("L", "U", &np1, &np2, &c_b17, &dwork[1], np, &d__[(nd2 + 1) * 
	    d_dim1 + 1], ldd, &c_b5, &dwork[np1 * *np + 1], np, (ftnlen)1, (
	    ftnlen)1);

/*     Compute [D1112' D1122']*inv(|D1111||D1111' D1121'| - */
/*                                 |D1121| */

/*                gamma^2*Inp1)*|D1112| + Inp2 . */
/*                              |D1122| */

#line 522 "SB10QD.f"
    dlaset_("Full", &np2, &np2, &c_b5, &c_b8, &dwork[np1 * (*np + 1) + 1], np,
	     (ftnlen)4);
#line 524 "SB10QD.f"
    mb01rx_("Left", "Upper", "Transpose", &np2, &np1, &c_b8, &c_b17, &dwork[
	    np1 * (*np + 1) + 1], np, &d__[(nd2 + 1) * d_dim1 + 1], ldd, &
	    dwork[np1 * *np + 1], np, &info2, (ftnlen)4, (ftnlen)5, (ftnlen)9)
	    ;

/*     Compute B1*D11' . */

#line 530 "SB10QD.f"
    dgemm_("N", "T", n, &np1, &m1, &c_b8, &b[b_offset], ldb, &d__[d_offset], 
	    ldd, &c_b5, &dwork[iw2], n, (ftnlen)1, (ftnlen)1);

/*     Compute B1*DD1' . */

#line 535 "SB10QD.f"
    dlacpy_("Full", n, &np2, &b[(nd2 + 1) * b_dim1 + 1], ldb, &dwork[iw2 + 
	    np1 * *n], n, (ftnlen)4);

/*     Compute B1*DD1'*inv(RT) in H . */

#line 540 "SB10QD.f"
    dsymm_("R", "U", n, np, &c_b8, &dwork[1], np, &dwork[iw2], n, &c_b5, &h__[
	    h_offset], ldh, (ftnlen)1, (ftnlen)1);

/*     Compute Ay = A - B1*DD1'*inv(RT)*C . */

#line 545 "SB10QD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[iwa], n, (ftnlen)4);
#line 546 "SB10QD.f"
    dgemm_("N", "N", n, n, np, &c_b17, &h__[h_offset], ldh, &c__[c_offset], 
	    ldc, &c_b8, &dwork[iwa], n, (ftnlen)1, (ftnlen)1);

/*     Compute Cy = B1*B1' - B1*DD1'*inv(RT)*DD1*B1' . */

#line 551 "SB10QD.f"
    if (nd2 == 0) {
#line 552 "SB10QD.f"
	dlaset_("U", n, n, &c_b5, &c_b5, &dwork[iwq], n, (ftnlen)1);
#line 553 "SB10QD.f"
    } else {
#line 554 "SB10QD.f"
	dsyrk_("U", "N", n, &m1, &c_b8, &b[b_offset], ldb, &c_b5, &dwork[iwq],
		 n, (ftnlen)1, (ftnlen)1);
#line 556 "SB10QD.f"
	mb01rx_("Right", "Upper", "Transpose", n, np, &c_b8, &c_b17, &dwork[
		iwq], n, &h__[h_offset], ldh, &dwork[iw2], n, &info2, (ftnlen)
		5, (ftnlen)5, (ftnlen)9);
#line 558 "SB10QD.f"
    }

/*     Compute Dy = C'*inv(RT)*C . */

#line 562 "SB10QD.f"
    iwrk = iw2;
#line 563 "SB10QD.f"
    i__1 = *n * *np;
#line 563 "SB10QD.f"
    mb01ru_("Upper", "Transpose", n, np, &c_b5, &c_b8, &dwork[iwg], n, &c__[
	    c_offset], ldc, &dwork[1], np, &dwork[iwrk], &i__1, &info2, (
	    ftnlen)5, (ftnlen)9);

/*     Solution of the Riccati equation Ay*Y + Y*Ay' + Cy - Y*Dy*Y = 0 . */
/*     Workspace:  need   NP*NP + 13*N*N + 12*N + 5; */
/*                 prefer larger. */

#line 570 "SB10QD.f"
    iwt = iw2;
#line 571 "SB10QD.f"
    iwv = iwt + nn;
#line 572 "SB10QD.f"
    iwr = iwv + nn;
#line 573 "SB10QD.f"
    iwi = iwr + n2;
#line 574 "SB10QD.f"
    iws = iwi + n2;
#line 575 "SB10QD.f"
    iwrk = iws + (nn << 2);

#line 577 "SB10QD.f"
    i__1 = *ldwork - iwrk + 1;
#line 577 "SB10QD.f"
    sb02rd_("All", "Continuous", "NotUsed", "Transpose", "Upper", "GeneralSc"\
	    "aling", "Stable", "NotFactored", "Original", n, &dwork[iwa], n, &
	    dwork[iwt], n, &dwork[iwv], n, &dwork[iwg], n, &dwork[iwq], n, &y[
	    y_offset], ldy, &sep, &xycond[2], &ferr, &dwork[iwr], &dwork[iwi],
	     &dwork[iws], &n2, &iwork[1], &dwork[iwrk], &i__1, &bwork[1], &
	    info2, (ftnlen)3, (ftnlen)10, (ftnlen)7, (ftnlen)9, (ftnlen)5, (
	    ftnlen)14, (ftnlen)6, (ftnlen)11, (ftnlen)8);
#line 584 "SB10QD.f"
    if (info2 > 0) {
#line 585 "SB10QD.f"
	*info = 3;
#line 586 "SB10QD.f"
	return 0;
#line 587 "SB10QD.f"
    }

/* Computing MAX */
#line 589 "SB10QD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 589 "SB10QD.f"
    lwamax = max(i__1,lwamax);

/*     Compute H = -|B1*DD1' + Y*C'|*inv(RT) . */

#line 593 "SB10QD.f"
    iwrk = iw2;
#line 594 "SB10QD.f"
    dgemm_("N", "T", n, np, n, &c_b8, &y[y_offset], ldy, &c__[c_offset], ldc, 
	    &c_b5, &dwork[iwrk], n, (ftnlen)1, (ftnlen)1);
#line 596 "SB10QD.f"
    dsymm_("R", "U", n, np, &c_b17, &dwork[1], np, &dwork[iwrk], n, &c_b17, &
	    h__[h_offset], ldh, (ftnlen)1, (ftnlen)1);

#line 599 "SB10QD.f"
    dwork[1] = (doublereal) lwamax;
#line 600 "SB10QD.f"
    return 0;
/* *** Last line of SB10QD *** */
} /* sb10qd_ */

