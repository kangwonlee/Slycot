#line 1 "SB10RD.f"
/* SB10RD.f -- translated by f2c (version 20100827).
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

#line 1 "SB10RD.f"
/* Table of constant values */

static doublereal c_b7 = 0.;
static doublereal c_b8 = 1.;
static doublereal c_b19 = -1.;

/* Subroutine */ int sb10rd_(integer *n, integer *m, integer *np, integer *
	ncon, integer *nmeas, doublereal *gamma, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *f, integer *ldf, 
	doublereal *h__, integer *ldh, doublereal *tu, integer *ldtu, 
	doublereal *ty, integer *ldty, doublereal *x, integer *ldx, 
	doublereal *y, integer *ldy, doublereal *ak, integer *ldak, 
	doublereal *bk, integer *ldbk, doublereal *ck, integer *ldck, 
	doublereal *dk, integer *lddk, integer *iwork, doublereal *dwork, 
	integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ak_dim1, ak_offset, b_dim1, b_offset, bk_dim1, 
	    bk_offset, c_dim1, c_offset, ck_dim1, ck_offset, d_dim1, d_offset,
	     dk_dim1, dk_offset, f_dim1, f_offset, h_dim1, h_offset, tu_dim1, 
	    tu_offset, ty_dim1, ty_offset, x_dim1, x_offset, y_dim1, y_offset,
	     i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10, 
	    i__11, i__12;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, m1, m2, ij, nd1, nd2, np1, np2, iw1, iw2, iw3, iw4,
	     id11, id12, id21, iwb, iwc;
    static doublereal eps;
    static integer iwrk, info2;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), dgemm_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen);
    static doublereal rcond;
    extern /* Subroutine */ int mb01rx_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dsyrk_(
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgecon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dgetrf_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), dlacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), dlaset_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), dgetri_(integer *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *), xerbla_(char *, integer *, ftnlen), dgetrs_(char *, 
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);
    static integer lwamax;
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int dpotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dsycon_(char *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
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

/*     To compute the matrices of an H-infinity (sub)optimal controller */

/*              | AK | BK | */
/*          K = |----|----|, */
/*              | CK | DK | */

/*     from the state feedback matrix F and output injection matrix H as */
/*     determined by the SLICOT Library routine SB10QD. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The column size of the matrix B.  M >= 0. */

/*     NP      (input) INTEGER */
/*             The row size of the matrix C.  NP >= 0. */

/*     NCON    (input) INTEGER */
/*             The number of control inputs (M2).  M >= NCON >= 0. */
/*             NP-NMEAS >= NCON. */

/*     NMEAS   (input) INTEGER */
/*             The number of measurements (NP2).  NP >= NMEAS >= 0. */
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

/*     F       (input) DOUBLE PRECISION array, dimension (LDF,N) */
/*             The leading M-by-N part of this array must contain the */
/*             state feedback matrix F. */

/*     LDF     INTEGER */
/*             The leading dimension of the array F.  LDF >= max(1,M). */

/*     H       (input) DOUBLE PRECISION array, dimension (LDH,NP) */
/*             The leading N-by-NP part of this array must contain the */
/*             output injection matrix H. */

/*     LDH     INTEGER */
/*             The leading dimension of the array H.  LDH >= max(1,N). */

/*     TU      (input) DOUBLE PRECISION array, dimension (LDTU,M2) */
/*             The leading M2-by-M2 part of this array must contain the */
/*             control transformation matrix TU, as obtained by the */
/*             SLICOT Library routine SB10PD. */

/*     LDTU    INTEGER */
/*             The leading dimension of the array TU.  LDTU >= max(1,M2). */

/*     TY      (input) DOUBLE PRECISION array, dimension (LDTY,NP2) */
/*             The leading NP2-by-NP2 part of this array must contain the */
/*             measurement transformation matrix TY, as obtained by the */
/*             SLICOT Library routine SB10PD. */

/*     LDTY    INTEGER */
/*             The leading dimension of the array TY. */
/*             LDTY >= max(1,NP2). */

/*     X       (input) DOUBLE PRECISION array, dimension (LDX,N) */
/*             The leading N-by-N part of this array must contain the */
/*             matrix X, solution of the X-Riccati equation, as obtained */
/*             by the SLICOT Library routine SB10QD. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X.  LDX >= max(1,N). */

/*     Y       (input) DOUBLE PRECISION array, dimension (LDY,N) */
/*             The leading N-by-N part of this array must contain the */
/*             matrix Y, solution of the Y-Riccati equation, as obtained */
/*             by the SLICOT Library routine SB10QD. */

/*     LDY     INTEGER */
/*             The leading dimension of the array Y.  LDY >= max(1,N). */

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

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK), where */
/*             LIWORK = max(2*(max(NP,M)-M2-NP2,M2,N),NP2) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal */
/*             LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= max(1, M2*NP2 + NP2*NP2 + M2*M2 + */
/*                           max(D1*D1 + max(2*D1, (D1+D2)*NP2), */
/*                               D2*D2 + max(2*D2, D2*M2), 3*N, */
/*                               N*(2*NP2 + M2) + */
/*                               max(2*N*M2, M2*NP2 + */
/*                                           max(M2*M2+3*M2, NP2*(2*NP2+ */
/*                                                  M2+max(NP2,N)))))) */
/*             where D1 = NP1 - M2, D2 = M1 - NP2, */
/*                  NP1 = NP - NP2, M1 = M - M2. */
/*             For good performance, LDWORK must generally be larger. */
/*             Denoting Q = max(M1,M2,NP1,NP2), an upper bound is */
/*             max( 1, Q*(3*Q + 3*N + max(2*N, 4*Q + max(Q, N)))). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the controller is not admissible (too small value */
/*                   of gamma); */
/*             = 2:  if the determinant of Im2 + Tu*D11HAT*Ty*D22 is zero. */

/*     METHOD */

/*     The routine implements the Glover's and Doyle's formulas [1],[2]. */

/*     REFERENCES */

/*     [1] Glover, K. and Doyle, J.C. */
/*         State-space formulae for all stabilizing controllers that */
/*         satisfy an Hinf norm bound and relations to risk sensitivity. */
/*         Systems and Control Letters, vol. 11, pp. 167-172, 1988. */

/*     [2] Balas, G.J., Doyle, J.C., Glover, K., Packard, A., and */
/*         Smith, R. */
/*         mu-Analysis and Synthesis Toolbox. */
/*         The MathWorks Inc., Natick, Mass., 1995. */

/*     NUMERICAL ASPECTS */

/*     The accuracy of the result depends on the condition numbers of the */
/*     input and output transformations. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 1998. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 1999, */
/*     Sept. 1999, Oct. 2001. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, H-infinity optimal control, robust */
/*     control. */

/*  ********************************************************************* */

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

#line 277 "SB10RD.f"
    /* Parameter adjustments */
#line 277 "SB10RD.f"
    a_dim1 = *lda;
#line 277 "SB10RD.f"
    a_offset = 1 + a_dim1;
#line 277 "SB10RD.f"
    a -= a_offset;
#line 277 "SB10RD.f"
    b_dim1 = *ldb;
#line 277 "SB10RD.f"
    b_offset = 1 + b_dim1;
#line 277 "SB10RD.f"
    b -= b_offset;
#line 277 "SB10RD.f"
    c_dim1 = *ldc;
#line 277 "SB10RD.f"
    c_offset = 1 + c_dim1;
#line 277 "SB10RD.f"
    c__ -= c_offset;
#line 277 "SB10RD.f"
    d_dim1 = *ldd;
#line 277 "SB10RD.f"
    d_offset = 1 + d_dim1;
#line 277 "SB10RD.f"
    d__ -= d_offset;
#line 277 "SB10RD.f"
    f_dim1 = *ldf;
#line 277 "SB10RD.f"
    f_offset = 1 + f_dim1;
#line 277 "SB10RD.f"
    f -= f_offset;
#line 277 "SB10RD.f"
    h_dim1 = *ldh;
#line 277 "SB10RD.f"
    h_offset = 1 + h_dim1;
#line 277 "SB10RD.f"
    h__ -= h_offset;
#line 277 "SB10RD.f"
    tu_dim1 = *ldtu;
#line 277 "SB10RD.f"
    tu_offset = 1 + tu_dim1;
#line 277 "SB10RD.f"
    tu -= tu_offset;
#line 277 "SB10RD.f"
    ty_dim1 = *ldty;
#line 277 "SB10RD.f"
    ty_offset = 1 + ty_dim1;
#line 277 "SB10RD.f"
    ty -= ty_offset;
#line 277 "SB10RD.f"
    x_dim1 = *ldx;
#line 277 "SB10RD.f"
    x_offset = 1 + x_dim1;
#line 277 "SB10RD.f"
    x -= x_offset;
#line 277 "SB10RD.f"
    y_dim1 = *ldy;
#line 277 "SB10RD.f"
    y_offset = 1 + y_dim1;
#line 277 "SB10RD.f"
    y -= y_offset;
#line 277 "SB10RD.f"
    ak_dim1 = *ldak;
#line 277 "SB10RD.f"
    ak_offset = 1 + ak_dim1;
#line 277 "SB10RD.f"
    ak -= ak_offset;
#line 277 "SB10RD.f"
    bk_dim1 = *ldbk;
#line 277 "SB10RD.f"
    bk_offset = 1 + bk_dim1;
#line 277 "SB10RD.f"
    bk -= bk_offset;
#line 277 "SB10RD.f"
    ck_dim1 = *ldck;
#line 277 "SB10RD.f"
    ck_offset = 1 + ck_dim1;
#line 277 "SB10RD.f"
    ck -= ck_offset;
#line 277 "SB10RD.f"
    dk_dim1 = *lddk;
#line 277 "SB10RD.f"
    dk_offset = 1 + dk_dim1;
#line 277 "SB10RD.f"
    dk -= dk_offset;
#line 277 "SB10RD.f"
    --iwork;
#line 277 "SB10RD.f"
    --dwork;
#line 277 "SB10RD.f"

#line 277 "SB10RD.f"
    /* Function Body */
#line 277 "SB10RD.f"
    m1 = *m - *ncon;
#line 278 "SB10RD.f"
    m2 = *ncon;
#line 279 "SB10RD.f"
    np1 = *np - *nmeas;
#line 280 "SB10RD.f"
    np2 = *nmeas;

#line 282 "SB10RD.f"
    *info = 0;
#line 283 "SB10RD.f"
    if (*n < 0) {
#line 284 "SB10RD.f"
	*info = -1;
#line 285 "SB10RD.f"
    } else if (*m < 0) {
#line 286 "SB10RD.f"
	*info = -2;
#line 287 "SB10RD.f"
    } else if (*np < 0) {
#line 288 "SB10RD.f"
	*info = -3;
#line 289 "SB10RD.f"
    } else if (*ncon < 0 || m1 < 0 || m2 > np1) {
#line 290 "SB10RD.f"
	*info = -4;
#line 291 "SB10RD.f"
    } else if (*nmeas < 0 || np1 < 0 || np2 > m1) {
#line 292 "SB10RD.f"
	*info = -5;
#line 293 "SB10RD.f"
    } else if (*gamma < 0.) {
#line 294 "SB10RD.f"
	*info = -6;
#line 295 "SB10RD.f"
    } else if (*lda < max(1,*n)) {
#line 296 "SB10RD.f"
	*info = -8;
#line 297 "SB10RD.f"
    } else if (*ldb < max(1,*n)) {
#line 298 "SB10RD.f"
	*info = -10;
#line 299 "SB10RD.f"
    } else if (*ldc < max(1,*np)) {
#line 300 "SB10RD.f"
	*info = -12;
#line 301 "SB10RD.f"
    } else if (*ldd < max(1,*np)) {
#line 302 "SB10RD.f"
	*info = -14;
#line 303 "SB10RD.f"
    } else if (*ldf < max(1,*m)) {
#line 304 "SB10RD.f"
	*info = -16;
#line 305 "SB10RD.f"
    } else if (*ldh < max(1,*n)) {
#line 306 "SB10RD.f"
	*info = -18;
#line 307 "SB10RD.f"
    } else if (*ldtu < max(1,m2)) {
#line 308 "SB10RD.f"
	*info = -20;
#line 309 "SB10RD.f"
    } else if (*ldty < max(1,np2)) {
#line 310 "SB10RD.f"
	*info = -22;
#line 311 "SB10RD.f"
    } else if (*ldx < max(1,*n)) {
#line 312 "SB10RD.f"
	*info = -24;
#line 313 "SB10RD.f"
    } else if (*ldy < max(1,*n)) {
#line 314 "SB10RD.f"
	*info = -26;
#line 315 "SB10RD.f"
    } else if (*ldak < max(1,*n)) {
#line 316 "SB10RD.f"
	*info = -28;
#line 317 "SB10RD.f"
    } else if (*ldbk < max(1,*n)) {
#line 318 "SB10RD.f"
	*info = -30;
#line 319 "SB10RD.f"
    } else if (*ldck < max(1,m2)) {
#line 320 "SB10RD.f"
	*info = -32;
#line 321 "SB10RD.f"
    } else if (*lddk < max(1,m2)) {
#line 322 "SB10RD.f"
	*info = -34;
#line 323 "SB10RD.f"
    } else {

/*        Compute workspace. */

#line 327 "SB10RD.f"
	nd1 = np1 - m2;
#line 328 "SB10RD.f"
	nd2 = m1 - np2;
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
#line 329 "SB10RD.f"
	i__5 = nd1 << 1, i__6 = (nd1 + nd2) * np2;
/* Computing MAX */
#line 329 "SB10RD.f"
	i__7 = nd2 << 1, i__8 = nd2 * m2;
/* Computing MAX */
/* Computing MAX */
#line 329 "SB10RD.f"
	i__11 = m2 * m2 + m2 * 3, i__12 = np2 * ((np2 << 1) + m2 + max(np2,*n)
		);
#line 329 "SB10RD.f"
	i__9 = (*n << 1) * m2, i__10 = m2 * np2 + max(i__11,i__12);
#line 329 "SB10RD.f"
	i__3 = nd1 * nd1 + max(i__5,i__6), i__4 = nd2 * nd2 + max(i__7,i__8), 
		i__3 = max(i__3,i__4), i__4 = *n * 3, i__3 = max(i__3,i__4), 
		i__4 = *n * ((np2 << 1) + m2) + max(i__9,i__10);
#line 329 "SB10RD.f"
	i__1 = 1, i__2 = m2 * np2 + np2 * np2 + m2 * m2 + max(i__3,i__4);
#line 329 "SB10RD.f"
	minwrk = max(i__1,i__2);
#line 336 "SB10RD.f"
	if (*ldwork < minwrk) {
#line 336 "SB10RD.f"
	    *info = -37;
#line 336 "SB10RD.f"
	}
#line 338 "SB10RD.f"
    }
#line 339 "SB10RD.f"
    if (*info != 0) {
#line 340 "SB10RD.f"
	i__1 = -(*info);
#line 340 "SB10RD.f"
	xerbla_("SB10RD", &i__1, (ftnlen)6);
#line 341 "SB10RD.f"
	return 0;
#line 342 "SB10RD.f"
    }

/*     Quick return if possible. */

#line 346 "SB10RD.f"
    if (*n == 0 || *m == 0 || *np == 0 || m1 == 0 || m2 == 0 || np1 == 0 || 
	    np2 == 0) {
#line 348 "SB10RD.f"
	dwork[1] = 1.;
#line 349 "SB10RD.f"
	return 0;
#line 350 "SB10RD.f"
    }

/*     Get the machine precision. */

#line 354 "SB10RD.f"
    eps = dlamch_("Epsilon", (ftnlen)7);

/*     Workspace usage. */

#line 358 "SB10RD.f"
    id11 = 1;
#line 359 "SB10RD.f"
    id21 = id11 + m2 * np2;
#line 360 "SB10RD.f"
    id12 = id21 + np2 * np2;
#line 361 "SB10RD.f"
    iw1 = id12 + m2 * m2;
#line 362 "SB10RD.f"
    iw2 = iw1 + nd1 * nd1;
#line 363 "SB10RD.f"
    iw3 = iw2 + nd1 * np2;
#line 364 "SB10RD.f"
    iwrk = iw2;

/*     Set D11HAT := -D1122 . */

#line 368 "SB10RD.f"
    ij = id11;
#line 369 "SB10RD.f"
    i__1 = np2;
#line 369 "SB10RD.f"
    for (j = 1; j <= i__1; ++j) {
#line 370 "SB10RD.f"
	i__2 = m2;
#line 370 "SB10RD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 371 "SB10RD.f"
	    dwork[ij] = -d__[nd1 + i__ + (nd2 + j) * d_dim1];
#line 372 "SB10RD.f"
	    ++ij;
#line 373 "SB10RD.f"
/* L10: */
#line 373 "SB10RD.f"
	}
#line 374 "SB10RD.f"
/* L20: */
#line 374 "SB10RD.f"
    }

/*     Set D21HAT := Inp2 . */

#line 378 "SB10RD.f"
    dlaset_("Upper", &np2, &np2, &c_b7, &c_b8, &dwork[id21], &np2, (ftnlen)5);

/*     Set D12HAT := Im2 . */

#line 382 "SB10RD.f"
    dlaset_("Lower", &m2, &m2, &c_b7, &c_b8, &dwork[id12], &m2, (ftnlen)5);

/*     Compute D11HAT, D21HAT, D12HAT . */

#line 386 "SB10RD.f"
    lwamax = 0;
#line 387 "SB10RD.f"
    if (nd1 > 0) {
#line 388 "SB10RD.f"
	if (nd2 == 0) {

/*           Compute D21HAT'*D21HAT = Inp2 - D1112'*D1112/gamma^2 . */

/* Computing 2nd power */
#line 392 "SB10RD.f"
	    d__2 = *gamma;
#line 392 "SB10RD.f"
	    d__1 = -1. / (d__2 * d__2);
#line 392 "SB10RD.f"
	    dsyrk_("U", "T", &np2, &nd1, &d__1, &d__[d_offset], ldd, &c_b8, &
		    dwork[id21], &np2, (ftnlen)1, (ftnlen)1);
#line 394 "SB10RD.f"
	} else {

/*           Compute gdum = gamma^2*Ind1 - D1111*D1111' . */

/* Computing 2nd power */
#line 398 "SB10RD.f"
	    d__2 = *gamma;
#line 398 "SB10RD.f"
	    d__1 = d__2 * d__2;
#line 398 "SB10RD.f"
	    dlaset_("U", &nd1, &nd1, &c_b7, &d__1, &dwork[iw1], &nd1, (ftnlen)
		    1);
#line 400 "SB10RD.f"
	    dsyrk_("U", "N", &nd1, &nd2, &c_b19, &d__[d_offset], ldd, &c_b8, &
		    dwork[iw1], &nd1, (ftnlen)1, (ftnlen)1);
#line 402 "SB10RD.f"
	    anorm = dlansy_("I", "U", &nd1, &dwork[iw1], &nd1, &dwork[iwrk], (
		    ftnlen)1, (ftnlen)1);
#line 404 "SB10RD.f"
	    i__1 = *ldwork - iwrk + 1;
#line 404 "SB10RD.f"
	    dsytrf_("U", &nd1, &dwork[iw1], &nd1, &iwork[1], &dwork[iwrk], &
		    i__1, &info2, (ftnlen)1);
#line 406 "SB10RD.f"
	    if (info2 > 0) {
#line 407 "SB10RD.f"
		*info = 1;
#line 408 "SB10RD.f"
		return 0;
#line 409 "SB10RD.f"
	    }
#line 410 "SB10RD.f"
	    lwamax = (integer) dwork[iwrk] + iwrk - 1;
#line 411 "SB10RD.f"
	    dsycon_("U", &nd1, &dwork[iw1], &nd1, &iwork[1], &anorm, &rcond, &
		    dwork[iwrk], &iwork[nd1 + 1], &info2, (ftnlen)1);

/*           Return if the matrix is singular to working precision. */

#line 416 "SB10RD.f"
	    if (rcond < eps) {
#line 417 "SB10RD.f"
		*info = 1;
#line 418 "SB10RD.f"
		return 0;
#line 419 "SB10RD.f"
	    }

/*           Compute inv(gdum)*D1112 . */

#line 423 "SB10RD.f"
	    dlacpy_("Full", &nd1, &np2, &d__[(nd2 + 1) * d_dim1 + 1], ldd, &
		    dwork[iw2], &nd1, (ftnlen)4);
#line 425 "SB10RD.f"
	    dsytrs_("U", &nd1, &np2, &dwork[iw1], &nd1, &iwork[1], &dwork[iw2]
		    , &nd1, &info2, (ftnlen)1);

/*           Compute D11HAT = -D1121*D1111'*inv(gdum)*D1112 - D1122 . */

#line 430 "SB10RD.f"
	    dgemm_("T", "N", &nd2, &np2, &nd1, &c_b8, &d__[d_offset], ldd, &
		    dwork[iw2], &nd1, &c_b7, &dwork[iw3], &nd2, (ftnlen)1, (
		    ftnlen)1);
#line 432 "SB10RD.f"
	    dgemm_("N", "N", &m2, &np2, &nd2, &c_b19, &d__[nd1 + 1 + d_dim1], 
		    ldd, &dwork[iw3], &nd2, &c_b8, &dwork[id11], &m2, (ftnlen)
		    1, (ftnlen)1);

/*           Compute D21HAT'*D21HAT = Inp2 - D1112'*inv(gdum)*D1112 . */

#line 437 "SB10RD.f"
	    mb01rx_("Left", "Upper", "Transpose", &np2, &nd1, &c_b8, &c_b19, &
		    dwork[id21], &np2, &d__[(nd2 + 1) * d_dim1 + 1], ldd, &
		    dwork[iw2], &nd1, &info2, (ftnlen)4, (ftnlen)5, (ftnlen)9)
		    ;

#line 441 "SB10RD.f"
	    iw2 = iw1 + nd2 * nd2;
#line 442 "SB10RD.f"
	    iwrk = iw2;

/*           Compute gdum = gamma^2*Ind2 - D1111'*D1111 . */

/* Computing 2nd power */
#line 446 "SB10RD.f"
	    d__2 = *gamma;
#line 446 "SB10RD.f"
	    d__1 = d__2 * d__2;
#line 446 "SB10RD.f"
	    dlaset_("L", &nd2, &nd2, &c_b7, &d__1, &dwork[iw1], &nd2, (ftnlen)
		    1);
#line 448 "SB10RD.f"
	    dsyrk_("L", "T", &nd2, &nd1, &c_b19, &d__[d_offset], ldd, &c_b8, &
		    dwork[iw1], &nd2, (ftnlen)1, (ftnlen)1);
#line 450 "SB10RD.f"
	    anorm = dlansy_("I", "L", &nd2, &dwork[iw1], &nd2, &dwork[iwrk], (
		    ftnlen)1, (ftnlen)1);
#line 452 "SB10RD.f"
	    i__1 = *ldwork - iwrk + 1;
#line 452 "SB10RD.f"
	    dsytrf_("L", &nd2, &dwork[iw1], &nd2, &iwork[1], &dwork[iwrk], &
		    i__1, &info2, (ftnlen)1);
#line 454 "SB10RD.f"
	    if (info2 > 0) {
#line 455 "SB10RD.f"
		*info = 1;
#line 456 "SB10RD.f"
		return 0;
#line 457 "SB10RD.f"
	    }
/* Computing MAX */
#line 458 "SB10RD.f"
	    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 458 "SB10RD.f"
	    lwamax = max(i__1,lwamax);
#line 459 "SB10RD.f"
	    dsycon_("L", &nd2, &dwork[iw1], &nd2, &iwork[1], &anorm, &rcond, &
		    dwork[iwrk], &iwork[nd2 + 1], &info2, (ftnlen)1);

/*           Return if the matrix is singular to working precision. */

#line 464 "SB10RD.f"
	    if (rcond < eps) {
#line 465 "SB10RD.f"
		*info = 1;
#line 466 "SB10RD.f"
		return 0;
#line 467 "SB10RD.f"
	    }

/*           Compute inv(gdum)*D1121' . */

#line 471 "SB10RD.f"
	    ma02ad_("Full", &m2, &nd2, &d__[nd1 + 1 + d_dim1], ldd, &dwork[
		    iw2], &nd2, (ftnlen)4);
#line 473 "SB10RD.f"
	    dsytrs_("L", &nd2, &m2, &dwork[iw1], &nd2, &iwork[1], &dwork[iw2],
		     &nd2, &info2, (ftnlen)1);

/*           Compute D12HAT*D12HAT' = Im2 - D1121*inv(gdum)*D1121' . */

#line 478 "SB10RD.f"
	    mb01rx_("Left", "Lower", "NoTranspose", &m2, &nd2, &c_b8, &c_b19, 
		    &dwork[id12], &m2, &d__[nd1 + 1 + d_dim1], ldd, &dwork[
		    iw2], &nd2, &info2, (ftnlen)4, (ftnlen)5, (ftnlen)11);
#line 481 "SB10RD.f"
	}
#line 482 "SB10RD.f"
    } else {
#line 483 "SB10RD.f"
	if (nd2 > 0) {

/*           Compute D12HAT*D12HAT' = Im2 - D1121*D1121'/gamma^2 . */

/* Computing 2nd power */
#line 487 "SB10RD.f"
	    d__2 = *gamma;
#line 487 "SB10RD.f"
	    d__1 = -1. / (d__2 * d__2);
#line 487 "SB10RD.f"
	    dsyrk_("L", "N", &m2, &nd2, &d__1, &d__[d_offset], ldd, &c_b8, &
		    dwork[id12], &m2, (ftnlen)1, (ftnlen)1);
#line 489 "SB10RD.f"
	}
#line 490 "SB10RD.f"
    }

/*     Compute D21HAT using Cholesky decomposition. */

#line 494 "SB10RD.f"
    dpotrf_("U", &np2, &dwork[id21], &np2, &info2, (ftnlen)1);
#line 495 "SB10RD.f"
    if (info2 > 0) {
#line 496 "SB10RD.f"
	*info = 1;
#line 497 "SB10RD.f"
	return 0;
#line 498 "SB10RD.f"
    }

/*     Compute D12HAT using Cholesky decomposition. */

#line 502 "SB10RD.f"
    dpotrf_("L", &m2, &dwork[id12], &m2, &info2, (ftnlen)1);
#line 503 "SB10RD.f"
    if (info2 > 0) {
#line 504 "SB10RD.f"
	*info = 1;
#line 505 "SB10RD.f"
	return 0;
#line 506 "SB10RD.f"
    }
/*             _ */
/*     Compute Z = In - Y*X/gamma^2 and its LU factorization in AK . */

#line 510 "SB10RD.f"
    iwrk = iw1;
#line 511 "SB10RD.f"
    dlaset_("Full", n, n, &c_b7, &c_b8, &ak[ak_offset], ldak, (ftnlen)4);
/* Computing 2nd power */
#line 512 "SB10RD.f"
    d__2 = *gamma;
#line 512 "SB10RD.f"
    d__1 = -1. / (d__2 * d__2);
#line 512 "SB10RD.f"
    dgemm_("N", "N", n, n, n, &d__1, &y[y_offset], ldy, &x[x_offset], ldx, &
	    c_b8, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);
#line 514 "SB10RD.f"
    anorm = dlange_("1", n, n, &ak[ak_offset], ldak, &dwork[iwrk], (ftnlen)1);
#line 515 "SB10RD.f"
    dgetrf_(n, n, &ak[ak_offset], ldak, &iwork[1], &info2);
#line 516 "SB10RD.f"
    if (info2 > 0) {
#line 517 "SB10RD.f"
	*info = 1;
#line 518 "SB10RD.f"
	return 0;
#line 519 "SB10RD.f"
    }
#line 520 "SB10RD.f"
    dgecon_("1", n, &ak[ak_offset], ldak, &anorm, &rcond, &dwork[iwrk], &
	    iwork[*n + 1], info, (ftnlen)1);

/*     Return if the matrix is singular to working precision. */

#line 525 "SB10RD.f"
    if (rcond < eps) {
#line 526 "SB10RD.f"
	*info = 1;
#line 527 "SB10RD.f"
	return 0;
#line 528 "SB10RD.f"
    }

#line 530 "SB10RD.f"
    iwb = iw1;
#line 531 "SB10RD.f"
    iwc = iwb + *n * np2;
#line 532 "SB10RD.f"
    iw1 = iwc + (m2 + np2) * *n;
#line 533 "SB10RD.f"
    iw2 = iw1 + *n * m2;

/*     Compute C2' + F12' in BK . */

#line 537 "SB10RD.f"
    i__1 = *n;
#line 537 "SB10RD.f"
    for (j = 1; j <= i__1; ++j) {
#line 538 "SB10RD.f"
	i__2 = np2;
#line 538 "SB10RD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 539 "SB10RD.f"
	    bk[j + i__ * bk_dim1] = c__[np1 + i__ + j * c_dim1] + f[nd2 + i__ 
		    + j * f_dim1];
#line 540 "SB10RD.f"
/* L30: */
#line 540 "SB10RD.f"
	}
#line 541 "SB10RD.f"
/* L40: */
#line 541 "SB10RD.f"
    }
/*                                                          _ */
/*     Compute the transpose of (C2 + F12)*Z , with Z = inv(Z) . */

#line 545 "SB10RD.f"
    dgetrs_("Transpose", n, &np2, &ak[ak_offset], ldak, &iwork[1], &bk[
	    bk_offset], ldbk, &info2, (ftnlen)9);

/*     Compute the transpose of F2*Z . */

#line 550 "SB10RD.f"
    ma02ad_("Full", &m2, n, &f[m1 + 1 + f_dim1], ldf, &dwork[iw1], n, (ftnlen)
	    4);
#line 551 "SB10RD.f"
    dgetrs_("Transpose", n, &m2, &ak[ak_offset], ldak, &iwork[1], &dwork[iw1],
	     n, &info2, (ftnlen)9);

/*     Compute the transpose of C1HAT = F2*Z - D11HAT*(C2 + F12)*Z . */

#line 556 "SB10RD.f"
    dgemm_("N", "T", n, &m2, &np2, &c_b19, &bk[bk_offset], ldbk, &dwork[id11],
	     &m2, &c_b8, &dwork[iw1], n, (ftnlen)1, (ftnlen)1);

/*     Compute CHAT . */

#line 561 "SB10RD.f"
    i__1 = m2 + np2;
#line 561 "SB10RD.f"
    dgemm_("N", "T", &m2, n, &m2, &c_b8, &tu[tu_offset], ldtu, &dwork[iw1], n,
	     &c_b7, &dwork[iwc], &i__1, (ftnlen)1, (ftnlen)1);
#line 563 "SB10RD.f"
    i__1 = m2 + np2;
#line 563 "SB10RD.f"
    ma02ad_("Full", n, &np2, &bk[bk_offset], ldbk, &dwork[iwc + m2], &i__1, (
	    ftnlen)4);
#line 564 "SB10RD.f"
    i__1 = m2 + np2;
#line 564 "SB10RD.f"
    dtrmm_("L", "U", "N", "N", &np2, n, &c_b19, &dwork[id21], &np2, &dwork[
	    iwc + m2], &i__1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Compute B2 + H12 . */

#line 569 "SB10RD.f"
    ij = iw2;
#line 570 "SB10RD.f"
    i__1 = m2;
#line 570 "SB10RD.f"
    for (j = 1; j <= i__1; ++j) {
#line 571 "SB10RD.f"
	i__2 = *n;
#line 571 "SB10RD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 572 "SB10RD.f"
	    dwork[ij] = b[i__ + (m1 + j) * b_dim1] + h__[i__ + (nd1 + j) * 
		    h_dim1];
#line 573 "SB10RD.f"
	    ++ij;
#line 574 "SB10RD.f"
/* L50: */
#line 574 "SB10RD.f"
	}
#line 575 "SB10RD.f"
/* L60: */
#line 575 "SB10RD.f"
    }

/*     Compute A + HC in AK . */

#line 579 "SB10RD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &ak[ak_offset], ldak, (ftnlen)4);
#line 580 "SB10RD.f"
    dgemm_("N", "N", n, n, np, &c_b8, &h__[h_offset], ldh, &c__[c_offset], 
	    ldc, &c_b8, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);

/*     Compute AHAT = A + HC + (B2 + H12)*C1HAT in AK . */

#line 585 "SB10RD.f"
    dgemm_("N", "T", n, n, &m2, &c_b8, &dwork[iw2], n, &dwork[iw1], n, &c_b8, 
	    &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);

/*     Compute B1HAT = -H2 + (B2 + H12)*D11HAT in BK . */

#line 590 "SB10RD.f"
    dlacpy_("Full", n, &np2, &h__[(np1 + 1) * h_dim1 + 1], ldh, &bk[bk_offset]
	    , ldbk, (ftnlen)4);
#line 591 "SB10RD.f"
    dgemm_("N", "N", n, &np2, &m2, &c_b8, &dwork[iw2], n, &dwork[id11], &m2, &
	    c_b19, &bk[bk_offset], ldbk, (ftnlen)1, (ftnlen)1);

/*     Compute the first block of BHAT, BHAT1 . */

#line 596 "SB10RD.f"
    dgemm_("N", "N", n, &np2, &np2, &c_b8, &bk[bk_offset], ldbk, &ty[
	    ty_offset], ldty, &c_b7, &dwork[iwb], n, (ftnlen)1, (ftnlen)1);

/*     Compute Tu*D11HAT . */

#line 601 "SB10RD.f"
    dgemm_("N", "N", &m2, &np2, &m2, &c_b8, &tu[tu_offset], ldtu, &dwork[id11]
	    , &m2, &c_b7, &dwork[iw1], &m2, (ftnlen)1, (ftnlen)1);

/*     Compute Tu*D11HAT*Ty in DK . */

#line 606 "SB10RD.f"
    dgemm_("N", "N", &m2, &np2, &np2, &c_b8, &dwork[iw1], &m2, &ty[ty_offset],
	     ldty, &c_b7, &dk[dk_offset], lddk, (ftnlen)1, (ftnlen)1);

/*     Compute P = Im2 + Tu*D11HAT*Ty*D22 and its condition. */

#line 611 "SB10RD.f"
    iw2 = iw1 + m2 * np2;
#line 612 "SB10RD.f"
    iwrk = iw2 + m2 * m2;
#line 613 "SB10RD.f"
    dlaset_("Full", &m2, &m2, &c_b7, &c_b8, &dwork[iw2], &m2, (ftnlen)4);
#line 614 "SB10RD.f"
    dgemm_("N", "N", &m2, &m2, &np2, &c_b8, &dk[dk_offset], lddk, &d__[np1 + 
	    1 + (m1 + 1) * d_dim1], ldd, &c_b8, &dwork[iw2], &m2, (ftnlen)1, (
	    ftnlen)1);
#line 616 "SB10RD.f"
    anorm = dlange_("1", &m2, &m2, &dwork[iw2], &m2, &dwork[iwrk], (ftnlen)1);
#line 617 "SB10RD.f"
    dgetrf_(&m2, &m2, &dwork[iw2], &m2, &iwork[1], &info2);
#line 618 "SB10RD.f"
    if (info2 > 0) {
#line 619 "SB10RD.f"
	*info = 2;
#line 620 "SB10RD.f"
	return 0;
#line 621 "SB10RD.f"
    }
#line 622 "SB10RD.f"
    dgecon_("1", &m2, &dwork[iw2], &m2, &anorm, &rcond, &dwork[iwrk], &iwork[
	    m2 + 1], &info2, (ftnlen)1);

/*     Return if the matrix is singular to working precision. */

#line 627 "SB10RD.f"
    if (rcond < eps) {
#line 628 "SB10RD.f"
	*info = 2;
#line 629 "SB10RD.f"
	return 0;
#line 630 "SB10RD.f"
    }

/*     Find the controller matrix CK, CK = inv(P)*CHAT(1:M2,:) . */

#line 634 "SB10RD.f"
    i__1 = m2 + np2;
#line 634 "SB10RD.f"
    dlacpy_("Full", &m2, n, &dwork[iwc], &i__1, &ck[ck_offset], ldck, (ftnlen)
	    4);
#line 635 "SB10RD.f"
    dgetrs_("NoTranspose", &m2, n, &dwork[iw2], &m2, &iwork[1], &ck[ck_offset]
	    , ldck, &info2, (ftnlen)11);

/*     Find the controller matrices AK, BK, and DK, exploiting the */
/*     special structure of the relations. */

/*     Compute Q = Inp2 + D22*Tu*D11HAT*Ty and its LU factorization. */

#line 643 "SB10RD.f"
    iw3 = iw2 + np2 * np2;
#line 644 "SB10RD.f"
    iw4 = iw3 + np2 * m2;
#line 645 "SB10RD.f"
    iwrk = iw4 + np2 * np2;
#line 646 "SB10RD.f"
    dlaset_("Full", &np2, &np2, &c_b7, &c_b8, &dwork[iw2], &np2, (ftnlen)4);
#line 647 "SB10RD.f"
    dgemm_("N", "N", &np2, &np2, &m2, &c_b8, &d__[np1 + 1 + (m1 + 1) * d_dim1]
	    , ldd, &dk[dk_offset], lddk, &c_b8, &dwork[iw2], &np2, (ftnlen)1, 
	    (ftnlen)1);
#line 649 "SB10RD.f"
    dgetrf_(&np2, &np2, &dwork[iw2], &np2, &iwork[1], &info2);
#line 650 "SB10RD.f"
    if (info2 > 0) {
#line 651 "SB10RD.f"
	*info = 2;
#line 652 "SB10RD.f"
	return 0;
#line 653 "SB10RD.f"
    }

/*     Compute A1 = inv(Q)*D22 and inv(Q) . */

#line 657 "SB10RD.f"
    dlacpy_("Full", &np2, &m2, &d__[np1 + 1 + (m1 + 1) * d_dim1], ldd, &dwork[
	    iw3], &np2, (ftnlen)4);
#line 659 "SB10RD.f"
    dgetrs_("NoTranspose", &np2, &m2, &dwork[iw2], &np2, &iwork[1], &dwork[
	    iw3], &np2, &info2, (ftnlen)11);
#line 661 "SB10RD.f"
    i__1 = *ldwork - iwrk + 1;
#line 661 "SB10RD.f"
    dgetri_(&np2, &dwork[iw2], &np2, &iwork[1], &dwork[iwrk], &i__1, &info2);
/* Computing MAX */
#line 663 "SB10RD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 663 "SB10RD.f"
    lwamax = max(i__1,lwamax);

/*     Compute A2 = ( inv(Ty) - inv(Q)*inv(Ty) - */
/*                    A1*Tu*D11HAT )*inv(D21HAT) . */

#line 668 "SB10RD.f"
    dlacpy_("Full", &np2, &np2, &ty[ty_offset], ldty, &dwork[iw4], &np2, (
	    ftnlen)4);
#line 669 "SB10RD.f"
    dgetrf_(&np2, &np2, &dwork[iw4], &np2, &iwork[1], &info2);
#line 670 "SB10RD.f"
    i__1 = *ldwork - iwrk + 1;
#line 670 "SB10RD.f"
    dgetri_(&np2, &dwork[iw4], &np2, &iwork[1], &dwork[iwrk], &i__1, &info2);

#line 673 "SB10RD.f"
    dlacpy_("Full", &np2, &np2, &dwork[iw4], &np2, &dwork[iwrk], &np2, (
	    ftnlen)4);
#line 675 "SB10RD.f"
    dgemm_("N", "N", &np2, &np2, &np2, &c_b19, &dwork[iw2], &np2, &dwork[iwrk]
	    , &np2, &c_b8, &dwork[iw4], &np2, (ftnlen)1, (ftnlen)1);
#line 677 "SB10RD.f"
    dgemm_("N", "N", &np2, &np2, &m2, &c_b19, &dwork[iw3], &np2, &dwork[iw1], 
	    &m2, &c_b8, &dwork[iw4], &np2, (ftnlen)1, (ftnlen)1);
#line 679 "SB10RD.f"
    dtrmm_("R", "U", "N", "N", &np2, &np2, &c_b8, &dwork[id21], &np2, &dwork[
	    iw4], &np2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Compute [ A1  A2 ]*CHAT . */

#line 684 "SB10RD.f"
    i__1 = m2 + np2;
#line 684 "SB10RD.f"
    i__2 = m2 + np2;
#line 684 "SB10RD.f"
    dgemm_("N", "N", &np2, n, &i__1, &c_b8, &dwork[iw3], &np2, &dwork[iwc], &
	    i__2, &c_b7, &dwork[iwrk], &np2, (ftnlen)1, (ftnlen)1);

/*     Compute AK := AHAT - BHAT1*[ A1  A2 ]*CHAT . */

#line 689 "SB10RD.f"
    dgemm_("N", "N", n, n, &np2, &c_b19, &dwork[iwb], n, &dwork[iwrk], &np2, &
	    c_b8, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);

/*     Compute BK := BHAT1*inv(Q) . */

#line 694 "SB10RD.f"
    dgemm_("N", "N", n, &np2, &np2, &c_b8, &dwork[iwb], n, &dwork[iw2], &np2, 
	    &c_b7, &bk[bk_offset], ldbk, (ftnlen)1, (ftnlen)1);

/*     Compute DK := Tu*D11HAT*Ty*inv(Q) . */

#line 699 "SB10RD.f"
    dgemm_("N", "N", &m2, &np2, &np2, &c_b8, &dk[dk_offset], lddk, &dwork[iw2]
	    , &np2, &c_b7, &dwork[iw3], &m2, (ftnlen)1, (ftnlen)1);
#line 701 "SB10RD.f"
    dlacpy_("Full", &m2, &np2, &dwork[iw3], &m2, &dk[dk_offset], lddk, (
	    ftnlen)4);

#line 703 "SB10RD.f"
    dwork[1] = (doublereal) lwamax;
#line 704 "SB10RD.f"
    return 0;
/* *** Last line of SB10RD *** */
} /* sb10rd_ */

