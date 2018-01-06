#line 1 "SB10HD.f"
/* SB10HD.f -- translated by f2c (version 20100827).
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

#line 1 "SB10HD.f"
/* Subroutine */ int sb10hd_(integer *n, integer *m, integer *np, integer *
	ncon, integer *nmeas, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *d__, integer 
	*ldd, doublereal *ak, integer *ldak, doublereal *bk, integer *ldbk, 
	doublereal *ck, integer *ldck, doublereal *dk, integer *lddk, 
	doublereal *rcond, doublereal *tol, integer *iwork, doublereal *dwork,
	 integer *ldwork, logical *bwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ak_dim1, ak_offset, b_dim1, b_offset, bk_dim1, 
	    bk_offset, c_dim1, c_offset, ck_dim1, ck_offset, d_dim1, d_offset,
	     dk_dim1, dk_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, 
	    i__8;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer m1, m2, np1, np2, iwc, iwd, iwf, iwh, iwy;
    static doublereal toll;
    static integer iwrk, iwtu, iwty, info2;
    extern /* Subroutine */ int sb10ud_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *), sb10vd_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    logical *, integer *), sb10wd_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *);
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

/*     To compute the matrices of the H2 optimal n-state controller */

/*              | AK | BK | */
/*          K = |----|----| */
/*              | CK | DK | */

/*     for the system */

/*                   | A  | B1  B2  |   | A | B | */
/*               P = |----|---------| = |---|---| , */
/*                   | C1 |  0  D12 |   | C | D | */
/*                   | C2 | D21 D22 | */

/*     where B2 has as column size the number of control inputs (NCON) */
/*     and C2 has as row size the number of measurements (NMEAS) being */
/*     provided to the controller. */

/*     It is assumed that */

/*     (A1) (A,B2) is stabilizable and (C2,A) is detectable, */

/*     (A2) The block D11 of D is zero, */

/*     (A3) D12 is full column rank and D21 is full row rank. */

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
/*             SLICOT Library routine SB10UD. Transformation matrices */
/*             whose reciprocal condition numbers are less than TOL are */
/*             not allowed. If TOL <= 0, then a default value equal to */
/*             sqrt(EPS) is used, where EPS is the relative machine */
/*             precision. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension max(2*N,N*N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal */
/*             LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= N*M + NP*(N+M) + M2*M2 + NP2*NP2 + */
/*                       max(max(M2 + NP1*NP1 + */
/*                               max(NP1*N,3*M2+NP1,5*M2), */
/*                               NP2 + M1*M1 + */
/*                               max(M1*N,3*NP2+M1,5*NP2), */
/*                               N*M2,NP2*N,NP2*M2,1), */
/*                               N*(14*N+12+M2+NP2)+5), */
/*             where M1 = M - M2 and NP1 = NP - NP2. */
/*             For good performance, LDWORK must generally be larger. */
/*             Denoting Q = max(M1,M2,NP1,NP2), an upper bound is */
/*             2*Q*(3*Q+2*N)+max(1,Q*(Q+max(N,5)+1),N*(14*N+12+2*Q)+5). */

/*     BWORK   LOGICAL array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the matrix D12 had not full column rank in */
/*                   respect to the tolerance TOL; */
/*             = 2:  if the matrix D21 had not full row rank in respect */
/*                   to the tolerance TOL; */
/*             = 3:  if the singular value decomposition (SVD) algorithm */
/*                   did not converge (when computing the SVD of one of */
/*                   the matrices D12 or D21). */
/*             = 4:  if the X-Riccati equation was not solved */
/*                   successfully; */
/*             = 5:  if the Y-Riccati equation was not solved */
/*                   successfully. */

/*     METHOD */

/*     The routine implements the formulas given in [1], [2]. */

/*     REFERENCES */

/*     [1] Zhou, K., Doyle, J.C., and Glover, K. */
/*         Robust and Optimal Control. */
/*         Prentice-Hall, Upper Saddle River, NJ, 1996. */

/*     [2] Balas, G.J., Doyle, J.C., Glover, K., Packard, A., and */
/*         Smith, R. */
/*         mu-Analysis and Synthesis Toolbox. */
/*         The MathWorks Inc., Natick, Mass., 1995. */

/*     NUMERICAL ASPECTS */

/*     The accuracy of the result depends on the condition numbers of the */
/*     input and output transformations and on the condition numbers of */
/*     the two Riccati equations, as given by the values of RCOND(1), */
/*     RCOND(2), RCOND(3) and RCOND(4), respectively. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, Oct. 1998. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 1999, */
/*     Sept. 1999, Jan. 2000, Feb. 2000. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, H2 optimal control, optimal regulator, */
/*     robust control. */

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

#line 266 "SB10HD.f"
    /* Parameter adjustments */
#line 266 "SB10HD.f"
    a_dim1 = *lda;
#line 266 "SB10HD.f"
    a_offset = 1 + a_dim1;
#line 266 "SB10HD.f"
    a -= a_offset;
#line 266 "SB10HD.f"
    b_dim1 = *ldb;
#line 266 "SB10HD.f"
    b_offset = 1 + b_dim1;
#line 266 "SB10HD.f"
    b -= b_offset;
#line 266 "SB10HD.f"
    c_dim1 = *ldc;
#line 266 "SB10HD.f"
    c_offset = 1 + c_dim1;
#line 266 "SB10HD.f"
    c__ -= c_offset;
#line 266 "SB10HD.f"
    d_dim1 = *ldd;
#line 266 "SB10HD.f"
    d_offset = 1 + d_dim1;
#line 266 "SB10HD.f"
    d__ -= d_offset;
#line 266 "SB10HD.f"
    ak_dim1 = *ldak;
#line 266 "SB10HD.f"
    ak_offset = 1 + ak_dim1;
#line 266 "SB10HD.f"
    ak -= ak_offset;
#line 266 "SB10HD.f"
    bk_dim1 = *ldbk;
#line 266 "SB10HD.f"
    bk_offset = 1 + bk_dim1;
#line 266 "SB10HD.f"
    bk -= bk_offset;
#line 266 "SB10HD.f"
    ck_dim1 = *ldck;
#line 266 "SB10HD.f"
    ck_offset = 1 + ck_dim1;
#line 266 "SB10HD.f"
    ck -= ck_offset;
#line 266 "SB10HD.f"
    dk_dim1 = *lddk;
#line 266 "SB10HD.f"
    dk_offset = 1 + dk_dim1;
#line 266 "SB10HD.f"
    dk -= dk_offset;
#line 266 "SB10HD.f"
    --rcond;
#line 266 "SB10HD.f"
    --iwork;
#line 266 "SB10HD.f"
    --dwork;
#line 266 "SB10HD.f"
    --bwork;
#line 266 "SB10HD.f"

#line 266 "SB10HD.f"
    /* Function Body */
#line 266 "SB10HD.f"
    m1 = *m - *ncon;
#line 267 "SB10HD.f"
    m2 = *ncon;
#line 268 "SB10HD.f"
    np1 = *np - *nmeas;
#line 269 "SB10HD.f"
    np2 = *nmeas;

#line 271 "SB10HD.f"
    *info = 0;
#line 272 "SB10HD.f"
    if (*n < 0) {
#line 273 "SB10HD.f"
	*info = -1;
#line 274 "SB10HD.f"
    } else if (*m < 0) {
#line 275 "SB10HD.f"
	*info = -2;
#line 276 "SB10HD.f"
    } else if (*np < 0) {
#line 277 "SB10HD.f"
	*info = -3;
#line 278 "SB10HD.f"
    } else if (*ncon < 0 || m1 < 0 || m2 > np1) {
#line 279 "SB10HD.f"
	*info = -4;
#line 280 "SB10HD.f"
    } else if (*nmeas < 0 || np1 < 0 || np2 > m1) {
#line 281 "SB10HD.f"
	*info = -5;
#line 282 "SB10HD.f"
    } else if (*lda < max(1,*n)) {
#line 283 "SB10HD.f"
	*info = -7;
#line 284 "SB10HD.f"
    } else if (*ldb < max(1,*n)) {
#line 285 "SB10HD.f"
	*info = -9;
#line 286 "SB10HD.f"
    } else if (*ldc < max(1,*np)) {
#line 287 "SB10HD.f"
	*info = -11;
#line 288 "SB10HD.f"
    } else if (*ldd < max(1,*np)) {
#line 289 "SB10HD.f"
	*info = -13;
#line 290 "SB10HD.f"
    } else if (*ldak < max(1,*n)) {
#line 291 "SB10HD.f"
	*info = -15;
#line 292 "SB10HD.f"
    } else if (*ldbk < max(1,*n)) {
#line 293 "SB10HD.f"
	*info = -17;
#line 294 "SB10HD.f"
    } else if (*ldck < max(1,m2)) {
#line 295 "SB10HD.f"
	*info = -19;
#line 296 "SB10HD.f"
    } else if (*lddk < max(1,m2)) {
#line 297 "SB10HD.f"
	*info = -21;
#line 298 "SB10HD.f"
    } else {

/*        Compute workspace. */

/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
#line 302 "SB10HD.f"
	i__5 = np1 * *n, i__6 = m2 * 3 + np1, i__5 = max(i__5,i__6), i__6 = 
		m2 * 5;
/* Computing MAX */
#line 302 "SB10HD.f"
	i__7 = m1 * *n, i__8 = np2 * 3 + m1, i__7 = max(i__7,i__8), i__8 = 
		np2 * 5;
#line 302 "SB10HD.f"
	i__3 = m2 + np1 * np1 + max(i__5,i__6), i__4 = np2 + m1 * m1 + max(
		i__7,i__8), i__3 = max(i__3,i__4), i__4 = *n * m2, i__3 = max(
		i__3,i__4), i__4 = np2 * *n, i__3 = max(i__3,i__4), i__4 = 
		np2 * m2, i__3 = max(i__3,i__4);
#line 302 "SB10HD.f"
	i__1 = max(i__3,1), i__2 = *n * (*n * 14 + 12 + m2 + np2) + 5;
#line 302 "SB10HD.f"
	minwrk = *n * *m + *np * (*n + *m) + m2 * m2 + np2 * np2 + max(i__1,
		i__2);
#line 309 "SB10HD.f"
	if (*ldwork < minwrk) {
#line 309 "SB10HD.f"
	    *info = -26;
#line 309 "SB10HD.f"
	}
#line 311 "SB10HD.f"
    }
#line 312 "SB10HD.f"
    if (*info != 0) {
#line 313 "SB10HD.f"
	i__1 = -(*info);
#line 313 "SB10HD.f"
	xerbla_("SB10HD", &i__1, (ftnlen)6);
#line 314 "SB10HD.f"
	return 0;
#line 315 "SB10HD.f"
    }

/*     Quick return if possible. */

#line 319 "SB10HD.f"
    if (*n == 0 || *m == 0 || *np == 0 || m1 == 0 || m2 == 0 || np1 == 0 || 
	    np2 == 0) {
#line 321 "SB10HD.f"
	rcond[1] = 1.;
#line 322 "SB10HD.f"
	rcond[2] = 1.;
#line 323 "SB10HD.f"
	rcond[3] = 1.;
#line 324 "SB10HD.f"
	rcond[4] = 1.;
#line 325 "SB10HD.f"
	dwork[1] = 1.;
#line 326 "SB10HD.f"
	return 0;
#line 327 "SB10HD.f"
    }

#line 329 "SB10HD.f"
    toll = *tol;
#line 330 "SB10HD.f"
    if (toll <= 0.) {

/*        Set the default value of the tolerance for rank tests. */

#line 334 "SB10HD.f"
	toll = sqrt(dlamch_("Epsilon", (ftnlen)7));
#line 335 "SB10HD.f"
    }

/*     Workspace usage. */

#line 339 "SB10HD.f"
    iwc = *n * *m + 1;
#line 340 "SB10HD.f"
    iwd = iwc + *np * *n;
#line 341 "SB10HD.f"
    iwtu = iwd + *np * *m;
#line 342 "SB10HD.f"
    iwty = iwtu + m2 * m2;
#line 343 "SB10HD.f"
    iwrk = iwty + np2 * np2;

#line 345 "SB10HD.f"
    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[1], n, (ftnlen)4);
#line 346 "SB10HD.f"
    dlacpy_("Full", np, n, &c__[c_offset], ldc, &dwork[iwc], np, (ftnlen)4);
#line 347 "SB10HD.f"
    dlacpy_("Full", np, m, &d__[d_offset], ldd, &dwork[iwd], np, (ftnlen)4);

/*     Transform the system so that D12 and D21 satisfy the formulas */
/*     in the computation of the H2 optimal controller. */

#line 352 "SB10HD.f"
    i__1 = *ldwork - iwrk + 1;
#line 352 "SB10HD.f"
    sb10ud_(n, m, np, ncon, nmeas, &dwork[1], n, &dwork[iwc], np, &dwork[iwd],
	     np, &dwork[iwtu], &m2, &dwork[iwty], &np2, &rcond[1], &toll, &
	    dwork[iwrk], &i__1, &info2);
#line 356 "SB10HD.f"
    if (info2 > 0) {
#line 357 "SB10HD.f"
	*info = info2;
#line 358 "SB10HD.f"
	return 0;
#line 359 "SB10HD.f"
    }
#line 360 "SB10HD.f"
    lwamax = (integer) dwork[iwrk] + iwrk - 1;

#line 362 "SB10HD.f"
    iwy = iwrk;
#line 363 "SB10HD.f"
    iwf = iwy + *n * *n;
#line 364 "SB10HD.f"
    iwh = iwf + m2 * *n;
#line 365 "SB10HD.f"
    iwrk = iwh + *n * np2;

/*     Compute the optimal state feedback and output injection matrices. */
/*     AK is used to store X. */

#line 370 "SB10HD.f"
    i__1 = *ldwork - iwrk + 1;
#line 370 "SB10HD.f"
    sb10vd_(n, m, np, ncon, nmeas, &a[a_offset], lda, &dwork[1], n, &dwork[
	    iwc], np, &dwork[iwf], &m2, &dwork[iwh], n, &ak[ak_offset], ldak, 
	    &dwork[iwy], n, &rcond[3], &iwork[1], &dwork[iwrk], &i__1, &bwork[
	    1], &info2);
#line 374 "SB10HD.f"
    if (info2 > 0) {
#line 375 "SB10HD.f"
	*info = info2 + 3;
#line 376 "SB10HD.f"
	return 0;
#line 377 "SB10HD.f"
    }
/* Computing MAX */
#line 378 "SB10HD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 378 "SB10HD.f"
    lwamax = max(i__1,lwamax);

/*     Compute the H2 optimal controller. */

#line 382 "SB10HD.f"
    sb10wd_(n, m, np, ncon, nmeas, &a[a_offset], lda, &dwork[1], n, &dwork[
	    iwc], np, &dwork[iwd], np, &dwork[iwf], &m2, &dwork[iwh], n, &
	    dwork[iwtu], &m2, &dwork[iwty], &np2, &ak[ak_offset], ldak, &bk[
	    bk_offset], ldbk, &ck[ck_offset], ldck, &dk[dk_offset], lddk, &
	    info2);

#line 387 "SB10HD.f"
    dwork[1] = (doublereal) lwamax;
#line 388 "SB10HD.f"
    return 0;
/* *** Last line of SB10HD *** */
} /* sb10hd_ */

