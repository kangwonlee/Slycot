#line 1 "SB10PD.f"
/* SB10PD.f -- translated by f2c (version 20100827).
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

#line 1 "SB10PD.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b27 = 1.;
static doublereal c_b28 = 0.;

/* Subroutine */ int sb10pd_(integer *n, integer *m, integer *np, integer *
	ncon, integer *nmeas, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *d__, integer 
	*ldd, doublereal *tu, integer *ldtu, doublereal *ty, integer *ldty, 
	doublereal *rcond, doublereal *tol, doublereal *dwork, integer *
	ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, tu_dim1, tu_offset, ty_dim1, ty_offset, i__1, i__2, 
	    i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, m1, m2, iq, nd1, nd2, np1, np2;
    static doublereal eps;
    static integer iext;
    static doublereal toll;
    static integer iwrk, info2;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen), dswap_(
	    integer *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dgesvd_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dlacpy_(char *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
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

/*     To reduce the matrices D12 and D21 of the linear time-invariant */
/*     system */

/*                   | A  | B1  B2  |   | A | B | */
/*               P = |----|---------| = |---|---| */
/*                   | C1 | D11 D12 |   | C | D | */
/*                   | C2 | D21 D22 | */

/*     to unit diagonal form, to transform the matrices B, C, and D11 to */
/*     satisfy the formulas in the computation of an H2 and H-infinity */
/*     (sub)optimal controllers and to check the rank conditions. */

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

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the system input matrix B. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the transformed system input matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading NP-by-N part of this array must */
/*             contain the system output matrix C. */
/*             On exit, the leading NP-by-N part of this array contains */
/*             the transformed system output matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= max(1,NP). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading NP-by-M part of this array must */
/*             contain the system input/output matrix D. The */
/*             NMEAS-by-NCON trailing submatrix D22 is not referenced. */
/*             On exit, the leading (NP-NMEAS)-by-(M-NCON) part of this */
/*             array contains the transformed submatrix D11. */
/*             The transformed submatrices D12 = [ 0  Im2 ]' and */
/*             D21 = [ 0  Inp2 ] are not stored. The corresponding part */
/*             of this array contains no useful information. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= max(1,NP). */

/*     TU      (output) DOUBLE PRECISION array, dimension (LDTU,M2) */
/*             The leading M2-by-M2 part of this array contains the */
/*             control transformation matrix TU. */

/*     LDTU    INTEGER */
/*             The leading dimension of the array TU.  LDTU >= max(1,M2). */

/*     TY      (output) DOUBLE PRECISION array, dimension (LDTY,NP2) */
/*             The leading NP2-by-NP2 part of this array contains the */
/*             measurement transformation matrix TY. */

/*     LDTY    INTEGER */
/*             The leading dimension of the array TY. */
/*             LDTY >= max(1,NP2). */

/*     RCOND   (output) DOUBLE PRECISION array, dimension (2) */
/*             RCOND(1) contains the reciprocal condition number of the */
/*                      control transformation matrix TU; */
/*             RCOND(2) contains the reciprocal condition number of the */
/*                      measurement transformation matrix TY. */
/*             RCOND is set even if INFO = 3 or INFO = 4; if INFO = 3, */
/*             then RCOND(2) was not computed, but it is set to 0. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             Tolerance used for controlling the accuracy of the applied */
/*             transformations. Transformation matrices TU and TY whose */
/*             reciprocal condition numbers are less than TOL are not */
/*             allowed. If TOL <= 0, then a default value equal to */
/*             sqrt(EPS) is used, where EPS is the relative machine */
/*             precision. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal */
/*             LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= MAX(1,LW1,LW2,LW3,LW4), where */
/*             LW1 = (N+NP1+1)*(N+M2) + MAX(3*(N+M2)+N+NP1,5*(N+M2)), */
/*             LW2 = (N+NP2)*(N+M1+1) + MAX(3*(N+NP2)+N+M1,5*(N+NP2)), */
/*             LW3 = M2 + NP1*NP1 + MAX(NP1*MAX(N,M1),3*M2+NP1,5*M2), */
/*             LW4 = NP2 + M1*M1 + MAX(MAX(N,NP1)*M1,3*NP2+M1,5*NP2), */
/*             with M1 = M - M2 and NP1 = NP - NP2. */
/*             For good performance, LDWORK must generally be larger. */
/*             Denoting Q = MAX(M1,M2,NP1,NP2), an upper bound is */
/*             MAX(1,(N+Q)*(N+Q+6),Q*(Q+MAX(N,Q,5)+1). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the matrix | A   B2  | had not full column rank */
/*                                 | C1  D12 | */
/*                   in respect to the tolerance EPS; */
/*             = 2:  if the matrix | A   B1  | had not full row rank in */
/*                                 | C2  D21 | */
/*                   respect to the tolerance EPS; */
/*             = 3:  if the matrix D12 had not full column rank in */
/*                   respect to the tolerance TOL; */
/*             = 4:  if the matrix D21 had not full row rank in respect */
/*                   to the tolerance TOL; */
/*             = 5:  if the singular value decomposition (SVD) algorithm */
/*                   did not converge (when computing the SVD of one of */
/*                   the matrices |A   B2 |, |A   B1 |, D12 or D21). */
/*                                |C1  D12|  |C2  D21| */

/*     METHOD */

/*     The routine performs the transformations described in [2]. */

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

/*     The precision of the transformations can be controlled by the */
/*     condition numbers of the matrices TU and TY as given by the */
/*     values of RCOND(1) and RCOND(2), respectively. An error return */
/*     with INFO = 3 or INFO = 4 will be obtained if the condition */
/*     number of TU or TY, respectively, would exceed 1/TOL. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 1998. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 1999, */
/*     Feb. 2000. */

/*     KEYWORDS */

/*     H-infinity optimal control, robust control, singular value */
/*     decomposition. */

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

#line 241 "SB10PD.f"
    /* Parameter adjustments */
#line 241 "SB10PD.f"
    a_dim1 = *lda;
#line 241 "SB10PD.f"
    a_offset = 1 + a_dim1;
#line 241 "SB10PD.f"
    a -= a_offset;
#line 241 "SB10PD.f"
    b_dim1 = *ldb;
#line 241 "SB10PD.f"
    b_offset = 1 + b_dim1;
#line 241 "SB10PD.f"
    b -= b_offset;
#line 241 "SB10PD.f"
    c_dim1 = *ldc;
#line 241 "SB10PD.f"
    c_offset = 1 + c_dim1;
#line 241 "SB10PD.f"
    c__ -= c_offset;
#line 241 "SB10PD.f"
    d_dim1 = *ldd;
#line 241 "SB10PD.f"
    d_offset = 1 + d_dim1;
#line 241 "SB10PD.f"
    d__ -= d_offset;
#line 241 "SB10PD.f"
    tu_dim1 = *ldtu;
#line 241 "SB10PD.f"
    tu_offset = 1 + tu_dim1;
#line 241 "SB10PD.f"
    tu -= tu_offset;
#line 241 "SB10PD.f"
    ty_dim1 = *ldty;
#line 241 "SB10PD.f"
    ty_offset = 1 + ty_dim1;
#line 241 "SB10PD.f"
    ty -= ty_offset;
#line 241 "SB10PD.f"
    --rcond;
#line 241 "SB10PD.f"
    --dwork;
#line 241 "SB10PD.f"

#line 241 "SB10PD.f"
    /* Function Body */
#line 241 "SB10PD.f"
    m1 = *m - *ncon;
#line 242 "SB10PD.f"
    m2 = *ncon;
#line 243 "SB10PD.f"
    np1 = *np - *nmeas;
#line 244 "SB10PD.f"
    np2 = *nmeas;

#line 246 "SB10PD.f"
    *info = 0;
#line 247 "SB10PD.f"
    if (*n < 0) {
#line 248 "SB10PD.f"
	*info = -1;
#line 249 "SB10PD.f"
    } else if (*m < 0) {
#line 250 "SB10PD.f"
	*info = -2;
#line 251 "SB10PD.f"
    } else if (*np < 0) {
#line 252 "SB10PD.f"
	*info = -3;
#line 253 "SB10PD.f"
    } else if (*ncon < 0 || m1 < 0 || m2 > np1) {
#line 254 "SB10PD.f"
	*info = -4;
#line 255 "SB10PD.f"
    } else if (*nmeas < 0 || np1 < 0 || np2 > m1) {
#line 256 "SB10PD.f"
	*info = -5;
#line 257 "SB10PD.f"
    } else if (*lda < max(1,*n)) {
#line 258 "SB10PD.f"
	*info = -7;
#line 259 "SB10PD.f"
    } else if (*ldb < max(1,*n)) {
#line 260 "SB10PD.f"
	*info = -9;
#line 261 "SB10PD.f"
    } else if (*ldc < max(1,*np)) {
#line 262 "SB10PD.f"
	*info = -11;
#line 263 "SB10PD.f"
    } else if (*ldd < max(1,*np)) {
#line 264 "SB10PD.f"
	*info = -13;
#line 265 "SB10PD.f"
    } else if (*ldtu < max(1,m2)) {
#line 266 "SB10PD.f"
	*info = -15;
#line 267 "SB10PD.f"
    } else if (*ldty < max(1,np2)) {
#line 268 "SB10PD.f"
	*info = -17;
#line 269 "SB10PD.f"
    } else {

/*        Compute workspace. */

/* Computing MAX */
/* Computing MAX */
#line 273 "SB10PD.f"
	i__3 = (*n + m2) * 3 + *n + np1, i__4 = (*n + m2) * 5;
/* Computing MAX */
#line 273 "SB10PD.f"
	i__5 = (*n + np2) * 3 + *n + m1, i__6 = (*n + np2) * 5;
/* Computing MAX */
#line 273 "SB10PD.f"
	i__7 = np1 * max(*n,m1), i__8 = m2 * 3 + np1, i__7 = max(i__7,i__8), 
		i__8 = m2 * 5;
/* Computing MAX */
#line 273 "SB10PD.f"
	i__9 = max(*n,np1) * m1, i__10 = np2 * 3 + m1, i__9 = max(i__9,i__10),
		 i__10 = np2 * 5;
#line 273 "SB10PD.f"
	i__1 = 1, i__2 = (*n + np1 + 1) * (*n + m2) + max(i__3,i__4), i__1 = 
		max(i__1,i__2), i__2 = (*n + np2) * (*n + m1 + 1) + max(i__5,
		i__6), i__1 = max(i__1,i__2), i__2 = m2 + np1 * np1 + max(
		i__7,i__8), i__1 = max(i__1,i__2), i__2 = np2 + m1 * m1 + max(
		i__9,i__10);
#line 273 "SB10PD.f"
	minwrk = max(i__1,i__2);
#line 282 "SB10PD.f"
	if (*ldwork < minwrk) {
#line 282 "SB10PD.f"
	    *info = -21;
#line 282 "SB10PD.f"
	}
#line 284 "SB10PD.f"
    }
#line 285 "SB10PD.f"
    if (*info != 0) {
#line 286 "SB10PD.f"
	i__1 = -(*info);
#line 286 "SB10PD.f"
	xerbla_("SB10PD", &i__1, (ftnlen)6);
#line 287 "SB10PD.f"
	return 0;
#line 288 "SB10PD.f"
    }

/*     Quick return if possible. */

#line 292 "SB10PD.f"
    if (*n == 0 || *m == 0 || *np == 0 || m1 == 0 || m2 == 0 || np1 == 0 || 
	    np2 == 0) {
#line 294 "SB10PD.f"
	rcond[1] = 1.;
#line 295 "SB10PD.f"
	rcond[2] = 1.;
#line 296 "SB10PD.f"
	dwork[1] = 1.;
#line 297 "SB10PD.f"
	return 0;
#line 298 "SB10PD.f"
    }

#line 300 "SB10PD.f"
    nd1 = np1 - m2;
#line 301 "SB10PD.f"
    nd2 = m1 - np2;
#line 302 "SB10PD.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 303 "SB10PD.f"
    toll = *tol;
#line 304 "SB10PD.f"
    if (toll <= 0.) {

/*        Set the default value of the tolerance for condition tests. */

#line 308 "SB10PD.f"
	toll = sqrt(eps);
#line 309 "SB10PD.f"
    }

/*     Determine if |A-jwI  B2 | has full column rank at w = 0. */
/*                  |  C1   D12| */
/*     Workspace:  need   (N+NP1+1)*(N+M2) + */
/*                        max(3*(N+M2)+N+NP1,5*(N+M2)); */
/*                 prefer larger. */

#line 317 "SB10PD.f"
    iext = *n + m2 + 1;
#line 318 "SB10PD.f"
    iwrk = iext + (*n + np1) * (*n + m2);
#line 319 "SB10PD.f"
    i__1 = *n + np1;
#line 319 "SB10PD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[iext], &i__1, (ftnlen)4);
#line 320 "SB10PD.f"
    i__1 = *n + np1;
#line 320 "SB10PD.f"
    dlacpy_("Full", &np1, n, &c__[c_offset], ldc, &dwork[iext + *n], &i__1, (
	    ftnlen)4);
#line 321 "SB10PD.f"
    i__1 = *n + np1;
#line 321 "SB10PD.f"
    dlacpy_("Full", n, &m2, &b[(m1 + 1) * b_dim1 + 1], ldb, &dwork[iext + (*n 
	    + np1) * *n], &i__1, (ftnlen)4);
#line 323 "SB10PD.f"
    i__1 = *n + np1;
#line 323 "SB10PD.f"
    dlacpy_("Full", &np1, &m2, &d__[(m1 + 1) * d_dim1 + 1], ldd, &dwork[iext 
	    + (*n + np1) * *n + *n], &i__1, (ftnlen)4);
#line 325 "SB10PD.f"
    i__1 = *n + np1;
#line 325 "SB10PD.f"
    i__2 = *n + m2;
#line 325 "SB10PD.f"
    i__3 = *n + np1;
#line 325 "SB10PD.f"
    i__4 = *ldwork - iwrk + 1;
#line 325 "SB10PD.f"
    dgesvd_("N", "N", &i__1, &i__2, &dwork[iext], &i__3, &dwork[1], &tu[
	    tu_offset], ldtu, &ty[ty_offset], ldty, &dwork[iwrk], &i__4, &
	    info2, (ftnlen)1, (ftnlen)1);
#line 328 "SB10PD.f"
    if (info2 != 0) {
#line 329 "SB10PD.f"
	*info = 5;
#line 330 "SB10PD.f"
	return 0;
#line 331 "SB10PD.f"
    }
#line 332 "SB10PD.f"
    if (dwork[*n + m2] / dwork[1] <= eps) {
#line 333 "SB10PD.f"
	*info = 1;
#line 334 "SB10PD.f"
	return 0;
#line 335 "SB10PD.f"
    }
#line 336 "SB10PD.f"
    lwamax = (integer) dwork[iwrk] + iwrk - 1;

/*     Determine if |A-jwI  B1 | has full row rank at w = 0. */
/*                  |  C2   D21| */
/*     Workspace:  need   (N+NP2)*(N+M1+1) + */
/*                        max(3*(N+NP2)+N+M1,5*(N+NP2)); */
/*                 prefer larger. */

#line 344 "SB10PD.f"
    iext = *n + np2 + 1;
#line 345 "SB10PD.f"
    iwrk = iext + (*n + np2) * (*n + m1);
#line 346 "SB10PD.f"
    i__1 = *n + np2;
#line 346 "SB10PD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[iext], &i__1, (ftnlen)4);
#line 347 "SB10PD.f"
    i__1 = *n + np2;
#line 347 "SB10PD.f"
    dlacpy_("Full", &np2, n, &c__[np1 + 1 + c_dim1], ldc, &dwork[iext + *n], &
	    i__1, (ftnlen)4);
#line 349 "SB10PD.f"
    i__1 = *n + np2;
#line 349 "SB10PD.f"
    dlacpy_("Full", n, &m1, &b[b_offset], ldb, &dwork[iext + (*n + np2) * *n],
	     &i__1, (ftnlen)4);
#line 351 "SB10PD.f"
    i__1 = *n + np2;
#line 351 "SB10PD.f"
    dlacpy_("Full", &np2, &m1, &d__[np1 + 1 + d_dim1], ldd, &dwork[iext + (*n 
	    + np2) * *n + *n], &i__1, (ftnlen)4);
#line 353 "SB10PD.f"
    i__1 = *n + np2;
#line 353 "SB10PD.f"
    i__2 = *n + m1;
#line 353 "SB10PD.f"
    i__3 = *n + np2;
#line 353 "SB10PD.f"
    i__4 = *ldwork - iwrk + 1;
#line 353 "SB10PD.f"
    dgesvd_("N", "N", &i__1, &i__2, &dwork[iext], &i__3, &dwork[1], &tu[
	    tu_offset], ldtu, &ty[ty_offset], ldty, &dwork[iwrk], &i__4, &
	    info2, (ftnlen)1, (ftnlen)1);
#line 356 "SB10PD.f"
    if (info2 != 0) {
#line 357 "SB10PD.f"
	*info = 5;
#line 358 "SB10PD.f"
	return 0;
#line 359 "SB10PD.f"
    }
#line 360 "SB10PD.f"
    if (dwork[*n + np2] / dwork[1] <= eps) {
#line 361 "SB10PD.f"
	*info = 2;
#line 362 "SB10PD.f"
	return 0;
#line 363 "SB10PD.f"
    }
/* Computing MAX */
#line 364 "SB10PD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 364 "SB10PD.f"
    lwamax = max(i__1,lwamax);

/*     Determine SVD of D12, D12 = U12 S12 V12', and check if D12 has */
/*     full column rank. V12' is stored in TU. */
/*     Workspace:  need   M2 + NP1*NP1 + max(3*M2+NP1,5*M2); */
/*                 prefer larger. */

#line 371 "SB10PD.f"
    iq = m2 + 1;
#line 372 "SB10PD.f"
    iwrk = iq + np1 * np1;

#line 374 "SB10PD.f"
    i__1 = *ldwork - iwrk + 1;
#line 374 "SB10PD.f"
    dgesvd_("A", "A", &np1, &m2, &d__[(m1 + 1) * d_dim1 + 1], ldd, &dwork[1], 
	    &dwork[iq], &np1, &tu[tu_offset], ldtu, &dwork[iwrk], &i__1, &
	    info2, (ftnlen)1, (ftnlen)1);
#line 377 "SB10PD.f"
    if (info2 != 0) {
#line 378 "SB10PD.f"
	*info = 5;
#line 379 "SB10PD.f"
	return 0;
#line 380 "SB10PD.f"
    }

#line 382 "SB10PD.f"
    rcond[1] = dwork[m2] / dwork[1];
#line 383 "SB10PD.f"
    if (rcond[1] <= toll) {
#line 384 "SB10PD.f"
	rcond[2] = 0.;
#line 385 "SB10PD.f"
	*info = 3;
#line 386 "SB10PD.f"
	return 0;
#line 387 "SB10PD.f"
    }
/* Computing MAX */
#line 388 "SB10PD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 388 "SB10PD.f"
    lwamax = max(i__1,lwamax);

/*     Determine Q12. */

#line 392 "SB10PD.f"
    if (nd1 > 0) {
#line 393 "SB10PD.f"
	dlacpy_("Full", &np1, &m2, &dwork[iq], &np1, &d__[(m1 + 1) * d_dim1 + 
		1], ldd, (ftnlen)4);
#line 395 "SB10PD.f"
	dlacpy_("Full", &np1, &nd1, &dwork[iq + np1 * m2], &np1, &dwork[iq], &
		np1, (ftnlen)4);
#line 397 "SB10PD.f"
	dlacpy_("Full", &np1, &m2, &d__[(m1 + 1) * d_dim1 + 1], ldd, &dwork[
		iq + np1 * nd1], &np1, (ftnlen)4);
#line 399 "SB10PD.f"
    }

/*     Determine Tu by transposing in-situ and scaling. */

#line 403 "SB10PD.f"
    i__1 = m2 - 1;
#line 403 "SB10PD.f"
    for (j = 1; j <= i__1; ++j) {
#line 404 "SB10PD.f"
	dswap_(&j, &tu[j + 1 + tu_dim1], ldtu, &tu[(j + 1) * tu_dim1 + 1], &
		c__1);
#line 405 "SB10PD.f"
/* L10: */
#line 405 "SB10PD.f"
    }

#line 407 "SB10PD.f"
    i__1 = m2;
#line 407 "SB10PD.f"
    for (j = 1; j <= i__1; ++j) {
#line 408 "SB10PD.f"
	d__1 = 1. / dwork[j];
#line 408 "SB10PD.f"
	dscal_(&m2, &d__1, &tu[j * tu_dim1 + 1], &c__1);
#line 409 "SB10PD.f"
/* L20: */
#line 409 "SB10PD.f"
    }

/*     Determine C1 =: Q12'*C1. */
/*     Workspace:  M2 + NP1*NP1 + NP1*N. */

#line 414 "SB10PD.f"
    dgemm_("T", "N", &np1, n, &np1, &c_b27, &dwork[iq], &np1, &c__[c_offset], 
	    ldc, &c_b28, &dwork[iwrk], &np1, (ftnlen)1, (ftnlen)1);
#line 416 "SB10PD.f"
    dlacpy_("Full", &np1, n, &dwork[iwrk], &np1, &c__[c_offset], ldc, (ftnlen)
	    4);
/* Computing MAX */
#line 417 "SB10PD.f"
    i__1 = iwrk + np1 * *n - 1;
#line 417 "SB10PD.f"
    lwamax = max(i__1,lwamax);

/*     Determine D11 =: Q12'*D11. */
/*     Workspace:  M2 + NP1*NP1 + NP1*M1. */

#line 422 "SB10PD.f"
    dgemm_("T", "N", &np1, &m1, &np1, &c_b27, &dwork[iq], &np1, &d__[d_offset]
	    , ldd, &c_b28, &dwork[iwrk], &np1, (ftnlen)1, (ftnlen)1);
#line 424 "SB10PD.f"
    dlacpy_("Full", &np1, &m1, &dwork[iwrk], &np1, &d__[d_offset], ldd, (
	    ftnlen)4);
/* Computing MAX */
#line 425 "SB10PD.f"
    i__1 = iwrk + np1 * m1 - 1;
#line 425 "SB10PD.f"
    lwamax = max(i__1,lwamax);

/*     Determine SVD of D21, D21 = U21 S21 V21', and check if D21 has */
/*     full row rank. U21 is stored in TY. */
/*     Workspace:  need   NP2 + M1*M1 + max(3*NP2+M1,5*NP2); */
/*                 prefer larger. */

#line 432 "SB10PD.f"
    iq = np2 + 1;
#line 433 "SB10PD.f"
    iwrk = iq + m1 * m1;

#line 435 "SB10PD.f"
    i__1 = *ldwork - iwrk + 1;
#line 435 "SB10PD.f"
    dgesvd_("A", "A", &np2, &m1, &d__[np1 + 1 + d_dim1], ldd, &dwork[1], &ty[
	    ty_offset], ldty, &dwork[iq], &m1, &dwork[iwrk], &i__1, &info2, (
	    ftnlen)1, (ftnlen)1);
#line 438 "SB10PD.f"
    if (info2 != 0) {
#line 439 "SB10PD.f"
	*info = 5;
#line 440 "SB10PD.f"
	return 0;
#line 441 "SB10PD.f"
    }

#line 443 "SB10PD.f"
    rcond[2] = dwork[np2] / dwork[1];
#line 444 "SB10PD.f"
    if (rcond[2] <= toll) {
#line 445 "SB10PD.f"
	*info = 4;
#line 446 "SB10PD.f"
	return 0;
#line 447 "SB10PD.f"
    }
/* Computing MAX */
#line 448 "SB10PD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 448 "SB10PD.f"
    lwamax = max(i__1,lwamax);

/*     Determine Q21. */

#line 452 "SB10PD.f"
    if (nd2 > 0) {
#line 453 "SB10PD.f"
	dlacpy_("Full", &np2, &m1, &dwork[iq], &m1, &d__[np1 + 1 + d_dim1], 
		ldd, (ftnlen)4);
#line 455 "SB10PD.f"
	dlacpy_("Full", &nd2, &m1, &dwork[iq + np2], &m1, &dwork[iq], &m1, (
		ftnlen)4);
#line 457 "SB10PD.f"
	dlacpy_("Full", &np2, &m1, &d__[np1 + 1 + d_dim1], ldd, &dwork[iq + 
		nd2], &m1, (ftnlen)4);
#line 459 "SB10PD.f"
    }

/*     Determine Ty by scaling and transposing in-situ. */

#line 463 "SB10PD.f"
    i__1 = np2;
#line 463 "SB10PD.f"
    for (j = 1; j <= i__1; ++j) {
#line 464 "SB10PD.f"
	d__1 = 1. / dwork[j];
#line 464 "SB10PD.f"
	dscal_(&np2, &d__1, &ty[j * ty_dim1 + 1], &c__1);
#line 465 "SB10PD.f"
/* L30: */
#line 465 "SB10PD.f"
    }

#line 467 "SB10PD.f"
    i__1 = np2 - 1;
#line 467 "SB10PD.f"
    for (j = 1; j <= i__1; ++j) {
#line 468 "SB10PD.f"
	dswap_(&j, &ty[j + 1 + ty_dim1], ldty, &ty[(j + 1) * ty_dim1 + 1], &
		c__1);
#line 469 "SB10PD.f"
/* L40: */
#line 469 "SB10PD.f"
    }

/*     Determine B1 =: B1*Q21'. */
/*     Workspace:  NP2 + M1*M1 + N*M1. */

#line 474 "SB10PD.f"
    dgemm_("N", "T", n, &m1, &m1, &c_b27, &b[b_offset], ldb, &dwork[iq], &m1, 
	    &c_b28, &dwork[iwrk], n, (ftnlen)1, (ftnlen)1);
#line 476 "SB10PD.f"
    dlacpy_("Full", n, &m1, &dwork[iwrk], n, &b[b_offset], ldb, (ftnlen)4);
/* Computing MAX */
#line 477 "SB10PD.f"
    i__1 = iwrk + *n * m1 - 1;
#line 477 "SB10PD.f"
    lwamax = max(i__1,lwamax);

/*     Determine D11 =: D11*Q21'. */
/*     Workspace:  NP2 + M1*M1 + NP1*M1. */

#line 482 "SB10PD.f"
    dgemm_("N", "T", &np1, &m1, &m1, &c_b27, &d__[d_offset], ldd, &dwork[iq], 
	    &m1, &c_b28, &dwork[iwrk], &np1, (ftnlen)1, (ftnlen)1);
#line 484 "SB10PD.f"
    dlacpy_("Full", &np1, &m1, &dwork[iwrk], &np1, &d__[d_offset], ldd, (
	    ftnlen)4);
/* Computing MAX */
#line 485 "SB10PD.f"
    i__1 = iwrk + np1 * m1 - 1;
#line 485 "SB10PD.f"
    lwamax = max(i__1,lwamax);

/*     Determine B2 =: B2*Tu. */
/*     Workspace:  N*M2. */

#line 490 "SB10PD.f"
    dgemm_("N", "N", n, &m2, &m2, &c_b27, &b[(m1 + 1) * b_dim1 + 1], ldb, &tu[
	    tu_offset], ldtu, &c_b28, &dwork[1], n, (ftnlen)1, (ftnlen)1);
#line 492 "SB10PD.f"
    dlacpy_("Full", n, &m2, &dwork[1], n, &b[(m1 + 1) * b_dim1 + 1], ldb, (
	    ftnlen)4);

/*     Determine C2 =: Ty*C2. */
/*     Workspace:  NP2*N. */

#line 497 "SB10PD.f"
    dgemm_("N", "N", &np2, n, &np2, &c_b27, &ty[ty_offset], ldty, &c__[np1 + 
	    1 + c_dim1], ldc, &c_b28, &dwork[1], &np2, (ftnlen)1, (ftnlen)1);
#line 499 "SB10PD.f"
    dlacpy_("Full", &np2, n, &dwork[1], &np2, &c__[np1 + 1 + c_dim1], ldc, (
	    ftnlen)4);

/* Computing MAX */
#line 501 "SB10PD.f"
    i__1 = *n * max(m2,np2);
#line 501 "SB10PD.f"
    lwamax = max(i__1,lwamax);
#line 502 "SB10PD.f"
    dwork[1] = (doublereal) lwamax;
#line 503 "SB10PD.f"
    return 0;
/* *** Last line of SB10PD *** */
} /* sb10pd_ */

