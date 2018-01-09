#line 1 "SB10UD.f"
/* SB10UD.f -- translated by f2c (version 20100827).
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

#line 1 "SB10UD.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b15 = 1.;
static doublereal c_b16 = 0.;

/* Subroutine */ int sb10ud_(integer *n, integer *m, integer *np, integer *
	ncon, integer *nmeas, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, doublereal *d__, integer *ldd, doublereal *tu, integer *
	ldtu, doublereal *ty, integer *ldty, doublereal *rcond, doublereal *
	tol, doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, c_dim1, c_offset, d_dim1, d_offset, tu_dim1, 
	    tu_offset, ty_dim1, ty_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, m1, m2, iq, nd1, nd2, np1, np2;
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
/*                   | C1 |  0  D12 |   | C | D | */
/*                   | C2 | D21 D22 | */

/*     to unit diagonal form, and to transform the matrices B and C to */
/*     satisfy the formulas in the computation of the H2 optimal */
/*     controller. */

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
/*             contain the system input/output matrix D. */
/*             The (NP-NMEAS)-by-(M-NCON) leading submatrix D11 is not */
/*             referenced. */
/*             On exit, the trailing NMEAS-by-NCON part (in the leading */
/*             NP-by-M part) of this array contains the transformed */
/*             submatrix D22. */
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
/*             RCOND is set even if INFO = 1 or INFO = 2; if INFO = 1, */
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
/*             LDWORK >= MAX( M2 + NP1*NP1 + MAX(NP1*N,3*M2+NP1,5*M2), */
/*                            NP2 + M1*M1  + MAX(M1*N,3*NP2+M1,5*NP2), */
/*                            N*M2, NP2*N, NP2*M2, 1 ) */
/*             where M1 = M - M2 and NP1 = NP - NP2. */
/*             For good performance, LDWORK must generally be larger. */
/*             Denoting Q = MAX(M1,M2,NP1,NP2), an upper bound is */
/*             MAX(1,Q*(Q+MAX(N,5)+1)). */

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
/*                   did not converge (when computing the SVD of D12 or */
/*                   D21). */

/*     METHOD */

/*     The routine performs the transformations described in [1], [2]. */

/*     REFERENCES */

/*     [1] Zhou, K., Doyle, J.C., and Glover, K. */
/*         Robust and Optimal Control. */
/*         Prentice-Hall, Upper Saddle River, NJ, 1996. */

/*     [2] Balas, G.J., Doyle, J.C., Glover, K., Packard, A., and */
/*         Smith, R. */
/*         mu-Analysis and Synthesis Toolbox. */
/*         The MathWorks Inc., Natick, Mass., 1995. */

/*     NUMERICAL ASPECTS */

/*     The precision of the transformations can be controlled by the */
/*     condition numbers of the matrices TU and TY as given by the */
/*     values of RCOND(1) and RCOND(2), respectively. An error return */
/*     with INFO = 1 or INFO = 2 will be obtained if the condition */
/*     number of TU or TY, respectively, would exceed 1/TOL. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 1998. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 1999, */
/*     Feb. 2000. */

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
/*     .. External Functions */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test input parameters. */

#line 226 "SB10UD.f"
    /* Parameter adjustments */
#line 226 "SB10UD.f"
    b_dim1 = *ldb;
#line 226 "SB10UD.f"
    b_offset = 1 + b_dim1;
#line 226 "SB10UD.f"
    b -= b_offset;
#line 226 "SB10UD.f"
    c_dim1 = *ldc;
#line 226 "SB10UD.f"
    c_offset = 1 + c_dim1;
#line 226 "SB10UD.f"
    c__ -= c_offset;
#line 226 "SB10UD.f"
    d_dim1 = *ldd;
#line 226 "SB10UD.f"
    d_offset = 1 + d_dim1;
#line 226 "SB10UD.f"
    d__ -= d_offset;
#line 226 "SB10UD.f"
    tu_dim1 = *ldtu;
#line 226 "SB10UD.f"
    tu_offset = 1 + tu_dim1;
#line 226 "SB10UD.f"
    tu -= tu_offset;
#line 226 "SB10UD.f"
    ty_dim1 = *ldty;
#line 226 "SB10UD.f"
    ty_offset = 1 + ty_dim1;
#line 226 "SB10UD.f"
    ty -= ty_offset;
#line 226 "SB10UD.f"
    --rcond;
#line 226 "SB10UD.f"
    --dwork;
#line 226 "SB10UD.f"

#line 226 "SB10UD.f"
    /* Function Body */
#line 226 "SB10UD.f"
    m1 = *m - *ncon;
#line 227 "SB10UD.f"
    m2 = *ncon;
#line 228 "SB10UD.f"
    np1 = *np - *nmeas;
#line 229 "SB10UD.f"
    np2 = *nmeas;

#line 231 "SB10UD.f"
    *info = 0;
#line 232 "SB10UD.f"
    if (*n < 0) {
#line 233 "SB10UD.f"
	*info = -1;
#line 234 "SB10UD.f"
    } else if (*m < 0) {
#line 235 "SB10UD.f"
	*info = -2;
#line 236 "SB10UD.f"
    } else if (*np < 0) {
#line 237 "SB10UD.f"
	*info = -3;
#line 238 "SB10UD.f"
    } else if (*ncon < 0 || m1 < 0 || m2 > np1) {
#line 239 "SB10UD.f"
	*info = -4;
#line 240 "SB10UD.f"
    } else if (*nmeas < 0 || np1 < 0 || np2 > m1) {
#line 241 "SB10UD.f"
	*info = -5;
#line 242 "SB10UD.f"
    } else if (*ldb < max(1,*n)) {
#line 243 "SB10UD.f"
	*info = -7;
#line 244 "SB10UD.f"
    } else if (*ldc < max(1,*np)) {
#line 245 "SB10UD.f"
	*info = -9;
#line 246 "SB10UD.f"
    } else if (*ldd < max(1,*np)) {
#line 247 "SB10UD.f"
	*info = -11;
#line 248 "SB10UD.f"
    } else if (*ldtu < max(1,m2)) {
#line 249 "SB10UD.f"
	*info = -13;
#line 250 "SB10UD.f"
    } else if (*ldty < max(1,np2)) {
#line 251 "SB10UD.f"
	*info = -15;
#line 252 "SB10UD.f"
    } else {

/*        Compute workspace. */

/* Computing MAX */
/* Computing MAX */
#line 256 "SB10UD.f"
	i__3 = np1 * *n, i__4 = m2 * 3 + np1, i__3 = max(i__3,i__4), i__4 = 
		m2 * 5;
/* Computing MAX */
#line 256 "SB10UD.f"
	i__5 = m1 * *n, i__6 = np2 * 3 + m1, i__5 = max(i__5,i__6), i__6 = 
		np2 * 5;
#line 256 "SB10UD.f"
	i__1 = 1, i__2 = m2 + np1 * np1 + max(i__3,i__4), i__1 = max(i__1,
		i__2), i__2 = np2 + m1 * m1 + max(i__5,i__6), i__1 = max(i__1,
		i__2), i__2 = *n * m2, i__1 = max(i__1,i__2), i__2 = np2 * *n,
		 i__1 = max(i__1,i__2), i__2 = np2 * m2;
#line 256 "SB10UD.f"
	minwrk = max(i__1,i__2);
#line 260 "SB10UD.f"
	if (*ldwork < minwrk) {
#line 260 "SB10UD.f"
	    *info = -19;
#line 260 "SB10UD.f"
	}
#line 262 "SB10UD.f"
    }
#line 263 "SB10UD.f"
    if (*info != 0) {
#line 264 "SB10UD.f"
	i__1 = -(*info);
#line 264 "SB10UD.f"
	xerbla_("SB10UD", &i__1, (ftnlen)6);
#line 265 "SB10UD.f"
	return 0;
#line 266 "SB10UD.f"
    }

/*     Quick return if possible. */

#line 270 "SB10UD.f"
    if (*n == 0 || *m == 0 || *np == 0 || m1 == 0 || m2 == 0 || np1 == 0 || 
	    np2 == 0) {
#line 272 "SB10UD.f"
	rcond[1] = 1.;
#line 273 "SB10UD.f"
	rcond[2] = 1.;
#line 274 "SB10UD.f"
	dwork[1] = 1.;
#line 275 "SB10UD.f"
	return 0;
#line 276 "SB10UD.f"
    }

#line 278 "SB10UD.f"
    nd1 = np1 - m2;
#line 279 "SB10UD.f"
    nd2 = m1 - np2;
#line 280 "SB10UD.f"
    toll = *tol;
#line 281 "SB10UD.f"
    if (toll <= 0.) {

/*        Set the default value of the tolerance for condition tests. */

#line 285 "SB10UD.f"
	toll = sqrt(dlamch_("Epsilon", (ftnlen)7));
#line 286 "SB10UD.f"
    }

/*     Determine SVD of D12, D12 = U12 S12 V12', and check if D12 has */
/*     full column rank. V12' is stored in TU. */
/*     Workspace:  need   M2 + NP1*NP1 + max(3*M2+NP1,5*M2); */
/*                 prefer larger. */

#line 293 "SB10UD.f"
    iq = m2 + 1;
#line 294 "SB10UD.f"
    iwrk = iq + np1 * np1;

#line 296 "SB10UD.f"
    i__1 = *ldwork - iwrk + 1;
#line 296 "SB10UD.f"
    dgesvd_("A", "A", &np1, &m2, &d__[(m1 + 1) * d_dim1 + 1], ldd, &dwork[1], 
	    &dwork[iq], &np1, &tu[tu_offset], ldtu, &dwork[iwrk], &i__1, &
	    info2, (ftnlen)1, (ftnlen)1);
#line 299 "SB10UD.f"
    if (info2 != 0) {
#line 300 "SB10UD.f"
	*info = 3;
#line 301 "SB10UD.f"
	return 0;
#line 302 "SB10UD.f"
    }

#line 304 "SB10UD.f"
    rcond[1] = dwork[m2] / dwork[1];
#line 305 "SB10UD.f"
    if (rcond[1] <= toll) {
#line 306 "SB10UD.f"
	rcond[2] = 0.;
#line 307 "SB10UD.f"
	*info = 1;
#line 308 "SB10UD.f"
	return 0;
#line 309 "SB10UD.f"
    }
#line 310 "SB10UD.f"
    lwamax = (integer) dwork[iwrk] + iwrk - 1;

/*     Determine Q12. */

#line 314 "SB10UD.f"
    if (nd1 > 0) {
#line 315 "SB10UD.f"
	dlacpy_("Full", &np1, &m2, &dwork[iq], &np1, &d__[(m1 + 1) * d_dim1 + 
		1], ldd, (ftnlen)4);
#line 317 "SB10UD.f"
	dlacpy_("Full", &np1, &nd1, &dwork[iq + np1 * m2], &np1, &dwork[iq], &
		np1, (ftnlen)4);
#line 319 "SB10UD.f"
	dlacpy_("Full", &np1, &m2, &d__[(m1 + 1) * d_dim1 + 1], ldd, &dwork[
		iq + np1 * nd1], &np1, (ftnlen)4);
#line 321 "SB10UD.f"
    }

/*     Determine Tu by transposing in-situ and scaling. */

#line 325 "SB10UD.f"
    i__1 = m2 - 1;
#line 325 "SB10UD.f"
    for (j = 1; j <= i__1; ++j) {
#line 326 "SB10UD.f"
	dswap_(&j, &tu[j + 1 + tu_dim1], ldtu, &tu[(j + 1) * tu_dim1 + 1], &
		c__1);
#line 327 "SB10UD.f"
/* L10: */
#line 327 "SB10UD.f"
    }

#line 329 "SB10UD.f"
    i__1 = m2;
#line 329 "SB10UD.f"
    for (j = 1; j <= i__1; ++j) {
#line 330 "SB10UD.f"
	d__1 = 1. / dwork[j];
#line 330 "SB10UD.f"
	dscal_(&m2, &d__1, &tu[j * tu_dim1 + 1], &c__1);
#line 331 "SB10UD.f"
/* L20: */
#line 331 "SB10UD.f"
    }

/*     Determine C1 =: Q12'*C1. */
/*     Workspace:  M2 + NP1*NP1 + NP1*N. */

#line 336 "SB10UD.f"
    dgemm_("T", "N", &np1, n, &np1, &c_b15, &dwork[iq], &np1, &c__[c_offset], 
	    ldc, &c_b16, &dwork[iwrk], &np1, (ftnlen)1, (ftnlen)1);
#line 338 "SB10UD.f"
    dlacpy_("Full", &np1, n, &dwork[iwrk], &np1, &c__[c_offset], ldc, (ftnlen)
	    4);
/* Computing MAX */
#line 339 "SB10UD.f"
    i__1 = iwrk + np1 * *n - 1;
#line 339 "SB10UD.f"
    lwamax = max(i__1,lwamax);

/*     Determine SVD of D21, D21 = U21 S21 V21', and check if D21 has */
/*     full row rank. U21 is stored in TY. */
/*     Workspace:  need   NP2 + M1*M1 + max(3*NP2+M1,5*NP2); */
/*                 prefer larger. */

#line 346 "SB10UD.f"
    iq = np2 + 1;
#line 347 "SB10UD.f"
    iwrk = iq + m1 * m1;

#line 349 "SB10UD.f"
    i__1 = *ldwork - iwrk + 1;
#line 349 "SB10UD.f"
    dgesvd_("A", "A", &np2, &m1, &d__[np1 + 1 + d_dim1], ldd, &dwork[1], &ty[
	    ty_offset], ldty, &dwork[iq], &m1, &dwork[iwrk], &i__1, &info2, (
	    ftnlen)1, (ftnlen)1);
#line 352 "SB10UD.f"
    if (info2 != 0) {
#line 353 "SB10UD.f"
	*info = 3;
#line 354 "SB10UD.f"
	return 0;
#line 355 "SB10UD.f"
    }

#line 357 "SB10UD.f"
    rcond[2] = dwork[np2] / dwork[1];
#line 358 "SB10UD.f"
    if (rcond[2] <= toll) {
#line 359 "SB10UD.f"
	*info = 2;
#line 360 "SB10UD.f"
	return 0;
#line 361 "SB10UD.f"
    }
/* Computing MAX */
#line 362 "SB10UD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 362 "SB10UD.f"
    lwamax = max(i__1,lwamax);

/*     Determine Q21. */

#line 366 "SB10UD.f"
    if (nd2 > 0) {
#line 367 "SB10UD.f"
	dlacpy_("Full", &np2, &m1, &dwork[iq], &m1, &d__[np1 + 1 + d_dim1], 
		ldd, (ftnlen)4);
#line 369 "SB10UD.f"
	dlacpy_("Full", &nd2, &m1, &dwork[iq + np2], &m1, &dwork[iq], &m1, (
		ftnlen)4);
#line 371 "SB10UD.f"
	dlacpy_("Full", &np2, &m1, &d__[np1 + 1 + d_dim1], ldd, &dwork[iq + 
		nd2], &m1, (ftnlen)4);
#line 373 "SB10UD.f"
    }

/*     Determine Ty by scaling and transposing in-situ. */

#line 377 "SB10UD.f"
    i__1 = np2;
#line 377 "SB10UD.f"
    for (j = 1; j <= i__1; ++j) {
#line 378 "SB10UD.f"
	d__1 = 1. / dwork[j];
#line 378 "SB10UD.f"
	dscal_(&np2, &d__1, &ty[j * ty_dim1 + 1], &c__1);
#line 379 "SB10UD.f"
/* L30: */
#line 379 "SB10UD.f"
    }

#line 381 "SB10UD.f"
    i__1 = np2 - 1;
#line 381 "SB10UD.f"
    for (j = 1; j <= i__1; ++j) {
#line 382 "SB10UD.f"
	dswap_(&j, &ty[j + 1 + ty_dim1], ldty, &ty[(j + 1) * ty_dim1 + 1], &
		c__1);
#line 383 "SB10UD.f"
/* L40: */
#line 383 "SB10UD.f"
    }

/*     Determine B1 =: B1*Q21'. */
/*     Workspace:  NP2 + M1*M1 + N*M1. */

#line 388 "SB10UD.f"
    dgemm_("N", "T", n, &m1, &m1, &c_b15, &b[b_offset], ldb, &dwork[iq], &m1, 
	    &c_b16, &dwork[iwrk], n, (ftnlen)1, (ftnlen)1);
#line 390 "SB10UD.f"
    dlacpy_("Full", n, &m1, &dwork[iwrk], n, &b[b_offset], ldb, (ftnlen)4);
/* Computing MAX */
#line 391 "SB10UD.f"
    i__1 = iwrk + *n * m1 - 1;
#line 391 "SB10UD.f"
    lwamax = max(i__1,lwamax);

/*     Determine B2 =: B2*Tu. */
/*     Workspace:  N*M2. */

#line 396 "SB10UD.f"
    dgemm_("N", "N", n, &m2, &m2, &c_b15, &b[(m1 + 1) * b_dim1 + 1], ldb, &tu[
	    tu_offset], ldtu, &c_b16, &dwork[1], n, (ftnlen)1, (ftnlen)1);
#line 398 "SB10UD.f"
    dlacpy_("Full", n, &m2, &dwork[1], n, &b[(m1 + 1) * b_dim1 + 1], ldb, (
	    ftnlen)4);

/*     Determine C2 =: Ty*C2. */
/*     Workspace:  NP2*N. */

#line 403 "SB10UD.f"
    dgemm_("N", "N", &np2, n, &np2, &c_b15, &ty[ty_offset], ldty, &c__[np1 + 
	    1 + c_dim1], ldc, &c_b16, &dwork[1], &np2, (ftnlen)1, (ftnlen)1);
#line 405 "SB10UD.f"
    dlacpy_("Full", &np2, n, &dwork[1], &np2, &c__[np1 + 1 + c_dim1], ldc, (
	    ftnlen)4);

/*     Determine D22 =: Ty*D22*Tu. */
/*     Workspace:  NP2*M2. */

#line 410 "SB10UD.f"
    dgemm_("N", "N", &np2, &m2, &np2, &c_b15, &ty[ty_offset], ldty, &d__[np1 
	    + 1 + (m1 + 1) * d_dim1], ldd, &c_b16, &dwork[1], &np2, (ftnlen)1,
	     (ftnlen)1);
#line 412 "SB10UD.f"
    dgemm_("N", "N", &np2, &m2, &m2, &c_b15, &dwork[1], &np2, &tu[tu_offset], 
	    ldtu, &c_b16, &d__[np1 + 1 + (m1 + 1) * d_dim1], ldd, (ftnlen)1, (
	    ftnlen)1);

/* Computing MAX */
#line 415 "SB10UD.f"
    i__1 = *n * max(m2,np2), i__2 = np2 * m2, i__1 = max(i__1,i__2);
#line 415 "SB10UD.f"
    lwamax = max(i__1,lwamax);
#line 416 "SB10UD.f"
    dwork[1] = (doublereal) lwamax;
#line 417 "SB10UD.f"
    return 0;
/* *** Last line of SB10UD *** */
} /* sb10ud_ */
