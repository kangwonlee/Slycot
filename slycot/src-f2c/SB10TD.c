#line 1 "SB10TD.f"
/* SB10TD.f -- translated by f2c (version 20100827).
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

#line 1 "SB10TD.f"
/* Table of constant values */

static doublereal c_b6 = 1.;
static doublereal c_b7 = 0.;
static doublereal c_b39 = -1.;

/* Subroutine */ int sb10td_(integer *n, integer *m, integer *np, integer *
	ncon, integer *nmeas, doublereal *d__, integer *ldd, doublereal *tu, 
	integer *ldtu, doublereal *ty, integer *ldty, doublereal *ak, integer 
	*ldak, doublereal *bk, integer *ldbk, doublereal *ck, integer *ldck, 
	doublereal *dk, integer *lddk, doublereal *rcond, doublereal *tol, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer ak_dim1, ak_offset, bk_dim1, bk_offset, ck_dim1, ck_offset, 
	    d_dim1, d_offset, dk_dim1, dk_offset, tu_dim1, tu_offset, ty_dim1,
	     ty_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer m1, m2, np1, np2;
    static doublereal toll;
    static integer iwrk, info2;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal anorm;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgecon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dgetrf_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), dlacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), dlaset_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen), dgetrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static integer minwrk;


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

/*     To compute the matrices of the H2 optimal discrete-time controller */

/*              | AK | BK | */
/*          K = |----|----|, */
/*              | CK | DK | */

/*     from the matrices of the controller for the normalized system, */
/*     as determined by the SLICOT Library routine SB10SD. */

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

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading NP-by-M part of this array must contain the */
/*             system input/output matrix D. Only the trailing */
/*             NMEAS-by-NCON submatrix D22 is used. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= max(1,NP). */

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

/*     AK      (input/output) DOUBLE PRECISION array, dimension (LDAK,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain controller state matrix for the normalized system */
/*             as obtained by the SLICOT Library routine SB10SD. */
/*             On exit, the leading N-by-N part of this array contains */
/*             controller state matrix AK. */

/*     LDAK    INTEGER */
/*             The leading dimension of the array AK.  LDAK >= max(1,N). */

/*     BK      (input/output) DOUBLE PRECISION array, dimension */
/*             (LDBK,NMEAS) */
/*             On entry, the leading N-by-NMEAS part of this array must */
/*             contain controller input matrix for the normalized system */
/*             as obtained by the SLICOT Library routine SB10SD. */
/*             On exit, the leading N-by-NMEAS part of this array */
/*             contains controller input matrix BK. */

/*     LDBK    INTEGER */
/*             The leading dimension of the array BK.  LDBK >= max(1,N). */

/*     CK      (input/output) DOUBLE PRECISION array, dimension (LDCK,N) */
/*             On entry, the leading NCON-by-N part of this array must */
/*             contain controller output matrix for the normalized */
/*             system as obtained by the SLICOT Library routine SB10SD. */
/*             On exit, the leading NCON-by-N part of this array contains */
/*             controller output matrix CK. */

/*     LDCK    INTEGER */
/*             The leading dimension of the array CK. */
/*             LDCK >= max(1,NCON). */

/*     DK      (input/output) DOUBLE PRECISION array, dimension */
/*             (LDDK,NMEAS) */
/*             On entry, the leading NCON-by-NMEAS part of this array */
/*             must contain controller matrix DK for the normalized */
/*             system as obtained by the SLICOT Library routine SB10SD. */
/*             On exit, the leading NCON-by-NMEAS part of this array */
/*             contains controller input/output matrix DK. */

/*     LDDK    INTEGER */
/*             The leading dimension of the array DK. */
/*             LDDK >= max(1,NCON). */

/*     RCOND   (output) DOUBLE PRECISION */
/*             RCOND contains an estimate of the reciprocal condition */
/*             number of the matrix Im2 + DKHAT*D22 which must be */
/*             inverted in the computation of the controller. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             Tolerance used in determining the nonsingularity of the */
/*             matrix which must be inverted. If TOL <= 0, then a default */
/*             value equal to sqrt(EPS) is used, where EPS is the */
/*             relative machine precision. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (2*M2) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= max(N*M2,N*NP2,M2*NP2,M2*M2+4*M2). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the matrix Im2 + DKHAT*D22 is singular, or the */
/*                   estimated condition number is larger than or equal */
/*                   to 1/TOL. */

/*     METHOD */

/*     The routine implements the formulas given in [1]. */

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
/*     input and output transformations and of the matrix Im2 + */
/*     DKHAT*D22. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, April 1999. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 1999, */
/*     Jan. 2000. */

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

#line 229 "SB10TD.f"
    /* Parameter adjustments */
#line 229 "SB10TD.f"
    d_dim1 = *ldd;
#line 229 "SB10TD.f"
    d_offset = 1 + d_dim1;
#line 229 "SB10TD.f"
    d__ -= d_offset;
#line 229 "SB10TD.f"
    tu_dim1 = *ldtu;
#line 229 "SB10TD.f"
    tu_offset = 1 + tu_dim1;
#line 229 "SB10TD.f"
    tu -= tu_offset;
#line 229 "SB10TD.f"
    ty_dim1 = *ldty;
#line 229 "SB10TD.f"
    ty_offset = 1 + ty_dim1;
#line 229 "SB10TD.f"
    ty -= ty_offset;
#line 229 "SB10TD.f"
    ak_dim1 = *ldak;
#line 229 "SB10TD.f"
    ak_offset = 1 + ak_dim1;
#line 229 "SB10TD.f"
    ak -= ak_offset;
#line 229 "SB10TD.f"
    bk_dim1 = *ldbk;
#line 229 "SB10TD.f"
    bk_offset = 1 + bk_dim1;
#line 229 "SB10TD.f"
    bk -= bk_offset;
#line 229 "SB10TD.f"
    ck_dim1 = *ldck;
#line 229 "SB10TD.f"
    ck_offset = 1 + ck_dim1;
#line 229 "SB10TD.f"
    ck -= ck_offset;
#line 229 "SB10TD.f"
    dk_dim1 = *lddk;
#line 229 "SB10TD.f"
    dk_offset = 1 + dk_dim1;
#line 229 "SB10TD.f"
    dk -= dk_offset;
#line 229 "SB10TD.f"
    --iwork;
#line 229 "SB10TD.f"
    --dwork;
#line 229 "SB10TD.f"

#line 229 "SB10TD.f"
    /* Function Body */
#line 229 "SB10TD.f"
    m1 = *m - *ncon;
#line 230 "SB10TD.f"
    m2 = *ncon;
#line 231 "SB10TD.f"
    np1 = *np - *nmeas;
#line 232 "SB10TD.f"
    np2 = *nmeas;

#line 234 "SB10TD.f"
    *info = 0;
#line 235 "SB10TD.f"
    if (*n < 0) {
#line 236 "SB10TD.f"
	*info = -1;
#line 237 "SB10TD.f"
    } else if (*m < 0) {
#line 238 "SB10TD.f"
	*info = -2;
#line 239 "SB10TD.f"
    } else if (*np < 0) {
#line 240 "SB10TD.f"
	*info = -3;
#line 241 "SB10TD.f"
    } else if (*ncon < 0 || m1 < 0 || m2 > np1) {
#line 242 "SB10TD.f"
	*info = -4;
#line 243 "SB10TD.f"
    } else if (*nmeas < 0 || np1 < 0 || np2 > m1) {
#line 244 "SB10TD.f"
	*info = -5;
#line 245 "SB10TD.f"
    } else if (*ldd < max(1,*np)) {
#line 246 "SB10TD.f"
	*info = -7;
#line 247 "SB10TD.f"
    } else if (*ldtu < max(1,m2)) {
#line 248 "SB10TD.f"
	*info = -9;
#line 249 "SB10TD.f"
    } else if (*ldty < max(1,np2)) {
#line 250 "SB10TD.f"
	*info = -11;
#line 251 "SB10TD.f"
    } else if (*ldak < max(1,*n)) {
#line 252 "SB10TD.f"
	*info = -13;
#line 253 "SB10TD.f"
    } else if (*ldbk < max(1,*n)) {
#line 254 "SB10TD.f"
	*info = -15;
#line 255 "SB10TD.f"
    } else if (*ldck < max(1,m2)) {
#line 256 "SB10TD.f"
	*info = -17;
#line 257 "SB10TD.f"
    } else if (*lddk < max(1,m2)) {
#line 258 "SB10TD.f"
	*info = -19;
#line 259 "SB10TD.f"
    } else {

/*        Compute workspace. */

/* Computing MAX */
#line 263 "SB10TD.f"
	i__1 = *n * m2, i__2 = *n * np2, i__1 = max(i__1,i__2), i__2 = m2 * 
		np2, i__1 = max(i__1,i__2), i__2 = m2 * (m2 + 4);
#line 263 "SB10TD.f"
	minwrk = max(i__1,i__2);
#line 264 "SB10TD.f"
	if (*ldwork < minwrk) {
#line 264 "SB10TD.f"
	    *info = -24;
#line 264 "SB10TD.f"
	}
#line 266 "SB10TD.f"
    }
#line 267 "SB10TD.f"
    if (*info != 0) {
#line 268 "SB10TD.f"
	i__1 = -(*info);
#line 268 "SB10TD.f"
	xerbla_("SB10TD", &i__1, (ftnlen)6);
#line 269 "SB10TD.f"
	return 0;
#line 270 "SB10TD.f"
    }

/*     Quick return if possible. */

#line 274 "SB10TD.f"
    if (*n == 0 || *m == 0 || *np == 0 || m1 == 0 || m2 == 0 || np1 == 0 || 
	    np2 == 0) {
#line 276 "SB10TD.f"
	*rcond = 1.;
#line 277 "SB10TD.f"
	return 0;
#line 278 "SB10TD.f"
    }

#line 280 "SB10TD.f"
    toll = *tol;
#line 281 "SB10TD.f"
    if (toll <= 0.) {

/*        Set the default value of the tolerance for nonsingularity test. */

#line 285 "SB10TD.f"
	toll = sqrt(dlamch_("Epsilon", (ftnlen)7));
#line 286 "SB10TD.f"
    }

/*     Find BKHAT . */

#line 290 "SB10TD.f"
    dgemm_("N", "N", n, &np2, &np2, &c_b6, &bk[bk_offset], ldbk, &ty[
	    ty_offset], ldty, &c_b7, &dwork[1], n, (ftnlen)1, (ftnlen)1);
#line 292 "SB10TD.f"
    dlacpy_("Full", n, &np2, &dwork[1], n, &bk[bk_offset], ldbk, (ftnlen)4);

/*     Find CKHAT . */

#line 296 "SB10TD.f"
    dgemm_("N", "N", &m2, n, &m2, &c_b6, &tu[tu_offset], ldtu, &ck[ck_offset],
	     ldck, &c_b7, &dwork[1], &m2, (ftnlen)1, (ftnlen)1);
#line 298 "SB10TD.f"
    dlacpy_("Full", &m2, n, &dwork[1], &m2, &ck[ck_offset], ldck, (ftnlen)4);

/*     Compute DKHAT . */

#line 302 "SB10TD.f"
    dgemm_("N", "N", &m2, &np2, &m2, &c_b6, &tu[tu_offset], ldtu, &dk[
	    dk_offset], lddk, &c_b7, &dwork[1], &m2, (ftnlen)1, (ftnlen)1);
#line 304 "SB10TD.f"
    dgemm_("N", "N", &m2, &np2, &np2, &c_b6, &dwork[1], &m2, &ty[ty_offset], 
	    ldty, &c_b7, &dk[dk_offset], lddk, (ftnlen)1, (ftnlen)1);

/*     Compute Im2 + DKHAT*D22 . */

#line 309 "SB10TD.f"
    iwrk = m2 * m2 + 1;
#line 310 "SB10TD.f"
    dlaset_("Full", &m2, &m2, &c_b7, &c_b6, &dwork[1], &m2, (ftnlen)4);
#line 311 "SB10TD.f"
    dgemm_("N", "N", &m2, &m2, &np2, &c_b6, &dk[dk_offset], lddk, &d__[np1 + 
	    1 + (m1 + 1) * d_dim1], ldd, &c_b6, &dwork[1], &m2, (ftnlen)1, (
	    ftnlen)1);
#line 313 "SB10TD.f"
    anorm = dlange_("1", &m2, &m2, &dwork[1], &m2, &dwork[iwrk], (ftnlen)1);
#line 314 "SB10TD.f"
    dgetrf_(&m2, &m2, &dwork[1], &m2, &iwork[1], &info2);
#line 315 "SB10TD.f"
    if (info2 > 0) {
#line 316 "SB10TD.f"
	*info = 1;
#line 317 "SB10TD.f"
	return 0;
#line 318 "SB10TD.f"
    }
#line 319 "SB10TD.f"
    dgecon_("1", &m2, &dwork[1], &m2, &anorm, rcond, &dwork[iwrk], &iwork[m2 
	    + 1], &info2, (ftnlen)1);

/*     Return if the matrix is singular to working precision. */

#line 324 "SB10TD.f"
    if (*rcond < toll) {
#line 325 "SB10TD.f"
	*info = 1;
#line 326 "SB10TD.f"
	return 0;
#line 327 "SB10TD.f"
    }

/*     Compute CK . */

#line 331 "SB10TD.f"
    dgetrs_("N", &m2, n, &dwork[1], &m2, &iwork[1], &ck[ck_offset], ldck, &
	    info2, (ftnlen)1);

/*     Compute DK . */

#line 335 "SB10TD.f"
    dgetrs_("N", &m2, &np2, &dwork[1], &m2, &iwork[1], &dk[dk_offset], lddk, &
	    info2, (ftnlen)1);

/*     Compute AK . */

#line 339 "SB10TD.f"
    dgemm_("N", "N", n, &m2, &np2, &c_b6, &bk[bk_offset], ldbk, &d__[np1 + 1 
	    + (m1 + 1) * d_dim1], ldd, &c_b7, &dwork[1], n, (ftnlen)1, (
	    ftnlen)1);
#line 341 "SB10TD.f"
    dgemm_("N", "N", n, n, &m2, &c_b39, &dwork[1], n, &ck[ck_offset], ldck, &
	    c_b6, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);

/*     Compute BK . */

#line 346 "SB10TD.f"
    dgemm_("N", "N", n, &np2, &m2, &c_b39, &dwork[1], n, &dk[dk_offset], lddk,
	     &c_b6, &bk[bk_offset], ldbk, (ftnlen)1, (ftnlen)1);
#line 348 "SB10TD.f"
    return 0;
/* *** Last line of SB10TD *** */
} /* sb10td_ */

