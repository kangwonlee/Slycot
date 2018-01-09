#line 1 "SB10LD.f"
/* SB10LD.f -- translated by f2c (version 20100827).
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

#line 1 "SB10LD.f"
/* Table of constant values */

static doublereal c_b5 = 0.;
static doublereal c_b6 = 1.;
static doublereal c_b9 = -1.;

/* Subroutine */ int sb10ld_(integer *n, integer *m, integer *np, integer *
	ncon, integer *nmeas, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *d__, integer 
	*ldd, doublereal *ak, integer *ldak, doublereal *bk, integer *ldbk, 
	doublereal *ck, integer *ldck, doublereal *dk, integer *lddk, 
	doublereal *ac, integer *ldac, doublereal *bc, integer *ldbc, 
	doublereal *cc, integer *ldcc, doublereal *dc, integer *lddc, integer 
	*iwork, doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ac_dim1, ac_offset, ak_dim1, ak_offset, b_dim1, 
	    b_offset, bc_dim1, bc_offset, bk_dim1, bk_offset, c_dim1, 
	    c_offset, cc_dim1, cc_offset, ck_dim1, ck_offset, d_dim1, 
	    d_offset, dc_dim1, dc_offset, dk_dim1, dk_offset, i__1;

    /* Local variables */
    static integer m1, m2, n2, np1, np2, iw2, iw3, iw4, iw5, iw6, iw7, iw8;
    static doublereal eps;
    static integer iwrk, info2;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal rcond, anorm;
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
	    integer *), xerbla_(char *, integer *, ftnlen);
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

/*     To compute the matrices of the closed-loop system */

/*              | AC | BC | */
/*          G = |----|----|, */
/*              | CC | DC | */

/*     from the matrices of the open-loop system */

/*               | A | B | */
/*           P = |---|---| */
/*               | C | D | */

/*     and the matrices of the controller */

/*              | AK | BK | */
/*          K = |----|----|. */
/*              | CK | DK | */

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

/*     AK      (input) DOUBLE PRECISION array, dimension (LDAK,N) */
/*             The leading N-by-N part of this array must contain the */
/*             controller state matrix AK. */

/*     LDAK    INTEGER */
/*             The leading dimension of the array AK.  LDAK >= max(1,N). */

/*     BK      (input) DOUBLE PRECISION array, dimension (LDBK,NMEAS) */
/*             The leading N-by-NMEAS part of this array must contain the */
/*             controller input matrix BK. */

/*     LDBK    INTEGER */
/*             The leading dimension of the array BK.  LDBK >= max(1,N). */

/*     CK      (input) DOUBLE PRECISION array, dimension (LDCK,N) */
/*             The leading NCON-by-N part of this array must contain the */
/*             controller output matrix CK. */

/*     LDCK    INTEGER */
/*             The leading dimension of the array CK. */
/*             LDCK >= max(1,NCON). */

/*     DK      (input) DOUBLE PRECISION array, dimension (LDDK,NMEAS) */
/*             The leading NCON-by-NMEAS part of this array must contain */
/*             the controller input/output matrix DK. */

/*     LDDK    INTEGER */
/*             The leading dimension of the array DK. */
/*             LDDK >= max(1,NCON). */

/*     AC      (output) DOUBLE PRECISION array, dimension (LDAC,2*N) */
/*             The leading 2*N-by-2*N part of this array contains the */
/*             closed-loop system state matrix AC. */

/*     LDAC    INTEGER */
/*             The leading dimension of the array AC. */
/*             LDAC >= max(1,2*N). */

/*     BC      (output) DOUBLE PRECISION array, dimension (LDBC,M-NCON) */
/*             The leading 2*N-by-(M-NCON) part of this array contains */
/*             the closed-loop system input matrix BC. */

/*     LDBC    INTEGER */
/*             The leading dimension of the array BC. */
/*             LDBC >= max(1,2*N). */

/*     CC      (output) DOUBLE PRECISION array, dimension (LDCC,2*N) */
/*             The leading (NP-NMEAS)-by-2*N part of this array contains */
/*             the closed-loop system output matrix CC. */

/*     LDCC    INTEGER */
/*             The leading dimension of the array CC. */
/*             LDCC >= max(1,NP-NMEAS). */

/*     DC      (output) DOUBLE PRECISION array, dimension (LDDC,M-NCON) */
/*             The leading (NP-NMEAS)-by-(M-NCON) part of this array */
/*             contains the closed-loop system input/output matrix DC. */

/*     LDDC    INTEGER */
/*             The leading dimension of the array DC. */
/*             LDDC >= max(1,NP-NMEAS). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension 2*max(NCON,NMEAS) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal */
/*             LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= 2*M*M+NP*NP+2*M*N+M*NP+2*N*NP. */
/*             For good performance, LDWORK must generally be larger. */

/*     Error Indicactor */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the matrix Inp2 - D22*DK is singular to working */
/*                   precision; */
/*             = 2:  if the matrix Im2 - DK*D22 is singular to working */
/*                   precision. */

/*     METHOD */

/*     The routine implements the formulas given in [1]. */

/*     REFERENCES */

/*     [1] Balas, G.J., Doyle, J.C., Glover, K., Packard, A., and */
/*         Smith, R. */
/*         mu-Analysis and Synthesis Toolbox. */
/*         The MathWorks Inc., Natick, Mass., 1995. */

/*     NUMERICAL ASPECTS */

/*     The accuracy of the result depends on the condition numbers of the */
/*     matrices  Inp2 - D22*DK  and  Im2 - DK*D22. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 1998. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 1999. */
/*     A. Markovski, Technical University, Sofia, April, 2003. */

/*     KEYWORDS */

/*     Closed loop systems, feedback control, robust control. */

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
/*     .. Executable Statements .. */

/*     Decode and Test input parameters. */

#line 244 "SB10LD.f"
    /* Parameter adjustments */
#line 244 "SB10LD.f"
    a_dim1 = *lda;
#line 244 "SB10LD.f"
    a_offset = 1 + a_dim1;
#line 244 "SB10LD.f"
    a -= a_offset;
#line 244 "SB10LD.f"
    b_dim1 = *ldb;
#line 244 "SB10LD.f"
    b_offset = 1 + b_dim1;
#line 244 "SB10LD.f"
    b -= b_offset;
#line 244 "SB10LD.f"
    c_dim1 = *ldc;
#line 244 "SB10LD.f"
    c_offset = 1 + c_dim1;
#line 244 "SB10LD.f"
    c__ -= c_offset;
#line 244 "SB10LD.f"
    d_dim1 = *ldd;
#line 244 "SB10LD.f"
    d_offset = 1 + d_dim1;
#line 244 "SB10LD.f"
    d__ -= d_offset;
#line 244 "SB10LD.f"
    ak_dim1 = *ldak;
#line 244 "SB10LD.f"
    ak_offset = 1 + ak_dim1;
#line 244 "SB10LD.f"
    ak -= ak_offset;
#line 244 "SB10LD.f"
    bk_dim1 = *ldbk;
#line 244 "SB10LD.f"
    bk_offset = 1 + bk_dim1;
#line 244 "SB10LD.f"
    bk -= bk_offset;
#line 244 "SB10LD.f"
    ck_dim1 = *ldck;
#line 244 "SB10LD.f"
    ck_offset = 1 + ck_dim1;
#line 244 "SB10LD.f"
    ck -= ck_offset;
#line 244 "SB10LD.f"
    dk_dim1 = *lddk;
#line 244 "SB10LD.f"
    dk_offset = 1 + dk_dim1;
#line 244 "SB10LD.f"
    dk -= dk_offset;
#line 244 "SB10LD.f"
    ac_dim1 = *ldac;
#line 244 "SB10LD.f"
    ac_offset = 1 + ac_dim1;
#line 244 "SB10LD.f"
    ac -= ac_offset;
#line 244 "SB10LD.f"
    bc_dim1 = *ldbc;
#line 244 "SB10LD.f"
    bc_offset = 1 + bc_dim1;
#line 244 "SB10LD.f"
    bc -= bc_offset;
#line 244 "SB10LD.f"
    cc_dim1 = *ldcc;
#line 244 "SB10LD.f"
    cc_offset = 1 + cc_dim1;
#line 244 "SB10LD.f"
    cc -= cc_offset;
#line 244 "SB10LD.f"
    dc_dim1 = *lddc;
#line 244 "SB10LD.f"
    dc_offset = 1 + dc_dim1;
#line 244 "SB10LD.f"
    dc -= dc_offset;
#line 244 "SB10LD.f"
    --iwork;
#line 244 "SB10LD.f"
    --dwork;
#line 244 "SB10LD.f"

#line 244 "SB10LD.f"
    /* Function Body */
#line 244 "SB10LD.f"
    n2 = *n << 1;
#line 245 "SB10LD.f"
    m1 = *m - *ncon;
#line 246 "SB10LD.f"
    m2 = *ncon;
#line 247 "SB10LD.f"
    np1 = *np - *nmeas;
#line 248 "SB10LD.f"
    np2 = *nmeas;

#line 250 "SB10LD.f"
    *info = 0;
#line 251 "SB10LD.f"
    if (*n < 0) {
#line 252 "SB10LD.f"
	*info = -1;
#line 253 "SB10LD.f"
    } else if (*m < 0) {
#line 254 "SB10LD.f"
	*info = -2;
#line 255 "SB10LD.f"
    } else if (*np < 0) {
#line 256 "SB10LD.f"
	*info = -3;
#line 257 "SB10LD.f"
    } else if (*ncon < 0 || m1 < 0 || m2 > np1) {
#line 258 "SB10LD.f"
	*info = -4;
#line 259 "SB10LD.f"
    } else if (*nmeas < 0 || np1 < 0 || np2 > m1) {
#line 260 "SB10LD.f"
	*info = -5;
#line 261 "SB10LD.f"
    } else if (*lda < max(1,*n)) {
#line 262 "SB10LD.f"
	*info = -7;
#line 263 "SB10LD.f"
    } else if (*ldb < max(1,*n)) {
#line 264 "SB10LD.f"
	*info = -9;
#line 265 "SB10LD.f"
    } else if (*ldc < max(1,*np)) {
#line 266 "SB10LD.f"
	*info = -11;
#line 267 "SB10LD.f"
    } else if (*ldd < max(1,*np)) {
#line 268 "SB10LD.f"
	*info = -13;
#line 269 "SB10LD.f"
    } else if (*ldak < max(1,*n)) {
#line 270 "SB10LD.f"
	*info = -15;
#line 271 "SB10LD.f"
    } else if (*ldbk < max(1,*n)) {
#line 272 "SB10LD.f"
	*info = -17;
#line 273 "SB10LD.f"
    } else if (*ldck < max(1,m2)) {
#line 274 "SB10LD.f"
	*info = -19;
#line 275 "SB10LD.f"
    } else if (*lddk < max(1,m2)) {
#line 276 "SB10LD.f"
	*info = -21;
#line 277 "SB10LD.f"
    } else if (*ldac < max(1,n2)) {
#line 278 "SB10LD.f"
	*info = -23;
#line 279 "SB10LD.f"
    } else if (*ldbc < max(1,n2)) {
#line 280 "SB10LD.f"
	*info = -25;
#line 281 "SB10LD.f"
    } else if (*ldcc < max(1,np1)) {
#line 282 "SB10LD.f"
	*info = -27;
#line 283 "SB10LD.f"
    } else if (*lddc < max(1,np1)) {
#line 284 "SB10LD.f"
	*info = -29;
#line 285 "SB10LD.f"
    } else {

/*        Compute workspace. */

#line 289 "SB10LD.f"
	minwrk = (*m << 1) * *m + *np * *np + (*m << 1) * *n + *m * *np + (*n 
		<< 1) * *np;
#line 290 "SB10LD.f"
	if (*ldwork < minwrk) {
#line 290 "SB10LD.f"
	    *info = -32;
#line 290 "SB10LD.f"
	}
#line 292 "SB10LD.f"
    }
#line 293 "SB10LD.f"
    if (*info != 0) {
#line 294 "SB10LD.f"
	i__1 = -(*info);
#line 294 "SB10LD.f"
	xerbla_("SB10LD", &i__1, (ftnlen)6);
#line 295 "SB10LD.f"
	return 0;
#line 296 "SB10LD.f"
    }

/*     Quick return if possible. */

#line 300 "SB10LD.f"
    if (*n == 0 || *m == 0 || *np == 0 || m1 == 0 || m2 == 0 || np1 == 0 || 
	    np2 == 0) {
#line 302 "SB10LD.f"
	dwork[1] = 1.;
#line 303 "SB10LD.f"
	return 0;
#line 304 "SB10LD.f"
    }

/*     Get the machine precision. */

#line 308 "SB10LD.f"
    eps = dlamch_("Epsilon", (ftnlen)7);

/*     Workspace usage. */

#line 312 "SB10LD.f"
    iw2 = np2 * np2 + 1;
#line 313 "SB10LD.f"
    iw3 = iw2 + m2 * m2;
#line 314 "SB10LD.f"
    iw4 = iw3 + np2 * *n;
#line 315 "SB10LD.f"
    iw5 = iw4 + m2 * *n;
#line 316 "SB10LD.f"
    iw6 = iw5 + np2 * m1;
#line 317 "SB10LD.f"
    iw7 = iw6 + m2 * m1;
#line 318 "SB10LD.f"
    iw8 = iw7 + m2 * *n;
#line 319 "SB10LD.f"
    iwrk = iw8 + np2 * *n;

/*     Compute inv(Inp2 - D22*DK) . */

#line 323 "SB10LD.f"
    dlaset_("Full", &np2, &np2, &c_b5, &c_b6, &dwork[1], &np2, (ftnlen)4);
#line 324 "SB10LD.f"
    dgemm_("N", "N", &np2, &np2, &m2, &c_b9, &d__[np1 + 1 + (m1 + 1) * d_dim1]
	    , ldd, &dk[dk_offset], lddk, &c_b6, &dwork[1], &np2, (ftnlen)1, (
	    ftnlen)1);
#line 326 "SB10LD.f"
    anorm = dlange_("1", &np2, &np2, &dwork[1], &np2, &dwork[iwrk], (ftnlen)1)
	    ;
#line 327 "SB10LD.f"
    dgetrf_(&np2, &np2, &dwork[1], &np2, &iwork[1], &info2);
#line 328 "SB10LD.f"
    if (info2 > 0) {
#line 329 "SB10LD.f"
	*info = 1;
#line 330 "SB10LD.f"
	return 0;
#line 331 "SB10LD.f"
    }
#line 332 "SB10LD.f"
    dgecon_("1", &np2, &dwork[1], &np2, &anorm, &rcond, &dwork[iwrk], &iwork[
	    np2 + 1], info, (ftnlen)1);
#line 334 "SB10LD.f"
    lwamax = (integer) dwork[iwrk] + iwrk - 1;

/*     Return if the matrix is singular to working precision. */

#line 338 "SB10LD.f"
    if (rcond < eps) {
#line 339 "SB10LD.f"
	*info = 1;
#line 340 "SB10LD.f"
	return 0;
#line 341 "SB10LD.f"
    }
#line 342 "SB10LD.f"
    i__1 = *ldwork - iwrk + 1;
#line 342 "SB10LD.f"
    dgetri_(&np2, &dwork[1], &np2, &iwork[1], &dwork[iwrk], &i__1, &info2);
/* Computing MAX */
#line 344 "SB10LD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 344 "SB10LD.f"
    lwamax = max(i__1,lwamax);

/*     Compute inv(Im2 - DK*D22) . */

#line 348 "SB10LD.f"
    dlaset_("Full", &m2, &m2, &c_b5, &c_b6, &dwork[iw2], &m2, (ftnlen)4);
#line 349 "SB10LD.f"
    dgemm_("N", "N", &m2, &m2, &np2, &c_b9, &dk[dk_offset], lddk, &d__[np1 + 
	    1 + (m1 + 1) * d_dim1], ldd, &c_b6, &dwork[iw2], &m2, (ftnlen)1, (
	    ftnlen)1);
#line 351 "SB10LD.f"
    anorm = dlange_("1", &m2, &m2, &dwork[iw2], &m2, &dwork[iwrk], (ftnlen)1);
#line 352 "SB10LD.f"
    dgetrf_(&m2, &m2, &dwork[iw2], &m2, &iwork[1], &info2);
#line 353 "SB10LD.f"
    if (info2 > 0) {
#line 354 "SB10LD.f"
	*info = 2;
#line 355 "SB10LD.f"
	return 0;
#line 356 "SB10LD.f"
    }
#line 357 "SB10LD.f"
    dgecon_("1", &m2, &dwork[iw2], &m2, &anorm, &rcond, &dwork[iwrk], &iwork[
	    m2 + 1], info, (ftnlen)1);
/* Computing MAX */
#line 359 "SB10LD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 359 "SB10LD.f"
    lwamax = max(i__1,lwamax);

/*     Return if the matrix is singular to working precision. */

#line 363 "SB10LD.f"
    if (rcond < eps) {
#line 364 "SB10LD.f"
	*info = 2;
#line 365 "SB10LD.f"
	return 0;
#line 366 "SB10LD.f"
    }
#line 367 "SB10LD.f"
    i__1 = *ldwork - iwrk + 1;
#line 367 "SB10LD.f"
    dgetri_(&m2, &dwork[iw2], &m2, &iwork[1], &dwork[iwrk], &i__1, &info2);
/* Computing MAX */
#line 369 "SB10LD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 369 "SB10LD.f"
    lwamax = max(i__1,lwamax);

/*     Compute inv(Inp2 - D22*DK)*C2 . */

#line 373 "SB10LD.f"
    dgemm_("N", "N", &np2, n, &np2, &c_b6, &dwork[1], &np2, &c__[np1 + 1 + 
	    c_dim1], ldc, &c_b5, &dwork[iw3], &np2, (ftnlen)1, (ftnlen)1);

/*     Compute DK*inv(Inp2 - D22*DK)*C2 . */

#line 378 "SB10LD.f"
    dgemm_("N", "N", &m2, n, &np2, &c_b6, &dk[dk_offset], lddk, &dwork[iw3], &
	    np2, &c_b5, &dwork[iw4], &m2, (ftnlen)1, (ftnlen)1);

/*     Compute inv(Inp2 - D22*DK)*D21 . */

#line 383 "SB10LD.f"
    dgemm_("N", "N", &np2, &m1, &np2, &c_b6, &dwork[1], &np2, &d__[np1 + 1 + 
	    d_dim1], ldd, &c_b5, &dwork[iw5], &np2, (ftnlen)1, (ftnlen)1);

/*     Compute DK*inv(Inp2 - D22*DK)*D21 . */

#line 388 "SB10LD.f"
    dgemm_("N", "N", &m2, &m1, &np2, &c_b6, &dk[dk_offset], lddk, &dwork[iw5],
	     &np2, &c_b5, &dwork[iw6], &m2, (ftnlen)1, (ftnlen)1);

/*     Compute inv(Im2 - DK*D22)*CK . */

#line 393 "SB10LD.f"
    dgemm_("N", "N", &m2, n, &m2, &c_b6, &dwork[iw2], &m2, &ck[ck_offset], 
	    ldck, &c_b5, &dwork[iw7], &m2, (ftnlen)1, (ftnlen)1);

/*     Compute D22*inv(Im2 - DK*D22)*CK . */

#line 398 "SB10LD.f"
    dgemm_("N", "N", &np2, n, &m2, &c_b6, &d__[np1 + 1 + (m1 + 1) * d_dim1], 
	    ldd, &dwork[iw7], &m2, &c_b5, &dwork[iw8], &np2, (ftnlen)1, (
	    ftnlen)1);

/*     Compute AC . */

#line 403 "SB10LD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &ac[ac_offset], ldac, (ftnlen)4);
#line 404 "SB10LD.f"
    dgemm_("N", "N", n, n, &m2, &c_b6, &b[(m1 + 1) * b_dim1 + 1], ldb, &dwork[
	    iw4], &m2, &c_b6, &ac[ac_offset], ldac, (ftnlen)1, (ftnlen)1);
#line 406 "SB10LD.f"
    dgemm_("N", "N", n, n, &m2, &c_b6, &b[(m1 + 1) * b_dim1 + 1], ldb, &dwork[
	    iw7], &m2, &c_b5, &ac[(*n + 1) * ac_dim1 + 1], ldac, (ftnlen)1, (
	    ftnlen)1);
#line 408 "SB10LD.f"
    dgemm_("N", "N", n, n, &np2, &c_b6, &bk[bk_offset], ldbk, &dwork[iw3], &
	    np2, &c_b5, &ac[*n + 1 + ac_dim1], ldac, (ftnlen)1, (ftnlen)1);
#line 410 "SB10LD.f"
    dlacpy_("Full", n, n, &ak[ak_offset], ldak, &ac[*n + 1 + (*n + 1) * 
	    ac_dim1], ldac, (ftnlen)4);
#line 411 "SB10LD.f"
    dgemm_("N", "N", n, n, &np2, &c_b6, &bk[bk_offset], ldbk, &dwork[iw8], &
	    np2, &c_b6, &ac[*n + 1 + (*n + 1) * ac_dim1], ldac, (ftnlen)1, (
	    ftnlen)1);

/*     Compute BC . */

#line 416 "SB10LD.f"
    dlacpy_("Full", n, &m1, &b[b_offset], ldb, &bc[bc_offset], ldbc, (ftnlen)
	    4);
#line 417 "SB10LD.f"
    dgemm_("N", "N", n, &m1, &m2, &c_b6, &b[(m1 + 1) * b_dim1 + 1], ldb, &
	    dwork[iw6], &m2, &c_b6, &bc[bc_offset], ldbc, (ftnlen)1, (ftnlen)
	    1);
#line 419 "SB10LD.f"
    dgemm_("N", "N", n, &m1, &np2, &c_b6, &bk[bk_offset], ldbk, &dwork[iw5], &
	    np2, &c_b5, &bc[*n + 1 + bc_dim1], ldbc, (ftnlen)1, (ftnlen)1);

/*     Compute CC . */

#line 424 "SB10LD.f"
    dlacpy_("Full", &np1, n, &c__[c_offset], ldc, &cc[cc_offset], ldcc, (
	    ftnlen)4);
#line 425 "SB10LD.f"
    dgemm_("N", "N", &np1, n, &m2, &c_b6, &d__[(m1 + 1) * d_dim1 + 1], ldd, &
	    dwork[iw4], &m2, &c_b6, &cc[cc_offset], ldcc, (ftnlen)1, (ftnlen)
	    1);
#line 427 "SB10LD.f"
    dgemm_("N", "N", &np1, n, &m2, &c_b6, &d__[(m1 + 1) * d_dim1 + 1], ldd, &
	    dwork[iw7], &m2, &c_b5, &cc[(*n + 1) * cc_dim1 + 1], ldcc, (
	    ftnlen)1, (ftnlen)1);

/*     Compute DC . */

#line 432 "SB10LD.f"
    dlacpy_("Full", &np1, &m1, &d__[d_offset], ldd, &dc[dc_offset], lddc, (
	    ftnlen)4);
#line 433 "SB10LD.f"
    dgemm_("N", "N", &np1, &m1, &m2, &c_b6, &d__[(m1 + 1) * d_dim1 + 1], ldd, 
	    &dwork[iw6], &m2, &c_b6, &dc[dc_offset], lddc, (ftnlen)1, (ftnlen)
	    1);

#line 436 "SB10LD.f"
    return 0;
/* *** Last line of SB10LD *** */
} /* sb10ld_ */

