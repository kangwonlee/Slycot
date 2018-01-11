#line 1 "AB08NX.f"
/* AB08NX.f -- translated by f2c (version 20100827).
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

#line 1 "AB08NX.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b22 = 0.;
static logical c_true = TRUE_;

/* Subroutine */ int ab08nx_(integer *n, integer *m, integer *p, integer *ro, 
	integer *sigma, doublereal *svlmax, doublereal *abcd, integer *ldabcd,
	 integer *ninfz, integer *infz, integer *kronl, integer *mu, integer *
	nu, integer *nkrol, doublereal *tol, integer *iwork, doublereal *
	dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer abcd_dim1, abcd_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static doublereal t;
    static integer i1, nb, ik, np, iz, mm1, ro1, mpm, tau, mnu, rank, itau;
    static doublereal sval[3];
    static integer irow;
    extern /* Subroutine */ int mb03oy_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *), mb03py_(
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *);
    static integer mntau, jwork;
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlapmt_(logical *, integer *, integer *, 
	    doublereal *, integer *, integer *), dlatzm_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, ftnlen), dormqr_(char *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), dormrq_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static logical lquery;
    static integer wrkopt;


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

/*     To extract from the (N+P)-by-(M+N) system */
/*                  ( B  A ) */
/*                  ( D  C ) */
/*     an (NU+MU)-by-(M+NU) "reduced" system */
/*                  ( B' A') */
/*                  ( D' C') */
/*     having the same transmission zeros but with D' of full row rank. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of state variables.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     RO      (input/output) INTEGER */
/*             On entry, */
/*             = P     for the original system; */
/*             = MAX(P-M, 0) for the pertransposed system. */
/*             On exit, RO contains the last computed rank. */

/*     SIGMA   (input/output) INTEGER */
/*             On entry, */
/*             = 0  for the original system; */
/*             = M  for the pertransposed system. */
/*             On exit, SIGMA contains the last computed value sigma in */
/*             the algorithm. */

/*     SVLMAX  (input) DOUBLE PRECISION */
/*             During each reduction step, the rank-revealing QR */
/*             factorization of a matrix stops when the estimated minimum */
/*             singular value is smaller than TOL * MAX(SVLMAX,EMSV), */
/*             where EMSV is the estimated maximum singular value. */
/*             SVLMAX >= 0. */

/*     ABCD    (input/output) DOUBLE PRECISION array, dimension */
/*             (LDABCD,M+N) */
/*             On entry, the leading (N+P)-by-(M+N) part of this array */
/*             must contain the compound input matrix of the system. */
/*             On exit, the leading (NU+MU)-by-(M+NU) part of this array */
/*             contains the reduced compound input matrix of the system. */

/*     LDABCD  INTEGER */
/*             The leading dimension of array ABCD. */
/*             LDABCD >= MAX(1,N+P). */

/*     NINFZ   (input/output) INTEGER */
/*             On entry, the currently computed number of infinite zeros. */
/*             It should be initialized to zero on the first call. */
/*             NINFZ >= 0. */
/*             On exit, the number of infinite zeros. */

/*     INFZ    (input/output) INTEGER array, dimension (N) */
/*             On entry, INFZ(i) must contain the current number of */
/*             infinite zeros of degree i, where i = 1,2,...,N, found in */
/*             the previous call(s) of the routine. It should be */
/*             initialized to zero on the first call. */
/*             On exit, INFZ(i) contains the number of infinite zeros of */
/*             degree i, where i = 1,2,...,N. */

/*     KRONL   (input/output) INTEGER array, dimension (N+1) */
/*             On entry, this array must contain the currently computed */
/*             left Kronecker (row) indices found in the previous call(s) */
/*             of the routine. It should be initialized to zero on the */
/*             first call. */
/*             On exit, the leading NKROL elements of this array contain */
/*             the left Kronecker (row) indices. */

/*     MU      (output) INTEGER */
/*             The normal rank of the transfer function matrix of the */
/*             original system. */

/*     NU      (output) INTEGER */
/*             The dimension of the reduced system matrix and the number */
/*             of (finite) invariant zeros if D' is invertible. */

/*     NKROL   (output) INTEGER */
/*             The number of left Kronecker indices. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             A tolerance used in rank decisions to determine the */
/*             effective rank, which is defined as the order of the */
/*             largest leading (or trailing) triangular submatrix in the */
/*             QR (or RQ) factorization with column (or row) pivoting */
/*             whose estimated condition number is less than 1/TOL. */
/*             NOTE that when SVLMAX > 0, the estimated ranks could be */
/*             less than those defined above (see SVLMAX). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (MAX(M,P)) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 1, MIN(P,M) + MAX(3*M-1,N), */
/*                               MIN(P,N) + MAX(3*P-1,N+P,N+M) ). */
/*             For optimum performance LDWORK should be larger. */

/*             If LDWORK = -1, then a workspace query is assumed; */
/*             the routine only calculates the optimal size of the */
/*             DWORK array, returns this value as the first entry of */
/*             the DWORK array, and no error message related to LDWORK */
/*             is issued by XERBLA. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     REFERENCES */

/*     [1] Svaricek, F. */
/*         Computation of the Structural Invariants of Linear */
/*         Multivariable Systems with an Extended Version of */
/*         the Program ZEROS. */
/*         System & Control Letters, 6, pp. 261-266, 1985. */

/*     [2] Emami-Naeini, A. and Van Dooren, P. */
/*         Computation of Zeros of Linear Multivariable Systems. */
/*         Automatica, 18, pp. 415-430, 1982. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996. */
/*     Supersedes Release 2.0 routine AB08BZ by F. Svaricek. */

/*     REVISIONS */

/*     V. Sima, Oct. 1997; Feb. 1998, Jan. 2009, Apr. 2009. */
/*     A. Varga, May 1999; May 2001. */

/*     KEYWORDS */

/*     Generalized eigenvalue problem, Kronecker indices, multivariable */
/*     system, orthogonal transformation, structural invariant. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 208 "AB08NX.f"
    /* Parameter adjustments */
#line 208 "AB08NX.f"
    abcd_dim1 = *ldabcd;
#line 208 "AB08NX.f"
    abcd_offset = 1 + abcd_dim1;
#line 208 "AB08NX.f"
    abcd -= abcd_offset;
#line 208 "AB08NX.f"
    --infz;
#line 208 "AB08NX.f"
    --kronl;
#line 208 "AB08NX.f"
    --iwork;
#line 208 "AB08NX.f"
    --dwork;
#line 208 "AB08NX.f"

#line 208 "AB08NX.f"
    /* Function Body */
#line 208 "AB08NX.f"
    np = *n + *p;
#line 209 "AB08NX.f"
    mpm = min(*p,*m);
#line 210 "AB08NX.f"
    *info = 0;
#line 211 "AB08NX.f"
    lquery = *ldwork == -1;

/*     Test the input scalar arguments. */

#line 215 "AB08NX.f"
    if (*n < 0) {
#line 216 "AB08NX.f"
	*info = -1;
#line 217 "AB08NX.f"
    } else if (*m < 0) {
#line 218 "AB08NX.f"
	*info = -2;
#line 219 "AB08NX.f"
    } else if (*p < 0) {
#line 220 "AB08NX.f"
	*info = -3;
#line 221 "AB08NX.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 221 "AB08NX.f"
	i__1 = *p - *m;
#line 221 "AB08NX.f"
	if (*ro != *p && *ro != max(i__1,0)) {
#line 222 "AB08NX.f"
	    *info = -4;
#line 223 "AB08NX.f"
	} else if (*sigma != 0 && *sigma != *m) {
#line 224 "AB08NX.f"
	    *info = -5;
#line 225 "AB08NX.f"
	} else if (*svlmax < 0.) {
#line 226 "AB08NX.f"
	    *info = -6;
#line 227 "AB08NX.f"
	} else if (*ldabcd < max(1,np)) {
#line 228 "AB08NX.f"
	    *info = -8;
#line 229 "AB08NX.f"
	} else if (*ninfz < 0) {
#line 230 "AB08NX.f"
	    *info = -9;
#line 231 "AB08NX.f"
	} else {
/* Computing MAX */
/* Computing MAX */
#line 232 "AB08NX.f"
	    i__3 = *m * 3 - 1;
/* Computing MAX */
#line 232 "AB08NX.f"
	    i__4 = *p * 3 - 1, i__4 = max(i__4,np), i__5 = *n + *m;
#line 232 "AB08NX.f"
	    i__1 = 1, i__2 = mpm + max(i__3,*n), i__1 = max(i__1,i__2), i__2 =
		     min(*p,*n) + max(i__4,i__5);
#line 232 "AB08NX.f"
	    jwork = max(i__1,i__2);
#line 234 "AB08NX.f"
	    if (lquery) {
#line 235 "AB08NX.f"
		if (*m > 0) {
/* Computing MIN */
#line 236 "AB08NX.f"
		    i__1 = 64, i__2 = ilaenv_(&c__1, "DORMQR", "LT", p, n, &
			    mpm, &c_n1, (ftnlen)6, (ftnlen)2);
#line 236 "AB08NX.f"
		    nb = min(i__1,i__2);
/* Computing MAX */
#line 238 "AB08NX.f"
		    i__1 = jwork, i__2 = mpm + max(1,*n) * nb;
#line 238 "AB08NX.f"
		    wrkopt = max(i__1,i__2);
#line 239 "AB08NX.f"
		} else {
#line 240 "AB08NX.f"
		    wrkopt = jwork;
#line 241 "AB08NX.f"
		}
/* Computing MIN */
#line 242 "AB08NX.f"
		i__3 = min(*p,*n);
#line 242 "AB08NX.f"
		i__1 = 64, i__2 = ilaenv_(&c__1, "DORMRQ", "RT", &np, n, &
			i__3, &c_n1, (ftnlen)6, (ftnlen)2);
#line 242 "AB08NX.f"
		nb = min(i__1,i__2);
/* Computing MAX */
#line 244 "AB08NX.f"
		i__1 = wrkopt, i__2 = min(*p,*n) + max(1,np) * nb;
#line 244 "AB08NX.f"
		wrkopt = max(i__1,i__2);
/* Computing MIN */
#line 245 "AB08NX.f"
		i__3 = *m + *n;
#line 245 "AB08NX.f"
		i__4 = min(*p,*n);
#line 245 "AB08NX.f"
		i__1 = 64, i__2 = ilaenv_(&c__1, "DORMRQ", "LN", n, &i__3, &
			i__4, &c_n1, (ftnlen)6, (ftnlen)2);
#line 245 "AB08NX.f"
		nb = min(i__1,i__2);
/* Computing MAX */
/* Computing MAX */
#line 247 "AB08NX.f"
		i__3 = 1, i__4 = *m + *n;
#line 247 "AB08NX.f"
		i__1 = wrkopt, i__2 = min(*p,*n) + max(i__3,i__4) * nb;
#line 247 "AB08NX.f"
		wrkopt = max(i__1,i__2);
#line 248 "AB08NX.f"
	    } else if (*ldwork < jwork) {
#line 249 "AB08NX.f"
		*info = -18;
#line 250 "AB08NX.f"
	    }
#line 251 "AB08NX.f"
	}
#line 251 "AB08NX.f"
    }

#line 253 "AB08NX.f"
    if (*info != 0) {

/*        Error return. */

#line 257 "AB08NX.f"
	i__1 = -(*info);
#line 257 "AB08NX.f"
	xerbla_("AB08NX", &i__1, (ftnlen)6);
#line 258 "AB08NX.f"
	return 0;
#line 259 "AB08NX.f"
    } else if (lquery) {
#line 260 "AB08NX.f"
	dwork[1] = (doublereal) wrkopt;
#line 261 "AB08NX.f"
	return 0;
#line 262 "AB08NX.f"
    }

#line 264 "AB08NX.f"
    *mu = *p;
#line 265 "AB08NX.f"
    *nu = *n;

#line 267 "AB08NX.f"
    iz = 0;
#line 268 "AB08NX.f"
    ik = 1;
#line 269 "AB08NX.f"
    mm1 = *m + 1;
#line 270 "AB08NX.f"
    itau = 1;
#line 271 "AB08NX.f"
    *nkrol = 0;
#line 272 "AB08NX.f"
    wrkopt = 1;

/*     Main reduction loop: */

/*            M   NU                  M     NU */
/*      NU  [ B   A ]           NU  [ B     A ] */
/*      MU  [ D   C ]  -->    SIGMA [ RD   C1 ]   (SIGMA = rank(D) = */
/*                             TAU  [ 0    C2 ]    row size of RD) */

/*                                    M   NU-RO  RO */
/*                            NU-RO [ B1   A11  A12 ] */
/*                     -->      RO  [ B2   A21  A22 ]  (RO = rank(C2) = */
/*                            SIGMA [ RD   C11  C12 ]   col size of LC) */
/*                             TAU  [ 0     0   LC  ] */

/*                                      M   NU-RO */
/*                            NU-RO [ B1   A11 ]     NU := NU - RO */
/*                                  [----------]     MU := RO + SIGMA */
/*                     -->      RO  [ B2   A21 ]      D := [B2;RD] */
/*                            SIGMA [ RD   C11 ]      C := [A21;C11] */

#line 293 "AB08NX.f"
L20:
#line 293 "AB08NX.f"
    if (*mu == 0) {
#line 293 "AB08NX.f"
	goto L80;
#line 293 "AB08NX.f"
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance.) */

#line 300 "AB08NX.f"
    ro1 = *ro;
#line 301 "AB08NX.f"
    mnu = *m + *nu;
#line 302 "AB08NX.f"
    if (*m > 0) {
#line 303 "AB08NX.f"
	if (*sigma != 0) {
#line 304 "AB08NX.f"
	    irow = *nu + 1;

/*           Compress rows of D.  First exploit triangular shape. */
/*           Workspace: need   M+N-1. */

#line 309 "AB08NX.f"
	    i__1 = *sigma;
#line 309 "AB08NX.f"
	    for (i1 = 1; i1 <= i__1; ++i1) {
#line 310 "AB08NX.f"
		i__2 = *ro + 1;
#line 310 "AB08NX.f"
		dlarfg_(&i__2, &abcd[irow + i1 * abcd_dim1], &abcd[irow + 1 + 
			i1 * abcd_dim1], &c__1, &t);
#line 311 "AB08NX.f"
		i__2 = *ro + 1;
#line 311 "AB08NX.f"
		i__3 = mnu - i1;
#line 311 "AB08NX.f"
		dlatzm_("L", &i__2, &i__3, &abcd[irow + 1 + i1 * abcd_dim1], &
			c__1, &t, &abcd[irow + (i1 + 1) * abcd_dim1], &abcd[
			irow + 1 + (i1 + 1) * abcd_dim1], ldabcd, &dwork[1], (
			ftnlen)1);
#line 314 "AB08NX.f"
		++irow;
#line 315 "AB08NX.f"
/* L40: */
#line 315 "AB08NX.f"
	    }
#line 316 "AB08NX.f"
	    i__1 = *ro + *sigma - 1;
#line 316 "AB08NX.f"
	    dlaset_("Lower", &i__1, sigma, &c_b22, &c_b22, &abcd[*nu + 2 + 
		    abcd_dim1], ldabcd, (ftnlen)5);
#line 318 "AB08NX.f"
	}

/*        Continue with Householder with column pivoting. */

/*        The rank of D is the number of (estimated) singular values */
/*        that are greater than TOL * MAX(SVLMAX,EMSV). This number */
/*        includes the singular values of the first SIGMA columns. */
/*        Integer workspace: need   M; */
/*        Workspace: need   min(RO1,M) + 3*M - 1.  RO1 <= P. */

#line 328 "AB08NX.f"
	if (*sigma < *m) {
#line 329 "AB08NX.f"
	    jwork = itau + min(ro1,*m);
#line 330 "AB08NX.f"
	    i1 = *sigma + 1;
#line 331 "AB08NX.f"
	    irow = *nu + i1;
#line 332 "AB08NX.f"
	    i__1 = *m - *sigma;
#line 332 "AB08NX.f"
	    mb03oy_(&ro1, &i__1, &abcd[irow + i1 * abcd_dim1], ldabcd, tol, 
		    svlmax, &rank, sval, &iwork[1], &dwork[itau], &dwork[
		    jwork], info);
/* Computing MAX */
#line 335 "AB08NX.f"
	    i__1 = wrkopt, i__2 = jwork + *m * 3 - 2;
#line 335 "AB08NX.f"
	    wrkopt = max(i__1,i__2);

/*           Apply the column permutations to matrices B and part of D. */

#line 339 "AB08NX.f"
	    i__1 = *nu + *sigma;
#line 339 "AB08NX.f"
	    i__2 = *m - *sigma;
#line 339 "AB08NX.f"
	    dlapmt_(&c_true, &i__1, &i__2, &abcd[i1 * abcd_dim1 + 1], ldabcd, 
		    &iwork[1]);

#line 342 "AB08NX.f"
	    if (rank > 0) {

/*              Apply the Householder transformations to the submatrix C. */
/*              Workspace: need   min(RO1,M) + NU; */
/*                         prefer min(RO1,M) + NU*NB. */

#line 348 "AB08NX.f"
		i__1 = *ldwork - jwork + 1;
#line 348 "AB08NX.f"
		dormqr_("Left", "Transpose", &ro1, nu, &rank, &abcd[irow + i1 
			* abcd_dim1], ldabcd, &dwork[itau], &abcd[irow + mm1 *
			 abcd_dim1], ldabcd, &dwork[jwork], &i__1, info, (
			ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 352 "AB08NX.f"
		i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 352 "AB08NX.f"
		wrkopt = max(i__1,i__2);
#line 353 "AB08NX.f"
		if (ro1 > 1) {
#line 353 "AB08NX.f"
		    i__1 = ro1 - 1;
/* Computing MIN */
#line 353 "AB08NX.f"
		    i__3 = ro1 - 1;
#line 353 "AB08NX.f"
		    i__2 = min(i__3,rank);
#line 353 "AB08NX.f"
		    dlaset_("Lower", &i__1, &i__2, &c_b22, &c_b22, &abcd[irow 
			    + 1 + i1 * abcd_dim1], ldabcd, (ftnlen)5);
#line 353 "AB08NX.f"
		}
#line 356 "AB08NX.f"
		ro1 -= rank;
#line 357 "AB08NX.f"
	    }
#line 358 "AB08NX.f"
	}
#line 359 "AB08NX.f"
    }

#line 361 "AB08NX.f"
    tau = ro1;
#line 362 "AB08NX.f"
    *sigma = *mu - tau;

/*     Determination of the orders of the infinite zeros. */

#line 366 "AB08NX.f"
    if (iz > 0) {
#line 367 "AB08NX.f"
	infz[iz] = infz[iz] + *ro - tau;
#line 368 "AB08NX.f"
	*ninfz += iz * (*ro - tau);
#line 369 "AB08NX.f"
    }
#line 370 "AB08NX.f"
    if (ro1 == 0) {
#line 370 "AB08NX.f"
	goto L80;
#line 370 "AB08NX.f"
    }
#line 372 "AB08NX.f"
    ++iz;

#line 374 "AB08NX.f"
    if (*nu <= 0) {
#line 375 "AB08NX.f"
	*mu = *sigma;
#line 376 "AB08NX.f"
	*nu = 0;
#line 377 "AB08NX.f"
	*ro = 0;
#line 378 "AB08NX.f"
    } else {

/*        Compress the columns of C2 using RQ factorization with row */
/*        pivoting, P * C2 = R * Q. */

#line 383 "AB08NX.f"
	i1 = *nu + *sigma + 1;
#line 384 "AB08NX.f"
	mntau = min(tau,*nu);
#line 385 "AB08NX.f"
	jwork = itau + mntau;

/*        The rank of C2 is the number of (estimated) singular values */
/*        greater than TOL * MAX(SVLMAX,EMSV). */
/*        Integer Workspace: need TAU; */
/*        Workspace: need min(TAU,NU) + 3*TAU - 1. */

#line 392 "AB08NX.f"
	mb03py_(&tau, nu, &abcd[i1 + mm1 * abcd_dim1], ldabcd, tol, svlmax, &
		rank, sval, &iwork[1], &dwork[itau], &dwork[jwork], info);
/* Computing MAX */
#line 394 "AB08NX.f"
	i__1 = wrkopt, i__2 = jwork + tau * 3 - 1;
#line 394 "AB08NX.f"
	wrkopt = max(i__1,i__2);
#line 395 "AB08NX.f"
	if (rank > 0) {
#line 396 "AB08NX.f"
	    irow = i1 + tau - rank;

/*           Apply Q' to the first NU columns of [A; C1] from the right. */
/*           Workspace: need   min(TAU,NU) + NU + SIGMA; SIGMA <= P; */
/*                      prefer min(TAU,NU) + (NU  + SIGMA)*NB. */

#line 402 "AB08NX.f"
	    i__1 = i1 - 1;
#line 402 "AB08NX.f"
	    i__2 = *ldwork - jwork + 1;
#line 402 "AB08NX.f"
	    dormrq_("Right", "Transpose", &i__1, nu, &rank, &abcd[irow + mm1 *
		     abcd_dim1], ldabcd, &dwork[mntau - rank + 1], &abcd[mm1 *
		     abcd_dim1 + 1], ldabcd, &dwork[jwork], &i__2, info, (
		    ftnlen)5, (ftnlen)9);
/* Computing MAX */
#line 406 "AB08NX.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 406 "AB08NX.f"
	    wrkopt = max(i__1,i__2);

/*           Apply Q to the first NU rows and M + NU columns of [ B  A ] */
/*           from the left. */
/*           Workspace: need   min(TAU,NU) + M + NU; */
/*                      prefer min(TAU,NU) + (M + NU)*NB. */

#line 413 "AB08NX.f"
	    i__1 = *ldwork - jwork + 1;
#line 413 "AB08NX.f"
	    dormrq_("Left", "NoTranspose", nu, &mnu, &rank, &abcd[irow + mm1 *
		     abcd_dim1], ldabcd, &dwork[mntau - rank + 1], &abcd[
		    abcd_offset], ldabcd, &dwork[jwork], &i__1, info, (ftnlen)
		    4, (ftnlen)11);
/* Computing MAX */
#line 417 "AB08NX.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 417 "AB08NX.f"
	    wrkopt = max(i__1,i__2);

#line 419 "AB08NX.f"
	    i__1 = *nu - rank;
#line 419 "AB08NX.f"
	    dlaset_("Full", &rank, &i__1, &c_b22, &c_b22, &abcd[irow + mm1 * 
		    abcd_dim1], ldabcd, (ftnlen)4);
#line 421 "AB08NX.f"
	    if (rank > 1) {
#line 421 "AB08NX.f"
		i__1 = rank - 1;
#line 421 "AB08NX.f"
		i__2 = rank - 1;
#line 421 "AB08NX.f"
		dlaset_("Lower", &i__1, &i__2, &c_b22, &c_b22, &abcd[irow + 1 
			+ (mm1 + *nu - rank) * abcd_dim1], ldabcd, (ftnlen)5);
#line 421 "AB08NX.f"
	    }
#line 424 "AB08NX.f"
	}

#line 426 "AB08NX.f"
	*ro = rank;
#line 427 "AB08NX.f"
    }

/*     Determine the left Kronecker indices (row indices). */

#line 431 "AB08NX.f"
    kronl[ik] = kronl[ik] + tau - *ro;
#line 432 "AB08NX.f"
    *nkrol += kronl[ik];
#line 433 "AB08NX.f"
    ++ik;

/*     C and D are updated to [A21 ; C11] and [B2 ; RD]. */

#line 437 "AB08NX.f"
    *nu -= *ro;
#line 438 "AB08NX.f"
    *mu = *sigma + *ro;
#line 439 "AB08NX.f"
    if (*ro != 0) {
#line 439 "AB08NX.f"
	goto L20;
#line 439 "AB08NX.f"
    }

#line 442 "AB08NX.f"
L80:
#line 443 "AB08NX.f"
    dwork[1] = (doublereal) wrkopt;
#line 444 "AB08NX.f"
    return 0;
/* *** Last line of AB08NX *** */
} /* ab08nx_ */

