#line 1 "AG08BY.f"
/* AG08BY.f -- translated by f2c (version 20100827).
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

#line 1 "AG08BY.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static logical c_true = TRUE_;
static doublereal c_b20 = 0.;
static integer c__2 = 2;

/* Subroutine */ int ag08by_(logical *first, integer *n, integer *m, integer *
	p, doublereal *svlmax, doublereal *abcd, integer *ldabcd, doublereal *
	e, integer *lde, integer *nr, integer *pr, integer *ninfz, integer *
	dinfz, integer *nkronl, integer *infz, integer *kronl, doublereal *
	tol, integer *iwork, doublereal *dwork, integer *ldwork, integer *
	info)
{
    /* System generated locals */
    integer abcd_dim1, abcd_offset, e_dim1, e_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c__;
    static integer i__, j, k;
    static doublereal s, t, c1, c2;
    static integer n1;
    static doublereal s1, s2;
    static integer nb, mn, pn, ro;
    static doublereal tt;
    static integer mn1, mp1, ro1, irc;
    static doublereal dum[1];
    static integer mpm, mui, mnr, icol, rank, itau, taui;
    static doublereal sval[3], smin, smax;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer irow;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static integer muim1, sigma;
    static doublereal rcond;
    extern /* Subroutine */ int mb03oy_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *);
    static integer ilast, jlast;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer ismin;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer ismax, mntau;
    extern /* Subroutine */ int dlaic1_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer jwork1, jwork2;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);
    static integer nblcks;
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dlapmt_(logical *, integer *, 
	    integer *, doublereal *, integer *, integer *), dlatzm_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, ftnlen), 
	    dormqr_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen);
    static doublereal sminpr, smaxpr;
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

/*     To extract from the (N+P)-by-(M+N) descriptor system pencil */

/*        S(lambda) = ( B   A - lambda*E  ) */
/*                    ( D        C        ) */

/*     with E nonsingular and upper triangular a */
/*     (NR+PR)-by-(M+NR) "reduced" descriptor system pencil */

/*                           ( Br  Ar-lambda*Er ) */
/*              Sr(lambda) = (                  ) */
/*                           ( Dr     Cr        ) */

/*     having the same finite Smith zeros as the pencil */
/*     S(lambda) but with Dr, a PR-by-M full row rank */
/*     left upper trapezoidal matrix, and Er, an NR-by-NR */
/*     upper triangular nonsingular matrix. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     FIRST   LOGICAL */
/*             Specifies if AG08BY is called first time or it is called */
/*             for an already reduced system, with D full column rank */
/*             with the last M rows in upper triangular form: */
/*             FIRST = .TRUE.,  first time called; */
/*             FIRST = .FALSE., not first time called. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of rows of matrix B, the number of columns of */
/*             matrix C and the order of square matrices A and E. */
/*             N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of matrices B and D.  M >= 0. */
/*             M <= P if FIRST = .FALSE. . */

/*     P       (input) INTEGER */
/*             The number of rows of matrices C and D.  P >= 0. */

/*     SVLMAX  (input) DOUBLE PRECISION */
/*             During each reduction step, the rank-revealing QR */
/*             factorization of a matrix stops when the estimated minimum */
/*             singular value is smaller than TOL * MAX(SVLMAX,EMSV), */
/*             where EMSV is the estimated maximum singular value. */
/*             SVLMAX >= 0. */

/*     ABCD    (input/output) DOUBLE PRECISION array, dimension */
/*             (LDABCD,M+N) */
/*             On entry, the leading (N+P)-by-(M+N) part of this array */
/*             must contain the compound matrix */
/*                      (  B   A  ) , */
/*                      (  D   C  ) */
/*             where A is an N-by-N matrix, B is an N-by-M matrix, */
/*             C is a P-by-N matrix and D is a P-by-M matrix. */
/*             If FIRST = .FALSE., then D must be a full column */
/*             rank matrix with the last M rows in upper triangular form. */
/*             On exit, the leading (NR+PR)-by-(M+NR) part of ABCD */
/*             contains the reduced compound matrix */
/*                       (  Br  Ar ) , */
/*                       (  Dr  Cr ) */
/*             where Ar is an NR-by-NR matrix, Br is an NR-by-M matrix, */
/*             Cr is a PR-by-NR matrix, Dr is a PR-by-M full row rank */
/*             left upper trapezoidal matrix with the first PR columns */
/*             in upper triangular form. */

/*     LDABCD  INTEGER */
/*             The leading dimension of array ABCD. */
/*             LDABCD >= MAX(1,N+P). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper triangular nonsingular matrix E. */
/*             On exit, the leading NR-by-NR part contains the reduced */
/*             upper triangular nonsingular matrix Er. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,N). */

/*     NR      (output) INTEGER */
/*             The order of the reduced matrices Ar and Er; also the */
/*             number of rows of the reduced matrix Br and the number */
/*             of columns of the reduced matrix Cr. */
/*             If Dr is invertible, NR is also the number of finite */
/*             Smith zeros. */

/*     PR      (output) INTEGER */
/*             The rank of the resulting matrix Dr; also the number of */
/*             rows of reduced matrices Cr and Dr. */

/*     NINFZ   (output) INTEGER */
/*             Number of infinite zeros.  NINFZ = 0 if FIRST = .FALSE. . */

/*     DINFZ   (output) INTEGER */
/*             The maximal multiplicity of infinite zeros. */
/*             DINFZ = 0 if FIRST = .FALSE. . */

/*     NKRONL  (output) INTEGER */
/*             The maximal dimension of left elementary Kronecker blocks. */

/*     INFZ    (output) INTEGER array, dimension (N) */
/*             INFZ(i) contains the number of infinite zeros of */
/*             degree i, where i = 1,2,...,DINFZ. */
/*             INFZ is not referenced if FIRST = .FALSE. . */

/*     KRONL   (output) INTEGER array, dimension (N+1) */
/*             KRONL(i) contains the number of left elementary Kronecker */
/*             blocks of dimension i-by-(i-1), where i = 1,2,...,NKRONL. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             A tolerance used in rank decisions to determine the */
/*             effective rank, which is defined as the order of the */
/*             largest leading (or trailing) triangular submatrix in the */
/*             QR (or RQ) factorization with column (or row) pivoting */
/*             whose estimated condition number is less than 1/TOL. */
/*             If the user sets TOL <= 0, then an implicitly computed, */
/*             default tolerance TOLDEF = (N+P)*(N+M)*EPS,  is used */
/*             instead, where EPS is the machine precision */
/*             (see LAPACK Library routine DLAMCH). */
/*             NOTE that when SVLMAX > 0, the estimated ranks could be */
/*             less than those defined above (see SVLMAX).  TOL <= 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M) */
/*             If FIRST = .FALSE., IWORK is not referenced. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 1, if P = 0; otherwise */
/*             LDWORK >= MAX( 1, N+M-1, MIN(P,M) + MAX(3*M-1,N), 5*P ), */
/*                                             if FIRST = .TRUE.; */
/*             LDWORK >= MAX( 1, N+M-1, 5*P ), if FIRST = .FALSE. . */
/*             The second term is not needed if M = 0. */
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

/*     METHOD */

/*     The subroutine is based on the reduction algorithm of [1]. */

/*     REFERENCES */

/*     [1] P. Misra, P. Van Dooren and A. Varga. */
/*         Computation of structural invariants of generalized */
/*         state-space systems. */
/*         Automatica, 30, pp. 1921-1936, 1994. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically backward stable and requires */
/*     0( (P+N)*(M+N)*N )  floating point operations. */

/*     FURTHER COMMENTS */

/*     The number of infinite zeros is computed as */

/*                   DINFZ */
/*        NINFZ =     Sum  (INFZ(i)*i) . */
/*                    i=1 */
/*     Note that each infinite zero of multiplicity k corresponds to */
/*     an infinite eigenvalue of multiplicity k+1. */
/*     The multiplicities of the infinite eigenvalues can be determined */
/*     from PR, DINFZ and INFZ(i), i = 1, ..., DINFZ, as follows: */

/*                     DINFZ */
/*     - there are PR - Sum (INFZ(i)) simple infinite eigenvalues; */
/*                      i=1 */

/*     - there are INFZ(i) infinite eigenvalues with multiplicity i+1, */
/*       for i = 1, ..., DINFZ. */

/*     The left Kronecker indices are: */

/*     [ 0  0 ...  0  | 1  1  ...  1 |  .... | NKRONL  ...  NKRONL ] */
/*     |<- KRONL(1) ->|<- KRONL(2) ->|       |<-  KRONL(NKRONL)  ->| */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     May 1999. Based on the RASP routine SRISEP. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Sep. 1999, */
/*     Jan. 2009, Apr. 2009. */
/*     A. Varga, DLR Oberpfaffenhofen, March 2002. */

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

/*     Test the input parameters. */

#line 276 "AG08BY.f"
    /* Parameter adjustments */
#line 276 "AG08BY.f"
    abcd_dim1 = *ldabcd;
#line 276 "AG08BY.f"
    abcd_offset = 1 + abcd_dim1;
#line 276 "AG08BY.f"
    abcd -= abcd_offset;
#line 276 "AG08BY.f"
    e_dim1 = *lde;
#line 276 "AG08BY.f"
    e_offset = 1 + e_dim1;
#line 276 "AG08BY.f"
    e -= e_offset;
#line 276 "AG08BY.f"
    --infz;
#line 276 "AG08BY.f"
    --kronl;
#line 276 "AG08BY.f"
    --iwork;
#line 276 "AG08BY.f"
    --dwork;
#line 276 "AG08BY.f"

#line 276 "AG08BY.f"
    /* Function Body */
#line 276 "AG08BY.f"
    lquery = *ldwork == -1;
#line 277 "AG08BY.f"
    *info = 0;
#line 278 "AG08BY.f"
    pn = *p + *n;
#line 279 "AG08BY.f"
    mn = *m + *n;
#line 280 "AG08BY.f"
    mpm = min(*p,*m);
#line 281 "AG08BY.f"
    if (*n < 0) {
#line 282 "AG08BY.f"
	*info = -2;
#line 283 "AG08BY.f"
    } else if (*m < 0 || ! (*first) && *m > *p) {
#line 284 "AG08BY.f"
	*info = -3;
#line 285 "AG08BY.f"
    } else if (*p < 0) {
#line 286 "AG08BY.f"
	*info = -4;
#line 287 "AG08BY.f"
    } else if (*svlmax < 0.) {
#line 288 "AG08BY.f"
	*info = -5;
#line 289 "AG08BY.f"
    } else if (*ldabcd < max(1,pn)) {
#line 290 "AG08BY.f"
	*info = -7;
#line 291 "AG08BY.f"
    } else if (*lde < max(1,*n)) {
#line 292 "AG08BY.f"
	*info = -9;
#line 293 "AG08BY.f"
    } else if (*tol > 1.) {
#line 294 "AG08BY.f"
	*info = -17;
#line 295 "AG08BY.f"
    } else {
/* Computing MAX */
#line 296 "AG08BY.f"
	i__1 = 1, i__2 = *p * 5;
#line 296 "AG08BY.f"
	wrkopt = max(i__1,i__2);
#line 297 "AG08BY.f"
	if (*p > 0) {
#line 298 "AG08BY.f"
	    if (*m > 0) {
/* Computing MAX */
#line 299 "AG08BY.f"
		i__1 = wrkopt, i__2 = mn - 1;
#line 299 "AG08BY.f"
		wrkopt = max(i__1,i__2);
#line 300 "AG08BY.f"
		if (*first) {
/* Computing MAX */
/* Computing MAX */
#line 301 "AG08BY.f"
		    i__3 = *m * 3 - 1;
#line 301 "AG08BY.f"
		    i__1 = wrkopt, i__2 = mpm + max(i__3,*n);
#line 301 "AG08BY.f"
		    wrkopt = max(i__1,i__2);
#line 302 "AG08BY.f"
		    if (lquery) {
/* Computing MIN */
#line 303 "AG08BY.f"
			i__1 = 64, i__2 = ilaenv_(&c__1, "DORMQR", "LT", p, n,
				 &mpm, &c_n1, (ftnlen)6, (ftnlen)2);
#line 303 "AG08BY.f"
			nb = min(i__1,i__2);
/* Computing MAX */
#line 305 "AG08BY.f"
			i__1 = wrkopt, i__2 = mpm + max(1,*n) * nb;
#line 305 "AG08BY.f"
			wrkopt = max(i__1,i__2);
#line 306 "AG08BY.f"
		    }
#line 307 "AG08BY.f"
		}
#line 308 "AG08BY.f"
	    }
#line 309 "AG08BY.f"
	}
#line 310 "AG08BY.f"
	if (*ldwork < wrkopt && ! lquery) {
#line 311 "AG08BY.f"
	    *info = -20;
#line 312 "AG08BY.f"
	}
#line 313 "AG08BY.f"
    }
#line 314 "AG08BY.f"
    if (*info != 0) {
#line 315 "AG08BY.f"
	i__1 = -(*info);
#line 315 "AG08BY.f"
	xerbla_("AG08BY", &i__1, (ftnlen)6);
#line 316 "AG08BY.f"
	return 0;
#line 317 "AG08BY.f"
    } else if (lquery) {
#line 318 "AG08BY.f"
	dwork[1] = (doublereal) wrkopt;
#line 319 "AG08BY.f"
	return 0;
#line 320 "AG08BY.f"
    }

/*     Initialize output variables. */

#line 324 "AG08BY.f"
    *pr = *p;
#line 325 "AG08BY.f"
    *nr = *n;
#line 326 "AG08BY.f"
    *dinfz = 0;
#line 327 "AG08BY.f"
    *ninfz = 0;
#line 328 "AG08BY.f"
    *nkronl = 0;

/*     Quick return if possible. */

#line 332 "AG08BY.f"
    if (*p == 0) {
#line 333 "AG08BY.f"
	dwork[1] = 1.;
#line 334 "AG08BY.f"
	return 0;
#line 335 "AG08BY.f"
    }
#line 336 "AG08BY.f"
    if (*n == 0 && *m == 0) {
#line 337 "AG08BY.f"
	*pr = 0;
#line 338 "AG08BY.f"
	*nkronl = 1;
#line 339 "AG08BY.f"
	kronl[1] = *p;
#line 340 "AG08BY.f"
	dwork[1] = 1.;
#line 341 "AG08BY.f"
	return 0;
#line 342 "AG08BY.f"
    }

#line 344 "AG08BY.f"
    rcond = *tol;
#line 345 "AG08BY.f"
    if (rcond <= 0.) {

/*        Use the default tolerance in rank determination. */

#line 349 "AG08BY.f"
	rcond = (doublereal) (pn * mn) * dlamch_("EPSILON", (ftnlen)7);
#line 350 "AG08BY.f"
    }

/*     The D matrix is (RO+SIGMA)-by-M, where RO = P - SIGMA and */
/*     SIGMA = 0 for FIRST = .TRUE. and SIGMA = M for FIRST = .FALSE.. */
/*     The leading (RO+SIGMA)-by-SIGMA submatrix of D has full column */
/*     rank, with the trailing SIGMA-by-SIGMA submatrix upper triangular. */

#line 357 "AG08BY.f"
    if (*first) {
#line 358 "AG08BY.f"
	sigma = 0;
#line 359 "AG08BY.f"
    } else {
#line 360 "AG08BY.f"
	sigma = *m;
#line 361 "AG08BY.f"
    }
#line 362 "AG08BY.f"
    ro = *p - sigma;
#line 363 "AG08BY.f"
    mp1 = *m + 1;
#line 364 "AG08BY.f"
    mui = 0;
#line 365 "AG08BY.f"
    dum[0] = 0.;

#line 367 "AG08BY.f"
    itau = 1;
#line 368 "AG08BY.f"
    jwork1 = itau + mpm;
#line 369 "AG08BY.f"
    ismin = (*p << 1) + 1;
#line 370 "AG08BY.f"
    ismax = ismin + *p;
#line 371 "AG08BY.f"
    jwork2 = ismax + *p;
#line 372 "AG08BY.f"
    nblcks = 0;
#line 373 "AG08BY.f"
    wrkopt = 1;

#line 375 "AG08BY.f"
L10:
#line 375 "AG08BY.f"
    if (*pr == 0) {
#line 375 "AG08BY.f"
	goto L90;
#line 375 "AG08BY.f"
    }

/*     (NR+1,ICOL+1) points to the current position of matrix D. */

#line 379 "AG08BY.f"
    ro1 = ro;
#line 380 "AG08BY.f"
    mnr = *m + *nr;
#line 381 "AG08BY.f"
    if (*m > 0) {

/*        Compress rows of D; first exploit the trapezoidal shape of the */
/*        (RO+SIGMA)-by-SIGMA matrix in the first SIGMA columns of D; */
/*        compress the first SIGMA columns without column pivoting: */

/*              ( x x x x x )       ( x x x x x ) */
/*              ( x x x x x )       ( 0 x x x x ) */
/*              ( x x x x x )  - >  ( 0 0 x x x ) */
/*              ( 0 x x x x )       ( 0 0 0 x x ) */
/*              ( 0 0 x x x )       ( 0 0 0 x x ) */

/*        where SIGMA = 3 and RO = 2. */
/*        Workspace: need maximum M+N-1. */

#line 396 "AG08BY.f"
	irow = *nr;
#line 397 "AG08BY.f"
	i__1 = sigma;
#line 397 "AG08BY.f"
	for (icol = 1; icol <= i__1; ++icol) {
#line 398 "AG08BY.f"
	    ++irow;
#line 399 "AG08BY.f"
	    i__2 = ro + 1;
#line 399 "AG08BY.f"
	    dlarfg_(&i__2, &abcd[irow + icol * abcd_dim1], &abcd[irow + 1 + 
		    icol * abcd_dim1], &c__1, &t);
#line 401 "AG08BY.f"
	    i__2 = ro + 1;
#line 401 "AG08BY.f"
	    i__3 = mnr - icol;
#line 401 "AG08BY.f"
	    dlatzm_("L", &i__2, &i__3, &abcd[irow + 1 + icol * abcd_dim1], &
		    c__1, &t, &abcd[irow + (icol + 1) * abcd_dim1], &abcd[
		    irow + 1 + (icol + 1) * abcd_dim1], ldabcd, &dwork[1], (
		    ftnlen)1);
#line 404 "AG08BY.f"
	    i__2 = *pr - icol;
#line 404 "AG08BY.f"
	    dcopy_(&i__2, dum, &c__0, &abcd[irow + 1 + icol * abcd_dim1], &
		    c__1);
#line 405 "AG08BY.f"
/* L20: */
#line 405 "AG08BY.f"
	}
/* Computing MAX */
#line 406 "AG08BY.f"
	i__1 = wrkopt, i__2 = mn - 1;
#line 406 "AG08BY.f"
	wrkopt = max(i__1,i__2);

#line 408 "AG08BY.f"
	if (*first) {

/*           Continue with Householder with column pivoting. */

/*              ( x x x x x )        ( x x x x x ) */
/*              ( 0 x x x x )        ( 0 x x x x ) */
/*              ( 0 0 x x x )  - >   ( 0 0 x x x ) */
/*              ( 0 0 0 x x )        ( 0 0 0 x x ) */
/*              ( 0 0 0 x x )        ( 0 0 0 0 0 ) */

/*           Real workspace:    need maximum min(P,M)+3*M-1; */
/*           Integer workspace: need maximum M. */

/* Computing MIN */
#line 421 "AG08BY.f"
	    i__1 = *nr + sigma + 1;
#line 421 "AG08BY.f"
	    irow = min(i__1,pn);
/* Computing MIN */
#line 422 "AG08BY.f"
	    i__1 = sigma + 1;
#line 422 "AG08BY.f"
	    icol = min(i__1,*m);
#line 423 "AG08BY.f"
	    i__1 = *m - sigma;
#line 423 "AG08BY.f"
	    mb03oy_(&ro1, &i__1, &abcd[irow + icol * abcd_dim1], ldabcd, &
		    rcond, svlmax, &rank, sval, &iwork[1], &dwork[itau], &
		    dwork[jwork1], info);
/* Computing MAX */
#line 426 "AG08BY.f"
	    i__1 = wrkopt, i__2 = jwork1 + *m * 3 - 2;
#line 426 "AG08BY.f"
	    wrkopt = max(i__1,i__2);

/*           Apply the column permutations to B and part of D. */

#line 430 "AG08BY.f"
	    i__1 = *nr + sigma;
#line 430 "AG08BY.f"
	    i__2 = *m - sigma;
#line 430 "AG08BY.f"
	    dlapmt_(&c_true, &i__1, &i__2, &abcd[icol * abcd_dim1 + 1], 
		    ldabcd, &iwork[1]);

#line 433 "AG08BY.f"
	    if (rank > 0) {

/*              Apply the Householder transformations to the submatrix C. */
/*              Workspace: need   maximum min(P,M) + N; */
/*                         prefer maximum min(P,M) + N*NB. */

#line 439 "AG08BY.f"
		i__1 = *ldwork - jwork1 + 1;
#line 439 "AG08BY.f"
		dormqr_("Left", "Transpose", &ro1, nr, &rank, &abcd[irow + 
			icol * abcd_dim1], ldabcd, &dwork[itau], &abcd[irow + 
			mp1 * abcd_dim1], ldabcd, &dwork[jwork1], &i__1, info,
			 (ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 443 "AG08BY.f"
		i__1 = wrkopt, i__2 = jwork1 + (integer) dwork[jwork1] - 1;
#line 443 "AG08BY.f"
		wrkopt = max(i__1,i__2);
#line 444 "AG08BY.f"
		i__1 = ro1 - 1;
/* Computing MIN */
#line 444 "AG08BY.f"
		i__3 = ro1 - 1;
#line 444 "AG08BY.f"
		i__2 = min(i__3,rank);
/* Computing MIN */
#line 444 "AG08BY.f"
		i__4 = irow + 1;
#line 444 "AG08BY.f"
		dlaset_("Lower", &i__1, &i__2, &c_b20, &c_b20, &abcd[min(i__4,
			pn) + icol * abcd_dim1], ldabcd, (ftnlen)5);
#line 446 "AG08BY.f"
		ro1 -= rank;
#line 447 "AG08BY.f"
	    }
#line 448 "AG08BY.f"
	}

/*        Terminate if Dr has maximal row rank. */

#line 452 "AG08BY.f"
	if (ro1 == 0) {
#line 452 "AG08BY.f"
	    goto L90;
#line 452 "AG08BY.f"
	}

#line 454 "AG08BY.f"
    }

/*     Update SIGMA. */

#line 458 "AG08BY.f"
    sigma = *pr - ro1;

#line 460 "AG08BY.f"
    ++nblcks;
#line 461 "AG08BY.f"
    taui = ro1;

/*     Compress the columns of current C to separate a TAUI-by-MUI */
/*     full column rank block. */

#line 466 "AG08BY.f"
    if (*nr == 0) {

/*        Finish for zero state dimension. */

#line 470 "AG08BY.f"
	*pr = sigma;
#line 471 "AG08BY.f"
	rank = 0;
#line 472 "AG08BY.f"
    } else {

/*        Perform RQ-decomposition with row pivoting on the current C */
/*        while keeping E upper triangular. */
/*        The current C is the TAUI-by-NR matrix delimited by rows */
/*        IRC+1 to IRC+TAUI and columns M+1 to M+NR of ABCD. */
/*        The rank of current C is computed in MUI. */
/*        Workspace: need maximum 5*P. */

#line 481 "AG08BY.f"
	irc = *nr + sigma;
#line 482 "AG08BY.f"
	n1 = *nr;
#line 483 "AG08BY.f"
	if (taui > 1) {

/*           Compute norms. */

#line 487 "AG08BY.f"
	    i__1 = taui;
#line 487 "AG08BY.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 488 "AG08BY.f"
		dwork[i__] = dnrm2_(nr, &abcd[irc + i__ + mp1 * abcd_dim1], 
			ldabcd);
#line 489 "AG08BY.f"
		dwork[*p + i__] = dwork[i__];
#line 490 "AG08BY.f"
/* L30: */
#line 490 "AG08BY.f"
	    }
#line 491 "AG08BY.f"
	}

#line 493 "AG08BY.f"
	rank = 0;
#line 494 "AG08BY.f"
	mntau = min(taui,*nr);

/*        ICOL and IROW will point to the current pivot position in C. */

#line 498 "AG08BY.f"
	ilast = *nr + *pr;
#line 499 "AG08BY.f"
	jlast = *m + *nr;
#line 500 "AG08BY.f"
	irow = ilast;
#line 501 "AG08BY.f"
	icol = jlast;
#line 502 "AG08BY.f"
	i__ = taui;
#line 503 "AG08BY.f"
L40:
#line 503 "AG08BY.f"
	if (rank < mntau) {
#line 504 "AG08BY.f"
	    mn1 = *m + n1;

/*           Pivot if necessary. */

#line 508 "AG08BY.f"
	    if (i__ != 1) {
#line 509 "AG08BY.f"
		j = idamax_(&i__, &dwork[1], &c__1);
#line 510 "AG08BY.f"
		if (j != i__) {
#line 511 "AG08BY.f"
		    dwork[j] = dwork[i__];
#line 512 "AG08BY.f"
		    dwork[*p + j] = dwork[*p + i__];
#line 513 "AG08BY.f"
		    dswap_(&n1, &abcd[irow + mp1 * abcd_dim1], ldabcd, &abcd[
			    irc + j + mp1 * abcd_dim1], ldabcd);
#line 515 "AG08BY.f"
		}
#line 516 "AG08BY.f"
	    }

/*           Zero elements left to ABCD(IROW,ICOL). */

#line 520 "AG08BY.f"
	    i__1 = n1 - 1;
#line 520 "AG08BY.f"
	    for (k = 1; k <= i__1; ++k) {
#line 521 "AG08BY.f"
		j = *m + k;

/*              Rotate columns J, J+1 to zero ABCD(IROW,J). */

#line 525 "AG08BY.f"
		t = abcd[irow + (j + 1) * abcd_dim1];
#line 526 "AG08BY.f"
		dlartg_(&t, &abcd[irow + j * abcd_dim1], &c__, &s, &abcd[irow 
			+ (j + 1) * abcd_dim1]);
#line 527 "AG08BY.f"
		abcd[irow + j * abcd_dim1] = 0.;
#line 528 "AG08BY.f"
		i__2 = irow - 1;
#line 528 "AG08BY.f"
		drot_(&i__2, &abcd[(j + 1) * abcd_dim1 + 1], &c__1, &abcd[j * 
			abcd_dim1 + 1], &c__1, &c__, &s);
#line 529 "AG08BY.f"
		i__2 = k + 1;
#line 529 "AG08BY.f"
		drot_(&i__2, &e[(k + 1) * e_dim1 + 1], &c__1, &e[k * e_dim1 + 
			1], &c__1, &c__, &s);

/*              Rotate rows K, K+1 to zero E(K+1,K). */

#line 533 "AG08BY.f"
		t = e[k + k * e_dim1];
#line 534 "AG08BY.f"
		dlartg_(&t, &e[k + 1 + k * e_dim1], &c__, &s, &e[k + k * 
			e_dim1]);
#line 535 "AG08BY.f"
		e[k + 1 + k * e_dim1] = 0.;
#line 536 "AG08BY.f"
		i__2 = n1 - k;
#line 536 "AG08BY.f"
		drot_(&i__2, &e[k + (k + 1) * e_dim1], lde, &e[k + 1 + (k + 1)
			 * e_dim1], lde, &c__, &s);
#line 537 "AG08BY.f"
		drot_(&mn1, &abcd[k + abcd_dim1], ldabcd, &abcd[k + 1 + 
			abcd_dim1], ldabcd, &c__, &s);
#line 539 "AG08BY.f"
/* L50: */
#line 539 "AG08BY.f"
	    }

#line 541 "AG08BY.f"
	    if (rank == 0) {

/*              Initialize; exit if matrix is zero (RANK = 0). */

#line 545 "AG08BY.f"
		smax = (d__1 = abcd[ilast + jlast * abcd_dim1], abs(d__1));
#line 546 "AG08BY.f"
		if (smax == 0.) {
#line 546 "AG08BY.f"
		    goto L80;
#line 546 "AG08BY.f"
		}
#line 547 "AG08BY.f"
		smin = smax;
#line 548 "AG08BY.f"
		smaxpr = smax;
#line 549 "AG08BY.f"
		sminpr = smin;
#line 550 "AG08BY.f"
		c1 = 1.;
#line 551 "AG08BY.f"
		c2 = 1.;
#line 552 "AG08BY.f"
	    } else {

/*              One step of incremental condition estimation. */

#line 556 "AG08BY.f"
		dcopy_(&rank, &abcd[irow + (icol + 1) * abcd_dim1], ldabcd, &
			dwork[jwork2], &c__1);
#line 558 "AG08BY.f"
		dlaic1_(&c__2, &rank, &dwork[ismin], &smin, &dwork[jwork2], &
			abcd[irow + icol * abcd_dim1], &sminpr, &s1, &c1);
#line 561 "AG08BY.f"
		dlaic1_(&c__1, &rank, &dwork[ismax], &smax, &dwork[jwork2], &
			abcd[irow + icol * abcd_dim1], &smaxpr, &s2, &c2);
/* Computing MAX */
#line 564 "AG08BY.f"
		i__1 = wrkopt, i__2 = *p * 5;
#line 564 "AG08BY.f"
		wrkopt = max(i__1,i__2);
#line 565 "AG08BY.f"
	    }

/*           Check the rank; finish the loop if rank loss occurs. */

#line 569 "AG08BY.f"
	    if (*svlmax * rcond <= smaxpr) {
#line 570 "AG08BY.f"
		if (*svlmax * rcond <= sminpr) {
#line 571 "AG08BY.f"
		    if (smaxpr * rcond <= sminpr) {

/*                    Finish the loop if last row. */

#line 575 "AG08BY.f"
			if (n1 == 0) {
#line 576 "AG08BY.f"
			    ++rank;
#line 577 "AG08BY.f"
			    goto L80;
#line 578 "AG08BY.f"
			}

#line 580 "AG08BY.f"
			if (n1 > 1) {

/*                       Update norms. */

#line 584 "AG08BY.f"
			    if (i__ - 1 > 1) {
#line 585 "AG08BY.f"
				i__1 = i__ - 1;
#line 585 "AG08BY.f"
				for (j = 1; j <= i__1; ++j) {
#line 586 "AG08BY.f"
				    if (dwork[j] != 0.) {
/* Computing 2nd power */
#line 587 "AG08BY.f"
					d__2 = (d__1 = abcd[irc + j + icol * 
						abcd_dim1], abs(d__1)) / 
						dwork[j];
#line 587 "AG08BY.f"
					t = 1. - d__2 * d__2;
#line 589 "AG08BY.f"
					t = max(t,0.);
/* Computing 2nd power */
#line 590 "AG08BY.f"
					d__1 = dwork[j] / dwork[*p + j];
#line 590 "AG08BY.f"
					tt = t * .05 * (d__1 * d__1) + 1.;
#line 592 "AG08BY.f"
					if (tt != 1.) {
#line 593 "AG08BY.f"
					    dwork[j] *= sqrt(t);
#line 594 "AG08BY.f"
					} else {
#line 595 "AG08BY.f"
					    i__2 = n1 - 1;
#line 595 "AG08BY.f"
					    dwork[j] = dnrm2_(&i__2, &abcd[
						    irc + j + mp1 * abcd_dim1]
						    , ldabcd);
#line 597 "AG08BY.f"
					    dwork[*p + j] = dwork[j];
#line 598 "AG08BY.f"
					}
#line 599 "AG08BY.f"
				    }
#line 600 "AG08BY.f"
/* L60: */
#line 600 "AG08BY.f"
				}
#line 601 "AG08BY.f"
			    }
#line 602 "AG08BY.f"
			}

#line 604 "AG08BY.f"
			i__1 = rank;
#line 604 "AG08BY.f"
			for (j = 1; j <= i__1; ++j) {
#line 605 "AG08BY.f"
			    dwork[ismin + j - 1] = s1 * dwork[ismin + j - 1];
#line 606 "AG08BY.f"
			    dwork[ismax + j - 1] = s2 * dwork[ismax + j - 1];
#line 607 "AG08BY.f"
/* L70: */
#line 607 "AG08BY.f"
			}

#line 609 "AG08BY.f"
			dwork[ismin + rank] = c1;
#line 610 "AG08BY.f"
			dwork[ismax + rank] = c2;
#line 611 "AG08BY.f"
			smin = sminpr;
#line 612 "AG08BY.f"
			smax = smaxpr;
#line 613 "AG08BY.f"
			++rank;
#line 614 "AG08BY.f"
			--icol;
#line 615 "AG08BY.f"
			--irow;
#line 616 "AG08BY.f"
			--n1;
#line 617 "AG08BY.f"
			--i__;
#line 618 "AG08BY.f"
			goto L40;
#line 619 "AG08BY.f"
		    }
#line 620 "AG08BY.f"
		}
#line 621 "AG08BY.f"
	    }
#line 622 "AG08BY.f"
	}
#line 623 "AG08BY.f"
    }

#line 625 "AG08BY.f"
L80:
#line 626 "AG08BY.f"
    mui = rank;
#line 627 "AG08BY.f"
    *nr -= mui;
#line 628 "AG08BY.f"
    *pr = sigma + mui;

/*     Set number of left Kronecker blocks of order (i-1)-by-i. */

#line 632 "AG08BY.f"
    kronl[nblcks] = taui - mui;

/*     Set number of infinite divisors of order i-1. */

#line 636 "AG08BY.f"
    if (*first && nblcks > 1) {
#line 636 "AG08BY.f"
	infz[nblcks - 1] = muim1 - taui;
#line 636 "AG08BY.f"
    }
#line 638 "AG08BY.f"
    muim1 = mui;
#line 639 "AG08BY.f"
    ro = mui;

/*     Continue reduction if rank of current C is positive. */

#line 643 "AG08BY.f"
    if (mui > 0) {
#line 643 "AG08BY.f"
	goto L10;
#line 643 "AG08BY.f"
    }

/*     Determine the maximal degree of infinite zeros and */
/*     the number of infinite zeros. */

#line 649 "AG08BY.f"
L90:
#line 650 "AG08BY.f"
    if (*first) {
#line 651 "AG08BY.f"
	if (mui == 0) {
/* Computing MAX */
#line 652 "AG08BY.f"
	    i__1 = 0, i__2 = nblcks - 1;
#line 652 "AG08BY.f"
	    *dinfz = max(i__1,i__2);
#line 653 "AG08BY.f"
	} else {
#line 654 "AG08BY.f"
	    *dinfz = nblcks;
#line 655 "AG08BY.f"
	    infz[nblcks] = mui;
#line 656 "AG08BY.f"
	}
#line 657 "AG08BY.f"
	k = *dinfz;
#line 658 "AG08BY.f"
	for (i__ = k; i__ >= 1; --i__) {
#line 659 "AG08BY.f"
	    if (infz[i__] != 0) {
#line 659 "AG08BY.f"
		goto L110;
#line 659 "AG08BY.f"
	    }
#line 660 "AG08BY.f"
	    --(*dinfz);
#line 661 "AG08BY.f"
/* L100: */
#line 661 "AG08BY.f"
	}
#line 662 "AG08BY.f"
L110:
#line 663 "AG08BY.f"
	i__1 = *dinfz;
#line 663 "AG08BY.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 664 "AG08BY.f"
	    *ninfz += infz[i__] * i__;
#line 665 "AG08BY.f"
/* L120: */
#line 665 "AG08BY.f"
	}
#line 666 "AG08BY.f"
    }

/*     Determine the maximal order of left elementary Kronecker blocks. */

#line 670 "AG08BY.f"
    *nkronl = nblcks;
#line 671 "AG08BY.f"
    for (i__ = nblcks; i__ >= 1; --i__) {
#line 672 "AG08BY.f"
	if (kronl[i__] != 0) {
#line 672 "AG08BY.f"
	    goto L140;
#line 672 "AG08BY.f"
	}
#line 673 "AG08BY.f"
	--(*nkronl);
#line 674 "AG08BY.f"
/* L130: */
#line 674 "AG08BY.f"
    }
#line 675 "AG08BY.f"
L140:

#line 677 "AG08BY.f"
    dwork[1] = (doublereal) wrkopt;
#line 678 "AG08BY.f"
    return 0;
/* *** Last line of AG08BY *** */
} /* ag08by_ */

