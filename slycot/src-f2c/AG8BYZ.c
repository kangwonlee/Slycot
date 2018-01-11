#line 1 "AG8BYZ.f"
/* AG8BYZ.f -- translated by f2c (version 20100827).
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

#line 1 "AG8BYZ.f"
/* Table of constant values */

static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static logical c_true = TRUE_;
static integer c__2 = 2;

/* Subroutine */ int ag8byz_(logical *first, integer *n, integer *m, integer *
	p, doublereal *svlmax, doublecomplex *abcd, integer *ldabcd, 
	doublecomplex *e, integer *lde, integer *nr, integer *pr, integer *
	ninfz, integer *dinfz, integer *nkronl, integer *infz, integer *kronl,
	 doublereal *tol, integer *iwork, doublereal *dwork, doublecomplex *
	zwork, integer *lzwork, integer *info)
{
    /* System generated locals */
    integer abcd_dim1, abcd_offset, e_dim1, e_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static doublereal c__;
    static integer i__, j, k;
    static doublecomplex s;
    static doublereal t;
    static doublecomplex c1, c2;
    static integer n1;
    static doublecomplex s1, s2;
    static integer nb;
    static doublecomplex tc;
    static integer mn, pn, ro;
    static doublereal tt;
    static integer mn1, mp1, ro1, irc;
    static doublecomplex dum[1];
    static integer mpm, mui, mnr, icol, rank, itau, taui;
    static doublereal sval[3], smin, smax;
    static integer irow;
    extern /* Subroutine */ int zrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static integer muim1, sigma;
    static doublereal rcond;
    static integer ilast, jlast, ismin, ismax, mntau;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), zlaic1_(integer *, 
	    integer *, doublecomplex *, doublereal *, doublecomplex *, 
	    doublecomplex *, doublereal *, doublecomplex *, doublecomplex *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);
    static integer jwork1, jwork2;
    extern /* Subroutine */ int mb3oyz_(integer *, integer *, doublecomplex *,
	     integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublecomplex *, doublereal *, doublecomplex *, 
	    integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    static integer nblcks;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zlarfg_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *), zlaset_(char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen), zlartg_(doublecomplex *, 
	    doublecomplex *, doublereal *, doublecomplex *, doublecomplex *), 
	    zlapmt_(logical *, integer *, integer *, doublecomplex *, integer 
	    *, integer *);
    static doublereal sminpr, smaxpr;
    static logical lquery;
    extern /* Subroutine */ int zlatzm_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, ftnlen);
    static integer wrkopt;
    extern /* Subroutine */ int zunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);


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
/*             Specifies if AG8BYZ is called first time or it is called */
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

/*     ABCD    (input/output) COMPLEX*16 array, dimension (LDABCD,M+N) */
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

/*     E       (input/output) COMPLEX*16 array, dimension (LDE,N) */
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
/*             LDWORK >= 2*MAX(M,P), if FIRST = .TRUE.; */
/*             LDWORK >= 2*P,        if FIRST = .FALSE. . */

/*     ZWORK   COMPLEX*16 array, dimension (LZWORK) */
/*             On exit, if INFO = 0, ZWORK(1) returns the optimal value */
/*             of LZWORK. */

/*     LZWORK  INTEGER */
/*             The length of the array ZWORK. */
/*             LZWORK >= 1, if P = 0; otherwise */
/*             LZWORK >= MAX( 1, N+M-1, MIN(P,M) + MAX(3*M-1,N), 3*P ), */
/*                                             if FIRST = .TRUE.; */
/*             LZWORK >= MAX( 1, N+M-1, 3*P ), if FIRST = .FALSE. . */
/*             The second term is not needed if M = 0. */
/*             For optimum performance LZWORK should be larger. */

/*             If LZWORK = -1, then a workspace query is assumed; */
/*             the routine only calculates the optimal size of the */
/*             ZWORK array, returns this value as the first entry of */
/*             the ZWORK array, and no error message related to LZWORK */
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
/*     May 1999. */
/*     Complex version: V. Sima, Research Institute for Informatics, */
/*     Bucharest, Nov. 2008. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2009. */

/*     KEYWORDS */

/*     Generalized eigenvalue problem, Kronecker indices, multivariable */
/*     system, unitary transformation, structural invariant. */

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

#line 284 "AG8BYZ.f"
    /* Parameter adjustments */
#line 284 "AG8BYZ.f"
    abcd_dim1 = *ldabcd;
#line 284 "AG8BYZ.f"
    abcd_offset = 1 + abcd_dim1;
#line 284 "AG8BYZ.f"
    abcd -= abcd_offset;
#line 284 "AG8BYZ.f"
    e_dim1 = *lde;
#line 284 "AG8BYZ.f"
    e_offset = 1 + e_dim1;
#line 284 "AG8BYZ.f"
    e -= e_offset;
#line 284 "AG8BYZ.f"
    --infz;
#line 284 "AG8BYZ.f"
    --kronl;
#line 284 "AG8BYZ.f"
    --iwork;
#line 284 "AG8BYZ.f"
    --dwork;
#line 284 "AG8BYZ.f"
    --zwork;
#line 284 "AG8BYZ.f"

#line 284 "AG8BYZ.f"
    /* Function Body */
#line 284 "AG8BYZ.f"
    lquery = *lzwork == -1;
#line 285 "AG8BYZ.f"
    *info = 0;
#line 286 "AG8BYZ.f"
    pn = *p + *n;
#line 287 "AG8BYZ.f"
    mn = *m + *n;
#line 288 "AG8BYZ.f"
    mpm = min(*p,*m);
#line 289 "AG8BYZ.f"
    if (*n < 0) {
#line 290 "AG8BYZ.f"
	*info = -2;
#line 291 "AG8BYZ.f"
    } else if (*m < 0 || ! (*first) && *m > *p) {
#line 292 "AG8BYZ.f"
	*info = -3;
#line 293 "AG8BYZ.f"
    } else if (*p < 0) {
#line 294 "AG8BYZ.f"
	*info = -4;
#line 295 "AG8BYZ.f"
    } else if (*svlmax < 0.) {
#line 296 "AG8BYZ.f"
	*info = -5;
#line 297 "AG8BYZ.f"
    } else if (*ldabcd < max(1,pn)) {
#line 298 "AG8BYZ.f"
	*info = -7;
#line 299 "AG8BYZ.f"
    } else if (*lde < max(1,*n)) {
#line 300 "AG8BYZ.f"
	*info = -9;
#line 301 "AG8BYZ.f"
    } else if (*tol > 1.) {
#line 302 "AG8BYZ.f"
	*info = -17;
#line 303 "AG8BYZ.f"
    } else {
/* Computing MAX */
#line 304 "AG8BYZ.f"
	i__1 = 1, i__2 = *p * 3;
#line 304 "AG8BYZ.f"
	wrkopt = max(i__1,i__2);
#line 305 "AG8BYZ.f"
	if (*p > 0) {
#line 306 "AG8BYZ.f"
	    if (*m > 0) {
/* Computing MAX */
#line 307 "AG8BYZ.f"
		i__1 = wrkopt, i__2 = mn - 1;
#line 307 "AG8BYZ.f"
		wrkopt = max(i__1,i__2);
#line 308 "AG8BYZ.f"
		if (*first) {
/* Computing MAX */
/* Computing MAX */
#line 309 "AG8BYZ.f"
		    i__3 = *m * 3 - 1;
#line 309 "AG8BYZ.f"
		    i__1 = wrkopt, i__2 = mpm + max(i__3,*n);
#line 309 "AG8BYZ.f"
		    wrkopt = max(i__1,i__2);
#line 310 "AG8BYZ.f"
		    if (lquery) {
/* Computing MIN */
#line 311 "AG8BYZ.f"
			i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMQR", "LC", p, n,
				 &mpm, &c_n1, (ftnlen)6, (ftnlen)2);
#line 311 "AG8BYZ.f"
			nb = min(i__1,i__2);
/* Computing MAX */
#line 313 "AG8BYZ.f"
			i__1 = wrkopt, i__2 = mpm + max(1,*n) * nb;
#line 313 "AG8BYZ.f"
			wrkopt = max(i__1,i__2);
#line 314 "AG8BYZ.f"
		    }
#line 315 "AG8BYZ.f"
		}
#line 316 "AG8BYZ.f"
	    }
#line 317 "AG8BYZ.f"
	}
#line 318 "AG8BYZ.f"
	if (*lzwork < wrkopt && ! lquery) {
#line 319 "AG8BYZ.f"
	    *info = -21;
#line 320 "AG8BYZ.f"
	}
#line 321 "AG8BYZ.f"
    }
#line 322 "AG8BYZ.f"
    if (*info != 0) {
#line 323 "AG8BYZ.f"
	i__1 = -(*info);
#line 323 "AG8BYZ.f"
	xerbla_("AG8BYZ", &i__1, (ftnlen)6);
#line 324 "AG8BYZ.f"
	return 0;
#line 325 "AG8BYZ.f"
    } else if (lquery) {
#line 326 "AG8BYZ.f"
	zwork[1].r = (doublereal) wrkopt, zwork[1].i = 0.;
#line 327 "AG8BYZ.f"
	return 0;
#line 328 "AG8BYZ.f"
    }

/*     Initialize output variables. */

#line 332 "AG8BYZ.f"
    *pr = *p;
#line 333 "AG8BYZ.f"
    *nr = *n;
#line 334 "AG8BYZ.f"
    *dinfz = 0;
#line 335 "AG8BYZ.f"
    *ninfz = 0;
#line 336 "AG8BYZ.f"
    *nkronl = 0;

/*     Quick return if possible. */

#line 340 "AG8BYZ.f"
    if (*p == 0) {
#line 341 "AG8BYZ.f"
	zwork[1].r = 1., zwork[1].i = 0.;
#line 342 "AG8BYZ.f"
	return 0;
#line 343 "AG8BYZ.f"
    }
#line 344 "AG8BYZ.f"
    if (*n == 0 && *m == 0) {
#line 345 "AG8BYZ.f"
	*pr = 0;
#line 346 "AG8BYZ.f"
	*nkronl = 1;
#line 347 "AG8BYZ.f"
	kronl[1] = *p;
#line 348 "AG8BYZ.f"
	zwork[1].r = 1., zwork[1].i = 0.;
#line 349 "AG8BYZ.f"
	return 0;
#line 350 "AG8BYZ.f"
    }

#line 352 "AG8BYZ.f"
    rcond = *tol;
#line 353 "AG8BYZ.f"
    if (rcond <= 0.) {

/*        Use the default tolerance in rank determination. */

#line 357 "AG8BYZ.f"
	rcond = (doublereal) (pn * mn) * dlamch_("EPSILON", (ftnlen)7);
#line 358 "AG8BYZ.f"
    }

/*     The D matrix is (RO+SIGMA)-by-M, where RO = P - SIGMA and */
/*     SIGMA = 0 for FIRST = .TRUE. and SIGMA = M for FIRST = .FALSE.. */
/*     The leading (RO+SIGMA)-by-SIGMA submatrix of D has full column */
/*     rank, with the trailing SIGMA-by-SIGMA submatrix upper triangular. */

#line 365 "AG8BYZ.f"
    if (*first) {
#line 366 "AG8BYZ.f"
	sigma = 0;
#line 367 "AG8BYZ.f"
    } else {
#line 368 "AG8BYZ.f"
	sigma = *m;
#line 369 "AG8BYZ.f"
    }
#line 370 "AG8BYZ.f"
    ro = *p - sigma;
#line 371 "AG8BYZ.f"
    mp1 = *m + 1;
#line 372 "AG8BYZ.f"
    mui = 0;
#line 373 "AG8BYZ.f"
    dum[0].r = 0., dum[0].i = 0.;

#line 375 "AG8BYZ.f"
    itau = 1;
#line 376 "AG8BYZ.f"
    jwork1 = itau + mpm;
#line 377 "AG8BYZ.f"
    ismin = 1;
#line 378 "AG8BYZ.f"
    ismax = ismin + *p;
#line 379 "AG8BYZ.f"
    jwork2 = ismax + *p;
#line 380 "AG8BYZ.f"
    nblcks = 0;
#line 381 "AG8BYZ.f"
    wrkopt = 1;

#line 383 "AG8BYZ.f"
L10:
#line 383 "AG8BYZ.f"
    if (*pr == 0) {
#line 383 "AG8BYZ.f"
	goto L90;
#line 383 "AG8BYZ.f"
    }

/*     (NR+1,ICOL+1) points to the current position of matrix D. */

#line 387 "AG8BYZ.f"
    ro1 = ro;
#line 388 "AG8BYZ.f"
    mnr = *m + *nr;
#line 389 "AG8BYZ.f"
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
/*        Complex workspace: need maximum M+N-1. */

#line 404 "AG8BYZ.f"
	irow = *nr;
#line 405 "AG8BYZ.f"
	i__1 = sigma;
#line 405 "AG8BYZ.f"
	for (icol = 1; icol <= i__1; ++icol) {
#line 406 "AG8BYZ.f"
	    ++irow;
#line 407 "AG8BYZ.f"
	    i__2 = ro + 1;
#line 407 "AG8BYZ.f"
	    zlarfg_(&i__2, &abcd[irow + icol * abcd_dim1], &abcd[irow + 1 + 
		    icol * abcd_dim1], &c__1, &tc);
#line 409 "AG8BYZ.f"
	    i__2 = ro + 1;
#line 409 "AG8BYZ.f"
	    i__3 = mnr - icol;
#line 409 "AG8BYZ.f"
	    d_cnjg(&z__1, &tc);
#line 409 "AG8BYZ.f"
	    zlatzm_("L", &i__2, &i__3, &abcd[irow + 1 + icol * abcd_dim1], &
		    c__1, &z__1, &abcd[irow + (icol + 1) * abcd_dim1], &abcd[
		    irow + 1 + (icol + 1) * abcd_dim1], ldabcd, &zwork[1], (
		    ftnlen)1);
#line 412 "AG8BYZ.f"
	    i__2 = *pr - icol;
#line 412 "AG8BYZ.f"
	    zcopy_(&i__2, dum, &c__0, &abcd[irow + 1 + icol * abcd_dim1], &
		    c__1);
#line 413 "AG8BYZ.f"
/* L20: */
#line 413 "AG8BYZ.f"
	}
/* Computing MAX */
#line 414 "AG8BYZ.f"
	i__1 = wrkopt, i__2 = mn - 1;
#line 414 "AG8BYZ.f"
	wrkopt = max(i__1,i__2);

#line 416 "AG8BYZ.f"
	if (*first) {

/*           Continue with Householder with column pivoting. */

/*              ( x x x x x )        ( x x x x x ) */
/*              ( 0 x x x x )        ( 0 x x x x ) */
/*              ( 0 0 x x x )  - >   ( 0 0 x x x ) */
/*              ( 0 0 0 x x )        ( 0 0 0 x x ) */
/*              ( 0 0 0 x x )        ( 0 0 0 0 0 ) */

/*           Real workspace:    need maximum 2*M; */
/*           Complex workspace: need maximum min(P,M)+3*M-1; */
/*           Integer workspace: need maximum M. */

/* Computing MIN */
#line 430 "AG8BYZ.f"
	    i__1 = *nr + sigma + 1;
#line 430 "AG8BYZ.f"
	    irow = min(i__1,pn);
/* Computing MIN */
#line 431 "AG8BYZ.f"
	    i__1 = sigma + 1;
#line 431 "AG8BYZ.f"
	    icol = min(i__1,*m);
#line 432 "AG8BYZ.f"
	    i__1 = *m - sigma;
#line 432 "AG8BYZ.f"
	    mb3oyz_(&ro1, &i__1, &abcd[irow + icol * abcd_dim1], ldabcd, &
		    rcond, svlmax, &rank, sval, &iwork[1], &zwork[itau], &
		    dwork[1], &zwork[jwork1], info);
/* Computing MAX */
#line 435 "AG8BYZ.f"
	    i__1 = wrkopt, i__2 = jwork1 + *m * 3 - 2;
#line 435 "AG8BYZ.f"
	    wrkopt = max(i__1,i__2);

/*           Apply the column permutations to B and part of D. */

#line 439 "AG8BYZ.f"
	    i__1 = *nr + sigma;
#line 439 "AG8BYZ.f"
	    i__2 = *m - sigma;
#line 439 "AG8BYZ.f"
	    zlapmt_(&c_true, &i__1, &i__2, &abcd[icol * abcd_dim1 + 1], 
		    ldabcd, &iwork[1]);

#line 442 "AG8BYZ.f"
	    if (rank > 0) {

/*              Apply the Householder transformations to the submatrix C. */
/*              Complex workspace: need   maximum min(P,M) + N; */
/*                                 prefer maximum min(P,M) + N*NB. */

#line 448 "AG8BYZ.f"
		i__1 = *lzwork - jwork1 + 1;
#line 448 "AG8BYZ.f"
		zunmqr_("Left", "ConjTranspose", &ro1, nr, &rank, &abcd[irow 
			+ icol * abcd_dim1], ldabcd, &zwork[itau], &abcd[irow 
			+ mp1 * abcd_dim1], ldabcd, &zwork[jwork1], &i__1, 
			info, (ftnlen)4, (ftnlen)13);
/* Computing MAX */
#line 452 "AG8BYZ.f"
		i__3 = jwork1;
#line 452 "AG8BYZ.f"
		i__1 = wrkopt, i__2 = jwork1 + (integer) zwork[i__3].r - 1;
#line 452 "AG8BYZ.f"
		wrkopt = max(i__1,i__2);
#line 453 "AG8BYZ.f"
		i__1 = ro1 - 1;
/* Computing MIN */
#line 453 "AG8BYZ.f"
		i__3 = ro1 - 1;
#line 453 "AG8BYZ.f"
		i__2 = min(i__3,rank);
/* Computing MIN */
#line 453 "AG8BYZ.f"
		i__4 = irow + 1;
#line 453 "AG8BYZ.f"
		zlaset_("Lower", &i__1, &i__2, &c_b2, &c_b2, &abcd[min(i__4,
			pn) + icol * abcd_dim1], ldabcd, (ftnlen)5);
#line 456 "AG8BYZ.f"
		ro1 -= rank;
#line 457 "AG8BYZ.f"
	    }
#line 458 "AG8BYZ.f"
	}

/*        Terminate if Dr has maximal row rank. */

#line 462 "AG8BYZ.f"
	if (ro1 == 0) {
#line 462 "AG8BYZ.f"
	    goto L90;
#line 462 "AG8BYZ.f"
	}

#line 464 "AG8BYZ.f"
    }

/*     Update SIGMA. */

#line 468 "AG8BYZ.f"
    sigma = *pr - ro1;

#line 470 "AG8BYZ.f"
    ++nblcks;
#line 471 "AG8BYZ.f"
    taui = ro1;

/*     Compress the columns of current C to separate a TAUI-by-MUI */
/*     full column rank block. */

#line 476 "AG8BYZ.f"
    if (*nr == 0) {

/*        Finish for zero state dimension. */

#line 480 "AG8BYZ.f"
	*pr = sigma;
#line 481 "AG8BYZ.f"
	rank = 0;
#line 482 "AG8BYZ.f"
    } else {

/*        Perform RQ-decomposition with row pivoting on the current C */
/*        while keeping E upper triangular. */
/*        The current C is the TAUI-by-NR matrix delimited by rows */
/*        IRC+1 to IRC+TAUI and columns M+1 to M+NR of ABCD. */
/*        The rank of current C is computed in MUI. */
/*        Real workspace:    need maximum 2*P; */
/*        Complex workspace: need maximum 3*P. */

#line 492 "AG8BYZ.f"
	irc = *nr + sigma;
#line 493 "AG8BYZ.f"
	n1 = *nr;
#line 494 "AG8BYZ.f"
	if (taui > 1) {

/*           Compute norms. */

#line 498 "AG8BYZ.f"
	    i__1 = taui;
#line 498 "AG8BYZ.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 499 "AG8BYZ.f"
		dwork[i__] = dznrm2_(nr, &abcd[irc + i__ + mp1 * abcd_dim1], 
			ldabcd);
#line 500 "AG8BYZ.f"
		dwork[*p + i__] = dwork[i__];
#line 501 "AG8BYZ.f"
/* L30: */
#line 501 "AG8BYZ.f"
	    }
#line 502 "AG8BYZ.f"
	}

#line 504 "AG8BYZ.f"
	rank = 0;
#line 505 "AG8BYZ.f"
	mntau = min(taui,*nr);

/*        ICOL and IROW will point to the current pivot position in C. */

#line 509 "AG8BYZ.f"
	ilast = *nr + *pr;
#line 510 "AG8BYZ.f"
	jlast = *m + *nr;
#line 511 "AG8BYZ.f"
	irow = ilast;
#line 512 "AG8BYZ.f"
	icol = jlast;
#line 513 "AG8BYZ.f"
	i__ = taui;
#line 514 "AG8BYZ.f"
L40:
#line 514 "AG8BYZ.f"
	if (rank < mntau) {
#line 515 "AG8BYZ.f"
	    mn1 = *m + n1;

/*           Pivot if necessary. */

#line 519 "AG8BYZ.f"
	    if (i__ != 1) {
#line 520 "AG8BYZ.f"
		j = idamax_(&i__, &dwork[1], &c__1);
#line 521 "AG8BYZ.f"
		if (j != i__) {
#line 522 "AG8BYZ.f"
		    dwork[j] = dwork[i__];
#line 523 "AG8BYZ.f"
		    dwork[*p + j] = dwork[*p + i__];
#line 524 "AG8BYZ.f"
		    zswap_(&n1, &abcd[irow + mp1 * abcd_dim1], ldabcd, &abcd[
			    irc + j + mp1 * abcd_dim1], ldabcd);
#line 526 "AG8BYZ.f"
		}
#line 527 "AG8BYZ.f"
	    }

/*           Zero elements left to ABCD(IROW,ICOL). */

#line 531 "AG8BYZ.f"
	    i__1 = n1 - 1;
#line 531 "AG8BYZ.f"
	    for (k = 1; k <= i__1; ++k) {
#line 532 "AG8BYZ.f"
		j = *m + k;

/*              Rotate columns J, J+1 to zero ABCD(IROW,J). */

#line 536 "AG8BYZ.f"
		i__2 = irow + (j + 1) * abcd_dim1;
#line 536 "AG8BYZ.f"
		tc.r = abcd[i__2].r, tc.i = abcd[i__2].i;
#line 537 "AG8BYZ.f"
		zlartg_(&tc, &abcd[irow + j * abcd_dim1], &c__, &s, &abcd[
			irow + (j + 1) * abcd_dim1]);
#line 538 "AG8BYZ.f"
		i__2 = irow + j * abcd_dim1;
#line 538 "AG8BYZ.f"
		abcd[i__2].r = 0., abcd[i__2].i = 0.;
#line 539 "AG8BYZ.f"
		i__2 = irow - 1;
#line 539 "AG8BYZ.f"
		zrot_(&i__2, &abcd[(j + 1) * abcd_dim1 + 1], &c__1, &abcd[j * 
			abcd_dim1 + 1], &c__1, &c__, &s);
#line 540 "AG8BYZ.f"
		i__2 = k + 1;
#line 540 "AG8BYZ.f"
		zrot_(&i__2, &e[(k + 1) * e_dim1 + 1], &c__1, &e[k * e_dim1 + 
			1], &c__1, &c__, &s);

/*              Rotate rows K, K+1 to zero E(K+1,K). */

#line 544 "AG8BYZ.f"
		i__2 = k + k * e_dim1;
#line 544 "AG8BYZ.f"
		tc.r = e[i__2].r, tc.i = e[i__2].i;
#line 545 "AG8BYZ.f"
		zlartg_(&tc, &e[k + 1 + k * e_dim1], &c__, &s, &e[k + k * 
			e_dim1]);
#line 546 "AG8BYZ.f"
		i__2 = k + 1 + k * e_dim1;
#line 546 "AG8BYZ.f"
		e[i__2].r = 0., e[i__2].i = 0.;
#line 547 "AG8BYZ.f"
		i__2 = n1 - k;
#line 547 "AG8BYZ.f"
		zrot_(&i__2, &e[k + (k + 1) * e_dim1], lde, &e[k + 1 + (k + 1)
			 * e_dim1], lde, &c__, &s);
#line 548 "AG8BYZ.f"
		zrot_(&mn1, &abcd[k + abcd_dim1], ldabcd, &abcd[k + 1 + 
			abcd_dim1], ldabcd, &c__, &s);
#line 550 "AG8BYZ.f"
/* L50: */
#line 550 "AG8BYZ.f"
	    }

#line 552 "AG8BYZ.f"
	    if (rank == 0) {

/*              Initialize; exit if matrix is zero (RANK = 0). */

#line 556 "AG8BYZ.f"
		smax = z_abs(&abcd[ilast + jlast * abcd_dim1]);
#line 557 "AG8BYZ.f"
		if (smax == 0.) {
#line 557 "AG8BYZ.f"
		    goto L80;
#line 557 "AG8BYZ.f"
		}
#line 558 "AG8BYZ.f"
		smin = smax;
#line 559 "AG8BYZ.f"
		smaxpr = smax;
#line 560 "AG8BYZ.f"
		sminpr = smin;
#line 561 "AG8BYZ.f"
		c1.r = 1., c1.i = 0.;
#line 562 "AG8BYZ.f"
		c2.r = 1., c2.i = 0.;
#line 563 "AG8BYZ.f"
	    } else {

/*              One step of incremental condition estimation. */
/*              Complex workspace: need maximum 3*P. */

#line 568 "AG8BYZ.f"
		zcopy_(&rank, &abcd[irow + (icol + 1) * abcd_dim1], ldabcd, &
			zwork[jwork2], &c__1);
#line 570 "AG8BYZ.f"
		zlaic1_(&c__2, &rank, &zwork[ismin], &smin, &zwork[jwork2], &
			abcd[irow + icol * abcd_dim1], &sminpr, &s1, &c1);
#line 573 "AG8BYZ.f"
		zlaic1_(&c__1, &rank, &zwork[ismax], &smax, &zwork[jwork2], &
			abcd[irow + icol * abcd_dim1], &smaxpr, &s2, &c2);
/* Computing MAX */
#line 576 "AG8BYZ.f"
		i__1 = wrkopt, i__2 = *p * 3;
#line 576 "AG8BYZ.f"
		wrkopt = max(i__1,i__2);
#line 577 "AG8BYZ.f"
	    }

/*           Check the rank; finish the loop if rank loss occurs. */

#line 581 "AG8BYZ.f"
	    if (*svlmax * rcond <= smaxpr) {
#line 582 "AG8BYZ.f"
		if (*svlmax * rcond <= sminpr) {
#line 583 "AG8BYZ.f"
		    if (smaxpr * rcond <= sminpr) {

/*                    Finish the loop if last row. */

#line 587 "AG8BYZ.f"
			if (n1 == 0) {
#line 588 "AG8BYZ.f"
			    ++rank;
#line 589 "AG8BYZ.f"
			    goto L80;
#line 590 "AG8BYZ.f"
			}

#line 592 "AG8BYZ.f"
			if (n1 > 1) {

/*                       Update norms. */

#line 596 "AG8BYZ.f"
			    if (i__ - 1 > 1) {
#line 597 "AG8BYZ.f"
				i__1 = i__ - 1;
#line 597 "AG8BYZ.f"
				for (j = 1; j <= i__1; ++j) {
#line 598 "AG8BYZ.f"
				    if (dwork[j] != 0.) {
/* Computing 2nd power */
#line 599 "AG8BYZ.f"
					d__1 = z_abs(&abcd[irc + j + icol * 
						abcd_dim1]) / dwork[j];
#line 599 "AG8BYZ.f"
					t = 1. - d__1 * d__1;
#line 601 "AG8BYZ.f"
					t = max(t,0.);
/* Computing 2nd power */
#line 602 "AG8BYZ.f"
					d__1 = dwork[j] / dwork[*p + j];
#line 602 "AG8BYZ.f"
					tt = t * .05 * (d__1 * d__1) + 1.;
#line 604 "AG8BYZ.f"
					if (tt != 1.) {
#line 605 "AG8BYZ.f"
					    dwork[j] *= sqrt(t);
#line 606 "AG8BYZ.f"
					} else {
#line 607 "AG8BYZ.f"
					    i__2 = n1 - 1;
#line 607 "AG8BYZ.f"
					    dwork[j] = dznrm2_(&i__2, &abcd[
						    irc + j + mp1 * abcd_dim1]
						    , ldabcd);
#line 609 "AG8BYZ.f"
					    dwork[*p + j] = dwork[j];
#line 610 "AG8BYZ.f"
					}
#line 611 "AG8BYZ.f"
				    }
#line 612 "AG8BYZ.f"
/* L60: */
#line 612 "AG8BYZ.f"
				}
#line 613 "AG8BYZ.f"
			    }
#line 614 "AG8BYZ.f"
			}

#line 616 "AG8BYZ.f"
			i__1 = rank;
#line 616 "AG8BYZ.f"
			for (j = 1; j <= i__1; ++j) {
#line 617 "AG8BYZ.f"
			    i__2 = ismin + j - 1;
#line 617 "AG8BYZ.f"
			    i__3 = ismin + j - 1;
#line 617 "AG8BYZ.f"
			    z__1.r = s1.r * zwork[i__3].r - s1.i * zwork[i__3]
				    .i, z__1.i = s1.r * zwork[i__3].i + s1.i *
				     zwork[i__3].r;
#line 617 "AG8BYZ.f"
			    zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
#line 618 "AG8BYZ.f"
			    i__2 = ismax + j - 1;
#line 618 "AG8BYZ.f"
			    i__3 = ismax + j - 1;
#line 618 "AG8BYZ.f"
			    z__1.r = s2.r * zwork[i__3].r - s2.i * zwork[i__3]
				    .i, z__1.i = s2.r * zwork[i__3].i + s2.i *
				     zwork[i__3].r;
#line 618 "AG8BYZ.f"
			    zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
#line 619 "AG8BYZ.f"
/* L70: */
#line 619 "AG8BYZ.f"
			}

#line 621 "AG8BYZ.f"
			i__1 = ismin + rank;
#line 621 "AG8BYZ.f"
			zwork[i__1].r = c1.r, zwork[i__1].i = c1.i;
#line 622 "AG8BYZ.f"
			i__1 = ismax + rank;
#line 622 "AG8BYZ.f"
			zwork[i__1].r = c2.r, zwork[i__1].i = c2.i;
#line 623 "AG8BYZ.f"
			smin = sminpr;
#line 624 "AG8BYZ.f"
			smax = smaxpr;
#line 625 "AG8BYZ.f"
			++rank;
#line 626 "AG8BYZ.f"
			--icol;
#line 627 "AG8BYZ.f"
			--irow;
#line 628 "AG8BYZ.f"
			--n1;
#line 629 "AG8BYZ.f"
			--i__;
#line 630 "AG8BYZ.f"
			goto L40;
#line 631 "AG8BYZ.f"
		    }
#line 632 "AG8BYZ.f"
		}
#line 633 "AG8BYZ.f"
	    }
#line 634 "AG8BYZ.f"
	}
#line 635 "AG8BYZ.f"
    }

#line 637 "AG8BYZ.f"
L80:
#line 638 "AG8BYZ.f"
    mui = rank;
#line 639 "AG8BYZ.f"
    *nr -= mui;
#line 640 "AG8BYZ.f"
    *pr = sigma + mui;

/*     Set number of left Kronecker blocks of order (i-1)-by-i. */

#line 644 "AG8BYZ.f"
    kronl[nblcks] = taui - mui;

/*     Set number of infinite divisors of order i-1. */

#line 648 "AG8BYZ.f"
    if (*first && nblcks > 1) {
#line 648 "AG8BYZ.f"
	infz[nblcks - 1] = muim1 - taui;
#line 648 "AG8BYZ.f"
    }
#line 650 "AG8BYZ.f"
    muim1 = mui;
#line 651 "AG8BYZ.f"
    ro = mui;

/*     Continue reduction if rank of current C is positive. */

#line 655 "AG8BYZ.f"
    if (mui > 0) {
#line 655 "AG8BYZ.f"
	goto L10;
#line 655 "AG8BYZ.f"
    }

/*     Determine the maximal degree of infinite zeros and */
/*     the number of infinite zeros. */

#line 661 "AG8BYZ.f"
L90:
#line 662 "AG8BYZ.f"
    if (*first) {
#line 663 "AG8BYZ.f"
	if (mui == 0) {
/* Computing MAX */
#line 664 "AG8BYZ.f"
	    i__1 = 0, i__2 = nblcks - 1;
#line 664 "AG8BYZ.f"
	    *dinfz = max(i__1,i__2);
#line 665 "AG8BYZ.f"
	} else {
#line 666 "AG8BYZ.f"
	    *dinfz = nblcks;
#line 667 "AG8BYZ.f"
	    infz[nblcks] = mui;
#line 668 "AG8BYZ.f"
	}
#line 669 "AG8BYZ.f"
	k = *dinfz;
#line 670 "AG8BYZ.f"
	for (i__ = k; i__ >= 1; --i__) {
#line 671 "AG8BYZ.f"
	    if (infz[i__] != 0) {
#line 671 "AG8BYZ.f"
		goto L110;
#line 671 "AG8BYZ.f"
	    }
#line 672 "AG8BYZ.f"
	    --(*dinfz);
#line 673 "AG8BYZ.f"
/* L100: */
#line 673 "AG8BYZ.f"
	}
#line 674 "AG8BYZ.f"
L110:
#line 675 "AG8BYZ.f"
	i__1 = *dinfz;
#line 675 "AG8BYZ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 676 "AG8BYZ.f"
	    *ninfz += infz[i__] * i__;
#line 677 "AG8BYZ.f"
/* L120: */
#line 677 "AG8BYZ.f"
	}
#line 678 "AG8BYZ.f"
    }

/*     Determine the maximal order of left elementary Kronecker blocks. */

#line 682 "AG8BYZ.f"
    *nkronl = nblcks;
#line 683 "AG8BYZ.f"
    for (i__ = nblcks; i__ >= 1; --i__) {
#line 684 "AG8BYZ.f"
	if (kronl[i__] != 0) {
#line 684 "AG8BYZ.f"
	    goto L140;
#line 684 "AG8BYZ.f"
	}
#line 685 "AG8BYZ.f"
	--(*nkronl);
#line 686 "AG8BYZ.f"
/* L130: */
#line 686 "AG8BYZ.f"
    }
#line 687 "AG8BYZ.f"
L140:

#line 689 "AG8BYZ.f"
    zwork[1].r = (doublereal) wrkopt, zwork[1].i = 0.;
#line 690 "AG8BYZ.f"
    return 0;
/* *** Last line of AG8BYZ *** */
} /* ag8byz_ */

