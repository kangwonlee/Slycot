#line 1 "TD03AD.f"
/* TD03AD.f -- translated by f2c (version 20100827).
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

#line 1 "TD03AD.f"
/* Table of constant values */

static doublereal c_b12 = 0.;
static integer c__1 = 1;

/* Subroutine */ int td03ad_(char *rowcol, char *leri, char *equil, integer *
	m, integer *p, integer *indexd, doublereal *dcoeff, integer *lddcoe, 
	doublereal *ucoeff, integer *lduco1, integer *lduco2, integer *nr, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *d__, integer *ldd, integer *indexp, 
	doublereal *pcoeff, integer *ldpco1, integer *ldpco2, doublereal *
	qcoeff, integer *ldqco1, integer *ldqco2, doublereal *vcoeff, integer 
	*ldvco1, integer *ldvco2, doublereal *tol, integer *iwork, doublereal 
	*dwork, integer *ldwork, integer *info, ftnlen rowcol_len, ftnlen 
	leri_len, ftnlen equil_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, dcoeff_dim1, dcoeff_offset, pcoeff_dim1, pcoeff_dim2, 
	    pcoeff_offset, qcoeff_dim1, qcoeff_dim2, qcoeff_offset, 
	    ucoeff_dim1, ucoeff_dim2, ucoeff_offset, vcoeff_dim1, vcoeff_dim2,
	     vcoeff_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, k, n;
    extern /* Subroutine */ int ab07md_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen), tb03ad_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen), tc01od_(char 
	    *, integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, ftnlen);
    static integer idual;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01xd_(char *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen), td03ay_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *);
    static logical lleri;
    static integer itemp, mplim;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer maxmp, jstop, mwork, pwork, kdcoef, kpcoef;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical lequil, lrowco;


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

/*     To find a relatively prime left or right polynomial matrix */
/*     representation for a proper transfer matrix T(s) given as either */
/*     row or column polynomial vectors over common denominator */
/*     polynomials, possibly with uncancelled common terms. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     ROWCOL  CHARACTER*1 */
/*             Indicates whether T(s) is to be factorized by rows or by */
/*             columns as follows: */
/*             = 'R':  T(s) is factorized by rows; */
/*             = 'C':  T(s) is factorized by columns. */

/*     LERI    CHARACTER*1 */
/*             Indicates whether a left or a right polynomial matrix */
/*             representation is required as follows: */
/*             = 'L':  A left polynomial matrix representation */
/*                     inv(P(s))*Q(s) is required; */
/*             = 'R':  A right polynomial matrix representation */
/*                     Q(s)*inv(P(s)) is required. */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to balance the triplet */
/*             (A,B,C), before computing a minimal state-space */
/*             representation, as follows: */
/*             = 'S':  Perform balancing (scaling); */
/*             = 'N':  Do not perform balancing. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     INDEXD  (input) INTEGER array, dimension (P), if ROWCOL = 'R', or */
/*                                    dimension (M), if ROWCOL = 'C'. */
/*             The leading pormd elements of this array must contain the */
/*             row degrees of the denominator polynomials in D(s). */
/*             pormd = P if the transfer matrix T(s) is given as row */
/*             polynomial vectors over denominator polynomials; */
/*             pormd = M if the transfer matrix T(s) is given as column */
/*             polynomial vectors over denominator polynomials. */

/*     DCOEFF  (input) DOUBLE PRECISION array, dimension (LDDCOE,kdcoef), */
/*             where kdcoef = MAX(INDEXD(I)) + 1. */
/*             The leading pormd-by-kdcoef part of this array must */
/*             contain the coefficients of each denominator polynomial. */
/*             DCOEFF(I,K) is the coefficient in s**(INDEXD(I)-K+1) of */
/*             the I-th denominator polynomial in D(s), where K = 1,2, */
/*             ...,kdcoef. */

/*     LDDCOE  INTEGER */
/*             The leading dimension of array DCOEFF. */
/*             LDDCOE >= MAX(1,P), if ROWCOL = 'R'; */
/*             LDDCOE >= MAX(1,M), if ROWCOL = 'C'. */

/*     UCOEFF  (input) DOUBLE PRECISION array, dimension */
/*             (LDUCO1,LDUCO2,kdcoef) */
/*             The leading P-by-M-by-kdcoef part of this array must */
/*             contain the coefficients of the numerator matrix U(s); */
/*             if ROWCOL = 'C', this array is modified internally but */
/*             restored on exit, and the remainder of the leading */
/*             MAX(M,P)-by-MAX(M,P)-by-kdcoef part is used as internal */
/*             workspace. */
/*             UCOEFF(I,J,K) is the coefficient in s**(INDEXD(iorj)-K+1) */
/*             of polynomial (I,J) of U(s), where K = 1,2,...,kdcoef; */
/*             iorj = I if T(s) is given as row polynomial vectors over */
/*             denominator polynomials; iorj = J if T(s) is given as */
/*             column polynomial vectors over denominator polynomials. */
/*             Thus for ROWCOL = 'R', U(s) = */
/*             diag(s**INDEXD(I))*(UCOEFF(.,.,1)+UCOEFF(.,.,2)/s+...). */

/*     LDUCO1  INTEGER */
/*             The leading dimension of array UCOEFF. */
/*             LDUCO1 >= MAX(1,P),   if ROWCOL = 'R'; */
/*             LDUCO1 >= MAX(1,M,P), if ROWCOL = 'C'. */

/*     LDUCO2  INTEGER */
/*             The second dimension of array UCOEFF. */
/*             LDUCO2 >= MAX(1,M),   if ROWCOL = 'R'; */
/*             LDUCO2 >= MAX(1,M,P), if ROWCOL = 'C'. */

/*     NR      (output) INTEGER */
/*             The order of the resulting minimal realization, i.e. the */
/*             order of the state dynamics matrix A. */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N), */
/*                      pormd */
/*             where N = SUM INDEXD(I) */
/*                       I=1 */
/*             The leading NR-by-NR part of this array contains the upper */
/*             block Hessenberg state dynamics matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,MAX(M,P)) */
/*             The leading NR-by-M part of this array contains the */
/*             input/state matrix B; the remainder of the leading */
/*             N-by-MAX(M,P) part is used as internal workspace. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-NR part of this array contains the */
/*             state/output matrix C; the remainder of the leading */
/*             MAX(M,P)-by-N part is used as internal workspace. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,M,P). */

/*     D       (output) DOUBLE PRECISION array, dimension (LDD,MAX(M,P)) */
/*             The leading P-by-M part of this array contains the direct */
/*             transmission matrix D; the remainder of the leading */
/*             MAX(M,P)-by-MAX(M,P) part is used as internal workspace. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,M,P). */

/*     INDEXP  (output) INTEGER array, dimension (P), if ROWCOL = 'R', or */
/*                                     dimension (M), if ROWCOL = 'C'. */
/*             The leading pormp elements of this array contain the */
/*             row (column if ROWCOL = 'C') degrees of the denominator */
/*             matrix P(s). */
/*             pormp = P if a left polynomial matrix representation */
/*             is requested; pormp = M if a right polynomial matrix */
/*             representation is requested. */
/*             These elements are ordered so that */
/*             INDEXP(1) >= INDEXP(2) >= ... >= INDEXP(pormp). */

/*     PCOEFF  (output) DOUBLE PRECISION array, dimension */
/*             (LDPCO1,LDPCO2,N+1) */
/*             The leading pormp-by-pormp-by-kpcoef part of this array */
/*             contains the coefficients of the denominator matrix P(s), */
/*             where kpcoef = MAX(INDEXP(I)) + 1. */
/*             PCOEFF(I,J,K) is the coefficient in s**(INDEXP(iorj)-K+1) */
/*             of polynomial (I,J) of P(s), where K = 1,2,...,kpcoef; */
/*             iorj = I if a left polynomial matrix representation is */
/*             requested; iorj = J if a right polynomial matrix */
/*             representation is requested. */
/*             Thus for a left polynomial matrix representation, P(s) = */
/*             diag(s**INDEXP(I))*(PCOEFF(.,.,1)+PCOEFF(.,.,2)/s+...). */

/*     LDPCO1  INTEGER */
/*             The leading dimension of array PCOEFF. */
/*             LDPCO1 >= MAX(1,P), if ROWCOL = 'R'; */
/*             LDPCO1 >= MAX(1,M), if ROWCOL = 'C'. */

/*     LDPCO2  INTEGER */
/*             The second dimension of array PCOEFF. */
/*             LDPCO2 >= MAX(1,P), if ROWCOL = 'R'; */
/*             LDPCO2 >= MAX(1,M), if ROWCOL = 'C'. */

/*     QCOEFF  (output) DOUBLE PRECISION array, dimension */
/*             (LDQCO1,LDQCO2,N+1) */
/*             The leading pormp-by-pormd-by-kpcoef part of this array */
/*             contains the coefficients of the numerator matrix Q(s). */
/*             QCOEFF(I,J,K) is defined as for PCOEFF(I,J,K). */

/*     LDQCO1  INTEGER */
/*             The leading dimension of array QCOEFF. */
/*             If LERI = 'L', LDQCO1 >= MAX(1,PM), */
/*                                      where PM = P, if ROWCOL = 'R'; */
/*                                            PM = M, if ROWCOL = 'C'. */
/*             If LERI = 'R', LDQCO1 >= MAX(1,M,P). */

/*     LDQCO2  INTEGER */
/*             The second dimension of array QCOEFF. */
/*             If LERI = 'L', LDQCO2 >= MAX(1,MP), */
/*                                      where MP = M, if ROWCOL = 'R'; */
/*                                            MP = P, if ROWCOL = 'C'. */
/*             If LERI = 'R', LDQCO2 >= MAX(1,M,P). */

/*     VCOEFF  (output) DOUBLE PRECISION array, dimension */
/*             (LDVCO1,LDVCO2,N+1) */
/*             The leading pormp-by-NR-by-kpcoef part of this array */
/*             contains the coefficients of the intermediate matrix */
/*             V(s) as produced by SLICOT Library routine TB03AD. */

/*     LDVCO1  INTEGER */
/*             The leading dimension of array VCOEFF. */
/*             LDVCO1 >= MAX(1,P), if ROWCOL = 'R'; */
/*             LDVCO1 >= MAX(1,M), if ROWCOL = 'C'. */

/*     LDVCO2  INTEGER */
/*             The second dimension of array VCOEFF.  LDVCO2 >= MAX(1,N). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in rank determination when */
/*             transforming (A, B, C). If the user sets TOL > 0, then */
/*             the given value of TOL is used as a lower bound for the */
/*             reciprocal condition number (see the description of the */
/*             argument RCOND in the SLICOT routine MB03OD);  a */
/*             (sub)matrix whose estimated condition number is less than */
/*             1/TOL is considered to be of full rank.  If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance */
/*             (determined by the SLICOT routine TB01UD) is used instead. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N+MAX(M,P)) */
/*             On exit, if INFO = 0, the first nonzero elements of */
/*             IWORK(1:N) return the orders of the diagonal blocks of A. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1, N + MAX(N, 3*M, 3*P), PM*(PM + 2)) */
/*             where  PM = P, if ROWCOL = 'R'; */
/*                    PM = M, if ROWCOL = 'C'. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i (i <= k = pormd), then i is the first */
/*                   integer I for which ABS( DCOEFF(I,1) ) is so small */
/*                   that the calculations would overflow (see SLICOT */
/*                   Library routine TD03AY); that is, the leading */
/*                   coefficient of a polynomial is nearly zero; no */
/*                   state-space representation or polynomial matrix */
/*                   representation is calculated; */
/*             = k+1:  if a singular matrix was encountered during the */
/*                   computation of V(s); */
/*             = k+2:  if a singular matrix was encountered during the */
/*                   computation of P(s). */

/*     METHOD */

/*     The method for transfer matrices factorized by rows will be */
/*     described here; T(s) factorized by columns is dealt with by */
/*     operating on the dual T'(s). The description for T(s) is actually */
/*     the left polynomial matrix representation */

/*          T(s) = inv(D(s))*U(s), */

/*     where D(s) is diagonal with its (I,I)-th polynomial element of */
/*     degree INDEXD(I). The first step is to check whether the leading */
/*     coefficient of any polynomial element of D(s) is approximately */
/*     zero, if so the routine returns with INFO > 0. Otherwise, */
/*     Wolovich's Observable Structure Theorem is used to construct a */
/*     state-space representation in observable companion form which is */
/*     equivalent to the above polynomial matrix representation. The */
/*     method is particularly easy here due to the diagonal form of D(s). */
/*     This state-space representation is not necessarily controllable */
/*     (as D(s) and U(s) are not necessarily relatively left prime), but */
/*     it is in theory completely observable; however, its observability */
/*     matrix may be poorly conditioned, so it is treated as a general */
/*     state-space representation and SLICOT Library routine TB03AD is */
/*     used to separate out a minimal realization for T(s) from it by */
/*     means of orthogonal similarity transformations and then to */
/*     calculate a relatively prime (left or right) polynomial matrix */
/*     representation which is equivalent to this. */

/*     REFERENCES */

/*     [1] Patel, R.V. */
/*         On Computing Matrix Fraction Descriptions and Canonical */
/*         Forms of Linear Time-Invariant Systems. */
/*         UMIST Control Systems Centre Report 489, 1980. */

/*     [2] Wolovich, W.A. */
/*         Linear Multivariable Systems, (Theorem 4.3.3). */
/*         Springer-Verlag, 1974. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1998. */
/*     Supersedes Release 3.0 routine TD01ND. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Coprime matrix fraction, elementary polynomial operations, */
/*     polynomial matrix, state-space representation, transfer matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 355 "TD03AD.f"
    /* Parameter adjustments */
#line 355 "TD03AD.f"
    --indexd;
#line 355 "TD03AD.f"
    dcoeff_dim1 = *lddcoe;
#line 355 "TD03AD.f"
    dcoeff_offset = 1 + dcoeff_dim1;
#line 355 "TD03AD.f"
    dcoeff -= dcoeff_offset;
#line 355 "TD03AD.f"
    ucoeff_dim1 = *lduco1;
#line 355 "TD03AD.f"
    ucoeff_dim2 = *lduco2;
#line 355 "TD03AD.f"
    ucoeff_offset = 1 + ucoeff_dim1 * (1 + ucoeff_dim2);
#line 355 "TD03AD.f"
    ucoeff -= ucoeff_offset;
#line 355 "TD03AD.f"
    a_dim1 = *lda;
#line 355 "TD03AD.f"
    a_offset = 1 + a_dim1;
#line 355 "TD03AD.f"
    a -= a_offset;
#line 355 "TD03AD.f"
    b_dim1 = *ldb;
#line 355 "TD03AD.f"
    b_offset = 1 + b_dim1;
#line 355 "TD03AD.f"
    b -= b_offset;
#line 355 "TD03AD.f"
    c_dim1 = *ldc;
#line 355 "TD03AD.f"
    c_offset = 1 + c_dim1;
#line 355 "TD03AD.f"
    c__ -= c_offset;
#line 355 "TD03AD.f"
    d_dim1 = *ldd;
#line 355 "TD03AD.f"
    d_offset = 1 + d_dim1;
#line 355 "TD03AD.f"
    d__ -= d_offset;
#line 355 "TD03AD.f"
    --indexp;
#line 355 "TD03AD.f"
    pcoeff_dim1 = *ldpco1;
#line 355 "TD03AD.f"
    pcoeff_dim2 = *ldpco2;
#line 355 "TD03AD.f"
    pcoeff_offset = 1 + pcoeff_dim1 * (1 + pcoeff_dim2);
#line 355 "TD03AD.f"
    pcoeff -= pcoeff_offset;
#line 355 "TD03AD.f"
    qcoeff_dim1 = *ldqco1;
#line 355 "TD03AD.f"
    qcoeff_dim2 = *ldqco2;
#line 355 "TD03AD.f"
    qcoeff_offset = 1 + qcoeff_dim1 * (1 + qcoeff_dim2);
#line 355 "TD03AD.f"
    qcoeff -= qcoeff_offset;
#line 355 "TD03AD.f"
    vcoeff_dim1 = *ldvco1;
#line 355 "TD03AD.f"
    vcoeff_dim2 = *ldvco2;
#line 355 "TD03AD.f"
    vcoeff_offset = 1 + vcoeff_dim1 * (1 + vcoeff_dim2);
#line 355 "TD03AD.f"
    vcoeff -= vcoeff_offset;
#line 355 "TD03AD.f"
    --iwork;
#line 355 "TD03AD.f"
    --dwork;
#line 355 "TD03AD.f"

#line 355 "TD03AD.f"
    /* Function Body */
#line 355 "TD03AD.f"
    *info = 0;
#line 356 "TD03AD.f"
    lrowco = lsame_(rowcol, "R", (ftnlen)1, (ftnlen)1);
#line 357 "TD03AD.f"
    lleri = lsame_(leri, "L", (ftnlen)1, (ftnlen)1);
#line 358 "TD03AD.f"
    lequil = lsame_(equil, "S", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 362 "TD03AD.f"
    maxmp = max(*m,*p);
#line 363 "TD03AD.f"
    mplim = max(1,maxmp);
#line 364 "TD03AD.f"
    if (lrowco) {

/*        Initialization for T(s) given as rows over common denominators. */

#line 368 "TD03AD.f"
	pwork = *p;
#line 369 "TD03AD.f"
	mwork = *m;
#line 370 "TD03AD.f"
    } else {

/*        Initialization for T(s) given as columns over common */
/*        denominators. */

#line 375 "TD03AD.f"
	pwork = *m;
#line 376 "TD03AD.f"
	mwork = *p;
#line 377 "TD03AD.f"
    }

#line 379 "TD03AD.f"
    if (! lrowco && ! lsame_(rowcol, "C", (ftnlen)1, (ftnlen)1)) {
#line 380 "TD03AD.f"
	*info = -1;
#line 381 "TD03AD.f"
    } else if (! lleri && ! lsame_(leri, "R", (ftnlen)1, (ftnlen)1)) {
#line 382 "TD03AD.f"
	*info = -2;
#line 383 "TD03AD.f"
    } else if (! lequil && ! lsame_(equil, "N", (ftnlen)1, (ftnlen)1)) {
#line 384 "TD03AD.f"
	*info = -3;
#line 385 "TD03AD.f"
    } else if (*m < 0) {
#line 386 "TD03AD.f"
	*info = -4;
#line 387 "TD03AD.f"
    } else if (*p < 0) {
#line 388 "TD03AD.f"
	*info = -5;
#line 389 "TD03AD.f"
    } else if (*lddcoe < max(1,pwork)) {
#line 390 "TD03AD.f"
	*info = -8;
#line 391 "TD03AD.f"
    } else if (*lduco1 < max(1,pwork) || ! lrowco && *lduco1 < mplim) {
#line 393 "TD03AD.f"
	*info = -10;
#line 394 "TD03AD.f"
    } else if (*lduco2 < max(1,mwork) || ! lrowco && *lduco2 < mplim) {
#line 396 "TD03AD.f"
	*info = -11;
#line 397 "TD03AD.f"
    }

#line 399 "TD03AD.f"
    n = 0;
#line 400 "TD03AD.f"
    if (*info == 0) {

/*        Calculate N, the order of the resulting state-space */
/*        representation, and the index kdcoef. */

#line 405 "TD03AD.f"
	kdcoef = 0;

#line 407 "TD03AD.f"
	i__1 = pwork;
#line 407 "TD03AD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 408 "TD03AD.f"
	    i__2 = kdcoef, i__3 = indexd[i__];
#line 408 "TD03AD.f"
	    kdcoef = max(i__2,i__3);
#line 409 "TD03AD.f"
	    n += indexd[i__];
#line 410 "TD03AD.f"
/* L10: */
#line 410 "TD03AD.f"
	}

#line 412 "TD03AD.f"
	++kdcoef;

#line 414 "TD03AD.f"
	if (*lda < max(1,n)) {
#line 415 "TD03AD.f"
	    *info = -14;
#line 416 "TD03AD.f"
	} else if (*ldb < max(1,n)) {
#line 417 "TD03AD.f"
	    *info = -16;
#line 418 "TD03AD.f"
	} else if (*ldc < mplim) {
#line 419 "TD03AD.f"
	    *info = -18;
#line 420 "TD03AD.f"
	} else if (*ldd < mplim) {
#line 421 "TD03AD.f"
	    *info = -20;
#line 422 "TD03AD.f"
	} else if (*ldpco1 < pwork) {
#line 423 "TD03AD.f"
	    *info = -23;
#line 424 "TD03AD.f"
	} else if (*ldpco2 < pwork) {
#line 425 "TD03AD.f"
	    *info = -24;
#line 426 "TD03AD.f"
	} else if (*ldqco1 < max(1,pwork) || ! lleri && *ldqco1 < mplim) {
#line 428 "TD03AD.f"
	    *info = -26;
#line 429 "TD03AD.f"
	} else if (*ldqco2 < max(1,mwork) || ! lleri && *ldqco2 < mplim) {
#line 431 "TD03AD.f"
	    *info = -27;
#line 432 "TD03AD.f"
	} else if (*ldvco1 < max(1,pwork)) {
#line 433 "TD03AD.f"
	    *info = -29;
#line 434 "TD03AD.f"
	} else if (*ldvco2 < max(1,n)) {
#line 435 "TD03AD.f"
	    *info = -30;

#line 437 "TD03AD.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
#line 437 "TD03AD.f"
	    i__3 = n, i__4 = maxmp * 3;
#line 437 "TD03AD.f"
	    i__1 = 1, i__2 = n + max(i__3,i__4), i__1 = max(i__1,i__2), i__2 =
		     pwork * (pwork + 2);
#line 437 "TD03AD.f"
	    if (*ldwork < max(i__1,i__2)) {
#line 439 "TD03AD.f"
		*info = -34;
#line 440 "TD03AD.f"
	    }
#line 440 "TD03AD.f"
	}
#line 441 "TD03AD.f"
    }

#line 443 "TD03AD.f"
    if (*info != 0) {

/*        Error return. */

#line 447 "TD03AD.f"
	i__1 = -(*info);
#line 447 "TD03AD.f"
	xerbla_("TD03AD", &i__1, (ftnlen)6);
#line 448 "TD03AD.f"
	return 0;
#line 449 "TD03AD.f"
    }

/*     Quick return if possible. */

/* Computing MAX */
#line 453 "TD03AD.f"
    i__1 = max(n,*m);
#line 453 "TD03AD.f"
    if (max(i__1,*p) == 0) {
#line 454 "TD03AD.f"
	*nr = 0;
#line 455 "TD03AD.f"
	dwork[1] = 1.;
#line 456 "TD03AD.f"
	return 0;
#line 457 "TD03AD.f"
    }

/*     IDUAL = 1 iff precisely ROWCOL = 'C' or (exclusively) LERI = 'R', */
/*     i.e. iff AB07MD call is required before TB03AD. */

#line 462 "TD03AD.f"
    idual = 0;
#line 463 "TD03AD.f"
    if (! lrowco) {
#line 463 "TD03AD.f"
	idual = 1;
#line 463 "TD03AD.f"
    }
#line 464 "TD03AD.f"
    if (! lleri) {
#line 464 "TD03AD.f"
	++idual;
#line 464 "TD03AD.f"
    }

#line 466 "TD03AD.f"
    if (! lrowco) {

/*        Initialize the remainder of the leading */
/*        MPLIM-by-MPLIM-by-KDCOEF part of U(s) to zero. */

#line 471 "TD03AD.f"
	if (*p < *m) {

#line 473 "TD03AD.f"
	    i__1 = kdcoef;
#line 473 "TD03AD.f"
	    for (k = 1; k <= i__1; ++k) {
#line 474 "TD03AD.f"
		i__2 = *m - *p;
#line 474 "TD03AD.f"
		dlacpy_("Full", &i__2, &mplim, &c_b12, &c_b12, &ucoeff[*p + 1 
			+ (k * ucoeff_dim2 + 1) * ucoeff_dim1], lduco1, (
			ftnlen)4);
#line 476 "TD03AD.f"
/* L20: */
#line 476 "TD03AD.f"
	    }

#line 478 "TD03AD.f"
	} else if (*p > *m) {

#line 480 "TD03AD.f"
	    i__1 = kdcoef;
#line 480 "TD03AD.f"
	    for (k = 1; k <= i__1; ++k) {
#line 481 "TD03AD.f"
		i__2 = *p - *m;
#line 481 "TD03AD.f"
		dlacpy_("Full", &mplim, &i__2, &c_b12, &c_b12, &ucoeff[(*m + 
			1 + k * ucoeff_dim2) * ucoeff_dim1 + 1], lduco1, (
			ftnlen)4);
#line 483 "TD03AD.f"
/* L30: */
#line 483 "TD03AD.f"
	    }

#line 485 "TD03AD.f"
	}

#line 487 "TD03AD.f"
	if (mplim != 1) {

/*           Non-scalar T(s) factorized by columns: transpose it */
/*           (i.e. U(s)). */

#line 492 "TD03AD.f"
	    jstop = mplim - 1;

#line 494 "TD03AD.f"
	    i__1 = kdcoef;
#line 494 "TD03AD.f"
	    for (k = 1; k <= i__1; ++k) {

#line 496 "TD03AD.f"
		i__2 = jstop;
#line 496 "TD03AD.f"
		for (j = 1; j <= i__2; ++j) {
#line 497 "TD03AD.f"
		    i__3 = mplim - j;
#line 497 "TD03AD.f"
		    dswap_(&i__3, &ucoeff[j + 1 + (j + k * ucoeff_dim2) * 
			    ucoeff_dim1], &c__1, &ucoeff[j + (j + 1 + k * 
			    ucoeff_dim2) * ucoeff_dim1], lduco1);
#line 499 "TD03AD.f"
/* L40: */
#line 499 "TD03AD.f"
		}

#line 501 "TD03AD.f"
/* L50: */
#line 501 "TD03AD.f"
	    }

#line 503 "TD03AD.f"
	}
#line 504 "TD03AD.f"
    }

/*     Construct non-minimal state-space representation (by Wolovich's */
/*     Structure Theorem) which has transfer matrix T(s) or T'(s) as */
/*     appropriate, */

#line 510 "TD03AD.f"
    td03ay_(&mwork, &pwork, &indexd[1], &dcoeff[dcoeff_offset], lddcoe, &
	    ucoeff[ucoeff_offset], lduco1, lduco2, &n, &a[a_offset], lda, &b[
	    b_offset], ldb, &c__[c_offset], ldc, &d__[d_offset], ldd, info);
#line 512 "TD03AD.f"
    if (*info > 0) {
#line 512 "TD03AD.f"
	return 0;
#line 512 "TD03AD.f"
    }

#line 515 "TD03AD.f"
    if (idual == 1) {

/*        and then obtain (MWORK x PWORK) dual of this system if */
/*        appropriate. */

#line 520 "TD03AD.f"
	ab07md_("D", &n, &mwork, &pwork, &a[a_offset], lda, &b[b_offset], ldb,
		 &c__[c_offset], ldc, &d__[d_offset], ldd, info, (ftnlen)1);
#line 522 "TD03AD.f"
	itemp = pwork;
#line 523 "TD03AD.f"
	pwork = mwork;
#line 524 "TD03AD.f"
	mwork = itemp;
#line 525 "TD03AD.f"
    }

/*     Find left polynomial matrix representation (and minimal */
/*     state-space representation en route) for the relevant state-space */
/*     representation ... */

#line 531 "TD03AD.f"
    tb03ad_("Left", equil, &n, &mwork, &pwork, &a[a_offset], lda, &b[b_offset]
	    , ldb, &c__[c_offset], ldc, &d__[d_offset], ldd, nr, &indexp[1], &
	    pcoeff[pcoeff_offset], ldpco1, ldpco2, &qcoeff[qcoeff_offset], 
	    ldqco1, ldqco2, &vcoeff[vcoeff_offset], ldvco1, ldvco2, tol, &
	    iwork[1], &dwork[1], ldwork, info, (ftnlen)4, (ftnlen)1);

#line 536 "TD03AD.f"
    if (*info > 0) {
#line 537 "TD03AD.f"
	*info = pwork + *info;
#line 538 "TD03AD.f"
	return 0;
#line 539 "TD03AD.f"
    }

#line 541 "TD03AD.f"
    if (! lleri) {

/*        and, if a right polynomial matrix representation is required, */
/*        transpose and reorder (to get a block upper Hessenberg */
/*        matrix A). */

#line 547 "TD03AD.f"
	k = iwork[1] - 1;
#line 548 "TD03AD.f"
	if (n >= 2) {
#line 548 "TD03AD.f"
	    k += iwork[2];
#line 548 "TD03AD.f"
	}
#line 550 "TD03AD.f"
	i__1 = *nr - 1;
#line 550 "TD03AD.f"
	tb01xd_("D", nr, &mwork, &pwork, &k, &i__1, &a[a_offset], lda, &b[
		b_offset], ldb, &c__[c_offset], ldc, &d__[d_offset], ldd, 
		info, (ftnlen)1);

#line 553 "TD03AD.f"
	kpcoef = 0;

#line 555 "TD03AD.f"
	i__1 = pwork;
#line 555 "TD03AD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 556 "TD03AD.f"
	    i__2 = kpcoef, i__3 = indexp[i__];
#line 556 "TD03AD.f"
	    kpcoef = max(i__2,i__3);
#line 557 "TD03AD.f"
/* L60: */
#line 557 "TD03AD.f"
	}

#line 559 "TD03AD.f"
	++kpcoef;
#line 560 "TD03AD.f"
	tc01od_("L", &mwork, &pwork, &kpcoef, &pcoeff[pcoeff_offset], ldpco1, 
		ldpco2, &qcoeff[qcoeff_offset], ldqco1, ldqco2, info, (ftnlen)
		1);
#line 562 "TD03AD.f"
    }

#line 564 "TD03AD.f"
    if (! lrowco && mplim != 1) {

/*        If non-scalar T(s) originally given by columns, */
/*        retranspose U(s). */

#line 569 "TD03AD.f"
	i__1 = kdcoef;
#line 569 "TD03AD.f"
	for (k = 1; k <= i__1; ++k) {

#line 571 "TD03AD.f"
	    i__2 = jstop;
#line 571 "TD03AD.f"
	    for (j = 1; j <= i__2; ++j) {
#line 572 "TD03AD.f"
		i__3 = mplim - j;
#line 572 "TD03AD.f"
		dswap_(&i__3, &ucoeff[j + 1 + (j + k * ucoeff_dim2) * 
			ucoeff_dim1], &c__1, &ucoeff[j + (j + 1 + k * 
			ucoeff_dim2) * ucoeff_dim1], lduco1);
#line 574 "TD03AD.f"
/* L70: */
#line 574 "TD03AD.f"
	    }

#line 576 "TD03AD.f"
/* L80: */
#line 576 "TD03AD.f"
	}

#line 578 "TD03AD.f"
    }
#line 579 "TD03AD.f"
    return 0;
/* *** Last line of TD03AD *** */
} /* td03ad_ */
