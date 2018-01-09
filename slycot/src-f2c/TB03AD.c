#line 1 "TB03AD.f"
/* TB03AD.f -- translated by f2c (version 20100827).
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

#line 1 "TB03AD.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b20 = 0.;
static doublereal c_b53 = 1.;
static integer c_n1 = -1;

/* Subroutine */ int tb03ad_(char *leri, char *equil, integer *n, integer *m, 
	integer *p, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, integer 
	*nr, integer *index, doublereal *pcoeff, integer *ldpco1, integer *
	ldpco2, doublereal *qcoeff, integer *ldqco1, integer *ldqco2, 
	doublereal *vcoeff, integer *ldvco1, integer *ldvco2, doublereal *tol,
	 integer *iwork, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen leri_len, ftnlen equil_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, pcoeff_dim1, pcoeff_dim2, pcoeff_offset, qcoeff_dim1, 
	    qcoeff_dim2, qcoeff_offset, vcoeff_dim1, vcoeff_dim2, 
	    vcoeff_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, ic, iz, ioff, joff, ncol, kmax, itau, nrow;
    extern /* Subroutine */ int ab07md_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen), ma02gd_(
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *), tb01id_(char *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), dgemm_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen), tc01od_(char *, integer 
	    *, integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, ftnlen), tb01ud_(
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb03ay_(integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *), tb01yd_(integer *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *);
    static integer mplim, ncont, maxmp;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer jwork, kwork, istop, kplus, mwork, pwork, indblk, irankc, 
	    kpcoef, nreflc;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dgetrf_(integer *, integer *, doublereal *, integer *, integer *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), dlaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen), xerbla_(char *, integer *, ftnlen);
    static doublereal maxred;
    static logical lleril;
    static integer ldwric;
    static logical llerir, lequil;
    static integer ifirst;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer istart, inplus, wrkopt;


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

/*     To find a relatively prime left polynomial matrix representation */
/*     inv(P(s))*Q(s) or right polynomial matrix representation */
/*     Q(s)*inv(P(s)) with the same transfer matrix T(s) as that of a */
/*     given state-space representation, i.e. */

/*        inv(P(s))*Q(s) = Q(s)*inv(P(s)) = T(s) = C*inv(s*I-A)*B + D. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     LERI    CHARACTER*1 */
/*             Indicates whether the left polynomial matrix */
/*             representation or the right polynomial matrix */
/*             representation is required as follows: */
/*             = 'L':  A left matrix fraction is required; */
/*             = 'R':  A right matrix fraction is required. */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to balance the triplet */
/*             (A,B,C), before computing a minimal state-space */
/*             representation, as follows: */
/*             = 'S':  Perform balancing (scaling); */
/*             = 'N':  Do not perform balancing. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the state-space representation, i.e. the */
/*             order of the original state dynamics matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original state dynamics matrix A. */
/*             On exit, the leading NR-by-NR part of this array contains */
/*             the upper block Hessenberg state dynamics matrix Amin of a */
/*             minimal realization for the original system. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDB,MAX(M,P)) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the original input/state matrix B; the remainder */
/*             of the leading N-by-MAX(M,P) part is used as internal */
/*             workspace. */
/*             On exit, the leading NR-by-M part of this array contains */
/*             the transformed input/state matrix Bmin. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the original state/output matrix C; the remainder */
/*             of the leading MAX(M,P)-by-N part is used as internal */
/*             workspace. */
/*             On exit, the leading P-by-NR part of this array contains */
/*             the transformed state/output matrix Cmin. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,M,P). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,MAX(M,P)) */
/*             The leading P-by-M part of this array must contain the */
/*             original direct transmission matrix D; the remainder of */
/*             the leading MAX(M,P)-by-MAX(M,P) part is used as internal */
/*             workspace. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,M,P). */

/*     NR      (output) INTEGER */
/*             The order of the minimal state-space representation */
/*             (Amin,Bmin,Cmin). */

/*     INDEX   (output) INTEGER array, dimension (P), if LERI = 'L', or */
/*                                     dimension (M), if LERI = 'R'. */
/*             If LERI = 'L', INDEX(I), I = 1,2,...,P, contains the */
/*             maximum degree of the polynomials in the I-th row of the */
/*             denominator matrix P(s) of the left polynomial matrix */
/*             representation. */
/*             These elements are ordered so that */
/*             INDEX(1) >= INDEX(2) >= ... >= INDEX(P). */
/*             If LERI = 'R', INDEX(I), I = 1,2,...,M, contains the */
/*             maximum degree of the polynomials in the I-th column of */
/*             the denominator matrix P(s) of the right polynomial */
/*             matrix representation. */
/*             These elements are ordered so that */
/*             INDEX(1) >= INDEX(2) >= ... >= INDEX(M). */

/*     PCOEFF  (output) DOUBLE PRECISION array, dimension */
/*             (LDPCO1,LDPCO2,N+1) */
/*             If LERI = 'L' then porm = P, otherwise porm = M. */
/*             The leading porm-by-porm-by-kpcoef part of this array */
/*             contains the coefficients of the denominator matrix P(s), */
/*             where kpcoef = MAX(INDEX(I)) + 1. */
/*             PCOEFF(I,J,K) is the coefficient in s**(INDEX(iorj)-K+1) */
/*             of polynomial (I,J) of P(s), where K = 1,2,...,kpcoef; if */
/*             LERI = 'L' then iorj = I, otherwise iorj = J. */
/*             Thus for LERI = 'L', P(s) = */
/*             diag(s**INDEX(I))*(PCOEFF(.,.,1)+PCOEFF(.,.,2)/s+...). */

/*     LDPCO1  INTEGER */
/*             The leading dimension of array PCOEFF. */
/*             LDPCO1 >= MAX(1,P), if LERI = 'L'; */
/*             LDPCO1 >= MAX(1,M), if LERI = 'R'. */

/*     LDPCO2  INTEGER */
/*             The second dimension of array PCOEFF. */
/*             LDPCO2 >= MAX(1,P), if LERI = 'L'; */
/*             LDPCO2 >= MAX(1,M), if LERI = 'R'. */

/*     QCOEFF  (output) DOUBLE PRECISION array, dimension */
/*             (LDQCO1,LDQCO2,N+1) */
/*             If LERI = 'L' then porp = M, otherwise porp = P. */
/*             If LERI = 'L', the leading porm-by-porp-by-kpcoef part */
/*             of this array contains the coefficients of the numerator */
/*             matrix Q(s). */
/*             If LERI = 'R', the leading porp-by-porm-by-kpcoef part */
/*             of this array contains the coefficients of the numerator */
/*             matrix Q(s). */
/*             QCOEFF(I,J,K) is defined as for PCOEFF(I,J,K). */

/*     LDQCO1  INTEGER */
/*             The leading dimension of array QCOEFF. */
/*             LDQCO1 >= MAX(1,P),   if LERI = 'L'; */
/*             LDQCO1 >= MAX(1,M,P), if LERI = 'R'. */

/*     LDQCO2  INTEGER */
/*             The second dimension of array QCOEFF. */
/*             LDQCO2 >= MAX(1,M),   if LERI = 'L'; */
/*             LDQCO2 >= MAX(1,M,P), if LERI = 'R'. */

/*     VCOEFF  (output) DOUBLE PRECISION array, dimension */
/*             (LDVCO1,LDVCO2,N+1) */
/*             The leading porm-by-NR-by-kpcoef part of this array */
/*             contains the coefficients of the intermediate matrix V(s). */
/*             VCOEFF(I,J,K) is defined as for PCOEFF(I,J,K). */

/*     LDVCO1  INTEGER */
/*             The leading dimension of array VCOEFF. */
/*             LDVCO1 >= MAX(1,P), if LERI = 'L'; */
/*             LDVCO1 >= MAX(1,M), if LERI = 'R'. */

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
/*             where  PM = P, if LERI = 'L'; */
/*                    PM = M, if LERI = 'R'. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if a singular matrix was encountered during the */
/*                   computation of V(s); */
/*             = 2:  if a singular matrix was encountered during the */
/*                   computation of P(s). */

/*     METHOD */

/*     The method for a left matrix fraction will be described here: */
/*     right matrix fractions are dealt with by constructing a left */
/*     fraction for the dual of the original system. The first step is to */
/*     obtain, by means of orthogonal similarity transformations, a */
/*     minimal state-space representation (Amin,Bmin,Cmin,D) for the */
/*     original system (A,B,C,D), where Amin is lower block Hessenberg */
/*     with all its superdiagonal blocks upper triangular and Cmin has */
/*     all but its first rank(C) columns zero.  The number and dimensions */
/*     of the blocks of Amin now immediately yield the row degrees of */
/*     P(s) with P(s) row proper: furthermore, the P-by-NR polynomial */
/*     matrix V(s) (playing a similar role to S(s) in Wolovich's */
/*     Structure Theorem) can be calculated a column block at a time, in */
/*     reverse order, from Amin. P(s) is then found as if it were the */
/*     O-th column block of V(s) (using Cmin as well as Amin), while */
/*     Q(s) = (V(s) * Bmin) + (P(s) * D). Finally, a special similarity */
/*     transformation is used to put Amin in an upper block Hessenberg */
/*     form. */

/*     REFERENCES */

/*     [1] Williams, T.W.C. */
/*         An Orthogonal Structure Theorem for Linear Systems. */
/*         Kingston Polytechnic Control Systems Research Group, */
/*         Internal Report 82/2, July 1982. */

/*     [2] Patel, R.V. */
/*         On Computing Matrix Fraction Descriptions and Canonical */
/*         Forms of Linear Time-Invariant Systems. */
/*         UMIST Control Systems Centre Report 489, 1980. */
/*         (Algorithms 1 and 2, extensively modified). */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, March 1998. */
/*     Supersedes Release 3.0 routine TB01SD. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2000. */

/*     KEYWORDS */

/*     Canonical form, coprime matrix fraction, dual system, elementary */
/*     polynomial operations, Hessenberg form, minimal realization, */
/*     orthogonal transformation, polynomial matrix, state-space */
/*     representation, transfer matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 309 "TB03AD.f"
    /* Parameter adjustments */
#line 309 "TB03AD.f"
    a_dim1 = *lda;
#line 309 "TB03AD.f"
    a_offset = 1 + a_dim1;
#line 309 "TB03AD.f"
    a -= a_offset;
#line 309 "TB03AD.f"
    b_dim1 = *ldb;
#line 309 "TB03AD.f"
    b_offset = 1 + b_dim1;
#line 309 "TB03AD.f"
    b -= b_offset;
#line 309 "TB03AD.f"
    c_dim1 = *ldc;
#line 309 "TB03AD.f"
    c_offset = 1 + c_dim1;
#line 309 "TB03AD.f"
    c__ -= c_offset;
#line 309 "TB03AD.f"
    d_dim1 = *ldd;
#line 309 "TB03AD.f"
    d_offset = 1 + d_dim1;
#line 309 "TB03AD.f"
    d__ -= d_offset;
#line 309 "TB03AD.f"
    --index;
#line 309 "TB03AD.f"
    pcoeff_dim1 = *ldpco1;
#line 309 "TB03AD.f"
    pcoeff_dim2 = *ldpco2;
#line 309 "TB03AD.f"
    pcoeff_offset = 1 + pcoeff_dim1 * (1 + pcoeff_dim2);
#line 309 "TB03AD.f"
    pcoeff -= pcoeff_offset;
#line 309 "TB03AD.f"
    qcoeff_dim1 = *ldqco1;
#line 309 "TB03AD.f"
    qcoeff_dim2 = *ldqco2;
#line 309 "TB03AD.f"
    qcoeff_offset = 1 + qcoeff_dim1 * (1 + qcoeff_dim2);
#line 309 "TB03AD.f"
    qcoeff -= qcoeff_offset;
#line 309 "TB03AD.f"
    vcoeff_dim1 = *ldvco1;
#line 309 "TB03AD.f"
    vcoeff_dim2 = *ldvco2;
#line 309 "TB03AD.f"
    vcoeff_offset = 1 + vcoeff_dim1 * (1 + vcoeff_dim2);
#line 309 "TB03AD.f"
    vcoeff -= vcoeff_offset;
#line 309 "TB03AD.f"
    --iwork;
#line 309 "TB03AD.f"
    --dwork;
#line 309 "TB03AD.f"

#line 309 "TB03AD.f"
    /* Function Body */
#line 309 "TB03AD.f"
    *info = 0;
#line 310 "TB03AD.f"
    lleril = lsame_(leri, "L", (ftnlen)1, (ftnlen)1);
#line 311 "TB03AD.f"
    llerir = lsame_(leri, "R", (ftnlen)1, (ftnlen)1);
#line 312 "TB03AD.f"
    lequil = lsame_(equil, "S", (ftnlen)1, (ftnlen)1);
#line 313 "TB03AD.f"
    maxmp = max(*m,*p);
#line 314 "TB03AD.f"
    mplim = max(1,maxmp);
#line 315 "TB03AD.f"
    if (llerir) {

/*        Initialization for right matrix fraction. */

#line 319 "TB03AD.f"
	pwork = *m;
#line 320 "TB03AD.f"
	mwork = *p;
#line 321 "TB03AD.f"
    } else {

/*        Initialization for left matrix fraction. */

#line 325 "TB03AD.f"
	pwork = *p;
#line 326 "TB03AD.f"
	mwork = *m;
#line 327 "TB03AD.f"
    }

/*     Test the input scalar arguments. */

#line 331 "TB03AD.f"
    if (! lleril && ! llerir) {
#line 332 "TB03AD.f"
	*info = -1;
#line 333 "TB03AD.f"
    } else if (! lequil && ! lsame_(equil, "N", (ftnlen)1, (ftnlen)1)) {
#line 334 "TB03AD.f"
	*info = -2;
#line 335 "TB03AD.f"
    } else if (*n < 0) {
#line 336 "TB03AD.f"
	*info = -3;
#line 337 "TB03AD.f"
    } else if (*m < 0) {
#line 338 "TB03AD.f"
	*info = -4;
#line 339 "TB03AD.f"
    } else if (*p < 0) {
#line 340 "TB03AD.f"
	*info = -5;
#line 341 "TB03AD.f"
    } else if (*lda < max(1,*n)) {
#line 342 "TB03AD.f"
	*info = -7;
#line 343 "TB03AD.f"
    } else if (*ldb < max(1,*n)) {
#line 344 "TB03AD.f"
	*info = -9;
#line 345 "TB03AD.f"
    } else if (*ldc < mplim) {
#line 346 "TB03AD.f"
	*info = -11;
#line 347 "TB03AD.f"
    } else if (*ldd < mplim) {
#line 348 "TB03AD.f"
	*info = -13;
#line 349 "TB03AD.f"
    } else if (*ldpco1 < max(1,pwork)) {
#line 350 "TB03AD.f"
	*info = -17;
#line 351 "TB03AD.f"
    } else if (*ldpco2 < max(1,pwork)) {
#line 352 "TB03AD.f"
	*info = -18;
#line 353 "TB03AD.f"
    } else if (*ldqco1 < max(1,pwork) || llerir && *ldqco1 < mplim) {
#line 355 "TB03AD.f"
	*info = -20;
#line 356 "TB03AD.f"
    } else if (*ldqco2 < max(1,mwork) || llerir && *ldqco2 < mplim) {
#line 358 "TB03AD.f"
	*info = -21;
#line 359 "TB03AD.f"
    } else if (*ldvco1 < max(1,pwork)) {
#line 360 "TB03AD.f"
	*info = -23;
#line 361 "TB03AD.f"
    } else if (*ldvco2 < max(1,*n)) {
#line 362 "TB03AD.f"
	*info = -24;
#line 363 "TB03AD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
#line 363 "TB03AD.f"
	i__3 = *n, i__4 = maxmp * 3;
#line 363 "TB03AD.f"
	i__1 = 1, i__2 = *n + max(i__3,i__4), i__1 = max(i__1,i__2), i__2 = 
		pwork * (pwork + 2);
#line 363 "TB03AD.f"
	if (*ldwork < max(i__1,i__2)) {
#line 365 "TB03AD.f"
	    *info = -28;
#line 366 "TB03AD.f"
	}
#line 366 "TB03AD.f"
    }

#line 368 "TB03AD.f"
    if (*info != 0) {

/*        Error return. */

#line 372 "TB03AD.f"
	i__1 = -(*info);
#line 372 "TB03AD.f"
	xerbla_("TB03AD", &i__1, (ftnlen)6);
#line 373 "TB03AD.f"
	return 0;
#line 374 "TB03AD.f"
    }

/*     Quick return if possible. */

/* Computing MAX */
#line 378 "TB03AD.f"
    i__1 = max(*n,*m);
#line 378 "TB03AD.f"
    if (max(i__1,*p) == 0) {
#line 379 "TB03AD.f"
	*nr = 0;
#line 380 "TB03AD.f"
	dwork[1] = 1.;
#line 381 "TB03AD.f"
	return 0;
#line 382 "TB03AD.f"
    }

#line 384 "TB03AD.f"
    if (llerir) {

/*        For right matrix fraction, obtain dual system. */

#line 388 "TB03AD.f"
	ab07md_("D", n, m, p, &a[a_offset], lda, &b[b_offset], ldb, &c__[
		c_offset], ldc, &d__[d_offset], ldd, info, (ftnlen)1);
#line 390 "TB03AD.f"
    }

/*     Obtain minimal realization, in canonical form, for this system. */
/*     Part of the code in SLICOT routine TB01PD is included in-line */
/*     here. (TB01PD cannot be directly used.) */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

/*     If required, balance the triplet (A,B,C) (default MAXRED). */
/*     Workspace: need N. */

#line 405 "TB03AD.f"
    if (lequil) {
#line 406 "TB03AD.f"
	maxred = 0.;
#line 407 "TB03AD.f"
	tb01id_("A", n, &mwork, &pwork, &maxred, &a[a_offset], lda, &b[
		b_offset], ldb, &c__[c_offset], ldc, &dwork[1], info, (ftnlen)
		1);
#line 409 "TB03AD.f"
    }

#line 411 "TB03AD.f"
    iz = 1;
#line 412 "TB03AD.f"
    itau = 1;
#line 413 "TB03AD.f"
    jwork = itau + *n;

/*     Separate out controllable subsystem (of order NCONT): */
/*     A <-- Z'*A*Z,  B <-- Z'*B,  C <-- C*Z. */

/*     Workspace: need   N + MAX(N, 3*MWORK, PWORK). */
/*                prefer larger. */

#line 421 "TB03AD.f"
    i__1 = *ldwork - jwork + 1;
#line 421 "TB03AD.f"
    tb01ud_("No Z", n, &mwork, &pwork, &a[a_offset], lda, &b[b_offset], ldb, &
	    c__[c_offset], ldc, &ncont, &indblk, &iwork[1], &dwork[iz], &c__1,
	     &dwork[itau], tol, &iwork[*n + 1], &dwork[jwork], &i__1, info, (
	    ftnlen)4);

#line 425 "TB03AD.f"
    wrkopt = (integer) dwork[jwork] + jwork - 1;

/*     Separate out the observable subsystem (of order NR): */
/*     Form the dual of the subsystem of order NCONT (which is */
/*     controllable), leaving rest as it is. */

#line 431 "TB03AD.f"
    ab07md_("Z", &ncont, &mwork, &pwork, &a[a_offset], lda, &b[b_offset], ldb,
	     &c__[c_offset], ldc, &dwork[1], &c__1, info, (ftnlen)1);

/*     And separate out the controllable part of this dual subsystem. */

/*     Workspace: need   NCONT + MAX(NCONT, 3*PWORK, MWORK). */
/*                prefer larger. */

#line 439 "TB03AD.f"
    i__1 = *ldwork - jwork + 1;
#line 439 "TB03AD.f"
    tb01ud_("No Z", &ncont, &pwork, &mwork, &a[a_offset], lda, &b[b_offset], 
	    ldb, &c__[c_offset], ldc, nr, &indblk, &iwork[1], &dwork[iz], &
	    c__1, &dwork[itau], tol, &iwork[*n + 1], &dwork[jwork], &i__1, 
	    info, (ftnlen)4);

/* Computing MAX */
#line 443 "TB03AD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 443 "TB03AD.f"
    wrkopt = max(i__1,i__2);

/*     Retranspose, giving controllable and observable (i.e. minimal) */
/*     part of original system. */

#line 448 "TB03AD.f"
    ab07md_("Z", nr, &pwork, &mwork, &a[a_offset], lda, &b[b_offset], ldb, &
	    c__[c_offset], ldc, &dwork[1], &c__1, info, (ftnlen)1);

/*     Annihilate the trailing components of IWORK(1:N). */

#line 453 "TB03AD.f"
    i__1 = *n;
#line 453 "TB03AD.f"
    for (i__ = indblk + 1; i__ <= i__1; ++i__) {
#line 454 "TB03AD.f"
	iwork[i__] = 0;
#line 455 "TB03AD.f"
/* L10: */
#line 455 "TB03AD.f"
    }

/*     Initialize polynomial matrices P(s), Q(s) and V(s) to zero. */

#line 459 "TB03AD.f"
    i__1 = *n + 1;
#line 459 "TB03AD.f"
    for (k = 1; k <= i__1; ++k) {
#line 460 "TB03AD.f"
	dlaset_("Full", &pwork, &pwork, &c_b20, &c_b20, &pcoeff[(k * 
		pcoeff_dim2 + 1) * pcoeff_dim1 + 1], ldpco1, (ftnlen)4);
#line 462 "TB03AD.f"
	dlaset_("Full", &pwork, &mwork, &c_b20, &c_b20, &qcoeff[(k * 
		qcoeff_dim2 + 1) * qcoeff_dim1 + 1], ldqco1, (ftnlen)4);
#line 464 "TB03AD.f"
	dlaset_("Full", &pwork, nr, &c_b20, &c_b20, &vcoeff[(k * vcoeff_dim2 
		+ 1) * vcoeff_dim1 + 1], ldvco1, (ftnlen)4);
#line 466 "TB03AD.f"
/* L20: */
#line 466 "TB03AD.f"
    }

/*     Finish initializing V(s), and set up row degrees of P(s). */

#line 470 "TB03AD.f"
    inplus = indblk + 1;
#line 471 "TB03AD.f"
    istart = 1;
#line 472 "TB03AD.f"
    joff = *nr;

#line 474 "TB03AD.f"
    i__1 = indblk;
#line 474 "TB03AD.f"
    for (k = 1; k <= i__1; ++k) {
#line 475 "TB03AD.f"
	kwork = inplus - k;
#line 476 "TB03AD.f"
	kplus = kwork + 1;
#line 477 "TB03AD.f"
	istop = iwork[kwork];
#line 478 "TB03AD.f"
	joff -= istop;

#line 480 "TB03AD.f"
	i__2 = istop;
#line 480 "TB03AD.f"
	for (i__ = istart; i__ <= i__2; ++i__) {
#line 481 "TB03AD.f"
	    index[i__] = kwork;
#line 482 "TB03AD.f"
	    vcoeff[i__ + (joff + i__ + kplus * vcoeff_dim2) * vcoeff_dim1] = 
		    1.;
#line 483 "TB03AD.f"
/* L30: */
#line 483 "TB03AD.f"
	}

#line 485 "TB03AD.f"
	istart = istop + 1;
#line 486 "TB03AD.f"
/* L40: */
#line 486 "TB03AD.f"
    }

/*     ISTART = IWORK(1)+1 now: if .LE. PWORK, set up final rows of P(s). */

#line 490 "TB03AD.f"
    i__1 = pwork;
#line 490 "TB03AD.f"
    for (i__ = istart; i__ <= i__1; ++i__) {
#line 491 "TB03AD.f"
	index[i__] = 0;
#line 492 "TB03AD.f"
	pcoeff[i__ + (i__ + pcoeff_dim2) * pcoeff_dim1] = 1.;
#line 493 "TB03AD.f"
/* L50: */
#line 493 "TB03AD.f"
    }

/*     Triangularize the superdiagonal blocks of Amin. */

#line 497 "TB03AD.f"
    nrow = iwork[indblk];
#line 498 "TB03AD.f"
    ioff = *nr - nrow;
#line 499 "TB03AD.f"
    kmax = indblk - 1;
#line 500 "TB03AD.f"
    itau = 1;
#line 501 "TB03AD.f"
    ifirst = 0;
#line 502 "TB03AD.f"
    if (indblk > 2) {
#line 502 "TB03AD.f"
	ifirst = ioff - iwork[kmax];
#line 502 "TB03AD.f"
    }

/*     QR decomposition of each superdiagonal block of A in turn */
/*     (done in reverse order to preserve upper triangular blocks in A). */

#line 507 "TB03AD.f"
    i__1 = kmax;
#line 507 "TB03AD.f"
    for (k = 1; k <= i__1; ++k) {

/*        Calculate dimensions of new block & its position in A. */

#line 511 "TB03AD.f"
	kwork = indblk - k;
#line 512 "TB03AD.f"
	ncol = nrow;
#line 513 "TB03AD.f"
	nrow = iwork[kwork];
#line 514 "TB03AD.f"
	joff = ioff;
#line 515 "TB03AD.f"
	ioff -= nrow;
#line 516 "TB03AD.f"
	nreflc = min(nrow,ncol);
#line 517 "TB03AD.f"
	jwork = itau + nreflc;
#line 518 "TB03AD.f"
	if (kwork >= 2) {
#line 518 "TB03AD.f"
	    ifirst -= iwork[kwork - 1];
#line 518 "TB03AD.f"
	}

/*        Find QR decomposition of this (full rank) block: */
/*        block = QR.  No pivoting is needed. */

/*        Workspace: need   MIN(NROW,NCOL) + NCOL; */
/*                   prefer MIN(NROW,NCOL) + NCOL*NB. */

#line 526 "TB03AD.f"
	i__2 = *ldwork - jwork + 1;
#line 526 "TB03AD.f"
	dgeqrf_(&nrow, &ncol, &a[ioff + 1 + (joff + 1) * a_dim1], lda, &dwork[
		itau], &dwork[jwork], &i__2, info);

/* Computing MAX */
#line 529 "TB03AD.f"
	i__2 = wrkopt, i__3 = (integer) dwork[jwork] + jwork - 1;
#line 529 "TB03AD.f"
	wrkopt = max(i__2,i__3);

/*        Premultiply appropriate row block of A by Q'. */

/*        Workspace: need   MIN(NROW,NCOL) + JOFF; */
/*                   prefer MIN(NROW,NCOL) + JOFF*NB. */

#line 536 "TB03AD.f"
	i__2 = *ldwork - jwork + 1;
#line 536 "TB03AD.f"
	dormqr_("Left", "Transpose", &nrow, &joff, &nreflc, &a[ioff + 1 + (
		joff + 1) * a_dim1], lda, &dwork[itau], &a[ioff + 1 + a_dim1],
		 lda, &dwork[jwork], &i__2, info, (ftnlen)4, (ftnlen)9);

/* Computing MAX */
#line 540 "TB03AD.f"
	i__2 = wrkopt, i__3 = (integer) dwork[jwork] + jwork - 1;
#line 540 "TB03AD.f"
	wrkopt = max(i__2,i__3);

/*        Premultiply appropriate row block of B by Q' also. */

/*        Workspace: need   MIN(NROW,NCOL) + MWORK; */
/*                   prefer MIN(NROW,NCOL) + MWORK*NB. */

#line 547 "TB03AD.f"
	i__2 = *ldwork - jwork + 1;
#line 547 "TB03AD.f"
	dormqr_("Left", "Transpose", &nrow, &mwork, &nreflc, &a[ioff + 1 + (
		joff + 1) * a_dim1], lda, &dwork[itau], &b[ioff + 1 + b_dim1],
		 ldb, &dwork[jwork], &i__2, info, (ftnlen)4, (ftnlen)9);

/* Computing MAX */
#line 551 "TB03AD.f"
	i__2 = wrkopt, i__3 = (integer) dwork[jwork] + jwork - 1;
#line 551 "TB03AD.f"
	wrkopt = max(i__2,i__3);

/*        And postmultiply the non-zero part of appropriate column */
/*        block of A by Q. */

/*        Workspace: need   MIN(NROW,NCOL) + NR; */
/*                   prefer MIN(NROW,NCOL) + NR*NB. */

#line 559 "TB03AD.f"
	i__2 = *nr - ifirst;
#line 559 "TB03AD.f"
	i__3 = *ldwork - jwork + 1;
#line 559 "TB03AD.f"
	dormqr_("Right", "No Transpose", &i__2, &nrow, &nreflc, &a[ioff + 1 + 
		(joff + 1) * a_dim1], lda, &dwork[itau], &a[ifirst + 1 + (
		ioff + 1) * a_dim1], lda, &dwork[jwork], &i__3, info, (ftnlen)
		5, (ftnlen)12);

/* Computing MAX */
#line 564 "TB03AD.f"
	i__2 = wrkopt, i__3 = (integer) dwork[jwork] + jwork - 1;
#line 564 "TB03AD.f"
	wrkopt = max(i__2,i__3);

/*        Annihilate the lower triangular part of the block in A. */

#line 568 "TB03AD.f"
	if (k != kmax && nrow > 1) {
#line 568 "TB03AD.f"
	    i__2 = nrow - 1;
#line 568 "TB03AD.f"
	    dlaset_("Lower", &i__2, &ncol, &c_b20, &c_b20, &a[ioff + 2 + (
		    joff + 1) * a_dim1], lda, (ftnlen)5);
#line 568 "TB03AD.f"
	}

#line 572 "TB03AD.f"
/* L60: */
#line 572 "TB03AD.f"
    }

/*     Finally: postmultiply non-zero columns of C by Q (K = KMAX). */

/*     Workspace: need   MIN(NROW,NCOL) + PWORK; */
/*                prefer MIN(NROW,NCOL) + PWORK*NB. */

#line 579 "TB03AD.f"
    i__1 = *ldwork - jwork + 1;
#line 579 "TB03AD.f"
    dormqr_("Right", "No Transpose", &pwork, &nrow, &nreflc, &a[ioff + 1 + (
	    joff + 1) * a_dim1], lda, &dwork[itau], &c__[c_offset], ldc, &
	    dwork[jwork], &i__1, info, (ftnlen)5, (ftnlen)12);

/* Computing MAX */
#line 583 "TB03AD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 583 "TB03AD.f"
    wrkopt = max(i__1,i__2);

/*     Annihilate the lower triangular part of the block in A. */

#line 587 "TB03AD.f"
    if (nrow > 1) {
#line 587 "TB03AD.f"
	i__1 = nrow - 1;
#line 587 "TB03AD.f"
	dlaset_("Lower", &i__1, &ncol, &c_b20, &c_b20, &a[ioff + 2 + (joff + 
		1) * a_dim1], lda, (ftnlen)5);
#line 587 "TB03AD.f"
    }

/*     Calculate the (PWORK x NR) polynomial matrix V(s) ... */

#line 593 "TB03AD.f"
    tb03ay_(nr, &a[a_offset], lda, &indblk, &iwork[1], &vcoeff[vcoeff_offset],
	     ldvco1, ldvco2, &pcoeff[pcoeff_offset], ldpco1, ldpco2, info);

#line 596 "TB03AD.f"
    if (*info != 0) {
#line 597 "TB03AD.f"
	*info = 1;
#line 598 "TB03AD.f"
	return 0;
#line 599 "TB03AD.f"
    } else {

/*        And then use this matrix to calculate P(s): first store */
/*        C1 from C. */

#line 604 "TB03AD.f"
	ic = 1;
#line 605 "TB03AD.f"
	irankc = iwork[1];
#line 606 "TB03AD.f"
	ldwric = max(1,pwork);
#line 607 "TB03AD.f"
	dlacpy_("Full", &pwork, &irankc, &c__[c_offset], ldc, &dwork[ic], &
		ldwric, (ftnlen)4);

#line 609 "TB03AD.f"
	if (irankc < pwork) {

/*           rank(C) .LT. PWORK: obtain QR decomposition of C1, */
/*           giving R and Q. */

/*           Workspace: need   PWORK*IRANKC + 2*IRANKC; */
/*                      prefer PWORK*IRANKC +   IRANKC + IRANKC*NB. */

#line 617 "TB03AD.f"
	    itau = ic + ldwric * irankc;
#line 618 "TB03AD.f"
	    jwork = itau + irankc;

#line 620 "TB03AD.f"
	    i__1 = *ldwork - jwork + 1;
#line 620 "TB03AD.f"
	    dgeqrf_(&pwork, &irankc, &dwork[ic], &ldwric, &dwork[itau], &
		    dwork[jwork], &i__1, info);

/* Computing MAX */
#line 623 "TB03AD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 623 "TB03AD.f"
	    wrkopt = max(i__1,i__2);

/*           First IRANKC rows of Pbar(s) are given by Wbar(s) * inv(R). */
/*           Check for zero diagonal elements of R. */

#line 628 "TB03AD.f"
	    i__1 = irankc;
#line 628 "TB03AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 629 "TB03AD.f"
		if (dwork[ic + (i__ - 1) * ldwric + i__ - 1] == 0.) {

/*                 Error return. */

#line 633 "TB03AD.f"
		    *info = 2;
#line 634 "TB03AD.f"
		    return 0;
#line 635 "TB03AD.f"
		}
#line 636 "TB03AD.f"
/* L70: */
#line 636 "TB03AD.f"
	    }

#line 638 "TB03AD.f"
	    nrow = irankc;

#line 640 "TB03AD.f"
	    i__1 = inplus;
#line 640 "TB03AD.f"
	    for (k = 1; k <= i__1; ++k) {
#line 641 "TB03AD.f"
		dtrsm_("Right", "Upper", "No Transpose", "Non-unit", &nrow, &
			irankc, &c_b53, &dwork[ic], &ldwric, &pcoeff[(k * 
			pcoeff_dim2 + 1) * pcoeff_dim1 + 1], ldpco1, (ftnlen)
			5, (ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 644 "TB03AD.f"
		nrow = iwork[k];
#line 645 "TB03AD.f"
/* L80: */
#line 645 "TB03AD.f"
	    }

/*           P(s) itself is now given by Pbar(s) * Q'. */

#line 649 "TB03AD.f"
	    nrow = pwork;

#line 651 "TB03AD.f"
	    i__1 = inplus;
#line 651 "TB03AD.f"
	    for (k = 1; k <= i__1; ++k) {

/*              Workspace: need   PWORK*IRANKC + IRANKC + NROW; */
/*                         prefer PWORK*IRANKC + IRANKC + NROW*NB. */

#line 656 "TB03AD.f"
		i__2 = *ldwork - jwork + 1;
#line 656 "TB03AD.f"
		dormqr_("Right", "Transpose", &nrow, &pwork, &irankc, &dwork[
			ic], &ldwric, &dwork[itau], &pcoeff[(k * pcoeff_dim2 
			+ 1) * pcoeff_dim1 + 1], ldpco1, &dwork[jwork], &i__2,
			 info, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
#line 660 "TB03AD.f"
		i__2 = wrkopt, i__3 = (integer) dwork[jwork] + jwork - 1;
#line 660 "TB03AD.f"
		wrkopt = max(i__2,i__3);
#line 661 "TB03AD.f"
		nrow = iwork[k];
#line 662 "TB03AD.f"
/* L90: */
#line 662 "TB03AD.f"
	    }

#line 664 "TB03AD.f"
	} else {

/*           Special case rank(C) = PWORK, full: */
/*           no QR decomposition (P(s)=Wbar(s)*inv(C1)). */

#line 669 "TB03AD.f"
	    dgetrf_(&pwork, &pwork, &dwork[ic], &ldwric, &iwork[*n + 1], info)
		    ;

#line 672 "TB03AD.f"
	    if (*info != 0) {

/*              Error return. */

#line 676 "TB03AD.f"
		*info = 2;
#line 677 "TB03AD.f"
		return 0;
#line 678 "TB03AD.f"
	    } else {

#line 680 "TB03AD.f"
		nrow = irankc;

/*              Workspace: need   PWORK*IRANKC + N. */

#line 684 "TB03AD.f"
		i__1 = inplus;
#line 684 "TB03AD.f"
		for (k = 1; k <= i__1; ++k) {
#line 685 "TB03AD.f"
		    dtrsm_("Right", "Upper", "No Transpose", "Non-unit", &
			    nrow, &pwork, &c_b53, &dwork[ic], &ldwric, &
			    pcoeff[(k * pcoeff_dim2 + 1) * pcoeff_dim1 + 1], 
			    ldpco1, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)
			    8);
#line 688 "TB03AD.f"
		    dtrsm_("Right", "Lower", "No Transpose", "Unit", &nrow, &
			    pwork, &c_b53, &dwork[ic], &ldwric, &pcoeff[(k * 
			    pcoeff_dim2 + 1) * pcoeff_dim1 + 1], ldpco1, (
			    ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
#line 691 "TB03AD.f"
		    ma02gd_(&nrow, &pcoeff[(k * pcoeff_dim2 + 1) * 
			    pcoeff_dim1 + 1], ldpco1, &c__1, &pwork, &iwork[*
			    n + 1], &c_n1);
#line 693 "TB03AD.f"
		    nrow = iwork[k];
#line 694 "TB03AD.f"
/* L100: */
#line 694 "TB03AD.f"
		}
#line 695 "TB03AD.f"
	    }
#line 696 "TB03AD.f"
	}

/*        Finally, Q(s) = V(s) * B + P(s) * D can now be evaluated. */

#line 700 "TB03AD.f"
	nrow = pwork;

#line 702 "TB03AD.f"
	i__1 = inplus;
#line 702 "TB03AD.f"
	for (k = 1; k <= i__1; ++k) {
#line 703 "TB03AD.f"
	    dgemm_("No transpose", "No transpose", &nrow, &mwork, nr, &c_b53, 
		    &vcoeff[(k * vcoeff_dim2 + 1) * vcoeff_dim1 + 1], ldvco1, 
		    &b[b_offset], ldb, &c_b20, &qcoeff[(k * qcoeff_dim2 + 1) *
		     qcoeff_dim1 + 1], ldqco1, (ftnlen)12, (ftnlen)12);
#line 706 "TB03AD.f"
	    dgemm_("No transpose", "No transpose", &nrow, &mwork, &pwork, &
		    c_b53, &pcoeff[(k * pcoeff_dim2 + 1) * pcoeff_dim1 + 1], 
		    ldpco1, &d__[d_offset], ldd, &c_b53, &qcoeff[(k * 
		    qcoeff_dim2 + 1) * qcoeff_dim1 + 1], ldqco1, (ftnlen)12, (
		    ftnlen)12);
#line 709 "TB03AD.f"
	    nrow = iwork[k];
#line 710 "TB03AD.f"
/* L110: */
#line 710 "TB03AD.f"
	}

#line 712 "TB03AD.f"
    }

#line 714 "TB03AD.f"
    if (llerir) {

/*        For right matrix fraction, return to original (dual of dual) */
/*        system. */

#line 719 "TB03AD.f"
	ab07md_("Z", nr, &mwork, &pwork, &a[a_offset], lda, &b[b_offset], ldb,
		 &c__[c_offset], ldc, &dwork[1], &c__1, info, (ftnlen)1);

/*        Also, obtain the dual of the polynomial matrix representation. */

#line 724 "TB03AD.f"
	kpcoef = 0;

#line 726 "TB03AD.f"
	i__1 = pwork;
#line 726 "TB03AD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 727 "TB03AD.f"
	    i__2 = kpcoef, i__3 = index[i__];
#line 727 "TB03AD.f"
	    kpcoef = max(i__2,i__3);
#line 728 "TB03AD.f"
/* L120: */
#line 728 "TB03AD.f"
	}

#line 730 "TB03AD.f"
	++kpcoef;
#line 731 "TB03AD.f"
	tc01od_("L", &mwork, &pwork, &kpcoef, &pcoeff[pcoeff_offset], ldpco1, 
		ldpco2, &qcoeff[qcoeff_offset], ldqco1, ldqco2, info, (ftnlen)
		1);
#line 733 "TB03AD.f"
    } else {

/*        Reorder the rows and columns of the system, to get an upper */
/*        block Hessenberg matrix A of the minimal system. */

#line 738 "TB03AD.f"
	tb01yd_(nr, m, p, &a[a_offset], lda, &b[b_offset], ldb, &c__[c_offset]
		, ldc, info);
#line 739 "TB03AD.f"
    }

/*     Set optimal workspace dimension. */

#line 743 "TB03AD.f"
    dwork[1] = (doublereal) wrkopt;
#line 744 "TB03AD.f"
    return 0;
/* *** Last line of TB03AD *** */
} /* tb03ad_ */

