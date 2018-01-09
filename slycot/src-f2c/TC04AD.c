#line 1 "TC04AD.f"
/* TC04AD.f -- translated by f2c (version 20100827).
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

#line 1 "TC04AD.f"
/* Table of constant values */

static doublereal c_b12 = 0.;
static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b26 = -1.;
static doublereal c_b27 = 1.;

/* Subroutine */ int tc04ad_(char *leri, integer *m, integer *p, integer *
	index, doublereal *pcoeff, integer *ldpco1, integer *ldpco2, 
	doublereal *qcoeff, integer *ldqco1, integer *ldqco2, integer *n, 
	doublereal *rcond, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *d__, integer 
	*ldd, integer *iwork, doublereal *dwork, integer *ldwork, integer *
	info, ftnlen leri_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, pcoeff_dim1, pcoeff_dim2, pcoeff_offset, qcoeff_dim1, 
	    qcoeff_dim2, qcoeff_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, ia, ja, jc, jw, ldw;
    extern /* Subroutine */ int ab07md_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen);
    static integer ibias;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     tc01od_(char *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lleri;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer jwork, mwork, kstop, pwork;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgecon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen);
    static integer kpcoef;
    extern /* Subroutine */ int dgetrf_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), dlacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), dlaset_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), dgetri_(integer *,
	     doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *), xerbla_(char *, integer *, ftnlen);
    static integer maxind, mindex;
    extern /* Subroutine */ int dgetrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static doublereal dwnorm;
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

/*     To find a state-space representation (A,B,C,D) with the same */
/*     transfer matrix T(s) as that of a given left or right polynomial */
/*     matrix representation, i.e. */

/*        C*inv(sI-A)*B + D = T(s) = inv(P(s))*Q(s) = Q(s)*inv(P(s)). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     LERI    CHARACTER*1 */
/*             Indicates whether a left polynomial matrix representation */
/*             or a right polynomial matrix representation is input as */
/*             follows: */
/*             = 'L':  A left matrix fraction is input; */
/*             = 'R':  A right matrix fraction is input. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     INDEX   (input) INTEGER array, dimension (MAX(M,P)) */
/*             If LERI = 'L', INDEX(I), I = 1,2,...,P, must contain the */
/*             maximum degree of the polynomials in the I-th row of the */
/*             denominator matrix P(s) of the given left polynomial */
/*             matrix representation. */
/*             If LERI = 'R', INDEX(I), I = 1,2,...,M, must contain the */
/*             maximum degree of the polynomials in the I-th column of */
/*             the denominator matrix P(s) of the given right polynomial */
/*             matrix representation. */

/*     PCOEFF  (input) DOUBLE PRECISION array, dimension */
/*             (LDPCO1,LDPCO2,kpcoef), where kpcoef = MAX(INDEX(I)) + 1. */
/*             If LERI = 'L' then porm = P, otherwise porm = M. */
/*             The leading porm-by-porm-by-kpcoef part of this array must */
/*             contain the coefficients of the denominator matrix P(s). */
/*             PCOEFF(I,J,K) is the coefficient in s**(INDEX(iorj)-K+1) */
/*             of polynomial (I,J) of P(s), where K = 1,2,...,kpcoef; if */
/*             LERI = 'L' then iorj = I, otherwise iorj = J. */
/*             Thus for LERI = 'L', P(s) = */
/*             diag(s**INDEX(I))*(PCOEFF(.,.,1)+PCOEFF(.,.,2)/s+...). */
/*             If LERI = 'R', PCOEFF is modified by the routine but */
/*             restored on exit. */

/*     LDPCO1  INTEGER */
/*             The leading dimension of array PCOEFF. */
/*             LDPCO1 >= MAX(1,P) if LERI = 'L', */
/*             LDPCO1 >= MAX(1,M) if LERI = 'R'. */

/*     LDPCO2  INTEGER */
/*             The second dimension of array PCOEFF. */
/*             LDPCO2 >= MAX(1,P) if LERI = 'L', */
/*             LDPCO2 >= MAX(1,M) if LERI = 'R'. */

/*     QCOEFF  (input) DOUBLE PRECISION array, dimension */
/*             (LDQCO1,LDQCO2,kpcoef) */
/*             If LERI = 'L' then porp = M, otherwise porp = P. */
/*             The leading porm-by-porp-by-kpcoef part of this array must */
/*             contain the coefficients of the numerator matrix Q(s). */
/*             QCOEFF(I,J,K) is defined as for PCOEFF(I,J,K). */
/*             If LERI = 'R', QCOEFF is modified by the routine but */
/*             restored on exit. */

/*     LDQCO1  INTEGER */
/*             The leading dimension of array QCOEFF. */
/*             LDQCO1 >= MAX(1,P)   if LERI = 'L', */
/*             LDQCO1 >= MAX(1,M,P) if LERI = 'R'. */

/*     LDQCO2  INTEGER */
/*             The second dimension of array QCOEFF. */
/*             LDQCO2 >= MAX(1,M)   if LERI = 'L', */
/*             LDQCO2 >= MAX(1,M,P) if LERI = 'R'. */

/*     N       (output) INTEGER */
/*             The order of the resulting state-space representation. */
/*                          porm */
/*             That is, N = SUM INDEX(I). */
/*                          I=1 */

/*     RCOND   (output) DOUBLE PRECISION */
/*             The estimated reciprocal of the condition number of the */
/*             leading row (if LERI = 'L') or the leading column (if */
/*             LERI = 'R') coefficient matrix of P(s). */
/*             If RCOND is nearly zero, P(s) is nearly row or column */
/*             non-proper. */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array contains the state */
/*             dynamics matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,MAX(M,P)) */
/*             The leading N-by-M part of this array contains the */
/*             input/state matrix B; the remainder of the leading */
/*             N-by-MAX(M,P) part is used as internal workspace. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array contains the */
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

/*     Workspace */

/*     IWORK   INTEGER array, dimension (2*MAX(M,P)) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,MAX(M,P)*(MAX(M,P)+4)). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if P(s) is not row (if LERI = 'L') or column */
/*                   (if LERI = 'R') proper. Consequently, no state-space */
/*                   representation is calculated. */

/*     METHOD */

/*     The method for a left matrix fraction will be described here; */
/*     right matrix fractions are dealt with by obtaining the dual left */
/*     polynomial matrix representation and constructing an equivalent */
/*     state-space representation for this. The first step is to check */
/*     if the denominator matrix P(s) is row proper; if it is not then */
/*     the routine returns with the Error Indicator (INFO) set to 1. */
/*     Otherwise, Wolovich's Observable  Structure Theorem is used to */
/*     construct a state-space representation (A,B,C,D) in observable */
/*     companion form. The sizes of the blocks of matrix A and matrix C */
/*     here are precisely the row degrees of P(s), while their */
/*     'non-trivial' columns are given easily from its coefficients. */
/*     Similarly, the matrix D is obtained from the leading coefficients */
/*     of P(s) and of the numerator matrix Q(s), while matrix B is given */
/*     by the relation Sbar(s)B = Q(s) - P(s)D, where Sbar(s) is a */
/*     polynomial matrix whose (j,k)(th) element is given by */

/*                  j-u(k-1)-1 */
/*               ( s           , j = u(k-1)+1,u(k-1)+2,....,u(k) */
/*     Sbar    = ( */
/*        j,k    (           0 , otherwise */

/*             k */
/*     u(k) = SUM d , k = 1,2,...,M and d ,d ,...,d  are the */
/*            i=1  i                     1  2      M */
/*     controllability indices. For convenience in solving this, C' and B */
/*     are initially set up to contain the coefficients of P(s) and Q(s), */
/*     respectively, stored by rows. */

/*     REFERENCES */

/*     [1] Wolovich, W.A. */
/*         Linear Multivariate Systems, (Theorem 4.3.3). */
/*         Springer-Verlag, 1974. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TC01BD by T.W.C.Williams, Kingston */
/*     Polytechnic, United Kingdom, March 1982. */

/*     REVISIONS */

/*     February 22, 1998 (changed the name of TC01ND). */
/*     May 12, 1998. */

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

#line 258 "TC04AD.f"
    /* Parameter adjustments */
#line 258 "TC04AD.f"
    --index;
#line 258 "TC04AD.f"
    pcoeff_dim1 = *ldpco1;
#line 258 "TC04AD.f"
    pcoeff_dim2 = *ldpco2;
#line 258 "TC04AD.f"
    pcoeff_offset = 1 + pcoeff_dim1 * (1 + pcoeff_dim2);
#line 258 "TC04AD.f"
    pcoeff -= pcoeff_offset;
#line 258 "TC04AD.f"
    qcoeff_dim1 = *ldqco1;
#line 258 "TC04AD.f"
    qcoeff_dim2 = *ldqco2;
#line 258 "TC04AD.f"
    qcoeff_offset = 1 + qcoeff_dim1 * (1 + qcoeff_dim2);
#line 258 "TC04AD.f"
    qcoeff -= qcoeff_offset;
#line 258 "TC04AD.f"
    a_dim1 = *lda;
#line 258 "TC04AD.f"
    a_offset = 1 + a_dim1;
#line 258 "TC04AD.f"
    a -= a_offset;
#line 258 "TC04AD.f"
    b_dim1 = *ldb;
#line 258 "TC04AD.f"
    b_offset = 1 + b_dim1;
#line 258 "TC04AD.f"
    b -= b_offset;
#line 258 "TC04AD.f"
    c_dim1 = *ldc;
#line 258 "TC04AD.f"
    c_offset = 1 + c_dim1;
#line 258 "TC04AD.f"
    c__ -= c_offset;
#line 258 "TC04AD.f"
    d_dim1 = *ldd;
#line 258 "TC04AD.f"
    d_offset = 1 + d_dim1;
#line 258 "TC04AD.f"
    d__ -= d_offset;
#line 258 "TC04AD.f"
    --iwork;
#line 258 "TC04AD.f"
    --dwork;
#line 258 "TC04AD.f"

#line 258 "TC04AD.f"
    /* Function Body */
#line 258 "TC04AD.f"
    *info = 0;
#line 259 "TC04AD.f"
    lleri = lsame_(leri, "L", (ftnlen)1, (ftnlen)1);
#line 260 "TC04AD.f"
    mindex = max(*m,*p);

/*     Test the input scalar arguments. */

#line 264 "TC04AD.f"
    if (! lleri && ! lsame_(leri, "R", (ftnlen)1, (ftnlen)1)) {
#line 265 "TC04AD.f"
	*info = -1;
#line 266 "TC04AD.f"
    } else if (*m < 0) {
#line 267 "TC04AD.f"
	*info = -2;
#line 268 "TC04AD.f"
    } else if (*p < 0) {
#line 269 "TC04AD.f"
	*info = -3;
#line 270 "TC04AD.f"
    } else if (lleri && *ldpco1 < max(1,*p) || ! lleri && *ldpco1 < max(1,*m))
	     {
#line 272 "TC04AD.f"
	*info = -6;
#line 273 "TC04AD.f"
    } else if (lleri && *ldpco2 < max(1,*p) || ! lleri && *ldpco2 < max(1,*m))
	     {
#line 275 "TC04AD.f"
	*info = -7;
#line 276 "TC04AD.f"
    } else if (lleri && *ldqco1 < max(1,*p) || ! lleri && *ldqco1 < max(1,
	    mindex)) {
#line 278 "TC04AD.f"
	*info = -9;
#line 279 "TC04AD.f"
    } else if (lleri && *ldqco2 < max(1,*m) || ! lleri && *ldqco2 < max(1,
	    mindex)) {
#line 281 "TC04AD.f"
	*info = -10;
#line 282 "TC04AD.f"
    }

#line 284 "TC04AD.f"
    *n = 0;
#line 285 "TC04AD.f"
    if (*info == 0) {
#line 286 "TC04AD.f"
	if (lleri) {
#line 287 "TC04AD.f"
	    pwork = *p;
#line 288 "TC04AD.f"
	    mwork = *m;
#line 289 "TC04AD.f"
	} else {
#line 290 "TC04AD.f"
	    pwork = *m;
#line 291 "TC04AD.f"
	    mwork = *p;
#line 292 "TC04AD.f"
	}

#line 294 "TC04AD.f"
	maxind = 0;
#line 295 "TC04AD.f"
	i__1 = pwork;
#line 295 "TC04AD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 296 "TC04AD.f"
	    *n += index[i__];
#line 297 "TC04AD.f"
	    if (index[i__] > maxind) {
#line 297 "TC04AD.f"
		maxind = index[i__];
#line 297 "TC04AD.f"
	    }
#line 298 "TC04AD.f"
/* L10: */
#line 298 "TC04AD.f"
	}
#line 299 "TC04AD.f"
	kpcoef = maxind + 1;
#line 300 "TC04AD.f"
    }

#line 302 "TC04AD.f"
    if (*lda < max(1,*n)) {
#line 303 "TC04AD.f"
	*info = -14;
#line 304 "TC04AD.f"
    } else if (*ldb < max(1,*n)) {
#line 305 "TC04AD.f"
	*info = -16;
#line 306 "TC04AD.f"
    } else if (*ldc < max(1,mindex)) {
#line 307 "TC04AD.f"
	*info = -18;
#line 308 "TC04AD.f"
    } else if (*ldd < max(1,mindex)) {
#line 309 "TC04AD.f"
	*info = -20;
#line 310 "TC04AD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 310 "TC04AD.f"
	i__1 = 1, i__2 = mindex * (mindex + 4);
#line 310 "TC04AD.f"
	if (*ldwork < max(i__1,i__2)) {
#line 311 "TC04AD.f"
	    *info = -23;
#line 312 "TC04AD.f"
	}
#line 312 "TC04AD.f"
    }

#line 314 "TC04AD.f"
    if (*info != 0) {

/*        Error return. */

#line 318 "TC04AD.f"
	i__1 = -(*info);
#line 318 "TC04AD.f"
	xerbla_("TC04AD", &i__1, (ftnlen)6);
#line 319 "TC04AD.f"
	return 0;
#line 320 "TC04AD.f"
    }

/*     Quick return if possible. */

#line 324 "TC04AD.f"
    if (*m == 0 || *p == 0) {
#line 325 "TC04AD.f"
	*n = 0;
#line 326 "TC04AD.f"
	*rcond = 1.;
#line 327 "TC04AD.f"
	dwork[1] = 1.;
#line 328 "TC04AD.f"
	return 0;
#line 329 "TC04AD.f"
    }

#line 331 "TC04AD.f"
    if (! lleri) {

/*        Initialization for right matrix fraction: obtain the dual */
/*        system. */

#line 336 "TC04AD.f"
	tc01od_("R", m, p, &kpcoef, &pcoeff[pcoeff_offset], ldpco1, ldpco2, &
		qcoeff[qcoeff_offset], ldqco1, ldqco2, info, (ftnlen)1);
#line 338 "TC04AD.f"
    }

/*     Store leading row coefficient matrix of P(s). */

#line 342 "TC04AD.f"
    ldw = max(1,pwork);
#line 343 "TC04AD.f"
    dlacpy_("Full", &pwork, &pwork, &pcoeff[pcoeff_offset], ldpco1, &dwork[1],
	     &ldw, (ftnlen)4);

/*     Check if P(s) is row proper: if not, exit. */

#line 347 "TC04AD.f"
    dwnorm = dlange_("1-norm", &pwork, &pwork, &dwork[1], &ldw, &dwork[1], (
	    ftnlen)6);

#line 349 "TC04AD.f"
    dgetrf_(&pwork, &pwork, &dwork[1], &ldw, &iwork[1], info);

/*     Workspace: need  PWORK*(PWORK + 4). */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 359 "TC04AD.f"
    jwork = ldw * pwork + 1;

#line 361 "TC04AD.f"
    dgecon_("1-norm", &pwork, &dwork[1], &ldw, &dwnorm, rcond, &dwork[jwork], 
	    &iwork[pwork + 1], info, (ftnlen)6);

/* Computing MAX */
#line 364 "TC04AD.f"
    i__1 = 1, i__2 = pwork * (pwork + 4);
#line 364 "TC04AD.f"
    wrkopt = max(i__1,i__2);

#line 366 "TC04AD.f"
    if (*rcond <= dlamch_("Epsilon", (ftnlen)7)) {

/*        Error return: P(s) is not row proper. */

#line 370 "TC04AD.f"
	*info = 1;
#line 371 "TC04AD.f"
	return 0;
#line 372 "TC04AD.f"
    } else {

/*        Calculate the order of equivalent state-space representation, */
/*        and initialize A. */

#line 377 "TC04AD.f"
	dlaset_("Full", n, n, &c_b12, &c_b12, &a[a_offset], lda, (ftnlen)4);

#line 379 "TC04AD.f"
	dwork[jwork] = 1.;
#line 380 "TC04AD.f"
	if (*n > 1) {
#line 380 "TC04AD.f"
	    i__1 = *n - 1;
#line 380 "TC04AD.f"
	    i__2 = *lda + 1;
#line 380 "TC04AD.f"
	    dcopy_(&i__1, &dwork[jwork], &c__0, &a[a_dim1 + 2], &i__2);
#line 380 "TC04AD.f"
	}

/*        Find the PWORK ordered 'non-trivial' columns row by row, */
/*        in PWORK row blocks, the I-th having INDEX(I) rows. */

#line 385 "TC04AD.f"
	ibias = 2;

#line 387 "TC04AD.f"
	i__1 = pwork;
#line 387 "TC04AD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 388 "TC04AD.f"
	    kstop = index[i__] + 1;
#line 389 "TC04AD.f"
	    if (kstop != 1) {
#line 390 "TC04AD.f"
		ibias += index[i__];

/*              These rows given from the lower coefficients of row I */
/*              of P(s). */

#line 395 "TC04AD.f"
		i__2 = kstop;
#line 395 "TC04AD.f"
		for (k = 2; k <= i__2; ++k) {
#line 396 "TC04AD.f"
		    ia = ibias - k;

#line 398 "TC04AD.f"
		    i__3 = pwork;
#line 398 "TC04AD.f"
		    for (j = 1; j <= i__3; ++j) {
#line 399 "TC04AD.f"
			dwork[jwork + j - 1] = -pcoeff[i__ + (j + k * 
				pcoeff_dim2) * pcoeff_dim1];
#line 400 "TC04AD.f"
/* L20: */
#line 400 "TC04AD.f"
		    }

#line 402 "TC04AD.f"
		    dgetrs_("Transpose", &pwork, &c__1, &dwork[1], &ldw, &
			    iwork[1], &dwork[jwork], &ldw, info, (ftnlen)9);

#line 405 "TC04AD.f"
		    ja = 0;

#line 407 "TC04AD.f"
		    i__3 = pwork;
#line 407 "TC04AD.f"
		    for (j = 1; j <= i__3; ++j) {
#line 408 "TC04AD.f"
			if (index[j] != 0) {
#line 409 "TC04AD.f"
			    ja += index[j];
#line 410 "TC04AD.f"
			    a[ia + ja * a_dim1] = dwork[jwork + j - 1];
#line 411 "TC04AD.f"
			}
#line 412 "TC04AD.f"
/* L30: */
#line 412 "TC04AD.f"
		    }

/*                 Also, set up B and C (temporarily) for use when */
/*                 finding B. */

#line 417 "TC04AD.f"
		    dcopy_(&mwork, &qcoeff[i__ + (k * qcoeff_dim2 + 1) * 
			    qcoeff_dim1], ldqco1, &b[ia + b_dim1], ldb);
#line 419 "TC04AD.f"
		    dcopy_(&pwork, &pcoeff[i__ + (k * pcoeff_dim2 + 1) * 
			    pcoeff_dim1], ldpco1, &c__[ia * c_dim1 + 1], &
			    c__1);
#line 420 "TC04AD.f"
/* L40: */
#line 420 "TC04AD.f"
		}

#line 422 "TC04AD.f"
	    }
#line 423 "TC04AD.f"
/* L50: */
#line 423 "TC04AD.f"
	}

/*        Calculate D from the leading coefficients of P and Q. */

#line 427 "TC04AD.f"
	dlacpy_("Full", &pwork, &mwork, &qcoeff[qcoeff_offset], ldqco1, &d__[
		d_offset], ldd, (ftnlen)4);

#line 429 "TC04AD.f"
	dgetrs_("No transpose", &pwork, &mwork, &dwork[1], &ldw, &iwork[1], &
		d__[d_offset], ldd, info, (ftnlen)12);

/*        For B and C as set up above, desired B = B - (C' * D). */

#line 434 "TC04AD.f"
	dgemm_("Transpose", "No transpose", n, &mwork, &pwork, &c_b26, &c__[
		c_offset], ldc, &d__[d_offset], ldd, &c_b27, &b[b_offset], 
		ldb, (ftnlen)9, (ftnlen)12);

/*        Finally, calculate C: zero, apart from ... */

#line 439 "TC04AD.f"
	dlaset_("Full", &pwork, n, &c_b12, &c_b12, &c__[c_offset], ldc, (
		ftnlen)4);

/*        PWORK ordered 'non-trivial' columns, equal to those */
/*        of inv(DWORK). */

/*        Workspace: need   PWORK*(PWORK + 1); */
/*                   prefer PWORK*PWORK + PWORK*NB. */

#line 447 "TC04AD.f"
	i__1 = *ldwork - jwork + 1;
#line 447 "TC04AD.f"
	dgetri_(&pwork, &dwork[1], &ldw, &iwork[1], &dwork[jwork], &i__1, 
		info);

/* Computing MAX */
#line 450 "TC04AD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 450 "TC04AD.f"
	wrkopt = max(i__1,i__2);
#line 451 "TC04AD.f"
	jc = 0;
#line 452 "TC04AD.f"
	jw = 1;

#line 454 "TC04AD.f"
	i__1 = pwork;
#line 454 "TC04AD.f"
	for (j = 1; j <= i__1; ++j) {
#line 455 "TC04AD.f"
	    if (index[j] != 0) {
#line 456 "TC04AD.f"
		jc += index[j];
#line 457 "TC04AD.f"
		dcopy_(&pwork, &dwork[jw], &c__1, &c__[jc * c_dim1 + 1], &
			c__1);
#line 458 "TC04AD.f"
	    }
#line 459 "TC04AD.f"
	    jw += ldw;
#line 460 "TC04AD.f"
/* L60: */
#line 460 "TC04AD.f"
	}

#line 462 "TC04AD.f"
    }

/*     For right matrix fraction, return to original (dual of dual) */
/*     system. */

#line 467 "TC04AD.f"
    if (! lleri) {
#line 468 "TC04AD.f"
	tc01od_("L", &mwork, &pwork, &kpcoef, &pcoeff[pcoeff_offset], ldpco1, 
		ldpco2, &qcoeff[qcoeff_offset], ldqco1, ldqco2, info, (ftnlen)
		1);

/*        Also, obtain dual of state-space representation. */

#line 473 "TC04AD.f"
	ab07md_("D", n, &mwork, &pwork, &a[a_offset], lda, &b[b_offset], ldb, 
		&c__[c_offset], ldc, &d__[d_offset], ldd, info, (ftnlen)1);
#line 475 "TC04AD.f"
    }

/*     Set optimal workspace dimension. */

#line 479 "TC04AD.f"
    dwork[1] = (doublereal) wrkopt;

#line 481 "TC04AD.f"
    return 0;
/* *** Last line of TC04AD *** */
} /* tc04ad_ */

