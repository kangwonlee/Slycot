#line 1 "MB05OD.f"
/* MB05OD.f -- translated by f2c (version 20100827).
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

#line 1 "MB05OD.f"
/* Table of constant values */

static doublereal c_b10 = 0.;
static doublereal c_b11 = 1.;
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b25 = 2.;

/* Subroutine */ int mb05od_(char *balanc, integer *n, integer *ndiag, 
	doublereal *delta, doublereal *a, integer *lda, integer *mdig, 
	integer *idig, integer *iwork, doublereal *dwork, integer *ldwork, 
	integer *iwarn, integer *info, ftnlen balanc_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_lg10(doublereal *), sqrt(doublereal), exp(doublereal), log(
	    doublereal), pow_di(doublereal *, integer *);

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal p, s, u, bd;
    static integer ij, ik;
    static doublereal fn, gn, ss, tr, xn;
    static integer im1;
    static doublereal sd2, big, eps, var, tmp1;
    static integer ndec, base;
    static doublereal eabs, rerl, temp, rerr, size;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal sum2d;
    extern /* Subroutine */ int mb04md_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *);
    static integer ifail;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static logical lbals;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal avgev, small;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int mb05oy_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer ndagm1, ndagm2, ndecm1, jwora1, jwora2, jwora3;
    static doublereal ovrth2;
    static char actbal[1];
    static integer jworv1, jworv2;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dgetrf_(integer *, integer *, 
	    doublereal *, integer *, integer *, integer *), dlacpy_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen);
    static doublereal eavgev, factor;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static doublereal maxred;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal underf;
    extern /* Subroutine */ int dgetrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static doublereal emnorm, vareps;
    static integer mpower;
    static doublereal ovrthr;


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

/*     To compute exp(A*delta) where A is a real N-by-N matrix and delta */
/*     is a scalar value. The routine also returns the minimal number of */
/*     accurate digits in the 1-norm of exp(A*delta) and the number of */
/*     accurate digits in the 1-norm of exp(A*delta) at 95% confidence */
/*     level. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     BALANC  CHARACTER*1 */
/*             Specifies whether or not a balancing transformation (done */
/*             by SLICOT Library routine MB04MD) is required, as follows: */
/*             = 'N', do not use balancing; */
/*             = 'S', use balancing (scaling). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     NDIAG   (input) INTEGER */
/*             The specified order of the diagonal Pade approximant. */
/*             In the absence of further information NDIAG should */
/*             be set to 9.  NDIAG should not exceed 15.  NDIAG >= 1. */

/*     DELTA   (input) DOUBLE PRECISION */
/*             The scalar value delta of the problem. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On input, the leading N-by-N part of this array must */
/*             contain the matrix A of the problem. (This is not needed */
/*             if DELTA = 0.) */
/*             On exit, if INFO = 0, the leading N-by-N part of this */
/*             array contains the solution matrix exp(A*delta). */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     MDIG    (output) INTEGER */
/*             The minimal number of accurate digits in the 1-norm of */
/*             exp(A*delta). */

/*     IDIG    (output) INTEGER */
/*             The number of accurate digits in the 1-norm of */
/*             exp(A*delta) at 95% confidence level. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= N*(2*N+NDIAG+1)+NDIAG, if N >  1. */
/*             LDWORK >= 1,                     if N <= 1. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  if MDIG = 0 and IDIG > 0, warning for possible */
/*                   inaccuracy (the exponential has been computed); */
/*             = 2:  if MDIG = 0 and IDIG = 0, warning for severe */
/*                   inaccuracy (the exponential has been computed); */
/*             = 3:  if balancing has been requested, but it failed to */
/*                   reduce the matrix norm and was not actually used. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the norm of matrix A*delta (after a possible */
/*                   balancing) is too large to obtain an accurate */
/*                   result; */
/*             = 2:  if the coefficient matrix (the denominator of the */
/*                   Pade approximant) is exactly singular; try a */
/*                   different value of NDIAG; */
/*             = 3:  if the solution exponential would overflow, possibly */
/*                   due to a too large value DELTA; the calculations */
/*                   stopped prematurely. This error is not likely to */
/*                   appear. */

/*     METHOD */

/*     The exponential of the matrix A is evaluated from a diagonal Pade */
/*     approximant. This routine is a modification of the subroutine */
/*     PADE, described in reference [1]. The routine implements an */
/*     algorithm which exploits the identity */

/*         (exp[(2**-m)*A]) ** (2**m) = exp(A), */

/*     where m is an integer determined by the algorithm, to improve the */
/*     accuracy for matrices with large norms. */

/*     REFERENCES */

/*     [1] Ward, R.C. */
/*         Numerical computation of the matrix exponential with accuracy */
/*         estimate. */
/*         SIAM J. Numer. Anal., 14, pp. 600-610, 1977. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997. */
/*     Supersedes Release 2.0 routine MB05CD by T.W.C. Williams, Kingston */
/*     Polytechnic, March 1982. */

/*     REVISIONS */

/*     June 14, 1997, April 25, 2003, December 12, 2004. */

/*     KEYWORDS */

/*     Continuous-time system, matrix algebra, matrix exponential, */
/*     matrix operations, Pade approximation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 187 "MB05OD.f"
    /* Parameter adjustments */
#line 187 "MB05OD.f"
    a_dim1 = *lda;
#line 187 "MB05OD.f"
    a_offset = 1 + a_dim1;
#line 187 "MB05OD.f"
    a -= a_offset;
#line 187 "MB05OD.f"
    --iwork;
#line 187 "MB05OD.f"
    --dwork;
#line 187 "MB05OD.f"

#line 187 "MB05OD.f"
    /* Function Body */
#line 187 "MB05OD.f"
    *iwarn = 0;
#line 188 "MB05OD.f"
    *info = 0;
#line 189 "MB05OD.f"
    lbals = lsame_(balanc, "S", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 193 "MB05OD.f"
    if (! (lsame_(balanc, "N", (ftnlen)1, (ftnlen)1) || lbals)) {
#line 194 "MB05OD.f"
	*info = -1;
#line 195 "MB05OD.f"
    } else if (*n < 0) {
#line 196 "MB05OD.f"
	*info = -2;
#line 197 "MB05OD.f"
    } else if (*ndiag < 1) {
#line 198 "MB05OD.f"
	*info = -3;
#line 199 "MB05OD.f"
    } else if (*lda < max(1,*n)) {
#line 200 "MB05OD.f"
	*info = -6;
#line 201 "MB05OD.f"
    } else if (*ldwork < 1 || *ldwork < *n * ((*n << 1) + *ndiag + 1) + *
	    ndiag && *n > 1) {
#line 204 "MB05OD.f"
	*info = -11;
#line 205 "MB05OD.f"
    }

#line 207 "MB05OD.f"
    if (*info != 0) {

/*        Error return. */

#line 211 "MB05OD.f"
	i__1 = -(*info);
#line 211 "MB05OD.f"
	xerbla_("MB05OD", &i__1, (ftnlen)6);
#line 212 "MB05OD.f"
	return 0;
#line 213 "MB05OD.f"
    }

/*     Quick return if possible. */

#line 217 "MB05OD.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 218 "MB05OD.f"
    d__1 = 1. / eps;
#line 218 "MB05OD.f"
    ndec = (integer) (d_lg10(&d__1) + 1.);

#line 220 "MB05OD.f"
    if (*n == 0) {
#line 221 "MB05OD.f"
	*mdig = ndec;
#line 222 "MB05OD.f"
	*idig = ndec;
#line 223 "MB05OD.f"
	return 0;
#line 224 "MB05OD.f"
    }

/*     Set some machine parameters. */

#line 228 "MB05OD.f"
    base = (integer) dlamch_("Base", (ftnlen)4);
#line 229 "MB05OD.f"
    ndecm1 = ndec - 1;
#line 230 "MB05OD.f"
    underf = dlamch_("Underflow", (ftnlen)9);
#line 231 "MB05OD.f"
    ovrthr = dlamch_("Overflow", (ftnlen)8);
#line 232 "MB05OD.f"
    ovrth2 = sqrt(ovrthr);

#line 234 "MB05OD.f"
    if (*delta == 0.) {

/*        The DELTA = 0 case. */

#line 238 "MB05OD.f"
	dlaset_("Full", n, n, &c_b10, &c_b11, &a[a_offset], lda, (ftnlen)4);
#line 239 "MB05OD.f"
	*mdig = ndecm1;
#line 240 "MB05OD.f"
	*idig = ndecm1;
#line 241 "MB05OD.f"
	return 0;
#line 242 "MB05OD.f"
    }

#line 244 "MB05OD.f"
    if (*n == 1) {

/*        The 1-by-1 case. */

#line 248 "MB05OD.f"
	a[a_dim1 + 1] = exp(a[a_dim1 + 1] * *delta);
#line 249 "MB05OD.f"
	*mdig = ndecm1;
#line 250 "MB05OD.f"
	*idig = ndecm1;
#line 251 "MB05OD.f"
	return 0;
#line 252 "MB05OD.f"
    }

/*     Set pointers for the workspace. */

#line 256 "MB05OD.f"
    jwora1 = 1;
#line 257 "MB05OD.f"
    jwora2 = jwora1 + *n * *n;
#line 258 "MB05OD.f"
    jwora3 = jwora2 + *n * *ndiag;
#line 259 "MB05OD.f"
    jworv1 = jwora3 + *n * *n;
#line 260 "MB05OD.f"
    jworv2 = jworv1 + *n;

/*     Compute Pade coefficients in DWORK(JWORV2:JWORV2+NDIAG-1). */

#line 264 "MB05OD.f"
    dwork[jworv2] = .5;

#line 266 "MB05OD.f"
    i__1 = *ndiag;
#line 266 "MB05OD.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 267 "MB05OD.f"
	im1 = i__ - 1;
#line 268 "MB05OD.f"
	dwork[jworv2 + im1] = dwork[jworv2 + i__ - 2] * (doublereal) (*ndiag 
		- im1) / (doublereal) (i__ * ((*ndiag << 1) - im1));
#line 270 "MB05OD.f"
/* L20: */
#line 270 "MB05OD.f"
    }

/* Computing 2nd power */
#line 272 "MB05OD.f"
    d__1 = eps;
/* Computing 2nd power */
#line 272 "MB05OD.f"
    d__2 = (doublereal) base;
#line 272 "MB05OD.f"
    vareps = d__1 * d__1 * ((d__2 * d__2 - 1.) / (log((doublereal) base) * 
	    24.));
#line 274 "MB05OD.f"
    xn = (doublereal) (*n);
#line 275 "MB05OD.f"
    tr = 0.;

/*     Apply a translation with the mean of the eigenvalues of A*DELTA. */

#line 279 "MB05OD.f"
    i__1 = *n;
#line 279 "MB05OD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 280 "MB05OD.f"
	dscal_(n, delta, &a[i__ * a_dim1 + 1], &c__1);
#line 281 "MB05OD.f"
	tr += a[i__ + i__ * a_dim1];
#line 282 "MB05OD.f"
/* L40: */
#line 282 "MB05OD.f"
    }

#line 284 "MB05OD.f"
    avgev = tr / xn;
#line 285 "MB05OD.f"
    if (avgev > log(ovrthr) || avgev < log(underf)) {
#line 285 "MB05OD.f"
	avgev = 0.;
#line 285 "MB05OD.f"
    }
#line 287 "MB05OD.f"
    if (avgev != 0.) {
#line 288 "MB05OD.f"
	anorm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[jwora1], (
		ftnlen)6);

#line 290 "MB05OD.f"
	i__1 = *n;
#line 290 "MB05OD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 291 "MB05OD.f"
	    a[i__ + i__ * a_dim1] -= avgev;
#line 292 "MB05OD.f"
/* L60: */
#line 292 "MB05OD.f"
	}

#line 294 "MB05OD.f"
	temp = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[jwora1], (
		ftnlen)6);
#line 295 "MB05OD.f"
	if (temp > anorm * .5) {

#line 297 "MB05OD.f"
	    i__1 = *n;
#line 297 "MB05OD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 298 "MB05OD.f"
		a[i__ + i__ * a_dim1] += avgev;
#line 299 "MB05OD.f"
/* L80: */
#line 299 "MB05OD.f"
	    }

#line 301 "MB05OD.f"
	    avgev = 0.;
#line 302 "MB05OD.f"
	}
#line 303 "MB05OD.f"
    }
#line 304 "MB05OD.f"
    *(unsigned char *)actbal = *(unsigned char *)balanc;
#line 305 "MB05OD.f"
    if (lbals) {

/*        Balancing (scaling) has been requested.  First, save A. */

#line 309 "MB05OD.f"
	dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[jwora1], n, (ftnlen)4)
		;
#line 310 "MB05OD.f"
	maxred = 200.;
#line 311 "MB05OD.f"
	mb04md_(n, &maxred, &a[a_offset], lda, &dwork[jworv1], info);
#line 312 "MB05OD.f"
	if (maxred < 1.) {

/*           Recover the matrix and reset DWORK(JWORV1,...,JWORV1+N-1) */
/*           to 1, as no reduction of the norm occured (unlikely event). */

#line 317 "MB05OD.f"
	    dlacpy_("Full", n, n, &dwork[jwora1], n, &a[a_offset], lda, (
		    ftnlen)4);
#line 318 "MB05OD.f"
	    *(unsigned char *)actbal = 'N';
#line 319 "MB05OD.f"
	    dwork[jworv1] = 1.;
#line 320 "MB05OD.f"
	    i__1 = *n - 1;
#line 320 "MB05OD.f"
	    dcopy_(&i__1, &dwork[jworv1], &c__0, &dwork[jworv1 + 1], &c__1);
#line 321 "MB05OD.f"
	    *iwarn = 3;
#line 322 "MB05OD.f"
	}
#line 323 "MB05OD.f"
    }

/*     Scale the matrix by 2**(-M), where M is the minimum integer */
/*     so that the resulted matrix has the 1-norm less than 0.5. */

#line 328 "MB05OD.f"
    anorm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[jwora1], (
	    ftnlen)6);
#line 329 "MB05OD.f"
    m = 0;
#line 330 "MB05OD.f"
    if (anorm >= .5) {
#line 331 "MB05OD.f"
	mpower = (integer) (log(ovrthr) / log(2.));
#line 332 "MB05OD.f"
	m = (integer) (log(anorm) / log(2.)) + 1;
#line 333 "MB05OD.f"
	if (m > mpower) {

/*           Error return: The norm of A*DELTA is too large. */

#line 337 "MB05OD.f"
	    *info = 1;
#line 338 "MB05OD.f"
	    return 0;
#line 339 "MB05OD.f"
	}
#line 340 "MB05OD.f"
	factor = pow_di(&c_b25, &m);
#line 341 "MB05OD.f"
	if (m + 1 < mpower) {
#line 342 "MB05OD.f"
	    ++m;
#line 343 "MB05OD.f"
	    factor *= 2.;
#line 344 "MB05OD.f"
	}

#line 346 "MB05OD.f"
	i__1 = *n;
#line 346 "MB05OD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 347 "MB05OD.f"
	    d__1 = 1. / factor;
#line 347 "MB05OD.f"
	    dscal_(n, &d__1, &a[i__ * a_dim1 + 1], &c__1);
#line 348 "MB05OD.f"
/* L120: */
#line 348 "MB05OD.f"
	}

#line 350 "MB05OD.f"
    }
#line 351 "MB05OD.f"
    ndagm1 = *ndiag - 1;
#line 352 "MB05OD.f"
    ndagm2 = ndagm1 - 1;
#line 353 "MB05OD.f"
    ij = 0;

/*     Compute the factors of the diagonal Pade approximant. */
/*     The loop 200 takes the accuracy requirements into account: */
/*     Pade coefficients decrease with K, so the calculations should */
/*     be performed in backward order, one column at a time. */
/*     (A BLAS 3 implementation in forward order, using DGEMM, could */
/*     possibly be less accurate.) */

#line 362 "MB05OD.f"
    i__1 = *n;
#line 362 "MB05OD.f"
    for (j = 1; j <= i__1; ++j) {
#line 363 "MB05OD.f"
	dgemv_("No transpose", n, n, &c_b11, &a[a_offset], lda, &a[j * a_dim1 
		+ 1], &c__1, &c_b10, &dwork[jwora2], &c__1, (ftnlen)12);
#line 365 "MB05OD.f"
	ik = 0;

#line 367 "MB05OD.f"
	i__2 = ndagm2;
#line 367 "MB05OD.f"
	for (k = 1; k <= i__2; ++k) {
#line 368 "MB05OD.f"
	    dgemv_("No transpose", n, n, &c_b11, &a[a_offset], lda, &dwork[
		    jwora2 + ik], &c__1, &c_b10, &dwork[jwora2 + ik + *n], &
		    c__1, (ftnlen)12);
#line 371 "MB05OD.f"
	    ik += *n;
#line 372 "MB05OD.f"
/* L140: */
#line 372 "MB05OD.f"
	}

#line 374 "MB05OD.f"
	i__2 = *n;
#line 374 "MB05OD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 375 "MB05OD.f"
	    s = 0.;
#line 376 "MB05OD.f"
	    u = 0.;
#line 377 "MB05OD.f"
	    ik = ndagm2 * *n + i__ - 1;

#line 379 "MB05OD.f"
	    for (k = ndagm1; k >= 1; --k) {
#line 380 "MB05OD.f"
		p = dwork[jworv2 + k] * dwork[jwora2 + ik];
#line 381 "MB05OD.f"
		ik -= *n;
#line 382 "MB05OD.f"
		s += p;
#line 383 "MB05OD.f"
		if ((k + 1) % 2 == 0) {
#line 384 "MB05OD.f"
		    u += p;
#line 385 "MB05OD.f"
		} else {
#line 386 "MB05OD.f"
		    u -= p;
#line 387 "MB05OD.f"
		}
#line 388 "MB05OD.f"
/* L160: */
#line 388 "MB05OD.f"
	    }

#line 390 "MB05OD.f"
	    p = dwork[jworv2] * a[i__ + j * a_dim1];
#line 391 "MB05OD.f"
	    s += p;
#line 392 "MB05OD.f"
	    u -= p;
#line 393 "MB05OD.f"
	    if (i__ == j) {
#line 394 "MB05OD.f"
		s += 1.;
#line 395 "MB05OD.f"
		u += 1.;
#line 396 "MB05OD.f"
	    }
#line 397 "MB05OD.f"
	    dwork[jwora3 + ij] = s;
#line 398 "MB05OD.f"
	    dwork[jwora1 + ij] = u;
#line 399 "MB05OD.f"
	    ++ij;
#line 400 "MB05OD.f"
/* L180: */
#line 400 "MB05OD.f"
	}

#line 402 "MB05OD.f"
/* L200: */
#line 402 "MB05OD.f"
    }

/*     Compute the exponential of the scaled matrix, using diagonal Pade */
/*     approximants.  As, in theory [1], the denominator of the Pade */
/*     approximant should be very well conditioned, no condition estimate */
/*     is computed. */

#line 409 "MB05OD.f"
    dgetrf_(n, n, &dwork[jwora1], n, &iwork[1], &ifail);
#line 410 "MB05OD.f"
    if (ifail > 0) {

/*        Error return: The matrix is exactly singular. */

#line 414 "MB05OD.f"
	*info = 2;
#line 415 "MB05OD.f"
	return 0;
#line 416 "MB05OD.f"
    }

#line 418 "MB05OD.f"
    dlacpy_("Full", n, n, &dwork[jwora3], n, &a[a_offset], lda, (ftnlen)4);
#line 419 "MB05OD.f"
    dgetrs_("No transpose", n, n, &dwork[jwora1], n, &iwork[1], &a[a_offset], 
	    lda, &ifail, (ftnlen)12);

/*     Prepare for the calculation of the accuracy estimates. */
/*     Note that ANORM here is in the range [1, e]. */

#line 425 "MB05OD.f"
    anorm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[jwora1], (
	    ftnlen)6);
#line 426 "MB05OD.f"
    if (anorm >= 1.) {
#line 427 "MB05OD.f"
	eabs = (xn * 19. + 47.) * (eps * anorm);
#line 428 "MB05OD.f"
    } else {
#line 429 "MB05OD.f"
	eabs = (xn * 19. + 47.) * eps * anorm;
#line 430 "MB05OD.f"
    }
#line 431 "MB05OD.f"
    if (m != 0) {
#line 432 "MB05OD.f"
	var = xn * vareps;
#line 433 "MB05OD.f"
	fn = xn * 4. / ((xn + 2.) * (xn + 1.));
/* Computing 2nd power */
#line 434 "MB05OD.f"
	d__1 = xn + 2.;
/* Computing 2nd power */
#line 434 "MB05OD.f"
	d__2 = xn + 1.;
#line 434 "MB05OD.f"
	gn = ((xn * 2. + 10.) * xn - 4.) / (d__1 * d__1 * (d__2 * d__2));

/*        Square-up the computed exponential matrix M times, with caution */
/*        for avoiding overflows. */

#line 440 "MB05OD.f"
	i__1 = m;
#line 440 "MB05OD.f"
	for (k = 1; k <= i__1; ++k) {
#line 441 "MB05OD.f"
	    if (anorm > ovrth2) {

/*              The solution could overflow. */

#line 445 "MB05OD.f"
		d__1 = 1. / anorm;
#line 445 "MB05OD.f"
		dgemm_("No transpose", "No transpose", n, n, n, &d__1, &a[
			a_offset], lda, &a[a_offset], lda, &c_b10, &dwork[
			jwora1], n, (ftnlen)12, (ftnlen)12);
#line 448 "MB05OD.f"
		s = dlange_("1-norm", n, n, &dwork[jwora1], n, &dwork[jwora1],
			 (ftnlen)6);
#line 450 "MB05OD.f"
		if (anorm <= ovrthr / s) {
#line 451 "MB05OD.f"
		    dlascl_("General", n, n, &c_b11, &anorm, n, n, &dwork[
			    jwora1], n, info, (ftnlen)7);
#line 453 "MB05OD.f"
		    temp = ovrthr;
#line 454 "MB05OD.f"
		} else {

/*                 Error return: The solution would overflow. */
/*                 This will not happen on most machines, due to the */
/*                 selection of M. */

#line 460 "MB05OD.f"
		    *info = 3;
#line 461 "MB05OD.f"
		    return 0;
#line 462 "MB05OD.f"
		}
#line 463 "MB05OD.f"
	    } else {
#line 464 "MB05OD.f"
		dgemm_("No transpose", "No transpose", n, n, n, &c_b11, &a[
			a_offset], lda, &a[a_offset], lda, &c_b10, &dwork[
			jwora1], n, (ftnlen)12, (ftnlen)12);
/* Computing 2nd power */
#line 466 "MB05OD.f"
		d__1 = anorm;
#line 466 "MB05OD.f"
		temp = d__1 * d__1;
#line 467 "MB05OD.f"
	    }
#line 468 "MB05OD.f"
	    if (eabs < 1.) {
#line 469 "MB05OD.f"
		eabs = (anorm * 2. + eabs) * eabs + xn * (eps * temp);
#line 470 "MB05OD.f"
	    } else if (eabs < sqrt(1. - xn * eps + ovrthr / temp) * anorm - 
		    anorm) {
/* Computing 2nd power */
#line 472 "MB05OD.f"
		d__1 = eabs;
#line 472 "MB05OD.f"
		eabs = xn * (eps * temp) + anorm * eabs * 2. + d__1 * d__1;
#line 473 "MB05OD.f"
	    } else {
#line 474 "MB05OD.f"
		eabs = ovrthr;
#line 475 "MB05OD.f"
	    }

#line 477 "MB05OD.f"
	    tmp1 = fn * var + gn * (temp * vareps);
#line 478 "MB05OD.f"
	    if (tmp1 > ovrthr / temp) {
#line 479 "MB05OD.f"
		var = ovrthr;
#line 480 "MB05OD.f"
	    } else {
#line 481 "MB05OD.f"
		var = tmp1 * temp;
#line 482 "MB05OD.f"
	    }

#line 484 "MB05OD.f"
	    dlacpy_("Full", n, n, &dwork[jwora1], n, &a[a_offset], lda, (
		    ftnlen)4);
#line 485 "MB05OD.f"
	    anorm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[jwora1],
		     (ftnlen)6);
#line 486 "MB05OD.f"
/* L220: */
#line 486 "MB05OD.f"
	}

#line 488 "MB05OD.f"
    } else {
#line 489 "MB05OD.f"
	var = xn * 12. * vareps;
#line 490 "MB05OD.f"
    }

/*     Apply back transformations, if balancing was effectively used. */

#line 494 "MB05OD.f"
    mb05oy_(actbal, n, &c__1, n, &a[a_offset], lda, &dwork[jworv1], info, (
	    ftnlen)1);
#line 495 "MB05OD.f"
    eavgev = exp(avgev);
#line 496 "MB05OD.f"
    emnorm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[jwora1], (
	    ftnlen)6);

/*     Compute auxiliary quantities needed for the accuracy estimates. */

#line 500 "MB05OD.f"
    big = 1.;
#line 501 "MB05OD.f"
    small = 1.;
#line 502 "MB05OD.f"
    if (lbals) {

/*        Compute norms of the diagonal scaling matrix and its inverse. */

#line 506 "MB05OD.f"
	i__1 = *n;
#line 506 "MB05OD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 507 "MB05OD.f"
	    u = dwork[jworv1 + i__ - 1];
#line 508 "MB05OD.f"
	    if (big < u) {
#line 508 "MB05OD.f"
		big = u;
#line 508 "MB05OD.f"
	    }
#line 509 "MB05OD.f"
	    if (small > u) {
#line 509 "MB05OD.f"
		small = u;
#line 509 "MB05OD.f"
	    }
#line 510 "MB05OD.f"
/* L240: */
#line 510 "MB05OD.f"
	}

#line 512 "MB05OD.f"
	sum2d = dnrm2_(n, &dwork[jworv1], &c__1);
#line 513 "MB05OD.f"
    } else {
#line 514 "MB05OD.f"
	sum2d = sqrt(xn);
#line 515 "MB05OD.f"
    }

/*     Update the exponential for the initial translation, and update the */
/*     auxiliary quantities needed for the accuracy estimates. */

#line 520 "MB05OD.f"
    sd2 = sqrt(xn * 8. * vareps) * anorm;
#line 521 "MB05OD.f"
    bd = sqrt(var);
#line 522 "MB05OD.f"
    ss = max(bd,sd2);
#line 523 "MB05OD.f"
    bd = min(bd,sd2);
/* Computing 2nd power */
#line 524 "MB05OD.f"
    d__1 = bd / ss;
#line 524 "MB05OD.f"
    sd2 = ss * sqrt(d__1 * d__1 + 1.);
#line 525 "MB05OD.f"
    if (sd2 <= 1.) {
#line 526 "MB05OD.f"
	sd2 = 2. / xn * sum2d * sd2;
#line 527 "MB05OD.f"
    } else if (sum2d / xn < ovrthr / 2. / sd2) {
#line 528 "MB05OD.f"
	sd2 = 2. / xn * sum2d * sd2;
#line 529 "MB05OD.f"
    } else {
#line 530 "MB05OD.f"
	sd2 = ovrthr;
#line 531 "MB05OD.f"
    }
#line 532 "MB05OD.f"
    if (lbals) {
#line 533 "MB05OD.f"
	size = 0.;
#line 534 "MB05OD.f"
    } else {
#line 535 "MB05OD.f"
	if (sd2 < ovrthr - emnorm) {
#line 536 "MB05OD.f"
	    size = emnorm + sd2;
#line 537 "MB05OD.f"
	} else {
#line 538 "MB05OD.f"
	    size = ovrthr;
#line 539 "MB05OD.f"
	}
#line 540 "MB05OD.f"
    }

#line 542 "MB05OD.f"
    i__1 = *n;
#line 542 "MB05OD.f"
    for (j = 1; j <= i__1; ++j) {
#line 543 "MB05OD.f"
	ss = dasum_(n, &a[j * a_dim1 + 1], &c__1);
#line 544 "MB05OD.f"
	dscal_(n, &eavgev, &a[j * a_dim1 + 1], &c__1);
#line 545 "MB05OD.f"
	if (lbals) {
#line 546 "MB05OD.f"
	    bd = dwork[jworv1 + j - 1];
/* Computing MAX */
#line 547 "MB05OD.f"
	    d__1 = size, d__2 = ss + sd2 / bd;
#line 547 "MB05OD.f"
	    size = max(d__1,d__2);
#line 548 "MB05OD.f"
	}
#line 549 "MB05OD.f"
/* L260: */
#line 549 "MB05OD.f"
    }

/*     Set the accuracy estimates and warning errors, if any. */

#line 553 "MB05OD.f"
    rerr = d_lg10(&big) + d_lg10(&eabs) - d_lg10(&small) - d_lg10(&emnorm) - 
	    d_lg10(&eps);
#line 555 "MB05OD.f"
    if (size > emnorm) {
#line 556 "MB05OD.f"
	d__1 = (size / emnorm - 1.) / eps;
#line 556 "MB05OD.f"
	rerl = d_lg10(&d__1);
#line 557 "MB05OD.f"
    } else {
#line 558 "MB05OD.f"
	rerl = 0.;
#line 559 "MB05OD.f"
    }
/* Computing MIN */
#line 560 "MB05OD.f"
    i__1 = ndec - (integer) (rerr + .5);
#line 560 "MB05OD.f"
    *mdig = min(i__1,ndecm1);
/* Computing MIN */
#line 561 "MB05OD.f"
    i__1 = ndec - (integer) (rerl + .5);
#line 561 "MB05OD.f"
    *idig = min(i__1,ndecm1);

#line 563 "MB05OD.f"
    if (*mdig <= 0) {
#line 564 "MB05OD.f"
	*mdig = 0;
#line 565 "MB05OD.f"
	*iwarn = 1;
#line 566 "MB05OD.f"
    }
#line 567 "MB05OD.f"
    if (*idig <= 0) {
#line 568 "MB05OD.f"
	*idig = 0;
#line 569 "MB05OD.f"
	*iwarn = 2;
#line 570 "MB05OD.f"
    }

#line 572 "MB05OD.f"
    return 0;
/* *** Last line of MB05OD *** */
} /* mb05od_ */

