#line 1 "SB08MD.f"
/* SB08MD.f -- translated by f2c (version 20100827).
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

#line 1 "SB08MD.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static doublereal c_b30 = -1.;

/* Subroutine */ int sb08md_(char *acona, integer *da, doublereal *a, 
	doublereal *res, doublereal *e, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen acona_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal w, a0;
    static integer i0, nc;
    static doublereal si;
    static integer lq;
    static doublereal mu;
    static integer da1;
    static doublereal xda;
    static integer lay;
    static doublereal eps, muj;
    static integer binc, ldif, lphi;
    static logical conv;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal signi, signj;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sb08my_(integer *, doublereal *, 
	    doublereal *, doublereal *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal signi0, simin1, sqrta0;
    static integer lambda;
    extern doublereal dlamch_(char *, ftnlen);
    static logical lacona;
    extern integer idamax_(integer *, doublereal *, integer *);
    static logical stable;
    static integer lphend, layend;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal tolphi, sqrtmj, sqrtmu;


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

/*     To compute a real polynomial E(s) such that */

/*        (a)  E(-s) * E(s) = A(-s) * A(s) and */
/*        (b)  E(s) is stable - that is, all the zeros of E(s) have */
/*             non-positive real parts, */

/*     which corresponds to computing the spectral factorization of the */
/*     real polynomial A(s) arising from continuous optimality problems. */

/*     The input polynomial may be supplied either in the form */

/*        A(s) = a(0) + a(1) * s + ... + a(DA) * s**DA */

/*     or as */

/*        B(s) = A(-s) * A(s) */
/*             = b(0) + b(1) * s**2  + ... + b(DA) * s**(2*DA)        (1) */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     ACONA   CHARACTER*1 */
/*             Indicates whether the coefficients of A(s) or B(s) = */
/*             A(-s) * A(s) are to be supplied as follows: */
/*             = 'A':  The coefficients of A(s) are to be supplied; */
/*             = 'B':  The coefficients of B(s) are to be supplied. */

/*     Input/Output Parameters */

/*     DA      (input) INTEGER */
/*             The degree of the polynomials A(s) and E(s).  DA >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (DA+1) */
/*             On entry, this array must contain either the coefficients */
/*             of the polynomial A(s) in increasing powers of s if */
/*             ACONA = 'A', or the coefficients of the polynomial B(s) in */
/*             increasing powers of s**2 (see equation (1)) if ACONA = */
/*             'B'. */
/*             On exit, this array contains the coefficients of the */
/*             polynomial B(s) in increasing powers of s**2. */

/*     RES     (output) DOUBLE PRECISION */
/*             An estimate of the accuracy with which the coefficients of */
/*             the polynomial E(s) have been computed (see also METHOD */
/*             and NUMERICAL ASPECTS). */

/*     E       (output) DOUBLE PRECISION array, dimension (DA+1) */
/*             The coefficients of the spectral factor E(s) in increasing */
/*             powers of s. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= 5*DA+5. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if on entry, A(I) = 0.0, for I = 1,2,...,DA+1. */
/*             = 2:  if on entry, ACONA = 'B' but the supplied */
/*                   coefficients of the polynomial B(s) are not the */
/*                   coefficients of A(-s) * A(s) for some real A(s); */
/*                   in this case, RES and E are unassigned; */
/*             = 3:  if the iterative process (see METHOD) has failed to */
/*                   converge in 30 iterations; */
/*             = 4:  if the last computed iterate (see METHOD) is */
/*                   unstable. If ACONA = 'B', then the supplied */
/*                   coefficients of the polynomial B(s) may not be the */
/*                   coefficients of A(-s) * A(s) for some real A(s). */

/*     METHOD */
/*         _                                               _ */
/*     Let A(s) be the conjugate polynomial of A(s), i.e., A(s) = A(-s). */

/*     The method used by the routine is based on applying the */
/*     Newton-Raphson iteration to the function */
/*               _       _ */
/*        F(e) = A * A - e * e, */

/*     which leads to the iteration formulae (see [1]): */

/*        _(i)   (i)  _(i)   (i)     _      ) */
/*        q   * x   + x   * q    = 2 A * A  ) */
/*                                          )   for i = 0, 1, 2,... */
/*         (i+1)    (i)   (i)               ) */
/*        q     = (q   + x   )/2            ) */

/*                    (0)         DA */
/*     Starting from q   = (1 + s)   (which has no zeros in the closed */
/*                                                  (1)   (2)   (3) */
/*     right half-plane), the sequence of iterates q   , q   , q   ,... */
/*     converges to a solution of F(e) = 0 which has no zeros in the */
/*     open right half-plane. */

/*     The iterates satisfy the following conditions: */

/*              (i) */
/*        (a)  q   is a stable polynomial (no zeros in the closed right */
/*             half-plane) and */

/*              (i)        (i-1) */
/*        (b)  q   (1) <= q     (1). */

/*                                       (i-1)                       (i) */
/*     The iterative process stops with q     , (where i <= 30)  if q */
/*     violates either (a) or (b), or if the condition */
/*                       _(i) (i)  _ */
/*        (c)  RES  = ||(q   q   - A A)|| < tol, */

/*     is satisfied, where || . || denotes the largest coefficient of */
/*                     _(i) (i)  _ */
/*     the polynomial (q   q   - A A) and tol is an estimate of the */
/*                                                    _(i)  (i) */
/*     rounding error in the computed coefficients of q    q   . If there */
/*     is no convergence after 30 iterations then the routine returns */
/*     with the Error Indicator (INFO) set to 3, and the value of RES may */
/*     indicate whether or not the last computed iterate is close to the */
/*     solution. */

/*     If ACONA = 'B', then it is possible that the equation e(-s) * */
/*     e(s) = B(s) has no real solution, which will be the case if A(1) */
/*     < 0 or if ( -1)**DA * A(DA+1) < 0. */

/*     REFERENCES */

/*     [1] Vostry, Z. */
/*         New Algorithm for Polynomial Spectral Factorization with */
/*         Quadratic Convergence II. */
/*         Kybernetika, 12, pp. 248-259, 1976. */

/*     NUMERICAL ASPECTS */

/*     The conditioning of the problem depends upon the distance of the */
/*     zeros of A(s) from the imaginary axis and on their multiplicity. */
/*     For a well-conditioned problem the accuracy of the computed */
/*     coefficients of E(s) is of the order of RES. However, for problems */
/*     with zeros near the imaginary axis or with multiple zeros, the */
/*     value of RES may be an overestimate of the true accuracy. */

/*     FURTHER COMMENTS */

/*     In order for the problem e(-s) * e(s) = B(s) to have a real */
/*     solution e(s), it is necessary and sufficient that B(j*omega) */
/*     >= 0 for any purely imaginary argument j*omega (see [1]). */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB08AD by A.J. Geurts. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Factorization, Laplace transform, optimal control, optimal */
/*     filtering, polynomial operations, spectral factorization, zeros. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 216 "SB08MD.f"
    /* Parameter adjustments */
#line 216 "SB08MD.f"
    --dwork;
#line 216 "SB08MD.f"
    --e;
#line 216 "SB08MD.f"
    --a;
#line 216 "SB08MD.f"

#line 216 "SB08MD.f"
    /* Function Body */
#line 216 "SB08MD.f"
    *info = 0;
#line 217 "SB08MD.f"
    lacona = lsame_(acona, "A", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 221 "SB08MD.f"
    if (! lacona && ! lsame_(acona, "B", (ftnlen)1, (ftnlen)1)) {
#line 222 "SB08MD.f"
	*info = -1;
#line 223 "SB08MD.f"
    } else if (*da < 0) {
#line 224 "SB08MD.f"
	*info = -2;
#line 225 "SB08MD.f"
    } else if (*ldwork < *da * 5 + 5) {
#line 226 "SB08MD.f"
	*info = -7;
#line 227 "SB08MD.f"
    }

#line 229 "SB08MD.f"
    if (*info != 0) {

/*        Error return. */

#line 233 "SB08MD.f"
	i__1 = -(*info);
#line 233 "SB08MD.f"
	xerbla_("SB08MD", &i__1, (ftnlen)6);
#line 234 "SB08MD.f"
	return 0;
#line 235 "SB08MD.f"
    }

#line 237 "SB08MD.f"
    if (! lacona) {
#line 238 "SB08MD.f"
	i__1 = *da + 1;
#line 238 "SB08MD.f"
	dcopy_(&i__1, &a[1], &c__1, &e[1], &c__1);
#line 239 "SB08MD.f"
    } else {
#line 240 "SB08MD.f"
	w = 0.;
#line 241 "SB08MD.f"
	sb08my_(da, &a[1], &e[1], &w);
#line 242 "SB08MD.f"
    }

/*     Reduce E such that the first and the last element are non-zero. */

#line 246 "SB08MD.f"
    da1 = *da + 1;

/*     WHILE ( DA1 >= 1 and E(DA1) = 0 ) DO */
#line 249 "SB08MD.f"
L20:
#line 249 "SB08MD.f"
    if (da1 >= 1) {
#line 250 "SB08MD.f"
	if (e[da1] == 0.) {
#line 251 "SB08MD.f"
	    --da1;
#line 252 "SB08MD.f"
	    goto L20;
#line 253 "SB08MD.f"
	}
#line 254 "SB08MD.f"
    }
/*     END WHILE 20 */

#line 257 "SB08MD.f"
    --da1;
#line 258 "SB08MD.f"
    if (da1 < 0) {
#line 259 "SB08MD.f"
	*info = 1;
#line 260 "SB08MD.f"
	return 0;
#line 261 "SB08MD.f"
    }

#line 263 "SB08MD.f"
    i0 = 1;

/*     WHILE ( E(I0) = 0 ) DO */
#line 266 "SB08MD.f"
L40:
#line 266 "SB08MD.f"
    if (e[i0] == 0.) {
#line 267 "SB08MD.f"
	++i0;
#line 268 "SB08MD.f"
	goto L40;
#line 269 "SB08MD.f"
    }
/*     END WHILE 40 */

#line 272 "SB08MD.f"
    --i0;
#line 273 "SB08MD.f"
    if (i0 != 0) {
#line 274 "SB08MD.f"
	if (i0 % 2 == 0) {
#line 275 "SB08MD.f"
	    signi0 = 1.;
#line 276 "SB08MD.f"
	} else {
#line 277 "SB08MD.f"
	    signi0 = -1.;
#line 278 "SB08MD.f"
	}

#line 280 "SB08MD.f"
	i__1 = da1 - i0 + 1;
#line 280 "SB08MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 281 "SB08MD.f"
	    e[i__] = signi0 * e[i__ + i0];
#line 282 "SB08MD.f"
/* L60: */
#line 282 "SB08MD.f"
	}

#line 284 "SB08MD.f"
	da1 -= i0;
#line 285 "SB08MD.f"
    }
#line 286 "SB08MD.f"
    if (da1 % 2 == 0) {
#line 287 "SB08MD.f"
	signi = 1.;
#line 288 "SB08MD.f"
    } else {
#line 289 "SB08MD.f"
	signi = -1.;
#line 290 "SB08MD.f"
    }
#line 291 "SB08MD.f"
    nc = da1 + 1;
#line 292 "SB08MD.f"
    if (e[1] < 0. || e[nc] * signi < 0.) {
#line 293 "SB08MD.f"
	*info = 2;
#line 294 "SB08MD.f"
	return 0;
#line 295 "SB08MD.f"
    }

/*     Initialization. */

#line 299 "SB08MD.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 300 "SB08MD.f"
    si = 1. / dlamch_("Safe minimum", (ftnlen)12);
#line 301 "SB08MD.f"
    lq = 1;
#line 302 "SB08MD.f"
    lay = lq + nc;
#line 303 "SB08MD.f"
    lambda = lay + nc;
#line 304 "SB08MD.f"
    lphi = lambda + nc;
#line 305 "SB08MD.f"
    ldif = lphi + nc;

#line 307 "SB08MD.f"
    a0 = e[1];
#line 308 "SB08MD.f"
    binc = 1;

/*     Computation of the starting polynomial and scaling of the input */
/*     polynomial. */

#line 313 "SB08MD.f"
    d__2 = a0 / (d__1 = e[nc], abs(d__1));
#line 313 "SB08MD.f"
    d__3 = 1. / (doublereal) da1;
#line 313 "SB08MD.f"
    mu = pow_dd(&d__2, &d__3);
#line 314 "SB08MD.f"
    muj = 1.;

#line 316 "SB08MD.f"
    i__1 = nc;
#line 316 "SB08MD.f"
    for (j = 1; j <= i__1; ++j) {
#line 317 "SB08MD.f"
	w = e[j] * muj / a0;
#line 318 "SB08MD.f"
	a[j] = w;
#line 319 "SB08MD.f"
	e[j] = (doublereal) binc;
#line 320 "SB08MD.f"
	dwork[lq + j - 1] = (doublereal) binc;
#line 321 "SB08MD.f"
	muj *= mu;
#line 322 "SB08MD.f"
	binc = binc * (nc - j) / j;
#line 323 "SB08MD.f"
/* L80: */
#line 323 "SB08MD.f"
    }

#line 325 "SB08MD.f"
    conv = FALSE_;
#line 326 "SB08MD.f"
    stable = TRUE_;

/*     The contents of the arrays is, cf [1], */

/*     E : the last computed stable polynomial q   ; */
/*                                              i-1 */
/*     DWORK(LAY+1,...,LAY+DA1-1)  : a'(1), ..., a'(DA1-1), these values */
/*                                   are changed during the computation */
/*                                   into y; */
/*          (LAMBDA+1,...,LAMBDA+DA1-2) : lambda(1), ..., lambda(DA1-2), */
/*                                        the factors of the Routh */
/*                                        stability test, (lambda(i) is */
/*                                        P(i) in [1]); */
/*          (LPHI+1,...,LPHI+DA1-1) : phi(1), ..., phi(DA1-1), the values */
/*                                    phi(i,j), see [1], scheme (11); */
/*          (LDIF,...,LDIF+DA1) : the coeffs of q (-s) * q (s) - b(s). */
/*                                               i        i */
/*     DWORK(LQ,...,LQ+DA1) : the last computed polynomial q . */
/*                                                          i */
#line 345 "SB08MD.f"
    i__ = 0;

/*     WHILE ( I < 30 and CONV = FALSE and STABLE = TRUE ) DO */
#line 348 "SB08MD.f"
L100:
#line 348 "SB08MD.f"
    if (i__ < 30 && ! conv && stable) {
#line 349 "SB08MD.f"
	++i__;
#line 350 "SB08MD.f"
	dcopy_(&nc, &a[1], &c__1, &dwork[lay], &c__1);
#line 351 "SB08MD.f"
	dcopy_(&nc, &dwork[lq], &c__1, &dwork[lphi], &c__1);
#line 352 "SB08MD.f"
	m = da1 / 2;
#line 353 "SB08MD.f"
	layend = lay + da1;
#line 354 "SB08MD.f"
	lphend = lphi + da1;
#line 355 "SB08MD.f"
	xda = a[nc] / dwork[lq + da1];

#line 357 "SB08MD.f"
	i__1 = m;
#line 357 "SB08MD.f"
	for (k = 1; k <= i__1; ++k) {
#line 358 "SB08MD.f"
	    dwork[lay + k] -= dwork[lphi + (k << 1)];
#line 359 "SB08MD.f"
	    dwork[layend - k] -= dwork[lphend - (k << 1)] * xda;
#line 360 "SB08MD.f"
/* L120: */
#line 360 "SB08MD.f"
	}

/*        Computation of lambda(k) and y(k). */

#line 364 "SB08MD.f"
	k = 1;

/*        WHILE ( K <= DA1 - 2 and STABLE = TRUE ) DO */
#line 367 "SB08MD.f"
L140:
#line 367 "SB08MD.f"
	if (k <= da1 - 2 && stable) {
#line 368 "SB08MD.f"
	    if (dwork[lphi + k] <= 0.) {
#line 368 "SB08MD.f"
		stable = FALSE_;
#line 368 "SB08MD.f"
	    }
#line 369 "SB08MD.f"
	    if (stable) {
#line 370 "SB08MD.f"
		w = dwork[lphi + k - 1] / dwork[lphi + k];
#line 371 "SB08MD.f"
		dwork[lambda + k] = w;
#line 372 "SB08MD.f"
		i__1 = (da1 - k) / 2;
#line 372 "SB08MD.f"
		d__1 = -w;
#line 372 "SB08MD.f"
		daxpy_(&i__1, &d__1, &dwork[lphi + k + 2], &c__2, &dwork[lphi 
			+ k + 1], &c__2);
#line 374 "SB08MD.f"
		w = dwork[lay + k] / dwork[lphi + k];
#line 375 "SB08MD.f"
		dwork[lay + k] = w;
#line 376 "SB08MD.f"
		i__1 = (da1 - k) / 2;
#line 376 "SB08MD.f"
		d__1 = -w;
#line 376 "SB08MD.f"
		daxpy_(&i__1, &d__1, &dwork[lphi + k + 2], &c__2, &dwork[lay 
			+ k + 1], &c__1);
#line 378 "SB08MD.f"
		++k;
#line 379 "SB08MD.f"
	    }
#line 380 "SB08MD.f"
	    goto L140;
#line 381 "SB08MD.f"
	}
/*        END WHILE 140 */

#line 384 "SB08MD.f"
	if (dwork[lphi + da1 - 1] <= 0.) {
#line 385 "SB08MD.f"
	    stable = FALSE_;
#line 386 "SB08MD.f"
	} else {
#line 387 "SB08MD.f"
	    dwork[lay + da1 - 1] /= dwork[lphi + da1 - 1];
#line 388 "SB08MD.f"
	}

/*        STABLE = The polynomial q    is stable. */
/*                                 i-1 */
#line 392 "SB08MD.f"
	if (stable) {

/*           Computation of x  and q . */
/*                           i      i */

#line 397 "SB08MD.f"
	    for (k = da1 - 2; k >= 1; --k) {
#line 398 "SB08MD.f"
		w = dwork[lambda + k];
#line 399 "SB08MD.f"
		i__1 = (da1 - k) / 2;
#line 399 "SB08MD.f"
		d__1 = -w;
#line 399 "SB08MD.f"
		daxpy_(&i__1, &d__1, &dwork[lay + k + 1], &c__2, &dwork[lay + 
			k], &c__2);
#line 401 "SB08MD.f"
/* L160: */
#line 401 "SB08MD.f"
	    }

#line 403 "SB08MD.f"
	    dwork[lay + da1] = xda;

#line 405 "SB08MD.f"
	    dcopy_(&nc, &dwork[lq], &c__1, &e[1], &c__1);
#line 406 "SB08MD.f"
	    simin1 = si;
#line 407 "SB08MD.f"
	    si = dwork[lq];
#line 408 "SB08MD.f"
	    signj = -1.;

#line 410 "SB08MD.f"
	    i__1 = da1;
#line 410 "SB08MD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 411 "SB08MD.f"
		w = (dwork[lq + j] + signj * dwork[lay + j]) * .5;
#line 412 "SB08MD.f"
		dwork[lq + j] = w;
#line 413 "SB08MD.f"
		si += w;
#line 414 "SB08MD.f"
		signj = -signj;
#line 415 "SB08MD.f"
/* L180: */
#line 415 "SB08MD.f"
	    }

#line 417 "SB08MD.f"
	    tolphi = eps;
#line 418 "SB08MD.f"
	    sb08my_(&da1, &e[1], &dwork[ldif], &tolphi);
#line 419 "SB08MD.f"
	    daxpy_(&nc, &c_b30, &a[1], &c__1, &dwork[ldif], &c__1);
#line 420 "SB08MD.f"
	    *res = (d__1 = dwork[idamax_(&nc, &dwork[ldif], &c__1) + ldif - 1]
		    , abs(d__1));

/*           Convergency test. */

#line 424 "SB08MD.f"
	    if (si > simin1 || *res < tolphi) {
#line 425 "SB08MD.f"
		conv = TRUE_;
#line 426 "SB08MD.f"
	    }
#line 427 "SB08MD.f"
	    goto L100;
#line 428 "SB08MD.f"
	}
#line 429 "SB08MD.f"
    }
/*     END WHILE 100 */

/*     Backscaling. */

#line 434 "SB08MD.f"
    mu = 1. / mu;
#line 435 "SB08MD.f"
    sqrta0 = sqrt(a0);
#line 436 "SB08MD.f"
    sqrtmu = sqrt(mu);
#line 437 "SB08MD.f"
    muj = 1.;
#line 438 "SB08MD.f"
    sqrtmj = 1.;

#line 440 "SB08MD.f"
    i__1 = nc;
#line 440 "SB08MD.f"
    for (j = 1; j <= i__1; ++j) {
#line 441 "SB08MD.f"
	e[j] = e[j] * sqrta0 * sqrtmj;
#line 442 "SB08MD.f"
	a[j] = a[j] * a0 * muj;
#line 443 "SB08MD.f"
	muj *= mu;
#line 444 "SB08MD.f"
	sqrtmj *= sqrtmu;
#line 445 "SB08MD.f"
/* L200: */
#line 445 "SB08MD.f"
    }

#line 447 "SB08MD.f"
    if (i0 != 0) {

#line 449 "SB08MD.f"
	for (j = nc; j >= 1; --j) {
#line 450 "SB08MD.f"
	    e[i0 + j] = e[j];
#line 451 "SB08MD.f"
	    a[i0 + j] = signi0 * a[j];
#line 452 "SB08MD.f"
/* L220: */
#line 452 "SB08MD.f"
	}

#line 454 "SB08MD.f"
	i__1 = i0;
#line 454 "SB08MD.f"
	for (j = 1; j <= i__1; ++j) {
#line 455 "SB08MD.f"
	    e[j] = 0.;
#line 456 "SB08MD.f"
	    a[j] = 0.;
#line 457 "SB08MD.f"
/* L240: */
#line 457 "SB08MD.f"
	}

#line 459 "SB08MD.f"
    }

#line 461 "SB08MD.f"
    if (! conv) {
#line 462 "SB08MD.f"
	if (stable) {
#line 463 "SB08MD.f"
	    *info = 3;
#line 464 "SB08MD.f"
	} else {
#line 465 "SB08MD.f"
	    *info = 4;
#line 466 "SB08MD.f"
	}
#line 467 "SB08MD.f"
    }

#line 469 "SB08MD.f"
    return 0;
/* *** Last line of SB08MD *** */
} /* sb08md_ */

