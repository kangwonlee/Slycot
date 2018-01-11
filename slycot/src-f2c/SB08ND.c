#line 1 "SB08ND.f"
/* SB08ND.f -- translated by f2c (version 20100827).
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

#line 1 "SB08ND.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b11 = 2.;
static integer c_n1 = -1;
static doublereal c_b24 = -1.;

/* Subroutine */ int sb08nd_(char *acona, integer *da, doublereal *a, 
	doublereal *res, doublereal *e, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen acona_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s, w, a0;
    static integer nc, lq;
    static doublereal sa0;
    static integer nck, lro;
    static doublereal res0;
    static integer leta;
    static logical conv;
    static doublereal tolq;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), sb08ny_(integer *, doublereal *, 
	    doublereal *, doublereal *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer lambda;
    static logical lacona;
    static integer lalpha;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical hurwtz;


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

/*     To compute a real polynomial E(z) such that */

/*        (a)  E(1/z) * E(z) = A(1/z) * A(z) and */
/*        (b)  E(z) is stable - that is, E(z) has no zeros with modulus */
/*             greater than 1, */

/*     which corresponds to computing the spectral factorization of the */
/*     real polynomial A(z) arising from discrete optimality problems. */

/*     The input polynomial may be supplied either in the form */

/*     A(z) = a(0) + a(1) * z + ... + a(DA) * z**DA */

/*     or as */

/*     B(z) = A(1/z) * A(z) */
/*          = b(0) + b(1) * (z + 1/z) + ... + b(DA) * (z**DA + 1/z**DA) */
/*                                                                    (1) */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     ACONA   CHARACTER*1 */
/*             Indicates whether the coefficients of A(z) or B(z) = */
/*             A(1/z) * A(z) are to be supplied as follows: */
/*             = 'A':  The coefficients of A(z) are to be supplied; */
/*             = 'B':  The coefficients of B(z) are to be supplied. */

/*     Input/Output Parameters */

/*     DA      (input) INTEGER */
/*             The degree of the polynomials A(z) and E(z).  DA >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (DA+1) */
/*             On entry, if ACONA = 'A', this array must contain the */
/*             coefficients of the polynomial A(z) in increasing powers */
/*             of z, and if ACONA = 'B', this array must contain the */
/*             coefficients b ,b ,...,b   of the polynomial B(z) in */
/*                           0  1      DA */
/*             equation (1). That is, A(i) = b    for i = 1,2,...,DA+1. */
/*                                            i-1 */
/*             On exit, this array contains the coefficients of the */
/*             polynomial B(z) in eqation (1). Specifically, A(i) */
/*             contains b   ,  for i = 1,2,...DA+1. */
/*                       i-1 */

/*     RES     (output) DOUBLE PRECISION */
/*             An estimate of the accuracy with which the coefficients of */
/*             the polynomial E(z) have been computed (see also METHOD */
/*             and NUMERICAL ASPECTS). */

/*     E       (output) DOUBLE PRECISION array, dimension (DA+1) */
/*             The coefficients of the spectral factor E(z) in increasing */
/*             powers of z. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= 5*DA+5. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 2:  if on entry, ACONA = 'B' but the supplied */
/*                   coefficients of the polynomial B(z) are not the */
/*                   coefficients of A(1/z) * A(z) for some real A(z); */
/*                   in this case, RES and E are unassigned; */
/*             = 3:  if the iterative process (see METHOD) has failed to */
/*                   converge in 30 iterations; */
/*             = 4:  if the last computed iterate (see METHOD) is */
/*                   unstable. If ACONA = 'B', then the supplied */
/*                   coefficients of the polynomial B(z) may not be the */
/*                   coefficients of A(1/z) * A(z) for some real A(z). */

/*     METHOD */
/*         _                                               _ */
/*     Let A(z) be the conjugate polynomial of A(z), i.e., A(z) = A(1/z). */

/*     The method used by the routine is based on applying the */
/*     Newton-Raphson iteration to the function */
/*               _       _ */
/*        F(e) = A * A - e * e, */

/*     which leads to the iteration formulae (see [1] and [2]) */

/*        _(i)   (i)  _(i)   (i)     _      ) */
/*        q   * x   + x   * q    = 2 A * A  ) */
/*                                          )   for i = 0, 1, 2,... */
/*         (i+1)    (i)   (i)               ) */
/*        q     = (q   + x   )/2            ) */

/*     The iteration starts from */

/*         (0)                                        DA */
/*        q   (z) = (b(0) + b(1) * z + ... + b(DA) * z  ) / SQRT( b(0)) */

/*     which is a Hurwitz polynomial that has no zeros in the closed unit */
/*                                            (i) */
/*     circle (see [2], Theorem 3). Then lim q   = e, the convergence is */
/*     uniform and e is a Hurwitz polynomial. */

/*     The iterates satisfy the following conditions: */
/*              (i) */
/*        (a)  q    has no zeros in the closed unit circle, */
/*              (i)     (i-1) */
/*        (b)  q    <= q     and */
/*              0       0 */
/*              DA   (i) 2    DA     2 */
/*        (c)  SUM (q   )  - SUM (A )  >= 0. */
/*             k=0   k       k=0   k */
/*                                     (i) */
/*     The iterative process stops if q    violates (a), (b) or (c), */
/*     or if the condition */
/*                       _(i) (i)  _ */
/*        (d)  RES  = ||(q   q   - A A)|| < tol, */

/*     is satisfied, where || . || denotes the largest coefficient of */
/*                     _(i) (i)  _ */
/*     the polynomial (q   q   - A A) and tol is an estimate of the */
/*                                                    _(i)  (i) */
/*     rounding error in the computed coefficients of q    q   . If */
/*                                            (i-1) */
/*     condition (a) or (b) is violated then q      is taken otherwise */
/*      (i) */
/*     q    is used. Thus the computed reciprocal polynomial E(z) = z**DA */
/*     * q(1/z) is stable. If there is no convergence after 30 iterations */
/*     then the routine returns with the Error Indicator (INFO) set to 3, */
/*     and the value of RES may indicate whether or not the last computed */
/*     iterate is close to the solution. */
/*                                               (0) */
/*     If ACONA = 'B', then it is possible that q    is not a Hurwitz */
/*     polynomial, in which case the equation e(1/z) * e(z) = B(z) has no */
/*     real solution (see [2], Theorem 3). */

/*     REFERENCES */

/*     [1] Kucera, V. */
/*         Discrete Linear Control, The polynomial Approach. */
/*         John Wiley & Sons, Chichester, 1979. */

/*     [2] Vostry, Z. */
/*         New Algorithm for Polynomial Spectral Factorization with */
/*         Quadratic Convergence I. */
/*         Kybernetika, 11, pp. 415-422, 1975. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB08BD by F. Delebecque and */
/*     A.J. Geurts. */

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

#line 219 "SB08ND.f"
    /* Parameter adjustments */
#line 219 "SB08ND.f"
    --dwork;
#line 219 "SB08ND.f"
    --e;
#line 219 "SB08ND.f"
    --a;
#line 219 "SB08ND.f"

#line 219 "SB08ND.f"
    /* Function Body */
#line 219 "SB08ND.f"
    *info = 0;
#line 220 "SB08ND.f"
    lacona = lsame_(acona, "A", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 224 "SB08ND.f"
    if (! lacona && ! lsame_(acona, "B", (ftnlen)1, (ftnlen)1)) {
#line 225 "SB08ND.f"
	*info = -1;
#line 226 "SB08ND.f"
    } else if (*da < 0) {
#line 227 "SB08ND.f"
	*info = -2;
#line 228 "SB08ND.f"
    } else if (*ldwork < *da * 5 + 5) {
#line 229 "SB08ND.f"
	*info = -7;
#line 230 "SB08ND.f"
    }

#line 232 "SB08ND.f"
    if (*info != 0) {

/*        Error return. */

#line 236 "SB08ND.f"
	i__1 = -(*info);
#line 236 "SB08ND.f"
	xerbla_("SB08ND", &i__1, (ftnlen)6);
#line 237 "SB08ND.f"
	return 0;
#line 238 "SB08ND.f"
    }

#line 240 "SB08ND.f"
    nc = *da + 1;
#line 241 "SB08ND.f"
    if (! lacona) {
#line 242 "SB08ND.f"
	if (a[1] <= 0.) {
#line 243 "SB08ND.f"
	    *info = 2;
#line 244 "SB08ND.f"
	    return 0;
#line 245 "SB08ND.f"
	}
#line 246 "SB08ND.f"
	dcopy_(&nc, &a[1], &c__1, &e[1], &c__1);
#line 247 "SB08ND.f"
    } else {
#line 248 "SB08ND.f"
	sb08ny_(da, &a[1], &e[1], &w);
#line 249 "SB08ND.f"
    }

/*     Initialization. */

#line 253 "SB08ND.f"
    lalpha = 1;
#line 254 "SB08ND.f"
    lro = lalpha + nc;
#line 255 "SB08ND.f"
    leta = lro + nc;
#line 256 "SB08ND.f"
    lambda = leta + nc;
#line 257 "SB08ND.f"
    lq = lambda + nc;

#line 259 "SB08ND.f"
    a0 = e[1];
#line 260 "SB08ND.f"
    sa0 = sqrt(a0);
#line 261 "SB08ND.f"
    s = 0.;

#line 263 "SB08ND.f"
    i__1 = nc;
#line 263 "SB08ND.f"
    for (j = 1; j <= i__1; ++j) {
#line 264 "SB08ND.f"
	w = e[j];
#line 265 "SB08ND.f"
	a[j] = w;
#line 266 "SB08ND.f"
	w /= sa0;
#line 267 "SB08ND.f"
	e[j] = w;
#line 268 "SB08ND.f"
	dwork[lq - 1 + j] = w;
/* Computing 2nd power */
#line 269 "SB08ND.f"
	d__1 = w;
#line 269 "SB08ND.f"
	s += d__1 * d__1;
#line 270 "SB08ND.f"
/* L20: */
#line 270 "SB08ND.f"
    }

#line 272 "SB08ND.f"
    res0 = s - a0;

/*     The contents of the arrays is, cf [1], Section 7.6, */

/*     E : the last computed Hurwitz polynomial q   ; */
/*                                               i-1 */
/*     DWORK(LALPHA,..,LALPHA+DA-K)  : alpha(k,0),...alpha(k,n-k); */
/*          (LRO,...,LRO+DA-K)       : alpha(k,n-k),...,alpha(k); */
/*          (LETA,...,LETA+DA)       : eta(0),...,eta(n); */
/*          (LAMBDA,...,LAMBDA+DA-1) : lambda(0),...,lambda(n-1) */

/*     DWORK(LQ,...,LQ+DA) : the last computed polynomial q . */
/*                                                         i */
#line 285 "SB08ND.f"
    i__ = 0;
#line 286 "SB08ND.f"
    conv = FALSE_;
#line 287 "SB08ND.f"
    hurwtz = TRUE_;

/*     WHILE ( I < 30 and CONV = FALSE and HURWTZ = TRUE ) DO */
#line 290 "SB08ND.f"
L40:
#line 290 "SB08ND.f"
    if (i__ < 30 && ! conv && hurwtz) {
#line 291 "SB08ND.f"
	++i__;
#line 292 "SB08ND.f"
	dcopy_(&nc, &a[1], &c__1, &dwork[leta], &c__1);
#line 293 "SB08ND.f"
	dscal_(&nc, &c_b11, &dwork[leta], &c__1);
#line 294 "SB08ND.f"
	dcopy_(&nc, &dwork[lq], &c__1, &dwork[lalpha], &c__1);

/*        Computation of lambda(k) and eta(k). */

#line 298 "SB08ND.f"
	k = 1;

/*        WHILE ( K <= DA and HURWTZ = TRUE ) DO */
#line 301 "SB08ND.f"
L60:
#line 301 "SB08ND.f"
	if (k <= *da && hurwtz) {
#line 302 "SB08ND.f"
	    nck = nc - k;
#line 303 "SB08ND.f"
	    i__1 = nck + 1;
#line 303 "SB08ND.f"
	    dcopy_(&i__1, &dwork[lalpha], &c_n1, &dwork[lro], &c__1);
#line 304 "SB08ND.f"
	    w = dwork[lalpha + nck] / dwork[lro + nck];
#line 305 "SB08ND.f"
	    if (abs(w) >= 1.) {
#line 305 "SB08ND.f"
		hurwtz = FALSE_;
#line 305 "SB08ND.f"
	    }
#line 306 "SB08ND.f"
	    if (hurwtz) {
#line 307 "SB08ND.f"
		dwork[lambda + k - 1] = w;
#line 308 "SB08ND.f"
		d__1 = -w;
#line 308 "SB08ND.f"
		daxpy_(&nck, &d__1, &dwork[lro], &c__1, &dwork[lalpha], &c__1)
			;
#line 309 "SB08ND.f"
		w = dwork[leta + nck] / dwork[lalpha];
#line 310 "SB08ND.f"
		dwork[leta + nck] = w;
#line 311 "SB08ND.f"
		i__1 = nck - 1;
#line 311 "SB08ND.f"
		d__1 = -w;
#line 311 "SB08ND.f"
		daxpy_(&i__1, &d__1, &dwork[lalpha + 1], &c_n1, &dwork[leta + 
			1], &c__1);
#line 313 "SB08ND.f"
		++k;
#line 314 "SB08ND.f"
	    }
#line 315 "SB08ND.f"
	    goto L60;
#line 316 "SB08ND.f"
	}
/*        END WHILE 60 */

/*        HURWTZ = The polynomial q    is a Hurwitz polynomial. */
/*                                 i-1 */
#line 321 "SB08ND.f"
	if (hurwtz) {
#line 322 "SB08ND.f"
	    dcopy_(&nc, &dwork[lq], &c__1, &e[1], &c__1);

/*           Accuracy test. */

#line 326 "SB08ND.f"
	    sb08ny_(da, &e[1], &dwork[lq], &tolq);
#line 327 "SB08ND.f"
	    daxpy_(&nc, &c_b24, &a[1], &c__1, &dwork[lq], &c__1);
#line 328 "SB08ND.f"
	    *res = (d__1 = dwork[idamax_(&nc, &dwork[lq], &c__1) + lq - 1], 
		    abs(d__1));
#line 329 "SB08ND.f"
	    conv = *res < tolq || res0 < 0.;

#line 331 "SB08ND.f"
	    if (! conv) {
#line 332 "SB08ND.f"
		dwork[leta] = dwork[leta] * .5 / dwork[lalpha];

/*              Computation of x  and q . */
/*                              i      i */
/*              DWORK(LETA,...,LETA+DA)   : eta(k,0),...,eta(k,n) */
/*                   (LRO,...,LRO+DA-K+1) : eta(k,n-k+1),...,eta(k,0) */

#line 339 "SB08ND.f"
		for (k = *da; k >= 1; --k) {
#line 340 "SB08ND.f"
		    nck = nc - k + 1;
#line 341 "SB08ND.f"
		    dcopy_(&nck, &dwork[leta], &c_n1, &dwork[lro], &c__1);
#line 342 "SB08ND.f"
		    w = dwork[lambda + k - 1];
#line 343 "SB08ND.f"
		    d__1 = -w;
#line 343 "SB08ND.f"
		    daxpy_(&nck, &d__1, &dwork[lro], &c__1, &dwork[leta], &
			    c__1);
#line 344 "SB08ND.f"
/* L80: */
#line 344 "SB08ND.f"
		}

#line 346 "SB08ND.f"
		s = 0.;

#line 348 "SB08ND.f"
		i__1 = *da;
#line 348 "SB08ND.f"
		for (j = 0; j <= i__1; ++j) {
#line 349 "SB08ND.f"
		    w = (dwork[leta + j] + e[j + 1]) * .5;
#line 350 "SB08ND.f"
		    dwork[lq + j] = w;
/* Computing 2nd power */
#line 351 "SB08ND.f"
		    d__1 = w;
#line 351 "SB08ND.f"
		    s += d__1 * d__1;
#line 352 "SB08ND.f"
/* L100: */
#line 352 "SB08ND.f"
		}

#line 354 "SB08ND.f"
		res0 = s - a0;

/*              Test on the monotonicity of q . */
/*                                           0 */
#line 358 "SB08ND.f"
		conv = dwork[lq] > e[1];
#line 359 "SB08ND.f"
		goto L40;
#line 360 "SB08ND.f"
	    }
#line 361 "SB08ND.f"
	}
#line 362 "SB08ND.f"
    }
/*     END WHILE 40 */

/*     Reverse the order of the coefficients in the array E. */

#line 367 "SB08ND.f"
    dswap_(&nc, &e[1], &c__1, &dwork[1], &c_n1);
#line 368 "SB08ND.f"
    dswap_(&nc, &dwork[1], &c__1, &e[1], &c__1);

#line 370 "SB08ND.f"
    if (! conv) {
#line 371 "SB08ND.f"
	if (hurwtz) {
#line 372 "SB08ND.f"
	    *info = 3;
#line 373 "SB08ND.f"
	} else if (i__ == 1) {
#line 374 "SB08ND.f"
	    *info = 2;
#line 375 "SB08ND.f"
	} else {
#line 376 "SB08ND.f"
	    *info = 4;
#line 377 "SB08ND.f"
	}
#line 378 "SB08ND.f"
    }

#line 380 "SB08ND.f"
    return 0;
/* *** Last line of SB08ND *** */
} /* sb08nd_ */

