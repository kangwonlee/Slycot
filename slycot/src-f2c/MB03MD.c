#line 1 "MB03MD.f"
/* MB03MD.f -- translated by f2c (version 20100827).
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

#line 1 "MB03MD.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mb03md_(integer *n, integer *l, doublereal *theta, 
	doublereal *q, doublereal *e, doublereal *q2, doublereal *e2, 
	doublereal *pivmin, doublereal *tol, doublereal *reltol, integer *
	iwarn, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    static doublereal h__;
    static integer i__;
    static doublereal y, z__, th;
    static integer num, numz;
    extern integer mb03nd_(integer *, doublereal *, doublereal *, doublereal *
	    , doublereal *, integer *);
    extern doublereal mb03my_(integer *, doublereal *, integer *), dlamch_(
	    char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

/*     To compute an upper bound THETA using a bisection method such that */
/*     the bidiagonal matrix */

/*              |q(1) e(1)  0    ...   0   | */
/*              | 0   q(2) e(2)        .   | */
/*          J = | .                    .   | */
/*              | .                  e(N-1)| */
/*              | 0   ...        ...  q(N) | */

/*     has precisely L singular values less than or equal to THETA plus */
/*     a given tolerance TOL. */

/*     This routine is mainly intended to be called only by other SLICOT */
/*     routines. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the bidiagonal matrix J.  N >= 0. */

/*     L       (input/output) INTEGER */
/*             On entry, L must contain the number of singular values */
/*             of J which must be less than or equal to the upper bound */
/*             computed by the routine.  0 <= L <= N. */
/*             On exit, L may be increased if the L-th smallest singular */
/*             value of J has multiplicity greater than 1. In this case, */
/*             L is increased by the number of singular values of J which */
/*             are larger than its L-th smallest one and approach the */
/*             L-th smallest singular value of J within a distance less */
/*             than TOL. */
/*             If L has been increased, then the routine returns with */
/*             IWARN set to 1. */

/*     THETA   (input/output) DOUBLE PRECISION */
/*             On entry, THETA must contain an initial estimate for the */
/*             upper bound to be computed. If THETA < 0.0 on entry, then */
/*             one of the following default values is used. */
/*             If L = 0, THETA is set to 0.0 irrespective of the input */
/*             value of THETA; if L = 1, then THETA is taken as */
/*             MIN(ABS(Q(i))), for i = 1,2,...,N; otherwise, THETA is */
/*             taken as ABS(Q(N-L+1)). */
/*             On exit, THETA contains the computed upper bound such that */
/*             the bidiagonal matrix J has precisely L singular values */
/*             less than or equal to THETA + TOL. */

/*     Q       (input) DOUBLE PRECISION array, dimension (N) */
/*             This array must contain the diagonal elements q(1), */
/*             q(2),...,q(N) of the bidiagonal matrix J. That is, */
/*             Q(i) = J(i,i) for i = 1,2,...,N. */

/*     E       (input) DOUBLE PRECISION array, dimension (N-1) */
/*             This array must contain the superdiagonal elements */
/*             e(1),e(2),...,e(N-1) of the bidiagonal matrix J. That is, */
/*             E(k) = J(k,k+1) for k = 1,2,...,N-1. */

/*     Q2      (input) DOUBLE PRECISION array, dimension (N) */
/*             This array must contain the squares of the diagonal */
/*             elements q(1),q(2),...,q(N) of the bidiagonal matrix J. */
/*             That is, Q2(i) = J(i,i)**2 for i = 1,2,...,N. */

/*     E2      (input) DOUBLE PRECISION array, dimension (N-1) */
/*             This array must contain the squares of the superdiagonal */
/*             elements e(1),e(2),...,e(N-1) of the bidiagonal matrix J. */
/*             That is, E2(k) = J(k,k+1)**2 for k = 1,2,...,N-1. */

/*     PIVMIN  (input) DOUBLE PRECISION */
/*             The minimum absolute value of a "pivot" in the Sturm */
/*             sequence loop. */
/*             PIVMIN >= max( max( |q(i)|, |e(k)| )**2*sf_min, sf_min ), */
/*             where i = 1,2,...,N, k = 1,2,...,N-1, and sf_min is at */
/*             least the smallest number that can divide one without */
/*             overflow (see LAPACK Library routine DLAMCH). */
/*             Note that this condition is not checked by the routine. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             This parameter defines the multiplicity of singular values */
/*             by considering all singular values within an interval of */
/*             length TOL as coinciding. TOL is used in checking how many */
/*             singular values are less than or equal to THETA. Also in */
/*             computing an appropriate upper bound THETA by a bisection */
/*             method, TOL is used as a stopping criterion defining the */
/*             minimum (absolute) subinterval width.  TOL >= 0. */

/*     RELTOL  DOUBLE PRECISION */
/*             This parameter specifies the minimum relative width of an */
/*             interval. When an interval is narrower than TOL, or than */
/*             RELTOL times the larger (in magnitude) endpoint, then it */
/*             is considered to be sufficiently small and bisection has */
/*             converged. */
/*             RELTOL >= BASE * EPS, where BASE is machine radix and EPS */
/*             is machine precision (see LAPACK Library routine DLAMCH). */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warnings; */
/*             = 1:  if the value of L has been increased as the L-th */
/*                   smallest singular value of J coincides with the */
/*                   (L+1)-th smallest one. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Let s(i), i = 1,2,...,N, be the N non-negative singular values of */
/*     the bidiagonal matrix J arranged so that s(1) >= ... >= s(N) >= 0. */
/*     The routine then computes an upper bound T such that s(N-L) > T >= */
/*     s(N-L+1) as follows (see [2]). */
/*     First, if the initial estimate of THETA is not specified by the */
/*     user then the routine initialises THETA to be an estimate which */
/*     is close to the requested value of THETA if s(N-L) >> s(N-L+1). */
/*     Second, a bisection method (see [1, 8.5]) is used which generates */
/*     a sequence of shrinking intervals [Y,Z] such that either THETA in */
/*     [Y,Z] was found (so that J has L singular values less than or */
/*     equal to THETA), or */

/*        (number of s(i) <= Y) < L < (number of s(i) <= Z). */

/*     This bisection method is applied to an associated 2N-by-2N */
/*     symmetric tridiagonal matrix T" whose eigenvalues (see [1]) are */
/*     given by s(1),s(2),...,s(N),-s(1),-s(2),...,-s(N). One of the */
/*     starting values for the bisection method is the initial value of */
/*     THETA. If this value is an upper bound, then the initial lower */
/*     bound is set to zero, else the initial upper bound is computed */
/*     from the Gershgorin Circle Theorem [1, Theorem 7.2-1], applied to */
/*     T". The computation of the "number of s(i) <= Y (or Z)" is */
/*     achieved by calling SLICOT Library routine MB03ND, which applies */
/*     Sylvester's Law of Inertia or equivalently Sturm sequences */
/*     [1, 8.5] to the associated matrix T". If */

/*        Z - Y <= MAX( TOL, PIVMIN, RELTOL*MAX( ABS( Y ), ABS( Z ) ) ) */

/*     at some stage of the bisection method, then at least two singular */
/*     values of J lie in the interval [Y,Z] within a distance less than */
/*     TOL from each other. In this case, s(N-L) and s(N-L+1) are assumed */
/*     to coincide, the upper bound T is set to the value of Z, the value */
/*     of L is increased and IWARN is set to 1. */

/*     REFERENCES */

/*     [1] Golub, G.H. and Van Loan, C.F. */
/*         Matrix Computations. */
/*         The Johns Hopkins University Press, Baltimore, Maryland, 1983. */

/*     [2] Van Huffel, S. and Vandewalle, J. */
/*         The Partial Total Least Squares Algorithm. */
/*         J. Comput. and Appl. Math., 21, pp. 333-341, 1988. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997. */
/*     Supersedes Release 2.0 routine MB03AD by S. Van Huffel, Katholieke */
/*     University, Leuven, Belgium. */

/*     REVISIONS */

/*     June 16, 1997, Oct. 26, 2003. */

/*     KEYWORDS */

/*     Bidiagonal matrix, singular values. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Test some input scalar arguments. */

#line 226 "MB03MD.f"
    /* Parameter adjustments */
#line 226 "MB03MD.f"
    --e2;
#line 226 "MB03MD.f"
    --q2;
#line 226 "MB03MD.f"
    --e;
#line 226 "MB03MD.f"
    --q;
#line 226 "MB03MD.f"

#line 226 "MB03MD.f"
    /* Function Body */
#line 226 "MB03MD.f"
    *iwarn = 0;
#line 227 "MB03MD.f"
    *info = 0;
#line 228 "MB03MD.f"
    if (*n < 0) {
#line 229 "MB03MD.f"
	*info = -1;
#line 230 "MB03MD.f"
    } else if (*l < 0 || *l > *n) {
#line 231 "MB03MD.f"
	*info = -2;
#line 232 "MB03MD.f"
    }

#line 234 "MB03MD.f"
    if (*info != 0) {

/*        Error return. */

#line 238 "MB03MD.f"
	i__1 = -(*info);
#line 238 "MB03MD.f"
	xerbla_("MB03MD", &i__1, (ftnlen)6);
#line 239 "MB03MD.f"
	return 0;
#line 240 "MB03MD.f"
    }

/*     Quick return if possible. */

#line 244 "MB03MD.f"
    if (*n == 0) {
#line 244 "MB03MD.f"
	return 0;
#line 244 "MB03MD.f"
    }

/*     Step 1: initialisation of THETA. */
/*             ----------------------- */
#line 249 "MB03MD.f"
    if (*l == 0) {
#line 249 "MB03MD.f"
	*theta = 0.;
#line 249 "MB03MD.f"
    }
#line 250 "MB03MD.f"
    if (*theta < 0.) {
#line 251 "MB03MD.f"
	if (*l == 1) {

/*           An upper bound which is close if S(N-1) >> S(N): */

#line 255 "MB03MD.f"
	    *theta = mb03my_(n, &q[1], &c__1);
#line 256 "MB03MD.f"
	    if (*n == 1) {
#line 256 "MB03MD.f"
		return 0;
#line 256 "MB03MD.f"
	    }
#line 258 "MB03MD.f"
	} else {

/*           An experimentally established estimate which is good if */
/*           S(N-L) >> S(N-L+1): */

#line 263 "MB03MD.f"
	    *theta = (d__1 = q[*n - *l + 1], abs(d__1));
#line 264 "MB03MD.f"
	}
#line 265 "MB03MD.f"
    }

/*     Step 2: Check quality of initial estimate THETA. */
/*             --------------------------------------- */
#line 269 "MB03MD.f"
    num = mb03nd_(n, theta, &q2[1], &e2[1], pivmin, info);
#line 270 "MB03MD.f"
    if (num == *l) {
#line 270 "MB03MD.f"
	return 0;
#line 270 "MB03MD.f"
    }

/*     Step 3: initialisation starting values for bisection method. */
/*             --------------------------------------------------- */
/*     Let S(i), i=1,...,N, be the singular values of J in decreasing */
/*     order. Then, the computed Y and Z will be such that */
/*     (number of S(i) <= Y) < L < (number of S(i) <= Z). */

#line 279 "MB03MD.f"
    if (num < *l) {
#line 280 "MB03MD.f"
	th = abs(q[1]);
#line 281 "MB03MD.f"
	z__ = 0.;
#line 282 "MB03MD.f"
	y = *theta;
#line 283 "MB03MD.f"
	numz = *n;

#line 285 "MB03MD.f"
	i__1 = *n - 1;
#line 285 "MB03MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 286 "MB03MD.f"
	    h__ = (d__1 = q[i__ + 1], abs(d__1));
/* Computing MAX */
#line 287 "MB03MD.f"
	    d__2 = max(th,h__) + (d__1 = e[i__], abs(d__1));
#line 287 "MB03MD.f"
	    z__ = max(d__2,z__);
#line 288 "MB03MD.f"
	    th = h__;
#line 289 "MB03MD.f"
/* L20: */
#line 289 "MB03MD.f"
	}

/*        Widen the Gershgorin interval a bit for machines with sloppy */
/*        arithmetic. */

#line 294 "MB03MD.f"
	z__ = z__ + abs(z__) * 2. * dlamch_("Epsilon", (ftnlen)7) * (
		doublereal) (*n) + *pivmin * 2.;
#line 296 "MB03MD.f"
    } else {
#line 297 "MB03MD.f"
	z__ = *theta;
#line 298 "MB03MD.f"
	y = 0.;
#line 299 "MB03MD.f"
	numz = num;
#line 300 "MB03MD.f"
    }

/*     Step 4: Bisection method for finding the upper bound on the L */
/*             smallest singular values of the bidiagonal. */
/*             ------------------------------------------ */
/*     A sequence of subintervals [Y,Z] is produced such that */
/*         (number of S(i) <= Y) < L < (number of S(i) <= Z). */
/*     NUM : number of S(i) <= TH, */
/*     NUMZ: number of S(i) <= Z. */

/*     WHILE ( ( NUM .NE. L ) .AND. */
/*             ( ( Z-Y ) .GT. MAX( TOL, PIVMIN, RELTOL*ABS( Z ) ) ) ) DO */
#line 312 "MB03MD.f"
L40:
/* Computing MAX */
/* Computing MAX */
#line 312 "MB03MD.f"
    d__4 = abs(y), d__5 = abs(z__);
#line 312 "MB03MD.f"
    d__2 = max(*tol,*pivmin), d__3 = *reltol * max(d__4,d__5);
#line 312 "MB03MD.f"
    if (num != *l && (d__1 = z__ - y, abs(d__1)) > max(d__2,d__3)) {
#line 316 "MB03MD.f"
	th = (y + z__) / 2.;
#line 317 "MB03MD.f"
	num = mb03nd_(n, &th, &q2[1], &e2[1], pivmin, info);
#line 318 "MB03MD.f"
	if (num < *l) {
#line 319 "MB03MD.f"
	    y = th;
#line 320 "MB03MD.f"
	} else {
#line 321 "MB03MD.f"
	    z__ = th;
#line 322 "MB03MD.f"
	    numz = num;
#line 323 "MB03MD.f"
	}
#line 324 "MB03MD.f"
	goto L40;
#line 325 "MB03MD.f"
    }
/*     END WHILE 40 */

/*     If NUM <> L and ( Z - Y ) <= TOL, then at least two singular */
/*     values of J lie in the interval [Y,Z] within a distance less than */
/*     TOL from each other. S(N-L) and S(N-L+1) are then assumed to */
/*     coincide. L is increased, and a warning is given. */

#line 333 "MB03MD.f"
    if (num != *l) {
#line 334 "MB03MD.f"
	*l = numz;
#line 335 "MB03MD.f"
	*theta = z__;
#line 336 "MB03MD.f"
	*iwarn = 1;
#line 337 "MB03MD.f"
    } else {
#line 338 "MB03MD.f"
	*theta = th;
#line 339 "MB03MD.f"
    }

#line 341 "MB03MD.f"
    return 0;
/* *** Last line of MB03MD *** */
} /* mb03md_ */
