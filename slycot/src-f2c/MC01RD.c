#line 1 "MC01RD.f"
/* MC01RD.f -- translated by f2c (version 20100827).
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

#line 1 "MC01RD.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;

/* Subroutine */ int mc01rd_(integer *dp1, integer *dp2, integer *dp3, 
	doublereal *alpha, doublereal *p1, doublereal *p2, doublereal *p3, 
	integer *info)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l, d1, d2, d3, e3, dmin__, dmax__;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer dsum;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *), daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), xerbla_(char *, integer *, 
	    ftnlen);


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

/*     To compute the coefficients of the polynomial */

/*        P(x) = P1(x) * P2(x) + alpha * P3(x), */

/*     where P1(x), P2(x) and P3(x) are given real polynomials and alpha */
/*     is a real scalar. */

/*     Each of the polynomials P1(x), P2(x) and P3(x) may be the zero */
/*     polynomial. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     DP1     (input) INTEGER */
/*             The degree of the polynomial P1(x).  DP1 >= -1. */

/*     DP2     (input) INTEGER */
/*             The degree of the polynomial P2(x).  DP2 >= -1. */

/*     DP3     (input/output) INTEGER */
/*             On entry, the degree of the polynomial P3(x).  DP3 >= -1. */
/*             On exit, the degree of the polynomial P(x). */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar value alpha of the problem. */

/*     P1      (input) DOUBLE PRECISION array, dimension (lenp1) */
/*             where lenp1 = DP1 + 1 if DP1 >= 0 and 1 otherwise. */
/*             If DP1 >= 0, then this array must contain the */
/*             coefficients of P1(x) in increasing powers of x. */
/*             If DP1 = -1, then P1(x) is taken to be the zero */
/*             polynomial, P1 is not referenced and can be supplied */
/*             as a dummy array. */

/*     P2      (input) DOUBLE PRECISION array, dimension (lenp2) */
/*             where lenp2 = DP2 + 1 if DP2 >= 0 and 1 otherwise. */
/*             If DP2 >= 0, then this array must contain the */
/*             coefficients of P2(x) in increasing powers of x. */
/*             If DP2 = -1, then P2(x) is taken to be the zero */
/*             polynomial, P2 is not referenced and can be supplied */
/*             as a dummy array. */

/*     P3      (input/output) DOUBLE PRECISION array, dimension (lenp3) */
/*             where lenp3 = MAX(DP1+DP2,DP3,0) + 1. */
/*             On entry, if DP3 >= 0, then this array must contain the */
/*             coefficients of P3(x) in increasing powers of x. */
/*             On entry, if DP3 = -1, then P3(x) is taken to be the zero */
/*             polynomial. */
/*             On exit, the leading (DP3+1) elements of this array */
/*             contain the coefficients of P(x) in increasing powers of x */
/*             unless DP3 = -1 on exit, in which case the coefficients of */
/*             P(x) (the zero polynomial) are not stored in the array. */
/*             This is the case, for instance, when ALPHA = 0.0 and */
/*             P1(x) or P2(x) is the zero polynomial. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Given real polynomials */

/*                DP1           i           DP2           i */
/*        P1(x) = SUM a(i+1) * x ,  P2(x) = SUM b(i+1) * x  and */
/*                i=0                       i=0 */

/*                DP3           i */
/*        P3(x) = SUM c(i+1) * x , */
/*                i=0 */

/*     the routine computes the coefficents of P(x) = P1(x) * P2(x) + */
/*                     DP3            i */
/*     alpha * P3(x) = SUM  d(i+1) * x  as follows. */
/*                     i=0 */

/*     Let e(i) = c(i) for 1 <= i <= DP3+1 and e(i) = 0 for i > DP3+1. */
/*     Then if DP1 >= DP2, */

/*                i */
/*        d(i) = SUM a(k) * b(i-k+1) + f(i), for i = 1, ..., DP2+1, */
/*               k=1 */

/*                 i */
/*        d(i)  = SUM a(k) * b(i-k+1) + f(i), for i = DP2+2, ..., DP1+1 */
/*               k=i-DP2 */

/*     and */
/*                DP1+1 */
/*        d(i)  = SUM a(k) * b(i-k+1) + f(i) for i = DP1+2,...,DP1+DP2+1, */
/*               k=i-DP2 */

/*     where f(i) = alpha * e(i). */

/*     Similar formulas hold for the case DP1 < DP2. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC01FD by C. Klimann and */
/*     A.J. Geurts. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary polynomial operations, polynomial operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

#line 164 "MC01RD.f"
    /* Parameter adjustments */
#line 164 "MC01RD.f"
    --p3;
#line 164 "MC01RD.f"
    --p2;
#line 164 "MC01RD.f"
    --p1;
#line 164 "MC01RD.f"

#line 164 "MC01RD.f"
    /* Function Body */
#line 164 "MC01RD.f"
    *info = 0;
#line 165 "MC01RD.f"
    if (*dp1 < -1) {
#line 166 "MC01RD.f"
	*info = -1;
#line 167 "MC01RD.f"
    } else if (*dp2 < -1) {
#line 168 "MC01RD.f"
	*info = -2;
#line 169 "MC01RD.f"
    } else if (*dp3 < -1) {
#line 170 "MC01RD.f"
	*info = -3;
#line 171 "MC01RD.f"
    }

#line 173 "MC01RD.f"
    if (*info != 0) {

/*        Error return. */

#line 177 "MC01RD.f"
	i__1 = -(*info);
#line 177 "MC01RD.f"
	xerbla_("MC01RD", &i__1, (ftnlen)6);
#line 178 "MC01RD.f"
	return 0;
#line 179 "MC01RD.f"
    }

/*     Computation of the exact degree of the polynomials, i.e., Di such */
/*     that either Di = -1 or Pi(Di+1) is non-zero. */

#line 184 "MC01RD.f"
    d1 = *dp1;
/*     WHILE ( D1 >= 0 and P1(D1+1) = 0 ) DO */
#line 186 "MC01RD.f"
L20:
#line 186 "MC01RD.f"
    if (d1 >= 0) {
#line 187 "MC01RD.f"
	if (p1[d1 + 1] == 0.) {
#line 188 "MC01RD.f"
	    --d1;
#line 189 "MC01RD.f"
	    goto L20;
#line 190 "MC01RD.f"
	}
#line 191 "MC01RD.f"
    }
/*     END WHILE 20 */
#line 193 "MC01RD.f"
    d2 = *dp2;
/*     WHILE ( D2 >= 0 and P2(D2+1) = 0 ) DO */
#line 195 "MC01RD.f"
L40:
#line 195 "MC01RD.f"
    if (d2 >= 0) {
#line 196 "MC01RD.f"
	if (p2[d2 + 1] == 0.) {
#line 197 "MC01RD.f"
	    --d2;
#line 198 "MC01RD.f"
	    goto L40;
#line 199 "MC01RD.f"
	}
#line 200 "MC01RD.f"
    }
/*     END WHILE 40 */
#line 202 "MC01RD.f"
    if (*alpha == 0.) {
#line 203 "MC01RD.f"
	d3 = -1;
#line 204 "MC01RD.f"
    } else {
#line 205 "MC01RD.f"
	d3 = *dp3;
#line 206 "MC01RD.f"
    }
/*     WHILE ( D3 >= 0 and P3(D3+1) = 0 ) DO */
#line 208 "MC01RD.f"
L60:
#line 208 "MC01RD.f"
    if (d3 >= 0) {
#line 209 "MC01RD.f"
	if (p3[d3 + 1] == 0.) {
#line 210 "MC01RD.f"
	    --d3;
#line 211 "MC01RD.f"
	    goto L60;
#line 212 "MC01RD.f"
	}
#line 213 "MC01RD.f"
    }
/*     END WHILE 60 */

/*     Computation of P3(x) := ALPHA * P3(x). */

#line 218 "MC01RD.f"
    i__1 = d3 + 1;
#line 218 "MC01RD.f"
    dscal_(&i__1, alpha, &p3[1], &c__1);

#line 220 "MC01RD.f"
    if (d1 == -1 || d2 == -1) {
#line 221 "MC01RD.f"
	*dp3 = d3;
#line 222 "MC01RD.f"
	return 0;
#line 223 "MC01RD.f"
    }

/*     P1(x) and P2(x) are non-zero polynomials. */

#line 227 "MC01RD.f"
    dsum = d1 + d2;
#line 228 "MC01RD.f"
    dmax__ = max(d1,d2);
#line 229 "MC01RD.f"
    dmin__ = dsum - dmax__;

#line 231 "MC01RD.f"
    if (d3 < dsum) {
#line 232 "MC01RD.f"
	p3[d3 + 2] = 0.;
#line 233 "MC01RD.f"
	i__1 = dsum - d3 - 1;
#line 233 "MC01RD.f"
	dcopy_(&i__1, &p3[d3 + 2], &c__0, &p3[d3 + 3], &c__1);
#line 234 "MC01RD.f"
	d3 = dsum;
#line 235 "MC01RD.f"
    }

#line 237 "MC01RD.f"
    if (d1 == 0 || d2 == 0) {

/*        D1 or D2 is zero. */

#line 241 "MC01RD.f"
	if (d1 != 0) {
#line 242 "MC01RD.f"
	    i__1 = d1 + 1;
#line 242 "MC01RD.f"
	    daxpy_(&i__1, &p2[1], &p1[1], &c__1, &p3[1], &c__1);
#line 243 "MC01RD.f"
	} else {
#line 244 "MC01RD.f"
	    i__1 = d2 + 1;
#line 244 "MC01RD.f"
	    daxpy_(&i__1, &p1[1], &p2[1], &c__1, &p3[1], &c__1);
#line 245 "MC01RD.f"
	}
#line 246 "MC01RD.f"
    } else {

/*        D1 and D2 are both nonzero. */

/*        First part of the computation. */

#line 252 "MC01RD.f"
	i__1 = dmin__ + 1;
#line 252 "MC01RD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 253 "MC01RD.f"
	    p3[i__] += ddot_(&i__, &p1[1], &c__1, &p2[1], &c_n1);
#line 254 "MC01RD.f"
/* L80: */
#line 254 "MC01RD.f"
	}

/*        Second part of the computation. */

#line 258 "MC01RD.f"
	i__1 = dmax__ + 1;
#line 258 "MC01RD.f"
	for (i__ = dmin__ + 2; i__ <= i__1; ++i__) {
#line 259 "MC01RD.f"
	    if (d1 > d2) {
#line 260 "MC01RD.f"
		k = i__ - d2;
#line 261 "MC01RD.f"
		i__2 = dmin__ + 1;
#line 261 "MC01RD.f"
		p3[i__] += ddot_(&i__2, &p1[k], &c__1, &p2[1], &c_n1);
#line 262 "MC01RD.f"
	    } else {
#line 263 "MC01RD.f"
		k = i__ - d1;
#line 264 "MC01RD.f"
		i__2 = dmin__ + 1;
#line 264 "MC01RD.f"
		p3[i__] += ddot_(&i__2, &p2[k], &c_n1, &p1[1], &c__1);
#line 265 "MC01RD.f"
	    }
#line 266 "MC01RD.f"
/* L100: */
#line 266 "MC01RD.f"
	}

/*        Third part of the computation. */

#line 270 "MC01RD.f"
	e3 = dsum + 2;

#line 272 "MC01RD.f"
	i__1 = dsum + 1;
#line 272 "MC01RD.f"
	for (i__ = dmax__ + 2; i__ <= i__1; ++i__) {
#line 273 "MC01RD.f"
	    j = e3 - i__;
#line 274 "MC01RD.f"
	    k = i__ - dmin__;
#line 275 "MC01RD.f"
	    l = i__ - dmax__;
#line 276 "MC01RD.f"
	    if (d1 > d2) {
#line 277 "MC01RD.f"
		p3[i__] += ddot_(&j, &p1[k], &c__1, &p2[l], &c_n1);
#line 278 "MC01RD.f"
	    } else {
#line 279 "MC01RD.f"
		p3[i__] += ddot_(&j, &p1[l], &c_n1, &p2[k], &c__1);
#line 280 "MC01RD.f"
	    }
#line 281 "MC01RD.f"
/* L120: */
#line 281 "MC01RD.f"
	}

#line 283 "MC01RD.f"
    }

/*     Computation of the exact degree of P3(x). */

/*     WHILE ( D3 >= 0 and P3(D3+1) = 0 ) DO */
#line 288 "MC01RD.f"
L140:
#line 288 "MC01RD.f"
    if (d3 >= 0) {
#line 289 "MC01RD.f"
	if (p3[d3 + 1] == 0.) {
#line 290 "MC01RD.f"
	    --d3;
#line 291 "MC01RD.f"
	    goto L140;
#line 292 "MC01RD.f"
	}
#line 293 "MC01RD.f"
    }
/*     END WHILE 140 */
#line 295 "MC01RD.f"
    *dp3 = d3;

#line 297 "MC01RD.f"
    return 0;
/* *** Last line of MC01RD *** */
} /* mc01rd_ */

