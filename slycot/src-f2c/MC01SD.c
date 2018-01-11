#line 1 "MC01SD.f"
/* MC01SD.f -- translated by f2c (version 20100827).
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

#line 1 "MC01SD.f"
/* Subroutine */ int mc01sd_(integer *dp, doublereal *p, integer *s, integer *
	t, doublereal *mant, integer *e, integer *iwork, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer i_dnnt(doublereal *);

    /* Local variables */
    static integer i__, j, m, v0, v1, lb, ub, dv, inc, beta;
    extern /* Subroutine */ int mc01sw_(doublereal *, integer *, doublereal *,
	     integer *);
    extern integer mc01sx_(integer *, integer *, integer *, doublereal *);
    extern /* Subroutine */ int mc01sy_(doublereal *, integer *, integer *, 
	    doublereal *, logical *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical ovflow;


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

/*     To scale the coefficients of the real polynomial P(x) such that */
/*     the coefficients of the scaled polynomial Q(x) = sP(tx) have */
/*     minimal variation, where s and t are real scalars. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     DP      (input) INTEGER */
/*             The degree of the polynomial P(x).  DP >= 0. */

/*     P       (input/output) DOUBLE PRECISION array, dimension (DP+1) */
/*             On entry, this array must contain the coefficients of P(x) */
/*             in increasing powers of x. */
/*             On exit, this array contains the coefficients of the */
/*             scaled polynomial Q(x) in increasing powers of x. */

/*     S       (output) INTEGER */
/*             The exponent of the floating-point representation of the */
/*             scaling factor s = BASE**S, where BASE is the base of the */
/*             machine representation of floating-point numbers (see */
/*             LAPACK Library Routine DLAMCH). */

/*     T       (output) INTEGER */
/*             The exponent of the floating-point representation of the */
/*             scaling factor t = BASE**T. */

/*     MANT    (output) DOUBLE PRECISION array, dimension (DP+1) */
/*             This array contains the mantissas of the standard */
/*             floating-point representation of the coefficients of the */
/*             scaled polynomial Q(x) in increasing powers of x. */

/*     E       (output) INTEGER array, dimension (DP+1) */
/*             This array contains the exponents of the standard */
/*             floating-point representation of the coefficients of the */
/*             scaled polynomial Q(x) in increasing powers of x. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (DP+1) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if on entry, P(x) is the zero polynomial. */

/*     METHOD */

/*     Define the variation of the coefficients of the real polynomial */

/*                                         2                DP */
/*        P(x) = p(0) + p(1) * x + p(2) * x  + ... + p(DP) x */

/*     whose non-zero coefficients can be represented as */
/*                          e(i) */
/*        p(i) = m(i) * BASE     (where 1 <= ABS(m(i)) < BASE) */

/*     by */

/*        V = max(e(i)) - min(e(i)), */

/*     where max and min are taken over the indices i for which p(i) is */
/*     non-zero. */
/*                                        DP         i    i */
/*     For the scaled polynomial P(cx) = SUM p(i) * c  * x  with */
/*                                       i=0 */
/*                j */
/*     c  = (BASE) , the variation V(j) is given by */

/*       V(j) = max(e(i) + j * i) - min(e(i) + j * i). */

/*     Using the fact that V(j) is a convex function of j, the routine */
/*     determines scaling factors s = (BASE)**S and t = (BASE)**T such */
/*     that the coefficients of the scaled polynomial Q(x) = sP(tx) */
/*     satisfy the following conditions: */

/*       (a) 1 <= q(0) < BASE and */

/*       (b) the variation of the coefficients of Q(x) is minimal. */

/*     Further details can be found in [1]. */

/*     REFERENCES */

/*     [1] Dunaway, D.K. */
/*         Calculation of Zeros of a Real Polynomial through */
/*         Factorization using Euclid's Algorithm. */
/*         SIAM J. Numer. Anal., 11, pp. 1087-1104, 1974. */

/*     NUMERICAL ASPECTS */

/*     Since the scaling is performed on the exponents of the floating- */
/*     point representation of the coefficients of P(x), no rounding */
/*     errors occur during the computation of the coefficients of Q(x). */

/*     FURTHER COMMENTS */

/*     The scaling factors s and t are BASE dependent. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC01GD by A.J. Geurts. */

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

#line 163 "MC01SD.f"
    /* Parameter adjustments */
#line 163 "MC01SD.f"
    --iwork;
#line 163 "MC01SD.f"
    --e;
#line 163 "MC01SD.f"
    --mant;
#line 163 "MC01SD.f"
    --p;
#line 163 "MC01SD.f"

#line 163 "MC01SD.f"
    /* Function Body */
#line 163 "MC01SD.f"
    if (*dp < 0) {
#line 164 "MC01SD.f"
	*info = -1;

/*        Error return. */

#line 168 "MC01SD.f"
	i__1 = -(*info);
#line 168 "MC01SD.f"
	xerbla_("MC01SD", &i__1, (ftnlen)6);
#line 169 "MC01SD.f"
	return 0;
#line 170 "MC01SD.f"
    }

#line 172 "MC01SD.f"
    *info = 0;
#line 173 "MC01SD.f"
    lb = 1;
/*     WHILE ( LB <= DP+1 and P(LB) = 0 ) DO */
#line 175 "MC01SD.f"
L20:
#line 175 "MC01SD.f"
    if (lb <= *dp + 1) {
#line 176 "MC01SD.f"
	if (p[lb] == 0.) {
#line 177 "MC01SD.f"
	    ++lb;
#line 178 "MC01SD.f"
	    goto L20;
#line 179 "MC01SD.f"
	}
#line 180 "MC01SD.f"
    }
/*     END WHILE 20 */

/*     LB = MIN( i: P(i) non-zero). */

#line 185 "MC01SD.f"
    if (lb == *dp + 2) {
#line 186 "MC01SD.f"
	*info = 1;
#line 187 "MC01SD.f"
	return 0;
#line 188 "MC01SD.f"
    }

#line 190 "MC01SD.f"
    ub = *dp + 1;
/*     WHILE ( P(UB) = 0 ) DO */
#line 192 "MC01SD.f"
L40:
#line 192 "MC01SD.f"
    if (p[ub] == 0.) {
#line 193 "MC01SD.f"
	--ub;
#line 194 "MC01SD.f"
	goto L40;
#line 195 "MC01SD.f"
    }
/*     END WHILE 40 */

/*     UB = MAX(i: P(i) non-zero). */

#line 200 "MC01SD.f"
    beta = (integer) dlamch_("Base", (ftnlen)4);

#line 202 "MC01SD.f"
    i__1 = *dp + 1;
#line 202 "MC01SD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 203 "MC01SD.f"
	mc01sw_(&p[i__], &beta, &mant[i__], &e[i__]);
#line 204 "MC01SD.f"
/* L60: */
#line 204 "MC01SD.f"
    }

/*     First prescaling. */

#line 208 "MC01SD.f"
    m = e[lb];
#line 209 "MC01SD.f"
    if (m != 0) {

#line 211 "MC01SD.f"
	i__1 = ub;
#line 211 "MC01SD.f"
	for (i__ = lb; i__ <= i__1; ++i__) {
#line 212 "MC01SD.f"
	    if (mant[i__] != 0.) {
#line 212 "MC01SD.f"
		e[i__] -= m;
#line 212 "MC01SD.f"
	    }
#line 213 "MC01SD.f"
/* L80: */
#line 213 "MC01SD.f"
	}

#line 215 "MC01SD.f"
    }
#line 216 "MC01SD.f"
    *s = -m;

/*     Second prescaling. */

#line 220 "MC01SD.f"
    if (ub > 1) {
#line 220 "MC01SD.f"
	d__1 = (doublereal) e[ub] / (doublereal) (ub - 1);
#line 220 "MC01SD.f"
	m = i_dnnt(&d__1);
#line 220 "MC01SD.f"
    }

#line 222 "MC01SD.f"
    i__1 = ub;
#line 222 "MC01SD.f"
    for (i__ = lb; i__ <= i__1; ++i__) {
#line 223 "MC01SD.f"
	if (mant[i__] != 0.) {
#line 223 "MC01SD.f"
	    e[i__] -= m * (i__ - 1);
#line 223 "MC01SD.f"
	}
#line 224 "MC01SD.f"
/* L100: */
#line 224 "MC01SD.f"
    }

#line 226 "MC01SD.f"
    *t = -m;

#line 228 "MC01SD.f"
    v0 = mc01sx_(&lb, &ub, &e[1], &mant[1]);
#line 229 "MC01SD.f"
    j = 1;

#line 231 "MC01SD.f"
    i__1 = ub;
#line 231 "MC01SD.f"
    for (i__ = lb; i__ <= i__1; ++i__) {
#line 232 "MC01SD.f"
	if (mant[i__] != 0.) {
#line 232 "MC01SD.f"
	    iwork[i__] = e[i__] + (i__ - 1);
#line 232 "MC01SD.f"
	}
#line 233 "MC01SD.f"
/* L120: */
#line 233 "MC01SD.f"
    }

#line 235 "MC01SD.f"
    v1 = mc01sx_(&lb, &ub, &iwork[1], &mant[1]);
#line 236 "MC01SD.f"
    dv = v1 - v0;
#line 237 "MC01SD.f"
    if (dv != 0) {
#line 238 "MC01SD.f"
	if (dv > 0) {
#line 239 "MC01SD.f"
	    j = 0;
#line 240 "MC01SD.f"
	    inc = -1;
#line 241 "MC01SD.f"
	    v1 = v0;
#line 242 "MC01SD.f"
	    dv = -dv;

#line 244 "MC01SD.f"
	    i__1 = ub;
#line 244 "MC01SD.f"
	    for (i__ = lb; i__ <= i__1; ++i__) {
#line 245 "MC01SD.f"
		iwork[i__] = e[i__];
#line 246 "MC01SD.f"
/* L130: */
#line 246 "MC01SD.f"
	    }

#line 248 "MC01SD.f"
	} else {
#line 249 "MC01SD.f"
	    inc = 1;
#line 250 "MC01SD.f"
	}
/*        WHILE ( DV < 0 ) DO */
#line 252 "MC01SD.f"
L140:
#line 252 "MC01SD.f"
	if (dv < 0) {
#line 253 "MC01SD.f"
	    v0 = v1;

#line 255 "MC01SD.f"
	    i__1 = ub;
#line 255 "MC01SD.f"
	    for (i__ = lb; i__ <= i__1; ++i__) {
#line 256 "MC01SD.f"
		e[i__] = iwork[i__];
#line 257 "MC01SD.f"
/* L150: */
#line 257 "MC01SD.f"
	    }

#line 259 "MC01SD.f"
	    j += inc;

#line 261 "MC01SD.f"
	    i__1 = ub;
#line 261 "MC01SD.f"
	    for (i__ = lb; i__ <= i__1; ++i__) {
#line 262 "MC01SD.f"
		iwork[i__] = e[i__] + inc * (i__ - 1);
#line 263 "MC01SD.f"
/* L160: */
#line 263 "MC01SD.f"
	    }

#line 265 "MC01SD.f"
	    v1 = mc01sx_(&lb, &ub, &iwork[1], &mant[1]);
#line 266 "MC01SD.f"
	    dv = v1 - v0;
#line 267 "MC01SD.f"
	    goto L140;
#line 268 "MC01SD.f"
	}
/*        END WHILE 140 */
#line 270 "MC01SD.f"
	*t = *t + j - inc;
#line 271 "MC01SD.f"
    }

/*     Evaluation of the output parameters. */

#line 275 "MC01SD.f"
    i__1 = ub;
#line 275 "MC01SD.f"
    for (i__ = lb; i__ <= i__1; ++i__) {
#line 276 "MC01SD.f"
	mc01sy_(&mant[i__], &e[i__], &beta, &p[i__], &ovflow);
#line 277 "MC01SD.f"
/* L180: */
#line 277 "MC01SD.f"
    }

#line 279 "MC01SD.f"
    return 0;
/* *** Last line of MC01SD *** */
} /* mc01sd_ */

