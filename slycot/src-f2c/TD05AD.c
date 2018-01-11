#line 1 "TD05AD.f"
/* TD05AD.f -- translated by f2c (version 20100827).
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

#line 1 "TD05AD.f"
/* Table of constant values */

static doublereal c_b13 = 90.;

/* Subroutine */ int td05ad_(char *unitf, char *output, integer *np1, integer 
	*mp1, doublereal *w, doublereal *a, doublereal *b, doublereal *valr, 
	doublereal *vali, integer *info, ftnlen unitf_len, ftnlen output_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double atan(doublereal), pow_di(doublereal *, integer *), d_imag(
	    doublecomplex *), d_sign(doublereal *, doublereal *), d_lg10(
	    doublereal *);

    /* Local variables */
    static doublereal g;
    static integer i__, m, n, m2, n2;
    static doublereal w2, wc, bimag, breal, timag;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal treal;
    static doublecomplex ztemp;
    static doublereal twopi;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static integer iphase;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern /* Double Complex */ VOID zladiv_(doublecomplex *, doublecomplex *,
	     doublecomplex *);
    static logical lunitf;
    static integer npzero, nzzero;
    static logical loutpu;


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

/*     Given a complex valued rational function of frequency (transfer */
/*     function) G(jW) this routine will calculate its complex value or */
/*     its magnitude and phase for a specified frequency value. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UNITF   CHARACTER*1 */
/*             Indicates the choice of frequency unit as follows: */
/*             = 'R':  Input frequency W in radians/second; */
/*             = 'H':  Input frequency W in hertz. */

/*     OUTPUT  CHARACTER*1 */
/*             Indicates the choice of co-ordinates for output as folows: */
/*             = 'C':  Cartesian co-ordinates (output real and imaginary */
/*                     parts of G(jW)); */
/*             = 'P':  Polar co-ordinates (output magnitude and phase */
/*                     of G(jW)). */

/*     Input/Output Parameters */

/*     NP1     (input) INTEGER */
/*             The order of the denominator + 1, i.e. N + 1.  NP1 >= 1. */

/*     MP1     (input) INTEGER */
/*             The order of the numerator + 1, i.e. M + 1.  MP1 >= 1. */

/*     W       (input) DOUBLE PRECISION */
/*             The frequency value W for which the transfer function is */
/*             to be evaluated. */

/*     A       (input) DOUBLE PRECISION array, dimension (NP1) */
/*             This array must contain the vector of denominator */
/*             coefficients in ascending order of powers. That is, A(i) */
/*             must contain the coefficient of (jW)**(i-1) for i = 1, */
/*             2,...,NP1. */

/*     B       (input) DOUBLE PRECISION array, dimension (MP1) */
/*             This array must contain the vector of numerator */
/*             coefficients in ascending order of powers. That is, B(i) */
/*             must contain the coefficient of (jW)**(i-1) for i = 1, */
/*             2,...,MP1. */

/*     VALR    (output) DOUBLE PRECISION */
/*             If OUTPUT = 'C', VALR contains the real part of G(jW). */
/*             If OUTPUT = 'P', VALR contains the magnitude of G(jW) */
/*                              in dBs. */

/*     VALI    (output) DOUBLE PRECISION */
/*             If OUTPUT = 'C', VALI contains the imaginary part of */
/*                              G(jW). */
/*             If OUTPUT = 'P', VALI contains the phase of G(jW) in */
/*                              degrees. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the frequency value W is a pole of G(jW), or all */
/*                   the coefficients of the A polynomial are zero. */

/*     METHOD */

/*     By substituting the values of A, B and W in the following */
/*     formula: */

/*            B(1)+B(2)*(jW)+B(3)*(jW)**2+...+B(MP1)*(jW)**(MP1-1) */
/*     G(jW) = ---------------------------------------------------. */
/*            A(1)+A(2)*(jW)+A(3)*(jW)**2+...+A(NP1)*(jW)**(NP1-1) */

/*     REFERENCES */

/*     None. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0(N+M) operations. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TD01AD by Control Systems Research */
/*     Group, Kingston Polytechnic, United Kingdom, March 1981. */

/*     REVISIONS */

/*     February 1997. */
/*     February 22, 1998 (changed the name of TD01MD). */

/*     KEYWORDS */

/*     Elementary polynomial operations, frequency response, matrix */
/*     fraction, polynomial matrix, state-space representation, transfer */
/*     matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 152 "TD05AD.f"
    /* Parameter adjustments */
#line 152 "TD05AD.f"
    --b;
#line 152 "TD05AD.f"
    --a;
#line 152 "TD05AD.f"

#line 152 "TD05AD.f"
    /* Function Body */
#line 152 "TD05AD.f"
    *info = 0;
#line 153 "TD05AD.f"
    lunitf = lsame_(unitf, "H", (ftnlen)1, (ftnlen)1);
#line 154 "TD05AD.f"
    loutpu = lsame_(output, "P", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 158 "TD05AD.f"
    if (! lunitf && ! lsame_(unitf, "R", (ftnlen)1, (ftnlen)1)) {
#line 159 "TD05AD.f"
	*info = -1;
#line 160 "TD05AD.f"
    } else if (! loutpu && ! lsame_(output, "C", (ftnlen)1, (ftnlen)1)) {
#line 161 "TD05AD.f"
	*info = -2;
#line 162 "TD05AD.f"
    } else if (*np1 < 1) {
#line 163 "TD05AD.f"
	*info = -3;
#line 164 "TD05AD.f"
    } else if (*mp1 < 1) {
#line 165 "TD05AD.f"
	*info = -4;
#line 166 "TD05AD.f"
    }

#line 168 "TD05AD.f"
    if (*info != 0) {

/*        Error return. */

#line 172 "TD05AD.f"
	i__1 = -(*info);
#line 172 "TD05AD.f"
	xerbla_("TD05AD", &i__1, (ftnlen)6);
#line 173 "TD05AD.f"
	return 0;
#line 174 "TD05AD.f"
    }

#line 176 "TD05AD.f"
    m = *mp1 - 1;
#line 177 "TD05AD.f"
    n = *np1 - 1;
#line 178 "TD05AD.f"
    wc = *w;
#line 179 "TD05AD.f"
    twopi = atan(1.) * 8.;
#line 180 "TD05AD.f"
    if (lunitf) {
#line 180 "TD05AD.f"
	wc *= twopi;
#line 180 "TD05AD.f"
    }
/* Computing 2nd power */
#line 181 "TD05AD.f"
    d__1 = wc;
#line 181 "TD05AD.f"
    w2 = d__1 * d__1;

/*     Determine the orders z (NZZERO) and p (NPZERO) of the factors */
/*     (jW)**k in the numerator and denominator polynomials, by counting */
/*     the zero trailing coefficients.  The value of G(jW) will then be */
/*     computed as (jW)**(z-p)*m(jW)/n(jW), for appropriate m and n. */

#line 188 "TD05AD.f"
    i__ = 0;

#line 190 "TD05AD.f"
L10:
#line 191 "TD05AD.f"
    ++i__;
#line 192 "TD05AD.f"
    if (i__ <= m) {
#line 193 "TD05AD.f"
	if (b[i__] == 0.) {
#line 193 "TD05AD.f"
	    goto L10;
#line 193 "TD05AD.f"
	}
#line 194 "TD05AD.f"
    }

#line 196 "TD05AD.f"
    nzzero = i__ - 1;
#line 197 "TD05AD.f"
    i__ = 0;

#line 199 "TD05AD.f"
L20:
#line 200 "TD05AD.f"
    ++i__;
#line 201 "TD05AD.f"
    if (i__ <= n) {
#line 202 "TD05AD.f"
	if (a[i__] == 0.) {
#line 202 "TD05AD.f"
	    goto L20;
#line 202 "TD05AD.f"
	}
#line 203 "TD05AD.f"
    }

#line 205 "TD05AD.f"
    npzero = i__ - 1;
#line 206 "TD05AD.f"
    iphase = nzzero - npzero;

#line 208 "TD05AD.f"
    m2 = (m - nzzero) % 2;

/*     Add real parts of the numerator m(jW). */

#line 212 "TD05AD.f"
    treal = b[*mp1 - m2];

#line 214 "TD05AD.f"
    i__1 = nzzero + 1;
#line 214 "TD05AD.f"
    for (i__ = m - 1 - m2; i__ >= i__1; i__ += -2) {
#line 215 "TD05AD.f"
	treal = b[i__] - w2 * treal;
#line 216 "TD05AD.f"
/* L30: */
#line 216 "TD05AD.f"
    }

/*     Add imaginary parts of the numerator m(jW). */

#line 220 "TD05AD.f"
    if (m == 0) {
#line 221 "TD05AD.f"
	timag = 0.;
#line 222 "TD05AD.f"
    } else {
#line 223 "TD05AD.f"
	timag = b[m + m2];

#line 225 "TD05AD.f"
	i__1 = nzzero + 2;
#line 225 "TD05AD.f"
	for (i__ = m + m2 - 2; i__ >= i__1; i__ += -2) {
#line 226 "TD05AD.f"
	    timag = b[i__] - w2 * timag;
#line 227 "TD05AD.f"
/* L40: */
#line 227 "TD05AD.f"
	}

#line 229 "TD05AD.f"
	timag *= wc;
#line 230 "TD05AD.f"
    }

#line 232 "TD05AD.f"
    n2 = (n - npzero) % 2;

/*     Add real parts of the denominator n(jW). */

#line 236 "TD05AD.f"
    breal = a[*np1 - n2];

#line 238 "TD05AD.f"
    i__1 = npzero + 1;
#line 238 "TD05AD.f"
    for (i__ = n - 1 - n2; i__ >= i__1; i__ += -2) {
#line 239 "TD05AD.f"
	breal = a[i__] - w2 * breal;
#line 240 "TD05AD.f"
/* L50: */
#line 240 "TD05AD.f"
    }

/*     Add imaginary parts of the denominator n(jW). */

#line 244 "TD05AD.f"
    if (n == 0) {
#line 245 "TD05AD.f"
	bimag = 0.;
#line 246 "TD05AD.f"
    } else {
#line 247 "TD05AD.f"
	bimag = a[n + n2];

#line 249 "TD05AD.f"
	i__1 = npzero + 2;
#line 249 "TD05AD.f"
	for (i__ = n + n2 - 2; i__ >= i__1; i__ += -2) {
#line 250 "TD05AD.f"
	    bimag = a[i__] - w2 * bimag;
#line 251 "TD05AD.f"
/* L60: */
#line 251 "TD05AD.f"
	}

#line 253 "TD05AD.f"
	bimag *= wc;
#line 254 "TD05AD.f"
    }

/* Computing MAX */
#line 256 "TD05AD.f"
    d__1 = abs(breal), d__2 = abs(bimag);
#line 256 "TD05AD.f"
    if (max(d__1,d__2) == 0. || *w == 0. && iphase < 0) {

/*        Error return:  The specified frequency W is a pole of G(jW), */
/*              or all the coefficients of the A polynomial are zero. */

#line 262 "TD05AD.f"
	*info = 1;
#line 263 "TD05AD.f"
    } else {

/*        Evaluate the complex number W**(z-p)*m(jW)/n(jW). */

#line 267 "TD05AD.f"
	z__2.r = treal, z__2.i = timag;
#line 267 "TD05AD.f"
	z__3.r = breal, z__3.i = bimag;
#line 267 "TD05AD.f"
	zladiv_(&z__1, &z__2, &z__3);
#line 267 "TD05AD.f"
	ztemp.r = z__1.r, ztemp.i = z__1.i;
#line 269 "TD05AD.f"
	*valr = ztemp.r * pow_di(&wc, &iphase);
#line 270 "TD05AD.f"
	*vali = d_imag(&ztemp) * pow_di(&wc, &iphase);

#line 272 "TD05AD.f"
	if (! loutpu) {

/*           Cartesian co-ordinates: Update the result for j**(z-p). */

#line 276 "TD05AD.f"
	    i__ = abs(iphase) % 4;
#line 277 "TD05AD.f"
	    if (iphase > 0 && i__ > 1 || iphase < 0 && (i__ == 1 || i__ == 2))
		     {
#line 279 "TD05AD.f"
		*valr = -(*valr);
#line 280 "TD05AD.f"
		*vali = -(*vali);
#line 281 "TD05AD.f"
	    }

#line 283 "TD05AD.f"
	    if (i__ % 2 != 0) {
#line 284 "TD05AD.f"
		g = *valr;
#line 285 "TD05AD.f"
		*valr = -(*vali);
#line 286 "TD05AD.f"
		*vali = g;
#line 287 "TD05AD.f"
	    }

#line 289 "TD05AD.f"
	} else {

/*           Polar co-ordinates: Compute the magnitude and phase. */

#line 293 "TD05AD.f"
	    g = dlapy2_(valr, vali);

#line 295 "TD05AD.f"
	    if (*valr == 0.) {
#line 296 "TD05AD.f"
		*vali = d_sign(&c_b13, vali);
#line 297 "TD05AD.f"
	    } else {
#line 298 "TD05AD.f"
		*vali = atan(*vali / *valr) / twopi * 360.;
#line 299 "TD05AD.f"
		if (*vali == 0. && nzzero == m && npzero == n && b[nzzero + 1]
			 * a[npzero + 1] < 0.) {
#line 299 "TD05AD.f"
		    *vali = 180.;
#line 299 "TD05AD.f"
		}
#line 302 "TD05AD.f"
	    }

#line 304 "TD05AD.f"
	    *valr = d_lg10(&g) * 20.;

#line 306 "TD05AD.f"
	    if (iphase != 0) {
#line 306 "TD05AD.f"
		*vali += (doublereal) (nzzero - npzero) * 90.;
#line 306 "TD05AD.f"
	    }
#line 308 "TD05AD.f"
	}

#line 310 "TD05AD.f"
    }

#line 312 "TD05AD.f"
    return 0;
/* *** Last line of TD05AD *** */
} /* td05ad_ */

