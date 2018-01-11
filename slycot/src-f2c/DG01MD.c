#line 1 "DG01MD.f"
/* DG01MD.f -- translated by f2c (version 20100827).
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

#line 1 "DG01MD.f"
/* Subroutine */ int dg01md_(char *indi, integer *n, doublereal *xr, 
	doublereal *xi, integer *info, ftnlen indi_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double atan(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublereal ti, wi, tr, wr, pi2;
    static logical lindi;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal whelp, wstpi, wstpr;
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

/*     To compute the discrete Fourier transform, or inverse transform, */
/*     of a complex signal. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     INDI    CHARACTER*1 */
/*             Indicates whether a Fourier transform or inverse Fourier */
/*             transform is to be performed as follows: */
/*             = 'D':  (Direct) Fourier transform; */
/*             = 'I':  Inverse Fourier transform. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of complex samples.  N must be a power of 2. */
/*             N >= 2. */

/*     XR      (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the real part of either */
/*             the complex signal z if INDI = 'D', or f(z) if INDI = 'I'. */
/*             On exit, this array contains either the real part of the */
/*             computed Fourier transform f(z) if INDI = 'D', or the */
/*             inverse Fourier transform z of f(z) if INDI = 'I'. */

/*     XI      (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the imaginary part of */
/*             either z if INDI = 'D', or f(z) if INDI = 'I'. */
/*             On exit, this array contains either the imaginary part of */
/*             f(z) if INDI = 'D', or z if INDI = 'I'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     If INDI = 'D', then the routine performs a discrete Fourier */
/*     transform on the complex signal Z(i), i = 1,2,...,N. If the result */
/*     is denoted by FZ(k), k = 1,2,...,N, then the relationship between */
/*     Z and FZ is given by the formula: */

/*                     N            ((k-1)*(i-1)) */
/*            FZ(k) = SUM ( Z(i) * V              ), */
/*                    i=1 */
/*                                     2 */
/*     where V = exp( -2*pi*j/N ) and j  = -1. */

/*     If INDI = 'I', then the routine performs an inverse discrete */
/*     Fourier transform on the complex signal FZ(k), k = 1,2,...,N. If */
/*     the result is denoted by Z(i), i = 1,2,...,N, then the */
/*     relationship between Z and FZ is given by the formula: */

/*                    N             ((k-1)*(i-1)) */
/*            Z(i) = SUM ( FZ(k) * W              ), */
/*                   k=1 */

/*     where W = exp( 2*pi*j/N ). */

/*     Note that a discrete Fourier transform, followed by an inverse */
/*     discrete Fourier transform, will result in a signal which is a */
/*     factor N larger than the original input signal. */

/*     REFERENCES */

/*     [1] Rabiner, L.R. and Rader, C.M. */
/*         Digital Signal Processing. */
/*         IEEE Press, 1972. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0( N*log(N) ) operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */
/*     Supersedes Release 2.0 routine DG01AD by R. Dekeyser, State */
/*     University of Gent, Belgium. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Complex signals, digital signal processing, fast Fourier */
/*     transform. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 139 "DG01MD.f"
    /* Parameter adjustments */
#line 139 "DG01MD.f"
    --xi;
#line 139 "DG01MD.f"
    --xr;
#line 139 "DG01MD.f"

#line 139 "DG01MD.f"
    /* Function Body */
#line 139 "DG01MD.f"
    *info = 0;
#line 140 "DG01MD.f"
    lindi = lsame_(indi, "D", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 144 "DG01MD.f"
    if (! lindi && ! lsame_(indi, "I", (ftnlen)1, (ftnlen)1)) {
#line 145 "DG01MD.f"
	*info = -1;
#line 146 "DG01MD.f"
    } else {
#line 147 "DG01MD.f"
	j = 0;
#line 148 "DG01MD.f"
	if (*n >= 2) {
#line 149 "DG01MD.f"
	    j = *n;
/*           WHILE ( MOD( J, 2 ).EQ.0 ) DO */
#line 151 "DG01MD.f"
L10:
#line 152 "DG01MD.f"
	    if (j % 2 == 0) {
#line 153 "DG01MD.f"
		j /= 2;
#line 154 "DG01MD.f"
		goto L10;
#line 155 "DG01MD.f"
	    }
/*           END WHILE 10 */
#line 157 "DG01MD.f"
	}
#line 158 "DG01MD.f"
	if (j != 1) {
#line 158 "DG01MD.f"
	    *info = -2;
#line 158 "DG01MD.f"
	}
#line 159 "DG01MD.f"
    }

#line 161 "DG01MD.f"
    if (*info != 0) {

/*        Error return. */

#line 165 "DG01MD.f"
	i__1 = -(*info);
#line 165 "DG01MD.f"
	xerbla_("DG01MD", &i__1, (ftnlen)6);
#line 166 "DG01MD.f"
	return 0;
#line 167 "DG01MD.f"
    }

/*     Inplace shuffling of data. */

#line 171 "DG01MD.f"
    j = 1;

#line 173 "DG01MD.f"
    i__1 = *n;
#line 173 "DG01MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 174 "DG01MD.f"
	if (j > i__) {
#line 175 "DG01MD.f"
	    tr = xr[i__];
#line 176 "DG01MD.f"
	    ti = xi[i__];
#line 177 "DG01MD.f"
	    xr[i__] = xr[j];
#line 178 "DG01MD.f"
	    xi[i__] = xi[j];
#line 179 "DG01MD.f"
	    xr[j] = tr;
#line 180 "DG01MD.f"
	    xi[j] = ti;
#line 181 "DG01MD.f"
	}
#line 182 "DG01MD.f"
	k = *n / 2;
/*        REPEAT */
#line 184 "DG01MD.f"
L20:
#line 184 "DG01MD.f"
	if (j > k) {
#line 185 "DG01MD.f"
	    j -= k;
#line 186 "DG01MD.f"
	    k /= 2;
#line 187 "DG01MD.f"
	    if (k >= 2) {
#line 187 "DG01MD.f"
		goto L20;
#line 187 "DG01MD.f"
	    }
#line 188 "DG01MD.f"
	}
/*        UNTIL ( K.LT.2 ) */
#line 190 "DG01MD.f"
	j += k;
#line 191 "DG01MD.f"
/* L30: */
#line 191 "DG01MD.f"
    }

/*     Transform by decimation in time. */

#line 195 "DG01MD.f"
    pi2 = atan(1.) * 8.;
#line 196 "DG01MD.f"
    if (lindi) {
#line 196 "DG01MD.f"
	pi2 = -pi2;
#line 196 "DG01MD.f"
    }

#line 198 "DG01MD.f"
    i__ = 1;

/*     WHILE ( I.LT.N ) DO */

#line 202 "DG01MD.f"
L40:
#line 202 "DG01MD.f"
    if (i__ < *n) {
#line 203 "DG01MD.f"
	l = i__ << 1;
#line 204 "DG01MD.f"
	whelp = pi2 / (doublereal) l;
#line 205 "DG01MD.f"
	wstpi = sin(whelp);
#line 206 "DG01MD.f"
	whelp = sin(whelp * .5);
#line 207 "DG01MD.f"
	wstpr = whelp * -2. * whelp;
#line 208 "DG01MD.f"
	wr = 1.;
#line 209 "DG01MD.f"
	wi = 0.;

#line 211 "DG01MD.f"
	i__1 = i__;
#line 211 "DG01MD.f"
	for (j = 1; j <= i__1; ++j) {

#line 213 "DG01MD.f"
	    i__2 = *n;
#line 213 "DG01MD.f"
	    i__3 = l;
#line 213 "DG01MD.f"
	    for (k = j; i__3 < 0 ? k >= i__2 : k <= i__2; k += i__3) {
#line 214 "DG01MD.f"
		m = k + i__;
#line 215 "DG01MD.f"
		tr = wr * xr[m] - wi * xi[m];
#line 216 "DG01MD.f"
		ti = wr * xi[m] + wi * xr[m];
#line 217 "DG01MD.f"
		xr[m] = xr[k] - tr;
#line 218 "DG01MD.f"
		xi[m] = xi[k] - ti;
#line 219 "DG01MD.f"
		xr[k] += tr;
#line 220 "DG01MD.f"
		xi[k] += ti;
#line 221 "DG01MD.f"
/* L50: */
#line 221 "DG01MD.f"
	    }

#line 223 "DG01MD.f"
	    whelp = wr;
#line 224 "DG01MD.f"
	    wr = wr + wr * wstpr - wi * wstpi;
#line 225 "DG01MD.f"
	    wi = wi + whelp * wstpi + wi * wstpr;
#line 226 "DG01MD.f"
/* L60: */
#line 226 "DG01MD.f"
	}

#line 228 "DG01MD.f"
	i__ = l;
#line 229 "DG01MD.f"
	goto L40;
/*        END WHILE 40 */
#line 231 "DG01MD.f"
    }

#line 233 "DG01MD.f"
    return 0;
/* *** Last line of DG01MD *** */
} /* dg01md_ */

