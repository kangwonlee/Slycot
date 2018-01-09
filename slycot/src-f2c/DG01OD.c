#line 1 "DG01OD.f"
/* DG01OD.f -- translated by f2c (version 20100827).
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

#line 1 "DG01OD.f"
/* Subroutine */ int dg01od_(char *scr, char *wght, integer *n, doublereal *a,
	 doublereal *w, integer *info, ftnlen scr_len, ftnlen wght_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double atan(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j, l, m, p1, p2, q1, q2, r1, r2, s1, s2;
    static doublereal t1, t2, cf, sf, th;
    static integer len;
    static logical lfwd, lscr;
    static integer wpos;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lwght;
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

/*     To compute the (scrambled) discrete Hartley transform of */
/*     a real signal. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SCR     CHARACTER*1 */
/*             Indicates whether the signal is scrambled on input or */
/*             on output as follows: */
/*             = 'N':  the signal is not scrambled at all; */
/*             = 'I':  the input signal is bit-reversed; */
/*             = 'O':  the output transform is bit-reversed. */

/*     WGHT    CHARACTER*1 */
/*             Indicates whether the precomputed weights are available */
/*             or not, as follows: */
/*             = 'A':  available; */
/*             = 'N':  not available. */
/*             Note that if N > 1 and WGHT = 'N' on entry, then WGHT is */
/*             set to 'A' on exit. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             Number of real samples. N must be a power of 2. */
/*             N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry with SCR = 'N' or SCR = 'O', this array must */
/*             contain the input signal. */
/*             On entry with SCR = 'I', this array must contain the */
/*             bit-reversed input signal. */
/*             On exit with SCR = 'N' or SCR = 'I', this array contains */
/*             the Hartley transform of the input signal. */
/*             On exit with SCR = 'O', this array contains the */
/*             bit-reversed Hartley transform. */

/*     W       (input/output) DOUBLE PRECISION array, */
/*                            dimension (N - LOG2(N)) */
/*             On entry with WGHT = 'A', this array must contain the long */
/*             weight vector computed by a previous call of this routine */
/*             with the same value of N. If WGHT = 'N', the contents of */
/*             this array on entry is ignored. */
/*             On exit, this array contains the long weight vector. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     This routine uses a Hartley butterfly algorithm as described */
/*     in [1]. */

/*     REFERENCES */

/*     [1] Van Loan, Charles. */
/*         Computational frameworks for the fast Fourier transform. */
/*         SIAM, 1992. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable and requires O(N log(N)) */
/*     floating point operations. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, April 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000. */

/*     KEYWORDS */

/*     Digital signal processing, fast Hartley transform, real signals. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 128 "DG01OD.f"
    /* Parameter adjustments */
#line 128 "DG01OD.f"
    --w;
#line 128 "DG01OD.f"
    --a;
#line 128 "DG01OD.f"

#line 128 "DG01OD.f"
    /* Function Body */
#line 128 "DG01OD.f"
    *info = 0;
#line 129 "DG01OD.f"
    lfwd = lsame_(scr, "N", (ftnlen)1, (ftnlen)1) || lsame_(scr, "I", (ftnlen)
	    1, (ftnlen)1);
#line 130 "DG01OD.f"
    lscr = lsame_(scr, "I", (ftnlen)1, (ftnlen)1) || lsame_(scr, "O", (ftnlen)
	    1, (ftnlen)1);
#line 131 "DG01OD.f"
    lwght = lsame_(wght, "A", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 135 "DG01OD.f"
    if (! (lfwd || lscr)) {
#line 136 "DG01OD.f"
	*info = -1;
#line 137 "DG01OD.f"
    } else if (! lwght && ! lsame_(wght, "N", (ftnlen)1, (ftnlen)1)) {
#line 138 "DG01OD.f"
	*info = -2;
#line 139 "DG01OD.f"
    } else {
#line 140 "DG01OD.f"
	m = 0;
#line 141 "DG01OD.f"
	j = 0;
#line 142 "DG01OD.f"
	if (*n >= 1) {
#line 143 "DG01OD.f"
	    j = *n;
/*           WHILE ( MOD( J, 2 ).EQ.0 ) DO */
#line 145 "DG01OD.f"
L10:
#line 146 "DG01OD.f"
	    if (j % 2 == 0) {
#line 147 "DG01OD.f"
		j /= 2;
#line 148 "DG01OD.f"
		++m;
#line 149 "DG01OD.f"
		goto L10;
#line 150 "DG01OD.f"
	    }
/*           END WHILE 10 */
#line 152 "DG01OD.f"
	    if (j != 1) {
#line 152 "DG01OD.f"
		*info = -3;
#line 152 "DG01OD.f"
	    }
#line 153 "DG01OD.f"
	} else if (*n < 0) {
#line 154 "DG01OD.f"
	    *info = -3;
#line 155 "DG01OD.f"
	}
#line 156 "DG01OD.f"
    }

#line 158 "DG01OD.f"
    if (*info != 0) {

/*        Error return. */

#line 162 "DG01OD.f"
	i__1 = -(*info);
#line 162 "DG01OD.f"
	xerbla_("DG01OD", &i__1, (ftnlen)6);
#line 163 "DG01OD.f"
	return 0;
#line 164 "DG01OD.f"
    }

/*     Quick return if possible. */

#line 168 "DG01OD.f"
    if (*n <= 1) {
#line 168 "DG01OD.f"
	return 0;
#line 168 "DG01OD.f"
    }

#line 171 "DG01OD.f"
    if (! lwght) {

/*        Compute the long weight vector via subvector scaling. */

#line 175 "DG01OD.f"
	r1 = 1;
#line 176 "DG01OD.f"
	len = 1;
#line 177 "DG01OD.f"
	th = atan(1.) * 4. / (doublereal) (*n);

#line 179 "DG01OD.f"
	i__1 = m - 2;
#line 179 "DG01OD.f"
	for (l = 1; l <= i__1; ++l) {
#line 180 "DG01OD.f"
	    len <<= 1;
#line 181 "DG01OD.f"
	    th *= 2.;
#line 182 "DG01OD.f"
	    cf = cos(th);
#line 183 "DG01OD.f"
	    sf = sin(th);
#line 184 "DG01OD.f"
	    w[r1] = cf;
#line 185 "DG01OD.f"
	    w[r1 + 1] = sf;
#line 186 "DG01OD.f"
	    r1 += 2;

#line 188 "DG01OD.f"
	    i__2 = len - 2;
#line 188 "DG01OD.f"
	    for (i__ = 1; i__ <= i__2; i__ += 2) {
#line 189 "DG01OD.f"
		w[r1] = cf * w[i__] - sf * w[i__ + 1];
#line 190 "DG01OD.f"
		w[r1 + 1] = sf * w[i__] + cf * w[i__ + 1];
#line 191 "DG01OD.f"
		r1 += 2;
#line 192 "DG01OD.f"
/* L20: */
#line 192 "DG01OD.f"
	    }

#line 194 "DG01OD.f"
/* L30: */
#line 194 "DG01OD.f"
	}

#line 196 "DG01OD.f"
	p1 = 3;
#line 197 "DG01OD.f"
	q1 = r1 - 2;

#line 199 "DG01OD.f"
	for (l = m - 2; l >= 1; --l) {

#line 201 "DG01OD.f"
	    i__1 = q1;
#line 201 "DG01OD.f"
	    for (i__ = p1; i__ <= i__1; i__ += 4) {
#line 202 "DG01OD.f"
		w[r1] = w[i__];
#line 203 "DG01OD.f"
		w[r1 + 1] = w[i__ + 1];
#line 204 "DG01OD.f"
		r1 += 2;
#line 205 "DG01OD.f"
/* L40: */
#line 205 "DG01OD.f"
	    }

#line 207 "DG01OD.f"
	    p1 = q1 + 4;
#line 208 "DG01OD.f"
	    q1 = r1 - 2;
#line 209 "DG01OD.f"
/* L50: */
#line 209 "DG01OD.f"
	}

#line 211 "DG01OD.f"
	*(unsigned char *)wght = 'A';

#line 213 "DG01OD.f"
    }

#line 215 "DG01OD.f"
    if (lfwd && ! lscr) {

/*        Inplace shuffling of data. */

#line 219 "DG01OD.f"
	j = 1;

#line 221 "DG01OD.f"
	i__1 = *n;
#line 221 "DG01OD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 222 "DG01OD.f"
	    if (j > i__) {
#line 223 "DG01OD.f"
		t1 = a[i__];
#line 224 "DG01OD.f"
		a[i__] = a[j];
#line 225 "DG01OD.f"
		a[j] = t1;
#line 226 "DG01OD.f"
	    }
#line 227 "DG01OD.f"
	    l = *n / 2;
/*           REPEAT */
#line 229 "DG01OD.f"
L60:
#line 229 "DG01OD.f"
	    if (j > l) {
#line 230 "DG01OD.f"
		j -= l;
#line 231 "DG01OD.f"
		l /= 2;
#line 232 "DG01OD.f"
		if (l >= 2) {
#line 232 "DG01OD.f"
		    goto L60;
#line 232 "DG01OD.f"
		}
#line 233 "DG01OD.f"
	    }
/*           UNTIL ( L.LT.2 ) */
#line 235 "DG01OD.f"
	    j += l;
#line 236 "DG01OD.f"
/* L70: */
#line 236 "DG01OD.f"
	}

#line 238 "DG01OD.f"
    }

#line 240 "DG01OD.f"
    if (lfwd) {

/*        Compute Hartley transform with butterfly operators. */

#line 244 "DG01OD.f"
	i__1 = *n;
#line 244 "DG01OD.f"
	for (j = 2; j <= i__1; j += 2) {
#line 245 "DG01OD.f"
	    t1 = a[j];
#line 246 "DG01OD.f"
	    a[j] = a[j - 1] - t1;
#line 247 "DG01OD.f"
	    a[j - 1] += t1;
#line 248 "DG01OD.f"
/* L110: */
#line 248 "DG01OD.f"
	}

#line 250 "DG01OD.f"
	len = 1;
#line 251 "DG01OD.f"
	wpos = *n - (m << 1) + 1;

#line 253 "DG01OD.f"
	i__1 = m - 1;
#line 253 "DG01OD.f"
	for (l = 1; l <= i__1; ++l) {
#line 254 "DG01OD.f"
	    len <<= 1;
#line 255 "DG01OD.f"
	    p2 = 1;
#line 256 "DG01OD.f"
	    q2 = len + 1;
#line 257 "DG01OD.f"
	    r2 = len / 2 + 1;
#line 258 "DG01OD.f"
	    s2 = r2 + q2 - 1;

#line 260 "DG01OD.f"
	    i__2 = *n / (len << 1) - 1;
#line 260 "DG01OD.f"
	    for (i__ = 0; i__ <= i__2; ++i__) {
#line 261 "DG01OD.f"
		t1 = a[q2];
#line 262 "DG01OD.f"
		a[q2] = a[p2] - t1;
#line 263 "DG01OD.f"
		a[p2] += t1;
#line 264 "DG01OD.f"
		t1 = a[s2];
#line 265 "DG01OD.f"
		a[s2] = a[r2] - t1;
#line 266 "DG01OD.f"
		a[r2] += t1;

#line 268 "DG01OD.f"
		p1 = p2 + 1;
#line 269 "DG01OD.f"
		q1 = p1 + len;
#line 270 "DG01OD.f"
		r1 = q1 - 2;
#line 271 "DG01OD.f"
		s1 = r1 + len;

#line 273 "DG01OD.f"
		i__3 = wpos + len - 3;
#line 273 "DG01OD.f"
		for (j = wpos; j <= i__3; j += 2) {
#line 274 "DG01OD.f"
		    cf = w[j];
#line 275 "DG01OD.f"
		    sf = w[j + 1];
#line 276 "DG01OD.f"
		    t1 = cf * a[q1] + sf * a[s1];
#line 277 "DG01OD.f"
		    t2 = -cf * a[s1] + sf * a[q1];
#line 278 "DG01OD.f"
		    a[q1] = a[p1] - t1;
#line 279 "DG01OD.f"
		    a[p1] += t1;
#line 280 "DG01OD.f"
		    a[s1] = a[r1] - t2;
#line 281 "DG01OD.f"
		    a[r1] += t2;
#line 282 "DG01OD.f"
		    ++p1;
#line 283 "DG01OD.f"
		    ++q1;
#line 284 "DG01OD.f"
		    --r1;
#line 285 "DG01OD.f"
		    --s1;
#line 286 "DG01OD.f"
/* L120: */
#line 286 "DG01OD.f"
		}

#line 288 "DG01OD.f"
		p2 += len << 1;
#line 289 "DG01OD.f"
		q2 += len << 1;
#line 290 "DG01OD.f"
		r2 += len << 1;
#line 291 "DG01OD.f"
		s2 += len << 1;
#line 292 "DG01OD.f"
/* L130: */
#line 292 "DG01OD.f"
	    }

#line 294 "DG01OD.f"
	    wpos = wpos - (len << 1) + 2;
#line 295 "DG01OD.f"
/* L140: */
#line 295 "DG01OD.f"
	}

#line 297 "DG01OD.f"
    } else {

/*        Compute Hartley transform with transposed butterfly operators. */

#line 301 "DG01OD.f"
	wpos = 1;
#line 302 "DG01OD.f"
	len = *n;

#line 304 "DG01OD.f"
	for (l = m - 1; l >= 1; --l) {
#line 305 "DG01OD.f"
	    len /= 2;
#line 306 "DG01OD.f"
	    p2 = 1;
#line 307 "DG01OD.f"
	    q2 = len + 1;
#line 308 "DG01OD.f"
	    r2 = len / 2 + 1;
#line 309 "DG01OD.f"
	    s2 = r2 + q2 - 1;

#line 311 "DG01OD.f"
	    i__1 = *n / (len << 1) - 1;
#line 311 "DG01OD.f"
	    for (i__ = 0; i__ <= i__1; ++i__) {
#line 312 "DG01OD.f"
		t1 = a[q2];
#line 313 "DG01OD.f"
		a[q2] = a[p2] - t1;
#line 314 "DG01OD.f"
		a[p2] += t1;
#line 315 "DG01OD.f"
		t1 = a[s2];
#line 316 "DG01OD.f"
		a[s2] = a[r2] - t1;
#line 317 "DG01OD.f"
		a[r2] += t1;

#line 319 "DG01OD.f"
		p1 = p2 + 1;
#line 320 "DG01OD.f"
		q1 = p1 + len;
#line 321 "DG01OD.f"
		r1 = q1 - 2;
#line 322 "DG01OD.f"
		s1 = r1 + len;

#line 324 "DG01OD.f"
		i__2 = wpos + len - 3;
#line 324 "DG01OD.f"
		for (j = wpos; j <= i__2; j += 2) {
#line 325 "DG01OD.f"
		    cf = w[j];
#line 326 "DG01OD.f"
		    sf = w[j + 1];
#line 327 "DG01OD.f"
		    t1 = a[p1] - a[q1];
#line 328 "DG01OD.f"
		    t2 = a[r1] - a[s1];
#line 329 "DG01OD.f"
		    a[p1] += a[q1];
#line 330 "DG01OD.f"
		    a[r1] += a[s1];
#line 331 "DG01OD.f"
		    a[q1] = cf * t1 + sf * t2;
#line 332 "DG01OD.f"
		    a[s1] = -cf * t2 + sf * t1;
#line 333 "DG01OD.f"
		    ++p1;
#line 334 "DG01OD.f"
		    ++q1;
#line 335 "DG01OD.f"
		    --r1;
#line 336 "DG01OD.f"
		    --s1;
#line 337 "DG01OD.f"
/* L210: */
#line 337 "DG01OD.f"
		}

#line 339 "DG01OD.f"
		p2 += len << 1;
#line 340 "DG01OD.f"
		q2 += len << 1;
#line 341 "DG01OD.f"
		r2 += len << 1;
#line 342 "DG01OD.f"
		s2 += len << 1;
#line 343 "DG01OD.f"
/* L220: */
#line 343 "DG01OD.f"
	    }

#line 345 "DG01OD.f"
	    wpos = wpos + len - 2;
#line 346 "DG01OD.f"
/* L230: */
#line 346 "DG01OD.f"
	}

#line 348 "DG01OD.f"
	i__1 = *n;
#line 348 "DG01OD.f"
	for (j = 2; j <= i__1; j += 2) {
#line 349 "DG01OD.f"
	    t1 = a[j];
#line 350 "DG01OD.f"
	    a[j] = a[j - 1] - t1;
#line 351 "DG01OD.f"
	    a[j - 1] += t1;
#line 352 "DG01OD.f"
/* L240: */
#line 352 "DG01OD.f"
	}

#line 354 "DG01OD.f"
    }
#line 355 "DG01OD.f"
    return 0;
/* *** Last line of DG01OD *** */
} /* dg01od_ */

