#line 1 "DE01PD.f"
/* DE01PD.f -- translated by f2c (version 20100827).
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

#line 1 "DE01PD.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int de01pd_(char *conv, char *wght, integer *n, doublereal *
	a, doublereal *b, doublereal *w, integer *info, ftnlen conv_len, 
	ftnlen wght_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j, l, m, p1, r1;
    static doublereal t1, t2, t3;
    static integer len;
    extern /* Subroutine */ int dg01od_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, ftnlen, ftnlen), dscal_(integer *, 
	    doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lconv, lwght;
    extern /* Subroutine */ int dladiv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), xerbla_(
	    char *, integer *, ftnlen);


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

/*     To compute the convolution or deconvolution of two real signals */
/*     A and B using the Hartley transform. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     CONV    CHARACTER*1 */
/*             Indicates whether convolution or deconvolution is to be */
/*             performed as follows: */
/*             = 'C':  Convolution; */
/*             = 'D':  Deconvolution. */

/*     WGHT    CHARACTER*1 */
/*             Indicates whether the precomputed weights are available */
/*             or not, as follows: */
/*             = 'A':  available; */
/*             = 'N':  not available. */
/*             Note that if N > 1 and WGHT = 'N' on entry, then WGHT is */
/*             set to 'A' on exit. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of samples.  N must be a power of 2.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the first signal. */
/*             On exit, this array contains the convolution (if */
/*             CONV = 'C') or deconvolution (if CONV = 'D') of the two */
/*             signals. */

/*     B       (input) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the second signal. */
/*             NOTE that this array is overwritten. */

/*     W       (input/output) DOUBLE PRECISION array, */
/*                            dimension (N - LOG2(N)) */
/*             On entry with WGHT = 'A', this array must contain the long */
/*             weight vector computed by a previous call of this routine */
/*             or of the SLICOT Library routine DG01OD.f, with the same */
/*             value of N. If WGHT = 'N', the contents of this array on */
/*             entry is ignored. */
/*             On exit, this array contains the long weight vector. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     This routine computes the convolution or deconvolution of two */
/*     real signals A and B using three scrambled Hartley transforms */
/*     (SLICOT Library routine DG01OD). */

/*     REFERENCES */

/*     [1] Van Loan, Charles. */
/*         Computational frameworks for the fast Fourier transform. */
/*         SIAM, 1992. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires O(N log(N)) floating point operations. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, April 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000. */

/*     KEYWORDS */

/*     Convolution, deconvolution, digital signal processing, */
/*     fast Hartley transform, real signals. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 127 "DE01PD.f"
    /* Parameter adjustments */
#line 127 "DE01PD.f"
    --w;
#line 127 "DE01PD.f"
    --b;
#line 127 "DE01PD.f"
    --a;
#line 127 "DE01PD.f"

#line 127 "DE01PD.f"
    /* Function Body */
#line 127 "DE01PD.f"
    *info = 0;
#line 128 "DE01PD.f"
    lconv = lsame_(conv, "C", (ftnlen)1, (ftnlen)1);
#line 129 "DE01PD.f"
    lwght = lsame_(wght, "A", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 133 "DE01PD.f"
    if (! lconv && ! lsame_(conv, "D", (ftnlen)1, (ftnlen)1)) {
#line 134 "DE01PD.f"
	*info = -1;
#line 135 "DE01PD.f"
    } else if (! lwght && ! lsame_(wght, "N", (ftnlen)1, (ftnlen)1)) {
#line 136 "DE01PD.f"
	*info = -2;
#line 137 "DE01PD.f"
    } else {
#line 138 "DE01PD.f"
	m = 0;
#line 139 "DE01PD.f"
	j = 0;
#line 140 "DE01PD.f"
	if (*n >= 1) {
#line 141 "DE01PD.f"
	    j = *n;
/*           WHILE ( MOD( J, 2 ).EQ.0 ) DO */
#line 143 "DE01PD.f"
L10:
#line 144 "DE01PD.f"
	    if (j % 2 == 0) {
#line 145 "DE01PD.f"
		j /= 2;
#line 146 "DE01PD.f"
		++m;
#line 147 "DE01PD.f"
		goto L10;
#line 148 "DE01PD.f"
	    }
/*           END WHILE 10 */
#line 150 "DE01PD.f"
	    if (j != 1) {
#line 150 "DE01PD.f"
		*info = -3;
#line 150 "DE01PD.f"
	    }
#line 151 "DE01PD.f"
	} else if (*n < 0) {
#line 152 "DE01PD.f"
	    *info = -3;
#line 153 "DE01PD.f"
	}
#line 154 "DE01PD.f"
    }

#line 156 "DE01PD.f"
    if (*info != 0) {

/*        Error return. */

#line 160 "DE01PD.f"
	i__1 = -(*info);
#line 160 "DE01PD.f"
	xerbla_("DE01PD", &i__1, (ftnlen)6);
#line 161 "DE01PD.f"
	return 0;
#line 162 "DE01PD.f"
    }

/*     Quick return if possible. */

#line 166 "DE01PD.f"
    if (*n <= 0) {
#line 167 "DE01PD.f"
	return 0;
#line 168 "DE01PD.f"
    } else if (*n == 1) {
#line 169 "DE01PD.f"
	if (lconv) {
#line 170 "DE01PD.f"
	    a[1] *= b[1];
#line 171 "DE01PD.f"
	} else {
#line 172 "DE01PD.f"
	    a[1] /= b[1];
#line 173 "DE01PD.f"
	}
#line 174 "DE01PD.f"
	return 0;
#line 175 "DE01PD.f"
    }

/*     Scrambled Hartley transforms of A and B. */

#line 179 "DE01PD.f"
    dg01od_("OutputScrambled", wght, n, &a[1], &w[1], info, (ftnlen)15, (
	    ftnlen)1);
#line 180 "DE01PD.f"
    dg01od_("OutputScrambled", wght, n, &b[1], &w[1], info, (ftnlen)15, (
	    ftnlen)1);

/*     Something similar to a Hadamard product/quotient. */

#line 184 "DE01PD.f"
    len = 1;
#line 185 "DE01PD.f"
    if (lconv) {
#line 186 "DE01PD.f"
	a[1] = a[1] * 2. * b[1];
#line 187 "DE01PD.f"
	a[2] = a[2] * 2. * b[2];

#line 189 "DE01PD.f"
	i__1 = m - 1;
#line 189 "DE01PD.f"
	for (l = 1; l <= i__1; ++l) {
#line 190 "DE01PD.f"
	    len <<= 1;
#line 191 "DE01PD.f"
	    r1 = len << 1;

#line 193 "DE01PD.f"
	    i__2 = len + len / 2;
#line 193 "DE01PD.f"
	    for (p1 = len + 1; p1 <= i__2; ++p1) {
#line 194 "DE01PD.f"
		t1 = b[p1] + b[r1];
#line 195 "DE01PD.f"
		t2 = b[p1] - b[r1];
#line 196 "DE01PD.f"
		t3 = t2 * a[p1];
#line 197 "DE01PD.f"
		a[p1] = t1 * a[p1] + t2 * a[r1];
#line 198 "DE01PD.f"
		a[r1] = t1 * a[r1] - t3;
#line 199 "DE01PD.f"
		--r1;
#line 200 "DE01PD.f"
/* L20: */
#line 200 "DE01PD.f"
	    }

#line 202 "DE01PD.f"
/* L30: */
#line 202 "DE01PD.f"
	}

#line 204 "DE01PD.f"
    } else {

#line 206 "DE01PD.f"
	a[1] = a[1] * .5 / b[1];
#line 207 "DE01PD.f"
	a[2] = a[2] * .5 / b[2];

#line 209 "DE01PD.f"
	i__1 = m - 1;
#line 209 "DE01PD.f"
	for (l = 1; l <= i__1; ++l) {
#line 210 "DE01PD.f"
	    len <<= 1;
#line 211 "DE01PD.f"
	    r1 = len << 1;

#line 213 "DE01PD.f"
	    i__2 = len + len / 2;
#line 213 "DE01PD.f"
	    for (p1 = len + 1; p1 <= i__2; ++p1) {
#line 214 "DE01PD.f"
		d__1 = b[p1] + b[r1];
#line 214 "DE01PD.f"
		d__2 = b[r1] - b[p1];
#line 214 "DE01PD.f"
		dladiv_(&a[p1], &a[r1], &d__1, &d__2, &t1, &t2);
#line 216 "DE01PD.f"
		a[p1] = t1;
#line 217 "DE01PD.f"
		a[r1] = t2;
#line 218 "DE01PD.f"
		--r1;
#line 219 "DE01PD.f"
/* L40: */
#line 219 "DE01PD.f"
	    }

#line 221 "DE01PD.f"
/* L50: */
#line 221 "DE01PD.f"
	}

#line 223 "DE01PD.f"
    }

/*     Transposed Hartley transform of A. */

#line 227 "DE01PD.f"
    dg01od_("InputScrambled", wght, n, &a[1], &w[1], info, (ftnlen)14, (
	    ftnlen)1);
#line 228 "DE01PD.f"
    if (lconv) {
#line 229 "DE01PD.f"
	d__1 = .5 / (doublereal) (*n);
#line 229 "DE01PD.f"
	dscal_(n, &d__1, &a[1], &c__1);
#line 230 "DE01PD.f"
    } else {
#line 231 "DE01PD.f"
	d__1 = 2. / (doublereal) (*n);
#line 231 "DE01PD.f"
	dscal_(n, &d__1, &a[1], &c__1);
#line 232 "DE01PD.f"
    }

#line 234 "DE01PD.f"
    return 0;
/* *** Last line of DE01PD *** */
} /* de01pd_ */

