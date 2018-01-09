#line 1 "DE01OD.f"
/* DE01OD.f -- translated by f2c (version 20100827).
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

#line 1 "DE01OD.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int de01od_(char *conv, integer *n, doublereal *a, 
	doublereal *b, integer *info, ftnlen conv_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j;
    static doublereal ac, bc, ci, as;
    static integer kj;
    static doublereal bs, cr, ast;
    static integer nd2p1;
    extern /* Subroutine */ int dg01md_(char *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), dscal_(integer *, doublereal *, 
	    doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lconv;
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
/*     A and B. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     CONV    CHARACTER*1 */
/*             Indicates whether convolution or deconvolution is to be */
/*             performed as follows: */
/*             = 'C':  Convolution; */
/*             = 'D':  Deconvolution. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of samples.  N must be a power of 2.  N >= 2. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the first signal. */
/*             On exit, this array contains the convolution (if */
/*             CONV = 'C') or deconvolution (if CONV = 'D') of the two */
/*             signals. */

/*     B       (input) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the second signal. */
/*             NOTE that this array is overwritten. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     This routine computes the convolution or deconvolution of two real */
/*     signals A and B using an FFT algorithm (SLICOT Library routine */
/*     DG01MD). */

/*     REFERENCES */

/*     [1] Rabiner, L.R. and Rader, C.M. */
/*         Digital Signal Processing. */
/*         IEEE Press, 1972. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0( N*log(N) ) operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */
/*     Supersedes Release 2.0 routine DE01CD by R. Dekeyser, State */
/*     University of Gent, Belgium. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Convolution, deconvolution, digital signal processing, fast */
/*     Fourier transform, real signals. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 112 "DE01OD.f"
    /* Parameter adjustments */
#line 112 "DE01OD.f"
    --b;
#line 112 "DE01OD.f"
    --a;
#line 112 "DE01OD.f"

#line 112 "DE01OD.f"
    /* Function Body */
#line 112 "DE01OD.f"
    *info = 0;
#line 113 "DE01OD.f"
    lconv = lsame_(conv, "C", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 117 "DE01OD.f"
    if (! lconv && ! lsame_(conv, "D", (ftnlen)1, (ftnlen)1)) {
#line 118 "DE01OD.f"
	*info = -1;
#line 119 "DE01OD.f"
    } else {
#line 120 "DE01OD.f"
	j = 0;
#line 121 "DE01OD.f"
	if (*n >= 2) {
#line 122 "DE01OD.f"
	    j = *n;
/*           WHILE ( MOD( J, 2 ).EQ.0 ) DO */
#line 124 "DE01OD.f"
L10:
#line 125 "DE01OD.f"
	    if (j % 2 == 0) {
#line 126 "DE01OD.f"
		j /= 2;
#line 127 "DE01OD.f"
		goto L10;
#line 128 "DE01OD.f"
	    }
/*           END WHILE 10 */
#line 130 "DE01OD.f"
	}
#line 131 "DE01OD.f"
	if (j != 1) {
#line 131 "DE01OD.f"
	    *info = -2;
#line 131 "DE01OD.f"
	}
#line 132 "DE01OD.f"
    }

#line 134 "DE01OD.f"
    if (*info != 0) {

/*        Error return. */

#line 138 "DE01OD.f"
	i__1 = -(*info);
#line 138 "DE01OD.f"
	xerbla_("DE01OD", &i__1, (ftnlen)6);
#line 139 "DE01OD.f"
	return 0;
#line 140 "DE01OD.f"
    }

/*     Fourier transform. */

#line 144 "DE01OD.f"
    dg01md_("Direct", n, &a[1], &b[1], info, (ftnlen)6);

#line 146 "DE01OD.f"
    if (lconv) {
#line 147 "DE01OD.f"
	ast = a[1] * b[1];
#line 148 "DE01OD.f"
    } else {
#line 149 "DE01OD.f"
	if (b[1] == 0.) {
#line 150 "DE01OD.f"
	    ast = 0.;
#line 151 "DE01OD.f"
	} else {
#line 152 "DE01OD.f"
	    ast = a[1] / b[1];
#line 153 "DE01OD.f"
	}
#line 154 "DE01OD.f"
    }

#line 156 "DE01OD.f"
    nd2p1 = *n / 2 + 1;
#line 157 "DE01OD.f"
    j = nd2p1;

#line 159 "DE01OD.f"
    i__1 = *n;
#line 159 "DE01OD.f"
    for (kj = nd2p1; kj <= i__1; ++kj) {

/*        Components of the transform of function A. */

#line 163 "DE01OD.f"
	ac = (a[j] + a[kj]) * .5;
#line 164 "DE01OD.f"
	as = (b[j] - b[kj]) * .5;

/*        Components of the transform of function B. */

#line 168 "DE01OD.f"
	bc = (b[kj] + b[j]) * .5;
#line 169 "DE01OD.f"
	bs = (a[kj] - a[j]) * .5;

/*        Deconvolution by complex division if CONV = 'D'; */
/*        Convolution by complex multiplication if CONV = 'C'. */

#line 174 "DE01OD.f"
	if (lconv) {
#line 175 "DE01OD.f"
	    cr = ac * bc - as * bs;
#line 176 "DE01OD.f"
	    ci = as * bc + ac * bs;
#line 177 "DE01OD.f"
	} else {
/* Computing MAX */
#line 178 "DE01OD.f"
	    d__1 = abs(bc), d__2 = abs(bs);
#line 178 "DE01OD.f"
	    if (max(d__1,d__2) == 0.) {
#line 179 "DE01OD.f"
		cr = 0.;
#line 180 "DE01OD.f"
		ci = 0.;
#line 181 "DE01OD.f"
	    } else {
#line 182 "DE01OD.f"
		dladiv_(&ac, &as, &bc, &bs, &cr, &ci);
#line 183 "DE01OD.f"
	    }
#line 184 "DE01OD.f"
	}

#line 186 "DE01OD.f"
	a[j] = cr;
#line 187 "DE01OD.f"
	b[j] = ci;
#line 188 "DE01OD.f"
	a[kj] = cr;
#line 189 "DE01OD.f"
	b[kj] = -ci;
#line 190 "DE01OD.f"
	--j;
#line 191 "DE01OD.f"
/* L20: */
#line 191 "DE01OD.f"
    }
#line 192 "DE01OD.f"
    a[1] = ast;
#line 193 "DE01OD.f"
    b[1] = 0.;

/*     Inverse Fourier transform. */

#line 197 "DE01OD.f"
    dg01md_("Inverse", n, &a[1], &b[1], info, (ftnlen)7);

#line 199 "DE01OD.f"
    d__1 = 1. / (doublereal) (*n);
#line 199 "DE01OD.f"
    dscal_(n, &d__1, &a[1], &c__1);

#line 201 "DE01OD.f"
    return 0;
/* *** Last line of DE01OD *** */
} /* de01od_ */

