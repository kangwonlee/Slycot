#line 1 "DK01MD.f"
/* DK01MD.f -- translated by f2c (version 20100827).
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

#line 1 "DK01MD.f"
/* Subroutine */ int dk01md_(char *type__, integer *n, doublereal *a, integer 
	*info, ftnlen type_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double atan(doublereal), cos(doublereal);

    /* Local variables */
    static integer i__, n1;
    static doublereal fn, buf, temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical mtype, ntype;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical mntype;


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

/*     To apply an anti-aliasing window to a real signal. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYPE    CHARACTER*1 */
/*             Indicates the type of window to be applied to the signal */
/*             as follows: */
/*             = 'M':  Hamming window; */
/*             = 'N':  Hann window; */
/*             = 'Q':  Quadratic window. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of samples.  N >= 1. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the signal to be */
/*             processed. */
/*             On exit, this array contains the windowing function. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     If TYPE = 'M', then a Hamming window is applied to A(1),...,A(N), */
/*     which yields */
/*       _ */
/*       A(i) = (0.54 + 0.46*cos(pi*(i-1)/(N-1)))*A(i), i = 1,2,...,N. */

/*     If TYPE = 'N', then a Hann window is applied to A(1),...,A(N), */
/*     which yields */
/*       _ */
/*       A(i) = 0.5*(1 + cos(pi*(i-1)/(N-1)))*A(i), i = 1,2,...,N. */

/*     If TYPE = 'Q', then a quadratic window is applied to A(1),..., */
/*     A(N), which yields */
/*       _ */
/*       A(i) = (1 - 2*((i-1)/(N-1))**2)*(1 - (i-1)/(N-1))*A(i), */
/*                                             i = 1,2,...,(N-1)/2+1; */
/*       _ */
/*       A(i) = 2*(1 - ((i-1)/(N-1))**3)*A(i), i = (N-1)/2+2,...,N. */

/*     REFERENCES */

/*     [1] Rabiner, L.R. and Rader, C.M. */
/*         Digital Signal Processing. */
/*         IEEE Press, 1972. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0( N ) operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */
/*     Supersedes Release 2.0 routine DK01AD by R. Dekeyser, State */
/*     University of Gent, Belgium. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Digital signal processing, Hamming window, Hann window, real */
/*     signals, windowing. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 122 "DK01MD.f"
    /* Parameter adjustments */
#line 122 "DK01MD.f"
    --a;
#line 122 "DK01MD.f"

#line 122 "DK01MD.f"
    /* Function Body */
#line 122 "DK01MD.f"
    *info = 0;
#line 123 "DK01MD.f"
    mtype = lsame_(type__, "M", (ftnlen)1, (ftnlen)1);
#line 124 "DK01MD.f"
    ntype = lsame_(type__, "N", (ftnlen)1, (ftnlen)1);
#line 125 "DK01MD.f"
    mntype = mtype || ntype;

/*     Test the input scalar arguments. */

#line 129 "DK01MD.f"
    if (! mntype && ! lsame_(type__, "Q", (ftnlen)1, (ftnlen)1)) {
#line 131 "DK01MD.f"
	*info = -1;
#line 132 "DK01MD.f"
    } else if (*n <= 0) {
#line 133 "DK01MD.f"
	*info = -2;
#line 134 "DK01MD.f"
    }

#line 136 "DK01MD.f"
    if (*info != 0) {

/*        Error return. */

#line 140 "DK01MD.f"
	i__1 = -(*info);
#line 140 "DK01MD.f"
	xerbla_("DK01MD", &i__1, (ftnlen)6);
#line 141 "DK01MD.f"
	return 0;
#line 142 "DK01MD.f"
    }

#line 144 "DK01MD.f"
    fn = (doublereal) (*n - 1);
#line 145 "DK01MD.f"
    if (mntype) {
#line 145 "DK01MD.f"
	temp = atan(1.) * 4. / fn;
#line 145 "DK01MD.f"
    }

#line 147 "DK01MD.f"
    if (mtype) {

/*        Hamming window. */

#line 151 "DK01MD.f"
	i__1 = *n;
#line 151 "DK01MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 152 "DK01MD.f"
	    a[i__] *= cos(temp * (doublereal) (i__ - 1)) * .46 + .54;
#line 153 "DK01MD.f"
/* L10: */
#line 153 "DK01MD.f"
	}

#line 155 "DK01MD.f"
    } else if (ntype) {

/*        Hann window. */

#line 159 "DK01MD.f"
	i__1 = *n;
#line 159 "DK01MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 160 "DK01MD.f"
	    a[i__] = a[i__] * .5 * (cos(temp * (doublereal) (i__ - 1)) + 1.);
#line 161 "DK01MD.f"
/* L20: */
#line 161 "DK01MD.f"
	}

#line 163 "DK01MD.f"
    } else {

/*        Quadratic window. */

#line 167 "DK01MD.f"
	n1 = (*n - 1) / 2 + 1;

#line 169 "DK01MD.f"
	i__1 = *n;
#line 169 "DK01MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 170 "DK01MD.f"
	    buf = (doublereal) (i__ - 1) / fn;
/* Computing 2nd power */
#line 171 "DK01MD.f"
	    d__1 = buf;
#line 171 "DK01MD.f"
	    temp = d__1 * d__1;
#line 172 "DK01MD.f"
	    if (i__ <= n1) {
#line 173 "DK01MD.f"
		a[i__] = a[i__] * (1. - temp * 2.) * (1. - buf);
#line 174 "DK01MD.f"
	    } else {
#line 175 "DK01MD.f"
		a[i__] = a[i__] * 2. * (1. - buf * temp);
#line 176 "DK01MD.f"
	    }
#line 177 "DK01MD.f"
/* L30: */
#line 177 "DK01MD.f"
	}

#line 179 "DK01MD.f"
    }

#line 181 "DK01MD.f"
    return 0;
/* *** Last line of DK01MD *** */
} /* dk01md_ */

