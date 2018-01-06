#line 1 "MC01ND.f"
/* MC01ND.f -- translated by f2c (version 20100827).
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

#line 1 "MC01ND.f"
/* Subroutine */ int mc01nd_(integer *dp, doublereal *xr, doublereal *xi, 
	doublereal *p, doublereal *vr, doublereal *vi, integer *info)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal t;
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

/*     To compute the value of the real polynomial P(x) at a given */
/*     complex point x = x0 using Horner's algorithm. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     DP      (input) INTEGER */
/*             The degree of the polynomial P(x).  DP >= 0. */

/*     XR      (input) DOUBLE PRECISION */
/*     XI      (input) DOUBLE PRECISION */
/*             The real and imaginary parts, respectively, of x0. */

/*     P       (input) DOUBLE PRECISION array, dimension (DP+1) */
/*             This array must contain the coefficients of the polynomial */
/*             P(x) in increasing powers of x. */

/*     VR      (output) DOUBLE PRECISION */
/*     VI      (output) DOUBLE PRECISION */
/*             The real and imaginary parts, respectively, of P(x0). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Given the real polynomial */
/*                                         2                   DP */
/*        P(x) = p(1) + p(2) * x + p(3) * x + ... + p(DP+1) * x  , */

/*     the routine computes the value of P(x0) using the recursion */

/*        q(DP+1) = p(DP+1), */
/*        q(i) = x0*q(i+1) + p(i) for i = DP, DP-1, ..., 1, */

/*     which is known as Horner's algorithm (see [1]). Then q(1) = P(x0). */

/*     REFERENCES */

/*     [1] STOER, J and BULIRSCH, R. */
/*         Introduction to Numerical Analysis. */
/*         Springer-Verlag. 1980. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires DP operations for real arguments and 4*DP */
/*     for complex arguments. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC01BD by Serge Steer. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary polynomial operations, polynomial operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

#line 108 "MC01ND.f"
    /* Parameter adjustments */
#line 108 "MC01ND.f"
    --p;
#line 108 "MC01ND.f"

#line 108 "MC01ND.f"
    /* Function Body */
#line 108 "MC01ND.f"
    if (*dp < 0) {
#line 109 "MC01ND.f"
	*info = -1;

/*        Error return. */

#line 113 "MC01ND.f"
	i__1 = -(*info);
#line 113 "MC01ND.f"
	xerbla_("MC01ND", &i__1, (ftnlen)6);
#line 114 "MC01ND.f"
	return 0;
#line 115 "MC01ND.f"
    }

#line 117 "MC01ND.f"
    *info = 0;
#line 118 "MC01ND.f"
    *vr = p[*dp + 1];
#line 119 "MC01ND.f"
    *vi = 0.;

#line 121 "MC01ND.f"
    if (*dp == 0) {
#line 121 "MC01ND.f"
	return 0;
#line 121 "MC01ND.f"
    }

#line 124 "MC01ND.f"
    if (*xi == 0.) {

/*        X real. */

#line 128 "MC01ND.f"
	for (i__ = *dp; i__ >= 1; --i__) {
#line 129 "MC01ND.f"
	    *vr = *vr * *xr + p[i__];
#line 130 "MC01ND.f"
/* L20: */
#line 130 "MC01ND.f"
	}

#line 132 "MC01ND.f"
    } else {

/*        X complex. */

#line 136 "MC01ND.f"
	for (i__ = *dp; i__ >= 1; --i__) {
#line 137 "MC01ND.f"
	    t = *vr * *xr - *vi * *xi + p[i__];
#line 138 "MC01ND.f"
	    *vi = *vi * *xr + *vr * *xi;
#line 139 "MC01ND.f"
	    *vr = t;
#line 140 "MC01ND.f"
/* L40: */
#line 140 "MC01ND.f"
	}

#line 142 "MC01ND.f"
    }

#line 144 "MC01ND.f"
    return 0;
/* *** Last line of MC01ND *** */
} /* mc01nd_ */

