#line 1 "DG01NY.f"
/* DG01NY.f -- translated by f2c (version 20100827).
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

#line 1 "DG01NY.f"
/* Subroutine */ int dg01ny_(char *indi, integer *n, doublereal *xr, 
	doublereal *xi, ftnlen indi_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double atan(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j, n2;
    static doublereal ai, bi, ar, br, wi, wr, pi2;
    static logical lindi;
    static doublereal helpi;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal helpr, whelp, wstpi, wstpr;


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

/*     For efficiency, no tests of the input scalar parameters are */
/*     performed. */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 45 "DG01NY.f"
    /* Parameter adjustments */
#line 45 "DG01NY.f"
    --xi;
#line 45 "DG01NY.f"
    --xr;
#line 45 "DG01NY.f"

#line 45 "DG01NY.f"
    /* Function Body */
#line 45 "DG01NY.f"
    lindi = lsame_(indi, "D", (ftnlen)1, (ftnlen)1);

/*     Initialisation. */

#line 49 "DG01NY.f"
    pi2 = atan(1.) * 8.;
#line 50 "DG01NY.f"
    if (lindi) {
#line 50 "DG01NY.f"
	pi2 = -pi2;
#line 50 "DG01NY.f"
    }

#line 52 "DG01NY.f"
    whelp = pi2 / (doublereal) (*n << 1);
#line 53 "DG01NY.f"
    wstpi = sin(whelp);
#line 54 "DG01NY.f"
    whelp = sin(whelp * .5);
#line 55 "DG01NY.f"
    wstpr = whelp * -2. * whelp;
#line 56 "DG01NY.f"
    wi = 0.;

#line 58 "DG01NY.f"
    if (lindi) {
#line 59 "DG01NY.f"
	wr = 1.;
#line 60 "DG01NY.f"
	xr[*n + 1] = xr[1];
#line 61 "DG01NY.f"
	xi[*n + 1] = xi[1];
#line 62 "DG01NY.f"
    } else {
#line 63 "DG01NY.f"
	wr = -1.;
#line 64 "DG01NY.f"
    }

/*     Recursion. */

#line 68 "DG01NY.f"
    n2 = *n / 2 + 1;
#line 69 "DG01NY.f"
    i__1 = n2;
#line 69 "DG01NY.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 70 "DG01NY.f"
	j = *n + 2 - i__;
#line 71 "DG01NY.f"
	ar = xr[i__] + xr[j];
#line 72 "DG01NY.f"
	ai = xi[i__] - xi[j];
#line 73 "DG01NY.f"
	br = xi[i__] + xi[j];
#line 74 "DG01NY.f"
	bi = xr[j] - xr[i__];
#line 75 "DG01NY.f"
	if (lindi) {
#line 76 "DG01NY.f"
	    ar *= .5;
#line 77 "DG01NY.f"
	    ai *= .5;
#line 78 "DG01NY.f"
	    br *= .5;
#line 79 "DG01NY.f"
	    bi *= .5;
#line 80 "DG01NY.f"
	}
#line 81 "DG01NY.f"
	helpr = wr * br - wi * bi;
#line 82 "DG01NY.f"
	helpi = wr * bi + wi * br;
#line 83 "DG01NY.f"
	xr[i__] = ar + helpr;
#line 84 "DG01NY.f"
	xi[i__] = ai + helpi;
#line 85 "DG01NY.f"
	xr[j] = ar - helpr;
#line 86 "DG01NY.f"
	xi[j] = helpi - ai;
#line 87 "DG01NY.f"
	whelp = wr;
#line 88 "DG01NY.f"
	wr = wr + wr * wstpr - wi * wstpi;
#line 89 "DG01NY.f"
	wi = wi + wi * wstpr + whelp * wstpi;
#line 90 "DG01NY.f"
/* L10: */
#line 90 "DG01NY.f"
    }

#line 92 "DG01NY.f"
    return 0;
/* *** Last line of DG01NY *** */
} /* dg01ny_ */

