#line 1 "MC01SW.f"
/* MC01SW.f -- translated by f2c (version 20100827).
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

#line 1 "MC01SW.f"
/* Subroutine */ int mc01sw_(doublereal *a, integer *b, doublereal *m, 
	integer *e)
{
    static doublereal db;


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

/*     To find the mantissa M and the exponent E of a real number A such */
/*     that */
/*        A = M * B**E */
/*        1 <= ABS( M ) < B */
/*     if A is non-zero. If A is zero, then M and E are set to 0. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     A       (input) DOUBLE PRECISION */
/*             The number whose mantissa and exponent are required. */

/*     B       (input) INTEGER */
/*             The base of the floating-point arithmetic. */

/*     M       (output) DOUBLE PRECISION */
/*             The mantissa of the floating-point representation of A. */

/*     E       (output) INTEGER */
/*             The exponent of the floating-point representation of A. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC01GZ by A.J. Geurts. */

/*     REVISIONS */

/*     - */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Local Scalars .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Quick return if possible. */

#line 74 "MC01SW.f"
    if (*a == 0.) {
#line 75 "MC01SW.f"
	*m = 0.;
#line 76 "MC01SW.f"
	*e = 0;
#line 77 "MC01SW.f"
	return 0;
#line 78 "MC01SW.f"
    }

/*     A non-zero. */

#line 82 "MC01SW.f"
    db = (doublereal) (*b);
#line 83 "MC01SW.f"
    *m = abs(*a);
#line 84 "MC01SW.f"
    *e = 0;
/*     WHILE ( M >= B ) DO */
#line 86 "MC01SW.f"
L20:
#line 86 "MC01SW.f"
    if (*m >= db) {
#line 87 "MC01SW.f"
	*m /= db;
#line 88 "MC01SW.f"
	++(*e);
#line 89 "MC01SW.f"
	goto L20;
#line 90 "MC01SW.f"
    }
/*     END WHILE 20 */
/*     WHILE ( M < 1 ) DO */
#line 93 "MC01SW.f"
L40:
#line 93 "MC01SW.f"
    if (*m < 1.) {
#line 94 "MC01SW.f"
	*m *= db;
#line 95 "MC01SW.f"
	--(*e);
#line 96 "MC01SW.f"
	goto L40;
#line 97 "MC01SW.f"
    }
/*     END WHILE 40 */

#line 100 "MC01SW.f"
    if (*a < 0.) {
#line 100 "MC01SW.f"
	*m = -(*m);
#line 100 "MC01SW.f"
    }

#line 102 "MC01SW.f"
    return 0;
/* *** Last line of MC01SW *** */
} /* mc01sw_ */

