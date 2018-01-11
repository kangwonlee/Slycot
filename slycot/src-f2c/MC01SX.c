#line 1 "MC01SX.f"
/* MC01SX.f -- translated by f2c (version 20100827).
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

#line 1 "MC01SX.f"
integer mc01sx_(integer *lb, integer *ub, integer *e, doublereal *mant)
{
    /* System generated locals */
    integer ret_val, i__1, i__2, i__3;

    /* Local variables */
    static integer j, mine, maxe;


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

/*     To compute the variation V of the exponents of a series of */
/*     non-zero floating-point numbers: a(j) = MANT(j) * beta**(E(j)), */
/*     where beta is the base of the machine representation of */
/*     floating-point numbers, i.e., */
/*     V = max(E(j)) - min(E(j)), j = LB,...,UB and MANT(j) non-zero. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC01GX by A.J. Geurts. */

/*     REVISIONS */

/*     - */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 54 "MC01SX.f"
    /* Parameter adjustments */
#line 54 "MC01SX.f"
    --mant;
#line 54 "MC01SX.f"
    --e;
#line 54 "MC01SX.f"

#line 54 "MC01SX.f"
    /* Function Body */
#line 54 "MC01SX.f"
    maxe = e[*lb];
#line 55 "MC01SX.f"
    mine = maxe;

#line 57 "MC01SX.f"
    i__1 = *ub;
#line 57 "MC01SX.f"
    for (j = *lb + 1; j <= i__1; ++j) {
#line 58 "MC01SX.f"
	if (mant[j] != 0.) {
/* Computing MAX */
#line 59 "MC01SX.f"
	    i__2 = maxe, i__3 = e[j];
#line 59 "MC01SX.f"
	    maxe = max(i__2,i__3);
/* Computing MIN */
#line 60 "MC01SX.f"
	    i__2 = mine, i__3 = e[j];
#line 60 "MC01SX.f"
	    mine = min(i__2,i__3);
#line 61 "MC01SX.f"
	}
#line 62 "MC01SX.f"
/* L20: */
#line 62 "MC01SX.f"
    }

#line 64 "MC01SX.f"
    ret_val = maxe - mine;

#line 66 "MC01SX.f"
    return ret_val;
/* *** Last line of MC01SX *** */
} /* mc01sx_ */

