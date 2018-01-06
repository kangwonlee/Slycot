#line 1 "MC01SY.f"
/* MC01SY.f -- translated by f2c (version 20100827).
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

#line 1 "MC01SY.f"
/* Subroutine */ int mc01sy_(doublereal *m, integer *e, integer *b, 
	doublereal *a, logical *ovflow)
{
    static integer et;
    static doublereal mt, base;
    static integer emin, emax, expon;
    extern doublereal dlamch_(char *, ftnlen);


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

/*     To find a real number A from its mantissa M and its exponent E, */
/*     i.e., */
/*        A = M * B**E. */
/*     M and E need not be the standard floating-point values. */
/*     If ABS(A) < B**(EMIN-1), i.e. the smallest positive model number, */
/*     then the routine returns A = 0. */
/*     If M = 0, then the routine returns A = 0 regardless of the value */
/*     of E. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) DOUBLE PRECISION */
/*             The mantissa of the floating-point representation of A. */

/*     E       (input) INTEGER */
/*             The exponent of the floating-point representation of A. */

/*     B       (input) INTEGER */
/*             The base of the floating-point arithmetic. */

/*     A       (output) DOUBLE PRECISION */
/*             The value of M * B**E. */

/*     OVFLOW  (output) LOGICAL */
/*             The value .TRUE., if ABS(M) * B**E >= B**EMAX (where EMAX */
/*             is the largest possible exponent) and .FALSE. otherwise. */
/*             A is not defined if OVFLOW = .TRUE.. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC01GY by A.J. Geurts. */

/*     REVISIONS */

/*     - */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 85 "MC01SY.f"
    *ovflow = FALSE_;

#line 87 "MC01SY.f"
    if (*m == 0. || *e == 0) {
#line 88 "MC01SY.f"
	*a = *m;
#line 89 "MC01SY.f"
	return 0;
#line 90 "MC01SY.f"
    }

/*     Determination of the mantissa MT and the exponent ET of the */
/*     standard floating-point representation. */

#line 95 "MC01SY.f"
    emin = (integer) dlamch_("Minimum exponent", (ftnlen)16);
#line 96 "MC01SY.f"
    emax = (integer) dlamch_("Largest exponent", (ftnlen)16);
#line 97 "MC01SY.f"
    mt = *m;
#line 98 "MC01SY.f"
    et = *e;
/*     WHILE ( ABS( MT ) >= B ) DO */
#line 100 "MC01SY.f"
L20:
#line 100 "MC01SY.f"
    if (abs(mt) >= (doublereal) (*b)) {
#line 101 "MC01SY.f"
	mt /= *b;
#line 102 "MC01SY.f"
	++et;
#line 103 "MC01SY.f"
	goto L20;
#line 104 "MC01SY.f"
    }
/*     END WHILE 20 */
/*     WHILE ( ABS( MT ) < 1 ) DO */
#line 107 "MC01SY.f"
L40:
#line 107 "MC01SY.f"
    if (abs(mt) < 1.) {
#line 108 "MC01SY.f"
	mt *= *b;
#line 109 "MC01SY.f"
	--et;
#line 110 "MC01SY.f"
	goto L40;
#line 111 "MC01SY.f"
    }
/*     END WHILE 40 */

#line 114 "MC01SY.f"
    if (et < emin) {
#line 115 "MC01SY.f"
	*a = 0.;
#line 116 "MC01SY.f"
	return 0;
#line 117 "MC01SY.f"
    }

#line 119 "MC01SY.f"
    if (et >= emax) {
#line 120 "MC01SY.f"
	*ovflow = TRUE_;
#line 121 "MC01SY.f"
	return 0;
#line 122 "MC01SY.f"
    }

/*     Computation of the value of A by the relation */
/*     M * B**E = A * (BASE)**EXPON */

#line 127 "MC01SY.f"
    expon = abs(et);
#line 128 "MC01SY.f"
    *a = mt;
#line 129 "MC01SY.f"
    base = (doublereal) (*b);
#line 130 "MC01SY.f"
    if (et < 0) {
#line 130 "MC01SY.f"
	base = 1. / base;
#line 130 "MC01SY.f"
    }
/*     WHILE ( not EXPON = 0 ) DO */
#line 132 "MC01SY.f"
L60:
#line 132 "MC01SY.f"
    if (expon != 0) {
#line 133 "MC01SY.f"
	if (expon % 2 == 0) {
#line 134 "MC01SY.f"
	    base *= base;
#line 135 "MC01SY.f"
	    expon /= 2;
#line 136 "MC01SY.f"
	} else {
#line 137 "MC01SY.f"
	    *a *= base;
#line 138 "MC01SY.f"
	    --expon;
#line 139 "MC01SY.f"
	}
#line 140 "MC01SY.f"
	goto L60;
#line 141 "MC01SY.f"
    }
/*     END WHILE 60 */

#line 144 "MC01SY.f"
    return 0;
/* *** Last line of MC01SY *** */
} /* mc01sy_ */

