#line 1 "SB02MW.f"
/* SB02MW.f -- translated by f2c (version 20100827).
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

#line 1 "SB02MW.f"
logical sb02mw_(doublereal *reig, doublereal *ieig)
{
    /* System generated locals */
    logical ret_val;

    /* Local variables */
    extern doublereal dlapy2_(doublereal *, doublereal *);


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

/*     To select the stable eigenvalues for solving the discrete-time */
/*     algebraic Riccati equation. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     REIG    (input) DOUBLE PRECISION */
/*             The real part of the current eigenvalue considered. */

/*     IEIG    (input) DOUBLE PRECISION */
/*             The imaginary part of the current eigenvalue considered. */

/*     METHOD */

/*     The function value SB02MW is set to .TRUE. for a stable */
/*     eigenvalue (i.e., with modulus less than one) and to .FALSE., */
/*     otherwise. */

/*     REFERENCES */

/*     None. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Algebraic Riccati equation, closed loop system, discrete-time */
/*     system, optimal regulator, Schur form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. External Functions .. */
/*     .. Executable Statements .. */

#line 75 "SB02MW.f"
    ret_val = dlapy2_(reig, ieig) < 1.;

#line 77 "SB02MW.f"
    return ret_val;
/* *** Last line of SB02MW *** */
} /* sb02mw_ */

