#line 1 "SB03OV.f"
/* SB03OV.f -- translated by f2c (version 20100827).
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

#line 1 "SB03OV.f"
/* Subroutine */ int sb03ov_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *s)
{
    static doublereal d__;
    extern doublereal dlapy3_(doublereal *, doublereal *, doublereal *);


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

/*     To construct a complex plane rotation such that, for a complex */
/*     number  a  and a real number  b, */

/*        ( conjg( c )  s )*( a ) = ( d ), */
/*        (       -s    c ) ( b )   ( 0 ) */

/*     where  d  is always real and is overwritten on  a,  so that on */
/*     return the imaginary part of  a  is zero.  b  is unaltered. */

/*     This routine has A and C declared as REAL, because it is intended */
/*     for use within a real Lyapunov solver and the REAL declarations */
/*     mean that a standard Fortran DOUBLE PRECISION version may be */
/*     readily constructed.  However A and C could safely be declared */
/*     COMPLEX in the calling program, although some systems may give a */
/*     type mismatch warning. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     A       (input/output) DOUBLE PRECISION array, dimension (2) */
/*             On entry, A(1) and A(2) must contain the real and */
/*             imaginary part, respectively, of the complex number a. */
/*             On exit, A(1) contains the real part of d, and A(2) is */
/*             set to zero. */

/*     B       (input) DOUBLE PRECISION */
/*             The real number b. */

/*     C       (output) DOUBLE PRECISION array, dimension (2) */
/*             C(1) and C(2) contain the real and imaginary part, */
/*             respectively, of the complex number c, the cosines of */
/*             the plane rotation. */

/*     S       (output) DOUBLE PRECISION */
/*             The real number s, the sines of the plane rotation. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB03CV by Sven Hammarling, */
/*     NAG Ltd., United Kingdom, May 1985. */

/*     REVISIONS */

/*     Dec. 1997. */

/*     KEYWORDS */

/*     Lyapunov equation, orthogonal transformation. */

/*     ***************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. Executable Statements .. */

#line 90 "SB03OV.f"
    /* Parameter adjustments */
#line 90 "SB03OV.f"
    --c__;
#line 90 "SB03OV.f"
    --a;
#line 90 "SB03OV.f"

#line 90 "SB03OV.f"
    /* Function Body */
#line 90 "SB03OV.f"
    d__ = dlapy3_(&a[1], &a[2], b);
#line 91 "SB03OV.f"
    if (d__ == 0.) {
#line 92 "SB03OV.f"
	c__[1] = 1.;
#line 93 "SB03OV.f"
	c__[2] = 0.;
#line 94 "SB03OV.f"
	*s = 0.;
#line 95 "SB03OV.f"
    } else {
#line 96 "SB03OV.f"
	c__[1] = a[1] / d__;
#line 97 "SB03OV.f"
	c__[2] = a[2] / d__;
#line 98 "SB03OV.f"
	*s = *b / d__;
#line 99 "SB03OV.f"
	a[1] = d__;
#line 100 "SB03OV.f"
	a[2] = 0.;
#line 101 "SB03OV.f"
    }

#line 103 "SB03OV.f"
    return 0;
/* *** Last line of SB03OV *** */
} /* sb03ov_ */

