#line 1 "SG03BY.f"
/* SG03BY.f -- translated by f2c (version 20100827).
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

#line 1 "SG03BY.f"
/* Subroutine */ int sg03by_(doublereal *xr, doublereal *xi, doublereal *yr, 
	doublereal *yi, doublereal *cr, doublereal *ci, doublereal *sr, 
	doublereal *si, doublereal *z__)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);


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

/*     To compute the parameters for the complex Givens rotation */

/*        (  CR-CI*I   SR-SI*I )   ( XR+XI*I )   ( Z ) */
/*        (                    ) * (         ) = (   ), */
/*        ( -SR-SI*I   CR+CI*I )   ( YR+YI*I )   ( 0 ) */

/*     where CR, CI, SR, SI, XR, XI, YR, YI are real numbers and I is the */
/*     imaginary unit, I = SQRT(-1). Z is a non-negative real number. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     XR, XI, (input) DOUBLE PRECISION */
/*     YR, YI  (input) DOUBLE PRECISION */
/*             The given real scalars XR, XI, YR, YI. */

/*     CR, CI, (output) DOUBLE PRECISION */
/*     SR, SI, (output) DOUBLE PRECISION */
/*     Z       (output) DOUBLE PRECISION */
/*             The computed real scalars CR, CI, SR, SI, Z, defining the */
/*             complex Givens rotation and Z. */

/*     NUMERICAL ASPECTS */

/*     The subroutine avoids unnecessary overflow. */

/*     FURTHER COMMENTS */

/*     In the interest of speed, this routine does not check the input */
/*     for errors. */

/*     CONTRIBUTOR */

/*     T. Penzl, Technical University Chemnitz, Germany, Aug. 1998. */

/*     REVISIONS */

/*     Sep. 1998 (V. Sima). */

/*     ****************************************************************** */

/*      .. Parameters .. */
/*      .. Scalar Arguments .. */
/*      .. Intrinsic Functions .. */
/*      .. Executable Statements .. */

/* Computing MAX */
#line 74 "SG03BY.f"
    d__1 = abs(*xr), d__2 = abs(*xi), d__1 = max(d__1,d__2), d__2 = abs(*yr), 
	    d__1 = max(d__1,d__2), d__2 = abs(*yi);
#line 74 "SG03BY.f"
    *z__ = max(d__1,d__2);

#line 76 "SG03BY.f"
    if (*z__ == 0.) {
#line 77 "SG03BY.f"
	*cr = 1.;
#line 78 "SG03BY.f"
	*ci = 0.;
#line 79 "SG03BY.f"
	*sr = 0.;
#line 80 "SG03BY.f"
	*si = 0.;
#line 81 "SG03BY.f"
    } else {
/* Computing 2nd power */
#line 82 "SG03BY.f"
	d__1 = *xr / *z__;
/* Computing 2nd power */
#line 82 "SG03BY.f"
	d__2 = *xi / *z__;
/* Computing 2nd power */
#line 82 "SG03BY.f"
	d__3 = *yr / *z__;
/* Computing 2nd power */
#line 82 "SG03BY.f"
	d__4 = *yi / *z__;
#line 82 "SG03BY.f"
	*z__ *= sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4);
#line 84 "SG03BY.f"
	*cr = *xr / *z__;
#line 85 "SG03BY.f"
	*ci = *xi / *z__;
#line 86 "SG03BY.f"
	*sr = *yr / *z__;
#line 87 "SG03BY.f"
	*si = *yi / *z__;
#line 88 "SG03BY.f"
    }

#line 90 "SG03BY.f"
    return 0;

/* *** Last line of SG03BY *** */
} /* sg03by_ */
