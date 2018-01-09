#line 1 "MB03QX.f"
/* MB03QX.f -- translated by f2c (version 20100827).
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

#line 1 "MB03QX.f"
/* Subroutine */ int mb03qx_(integer *n, doublereal *t, integer *ldt, 
	doublereal *wr, doublereal *wi, integer *info)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1;

    /* Local variables */
    static integer i__, i1;
    static doublereal a11, a12, a21, a22, cs, sn;
    static integer inext;
    extern /* Subroutine */ int dlanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
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

/*     To compute the eigenvalues of an upper quasi-triangular matrix. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix T.  N >= 0. */

/*     T       (input) DOUBLE PRECISION array, dimension(LDT,N) */
/*             The upper quasi-triangular matrix T. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T.  LDT >= max(1,N). */

/*     WR, WI  (output) DOUBLE PRECISION arrays, dimension (N) */
/*             The real and imaginary parts, respectively, of the */
/*             eigenvalues of T. The eigenvalues are stored in the same */
/*             order as on the diagonal of T. If T(i:i+1,i:i+1) is a */
/*             2-by-2 diagonal block with complex conjugated eigenvalues */
/*             then WI(i) > 0 and WI(i+1) = -WI(i). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen, */
/*     March 1998. Based on the RASP routine SEIG. */

/*     ****************************************************************** */
/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 74 "MB03QX.f"
    /* Parameter adjustments */
#line 74 "MB03QX.f"
    t_dim1 = *ldt;
#line 74 "MB03QX.f"
    t_offset = 1 + t_dim1;
#line 74 "MB03QX.f"
    t -= t_offset;
#line 74 "MB03QX.f"
    --wr;
#line 74 "MB03QX.f"
    --wi;
#line 74 "MB03QX.f"

#line 74 "MB03QX.f"
    /* Function Body */
#line 74 "MB03QX.f"
    *info = 0;

/*     Test the input scalar arguments. */

#line 78 "MB03QX.f"
    if (*n < 0) {
#line 79 "MB03QX.f"
	*info = -1;
#line 80 "MB03QX.f"
    } else if (*ldt < max(1,*n)) {
#line 81 "MB03QX.f"
	*info = -3;
#line 82 "MB03QX.f"
    }

#line 84 "MB03QX.f"
    if (*info != 0) {

/*        Error return. */

#line 88 "MB03QX.f"
	i__1 = -(*info);
#line 88 "MB03QX.f"
	xerbla_("MB03QX", &i__1, (ftnlen)6);
#line 89 "MB03QX.f"
	return 0;
#line 90 "MB03QX.f"
    }

#line 92 "MB03QX.f"
    inext = 1;
#line 93 "MB03QX.f"
    i__1 = *n;
#line 93 "MB03QX.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 94 "MB03QX.f"
	if (i__ < inext) {
#line 94 "MB03QX.f"
	    goto L10;
#line 94 "MB03QX.f"
	}
#line 96 "MB03QX.f"
	if (i__ != *n) {
#line 97 "MB03QX.f"
	    if (t[i__ + 1 + i__ * t_dim1] != 0.) {

/*              A pair of eigenvalues. */

#line 101 "MB03QX.f"
		inext = i__ + 2;
#line 102 "MB03QX.f"
		i1 = i__ + 1;
#line 103 "MB03QX.f"
		a11 = t[i__ + i__ * t_dim1];
#line 104 "MB03QX.f"
		a12 = t[i__ + i1 * t_dim1];
#line 105 "MB03QX.f"
		a21 = t[i1 + i__ * t_dim1];
#line 106 "MB03QX.f"
		a22 = t[i1 + i1 * t_dim1];
#line 107 "MB03QX.f"
		dlanv2_(&a11, &a12, &a21, &a22, &wr[i__], &wi[i__], &wr[i1], &
			wi[i1], &cs, &sn);
#line 109 "MB03QX.f"
		goto L10;
#line 110 "MB03QX.f"
	    }
#line 111 "MB03QX.f"
	}

/*        Simple eigenvalue. */

#line 115 "MB03QX.f"
	inext = i__ + 1;
#line 116 "MB03QX.f"
	wr[i__] = t[i__ + i__ * t_dim1];
#line 117 "MB03QX.f"
	wi[i__] = 0.;
#line 118 "MB03QX.f"
L10:
#line 118 "MB03QX.f"
	;
#line 118 "MB03QX.f"
    }

#line 120 "MB03QX.f"
    return 0;
/* *** Last line of MB03QX *** */
} /* mb03qx_ */

