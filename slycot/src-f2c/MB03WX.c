#line 1 "MB03WX.f"
/* MB03WX.f -- translated by f2c (version 20100827).
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

#line 1 "MB03WX.f"
/* Subroutine */ int mb03wx_(integer *n, integer *p, doublereal *t, integer *
	ldt1, integer *ldt2, doublereal *wr, doublereal *wi, integer *info)
{
    /* System generated locals */
    integer t_dim1, t_dim2, t_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, i1;
    static doublereal a11, a12, a21, a22, cs, t11, t12, t22, sn;
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

/*     To compute the eigenvalues of a product of matrices, */
/*     T = T_1*T_2*...*T_p, where T_1 is an upper quasi-triangular */
/*     matrix and T_2, ..., T_p are upper triangular matrices. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix T.  N >= 0. */

/*     P       (input) INTEGER */
/*             The number of matrices in the product T_1*T_2*...*T_p. */
/*             P >= 1. */

/*     T       (input) DOUBLE PRECISION array, dimension (LDT1,LDT2,P) */
/*             The leading N-by-N part of T(*,*,1) must contain the upper */
/*             quasi-triangular matrix T_1 and the leading N-by-N part of */
/*             T(*,*,j) for j > 1 must contain the upper-triangular */
/*             matrix T_j, j = 2, ..., p. */
/*             The elements below the subdiagonal of T(*,*,1) and below */
/*             the diagonal of T(*,*,j), j = 2, ..., p, are not */
/*             referenced. */

/*     LDT1    INTEGER */
/*             The first leading dimension of the array T. */
/*             LDT1 >= max(1,N). */

/*     LDT2    INTEGER */
/*             The second leading dimension of the array T. */
/*             LDT2 >= max(1,N). */

/*     WR, WI  (output) DOUBLE PRECISION arrays, dimension (N) */
/*             The real and imaginary parts, respectively, of the */
/*             eigenvalues of T. The eigenvalues are stored in the same */
/*             order as on the diagonal of T_1. If T(i:i+1,i:i+1,1) is a */
/*             2-by-2 diagonal block with complex conjugated eigenvalues */
/*             then WI(i) > 0 and WI(i+1) = -WI(i). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, February 1999. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Eigenvalue, eigenvalue decomposition, periodic systems, */
/*     real Schur form, triangular form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 100 "MB03WX.f"
    /* Parameter adjustments */
#line 100 "MB03WX.f"
    t_dim1 = *ldt1;
#line 100 "MB03WX.f"
    t_dim2 = *ldt2;
#line 100 "MB03WX.f"
    t_offset = 1 + t_dim1 * (1 + t_dim2);
#line 100 "MB03WX.f"
    t -= t_offset;
#line 100 "MB03WX.f"
    --wr;
#line 100 "MB03WX.f"
    --wi;
#line 100 "MB03WX.f"

#line 100 "MB03WX.f"
    /* Function Body */
#line 100 "MB03WX.f"
    *info = 0;

/*     Test the input scalar arguments. */

#line 104 "MB03WX.f"
    if (*n < 0) {
#line 105 "MB03WX.f"
	*info = -1;
#line 106 "MB03WX.f"
    } else if (*p < 1) {
#line 107 "MB03WX.f"
	*info = -2;
#line 108 "MB03WX.f"
    } else if (*ldt1 < max(1,*n)) {
#line 109 "MB03WX.f"
	*info = -4;
#line 110 "MB03WX.f"
    } else if (*ldt2 < max(1,*n)) {
#line 111 "MB03WX.f"
	*info = -5;
#line 112 "MB03WX.f"
    }

#line 114 "MB03WX.f"
    if (*info != 0) {

/*        Error return. */

#line 118 "MB03WX.f"
	i__1 = -(*info);
#line 118 "MB03WX.f"
	xerbla_("MB03WX", &i__1, (ftnlen)6);
#line 119 "MB03WX.f"
	return 0;
#line 120 "MB03WX.f"
    }

#line 122 "MB03WX.f"
    inext = 1;
#line 123 "MB03WX.f"
    i__1 = *n;
#line 123 "MB03WX.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 124 "MB03WX.f"
	if (i__ < inext) {
#line 124 "MB03WX.f"
	    goto L30;
#line 124 "MB03WX.f"
	}
#line 126 "MB03WX.f"
	if (i__ != *n) {
#line 127 "MB03WX.f"
	    if (t[i__ + 1 + (i__ + t_dim2) * t_dim1] != 0.) {

/*              A pair of eigenvalues. First compute the corresponding */
/*              elements of T(I:I+1,I:I+1). */

#line 132 "MB03WX.f"
		inext = i__ + 2;
#line 133 "MB03WX.f"
		i1 = i__ + 1;
#line 134 "MB03WX.f"
		t11 = 1.;
#line 135 "MB03WX.f"
		t12 = 0.;
#line 136 "MB03WX.f"
		t22 = 1.;

#line 138 "MB03WX.f"
		i__2 = *p;
#line 138 "MB03WX.f"
		for (j = 2; j <= i__2; ++j) {
#line 139 "MB03WX.f"
		    t22 *= t[i1 + (i1 + j * t_dim2) * t_dim1];
#line 140 "MB03WX.f"
		    t12 = t11 * t[i__ + (i1 + j * t_dim2) * t_dim1] + t12 * t[
			    i1 + (i1 + j * t_dim2) * t_dim1];
#line 141 "MB03WX.f"
		    t11 *= t[i__ + (i__ + j * t_dim2) * t_dim1];
#line 142 "MB03WX.f"
/* L10: */
#line 142 "MB03WX.f"
		}

#line 144 "MB03WX.f"
		a11 = t[i__ + (i__ + t_dim2) * t_dim1] * t11;
#line 145 "MB03WX.f"
		a12 = t[i__ + (i__ + t_dim2) * t_dim1] * t12 + t[i__ + (i1 + 
			t_dim2) * t_dim1] * t22;
#line 146 "MB03WX.f"
		a21 = t[i1 + (i__ + t_dim2) * t_dim1] * t11;
#line 147 "MB03WX.f"
		a22 = t[i1 + (i__ + t_dim2) * t_dim1] * t12 + t[i1 + (i1 + 
			t_dim2) * t_dim1] * t22;

#line 149 "MB03WX.f"
		dlanv2_(&a11, &a12, &a21, &a22, &wr[i__], &wi[i__], &wr[i1], &
			wi[i1], &cs, &sn);
#line 151 "MB03WX.f"
		goto L30;
#line 152 "MB03WX.f"
	    }
#line 153 "MB03WX.f"
	}

/*        Simple eigenvalue. Compute the corresponding element of T(I,I). */

#line 157 "MB03WX.f"
	inext = i__ + 1;
#line 158 "MB03WX.f"
	t11 = 1.;

#line 160 "MB03WX.f"
	i__2 = *p;
#line 160 "MB03WX.f"
	for (j = 1; j <= i__2; ++j) {
#line 161 "MB03WX.f"
	    t11 *= t[i__ + (i__ + j * t_dim2) * t_dim1];
#line 162 "MB03WX.f"
/* L20: */
#line 162 "MB03WX.f"
	}

#line 164 "MB03WX.f"
	wr[i__] = t11;
#line 165 "MB03WX.f"
	wi[i__] = 0.;
#line 166 "MB03WX.f"
L30:
#line 166 "MB03WX.f"
	;
#line 166 "MB03WX.f"
    }

#line 168 "MB03WX.f"
    return 0;
/* *** Last line of MB03WX *** */
} /* mb03wx_ */

