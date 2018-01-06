#line 1 "SB04MY.f"
/* SB04MY.f -- translated by f2c (version 20100827).
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

#line 1 "SB04MY.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int sb04my_(integer *n, integer *m, integer *ind, doublereal 
	*a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, doublereal *d__, integer *ipr, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, i2, k1, k2, m1;
    extern /* Subroutine */ int sb04mw_(integer *, doublereal *, integer *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *), daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


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

/*     To construct and solve a linear algebraic system of order M whose */
/*     coefficient matrix is in upper Hessenberg form. Such systems */
/*     appear when solving Sylvester equations using the Hessenberg-Schur */
/*     method. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix B.  N >= 0. */

/*     M       (input) INTEGER */
/*             The order of the matrix A.  M >= 0. */

/*     IND     (input) INTEGER */
/*             The index of the column in C to be computed.  IND >= 1. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,M) */
/*             The leading M-by-M part of this array must contain an */
/*             upper Hessenberg matrix. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,N) */
/*             The leading N-by-N part of this array must contain a */
/*             matrix in real Schur form. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the coefficient matrix C of the equation. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the matrix C with column IND updated. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,M). */

/*     Workspace */

/*     D       DOUBLE PRECISION array, dimension (M*(M+1)/2+2*M) */

/*     IPR     INTEGER array, dimension (2*M) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             > 0:  if INFO = IND, a singular matrix was encountered. */

/*     METHOD */

/*     A special linear algebraic system of order M, with coefficient */
/*     matrix in upper Hessenberg form is constructed and solved. The */
/*     coefficient matrix is stored compactly, row-wise. */

/*     REFERENCES */

/*     [1] Golub, G.H., Nash, S. and Van Loan, C.F. */
/*         A Hessenberg-Schur method for the problem AX + XB = C. */
/*         IEEE Trans. Auto. Contr., AC-24, pp. 909-913, 1979. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Sep. 1997. */
/*     Supersedes Release 2.0 routine SB04AY by G. Golub, S. Nash, and */
/*     C. Van Loan, Stanford University, California, United States of */
/*     America, January 1982. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Hessenberg form, orthogonal transformation, real Schur form, */
/*     Sylvester equation. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 124 "SB04MY.f"
    /* Parameter adjustments */
#line 124 "SB04MY.f"
    a_dim1 = *lda;
#line 124 "SB04MY.f"
    a_offset = 1 + a_dim1;
#line 124 "SB04MY.f"
    a -= a_offset;
#line 124 "SB04MY.f"
    b_dim1 = *ldb;
#line 124 "SB04MY.f"
    b_offset = 1 + b_dim1;
#line 124 "SB04MY.f"
    b -= b_offset;
#line 124 "SB04MY.f"
    c_dim1 = *ldc;
#line 124 "SB04MY.f"
    c_offset = 1 + c_dim1;
#line 124 "SB04MY.f"
    c__ -= c_offset;
#line 124 "SB04MY.f"
    --d__;
#line 124 "SB04MY.f"
    --ipr;
#line 124 "SB04MY.f"

#line 124 "SB04MY.f"
    /* Function Body */
#line 124 "SB04MY.f"
    i__1 = *n;
#line 124 "SB04MY.f"
    for (i__ = *ind + 1; i__ <= i__1; ++i__) {
#line 125 "SB04MY.f"
	d__1 = -b[*ind + i__ * b_dim1];
#line 125 "SB04MY.f"
	daxpy_(m, &d__1, &c__[i__ * c_dim1 + 1], &c__1, &c__[*ind * c_dim1 + 
		1], &c__1);
#line 126 "SB04MY.f"
/* L20: */
#line 126 "SB04MY.f"
    }

#line 128 "SB04MY.f"
    m1 = *m + 1;
#line 129 "SB04MY.f"
    i2 = *m * m1 / 2 + m1;
#line 130 "SB04MY.f"
    k2 = 1;
#line 131 "SB04MY.f"
    k = *m;

/*     Construct the linear algebraic system of order M. */

#line 135 "SB04MY.f"
    i__1 = *m;
#line 135 "SB04MY.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 136 "SB04MY.f"
	j = m1 - k;
#line 137 "SB04MY.f"
	dcopy_(&k, &a[i__ + j * a_dim1], lda, &d__[k2], &c__1);
#line 138 "SB04MY.f"
	k1 = k2;
#line 139 "SB04MY.f"
	k2 += k;
#line 140 "SB04MY.f"
	if (i__ > 1) {
#line 141 "SB04MY.f"
	    ++k1;
#line 142 "SB04MY.f"
	    --k;
#line 143 "SB04MY.f"
	}
#line 144 "SB04MY.f"
	d__[k1] += b[*ind + *ind * b_dim1];

/*        Store the right hand side. */

#line 148 "SB04MY.f"
	d__[i2] = c__[i__ + *ind * c_dim1];
#line 149 "SB04MY.f"
	++i2;
#line 150 "SB04MY.f"
/* L40: */
#line 150 "SB04MY.f"
    }

/*     Solve the linear algebraic system and store the solution in C. */

#line 154 "SB04MY.f"
    sb04mw_(m, &d__[1], &ipr[1], info);

#line 156 "SB04MY.f"
    if (*info != 0) {
#line 157 "SB04MY.f"
	*info = *ind;
#line 158 "SB04MY.f"
    } else {

#line 160 "SB04MY.f"
	i__1 = *m;
#line 160 "SB04MY.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 161 "SB04MY.f"
	    c__[i__ + *ind * c_dim1] = d__[ipr[i__]];
#line 162 "SB04MY.f"
/* L60: */
#line 162 "SB04MY.f"
	}

#line 164 "SB04MY.f"
    }

#line 166 "SB04MY.f"
    return 0;
/* *** Last line of SB04MY *** */
} /* sb04my_ */

