#line 1 "SB04MU.f"
/* SB04MU.f -- translated by f2c (version 20100827).
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

#line 1 "SB04MU.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int sb04mu_(integer *n, integer *m, integer *ind, doublereal 
	*a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, doublereal *d__, integer *ipr, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, i2, k1, k2, m2, ind1;
    static doublereal temp;
    extern /* Subroutine */ int sb04mr_(integer *, doublereal *, integer *, 
	    integer *), daxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *);


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

/*     To construct and solve a linear algebraic system of order 2*M */
/*     whose coefficient matrix has zeros below the second subdiagonal. */
/*     Such systems appear when solving continuous-time Sylvester */
/*     equations using the Hessenberg-Schur method. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix B.  N >= 0. */

/*     M       (input) INTEGER */
/*             The order of the matrix A.  M >= 0. */

/*     IND     (input) INTEGER */
/*             IND and IND - 1 specify the indices of the columns in C */
/*             to be computed.  IND > 1. */

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
/*             the matrix C with columns IND-1 and IND updated. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,M). */

/*     Workspace */

/*     D       DOUBLE PRECISION array, dimension (2*M*M+7*M) */

/*     IPR     INTEGER array, dimension (4*M) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             > 0:  if INFO = IND, a singular matrix was encountered. */

/*     METHOD */

/*     A special linear algebraic system of order 2*M, whose coefficient */
/*     matrix has zeros below the second subdiagonal is constructed and */
/*     solved. The coefficient matrix is stored compactly, row-wise. */

/*     REFERENCES */

/*     [1] Golub, G.H., Nash, S. and Van Loan, C.F. */
/*         A Hessenberg-Schur method for the problem AX + XB = C. */
/*         IEEE Trans. Auto. Contr., AC-24, pp. 909-913, 1979. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Sep. 1997. */
/*     Supersedes Release 2.0 routine SB04AU by G. Golub, S. Nash, and */
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

#line 128 "SB04MU.f"
    /* Parameter adjustments */
#line 128 "SB04MU.f"
    a_dim1 = *lda;
#line 128 "SB04MU.f"
    a_offset = 1 + a_dim1;
#line 128 "SB04MU.f"
    a -= a_offset;
#line 128 "SB04MU.f"
    b_dim1 = *ldb;
#line 128 "SB04MU.f"
    b_offset = 1 + b_dim1;
#line 128 "SB04MU.f"
    b -= b_offset;
#line 128 "SB04MU.f"
    c_dim1 = *ldc;
#line 128 "SB04MU.f"
    c_offset = 1 + c_dim1;
#line 128 "SB04MU.f"
    c__ -= c_offset;
#line 128 "SB04MU.f"
    --d__;
#line 128 "SB04MU.f"
    --ipr;
#line 128 "SB04MU.f"

#line 128 "SB04MU.f"
    /* Function Body */
#line 128 "SB04MU.f"
    ind1 = *ind - 1;

#line 130 "SB04MU.f"
    i__1 = *n;
#line 130 "SB04MU.f"
    for (i__ = *ind + 1; i__ <= i__1; ++i__) {
#line 131 "SB04MU.f"
	d__1 = -b[ind1 + i__ * b_dim1];
#line 131 "SB04MU.f"
	daxpy_(m, &d__1, &c__[i__ * c_dim1 + 1], &c__1, &c__[ind1 * c_dim1 + 
		1], &c__1);
#line 132 "SB04MU.f"
	d__1 = -b[*ind + i__ * b_dim1];
#line 132 "SB04MU.f"
	daxpy_(m, &d__1, &c__[i__ * c_dim1 + 1], &c__1, &c__[*ind * c_dim1 + 
		1], &c__1);
#line 133 "SB04MU.f"
/* L20: */
#line 133 "SB04MU.f"
    }

/*     Construct the linear algebraic system of order 2*M. */

#line 137 "SB04MU.f"
    k1 = -1;
#line 138 "SB04MU.f"
    m2 = *m << 1;
#line 139 "SB04MU.f"
    i2 = *m * (m2 + 5);
#line 140 "SB04MU.f"
    k = m2;

#line 142 "SB04MU.f"
    i__1 = *m;
#line 142 "SB04MU.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/* Computing MAX */
#line 144 "SB04MU.f"
	i__2 = 1, i__3 = i__ - 1;
#line 144 "SB04MU.f"
	i__4 = *m;
#line 144 "SB04MU.f"
	for (j = max(i__2,i__3); j <= i__4; ++j) {
#line 145 "SB04MU.f"
	    k1 += 2;
#line 146 "SB04MU.f"
	    k2 = k1 + k;
#line 147 "SB04MU.f"
	    temp = a[i__ + j * a_dim1];
#line 148 "SB04MU.f"
	    if (i__ != j) {
#line 149 "SB04MU.f"
		d__[k1] = temp;
#line 150 "SB04MU.f"
		d__[k1 + 1] = 0.;
#line 151 "SB04MU.f"
		if (j > i__) {
#line 151 "SB04MU.f"
		    d__[k2] = 0.;
#line 151 "SB04MU.f"
		}
#line 152 "SB04MU.f"
		d__[k2 + 1] = temp;
#line 153 "SB04MU.f"
	    } else {
#line 154 "SB04MU.f"
		d__[k1] = temp + b[ind1 + ind1 * b_dim1];
#line 155 "SB04MU.f"
		d__[k1 + 1] = b[ind1 + *ind * b_dim1];
#line 156 "SB04MU.f"
		d__[k2] = b[*ind + ind1 * b_dim1];
#line 157 "SB04MU.f"
		d__[k2 + 1] = temp + b[*ind + *ind * b_dim1];
#line 158 "SB04MU.f"
	    }
#line 159 "SB04MU.f"
/* L40: */
#line 159 "SB04MU.f"
	}

#line 161 "SB04MU.f"
	k1 = k2;
#line 162 "SB04MU.f"
	k -= min(2,i__);

/*        Store the right hand side. */

#line 166 "SB04MU.f"
	i2 += 2;
#line 167 "SB04MU.f"
	d__[i2] = c__[i__ + *ind * c_dim1];
#line 168 "SB04MU.f"
	d__[i2 - 1] = c__[i__ + ind1 * c_dim1];
#line 169 "SB04MU.f"
/* L60: */
#line 169 "SB04MU.f"
    }

/*     Solve the linear algebraic system and store the solution in C. */

#line 173 "SB04MU.f"
    sb04mr_(&m2, &d__[1], &ipr[1], info);

#line 175 "SB04MU.f"
    if (*info != 0) {
#line 176 "SB04MU.f"
	*info = *ind;
#line 177 "SB04MU.f"
    } else {
#line 178 "SB04MU.f"
	i2 = 0;

#line 180 "SB04MU.f"
	i__1 = *m;
#line 180 "SB04MU.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 181 "SB04MU.f"
	    i2 += 2;
#line 182 "SB04MU.f"
	    c__[i__ + ind1 * c_dim1] = d__[ipr[i2 - 1]];
#line 183 "SB04MU.f"
	    c__[i__ + *ind * c_dim1] = d__[ipr[i2]];
#line 184 "SB04MU.f"
/* L80: */
#line 184 "SB04MU.f"
	}

#line 186 "SB04MU.f"
    }

#line 188 "SB04MU.f"
    return 0;
/* *** Last line of SB04MU *** */
} /* sb04mu_ */

