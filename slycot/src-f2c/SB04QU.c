#line 1 "SB04QU.f"
/* SB04QU.f -- translated by f2c (version 20100827).
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

#line 1 "SB04QU.f"
/* Table of constant values */

static integer c__0 = 0;
static integer c__1 = 1;

/* Subroutine */ int sb04qu_(integer *n, integer *m, integer *ind, doublereal 
	*a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, doublereal *d__, integer *ipr, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, j, k, i2, k1, k2, m2;
    static doublereal dum[1];
    static integer ind1;
    static doublereal temp;
    extern /* Subroutine */ int sb04qr_(integer *, doublereal *, integer *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *), daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dtrmv_(char *, char *, char *
	    , integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);


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
/*     whose coefficient matrix has zeros below the third subdiagonal, */
/*     and zero elements on the third subdiagonal with even column */
/*     indices. Such systems appear when solving discrete-time Sylvester */
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

/*     D       DOUBLE PRECISION array, dimension (2*M*M+8*M) */

/*     IPR     INTEGER array, dimension (4*M) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             > 0:  if INFO = IND, a singular matrix was encountered. */

/*     METHOD */

/*     A special linear algebraic system of order 2*M, whose coefficient */
/*     matrix has zeros below the third subdiagonal and zero elements on */
/*     the third subdiagonal with even column indices, is constructed and */
/*     solved. The coefficient matrix is stored compactly, row-wise. */

/*     REFERENCES */

/*     [1] Golub, G.H., Nash, S. and Van Loan, C.F. */
/*         A Hessenberg-Schur method for the problem AX + XB = C. */
/*         IEEE Trans. Auto. Contr., AC-24, pp. 909-913, 1979. */

/*     [2] Sima, V. */
/*         Algorithms for Linear-quadratic Optimization. */
/*         Marcel Dekker, Inc., New York, 1996. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTORS */

/*     D. Sima, University of Bucharest, May 2000. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Hessenberg form, orthogonal transformation, real Schur form, */
/*     Sylvester equation. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 133 "SB04QU.f"
    /* Parameter adjustments */
#line 133 "SB04QU.f"
    a_dim1 = *lda;
#line 133 "SB04QU.f"
    a_offset = 1 + a_dim1;
#line 133 "SB04QU.f"
    a -= a_offset;
#line 133 "SB04QU.f"
    b_dim1 = *ldb;
#line 133 "SB04QU.f"
    b_offset = 1 + b_dim1;
#line 133 "SB04QU.f"
    b -= b_offset;
#line 133 "SB04QU.f"
    c_dim1 = *ldc;
#line 133 "SB04QU.f"
    c_offset = 1 + c_dim1;
#line 133 "SB04QU.f"
    c__ -= c_offset;
#line 133 "SB04QU.f"
    --d__;
#line 133 "SB04QU.f"
    --ipr;
#line 133 "SB04QU.f"

#line 133 "SB04QU.f"
    /* Function Body */
#line 133 "SB04QU.f"
    ind1 = *ind - 1;

#line 135 "SB04QU.f"
    if (*ind < *n) {
#line 136 "SB04QU.f"
	dum[0] = 0.;
#line 137 "SB04QU.f"
	dcopy_(m, dum, &c__0, &d__[1], &c__1);
#line 138 "SB04QU.f"
	i__1 = *n;
#line 138 "SB04QU.f"
	for (i__ = *ind + 1; i__ <= i__1; ++i__) {
#line 139 "SB04QU.f"
	    daxpy_(m, &b[ind1 + i__ * b_dim1], &c__[i__ * c_dim1 + 1], &c__1, 
		    &d__[1], &c__1);
#line 140 "SB04QU.f"
/* L10: */
#line 140 "SB04QU.f"
	}

#line 142 "SB04QU.f"
	i__1 = *m;
#line 142 "SB04QU.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 143 "SB04QU.f"
	    c__[i__ + ind1 * c_dim1] -= a[i__ + (i__ - 1) * a_dim1] * d__[i__ 
		    - 1];
#line 144 "SB04QU.f"
/* L20: */
#line 144 "SB04QU.f"
	}
#line 145 "SB04QU.f"
	dtrmv_("Upper", "No Transpose", "Non Unit", m, &a[a_offset], lda, &
		d__[1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 147 "SB04QU.f"
	i__1 = *m;
#line 147 "SB04QU.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 148 "SB04QU.f"
	    c__[i__ + ind1 * c_dim1] -= d__[i__];
#line 149 "SB04QU.f"
/* L30: */
#line 149 "SB04QU.f"
	}

#line 151 "SB04QU.f"
	dcopy_(m, dum, &c__0, &d__[1], &c__1);
#line 152 "SB04QU.f"
	i__1 = *n;
#line 152 "SB04QU.f"
	for (i__ = *ind + 1; i__ <= i__1; ++i__) {
#line 153 "SB04QU.f"
	    daxpy_(m, &b[*ind + i__ * b_dim1], &c__[i__ * c_dim1 + 1], &c__1, 
		    &d__[1], &c__1);
#line 154 "SB04QU.f"
/* L40: */
#line 154 "SB04QU.f"
	}

#line 156 "SB04QU.f"
	i__1 = *m;
#line 156 "SB04QU.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 157 "SB04QU.f"
	    c__[i__ + *ind * c_dim1] -= a[i__ + (i__ - 1) * a_dim1] * d__[i__ 
		    - 1];
#line 158 "SB04QU.f"
/* L50: */
#line 158 "SB04QU.f"
	}
#line 159 "SB04QU.f"
	dtrmv_("Upper", "No Transpose", "Non Unit", m, &a[a_offset], lda, &
		d__[1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 161 "SB04QU.f"
	i__1 = *m;
#line 161 "SB04QU.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 162 "SB04QU.f"
	    c__[i__ + *ind * c_dim1] -= d__[i__];
#line 163 "SB04QU.f"
/* L60: */
#line 163 "SB04QU.f"
	}
#line 164 "SB04QU.f"
    }

/*     Construct the linear algebraic system of order 2*M. */

#line 168 "SB04QU.f"
    k1 = -1;
#line 169 "SB04QU.f"
    m2 = *m << 1;
#line 170 "SB04QU.f"
    i2 = m2 * (*m + 3);
#line 171 "SB04QU.f"
    k = m2;

#line 173 "SB04QU.f"
    i__1 = *m;
#line 173 "SB04QU.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/* Computing MAX */
#line 175 "SB04QU.f"
	i__2 = 1, i__3 = i__ - 1;
#line 175 "SB04QU.f"
	i__4 = *m;
#line 175 "SB04QU.f"
	for (j = max(i__2,i__3); j <= i__4; ++j) {
#line 176 "SB04QU.f"
	    k1 += 2;
#line 177 "SB04QU.f"
	    k2 = k1 + k;
#line 178 "SB04QU.f"
	    temp = a[i__ + j * a_dim1];
#line 179 "SB04QU.f"
	    d__[k1] = temp * b[ind1 + ind1 * b_dim1];
#line 180 "SB04QU.f"
	    d__[k1 + 1] = temp * b[ind1 + *ind * b_dim1];
#line 181 "SB04QU.f"
	    d__[k2] = temp * b[*ind + ind1 * b_dim1];
#line 182 "SB04QU.f"
	    d__[k2 + 1] = temp * b[*ind + *ind * b_dim1];
#line 183 "SB04QU.f"
	    if (i__ == j) {
#line 184 "SB04QU.f"
		d__[k1] += 1.;
#line 185 "SB04QU.f"
		d__[k2 + 1] += 1.;
#line 186 "SB04QU.f"
	    }
#line 187 "SB04QU.f"
/* L70: */
#line 187 "SB04QU.f"
	}

#line 189 "SB04QU.f"
	k1 = k2;
#line 190 "SB04QU.f"
	if (i__ > 1) {
#line 190 "SB04QU.f"
	    k += -2;
#line 190 "SB04QU.f"
	}

/*        Store the right hand side. */

#line 194 "SB04QU.f"
	i2 += 2;
#line 195 "SB04QU.f"
	d__[i2] = c__[i__ + *ind * c_dim1];
#line 196 "SB04QU.f"
	d__[i2 - 1] = c__[i__ + ind1 * c_dim1];
#line 197 "SB04QU.f"
/* L80: */
#line 197 "SB04QU.f"
    }

/*     Solve the linear algebraic system and store the solution in C. */

#line 201 "SB04QU.f"
    sb04qr_(&m2, &d__[1], &ipr[1], info);

#line 203 "SB04QU.f"
    if (*info != 0) {
#line 204 "SB04QU.f"
	*info = *ind;
#line 205 "SB04QU.f"
    } else {
#line 206 "SB04QU.f"
	i2 = 0;

#line 208 "SB04QU.f"
	i__1 = *m;
#line 208 "SB04QU.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 209 "SB04QU.f"
	    i2 += 2;
#line 210 "SB04QU.f"
	    c__[i__ + ind1 * c_dim1] = d__[ipr[i2 - 1]];
#line 211 "SB04QU.f"
	    c__[i__ + *ind * c_dim1] = d__[ipr[i2]];
#line 212 "SB04QU.f"
/* L90: */
#line 212 "SB04QU.f"
	}

#line 214 "SB04QU.f"
    }

#line 216 "SB04QU.f"
    return 0;
/* *** Last line of SB04QU *** */
} /* sb04qu_ */

