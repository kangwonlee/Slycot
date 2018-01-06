#line 1 "TB01YD.f"
/* TB01YD.f -- translated by f2c (version 20100827).
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

#line 1 "TB01YD.f"
/* Table of constant values */

static integer c_n1 = -1;
static integer c__1 = 1;

/* Subroutine */ int tb01yd_(integer *n, integer *m, integer *p, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1;

    /* Local variables */
    static integer j, nby2;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), xerbla_(char *, integer *, ftnlen);


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

/*     To apply a special similarity transformation to a system given as */
/*     a triple (A,B,C), */

/*        A <-- P * A * P,  B <-- P * B,  C <-- C * P, */

/*     where P is a matrix with 1 on the secondary diagonal, and with 0 */
/*     in the other entries. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A, the number of rows of matrix B */
/*             and the number of columns of matrix C. */
/*             N represents the dimension of the state vector.  N >= 0. */

/*     M       (input) INTEGER. */
/*             The number of columns of matrix B. */
/*             M represents the dimension of input vector.  M >= 0. */

/*     P       (input) INTEGER. */
/*             The number of rows of matrix C. */
/*             P represents the dimension of output vector.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the system state matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the transformed matrix P*A*P. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the system input matrix B. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the transformed matrix P*B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= MAX(1,N) if M > 0. */
/*             LDB >= 1        if M = 0. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the system output matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the transformed matrix C*P. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= MAX(1,P). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit. */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The rows and/or columns of the matrices of the triplet (A,B,C) */
/*     are swapped in a special way. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1998. */


/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2004. */

/*     KEYWORDS */

/*     Matrix algebra, matrix operations, similarity transformation. */

/*  ********************************************************************* */

/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the scalar input arguments. */

#line 128 "TB01YD.f"
    /* Parameter adjustments */
#line 128 "TB01YD.f"
    a_dim1 = *lda;
#line 128 "TB01YD.f"
    a_offset = 1 + a_dim1;
#line 128 "TB01YD.f"
    a -= a_offset;
#line 128 "TB01YD.f"
    b_dim1 = *ldb;
#line 128 "TB01YD.f"
    b_offset = 1 + b_dim1;
#line 128 "TB01YD.f"
    b -= b_offset;
#line 128 "TB01YD.f"
    c_dim1 = *ldc;
#line 128 "TB01YD.f"
    c_offset = 1 + c_dim1;
#line 128 "TB01YD.f"
    c__ -= c_offset;
#line 128 "TB01YD.f"

#line 128 "TB01YD.f"
    /* Function Body */
#line 128 "TB01YD.f"
    *info = 0;

#line 130 "TB01YD.f"
    if (*n < 0) {
#line 131 "TB01YD.f"
	*info = -1;
#line 132 "TB01YD.f"
    } else if (*m < 0) {
#line 133 "TB01YD.f"
	*info = -2;
#line 134 "TB01YD.f"
    } else if (*p < 0) {
#line 135 "TB01YD.f"
	*info = -3;
#line 136 "TB01YD.f"
    } else if (*lda < max(1,*n)) {
#line 137 "TB01YD.f"
	*info = -5;
#line 138 "TB01YD.f"
    } else if (*ldb < 1 || *m > 0 && *ldb < *n) {
#line 139 "TB01YD.f"
	*info = -7;
#line 140 "TB01YD.f"
    } else if (*ldc < max(1,*p)) {
#line 141 "TB01YD.f"
	*info = -9;
#line 142 "TB01YD.f"
    }

#line 144 "TB01YD.f"
    if (*info != 0) {

/*        Error return. */

#line 148 "TB01YD.f"
	i__1 = -(*info);
#line 148 "TB01YD.f"
	xerbla_("TB01YD", &i__1, (ftnlen)6);
#line 149 "TB01YD.f"
	return 0;
#line 150 "TB01YD.f"
    }

#line 152 "TB01YD.f"
    if (*n <= 1) {
#line 152 "TB01YD.f"
	return 0;
#line 152 "TB01YD.f"
    }

/*     Transform the matrix A. */

#line 157 "TB01YD.f"
    nby2 = *n / 2;

#line 159 "TB01YD.f"
    i__1 = nby2;
#line 159 "TB01YD.f"
    for (j = 1; j <= i__1; ++j) {
#line 160 "TB01YD.f"
	dswap_(n, &a[j * a_dim1 + 1], &c_n1, &a[(*n - j + 1) * a_dim1 + 1], &
		c__1);
#line 161 "TB01YD.f"
/* L10: */
#line 161 "TB01YD.f"
    }

#line 163 "TB01YD.f"
    if (*n % 2 != 0 && *n > 2) {
#line 163 "TB01YD.f"
	dswap_(&nby2, &a[nby2 + 2 + (nby2 + 1) * a_dim1], &c_n1, &a[(nby2 + 1)
		 * a_dim1 + 1], &c__1);
#line 163 "TB01YD.f"
    }

#line 166 "TB01YD.f"
    if (*m > 0) {

/*        Transform the matrix B. */

#line 170 "TB01YD.f"
	i__1 = nby2;
#line 170 "TB01YD.f"
	for (j = 1; j <= i__1; ++j) {
#line 171 "TB01YD.f"
	    dswap_(m, &b[j + b_dim1], ldb, &b[*n - j + 1 + b_dim1], ldb);
#line 172 "TB01YD.f"
/* L20: */
#line 172 "TB01YD.f"
	}

#line 174 "TB01YD.f"
    }

#line 176 "TB01YD.f"
    if (*p > 0) {

/*        Transform the matrix C. */

#line 180 "TB01YD.f"
	i__1 = nby2;
#line 180 "TB01YD.f"
	for (j = 1; j <= i__1; ++j) {
#line 181 "TB01YD.f"
	    dswap_(p, &c__[j * c_dim1 + 1], &c__1, &c__[(*n - j + 1) * c_dim1 
		    + 1], &c__1);
#line 182 "TB01YD.f"
/* L30: */
#line 182 "TB01YD.f"
	}

#line 184 "TB01YD.f"
    }

#line 186 "TB01YD.f"
    return 0;
/* *** Last line of TB01YD *** */
} /* tb01yd_ */

