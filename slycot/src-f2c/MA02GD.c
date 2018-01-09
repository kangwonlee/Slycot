#line 1 "MA02GD.f"
/* MA02GD.f -- translated by f2c (version 20100827).
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

#line 1 "MA02GD.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int ma02gd_(integer *n, doublereal *a, integer *lda, integer 
	*k1, integer *k2, integer *ipiv, integer *incx)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer j, jp, jx;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);


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

/*     To perform a series of column interchanges on the matrix A. */
/*     One column interchange is initiated for each of columns K1 through */
/*     K2 of A. This is useful for solving linear systems X*A = B, when */
/*     the matrix A has already been factored by LAPACK Library routine */
/*     DGETRF. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of rows of the matrix A.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,*) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the matrix A to which the column interchanges will */
/*             be applied, where M is the largest element of IPIV(K), for */
/*             K = K1, ..., K2. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the permuted matrix. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     K1      (input) INTEGER */
/*             The first element of IPIV for which a column interchange */
/*             will be done. */

/*     K2      (input) INTEGER */
/*             The last element of IPIV for which a column interchange */
/*             will be done. */

/*     IPIV    (input) INTEGER array, dimension (K1+(K2-K1)*abs(INCX)) */
/*             The vector of interchanging (pivot) indices.  Only the */
/*             elements in positions K1 through K2 of IPIV are accessed. */
/*             IPIV(K) = L implies columns K and L are to be */
/*             interchanged. */

/*     INCX    (input) INTEGER */
/*             The increment between successive values of IPIV. */
/*             If INCX is negative, the interchanges are applied in */
/*             reverse order. */

/*     METHOD */

/*     The columns IPIV(K) and K are swapped for K = K1, ..., K2, for */
/*     INCX = 1 (and similarly, for INCX <> 1). */

/*     FURTHER COMMENTS */

/*     This routine is the column-oriented counterpart of the LAPACK */
/*     Library routine DLASWP. The LAPACK Library routine DLAPMT cannot */
/*     be used in this context. To solve the system X*A = B, where A and */
/*     B are N-by-N and M-by-N, respectively, the following statements */
/*     can be used: */

/*         CALL DGETRF( N, N, A, LDA, IPIV, INFO ) */
/*         CALL DTRSM( 'R', 'U', 'N', 'N', M, N, ONE, A, LDA, B, LDB ) */
/*         CALL DTRSM( 'R', 'L', 'N', 'U', M, N, ONE, A, LDA, B, LDB ) */
/*         CALL MA02GD( M, B, LDB, 1, N, IPIV, -1 ) */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2000. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2008. */

/*     KEYWORDS */

/*     Elementary matrix operations, linear algebra. */

/*    ****************************************************************** */

/*      .. Scalar Arguments .. */
/*      .. */
/*      .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible. */

#line 115 "MA02GD.f"
    /* Parameter adjustments */
#line 115 "MA02GD.f"
    a_dim1 = *lda;
#line 115 "MA02GD.f"
    a_offset = 1 + a_dim1;
#line 115 "MA02GD.f"
    a -= a_offset;
#line 115 "MA02GD.f"
    --ipiv;
#line 115 "MA02GD.f"

#line 115 "MA02GD.f"
    /* Function Body */
#line 115 "MA02GD.f"
    if (*incx == 0 || *n == 0) {
#line 115 "MA02GD.f"
	return 0;
#line 115 "MA02GD.f"
    }

/*     Interchange column J with column IPIV(J) for each of columns K1 */
/*     through K2. */

#line 121 "MA02GD.f"
    if (*incx > 0) {
#line 122 "MA02GD.f"
	jx = *k1;
#line 123 "MA02GD.f"
    } else {
#line 124 "MA02GD.f"
	jx = (1 - *k2) * *incx + 1;
#line 125 "MA02GD.f"
    }

#line 127 "MA02GD.f"
    if (*incx == 1) {

#line 129 "MA02GD.f"
	i__1 = *k2;
#line 129 "MA02GD.f"
	for (j = *k1; j <= i__1; ++j) {
#line 130 "MA02GD.f"
	    jp = ipiv[j];
#line 131 "MA02GD.f"
	    if (jp != j) {
#line 131 "MA02GD.f"
		dswap_(n, &a[j * a_dim1 + 1], &c__1, &a[jp * a_dim1 + 1], &
			c__1);
#line 131 "MA02GD.f"
	    }
#line 133 "MA02GD.f"
/* L10: */
#line 133 "MA02GD.f"
	}

#line 135 "MA02GD.f"
    } else if (*incx > 1) {

#line 137 "MA02GD.f"
	i__1 = *k2;
#line 137 "MA02GD.f"
	for (j = *k1; j <= i__1; ++j) {
#line 138 "MA02GD.f"
	    jp = ipiv[jx];
#line 139 "MA02GD.f"
	    if (jp != j) {
#line 139 "MA02GD.f"
		dswap_(n, &a[j * a_dim1 + 1], &c__1, &a[jp * a_dim1 + 1], &
			c__1);
#line 139 "MA02GD.f"
	    }
#line 141 "MA02GD.f"
	    jx += *incx;
#line 142 "MA02GD.f"
/* L20: */
#line 142 "MA02GD.f"
	}

#line 144 "MA02GD.f"
    } else if (*incx < 0) {

#line 146 "MA02GD.f"
	i__1 = *k1;
#line 146 "MA02GD.f"
	for (j = *k2; j >= i__1; --j) {
#line 147 "MA02GD.f"
	    jp = ipiv[jx];
#line 148 "MA02GD.f"
	    if (jp != j) {
#line 148 "MA02GD.f"
		dswap_(n, &a[j * a_dim1 + 1], &c__1, &a[jp * a_dim1 + 1], &
			c__1);
#line 148 "MA02GD.f"
	    }
#line 150 "MA02GD.f"
	    jx += *incx;
#line 151 "MA02GD.f"
/* L30: */
#line 151 "MA02GD.f"
	}

#line 153 "MA02GD.f"
    }

#line 155 "MA02GD.f"
    return 0;

/* *** Last line of MA02GD *** */
} /* ma02gd_ */

