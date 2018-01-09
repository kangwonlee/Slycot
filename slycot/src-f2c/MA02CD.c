#line 1 "MA02CD.f"
/* MA02CD.f -- translated by f2c (version 20100827).
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

#line 1 "MA02CD.f"
/* Subroutine */ int ma02cd_(integer *n, integer *kl, integer *ku, doublereal 
	*a, integer *lda)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, i1, lda1;
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

/*     To compute the pertranspose of a central band of a square matrix. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the square matrix A.  N >= 0. */

/*     KL      (input) INTEGER */
/*             The number of subdiagonals of A to be pertransposed. */
/*             0 <= KL <= N-1. */

/*     KU      (input) INTEGER */
/*             The number of superdiagonals of A to be pertransposed. */
/*             0 <= KU <= N-1. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain a square matrix whose central band formed from */
/*             the KL subdiagonals, the main diagonal and the KU */
/*             superdiagonals will be pertransposed. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the matrix A with its central band (the KL subdiagonals, */
/*             the main diagonal and the KU superdiagonals) pertransposed */
/*             (that is the elements of each antidiagonal appear in */
/*             reversed order). This is equivalent to forming P*B'*P, */
/*             where B is the matrix formed from the central band of A */
/*             and P is a permutation matrix with ones down the secondary */
/*             diagonal. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, March 1998. */
/*     Based on the RASP routine DMPTR. */

/*     REVISIONS */

/*     A. Varga, December 2001. */
/*     V. Sima, March 2004. */

/*    ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Quick return if possible. */

#line 84 "MA02CD.f"
    /* Parameter adjustments */
#line 84 "MA02CD.f"
    a_dim1 = *lda;
#line 84 "MA02CD.f"
    a_offset = 1 + a_dim1;
#line 84 "MA02CD.f"
    a -= a_offset;
#line 84 "MA02CD.f"

#line 84 "MA02CD.f"
    /* Function Body */
#line 84 "MA02CD.f"
    if (*n <= 1) {
#line 84 "MA02CD.f"
	return 0;
#line 84 "MA02CD.f"
    }

#line 87 "MA02CD.f"
    lda1 = *lda + 1;

/*     Pertranspose the KL subdiagonals. */

/* Computing MIN */
#line 91 "MA02CD.f"
    i__2 = *kl, i__3 = *n - 2;
#line 91 "MA02CD.f"
    i__1 = min(i__2,i__3);
#line 91 "MA02CD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 92 "MA02CD.f"
	i1 = (*n - i__) / 2;
#line 93 "MA02CD.f"
	if (i1 > 0) {
#line 93 "MA02CD.f"
	    i__2 = -lda1;
#line 93 "MA02CD.f"
	    dswap_(&i1, &a[i__ + 1 + a_dim1], &lda1, &a[*n - i1 + 1 + (*n - 
		    i1 + 1 - i__) * a_dim1], &i__2);
#line 93 "MA02CD.f"
	}
#line 95 "MA02CD.f"
/* L10: */
#line 95 "MA02CD.f"
    }

/*     Pertranspose the KU superdiagonals. */

/* Computing MIN */
#line 99 "MA02CD.f"
    i__2 = *ku, i__3 = *n - 2;
#line 99 "MA02CD.f"
    i__1 = min(i__2,i__3);
#line 99 "MA02CD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 100 "MA02CD.f"
	i1 = (*n - i__) / 2;
#line 101 "MA02CD.f"
	if (i1 > 0) {
#line 101 "MA02CD.f"
	    i__2 = -lda1;
#line 101 "MA02CD.f"
	    dswap_(&i1, &a[(i__ + 1) * a_dim1 + 1], &lda1, &a[*n - i1 + 1 - 
		    i__ + (*n - i1 + 1) * a_dim1], &i__2);
#line 101 "MA02CD.f"
	}
#line 103 "MA02CD.f"
/* L20: */
#line 103 "MA02CD.f"
    }

/*     Pertranspose the diagonal. */

#line 107 "MA02CD.f"
    i1 = *n / 2;
#line 108 "MA02CD.f"
    if (i1 > 0) {
#line 108 "MA02CD.f"
	i__1 = -lda1;
#line 108 "MA02CD.f"
	dswap_(&i1, &a[a_dim1 + 1], &lda1, &a[*n - i1 + 1 + (*n - i1 + 1) * 
		a_dim1], &i__1);
#line 108 "MA02CD.f"
    }

#line 111 "MA02CD.f"
    return 0;
/* *** Last line of MA02CD *** */
} /* ma02cd_ */
