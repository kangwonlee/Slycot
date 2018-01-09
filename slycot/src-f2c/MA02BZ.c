#line 1 "MA02BZ.f"
/* MA02BZ.f -- translated by f2c (version 20100827).
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

#line 1 "MA02BZ.f"
/* Table of constant values */

static integer c_n1 = -1;
static integer c__1 = 1;

/* Subroutine */ int ma02bz_(char *side, integer *m, integer *n, 
	doublecomplex *a, integer *lda, ftnlen side_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k, m2, n2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical bsides;


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

/*     To reverse the order of rows and/or columns of a given matrix A */
/*     by pre-multiplying and/or post-multiplying it, respectively, with */
/*     a permutation matrix P, where P is a square matrix of appropriate */
/*     order, with ones down the secondary diagonal. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SIDE    CHARACTER*1 */
/*             Specifies the operation to be performed, as follows: */
/*             = 'L': the order of rows of A is to be reversed by */
/*                    pre-multiplying A with P; */
/*             = 'R': the order of columns of A is to be reversed by */
/*                    post-multiplying A with P; */
/*             = 'B': both the order of rows and the order of columns */
/*                    of A is to be reversed by pre-multiplying and */
/*                    post-multiplying A with P. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     A       (input/output) COMPLEX*16 array, dimension (LDA,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the given matrix whose rows and/or columns are to */
/*             be permuted. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the matrix P*A if SIDE = 'L', or A*P if SIDE = 'R', or */
/*             P*A*P if SIDE = 'B'. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, March 1998. */
/*     Complex version: V. Sima, Research Institute for Informatics, */
/*     Bucharest, Nov. 2008. */

/*     REVISIONS */

/*     - */

/*    ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

#line 89 "MA02BZ.f"
    /* Parameter adjustments */
#line 89 "MA02BZ.f"
    a_dim1 = *lda;
#line 89 "MA02BZ.f"
    a_offset = 1 + a_dim1;
#line 89 "MA02BZ.f"
    a -= a_offset;
#line 89 "MA02BZ.f"

#line 89 "MA02BZ.f"
    /* Function Body */
#line 89 "MA02BZ.f"
    bsides = lsame_(side, "B", (ftnlen)1, (ftnlen)1);

#line 91 "MA02BZ.f"
    if ((lsame_(side, "L", (ftnlen)1, (ftnlen)1) || bsides) && *m > 1) {

/*        Compute P*A. */

#line 95 "MA02BZ.f"
	m2 = *m / 2;
#line 96 "MA02BZ.f"
	k = *m - m2 + 1;
#line 97 "MA02BZ.f"
	i__1 = *n;
#line 97 "MA02BZ.f"
	for (j = 1; j <= i__1; ++j) {
#line 98 "MA02BZ.f"
	    zswap_(&m2, &a[j * a_dim1 + 1], &c_n1, &a[k + j * a_dim1], &c__1);
#line 99 "MA02BZ.f"
/* L10: */
#line 99 "MA02BZ.f"
	}
#line 100 "MA02BZ.f"
    }
#line 101 "MA02BZ.f"
    if ((lsame_(side, "R", (ftnlen)1, (ftnlen)1) || bsides) && *n > 1) {

/*        Compute A*P. */

#line 105 "MA02BZ.f"
	n2 = *n / 2;
#line 106 "MA02BZ.f"
	k = *n - n2 + 1;
#line 107 "MA02BZ.f"
	i__1 = *m;
#line 107 "MA02BZ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 108 "MA02BZ.f"
	    i__2 = -(*lda);
#line 108 "MA02BZ.f"
	    zswap_(&n2, &a[i__ + a_dim1], &i__2, &a[i__ + k * a_dim1], lda);
#line 109 "MA02BZ.f"
/* L20: */
#line 109 "MA02BZ.f"
	}
#line 110 "MA02BZ.f"
    }

#line 112 "MA02BZ.f"
    return 0;
/* *** Last line of MA02BZ *** */
} /* ma02bz_ */

