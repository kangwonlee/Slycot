#line 1 "MB04OX.f"
/* MB04OX.f -- translated by f2c (version 20100827).
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

#line 1 "MB04OX.f"
/* Subroutine */ int mb04ox_(integer *n, doublereal *a, integer *lda, 
	doublereal *x, integer *incx)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__;
    static doublereal ci, si;
    static integer ix;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dlartg_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);


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

/*     To perform the QR factorization */

/*        (U ) = Q*(R), */
/*        (x')     (0) */

/*     where U and R are n-by-n upper triangular matrices, x is an */
/*     n element vector and Q is an (n+1)-by-(n+1) orthogonal matrix. */

/*     U must be supplied in the n-by-n upper triangular part of the */
/*     array A and this is overwritten by R. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N      (input) INTEGER */
/*            The number of elements of X and the order of the square */
/*            matrix A.  N >= 0. */

/*     A      (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*            On entry, the leading N-by-N upper triangular part of this */
/*            array must contain the upper triangular matrix U. */
/*            On exit, the leading N-by-N upper triangular part of this */
/*            array contains the upper triangular matrix R. */
/*            The strict lower triangle of A is not referenced. */

/*     LDA    INTEGER */
/*            The leading dimension of the array A.  LDA >= max(1,N). */

/*     X      (input/output) DOUBLE PRECISION array, dimension */
/*            (1+(N-1)*INCX) */
/*            On entry, the incremented array X must contain the */
/*            vector x. On exit, the content of X is changed. */

/*     INCX   (input) INTEGER. */
/*            Specifies the increment for the elements of X.  INCX > 0. */

/*     METHOD */

/*     The matrix Q is formed as a sequence of plane rotations in planes */
/*     (1, n+1), (2, n+1), ..., (n, n+1), the rotation in the (j, n+1)th */
/*     plane, Q(j), being chosen to annihilate the jth element of x. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, July 1998. */
/*     Based on the RASP routine DUTUPD. */

/*     REVISIONS */

/*     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */

/*     .. Executable Statements .. */

/*     For efficiency reasons, the parameters are not checked. */

#line 92 "MB04OX.f"
    /* Parameter adjustments */
#line 92 "MB04OX.f"
    a_dim1 = *lda;
#line 92 "MB04OX.f"
    a_offset = 1 + a_dim1;
#line 92 "MB04OX.f"
    a -= a_offset;
#line 92 "MB04OX.f"
    --x;
#line 92 "MB04OX.f"

#line 92 "MB04OX.f"
    /* Function Body */
#line 92 "MB04OX.f"
    ix = 1;

#line 94 "MB04OX.f"
    i__1 = *n - 1;
#line 94 "MB04OX.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 95 "MB04OX.f"
	dlartg_(&a[i__ + i__ * a_dim1], &x[ix], &ci, &si, &temp);
#line 96 "MB04OX.f"
	a[i__ + i__ * a_dim1] = temp;
#line 97 "MB04OX.f"
	ix += *incx;
#line 98 "MB04OX.f"
	i__2 = *n - i__;
#line 98 "MB04OX.f"
	drot_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda, &x[ix], incx, &ci, &
		si);
#line 99 "MB04OX.f"
/* L20: */
#line 99 "MB04OX.f"
    }

#line 101 "MB04OX.f"
    dlartg_(&a[*n + *n * a_dim1], &x[ix], &ci, &si, &temp);
#line 102 "MB04OX.f"
    a[*n + *n * a_dim1] = temp;

#line 104 "MB04OX.f"
    return 0;
/* *** Last line of MB04OX *** */
} /* mb04ox_ */

