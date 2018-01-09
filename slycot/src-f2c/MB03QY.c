#line 1 "MB03QY.f"
/* MB03QY.f -- translated by f2c (version 20100827).
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

#line 1 "MB03QY.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mb03qy_(integer *n, integer *l, doublereal *a, integer *
	lda, doublereal *u, integer *ldu, doublereal *e1, doublereal *e2, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, u_dim1, u_offset, i__1;

    /* Local variables */
    static integer l1;
    static doublereal cs, sn, ew1, ew2;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dlanv2_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), xerbla_(char *, integer *, ftnlen);


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

/*     To compute the eigenvalues of a selected 2-by-2 diagonal block */
/*     of an upper quasi-triangular matrix, to reduce the selected block */
/*     to the standard form and to split the block in the case of real */
/*     eigenvalues by constructing an orthogonal transformation UT. */
/*     This transformation is applied to A (by similarity) and to */
/*     another matrix U from the right. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A and UT.  N >= 2. */

/*     L       (input) INTEGER */
/*             Specifies the position of the block.  1 <= L < N. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper quasi-triangular matrix A whose */
/*             selected 2-by-2 diagonal block is to be processed. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the upper quasi-triangular matrix A after its selected */
/*             block has been splitt and/or put in the LAPACK standard */
/*             form. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= N. */

/*     U       (input/output) DOUBLE PRECISION array, dimension (LDU,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain a transformation matrix U. */
/*             On exit, the leading N-by-N part of this array contains */
/*             U*UT, where UT is the transformation matrix used to */
/*             split and/or standardize the selected block. */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= N. */

/*     E1, E2  (output) DOUBLE PRECISION */
/*             E1 and E2 contain either the real eigenvalues or the real */
/*             and positive imaginary parts, respectively, of the complex */
/*             eigenvalues of the selected 2-by-2 diagonal block of A. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Let A1 = ( A(L,L)    A(L,L+1)   ) */
/*              ( A(L+1,L)  A(L+1,L+1) ) */
/*     be the specified 2-by-2 diagonal block of matrix A. */
/*     If the eigenvalues of A1 are complex, then they are computed and */
/*     stored in E1 and E2, where the real part is stored in E1 and the */
/*     positive imaginary part in E2. The 2-by-2 block is reduced if */
/*     necessary to the standard form, such that A(L,L) = A(L+1,L+1), and */
/*     A(L,L+1) and A(L+1,L) have oposite signs. If the eigenvalues are */
/*     real, the 2-by-2 block is reduced to an upper triangular form such */
/*     that ABS(A(L,L)) >= ABS(A(L+1,L+1)). */
/*     In both cases, an orthogonal rotation U1' is constructed such that */
/*     U1'*A1*U1 has the appropriate form. Let UT be an extension of U1 */
/*     to an N-by-N orthogonal matrix, using identity submatrices. Then A */
/*     is replaced by UT'*A*UT and the contents of array U is U * UT. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen, */
/*     March 1998. Based on the RASP routine SPLITB. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Eigenvalues, orthogonal transformation, real Schur form, */
/*     similarity transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

#line 122 "MB03QY.f"
    /* Parameter adjustments */
#line 122 "MB03QY.f"
    a_dim1 = *lda;
#line 122 "MB03QY.f"
    a_offset = 1 + a_dim1;
#line 122 "MB03QY.f"
    a -= a_offset;
#line 122 "MB03QY.f"
    u_dim1 = *ldu;
#line 122 "MB03QY.f"
    u_offset = 1 + u_dim1;
#line 122 "MB03QY.f"
    u -= u_offset;
#line 122 "MB03QY.f"

#line 122 "MB03QY.f"
    /* Function Body */
#line 122 "MB03QY.f"
    *info = 0;

/*     Test the input scalar arguments. */

#line 126 "MB03QY.f"
    if (*n < 2) {
#line 127 "MB03QY.f"
	*info = -1;
#line 128 "MB03QY.f"
    } else if (*l < 1 || *l >= *n) {
#line 129 "MB03QY.f"
	*info = -2;
#line 130 "MB03QY.f"
    } else if (*lda < *n) {
#line 131 "MB03QY.f"
	*info = -4;
#line 132 "MB03QY.f"
    } else if (*ldu < *n) {
#line 133 "MB03QY.f"
	*info = -6;
#line 134 "MB03QY.f"
    }

#line 136 "MB03QY.f"
    if (*info != 0) {

/*        Error return. */

#line 140 "MB03QY.f"
	i__1 = -(*info);
#line 140 "MB03QY.f"
	xerbla_("MB03QY", &i__1, (ftnlen)6);
#line 141 "MB03QY.f"
	return 0;
#line 142 "MB03QY.f"
    }

/*     Compute the eigenvalues and the elements of the Givens */
/*     transformation. */

#line 147 "MB03QY.f"
    l1 = *l + 1;
#line 148 "MB03QY.f"
    dlanv2_(&a[*l + *l * a_dim1], &a[*l + l1 * a_dim1], &a[l1 + *l * a_dim1], 
	    &a[l1 + l1 * a_dim1], e1, e2, &ew1, &ew2, &cs, &sn);
#line 150 "MB03QY.f"
    if (*e2 == 0.) {
#line 150 "MB03QY.f"
	*e2 = ew1;
#line 150 "MB03QY.f"
    }

/*     Apply the transformation to A. */

#line 154 "MB03QY.f"
    if (l1 < *n) {
#line 154 "MB03QY.f"
	i__1 = *n - l1;
#line 154 "MB03QY.f"
	drot_(&i__1, &a[*l + (l1 + 1) * a_dim1], lda, &a[l1 + (l1 + 1) * 
		a_dim1], lda, &cs, &sn);
#line 154 "MB03QY.f"
    }
#line 156 "MB03QY.f"
    i__1 = *l - 1;
#line 156 "MB03QY.f"
    drot_(&i__1, &a[*l * a_dim1 + 1], &c__1, &a[l1 * a_dim1 + 1], &c__1, &cs, 
	    &sn);

/*     Accumulate the transformation in U. */

#line 160 "MB03QY.f"
    drot_(n, &u[*l * u_dim1 + 1], &c__1, &u[l1 * u_dim1 + 1], &c__1, &cs, &sn)
	    ;

#line 162 "MB03QY.f"
    return 0;
/* *** Last line of MB03QY *** */
} /* mb03qy_ */
