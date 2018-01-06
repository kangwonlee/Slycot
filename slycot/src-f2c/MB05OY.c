#line 1 "MB05OY.f"
/* MB05OY.f -- translated by f2c (version 20100827).
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

#line 1 "MB05OY.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mb05oy_(char *job, integer *n, integer *low, integer *
	igh, doublereal *a, integer *lda, doublereal *scale, integer *info, 
	ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, ii;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
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

/*     To restore a matrix after it has been transformed by applying */
/*     balancing transformations (permutations and scalings), as */
/*     determined by LAPACK Library routine DGEBAL. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the type of backward transformation required, */
/*             as follows: */
/*             = 'N', do nothing, return immediately; */
/*             = 'P', do backward transformation for permutation only; */
/*             = 'S', do backward transformation for scaling only; */
/*             = 'B', do backward transformations for both permutation */
/*                    and scaling. */
/*             JOB must be the same as the argument JOB supplied */
/*             to DGEBAL. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     LOW     (input) INTEGER */
/*     IGH     (input) INTEGER */
/*             The integers LOW and IGH determined by DGEBAL. */
/*             1 <= LOW <= IGH <= N, if N > 0; LOW=1 and IGH=0, if N=0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix to be back-transformed. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the transformed matrix. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     SCALE   (input) DOUBLE PRECISION array, dimension (N) */
/*             Details of the permutation and scaling factors, as */
/*             returned by DGEBAL. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Let P be a permutation matrix, and D a diagonal matrix of scaling */
/*     factors, both of order N. The routine computes */
/*                     -1 */
/*        A <-- P D A D  P'. */

/*     where the permutation and scaling factors are encoded in the */
/*     array SCALE. */

/*     REFERENCES */

/*     None. */

/*     NUMERICAL ASPECTS */
/*                               2 */
/*     The algorithm requires O(N ) operations. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997. */
/*     Supersedes Release 2.0 routine MB05CY. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix algebra, matrix operations. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

#line 125 "MB05OY.f"
    /* Parameter adjustments */
#line 125 "MB05OY.f"
    a_dim1 = *lda;
#line 125 "MB05OY.f"
    a_offset = 1 + a_dim1;
#line 125 "MB05OY.f"
    a -= a_offset;
#line 125 "MB05OY.f"
    --scale;
#line 125 "MB05OY.f"

#line 125 "MB05OY.f"
    /* Function Body */
#line 125 "MB05OY.f"
    *info = 0;
#line 126 "MB05OY.f"
    if (! lsame_(job, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(job, "P", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(job, "S", (ftnlen)1, (ftnlen)1) 
	    && ! lsame_(job, "B", (ftnlen)1, (ftnlen)1)) {
#line 128 "MB05OY.f"
	*info = -1;
#line 129 "MB05OY.f"
    } else if (*n < 0) {
#line 130 "MB05OY.f"
	*info = -2;
#line 131 "MB05OY.f"
    } else if (*low < 1 || *low > max(1,*n)) {
#line 132 "MB05OY.f"
	*info = -3;
#line 133 "MB05OY.f"
    } else if (*igh < min(*low,*n) || *igh > *n) {
#line 134 "MB05OY.f"
	*info = -4;
#line 135 "MB05OY.f"
    } else if (*lda < max(1,*n)) {
#line 136 "MB05OY.f"
	*info = -6;
#line 137 "MB05OY.f"
    }

#line 139 "MB05OY.f"
    if (*info != 0) {
#line 140 "MB05OY.f"
	i__1 = -(*info);
#line 140 "MB05OY.f"
	xerbla_("MB05OY", &i__1, (ftnlen)6);
#line 141 "MB05OY.f"
	return 0;
#line 142 "MB05OY.f"
    }

/*     Quick return if possible. */

#line 146 "MB05OY.f"
    if (*n == 0 || lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 146 "MB05OY.f"
	return 0;
#line 146 "MB05OY.f"
    }

#line 149 "MB05OY.f"
    if (! lsame_(job, "P", (ftnlen)1, (ftnlen)1) && *igh != *low) {

#line 151 "MB05OY.f"
	i__1 = *igh;
#line 151 "MB05OY.f"
	for (i__ = *low; i__ <= i__1; ++i__) {
#line 152 "MB05OY.f"
	    dscal_(n, &scale[i__], &a[i__ + a_dim1], lda);
#line 153 "MB05OY.f"
/* L20: */
#line 153 "MB05OY.f"
	}

#line 155 "MB05OY.f"
	i__1 = *igh;
#line 155 "MB05OY.f"
	for (j = *low; j <= i__1; ++j) {
#line 156 "MB05OY.f"
	    d__1 = 1. / scale[j];
#line 156 "MB05OY.f"
	    dscal_(n, &d__1, &a[j * a_dim1 + 1], &c__1);
#line 157 "MB05OY.f"
/* L40: */
#line 157 "MB05OY.f"
	}

#line 159 "MB05OY.f"
    }

#line 161 "MB05OY.f"
    if (! lsame_(job, "S", (ftnlen)1, (ftnlen)1)) {

#line 163 "MB05OY.f"
	i__1 = *n;
#line 163 "MB05OY.f"
	for (ii = 1; ii <= i__1; ++ii) {
#line 164 "MB05OY.f"
	    i__ = ii;
#line 165 "MB05OY.f"
	    if (i__ < *low || i__ > *igh) {
#line 166 "MB05OY.f"
		if (i__ < *low) {
#line 166 "MB05OY.f"
		    i__ = *low - ii;
#line 166 "MB05OY.f"
		}
#line 167 "MB05OY.f"
		k = (integer) scale[i__];
#line 168 "MB05OY.f"
		if (k != i__) {
#line 169 "MB05OY.f"
		    dswap_(n, &a[i__ + a_dim1], lda, &a[k + a_dim1], lda);
#line 170 "MB05OY.f"
		    dswap_(n, &a[i__ * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
			     &c__1);
#line 171 "MB05OY.f"
		}
#line 172 "MB05OY.f"
	    }
#line 173 "MB05OY.f"
/* L60: */
#line 173 "MB05OY.f"
	}

#line 175 "MB05OY.f"
    }

#line 177 "MB05OY.f"
    return 0;
/* *** Last line of MB05OY *** */
} /* mb05oy_ */

