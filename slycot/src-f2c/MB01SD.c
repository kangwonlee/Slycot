#line 1 "MB01SD.f"
/* MB01SD.f -- translated by f2c (version 20100827).
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

#line 1 "MB01SD.f"
/* Subroutine */ int mb01sd_(char *jobs, integer *m, integer *n, doublereal *
	a, integer *lda, doublereal *r__, doublereal *c__, ftnlen jobs_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal cj;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);


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

/*     To scale a general M-by-N matrix A using the row and column */
/*     scaling factors in the vectors R and C. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBS    CHARACTER*1 */
/*             Specifies the scaling operation to be done, as follows: */
/*             = 'R':  row scaling, i.e., A will be premultiplied */
/*                     by diag(R); */
/*             = 'C':  column scaling, i.e., A will be postmultiplied */
/*                     by diag(C); */
/*             = 'B':  both row and column scaling, i.e., A will be */
/*                     replaced by diag(R) * A * diag(C). */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the M-by-N matrix A. */
/*             On exit, the scaled matrix.  See JOBS for the form of the */
/*             scaled matrix. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     R       (input) DOUBLE PRECISION array, dimension (M) */
/*             The row scale factors for A. */
/*             R is not referenced if JOBS = 'C'. */

/*     C       (input) DOUBLE PRECISION array, dimension (N) */
/*             The column scale factors for A. */
/*             C is not referenced if JOBS = 'R'. */


/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, April 1998. */
/*     Based on the RASP routine DMSCAL. */

/*    ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. Executable Statements .. */

/*     Quick return if possible. */

#line 87 "MB01SD.f"
    /* Parameter adjustments */
#line 87 "MB01SD.f"
    a_dim1 = *lda;
#line 87 "MB01SD.f"
    a_offset = 1 + a_dim1;
#line 87 "MB01SD.f"
    a -= a_offset;
#line 87 "MB01SD.f"
    --r__;
#line 87 "MB01SD.f"
    --c__;
#line 87 "MB01SD.f"

#line 87 "MB01SD.f"
    /* Function Body */
#line 87 "MB01SD.f"
    if (*m == 0 || *n == 0) {
#line 87 "MB01SD.f"
	return 0;
#line 87 "MB01SD.f"
    }

#line 90 "MB01SD.f"
    if (lsame_(jobs, "C", (ftnlen)1, (ftnlen)1)) {

/*        Column scaling, no row scaling. */

#line 94 "MB01SD.f"
	i__1 = *n;
#line 94 "MB01SD.f"
	for (j = 1; j <= i__1; ++j) {
#line 95 "MB01SD.f"
	    cj = c__[j];
#line 96 "MB01SD.f"
	    i__2 = *m;
#line 96 "MB01SD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 97 "MB01SD.f"
		a[i__ + j * a_dim1] = cj * a[i__ + j * a_dim1];
#line 98 "MB01SD.f"
/* L10: */
#line 98 "MB01SD.f"
	    }
#line 99 "MB01SD.f"
/* L20: */
#line 99 "MB01SD.f"
	}
#line 100 "MB01SD.f"
    } else if (lsame_(jobs, "R", (ftnlen)1, (ftnlen)1)) {

/*        Row scaling, no column scaling. */

#line 104 "MB01SD.f"
	i__1 = *n;
#line 104 "MB01SD.f"
	for (j = 1; j <= i__1; ++j) {
#line 105 "MB01SD.f"
	    i__2 = *m;
#line 105 "MB01SD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 106 "MB01SD.f"
		a[i__ + j * a_dim1] = r__[i__] * a[i__ + j * a_dim1];
#line 107 "MB01SD.f"
/* L30: */
#line 107 "MB01SD.f"
	    }
#line 108 "MB01SD.f"
/* L40: */
#line 108 "MB01SD.f"
	}
#line 109 "MB01SD.f"
    } else if (lsame_(jobs, "B", (ftnlen)1, (ftnlen)1)) {

/*        Row and column scaling. */

#line 113 "MB01SD.f"
	i__1 = *n;
#line 113 "MB01SD.f"
	for (j = 1; j <= i__1; ++j) {
#line 114 "MB01SD.f"
	    cj = c__[j];
#line 115 "MB01SD.f"
	    i__2 = *m;
#line 115 "MB01SD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 116 "MB01SD.f"
		a[i__ + j * a_dim1] = cj * r__[i__] * a[i__ + j * a_dim1];
#line 117 "MB01SD.f"
/* L50: */
#line 117 "MB01SD.f"
	    }
#line 118 "MB01SD.f"
/* L60: */
#line 118 "MB01SD.f"
	}
#line 119 "MB01SD.f"
    }

#line 121 "MB01SD.f"
    return 0;
/* *** Last line of MB01SD *** */
} /* mb01sd_ */
