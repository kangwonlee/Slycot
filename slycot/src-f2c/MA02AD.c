#line 1 "MA02AD.f"
/* MA02AD.f -- translated by f2c (version 20100827).
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

#line 1 "MA02AD.f"
/* Subroutine */ int ma02ad_(char *job, integer *m, integer *n, doublereal *a,
	 integer *lda, doublereal *b, integer *ldb, ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
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

/*     To transpose all or part of a two-dimensional matrix A into */
/*     another matrix B. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the part of the matrix A to be transposed into B */
/*             as follows: */
/*             = 'U': Upper triangular part; */
/*             = 'L': Lower triangular part; */
/*             Otherwise:  All of the matrix A. */

/*     Input/Output Parameters */

/*     M      (input) INTEGER */
/*            The number of rows of the matrix A.  M >= 0. */

/*     N      (input) INTEGER */
/*            The number of columns of the matrix A.  N >= 0. */

/*     A      (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*            The m-by-n matrix A.  If JOB = 'U', only the upper */
/*            triangle or trapezoid is accessed; if JOB = 'L', only the */
/*            lower triangle or trapezoid is accessed. */

/*     LDA    INTEGER */
/*            The leading dimension of the array A.  LDA >= max(1,M). */

/*     B      (output) DOUBLE PRECISION array, dimension (LDB,M) */
/*            B = A' in the locations specified by JOB. */

/*     LDB    INTEGER */
/*            The leading dimension of the array B.  LDB >= max(1,N). */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, March 1998. */
/*     Based on the RASP routine DMTRA. */

/*     REVISIONS */

/*     - */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

#line 86 "MA02AD.f"
    /* Parameter adjustments */
#line 86 "MA02AD.f"
    a_dim1 = *lda;
#line 86 "MA02AD.f"
    a_offset = 1 + a_dim1;
#line 86 "MA02AD.f"
    a -= a_offset;
#line 86 "MA02AD.f"
    b_dim1 = *ldb;
#line 86 "MA02AD.f"
    b_offset = 1 + b_dim1;
#line 86 "MA02AD.f"
    b -= b_offset;
#line 86 "MA02AD.f"

#line 86 "MA02AD.f"
    /* Function Body */
#line 86 "MA02AD.f"
    if (lsame_(job, "U", (ftnlen)1, (ftnlen)1)) {
#line 87 "MA02AD.f"
	i__1 = *n;
#line 87 "MA02AD.f"
	for (j = 1; j <= i__1; ++j) {
#line 88 "MA02AD.f"
	    i__2 = min(j,*m);
#line 88 "MA02AD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 89 "MA02AD.f"
		b[j + i__ * b_dim1] = a[i__ + j * a_dim1];
#line 90 "MA02AD.f"
/* L10: */
#line 90 "MA02AD.f"
	    }
#line 91 "MA02AD.f"
/* L20: */
#line 91 "MA02AD.f"
	}
#line 92 "MA02AD.f"
    } else if (lsame_(job, "L", (ftnlen)1, (ftnlen)1)) {
#line 93 "MA02AD.f"
	i__1 = *n;
#line 93 "MA02AD.f"
	for (j = 1; j <= i__1; ++j) {
#line 94 "MA02AD.f"
	    i__2 = *m;
#line 94 "MA02AD.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
#line 95 "MA02AD.f"
		b[j + i__ * b_dim1] = a[i__ + j * a_dim1];
#line 96 "MA02AD.f"
/* L30: */
#line 96 "MA02AD.f"
	    }
#line 97 "MA02AD.f"
/* L40: */
#line 97 "MA02AD.f"
	}
#line 98 "MA02AD.f"
    } else {
#line 99 "MA02AD.f"
	i__1 = *n;
#line 99 "MA02AD.f"
	for (j = 1; j <= i__1; ++j) {
#line 100 "MA02AD.f"
	    i__2 = *m;
#line 100 "MA02AD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 101 "MA02AD.f"
		b[j + i__ * b_dim1] = a[i__ + j * a_dim1];
#line 102 "MA02AD.f"
/* L50: */
#line 102 "MA02AD.f"
	    }
#line 103 "MA02AD.f"
/* L60: */
#line 103 "MA02AD.f"
	}
#line 104 "MA02AD.f"
    }

#line 106 "MA02AD.f"
    return 0;
/* *** Last line of MA02AD *** */
} /* ma02ad_ */

