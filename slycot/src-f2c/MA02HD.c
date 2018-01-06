#line 1 "MA02HD.f"
/* MA02HD.f -- translated by f2c (version 20100827).
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

#line 1 "MA02HD.f"
logical ma02hd_(char *job, integer *m, integer *n, doublereal *diag, 
	doublereal *a, integer *lda, ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    logical ret_val;

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

/*     To check if A = DIAG*I, where I is an M-by-N matrix with ones on */
/*     the diagonal and zeros elsewhere. */

/*     FUNCTION VALUE */

/*     MA02HD  LOGICAL */
/*             The function value is set to .TRUE. if A = DIAG*I, and to */
/*             .FALSE., otherwise. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the part of the matrix A to be checked out, */
/*             as follows: */
/*             = 'U': Upper triangular/trapezoidal part; */
/*             = 'L': Lower triangular/trapezoidal part. */
/*             Otherwise:  All of the matrix A. */

/*     Input/Output Parameters */

/*     M      (input) INTEGER */
/*            The number of rows of the matrix A.  M >= 0. */

/*     N      (input) INTEGER */
/*            The number of columns of the matrix A.  N >= 0. */

/*     DIAG   (input) DOUBLE PRECISION */
/*            The scalar DIAG. */

/*     A      (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*            The leading M-by-N part of this array must contain the */
/*            matrix A.  If JOB = 'U', only the upper triangle or */
/*            trapezoid is accessed; if JOB = 'L', only the lower */
/*            triangle or trapezoid is accessed. */

/*     LDA    INTEGER */
/*            The leading dimension of the array A.  LDA >= max(1,M). */

/*     METHOD */

/*     The routine returns immediately after detecting a diagonal element */
/*     which differs from DIAG, or a nonzero off-diagonal element in the */
/*     searched part of A. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 2001. */
/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Jan. 2003. */

/*     KEYWORDS */

/*     Elementary operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Do not check parameters, for efficiency. */

#line 105 "MA02HD.f"
    /* Parameter adjustments */
#line 105 "MA02HD.f"
    a_dim1 = *lda;
#line 105 "MA02HD.f"
    a_offset = 1 + a_dim1;
#line 105 "MA02HD.f"
    a -= a_offset;
#line 105 "MA02HD.f"

#line 105 "MA02HD.f"
    /* Function Body */
#line 105 "MA02HD.f"
    if (lsame_(job, "U", (ftnlen)1, (ftnlen)1)) {

#line 107 "MA02HD.f"
	i__1 = *n;
#line 107 "MA02HD.f"
	for (j = 1; j <= i__1; ++j) {

/* Computing MIN */
#line 109 "MA02HD.f"
	    i__3 = j - 1;
#line 109 "MA02HD.f"
	    i__2 = min(i__3,*m);
#line 109 "MA02HD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 110 "MA02HD.f"
		if (a[i__ + j * a_dim1] != 0.) {
#line 111 "MA02HD.f"
		    ret_val = FALSE_;
#line 112 "MA02HD.f"
		    return ret_val;
#line 113 "MA02HD.f"
		}
#line 114 "MA02HD.f"
/* L10: */
#line 114 "MA02HD.f"
	    }

#line 116 "MA02HD.f"
	    if (j <= *m) {
#line 117 "MA02HD.f"
		if (a[j + j * a_dim1] != *diag) {
#line 118 "MA02HD.f"
		    ret_val = FALSE_;
#line 119 "MA02HD.f"
		    return ret_val;
#line 120 "MA02HD.f"
		}
#line 121 "MA02HD.f"
	    }
#line 122 "MA02HD.f"
/* L20: */
#line 122 "MA02HD.f"
	}

#line 124 "MA02HD.f"
    } else if (lsame_(job, "L", (ftnlen)1, (ftnlen)1)) {

#line 126 "MA02HD.f"
	i__1 = min(*m,*n);
#line 126 "MA02HD.f"
	for (j = 1; j <= i__1; ++j) {
#line 127 "MA02HD.f"
	    if (a[j + j * a_dim1] != *diag) {
#line 128 "MA02HD.f"
		ret_val = FALSE_;
#line 129 "MA02HD.f"
		return ret_val;
#line 130 "MA02HD.f"
	    }

#line 132 "MA02HD.f"
	    if (j != *m) {

/* Computing MIN */
#line 134 "MA02HD.f"
		i__2 = j + 1;
#line 134 "MA02HD.f"
		i__3 = *m;
#line 134 "MA02HD.f"
		for (i__ = min(i__2,*m); i__ <= i__3; ++i__) {
#line 135 "MA02HD.f"
		    if (a[i__ + j * a_dim1] != 0.) {
#line 136 "MA02HD.f"
			ret_val = FALSE_;
#line 137 "MA02HD.f"
			return ret_val;
#line 138 "MA02HD.f"
		    }
#line 139 "MA02HD.f"
/* L30: */
#line 139 "MA02HD.f"
		}

#line 141 "MA02HD.f"
	    }
#line 142 "MA02HD.f"
/* L40: */
#line 142 "MA02HD.f"
	}

#line 144 "MA02HD.f"
    } else {

#line 146 "MA02HD.f"
	i__1 = *n;
#line 146 "MA02HD.f"
	for (j = 1; j <= i__1; ++j) {

/* Computing MIN */
#line 148 "MA02HD.f"
	    i__2 = j - 1;
#line 148 "MA02HD.f"
	    i__3 = min(i__2,*m);
#line 148 "MA02HD.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 149 "MA02HD.f"
		if (a[i__ + j * a_dim1] != 0.) {
#line 150 "MA02HD.f"
		    ret_val = FALSE_;
#line 151 "MA02HD.f"
		    return ret_val;
#line 152 "MA02HD.f"
		}
#line 153 "MA02HD.f"
/* L50: */
#line 153 "MA02HD.f"
	    }

#line 155 "MA02HD.f"
	    if (j <= *m) {
#line 156 "MA02HD.f"
		if (a[j + j * a_dim1] != *diag) {
#line 157 "MA02HD.f"
		    ret_val = FALSE_;
#line 158 "MA02HD.f"
		    return ret_val;
#line 159 "MA02HD.f"
		}
#line 160 "MA02HD.f"
	    }

#line 162 "MA02HD.f"
	    if (j < *m) {

/* Computing MIN */
#line 164 "MA02HD.f"
		i__3 = j + 1;
#line 164 "MA02HD.f"
		i__2 = *m;
#line 164 "MA02HD.f"
		for (i__ = min(i__3,*m); i__ <= i__2; ++i__) {
#line 165 "MA02HD.f"
		    if (a[i__ + j * a_dim1] != 0.) {
#line 166 "MA02HD.f"
			ret_val = FALSE_;
#line 167 "MA02HD.f"
			return ret_val;
#line 168 "MA02HD.f"
		    }
#line 169 "MA02HD.f"
/* L60: */
#line 169 "MA02HD.f"
		}

#line 171 "MA02HD.f"
	    }
#line 172 "MA02HD.f"
/* L70: */
#line 172 "MA02HD.f"
	}

#line 174 "MA02HD.f"
    }

#line 176 "MA02HD.f"
    ret_val = TRUE_;

#line 178 "MA02HD.f"
    return ret_val;
/* *** Last line of MA02HD *** */
} /* ma02hd_ */

