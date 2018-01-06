#line 1 "MB03NY.f"
/* MB03NY.f -- translated by f2c (version 20100827).
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

#line 1 "MB03NY.f"
/* Table of constant values */

static integer c__1 = 1;

doublereal mb03ny_(integer *n, doublereal *omega, doublereal *a, integer *lda,
	 doublereal *s, doublereal *dwork, integer *ldwork, doublecomplex *
	cwork, integer *lcwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__, j, ic;
    static doublereal dummy[1]	/* was [1][1] */;
    extern /* Subroutine */ int dgesvd_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen), zgesvd_(char 
	    *, char *, integer *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    static doublecomplex zdummy[1]	/* was [1][1] */;


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

/*     To compute the smallest singular value of A - jwI. */

/*     FUNCTION VALUE */

/*     MB03NY  DOUBLE PRECISION */
/*             The smallest singular value of A - jwI (if INFO = 0). */
/*             If N = 0, the function value is set to zero. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the the matrix A.  N >= 0. */

/*     OMEGA   (input) DOUBLE PRECISION */
/*             The constant factor of A - jwI. */

/*     A       (input/workspace) DOUBLE PRECISION array, dimension */
/*             (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A. */
/*             On exit, if OMEGA = 0, the contents of this array are */
/*             destroyed. Otherwise, this array is unchanged. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     S       (output) DOUBLE PRECISION array, dimension (N) */
/*             The singular values of A - jwI in decreasing order. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX( 1, 5*N ). */
/*             For optimum performance LDWORK should be larger. */

/*     CWORK   COMPLEX*16 array, dimension (LCWORK) */
/*             On exit, if INFO = 0 and OMEGA <> 0, CWORK(1) returns the */
/*             optimal value of LCWORK. */
/*             If OMEGA is zero, this array is not referenced. */

/*     LCWORK  INTEGER */
/*             The length of the array CWORK. */
/*             LCWORK >= 1,                 if OMEGA =  0; */
/*             LCWORK >= MAX( 1, N*N+3*N ), if OMEGA <> 0. */
/*             For optimum performance LCWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 2:  The SVD algorithm (in either LAPACK Library routine */
/*                   DGESVD or ZGESVD) fails to converge; this error is */
/*                   very rare. */

/*     METHOD */

/*     This procedure simply constructs the matrix A - jwI, and calls */
/*     ZGESVD if w is not zero, or DGESVD if w = 0. */

/*     FURTHER COMMENTS */

/*     This routine is not very efficient because it computes all */
/*     singular values, but it is very accurate. The routine is intended */
/*     to be called only from the SLICOT Library routine AB13FD. */

/*     CONTRIBUTOR */

/*     R. Byers, the routine SIGMIN (January, 1995). */

/*     REVISIONS */

/*     Release 4.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1999. */

/*     REVISIONS */

/*     Oct. 2001, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Apr. 2002, V. Sima. */

/*     KEYWORDS */

/*     singular values. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

#line 141 "MB03NY.f"
    /* Parameter adjustments */
#line 141 "MB03NY.f"
    a_dim1 = *lda;
#line 141 "MB03NY.f"
    a_offset = 1 + a_dim1;
#line 141 "MB03NY.f"
    a -= a_offset;
#line 141 "MB03NY.f"
    --s;
#line 141 "MB03NY.f"
    --dwork;
#line 141 "MB03NY.f"
    --cwork;
#line 141 "MB03NY.f"

#line 141 "MB03NY.f"
    /* Function Body */
#line 141 "MB03NY.f"
    *info = 0;

#line 143 "MB03NY.f"
    if (*n < 0) {
#line 144 "MB03NY.f"
	*info = -1;
#line 145 "MB03NY.f"
    } else if (*lda < max(1,*n)) {
#line 146 "MB03NY.f"
	*info = -4;
#line 147 "MB03NY.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 147 "MB03NY.f"
	i__1 = 1, i__2 = *n * 5;
#line 147 "MB03NY.f"
	if (*ldwork < max(i__1,i__2)) {
#line 148 "MB03NY.f"
	    *info = -7;
#line 149 "MB03NY.f"
	} else if (*lcwork < 1 || *omega != 0. && *lcwork < *n * *n + *n * 3) 
		{
#line 151 "MB03NY.f"
	    *info = -9;
#line 152 "MB03NY.f"
	}
#line 152 "MB03NY.f"
    }

#line 154 "MB03NY.f"
    if (*info != 0) {

/*        Error return. */

#line 158 "MB03NY.f"
	i__1 = -(*info);
#line 158 "MB03NY.f"
	xerbla_("MB03NY", &i__1, (ftnlen)6);
#line 159 "MB03NY.f"
	return ret_val;
#line 160 "MB03NY.f"
    }

/*     Quick return if possible. */

#line 164 "MB03NY.f"
    if (*n == 0) {
#line 165 "MB03NY.f"
	ret_val = 0.;
#line 166 "MB03NY.f"
	dwork[1] = 1.;
#line 167 "MB03NY.f"
	if (*omega != 0.) {
#line 167 "MB03NY.f"
	    cwork[1].r = 1., cwork[1].i = 0.;
#line 167 "MB03NY.f"
	}
#line 169 "MB03NY.f"
	return ret_val;
#line 170 "MB03NY.f"
    }

#line 172 "MB03NY.f"
    if (*omega == 0.) {

/*        OMEGA = 0 allows real SVD. */

#line 176 "MB03NY.f"
	dgesvd_("No vectors", "No vectors", n, n, &a[a_offset], n, &s[1], 
		dummy, &c__1, dummy, &c__1, &dwork[1], ldwork, info, (ftnlen)
		10, (ftnlen)10);
#line 178 "MB03NY.f"
	if (*info != 0) {
#line 179 "MB03NY.f"
	    *info = 2;
#line 180 "MB03NY.f"
	    return ret_val;
#line 181 "MB03NY.f"
	}
#line 182 "MB03NY.f"
    } else {

/*        General case, that is complex SVD. */

#line 186 "MB03NY.f"
	ic = 1;
#line 187 "MB03NY.f"
	i__1 = *n;
#line 187 "MB03NY.f"
	for (j = 1; j <= i__1; ++j) {
#line 188 "MB03NY.f"
	    i__2 = *n;
#line 188 "MB03NY.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 189 "MB03NY.f"
		i__3 = ic;
#line 189 "MB03NY.f"
		i__4 = i__ + j * a_dim1;
#line 189 "MB03NY.f"
		cwork[i__3].r = a[i__4], cwork[i__3].i = 0.;
#line 190 "MB03NY.f"
		++ic;
#line 191 "MB03NY.f"
/* L10: */
#line 191 "MB03NY.f"
	    }
#line 192 "MB03NY.f"
	    i__2 = (j - 1) * *n + j;
#line 192 "MB03NY.f"
	    i__3 = (j - 1) * *n + j;
#line 192 "MB03NY.f"
	    z__2.r = *omega * 0., z__2.i = *omega * 1.;
#line 192 "MB03NY.f"
	    z__1.r = cwork[i__3].r - z__2.r, z__1.i = cwork[i__3].i - z__2.i;
#line 192 "MB03NY.f"
	    cwork[i__2].r = z__1.r, cwork[i__2].i = z__1.i;
#line 193 "MB03NY.f"
/* L20: */
#line 193 "MB03NY.f"
	}
#line 194 "MB03NY.f"
	i__1 = *lcwork - *n * *n;
#line 194 "MB03NY.f"
	zgesvd_("No vectors", "No vectors", n, n, &cwork[1], n, &s[1], zdummy,
		 &c__1, zdummy, &c__1, &cwork[*n * *n + 1], &i__1, &dwork[1], 
		info, (ftnlen)10, (ftnlen)10);
#line 197 "MB03NY.f"
	if (*info != 0) {
#line 198 "MB03NY.f"
	    *info = 2;
#line 199 "MB03NY.f"
	    return ret_val;
#line 200 "MB03NY.f"
	}
#line 201 "MB03NY.f"
	i__1 = *n * *n + 1;
#line 201 "MB03NY.f"
	d__1 = (doublereal) (*n * *n);
#line 201 "MB03NY.f"
	z__2.r = d__1 * 1., z__2.i = d__1 * 0.;
#line 201 "MB03NY.f"
	z__1.r = cwork[i__1].r + z__2.r, z__1.i = cwork[i__1].i + z__2.i;
#line 201 "MB03NY.f"
	cwork[1].r = z__1.r, cwork[1].i = z__1.i;
#line 202 "MB03NY.f"
	dwork[1] = (doublereal) (*n * 5);
#line 203 "MB03NY.f"
    }

#line 205 "MB03NY.f"
    ret_val = s[*n];

/* *** Last line of MB03NY *** */
#line 208 "MB03NY.f"
    return ret_val;
} /* mb03ny_ */

