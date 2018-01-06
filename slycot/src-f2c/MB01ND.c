#line 1 "MB01ND.f"
/* MB01ND.f -- translated by f2c (version 20100827).
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

#line 1 "MB01ND.f"
/* Subroutine */ int mb01nd_(char *uplo, integer *n, doublereal *alpha, 
	doublereal *x, integer *incx, doublereal *y, integer *incy, 
	doublereal *a, integer *lda, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ix, iy, jx, jy, kx, ky, info;
    static doublereal temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

/*     To perform the skew-symmetric rank 2 operation */

/*          A := alpha*x*y' - alpha*y*x' + A, */

/*     where alpha is a scalar, x and y are vectors of length n and A is */
/*     an n-by-n skew-symmetric matrix. */

/*     This is a modified version of the vanilla implemented BLAS */
/*     routine DSYR2 written by Jack Dongarra, Jeremy Du Croz, */
/*     Sven Hammarling, and Richard Hanson. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPLO    CHARACTER*1 */
/*             Specifies whether the upper or lower triangular part of */
/*             the array A is to be referenced as follows: */
/*             = 'U':  only the strictly upper triangular part of A is to */
/*                     be referenced; */
/*             = 'L':  only the strictly lower triangular part of A is to */
/*                     be referenced. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar alpha. If alpha is zero X and Y are not */
/*             referenced. */

/*     X       (input) DOUBLE PRECISION array, dimension */
/*             ( 1 + ( N - 1 )*abs( INCX ) ). */
/*             On entry, elements 1, INCX+1, .., ( N - 1 )*INCX + 1 of */
/*             this array must contain the elements of the vector X. */

/*     INCX    (input) INTEGER */
/*             The increment for the elements of X. IF INCX < 0 then the */
/*             elements of X are accessed in reversed order.  INCX <> 0. */

/*     Y       (input) DOUBLE PRECISION array, dimension */
/*             ( 1 + ( N - 1 )*abs( INCY ) ). */
/*             On entry, elements 1, INCY+1, .., ( N - 1 )*INCY + 1 of */
/*             this array must contain the elements of the vector Y. */

/*     INCY    (input) INTEGER */
/*             The increment for the elements of Y. IF INCY < 0 then the */
/*             elements of Y are accessed in reversed order.  INCY <> 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry with UPLO = 'U', the leading N-by-N part of this */
/*             array must contain the strictly upper triangular part of */
/*             the matrix A. The lower triangular part of this array is */
/*             not referenced. */
/*             On entry with UPLO = 'L', the leading N-by-N part of this */
/*             array must contain the strictly lower triangular part of */
/*             the matrix A. The upper triangular part of this array is */
/*             not referenced. */
/*             On exit with UPLO = 'U', the leading N-by-N part of this */
/*             array contains the strictly upper triangular part of the */
/*             updated matrix A. */
/*             On exit with UPLO = 'L', the leading N-by-N part of this */
/*             array contains the strictly lower triangular part of the */
/*             updated matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N) */

/*     NUMERICAL ASPECTS */

/*     Though being almost identical with the vanilla implementation */
/*     of the BLAS routine DSYR2 the performance of this routine could */
/*     be significantly lower in the case of vendor supplied, highly */
/*     optimized BLAS. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, May 2008 (SLICOT version of the HAPACK routine DSKR2). */

/*     KEYWORDS */

/*     Elementary matrix operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 138 "MB01ND.f"
    /* Parameter adjustments */
#line 138 "MB01ND.f"
    --x;
#line 138 "MB01ND.f"
    --y;
#line 138 "MB01ND.f"
    a_dim1 = *lda;
#line 138 "MB01ND.f"
    a_offset = 1 + a_dim1;
#line 138 "MB01ND.f"
    a -= a_offset;
#line 138 "MB01ND.f"

#line 138 "MB01ND.f"
    /* Function Body */
#line 138 "MB01ND.f"
    info = 0;
#line 139 "MB01ND.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 141 "MB01ND.f"
	info = 1;
#line 142 "MB01ND.f"
    } else if (*n < 0) {
#line 143 "MB01ND.f"
	info = 2;
#line 144 "MB01ND.f"
    } else if (*incx == 0) {
#line 145 "MB01ND.f"
	info = 5;
#line 146 "MB01ND.f"
    } else if (*incy == 0) {
#line 147 "MB01ND.f"
	info = 7;
#line 148 "MB01ND.f"
    } else if (*lda < max(1,*n)) {
#line 149 "MB01ND.f"
	info = 9;
#line 150 "MB01ND.f"
    }

#line 152 "MB01ND.f"
    if (info != 0) {
#line 153 "MB01ND.f"
	xerbla_("MB01ND", &info, (ftnlen)6);
#line 154 "MB01ND.f"
	return 0;
#line 155 "MB01ND.f"
    }

/*     Quick return if possible. */

#line 159 "MB01ND.f"
    if (*n == 0 || *alpha == 0.) {
#line 159 "MB01ND.f"
	return 0;
#line 159 "MB01ND.f"
    }

/*     Set up the start points in X and Y if the increments are not both */
/*     unity. */

#line 165 "MB01ND.f"
    if (*incx != 1 || *incy != 1) {
#line 166 "MB01ND.f"
	if (*incx > 0) {
#line 167 "MB01ND.f"
	    kx = 1;
#line 168 "MB01ND.f"
	} else {
#line 169 "MB01ND.f"
	    kx = 1 - (*n - 1) * *incx;
#line 170 "MB01ND.f"
	}
#line 171 "MB01ND.f"
	if (*incy > 0) {
#line 172 "MB01ND.f"
	    ky = 1;
#line 173 "MB01ND.f"
	} else {
#line 174 "MB01ND.f"
	    ky = 1 - (*n - 1) * *incy;
#line 175 "MB01ND.f"
	}
#line 176 "MB01ND.f"
	jx = kx;
#line 177 "MB01ND.f"
	jy = ky;
#line 178 "MB01ND.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

#line 184 "MB01ND.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form A when A is stored in the upper triangle. */

#line 188 "MB01ND.f"
	if (*incx == 1 && *incy == 1) {
#line 189 "MB01ND.f"
	    i__1 = *n;
#line 189 "MB01ND.f"
	    for (j = 2; j <= i__1; ++j) {
#line 190 "MB01ND.f"
		if (x[j] != 0. || y[j] != 0.) {
#line 191 "MB01ND.f"
		    temp1 = *alpha * y[j];
#line 192 "MB01ND.f"
		    temp2 = *alpha * x[j];
#line 193 "MB01ND.f"
		    i__2 = j - 1;
#line 193 "MB01ND.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 194 "MB01ND.f"
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[i__] * 
				temp1 - y[i__] * temp2;
#line 195 "MB01ND.f"
/* L10: */
#line 195 "MB01ND.f"
		    }
#line 196 "MB01ND.f"
		}
#line 197 "MB01ND.f"
/* L20: */
#line 197 "MB01ND.f"
	    }
#line 198 "MB01ND.f"
	} else {
#line 199 "MB01ND.f"
	    i__1 = *n;
#line 199 "MB01ND.f"
	    for (j = 2; j <= i__1; ++j) {
#line 200 "MB01ND.f"
		if (x[jx] != 0. || y[jy] != 0.) {
#line 201 "MB01ND.f"
		    temp1 = *alpha * y[jy];
#line 202 "MB01ND.f"
		    temp2 = *alpha * x[jx];
#line 203 "MB01ND.f"
		    ix = kx;
#line 204 "MB01ND.f"
		    iy = ky;
#line 205 "MB01ND.f"
		    i__2 = j - 1;
#line 205 "MB01ND.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 206 "MB01ND.f"
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[ix] * 
				temp1 - y[iy] * temp2;
#line 207 "MB01ND.f"
			ix += *incx;
#line 208 "MB01ND.f"
			iy += *incy;
#line 209 "MB01ND.f"
/* L30: */
#line 209 "MB01ND.f"
		    }
#line 210 "MB01ND.f"
		}
#line 211 "MB01ND.f"
		jx += *incx;
#line 212 "MB01ND.f"
		jy += *incy;
#line 213 "MB01ND.f"
/* L40: */
#line 213 "MB01ND.f"
	    }
#line 214 "MB01ND.f"
	}
#line 215 "MB01ND.f"
    } else {

/*        Form A when A is stored in the lower triangle. */

#line 219 "MB01ND.f"
	if (*incx == 1 && *incy == 1) {
#line 220 "MB01ND.f"
	    i__1 = *n - 1;
#line 220 "MB01ND.f"
	    for (j = 1; j <= i__1; ++j) {
#line 221 "MB01ND.f"
		if (x[j] != 0. || y[j] != 0.) {
#line 222 "MB01ND.f"
		    temp1 = *alpha * y[j];
#line 223 "MB01ND.f"
		    temp2 = *alpha * x[j];
#line 224 "MB01ND.f"
		    i__2 = *n;
#line 224 "MB01ND.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 225 "MB01ND.f"
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[i__] * 
				temp1 - y[i__] * temp2;
#line 226 "MB01ND.f"
/* L50: */
#line 226 "MB01ND.f"
		    }
#line 227 "MB01ND.f"
		}
#line 228 "MB01ND.f"
/* L60: */
#line 228 "MB01ND.f"
	    }
#line 229 "MB01ND.f"
	} else {
#line 230 "MB01ND.f"
	    i__1 = *n - 1;
#line 230 "MB01ND.f"
	    for (j = 1; j <= i__1; ++j) {
#line 231 "MB01ND.f"
		if (x[jx] != 0. || y[jy] != 0.) {
#line 232 "MB01ND.f"
		    temp1 = *alpha * y[jy];
#line 233 "MB01ND.f"
		    temp2 = *alpha * x[jx];
#line 234 "MB01ND.f"
		    ix = jx;
#line 235 "MB01ND.f"
		    iy = jy;
#line 236 "MB01ND.f"
		    i__2 = *n;
#line 236 "MB01ND.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 237 "MB01ND.f"
			a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[ix] * 
				temp1 - y[iy] * temp2;
#line 238 "MB01ND.f"
			ix += *incx;
#line 239 "MB01ND.f"
			iy += *incy;
#line 240 "MB01ND.f"
/* L70: */
#line 240 "MB01ND.f"
		    }
#line 241 "MB01ND.f"
		}
#line 242 "MB01ND.f"
		jx += *incx;
#line 243 "MB01ND.f"
		jy += *incy;
#line 244 "MB01ND.f"
/* L80: */
#line 244 "MB01ND.f"
	    }
#line 245 "MB01ND.f"
	}
#line 246 "MB01ND.f"
    }
#line 247 "MB01ND.f"
    return 0;
/* *** Last line of MB01ND *** */
} /* mb01nd_ */

