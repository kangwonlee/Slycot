#line 1 "MB01MD.f"
/* MB01MD.f -- translated by f2c (version 20100827).
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

#line 1 "MB01MD.f"
/* Subroutine */ int mb01md_(char *uplo, integer *n, doublereal *alpha, 
	doublereal *a, integer *lda, doublereal *x, integer *incx, doublereal 
	*beta, doublereal *y, integer *incy, ftnlen uplo_len)
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

/*     To perform the matrix-vector operation */

/*        y := alpha*A*x + beta*y, */

/*     where alpha and beta are scalars, x and y are vectors of length */
/*     n and A is an n-by-n skew-symmetric matrix. */

/*     This is a modified version of the vanilla implemented BLAS */
/*     routine DSYMV written by Jack Dongarra, Jeremy Du Croz, */
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
/*             The scalar alpha. If alpha is zero the array A is not */
/*             referenced. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry with UPLO = 'U', the leading N-by-N part of this */
/*             array must contain the strictly upper triangular part of */
/*             the matrix A. The lower triangular part of this array is */
/*             not referenced. */
/*             On entry with UPLO = 'L', the leading N-by-N part of this */
/*             array must contain the strictly lower triangular part of */
/*             the matrix A. The upper triangular part of this array is */
/*             not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N) */

/*     X       (input) DOUBLE PRECISION array, dimension */
/*             ( 1 + ( N - 1 )*abs( INCX ) ). */
/*             On entry, elements 1, INCX+1, .., ( N - 1 )*INCX + 1 of */
/*             this array must contain the elements of the vector X. */

/*     INCX    (input) INTEGER */
/*             The increment for the elements of X. IF INCX < 0 then the */
/*             elements of X are accessed in reversed order.  INCX <> 0. */

/*     BETA    (input) DOUBLE PRECISION */
/*             The scalar beta. If beta is zero then Y need not be set on */
/*             input. */

/*     Y       (input/output) DOUBLE PRECISION array, dimension */
/*             ( 1 + ( N - 1 )*abs( INCY ) ). */
/*             On entry, elements 1, INCY+1, .., ( N - 1 )*INCY + 1 of */
/*             this array must contain the elements of the vector Y. */
/*             On exit, elements 1, INCY+1, .., ( N - 1 )*INCY + 1 of */
/*             this array contain the updated elements of the vector Y. */

/*     INCY    (input) INTEGER */
/*             The increment for the elements of Y. IF INCY < 0 then the */
/*             elements of Y are accessed in reversed order.  INCY <> 0. */

/*     NUMERICAL ASPECTS */

/*     Though being almost identical with the vanilla implementation */
/*     of the BLAS routine DSYMV the performance of this routine could */
/*     be significantly lower in the case of vendor supplied, highly */
/*     optimized BLAS. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, May 2008 (SLICOT version of the HAPACK routine DSKMV). */

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

#line 139 "MB01MD.f"
    /* Parameter adjustments */
#line 139 "MB01MD.f"
    a_dim1 = *lda;
#line 139 "MB01MD.f"
    a_offset = 1 + a_dim1;
#line 139 "MB01MD.f"
    a -= a_offset;
#line 139 "MB01MD.f"
    --x;
#line 139 "MB01MD.f"
    --y;
#line 139 "MB01MD.f"

#line 139 "MB01MD.f"
    /* Function Body */
#line 139 "MB01MD.f"
    info = 0;
#line 140 "MB01MD.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 141 "MB01MD.f"
	info = 1;
#line 142 "MB01MD.f"
    } else if (*n < 0) {
#line 143 "MB01MD.f"
	info = 2;
#line 144 "MB01MD.f"
    } else if (*lda < max(1,*n)) {
#line 145 "MB01MD.f"
	info = 5;
#line 146 "MB01MD.f"
    } else if (*incx == 0) {
#line 147 "MB01MD.f"
	info = 7;
#line 148 "MB01MD.f"
    } else if (*incy == 0) {
#line 149 "MB01MD.f"
	info = 10;
#line 150 "MB01MD.f"
    }
#line 151 "MB01MD.f"
    if (info != 0) {
#line 152 "MB01MD.f"
	xerbla_("MB01MD", &info, (ftnlen)6);
#line 153 "MB01MD.f"
	return 0;
#line 154 "MB01MD.f"
    }

/*     Quick return if possible. */

#line 158 "MB01MD.f"
    if (*n == 0 || *alpha == 0. && *beta == 1.) {
#line 158 "MB01MD.f"
	return 0;
#line 158 "MB01MD.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 163 "MB01MD.f"
    if (*incx > 0) {
#line 164 "MB01MD.f"
	kx = 1;
#line 165 "MB01MD.f"
    } else {
#line 166 "MB01MD.f"
	kx = 1 - (*n - 1) * *incx;
#line 167 "MB01MD.f"
    }
#line 168 "MB01MD.f"
    if (*incy > 0) {
#line 169 "MB01MD.f"
	ky = 1;
#line 170 "MB01MD.f"
    } else {
#line 171 "MB01MD.f"
	ky = 1 - (*n - 1) * *incy;
#line 172 "MB01MD.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

/*     First form  y := beta*y. */

#line 180 "MB01MD.f"
    if (*beta != 1.) {
#line 181 "MB01MD.f"
	if (*incy == 1) {
#line 182 "MB01MD.f"
	    if (*beta == 0.) {
#line 183 "MB01MD.f"
		i__1 = *n;
#line 183 "MB01MD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 184 "MB01MD.f"
		    y[i__] = 0.;
#line 185 "MB01MD.f"
/* L10: */
#line 185 "MB01MD.f"
		}
#line 186 "MB01MD.f"
	    } else {
#line 187 "MB01MD.f"
		i__1 = *n;
#line 187 "MB01MD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 188 "MB01MD.f"
		    y[i__] = *beta * y[i__];
#line 189 "MB01MD.f"
/* L20: */
#line 189 "MB01MD.f"
		}
#line 190 "MB01MD.f"
	    }
#line 191 "MB01MD.f"
	} else {
#line 192 "MB01MD.f"
	    iy = ky;
#line 193 "MB01MD.f"
	    if (*beta == 0.) {
#line 194 "MB01MD.f"
		i__1 = *n;
#line 194 "MB01MD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 195 "MB01MD.f"
		    y[iy] = 0.;
#line 196 "MB01MD.f"
		    iy += *incy;
#line 197 "MB01MD.f"
/* L30: */
#line 197 "MB01MD.f"
		}
#line 198 "MB01MD.f"
	    } else {
#line 199 "MB01MD.f"
		i__1 = *n;
#line 199 "MB01MD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 200 "MB01MD.f"
		    y[iy] = *beta * y[iy];
#line 201 "MB01MD.f"
		    iy += *incy;
#line 202 "MB01MD.f"
/* L40: */
#line 202 "MB01MD.f"
		}
#line 203 "MB01MD.f"
	    }
#line 204 "MB01MD.f"
	}
#line 205 "MB01MD.f"
    }

/*     Quick return if possible. */

#line 209 "MB01MD.f"
    if (*alpha == 0.) {
#line 209 "MB01MD.f"
	return 0;
#line 209 "MB01MD.f"
    }
#line 211 "MB01MD.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form y when A is stored in upper triangle. */

#line 215 "MB01MD.f"
	if (*incx == 1 && *incy == 1) {
#line 216 "MB01MD.f"
	    i__1 = *n;
#line 216 "MB01MD.f"
	    for (j = 2; j <= i__1; ++j) {
#line 217 "MB01MD.f"
		temp1 = *alpha * x[j];
#line 218 "MB01MD.f"
		temp2 = 0.;
#line 219 "MB01MD.f"
		i__2 = j - 1;
#line 219 "MB01MD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 220 "MB01MD.f"
		    y[i__] += temp1 * a[i__ + j * a_dim1];
#line 221 "MB01MD.f"
		    temp2 += a[i__ + j * a_dim1] * x[i__];
#line 222 "MB01MD.f"
/* L50: */
#line 222 "MB01MD.f"
		}
#line 223 "MB01MD.f"
		y[j] -= *alpha * temp2;
#line 224 "MB01MD.f"
/* L60: */
#line 224 "MB01MD.f"
	    }
#line 225 "MB01MD.f"
	} else {
#line 226 "MB01MD.f"
	    jx = kx + *incx;
#line 227 "MB01MD.f"
	    jy = ky + *incy;
#line 228 "MB01MD.f"
	    i__1 = *n;
#line 228 "MB01MD.f"
	    for (j = 2; j <= i__1; ++j) {
#line 229 "MB01MD.f"
		temp1 = *alpha * x[jx];
#line 230 "MB01MD.f"
		temp2 = 0.;
#line 231 "MB01MD.f"
		ix = kx;
#line 232 "MB01MD.f"
		iy = ky;
#line 233 "MB01MD.f"
		i__2 = j - 1;
#line 233 "MB01MD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 234 "MB01MD.f"
		    y[iy] += temp1 * a[i__ + j * a_dim1];
#line 235 "MB01MD.f"
		    temp2 += a[i__ + j * a_dim1] * x[ix];
#line 236 "MB01MD.f"
		    ix += *incx;
#line 237 "MB01MD.f"
		    iy += *incy;
#line 238 "MB01MD.f"
/* L70: */
#line 238 "MB01MD.f"
		}
#line 239 "MB01MD.f"
		y[jy] -= *alpha * temp2;
#line 240 "MB01MD.f"
		jx += *incx;
#line 241 "MB01MD.f"
		jy += *incy;
#line 242 "MB01MD.f"
/* L80: */
#line 242 "MB01MD.f"
	    }
#line 243 "MB01MD.f"
	}
#line 244 "MB01MD.f"
    } else {

/*        Form y when A is stored in lower triangle. */

#line 248 "MB01MD.f"
	if (*incx == 1 && *incy == 1) {
#line 249 "MB01MD.f"
	    i__1 = *n - 1;
#line 249 "MB01MD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 250 "MB01MD.f"
		temp1 = *alpha * x[j];
#line 251 "MB01MD.f"
		temp2 = 0.;
#line 252 "MB01MD.f"
		i__2 = *n;
#line 252 "MB01MD.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 253 "MB01MD.f"
		    y[i__] += temp1 * a[i__ + j * a_dim1];
#line 254 "MB01MD.f"
		    temp2 += a[i__ + j * a_dim1] * x[i__];
#line 255 "MB01MD.f"
/* L90: */
#line 255 "MB01MD.f"
		}
#line 256 "MB01MD.f"
		y[j] -= *alpha * temp2;
#line 257 "MB01MD.f"
/* L100: */
#line 257 "MB01MD.f"
	    }
#line 258 "MB01MD.f"
	} else {
#line 259 "MB01MD.f"
	    jx = kx;
#line 260 "MB01MD.f"
	    jy = ky;
#line 261 "MB01MD.f"
	    i__1 = *n - 1;
#line 261 "MB01MD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 262 "MB01MD.f"
		temp1 = *alpha * x[jx];
#line 263 "MB01MD.f"
		temp2 = 0.;
#line 264 "MB01MD.f"
		ix = jx;
#line 265 "MB01MD.f"
		iy = jy;
#line 266 "MB01MD.f"
		i__2 = *n;
#line 266 "MB01MD.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 267 "MB01MD.f"
		    ix += *incx;
#line 268 "MB01MD.f"
		    iy += *incy;
#line 269 "MB01MD.f"
		    y[iy] += temp1 * a[i__ + j * a_dim1];
#line 270 "MB01MD.f"
		    temp2 += a[i__ + j * a_dim1] * x[ix];
#line 271 "MB01MD.f"
/* L110: */
#line 271 "MB01MD.f"
		}
#line 272 "MB01MD.f"
		y[jy] -= *alpha * temp2;
#line 273 "MB01MD.f"
		jx += *incx;
#line 274 "MB01MD.f"
		jy += *incy;
#line 275 "MB01MD.f"
/* L120: */
#line 275 "MB01MD.f"
	    }
#line 276 "MB01MD.f"
	}
#line 277 "MB01MD.f"
    }
/* *** Last line of MB01MD *** */
#line 279 "MB01MD.f"
    return 0;
} /* mb01md_ */

