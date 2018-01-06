#line 1 "MB01ZD.f"
/* MB01ZD.f -- translated by f2c (version 20100827).
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

#line 1 "MB01ZD.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mb01zd_(char *side, char *uplo, char *transt, char *diag,
	 integer *m, integer *n, integer *l, doublereal *alpha, doublereal *t,
	 integer *ldt, doublereal *h__, integer *ldh, integer *info, ftnlen 
	side_len, ftnlen uplo_len, ftnlen transt_len, ftnlen diag_len)
{
    /* System generated locals */
    integer h_dim1, h_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, i1, i2, m2;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical trans, upper;
    static integer nrowt;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical nounit;


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

/*     To compute the matrix product */

/*        H := alpha*op( T )*H,   or   H := alpha*H*op( T ), */

/*     where alpha is a scalar, H is an m-by-n upper or lower */
/*     Hessenberg-like matrix (with l nonzero subdiagonals or */
/*     superdiagonals, respectively), T is a unit, or non-unit, */
/*     upper or lower triangular matrix, and op( T ) is one of */

/*        op( T ) = T   or   op( T ) = T'. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SIDE    CHARACTER*1 */
/*             Specifies whether the triangular matrix T appears on the */
/*             left or right in the matrix product, as follows: */
/*             = 'L':  the product alpha*op( T )*H is computed; */
/*             = 'R':  the product alpha*H*op( T ) is computed. */

/*     UPLO    CHARACTER*1 */
/*             Specifies the form of the matrices T and H, as follows: */
/*             = 'U':  the matrix T is upper triangular and the matrix H */
/*                     is upper Hessenberg-like; */
/*             = 'L':  the matrix T is lower triangular and the matrix H */
/*                     is lower Hessenberg-like. */

/*     TRANST  CHARACTER*1 */
/*             Specifies the form of op( T ) to be used, as follows: */
/*             = 'N':  op( T ) = T; */
/*             = 'T':  op( T ) = T'; */
/*             = 'C':  op( T ) = T'. */

/*     DIAG    CHARACTER*1. */
/*             Specifies whether or not T is unit triangular, as follows: */
/*             = 'U':  the matrix T is assumed to be unit triangular; */
/*             = 'N':  the matrix T is not assumed to be unit triangular. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of H.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of H.  N >= 0. */

/*     L       (input) INTEGER */
/*             If UPLO = 'U', matrix H has L nonzero subdiagonals. */
/*             If UPLO = 'L', matrix H has L nonzero superdiagonals. */
/*             MAX(0,M-1) >= L >= 0, if UPLO = 'U'; */
/*             MAX(0,N-1) >= L >= 0, if UPLO = 'L'. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar alpha. When alpha is zero then T is not */
/*             referenced and H need not be set before entry. */

/*     T       (input) DOUBLE PRECISION array, dimension (LDT,k), where */
/*             k is m when SIDE = 'L' and is n when SIDE = 'R'. */
/*             If UPLO = 'U', the leading k-by-k upper triangular part */
/*             of this array must contain the upper triangular matrix T */
/*             and the strictly lower triangular part is not referenced. */
/*             If UPLO = 'L', the leading k-by-k lower triangular part */
/*             of this array must contain the lower triangular matrix T */
/*             and the strictly upper triangular part is not referenced. */
/*             Note that when DIAG = 'U', the diagonal elements of T are */
/*             not referenced either, but are assumed to be unity. */

/*     LDT     INTEGER */
/*             The leading dimension of array T. */
/*             LDT >= MAX(1,M), if SIDE = 'L'; */
/*             LDT >= MAX(1,N), if SIDE = 'R'. */

/*     H       (input/output) DOUBLE PRECISION array, dimension (LDH,N) */
/*             On entry, if UPLO = 'U', the leading M-by-N upper */
/*             Hessenberg part of this array must contain the upper */
/*             Hessenberg-like matrix H. */
/*             On entry, if UPLO = 'L', the leading M-by-N lower */
/*             Hessenberg part of this array must contain the lower */
/*             Hessenberg-like matrix H. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the matrix product alpha*op( T )*H, if SIDE = 'L', */
/*             or alpha*H*op( T ), if SIDE = 'R'. If TRANST = 'N', this */
/*             product has the same pattern as the given matrix H; */
/*             the elements below the L-th subdiagonal (if UPLO = 'U'), */
/*             or above the L-th superdiagonal (if UPLO = 'L'), are not */
/*             referenced in this case. If TRANST = 'T', the elements */
/*             below the (N+L)-th row (if UPLO = 'U', SIDE = 'R', and */
/*             M > N+L), or at the right of the (M+L)-th column */
/*             (if UPLO = 'L', SIDE = 'L', and N > M+L), are not set to */
/*             zero nor referenced. */

/*     LDH     INTEGER */
/*             The leading dimension of array H.  LDH >= max(1,M). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The calculations are efficiently performed taking the problem */
/*     structure into account. */

/*     FURTHER COMMENTS */

/*     The matrix H may have the following patterns, when m = 7, n = 6, */
/*     and l = 2 are used for illustration: */

/*               UPLO = 'U'                    UPLO = 'L' */

/*            [ x x x x x x ]               [ x x x 0 0 0 ] */
/*            [ x x x x x x ]               [ x x x x 0 0 ] */
/*            [ x x x x x x ]               [ x x x x x 0 ] */
/*        H = [ 0 x x x x x ],          H = [ x x x x x x ]. */
/*            [ 0 0 x x x x ]               [ x x x x x x ] */
/*            [ 0 0 0 x x x ]               [ x x x x x x ] */
/*            [ 0 0 0 0 x x ]               [ x x x x x x ] */

/*     The products T*H or H*T have the same pattern as H, but the */
/*     products T'*H or H*T' may be full matrices. */

/*     If m = n, the matrix H is upper or lower triangular, for l = 0, */
/*     and upper or lower Hessenberg, for l = 1. */

/*     This routine is a specialization of the BLAS 3 routine DTRMM. */
/*     BLAS 1 calls are used when appropriate, instead of in-line code, */
/*     in order to increase the efficiency. If the matrix H is full, or */
/*     its zero triangle has small order, an optimized DTRMM code could */
/*     be faster than MB01ZD. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Nov. 2000. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

#line 204 "MB01ZD.f"
    /* Parameter adjustments */
#line 204 "MB01ZD.f"
    t_dim1 = *ldt;
#line 204 "MB01ZD.f"
    t_offset = 1 + t_dim1;
#line 204 "MB01ZD.f"
    t -= t_offset;
#line 204 "MB01ZD.f"
    h_dim1 = *ldh;
#line 204 "MB01ZD.f"
    h_offset = 1 + h_dim1;
#line 204 "MB01ZD.f"
    h__ -= h_offset;
#line 204 "MB01ZD.f"

#line 204 "MB01ZD.f"
    /* Function Body */
#line 204 "MB01ZD.f"
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 205 "MB01ZD.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 206 "MB01ZD.f"
    trans = lsame_(transt, "T", (ftnlen)1, (ftnlen)1) || lsame_(transt, "C", (
	    ftnlen)1, (ftnlen)1);
#line 207 "MB01ZD.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 208 "MB01ZD.f"
    if (lside) {
#line 209 "MB01ZD.f"
	nrowt = *m;
#line 210 "MB01ZD.f"
    } else {
#line 211 "MB01ZD.f"
	nrowt = *n;
#line 212 "MB01ZD.f"
    }

#line 214 "MB01ZD.f"
    if (upper) {
#line 215 "MB01ZD.f"
	m2 = *m;
#line 216 "MB01ZD.f"
    } else {
#line 217 "MB01ZD.f"
	m2 = *n;
#line 218 "MB01ZD.f"
    }

#line 220 "MB01ZD.f"
    *info = 0;
#line 221 "MB01ZD.f"
    if (! (lside || lsame_(side, "R", (ftnlen)1, (ftnlen)1))) {
#line 222 "MB01ZD.f"
	*info = -1;
#line 223 "MB01ZD.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 224 "MB01ZD.f"
	*info = -2;
#line 225 "MB01ZD.f"
    } else if (! (trans || lsame_(transt, "N", (ftnlen)1, (ftnlen)1))) {
#line 226 "MB01ZD.f"
	*info = -3;
#line 227 "MB01ZD.f"
    } else if (! (nounit || lsame_(diag, "U", (ftnlen)1, (ftnlen)1))) {
#line 228 "MB01ZD.f"
	*info = -4;
#line 229 "MB01ZD.f"
    } else if (*m < 0) {
#line 230 "MB01ZD.f"
	*info = -5;
#line 231 "MB01ZD.f"
    } else if (*n < 0) {
#line 232 "MB01ZD.f"
	*info = -6;
#line 233 "MB01ZD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 233 "MB01ZD.f"
	i__1 = 0, i__2 = m2 - 1;
#line 233 "MB01ZD.f"
	if (*l < 0 || *l > max(i__1,i__2)) {
#line 234 "MB01ZD.f"
	    *info = -7;
#line 235 "MB01ZD.f"
	} else if (*ldt < max(1,nrowt)) {
#line 236 "MB01ZD.f"
	    *info = -10;
#line 237 "MB01ZD.f"
	} else if (*ldh < max(1,*m)) {
#line 238 "MB01ZD.f"
	    *info = -12;
#line 239 "MB01ZD.f"
	}
#line 239 "MB01ZD.f"
    }

#line 241 "MB01ZD.f"
    if (*info != 0) {

/*        Error return. */

#line 245 "MB01ZD.f"
	i__1 = -(*info);
#line 245 "MB01ZD.f"
	xerbla_("MB01ZD", &i__1, (ftnlen)6);
#line 246 "MB01ZD.f"
	return 0;
#line 247 "MB01ZD.f"
    }

/*     Quick return, if possible. */

#line 251 "MB01ZD.f"
    if (min(*m,*n) == 0) {
#line 251 "MB01ZD.f"
	return 0;
#line 251 "MB01ZD.f"
    }

/*     Also, when alpha = 0. */

#line 256 "MB01ZD.f"
    if (*alpha == 0.) {

#line 258 "MB01ZD.f"
	i__1 = *n;
#line 258 "MB01ZD.f"
	for (j = 1; j <= i__1; ++j) {
#line 259 "MB01ZD.f"
	    if (upper) {
#line 260 "MB01ZD.f"
		i1 = 1;
/* Computing MIN */
#line 261 "MB01ZD.f"
		i__2 = j + *l;
#line 261 "MB01ZD.f"
		i2 = min(i__2,*m);
#line 262 "MB01ZD.f"
	    } else {
/* Computing MAX */
#line 263 "MB01ZD.f"
		i__2 = 1, i__3 = j - *l;
#line 263 "MB01ZD.f"
		i1 = max(i__2,i__3);
#line 264 "MB01ZD.f"
		i2 = *m;
#line 265 "MB01ZD.f"
	    }

#line 267 "MB01ZD.f"
	    i__2 = i2;
#line 267 "MB01ZD.f"
	    for (i__ = i1; i__ <= i__2; ++i__) {
#line 268 "MB01ZD.f"
		h__[i__ + j * h_dim1] = 0.;
#line 269 "MB01ZD.f"
/* L10: */
#line 269 "MB01ZD.f"
	    }

#line 271 "MB01ZD.f"
/* L20: */
#line 271 "MB01ZD.f"
	}

#line 273 "MB01ZD.f"
	return 0;
#line 274 "MB01ZD.f"
    }

/*     Start the operations. */

#line 278 "MB01ZD.f"
    if (lside) {
#line 279 "MB01ZD.f"
	if (! trans) {

/*           Form  H := alpha*T*H. */

#line 283 "MB01ZD.f"
	    if (upper) {

#line 285 "MB01ZD.f"
		i__1 = *n;
#line 285 "MB01ZD.f"
		for (j = 1; j <= i__1; ++j) {

/* Computing MIN */
#line 287 "MB01ZD.f"
		    i__3 = j + *l;
#line 287 "MB01ZD.f"
		    i__2 = min(i__3,*m);
#line 287 "MB01ZD.f"
		    for (k = 1; k <= i__2; ++k) {
#line 288 "MB01ZD.f"
			if (h__[k + j * h_dim1] != 0.) {
#line 289 "MB01ZD.f"
			    temp = *alpha * h__[k + j * h_dim1];
#line 290 "MB01ZD.f"
			    i__3 = k - 1;
#line 290 "MB01ZD.f"
			    daxpy_(&i__3, &temp, &t[k * t_dim1 + 1], &c__1, &
				    h__[j * h_dim1 + 1], &c__1);
#line 292 "MB01ZD.f"
			    if (nounit) {
#line 292 "MB01ZD.f"
				temp *= t[k + k * t_dim1];
#line 292 "MB01ZD.f"
			    }
#line 294 "MB01ZD.f"
			    h__[k + j * h_dim1] = temp;
#line 295 "MB01ZD.f"
			}
#line 296 "MB01ZD.f"
/* L30: */
#line 296 "MB01ZD.f"
		    }

#line 298 "MB01ZD.f"
/* L40: */
#line 298 "MB01ZD.f"
		}

#line 300 "MB01ZD.f"
	    } else {

#line 302 "MB01ZD.f"
		i__1 = *n;
#line 302 "MB01ZD.f"
		for (j = 1; j <= i__1; ++j) {

/* Computing MAX */
#line 304 "MB01ZD.f"
		    i__3 = 1, i__4 = j - *l;
#line 304 "MB01ZD.f"
		    i__2 = max(i__3,i__4);
#line 304 "MB01ZD.f"
		    for (k = *m; k >= i__2; --k) {
#line 305 "MB01ZD.f"
			if (h__[k + j * h_dim1] != 0.) {
#line 306 "MB01ZD.f"
			    temp = *alpha * h__[k + j * h_dim1];
#line 307 "MB01ZD.f"
			    h__[k + j * h_dim1] = temp;
#line 308 "MB01ZD.f"
			    if (nounit) {
#line 308 "MB01ZD.f"
				h__[k + j * h_dim1] *= t[k + k * t_dim1];
#line 308 "MB01ZD.f"
			    }
#line 310 "MB01ZD.f"
			    i__3 = *m - k;
#line 310 "MB01ZD.f"
			    daxpy_(&i__3, &temp, &t[k + 1 + k * t_dim1], &
				    c__1, &h__[k + 1 + j * h_dim1], &c__1);
#line 312 "MB01ZD.f"
			}
#line 313 "MB01ZD.f"
/* L50: */
#line 313 "MB01ZD.f"
		    }

#line 315 "MB01ZD.f"
/* L60: */
#line 315 "MB01ZD.f"
		}

#line 317 "MB01ZD.f"
	    }

#line 319 "MB01ZD.f"
	} else {

/*           Form  H := alpha*T'*H. */

#line 323 "MB01ZD.f"
	    if (upper) {

#line 325 "MB01ZD.f"
		i__1 = *n;
#line 325 "MB01ZD.f"
		for (j = 1; j <= i__1; ++j) {
#line 326 "MB01ZD.f"
		    i1 = j + *l;

#line 328 "MB01ZD.f"
		    for (i__ = *m; i__ >= 1; --i__) {
#line 329 "MB01ZD.f"
			if (i__ > i1) {
#line 330 "MB01ZD.f"
			    temp = ddot_(&i1, &t[i__ * t_dim1 + 1], &c__1, &
				    h__[j * h_dim1 + 1], &c__1);
#line 331 "MB01ZD.f"
			} else {
#line 332 "MB01ZD.f"
			    temp = h__[i__ + j * h_dim1];
#line 333 "MB01ZD.f"
			    if (nounit) {
#line 333 "MB01ZD.f"
				temp *= t[i__ + i__ * t_dim1];
#line 333 "MB01ZD.f"
			    }
#line 335 "MB01ZD.f"
			    i__2 = i__ - 1;
#line 335 "MB01ZD.f"
			    temp += ddot_(&i__2, &t[i__ * t_dim1 + 1], &c__1, 
				    &h__[j * h_dim1 + 1], &c__1);
#line 337 "MB01ZD.f"
			}
#line 338 "MB01ZD.f"
			h__[i__ + j * h_dim1] = *alpha * temp;
#line 339 "MB01ZD.f"
/* L70: */
#line 339 "MB01ZD.f"
		    }

#line 341 "MB01ZD.f"
/* L80: */
#line 341 "MB01ZD.f"
		}

#line 343 "MB01ZD.f"
	    } else {

/* Computing MIN */
#line 345 "MB01ZD.f"
		i__2 = *m + *l;
#line 345 "MB01ZD.f"
		i__1 = min(i__2,*n);
#line 345 "MB01ZD.f"
		for (j = 1; j <= i__1; ++j) {
#line 346 "MB01ZD.f"
		    i1 = j - *l;

#line 348 "MB01ZD.f"
		    i__2 = *m;
#line 348 "MB01ZD.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 349 "MB01ZD.f"
			if (i__ < i1) {
#line 350 "MB01ZD.f"
			    i__3 = *m - i1 + 1;
#line 350 "MB01ZD.f"
			    temp = ddot_(&i__3, &t[i1 + i__ * t_dim1], &c__1, 
				    &h__[i1 + j * h_dim1], &c__1);
#line 352 "MB01ZD.f"
			} else {
#line 353 "MB01ZD.f"
			    temp = h__[i__ + j * h_dim1];
#line 354 "MB01ZD.f"
			    if (nounit) {
#line 354 "MB01ZD.f"
				temp *= t[i__ + i__ * t_dim1];
#line 354 "MB01ZD.f"
			    }
#line 356 "MB01ZD.f"
			    i__3 = *m - i__;
#line 356 "MB01ZD.f"
			    temp += ddot_(&i__3, &t[i__ + 1 + i__ * t_dim1], &
				    c__1, &h__[i__ + 1 + j * h_dim1], &c__1);
#line 358 "MB01ZD.f"
			}
#line 359 "MB01ZD.f"
			h__[i__ + j * h_dim1] = *alpha * temp;
#line 360 "MB01ZD.f"
/* L90: */
#line 360 "MB01ZD.f"
		    }

#line 362 "MB01ZD.f"
/* L100: */
#line 362 "MB01ZD.f"
		}

#line 364 "MB01ZD.f"
	    }

#line 366 "MB01ZD.f"
	}

#line 368 "MB01ZD.f"
    } else {

#line 370 "MB01ZD.f"
	if (! trans) {

/*           Form  H := alpha*H*T. */

#line 374 "MB01ZD.f"
	    if (upper) {

#line 376 "MB01ZD.f"
		for (j = *n; j >= 1; --j) {
/* Computing MIN */
#line 377 "MB01ZD.f"
		    i__1 = j + *l;
#line 377 "MB01ZD.f"
		    i2 = min(i__1,*m);
#line 378 "MB01ZD.f"
		    temp = *alpha;
#line 379 "MB01ZD.f"
		    if (nounit) {
#line 379 "MB01ZD.f"
			temp *= t[j + j * t_dim1];
#line 379 "MB01ZD.f"
		    }
#line 381 "MB01ZD.f"
		    dscal_(&i2, &temp, &h__[j * h_dim1 + 1], &c__1);

#line 383 "MB01ZD.f"
		    i__1 = j - 1;
#line 383 "MB01ZD.f"
		    for (k = 1; k <= i__1; ++k) {
#line 384 "MB01ZD.f"
			d__1 = *alpha * t[k + j * t_dim1];
#line 384 "MB01ZD.f"
			daxpy_(&i2, &d__1, &h__[k * h_dim1 + 1], &c__1, &h__[
				j * h_dim1 + 1], &c__1);
#line 386 "MB01ZD.f"
/* L110: */
#line 386 "MB01ZD.f"
		    }

#line 388 "MB01ZD.f"
/* L120: */
#line 388 "MB01ZD.f"
		}

#line 390 "MB01ZD.f"
	    } else {

#line 392 "MB01ZD.f"
		i__1 = *n;
#line 392 "MB01ZD.f"
		for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 393 "MB01ZD.f"
		    i__2 = 1, i__3 = j - *l;
#line 393 "MB01ZD.f"
		    i1 = max(i__2,i__3);
#line 394 "MB01ZD.f"
		    temp = *alpha;
#line 395 "MB01ZD.f"
		    if (nounit) {
#line 395 "MB01ZD.f"
			temp *= t[j + j * t_dim1];
#line 395 "MB01ZD.f"
		    }
#line 397 "MB01ZD.f"
		    i__2 = *m - i1 + 1;
#line 397 "MB01ZD.f"
		    dscal_(&i__2, &temp, &h__[i1 + j * h_dim1], &c__1);

#line 399 "MB01ZD.f"
		    i__2 = *n;
#line 399 "MB01ZD.f"
		    for (k = j + 1; k <= i__2; ++k) {
#line 400 "MB01ZD.f"
			i__3 = *m - i1 + 1;
#line 400 "MB01ZD.f"
			d__1 = *alpha * t[k + j * t_dim1];
#line 400 "MB01ZD.f"
			daxpy_(&i__3, &d__1, &h__[i1 + k * h_dim1], &c__1, &
				h__[i1 + j * h_dim1], &c__1);
#line 402 "MB01ZD.f"
/* L130: */
#line 402 "MB01ZD.f"
		    }

#line 404 "MB01ZD.f"
/* L140: */
#line 404 "MB01ZD.f"
		}

#line 406 "MB01ZD.f"
	    }

#line 408 "MB01ZD.f"
	} else {

/*           Form  H := alpha*H*T'. */

#line 412 "MB01ZD.f"
	    if (upper) {
/* Computing MIN */
#line 413 "MB01ZD.f"
		i__1 = *n + *l;
#line 413 "MB01ZD.f"
		m2 = min(i__1,*m);

#line 415 "MB01ZD.f"
		i__1 = *n;
#line 415 "MB01ZD.f"
		for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
#line 416 "MB01ZD.f"
		    i__2 = k + *l;
#line 416 "MB01ZD.f"
		    i1 = min(i__2,*m);
/* Computing MIN */
#line 417 "MB01ZD.f"
		    i__2 = k + *l;
#line 417 "MB01ZD.f"
		    i2 = min(i__2,m2);

#line 419 "MB01ZD.f"
		    i__2 = k - 1;
#line 419 "MB01ZD.f"
		    for (j = 1; j <= i__2; ++j) {
#line 420 "MB01ZD.f"
			if (t[j + k * t_dim1] != 0.) {
#line 421 "MB01ZD.f"
			    temp = *alpha * t[j + k * t_dim1];
#line 422 "MB01ZD.f"
			    daxpy_(&i1, &temp, &h__[k * h_dim1 + 1], &c__1, &
				    h__[j * h_dim1 + 1], &c__1);

#line 425 "MB01ZD.f"
			    i__3 = i2;
#line 425 "MB01ZD.f"
			    for (i__ = i1 + 1; i__ <= i__3; ++i__) {
#line 426 "MB01ZD.f"
				h__[i__ + j * h_dim1] = temp * h__[i__ + k * 
					h_dim1];
#line 427 "MB01ZD.f"
/* L150: */
#line 427 "MB01ZD.f"
			    }

#line 429 "MB01ZD.f"
			}
#line 430 "MB01ZD.f"
/* L160: */
#line 430 "MB01ZD.f"
		    }

#line 432 "MB01ZD.f"
		    temp = *alpha;
#line 433 "MB01ZD.f"
		    if (nounit) {
#line 433 "MB01ZD.f"
			temp *= t[k + k * t_dim1];
#line 433 "MB01ZD.f"
		    }
#line 435 "MB01ZD.f"
		    if (temp != 1.) {
#line 435 "MB01ZD.f"
			dscal_(&i2, &temp, &h__[k * h_dim1 + 1], &c__1);
#line 435 "MB01ZD.f"
		    }
#line 437 "MB01ZD.f"
/* L170: */
#line 437 "MB01ZD.f"
		}

#line 439 "MB01ZD.f"
	    } else {

#line 441 "MB01ZD.f"
		for (k = *n; k >= 1; --k) {
/* Computing MAX */
#line 442 "MB01ZD.f"
		    i__1 = 1, i__2 = k - *l;
#line 442 "MB01ZD.f"
		    i1 = max(i__1,i__2);
/* Computing MAX */
#line 443 "MB01ZD.f"
		    i__1 = 1, i__2 = k - *l + 1;
#line 443 "MB01ZD.f"
		    i2 = max(i__1,i__2);
/* Computing MIN */
#line 444 "MB01ZD.f"
		    i__1 = *m, i__2 = i2 - 1;
#line 444 "MB01ZD.f"
		    m2 = min(i__1,i__2);

#line 446 "MB01ZD.f"
		    i__1 = *n;
#line 446 "MB01ZD.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 447 "MB01ZD.f"
			if (t[j + k * t_dim1] != 0.) {
#line 448 "MB01ZD.f"
			    temp = *alpha * t[j + k * t_dim1];
#line 449 "MB01ZD.f"
			    i__2 = *m - i2 + 1;
#line 449 "MB01ZD.f"
			    daxpy_(&i__2, &temp, &h__[i2 + k * h_dim1], &c__1,
				     &h__[i2 + j * h_dim1], &c__1);

#line 452 "MB01ZD.f"
			    i__2 = m2;
#line 452 "MB01ZD.f"
			    for (i__ = i1; i__ <= i__2; ++i__) {
#line 453 "MB01ZD.f"
				h__[i__ + j * h_dim1] = temp * h__[i__ + k * 
					h_dim1];
#line 454 "MB01ZD.f"
/* L180: */
#line 454 "MB01ZD.f"
			    }

#line 456 "MB01ZD.f"
			}
#line 457 "MB01ZD.f"
/* L190: */
#line 457 "MB01ZD.f"
		    }

#line 459 "MB01ZD.f"
		    temp = *alpha;
#line 460 "MB01ZD.f"
		    if (nounit) {
#line 460 "MB01ZD.f"
			temp *= t[k + k * t_dim1];
#line 460 "MB01ZD.f"
		    }
#line 462 "MB01ZD.f"
		    if (temp != 1.) {
#line 462 "MB01ZD.f"
			i__1 = *m - i1 + 1;
#line 462 "MB01ZD.f"
			dscal_(&i__1, &temp, &h__[i1 + k * h_dim1], &c__1);
#line 462 "MB01ZD.f"
		    }
#line 464 "MB01ZD.f"
/* L200: */
#line 464 "MB01ZD.f"
		}

#line 466 "MB01ZD.f"
	    }

#line 468 "MB01ZD.f"
	}

#line 470 "MB01ZD.f"
    }

#line 472 "MB01ZD.f"
    return 0;

/* *** Last line of MB01ZD *** */
} /* mb01zd_ */

