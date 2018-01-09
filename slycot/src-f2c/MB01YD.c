#line 1 "MB01YD.f"
/* MB01YD.f -- translated by f2c (version 20100827).
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

#line 1 "MB01YD.f"
/* Table of constant values */

static doublereal c_b8 = 0.;
static integer c__0 = 0;
static doublereal c_b12 = 1.;
static integer c__1 = 1;

/* Subroutine */ int mb01yd_(char *uplo, char *trans, integer *n, integer *k, 
	integer *l, doublereal *alpha, doublereal *beta, doublereal *a, 
	integer *lda, doublereal *c__, integer *ldc, integer *info, ftnlen 
	uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, m;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer ncola;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dlaset_(char *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical transp;


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

/*     To perform the symmetric rank k operations */

/*        C := alpha*op( A )*op( A )' + beta*C, */

/*     where alpha and beta are scalars, C is an n-by-n symmetric matrix, */
/*     op( A ) is an n-by-k matrix, and op( A ) is one of */

/*        op( A ) = A   or   op( A ) = A'. */

/*     The matrix A has l nonzero codiagonals, either upper or lower. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPLO    CHARACTER*1 */
/*             Specifies which triangle of the symmetric matrix C */
/*             is given and computed, as follows: */
/*             = 'U':  the upper triangular part is given/computed; */
/*             = 'L':  the lower triangular part is given/computed. */
/*             UPLO also defines the pattern of the matrix A (see below). */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of op( A ) to be used, as follows: */
/*             = 'N':  op( A ) = A; */
/*             = 'T':  op( A ) = A'; */
/*             = 'C':  op( A ) = A'. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix C.  N >= 0. */

/*     K       (input) INTEGER */
/*             The number of columns of the matrix op( A ).  K >= 0. */

/*     L       (input) INTEGER */
/*             If UPLO = 'U', matrix A has L nonzero subdiagonals. */
/*             If UPLO = 'L', matrix A has L nonzero superdiagonals. */
/*             MAX(0,NR-1) >= L >= 0, if UPLO = 'U', */
/*             MAX(0,NC-1) >= L >= 0, if UPLO = 'L', */
/*             where NR and NC are the numbers of rows and columns of the */
/*             matrix A, respectively. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar alpha. When alpha is zero then the array A is */
/*             not referenced. */

/*     BETA    (input) DOUBLE PRECISION */
/*             The scalar beta. When beta is zero then the array C need */
/*             not be set before entry. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,NC), where */
/*             NC is K when TRANS = 'N', and is N otherwise. */
/*             If TRANS = 'N', the leading N-by-K part of this array must */
/*             contain the matrix A, otherwise the leading K-by-N part of */
/*             this array must contain the matrix A. */
/*             If UPLO = 'U', only the upper triangular part and the */
/*             first L subdiagonals are referenced, and the remaining */
/*             subdiagonals are assumed to be zero. */
/*             If UPLO = 'L', only the lower triangular part and the */
/*             first L superdiagonals are referenced, and the remaining */
/*             superdiagonals are assumed to be zero. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= max(1,NR), */
/*             where NR = N, if TRANS = 'N', and NR = K, otherwise. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry with UPLO = 'U', the leading N-by-N upper */
/*             triangular part of this array must contain the upper */
/*             triangular part of the symmetric matrix C. */
/*             On entry with UPLO = 'L', the leading N-by-N lower */
/*             triangular part of this array must contain the lower */
/*             triangular part of the symmetric matrix C. */
/*             On exit, the leading N-by-N upper triangular part (if */
/*             UPLO = 'U'), or lower triangular part (if UPLO = 'L'), of */
/*             this array contains the corresponding triangular part of */
/*             the updated matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The calculations are efficiently performed taking the symmetry */
/*     and structure into account. */

/*     FURTHER COMMENTS */

/*     The matrix A may have the following patterns, when n = 7, k = 5, */
/*     and l = 2 are used for illustration: */

/*     UPLO = 'U', TRANS = 'N'         UPLO = 'L', TRANS = 'N' */

/*            [ x x x x x ]                   [ x x x 0 0 ] */
/*            [ x x x x x ]                   [ x x x x 0 ] */
/*            [ x x x x x ]                   [ x x x x x ] */
/*        A = [ 0 x x x x ],              A = [ x x x x x ], */
/*            [ 0 0 x x x ]                   [ x x x x x ] */
/*            [ 0 0 0 x x ]                   [ x x x x x ] */
/*            [ 0 0 0 0 x ]                   [ x x x x x ] */

/*     UPLO = 'U', TRANS = 'T'         UPLO = 'L', TRANS = 'T' */

/*            [ x x x x x x x ]               [ x x x 0 0 0 0 ] */
/*            [ x x x x x x x ]               [ x x x x 0 0 0 ] */
/*        A = [ x x x x x x x ],          A = [ x x x x x 0 0 ]. */
/*            [ 0 x x x x x x ]               [ x x x x x x 0 ] */
/*            [ 0 0 x x x x x ]               [ x x x x x x x ] */

/*     If N = K, the matrix A is upper or lower triangular, for L = 0, */
/*     and upper or lower Hessenberg, for L = 1. */

/*     This routine is a specialization of the BLAS 3 routine DSYRK. */
/*     BLAS 1 calls are used when appropriate, instead of in-line code, */
/*     in order to increase the efficiency. If the matrix A is full, or */
/*     its zero triangle has small order, an optimized DSYRK code could */
/*     be faster than MB01YD. */

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

#line 197 "MB01YD.f"
    /* Parameter adjustments */
#line 197 "MB01YD.f"
    a_dim1 = *lda;
#line 197 "MB01YD.f"
    a_offset = 1 + a_dim1;
#line 197 "MB01YD.f"
    a -= a_offset;
#line 197 "MB01YD.f"
    c_dim1 = *ldc;
#line 197 "MB01YD.f"
    c_offset = 1 + c_dim1;
#line 197 "MB01YD.f"
    c__ -= c_offset;
#line 197 "MB01YD.f"

#line 197 "MB01YD.f"
    /* Function Body */
#line 197 "MB01YD.f"
    *info = 0;
#line 198 "MB01YD.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 199 "MB01YD.f"
    transp = lsame_(trans, "T", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1);

#line 201 "MB01YD.f"
    if (transp) {
#line 202 "MB01YD.f"
	nrowa = *k;
#line 203 "MB01YD.f"
	ncola = *n;
#line 204 "MB01YD.f"
    } else {
#line 205 "MB01YD.f"
	nrowa = *n;
#line 206 "MB01YD.f"
	ncola = *k;
#line 207 "MB01YD.f"
    }

#line 209 "MB01YD.f"
    if (upper) {
#line 210 "MB01YD.f"
	m = nrowa;
#line 211 "MB01YD.f"
    } else {
#line 212 "MB01YD.f"
	m = ncola;
#line 213 "MB01YD.f"
    }

#line 215 "MB01YD.f"
    if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 216 "MB01YD.f"
	*info = -1;
#line 217 "MB01YD.f"
    } else if (! (transp || lsame_(trans, "N", (ftnlen)1, (ftnlen)1))) {
#line 218 "MB01YD.f"
	*info = -2;
#line 219 "MB01YD.f"
    } else if (*n < 0) {
#line 220 "MB01YD.f"
	*info = -3;
#line 221 "MB01YD.f"
    } else if (*k < 0) {
#line 222 "MB01YD.f"
	*info = -4;
#line 223 "MB01YD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 223 "MB01YD.f"
	i__1 = 0, i__2 = m - 1;
#line 223 "MB01YD.f"
	if (*l < 0 || *l > max(i__1,i__2)) {
#line 224 "MB01YD.f"
	    *info = -5;
#line 225 "MB01YD.f"
	} else if (*lda < max(1,nrowa)) {
#line 226 "MB01YD.f"
	    *info = -9;
#line 227 "MB01YD.f"
	} else if (*ldc < max(1,*n)) {
#line 228 "MB01YD.f"
	    *info = -11;
#line 229 "MB01YD.f"
	}
#line 229 "MB01YD.f"
    }

#line 231 "MB01YD.f"
    if (*info != 0) {

/*        Error return. */

#line 235 "MB01YD.f"
	i__1 = -(*info);
#line 235 "MB01YD.f"
	xerbla_("MB01YD", &i__1, (ftnlen)6);
#line 236 "MB01YD.f"
	return 0;
#line 237 "MB01YD.f"
    }

/*     Quick return, if possible. */

#line 241 "MB01YD.f"
    if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
#line 241 "MB01YD.f"
	return 0;
#line 241 "MB01YD.f"
    }

#line 245 "MB01YD.f"
    if (*alpha == 0.) {
#line 246 "MB01YD.f"
	if (*beta == 0.) {

/*           Special case when both alpha = 0 and beta = 0. */

#line 250 "MB01YD.f"
	    dlaset_(uplo, n, n, &c_b8, &c_b8, &c__[c_offset], ldc, (ftnlen)1);
#line 251 "MB01YD.f"
	} else {

/*           Special case alpha = 0. */

#line 255 "MB01YD.f"
	    dlascl_(uplo, &c__0, &c__0, &c_b12, beta, n, n, &c__[c_offset], 
		    ldc, info, (ftnlen)1);
#line 256 "MB01YD.f"
	}
#line 257 "MB01YD.f"
	return 0;
#line 258 "MB01YD.f"
    }

/*     General case: alpha <> 0. */

#line 262 "MB01YD.f"
    if (! transp) {

/*        Form  C := alpha*A*A' + beta*C. */

#line 266 "MB01YD.f"
	if (upper) {

#line 268 "MB01YD.f"
	    i__1 = *n;
#line 268 "MB01YD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 269 "MB01YD.f"
		if (*beta == 0.) {

#line 271 "MB01YD.f"
		    i__2 = j;
#line 271 "MB01YD.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 272 "MB01YD.f"
			c__[i__ + j * c_dim1] = 0.;
#line 273 "MB01YD.f"
/* L10: */
#line 273 "MB01YD.f"
		    }

#line 275 "MB01YD.f"
		} else if (*beta != 1.) {
#line 276 "MB01YD.f"
		    dscal_(&j, beta, &c__[j * c_dim1 + 1], &c__1);
#line 277 "MB01YD.f"
		}

/* Computing MAX */
#line 279 "MB01YD.f"
		i__2 = 1, i__3 = j - *l;
#line 279 "MB01YD.f"
		i__4 = *k;
#line 279 "MB01YD.f"
		for (m = max(i__2,i__3); m <= i__4; ++m) {
/* Computing MIN */
#line 280 "MB01YD.f"
		    i__3 = j, i__5 = *l + m;
#line 280 "MB01YD.f"
		    i__2 = min(i__3,i__5);
#line 280 "MB01YD.f"
		    d__1 = *alpha * a[j + m * a_dim1];
#line 280 "MB01YD.f"
		    daxpy_(&i__2, &d__1, &a[m * a_dim1 + 1], &c__1, &c__[j * 
			    c_dim1 + 1], &c__1);
#line 282 "MB01YD.f"
/* L20: */
#line 282 "MB01YD.f"
		}

#line 284 "MB01YD.f"
/* L30: */
#line 284 "MB01YD.f"
	    }

#line 286 "MB01YD.f"
	} else {

#line 288 "MB01YD.f"
	    i__1 = *n;
#line 288 "MB01YD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 289 "MB01YD.f"
		if (*beta == 0.) {

#line 291 "MB01YD.f"
		    i__4 = *n;
#line 291 "MB01YD.f"
		    for (i__ = j; i__ <= i__4; ++i__) {
#line 292 "MB01YD.f"
			c__[i__ + j * c_dim1] = 0.;
#line 293 "MB01YD.f"
/* L40: */
#line 293 "MB01YD.f"
		    }

#line 295 "MB01YD.f"
		} else if (*beta != 1.) {
#line 296 "MB01YD.f"
		    i__4 = *n - j + 1;
#line 296 "MB01YD.f"
		    dscal_(&i__4, beta, &c__[j + j * c_dim1], &c__1);
#line 297 "MB01YD.f"
		}

/* Computing MIN */
#line 299 "MB01YD.f"
		i__2 = j + *l;
#line 299 "MB01YD.f"
		i__4 = min(i__2,*k);
#line 299 "MB01YD.f"
		for (m = 1; m <= i__4; ++m) {
#line 300 "MB01YD.f"
		    i__2 = *n - j + 1;
#line 300 "MB01YD.f"
		    d__1 = *alpha * a[j + m * a_dim1];
#line 300 "MB01YD.f"
		    daxpy_(&i__2, &d__1, &a[j + m * a_dim1], &c__1, &c__[j + 
			    j * c_dim1], &c__1);
#line 302 "MB01YD.f"
/* L50: */
#line 302 "MB01YD.f"
		}

#line 304 "MB01YD.f"
/* L60: */
#line 304 "MB01YD.f"
	    }

#line 306 "MB01YD.f"
	}

#line 308 "MB01YD.f"
    } else {

/*        Form  C := alpha*A'*A + beta*C. */

#line 312 "MB01YD.f"
	if (upper) {

#line 314 "MB01YD.f"
	    i__1 = *n;
#line 314 "MB01YD.f"
	    for (j = 1; j <= i__1; ++j) {

#line 316 "MB01YD.f"
		i__4 = j;
#line 316 "MB01YD.f"
		for (i__ = 1; i__ <= i__4; ++i__) {
/* Computing MIN */
#line 317 "MB01YD.f"
		    i__3 = j + *l;
#line 317 "MB01YD.f"
		    i__2 = min(i__3,*k);
#line 317 "MB01YD.f"
		    temp = *alpha * ddot_(&i__2, &a[i__ * a_dim1 + 1], &c__1, 
			    &a[j * a_dim1 + 1], &c__1);
#line 319 "MB01YD.f"
		    if (*beta == 0.) {
#line 320 "MB01YD.f"
			c__[i__ + j * c_dim1] = temp;
#line 321 "MB01YD.f"
		    } else {
#line 322 "MB01YD.f"
			c__[i__ + j * c_dim1] = temp + *beta * c__[i__ + j * 
				c_dim1];
#line 323 "MB01YD.f"
		    }
#line 324 "MB01YD.f"
/* L70: */
#line 324 "MB01YD.f"
		}

#line 326 "MB01YD.f"
/* L80: */
#line 326 "MB01YD.f"
	    }

#line 328 "MB01YD.f"
	} else {

#line 330 "MB01YD.f"
	    i__1 = *n;
#line 330 "MB01YD.f"
	    for (j = 1; j <= i__1; ++j) {

#line 332 "MB01YD.f"
		i__4 = *n;
#line 332 "MB01YD.f"
		for (i__ = j; i__ <= i__4; ++i__) {
/* Computing MAX */
#line 333 "MB01YD.f"
		    i__2 = 1, i__3 = i__ - *l;
#line 333 "MB01YD.f"
		    m = max(i__2,i__3);
#line 334 "MB01YD.f"
		    i__2 = *k - m + 1;
#line 334 "MB01YD.f"
		    temp = *alpha * ddot_(&i__2, &a[m + i__ * a_dim1], &c__1, 
			    &a[m + j * a_dim1], &c__1);
#line 336 "MB01YD.f"
		    if (*beta == 0.) {
#line 337 "MB01YD.f"
			c__[i__ + j * c_dim1] = temp;
#line 338 "MB01YD.f"
		    } else {
#line 339 "MB01YD.f"
			c__[i__ + j * c_dim1] = temp + *beta * c__[i__ + j * 
				c_dim1];
#line 340 "MB01YD.f"
		    }
#line 341 "MB01YD.f"
/* L90: */
#line 341 "MB01YD.f"
		}

#line 343 "MB01YD.f"
/* L100: */
#line 343 "MB01YD.f"
	    }

#line 345 "MB01YD.f"
	}

#line 347 "MB01YD.f"
    }

#line 349 "MB01YD.f"
    return 0;

/* *** Last line of MB01YD *** */
} /* mb01yd_ */

