#line 1 "MB01RX.f"
/* MB01RX.f -- translated by f2c (version 20100827).
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

#line 1 "MB01RX.f"
/* Table of constant values */

static doublereal c_b10 = 0.;
static integer c__0 = 0;
static doublereal c_b14 = 1.;
static integer c__1 = 1;

/* Subroutine */ int mb01rx_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, doublereal *alpha, doublereal *beta, doublereal *r__, 
	integer *ldr, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, integer *info, ftnlen side_len, ftnlen uplo_len, ftnlen 
	trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, r_dim1, r_offset, i__1, i__2;

    /* Local variables */
    static integer j;
    static logical lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static logical luplo;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dlaset_(char *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical ltrans;


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

/*     To compute either the upper or lower triangular part of one of the */
/*     matrix formulas */
/*        _ */
/*        R = alpha*R + beta*op( A )*B,                               (1) */
/*        _ */
/*        R = alpha*R + beta*B*op( A ),                               (2) */
/*                                             _ */
/*     where alpha and beta are scalars, R and R are m-by-m matrices, */
/*     op( A ) and B are m-by-n and n-by-m matrices for (1), or n-by-m */
/*     and m-by-n matrices for (2), respectively, and op( A ) is one of */

/*        op( A ) = A   or   op( A ) = A',  the transpose of A. */

/*     The result is overwritten on R. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SIDE    CHARACTER*1 */
/*             Specifies whether the matrix A appears on the left or */
/*             right in the matrix product as follows: */
/*                     _ */
/*             = 'L':  R = alpha*R + beta*op( A )*B; */
/*                     _ */
/*             = 'R':  R = alpha*R + beta*B*op( A ). */

/*     UPLO    CHARACTER*1                               _ */
/*             Specifies which triangles of the matrices R and R are */
/*             computed and given, respectively, as follows: */
/*             = 'U':  the upper triangular part; */
/*             = 'L':  the lower triangular part. */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of op( A ) to be used in the matrix */
/*             multiplication as follows: */
/*             = 'N':  op( A ) = A; */
/*             = 'T':  op( A ) = A'; */
/*             = 'C':  op( A ) = A'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER           _ */
/*             The order of the matrices R and R, the number of rows of */
/*             the matrix op( A ) and the number of columns of the */
/*             matrix B, for SIDE = 'L', or the number of rows of the */
/*             matrix B and the number of columns of the matrix op( A ), */
/*             for SIDE = 'R'.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of rows of the matrix B and the number of */
/*             columns of the matrix op( A ), for SIDE = 'L', or the */
/*             number of rows of the matrix op( A ) and the number of */
/*             columns of the matrix B, for SIDE = 'R'.  N >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar alpha. When alpha is zero then R need not be */
/*             set before entry. */

/*     BETA    (input) DOUBLE PRECISION */
/*             The scalar beta. When beta is zero then A and B are not */
/*             referenced. */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR,M) */
/*             On entry with UPLO = 'U', the leading M-by-M upper */
/*             triangular part of this array must contain the upper */
/*             triangular part of the matrix R; the strictly lower */
/*             triangular part of the array is not referenced. */
/*             On entry with UPLO = 'L', the leading M-by-M lower */
/*             triangular part of this array must contain the lower */
/*             triangular part of the matrix R; the strictly upper */
/*             triangular part of the array is not referenced. */
/*             On exit, the leading M-by-M upper triangular part (if */
/*             UPLO = 'U'), or lower triangular part (if UPLO = 'L') of */
/*             this array contains the corresponding triangular part of */
/*                                 _ */
/*             the computed matrix R. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,M). */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,k), where */
/*             k = N  when  SIDE = 'L', and TRANS = 'N', or */
/*                          SIDE = 'R', and TRANS = 'T'; */
/*             k = M  when  SIDE = 'R', and TRANS = 'N', or */
/*                          SIDE = 'L', and TRANS = 'T'. */
/*             On entry, if SIDE = 'L', and TRANS = 'N', or */
/*                          SIDE = 'R', and TRANS = 'T', */
/*             the leading M-by-N part of this array must contain the */
/*             matrix A. */
/*             On entry, if SIDE = 'R', and TRANS = 'N', or */
/*                          SIDE = 'L', and TRANS = 'T', */
/*             the leading N-by-M part of this array must contain the */
/*             matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,l), where */
/*             l = M  when  SIDE = 'L', and TRANS = 'N', or */
/*                          SIDE = 'R', and TRANS = 'T'; */
/*             l = N  when  SIDE = 'R', and TRANS = 'N', or */
/*                          SIDE = 'L', and TRANS = 'T'. */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,p), where */
/*             p = M  when  SIDE = 'L'; */
/*             p = N  when  SIDE = 'R'. */
/*             On entry, the leading N-by-M part, if SIDE = 'L', or */
/*             M-by-N part, if SIDE = 'R', of this array must contain the */
/*             matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B. */
/*             LDB >= MAX(1,N), if SIDE = 'L'; */
/*             LDB >= MAX(1,M), if SIDE = 'R'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The matrix expression is evaluated taking the triangular */
/*     structure into account. BLAS 2 operations are used. A block */
/*     algorithm can be easily constructed; it can use BLAS 3 GEMM */
/*     operations for most computations, and calls of this BLAS 2 */
/*     algorithm for computing the triangles. */

/*     FURTHER COMMENTS */

/*     The main application of this routine is when the result should */
/*     be a symmetric matrix, e.g., when B = X*op( A )', for (1), or */
/*     B = op( A )'*X, for (2), where B is already available and X = X'. */

/*     CONTRIBUTORS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1999. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2004. */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix algebra, matrix operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

#line 196 "MB01RX.f"
    /* Parameter adjustments */
#line 196 "MB01RX.f"
    r_dim1 = *ldr;
#line 196 "MB01RX.f"
    r_offset = 1 + r_dim1;
#line 196 "MB01RX.f"
    r__ -= r_offset;
#line 196 "MB01RX.f"
    a_dim1 = *lda;
#line 196 "MB01RX.f"
    a_offset = 1 + a_dim1;
#line 196 "MB01RX.f"
    a -= a_offset;
#line 196 "MB01RX.f"
    b_dim1 = *ldb;
#line 196 "MB01RX.f"
    b_offset = 1 + b_dim1;
#line 196 "MB01RX.f"
    b -= b_offset;
#line 196 "MB01RX.f"

#line 196 "MB01RX.f"
    /* Function Body */
#line 196 "MB01RX.f"
    *info = 0;
#line 197 "MB01RX.f"
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 198 "MB01RX.f"
    luplo = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 199 "MB01RX.f"
    ltrans = lsame_(trans, "T", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1);

#line 201 "MB01RX.f"
    if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 202 "MB01RX.f"
	*info = -1;
#line 203 "MB01RX.f"
    } else if (! luplo && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 204 "MB01RX.f"
	*info = -2;
#line 205 "MB01RX.f"
    } else if (! ltrans && ! lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 206 "MB01RX.f"
	*info = -3;
#line 207 "MB01RX.f"
    } else if (*m < 0) {
#line 208 "MB01RX.f"
	*info = -4;
#line 209 "MB01RX.f"
    } else if (*n < 0) {
#line 210 "MB01RX.f"
	*info = -5;
#line 211 "MB01RX.f"
    } else if (*ldr < max(1,*m)) {
#line 212 "MB01RX.f"
	*info = -9;
#line 213 "MB01RX.f"
    } else if (*lda < 1 || (lside && ! ltrans || ! lside && ltrans) && *lda < 
	    *m || (lside && ltrans || ! lside && ! ltrans) && *lda < *n) {
#line 218 "MB01RX.f"
	*info = -11;
#line 219 "MB01RX.f"
    } else if (*ldb < 1 || lside && *ldb < *n || ! lside && *ldb < *m) {
#line 222 "MB01RX.f"
	*info = -13;
#line 223 "MB01RX.f"
    }

#line 225 "MB01RX.f"
    if (*info != 0) {

/*        Error return. */

#line 229 "MB01RX.f"
	i__1 = -(*info);
#line 229 "MB01RX.f"
	xerbla_("MB01RX", &i__1, (ftnlen)6);
#line 230 "MB01RX.f"
	return 0;
#line 231 "MB01RX.f"
    }

/*     Quick return if possible. */

#line 235 "MB01RX.f"
    if (*m == 0) {
#line 235 "MB01RX.f"
	return 0;
#line 235 "MB01RX.f"
    }

#line 238 "MB01RX.f"
    if (*beta == 0. || *n == 0) {
#line 239 "MB01RX.f"
	if (*alpha == 0.) {

/*           Special case alpha = 0. */

#line 243 "MB01RX.f"
	    dlaset_(uplo, m, m, &c_b10, &c_b10, &r__[r_offset], ldr, (ftnlen)
		    1);
#line 244 "MB01RX.f"
	} else {

/*           Special case beta = 0 or N = 0. */

#line 248 "MB01RX.f"
	    if (*alpha != 1.) {
#line 248 "MB01RX.f"
		dlascl_(uplo, &c__0, &c__0, &c_b14, alpha, m, m, &r__[
			r_offset], ldr, info, (ftnlen)1);
#line 248 "MB01RX.f"
	    }
#line 250 "MB01RX.f"
	}
#line 251 "MB01RX.f"
	return 0;
#line 252 "MB01RX.f"
    }

/*     General case: beta <> 0. */
/*     Compute the required triangle of (1) or (2) using BLAS 2 */
/*     operations. */

#line 258 "MB01RX.f"
    if (lside) {
#line 259 "MB01RX.f"
	if (luplo) {
#line 260 "MB01RX.f"
	    if (ltrans) {
#line 261 "MB01RX.f"
		i__1 = *m;
#line 261 "MB01RX.f"
		for (j = 1; j <= i__1; ++j) {
#line 262 "MB01RX.f"
		    dgemv_(trans, n, &j, beta, &a[a_offset], lda, &b[j * 
			    b_dim1 + 1], &c__1, alpha, &r__[j * r_dim1 + 1], &
			    c__1, (ftnlen)1);
#line 264 "MB01RX.f"
/* L10: */
#line 264 "MB01RX.f"
		}
#line 265 "MB01RX.f"
	    } else {
#line 266 "MB01RX.f"
		i__1 = *m;
#line 266 "MB01RX.f"
		for (j = 1; j <= i__1; ++j) {
#line 267 "MB01RX.f"
		    dgemv_(trans, &j, n, beta, &a[a_offset], lda, &b[j * 
			    b_dim1 + 1], &c__1, alpha, &r__[j * r_dim1 + 1], &
			    c__1, (ftnlen)1);
#line 269 "MB01RX.f"
/* L20: */
#line 269 "MB01RX.f"
		}
#line 270 "MB01RX.f"
	    }
#line 271 "MB01RX.f"
	} else {
#line 272 "MB01RX.f"
	    if (ltrans) {
#line 273 "MB01RX.f"
		i__1 = *m;
#line 273 "MB01RX.f"
		for (j = 1; j <= i__1; ++j) {
#line 274 "MB01RX.f"
		    i__2 = *m - j + 1;
#line 274 "MB01RX.f"
		    dgemv_(trans, n, &i__2, beta, &a[j * a_dim1 + 1], lda, &b[
			    j * b_dim1 + 1], &c__1, alpha, &r__[j + j * 
			    r_dim1], &c__1, (ftnlen)1);
#line 276 "MB01RX.f"
/* L30: */
#line 276 "MB01RX.f"
		}
#line 277 "MB01RX.f"
	    } else {
#line 278 "MB01RX.f"
		i__1 = *m;
#line 278 "MB01RX.f"
		for (j = 1; j <= i__1; ++j) {
#line 279 "MB01RX.f"
		    i__2 = *m - j + 1;
#line 279 "MB01RX.f"
		    dgemv_(trans, &i__2, n, beta, &a[j + a_dim1], lda, &b[j * 
			    b_dim1 + 1], &c__1, alpha, &r__[j + j * r_dim1], &
			    c__1, (ftnlen)1);
#line 281 "MB01RX.f"
/* L40: */
#line 281 "MB01RX.f"
		}
#line 282 "MB01RX.f"
	    }
#line 283 "MB01RX.f"
	}

#line 285 "MB01RX.f"
    } else {
#line 286 "MB01RX.f"
	if (luplo) {
#line 287 "MB01RX.f"
	    if (ltrans) {
#line 288 "MB01RX.f"
		i__1 = *m;
#line 288 "MB01RX.f"
		for (j = 1; j <= i__1; ++j) {
#line 289 "MB01RX.f"
		    dgemv_("NoTranspose", &j, n, beta, &b[b_offset], ldb, &a[
			    j + a_dim1], lda, alpha, &r__[j * r_dim1 + 1], &
			    c__1, (ftnlen)11);
#line 291 "MB01RX.f"
/* L50: */
#line 291 "MB01RX.f"
		}
#line 292 "MB01RX.f"
	    } else {
#line 293 "MB01RX.f"
		i__1 = *m;
#line 293 "MB01RX.f"
		for (j = 1; j <= i__1; ++j) {
#line 294 "MB01RX.f"
		    dgemv_("NoTranspose", &j, n, beta, &b[b_offset], ldb, &a[
			    j * a_dim1 + 1], &c__1, alpha, &r__[j * r_dim1 + 
			    1], &c__1, (ftnlen)11);
#line 296 "MB01RX.f"
/* L60: */
#line 296 "MB01RX.f"
		}
#line 297 "MB01RX.f"
	    }
#line 298 "MB01RX.f"
	} else {
#line 299 "MB01RX.f"
	    if (ltrans) {
#line 300 "MB01RX.f"
		i__1 = *m;
#line 300 "MB01RX.f"
		for (j = 1; j <= i__1; ++j) {
#line 301 "MB01RX.f"
		    i__2 = *m - j + 1;
#line 301 "MB01RX.f"
		    dgemv_("NoTranspose", &i__2, n, beta, &b[j + b_dim1], ldb,
			     &a[j + a_dim1], lda, alpha, &r__[j + j * r_dim1],
			     &c__1, (ftnlen)11);
#line 303 "MB01RX.f"
/* L70: */
#line 303 "MB01RX.f"
		}
#line 304 "MB01RX.f"
	    } else {
#line 305 "MB01RX.f"
		i__1 = *m;
#line 305 "MB01RX.f"
		for (j = 1; j <= i__1; ++j) {
#line 306 "MB01RX.f"
		    i__2 = *m - j + 1;
#line 306 "MB01RX.f"
		    dgemv_("NoTranspose", &i__2, n, beta, &b[j + b_dim1], ldb,
			     &a[j * a_dim1 + 1], &c__1, alpha, &r__[j + j * 
			    r_dim1], &c__1, (ftnlen)11);
#line 308 "MB01RX.f"
/* L80: */
#line 308 "MB01RX.f"
		}
#line 309 "MB01RX.f"
	    }
#line 310 "MB01RX.f"
	}
#line 311 "MB01RX.f"
    }

#line 313 "MB01RX.f"
    return 0;
/* *** Last line of MB01RX *** */
} /* mb01rx_ */

