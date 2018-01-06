#line 1 "MB01RY.f"
/* MB01RY.f -- translated by f2c (version 20100827).
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

#line 1 "MB01RY.f"
/* Table of constant values */

static doublereal c_b10 = 0.;
static integer c__0 = 0;
static doublereal c_b14 = 1.;
static integer c__1 = 1;

/* Subroutine */ int mb01ry_(char *side, char *uplo, char *trans, integer *m, 
	doublereal *alpha, doublereal *beta, doublereal *r__, integer *ldr, 
	doublereal *h__, integer *ldh, doublereal *b, integer *ldb, 
	doublereal *dwork, integer *info, ftnlen side_len, ftnlen uplo_len, 
	ftnlen trans_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, h_dim1, h_offset, r_dim1, r_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, j;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), dswap_(integer 
	    *, doublereal *, integer *, doublereal *, integer *);
    static logical luplo;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen), dlascl_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
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
/*        R = alpha*R + beta*op( H )*B,                               (1) */
/*        _ */
/*        R = alpha*R + beta*B*op( H ),                               (2) */
/*                                                    _ */
/*     where alpha and beta are scalars, H, B, R, and R are m-by-m */
/*     matrices, H is an upper Hessenberg matrix, and op( H ) is one of */

/*        op( H ) = H   or   op( H ) = H',  the transpose of H. */

/*     The result is overwritten on R. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SIDE    CHARACTER*1 */
/*             Specifies whether the Hessenberg matrix H appears on the */
/*             left or right in the matrix product as follows: */
/*                     _ */
/*             = 'L':  R = alpha*R + beta*op( H )*B; */
/*                     _ */
/*             = 'R':  R = alpha*R + beta*B*op( H ). */

/*     UPLO    CHARACTER*1                               _ */
/*             Specifies which triangles of the matrices R and R are */
/*             computed and given, respectively, as follows: */
/*             = 'U':  the upper triangular part; */
/*             = 'L':  the lower triangular part. */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of op( H ) to be used in the matrix */
/*             multiplication as follows: */
/*             = 'N':  op( H ) = H; */
/*             = 'T':  op( H ) = H'; */
/*             = 'C':  op( H ) = H'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER           _ */
/*             The order of the matrices R, R, H and B.  M >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar alpha. When alpha is zero then R need not be */
/*             set before entry. */

/*     BETA    (input) DOUBLE PRECISION */
/*             The scalar beta. When beta is zero then H and B are not */
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

/*     H       (input) DOUBLE PRECISION array, dimension (LDH,M) */
/*             On entry, the leading M-by-M upper Hessenberg part of */
/*             this array must contain the upper Hessenberg part of the */
/*             matrix H. */
/*             The elements below the subdiagonal are not referenced, */
/*             except possibly for those in the first column, which */
/*             could be overwritten, but are restored on exit. */

/*     LDH     INTEGER */
/*             The leading dimension of array H.  LDH >= MAX(1,M). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading M-by-M part of this array must */
/*             contain the matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,M). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             LDWORK >= M, if  beta <> 0 and SIDE = 'L'; */
/*             LDWORK >= 0, if  beta =  0 or  SIDE = 'R'. */
/*             This array is not referenced when beta = 0 or SIDE = 'R'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The matrix expression is efficiently evaluated taking the */
/*     Hessenberg/triangular structure into account. BLAS 2 operations */
/*     are used. A block algorithm can be constructed; it can use BLAS 3 */
/*     GEMM operations for most computations, and calls of this BLAS 2 */
/*     algorithm for computing the triangles. */

/*     FURTHER COMMENTS */

/*     The main application of this routine is when the result should */
/*     be a symmetric matrix, e.g., when B = X*op( H )', for (1), or */
/*     B = op( H )'*X, for (2), where B is already available and X = X'. */

/*     CONTRIBUTORS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1999. */

/*     REVISIONS */

/*     - */

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

#line 179 "MB01RY.f"
    /* Parameter adjustments */
#line 179 "MB01RY.f"
    r_dim1 = *ldr;
#line 179 "MB01RY.f"
    r_offset = 1 + r_dim1;
#line 179 "MB01RY.f"
    r__ -= r_offset;
#line 179 "MB01RY.f"
    h_dim1 = *ldh;
#line 179 "MB01RY.f"
    h_offset = 1 + h_dim1;
#line 179 "MB01RY.f"
    h__ -= h_offset;
#line 179 "MB01RY.f"
    b_dim1 = *ldb;
#line 179 "MB01RY.f"
    b_offset = 1 + b_dim1;
#line 179 "MB01RY.f"
    b -= b_offset;
#line 179 "MB01RY.f"
    --dwork;
#line 179 "MB01RY.f"

#line 179 "MB01RY.f"
    /* Function Body */
#line 179 "MB01RY.f"
    *info = 0;
#line 180 "MB01RY.f"
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 181 "MB01RY.f"
    luplo = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 182 "MB01RY.f"
    ltrans = lsame_(trans, "T", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1);

#line 184 "MB01RY.f"
    if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 185 "MB01RY.f"
	*info = -1;
#line 186 "MB01RY.f"
    } else if (! luplo && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 187 "MB01RY.f"
	*info = -2;
#line 188 "MB01RY.f"
    } else if (! ltrans && ! lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 189 "MB01RY.f"
	*info = -3;
#line 190 "MB01RY.f"
    } else if (*m < 0) {
#line 191 "MB01RY.f"
	*info = -4;
#line 192 "MB01RY.f"
    } else if (*ldr < max(1,*m)) {
#line 193 "MB01RY.f"
	*info = -8;
#line 194 "MB01RY.f"
    } else if (*ldh < max(1,*m)) {
#line 195 "MB01RY.f"
	*info = -10;
#line 196 "MB01RY.f"
    } else if (*ldb < max(1,*m)) {
#line 197 "MB01RY.f"
	*info = -12;
#line 198 "MB01RY.f"
    }

#line 200 "MB01RY.f"
    if (*info != 0) {

/*        Error return. */

#line 204 "MB01RY.f"
	i__1 = -(*info);
#line 204 "MB01RY.f"
	xerbla_("MB01RY", &i__1, (ftnlen)6);
#line 205 "MB01RY.f"
	return 0;
#line 206 "MB01RY.f"
    }

/*     Quick return if possible. */

#line 210 "MB01RY.f"
    if (*m == 0) {
#line 210 "MB01RY.f"
	return 0;
#line 210 "MB01RY.f"
    }

#line 213 "MB01RY.f"
    if (*beta == 0.) {
#line 214 "MB01RY.f"
	if (*alpha == 0.) {

/*           Special case when both alpha = 0 and beta = 0. */

#line 218 "MB01RY.f"
	    dlaset_(uplo, m, m, &c_b10, &c_b10, &r__[r_offset], ldr, (ftnlen)
		    1);
#line 219 "MB01RY.f"
	} else {

/*           Special case beta = 0. */

#line 223 "MB01RY.f"
	    if (*alpha != 1.) {
#line 223 "MB01RY.f"
		dlascl_(uplo, &c__0, &c__0, &c_b14, alpha, m, m, &r__[
			r_offset], ldr, info, (ftnlen)1);
#line 223 "MB01RY.f"
	    }
#line 225 "MB01RY.f"
	}
#line 226 "MB01RY.f"
	return 0;
#line 227 "MB01RY.f"
    }

/*     General case: beta <> 0. */
/*     Compute the required triangle of (1) or (2) using BLAS 2 */
/*     operations. */

#line 233 "MB01RY.f"
    if (lside) {

/*        To avoid repeated references to the subdiagonal elements of H, */
/*        these are swapped with the corresponding elements of H in the */
/*        first column, and are finally restored. */

#line 239 "MB01RY.f"
	if (*m > 2) {
#line 239 "MB01RY.f"
	    i__1 = *m - 2;
#line 239 "MB01RY.f"
	    i__2 = *ldh + 1;
#line 239 "MB01RY.f"
	    dswap_(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3], &
		    c__1);
#line 239 "MB01RY.f"
	}

#line 242 "MB01RY.f"
	if (luplo) {
#line 243 "MB01RY.f"
	    if (ltrans) {

#line 245 "MB01RY.f"
		i__1 = *m;
#line 245 "MB01RY.f"
		for (j = 1; j <= i__1; ++j) {

/*                 Multiply the transposed upper triangle of the leading */
/*                 j-by-j submatrix of H by the leading part of the j-th */
/*                 column of B. */

#line 251 "MB01RY.f"
		    dcopy_(&j, &b[j * b_dim1 + 1], &c__1, &dwork[1], &c__1);
#line 252 "MB01RY.f"
		    dtrmv_("Upper", trans, "Non-unit", &j, &h__[h_offset], 
			    ldh, &dwork[1], &c__1, (ftnlen)5, (ftnlen)1, (
			    ftnlen)8);

/*                 Add the contribution of the subdiagonal of H to */
/*                 the j-th column of the product. */

/* Computing MIN */
#line 258 "MB01RY.f"
		    i__3 = j, i__4 = *m - 1;
#line 258 "MB01RY.f"
		    i__2 = min(i__3,i__4);
#line 258 "MB01RY.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 259 "MB01RY.f"
			r__[i__ + j * r_dim1] = *alpha * r__[i__ + j * r_dim1]
				 + *beta * (dwork[i__] + h__[i__ + 1 + h_dim1]
				 * b[i__ + 1 + j * b_dim1]);
#line 261 "MB01RY.f"
/* L10: */
#line 261 "MB01RY.f"
		    }

#line 263 "MB01RY.f"
/* L20: */
#line 263 "MB01RY.f"
		}

#line 265 "MB01RY.f"
		r__[*m + *m * r_dim1] = *alpha * r__[*m + *m * r_dim1] + *
			beta * dwork[*m];

#line 267 "MB01RY.f"
	    } else {

#line 269 "MB01RY.f"
		i__1 = *m;
#line 269 "MB01RY.f"
		for (j = 1; j <= i__1; ++j) {

/*                 Multiply the upper triangle of the leading j-by-j */
/*                 submatrix of H by the leading part of the j-th column */
/*                 of B. */

#line 275 "MB01RY.f"
		    dcopy_(&j, &b[j * b_dim1 + 1], &c__1, &dwork[1], &c__1);
#line 276 "MB01RY.f"
		    dtrmv_("Upper", trans, "Non-unit", &j, &h__[h_offset], 
			    ldh, &dwork[1], &c__1, (ftnlen)5, (ftnlen)1, (
			    ftnlen)8);
#line 278 "MB01RY.f"
		    if (j < *m) {

/*                    Multiply the remaining right part of the leading */
/*                    j-by-M submatrix of H by the trailing part of the */
/*                    j-th column of B. */

#line 284 "MB01RY.f"
			i__2 = *m - j;
#line 284 "MB01RY.f"
			dgemv_(trans, &j, &i__2, beta, &h__[(j + 1) * h_dim1 
				+ 1], ldh, &b[j + 1 + j * b_dim1], &c__1, 
				alpha, &r__[j * r_dim1 + 1], &c__1, (ftnlen)1)
				;
#line 286 "MB01RY.f"
		    } else {
#line 287 "MB01RY.f"
			dscal_(m, alpha, &r__[*m * r_dim1 + 1], &c__1);
#line 288 "MB01RY.f"
		    }

/*                 Add the contribution of the subdiagonal of H to */
/*                 the j-th column of the product. */

#line 293 "MB01RY.f"
		    r__[j * r_dim1 + 1] += *beta * dwork[1];

#line 295 "MB01RY.f"
		    i__2 = j;
#line 295 "MB01RY.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 296 "MB01RY.f"
			r__[i__ + j * r_dim1] += *beta * (dwork[i__] + h__[
				i__ + h_dim1] * b[i__ - 1 + j * b_dim1]);
#line 298 "MB01RY.f"
/* L30: */
#line 298 "MB01RY.f"
		    }

#line 300 "MB01RY.f"
/* L40: */
#line 300 "MB01RY.f"
		}

#line 302 "MB01RY.f"
	    }

#line 304 "MB01RY.f"
	} else {

#line 306 "MB01RY.f"
	    if (ltrans) {

#line 308 "MB01RY.f"
		for (j = *m; j >= 1; --j) {

/*                 Multiply the transposed upper triangle of the trailing */
/*                 (M-j+1)-by-(M-j+1) submatrix of H by the trailing part */
/*                 of the j-th column of B. */

#line 314 "MB01RY.f"
		    i__1 = *m - j + 1;
#line 314 "MB01RY.f"
		    dcopy_(&i__1, &b[j + j * b_dim1], &c__1, &dwork[j], &c__1)
			    ;
#line 315 "MB01RY.f"
		    i__1 = *m - j + 1;
#line 315 "MB01RY.f"
		    dtrmv_("Upper", trans, "Non-unit", &i__1, &h__[j + j * 
			    h_dim1], ldh, &dwork[j], &c__1, (ftnlen)5, (
			    ftnlen)1, (ftnlen)8);
#line 317 "MB01RY.f"
		    if (j > 1) {

/*                    Multiply the remaining left part of the trailing */
/*                    (M-j+1)-by-(j-1) submatrix of H' by the leading */
/*                    part of the j-th column of B. */

#line 323 "MB01RY.f"
			i__1 = j - 1;
#line 323 "MB01RY.f"
			i__2 = *m - j + 1;
#line 323 "MB01RY.f"
			dgemv_(trans, &i__1, &i__2, beta, &h__[j * h_dim1 + 1]
				, ldh, &b[j * b_dim1 + 1], &c__1, alpha, &r__[
				j + j * r_dim1], &c__1, (ftnlen)1);
#line 326 "MB01RY.f"
		    } else {
#line 327 "MB01RY.f"
			dscal_(m, alpha, &r__[r_dim1 + 1], &c__1);
#line 328 "MB01RY.f"
		    }

/*                 Add the contribution of the subdiagonal of H to */
/*                 the j-th column of the product. */

#line 333 "MB01RY.f"
		    i__1 = *m - 1;
#line 333 "MB01RY.f"
		    for (i__ = j; i__ <= i__1; ++i__) {
#line 334 "MB01RY.f"
			r__[i__ + j * r_dim1] += *beta * (dwork[i__] + h__[
				i__ + 1 + h_dim1] * b[i__ + 1 + j * b_dim1]);
#line 336 "MB01RY.f"
/* L50: */
#line 336 "MB01RY.f"
		    }

#line 338 "MB01RY.f"
		    r__[*m + j * r_dim1] += *beta * dwork[*m];
#line 339 "MB01RY.f"
/* L60: */
#line 339 "MB01RY.f"
		}

#line 341 "MB01RY.f"
	    } else {

#line 343 "MB01RY.f"
		for (j = *m; j >= 1; --j) {

/*                 Multiply the upper triangle of the trailing */
/*                 (M-j+1)-by-(M-j+1) submatrix of H by the trailing */
/*                 part of the j-th column of B. */

#line 349 "MB01RY.f"
		    i__1 = *m - j + 1;
#line 349 "MB01RY.f"
		    dcopy_(&i__1, &b[j + j * b_dim1], &c__1, &dwork[j], &c__1)
			    ;
#line 350 "MB01RY.f"
		    i__1 = *m - j + 1;
#line 350 "MB01RY.f"
		    dtrmv_("Upper", trans, "Non-unit", &i__1, &h__[j + j * 
			    h_dim1], ldh, &dwork[j], &c__1, (ftnlen)5, (
			    ftnlen)1, (ftnlen)8);

/*                 Add the contribution of the subdiagonal of H to */
/*                 the j-th column of the product. */

#line 356 "MB01RY.f"
		    i__1 = *m;
#line 356 "MB01RY.f"
		    for (i__ = max(j,2); i__ <= i__1; ++i__) {
#line 357 "MB01RY.f"
			r__[i__ + j * r_dim1] = *alpha * r__[i__ + j * r_dim1]
				 + *beta * (dwork[i__] + h__[i__ + h_dim1] * 
				b[i__ - 1 + j * b_dim1]);
#line 359 "MB01RY.f"
/* L70: */
#line 359 "MB01RY.f"
		    }

#line 361 "MB01RY.f"
/* L80: */
#line 361 "MB01RY.f"
		}

#line 363 "MB01RY.f"
		r__[r_dim1 + 1] = *alpha * r__[r_dim1 + 1] + *beta * dwork[1];

#line 365 "MB01RY.f"
	    }
#line 366 "MB01RY.f"
	}

#line 368 "MB01RY.f"
	if (*m > 2) {
#line 368 "MB01RY.f"
	    i__1 = *m - 2;
#line 368 "MB01RY.f"
	    i__2 = *ldh + 1;
#line 368 "MB01RY.f"
	    dswap_(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3], &
		    c__1);
#line 368 "MB01RY.f"
	}

#line 371 "MB01RY.f"
    } else {

/*        Row-wise calculations are used for H, if SIDE = 'R' and */
/*        TRANS = 'T'. */

#line 376 "MB01RY.f"
	if (luplo) {
#line 377 "MB01RY.f"
	    if (ltrans) {
#line 378 "MB01RY.f"
		r__[r_dim1 + 1] = *alpha * r__[r_dim1 + 1] + *beta * ddot_(m, 
			&b[b_offset], ldb, &h__[h_offset], ldh);

#line 381 "MB01RY.f"
		i__1 = *m;
#line 381 "MB01RY.f"
		for (j = 2; j <= i__1; ++j) {
#line 382 "MB01RY.f"
		    i__2 = *m - j + 2;
#line 382 "MB01RY.f"
		    dgemv_("NoTranspose", &j, &i__2, beta, &b[(j - 1) * 
			    b_dim1 + 1], ldb, &h__[j + (j - 1) * h_dim1], ldh,
			     alpha, &r__[j * r_dim1 + 1], &c__1, (ftnlen)11);
#line 385 "MB01RY.f"
/* L90: */
#line 385 "MB01RY.f"
		}

#line 387 "MB01RY.f"
	    } else {

#line 389 "MB01RY.f"
		i__1 = *m - 1;
#line 389 "MB01RY.f"
		for (j = 1; j <= i__1; ++j) {
#line 390 "MB01RY.f"
		    i__2 = j + 1;
#line 390 "MB01RY.f"
		    dgemv_("NoTranspose", &j, &i__2, beta, &b[b_offset], ldb, 
			    &h__[j * h_dim1 + 1], &c__1, alpha, &r__[j * 
			    r_dim1 + 1], &c__1, (ftnlen)11);
#line 392 "MB01RY.f"
/* L100: */
#line 392 "MB01RY.f"
		}

#line 394 "MB01RY.f"
		dgemv_("NoTranspose", m, m, beta, &b[b_offset], ldb, &h__[*m *
			 h_dim1 + 1], &c__1, alpha, &r__[*m * r_dim1 + 1], &
			c__1, (ftnlen)11);

#line 397 "MB01RY.f"
	    }

#line 399 "MB01RY.f"
	} else {

#line 401 "MB01RY.f"
	    if (ltrans) {

#line 403 "MB01RY.f"
		dgemv_("NoTranspose", m, m, beta, &b[b_offset], ldb, &h__[
			h_offset], ldh, alpha, &r__[r_dim1 + 1], &c__1, (
			ftnlen)11);

#line 406 "MB01RY.f"
		i__1 = *m;
#line 406 "MB01RY.f"
		for (j = 2; j <= i__1; ++j) {
#line 407 "MB01RY.f"
		    i__2 = *m - j + 1;
#line 407 "MB01RY.f"
		    i__3 = *m - j + 2;
#line 407 "MB01RY.f"
		    dgemv_("NoTranspose", &i__2, &i__3, beta, &b[j + (j - 1) *
			     b_dim1], ldb, &h__[j + (j - 1) * h_dim1], ldh, 
			    alpha, &r__[j + j * r_dim1], &c__1, (ftnlen)11);
#line 410 "MB01RY.f"
/* L110: */
#line 410 "MB01RY.f"
		}

#line 412 "MB01RY.f"
	    } else {

#line 414 "MB01RY.f"
		i__1 = *m - 1;
#line 414 "MB01RY.f"
		for (j = 1; j <= i__1; ++j) {
#line 415 "MB01RY.f"
		    i__2 = *m - j + 1;
#line 415 "MB01RY.f"
		    i__3 = j + 1;
#line 415 "MB01RY.f"
		    dgemv_("NoTranspose", &i__2, &i__3, beta, &b[j + b_dim1], 
			    ldb, &h__[j * h_dim1 + 1], &c__1, alpha, &r__[j + 
			    j * r_dim1], &c__1, (ftnlen)11);
#line 418 "MB01RY.f"
/* L120: */
#line 418 "MB01RY.f"
		}

#line 420 "MB01RY.f"
		r__[*m + *m * r_dim1] = *alpha * r__[*m + *m * r_dim1] + *
			beta * ddot_(m, &b[*m + b_dim1], ldb, &h__[*m * 
			h_dim1 + 1], &c__1);

#line 423 "MB01RY.f"
	    }
#line 424 "MB01RY.f"
	}
#line 425 "MB01RY.f"
    }

#line 427 "MB01RY.f"
    return 0;
/* *** Last line of MB01RY *** */
} /* mb01ry_ */

