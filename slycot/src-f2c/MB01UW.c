#line 1 "MB01UW.f"
/* MB01UW.f -- translated by f2c (version 20100827).
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

#line 1 "MB01UW.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b10 = 0.;
static doublereal c_b33 = 1.;
static integer c__0 = 0;

/* Subroutine */ int mb01uw_(char *side, char *trans, integer *m, integer *n, 
	doublereal *alpha, doublereal *h__, integer *ldh, doublereal *a, 
	integer *lda, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, h_dim1, h_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, jw;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *), dtrmv_(char *, char *, char *, integer *, doublereal *
	    , integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen), 
	    dlascl_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
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

/*     To compute one of the matrix products */

/*        A : = alpha*op( H ) * A, or A : = alpha*A * op( H ), */

/*     where alpha is a scalar, A is an m-by-n matrix, H is an upper */
/*     Hessenberg matrix, and op( H ) is one of */

/*        op( H ) = H   or   op( H ) = H',  the transpose of H. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SIDE    CHARACTER*1 */
/*             Specifies whether the Hessenberg matrix H appears on the */
/*             left or right in the matrix product as follows: */
/*             = 'L':  A := alpha*op( H ) * A; */
/*             = 'R':  A := alpha*A * op( H ). */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of op( H ) to be used in the matrix */
/*             multiplication as follows: */
/*             = 'N':  op( H ) = H; */
/*             = 'T':  op( H ) = H'; */
/*             = 'C':  op( H ) = H'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar alpha. When alpha is zero then H is not */
/*             referenced and A need not be set before entry. */

/*     H       (input) DOUBLE PRECISION array, dimension (LDH,k) */
/*             where k is M when SIDE = 'L' and is N when SIDE = 'R'. */
/*             On entry with SIDE = 'L', the leading M-by-M upper */
/*             Hessenberg part of this array must contain the upper */
/*             Hessenberg matrix H. */
/*             On entry with SIDE = 'R', the leading N-by-N upper */
/*             Hessenberg part of this array must contain the upper */
/*             Hessenberg matrix H. */
/*             The elements below the subdiagonal are not referenced, */
/*             except possibly for those in the first column, which */
/*             could be overwritten, but are restored on exit. */

/*     LDH     INTEGER */
/*             The leading dimension of the array H.  LDH >= max(1,k), */
/*             where k is M when SIDE = 'L' and is N when SIDE = 'R'. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the matrix A. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the computed product. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, alpha <> 0, and LDWORK >= M*N > 0, */
/*             DWORK contains a copy of the matrix A, having the leading */
/*             dimension M. */
/*             This array is not referenced when alpha = 0. */

/*     LDWORK  The length of the array DWORK. */
/*             LDWORK >= 0,   if  alpha =  0 or MIN(M,N) = 0; */
/*             LDWORK >= M-1, if  SIDE  = 'L'; */
/*             LDWORK >= N-1, if  SIDE  = 'R'. */
/*             For maximal efficiency LDWORK should be at least M*N. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The required matrix product is computed in two steps. In the first */
/*     step, the upper triangle of H is used; in the second step, the */
/*     contribution of the subdiagonal is added. If the workspace can */
/*     accomodate a copy of A, a fast BLAS 3 DTRMM operation is used in */
/*     the first step. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, January 1999. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2004. */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix operations. */

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

#line 155 "MB01UW.f"
    /* Parameter adjustments */
#line 155 "MB01UW.f"
    h_dim1 = *ldh;
#line 155 "MB01UW.f"
    h_offset = 1 + h_dim1;
#line 155 "MB01UW.f"
    h__ -= h_offset;
#line 155 "MB01UW.f"
    a_dim1 = *lda;
#line 155 "MB01UW.f"
    a_offset = 1 + a_dim1;
#line 155 "MB01UW.f"
    a -= a_offset;
#line 155 "MB01UW.f"
    --dwork;
#line 155 "MB01UW.f"

#line 155 "MB01UW.f"
    /* Function Body */
#line 155 "MB01UW.f"
    *info = 0;
#line 156 "MB01UW.f"
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 157 "MB01UW.f"
    ltrans = lsame_(trans, "T", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1);

#line 159 "MB01UW.f"
    if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 160 "MB01UW.f"
	*info = -1;
#line 161 "MB01UW.f"
    } else if (! ltrans && ! lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 162 "MB01UW.f"
	*info = -2;
#line 163 "MB01UW.f"
    } else if (*m < 0) {
#line 164 "MB01UW.f"
	*info = -3;
#line 165 "MB01UW.f"
    } else if (*n < 0) {
#line 166 "MB01UW.f"
	*info = -4;
#line 167 "MB01UW.f"
    } else if (*ldh < 1 || lside && *ldh < *m || ! lside && *ldh < *n) {
#line 169 "MB01UW.f"
	*info = -7;
#line 170 "MB01UW.f"
    } else if (*lda < max(1,*m)) {
#line 171 "MB01UW.f"
	*info = -9;
#line 172 "MB01UW.f"
    } else if (*ldwork < 0 || *alpha != 0. && min(*m,*n) > 0 && (lside && *
	    ldwork < *m - 1 || ! lside && *ldwork < *n - 1)) {
#line 176 "MB01UW.f"
	*info = -11;
#line 177 "MB01UW.f"
    }

#line 179 "MB01UW.f"
    if (*info != 0) {

/*        Error return. */

#line 183 "MB01UW.f"
	i__1 = -(*info);
#line 183 "MB01UW.f"
	xerbla_("MB01UW", &i__1, (ftnlen)6);
#line 184 "MB01UW.f"
	return 0;
#line 185 "MB01UW.f"
    }

/*     Quick return, if possible. */

#line 189 "MB01UW.f"
    if (min(*m,*n) == 0) {
#line 190 "MB01UW.f"
	return 0;
#line 191 "MB01UW.f"
    } else if (lside) {
#line 192 "MB01UW.f"
	if (*m == 1) {
#line 193 "MB01UW.f"
	    d__1 = *alpha * h__[h_dim1 + 1];
#line 193 "MB01UW.f"
	    dscal_(n, &d__1, &a[a_offset], lda);
#line 194 "MB01UW.f"
	    return 0;
#line 195 "MB01UW.f"
	}
#line 196 "MB01UW.f"
    } else {
#line 197 "MB01UW.f"
	if (*n == 1) {
#line 198 "MB01UW.f"
	    d__1 = *alpha * h__[h_dim1 + 1];
#line 198 "MB01UW.f"
	    dscal_(m, &d__1, &a[a_offset], &c__1);
#line 199 "MB01UW.f"
	    return 0;
#line 200 "MB01UW.f"
	}
#line 201 "MB01UW.f"
    }

#line 203 "MB01UW.f"
    if (*alpha == 0.) {

/*        Set A to zero and return. */

#line 207 "MB01UW.f"
	dlaset_("Full", m, n, &c_b10, &c_b10, &a[a_offset], lda, (ftnlen)4);
#line 208 "MB01UW.f"
	return 0;
#line 209 "MB01UW.f"
    }

#line 211 "MB01UW.f"
    if (*ldwork >= *m * *n) {

/*        Enough workspace for a fast BLAS 3 calculation. */
/*        Save A in the workspace and compute one of the matrix products */
/*          A : = alpha*op( triu( H ) ) * A, or */
/*          A : = alpha*A * op( triu( H ) ), */
/*        involving the upper triangle of H. */

#line 219 "MB01UW.f"
	dlacpy_("Full", m, n, &a[a_offset], lda, &dwork[1], m, (ftnlen)4);
#line 220 "MB01UW.f"
	dtrmm_(side, "Upper", trans, "Non-unit", m, n, alpha, &h__[h_offset], 
		ldh, &a[a_offset], lda, (ftnlen)1, (ftnlen)5, (ftnlen)1, (
		ftnlen)8);

/*        Add the contribution of the subdiagonal of H. */
/*        If SIDE = 'L', the subdiagonal of H is swapped with the */
/*        corresponding elements in the first column of H, and the */
/*        calculations are organized for column operations. */

#line 228 "MB01UW.f"
	if (lside) {
#line 229 "MB01UW.f"
	    if (*m > 2) {
#line 229 "MB01UW.f"
		i__1 = *m - 2;
#line 229 "MB01UW.f"
		i__2 = *ldh + 1;
#line 229 "MB01UW.f"
		dswap_(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3]
			, &c__1);
#line 229 "MB01UW.f"
	    }
#line 231 "MB01UW.f"
	    if (ltrans) {
#line 232 "MB01UW.f"
		jw = 1;
#line 233 "MB01UW.f"
		i__1 = *n;
#line 233 "MB01UW.f"
		for (j = 1; j <= i__1; ++j) {
#line 234 "MB01UW.f"
		    ++jw;
#line 235 "MB01UW.f"
		    i__2 = *m - 1;
#line 235 "MB01UW.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 236 "MB01UW.f"
			a[i__ + j * a_dim1] += *alpha * h__[i__ + 1 + h_dim1] 
				* dwork[jw];
#line 238 "MB01UW.f"
			++jw;
#line 239 "MB01UW.f"
/* L10: */
#line 239 "MB01UW.f"
		    }
#line 240 "MB01UW.f"
/* L20: */
#line 240 "MB01UW.f"
		}
#line 241 "MB01UW.f"
	    } else {
#line 242 "MB01UW.f"
		jw = 0;
#line 243 "MB01UW.f"
		i__1 = *n;
#line 243 "MB01UW.f"
		for (j = 1; j <= i__1; ++j) {
#line 244 "MB01UW.f"
		    ++jw;
#line 245 "MB01UW.f"
		    i__2 = *m;
#line 245 "MB01UW.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 246 "MB01UW.f"
			a[i__ + j * a_dim1] += *alpha * h__[i__ + h_dim1] * 
				dwork[jw];
#line 248 "MB01UW.f"
			++jw;
#line 249 "MB01UW.f"
/* L30: */
#line 249 "MB01UW.f"
		    }
#line 250 "MB01UW.f"
/* L40: */
#line 250 "MB01UW.f"
		}
#line 251 "MB01UW.f"
	    }
#line 252 "MB01UW.f"
	    if (*m > 2) {
#line 252 "MB01UW.f"
		i__1 = *m - 2;
#line 252 "MB01UW.f"
		i__2 = *ldh + 1;
#line 252 "MB01UW.f"
		dswap_(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3]
			, &c__1);
#line 252 "MB01UW.f"
	    }

#line 255 "MB01UW.f"
	} else {

#line 257 "MB01UW.f"
	    if (ltrans) {
#line 258 "MB01UW.f"
		jw = 1;
#line 259 "MB01UW.f"
		i__1 = *n - 1;
#line 259 "MB01UW.f"
		for (j = 1; j <= i__1; ++j) {
#line 260 "MB01UW.f"
		    if (h__[j + 1 + j * h_dim1] != 0.) {
#line 260 "MB01UW.f"
			d__1 = *alpha * h__[j + 1 + j * h_dim1];
#line 260 "MB01UW.f"
			daxpy_(m, &d__1, &dwork[jw], &c__1, &a[(j + 1) * 
				a_dim1 + 1], &c__1);
#line 260 "MB01UW.f"
		    }
#line 263 "MB01UW.f"
		    jw += *m;
#line 264 "MB01UW.f"
/* L50: */
#line 264 "MB01UW.f"
		}
#line 265 "MB01UW.f"
	    } else {
#line 266 "MB01UW.f"
		jw = *m + 1;
#line 267 "MB01UW.f"
		i__1 = *n - 1;
#line 267 "MB01UW.f"
		for (j = 1; j <= i__1; ++j) {
#line 268 "MB01UW.f"
		    if (h__[j + 1 + j * h_dim1] != 0.) {
#line 268 "MB01UW.f"
			d__1 = *alpha * h__[j + 1 + j * h_dim1];
#line 268 "MB01UW.f"
			daxpy_(m, &d__1, &dwork[jw], &c__1, &a[j * a_dim1 + 1]
				, &c__1);
#line 268 "MB01UW.f"
		    }
#line 271 "MB01UW.f"
		    jw += *m;
#line 272 "MB01UW.f"
/* L60: */
#line 272 "MB01UW.f"
		}
#line 273 "MB01UW.f"
	    }
#line 274 "MB01UW.f"
	}

#line 276 "MB01UW.f"
    } else {

/*        Use a BLAS 2 calculation. */

#line 280 "MB01UW.f"
	if (lside) {
#line 281 "MB01UW.f"
	    if (*m > 2) {
#line 281 "MB01UW.f"
		i__1 = *m - 2;
#line 281 "MB01UW.f"
		i__2 = *ldh + 1;
#line 281 "MB01UW.f"
		dswap_(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3]
			, &c__1);
#line 281 "MB01UW.f"
	    }
#line 283 "MB01UW.f"
	    if (ltrans) {
#line 284 "MB01UW.f"
		i__1 = *n;
#line 284 "MB01UW.f"
		for (j = 1; j <= i__1; ++j) {

/*                 Compute the contribution of the subdiagonal of H to */
/*                 the j-th column of the product. */

#line 289 "MB01UW.f"
		    i__2 = *m - 1;
#line 289 "MB01UW.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 290 "MB01UW.f"
			dwork[i__] = h__[i__ + 1 + h_dim1] * a[i__ + 1 + j * 
				a_dim1];
#line 291 "MB01UW.f"
/* L70: */
#line 291 "MB01UW.f"
		    }

/*                 Multiply the upper triangle of H by the j-th column */
/*                 of A, and add to the above result. */

#line 296 "MB01UW.f"
		    dtrmv_("Upper", trans, "Non-unit", m, &h__[h_offset], ldh,
			     &a[j * a_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)1, 
			    (ftnlen)8);
#line 298 "MB01UW.f"
		    i__2 = *m - 1;
#line 298 "MB01UW.f"
		    daxpy_(&i__2, &c_b33, &dwork[1], &c__1, &a[j * a_dim1 + 1]
			    , &c__1);
#line 299 "MB01UW.f"
/* L80: */
#line 299 "MB01UW.f"
		}

#line 301 "MB01UW.f"
	    } else {
#line 302 "MB01UW.f"
		i__1 = *n;
#line 302 "MB01UW.f"
		for (j = 1; j <= i__1; ++j) {

/*                 Compute the contribution of the subdiagonal of H to */
/*                 the j-th column of the product. */

#line 307 "MB01UW.f"
		    i__2 = *m - 1;
#line 307 "MB01UW.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 308 "MB01UW.f"
			dwork[i__] = h__[i__ + 1 + h_dim1] * a[i__ + j * 
				a_dim1];
#line 309 "MB01UW.f"
/* L90: */
#line 309 "MB01UW.f"
		    }

/*                 Multiply the upper triangle of H by the j-th column */
/*                 of A, and add to the above result. */

#line 314 "MB01UW.f"
		    dtrmv_("Upper", trans, "Non-unit", m, &h__[h_offset], ldh,
			     &a[j * a_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)1, 
			    (ftnlen)8);
#line 316 "MB01UW.f"
		    i__2 = *m - 1;
#line 316 "MB01UW.f"
		    daxpy_(&i__2, &c_b33, &dwork[1], &c__1, &a[j * a_dim1 + 2]
			    , &c__1);
#line 317 "MB01UW.f"
/* L100: */
#line 317 "MB01UW.f"
		}
#line 318 "MB01UW.f"
	    }
#line 319 "MB01UW.f"
	    if (*m > 2) {
#line 319 "MB01UW.f"
		i__1 = *m - 2;
#line 319 "MB01UW.f"
		i__2 = *ldh + 1;
#line 319 "MB01UW.f"
		dswap_(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3]
			, &c__1);
#line 319 "MB01UW.f"
	    }

#line 322 "MB01UW.f"
	} else {

/*           Below, row-wise calculations are used for A. */

#line 326 "MB01UW.f"
	    if (*n > 2) {
#line 326 "MB01UW.f"
		i__1 = *n - 2;
#line 326 "MB01UW.f"
		i__2 = *ldh + 1;
#line 326 "MB01UW.f"
		dswap_(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3]
			, &c__1);
#line 326 "MB01UW.f"
	    }
#line 328 "MB01UW.f"
	    if (ltrans) {
#line 329 "MB01UW.f"
		i__1 = *m;
#line 329 "MB01UW.f"
		for (i__ = 1; i__ <= i__1; ++i__) {

/*                 Compute the contribution of the subdiagonal of H to */
/*                 the i-th row of the product. */

#line 334 "MB01UW.f"
		    i__2 = *n - 1;
#line 334 "MB01UW.f"
		    for (j = 1; j <= i__2; ++j) {
#line 335 "MB01UW.f"
			dwork[j] = a[i__ + j * a_dim1] * h__[j + 1 + h_dim1];
#line 336 "MB01UW.f"
/* L110: */
#line 336 "MB01UW.f"
		    }

/*                 Multiply the i-th row of A by the upper triangle of H, */
/*                 and add to the above result. */

#line 341 "MB01UW.f"
		    dtrmv_("Upper", "NoTranspose", "Non-unit", n, &h__[
			    h_offset], ldh, &a[i__ + a_dim1], lda, (ftnlen)5, 
			    (ftnlen)11, (ftnlen)8);
#line 343 "MB01UW.f"
		    i__2 = *n - 1;
#line 343 "MB01UW.f"
		    daxpy_(&i__2, &c_b33, &dwork[1], &c__1, &a[i__ + (a_dim1 
			    << 1)], lda);
#line 344 "MB01UW.f"
/* L120: */
#line 344 "MB01UW.f"
		}

#line 346 "MB01UW.f"
	    } else {
#line 347 "MB01UW.f"
		i__1 = *m;
#line 347 "MB01UW.f"
		for (i__ = 1; i__ <= i__1; ++i__) {

/*                 Compute the contribution of the subdiagonal of H to */
/*                 the i-th row of the product. */

#line 352 "MB01UW.f"
		    i__2 = *n - 1;
#line 352 "MB01UW.f"
		    for (j = 1; j <= i__2; ++j) {
#line 353 "MB01UW.f"
			dwork[j] = a[i__ + (j + 1) * a_dim1] * h__[j + 1 + 
				h_dim1];
#line 354 "MB01UW.f"
/* L130: */
#line 354 "MB01UW.f"
		    }

/*                 Multiply the i-th row of A by the upper triangle of H, */
/*                 and add to the above result. */

#line 359 "MB01UW.f"
		    dtrmv_("Upper", "Transpose", "Non-unit", n, &h__[h_offset]
			    , ldh, &a[i__ + a_dim1], lda, (ftnlen)5, (ftnlen)
			    9, (ftnlen)8);
#line 361 "MB01UW.f"
		    i__2 = *n - 1;
#line 361 "MB01UW.f"
		    daxpy_(&i__2, &c_b33, &dwork[1], &c__1, &a[i__ + a_dim1], 
			    lda);
#line 362 "MB01UW.f"
/* L140: */
#line 362 "MB01UW.f"
		}
#line 363 "MB01UW.f"
	    }
#line 364 "MB01UW.f"
	    if (*n > 2) {
#line 364 "MB01UW.f"
		i__1 = *n - 2;
#line 364 "MB01UW.f"
		i__2 = *ldh + 1;
#line 364 "MB01UW.f"
		dswap_(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3]
			, &c__1);
#line 364 "MB01UW.f"
	    }

#line 367 "MB01UW.f"
	}

/*        Scale the result by alpha. */

#line 371 "MB01UW.f"
	if (*alpha != 1.) {
#line 371 "MB01UW.f"
	    dlascl_("General", &c__0, &c__0, &c_b33, alpha, m, n, &a[a_offset]
		    , lda, info, (ftnlen)7);
#line 371 "MB01UW.f"
	}
#line 374 "MB01UW.f"
    }
#line 375 "MB01UW.f"
    return 0;
/* *** Last line of MB01UW *** */
} /* mb01uw_ */

