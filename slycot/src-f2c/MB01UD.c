#line 1 "MB01UD.f"
/* MB01UD.f -- translated by f2c (version 20100827).
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

#line 1 "MB01UD.f"
/* Table of constant values */

static doublereal c_b9 = 0.;
static integer c__1 = 1;

/* Subroutine */ int mb01ud_(char *side, char *trans, integer *m, integer *n, 
	doublereal *alpha, doublereal *h__, integer *ldh, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, integer *info, ftnlen 
	side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, h_dim1, h_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static logical lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), dlaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen), xerbla_(char *, integer *, ftnlen);
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

/*        B = alpha*op( H ) * A, or B = alpha*A * op( H ), */

/*     where alpha is a scalar, A and B are m-by-n matrices, H is an */
/*     upper Hessenberg matrix, and op( H ) is one of */

/*        op( H ) = H   or   op( H ) = H',  the transpose of H. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SIDE    CHARACTER*1 */
/*             Specifies whether the Hessenberg matrix H appears on the */
/*             left or right in the matrix product as follows: */
/*             = 'L':  B = alpha*op( H ) * A; */
/*             = 'R':  B = alpha*A * op( H ). */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of op( H ) to be used in the matrix */
/*             multiplication as follows: */
/*             = 'N':  op( H ) = H; */
/*             = 'T':  op( H ) = H'; */
/*             = 'C':  op( H ) = H'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrices A and B.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrices A and B.  N >= 0. */

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

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading M-by-N part of this array must contain the */
/*             matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             The leading M-by-N part of this array contains the */
/*             computed product. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,M). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The required matrix product is computed in two steps. In the first */
/*     step, the upper triangle of H is used; in the second step, the */
/*     contribution of the subdiagonal is added. A fast BLAS 3 DTRMM */
/*     operation is used in the first step. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, January 1999. */

/*     REVISIONS */

/*     - */

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

#line 144 "MB01UD.f"
    /* Parameter adjustments */
#line 144 "MB01UD.f"
    h_dim1 = *ldh;
#line 144 "MB01UD.f"
    h_offset = 1 + h_dim1;
#line 144 "MB01UD.f"
    h__ -= h_offset;
#line 144 "MB01UD.f"
    a_dim1 = *lda;
#line 144 "MB01UD.f"
    a_offset = 1 + a_dim1;
#line 144 "MB01UD.f"
    a -= a_offset;
#line 144 "MB01UD.f"
    b_dim1 = *ldb;
#line 144 "MB01UD.f"
    b_offset = 1 + b_dim1;
#line 144 "MB01UD.f"
    b -= b_offset;
#line 144 "MB01UD.f"

#line 144 "MB01UD.f"
    /* Function Body */
#line 144 "MB01UD.f"
    *info = 0;
#line 145 "MB01UD.f"
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 146 "MB01UD.f"
    ltrans = lsame_(trans, "T", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1);

#line 148 "MB01UD.f"
    if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 149 "MB01UD.f"
	*info = -1;
#line 150 "MB01UD.f"
    } else if (! ltrans && ! lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 151 "MB01UD.f"
	*info = -2;
#line 152 "MB01UD.f"
    } else if (*m < 0) {
#line 153 "MB01UD.f"
	*info = -3;
#line 154 "MB01UD.f"
    } else if (*n < 0) {
#line 155 "MB01UD.f"
	*info = -4;
#line 156 "MB01UD.f"
    } else if (*ldh < 1 || lside && *ldh < *m || ! lside && *ldh < *n) {
#line 158 "MB01UD.f"
	*info = -7;
#line 159 "MB01UD.f"
    } else if (*lda < max(1,*m)) {
#line 160 "MB01UD.f"
	*info = -9;
#line 161 "MB01UD.f"
    } else if (*ldb < max(1,*m)) {
#line 162 "MB01UD.f"
	*info = -11;
#line 163 "MB01UD.f"
    }

#line 165 "MB01UD.f"
    if (*info != 0) {

/*        Error return. */

#line 169 "MB01UD.f"
	i__1 = -(*info);
#line 169 "MB01UD.f"
	xerbla_("MB01UD", &i__1, (ftnlen)6);
#line 170 "MB01UD.f"
	return 0;
#line 171 "MB01UD.f"
    }

/*     Quick return, if possible. */

#line 175 "MB01UD.f"
    if (min(*m,*n) == 0) {
#line 175 "MB01UD.f"
	return 0;
#line 175 "MB01UD.f"
    }

#line 178 "MB01UD.f"
    if (*alpha == 0.) {

/*        Set B to zero and return. */

#line 182 "MB01UD.f"
	dlaset_("Full", m, n, &c_b9, &c_b9, &b[b_offset], ldb, (ftnlen)4);
#line 183 "MB01UD.f"
	return 0;
#line 184 "MB01UD.f"
    }

/*     Copy A in B and compute one of the matrix products */
/*       B = alpha*op( triu( H ) ) * A, or */
/*       B = alpha*A * op( triu( H ) ), */
/*     involving the upper triangle of H. */

#line 191 "MB01UD.f"
    dlacpy_("Full", m, n, &a[a_offset], lda, &b[b_offset], ldb, (ftnlen)4);
#line 192 "MB01UD.f"
    dtrmm_(side, "Upper", trans, "Non-unit", m, n, alpha, &h__[h_offset], ldh,
	     &b[b_offset], ldb, (ftnlen)1, (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*     Add the contribution of the subdiagonal of H. */
/*     If SIDE = 'L', the subdiagonal of H is swapped with the */
/*     corresponding elements in the first column of H, and the */
/*     calculations are organized for column operations. */

#line 200 "MB01UD.f"
    if (lside) {
#line 201 "MB01UD.f"
	if (*m > 2) {
#line 201 "MB01UD.f"
	    i__1 = *m - 2;
#line 201 "MB01UD.f"
	    i__2 = *ldh + 1;
#line 201 "MB01UD.f"
	    dswap_(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3], &
		    c__1);
#line 201 "MB01UD.f"
	}
#line 203 "MB01UD.f"
	if (ltrans) {
#line 204 "MB01UD.f"
	    i__1 = *n;
#line 204 "MB01UD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 205 "MB01UD.f"
		i__2 = *m - 1;
#line 205 "MB01UD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 206 "MB01UD.f"
		    b[i__ + j * b_dim1] += *alpha * h__[i__ + 1 + h_dim1] * a[
			    i__ + 1 + j * a_dim1];
#line 207 "MB01UD.f"
/* L10: */
#line 207 "MB01UD.f"
		}
#line 208 "MB01UD.f"
/* L20: */
#line 208 "MB01UD.f"
	    }
#line 209 "MB01UD.f"
	} else {
#line 210 "MB01UD.f"
	    i__1 = *n;
#line 210 "MB01UD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 211 "MB01UD.f"
		i__2 = *m;
#line 211 "MB01UD.f"
		for (i__ = 2; i__ <= i__2; ++i__) {
#line 212 "MB01UD.f"
		    b[i__ + j * b_dim1] += *alpha * h__[i__ + h_dim1] * a[i__ 
			    - 1 + j * a_dim1];
#line 213 "MB01UD.f"
/* L30: */
#line 213 "MB01UD.f"
		}
#line 214 "MB01UD.f"
/* L40: */
#line 214 "MB01UD.f"
	    }
#line 215 "MB01UD.f"
	}
#line 216 "MB01UD.f"
	if (*m > 2) {
#line 216 "MB01UD.f"
	    i__1 = *m - 2;
#line 216 "MB01UD.f"
	    i__2 = *ldh + 1;
#line 216 "MB01UD.f"
	    dswap_(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3], &
		    c__1);
#line 216 "MB01UD.f"
	}

#line 219 "MB01UD.f"
    } else {

#line 221 "MB01UD.f"
	if (ltrans) {
#line 222 "MB01UD.f"
	    i__1 = *n - 1;
#line 222 "MB01UD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 223 "MB01UD.f"
		if (h__[j + 1 + j * h_dim1] != 0.) {
#line 223 "MB01UD.f"
		    d__1 = *alpha * h__[j + 1 + j * h_dim1];
#line 223 "MB01UD.f"
		    daxpy_(m, &d__1, &a[j * a_dim1 + 1], &c__1, &b[(j + 1) * 
			    b_dim1 + 1], &c__1);
#line 223 "MB01UD.f"
		}
#line 226 "MB01UD.f"
/* L50: */
#line 226 "MB01UD.f"
	    }
#line 227 "MB01UD.f"
	} else {
#line 228 "MB01UD.f"
	    i__1 = *n - 1;
#line 228 "MB01UD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 229 "MB01UD.f"
		if (h__[j + 1 + j * h_dim1] != 0.) {
#line 229 "MB01UD.f"
		    d__1 = *alpha * h__[j + 1 + j * h_dim1];
#line 229 "MB01UD.f"
		    daxpy_(m, &d__1, &a[(j + 1) * a_dim1 + 1], &c__1, &b[j * 
			    b_dim1 + 1], &c__1);
#line 229 "MB01UD.f"
		}
#line 232 "MB01UD.f"
/* L60: */
#line 232 "MB01UD.f"
	    }
#line 233 "MB01UD.f"
	}
#line 234 "MB01UD.f"
    }

#line 236 "MB01UD.f"
    return 0;
/* *** Last line of MB01UD *** */
} /* mb01ud_ */
