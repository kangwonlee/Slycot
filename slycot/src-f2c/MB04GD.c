#line 1 "MB04GD.f"
/* MB04GD.f -- translated by f2c (version 20100827).
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

#line 1 "MB04GD.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mb04gd_(integer *m, integer *n, doublereal *a, integer *
	lda, integer *jpvt, doublereal *tau, doublereal *dwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, ma;
    static doublereal aii;
    static integer mki, nki, pvt;
    static doublereal temp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal temp2;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    static integer nfree, itemp;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dgerq2_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *), 
	    dormr2_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen), dlarfg_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);
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

/*     To compute an RQ factorization with row pivoting of a */
/*     real m-by-n matrix A: P*A = R*Q. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the m-by-n matrix A. */
/*             On exit, */
/*             if m <= n, the upper triangle of the subarray */
/*             A(1:m,n-m+1:n) contains the m-by-m upper triangular */
/*             matrix R; */
/*             if m >= n, the elements on and above the (m-n)-th */
/*             subdiagonal contain the m-by-n upper trapezoidal matrix R; */
/*             the remaining elements, with the array TAU, represent the */
/*             orthogonal matrix Q as a product of min(m,n) elementary */
/*             reflectors (see METHOD). */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. LDA >= max(1,M). */

/*     JPVT    (input/output) INTEGER array, dimension (M) */
/*             On entry, if JPVT(i) .ne. 0, the i-th row of A is permuted */
/*             to the bottom of P*A (a trailing row); if JPVT(i) = 0, */
/*             the i-th row of A is a free row. */
/*             On exit, if JPVT(i) = k, then the i-th row of P*A */
/*             was the k-th row of A. */

/*     TAU     (output) DOUBLE PRECISION array, dimension (min(M,N)) */
/*             The scalar factors of the elementary reflectors. */

/*     Workspace */

/*     DWORK    DOUBLE PRECISION array, dimension (3*M) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The matrix Q is represented as a product of elementary reflectors */

/*        Q = H(1) H(2) . . . H(k), where k = min(m,n). */

/*     Each H(i) has the form */

/*        H = I - tau * v * v' */

/*     where tau is a real scalar, and v is a real vector with */
/*     v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit */
/*     in A(m-k+i,1:n-k+i-1), and tau in TAU(i). */

/*     The matrix P is represented in jpvt as follows: If */
/*        jpvt(j) = i */
/*     then the jth row of P is the ith canonical unit vector. */

/*     REFERENCES */

/*     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J., */
/*         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A., */
/*         Ostrouchov, S., and Sorensen, D. */
/*         LAPACK Users' Guide: Second Edition. */
/*         SIAM, Philadelphia, 1995. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Sep. 1997. */
/*     Based on LAPACK Library routines DGEQPF and DGERQ2. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Factorization, matrix algebra, matrix operations, orthogonal */
/*     transformation, triangular form. */

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

#line 148 "MB04GD.f"
    /* Parameter adjustments */
#line 148 "MB04GD.f"
    a_dim1 = *lda;
#line 148 "MB04GD.f"
    a_offset = 1 + a_dim1;
#line 148 "MB04GD.f"
    a -= a_offset;
#line 148 "MB04GD.f"
    --jpvt;
#line 148 "MB04GD.f"
    --tau;
#line 148 "MB04GD.f"
    --dwork;
#line 148 "MB04GD.f"

#line 148 "MB04GD.f"
    /* Function Body */
#line 148 "MB04GD.f"
    *info = 0;
#line 149 "MB04GD.f"
    if (*m < 0) {
#line 150 "MB04GD.f"
	*info = -1;
#line 151 "MB04GD.f"
    } else if (*n < 0) {
#line 152 "MB04GD.f"
	*info = -2;
#line 153 "MB04GD.f"
    } else if (*lda < max(1,*m)) {
#line 154 "MB04GD.f"
	*info = -4;
#line 155 "MB04GD.f"
    }
#line 156 "MB04GD.f"
    if (*info != 0) {
#line 157 "MB04GD.f"
	i__1 = -(*info);
#line 157 "MB04GD.f"
	xerbla_("MB04GD", &i__1, (ftnlen)6);
#line 158 "MB04GD.f"
	return 0;
#line 159 "MB04GD.f"
    }

#line 161 "MB04GD.f"
    k = min(*m,*n);

/*     Move non-free rows bottom. */

#line 165 "MB04GD.f"
    itemp = *m;
#line 166 "MB04GD.f"
    for (i__ = *m; i__ >= 1; --i__) {
#line 167 "MB04GD.f"
	if (jpvt[i__] != 0) {
#line 168 "MB04GD.f"
	    if (i__ != itemp) {
#line 169 "MB04GD.f"
		dswap_(n, &a[i__ + a_dim1], lda, &a[itemp + a_dim1], lda);
#line 170 "MB04GD.f"
		jpvt[i__] = jpvt[itemp];
#line 171 "MB04GD.f"
		jpvt[itemp] = i__;
#line 172 "MB04GD.f"
	    } else {
#line 173 "MB04GD.f"
		jpvt[i__] = i__;
#line 174 "MB04GD.f"
	    }
#line 175 "MB04GD.f"
	    --itemp;
#line 176 "MB04GD.f"
	} else {
#line 177 "MB04GD.f"
	    jpvt[i__] = i__;
#line 178 "MB04GD.f"
	}
#line 179 "MB04GD.f"
/* L10: */
#line 179 "MB04GD.f"
    }
#line 180 "MB04GD.f"
    nfree = *m - itemp;

/*     Compute the RQ factorization and update remaining rows. */

#line 184 "MB04GD.f"
    if (nfree > 0) {
#line 185 "MB04GD.f"
	ma = min(nfree,*n);
#line 186 "MB04GD.f"
	dgerq2_(&ma, n, &a[*m - ma + 1 + a_dim1], lda, &tau[k - ma + 1], &
		dwork[1], info);
#line 188 "MB04GD.f"
	i__1 = *m - ma;
#line 188 "MB04GD.f"
	dormr2_("Right", "Transpose", &i__1, n, &ma, &a[*m - ma + 1 + a_dim1],
		 lda, &tau[k - ma + 1], &a[a_offset], lda, &dwork[1], info, (
		ftnlen)5, (ftnlen)9);
#line 190 "MB04GD.f"
    }

#line 192 "MB04GD.f"
    if (nfree < k) {

/*        Initialize partial row norms. The first ITEMP elements of */
/*        DWORK store the exact row norms. (Here, ITEMP is the number of */
/*        free rows, which have been permuted to be the first ones.) */

#line 198 "MB04GD.f"
	i__1 = itemp;
#line 198 "MB04GD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 199 "MB04GD.f"
	    i__2 = *n - nfree;
#line 199 "MB04GD.f"
	    dwork[i__] = dnrm2_(&i__2, &a[i__ + a_dim1], lda);
#line 200 "MB04GD.f"
	    dwork[*m + i__] = dwork[i__];
#line 201 "MB04GD.f"
/* L20: */
#line 201 "MB04GD.f"
	}

/*        Compute factorization. */

#line 205 "MB04GD.f"
	for (i__ = k - nfree; i__ >= 1; --i__) {

/*           Determine ith pivot row and swap if necessary. */

#line 209 "MB04GD.f"
	    mki = *m - k + i__;
#line 210 "MB04GD.f"
	    nki = *n - k + i__;
#line 211 "MB04GD.f"
	    pvt = idamax_(&mki, &dwork[1], &c__1);

#line 213 "MB04GD.f"
	    if (pvt != mki) {
#line 214 "MB04GD.f"
		dswap_(n, &a[pvt + a_dim1], lda, &a[mki + a_dim1], lda);
#line 215 "MB04GD.f"
		itemp = jpvt[pvt];
#line 216 "MB04GD.f"
		jpvt[pvt] = jpvt[mki];
#line 217 "MB04GD.f"
		jpvt[mki] = itemp;
#line 218 "MB04GD.f"
		dwork[pvt] = dwork[mki];
#line 219 "MB04GD.f"
		dwork[*m + pvt] = dwork[*m + mki];
#line 220 "MB04GD.f"
	    }

/*           Generate elementary reflector H(i) to annihilate */
/*           A(m-k+i,1:n-k+i-1), k = min(m,n). */

#line 225 "MB04GD.f"
	    dlarfg_(&nki, &a[mki + nki * a_dim1], &a[mki + a_dim1], lda, &tau[
		    i__]);

/*           Apply H(i) to A(1:m-k+i-1,1:n-k+i) from the right. */

#line 230 "MB04GD.f"
	    aii = a[mki + nki * a_dim1];
#line 231 "MB04GD.f"
	    a[mki + nki * a_dim1] = 1.;
#line 232 "MB04GD.f"
	    i__1 = mki - 1;
#line 232 "MB04GD.f"
	    dlarf_("Right", &i__1, &nki, &a[mki + a_dim1], lda, &tau[i__], &a[
		    a_offset], lda, &dwork[(*m << 1) + 1], (ftnlen)5);
#line 234 "MB04GD.f"
	    a[mki + nki * a_dim1] = aii;

/*           Update partial row norms. */

#line 238 "MB04GD.f"
	    i__1 = mki - 1;
#line 238 "MB04GD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 239 "MB04GD.f"
		if (dwork[j] != 0.) {
/* Computing 2nd power */
#line 240 "MB04GD.f"
		    d__2 = (d__1 = a[j + nki * a_dim1], abs(d__1)) / dwork[j];
#line 240 "MB04GD.f"
		    temp = 1. - d__2 * d__2;
#line 241 "MB04GD.f"
		    temp = max(temp,0.);
/* Computing 2nd power */
#line 242 "MB04GD.f"
		    d__1 = dwork[j] / dwork[*m + j];
#line 242 "MB04GD.f"
		    temp2 = temp * .05 * (d__1 * d__1) + 1.;
#line 244 "MB04GD.f"
		    if (temp2 == 1.) {
#line 245 "MB04GD.f"
			i__2 = nki - 1;
#line 245 "MB04GD.f"
			dwork[j] = dnrm2_(&i__2, &a[j + a_dim1], lda);
#line 246 "MB04GD.f"
			dwork[*m + j] = dwork[j];
#line 247 "MB04GD.f"
		    } else {
#line 248 "MB04GD.f"
			dwork[j] *= sqrt(temp);
#line 249 "MB04GD.f"
		    }
#line 250 "MB04GD.f"
		}
#line 251 "MB04GD.f"
/* L30: */
#line 251 "MB04GD.f"
	    }

#line 253 "MB04GD.f"
/* L40: */
#line 253 "MB04GD.f"
	}
#line 254 "MB04GD.f"
    }

#line 256 "MB04GD.f"
    return 0;
/* *** Last line of MB04GD *** */
} /* mb04gd_ */

