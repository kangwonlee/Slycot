#line 1 "SB03MY.f"
/* SB03MY.f -- translated by f2c (version 20100827).
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

#line 1 "SB03MY.f"
/* Table of constant values */

static integer c__1 = 1;
static logical c_true = TRUE_;
static integer c__2 = 2;
static doublereal c_b25 = 1.;
static doublereal c_b29 = 0.;
static logical c_false = FALSE_;

/* Subroutine */ int sb03my_(char *trana, integer *n, doublereal *a, integer *
	lda, doublereal *c__, integer *ldc, doublereal *scale, integer *info, 
	ftnlen trana_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j, k, l;
    static doublereal x[4]	/* was [2][2] */;
    static integer k1, k2, l1, l2;
    static doublereal a11, db, da11, vec[4]	/* was [2][2] */, dum[1], eps;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ierr;
    static doublereal smin;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sb03mw_(logical *, logical *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer knext, lnext;
    static doublereal xnorm;
    extern /* Subroutine */ int dlaln2_(logical *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *);
    static integer mink1n, mink2n, minl1n, minl2n;
    extern /* Subroutine */ int dlasy2_(logical *, logical *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dlabad_(doublereal *, 
	    doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal scaloc;
    extern doublereal dlanhs_(char *, integer *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static logical notrna, lupper;
    static doublereal smlnum;


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

/*     To solve the real Lyapunov matrix equation */

/*            op(A)'*X + X*op(A) = scale*C */

/*     where op(A) = A or A' (A**T), A is upper quasi-triangular and C is */
/*     symmetric (C = C'). (A' denotes the transpose of the matrix A.) */
/*     A is N-by-N, the right hand side C and the solution X are N-by-N, */
/*     and scale is an output scale factor, set less than or equal to 1 */
/*     to avoid overflow in X. The solution matrix X is overwritten */
/*     onto C. */

/*     A must be in Schur canonical form (as returned by LAPACK routines */
/*     DGEES or DHSEQR), that is, block upper triangular with 1-by-1 and */
/*     2-by-2 diagonal blocks; each 2-by-2 diagonal block has its */
/*     diagonal elements equal and its off-diagonal elements of opposite */
/*     sign. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op(A) to be used, as follows: */
/*             = 'N':  op(A) = A    (No transpose); */
/*             = 'T':  op(A) = A**T (Transpose); */
/*             = 'C':  op(A) = A**T (Conjugate transpose = Transpose). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, X, and C.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             upper quasi-triangular matrix A, in Schur canonical form. */
/*             The part of A below the first sub-diagonal is not */
/*             referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the symmetric matrix C. */
/*             On exit, if INFO >= 0, the leading N-by-N part of this */
/*             array contains the symmetric solution matrix X. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,N). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor, scale, set less than or equal to 1 to */
/*             prevent the solution overflowing. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if A and -A have common or very close eigenvalues; */
/*                   perturbed values were used to solve the equation */
/*                   (but the matrix A is unchanged). */

/*     METHOD */

/*     Bartels-Stewart algorithm is used. A set of equivalent linear */
/*     algebraic systems of equations of order at most four are formed */
/*     and solved using Gaussian elimination with complete pivoting. */

/*     REFERENCES */

/*     [1] Bartels, R.H. and Stewart, G.W.  T */
/*         Solution of the matrix equation A X + XB = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, May 1997. */
/*     Supersedes Release 2.0 routine SB03AY by Control Systems Research */
/*     Group, Kingston Polytechnic, United Kingdom, October 1982. */
/*     Based on DTRLYP by P. Petkov, Tech. University of Sofia, September */
/*     1993. */

/*     REVISIONS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, May 1999. */

/*     KEYWORDS */

/*     Continuous-time system, Lyapunov equation, matrix algebra, real */
/*     Schur form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test input parameters. */

#line 159 "SB03MY.f"
    /* Parameter adjustments */
#line 159 "SB03MY.f"
    a_dim1 = *lda;
#line 159 "SB03MY.f"
    a_offset = 1 + a_dim1;
#line 159 "SB03MY.f"
    a -= a_offset;
#line 159 "SB03MY.f"
    c_dim1 = *ldc;
#line 159 "SB03MY.f"
    c_offset = 1 + c_dim1;
#line 159 "SB03MY.f"
    c__ -= c_offset;
#line 159 "SB03MY.f"

#line 159 "SB03MY.f"
    /* Function Body */
#line 159 "SB03MY.f"
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
#line 160 "SB03MY.f"
    lupper = TRUE_;

#line 162 "SB03MY.f"
    *info = 0;
#line 163 "SB03MY.f"
    if (! notrna && ! lsame_(trana, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trana, "C", (ftnlen)1, (ftnlen)1)) {
#line 165 "SB03MY.f"
	*info = -1;
#line 166 "SB03MY.f"
    } else if (*n < 0) {
#line 167 "SB03MY.f"
	*info = -2;
#line 168 "SB03MY.f"
    } else if (*lda < max(1,*n)) {
#line 169 "SB03MY.f"
	*info = -4;
#line 170 "SB03MY.f"
    } else if (*ldc < max(1,*n)) {
#line 171 "SB03MY.f"
	*info = -6;
#line 172 "SB03MY.f"
    }

#line 174 "SB03MY.f"
    if (*info != 0) {
#line 175 "SB03MY.f"
	i__1 = -(*info);
#line 175 "SB03MY.f"
	xerbla_("SB03MY", &i__1, (ftnlen)6);
#line 176 "SB03MY.f"
	return 0;
#line 177 "SB03MY.f"
    }

#line 179 "SB03MY.f"
    *scale = 1.;

/*     Quick return if possible. */

#line 183 "SB03MY.f"
    if (*n == 0) {
#line 183 "SB03MY.f"
	return 0;
#line 183 "SB03MY.f"
    }

/*     Set constants to control overflow. */

#line 188 "SB03MY.f"
    eps = dlamch_("P", (ftnlen)1);
#line 189 "SB03MY.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 190 "SB03MY.f"
    bignum = 1. / smlnum;
#line 191 "SB03MY.f"
    dlabad_(&smlnum, &bignum);
#line 192 "SB03MY.f"
    smlnum = smlnum * (doublereal) (*n * *n) / eps;
#line 193 "SB03MY.f"
    bignum = 1. / smlnum;

/* Computing MAX */
#line 195 "SB03MY.f"
    d__1 = smlnum, d__2 = eps * dlanhs_("Max", n, &a[a_offset], lda, dum, (
	    ftnlen)3);
#line 195 "SB03MY.f"
    smin = max(d__1,d__2);

#line 197 "SB03MY.f"
    if (notrna) {

/*        Solve    A'*X + X*A = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        upper-left corner column by column by */

/*          A(K,K)'*X(K,L) + X(K,L)*A(L,L) = C(K,L) - R(K,L), */

/*        where */
/*                   K-1                    L-1 */
/*          R(K,L) = SUM [A(I,K)'*X(I,L)] + SUM [X(K,J)*A(J,L)]. */
/*                   I=1                    J=1 */

/*        Start column loop (index = L). */
/*        L1 (L2): column index of the first (last) row of X(K,L). */

#line 214 "SB03MY.f"
	lnext = 1;

#line 216 "SB03MY.f"
	i__1 = *n;
#line 216 "SB03MY.f"
	for (l = 1; l <= i__1; ++l) {
#line 217 "SB03MY.f"
	    if (l < lnext) {
#line 217 "SB03MY.f"
		goto L60;
#line 217 "SB03MY.f"
	    }
#line 219 "SB03MY.f"
	    l1 = l;
#line 220 "SB03MY.f"
	    l2 = l;
#line 221 "SB03MY.f"
	    if (l < *n) {
#line 222 "SB03MY.f"
		if (a[l + 1 + l * a_dim1] != 0.) {
#line 222 "SB03MY.f"
		    ++l2;
#line 222 "SB03MY.f"
		}
#line 224 "SB03MY.f"
		lnext = l2 + 1;
#line 225 "SB03MY.f"
	    }

/*           Start row loop (index = K). */
/*           K1 (K2): row index of the first (last) row of X(K,L). */

#line 230 "SB03MY.f"
	    knext = l;

#line 232 "SB03MY.f"
	    i__2 = *n;
#line 232 "SB03MY.f"
	    for (k = l; k <= i__2; ++k) {
#line 233 "SB03MY.f"
		if (k < knext) {
#line 233 "SB03MY.f"
		    goto L50;
#line 233 "SB03MY.f"
		}
#line 235 "SB03MY.f"
		k1 = k;
#line 236 "SB03MY.f"
		k2 = k;
#line 237 "SB03MY.f"
		if (k < *n) {
#line 238 "SB03MY.f"
		    if (a[k + 1 + k * a_dim1] != 0.) {
#line 238 "SB03MY.f"
			++k2;
#line 238 "SB03MY.f"
		    }
#line 240 "SB03MY.f"
		    knext = k2 + 1;
#line 241 "SB03MY.f"
		}

#line 243 "SB03MY.f"
		if (l1 == l2 && k1 == k2) {
#line 244 "SB03MY.f"
		    i__3 = k1 - 1;
#line 244 "SB03MY.f"
		    i__4 = l1 - 1;
#line 244 "SB03MY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (ddot_(&i__3, &a[k1 * 
			    a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1) 
			    + ddot_(&i__4, &c__[k1 + c_dim1], ldc, &a[l1 * 
			    a_dim1 + 1], &c__1));
#line 247 "SB03MY.f"
		    scaloc = 1.;

#line 249 "SB03MY.f"
		    a11 = a[k1 + k1 * a_dim1] + a[l1 + l1 * a_dim1];
#line 250 "SB03MY.f"
		    da11 = abs(a11);
#line 251 "SB03MY.f"
		    if (da11 <= smin) {
#line 252 "SB03MY.f"
			a11 = smin;
#line 253 "SB03MY.f"
			da11 = smin;
#line 254 "SB03MY.f"
			*info = 1;
#line 255 "SB03MY.f"
		    }
#line 256 "SB03MY.f"
		    db = abs(vec[0]);
#line 257 "SB03MY.f"
		    if (da11 < 1. && db > 1.) {
#line 258 "SB03MY.f"
			if (db > bignum * da11) {
#line 258 "SB03MY.f"
			    scaloc = 1. / db;
#line 258 "SB03MY.f"
			}
#line 260 "SB03MY.f"
		    }
#line 261 "SB03MY.f"
		    x[0] = vec[0] * scaloc / a11;

#line 263 "SB03MY.f"
		    if (scaloc != 1.) {

#line 265 "SB03MY.f"
			i__3 = *n;
#line 265 "SB03MY.f"
			for (j = 1; j <= i__3; ++j) {
#line 266 "SB03MY.f"
			    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 267 "SB03MY.f"
/* L10: */
#line 267 "SB03MY.f"
			}

#line 269 "SB03MY.f"
			*scale *= scaloc;
#line 270 "SB03MY.f"
		    }
#line 271 "SB03MY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 272 "SB03MY.f"
		    if (k1 != l1) {
#line 273 "SB03MY.f"
			c__[l1 + k1 * c_dim1] = x[0];
#line 274 "SB03MY.f"
		    }

#line 276 "SB03MY.f"
		} else if (l1 == l2 && k1 != k2) {

#line 278 "SB03MY.f"
		    i__3 = k1 - 1;
#line 278 "SB03MY.f"
		    i__4 = l1 - 1;
#line 278 "SB03MY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (ddot_(&i__3, &a[k1 * 
			    a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1) 
			    + ddot_(&i__4, &c__[k1 + c_dim1], ldc, &a[l1 * 
			    a_dim1 + 1], &c__1));

#line 282 "SB03MY.f"
		    i__3 = k1 - 1;
#line 282 "SB03MY.f"
		    i__4 = l1 - 1;
#line 282 "SB03MY.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (ddot_(&i__3, &a[k2 * 
			    a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1) 
			    + ddot_(&i__4, &c__[k2 + c_dim1], ldc, &a[l1 * 
			    a_dim1 + 1], &c__1));

#line 286 "SB03MY.f"
		    d__1 = -a[l1 + l1 * a_dim1];
#line 286 "SB03MY.f"
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b25, &a[k1 + k1 *
			     a_dim1], lda, &c_b25, &c_b25, vec, &c__2, &d__1, 
			    &c_b29, x, &c__2, &scaloc, &xnorm, &ierr);
#line 289 "SB03MY.f"
		    if (ierr != 0) {
#line 289 "SB03MY.f"
			*info = 1;
#line 289 "SB03MY.f"
		    }

#line 292 "SB03MY.f"
		    if (scaloc != 1.) {

#line 294 "SB03MY.f"
			i__3 = *n;
#line 294 "SB03MY.f"
			for (j = 1; j <= i__3; ++j) {
#line 295 "SB03MY.f"
			    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 296 "SB03MY.f"
/* L20: */
#line 296 "SB03MY.f"
			}

#line 298 "SB03MY.f"
			*scale *= scaloc;
#line 299 "SB03MY.f"
		    }
#line 300 "SB03MY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 301 "SB03MY.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 302 "SB03MY.f"
		    c__[l1 + k1 * c_dim1] = x[0];
#line 303 "SB03MY.f"
		    c__[l1 + k2 * c_dim1] = x[1];

#line 305 "SB03MY.f"
		} else if (l1 != l2 && k1 == k2) {

#line 307 "SB03MY.f"
		    i__3 = k1 - 1;
#line 307 "SB03MY.f"
		    i__4 = l1 - 1;
#line 307 "SB03MY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (ddot_(&i__3, &a[k1 * 
			    a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1) 
			    + ddot_(&i__4, &c__[k1 + c_dim1], ldc, &a[l1 * 
			    a_dim1 + 1], &c__1));

#line 311 "SB03MY.f"
		    i__3 = k1 - 1;
#line 311 "SB03MY.f"
		    i__4 = l1 - 1;
#line 311 "SB03MY.f"
		    vec[1] = c__[k1 + l2 * c_dim1] - (ddot_(&i__3, &a[k1 * 
			    a_dim1 + 1], &c__1, &c__[l2 * c_dim1 + 1], &c__1) 
			    + ddot_(&i__4, &c__[k1 + c_dim1], ldc, &a[l2 * 
			    a_dim1 + 1], &c__1));

#line 315 "SB03MY.f"
		    d__1 = -a[k1 + k1 * a_dim1];
#line 315 "SB03MY.f"
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b25, &a[l1 + l1 *
			     a_dim1], lda, &c_b25, &c_b25, vec, &c__2, &d__1, 
			    &c_b29, x, &c__2, &scaloc, &xnorm, &ierr);
#line 318 "SB03MY.f"
		    if (ierr != 0) {
#line 318 "SB03MY.f"
			*info = 1;
#line 318 "SB03MY.f"
		    }

#line 321 "SB03MY.f"
		    if (scaloc != 1.) {

#line 323 "SB03MY.f"
			i__3 = *n;
#line 323 "SB03MY.f"
			for (j = 1; j <= i__3; ++j) {
#line 324 "SB03MY.f"
			    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 325 "SB03MY.f"
/* L30: */
#line 325 "SB03MY.f"
			}

#line 327 "SB03MY.f"
			*scale *= scaloc;
#line 328 "SB03MY.f"
		    }
#line 329 "SB03MY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 330 "SB03MY.f"
		    c__[k1 + l2 * c_dim1] = x[1];
#line 331 "SB03MY.f"
		    c__[l1 + k1 * c_dim1] = x[0];
#line 332 "SB03MY.f"
		    c__[l2 + k1 * c_dim1] = x[1];

#line 334 "SB03MY.f"
		} else if (l1 != l2 && k1 != k2) {

#line 336 "SB03MY.f"
		    i__3 = k1 - 1;
#line 336 "SB03MY.f"
		    i__4 = l1 - 1;
#line 336 "SB03MY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (ddot_(&i__3, &a[k1 * 
			    a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1) 
			    + ddot_(&i__4, &c__[k1 + c_dim1], ldc, &a[l1 * 
			    a_dim1 + 1], &c__1));

#line 340 "SB03MY.f"
		    i__3 = k1 - 1;
#line 340 "SB03MY.f"
		    i__4 = l1 - 1;
#line 340 "SB03MY.f"
		    vec[2] = c__[k1 + l2 * c_dim1] - (ddot_(&i__3, &a[k1 * 
			    a_dim1 + 1], &c__1, &c__[l2 * c_dim1 + 1], &c__1) 
			    + ddot_(&i__4, &c__[k1 + c_dim1], ldc, &a[l2 * 
			    a_dim1 + 1], &c__1));

#line 344 "SB03MY.f"
		    i__3 = k1 - 1;
#line 344 "SB03MY.f"
		    i__4 = l1 - 1;
#line 344 "SB03MY.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (ddot_(&i__3, &a[k2 * 
			    a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1) 
			    + ddot_(&i__4, &c__[k2 + c_dim1], ldc, &a[l1 * 
			    a_dim1 + 1], &c__1));

#line 348 "SB03MY.f"
		    i__3 = k1 - 1;
#line 348 "SB03MY.f"
		    i__4 = l1 - 1;
#line 348 "SB03MY.f"
		    vec[3] = c__[k2 + l2 * c_dim1] - (ddot_(&i__3, &a[k2 * 
			    a_dim1 + 1], &c__1, &c__[l2 * c_dim1 + 1], &c__1) 
			    + ddot_(&i__4, &c__[k2 + c_dim1], ldc, &a[l2 * 
			    a_dim1 + 1], &c__1));

#line 352 "SB03MY.f"
		    if (k1 == l1) {
#line 353 "SB03MY.f"
			sb03mw_(&c_false, &lupper, &a[k1 + k1 * a_dim1], lda, 
				vec, &c__2, &scaloc, x, &c__2, &xnorm, &ierr);
#line 355 "SB03MY.f"
			if (lupper) {
#line 356 "SB03MY.f"
			    x[1] = x[2];
#line 357 "SB03MY.f"
			} else {
#line 358 "SB03MY.f"
			    x[2] = x[1];
#line 359 "SB03MY.f"
			}
#line 360 "SB03MY.f"
		    } else {
#line 361 "SB03MY.f"
			dlasy2_(&c_true, &c_false, &c__1, &c__2, &c__2, &a[k1 
				+ k1 * a_dim1], lda, &a[l1 + l1 * a_dim1], 
				lda, vec, &c__2, &scaloc, x, &c__2, &xnorm, &
				ierr);
#line 364 "SB03MY.f"
		    }
#line 365 "SB03MY.f"
		    if (ierr != 0) {
#line 365 "SB03MY.f"
			*info = 1;
#line 365 "SB03MY.f"
		    }

#line 368 "SB03MY.f"
		    if (scaloc != 1.) {

#line 370 "SB03MY.f"
			i__3 = *n;
#line 370 "SB03MY.f"
			for (j = 1; j <= i__3; ++j) {
#line 371 "SB03MY.f"
			    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 372 "SB03MY.f"
/* L40: */
#line 372 "SB03MY.f"
			}

#line 374 "SB03MY.f"
			*scale *= scaloc;
#line 375 "SB03MY.f"
		    }
#line 376 "SB03MY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 377 "SB03MY.f"
		    c__[k1 + l2 * c_dim1] = x[2];
#line 378 "SB03MY.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 379 "SB03MY.f"
		    c__[k2 + l2 * c_dim1] = x[3];
#line 380 "SB03MY.f"
		    if (k1 != l1) {
#line 381 "SB03MY.f"
			c__[l1 + k1 * c_dim1] = x[0];
#line 382 "SB03MY.f"
			c__[l2 + k1 * c_dim1] = x[2];
#line 383 "SB03MY.f"
			c__[l1 + k2 * c_dim1] = x[1];
#line 384 "SB03MY.f"
			c__[l2 + k2 * c_dim1] = x[3];
#line 385 "SB03MY.f"
		    }
#line 386 "SB03MY.f"
		}

#line 388 "SB03MY.f"
L50:
#line 388 "SB03MY.f"
		;
#line 388 "SB03MY.f"
	    }

#line 390 "SB03MY.f"
L60:
#line 390 "SB03MY.f"
	    ;
#line 390 "SB03MY.f"
	}

#line 392 "SB03MY.f"
    } else {

/*        Solve    A*X + X*A' = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        bottom-right corner column by column by */

/*            A(K,K)*X(K,L) + X(K,L)*A(L,L)' = C(K,L) - R(K,L), */

/*        where */
/*                      N                     N */
/*            R(K,L) = SUM [A(K,I)*X(I,L)] + SUM [X(K,J)*A(L,J)']. */
/*                    I=K+1                 J=L+1 */

/*        Start column loop (index = L). */
/*        L1 (L2): column index of the first (last) row of X(K,L). */

#line 409 "SB03MY.f"
	lnext = *n;

#line 411 "SB03MY.f"
	for (l = *n; l >= 1; --l) {
#line 412 "SB03MY.f"
	    if (l > lnext) {
#line 412 "SB03MY.f"
		goto L120;
#line 412 "SB03MY.f"
	    }
#line 414 "SB03MY.f"
	    l1 = l;
#line 415 "SB03MY.f"
	    l2 = l;
#line 416 "SB03MY.f"
	    if (l > 1) {
#line 417 "SB03MY.f"
		if (a[l + (l - 1) * a_dim1] != 0.) {
#line 417 "SB03MY.f"
		    --l1;
#line 417 "SB03MY.f"
		}
#line 419 "SB03MY.f"
		lnext = l1 - 1;
#line 420 "SB03MY.f"
	    }
/* Computing MIN */
#line 421 "SB03MY.f"
	    i__1 = l1 + 1;
#line 421 "SB03MY.f"
	    minl1n = min(i__1,*n);
/* Computing MIN */
#line 422 "SB03MY.f"
	    i__1 = l2 + 1;
#line 422 "SB03MY.f"
	    minl2n = min(i__1,*n);

/*           Start row loop (index = K). */
/*           K1 (K2): row index of the first (last) row of X(K,L). */

#line 427 "SB03MY.f"
	    knext = l;

#line 429 "SB03MY.f"
	    for (k = l; k >= 1; --k) {
#line 430 "SB03MY.f"
		if (k > knext) {
#line 430 "SB03MY.f"
		    goto L110;
#line 430 "SB03MY.f"
		}
#line 432 "SB03MY.f"
		k1 = k;
#line 433 "SB03MY.f"
		k2 = k;
#line 434 "SB03MY.f"
		if (k > 1) {
#line 435 "SB03MY.f"
		    if (a[k + (k - 1) * a_dim1] != 0.) {
#line 435 "SB03MY.f"
			--k1;
#line 435 "SB03MY.f"
		    }
#line 437 "SB03MY.f"
		    knext = k1 - 1;
#line 438 "SB03MY.f"
		}
/* Computing MIN */
#line 439 "SB03MY.f"
		i__1 = k1 + 1;
#line 439 "SB03MY.f"
		mink1n = min(i__1,*n);
/* Computing MIN */
#line 440 "SB03MY.f"
		i__1 = k2 + 1;
#line 440 "SB03MY.f"
		mink2n = min(i__1,*n);

#line 442 "SB03MY.f"
		if (l1 == l2 && k1 == k2) {
#line 443 "SB03MY.f"
		    i__1 = *n - k1;
#line 443 "SB03MY.f"
		    i__2 = *n - l1;
#line 443 "SB03MY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (ddot_(&i__1, &a[k1 + 
			    mink1n * a_dim1], lda, &c__[mink1n + l1 * c_dim1],
			     &c__1) + ddot_(&i__2, &c__[k1 + minl1n * c_dim1],
			     ldc, &a[l1 + minl1n * a_dim1], lda));
#line 448 "SB03MY.f"
		    scaloc = 1.;

#line 450 "SB03MY.f"
		    a11 = a[k1 + k1 * a_dim1] + a[l1 + l1 * a_dim1];
#line 451 "SB03MY.f"
		    da11 = abs(a11);
#line 452 "SB03MY.f"
		    if (da11 <= smin) {
#line 453 "SB03MY.f"
			a11 = smin;
#line 454 "SB03MY.f"
			da11 = smin;
#line 455 "SB03MY.f"
			*info = 1;
#line 456 "SB03MY.f"
		    }
#line 457 "SB03MY.f"
		    db = abs(vec[0]);
#line 458 "SB03MY.f"
		    if (da11 < 1. && db > 1.) {
#line 459 "SB03MY.f"
			if (db > bignum * da11) {
#line 459 "SB03MY.f"
			    scaloc = 1. / db;
#line 459 "SB03MY.f"
			}
#line 461 "SB03MY.f"
		    }
#line 462 "SB03MY.f"
		    x[0] = vec[0] * scaloc / a11;

#line 464 "SB03MY.f"
		    if (scaloc != 1.) {

#line 466 "SB03MY.f"
			i__1 = *n;
#line 466 "SB03MY.f"
			for (j = 1; j <= i__1; ++j) {
#line 467 "SB03MY.f"
			    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 468 "SB03MY.f"
/* L70: */
#line 468 "SB03MY.f"
			}

#line 470 "SB03MY.f"
			*scale *= scaloc;
#line 471 "SB03MY.f"
		    }
#line 472 "SB03MY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 473 "SB03MY.f"
		    if (k1 != l1) {
#line 474 "SB03MY.f"
			c__[l1 + k1 * c_dim1] = x[0];
#line 475 "SB03MY.f"
		    }

#line 477 "SB03MY.f"
		} else if (l1 == l2 && k1 != k2) {

#line 479 "SB03MY.f"
		    i__1 = *n - k2;
#line 479 "SB03MY.f"
		    i__2 = *n - l2;
#line 479 "SB03MY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (ddot_(&i__1, &a[k1 + 
			    mink2n * a_dim1], lda, &c__[mink2n + l1 * c_dim1],
			     &c__1) + ddot_(&i__2, &c__[k1 + minl2n * c_dim1],
			     ldc, &a[l1 + minl2n * a_dim1], lda));

#line 485 "SB03MY.f"
		    i__1 = *n - k2;
#line 485 "SB03MY.f"
		    i__2 = *n - l2;
#line 485 "SB03MY.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (ddot_(&i__1, &a[k2 + 
			    mink2n * a_dim1], lda, &c__[mink2n + l1 * c_dim1],
			     &c__1) + ddot_(&i__2, &c__[k2 + minl2n * c_dim1],
			     ldc, &a[l1 + minl2n * a_dim1], lda));

#line 491 "SB03MY.f"
		    d__1 = -a[l1 + l1 * a_dim1];
#line 491 "SB03MY.f"
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b25, &a[k1 + k1 
			    * a_dim1], lda, &c_b25, &c_b25, vec, &c__2, &d__1,
			     &c_b29, x, &c__2, &scaloc, &xnorm, &ierr);
#line 494 "SB03MY.f"
		    if (ierr != 0) {
#line 494 "SB03MY.f"
			*info = 1;
#line 494 "SB03MY.f"
		    }

#line 497 "SB03MY.f"
		    if (scaloc != 1.) {

#line 499 "SB03MY.f"
			i__1 = *n;
#line 499 "SB03MY.f"
			for (j = 1; j <= i__1; ++j) {
#line 500 "SB03MY.f"
			    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 501 "SB03MY.f"
/* L80: */
#line 501 "SB03MY.f"
			}

#line 503 "SB03MY.f"
			*scale *= scaloc;
#line 504 "SB03MY.f"
		    }
#line 505 "SB03MY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 506 "SB03MY.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 507 "SB03MY.f"
		    c__[l1 + k1 * c_dim1] = x[0];
#line 508 "SB03MY.f"
		    c__[l1 + k2 * c_dim1] = x[1];

#line 510 "SB03MY.f"
		} else if (l1 != l2 && k1 == k2) {

#line 512 "SB03MY.f"
		    i__1 = *n - k1;
#line 512 "SB03MY.f"
		    i__2 = *n - l2;
#line 512 "SB03MY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (ddot_(&i__1, &a[k1 + 
			    mink1n * a_dim1], lda, &c__[mink1n + l1 * c_dim1],
			     &c__1) + ddot_(&i__2, &c__[k1 + minl2n * c_dim1],
			     ldc, &a[l1 + minl2n * a_dim1], lda));

#line 518 "SB03MY.f"
		    i__1 = *n - k1;
#line 518 "SB03MY.f"
		    i__2 = *n - l2;
#line 518 "SB03MY.f"
		    vec[1] = c__[k1 + l2 * c_dim1] - (ddot_(&i__1, &a[k1 + 
			    mink1n * a_dim1], lda, &c__[mink1n + l2 * c_dim1],
			     &c__1) + ddot_(&i__2, &c__[k1 + minl2n * c_dim1],
			     ldc, &a[l2 + minl2n * a_dim1], lda));

#line 524 "SB03MY.f"
		    d__1 = -a[k1 + k1 * a_dim1];
#line 524 "SB03MY.f"
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b25, &a[l1 + l1 
			    * a_dim1], lda, &c_b25, &c_b25, vec, &c__2, &d__1,
			     &c_b29, x, &c__2, &scaloc, &xnorm, &ierr);
#line 527 "SB03MY.f"
		    if (ierr != 0) {
#line 527 "SB03MY.f"
			*info = 1;
#line 527 "SB03MY.f"
		    }

#line 530 "SB03MY.f"
		    if (scaloc != 1.) {

#line 532 "SB03MY.f"
			i__1 = *n;
#line 532 "SB03MY.f"
			for (j = 1; j <= i__1; ++j) {
#line 533 "SB03MY.f"
			    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 534 "SB03MY.f"
/* L90: */
#line 534 "SB03MY.f"
			}

#line 536 "SB03MY.f"
			*scale *= scaloc;
#line 537 "SB03MY.f"
		    }
#line 538 "SB03MY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 539 "SB03MY.f"
		    c__[k1 + l2 * c_dim1] = x[1];
#line 540 "SB03MY.f"
		    c__[l1 + k1 * c_dim1] = x[0];
#line 541 "SB03MY.f"
		    c__[l2 + k1 * c_dim1] = x[1];

#line 543 "SB03MY.f"
		} else if (l1 != l2 && k1 != k2) {

#line 545 "SB03MY.f"
		    i__1 = *n - k2;
#line 545 "SB03MY.f"
		    i__2 = *n - l2;
#line 545 "SB03MY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (ddot_(&i__1, &a[k1 + 
			    mink2n * a_dim1], lda, &c__[mink2n + l1 * c_dim1],
			     &c__1) + ddot_(&i__2, &c__[k1 + minl2n * c_dim1],
			     ldc, &a[l1 + minl2n * a_dim1], lda));

#line 551 "SB03MY.f"
		    i__1 = *n - k2;
#line 551 "SB03MY.f"
		    i__2 = *n - l2;
#line 551 "SB03MY.f"
		    vec[2] = c__[k1 + l2 * c_dim1] - (ddot_(&i__1, &a[k1 + 
			    mink2n * a_dim1], lda, &c__[mink2n + l2 * c_dim1],
			     &c__1) + ddot_(&i__2, &c__[k1 + minl2n * c_dim1],
			     ldc, &a[l2 + minl2n * a_dim1], lda));

#line 557 "SB03MY.f"
		    i__1 = *n - k2;
#line 557 "SB03MY.f"
		    i__2 = *n - l2;
#line 557 "SB03MY.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (ddot_(&i__1, &a[k2 + 
			    mink2n * a_dim1], lda, &c__[mink2n + l1 * c_dim1],
			     &c__1) + ddot_(&i__2, &c__[k2 + minl2n * c_dim1],
			     ldc, &a[l1 + minl2n * a_dim1], lda));

#line 563 "SB03MY.f"
		    i__1 = *n - k2;
#line 563 "SB03MY.f"
		    i__2 = *n - l2;
#line 563 "SB03MY.f"
		    vec[3] = c__[k2 + l2 * c_dim1] - (ddot_(&i__1, &a[k2 + 
			    mink2n * a_dim1], lda, &c__[mink2n + l2 * c_dim1],
			     &c__1) + ddot_(&i__2, &c__[k2 + minl2n * c_dim1],
			     ldc, &a[l2 + minl2n * a_dim1], lda));

#line 569 "SB03MY.f"
		    if (k1 == l1) {
#line 570 "SB03MY.f"
			sb03mw_(&c_true, &lupper, &a[k1 + k1 * a_dim1], lda, 
				vec, &c__2, &scaloc, x, &c__2, &xnorm, &ierr);
#line 572 "SB03MY.f"
			if (lupper) {
#line 573 "SB03MY.f"
			    x[1] = x[2];
#line 574 "SB03MY.f"
			} else {
#line 575 "SB03MY.f"
			    x[2] = x[1];
#line 576 "SB03MY.f"
			}
#line 577 "SB03MY.f"
		    } else {
#line 578 "SB03MY.f"
			dlasy2_(&c_false, &c_true, &c__1, &c__2, &c__2, &a[k1 
				+ k1 * a_dim1], lda, &a[l1 + l1 * a_dim1], 
				lda, vec, &c__2, &scaloc, x, &c__2, &xnorm, &
				ierr);
#line 581 "SB03MY.f"
		    }
#line 582 "SB03MY.f"
		    if (ierr != 0) {
#line 582 "SB03MY.f"
			*info = 1;
#line 582 "SB03MY.f"
		    }

#line 585 "SB03MY.f"
		    if (scaloc != 1.) {

#line 587 "SB03MY.f"
			i__1 = *n;
#line 587 "SB03MY.f"
			for (j = 1; j <= i__1; ++j) {
#line 588 "SB03MY.f"
			    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 589 "SB03MY.f"
/* L100: */
#line 589 "SB03MY.f"
			}

#line 591 "SB03MY.f"
			*scale *= scaloc;
#line 592 "SB03MY.f"
		    }
#line 593 "SB03MY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 594 "SB03MY.f"
		    c__[k1 + l2 * c_dim1] = x[2];
#line 595 "SB03MY.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 596 "SB03MY.f"
		    c__[k2 + l2 * c_dim1] = x[3];
#line 597 "SB03MY.f"
		    if (k1 != l1) {
#line 598 "SB03MY.f"
			c__[l1 + k1 * c_dim1] = x[0];
#line 599 "SB03MY.f"
			c__[l2 + k1 * c_dim1] = x[2];
#line 600 "SB03MY.f"
			c__[l1 + k2 * c_dim1] = x[1];
#line 601 "SB03MY.f"
			c__[l2 + k2 * c_dim1] = x[3];
#line 602 "SB03MY.f"
		    }
#line 603 "SB03MY.f"
		}

#line 605 "SB03MY.f"
L110:
#line 605 "SB03MY.f"
		;
#line 605 "SB03MY.f"
	    }

#line 607 "SB03MY.f"
L120:
#line 607 "SB03MY.f"
	    ;
#line 607 "SB03MY.f"
	}

#line 609 "SB03MY.f"
    }

#line 611 "SB03MY.f"
    return 0;
/* *** Last line of SB03MY *** */
} /* sb03my_ */

