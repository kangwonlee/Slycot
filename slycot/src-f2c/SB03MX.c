#line 1 "SB03MX.f"
/* SB03MX.f -- translated by f2c (version 20100827).
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

#line 1 "SB03MX.f"
/* Table of constant values */

static doublereal c_b11 = 1.;
static integer c__1 = 1;
static doublereal c_b13 = 0.;
static logical c_true = TRUE_;
static integer c__2 = 2;
static logical c_false = FALSE_;
static integer c_n1 = -1;

/* Subroutine */ int sb03mx_(char *trana, integer *n, doublereal *a, integer *
	lda, doublereal *c__, integer *ldc, doublereal *scale, doublereal *
	dwork, integer *info, ftnlen trana_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j, k, l;
    static doublereal x[4]	/* was [2][2] */;
    static integer k1, k2, l1, l2;
    static doublereal a11, db, p11, p12, p21, p22;
    static integer np1;
    static doublereal da11, vec[4]	/* was [2][2] */, eps;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ierr;
    static doublereal smin;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sb03mv_(logical *, logical *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), sb04px_(logical *, logical *,
	     integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer knext, lnext;
    extern /* Subroutine */ int dsymv_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    static doublereal xnorm;
    extern /* Subroutine */ int dlaln2_(logical *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *);
    static integer mink1n, mink2n, minl1n, minl2n;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
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

/*     To solve the real discrete Lyapunov matrix equation */

/*            op(A)'*X*op(A) - X = scale*C */

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

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if A has almost reciprocal eigenvalues; perturbed */
/*                   values were used to solve the equation (but the */
/*                   matrix A is unchanged). */

/*     METHOD */

/*     A discrete-time version of the Bartels-Stewart algorithm is used. */
/*     A set of equivalent linear algebraic systems of equations of order */
/*     at most four are formed and solved using Gaussian elimination with */
/*     complete pivoting. */

/*     REFERENCES */

/*     [1] Barraud, A.Y.                   T */
/*         A numerical algorithm to solve A XA - X = Q. */
/*         IEEE Trans. Auto. Contr., AC-22, pp. 883-885, 1977. */

/*     [2] Bartels, R.H. and Stewart, G.W.  T */
/*         Solution of the matrix equation A X + XB = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, May 1997. */
/*     Supersedes Release 2.0 routine SB03AZ by Control Systems Research */
/*     Group, Kingston Polytechnic, United Kingdom, October 1982. */
/*     Based on DTRLPD by P. Petkov, Tech. University of Sofia, September */
/*     1993. */

/*     REVISIONS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, May 1999. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000. */
/*     A. Varga, DLR Oberpfaffenhofen, March 2002. */

/*     KEYWORDS */

/*     Discrete-time system, Lyapunov equation, matrix algebra, real */
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

#line 171 "SB03MX.f"
    /* Parameter adjustments */
#line 171 "SB03MX.f"
    a_dim1 = *lda;
#line 171 "SB03MX.f"
    a_offset = 1 + a_dim1;
#line 171 "SB03MX.f"
    a -= a_offset;
#line 171 "SB03MX.f"
    c_dim1 = *ldc;
#line 171 "SB03MX.f"
    c_offset = 1 + c_dim1;
#line 171 "SB03MX.f"
    c__ -= c_offset;
#line 171 "SB03MX.f"
    --dwork;
#line 171 "SB03MX.f"

#line 171 "SB03MX.f"
    /* Function Body */
#line 171 "SB03MX.f"
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
#line 172 "SB03MX.f"
    lupper = TRUE_;

#line 174 "SB03MX.f"
    *info = 0;
#line 175 "SB03MX.f"
    if (! notrna && ! lsame_(trana, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trana, "C", (ftnlen)1, (ftnlen)1)) {
#line 177 "SB03MX.f"
	*info = -1;
#line 178 "SB03MX.f"
    } else if (*n < 0) {
#line 179 "SB03MX.f"
	*info = -2;
#line 180 "SB03MX.f"
    } else if (*lda < max(1,*n)) {
#line 181 "SB03MX.f"
	*info = -4;
#line 182 "SB03MX.f"
    } else if (*ldc < max(1,*n)) {
#line 183 "SB03MX.f"
	*info = -6;
#line 184 "SB03MX.f"
    }

#line 186 "SB03MX.f"
    if (*info != 0) {
#line 187 "SB03MX.f"
	i__1 = -(*info);
#line 187 "SB03MX.f"
	xerbla_("SB03MX", &i__1, (ftnlen)6);
#line 188 "SB03MX.f"
	return 0;
#line 189 "SB03MX.f"
    }

#line 191 "SB03MX.f"
    *scale = 1.;

/*     Quick return if possible. */

#line 195 "SB03MX.f"
    if (*n == 0) {
#line 195 "SB03MX.f"
	return 0;
#line 195 "SB03MX.f"
    }

/*     Set constants to control overflow. */

#line 200 "SB03MX.f"
    eps = dlamch_("P", (ftnlen)1);
#line 201 "SB03MX.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 202 "SB03MX.f"
    bignum = 1. / smlnum;
#line 203 "SB03MX.f"
    dlabad_(&smlnum, &bignum);
#line 204 "SB03MX.f"
    smlnum = smlnum * (doublereal) (*n * *n) / eps;
#line 205 "SB03MX.f"
    bignum = 1. / smlnum;

/* Computing MAX */
#line 207 "SB03MX.f"
    d__1 = smlnum, d__2 = eps * dlanhs_("Max", n, &a[a_offset], lda, &dwork[1]
	    , (ftnlen)3);
#line 207 "SB03MX.f"
    smin = max(d__1,d__2);
#line 208 "SB03MX.f"
    np1 = *n + 1;

#line 210 "SB03MX.f"
    if (notrna) {

/*        Solve    A'*X*A - X = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        upper-left corner column by column by */

/*          A(K,K)'*X(K,L)*A(L,L) - X(K,L) = C(K,L) - R(K,L), */

/*        where */
/*                    K           L-1 */
/*          R(K,L) = SUM {A(I,K)'*SUM [X(I,J)*A(J,L)]} + */
/*                   I=1          J=1 */

/*                    K-1 */
/*                   {SUM [A(I,K)'*X(I,L)]}*A(L,L). */
/*                    I=1 */

/*        Start column loop (index = L). */
/*        L1 (L2): column index of the first (last) row of X(K,L). */

#line 231 "SB03MX.f"
	lnext = 1;

#line 233 "SB03MX.f"
	i__1 = *n;
#line 233 "SB03MX.f"
	for (l = 1; l <= i__1; ++l) {
#line 234 "SB03MX.f"
	    if (l < lnext) {
#line 234 "SB03MX.f"
		goto L60;
#line 234 "SB03MX.f"
	    }
#line 236 "SB03MX.f"
	    l1 = l;
#line 237 "SB03MX.f"
	    l2 = l;
#line 238 "SB03MX.f"
	    if (l < *n) {
#line 239 "SB03MX.f"
		if (a[l + 1 + l * a_dim1] != 0.) {
#line 239 "SB03MX.f"
		    ++l2;
#line 239 "SB03MX.f"
		}
#line 241 "SB03MX.f"
		lnext = l2 + 1;
#line 242 "SB03MX.f"
	    }

/*           Start row loop (index = K). */
/*           K1 (K2): row index of the first (last) row of X(K,L). */

#line 247 "SB03MX.f"
	    dwork[l1] = 0.;
#line 248 "SB03MX.f"
	    dwork[*n + l1] = 0.;
#line 249 "SB03MX.f"
	    i__2 = l1 - 1;
#line 249 "SB03MX.f"
	    dsymv_("Lower", &i__2, &c_b11, &c__[c_offset], ldc, &a[l1 * 
		    a_dim1 + 1], &c__1, &c_b13, &dwork[1], &c__1, (ftnlen)5);
#line 251 "SB03MX.f"
	    i__2 = l1 - 1;
#line 251 "SB03MX.f"
	    dsymv_("Lower", &i__2, &c_b11, &c__[c_offset], ldc, &a[l2 * 
		    a_dim1 + 1], &c__1, &c_b13, &dwork[np1], &c__1, (ftnlen)5)
		    ;

#line 254 "SB03MX.f"
	    knext = l;

#line 256 "SB03MX.f"
	    i__2 = *n;
#line 256 "SB03MX.f"
	    for (k = l; k <= i__2; ++k) {
#line 257 "SB03MX.f"
		if (k < knext) {
#line 257 "SB03MX.f"
		    goto L50;
#line 257 "SB03MX.f"
		}
#line 259 "SB03MX.f"
		k1 = k;
#line 260 "SB03MX.f"
		k2 = k;
#line 261 "SB03MX.f"
		if (k < *n) {
#line 262 "SB03MX.f"
		    if (a[k + 1 + k * a_dim1] != 0.) {
#line 262 "SB03MX.f"
			++k2;
#line 262 "SB03MX.f"
		    }
#line 264 "SB03MX.f"
		    knext = k2 + 1;
#line 265 "SB03MX.f"
		}

#line 267 "SB03MX.f"
		if (l1 == l2 && k1 == k2) {
#line 268 "SB03MX.f"
		    i__3 = l1 - 1;
#line 268 "SB03MX.f"
		    dwork[k1] = ddot_(&i__3, &c__[k1 + c_dim1], ldc, &a[l1 * 
			    a_dim1 + 1], &c__1);

#line 271 "SB03MX.f"
		    i__3 = k1 - 1;
#line 271 "SB03MX.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (ddot_(&k1, &a[k1 * 
			    a_dim1 + 1], &c__1, &dwork[1], &c__1) + a[l1 + l1 
			    * a_dim1] * ddot_(&i__3, &a[k1 * a_dim1 + 1], &
			    c__1, &c__[l1 * c_dim1 + 1], &c__1));
#line 274 "SB03MX.f"
		    scaloc = 1.;

#line 276 "SB03MX.f"
		    a11 = a[k1 + k1 * a_dim1] * a[l1 + l1 * a_dim1] - 1.;
#line 277 "SB03MX.f"
		    da11 = abs(a11);
#line 278 "SB03MX.f"
		    if (da11 <= smin) {
#line 279 "SB03MX.f"
			a11 = smin;
#line 280 "SB03MX.f"
			da11 = smin;
#line 281 "SB03MX.f"
			*info = 1;
#line 282 "SB03MX.f"
		    }
#line 283 "SB03MX.f"
		    db = abs(vec[0]);
#line 284 "SB03MX.f"
		    if (da11 < 1. && db > 1.) {
#line 285 "SB03MX.f"
			if (db > bignum * da11) {
#line 285 "SB03MX.f"
			    scaloc = 1. / db;
#line 285 "SB03MX.f"
			}
#line 287 "SB03MX.f"
		    }
#line 288 "SB03MX.f"
		    x[0] = vec[0] * scaloc / a11;

#line 290 "SB03MX.f"
		    if (scaloc != 1.) {

#line 292 "SB03MX.f"
			i__3 = *n;
#line 292 "SB03MX.f"
			for (j = 1; j <= i__3; ++j) {
#line 293 "SB03MX.f"
			    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 294 "SB03MX.f"
/* L10: */
#line 294 "SB03MX.f"
			}

#line 296 "SB03MX.f"
			dscal_(n, &scaloc, &dwork[1], &c__1);
#line 297 "SB03MX.f"
			*scale *= scaloc;
#line 298 "SB03MX.f"
		    }
#line 299 "SB03MX.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 300 "SB03MX.f"
		    if (k1 != l1) {
#line 301 "SB03MX.f"
			c__[l1 + k1 * c_dim1] = x[0];
#line 302 "SB03MX.f"
		    }

#line 304 "SB03MX.f"
		} else if (l1 == l2 && k1 != k2) {

#line 306 "SB03MX.f"
		    i__3 = l1 - 1;
#line 306 "SB03MX.f"
		    dwork[k1] = ddot_(&i__3, &c__[k1 + c_dim1], ldc, &a[l1 * 
			    a_dim1 + 1], &c__1);
#line 308 "SB03MX.f"
		    i__3 = l1 - 1;
#line 308 "SB03MX.f"
		    dwork[k2] = ddot_(&i__3, &c__[k2 + c_dim1], ldc, &a[l1 * 
			    a_dim1 + 1], &c__1);

#line 311 "SB03MX.f"
		    i__3 = k1 - 1;
#line 311 "SB03MX.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (ddot_(&k2, &a[k1 * 
			    a_dim1 + 1], &c__1, &dwork[1], &c__1) + a[l1 + l1 
			    * a_dim1] * ddot_(&i__3, &a[k1 * a_dim1 + 1], &
			    c__1, &c__[l1 * c_dim1 + 1], &c__1));

#line 315 "SB03MX.f"
		    i__3 = k1 - 1;
#line 315 "SB03MX.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (ddot_(&k2, &a[k2 * 
			    a_dim1 + 1], &c__1, &dwork[1], &c__1) + a[l1 + l1 
			    * a_dim1] * ddot_(&i__3, &a[k2 * a_dim1 + 1], &
			    c__1, &c__[l1 * c_dim1 + 1], &c__1));

#line 319 "SB03MX.f"
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &a[l1 + l1 * a_dim1]
			    , &a[k1 + k1 * a_dim1], lda, &c_b11, &c_b11, vec, 
			    &c__2, &c_b11, &c_b13, x, &c__2, &scaloc, &xnorm, 
			    &ierr);
#line 322 "SB03MX.f"
		    if (ierr != 0) {
#line 322 "SB03MX.f"
			*info = 1;
#line 322 "SB03MX.f"
		    }

#line 325 "SB03MX.f"
		    if (scaloc != 1.) {

#line 327 "SB03MX.f"
			i__3 = *n;
#line 327 "SB03MX.f"
			for (j = 1; j <= i__3; ++j) {
#line 328 "SB03MX.f"
			    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 329 "SB03MX.f"
/* L20: */
#line 329 "SB03MX.f"
			}

#line 331 "SB03MX.f"
			dscal_(n, &scaloc, &dwork[1], &c__1);
#line 332 "SB03MX.f"
			*scale *= scaloc;
#line 333 "SB03MX.f"
		    }
#line 334 "SB03MX.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 335 "SB03MX.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 336 "SB03MX.f"
		    c__[l1 + k1 * c_dim1] = x[0];
#line 337 "SB03MX.f"
		    c__[l1 + k2 * c_dim1] = x[1];

#line 339 "SB03MX.f"
		} else if (l1 != l2 && k1 == k2) {

#line 341 "SB03MX.f"
		    i__3 = l1 - 1;
#line 341 "SB03MX.f"
		    dwork[k1] = ddot_(&i__3, &c__[k1 + c_dim1], ldc, &a[l1 * 
			    a_dim1 + 1], &c__1);
#line 343 "SB03MX.f"
		    i__3 = l1 - 1;
#line 343 "SB03MX.f"
		    dwork[*n + k1] = ddot_(&i__3, &c__[k1 + c_dim1], ldc, &a[
			    l2 * a_dim1 + 1], &c__1);
#line 345 "SB03MX.f"
		    i__3 = k1 - 1;
#line 345 "SB03MX.f"
		    p11 = ddot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 346 "SB03MX.f"
		    i__3 = k1 - 1;
#line 346 "SB03MX.f"
		    p12 = ddot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);

#line 348 "SB03MX.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (ddot_(&k1, &a[k1 * 
			    a_dim1 + 1], &c__1, &dwork[1], &c__1) + p11 * a[
			    l1 + l1 * a_dim1] + p12 * a[l2 + l1 * a_dim1]);

#line 352 "SB03MX.f"
		    vec[1] = c__[k1 + l2 * c_dim1] - (ddot_(&k1, &a[k1 * 
			    a_dim1 + 1], &c__1, &dwork[np1], &c__1) + p11 * a[
			    l1 + l2 * a_dim1] + p12 * a[l2 + l2 * a_dim1]);

#line 356 "SB03MX.f"
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &a[k1 + k1 * a_dim1]
			    , &a[l1 + l1 * a_dim1], lda, &c_b11, &c_b11, vec, 
			    &c__2, &c_b11, &c_b13, x, &c__2, &scaloc, &xnorm, 
			    &ierr);
#line 359 "SB03MX.f"
		    if (ierr != 0) {
#line 359 "SB03MX.f"
			*info = 1;
#line 359 "SB03MX.f"
		    }

#line 362 "SB03MX.f"
		    if (scaloc != 1.) {

#line 364 "SB03MX.f"
			i__3 = *n;
#line 364 "SB03MX.f"
			for (j = 1; j <= i__3; ++j) {
#line 365 "SB03MX.f"
			    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 366 "SB03MX.f"
/* L30: */
#line 366 "SB03MX.f"
			}

#line 368 "SB03MX.f"
			dscal_(n, &scaloc, &dwork[1], &c__1);
#line 369 "SB03MX.f"
			dscal_(n, &scaloc, &dwork[np1], &c__1);
#line 370 "SB03MX.f"
			*scale *= scaloc;
#line 371 "SB03MX.f"
		    }
#line 372 "SB03MX.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 373 "SB03MX.f"
		    c__[k1 + l2 * c_dim1] = x[1];
#line 374 "SB03MX.f"
		    c__[l1 + k1 * c_dim1] = x[0];
#line 375 "SB03MX.f"
		    c__[l2 + k1 * c_dim1] = x[1];

#line 377 "SB03MX.f"
		} else if (l1 != l2 && k1 != k2) {

#line 379 "SB03MX.f"
		    i__3 = l1 - 1;
#line 379 "SB03MX.f"
		    dwork[k1] = ddot_(&i__3, &c__[k1 + c_dim1], ldc, &a[l1 * 
			    a_dim1 + 1], &c__1);
#line 381 "SB03MX.f"
		    i__3 = l1 - 1;
#line 381 "SB03MX.f"
		    dwork[k2] = ddot_(&i__3, &c__[k2 + c_dim1], ldc, &a[l1 * 
			    a_dim1 + 1], &c__1);
#line 383 "SB03MX.f"
		    i__3 = l1 - 1;
#line 383 "SB03MX.f"
		    dwork[*n + k1] = ddot_(&i__3, &c__[k1 + c_dim1], ldc, &a[
			    l2 * a_dim1 + 1], &c__1);
#line 385 "SB03MX.f"
		    i__3 = l1 - 1;
#line 385 "SB03MX.f"
		    dwork[*n + k2] = ddot_(&i__3, &c__[k2 + c_dim1], ldc, &a[
			    l2 * a_dim1 + 1], &c__1);
#line 387 "SB03MX.f"
		    i__3 = k1 - 1;
#line 387 "SB03MX.f"
		    p11 = ddot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 388 "SB03MX.f"
		    i__3 = k1 - 1;
#line 388 "SB03MX.f"
		    p12 = ddot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
#line 389 "SB03MX.f"
		    i__3 = k1 - 1;
#line 389 "SB03MX.f"
		    p21 = ddot_(&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 390 "SB03MX.f"
		    i__3 = k1 - 1;
#line 390 "SB03MX.f"
		    p22 = ddot_(&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);

#line 392 "SB03MX.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (ddot_(&k2, &a[k1 * 
			    a_dim1 + 1], &c__1, &dwork[1], &c__1) + p11 * a[
			    l1 + l1 * a_dim1] + p12 * a[l2 + l1 * a_dim1]);

#line 396 "SB03MX.f"
		    vec[2] = c__[k1 + l2 * c_dim1] - (ddot_(&k2, &a[k1 * 
			    a_dim1 + 1], &c__1, &dwork[np1], &c__1) + p11 * a[
			    l1 + l2 * a_dim1] + p12 * a[l2 + l2 * a_dim1]);

#line 400 "SB03MX.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (ddot_(&k2, &a[k2 * 
			    a_dim1 + 1], &c__1, &dwork[1], &c__1) + p21 * a[
			    l1 + l1 * a_dim1] + p22 * a[l2 + l1 * a_dim1]);

#line 404 "SB03MX.f"
		    vec[3] = c__[k2 + l2 * c_dim1] - (ddot_(&k2, &a[k2 * 
			    a_dim1 + 1], &c__1, &dwork[np1], &c__1) + p21 * a[
			    l1 + l2 * a_dim1] + p22 * a[l2 + l2 * a_dim1]);

#line 408 "SB03MX.f"
		    if (k1 == l1) {
#line 409 "SB03MX.f"
			sb03mv_(&c_false, &lupper, &a[k1 + k1 * a_dim1], lda, 
				vec, &c__2, &scaloc, x, &c__2, &xnorm, &ierr);
#line 411 "SB03MX.f"
			if (lupper) {
#line 412 "SB03MX.f"
			    x[1] = x[2];
#line 413 "SB03MX.f"
			} else {
#line 414 "SB03MX.f"
			    x[2] = x[1];
#line 415 "SB03MX.f"
			}
#line 416 "SB03MX.f"
		    } else {
#line 417 "SB03MX.f"
			sb04px_(&c_true, &c_false, &c_n1, &c__2, &c__2, &a[k1 
				+ k1 * a_dim1], lda, &a[l1 + l1 * a_dim1], 
				lda, vec, &c__2, &scaloc, x, &c__2, &xnorm, &
				ierr);
#line 420 "SB03MX.f"
		    }
#line 421 "SB03MX.f"
		    if (ierr != 0) {
#line 421 "SB03MX.f"
			*info = 1;
#line 421 "SB03MX.f"
		    }

#line 424 "SB03MX.f"
		    if (scaloc != 1.) {

#line 426 "SB03MX.f"
			i__3 = *n;
#line 426 "SB03MX.f"
			for (j = 1; j <= i__3; ++j) {
#line 427 "SB03MX.f"
			    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 428 "SB03MX.f"
/* L40: */
#line 428 "SB03MX.f"
			}

#line 430 "SB03MX.f"
			dscal_(n, &scaloc, &dwork[1], &c__1);
#line 431 "SB03MX.f"
			dscal_(n, &scaloc, &dwork[np1], &c__1);
#line 432 "SB03MX.f"
			*scale *= scaloc;
#line 433 "SB03MX.f"
		    }
#line 434 "SB03MX.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 435 "SB03MX.f"
		    c__[k1 + l2 * c_dim1] = x[2];
#line 436 "SB03MX.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 437 "SB03MX.f"
		    c__[k2 + l2 * c_dim1] = x[3];
#line 438 "SB03MX.f"
		    if (k1 != l1) {
#line 439 "SB03MX.f"
			c__[l1 + k1 * c_dim1] = x[0];
#line 440 "SB03MX.f"
			c__[l2 + k1 * c_dim1] = x[2];
#line 441 "SB03MX.f"
			c__[l1 + k2 * c_dim1] = x[1];
#line 442 "SB03MX.f"
			c__[l2 + k2 * c_dim1] = x[3];
#line 443 "SB03MX.f"
		    }
#line 444 "SB03MX.f"
		}

#line 446 "SB03MX.f"
L50:
#line 446 "SB03MX.f"
		;
#line 446 "SB03MX.f"
	    }

#line 448 "SB03MX.f"
L60:
#line 448 "SB03MX.f"
	    ;
#line 448 "SB03MX.f"
	}

#line 450 "SB03MX.f"
    } else {

/*        Solve    A*X*A' - X = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        bottom-right corner column by column by */

/*            A(K,K)*X(K,L)*A(L,L)' - X(K,L) = C(K,L) - R(K,L), */

/*        where */

/*                    N            N */
/*          R(K,L) = SUM {A(K,I)* SUM [X(I,J)*A(L,J)']} + */
/*                   I=K         J=L+1 */

/*                      N */
/*                   { SUM [A(K,J)*X(J,L)]}*A(L,L)' */
/*                    J=K+1 */

/*        Start column loop (index = L) */
/*        L1 (L2): column index of the first (last) row of X(K,L) */

#line 472 "SB03MX.f"
	lnext = *n;

#line 474 "SB03MX.f"
	for (l = *n; l >= 1; --l) {
#line 475 "SB03MX.f"
	    if (l > lnext) {
#line 475 "SB03MX.f"
		goto L120;
#line 475 "SB03MX.f"
	    }
#line 477 "SB03MX.f"
	    l1 = l;
#line 478 "SB03MX.f"
	    l2 = l;
#line 479 "SB03MX.f"
	    if (l > 1) {
#line 480 "SB03MX.f"
		if (a[l + (l - 1) * a_dim1] != 0.) {
#line 481 "SB03MX.f"
		    --l1;
#line 482 "SB03MX.f"
		    dwork[l1] = 0.;
#line 483 "SB03MX.f"
		    dwork[*n + l1] = 0.;
#line 484 "SB03MX.f"
		}
#line 485 "SB03MX.f"
		lnext = l1 - 1;
#line 486 "SB03MX.f"
	    }
/* Computing MIN */
#line 487 "SB03MX.f"
	    i__1 = l1 + 1;
#line 487 "SB03MX.f"
	    minl1n = min(i__1,*n);
/* Computing MIN */
#line 488 "SB03MX.f"
	    i__1 = l2 + 1;
#line 488 "SB03MX.f"
	    minl2n = min(i__1,*n);

/*           Start row loop (index = K) */
/*           K1 (K2): row index of the first (last) row of X(K,L) */

#line 493 "SB03MX.f"
	    if (l2 < *n) {
#line 494 "SB03MX.f"
		i__1 = *n - l2;
#line 494 "SB03MX.f"
		dsymv_("Upper", &i__1, &c_b11, &c__[l2 + 1 + (l2 + 1) * 
			c_dim1], ldc, &a[l1 + (l2 + 1) * a_dim1], lda, &c_b13,
			 &dwork[l2 + 1], &c__1, (ftnlen)5);
#line 496 "SB03MX.f"
		i__1 = *n - l2;
#line 496 "SB03MX.f"
		dsymv_("Upper", &i__1, &c_b11, &c__[l2 + 1 + (l2 + 1) * 
			c_dim1], ldc, &a[l2 + (l2 + 1) * a_dim1], lda, &c_b13,
			 &dwork[np1 + l2], &c__1, (ftnlen)5);
#line 498 "SB03MX.f"
	    }

#line 500 "SB03MX.f"
	    knext = l;

#line 502 "SB03MX.f"
	    for (k = l; k >= 1; --k) {
#line 503 "SB03MX.f"
		if (k > knext) {
#line 503 "SB03MX.f"
		    goto L110;
#line 503 "SB03MX.f"
		}
#line 505 "SB03MX.f"
		k1 = k;
#line 506 "SB03MX.f"
		k2 = k;
#line 507 "SB03MX.f"
		if (k > 1) {
#line 508 "SB03MX.f"
		    if (a[k + (k - 1) * a_dim1] != 0.) {
#line 508 "SB03MX.f"
			--k1;
#line 508 "SB03MX.f"
		    }
#line 510 "SB03MX.f"
		    knext = k1 - 1;
#line 511 "SB03MX.f"
		}
/* Computing MIN */
#line 512 "SB03MX.f"
		i__1 = k1 + 1;
#line 512 "SB03MX.f"
		mink1n = min(i__1,*n);
/* Computing MIN */
#line 513 "SB03MX.f"
		i__1 = k2 + 1;
#line 513 "SB03MX.f"
		mink2n = min(i__1,*n);

#line 515 "SB03MX.f"
		if (l1 == l2 && k1 == k2) {
#line 516 "SB03MX.f"
		    i__1 = *n - l1;
#line 516 "SB03MX.f"
		    dwork[k1] = ddot_(&i__1, &c__[k1 + minl1n * c_dim1], ldc, 
			    &a[l1 + minl1n * a_dim1], lda);

#line 519 "SB03MX.f"
		    i__1 = *n - k1 + 1;
#line 519 "SB03MX.f"
		    i__2 = *n - k1;
#line 519 "SB03MX.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (ddot_(&i__1, &a[k1 + k1 
			    * a_dim1], lda, &dwork[k1], &c__1) + ddot_(&i__2, 
			    &a[k1 + mink1n * a_dim1], lda, &c__[mink1n + l1 * 
			    c_dim1], &c__1) * a[l1 + l1 * a_dim1]);
#line 523 "SB03MX.f"
		    scaloc = 1.;

#line 525 "SB03MX.f"
		    a11 = a[k1 + k1 * a_dim1] * a[l1 + l1 * a_dim1] - 1.;
#line 526 "SB03MX.f"
		    da11 = abs(a11);
#line 527 "SB03MX.f"
		    if (da11 <= smin) {
#line 528 "SB03MX.f"
			a11 = smin;
#line 529 "SB03MX.f"
			da11 = smin;
#line 530 "SB03MX.f"
			*info = 1;
#line 531 "SB03MX.f"
		    }
#line 532 "SB03MX.f"
		    db = abs(vec[0]);
#line 533 "SB03MX.f"
		    if (da11 < 1. && db > 1.) {
#line 534 "SB03MX.f"
			if (db > bignum * da11) {
#line 534 "SB03MX.f"
			    scaloc = 1. / db;
#line 534 "SB03MX.f"
			}
#line 536 "SB03MX.f"
		    }
#line 537 "SB03MX.f"
		    x[0] = vec[0] * scaloc / a11;

#line 539 "SB03MX.f"
		    if (scaloc != 1.) {

#line 541 "SB03MX.f"
			i__1 = *n;
#line 541 "SB03MX.f"
			for (j = 1; j <= i__1; ++j) {
#line 542 "SB03MX.f"
			    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 543 "SB03MX.f"
/* L70: */
#line 543 "SB03MX.f"
			}

#line 545 "SB03MX.f"
			dscal_(n, &scaloc, &dwork[1], &c__1);
#line 546 "SB03MX.f"
			*scale *= scaloc;
#line 547 "SB03MX.f"
		    }
#line 548 "SB03MX.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 549 "SB03MX.f"
		    if (k1 != l1) {
#line 550 "SB03MX.f"
			c__[l1 + k1 * c_dim1] = x[0];
#line 551 "SB03MX.f"
		    }

#line 553 "SB03MX.f"
		} else if (l1 == l2 && k1 != k2) {

#line 555 "SB03MX.f"
		    i__1 = *n - l1;
#line 555 "SB03MX.f"
		    dwork[k1] = ddot_(&i__1, &c__[k1 + minl1n * c_dim1], ldc, 
			    &a[l1 + minl1n * a_dim1], lda);
#line 557 "SB03MX.f"
		    i__1 = *n - l1;
#line 557 "SB03MX.f"
		    dwork[k2] = ddot_(&i__1, &c__[k2 + minl1n * c_dim1], ldc, 
			    &a[l1 + minl1n * a_dim1], lda);

#line 560 "SB03MX.f"
		    i__1 = np1 - k1;
#line 560 "SB03MX.f"
		    i__2 = *n - k2;
#line 560 "SB03MX.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (ddot_(&i__1, &a[k1 + k1 
			    * a_dim1], lda, &dwork[k1], &c__1) + ddot_(&i__2, 
			    &a[k1 + mink2n * a_dim1], lda, &c__[mink2n + l1 * 
			    c_dim1], &c__1) * a[l1 + l1 * a_dim1]);

#line 565 "SB03MX.f"
		    i__1 = np1 - k1;
#line 565 "SB03MX.f"
		    i__2 = *n - k2;
#line 565 "SB03MX.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (ddot_(&i__1, &a[k2 + k1 
			    * a_dim1], lda, &dwork[k1], &c__1) + ddot_(&i__2, 
			    &a[k2 + mink2n * a_dim1], lda, &c__[mink2n + l1 * 
			    c_dim1], &c__1) * a[l1 + l1 * a_dim1]);

#line 570 "SB03MX.f"
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &a[l1 + l1 * 
			    a_dim1], &a[k1 + k1 * a_dim1], lda, &c_b11, &
			    c_b11, vec, &c__2, &c_b11, &c_b13, x, &c__2, &
			    scaloc, &xnorm, &ierr);
#line 573 "SB03MX.f"
		    if (ierr != 0) {
#line 573 "SB03MX.f"
			*info = 1;
#line 573 "SB03MX.f"
		    }

#line 576 "SB03MX.f"
		    if (scaloc != 1.) {

#line 578 "SB03MX.f"
			i__1 = *n;
#line 578 "SB03MX.f"
			for (j = 1; j <= i__1; ++j) {
#line 579 "SB03MX.f"
			    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 580 "SB03MX.f"
/* L80: */
#line 580 "SB03MX.f"
			}

#line 582 "SB03MX.f"
			dscal_(n, &scaloc, &dwork[1], &c__1);
#line 583 "SB03MX.f"
			*scale *= scaloc;
#line 584 "SB03MX.f"
		    }
#line 585 "SB03MX.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 586 "SB03MX.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 587 "SB03MX.f"
		    c__[l1 + k1 * c_dim1] = x[0];
#line 588 "SB03MX.f"
		    c__[l1 + k2 * c_dim1] = x[1];

#line 590 "SB03MX.f"
		} else if (l1 != l2 && k1 == k2) {

#line 592 "SB03MX.f"
		    i__1 = *n - l2;
#line 592 "SB03MX.f"
		    dwork[k1] = ddot_(&i__1, &c__[k1 + minl2n * c_dim1], ldc, 
			    &a[l1 + minl2n * a_dim1], lda);
#line 594 "SB03MX.f"
		    i__1 = *n - l2;
#line 594 "SB03MX.f"
		    dwork[*n + k1] = ddot_(&i__1, &c__[k1 + minl2n * c_dim1], 
			    ldc, &a[l2 + minl2n * a_dim1], lda);
#line 596 "SB03MX.f"
		    i__1 = *n - k1;
#line 596 "SB03MX.f"
		    p11 = ddot_(&i__1, &a[k1 + mink1n * a_dim1], lda, &c__[
			    mink1n + l1 * c_dim1], &c__1);
#line 598 "SB03MX.f"
		    i__1 = *n - k1;
#line 598 "SB03MX.f"
		    p12 = ddot_(&i__1, &a[k1 + mink1n * a_dim1], lda, &c__[
			    mink1n + l2 * c_dim1], &c__1);

#line 601 "SB03MX.f"
		    i__1 = np1 - k1;
#line 601 "SB03MX.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (ddot_(&i__1, &a[k1 + k1 
			    * a_dim1], lda, &dwork[k1], &c__1) + p11 * a[l1 + 
			    l1 * a_dim1] + p12 * a[l1 + l2 * a_dim1]);

#line 605 "SB03MX.f"
		    i__1 = np1 - k1;
#line 605 "SB03MX.f"
		    vec[1] = c__[k1 + l2 * c_dim1] - (ddot_(&i__1, &a[k1 + k1 
			    * a_dim1], lda, &dwork[*n + k1], &c__1) + p11 * a[
			    l2 + l1 * a_dim1] + p12 * a[l2 + l2 * a_dim1]);

#line 609 "SB03MX.f"
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &a[k1 + k1 * 
			    a_dim1], &a[l1 + l1 * a_dim1], lda, &c_b11, &
			    c_b11, vec, &c__2, &c_b11, &c_b13, x, &c__2, &
			    scaloc, &xnorm, &ierr);
#line 612 "SB03MX.f"
		    if (ierr != 0) {
#line 612 "SB03MX.f"
			*info = 1;
#line 612 "SB03MX.f"
		    }

#line 615 "SB03MX.f"
		    if (scaloc != 1.) {

#line 617 "SB03MX.f"
			i__1 = *n;
#line 617 "SB03MX.f"
			for (j = 1; j <= i__1; ++j) {
#line 618 "SB03MX.f"
			    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 619 "SB03MX.f"
/* L90: */
#line 619 "SB03MX.f"
			}

#line 621 "SB03MX.f"
			dscal_(n, &scaloc, &dwork[1], &c__1);
#line 622 "SB03MX.f"
			dscal_(n, &scaloc, &dwork[np1], &c__1);
#line 623 "SB03MX.f"
			*scale *= scaloc;
#line 624 "SB03MX.f"
		    }
#line 625 "SB03MX.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 626 "SB03MX.f"
		    c__[k1 + l2 * c_dim1] = x[1];
#line 627 "SB03MX.f"
		    c__[l1 + k1 * c_dim1] = x[0];
#line 628 "SB03MX.f"
		    c__[l2 + k1 * c_dim1] = x[1];

#line 630 "SB03MX.f"
		} else if (l1 != l2 && k1 != k2) {

#line 632 "SB03MX.f"
		    i__1 = *n - l2;
#line 632 "SB03MX.f"
		    dwork[k1] = ddot_(&i__1, &c__[k1 + minl2n * c_dim1], ldc, 
			    &a[l1 + minl2n * a_dim1], lda);
#line 634 "SB03MX.f"
		    i__1 = *n - l2;
#line 634 "SB03MX.f"
		    dwork[k2] = ddot_(&i__1, &c__[k2 + minl2n * c_dim1], ldc, 
			    &a[l1 + minl2n * a_dim1], lda);
#line 636 "SB03MX.f"
		    i__1 = *n - l2;
#line 636 "SB03MX.f"
		    dwork[*n + k1] = ddot_(&i__1, &c__[k1 + minl2n * c_dim1], 
			    ldc, &a[l2 + minl2n * a_dim1], lda);
#line 638 "SB03MX.f"
		    i__1 = *n - l2;
#line 638 "SB03MX.f"
		    dwork[*n + k2] = ddot_(&i__1, &c__[k2 + minl2n * c_dim1], 
			    ldc, &a[l2 + minl2n * a_dim1], lda);
#line 640 "SB03MX.f"
		    i__1 = *n - k2;
#line 640 "SB03MX.f"
		    p11 = ddot_(&i__1, &a[k1 + mink2n * a_dim1], lda, &c__[
			    mink2n + l1 * c_dim1], &c__1);
#line 642 "SB03MX.f"
		    i__1 = *n - k2;
#line 642 "SB03MX.f"
		    p12 = ddot_(&i__1, &a[k1 + mink2n * a_dim1], lda, &c__[
			    mink2n + l2 * c_dim1], &c__1);
#line 644 "SB03MX.f"
		    i__1 = *n - k2;
#line 644 "SB03MX.f"
		    p21 = ddot_(&i__1, &a[k2 + mink2n * a_dim1], lda, &c__[
			    mink2n + l1 * c_dim1], &c__1);
#line 646 "SB03MX.f"
		    i__1 = *n - k2;
#line 646 "SB03MX.f"
		    p22 = ddot_(&i__1, &a[k2 + mink2n * a_dim1], lda, &c__[
			    mink2n + l2 * c_dim1], &c__1);

#line 649 "SB03MX.f"
		    i__1 = np1 - k1;
#line 649 "SB03MX.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (ddot_(&i__1, &a[k1 + k1 
			    * a_dim1], lda, &dwork[k1], &c__1) + p11 * a[l1 + 
			    l1 * a_dim1] + p12 * a[l1 + l2 * a_dim1]);

#line 653 "SB03MX.f"
		    i__1 = np1 - k1;
#line 653 "SB03MX.f"
		    vec[2] = c__[k1 + l2 * c_dim1] - (ddot_(&i__1, &a[k1 + k1 
			    * a_dim1], lda, &dwork[*n + k1], &c__1) + p11 * a[
			    l2 + l1 * a_dim1] + p12 * a[l2 + l2 * a_dim1]);

#line 657 "SB03MX.f"
		    i__1 = np1 - k1;
#line 657 "SB03MX.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (ddot_(&i__1, &a[k2 + k1 
			    * a_dim1], lda, &dwork[k1], &c__1) + p21 * a[l1 + 
			    l1 * a_dim1] + p22 * a[l1 + l2 * a_dim1]);

#line 661 "SB03MX.f"
		    i__1 = np1 - k1;
#line 661 "SB03MX.f"
		    vec[3] = c__[k2 + l2 * c_dim1] - (ddot_(&i__1, &a[k2 + k1 
			    * a_dim1], lda, &dwork[*n + k1], &c__1) + p21 * a[
			    l2 + l1 * a_dim1] + p22 * a[l2 + l2 * a_dim1]);

#line 665 "SB03MX.f"
		    if (k1 == l1) {
#line 666 "SB03MX.f"
			sb03mv_(&c_true, &lupper, &a[k1 + k1 * a_dim1], lda, 
				vec, &c__2, &scaloc, x, &c__2, &xnorm, &ierr);
#line 668 "SB03MX.f"
			if (lupper) {
#line 669 "SB03MX.f"
			    x[1] = x[2];
#line 670 "SB03MX.f"
			} else {
#line 671 "SB03MX.f"
			    x[2] = x[1];
#line 672 "SB03MX.f"
			}
#line 673 "SB03MX.f"
		    } else {
#line 674 "SB03MX.f"
			sb04px_(&c_false, &c_true, &c_n1, &c__2, &c__2, &a[k1 
				+ k1 * a_dim1], lda, &a[l1 + l1 * a_dim1], 
				lda, vec, &c__2, &scaloc, x, &c__2, &xnorm, &
				ierr);
#line 677 "SB03MX.f"
		    }
#line 678 "SB03MX.f"
		    if (ierr != 0) {
#line 678 "SB03MX.f"
			*info = 1;
#line 678 "SB03MX.f"
		    }

#line 681 "SB03MX.f"
		    if (scaloc != 1.) {

#line 683 "SB03MX.f"
			i__1 = *n;
#line 683 "SB03MX.f"
			for (j = 1; j <= i__1; ++j) {
#line 684 "SB03MX.f"
			    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 685 "SB03MX.f"
/* L100: */
#line 685 "SB03MX.f"
			}

#line 687 "SB03MX.f"
			dscal_(n, &scaloc, &dwork[1], &c__1);
#line 688 "SB03MX.f"
			dscal_(n, &scaloc, &dwork[np1], &c__1);
#line 689 "SB03MX.f"
			*scale *= scaloc;
#line 690 "SB03MX.f"
		    }
#line 691 "SB03MX.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 692 "SB03MX.f"
		    c__[k1 + l2 * c_dim1] = x[2];
#line 693 "SB03MX.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 694 "SB03MX.f"
		    c__[k2 + l2 * c_dim1] = x[3];
#line 695 "SB03MX.f"
		    if (k1 != l1) {
#line 696 "SB03MX.f"
			c__[l1 + k1 * c_dim1] = x[0];
#line 697 "SB03MX.f"
			c__[l2 + k1 * c_dim1] = x[2];
#line 698 "SB03MX.f"
			c__[l1 + k2 * c_dim1] = x[1];
#line 699 "SB03MX.f"
			c__[l2 + k2 * c_dim1] = x[3];
#line 700 "SB03MX.f"
		    }
#line 701 "SB03MX.f"
		}

#line 703 "SB03MX.f"
L110:
#line 703 "SB03MX.f"
		;
#line 703 "SB03MX.f"
	    }

#line 705 "SB03MX.f"
L120:
#line 705 "SB03MX.f"
	    ;
#line 705 "SB03MX.f"
	}

#line 707 "SB03MX.f"
    }

#line 709 "SB03MX.f"
    return 0;
/* *** Last line of SB03MX *** */
} /* sb03mx_ */

