#line 1 "SB04PY.f"
/* SB04PY.f -- translated by f2c (version 20100827).
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

#line 1 "SB04PY.f"
/* Table of constant values */

static integer c__1 = 1;
static logical c_false = FALSE_;
static integer c__2 = 2;
static doublereal c_b28 = 1.;
static doublereal c_b31 = 0.;
static logical c_true = TRUE_;

/* Subroutine */ int sb04py_(char *trana, char *tranb, integer *isgn, integer 
	*m, integer *n, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *scale, doublereal *
	dwork, integer *info, ftnlen trana_len, ftnlen tranb_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j, k, l;
    static doublereal x[4]	/* was [2][2] */;
    static integer k1, k2, l1, l2;
    static doublereal a11, db, p11, p12, p21, p22, da11, vec[4]	/* was [2][2] 
	    */, dum[1], eps, sgn;
    static integer mnk1, mnk2, mnl1, mnl2;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ierr;
    static doublereal smin, sumr;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sb04px_(logical *, logical *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer knext, lnext;
    static doublereal xnorm;
    extern /* Subroutine */ int dlaln2_(logical *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *),
	     dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    static doublereal scaloc;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static logical notrna, notrnb;
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

/*     To solve for X the discrete-time Sylvester equation */

/*        op(A)*X*op(B) + ISGN*X = scale*C, */

/*     where op(A) = A or A**T, A and B are both upper quasi-triangular, */
/*     and ISGN = 1 or -1. A is M-by-M and B is N-by-N; the right hand */
/*     side C and the solution X are M-by-N; and scale is an output scale */
/*     factor, set less than or equal to 1 to avoid overflow in X. The */
/*     solution matrix X is overwritten onto C. */

/*     A and B must be in Schur canonical form (as returned by LAPACK */
/*     Library routine DHSEQR), that is, block upper triangular with */
/*     1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block has */
/*     its diagonal elements equal and its off-diagonal elements of */
/*     opposite sign. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op(A) to be used, as follows: */
/*             = 'N':  op(A) = A    (No transpose); */
/*             = 'T':  op(A) = A**T (Transpose); */
/*             = 'C':  op(A) = A**T (Conjugate transpose = Transpose). */

/*     TRANB   CHARACTER*1 */
/*             Specifies the form of op(B) to be used, as follows: */
/*             = 'N':  op(B) = B    (No transpose); */
/*             = 'T':  op(B) = B**T (Transpose); */
/*             = 'C':  op(B) = B**T (Conjugate transpose = Transpose). */

/*     ISGN    INTEGER */
/*             Specifies the sign of the equation as described before. */
/*             ISGN may only be 1 or -1. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of the matrix A, and the number of rows in the */
/*             matrices X and C.  M >= 0. */

/*     N       (input) INTEGER */
/*             The order of the matrix B, and the number of columns in */
/*             the matrices X and C.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,M) */
/*             The leading M-by-M part of this array must contain the */
/*             upper quasi-triangular matrix A, in Schur canonical form. */
/*             The part of A below the first sub-diagonal is not */
/*             referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,N) */
/*             The leading N-by-N part of this array must contain the */
/*             upper quasi-triangular matrix B, in Schur canonical form. */
/*             The part of B below the first sub-diagonal is not */
/*             referenced. */

/*     LDB     (input) INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the right hand side matrix C. */
/*             On exit, if INFO >= 0, the leading M-by-N part of this */
/*             array contains the solution matrix X. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,M). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor, scale, set less than or equal to 1 to */
/*             prevent the solution overflowing. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (2*M) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  A and -ISGN*B have almost reciprocal eigenvalues; */
/*                   perturbed values were used to solve the equation */
/*                   (but the matrices A and B are unchanged). */

/*     METHOD */

/*     The solution matrix X is computed column-wise via a back */
/*     substitution scheme, an extension and refinement of the algorithm */
/*     in [1], similar to that used in [2] for continuous-time Sylvester */
/*     equations. A set of equivalent linear algebraic systems of */
/*     equations of order at most four are formed and solved using */
/*     Gaussian elimination with complete pivoting. */

/*     REFERENCES */

/*     [1] Bartels, R.H. and Stewart, G.W.  T */
/*         Solution of the matrix equation A X + XB = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     [2] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J., */
/*         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A., */
/*         Ostrouchov, S., and Sorensen, D. */
/*         LAPACK Users' Guide: Second Edition. */
/*         SIAM, Philadelphia, 1995. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is stable and reliable, since Gaussian elimination */
/*     with complete pivoting is used. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, March 2000. */
/*     D. Sima, University of Bucharest, April 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000. */
/*     Partly based on the routine SYLSV, A. Varga, 1992. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Discrete-time system, matrix algebra, Sylvester equation. */

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

/*     Decode and Test input parameters */

#line 196 "SB04PY.f"
    /* Parameter adjustments */
#line 196 "SB04PY.f"
    a_dim1 = *lda;
#line 196 "SB04PY.f"
    a_offset = 1 + a_dim1;
#line 196 "SB04PY.f"
    a -= a_offset;
#line 196 "SB04PY.f"
    b_dim1 = *ldb;
#line 196 "SB04PY.f"
    b_offset = 1 + b_dim1;
#line 196 "SB04PY.f"
    b -= b_offset;
#line 196 "SB04PY.f"
    c_dim1 = *ldc;
#line 196 "SB04PY.f"
    c_offset = 1 + c_dim1;
#line 196 "SB04PY.f"
    c__ -= c_offset;
#line 196 "SB04PY.f"
    --dwork;
#line 196 "SB04PY.f"

#line 196 "SB04PY.f"
    /* Function Body */
#line 196 "SB04PY.f"
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
#line 197 "SB04PY.f"
    notrnb = lsame_(tranb, "N", (ftnlen)1, (ftnlen)1);

#line 199 "SB04PY.f"
    *info = 0;
#line 200 "SB04PY.f"
    if (! notrna && ! lsame_(trana, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trana, "C", (ftnlen)1, (ftnlen)1)) {
#line 202 "SB04PY.f"
	*info = -1;
#line 203 "SB04PY.f"
    } else if (! notrnb && ! lsame_(tranb, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(tranb, "C", (ftnlen)1, (ftnlen)1)) {
#line 205 "SB04PY.f"
	*info = -2;
#line 206 "SB04PY.f"
    } else if (*isgn != 1 && *isgn != -1) {
#line 207 "SB04PY.f"
	*info = -3;
#line 208 "SB04PY.f"
    } else if (*m < 0) {
#line 209 "SB04PY.f"
	*info = -4;
#line 210 "SB04PY.f"
    } else if (*n < 0) {
#line 211 "SB04PY.f"
	*info = -5;
#line 212 "SB04PY.f"
    } else if (*lda < max(1,*m)) {
#line 213 "SB04PY.f"
	*info = -7;
#line 214 "SB04PY.f"
    } else if (*ldb < max(1,*n)) {
#line 215 "SB04PY.f"
	*info = -9;
#line 216 "SB04PY.f"
    } else if (*ldc < max(1,*m)) {
#line 217 "SB04PY.f"
	*info = -11;
#line 218 "SB04PY.f"
    }
#line 219 "SB04PY.f"
    if (*info != 0) {
#line 220 "SB04PY.f"
	i__1 = -(*info);
#line 220 "SB04PY.f"
	xerbla_("SB04PY", &i__1, (ftnlen)6);
#line 221 "SB04PY.f"
	return 0;
#line 222 "SB04PY.f"
    }

/*     Quick return if possible. */

#line 226 "SB04PY.f"
    *scale = 1.;
#line 227 "SB04PY.f"
    if (*m == 0 || *n == 0) {
#line 227 "SB04PY.f"
	return 0;
#line 227 "SB04PY.f"
    }

/*     Set constants to control overflow. */

#line 232 "SB04PY.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 233 "SB04PY.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12);
#line 234 "SB04PY.f"
    bignum = 1. / smlnum;
#line 235 "SB04PY.f"
    dlabad_(&smlnum, &bignum);
#line 236 "SB04PY.f"
    smlnum = smlnum * (doublereal) (*m * *n) / eps;
#line 237 "SB04PY.f"
    bignum = 1. / smlnum;

/* Computing MAX */
#line 239 "SB04PY.f"
    d__1 = smlnum, d__2 = eps * dlange_("M", m, m, &a[a_offset], lda, dum, (
	    ftnlen)1), d__1 = max(d__1,d__2), d__2 = eps * dlange_("M", n, n, 
	    &b[b_offset], ldb, dum, (ftnlen)1);
#line 239 "SB04PY.f"
    smin = max(d__1,d__2);

#line 242 "SB04PY.f"
    sgn = (doublereal) (*isgn);

#line 244 "SB04PY.f"
    if (notrna && notrnb) {

/*        Solve    A*X*B + ISGN*X = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        bottom-left corner column by column by */

/*           A(K,K)*X(K,L)*B(L,L) + ISGN*X(K,L) = C(K,L) - R(K,L) */

/*        where */
/*                       M */
/*           R(K,L) = { SUM [A(K,J)*X(J,L)] } * B(L,L) + */
/*                     J=K+1 */
/*                       M             L-1 */
/*                      SUM { A(K,J) * SUM [X(J,I)*B(I,L)] }. */
/*                      J=K            I=1 */

/*        Start column loop (index = L) */
/*        L1 (L2) : column index of the first (last) row of X(K,L). */

#line 264 "SB04PY.f"
	lnext = 1;

#line 266 "SB04PY.f"
	i__1 = *n;
#line 266 "SB04PY.f"
	for (l = 1; l <= i__1; ++l) {
#line 267 "SB04PY.f"
	    if (l < lnext) {
#line 267 "SB04PY.f"
		goto L60;
#line 267 "SB04PY.f"
	    }
#line 269 "SB04PY.f"
	    l1 = l;
#line 270 "SB04PY.f"
	    if (l == *n) {
#line 271 "SB04PY.f"
		l2 = l;
#line 272 "SB04PY.f"
	    } else {
#line 273 "SB04PY.f"
		if (b[l + 1 + l * b_dim1] != 0.) {
#line 274 "SB04PY.f"
		    l2 = l + 1;
#line 275 "SB04PY.f"
		} else {
#line 276 "SB04PY.f"
		    l2 = l;
#line 277 "SB04PY.f"
		}
#line 278 "SB04PY.f"
		lnext = l2 + 1;
#line 279 "SB04PY.f"
	    }

/*           Start row loop (index = K) */
/*           K1 (K2): row index of the first (last) row of X(K,L). */

#line 284 "SB04PY.f"
	    knext = *m;

#line 286 "SB04PY.f"
	    for (k = *m; k >= 1; --k) {
#line 287 "SB04PY.f"
		if (k > knext) {
#line 287 "SB04PY.f"
		    goto L50;
#line 287 "SB04PY.f"
		}
#line 289 "SB04PY.f"
		k2 = k;
#line 290 "SB04PY.f"
		if (k == 1) {
#line 291 "SB04PY.f"
		    k1 = k;
#line 292 "SB04PY.f"
		} else {
#line 293 "SB04PY.f"
		    if (a[k + (k - 1) * a_dim1] != 0.) {
#line 294 "SB04PY.f"
			k1 = k - 1;
#line 295 "SB04PY.f"
		    } else {
#line 296 "SB04PY.f"
			k1 = k;
#line 297 "SB04PY.f"
		    }
#line 298 "SB04PY.f"
		    knext = k1 - 1;
#line 299 "SB04PY.f"
		}

/* Computing MIN */
#line 301 "SB04PY.f"
		i__2 = k1 + 1;
#line 301 "SB04PY.f"
		mnk1 = min(i__2,*m);
/* Computing MIN */
#line 302 "SB04PY.f"
		i__2 = k2 + 1;
#line 302 "SB04PY.f"
		mnk2 = min(i__2,*m);
#line 303 "SB04PY.f"
		i__2 = *m - k2;
#line 303 "SB04PY.f"
		p11 = ddot_(&i__2, &a[k1 + mnk2 * a_dim1], lda, &c__[mnk2 + 
			l1 * c_dim1], &c__1);
#line 304 "SB04PY.f"
		i__2 = l1 - 1;
#line 304 "SB04PY.f"
		dwork[k1] = ddot_(&i__2, &c__[k1 + c_dim1], ldc, &b[l1 * 
			b_dim1 + 1], &c__1);

#line 307 "SB04PY.f"
		if (l1 == l2 && k1 == k2) {

#line 309 "SB04PY.f"
		    i__2 = *m - k1 + 1;
#line 309 "SB04PY.f"
		    sumr = ddot_(&i__2, &a[k1 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
#line 311 "SB04PY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1]);
#line 312 "SB04PY.f"
		    scaloc = 1.;

#line 314 "SB04PY.f"
		    a11 = a[k1 + k1 * a_dim1] * b[l1 + l1 * b_dim1] + sgn;
#line 315 "SB04PY.f"
		    da11 = abs(a11);
#line 316 "SB04PY.f"
		    if (da11 <= smin) {
#line 317 "SB04PY.f"
			a11 = smin;
#line 318 "SB04PY.f"
			da11 = smin;
#line 319 "SB04PY.f"
			*info = 1;
#line 320 "SB04PY.f"
		    }
#line 321 "SB04PY.f"
		    db = abs(vec[0]);
#line 322 "SB04PY.f"
		    if (da11 < 1. && db > 1.) {
#line 323 "SB04PY.f"
			if (db > bignum * da11) {
#line 323 "SB04PY.f"
			    scaloc = 1. / db;
#line 323 "SB04PY.f"
			}
#line 325 "SB04PY.f"
		    }
#line 326 "SB04PY.f"
		    x[0] = vec[0] * scaloc / a11;

#line 328 "SB04PY.f"
		    if (scaloc != 1.) {

#line 330 "SB04PY.f"
			i__2 = *n;
#line 330 "SB04PY.f"
			for (j = 1; j <= i__2; ++j) {
#line 331 "SB04PY.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 332 "SB04PY.f"
/* L10: */
#line 332 "SB04PY.f"
			}

#line 334 "SB04PY.f"
			i__2 = *m - k1 + 1;
#line 334 "SB04PY.f"
			dscal_(&i__2, &scaloc, &dwork[k1], &c__1);
#line 335 "SB04PY.f"
			*scale *= scaloc;
#line 336 "SB04PY.f"
		    }
#line 337 "SB04PY.f"
		    c__[k1 + l1 * c_dim1] = x[0];

#line 339 "SB04PY.f"
		} else if (l1 == l2 && k1 != k2) {

#line 341 "SB04PY.f"
		    i__2 = *m - k2;
#line 341 "SB04PY.f"
		    p21 = ddot_(&i__2, &a[k2 + mnk2 * a_dim1], lda, &c__[mnk2 
			    + l1 * c_dim1], &c__1);
#line 343 "SB04PY.f"
		    i__2 = l1 - 1;
#line 343 "SB04PY.f"
		    dwork[k2] = ddot_(&i__2, &c__[k2 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 345 "SB04PY.f"
		    i__2 = *m - k1 + 1;
#line 345 "SB04PY.f"
		    sumr = ddot_(&i__2, &a[k1 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
#line 347 "SB04PY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1]);

#line 349 "SB04PY.f"
		    i__2 = *m - k1 + 1;
#line 349 "SB04PY.f"
		    sumr = ddot_(&i__2, &a[k2 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
#line 351 "SB04PY.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (sumr + p21 * b[l1 + l1 *
			     b_dim1]);

#line 353 "SB04PY.f"
		    d__1 = -sgn;
#line 353 "SB04PY.f"
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &b[l1 + l1 * 
			    b_dim1], &a[k1 + k1 * a_dim1], lda, &c_b28, &
			    c_b28, vec, &c__2, &d__1, &c_b31, x, &c__2, &
			    scaloc, &xnorm, &ierr);
#line 356 "SB04PY.f"
		    if (ierr != 0) {
#line 356 "SB04PY.f"
			*info = 1;
#line 356 "SB04PY.f"
		    }

#line 359 "SB04PY.f"
		    if (scaloc != 1.) {

#line 361 "SB04PY.f"
			i__2 = *n;
#line 361 "SB04PY.f"
			for (j = 1; j <= i__2; ++j) {
#line 362 "SB04PY.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 363 "SB04PY.f"
/* L20: */
#line 363 "SB04PY.f"
			}

#line 365 "SB04PY.f"
			i__2 = *m - k1 + 1;
#line 365 "SB04PY.f"
			dscal_(&i__2, &scaloc, &dwork[k1], &c__1);
#line 366 "SB04PY.f"
			*scale *= scaloc;
#line 367 "SB04PY.f"
		    }
#line 368 "SB04PY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 369 "SB04PY.f"
		    c__[k2 + l1 * c_dim1] = x[1];

#line 371 "SB04PY.f"
		} else if (l1 != l2 && k1 == k2) {

#line 373 "SB04PY.f"
		    i__2 = *m - k1;
#line 373 "SB04PY.f"
		    p12 = ddot_(&i__2, &a[k1 + mnk1 * a_dim1], lda, &c__[mnk1 
			    + l2 * c_dim1], &c__1);
#line 375 "SB04PY.f"
		    i__2 = *m - k1 + 1;
#line 375 "SB04PY.f"
		    sumr = ddot_(&i__2, &a[k1 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
#line 377 "SB04PY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1] + p12 * b[l2 + l1 * b_dim1]);

#line 380 "SB04PY.f"
		    i__2 = l1 - 1;
#line 380 "SB04PY.f"
		    dwork[k1 + *m] = ddot_(&i__2, &c__[k1 + c_dim1], ldc, &b[
			    l2 * b_dim1 + 1], &c__1);
#line 382 "SB04PY.f"
		    i__2 = *m - k1 + 1;
#line 382 "SB04PY.f"
		    sumr = ddot_(&i__2, &a[k1 + k1 * a_dim1], lda, &dwork[k1 
			    + *m], &c__1);
#line 384 "SB04PY.f"
		    vec[1] = c__[k1 + l2 * c_dim1] - (sumr + p11 * b[l1 + l2 *
			     b_dim1] + p12 * b[l2 + l2 * b_dim1]);

#line 387 "SB04PY.f"
		    d__1 = -sgn;
#line 387 "SB04PY.f"
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &a[k1 + k1 * a_dim1]
			    , &b[l1 + l1 * b_dim1], ldb, &c_b28, &c_b28, vec, 
			    &c__2, &d__1, &c_b31, x, &c__2, &scaloc, &xnorm, &
			    ierr);
#line 390 "SB04PY.f"
		    if (ierr != 0) {
#line 390 "SB04PY.f"
			*info = 1;
#line 390 "SB04PY.f"
		    }

#line 393 "SB04PY.f"
		    if (scaloc != 1.) {

#line 395 "SB04PY.f"
			i__2 = *n;
#line 395 "SB04PY.f"
			for (j = 1; j <= i__2; ++j) {
#line 396 "SB04PY.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 397 "SB04PY.f"
/* L30: */
#line 397 "SB04PY.f"
			}

#line 399 "SB04PY.f"
			i__2 = *m - k1 + 1;
#line 399 "SB04PY.f"
			dscal_(&i__2, &scaloc, &dwork[k1], &c__1);
#line 400 "SB04PY.f"
			i__2 = *m - k1 + 1;
#line 400 "SB04PY.f"
			dscal_(&i__2, &scaloc, &dwork[k1 + *m], &c__1);
#line 401 "SB04PY.f"
			*scale *= scaloc;
#line 402 "SB04PY.f"
		    }
#line 403 "SB04PY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 404 "SB04PY.f"
		    c__[k1 + l2 * c_dim1] = x[1];

#line 406 "SB04PY.f"
		} else if (l1 != l2 && k1 != k2) {

#line 408 "SB04PY.f"
		    i__2 = *m - k2;
#line 408 "SB04PY.f"
		    p21 = ddot_(&i__2, &a[k2 + mnk2 * a_dim1], lda, &c__[mnk2 
			    + l1 * c_dim1], &c__1);
#line 410 "SB04PY.f"
		    i__2 = *m - k2;
#line 410 "SB04PY.f"
		    p12 = ddot_(&i__2, &a[k1 + mnk2 * a_dim1], lda, &c__[mnk2 
			    + l2 * c_dim1], &c__1);
#line 412 "SB04PY.f"
		    i__2 = *m - k2;
#line 412 "SB04PY.f"
		    p22 = ddot_(&i__2, &a[k2 + mnk2 * a_dim1], lda, &c__[mnk2 
			    + l2 * c_dim1], &c__1);

#line 415 "SB04PY.f"
		    i__2 = l1 - 1;
#line 415 "SB04PY.f"
		    dwork[k2] = ddot_(&i__2, &c__[k2 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 417 "SB04PY.f"
		    i__2 = *m - k1 + 1;
#line 417 "SB04PY.f"
		    sumr = ddot_(&i__2, &a[k1 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
#line 419 "SB04PY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1] + p12 * b[l2 + l1 * b_dim1]);

#line 422 "SB04PY.f"
		    i__2 = l1 - 1;
#line 422 "SB04PY.f"
		    dwork[k1 + *m] = ddot_(&i__2, &c__[k1 + c_dim1], ldc, &b[
			    l2 * b_dim1 + 1], &c__1);
#line 424 "SB04PY.f"
		    i__2 = l1 - 1;
#line 424 "SB04PY.f"
		    dwork[k2 + *m] = ddot_(&i__2, &c__[k2 + c_dim1], ldc, &b[
			    l2 * b_dim1 + 1], &c__1);
#line 426 "SB04PY.f"
		    i__2 = *m - k1 + 1;
#line 426 "SB04PY.f"
		    sumr = ddot_(&i__2, &a[k1 + k1 * a_dim1], lda, &dwork[k1 
			    + *m], &c__1);
#line 428 "SB04PY.f"
		    vec[2] = c__[k1 + l2 * c_dim1] - (sumr + p11 * b[l1 + l2 *
			     b_dim1] + p12 * b[l2 + l2 * b_dim1]);

#line 431 "SB04PY.f"
		    i__2 = *m - k1 + 1;
#line 431 "SB04PY.f"
		    sumr = ddot_(&i__2, &a[k2 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
#line 433 "SB04PY.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (sumr + p21 * b[l1 + l1 *
			     b_dim1] + p22 * b[l2 + l1 * b_dim1]);

#line 436 "SB04PY.f"
		    i__2 = *m - k1 + 1;
#line 436 "SB04PY.f"
		    sumr = ddot_(&i__2, &a[k2 + k1 * a_dim1], lda, &dwork[k1 
			    + *m], &c__1);
#line 438 "SB04PY.f"
		    vec[3] = c__[k2 + l2 * c_dim1] - (sumr + p21 * b[l1 + l2 *
			     b_dim1] + p22 * b[l2 + l2 * b_dim1]);

#line 441 "SB04PY.f"
		    sb04px_(&c_false, &c_false, isgn, &c__2, &c__2, &a[k1 + 
			    k1 * a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec,
			     &c__2, &scaloc, x, &c__2, &xnorm, &ierr);
#line 444 "SB04PY.f"
		    if (ierr != 0) {
#line 444 "SB04PY.f"
			*info = 1;
#line 444 "SB04PY.f"
		    }

#line 447 "SB04PY.f"
		    if (scaloc != 1.) {

#line 449 "SB04PY.f"
			i__2 = *n;
#line 449 "SB04PY.f"
			for (j = 1; j <= i__2; ++j) {
#line 450 "SB04PY.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 451 "SB04PY.f"
/* L40: */
#line 451 "SB04PY.f"
			}

#line 453 "SB04PY.f"
			i__2 = *m - k1 + 1;
#line 453 "SB04PY.f"
			dscal_(&i__2, &scaloc, &dwork[k1], &c__1);
#line 454 "SB04PY.f"
			i__2 = *m - k1 + 1;
#line 454 "SB04PY.f"
			dscal_(&i__2, &scaloc, &dwork[k1 + *m], &c__1);
#line 455 "SB04PY.f"
			*scale *= scaloc;
#line 456 "SB04PY.f"
		    }
#line 457 "SB04PY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 458 "SB04PY.f"
		    c__[k1 + l2 * c_dim1] = x[2];
#line 459 "SB04PY.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 460 "SB04PY.f"
		    c__[k2 + l2 * c_dim1] = x[3];
#line 461 "SB04PY.f"
		}

#line 463 "SB04PY.f"
L50:
#line 463 "SB04PY.f"
		;
#line 463 "SB04PY.f"
	    }

#line 465 "SB04PY.f"
L60:
#line 465 "SB04PY.f"
	    ;
#line 465 "SB04PY.f"
	}

#line 467 "SB04PY.f"
    } else if (! notrna && notrnb) {

/*        Solve     A'*X*B + ISGN*X = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        upper-left corner column by column by */

/*         A(K,K)'*X(K,L)*B(L,L) + ISGN*X(K,L) = C(K,L) - R(K,L) */

/*        where */
/*                      K-1 */
/*           R(K,L) = { SUM [A(J,K)'*X(J,L)] } * B(L,L) + */
/*                      J=1 */
/*                       K              L-1 */
/*                      SUM A(J,K)' * { SUM [X(J,I)*B(I,L)] }. */
/*                      J=1             I=1 */

/*        Start column loop (index = L) */
/*        L1 (L2): column index of the first (last) row of X(K,L). */

#line 487 "SB04PY.f"
	lnext = 1;

#line 489 "SB04PY.f"
	i__1 = *n;
#line 489 "SB04PY.f"
	for (l = 1; l <= i__1; ++l) {
#line 490 "SB04PY.f"
	    if (l < lnext) {
#line 490 "SB04PY.f"
		goto L120;
#line 490 "SB04PY.f"
	    }
#line 492 "SB04PY.f"
	    l1 = l;
#line 493 "SB04PY.f"
	    if (l == *n) {
#line 494 "SB04PY.f"
		l2 = l;
#line 495 "SB04PY.f"
	    } else {
#line 496 "SB04PY.f"
		if (b[l + 1 + l * b_dim1] != 0.) {
#line 497 "SB04PY.f"
		    l2 = l + 1;
#line 498 "SB04PY.f"
		} else {
#line 499 "SB04PY.f"
		    l2 = l;
#line 500 "SB04PY.f"
		}
#line 501 "SB04PY.f"
		lnext = l2 + 1;
#line 502 "SB04PY.f"
	    }

/*           Start row loop (index = K) */
/*           K1 (K2): row index of the first (last) row of X(K,L). */

#line 507 "SB04PY.f"
	    knext = 1;

#line 509 "SB04PY.f"
	    i__2 = *m;
#line 509 "SB04PY.f"
	    for (k = 1; k <= i__2; ++k) {
#line 510 "SB04PY.f"
		if (k < knext) {
#line 510 "SB04PY.f"
		    goto L110;
#line 510 "SB04PY.f"
		}
#line 512 "SB04PY.f"
		k1 = k;
#line 513 "SB04PY.f"
		if (k == *m) {
#line 514 "SB04PY.f"
		    k2 = k;
#line 515 "SB04PY.f"
		} else {
#line 516 "SB04PY.f"
		    if (a[k + 1 + k * a_dim1] != 0.) {
#line 517 "SB04PY.f"
			k2 = k + 1;
#line 518 "SB04PY.f"
		    } else {
#line 519 "SB04PY.f"
			k2 = k;
#line 520 "SB04PY.f"
		    }
#line 521 "SB04PY.f"
		    knext = k2 + 1;
#line 522 "SB04PY.f"
		}

#line 524 "SB04PY.f"
		i__3 = k1 - 1;
#line 524 "SB04PY.f"
		p11 = ddot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			c_dim1 + 1], &c__1);
#line 525 "SB04PY.f"
		i__3 = l1 - 1;
#line 525 "SB04PY.f"
		dwork[k1] = ddot_(&i__3, &c__[k1 + c_dim1], ldc, &b[l1 * 
			b_dim1 + 1], &c__1);

#line 528 "SB04PY.f"
		if (l1 == l2 && k1 == k2) {

#line 530 "SB04PY.f"
		    sumr = ddot_(&k1, &a[k1 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
#line 531 "SB04PY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1]);
#line 532 "SB04PY.f"
		    scaloc = 1.;

#line 534 "SB04PY.f"
		    a11 = a[k1 + k1 * a_dim1] * b[l1 + l1 * b_dim1] + sgn;
#line 535 "SB04PY.f"
		    da11 = abs(a11);
#line 536 "SB04PY.f"
		    if (da11 <= smin) {
#line 537 "SB04PY.f"
			a11 = smin;
#line 538 "SB04PY.f"
			da11 = smin;
#line 539 "SB04PY.f"
			*info = 1;
#line 540 "SB04PY.f"
		    }
#line 541 "SB04PY.f"
		    db = abs(vec[0]);
#line 542 "SB04PY.f"
		    if (da11 < 1. && db > 1.) {
#line 543 "SB04PY.f"
			if (db > bignum * da11) {
#line 543 "SB04PY.f"
			    scaloc = 1. / db;
#line 543 "SB04PY.f"
			}
#line 545 "SB04PY.f"
		    }
#line 546 "SB04PY.f"
		    x[0] = vec[0] * scaloc / a11;

#line 548 "SB04PY.f"
		    if (scaloc != 1.) {

#line 550 "SB04PY.f"
			i__3 = *n;
#line 550 "SB04PY.f"
			for (j = 1; j <= i__3; ++j) {
#line 551 "SB04PY.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 552 "SB04PY.f"
/* L70: */
#line 552 "SB04PY.f"
			}

#line 554 "SB04PY.f"
			dscal_(&k1, &scaloc, &dwork[1], &c__1);
#line 555 "SB04PY.f"
			*scale *= scaloc;
#line 556 "SB04PY.f"
		    }
#line 557 "SB04PY.f"
		    c__[k1 + l1 * c_dim1] = x[0];

#line 559 "SB04PY.f"
		} else if (l1 == l2 && k1 != k2) {

#line 561 "SB04PY.f"
		    i__3 = k1 - 1;
#line 561 "SB04PY.f"
		    p21 = ddot_(&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 562 "SB04PY.f"
		    i__3 = l1 - 1;
#line 562 "SB04PY.f"
		    dwork[k2] = ddot_(&i__3, &c__[k2 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 564 "SB04PY.f"
		    sumr = ddot_(&k2, &a[k1 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
#line 565 "SB04PY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1]);

#line 567 "SB04PY.f"
		    sumr = ddot_(&k2, &a[k2 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
#line 568 "SB04PY.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (sumr + p21 * b[l1 + l1 *
			     b_dim1]);

#line 570 "SB04PY.f"
		    d__1 = -sgn;
#line 570 "SB04PY.f"
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &b[l1 + l1 * b_dim1]
			    , &a[k1 + k1 * a_dim1], lda, &c_b28, &c_b28, vec, 
			    &c__2, &d__1, &c_b31, x, &c__2, &scaloc, &xnorm, &
			    ierr);
#line 573 "SB04PY.f"
		    if (ierr != 0) {
#line 573 "SB04PY.f"
			*info = 1;
#line 573 "SB04PY.f"
		    }

#line 576 "SB04PY.f"
		    if (scaloc != 1.) {

#line 578 "SB04PY.f"
			i__3 = *n;
#line 578 "SB04PY.f"
			for (j = 1; j <= i__3; ++j) {
#line 579 "SB04PY.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 580 "SB04PY.f"
/* L80: */
#line 580 "SB04PY.f"
			}

#line 582 "SB04PY.f"
			dscal_(&k2, &scaloc, &dwork[1], &c__1);
#line 583 "SB04PY.f"
			*scale *= scaloc;
#line 584 "SB04PY.f"
		    }
#line 585 "SB04PY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 586 "SB04PY.f"
		    c__[k2 + l1 * c_dim1] = x[1];

#line 588 "SB04PY.f"
		} else if (l1 != l2 && k1 == k2) {

#line 590 "SB04PY.f"
		    i__3 = k1 - 1;
#line 590 "SB04PY.f"
		    p12 = ddot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
#line 591 "SB04PY.f"
		    sumr = ddot_(&k1, &a[k1 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
#line 592 "SB04PY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1] + p12 * b[l2 + l1 * b_dim1]);

#line 595 "SB04PY.f"
		    i__3 = l1 - 1;
#line 595 "SB04PY.f"
		    dwork[k1 + *m] = ddot_(&i__3, &c__[k1 + c_dim1], ldc, &b[
			    l2 * b_dim1 + 1], &c__1);
#line 597 "SB04PY.f"
		    sumr = ddot_(&k1, &a[k1 * a_dim1 + 1], &c__1, &dwork[*m + 
			    1], &c__1);
#line 598 "SB04PY.f"
		    vec[1] = c__[k1 + l2 * c_dim1] - (sumr + p11 * b[l1 + l2 *
			     b_dim1] + p12 * b[l2 + l2 * b_dim1]);

#line 601 "SB04PY.f"
		    d__1 = -sgn;
#line 601 "SB04PY.f"
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &a[k1 + k1 * a_dim1]
			    , &b[l1 + l1 * b_dim1], ldb, &c_b28, &c_b28, vec, 
			    &c__2, &d__1, &c_b31, x, &c__2, &scaloc, &xnorm, &
			    ierr);
#line 604 "SB04PY.f"
		    if (ierr != 0) {
#line 604 "SB04PY.f"
			*info = 1;
#line 604 "SB04PY.f"
		    }

#line 607 "SB04PY.f"
		    if (scaloc != 1.) {

#line 609 "SB04PY.f"
			i__3 = *n;
#line 609 "SB04PY.f"
			for (j = 1; j <= i__3; ++j) {
#line 610 "SB04PY.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 611 "SB04PY.f"
/* L90: */
#line 611 "SB04PY.f"
			}

#line 613 "SB04PY.f"
			dscal_(&k1, &scaloc, &dwork[1], &c__1);
#line 614 "SB04PY.f"
			dscal_(&k1, &scaloc, &dwork[*m + 1], &c__1);
#line 615 "SB04PY.f"
			*scale *= scaloc;
#line 616 "SB04PY.f"
		    }
#line 617 "SB04PY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 618 "SB04PY.f"
		    c__[k1 + l2 * c_dim1] = x[1];

#line 620 "SB04PY.f"
		} else if (l1 != l2 && k1 != k2) {

#line 622 "SB04PY.f"
		    i__3 = k1 - 1;
#line 622 "SB04PY.f"
		    p21 = ddot_(&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 623 "SB04PY.f"
		    i__3 = k1 - 1;
#line 623 "SB04PY.f"
		    p12 = ddot_(&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
#line 624 "SB04PY.f"
		    i__3 = k1 - 1;
#line 624 "SB04PY.f"
		    p22 = ddot_(&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);

#line 626 "SB04PY.f"
		    i__3 = l1 - 1;
#line 626 "SB04PY.f"
		    dwork[k2] = ddot_(&i__3, &c__[k2 + c_dim1], ldc, &b[l1 * 
			    b_dim1 + 1], &c__1);
#line 628 "SB04PY.f"
		    sumr = ddot_(&k2, &a[k1 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
#line 629 "SB04PY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1] + p12 * b[l2 + l1 * b_dim1]);

#line 632 "SB04PY.f"
		    sumr = ddot_(&k2, &a[k2 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
#line 633 "SB04PY.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (sumr + p21 * b[l1 + l1 *
			     b_dim1] + p22 * b[l2 + l1 * b_dim1]);

#line 636 "SB04PY.f"
		    i__3 = l1 - 1;
#line 636 "SB04PY.f"
		    dwork[k1 + *m] = ddot_(&i__3, &c__[k1 + c_dim1], ldc, &b[
			    l2 * b_dim1 + 1], &c__1);
#line 638 "SB04PY.f"
		    i__3 = l1 - 1;
#line 638 "SB04PY.f"
		    dwork[k2 + *m] = ddot_(&i__3, &c__[k2 + c_dim1], ldc, &b[
			    l2 * b_dim1 + 1], &c__1);
#line 640 "SB04PY.f"
		    sumr = ddot_(&k2, &a[k1 * a_dim1 + 1], &c__1, &dwork[*m + 
			    1], &c__1);
#line 641 "SB04PY.f"
		    vec[2] = c__[k1 + l2 * c_dim1] - (sumr + p11 * b[l1 + l2 *
			     b_dim1] + p12 * b[l2 + l2 * b_dim1]);

#line 644 "SB04PY.f"
		    sumr = ddot_(&k2, &a[k2 * a_dim1 + 1], &c__1, &dwork[*m + 
			    1], &c__1);
#line 645 "SB04PY.f"
		    vec[3] = c__[k2 + l2 * c_dim1] - (sumr + p21 * b[l1 + l2 *
			     b_dim1] + p22 * b[l2 + l2 * b_dim1]);

#line 648 "SB04PY.f"
		    sb04px_(&c_true, &c_false, isgn, &c__2, &c__2, &a[k1 + k1 
			    * a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec, &
			    c__2, &scaloc, x, &c__2, &xnorm, &ierr);
#line 651 "SB04PY.f"
		    if (ierr != 0) {
#line 651 "SB04PY.f"
			*info = 1;
#line 651 "SB04PY.f"
		    }

#line 654 "SB04PY.f"
		    if (scaloc != 1.) {

#line 656 "SB04PY.f"
			i__3 = *n;
#line 656 "SB04PY.f"
			for (j = 1; j <= i__3; ++j) {
#line 657 "SB04PY.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 658 "SB04PY.f"
/* L100: */
#line 658 "SB04PY.f"
			}

#line 660 "SB04PY.f"
			dscal_(&k2, &scaloc, &dwork[1], &c__1);
#line 661 "SB04PY.f"
			dscal_(&k2, &scaloc, &dwork[*m + 1], &c__1);
#line 662 "SB04PY.f"
			*scale *= scaloc;
#line 663 "SB04PY.f"
		    }
#line 664 "SB04PY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 665 "SB04PY.f"
		    c__[k1 + l2 * c_dim1] = x[2];
#line 666 "SB04PY.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 667 "SB04PY.f"
		    c__[k2 + l2 * c_dim1] = x[3];
#line 668 "SB04PY.f"
		}

#line 670 "SB04PY.f"
L110:
#line 670 "SB04PY.f"
		;
#line 670 "SB04PY.f"
	    }

#line 672 "SB04PY.f"
L120:
#line 672 "SB04PY.f"
	    ;
#line 672 "SB04PY.f"
	}

#line 674 "SB04PY.f"
    } else if (! notrna && ! notrnb) {

/*        Solve    A'*X*B' + ISGN*X = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        top-right corner column by column by */

/*           A(K,K)'*X(K,L)*B(L,L)' + ISGN*X(K,L) = C(K,L) - R(K,L) */

/*        where */
/*                      K-1 */
/*           R(K,L) = { SUM [A(J,K)'*X(J,L)] } * B(L,L)' + */
/*                      J=1 */
/*                       K               N */
/*                      SUM A(J,K)' * { SUM [X(J,I)*B(L,I)'] }. */
/*                      J=1            I=L+1 */

/*        Start column loop (index = L) */
/*        L1 (L2): column index of the first (last) row of X(K,L). */

#line 694 "SB04PY.f"
	lnext = *n;

#line 696 "SB04PY.f"
	for (l = *n; l >= 1; --l) {
#line 697 "SB04PY.f"
	    if (l > lnext) {
#line 697 "SB04PY.f"
		goto L180;
#line 697 "SB04PY.f"
	    }
#line 699 "SB04PY.f"
	    l2 = l;
#line 700 "SB04PY.f"
	    if (l == 1) {
#line 701 "SB04PY.f"
		l1 = l;
#line 702 "SB04PY.f"
	    } else {
#line 703 "SB04PY.f"
		if (b[l + (l - 1) * b_dim1] != 0.) {
#line 704 "SB04PY.f"
		    l1 = l - 1;
#line 705 "SB04PY.f"
		} else {
#line 706 "SB04PY.f"
		    l1 = l;
#line 707 "SB04PY.f"
		}
#line 708 "SB04PY.f"
		lnext = l1 - 1;
#line 709 "SB04PY.f"
	    }

/*           Start row loop (index = K) */
/*           K1 (K2): row index of the first (last) row of X(K,L). */

#line 714 "SB04PY.f"
	    knext = 1;

#line 716 "SB04PY.f"
	    i__1 = *m;
#line 716 "SB04PY.f"
	    for (k = 1; k <= i__1; ++k) {
#line 717 "SB04PY.f"
		if (k < knext) {
#line 717 "SB04PY.f"
		    goto L170;
#line 717 "SB04PY.f"
		}
#line 719 "SB04PY.f"
		k1 = k;
#line 720 "SB04PY.f"
		if (k == *m) {
#line 721 "SB04PY.f"
		    k2 = k;
#line 722 "SB04PY.f"
		} else {
#line 723 "SB04PY.f"
		    if (a[k + 1 + k * a_dim1] != 0.) {
#line 724 "SB04PY.f"
			k2 = k + 1;
#line 725 "SB04PY.f"
		    } else {
#line 726 "SB04PY.f"
			k2 = k;
#line 727 "SB04PY.f"
		    }
#line 728 "SB04PY.f"
		    knext = k2 + 1;
#line 729 "SB04PY.f"
		}

/* Computing MIN */
#line 731 "SB04PY.f"
		i__2 = l1 + 1;
#line 731 "SB04PY.f"
		mnl1 = min(i__2,*n);
/* Computing MIN */
#line 732 "SB04PY.f"
		i__2 = l2 + 1;
#line 732 "SB04PY.f"
		mnl2 = min(i__2,*n);
#line 733 "SB04PY.f"
		i__2 = k1 - 1;
#line 733 "SB04PY.f"
		p11 = ddot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * 
			c_dim1 + 1], &c__1);
#line 734 "SB04PY.f"
		i__2 = *n - l2;
#line 734 "SB04PY.f"
		dwork[k1] = ddot_(&i__2, &c__[k1 + mnl2 * c_dim1], ldc, &b[l1 
			+ mnl2 * b_dim1], ldb);

#line 737 "SB04PY.f"
		if (l1 == l2 && k1 == k2) {
#line 738 "SB04PY.f"
		    sumr = ddot_(&k1, &a[k1 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
#line 739 "SB04PY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1]);
#line 740 "SB04PY.f"
		    scaloc = 1.;

#line 742 "SB04PY.f"
		    a11 = a[k1 + k1 * a_dim1] * b[l1 + l1 * b_dim1] + sgn;
#line 743 "SB04PY.f"
		    da11 = abs(a11);
#line 744 "SB04PY.f"
		    if (da11 <= smin) {
#line 745 "SB04PY.f"
			a11 = smin;
#line 746 "SB04PY.f"
			da11 = smin;
#line 747 "SB04PY.f"
			*info = 1;
#line 748 "SB04PY.f"
		    }
#line 749 "SB04PY.f"
		    db = abs(vec[0]);
#line 750 "SB04PY.f"
		    if (da11 < 1. && db > 1.) {
#line 751 "SB04PY.f"
			if (db > bignum * da11) {
#line 751 "SB04PY.f"
			    scaloc = 1. / db;
#line 751 "SB04PY.f"
			}
#line 753 "SB04PY.f"
		    }
#line 754 "SB04PY.f"
		    x[0] = vec[0] * scaloc / a11;

#line 756 "SB04PY.f"
		    if (scaloc != 1.) {

#line 758 "SB04PY.f"
			i__2 = *n;
#line 758 "SB04PY.f"
			for (j = 1; j <= i__2; ++j) {
#line 759 "SB04PY.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 760 "SB04PY.f"
/* L130: */
#line 760 "SB04PY.f"
			}

#line 762 "SB04PY.f"
			dscal_(&k1, &scaloc, &dwork[1], &c__1);
#line 763 "SB04PY.f"
			*scale *= scaloc;
#line 764 "SB04PY.f"
		    }
#line 765 "SB04PY.f"
		    c__[k1 + l1 * c_dim1] = x[0];

#line 767 "SB04PY.f"
		} else if (l1 == l2 && k1 != k2) {

#line 769 "SB04PY.f"
		    i__2 = k1 - 1;
#line 769 "SB04PY.f"
		    p21 = ddot_(&i__2, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 770 "SB04PY.f"
		    i__2 = *n - l1;
#line 770 "SB04PY.f"
		    dwork[k2] = ddot_(&i__2, &c__[k2 + mnl1 * c_dim1], ldc, &
			    b[l1 + mnl1 * b_dim1], ldb);
#line 772 "SB04PY.f"
		    sumr = ddot_(&k2, &a[k1 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
#line 773 "SB04PY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1]);

#line 775 "SB04PY.f"
		    sumr = ddot_(&k2, &a[k2 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
#line 776 "SB04PY.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (sumr + p21 * b[l1 + l1 *
			     b_dim1]);

#line 778 "SB04PY.f"
		    d__1 = -sgn;
#line 778 "SB04PY.f"
		    dlaln2_(&c_true, &c__2, &c__1, &smin, &b[l1 + l1 * b_dim1]
			    , &a[k1 + k1 * a_dim1], lda, &c_b28, &c_b28, vec, 
			    &c__2, &d__1, &c_b31, x, &c__2, &scaloc, &xnorm, &
			    ierr);
#line 781 "SB04PY.f"
		    if (ierr != 0) {
#line 781 "SB04PY.f"
			*info = 1;
#line 781 "SB04PY.f"
		    }

#line 784 "SB04PY.f"
		    if (scaloc != 1.) {

#line 786 "SB04PY.f"
			i__2 = *n;
#line 786 "SB04PY.f"
			for (j = 1; j <= i__2; ++j) {
#line 787 "SB04PY.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 788 "SB04PY.f"
/* L140: */
#line 788 "SB04PY.f"
			}

#line 790 "SB04PY.f"
			dscal_(&k2, &scaloc, &dwork[1], &c__1);
#line 791 "SB04PY.f"
			*scale *= scaloc;
#line 792 "SB04PY.f"
		    }
#line 793 "SB04PY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 794 "SB04PY.f"
		    c__[k2 + l1 * c_dim1] = x[1];

#line 796 "SB04PY.f"
		} else if (l1 != l2 && k1 == k2) {

#line 798 "SB04PY.f"
		    i__2 = k1 - 1;
#line 798 "SB04PY.f"
		    p12 = ddot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
#line 799 "SB04PY.f"
		    sumr = ddot_(&k1, &a[k1 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
#line 800 "SB04PY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1] + p12 * b[l1 + l2 * b_dim1]);

#line 803 "SB04PY.f"
		    i__2 = *n - l2;
#line 803 "SB04PY.f"
		    dwork[k1 + *m] = ddot_(&i__2, &c__[k1 + mnl2 * c_dim1], 
			    ldc, &b[l2 + mnl2 * b_dim1], ldb);
#line 805 "SB04PY.f"
		    sumr = ddot_(&k1, &a[k1 * a_dim1 + 1], &c__1, &dwork[*m + 
			    1], &c__1);
#line 806 "SB04PY.f"
		    vec[1] = c__[k1 + l2 * c_dim1] - (sumr + p11 * b[l2 + l1 *
			     b_dim1] + p12 * b[l2 + l2 * b_dim1]);

#line 809 "SB04PY.f"
		    d__1 = -sgn;
#line 809 "SB04PY.f"
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &a[k1 + k1 * 
			    a_dim1], &b[l1 + l1 * b_dim1], ldb, &c_b28, &
			    c_b28, vec, &c__2, &d__1, &c_b31, x, &c__2, &
			    scaloc, &xnorm, &ierr);
#line 812 "SB04PY.f"
		    if (ierr != 0) {
#line 812 "SB04PY.f"
			*info = 1;
#line 812 "SB04PY.f"
		    }

#line 815 "SB04PY.f"
		    if (scaloc != 1.) {

#line 817 "SB04PY.f"
			i__2 = *n;
#line 817 "SB04PY.f"
			for (j = 1; j <= i__2; ++j) {
#line 818 "SB04PY.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 819 "SB04PY.f"
/* L150: */
#line 819 "SB04PY.f"
			}

#line 821 "SB04PY.f"
			dscal_(&k1, &scaloc, &dwork[1], &c__1);
#line 822 "SB04PY.f"
			dscal_(&k1, &scaloc, &dwork[*m + 1], &c__1);
#line 823 "SB04PY.f"
			*scale *= scaloc;
#line 824 "SB04PY.f"
		    }
#line 825 "SB04PY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 826 "SB04PY.f"
		    c__[k1 + l2 * c_dim1] = x[1];

#line 828 "SB04PY.f"
		} else if (l1 != l2 && k1 != k2) {

#line 830 "SB04PY.f"
		    i__2 = k1 - 1;
#line 830 "SB04PY.f"
		    p21 = ddot_(&i__2, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * 
			    c_dim1 + 1], &c__1);
#line 831 "SB04PY.f"
		    i__2 = k1 - 1;
#line 831 "SB04PY.f"
		    p12 = ddot_(&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);
#line 832 "SB04PY.f"
		    i__2 = k1 - 1;
#line 832 "SB04PY.f"
		    p22 = ddot_(&i__2, &a[k2 * a_dim1 + 1], &c__1, &c__[l2 * 
			    c_dim1 + 1], &c__1);

#line 834 "SB04PY.f"
		    i__2 = *n - l2;
#line 834 "SB04PY.f"
		    dwork[k2] = ddot_(&i__2, &c__[k2 + mnl2 * c_dim1], ldc, &
			    b[l1 + mnl2 * b_dim1], ldb);
#line 836 "SB04PY.f"
		    sumr = ddot_(&k2, &a[k1 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
#line 837 "SB04PY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1] + p12 * b[l1 + l2 * b_dim1]);

#line 840 "SB04PY.f"
		    sumr = ddot_(&k2, &a[k2 * a_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
#line 841 "SB04PY.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (sumr + p21 * b[l1 + l1 *
			     b_dim1] + p22 * b[l1 + l2 * b_dim1]);

#line 844 "SB04PY.f"
		    i__2 = *n - l2;
#line 844 "SB04PY.f"
		    dwork[k1 + *m] = ddot_(&i__2, &c__[k1 + mnl2 * c_dim1], 
			    ldc, &b[l2 + mnl2 * b_dim1], ldb);
#line 846 "SB04PY.f"
		    i__2 = *n - l2;
#line 846 "SB04PY.f"
		    dwork[k2 + *m] = ddot_(&i__2, &c__[k2 + mnl2 * c_dim1], 
			    ldc, &b[l2 + mnl2 * b_dim1], ldb);
#line 848 "SB04PY.f"
		    sumr = ddot_(&k2, &a[k1 * a_dim1 + 1], &c__1, &dwork[*m + 
			    1], &c__1);
#line 849 "SB04PY.f"
		    vec[2] = c__[k1 + l2 * c_dim1] - (sumr + p11 * b[l2 + l1 *
			     b_dim1] + p12 * b[l2 + l2 * b_dim1]);

#line 852 "SB04PY.f"
		    sumr = ddot_(&k2, &a[k2 * a_dim1 + 1], &c__1, &dwork[*m + 
			    1], &c__1);
#line 853 "SB04PY.f"
		    vec[3] = c__[k2 + l2 * c_dim1] - (sumr + p21 * b[l2 + l1 *
			     b_dim1] + p22 * b[l2 + l2 * b_dim1]);

#line 856 "SB04PY.f"
		    sb04px_(&c_true, &c_true, isgn, &c__2, &c__2, &a[k1 + k1 *
			     a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec, &
			    c__2, &scaloc, x, &c__2, &xnorm, &ierr);
#line 859 "SB04PY.f"
		    if (ierr != 0) {
#line 859 "SB04PY.f"
			*info = 1;
#line 859 "SB04PY.f"
		    }

#line 862 "SB04PY.f"
		    if (scaloc != 1.) {

#line 864 "SB04PY.f"
			i__2 = *n;
#line 864 "SB04PY.f"
			for (j = 1; j <= i__2; ++j) {
#line 865 "SB04PY.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 866 "SB04PY.f"
/* L160: */
#line 866 "SB04PY.f"
			}

#line 868 "SB04PY.f"
			dscal_(&k2, &scaloc, &dwork[1], &c__1);
#line 869 "SB04PY.f"
			dscal_(&k2, &scaloc, &dwork[*m + 1], &c__1);
#line 870 "SB04PY.f"
			*scale *= scaloc;
#line 871 "SB04PY.f"
		    }
#line 872 "SB04PY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 873 "SB04PY.f"
		    c__[k1 + l2 * c_dim1] = x[2];
#line 874 "SB04PY.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 875 "SB04PY.f"
		    c__[k2 + l2 * c_dim1] = x[3];
#line 876 "SB04PY.f"
		}

#line 878 "SB04PY.f"
L170:
#line 878 "SB04PY.f"
		;
#line 878 "SB04PY.f"
	    }

#line 880 "SB04PY.f"
L180:
#line 880 "SB04PY.f"
	    ;
#line 880 "SB04PY.f"
	}

#line 882 "SB04PY.f"
    } else {

/*        Solve    A*X*B' + ISGN*X = scale*C. */

/*        The (K,L)th block of X is determined starting from */
/*        bottom-right corner column by column by */

/*            A(K,K)*X(K,L)*B(L,L)' + ISGN*X(K,L) = C(K,L) - R(K,L) */

/*        where */
/*                       M */
/*           R(K,L) = { SUM [A(K,J)*X(J,L)] } * B(L,L)' + */
/*                     J=K+1 */
/*                       M              N */
/*                      SUM { A(K,J) * SUM [X(J,I)*B(L,I)'] }. */
/*                      J=K           I=L+1 */

/*        Start column loop (index = L) */
/*        L1 (L2): column index of the first (last) row of X(K,L). */

#line 902 "SB04PY.f"
	lnext = *n;

#line 904 "SB04PY.f"
	for (l = *n; l >= 1; --l) {
#line 905 "SB04PY.f"
	    if (l > lnext) {
#line 905 "SB04PY.f"
		goto L240;
#line 905 "SB04PY.f"
	    }
#line 907 "SB04PY.f"
	    l2 = l;
#line 908 "SB04PY.f"
	    if (l == 1) {
#line 909 "SB04PY.f"
		l1 = l;
#line 910 "SB04PY.f"
	    } else {
#line 911 "SB04PY.f"
		if (b[l + (l - 1) * b_dim1] != 0.) {
#line 912 "SB04PY.f"
		    l1 = l - 1;
#line 913 "SB04PY.f"
		} else {
#line 914 "SB04PY.f"
		    l1 = l;
#line 915 "SB04PY.f"
		}
#line 916 "SB04PY.f"
		lnext = l1 - 1;
#line 917 "SB04PY.f"
	    }

/*           Start row loop (index = K) */
/*           K1 (K2): row index of the first (last) row of X(K,L). */

#line 922 "SB04PY.f"
	    knext = *m;

#line 924 "SB04PY.f"
	    for (k = *m; k >= 1; --k) {
#line 925 "SB04PY.f"
		if (k > knext) {
#line 925 "SB04PY.f"
		    goto L230;
#line 925 "SB04PY.f"
		}
#line 927 "SB04PY.f"
		k2 = k;
#line 928 "SB04PY.f"
		if (k == 1) {
#line 929 "SB04PY.f"
		    k1 = k;
#line 930 "SB04PY.f"
		} else {
#line 931 "SB04PY.f"
		    if (a[k + (k - 1) * a_dim1] != 0.) {
#line 932 "SB04PY.f"
			k1 = k - 1;
#line 933 "SB04PY.f"
		    } else {
#line 934 "SB04PY.f"
			k1 = k;
#line 935 "SB04PY.f"
		    }
#line 936 "SB04PY.f"
		    knext = k1 - 1;
#line 937 "SB04PY.f"
		}

/* Computing MIN */
#line 939 "SB04PY.f"
		i__1 = k1 + 1;
#line 939 "SB04PY.f"
		mnk1 = min(i__1,*m);
/* Computing MIN */
#line 940 "SB04PY.f"
		i__1 = k2 + 1;
#line 940 "SB04PY.f"
		mnk2 = min(i__1,*m);
/* Computing MIN */
#line 941 "SB04PY.f"
		i__1 = l1 + 1;
#line 941 "SB04PY.f"
		mnl1 = min(i__1,*n);
/* Computing MIN */
#line 942 "SB04PY.f"
		i__1 = l2 + 1;
#line 942 "SB04PY.f"
		mnl2 = min(i__1,*n);
#line 943 "SB04PY.f"
		i__1 = *m - k2;
#line 943 "SB04PY.f"
		p11 = ddot_(&i__1, &a[k1 + mnk2 * a_dim1], lda, &c__[mnk2 + 
			l1 * c_dim1], &c__1);
#line 944 "SB04PY.f"
		i__1 = *n - l2;
#line 944 "SB04PY.f"
		dwork[k1] = ddot_(&i__1, &c__[k1 + mnl2 * c_dim1], ldc, &b[l1 
			+ mnl2 * b_dim1], ldb);

#line 947 "SB04PY.f"
		if (l1 == l2 && k1 == k2) {

#line 949 "SB04PY.f"
		    i__1 = *m - k1 + 1;
#line 949 "SB04PY.f"
		    sumr = ddot_(&i__1, &a[k1 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
#line 951 "SB04PY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1]);
#line 952 "SB04PY.f"
		    scaloc = 1.;

#line 954 "SB04PY.f"
		    a11 = a[k1 + k1 * a_dim1] * b[l1 + l1 * b_dim1] + sgn;
#line 955 "SB04PY.f"
		    da11 = abs(a11);
#line 956 "SB04PY.f"
		    if (da11 <= smin) {
#line 957 "SB04PY.f"
			a11 = smin;
#line 958 "SB04PY.f"
			da11 = smin;
#line 959 "SB04PY.f"
			*info = 1;
#line 960 "SB04PY.f"
		    }
#line 961 "SB04PY.f"
		    db = abs(vec[0]);
#line 962 "SB04PY.f"
		    if (da11 < 1. && db > 1.) {
#line 963 "SB04PY.f"
			if (db > bignum * da11) {
#line 963 "SB04PY.f"
			    scaloc = 1. / db;
#line 963 "SB04PY.f"
			}
#line 965 "SB04PY.f"
		    }
#line 966 "SB04PY.f"
		    x[0] = vec[0] * scaloc / a11;

#line 968 "SB04PY.f"
		    if (scaloc != 1.) {

#line 970 "SB04PY.f"
			i__1 = *n;
#line 970 "SB04PY.f"
			for (j = 1; j <= i__1; ++j) {
#line 971 "SB04PY.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 972 "SB04PY.f"
/* L190: */
#line 972 "SB04PY.f"
			}

#line 974 "SB04PY.f"
			i__1 = *m - k1 + 1;
#line 974 "SB04PY.f"
			dscal_(&i__1, &scaloc, &dwork[k1], &c__1);
#line 975 "SB04PY.f"
			*scale *= scaloc;
#line 976 "SB04PY.f"
		    }
#line 977 "SB04PY.f"
		    c__[k1 + l1 * c_dim1] = x[0];

#line 979 "SB04PY.f"
		} else if (l1 == l2 && k1 != k2) {

#line 981 "SB04PY.f"
		    i__1 = *m - k2;
#line 981 "SB04PY.f"
		    p21 = ddot_(&i__1, &a[k2 + mnk2 * a_dim1], lda, &c__[mnk2 
			    + l1 * c_dim1], &c__1);
#line 983 "SB04PY.f"
		    i__1 = *n - l1;
#line 983 "SB04PY.f"
		    dwork[k2] = ddot_(&i__1, &c__[k2 + mnl1 * c_dim1], ldc, &
			    b[l1 + mnl1 * b_dim1], ldb);
#line 985 "SB04PY.f"
		    i__1 = *m - k1 + 1;
#line 985 "SB04PY.f"
		    sumr = ddot_(&i__1, &a[k1 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
#line 987 "SB04PY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1]);

#line 989 "SB04PY.f"
		    i__1 = *m - k1 + 1;
#line 989 "SB04PY.f"
		    sumr = ddot_(&i__1, &a[k2 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
#line 991 "SB04PY.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (sumr + p21 * b[l1 + l1 *
			     b_dim1]);

#line 993 "SB04PY.f"
		    d__1 = -sgn;
#line 993 "SB04PY.f"
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &b[l1 + l1 * 
			    b_dim1], &a[k1 + k1 * a_dim1], lda, &c_b28, &
			    c_b28, vec, &c__2, &d__1, &c_b31, x, &c__2, &
			    scaloc, &xnorm, &ierr);
#line 996 "SB04PY.f"
		    if (ierr != 0) {
#line 996 "SB04PY.f"
			*info = 1;
#line 996 "SB04PY.f"
		    }

#line 999 "SB04PY.f"
		    if (scaloc != 1.) {

#line 1001 "SB04PY.f"
			i__1 = *n;
#line 1001 "SB04PY.f"
			for (j = 1; j <= i__1; ++j) {
#line 1002 "SB04PY.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 1003 "SB04PY.f"
/* L200: */
#line 1003 "SB04PY.f"
			}

#line 1005 "SB04PY.f"
			i__1 = *m - k1 + 1;
#line 1005 "SB04PY.f"
			dscal_(&i__1, &scaloc, &dwork[k1], &c__1);
#line 1006 "SB04PY.f"
			*scale *= scaloc;
#line 1007 "SB04PY.f"
		    }
#line 1008 "SB04PY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 1009 "SB04PY.f"
		    c__[k2 + l1 * c_dim1] = x[1];

#line 1011 "SB04PY.f"
		} else if (l1 != l2 && k1 == k2) {

#line 1013 "SB04PY.f"
		    i__1 = *m - k1;
#line 1013 "SB04PY.f"
		    p12 = ddot_(&i__1, &a[k1 + mnk1 * a_dim1], lda, &c__[mnk1 
			    + l2 * c_dim1], &c__1);
#line 1015 "SB04PY.f"
		    i__1 = *m - k1 + 1;
#line 1015 "SB04PY.f"
		    sumr = ddot_(&i__1, &a[k1 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
#line 1017 "SB04PY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1] + p12 * b[l1 + l2 * b_dim1]);

#line 1020 "SB04PY.f"
		    i__1 = *n - l2;
#line 1020 "SB04PY.f"
		    dwork[k1 + *m] = ddot_(&i__1, &c__[k1 + mnl2 * c_dim1], 
			    ldc, &b[l2 + mnl2 * b_dim1], ldb);
#line 1022 "SB04PY.f"
		    i__1 = *m - k1 + 1;
#line 1022 "SB04PY.f"
		    sumr = ddot_(&i__1, &a[k1 + k1 * a_dim1], lda, &dwork[k1 
			    + *m], &c__1);
#line 1024 "SB04PY.f"
		    vec[1] = c__[k1 + l2 * c_dim1] - (sumr + p11 * b[l2 + l1 *
			     b_dim1] + p12 * b[l2 + l2 * b_dim1]);

#line 1027 "SB04PY.f"
		    d__1 = -sgn;
#line 1027 "SB04PY.f"
		    dlaln2_(&c_false, &c__2, &c__1, &smin, &a[k1 + k1 * 
			    a_dim1], &b[l1 + l1 * b_dim1], ldb, &c_b28, &
			    c_b28, vec, &c__2, &d__1, &c_b31, x, &c__2, &
			    scaloc, &xnorm, &ierr);
#line 1030 "SB04PY.f"
		    if (ierr != 0) {
#line 1030 "SB04PY.f"
			*info = 1;
#line 1030 "SB04PY.f"
		    }

#line 1033 "SB04PY.f"
		    if (scaloc != 1.) {

#line 1035 "SB04PY.f"
			i__1 = *n;
#line 1035 "SB04PY.f"
			for (j = 1; j <= i__1; ++j) {
#line 1036 "SB04PY.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 1037 "SB04PY.f"
/* L210: */
#line 1037 "SB04PY.f"
			}

#line 1039 "SB04PY.f"
			i__1 = *m - k1 + 1;
#line 1039 "SB04PY.f"
			dscal_(&i__1, &scaloc, &dwork[k1], &c__1);
#line 1040 "SB04PY.f"
			i__1 = *m - k1 + 1;
#line 1040 "SB04PY.f"
			dscal_(&i__1, &scaloc, &dwork[k1 + *m], &c__1);
#line 1041 "SB04PY.f"
			*scale *= scaloc;
#line 1042 "SB04PY.f"
		    }
#line 1043 "SB04PY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 1044 "SB04PY.f"
		    c__[k1 + l2 * c_dim1] = x[1];

#line 1046 "SB04PY.f"
		} else if (l1 != l2 && k1 != k2) {

#line 1048 "SB04PY.f"
		    i__1 = *m - k2;
#line 1048 "SB04PY.f"
		    p21 = ddot_(&i__1, &a[k2 + mnk2 * a_dim1], lda, &c__[mnk2 
			    + l1 * c_dim1], &c__1);
#line 1050 "SB04PY.f"
		    i__1 = *m - k2;
#line 1050 "SB04PY.f"
		    p12 = ddot_(&i__1, &a[k1 + mnk2 * a_dim1], lda, &c__[mnk2 
			    + l2 * c_dim1], &c__1);
#line 1052 "SB04PY.f"
		    i__1 = *m - k2;
#line 1052 "SB04PY.f"
		    p22 = ddot_(&i__1, &a[k2 + mnk2 * a_dim1], lda, &c__[mnk2 
			    + l2 * c_dim1], &c__1);

#line 1055 "SB04PY.f"
		    i__1 = *n - l2;
#line 1055 "SB04PY.f"
		    dwork[k2] = ddot_(&i__1, &c__[k2 + mnl2 * c_dim1], ldc, &
			    b[l1 + mnl2 * b_dim1], ldb);
#line 1057 "SB04PY.f"
		    i__1 = *m - k1 + 1;
#line 1057 "SB04PY.f"
		    sumr = ddot_(&i__1, &a[k1 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
#line 1059 "SB04PY.f"
		    vec[0] = c__[k1 + l1 * c_dim1] - (sumr + p11 * b[l1 + l1 *
			     b_dim1] + p12 * b[l1 + l2 * b_dim1]);

#line 1062 "SB04PY.f"
		    i__1 = *m - k1 + 1;
#line 1062 "SB04PY.f"
		    sumr = ddot_(&i__1, &a[k2 + k1 * a_dim1], lda, &dwork[k1],
			     &c__1);
#line 1064 "SB04PY.f"
		    vec[1] = c__[k2 + l1 * c_dim1] - (sumr + p21 * b[l1 + l1 *
			     b_dim1] + p22 * b[l1 + l2 * b_dim1]);

#line 1067 "SB04PY.f"
		    i__1 = *n - l2;
#line 1067 "SB04PY.f"
		    dwork[k1 + *m] = ddot_(&i__1, &c__[k1 + mnl2 * c_dim1], 
			    ldc, &b[l2 + mnl2 * b_dim1], ldb);
#line 1069 "SB04PY.f"
		    i__1 = *n - l2;
#line 1069 "SB04PY.f"
		    dwork[k2 + *m] = ddot_(&i__1, &c__[k2 + mnl2 * c_dim1], 
			    ldc, &b[l2 + mnl2 * b_dim1], ldb);
#line 1071 "SB04PY.f"
		    i__1 = *m - k1 + 1;
#line 1071 "SB04PY.f"
		    sumr = ddot_(&i__1, &a[k1 + k1 * a_dim1], lda, &dwork[k1 
			    + *m], &c__1);
#line 1073 "SB04PY.f"
		    vec[2] = c__[k1 + l2 * c_dim1] - (sumr + p11 * b[l2 + l1 *
			     b_dim1] + p12 * b[l2 + l2 * b_dim1]);

#line 1076 "SB04PY.f"
		    i__1 = *m - k1 + 1;
#line 1076 "SB04PY.f"
		    sumr = ddot_(&i__1, &a[k2 + k1 * a_dim1], lda, &dwork[k1 
			    + *m], &c__1);
#line 1078 "SB04PY.f"
		    vec[3] = c__[k2 + l2 * c_dim1] - (sumr + p21 * b[l2 + l1 *
			     b_dim1] + p22 * b[l2 + l2 * b_dim1]);

#line 1081 "SB04PY.f"
		    sb04px_(&c_false, &c_true, isgn, &c__2, &c__2, &a[k1 + k1 
			    * a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec, &
			    c__2, &scaloc, x, &c__2, &xnorm, &ierr);
#line 1084 "SB04PY.f"
		    if (ierr != 0) {
#line 1084 "SB04PY.f"
			*info = 1;
#line 1084 "SB04PY.f"
		    }

#line 1087 "SB04PY.f"
		    if (scaloc != 1.) {

#line 1089 "SB04PY.f"
			i__1 = *n;
#line 1089 "SB04PY.f"
			for (j = 1; j <= i__1; ++j) {
#line 1090 "SB04PY.f"
			    dscal_(m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 1091 "SB04PY.f"
/* L220: */
#line 1091 "SB04PY.f"
			}

#line 1093 "SB04PY.f"
			i__1 = *m - k1 + 1;
#line 1093 "SB04PY.f"
			dscal_(&i__1, &scaloc, &dwork[k1], &c__1);
#line 1094 "SB04PY.f"
			i__1 = *m - k1 + 1;
#line 1094 "SB04PY.f"
			dscal_(&i__1, &scaloc, &dwork[k1 + *m], &c__1);
#line 1095 "SB04PY.f"
			*scale *= scaloc;
#line 1096 "SB04PY.f"
		    }
#line 1097 "SB04PY.f"
		    c__[k1 + l1 * c_dim1] = x[0];
#line 1098 "SB04PY.f"
		    c__[k1 + l2 * c_dim1] = x[2];
#line 1099 "SB04PY.f"
		    c__[k2 + l1 * c_dim1] = x[1];
#line 1100 "SB04PY.f"
		    c__[k2 + l2 * c_dim1] = x[3];
#line 1101 "SB04PY.f"
		}

#line 1103 "SB04PY.f"
L230:
#line 1103 "SB04PY.f"
		;
#line 1103 "SB04PY.f"
	    }

#line 1105 "SB04PY.f"
L240:
#line 1105 "SB04PY.f"
	    ;
#line 1105 "SB04PY.f"
	}

#line 1107 "SB04PY.f"
    }

#line 1109 "SB04PY.f"
    return 0;
/* *** Last line of SB04PY *** */
} /* sb04py_ */

