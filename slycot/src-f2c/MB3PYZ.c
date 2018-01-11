#line 1 "MB3PYZ.f"
/* MB3PYZ.f -- translated by f2c (version 20100827).
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

#line 1 "MB3PYZ.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int mb3pyz_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *rcond, doublereal *svlmax, integer *rank, 
	doublereal *sval, integer *jpvt, doublecomplex *tau, doublereal *
	dwork, doublecomplex *zwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublecomplex c1, c2, s1, s2, aii;
    static integer mki, nki, pvt;
    static doublereal smin, temp, smax, temp2;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static integer itemp, ismin;
    extern /* Subroutine */ int zlarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen);
    static integer ismax, jwork;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), zlaic1_(integer *, 
	    integer *, doublecomplex *, doublereal *, doublecomplex *, 
	    doublecomplex *, doublereal *, doublecomplex *, doublecomplex *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zlarfg_(
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *), zlacgv_(integer *, doublecomplex *, integer *);
    static doublereal sminpr, smaxpr;


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

/*     To compute a rank-revealing RQ factorization of a complex general */
/*     M-by-N matrix  A,  which may be rank-deficient, and estimate its */
/*     effective rank using incremental condition estimation. */

/*     The routine uses a truncated RQ factorization with row pivoting: */
/*                                   [ R11 R12 ] */
/*        P * A = R * Q,  where  R = [         ], */
/*                                   [  0  R22 ] */
/*     with R22 defined as the largest trailing upper triangular */
/*     submatrix whose estimated condition number is less than 1/RCOND. */
/*     The order of R22, RANK, is the effective rank of A.  Condition */
/*     estimation is performed during the RQ factorization process. */
/*     Matrix R11 is full (but of small norm), or empty. */

/*     MB3PYZ  does not perform any scaling of the matrix A. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     A       (input/output) COMPLEX*16 array, dimension ( LDA, N ) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the given matrix A. */
/*             On exit, the upper triangle of the subarray */
/*             A(M-RANK+1:M,N-RANK+1:N) contains the RANK-by-RANK upper */
/*             triangular matrix R22;  the remaining elements in the last */
/*             RANK  rows, with the array TAU, represent the unitary */
/*             matrix Q as a product of  RANK  elementary reflectors */
/*             (see METHOD).  The first  M-RANK  rows contain the result */
/*             of the RQ factorization process used. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     RCOND   (input) DOUBLE PRECISION */
/*             RCOND is used to determine the effective rank of A, which */
/*             is defined as the order of the largest trailing triangular */
/*             submatrix R22 in the RQ factorization with pivoting of A, */
/*             whose estimated condition number is less than 1/RCOND. */
/*             0 <= RCOND <= 1. */
/*             NOTE that when SVLMAX > 0, the estimated rank could be */
/*             less than that defined above (see SVLMAX). */

/*     SVLMAX  (input) DOUBLE PRECISION */
/*             If A is a submatrix of another matrix B, and the rank */
/*             decision should be related to that matrix, then SVLMAX */
/*             should be an estimate of the largest singular value of B */
/*             (for instance, the Frobenius norm of B).  If this is not */
/*             the case, the input value SVLMAX = 0 should work. */
/*             SVLMAX >= 0. */

/*     RANK    (output) INTEGER */
/*             The effective (estimated) rank of A, i.e., the order of */
/*             the submatrix R22. */

/*     SVAL    (output) DOUBLE PRECISION array, dimension ( 3 ) */
/*             The estimates of some of the singular values of the */
/*             triangular factor R: */
/*             SVAL(1): largest singular value of */
/*                      R(M-RANK+1:M,N-RANK+1:N); */
/*             SVAL(2): smallest singular value of */
/*                      R(M-RANK+1:M,N-RANK+1:N); */
/*             SVAL(3): smallest singular value of R(M-RANK:M,N-RANK:N), */
/*                      if RANK < MIN( M, N ), or of */
/*                      R(M-RANK+1:M,N-RANK+1:N), otherwise. */
/*             If the triangular factorization is a rank-revealing one */
/*             (which will be the case if the trailing rows were well- */
/*             conditioned), then SVAL(1) will also be an estimate for */
/*             the largest singular value of A, and SVAL(2) and SVAL(3) */
/*             will be estimates for the RANK-th and (RANK+1)-st singular */
/*             values of A, respectively. */
/*             By examining these values, one can confirm that the rank */
/*             is well defined with respect to the chosen value of RCOND. */
/*             The ratio SVAL(1)/SVAL(2) is an estimate of the condition */
/*             number of R(M-RANK+1:M,N-RANK+1:N). */

/*     JPVT    (output) INTEGER array, dimension ( M ) */
/*             If JPVT(i) = k, then the i-th row of P*A was the k-th row */
/*             of A. */

/*     TAU     (output) COMPLEX*16 array, dimension ( MIN( M, N ) ) */
/*             The trailing  RANK  elements of TAU contain the scalar */
/*             factors of the elementary reflectors. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension ( 2*M ) */

/*     ZWORK   COMPLEX*16 array, dimension ( 3*M-1 ) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine computes a truncated RQ factorization with row */
/*     pivoting of A,  P * A = R * Q,  with  R  defined above, and, */
/*     during this process, finds the largest trailing submatrix whose */
/*     estimated condition number is less than 1/RCOND, taking the */
/*     possible positive value of SVLMAX into account.  This is performed */
/*     using an adaptation of the LAPACK incremental condition estimation */
/*     scheme and a slightly modified rank decision test.  The */
/*     factorization process stops when  RANK  has been determined. */

/*     The matrix Q is represented as a product of elementary reflectors */

/*        Q = H(k-rank+1)' H(k-rank+2)' . . . H(k)', where k = min(m,n). */

/*     Each H(i) has the form */

/*        H = I - tau * v * v' */

/*     where tau is a complex scalar, and v is a complex vector with */
/*     v(n-k+i+1:n) = 0 and v(n-k+i) = 1; conjg(v(1:n-k+i-1)) is stored */
/*     on exit in A(m-k+i,1:n-k+i-1), and tau in TAU(i). */

/*     The matrix P is represented in jpvt as follows: If */
/*        jpvt(j) = i */
/*     then the jth row of P is the ith canonical unit vector. */

/*     REFERENCES */

/*     [1] Bischof, C.H. and P. Tang. */
/*         Generalizing Incremental Condition Estimation. */
/*         LAPACK Working Notes 32, Mathematics and Computer Science */
/*         Division, Argonne National Laboratory, UT, CS-91-132, */
/*         May 1991. */

/*     [2] Bischof, C.H. and P. Tang. */
/*         Robust Incremental Condition Estimation. */
/*         LAPACK Working Notes 33, Mathematics and Computer Science */
/*         Division, Argonne National Laboratory, UT, CS-91-133, */
/*         May 1991. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1998. */
/*     Complex version: V. Sima, Research Institute for Informatics, */
/*     Bucharest, Nov. 2008. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Eigenvalue problem, matrix operations, unitary transformation, */
/*     singular values. */

/*    ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
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

#line 225 "MB3PYZ.f"
    /* Parameter adjustments */
#line 225 "MB3PYZ.f"
    a_dim1 = *lda;
#line 225 "MB3PYZ.f"
    a_offset = 1 + a_dim1;
#line 225 "MB3PYZ.f"
    a -= a_offset;
#line 225 "MB3PYZ.f"
    --sval;
#line 225 "MB3PYZ.f"
    --jpvt;
#line 225 "MB3PYZ.f"
    --tau;
#line 225 "MB3PYZ.f"
    --dwork;
#line 225 "MB3PYZ.f"
    --zwork;
#line 225 "MB3PYZ.f"

#line 225 "MB3PYZ.f"
    /* Function Body */
#line 225 "MB3PYZ.f"
    *info = 0;
#line 226 "MB3PYZ.f"
    if (*m < 0) {
#line 227 "MB3PYZ.f"
	*info = -1;
#line 228 "MB3PYZ.f"
    } else if (*n < 0) {
#line 229 "MB3PYZ.f"
	*info = -2;
#line 230 "MB3PYZ.f"
    } else if (*lda < max(1,*m)) {
#line 231 "MB3PYZ.f"
	*info = -4;
#line 232 "MB3PYZ.f"
    } else if (*rcond < 0. || *rcond > 1.) {
#line 233 "MB3PYZ.f"
	*info = -5;
#line 234 "MB3PYZ.f"
    } else if (*svlmax < 0.) {
#line 235 "MB3PYZ.f"
	*info = -6;
#line 236 "MB3PYZ.f"
    }

#line 238 "MB3PYZ.f"
    if (*info != 0) {
#line 239 "MB3PYZ.f"
	i__1 = -(*info);
#line 239 "MB3PYZ.f"
	xerbla_("MB3PYZ", &i__1, (ftnlen)6);
#line 240 "MB3PYZ.f"
	return 0;
#line 241 "MB3PYZ.f"
    }

/*     Quick return if possible. */

#line 245 "MB3PYZ.f"
    k = min(*m,*n);
#line 246 "MB3PYZ.f"
    if (k == 0) {
#line 247 "MB3PYZ.f"
	*rank = 0;
#line 248 "MB3PYZ.f"
	sval[1] = 0.;
#line 249 "MB3PYZ.f"
	sval[2] = 0.;
#line 250 "MB3PYZ.f"
	sval[3] = 0.;
#line 251 "MB3PYZ.f"
	return 0;
#line 252 "MB3PYZ.f"
    }

#line 254 "MB3PYZ.f"
    ismin = 1;
#line 255 "MB3PYZ.f"
    ismax = ismin + *m;
#line 256 "MB3PYZ.f"
    jwork = ismax + *m;

/*     Initialize partial row norms and pivoting vector. The first m */
/*     elements of DWORK store the exact row norms. */

#line 261 "MB3PYZ.f"
    i__1 = *m;
#line 261 "MB3PYZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 262 "MB3PYZ.f"
	dwork[i__] = dznrm2_(n, &a[i__ + a_dim1], lda);
#line 263 "MB3PYZ.f"
	dwork[*m + i__] = dwork[i__];
#line 264 "MB3PYZ.f"
	jpvt[i__] = i__;
#line 265 "MB3PYZ.f"
/* L10: */
#line 265 "MB3PYZ.f"
    }

/*     Compute factorization and determine RANK using incremental */
/*     condition estimation. */

#line 270 "MB3PYZ.f"
    *rank = 0;

#line 272 "MB3PYZ.f"
L20:
#line 273 "MB3PYZ.f"
    if (*rank < k) {
#line 274 "MB3PYZ.f"
	i__ = k - *rank;

/*        Determine ith pivot row and swap if necessary. */

#line 278 "MB3PYZ.f"
	mki = *m - *rank;
#line 279 "MB3PYZ.f"
	nki = *n - *rank;
#line 280 "MB3PYZ.f"
	pvt = idamax_(&mki, &dwork[1], &c__1);

#line 282 "MB3PYZ.f"
	if (pvt != mki) {
#line 283 "MB3PYZ.f"
	    zswap_(n, &a[pvt + a_dim1], lda, &a[mki + a_dim1], lda);
#line 284 "MB3PYZ.f"
	    itemp = jpvt[pvt];
#line 285 "MB3PYZ.f"
	    jpvt[pvt] = jpvt[mki];
#line 286 "MB3PYZ.f"
	    jpvt[mki] = itemp;
#line 287 "MB3PYZ.f"
	    dwork[pvt] = dwork[mki];
#line 288 "MB3PYZ.f"
	    dwork[*m + pvt] = dwork[*m + mki];
#line 289 "MB3PYZ.f"
	}

#line 291 "MB3PYZ.f"
	if (nki > 1) {

/*           Save A(m-k+i,n-k+i) and generate elementary reflector H(i) */
/*           to annihilate A(m-k+i,1:n-k+i-1), k = min(m,n). */
/*               A(m-k+i,1:n-k+i) * H(tau,v)        = [0 , *]         <=> */
/*               H(conj(tau),v) A(m-k+i,1:n-k+i)^H  = [0 ; *], */
/*           using H(tau,v)^H = H(conj(tau),v). */

#line 299 "MB3PYZ.f"
	    zlacgv_(&nki, &a[mki + a_dim1], lda);
#line 300 "MB3PYZ.f"
	    i__1 = mki + nki * a_dim1;
#line 300 "MB3PYZ.f"
	    aii.r = a[i__1].r, aii.i = a[i__1].i;
#line 301 "MB3PYZ.f"
	    zlarfg_(&nki, &a[mki + nki * a_dim1], &a[mki + a_dim1], lda, &tau[
		    i__]);
#line 303 "MB3PYZ.f"
	}

#line 305 "MB3PYZ.f"
	if (*rank == 0) {

/*           Initialize; exit if matrix is zero (RANK = 0). */

#line 309 "MB3PYZ.f"
	    smax = z_abs(&a[*m + *n * a_dim1]);
#line 310 "MB3PYZ.f"
	    if (smax == 0.) {
#line 311 "MB3PYZ.f"
		sval[1] = 0.;
#line 312 "MB3PYZ.f"
		sval[2] = 0.;
#line 313 "MB3PYZ.f"
		sval[3] = 0.;
#line 314 "MB3PYZ.f"
		return 0;
#line 315 "MB3PYZ.f"
	    }
#line 316 "MB3PYZ.f"
	    smin = smax;
#line 317 "MB3PYZ.f"
	    smaxpr = smax;
#line 318 "MB3PYZ.f"
	    sminpr = smin;
#line 319 "MB3PYZ.f"
	    c1.r = 1., c1.i = 0.;
#line 320 "MB3PYZ.f"
	    c2.r = 1., c2.i = 0.;
#line 321 "MB3PYZ.f"
	} else {

/*           One step of incremental condition estimation. */

#line 325 "MB3PYZ.f"
	    zcopy_(rank, &a[mki + (nki + 1) * a_dim1], lda, &zwork[jwork], &
		    c__1);
#line 326 "MB3PYZ.f"
	    zlaic1_(&c__2, rank, &zwork[ismin], &smin, &zwork[jwork], &a[mki 
		    + nki * a_dim1], &sminpr, &s1, &c1);
#line 328 "MB3PYZ.f"
	    zlaic1_(&c__1, rank, &zwork[ismax], &smax, &zwork[jwork], &a[mki 
		    + nki * a_dim1], &smaxpr, &s2, &c2);
#line 330 "MB3PYZ.f"
	}

#line 332 "MB3PYZ.f"
	if (*svlmax * *rcond <= smaxpr) {
#line 333 "MB3PYZ.f"
	    if (*svlmax * *rcond <= sminpr) {
#line 334 "MB3PYZ.f"
		if (smaxpr * *rcond <= sminpr) {

#line 336 "MB3PYZ.f"
		    if (mki > 1) {

/*                    Continue factorization, as rank is at least RANK. */
/*                    Apply H(i) to A(1:m-k+i-1,1:n-k+i) from the right. */

#line 341 "MB3PYZ.f"
			i__1 = mki + nki * a_dim1;
#line 341 "MB3PYZ.f"
			aii.r = a[i__1].r, aii.i = a[i__1].i;
#line 342 "MB3PYZ.f"
			i__1 = mki + nki * a_dim1;
#line 342 "MB3PYZ.f"
			a[i__1].r = 1., a[i__1].i = 0.;
#line 343 "MB3PYZ.f"
			i__1 = mki - 1;
#line 343 "MB3PYZ.f"
			zlarf_("Right", &i__1, &nki, &a[mki + a_dim1], lda, &
				tau[i__], &a[a_offset], lda, &zwork[jwork], (
				ftnlen)5);
#line 345 "MB3PYZ.f"
			i__1 = mki + nki * a_dim1;
#line 345 "MB3PYZ.f"
			a[i__1].r = aii.r, a[i__1].i = aii.i;

/*                    Update partial row norms. */

#line 349 "MB3PYZ.f"
			i__1 = mki - 1;
#line 349 "MB3PYZ.f"
			for (j = 1; j <= i__1; ++j) {
#line 350 "MB3PYZ.f"
			    if (dwork[j] != 0.) {
/* Computing 2nd power */
#line 351 "MB3PYZ.f"
				d__1 = z_abs(&a[j + nki * a_dim1]) / dwork[j];
#line 351 "MB3PYZ.f"
				temp = 1. - d__1 * d__1;
#line 353 "MB3PYZ.f"
				temp = max(temp,0.);
/* Computing 2nd power */
#line 354 "MB3PYZ.f"
				d__1 = dwork[j] / dwork[*m + j];
#line 354 "MB3PYZ.f"
				temp2 = temp * .05 * (d__1 * d__1) + 1.;
#line 356 "MB3PYZ.f"
				if (temp2 == 1.) {
#line 357 "MB3PYZ.f"
				    i__2 = nki - 1;
#line 357 "MB3PYZ.f"
				    dwork[j] = dznrm2_(&i__2, &a[j + a_dim1], 
					    lda);
#line 359 "MB3PYZ.f"
				    dwork[*m + j] = dwork[j];
#line 360 "MB3PYZ.f"
				} else {
#line 361 "MB3PYZ.f"
				    dwork[j] *= sqrt(temp);
#line 362 "MB3PYZ.f"
				}
#line 363 "MB3PYZ.f"
			    }
#line 364 "MB3PYZ.f"
/* L30: */
#line 364 "MB3PYZ.f"
			}

#line 366 "MB3PYZ.f"
		    }

#line 368 "MB3PYZ.f"
		    i__1 = *rank;
#line 368 "MB3PYZ.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 369 "MB3PYZ.f"
			i__2 = ismin + i__ - 1;
#line 369 "MB3PYZ.f"
			i__3 = ismin + i__ - 1;
#line 369 "MB3PYZ.f"
			z__1.r = s1.r * zwork[i__3].r - s1.i * zwork[i__3].i, 
				z__1.i = s1.r * zwork[i__3].i + s1.i * zwork[
				i__3].r;
#line 369 "MB3PYZ.f"
			zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
#line 370 "MB3PYZ.f"
			i__2 = ismax + i__ - 1;
#line 370 "MB3PYZ.f"
			i__3 = ismax + i__ - 1;
#line 370 "MB3PYZ.f"
			z__1.r = s2.r * zwork[i__3].r - s2.i * zwork[i__3].i, 
				z__1.i = s2.r * zwork[i__3].i + s2.i * zwork[
				i__3].r;
#line 370 "MB3PYZ.f"
			zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
#line 371 "MB3PYZ.f"
/* L40: */
#line 371 "MB3PYZ.f"
		    }

#line 373 "MB3PYZ.f"
		    i__1 = ismin + *rank;
#line 373 "MB3PYZ.f"
		    zwork[i__1].r = c1.r, zwork[i__1].i = c1.i;
#line 374 "MB3PYZ.f"
		    i__1 = ismax + *rank;
#line 374 "MB3PYZ.f"
		    zwork[i__1].r = c2.r, zwork[i__1].i = c2.i;
#line 375 "MB3PYZ.f"
		    smin = sminpr;
#line 376 "MB3PYZ.f"
		    smax = smaxpr;
#line 377 "MB3PYZ.f"
		    ++(*rank);
#line 378 "MB3PYZ.f"
		    i__1 = nki - 1;
#line 378 "MB3PYZ.f"
		    zlacgv_(&i__1, &a[mki + a_dim1], lda);
#line 379 "MB3PYZ.f"
		    goto L20;
#line 380 "MB3PYZ.f"
		}
#line 381 "MB3PYZ.f"
	    }
#line 382 "MB3PYZ.f"
	}
#line 383 "MB3PYZ.f"
    }

/*     Restore the changed part of the (M-RANK)-th row and set SVAL. */

#line 387 "MB3PYZ.f"
    if (*rank < k && nki > 1) {
#line 388 "MB3PYZ.f"
	i__1 = nki - 1;
#line 388 "MB3PYZ.f"
	zlacgv_(&i__1, &a[mki + a_dim1], lda);
#line 389 "MB3PYZ.f"
	i__1 = nki - 1;
#line 389 "MB3PYZ.f"
	i__2 = mki + nki * a_dim1;
#line 389 "MB3PYZ.f"
	z__2.r = -a[i__2].r, z__2.i = -a[i__2].i;
#line 389 "MB3PYZ.f"
	i__3 = i__;
#line 389 "MB3PYZ.f"
	z__1.r = z__2.r * tau[i__3].r - z__2.i * tau[i__3].i, z__1.i = z__2.r 
		* tau[i__3].i + z__2.i * tau[i__3].r;
#line 389 "MB3PYZ.f"
	zscal_(&i__1, &z__1, &a[mki + a_dim1], lda);
#line 390 "MB3PYZ.f"
	i__1 = mki + nki * a_dim1;
#line 390 "MB3PYZ.f"
	a[i__1].r = aii.r, a[i__1].i = aii.i;
#line 391 "MB3PYZ.f"
    }
#line 392 "MB3PYZ.f"
    sval[1] = smax;
#line 393 "MB3PYZ.f"
    sval[2] = smin;
#line 394 "MB3PYZ.f"
    sval[3] = sminpr;

#line 396 "MB3PYZ.f"
    return 0;
/* *** Last line of MB3PYZ *** */
} /* mb3pyz_ */

