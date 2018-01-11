#line 1 "MB03PY.f"
/* MB03PY.f -- translated by f2c (version 20100827).
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

#line 1 "MB03PY.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int mb03py_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *rcond, doublereal *svlmax, integer *rank, doublereal 
	*sval, integer *jpvt, doublereal *tau, doublereal *dwork, integer *
	info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal c1, c2, s1, s2, aii;
    static integer mki, nki, pvt;
    static doublereal smin, temp, smax;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal temp2;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dlarf_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    ftnlen);
    static integer itemp, ismin;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer ismax;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer jwork;
    extern /* Subroutine */ int dlaic1_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), dlarfg_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
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

/*     To compute a rank-revealing RQ factorization of a real general */
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

/*     MB03PY  does not perform any scaling of the matrix A. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension */
/*             ( LDA, N ) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the given matrix A. */
/*             On exit, the upper triangle of the subarray */
/*             A(M-RANK+1:M,N-RANK+1:N) contains the RANK-by-RANK upper */
/*             triangular matrix R22;  the remaining elements in the last */
/*             RANK  rows, with the array TAU, represent the orthogonal */
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

/*     TAU     (output) DOUBLE PRECISION array, dimension ( MIN( M, N ) ) */
/*             The trailing  RANK  elements of TAU contain the scalar */
/*             factors of the elementary reflectors. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension ( 3*M-1 ) */

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

/*        Q = H(k-rank+1) H(k-rank+2) . . . H(k), where k = min(m,n). */

/*     Each H(i) has the form */

/*        H = I - tau * v * v' */

/*     where tau is a real scalar, and v is a real vector with */
/*     v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit */
/*     in A(m-k+i,1:n-k+i-1), and tau in TAU(i). */

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

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2001, */
/*     Jan. 2009. */

/*     KEYWORDS */

/*     Eigenvalue problem, matrix operations, orthogonal transformation, */
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

#line 220 "MB03PY.f"
    /* Parameter adjustments */
#line 220 "MB03PY.f"
    a_dim1 = *lda;
#line 220 "MB03PY.f"
    a_offset = 1 + a_dim1;
#line 220 "MB03PY.f"
    a -= a_offset;
#line 220 "MB03PY.f"
    --sval;
#line 220 "MB03PY.f"
    --jpvt;
#line 220 "MB03PY.f"
    --tau;
#line 220 "MB03PY.f"
    --dwork;
#line 220 "MB03PY.f"

#line 220 "MB03PY.f"
    /* Function Body */
#line 220 "MB03PY.f"
    *info = 0;
#line 221 "MB03PY.f"
    if (*m < 0) {
#line 222 "MB03PY.f"
	*info = -1;
#line 223 "MB03PY.f"
    } else if (*n < 0) {
#line 224 "MB03PY.f"
	*info = -2;
#line 225 "MB03PY.f"
    } else if (*lda < max(1,*m)) {
#line 226 "MB03PY.f"
	*info = -4;
#line 227 "MB03PY.f"
    } else if (*rcond < 0. || *rcond > 1.) {
#line 228 "MB03PY.f"
	*info = -5;
#line 229 "MB03PY.f"
    } else if (*svlmax < 0.) {
#line 230 "MB03PY.f"
	*info = -6;
#line 231 "MB03PY.f"
    }

#line 233 "MB03PY.f"
    if (*info != 0) {
#line 234 "MB03PY.f"
	i__1 = -(*info);
#line 234 "MB03PY.f"
	xerbla_("MB03PY", &i__1, (ftnlen)6);
#line 235 "MB03PY.f"
	return 0;
#line 236 "MB03PY.f"
    }

/*     Quick return if possible. */

#line 240 "MB03PY.f"
    k = min(*m,*n);
#line 241 "MB03PY.f"
    if (k == 0) {
#line 242 "MB03PY.f"
	*rank = 0;
#line 243 "MB03PY.f"
	sval[1] = 0.;
#line 244 "MB03PY.f"
	sval[2] = 0.;
#line 245 "MB03PY.f"
	sval[3] = 0.;
#line 246 "MB03PY.f"
	return 0;
#line 247 "MB03PY.f"
    }

#line 249 "MB03PY.f"
    ismin = *m;
#line 250 "MB03PY.f"
    ismax = ismin + *m;
#line 251 "MB03PY.f"
    jwork = ismax + 1;

/*     Initialize partial row norms and pivoting vector. The first m */
/*     elements of DWORK store the exact row norms. The already used */
/*     trailing part is then overwritten by the condition estimator. */

#line 257 "MB03PY.f"
    i__1 = *m;
#line 257 "MB03PY.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 258 "MB03PY.f"
	dwork[i__] = dnrm2_(n, &a[i__ + a_dim1], lda);
#line 259 "MB03PY.f"
	dwork[*m + i__] = dwork[i__];
#line 260 "MB03PY.f"
	jpvt[i__] = i__;
#line 261 "MB03PY.f"
/* L10: */
#line 261 "MB03PY.f"
    }

/*     Compute factorization and determine RANK using incremental */
/*     condition estimation. */

#line 266 "MB03PY.f"
    *rank = 0;

#line 268 "MB03PY.f"
L20:
#line 269 "MB03PY.f"
    if (*rank < k) {
#line 270 "MB03PY.f"
	i__ = k - *rank;

/*        Determine ith pivot row and swap if necessary. */

#line 274 "MB03PY.f"
	mki = *m - *rank;
#line 275 "MB03PY.f"
	nki = *n - *rank;
#line 276 "MB03PY.f"
	pvt = idamax_(&mki, &dwork[1], &c__1);

#line 278 "MB03PY.f"
	if (pvt != mki) {
#line 279 "MB03PY.f"
	    dswap_(n, &a[pvt + a_dim1], lda, &a[mki + a_dim1], lda);
#line 280 "MB03PY.f"
	    itemp = jpvt[pvt];
#line 281 "MB03PY.f"
	    jpvt[pvt] = jpvt[mki];
#line 282 "MB03PY.f"
	    jpvt[mki] = itemp;
#line 283 "MB03PY.f"
	    dwork[pvt] = dwork[mki];
#line 284 "MB03PY.f"
	    dwork[*m + pvt] = dwork[*m + mki];
#line 285 "MB03PY.f"
	}

#line 287 "MB03PY.f"
	if (nki > 1) {

/*           Save A(m-k+i,n-k+i) and generate elementary reflector H(i) */
/*           to annihilate A(m-k+i,1:n-k+i-1), k = min(m,n). */

#line 292 "MB03PY.f"
	    aii = a[mki + nki * a_dim1];
#line 293 "MB03PY.f"
	    dlarfg_(&nki, &a[mki + nki * a_dim1], &a[mki + a_dim1], lda, &tau[
		    i__]);
#line 295 "MB03PY.f"
	}

#line 297 "MB03PY.f"
	if (*rank == 0) {

/*           Initialize; exit if matrix is zero (RANK = 0). */

#line 301 "MB03PY.f"
	    smax = (d__1 = a[*m + *n * a_dim1], abs(d__1));
#line 302 "MB03PY.f"
	    if (smax == 0.) {
#line 303 "MB03PY.f"
		sval[1] = 0.;
#line 304 "MB03PY.f"
		sval[2] = 0.;
#line 305 "MB03PY.f"
		sval[3] = 0.;
#line 306 "MB03PY.f"
		return 0;
#line 307 "MB03PY.f"
	    }
#line 308 "MB03PY.f"
	    smin = smax;
#line 309 "MB03PY.f"
	    smaxpr = smax;
#line 310 "MB03PY.f"
	    sminpr = smin;
#line 311 "MB03PY.f"
	    c1 = 1.;
#line 312 "MB03PY.f"
	    c2 = 1.;
#line 313 "MB03PY.f"
	} else {

/*           One step of incremental condition estimation. */

#line 317 "MB03PY.f"
	    dcopy_(rank, &a[mki + (nki + 1) * a_dim1], lda, &dwork[jwork], &
		    c__1);
#line 318 "MB03PY.f"
	    dlaic1_(&c__2, rank, &dwork[ismin], &smin, &dwork[jwork], &a[mki 
		    + nki * a_dim1], &sminpr, &s1, &c1);
#line 320 "MB03PY.f"
	    dlaic1_(&c__1, rank, &dwork[ismax], &smax, &dwork[jwork], &a[mki 
		    + nki * a_dim1], &smaxpr, &s2, &c2);
#line 322 "MB03PY.f"
	}

#line 324 "MB03PY.f"
	if (*svlmax * *rcond <= smaxpr) {
#line 325 "MB03PY.f"
	    if (*svlmax * *rcond <= sminpr) {
#line 326 "MB03PY.f"
		if (smaxpr * *rcond <= sminpr) {

#line 328 "MB03PY.f"
		    if (mki > 1) {

/*                    Continue factorization, as rank is at least RANK. */
/*                    Apply H(i) to A(1:m-k+i-1,1:n-k+i) from the right. */

#line 333 "MB03PY.f"
			aii = a[mki + nki * a_dim1];
#line 334 "MB03PY.f"
			a[mki + nki * a_dim1] = 1.;
#line 335 "MB03PY.f"
			i__1 = mki - 1;
#line 335 "MB03PY.f"
			dlarf_("Right", &i__1, &nki, &a[mki + a_dim1], lda, &
				tau[i__], &a[a_offset], lda, &dwork[jwork], (
				ftnlen)5);
#line 337 "MB03PY.f"
			a[mki + nki * a_dim1] = aii;

/*                    Update partial row norms. */

#line 341 "MB03PY.f"
			i__1 = mki - 1;
#line 341 "MB03PY.f"
			for (j = 1; j <= i__1; ++j) {
#line 342 "MB03PY.f"
			    if (dwork[j] != 0.) {
/* Computing 2nd power */
#line 343 "MB03PY.f"
				d__2 = (d__1 = a[j + nki * a_dim1], abs(d__1))
					 / dwork[j];
#line 343 "MB03PY.f"
				temp = 1. - d__2 * d__2;
#line 345 "MB03PY.f"
				temp = max(temp,0.);
/* Computing 2nd power */
#line 346 "MB03PY.f"
				d__1 = dwork[j] / dwork[*m + j];
#line 346 "MB03PY.f"
				temp2 = temp * .05 * (d__1 * d__1) + 1.;
#line 348 "MB03PY.f"
				if (temp2 == 1.) {
#line 349 "MB03PY.f"
				    i__2 = nki - 1;
#line 349 "MB03PY.f"
				    dwork[j] = dnrm2_(&i__2, &a[j + a_dim1], 
					    lda);
#line 351 "MB03PY.f"
				    dwork[*m + j] = dwork[j];
#line 352 "MB03PY.f"
				} else {
#line 353 "MB03PY.f"
				    dwork[j] *= sqrt(temp);
#line 354 "MB03PY.f"
				}
#line 355 "MB03PY.f"
			    }
#line 356 "MB03PY.f"
/* L30: */
#line 356 "MB03PY.f"
			}

#line 358 "MB03PY.f"
		    }

#line 360 "MB03PY.f"
		    i__1 = *rank;
#line 360 "MB03PY.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 361 "MB03PY.f"
			dwork[ismin + i__ - 1] = s1 * dwork[ismin + i__ - 1];
#line 362 "MB03PY.f"
			dwork[ismax + i__ - 1] = s2 * dwork[ismax + i__ - 1];
#line 363 "MB03PY.f"
/* L40: */
#line 363 "MB03PY.f"
		    }

#line 365 "MB03PY.f"
		    if (*rank > 0) {
#line 366 "MB03PY.f"
			--ismin;
#line 367 "MB03PY.f"
			--ismax;
#line 368 "MB03PY.f"
		    }
#line 369 "MB03PY.f"
		    dwork[ismin] = c1;
#line 370 "MB03PY.f"
		    dwork[ismax] = c2;
#line 371 "MB03PY.f"
		    smin = sminpr;
#line 372 "MB03PY.f"
		    smax = smaxpr;
#line 373 "MB03PY.f"
		    ++(*rank);
#line 374 "MB03PY.f"
		    goto L20;
#line 375 "MB03PY.f"
		}
#line 376 "MB03PY.f"
	    }
#line 377 "MB03PY.f"
	}
#line 378 "MB03PY.f"
    }

/*     Restore the changed part of the (M-RANK)-th row and set SVAL. */

#line 382 "MB03PY.f"
    if (*rank < k && nki > 1) {
#line 383 "MB03PY.f"
	i__1 = nki - 1;
#line 383 "MB03PY.f"
	d__1 = -a[mki + nki * a_dim1] * tau[i__];
#line 383 "MB03PY.f"
	dscal_(&i__1, &d__1, &a[mki + a_dim1], lda);
#line 384 "MB03PY.f"
	a[mki + nki * a_dim1] = aii;
#line 385 "MB03PY.f"
    }
#line 386 "MB03PY.f"
    sval[1] = smax;
#line 387 "MB03PY.f"
    sval[2] = smin;
#line 388 "MB03PY.f"
    sval[3] = sminpr;

#line 390 "MB03PY.f"
    return 0;
/* *** Last line of MB03PY *** */
} /* mb03py_ */

