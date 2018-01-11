#line 1 "MB03OY.f"
/* MB03OY.f -- translated by f2c (version 20100827).
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

#line 1 "MB03OY.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int mb03oy_(integer *m, integer *n, doublereal *a, integer *
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
    static integer i__, j;
    static doublereal c1, c2, s1, s2;
    static integer mn;
    static doublereal aii;
    static integer pvt;
    static doublereal smin, temp, smax;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal temp2;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dlarf_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    ftnlen);
    static integer itemp, ismin;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer ismax;
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

/*     To compute a rank-revealing QR factorization of a real general */
/*     M-by-N matrix  A,  which may be rank-deficient, and estimate its */
/*     effective rank using incremental condition estimation. */

/*     The routine uses a truncated QR factorization with column pivoting */
/*                                   [ R11 R12 ] */
/*        A * P = Q * R,  where  R = [         ], */
/*                                   [  0  R22 ] */
/*     with R11 defined as the largest leading upper triangular submatrix */
/*     whose estimated condition number is less than 1/RCOND.  The order */
/*     of R11, RANK, is the effective rank of A.  Condition estimation is */
/*     performed during the QR factorization process.  Matrix R22 is full */
/*     (but of small norm), or empty. */

/*     MB03OY  does not perform any scaling of the matrix A. */

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
/*             On exit, the leading RANK-by-RANK upper triangular part */
/*             of A contains the triangular factor R11, and the elements */
/*             below the diagonal in the first  RANK  columns, with the */
/*             array TAU, represent the orthogonal matrix Q as a product */
/*             of  RANK  elementary reflectors. */
/*             The remaining  N-RANK  columns contain the result of the */
/*             QR factorization process used. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     RCOND   (input) DOUBLE PRECISION */
/*             RCOND is used to determine the effective rank of A, which */
/*             is defined as the order of the largest leading triangular */
/*             submatrix R11 in the QR factorization with pivoting of A, */
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
/*             the submatrix R11. */

/*     SVAL    (output) DOUBLE PRECISION array, dimension ( 3 ) */
/*             The estimates of some of the singular values of the */
/*             triangular factor R: */
/*             SVAL(1): largest singular value of R(1:RANK,1:RANK); */
/*             SVAL(2): smallest singular value of R(1:RANK,1:RANK); */
/*             SVAL(3): smallest singular value of R(1:RANK+1,1:RANK+1), */
/*                      if RANK < MIN( M, N ), or of R(1:RANK,1:RANK), */
/*                      otherwise. */
/*             If the triangular factorization is a rank-revealing one */
/*             (which will be the case if the leading columns were well- */
/*             conditioned), then SVAL(1) will also be an estimate for */
/*             the largest singular value of A, and SVAL(2) and SVAL(3) */
/*             will be estimates for the RANK-th and (RANK+1)-st singular */
/*             values of A, respectively. */
/*             By examining these values, one can confirm that the rank */
/*             is well defined with respect to the chosen value of RCOND. */
/*             The ratio SVAL(1)/SVAL(2) is an estimate of the condition */
/*             number of R(1:RANK,1:RANK). */

/*     JPVT    (output) INTEGER array, dimension ( N ) */
/*             If JPVT(i) = k, then the i-th column of A*P was the k-th */
/*             column of A. */

/*     TAU     (output) DOUBLE PRECISION array, dimension ( MIN( M, N ) ) */
/*             The leading  RANK  elements of TAU contain the scalar */
/*             factors of the elementary reflectors. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension ( 3*N-1 ) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine computes a truncated QR factorization with column */
/*     pivoting of A,  A * P = Q * R,  with  R  defined above, and, */
/*     during this process, finds the largest leading submatrix whose */
/*     estimated condition number is less than 1/RCOND, taking the */
/*     possible positive value of SVLMAX into account.  This is performed */
/*     using the LAPACK incremental condition estimation scheme and a */
/*     slightly modified rank decision test.  The factorization process */
/*     stops when  RANK  has been determined. */

/*     The matrix Q is represented as a product of elementary reflectors */

/*        Q = H(1) H(2) . . . H(k), where k = rank <= min(m,n). */

/*     Each H(i) has the form */

/*        H = I - tau * v * v' */

/*     where tau is a real scalar, and v is a real vector with */
/*     v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in */
/*     A(i+1:m,i), and tau in TAU(i). */

/*     The matrix P is represented in jpvt as follows: If */
/*        jpvt(j) = i */
/*     then the jth column of P is the ith canonical unit vector. */

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

/*     V. Sima, Research Institute for Informatics, Bucharest, Jan. 2009. */

/*     KEYWORDS */

/*     Eigenvalue problem, matrix operations, orthogonal transformation, */
/*     singular values. */

/*    ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

#line 214 "MB03OY.f"
    /* Parameter adjustments */
#line 214 "MB03OY.f"
    a_dim1 = *lda;
#line 214 "MB03OY.f"
    a_offset = 1 + a_dim1;
#line 214 "MB03OY.f"
    a -= a_offset;
#line 214 "MB03OY.f"
    --sval;
#line 214 "MB03OY.f"
    --jpvt;
#line 214 "MB03OY.f"
    --tau;
#line 214 "MB03OY.f"
    --dwork;
#line 214 "MB03OY.f"

#line 214 "MB03OY.f"
    /* Function Body */
#line 214 "MB03OY.f"
    *info = 0;
#line 215 "MB03OY.f"
    if (*m < 0) {
#line 216 "MB03OY.f"
	*info = -1;
#line 217 "MB03OY.f"
    } else if (*n < 0) {
#line 218 "MB03OY.f"
	*info = -2;
#line 219 "MB03OY.f"
    } else if (*lda < max(1,*m)) {
#line 220 "MB03OY.f"
	*info = -4;
#line 221 "MB03OY.f"
    } else if (*rcond < 0. || *rcond > 1.) {
#line 222 "MB03OY.f"
	*info = -5;
#line 223 "MB03OY.f"
    } else if (*svlmax < 0.) {
#line 224 "MB03OY.f"
	*info = -6;
#line 225 "MB03OY.f"
    }

#line 227 "MB03OY.f"
    if (*info != 0) {
#line 228 "MB03OY.f"
	i__1 = -(*info);
#line 228 "MB03OY.f"
	xerbla_("MB03OY", &i__1, (ftnlen)6);
#line 229 "MB03OY.f"
	return 0;
#line 230 "MB03OY.f"
    }

/*     Quick return if possible. */

#line 234 "MB03OY.f"
    mn = min(*m,*n);
#line 235 "MB03OY.f"
    if (mn == 0) {
#line 236 "MB03OY.f"
	*rank = 0;
#line 237 "MB03OY.f"
	sval[1] = 0.;
#line 238 "MB03OY.f"
	sval[2] = 0.;
#line 239 "MB03OY.f"
	sval[3] = 0.;
#line 240 "MB03OY.f"
	return 0;
#line 241 "MB03OY.f"
    }

#line 243 "MB03OY.f"
    ismin = 1;
#line 244 "MB03OY.f"
    ismax = ismin + *n;

/*     Initialize partial column norms and pivoting vector. The first n */
/*     elements of DWORK store the exact column norms. The already used */
/*     leading part is then overwritten by the condition estimator. */

#line 250 "MB03OY.f"
    i__1 = *n;
#line 250 "MB03OY.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 251 "MB03OY.f"
	dwork[i__] = dnrm2_(m, &a[i__ * a_dim1 + 1], &c__1);
#line 252 "MB03OY.f"
	dwork[*n + i__] = dwork[i__];
#line 253 "MB03OY.f"
	jpvt[i__] = i__;
#line 254 "MB03OY.f"
/* L10: */
#line 254 "MB03OY.f"
    }

/*     Compute factorization and determine RANK using incremental */
/*     condition estimation. */

#line 259 "MB03OY.f"
    *rank = 0;

#line 261 "MB03OY.f"
L20:
#line 262 "MB03OY.f"
    if (*rank < mn) {
#line 263 "MB03OY.f"
	i__ = *rank + 1;

/*        Determine ith pivot column and swap if necessary. */

#line 267 "MB03OY.f"
	i__1 = *n - i__ + 1;
#line 267 "MB03OY.f"
	pvt = i__ - 1 + idamax_(&i__1, &dwork[i__], &c__1);

#line 269 "MB03OY.f"
	if (pvt != i__) {
#line 270 "MB03OY.f"
	    dswap_(m, &a[pvt * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], &
		    c__1);
#line 271 "MB03OY.f"
	    itemp = jpvt[pvt];
#line 272 "MB03OY.f"
	    jpvt[pvt] = jpvt[i__];
#line 273 "MB03OY.f"
	    jpvt[i__] = itemp;
#line 274 "MB03OY.f"
	    dwork[pvt] = dwork[i__];
#line 275 "MB03OY.f"
	    dwork[*n + pvt] = dwork[*n + i__];
#line 276 "MB03OY.f"
	}

/*        Save A(I,I) and generate elementary reflector H(i). */

#line 280 "MB03OY.f"
	if (i__ < *m) {
#line 281 "MB03OY.f"
	    aii = a[i__ + i__ * a_dim1];
#line 282 "MB03OY.f"
	    i__1 = *m - i__ + 1;
#line 282 "MB03OY.f"
	    dlarfg_(&i__1, &a[i__ + i__ * a_dim1], &a[i__ + 1 + i__ * a_dim1],
		     &c__1, &tau[i__]);
#line 283 "MB03OY.f"
	} else {
#line 284 "MB03OY.f"
	    tau[*m] = 0.;
#line 285 "MB03OY.f"
	}

#line 287 "MB03OY.f"
	if (*rank == 0) {

/*           Initialize; exit if matrix is zero (RANK = 0). */

#line 291 "MB03OY.f"
	    smax = (d__1 = a[a_dim1 + 1], abs(d__1));
#line 292 "MB03OY.f"
	    if (smax == 0.) {
#line 293 "MB03OY.f"
		sval[1] = 0.;
#line 294 "MB03OY.f"
		sval[2] = 0.;
#line 295 "MB03OY.f"
		sval[3] = 0.;
#line 296 "MB03OY.f"
		return 0;
#line 297 "MB03OY.f"
	    }
#line 298 "MB03OY.f"
	    smin = smax;
#line 299 "MB03OY.f"
	    smaxpr = smax;
#line 300 "MB03OY.f"
	    sminpr = smin;
#line 301 "MB03OY.f"
	    c1 = 1.;
#line 302 "MB03OY.f"
	    c2 = 1.;
#line 303 "MB03OY.f"
	} else {

/*           One step of incremental condition estimation. */

#line 307 "MB03OY.f"
	    dlaic1_(&c__2, rank, &dwork[ismin], &smin, &a[i__ * a_dim1 + 1], &
		    a[i__ + i__ * a_dim1], &sminpr, &s1, &c1);
#line 309 "MB03OY.f"
	    dlaic1_(&c__1, rank, &dwork[ismax], &smax, &a[i__ * a_dim1 + 1], &
		    a[i__ + i__ * a_dim1], &smaxpr, &s2, &c2);
#line 311 "MB03OY.f"
	}

#line 313 "MB03OY.f"
	if (*svlmax * *rcond <= smaxpr) {
#line 314 "MB03OY.f"
	    if (*svlmax * *rcond <= sminpr) {
#line 315 "MB03OY.f"
		if (smaxpr * *rcond <= sminpr) {

/*                 Continue factorization, as rank is at least RANK. */

#line 319 "MB03OY.f"
		    if (i__ < *n) {

/*                    Apply H(i) to A(i:m,i+1:n) from the left. */

#line 323 "MB03OY.f"
			aii = a[i__ + i__ * a_dim1];
#line 324 "MB03OY.f"
			a[i__ + i__ * a_dim1] = 1.;
#line 325 "MB03OY.f"
			i__1 = *m - i__ + 1;
#line 325 "MB03OY.f"
			i__2 = *n - i__;
#line 325 "MB03OY.f"
			dlarf_("Left", &i__1, &i__2, &a[i__ + i__ * a_dim1], &
				c__1, &tau[i__], &a[i__ + (i__ + 1) * a_dim1],
				 lda, &dwork[(*n << 1) + 1], (ftnlen)4);
#line 328 "MB03OY.f"
			a[i__ + i__ * a_dim1] = aii;
#line 329 "MB03OY.f"
		    }

/*                 Update partial column norms. */

#line 333 "MB03OY.f"
		    i__1 = *n;
#line 333 "MB03OY.f"
		    for (j = i__ + 1; j <= i__1; ++j) {
#line 334 "MB03OY.f"
			if (dwork[j] != 0.) {
/* Computing 2nd power */
#line 335 "MB03OY.f"
			    d__2 = (d__1 = a[i__ + j * a_dim1], abs(d__1)) / 
				    dwork[j];
#line 335 "MB03OY.f"
			    temp = 1. - d__2 * d__2;
#line 337 "MB03OY.f"
			    temp = max(temp,0.);
/* Computing 2nd power */
#line 338 "MB03OY.f"
			    d__1 = dwork[j] / dwork[*n + j];
#line 338 "MB03OY.f"
			    temp2 = temp * .05 * (d__1 * d__1) + 1.;
#line 340 "MB03OY.f"
			    if (temp2 == 1.) {
#line 341 "MB03OY.f"
				if (*m - i__ > 0) {
#line 342 "MB03OY.f"
				    i__2 = *m - i__;
#line 342 "MB03OY.f"
				    dwork[j] = dnrm2_(&i__2, &a[i__ + 1 + j * 
					    a_dim1], &c__1);
#line 343 "MB03OY.f"
				    dwork[*n + j] = dwork[j];
#line 344 "MB03OY.f"
				} else {
#line 345 "MB03OY.f"
				    dwork[j] = 0.;
#line 346 "MB03OY.f"
				    dwork[*n + j] = 0.;
#line 347 "MB03OY.f"
				}
#line 348 "MB03OY.f"
			    } else {
#line 349 "MB03OY.f"
				dwork[j] *= sqrt(temp);
#line 350 "MB03OY.f"
			    }
#line 351 "MB03OY.f"
			}
#line 352 "MB03OY.f"
/* L30: */
#line 352 "MB03OY.f"
		    }

#line 354 "MB03OY.f"
		    i__1 = *rank;
#line 354 "MB03OY.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 355 "MB03OY.f"
			dwork[ismin + i__ - 1] = s1 * dwork[ismin + i__ - 1];
#line 356 "MB03OY.f"
			dwork[ismax + i__ - 1] = s2 * dwork[ismax + i__ - 1];
#line 357 "MB03OY.f"
/* L40: */
#line 357 "MB03OY.f"
		    }

#line 359 "MB03OY.f"
		    dwork[ismin + *rank] = c1;
#line 360 "MB03OY.f"
		    dwork[ismax + *rank] = c2;
#line 361 "MB03OY.f"
		    smin = sminpr;
#line 362 "MB03OY.f"
		    smax = smaxpr;
#line 363 "MB03OY.f"
		    ++(*rank);
#line 364 "MB03OY.f"
		    goto L20;
#line 365 "MB03OY.f"
		}
#line 366 "MB03OY.f"
	    }
#line 367 "MB03OY.f"
	}
#line 368 "MB03OY.f"
    }

/*     Restore the changed part of the (RANK+1)-th column and set SVAL. */

#line 372 "MB03OY.f"
    if (*rank < *n) {
#line 373 "MB03OY.f"
	if (i__ < *m) {
#line 374 "MB03OY.f"
	    i__1 = *m - i__;
#line 374 "MB03OY.f"
	    d__1 = -a[i__ + i__ * a_dim1] * tau[i__];
#line 374 "MB03OY.f"
	    dscal_(&i__1, &d__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
#line 375 "MB03OY.f"
	    a[i__ + i__ * a_dim1] = aii;
#line 376 "MB03OY.f"
	}
#line 377 "MB03OY.f"
    }
#line 378 "MB03OY.f"
    if (*rank == 0) {
#line 379 "MB03OY.f"
	smin = 0.;
#line 380 "MB03OY.f"
	sminpr = 0.;
#line 381 "MB03OY.f"
    }
#line 382 "MB03OY.f"
    sval[1] = smax;
#line 383 "MB03OY.f"
    sval[2] = smin;
#line 384 "MB03OY.f"
    sval[3] = sminpr;

#line 386 "MB03OY.f"
    return 0;
/* *** Last line of MB03OY *** */
} /* mb03oy_ */

