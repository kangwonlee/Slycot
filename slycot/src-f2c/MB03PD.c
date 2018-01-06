#line 1 "MB03PD.f"
/* MB03PD.f -- translated by f2c (version 20100827).
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

#line 1 "MB03PD.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int mb03pd_(char *jobrq, integer *m, integer *n, doublereal *
	a, integer *lda, integer *jpvt, doublereal *rcond, doublereal *svlmax,
	 doublereal *tau, integer *rank, doublereal *sval, doublereal *dwork, 
	integer *info, ftnlen jobrq_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal c1, c2, s1, s2;
    static integer mn;
    static doublereal smin, smax;
    extern /* Subroutine */ int mb04gd_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer ismin, ismax, jwork;
    extern /* Subroutine */ int dlaic1_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), xerbla_(char *, integer *, ftnlen);
    static logical ljobrq;
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

/*     To compute (optionally) a rank-revealing RQ factorization of a */
/*     real general M-by-N matrix  A,  which may be rank-deficient, */
/*     and estimate its effective rank using incremental condition */
/*     estimation. */

/*     The routine uses an RQ factorization with row pivoting: */
/*        P * A = R * Q,  where  R = [ R11 R12 ], */
/*                                   [  0  R22 ] */
/*     with R22 defined as the largest trailing submatrix whose estimated */
/*     condition number is less than 1/RCOND.  The order of R22, RANK, */
/*     is the effective rank of A. */

/*     MB03PD  does not perform any scaling of the matrix A. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBRQ   CHARACTER*1 */
/*             = 'R':  Perform an RQ factorization with row pivoting; */
/*             = 'N':  Do not perform the RQ factorization (but assume */
/*                     that it has been done outside). */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension */
/*             ( LDA, N ) */
/*             On entry with JOBRQ = 'R', the leading M-by-N part of this */
/*             array must contain the given matrix A. */
/*             On exit with JOBRQ = 'R', */
/*             if M <= N, the upper triangle of the subarray */
/*             A(1:M,N-M+1:N) contains the M-by-M upper triangular */
/*             matrix R; */
/*             if M >= N, the elements on and above the (M-N)-th */
/*             subdiagonal contain the M-by-N upper trapezoidal matrix R; */
/*             the remaining elements, with the array TAU, represent the */
/*             orthogonal matrix Q as a product of min(M,N) elementary */
/*             reflectors (see METHOD). */
/*             On entry and on exit with JOBRQ = 'N', */
/*             if M <= N, the upper triangle of the subarray */
/*             A(1:M,N-M+1:N) must contain the M-by-M upper triangular */
/*             matrix R; */
/*             if M >= N, the elements on and above the (M-N)-th */
/*             subdiagonal must contain the M-by-N upper trapezoidal */
/*             matrix R; */
/*             the remaining elements are not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     JPVT    (input/output) INTEGER array, dimension ( M ) */
/*             On entry with JOBRQ = 'R', if JPVT(i) <> 0, the i-th row */
/*             of A is a final row, otherwise it is a free row. Before */
/*             the RQ factorization of A, all final rows are permuted */
/*             to the trailing positions; only the remaining free rows */
/*             are moved as a result of row pivoting during the */
/*             factorization.  For rank determination it is preferable */
/*             that all rows be free. */
/*             On exit with JOBRQ = 'R', if JPVT(i) = k, then the i-th */
/*             row of P*A was the k-th row of A. */
/*             Array JPVT is not referenced when JOBRQ = 'N'. */

/*     RCOND   (input) DOUBLE PRECISION */
/*             RCOND is used to determine the effective rank of A, which */
/*             is defined as the order of the largest trailing triangular */
/*             submatrix R22 in the RQ factorization with pivoting of A, */
/*             whose estimated condition number is less than 1/RCOND. */
/*             RCOND >= 0. */
/*             NOTE that when SVLMAX > 0, the estimated rank could be */
/*             less than that defined above (see SVLMAX). */

/*     SVLMAX  (input) DOUBLE PRECISION */
/*             If A is a submatrix of another matrix B, and the rank */
/*             decision should be related to that matrix, then SVLMAX */
/*             should be an estimate of the largest singular value of B */
/*             (for instance, the Frobenius norm of B).  If this is not */
/*             the case, the input value SVLMAX = 0 should work. */
/*             SVLMAX >= 0. */

/*     TAU     (output) DOUBLE PRECISION array, dimension ( MIN( M, N ) ) */
/*             On exit with JOBRQ = 'R', the leading min(M,N) elements of */
/*             TAU contain the scalar factors of the elementary */
/*             reflectors. */
/*             Array TAU is not referenced when JOBRQ = 'N'. */

/*     RANK    (output) INTEGER */
/*             The effective (estimated) rank of A, i.e. the order of */
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

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension ( LDWORK ) */
/*             where LDWORK = max( 1, 3*M ),           if JOBRQ = 'R'; */
/*                   LDWORK = max( 1, 3*min( M, N ) ), if JOBRQ = 'N'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine computes or uses an RQ factorization with row */
/*     pivoting of A,  P * A = R * Q,  with  R  defined above, and then */
/*     finds the largest trailing submatrix whose estimated condition */
/*     number is less than 1/RCOND, taking the possible positive value of */
/*     SVLMAX into account.  This is performed using an adaptation of the */
/*     LAPACK incremental condition estimation scheme and a slightly */
/*     modified rank decision test. */

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

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996. */

/*     REVISIONS */

/*     Nov. 1997 */

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
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 239 "MB03PD.f"
    /* Parameter adjustments */
#line 239 "MB03PD.f"
    a_dim1 = *lda;
#line 239 "MB03PD.f"
    a_offset = 1 + a_dim1;
#line 239 "MB03PD.f"
    a -= a_offset;
#line 239 "MB03PD.f"
    --jpvt;
#line 239 "MB03PD.f"
    --tau;
#line 239 "MB03PD.f"
    --sval;
#line 239 "MB03PD.f"
    --dwork;
#line 239 "MB03PD.f"

#line 239 "MB03PD.f"
    /* Function Body */
#line 239 "MB03PD.f"
    ljobrq = lsame_(jobrq, "R", (ftnlen)1, (ftnlen)1);
#line 240 "MB03PD.f"
    mn = min(*m,*n);

/*     Test the input scalar arguments. */

#line 244 "MB03PD.f"
    *info = 0;
#line 245 "MB03PD.f"
    if (! ljobrq && ! lsame_(jobrq, "N", (ftnlen)1, (ftnlen)1)) {
#line 246 "MB03PD.f"
	*info = -1;
#line 247 "MB03PD.f"
    } else if (*m < 0) {
#line 248 "MB03PD.f"
	*info = -2;
#line 249 "MB03PD.f"
    } else if (*n < 0) {
#line 250 "MB03PD.f"
	*info = -3;
#line 251 "MB03PD.f"
    } else if (*lda < max(1,*m)) {
#line 252 "MB03PD.f"
	*info = -5;
#line 253 "MB03PD.f"
    } else if (*rcond < 0.) {
#line 254 "MB03PD.f"
	*info = -7;
#line 255 "MB03PD.f"
    } else if (*svlmax < 0.) {
#line 256 "MB03PD.f"
	*info = -8;
#line 257 "MB03PD.f"
    }

#line 259 "MB03PD.f"
    if (*info != 0) {
#line 260 "MB03PD.f"
	i__1 = -(*info);
#line 260 "MB03PD.f"
	xerbla_("MB03PD", &i__1, (ftnlen)6);
#line 261 "MB03PD.f"
	return 0;
#line 262 "MB03PD.f"
    }

/*     Quick return if possible. */

#line 266 "MB03PD.f"
    if (mn == 0) {
#line 267 "MB03PD.f"
	*rank = 0;
#line 268 "MB03PD.f"
	sval[1] = 0.;
#line 269 "MB03PD.f"
	sval[2] = 0.;
#line 270 "MB03PD.f"
	sval[3] = 0.;
#line 271 "MB03PD.f"
	return 0;
#line 272 "MB03PD.f"
    }

#line 274 "MB03PD.f"
    if (ljobrq) {

/*        Compute RQ factorization with row pivoting of A: */
/*           P * A = R * Q */
/*        Workspace 3*M. Details of Householder rotations stored in TAU. */

#line 280 "MB03PD.f"
	mb04gd_(m, n, &a[a_offset], lda, &jpvt[1], &tau[1], &dwork[1], info);
#line 281 "MB03PD.f"
    }

/*     Determine RANK using incremental condition estimation. */
/*        Workspace 3*min(M,N). */

#line 286 "MB03PD.f"
    smax = (d__1 = a[*m + *n * a_dim1], abs(d__1));
#line 287 "MB03PD.f"
    if (smax == 0. || *svlmax * *rcond > smax) {
#line 288 "MB03PD.f"
	*rank = 0;
#line 289 "MB03PD.f"
	sval[1] = smax;
#line 290 "MB03PD.f"
	sval[2] = 0.;
#line 291 "MB03PD.f"
	sval[3] = 0.;
#line 292 "MB03PD.f"
    } else {
#line 293 "MB03PD.f"
	ismin = mn;
#line 294 "MB03PD.f"
	ismax = mn << 1;
#line 295 "MB03PD.f"
	jwork = ismax + 1;
#line 296 "MB03PD.f"
	dwork[ismin] = 1.;
#line 297 "MB03PD.f"
	dwork[ismax] = 1.;
#line 298 "MB03PD.f"
	*rank = 1;
#line 299 "MB03PD.f"
	smin = smax;
#line 300 "MB03PD.f"
	sminpr = smin;

#line 302 "MB03PD.f"
L10:
#line 303 "MB03PD.f"
	if (*rank < mn) {
#line 304 "MB03PD.f"
	    dcopy_(rank, &a[*m - *rank + (*n - *rank + 1) * a_dim1], lda, &
		    dwork[jwork], &c__1);
#line 306 "MB03PD.f"
	    dlaic1_(&c__2, rank, &dwork[ismin], &smin, &dwork[jwork], &a[*m - 
		    *rank + (*n - *rank) * a_dim1], &sminpr, &s1, &c1);
#line 309 "MB03PD.f"
	    dlaic1_(&c__1, rank, &dwork[ismax], &smax, &dwork[jwork], &a[*m - 
		    *rank + (*n - *rank) * a_dim1], &smaxpr, &s2, &c2);

#line 313 "MB03PD.f"
	    if (*svlmax * *rcond <= smaxpr) {
#line 314 "MB03PD.f"
		if (*svlmax * *rcond <= sminpr) {
#line 315 "MB03PD.f"
		    if (smaxpr * *rcond <= sminpr) {
#line 316 "MB03PD.f"
			i__1 = *rank;
#line 316 "MB03PD.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 317 "MB03PD.f"
			    dwork[ismin + i__ - 1] = s1 * dwork[ismin + i__ - 
				    1];
#line 318 "MB03PD.f"
			    dwork[ismax + i__ - 1] = s2 * dwork[ismax + i__ - 
				    1];
#line 319 "MB03PD.f"
/* L20: */
#line 319 "MB03PD.f"
			}
#line 320 "MB03PD.f"
			--ismin;
#line 321 "MB03PD.f"
			--ismax;
#line 322 "MB03PD.f"
			dwork[ismin] = c1;
#line 323 "MB03PD.f"
			dwork[ismax] = c2;
#line 324 "MB03PD.f"
			smin = sminpr;
#line 325 "MB03PD.f"
			smax = smaxpr;
#line 326 "MB03PD.f"
			++(*rank);
#line 327 "MB03PD.f"
			goto L10;
#line 328 "MB03PD.f"
		    }
#line 329 "MB03PD.f"
		}
#line 330 "MB03PD.f"
	    }
#line 331 "MB03PD.f"
	}
#line 332 "MB03PD.f"
	sval[1] = smax;
#line 333 "MB03PD.f"
	sval[2] = smin;
#line 334 "MB03PD.f"
	sval[3] = sminpr;
#line 335 "MB03PD.f"
    }

#line 337 "MB03PD.f"
    return 0;
/* *** Last line of MB03PD *** */
} /* mb03pd_ */

