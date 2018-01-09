#line 1 "MB02QD.f"
/* MB02QD.f -- translated by f2c (version 20100827).
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

#line 1 "MB02QD.f"
/* Table of constant values */

static integer c__0 = 0;
static doublereal c_b17 = 0.;
static doublereal c_b33 = 1.;

/* Subroutine */ int mb02qd_(char *job, char *iniper, integer *m, integer *n, 
	integer *nrhs, doublereal *rcond, doublereal *svlmax, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, doublereal *y, integer *
	jpvt, integer *rank, doublereal *sval, doublereal *dwork, integer *
	ldwork, integer *info, ftnlen job_len, ftnlen iniper_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal t1, t2;
    static integer mn;
    static doublereal anrm, bnrm;
    extern /* Subroutine */ int mb03od_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen);
    static integer iascl, ibscl;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb03oy_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *), dtrsm_(char *, 
	    char *, char *, char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen), dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dlacpy_(char *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static doublereal bignum;
    static logical leasts;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    static doublereal smlnum;
    static logical permut;
    extern /* Subroutine */ int dormrz_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dtzrzf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);


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

/*     To compute a solution, optionally corresponding to specified free */
/*     elements, to a real linear least squares problem: */

/*         minimize || A * X - B || */

/*     using a complete orthogonal factorization of the M-by-N matrix A, */
/*     which may be rank-deficient. */

/*     Several right hand side vectors b and solution vectors x can be */
/*     handled in a single call; they are stored as the columns of the */
/*     M-by-NRHS right hand side matrix B and the N-by-NRHS solution */
/*     matrix X. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies whether or not a standard least squares solution */
/*             must be computed, as follows: */
/*             = 'L':  Compute a standard least squares solution (Y = 0); */
/*             = 'F':  Compute a solution with specified free elements */
/*                     (given in Y). */

/*     INIPER  CHARACTER*1 */
/*             Specifies whether an initial column permutation, defined */
/*             by JPVT, must be performed, as follows: */
/*             = 'P':  Perform an initial column permutation; */
/*             = 'N':  Do not perform an initial column permutation. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     NRHS    (input) INTEGER */
/*             The number of right hand sides, i.e., the number of */
/*             columns of the matrices B and X.  NRHS >= 0. */

/*     RCOND   (input) DOUBLE PRECISION */
/*             RCOND is used to determine the effective rank of A, which */
/*             is defined as the order of the largest leading triangular */
/*             submatrix R11 in the QR factorization with pivoting of A, */
/*             whose estimated condition number is less than 1/RCOND. */
/*             0 <= RCOND <= 1. */

/*     SVLMAX  (input) DOUBLE PRECISION */
/*             If A is a submatrix of another matrix C, and the rank */
/*             decision should be related to that matrix, then SVLMAX */
/*             should be an estimate of the largest singular value of C */
/*             (for instance, the Frobenius norm of C).  If this is not */
/*             the case, the input value SVLMAX = 0 should work. */
/*             SVLMAX >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the given matrix A. */
/*             On exit, the leading M-by-N part of this array contains */
/*             details of its complete orthogonal factorization: */
/*             the leading RANK-by-RANK upper triangular part contains */
/*             the upper triangular factor T11 (see METHOD); */
/*             the elements below the diagonal, with the entries 2 to */
/*             min(M,N)+1 of the array DWORK, represent the orthogonal */
/*             matrix Q as a product of min(M,N) elementary reflectors */
/*             (see METHOD); */
/*             the elements of the subarray A(1:RANK,RANK+1:N), with the */
/*             next RANK entries of the array DWORK, represent the */
/*             orthogonal matrix Z as a product of RANK elementary */
/*             reflectors (see METHOD). */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     B       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDB,NRHS) */
/*             On entry, the leading M-by-NRHS part of this array must */
/*             contain the right hand side matrix B. */
/*             On exit, the leading N-by-NRHS part of this array contains */
/*             the solution matrix X. */
/*             If M >= N and RANK = N, the residual sum-of-squares for */
/*             the solution in the i-th column is given by the sum of */
/*             squares of elements N+1:M in that column. */
/*             If NRHS = 0, this array is not referenced, and the routine */
/*             returns the effective rank of A, and its QR factorization. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,M,N). */

/*     Y       (input) DOUBLE PRECISION array, dimension ( N*NRHS ) */
/*             If JOB = 'F', the elements Y(1:(N-RANK)*NRHS) are used as */
/*             free elements in computing the solution (see METHOD). */
/*             The remaining elements are not referenced. */
/*             If JOB = 'L', or NRHS = 0, this array is not referenced. */

/*     JPVT    (input/output) INTEGER array, dimension (N) */
/*             On entry with INIPER = 'P', if JPVT(i) <> 0, the i-th */
/*             column of A is an initial column, otherwise it is a free */
/*             column.  Before the QR factorization of A, all initial */
/*             columns are permuted to the leading positions; only the */
/*             remaining free columns are moved as a result of column */
/*             pivoting during the factorization. */
/*             If INIPER = 'N', JPVT need not be set on entry. */
/*             On exit, if JPVT(i) = k, then the i-th column of A*P */
/*             was the k-th column of A. */

/*     RANK    (output) INTEGER */
/*             The effective rank of A, i.e., the order of the submatrix */
/*             R11.  This is the same as the order of the submatrix T11 */
/*             in the complete orthogonal factorization of A. */

/*     SVAL    (output) DOUBLE PRECISION array, dimension ( 3 ) */
/*             The estimates of some of the singular values of the */
/*             triangular factor R11: */
/*             SVAL(1): largest singular value of  R(1:RANK,1:RANK); */
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

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension LDWORK */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK, and the entries 2 to min(M,N) + RANK + 1 */
/*             contain the scalar factors of the elementary reflectors */
/*             used in the complete orthogonal factorization of A. */
/*             Among the entries 2 to min(M,N) + 1, only the first RANK */
/*             elements are useful, if INIPER = 'N'. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= max( min(M,N)+3*N+1, 2*min(M,N)+NRHS ) */
/*             For optimum performance LDWORK should be larger. */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     If INIPER = 'P', the routine first computes a QR factorization */
/*     with column pivoting: */
/*         A * P = Q * [ R11 R12 ] */
/*                     [  0  R22 ] */
/*     with R11 defined as the largest leading submatrix whose estimated */
/*     condition number is less than 1/RCOND.  The order of R11, RANK, */
/*     is the effective rank of A. */
/*     If INIPER = 'N', the effective rank is estimated during a */
/*     truncated QR factorization (with column pivoting) process, and */
/*     the submatrix R22 is not upper triangular, but full and of small */
/*     norm. (See SLICOT Library routines MB03OD or MB03OY, respectively, */
/*     for further details.) */

/*     Then, R22 is considered to be negligible, and R12 is annihilated */
/*     by orthogonal transformations from the right, arriving at the */
/*     complete orthogonal factorization: */
/*        A * P = Q * [ T11 0 ] * Z */
/*                    [  0  0 ] */
/*     The solution is then */
/*        X = P * Z' [ inv(T11)*Q1'*B ] */
/*                   [        Y       ] */
/*     where Q1 consists of the first RANK columns of Q, and Y contains */
/*     free elements (if JOB = 'F'), or is zero (if JOB = 'L'). */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     FURTHER COMMENTS */

/*     Significant gain in efficiency is possible for small-rank problems */
/*     using truncated QR factorization (option INIPER = 'N'). */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, Technical University of Sofia, Oct. 1998, */
/*     modification of the LAPACK routine DGELSX. */
/*     V. Sima, Katholieke Universiteit Leuven, Jan. 1999, SLICOT Library */
/*     version. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2005. */

/*     KEYWORDS */

/*     Least squares problems, QR factorization. */

/*     ****************************************************************** */

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

#line 264 "MB02QD.f"
    /* Parameter adjustments */
#line 264 "MB02QD.f"
    a_dim1 = *lda;
#line 264 "MB02QD.f"
    a_offset = 1 + a_dim1;
#line 264 "MB02QD.f"
    a -= a_offset;
#line 264 "MB02QD.f"
    b_dim1 = *ldb;
#line 264 "MB02QD.f"
    b_offset = 1 + b_dim1;
#line 264 "MB02QD.f"
    b -= b_offset;
#line 264 "MB02QD.f"
    --y;
#line 264 "MB02QD.f"
    --jpvt;
#line 264 "MB02QD.f"
    --sval;
#line 264 "MB02QD.f"
    --dwork;
#line 264 "MB02QD.f"

#line 264 "MB02QD.f"
    /* Function Body */
#line 264 "MB02QD.f"
    mn = min(*m,*n);
#line 265 "MB02QD.f"
    leasts = lsame_(job, "L", (ftnlen)1, (ftnlen)1);
#line 266 "MB02QD.f"
    permut = lsame_(iniper, "P", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 270 "MB02QD.f"
    *info = 0;
/* Computing MAX */
#line 271 "MB02QD.f"
    i__1 = mn + *n * 3 + 1, i__2 = (mn << 1) + *nrhs;
#line 271 "MB02QD.f"
    minwrk = max(i__1,i__2);
#line 272 "MB02QD.f"
    if (! (leasts || lsame_(job, "F", (ftnlen)1, (ftnlen)1))) {
#line 273 "MB02QD.f"
	*info = -1;
#line 274 "MB02QD.f"
    } else if (! (permut || lsame_(iniper, "N", (ftnlen)1, (ftnlen)1))) {
#line 275 "MB02QD.f"
	*info = -2;
#line 276 "MB02QD.f"
    } else if (*m < 0) {
#line 277 "MB02QD.f"
	*info = -3;
#line 278 "MB02QD.f"
    } else if (*n < 0) {
#line 279 "MB02QD.f"
	*info = -4;
#line 280 "MB02QD.f"
    } else if (*nrhs < 0) {
#line 281 "MB02QD.f"
	*info = -5;
#line 282 "MB02QD.f"
    } else if (*rcond < 0. || *rcond > 1.) {
#line 283 "MB02QD.f"
	*info = -6;
#line 284 "MB02QD.f"
    } else if (*svlmax < 0.) {
#line 285 "MB02QD.f"
	*info = -7;
#line 286 "MB02QD.f"
    } else if (*lda < max(1,*m)) {
#line 287 "MB02QD.f"
	*info = -9;
#line 288 "MB02QD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 288 "MB02QD.f"
	i__1 = max(1,*m);
#line 288 "MB02QD.f"
	if (*ldb < max(i__1,*n)) {
#line 289 "MB02QD.f"
	    *info = -11;
#line 290 "MB02QD.f"
	} else if (*ldwork < minwrk) {
#line 291 "MB02QD.f"
	    *info = -17;
#line 292 "MB02QD.f"
	}
#line 292 "MB02QD.f"
    }

#line 294 "MB02QD.f"
    if (*info != 0) {
#line 295 "MB02QD.f"
	i__1 = -(*info);
#line 295 "MB02QD.f"
	xerbla_("MB02QD", &i__1, (ftnlen)6);
#line 296 "MB02QD.f"
	return 0;
#line 297 "MB02QD.f"
    }

/*     Quick return if possible. */

#line 301 "MB02QD.f"
    if (mn == 0) {
#line 302 "MB02QD.f"
	*rank = 0;
#line 303 "MB02QD.f"
	dwork[1] = 1.;
#line 304 "MB02QD.f"
	return 0;
#line 305 "MB02QD.f"
    }

/*     Get machine parameters. */

#line 309 "MB02QD.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12) / dlamch_("Precision", (
	    ftnlen)9);
#line 310 "MB02QD.f"
    bignum = 1. / smlnum;
#line 311 "MB02QD.f"
    dlabad_(&smlnum, &bignum);

/*     Scale A, B if max entries outside range [SMLNUM,BIGNUM]. */

#line 315 "MB02QD.f"
    anrm = dlange_("M", m, n, &a[a_offset], lda, &dwork[1], (ftnlen)1);
#line 316 "MB02QD.f"
    iascl = 0;
#line 317 "MB02QD.f"
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM. */

#line 321 "MB02QD.f"
	dlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 322 "MB02QD.f"
	iascl = 1;
#line 323 "MB02QD.f"
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM. */

#line 327 "MB02QD.f"
	dlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, 
		info, (ftnlen)1);
#line 328 "MB02QD.f"
	iascl = 2;
#line 329 "MB02QD.f"
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

#line 333 "MB02QD.f"
	if (*nrhs > 0) {
#line 333 "MB02QD.f"
	    i__1 = max(*m,*n);
#line 333 "MB02QD.f"
	    dlaset_("Full", &i__1, nrhs, &c_b17, &c_b17, &b[b_offset], ldb, (
		    ftnlen)4);
#line 333 "MB02QD.f"
	}
#line 335 "MB02QD.f"
	*rank = 0;
#line 336 "MB02QD.f"
	dwork[1] = 1.;
#line 337 "MB02QD.f"
	return 0;
#line 338 "MB02QD.f"
    }

#line 340 "MB02QD.f"
    if (*nrhs > 0) {
#line 341 "MB02QD.f"
	bnrm = dlange_("M", m, nrhs, &b[b_offset], ldb, &dwork[1], (ftnlen)1);
#line 342 "MB02QD.f"
	ibscl = 0;
#line 343 "MB02QD.f"
	if (bnrm > 0. && bnrm < smlnum) {

/*           Scale matrix norm up to SMLNUM. */

#line 347 "MB02QD.f"
	    dlascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], 
		    ldb, info, (ftnlen)1);
#line 349 "MB02QD.f"
	    ibscl = 1;
#line 350 "MB02QD.f"
	} else if (bnrm > bignum) {

/*           Scale matrix norm down to BIGNUM. */

#line 354 "MB02QD.f"
	    dlascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], 
		    ldb, info, (ftnlen)1);
#line 356 "MB02QD.f"
	    ibscl = 2;
#line 357 "MB02QD.f"
	}
#line 358 "MB02QD.f"
    }

/*     Compute a rank-revealing QR factorization of A and estimate its */
/*     effective rank using incremental condition estimation: */
/*        A * P = Q * R. */
/*     Workspace need   min(M,N)+3*N+1; */
/*               prefer min(M,N)+2*N+N*NB. */
/*     Details of Householder transformations stored in DWORK(1:MN). */
/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*      minimal amount of workspace needed at that point in the code, */
/*      as well as the preferred amount for good performance. */
/*      NB refers to the optimal block size for the immediately */
/*      following subroutine, as returned by ILAENV.) */

#line 372 "MB02QD.f"
    maxwrk = minwrk;
#line 373 "MB02QD.f"
    if (permut) {
#line 374 "MB02QD.f"
	i__1 = *ldwork - mn;
#line 374 "MB02QD.f"
	mb03od_("Q", m, n, &a[a_offset], lda, &jpvt[1], rcond, svlmax, &dwork[
		1], rank, &sval[1], &dwork[mn + 1], &i__1, info, (ftnlen)1);
/* Computing MAX */
#line 377 "MB02QD.f"
	i__1 = maxwrk, i__2 = (integer) dwork[mn + 1] + mn;
#line 377 "MB02QD.f"
	maxwrk = max(i__1,i__2);
#line 378 "MB02QD.f"
    } else {
#line 379 "MB02QD.f"
	mb03oy_(m, n, &a[a_offset], lda, rcond, svlmax, rank, &sval[1], &jpvt[
		1], &dwork[1], &dwork[mn + 1], info);
#line 381 "MB02QD.f"
    }

/*     Logically partition R = [ R11 R12 ] */
/*                             [  0  R22 ], */
/*     where R11 = R(1:RANK,1:RANK). */

/*     [R11,R12] = [ T11, 0 ] * Z. */

/*     Details of Householder transformations stored in DWORK(MN+1:2*MN). */
/*     Workspace need   3*min(M,N); */
/*               prefer 2*min(M,N)+min(M,N)*NB. */

#line 393 "MB02QD.f"
    if (*rank < *n) {
#line 394 "MB02QD.f"
	i__1 = *ldwork - (mn << 1);
#line 394 "MB02QD.f"
	dtzrzf_(rank, n, &a[a_offset], lda, &dwork[mn + 1], &dwork[(mn << 1) 
		+ 1], &i__1, info);
/* Computing MAX */
#line 396 "MB02QD.f"
	i__1 = maxwrk, i__2 = (integer) dwork[(mn << 1) + 1] + (mn << 1);
#line 396 "MB02QD.f"
	maxwrk = max(i__1,i__2);
#line 397 "MB02QD.f"
    }

#line 399 "MB02QD.f"
    if (*nrhs > 0) {

/*        B(1:M,1:NRHS) := Q' * B(1:M,1:NRHS). */

/*        Workspace: need   2*min(M,N)+NRHS; */
/*                   prefer   min(M,N)+NRHS*NB. */

#line 406 "MB02QD.f"
	i__1 = *ldwork - (mn << 1);
#line 406 "MB02QD.f"
	dormqr_("Left", "Transpose", m, nrhs, &mn, &a[a_offset], lda, &dwork[
		1], &b[b_offset], ldb, &dwork[(mn << 1) + 1], &i__1, info, (
		ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 409 "MB02QD.f"
	i__1 = maxwrk, i__2 = (integer) dwork[(mn << 1) + 1] + (mn << 1);
#line 409 "MB02QD.f"
	maxwrk = max(i__1,i__2);

/*        B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS). */

#line 413 "MB02QD.f"
	dtrsm_("Left", "Upper", "No transpose", "Non-unit", rank, nrhs, &
		c_b33, &a[a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);

#line 416 "MB02QD.f"
	if (*rank < *n) {

/*           Set B(RANK+1:N,1:NRHS). */

#line 420 "MB02QD.f"
	    if (leasts) {
#line 421 "MB02QD.f"
		i__1 = *n - *rank;
#line 421 "MB02QD.f"
		dlaset_("Full", &i__1, nrhs, &c_b17, &c_b17, &b[*rank + 1 + 
			b_dim1], ldb, (ftnlen)4);
#line 423 "MB02QD.f"
	    } else {
#line 424 "MB02QD.f"
		i__1 = *n - *rank;
#line 424 "MB02QD.f"
		i__2 = *n - *rank;
#line 424 "MB02QD.f"
		dlacpy_("Full", &i__1, nrhs, &y[1], &i__2, &b[*rank + 1 + 
			b_dim1], ldb, (ftnlen)4);
#line 426 "MB02QD.f"
	    }

/*           B(1:N,1:NRHS) := Z' * B(1:N,1:NRHS). */

/*           Workspace need   2*min(M,N)+NRHS; */
/*                     prefer 2*min(M,N)+NRHS*NB. */

#line 433 "MB02QD.f"
	    i__1 = *n - *rank;
#line 433 "MB02QD.f"
	    i__2 = *ldwork - (mn << 1);
#line 433 "MB02QD.f"
	    dormrz_("Left", "Transpose", n, nrhs, rank, &i__1, &a[a_offset], 
		    lda, &dwork[mn + 1], &b[b_offset], ldb, &dwork[(mn << 1) 
		    + 1], &i__2, info, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 436 "MB02QD.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[(mn << 1) + 1] + (mn << 1);
#line 436 "MB02QD.f"
	    maxwrk = max(i__1,i__2);
#line 437 "MB02QD.f"
	}

/*        Additional workspace: NRHS. */

/*        B(1:N,1:NRHS) := P * B(1:N,1:NRHS). */

#line 443 "MB02QD.f"
	i__1 = *nrhs;
#line 443 "MB02QD.f"
	for (j = 1; j <= i__1; ++j) {
#line 444 "MB02QD.f"
	    i__2 = *n;
#line 444 "MB02QD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 445 "MB02QD.f"
		dwork[(mn << 1) + i__] = 1.;
#line 446 "MB02QD.f"
/* L20: */
#line 446 "MB02QD.f"
	    }
#line 447 "MB02QD.f"
	    i__2 = *n;
#line 447 "MB02QD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 448 "MB02QD.f"
		if (dwork[(mn << 1) + i__] == 1.) {
#line 449 "MB02QD.f"
		    if (jpvt[i__] != i__) {
#line 450 "MB02QD.f"
			k = i__;
#line 451 "MB02QD.f"
			t1 = b[k + j * b_dim1];
#line 452 "MB02QD.f"
			t2 = b[jpvt[k] + j * b_dim1];
#line 453 "MB02QD.f"
L30:
#line 454 "MB02QD.f"
			b[jpvt[k] + j * b_dim1] = t1;
#line 455 "MB02QD.f"
			dwork[(mn << 1) + k] = 0.;
#line 456 "MB02QD.f"
			t1 = t2;
#line 457 "MB02QD.f"
			k = jpvt[k];
#line 458 "MB02QD.f"
			t2 = b[jpvt[k] + j * b_dim1];
#line 459 "MB02QD.f"
			if (jpvt[k] != i__) {
#line 459 "MB02QD.f"
			    goto L30;
#line 459 "MB02QD.f"
			}
#line 461 "MB02QD.f"
			b[i__ + j * b_dim1] = t1;
#line 462 "MB02QD.f"
			dwork[(mn << 1) + k] = 0.;
#line 463 "MB02QD.f"
		    }
#line 464 "MB02QD.f"
		}
#line 465 "MB02QD.f"
/* L40: */
#line 465 "MB02QD.f"
	    }
#line 466 "MB02QD.f"
/* L50: */
#line 466 "MB02QD.f"
	}

/*        Undo scaling for B. */

#line 470 "MB02QD.f"
	if (ibscl == 1) {
#line 471 "MB02QD.f"
	    dlascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], 
		    ldb, info, (ftnlen)1);
#line 473 "MB02QD.f"
	} else if (ibscl == 2) {
#line 474 "MB02QD.f"
	    dlascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], 
		    ldb, info, (ftnlen)1);
#line 476 "MB02QD.f"
	}
#line 477 "MB02QD.f"
    }

/*     Undo scaling for A. */

#line 481 "MB02QD.f"
    if (iascl == 1) {
#line 482 "MB02QD.f"
	if (*nrhs > 0) {
#line 482 "MB02QD.f"
	    dlascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], 
		    ldb, info, (ftnlen)1);
#line 482 "MB02QD.f"
	}
#line 485 "MB02QD.f"
	dlascl_("U", &c__0, &c__0, &smlnum, &anrm, rank, rank, &a[a_offset], 
		lda, info, (ftnlen)1);
#line 487 "MB02QD.f"
    } else if (iascl == 2) {
#line 488 "MB02QD.f"
	if (*nrhs > 0) {
#line 488 "MB02QD.f"
	    dlascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], 
		    ldb, info, (ftnlen)1);
#line 488 "MB02QD.f"
	}
#line 491 "MB02QD.f"
	dlascl_("U", &c__0, &c__0, &bignum, &anrm, rank, rank, &a[a_offset], 
		lda, info, (ftnlen)1);
#line 493 "MB02QD.f"
    }

#line 495 "MB02QD.f"
    for (i__ = mn + *rank; i__ >= 1; --i__) {
#line 496 "MB02QD.f"
	dwork[i__ + 1] = dwork[i__];
#line 497 "MB02QD.f"
/* L60: */
#line 497 "MB02QD.f"
    }

#line 499 "MB02QD.f"
    dwork[1] = (doublereal) maxwrk;
#line 500 "MB02QD.f"
    return 0;
/* *** Last line of MB02QD *** */
} /* mb02qd_ */

