#line 1 "SB02MT.f"
/* SB02MT.f -- translated by f2c (version 20100827).
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

#line 1 "SB02MT.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b28 = 1.;
static doublereal c_b31 = 0.;
static doublereal c_b37 = -1.;

/* Subroutine */ int sb02mt_(char *jobg, char *jobl, char *fact, char *uplo, 
	integer *n, integer *m, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *q, integer *ldq, doublereal *r__, integer *
	ldr, doublereal *l, integer *ldl, integer *ipiv, integer *oufact, 
	doublereal *g, integer *ldg, integer *iwork, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen jobg_len, ftnlen jobl_len, 
	ftnlen fact_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, g_dim1, g_offset, l_dim1, 
	    l_offset, q_dim1, q_offset, r_dim1, r_offset, i__1, i__2, i__3, 
	    i__4;

    /* Local variables */
    static integer i__, j;
    static doublereal eps;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static logical ljobg;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical ljobl;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal rcond;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static char trans[1];
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dsyrk_(
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal rnorm;
    extern doublereal dlamch_(char *, ftnlen);
    static logical lfacta, lfactc;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical lfactu;
    extern /* Subroutine */ int dpocon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen);
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int dpotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dsycon_(char *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen), dsytrf_(char *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static logical luplou;
    static integer wrkopt;
    extern /* Subroutine */ int dsytrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);


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

/*     To compute the following matrices */

/*                -1 */
/*         G = B*R  *B', */

/*         -          -1 */
/*         A = A - B*R  *L', */

/*         -          -1 */
/*         Q = Q - L*R  *L', */

/*     where A, B, Q, R, L, and G are N-by-N, N-by-M, N-by-N, M-by-M, */
/*     N-by-M, and N-by-N matrices, respectively, with Q, R and G */
/*     symmetric matrices. */

/*     When R is well-conditioned with respect to inversion, standard */
/*     algorithms for solving linear-quadratic optimization problems will */
/*     then also solve optimization problems with coupling weighting */
/*     matrix L. Moreover, a gain in efficiency is possible using matrix */
/*     G in the deflating subspace algorithms (see SLICOT Library routine */
/*     SB02OD). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBG    CHARACTER*1 */
/*             Specifies whether or not the matrix G is to be computed, */
/*             as follows: */
/*             = 'G':  Compute G; */
/*             = 'N':  Do not compute G. */

/*     JOBL    CHARACTER*1 */
/*             Specifies whether or not the matrix L is zero, as follows: */
/*             = 'Z':  L is zero; */
/*             = 'N':  L is nonzero. */

/*     FACT    CHARACTER*1 */
/*             Specifies how the matrix R is given (factored or not), as */
/*             follows: */
/*             = 'N':  Array R contains the matrix R; */
/*             = 'C':  Array R contains the Cholesky factor of R; */
/*             = 'U':  Array R contains the symmetric indefinite UdU' or */
/*                     LdL' factorization of R. */

/*     UPLO    CHARACTER*1 */
/*             Specifies which triangle of the matrices R and Q (if */
/*             JOBL = 'N') is stored, as follows: */
/*             = 'U':  Upper triangle is stored; */
/*             = 'L':  Lower triangle is stored. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, Q, and G, and the number of */
/*             rows of the matrices B and L.  N >= 0. */

/*     M       (input) INTEGER */
/*             The order of the matrix R, and the number of columns of */
/*             the matrices B and L.  M >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, if JOBL = 'N', the leading N-by-N part of this */
/*             array must contain the matrix A. */
/*             On exit, if JOBL = 'N', and INFO = 0, the leading N-by-N */
/*                                                    -          -1 */
/*             part of this array contains the matrix A = A - B*R  L'. */
/*             If JOBL = 'Z', this array is not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of array A. */
/*             LDA >= MAX(1,N) if JOBL = 'N'; */
/*             LDA >= 1        if JOBL = 'Z'. */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the matrix B. */
/*             On exit, if OUFACT = 1, and INFO = 0, the leading N-by-M */
/*                                                             -1 */
/*             part of this array contains the matrix B*chol(R)  . */
/*             On exit, B is unchanged if OUFACT = 2 (hence also when */
/*             FACT = 'U'). */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             On entry, if JOBL = 'N', the leading N-by-N upper */
/*             triangular part (if UPLO = 'U') or lower triangular part */
/*             (if UPLO = 'L') of this array must contain the upper */
/*             triangular part or lower triangular part, respectively, of */
/*             the symmetric matrix Q. The stricly lower triangular part */
/*             (if UPLO = 'U') or stricly upper triangular part (if */
/*             UPLO = 'L') is not referenced. */
/*             On exit, if JOBL = 'N' and INFO = 0, the leading N-by-N */
/*             upper triangular part (if UPLO = 'U') or lower triangular */
/*             part (if UPLO = 'L') of this array contains the upper */
/*             triangular part or lower triangular part, respectively, of */
/*                                  -          -1 */
/*             the symmetric matrix Q = Q - L*R  *L'. */
/*             If JOBL = 'Z', this array is not referenced. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. */
/*             LDQ >= MAX(1,N) if JOBL = 'N'; */
/*             LDQ >= 1        if JOBL = 'Z'. */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR,M) */
/*             On entry, if FACT = 'N', the leading M-by-M upper */
/*             triangular part (if UPLO = 'U') or lower triangular part */
/*             (if UPLO = 'L') of this array must contain the upper */
/*             triangular part or lower triangular part, respectively, */
/*             of the symmetric input weighting matrix R. */
/*             On entry, if FACT = 'C', the leading M-by-M upper */
/*             triangular part (if UPLO = 'U') or lower triangular part */
/*             (if UPLO = 'L') of this array must contain the Cholesky */
/*             factor of the positive definite input weighting matrix R */
/*             (as produced by LAPACK routine DPOTRF). */
/*             On entry, if FACT = 'U', the leading M-by-M upper */
/*             triangular part (if UPLO = 'U') or lower triangular part */
/*             (if UPLO = 'L') of this array must contain the factors of */
/*             the UdU' or LdL' factorization, respectively, of the */
/*             symmetric indefinite input weighting matrix R (as produced */
/*             by LAPACK routine DSYTRF). */
/*             If FACT = 'N', the stricly lower triangular part (if UPLO */
/*             = 'U') or stricly upper triangular part (if UPLO = 'L') of */
/*             this array is used as workspace. */
/*             On exit, if OUFACT = 1, and INFO = 0 (or INFO = M+1), */
/*             the leading M-by-M upper triangular part (if UPLO = 'U') */
/*             or lower triangular part (if UPLO = 'L') of this array */
/*             contains the Cholesky factor of the given input weighting */
/*             matrix. */
/*             On exit, if OUFACT = 2, and INFO = 0 (or INFO = M+1), */
/*             the leading M-by-M upper triangular part (if UPLO = 'U') */
/*             or lower triangular part (if UPLO = 'L') of this array */
/*             contains the factors of the UdU' or LdL' factorization, */
/*             respectively, of the given input weighting matrix. */
/*             On exit R is unchanged if FACT = 'C' or 'U'. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,M). */

/*     L       (input/output) DOUBLE PRECISION array, dimension (LDL,M) */
/*             On entry, if JOBL = 'N', the leading N-by-M part of this */
/*             array must contain the matrix L. */
/*             On exit, if JOBL = 'N', OUFACT = 1, and INFO = 0, the */
/*             leading N-by-M part of this array contains the matrix */
/*                      -1 */
/*             L*chol(R)  . */
/*             On exit, L is unchanged if OUFACT = 2 (hence also when */
/*             FACT = 'U'). */
/*             L is not referenced if JOBL = 'Z'. */

/*     LDL     INTEGER */
/*             The leading dimension of array L. */
/*             LDL >= MAX(1,N) if JOBL = 'N'; */
/*             LDL >= 1        if JOBL = 'Z'. */

/*     IPIV    (input/output) INTEGER array, dimension (M) */
/*             On entry, if FACT = 'U', this array must contain details */
/*             of the interchanges performed and the block structure of */
/*             the d factor in the UdU' or LdL' factorization of matrix R */
/*             (as produced by LAPACK routine DSYTRF). */
/*             On exit, if OUFACT = 2, this array contains details of */
/*             the interchanges performed and the block structure of the */
/*             d factor in the UdU' or LdL' factorization of matrix R, */
/*             as produced by LAPACK routine DSYTRF. */
/*             This array is not referenced if FACT = 'C'. */

/*     OUFACT  (output) INTEGER */
/*             Information about the factorization finally used. */
/*             OUFACT = 1:  Cholesky factorization of R has been used; */
/*             OUFACT = 2:  UdU' (if UPLO = 'U') or LdL' (if UPLO = 'L') */
/*                          factorization of R has been used. */

/*     G       (output) DOUBLE PRECISION array, dimension (LDG,N) */
/*             If JOBG = 'G', and INFO = 0, the leading N-by-N upper */
/*             triangular part (if UPLO = 'U') or lower triangular part */
/*             (if UPLO = 'L') of this array contains the upper */
/*             triangular part (if UPLO = 'U') or lower triangular part */
/*                                                                 -1 */
/*             (if UPLO = 'L'), respectively, of the matrix G = B*R  B'. */
/*             If JOBG = 'N', this array is not referenced. */

/*     LDG     INTEGER */
/*             The leading dimension of array G. */
/*             LDG >= MAX(1,N) if JOBG = 'G', */
/*             LDG >= 1        if JOBG = 'N'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK; if FACT = 'N', DWORK(2) contains the reciprocal */
/*             condition number of the given matrix R. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 1              if FACT = 'C'; */
/*             LDWORK >= MAX(2,3*M,N*M) if FACT = 'N'; */
/*             LDWORK >= MAX(1,N*M)     if FACT = 'U'. */
/*             For optimum performance LDWORK should be larger than 3*M, */
/*             if FACT = 'N'. */
/*             The N*M workspace is not needed for FACT = 'N', if matrix */
/*             R is positive definite. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = i:  if the i-th element (1 <= i <= M) of the d factor is */
/*                   exactly zero; the UdU' (or LdL') factorization has */
/*                   been completed, but the block diagonal matrix d is */
/*                   exactly singular; */
/*             = M+1:  if the matrix R is numerically singular. */

/*     METHOD */
/*                            -     - */
/*     The matrices G, and/or A and Q are evaluated using the given or */
/*     computed symmetric factorization of R. */

/*     NUMERICAL ASPECTS */

/*     The routine should not be used when R is ill-conditioned. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Sep. 1997. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Algebraic Riccati equation, closed loop system, continuous-time */
/*     system, discrete-time system, optimal regulator, Schur form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 296 "SB02MT.f"
    /* Parameter adjustments */
#line 296 "SB02MT.f"
    a_dim1 = *lda;
#line 296 "SB02MT.f"
    a_offset = 1 + a_dim1;
#line 296 "SB02MT.f"
    a -= a_offset;
#line 296 "SB02MT.f"
    b_dim1 = *ldb;
#line 296 "SB02MT.f"
    b_offset = 1 + b_dim1;
#line 296 "SB02MT.f"
    b -= b_offset;
#line 296 "SB02MT.f"
    q_dim1 = *ldq;
#line 296 "SB02MT.f"
    q_offset = 1 + q_dim1;
#line 296 "SB02MT.f"
    q -= q_offset;
#line 296 "SB02MT.f"
    r_dim1 = *ldr;
#line 296 "SB02MT.f"
    r_offset = 1 + r_dim1;
#line 296 "SB02MT.f"
    r__ -= r_offset;
#line 296 "SB02MT.f"
    l_dim1 = *ldl;
#line 296 "SB02MT.f"
    l_offset = 1 + l_dim1;
#line 296 "SB02MT.f"
    l -= l_offset;
#line 296 "SB02MT.f"
    --ipiv;
#line 296 "SB02MT.f"
    g_dim1 = *ldg;
#line 296 "SB02MT.f"
    g_offset = 1 + g_dim1;
#line 296 "SB02MT.f"
    g -= g_offset;
#line 296 "SB02MT.f"
    --iwork;
#line 296 "SB02MT.f"
    --dwork;
#line 296 "SB02MT.f"

#line 296 "SB02MT.f"
    /* Function Body */
#line 296 "SB02MT.f"
    *info = 0;
#line 297 "SB02MT.f"
    ljobg = lsame_(jobg, "G", (ftnlen)1, (ftnlen)1);
#line 298 "SB02MT.f"
    ljobl = lsame_(jobl, "N", (ftnlen)1, (ftnlen)1);
#line 299 "SB02MT.f"
    lfactc = lsame_(fact, "C", (ftnlen)1, (ftnlen)1);
#line 300 "SB02MT.f"
    lfactu = lsame_(fact, "U", (ftnlen)1, (ftnlen)1);
#line 301 "SB02MT.f"
    luplou = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 302 "SB02MT.f"
    lfacta = lfactc || lfactu;

/*     Test the input scalar arguments. */

#line 306 "SB02MT.f"
    if (! ljobg && ! lsame_(jobg, "N", (ftnlen)1, (ftnlen)1)) {
#line 307 "SB02MT.f"
	*info = -1;
#line 308 "SB02MT.f"
    } else if (! ljobl && ! lsame_(jobl, "Z", (ftnlen)1, (ftnlen)1)) {
#line 309 "SB02MT.f"
	*info = -2;
#line 310 "SB02MT.f"
    } else if (! lfacta && ! lsame_(fact, "N", (ftnlen)1, (ftnlen)1)) {
#line 311 "SB02MT.f"
	*info = -3;
#line 312 "SB02MT.f"
    } else if (! luplou && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 313 "SB02MT.f"
	*info = -4;
#line 314 "SB02MT.f"
    } else if (*n < 0) {
#line 315 "SB02MT.f"
	*info = -5;
#line 316 "SB02MT.f"
    } else if (*m < 0) {
#line 317 "SB02MT.f"
	*info = -6;
#line 318 "SB02MT.f"
    } else if (*lda < 1 || ljobl && *lda < *n) {
#line 319 "SB02MT.f"
	*info = -8;
#line 320 "SB02MT.f"
    } else if (*ldb < max(1,*n)) {
#line 321 "SB02MT.f"
	*info = -10;
#line 322 "SB02MT.f"
    } else if (*ldq < 1 || ljobl && *ldq < *n) {
#line 323 "SB02MT.f"
	*info = -12;
#line 324 "SB02MT.f"
    } else if (*ldr < max(1,*m)) {
#line 325 "SB02MT.f"
	*info = -14;
#line 326 "SB02MT.f"
    } else if (*ldl < 1 || ljobl && *ldl < *n) {
#line 327 "SB02MT.f"
	*info = -16;
#line 328 "SB02MT.f"
    } else if (*ldg < 1 || ljobg && *ldg < *n) {
#line 329 "SB02MT.f"
	*info = -20;
#line 330 "SB02MT.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 330 "SB02MT.f"
	i__1 = 1, i__2 = *n * *m;
/* Computing MAX */
#line 330 "SB02MT.f"
	i__3 = 2, i__4 = *n * *m, i__3 = max(i__3,i__4), i__4 = *m * 3;
#line 330 "SB02MT.f"
	if (lfactc && *ldwork < 1 || lfactu && *ldwork < max(i__1,i__2) || ! 
		lfacta && *ldwork < max(i__3,i__4)) {
#line 333 "SB02MT.f"
	    *info = -23;
#line 334 "SB02MT.f"
	}
#line 334 "SB02MT.f"
    }

#line 336 "SB02MT.f"
    if (*info != 0) {

/*        Error return. */

#line 340 "SB02MT.f"
	i__1 = -(*info);
#line 340 "SB02MT.f"
	xerbla_("SB02MT", &i__1, (ftnlen)6);
#line 341 "SB02MT.f"
	return 0;
#line 342 "SB02MT.f"
    }

#line 344 "SB02MT.f"
    if (lfactc) {
#line 345 "SB02MT.f"
	*oufact = 1;
#line 346 "SB02MT.f"
    } else if (lfactu) {
#line 347 "SB02MT.f"
	*oufact = 2;
#line 348 "SB02MT.f"
    }

/*     Quick return if possible. */

#line 352 "SB02MT.f"
    if (*n == 0 || *m == 0 || ! (ljobl || ljobg)) {
#line 353 "SB02MT.f"
	dwork[1] = 1.;
#line 354 "SB02MT.f"
	if (! lfacta) {
#line 354 "SB02MT.f"
	    dwork[2] = 1.;
#line 354 "SB02MT.f"
	}
#line 355 "SB02MT.f"
	return 0;
#line 356 "SB02MT.f"
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of workspace needed at that point in the code, */
/*     as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 364 "SB02MT.f"
    wrkopt = 1;

/*     Set relative machine precision. */

#line 368 "SB02MT.f"
    eps = dlamch_("Epsilon", (ftnlen)7);

#line 370 "SB02MT.f"
    if (! lfacta) {

/*        Compute the norm of the matrix R, which is not factored. */
/*        Then save the given triangle of R in the other strict triangle */
/*        and the diagonal in the workspace, and try Cholesky */
/*        factorization. */
/*        Workspace: need M. */

#line 378 "SB02MT.f"
	rnorm = dlansy_("1-norm", uplo, m, &r__[r_offset], ldr, &dwork[1], (
		ftnlen)6, (ftnlen)1);
#line 379 "SB02MT.f"
	i__1 = *ldr + 1;
#line 379 "SB02MT.f"
	dcopy_(m, &r__[r_offset], &i__1, &dwork[1], &c__1);
#line 380 "SB02MT.f"
	if (luplou) {

#line 382 "SB02MT.f"
	    i__1 = *m;
#line 382 "SB02MT.f"
	    for (j = 2; j <= i__1; ++j) {
#line 383 "SB02MT.f"
		i__2 = j - 1;
#line 383 "SB02MT.f"
		dcopy_(&i__2, &r__[j * r_dim1 + 1], &c__1, &r__[j + r_dim1], 
			ldr);
#line 384 "SB02MT.f"
/* L20: */
#line 384 "SB02MT.f"
	    }

#line 386 "SB02MT.f"
	} else {

#line 388 "SB02MT.f"
	    i__1 = *m;
#line 388 "SB02MT.f"
	    for (j = 2; j <= i__1; ++j) {
#line 389 "SB02MT.f"
		i__2 = j - 1;
#line 389 "SB02MT.f"
		dcopy_(&i__2, &r__[j + r_dim1], ldr, &r__[j * r_dim1 + 1], &
			c__1);
#line 390 "SB02MT.f"
/* L40: */
#line 390 "SB02MT.f"
	    }

#line 392 "SB02MT.f"
	}
#line 393 "SB02MT.f"
	dpotrf_(uplo, m, &r__[r_offset], ldr, info, (ftnlen)1);
#line 394 "SB02MT.f"
	if (*info == 0) {

/*           Compute the reciprocal of the condition number of R. */
/*           Workspace: need 3*M. */

#line 399 "SB02MT.f"
	    dpocon_(uplo, m, &r__[r_offset], ldr, &rnorm, &rcond, &dwork[1], &
		    iwork[1], info, (ftnlen)1);

/*           Return if the matrix is singular to working precision. */

#line 404 "SB02MT.f"
	    *oufact = 1;
#line 405 "SB02MT.f"
	    dwork[2] = rcond;
#line 406 "SB02MT.f"
	    if (rcond < eps) {
#line 407 "SB02MT.f"
		*info = *m + 1;
#line 408 "SB02MT.f"
		return 0;
#line 409 "SB02MT.f"
	    }
/* Computing MAX */
#line 410 "SB02MT.f"
	    i__1 = wrkopt, i__2 = *m * 3;
#line 410 "SB02MT.f"
	    wrkopt = max(i__1,i__2);
#line 411 "SB02MT.f"
	} else {

/*           Use UdU' or LdL' factorization, first restoring the saved */
/*           triangle. */

#line 416 "SB02MT.f"
	    i__1 = *ldr + 1;
#line 416 "SB02MT.f"
	    dcopy_(m, &dwork[1], &c__1, &r__[r_offset], &i__1);
#line 417 "SB02MT.f"
	    if (luplou) {

#line 419 "SB02MT.f"
		i__1 = *m;
#line 419 "SB02MT.f"
		for (j = 2; j <= i__1; ++j) {
#line 420 "SB02MT.f"
		    i__2 = j - 1;
#line 420 "SB02MT.f"
		    dcopy_(&i__2, &r__[j + r_dim1], ldr, &r__[j * r_dim1 + 1],
			     &c__1);
#line 421 "SB02MT.f"
/* L60: */
#line 421 "SB02MT.f"
		}

#line 423 "SB02MT.f"
	    } else {

#line 425 "SB02MT.f"
		i__1 = *m;
#line 425 "SB02MT.f"
		for (j = 2; j <= i__1; ++j) {
#line 426 "SB02MT.f"
		    i__2 = j - 1;
#line 426 "SB02MT.f"
		    dcopy_(&i__2, &r__[j * r_dim1 + 1], &c__1, &r__[j + 
			    r_dim1], ldr);
#line 427 "SB02MT.f"
/* L80: */
#line 427 "SB02MT.f"
		}

#line 429 "SB02MT.f"
	    }

/*           Compute the UdU' or LdL' factorization. */
/*           Workspace: need   1, */
/*                      prefer M*NB. */

#line 435 "SB02MT.f"
	    dsytrf_(uplo, m, &r__[r_offset], ldr, &ipiv[1], &dwork[1], ldwork,
		     info, (ftnlen)1);
#line 436 "SB02MT.f"
	    *oufact = 2;
#line 437 "SB02MT.f"
	    if (*info > 0) {
#line 438 "SB02MT.f"
		dwork[2] = 1.;
#line 439 "SB02MT.f"
		return 0;
#line 440 "SB02MT.f"
	    }
/* Computing MAX */
#line 441 "SB02MT.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 441 "SB02MT.f"
	    wrkopt = max(i__1,i__2);

/*           Compute the reciprocal of the condition number of R. */
/*           Workspace: need 2*M. */

#line 446 "SB02MT.f"
	    dsycon_(uplo, m, &r__[r_offset], ldr, &ipiv[1], &rnorm, &rcond, &
		    dwork[1], &iwork[1], info, (ftnlen)1);

/*           Return if the matrix is singular to working precision. */

#line 451 "SB02MT.f"
	    dwork[2] = rcond;
#line 452 "SB02MT.f"
	    if (rcond < eps) {
#line 453 "SB02MT.f"
		*info = *m + 1;
#line 454 "SB02MT.f"
		return 0;
#line 455 "SB02MT.f"
	    }
#line 456 "SB02MT.f"
	}
#line 457 "SB02MT.f"
    }

#line 459 "SB02MT.f"
    if (*oufact == 1) {

/*        Solve positive definite linear system(s). */

#line 463 "SB02MT.f"
	if (luplou) {
#line 464 "SB02MT.f"
	    *(unsigned char *)trans = 'N';
#line 465 "SB02MT.f"
	} else {
#line 466 "SB02MT.f"
	    *(unsigned char *)trans = 'T';
#line 467 "SB02MT.f"
	}

/*        Solve the system X*U = B, overwriting B with X. */

#line 471 "SB02MT.f"
	dtrsm_("Right", uplo, trans, "Non-unit", n, m, &c_b28, &r__[r_offset],
		 ldr, &b[b_offset], ldb, (ftnlen)5, (ftnlen)1, (ftnlen)1, (
		ftnlen)8);

#line 474 "SB02MT.f"
	if (ljobg) {
/*                                      -1 */
/*           Compute the matrix  G = B*R  *B', multiplying X*X' in G. */

#line 478 "SB02MT.f"
	    dsyrk_(uplo, "No transpose", n, m, &c_b28, &b[b_offset], ldb, &
		    c_b31, &g[g_offset], ldg, (ftnlen)1, (ftnlen)12);
#line 480 "SB02MT.f"
	}

#line 482 "SB02MT.f"
	if (ljobl) {

/*           Update matrices A and Q. */

/*           Solve the system Y*U = L, overwriting L with Y. */

#line 488 "SB02MT.f"
	    dtrsm_("Right", uplo, trans, "Non-unit", n, m, &c_b28, &r__[
		    r_offset], ldr, &l[l_offset], ldl, (ftnlen)5, (ftnlen)1, (
		    ftnlen)1, (ftnlen)8);

/*           Compute A <- A - X*Y'. */

#line 493 "SB02MT.f"
	    dgemm_("No transpose", "Transpose", n, n, m, &c_b37, &b[b_offset],
		     ldb, &l[l_offset], ldl, &c_b28, &a[a_offset], lda, (
		    ftnlen)12, (ftnlen)9);

/*           Compute Q <- Q - Y*Y'. */

#line 498 "SB02MT.f"
	    dsyrk_(uplo, "No transpose", n, m, &c_b37, &l[l_offset], ldl, &
		    c_b28, &q[q_offset], ldq, (ftnlen)1, (ftnlen)12);
#line 500 "SB02MT.f"
	}
#line 501 "SB02MT.f"
    } else {

/*        Solve indefinite linear system(s). */

/*        Solve the system UdU'*X = B' (or LdL'*X = B'). */
/*        Workspace: need N*M. */

#line 508 "SB02MT.f"
	i__1 = *m;
#line 508 "SB02MT.f"
	for (j = 1; j <= i__1; ++j) {
#line 509 "SB02MT.f"
	    dcopy_(n, &b[j * b_dim1 + 1], &c__1, &dwork[j], m);
#line 510 "SB02MT.f"
/* L100: */
#line 510 "SB02MT.f"
	}

#line 512 "SB02MT.f"
	dsytrs_(uplo, m, n, &r__[r_offset], ldr, &ipiv[1], &dwork[1], m, info,
		 (ftnlen)1);

#line 514 "SB02MT.f"
	if (ljobg) {
/*                                                    -1 */
/*           Compute a triangle of the matrix  G = B*R  *B' = B*X. */

#line 518 "SB02MT.f"
	    if (luplou) {
#line 519 "SB02MT.f"
		i__ = 1;

#line 521 "SB02MT.f"
		i__1 = *n;
#line 521 "SB02MT.f"
		for (j = 1; j <= i__1; ++j) {
#line 522 "SB02MT.f"
		    dgemv_("No transpose", &j, m, &c_b28, &b[b_offset], ldb, &
			    dwork[i__], &c__1, &c_b31, &g[j * g_dim1 + 1], &
			    c__1, (ftnlen)12);
#line 524 "SB02MT.f"
		    i__ += *m;
#line 525 "SB02MT.f"
/* L120: */
#line 525 "SB02MT.f"
		}

#line 527 "SB02MT.f"
	    } else {

#line 529 "SB02MT.f"
		i__1 = *n;
#line 529 "SB02MT.f"
		for (j = 1; j <= i__1; ++j) {
#line 530 "SB02MT.f"
		    dgemv_("Transpose", m, &j, &c_b28, &dwork[1], m, &b[j + 
			    b_dim1], ldb, &c_b31, &g[j + g_dim1], ldg, (
			    ftnlen)9);
#line 532 "SB02MT.f"
/* L140: */
#line 532 "SB02MT.f"
		}

#line 534 "SB02MT.f"
	    }
#line 535 "SB02MT.f"
	}

#line 537 "SB02MT.f"
	if (ljobl) {

/*           Update matrices A and Q. */

/*           Solve the system UdU'*Y = L' (or LdL'*Y = L'). */

#line 543 "SB02MT.f"
	    i__1 = *m;
#line 543 "SB02MT.f"
	    for (j = 1; j <= i__1; ++j) {
#line 544 "SB02MT.f"
		dcopy_(n, &l[j * l_dim1 + 1], &c__1, &dwork[j], m);
#line 545 "SB02MT.f"
/* L160: */
#line 545 "SB02MT.f"
	    }

#line 547 "SB02MT.f"
	    dsytrs_(uplo, m, n, &r__[r_offset], ldr, &ipiv[1], &dwork[1], m, 
		    info, (ftnlen)1);

/*           A <- A - B*Y. */

#line 551 "SB02MT.f"
	    dgemm_("No transpose", "No transpose", n, n, m, &c_b37, &b[
		    b_offset], ldb, &dwork[1], m, &c_b28, &a[a_offset], lda, (
		    ftnlen)12, (ftnlen)12);
/*                                            -          -1 */
/*           Compute a triangle of the matrix Q = Q - L*R  *L' = Q - L*Y. */

#line 556 "SB02MT.f"
	    if (luplou) {
#line 557 "SB02MT.f"
		i__ = 1;

#line 559 "SB02MT.f"
		i__1 = *n;
#line 559 "SB02MT.f"
		for (j = 1; j <= i__1; ++j) {
#line 560 "SB02MT.f"
		    dgemv_("No transpose", &j, m, &c_b37, &l[l_offset], ldl, &
			    dwork[i__], &c__1, &c_b28, &q[j * q_dim1 + 1], &
			    c__1, (ftnlen)12);
#line 562 "SB02MT.f"
		    i__ += *m;
#line 563 "SB02MT.f"
/* L180: */
#line 563 "SB02MT.f"
		}

#line 565 "SB02MT.f"
	    } else {

#line 567 "SB02MT.f"
		i__1 = *n;
#line 567 "SB02MT.f"
		for (j = 1; j <= i__1; ++j) {
#line 568 "SB02MT.f"
		    dgemv_("Transpose", m, &j, &c_b37, &dwork[1], m, &l[j + 
			    l_dim1], ldl, &c_b28, &q[j + q_dim1], ldq, (
			    ftnlen)9);
#line 570 "SB02MT.f"
/* L200: */
#line 570 "SB02MT.f"
		}

#line 572 "SB02MT.f"
	    }
#line 573 "SB02MT.f"
	}
#line 574 "SB02MT.f"
    }

#line 576 "SB02MT.f"
    dwork[1] = (doublereal) wrkopt;
#line 577 "SB02MT.f"
    if (! lfacta) {
#line 577 "SB02MT.f"
	dwork[2] = rcond;
#line 577 "SB02MT.f"
    }

/* *** Last line of SB02MT *** */
#line 580 "SB02MT.f"
    return 0;
} /* sb02mt_ */

