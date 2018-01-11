#line 1 "MB02PD.f"
/* MB02PD.f -- translated by f2c (version 20100827).
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

#line 1 "MB02PD.f"
/* Subroutine */ int mb02pd_(char *fact, char *trans, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	integer *ipiv, char *equed, doublereal *r__, doublereal *c__, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	rcond, doublereal *ferr, doublereal *berr, integer *iwork, doublereal 
	*dwork, integer *info, ftnlen fact_len, ftnlen trans_len, ftnlen 
	equed_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, 
	    x_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal amax;
    static char norm[1];
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal rcmin, rcmax, anorm;
    static logical equil;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlaqge_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, char *, ftnlen), dgecon_(char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, ftnlen);
    static doublereal colcnd;
    static logical nofact;
    extern /* Subroutine */ int dgeequ_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *), dgerfs_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen), 
	    dgetrf_(integer *, integer *, doublereal *, integer *, integer *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    static doublereal bignum;
    extern doublereal dlantr_(char *, char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    static integer infequ;
    static logical colequ;
    extern /* Subroutine */ int dgetrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static doublereal rowcnd;
    static logical notran;
    static doublereal smlnum;
    static logical rowequ;
    static doublereal rpvgrw;


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

/*     To solve (if well-conditioned) the matrix equations */

/*        op( A )*X = B, */

/*     where X and B are N-by-NRHS matrices, A is an N-by-N matrix and */
/*     op( A ) is one of */

/*        op( A ) = A   or   op( A ) = A'. */

/*     Error bounds on the solution and a condition estimate are also */
/*     provided. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     FACT    CHARACTER*1 */
/*             Specifies whether or not the factored form of the matrix A */
/*             is supplied on entry, and if not, whether the matrix A */
/*             should be equilibrated before it is factored. */
/*             = 'F':  On entry, AF and IPIV contain the factored form */
/*                     of A. If EQUED is not 'N', the matrix A has been */
/*                     equilibrated with scaling factors given by R */
/*                     and C. A, AF, and IPIV are not modified. */
/*             = 'N':  The matrix A will be copied to AF and factored. */
/*             = 'E':  The matrix A will be equilibrated if necessary, */
/*                     then copied to AF and factored. */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of the system of equations as follows: */
/*             = 'N':  A * X = B     (No transpose); */
/*             = 'T':  A**T * X = B  (Transpose); */
/*             = 'C':  A**H * X = B  (Transpose). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of linear equations, i.e., the order of the */
/*             matrix A.  N >= 0. */

/*     NRHS    (input) INTEGER */
/*             The number of right hand sides, i.e., the number of */
/*             columns of the matrices B and X.  NRHS >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A.  If FACT = 'F' and EQUED is not 'N', */
/*             then A must have been equilibrated by the scaling factors */
/*             in R and/or C.  A is not modified if FACT = 'F' or 'N', */
/*             or if FACT = 'E' and EQUED = 'N' on exit. */
/*             On exit, if EQUED .NE. 'N', the leading N-by-N part of */
/*             this array contains the matrix A scaled as follows: */
/*             EQUED = 'R':  A := diag(R) * A; */
/*             EQUED = 'C':  A := A * diag(C); */
/*             EQUED = 'B':  A := diag(R) * A * diag(C). */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     AF      (input or output) DOUBLE PRECISION array, dimension */
/*             (LDAF,N) */
/*             If FACT = 'F', then AF is an input argument and on entry */
/*             the leading N-by-N part of this array must contain the */
/*             factors L and U from the factorization A = P*L*U as */
/*             computed by DGETRF.  If EQUED .NE. 'N', then AF is the */
/*             factored form of the equilibrated matrix A. */
/*             If FACT = 'N', then AF is an output argument and on exit */
/*             the leading N-by-N part of this array contains the factors */
/*             L and U from the factorization A = P*L*U of the original */
/*             matrix A. */
/*             If FACT = 'E', then AF is an output argument and on exit */
/*             the leading N-by-N part of this array contains the factors */
/*             L and U from the factorization A = P*L*U of the */
/*             equilibrated matrix A (see the description of A for the */
/*             form of the equilibrated matrix). */

/*     LDAF    (input) INTEGER */
/*             The leading dimension of the array AF.  LDAF >= max(1,N). */

/*     IPIV    (input or output) INTEGER array, dimension (N) */
/*             If FACT = 'F', then IPIV is an input argument and on entry */
/*             it must contain the pivot indices from the factorization */
/*             A = P*L*U as computed by DGETRF; row i of the matrix was */
/*             interchanged with row IPIV(i). */
/*             If FACT = 'N', then IPIV is an output argument and on exit */
/*             it contains the pivot indices from the factorization */
/*             A = P*L*U of the original matrix A. */
/*             If FACT = 'E', then IPIV is an output argument and on exit */
/*             it contains the pivot indices from the factorization */
/*             A = P*L*U of the equilibrated matrix A. */

/*     EQUED   (input or output) CHARACTER*1 */
/*             Specifies the form of equilibration that was done as */
/*             follows: */
/*             = 'N':  No equilibration (always true if FACT = 'N'); */
/*             = 'R':  Row equilibration, i.e., A has been premultiplied */
/*                     by diag(R); */
/*             = 'C':  Column equilibration, i.e., A has been */
/*                     postmultiplied by diag(C); */
/*             = 'B':  Both row and column equilibration, i.e., A has */
/*                     been replaced by diag(R) * A * diag(C). */
/*             EQUED is an input argument if FACT = 'F'; otherwise, it is */
/*             an output argument. */

/*     R       (input or output) DOUBLE PRECISION array, dimension (N) */
/*             The row scale factors for A.  If EQUED = 'R' or 'B', A is */
/*             multiplied on the left by diag(R); if EQUED = 'N' or 'C', */
/*             R is not accessed.  R is an input argument if FACT = 'F'; */
/*             otherwise, R is an output argument.  If FACT = 'F' and */
/*             EQUED = 'R' or 'B', each element of R must be positive. */

/*     C       (input or output) DOUBLE PRECISION array, dimension (N) */
/*             The column scale factors for A.  If EQUED = 'C' or 'B', */
/*             A is multiplied on the right by diag(C); if EQUED = 'N' */
/*             or 'R', C is not accessed.  C is an input argument if */
/*             FACT = 'F'; otherwise, C is an output argument.  If */
/*             FACT = 'F' and EQUED = 'C' or 'B', each element of C must */
/*             be positive. */

/*     B       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDB,NRHS) */
/*             On entry, the leading N-by-NRHS part of this array must */
/*             contain the right-hand side matrix B. */
/*             On exit, */
/*             if EQUED = 'N', B is not modified; */
/*             if TRANS = 'N' and EQUED = 'R' or 'B', the leading */
/*             N-by-NRHS part of this array contains diag(R)*B; */
/*             if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', the leading */
/*             N-by-NRHS part of this array contains diag(C)*B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,NRHS) */
/*             If INFO = 0 or INFO = N+1, the leading N-by-NRHS part of */
/*             this array contains the solution matrix X to the original */
/*             system of equations.  Note that A and B are modified on */
/*             exit if EQUED .NE. 'N', and the solution to the */
/*             equilibrated system is inv(diag(C))*X if TRANS = 'N' and */
/*             EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or */
/*             'C' and EQUED = 'R' or 'B'. */

/*     LDX     (input) INTEGER */
/*             The leading dimension of the array X.  LDX >= max(1,N). */

/*     RCOND   (output) DOUBLE PRECISION */
/*             The estimate of the reciprocal condition number of the */
/*             matrix A after equilibration (if done).  If RCOND is less */
/*             than the machine precision (in particular, if RCOND = 0), */
/*             the matrix is singular to working precision.  This */
/*             condition is indicated by a return code of INFO > 0. */
/*             For efficiency reasons, RCOND is computed only when the */
/*             matrix A is factored, i.e., for FACT = 'N' or 'E'.  For */
/*             FACT = 'F', RCOND is not used, but it is assumed that it */
/*             has been computed and checked before the routine call. */

/*     FERR    (output) DOUBLE PRECISION array, dimension (NRHS) */
/*             The estimated forward error bound for each solution vector */
/*             X(j) (the j-th column of the solution matrix X). */
/*             If XTRUE is the true solution corresponding to X(j), */
/*             FERR(j) is an estimated upper bound for the magnitude of */
/*             the largest element in (X(j) - XTRUE) divided by the */
/*             magnitude of the largest element in X(j).  The estimate */
/*             is as reliable as the estimate for RCOND, and is almost */
/*             always a slight overestimate of the true error. */

/*     BERR    (output) DOUBLE PRECISION array, dimension (NRHS) */
/*             The componentwise relative backward error of each solution */
/*             vector X(j) (i.e., the smallest relative change in */
/*             any element of A or B that makes X(j) an exact solution). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     DWORK   DOUBLE PRECISION array, dimension (4*N) */
/*             On exit, DWORK(1) contains the reciprocal pivot growth */
/*             factor norm(A)/norm(U). The "max absolute element" norm is */
/*             used. If DWORK(1) is much less than 1, then the stability */
/*             of the LU factorization of the (equilibrated) matrix A */
/*             could be poor. This also means that the solution X, */
/*             condition estimator RCOND, and forward error bound FERR */
/*             could be unreliable. If factorization fails with */
/*             0 < INFO <= N, then DWORK(1) contains the reciprocal pivot */
/*             growth factor for the leading INFO columns of A. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, and i is */
/*                   <= N:  U(i,i) is exactly zero.  The factorization */
/*                          has been completed, but the factor U is */
/*                          exactly singular, so the solution and error */
/*                          bounds could not be computed. RCOND = 0 is */
/*                          returned. */
/*                   = N+1: U is nonsingular, but RCOND is less than */
/*                          machine precision, meaning that the matrix is */
/*                          singular to working precision.  Nevertheless, */
/*                          the solution and error bounds are computed */
/*                          because there are a number of situations */
/*                          where the computed solution can be more */
/*                          accurate than the value of RCOND would */
/*                          suggest. */
/*             The positive values for INFO are set only when the */
/*             matrix A is factored, i.e., for FACT = 'N' or 'E'. */

/*     METHOD */

/*     The following steps are performed: */

/*     1. If FACT = 'E', real scaling factors are computed to equilibrate */
/*        the system: */

/*        TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B */
/*        TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B */
/*        TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B */

/*        Whether or not the system will be equilibrated depends on the */
/*        scaling of the matrix A, but if equilibration is used, A is */
/*        overwritten by diag(R)*A*diag(C) and B by diag(R)*B */
/*        (if TRANS='N') or diag(C)*B (if TRANS = 'T' or 'C'). */

/*     2. If FACT = 'N' or 'E', the LU decomposition is used to factor */
/*        the matrix A (after equilibration if FACT = 'E') as */
/*           A = P * L * U, */
/*        where P is a permutation matrix, L is a unit lower triangular */
/*        matrix, and U is upper triangular. */

/*     3. If some U(i,i)=0, so that U is exactly singular, then the */
/*        routine returns with INFO = i. Otherwise, the factored form */
/*        of A is used to estimate the condition number of the matrix A. */
/*        If the reciprocal of the condition number is less than machine */
/*        precision, INFO = N+1 is returned as a warning, but the routine */
/*        still goes on to solve for X and compute error bounds as */
/*        described below. */

/*     4. The system of equations is solved for X using the factored form */
/*        of A. */

/*     5. Iterative refinement is applied to improve the computed */
/*        solution matrix and calculate error bounds and backward error */
/*        estimates for it. */

/*     6. If equilibration was used, the matrix X is premultiplied by */
/*        diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so */
/*        that it solves the original system before equilibration. */

/*     REFERENCES */

/*     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J., */
/*         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A., */
/*         Ostrouchov, S., Sorensen, D. */
/*         LAPACK Users' Guide: Second Edition, SIAM, Philadelphia, 1995. */

/*     FURTHER COMMENTS */

/*     This is a simplified version of the LAPACK Library routine DGESVX, */
/*     useful when several sets of matrix equations with the same */
/*     coefficient matrix  A and/or A'  should be solved. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 1999. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Condition number, matrix algebra, matrix operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
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
/*     .. Save Statement .. */
/*     .. */
/*     .. Executable Statements .. */

#line 344 "MB02PD.f"
    /* Parameter adjustments */
#line 344 "MB02PD.f"
    a_dim1 = *lda;
#line 344 "MB02PD.f"
    a_offset = 1 + a_dim1;
#line 344 "MB02PD.f"
    a -= a_offset;
#line 344 "MB02PD.f"
    af_dim1 = *ldaf;
#line 344 "MB02PD.f"
    af_offset = 1 + af_dim1;
#line 344 "MB02PD.f"
    af -= af_offset;
#line 344 "MB02PD.f"
    --ipiv;
#line 344 "MB02PD.f"
    --r__;
#line 344 "MB02PD.f"
    --c__;
#line 344 "MB02PD.f"
    b_dim1 = *ldb;
#line 344 "MB02PD.f"
    b_offset = 1 + b_dim1;
#line 344 "MB02PD.f"
    b -= b_offset;
#line 344 "MB02PD.f"
    x_dim1 = *ldx;
#line 344 "MB02PD.f"
    x_offset = 1 + x_dim1;
#line 344 "MB02PD.f"
    x -= x_offset;
#line 344 "MB02PD.f"
    --ferr;
#line 344 "MB02PD.f"
    --berr;
#line 344 "MB02PD.f"
    --iwork;
#line 344 "MB02PD.f"
    --dwork;
#line 344 "MB02PD.f"

#line 344 "MB02PD.f"
    /* Function Body */
#line 344 "MB02PD.f"
    *info = 0;
#line 345 "MB02PD.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 346 "MB02PD.f"
    equil = lsame_(fact, "E", (ftnlen)1, (ftnlen)1);
#line 347 "MB02PD.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 348 "MB02PD.f"
    if (nofact || equil) {
#line 349 "MB02PD.f"
	*(unsigned char *)equed = 'N';
#line 350 "MB02PD.f"
	rowequ = FALSE_;
#line 351 "MB02PD.f"
	colequ = FALSE_;
#line 352 "MB02PD.f"
    } else {
#line 353 "MB02PD.f"
	rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
		"B", (ftnlen)1, (ftnlen)1);
#line 354 "MB02PD.f"
	colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
		"B", (ftnlen)1, (ftnlen)1);
#line 355 "MB02PD.f"
	smlnum = dlamch_("Safe minimum", (ftnlen)12);
#line 356 "MB02PD.f"
	bignum = 1. / smlnum;
#line 357 "MB02PD.f"
    }

/*     Test the input parameters. */

#line 361 "MB02PD.f"
    if (! nofact && ! equil && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 363 "MB02PD.f"
	*info = -1;
#line 364 "MB02PD.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 366 "MB02PD.f"
	*info = -2;
#line 367 "MB02PD.f"
    } else if (*n < 0) {
#line 368 "MB02PD.f"
	*info = -3;
#line 369 "MB02PD.f"
    } else if (*nrhs < 0) {
#line 370 "MB02PD.f"
	*info = -4;
#line 371 "MB02PD.f"
    } else if (*lda < max(1,*n)) {
#line 372 "MB02PD.f"
	*info = -6;
#line 373 "MB02PD.f"
    } else if (*ldaf < max(1,*n)) {
#line 374 "MB02PD.f"
	*info = -8;
#line 375 "MB02PD.f"
    } else if (lsame_(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rowequ || colequ 
	    || lsame_(equed, "N", (ftnlen)1, (ftnlen)1))) {
#line 377 "MB02PD.f"
	*info = -10;
#line 378 "MB02PD.f"
    } else {
#line 379 "MB02PD.f"
	if (rowequ) {
#line 380 "MB02PD.f"
	    rcmin = bignum;
#line 381 "MB02PD.f"
	    rcmax = 0.;
#line 382 "MB02PD.f"
	    i__1 = *n;
#line 382 "MB02PD.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 383 "MB02PD.f"
		d__1 = rcmin, d__2 = r__[j];
#line 383 "MB02PD.f"
		rcmin = min(d__1,d__2);
/* Computing MAX */
#line 384 "MB02PD.f"
		d__1 = rcmax, d__2 = r__[j];
#line 384 "MB02PD.f"
		rcmax = max(d__1,d__2);
#line 385 "MB02PD.f"
/* L10: */
#line 385 "MB02PD.f"
	    }
#line 386 "MB02PD.f"
	    if (rcmin <= 0.) {
#line 387 "MB02PD.f"
		*info = -11;
#line 388 "MB02PD.f"
	    } else if (*n > 0) {
#line 389 "MB02PD.f"
		rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 390 "MB02PD.f"
	    } else {
#line 391 "MB02PD.f"
		rowcnd = 1.;
#line 392 "MB02PD.f"
	    }
#line 393 "MB02PD.f"
	}
#line 394 "MB02PD.f"
	if (colequ && *info == 0) {
#line 395 "MB02PD.f"
	    rcmin = bignum;
#line 396 "MB02PD.f"
	    rcmax = 0.;
#line 397 "MB02PD.f"
	    i__1 = *n;
#line 397 "MB02PD.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 398 "MB02PD.f"
		d__1 = rcmin, d__2 = c__[j];
#line 398 "MB02PD.f"
		rcmin = min(d__1,d__2);
/* Computing MAX */
#line 399 "MB02PD.f"
		d__1 = rcmax, d__2 = c__[j];
#line 399 "MB02PD.f"
		rcmax = max(d__1,d__2);
#line 400 "MB02PD.f"
/* L20: */
#line 400 "MB02PD.f"
	    }
#line 401 "MB02PD.f"
	    if (rcmin <= 0.) {
#line 402 "MB02PD.f"
		*info = -12;
#line 403 "MB02PD.f"
	    } else if (*n > 0) {
#line 404 "MB02PD.f"
		colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 405 "MB02PD.f"
	    } else {
#line 406 "MB02PD.f"
		colcnd = 1.;
#line 407 "MB02PD.f"
	    }
#line 408 "MB02PD.f"
	}
#line 409 "MB02PD.f"
	if (*info == 0) {
#line 410 "MB02PD.f"
	    if (*ldb < max(1,*n)) {
#line 411 "MB02PD.f"
		*info = -14;
#line 412 "MB02PD.f"
	    } else if (*ldx < max(1,*n)) {
#line 413 "MB02PD.f"
		*info = -16;
#line 414 "MB02PD.f"
	    }
#line 415 "MB02PD.f"
	}
#line 416 "MB02PD.f"
    }

#line 418 "MB02PD.f"
    if (*info != 0) {
#line 419 "MB02PD.f"
	i__1 = -(*info);
#line 419 "MB02PD.f"
	xerbla_("MB02PD", &i__1, (ftnlen)6);
#line 420 "MB02PD.f"
	return 0;
#line 421 "MB02PD.f"
    }

#line 423 "MB02PD.f"
    if (equil) {

/*        Compute row and column scalings to equilibrate the matrix A. */

#line 427 "MB02PD.f"
	dgeequ_(n, n, &a[a_offset], lda, &r__[1], &c__[1], &rowcnd, &colcnd, &
		amax, &infequ);
#line 428 "MB02PD.f"
	if (infequ == 0) {

/*           Equilibrate the matrix. */

#line 432 "MB02PD.f"
	    dlaqge_(n, n, &a[a_offset], lda, &r__[1], &c__[1], &rowcnd, &
		    colcnd, &amax, equed, (ftnlen)1);
#line 434 "MB02PD.f"
	    rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
#line 435 "MB02PD.f"
	    colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
#line 436 "MB02PD.f"
	}
#line 437 "MB02PD.f"
    }

/*     Scale the right hand side. */

#line 441 "MB02PD.f"
    if (notran) {
#line 442 "MB02PD.f"
	if (rowequ) {
#line 443 "MB02PD.f"
	    i__1 = *nrhs;
#line 443 "MB02PD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 444 "MB02PD.f"
		i__2 = *n;
#line 444 "MB02PD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 445 "MB02PD.f"
		    b[i__ + j * b_dim1] = r__[i__] * b[i__ + j * b_dim1];
#line 446 "MB02PD.f"
/* L30: */
#line 446 "MB02PD.f"
		}
#line 447 "MB02PD.f"
/* L40: */
#line 447 "MB02PD.f"
	    }
#line 448 "MB02PD.f"
	}
#line 449 "MB02PD.f"
    } else if (colequ) {
#line 450 "MB02PD.f"
	i__1 = *nrhs;
#line 450 "MB02PD.f"
	for (j = 1; j <= i__1; ++j) {
#line 451 "MB02PD.f"
	    i__2 = *n;
#line 451 "MB02PD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 452 "MB02PD.f"
		b[i__ + j * b_dim1] = c__[i__] * b[i__ + j * b_dim1];
#line 453 "MB02PD.f"
/* L50: */
#line 453 "MB02PD.f"
	    }
#line 454 "MB02PD.f"
/* L60: */
#line 454 "MB02PD.f"
	}
#line 455 "MB02PD.f"
    }

#line 457 "MB02PD.f"
    if (nofact || equil) {

/*        Compute the LU factorization of A. */

#line 461 "MB02PD.f"
	dlacpy_("Full", n, n, &a[a_offset], lda, &af[af_offset], ldaf, (
		ftnlen)4);
#line 462 "MB02PD.f"
	dgetrf_(n, n, &af[af_offset], ldaf, &ipiv[1], info);

/*        Return if INFO is non-zero. */

#line 466 "MB02PD.f"
	if (*info != 0) {
#line 467 "MB02PD.f"
	    if (*info > 0) {

/*              Compute the reciprocal pivot growth factor of the */
/*              leading rank-deficient INFO columns of A. */

#line 472 "MB02PD.f"
		rpvgrw = dlantr_("M", "U", "N", info, info, &af[af_offset], 
			ldaf, &dwork[1], (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 474 "MB02PD.f"
		if (rpvgrw == 0.) {
#line 475 "MB02PD.f"
		    rpvgrw = 1.;
#line 476 "MB02PD.f"
		} else {
#line 477 "MB02PD.f"
		    rpvgrw = dlange_("M", n, info, &a[a_offset], lda, &dwork[
			    1], (ftnlen)1) / rpvgrw;
#line 479 "MB02PD.f"
		}
#line 480 "MB02PD.f"
		dwork[1] = rpvgrw;
#line 481 "MB02PD.f"
		*rcond = 0.;
#line 482 "MB02PD.f"
	    }
#line 483 "MB02PD.f"
	    return 0;
#line 484 "MB02PD.f"
	}

/*        Compute the norm of the matrix A and the */
/*        reciprocal pivot growth factor RPVGRW. */

#line 489 "MB02PD.f"
	if (notran) {
#line 490 "MB02PD.f"
	    *(unsigned char *)norm = '1';
#line 491 "MB02PD.f"
	} else {
#line 492 "MB02PD.f"
	    *(unsigned char *)norm = 'I';
#line 493 "MB02PD.f"
	}
#line 494 "MB02PD.f"
	anorm = dlange_(norm, n, n, &a[a_offset], lda, &dwork[1], (ftnlen)1);
#line 495 "MB02PD.f"
	rpvgrw = dlantr_("M", "U", "N", n, n, &af[af_offset], ldaf, &dwork[1],
		 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 496 "MB02PD.f"
	if (rpvgrw == 0.) {
#line 497 "MB02PD.f"
	    rpvgrw = 1.;
#line 498 "MB02PD.f"
	} else {
#line 499 "MB02PD.f"
	    rpvgrw = dlange_("M", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)
		    1) / rpvgrw;
#line 500 "MB02PD.f"
	}

/*        Compute the reciprocal of the condition number of A. */

#line 504 "MB02PD.f"
	dgecon_(norm, n, &af[af_offset], ldaf, &anorm, rcond, &dwork[1], &
		iwork[1], info, (ftnlen)1);

/*        Set INFO = N+1 if the matrix is singular to working precision. */

#line 509 "MB02PD.f"
	if (*rcond < dlamch_("Epsilon", (ftnlen)7)) {
#line 509 "MB02PD.f"
	    *info = *n + 1;
#line 509 "MB02PD.f"
	}
#line 511 "MB02PD.f"
    }

/*     Compute the solution matrix X. */

#line 515 "MB02PD.f"
    dlacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 516 "MB02PD.f"
    dgetrs_(trans, n, nrhs, &af[af_offset], ldaf, &ipiv[1], &x[x_offset], ldx,
	     info, (ftnlen)1);

/*     Use iterative refinement to improve the computed solution and */
/*     compute error bounds and backward error estimates for it. */

#line 521 "MB02PD.f"
    dgerfs_(trans, n, nrhs, &a[a_offset], lda, &af[af_offset], ldaf, &ipiv[1],
	     &b[b_offset], ldb, &x[x_offset], ldx, &ferr[1], &berr[1], &dwork[
	    1], &iwork[1], info, (ftnlen)1);

/*     Transform the solution matrix X to a solution of the original */
/*     system. */

#line 527 "MB02PD.f"
    if (notran) {
#line 528 "MB02PD.f"
	if (colequ) {
#line 529 "MB02PD.f"
	    i__1 = *nrhs;
#line 529 "MB02PD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 530 "MB02PD.f"
		i__2 = *n;
#line 530 "MB02PD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 531 "MB02PD.f"
		    x[i__ + j * x_dim1] = c__[i__] * x[i__ + j * x_dim1];
#line 532 "MB02PD.f"
/* L70: */
#line 532 "MB02PD.f"
		}
#line 533 "MB02PD.f"
/* L80: */
#line 533 "MB02PD.f"
	    }
#line 534 "MB02PD.f"
	    i__1 = *nrhs;
#line 534 "MB02PD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 535 "MB02PD.f"
		ferr[j] /= colcnd;
#line 536 "MB02PD.f"
/* L90: */
#line 536 "MB02PD.f"
	    }
#line 537 "MB02PD.f"
	}
#line 538 "MB02PD.f"
    } else if (rowequ) {
#line 539 "MB02PD.f"
	i__1 = *nrhs;
#line 539 "MB02PD.f"
	for (j = 1; j <= i__1; ++j) {
#line 540 "MB02PD.f"
	    i__2 = *n;
#line 540 "MB02PD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 541 "MB02PD.f"
		x[i__ + j * x_dim1] = r__[i__] * x[i__ + j * x_dim1];
#line 542 "MB02PD.f"
/* L100: */
#line 542 "MB02PD.f"
	    }
#line 543 "MB02PD.f"
/* L110: */
#line 543 "MB02PD.f"
	}
#line 544 "MB02PD.f"
	i__1 = *nrhs;
#line 544 "MB02PD.f"
	for (j = 1; j <= i__1; ++j) {
#line 545 "MB02PD.f"
	    ferr[j] /= rowcnd;
#line 546 "MB02PD.f"
/* L120: */
#line 546 "MB02PD.f"
	}
#line 547 "MB02PD.f"
    }

#line 549 "MB02PD.f"
    dwork[1] = rpvgrw;
#line 550 "MB02PD.f"
    return 0;

/* *** Last line of MB02PD *** */
} /* mb02pd_ */

