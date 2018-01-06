#line 1 "SB03SD.f"
/* SB03SD.f -- translated by f2c (version 20100827).
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

#line 1 "SB03SD.f"
/* Table of constant values */

static doublereal c_b21 = 0.;
static doublereal c_b23 = -1.;
static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b39 = 1.;

/* Subroutine */ int sb03sd_(char *job, char *fact, char *trana, char *uplo, 
	char *lyapun, integer *n, doublereal *scale, doublereal *a, integer *
	lda, doublereal *t, integer *ldt, doublereal *u, integer *ldu, 
	doublereal *c__, integer *ldc, doublereal *x, integer *ldx, 
	doublereal *sepd, doublereal *rcond, doublereal *ferr, integer *iwork,
	 doublereal *dwork, integer *ldwork, integer *info, ftnlen job_len, 
	ftnlen fact_len, ftnlen trana_len, ftnlen uplo_len, ftnlen lyapun_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, t_dim1, t_offset, u_dim1, 
	    u_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, nn, ldw;
    static doublereal eps;
    static logical jobb, jobc, jobe;
    static integer iabs;
    static char sjob[1];
    static integer ixma, sdim, ires;
    static doublereal epsn, temp, tmax;
    static integer iwrk;
    extern /* Subroutine */ int ma02ed_(char *, integer *, doublereal *, 
	    integer *, ftnlen), dgees_(char *, char *, L_fp, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen), dgemm_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen), mb01ud_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal denom;
    extern /* Subroutine */ int mb01ru_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), mb01rx_(char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), mb01ry_(char *, char *, char *
	    , integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static doublereal anorm, cnorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sb03sx_(char *, char *, char *, integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), sb03sy_(char *, 
	    char *, char *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static logical bwork[1];
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical lower;
    static doublereal xnorm;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen), 
	    dlanhs_(char *, integer *, doublereal *, integer *, doublereal *, 
	    ftnlen);
    static logical nofact;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    extern logical select_();
    static logical update;
    static char tranat[1];
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static logical notrna;
    static doublereal xanorm, thnorm;
    static integer wrkopt;


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

/*     To estimate the conditioning and compute an error bound on the */
/*     solution of the real discrete-time Lyapunov matrix equation */

/*         op(A)'*X*op(A) - X = scale*C */

/*     where op(A) = A or A' (A**T) and C is symmetric (C = C**T). The */
/*     matrix A is N-by-N, the right hand side C and the solution X are */
/*     N-by-N symmetric matrices, and scale is a given scale factor. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the computation to be performed, as follows: */
/*             = 'C':  Compute the reciprocal condition number only; */
/*             = 'E':  Compute the error bound only; */
/*             = 'B':  Compute both the reciprocal condition number and */
/*                     the error bound. */

/*     FACT    CHARACTER*1 */
/*             Specifies whether or not the real Schur factorization */
/*             of the matrix A is supplied on entry, as follows: */
/*             = 'F':  On entry, T and U (if LYAPUN = 'O') contain the */
/*                     factors from the real Schur factorization of the */
/*                     matrix A; */
/*             = 'N':  The Schur factorization of A will be computed */
/*                     and the factors will be stored in T and U (if */
/*                     LYAPUN = 'O'). */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op(A) to be used, as follows: */
/*             = 'N':  op(A) = A    (No transpose); */
/*             = 'T':  op(A) = A**T (Transpose); */
/*             = 'C':  op(A) = A**T (Conjugate transpose = Transpose). */

/*     UPLO    CHARACTER*1 */
/*             Specifies which part of the symmetric matrix C is to be */
/*             used, as follows: */
/*             = 'U':  Upper triangular part; */
/*             = 'L':  Lower triangular part. */

/*     LYAPUN  CHARACTER*1 */
/*             Specifies whether or not the original Lyapunov equations */
/*             should be solved in the iterative estimation process, */
/*             as follows: */
/*             = 'O':  Solve the original Lyapunov equations, updating */
/*                     the right-hand sides and solutions with the */
/*                     matrix U, e.g., X <-- U'*X*U; */
/*             = 'R':  Solve reduced Lyapunov equations only, without */
/*                     updating the right-hand sides and solutions. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, X and C.  N >= 0. */

/*     SCALE   (input) DOUBLE PRECISION */
/*             The scale factor, scale, set by a Lyapunov solver. */
/*             0 <= SCALE <= 1. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             If FACT = 'N' or LYAPUN = 'O', the leading N-by-N part of */
/*             this array must contain the original matrix A. */
/*             If FACT = 'F' and LYAPUN = 'R', A is not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. */
/*             LDA >= MAX(1,N), if FACT = 'N' or  LYAPUN = 'O'; */
/*             LDA >= 1,        if FACT = 'F' and LYAPUN = 'R'. */

/*     T       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDT,N) */
/*             If FACT = 'F', then on entry the leading N-by-N upper */
/*             Hessenberg part of this array must contain the upper */
/*             quasi-triangular matrix T in Schur canonical form from a */
/*             Schur factorization of A. */
/*             If FACT = 'N', then this array need not be set on input. */
/*             On exit, (if INFO = 0 or INFO = N+1, for FACT = 'N') the */
/*             leading N-by-N upper Hessenberg part of this array */
/*             contains the upper quasi-triangular matrix T in Schur */
/*             canonical form from a Schur factorization of A. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T.  LDT >= MAX(1,N). */

/*     U       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDU,N) */
/*             If LYAPUN = 'O' and FACT = 'F', then U is an input */
/*             argument and on entry, the leading N-by-N part of this */
/*             array must contain the orthogonal matrix U from a real */
/*             Schur factorization of A. */
/*             If LYAPUN = 'O' and FACT = 'N', then U is an output */
/*             argument and on exit, if INFO = 0 or INFO = N+1, it */
/*             contains the orthogonal N-by-N matrix from a real Schur */
/*             factorization of A. */
/*             If LYAPUN = 'R', the array U is not referenced. */

/*     LDU     INTEGER */
/*             The leading dimension of the array U. */
/*             LDU >= 1,        if LYAPUN = 'R'; */
/*             LDU >= MAX(1,N), if LYAPUN = 'O'. */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             If UPLO = 'U', the leading N-by-N upper triangular part of */
/*             this array must contain the upper triangular part of the */
/*             matrix C of the original Lyapunov equation (with */
/*             matrix A), if LYAPUN = 'O', or of the reduced Lyapunov */
/*             equation (with matrix T), if LYAPUN = 'R'. */
/*             If UPLO = 'L', the leading N-by-N lower triangular part of */
/*             this array must contain the lower triangular part of the */
/*             matrix C of the original Lyapunov equation (with */
/*             matrix A), if LYAPUN = 'O', or of the reduced Lyapunov */
/*             equation (with matrix T), if LYAPUN = 'R'. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= MAX(1,N). */

/*     X       (input) DOUBLE PRECISION array, dimension (LDX,N) */
/*             The leading N-by-N part of this array must contain the */
/*             symmetric solution matrix X of the original Lyapunov */
/*             equation (with matrix A), if LYAPUN = 'O', or of the */
/*             reduced Lyapunov equation (with matrix T), if */
/*             LYAPUN = 'R'. */
/*             The array X is modified internally, but restored on exit. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X.  LDX >= MAX(1,N). */

/*     SEPD    (output) DOUBLE PRECISION */
/*             If JOB = 'C' or JOB = 'B', the estimated quantity */
/*             sepd(op(A),op(A)'). */
/*             If N = 0, or X = 0, or JOB = 'E', SEPD is not referenced. */

/*     RCOND   (output) DOUBLE PRECISION */
/*             If JOB = 'C' or JOB = 'B', an estimate of the reciprocal */
/*             condition number of the discrete-time Lyapunov equation. */
/*             If N = 0 or X = 0, RCOND is set to 1 or 0, respectively. */
/*             If JOB = 'E', RCOND is not referenced. */

/*     FERR    (output) DOUBLE PRECISION */
/*             If JOB = 'E' or JOB = 'B', an estimated forward error */
/*             bound for the solution X. If XTRUE is the true solution, */
/*             FERR bounds the magnitude of the largest entry in */
/*             (X - XTRUE) divided by the magnitude of the largest entry */
/*             in X. */
/*             If N = 0 or X = 0, FERR is set to 0. */
/*             If JOB = 'C', FERR is not referenced. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N*N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0 or INFO = N+1, DWORK(1) returns the */
/*             optimal value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= 1,                            if N = 0; else, */
/*             LDWORK >= MAX(3,2*N*N) + N*N,           if JOB  = 'C', */
/*                                                        FACT = 'F'; */
/*             LDWORK >= MAX(MAX(3,2*N*N) + N*N, 5*N), if JOB  = 'C', */
/*                                                        FACT = 'N'; */
/*             LDWORK >= MAX(3,2*N*N) + N*N + 2*N,     if JOB  = 'E', or */
/*                                                        JOB  = 'B'. */
/*             For optimum performance LDWORK should sometimes be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, i <= N, the QR algorithm failed to */
/*                   complete the reduction to Schur canonical form (see */
/*                   LAPACK Library routine DGEES); on exit, the matrix */
/*                   T(i+1:N,i+1:N) contains the partially converged */
/*                   Schur form, and DWORK(i+1:N) and DWORK(N+i+1:2*N) */
/*                   contain the real and imaginary parts, respectively, */
/*                   of the converged eigenvalues; this error is unlikely */
/*                   to appear; */
/*             = N+1:  if the matrix T has almost reciprocal eigenvalues; */
/*                   perturbed values were used to solve Lyapunov */
/*                   equations, but the matrix T, if given (for */
/*                   FACT = 'F'), is unchanged. */

/*     METHOD */

/*     The condition number of the discrete-time Lyapunov equation is */
/*     estimated as */

/*     cond = (norm(Theta)*norm(A) + norm(inv(Omega))*norm(C))/norm(X), */

/*     where Omega and Theta are linear operators defined by */

/*     Omega(W) = op(A)'*W*op(A) - W, */
/*     Theta(W) = inv(Omega(op(W)'*X*op(A) + op(A)'*X*op(W))). */

/*     The routine estimates the quantities */

/*     sepd(op(A),op(A)') = 1 / norm(inv(Omega)) */

/*     and norm(Theta) using 1-norm condition estimators. */

/*     The forward error bound is estimated using a practical error bound */
/*     similar to the one proposed in [1]. */

/*     REFERENCES */

/*     [1] Higham, N.J. */
/*         Perturbation theory and backward error for AX-XB=C. */
/*         BIT, vol. 33, pp. 124-136, 1993. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */
/*     The accuracy of the estimates obtained depends on the solution */
/*     accuracy and on the properties of the 1-norm estimator. */

/*     FURTHER COMMENTS */

/*     The option LYAPUN = 'R' may occasionally produce slightly worse */
/*     or better estimates, and it is much faster than the option 'O'. */
/*     When SEPD is computed and it is zero, the routine returns */
/*     immediately, with RCOND and FERR (if requested) set to 0 and 1, */
/*     respectively. In this case, the equation is singular. */

/*     CONTRIBUTORS */

/*     P. Petkov, Tech. University of Sofia, December 1998. */
/*     V. Sima, Katholieke Univ. Leuven, Belgium, February 1999. */

/*     REVISIONS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, March 2003. */

/*     KEYWORDS */

/*     Lyapunov equation, orthogonal transformation, real Schur form. */

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

#line 311 "SB03SD.f"
    /* Parameter adjustments */
#line 311 "SB03SD.f"
    a_dim1 = *lda;
#line 311 "SB03SD.f"
    a_offset = 1 + a_dim1;
#line 311 "SB03SD.f"
    a -= a_offset;
#line 311 "SB03SD.f"
    t_dim1 = *ldt;
#line 311 "SB03SD.f"
    t_offset = 1 + t_dim1;
#line 311 "SB03SD.f"
    t -= t_offset;
#line 311 "SB03SD.f"
    u_dim1 = *ldu;
#line 311 "SB03SD.f"
    u_offset = 1 + u_dim1;
#line 311 "SB03SD.f"
    u -= u_offset;
#line 311 "SB03SD.f"
    c_dim1 = *ldc;
#line 311 "SB03SD.f"
    c_offset = 1 + c_dim1;
#line 311 "SB03SD.f"
    c__ -= c_offset;
#line 311 "SB03SD.f"
    x_dim1 = *ldx;
#line 311 "SB03SD.f"
    x_offset = 1 + x_dim1;
#line 311 "SB03SD.f"
    x -= x_offset;
#line 311 "SB03SD.f"
    --iwork;
#line 311 "SB03SD.f"
    --dwork;
#line 311 "SB03SD.f"

#line 311 "SB03SD.f"
    /* Function Body */
#line 311 "SB03SD.f"
    jobc = lsame_(job, "C", (ftnlen)1, (ftnlen)1);
#line 312 "SB03SD.f"
    jobe = lsame_(job, "E", (ftnlen)1, (ftnlen)1);
#line 313 "SB03SD.f"
    jobb = lsame_(job, "B", (ftnlen)1, (ftnlen)1);
#line 314 "SB03SD.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 315 "SB03SD.f"
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
#line 316 "SB03SD.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 317 "SB03SD.f"
    update = lsame_(lyapun, "O", (ftnlen)1, (ftnlen)1);

#line 319 "SB03SD.f"
    nn = *n * *n;
/* Computing MAX */
#line 320 "SB03SD.f"
    i__1 = 3, i__2 = nn << 1;
#line 320 "SB03SD.f"
    ldw = max(i__1,i__2) + nn;

#line 322 "SB03SD.f"
    *info = 0;
#line 323 "SB03SD.f"
    if (! (jobb || jobc || jobe)) {
#line 324 "SB03SD.f"
	*info = -1;
#line 325 "SB03SD.f"
    } else if (! (nofact || lsame_(fact, "F", (ftnlen)1, (ftnlen)1))) {
#line 326 "SB03SD.f"
	*info = -2;
#line 327 "SB03SD.f"
    } else if (! (notrna || lsame_(trana, "T", (ftnlen)1, (ftnlen)1) || 
	    lsame_(trana, "C", (ftnlen)1, (ftnlen)1))) {
#line 329 "SB03SD.f"
	*info = -3;
#line 330 "SB03SD.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 331 "SB03SD.f"
	*info = -4;
#line 332 "SB03SD.f"
    } else if (! (update || lsame_(lyapun, "R", (ftnlen)1, (ftnlen)1))) {
#line 333 "SB03SD.f"
	*info = -5;
#line 334 "SB03SD.f"
    } else if (*n < 0) {
#line 335 "SB03SD.f"
	*info = -6;
#line 336 "SB03SD.f"
    } else if (*scale < 0. || *scale > 1.) {
#line 337 "SB03SD.f"
	*info = -7;
#line 338 "SB03SD.f"
    } else if (*lda < 1 || *lda < *n && (update || nofact)) {
#line 340 "SB03SD.f"
	*info = -9;
#line 341 "SB03SD.f"
    } else if (*ldt < max(1,*n)) {
#line 342 "SB03SD.f"
	*info = -11;
#line 343 "SB03SD.f"
    } else if (*ldu < 1 || *ldu < *n && update) {
#line 344 "SB03SD.f"
	*info = -13;
#line 345 "SB03SD.f"
    } else if (*ldc < max(1,*n)) {
#line 346 "SB03SD.f"
	*info = -15;
#line 347 "SB03SD.f"
    } else if (*ldx < max(1,*n)) {
#line 348 "SB03SD.f"
	*info = -17;
#line 349 "SB03SD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 349 "SB03SD.f"
	i__1 = ldw, i__2 = *n * 5;
#line 349 "SB03SD.f"
	if (*ldwork < 1 || *ldwork < ldw && jobc && ! nofact || *ldwork < max(
		i__1,i__2) && jobc && nofact || *ldwork < ldw + (*n << 1) && !
		 jobc) {
#line 353 "SB03SD.f"
	    *info = -23;
#line 354 "SB03SD.f"
	}
#line 354 "SB03SD.f"
    }

#line 356 "SB03SD.f"
    if (*info != 0) {
#line 357 "SB03SD.f"
	i__1 = -(*info);
#line 357 "SB03SD.f"
	xerbla_("SB03SD", &i__1, (ftnlen)6);
#line 358 "SB03SD.f"
	return 0;
#line 359 "SB03SD.f"
    }

/*     Quick return if possible. */

#line 363 "SB03SD.f"
    if (*n == 0) {
#line 364 "SB03SD.f"
	if (! jobe) {
#line 364 "SB03SD.f"
	    *rcond = 1.;
#line 364 "SB03SD.f"
	}
#line 366 "SB03SD.f"
	if (! jobc) {
#line 366 "SB03SD.f"
	    *ferr = 0.;
#line 366 "SB03SD.f"
	}
#line 368 "SB03SD.f"
	dwork[1] = 1.;
#line 369 "SB03SD.f"
	return 0;
#line 370 "SB03SD.f"
    }

/*     Compute the 1-norm of the matrix X. */

#line 374 "SB03SD.f"
    xnorm = dlansy_("1-norm", uplo, n, &x[x_offset], ldx, &dwork[1], (ftnlen)
	    6, (ftnlen)1);
#line 375 "SB03SD.f"
    if (xnorm == 0.) {

/*        The solution is zero. */

#line 379 "SB03SD.f"
	if (! jobe) {
#line 379 "SB03SD.f"
	    *rcond = 0.;
#line 379 "SB03SD.f"
	}
#line 381 "SB03SD.f"
	if (! jobc) {
#line 381 "SB03SD.f"
	    *ferr = 0.;
#line 381 "SB03SD.f"
	}
#line 383 "SB03SD.f"
	dwork[1] = (doublereal) (*n);
#line 384 "SB03SD.f"
	return 0;
#line 385 "SB03SD.f"
    }

/*     Compute the 1-norm of A or T. */

#line 389 "SB03SD.f"
    if (nofact || update) {
#line 390 "SB03SD.f"
	anorm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)
		6);
#line 391 "SB03SD.f"
    } else {
#line 392 "SB03SD.f"
	anorm = dlanhs_("1-norm", n, &t[t_offset], ldt, &dwork[1], (ftnlen)6);
#line 393 "SB03SD.f"
    }

/*     For the special case A = I, set SEPD and RCOND to 0. */
/*     For the special case A = 0, set SEPD and RCOND to 1. */
/*     A quick test is used in general. */

#line 399 "SB03SD.f"
    if (anorm == 1.) {
#line 400 "SB03SD.f"
	if (nofact || update) {
#line 401 "SB03SD.f"
	    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[1], n, (ftnlen)4);
#line 402 "SB03SD.f"
	} else {
#line 403 "SB03SD.f"
	    dlacpy_("Full", n, n, &t[t_offset], ldt, &dwork[1], n, (ftnlen)4);
#line 404 "SB03SD.f"
	    if (*n > 2) {
#line 404 "SB03SD.f"
		i__1 = *n - 2;
#line 404 "SB03SD.f"
		i__2 = *n - 2;
#line 404 "SB03SD.f"
		dlaset_("Lower", &i__1, &i__2, &c_b21, &c_b21, &dwork[3], n, (
			ftnlen)5);
#line 404 "SB03SD.f"
	    }
#line 407 "SB03SD.f"
	}
#line 408 "SB03SD.f"
	dwork[nn + 1] = 1.;
#line 409 "SB03SD.f"
	i__1 = *n + 1;
#line 409 "SB03SD.f"
	daxpy_(n, &c_b23, &dwork[nn + 1], &c__0, &dwork[1], &i__1);
#line 410 "SB03SD.f"
	if (dlange_("Max", n, n, &dwork[1], n, &dwork[1], (ftnlen)3) == 0.) {
#line 411 "SB03SD.f"
	    if (! jobe) {
#line 412 "SB03SD.f"
		*sepd = 0.;
#line 413 "SB03SD.f"
		*rcond = 0.;
#line 414 "SB03SD.f"
	    }
#line 415 "SB03SD.f"
	    if (! jobc) {
#line 415 "SB03SD.f"
		*ferr = 1.;
#line 415 "SB03SD.f"
	    }
#line 417 "SB03SD.f"
	    dwork[1] = (doublereal) (nn + 1);
#line 418 "SB03SD.f"
	    return 0;
#line 419 "SB03SD.f"
	}

#line 421 "SB03SD.f"
    } else if (anorm == 0.) {
#line 422 "SB03SD.f"
	if (! jobe) {
#line 423 "SB03SD.f"
	    *sepd = 1.;
#line 424 "SB03SD.f"
	    *rcond = 1.;
#line 425 "SB03SD.f"
	}
#line 426 "SB03SD.f"
	if (jobc) {
#line 427 "SB03SD.f"
	    dwork[1] = (doublereal) (*n);
#line 428 "SB03SD.f"
	    return 0;
#line 429 "SB03SD.f"
	} else {

/*           Set FERR for the special case A = 0. */

#line 433 "SB03SD.f"
	    dlacpy_(uplo, n, n, &x[x_offset], ldx, &dwork[1], n, (ftnlen)1);

#line 435 "SB03SD.f"
	    if (lower) {
#line 436 "SB03SD.f"
		i__1 = *n;
#line 436 "SB03SD.f"
		for (j = 1; j <= i__1; ++j) {
#line 437 "SB03SD.f"
		    i__2 = *n - j + 1;
#line 437 "SB03SD.f"
		    daxpy_(&i__2, scale, &c__[j + j * c_dim1], &c__1, &dwork[(
			    j - 1) * *n + j], &c__1);
#line 439 "SB03SD.f"
/* L10: */
#line 439 "SB03SD.f"
		}
#line 440 "SB03SD.f"
	    } else {
#line 441 "SB03SD.f"
		i__1 = *n;
#line 441 "SB03SD.f"
		for (j = 1; j <= i__1; ++j) {
#line 442 "SB03SD.f"
		    daxpy_(&j, scale, &c__[j * c_dim1 + 1], &c__1, &dwork[(j 
			    - 1) * *n + 1], &c__1);
#line 444 "SB03SD.f"
/* L20: */
#line 444 "SB03SD.f"
		}
#line 445 "SB03SD.f"
	    }

/* Computing MIN */
#line 447 "SB03SD.f"
	    d__1 = 1., d__2 = dlansy_("1-norm", uplo, n, &dwork[1], n, &dwork[
		    nn + 1], (ftnlen)6, (ftnlen)1) / xnorm;
#line 447 "SB03SD.f"
	    *ferr = min(d__1,d__2);
#line 449 "SB03SD.f"
	    dwork[1] = (doublereal) (nn + *n);
#line 450 "SB03SD.f"
	    return 0;
#line 451 "SB03SD.f"
	}
#line 452 "SB03SD.f"
    }

/*     General case. */

#line 456 "SB03SD.f"
    cnorm = dlansy_("1-norm", uplo, n, &c__[c_offset], ldc, &dwork[1], (
	    ftnlen)6, (ftnlen)1);

/*     Workspace usage. */

#line 460 "SB03SD.f"
    iabs = nn;
/* Computing MAX */
#line 461 "SB03SD.f"
    i__1 = 3, i__2 = nn << 1;
#line 461 "SB03SD.f"
    ixma = max(i__1,i__2);
#line 462 "SB03SD.f"
    ires = ixma;
#line 463 "SB03SD.f"
    iwrk = ixma + nn;
#line 464 "SB03SD.f"
    wrkopt = 0;

#line 466 "SB03SD.f"
    if (nofact) {

/*        Compute the Schur factorization of A, A = U*T*U'. */
/*        Workspace:  need   5*N; */
/*                    prefer larger. */
/*        (Note: Comments in the code beginning "Workspace:" describe the */
/*        minimal amount of real workspace needed at that point in the */
/*        code, as well as the preferred amount for good performance.) */

#line 475 "SB03SD.f"
	dlacpy_("Full", n, n, &a[a_offset], lda, &t[t_offset], ldt, (ftnlen)4)
		;
#line 476 "SB03SD.f"
	if (update) {
#line 477 "SB03SD.f"
	    *(unsigned char *)sjob = 'V';
#line 478 "SB03SD.f"
	} else {
#line 479 "SB03SD.f"
	    *(unsigned char *)sjob = 'N';
#line 480 "SB03SD.f"
	}
#line 481 "SB03SD.f"
	i__1 = *ldwork - (*n << 1);
#line 481 "SB03SD.f"
	dgees_(sjob, "Not ordered", (L_fp)select_, n, &t[t_offset], ldt, &
		sdim, &dwork[1], &dwork[*n + 1], &u[u_offset], ldu, &dwork[(*
		n << 1) + 1], &i__1, bwork, info, (ftnlen)1, (ftnlen)11);
#line 484 "SB03SD.f"
	if (*info > 0) {
#line 484 "SB03SD.f"
	    return 0;
#line 484 "SB03SD.f"
	}
#line 486 "SB03SD.f"
	wrkopt = (integer) dwork[(*n << 1) + 1] + (*n << 1);
#line 487 "SB03SD.f"
    }

/*     Compute X*op(A) or X*op(T). */

#line 491 "SB03SD.f"
    if (update) {
#line 492 "SB03SD.f"
	dgemm_("NoTranspose", trana, n, n, n, &c_b39, &x[x_offset], ldx, &a[
		a_offset], lda, &c_b21, &dwork[ixma + 1], n, (ftnlen)11, (
		ftnlen)1);
#line 494 "SB03SD.f"
    } else {
#line 495 "SB03SD.f"
	mb01ud_("Right", trana, n, n, &c_b39, &t[t_offset], ldt, &x[x_offset],
		 ldx, &dwork[ixma + 1], n, info, (ftnlen)5, (ftnlen)1);
#line 497 "SB03SD.f"
    }

#line 499 "SB03SD.f"
    if (! jobe) {

/*        Estimate sepd(op(A),op(A)') = sepd(op(T),op(T)') and */
/*        norm(Theta). */
/*        Workspace max(3,2*N*N) + N*N. */

#line 505 "SB03SD.f"
	sb03sy_("Both", trana, lyapun, n, &t[t_offset], ldt, &u[u_offset], 
		ldu, &dwork[ixma + 1], n, sepd, &thnorm, &iwork[1], &dwork[1],
		 &ixma, info, (ftnlen)4, (ftnlen)1, (ftnlen)1);

/* Computing MAX */
/* Computing MAX */
#line 509 "SB03SD.f"
	i__3 = 3, i__4 = nn << 1;
#line 509 "SB03SD.f"
	i__1 = wrkopt, i__2 = max(i__3,i__4) + nn;
#line 509 "SB03SD.f"
	wrkopt = max(i__1,i__2);

/*        Return if the equation is singular. */

#line 513 "SB03SD.f"
	if (*sepd == 0.) {
#line 514 "SB03SD.f"
	    *rcond = 0.;
#line 515 "SB03SD.f"
	    if (jobb) {
#line 515 "SB03SD.f"
		*ferr = 1.;
#line 515 "SB03SD.f"
	    }
#line 517 "SB03SD.f"
	    dwork[1] = (doublereal) wrkopt;
#line 518 "SB03SD.f"
	    return 0;
#line 519 "SB03SD.f"
	}

/*        Estimate the reciprocal condition number. */

/* Computing MAX */
#line 523 "SB03SD.f"
	d__1 = max(*sepd,xnorm);
#line 523 "SB03SD.f"
	tmax = max(d__1,anorm);
#line 524 "SB03SD.f"
	if (tmax <= 1.) {
#line 525 "SB03SD.f"
	    temp = *sepd * xnorm;
#line 526 "SB03SD.f"
	    denom = *scale * cnorm + *sepd * anorm * thnorm;
#line 527 "SB03SD.f"
	} else {
#line 528 "SB03SD.f"
	    temp = *sepd / tmax * (xnorm / tmax);
#line 529 "SB03SD.f"
	    denom = *scale / tmax * (cnorm / tmax) + *sepd / tmax * (anorm / 
		    tmax) * thnorm;
#line 531 "SB03SD.f"
	}
#line 532 "SB03SD.f"
	if (temp >= denom) {
#line 533 "SB03SD.f"
	    *rcond = 1.;
#line 534 "SB03SD.f"
	} else {
#line 535 "SB03SD.f"
	    *rcond = temp / denom;
#line 536 "SB03SD.f"
	}
#line 537 "SB03SD.f"
    }

#line 539 "SB03SD.f"
    if (! jobc) {

/*        Form a triangle of the residual matrix */
/*        R = scale*C + X - op(A)'*X*op(A), or */
/*        R = scale*C + X - op(T)'*X*op(T), */
/*        exploiting the symmetry. For memory savings, R is formed in the */
/*        leading N-by-N upper/lower triangular part of DWORK, and it is */
/*        finally moved in the location where X*op(A) or X*op(T) was */
/*        stored, freeing workspace for the SB03SX call. */

#line 549 "SB03SD.f"
	if (notrna) {
#line 550 "SB03SD.f"
	    *(unsigned char *)tranat = 'T';
#line 551 "SB03SD.f"
	} else {
#line 552 "SB03SD.f"
	    *(unsigned char *)tranat = 'N';
#line 553 "SB03SD.f"
	}

#line 555 "SB03SD.f"
	dlacpy_(uplo, n, n, &c__[c_offset], ldc, &dwork[1], n, (ftnlen)1);

#line 557 "SB03SD.f"
	if (update) {
#line 558 "SB03SD.f"
	    mb01rx_("Left", uplo, tranat, n, n, scale, &c_b23, &dwork[1], n, &
		    a[a_offset], lda, &dwork[ixma + 1], n, info, (ftnlen)4, (
		    ftnlen)1, (ftnlen)1);
#line 560 "SB03SD.f"
	} else {
#line 561 "SB03SD.f"
	    mb01ry_("Left", uplo, tranat, n, scale, &c_b23, &dwork[1], n, &t[
		    t_offset], ldt, &dwork[ixma + 1], n, &dwork[iwrk + 1], 
		    info, (ftnlen)4, (ftnlen)1, (ftnlen)1);
#line 564 "SB03SD.f"
	}

#line 566 "SB03SD.f"
	if (lower) {
#line 567 "SB03SD.f"
	    i__1 = *n;
#line 567 "SB03SD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 568 "SB03SD.f"
		i__2 = *n - j + 1;
#line 568 "SB03SD.f"
		daxpy_(&i__2, &c_b39, &x[j + j * x_dim1], &c__1, &dwork[(j - 
			1) * *n + j], &c__1);
#line 570 "SB03SD.f"
/* L30: */
#line 570 "SB03SD.f"
	    }
#line 571 "SB03SD.f"
	} else {
#line 572 "SB03SD.f"
	    i__1 = *n;
#line 572 "SB03SD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 573 "SB03SD.f"
		daxpy_(&j, &c_b39, &x[j * x_dim1 + 1], &c__1, &dwork[(j - 1) *
			 *n + 1], &c__1);
#line 574 "SB03SD.f"
/* L40: */
#line 574 "SB03SD.f"
	    }
#line 575 "SB03SD.f"
	}

#line 577 "SB03SD.f"
	dlacpy_(uplo, n, n, &dwork[1], n, &dwork[ires + 1], n, (ftnlen)1);

/*        Get the machine precision. */

#line 581 "SB03SD.f"
	eps = dlamch_("Epsilon", (ftnlen)7);
#line 582 "SB03SD.f"
	epsn = eps * (doublereal) ((*n << 1) + 2);

/*        Add to abs(R) a term that takes account of rounding errors in */
/*        forming R: */
/*          abs(R) := abs(R) + EPS*(3*scale*abs(C) + 3*abs(X) + */
/*                    2*(n+1)*abs(op(A))'*abs(X)*abs(op(A))), or */
/*          abs(R) := abs(R) + EPS*(3*scale*abs(C) + 3*abs(X) + */
/*                    2*(n+1)*abs(op(T))'*abs(X)*abs(op(T))), */
/*        where EPS is the machine precision. */
/*        Workspace max(3,2*N*N) + N*N + 2*N. */
/*        Note that the lower or upper triangular part of X specified by */
/*        UPLO is used as workspace, but it is finally restored. */

#line 595 "SB03SD.f"
	if (update) {
#line 596 "SB03SD.f"
	    i__1 = *n;
#line 596 "SB03SD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 597 "SB03SD.f"
		i__2 = *n;
#line 597 "SB03SD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 598 "SB03SD.f"
		    dwork[iabs + (j - 1) * *n + i__] = (d__1 = a[i__ + j * 
			    a_dim1], abs(d__1));
#line 599 "SB03SD.f"
/* L50: */
#line 599 "SB03SD.f"
		}
#line 600 "SB03SD.f"
/* L60: */
#line 600 "SB03SD.f"
	    }
#line 601 "SB03SD.f"
	} else {
#line 602 "SB03SD.f"
	    i__1 = *n;
#line 602 "SB03SD.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 603 "SB03SD.f"
		i__3 = j + 1;
#line 603 "SB03SD.f"
		i__2 = min(i__3,*n);
#line 603 "SB03SD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 604 "SB03SD.f"
		    dwork[iabs + (j - 1) * *n + i__] = (d__1 = t[i__ + j * 
			    t_dim1], abs(d__1));
#line 605 "SB03SD.f"
/* L70: */
#line 605 "SB03SD.f"
		}
#line 606 "SB03SD.f"
/* L80: */
#line 606 "SB03SD.f"
	    }
#line 607 "SB03SD.f"
	}

#line 609 "SB03SD.f"
	i__1 = *ldx + 1;
#line 609 "SB03SD.f"
	dcopy_(n, &x[x_offset], &i__1, &dwork[iwrk + 1], &c__1);

#line 611 "SB03SD.f"
	if (lower) {
#line 612 "SB03SD.f"
	    i__1 = *n;
#line 612 "SB03SD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 613 "SB03SD.f"
		i__2 = *n;
#line 613 "SB03SD.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 614 "SB03SD.f"
		    temp = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 615 "SB03SD.f"
		    x[i__ + j * x_dim1] = temp;
#line 616 "SB03SD.f"
		    dwork[ires + (j - 1) * *n + i__] = (d__2 = dwork[ires + (
			    j - 1) * *n + i__], abs(d__2)) + eps * 3. * (*
			    scale * (d__1 = c__[i__ + j * c_dim1], abs(d__1)) 
			    + temp);
#line 619 "SB03SD.f"
/* L90: */
#line 619 "SB03SD.f"
		}
#line 620 "SB03SD.f"
/* L100: */
#line 620 "SB03SD.f"
	    }
#line 621 "SB03SD.f"
	} else {
#line 622 "SB03SD.f"
	    i__1 = *n;
#line 622 "SB03SD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 623 "SB03SD.f"
		i__2 = j;
#line 623 "SB03SD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 624 "SB03SD.f"
		    temp = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 625 "SB03SD.f"
		    x[i__ + j * x_dim1] = temp;
#line 626 "SB03SD.f"
		    dwork[ires + (j - 1) * *n + i__] = (d__2 = dwork[ires + (
			    j - 1) * *n + i__], abs(d__2)) + eps * 3. * (*
			    scale * (d__1 = c__[i__ + j * c_dim1], abs(d__1)) 
			    + temp);
#line 629 "SB03SD.f"
/* L110: */
#line 629 "SB03SD.f"
		}
#line 630 "SB03SD.f"
/* L120: */
#line 630 "SB03SD.f"
	    }
#line 631 "SB03SD.f"
	}

#line 633 "SB03SD.f"
	if (update) {
#line 634 "SB03SD.f"
	    mb01ru_(uplo, tranat, n, n, &c_b39, &epsn, &dwork[ires + 1], n, &
		    dwork[iabs + 1], n, &x[x_offset], ldx, &dwork[1], &nn, 
		    info, (ftnlen)1, (ftnlen)1);
#line 637 "SB03SD.f"
	} else {

/*           Compute W = abs(X)*abs(op(T)), and then premultiply by */
/*           abs(T)' and add in the result. */

#line 642 "SB03SD.f"
	    mb01ud_("Right", trana, n, n, &c_b39, &dwork[iabs + 1], n, &x[
		    x_offset], ldx, &dwork[1], n, info, (ftnlen)5, (ftnlen)1);
#line 644 "SB03SD.f"
	    mb01ry_("Left", uplo, tranat, n, &c_b39, &epsn, &dwork[ires + 1], 
		    n, &dwork[iabs + 1], n, &dwork[1], n, &dwork[iwrk + *n + 
		    1], info, (ftnlen)4, (ftnlen)1, (ftnlen)1);
#line 647 "SB03SD.f"
	}

/* Computing MAX */
/* Computing MAX */
#line 649 "SB03SD.f"
	i__3 = 3, i__4 = nn << 1;
#line 649 "SB03SD.f"
	i__1 = wrkopt, i__2 = max(i__3,i__4) + nn + (*n << 1);
#line 649 "SB03SD.f"
	wrkopt = max(i__1,i__2);

/*        Restore X. */

#line 653 "SB03SD.f"
	i__1 = *ldx + 1;
#line 653 "SB03SD.f"
	dcopy_(n, &dwork[iwrk + 1], &c__1, &x[x_offset], &i__1);
#line 654 "SB03SD.f"
	if (lower) {
#line 655 "SB03SD.f"
	    ma02ed_("Upper", n, &x[x_offset], ldx, (ftnlen)5);
#line 656 "SB03SD.f"
	} else {
#line 657 "SB03SD.f"
	    ma02ed_("Lower", n, &x[x_offset], ldx, (ftnlen)5);
#line 658 "SB03SD.f"
	}

/*        Compute forward error bound, using matrix norm estimator. */
/*        Workspace max(3,2*N*N) + N*N. */

#line 663 "SB03SD.f"
	xanorm = dlansy_("Max", uplo, n, &x[x_offset], ldx, &dwork[1], (
		ftnlen)3, (ftnlen)1);

#line 665 "SB03SD.f"
	sb03sx_(trana, uplo, lyapun, n, &xanorm, &t[t_offset], ldt, &u[
		u_offset], ldu, &dwork[ires + 1], n, ferr, &iwork[1], &dwork[
		1], &ires, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 668 "SB03SD.f"
    }

#line 670 "SB03SD.f"
    dwork[1] = (doublereal) wrkopt;
#line 671 "SB03SD.f"
    return 0;

/* *** Last line of SB03SD *** */
} /* sb03sd_ */

