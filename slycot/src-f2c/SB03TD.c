#line 1 "SB03TD.f"
/* SB03TD.f -- translated by f2c (version 20100827).
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

#line 1 "SB03TD.f"
/* Table of constant values */

static doublereal c_b23 = 0.;
static doublereal c_b24 = 1.;
static doublereal c_b25 = .5;

/* Subroutine */ int sb03td_(char *job, char *fact, char *trana, char *uplo, 
	char *lyapun, integer *n, doublereal *scale, doublereal *a, integer *
	lda, doublereal *t, integer *ldt, doublereal *u, integer *ldu, 
	doublereal *c__, integer *ldc, doublereal *x, integer *ldx, 
	doublereal *sep, doublereal *rcond, doublereal *ferr, doublereal *wr, 
	doublereal *wi, integer *iwork, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen job_len, ftnlen fact_len, ftnlen trana_len, 
	ftnlen uplo_len, ftnlen lyapun_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, t_dim1, t_offset, u_dim1, 
	    u_offset, x_dim1, x_offset, i__1, i__2;

    /* Local variables */
    static integer nn, ldw;
    static logical joba, jobc, jobe;
    static char jobl[1];
    static integer sdim;
    static logical jobs;
    static char sjob[1];
    static logical jobx;
    extern /* Subroutine */ int ma02ed_(char *, integer *, doublereal *, 
	    integer *, ftnlen);
    static char cfact[1];
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgees_(char *, char *, L_fp, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, logical *, integer *, ftnlen, 
	    ftnlen), sb03qd_(char *, char *, char *, char *, char *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb01ru_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), sb03my_(char *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), sb03qy_(char *, char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static logical bwork[1], lower, nofact;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern logical select_();
    static logical update, notrna;
    static doublereal thnorm;


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

/*     To solve the real continuous-time Lyapunov matrix equation */

/*            op(A)'*X + X*op(A) = scale*C, */

/*     estimate the conditioning, and compute an error bound on the */
/*     solution X, where op(A) = A or A' (A**T), the matrix A is N-by-N, */
/*     the right hand side C and the solution X are N-by-N symmetric */
/*     matrices (C = C', X = X'), and scale is an output scale factor, */
/*     set less than or equal to 1 to avoid overflow in X. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the computation to be performed, as follows: */
/*             = 'X':  Compute the solution only; */
/*             = 'S':  Compute the separation only; */
/*             = 'C':  Compute the reciprocal condition number only; */
/*             = 'E':  Compute the error bound only; */
/*             = 'A':  Compute all: the solution, separation, reciprocal */
/*                     condition number, and the error bound. */

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
/*             Specifies whether or not the original or "reduced" */
/*             Lyapunov equations should be solved, as follows: */
/*             = 'O':  Solve the original Lyapunov equations, updating */
/*                     the right-hand sides and solutions with the */
/*                     matrix U, e.g., X <-- U'*X*U; */
/*             = 'R':  Solve reduced Lyapunov equations only, without */
/*                     updating the right-hand sides and solutions. */
/*                     This means that a real Schur form T of A appears */
/*                     in the equation, instead of A. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, X, and C.  N >= 0. */

/*     SCALE   (input or output) DOUBLE PRECISION */
/*             If JOB = 'C' or JOB = 'E', SCALE is an input argument: */
/*             the scale factor, set by a Lyapunov solver. */
/*             0 <= SCALE <= 1. */
/*             If JOB = 'X' or JOB = 'A', SCALE is an output argument: */
/*             the scale factor, scale, set less than or equal to 1 to */
/*             prevent the solution overflowing. */
/*             If JOB = 'S', this argument is not used. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             If FACT = 'N' or (LYAPUN = 'O' and JOB <> 'X'), the */
/*             leading N-by-N part of this array must contain the */
/*             original matrix A. */
/*             If FACT = 'F' and (LYAPUN = 'R' or JOB = 'X'), A is */
/*             not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. */
/*             LDA >= MAX(1,N), if FACT = 'N' or LYAPUN = 'O' and */
/*                                               JOB <> 'X'; */
/*             LDA >= 1,        otherwise. */

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
/*             The contents of array T is not modified if FACT = 'F'. */

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
/*             If JOB <> 'S' and UPLO = 'U', the leading N-by-N upper */
/*             triangular part of this array must contain the upper */
/*             triangular part of the matrix C of the original Lyapunov */
/*             equation (with matrix A), if LYAPUN = 'O', or of the */
/*             reduced Lyapunov equation (with matrix T), if */
/*             LYAPUN = 'R'. */
/*             If JOB <> 'S' and UPLO = 'L', the leading N-by-N lower */
/*             triangular part of this array must contain the lower */
/*             triangular part of the matrix C of the original Lyapunov */
/*             equation (with matrix A), if LYAPUN = 'O', or of the */
/*             reduced Lyapunov equation (with matrix T), if */
/*             LYAPUN = 'R'. */
/*             The remaining strictly triangular part of this array is */
/*             used as workspace. */
/*             If JOB = 'X', then this array may be identified with X */
/*             in the call of this routine. */
/*             If JOB = 'S', the array C is not referenced. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C. */
/*             LDC >= 1,        if JOB = 'S'; */
/*             LDC >= MAX(1,N), otherwise. */

/*     X       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDX,N) */
/*             If JOB = 'C' or 'E', then X is an input argument and on */
/*             entry, the leading N-by-N part of this array must contain */
/*             the symmetric solution matrix X of the original Lyapunov */
/*             equation (with matrix A), if LYAPUN = 'O', or of the */
/*             reduced Lyapunov equation (with matrix T), if */
/*             LYAPUN = 'R'. */
/*             If JOB = 'X' or 'A', then X is an output argument and on */
/*             exit, if INFO = 0 or INFO = N+1, the leading N-by-N part */
/*             of this array contains the symmetric solution matrix X of */
/*             of the original Lyapunov equation (with matrix A), if */
/*             LYAPUN = 'O', or of the reduced Lyapunov equation (with */
/*             matrix T), if LYAPUN = 'R'. */
/*             If JOB = 'S', the array X is not referenced. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X. */
/*             LDX >= 1,        if JOB = 'S'; */
/*             LDX >= MAX(1,N), otherwise. */

/*     SEP     (output) DOUBLE PRECISION */
/*             If JOB = 'S' or JOB = 'C' or JOB = 'A', and INFO = 0 or */
/*             INFO = N+1, SEP contains the estimated separation of the */
/*             matrices op(A) and -op(A)', sep(op(A),-op(A)'). */
/*             If N = 0, or X = 0, or JOB = 'X' or JOB = 'E', SEP is not */
/*             referenced. */

/*     RCOND   (output) DOUBLE PRECISION */
/*             If JOB = 'C' or JOB = 'A', an estimate of the reciprocal */
/*             condition number of the continuous-time Lyapunov equation. */
/*             If N = 0 or X = 0, RCOND is set to 1 or 0, respectively. */
/*             If JOB = 'X' or JOB = 'S' or JOB = 'E', RCOND is not */
/*             referenced. */

/*     FERR    (output) DOUBLE PRECISION */
/*             If JOB = 'E' or JOB = 'A', and INFO = 0 or INFO = N+1, */
/*             FERR contains an estimated forward error bound for the */
/*             solution X. If XTRUE is the true solution, FERR bounds the */
/*             relative error in the computed solution, measured in the */
/*             Frobenius norm:  norm(X - XTRUE)/norm(XTRUE). */
/*             If N = 0 or X = 0, FERR is set to 0. */
/*             If JOB = 'X' or JOB = 'S' or JOB = 'C', FERR is not */
/*             referenced. */

/*     WR      (output) DOUBLE PRECISION array, dimension (N) */
/*     WI      (output) DOUBLE PRECISION array, dimension (N) */
/*             If FACT = 'N', and INFO = 0 or INFO = N+1, WR and WI */
/*             contain the real and imaginary parts, respectively, of the */
/*             eigenvalues of A. */
/*             If FACT = 'F', WR and WI are not referenced. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N*N) */
/*             This array is not referenced if JOB = 'X'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0 or INFO = N+1, DWORK(1) returns the */
/*             optimal value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             If JOB = 'X', then */
/*             LDWORK >= MAX(1,N*N),           if FACT = 'F'; */
/*             LDWORK >= MAX(1,MAX(N*N,3*N)),  if FACT = 'N'. */
/*             If JOB = 'S' or JOB = 'C', then */
/*             LDWORK >= MAX(1,2*N*N),         if FACT = 'F'; */
/*             LDWORK >= MAX(1,2*N*N,3*N),     if FACT = 'N'. */
/*             If JOB = 'E', or JOB = 'A', and LYAPUN  = 'O', then */
/*             LDWORK >= MAX(1,3*N*N); */
/*             If JOB = 'E', or JOB = 'A', and LYAPUN  = 'R', then */
/*             LDWORK >= MAX(1,3*N*N+N-1). */
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
/*                   Schur form, and the elements i+1:n of WR and WI */
/*                   contain the real and imaginary parts, respectively, */
/*                   of the converged eigenvalues; this error is unlikely */
/*                   to appear; */
/*             = N+1:  if the matrices T and -T' have common or very */
/*                   close eigenvalues; perturbed values were used to */
/*                   solve Lyapunov equations, but the matrix T, if given */
/*                   (for FACT = 'F'), is unchanged. */

/*     METHOD */

/*     After reducing matrix A to real Schur canonical form (if needed), */
/*     the Bartels-Stewart algorithm is used. A set of equivalent linear */
/*     algebraic systems of equations of order at most four are formed */
/*     and solved using Gaussian elimination with complete pivoting. */

/*     The condition number of the continuous-time Lyapunov equation is */
/*     estimated as */

/*     cond = (norm(Theta)*norm(A) + norm(inv(Omega))*norm(C))/norm(X), */

/*     where Omega and Theta are linear operators defined by */

/*     Omega(W) = op(A)'*W + W*op(A), */
/*     Theta(W) = inv(Omega(op(W)'*X + X*op(W))). */

/*     The routine estimates the quantities */

/*     sep(op(A),-op(A)') = 1 / norm(inv(Omega)) */

/*     and norm(Theta) using 1-norm condition estimators. */

/*     The forward error bound is estimated using a practical error bound */
/*     similar to the one proposed in [2]. */

/*     REFERENCES */

/*     [1] Bartels, R.H. and Stewart, G.W.  T */
/*         Solution of the matrix equation A X + XB = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     [2] Higham, N.J. */
/*         Perturbation theory and backward error for AX-XB=C. */
/*         BIT, vol. 33, pp. 124-136, 1993. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */
/*     The accuracy of the estimates obtained depends on the solution */
/*     accuracy and on the properties of the 1-norm estimator. */

/*     FURTHER COMMENTS */

/*     The separation of op(A) and -op(A)' can also be defined as */

/*            sep( op(A), -op(A)' ) = sigma_min( T ), */

/*     where sigma_min(T) is the smallest singular value of the */
/*     N*N-by-N*N matrix */

/*        T = kprod( I(N), op(A)' ) + kprod( op(A)', I(N) ). */

/*     I(N) is an N-by-N identity matrix, and kprod denotes the Kronecker */
/*     product. The routine estimates sigma_min(T) by the reciprocal of */
/*     an estimate of the 1-norm of inverse(T). The true reciprocal */
/*     1-norm of inverse(T) cannot differ from sigma_min(T) by more */
/*     than a factor of N. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, February 1999. */
/*     This is an extended and improved version of Release 3.0 routine */
/*     SB03RD. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2004. */

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

/*     Decode option parameters. */

#line 375 "SB03TD.f"
    /* Parameter adjustments */
#line 375 "SB03TD.f"
    a_dim1 = *lda;
#line 375 "SB03TD.f"
    a_offset = 1 + a_dim1;
#line 375 "SB03TD.f"
    a -= a_offset;
#line 375 "SB03TD.f"
    t_dim1 = *ldt;
#line 375 "SB03TD.f"
    t_offset = 1 + t_dim1;
#line 375 "SB03TD.f"
    t -= t_offset;
#line 375 "SB03TD.f"
    u_dim1 = *ldu;
#line 375 "SB03TD.f"
    u_offset = 1 + u_dim1;
#line 375 "SB03TD.f"
    u -= u_offset;
#line 375 "SB03TD.f"
    c_dim1 = *ldc;
#line 375 "SB03TD.f"
    c_offset = 1 + c_dim1;
#line 375 "SB03TD.f"
    c__ -= c_offset;
#line 375 "SB03TD.f"
    x_dim1 = *ldx;
#line 375 "SB03TD.f"
    x_offset = 1 + x_dim1;
#line 375 "SB03TD.f"
    x -= x_offset;
#line 375 "SB03TD.f"
    --wr;
#line 375 "SB03TD.f"
    --wi;
#line 375 "SB03TD.f"
    --iwork;
#line 375 "SB03TD.f"
    --dwork;
#line 375 "SB03TD.f"

#line 375 "SB03TD.f"
    /* Function Body */
#line 375 "SB03TD.f"
    jobx = lsame_(job, "X", (ftnlen)1, (ftnlen)1);
#line 376 "SB03TD.f"
    jobs = lsame_(job, "S", (ftnlen)1, (ftnlen)1);
#line 377 "SB03TD.f"
    jobc = lsame_(job, "C", (ftnlen)1, (ftnlen)1);
#line 378 "SB03TD.f"
    jobe = lsame_(job, "E", (ftnlen)1, (ftnlen)1);
#line 379 "SB03TD.f"
    joba = lsame_(job, "A", (ftnlen)1, (ftnlen)1);
#line 380 "SB03TD.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 381 "SB03TD.f"
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
#line 382 "SB03TD.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 383 "SB03TD.f"
    update = lsame_(lyapun, "O", (ftnlen)1, (ftnlen)1);

/*     Compute workspace. */

#line 387 "SB03TD.f"
    nn = *n * *n;
#line 388 "SB03TD.f"
    if (jobx) {
#line 389 "SB03TD.f"
	ldw = nn;
#line 390 "SB03TD.f"
    } else if (jobs || jobc) {
#line 391 "SB03TD.f"
	ldw = nn << 1;
#line 392 "SB03TD.f"
    } else {
#line 393 "SB03TD.f"
	ldw = nn * 3;
#line 394 "SB03TD.f"
    }
#line 395 "SB03TD.f"
    if ((jobe || joba) && ! update) {
#line 395 "SB03TD.f"
	ldw = ldw + *n - 1;
#line 395 "SB03TD.f"
    }
#line 397 "SB03TD.f"
    if (nofact) {
/* Computing MAX */
#line 397 "SB03TD.f"
	i__1 = ldw, i__2 = *n * 3;
#line 397 "SB03TD.f"
	ldw = max(i__1,i__2);
#line 397 "SB03TD.f"
    }

/*     Test the scalar input parameters. */

#line 402 "SB03TD.f"
    *info = 0;
#line 403 "SB03TD.f"
    if (! (jobx || jobs || jobc || jobe || joba)) {
#line 404 "SB03TD.f"
	*info = -1;
#line 405 "SB03TD.f"
    } else if (! (nofact || lsame_(fact, "F", (ftnlen)1, (ftnlen)1))) {
#line 406 "SB03TD.f"
	*info = -2;
#line 407 "SB03TD.f"
    } else if (! (notrna || lsame_(trana, "T", (ftnlen)1, (ftnlen)1) || 
	    lsame_(trana, "C", (ftnlen)1, (ftnlen)1))) {
#line 409 "SB03TD.f"
	*info = -3;
#line 410 "SB03TD.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 411 "SB03TD.f"
	*info = -4;
#line 412 "SB03TD.f"
    } else if (! (update || lsame_(lyapun, "R", (ftnlen)1, (ftnlen)1))) {
#line 413 "SB03TD.f"
	*info = -5;
#line 414 "SB03TD.f"
    } else if (*n < 0) {
#line 415 "SB03TD.f"
	*info = -6;
#line 416 "SB03TD.f"
    } else if ((jobc || jobe) && (*scale < 0. || *scale > 1.)) {
#line 418 "SB03TD.f"
	*info = -7;
#line 419 "SB03TD.f"
    } else if (*lda < 1 || *lda < *n && (update && ! jobx || nofact)) {
#line 422 "SB03TD.f"
	*info = -9;
#line 423 "SB03TD.f"
    } else if (*ldt < max(1,*n)) {
#line 424 "SB03TD.f"
	*info = -11;
#line 425 "SB03TD.f"
    } else if (*ldu < 1 || *ldu < *n && update) {
#line 426 "SB03TD.f"
	*info = -13;
#line 427 "SB03TD.f"
    } else if (*ldc < 1 || ! jobs && *ldc < *n) {
#line 428 "SB03TD.f"
	*info = -15;
#line 429 "SB03TD.f"
    } else if (*ldx < 1 || ! jobs && *ldx < *n) {
#line 430 "SB03TD.f"
	*info = -17;
#line 431 "SB03TD.f"
    } else if (*ldwork < 1 || *ldwork < ldw) {
#line 432 "SB03TD.f"
	*info = -25;
#line 433 "SB03TD.f"
    }

#line 435 "SB03TD.f"
    if (*info != 0) {
#line 436 "SB03TD.f"
	i__1 = -(*info);
#line 436 "SB03TD.f"
	xerbla_("SB03TD", &i__1, (ftnlen)6);
#line 437 "SB03TD.f"
	return 0;
#line 438 "SB03TD.f"
    }

/*     Quick return if possible. */

#line 442 "SB03TD.f"
    if (*n == 0) {
#line 443 "SB03TD.f"
	if (jobx || joba) {
#line 443 "SB03TD.f"
	    *scale = 1.;
#line 443 "SB03TD.f"
	}
#line 445 "SB03TD.f"
	if (jobc || joba) {
#line 445 "SB03TD.f"
	    *rcond = 1.;
#line 445 "SB03TD.f"
	}
#line 447 "SB03TD.f"
	if (jobe || joba) {
#line 447 "SB03TD.f"
	    *ferr = 0.;
#line 447 "SB03TD.f"
	}
#line 449 "SB03TD.f"
	dwork[1] = 1.;
#line 450 "SB03TD.f"
	return 0;
#line 451 "SB03TD.f"
    }

#line 453 "SB03TD.f"
    if (nofact) {

/*        Compute the Schur factorization of A. */
/*        Workspace:  need   3*N; */
/*                    prefer larger. */

#line 459 "SB03TD.f"
	dlacpy_("Full", n, n, &a[a_offset], lda, &t[t_offset], ldt, (ftnlen)4)
		;
#line 460 "SB03TD.f"
	if (update) {
#line 461 "SB03TD.f"
	    *(unsigned char *)sjob = 'V';
#line 462 "SB03TD.f"
	} else {
#line 463 "SB03TD.f"
	    *(unsigned char *)sjob = 'N';
#line 464 "SB03TD.f"
	}
#line 465 "SB03TD.f"
	dgees_(sjob, "Not ordered", (L_fp)select_, n, &t[t_offset], ldt, &
		sdim, &wr[1], &wi[1], &u[u_offset], ldu, &dwork[1], ldwork, 
		bwork, info, (ftnlen)1, (ftnlen)11);
#line 467 "SB03TD.f"
	if (*info > 0) {
#line 467 "SB03TD.f"
	    return 0;
#line 467 "SB03TD.f"
	}
#line 469 "SB03TD.f"
	*(unsigned char *)cfact = 'F';
#line 470 "SB03TD.f"
    } else {
#line 471 "SB03TD.f"
	*(unsigned char *)cfact = *(unsigned char *)fact;
#line 472 "SB03TD.f"
    }

#line 474 "SB03TD.f"
    if (jobx || joba) {

/*        Copy the right-hand side in X. */

#line 478 "SB03TD.f"
	dlacpy_(uplo, n, n, &c__[c_offset], ldc, &x[x_offset], ldx, (ftnlen)1)
		;

#line 480 "SB03TD.f"
	if (update) {

/*           Transform the right-hand side. */
/*           Workspace:  need   N*N. */

#line 485 "SB03TD.f"
	    mb01ru_(uplo, "Transpose", n, n, &c_b23, &c_b24, &x[x_offset], 
		    ldx, &u[u_offset], ldu, &x[x_offset], ldx, &dwork[1], 
		    ldwork, info, (ftnlen)1, (ftnlen)9);
#line 487 "SB03TD.f"
	    i__1 = *ldx + 1;
#line 487 "SB03TD.f"
	    dscal_(n, &c_b25, &x[x_offset], &i__1);
#line 488 "SB03TD.f"
	}

/*        Fill in the remaining triangle of X. */

#line 492 "SB03TD.f"
	ma02ed_(uplo, n, &x[x_offset], ldx, (ftnlen)1);

/*        Solve the transformed equation. */

#line 496 "SB03TD.f"
	sb03my_(trana, n, &t[t_offset], ldt, &x[x_offset], ldx, scale, info, (
		ftnlen)1);
#line 497 "SB03TD.f"
	if (*info > 0) {
#line 497 "SB03TD.f"
	    *info = *n + 1;
#line 497 "SB03TD.f"
	}

#line 500 "SB03TD.f"
	if (update) {

/*           Transform back the solution. */

#line 504 "SB03TD.f"
	    mb01ru_(uplo, "No transpose", n, n, &c_b23, &c_b24, &x[x_offset], 
		    ldx, &u[u_offset], ldu, &x[x_offset], ldx, &dwork[1], 
		    ldwork, info, (ftnlen)1, (ftnlen)12);
#line 506 "SB03TD.f"
	    i__1 = *ldx + 1;
#line 506 "SB03TD.f"
	    dscal_(n, &c_b25, &x[x_offset], &i__1);

/*           Fill in the remaining triangle of X. */

#line 510 "SB03TD.f"
	    ma02ed_(uplo, n, &x[x_offset], ldx, (ftnlen)1);
#line 511 "SB03TD.f"
	}
#line 512 "SB03TD.f"
    }

#line 514 "SB03TD.f"
    if (jobs) {

/*        Estimate sep(op(A),-op(A)'). */
/*        Workspace:  2*N*N. */

#line 519 "SB03TD.f"
	sb03qy_("Separation", trana, lyapun, n, &t[t_offset], ldt, &u[
		u_offset], ldu, &x[x_offset], ldx, sep, &thnorm, &iwork[1], &
		dwork[1], ldwork, info, (ftnlen)10, (ftnlen)1, (ftnlen)1);

#line 522 "SB03TD.f"
    } else if (! jobx) {

/*        Estimate the reciprocal condition and/or the error bound. */
/*        Workspace:  2*N*N, if JOB = 'C'; */
/*                    3*N*N + a*(N-1), where: */
/*                    a = 1, if JOB = 'E' or JOB = 'A', and LYAPUN = 'R'; */
/*                    a = 0, otherwise. */

#line 530 "SB03TD.f"
	if (joba) {
#line 531 "SB03TD.f"
	    *(unsigned char *)jobl = 'B';
#line 532 "SB03TD.f"
	} else {
#line 533 "SB03TD.f"
	    *(unsigned char *)jobl = *(unsigned char *)job;
#line 534 "SB03TD.f"
	}
#line 535 "SB03TD.f"
	sb03qd_(jobl, cfact, trana, uplo, lyapun, n, scale, &a[a_offset], lda,
		 &t[t_offset], ldt, &u[u_offset], ldu, &c__[c_offset], ldc, &
		x[x_offset], ldx, sep, rcond, ferr, &iwork[1], &dwork[1], 
		ldwork, info, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
/* Computing MAX */
#line 538 "SB03TD.f"
	i__1 = ldw, i__2 = (integer) dwork[1];
#line 538 "SB03TD.f"
	ldw = max(i__1,i__2);
#line 539 "SB03TD.f"
    }

#line 541 "SB03TD.f"
    dwork[1] = (doublereal) ldw;

#line 543 "SB03TD.f"
    return 0;
/* *** Last line of SB03TD *** */
} /* sb03td_ */

