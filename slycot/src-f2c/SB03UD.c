#line 1 "SB03UD.f"
/* SB03UD.f -- translated by f2c (version 20100827).
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

#line 1 "SB03UD.f"
/* Table of constant values */

static doublereal c_b23 = 0.;
static doublereal c_b24 = 1.;
static doublereal c_b25 = .5;
static integer c__1 = 1;

/* Subroutine */ int sb03ud_(char *job, char *fact, char *trana, char *uplo, 
	char *lyapun, integer *n, doublereal *scale, doublereal *a, integer *
	lda, doublereal *t, integer *ldt, doublereal *u, integer *ldu, 
	doublereal *c__, integer *ldc, doublereal *x, integer *ldx, 
	doublereal *sepd, doublereal *rcond, doublereal *ferr, doublereal *wr,
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
    static logical jobs;
    static char sjob[1];
    static integer sdim;
    static logical jobx;
    extern /* Subroutine */ int ma02ed_(char *, integer *, doublereal *, 
	    integer *, ftnlen);
    static char cfact[1];
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgees_(char *, char *, L_fp, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, logical *, integer *, ftnlen, 
	    ftnlen), sb03sd_(char *, char *, char *, char *, char *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb01ru_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), sb03mx_(char *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen), sb03sy_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
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

/*     To solve the real discrete-time Lyapunov matrix equation */

/*            op(A)'*X*op(A) - X = scale*C, */

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

/*     SEPD    (output) DOUBLE PRECISION */
/*             If JOB = 'S' or JOB = 'C' or JOB = 'A', and INFO = 0 or */
/*             INFO = N+1, SEPD contains the estimated separation of the */
/*             matrices op(A) and op(A)', sepd(op(A),op(A)'). */
/*             If N = 0, or X = 0, or JOB = 'X' or JOB = 'E', SEPD is not */
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
/*             LDWORK >= MAX(1,N*N,2*N),       if FACT = 'F'; */
/*             LDWORK >= MAX(1,N*N,3*N),       if FACT = 'N'. */
/*             If JOB = 'S', then */
/*             LDWORK >= MAX(3,2*N*N). */
/*             If JOB = 'C', then */
/*             LDWORK >= MAX(3,2*N*N) + N*N. */
/*             If JOB = 'E', or JOB = 'A', then */
/*             LDWORK >= MAX(3,2*N*N) + N*N + 2*N. */
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
/*             = N+1:  if the matrix T has almost reciprocal eigenvalues; */
/*                   perturbed values were used to solve Lyapunov */
/*                   equations, but the matrix T, if given (for */
/*                   FACT = 'F'), is unchanged. */

/*     METHOD */

/*     After reducing matrix A to real Schur canonical form (if needed), */
/*     a discrete-time version of the Bartels-Stewart algorithm is used. */
/*     A set of equivalent linear algebraic systems of equations of order */
/*     at most four are formed and solved using Gaussian elimination with */
/*     complete pivoting. */

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
/*     similar to the one proposed in [3]. */

/*     REFERENCES */

/*     [1] Barraud, A.Y.                   T */
/*         A numerical algorithm to solve A XA - X = Q. */
/*         IEEE Trans. Auto. Contr., AC-22, pp. 883-885, 1977. */

/*     [2] Bartels, R.H. and Stewart, G.W.  T */
/*         Solution of the matrix equation A X + XB = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     [3] Higham, N.J. */
/*         Perturbation theory and backward error for AX-XB=C. */
/*         BIT, vol. 33, pp. 124-136, 1993. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */
/*     The accuracy of the estimates obtained depends on the solution */
/*     accuracy and on the properties of the 1-norm estimator. */

/*     FURTHER COMMENTS */

/*     The "separation" sepd of op(A) and op(A)' can also be defined as */

/*            sepd( op(A), op(A)' ) = sigma_min( T ), */

/*     where sigma_min(T) is the smallest singular value of the */
/*     N*N-by-N*N matrix */

/*        T = kprod( op(A)', op(A)' ) - I(N**2). */

/*     I(N**2) is an N*N-by-N*N identity matrix, and kprod denotes the */
/*     Kronecker product. The routine estimates sigma_min(T) by the */
/*     reciprocal of an estimate of the 1-norm of inverse(T). The true */
/*     reciprocal 1-norm of inverse(T) cannot differ from sigma_min(T) by */
/*     more than a factor of N. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, February 1999. */
/*     This is an extended and improved version of Release 3.0 routine */
/*     SB03PD. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2004. */

/*     KEYWORDS */

/*     Lyapunov equation, orthogonal transformation, real Schur form, */
/*     Sylvester equation. */

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

#line 380 "SB03UD.f"
    /* Parameter adjustments */
#line 380 "SB03UD.f"
    a_dim1 = *lda;
#line 380 "SB03UD.f"
    a_offset = 1 + a_dim1;
#line 380 "SB03UD.f"
    a -= a_offset;
#line 380 "SB03UD.f"
    t_dim1 = *ldt;
#line 380 "SB03UD.f"
    t_offset = 1 + t_dim1;
#line 380 "SB03UD.f"
    t -= t_offset;
#line 380 "SB03UD.f"
    u_dim1 = *ldu;
#line 380 "SB03UD.f"
    u_offset = 1 + u_dim1;
#line 380 "SB03UD.f"
    u -= u_offset;
#line 380 "SB03UD.f"
    c_dim1 = *ldc;
#line 380 "SB03UD.f"
    c_offset = 1 + c_dim1;
#line 380 "SB03UD.f"
    c__ -= c_offset;
#line 380 "SB03UD.f"
    x_dim1 = *ldx;
#line 380 "SB03UD.f"
    x_offset = 1 + x_dim1;
#line 380 "SB03UD.f"
    x -= x_offset;
#line 380 "SB03UD.f"
    --wr;
#line 380 "SB03UD.f"
    --wi;
#line 380 "SB03UD.f"
    --iwork;
#line 380 "SB03UD.f"
    --dwork;
#line 380 "SB03UD.f"

#line 380 "SB03UD.f"
    /* Function Body */
#line 380 "SB03UD.f"
    jobx = lsame_(job, "X", (ftnlen)1, (ftnlen)1);
#line 381 "SB03UD.f"
    jobs = lsame_(job, "S", (ftnlen)1, (ftnlen)1);
#line 382 "SB03UD.f"
    jobc = lsame_(job, "C", (ftnlen)1, (ftnlen)1);
#line 383 "SB03UD.f"
    jobe = lsame_(job, "E", (ftnlen)1, (ftnlen)1);
#line 384 "SB03UD.f"
    joba = lsame_(job, "A", (ftnlen)1, (ftnlen)1);
#line 385 "SB03UD.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 386 "SB03UD.f"
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
#line 387 "SB03UD.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 388 "SB03UD.f"
    update = lsame_(lyapun, "O", (ftnlen)1, (ftnlen)1);

/*     Compute workspace. */

#line 392 "SB03UD.f"
    nn = *n * *n;
#line 393 "SB03UD.f"
    if (jobx) {
#line 394 "SB03UD.f"
	if (nofact) {
/* Computing MAX */
#line 395 "SB03UD.f"
	    i__1 = max(1,nn), i__2 = *n * 3;
#line 395 "SB03UD.f"
	    ldw = max(i__1,i__2);
#line 396 "SB03UD.f"
	} else {
/* Computing MAX */
#line 397 "SB03UD.f"
	    i__1 = max(1,nn), i__2 = *n << 1;
#line 397 "SB03UD.f"
	    ldw = max(i__1,i__2);
#line 398 "SB03UD.f"
	}
#line 399 "SB03UD.f"
    } else if (jobs) {
/* Computing MAX */
#line 400 "SB03UD.f"
	i__1 = 3, i__2 = nn << 1;
#line 400 "SB03UD.f"
	ldw = max(i__1,i__2);
#line 401 "SB03UD.f"
    } else if (jobc) {
/* Computing MAX */
#line 402 "SB03UD.f"
	i__1 = 3, i__2 = nn << 1;
#line 402 "SB03UD.f"
	ldw = max(i__1,i__2) + nn;
#line 403 "SB03UD.f"
    } else {
/* Computing MAX */
#line 404 "SB03UD.f"
	i__1 = 3, i__2 = nn << 1;
#line 404 "SB03UD.f"
	ldw = max(i__1,i__2) + nn + (*n << 1);
#line 405 "SB03UD.f"
    }

/*     Test the scalar input parameters. */

#line 409 "SB03UD.f"
    *info = 0;
#line 410 "SB03UD.f"
    if (! (jobx || jobs || jobc || jobe || joba)) {
#line 411 "SB03UD.f"
	*info = -1;
#line 412 "SB03UD.f"
    } else if (! (nofact || lsame_(fact, "F", (ftnlen)1, (ftnlen)1))) {
#line 413 "SB03UD.f"
	*info = -2;
#line 414 "SB03UD.f"
    } else if (! (notrna || lsame_(trana, "T", (ftnlen)1, (ftnlen)1) || 
	    lsame_(trana, "C", (ftnlen)1, (ftnlen)1))) {
#line 416 "SB03UD.f"
	*info = -3;
#line 417 "SB03UD.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 418 "SB03UD.f"
	*info = -4;
#line 419 "SB03UD.f"
    } else if (! (update || lsame_(lyapun, "R", (ftnlen)1, (ftnlen)1))) {
#line 420 "SB03UD.f"
	*info = -5;
#line 421 "SB03UD.f"
    } else if (*n < 0) {
#line 422 "SB03UD.f"
	*info = -6;
#line 423 "SB03UD.f"
    } else if ((jobc || jobe) && (*scale < 0. || *scale > 1.)) {
#line 425 "SB03UD.f"
	*info = -7;
#line 426 "SB03UD.f"
    } else if (*lda < 1 || *lda < *n && (update && ! jobx || nofact)) {
#line 429 "SB03UD.f"
	*info = -9;
#line 430 "SB03UD.f"
    } else if (*ldt < max(1,*n)) {
#line 431 "SB03UD.f"
	*info = -11;
#line 432 "SB03UD.f"
    } else if (*ldu < 1 || *ldu < *n && update) {
#line 433 "SB03UD.f"
	*info = -13;
#line 434 "SB03UD.f"
    } else if (*ldc < 1 || ! jobs && *ldc < *n) {
#line 435 "SB03UD.f"
	*info = -15;
#line 436 "SB03UD.f"
    } else if (*ldx < 1 || ! jobs && *ldx < *n) {
#line 437 "SB03UD.f"
	*info = -17;
#line 438 "SB03UD.f"
    } else if (*ldwork < ldw) {
#line 439 "SB03UD.f"
	*info = -25;
#line 440 "SB03UD.f"
    }

#line 442 "SB03UD.f"
    if (*info != 0) {
#line 443 "SB03UD.f"
	i__1 = -(*info);
#line 443 "SB03UD.f"
	xerbla_("SB03UD", &i__1, (ftnlen)6);
#line 444 "SB03UD.f"
	return 0;
#line 445 "SB03UD.f"
    }

/*     Quick return if possible. */

#line 449 "SB03UD.f"
    if (*n == 0) {
#line 450 "SB03UD.f"
	if (jobx || joba) {
#line 450 "SB03UD.f"
	    *scale = 1.;
#line 450 "SB03UD.f"
	}
#line 452 "SB03UD.f"
	if (jobc || joba) {
#line 452 "SB03UD.f"
	    *rcond = 1.;
#line 452 "SB03UD.f"
	}
#line 454 "SB03UD.f"
	if (jobe || joba) {
#line 454 "SB03UD.f"
	    *ferr = 0.;
#line 454 "SB03UD.f"
	}
#line 456 "SB03UD.f"
	dwork[1] = 1.;
#line 457 "SB03UD.f"
	return 0;
#line 458 "SB03UD.f"
    }

#line 460 "SB03UD.f"
    if (nofact) {

/*        Compute the Schur factorization of A. */
/*        Workspace:  need   3*N; */
/*                    prefer larger. */

#line 466 "SB03UD.f"
	dlacpy_("Full", n, n, &a[a_offset], lda, &t[t_offset], ldt, (ftnlen)4)
		;
#line 467 "SB03UD.f"
	if (update) {
#line 468 "SB03UD.f"
	    *(unsigned char *)sjob = 'V';
#line 469 "SB03UD.f"
	} else {
#line 470 "SB03UD.f"
	    *(unsigned char *)sjob = 'N';
#line 471 "SB03UD.f"
	}
#line 472 "SB03UD.f"
	dgees_(sjob, "Not ordered", (L_fp)select_, n, &t[t_offset], ldt, &
		sdim, &wr[1], &wi[1], &u[u_offset], ldu, &dwork[1], ldwork, 
		bwork, info, (ftnlen)1, (ftnlen)11);
#line 474 "SB03UD.f"
	if (*info > 0) {
#line 474 "SB03UD.f"
	    return 0;
#line 474 "SB03UD.f"
	}
/* Computing MAX */
#line 476 "SB03UD.f"
	i__1 = ldw, i__2 = (integer) dwork[1];
#line 476 "SB03UD.f"
	ldw = max(i__1,i__2);
#line 477 "SB03UD.f"
	*(unsigned char *)cfact = 'F';
#line 478 "SB03UD.f"
    } else {
#line 479 "SB03UD.f"
	*(unsigned char *)cfact = *(unsigned char *)fact;
#line 480 "SB03UD.f"
    }

#line 482 "SB03UD.f"
    if (jobx || joba) {

/*        Copy the right-hand side in X. */

#line 486 "SB03UD.f"
	dlacpy_(uplo, n, n, &c__[c_offset], ldc, &x[x_offset], ldx, (ftnlen)1)
		;

#line 488 "SB03UD.f"
	if (update) {

/*           Transform the right-hand side. */
/*           Workspace:  need   N*N. */

#line 493 "SB03UD.f"
	    mb01ru_(uplo, "Transpose", n, n, &c_b23, &c_b24, &x[x_offset], 
		    ldx, &u[u_offset], ldu, &x[x_offset], ldx, &dwork[1], 
		    ldwork, info, (ftnlen)1, (ftnlen)9);
#line 495 "SB03UD.f"
	    i__1 = *ldx + 1;
#line 495 "SB03UD.f"
	    dscal_(n, &c_b25, &x[x_offset], &i__1);
#line 496 "SB03UD.f"
	}

/*        Fill in the remaining triangle of X. */

#line 500 "SB03UD.f"
	ma02ed_(uplo, n, &x[x_offset], ldx, (ftnlen)1);

/*        Solve the transformed equation. */
/*        Workspace:  2*N. */

#line 505 "SB03UD.f"
	sb03mx_(trana, n, &t[t_offset], ldt, &x[x_offset], ldx, scale, &dwork[
		1], info, (ftnlen)1);
#line 506 "SB03UD.f"
	if (*info > 0) {
#line 506 "SB03UD.f"
	    *info = *n + 1;
#line 506 "SB03UD.f"
	}

#line 509 "SB03UD.f"
	if (update) {

/*           Transform back the solution. */

#line 513 "SB03UD.f"
	    mb01ru_(uplo, "No transpose", n, n, &c_b23, &c_b24, &x[x_offset], 
		    ldx, &u[u_offset], ldu, &x[x_offset], ldx, &dwork[1], 
		    ldwork, info, (ftnlen)1, (ftnlen)12);
#line 515 "SB03UD.f"
	    i__1 = *ldx + 1;
#line 515 "SB03UD.f"
	    dscal_(n, &c_b25, &x[x_offset], &i__1);

/*           Fill in the remaining triangle of X. */

#line 519 "SB03UD.f"
	    ma02ed_(uplo, n, &x[x_offset], ldx, (ftnlen)1);
#line 520 "SB03UD.f"
	}
#line 521 "SB03UD.f"
    }

#line 523 "SB03UD.f"
    if (jobs) {

/*        Estimate sepd(op(A),op(A)'). */
/*        Workspace:  MAX(3,2*N*N). */

#line 528 "SB03UD.f"
	sb03sy_("Separation", trana, lyapun, n, &t[t_offset], ldt, &u[
		u_offset], ldu, &dwork[1], &c__1, sepd, &thnorm, &iwork[1], &
		dwork[1], ldwork, info, (ftnlen)10, (ftnlen)1, (ftnlen)1);

#line 532 "SB03UD.f"
    } else if (! jobx) {

/*        Estimate the reciprocal condition and/or the error bound. */
/*        Workspace:  MAX(3,2*N*N) + N*N + a*N, where: */
/*                    a = 2, if JOB = 'E' or JOB = 'A'; */
/*                    a = 0, otherwise. */

#line 539 "SB03UD.f"
	if (joba) {
#line 540 "SB03UD.f"
	    *(unsigned char *)jobl = 'B';
#line 541 "SB03UD.f"
	} else {
#line 542 "SB03UD.f"
	    *(unsigned char *)jobl = *(unsigned char *)job;
#line 543 "SB03UD.f"
	}
#line 544 "SB03UD.f"
	sb03sd_(jobl, cfact, trana, uplo, lyapun, n, scale, &a[a_offset], lda,
		 &t[t_offset], ldt, &u[u_offset], ldu, &c__[c_offset], ldc, &
		x[x_offset], ldx, sepd, rcond, ferr, &iwork[1], &dwork[1], 
		ldwork, info, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
/* Computing MAX */
#line 547 "SB03UD.f"
	i__1 = ldw, i__2 = (integer) dwork[1];
#line 547 "SB03UD.f"
	ldw = max(i__1,i__2);
#line 548 "SB03UD.f"
    }

#line 550 "SB03UD.f"
    dwork[1] = (doublereal) ldw;

#line 552 "SB03UD.f"
    return 0;
/* *** Last line of SB03UD *** */
} /* sb03ud_ */

