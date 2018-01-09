#line 1 "SB02QD.f"
/* SB02QD.f -- translated by f2c (version 20100827).
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

#line 1 "SB02QD.f"
/* Table of constant values */

static doublereal c_b18 = -1.;
static doublereal c_b19 = 1.;
static integer c__1 = 1;
static doublereal c_b41 = 0.;
static doublereal c_b43 = .5;

/* Subroutine */ int sb02qd_(char *job, char *fact, char *trana, char *uplo, 
	char *lyapun, integer *n, doublereal *a, integer *lda, doublereal *t, 
	integer *ldt, doublereal *u, integer *ldu, doublereal *g, integer *
	ldg, doublereal *q, integer *ldq, doublereal *x, integer *ldx, 
	doublereal *sep, doublereal *rcond, doublereal *ferr, integer *iwork, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen job_len, 
	ftnlen fact_len, ftnlen trana_len, ftnlen uplo_len, ftnlen lyapun_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, g_dim1, g_offset, q_dim1, q_offset, t_dim1, 
	    t_offset, u_dim1, u_offset, x_dim1, x_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, jj, nn, jx;
    static doublereal sig;
    static integer lwa, ldw;
    static doublereal eps, est;
    static logical jobb, jobc, jobe;
    static integer iabs, kase, sdim;
    static char sjob[1];
    static integer ires, ixbs;
    static doublereal epsn, temp;
    static integer itmp;
    static doublereal tmax;
    static char loup[1];
    static integer info2;
    extern /* Subroutine */ int ma02ed_(char *, integer *, doublereal *, 
	    integer *, ftnlen), dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal scale;
    extern /* Subroutine */ int dgees_(char *, char *, L_fp, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen), mb01ud_(char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal denom;
    extern /* Subroutine */ int mb01ru_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int sb03my_(char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), dcopy_(integer *, doublereal *, integer *, doublereal *, 
	    integer *), sb03qx_(char *, char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), sb03qy_(char *, char *, char *
	    , integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    static doublereal gnorm;
    static logical bwork[1];
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical lower;
    extern /* Subroutine */ int dsymm_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal qnorm, xnorm;
    extern /* Subroutine */ int dsyr2k_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static logical needac;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlacon_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *);
    extern doublereal dlanhs_(char *, integer *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    static logical nofact;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern logical select_();
    static doublereal bignum;
    static logical update;
    static char tranat[1];
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static logical notrna;
    static doublereal pinorm, xanorm, thnorm;
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
/*     solution of the real continuous-time matrix algebraic Riccati */
/*     equation */

/*         op(A)'*X + X*op(A) + Q - X*G*X = 0,                        (1) */

/*     where op(A) = A or A' (A**T) and Q, G are symmetric (Q = Q**T, */
/*     G = G**T). The matrices A, Q and G are N-by-N and the solution X */
/*     is N-by-N. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the computation to be performed, as follows: */
/*             = 'C':  Compute the reciprocal condition number only; */
/*             = 'E':  Compute the error bound only; */
/*             = 'B':  Compute both the reciprocal condition number and */
/*                     the error bound. */

/*     FACT    CHARACTER*1 */
/*             Specifies whether or not the real Schur factorization of */
/*             the matrix Ac = A - G*X (if TRANA = 'N') or Ac = A - X*G */
/*             (if TRANA = 'T' or 'C') is supplied on entry, as follows: */
/*             = 'F':  On entry, T and U (if LYAPUN = 'O') contain the */
/*                     factors from the real Schur factorization of the */
/*                     matrix Ac; */
/*             = 'N':  The Schur factorization of Ac will be computed */
/*                     and the factors will be stored in T and U (if */
/*                     LYAPUN = 'O'). */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op(A) to be used, as follows: */
/*             = 'N':  op(A) = A    (No transpose); */
/*             = 'T':  op(A) = A**T (Transpose); */
/*             = 'C':  op(A) = A**T (Conjugate transpose = Transpose). */

/*     UPLO    CHARACTER*1 */
/*             Specifies which part of the symmetric matrices Q and G is */
/*             to be used, as follows: */
/*             = 'U':  Upper triangular part; */
/*             = 'L':  Lower triangular part. */

/*     LYAPUN  CHARACTER*1 */
/*             Specifies whether or not the original Lyapunov equations */
/*             should be solved in the iterative estimation process, */
/*             as follows: */
/*             = 'O':  Solve the original Lyapunov equations, updating */
/*                     the right-hand sides and solutions with the */
/*                     matrix U, e.g., RHS <-- U'*RHS*U; */
/*             = 'R':  Solve reduced Lyapunov equations only, without */
/*                     updating the right-hand sides and solutions. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, X, Q, and G.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             If FACT = 'N' or LYAPUN = 'O', the leading N-by-N part of */
/*             this array must contain the matrix A. */
/*             If FACT = 'F' and LYAPUN = 'R', A is not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. */
/*             LDA >= max(1,N), if FACT = 'N' or  LYAPUN = 'O'; */
/*             LDA >= 1,        if FACT = 'F' and LYAPUN = 'R'. */

/*     T       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDT,N) */
/*             If FACT = 'F', then T is an input argument and on entry, */
/*             the leading N-by-N upper Hessenberg part of this array */
/*             must contain the upper quasi-triangular matrix T in Schur */
/*             canonical form from a Schur factorization of Ac (see */
/*             argument FACT). */
/*             If FACT = 'N', then T is an output argument and on exit, */
/*             if INFO = 0 or INFO = N+1, the leading N-by-N upper */
/*             Hessenberg part of this array contains the upper quasi- */
/*             triangular matrix T in Schur canonical form from a Schur */
/*             factorization of Ac (see argument FACT). */

/*     LDT     INTEGER */
/*             The leading dimension of the array T.  LDT >= max(1,N). */

/*     U       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDU,N) */
/*             If LYAPUN = 'O' and FACT = 'F', then U is an input */
/*             argument and on entry, the leading N-by-N part of this */
/*             array must contain the orthogonal matrix U from a real */
/*             Schur factorization of Ac (see argument FACT). */
/*             If LYAPUN = 'O' and FACT = 'N', then U is an output */
/*             argument and on exit, if INFO = 0 or INFO = N+1, it */
/*             contains the orthogonal N-by-N matrix from a real Schur */
/*             factorization of Ac (see argument FACT). */
/*             If LYAPUN = 'R', the array U is not referenced. */

/*     LDU     INTEGER */
/*             The leading dimension of the array U. */
/*             LDU >= 1,        if LYAPUN = 'R'; */
/*             LDU >= MAX(1,N), if LYAPUN = 'O'. */

/*     G       (input) DOUBLE PRECISION array, dimension (LDG,N) */
/*             If UPLO = 'U', the leading N-by-N upper triangular part of */
/*             this array must contain the upper triangular part of the */
/*             matrix G. */
/*             If UPLO = 'L', the leading N-by-N lower triangular part of */
/*             this array must contain the lower triangular part of the */
/*             matrix G.                     _ */
/*             Matrix G should correspond to G in the "reduced" Riccati */
/*             equation (with matrix T, instead of A), if LYAPUN = 'R'. */
/*             See METHOD. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G.  LDG >= max(1,N). */

/*     Q       (input) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             If UPLO = 'U', the leading N-by-N upper triangular part of */
/*             this array must contain the upper triangular part of the */
/*             matrix Q. */
/*             If UPLO = 'L', the leading N-by-N lower triangular part of */
/*             this array must contain the lower triangular part of the */
/*             matrix Q.                     _ */
/*             Matrix Q should correspond to Q in the "reduced" Riccati */
/*             equation (with matrix T, instead of A), if LYAPUN = 'R'. */
/*             See METHOD. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q.  LDQ >= max(1,N). */

/*     X       (input) DOUBLE PRECISION array, dimension (LDX,N) */
/*             The leading N-by-N part of this array must contain the */
/*             symmetric solution matrix of the original Riccati */
/*             equation (with matrix A), if LYAPUN = 'O', or of the */
/*             "reduced" Riccati equation (with matrix T), if */
/*             LYAPUN = 'R'. See METHOD. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X.  LDX >= max(1,N). */

/*     SEP     (output) DOUBLE PRECISION */
/*             If JOB = 'C' or JOB = 'B', the estimated quantity */
/*             sep(op(Ac),-op(Ac)'). */
/*             If N = 0, or X = 0, or JOB = 'E', SEP is not referenced. */

/*     RCOND   (output) DOUBLE PRECISION */
/*             If JOB = 'C' or JOB = 'B', an estimate of the reciprocal */
/*             condition number of the continuous-time Riccati equation. */
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
/*             Let LWA = N*N, if LYAPUN = 'O' and JOB = 'E' or 'B'; */
/*                 LWA = 0,   otherwise. */
/*             If FACT = 'N', then */
/*                LDWORK  = MAX(1, 5*N, 2*N*N),        if JOB = 'C'; */
/*                LDWORK  = MAX(1, LWA + 5*N, 4*N*N ), if JOB = 'E', 'B'. */
/*             If FACT = 'F', then */
/*                LDWORK  = MAX(1, 2*N*N),  if JOB = 'C'; */
/*                LDWORK  = MAX(1, 4*N*N ), if JOB = 'E' or 'B'. */
/*             For good performance, LDWORK must generally be larger. */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, i <= N, the QR algorithm failed to */
/*                   complete the reduction of the matrix Ac to Schur */
/*                   canonical form (see LAPACK Library routine DGEES); */
/*                   on exit, the matrix T(i+1:N,i+1:N) contains the */
/*                   partially converged Schur form, and DWORK(i+1:N) and */
/*                   DWORK(N+i+1:2*N) contain the real and imaginary */
/*                   parts, respectively, of the converged eigenvalues; */
/*                   this error is unlikely to appear; */
/*             = N+1:  if the matrices T and -T' have common or very */
/*                   close eigenvalues; perturbed values were used to */
/*                   solve Lyapunov equations, but the matrix T, if given */
/*                   (for FACT = 'F'), is unchanged. */

/*     METHOD */

/*     The condition number of the Riccati equation is estimated as */

/*     cond = ( norm(Theta)*norm(A) + norm(inv(Omega))*norm(Q) + */
/*                 norm(Pi)*norm(G) ) / norm(X), */

/*     where Omega, Theta and Pi are linear operators defined by */

/*     Omega(W) = op(Ac)'*W + W*op(Ac), */
/*     Theta(W) = inv(Omega(op(W)'*X + X*op(W))), */
/*        Pi(W) = inv(Omega(X*W*X)), */

/*     and Ac = A - G*X (if TRANA = 'N') or Ac = A - X*G (if TRANA = 'T' */
/*     or 'C'). Note that the Riccati equation (1) is equivalent to */
/*                _   _         _   _ _ _ */
/*         op(T)'*X + X*op(T) + Q + X*G*X = 0,                        (2) */
/*           _           _               _ */
/*     where X = U'*X*U, Q = U'*Q*U, and G = U'*G*U, with U the */
/*     orthogonal matrix reducing Ac to a real Schur form, T = U'*Ac*U. */

/*     The routine estimates the quantities */

/*     sep(op(Ac),-op(Ac)') = 1 / norm(inv(Omega)), */

/*     norm(Theta) and norm(Pi) using 1-norm condition estimator. */

/*     The forward error bound is estimated using a practical error bound */
/*     similar to the one proposed in [2]. */

/*     REFERENCES */

/*     [1] Ghavimi, A.R. and Laub, A.J. */
/*         Backward error, sensitivity, and refinement of computed */
/*         solutions of algebraic Riccati equations. */
/*         Numerical Linear Algebra with Applications, vol. 2, pp. 29-49, */
/*         1995. */

/*     [2] Higham, N.J. */
/*         Perturbation theory and backward error for AX-XB=C. */
/*         BIT, vol. 33, pp. 124-136, 1993. */

/*     [3] Petkov, P.Hr., Konstantinov, M.M., and Mehrmann, V. */
/*         DGRSVX and DMSRIC: Fortran 77 subroutines for solving */
/*         continuous-time matrix algebraic Riccati equations with */
/*         condition and accuracy estimates. */
/*         Preprint SFB393/98-16, Fak. f. Mathematik, Tech. Univ. */
/*         Chemnitz, May 1998. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */
/*     The accuracy of the estimates obtained depends on the solution */
/*     accuracy and on the properties of the 1-norm estimator. */

/*     FURTHER COMMENTS */

/*     The option LYAPUN = 'R' may occasionally produce slightly worse */
/*     or better estimates, and it is much faster than the option 'O'. */
/*     When SEP is computed and it is zero, the routine returns */
/*     immediately, with RCOND and FERR (if requested) set to 0 and 1, */
/*     respectively. In this case, the equation is singular. */

/*     CONTRIBUTOR */

/*     P.Hr. Petkov, Technical University of Sofia, December 1998. */
/*     V. Sima, Katholieke Univ. Leuven, Belgium, February 1999. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2004. */

/*     KEYWORDS */

/*     Conditioning, error estimates, orthogonal transformation, */
/*     real Schur form, Riccati equation. */

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

#line 349 "SB02QD.f"
    /* Parameter adjustments */
#line 349 "SB02QD.f"
    a_dim1 = *lda;
#line 349 "SB02QD.f"
    a_offset = 1 + a_dim1;
#line 349 "SB02QD.f"
    a -= a_offset;
#line 349 "SB02QD.f"
    t_dim1 = *ldt;
#line 349 "SB02QD.f"
    t_offset = 1 + t_dim1;
#line 349 "SB02QD.f"
    t -= t_offset;
#line 349 "SB02QD.f"
    u_dim1 = *ldu;
#line 349 "SB02QD.f"
    u_offset = 1 + u_dim1;
#line 349 "SB02QD.f"
    u -= u_offset;
#line 349 "SB02QD.f"
    g_dim1 = *ldg;
#line 349 "SB02QD.f"
    g_offset = 1 + g_dim1;
#line 349 "SB02QD.f"
    g -= g_offset;
#line 349 "SB02QD.f"
    q_dim1 = *ldq;
#line 349 "SB02QD.f"
    q_offset = 1 + q_dim1;
#line 349 "SB02QD.f"
    q -= q_offset;
#line 349 "SB02QD.f"
    x_dim1 = *ldx;
#line 349 "SB02QD.f"
    x_offset = 1 + x_dim1;
#line 349 "SB02QD.f"
    x -= x_offset;
#line 349 "SB02QD.f"
    --iwork;
#line 349 "SB02QD.f"
    --dwork;
#line 349 "SB02QD.f"

#line 349 "SB02QD.f"
    /* Function Body */
#line 349 "SB02QD.f"
    jobc = lsame_(job, "C", (ftnlen)1, (ftnlen)1);
#line 350 "SB02QD.f"
    jobe = lsame_(job, "E", (ftnlen)1, (ftnlen)1);
#line 351 "SB02QD.f"
    jobb = lsame_(job, "B", (ftnlen)1, (ftnlen)1);
#line 352 "SB02QD.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 353 "SB02QD.f"
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
#line 354 "SB02QD.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 355 "SB02QD.f"
    update = lsame_(lyapun, "O", (ftnlen)1, (ftnlen)1);

#line 357 "SB02QD.f"
    needac = update && ! jobc;

#line 359 "SB02QD.f"
    nn = *n * *n;
#line 360 "SB02QD.f"
    if (needac) {
#line 361 "SB02QD.f"
	lwa = nn;
#line 362 "SB02QD.f"
    } else {
#line 363 "SB02QD.f"
	lwa = 0;
#line 364 "SB02QD.f"
    }

#line 366 "SB02QD.f"
    if (nofact) {
#line 367 "SB02QD.f"
	if (jobc) {
/* Computing MAX */
#line 368 "SB02QD.f"
	    i__1 = *n * 5, i__2 = nn << 1;
#line 368 "SB02QD.f"
	    ldw = max(i__1,i__2);
#line 369 "SB02QD.f"
	} else {
/* Computing MAX */
#line 370 "SB02QD.f"
	    i__1 = lwa + *n * 5, i__2 = nn << 2;
#line 370 "SB02QD.f"
	    ldw = max(i__1,i__2);
#line 371 "SB02QD.f"
	}
#line 372 "SB02QD.f"
    } else {
#line 373 "SB02QD.f"
	if (jobc) {
#line 374 "SB02QD.f"
	    ldw = nn << 1;
#line 375 "SB02QD.f"
	} else {
#line 376 "SB02QD.f"
	    ldw = nn << 2;
#line 377 "SB02QD.f"
	}
#line 378 "SB02QD.f"
    }

#line 380 "SB02QD.f"
    *info = 0;
#line 381 "SB02QD.f"
    if (! (jobb || jobc || jobe)) {
#line 382 "SB02QD.f"
	*info = -1;
#line 383 "SB02QD.f"
    } else if (! (nofact || lsame_(fact, "F", (ftnlen)1, (ftnlen)1))) {
#line 384 "SB02QD.f"
	*info = -2;
#line 385 "SB02QD.f"
    } else if (! (notrna || lsame_(trana, "T", (ftnlen)1, (ftnlen)1) || 
	    lsame_(trana, "C", (ftnlen)1, (ftnlen)1))) {
#line 387 "SB02QD.f"
	*info = -3;
#line 388 "SB02QD.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 389 "SB02QD.f"
	*info = -4;
#line 390 "SB02QD.f"
    } else if (! (update || lsame_(lyapun, "R", (ftnlen)1, (ftnlen)1))) {
#line 391 "SB02QD.f"
	*info = -5;
#line 392 "SB02QD.f"
    } else if (*n < 0) {
#line 393 "SB02QD.f"
	*info = -6;
#line 394 "SB02QD.f"
    } else if (*lda < 1 || *lda < *n && (update || nofact)) {
#line 396 "SB02QD.f"
	*info = -8;
#line 397 "SB02QD.f"
    } else if (*ldt < max(1,*n)) {
#line 398 "SB02QD.f"
	*info = -10;
#line 399 "SB02QD.f"
    } else if (*ldu < 1 || *ldu < *n && update) {
#line 400 "SB02QD.f"
	*info = -12;
#line 401 "SB02QD.f"
    } else if (*ldg < max(1,*n)) {
#line 402 "SB02QD.f"
	*info = -14;
#line 403 "SB02QD.f"
    } else if (*ldq < max(1,*n)) {
#line 404 "SB02QD.f"
	*info = -16;
#line 405 "SB02QD.f"
    } else if (*ldx < max(1,*n)) {
#line 406 "SB02QD.f"
	*info = -18;
#line 407 "SB02QD.f"
    } else if (*ldwork < max(1,ldw)) {
#line 408 "SB02QD.f"
	*info = -24;
#line 409 "SB02QD.f"
    }

#line 411 "SB02QD.f"
    if (*info != 0) {
#line 412 "SB02QD.f"
	i__1 = -(*info);
#line 412 "SB02QD.f"
	xerbla_("SB02QD", &i__1, (ftnlen)6);
#line 413 "SB02QD.f"
	return 0;
#line 414 "SB02QD.f"
    }

/*     Quick return if possible. */

#line 418 "SB02QD.f"
    if (*n == 0) {
#line 419 "SB02QD.f"
	if (! jobe) {
#line 419 "SB02QD.f"
	    *rcond = 1.;
#line 419 "SB02QD.f"
	}
#line 421 "SB02QD.f"
	if (! jobc) {
#line 421 "SB02QD.f"
	    *ferr = 0.;
#line 421 "SB02QD.f"
	}
#line 423 "SB02QD.f"
	dwork[1] = 1.;
#line 424 "SB02QD.f"
	return 0;
#line 425 "SB02QD.f"
    }

/*     Compute the 1-norm of the matrix X. */

#line 429 "SB02QD.f"
    xnorm = dlansy_("1-norm", uplo, n, &x[x_offset], ldx, &dwork[1], (ftnlen)
	    6, (ftnlen)1);
#line 430 "SB02QD.f"
    if (xnorm == 0.) {

/*        The solution is zero. */

#line 434 "SB02QD.f"
	if (! jobe) {
#line 434 "SB02QD.f"
	    *rcond = 0.;
#line 434 "SB02QD.f"
	}
#line 436 "SB02QD.f"
	if (! jobc) {
#line 436 "SB02QD.f"
	    *ferr = 0.;
#line 436 "SB02QD.f"
	}
#line 438 "SB02QD.f"
	dwork[1] = (doublereal) (*n);
#line 439 "SB02QD.f"
	return 0;
#line 440 "SB02QD.f"
    }

/*     Workspace usage. */

#line 444 "SB02QD.f"
    ixbs = 0;
#line 445 "SB02QD.f"
    itmp = ixbs + nn;
#line 446 "SB02QD.f"
    iabs = itmp + nn;
#line 447 "SB02QD.f"
    ires = iabs + nn;

/*     Workspace:  LWR, where */
/*                 LWR = N*N, if LYAPUN = 'O' and JOB = 'E' or 'B', or */
/*                               FACT = 'N', */
/*                 LWR = 0,   otherwise. */

#line 454 "SB02QD.f"
    if (needac || nofact) {

#line 456 "SB02QD.f"
	dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[1], n, (ftnlen)4);
#line 457 "SB02QD.f"
	if (notrna) {

/*           Compute Ac = A - G*X. */

#line 461 "SB02QD.f"
	    dsymm_("Left", uplo, n, n, &c_b18, &g[g_offset], ldg, &x[x_offset]
		    , ldx, &c_b19, &dwork[1], n, (ftnlen)4, (ftnlen)1);
#line 463 "SB02QD.f"
	} else {

/*           Compute Ac = A - X*G. */

#line 467 "SB02QD.f"
	    dsymm_("Right", uplo, n, n, &c_b18, &g[g_offset], ldg, &x[
		    x_offset], ldx, &c_b19, &dwork[1], n, (ftnlen)5, (ftnlen)
		    1);
#line 469 "SB02QD.f"
	}

#line 471 "SB02QD.f"
	wrkopt = (integer) ((doublereal) nn);
#line 472 "SB02QD.f"
	if (nofact) {
#line 472 "SB02QD.f"
	    dlacpy_("Full", n, n, &dwork[1], n, &t[t_offset], ldt, (ftnlen)4);
#line 472 "SB02QD.f"
	}
#line 474 "SB02QD.f"
    } else {
#line 475 "SB02QD.f"
	wrkopt = (integer) ((doublereal) (*n));
#line 476 "SB02QD.f"
    }

#line 478 "SB02QD.f"
    if (nofact) {

/*        Compute the Schur factorization of Ac, Ac = U*T*U'. */
/*        Workspace:  need   LWA + 5*N; */
/*                    prefer larger; */
/*                    LWA = N*N, if LYAPUN = 'O' and JOB = 'E' or 'B'; */
/*                    LWA = 0,   otherwise. */
/*        (Note: Comments in the code beginning "Workspace:" describe the */
/*        minimal amount of real workspace needed at that point in the */
/*        code, as well as the preferred amount for good performance.) */

#line 489 "SB02QD.f"
	if (update) {
#line 490 "SB02QD.f"
	    *(unsigned char *)sjob = 'V';
#line 491 "SB02QD.f"
	} else {
#line 492 "SB02QD.f"
	    *(unsigned char *)sjob = 'N';
#line 493 "SB02QD.f"
	}
#line 494 "SB02QD.f"
	i__1 = *ldwork - lwa - (*n << 1);
#line 494 "SB02QD.f"
	dgees_(sjob, "Not ordered", (L_fp)select_, n, &t[t_offset], ldt, &
		sdim, &dwork[lwa + 1], &dwork[lwa + *n + 1], &u[u_offset], 
		ldu, &dwork[lwa + (*n << 1) + 1], &i__1, bwork, info, (ftnlen)
		1, (ftnlen)11);
#line 497 "SB02QD.f"
	if (*info > 0) {
#line 498 "SB02QD.f"
	    if (lwa > 0) {
#line 498 "SB02QD.f"
		i__1 = *n << 1;
#line 498 "SB02QD.f"
		dcopy_(&i__1, &dwork[lwa + 1], &c__1, &dwork[1], &c__1);
#line 498 "SB02QD.f"
	    }
#line 500 "SB02QD.f"
	    return 0;
#line 501 "SB02QD.f"
	}

/* Computing MAX */
#line 503 "SB02QD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[lwa + (*n << 1) + 1] + lwa + (*
		n << 1);
#line 503 "SB02QD.f"
	wrkopt = max(i__1,i__2);
#line 504 "SB02QD.f"
    }
#line 505 "SB02QD.f"
    if (needac) {
#line 505 "SB02QD.f"
	dlacpy_("Full", n, n, &dwork[1], n, &dwork[iabs + 1], n, (ftnlen)4);
#line 505 "SB02QD.f"
    }

#line 508 "SB02QD.f"
    if (notrna) {
#line 509 "SB02QD.f"
	*(unsigned char *)tranat = 'T';
#line 510 "SB02QD.f"
    } else {
#line 511 "SB02QD.f"
	*(unsigned char *)tranat = 'N';
#line 512 "SB02QD.f"
    }

#line 514 "SB02QD.f"
    if (! jobe) {

/*        Estimate sep(op(Ac),-op(Ac)') = sep(op(T),-op(T)') and */
/*        norm(Theta). */
/*        Workspace LWA + 2*N*N. */

#line 520 "SB02QD.f"
	sb03qy_("Both", trana, lyapun, n, &t[t_offset], ldt, &u[u_offset], 
		ldu, &x[x_offset], ldx, sep, &thnorm, &iwork[1], &dwork[1], 
		ldwork, info, (ftnlen)4, (ftnlen)1, (ftnlen)1);

/* Computing MAX */
#line 523 "SB02QD.f"
	i__1 = wrkopt, i__2 = lwa + (nn << 1);
#line 523 "SB02QD.f"
	wrkopt = max(i__1,i__2);

/*        Return if the equation is singular. */

#line 527 "SB02QD.f"
	if (*sep == 0.) {
#line 528 "SB02QD.f"
	    *rcond = 0.;
#line 529 "SB02QD.f"
	    if (jobb) {
#line 529 "SB02QD.f"
		*ferr = 1.;
#line 529 "SB02QD.f"
	    }
#line 531 "SB02QD.f"
	    dwork[1] = (doublereal) wrkopt;
#line 532 "SB02QD.f"
	    return 0;
#line 533 "SB02QD.f"
	}

/*        Estimate norm(Pi). */
/*        Workspace LWA + 2*N*N. */

#line 538 "SB02QD.f"
	kase = 0;

/*        REPEAT */
#line 541 "SB02QD.f"
L10:
#line 542 "SB02QD.f"
	dlacon_(&nn, &dwork[itmp + 1], &dwork[1], &iwork[1], &est, &kase);
#line 543 "SB02QD.f"
	if (kase != 0) {

/*           Select the triangular part of symmetric matrix to be used. */

#line 547 "SB02QD.f"
	    if (dlansy_("1-norm", "Upper", n, &dwork[1], n, &dwork[itmp + 1], 
		    (ftnlen)6, (ftnlen)5) >= dlansy_("1-norm", "Lower", n, &
		    dwork[1], n, &dwork[itmp + 1], (ftnlen)6, (ftnlen)5)) {
#line 551 "SB02QD.f"
		*(unsigned char *)loup = 'U';
#line 552 "SB02QD.f"
	    } else {
#line 553 "SB02QD.f"
		*(unsigned char *)loup = 'L';
#line 554 "SB02QD.f"
	    }

/*           Compute RHS = X*W*X. */

#line 558 "SB02QD.f"
	    mb01ru_(loup, "No Transpose", n, n, &c_b41, &c_b19, &dwork[1], n, 
		    &x[x_offset], ldx, &dwork[1], n, &dwork[itmp + 1], &nn, &
		    info2, (ftnlen)1, (ftnlen)12);
#line 561 "SB02QD.f"
	    i__1 = *n + 1;
#line 561 "SB02QD.f"
	    dscal_(n, &c_b43, &dwork[1], &i__1);

#line 563 "SB02QD.f"
	    if (update) {

/*              Transform the right-hand side: RHS := U'*RHS*U. */

#line 567 "SB02QD.f"
		mb01ru_(loup, "Transpose", n, n, &c_b41, &c_b19, &dwork[1], n,
			 &u[u_offset], ldu, &dwork[1], n, &dwork[itmp + 1], &
			nn, &info2, (ftnlen)1, (ftnlen)9);
#line 570 "SB02QD.f"
		i__1 = *n + 1;
#line 570 "SB02QD.f"
		dscal_(n, &c_b43, &dwork[1], &i__1);
#line 571 "SB02QD.f"
	    }

/*           Fill in the remaining triangle of the symmetric matrix. */

#line 575 "SB02QD.f"
	    ma02ed_(loup, n, &dwork[1], n, (ftnlen)1);

#line 577 "SB02QD.f"
	    if (kase == 1) {

/*              Solve op(T)'*Y + Y*op(T) = scale*RHS. */

#line 581 "SB02QD.f"
		sb03my_(trana, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
			info2, (ftnlen)1);
#line 582 "SB02QD.f"
	    } else {

/*              Solve op(T)*W + W*op(T)' = scale*RHS. */

#line 586 "SB02QD.f"
		sb03my_(tranat, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
			info2, (ftnlen)1);
#line 587 "SB02QD.f"
	    }

#line 589 "SB02QD.f"
	    if (update) {

/*              Transform back to obtain the solution: Z := U*Z*U', with */
/*              Z = Y or Z = W. */

#line 594 "SB02QD.f"
		mb01ru_(loup, "No transpose", n, n, &c_b41, &c_b19, &dwork[1],
			 n, &u[u_offset], ldu, &dwork[1], n, &dwork[itmp + 1],
			 &nn, &info2, (ftnlen)1, (ftnlen)12);
#line 597 "SB02QD.f"
		i__1 = *n + 1;
#line 597 "SB02QD.f"
		dscal_(n, &c_b43, &dwork[1], &i__1);

/*              Fill in the remaining triangle of the symmetric matrix. */

#line 601 "SB02QD.f"
		ma02ed_(loup, n, &dwork[1], n, (ftnlen)1);
#line 602 "SB02QD.f"
	    }
#line 603 "SB02QD.f"
	    goto L10;
#line 604 "SB02QD.f"
	}
/*        UNTIL KASE = 0 */

#line 607 "SB02QD.f"
	if (est < scale) {
#line 608 "SB02QD.f"
	    pinorm = est / scale;
#line 609 "SB02QD.f"
	} else {
#line 610 "SB02QD.f"
	    bignum = 1. / dlamch_("Safe minimum", (ftnlen)12);
#line 611 "SB02QD.f"
	    if (est < scale * bignum) {
#line 612 "SB02QD.f"
		pinorm = est / scale;
#line 613 "SB02QD.f"
	    } else {
#line 614 "SB02QD.f"
		pinorm = bignum;
#line 615 "SB02QD.f"
	    }
#line 616 "SB02QD.f"
	}

/*        Compute the 1-norm of A or T. */

#line 620 "SB02QD.f"
	if (update) {
#line 621 "SB02QD.f"
	    anorm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[1], (
		    ftnlen)6);
#line 622 "SB02QD.f"
	} else {
#line 623 "SB02QD.f"
	    anorm = dlanhs_("1-norm", n, &t[t_offset], ldt, &dwork[1], (
		    ftnlen)6);
#line 624 "SB02QD.f"
	}

/*        Compute the 1-norms of the matrices Q and G. */

#line 628 "SB02QD.f"
	qnorm = dlansy_("1-norm", uplo, n, &q[q_offset], ldq, &dwork[1], (
		ftnlen)6, (ftnlen)1);
#line 629 "SB02QD.f"
	gnorm = dlansy_("1-norm", uplo, n, &g[g_offset], ldg, &dwork[1], (
		ftnlen)6, (ftnlen)1);

/*        Estimate the reciprocal condition number. */

/* Computing MAX */
#line 633 "SB02QD.f"
	d__1 = max(*sep,xnorm), d__1 = max(d__1,anorm);
#line 633 "SB02QD.f"
	tmax = max(d__1,gnorm);
#line 634 "SB02QD.f"
	if (tmax <= 1.) {
#line 635 "SB02QD.f"
	    temp = *sep * xnorm;
#line 636 "SB02QD.f"
	    denom = qnorm + *sep * anorm * thnorm + *sep * gnorm * pinorm;
#line 638 "SB02QD.f"
	} else {
#line 639 "SB02QD.f"
	    temp = *sep / tmax * (xnorm / tmax);
#line 640 "SB02QD.f"
	    denom = 1. / tmax * (qnorm / tmax) + *sep / tmax * (anorm / tmax) 
		    * thnorm + *sep / tmax * (gnorm / tmax) * pinorm;
#line 643 "SB02QD.f"
	}
#line 644 "SB02QD.f"
	if (temp >= denom) {
#line 645 "SB02QD.f"
	    *rcond = 1.;
#line 646 "SB02QD.f"
	} else {
#line 647 "SB02QD.f"
	    *rcond = temp / denom;
#line 648 "SB02QD.f"
	}
#line 649 "SB02QD.f"
    }

#line 651 "SB02QD.f"
    if (! jobc) {

/*        Form a triangle of the residual matrix */
/*          R = op(A)'*X + X*op(A) + Q - X*G*X, */
/*        or           _   _         _   _ _ _ */
/*          R = op(T)'*X + X*op(T) + Q + X*G*X, */
/*        exploiting the symmetry. */
/*        Workspace 4*N*N. */

#line 660 "SB02QD.f"
	if (update) {
#line 661 "SB02QD.f"
	    dlacpy_(uplo, n, n, &q[q_offset], ldq, &dwork[ires + 1], n, (
		    ftnlen)1);
#line 662 "SB02QD.f"
	    dsyr2k_(uplo, tranat, n, n, &c_b19, &a[a_offset], lda, &x[
		    x_offset], ldx, &c_b19, &dwork[ires + 1], n, (ftnlen)1, (
		    ftnlen)1);
#line 664 "SB02QD.f"
	    sig = -1.;
#line 665 "SB02QD.f"
	} else {
#line 666 "SB02QD.f"
	    mb01ud_("Right", trana, n, n, &c_b19, &t[t_offset], ldt, &x[
		    x_offset], ldx, &dwork[ires + 1], n, &info2, (ftnlen)5, (
		    ftnlen)1);
#line 668 "SB02QD.f"
	    jj = ires + 1;
#line 669 "SB02QD.f"
	    if (lower) {
#line 670 "SB02QD.f"
		i__1 = *n;
#line 670 "SB02QD.f"
		for (j = 1; j <= i__1; ++j) {
#line 671 "SB02QD.f"
		    i__2 = *n - j + 1;
#line 671 "SB02QD.f"
		    daxpy_(&i__2, &c_b19, &dwork[jj], n, &dwork[jj], &c__1);
#line 673 "SB02QD.f"
		    i__2 = *n - j + 1;
#line 673 "SB02QD.f"
		    daxpy_(&i__2, &c_b19, &q[j + j * q_dim1], &c__1, &dwork[
			    jj], &c__1);
#line 674 "SB02QD.f"
		    jj = jj + *n + 1;
#line 675 "SB02QD.f"
/* L20: */
#line 675 "SB02QD.f"
		}
#line 676 "SB02QD.f"
	    } else {
#line 677 "SB02QD.f"
		i__1 = *n;
#line 677 "SB02QD.f"
		for (j = 1; j <= i__1; ++j) {
#line 678 "SB02QD.f"
		    daxpy_(&j, &c_b19, &dwork[ires + j], n, &dwork[jj], &c__1)
			    ;
#line 680 "SB02QD.f"
		    daxpy_(&j, &c_b19, &q[j * q_dim1 + 1], &c__1, &dwork[jj], 
			    &c__1);
#line 681 "SB02QD.f"
		    jj += *n;
#line 682 "SB02QD.f"
/* L30: */
#line 682 "SB02QD.f"
		}
#line 683 "SB02QD.f"
	    }
#line 684 "SB02QD.f"
	    sig = 1.;
#line 685 "SB02QD.f"
	}
#line 686 "SB02QD.f"
	mb01ru_(uplo, tranat, n, n, &c_b19, &sig, &dwork[ires + 1], n, &x[
		x_offset], ldx, &g[g_offset], ldg, &dwork[itmp + 1], &nn, &
		info2, (ftnlen)1, (ftnlen)1);

/*        Get the machine precision. */

#line 691 "SB02QD.f"
	eps = dlamch_("Epsilon", (ftnlen)7);
#line 692 "SB02QD.f"
	epsn = eps * (doublereal) (*n + 4);
#line 693 "SB02QD.f"
	temp = eps * 4.;

/*        Add to abs(R) a term that takes account of rounding errors in */
/*        forming R: */
/*         abs(R) := abs(R) + EPS*(4*abs(Q) + (n+4)*(abs(op(Ac))'*abs(X) */
/*                 + abs(X)*abs(op(Ac))) + 2*(n+1)*abs(X)*abs(G)*abs(X)), */
/*        or                             _                           _ */
/*         abs(R) := abs(R) + EPS*(4*abs(Q) + (n+4)*(abs(op(T))'*abs(X) */
/*                       _                            _      _      _ */
/*                 + abs(X)*abs(op(T))) + 2*(n+1)*abs(X)*abs(G)*abs(X)), */
/*        where EPS is the machine precision. */

#line 705 "SB02QD.f"
	i__1 = *n;
#line 705 "SB02QD.f"
	for (j = 1; j <= i__1; ++j) {
#line 706 "SB02QD.f"
	    i__2 = *n;
#line 706 "SB02QD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 707 "SB02QD.f"
		dwork[ixbs + (j - 1) * *n + i__] = (d__1 = x[i__ + j * x_dim1]
			, abs(d__1));
#line 708 "SB02QD.f"
/* L40: */
#line 708 "SB02QD.f"
	    }
#line 709 "SB02QD.f"
/* L50: */
#line 709 "SB02QD.f"
	}

#line 711 "SB02QD.f"
	if (lower) {
#line 712 "SB02QD.f"
	    i__1 = *n;
#line 712 "SB02QD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 713 "SB02QD.f"
		i__2 = *n;
#line 713 "SB02QD.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 714 "SB02QD.f"
		    dwork[ires + (j - 1) * *n + i__] = temp * (d__1 = q[i__ + 
			    j * q_dim1], abs(d__1)) + (d__2 = dwork[ires + (j 
			    - 1) * *n + i__], abs(d__2));
#line 716 "SB02QD.f"
/* L60: */
#line 716 "SB02QD.f"
		}
#line 717 "SB02QD.f"
/* L70: */
#line 717 "SB02QD.f"
	    }
#line 718 "SB02QD.f"
	} else {
#line 719 "SB02QD.f"
	    i__1 = *n;
#line 719 "SB02QD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 720 "SB02QD.f"
		i__2 = j;
#line 720 "SB02QD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 721 "SB02QD.f"
		    dwork[ires + (j - 1) * *n + i__] = temp * (d__1 = q[i__ + 
			    j * q_dim1], abs(d__1)) + (d__2 = dwork[ires + (j 
			    - 1) * *n + i__], abs(d__2));
#line 723 "SB02QD.f"
/* L80: */
#line 723 "SB02QD.f"
		}
#line 724 "SB02QD.f"
/* L90: */
#line 724 "SB02QD.f"
	    }
#line 725 "SB02QD.f"
	}

#line 727 "SB02QD.f"
	if (update) {

#line 729 "SB02QD.f"
	    i__1 = *n;
#line 729 "SB02QD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 730 "SB02QD.f"
		i__2 = *n;
#line 730 "SB02QD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 731 "SB02QD.f"
		    dwork[iabs + (j - 1) * *n + i__] = (d__1 = dwork[iabs + (
			    j - 1) * *n + i__], abs(d__1));
#line 733 "SB02QD.f"
/* L100: */
#line 733 "SB02QD.f"
		}
#line 734 "SB02QD.f"
/* L110: */
#line 734 "SB02QD.f"
	    }

#line 736 "SB02QD.f"
	    dsyr2k_(uplo, tranat, n, n, &epsn, &dwork[iabs + 1], n, &dwork[
		    ixbs + 1], n, &c_b19, &dwork[ires + 1], n, (ftnlen)1, (
		    ftnlen)1);
#line 738 "SB02QD.f"
	} else {

#line 740 "SB02QD.f"
	    i__1 = *n;
#line 740 "SB02QD.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 741 "SB02QD.f"
		i__3 = j + 1;
#line 741 "SB02QD.f"
		i__2 = min(i__3,*n);
#line 741 "SB02QD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 742 "SB02QD.f"
		    dwork[iabs + (j - 1) * *n + i__] = (d__1 = t[i__ + j * 
			    t_dim1], abs(d__1));
#line 743 "SB02QD.f"
/* L120: */
#line 743 "SB02QD.f"
		}
#line 744 "SB02QD.f"
/* L130: */
#line 744 "SB02QD.f"
	    }

#line 746 "SB02QD.f"
	    mb01ud_("Left", tranat, n, n, &epsn, &dwork[iabs + 1], n, &dwork[
		    ixbs + 1], n, &dwork[itmp + 1], n, &info2, (ftnlen)4, (
		    ftnlen)1);
#line 748 "SB02QD.f"
	    jj = ires + 1;
#line 749 "SB02QD.f"
	    jx = itmp + 1;
#line 750 "SB02QD.f"
	    if (lower) {
#line 751 "SB02QD.f"
		i__1 = *n;
#line 751 "SB02QD.f"
		for (j = 1; j <= i__1; ++j) {
#line 752 "SB02QD.f"
		    i__2 = *n - j + 1;
#line 752 "SB02QD.f"
		    daxpy_(&i__2, &c_b19, &dwork[jx], n, &dwork[jx], &c__1);
#line 754 "SB02QD.f"
		    i__2 = *n - j + 1;
#line 754 "SB02QD.f"
		    daxpy_(&i__2, &c_b19, &dwork[jx], &c__1, &dwork[jj], &
			    c__1);
#line 756 "SB02QD.f"
		    jj = jj + *n + 1;
#line 757 "SB02QD.f"
		    jx = jx + *n + 1;
#line 758 "SB02QD.f"
/* L140: */
#line 758 "SB02QD.f"
		}
#line 759 "SB02QD.f"
	    } else {
#line 760 "SB02QD.f"
		i__1 = *n;
#line 760 "SB02QD.f"
		for (j = 1; j <= i__1; ++j) {
#line 761 "SB02QD.f"
		    daxpy_(&j, &c_b19, &dwork[itmp + j], n, &dwork[jx], &c__1)
			    ;
#line 763 "SB02QD.f"
		    daxpy_(&j, &c_b19, &dwork[jx], &c__1, &dwork[jj], &c__1);
#line 764 "SB02QD.f"
		    jj += *n;
#line 765 "SB02QD.f"
		    jx += *n;
#line 766 "SB02QD.f"
/* L150: */
#line 766 "SB02QD.f"
		}
#line 767 "SB02QD.f"
	    }
#line 768 "SB02QD.f"
	}

#line 770 "SB02QD.f"
	if (lower) {
#line 771 "SB02QD.f"
	    i__1 = *n;
#line 771 "SB02QD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 772 "SB02QD.f"
		i__2 = *n;
#line 772 "SB02QD.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 773 "SB02QD.f"
		    dwork[iabs + (j - 1) * *n + i__] = (d__1 = g[i__ + j * 
			    g_dim1], abs(d__1));
#line 774 "SB02QD.f"
/* L160: */
#line 774 "SB02QD.f"
		}
#line 775 "SB02QD.f"
/* L170: */
#line 775 "SB02QD.f"
	    }
#line 776 "SB02QD.f"
	} else {
#line 777 "SB02QD.f"
	    i__1 = *n;
#line 777 "SB02QD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 778 "SB02QD.f"
		i__2 = j;
#line 778 "SB02QD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 779 "SB02QD.f"
		    dwork[iabs + (j - 1) * *n + i__] = (d__1 = g[i__ + j * 
			    g_dim1], abs(d__1));
#line 780 "SB02QD.f"
/* L180: */
#line 780 "SB02QD.f"
		}
#line 781 "SB02QD.f"
/* L190: */
#line 781 "SB02QD.f"
	    }
#line 782 "SB02QD.f"
	}

#line 784 "SB02QD.f"
	d__1 = eps * (doublereal) (*n + 1 << 1);
#line 784 "SB02QD.f"
	mb01ru_(uplo, trana, n, n, &c_b19, &d__1, &dwork[ires + 1], n, &dwork[
		ixbs + 1], n, &dwork[iabs + 1], n, &dwork[itmp + 1], &nn, &
		info2, (ftnlen)1, (ftnlen)1);

/* Computing MAX */
#line 788 "SB02QD.f"
	i__1 = wrkopt, i__2 = nn << 2;
#line 788 "SB02QD.f"
	wrkopt = max(i__1,i__2);

/*        Compute forward error bound, using matrix norm estimator. */
/*        Workspace 4*N*N. */

#line 793 "SB02QD.f"
	xanorm = dlansy_("Max", uplo, n, &x[x_offset], ldx, &dwork[1], (
		ftnlen)3, (ftnlen)1);

#line 795 "SB02QD.f"
	sb03qx_(trana, uplo, lyapun, n, &xanorm, &t[t_offset], ldt, &u[
		u_offset], ldu, &dwork[ires + 1], n, ferr, &iwork[1], &dwork[
		1], &ires, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 798 "SB02QD.f"
    }

#line 800 "SB02QD.f"
    dwork[1] = (doublereal) wrkopt;
#line 801 "SB02QD.f"
    return 0;

/* *** Last line of SB02QD *** */
} /* sb02qd_ */

