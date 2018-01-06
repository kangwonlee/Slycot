#line 1 "SB02SD.f"
/* SB02SD.f -- translated by f2c (version 20100827).
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

#line 1 "SB02SD.f"
/* Table of constant values */

static doublereal c_b17 = 0.;
static doublereal c_b18 = 1.;
static integer c__1 = 1;
static doublereal c_b51 = .5;
static doublereal c_b66 = -1.;

/* Subroutine */ int sb02sd_(char *job, char *fact, char *trana, char *uplo, 
	char *lyapun, integer *n, doublereal *a, integer *lda, doublereal *t, 
	integer *ldt, doublereal *u, integer *ldu, doublereal *g, integer *
	ldg, doublereal *q, integer *ldq, doublereal *x, integer *ldx, 
	doublereal *sepd, doublereal *rcond, doublereal *ferr, integer *iwork,
	 doublereal *dwork, integer *ldwork, integer *info, ftnlen job_len, 
	ftnlen fact_len, ftnlen trana_len, ftnlen uplo_len, ftnlen lyapun_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, g_dim1, g_offset, q_dim1, q_offset, t_dim1, 
	    t_offset, u_dim1, u_offset, x_dim1, x_offset, i__1, i__2, i__3, 
	    i__4;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, jj, nn, lwa, ldw;
    static doublereal eps, est;
    static integer lwr;
    static logical jobb, jobc, jobe;
    static integer iabs, kase;
    static char sjob[1];
    static integer ixma, sdim, ires, ixbs;
    static doublereal epsn, temp, tmax;
    static integer iwrk;
    static doublereal epst;
    static char loup[1];
    static integer info2;
    extern /* Subroutine */ int ma02ed_(char *, integer *, doublereal *, 
	    integer *, ftnlen), dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal scale;
    extern /* Subroutine */ int dgees_(char *, char *, L_fp, integer *, 
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
	    integer *, ftnlen, ftnlen), dgesv_(integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *), mb01rx_(char *, char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen), mb01ry_(char *, char *, char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int sb03mx_(char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen), dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static doublereal gnorm;
    extern /* Subroutine */ int sb03sx_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    static logical bwork[1];
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), sb03sy_(char *, char *, char 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    static logical lower;
    extern /* Subroutine */ int dsymm_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal qnorm, xnorm;
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
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
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
/*     solution of the real discrete-time matrix algebraic Riccati */
/*     equation (see FURTHER COMMENTS) */
/*                                 -1 */
/*         X = op(A)'*X*(I_n + G*X)  *op(A) + Q,                      (1) */

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
/*             the matrix Ac = inv(I_n + G*X)*A (if TRANA = 'N'), or */
/*             Ac = A*inv(I_n + X*G) (if TRANA = 'T' or 'C'), is supplied */
/*             on entry, as follows: */
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

/*     SEPD    (output) DOUBLE PRECISION */
/*             If JOB = 'C' or JOB = 'B', the estimated quantity */
/*             sepd(op(Ac),op(Ac)'). */
/*             If N = 0, or X = 0, or JOB = 'E', SEPD is not referenced. */

/*     RCOND   (output) DOUBLE PRECISION */
/*             If JOB = 'C' or JOB = 'B', an estimate of the reciprocal */
/*             condition number of the discrete-time Riccati equation. */
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
/*             Let LWA = N*N, if LYAPUN = 'O'; */
/*                 LWA = 0,   otherwise, */
/*             and LWN = N,   if LYAPUN = 'R' and JOB = 'E' or 'B'; */
/*                 LWN = 0,   otherwise. */
/*             If FACT = 'N', then */
/*                LDWORK  = MAX(LWA + 5*N, MAX(3,2*N*N) + N*N), */
/*                                                 if JOB = 'C'; */
/*                LDWORK  = MAX(LWA + 5*N, MAX(3,2*N*N) + 2*N*N + LWN), */
/*                                                 if JOB = 'E' or 'B'. */
/*             If FACT = 'F', then */
/*                LDWORK  = MAX(3,2*N*N) + N*N,    if JOB = 'C'; */
/*                LDWORK  = MAX(3,2*N*N) + 2*N*N + LWN, */
/*                                                 if JOB = 'E' or 'B'. */
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
/*             = N+1:  if T has almost reciprocal eigenvalues; perturbed */
/*                   values were used to solve Lyapunov equations, but */
/*                   the matrix T, if given (for FACT = 'F'), is */
/*                   unchanged. */

/*     METHOD */

/*     The condition number of the Riccati equation is estimated as */

/*     cond = ( norm(Theta)*norm(A) + norm(inv(Omega))*norm(Q) + */
/*                 norm(Pi)*norm(G) ) / norm(X), */

/*     where Omega, Theta and Pi are linear operators defined by */

/*     Omega(W) = op(Ac)'*W*op(Ac) - W, */
/*     Theta(W) = inv(Omega(op(W)'*X*op(Ac) + op(Ac)'X*op(W))), */
/*        Pi(W) = inv(Omega(op(Ac)'*X*W*X*op(Ac))), */

/*     and Ac = inv(I_n + G*X)*A (if TRANA = 'N'), or */
/*         Ac = A*inv(I_n + X*G) (if TRANA = 'T' or 'C'). */

/*     Note that the Riccati equation (1) is equivalent to */

/*         X = op(Ac)'*X*op(Ac) + op(Ac)'*X*G*X*op(Ac) + Q,           (2) */

/*     and to */
/*         _          _                _ _ _         _ */
/*         X = op(T)'*X*op(T) + op(T)'*X*G*X*op(T) + Q,               (3) */
/*           _           _               _ */
/*     where X = U'*X*U, Q = U'*Q*U, and G = U'*G*U, with U the */
/*     orthogonal matrix reducing Ac to a real Schur form, T = U'*Ac*U. */

/*     The routine estimates the quantities */

/*     sepd(op(Ac),op(Ac)') = 1 / norm(inv(Omega)), */

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
/*     When SEPD is computed and it is zero, the routine returns */
/*     immediately, with RCOND and FERR (if requested) set to 0 and 1, */
/*     respectively. In this case, the equation is singular. */

/*     Let B be an N-by-M matrix (if TRANA = 'N') or an M-by-N matrix */
/*     (if TRANA = 'T' or 'C'), let R be an M-by-M symmetric positive */
/*     definite matrix (R = R**T), and denote G = op(B)*inv(R)*op(B)'. */
/*     Then, the Riccati equation (1) is equivalent to the standard */
/*     discrete-time matrix algebraic Riccati equation */

/*         X = op(A)'*X*op(A) -                                       (4) */
/*                                                -1 */
/*             op(A)'*X*op(B)*(R + op(B)'*X*op(B))  *op(B)'*X*op(A) + Q. */

/*     By symmetry, the equation (1) is also equivalent to */
/*                               -1 */
/*         X = op(A)'*(I_n + X*G)  *X*op(A) + Q. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, and */
/*     P.Hr. Petkov, Technical University of Sofia, March 1999. */

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

#line 376 "SB02SD.f"
    /* Parameter adjustments */
#line 376 "SB02SD.f"
    a_dim1 = *lda;
#line 376 "SB02SD.f"
    a_offset = 1 + a_dim1;
#line 376 "SB02SD.f"
    a -= a_offset;
#line 376 "SB02SD.f"
    t_dim1 = *ldt;
#line 376 "SB02SD.f"
    t_offset = 1 + t_dim1;
#line 376 "SB02SD.f"
    t -= t_offset;
#line 376 "SB02SD.f"
    u_dim1 = *ldu;
#line 376 "SB02SD.f"
    u_offset = 1 + u_dim1;
#line 376 "SB02SD.f"
    u -= u_offset;
#line 376 "SB02SD.f"
    g_dim1 = *ldg;
#line 376 "SB02SD.f"
    g_offset = 1 + g_dim1;
#line 376 "SB02SD.f"
    g -= g_offset;
#line 376 "SB02SD.f"
    q_dim1 = *ldq;
#line 376 "SB02SD.f"
    q_offset = 1 + q_dim1;
#line 376 "SB02SD.f"
    q -= q_offset;
#line 376 "SB02SD.f"
    x_dim1 = *ldx;
#line 376 "SB02SD.f"
    x_offset = 1 + x_dim1;
#line 376 "SB02SD.f"
    x -= x_offset;
#line 376 "SB02SD.f"
    --iwork;
#line 376 "SB02SD.f"
    --dwork;
#line 376 "SB02SD.f"

#line 376 "SB02SD.f"
    /* Function Body */
#line 376 "SB02SD.f"
    jobc = lsame_(job, "C", (ftnlen)1, (ftnlen)1);
#line 377 "SB02SD.f"
    jobe = lsame_(job, "E", (ftnlen)1, (ftnlen)1);
#line 378 "SB02SD.f"
    jobb = lsame_(job, "B", (ftnlen)1, (ftnlen)1);
#line 379 "SB02SD.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 380 "SB02SD.f"
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
#line 381 "SB02SD.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 382 "SB02SD.f"
    update = lsame_(lyapun, "O", (ftnlen)1, (ftnlen)1);

#line 384 "SB02SD.f"
    needac = update && ! jobc;

#line 386 "SB02SD.f"
    nn = *n * *n;
#line 387 "SB02SD.f"
    if (update) {
#line 388 "SB02SD.f"
	lwa = nn;
#line 389 "SB02SD.f"
    } else {
#line 390 "SB02SD.f"
	lwa = 0;
#line 391 "SB02SD.f"
    }

#line 393 "SB02SD.f"
    if (jobc) {
/* Computing MAX */
#line 394 "SB02SD.f"
	i__1 = 3, i__2 = nn << 1;
#line 394 "SB02SD.f"
	ldw = max(i__1,i__2) + nn;
#line 395 "SB02SD.f"
    } else {
/* Computing MAX */
#line 396 "SB02SD.f"
	i__1 = 3, i__2 = nn << 1;
#line 396 "SB02SD.f"
	ldw = max(i__1,i__2) + (nn << 1);
#line 397 "SB02SD.f"
	if (! update) {
#line 397 "SB02SD.f"
	    ldw += *n;
#line 397 "SB02SD.f"
	}
#line 399 "SB02SD.f"
    }
#line 400 "SB02SD.f"
    if (nofact) {
/* Computing MAX */
#line 400 "SB02SD.f"
	i__1 = lwa + *n * 5;
#line 400 "SB02SD.f"
	ldw = max(i__1,ldw);
#line 400 "SB02SD.f"
    }

#line 403 "SB02SD.f"
    *info = 0;
#line 404 "SB02SD.f"
    if (! (jobb || jobc || jobe)) {
#line 405 "SB02SD.f"
	*info = -1;
#line 406 "SB02SD.f"
    } else if (! (nofact || lsame_(fact, "F", (ftnlen)1, (ftnlen)1))) {
#line 407 "SB02SD.f"
	*info = -2;
#line 408 "SB02SD.f"
    } else if (! (notrna || lsame_(trana, "T", (ftnlen)1, (ftnlen)1) || 
	    lsame_(trana, "C", (ftnlen)1, (ftnlen)1))) {
#line 410 "SB02SD.f"
	*info = -3;
#line 411 "SB02SD.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 412 "SB02SD.f"
	*info = -4;
#line 413 "SB02SD.f"
    } else if (! (update || lsame_(lyapun, "R", (ftnlen)1, (ftnlen)1))) {
#line 414 "SB02SD.f"
	*info = -5;
#line 415 "SB02SD.f"
    } else if (*n < 0) {
#line 416 "SB02SD.f"
	*info = -6;
#line 417 "SB02SD.f"
    } else if (*lda < 1 || *lda < *n && (update || nofact)) {
#line 419 "SB02SD.f"
	*info = -8;
#line 420 "SB02SD.f"
    } else if (*ldt < max(1,*n)) {
#line 421 "SB02SD.f"
	*info = -10;
#line 422 "SB02SD.f"
    } else if (*ldu < 1 || *ldu < *n && update) {
#line 423 "SB02SD.f"
	*info = -12;
#line 424 "SB02SD.f"
    } else if (*ldg < max(1,*n)) {
#line 425 "SB02SD.f"
	*info = -14;
#line 426 "SB02SD.f"
    } else if (*ldq < max(1,*n)) {
#line 427 "SB02SD.f"
	*info = -16;
#line 428 "SB02SD.f"
    } else if (*ldx < max(1,*n)) {
#line 429 "SB02SD.f"
	*info = -18;
#line 430 "SB02SD.f"
    } else if (*ldwork < ldw) {
#line 431 "SB02SD.f"
	*info = -24;
#line 432 "SB02SD.f"
    }

#line 434 "SB02SD.f"
    if (*info != 0) {
#line 435 "SB02SD.f"
	i__1 = -(*info);
#line 435 "SB02SD.f"
	xerbla_("SB02SD", &i__1, (ftnlen)6);
#line 436 "SB02SD.f"
	return 0;
#line 437 "SB02SD.f"
    }

/*     Quick return if possible. */

#line 441 "SB02SD.f"
    if (*n == 0) {
#line 442 "SB02SD.f"
	if (! jobe) {
#line 442 "SB02SD.f"
	    *rcond = 1.;
#line 442 "SB02SD.f"
	}
#line 444 "SB02SD.f"
	if (! jobc) {
#line 444 "SB02SD.f"
	    *ferr = 0.;
#line 444 "SB02SD.f"
	}
#line 446 "SB02SD.f"
	dwork[1] = 1.;
#line 447 "SB02SD.f"
	return 0;
#line 448 "SB02SD.f"
    }

/*     Compute the 1-norm of the matrix X. */

#line 452 "SB02SD.f"
    xnorm = dlansy_("1-norm", uplo, n, &x[x_offset], ldx, &dwork[1], (ftnlen)
	    6, (ftnlen)1);
#line 453 "SB02SD.f"
    if (xnorm == 0.) {

/*        The solution is zero. */

#line 457 "SB02SD.f"
	if (! jobe) {
#line 457 "SB02SD.f"
	    *rcond = 0.;
#line 457 "SB02SD.f"
	}
#line 459 "SB02SD.f"
	if (! jobc) {
#line 459 "SB02SD.f"
	    *ferr = 0.;
#line 459 "SB02SD.f"
	}
#line 461 "SB02SD.f"
	dwork[1] = (doublereal) (*n);
#line 462 "SB02SD.f"
	return 0;
#line 463 "SB02SD.f"
    }

/*     Workspace usage. */

#line 467 "SB02SD.f"
    ires = 0;
#line 468 "SB02SD.f"
    ixbs = ires + nn;
/* Computing MAX */
#line 469 "SB02SD.f"
    i__1 = 3, i__2 = nn << 1;
#line 469 "SB02SD.f"
    ixma = max(i__1,i__2);
#line 470 "SB02SD.f"
    iabs = ixma + nn;
#line 471 "SB02SD.f"
    iwrk = iabs + nn;

/*     Workspace:  LWK, where */
/*                 LWK = 2*N*N, if LYAPUN = 'O', or FACT = 'N', */
/*                 LWK = N,     otherwise. */

#line 477 "SB02SD.f"
    if (update || nofact) {

#line 479 "SB02SD.f"
	dlaset_("Full", n, n, &c_b17, &c_b18, &dwork[ixbs + 1], n, (ftnlen)4);
#line 480 "SB02SD.f"
	dsymm_("Left", uplo, n, n, &c_b18, &g[g_offset], ldg, &x[x_offset], 
		ldx, &c_b18, &dwork[ixbs + 1], n, (ftnlen)4, (ftnlen)1);
#line 482 "SB02SD.f"
	if (notrna) {
/*                                   -1 */
/*           Compute Ac = (I_n + G*X)  *A. */

#line 486 "SB02SD.f"
	    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[1], n, (ftnlen)4);
#line 487 "SB02SD.f"
	    dgesv_(n, n, &dwork[ixbs + 1], n, &iwork[1], &dwork[1], n, &info2)
		    ;
#line 489 "SB02SD.f"
	} else {
/*                                     -1 */
/*           Compute Ac = A*(I_n + X*G)  . */

#line 493 "SB02SD.f"
	    i__1 = *n;
#line 493 "SB02SD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 494 "SB02SD.f"
		dcopy_(n, &a[j * a_dim1 + 1], &c__1, &dwork[j], n);
#line 495 "SB02SD.f"
/* L10: */
#line 495 "SB02SD.f"
	    }
#line 496 "SB02SD.f"
	    dgesv_(n, n, &dwork[ixbs + 1], n, &iwork[1], &dwork[1], n, &info2)
		    ;
#line 498 "SB02SD.f"
	    i__1 = *n;
#line 498 "SB02SD.f"
	    for (j = 2; j <= i__1; ++j) {
#line 499 "SB02SD.f"
		i__2 = j - 1;
#line 499 "SB02SD.f"
		dswap_(&i__2, &dwork[(j - 1) * *n + 1], &c__1, &dwork[j], n);
#line 500 "SB02SD.f"
/* L20: */
#line 500 "SB02SD.f"
	    }
#line 501 "SB02SD.f"
	}

#line 503 "SB02SD.f"
	wrkopt = (integer) ((doublereal) (nn << 1));
#line 504 "SB02SD.f"
	if (nofact) {
#line 504 "SB02SD.f"
	    dlacpy_("Full", n, n, &dwork[1], n, &t[t_offset], ldt, (ftnlen)4);
#line 504 "SB02SD.f"
	}
#line 506 "SB02SD.f"
    } else {
#line 507 "SB02SD.f"
	wrkopt = (integer) ((doublereal) (*n));
#line 508 "SB02SD.f"
    }

#line 510 "SB02SD.f"
    if (nofact) {

/*        Compute the Schur factorization of Ac, Ac = U*T*U'. */
/*        Workspace:  need   LWA + 5*N; */
/*                    prefer larger; */
/*                    LWA = N*N, if LYAPUN = 'O'; */
/*                    LWA = 0,   otherwise. */
/*        (Note: Comments in the code beginning "Workspace:" describe the */
/*        minimal amount of real workspace needed at that point in the */
/*        code, as well as the preferred amount for good performance.) */

#line 521 "SB02SD.f"
	if (update) {
#line 522 "SB02SD.f"
	    *(unsigned char *)sjob = 'V';
#line 523 "SB02SD.f"
	} else {
#line 524 "SB02SD.f"
	    *(unsigned char *)sjob = 'N';
#line 525 "SB02SD.f"
	}
#line 526 "SB02SD.f"
	i__1 = *ldwork - lwa - (*n << 1);
#line 526 "SB02SD.f"
	dgees_(sjob, "Not ordered", (L_fp)select_, n, &t[t_offset], ldt, &
		sdim, &dwork[lwa + 1], &dwork[lwa + *n + 1], &u[u_offset], 
		ldu, &dwork[lwa + (*n << 1) + 1], &i__1, bwork, info, (ftnlen)
		1, (ftnlen)11);
#line 529 "SB02SD.f"
	if (*info > 0) {
#line 530 "SB02SD.f"
	    if (lwa > 0) {
#line 530 "SB02SD.f"
		i__1 = *n << 1;
#line 530 "SB02SD.f"
		dcopy_(&i__1, &dwork[lwa + 1], &c__1, &dwork[1], &c__1);
#line 530 "SB02SD.f"
	    }
#line 532 "SB02SD.f"
	    return 0;
#line 533 "SB02SD.f"
	}

/* Computing MAX */
#line 535 "SB02SD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[lwa + (*n << 1) + 1] + lwa + (*
		n << 1);
#line 535 "SB02SD.f"
	wrkopt = max(i__1,i__2);
#line 536 "SB02SD.f"
    }
#line 537 "SB02SD.f"
    if (needac) {
#line 538 "SB02SD.f"
	dlacpy_("Full", n, n, &dwork[1], n, &dwork[iabs + 1], n, (ftnlen)4);
#line 539 "SB02SD.f"
	lwr = nn;
#line 540 "SB02SD.f"
    } else {
#line 541 "SB02SD.f"
	lwr = 0;
#line 542 "SB02SD.f"
    }

#line 544 "SB02SD.f"
    if (notrna) {
#line 545 "SB02SD.f"
	*(unsigned char *)tranat = 'T';
#line 546 "SB02SD.f"
    } else {
#line 547 "SB02SD.f"
	*(unsigned char *)tranat = 'N';
#line 548 "SB02SD.f"
    }
/*                         _ */
/*     Compute X*op(Ac) or X*op(T). */

#line 552 "SB02SD.f"
    if (update) {
#line 553 "SB02SD.f"
	dgemm_("NoTranspose", trana, n, n, n, &c_b18, &x[x_offset], ldx, &
		dwork[1], n, &c_b17, &dwork[ixma + 1], n, (ftnlen)11, (ftnlen)
		1);
#line 555 "SB02SD.f"
    } else {
#line 556 "SB02SD.f"
	mb01ud_("Right", trana, n, n, &c_b18, &t[t_offset], ldt, &x[x_offset],
		 ldx, &dwork[ixma + 1], n, &info2, (ftnlen)5, (ftnlen)1);
#line 558 "SB02SD.f"
    }

#line 560 "SB02SD.f"
    if (! jobe) {

/*        Estimate sepd(op(Ac),op(Ac)') = sepd(op(T),op(T)') and */
/*        norm(Theta). */
/*        Workspace LWR + MAX(3,2*N*N) + N*N, where */
/*                  LWR = N*N, if LYAPUN = 'O' and JOB = 'B', */
/*                  LWR = 0,   otherwise. */

#line 568 "SB02SD.f"
	sb03sy_("Both", trana, lyapun, n, &t[t_offset], ldt, &u[u_offset], 
		ldu, &dwork[ixma + 1], n, sepd, &thnorm, &iwork[1], &dwork[1],
		 &ixma, info, (ftnlen)4, (ftnlen)1, (ftnlen)1);

/* Computing MAX */
/* Computing MAX */
#line 572 "SB02SD.f"
	i__3 = 3, i__4 = nn << 1;
#line 572 "SB02SD.f"
	i__1 = wrkopt, i__2 = lwr + max(i__3,i__4) + nn;
#line 572 "SB02SD.f"
	wrkopt = max(i__1,i__2);

/*        Return if the equation is singular. */

#line 576 "SB02SD.f"
	if (*sepd == 0.) {
#line 577 "SB02SD.f"
	    *rcond = 0.;
#line 578 "SB02SD.f"
	    if (jobb) {
#line 578 "SB02SD.f"
		*ferr = 1.;
#line 578 "SB02SD.f"
	    }
#line 580 "SB02SD.f"
	    dwork[1] = (doublereal) wrkopt;
#line 581 "SB02SD.f"
	    return 0;
#line 582 "SB02SD.f"
	}

/*        Estimate norm(Pi). */
/*        Workspace LWR + MAX(3,2*N*N) + N*N. */

#line 587 "SB02SD.f"
	kase = 0;

/*        REPEAT */
#line 590 "SB02SD.f"
L30:
#line 591 "SB02SD.f"
	dlacon_(&nn, &dwork[ixbs + 1], &dwork[1], &iwork[1], &est, &kase);
#line 592 "SB02SD.f"
	if (kase != 0) {

/*           Select the triangular part of symmetric matrix to be used. */

#line 596 "SB02SD.f"
	    if (dlansy_("1-norm", "Upper", n, &dwork[1], n, &dwork[ixbs + 1], 
		    (ftnlen)6, (ftnlen)5) >= dlansy_("1-norm", "Lower", n, &
		    dwork[1], n, &dwork[ixbs + 1], (ftnlen)6, (ftnlen)5)) {
#line 600 "SB02SD.f"
		*(unsigned char *)loup = 'U';
#line 601 "SB02SD.f"
	    } else {
#line 602 "SB02SD.f"
		*(unsigned char *)loup = 'L';
#line 603 "SB02SD.f"
	    }
/*                                                        _   _ */
/*           Compute RHS = op(Ac)'*X*W*X*op(Ac) or op(T)'*X*W*X*op(T). */

#line 607 "SB02SD.f"
	    mb01ru_(loup, tranat, n, n, &c_b17, &c_b18, &dwork[1], n, &dwork[
		    ixma + 1], n, &dwork[1], n, &dwork[ixbs + 1], &nn, &info2,
		     (ftnlen)1, (ftnlen)1);
#line 610 "SB02SD.f"
	    i__1 = *n + 1;
#line 610 "SB02SD.f"
	    dscal_(n, &c_b51, &dwork[1], &i__1);

#line 612 "SB02SD.f"
	    if (update) {

/*              Transform the right-hand side: RHS := U'*RHS*U. */

#line 616 "SB02SD.f"
		mb01ru_(loup, "Transpose", n, n, &c_b17, &c_b18, &dwork[1], n,
			 &u[u_offset], ldu, &dwork[1], n, &dwork[ixbs + 1], &
			nn, &info2, (ftnlen)1, (ftnlen)9);
#line 619 "SB02SD.f"
		i__1 = *n + 1;
#line 619 "SB02SD.f"
		dscal_(n, &c_b51, &dwork[1], &i__1);
#line 620 "SB02SD.f"
	    }

/*           Fill in the remaining triangle of the symmetric matrix. */

#line 624 "SB02SD.f"
	    ma02ed_(loup, n, &dwork[1], n, (ftnlen)1);

#line 626 "SB02SD.f"
	    if (kase == 1) {

/*              Solve op(T)'*Y*op(T) - Y = scale*RHS. */

#line 630 "SB02SD.f"
		sb03mx_(trana, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
			dwork[ixbs + 1], &info2, (ftnlen)1);
#line 632 "SB02SD.f"
	    } else {

/*              Solve op(T)*W*op(T)' - W = scale*RHS. */

#line 636 "SB02SD.f"
		sb03mx_(tranat, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
			dwork[ixbs + 1], &info2, (ftnlen)1);
#line 638 "SB02SD.f"
	    }

#line 640 "SB02SD.f"
	    if (update) {

/*              Transform back to obtain the solution: Z := U*Z*U', with */
/*              Z = Y or Z = W. */

#line 645 "SB02SD.f"
		mb01ru_(loup, "No transpose", n, n, &c_b17, &c_b18, &dwork[1],
			 n, &u[u_offset], ldu, &dwork[1], n, &dwork[ixbs + 1],
			 &nn, &info2, (ftnlen)1, (ftnlen)12);
#line 648 "SB02SD.f"
		i__1 = *n + 1;
#line 648 "SB02SD.f"
		dscal_(n, &c_b51, &dwork[1], &i__1);

/*              Fill in the remaining triangle of the symmetric matrix. */

#line 652 "SB02SD.f"
		ma02ed_(loup, n, &dwork[1], n, (ftnlen)1);
#line 653 "SB02SD.f"
	    }
#line 654 "SB02SD.f"
	    goto L30;
#line 655 "SB02SD.f"
	}
/*        UNTIL KASE = 0 */

#line 658 "SB02SD.f"
	if (est < scale) {
#line 659 "SB02SD.f"
	    pinorm = est / scale;
#line 660 "SB02SD.f"
	} else {
#line 661 "SB02SD.f"
	    bignum = 1. / dlamch_("Safe minimum", (ftnlen)12);
#line 662 "SB02SD.f"
	    if (est < scale * bignum) {
#line 663 "SB02SD.f"
		pinorm = est / scale;
#line 664 "SB02SD.f"
	    } else {
#line 665 "SB02SD.f"
		pinorm = bignum;
#line 666 "SB02SD.f"
	    }
#line 667 "SB02SD.f"
	}

/*        Compute the 1-norm of A or T. */

#line 671 "SB02SD.f"
	if (update) {
#line 672 "SB02SD.f"
	    anorm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[1], (
		    ftnlen)6);
#line 673 "SB02SD.f"
	} else {
#line 674 "SB02SD.f"
	    anorm = dlanhs_("1-norm", n, &t[t_offset], ldt, &dwork[1], (
		    ftnlen)6);
#line 675 "SB02SD.f"
	}

/*        Compute the 1-norms of the matrices Q and G. */

#line 679 "SB02SD.f"
	qnorm = dlansy_("1-norm", uplo, n, &q[q_offset], ldq, &dwork[1], (
		ftnlen)6, (ftnlen)1);
#line 680 "SB02SD.f"
	gnorm = dlansy_("1-norm", uplo, n, &g[g_offset], ldg, &dwork[1], (
		ftnlen)6, (ftnlen)1);

/*        Estimate the reciprocal condition number. */

/* Computing MAX */
#line 684 "SB02SD.f"
	d__1 = max(*sepd,xnorm), d__1 = max(d__1,anorm);
#line 684 "SB02SD.f"
	tmax = max(d__1,gnorm);
#line 685 "SB02SD.f"
	if (tmax <= 1.) {
#line 686 "SB02SD.f"
	    temp = *sepd * xnorm;
#line 687 "SB02SD.f"
	    denom = qnorm + *sepd * anorm * thnorm + *sepd * gnorm * pinorm;
#line 689 "SB02SD.f"
	} else {
#line 690 "SB02SD.f"
	    temp = *sepd / tmax * (xnorm / tmax);
#line 691 "SB02SD.f"
	    denom = 1. / tmax * (qnorm / tmax) + *sepd / tmax * (anorm / tmax)
		     * thnorm + *sepd / tmax * (gnorm / tmax) * pinorm;
#line 694 "SB02SD.f"
	}
#line 695 "SB02SD.f"
	if (temp >= denom) {
#line 696 "SB02SD.f"
	    *rcond = 1.;
#line 697 "SB02SD.f"
	} else {
#line 698 "SB02SD.f"
	    *rcond = temp / denom;
#line 699 "SB02SD.f"
	}
#line 700 "SB02SD.f"
    }

#line 702 "SB02SD.f"
    if (! jobc) {

/*        Form a triangle of the residual matrix */
/*          R = op(Ac)'*X*op(Ac) + op(Ac)'*X*G*X*op(Ac) + Q - X, */
/*        or           _                _ _ _         _   _ */
/*          R = op(T)'*X*op(T) + op(T)'*X*G*X*op(T) + Q - X, */
/*        exploiting the symmetry. Actually, the equivalent formula */
/*          R = op(A)'*X*op(Ac) + Q - X */
/*        is used in the first case. */
/*        Workspace MAX(3,2*N*N) + 2*N*N,     if LYAPUN = 'O'; */
/*                  MAX(3,2*N*N) + 2*N*N + N, if LYAPUN = 'R'. */

#line 714 "SB02SD.f"
	dlacpy_(uplo, n, n, &q[q_offset], ldq, &dwork[ires + 1], n, (ftnlen)1)
		;
#line 715 "SB02SD.f"
	jj = ires + 1;
#line 716 "SB02SD.f"
	if (lower) {
#line 717 "SB02SD.f"
	    i__1 = *n;
#line 717 "SB02SD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 718 "SB02SD.f"
		i__2 = *n - j + 1;
#line 718 "SB02SD.f"
		daxpy_(&i__2, &c_b66, &x[j + j * x_dim1], &c__1, &dwork[jj], &
			c__1);
#line 719 "SB02SD.f"
		jj = jj + *n + 1;
#line 720 "SB02SD.f"
/* L40: */
#line 720 "SB02SD.f"
	    }
#line 721 "SB02SD.f"
	} else {
#line 722 "SB02SD.f"
	    i__1 = *n;
#line 722 "SB02SD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 723 "SB02SD.f"
		daxpy_(&j, &c_b66, &x[j * x_dim1 + 1], &c__1, &dwork[jj], &
			c__1);
#line 724 "SB02SD.f"
		jj += *n;
#line 725 "SB02SD.f"
/* L50: */
#line 725 "SB02SD.f"
	    }
#line 726 "SB02SD.f"
	}

#line 728 "SB02SD.f"
	if (update) {
#line 729 "SB02SD.f"
	    mb01rx_("Left", uplo, tranat, n, n, &c_b18, &c_b18, &dwork[ires + 
		    1], n, &a[a_offset], lda, &dwork[ixma + 1], n, &info2, (
		    ftnlen)4, (ftnlen)1, (ftnlen)1);
#line 732 "SB02SD.f"
	} else {
#line 733 "SB02SD.f"
	    mb01ry_("Left", uplo, tranat, n, &c_b18, &c_b18, &dwork[ires + 1],
		     n, &t[t_offset], ldt, &dwork[ixma + 1], n, &dwork[iwrk + 
		    1], &info2, (ftnlen)4, (ftnlen)1, (ftnlen)1);
#line 736 "SB02SD.f"
	    dsymm_("Left", uplo, n, n, &c_b18, &g[g_offset], ldg, &dwork[ixma 
		    + 1], n, &c_b17, &dwork[ixbs + 1], n, (ftnlen)4, (ftnlen)
		    1);
#line 738 "SB02SD.f"
	    mb01rx_("Left", uplo, "Transpose", n, n, &c_b18, &c_b18, &dwork[
		    ires + 1], n, &dwork[ixma + 1], n, &dwork[ixbs + 1], n, &
		    info2, (ftnlen)4, (ftnlen)1, (ftnlen)9);
#line 741 "SB02SD.f"
	}

/*        Get the machine precision. */

#line 745 "SB02SD.f"
	eps = dlamch_("Epsilon", (ftnlen)7);
#line 746 "SB02SD.f"
	epsn = eps * (doublereal) (*n + 4);
#line 747 "SB02SD.f"
	epst = eps * (doublereal) (*n + 1 << 1);
#line 748 "SB02SD.f"
	temp = eps * 4.;

/*        Add to abs(R) a term that takes account of rounding errors in */
/*        forming R: */
/*         abs(R) := abs(R) + EPS*(4*abs(Q) + 4*abs(X) + */
/*                   (n+4)*abs(op(Ac))'*abs(X)*abs(op(Ac)) + 2*(n+1)* */
/*                   abs(op(Ac))'*abs(X)*abs(G)*abs(X)*abs(op(Ac))), */
/*        or                             _          _ */
/*         abs(R) := abs(R) + EPS*(4*abs(Q) + 4*abs(X) + */
/*                                         _ */
/*                   (n+4)*abs(op(T))'*abs(X)*abs(op(T)) + */
/*                                         _      _      _ */
/*                 2*(n+1)*abs(op(T))'*abs(X)*abs(G)*abs(X)*abs(op(T))), */
/*        where EPS is the machine precision. */

#line 763 "SB02SD.f"
	i__1 = *n;
#line 763 "SB02SD.f"
	for (j = 1; j <= i__1; ++j) {
#line 764 "SB02SD.f"
	    i__2 = *n;
#line 764 "SB02SD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 765 "SB02SD.f"
		dwork[ixbs + (j - 1) * *n + i__] = (d__1 = x[i__ + j * x_dim1]
			, abs(d__1));
#line 766 "SB02SD.f"
/* L60: */
#line 766 "SB02SD.f"
	    }
#line 767 "SB02SD.f"
/* L70: */
#line 767 "SB02SD.f"
	}

#line 769 "SB02SD.f"
	if (lower) {
#line 770 "SB02SD.f"
	    i__1 = *n;
#line 770 "SB02SD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 771 "SB02SD.f"
		i__2 = *n;
#line 771 "SB02SD.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 772 "SB02SD.f"
		    dwork[ires + (j - 1) * *n + i__] = temp * ((d__1 = q[i__ 
			    + j * q_dim1], abs(d__1)) + (d__2 = x[i__ + j * 
			    x_dim1], abs(d__2))) + (d__3 = dwork[ires + (j - 
			    1) * *n + i__], abs(d__3));
#line 775 "SB02SD.f"
/* L80: */
#line 775 "SB02SD.f"
		}
#line 776 "SB02SD.f"
/* L90: */
#line 776 "SB02SD.f"
	    }
#line 777 "SB02SD.f"
	} else {
#line 778 "SB02SD.f"
	    i__1 = *n;
#line 778 "SB02SD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 779 "SB02SD.f"
		i__2 = j;
#line 779 "SB02SD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 780 "SB02SD.f"
		    dwork[ires + (j - 1) * *n + i__] = temp * ((d__1 = q[i__ 
			    + j * q_dim1], abs(d__1)) + (d__2 = x[i__ + j * 
			    x_dim1], abs(d__2))) + (d__3 = dwork[ires + (j - 
			    1) * *n + i__], abs(d__3));
#line 783 "SB02SD.f"
/* L100: */
#line 783 "SB02SD.f"
		}
#line 784 "SB02SD.f"
/* L110: */
#line 784 "SB02SD.f"
	    }
#line 785 "SB02SD.f"
	}

#line 787 "SB02SD.f"
	if (update) {

#line 789 "SB02SD.f"
	    i__1 = *n;
#line 789 "SB02SD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 790 "SB02SD.f"
		i__2 = *n;
#line 790 "SB02SD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 791 "SB02SD.f"
		    dwork[iabs + (j - 1) * *n + i__] = (d__1 = dwork[iabs + (
			    j - 1) * *n + i__], abs(d__1));
#line 793 "SB02SD.f"
/* L120: */
#line 793 "SB02SD.f"
		}
#line 794 "SB02SD.f"
/* L130: */
#line 794 "SB02SD.f"
	    }

#line 796 "SB02SD.f"
	    dgemm_("NoTranspose", trana, n, n, n, &c_b18, &dwork[ixbs + 1], n,
		     &dwork[iabs + 1], n, &c_b17, &dwork[ixma + 1], n, (
		    ftnlen)11, (ftnlen)1);
#line 799 "SB02SD.f"
	    mb01rx_("Left", uplo, tranat, n, n, &c_b18, &epsn, &dwork[ires + 
		    1], n, &dwork[iabs + 1], n, &dwork[ixma + 1], n, &info2, (
		    ftnlen)4, (ftnlen)1, (ftnlen)1);
#line 802 "SB02SD.f"
	} else {

#line 804 "SB02SD.f"
	    i__1 = *n;
#line 804 "SB02SD.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 805 "SB02SD.f"
		i__3 = j + 1;
#line 805 "SB02SD.f"
		i__2 = min(i__3,*n);
#line 805 "SB02SD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 806 "SB02SD.f"
		    dwork[iabs + (j - 1) * *n + i__] = (d__1 = t[i__ + j * 
			    t_dim1], abs(d__1));
#line 807 "SB02SD.f"
/* L140: */
#line 807 "SB02SD.f"
		}
#line 808 "SB02SD.f"
/* L150: */
#line 808 "SB02SD.f"
	    }

#line 810 "SB02SD.f"
	    mb01ud_("Right", trana, n, n, &c_b18, &dwork[iabs + 1], n, &dwork[
		    ixbs + 1], n, &dwork[ixma + 1], n, &info2, (ftnlen)5, (
		    ftnlen)1);
#line 812 "SB02SD.f"
	    mb01ry_("Left", uplo, tranat, n, &c_b18, &epsn, &dwork[ires + 1], 
		    n, &dwork[iabs + 1], n, &dwork[ixma + 1], n, &dwork[iwrk 
		    + 1], &info2, (ftnlen)4, (ftnlen)1, (ftnlen)1);
#line 815 "SB02SD.f"
	}

#line 817 "SB02SD.f"
	if (lower) {
#line 818 "SB02SD.f"
	    i__1 = *n;
#line 818 "SB02SD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 819 "SB02SD.f"
		i__2 = *n;
#line 819 "SB02SD.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 820 "SB02SD.f"
		    dwork[iabs + (j - 1) * *n + i__] = (d__1 = g[i__ + j * 
			    g_dim1], abs(d__1));
#line 821 "SB02SD.f"
/* L160: */
#line 821 "SB02SD.f"
		}
#line 822 "SB02SD.f"
/* L170: */
#line 822 "SB02SD.f"
	    }
#line 823 "SB02SD.f"
	} else {
#line 824 "SB02SD.f"
	    i__1 = *n;
#line 824 "SB02SD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 825 "SB02SD.f"
		i__2 = j;
#line 825 "SB02SD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 826 "SB02SD.f"
		    dwork[iabs + (j - 1) * *n + i__] = (d__1 = g[i__ + j * 
			    g_dim1], abs(d__1));
#line 827 "SB02SD.f"
/* L180: */
#line 827 "SB02SD.f"
		}
#line 828 "SB02SD.f"
/* L190: */
#line 828 "SB02SD.f"
	    }
#line 829 "SB02SD.f"
	}

#line 831 "SB02SD.f"
	if (update) {
#line 832 "SB02SD.f"
	    mb01ru_(uplo, tranat, n, n, &c_b18, &epst, &dwork[ires + 1], n, &
		    dwork[ixma + 1], n, &dwork[iabs + 1], n, &dwork[ixbs + 1],
		     &nn, &info2, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
/* Computing MAX */
#line 835 "SB02SD.f"
	    i__3 = 3, i__4 = nn << 1;
#line 835 "SB02SD.f"
	    i__1 = wrkopt, i__2 = max(i__3,i__4) + (nn << 1);
#line 835 "SB02SD.f"
	    wrkopt = max(i__1,i__2);
#line 836 "SB02SD.f"
	} else {
#line 837 "SB02SD.f"
	    dsymm_("Left", uplo, n, n, &c_b18, &dwork[iabs + 1], n, &dwork[
		    ixma + 1], n, &c_b17, &dwork[ixbs + 1], n, (ftnlen)4, (
		    ftnlen)1);
#line 839 "SB02SD.f"
	    mb01ry_("Left", uplo, tranat, n, &c_b18, &epst, &dwork[ires + 1], 
		    n, &dwork[ixma + 1], n, &dwork[ixbs + 1], n, &dwork[iwrk 
		    + 1], &info2, (ftnlen)4, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
/* Computing MAX */
#line 842 "SB02SD.f"
	    i__3 = 3, i__4 = nn << 1;
#line 842 "SB02SD.f"
	    i__1 = wrkopt, i__2 = max(i__3,i__4) + (nn << 1) + *n;
#line 842 "SB02SD.f"
	    wrkopt = max(i__1,i__2);
#line 843 "SB02SD.f"
	}

/*        Compute forward error bound, using matrix norm estimator. */
/*        Workspace MAX(3,2*N*N) + N*N. */

#line 848 "SB02SD.f"
	xanorm = dlansy_("Max", uplo, n, &x[x_offset], ldx, &dwork[1], (
		ftnlen)3, (ftnlen)1);

#line 850 "SB02SD.f"
	sb03sx_(trana, uplo, lyapun, n, &xanorm, &t[t_offset], ldt, &u[
		u_offset], ldu, &dwork[ires + 1], n, ferr, &iwork[1], &dwork[
		ixbs + 1], &ixma, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 853 "SB02SD.f"
    }

#line 855 "SB02SD.f"
    dwork[1] = (doublereal) wrkopt;
#line 856 "SB02SD.f"
    return 0;

/* *** Last line of SB02SD *** */
} /* sb02sd_ */

