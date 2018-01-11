#line 1 "SB02PD.f"
/* SB02PD.f -- translated by f2c (version 20100827).
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

#line 1 "SB02PD.f"
/* Table of constant values */

static integer c__0 = 0;
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b57 = -.5;
static doublereal c_b66 = 0.;
static doublereal c_b67 = 1.;
static doublereal c_b81 = -1.;

/* Subroutine */ int sb02pd_(char *job, char *trana, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *g, integer *ldg, doublereal *
	q, integer *ldq, doublereal *x, integer *ldx, doublereal *rcond, 
	doublereal *ferr, doublereal *wr, doublereal *wi, integer *iwork, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen job_len, 
	ftnlen trana_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, g_dim1, g_offset, q_dim1, q_offset, x_dim1, 
	    x_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, n2, ib, ic, ij, ji, ir, it, iu, ij1, ij2, iaf;
    static logical all;
    static integer ibr, ini, ifr;
    static doublereal eps, sep, tol;
    static integer isv, iscl, sdim, itau, iter;
    static doublereal conv, temp;
    static integer iwrk;
    static char loup[1];
    static integer info2;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    ma02ed_(char *, integer *, doublereal *, integer *, ftnlen), 
	    dscal_(integer *, doublereal *, doublereal *, integer *);
    static doublereal scale;
    extern /* Subroutine */ int dgees_(char *, char *, L_fp, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen), sb02qd_(char *, char *, char *, char *
	    , char *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static char equed[1];
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static doublereal hnorm;
    static logical bwork[1], lower;
    extern /* Subroutine */ int dsymm_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen), dgeqp3_(
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);
    static doublereal gnorm2, qnorm2;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dlacpy_(char *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern logical select_();
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lwamax;
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static doublereal hinnrm;
    extern /* Subroutine */ int dgesvx_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, char 
	    *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    static logical notrna;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer minwrk;
    extern /* Subroutine */ int dsytrf_(char *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     dsytri_(char *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, ftnlen);


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

/*     To solve the real continuous-time matrix algebraic Riccati */
/*     equation */

/*        op(A)'*X + X*op(A) + Q - X*G*X = 0, */

/*     where op(A) = A or A' = A**T and G, Q are symmetric (G = G**T, */
/*     Q = Q**T). The matrices A, G and Q are N-by-N and the solution X */
/*     is an N-by-N symmetric matrix. */

/*     An error bound on the solution and a condition estimate are also */
/*     optionally provided. */

/*     It is assumed that the matrices A, G and Q are such that the */
/*     corresponding Hamiltonian matrix has N eigenvalues with negative */
/*     real parts. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the computation to be performed, as follows: */
/*             = 'X':  Compute the solution only; */
/*             = 'A':  Compute all: the solution, reciprocal condition */
/*                     number, and the error bound. */

/*     TRANA   CHARACTER*1 */
/*             Specifies the option op(A): */
/*             = 'N':  op(A) = A    (No transpose); */
/*             = 'T':  op(A) = A**T (Transpose); */
/*             = 'C':  op(A) = A**T (Conjugate transpose = Transpose). */

/*     UPLO    CHARACTER*1 */
/*             Specifies which triangle of the matrices G and Q is */
/*             stored, as follows: */
/*             = 'U':  Upper triangles of G and Q are stored; */
/*             = 'L':  Lower triangles of G and Q are stored. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, G, Q, and X.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             coefficient matrix A of the equation. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     G       (input) DOUBLE PRECISION array, dimension (LDG,N) */
/*             If UPLO = 'U', the leading N-by-N upper triangular part of */
/*             this array must contain the upper triangular part of the */
/*             matrix G. */
/*             If UPLO = 'L', the leading N-by-N lower triangular part of */
/*             this array must contain the lower triangular part of the */
/*             matrix G. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G.  LDG >= max(1,N). */

/*     Q       (input) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             If UPLO = 'U', the leading N-by-N upper triangular part of */
/*             this array must contain the upper triangular part of the */
/*             matrix Q. */
/*             If UPLO = 'L', the leading N-by-N lower triangular part of */
/*             this array must contain the lower triangular part of the */
/*             matrix Q. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q.  LDQ >= max(1,N). */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             If INFO = 0, INFO = 2, or INFO = 4, the leading N-by-N */
/*             part of this array contains the symmetric solution matrix */
/*             X of the algebraic Riccati equation. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X.  LDX >= max(1,N). */

/*     RCOND   (output) DOUBLE PRECISION */
/*             If JOB = 'A', the estimate of the reciprocal condition */
/*             number of the Riccati equation. */

/*     FERR    (output) DOUBLE PRECISION */
/*             If JOB = 'A', the estimated forward error bound for the */
/*             solution X. If XTRUE is the true solution, FERR bounds the */
/*             magnitude of the largest entry in (X - XTRUE) divided by */
/*             the magnitude of the largest entry in X. */

/*     WR      (output) DOUBLE PRECISION array, dimension (N) */
/*     WI      (output) DOUBLE PRECISION array, dimension (N) */
/*             If JOB = 'A' and TRANA = 'N', WR and WI contain the real */
/*             and imaginary parts, respectively, of the eigenvalues of */
/*             the matrix A - G*X, i.e., the closed-loop system poles. */
/*             If JOB = 'A' and TRANA = 'T' or 'C', WR and WI contain the */
/*             real and imaginary parts, respectively, of the eigenvalues */
/*             of the matrix A - X*G, i.e., the closed-loop system poles. */
/*             If JOB = 'X', these arrays are not referenced. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK), where */
/*             LIWORK >= 2*N,          if JOB = 'X'; */
/*             LIWORK >= max(2*N,N*N), if JOB = 'A'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0 or INFO = 2, DWORK(1) contains the */
/*             optimal value of LDWORK. If JOB = 'A', then DWORK(2:N*N+1) */
/*             and DWORK(N*N+2:2*N*N+1) contain a real Schur form of the */
/*             closed-loop system matrix, Ac = A - G*X (if TRANA = 'N') */
/*             or Ac = A - X*G (if TRANA = 'T' or 'C'), and the */
/*             orthogonal matrix which reduced Ac to real Schur form, */
/*             respectively. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= 4*N*N + 8*N + 1,               if JOB = 'X'; */
/*             LDWORK >= max( 4*N*N + 8*N, 6*N*N ) + 1, if JOB = 'A'. */
/*             For good performance, LDWORK should be larger, e.g., */
/*             LDWORK >= 4*N*N + 6*N +( 2*N+1 )*NB,     if JOB = 'X', */
/*             where NB is the optimal blocksize. */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the Hamiltonian matrix has eigenvalues on the */
/*                   imaginary axis, so the solution and error bounds */
/*                   could not be computed; */
/*             = 2:  the iteration for the matrix sign function failed to */
/*                   converge after 50 iterations, but an approximate */
/*                   solution and error bounds (if JOB = 'A') have been */
/*                   computed; */
/*             = 3:  the system of linear equations for the solution is */
/*                   singular to working precision, so the solution and */
/*                   error bounds could not be computed; */
/*             = 4:  the matrix A-G*X (or A-X*G) cannot be reduced to */
/*                   Schur canonical form and condition number estimate */
/*                   and forward error estimate have not been computed. */

/*     METHOD */

/*     The Riccati equation is solved by the matrix sign function */
/*     approach [1], [2], implementing a scaling which enhances the */
/*     numerical stability [4]. */

/*     REFERENCES */

/*     [1] Bai, Z., Demmel, J., Dongarra, J., Petitet, A., Robinson, H., */
/*         and Stanley, K. */
/*         The spectral decomposition of nonsymmetric matrices on */
/*         distributed memory parallel computers. */
/*         SIAM J. Sci. Comput., vol. 18, pp. 1446-1461, 1997. */

/*     [2] Byers, R., He, C., and Mehrmann, V. */
/*         The matrix sign function method and the computation of */
/*         invariant subspaces. */
/*         SIAM J. Matrix Anal. Appl., vol. 18, pp. 615-632, 1997. */

/*     [3] Higham, N.J. */
/*         Perturbation theory and backward error for AX-XB=C. */
/*         BIT, vol. 33, pp. 124-136, 1993. */

/*     [4] Petkov, P.Hr., Konstantinov, M.M., and Mehrmann, V., */
/*         DGRSVX and DMSRIC: Fortran 77 subroutines for solving */
/*         continuous-time matrix algebraic Riccati equations with */
/*         condition and accuracy estimates. */
/*         Preprint SFB393/98-16, Fak. f. Mathematik, Technical */
/*         University Chemnitz, May 1998. */

/*     NUMERICAL ASPECTS */

/*     The solution accuracy can be controlled by the output parameter */
/*     FERR. */

/*     FURTHER COMMENTS */

/*     The condition number of the Riccati equation is estimated as */

/*     cond = ( norm(Theta)*norm(A) + norm(inv(Omega))*norm(Q) + */
/*                 norm(Pi)*norm(G) ) / norm(X), */

/*     where Omega, Theta and Pi are linear operators defined by */

/*     Omega(W) = op(Ac)'*W + W*op(Ac), */
/*     Theta(W) = inv(Omega(op(W)'*X + X*op(W))), */
/*        Pi(W) = inv(Omega(X*W*X)), */

/*     and the matrix Ac (the closed-loop system matrix) is given by */
/*        Ac = A - G*X, if TRANA = 'N', or */
/*        Ac = A - X*G, if TRANA = 'T' or 'C'. */

/*     The program estimates the quantities */

/*     sep(op(Ac),-op(Ac)') = 1 / norm(inv(Omega)), */

/*     norm(Theta) and norm(Pi) using 1-norm condition estimator. */

/*     The forward error bound is estimated using a practical error bound */
/*     similar to the one proposed in [3]. */

/*     CONTRIBUTOR */

/*     P. Petkov, Tech. University of Sofia, March 2000. */

/*     REVISIONS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, June 2000. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, continuous-time system, */
/*     optimal control, optimal regulator. */

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

#line 291 "SB02PD.f"
    /* Parameter adjustments */
#line 291 "SB02PD.f"
    a_dim1 = *lda;
#line 291 "SB02PD.f"
    a_offset = 1 + a_dim1;
#line 291 "SB02PD.f"
    a -= a_offset;
#line 291 "SB02PD.f"
    g_dim1 = *ldg;
#line 291 "SB02PD.f"
    g_offset = 1 + g_dim1;
#line 291 "SB02PD.f"
    g -= g_offset;
#line 291 "SB02PD.f"
    q_dim1 = *ldq;
#line 291 "SB02PD.f"
    q_offset = 1 + q_dim1;
#line 291 "SB02PD.f"
    q -= q_offset;
#line 291 "SB02PD.f"
    x_dim1 = *ldx;
#line 291 "SB02PD.f"
    x_offset = 1 + x_dim1;
#line 291 "SB02PD.f"
    x -= x_offset;
#line 291 "SB02PD.f"
    --wr;
#line 291 "SB02PD.f"
    --wi;
#line 291 "SB02PD.f"
    --iwork;
#line 291 "SB02PD.f"
    --dwork;
#line 291 "SB02PD.f"

#line 291 "SB02PD.f"
    /* Function Body */
#line 291 "SB02PD.f"
    all = lsame_(job, "A", (ftnlen)1, (ftnlen)1);
#line 292 "SB02PD.f"
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
#line 293 "SB02PD.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);

#line 295 "SB02PD.f"
    *info = 0;
#line 296 "SB02PD.f"
    if (! all && ! lsame_(job, "X", (ftnlen)1, (ftnlen)1)) {
#line 297 "SB02PD.f"
	*info = -1;
#line 298 "SB02PD.f"
    } else if (! lsame_(trana, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trana, 
	    "C", (ftnlen)1, (ftnlen)1) && ! notrna) {
#line 300 "SB02PD.f"
	*info = -2;
#line 301 "SB02PD.f"
    } else if (! lower && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 302 "SB02PD.f"
	*info = -3;
#line 303 "SB02PD.f"
    } else if (*n < 0) {
#line 304 "SB02PD.f"
	*info = -4;
#line 305 "SB02PD.f"
    } else if (*lda < max(1,*n)) {
#line 306 "SB02PD.f"
	*info = -6;
#line 307 "SB02PD.f"
    } else if (*ldg < max(1,*n)) {
#line 308 "SB02PD.f"
	*info = -8;
#line 309 "SB02PD.f"
    } else if (*ldq < max(1,*n)) {
#line 310 "SB02PD.f"
	*info = -10;
#line 311 "SB02PD.f"
    } else if (*ldx < max(1,*n)) {
#line 312 "SB02PD.f"
	*info = -12;
#line 313 "SB02PD.f"
    } else {

/*        Compute workspace. */

#line 317 "SB02PD.f"
	if (all) {
/* Computing MAX */
#line 318 "SB02PD.f"
	    i__1 = (*n << 2) * *n + (*n << 3) + 1, i__2 = *n * 6 * *n;
#line 318 "SB02PD.f"
	    minwrk = max(i__1,i__2);
#line 319 "SB02PD.f"
	} else {
#line 320 "SB02PD.f"
	    minwrk = (*n << 2) * *n + (*n << 3) + 1;
#line 321 "SB02PD.f"
	}
#line 322 "SB02PD.f"
	if (*ldwork < minwrk) {
#line 323 "SB02PD.f"
	    *info = -19;
#line 324 "SB02PD.f"
	}
#line 325 "SB02PD.f"
    }
#line 326 "SB02PD.f"
    if (*info != 0) {
#line 327 "SB02PD.f"
	i__1 = -(*info);
#line 327 "SB02PD.f"
	xerbla_("SB02PD", &i__1, (ftnlen)6);
#line 328 "SB02PD.f"
	return 0;
#line 329 "SB02PD.f"
    }

/*     Quick return if possible. */

#line 333 "SB02PD.f"
    if (*n == 0) {
#line 334 "SB02PD.f"
	if (all) {
#line 335 "SB02PD.f"
	    *rcond = 1.;
#line 336 "SB02PD.f"
	    *ferr = 0.;
#line 337 "SB02PD.f"
	}
#line 338 "SB02PD.f"
	dwork[1] = 1.;
#line 339 "SB02PD.f"
	return 0;
#line 340 "SB02PD.f"
    }

/*     Set tol. */

#line 344 "SB02PD.f"
    eps = dlamch_("P", (ftnlen)1);
#line 345 "SB02PD.f"
    tol = (doublereal) (*n) * 10. * eps;

/*     Compute the square-roots of the norms of the matrices Q and G . */

#line 349 "SB02PD.f"
    qnorm2 = sqrt(dlansy_("1", uplo, n, &q[q_offset], ldq, &dwork[1], (ftnlen)
	    1, (ftnlen)1));
#line 350 "SB02PD.f"
    gnorm2 = sqrt(dlansy_("1", uplo, n, &g[g_offset], ldg, &dwork[1], (ftnlen)
	    1, (ftnlen)1));

#line 352 "SB02PD.f"
    n2 = *n << 1;

/*     Construct the lower (if UPLO = 'L') or upper (if UPLO = 'U') */
/*     triangle of the symmetric block-permuted Hamiltonian matrix. */
/*     During iteration, both the current iterate corresponding to the */
/*     Hamiltonian matrix, and its inverse are needed. To reduce the */
/*     workspace length, the transpose of the triangle specified by UPLO */
/*     of the current iterate H is saved in the opposite triangle, */
/*     suitably shifted with one column, and then the inverse of H */
/*     overwrites H. The triangles of the saved iterate and its inverse */
/*     are stored together in an 2*N-by-(2*N+1) matrix. For instance, if */
/*     UPLO = 'U', then the upper triangle is built starting from the */
/*     location 2*N+1 of the array DWORK, so that its transpose can be */
/*     stored in the lower triangle of DWORK. */
/*     Workspace: need   4*N*N,        if UPLO = 'L'; */
/*                       4*N*N + 2*N,  if UPLO = 'U'. */

#line 369 "SB02PD.f"
    if (lower) {
#line 370 "SB02PD.f"
	ini = 0;
#line 371 "SB02PD.f"
	isv = n2;
#line 372 "SB02PD.f"
	*(unsigned char *)loup = 'U';

#line 374 "SB02PD.f"
	i__1 = *n;
#line 374 "SB02PD.f"
	for (j = 1; j <= i__1; ++j) {
#line 375 "SB02PD.f"
	    ij = (j - 1) * n2 + j;

#line 377 "SB02PD.f"
	    i__2 = *n;
#line 377 "SB02PD.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
#line 378 "SB02PD.f"
		dwork[ij] = -q[i__ + j * q_dim1];
#line 379 "SB02PD.f"
		++ij;
#line 380 "SB02PD.f"
/* L10: */
#line 380 "SB02PD.f"
	    }

#line 382 "SB02PD.f"
	    if (notrna) {

#line 384 "SB02PD.f"
		i__2 = *n;
#line 384 "SB02PD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 385 "SB02PD.f"
		    dwork[ij] = -a[i__ + j * a_dim1];
#line 386 "SB02PD.f"
		    ++ij;
#line 387 "SB02PD.f"
/* L20: */
#line 387 "SB02PD.f"
		}

#line 389 "SB02PD.f"
	    } else {

#line 391 "SB02PD.f"
		i__2 = *n;
#line 391 "SB02PD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 392 "SB02PD.f"
		    dwork[ij] = -a[j + i__ * a_dim1];
#line 393 "SB02PD.f"
		    ++ij;
#line 394 "SB02PD.f"
/* L30: */
#line 394 "SB02PD.f"
		}

#line 396 "SB02PD.f"
	    }
#line 397 "SB02PD.f"
/* L40: */
#line 397 "SB02PD.f"
	}

#line 399 "SB02PD.f"
	i__1 = *n;
#line 399 "SB02PD.f"
	for (j = 1; j <= i__1; ++j) {
#line 400 "SB02PD.f"
	    ij = (*n + j - 1) * n2 + *n + j;

#line 402 "SB02PD.f"
	    i__2 = *n;
#line 402 "SB02PD.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
#line 403 "SB02PD.f"
		dwork[ij] = g[i__ + j * g_dim1];
#line 404 "SB02PD.f"
		++ij;
#line 405 "SB02PD.f"
/* L50: */
#line 405 "SB02PD.f"
	    }

#line 407 "SB02PD.f"
/* L60: */
#line 407 "SB02PD.f"
	}

#line 409 "SB02PD.f"
    } else {
#line 410 "SB02PD.f"
	ini = n2;
#line 411 "SB02PD.f"
	isv = 0;
#line 412 "SB02PD.f"
	*(unsigned char *)loup = 'L';

#line 414 "SB02PD.f"
	i__1 = *n;
#line 414 "SB02PD.f"
	for (j = 1; j <= i__1; ++j) {
#line 415 "SB02PD.f"
	    ij = j * n2 + 1;

#line 417 "SB02PD.f"
	    i__2 = j;
#line 417 "SB02PD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 418 "SB02PD.f"
		dwork[ij] = -q[i__ + j * q_dim1];
#line 419 "SB02PD.f"
		++ij;
#line 420 "SB02PD.f"
/* L70: */
#line 420 "SB02PD.f"
	    }

#line 422 "SB02PD.f"
/* L80: */
#line 422 "SB02PD.f"
	}

#line 424 "SB02PD.f"
	i__1 = *n;
#line 424 "SB02PD.f"
	for (j = 1; j <= i__1; ++j) {
#line 425 "SB02PD.f"
	    ij = (*n + j) * n2 + 1;

#line 427 "SB02PD.f"
	    if (notrna) {

#line 429 "SB02PD.f"
		i__2 = *n;
#line 429 "SB02PD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 430 "SB02PD.f"
		    dwork[ij] = -a[j + i__ * a_dim1];
#line 431 "SB02PD.f"
		    ++ij;
#line 432 "SB02PD.f"
/* L90: */
#line 432 "SB02PD.f"
		}

#line 434 "SB02PD.f"
	    } else {

#line 436 "SB02PD.f"
		i__2 = *n;
#line 436 "SB02PD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 437 "SB02PD.f"
		    dwork[ij] = -a[i__ + j * a_dim1];
#line 438 "SB02PD.f"
		    ++ij;
#line 439 "SB02PD.f"
/* L100: */
#line 439 "SB02PD.f"
		}

#line 441 "SB02PD.f"
	    }

#line 443 "SB02PD.f"
	    i__2 = j;
#line 443 "SB02PD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 444 "SB02PD.f"
		dwork[ij] = g[i__ + j * g_dim1];
#line 445 "SB02PD.f"
		++ij;
#line 446 "SB02PD.f"
/* L110: */
#line 446 "SB02PD.f"
	    }

#line 448 "SB02PD.f"
/* L120: */
#line 448 "SB02PD.f"
	}

#line 450 "SB02PD.f"
    }

/*     Block-scaling. */

#line 454 "SB02PD.f"
    iscl = 0;
#line 455 "SB02PD.f"
    if (qnorm2 > gnorm2 && gnorm2 > 0.) {
#line 456 "SB02PD.f"
	dlascl_(uplo, &c__0, &c__0, &qnorm2, &gnorm2, n, n, &dwork[ini + 1], &
		n2, &info2, (ftnlen)1);
#line 458 "SB02PD.f"
	dlascl_(uplo, &c__0, &c__0, &gnorm2, &qnorm2, n, n, &dwork[n2 * *n + *
		n + ini + 1], &n2, &info2, (ftnlen)1);
#line 460 "SB02PD.f"
	iscl = 1;
#line 461 "SB02PD.f"
    }

/*     Workspace usage. */

#line 465 "SB02PD.f"
    itau = n2 * n2;
#line 466 "SB02PD.f"
    iwrk = itau + n2;

#line 468 "SB02PD.f"
    lwamax = n2 * ilaenv_(&c__1, "DSYTRF", uplo, &n2, &c_n1, &c_n1, &c_n1, (
	    ftnlen)6, (ftnlen)1);

/*     Compute the matrix sign function. */

#line 472 "SB02PD.f"
    for (iter = 1; iter <= 50; ++iter) {

/*        Save the transpose of the corresponding triangle of the */
/*        current iterate in the free locations of the shifted opposite */
/*        triangle. */
/*        Workspace: need   4*N*N + 2*N. */

#line 479 "SB02PD.f"
	if (lower) {

#line 481 "SB02PD.f"
	    i__1 = n2;
#line 481 "SB02PD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 482 "SB02PD.f"
		dcopy_(&i__, &dwork[i__], &n2, &dwork[i__ * n2 + 1], &c__1);
#line 483 "SB02PD.f"
/* L130: */
#line 483 "SB02PD.f"
	    }

#line 485 "SB02PD.f"
	} else {

#line 487 "SB02PD.f"
	    i__1 = n2;
#line 487 "SB02PD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 488 "SB02PD.f"
		dcopy_(&i__, &dwork[i__ * n2 + 1], &c__1, &dwork[i__], &n2);
#line 489 "SB02PD.f"
/* L140: */
#line 489 "SB02PD.f"
	    }

#line 491 "SB02PD.f"
	}

/*        Store the norm of the Hamiltonian matrix. */

#line 495 "SB02PD.f"
	hnorm = dlansy_("F", uplo, &n2, &dwork[ini + 1], &n2, &dwork[1], (
		ftnlen)1, (ftnlen)1);

/*        Compute the inverse of the block-permuted Hamiltonian matrix. */
/*        Workspace: need   4*N*N + 2*N + 1; */
/*                   prefer 4*N*N + 2*N + 2*N*NB. */

#line 501 "SB02PD.f"
	i__1 = *ldwork - iwrk;
#line 501 "SB02PD.f"
	dsytrf_(uplo, &n2, &dwork[ini + 1], &n2, &iwork[1], &dwork[iwrk + 1], 
		&i__1, &info2, (ftnlen)1);
#line 503 "SB02PD.f"
	if (info2 > 0) {
#line 504 "SB02PD.f"
	    *info = 1;
#line 505 "SB02PD.f"
	    return 0;
#line 506 "SB02PD.f"
	}

/*        Workspace: need   4*N*N + 4*N. */

#line 510 "SB02PD.f"
	dsytri_(uplo, &n2, &dwork[ini + 1], &n2, &iwork[1], &dwork[iwrk + 1], 
		&info2, (ftnlen)1);

/*        Block-permutation of the inverse matrix. */

#line 515 "SB02PD.f"
	if (lower) {

#line 517 "SB02PD.f"
	    i__1 = *n;
#line 517 "SB02PD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 518 "SB02PD.f"
		ij2 = (*n + j - 1) * n2 + *n + j;

#line 520 "SB02PD.f"
		i__2 = (j - 1) * n2 + *n;
#line 520 "SB02PD.f"
		for (ij1 = (j - 1) * n2 + j; ij1 <= i__2; ++ij1) {
#line 521 "SB02PD.f"
		    temp = dwork[ij1];
#line 522 "SB02PD.f"
		    dwork[ij1] = -dwork[ij2];
#line 523 "SB02PD.f"
		    dwork[ij2] = -temp;
#line 524 "SB02PD.f"
		    ++ij2;
#line 525 "SB02PD.f"
/* L150: */
#line 525 "SB02PD.f"
		}

#line 527 "SB02PD.f"
		i__2 = j - 1;
#line 527 "SB02PD.f"
		dswap_(&i__2, &dwork[*n + j], &n2, &dwork[(j - 1) * n2 + *n + 
			1], &c__1);
#line 529 "SB02PD.f"
/* L160: */
#line 529 "SB02PD.f"
	    }

#line 531 "SB02PD.f"
	} else {

#line 533 "SB02PD.f"
	    i__1 = *n;
#line 533 "SB02PD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 534 "SB02PD.f"
		ij2 = (*n + j) * n2 + *n + 1;

#line 536 "SB02PD.f"
		i__2 = j * n2 + j;
#line 536 "SB02PD.f"
		for (ij1 = j * n2 + 1; ij1 <= i__2; ++ij1) {
#line 537 "SB02PD.f"
		    temp = dwork[ij1];
#line 538 "SB02PD.f"
		    dwork[ij1] = -dwork[ij2];
#line 539 "SB02PD.f"
		    dwork[ij2] = -temp;
#line 540 "SB02PD.f"
		    ++ij2;
#line 541 "SB02PD.f"
/* L170: */
#line 541 "SB02PD.f"
		}

#line 543 "SB02PD.f"
		i__2 = j - 1;
#line 543 "SB02PD.f"
		dswap_(&i__2, &dwork[(*n + 1) * n2 + j], &n2, &dwork[(*n + j) 
			* n2 + 1], &c__1);
#line 545 "SB02PD.f"
/* L180: */
#line 545 "SB02PD.f"
	    }

#line 547 "SB02PD.f"
	}

/*        Scale the Hamiltonian matrix and its inverse and compute */
/*        the next iterate. */

#line 552 "SB02PD.f"
	hinnrm = dlansy_("F", uplo, &n2, &dwork[ini + 1], &n2, &dwork[1], (
		ftnlen)1, (ftnlen)1);
#line 553 "SB02PD.f"
	scale = sqrt(hinnrm / hnorm);

#line 555 "SB02PD.f"
	if (lower) {

#line 557 "SB02PD.f"
	    i__1 = n2;
#line 557 "SB02PD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 558 "SB02PD.f"
		ji = (j - 1) * n2 + j;

#line 560 "SB02PD.f"
		i__2 = j * n2;
#line 560 "SB02PD.f"
		for (ij = ji; ij <= i__2; ++ij) {
#line 561 "SB02PD.f"
		    ji += n2;
#line 562 "SB02PD.f"
		    dwork[ij] = (dwork[ij] / scale + dwork[ji] * scale) / 2.;
#line 564 "SB02PD.f"
		    dwork[ji] -= dwork[ij];
#line 565 "SB02PD.f"
/* L190: */
#line 565 "SB02PD.f"
		}

#line 567 "SB02PD.f"
/* L200: */
#line 567 "SB02PD.f"
	    }

#line 569 "SB02PD.f"
	} else {

#line 571 "SB02PD.f"
	    i__1 = n2;
#line 571 "SB02PD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 572 "SB02PD.f"
		ji = j;

#line 574 "SB02PD.f"
		i__2 = j * n2 + j;
#line 574 "SB02PD.f"
		for (ij = j * n2 + 1; ij <= i__2; ++ij) {
#line 575 "SB02PD.f"
		    dwork[ij] = (dwork[ij] / scale + dwork[ji] * scale) / 2.;
#line 577 "SB02PD.f"
		    dwork[ji] -= dwork[ij];
#line 578 "SB02PD.f"
		    ji += n2;
#line 579 "SB02PD.f"
/* L210: */
#line 579 "SB02PD.f"
		}

#line 581 "SB02PD.f"
/* L220: */
#line 581 "SB02PD.f"
	    }

#line 583 "SB02PD.f"
	}

/*        Test for convergence. */

#line 587 "SB02PD.f"
	conv = dlansy_("F", loup, &n2, &dwork[isv + 1], &n2, &dwork[1], (
		ftnlen)1, (ftnlen)1);
#line 588 "SB02PD.f"
	if (conv <= tol * hnorm) {
#line 588 "SB02PD.f"
	    goto L240;
#line 588 "SB02PD.f"
	}
#line 589 "SB02PD.f"
/* L230: */
#line 589 "SB02PD.f"
    }

/*     No convergence after MAXIT iterations, but an approximate solution */
/*     has been found. */

#line 594 "SB02PD.f"
    *info = 2;

#line 596 "SB02PD.f"
L240:

/*     If UPLO = 'U', shift the upper triangle one column to the left. */

#line 600 "SB02PD.f"
    if (! lower) {
#line 600 "SB02PD.f"
	dlacpy_("U", &n2, &n2, &dwork[ini + 1], &n2, &dwork[1], &n2, (ftnlen)
		1);
#line 600 "SB02PD.f"
    }

/*     Divide the triangle elements by -2 and then fill-in the other */
/*     triangle by symmetry. */

#line 606 "SB02PD.f"
    if (lower) {

#line 608 "SB02PD.f"
	i__1 = n2;
#line 608 "SB02PD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 609 "SB02PD.f"
	    i__2 = n2 - i__ + 1;
#line 609 "SB02PD.f"
	    dscal_(&i__2, &c_b57, &dwork[(i__ - 1) * n2 + i__], &c__1);
#line 610 "SB02PD.f"
/* L250: */
#line 610 "SB02PD.f"
	}

#line 612 "SB02PD.f"
    } else {

#line 614 "SB02PD.f"
	i__1 = n2;
#line 614 "SB02PD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 615 "SB02PD.f"
	    dscal_(&i__, &c_b57, &dwork[(i__ - 1) * n2 + 1], &c__1);
#line 616 "SB02PD.f"
/* L260: */
#line 616 "SB02PD.f"
	}

#line 618 "SB02PD.f"
    }
#line 619 "SB02PD.f"
    ma02ed_(uplo, &n2, &dwork[1], &n2, (ftnlen)1);

/*     Back block-permutation. */

#line 623 "SB02PD.f"
    i__1 = n2;
#line 623 "SB02PD.f"
    for (j = 1; j <= i__1; ++j) {

#line 625 "SB02PD.f"
	i__2 = (j - 1) * n2 + *n;
#line 625 "SB02PD.f"
	for (i__ = (j - 1) * n2 + 1; i__ <= i__2; ++i__) {
#line 626 "SB02PD.f"
	    temp = dwork[i__];
#line 627 "SB02PD.f"
	    dwork[i__] = -dwork[i__ + *n];
#line 628 "SB02PD.f"
	    dwork[i__ + *n] = temp;
#line 629 "SB02PD.f"
/* L270: */
#line 629 "SB02PD.f"
	}

#line 631 "SB02PD.f"
/* L280: */
#line 631 "SB02PD.f"
    }

/*     Compute the QR decomposition of the projector onto the stable */
/*     invariant subspace. */
/*     Workspace: need   4*N*N + 8*N + 1. */
/*                prefer 4*N*N + 6*N + ( 2*N+1 )*NB. */

#line 638 "SB02PD.f"
    i__1 = n2;
#line 638 "SB02PD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 639 "SB02PD.f"
	iwork[i__] = 0;
#line 640 "SB02PD.f"
	dwork[(i__ - 1) * n2 + i__] += .5;
#line 641 "SB02PD.f"
/* L290: */
#line 641 "SB02PD.f"
    }

#line 643 "SB02PD.f"
    i__1 = *ldwork - iwrk;
#line 643 "SB02PD.f"
    dgeqp3_(&n2, &n2, &dwork[1], &n2, &iwork[1], &dwork[itau + 1], &dwork[
	    iwrk + 1], &i__1, &info2);
/* Computing MAX */
#line 645 "SB02PD.f"
    i__1 = (integer) dwork[iwrk + 1];
#line 645 "SB02PD.f"
    lwamax = max(i__1,lwamax);

/*     Accumulate the orthogonal transformations. Note that only the */
/*     first N columns of the array DWORK, returned by DGEQP3, are */
/*     needed, so that the last N columns of DWORK are used to get the */
/*     orthogonal basis for the stable invariant subspace. */
/*     Workspace: need   4*N*N + 3*N. */
/*                prefer 4*N*N + 2*N + N*NB. */

#line 654 "SB02PD.f"
    ib = *n * *n;
#line 655 "SB02PD.f"
    iaf = n2 * *n;
#line 656 "SB02PD.f"
    dlaset_("F", &n2, n, &c_b66, &c_b67, &dwork[iaf + 1], &n2, (ftnlen)1);
#line 657 "SB02PD.f"
    i__1 = *ldwork - iwrk;
#line 657 "SB02PD.f"
    dormqr_("L", "N", &n2, n, n, &dwork[1], &n2, &dwork[itau + 1], &dwork[iaf 
	    + 1], &n2, &dwork[iwrk + 1], &i__1, &info2, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 660 "SB02PD.f"
    i__1 = (integer) dwork[iwrk + 1];
#line 660 "SB02PD.f"
    lwamax = iwrk + max(i__1,lwamax);

/*     Store the matrices V11 and V21' . */

#line 664 "SB02PD.f"
    dlacpy_("F", n, n, &dwork[iaf + 1], &n2, &dwork[1], n, (ftnlen)1);
#line 665 "SB02PD.f"
    ma02ad_("F", n, n, &dwork[iaf + *n + 1], &n2, &dwork[ib + 1], n, (ftnlen)
	    1);

#line 667 "SB02PD.f"
    ir = iaf + ib;
#line 668 "SB02PD.f"
    ic = ir + *n;
#line 669 "SB02PD.f"
    ifr = ic + *n;
#line 670 "SB02PD.f"
    ibr = ifr + *n;
#line 671 "SB02PD.f"
    iwrk = ibr + *n;

/*     Compute the solution matrix X . */
/*     Workspace: need   3*N*N + 8*N. */

#line 676 "SB02PD.f"
    dgesvx_("E", "T", n, n, &dwork[1], n, &dwork[iaf + 1], n, &iwork[1], 
	    equed, &dwork[ir + 1], &dwork[ic + 1], &dwork[ib + 1], n, &x[
	    x_offset], ldx, rcond, &dwork[ifr + 1], &dwork[ibr + 1], &dwork[
	    iwrk + 1], &iwork[*n + 1], &info2, (ftnlen)1, (ftnlen)1, (ftnlen)
	    1);
#line 681 "SB02PD.f"
    if (info2 > 0) {
#line 682 "SB02PD.f"
	*info = 3;
#line 683 "SB02PD.f"
	return 0;
#line 684 "SB02PD.f"
    }

/*     Symmetrize the solution. */

#line 688 "SB02PD.f"
    i__1 = *n - 1;
#line 688 "SB02PD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

#line 690 "SB02PD.f"
	i__2 = *n;
#line 690 "SB02PD.f"
	for (j = i__ + 1; j <= i__2; ++j) {
#line 691 "SB02PD.f"
	    temp = (x[i__ + j * x_dim1] + x[j + i__ * x_dim1]) / 2.;
#line 692 "SB02PD.f"
	    x[i__ + j * x_dim1] = temp;
#line 693 "SB02PD.f"
	    x[j + i__ * x_dim1] = temp;
#line 694 "SB02PD.f"
/* L300: */
#line 694 "SB02PD.f"
	}

#line 696 "SB02PD.f"
/* L310: */
#line 696 "SB02PD.f"
    }

/*     Undo scaling for the solution matrix. */

#line 700 "SB02PD.f"
    if (iscl == 1) {
#line 701 "SB02PD.f"
	dlascl_("G", &c__0, &c__0, &gnorm2, &qnorm2, n, n, &x[x_offset], ldx, 
		&info2, (ftnlen)1);
#line 702 "SB02PD.f"
    }

#line 704 "SB02PD.f"
    if (all) {

/*        Compute the estimates of the reciprocal condition number and */
/*        error bound. */
/*        Workspace usage. */

#line 710 "SB02PD.f"
	it = 1;
#line 711 "SB02PD.f"
	iu = it + *n * *n;
#line 712 "SB02PD.f"
	iwrk = iu + *n * *n;

#line 714 "SB02PD.f"
	dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[it + 1], n, (ftnlen)4)
		;
#line 715 "SB02PD.f"
	if (notrna) {

/*           Compute Ac = A-G*X . */

#line 719 "SB02PD.f"
	    dsymm_("L", uplo, n, n, &c_b81, &g[g_offset], ldg, &x[x_offset], 
		    ldx, &c_b67, &dwork[it + 1], n, (ftnlen)1, (ftnlen)1);
#line 721 "SB02PD.f"
	} else {

/*           Compute Ac = A-X*G . */

#line 725 "SB02PD.f"
	    dsymm_("R", uplo, n, n, &c_b81, &g[g_offset], ldg, &x[x_offset], 
		    ldx, &c_b67, &dwork[it + 1], n, (ftnlen)1, (ftnlen)1);
#line 727 "SB02PD.f"
	}

/*        Compute the Schur factorization of Ac . */
/*        Workspace: need   2*N*N + 5*N + 1; */
/*                   prefer larger. */

#line 733 "SB02PD.f"
	i__1 = *ldwork - iwrk;
#line 733 "SB02PD.f"
	dgees_("V", "N", (L_fp)select_, n, &dwork[it + 1], n, &sdim, &wr[1], &
		wi[1], &dwork[iu + 1], n, &dwork[iwrk + 1], &i__1, bwork, &
		info2, (ftnlen)1, (ftnlen)1);
#line 736 "SB02PD.f"
	if (info2 > 0) {
#line 737 "SB02PD.f"
	    *info = 4;
#line 738 "SB02PD.f"
	    return 0;
#line 739 "SB02PD.f"
	}
/* Computing MAX */
#line 740 "SB02PD.f"
	i__1 = (integer) dwork[iwrk + 1];
#line 740 "SB02PD.f"
	lwamax = iwrk + max(i__1,lwamax);

/*        Estimate the reciprocal condition number and the forward error. */
/*        Workspace: need   6*N*N + 1; */
/*                   prefer larger. */

#line 746 "SB02PD.f"
	i__1 = *ldwork - iwrk;
#line 746 "SB02PD.f"
	sb02qd_("B", "F", trana, uplo, "O", n, &a[a_offset], lda, &dwork[it + 
		1], n, &dwork[iu + 1], n, &g[g_offset], ldg, &q[q_offset], 
		ldq, &x[x_offset], ldx, &sep, rcond, ferr, &iwork[1], &dwork[
		iwrk + 1], &i__1, &info2, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 750 "SB02PD.f"
	i__1 = (integer) dwork[iwrk + 1];
#line 750 "SB02PD.f"
	lwamax = iwrk + max(i__1,lwamax);
#line 751 "SB02PD.f"
    }

#line 753 "SB02PD.f"
    dwork[1] = (doublereal) lwamax;
#line 754 "SB02PD.f"
    return 0;
/* *** Last line of SB02PD */
} /* sb02pd_ */

