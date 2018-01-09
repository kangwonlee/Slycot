#line 1 "SG03AD.f"
/* SG03AD.f -- translated by f2c (version 20100827).
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

#line 1 "SG03AD.f"
/* Table of constant values */

static doublereal c_b20 = 0.;
static doublereal c_b21 = 1.;
static integer c__1 = 1;

/* Subroutine */ int sg03ad_(char *dico, char *job, char *fact, char *trans, 
	char *uplo, integer *n, doublereal *a, integer *lda, doublereal *e, 
	integer *lde, doublereal *q, integer *ldq, doublereal *z__, integer *
	ldz, doublereal *x, integer *ldx, doublereal *scale, doublereal *sep, 
	doublereal *ferr, doublereal *alphar, doublereal *alphai, doublereal *
	beta, integer *iwork, doublereal *dwork, integer *ldwork, integer *
	info, ftnlen dico_len, ftnlen job_len, ftnlen fact_len, ftnlen 
	trans_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, e_dim1, e_offset, q_dim1, q_offset, x_dim1, 
	    x_offset, z_dim1, z_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__;
    static doublereal eps, est;
    static integer kase, info1;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int mb01rd_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), dgegs_(char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sg03ax_(char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), sg03ay_(char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), mb01rw_(char *, char 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal norma;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal norme;
    static logical wantx;
    static doublereal scale1;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacon_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *);
    static logical isfact;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical isdisc, wantbh;
    static char etrans[1];
    static logical istran;
    static integer minwrk;
    static logical wantsp, isuppr;
    static integer optwrk;


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

/*     To solve for X either the generalized continuous-time Lyapunov */
/*     equation */

/*             T                T */
/*        op(A)  X op(E) + op(E)  X op(A) = SCALE * Y,                (1) */

/*     or the generalized discrete-time Lyapunov equation */

/*             T                T */
/*        op(A)  X op(A) - op(E)  X op(E) = SCALE * Y,                (2) */

/*     where op(M) is either M or M**T for M = A, E and the right hand */
/*     side Y is symmetric. A, E, Y, and the solution X are N-by-N */
/*     matrices. SCALE is an output scale factor, set to avoid overflow */
/*     in X. */

/*     Estimates of the separation and the relative forward error norm */
/*     are provided. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies which type of the equation is considered: */
/*             = 'C':  Continuous-time equation (1); */
/*             = 'D':  Discrete-time equation (2). */

/*     JOB     CHARACTER*1 */
/*             Specifies if the solution is to be computed and if the */
/*             separation is to be estimated: */
/*             = 'X':  Compute the solution only; */
/*             = 'S':  Estimate the separation only; */
/*             = 'B':  Compute the solution and estimate the separation. */

/*     FACT    CHARACTER*1 */
/*             Specifies whether the generalized real Schur */
/*             factorization of the pencil A - lambda * E is supplied */
/*             on entry or not: */
/*             = 'N':  Factorization is not supplied; */
/*             = 'F':  Factorization is supplied. */

/*     TRANS   CHARACTER*1 */
/*             Specifies whether the transposed equation is to be solved */
/*             or not: */
/*             = 'N':  op(A) = A,    op(E) = E; */
/*             = 'T':  op(A) = A**T, op(E) = E**T. */

/*     UPLO    CHARACTER*1 */
/*             Specifies whether the lower or the upper triangle of the */
/*             array X is needed on input: */
/*             = 'L':  Only the lower triangle is needed on input; */
/*             = 'U':  Only the upper triangle is needed on input. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, if FACT = 'F', then the leading N-by-N upper */
/*             Hessenberg part of this array must contain the */
/*             generalized Schur factor A_s of the matrix A (see */
/*             definition (3) in section METHOD). A_s must be an upper */
/*             quasitriangular matrix. The elements below the upper */
/*             Hessenberg part of the array A are not referenced. */
/*             If FACT = 'N', then the leading N-by-N part of this */
/*             array must contain the matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the generalized Schur factor A_s of the matrix A. (A_s is */
/*             an upper quasitriangular matrix.) */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, if FACT = 'F', then the leading N-by-N upper */
/*             triangular part of this array must contain the */
/*             generalized Schur factor E_s of the matrix E (see */
/*             definition (4) in section METHOD). The elements below the */
/*             upper triangular part of the array E are not referenced. */
/*             If FACT = 'N', then the leading N-by-N part of this */
/*             array must contain the coefficient matrix E of the */
/*             equation. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the generalized Schur factor E_s of the matrix E. (E_s is */
/*             an upper triangular matrix.) */

/*     LDE     INTEGER */
/*             The leading dimension of the array E.  LDE >= MAX(1,N). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             On entry, if FACT = 'F', then the leading N-by-N part of */
/*             this array must contain the orthogonal matrix Q from */
/*             the generalized Schur factorization (see definitions (3) */
/*             and (4) in section METHOD). */
/*             If FACT = 'N', Q need not be set on entry. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the orthogonal matrix Q from the generalized Schur */
/*             factorization. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q.  LDQ >= MAX(1,N). */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             On entry, if FACT = 'F', then the leading N-by-N part of */
/*             this array must contain the orthogonal matrix Z from */
/*             the generalized Schur factorization (see definitions (3) */
/*             and (4) in section METHOD). */
/*             If FACT = 'N', Z need not be set on entry. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the orthogonal matrix Z from the generalized Schur */
/*             factorization. */

/*     LDZ     INTEGER */
/*             The leading dimension of the array Z.  LDZ >= MAX(1,N). */

/*     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             On entry, if JOB = 'B' or 'X', then the leading N-by-N */
/*             part of this array must contain the right hand side matrix */
/*             Y of the equation. Either the lower or the upper */
/*             triangular part of this array is needed (see mode */
/*             parameter UPLO). */
/*             If JOB = 'S', X is not referenced. */
/*             On exit, if JOB = 'B' or 'X', and INFO = 0, 3, or 4, then */
/*             the leading N-by-N part of this array contains the */
/*             solution matrix X of the equation. */
/*             If JOB = 'S', X is not referenced. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X.  LDX >= MAX(1,N). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor set to avoid overflow in X. */
/*             (0 < SCALE <= 1) */

/*     SEP     (output) DOUBLE PRECISION */
/*             If JOB = 'S' or JOB = 'B', and INFO = 0, 3, or 4, then */
/*             SEP contains an estimate of the separation of the */
/*             Lyapunov operator. */

/*     FERR    (output) DOUBLE PRECISION */
/*             If JOB = 'B', and INFO = 0, 3, or 4, then FERR contains an */
/*             estimated forward error bound for the solution X. If XTRUE */
/*             is the true solution, FERR estimates the relative error */
/*             in the computed solution, measured in the Frobenius norm: */
/*             norm(X - XTRUE) / norm(XTRUE) */

/*     ALPHAR  (output) DOUBLE PRECISION array, dimension (N) */
/*     ALPHAI  (output) DOUBLE PRECISION array, dimension (N) */
/*     BETA    (output) DOUBLE PRECISION array, dimension (N) */
/*             If FACT = 'N' and INFO = 0, 3, or 4, then */
/*             (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, are the */
/*             eigenvalues of the matrix pencil A - lambda * E. */
/*             If FACT = 'F', ALPHAR, ALPHAI, and BETA are not */
/*             referenced. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N**2) */
/*             IWORK is not referenced if JOB = 'X'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. The following table */
/*             contains the minimal work space requirements depending */
/*             on the choice of JOB and FACT. */

/*                    JOB        FACT    |  LDWORK */
/*                    -------------------+------------------- */
/*                    'X'        'F'     |  MAX(1,N) */
/*                    'X'        'N'     |  MAX(1,4*N) */
/*                    'B', 'S'   'F'     |  MAX(1,2*N**2) */
/*                    'B', 'S'   'N'     |  MAX(1,2*N**2,4*N) */

/*             For optimum performance, LDWORK should be larger. */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  FACT = 'F' and the matrix contained in the upper */
/*                   Hessenberg part of the array A is not in upper */
/*                   quasitriangular form; */
/*             = 2:  FACT = 'N' and the pencil A - lambda * E cannot be */
/*                   reduced to generalized Schur form: LAPACK routine */
/*                   DGEGS has failed to converge; */
/*             = 3:  DICO = 'D' and the pencil A - lambda * E has a */
/*                   pair of reciprocal eigenvalues. That is, lambda_i = */
/*                   1/lambda_j for some i and j, where lambda_i and */
/*                   lambda_j are eigenvalues of A - lambda * E. Hence, */
/*                   equation (2) is singular;  perturbed values were */
/*                   used to solve the equation (but the matrices A and */
/*                   E are unchanged); */
/*             = 4:  DICO = 'C' and the pencil A - lambda * E has a */
/*                   degenerate pair of eigenvalues. That is, lambda_i = */
/*                   -lambda_j for some i and j, where lambda_i and */
/*                   lambda_j are eigenvalues of A - lambda * E. Hence, */
/*                   equation (1) is singular;  perturbed values were */
/*                   used to solve the equation (but the matrices A and */
/*                   E are unchanged). */

/*     METHOD */

/*     A straightforward generalization [3] of the method proposed by */
/*     Bartels and Stewart [1] is utilized to solve (1) or (2). */

/*     First the pencil A - lambda * E is reduced to real generalized */
/*     Schur form A_s - lambda * E_s by means of orthogonal */
/*     transformations (QZ-algorithm): */

/*        A_s = Q**T * A * Z   (upper quasitriangular)                (3) */

/*        E_s = Q**T * E * Z   (upper triangular).                    (4) */

/*     If FACT = 'F', this step is omitted. Assuming SCALE = 1 and */
/*     defining */

/*              ( Z**T * Y * Z   :   TRANS = 'N' */
/*        Y_s = < */
/*              ( Q**T * Y * Q   :   TRANS = 'T' */


/*              ( Q**T * X * Q    if TRANS = 'N' */
/*        X_s = <                                                     (5) */
/*              ( Z**T * X * Z    if TRANS = 'T' */

/*     leads to the reduced Lyapunov equation */

/*               T                      T */
/*        op(A_s)  X_s op(E_s) + op(E_s)  X_s op(A_s) = Y_s,          (6) */

/*     or */
/*               T                      T */
/*        op(A_s)  X_s op(A_s) - op(E_s)  X_s op(E_s) = Y_s,          (7) */

/*     which are equivalent to (1) or (2), respectively. The solution X_s */
/*     of (6) or (7) is computed via block back substitution (if TRANS = */
/*     'N') or block forward substitution (if TRANS = 'T'), where the */
/*     block order is at most 2. (See [1] and [3] for details.) */
/*     Equation (5) yields the solution matrix X. */

/*     For fast computation the estimates of the separation and the */
/*     forward error are gained from (6) or (7) rather than (1) or */
/*     (2), respectively. We consider (6) and (7) as special cases of the */
/*     generalized Sylvester equation */

/*        R * X * S + U * X * V = Y,                                  (8) */

/*     whose separation is defined as follows */

/*        sep = sep(R,S,U,V) =   min   || R * X * S + U * X * V || . */
/*                            ||X|| = 1                           F */
/*                                 F */

/*     Equation (8) is equivalent to the system of linear equations */

/*        K * vec(X) = (kron(S**T,R) + kron(V**T,U)) * vec(X) = vec(Y), */

/*     where kron is the Kronecker product of two matrices and vec */
/*     is the mapping that stacks the columns of a matrix. If K is */
/*     nonsingular then */

/*        sep = 1 / ||K**(-1)|| . */
/*                             2 */

/*     We estimate ||K**(-1)|| by a method devised by Higham [2]. Note */
/*     that this method yields an estimation for the 1-norm but we use it */
/*     as an approximation for the 2-norm. Estimates for the forward */
/*     error norm are provided by */

/*        FERR = 2 * EPS * ||A_s||  * ||E_s||  / sep */
/*                                F          F */

/*     in the continuous-time case (1) and */

/*        FERR = EPS * ( ||A_s|| **2 + ||E_s|| **2 ) / sep */
/*                              F             F */

/*     in the discrete-time case (2). */
/*     The reciprocal condition number, RCOND, of the Lyapunov equation */
/*     can be estimated by FERR/EPS. */

/*     REFERENCES */

/*     [1] Bartels, R.H., Stewart, G.W. */
/*         Solution of the equation A X + X B = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     [2] Higham, N.J. */
/*         FORTRAN codes for estimating the one-norm of a real or complex */
/*         matrix, with applications to condition estimation. */
/*         A.C.M. Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, 1988. */

/*     [3] Penzl, T. */
/*         Numerical solution of generalized Lyapunov equations. */
/*         Advances in Comp. Math., vol. 8, pp. 33-48, 1998. */

/*     NUMERICAL ASPECTS */

/*     The number of flops required by the routine is given by the */
/*     following table. Note that we count a single floating point */
/*     arithmetic operation as one flop. c is an integer number of modest */
/*     size (say 4 or 5). */

/*                   |  FACT = 'F'            FACT = 'N' */
/*        -----------+------------------------------------------ */
/*        JOB = 'B'  |  (26+8*c)/3 * N**3     (224+8*c)/3 * N**3 */
/*        JOB = 'S'  |  8*c/3 * N**3          (198+8*c)/3 * N**3 */
/*        JOB = 'X'  |  26/3 * N**3           224/3 * N**3 */

/*     The algorithm is backward stable if the eigenvalues of the pencil */
/*     A - lambda * E are real. Otherwise, linear systems of order at */
/*     most 4 are involved into the computation. These systems are solved */
/*     by Gauss elimination with complete pivoting. The loss of stability */
/*     of the Gauss elimination with complete pivoting is rarely */
/*     encountered in practice. */

/*     The Lyapunov equation may be very ill-conditioned. In particular, */
/*     if DICO = 'D' and the pencil A - lambda * E has a pair of almost */
/*     reciprocal eigenvalues, or DICO = 'C' and the pencil has an almost */
/*     degenerate pair of eigenvalues, then the Lyapunov equation will be */
/*     ill-conditioned. Perturbed values were used to solve the equation. */
/*     Ill-conditioning can be detected by a very small value of the */
/*     reciprocal condition number RCOND. */

/*     CONTRIBUTOR */

/*     T. Penzl, Technical University Chemnitz, Germany, Aug. 1998. */

/*     REVISIONS */

/*     Sep. 1998 (V. Sima). */
/*     Dec. 1998 (V. Sima). */

/*     KEYWORDS */

/*     Lyapunov equation */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Decode input parameters. */

#line 403 "SG03AD.f"
    /* Parameter adjustments */
#line 403 "SG03AD.f"
    a_dim1 = *lda;
#line 403 "SG03AD.f"
    a_offset = 1 + a_dim1;
#line 403 "SG03AD.f"
    a -= a_offset;
#line 403 "SG03AD.f"
    e_dim1 = *lde;
#line 403 "SG03AD.f"
    e_offset = 1 + e_dim1;
#line 403 "SG03AD.f"
    e -= e_offset;
#line 403 "SG03AD.f"
    q_dim1 = *ldq;
#line 403 "SG03AD.f"
    q_offset = 1 + q_dim1;
#line 403 "SG03AD.f"
    q -= q_offset;
#line 403 "SG03AD.f"
    z_dim1 = *ldz;
#line 403 "SG03AD.f"
    z_offset = 1 + z_dim1;
#line 403 "SG03AD.f"
    z__ -= z_offset;
#line 403 "SG03AD.f"
    x_dim1 = *ldx;
#line 403 "SG03AD.f"
    x_offset = 1 + x_dim1;
#line 403 "SG03AD.f"
    x -= x_offset;
#line 403 "SG03AD.f"
    --alphar;
#line 403 "SG03AD.f"
    --alphai;
#line 403 "SG03AD.f"
    --beta;
#line 403 "SG03AD.f"
    --iwork;
#line 403 "SG03AD.f"
    --dwork;
#line 403 "SG03AD.f"

#line 403 "SG03AD.f"
    /* Function Body */
#line 403 "SG03AD.f"
    isdisc = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 404 "SG03AD.f"
    wantx = lsame_(job, "X", (ftnlen)1, (ftnlen)1);
#line 405 "SG03AD.f"
    wantsp = lsame_(job, "S", (ftnlen)1, (ftnlen)1);
#line 406 "SG03AD.f"
    wantbh = lsame_(job, "B", (ftnlen)1, (ftnlen)1);
#line 407 "SG03AD.f"
    isfact = lsame_(fact, "F", (ftnlen)1, (ftnlen)1);
#line 408 "SG03AD.f"
    istran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 409 "SG03AD.f"
    isuppr = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

#line 413 "SG03AD.f"
    if (! (isdisc || lsame_(dico, "C", (ftnlen)1, (ftnlen)1))) {
#line 414 "SG03AD.f"
	*info = -1;
#line 415 "SG03AD.f"
    } else if (! (wantx || wantsp || wantbh)) {
#line 416 "SG03AD.f"
	*info = -2;
#line 417 "SG03AD.f"
    } else if (! (isfact || lsame_(fact, "N", (ftnlen)1, (ftnlen)1))) {
#line 418 "SG03AD.f"
	*info = -3;
#line 419 "SG03AD.f"
    } else if (! (istran || lsame_(trans, "N", (ftnlen)1, (ftnlen)1))) {
#line 420 "SG03AD.f"
	*info = -4;
#line 421 "SG03AD.f"
    } else if (! (isuppr || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 422 "SG03AD.f"
	*info = -5;
#line 423 "SG03AD.f"
    } else if (*n < 0) {
#line 424 "SG03AD.f"
	*info = -6;
#line 425 "SG03AD.f"
    } else if (*lda < max(1,*n)) {
#line 426 "SG03AD.f"
	*info = -8;
#line 427 "SG03AD.f"
    } else if (*lde < max(1,*n)) {
#line 428 "SG03AD.f"
	*info = -10;
#line 429 "SG03AD.f"
    } else if (*ldq < max(1,*n)) {
#line 430 "SG03AD.f"
	*info = -12;
#line 431 "SG03AD.f"
    } else if (*ldz < max(1,*n)) {
#line 432 "SG03AD.f"
	*info = -14;
#line 433 "SG03AD.f"
    } else if (*ldx < max(1,*n)) {
#line 434 "SG03AD.f"
	*info = -16;
#line 435 "SG03AD.f"
    } else {
#line 436 "SG03AD.f"
	*info = 0;
#line 437 "SG03AD.f"
    }
#line 438 "SG03AD.f"
    if (*info == 0) {

/*        Compute minimal workspace. */

#line 442 "SG03AD.f"
	if (wantx) {
#line 443 "SG03AD.f"
	    if (isfact) {
#line 444 "SG03AD.f"
		minwrk = max(*n,1);
#line 445 "SG03AD.f"
	    } else {
/* Computing MAX */
#line 446 "SG03AD.f"
		i__1 = *n << 2;
#line 446 "SG03AD.f"
		minwrk = max(i__1,1);
#line 447 "SG03AD.f"
	    }
#line 448 "SG03AD.f"
	} else {
#line 449 "SG03AD.f"
	    if (isfact) {
/* Computing MAX */
#line 450 "SG03AD.f"
		i__1 = (*n << 1) * *n;
#line 450 "SG03AD.f"
		minwrk = max(i__1,1);
#line 451 "SG03AD.f"
	    } else {
/* Computing MAX */
#line 452 "SG03AD.f"
		i__1 = (*n << 1) * *n, i__2 = *n << 2, i__1 = max(i__1,i__2);
#line 452 "SG03AD.f"
		minwrk = max(i__1,1);
#line 453 "SG03AD.f"
	    }
#line 454 "SG03AD.f"
	}
#line 455 "SG03AD.f"
	if (minwrk > *ldwork) {
#line 456 "SG03AD.f"
	    *info = -25;
#line 457 "SG03AD.f"
	}
#line 458 "SG03AD.f"
    }
#line 459 "SG03AD.f"
    if (*info != 0) {
#line 460 "SG03AD.f"
	i__1 = -(*info);
#line 460 "SG03AD.f"
	xerbla_("SG03AD", &i__1, (ftnlen)6);
#line 461 "SG03AD.f"
	return 0;
#line 462 "SG03AD.f"
    }

/*     Quick return if possible. */

#line 466 "SG03AD.f"
    if (*n == 0) {
#line 467 "SG03AD.f"
	*scale = 1.;
#line 468 "SG03AD.f"
	if (! wantx) {
#line 468 "SG03AD.f"
	    *sep = 0.;
#line 468 "SG03AD.f"
	}
#line 469 "SG03AD.f"
	if (wantbh) {
#line 469 "SG03AD.f"
	    *ferr = 0.;
#line 469 "SG03AD.f"
	}
#line 470 "SG03AD.f"
	dwork[1] = 1.;
#line 471 "SG03AD.f"
	return 0;
#line 472 "SG03AD.f"
    }

#line 474 "SG03AD.f"
    if (isfact) {

/*        Make sure the upper Hessenberg part of A is quasitriangular. */

#line 478 "SG03AD.f"
	i__1 = *n - 2;
#line 478 "SG03AD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 479 "SG03AD.f"
	    if (a[i__ + 1 + i__ * a_dim1] != 0. && a[i__ + 2 + (i__ + 1) * 
		    a_dim1] != 0.) {
#line 480 "SG03AD.f"
		*info = 1;
#line 481 "SG03AD.f"
		return 0;
#line 482 "SG03AD.f"
	    }
#line 483 "SG03AD.f"
/* L20: */
#line 483 "SG03AD.f"
	}
#line 484 "SG03AD.f"
    }

#line 486 "SG03AD.f"
    if (! isfact) {

/*        Reduce A - lambda * E to generalized Schur form. */

/*           A := Q**T * A * Z   (upper quasitriangular) */
/*           E := Q**T * E * Z   (upper triangular) */

/*        ( Workspace: >= MAX(1,4*N) ) */

#line 495 "SG03AD.f"
	dgegs_("Vectors", "Vectors", n, &a[a_offset], lda, &e[e_offset], lde, 
		&alphar[1], &alphai[1], &beta[1], &q[q_offset], ldq, &z__[
		z_offset], ldz, &dwork[1], ldwork, &info1, (ftnlen)7, (ftnlen)
		7);
#line 498 "SG03AD.f"
	if (info1 != 0) {
#line 499 "SG03AD.f"
	    *info = 2;
#line 500 "SG03AD.f"
	    return 0;
#line 501 "SG03AD.f"
	}
#line 502 "SG03AD.f"
	optwrk = (integer) dwork[1];
#line 503 "SG03AD.f"
    } else {
#line 504 "SG03AD.f"
	optwrk = minwrk;
#line 505 "SG03AD.f"
    }

#line 507 "SG03AD.f"
    if (wantbh || wantx) {

/*        Transform right hand side. */

/*           X := Z**T * X * Z  or  X := Q**T * X * Q */

/*        Use BLAS 3 if there is enough workspace. Otherwise, use BLAS 2. */

/*        ( Workspace: >= N ) */

#line 517 "SG03AD.f"
	if (*ldwork < *n * *n) {
#line 518 "SG03AD.f"
	    if (istran) {
#line 519 "SG03AD.f"
		mb01rw_(uplo, "Transpose", n, n, &x[x_offset], ldx, &q[
			q_offset], ldq, &dwork[1], &info1, (ftnlen)1, (ftnlen)
			9);
#line 521 "SG03AD.f"
	    } else {
#line 522 "SG03AD.f"
		mb01rw_(uplo, "Transpose", n, n, &x[x_offset], ldx, &z__[
			z_offset], ldz, &dwork[1], &info1, (ftnlen)1, (ftnlen)
			9);
#line 524 "SG03AD.f"
	    }
#line 525 "SG03AD.f"
	} else {
#line 526 "SG03AD.f"
	    if (istran) {
#line 527 "SG03AD.f"
		mb01rd_(uplo, "Transpose", n, n, &c_b20, &c_b21, &x[x_offset],
			 ldx, &q[q_offset], ldq, &x[x_offset], ldx, &dwork[1],
			 ldwork, info, (ftnlen)1, (ftnlen)9);
#line 529 "SG03AD.f"
	    } else {
#line 530 "SG03AD.f"
		mb01rd_(uplo, "Transpose", n, n, &c_b20, &c_b21, &x[x_offset],
			 ldx, &z__[z_offset], ldz, &x[x_offset], ldx, &dwork[
			1], ldwork, info, (ftnlen)1, (ftnlen)9);
#line 532 "SG03AD.f"
	    }
#line 533 "SG03AD.f"
	}
#line 534 "SG03AD.f"
	if (! isuppr) {
#line 535 "SG03AD.f"
	    i__1 = *n - 1;
#line 535 "SG03AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 536 "SG03AD.f"
		i__2 = *n - i__;
#line 536 "SG03AD.f"
		dcopy_(&i__2, &x[i__ + 1 + i__ * x_dim1], &c__1, &x[i__ + (
			i__ + 1) * x_dim1], ldx);
#line 537 "SG03AD.f"
/* L40: */
#line 537 "SG03AD.f"
	    }
#line 538 "SG03AD.f"
	}
/* Computing MAX */
#line 539 "SG03AD.f"
	i__1 = optwrk, i__2 = *n * *n;
#line 539 "SG03AD.f"
	optwrk = max(i__1,i__2);

/*        Solve reduced generalized Lyapunov equation. */

#line 543 "SG03AD.f"
	if (isdisc) {
#line 544 "SG03AD.f"
	    sg03ax_(trans, n, &a[a_offset], lda, &e[e_offset], lde, &x[
		    x_offset], ldx, scale, &info1, (ftnlen)1);
#line 545 "SG03AD.f"
	    if (info1 != 0) {
#line 545 "SG03AD.f"
		*info = 3;
#line 545 "SG03AD.f"
	    }
#line 547 "SG03AD.f"
	} else {
#line 548 "SG03AD.f"
	    sg03ay_(trans, n, &a[a_offset], lda, &e[e_offset], lde, &x[
		    x_offset], ldx, scale, &info1, (ftnlen)1);
#line 549 "SG03AD.f"
	    if (info1 != 0) {
#line 549 "SG03AD.f"
		*info = 4;
#line 549 "SG03AD.f"
	    }
#line 551 "SG03AD.f"
	}

/*        Transform the solution matrix back. */

/*           X := Q * X * Q**T  or  X := Z * X * Z**T. */

/*        Use BLAS 3 if there is enough workspace. Otherwise, use BLAS 2. */

/*        ( Workspace: >= N ) */

#line 561 "SG03AD.f"
	if (*ldwork < *n * *n) {
#line 562 "SG03AD.f"
	    if (istran) {
#line 563 "SG03AD.f"
		mb01rw_("Upper", "NoTranspose", n, n, &x[x_offset], ldx, &z__[
			z_offset], ldz, &dwork[1], &info1, (ftnlen)5, (ftnlen)
			11);
#line 565 "SG03AD.f"
	    } else {
#line 566 "SG03AD.f"
		mb01rw_("Upper", "NoTranspose", n, n, &x[x_offset], ldx, &q[
			q_offset], ldq, &dwork[1], &info1, (ftnlen)5, (ftnlen)
			11);
#line 568 "SG03AD.f"
	    }
#line 569 "SG03AD.f"
	} else {
#line 570 "SG03AD.f"
	    if (istran) {
#line 571 "SG03AD.f"
		mb01rd_("Upper", "NoTranspose", n, n, &c_b20, &c_b21, &x[
			x_offset], ldx, &z__[z_offset], ldz, &x[x_offset], 
			ldx, &dwork[1], ldwork, info, (ftnlen)5, (ftnlen)11);
#line 573 "SG03AD.f"
	    } else {
#line 574 "SG03AD.f"
		mb01rd_("Upper", "NoTranspose", n, n, &c_b20, &c_b21, &x[
			x_offset], ldx, &q[q_offset], ldq, &x[x_offset], ldx, 
			&dwork[1], ldwork, info, (ftnlen)5, (ftnlen)11);
#line 576 "SG03AD.f"
	    }
#line 577 "SG03AD.f"
	}
#line 578 "SG03AD.f"
	i__1 = *n - 1;
#line 578 "SG03AD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 579 "SG03AD.f"
	    i__2 = *n - i__;
#line 579 "SG03AD.f"
	    dcopy_(&i__2, &x[i__ + (i__ + 1) * x_dim1], ldx, &x[i__ + 1 + i__ 
		    * x_dim1], &c__1);
#line 580 "SG03AD.f"
/* L60: */
#line 580 "SG03AD.f"
	}
#line 581 "SG03AD.f"
    }

#line 583 "SG03AD.f"
    if (wantbh || wantsp) {

/*        Estimate the 1-norm of the inverse Kronecker product matrix */
/*        belonging to the reduced generalized Lyapunov equation. */

/*        ( Workspace: 2*N*N ) */

#line 590 "SG03AD.f"
	est = 0.;
#line 591 "SG03AD.f"
	kase = 0;
#line 592 "SG03AD.f"
L80:
#line 593 "SG03AD.f"
	i__1 = *n * *n;
#line 593 "SG03AD.f"
	dlacon_(&i__1, &dwork[*n * *n + 1], &dwork[1], &iwork[1], &est, &kase)
		;
#line 594 "SG03AD.f"
	if (kase != 0) {
#line 595 "SG03AD.f"
	    if (kase == 1 && ! istran || kase != 1 && istran) {
#line 597 "SG03AD.f"
		*(unsigned char *)etrans = 'N';
#line 598 "SG03AD.f"
	    } else {
#line 599 "SG03AD.f"
		*(unsigned char *)etrans = 'T';
#line 600 "SG03AD.f"
	    }
#line 601 "SG03AD.f"
	    if (isdisc) {
#line 602 "SG03AD.f"
		sg03ax_(etrans, n, &a[a_offset], lda, &e[e_offset], lde, &
			dwork[1], n, &scale1, &info1, (ftnlen)1);
#line 604 "SG03AD.f"
		if (info1 != 0) {
#line 604 "SG03AD.f"
		    *info = 3;
#line 604 "SG03AD.f"
		}
#line 606 "SG03AD.f"
	    } else {
#line 607 "SG03AD.f"
		sg03ay_(etrans, n, &a[a_offset], lda, &e[e_offset], lde, &
			dwork[1], n, &scale1, &info1, (ftnlen)1);
#line 609 "SG03AD.f"
		if (info1 != 0) {
#line 609 "SG03AD.f"
		    *info = 4;
#line 609 "SG03AD.f"
		}
#line 611 "SG03AD.f"
	    }
#line 612 "SG03AD.f"
	    goto L80;
#line 613 "SG03AD.f"
	}
#line 614 "SG03AD.f"
	*sep = scale1 / est;
#line 615 "SG03AD.f"
    }

/*     Estimate the relative forward error. */

/*     ( Workspace: 2*N ) */

#line 621 "SG03AD.f"
    if (wantbh) {
#line 622 "SG03AD.f"
	eps = dlamch_("Precision", (ftnlen)9);
#line 623 "SG03AD.f"
	i__1 = *n;
#line 623 "SG03AD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
#line 624 "SG03AD.f"
	    i__3 = i__ + 1;
#line 624 "SG03AD.f"
	    i__2 = min(i__3,*n);
#line 624 "SG03AD.f"
	    dwork[i__] = dnrm2_(&i__2, &a[i__ * a_dim1 + 1], &c__1);
#line 625 "SG03AD.f"
	    dwork[*n + i__] = dnrm2_(&i__, &e[i__ * e_dim1 + 1], &c__1);
#line 626 "SG03AD.f"
/* L100: */
#line 626 "SG03AD.f"
	}
#line 627 "SG03AD.f"
	norma = dnrm2_(n, &dwork[1], &c__1);
#line 628 "SG03AD.f"
	norme = dnrm2_(n, &dwork[*n + 1], &c__1);
#line 629 "SG03AD.f"
	if (isdisc) {
/* Computing 2nd power */
#line 630 "SG03AD.f"
	    d__1 = norma;
/* Computing 2nd power */
#line 630 "SG03AD.f"
	    d__2 = norme;
#line 630 "SG03AD.f"
	    *ferr = (d__1 * d__1 + d__2 * d__2) * eps / *sep;
#line 631 "SG03AD.f"
	} else {
#line 632 "SG03AD.f"
	    *ferr = norma * 2. * norme * eps / *sep;
#line 633 "SG03AD.f"
	}
#line 634 "SG03AD.f"
    }

#line 636 "SG03AD.f"
    dwork[1] = (doublereal) max(optwrk,minwrk);
#line 637 "SG03AD.f"
    return 0;
/* *** Last line of SG03AD *** */
} /* sg03ad_ */

