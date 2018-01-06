#line 1 "NF01BP.f"
/* NF01BP.f -- translated by f2c (version 20100827).
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

#line 1 "NF01BP.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b53 = 1.;

/* Subroutine */ int nf01bp_(char *cond, integer *n, integer *ipar, integer *
	lipar, doublereal *r__, integer *ldr, integer *ipvt, doublereal *diag,
	 doublereal *qtb, doublereal *delta, doublereal *par, integer *ranks, 
	doublereal *x, doublereal *rx, doublereal *tol, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen cond_len)
{
    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, l, n2, bn;
    static doublereal fp;
    static integer jw, st, bsm, bsn, lds;
    static doublereal sum, parc;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ibsn, rank;
    static doublereal parl;
    static logical sing;
    static integer iter;
    static doublereal temp, paru;
    static integer nths;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static logical badrk;
    extern /* Subroutine */ int nf01bq_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen);
    static logical econd;
    extern /* Subroutine */ int nf01br_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen), 
	    md03by_(char *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , integer *, ftnlen);
    static char condl[1];
    static logical ncond;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal dwarf;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal dmino;
    static logical ucond;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal gnorm;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal dxnorm;


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

/*     To determine a value for the Levenberg-Marquardt parameter PAR */
/*     such that if x solves the system */

/*           J*x = b ,     sqrt(PAR)*D*x = 0 , */

/*     in the least squares sense, where J is an m-by-n matrix, D is an */
/*     n-by-n nonsingular diagonal matrix, and b is an m-vector, and if */
/*     DELTA is a positive number, DXNORM is the Euclidean norm of D*x, */
/*     then either PAR is zero and */

/*           ( DXNORM - DELTA ) .LE. 0.1*DELTA , */

/*     or PAR is positive and */

/*           ABS( DXNORM - DELTA ) .LE. 0.1*DELTA . */

/*     The matrix J is the current Jacobian matrix of a nonlinear least */
/*     squares problem, provided in a compressed form by SLICOT Library */
/*     routine NF01BD. It is assumed that a block QR factorization, with */
/*     column pivoting, of J is available, that is, J*P = Q*R, where P is */
/*     a permutation matrix, Q has orthogonal columns, and R is an upper */
/*     triangular matrix with diagonal elements of nonincreasing */
/*     magnitude for each block, as returned by SLICOT Library */
/*     routine NF01BS. The routine NF01BP needs the upper triangle of R */
/*     in compressed form, the permutation matrix P, and the first */
/*     n components of Q'*b (' denotes the transpose). On output, */
/*     NF01BP also provides a compressed representation of an upper */
/*     triangular matrix S, such that */

/*           P'*(J'*J + PAR*D*D)*P = S'*S . */

/*     Matrix S is used in the solution process. The matrix R has the */
/*     following structure */

/*         /   R_1    0    ..   0   |   L_1   \ */
/*         |    0    R_2   ..   0   |   L_2   | */
/*         |    :     :    ..   :   |    :    | , */
/*         |    0     0    ..  R_l  |   L_l   | */
/*         \    0     0    ..   0   |  R_l+1  / */

/*     where the submatrices R_k, k = 1:l, have the same order BSN, */
/*     and R_k, k = 1:l+1, are square and upper triangular. This matrix */
/*     is stored in the compressed form */

/*              /   R_1  |   L_1   \ */
/*              |   R_2  |   L_2   | */
/*       Rc =   |    :   |    :    | , */
/*              |   R_l  |   L_l   | */
/*              \    X   |  R_l+1  / */

/*     where the submatrix X is irrelevant. The matrix S has the same */
/*     structure as R, and its diagonal blocks are denoted by S_k, */
/*     k = 1:l+1. */

/*     If l <= 1, then the full upper triangle of the matrix R is stored. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     COND    CHARACTER*1 */
/*             Specifies whether the condition of the diagonal blocks R_k */
/*             and S_k of the matrices R and S should be estimated, */
/*             as follows: */
/*             = 'E' :  use incremental condition estimation for each */
/*                      diagonal block of R_k and S_k to find its */
/*                      numerical rank; */
/*             = 'N' :  do not use condition estimation, but check the */
/*                      diagonal entries of R_k and S_k for zero values; */
/*             = 'U' :  use the ranks already stored in RANKS (for R). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix R.  N = BN*BSN + ST >= 0. */
/*             (See parameter description below.) */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             The integer parameters describing the structure of the */
/*             matrix R, as follows: */
/*             IPAR(1) must contain ST, the number of columns of the */
/*                     submatrices L_k and the order of R_l+1.  ST >= 0. */
/*             IPAR(2) must contain BN, the number of blocks, l, in the */
/*                     block diagonal part of R.  BN >= 0. */
/*             IPAR(3) must contain BSM, the number of rows of the blocks */
/*                     R_k, k = 1:l.  BSM >= 0. */
/*             IPAR(4) must contain BSN, the number of columns of the */
/*                     blocks R_k, k = 1:l.  BSN >= 0. */
/*             BSM is not used by this routine, but assumed equal to BSN. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 4. */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR, NC) */
/*             where NC = N if BN <= 1, and NC = BSN+ST, if BN > 1. */
/*             On entry, the leading N-by-NC part of this array must */
/*             contain the (compressed) representation (Rc) of the upper */
/*             triangular matrix R. If BN > 1, the submatrix X in Rc is */
/*             not referenced. The zero strict lower triangles of R_k, */
/*             k = 1:l+1, need not be set. If BN <= 1 or BSN = 0, then */
/*             the full upper triangle of R must be stored. */
/*             On exit, the full upper triangles of R_k, k = 1:l+1, and */
/*             L_k, k = 1:l, are unaltered, and the strict lower */
/*             triangles of R_k, k = 1:l+1, contain the corresponding */
/*             strict upper triangles (transposed) of the upper */
/*             triangular matrix S. */
/*             If BN <= 1 or BSN = 0, then the transpose of the strict */
/*             upper triangle of S is stored in the strict lower triangle */
/*             of R. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,N). */

/*     IPVT    (input) INTEGER array, dimension (N) */
/*             This array must define the permutation matrix P such that */
/*             J*P = Q*R. Column j of P is column IPVT(j) of the identity */
/*             matrix. */

/*     DIAG    (input) DOUBLE PRECISION array, dimension (N) */
/*             This array must contain the diagonal elements of the */
/*             matrix D.  DIAG(I) <> 0, I = 1,...,N. */

/*     QTB     (input) DOUBLE PRECISION array, dimension (N) */
/*             This array must contain the first n elements of the */
/*             vector Q'*b. */

/*     DELTA   (input) DOUBLE PRECISION */
/*             An upper bound on the Euclidean norm of D*x.  DELTA > 0. */

/*     PAR     (input/output) DOUBLE PRECISION */
/*             On entry, PAR must contain an initial estimate of the */
/*             Levenberg-Marquardt parameter.  PAR >= 0. */
/*             On exit, it contains the final estimate of this parameter. */

/*     RANKS   (input or output) INTEGER array, dimension (r), where */
/*             r = BN + 1,  if ST > 0, BSN > 0, and BN > 1; */
/*             r = BN,      if ST = 0 and BSN > 0; */
/*             r = 1,       if ST > 0 and ( BSN = 0 or BN <= 1 ); */
/*             r = 0,       if ST = 0 and BSN = 0. */
/*             On entry, if COND = 'U' and N > 0, this array must contain */
/*             the numerical ranks of the submatrices R_k, k = 1:l(+1). */
/*             On exit, if N > 0, this array contains the numerical ranks */
/*             of the submatrices S_k, k = 1:l(+1). */

/*     X       (output) DOUBLE PRECISION array, dimension (N) */
/*             This array contains the least squares solution of the */
/*             system J*x = b, sqrt(PAR)*D*x = 0. */

/*     RX      (output) DOUBLE PRECISION array, dimension (N) */
/*             This array contains the matrix-vector product -R*P'*x. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If COND = 'E', the tolerance to be used for finding the */
/*             ranks of the submatrices R_k and S_k. If the user sets */
/*             TOL > 0, then the given value of TOL is used as a lower */
/*             bound for the reciprocal condition number;  a (sub)matrix */
/*             whose estimated condition number is less than 1/TOL is */
/*             considered to be of full rank.  If the user sets TOL <= 0, */
/*             then an implicitly computed, default tolerance, defined by */
/*             TOLDEF = N*EPS,  is used instead, where EPS is the machine */
/*             precision (see LAPACK Library routine DLAMCH). */
/*             This parameter is not relevant if COND = 'U' or 'N'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, the first N elements of this array contain the */
/*             diagonal elements of the upper triangular matrix S. */
/*             If BN > 1 and BSN > 0, the elements N+1 : N+ST*(N-ST) */
/*             contain the submatrix (S(1:N-ST,N-ST+1:N))' of the */
/*             matrix S. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 2*N,              if BN <= 1 or  BSN = 0 and */
/*                                                        COND <> 'E'; */
/*             LDWORK >= 4*N,              if BN <= 1 or  BSN = 0 and */
/*                                                        COND =  'E'; */
/*             LDWORK >= ST*(N-ST) + 2*N,  if BN >  1 and BSN > 0 and */
/*                                                        COND <> 'E'; */
/*             LDWORK >= ST*(N-ST) + 2*N + 2*MAX(BSN,ST), */
/*                                         if BN >  1 and BSN > 0 and */
/*                                                        COND =  'E'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The algorithm computes the Gauss-Newton direction. An approximate */
/*     basic least squares solution is found if the Jacobian is rank */
/*     deficient. The computations exploit the special structure and */
/*     storage scheme of the matrix R. If one or more of the submatrices */
/*     R_k or S_k, k = 1:l+1, is singular, then the computed result is */
/*     not the basic least squares solution for the whole problem, but a */
/*     concatenation of (least squares) solutions of the individual */
/*     subproblems involving R_k or S_k, k = 1:l+1 (with adapted right */
/*     hand sides). */

/*     If the Gauss-Newton direction is not acceptable, then an iterative */
/*     algorithm obtains improved lower and upper bounds for the */
/*     Levenberg-Marquardt parameter PAR. Only a few iterations are */
/*     generally needed for convergence of the algorithm. If, however, */
/*     the limit of ITMAX = 10 iterations is reached, then the output PAR */
/*     will contain the best value obtained so far. If the Gauss-Newton */
/*     step is acceptable, it is stored in x, and PAR is set to zero, */
/*     hence S = R. */

/*     REFERENCES */

/*     [1] More, J.J., Garbow, B.S, and Hillstrom, K.E. */
/*         User's Guide for MINPACK-1. */
/*         Applied Math. Division, Argonne National Laboratory, Argonne, */
/*         Illinois, Report ANL-80-74, 1980. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0(N*(BSN+ST)) operations and is backward */
/*     stable, if R is nonsingular. */

/*     FURTHER COMMENTS */

/*     This routine is a structure-exploiting, LAPACK-based modification */
/*     of LMPAR from the MINPACK package [1], and with optional condition */
/*     estimation. The option COND = 'U' is useful when dealing with */
/*     several right-hand side vectors, but RANKS array should be reset. */
/*     If COND = 'E', but the matrix S is guaranteed to be nonsingular */
/*     and well conditioned relative to TOL, i.e., rank(R) = N, and */
/*     min(DIAG) > 0, then its condition is not estimated. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     V. Sima, Feb. 2004. */

/*     KEYWORDS */

/*     Linear system of equations, matrix operations, plane rotations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

#line 309 "NF01BP.f"
    /* Parameter adjustments */
#line 309 "NF01BP.f"
    --ipar;
#line 309 "NF01BP.f"
    r_dim1 = *ldr;
#line 309 "NF01BP.f"
    r_offset = 1 + r_dim1;
#line 309 "NF01BP.f"
    r__ -= r_offset;
#line 309 "NF01BP.f"
    --ipvt;
#line 309 "NF01BP.f"
    --diag;
#line 309 "NF01BP.f"
    --qtb;
#line 309 "NF01BP.f"
    --ranks;
#line 309 "NF01BP.f"
    --x;
#line 309 "NF01BP.f"
    --rx;
#line 309 "NF01BP.f"
    --dwork;
#line 309 "NF01BP.f"

#line 309 "NF01BP.f"
    /* Function Body */
#line 309 "NF01BP.f"
    econd = lsame_(cond, "E", (ftnlen)1, (ftnlen)1);
#line 310 "NF01BP.f"
    ncond = lsame_(cond, "N", (ftnlen)1, (ftnlen)1);
#line 311 "NF01BP.f"
    ucond = lsame_(cond, "U", (ftnlen)1, (ftnlen)1);
#line 312 "NF01BP.f"
    *info = 0;
#line 313 "NF01BP.f"
    n2 = *n << 1;
#line 314 "NF01BP.f"
    if (! (econd || ncond || ucond)) {
#line 315 "NF01BP.f"
	*info = -1;
#line 316 "NF01BP.f"
    } else if (*n < 0) {
#line 317 "NF01BP.f"
	*info = -2;
#line 318 "NF01BP.f"
    } else if (*lipar < 4) {
#line 319 "NF01BP.f"
	*info = -4;
#line 320 "NF01BP.f"
    } else if (*ldr < max(1,*n)) {
#line 321 "NF01BP.f"
	*info = -6;
#line 322 "NF01BP.f"
    } else if (*delta <= 0.) {
#line 323 "NF01BP.f"
	*info = -10;
#line 324 "NF01BP.f"
    } else if (*par < 0.) {
#line 325 "NF01BP.f"
	*info = -11;
#line 326 "NF01BP.f"
    } else {
#line 327 "NF01BP.f"
	st = ipar[1];
#line 328 "NF01BP.f"
	bn = ipar[2];
#line 329 "NF01BP.f"
	bsm = ipar[3];
#line 330 "NF01BP.f"
	bsn = ipar[4];
#line 331 "NF01BP.f"
	nths = bn * bsn;
/* Computing MIN */
#line 332 "NF01BP.f"
	i__1 = min(st,bn), i__1 = min(i__1,bsm);
#line 332 "NF01BP.f"
	if (min(i__1,bsn) < 0) {
#line 333 "NF01BP.f"
	    *info = -3;
#line 334 "NF01BP.f"
	} else if (*n != nths + st) {
#line 335 "NF01BP.f"
	    *info = -2;
#line 336 "NF01BP.f"
	} else {
#line 337 "NF01BP.f"
	    if (*n > 0) {
#line 337 "NF01BP.f"
		dmino = diag[1];
#line 337 "NF01BP.f"
	    }
#line 339 "NF01BP.f"
	    sing = FALSE_;

#line 341 "NF01BP.f"
	    i__1 = *n;
#line 341 "NF01BP.f"
	    for (j = 1; j <= i__1; ++j) {
#line 342 "NF01BP.f"
		if (diag[j] < dmino) {
#line 342 "NF01BP.f"
		    dmino = diag[j];
#line 342 "NF01BP.f"
		}
#line 344 "NF01BP.f"
		sing = sing || diag[j] == 0.;
#line 345 "NF01BP.f"
/* L10: */
#line 345 "NF01BP.f"
	    }

#line 347 "NF01BP.f"
	    if (sing) {
#line 348 "NF01BP.f"
		*info = -8;
#line 349 "NF01BP.f"
	    } else if (ucond) {
#line 350 "NF01BP.f"
		badrk = FALSE_;
#line 351 "NF01BP.f"
		if (bn <= 1 || bsn == 0) {
#line 352 "NF01BP.f"
		    if (*n > 0) {
#line 352 "NF01BP.f"
			badrk = ranks[1] < 0 || ranks[1] > *n;
#line 352 "NF01BP.f"
		    }
#line 354 "NF01BP.f"
		} else {
#line 355 "NF01BP.f"
		    rank = 0;

#line 357 "NF01BP.f"
		    i__1 = bn;
#line 357 "NF01BP.f"
		    for (k = 1; k <= i__1; ++k) {
#line 358 "NF01BP.f"
			badrk = badrk || ranks[k] < 0 || ranks[k] > bsn;
#line 360 "NF01BP.f"
			rank += ranks[k];
#line 361 "NF01BP.f"
/* L20: */
#line 361 "NF01BP.f"
		    }

#line 363 "NF01BP.f"
		    if (st > 0) {
#line 364 "NF01BP.f"
			badrk = badrk || ranks[bn + 1] < 0 || ranks[bn + 1] > 
				st;
#line 366 "NF01BP.f"
			rank += ranks[bn + 1];
#line 367 "NF01BP.f"
		    }
#line 368 "NF01BP.f"
		}
#line 369 "NF01BP.f"
		if (badrk) {
#line 369 "NF01BP.f"
		    *info = -12;
#line 369 "NF01BP.f"
		}
#line 371 "NF01BP.f"
	    } else {
#line 372 "NF01BP.f"
		jw = n2;
#line 373 "NF01BP.f"
		if (bn <= 1 || bsn == 0) {
#line 374 "NF01BP.f"
		    if (econd) {
#line 374 "NF01BP.f"
			jw = *n << 2;
#line 374 "NF01BP.f"
		    }
#line 376 "NF01BP.f"
		} else {
#line 377 "NF01BP.f"
		    jw = st * nths + jw;
#line 378 "NF01BP.f"
		    if (econd) {
#line 378 "NF01BP.f"
			jw = (max(bsn,st) << 1) + jw;
#line 378 "NF01BP.f"
		    }
#line 380 "NF01BP.f"
		}
#line 381 "NF01BP.f"
		if (*ldwork < jw) {
#line 381 "NF01BP.f"
		    *info = -17;
#line 381 "NF01BP.f"
		}
#line 383 "NF01BP.f"
	    }
#line 384 "NF01BP.f"
	}
#line 385 "NF01BP.f"
    }

/*     Return if there are illegal arguments. */

#line 389 "NF01BP.f"
    if (*info != 0) {
#line 390 "NF01BP.f"
	i__1 = -(*info);
#line 390 "NF01BP.f"
	xerbla_("NF01BP", &i__1, (ftnlen)6);
#line 391 "NF01BP.f"
	return 0;
#line 392 "NF01BP.f"
    }

/*     Quick return if possible. */

#line 396 "NF01BP.f"
    if (*n == 0) {
#line 397 "NF01BP.f"
	*par = 0.;
#line 398 "NF01BP.f"
	return 0;
#line 399 "NF01BP.f"
    }

#line 401 "NF01BP.f"
    if (bn <= 1 || bsn == 0) {

/*        Special case: R is just an upper triangular matrix. */
/*        Workspace: 4*N, if COND =  'E'; */
/*                   2*N, if COND <> 'E'. */

#line 407 "NF01BP.f"
	md03by_(cond, n, &r__[r_offset], ldr, &ipvt[1], &diag[1], &qtb[1], 
		delta, par, &ranks[1], &x[1], &rx[1], tol, &dwork[1], ldwork, 
		info, (ftnlen)1);
#line 409 "NF01BP.f"
	return 0;
#line 410 "NF01BP.f"
    }

/*     General case: l > 1 and BSN > 0. */
/*     DWARF is the smallest positive magnitude. */

#line 415 "NF01BP.f"
    dwarf = dlamch_("Underflow", (ftnlen)9);

/*     Compute and store in x the Gauss-Newton direction. If the */
/*     Jacobian is rank-deficient, obtain a least squares solution. */
/*     The array RX is used as workspace. */
/*     Workspace: 2*MAX(BSN,ST), if COND =  'E'; */
/*                0,             if COND <> 'E'. */

#line 423 "NF01BP.f"
    dcopy_(n, &qtb[1], &c__1, &rx[1], &c__1);
#line 424 "NF01BP.f"
    nf01br_(cond, "Upper", "No transpose", n, &ipar[1], lipar, &r__[r_offset],
	     ldr, &dwork[1], &dwork[1], &c__1, &rx[1], &ranks[1], tol, &dwork[
	    1], ldwork, info, (ftnlen)1, (ftnlen)5, (ftnlen)12);

#line 428 "NF01BP.f"
    i__1 = *n;
#line 428 "NF01BP.f"
    for (j = 1; j <= i__1; ++j) {
#line 429 "NF01BP.f"
	l = ipvt[j];
#line 430 "NF01BP.f"
	x[l] = rx[j];
#line 431 "NF01BP.f"
/* L30: */
#line 431 "NF01BP.f"
    }

/*     Initialize the iteration counter. */
/*     Evaluate the function at the origin, and test */
/*     for acceptance of the Gauss-Newton direction. */

#line 437 "NF01BP.f"
    iter = 0;

#line 439 "NF01BP.f"
    i__1 = *n;
#line 439 "NF01BP.f"
    for (j = 1; j <= i__1; ++j) {
#line 440 "NF01BP.f"
	dwork[j] = diag[j] * x[j];
#line 441 "NF01BP.f"
/* L40: */
#line 441 "NF01BP.f"
    }

#line 443 "NF01BP.f"
    dxnorm = dnrm2_(n, &dwork[1], &c__1);
#line 444 "NF01BP.f"
    fp = dxnorm - *delta;
#line 445 "NF01BP.f"
    if (fp > *delta * .1) {

/*        Set an appropriate option for estimating the condition of */
/*        the matrix S. */

#line 450 "NF01BP.f"
	lds = max(1,st);
#line 451 "NF01BP.f"
	jw = n2 + st * nths;
#line 452 "NF01BP.f"
	if (ucond) {
#line 453 "NF01BP.f"
	    if (*ldwork >= jw + (max(bsn,st) << 1)) {
#line 454 "NF01BP.f"
		*(unsigned char *)condl = 'E';
#line 455 "NF01BP.f"
		toldef = (doublereal) (*n) * dlamch_("Epsilon", (ftnlen)7);
#line 456 "NF01BP.f"
	    } else {
#line 457 "NF01BP.f"
		*(unsigned char *)condl = 'N';
#line 458 "NF01BP.f"
		toldef = *tol;
#line 459 "NF01BP.f"
	    }
#line 460 "NF01BP.f"
	} else {
#line 461 "NF01BP.f"
	    rank = 0;

#line 463 "NF01BP.f"
	    i__1 = bn;
#line 463 "NF01BP.f"
	    for (k = 1; k <= i__1; ++k) {
#line 464 "NF01BP.f"
		rank += ranks[k];
#line 465 "NF01BP.f"
/* L50: */
#line 465 "NF01BP.f"
	    }

#line 467 "NF01BP.f"
	    if (st > 0) {
#line 467 "NF01BP.f"
		rank += ranks[bn + 1];
#line 467 "NF01BP.f"
	    }
#line 469 "NF01BP.f"
	    *(unsigned char *)condl = *(unsigned char *)cond;
#line 470 "NF01BP.f"
	    toldef = *tol;
#line 471 "NF01BP.f"
	}

/*        If the Jacobian is not rank deficient, the Newton */
/*        step provides a lower bound, PARL, for the zero of */
/*        the function. Otherwise set this bound to zero. */

#line 477 "NF01BP.f"
	if (rank == *n) {

#line 479 "NF01BP.f"
	    i__1 = *n;
#line 479 "NF01BP.f"
	    for (j = 1; j <= i__1; ++j) {
#line 480 "NF01BP.f"
		l = ipvt[j];
#line 481 "NF01BP.f"
		rx[j] = diag[l] * (dwork[l] / dxnorm);
#line 482 "NF01BP.f"
/* L60: */
#line 482 "NF01BP.f"
	    }

#line 484 "NF01BP.f"
	    nf01br_("Use ranks", "Upper", "Transpose", n, &ipar[1], lipar, &
		    r__[r_offset], ldr, &dwork[1], &dwork[1], &c__1, &rx[1], &
		    ranks[1], tol, &dwork[1], ldwork, info, (ftnlen)9, (
		    ftnlen)5, (ftnlen)9);
#line 487 "NF01BP.f"
	    temp = dnrm2_(n, &rx[1], &c__1);
#line 488 "NF01BP.f"
	    parl = fp / *delta / temp / temp;

/*           For efficiency, use CONDL = 'U', if possible. */

#line 492 "NF01BP.f"
	    if (! lsame_(condl, "U", (ftnlen)1, (ftnlen)1) && dmino > 0.) {
#line 492 "NF01BP.f"
		*(unsigned char *)condl = 'U';
#line 492 "NF01BP.f"
	    }
#line 494 "NF01BP.f"
	} else {
#line 495 "NF01BP.f"
	    parl = 0.;
#line 496 "NF01BP.f"
	}

#line 498 "NF01BP.f"
	ibsn = 0;
#line 499 "NF01BP.f"
	k = 1;

/*        Calculate an upper bound, PARU, for the zero of the function. */

#line 503 "NF01BP.f"
	i__1 = *n;
#line 503 "NF01BP.f"
	for (j = 1; j <= i__1; ++j) {
#line 504 "NF01BP.f"
	    ++ibsn;
#line 505 "NF01BP.f"
	    if (j < nths) {
#line 506 "NF01BP.f"
		sum = ddot_(&ibsn, &r__[k + ibsn * r_dim1], &c__1, &qtb[k], &
			c__1);
#line 507 "NF01BP.f"
		if (ibsn == bsn) {
#line 508 "NF01BP.f"
		    ibsn = 0;
#line 509 "NF01BP.f"
		    k += bsn;
#line 510 "NF01BP.f"
		}
#line 511 "NF01BP.f"
	    } else if (j == nths) {
#line 512 "NF01BP.f"
		sum = ddot_(&ibsn, &r__[k + ibsn * r_dim1], &c__1, &qtb[k], &
			c__1);
#line 513 "NF01BP.f"
	    } else {
#line 514 "NF01BP.f"
		sum = ddot_(&j, &r__[ibsn * r_dim1 + 1], &c__1, &qtb[1], &
			c__1);
#line 515 "NF01BP.f"
	    }
#line 516 "NF01BP.f"
	    l = ipvt[j];
#line 517 "NF01BP.f"
	    rx[j] = sum / diag[l];
#line 518 "NF01BP.f"
/* L70: */
#line 518 "NF01BP.f"
	}

#line 520 "NF01BP.f"
	gnorm = dnrm2_(n, &rx[1], &c__1);
#line 521 "NF01BP.f"
	paru = gnorm / *delta;
#line 522 "NF01BP.f"
	if (paru == 0.) {
#line 522 "NF01BP.f"
	    paru = dwarf / min(*delta,.1) / .001;
#line 522 "NF01BP.f"
	}

/*        If the input PAR lies outside of the interval (PARL,PARU), */
/*        set PAR to the closer endpoint. */

#line 528 "NF01BP.f"
	*par = max(*par,parl);
#line 529 "NF01BP.f"
	*par = min(*par,paru);
#line 530 "NF01BP.f"
	if (*par == 0.) {
#line 530 "NF01BP.f"
	    *par = gnorm / dxnorm;
#line 530 "NF01BP.f"
	}

/*        Beginning of an iteration. */

#line 535 "NF01BP.f"
L80:
#line 536 "NF01BP.f"
	++iter;

/*           Evaluate the function at the current value of PAR. */

#line 540 "NF01BP.f"
	if (*par == 0.) {
/* Computing MAX */
#line 540 "NF01BP.f"
	    d__1 = dwarf, d__2 = paru * .001;
#line 540 "NF01BP.f"
	    *par = max(d__1,d__2);
#line 540 "NF01BP.f"
	}
#line 542 "NF01BP.f"
	temp = sqrt(*par);

#line 544 "NF01BP.f"
	i__1 = *n;
#line 544 "NF01BP.f"
	for (j = 1; j <= i__1; ++j) {
#line 545 "NF01BP.f"
	    rx[j] = temp * diag[j];
#line 546 "NF01BP.f"
/* L90: */
#line 546 "NF01BP.f"
	}

/*           Solve the system J*x = b , sqrt(PAR)*D*x = 0 , in a least */
/*           square sense. */
/*           The first N elements of DWORK contain the diagonal elements */
/*           of the upper triangular matrix S, and the next N elements */
/*           contain the the vector z, so that x = P*z (see NF01BQ). */
/*           The vector z is not preserved, to reduce the workspace. */
/*           The elements 2*N+1 : 2*N+ST*(N-ST) contain the */
/*           submatrix (S(1:N-ST,N-ST+1:N))' of the matrix S. */
/*           Workspace: ST*(N-ST) + 2*N,                 if CONDL <> 'E'; */
/*                      ST*(N-ST) + 2*N + 2*MAX(BSN,ST), if CONDL =  'E'. */

#line 559 "NF01BP.f"
	nf01bq_(condl, n, &ipar[1], lipar, &r__[r_offset], ldr, &ipvt[1], &rx[
		1], &qtb[1], &ranks[1], &x[1], &toldef, &dwork[1], ldwork, 
		info, (ftnlen)1);

#line 562 "NF01BP.f"
	i__1 = *n;
#line 562 "NF01BP.f"
	for (j = 1; j <= i__1; ++j) {
#line 563 "NF01BP.f"
	    dwork[*n + j] = diag[j] * x[j];
#line 564 "NF01BP.f"
/* L100: */
#line 564 "NF01BP.f"
	}

#line 566 "NF01BP.f"
	dxnorm = dnrm2_(n, &dwork[*n + 1], &c__1);
#line 567 "NF01BP.f"
	temp = fp;
#line 568 "NF01BP.f"
	fp = dxnorm - *delta;

/*           If the function is small enough, accept the current value */
/*           of PAR. Also test for the exceptional cases where PARL */
/*           is zero or the number of iterations has reached ITMAX. */

#line 574 "NF01BP.f"
	if (abs(fp) > *delta * .1 && (parl != 0. || fp > temp || temp >= 0.) 
		&& iter < 10) {

/*              Compute the Newton correction. */

#line 580 "NF01BP.f"
	    i__1 = *n;
#line 580 "NF01BP.f"
	    for (j = 1; j <= i__1; ++j) {
#line 581 "NF01BP.f"
		l = ipvt[j];
#line 582 "NF01BP.f"
		rx[j] = diag[l] * (dwork[*n + l] / dxnorm);
#line 583 "NF01BP.f"
/* L110: */
#line 583 "NF01BP.f"
	    }

#line 585 "NF01BP.f"
	    i__1 = *ldwork - jw;
#line 585 "NF01BP.f"
	    nf01br_("Use ranks", "Lower", "Transpose", n, &ipar[1], lipar, &
		    r__[r_offset], ldr, &dwork[1], &dwork[n2 + 1], &lds, &rx[
		    1], &ranks[1], tol, &dwork[jw], &i__1, info, (ftnlen)9, (
		    ftnlen)5, (ftnlen)9);
#line 588 "NF01BP.f"
	    temp = dnrm2_(n, &rx[1], &c__1);
#line 589 "NF01BP.f"
	    parc = fp / *delta / temp / temp;

/*              Depending on the sign of the function, update PARL */
/*              or PARU. */

#line 594 "NF01BP.f"
	    if (fp > 0.) {
#line 595 "NF01BP.f"
		parl = max(parl,*par);
#line 596 "NF01BP.f"
	    } else if (fp < 0.) {
#line 597 "NF01BP.f"
		paru = min(paru,*par);
#line 598 "NF01BP.f"
	    }

/*              Compute an improved estimate for PAR. */

/* Computing MAX */
#line 602 "NF01BP.f"
	    d__1 = parl, d__2 = *par + parc;
#line 602 "NF01BP.f"
	    *par = max(d__1,d__2);

/*              End of an iteration. */

#line 606 "NF01BP.f"
	    goto L80;
#line 607 "NF01BP.f"
	}
#line 608 "NF01BP.f"
    }

/*     Compute -R*P'*x = -R*z. */

#line 612 "NF01BP.f"
    i__1 = *n;
#line 612 "NF01BP.f"
    for (j = 1; j <= i__1; ++j) {
#line 613 "NF01BP.f"
	l = ipvt[j];
#line 614 "NF01BP.f"
	rx[j] = -x[l];
#line 615 "NF01BP.f"
/* L120: */
#line 615 "NF01BP.f"
    }

#line 617 "NF01BP.f"
    i__1 = nths;
#line 617 "NF01BP.f"
    i__2 = bsn;
#line 617 "NF01BP.f"
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 618 "NF01BP.f"
	dtrmv_("Upper", "NoTranspose", "NonUnit", &bsn, &r__[i__ + r_dim1], 
		ldr, &rx[i__], &c__1, (ftnlen)5, (ftnlen)11, (ftnlen)7);
#line 620 "NF01BP.f"
/* L130: */
#line 620 "NF01BP.f"
    }

#line 622 "NF01BP.f"
    if (st > 0) {
#line 623 "NF01BP.f"
	dgemv_("NoTranspose", &nths, &st, &c_b53, &r__[(bsn + 1) * r_dim1 + 1]
		, ldr, &rx[nths + 1], &c__1, &c_b53, &rx[1], &c__1, (ftnlen)
		11);
#line 625 "NF01BP.f"
	dtrmv_("Upper", "NoTranspose", "NonUnit", &st, &r__[nths + 1 + (bsn + 
		1) * r_dim1], ldr, &rx[nths + 1], &c__1, (ftnlen)5, (ftnlen)
		11, (ftnlen)7);
#line 627 "NF01BP.f"
    }

/*     Termination. If PAR = 0, set S. */

#line 631 "NF01BP.f"
    if (iter == 0) {
#line 632 "NF01BP.f"
	*par = 0.;
#line 633 "NF01BP.f"
	i__ = 1;

#line 635 "NF01BP.f"
	i__2 = bn;
#line 635 "NF01BP.f"
	for (k = 1; k <= i__2; ++k) {

#line 637 "NF01BP.f"
	    i__1 = bsn;
#line 637 "NF01BP.f"
	    for (j = 1; j <= i__1; ++j) {
#line 638 "NF01BP.f"
		dwork[i__] = r__[i__ + j * r_dim1];
#line 639 "NF01BP.f"
		i__3 = bsn - j + 1;
#line 639 "NF01BP.f"
		dcopy_(&i__3, &r__[i__ + j * r_dim1], ldr, &r__[i__ + j * 
			r_dim1], &c__1);
#line 640 "NF01BP.f"
		++i__;
#line 641 "NF01BP.f"
/* L140: */
#line 641 "NF01BP.f"
	    }

#line 643 "NF01BP.f"
/* L150: */
#line 643 "NF01BP.f"
	}

#line 645 "NF01BP.f"
	if (st > 0) {

#line 647 "NF01BP.f"
	    i__2 = bsn + st;
#line 647 "NF01BP.f"
	    for (j = bsn + 1; j <= i__2; ++j) {
#line 648 "NF01BP.f"
		dcopy_(&nths, &r__[j * r_dim1 + 1], &c__1, &dwork[*n + j - 
			bsn], &st);
#line 649 "NF01BP.f"
		dwork[i__] = r__[i__ + j * r_dim1];
#line 650 "NF01BP.f"
		i__1 = bsn + st - j + 1;
#line 650 "NF01BP.f"
		dcopy_(&i__1, &r__[i__ + j * r_dim1], ldr, &r__[i__ + j * 
			r_dim1], &c__1);
#line 651 "NF01BP.f"
		++i__;
#line 652 "NF01BP.f"
/* L160: */
#line 652 "NF01BP.f"
	    }

#line 654 "NF01BP.f"
	}
#line 655 "NF01BP.f"
    } else {

#line 657 "NF01BP.f"
	i__2 = *n + st * nths;
#line 657 "NF01BP.f"
	for (k = *n + 1; k <= i__2; ++k) {
#line 658 "NF01BP.f"
	    dwork[k] = dwork[k + *n];
#line 659 "NF01BP.f"
/* L170: */
#line 659 "NF01BP.f"
	}

#line 661 "NF01BP.f"
    }

#line 663 "NF01BP.f"
    return 0;

/* *** Last line of NF01BP *** */
} /* nf01bp_ */

