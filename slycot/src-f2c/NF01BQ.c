#line 1 "NF01BQ.f"
/* NF01BQ.f -- translated by f2c (version 20100827).
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

#line 1 "NF01BQ.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;

/* Subroutine */ int nf01bq_(char *cond, integer *n, integer *ipar, integer *
	lipar, doublereal *r__, integer *ldr, integer *ipvt, doublereal *diag,
	 doublereal *qtb, integer *ranks, doublereal *x, doublereal *tol, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen cond_len)
{
    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, l, ib, bn, kf, nc, is, jw, st, itc, bsm, bsn, 
	    itr, ibsn, nths;
    static logical econd;
    extern /* Subroutine */ int nf01br_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen), 
	    mb02yd_(char *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb04ow_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *), dcopy_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static doublereal qtbpj;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), xerbla_(char *, integer *, ftnlen);


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

/*     To determine a vector x which solves the system of linear */
/*     equations */

/*           J*x = b ,     D*x = 0 , */

/*     in the least squares sense, where J is an m-by-n matrix, */
/*     D is an n-by-n diagonal matrix, and b is an m-vector. The matrix J */
/*     is the current Jacobian of a nonlinear least squares problem, */
/*     provided in a compressed form by SLICOT Library routine NF01BD. */
/*     It is assumed that a block QR factorization, with column pivoting, */
/*     of J is available, that is, J*P = Q*R, where P is a permutation */
/*     matrix, Q has orthogonal columns, and R is an upper triangular */
/*     matrix with diagonal elements of nonincreasing magnitude for each */
/*     block, as returned by SLICOT Library routine NF01BS. The routine */
/*     NF01BQ needs the upper triangle of R in compressed form, the */
/*     permutation matrix P, and the first n components of Q'*b */
/*     (' denotes the transpose). The system J*x = b, D*x = 0, is then */
/*     equivalent to */

/*           R*z = Q'*b ,  P'*D*P*z = 0 ,                             (1) */

/*     where x = P*z. If this system does not have full rank, then an */
/*     approximate least squares solution is obtained (see METHOD). */
/*     On output, NF01BQ also provides an upper triangular matrix S */
/*     such that */

/*           P'*(J'*J + D*D)*P = S'*S . */

/*     The system (1) is equivalent to S*z = c , where c contains the */
/*     first n components of the vector obtained by applying to */
/*     [ (Q'*b)'  0 ]' the transformations which triangularized */
/*     [ R'  P'*D*P ]', getting S. */

/*     The matrix R has the following structure */

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
/*             Specifies whether the condition of the matrices S_k should */
/*             be estimated, as follows: */
/*             = 'E' :  use incremental condition estimation and store */
/*                      the numerical rank of S_k in the array entry */
/*                      RANKS(k), for k = 1:l+1; */
/*             = 'N' :  do not use condition estimation, but check the */
/*                      diagonal entries of S_k for zero values; */
/*             = 'U' :  use the ranks already stored in RANKS(1:l+1). */

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
/*             The leading dimension of the array R.  LDR >= MAX(1,N). */

/*     IPVT    (input) INTEGER array, dimension (N) */
/*             This array must define the permutation matrix P such that */
/*             J*P = Q*R. Column j of P is column IPVT(j) of the identity */
/*             matrix. */

/*     DIAG    (input) DOUBLE PRECISION array, dimension (N) */
/*             This array must contain the diagonal elements of the */
/*             matrix D. */

/*     QTB     (input) DOUBLE PRECISION array, dimension (N) */
/*             This array must contain the first n elements of the */
/*             vector Q'*b. */

/*     RANKS   (input or output) INTEGER array, dimension (r), where */
/*             r = BN + 1,  if ST > 0, BSN > 0, and BN > 1; */
/*             r = BN,      if ST = 0 and BSN > 0; */
/*             r = 1,       if ST > 0 and ( BSN = 0 or BN <= 1 ); */
/*             r = 0,       if ST = 0 and BSN = 0. */
/*             On entry, if COND = 'U' and N > 0, this array must contain */
/*             the numerical ranks of the submatrices S_k, k = 1:l(+1). */
/*             On exit, if COND = 'E' or 'N' and N > 0, this array */
/*             contains the numerical ranks of the submatrices S_k, */
/*             k = 1:l(+1), estimated according to the value of COND. */

/*     X       (output) DOUBLE PRECISION array, dimension (N) */
/*             This array contains the least squares solution of the */
/*             system J*x = b, D*x = 0. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If COND = 'E', the tolerance to be used for finding the */
/*             ranks of the submatrices S_k. If the user sets TOL > 0, */
/*             then the given value of TOL is used as a lower bound for */
/*             the reciprocal condition number;  a (sub)matrix whose */
/*             estimated condition number is less than 1/TOL is */
/*             considered to be of full rank.  If the user sets TOL <= 0, */
/*             then an implicitly computed, default tolerance, defined by */
/*             TOLDEF = N*EPS,  is used instead, where EPS is the machine */
/*             precision (see LAPACK Library routine DLAMCH). */
/*             This parameter is not relevant if COND = 'U' or 'N'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, the first N elements of this array contain the */
/*             diagonal elements of the upper triangular matrix S, and */
/*             the next N elements contain the solution z. */
/*             If BN > 1 and BSN > 0, the elements 2*N+1 : 2*N+ST*(N-ST) */
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
/*                                                        COND = 'E'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Standard plane rotations are used to annihilate the elements of */
/*     the diagonal matrix D, updating the upper triangular matrix R */
/*     and the first n elements of the vector Q'*b. A basic least squares */
/*     solution is computed. The computations exploit the special */
/*     structure and storage scheme of the matrix R. If one or more of */
/*     the submatrices S_k, k = 1:l+1, is singular, then the computed */
/*     result is not the basic least squares solution for the whole */
/*     problem, but a concatenation of (least squares) solutions of the */
/*     individual subproblems involving R_k, k = 1:l+1 (with adapted */
/*     right hand sides). */

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
/*     of QRSOLV from the MINPACK package [1], and with optional */
/*     condition estimation. */
/*     The option COND = 'U' is useful when dealing with several */
/*     right-hand side vectors. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     - */

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

#line 282 "NF01BQ.f"
    /* Parameter adjustments */
#line 282 "NF01BQ.f"
    --ipar;
#line 282 "NF01BQ.f"
    r_dim1 = *ldr;
#line 282 "NF01BQ.f"
    r_offset = 1 + r_dim1;
#line 282 "NF01BQ.f"
    r__ -= r_offset;
#line 282 "NF01BQ.f"
    --ipvt;
#line 282 "NF01BQ.f"
    --diag;
#line 282 "NF01BQ.f"
    --qtb;
#line 282 "NF01BQ.f"
    --ranks;
#line 282 "NF01BQ.f"
    --x;
#line 282 "NF01BQ.f"
    --dwork;
#line 282 "NF01BQ.f"

#line 282 "NF01BQ.f"
    /* Function Body */
#line 282 "NF01BQ.f"
    econd = lsame_(cond, "E", (ftnlen)1, (ftnlen)1);
#line 283 "NF01BQ.f"
    *info = 0;
#line 284 "NF01BQ.f"
    if (! (econd || lsame_(cond, "N", (ftnlen)1, (ftnlen)1) || lsame_(cond, 
	    "U", (ftnlen)1, (ftnlen)1))) {
#line 286 "NF01BQ.f"
	*info = -1;
#line 287 "NF01BQ.f"
    } else if (*n < 0) {
#line 288 "NF01BQ.f"
	*info = -2;
#line 289 "NF01BQ.f"
    } else if (*lipar < 4) {
#line 290 "NF01BQ.f"
	*info = -4;
#line 291 "NF01BQ.f"
    } else {
#line 292 "NF01BQ.f"
	st = ipar[1];
#line 293 "NF01BQ.f"
	bn = ipar[2];
#line 294 "NF01BQ.f"
	bsm = ipar[3];
#line 295 "NF01BQ.f"
	bsn = ipar[4];
#line 296 "NF01BQ.f"
	nths = bn * bsn;
/* Computing MIN */
#line 297 "NF01BQ.f"
	i__1 = min(st,bn), i__1 = min(i__1,bsm);
#line 297 "NF01BQ.f"
	if (min(i__1,bsn) < 0) {
#line 298 "NF01BQ.f"
	    *info = -3;
#line 299 "NF01BQ.f"
	} else if (*n != nths + st) {
#line 300 "NF01BQ.f"
	    *info = -2;
#line 301 "NF01BQ.f"
	} else if (*ldr < max(1,*n)) {
#line 302 "NF01BQ.f"
	    *info = -6;
#line 303 "NF01BQ.f"
	} else {
#line 304 "NF01BQ.f"
	    jw = *n << 1;
#line 305 "NF01BQ.f"
	    if (bn <= 1 || bsn == 0) {
#line 306 "NF01BQ.f"
		if (econd) {
#line 306 "NF01BQ.f"
		    jw = *n << 2;
#line 306 "NF01BQ.f"
		}
#line 308 "NF01BQ.f"
	    } else {
#line 309 "NF01BQ.f"
		jw = st * nths + jw;
#line 310 "NF01BQ.f"
		if (econd) {
#line 310 "NF01BQ.f"
		    jw = (max(bsn,st) << 1) + jw;
#line 310 "NF01BQ.f"
		}
#line 312 "NF01BQ.f"
	    }
#line 313 "NF01BQ.f"
	    if (*ldwork < jw) {
#line 313 "NF01BQ.f"
		*info = -14;
#line 313 "NF01BQ.f"
	    }
#line 315 "NF01BQ.f"
	}
#line 316 "NF01BQ.f"
    }

/*     Return if there are illegal arguments. */

#line 320 "NF01BQ.f"
    if (*info != 0) {
#line 321 "NF01BQ.f"
	i__1 = -(*info);
#line 321 "NF01BQ.f"
	xerbla_("NF01BQ", &i__1, (ftnlen)6);
#line 322 "NF01BQ.f"
	return 0;
#line 323 "NF01BQ.f"
    }

/*     Quick return if possible. */

#line 327 "NF01BQ.f"
    if (*n == 0) {
#line 327 "NF01BQ.f"
	return 0;
#line 327 "NF01BQ.f"
    }

#line 330 "NF01BQ.f"
    if (bn <= 1 || bsn == 0) {

/*        Special case: R is an upper triangular matrix. */
/*        Workspace: 4*N, if COND =  'E'; */
/*                   2*N, if COND <> 'E'. */

#line 336 "NF01BQ.f"
	mb02yd_(cond, n, &r__[r_offset], ldr, &ipvt[1], &diag[1], &qtb[1], &
		ranks[1], &x[1], tol, &dwork[1], ldwork, info, (ftnlen)1);
#line 338 "NF01BQ.f"
	return 0;
#line 339 "NF01BQ.f"
    }

/*     General case: BN > 1 and BSN > 0. */
/*     Copy R and Q'*b to preserve input and initialize S. */
/*     In particular, save the diagonal elements of R in X. */

#line 345 "NF01BQ.f"
    ib = *n + 1;
#line 346 "NF01BQ.f"
    is = ib + *n;
#line 347 "NF01BQ.f"
    jw = is + st * nths;
#line 348 "NF01BQ.f"
    i__ = 1;
#line 349 "NF01BQ.f"
    l = is;
#line 350 "NF01BQ.f"
    nc = bsn + st;
#line 351 "NF01BQ.f"
    kf = nc;

#line 353 "NF01BQ.f"
    i__1 = bn;
#line 353 "NF01BQ.f"
    for (k = 1; k <= i__1; ++k) {

#line 355 "NF01BQ.f"
	i__2 = bsn;
#line 355 "NF01BQ.f"
	for (j = 1; j <= i__2; ++j) {
#line 356 "NF01BQ.f"
	    x[i__] = r__[i__ + j * r_dim1];
#line 357 "NF01BQ.f"
	    i__3 = bsn - j + 1;
#line 357 "NF01BQ.f"
	    dcopy_(&i__3, &r__[i__ + j * r_dim1], ldr, &r__[i__ + j * r_dim1],
		     &c__1);
#line 358 "NF01BQ.f"
	    ++i__;
#line 359 "NF01BQ.f"
/* L10: */
#line 359 "NF01BQ.f"
	}

#line 361 "NF01BQ.f"
/* L20: */
#line 361 "NF01BQ.f"
    }

/*     DWORK(IS) contains a copy of [ L_1' ... L_l' ]. */
/*     Workspace:  ST*(N-ST)+2*N; */

#line 366 "NF01BQ.f"
    i__1 = nc;
#line 366 "NF01BQ.f"
    for (j = bsn + 1; j <= i__1; ++j) {
#line 367 "NF01BQ.f"
	dcopy_(&nths, &r__[j * r_dim1 + 1], &c__1, &dwork[l], &st);
#line 368 "NF01BQ.f"
	x[i__] = r__[i__ + j * r_dim1];
#line 369 "NF01BQ.f"
	i__2 = nc - j + 1;
#line 369 "NF01BQ.f"
	dcopy_(&i__2, &r__[i__ + j * r_dim1], ldr, &r__[i__ + j * r_dim1], &
		c__1);
#line 370 "NF01BQ.f"
	++i__;
#line 371 "NF01BQ.f"
	++l;
#line 372 "NF01BQ.f"
/* L30: */
#line 372 "NF01BQ.f"
    }

#line 374 "NF01BQ.f"
    dcopy_(n, &qtb[1], &c__1, &dwork[ib], &c__1);
#line 375 "NF01BQ.f"
    if (st > 0) {
#line 376 "NF01BQ.f"
	itr = nths + 1;
#line 377 "NF01BQ.f"
	itc = bsn + 1;
#line 378 "NF01BQ.f"
    } else {
#line 379 "NF01BQ.f"
	itr = 1;
#line 380 "NF01BQ.f"
	itc = 1;
#line 381 "NF01BQ.f"
    }
#line 382 "NF01BQ.f"
    ibsn = 0;

/*     Eliminate the diagonal matrix D using Givens rotations. */

#line 386 "NF01BQ.f"
    i__1 = *n;
#line 386 "NF01BQ.f"
    for (j = 1; j <= i__1; ++j) {
#line 387 "NF01BQ.f"
	++ibsn;
#line 388 "NF01BQ.f"
	i__ = ibsn;

/*        Prepare the row of D to be eliminated, locating the */
/*        diagonal element using P from the QR factorization. */

#line 393 "NF01BQ.f"
	l = ipvt[j];
#line 394 "NF01BQ.f"
	if (diag[l] != 0.) {
#line 395 "NF01BQ.f"
	    qtbpj = 0.;
#line 396 "NF01BQ.f"
	    dwork[j] = diag[l];

/* Computing MIN */
#line 398 "NF01BQ.f"
	    i__3 = j + kf - 1;
#line 398 "NF01BQ.f"
	    i__2 = min(i__3,*n);
#line 398 "NF01BQ.f"
	    for (k = j + 1; k <= i__2; ++k) {
#line 399 "NF01BQ.f"
		dwork[k] = 0.;
#line 400 "NF01BQ.f"
/* L40: */
#line 400 "NF01BQ.f"
	    }

/*           The transformations to eliminate the row of D modify only */
/*           a single element of Q'*b beyond the first n, which is */
/*           initially zero. */

#line 406 "NF01BQ.f"
	    if (j < nths) {
#line 407 "NF01BQ.f"
		i__2 = bsn - ibsn + 1;
#line 407 "NF01BQ.f"
		mb04ow_(&i__2, &st, &c__1, &r__[j + ibsn * r_dim1], ldr, &r__[
			itr + itc * r_dim1], ldr, &dwork[j], &c__1, &dwork[ib 
			+ j - 1], &bsn, &dwork[ib + nths], &st, &qtbpj, &c__1)
			;
#line 410 "NF01BQ.f"
		if (ibsn == bsn) {
#line 410 "NF01BQ.f"
		    ibsn = 0;
#line 410 "NF01BQ.f"
		}
#line 412 "NF01BQ.f"
	    } else if (j == nths) {
#line 413 "NF01BQ.f"
		mb04ow_(&c__1, &st, &c__1, &r__[j + ibsn * r_dim1], ldr, &r__[
			itr + itc * r_dim1], ldr, &dwork[j], &c__1, &dwork[ib 
			+ j - 1], &bsn, &dwork[ib + nths], &st, &qtbpj, &c__1)
			;
#line 416 "NF01BQ.f"
		kf = st;
#line 417 "NF01BQ.f"
	    } else {
#line 418 "NF01BQ.f"
		i__2 = *n - j + 1;
#line 418 "NF01BQ.f"
		mb04ow_(&c__0, &i__2, &c__1, &r__[j + ibsn * r_dim1], ldr, &
			r__[j + ibsn * r_dim1], ldr, &dwork[j], &c__1, &dwork[
			ib + j - 1], &c__1, &dwork[ib + j - 1], &st, &qtbpj, &
			c__1);
#line 421 "NF01BQ.f"
	    }
#line 422 "NF01BQ.f"
	} else {
#line 423 "NF01BQ.f"
	    if (j < nths) {
#line 424 "NF01BQ.f"
		if (ibsn == bsn) {
#line 424 "NF01BQ.f"
		    ibsn = 0;
#line 424 "NF01BQ.f"
		}
#line 426 "NF01BQ.f"
	    } else if (j == nths) {
#line 427 "NF01BQ.f"
		kf = st;
#line 428 "NF01BQ.f"
	    }
#line 429 "NF01BQ.f"
	}

/*        Store the diagonal element of S. */

#line 433 "NF01BQ.f"
	dwork[j] = r__[j + i__ * r_dim1];
#line 434 "NF01BQ.f"
/* L50: */
#line 434 "NF01BQ.f"
    }

/*     Solve the triangular system for z. If the system is singular, */
/*     then obtain an approximate least squares solution. */
/*     Additional workspace:   2*MAX(BSN,ST), if COND =  'E'; */
/*                             0,             if COND <> 'E'. */

#line 441 "NF01BQ.f"
    i__1 = *ldwork - jw + 1;
#line 441 "NF01BQ.f"
    nf01br_(cond, "Upper", "NoTranspose", n, &ipar[1], lipar, &r__[r_offset], 
	    ldr, &dwork[1], &dwork[is], &c__1, &dwork[ib], &ranks[1], tol, &
	    dwork[jw], &i__1, info, (ftnlen)1, (ftnlen)5, (ftnlen)11);
#line 444 "NF01BQ.f"
    i__ = 1;

/*     Restore the diagonal elements of R from X and interchange */
/*     the upper and lower triangular parts of R. */

#line 449 "NF01BQ.f"
    i__1 = bn;
#line 449 "NF01BQ.f"
    for (k = 1; k <= i__1; ++k) {

#line 451 "NF01BQ.f"
	i__2 = bsn;
#line 451 "NF01BQ.f"
	for (j = 1; j <= i__2; ++j) {
#line 452 "NF01BQ.f"
	    r__[i__ + j * r_dim1] = x[i__];
#line 453 "NF01BQ.f"
	    i__3 = bsn - j + 1;
#line 453 "NF01BQ.f"
	    dswap_(&i__3, &r__[i__ + j * r_dim1], ldr, &r__[i__ + j * r_dim1],
		     &c__1);
#line 454 "NF01BQ.f"
	    ++i__;
#line 455 "NF01BQ.f"
/* L60: */
#line 455 "NF01BQ.f"
	}

#line 457 "NF01BQ.f"
/* L70: */
#line 457 "NF01BQ.f"
    }

#line 459 "NF01BQ.f"
    i__1 = nc;
#line 459 "NF01BQ.f"
    for (j = bsn + 1; j <= i__1; ++j) {
#line 460 "NF01BQ.f"
	dswap_(&nths, &r__[j * r_dim1 + 1], &c__1, &dwork[is], &st);
#line 461 "NF01BQ.f"
	r__[i__ + j * r_dim1] = x[i__];
#line 462 "NF01BQ.f"
	i__2 = nc - j + 1;
#line 462 "NF01BQ.f"
	dswap_(&i__2, &r__[i__ + j * r_dim1], ldr, &r__[i__ + j * r_dim1], &
		c__1);
#line 463 "NF01BQ.f"
	++i__;
#line 464 "NF01BQ.f"
	++is;
#line 465 "NF01BQ.f"
/* L80: */
#line 465 "NF01BQ.f"
    }

/*     Permute the components of z back to components of x. */

#line 469 "NF01BQ.f"
    i__1 = *n;
#line 469 "NF01BQ.f"
    for (j = 1; j <= i__1; ++j) {
#line 470 "NF01BQ.f"
	l = ipvt[j];
#line 471 "NF01BQ.f"
	x[l] = dwork[*n + j];
#line 472 "NF01BQ.f"
/* L90: */
#line 472 "NF01BQ.f"
    }

#line 474 "NF01BQ.f"
    return 0;

/* *** Last line of NF01BQ *** */
} /* nf01bq_ */

