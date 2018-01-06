#line 1 "NF01BR.f"
/* NF01BR.f -- translated by f2c (version 20100827).
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

#line 1 "NF01BR.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b19 = 0.;
static integer c__0 = 0;
static doublereal c_b56 = -1.;
static doublereal c_b58 = 1.;

/* Subroutine */ int nf01br_(char *cond, char *uplo, char *trans, integer *n, 
	integer *ipar, integer *lipar, doublereal *r__, integer *ldr, 
	doublereal *sdiag, doublereal *s, integer *lds, doublereal *b, 
	integer *ranks, doublereal *tol, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen cond_len, ftnlen uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    integer r_dim1, r_offset, s_dim1, s_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, l, i1, bn, nc, st, bsm, bsn;
    static doublereal dum[3];
    static integer rank;
    static logical full;
    static integer nths;
    extern /* Subroutine */ int mb03od_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen);
    static logical econd, ncond;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), dswap_(integer 
	    *, doublereal *, integer *, doublereal *, integer *);
    static logical lower, tranr;
    static char uplol[1];
    extern /* Subroutine */ int dtrsv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static char transl[1];


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

/*     To solve one of the systems of linear equations */

/*           R*x = b ,  or  R'*x = b , */

/*     in the least squares sense, where R is an n-by-n block upper */
/*     triangular matrix, with the structure */

/*         /   R_1    0    ..   0   |   L_1   \ */
/*         |    0    R_2   ..   0   |   L_2   | */
/*         |    :     :    ..   :   |    :    | , */
/*         |    0     0    ..  R_l  |   L_l   | */
/*         \    0     0    ..   0   |  R_l+1  / */

/*     with the upper triangular submatrices R_k, k = 1:l+1, square, and */
/*     the first l of the same order, BSN. The diagonal elements of each */
/*     block R_k have nonincreasing magnitude. The matrix R is stored in */
/*     the compressed form, as returned by SLICOT Library routine NF01BS, */

/*              /   R_1  |   L_1   \ */
/*              |   R_2  |   L_2   | */
/*       Rc =   |    :   |    :    | , */
/*              |   R_l  |   L_l   | */
/*              \    X   |  R_l+1  / */

/*     where the submatrix X is irrelevant. If the matrix R does not have */
/*     full rank, then a least squares solution is obtained. If l <= 1, */
/*     then R is an upper triangular matrix and its full upper triangle */
/*     is stored. */

/*     Optionally, the transpose of the matrix R can be stored in the */
/*     strict lower triangles of the submatrices R_k, k = 1:l+1, and in */
/*     the arrays SDIAG and S, as described at the parameter UPLO below. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     COND    CHARACTER*1 */
/*             Specifies whether the condition of submatrices R_k should */
/*             be estimated, as follows: */
/*             = 'E' :  use incremental condition estimation and store */
/*                      the numerical rank of R_k in the array entry */
/*                      RANKS(k), for k = 1:l+1; */
/*             = 'N' :  do not use condition estimation, but check the */
/*                      diagonal entries of R_k for zero values; */
/*             = 'U' :  use the ranks already stored in RANKS(1:l+1). */

/*     UPLO    CHARACTER*1 */
/*             Specifies the storage scheme for the matrix R, as follows: */
/*             = 'U' :  the upper triangular part is stored as in Rc; */
/*             = 'L' :  the lower triangular part is stored, namely, */
/*                      - the transpose of the strict upper triangle of */
/*                        R_k is stored in the strict lower triangle of */
/*                        R_k, for k = 1:l+1; */
/*                      - the diagonal elements of R_k, k = 1:l+1, are */
/*                        stored in the array SDIAG; */
/*                      - the transpose of the last block column in R */
/*                        (without R_l+1) is stored in the array S. */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of the system of equations, as follows: */
/*             = 'N':  R*x  = b  (No transpose); */
/*             = 'T':  R'*x = b  (Transpose); */
/*             = 'C':  R'*x = b  (Transpose). */

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

/*     R       (input) DOUBLE PRECISION array, dimension (LDR, NC) */
/*             where NC = N if BN <= 1, and NC = BSN+ST, if BN > 1. */
/*             If UPLO = 'U', the leading N-by-NC part of this array must */
/*             contain the (compressed) representation (Rc) of the upper */
/*             triangular matrix R. The submatrix X in Rc and the strict */
/*             lower triangular parts of the diagonal blocks R_k, */
/*             k = 1:l+1, are not referenced. If BN <= 1 or BSN = 0, then */
/*             the full upper triangle of R must be stored. */
/*             If UPLO = 'L', BN > 1 and BSN > 0, the leading */
/*             (N-ST)-by-BSN part of this array must contain the */
/*             transposes of the strict upper triangles of R_k, k = 1:l, */
/*             stored in the strict lower triangles of R_k, and the */
/*             strict lower triangle of R_l+1 must contain the transpose */
/*             of the strict upper triangle of R_l+1. The submatrix X */
/*             in Rc is not referenced. The diagonal elements of R_k, */
/*             and, if COND = 'E', the upper triangular parts of R_k, */
/*             k = 1:l+1, are modified internally, but are restored */
/*             on exit. */
/*             If UPLO = 'L' and BN <= 1 or BSN = 0, the leading N-by-N */
/*             strict lower triangular part of this array must contain */
/*             the transpose of the strict upper triangular part of R. */
/*             The diagonal elements and, if COND = 'E', the upper */
/*             triangular elements are modified internally, but are */
/*             restored on exit. */

/*     LDR     INTEGER */
/*             The leading dimension of the array R.  LDR >= MAX(1,N). */

/*     SDIAG   (input) DOUBLE PRECISION array, dimension (N) */
/*             If UPLO = 'L', this array must contain the diagonal */
/*             entries of R_k, k = 1:l+1. This array is modified */
/*             internally, but is restored on exit. */
/*             This parameter is not referenced if UPLO = 'U'. */

/*     S       (input) DOUBLE PRECISION array, dimension (LDS,N-ST) */
/*             If UPLO = 'L', BN > 1, and BSN > 0, the leading */
/*             ST-by-(N-ST) part of this array must contain the transpose */
/*             of the rectangular part of the last block column in R, */
/*             that is [ L_1' L_2' ... L_l' ] . If COND = 'E', S is */
/*             modified internally, but is restored on exit. */
/*             This parameter is not referenced if UPLO = 'U', or */
/*             BN <= 1, or BSN = 0. */

/*     LDS     INTEGER */
/*             The leading dimension of the array S. */
/*             LDS >= 1,         if UPLO = 'U', or BN <= 1, or BSN = 0; */
/*             LDS >= MAX(1,ST), if UPLO = 'L', BN > 1, and BSN > 0. */

/*     B       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the right hand side */
/*             vector b. */
/*             On exit, this array contains the (least squares) solution */
/*             of the system R*x = b or R'*x = b. */

/*     RANKS   (input or output) INTEGER array, dimension (r), where */
/*             r = BN + 1,  if ST > 0, BSN > 0, and BN > 1; */
/*             r = BN,      if ST = 0 and BSN > 0; */
/*             r = 1,       if ST > 0 and ( BSN = 0 or BN <= 1 ); */
/*             r = 0,       if ST = 0 and BSN = 0. */
/*             On entry, if COND = 'U' and N > 0, this array must contain */
/*             the numerical ranks of the submatrices R_k, k = 1:l(+1). */
/*             On exit, if COND = 'E' or 'N' and N > 0, this array */
/*             contains the numerical ranks of the submatrices R_k, */
/*             k = 1:l(+1), estimated according to the value of COND. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If COND = 'E', the tolerance to be used for finding the */
/*             ranks of the submatrices R_k. If the user sets TOL > 0, */
/*             then the given value of TOL is used as a lower bound for */
/*             the reciprocal condition number;  a (sub)matrix whose */
/*             estimated condition number is less than 1/TOL is */
/*             considered to be of full rank. If the user sets TOL <= 0, */
/*             then an implicitly computed, default tolerance, defined by */
/*             TOLDEF = N*EPS,  is used instead, where EPS is the machine */
/*             precision (see LAPACK Library routine DLAMCH). */
/*             This parameter is not relevant if COND = 'U' or 'N'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             Denote Full = ( BN <= 1 or  BSN = 0 ); */
/*                    Comp = ( BN >  1 and BSN > 0 ). */
/*             LDWORK >= 2*N,           if Full and COND = 'E'; */
/*             LDWORK >= 2*MAX(BSN,ST), if Comp and COND = 'E'; */
/*             LDWORK >= 0,   in the remaining cases. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Block back or forward substitution is used (depending on TRANS */
/*     and UPLO), exploiting the special structure and storage scheme of */
/*     the matrix R. If a submatrix R_k, k = 1:l+1, is singular, a local */
/*     basic least squares solution is computed. Therefore, the returned */
/*     result is not the basic least squares solution for the whole */
/*     problem, but a concatenation of (least squares) solutions of the */
/*     individual subproblems involving R_k, k = 1:l+1 (with adapted */
/*     right hand sides). */

/*     NUMERICAL ASPECTS */
/*                                    2    2 */
/*     The algorithm requires 0(BN*BSN + ST + N*ST) operations and is */
/*     backward stable, if R is nonsingular. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2005. */

/*     KEYWORDS */

/*     Linear system of equations, matrix operations, plane rotations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

#line 270 "NF01BR.f"
    /* Parameter adjustments */
#line 270 "NF01BR.f"
    --ipar;
#line 270 "NF01BR.f"
    r_dim1 = *ldr;
#line 270 "NF01BR.f"
    r_offset = 1 + r_dim1;
#line 270 "NF01BR.f"
    r__ -= r_offset;
#line 270 "NF01BR.f"
    --sdiag;
#line 270 "NF01BR.f"
    s_dim1 = *lds;
#line 270 "NF01BR.f"
    s_offset = 1 + s_dim1;
#line 270 "NF01BR.f"
    s -= s_offset;
#line 270 "NF01BR.f"
    --b;
#line 270 "NF01BR.f"
    --ranks;
#line 270 "NF01BR.f"
    --dwork;
#line 270 "NF01BR.f"

#line 270 "NF01BR.f"
    /* Function Body */
#line 270 "NF01BR.f"
    econd = lsame_(cond, "E", (ftnlen)1, (ftnlen)1);
#line 271 "NF01BR.f"
    ncond = lsame_(cond, "N", (ftnlen)1, (ftnlen)1);
#line 272 "NF01BR.f"
    lower = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 273 "NF01BR.f"
    tranr = lsame_(trans, "T", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1);

#line 275 "NF01BR.f"
    *info = 0;
#line 276 "NF01BR.f"
    if (! (econd || ncond || lsame_(cond, "U", (ftnlen)1, (ftnlen)1))) {
#line 277 "NF01BR.f"
	*info = -1;
#line 278 "NF01BR.f"
    } else if (! (lower || lsame_(uplo, "U", (ftnlen)1, (ftnlen)1))) {
#line 279 "NF01BR.f"
	*info = -2;
#line 280 "NF01BR.f"
    } else if (! (tranr || lsame_(trans, "N", (ftnlen)1, (ftnlen)1))) {
#line 281 "NF01BR.f"
	*info = -3;
#line 282 "NF01BR.f"
    } else if (*n < 0) {
#line 283 "NF01BR.f"
	*info = -4;
#line 284 "NF01BR.f"
    } else if (*lipar < 4) {
#line 285 "NF01BR.f"
	*info = -6;
#line 286 "NF01BR.f"
    } else {
#line 287 "NF01BR.f"
	st = ipar[1];
#line 288 "NF01BR.f"
	bn = ipar[2];
#line 289 "NF01BR.f"
	bsm = ipar[3];
#line 290 "NF01BR.f"
	bsn = ipar[4];
#line 291 "NF01BR.f"
	nths = bn * bsn;
#line 292 "NF01BR.f"
	full = bn <= 1 || bsn == 0;
/* Computing MIN */
#line 293 "NF01BR.f"
	i__1 = min(st,bn), i__1 = min(i__1,bsm);
#line 293 "NF01BR.f"
	if (min(i__1,bsn) < 0) {
#line 294 "NF01BR.f"
	    *info = -5;
#line 295 "NF01BR.f"
	} else if (*n != nths + st) {
#line 296 "NF01BR.f"
	    *info = -4;
#line 297 "NF01BR.f"
	} else if (*ldr < max(1,*n)) {
#line 298 "NF01BR.f"
	    *info = -8;
#line 299 "NF01BR.f"
	} else if (*lds < 1 || lower && ! full && *lds < st) {
#line 301 "NF01BR.f"
	    *info = -11;
#line 302 "NF01BR.f"
	} else {
#line 303 "NF01BR.f"
	    if (econd) {
#line 304 "NF01BR.f"
		if (full) {
#line 305 "NF01BR.f"
		    l = *n << 1;
#line 306 "NF01BR.f"
		} else {
#line 307 "NF01BR.f"
		    l = max(bsn,st) << 1;
#line 308 "NF01BR.f"
		}
#line 309 "NF01BR.f"
	    } else {
#line 310 "NF01BR.f"
		l = 0;
#line 311 "NF01BR.f"
	    }
#line 312 "NF01BR.f"
	    if (*ldwork < l) {
#line 312 "NF01BR.f"
		*info = -16;
#line 312 "NF01BR.f"
	    }
#line 314 "NF01BR.f"
	}
#line 315 "NF01BR.f"
    }

/*     Return if there are illegal arguments. */

#line 319 "NF01BR.f"
    if (*info != 0) {
#line 320 "NF01BR.f"
	i__1 = -(*info);
#line 320 "NF01BR.f"
	xerbla_("NF01BR", &i__1, (ftnlen)6);
#line 321 "NF01BR.f"
	return 0;
#line 322 "NF01BR.f"
    }

/*     Quick return if possible. */

#line 326 "NF01BR.f"
    if (*n == 0) {
#line 326 "NF01BR.f"
	return 0;
#line 326 "NF01BR.f"
    }

#line 329 "NF01BR.f"
    if (econd) {
#line 330 "NF01BR.f"
	toldef = *tol;
#line 331 "NF01BR.f"
	if (toldef <= 0.) {

/*           Use the default tolerance in rank determination. */

#line 335 "NF01BR.f"
	    toldef = (doublereal) (*n) * dlamch_("Epsilon", (ftnlen)7);
#line 336 "NF01BR.f"
	}
#line 337 "NF01BR.f"
    }

#line 339 "NF01BR.f"
    nc = bsn + st;
#line 340 "NF01BR.f"
    if (full) {

/*        Special case: l <= 1 or BSN = 0; R is just an upper triangular */
/*        matrix. */

#line 345 "NF01BR.f"
	if (lower) {

/*           Swap the diagonal elements of R and the elements of SDIAG */
/*           and, if COND = 'E', swap the upper and lower triangular */
/*           parts of R, in order to find the numerical rank. */

#line 351 "NF01BR.f"
	    i__1 = *ldr + 1;
#line 351 "NF01BR.f"
	    dswap_(n, &r__[r_offset], &i__1, &sdiag[1], &c__1);
#line 352 "NF01BR.f"
	    if (econd) {
#line 353 "NF01BR.f"
		*(unsigned char *)uplol = 'U';
#line 354 "NF01BR.f"
		*(unsigned char *)transl = *(unsigned char *)trans;

#line 356 "NF01BR.f"
		i__1 = *n;
#line 356 "NF01BR.f"
		for (j = 1; j <= i__1; ++j) {
#line 357 "NF01BR.f"
		    i__2 = *n - j + 1;
#line 357 "NF01BR.f"
		    dswap_(&i__2, &r__[j + j * r_dim1], ldr, &r__[j + j * 
			    r_dim1], &c__1);
#line 358 "NF01BR.f"
/* L10: */
#line 358 "NF01BR.f"
		}

#line 360 "NF01BR.f"
	    } else {
#line 361 "NF01BR.f"
		*(unsigned char *)uplol = *(unsigned char *)uplo;
#line 362 "NF01BR.f"
		if (tranr) {
#line 363 "NF01BR.f"
		    *(unsigned char *)transl = 'N';
#line 364 "NF01BR.f"
		} else {
#line 365 "NF01BR.f"
		    *(unsigned char *)transl = 'T';
#line 366 "NF01BR.f"
		}
#line 367 "NF01BR.f"
	    }
#line 368 "NF01BR.f"
	} else {
#line 369 "NF01BR.f"
	    *(unsigned char *)uplol = *(unsigned char *)uplo;
#line 370 "NF01BR.f"
	    *(unsigned char *)transl = *(unsigned char *)trans;
#line 371 "NF01BR.f"
	}

#line 373 "NF01BR.f"
	if (econd) {

/*           Estimate the reciprocal condition number and set the rank. */
/*           Workspace: 2*N. */

#line 378 "NF01BR.f"
	    mb03od_("No QR", n, n, &r__[r_offset], ldr, &ipar[1], &toldef, &
		    c_b19, &dwork[1], &rank, dum, &dwork[1], ldwork, info, (
		    ftnlen)5);
#line 380 "NF01BR.f"
	    ranks[1] = rank;

#line 382 "NF01BR.f"
	} else if (ncond) {

/*           Determine rank(R) by checking zero diagonal entries. */

#line 386 "NF01BR.f"
	    rank = *n;

#line 388 "NF01BR.f"
	    i__1 = *n;
#line 388 "NF01BR.f"
	    for (j = 1; j <= i__1; ++j) {
#line 389 "NF01BR.f"
		if (r__[j + j * r_dim1] == 0. && rank == *n) {
#line 389 "NF01BR.f"
		    rank = j - 1;
#line 389 "NF01BR.f"
		}
#line 391 "NF01BR.f"
/* L20: */
#line 391 "NF01BR.f"
	    }

#line 393 "NF01BR.f"
	    ranks[1] = rank;

#line 395 "NF01BR.f"
	} else {

/*           Use the stored rank. */

#line 399 "NF01BR.f"
	    rank = ranks[1];
#line 400 "NF01BR.f"
	}

/*        Solve R*x = b, or R'*x = b using back or forward substitution. */

#line 404 "NF01BR.f"
	dum[0] = 0.;
#line 405 "NF01BR.f"
	if (rank < *n) {
#line 405 "NF01BR.f"
	    i__1 = *n - rank;
#line 405 "NF01BR.f"
	    dcopy_(&i__1, dum, &c__0, &b[rank + 1], &c__1);
#line 405 "NF01BR.f"
	}
#line 407 "NF01BR.f"
	dtrsv_(uplol, transl, "NonUnit", &rank, &r__[r_offset], ldr, &b[1], &
		c__1, (ftnlen)1, (ftnlen)1, (ftnlen)7);

#line 409 "NF01BR.f"
	if (lower) {

/*           Swap the diagonal elements of R and the elements of SDIAG */
/*           and, if COND = 'E', swap back the upper and lower triangular */
/*           parts of R. */

#line 415 "NF01BR.f"
	    i__1 = *ldr + 1;
#line 415 "NF01BR.f"
	    dswap_(n, &r__[r_offset], &i__1, &sdiag[1], &c__1);
#line 416 "NF01BR.f"
	    if (econd) {

#line 418 "NF01BR.f"
		i__1 = *n;
#line 418 "NF01BR.f"
		for (j = 1; j <= i__1; ++j) {
#line 419 "NF01BR.f"
		    i__2 = *n - j + 1;
#line 419 "NF01BR.f"
		    dswap_(&i__2, &r__[j + j * r_dim1], ldr, &r__[j + j * 
			    r_dim1], &c__1);
#line 420 "NF01BR.f"
/* L30: */
#line 420 "NF01BR.f"
		}

#line 422 "NF01BR.f"
	    }

#line 424 "NF01BR.f"
	}
#line 425 "NF01BR.f"
	return 0;
#line 426 "NF01BR.f"
    }

/*     General case: l > 1 and BSN > 0. */

#line 430 "NF01BR.f"
    i__ = 1;
#line 431 "NF01BR.f"
    l = bn;
#line 432 "NF01BR.f"
    if (econd) {

/*        Estimate the reciprocal condition numbers and set the ranks. */

#line 436 "NF01BR.f"
	if (lower) {

/*           Swap the diagonal elements of R and the elements of SDIAG */
/*           and swap the upper and lower triangular parts of R, in order */
/*           to find the numerical rank. Swap S and the transpose of the */
/*           rectangular part of the last block column of R. */

#line 443 "NF01BR.f"
	    i__1 = bn;
#line 443 "NF01BR.f"
	    for (k = 1; k <= i__1; ++k) {
#line 444 "NF01BR.f"
		i__2 = *ldr + 1;
#line 444 "NF01BR.f"
		dswap_(&bsn, &r__[i__ + r_dim1], &i__2, &sdiag[i__], &c__1);

#line 446 "NF01BR.f"
		i__2 = bsn;
#line 446 "NF01BR.f"
		for (j = 1; j <= i__2; ++j) {
#line 447 "NF01BR.f"
		    i__3 = bsn - j + 1;
#line 447 "NF01BR.f"
		    dswap_(&i__3, &r__[i__ + j * r_dim1], ldr, &r__[i__ + j * 
			    r_dim1], &c__1);
#line 448 "NF01BR.f"
		    ++i__;
#line 449 "NF01BR.f"
/* L40: */
#line 449 "NF01BR.f"
		}

#line 451 "NF01BR.f"
/* L50: */
#line 451 "NF01BR.f"
	    }

#line 453 "NF01BR.f"
	    if (st > 0) {
#line 454 "NF01BR.f"
		i__1 = *ldr + 1;
#line 454 "NF01BR.f"
		dswap_(&st, &r__[i__ + (bsn + 1) * r_dim1], &i__1, &sdiag[i__]
			, &c__1);

#line 456 "NF01BR.f"
		i__1 = nc;
#line 456 "NF01BR.f"
		for (j = bsn + 1; j <= i__1; ++j) {
#line 457 "NF01BR.f"
		    dswap_(&nths, &r__[j * r_dim1 + 1], &c__1, &s[j - bsn + 
			    s_dim1], lds);
#line 458 "NF01BR.f"
		    i__2 = nc - j + 1;
#line 458 "NF01BR.f"
		    dswap_(&i__2, &r__[i__ + j * r_dim1], ldr, &r__[i__ + j * 
			    r_dim1], &c__1);
#line 459 "NF01BR.f"
		    ++i__;
#line 460 "NF01BR.f"
/* L60: */
#line 460 "NF01BR.f"
		}

#line 462 "NF01BR.f"
	    }

#line 464 "NF01BR.f"
	}

#line 466 "NF01BR.f"
	i1 = 1;

/*        Determine rank(R_k) using incremental condition estimation. */
/*        Workspace 2*MAX(BSN,ST). */

#line 471 "NF01BR.f"
	i__1 = bn;
#line 471 "NF01BR.f"
	for (k = 1; k <= i__1; ++k) {
#line 472 "NF01BR.f"
	    mb03od_("No QR", &bsn, &bsn, &r__[i1 + r_dim1], ldr, &ipar[1], &
		    toldef, &c_b19, &dwork[1], &ranks[k], dum, &dwork[1], 
		    ldwork, info, (ftnlen)5);
#line 475 "NF01BR.f"
	    i1 += bsn;
#line 476 "NF01BR.f"
/* L70: */
#line 476 "NF01BR.f"
	}

#line 478 "NF01BR.f"
	if (st > 0) {
#line 479 "NF01BR.f"
	    ++l;
#line 480 "NF01BR.f"
	    mb03od_("No QR", &st, &st, &r__[i1 + (bsn + 1) * r_dim1], ldr, &
		    ipar[1], &toldef, &c_b19, &dwork[1], &ranks[l], dum, &
		    dwork[1], ldwork, info, (ftnlen)5);
#line 483 "NF01BR.f"
	}

#line 485 "NF01BR.f"
    } else if (ncond) {

/*        Determine rank(R_k) by checking zero diagonal entries. */

#line 489 "NF01BR.f"
	if (lower) {

#line 491 "NF01BR.f"
	    i__1 = bn;
#line 491 "NF01BR.f"
	    for (k = 1; k <= i__1; ++k) {
#line 492 "NF01BR.f"
		rank = bsn;

#line 494 "NF01BR.f"
		i__2 = bsn;
#line 494 "NF01BR.f"
		for (j = 1; j <= i__2; ++j) {
#line 495 "NF01BR.f"
		    if (sdiag[i__] == 0. && rank == bsn) {
#line 495 "NF01BR.f"
			rank = j - 1;
#line 495 "NF01BR.f"
		    }
#line 497 "NF01BR.f"
		    ++i__;
#line 498 "NF01BR.f"
/* L80: */
#line 498 "NF01BR.f"
		}

#line 500 "NF01BR.f"
		ranks[k] = rank;
#line 501 "NF01BR.f"
/* L90: */
#line 501 "NF01BR.f"
	    }

#line 503 "NF01BR.f"
	    if (st > 0) {
#line 504 "NF01BR.f"
		++l;
#line 505 "NF01BR.f"
		rank = st;

#line 507 "NF01BR.f"
		i__1 = st;
#line 507 "NF01BR.f"
		for (j = 1; j <= i__1; ++j) {
#line 508 "NF01BR.f"
		    if (sdiag[i__] == 0. && rank == st) {
#line 508 "NF01BR.f"
			rank = j - 1;
#line 508 "NF01BR.f"
		    }
#line 510 "NF01BR.f"
		    ++i__;
#line 511 "NF01BR.f"
/* L100: */
#line 511 "NF01BR.f"
		}

#line 513 "NF01BR.f"
		ranks[l] = rank;
#line 514 "NF01BR.f"
	    }

#line 516 "NF01BR.f"
	} else {

#line 518 "NF01BR.f"
	    i__1 = bn;
#line 518 "NF01BR.f"
	    for (k = 1; k <= i__1; ++k) {
#line 519 "NF01BR.f"
		rank = bsn;

#line 521 "NF01BR.f"
		i__2 = bsn;
#line 521 "NF01BR.f"
		for (j = 1; j <= i__2; ++j) {
#line 522 "NF01BR.f"
		    if (r__[i__ + j * r_dim1] == 0. && rank == bsn) {
#line 522 "NF01BR.f"
			rank = j - 1;
#line 522 "NF01BR.f"
		    }
#line 524 "NF01BR.f"
		    ++i__;
#line 525 "NF01BR.f"
/* L110: */
#line 525 "NF01BR.f"
		}

#line 527 "NF01BR.f"
		ranks[k] = rank;
#line 528 "NF01BR.f"
/* L120: */
#line 528 "NF01BR.f"
	    }

#line 530 "NF01BR.f"
	    if (st > 0) {
#line 531 "NF01BR.f"
		++l;
#line 532 "NF01BR.f"
		rank = st;

#line 534 "NF01BR.f"
		i__1 = nc;
#line 534 "NF01BR.f"
		for (j = bsn + 1; j <= i__1; ++j) {
#line 535 "NF01BR.f"
		    if (r__[i__ + j * r_dim1] == 0. && rank == st) {
#line 535 "NF01BR.f"
			rank = j - bsn - 1;
#line 535 "NF01BR.f"
		    }
#line 537 "NF01BR.f"
		    ++i__;
#line 538 "NF01BR.f"
/* L130: */
#line 538 "NF01BR.f"
		}

#line 540 "NF01BR.f"
		ranks[l] = rank;
#line 541 "NF01BR.f"
	    }
#line 542 "NF01BR.f"
	}

#line 544 "NF01BR.f"
    } else {

/*        Set the number of elements of RANKS. Then use the stored ranks. */

#line 548 "NF01BR.f"
	if (st > 0) {
#line 548 "NF01BR.f"
	    ++l;
#line 548 "NF01BR.f"
	}
#line 550 "NF01BR.f"
    }

/*     Solve the triangular system for x. If the system is singular, */
/*     then obtain a basic least squares solution. */

#line 555 "NF01BR.f"
    dum[0] = 0.;
#line 556 "NF01BR.f"
    if (lower && ! econd) {

#line 558 "NF01BR.f"
	if (! tranr) {

/*           Solve R*x = b using back substitution, with R' stored in */
/*           the arrays R, SDIAG and S. Swap diag(R) and SDIAG. */

#line 563 "NF01BR.f"
	    i1 = nths + 1;
#line 564 "NF01BR.f"
	    if (st > 0) {
#line 565 "NF01BR.f"
		rank = ranks[l];
#line 566 "NF01BR.f"
		if (rank < st) {
#line 566 "NF01BR.f"
		    i__1 = st - rank;
#line 566 "NF01BR.f"
		    dcopy_(&i__1, dum, &c__0, &b[i1 + rank], &c__1);
#line 566 "NF01BR.f"
		}
#line 568 "NF01BR.f"
		i__1 = *ldr + 1;
#line 568 "NF01BR.f"
		dswap_(&st, &r__[i1 + (bsn + 1) * r_dim1], &i__1, &sdiag[i1], 
			&c__1);
#line 569 "NF01BR.f"
		dtrsv_("Lower", "Transpose", "NonUnit", &rank, &r__[i1 + (bsn 
			+ 1) * r_dim1], ldr, &b[i1], &c__1, (ftnlen)5, (
			ftnlen)9, (ftnlen)7);
#line 571 "NF01BR.f"
		i__1 = *ldr + 1;
#line 571 "NF01BR.f"
		dswap_(&st, &r__[i1 + (bsn + 1) * r_dim1], &i__1, &sdiag[i1], 
			&c__1);
#line 572 "NF01BR.f"
		dgemv_("Transpose", &st, &nths, &c_b56, &s[s_offset], lds, &b[
			nths + 1], &c__1, &c_b58, &b[1], &c__1, (ftnlen)9);
#line 574 "NF01BR.f"
	    }

#line 576 "NF01BR.f"
	    for (k = bn; k >= 1; --k) {
#line 577 "NF01BR.f"
		i1 -= bsn;
#line 578 "NF01BR.f"
		rank = ranks[k];
#line 579 "NF01BR.f"
		if (rank < bsn) {
#line 579 "NF01BR.f"
		    i__1 = bsn - rank;
#line 579 "NF01BR.f"
		    dcopy_(&i__1, dum, &c__0, &b[i1 + rank], &c__1);
#line 579 "NF01BR.f"
		}
#line 581 "NF01BR.f"
		i__1 = *ldr + 1;
#line 581 "NF01BR.f"
		dswap_(&bsn, &r__[i1 + r_dim1], &i__1, &sdiag[i1], &c__1);
#line 582 "NF01BR.f"
		dtrsv_("Lower", "Transpose", "NonUnit", &rank, &r__[i1 + 
			r_dim1], ldr, &b[i1], &c__1, (ftnlen)5, (ftnlen)9, (
			ftnlen)7);
#line 584 "NF01BR.f"
		i__1 = *ldr + 1;
#line 584 "NF01BR.f"
		dswap_(&bsn, &r__[i1 + r_dim1], &i__1, &sdiag[i1], &c__1);
#line 585 "NF01BR.f"
/* L140: */
#line 585 "NF01BR.f"
	    }

#line 587 "NF01BR.f"
	} else {

/*           Solve R'*x = b using forward substitution, with R' stored in */
/*           the arrays R, SDIAG and S. Swap diag(R) and SDIAG. */

#line 592 "NF01BR.f"
	    i1 = 1;
#line 593 "NF01BR.f"
	    if (tranr) {
#line 594 "NF01BR.f"
		*(unsigned char *)transl = 'N';
#line 595 "NF01BR.f"
	    } else {
#line 596 "NF01BR.f"
		*(unsigned char *)transl = 'T';
#line 597 "NF01BR.f"
	    }

#line 599 "NF01BR.f"
	    i__1 = bn;
#line 599 "NF01BR.f"
	    for (k = 1; k <= i__1; ++k) {
#line 600 "NF01BR.f"
		rank = ranks[k];
#line 601 "NF01BR.f"
		if (rank < bsn) {
#line 601 "NF01BR.f"
		    i__2 = bsn - rank;
#line 601 "NF01BR.f"
		    dcopy_(&i__2, dum, &c__0, &b[i1 + rank], &c__1);
#line 601 "NF01BR.f"
		}
#line 603 "NF01BR.f"
		i__2 = *ldr + 1;
#line 603 "NF01BR.f"
		dswap_(&bsn, &r__[i1 + r_dim1], &i__2, &sdiag[i1], &c__1);
#line 604 "NF01BR.f"
		dtrsv_("Lower", transl, "NonUnit", &rank, &r__[i1 + r_dim1], 
			ldr, &b[i1], &c__1, (ftnlen)5, (ftnlen)1, (ftnlen)7);
#line 606 "NF01BR.f"
		i__2 = *ldr + 1;
#line 606 "NF01BR.f"
		dswap_(&bsn, &r__[i1 + r_dim1], &i__2, &sdiag[i1], &c__1);
#line 607 "NF01BR.f"
		i1 += bsn;
#line 608 "NF01BR.f"
/* L150: */
#line 608 "NF01BR.f"
	    }

#line 610 "NF01BR.f"
	    if (st > 0) {
#line 611 "NF01BR.f"
		rank = ranks[l];
#line 612 "NF01BR.f"
		if (rank < st) {
#line 612 "NF01BR.f"
		    i__1 = st - rank;
#line 612 "NF01BR.f"
		    dcopy_(&i__1, dum, &c__0, &b[i1 + rank], &c__1);
#line 612 "NF01BR.f"
		}
#line 614 "NF01BR.f"
		dgemv_("NoTranspose", &st, &nths, &c_b56, &s[s_offset], lds, &
			b[1], &c__1, &c_b58, &b[i1], &c__1, (ftnlen)11);
#line 616 "NF01BR.f"
		i__1 = *ldr + 1;
#line 616 "NF01BR.f"
		dswap_(&st, &r__[i1 + (bsn + 1) * r_dim1], &i__1, &sdiag[i1], 
			&c__1);
#line 617 "NF01BR.f"
		dtrsv_("Lower", transl, "NonUnit", &rank, &r__[i1 + (bsn + 1) 
			* r_dim1], ldr, &b[i1], &c__1, (ftnlen)5, (ftnlen)1, (
			ftnlen)7);
#line 619 "NF01BR.f"
		i__1 = *ldr + 1;
#line 619 "NF01BR.f"
		dswap_(&st, &r__[i1 + (bsn + 1) * r_dim1], &i__1, &sdiag[i1], 
			&c__1);
#line 620 "NF01BR.f"
	    }

#line 622 "NF01BR.f"
	}

#line 624 "NF01BR.f"
    } else {

#line 626 "NF01BR.f"
	if (! tranr) {

/*           Solve R*x = b using back substitution. */

#line 630 "NF01BR.f"
	    i1 = nths + 1;
#line 631 "NF01BR.f"
	    if (st > 0) {
#line 632 "NF01BR.f"
		rank = ranks[l];
#line 633 "NF01BR.f"
		if (rank < st) {
#line 633 "NF01BR.f"
		    i__1 = st - rank;
#line 633 "NF01BR.f"
		    dcopy_(&i__1, dum, &c__0, &b[i1 + rank], &c__1);
#line 633 "NF01BR.f"
		}
#line 635 "NF01BR.f"
		dtrsv_("Upper", trans, "NonUnit", &rank, &r__[i1 + (bsn + 1) *
			 r_dim1], ldr, &b[i1], &c__1, (ftnlen)5, (ftnlen)1, (
			ftnlen)7);
#line 637 "NF01BR.f"
		dgemv_(trans, &nths, &st, &c_b56, &r__[(bsn + 1) * r_dim1 + 1]
			, ldr, &b[nths + 1], &c__1, &c_b58, &b[1], &c__1, (
			ftnlen)1);
#line 639 "NF01BR.f"
	    }

#line 641 "NF01BR.f"
	    for (k = bn; k >= 1; --k) {
#line 642 "NF01BR.f"
		i1 -= bsn;
#line 643 "NF01BR.f"
		rank = ranks[k];
#line 644 "NF01BR.f"
		if (rank < bsn) {
#line 644 "NF01BR.f"
		    i__1 = bsn - rank;
#line 644 "NF01BR.f"
		    dcopy_(&i__1, dum, &c__0, &b[i1 + rank], &c__1);
#line 644 "NF01BR.f"
		}
#line 646 "NF01BR.f"
		dtrsv_("Upper", trans, "NonUnit", &rank, &r__[i1 + r_dim1], 
			ldr, &b[i1], &c__1, (ftnlen)5, (ftnlen)1, (ftnlen)7);
#line 648 "NF01BR.f"
/* L160: */
#line 648 "NF01BR.f"
	    }

#line 650 "NF01BR.f"
	} else {

/*           Solve R'*x = b using forward substitution. */

#line 654 "NF01BR.f"
	    i1 = 1;

#line 656 "NF01BR.f"
	    i__1 = bn;
#line 656 "NF01BR.f"
	    for (k = 1; k <= i__1; ++k) {
#line 657 "NF01BR.f"
		rank = ranks[k];
#line 658 "NF01BR.f"
		if (rank < bsn) {
#line 658 "NF01BR.f"
		    i__2 = bsn - rank;
#line 658 "NF01BR.f"
		    dcopy_(&i__2, dum, &c__0, &b[i1 + rank], &c__1);
#line 658 "NF01BR.f"
		}
#line 660 "NF01BR.f"
		dtrsv_("Upper", trans, "NonUnit", &rank, &r__[i1 + r_dim1], 
			ldr, &b[i1], &c__1, (ftnlen)5, (ftnlen)1, (ftnlen)7);
#line 662 "NF01BR.f"
		i1 += bsn;
#line 663 "NF01BR.f"
/* L170: */
#line 663 "NF01BR.f"
	    }

#line 665 "NF01BR.f"
	    if (st > 0) {
#line 666 "NF01BR.f"
		rank = ranks[l];
#line 667 "NF01BR.f"
		if (rank < st) {
#line 667 "NF01BR.f"
		    i__1 = st - rank;
#line 667 "NF01BR.f"
		    dcopy_(&i__1, dum, &c__0, &b[i1 + rank], &c__1);
#line 667 "NF01BR.f"
		}
#line 669 "NF01BR.f"
		dgemv_(trans, &nths, &st, &c_b56, &r__[(bsn + 1) * r_dim1 + 1]
			, ldr, &b[1], &c__1, &c_b58, &b[i1], &c__1, (ftnlen)1)
			;
#line 671 "NF01BR.f"
		dtrsv_("Upper", trans, "NonUnit", &rank, &r__[i1 + (bsn + 1) *
			 r_dim1], ldr, &b[i1], &c__1, (ftnlen)5, (ftnlen)1, (
			ftnlen)7);
#line 673 "NF01BR.f"
	    }

#line 675 "NF01BR.f"
	}
#line 676 "NF01BR.f"
    }

#line 678 "NF01BR.f"
    if (econd && lower) {
#line 679 "NF01BR.f"
	i__ = 1;

/*        If COND = 'E' and UPLO = 'L', swap the diagonal elements of R */
/*        and the elements of SDIAG and swap back the upper and lower */
/*        triangular parts of R, including the part corresponding to S. */

#line 685 "NF01BR.f"
	i__1 = bn;
#line 685 "NF01BR.f"
	for (k = 1; k <= i__1; ++k) {
#line 686 "NF01BR.f"
	    i__2 = *ldr + 1;
#line 686 "NF01BR.f"
	    dswap_(&bsn, &r__[i__ + r_dim1], &i__2, &sdiag[i__], &c__1);

#line 688 "NF01BR.f"
	    i__2 = bsn;
#line 688 "NF01BR.f"
	    for (j = 1; j <= i__2; ++j) {
#line 689 "NF01BR.f"
		i__3 = bsn - j + 1;
#line 689 "NF01BR.f"
		dswap_(&i__3, &r__[i__ + j * r_dim1], ldr, &r__[i__ + j * 
			r_dim1], &c__1);
#line 690 "NF01BR.f"
		++i__;
#line 691 "NF01BR.f"
/* L180: */
#line 691 "NF01BR.f"
	    }

#line 693 "NF01BR.f"
/* L190: */
#line 693 "NF01BR.f"
	}

#line 695 "NF01BR.f"
	if (st > 0) {
#line 696 "NF01BR.f"
	    i__1 = *ldr + 1;
#line 696 "NF01BR.f"
	    dswap_(&st, &r__[i__ + (bsn + 1) * r_dim1], &i__1, &sdiag[i__], &
		    c__1);

#line 698 "NF01BR.f"
	    i__1 = nc;
#line 698 "NF01BR.f"
	    for (j = bsn + 1; j <= i__1; ++j) {
#line 699 "NF01BR.f"
		dswap_(&nths, &r__[j * r_dim1 + 1], &c__1, &s[j - bsn + 
			s_dim1], lds);
#line 700 "NF01BR.f"
		i__2 = nc - j + 1;
#line 700 "NF01BR.f"
		dswap_(&i__2, &r__[i__ + j * r_dim1], ldr, &r__[i__ + j * 
			r_dim1], &c__1);
#line 701 "NF01BR.f"
		++i__;
#line 702 "NF01BR.f"
/* L200: */
#line 702 "NF01BR.f"
	    }

#line 704 "NF01BR.f"
	}

#line 706 "NF01BR.f"
    }

#line 708 "NF01BR.f"
    return 0;

/* *** Last line of NF01BR *** */
} /* nf01br_ */

