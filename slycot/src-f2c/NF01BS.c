#line 1 "NF01BS.f"
/* NF01BS.f -- translated by f2c (version 20100827).
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

#line 1 "NF01BS.f"
/* Table of constant values */

static integer c__1 = 1;
static logical c_true = TRUE_;

/* Subroutine */ int nf01bs_(integer *n, integer *ipar, integer *lipar, 
	doublereal *fnorm, doublereal *j, integer *ldj, doublereal *e, 
	doublereal *jnorms, doublereal *gnorm, integer *ipvt, doublereal *
	dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, k, l, m, bn, jl, st, bsm, bsn, jlm, mmn;
    static doublereal sum;
    static integer ibsm, ibsn;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer itau, nths;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int md03bx_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    static integer ibsni;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static integer jwork;
    extern /* Subroutine */ int dgeqp3_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen), dlapmt_(logical *, integer *, integer *, 
	    doublereal *, integer *, integer *), dormqr_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
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

/*     To compute the QR factorization of the Jacobian matrix J, as */
/*     received in compressed form from SLICOT Library routine NF01BD, */

/*            /  dy(1)/dwb(1)  |  dy(1)/ dtheta  \ */
/*       Jc = |       :        |       :         | , */
/*            \  dy(L)/dwb(L)  |  dy(L)/ dtheta  / */

/*     and to apply the transformation Q on the error vector e (in-situ). */
/*     The factorization is J*P = Q*R, where Q is a matrix with */
/*     orthogonal columns, P a permutation matrix, and R an upper */
/*     trapezoidal matrix with diagonal elements of nonincreasing */
/*     magnitude for each block column (see below). The 1-norm of the */
/*     scaled gradient is also returned. */

/*     Actually, the Jacobian J has the block form */

/*       dy(1)/dwb(1)       0         .....       0        dy(1)/dtheta */
/*            0        dy(2)/dwb(2)   .....       0        dy(2)/dtheta */
/*          .....         .....       .....     .....         ..... */
/*            0           .....         0    dy(L)/dwb(L)  dy(L)/dtheta */

/*     but the zero blocks are omitted. The diagonal blocks have the */
/*     same size and correspond to the nonlinear part. The last block */
/*     column corresponds to the linear part. It is assumed that the */
/*     Jacobian matrix has at least as many rows as columns. The linear */
/*     or nonlinear parts can be empty. If L <= 1, the Jacobian is */
/*     represented as a full matrix. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of columns of the Jacobian matrix J. */
/*             N = BN*BSN + ST >= 0.  (See parameter description below.) */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             The integer parameters describing the structure of the */
/*             matrix J, as follows: */
/*             IPAR(1) must contain ST, the number of parameters */
/*                     corresponding to the linear part.  ST >= 0. */
/*             IPAR(2) must contain BN, the number of blocks, BN = L, */
/*                     for the parameters corresponding to the nonlinear */
/*                     part.  BN >= 0. */
/*             IPAR(3) must contain BSM, the number of rows of the blocks */
/*                     J_k = dy(k)/dwb(k), k = 1:BN, if BN > 0, or the */
/*                     number of rows of the matrix J, if BN <= 1. */
/*                     BN*BSM >= N, if BN > 0; */
/*                     BSM >= N,    if BN = 0. */
/*             IPAR(4) must contain BSN, the number of columns of the */
/*                     blocks J_k, k = 1:BN.  BSN >= 0. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 4. */

/*     FNORM   (input) DOUBLE PRECISION */
/*             The Euclidean norm of the vector e.  FNORM >= 0. */

/*     J       (input/output) DOUBLE PRECISION array, dimension (LDJ, NC) */
/*             where NC = N if BN <= 1, and NC = BSN+ST, if BN > 1. */
/*             On entry, the leading NR-by-NC part of this array must */
/*             contain the (compressed) representation (Jc) of the */
/*             Jacobian matrix J, where NR = BSM if BN <= 1, and */
/*             NR = BN*BSM, if BN > 1. */
/*             On exit, the leading N-by-NC part of this array contains */
/*             a (compressed) representation of the upper triangular */
/*             factor R of the Jacobian matrix. The matrix R has the same */
/*             structure as the Jacobian matrix J, but with an additional */
/*             diagonal block. Note that for efficiency of the later */
/*             calculations, the matrix R is delivered with the leading */
/*             dimension MAX(1,N), possibly much smaller than the value */
/*             of LDJ on entry. */

/*     LDJ     (input/output) INTEGER */
/*             The leading dimension of array J. */
/*             On entry, LDJ >= MAX(1,NR). */
/*             On exit,  LDJ >= MAX(1,N). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (NR) */
/*             On entry, this array contains the vector e, */
/*             e = vec( Y - y ), where Y is set of output samples, and */
/*             vec denotes the concatenation of the columns of a matrix. */
/*             On exit, this array contains the updated vector Z*Q'*e, */
/*             where Z is the block row permutation matrix used in the */
/*             QR factorization of J (see METHOD). */

/*     JNORMS  (output) DOUBLE PRECISION array, dimension (N) */
/*             This array contains the Euclidean norms of the columns */
/*             of the Jacobian matrix, considered in the initial order. */

/*     GNORM   (output) DOUBLE PRECISION */
/*             If FNORM > 0, the 1-norm of the scaled vector J'*e/FNORM, */
/*             with each element i further divided by JNORMS(i) (if */
/*             JNORMS(i) is nonzero). */
/*             If FNORM = 0, the returned value of GNORM is 0. */

/*     IPVT    (output) INTEGER array, dimension (N) */
/*             This array defines the permutation matrix P such that */
/*             J*P = Q*R. Column j of P is column IPVT(j) of the identity */
/*             matrix. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 1,      if N = 0 or BN <= 1 and BSM = N = 1; */
/*                               otherwise, */
/*             LDWORK >= 4*N+1,  if BN <= 1 or  BSN = 0; */
/*             LDWORK >= JWORK,  if BN >  1 and BSN > 0, where JWORK is */
/*                               given by the following procedure: */
/*              JWORK  = BSN + MAX(3*BSN+1,ST); */
/*              JWORK  = MAX(JWORK,4*ST+1),         if BSM > BSN; */
/*              JWORK  = MAX(JWORK,(BSM-BSN)*(BN-1)), */
/*                                                  if BSN < BSM < 2*BSN. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     A QR factorization with column pivoting of the matrix J is */
/*     computed, J*P = Q*R. */

/*     If l = L > 1, the R factor of the QR factorization has the same */
/*     structure as the Jacobian, but with an additional diagonal block. */
/*     Denote */

/*         /   J_1    0    ..   0   |  L_1  \ */
/*         |    0    J_2   ..   0   |  L_2  | */
/*     J = |    :     :    ..   :   |   :   | . */
/*         |    :     :    ..   :   |   :   | */
/*         \    0     0    ..  J_l  |  L_l  / */

/*     The algorithm consists in two phases. In the first phase, the */
/*     algorithm uses QR factorizations with column pivoting for each */
/*     block J_k, k = 1:l, and applies the orthogonal matrix Q'_k to the */
/*     corresponding part of the last block column and of e. After all */
/*     block rows have been processed, the block rows are interchanged */
/*     so that the zeroed submatrices in the first l block columns are */
/*     moved to the bottom part. The same block row permutation Z is */
/*     also applied to the vector e. At the end of the first phase, */
/*     the structure of the processed matrix J is */

/*         /   R_1    0    ..   0   |  L^1_1  \ */
/*         |    0    R_2   ..   0   |  L^1_2  | */
/*         |    :     :    ..   :   |    :    | . */
/*         |    :     :    ..   :   |    :    | */
/*         |    0     0    ..  R_l  |  L^1_l  | */
/*         |    0     0    ..   0   |  L^2_1  | */
/*         |    :     :    ..   :   |    :    | */
/*         \    0     0    ..   0   |  L^2_l  / */

/*     In the second phase, the submatrix L^2_1:l is triangularized */
/*     using an additional QR factorization with pivoting. (The columns */
/*     of L^1_1:l are also permuted accordingly.) Therefore, the column */
/*     pivoting is restricted to each such local block column. */

/*     If l <= 1, the matrix J is triangularized in one phase, by one */
/*     QR factorization with pivoting. In this case, the column */
/*     pivoting is global. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     Feb. 22, 2004. */

/*     KEYWORDS */

/*     Elementary matrix operations, Jacobian matrix, matrix algebra, */
/*     matrix operations, Wiener system. */

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

#line 233 "NF01BS.f"
    /* Parameter adjustments */
#line 233 "NF01BS.f"
    --dwork;
#line 233 "NF01BS.f"
    --ipvt;
#line 233 "NF01BS.f"
    --jnorms;
#line 233 "NF01BS.f"
    --e;
#line 233 "NF01BS.f"
    --j;
#line 233 "NF01BS.f"
    --ipar;
#line 233 "NF01BS.f"

#line 233 "NF01BS.f"
    /* Function Body */
#line 233 "NF01BS.f"
    *info = 0;
#line 234 "NF01BS.f"
    if (*n < 0) {
#line 235 "NF01BS.f"
	*info = -1;
#line 236 "NF01BS.f"
    } else if (*lipar < 4) {
#line 237 "NF01BS.f"
	*info = -3;
#line 238 "NF01BS.f"
    } else if (*fnorm < 0.) {
#line 239 "NF01BS.f"
	*info = -4;
#line 240 "NF01BS.f"
    } else if (*ldj < max(1,*n)) {
#line 241 "NF01BS.f"
	*info = -6;
#line 242 "NF01BS.f"
    } else {
#line 243 "NF01BS.f"
	st = ipar[1];
#line 244 "NF01BS.f"
	bn = ipar[2];
#line 245 "NF01BS.f"
	bsm = ipar[3];
#line 246 "NF01BS.f"
	bsn = ipar[4];
#line 247 "NF01BS.f"
	nths = bn * bsn;
#line 248 "NF01BS.f"
	mmn = bsm - bsn;
#line 249 "NF01BS.f"
	if (bn > 0) {
#line 250 "NF01BS.f"
	    m = bn * bsm;
#line 251 "NF01BS.f"
	} else {
#line 252 "NF01BS.f"
	    m = *n;
#line 253 "NF01BS.f"
	}
/* Computing MIN */
#line 254 "NF01BS.f"
	i__1 = min(st,bn), i__1 = min(i__1,bsm);
#line 254 "NF01BS.f"
	if (min(i__1,bsn) < 0) {
#line 255 "NF01BS.f"
	    *info = -2;
#line 256 "NF01BS.f"
	} else if (*n != nths + st) {
#line 257 "NF01BS.f"
	    *info = -1;
#line 258 "NF01BS.f"
	} else if (m < *n) {
#line 259 "NF01BS.f"
	    *info = -2;
#line 260 "NF01BS.f"
	} else if (*ldj < max(1,m)) {
#line 261 "NF01BS.f"
	    *info = -6;
#line 262 "NF01BS.f"
	} else {
#line 263 "NF01BS.f"
	    if (*n == 0) {
#line 264 "NF01BS.f"
		jwork = 1;
#line 265 "NF01BS.f"
	    } else if (bn <= 1 || bsn == 0) {
#line 266 "NF01BS.f"
		if (bn <= 1 && bsm == 1 && *n == 1) {
#line 267 "NF01BS.f"
		    jwork = 1;
#line 268 "NF01BS.f"
		} else {
#line 269 "NF01BS.f"
		    jwork = (*n << 2) + 1;
#line 270 "NF01BS.f"
		}
#line 271 "NF01BS.f"
	    } else {
/* Computing MAX */
#line 272 "NF01BS.f"
		i__1 = bsn * 3 + 1;
#line 272 "NF01BS.f"
		jwork = bsn + max(i__1,st);
#line 273 "NF01BS.f"
		if (bsm > bsn) {
/* Computing MAX */
#line 274 "NF01BS.f"
		    i__1 = jwork, i__2 = (st << 2) + 1;
#line 274 "NF01BS.f"
		    jwork = max(i__1,i__2);
#line 275 "NF01BS.f"
		    if (bsm < bsn << 1) {
/* Computing MAX */
#line 275 "NF01BS.f"
			i__1 = jwork, i__2 = mmn * (bn - 1);
#line 275 "NF01BS.f"
			jwork = max(i__1,i__2);
#line 275 "NF01BS.f"
		    }
#line 277 "NF01BS.f"
		}
#line 278 "NF01BS.f"
	    }
#line 279 "NF01BS.f"
	    if (*ldwork < jwork) {
#line 279 "NF01BS.f"
		*info = -12;
#line 279 "NF01BS.f"
	    }
#line 281 "NF01BS.f"
	}
#line 282 "NF01BS.f"
    }

#line 284 "NF01BS.f"
    if (*info != 0) {

/*        Error return. */

#line 288 "NF01BS.f"
	i__1 = -(*info);
#line 288 "NF01BS.f"
	xerbla_("NF01BS", &i__1, (ftnlen)6);
#line 289 "NF01BS.f"
	return 0;
#line 290 "NF01BS.f"
    }

/*     Quick return if possible. */

#line 294 "NF01BS.f"
    *gnorm = 0.;
#line 295 "NF01BS.f"
    if (*n == 0) {
#line 296 "NF01BS.f"
	*ldj = 1;
#line 297 "NF01BS.f"
	dwork[1] = 1.;
#line 298 "NF01BS.f"
	return 0;
#line 299 "NF01BS.f"
    }

#line 301 "NF01BS.f"
    if (bn <= 1 || bsn == 0) {

/*        Special case, l <= 1 or BSN = 0: the Jacobian is represented */
/*        as a full matrix. */
/*        (Note: Comments in the code beginning "Workspace:" describe the */
/*        minimal amount of real workspace needed at that point in the */
/*        code, as well as the preferred amount for good performance. */
/*        NB refers to the optimal block size for the immediately */
/*        following subroutine, as returned by ILAENV.) */

/*        Workspace: need:    4*N + 1; */
/*                   prefer:  3*N + ( N+1 )*NB. */

#line 314 "NF01BS.f"
	md03bx_(&m, n, fnorm, &j[1], ldj, &e[1], &jnorms[1], gnorm, &ipvt[1], 
		&dwork[1], ldwork, info);
#line 316 "NF01BS.f"
	return 0;
#line 317 "NF01BS.f"
    }

/*     General case: l > 1 and BSN > 0. */
/*     Initialize the column pivoting indices. */

#line 322 "NF01BS.f"
    i__1 = *n;
#line 322 "NF01BS.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 323 "NF01BS.f"
	ipvt[i__] = 0;
#line 324 "NF01BS.f"
/* L10: */
#line 324 "NF01BS.f"
    }

/*     Compute the QR factorization with pivoting of J. */
/*     Pivoting is done separately on each block column of J. */

#line 329 "NF01BS.f"
    wrkopt = 1;
#line 330 "NF01BS.f"
    ibsn = 1;
#line 331 "NF01BS.f"
    jl = *ldj * bsn + 1;
#line 332 "NF01BS.f"
    jwork = bsn + 1;

#line 334 "NF01BS.f"
    i__1 = m;
#line 334 "NF01BS.f"
    i__2 = bsm;
#line 334 "NF01BS.f"
    for (ibsm = 1; i__2 < 0 ? ibsm >= i__1 : ibsm <= i__1; ibsm += i__2) {

/*        Compute the QR factorization with pivoting of J_k, and apply Q' */
/*        to the corresponding part of the last block-column and of e. */
/*        Workspace: need:    4*BSN + 1; */
/*                   prefer:  3*BSN + ( BSN+1 )*NB. */

#line 341 "NF01BS.f"
	i__3 = *ldwork - jwork + 1;
#line 341 "NF01BS.f"
	dgeqp3_(&bsm, &bsn, &j[ibsm], ldj, &ipvt[ibsn], &dwork[1], &dwork[
		jwork], &i__3, info);
/* Computing MAX */
#line 343 "NF01BS.f"
	i__3 = wrkopt, i__4 = (integer) dwork[jwork] + jwork - 1;
#line 343 "NF01BS.f"
	wrkopt = max(i__3,i__4);
#line 344 "NF01BS.f"
	if (ibsm > 1) {

/*           Adjust the column pivoting indices. */

#line 348 "NF01BS.f"
	    i__3 = ibsn + bsn - 1;
#line 348 "NF01BS.f"
	    for (i__ = ibsn; i__ <= i__3; ++i__) {
#line 349 "NF01BS.f"
		ipvt[i__] = ipvt[i__] + ibsn - 1;
#line 350 "NF01BS.f"
/* L20: */
#line 350 "NF01BS.f"
	    }

#line 352 "NF01BS.f"
	}

#line 354 "NF01BS.f"
	if (st > 0) {

/*           Workspace: need:    BSN + ST; */
/*                      prefer:  BSN + ST*NB. */

#line 359 "NF01BS.f"
	    i__3 = *ldwork - jwork + 1;
#line 359 "NF01BS.f"
	    dormqr_("Left", "Transpose", &bsm, &st, &bsn, &j[ibsm], ldj, &
		    dwork[1], &j[jl], ldj, &dwork[jwork], &i__3, info, (
		    ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 362 "NF01BS.f"
	    i__3 = wrkopt, i__4 = (integer) dwork[jwork] + jwork - 1;
#line 362 "NF01BS.f"
	    wrkopt = max(i__3,i__4);
#line 363 "NF01BS.f"
	}

/*        Workspace: need:    BSN + 1; */
/*                   prefer:  BSN + NB. */

#line 368 "NF01BS.f"
	i__3 = *ldwork - jwork + 1;
#line 368 "NF01BS.f"
	dormqr_("Left", "Transpose", &bsm, &c__1, &bsn, &j[ibsm], ldj, &dwork[
		1], &e[ibsm], &bsm, &dwork[jwork], &i__3, info, (ftnlen)4, (
		ftnlen)9);
#line 371 "NF01BS.f"
	jl += bsm;
#line 372 "NF01BS.f"
	ibsn += bsn;
#line 373 "NF01BS.f"
/* L30: */
#line 373 "NF01BS.f"
    }

#line 375 "NF01BS.f"
    if (mmn > 0) {

/*        Case BSM > BSN. */
/*        Compute the original column norms for the first block column */
/*        of Jc. */
/*        Permute the rows of the first block column to move the zeroed */
/*        submatrices to the bottom. In the same loops, reshape the */
/*        first block column of R to have the leading dimension N. */

#line 384 "NF01BS.f"
	l = ipvt[1];
#line 385 "NF01BS.f"
	jnorms[l] = abs(j[1]);
#line 386 "NF01BS.f"
	ibsm = bsm + 1;
#line 387 "NF01BS.f"
	ibsn = bsn + 1;

#line 389 "NF01BS.f"
	i__2 = bn - 1;
#line 389 "NF01BS.f"
	for (k = 1; k <= i__2; ++k) {
#line 390 "NF01BS.f"
	    j[ibsn] = j[ibsm];
#line 391 "NF01BS.f"
	    l = ipvt[ibsn];
#line 392 "NF01BS.f"
	    jnorms[l] = (d__1 = j[ibsn], abs(d__1));
#line 393 "NF01BS.f"
	    ibsm += bsm;
#line 394 "NF01BS.f"
	    ibsn += bsn;
#line 395 "NF01BS.f"
/* L40: */
#line 395 "NF01BS.f"
	}

#line 397 "NF01BS.f"
	ibsn += st;

#line 399 "NF01BS.f"
	i__2 = bsn;
#line 399 "NF01BS.f"
	for (i__ = 2; i__ <= i__2; ++i__) {
#line 400 "NF01BS.f"
	    ibsm = (i__ - 1) * *ldj + 1;
#line 401 "NF01BS.f"
	    jl = i__;

#line 403 "NF01BS.f"
	    i__1 = bn;
#line 403 "NF01BS.f"
	    for (k = 1; k <= i__1; ++k) {

#line 405 "NF01BS.f"
		i__3 = i__ - 1;
#line 405 "NF01BS.f"
		for (l = 0; l <= i__3; ++l) {
#line 406 "NF01BS.f"
		    j[ibsn + l] = j[ibsm + l];
#line 407 "NF01BS.f"
/* L45: */
#line 407 "NF01BS.f"
		}

#line 409 "NF01BS.f"
		l = ipvt[jl];
#line 410 "NF01BS.f"
		jnorms[l] = dnrm2_(&i__, &j[ibsn], &c__1);
#line 411 "NF01BS.f"
		ibsm += bsm;
#line 412 "NF01BS.f"
		ibsn += bsn;
#line 413 "NF01BS.f"
		jl += bsn;
#line 414 "NF01BS.f"
/* L50: */
#line 414 "NF01BS.f"
	    }

#line 416 "NF01BS.f"
	    ibsn += st;
#line 417 "NF01BS.f"
/* L60: */
#line 417 "NF01BS.f"
	}

/*        Permute the rows of the second block column of Jc and of */
/*        the vector e. */

#line 422 "NF01BS.f"
	jl = *ldj * bsn;
#line 423 "NF01BS.f"
	if (bsm >= bsn << 1) {

/*           A swap operation can be used. */

#line 427 "NF01BS.f"
	    i__2 = st;
#line 427 "NF01BS.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 428 "NF01BS.f"
		ibsn = bsn + 1;

#line 430 "NF01BS.f"
		i__1 = m;
#line 430 "NF01BS.f"
		i__3 = bsm;
#line 430 "NF01BS.f"
		for (ibsm = bsm + 1; i__3 < 0 ? ibsm >= i__1 : ibsm <= i__1; 
			ibsm += i__3) {
#line 431 "NF01BS.f"
		    dswap_(&mmn, &j[jl + ibsm], &c__1, &j[jl + ibsn], &c__1);
#line 432 "NF01BS.f"
		    ibsn += bsn;
#line 433 "NF01BS.f"
/* L70: */
#line 433 "NF01BS.f"
		}

#line 435 "NF01BS.f"
		jl += *ldj;
#line 436 "NF01BS.f"
/* L80: */
#line 436 "NF01BS.f"
	    }

/*           Permute the rows of e. */

#line 440 "NF01BS.f"
	    ibsn = bsn + 1;

#line 442 "NF01BS.f"
	    i__2 = m;
#line 442 "NF01BS.f"
	    i__3 = bsm;
#line 442 "NF01BS.f"
	    for (ibsm = bsm + 1; i__3 < 0 ? ibsm >= i__2 : ibsm <= i__2; ibsm 
		    += i__3) {
#line 443 "NF01BS.f"
		dswap_(&mmn, &e[ibsm], &c__1, &e[ibsn], &c__1);
#line 444 "NF01BS.f"
		ibsn += bsn;
#line 445 "NF01BS.f"
/* L90: */
#line 445 "NF01BS.f"
	    }

#line 447 "NF01BS.f"
	} else {

/*           A swap operation cannot be used. */
/*           Workspace: need:    ( BSM-BSN )*( BN-1 ). */

#line 452 "NF01BS.f"
	    i__3 = st;
#line 452 "NF01BS.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 453 "NF01BS.f"
		ibsn = bsn + 1;
#line 454 "NF01BS.f"
		jlm = jl + ibsn;
#line 455 "NF01BS.f"
		jwork = 1;

#line 457 "NF01BS.f"
		i__2 = m;
#line 457 "NF01BS.f"
		i__1 = bsm;
#line 457 "NF01BS.f"
		for (ibsm = bsm + 1; i__1 < 0 ? ibsm >= i__2 : ibsm <= i__2; 
			ibsm += i__1) {
#line 458 "NF01BS.f"
		    dcopy_(&mmn, &j[jlm], &c__1, &dwork[jwork], &c__1);

#line 460 "NF01BS.f"
		    i__4 = jl + bsn - 1;
#line 460 "NF01BS.f"
		    for (k = jl; k <= i__4; ++k) {
#line 461 "NF01BS.f"
			j[ibsn + k] = j[ibsm + k];
#line 462 "NF01BS.f"
/* L105: */
#line 462 "NF01BS.f"
		    }

#line 464 "NF01BS.f"
		    jlm += bsm;
#line 465 "NF01BS.f"
		    ibsn += bsn;
#line 466 "NF01BS.f"
		    jwork += mmn;
#line 467 "NF01BS.f"
/* L100: */
#line 467 "NF01BS.f"
		}

#line 469 "NF01BS.f"
		i__1 = mmn * (bn - 1);
#line 469 "NF01BS.f"
		dcopy_(&i__1, &dwork[1], &c__1, &j[jl + ibsn], &c__1);
#line 470 "NF01BS.f"
		jl += *ldj;
#line 471 "NF01BS.f"
/* L110: */
#line 471 "NF01BS.f"
	    }

/*           Permute the rows of e. */

#line 475 "NF01BS.f"
	    ibsn = bsn + 1;
#line 476 "NF01BS.f"
	    jlm = ibsn;
#line 477 "NF01BS.f"
	    jwork = 1;

#line 479 "NF01BS.f"
	    i__3 = m;
#line 479 "NF01BS.f"
	    i__1 = bsm;
#line 479 "NF01BS.f"
	    for (ibsm = bsm + 1; i__1 < 0 ? ibsm >= i__3 : ibsm <= i__3; ibsm 
		    += i__1) {
#line 480 "NF01BS.f"
		dcopy_(&mmn, &e[jlm], &c__1, &dwork[jwork], &c__1);

#line 482 "NF01BS.f"
		i__2 = bsn - 1;
#line 482 "NF01BS.f"
		for (k = 0; k <= i__2; ++k) {
#line 483 "NF01BS.f"
		    e[ibsn + k] = e[ibsm + k];
#line 484 "NF01BS.f"
/* L115: */
#line 484 "NF01BS.f"
		}

#line 486 "NF01BS.f"
		jlm += bsm;
#line 487 "NF01BS.f"
		ibsn += bsn;
#line 488 "NF01BS.f"
		jwork += mmn;
#line 489 "NF01BS.f"
/* L120: */
#line 489 "NF01BS.f"
	    }

#line 491 "NF01BS.f"
	    i__1 = mmn * (bn - 1);
#line 491 "NF01BS.f"
	    dcopy_(&i__1, &dwork[1], &c__1, &e[ibsn], &c__1);
#line 492 "NF01BS.f"
	}

#line 494 "NF01BS.f"
	if (st > 0) {

/*           Compute the QR factorization with pivoting of the submatrix */
/*           L^2_1:l, and apply Q' to the corresponding part of e. */

/*           Workspace: need:    4*ST + 1; */
/*                      prefer:  3*ST + ( ST+1 )*NB. */

#line 502 "NF01BS.f"
	    jl = (*ldj + bn) * bsn + 1;
#line 503 "NF01BS.f"
	    itau = 1;
#line 504 "NF01BS.f"
	    jwork = itau + st;
#line 505 "NF01BS.f"
	    i__1 = mmn * bn;
#line 505 "NF01BS.f"
	    i__3 = *ldwork - jwork + 1;
#line 505 "NF01BS.f"
	    dgeqp3_(&i__1, &st, &j[jl], ldj, &ipvt[nths + 1], &dwork[itau], &
		    dwork[jwork], &i__3, info);
/* Computing MAX */
#line 508 "NF01BS.f"
	    i__1 = wrkopt, i__3 = (integer) dwork[jwork] + jwork - 1;
#line 508 "NF01BS.f"
	    wrkopt = max(i__1,i__3);

/*           Permute columns of the upper part of the second block */
/*           column of Jc. */

#line 513 "NF01BS.f"
	    dlapmt_(&c_true, &nths, &st, &j[jl - nths], ldj, &ipvt[nths + 1]);

/*           Adjust the column pivoting indices. */

#line 518 "NF01BS.f"
	    i__1 = *n;
#line 518 "NF01BS.f"
	    for (i__ = nths + 1; i__ <= i__1; ++i__) {
#line 519 "NF01BS.f"
		ipvt[i__] += nths;
#line 520 "NF01BS.f"
/* L130: */
#line 520 "NF01BS.f"
	    }

/*           Workspace: need:    ST + 1; */
/*                      prefer:  ST + NB. */

#line 525 "NF01BS.f"
	    i__1 = mmn * bn;
#line 525 "NF01BS.f"
	    i__3 = *ldwork - jwork + 1;
#line 525 "NF01BS.f"
	    dormqr_("Left", "Transpose", &i__1, &c__1, &st, &j[jl], ldj, &
		    dwork[itau], &e[ibsn], ldj, &dwork[jwork], &i__3, info, (
		    ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 528 "NF01BS.f"
	    i__1 = wrkopt, i__3 = (integer) dwork[jwork] + jwork - 1;
#line 528 "NF01BS.f"
	    wrkopt = max(i__1,i__3);

/*           Reshape the second block column of R to have the leading */
/*           dimension N. */

#line 533 "NF01BS.f"
	    ibsn = *n * bsn + 1;
#line 534 "NF01BS.f"
	    dlacpy_("Full", n, &st, &j[*ldj * bsn + 1], ldj, &j[ibsn], n, (
		    ftnlen)4);

/*           Compute the original column norms for the second block */
/*           column. */

#line 539 "NF01BS.f"
	    i__1 = *n;
#line 539 "NF01BS.f"
	    for (i__ = nths + 1; i__ <= i__1; ++i__) {
#line 540 "NF01BS.f"
		l = ipvt[i__];
#line 541 "NF01BS.f"
		jnorms[l] = dnrm2_(&i__, &j[ibsn], &c__1);
#line 542 "NF01BS.f"
		ibsn += *n;
#line 543 "NF01BS.f"
/* L140: */
#line 543 "NF01BS.f"
	    }

#line 545 "NF01BS.f"
	}

#line 547 "NF01BS.f"
    } else {

/*        Case BSM = BSN. */
/*        Compute the original column norms for the first block column */
/*        of Jc. */

#line 553 "NF01BS.f"
	ibsn = 1;

#line 555 "NF01BS.f"
	i__1 = bsn;
#line 555 "NF01BS.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 556 "NF01BS.f"
	    jl = i__;

#line 558 "NF01BS.f"
	    i__3 = bn;
#line 558 "NF01BS.f"
	    for (k = 1; k <= i__3; ++k) {
#line 559 "NF01BS.f"
		l = ipvt[jl];
#line 560 "NF01BS.f"
		jnorms[l] = dnrm2_(&i__, &j[ibsn], &c__1);
#line 561 "NF01BS.f"
		ibsn += bsn;
#line 562 "NF01BS.f"
		jl += bsn;
#line 563 "NF01BS.f"
/* L150: */
#line 563 "NF01BS.f"
	    }

#line 565 "NF01BS.f"
	    ibsn += st;
#line 566 "NF01BS.f"
/* L160: */
#line 566 "NF01BS.f"
	}

#line 568 "NF01BS.f"
	i__1 = *n;
#line 568 "NF01BS.f"
	for (i__ = nths + 1; i__ <= i__1; ++i__) {
#line 569 "NF01BS.f"
	    ipvt[i__] = i__;
#line 570 "NF01BS.f"
/* L170: */
#line 570 "NF01BS.f"
	}

#line 572 "NF01BS.f"
    }

/*     Compute the norm of the scaled gradient. */

#line 576 "NF01BS.f"
    if (*fnorm != 0.) {

#line 578 "NF01BS.f"
	i__1 = nths;
#line 578 "NF01BS.f"
	i__3 = bsn;
#line 578 "NF01BS.f"
	for (ibsn = 1; i__3 < 0 ? ibsn >= i__1 : ibsn <= i__1; ibsn += i__3) {
#line 579 "NF01BS.f"
	    ibsni = ibsn;

#line 581 "NF01BS.f"
	    i__2 = bsn;
#line 581 "NF01BS.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 582 "NF01BS.f"
		l = ipvt[ibsn + i__ - 1];
#line 583 "NF01BS.f"
		if (jnorms[l] != 0.) {
#line 584 "NF01BS.f"
		    sum = ddot_(&i__, &j[ibsni], &c__1, &e[ibsn], &c__1) / *
			    fnorm;
/* Computing MAX */
#line 585 "NF01BS.f"
		    d__2 = *gnorm, d__3 = (d__1 = sum / jnorms[l], abs(d__1));
#line 585 "NF01BS.f"
		    *gnorm = max(d__2,d__3);
#line 586 "NF01BS.f"
		}
#line 587 "NF01BS.f"
		ibsni += *n;
#line 588 "NF01BS.f"
/* L180: */
#line 588 "NF01BS.f"
	    }

#line 590 "NF01BS.f"
/* L190: */
#line 590 "NF01BS.f"
	}

#line 592 "NF01BS.f"
	ibsni = *n * bsn + 1;

#line 594 "NF01BS.f"
	i__3 = *n;
#line 594 "NF01BS.f"
	for (i__ = nths + 1; i__ <= i__3; ++i__) {
#line 595 "NF01BS.f"
	    l = ipvt[i__];
#line 596 "NF01BS.f"
	    if (jnorms[l] != 0.) {
#line 597 "NF01BS.f"
		sum = ddot_(&i__, &j[ibsni], &c__1, &e[1], &c__1) / *fnorm;
/* Computing MAX */
#line 598 "NF01BS.f"
		d__2 = *gnorm, d__3 = (d__1 = sum / jnorms[l], abs(d__1));
#line 598 "NF01BS.f"
		*gnorm = max(d__2,d__3);
#line 599 "NF01BS.f"
	    }
#line 600 "NF01BS.f"
	    ibsni += *n;
#line 601 "NF01BS.f"
/* L200: */
#line 601 "NF01BS.f"
	}

#line 603 "NF01BS.f"
    }

#line 605 "NF01BS.f"
    *ldj = *n;
#line 606 "NF01BS.f"
    dwork[1] = (doublereal) wrkopt;
#line 607 "NF01BS.f"
    return 0;

/* *** Last line of NF01BS *** */
} /* nf01bs_ */

