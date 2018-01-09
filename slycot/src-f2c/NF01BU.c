#line 1 "NF01BU.f"
/* NF01BU.f -- translated by f2c (version 20100827).
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

#line 1 "NF01BU.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b9 = 0.;
static doublereal c_b11 = 1.;
static integer c__0 = 0;

/* Subroutine */ int nf01bu_(char *stor, char *uplo, integer *n, integer *
	ipar, integer *lipar, doublereal *dpar, integer *ldpar, doublereal *j,
	 integer *ldj, doublereal *jtj, integer *ldjtj, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen stor_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer j_dim1, j_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static doublereal c__;
    static integer k, m, i1, bn, ii, jl, st, bsm, bsn;
    static doublereal tmp[1];
    static integer ibsm, ibsn, nbsn;
    static logical full;
    static integer itmp[1], nths;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     nf01bv_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static logical upper;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen), dlaset_(char *, integer *, integer *,
	     doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);


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

/*     To compute the matrix J'*J + c*I, for the Jacobian J as received */
/*     from SLICOT Library routine NF01BD: */

/*          /  dy(1)/dwb(1)  |  dy(1)/dtheta  \ */
/*     Jc = |       :        |       :        | . */
/*          \  dy(L)/dwb(L)  |  dy(L)/dtheta  / */

/*     This is a compressed representation of the actual structure */

/*         /   J_1    0    ..   0   |  L_1  \ */
/*         |    0    J_2   ..   0   |  L_2  | */
/*     J = |    :     :    ..   :   |   :   | . */
/*         |    :     :    ..   :   |   :   | */
/*         \    0     0    ..  J_L  |  L_L  / */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     STOR    CHARACTER*1 */
/*             Specifies the storage scheme for the symmetric */
/*             matrix J'*J + c*I, as follows: */
/*             = 'F' :  full storage is used; */
/*             = 'P' :  packed storage is used. */

/*     UPLO    CHARACTER*1 */
/*             Specifies which part of the matrix J'*J + c*I is stored, */
/*             as follows: */
/*             = 'U' :  the upper triagular part is stored; */
/*             = 'L' :  the lower triagular part is stored. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix J'*J + c*I. */
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
/*             IPAR(4) must contain BSN, the number of columns of the */
/*                     blocks J_k, k = 1:BN.  BSN >= 0. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 4. */

/*     DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR) */
/*             The real parameters needed for solving the problem. */
/*             The entry DPAR(1) must contain the real scalar c. */

/*     LDPAR   (input) INTEGER */
/*             The length of the array DPAR.  LDPAR >= 1. */

/*     J       (input) DOUBLE PRECISION array, dimension (LDJ, NC) */
/*             where NC = N if BN <= 1, and NC = BSN+ST, if BN > 1. */
/*             The leading NR-by-NC part of this array must contain */
/*             the (compressed) representation (Jc) of the Jacobian */
/*             matrix J, where NR = BSM if BN <= 1, and NR = BN*BSM, */
/*             if BN > 1. */

/*     LDJ     (input) INTEGER */
/*             The leading dimension of array J.  LDJ >= MAX(1,NR). */

/*     JTJ     (output) DOUBLE PRECISION array, */
/*                      dimension (LDJTJ,N),    if STOR = 'F', */
/*                      dimension (N*(N+1)/2),  if STOR = 'P'. */
/*             The leading N-by-N (if STOR = 'F'), or N*(N+1)/2 (if */
/*             STOR = 'P') part of this array contains the upper or */
/*             lower triangle of the matrix J'*J + c*I, depending on */
/*             UPLO = 'U', or UPLO = 'L', respectively, stored either as */
/*             a two-dimensional, or one-dimensional array, depending */
/*             on STOR. */

/*     LDJTJ   INTEGER */
/*             The leading dimension of the array JTJ. */
/*             LDJTJ >= MAX(1,N), if STOR = 'F'. */
/*             LDJTJ >= 1,        if STOR = 'P'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             Currently, this array is not used. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= 0. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The matrix product is computed columnn-wise, exploiting the */
/*     symmetry. BLAS 3 routines DGEMM and DSYRK are used if STOR = 'F', */
/*     and BLAS 2 routine DGEMV is used if STOR = 'P'. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2001. */

/*     REVISIONS */

/*     V. Sima, Dec. 2001, Mar. 2002. */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix algebra, matrix operations, */
/*     Wiener system. */

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

#line 174 "NF01BU.f"
    /* Parameter adjustments */
#line 174 "NF01BU.f"
    --ipar;
#line 174 "NF01BU.f"
    --dpar;
#line 174 "NF01BU.f"
    j_dim1 = *ldj;
#line 174 "NF01BU.f"
    j_offset = 1 + j_dim1;
#line 174 "NF01BU.f"
    j -= j_offset;
#line 174 "NF01BU.f"
    --jtj;
#line 174 "NF01BU.f"
    --dwork;
#line 174 "NF01BU.f"

#line 174 "NF01BU.f"
    /* Function Body */
#line 174 "NF01BU.f"
    *info = 0;

#line 176 "NF01BU.f"
    full = lsame_(stor, "F", (ftnlen)1, (ftnlen)1);
#line 177 "NF01BU.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 179 "NF01BU.f"
    if (! (full || lsame_(stor, "P", (ftnlen)1, (ftnlen)1))) {
#line 180 "NF01BU.f"
	*info = -1;
#line 181 "NF01BU.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 182 "NF01BU.f"
	*info = -2;
#line 183 "NF01BU.f"
    } else if (*n < 0) {
#line 184 "NF01BU.f"
	*info = -3;
#line 185 "NF01BU.f"
    } else if (*lipar < 4) {
#line 186 "NF01BU.f"
	*info = -5;
#line 187 "NF01BU.f"
    } else if (*ldpar < 1) {
#line 188 "NF01BU.f"
	*info = -7;
#line 189 "NF01BU.f"
    } else if (*ldjtj < 1 || full && *ldjtj < *n) {
#line 190 "NF01BU.f"
	*info = -11;
#line 191 "NF01BU.f"
    } else if (*ldwork < 0) {
#line 192 "NF01BU.f"
	*info = -13;
#line 193 "NF01BU.f"
    } else {
#line 194 "NF01BU.f"
	st = ipar[1];
#line 195 "NF01BU.f"
	bn = ipar[2];
#line 196 "NF01BU.f"
	bsm = ipar[3];
#line 197 "NF01BU.f"
	bsn = ipar[4];
#line 198 "NF01BU.f"
	nths = bn * bsn;
#line 199 "NF01BU.f"
	if (bn > 1) {
#line 200 "NF01BU.f"
	    m = bn * bsm;
#line 201 "NF01BU.f"
	} else {
#line 202 "NF01BU.f"
	    m = bsm;
#line 203 "NF01BU.f"
	}
/* Computing MIN */
#line 204 "NF01BU.f"
	i__1 = min(st,bn), i__1 = min(i__1,bsm);
#line 204 "NF01BU.f"
	if (min(i__1,bsn) < 0) {
#line 205 "NF01BU.f"
	    *info = -4;
#line 206 "NF01BU.f"
	} else if (*n != nths + st) {
#line 207 "NF01BU.f"
	    *info = -3;
#line 208 "NF01BU.f"
	} else if (*ldj < max(1,m)) {
#line 209 "NF01BU.f"
	    *info = -9;
#line 210 "NF01BU.f"
	}
#line 211 "NF01BU.f"
    }

/*     Return if there are illegal arguments. */

#line 215 "NF01BU.f"
    if (*info != 0) {
#line 216 "NF01BU.f"
	i__1 = -(*info);
#line 216 "NF01BU.f"
	xerbla_("NF01BU", &i__1, (ftnlen)6);
#line 217 "NF01BU.f"
	return 0;
#line 218 "NF01BU.f"
    }

/*     Quick return if possible. */

#line 222 "NF01BU.f"
    if (*n == 0) {
#line 222 "NF01BU.f"
	return 0;
#line 222 "NF01BU.f"
    }

#line 225 "NF01BU.f"
    c__ = dpar[1];

#line 227 "NF01BU.f"
    if (bn <= 1 || bsn == 0 || bsm == 0) {

/*        Special case, l <= 1 or BSN = 0 or BSM = 0: the Jacobian is */
/*        represented as a full matrix. */

#line 232 "NF01BU.f"
	itmp[0] = m;
#line 233 "NF01BU.f"
	nf01bv_(stor, uplo, n, itmp, &c__1, &dpar[1], &c__1, &j[j_offset], 
		ldj, &jtj[1], ldjtj, &dwork[1], ldwork, info, (ftnlen)1, (
		ftnlen)1);
#line 235 "NF01BU.f"
	return 0;
#line 236 "NF01BU.f"
    }

/*     General case: l > 1, BSN > 0, BSM > 0. */

#line 240 "NF01BU.f"
    jl = bsn + 1;

#line 242 "NF01BU.f"
    if (full) {

#line 244 "NF01BU.f"
	nbsn = *n * bsn;

#line 246 "NF01BU.f"
	if (upper) {

/*           Compute the leading upper triangular part (full storage). */

#line 250 "NF01BU.f"
	    dlaset_(uplo, &bsn, &bsn, &c_b9, &c__, &jtj[1], ldjtj, (ftnlen)1);
#line 251 "NF01BU.f"
	    dsyrk_(uplo, "Transpose", &bsn, &bsm, &c_b11, &j[j_offset], ldj, &
		    c_b11, &jtj[1], ldjtj, (ftnlen)1, (ftnlen)9);
#line 253 "NF01BU.f"
	    ibsn = bsn;
#line 254 "NF01BU.f"
	    i1 = nbsn + 1;

#line 256 "NF01BU.f"
	    i__1 = m;
#line 256 "NF01BU.f"
	    i__2 = bsm;
#line 256 "NF01BU.f"
	    for (ibsm = bsm + 1; i__2 < 0 ? ibsm >= i__1 : ibsm <= i__1; ibsm 
		    += i__2) {
#line 257 "NF01BU.f"
		ii = i1 + ibsn;
#line 258 "NF01BU.f"
		dlaset_("Full", &ibsn, &bsn, &c_b9, &c_b9, &jtj[i1], ldjtj, (
			ftnlen)4);
#line 260 "NF01BU.f"
		i1 += nbsn;
#line 261 "NF01BU.f"
		dlaset_(uplo, &bsn, &bsn, &c_b9, &c__, &jtj[ii], ldjtj, (
			ftnlen)1);
#line 262 "NF01BU.f"
		dsyrk_(uplo, "Transpose", &bsn, &bsm, &c_b11, &j[ibsm + 
			j_dim1], ldj, &c_b11, &jtj[ii], ldjtj, (ftnlen)1, (
			ftnlen)9);
#line 264 "NF01BU.f"
		ibsn += bsn;
#line 265 "NF01BU.f"
/* L10: */
#line 265 "NF01BU.f"
	    }

#line 267 "NF01BU.f"
	    if (st > 0) {

/*              Compute the last block column. */

#line 271 "NF01BU.f"
		i__2 = m;
#line 271 "NF01BU.f"
		i__1 = bsm;
#line 271 "NF01BU.f"
		for (ibsm = 1; i__1 < 0 ? ibsm >= i__2 : ibsm <= i__2; ibsm +=
			 i__1) {
#line 272 "NF01BU.f"
		    dgemm_("Transpose", "NoTranspose", &bsn, &st, &bsm, &
			    c_b11, &j[ibsm + j_dim1], ldj, &j[ibsm + jl * 
			    j_dim1], ldj, &c_b9, &jtj[i1], ldjtj, (ftnlen)9, (
			    ftnlen)11);
#line 275 "NF01BU.f"
		    i1 += bsn;
#line 276 "NF01BU.f"
/* L20: */
#line 276 "NF01BU.f"
		}

#line 278 "NF01BU.f"
		dlaset_(uplo, &st, &st, &c_b9, &c__, &jtj[i1], ldjtj, (ftnlen)
			1);
#line 279 "NF01BU.f"
		dsyrk_(uplo, "Transpose", &st, &m, &c_b11, &j[jl * j_dim1 + 1]
			, ldj, &c_b11, &jtj[i1], ldjtj, (ftnlen)1, (ftnlen)9);
#line 281 "NF01BU.f"
	    }

#line 283 "NF01BU.f"
	} else {

/*           Compute the leading lower triangular part (full storage). */

#line 287 "NF01BU.f"
	    ibsn = nths;
#line 288 "NF01BU.f"
	    ii = 1;

#line 290 "NF01BU.f"
	    i__1 = m;
#line 290 "NF01BU.f"
	    i__2 = bsm;
#line 290 "NF01BU.f"
	    for (ibsm = 1; i__2 < 0 ? ibsm >= i__1 : ibsm <= i__1; ibsm += 
		    i__2) {
#line 291 "NF01BU.f"
		i1 = ii + bsn;
#line 292 "NF01BU.f"
		dlaset_(uplo, &bsn, &bsn, &c_b9, &c__, &jtj[ii], ldjtj, (
			ftnlen)1);
#line 293 "NF01BU.f"
		dsyrk_(uplo, "Transpose", &bsn, &bsm, &c_b11, &j[ibsm + 
			j_dim1], ldj, &c_b11, &jtj[ii], ldjtj, (ftnlen)1, (
			ftnlen)9);
#line 295 "NF01BU.f"
		ibsn -= bsn;
#line 296 "NF01BU.f"
		dlaset_("Full", &ibsn, &bsn, &c_b9, &c_b9, &jtj[i1], ldjtj, (
			ftnlen)4);
#line 298 "NF01BU.f"
		ii = i1 + nbsn;
#line 299 "NF01BU.f"
		if (st > 0) {
#line 299 "NF01BU.f"
		    dgemm_("Transpose", "NoTranspose", &st, &bsn, &bsm, &
			    c_b11, &j[ibsm + jl * j_dim1], ldj, &j[ibsm + 
			    j_dim1], ldj, &c_b9, &jtj[i1 + ibsn], ldjtj, (
			    ftnlen)9, (ftnlen)11);
#line 299 "NF01BU.f"
		}
#line 303 "NF01BU.f"
/* L30: */
#line 303 "NF01BU.f"
	    }

#line 305 "NF01BU.f"
	    if (st > 0) {

/*              Compute the last diagonal block. */

#line 309 "NF01BU.f"
		dlaset_(uplo, &st, &st, &c_b9, &c__, &jtj[ii], ldjtj, (ftnlen)
			1);
#line 310 "NF01BU.f"
		dsyrk_(uplo, "Transpose", &st, &m, &c_b11, &j[jl * j_dim1 + 1]
			, ldj, &c_b11, &jtj[ii], ldjtj, (ftnlen)1, (ftnlen)9);
#line 312 "NF01BU.f"
	    }

#line 314 "NF01BU.f"
	}

#line 316 "NF01BU.f"
    } else {

#line 318 "NF01BU.f"
	tmp[0] = 0.;

#line 320 "NF01BU.f"
	if (upper) {

/*           Compute the leading upper triangular part (packed storage). */

#line 324 "NF01BU.f"
	    ibsn = 0;
#line 325 "NF01BU.f"
	    i1 = 1;

#line 327 "NF01BU.f"
	    i__2 = m;
#line 327 "NF01BU.f"
	    i__1 = bsm;
#line 327 "NF01BU.f"
	    for (ibsm = 1; i__1 < 0 ? ibsm >= i__2 : ibsm <= i__2; ibsm += 
		    i__1) {

#line 329 "NF01BU.f"
		i__3 = bsn;
#line 329 "NF01BU.f"
		for (k = 1; k <= i__3; ++k) {
#line 330 "NF01BU.f"
		    ii = i1 + ibsn;
#line 331 "NF01BU.f"
		    dcopy_(&ibsn, tmp, &c__0, &jtj[i1], &c__1);
#line 332 "NF01BU.f"
		    dgemv_("Transpose", &bsm, &k, &c_b11, &j[ibsm + j_dim1], 
			    ldj, &j[ibsm + k * j_dim1], &c__1, &c_b9, &jtj[ii]
			    , &c__1, (ftnlen)9);
#line 334 "NF01BU.f"
		    i1 = ii + k;
#line 335 "NF01BU.f"
		    jtj[i1 - 1] += c__;
#line 336 "NF01BU.f"
/* L40: */
#line 336 "NF01BU.f"
		}

#line 338 "NF01BU.f"
		ibsn += bsn;
#line 339 "NF01BU.f"
/* L50: */
#line 339 "NF01BU.f"
	    }

/*           Compute the last block column. */

#line 343 "NF01BU.f"
	    i__1 = st;
#line 343 "NF01BU.f"
	    for (k = 1; k <= i__1; ++k) {

#line 345 "NF01BU.f"
		i__2 = m;
#line 345 "NF01BU.f"
		i__3 = bsm;
#line 345 "NF01BU.f"
		for (ibsm = 1; i__3 < 0 ? ibsm >= i__2 : ibsm <= i__2; ibsm +=
			 i__3) {
#line 346 "NF01BU.f"
		    dgemv_("Transpose", &bsm, &bsn, &c_b11, &j[ibsm + j_dim1],
			     ldj, &j[ibsm + (bsn + k) * j_dim1], &c__1, &c_b9,
			     &jtj[i1], &c__1, (ftnlen)9);
#line 348 "NF01BU.f"
		    i1 += bsn;
#line 349 "NF01BU.f"
/* L60: */
#line 349 "NF01BU.f"
		}

#line 351 "NF01BU.f"
		dgemv_("Transpose", &m, &k, &c_b11, &j[jl * j_dim1 + 1], ldj, 
			&j[(bsn + k) * j_dim1 + 1], &c__1, &c_b9, &jtj[i1], &
			c__1, (ftnlen)9);
#line 353 "NF01BU.f"
		i1 += k;
#line 354 "NF01BU.f"
		jtj[i1 - 1] += c__;
#line 355 "NF01BU.f"
/* L70: */
#line 355 "NF01BU.f"
	    }

#line 357 "NF01BU.f"
	} else {

/*           Compute the leading lower triangular part (packed storage). */

#line 361 "NF01BU.f"
	    ibsn = nths;
#line 362 "NF01BU.f"
	    ii = 1;

#line 364 "NF01BU.f"
	    i__1 = m;
#line 364 "NF01BU.f"
	    i__3 = bsm;
#line 364 "NF01BU.f"
	    for (ibsm = 1; i__3 < 0 ? ibsm >= i__1 : ibsm <= i__1; ibsm += 
		    i__3) {
#line 365 "NF01BU.f"
		ibsn -= bsn;

#line 367 "NF01BU.f"
		i__2 = bsn;
#line 367 "NF01BU.f"
		for (k = 1; k <= i__2; ++k) {
#line 368 "NF01BU.f"
		    i1 = ii + bsn - k + 1;
#line 369 "NF01BU.f"
		    dcopy_(&ibsn, tmp, &c__0, &jtj[i1], &c__1);
#line 370 "NF01BU.f"
		    i__4 = bsn - k + 1;
#line 370 "NF01BU.f"
		    dgemv_("Transpose", &bsm, &i__4, &c_b11, &j[ibsm + k * 
			    j_dim1], ldj, &j[ibsm + k * j_dim1], &c__1, &c_b9,
			     &jtj[ii], &c__1, (ftnlen)9);
#line 372 "NF01BU.f"
		    jtj[ii] += c__;
#line 373 "NF01BU.f"
		    i1 += ibsn;
#line 374 "NF01BU.f"
		    ii = i1 + st;
#line 375 "NF01BU.f"
		    if (st > 0) {
#line 375 "NF01BU.f"
			dgemv_("Transpose", &bsm, &st, &c_b11, &j[ibsm + jl * 
				j_dim1], ldj, &j[ibsm + k * j_dim1], &c__1, &
				c_b9, &jtj[i1], &c__1, (ftnlen)9);
#line 375 "NF01BU.f"
		    }
#line 378 "NF01BU.f"
/* L80: */
#line 378 "NF01BU.f"
		}

#line 380 "NF01BU.f"
/* L90: */
#line 380 "NF01BU.f"
	    }

/*           Compute the last diagonal block. */

#line 384 "NF01BU.f"
	    i__3 = st;
#line 384 "NF01BU.f"
	    for (k = 1; k <= i__3; ++k) {
#line 385 "NF01BU.f"
		i__1 = st - k + 1;
#line 385 "NF01BU.f"
		dgemv_("Transpose", &m, &i__1, &c_b11, &j[(bsn + k) * j_dim1 
			+ 1], ldj, &j[(bsn + k) * j_dim1 + 1], &c__1, &c_b9, &
			jtj[ii], &c__1, (ftnlen)9);
#line 387 "NF01BU.f"
		jtj[ii] += c__;
#line 388 "NF01BU.f"
		ii = ii + st - k + 1;
#line 389 "NF01BU.f"
/* L100: */
#line 389 "NF01BU.f"
	    }

#line 391 "NF01BU.f"
	}

#line 393 "NF01BU.f"
    }

#line 395 "NF01BU.f"
    return 0;

/* *** Last line of NF01BU *** */
} /* nf01bu_ */

