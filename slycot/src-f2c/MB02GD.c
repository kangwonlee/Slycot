#line 1 "MB02GD.f"
/* MB02GD.f -- translated by f2c (version 20100827).
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

#line 1 "MB02GD.f"
/* Table of constant values */

static doublereal c_b12 = 1.;
static integer c__1 = 1;
static doublereal c_b24 = 0.;
static integer c_n1 = -1;
static integer c__2 = 2;

/* Subroutine */ int mb02gd_(char *typet, char *triu, integer *k, integer *n, 
	integer *nl, integer *p, integer *s, doublereal *t, integer *ldt, 
	doublereal *rb, integer *ldrb, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen typet_len, ftnlen triu_len)
{
    /* System generated locals */
    integer rb_dim1, rb_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4, 
	    i__5;

    /* Local variables */
    static integer i__, j, nb, jj, kk, len;
    static doublereal dum[1];
    static integer pre, pdw, rnk, len2, head, lenr, ierr;
    static logical ltri;
    static integer ipvt[1], posr, sizr, stps;
    extern /* Subroutine */ int mb02cu_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen), 
	    mb02cv_(char *, char *, integer *, integer *, integer *, integer *
	    , integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical isrow;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dpotrf_(
	    char *, integer *, doublereal *, integer *, integer *, ftnlen);
    static integer wrkmin;
    static char struct__[1];
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

/*     To compute the Cholesky factor of a banded symmetric positive */
/*     definite (s.p.d.) block Toeplitz matrix, defined by either its */
/*     first block row, or its first block column, depending on the */
/*     routine parameter TYPET. */

/*     By subsequent calls of this routine the Cholesky factor can be */
/*     computed block column by block column. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYPET   CHARACTER*1 */
/*             Specifies the type of T, as follows: */
/*             = 'R':  T contains the first block row of an s.p.d. block */
/*                     Toeplitz matrix; the Cholesky factor is upper */
/*                     triangular; */
/*             = 'C':  T contains the first block column of an s.p.d. */
/*                     block Toeplitz matrix; the Cholesky factor is */
/*                     lower triangular. This choice results in a column */
/*                     oriented algorithm which is usually faster. */
/*             Note:   in the sequel, the notation x / y means that */
/*                     x corresponds to TYPET = 'R' and y corresponds to */
/*                     TYPET = 'C'. */

/*     TRIU    CHARACTER*1 */
/*             Specifies the structure of the last block in T, as */
/*             follows: */
/*             = 'N':  the last block has no special structure; */
/*             = 'T':  the last block is lower / upper triangular. */

/*     Input/Output Parameters */

/*     K       (input)  INTEGER */
/*             The number of rows / columns in T, which should be equal */
/*             to the blocksize.  K >= 0. */

/*     N       (input)  INTEGER */
/*             The number of blocks in T.  N >= 1. */
/*             If TRIU = 'N',   N >= 1; */
/*             if TRIU = 'T',   N >= 2. */

/*     NL      (input)  INTEGER */
/*             The lower block bandwidth, i.e., NL + 1 is the number of */
/*             nonzero blocks in the first block column of the block */
/*             Toeplitz matrix. */
/*             If TRIU = 'N',   0 <= NL < N; */
/*             if TRIU = 'T',   1 <= NL < N. */

/*     P       (input)  INTEGER */
/*             The number of previously computed block rows / columns of */
/*             the Cholesky factor.  0 <= P <= N. */

/*     S       (input)  INTEGER */
/*             The number of block rows / columns of the Cholesky factor */
/*             to compute.  0 <= S <= N - P. */

/*     T       (input/output)  DOUBLE PRECISION array, dimension */
/*             (LDT,(NL+1)*K) / (LDT,K) */
/*             On entry, if P = 0, the leading K-by-(NL+1)*K / */
/*             (NL+1)*K-by-K part of this array must contain the first */
/*             block row / column of an s.p.d. block Toeplitz matrix. */
/*             On entry, if P > 0, the leading K-by-(NL+1)*K / */
/*             (NL+1)*K-by-K part of this array must contain the P-th */
/*             block row / column of the Cholesky factor. */
/*             On exit, if INFO = 0, then the leading K-by-(NL+1)*K / */
/*             (NL+1)*K-by-K part of this array contains the (P+S)-th */
/*             block row / column of the Cholesky factor. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T. */
/*             LDT >= MAX(1,K) / MAX(1,(NL+1)*K). */

/*     RB      (input/output)  DOUBLE PRECISION array, dimension */
/*             (LDRB,MIN(P+NL+S,N)*K) / (LDRB,MIN(P+S,N)*K) */
/*             On entry, if TYPET = 'R'  and  TRIU = 'N'  and  P > 0, */
/*             the leading (NL+1)*K-by-MIN(NL,N-P)*K part of this array */
/*             must contain the (P*K+1)-st to ((P+NL)*K)-th columns */
/*             of the upper Cholesky factor in banded format from a */
/*             previous call of this routine. */
/*             On entry, if TYPET = 'R'  and  TRIU = 'T'  and  P > 0, */
/*             the leading (NL*K+1)-by-MIN(NL,N-P)*K part of this array */
/*             must contain the (P*K+1)-st to (MIN(P+NL,N)*K)-th columns */
/*             of the upper Cholesky factor in banded format from a */
/*             previous call of this routine. */
/*             On exit, if TYPET = 'R'  and  TRIU = 'N', the leading */
/*             (NL+1)*K-by-MIN(NL+S,N-P)*K part of this array contains */
/*             the (P*K+1)-st to (MIN(P+NL+S,N)*K)-th columns of the */
/*             upper Cholesky factor in banded format. */
/*             On exit, if TYPET = 'R'  and  TRIU = 'T', the leading */
/*             (NL*K+1)-by-MIN(NL+S,N-P)*K part of this array contains */
/*             the (P*K+1)-st to (MIN(P+NL+S,N)*K)-th columns of the */
/*             upper Cholesky factor in banded format. */
/*             On exit, if TYPET = 'C'  and  TRIU = 'N', the leading */
/*             (NL+1)*K-by-MIN(S,N-P)*K part of this array contains */
/*             the (P*K+1)-st to (MIN(P+S,N)*K)-th columns of the lower */
/*             Cholesky factor in banded format. */
/*             On exit, if TYPET = 'C'  and  TRIU = 'T', the leading */
/*             (NL*K+1)-by-MIN(S,N-P)*K part of this array contains */
/*             the (P*K+1)-st to (MIN(P+S,N)*K)-th columns of the lower */
/*             Cholesky factor in banded format. */
/*             For further details regarding the band storage scheme see */
/*             the documentation of the LAPACK routine DPBTF2. */

/*     LDRB    INTEGER */
/*             The leading dimension of the array RB. */
/*             If TRIU = 'N',   LDRB >= MAX( (NL+1)*K,1 ); */
/*             if TRIU = 'T',   LDRB >= NL*K+1. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -13,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */
/*             The first 1 + ( NL + 1 )*K*K elements of DWORK should be */
/*             preserved during successive calls of the routine. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 1 + ( NL + 1 )*K*K + NL*K. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the reduction algorithm failed. The Toeplitz matrix */
/*                   associated with T is not (numerically) positive */
/*                   definite. */

/*     METHOD */

/*     Householder transformations and modified hyperbolic rotations */
/*     are used in the Schur algorithm [1], [2]. */

/*     REFERENCES */

/*     [1] Kailath, T. and Sayed, A. */
/*         Fast Reliable Algorithms for Matrices with Structure. */
/*         SIAM Publications, Philadelphia, 1999. */

/*     [2] Kressner, D. and Van Dooren, P. */
/*         Factorizations and linear system solvers for matrices with */
/*         Toeplitz structure. */
/*         SLICOT Working Note 2000-2, 2000. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is numerically stable. */
/*                                3 */
/*     The algorithm requires O( K *N*NL ) floating point operations. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, May 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, June 2001, */
/*     Mar. 2004. */

/*     KEYWORDS */

/*     Elementary matrix operations, Householder transformation, matrix */
/*     operations, Toeplitz matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

#line 227 "MB02GD.f"
    /* Parameter adjustments */
#line 227 "MB02GD.f"
    t_dim1 = *ldt;
#line 227 "MB02GD.f"
    t_offset = 1 + t_dim1;
#line 227 "MB02GD.f"
    t -= t_offset;
#line 227 "MB02GD.f"
    rb_dim1 = *ldrb;
#line 227 "MB02GD.f"
    rb_offset = 1 + rb_dim1;
#line 227 "MB02GD.f"
    rb -= rb_offset;
#line 227 "MB02GD.f"
    --dwork;
#line 227 "MB02GD.f"

#line 227 "MB02GD.f"
    /* Function Body */
#line 227 "MB02GD.f"
    *info = 0;
#line 228 "MB02GD.f"
    ltri = lsame_(triu, "T", (ftnlen)1, (ftnlen)1);
#line 229 "MB02GD.f"
    lenr = (*nl + 1) * *k;
#line 230 "MB02GD.f"
    if (ltri) {
#line 231 "MB02GD.f"
	sizr = *nl * *k + 1;
#line 232 "MB02GD.f"
    } else {
#line 233 "MB02GD.f"
	sizr = lenr;
#line 234 "MB02GD.f"
    }
#line 235 "MB02GD.f"
    isrow = lsame_(typet, "R", (ftnlen)1, (ftnlen)1);
#line 236 "MB02GD.f"
    wrkmin = (lenr + *nl) * *k + 1;

/*     Check the scalar input parameters. */

#line 240 "MB02GD.f"
    if (! (isrow || lsame_(typet, "C", (ftnlen)1, (ftnlen)1))) {
#line 241 "MB02GD.f"
	*info = -1;
#line 242 "MB02GD.f"
    } else if (! (ltri || lsame_(triu, "N", (ftnlen)1, (ftnlen)1))) {
#line 243 "MB02GD.f"
	*info = -2;
#line 244 "MB02GD.f"
    } else if (*k < 0) {
#line 245 "MB02GD.f"
	*info = -3;
#line 246 "MB02GD.f"
    } else if (ltri && *n < 2 || ! ltri && *n < 1) {
#line 248 "MB02GD.f"
	*info = -4;
#line 249 "MB02GD.f"
    } else if (*nl >= *n || ltri && *nl < 1 || ! ltri && *nl < 0) {
#line 251 "MB02GD.f"
	*info = -5;
#line 252 "MB02GD.f"
    } else if (*p < 0 || *p > *n) {
#line 253 "MB02GD.f"
	*info = -6;
#line 254 "MB02GD.f"
    } else if (*s < 0 || *s > *n - *p) {
#line 255 "MB02GD.f"
	*info = -7;
#line 256 "MB02GD.f"
    } else if (isrow && *ldt < max(1,*k) || ! isrow && *ldt < max(1,lenr)) {
#line 259 "MB02GD.f"
	*info = -9;
#line 260 "MB02GD.f"
    } else if (ltri && *ldrb < sizr || ! ltri && *ldrb < max(1,lenr)) {
#line 263 "MB02GD.f"
	*info = -11;
#line 264 "MB02GD.f"
    } else if (*ldwork < wrkmin) {
#line 265 "MB02GD.f"
	dwork[1] = (doublereal) wrkmin;
#line 266 "MB02GD.f"
	*info = -13;
#line 267 "MB02GD.f"
    }

/*     Return if there were illegal values. */

#line 271 "MB02GD.f"
    if (*info != 0) {
#line 272 "MB02GD.f"
	i__1 = -(*info);
#line 272 "MB02GD.f"
	xerbla_("MB02GD", &i__1, (ftnlen)6);
#line 273 "MB02GD.f"
	return 0;
#line 274 "MB02GD.f"
    }

/*     Quick return if possible. */

#line 278 "MB02GD.f"
    if (*s * *k == 0) {
#line 279 "MB02GD.f"
	dwork[1] = 1.;
#line 280 "MB02GD.f"
	return 0;
#line 281 "MB02GD.f"
    }

/*     Compute the generator if P = 0. */

#line 285 "MB02GD.f"
    if (*p == 0) {
#line 286 "MB02GD.f"
	if (isrow) {
#line 287 "MB02GD.f"
	    dpotrf_("Upper", k, &t[t_offset], ldt, &ierr, (ftnlen)5);
#line 288 "MB02GD.f"
	    if (ierr != 0) {

/*              Error return:  The matrix is not positive definite. */

#line 292 "MB02GD.f"
		*info = 1;
#line 293 "MB02GD.f"
		return 0;
#line 294 "MB02GD.f"
	    }
#line 295 "MB02GD.f"
	    if (*nl > 0) {
#line 295 "MB02GD.f"
		i__1 = *nl * *k;
#line 295 "MB02GD.f"
		dtrsm_("Left", "Upper", "Transpose", "NonUnit", k, &i__1, &
			c_b12, &t[t_offset], ldt, &t[(*k + 1) * t_dim1 + 1], 
			ldt, (ftnlen)4, (ftnlen)5, (ftnlen)9, (ftnlen)7);
#line 295 "MB02GD.f"
	    }

/*           Copy the first block row to RB. */

#line 301 "MB02GD.f"
	    if (ltri) {

#line 303 "MB02GD.f"
		i__1 = lenr - *k;
#line 303 "MB02GD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 304 "MB02GD.f"
		    i__2 = min(i__,*k);
/* Computing MAX */
#line 304 "MB02GD.f"
		    i__3 = sizr - i__ + 1;
#line 304 "MB02GD.f"
		    dcopy_(&i__2, &t[i__ * t_dim1 + 1], &c__1, &rb[max(i__3,1)
			     + i__ * rb_dim1], &c__1);
#line 306 "MB02GD.f"
/* L10: */
#line 306 "MB02GD.f"
		}

#line 308 "MB02GD.f"
		for (i__ = *k; i__ >= 1; --i__) {
#line 309 "MB02GD.f"
		    dcopy_(&i__, &t[*k - i__ + 1 + (lenr - i__ + 1) * t_dim1],
			     &c__1, &rb[(lenr - i__ + 1) * rb_dim1 + 1], &
			    c__1);
#line 311 "MB02GD.f"
/* L20: */
#line 311 "MB02GD.f"
		}

#line 313 "MB02GD.f"
	    } else {

#line 315 "MB02GD.f"
		i__1 = lenr;
#line 315 "MB02GD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 316 "MB02GD.f"
		    i__2 = min(i__,*k);
/* Computing MAX */
#line 316 "MB02GD.f"
		    i__3 = sizr - i__ + 1;
#line 316 "MB02GD.f"
		    dcopy_(&i__2, &t[i__ * t_dim1 + 1], &c__1, &rb[max(i__3,1)
			     + i__ * rb_dim1], &c__1);
#line 318 "MB02GD.f"
/* L30: */
#line 318 "MB02GD.f"
		}

#line 320 "MB02GD.f"
	    }

/*           Quick return if N = 1. */

#line 324 "MB02GD.f"
	    if (*n == 1) {
#line 325 "MB02GD.f"
		dwork[1] = 1.;
#line 326 "MB02GD.f"
		return 0;
#line 327 "MB02GD.f"
	    }

#line 329 "MB02GD.f"
	    i__1 = *nl * *k;
#line 329 "MB02GD.f"
	    dlacpy_("All", k, &i__1, &t[(*k + 1) * t_dim1 + 1], ldt, &dwork[2]
		    , k, (ftnlen)3);
#line 330 "MB02GD.f"
	    dlaset_("All", k, k, &c_b24, &c_b24, &dwork[*nl * *k * *k + 2], k,
		     (ftnlen)3);
#line 331 "MB02GD.f"
	    posr = *k + 1;
#line 332 "MB02GD.f"
	} else {
#line 333 "MB02GD.f"
	    dpotrf_("Lower", k, &t[t_offset], ldt, &ierr, (ftnlen)5);
#line 334 "MB02GD.f"
	    if (ierr != 0) {

/*              Error return:  The matrix is not positive definite. */

#line 338 "MB02GD.f"
		*info = 1;
#line 339 "MB02GD.f"
		return 0;
#line 340 "MB02GD.f"
	    }
#line 341 "MB02GD.f"
	    if (*nl > 0) {
#line 341 "MB02GD.f"
		i__1 = *nl * *k;
#line 341 "MB02GD.f"
		dtrsm_("Right", "Lower", "Transpose", "NonUnit", &i__1, k, &
			c_b12, &t[t_offset], ldt, &t[*k + 1 + t_dim1], ldt, (
			ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)7);
#line 341 "MB02GD.f"
	    }

/*           Copy the first block column to RB. */

#line 347 "MB02GD.f"
	    posr = 1;
#line 348 "MB02GD.f"
	    if (ltri) {

#line 350 "MB02GD.f"
		i__1 = *k;
#line 350 "MB02GD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 351 "MB02GD.f"
		    dcopy_(&sizr, &t[i__ + i__ * t_dim1], &c__1, &rb[posr * 
			    rb_dim1 + 1], &c__1);
#line 352 "MB02GD.f"
		    ++posr;
#line 353 "MB02GD.f"
/* L40: */
#line 353 "MB02GD.f"
		}

#line 355 "MB02GD.f"
	    } else {

#line 357 "MB02GD.f"
		i__1 = *k;
#line 357 "MB02GD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 358 "MB02GD.f"
		    i__2 = lenr - i__ + 1;
#line 358 "MB02GD.f"
		    dcopy_(&i__2, &t[i__ + i__ * t_dim1], &c__1, &rb[posr * 
			    rb_dim1 + 1], &c__1);
#line 359 "MB02GD.f"
		    if (lenr < *n * *k && i__ > 1) {
#line 360 "MB02GD.f"
			i__2 = i__ - 1;
#line 360 "MB02GD.f"
			dlaset_("All", &i__2, &c__1, &c_b24, &c_b24, &rb[lenr 
				- i__ + 2 + posr * rb_dim1], ldrb, (ftnlen)3);
#line 362 "MB02GD.f"
		    }
#line 363 "MB02GD.f"
		    ++posr;
#line 364 "MB02GD.f"
/* L50: */
#line 364 "MB02GD.f"
		}

#line 366 "MB02GD.f"
	    }

/*           Quick return if N = 1. */

#line 370 "MB02GD.f"
	    if (*n == 1) {
#line 371 "MB02GD.f"
		dwork[1] = 1.;
#line 372 "MB02GD.f"
		return 0;
#line 373 "MB02GD.f"
	    }

#line 375 "MB02GD.f"
	    i__1 = *nl * *k;
#line 375 "MB02GD.f"
	    dlacpy_("All", &i__1, k, &t[*k + 1 + t_dim1], ldt, &dwork[2], &
		    lenr, (ftnlen)3);
#line 376 "MB02GD.f"
	    dlaset_("All", k, k, &c_b24, &c_b24, &dwork[*nl * *k + 2], &lenr, 
		    (ftnlen)3);
#line 377 "MB02GD.f"
	}
#line 378 "MB02GD.f"
	pre = 1;
#line 379 "MB02GD.f"
	stps = *s - 1;
#line 380 "MB02GD.f"
    } else {
#line 381 "MB02GD.f"
	pre = *p;
#line 382 "MB02GD.f"
	stps = *s;
#line 383 "MB02GD.f"
	posr = 1;
#line 384 "MB02GD.f"
    }

#line 386 "MB02GD.f"
    pdw = lenr * *k + 1;
#line 387 "MB02GD.f"
    head = (pre - 1) * *k % lenr;

/*     Determine block size for the involved block Householder */
/*     transformations. */

#line 392 "MB02GD.f"
    if (isrow) {
/* Computing MIN */
#line 393 "MB02GD.f"
	i__1 = ilaenv_(&c__1, "DGEQRF", " ", k, &lenr, &c_n1, &c_n1, (ftnlen)
		6, (ftnlen)1);
#line 393 "MB02GD.f"
	nb = min(i__1,*k);
#line 394 "MB02GD.f"
    } else {
/* Computing MIN */
#line 395 "MB02GD.f"
	i__1 = ilaenv_(&c__1, "DGELQF", " ", &lenr, k, &c_n1, &c_n1, (ftnlen)
		6, (ftnlen)1);
#line 395 "MB02GD.f"
	nb = min(i__1,*k);
#line 396 "MB02GD.f"
    }
#line 397 "MB02GD.f"
    kk = pdw + (*k << 2);
#line 398 "MB02GD.f"
    wrkopt = kk + lenr * nb;
#line 399 "MB02GD.f"
    kk = *ldwork - kk;
#line 400 "MB02GD.f"
    if (kk < lenr * nb) {
#line 400 "MB02GD.f"
	nb = kk / lenr;
#line 400 "MB02GD.f"
    }
#line 401 "MB02GD.f"
    if (isrow) {
/* Computing MAX */
#line 402 "MB02GD.f"
	i__1 = 2, i__2 = ilaenv_(&c__2, "DGEQRF", " ", k, &lenr, &c_n1, &c_n1,
		 (ftnlen)6, (ftnlen)1);
#line 402 "MB02GD.f"
	nbmin = max(i__1,i__2);
#line 403 "MB02GD.f"
    } else {
/* Computing MAX */
#line 404 "MB02GD.f"
	i__1 = 2, i__2 = ilaenv_(&c__2, "DGELQF", " ", &lenr, k, &c_n1, &c_n1,
		 (ftnlen)6, (ftnlen)1);
#line 404 "MB02GD.f"
	nbmin = max(i__1,i__2);
#line 405 "MB02GD.f"
    }
#line 406 "MB02GD.f"
    if (nb < nbmin) {
#line 406 "MB02GD.f"
	nb = 0;
#line 406 "MB02GD.f"
    }

/*     Generator reduction process. */

#line 410 "MB02GD.f"
    if (isrow) {

#line 412 "MB02GD.f"
	i__1 = pre + stps - 1;
#line 412 "MB02GD.f"
	for (i__ = pre; i__ <= i__1; ++i__) {
#line 413 "MB02GD.f"
	    i__2 = *ldwork - pdw - (*k << 2);
#line 413 "MB02GD.f"
	    mb02cu_("Row", k, k, k, &nb, &t[t_offset], ldt, dum, &c__1, &
		    dwork[head * *k + 2], k, &rnk, ipvt, &dwork[pdw + 1], &
		    c_b24, &dwork[pdw + (*k << 2) + 1], &i__2, &ierr, (ftnlen)
		    3);

#line 417 "MB02GD.f"
	    if (ierr != 0) {

/*              Error return:  The positive definiteness is (numerically) */
/*                             not satisfied. */

#line 422 "MB02GD.f"
		*info = 1;
#line 423 "MB02GD.f"
		return 0;
#line 424 "MB02GD.f"
	    }

/* Computing MAX */
/* Computing MIN */
#line 426 "MB02GD.f"
	    i__3 = (*n - i__) * *k - *k, i__4 = lenr - head - *k;
#line 426 "MB02GD.f"
	    i__2 = min(i__3,i__4);
#line 426 "MB02GD.f"
	    len = max(i__2,0);
/* Computing MAX */
/* Computing MIN */
#line 427 "MB02GD.f"
	    i__3 = (*n - i__) * *k - len - *k;
#line 427 "MB02GD.f"
	    i__2 = min(i__3,head);
#line 427 "MB02GD.f"
	    len2 = max(i__2,0);
#line 428 "MB02GD.f"
	    if (len == lenr - *k) {
#line 429 "MB02GD.f"
		*(unsigned char *)struct__ = *(unsigned char *)triu;
#line 430 "MB02GD.f"
	    } else {
#line 431 "MB02GD.f"
		*(unsigned char *)struct__ = 'N';
#line 432 "MB02GD.f"
	    }
#line 433 "MB02GD.f"
	    i__2 = *ldwork - pdw - (*k << 2);
#line 433 "MB02GD.f"
	    mb02cv_("Row", struct__, k, &len, k, k, &nb, &c_n1, dum, &c__1, 
		    dum, &c__1, &dwork[head * *k + 2], k, &t[(*k + 1) * 
		    t_dim1 + 1], ldt, dum, &c__1, &dwork[(head + *k) * *k + 2]
		    , k, &dwork[pdw + 1], &dwork[pdw + (*k << 2) + 1], &i__2, 
		    &ierr, (ftnlen)3, (ftnlen)1);

#line 438 "MB02GD.f"
	    if ((*n - i__) * *k >= lenr) {
#line 439 "MB02GD.f"
		*(unsigned char *)struct__ = *(unsigned char *)triu;
#line 440 "MB02GD.f"
	    } else {
#line 441 "MB02GD.f"
		*(unsigned char *)struct__ = 'N';
#line 442 "MB02GD.f"
	    }
#line 443 "MB02GD.f"
	    i__2 = *ldwork - pdw - (*k << 2);
#line 443 "MB02GD.f"
	    mb02cv_("Row", struct__, k, &len2, k, k, &nb, &c_n1, dum, &c__1, 
		    dum, &c__1, &dwork[head * *k + 2], k, &t[(*k + len + 1) * 
		    t_dim1 + 1], ldt, dum, &c__1, &dwork[2], k, &dwork[pdw + 
		    1], &dwork[pdw + (*k << 2) + 1], &i__2, &ierr, (ftnlen)3, 
		    (ftnlen)1);

#line 448 "MB02GD.f"
	    dlaset_("All", k, k, &c_b24, &c_b24, &dwork[head * *k + 2], k, (
		    ftnlen)3);

/*           Copy current block row to RB. */

#line 452 "MB02GD.f"
	    if (ltri) {

/* Computing MIN */
#line 454 "MB02GD.f"
		i__3 = len + len2 + *k, i__4 = lenr - *k;
#line 454 "MB02GD.f"
		i__2 = min(i__3,i__4);
#line 454 "MB02GD.f"
		for (j = 1; j <= i__2; ++j) {
#line 455 "MB02GD.f"
		    i__3 = min(j,*k);
/* Computing MAX */
#line 455 "MB02GD.f"
		    i__4 = sizr - j + 1;
#line 455 "MB02GD.f"
		    dcopy_(&i__3, &t[j * t_dim1 + 1], &c__1, &rb[max(i__4,1) 
			    + (posr + j - 1) * rb_dim1], &c__1);
#line 457 "MB02GD.f"
/* L60: */
#line 457 "MB02GD.f"
		}

#line 459 "MB02GD.f"
		if (len + len2 + *k >= lenr) {

#line 461 "MB02GD.f"
		    for (jj = *k; jj >= 1; --jj) {
#line 462 "MB02GD.f"
			dcopy_(&jj, &t[*k - jj + 1 + (lenr - jj + 1) * t_dim1]
				, &c__1, &rb[(posr + lenr - jj) * rb_dim1 + 1]
				, &c__1);
#line 464 "MB02GD.f"
/* L70: */
#line 464 "MB02GD.f"
		    }

#line 466 "MB02GD.f"
		}
#line 467 "MB02GD.f"
		posr += *k;

#line 469 "MB02GD.f"
	    } else {

#line 471 "MB02GD.f"
		i__2 = len + len2 + *k;
#line 471 "MB02GD.f"
		for (j = 1; j <= i__2; ++j) {
#line 472 "MB02GD.f"
		    i__3 = min(j,*k);
/* Computing MAX */
#line 472 "MB02GD.f"
		    i__4 = sizr - j + 1;
#line 472 "MB02GD.f"
		    dcopy_(&i__3, &t[j * t_dim1 + 1], &c__1, &rb[max(i__4,1) 
			    + (posr + j - 1) * rb_dim1], &c__1);
#line 474 "MB02GD.f"
		    if (j > lenr - *k) {
#line 475 "MB02GD.f"
			i__3 = sizr - j;
#line 475 "MB02GD.f"
			dlaset_("All", &i__3, &c__1, &c_b24, &c_b24, &rb[(
				posr + j - 1) * rb_dim1 + 1], &c__1, (ftnlen)
				3);
#line 477 "MB02GD.f"
		    }
#line 478 "MB02GD.f"
/* L80: */
#line 478 "MB02GD.f"
		}

#line 480 "MB02GD.f"
		posr += *k;
#line 481 "MB02GD.f"
	    }
#line 482 "MB02GD.f"
	    head = (head + *k) % lenr;
#line 483 "MB02GD.f"
/* L90: */
#line 483 "MB02GD.f"
	}

#line 485 "MB02GD.f"
    } else {

#line 487 "MB02GD.f"
	i__1 = pre + stps - 1;
#line 487 "MB02GD.f"
	for (i__ = pre; i__ <= i__1; ++i__) {

#line 489 "MB02GD.f"
	    i__2 = *ldwork - pdw - (*k << 2);
#line 489 "MB02GD.f"
	    mb02cu_("Column", k, k, k, &nb, &t[t_offset], ldt, dum, &c__1, &
		    dwork[head + 2], &lenr, &rnk, ipvt, &dwork[pdw + 1], &
		    c_b24, &dwork[pdw + (*k << 2) + 1], &i__2, &ierr, (ftnlen)
		    6);

#line 493 "MB02GD.f"
	    if (ierr != 0) {

/*              Error return:  The positive definiteness is (numerically) */
/*                             not satisfied. */

#line 498 "MB02GD.f"
		*info = 1;
#line 499 "MB02GD.f"
		return 0;
#line 500 "MB02GD.f"
	    }

/* Computing MAX */
/* Computing MIN */
#line 502 "MB02GD.f"
	    i__3 = (*n - i__) * *k - *k, i__4 = lenr - head - *k;
#line 502 "MB02GD.f"
	    i__2 = min(i__3,i__4);
#line 502 "MB02GD.f"
	    len = max(i__2,0);
/* Computing MAX */
/* Computing MIN */
#line 503 "MB02GD.f"
	    i__3 = (*n - i__) * *k - len - *k;
#line 503 "MB02GD.f"
	    i__2 = min(i__3,head);
#line 503 "MB02GD.f"
	    len2 = max(i__2,0);
#line 504 "MB02GD.f"
	    if (len == lenr - *k) {
#line 505 "MB02GD.f"
		*(unsigned char *)struct__ = *(unsigned char *)triu;
#line 506 "MB02GD.f"
	    } else {
#line 507 "MB02GD.f"
		*(unsigned char *)struct__ = 'N';
#line 508 "MB02GD.f"
	    }
#line 509 "MB02GD.f"
	    i__2 = *ldwork - pdw - (*k << 2);
#line 509 "MB02GD.f"
	    mb02cv_("Column", struct__, k, &len, k, k, &nb, &c_n1, dum, &c__1,
		     dum, &c__1, &dwork[head + 2], &lenr, &t[*k + 1 + t_dim1],
		     ldt, dum, &c__1, &dwork[head + *k + 2], &lenr, &dwork[
		    pdw + 1], &dwork[pdw + (*k << 2) + 1], &i__2, &ierr, (
		    ftnlen)6, (ftnlen)1);

#line 514 "MB02GD.f"
	    if ((*n - i__) * *k >= lenr) {
#line 515 "MB02GD.f"
		*(unsigned char *)struct__ = *(unsigned char *)triu;
#line 516 "MB02GD.f"
	    } else {
#line 517 "MB02GD.f"
		*(unsigned char *)struct__ = 'N';
#line 518 "MB02GD.f"
	    }
#line 519 "MB02GD.f"
	    i__2 = *ldwork - pdw - (*k << 2);
#line 519 "MB02GD.f"
	    mb02cv_("Column", struct__, k, &len2, k, k, &nb, &c_n1, dum, &
		    c__1, dum, &c__1, &dwork[head + 2], &lenr, &t[*k + len + 
		    1 + t_dim1], ldt, dum, &c__1, &dwork[2], &lenr, &dwork[
		    pdw + 1], &dwork[pdw + (*k << 2) + 1], &i__2, &ierr, (
		    ftnlen)6, (ftnlen)1);

#line 524 "MB02GD.f"
	    dlaset_("All", k, k, &c_b24, &c_b24, &dwork[head + 2], &lenr, (
		    ftnlen)3);

/*           Copy current block column to RB. */

#line 528 "MB02GD.f"
	    if (ltri) {

#line 530 "MB02GD.f"
		i__2 = *k;
#line 530 "MB02GD.f"
		for (j = 1; j <= i__2; ++j) {
/* Computing MIN */
#line 531 "MB02GD.f"
		    i__4 = sizr, i__5 = (*n - i__) * *k - j + 1;
#line 531 "MB02GD.f"
		    i__3 = min(i__4,i__5);
#line 531 "MB02GD.f"
		    dcopy_(&i__3, &t[j + j * t_dim1], &c__1, &rb[posr * 
			    rb_dim1 + 1], &c__1);
#line 533 "MB02GD.f"
		    ++posr;
#line 534 "MB02GD.f"
/* L100: */
#line 534 "MB02GD.f"
		}

#line 536 "MB02GD.f"
	    } else {

#line 538 "MB02GD.f"
		i__2 = *k;
#line 538 "MB02GD.f"
		for (j = 1; j <= i__2; ++j) {
/* Computing MIN */
#line 539 "MB02GD.f"
		    i__4 = sizr - j + 1, i__5 = (*n - i__) * *k - j + 1;
#line 539 "MB02GD.f"
		    i__3 = min(i__4,i__5);
#line 539 "MB02GD.f"
		    dcopy_(&i__3, &t[j + j * t_dim1], &c__1, &rb[posr * 
			    rb_dim1 + 1], &c__1);
#line 541 "MB02GD.f"
		    if (lenr < (*n - i__) * *k) {
#line 542 "MB02GD.f"
			i__3 = j - 1;
/* Computing MIN */
#line 542 "MB02GD.f"
			i__4 = sizr - j + 1, i__5 = (*n - i__) * *k - j + 1;
#line 542 "MB02GD.f"
			dlaset_("All", &i__3, &c__1, &c_b24, &c_b24, &rb[min(
				i__4,i__5) + 1 + posr * rb_dim1], ldrb, (
				ftnlen)3);
#line 545 "MB02GD.f"
		    }
#line 546 "MB02GD.f"
		    ++posr;
#line 547 "MB02GD.f"
/* L110: */
#line 547 "MB02GD.f"
		}

#line 549 "MB02GD.f"
	    }
#line 550 "MB02GD.f"
	    head = (head + *k) % lenr;
#line 551 "MB02GD.f"
/* L120: */
#line 551 "MB02GD.f"
	}

#line 553 "MB02GD.f"
    }
#line 554 "MB02GD.f"
    dwork[1] = (doublereal) wrkopt;
#line 555 "MB02GD.f"
    return 0;

/* *** Last line of MB02GD *** */
} /* mb02gd_ */

