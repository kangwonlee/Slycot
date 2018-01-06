#line 1 "MB02HD.f"
/* MB02HD.f -- translated by f2c (version 20100827).
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

#line 1 "MB02HD.f"
/* Table of constant values */

static doublereal c_b10 = 0.;
static doublereal c_b15 = 1.;
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* Subroutine */ int mb02hd_(char *triu, integer *k, integer *l, integer *m, 
	integer *ml, integer *n, integer *nu, integer *p, integer *s, 
	doublereal *tc, integer *ldtc, doublereal *tr, integer *ldtr, 
	doublereal *rb, integer *ldrb, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen triu_len)
{
    /* System generated locals */
    integer rb_dim1, rb_offset, tc_dim1, tc_offset, tr_dim1, tr_offset, i__1, 
	    i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j, x, nb, kk, pt, pdc, len, pdr, pre, pfr, pdw, rnk, 
	    pnr, col2, len2, head, lenc, lenl, lenr, ierr;
    static logical ltri;
    static integer ipvt[1], posr, sizr, stps;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    mb02cu_(char *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen), dgemm_(char *, char *
	    , integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen), mb02cv_(char *, char *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dgeqrf_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), dlaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dorgqr_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);
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

/*     To compute, for a banded K*M-by-L*N block Toeplitz matrix T with */
/*     block size (K,L), specified by the nonzero blocks of its first */
/*     block column TC and row TR, a LOWER triangular matrix R (in band */
/*     storage scheme) such that */
/*                          T          T */
/*                         T  T  =  R R .                             (1) */

/*     It is assumed that the first MIN(M*K, N*L) columns of T are */
/*     linearly independent. */

/*     By subsequent calls of this routine, the matrix R can be computed */
/*     block column by block column. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRIU    CHARACTER*1 */
/*             Specifies the structure, if any, of the last blocks in TC */
/*             and TR, as follows: */
/*             = 'N':  TC and TR have no special structure; */
/*             = 'T':  TC and TR are upper and lower triangular, */
/*                     respectively. Depending on the block sizes, two */
/*                     different shapes of the last blocks in TC and TR */
/*                     are possible, as illustrated below: */

/*                     1)    TC       TR     2)   TC         TR */

/*                          x x x    x 0 0      x x x x    x 0 0 0 */
/*                          0 x x    x x 0      0 x x x    x x 0 0 */
/*                          0 0 x    x x x      0 0 x x    x x x 0 */
/*                          0 0 0    x x x */

/*     Input/Output Parameters */

/*     K       (input) INTEGER */
/*             The number of rows in the blocks of T.  K >= 0. */

/*     L       (input) INTEGER */
/*             The number of columns in the blocks of T.  L >= 0. */

/*     M       (input) INTEGER */
/*             The number of blocks in the first block column of T. */
/*             M >= 1. */

/*     ML      (input) INTEGER */
/*             The lower block bandwidth, i.e., ML + 1 is the number of */
/*             nonzero blocks in the first block column of T. */
/*             0 <= ML < M and (ML + 1)*K >= L and */
/*             if ( M*K <= N*L ),  ML >= M - INT( ( M*K - 1 )/L ) - 1; */
/*                                 ML >= M - INT( M*K/L ) or */
/*                                 MOD( M*K, L ) >= K; */
/*             if ( M*K >= N*L ),  ML*K >= N*( L - K ). */

/*     N       (input) INTEGER */
/*             The number of blocks in the first block row of T. */
/*             N >= 1. */

/*     NU      (input) INTEGER */
/*             The upper block bandwidth, i.e., NU + 1 is the number of */
/*             nonzero blocks in the first block row of T. */
/*             If TRIU = 'N',   0 <= NU < N and */
/*                              (M + NU)*L >= MIN( M*K, N*L ); */
/*             if TRIU = 'T',   MAX(1-ML,0) <= NU < N and */
/*                              (M + NU)*L >= MIN( M*K, N*L ). */

/*     P       (input)  INTEGER */
/*             The number of previously computed block columns of R. */
/*             P*L < MIN( M*K,N*L ) + L and P >= 0. */

/*     S       (input)  INTEGER */
/*             The number of block columns of R to compute. */
/*             (P+S)*L < MIN( M*K,N*L ) + L and S >= 0. */

/*     TC      (input)  DOUBLE PRECISION array, dimension (LDTC,L) */
/*             On entry, if P = 0, the leading (ML+1)*K-by-L part of this */
/*             array must contain the nonzero blocks in the first block */
/*             column of T. */

/*     LDTC    INTEGER */
/*             The leading dimension of the array TC. */
/*             LDTC >= MAX(1,(ML+1)*K),  if P = 0. */

/*     TR      (input)  DOUBLE PRECISION array, dimension (LDTR,NU*L) */
/*             On entry, if P = 0, the leading K-by-NU*L part of this */
/*             array must contain the 2nd to the (NU+1)-st blocks of */
/*             the first block row of T. */

/*     LDTR    INTEGER */
/*             The leading dimension of the array TR. */
/*             LDTR >= MAX(1,K),  if P = 0. */

/*     RB      (output)  DOUBLE PRECISION array, dimension */
/*             (LDRB,MIN( S*L,MIN( M*K,N*L )-P*L )) */
/*             On exit, if INFO = 0 and TRIU = 'N', the leading */
/*             MIN( ML+NU+1,N )*L-by-MIN( S*L,MIN( M*K,N*L )-P*L ) part */
/*             of this array contains the (P+1)-th to (P+S)-th block */
/*             column of the lower R factor (1) in band storage format. */
/*             On exit, if INFO = 0 and TRIU = 'T', the leading */
/*             MIN( (ML+NU)*L+1,N*L )-by-MIN( S*L,MIN( M*K,N*L )-P*L ) */
/*             part of this array contains the (P+1)-th to (P+S)-th block */
/*             column of the lower R factor (1) in band storage format. */
/*             For further details regarding the band storage scheme see */
/*             the documentation of the LAPACK routine DPBTF2. */

/*     LDRB    INTEGER */
/*             The leading dimension of the array RB. */
/*             LDRB >= MAX( MIN( ML+NU+1,N )*L,1 ),      if TRIU = 'N'; */
/*             LDRB >= MAX( MIN( (ML+NU)*L+1,N*L ),1 ),  if TRIU = 'T'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -17,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */
/*             The first 1 + 2*MIN( ML+NU+1,N )*L*(K+L) elements of DWORK */
/*             should be preserved during successive calls of the routine. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             Let x = MIN( ML+NU+1,N ), then */
/*             LDWORK >= 1 + MAX( x*L*L + (2*NU+1)*L*K, */
/*                                2*x*L*(K+L) + (6+x)*L ),  if P = 0; */
/*             LDWORK >= 1 + 2*x*L*(K+L) + (6+x)*L,         if P > 0. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the full rank condition for the first MIN(M*K, N*L) */
/*                   columns of T is (numerically) violated. */

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

/*     The implemented method yields a factor R which has comparable */
/*     accuracy with the Cholesky factor of T^T * T. */
/*     The algorithm requires */
/*               2                                  2 */
/*           O( L *K*N*( ML + NU ) + N*( ML + NU )*L *( L + K ) ) */

/*     floating point operations. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, May 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, June 2001. */
/*     D. Kressner, Technical Univ. Berlin, Germany, July 2002. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2004. */

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

#line 237 "MB02HD.f"
    /* Parameter adjustments */
#line 237 "MB02HD.f"
    tc_dim1 = *ldtc;
#line 237 "MB02HD.f"
    tc_offset = 1 + tc_dim1;
#line 237 "MB02HD.f"
    tc -= tc_offset;
#line 237 "MB02HD.f"
    tr_dim1 = *ldtr;
#line 237 "MB02HD.f"
    tr_offset = 1 + tr_dim1;
#line 237 "MB02HD.f"
    tr -= tr_offset;
#line 237 "MB02HD.f"
    rb_dim1 = *ldrb;
#line 237 "MB02HD.f"
    rb_offset = 1 + rb_dim1;
#line 237 "MB02HD.f"
    rb -= rb_offset;
#line 237 "MB02HD.f"
    --dwork;
#line 237 "MB02HD.f"

#line 237 "MB02HD.f"
    /* Function Body */
#line 237 "MB02HD.f"
    *info = 0;
#line 238 "MB02HD.f"
    ltri = lsame_(triu, "T", (ftnlen)1, (ftnlen)1);
/* Computing MIN */
#line 239 "MB02HD.f"
    i__1 = *ml + *nu + 1;
#line 239 "MB02HD.f"
    x = min(i__1,*n);
#line 240 "MB02HD.f"
    lenr = x * *l;
#line 241 "MB02HD.f"
    if (ltri) {
/* Computing MIN */
#line 242 "MB02HD.f"
	i__1 = (*ml + *nu) * *l + 1, i__2 = *n * *l;
#line 242 "MB02HD.f"
	sizr = min(i__1,i__2);
#line 243 "MB02HD.f"
    } else {
#line 244 "MB02HD.f"
	sizr = lenr;
#line 245 "MB02HD.f"
    }
#line 246 "MB02HD.f"
    if (*p == 0) {
/* Computing MAX */
#line 247 "MB02HD.f"
	i__1 = lenr * *l + ((*nu << 1) + 1) * *l * *k, i__2 = (lenr << 1) * (*
		k + *l) + (x + 6) * *l;
#line 247 "MB02HD.f"
	wrkmin = max(i__1,i__2) + 1;
#line 249 "MB02HD.f"
    } else {
#line 250 "MB02HD.f"
	wrkmin = (lenr << 1) * (*k + *l) + 1 + (x + 6) * *l;
#line 251 "MB02HD.f"
    }
#line 252 "MB02HD.f"
    posr = 1;

/*     Check the scalar input parameters. */

#line 256 "MB02HD.f"
    if (! (ltri || lsame_(triu, "N", (ftnlen)1, (ftnlen)1))) {
#line 257 "MB02HD.f"
	*info = -1;
#line 258 "MB02HD.f"
    } else if (*k < 0) {
#line 259 "MB02HD.f"
	*info = -2;
#line 260 "MB02HD.f"
    } else if (*l < 0) {
#line 261 "MB02HD.f"
	*info = -3;
#line 262 "MB02HD.f"
    } else if (*m < 1) {
#line 263 "MB02HD.f"
	*info = -4;
#line 264 "MB02HD.f"
    } else if (*ml >= *m || (*ml + 1) * *k < *l || *m * *k <= *n * *l && (*ml 
	    < *m - (*m * *k - 1) / *l - 1 || *ml < *m - *m * *k / *l && *m * *
	    k % *l < *k) || *m * *k >= *n * *l && *ml * *k < *n * (*l - *k)) {
#line 268 "MB02HD.f"
	*info = -5;
#line 269 "MB02HD.f"
    } else if (*n < 1) {
#line 270 "MB02HD.f"
	*info = -6;
#line 271 "MB02HD.f"
    } else /* if(complicated condition) */ {
/* Computing MIN */
#line 271 "MB02HD.f"
	i__1 = *m * *k, i__2 = *n * *l;
#line 271 "MB02HD.f"
	if (*nu >= *n || *nu < 0 || ltri && *nu < 1 - *ml || (*m + *nu) * *l <
		 min(i__1,i__2)) {
#line 273 "MB02HD.f"
	    *info = -7;
#line 274 "MB02HD.f"
	} else /* if(complicated condition) */ {
/* Computing MIN */
#line 274 "MB02HD.f"
	    i__1 = *m * *k, i__2 = *n * *l;
#line 274 "MB02HD.f"
	    if (*p < 0 || *p * *l - *l >= min(i__1,i__2)) {
#line 275 "MB02HD.f"
		*info = -8;
#line 276 "MB02HD.f"
	    } else /* if(complicated condition) */ {
/* Computing MIN */
#line 276 "MB02HD.f"
		i__1 = *m * *k, i__2 = *n * *l;
#line 276 "MB02HD.f"
		if (*s < 0 || (*p + *s - 1) * *l >= min(i__1,i__2)) {
#line 277 "MB02HD.f"
		    *info = -9;
#line 278 "MB02HD.f"
		} else /* if(complicated condition) */ {
/* Computing MAX */
#line 278 "MB02HD.f"
		    i__1 = 1, i__2 = (*ml + 1) * *k;
#line 278 "MB02HD.f"
		    if (*p == 0 && *ldtc < max(i__1,i__2)) {
#line 279 "MB02HD.f"
			*info = -11;
#line 280 "MB02HD.f"
		    } else if (*p == 0 && *ldtr < max(1,*k)) {
#line 281 "MB02HD.f"
			*info = -13;
#line 282 "MB02HD.f"
		    } else if (*ldrb < max(sizr,1)) {
#line 283 "MB02HD.f"
			*info = 15;
#line 284 "MB02HD.f"
		    } else if (*ldwork < wrkmin) {
#line 285 "MB02HD.f"
			dwork[1] = (doublereal) wrkmin;
#line 286 "MB02HD.f"
			*info = -17;
#line 287 "MB02HD.f"
		    }
#line 287 "MB02HD.f"
		}
#line 287 "MB02HD.f"
	    }
#line 287 "MB02HD.f"
	}
#line 287 "MB02HD.f"
    }

/*     Return if there were illegal values. */

#line 291 "MB02HD.f"
    if (*info != 0) {
#line 292 "MB02HD.f"
	i__1 = -(*info);
#line 292 "MB02HD.f"
	xerbla_("MB02HD", &i__1, (ftnlen)6);
#line 293 "MB02HD.f"
	return 0;
#line 294 "MB02HD.f"
    }

/*     Quick return if possible. */

#line 298 "MB02HD.f"
    if (*l * *k * *s == 0) {
#line 299 "MB02HD.f"
	dwork[1] = 1.;
#line 300 "MB02HD.f"
	return 0;
#line 301 "MB02HD.f"
    }

#line 303 "MB02HD.f"
    wrkopt = 1;

/*     Compute the generator if P = 0. */

#line 307 "MB02HD.f"
    if (*p == 0) {

/*        1st column of the generator. */

#line 311 "MB02HD.f"
	lenc = (*ml + 1) * *k;
/* Computing MAX */
/* Computing MIN */
#line 312 "MB02HD.f"
	i__2 = *nu, i__3 = *n - *m;
#line 312 "MB02HD.f"
	i__1 = *ml + 1 + min(i__2,i__3);
#line 312 "MB02HD.f"
	lenl = max(i__1,0);
#line 313 "MB02HD.f"
	pdc = lenr * *l + 1;
#line 314 "MB02HD.f"
	pdw = pdc + lenc * *l;

/*        QR decomposition of the nonzero blocks in TC. */

#line 318 "MB02HD.f"
	dlacpy_("All", &lenc, l, &tc[tc_offset], ldtc, &dwork[pdc + 1], &lenc,
		 (ftnlen)3);
#line 319 "MB02HD.f"
	i__1 = *ldwork - pdw - *l;
#line 319 "MB02HD.f"
	dgeqrf_(&lenc, l, &dwork[pdc + 1], &lenc, &dwork[pdw + 1], &dwork[pdw 
		+ *l + 1], &i__1, &ierr);
/* Computing MAX */
#line 321 "MB02HD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + *l + 1] + pdw + *l;
#line 321 "MB02HD.f"
	wrkopt = max(i__1,i__2);

/*        The R factor is the transposed of the first block in the */
/*        generator. */

#line 326 "MB02HD.f"
	ma02ad_("Upper part", l, l, &dwork[pdc + 1], &lenc, &dwork[2], &lenr, 
		(ftnlen)10);

/*        Get the first block column of the Q factor. */

#line 331 "MB02HD.f"
	i__1 = *ldwork - pdw - *l;
#line 331 "MB02HD.f"
	dorgqr_(&lenc, l, l, &dwork[pdc + 1], &lenc, &dwork[pdw + 1], &dwork[
		pdw + *l + 1], &i__1, &ierr);
/* Computing MAX */
#line 333 "MB02HD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + *l + 1] + pdw + *l;
#line 333 "MB02HD.f"
	wrkopt = max(i__1,i__2);

/*        Construct a flipped copy of TC for faster multiplication. */

#line 337 "MB02HD.f"
	pt = lenc - (*k << 1) + 1;

#line 339 "MB02HD.f"
	i__1 = pdw + *ml * *k * *l;
#line 339 "MB02HD.f"
	i__2 = *k * *l;
#line 339 "MB02HD.f"
	for (i__ = pdw + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
		 {
#line 340 "MB02HD.f"
	    dlacpy_("All", k, l, &tc[pt + tc_dim1], ldtc, &dwork[i__], k, (
		    ftnlen)3);
#line 341 "MB02HD.f"
	    pt -= *k;
#line 342 "MB02HD.f"
/* L10: */
#line 342 "MB02HD.f"
	}

/*        Multiply T^T with the first block column of Q. */

#line 346 "MB02HD.f"
	pdw = i__;
#line 347 "MB02HD.f"
	pdr = *l + 2;
#line 348 "MB02HD.f"
	len = *nu * *l;
#line 349 "MB02HD.f"
	i__2 = lenr - *l;
#line 349 "MB02HD.f"
	dlaset_("All", &i__2, l, &c_b10, &c_b10, &dwork[pdr], &lenr, (ftnlen)
		3);

#line 351 "MB02HD.f"
	i__2 = *ml + 1;
#line 351 "MB02HD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MIN */
#line 352 "MB02HD.f"
	    i__3 = i__ - 1, i__4 = *n - 1;
#line 352 "MB02HD.f"
	    i__1 = min(i__3,i__4) * *l;
#line 352 "MB02HD.f"
	    dgemm_("Transpose", "NonTranspose", &i__1, l, k, &c_b15, &dwork[
		    pdw], k, &dwork[pdc + 1], &lenc, &c_b15, &dwork[pdr], &
		    lenr, (ftnlen)9, (ftnlen)12);
#line 355 "MB02HD.f"
	    if (len > 0) {
#line 356 "MB02HD.f"
		dgemm_("Transpose", "NonTranspose", &len, l, k, &c_b15, &tr[
			tr_offset], ldtr, &dwork[pdc + 1], &lenc, &c_b15, &
			dwork[pdr + (i__ - 1) * *l], &lenr, (ftnlen)9, (
			ftnlen)12);
#line 359 "MB02HD.f"
	    }
#line 360 "MB02HD.f"
	    pdw -= *k * *l;
#line 361 "MB02HD.f"
	    pdc += *k;
#line 362 "MB02HD.f"
	    if (i__ >= *n - *nu) {
#line 362 "MB02HD.f"
		len -= *l;
#line 362 "MB02HD.f"
	    }
#line 363 "MB02HD.f"
/* L20: */
#line 363 "MB02HD.f"
	}

/*        Copy the first block column to R. */

#line 367 "MB02HD.f"
	if (ltri) {

#line 369 "MB02HD.f"
	    i__2 = *l;
#line 369 "MB02HD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MIN */
#line 370 "MB02HD.f"
		i__3 = sizr, i__4 = *n * *l - i__ + 1;
#line 370 "MB02HD.f"
		i__1 = min(i__3,i__4);
#line 370 "MB02HD.f"
		dcopy_(&i__1, &dwork[(i__ - 1) * lenr + i__ + 1], &c__1, &rb[
			posr * rb_dim1 + 1], &c__1);
#line 373 "MB02HD.f"
		++posr;
#line 374 "MB02HD.f"
/* L30: */
#line 374 "MB02HD.f"
	    }

#line 376 "MB02HD.f"
	} else {

#line 378 "MB02HD.f"
	    i__2 = *l;
#line 378 "MB02HD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 379 "MB02HD.f"
		i__1 = lenr - i__ + 1;
#line 379 "MB02HD.f"
		dcopy_(&i__1, &dwork[(i__ - 1) * lenr + i__ + 1], &c__1, &rb[
			posr * rb_dim1 + 1], &c__1);
#line 381 "MB02HD.f"
		if (lenr < *n * *l && i__ > 1) {
#line 382 "MB02HD.f"
		    i__1 = i__ - 1;
#line 382 "MB02HD.f"
		    dlaset_("All", &i__1, &c__1, &c_b10, &c_b10, &rb[lenr - 
			    i__ + 2 + posr * rb_dim1], ldrb, (ftnlen)3);
#line 384 "MB02HD.f"
		}
#line 385 "MB02HD.f"
		++posr;
#line 386 "MB02HD.f"
/* L40: */
#line 386 "MB02HD.f"
	    }

#line 388 "MB02HD.f"
	}

/*        Quick return if N = 1. */

#line 392 "MB02HD.f"
	if (*n == 1) {
#line 393 "MB02HD.f"
	    dwork[1] = (doublereal) wrkopt;
#line 394 "MB02HD.f"
	    return 0;
#line 395 "MB02HD.f"
	}

/*        2nd column of the generator. */

#line 399 "MB02HD.f"
	pdr = lenr * *l + 1;
#line 400 "MB02HD.f"
	i__2 = *nu * *l;
#line 400 "MB02HD.f"
	ma02ad_("All", k, &i__2, &tr[tr_offset], ldtr, &dwork[pdr + 1], &lenr,
		 (ftnlen)3);
#line 401 "MB02HD.f"
	i__2 = lenr - *nu * *l;
#line 401 "MB02HD.f"
	dlaset_("All", &i__2, k, &c_b10, &c_b10, &dwork[pdr + *nu * *l + 1], &
		lenr, (ftnlen)3);

/*        3rd column of the generator. */

#line 406 "MB02HD.f"
	pnr = pdr + lenr * *k;
#line 407 "MB02HD.f"
	i__2 = lenr - *l;
#line 407 "MB02HD.f"
	dlacpy_("All", &i__2, l, &dwork[*l + 2], &lenr, &dwork[pnr + 1], &
		lenr, (ftnlen)3);
#line 409 "MB02HD.f"
	dlaset_("All", l, l, &c_b10, &c_b10, &dwork[pnr + lenr - *l + 1], &
		lenr, (ftnlen)3);

/*        4th column of the generator. */

#line 414 "MB02HD.f"
	pfr = pnr + lenr * *l;

#line 416 "MB02HD.f"
	pdw = pfr + (*m - *ml - 1) * *l % lenr;
#line 417 "MB02HD.f"
	pt = *ml * *k + 1;
/* Computing MIN */
#line 418 "MB02HD.f"
	i__1 = *ml + 1;
#line 418 "MB02HD.f"
	i__2 = min(i__1,lenl);
#line 418 "MB02HD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 419 "MB02HD.f"
	    ma02ad_("All", k, l, &tc[pt + tc_dim1], ldtc, &dwork[pdw + 1], &
		    lenr, (ftnlen)3);
#line 421 "MB02HD.f"
	    pt -= *k;
#line 422 "MB02HD.f"
	    pdw = pfr + (pdw + *l - pfr) % lenr;
#line 423 "MB02HD.f"
/* L50: */
#line 423 "MB02HD.f"
	}
#line 424 "MB02HD.f"
	pt = 1;
#line 425 "MB02HD.f"
	i__2 = lenl;
#line 425 "MB02HD.f"
	for (i__ = *ml + 2; i__ <= i__2; ++i__) {
#line 426 "MB02HD.f"
	    ma02ad_("All", k, l, &tr[pt * tr_dim1 + 1], ldtr, &dwork[pdw + 1],
		     &lenr, (ftnlen)3);
#line 428 "MB02HD.f"
	    pt += *l;
#line 429 "MB02HD.f"
	    pdw = pfr + (pdw + *l - pfr) % lenr;
#line 430 "MB02HD.f"
/* L60: */
#line 430 "MB02HD.f"
	}
#line 431 "MB02HD.f"
	pre = 1;
#line 432 "MB02HD.f"
	stps = *s - 1;
#line 433 "MB02HD.f"
    } else {
#line 434 "MB02HD.f"
	pdr = lenr * *l + 1;
#line 435 "MB02HD.f"
	pnr = pdr + lenr * *k;
#line 436 "MB02HD.f"
	pfr = pnr + lenr * *l;
#line 437 "MB02HD.f"
	pre = *p;
#line 438 "MB02HD.f"
	stps = *s;
#line 439 "MB02HD.f"
    }

#line 441 "MB02HD.f"
    pdw = pfr + lenr * *k;
#line 442 "MB02HD.f"
    head = (pre - 1) * *l % lenr;

/*     Determine block size for the involved block Householder */
/*     transformations. */

/* Computing MIN */
#line 447 "MB02HD.f"
    i__2 = ilaenv_(&c__1, "DGELQF", " ", &lenr, l, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 447 "MB02HD.f"
    nb = min(i__2,*l);
#line 448 "MB02HD.f"
    kk = pdw + *l * 6;
/* Computing MAX */
#line 449 "MB02HD.f"
    i__2 = wrkopt, i__1 = kk + lenr * nb;
#line 449 "MB02HD.f"
    wrkopt = max(i__2,i__1);
#line 450 "MB02HD.f"
    kk = *ldwork - kk;
#line 451 "MB02HD.f"
    if (kk < lenr * nb) {
#line 451 "MB02HD.f"
	nb = kk / lenr;
#line 451 "MB02HD.f"
    }
/* Computing MAX */
#line 452 "MB02HD.f"
    i__2 = 2, i__1 = ilaenv_(&c__2, "DGELQF", " ", &lenr, l, &c_n1, &c_n1, (
	    ftnlen)6, (ftnlen)1);
#line 452 "MB02HD.f"
    nbmin = max(i__2,i__1);
#line 453 "MB02HD.f"
    if (nb < nbmin) {
#line 453 "MB02HD.f"
	nb = 0;
#line 453 "MB02HD.f"
    }

/*     Generator reduction process. */

#line 457 "MB02HD.f"
    i__2 = pre + stps - 1;
#line 457 "MB02HD.f"
    for (i__ = pre; i__ <= i__2; ++i__) {

/*        The 4th generator column is not used in the first (M-ML) steps. */

#line 461 "MB02HD.f"
	if (i__ < *m - *ml) {
#line 462 "MB02HD.f"
	    col2 = *l;
#line 463 "MB02HD.f"
	} else {
#line 464 "MB02HD.f"
	    col2 = *k + *l;
#line 465 "MB02HD.f"
	}

/* Computing MIN */
#line 467 "MB02HD.f"
	i__1 = *l, i__3 = *m * *k - i__ * *l;
#line 467 "MB02HD.f"
	kk = min(i__1,i__3);
#line 468 "MB02HD.f"
	i__1 = kk + *k;
#line 468 "MB02HD.f"
	i__3 = *ldwork - pdw - *l * 6;
#line 468 "MB02HD.f"
	mb02cu_("Column", &kk, &i__1, &col2, &nb, &dwork[2], &lenr, &dwork[
		pdr + head + 1], &lenr, &dwork[pnr + head + 1], &lenr, &rnk, 
		ipvt, &dwork[pdw + 1], &c_b10, &dwork[pdw + *l * 6 + 1], &
		i__3, &ierr, (ftnlen)6);
#line 472 "MB02HD.f"
	if (ierr != 0) {

/*           Error return:  The rank condition is (numerically) not */
/*                          satisfied. */

#line 477 "MB02HD.f"
	    *info = 1;
#line 478 "MB02HD.f"
	    return 0;
#line 479 "MB02HD.f"
	}

/* Computing MAX */
/* Computing MIN */
#line 481 "MB02HD.f"
	i__3 = (*n - i__) * *l - kk, i__4 = lenr - head - kk;
#line 481 "MB02HD.f"
	i__1 = min(i__3,i__4);
#line 481 "MB02HD.f"
	len = max(i__1,0);
/* Computing MAX */
/* Computing MIN */
#line 482 "MB02HD.f"
	i__3 = (*n - i__) * *l - len - kk;
#line 482 "MB02HD.f"
	i__1 = min(i__3,head);
#line 482 "MB02HD.f"
	len2 = max(i__1,0);
#line 483 "MB02HD.f"
	if (len == lenr - kk) {
#line 484 "MB02HD.f"
	    *(unsigned char *)struct__ = *(unsigned char *)triu;
#line 485 "MB02HD.f"
	} else {
#line 486 "MB02HD.f"
	    *(unsigned char *)struct__ = 'N';
#line 487 "MB02HD.f"
	}
#line 488 "MB02HD.f"
	i__1 = kk + *k;
#line 488 "MB02HD.f"
	i__3 = *ldwork - pdw - *l * 6;
#line 488 "MB02HD.f"
	mb02cv_("Column", struct__, &kk, &len, &i__1, &col2, &nb, &c_n1, &
		dwork[2], &lenr, &dwork[pdr + head + 1], &lenr, &dwork[pnr + 
		head + 1], &lenr, &dwork[kk + 2], &lenr, &dwork[pdr + head + 
		kk + 1], &lenr, &dwork[pnr + head + kk + 1], &lenr, &dwork[
		pdw + 1], &dwork[pdw + *l * 6 + 1], &i__3, &ierr, (ftnlen)6, (
		ftnlen)1);

#line 495 "MB02HD.f"
	if ((*n - i__) * *l >= lenr) {
#line 496 "MB02HD.f"
	    *(unsigned char *)struct__ = *(unsigned char *)triu;
#line 497 "MB02HD.f"
	} else {
#line 498 "MB02HD.f"
	    *(unsigned char *)struct__ = 'N';
#line 499 "MB02HD.f"
	}

#line 501 "MB02HD.f"
	i__1 = kk + *k;
#line 501 "MB02HD.f"
	i__3 = *ldwork - pdw - *l * 6;
#line 501 "MB02HD.f"
	mb02cv_("Column", struct__, &kk, &len2, &i__1, &col2, &nb, &c_n1, &
		dwork[2], &lenr, &dwork[pdr + head + 1], &lenr, &dwork[pnr + 
		head + 1], &lenr, &dwork[kk + len + 2], &lenr, &dwork[pdr + 1]
		, &lenr, &dwork[pnr + 1], &lenr, &dwork[pdw + 1], &dwork[pdw 
		+ *l * 6 + 1], &i__3, &ierr, (ftnlen)6, (ftnlen)1);

#line 508 "MB02HD.f"
	i__1 = *k + col2;
#line 508 "MB02HD.f"
	dlaset_("All", l, &i__1, &c_b10, &c_b10, &dwork[pdr + head + 1], &
		lenr, (ftnlen)3);

/*        Copy current block column to R. */

#line 513 "MB02HD.f"
	if (ltri) {

#line 515 "MB02HD.f"
	    i__1 = kk;
#line 515 "MB02HD.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 516 "MB02HD.f"
		i__4 = sizr, i__5 = (*n - i__) * *l - j + 1;
#line 516 "MB02HD.f"
		i__3 = min(i__4,i__5);
#line 516 "MB02HD.f"
		dcopy_(&i__3, &dwork[(j - 1) * lenr + j + 1], &c__1, &rb[posr 
			* rb_dim1 + 1], &c__1);
#line 519 "MB02HD.f"
		++posr;
#line 520 "MB02HD.f"
/* L70: */
#line 520 "MB02HD.f"
	    }

#line 522 "MB02HD.f"
	} else {

#line 524 "MB02HD.f"
	    i__1 = kk;
#line 524 "MB02HD.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 525 "MB02HD.f"
		i__4 = sizr - j + 1, i__5 = (*n - i__) * *l - j + 1;
#line 525 "MB02HD.f"
		i__3 = min(i__4,i__5);
#line 525 "MB02HD.f"
		dcopy_(&i__3, &dwork[(j - 1) * lenr + j + 1], &c__1, &rb[posr 
			* rb_dim1 + 1], &c__1);
#line 528 "MB02HD.f"
		if (lenr < (*n - i__) * *l && j > 1) {
#line 529 "MB02HD.f"
		    i__3 = j - 1;
/* Computing MIN */
#line 529 "MB02HD.f"
		    i__4 = sizr - j + 1, i__5 = (*n - i__) * *l - j + 1;
#line 529 "MB02HD.f"
		    dlaset_("All", &i__3, &c__1, &c_b10, &c_b10, &rb[min(i__4,
			    i__5) + 1 + posr * rb_dim1], ldrb, (ftnlen)3);
#line 532 "MB02HD.f"
		}
#line 533 "MB02HD.f"
		++posr;
#line 534 "MB02HD.f"
/* L80: */
#line 534 "MB02HD.f"
	    }

#line 536 "MB02HD.f"
	}

#line 538 "MB02HD.f"
	head = (head + *l) % lenr;
#line 539 "MB02HD.f"
/* L90: */
#line 539 "MB02HD.f"
    }

#line 541 "MB02HD.f"
    dwork[1] = (doublereal) wrkopt;
#line 542 "MB02HD.f"
    return 0;

/* *** Last line of MB02HD *** */
} /* mb02hd_ */

