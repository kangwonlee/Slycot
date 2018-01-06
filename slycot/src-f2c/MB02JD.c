#line 1 "MB02JD.f"
/* MB02JD.f -- translated by f2c (version 20100827).
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

#line 1 "MB02JD.f"
/* Table of constant values */

static doublereal c_b10 = 1.;
static doublereal c_b11 = 0.;
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* Subroutine */ int mb02jd_(char *job, integer *k, integer *l, integer *m, 
	integer *n, integer *p, integer *s, doublereal *tc, integer *ldtc, 
	doublereal *tr, integer *ldtr, doublereal *q, integer *ldq, 
	doublereal *r__, integer *ldr, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen job_len)
{
    /* System generated locals */
    integer q_dim1, q_offset, r_dim1, r_offset, tc_dim1, tc_offset, tr_dim1, 
	    tr_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9;

    /* Local variables */
    static integer i__, nb, kk, pt, len, pdq, pre, pdw, rnk, pnq, pnr, colr, 
	    ierr, shfr, ipvt[1], stps;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    mb02kd_(char *, char *, integer *, integer *, integer *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    mb02cu_(char *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen), mb02cv_(char *, char 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin;
    static logical compq;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dorgqr_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);
    static integer wrkmin, wrkopt;


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

/*     To compute a lower triangular matrix R and a matrix Q with */
/*     Q^T Q = I such that */
/*                                    T */
/*                           T  =  Q R , */

/*     where T is a K*M-by-L*N block Toeplitz matrix with blocks of size */
/*     (K,L). The first column of T will be denoted by TC and the first */
/*     row by TR. It is assumed that the first MIN(M*K, N*L) columns of T */
/*     have full rank. */

/*     By subsequent calls of this routine the factors Q and R can be */
/*     computed block column by block column. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the output of the routine as follows: */
/*             = 'Q':  computes Q and R; */
/*             = 'R':  only computes R. */

/*     Input/Output Parameters */

/*     K       (input)  INTEGER */
/*             The number of rows in one block of T.  K >= 0. */

/*     L       (input)  INTEGER */
/*             The number of columns in one block of T.  L >= 0. */

/*     M       (input)  INTEGER */
/*             The number of blocks in one block column of T.  M >= 0. */

/*     N       (input)  INTEGER */
/*             The number of blocks in one block row of T.  N >= 0. */

/*     P       (input)  INTEGER */
/*             The number of previously computed block columns of R. */
/*             P*L < MIN( M*K,N*L ) + L and P >= 0. */

/*     S       (input)  INTEGER */
/*             The number of block columns of R to compute. */
/*             (P+S)*L < MIN( M*K,N*L ) + L and S >= 0. */

/*     TC      (input) DOUBLE PRECISION array, dimension (LDTC, L) */
/*             On entry, if P = 0, the leading M*K-by-L part of this */
/*             array must contain the first block column of T. */

/*     LDTC    INTEGER */
/*             The leading dimension of the array TC. */
/*             LDTC >= MAX(1,M*K). */

/*     TR      (input)  DOUBLE PRECISION array, dimension (LDTR,(N-1)*L) */
/*             On entry, if P = 0, the leading K-by-(N-1)*L part of this */
/*             array must contain the first block row of T without the */
/*             leading K-by-L block. */

/*     LDTR    INTEGER */
/*             The leading dimension of the array TR. */
/*             LDTR >= MAX(1,K). */

/*     Q       (input/output)  DOUBLE PRECISION array, dimension */
/*                             (LDQ,MIN( S*L, MIN( M*K,N*L )-P*L )) */
/*             On entry, if JOB = 'Q'  and  P > 0, the leading M*K-by-L */
/*             part of this array must contain the last block column of Q */
/*             from a previous call of this routine. */
/*             On exit, if JOB = 'Q'  and  INFO = 0, the leading */
/*             M*K-by-MIN( S*L, MIN( M*K,N*L )-P*L ) part of this array */
/*             contains the P-th to (P+S)-th block columns of the factor */
/*             Q. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q. */
/*             LDQ >= MAX(1,M*K), if JOB = 'Q'; */
/*             LDQ >= 1,          if JOB = 'R'. */

/*     R       (input/output)  DOUBLE PRECISION array, dimension */
/*                             (LDR,MIN( S*L, MIN( M*K,N*L )-P*L )) */
/*             On entry, if P > 0, the leading (N-P+1)*L-by-L */
/*             part of this array must contain the nozero part of the */
/*             last block column of R from a previous call of this */
/*             routine. */
/*             One exit, if INFO = 0, the leading */
/*             MIN( N, N-P+1 )*L-by-MIN( S*L, MIN( M*K,N*L )-P*L ) */
/*             part of this array contains the nonzero parts of the P-th */
/*             to (P+S)-th block columns of the lower triangular */
/*             factor R. */
/*             Note that elements in the strictly upper triangular part */
/*             will not be referenced. */

/*     LDR     INTEGER */
/*             The leading dimension of the array R. */
/*             LDR >= MAX( 1, MIN( N, N-P+1 )*L ) */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */
/*             On exit, if INFO = -17,  DWORK(1) returns the minimum */
/*             value of LDWORK. */
/*             If JOB = 'Q', the first 1 + ( (N-1)*L + M*K )*( 2*K + L ) */
/*             elements of DWORK should be preserved during successive */
/*             calls of the routine. */
/*             If JOB = 'R', the first 1 + (N-1)*L*( 2*K + L ) elements */
/*             of DWORK should be preserved during successive calls of */
/*             the routine. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             JOB = 'Q': */
/*                LDWORK >= 1 + ( M*K + ( N - 1 )*L )*( L + 2*K ) + 6*L */
/*                            + MAX( M*K,( N - MAX( 1,P )*L ) ); */
/*             JOB = 'R': */
/*                If P = 0, */
/*                   LDWORK >= MAX( 1 + ( N - 1 )*L*( L + 2*K ) + 6*L */
/*                                    + (N-1)*L, M*K*( L + 1 ) + L ); */
/*                If P > 0, */
/*                   LDWORK >= 1 + (N-1)*L*( L + 2*K ) + 6*L + (N-P)*L. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the full rank condition for the first MIN(M*K, N*L) */
/*                   columns of T is (numerically) violated. */

/*     METHOD */

/*     Block Householder transformations and modified hyperbolic */
/*     rotations are used in the Schur algorithm [1], [2]. */

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
/*     accuracy with the Cholesky factor of T^T * T. Q is implicitly */
/*     computed from the formula Q = T * inv(R^T R) * R, i.e., for ill */
/*     conditioned problems this factor is of very limited value. */
/*                                 2 */
/*     The algorithm requires 0(K*L *M*N) floating point operations. */

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

#line 227 "MB02JD.f"
    /* Parameter adjustments */
#line 227 "MB02JD.f"
    tc_dim1 = *ldtc;
#line 227 "MB02JD.f"
    tc_offset = 1 + tc_dim1;
#line 227 "MB02JD.f"
    tc -= tc_offset;
#line 227 "MB02JD.f"
    tr_dim1 = *ldtr;
#line 227 "MB02JD.f"
    tr_offset = 1 + tr_dim1;
#line 227 "MB02JD.f"
    tr -= tr_offset;
#line 227 "MB02JD.f"
    q_dim1 = *ldq;
#line 227 "MB02JD.f"
    q_offset = 1 + q_dim1;
#line 227 "MB02JD.f"
    q -= q_offset;
#line 227 "MB02JD.f"
    r_dim1 = *ldr;
#line 227 "MB02JD.f"
    r_offset = 1 + r_dim1;
#line 227 "MB02JD.f"
    r__ -= r_offset;
#line 227 "MB02JD.f"
    --dwork;
#line 227 "MB02JD.f"

#line 227 "MB02JD.f"
    /* Function Body */
#line 227 "MB02JD.f"
    *info = 0;
#line 228 "MB02JD.f"
    compq = lsame_(job, "Q", (ftnlen)1, (ftnlen)1);
#line 229 "MB02JD.f"
    if (compq) {
/* Computing MAX */
#line 230 "MB02JD.f"
	i__1 = *m * *k, i__2 = (*n - max(1,*p)) * *l;
#line 230 "MB02JD.f"
	wrkmin = (*m * *k + (*n - 1) * *l) * (*l + (*k << 1)) + 1 + *l * 6 + 
		max(i__1,i__2);
#line 232 "MB02JD.f"
    } else {
#line 233 "MB02JD.f"
	wrkmin = (*n - 1) * *l * (*l + (*k << 1)) + 1 + *l * 6 + (*n - max(*p,
		1)) * *l;
#line 235 "MB02JD.f"
	if (*p == 0) {
/* Computing MAX */
#line 236 "MB02JD.f"
	    i__1 = wrkmin, i__2 = *m * *k * (*l + 1) + *l;
#line 236 "MB02JD.f"
	    wrkmin = max(i__1,i__2);
#line 237 "MB02JD.f"
	}
#line 238 "MB02JD.f"
    }

/*     Check the scalar input parameters. */

#line 242 "MB02JD.f"
    if (! (compq || lsame_(job, "R", (ftnlen)1, (ftnlen)1))) {
#line 243 "MB02JD.f"
	*info = -1;
#line 244 "MB02JD.f"
    } else if (*k < 0) {
#line 245 "MB02JD.f"
	*info = -2;
#line 246 "MB02JD.f"
    } else if (*l < 0) {
#line 247 "MB02JD.f"
	*info = -3;
#line 248 "MB02JD.f"
    } else if (*m < 0) {
#line 249 "MB02JD.f"
	*info = -4;
#line 250 "MB02JD.f"
    } else if (*n < 0) {
#line 251 "MB02JD.f"
	*info = -5;
#line 252 "MB02JD.f"
    } else /* if(complicated condition) */ {
/* Computing MIN */
#line 252 "MB02JD.f"
	i__1 = *m * *k, i__2 = *n * *l;
#line 252 "MB02JD.f"
	if (*p * *l >= min(i__1,i__2) + *l || *p < 0) {
#line 253 "MB02JD.f"
	    *info = -6;
#line 254 "MB02JD.f"
	} else /* if(complicated condition) */ {
/* Computing MIN */
#line 254 "MB02JD.f"
	    i__1 = *m * *k, i__2 = *n * *l;
#line 254 "MB02JD.f"
	    if ((*p + *s) * *l >= min(i__1,i__2) + *l || *s < 0) {
#line 255 "MB02JD.f"
		*info = -7;
#line 256 "MB02JD.f"
	    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 256 "MB02JD.f"
		i__1 = 1, i__2 = *m * *k;
#line 256 "MB02JD.f"
		if (*ldtc < max(i__1,i__2)) {
#line 257 "MB02JD.f"
		    *info = -9;
#line 258 "MB02JD.f"
		} else if (*ldtr < max(1,*k)) {
#line 259 "MB02JD.f"
		    *info = -11;
#line 260 "MB02JD.f"
		} else if (*ldq < 1 || compq && *ldq < *m * *k) {
#line 261 "MB02JD.f"
		    *info = -13;
#line 262 "MB02JD.f"
		} else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MIN */
#line 262 "MB02JD.f"
		    i__3 = *n, i__4 = *n - *p + 1;
#line 262 "MB02JD.f"
		    i__1 = 1, i__2 = min(i__3,i__4) * *l;
#line 262 "MB02JD.f"
		    if (*ldr < max(i__1,i__2)) {
#line 263 "MB02JD.f"
			*info = -15;
#line 264 "MB02JD.f"
		    } else if (*ldwork < wrkmin) {
#line 265 "MB02JD.f"
			dwork[1] = (doublereal) wrkmin;
#line 266 "MB02JD.f"
			*info = -17;
#line 267 "MB02JD.f"
		    }
#line 267 "MB02JD.f"
		}
#line 267 "MB02JD.f"
	    }
#line 267 "MB02JD.f"
	}
#line 267 "MB02JD.f"
    }

/*     Return if there were illegal values. */

#line 271 "MB02JD.f"
    if (*info != 0) {
#line 272 "MB02JD.f"
	i__1 = -(*info);
#line 272 "MB02JD.f"
	xerbla_("MB02JD", &i__1, (ftnlen)6);
#line 273 "MB02JD.f"
	return 0;
#line 274 "MB02JD.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 278 "MB02JD.f"
    i__1 = min(*m,*n), i__2 = *k * *l, i__1 = min(i__1,i__2);
#line 278 "MB02JD.f"
    if (min(i__1,*s) == 0) {
#line 279 "MB02JD.f"
	dwork[1] = 1.;
#line 280 "MB02JD.f"
	return 0;
#line 281 "MB02JD.f"
    }

/*     Catch M*K <= L. */

#line 285 "MB02JD.f"
    wrkopt = 1;
#line 286 "MB02JD.f"
    if (*m * *k <= *l) {
#line 287 "MB02JD.f"
	i__1 = *m * *k;
#line 287 "MB02JD.f"
	i__2 = *m * *k;
#line 287 "MB02JD.f"
	dlacpy_("All", &i__1, l, &tc[tc_offset], ldtc, &dwork[1], &i__2, (
		ftnlen)3);
#line 288 "MB02JD.f"
	pdw = *m * *k * *l + 1;
#line 289 "MB02JD.f"
	i__1 = *m * *k;
#line 289 "MB02JD.f"
	i__2 = *m * *k;
#line 289 "MB02JD.f"
	i__3 = *ldwork - pdw - *m * *k + 1;
#line 289 "MB02JD.f"
	dgeqrf_(&i__1, l, &dwork[1], &i__2, &dwork[pdw], &dwork[pdw + *m * *k]
		, &i__3, &ierr);
/* Computing MAX */
#line 291 "MB02JD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + *m * *k] + pdw + *m * *k 
		- 1;
#line 291 "MB02JD.f"
	wrkopt = max(i__1,i__2);
#line 292 "MB02JD.f"
	i__1 = *m * *k;
#line 292 "MB02JD.f"
	i__2 = *m * *k;
#line 292 "MB02JD.f"
	ma02ad_("Upper part", &i__1, l, &dwork[1], &i__2, &r__[r_offset], ldr,
		 (ftnlen)10);
#line 293 "MB02JD.f"
	i__1 = *m * *k;
#line 293 "MB02JD.f"
	i__2 = *m * *k;
#line 293 "MB02JD.f"
	i__3 = *m * *k;
#line 293 "MB02JD.f"
	i__4 = *m * *k;
#line 293 "MB02JD.f"
	i__5 = *ldwork - pdw - *m * *k + 1;
#line 293 "MB02JD.f"
	dorgqr_(&i__1, &i__2, &i__3, &dwork[1], &i__4, &dwork[pdw], &dwork[
		pdw + *m * *k], &i__5, &ierr);
/* Computing MAX */
#line 295 "MB02JD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + *m * *k] + pdw + *m * *k 
		- 1;
#line 295 "MB02JD.f"
	wrkopt = max(i__1,i__2);
#line 296 "MB02JD.f"
	if (compq) {
#line 297 "MB02JD.f"
	    i__1 = *m * *k;
#line 297 "MB02JD.f"
	    i__2 = *m * *k;
#line 297 "MB02JD.f"
	    i__3 = *m * *k;
#line 297 "MB02JD.f"
	    dlacpy_("All", &i__1, &i__2, &dwork[1], &i__3, &q[q_offset], ldq, 
		    (ftnlen)3);
#line 298 "MB02JD.f"
	}
#line 299 "MB02JD.f"
	pdw = *m * *k * *m * *k + 1;
#line 300 "MB02JD.f"
	if (*n > 1) {
#line 301 "MB02JD.f"
	    i__1 = *n - 1;
#line 301 "MB02JD.f"
	    i__2 = *m * *k;
#line 301 "MB02JD.f"
	    i__3 = *m * *k;
#line 301 "MB02JD.f"
	    i__4 = *ldwork - pdw + 1;
#line 301 "MB02JD.f"
	    mb02kd_("Row", "Transpose", k, l, m, &i__1, &i__2, &c_b10, &c_b11,
		     &tc[tc_offset], ldtc, &tr[tr_offset], ldtr, &dwork[1], &
		    i__3, &r__[*l + 1 + r_dim1], ldr, &dwork[pdw], &i__4, &
		    ierr, (ftnlen)3, (ftnlen)9);
#line 304 "MB02JD.f"
	}
/* Computing MAX */
#line 305 "MB02JD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
#line 305 "MB02JD.f"
	wrkopt = max(i__1,i__2);
#line 306 "MB02JD.f"
	dwork[1] = (doublereal) wrkopt;
#line 307 "MB02JD.f"
	return 0;
#line 308 "MB02JD.f"
    }

/*     Compute the generator if P = 0. */

#line 312 "MB02JD.f"
    if (*p == 0) {

/*        1st column of the generator. */

#line 316 "MB02JD.f"
	if (compq) {
#line 317 "MB02JD.f"
	    i__1 = *m * *k;
#line 317 "MB02JD.f"
	    dlacpy_("All", &i__1, l, &tc[tc_offset], ldtc, &q[q_offset], ldq, 
		    (ftnlen)3);
#line 318 "MB02JD.f"
	    i__1 = *m * *k;
#line 318 "MB02JD.f"
	    i__2 = *ldwork - *l;
#line 318 "MB02JD.f"
	    dgeqrf_(&i__1, l, &q[q_offset], ldq, &dwork[1], &dwork[*l + 1], &
		    i__2, &ierr);
/* Computing MAX */
#line 320 "MB02JD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[*l + 1] + *l;
#line 320 "MB02JD.f"
	    wrkopt = max(i__1,i__2);
#line 321 "MB02JD.f"
	    ma02ad_("Upper part", l, l, &q[q_offset], ldq, &r__[r_offset], 
		    ldr, (ftnlen)10);
#line 322 "MB02JD.f"
	    i__1 = *m * *k;
#line 322 "MB02JD.f"
	    i__2 = *ldwork - *l;
#line 322 "MB02JD.f"
	    dorgqr_(&i__1, l, l, &q[q_offset], ldq, &dwork[1], &dwork[*l + 1],
		     &i__2, &ierr);
/* Computing MAX */
#line 324 "MB02JD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[*l + 1] + *l;
#line 324 "MB02JD.f"
	    wrkopt = max(i__1,i__2);
#line 325 "MB02JD.f"
	    if (*n > 1) {
#line 326 "MB02JD.f"
		i__1 = *n - 1;
#line 326 "MB02JD.f"
		mb02kd_("Row", "Transpose", k, l, m, &i__1, l, &c_b10, &c_b11,
			 &tc[tc_offset], ldtc, &tr[tr_offset], ldtr, &q[
			q_offset], ldq, &r__[*l + 1 + r_dim1], ldr, &dwork[1],
			 ldwork, &ierr, (ftnlen)3, (ftnlen)9);
#line 329 "MB02JD.f"
	    }
/* Computing MAX */
#line 330 "MB02JD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 330 "MB02JD.f"
	    wrkopt = max(i__1,i__2);
#line 331 "MB02JD.f"
	} else {
#line 332 "MB02JD.f"
	    pdw = *m * *k * *l + 1;
#line 333 "MB02JD.f"
	    i__1 = *m * *k;
#line 333 "MB02JD.f"
	    i__2 = *m * *k;
#line 333 "MB02JD.f"
	    dlacpy_("All", &i__1, l, &tc[tc_offset], ldtc, &dwork[1], &i__2, (
		    ftnlen)3);
#line 334 "MB02JD.f"
	    i__1 = *m * *k;
#line 334 "MB02JD.f"
	    i__2 = *m * *k;
#line 334 "MB02JD.f"
	    i__3 = *ldwork - pdw - *l + 1;
#line 334 "MB02JD.f"
	    dgeqrf_(&i__1, l, &dwork[1], &i__2, &dwork[pdw], &dwork[pdw + *l],
		     &i__3, &ierr);
/* Computing MAX */
#line 336 "MB02JD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[pdw + *l] + pdw + *l - 1;
#line 336 "MB02JD.f"
	    wrkopt = max(i__1,i__2);
#line 337 "MB02JD.f"
	    i__1 = *m * *k;
#line 337 "MB02JD.f"
	    ma02ad_("Upper part", l, l, &dwork[1], &i__1, &r__[r_offset], ldr,
		     (ftnlen)10);
#line 338 "MB02JD.f"
	    i__1 = *m * *k;
#line 338 "MB02JD.f"
	    i__2 = *m * *k;
#line 338 "MB02JD.f"
	    i__3 = *ldwork - pdw - *l + 1;
#line 338 "MB02JD.f"
	    dorgqr_(&i__1, l, l, &dwork[1], &i__2, &dwork[pdw], &dwork[pdw + *
		    l], &i__3, &ierr);
/* Computing MAX */
#line 340 "MB02JD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[pdw + *l] + pdw + *l - 1;
#line 340 "MB02JD.f"
	    wrkopt = max(i__1,i__2);
#line 341 "MB02JD.f"
	    if (*n > 1) {
#line 342 "MB02JD.f"
		i__1 = *n - 1;
#line 342 "MB02JD.f"
		i__2 = *m * *k;
#line 342 "MB02JD.f"
		i__3 = *ldwork - pdw + 1;
#line 342 "MB02JD.f"
		mb02kd_("Row", "Transpose", k, l, m, &i__1, l, &c_b10, &c_b11,
			 &tc[tc_offset], ldtc, &tr[tr_offset], ldtr, &dwork[1]
			, &i__2, &r__[*l + 1 + r_dim1], ldr, &dwork[pdw], &
			i__3, &ierr, (ftnlen)3, (ftnlen)9);
#line 346 "MB02JD.f"
	    }
/* Computing MAX */
#line 347 "MB02JD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
#line 347 "MB02JD.f"
	    wrkopt = max(i__1,i__2);
#line 348 "MB02JD.f"
	}

/*        Quick return if N = 1. */

#line 352 "MB02JD.f"
	if (*n == 1) {
#line 353 "MB02JD.f"
	    dwork[1] = (doublereal) wrkopt;
#line 354 "MB02JD.f"
	    return 0;
#line 355 "MB02JD.f"
	}

/*        2nd column of the generator. */

#line 359 "MB02JD.f"
	pnr = (*n - 1) * *l * *k + 2;
#line 360 "MB02JD.f"
	i__1 = (*n - 1) * *l;
#line 360 "MB02JD.f"
	i__2 = (*n - 1) * *l;
#line 360 "MB02JD.f"
	ma02ad_("All", k, &i__1, &tr[tr_offset], ldtr, &dwork[2], &i__2, (
		ftnlen)3);

/*        3rd and 4th column of the generator. */

#line 364 "MB02JD.f"
	i__1 = (*n - 1) * *l;
#line 364 "MB02JD.f"
	i__2 = (*n - 1) * *l;
#line 364 "MB02JD.f"
	dlacpy_("All", &i__1, l, &r__[*l + 1 + r_dim1], ldr, &dwork[pnr], &
		i__2, (ftnlen)3);
#line 366 "MB02JD.f"
	pt = (*m - 1) * *k + 1;
#line 367 "MB02JD.f"
	pdw = pnr + (*n - 1) * *l * *l;

/* Computing MIN */
#line 369 "MB02JD.f"
	i__2 = *m, i__3 = *n - 1;
#line 369 "MB02JD.f"
	i__1 = min(i__2,i__3);
#line 369 "MB02JD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 370 "MB02JD.f"
	    i__2 = (*n - 1) * *l;
#line 370 "MB02JD.f"
	    ma02ad_("All", k, l, &tc[pt + tc_dim1], ldtc, &dwork[pdw], &i__2, 
		    (ftnlen)3);
#line 372 "MB02JD.f"
	    pt -= *k;
#line 373 "MB02JD.f"
	    pdw += *l;
#line 374 "MB02JD.f"
/* L10: */
#line 374 "MB02JD.f"
	}

#line 376 "MB02JD.f"
	pt = 1;

#line 378 "MB02JD.f"
	i__1 = *n - 1;
#line 378 "MB02JD.f"
	for (i__ = *m + 1; i__ <= i__1; ++i__) {
#line 379 "MB02JD.f"
	    i__2 = (*n - 1) * *l;
#line 379 "MB02JD.f"
	    ma02ad_("All", k, l, &tr[pt * tr_dim1 + 1], ldtr, &dwork[pdw], &
		    i__2, (ftnlen)3);
#line 381 "MB02JD.f"
	    pt += *l;
#line 382 "MB02JD.f"
	    pdw += *l;
#line 383 "MB02JD.f"
/* L20: */
#line 383 "MB02JD.f"
	}

#line 385 "MB02JD.f"
	if (compq) {
#line 386 "MB02JD.f"
	    pdq = ((*k << 1) + *l) * (*n - 1) * *l + 2;
#line 387 "MB02JD.f"
	    pdw = ((*k << 1) + *l) * ((*n - 1) * *l + *m * *k) + 2;
#line 388 "MB02JD.f"
	    pnq = pdq + *m * *k * *k;
#line 389 "MB02JD.f"
	    i__1 = *m * *k;
#line 389 "MB02JD.f"
	    dlaset_("All", k, k, &c_b11, &c_b10, &dwork[pdq], &i__1, (ftnlen)
		    3);
#line 390 "MB02JD.f"
	    i__1 = (*m - 1) * *k;
#line 390 "MB02JD.f"
	    i__2 = *m * *k;
#line 390 "MB02JD.f"
	    dlaset_("All", &i__1, k, &c_b11, &c_b11, &dwork[pdq + *k], &i__2, 
		    (ftnlen)3);
#line 392 "MB02JD.f"
	    i__1 = *m * *k;
#line 392 "MB02JD.f"
	    i__2 = *m * *k;
#line 392 "MB02JD.f"
	    dlacpy_("All", &i__1, l, &q[q_offset], ldq, &dwork[pnq], &i__2, (
		    ftnlen)3);
#line 393 "MB02JD.f"
	    i__1 = *m * *k;
#line 393 "MB02JD.f"
	    i__2 = *m * *k;
#line 393 "MB02JD.f"
	    dlaset_("All", &i__1, k, &c_b11, &c_b11, &dwork[pnq + *m * *l * *
		    k], &i__2, (ftnlen)3);
#line 395 "MB02JD.f"
	} else {
#line 396 "MB02JD.f"
	    pdw = ((*k << 1) + *l) * (*n - 1) * *l + 2;
#line 397 "MB02JD.f"
	}
#line 398 "MB02JD.f"
	pre = 1;
#line 399 "MB02JD.f"
	stps = *s - 1;
#line 400 "MB02JD.f"
    } else {

/*        Set workspace pointers. */

#line 404 "MB02JD.f"
	pnr = (*n - 1) * *l * *k + 2;
#line 405 "MB02JD.f"
	if (compq) {
#line 406 "MB02JD.f"
	    pdq = ((*k << 1) + *l) * (*n - 1) * *l + 2;
#line 407 "MB02JD.f"
	    pdw = ((*k << 1) + *l) * ((*n - 1) * *l + *m * *k) + 2;
#line 408 "MB02JD.f"
	    pnq = pdq + *m * *k * *k;
#line 409 "MB02JD.f"
	} else {
#line 410 "MB02JD.f"
	    pdw = ((*k << 1) + *l) * (*n - 1) * *l + 2;
#line 411 "MB02JD.f"
	}
#line 412 "MB02JD.f"
	pre = *p;
#line 413 "MB02JD.f"
	stps = *s;
#line 414 "MB02JD.f"
    }

/*     Determine suitable size for the block Housholder reflectors. */

#line 418 "MB02JD.f"
    if (compq) {
/* Computing MAX */
#line 419 "MB02JD.f"
	i__1 = *l + *m * *k, i__2 = (*n - pre + 1) * *l;
#line 419 "MB02JD.f"
	len = max(i__1,i__2);
#line 420 "MB02JD.f"
    } else {
#line 421 "MB02JD.f"
	len = (*n - pre + 1) * *l;
#line 422 "MB02JD.f"
    }
/* Computing MIN */
#line 423 "MB02JD.f"
    i__1 = ilaenv_(&c__1, "DGELQF", " ", &len, l, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 423 "MB02JD.f"
    nb = min(i__1,*l);
#line 424 "MB02JD.f"
    kk = pdw + *l * 6 - 1;
/* Computing MAX */
#line 425 "MB02JD.f"
    i__1 = wrkopt, i__2 = kk + len * nb;
#line 425 "MB02JD.f"
    wrkopt = max(i__1,i__2);
#line 426 "MB02JD.f"
    kk = *ldwork - kk;
#line 427 "MB02JD.f"
    if (kk < len * nb) {
#line 427 "MB02JD.f"
	nb = kk / len;
#line 427 "MB02JD.f"
    }
/* Computing MAX */
#line 428 "MB02JD.f"
    i__1 = 2, i__2 = ilaenv_(&c__2, "DGELQF", " ", &len, l, &c_n1, &c_n1, (
	    ftnlen)6, (ftnlen)1);
#line 428 "MB02JD.f"
    nbmin = max(i__1,i__2);
#line 429 "MB02JD.f"
    if (nb < nbmin) {
#line 429 "MB02JD.f"
	nb = 0;
#line 429 "MB02JD.f"
    }
#line 430 "MB02JD.f"
    colr = *l + 1;

/*     Generator reduction process. */

#line 434 "MB02JD.f"
    len = (*n - pre) * *l;
#line 435 "MB02JD.f"
    shfr = (pre - 1) * *l;
#line 436 "MB02JD.f"
    i__1 = pre + stps - 1;
#line 436 "MB02JD.f"
    for (i__ = pre; i__ <= i__1; ++i__) {

/*        IF M*K < N*L the last block might have less than L columns. */

/* Computing MIN */
#line 440 "MB02JD.f"
	i__2 = *l, i__3 = *m * *k - i__ * *l;
#line 440 "MB02JD.f"
	kk = min(i__2,i__3);
#line 441 "MB02JD.f"
	dlacpy_("Lower", &len, &kk, &r__[colr - *l + (colr - *l) * r_dim1], 
		ldr, &r__[colr + colr * r_dim1], ldr, (ftnlen)5);
#line 443 "MB02JD.f"
	i__2 = kk + *k;
#line 443 "MB02JD.f"
	i__3 = *l + *k;
#line 443 "MB02JD.f"
	i__4 = (*n - 1) * *l;
#line 443 "MB02JD.f"
	i__5 = (*n - 1) * *l;
#line 443 "MB02JD.f"
	i__6 = *ldwork - pdw - *l * 6 + 1;
#line 443 "MB02JD.f"
	mb02cu_("Column", &kk, &i__2, &i__3, &nb, &r__[colr + colr * r_dim1], 
		ldr, &dwork[shfr + 2], &i__4, &dwork[pnr + shfr], &i__5, &rnk,
		 ipvt, &dwork[pdw], &c_b11, &dwork[pdw + *l * 6], &i__6, &
		ierr, (ftnlen)6);
#line 447 "MB02JD.f"
	if (ierr != 0) {

/*           Error return:  The rank condition is (numerically) not */
/*                          satisfied. */

#line 452 "MB02JD.f"
	    *info = 1;
#line 453 "MB02JD.f"
	    return 0;
#line 454 "MB02JD.f"
	}
#line 455 "MB02JD.f"
	if (len > kk) {
#line 456 "MB02JD.f"
	    i__2 = len - kk;
#line 456 "MB02JD.f"
	    i__3 = kk + *k;
#line 456 "MB02JD.f"
	    i__4 = *l + *k;
#line 456 "MB02JD.f"
	    i__5 = (*n - 1) * *l;
#line 456 "MB02JD.f"
	    i__6 = (*n - 1) * *l;
#line 456 "MB02JD.f"
	    i__7 = (*n - 1) * *l;
#line 456 "MB02JD.f"
	    i__8 = (*n - 1) * *l;
#line 456 "MB02JD.f"
	    i__9 = *ldwork - pdw - *l * 6 + 1;
#line 456 "MB02JD.f"
	    mb02cv_("Column", "NoStructure", &kk, &i__2, &i__3, &i__4, &nb, &
		    c_n1, &r__[colr + colr * r_dim1], ldr, &dwork[shfr + 2], &
		    i__5, &dwork[pnr + shfr], &i__6, &r__[colr + kk + colr * 
		    r_dim1], ldr, &dwork[shfr + kk + 2], &i__7, &dwork[pnr + 
		    shfr + kk], &i__8, &dwork[pdw], &dwork[pdw + *l * 6], &
		    i__9, &ierr, (ftnlen)6, (ftnlen)11);
#line 463 "MB02JD.f"
	}
#line 464 "MB02JD.f"
	if (compq) {
#line 465 "MB02JD.f"
	    dlaset_("All", k, &kk, &c_b11, &c_b11, &q[colr * q_dim1 + 1], ldq,
		     (ftnlen)3);
#line 466 "MB02JD.f"
	    if (*m > 1) {
#line 467 "MB02JD.f"
		i__2 = (*m - 1) * *k;
#line 467 "MB02JD.f"
		dlacpy_("All", &i__2, &kk, &q[(colr - *l) * q_dim1 + 1], ldq, 
			&q[*k + 1 + colr * q_dim1], ldq, (ftnlen)3);
#line 469 "MB02JD.f"
	    }
#line 470 "MB02JD.f"
	    i__2 = *m * *k;
#line 470 "MB02JD.f"
	    i__3 = kk + *k;
#line 470 "MB02JD.f"
	    i__4 = *l + *k;
#line 470 "MB02JD.f"
	    i__5 = (*n - 1) * *l;
#line 470 "MB02JD.f"
	    i__6 = (*n - 1) * *l;
#line 470 "MB02JD.f"
	    i__7 = *m * *k;
#line 470 "MB02JD.f"
	    i__8 = *m * *k;
#line 470 "MB02JD.f"
	    i__9 = *ldwork - pdw - *l * 6 + 1;
#line 470 "MB02JD.f"
	    mb02cv_("Column", "NoStructure", &kk, &i__2, &i__3, &i__4, &nb, &
		    c_n1, &r__[colr + colr * r_dim1], ldr, &dwork[shfr + 2], &
		    i__5, &dwork[pnr + shfr], &i__6, &q[colr * q_dim1 + 1], 
		    ldq, &dwork[pdq], &i__7, &dwork[pnq], &i__8, &dwork[pdw], 
		    &dwork[pdw + *l * 6], &i__9, &ierr, (ftnlen)6, (ftnlen)11)
		    ;
#line 476 "MB02JD.f"
	}
#line 477 "MB02JD.f"
	len -= *l;
#line 478 "MB02JD.f"
	colr += *l;
#line 479 "MB02JD.f"
	shfr += *l;
#line 480 "MB02JD.f"
/* L30: */
#line 480 "MB02JD.f"
    }

#line 482 "MB02JD.f"
    dwork[1] = (doublereal) wrkopt;
#line 483 "MB02JD.f"
    return 0;

/* *** Last line of MB02JD *** */
} /* mb02jd_ */

