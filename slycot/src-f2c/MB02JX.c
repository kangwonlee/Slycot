#line 1 "MB02JX.f"
/* MB02JX.f -- translated by f2c (version 20100827).
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

#line 1 "MB02JX.f"
/* Table of constant values */

static doublereal c_b10 = 1.;
static doublereal c_b11 = 0.;
static integer c__1 = 1;
static integer c_n1 = -1;

/* Subroutine */ int mb02jx_(char *job, integer *k, integer *l, integer *m, 
	integer *n, doublereal *tc, integer *ldtc, doublereal *tr, integer *
	ldtr, integer *rnk, doublereal *q, integer *ldq, doublereal *r__, 
	integer *ldr, integer *jpvt, doublereal *tol1, doublereal *tol2, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen job_len)
{
    /* System generated locals */
    integer q_dim1, q_offset, r_dim1, r_offset, tc_dim1, tc_offset, tr_dim1, 
	    tr_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, jj, kk, mk, pp, pt, gap, len, pdp, pdq, nzc, pdw;
    static doublereal nrm;
    static integer pnq, pnr, ppr, rdef, rrdf, ierr;
    static logical last;
    static doublereal temp;
    static integer rrnk;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal ltol1, ltol2;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    mb02kd_(char *, char *, integer *, integer *, integer *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    dscal_(integer *, doublereal *, doublereal *, integer *), mb02cu_(
	    char *, integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, ftnlen), mb02cv_(char *, char *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static integer cpcol;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical compq;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer jwork;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen), dgeqp3_(integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen), dorgqr_(integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *);
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

/*     To compute a low rank QR factorization with column pivoting of a */
/*     K*M-by-L*N block Toeplitz matrix T with blocks of size (K,L); */
/*     specifically, */
/*                                     T */
/*                           T P =  Q R , */

/*     where R is lower trapezoidal, P is a block permutation matrix */
/*     and Q^T Q = I. The number of columns in R is equivalent to the */
/*     numerical rank of T with respect to the given tolerance TOL1. */
/*     Note that the pivoting scheme is local, i.e., only columns */
/*     belonging to the same block in T are permuted. */

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

/*     TC      (input) DOUBLE PRECISION array, dimension (LDTC, L) */
/*             The leading M*K-by-L part of this array must contain */
/*             the first block column of T. */

/*     LDTC    INTEGER */
/*             The leading dimension of the array TC. */
/*             LDTC >= MAX(1,M*K). */

/*     TR      (input)  DOUBLE PRECISION array, dimension (LDTR,(N-1)*L) */
/*             The leading K-by-(N-1)*L part of this array must contain */
/*             the first block row of T without the leading K-by-L */
/*             block. */

/*     LDTR    INTEGER */
/*             The leading dimension of the array TR.  LDTR >= MAX(1,K). */

/*     RNK     (output)  INTEGER */
/*             The number of columns in R, which is equivalent to the */
/*             numerical rank of T. */

/*     Q       (output)  DOUBLE PRECISION array, dimension (LDQ,RNK) */
/*             If JOB = 'Q', then the leading M*K-by-RNK part of this */
/*             array contains the factor Q. */
/*             If JOB = 'R', then this array is not referenced. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q. */
/*             LDQ >= MAX(1,M*K),  if JOB = 'Q'; */
/*             LDQ >= 1,           if JOB = 'R'. */

/*     R       (output)  DOUBLE PRECISION array, dimension (LDR,RNK) */
/*             The leading N*L-by-RNK part of this array contains the */
/*             lower trapezoidal factor R. */

/*     LDR     INTEGER */
/*             The leading dimension of the array R. */
/*             LDR >= MAX(1,N*L) */

/*     JPVT    (output)  INTEGER array, dimension (MIN(M*K,N*L)) */
/*             This array records the column pivoting performed. */
/*             If JPVT(j) = k, then the j-th column of T*P was */
/*             the k-th column of T. */

/*     Tolerances */

/*     TOL1    DOUBLE PRECISION */
/*             If TOL1 >= 0.0, the user supplied diagonal tolerance; */
/*             if TOL1 < 0.0, a default diagonal tolerance is used. */

/*     TOL2    DOUBLE PRECISION */
/*             If TOL2 >= 0.0, the user supplied offdiagonal tolerance; */
/*             if TOL2 < 0.0, a default offdiagonal tolerance is used. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK;  DWORK(2) and DWORK(3) return the used values */
/*             for TOL1 and TOL2, respectively. */
/*             On exit, if INFO = -19,  DWORK(1) returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 3, ( M*K + ( N - 1 )*L )*( L + 2*K ) + 9*L */
/*                                 + MAX(M*K,(N-1)*L) ),    if JOB = 'Q'; */
/*             LDWORK >= MAX( 3, ( N - 1 )*L*( L + 2*K + 1 ) + 9*L, */
/*                                 M*K*( L + 1 ) + L ),     if JOB = 'R'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  due to perturbations induced by roundoff errors, or */
/*                   removal of nearly linearly dependent columns of the */
/*                   generator, the Schur algorithm encountered a */
/*                   situation where a diagonal element in the negative */
/*                   generator is larger in magnitude than the */
/*                   corresponding diagonal element in the positive */
/*                   generator (modulo TOL1); */
/*             = 2:  due to perturbations induced by roundoff errors, or */
/*                   removal of nearly linearly dependent columns of the */
/*                   generator, the Schur algorithm encountered a */
/*                   situation where diagonal elements in the positive */
/*                   and negative generator are equal in magnitude */
/*                   (modulo TOL1), but the offdiagonal elements suggest */
/*                   that these columns are not linearly dependent */
/*                   (modulo TOL2*ABS(diagonal element)). */

/*     METHOD */

/*     Householder transformations and modified hyperbolic rotations */
/*     are used in the Schur algorithm [1], [2]. */
/*     If, during the process, the hyperbolic norm of a row in the */
/*     leading part of the generator is found to be less than or equal */
/*     to TOL1, then this row is not reduced. If the difference of the */
/*     corresponding columns has a norm less than or equal to TOL2 times */
/*     the magnitude of the leading element, then this column is removed */
/*     from the generator, as well as from R. Otherwise, the algorithm */
/*     breaks down. TOL1 is set to norm(TC)*sqrt(eps) and TOL2 is set */
/*     to N*L*sqrt(eps) by default. */
/*     If M*K > L, the columns of T are permuted so that the diagonal */
/*     elements in one block column of R have decreasing magnitudes. */

/*     REFERENCES */

/*     [1] Kailath, T. and Sayed, A. */
/*         Fast Reliable Algorithms for Matrices with Structure. */
/*         SIAM Publications, Philadelphia, 1999. */

/*     [2] Kressner, D. and Van Dooren, P. */
/*         Factorizations and linear system solvers for matrices with */
/*         Toeplitz structure. */
/*         SLICOT Working Note 2000-2, 2000. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0(K*RNK*L*M*N) floating point operations. */

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
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

#line 231 "MB02JX.f"
    /* Parameter adjustments */
#line 231 "MB02JX.f"
    tc_dim1 = *ldtc;
#line 231 "MB02JX.f"
    tc_offset = 1 + tc_dim1;
#line 231 "MB02JX.f"
    tc -= tc_offset;
#line 231 "MB02JX.f"
    tr_dim1 = *ldtr;
#line 231 "MB02JX.f"
    tr_offset = 1 + tr_dim1;
#line 231 "MB02JX.f"
    tr -= tr_offset;
#line 231 "MB02JX.f"
    q_dim1 = *ldq;
#line 231 "MB02JX.f"
    q_offset = 1 + q_dim1;
#line 231 "MB02JX.f"
    q -= q_offset;
#line 231 "MB02JX.f"
    r_dim1 = *ldr;
#line 231 "MB02JX.f"
    r_offset = 1 + r_dim1;
#line 231 "MB02JX.f"
    r__ -= r_offset;
#line 231 "MB02JX.f"
    --jpvt;
#line 231 "MB02JX.f"
    --dwork;
#line 231 "MB02JX.f"

#line 231 "MB02JX.f"
    /* Function Body */
#line 231 "MB02JX.f"
    *info = 0;
#line 232 "MB02JX.f"
    wrkopt = 3;
#line 233 "MB02JX.f"
    mk = *m * *k;
#line 234 "MB02JX.f"
    compq = lsame_(job, "Q", (ftnlen)1, (ftnlen)1);
#line 235 "MB02JX.f"
    if (compq) {
/* Computing MAX */
/* Computing MAX */
#line 236 "MB02JX.f"
	i__3 = mk, i__4 = (*n - 1) * *l;
#line 236 "MB02JX.f"
	i__1 = 3, i__2 = (mk + (*n - 1) * *l) * (*l + (*k << 1)) + *l * 9 + 
		max(i__3,i__4);
#line 236 "MB02JX.f"
	wrkmin = max(i__1,i__2);
#line 238 "MB02JX.f"
    } else {
/* Computing MAX */
/* Computing MAX */
#line 239 "MB02JX.f"
	i__3 = (*n - 1) * *l * (*l + (*k << 1) + 1) + *l * 9, i__4 = mk * (*l 
		+ 1) + *l;
#line 239 "MB02JX.f"
	i__1 = 3, i__2 = max(i__3,i__4);
#line 239 "MB02JX.f"
	wrkmin = max(i__1,i__2);
#line 241 "MB02JX.f"
    }

/*     Check the scalar input parameters. */

#line 245 "MB02JX.f"
    if (! (compq || lsame_(job, "R", (ftnlen)1, (ftnlen)1))) {
#line 246 "MB02JX.f"
	*info = -1;
#line 247 "MB02JX.f"
    } else if (*k < 0) {
#line 248 "MB02JX.f"
	*info = -2;
#line 249 "MB02JX.f"
    } else if (*l < 0) {
#line 250 "MB02JX.f"
	*info = -3;
#line 251 "MB02JX.f"
    } else if (*m < 0) {
#line 252 "MB02JX.f"
	*info = -4;
#line 253 "MB02JX.f"
    } else if (*n < 0) {
#line 254 "MB02JX.f"
	*info = -5;
#line 255 "MB02JX.f"
    } else if (*ldtc < max(1,mk)) {
#line 256 "MB02JX.f"
	*info = -7;
#line 257 "MB02JX.f"
    } else if (*ldtr < max(1,*k)) {
#line 258 "MB02JX.f"
	*info = -9;
#line 259 "MB02JX.f"
    } else if (*ldq < 1 || compq && *ldq < mk) {
#line 260 "MB02JX.f"
	*info = -12;
#line 261 "MB02JX.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 261 "MB02JX.f"
	i__1 = 1, i__2 = *n * *l;
#line 261 "MB02JX.f"
	if (*ldr < max(i__1,i__2)) {
#line 262 "MB02JX.f"
	    *info = -14;
#line 263 "MB02JX.f"
	} else if (*ldwork < wrkmin) {
#line 264 "MB02JX.f"
	    dwork[1] = (doublereal) wrkmin;
#line 265 "MB02JX.f"
	    *info = -19;
#line 266 "MB02JX.f"
	}
#line 266 "MB02JX.f"
    }

/*     Return if there were illegal values. */

#line 270 "MB02JX.f"
    if (*info != 0) {
#line 271 "MB02JX.f"
	i__1 = -(*info);
#line 271 "MB02JX.f"
	xerbla_("MB02JX", &i__1, (ftnlen)6);
#line 272 "MB02JX.f"
	return 0;
#line 273 "MB02JX.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 277 "MB02JX.f"
    i__1 = min(*m,*n), i__1 = min(i__1,*k);
#line 277 "MB02JX.f"
    if (min(i__1,*l) == 0) {
#line 278 "MB02JX.f"
	*rnk = 0;
#line 279 "MB02JX.f"
	dwork[1] = (doublereal) wrkopt;
#line 280 "MB02JX.f"
	dwork[2] = 0.;
#line 281 "MB02JX.f"
	dwork[3] = 0.;
#line 282 "MB02JX.f"
	return 0;
#line 283 "MB02JX.f"
    }

#line 285 "MB02JX.f"
    wrkopt = wrkmin;

#line 287 "MB02JX.f"
    if (mk <= *l) {

/*        Catch M*K <= L. */

#line 291 "MB02JX.f"
	dlacpy_("All", &mk, l, &tc[tc_offset], ldtc, &dwork[1], &mk, (ftnlen)
		3);
#line 292 "MB02JX.f"
	pdw = mk * *l + 1;
#line 293 "MB02JX.f"
	jwork = pdw + mk;
#line 294 "MB02JX.f"
	i__1 = *ldwork - jwork + 1;
#line 294 "MB02JX.f"
	dgeqrf_(&mk, l, &dwork[1], &mk, &dwork[pdw], &dwork[jwork], &i__1, &
		ierr);
/* Computing MAX */
#line 296 "MB02JX.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 296 "MB02JX.f"
	wrkopt = max(i__1,i__2);
#line 297 "MB02JX.f"
	ma02ad_("Upper part", &mk, l, &dwork[1], &mk, &r__[r_offset], ldr, (
		ftnlen)10);
#line 298 "MB02JX.f"
	i__1 = *ldwork - jwork + 1;
#line 298 "MB02JX.f"
	dorgqr_(&mk, &mk, &mk, &dwork[1], &mk, &dwork[pdw], &dwork[jwork], &
		i__1, &ierr);
/* Computing MAX */
#line 300 "MB02JX.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 300 "MB02JX.f"
	wrkopt = max(i__1,i__2);
#line 301 "MB02JX.f"
	if (compq) {
#line 301 "MB02JX.f"
	    dlacpy_("All", &mk, &mk, &dwork[1], &mk, &q[q_offset], ldq, (
		    ftnlen)3);
#line 301 "MB02JX.f"
	}
#line 303 "MB02JX.f"
	pdw = mk * mk + 1;
#line 304 "MB02JX.f"
	if (*n > 1) {
#line 305 "MB02JX.f"
	    i__1 = *n - 1;
#line 305 "MB02JX.f"
	    i__2 = *ldwork - pdw + 1;
#line 305 "MB02JX.f"
	    mb02kd_("Row", "Transpose", k, l, m, &i__1, &mk, &c_b10, &c_b11, &
		    tc[tc_offset], ldtc, &tr[tr_offset], ldtr, &dwork[1], &mk,
		     &r__[*l + 1 + r_dim1], ldr, &dwork[pdw], &i__2, &ierr, (
		    ftnlen)3, (ftnlen)9);
#line 308 "MB02JX.f"
	}
/* Computing MAX */
#line 309 "MB02JX.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
#line 309 "MB02JX.f"
	wrkopt = max(i__1,i__2);

#line 311 "MB02JX.f"
	i__1 = mk;
#line 311 "MB02JX.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 312 "MB02JX.f"
	    jpvt[i__] = i__;
#line 313 "MB02JX.f"
/* L10: */
#line 313 "MB02JX.f"
	}

#line 315 "MB02JX.f"
	*rnk = mk;
#line 316 "MB02JX.f"
	dwork[1] = (doublereal) wrkopt;
#line 317 "MB02JX.f"
	dwork[2] = 0.;
#line 318 "MB02JX.f"
	dwork[3] = 0.;
#line 319 "MB02JX.f"
	return 0;
#line 320 "MB02JX.f"
    }

/*     Compute the generator: */

/*     1st column of the generator. */

#line 326 "MB02JX.f"
    i__1 = *l;
#line 326 "MB02JX.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 327 "MB02JX.f"
	jpvt[i__] = 0;
#line 328 "MB02JX.f"
/* L20: */
#line 328 "MB02JX.f"
    }

#line 330 "MB02JX.f"
    ltol1 = *tol1;
#line 331 "MB02JX.f"
    ltol2 = *tol2;

#line 333 "MB02JX.f"
    if (compq) {
#line 334 "MB02JX.f"
	dlacpy_("All", &mk, l, &tc[tc_offset], ldtc, &q[q_offset], ldq, (
		ftnlen)3);
#line 335 "MB02JX.f"
	i__1 = *ldwork - *l;
#line 335 "MB02JX.f"
	dgeqp3_(&mk, l, &q[q_offset], ldq, &jpvt[1], &dwork[1], &dwork[*l + 1]
		, &i__1, &ierr);
/* Computing MAX */
#line 337 "MB02JX.f"
	i__1 = wrkopt, i__2 = (integer) dwork[*l + 1] + *l;
#line 337 "MB02JX.f"
	wrkopt = max(i__1,i__2);

#line 339 "MB02JX.f"
	if (ltol1 < 0.) {

/*           Compute default tolerance LTOL1. */

/*           Estimate the 2-norm of the first block column of the */
/*           matrix with 5 power iterations. */

#line 346 "MB02JX.f"
	    temp = 1. / sqrt((doublereal) (*l));
#line 347 "MB02JX.f"
	    dlaset_("All", l, &c__1, &temp, &temp, &dwork[*l + 1], &c__1, (
		    ftnlen)3);

#line 349 "MB02JX.f"
	    for (i__ = 1; i__ <= 5; ++i__) {
#line 350 "MB02JX.f"
		dtrmv_("Upper", "NonTranspose", "NonUnit", l, &q[q_offset], 
			ldq, &dwork[*l + 1], &c__1, (ftnlen)5, (ftnlen)12, (
			ftnlen)7);
#line 352 "MB02JX.f"
		dtrmv_("Upper", "Transpose", "NonUnit", l, &q[q_offset], ldq, 
			&dwork[*l + 1], &c__1, (ftnlen)5, (ftnlen)9, (ftnlen)
			7);
#line 354 "MB02JX.f"
		nrm = dnrm2_(l, &dwork[*l + 1], &c__1);
#line 355 "MB02JX.f"
		d__1 = 1. / nrm;
#line 355 "MB02JX.f"
		dscal_(l, &d__1, &dwork[*l + 1], &c__1);
#line 356 "MB02JX.f"
/* L30: */
#line 356 "MB02JX.f"
	    }

#line 358 "MB02JX.f"
	    ltol1 = sqrt(nrm * dlamch_("Epsilon", (ftnlen)7));
#line 359 "MB02JX.f"
	}

#line 361 "MB02JX.f"
	i__ = *l;

#line 363 "MB02JX.f"
L40:
#line 364 "MB02JX.f"
	if ((d__1 = q[i__ + i__ * q_dim1], abs(d__1)) <= ltol1) {
#line 365 "MB02JX.f"
	    --i__;
#line 366 "MB02JX.f"
	    if (i__ > 0) {
#line 366 "MB02JX.f"
		goto L40;
#line 366 "MB02JX.f"
	    }
#line 367 "MB02JX.f"
	}

#line 369 "MB02JX.f"
	rrnk = i__;
#line 370 "MB02JX.f"
	rrdf = *l - rrnk;
#line 371 "MB02JX.f"
	ma02ad_("Upper", &rrnk, l, &q[q_offset], ldq, &r__[r_offset], ldr, (
		ftnlen)5);
#line 372 "MB02JX.f"
	if (rrnk > 1) {
#line 372 "MB02JX.f"
	    i__1 = *l - 1;
#line 372 "MB02JX.f"
	    i__2 = rrnk - 1;
#line 372 "MB02JX.f"
	    dlaset_("Upper", &i__1, &i__2, &c_b11, &c_b11, &r__[(r_dim1 << 1) 
		    + 1], ldr, (ftnlen)5);
#line 372 "MB02JX.f"
	}
#line 374 "MB02JX.f"
	i__1 = *ldwork - *l;
#line 374 "MB02JX.f"
	dorgqr_(&mk, l, &rrnk, &q[q_offset], ldq, &dwork[1], &dwork[*l + 1], &
		i__1, &ierr);
/* Computing MAX */
#line 376 "MB02JX.f"
	i__1 = wrkopt, i__2 = (integer) dwork[*l + 1] + *l;
#line 376 "MB02JX.f"
	wrkopt = max(i__1,i__2);
#line 377 "MB02JX.f"
	if (*n > 1) {
#line 378 "MB02JX.f"
	    i__1 = *n - 1;
#line 378 "MB02JX.f"
	    mb02kd_("Row", "Transpose", k, l, m, &i__1, &rrnk, &c_b10, &c_b11,
		     &tc[tc_offset], ldtc, &tr[tr_offset], ldtr, &q[q_offset],
		     ldq, &r__[*l + 1 + r_dim1], ldr, &dwork[1], ldwork, &
		    ierr, (ftnlen)3, (ftnlen)9);
/* Computing MAX */
#line 381 "MB02JX.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 381 "MB02JX.f"
	    wrkopt = max(i__1,i__2);
#line 382 "MB02JX.f"
	}

#line 384 "MB02JX.f"
    } else {

#line 386 "MB02JX.f"
	pdw = mk * *l + 1;
#line 387 "MB02JX.f"
	jwork = pdw + *l;
#line 388 "MB02JX.f"
	dlacpy_("All", &mk, l, &tc[tc_offset], ldtc, &dwork[1], &mk, (ftnlen)
		3);
#line 389 "MB02JX.f"
	i__1 = *ldwork - jwork + 1;
#line 389 "MB02JX.f"
	dgeqp3_(&mk, l, &dwork[1], &mk, &jpvt[1], &dwork[pdw], &dwork[jwork], 
		&i__1, &ierr);
/* Computing MAX */
#line 391 "MB02JX.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 391 "MB02JX.f"
	wrkopt = max(i__1,i__2);

#line 393 "MB02JX.f"
	if (ltol1 < 0.) {

/*           Compute default tolerance LTOL1. */

/*           Estimate the 2-norm of the first block column of the */
/*           matrix with 5 power iterations. */

#line 400 "MB02JX.f"
	    temp = 1. / sqrt((doublereal) (*l));
#line 401 "MB02JX.f"
	    dlaset_("All", l, &c__1, &temp, &temp, &dwork[jwork], &c__1, (
		    ftnlen)3);

#line 403 "MB02JX.f"
	    for (i__ = 1; i__ <= 5; ++i__) {
#line 404 "MB02JX.f"
		dtrmv_("Upper", "NonTranspose", "NonUnit", l, &dwork[1], &mk, 
			&dwork[jwork], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)
			7);
#line 406 "MB02JX.f"
		dtrmv_("Upper", "Transpose", "NonUnit", l, &dwork[1], &mk, &
			dwork[jwork], &c__1, (ftnlen)5, (ftnlen)9, (ftnlen)7);
#line 408 "MB02JX.f"
		nrm = dnrm2_(l, &dwork[jwork], &c__1);
#line 409 "MB02JX.f"
		d__1 = 1. / nrm;
#line 409 "MB02JX.f"
		dscal_(l, &d__1, &dwork[jwork], &c__1);
#line 410 "MB02JX.f"
/* L50: */
#line 410 "MB02JX.f"
	    }

#line 412 "MB02JX.f"
	    ltol1 = sqrt(nrm * dlamch_("Epsilon", (ftnlen)7));
#line 413 "MB02JX.f"
	}

#line 415 "MB02JX.f"
	rrnk = *l;
#line 416 "MB02JX.f"
	i__ = (*l - 1) * mk + *l;

#line 418 "MB02JX.f"
L60:
#line 419 "MB02JX.f"
	if ((d__1 = dwork[i__], abs(d__1)) <= ltol1) {
#line 420 "MB02JX.f"
	    --rrnk;
#line 421 "MB02JX.f"
	    i__ = i__ - mk - 1;
#line 422 "MB02JX.f"
	    if (i__ > 0) {
#line 422 "MB02JX.f"
		goto L60;
#line 422 "MB02JX.f"
	    }
#line 423 "MB02JX.f"
	}

#line 425 "MB02JX.f"
	rrdf = *l - rrnk;
#line 426 "MB02JX.f"
	ma02ad_("Upper part", &rrnk, l, &dwork[1], &mk, &r__[r_offset], ldr, (
		ftnlen)10);
#line 427 "MB02JX.f"
	if (rrnk > 1) {
#line 427 "MB02JX.f"
	    i__1 = *l - 1;
#line 427 "MB02JX.f"
	    i__2 = rrnk - 1;
#line 427 "MB02JX.f"
	    dlaset_("Upper", &i__1, &i__2, &c_b11, &c_b11, &r__[(r_dim1 << 1) 
		    + 1], ldr, (ftnlen)5);
#line 427 "MB02JX.f"
	}
#line 429 "MB02JX.f"
	i__1 = *ldwork - jwork + 1;
#line 429 "MB02JX.f"
	dorgqr_(&mk, l, &rrnk, &dwork[1], &mk, &dwork[pdw], &dwork[jwork], &
		i__1, &ierr);
/* Computing MAX */
#line 431 "MB02JX.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 431 "MB02JX.f"
	wrkopt = max(i__1,i__2);
#line 432 "MB02JX.f"
	if (*n > 1) {
#line 433 "MB02JX.f"
	    i__1 = *n - 1;
#line 433 "MB02JX.f"
	    i__2 = *ldwork - pdw + 1;
#line 433 "MB02JX.f"
	    mb02kd_("Row", "Transpose", k, l, m, &i__1, &rrnk, &c_b10, &c_b11,
		     &tc[tc_offset], ldtc, &tr[tr_offset], ldtr, &dwork[1], &
		    mk, &r__[*l + 1 + r_dim1], ldr, &dwork[pdw], &i__2, &ierr,
		     (ftnlen)3, (ftnlen)9);
/* Computing MAX */
#line 436 "MB02JX.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
#line 436 "MB02JX.f"
	    wrkopt = max(i__1,i__2);
#line 437 "MB02JX.f"
	}
#line 438 "MB02JX.f"
    }

/*     Quick return if N = 1. */

#line 442 "MB02JX.f"
    if (*n == 1) {
#line 443 "MB02JX.f"
	*rnk = rrnk;
#line 444 "MB02JX.f"
	dwork[1] = (doublereal) wrkopt;
#line 445 "MB02JX.f"
	dwork[2] = ltol1;
#line 446 "MB02JX.f"
	dwork[3] = 0.;
#line 447 "MB02JX.f"
	return 0;
#line 448 "MB02JX.f"
    }

/*     Compute default tolerance LTOL2. */

#line 452 "MB02JX.f"
    if (ltol2 < 0.) {
#line 452 "MB02JX.f"
	ltol2 = (doublereal) (*n * *l) * sqrt(dlamch_("Epsilon", (ftnlen)7));
#line 452 "MB02JX.f"
    }

#line 455 "MB02JX.f"
    i__1 = *l;
#line 455 "MB02JX.f"
    for (j = 1; j <= i__1; ++j) {
#line 456 "MB02JX.f"
	dcopy_(&rrnk, &r__[j + r_dim1], ldr, &r__[*l + jpvt[j] + (rrnk + 1) * 
		r_dim1], ldr);
#line 457 "MB02JX.f"
/* L70: */
#line 457 "MB02JX.f"
    }

#line 459 "MB02JX.f"
    if (*n > 2) {
#line 459 "MB02JX.f"
	i__1 = (*n - 2) * *l;
#line 459 "MB02JX.f"
	dlacpy_("All", &i__1, &rrnk, &r__[*l + 1 + r_dim1], ldr, &r__[(*l << 
		1) + 1 + (rrnk + 1) * r_dim1], ldr, (ftnlen)3);
#line 459 "MB02JX.f"
    }

/*     2nd column of the generator. */

#line 465 "MB02JX.f"
    if (rrdf > 0) {
#line 465 "MB02JX.f"
	i__1 = min(rrdf,*k);
#line 465 "MB02JX.f"
	i__2 = (*n - 1) * *l;
#line 465 "MB02JX.f"
	ma02ad_("All", &i__1, &i__2, &tr[tr_offset], ldtr, &r__[*l + 1 + ((
		rrnk << 1) + 1) * r_dim1], ldr, (ftnlen)3);
#line 465 "MB02JX.f"
    }
#line 468 "MB02JX.f"
    if (*k > rrdf) {
#line 468 "MB02JX.f"
	i__1 = *k - rrdf;
#line 468 "MB02JX.f"
	i__2 = (*n - 1) * *l;
#line 468 "MB02JX.f"
	i__3 = (*n - 1) * *l;
#line 468 "MB02JX.f"
	ma02ad_("All", &i__1, &i__2, &tr[rrdf + 1 + tr_dim1], ldtr, &dwork[1],
		 &i__3, (ftnlen)3);
#line 468 "MB02JX.f"
    }

/*     3rd column of the generator. */

/* Computing MAX */
#line 474 "MB02JX.f"
    i__1 = 0, i__2 = *k - rrdf;
#line 474 "MB02JX.f"
    pnr = (*n - 1) * *l * max(i__1,i__2) + 1;
#line 475 "MB02JX.f"
    i__1 = (*n - 1) * *l;
#line 475 "MB02JX.f"
    i__2 = (*n - 1) * *l;
#line 475 "MB02JX.f"
    dlacpy_("All", &i__1, &rrnk, &r__[*l + 1 + r_dim1], ldr, &dwork[pnr], &
	    i__2, (ftnlen)3);

/*     4th column of the generator. */

#line 480 "MB02JX.f"
    pdw = pnr + (*n - 1) * *l * rrnk;
#line 481 "MB02JX.f"
    pt = (*m - 1) * *k + 1;

/* Computing MIN */
#line 483 "MB02JX.f"
    i__2 = *m, i__3 = *n - 1;
#line 483 "MB02JX.f"
    i__1 = min(i__2,i__3);
#line 483 "MB02JX.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 484 "MB02JX.f"
	i__2 = (*n - 1) * *l;
#line 484 "MB02JX.f"
	ma02ad_("All", k, l, &tc[pt + tc_dim1], ldtc, &dwork[pdw], &i__2, (
		ftnlen)3);
#line 485 "MB02JX.f"
	pt -= *k;
#line 486 "MB02JX.f"
	pdw += *l;
#line 487 "MB02JX.f"
/* L80: */
#line 487 "MB02JX.f"
    }

#line 489 "MB02JX.f"
    pt = 1;

#line 491 "MB02JX.f"
    i__1 = *n - 1;
#line 491 "MB02JX.f"
    for (i__ = *m + 1; i__ <= i__1; ++i__) {
#line 492 "MB02JX.f"
	i__2 = (*n - 1) * *l;
#line 492 "MB02JX.f"
	ma02ad_("All", k, l, &tr[pt * tr_dim1 + 1], ldtr, &dwork[pdw], &i__2, 
		(ftnlen)3);
#line 493 "MB02JX.f"
	pt += *l;
#line 494 "MB02JX.f"
	pdw += *l;
#line 495 "MB02JX.f"
/* L90: */
#line 495 "MB02JX.f"
    }

#line 497 "MB02JX.f"
    if (compq) {
#line 498 "MB02JX.f"
	pdq = pnr + (*n - 1) * *l * (rrnk + *k);
/* Computing MAX */
#line 499 "MB02JX.f"
	i__1 = 0, i__2 = *k - rrdf;
#line 499 "MB02JX.f"
	pnq = pdq + mk * max(i__1,i__2);
#line 500 "MB02JX.f"
	pdw = pnq + mk * (rrnk + *k);
#line 501 "MB02JX.f"
	dlacpy_("All", &mk, &rrnk, &q[q_offset], ldq, &dwork[pnq], &mk, (
		ftnlen)3);
#line 502 "MB02JX.f"
	if (*m > 1) {
#line 502 "MB02JX.f"
	    i__1 = (*m - 1) * *k;
#line 502 "MB02JX.f"
	    dlacpy_("All", &i__1, &rrnk, &q[q_offset], ldq, &q[*k + 1 + (rrnk 
		    + 1) * q_dim1], ldq, (ftnlen)3);
#line 502 "MB02JX.f"
	}
#line 505 "MB02JX.f"
	dlaset_("All", k, &rrnk, &c_b11, &c_b11, &q[(rrnk + 1) * q_dim1 + 1], 
		ldq, (ftnlen)3);
#line 506 "MB02JX.f"
	if (rrdf > 0) {
#line 506 "MB02JX.f"
	    dlaset_("All", &mk, &rrdf, &c_b11, &c_b10, &q[((rrnk << 1) + 1) * 
		    q_dim1 + 1], ldq, (ftnlen)3);
#line 506 "MB02JX.f"
	}
/* Computing MAX */
#line 509 "MB02JX.f"
	i__2 = 0, i__3 = *k - rrdf;
#line 509 "MB02JX.f"
	i__1 = max(i__2,i__3);
#line 509 "MB02JX.f"
	dlaset_("All", &rrdf, &i__1, &c_b11, &c_b11, &dwork[pdq], &mk, (
		ftnlen)3);
#line 511 "MB02JX.f"
	i__1 = *m * *k - rrdf;
/* Computing MAX */
#line 511 "MB02JX.f"
	i__3 = 0, i__4 = *k - rrdf;
#line 511 "MB02JX.f"
	i__2 = max(i__3,i__4);
#line 511 "MB02JX.f"
	dlaset_("All", &i__1, &i__2, &c_b11, &c_b10, &dwork[pdq + rrdf], &mk, 
		(ftnlen)3);
#line 513 "MB02JX.f"
	dlaset_("All", &mk, k, &c_b11, &c_b11, &dwork[pnq + mk * rrnk], &mk, (
		ftnlen)3);
#line 514 "MB02JX.f"
    } else {
#line 515 "MB02JX.f"
	pdw = pnr + (*n - 1) * *l * (rrnk + *k);
#line 516 "MB02JX.f"
    }
#line 517 "MB02JX.f"
    ppr = 1;
#line 518 "MB02JX.f"
    *rnk = rrnk;
#line 519 "MB02JX.f"
    rdef = rrdf;
#line 520 "MB02JX.f"
    len = *n * *l;
/* Computing MIN */
#line 521 "MB02JX.f"
    i__1 = *n * *l;
#line 521 "MB02JX.f"
    gap = *n * *l - min(i__1,mk);

/*     KK is the number of columns in the leading part of the */
/*     generator. After sufficiently many rank drops or if */
/*     M*K < N*L it may be less than L. */

/* Computing MIN */
#line 527 "MB02JX.f"
    i__1 = *l + *k - rdef;
#line 527 "MB02JX.f"
    kk = min(i__1,*l);
/* Computing MIN */
#line 528 "MB02JX.f"
    i__1 = kk, i__2 = mk - *l;
#line 528 "MB02JX.f"
    kk = min(i__1,i__2);

/*     Generator reduction process. */

/* Computing MIN */
#line 532 "MB02JX.f"
    i__2 = mk, i__3 = *n * *l;
#line 532 "MB02JX.f"
    i__1 = min(i__2,i__3);
#line 532 "MB02JX.f"
    i__4 = *l;
#line 532 "MB02JX.f"
    for (i__ = *l + 1; i__4 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__4) {
/* Computing MIN */
#line 533 "MB02JX.f"
	i__2 = mk, i__3 = *n * *l;
#line 533 "MB02JX.f"
	if (i__ + *l <= min(i__2,i__3)) {
#line 534 "MB02JX.f"
	    last = FALSE_;
#line 535 "MB02JX.f"
	} else {
#line 536 "MB02JX.f"
	    last = TRUE_;
#line 537 "MB02JX.f"
	}
/* Computing MAX */
#line 538 "MB02JX.f"
	i__2 = *k - rdef;
#line 538 "MB02JX.f"
	pp = kk + max(i__2,0);
#line 539 "MB02JX.f"
	len -= *l;
#line 540 "MB02JX.f"
	i__2 = *l + *k - rdef;
#line 540 "MB02JX.f"
	i__3 = (*n - 1) * *l;
#line 540 "MB02JX.f"
	i__5 = (*n - 1) * *l;
#line 540 "MB02JX.f"
	i__6 = *ldwork - pdw - *l * 5 + 1;
#line 540 "MB02JX.f"
	mb02cu_("Deficient", &kk, &pp, &i__2, &c_n1, &r__[i__ + (*rnk + 1) * 
		r_dim1], ldr, &dwork[ppr], &i__3, &dwork[pnr], &i__5, &rrnk, &
		jpvt[i__], &dwork[pdw], &ltol1, &dwork[pdw + *l * 5], &i__6, &
		ierr, (ftnlen)9);
#line 544 "MB02JX.f"
	if (ierr != 0) {

/*           Error return:  The current generator is indefinite. */

#line 548 "MB02JX.f"
	    *info = 1;
#line 549 "MB02JX.f"
	    return 0;
#line 550 "MB02JX.f"
	}

/*        Apply pivoting to other columns of R. */

#line 554 "MB02JX.f"
	pdp = pdw + *l * 6 - i__;

#line 556 "MB02JX.f"
	i__2 = i__ + kk - 1;
#line 556 "MB02JX.f"
	for (j = i__; j <= i__2; ++j) {
#line 557 "MB02JX.f"
	    jpvt[j] = jpvt[j] + i__ - 1;
#line 558 "MB02JX.f"
	    dwork[pdp + jpvt[j]] = (doublereal) j;
#line 559 "MB02JX.f"
/* L100: */
#line 559 "MB02JX.f"
	}

#line 561 "MB02JX.f"
	i__2 = i__ + kk - 1;
#line 561 "MB02JX.f"
	for (j = i__; j <= i__2; ++j) {
#line 562 "MB02JX.f"
	    temp = (doublereal) j;
#line 563 "MB02JX.f"
	    jj = j - 1;

#line 565 "MB02JX.f"
L110:
#line 566 "MB02JX.f"
	    ++jj;
#line 567 "MB02JX.f"
	    if (dwork[pdp + jj] != temp) {
#line 567 "MB02JX.f"
		goto L110;
#line 567 "MB02JX.f"
	    }

#line 569 "MB02JX.f"
	    if (jj != j) {
#line 570 "MB02JX.f"
		dwork[pdp + jj] = dwork[pdp + j];
#line 571 "MB02JX.f"
		dswap_(rnk, &r__[j + r_dim1], ldr, &r__[jj + r_dim1], ldr);
#line 572 "MB02JX.f"
	    }
#line 573 "MB02JX.f"
/* L120: */
#line 573 "MB02JX.f"
	}

#line 575 "MB02JX.f"
	i__2 = i__ + *l - 1;
#line 575 "MB02JX.f"
	for (j = i__ + kk; j <= i__2; ++j) {
#line 576 "MB02JX.f"
	    jpvt[j] = j;
#line 577 "MB02JX.f"
/* L130: */
#line 577 "MB02JX.f"
	}

/*        Apply reduction to other rows of R. */

#line 581 "MB02JX.f"
	if (len > kk) {
#line 582 "MB02JX.f"
	    i__2 = len - kk;
#line 582 "MB02JX.f"
	    i__3 = *l + *k - rdef;
#line 582 "MB02JX.f"
	    i__5 = (*n - 1) * *l;
#line 582 "MB02JX.f"
	    i__6 = (*n - 1) * *l;
#line 582 "MB02JX.f"
	    i__7 = (*n - 1) * *l;
#line 582 "MB02JX.f"
	    i__8 = (*n - 1) * *l;
#line 582 "MB02JX.f"
	    i__9 = *ldwork - pdw - *l * 5 + 1;
#line 582 "MB02JX.f"
	    mb02cv_("Deficient", "NoStructure", &kk, &i__2, &pp, &i__3, &c_n1,
		     &rrnk, &r__[i__ + (*rnk + 1) * r_dim1], ldr, &dwork[ppr],
		     &i__5, &dwork[pnr], &i__6, &r__[i__ + kk + (*rnk + 1) * 
		    r_dim1], ldr, &dwork[ppr + kk], &i__7, &dwork[pnr + kk], &
		    i__8, &dwork[pdw], &dwork[pdw + *l * 5], &i__9, &ierr, (
		    ftnlen)9, (ftnlen)11);
#line 588 "MB02JX.f"
	}

/*        Apply reduction to Q. */

#line 592 "MB02JX.f"
	if (compq) {
#line 593 "MB02JX.f"
	    i__2 = *l + *k - rdef;
#line 593 "MB02JX.f"
	    i__3 = (*n - 1) * *l;
#line 593 "MB02JX.f"
	    i__5 = (*n - 1) * *l;
#line 593 "MB02JX.f"
	    i__6 = *ldwork - pdw - *l * 5 + 1;
#line 593 "MB02JX.f"
	    mb02cv_("Deficient", "NoStructure", &kk, &mk, &pp, &i__2, &c_n1, &
		    rrnk, &r__[i__ + (*rnk + 1) * r_dim1], ldr, &dwork[ppr], &
		    i__3, &dwork[pnr], &i__5, &q[(*rnk + 1) * q_dim1 + 1], 
		    ldq, &dwork[pdq], &mk, &dwork[pnq], &mk, &dwork[pdw], &
		    dwork[pdw + *l * 5], &i__6, &ierr, (ftnlen)9, (ftnlen)11);
#line 599 "MB02JX.f"
	}

/*        Inspection of the rank deficient columns: */
/*        Look for small diagonal entries. */

#line 604 "MB02JX.f"
	nzc = 0;

#line 606 "MB02JX.f"
	i__2 = rrnk + 1;
#line 606 "MB02JX.f"
	for (j = kk; j >= i__2; --j) {
#line 607 "MB02JX.f"
	    if ((d__1 = r__[i__ + j - 1 + (*rnk + j) * r_dim1], abs(d__1)) <= 
		    ltol1) {
#line 607 "MB02JX.f"
		++nzc;
#line 607 "MB02JX.f"
	    }
#line 608 "MB02JX.f"
/* L140: */
#line 608 "MB02JX.f"
	}

/*        The last NZC columns of the generator cannot be removed. */
/*        Now, decide whether for the other rank deficient columns */
/*        it is safe to remove. */

#line 614 "MB02JX.f"
	pt = pnr;

#line 616 "MB02JX.f"
	i__2 = kk - nzc;
#line 616 "MB02JX.f"
	for (j = rrnk + 1; j <= i__2; ++j) {
#line 617 "MB02JX.f"
	    temp = r__[i__ + j - 1 + (*rnk + j) * r_dim1];
#line 618 "MB02JX.f"
	    i__3 = len - j - gap;
#line 618 "MB02JX.f"
	    dscal_(&i__3, &temp, &r__[i__ + j + (*rnk + j) * r_dim1], &c__1);
#line 619 "MB02JX.f"
	    i__3 = len - j - gap;
#line 619 "MB02JX.f"
	    d__1 = -dwork[pt + j - 1];
#line 619 "MB02JX.f"
	    daxpy_(&i__3, &d__1, &dwork[pt + j], &c__1, &r__[i__ + j + (*rnk 
		    + j) * r_dim1], &c__1);
#line 621 "MB02JX.f"
	    i__3 = len - j - gap;
#line 621 "MB02JX.f"
	    if (dnrm2_(&i__3, &r__[i__ + j + (*rnk + j) * r_dim1], &c__1) > 
		    ltol2 * abs(temp)) {

/*              Unlucky case: */
/*              It is neither advisable to remove the whole column nor */
/*              possible to remove the diagonal entries by Hyperbolic */
/*              rotations. */

#line 629 "MB02JX.f"
		*info = 2;
#line 630 "MB02JX.f"
		return 0;
#line 631 "MB02JX.f"
	    }
#line 632 "MB02JX.f"
	    pt += (*n - 1) * *l;
#line 633 "MB02JX.f"
/* L150: */
#line 633 "MB02JX.f"
	}

/*        Annihilate unwanted elements in the factor R. */

#line 637 "MB02JX.f"
	rrdf = kk - rrnk;
#line 638 "MB02JX.f"
	i__2 = i__ - 1;
#line 638 "MB02JX.f"
	dlaset_("All", &i__2, &rrnk, &c_b11, &c_b11, &r__[(*rnk + 1) * r_dim1 
		+ 1], ldr, (ftnlen)3);
#line 639 "MB02JX.f"
	i__2 = *l - 1;
#line 639 "MB02JX.f"
	i__3 = rrnk - 1;
#line 639 "MB02JX.f"
	dlaset_("Upper", &i__2, &i__3, &c_b11, &c_b11, &r__[i__ + (*rnk + 2) *
		 r_dim1], ldr, (ftnlen)5);

/*        Construct the generator for the next step. */

#line 644 "MB02JX.f"
	if (! last) {

/*           Compute KK for the next step. */

/* Computing MIN */
#line 648 "MB02JX.f"
	    i__2 = *l + *k - rdef - rrdf + nzc;
#line 648 "MB02JX.f"
	    kk = min(i__2,*l);
/* Computing MIN */
#line 649 "MB02JX.f"
	    i__2 = kk, i__3 = mk - i__ - *l + 1;
#line 649 "MB02JX.f"
	    kk = min(i__2,i__3);

#line 651 "MB02JX.f"
	    if (kk <= 0) {
#line 652 "MB02JX.f"
		*rnk += rrnk;
#line 653 "MB02JX.f"
		goto L200;
#line 654 "MB02JX.f"
	    }

#line 656 "MB02JX.f"
	    dlaset_("All", l, &rrdf, &c_b11, &c_b11, &r__[i__ + (*rnk + rrnk 
		    + 1) * r_dim1], ldr, (ftnlen)3);

/*           The columns with small diagonal entries form parts of the */
/*           new positive generator. */

#line 662 "MB02JX.f"
	    if (rrdf - nzc > 0 && nzc > 0) {
#line 663 "MB02JX.f"
		cpcol = min(nzc,kk);

#line 665 "MB02JX.f"
		i__2 = *rnk + rrnk + cpcol;
#line 665 "MB02JX.f"
		for (j = *rnk + rrnk + 1; j <= i__2; ++j) {
#line 666 "MB02JX.f"
		    i__3 = len - *l;
#line 666 "MB02JX.f"
		    dcopy_(&i__3, &r__[i__ + *l + (j + rrdf - nzc) * r_dim1], 
			    &c__1, &r__[i__ + *l + j * r_dim1], &c__1);
#line 668 "MB02JX.f"
/* L160: */
#line 668 "MB02JX.f"
		}

#line 670 "MB02JX.f"
	    }

/*           Construct the leading parts of the positive generator. */

/* Computing MIN */
#line 674 "MB02JX.f"
	    i__2 = rrnk, i__3 = kk - nzc;
#line 674 "MB02JX.f"
	    cpcol = min(i__2,i__3);
#line 675 "MB02JX.f"
	    if (cpcol > 0) {

#line 677 "MB02JX.f"
		i__2 = i__ + *l - 1;
#line 677 "MB02JX.f"
		for (j = i__; j <= i__2; ++j) {
#line 678 "MB02JX.f"
		    dcopy_(&cpcol, &r__[j + (*rnk + 1) * r_dim1], ldr, &r__[
			    jpvt[j] + *l + (*rnk + rrnk + nzc + 1) * r_dim1], 
			    ldr);
#line 680 "MB02JX.f"
/* L170: */
#line 680 "MB02JX.f"
		}

#line 682 "MB02JX.f"
		if (len > *l << 1) {
#line 683 "MB02JX.f"
		    i__2 = len - (*l << 1);
#line 683 "MB02JX.f"
		    dlacpy_("All", &i__2, &cpcol, &r__[i__ + *l + (*rnk + 1) *
			     r_dim1], ldr, &r__[i__ + (*l << 1) + (*rnk + 
			    rrnk + nzc + 1) * r_dim1], ldr, (ftnlen)3);
#line 685 "MB02JX.f"
		}
#line 686 "MB02JX.f"
	    }
#line 687 "MB02JX.f"
	    ppr += *l;

/*           Refill the leading parts of the positive generator. */

/* Computing MIN */
#line 691 "MB02JX.f"
	    i__2 = *k - rdef, i__3 = kk - rrnk - nzc;
#line 691 "MB02JX.f"
	    cpcol = min(i__2,i__3);
#line 692 "MB02JX.f"
	    if (cpcol > 0) {
#line 693 "MB02JX.f"
		i__2 = len - *l;
#line 693 "MB02JX.f"
		i__3 = (*n - 1) * *l;
#line 693 "MB02JX.f"
		dlacpy_("All", &i__2, &cpcol, &dwork[ppr], &i__3, &r__[i__ + *
			l + (*rnk + (rrnk << 1) + nzc + 1) * r_dim1], ldr, (
			ftnlen)3);
#line 695 "MB02JX.f"
		ppr += cpcol * (*n - 1) * *l;
#line 696 "MB02JX.f"
	    }
#line 697 "MB02JX.f"
	    pnr = pnr + (rrdf - nzc) * (*n - 1) * *l + *l;

/*           Do the same things for Q. */

#line 701 "MB02JX.f"
	    if (compq) {
#line 702 "MB02JX.f"
		if (rrdf - nzc > 0 && nzc > 0) {
#line 703 "MB02JX.f"
		    cpcol = min(nzc,kk);

#line 705 "MB02JX.f"
		    i__2 = *rnk + rrnk + cpcol;
#line 705 "MB02JX.f"
		    for (j = *rnk + rrnk + 1; j <= i__2; ++j) {
#line 706 "MB02JX.f"
			dcopy_(&mk, &q[(j + rrdf - nzc) * q_dim1 + 1], &c__1, 
				&q[j * q_dim1 + 1], &c__1);
#line 707 "MB02JX.f"
/* L180: */
#line 707 "MB02JX.f"
		    }

#line 709 "MB02JX.f"
		}
/* Computing MIN */
#line 710 "MB02JX.f"
		i__2 = rrnk, i__3 = kk - nzc;
#line 710 "MB02JX.f"
		cpcol = min(i__2,i__3);
#line 711 "MB02JX.f"
		if (cpcol > 0) {
#line 712 "MB02JX.f"
		    dlaset_("All", k, &cpcol, &c_b11, &c_b11, &q[(*rnk + rrnk 
			    + nzc + 1) * q_dim1 + 1], ldq, (ftnlen)3);
#line 714 "MB02JX.f"
		    if (*m > 1) {
#line 714 "MB02JX.f"
			i__2 = (*m - 1) * *k;
#line 714 "MB02JX.f"
			dlacpy_("All", &i__2, &cpcol, &q[(*rnk + 1) * q_dim1 
				+ 1], ldq, &q[*k + 1 + (*rnk + rrnk + nzc + 1)
				 * q_dim1], ldq, (ftnlen)3);
#line 714 "MB02JX.f"
		    }
#line 717 "MB02JX.f"
		}
/* Computing MIN */
#line 718 "MB02JX.f"
		i__2 = *k - rdef, i__3 = kk - rrnk - nzc;
#line 718 "MB02JX.f"
		cpcol = min(i__2,i__3);
#line 719 "MB02JX.f"
		if (cpcol > 0) {
#line 720 "MB02JX.f"
		    dlacpy_("All", &mk, &cpcol, &dwork[pdq], &mk, &q[(*rnk + (
			    rrnk << 1) + nzc + 1) * q_dim1 + 1], ldq, (ftnlen)
			    3);
#line 722 "MB02JX.f"
		    pdq += cpcol * mk;
#line 723 "MB02JX.f"
		}
#line 724 "MB02JX.f"
		pnq += (rrdf - nzc) * mk;
#line 725 "MB02JX.f"
	    }
#line 726 "MB02JX.f"
	}
#line 727 "MB02JX.f"
	*rnk += rrnk;
#line 728 "MB02JX.f"
	rdef = rdef + rrdf - nzc;
#line 729 "MB02JX.f"
/* L190: */
#line 729 "MB02JX.f"
    }

#line 731 "MB02JX.f"
L200:
#line 732 "MB02JX.f"
    dwork[1] = (doublereal) wrkopt;
#line 733 "MB02JX.f"
    dwork[2] = ltol1;
#line 734 "MB02JX.f"
    dwork[3] = ltol2;

/* *** Last line of MB02JX *** */
#line 737 "MB02JX.f"
    return 0;
} /* mb02jx_ */

