#line 1 "MB02CD.f"
/* MB02CD.f -- translated by f2c (version 20100827).
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

#line 1 "MB02CD.f"
/* Table of constant values */

static doublereal c_b17 = 1.;
static doublereal c_b19 = 0.;
static integer c__1 = 1;

/* Subroutine */ int mb02cd_(char *job, char *typet, integer *k, integer *n, 
	doublereal *t, integer *ldt, doublereal *g, integer *ldg, doublereal *
	r__, integer *ldr, doublereal *l, integer *ldl, doublereal *cs, 
	integer *lcs, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen job_len, ftnlen typet_len)
{
    /* System generated locals */
    integer g_dim1, g_offset, l_dim1, l_offset, r_dim1, r_offset, t_dim1, 
	    t_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, ierr;
    extern /* Subroutine */ int mb02cx_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen), mb02cy_(
	    char *, char *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical compg, compl, compr;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical isrow;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), dpotrf_(char *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static integer starti, maxwrk, startr, startt;


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

/*     To compute the Cholesky factor and the generator and/or the */
/*     Cholesky factor of the inverse of a symmetric positive definite */
/*     (s.p.d.) block Toeplitz matrix T, defined by either its first */
/*     block row, or its first block column, depending on the routine */
/*     parameter TYPET. Transformation information is stored. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the output of the routine, as follows: */
/*             = 'G':  only computes the generator G of the inverse; */
/*             = 'R':  computes the generator G of the inverse and the */
/*                     Cholesky factor R of T, i.e., if TYPET = 'R', */
/*                     then R'*R = T, and if TYPET = 'C', then R*R' = T; */
/*             = 'L':  computes the generator G and the Cholesky factor L */
/*                     of the inverse, i.e., if TYPET = 'R', then */
/*                     L'*L = inv(T), and if TYPET = 'C', then */
/*                     L*L' = inv(T); */
/*             = 'A':  computes the generator G, the Cholesky factor L */
/*                     of the inverse and the Cholesky factor R of T; */
/*             = 'O':  only computes the Cholesky factor R of T. */

/*     TYPET   CHARACTER*1 */
/*             Specifies the type of T, as follows: */
/*             = 'R':  T contains the first block row of an s.p.d. block */
/*                     Toeplitz matrix; if demanded, the Cholesky factors */
/*                     R and L are upper and lower triangular, */
/*                     respectively, and G contains the transposed */
/*                     generator of the inverse; */
/*             = 'C':  T contains the first block column of an s.p.d. */
/*                     block Toeplitz matrix; if demanded, the Cholesky */
/*                     factors R and L are lower and upper triangular, */
/*                     respectively, and G contains the generator of the */
/*                     inverse. This choice results in a column oriented */
/*                     algorithm which is usually faster. */
/*             Note:   in the sequel, the notation x / y means that */
/*                     x corresponds to TYPET = 'R' and y corresponds to */
/*                     TYPET = 'C'. */

/*     Input/Output Parameters */

/*     K       (input)  INTEGER */
/*             The number of rows / columns in T, which should be equal */
/*             to the blocksize.  K >= 0. */

/*     N       (input)  INTEGER */
/*             The number of blocks in T.  N >= 0. */

/*     T       (input/output)  DOUBLE PRECISION array, dimension */
/*             (LDT,N*K) / (LDT,K) */
/*             On entry, the leading K-by-N*K / N*K-by-K part of this */
/*             array must contain the first block row / column of an */
/*             s.p.d. block Toeplitz matrix. */
/*             On exit, if INFO = 0, then the leading K-by-N*K / N*K-by-K */
/*             part of this array contains, in the first K-by-K block, */
/*             the upper / lower Cholesky factor of T(1:K,1:K), and in */
/*             the remaining part, the Householder transformations */
/*             applied during the process. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T. */
/*             LDT >= MAX(1,K),    if TYPET = 'R'; */
/*             LDT >= MAX(1,N*K),  if TYPET = 'C'. */

/*     G       (output)  DOUBLE PRECISION array, dimension */
/*             (LDG,N*K) / (LDG,2*K) */
/*             If INFO = 0 and JOB = 'G', 'R', 'L', or 'A', the leading */
/*             2*K-by-N*K / N*K-by-2*K part of this array contains, in */
/*             the first K-by-K block of the second block row / column, */
/*             the lower right block of L (necessary for updating */
/*             factorizations in SLICOT Library routine MB02DD), and */
/*             in the remaining part, the generator of the inverse of T. */
/*             Actually, to obtain a generator one has to set */
/*                 G(K+1:2*K, 1:K) = 0,    if TYPET = 'R'; */
/*                 G(1:K, K+1:2*K) = 0,    if TYPET = 'C'. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G. */
/*             LDG >= MAX(1,2*K),  if TYPET = 'R' and */
/*                                    JOB = 'G', 'R', 'L', or 'A'; */
/*             LDG >= MAX(1,N*K),  if TYPET = 'C' and */
/*                                    JOB = 'G', 'R', 'L', or 'A'; */
/*             LDG >= 1,           if JOB = 'O'. */

/*     R       (output)  DOUBLE PRECISION array, dimension (LDR,N*K) */
/*             If INFO = 0 and JOB = 'R', 'A', or 'O', then the leading */
/*             N*K-by-N*K part of this array contains the upper / lower */
/*             Cholesky factor of T. */
/*             The elements in the strictly lower / upper triangular part */
/*             are not referenced. */

/*     LDR     INTEGER */
/*             The leading dimension of the array R. */
/*             LDR >= MAX(1,N*K),  if JOB = 'R', 'A', or 'O'; */
/*             LDR >= 1,           if JOB = 'G', or 'L'. */

/*     L       (output)  DOUBLE PRECISION array, dimension (LDL,N*K) */
/*             If INFO = 0 and JOB = 'L', or 'A', then the leading */
/*             N*K-by-N*K part of this array contains the lower / upper */
/*             Cholesky factor of the inverse of T. */
/*             The elements in the strictly upper / lower triangular part */
/*             are not referenced. */

/*     LDL     INTEGER */
/*             The leading dimension of the array L. */
/*             LDL >= MAX(1,N*K),  if JOB = 'L', or 'A'; */
/*             LDL >= 1,           if JOB = 'G', 'R', or 'O'. */

/*     CS      (output)  DOUBLE PRECISION array, dimension (LCS) */
/*             If INFO = 0, then the leading 3*(N-1)*K part of this */
/*             array contains information about the hyperbolic rotations */
/*             and Householder transformations applied during the */
/*             process. This information is needed for updating the */
/*             factorizations in SLICOT Library routine MB02DD. */

/*     LCS     INTEGER */
/*             The length of the array CS.  LCS >= 3*(N-1)*K. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -16,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,(N-1)*K). */
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
/*                               3 2 */
/*     The algorithm requires 0(K N ) floating point operations. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Chemnitz, Germany, June 2000. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, July 2000, */
/*     February 2004. */

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

#line 230 "MB02CD.f"
    /* Parameter adjustments */
#line 230 "MB02CD.f"
    t_dim1 = *ldt;
#line 230 "MB02CD.f"
    t_offset = 1 + t_dim1;
#line 230 "MB02CD.f"
    t -= t_offset;
#line 230 "MB02CD.f"
    g_dim1 = *ldg;
#line 230 "MB02CD.f"
    g_offset = 1 + g_dim1;
#line 230 "MB02CD.f"
    g -= g_offset;
#line 230 "MB02CD.f"
    r_dim1 = *ldr;
#line 230 "MB02CD.f"
    r_offset = 1 + r_dim1;
#line 230 "MB02CD.f"
    r__ -= r_offset;
#line 230 "MB02CD.f"
    l_dim1 = *ldl;
#line 230 "MB02CD.f"
    l_offset = 1 + l_dim1;
#line 230 "MB02CD.f"
    l -= l_offset;
#line 230 "MB02CD.f"
    --cs;
#line 230 "MB02CD.f"
    --dwork;
#line 230 "MB02CD.f"

#line 230 "MB02CD.f"
    /* Function Body */
#line 230 "MB02CD.f"
    *info = 0;
#line 231 "MB02CD.f"
    compl = lsame_(job, "L", (ftnlen)1, (ftnlen)1) || lsame_(job, "A", (
	    ftnlen)1, (ftnlen)1);
#line 232 "MB02CD.f"
    compg = lsame_(job, "G", (ftnlen)1, (ftnlen)1) || lsame_(job, "R", (
	    ftnlen)1, (ftnlen)1) || compl;
#line 233 "MB02CD.f"
    compr = lsame_(job, "R", (ftnlen)1, (ftnlen)1) || lsame_(job, "A", (
	    ftnlen)1, (ftnlen)1) || lsame_(job, "O", (ftnlen)1, (ftnlen)1);
#line 235 "MB02CD.f"
    isrow = lsame_(typet, "R", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

#line 239 "MB02CD.f"
    if (! (compg || compr)) {
#line 240 "MB02CD.f"
	*info = -1;
#line 241 "MB02CD.f"
    } else if (! (isrow || lsame_(typet, "C", (ftnlen)1, (ftnlen)1))) {
#line 242 "MB02CD.f"
	*info = -2;
#line 243 "MB02CD.f"
    } else if (*k < 0) {
#line 244 "MB02CD.f"
	*info = -3;
#line 245 "MB02CD.f"
    } else if (*n < 0) {
#line 246 "MB02CD.f"
	*info = -4;
#line 247 "MB02CD.f"
    } else if (*ldt < 1 || isrow && *ldt < *k || ! isrow && *ldt < *n * *k) {
#line 249 "MB02CD.f"
	*info = -6;
#line 250 "MB02CD.f"
    } else if (*ldg < 1 || compg && (isrow && *ldg < *k << 1 || ! isrow && *
	    ldg < *n * *k)) {
#line 253 "MB02CD.f"
	*info = -8;
#line 254 "MB02CD.f"
    } else if (*ldr < 1 || compr && *ldr < *n * *k) {
#line 255 "MB02CD.f"
	*info = -10;
#line 256 "MB02CD.f"
    } else if (*ldl < 1 || compl && *ldl < *n * *k) {
#line 257 "MB02CD.f"
	*info = -12;
#line 258 "MB02CD.f"
    } else if (*lcs < (*n - 1) * 3 * *k) {
#line 259 "MB02CD.f"
	*info = -14;
#line 260 "MB02CD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 260 "MB02CD.f"
	i__1 = 1, i__2 = (*n - 1) * *k;
#line 260 "MB02CD.f"
	if (*ldwork < max(i__1,i__2)) {
/* Computing MAX */
#line 261 "MB02CD.f"
	    i__1 = 1, i__2 = (*n - 1) * *k;
#line 261 "MB02CD.f"
	    dwork[1] = (doublereal) max(i__1,i__2);
#line 262 "MB02CD.f"
	    *info = -16;
#line 263 "MB02CD.f"
	}
#line 263 "MB02CD.f"
    }

/*     Return if there were illegal values. */

#line 267 "MB02CD.f"
    if (*info != 0) {
#line 268 "MB02CD.f"
	i__1 = -(*info);
#line 268 "MB02CD.f"
	xerbla_("MB02CD", &i__1, (ftnlen)6);
#line 269 "MB02CD.f"
	return 0;
#line 270 "MB02CD.f"
    }

/*     Quick return if possible. */

#line 274 "MB02CD.f"
    if (min(*k,*n) == 0) {
#line 275 "MB02CD.f"
	dwork[1] = 1.;
#line 276 "MB02CD.f"
	return 0;
#line 277 "MB02CD.f"
    }

#line 279 "MB02CD.f"
    maxwrk = 1;
#line 280 "MB02CD.f"
    if (isrow) {

/*        T is the first block row of a block Toeplitz matrix. */
/*        Bring T to proper form by triangularizing its first block. */

#line 285 "MB02CD.f"
	dpotrf_("Upper", k, &t[t_offset], ldt, &ierr, (ftnlen)5);
#line 286 "MB02CD.f"
	if (ierr != 0) {

/*           Error return:  The matrix is not positive definite. */

#line 290 "MB02CD.f"
	    *info = 1;
#line 291 "MB02CD.f"
	    return 0;
#line 292 "MB02CD.f"
	}

#line 294 "MB02CD.f"
	if (*n > 1) {
#line 294 "MB02CD.f"
	    i__1 = (*n - 1) * *k;
#line 294 "MB02CD.f"
	    dtrsm_("Left", "Upper", "Transpose", "NonUnit", k, &i__1, &c_b17, 
		    &t[t_offset], ldt, &t[(*k + 1) * t_dim1 + 1], ldt, (
		    ftnlen)4, (ftnlen)5, (ftnlen)9, (ftnlen)7);
#line 294 "MB02CD.f"
	}

/*        Initialize the output matrices. */

#line 300 "MB02CD.f"
	if (compg) {
#line 301 "MB02CD.f"
	    i__1 = *k << 1;
#line 301 "MB02CD.f"
	    i__2 = *n * *k;
#line 301 "MB02CD.f"
	    dlaset_("All", &i__1, &i__2, &c_b19, &c_b19, &g[g_offset], ldg, (
		    ftnlen)3);
#line 302 "MB02CD.f"
	    i__1 = *ldg + 1;
#line 302 "MB02CD.f"
	    dlaset_("All", &c__1, k, &c_b17, &c_b17, &g[*k + 1 + g_dim1], &
		    i__1, (ftnlen)3);
#line 303 "MB02CD.f"
	    dtrsm_("Left", "Upper", "Transpose", "NonUnit", k, k, &c_b17, &t[
		    t_offset], ldt, &g[*k + 1 + g_dim1], ldg, (ftnlen)4, (
		    ftnlen)5, (ftnlen)9, (ftnlen)7);
#line 305 "MB02CD.f"
	    if (*n > 1) {
#line 305 "MB02CD.f"
		i__1 = (*n - 1) * *k;
#line 305 "MB02CD.f"
		dlacpy_("Upper", k, &i__1, &t[t_offset], ldt, &g[*k + 1 + (*k 
			+ 1) * g_dim1], ldg, (ftnlen)5);
#line 305 "MB02CD.f"
	    }
#line 308 "MB02CD.f"
	    dlacpy_("Lower", k, k, &g[*k + 1 + g_dim1], ldg, &g[g_offset], 
		    ldg, (ftnlen)5);
#line 309 "MB02CD.f"
	}

#line 311 "MB02CD.f"
	if (compl) {
#line 312 "MB02CD.f"
	    dlacpy_("Lower", k, k, &g[*k + 1 + g_dim1], ldg, &l[l_offset], 
		    ldl, (ftnlen)5);
#line 313 "MB02CD.f"
	}

#line 315 "MB02CD.f"
	if (compr) {
#line 316 "MB02CD.f"
	    i__1 = *n * *k;
#line 316 "MB02CD.f"
	    dlacpy_("Upper", k, &i__1, &t[t_offset], ldt, &r__[r_offset], ldr,
		     (ftnlen)5);
#line 317 "MB02CD.f"
	}

/*        Processing the generator. */

#line 321 "MB02CD.f"
	if (compg) {

/*           Here we use G as working array for holding the generator. */
/*           T contains the second row of the generator. */
/*           G contains in its first block row the second row of the */
/*           inverse generator. */
/*           The second block row of G is partitioned as follows: */

/*           [ First block of the inverse generator, ... */
/*             First row of the generator, ... */
/*             The rest of the blocks of the inverse generator ] */

/*           The reason for the odd partitioning is that the first block */
/*           of the inverse generator will be thrown out at the end and */
/*           we want to avoid reordering. */

/*           (N-1)*K locations of DWORK are used by SLICOT Library */
/*           routine MB02CY. */

#line 340 "MB02CD.f"
	    i__1 = *n;
#line 340 "MB02CD.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 341 "MB02CD.f"
		startr = (i__ - 1) * *k + 1;
#line 342 "MB02CD.f"
		starti = (*n - i__ + 1) * *k + 1;
#line 343 "MB02CD.f"
		startt = (i__ - 2) * 3 * *k + 1;

/*              Transformations acting on the generator: */

#line 347 "MB02CD.f"
		i__2 = *k * 3;
#line 347 "MB02CD.f"
		mb02cx_("Row", k, k, k, &g[*k + 1 + (*k + 1) * g_dim1], ldg, &
			t[startr * t_dim1 + 1], ldt, &cs[startt], &i__2, &
			dwork[1], ldwork, &ierr, (ftnlen)3);

#line 351 "MB02CD.f"
		if (ierr != 0) {

/*                 Error return:  The matrix is not positive definite. */

#line 355 "MB02CD.f"
		    *info = 1;
#line 356 "MB02CD.f"
		    return 0;
#line 357 "MB02CD.f"
		}

/* Computing MAX */
#line 359 "MB02CD.f"
		i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 359 "MB02CD.f"
		maxwrk = max(i__2,i__3);
#line 360 "MB02CD.f"
		if (*n > i__) {
#line 361 "MB02CD.f"
		    i__2 = (*n - i__) * *k;
#line 361 "MB02CD.f"
		    i__3 = *k * 3;
#line 361 "MB02CD.f"
		    mb02cy_("Row", "NoStructure", k, k, &i__2, k, &g[*k + 1 + 
			    ((*k << 1) + 1) * g_dim1], ldg, &t[(startr + *k) *
			     t_dim1 + 1], ldt, &t[startr * t_dim1 + 1], ldt, &
			    cs[startt], &i__3, &dwork[1], ldwork, &ierr, (
			    ftnlen)3, (ftnlen)11);
/* Computing MAX */
#line 365 "MB02CD.f"
		    i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 365 "MB02CD.f"
		    maxwrk = max(i__2,i__3);
#line 366 "MB02CD.f"
		}

#line 368 "MB02CD.f"
		if (compr) {
#line 369 "MB02CD.f"
		    i__2 = (*n - i__ + 1) * *k;
#line 369 "MB02CD.f"
		    dlacpy_("Upper", k, &i__2, &g[*k + 1 + (*k + 1) * g_dim1],
			     ldg, &r__[startr + startr * r_dim1], ldr, (
			    ftnlen)5);
#line 371 "MB02CD.f"
		}

/*              Transformations acting on the inverse generator: */

#line 375 "MB02CD.f"
		dlaset_("All", k, k, &c_b19, &c_b19, &g[*k + 1 + starti * 
			g_dim1], ldg, (ftnlen)3);
#line 377 "MB02CD.f"
		i__2 = *k * 3;
#line 377 "MB02CD.f"
		mb02cy_("Row", "Triangular", k, k, k, k, &g[*k + 1 + g_dim1], 
			ldg, &g[startr * g_dim1 + 1], ldg, &t[startr * t_dim1 
			+ 1], ldt, &cs[startt], &i__2, &dwork[1], ldwork, &
			ierr, (ftnlen)3, (ftnlen)10);
/* Computing MAX */
#line 380 "MB02CD.f"
		i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 380 "MB02CD.f"
		maxwrk = max(i__2,i__3);

#line 382 "MB02CD.f"
		i__2 = (i__ - 1) * *k;
#line 382 "MB02CD.f"
		i__3 = *k * 3;
#line 382 "MB02CD.f"
		mb02cy_("Row", "NoStructure", k, k, &i__2, k, &g[*k + 1 + 
			starti * g_dim1], ldg, &g[g_offset], ldg, &t[startr * 
			t_dim1 + 1], ldt, &cs[startt], &i__3, &dwork[1], 
			ldwork, &ierr, (ftnlen)3, (ftnlen)11);
/* Computing MAX */
#line 385 "MB02CD.f"
		i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 385 "MB02CD.f"
		maxwrk = max(i__2,i__3);

#line 387 "MB02CD.f"
		if (compl) {
#line 388 "MB02CD.f"
		    i__2 = (i__ - 1) * *k;
#line 388 "MB02CD.f"
		    dlacpy_("All", k, &i__2, &g[*k + 1 + starti * g_dim1], 
			    ldg, &l[startr + l_dim1], ldl, (ftnlen)3);
#line 390 "MB02CD.f"
		    dlacpy_("Lower", k, k, &g[*k + 1 + g_dim1], ldg, &l[
			    startr + ((i__ - 1) * *k + 1) * l_dim1], ldl, (
			    ftnlen)5);
#line 392 "MB02CD.f"
		}
#line 393 "MB02CD.f"
/* L10: */
#line 393 "MB02CD.f"
	    }

#line 395 "MB02CD.f"
	} else {

/*           Here R is used as working array for holding the generator. */
/*           Again, T contains the second row of the generator. */
/*           The current row of R contains the first row of the */
/*           generator. */

#line 402 "MB02CD.f"
	    if (*n > 1) {
#line 402 "MB02CD.f"
		i__1 = (*n - 1) * *k;
#line 402 "MB02CD.f"
		dlacpy_("Upper", k, &i__1, &t[t_offset], ldt, &r__[*k + 1 + (*
			k + 1) * r_dim1], ldr, (ftnlen)5);
#line 402 "MB02CD.f"
	    }

#line 406 "MB02CD.f"
	    i__1 = *n;
#line 406 "MB02CD.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 407 "MB02CD.f"
		startr = (i__ - 1) * *k + 1;
#line 408 "MB02CD.f"
		startt = (i__ - 2) * 3 * *k + 1;
#line 409 "MB02CD.f"
		i__2 = *k * 3;
#line 409 "MB02CD.f"
		mb02cx_("Row", k, k, k, &r__[startr + startr * r_dim1], ldr, &
			t[startr * t_dim1 + 1], ldt, &cs[startt], &i__2, &
			dwork[1], ldwork, &ierr, (ftnlen)3);
#line 412 "MB02CD.f"
		if (ierr != 0) {

/*                 Error return:  The matrix is not positive definite. */

#line 416 "MB02CD.f"
		    *info = 1;
#line 417 "MB02CD.f"
		    return 0;
#line 418 "MB02CD.f"
		}

/* Computing MAX */
#line 420 "MB02CD.f"
		i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 420 "MB02CD.f"
		maxwrk = max(i__2,i__3);
#line 421 "MB02CD.f"
		if (*n > i__) {
#line 422 "MB02CD.f"
		    i__2 = (*n - i__) * *k;
#line 422 "MB02CD.f"
		    i__3 = *k * 3;
#line 422 "MB02CD.f"
		    mb02cy_("Row", "NoStructure", k, k, &i__2, k, &r__[startr 
			    + (startr + *k) * r_dim1], ldr, &t[(startr + *k) *
			     t_dim1 + 1], ldt, &t[startr * t_dim1 + 1], ldt, &
			    cs[startt], &i__3, &dwork[1], ldwork, &ierr, (
			    ftnlen)3, (ftnlen)11);
/* Computing MAX */
#line 426 "MB02CD.f"
		    i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 426 "MB02CD.f"
		    maxwrk = max(i__2,i__3);

#line 428 "MB02CD.f"
		    i__2 = (*n - i__) * *k;
#line 428 "MB02CD.f"
		    dlacpy_("Upper", k, &i__2, &r__[startr + startr * r_dim1],
			     ldr, &r__[startr + *k + (startr + *k) * r_dim1], 
			    ldr, (ftnlen)5);
#line 430 "MB02CD.f"
		}
#line 431 "MB02CD.f"
/* L20: */
#line 431 "MB02CD.f"
	    }

#line 433 "MB02CD.f"
	}

#line 435 "MB02CD.f"
    } else {

/*        T is the first block column of a block Toeplitz matrix. */
/*        Bring T to proper form by triangularizing its first block. */

#line 440 "MB02CD.f"
	dpotrf_("Lower", k, &t[t_offset], ldt, &ierr, (ftnlen)5);
#line 441 "MB02CD.f"
	if (ierr != 0) {

/*           Error return:  The matrix is not positive definite. */

#line 445 "MB02CD.f"
	    *info = 1;
#line 446 "MB02CD.f"
	    return 0;
#line 447 "MB02CD.f"
	}

#line 449 "MB02CD.f"
	if (*n > 1) {
#line 449 "MB02CD.f"
	    i__1 = (*n - 1) * *k;
#line 449 "MB02CD.f"
	    dtrsm_("Right", "Lower", "Transpose", "NonUnit", &i__1, k, &c_b17,
		     &t[t_offset], ldt, &t[*k + 1 + t_dim1], ldt, (ftnlen)5, (
		    ftnlen)5, (ftnlen)9, (ftnlen)7);
#line 449 "MB02CD.f"
	}

/*        Initialize the output matrices. */

#line 455 "MB02CD.f"
	if (compg) {
#line 456 "MB02CD.f"
	    i__1 = *n * *k;
#line 456 "MB02CD.f"
	    i__2 = *k << 1;
#line 456 "MB02CD.f"
	    dlaset_("All", &i__1, &i__2, &c_b19, &c_b19, &g[g_offset], ldg, (
		    ftnlen)3);
#line 457 "MB02CD.f"
	    i__1 = *ldg + 1;
#line 457 "MB02CD.f"
	    dlaset_("All", &c__1, k, &c_b17, &c_b17, &g[(*k + 1) * g_dim1 + 1]
		    , &i__1, (ftnlen)3);
#line 458 "MB02CD.f"
	    dtrsm_("Right", "Lower", "Transpose", "NonUnit", k, k, &c_b17, &t[
		    t_offset], ldt, &g[(*k + 1) * g_dim1 + 1], ldg, (ftnlen)5,
		     (ftnlen)5, (ftnlen)9, (ftnlen)7);
#line 460 "MB02CD.f"
	    if (*n > 1) {
#line 460 "MB02CD.f"
		i__1 = (*n - 1) * *k;
#line 460 "MB02CD.f"
		dlacpy_("Lower", &i__1, k, &t[t_offset], ldt, &g[*k + 1 + (*k 
			+ 1) * g_dim1], ldg, (ftnlen)5);
#line 460 "MB02CD.f"
	    }
#line 463 "MB02CD.f"
	    dlacpy_("Upper", k, k, &g[(*k + 1) * g_dim1 + 1], ldg, &g[
		    g_offset], ldg, (ftnlen)5);
#line 464 "MB02CD.f"
	}

#line 466 "MB02CD.f"
	if (compl) {
#line 467 "MB02CD.f"
	    dlacpy_("Upper", k, k, &g[(*k + 1) * g_dim1 + 1], ldg, &l[
		    l_offset], ldl, (ftnlen)5);
#line 468 "MB02CD.f"
	}

#line 470 "MB02CD.f"
	if (compr) {
#line 471 "MB02CD.f"
	    i__1 = *n * *k;
#line 471 "MB02CD.f"
	    dlacpy_("Lower", &i__1, k, &t[t_offset], ldt, &r__[r_offset], ldr,
		     (ftnlen)5);
#line 472 "MB02CD.f"
	}

/*        Processing the generator. */

#line 476 "MB02CD.f"
	if (compg) {

/*           Here we use G as working array for holding the generator. */
/*           T contains the second column of the generator. */
/*           G contains in its first block column the second column of */
/*           the inverse generator. */
/*           The second block column of G is partitioned as follows: */

/*           [ First block of the inverse generator; ... */
/*             First column of the generator; ... */
/*             The rest of the blocks of the inverse generator ] */

/*           The reason for the odd partitioning is that the first block */
/*           of the inverse generator will be thrown out at the end and */
/*           we want to avoid reordering. */

/*           (N-1)*K locations of DWORK are used by SLICOT Library */
/*           routine MB02CY. */

#line 495 "MB02CD.f"
	    i__1 = *n;
#line 495 "MB02CD.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 496 "MB02CD.f"
		startr = (i__ - 1) * *k + 1;
#line 497 "MB02CD.f"
		starti = (*n - i__ + 1) * *k + 1;
#line 498 "MB02CD.f"
		startt = (i__ - 2) * 3 * *k + 1;

/*              Transformations acting on the generator: */

#line 502 "MB02CD.f"
		i__2 = *k * 3;
#line 502 "MB02CD.f"
		mb02cx_("Column", k, k, k, &g[*k + 1 + (*k + 1) * g_dim1], 
			ldg, &t[startr + t_dim1], ldt, &cs[startt], &i__2, &
			dwork[1], ldwork, &ierr, (ftnlen)6);

#line 506 "MB02CD.f"
		if (ierr != 0) {

/*                 Error return:  The matrix is not positive definite. */

#line 510 "MB02CD.f"
		    *info = 1;
#line 511 "MB02CD.f"
		    return 0;
#line 512 "MB02CD.f"
		}

/* Computing MAX */
#line 514 "MB02CD.f"
		i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 514 "MB02CD.f"
		maxwrk = max(i__2,i__3);
#line 515 "MB02CD.f"
		if (*n > i__) {
#line 516 "MB02CD.f"
		    i__2 = (*n - i__) * *k;
#line 516 "MB02CD.f"
		    i__3 = *k * 3;
#line 516 "MB02CD.f"
		    mb02cy_("Column", "NoStructure", k, k, &i__2, k, &g[(*k <<
			     1) + 1 + (*k + 1) * g_dim1], ldg, &t[startr + *k 
			    + t_dim1], ldt, &t[startr + t_dim1], ldt, &cs[
			    startt], &i__3, &dwork[1], ldwork, &ierr, (ftnlen)
			    6, (ftnlen)11);
/* Computing MAX */
#line 520 "MB02CD.f"
		    i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 520 "MB02CD.f"
		    maxwrk = max(i__2,i__3);
#line 521 "MB02CD.f"
		}

#line 523 "MB02CD.f"
		if (compr) {
#line 524 "MB02CD.f"
		    i__2 = (*n - i__ + 1) * *k;
#line 524 "MB02CD.f"
		    dlacpy_("Lower", &i__2, k, &g[*k + 1 + (*k + 1) * g_dim1],
			     ldg, &r__[startr + startr * r_dim1], ldr, (
			    ftnlen)5);
#line 526 "MB02CD.f"
		}

/*              Transformations acting on the inverse generator: */

#line 530 "MB02CD.f"
		dlaset_("All", k, k, &c_b19, &c_b19, &g[starti + (*k + 1) * 
			g_dim1], ldg, (ftnlen)3);
#line 532 "MB02CD.f"
		i__2 = *k * 3;
#line 532 "MB02CD.f"
		mb02cy_("Column", "Triangular", k, k, k, k, &g[(*k + 1) * 
			g_dim1 + 1], ldg, &g[startr + g_dim1], ldg, &t[startr 
			+ t_dim1], ldt, &cs[startt], &i__2, &dwork[1], ldwork,
			 &ierr, (ftnlen)6, (ftnlen)10);
/* Computing MAX */
#line 536 "MB02CD.f"
		i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 536 "MB02CD.f"
		maxwrk = max(i__2,i__3);

#line 538 "MB02CD.f"
		i__2 = (i__ - 1) * *k;
#line 538 "MB02CD.f"
		i__3 = *k * 3;
#line 538 "MB02CD.f"
		mb02cy_("Column", "NoStructure", k, k, &i__2, k, &g[starti + (
			*k + 1) * g_dim1], ldg, &g[g_offset], ldg, &t[startr 
			+ t_dim1], ldt, &cs[startt], &i__3, &dwork[1], ldwork,
			 &ierr, (ftnlen)6, (ftnlen)11);
/* Computing MAX */
#line 541 "MB02CD.f"
		i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 541 "MB02CD.f"
		maxwrk = max(i__2,i__3);

#line 543 "MB02CD.f"
		if (compl) {
#line 544 "MB02CD.f"
		    i__2 = (i__ - 1) * *k;
#line 544 "MB02CD.f"
		    dlacpy_("All", &i__2, k, &g[starti + (*k + 1) * g_dim1], 
			    ldg, &l[startr * l_dim1 + 1], ldl, (ftnlen)3);
#line 546 "MB02CD.f"
		    dlacpy_("Upper", k, k, &g[(*k + 1) * g_dim1 + 1], ldg, &l[
			    (i__ - 1) * *k + 1 + startr * l_dim1], ldl, (
			    ftnlen)5);
#line 548 "MB02CD.f"
		}
#line 549 "MB02CD.f"
/* L30: */
#line 549 "MB02CD.f"
	    }

#line 551 "MB02CD.f"
	} else {

/*           Here R is used as working array for holding the generator. */
/*           Again, T contains the second column of the generator. */
/*           The current column of R contains the first column of the */
/*           generator. */

#line 558 "MB02CD.f"
	    if (*n > 1) {
#line 558 "MB02CD.f"
		i__1 = (*n - 1) * *k;
#line 558 "MB02CD.f"
		dlacpy_("Lower", &i__1, k, &t[t_offset], ldt, &r__[*k + 1 + (*
			k + 1) * r_dim1], ldr, (ftnlen)5);
#line 558 "MB02CD.f"
	    }

#line 562 "MB02CD.f"
	    i__1 = *n;
#line 562 "MB02CD.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 563 "MB02CD.f"
		startr = (i__ - 1) * *k + 1;
#line 564 "MB02CD.f"
		startt = (i__ - 2) * 3 * *k + 1;
#line 565 "MB02CD.f"
		i__2 = *k * 3;
#line 565 "MB02CD.f"
		mb02cx_("Column", k, k, k, &r__[startr + startr * r_dim1], 
			ldr, &t[startr + t_dim1], ldt, &cs[startt], &i__2, &
			dwork[1], ldwork, &ierr, (ftnlen)6);
#line 568 "MB02CD.f"
		if (ierr != 0) {

/*                 Error return:  The matrix is not positive definite. */

#line 572 "MB02CD.f"
		    *info = 1;
#line 573 "MB02CD.f"
		    return 0;
#line 574 "MB02CD.f"
		}

/* Computing MAX */
#line 576 "MB02CD.f"
		i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 576 "MB02CD.f"
		maxwrk = max(i__2,i__3);
#line 577 "MB02CD.f"
		if (*n > i__) {
#line 578 "MB02CD.f"
		    i__2 = (*n - i__) * *k;
#line 578 "MB02CD.f"
		    i__3 = *k * 3;
#line 578 "MB02CD.f"
		    mb02cy_("Column", "NoStructure", k, k, &i__2, k, &r__[
			    startr + *k + startr * r_dim1], ldr, &t[startr + *
			    k + t_dim1], ldt, &t[startr + t_dim1], ldt, &cs[
			    startt], &i__3, &dwork[1], ldwork, &ierr, (ftnlen)
			    6, (ftnlen)11);
/* Computing MAX */
#line 582 "MB02CD.f"
		    i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 582 "MB02CD.f"
		    maxwrk = max(i__2,i__3);

#line 584 "MB02CD.f"
		    i__2 = (*n - i__) * *k;
#line 584 "MB02CD.f"
		    dlacpy_("Lower", &i__2, k, &r__[startr + startr * r_dim1],
			     ldr, &r__[startr + *k + (startr + *k) * r_dim1], 
			    ldr, (ftnlen)5);
#line 586 "MB02CD.f"
		}
#line 587 "MB02CD.f"
/* L40: */
#line 587 "MB02CD.f"
	    }

#line 589 "MB02CD.f"
	}
#line 590 "MB02CD.f"
    }

#line 592 "MB02CD.f"
    dwork[1] = (doublereal) maxwrk;

#line 594 "MB02CD.f"
    return 0;

/* *** Last line of MB02CD *** */
} /* mb02cd_ */

