#line 1 "MB02DD.f"
/* MB02DD.f -- translated by f2c (version 20100827).
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

#line 1 "MB02DD.f"
/* Table of constant values */

static doublereal c_b12 = 1.;
static doublereal c_b14 = 0.;
static integer c__1 = 1;

/* Subroutine */ int mb02dd_(char *job, char *typet, integer *k, integer *m, 
	integer *n, doublereal *ta, integer *ldta, doublereal *t, integer *
	ldt, doublereal *g, integer *ldg, doublereal *r__, integer *ldr, 
	doublereal *l, integer *ldl, doublereal *cs, integer *lcs, doublereal 
	*dwork, integer *ldwork, integer *info, ftnlen job_len, ftnlen 
	typet_len)
{
    /* System generated locals */
    integer g_dim1, g_offset, l_dim1, l_offset, r_dim1, r_offset, t_dim1, 
	    t_offset, ta_dim1, ta_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, ierr;
    extern /* Subroutine */ int mb02cx_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen), mb02cy_(
	    char *, char *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical compg, compl;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical isrow;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
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

/*     To update the Cholesky factor and the generator and/or the */
/*     Cholesky factor of the inverse of a symmetric positive definite */
/*     (s.p.d.) block Toeplitz matrix T, given the information from */
/*     a previous factorization and additional blocks in TA of its first */
/*     block row, or its first block column, depending on the routine */
/*     parameter TYPET. Transformation information is stored. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the output of the routine, as follows: */
/*             = 'R':  updates the generator G of the inverse and */
/*                     computes the new columns / rows for the Cholesky */
/*                     factor R of T; */
/*             = 'A':  updates the generator G, computes the new */
/*                     columns / rows for the Cholesky factor R of T and */
/*                     the new rows / columns for the Cholesky factor L */
/*                     of the inverse; */
/*             = 'O':  only computes the new columns / rows for the */
/*                     Cholesky factor R of T. */

/*     TYPET   CHARACTER*1 */
/*             Specifies the type of T, as follows: */
/*             = 'R':  the first block row of an s.p.d. block Toeplitz */
/*                     matrix was/is defined; if demanded, the Cholesky */
/*                     factors R and L are upper and lower triangular, */
/*                     respectively, and G contains the transposed */
/*                     generator of the inverse; */
/*             = 'C':  the first block column of an s.p.d. block Toeplitz */
/*                     matrix was/is defined; if demanded, the Cholesky */
/*                     factors R and L are lower and upper triangular, */
/*                     respectively, and G contains the generator of the */
/*                     inverse. This choice results in a column oriented */
/*                     algorithm which is usually faster. */
/*             Note:   in this routine, the notation x / y means that */
/*                     x corresponds to TYPET = 'R' and y corresponds to */
/*                     TYPET = 'C'. */

/*     Input/Output Parameters */

/*     K       (input)  INTEGER */
/*             The number of rows / columns in T, which should be equal */
/*             to the blocksize.  K >= 0. */

/*     M       (input)  INTEGER */
/*             The number of blocks in TA.  M >= 0. */

/*     N       (input)  INTEGER */
/*             The number of blocks in T.  N >= 0. */

/*     TA      (input/output)  DOUBLE PRECISION array, dimension */
/*             (LDTA,M*K) / (LDTA,K) */
/*             On entry, the leading K-by-M*K / M*K-by-K part of this */
/*             array must contain the (N+1)-th to (N+M)-th blocks in the */
/*             first block row / column of an s.p.d. block Toeplitz */
/*             matrix. */
/*             On exit, if INFO = 0, the leading K-by-M*K / M*K-by-K part */
/*             of this array contains information on the Householder */
/*             transformations used, such that the array */

/*                        [ T  TA ]    /    [ T  ] */
/*                                          [ TA ] */

/*             serves as the new transformation matrix T for further */
/*             applications of this routine. */

/*     LDTA    INTEGER */
/*             The leading dimension of the array TA. */
/*             LDTA >= MAX(1,K),   if TYPET = 'R'; */
/*             LDTA >= MAX(1,M*K), if TYPET = 'C'. */

/*     T       (input)  DOUBLE PRECISION array, dimension (LDT,N*K) / */
/*             (LDT,K) */
/*             The leading K-by-N*K / N*K-by-K part of this array must */
/*             contain transformation information generated by the SLICOT */
/*             Library routine MB02CD, i.e., in the first K-by-K block, */
/*             the upper / lower Cholesky factor of T(1:K,1:K), and in */
/*             the remaining part, the Householder transformations */
/*             applied during the initial factorization process. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T. */
/*             LDT >= MAX(1,K),    if TYPET = 'R'; */
/*             LDT >= MAX(1,N*K),  if TYPET = 'C'. */

/*     G       (input/output)  DOUBLE PRECISION array, dimension */
/*             (LDG,( N + M )*K) / (LDG,2*K) */
/*             On entry, if JOB = 'R', or 'A', then the leading */
/*             2*K-by-N*K / N*K-by-2*K part of this array must contain, */
/*             in the first K-by-K block of the second block row / */
/*             column, the lower right block of the Cholesky factor of */
/*             the inverse of T, and in the remaining part, the generator */
/*             of the inverse of T. */
/*             On exit, if INFO = 0 and JOB = 'R', or 'A', then the */
/*             leading 2*K-by-( N + M )*K / ( N + M )*K-by-2*K part of */
/*             this array contains the same information as on entry, now */
/*             for the updated Toeplitz matrix. Actually, to obtain a */
/*             generator of the inverse one has to set */
/*               G(K+1:2*K, 1:K) = 0,    if TYPET = 'R'; */
/*               G(1:K, K+1:2*K) = 0,    if TYPET = 'C'. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G. */
/*             LDG >= MAX(1,2*K),  if TYPET = 'R' and JOB = 'R', or 'A'; */
/*             LDG >= MAX(1,( N + M )*K), */
/*                                 if TYPET = 'C' and JOB = 'R', or 'A'; */
/*             LDG >= 1,           if JOB = 'O'. */

/*     R       (input/output)  DOUBLE PRECISION array, dimension */
/*             (LDR,M*K) / (LDR,( N + M )*K) */
/*             On input, the leading N*K-by-K part of R(K+1,1) / */
/*             K-by-N*K part of R(1,K+1) contains the last block column / */
/*             row of the previous Cholesky factor R. */
/*             On exit, if INFO = 0, then the leading */
/*             ( N + M )*K-by-M*K / M*K-by-( N + M )*K part of this */
/*             array contains the last M*K columns / rows of the upper / */
/*             lower Cholesky factor of T. The elements in the strictly */
/*             lower / upper triangular part are not referenced. */

/*     LDR     INTEGER */
/*             The leading dimension of the array R. */
/*             LDR >= MAX(1, ( N + M )*K), if TYPET = 'R'; */
/*             LDR >= MAX(1, M*K),         if TYPET = 'C'. */

/*     L       (output)  DOUBLE PRECISION array, dimension */
/*             (LDL,( N + M )*K) / (LDL,M*K) */
/*             If INFO = 0 and JOB = 'A', then the leading */
/*             M*K-by-( N + M )*K / ( N + M )*K-by-M*K part of this */
/*             array contains the last M*K rows / columns of the lower / */
/*             upper Cholesky factor of the inverse of T. The elements */
/*             in the strictly upper / lower triangular part are not */
/*             referenced. */

/*     LDL     INTEGER */
/*             The leading dimension of the array L. */
/*             LDL >= MAX(1, M*K),         if TYPET = 'R' and JOB = 'A'; */
/*             LDL >= MAX(1, ( N + M )*K), if TYPET = 'C' and JOB = 'A'; */
/*             LDL >= 1,                   if JOB = 'R', or 'O'. */

/*     CS      (input/output)  DOUBLE PRECISION array, dimension (LCS) */
/*             On input, the leading 3*(N-1)*K part of this array must */
/*             contain the necessary information about the hyperbolic */
/*             rotations and Householder transformations applied */
/*             previously by SLICOT Library routine MB02CD. */
/*             On exit, if INFO = 0, then the leading 3*(N+M-1)*K part of */
/*             this array contains information about all the hyperbolic */
/*             rotations and Householder transformations applied during */
/*             the whole process. */

/*     LCS     INTEGER */
/*             The length of the array CS.  LCS >= 3*(N+M-1)*K. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -19,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,(N+M-1)*K). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the reduction algorithm failed. The block Toeplitz */
/*                   matrix associated with [ T  TA ] / [ T'  TA' ]' is */
/*                   not (numerically) positive definite. */

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
/*                               3         2 */
/*     The algorithm requires 0(K ( N M + M ) ) floating point */
/*     operations. */

/*     FURTHER COMMENTS */

/*     For min(K,N,M) = 0, the routine sets DWORK(1) = 1 and returns. */
/*     Although the calculations could still be performed when N = 0, */
/*     but min(K,M) > 0, this case is not considered as an "update". */
/*     SLICOT Library routine MB02CD should be called with the argument */
/*     M instead of N. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Chemnitz, Germany, December 2000. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000, */
/*     Feb. 2004. */

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

#line 275 "MB02DD.f"
    /* Parameter adjustments */
#line 275 "MB02DD.f"
    ta_dim1 = *ldta;
#line 275 "MB02DD.f"
    ta_offset = 1 + ta_dim1;
#line 275 "MB02DD.f"
    ta -= ta_offset;
#line 275 "MB02DD.f"
    t_dim1 = *ldt;
#line 275 "MB02DD.f"
    t_offset = 1 + t_dim1;
#line 275 "MB02DD.f"
    t -= t_offset;
#line 275 "MB02DD.f"
    g_dim1 = *ldg;
#line 275 "MB02DD.f"
    g_offset = 1 + g_dim1;
#line 275 "MB02DD.f"
    g -= g_offset;
#line 275 "MB02DD.f"
    r_dim1 = *ldr;
#line 275 "MB02DD.f"
    r_offset = 1 + r_dim1;
#line 275 "MB02DD.f"
    r__ -= r_offset;
#line 275 "MB02DD.f"
    l_dim1 = *ldl;
#line 275 "MB02DD.f"
    l_offset = 1 + l_dim1;
#line 275 "MB02DD.f"
    l -= l_offset;
#line 275 "MB02DD.f"
    --cs;
#line 275 "MB02DD.f"
    --dwork;
#line 275 "MB02DD.f"

#line 275 "MB02DD.f"
    /* Function Body */
#line 275 "MB02DD.f"
    *info = 0;
#line 276 "MB02DD.f"
    compl = lsame_(job, "A", (ftnlen)1, (ftnlen)1);
#line 277 "MB02DD.f"
    compg = lsame_(job, "R", (ftnlen)1, (ftnlen)1) || compl;
#line 278 "MB02DD.f"
    isrow = lsame_(typet, "R", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

#line 282 "MB02DD.f"
    if (! (compg || lsame_(job, "O", (ftnlen)1, (ftnlen)1))) {
#line 283 "MB02DD.f"
	*info = -1;
#line 284 "MB02DD.f"
    } else if (! (isrow || lsame_(typet, "C", (ftnlen)1, (ftnlen)1))) {
#line 285 "MB02DD.f"
	*info = -2;
#line 286 "MB02DD.f"
    } else if (*k < 0) {
#line 287 "MB02DD.f"
	*info = -3;
#line 288 "MB02DD.f"
    } else if (*m < 0) {
#line 289 "MB02DD.f"
	*info = -4;
#line 290 "MB02DD.f"
    } else if (*n < 0) {
#line 291 "MB02DD.f"
	*info = -5;
#line 292 "MB02DD.f"
    } else if (*ldta < 1 || isrow && *ldta < *k || ! isrow && *ldta < *m * *k)
	     {
#line 294 "MB02DD.f"
	*info = -7;
#line 295 "MB02DD.f"
    } else if (*ldt < 1 || isrow && *ldt < *k || ! isrow && *ldt < *n * *k) {
#line 297 "MB02DD.f"
	*info = -9;
#line 298 "MB02DD.f"
    } else if (compg && (isrow && *ldg < *k << 1 || ! isrow && *ldg < (*n + *
	    m) * *k) || *ldg < 1) {
#line 301 "MB02DD.f"
	*info = -11;
#line 302 "MB02DD.f"
    } else if (isrow && *ldr < (*n + *m) * *k || ! isrow && *ldr < *m * *k || 
	    *ldr < 1) {
#line 305 "MB02DD.f"
	*info = -13;
#line 306 "MB02DD.f"
    } else if (compl && (isrow && *ldl < *m * *k || ! isrow && *ldl < (*n + *
	    m) * *k) || *ldl < 1) {
#line 309 "MB02DD.f"
	*info = -15;
#line 310 "MB02DD.f"
    } else if (*lcs < (*n + *m - 1) * 3 * *k) {
#line 311 "MB02DD.f"
	*info = -17;
#line 312 "MB02DD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 312 "MB02DD.f"
	i__1 = 1, i__2 = (*n + *m - 1) * *k;
#line 312 "MB02DD.f"
	if (*ldwork < max(i__1,i__2)) {
/* Computing MAX */
#line 313 "MB02DD.f"
	    i__1 = 1, i__2 = (*n + *m - 1) * *k;
#line 313 "MB02DD.f"
	    dwork[1] = (doublereal) max(i__1,i__2);
#line 314 "MB02DD.f"
	    *info = -19;
#line 315 "MB02DD.f"
	}
#line 315 "MB02DD.f"
    }

/*     Return if there were illegal values. */

#line 319 "MB02DD.f"
    if (*info != 0) {
#line 320 "MB02DD.f"
	i__1 = -(*info);
#line 320 "MB02DD.f"
	xerbla_("MB02DD", &i__1, (ftnlen)6);
#line 321 "MB02DD.f"
	return 0;
#line 322 "MB02DD.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 326 "MB02DD.f"
    i__1 = min(*k,*n);
#line 326 "MB02DD.f"
    if (min(i__1,*m) == 0) {
#line 327 "MB02DD.f"
	dwork[1] = 1.;
#line 328 "MB02DD.f"
	return 0;
#line 329 "MB02DD.f"
    }

#line 331 "MB02DD.f"
    maxwrk = 1;
#line 332 "MB02DD.f"
    if (isrow) {

/*        Apply Cholesky factor of T(1:K, 1:K) on TA. */

#line 336 "MB02DD.f"
	i__1 = *m * *k;
#line 336 "MB02DD.f"
	dtrsm_("Left", "Upper", "Transpose", "NonUnit", k, &i__1, &c_b12, &t[
		t_offset], ldt, &ta[ta_offset], ldta, (ftnlen)4, (ftnlen)5, (
		ftnlen)9, (ftnlen)7);

/*        Initialize the output matrices. */

#line 341 "MB02DD.f"
	if (compg) {
#line 342 "MB02DD.f"
	    i__1 = *m * *k;
#line 342 "MB02DD.f"
	    dlaset_("All", k, &i__1, &c_b14, &c_b14, &g[(*n * *k + 1) * 
		    g_dim1 + 1], ldg, (ftnlen)3);
#line 343 "MB02DD.f"
	    if (*m >= *n - 1 && *n > 1) {
#line 344 "MB02DD.f"
		i__1 = (*n - 1) * *k;
#line 344 "MB02DD.f"
		dlacpy_("All", k, &i__1, &g[*k + 1 + (*k + 1) * g_dim1], ldg, 
			&g[*k + 1 + (*k * (*m + 1) + 1) * g_dim1], ldg, (
			ftnlen)3);
#line 346 "MB02DD.f"
	    } else {
#line 347 "MB02DD.f"
		i__1 = *k + 1;
#line 347 "MB02DD.f"
		for (i__ = *n * *k; i__ >= i__1; --i__) {
#line 348 "MB02DD.f"
		    dcopy_(k, &g[*k + 1 + i__ * g_dim1], &c__1, &g[*k + 1 + (*
			    m * *k + i__) * g_dim1], &c__1);
#line 349 "MB02DD.f"
/* L10: */
#line 349 "MB02DD.f"
		}
#line 350 "MB02DD.f"
	    }
#line 351 "MB02DD.f"
	    i__1 = *m * *k;
#line 351 "MB02DD.f"
	    dlaset_("All", k, &i__1, &c_b14, &c_b14, &g[*k + 1 + (*k + 1) * 
		    g_dim1], ldg, (ftnlen)3);
#line 352 "MB02DD.f"
	}

#line 354 "MB02DD.f"
	i__1 = *m * *k;
#line 354 "MB02DD.f"
	dlacpy_("All", k, &i__1, &ta[ta_offset], ldta, &r__[r_offset], ldr, (
		ftnlen)3);

/*        Apply the stored transformations on the new columns. */

#line 358 "MB02DD.f"
	i__1 = *n;
#line 358 "MB02DD.f"
	for (i__ = 2; i__ <= i__1; ++i__) {

/*           Copy the last M-1 blocks of the positive generator together; */
/*           the last M blocks of the negative generator are contained */
/*           in TA. */

#line 364 "MB02DD.f"
	    startr = (i__ - 1) * *k + 1;
#line 365 "MB02DD.f"
	    startt = (i__ - 2) * 3 * *k + 1;
#line 366 "MB02DD.f"
	    i__2 = (*m - 1) * *k;
#line 366 "MB02DD.f"
	    dlacpy_("All", k, &i__2, &r__[startr - *k + r_dim1], ldr, &r__[
		    startr + (*k + 1) * r_dim1], ldr, (ftnlen)3);

/*           Apply the transformations stored in T on the generator. */

#line 371 "MB02DD.f"
	    i__2 = *m * *k;
#line 371 "MB02DD.f"
	    i__3 = *k * 3;
#line 371 "MB02DD.f"
	    mb02cy_("Row", "NoStructure", k, k, &i__2, k, &r__[startr + 
		    r_dim1], ldr, &ta[ta_offset], ldta, &t[startr * t_dim1 + 
		    1], ldt, &cs[startt], &i__3, &dwork[1], ldwork, &ierr, (
		    ftnlen)3, (ftnlen)11);
/* Computing MAX */
#line 374 "MB02DD.f"
	    i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 374 "MB02DD.f"
	    maxwrk = max(i__2,i__3);
#line 375 "MB02DD.f"
/* L20: */
#line 375 "MB02DD.f"
	}

/*        Now, we have "normality" and can apply further M Schur steps. */

#line 379 "MB02DD.f"
	i__1 = *m;
#line 379 "MB02DD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Copy the first M-I+1 blocks of the positive generator */
/*           together; the first M-I+1 blocks of the negative generator */
/*           are contained in TA. */

#line 385 "MB02DD.f"
	    startt = (*n + i__ - 2) * 3 * *k + 1;
#line 386 "MB02DD.f"
	    starti = (*m - i__ + 1) * *k + 1;
#line 387 "MB02DD.f"
	    startr = (*n + i__ - 1) * *k + 1;
#line 388 "MB02DD.f"
	    if (i__ == 1) {
#line 389 "MB02DD.f"
		i__2 = (*m - 1) * *k;
#line 389 "MB02DD.f"
		dlacpy_("All", k, &i__2, &r__[startr - *k + r_dim1], ldr, &
			r__[startr + (*k + 1) * r_dim1], ldr, (ftnlen)3);
#line 391 "MB02DD.f"
	    } else {
#line 392 "MB02DD.f"
		i__2 = (*m - i__ + 1) * *k;
#line 392 "MB02DD.f"
		dlacpy_("Upper", k, &i__2, &r__[startr - *k + ((i__ - 2) * *k 
			+ 1) * r_dim1], ldr, &r__[startr + ((i__ - 1) * *k + 
			1) * r_dim1], ldr, (ftnlen)5);
#line 395 "MB02DD.f"
	    }

/*           Reduce the generator to proper form. */

#line 399 "MB02DD.f"
	    i__2 = *k * 3;
#line 399 "MB02DD.f"
	    mb02cx_("Row", k, k, k, &r__[startr + ((i__ - 1) * *k + 1) * 
		    r_dim1], ldr, &ta[((i__ - 1) * *k + 1) * ta_dim1 + 1], 
		    ldta, &cs[startt], &i__2, &dwork[1], ldwork, &ierr, (
		    ftnlen)3);
#line 402 "MB02DD.f"
	    if (ierr != 0) {

/*              Error return:  The matrix is not positive definite. */

#line 406 "MB02DD.f"
		*info = 1;
#line 407 "MB02DD.f"
		return 0;
#line 408 "MB02DD.f"
	    }

/* Computing MAX */
#line 410 "MB02DD.f"
	    i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 410 "MB02DD.f"
	    maxwrk = max(i__2,i__3);
#line 411 "MB02DD.f"
	    if (*m > i__) {
#line 412 "MB02DD.f"
		i__2 = (*m - i__) * *k;
#line 412 "MB02DD.f"
		i__3 = *k * 3;
#line 412 "MB02DD.f"
		mb02cy_("Row", "NoStructure", k, k, &i__2, k, &r__[startr + (
			i__ * *k + 1) * r_dim1], ldr, &ta[(i__ * *k + 1) * 
			ta_dim1 + 1], ldta, &ta[((i__ - 1) * *k + 1) * 
			ta_dim1 + 1], ldta, &cs[startt], &i__3, &dwork[1], 
			ldwork, &ierr, (ftnlen)3, (ftnlen)11);
/* Computing MAX */
#line 416 "MB02DD.f"
		i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 416 "MB02DD.f"
		maxwrk = max(i__2,i__3);
#line 417 "MB02DD.f"
	    }

#line 419 "MB02DD.f"
	    if (compg) {

/*              Transformations acting on the inverse generator: */

#line 423 "MB02DD.f"
		i__2 = *k * 3;
#line 423 "MB02DD.f"
		mb02cy_("Row", "Triangular", k, k, k, k, &g[*k + 1 + g_dim1], 
			ldg, &g[startr * g_dim1 + 1], ldg, &ta[((i__ - 1) * *
			k + 1) * ta_dim1 + 1], ldta, &cs[startt], &i__2, &
			dwork[1], ldwork, &ierr, (ftnlen)3, (ftnlen)10);
/* Computing MAX */
#line 426 "MB02DD.f"
		i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 426 "MB02DD.f"
		maxwrk = max(i__2,i__3);

#line 428 "MB02DD.f"
		i__2 = (*n + i__ - 1) * *k;
#line 428 "MB02DD.f"
		i__3 = *k * 3;
#line 428 "MB02DD.f"
		mb02cy_("Row", "NoStructure", k, k, &i__2, k, &g[*k + 1 + 
			starti * g_dim1], ldg, &g[g_offset], ldg, &ta[((i__ - 
			1) * *k + 1) * ta_dim1 + 1], ldta, &cs[startt], &i__3,
			 &dwork[1], ldwork, &ierr, (ftnlen)3, (ftnlen)11);
/* Computing MAX */
#line 431 "MB02DD.f"
		i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 431 "MB02DD.f"
		maxwrk = max(i__2,i__3);

#line 433 "MB02DD.f"
		if (compl) {
#line 434 "MB02DD.f"
		    i__2 = (*n + i__ - 1) * *k;
#line 434 "MB02DD.f"
		    dlacpy_("All", k, &i__2, &g[*k + 1 + starti * g_dim1], 
			    ldg, &l[(i__ - 1) * *k + 1 + l_dim1], ldl, (
			    ftnlen)3);
#line 436 "MB02DD.f"
		    dlacpy_("Lower", k, k, &g[*k + 1 + g_dim1], ldg, &l[(i__ 
			    - 1) * *k + 1 + startr * l_dim1], ldl, (ftnlen)5);
#line 438 "MB02DD.f"
		}

#line 440 "MB02DD.f"
	    }
#line 441 "MB02DD.f"
/* L30: */
#line 441 "MB02DD.f"
	}

#line 443 "MB02DD.f"
    } else {

/*        Apply Cholesky factor of T(1:K, 1:K) on TA. */

#line 447 "MB02DD.f"
	i__1 = *m * *k;
#line 447 "MB02DD.f"
	dtrsm_("Right", "Lower", "Transpose", "NonUnit", &i__1, k, &c_b12, &t[
		t_offset], ldt, &ta[ta_offset], ldta, (ftnlen)5, (ftnlen)5, (
		ftnlen)9, (ftnlen)7);

/*        Initialize the output matrices. */

#line 452 "MB02DD.f"
	if (compg) {
#line 453 "MB02DD.f"
	    i__1 = *m * *k;
#line 453 "MB02DD.f"
	    dlaset_("All", &i__1, k, &c_b14, &c_b14, &g[*n * *k + 1 + g_dim1],
		     ldg, (ftnlen)3);
#line 454 "MB02DD.f"
	    if (*m >= *n - 1 && *n > 1) {
#line 455 "MB02DD.f"
		i__1 = (*n - 1) * *k;
#line 455 "MB02DD.f"
		dlacpy_("All", &i__1, k, &g[*k + 1 + (*k + 1) * g_dim1], ldg, 
			&g[*k * (*m + 1) + 1 + (*k + 1) * g_dim1], ldg, (
			ftnlen)3);
#line 457 "MB02DD.f"
	    } else {
#line 458 "MB02DD.f"
		i__1 = *k;
#line 458 "MB02DD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 459 "MB02DD.f"
		    i__2 = *k + 1;
#line 459 "MB02DD.f"
		    for (j = *n * *k; j >= i__2; --j) {
#line 460 "MB02DD.f"
			g[j + *m * *k + (*k + i__) * g_dim1] = g[j + (*k + 
				i__) * g_dim1];
#line 461 "MB02DD.f"
/* L35: */
#line 461 "MB02DD.f"
		    }
#line 462 "MB02DD.f"
/* L40: */
#line 462 "MB02DD.f"
		}
#line 463 "MB02DD.f"
	    }
#line 464 "MB02DD.f"
	    i__1 = *m * *k;
#line 464 "MB02DD.f"
	    dlaset_("All", &i__1, k, &c_b14, &c_b14, &g[*k + 1 + (*k + 1) * 
		    g_dim1], ldg, (ftnlen)3);
#line 465 "MB02DD.f"
	}

#line 467 "MB02DD.f"
	i__1 = *m * *k;
#line 467 "MB02DD.f"
	dlacpy_("All", &i__1, k, &ta[ta_offset], ldta, &r__[r_offset], ldr, (
		ftnlen)3);

/*        Apply the stored transformations on the new rows. */

#line 471 "MB02DD.f"
	i__1 = *n;
#line 471 "MB02DD.f"
	for (i__ = 2; i__ <= i__1; ++i__) {

/*           Copy the last M-1 blocks of the positive generator together; */
/*           the last M blocks of the negative generator are contained */
/*           in TA. */

#line 477 "MB02DD.f"
	    startr = (i__ - 1) * *k + 1;
#line 478 "MB02DD.f"
	    startt = (i__ - 2) * 3 * *k + 1;
#line 479 "MB02DD.f"
	    i__2 = (*m - 1) * *k;
#line 479 "MB02DD.f"
	    dlacpy_("All", &i__2, k, &r__[(startr - *k) * r_dim1 + 1], ldr, &
		    r__[*k + 1 + startr * r_dim1], ldr, (ftnlen)3);

/*           Apply the transformations stored in T on the generator. */

#line 484 "MB02DD.f"
	    i__2 = *m * *k;
#line 484 "MB02DD.f"
	    i__3 = *k * 3;
#line 484 "MB02DD.f"
	    mb02cy_("Column", "NoStructure", k, k, &i__2, k, &r__[startr * 
		    r_dim1 + 1], ldr, &ta[ta_offset], ldta, &t[startr + 
		    t_dim1], ldt, &cs[startt], &i__3, &dwork[1], ldwork, &
		    ierr, (ftnlen)6, (ftnlen)11);
/* Computing MAX */
#line 487 "MB02DD.f"
	    i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 487 "MB02DD.f"
	    maxwrk = max(i__2,i__3);
#line 488 "MB02DD.f"
/* L50: */
#line 488 "MB02DD.f"
	}

/*        Now, we have "normality" and can apply further M Schur steps. */

#line 492 "MB02DD.f"
	i__1 = *m;
#line 492 "MB02DD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Copy the first M-I+1 blocks of the positive generator */
/*           together; the first M-I+1 blocks of the negative generator */
/*           are contained in TA. */

#line 498 "MB02DD.f"
	    startt = (*n + i__ - 2) * 3 * *k + 1;
#line 499 "MB02DD.f"
	    starti = (*m - i__ + 1) * *k + 1;
#line 500 "MB02DD.f"
	    startr = (*n + i__ - 1) * *k + 1;
#line 501 "MB02DD.f"
	    if (i__ == 1) {
#line 502 "MB02DD.f"
		i__2 = (*m - 1) * *k;
#line 502 "MB02DD.f"
		dlacpy_("All", &i__2, k, &r__[(startr - *k) * r_dim1 + 1], 
			ldr, &r__[*k + 1 + startr * r_dim1], ldr, (ftnlen)3);
#line 504 "MB02DD.f"
	    } else {
#line 505 "MB02DD.f"
		i__2 = (*m - i__ + 1) * *k;
#line 505 "MB02DD.f"
		dlacpy_("Lower", &i__2, k, &r__[(i__ - 2) * *k + 1 + (startr 
			- *k) * r_dim1], ldr, &r__[(i__ - 1) * *k + 1 + 
			startr * r_dim1], ldr, (ftnlen)5);
#line 508 "MB02DD.f"
	    }

/*           Reduce the generator to proper form. */

#line 512 "MB02DD.f"
	    i__2 = *k * 3;
#line 512 "MB02DD.f"
	    mb02cx_("Column", k, k, k, &r__[(i__ - 1) * *k + 1 + startr * 
		    r_dim1], ldr, &ta[(i__ - 1) * *k + 1 + ta_dim1], ldta, &
		    cs[startt], &i__2, &dwork[1], ldwork, &ierr, (ftnlen)6);
#line 515 "MB02DD.f"
	    if (ierr != 0) {

/*              Error return:  The matrix is not positive definite. */

#line 519 "MB02DD.f"
		*info = 1;
#line 520 "MB02DD.f"
		return 0;
#line 521 "MB02DD.f"
	    }

/* Computing MAX */
#line 523 "MB02DD.f"
	    i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 523 "MB02DD.f"
	    maxwrk = max(i__2,i__3);
#line 524 "MB02DD.f"
	    if (*m > i__) {
#line 525 "MB02DD.f"
		i__2 = (*m - i__) * *k;
#line 525 "MB02DD.f"
		i__3 = *k * 3;
#line 525 "MB02DD.f"
		mb02cy_("Column", "NoStructure", k, k, &i__2, k, &r__[i__ * *
			k + 1 + startr * r_dim1], ldr, &ta[i__ * *k + 1 + 
			ta_dim1], ldta, &ta[(i__ - 1) * *k + 1 + ta_dim1], 
			ldta, &cs[startt], &i__3, &dwork[1], ldwork, &ierr, (
			ftnlen)6, (ftnlen)11);
/* Computing MAX */
#line 529 "MB02DD.f"
		i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 529 "MB02DD.f"
		maxwrk = max(i__2,i__3);
#line 530 "MB02DD.f"
	    }

#line 532 "MB02DD.f"
	    if (compg) {

/*              Transformations acting on the inverse generator: */

#line 536 "MB02DD.f"
		i__2 = *k * 3;
#line 536 "MB02DD.f"
		mb02cy_("Column", "Triangular", k, k, k, k, &g[(*k + 1) * 
			g_dim1 + 1], ldg, &g[startr + g_dim1], ldg, &ta[(i__ 
			- 1) * *k + 1 + ta_dim1], ldta, &cs[startt], &i__2, &
			dwork[1], ldwork, &ierr, (ftnlen)6, (ftnlen)10);
/* Computing MAX */
#line 540 "MB02DD.f"
		i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 540 "MB02DD.f"
		maxwrk = max(i__2,i__3);

#line 542 "MB02DD.f"
		i__2 = (*n + i__ - 1) * *k;
#line 542 "MB02DD.f"
		i__3 = *k * 3;
#line 542 "MB02DD.f"
		mb02cy_("Column", "NoStructure", k, k, &i__2, k, &g[starti + (
			*k + 1) * g_dim1], ldg, &g[g_offset], ldg, &ta[(i__ - 
			1) * *k + 1 + ta_dim1], ldta, &cs[startt], &i__3, &
			dwork[1], ldwork, &ierr, (ftnlen)6, (ftnlen)11);
/* Computing MAX */
#line 545 "MB02DD.f"
		i__2 = maxwrk, i__3 = (integer) dwork[1];
#line 545 "MB02DD.f"
		maxwrk = max(i__2,i__3);

#line 547 "MB02DD.f"
		if (compl) {
#line 548 "MB02DD.f"
		    i__2 = (*n + i__ - 1) * *k;
#line 548 "MB02DD.f"
		    dlacpy_("All", &i__2, k, &g[starti + (*k + 1) * g_dim1], 
			    ldg, &l[((i__ - 1) * *k + 1) * l_dim1 + 1], ldl, (
			    ftnlen)3);
#line 550 "MB02DD.f"
		    dlacpy_("Upper", k, k, &g[(*k + 1) * g_dim1 + 1], ldg, &l[
			    startr + ((i__ - 1) * *k + 1) * l_dim1], ldl, (
			    ftnlen)5);
#line 552 "MB02DD.f"
		}

#line 554 "MB02DD.f"
	    }
#line 555 "MB02DD.f"
/* L60: */
#line 555 "MB02DD.f"
	}

#line 557 "MB02DD.f"
    }

#line 559 "MB02DD.f"
    dwork[1] = (doublereal) maxwrk;

#line 561 "MB02DD.f"
    return 0;

/* *** Last line of MB02DD *** */
} /* mb02dd_ */

