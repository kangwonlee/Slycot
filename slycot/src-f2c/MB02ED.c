#line 1 "MB02ED.f"
/* MB02ED.f -- translated by f2c (version 20100827).
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

#line 1 "MB02ED.f"
/* Table of constant values */

static doublereal c_b10 = 1.;
static doublereal c_b19 = -1.;
static doublereal c_b21 = 0.;

/* Subroutine */ int mb02ed_(char *typet, integer *k, integer *n, integer *
	nrhs, doublereal *t, integer *ldt, doublereal *b, integer *ldb, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen typet_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;

    /* Local variables */
    static integer i__, ierr;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     mb02cx_(char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen), mb02cy_(char *, char 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dtrsm_(
	    char *, char *, char *, char *, integer *, integer *, doublereal *
	    , doublereal *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static logical isrow;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), dpotrf_(char *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static integer starth, starti, maxwrk, startn, startr, startt;


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

/*     To solve a system of linear equations  T*X = B  or  X*T = B  with */
/*     a symmetric positive definite (s.p.d.) block Toeplitz matrix T. */
/*     T is defined either by its first block row or its first block */
/*     column, depending on the parameter TYPET. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYPET   CHARACTER*1 */
/*             Specifies the type of T, as follows: */
/*             = 'R':  T contains the first block row of an s.p.d. block */
/*                     Toeplitz matrix, and the system X*T = B is solved; */
/*             = 'C':  T contains the first block column of an s.p.d. */
/*                     block Toeplitz matrix, and the system T*X = B is */
/*                     solved. */
/*             Note:   in the sequel, the notation x / y means that */
/*                     x corresponds to TYPET = 'R' and y corresponds to */
/*                     TYPET = 'C'. */

/*     Input/Output Parameters */

/*     K       (input)  INTEGER */
/*             The number of rows / columns in T, which should be equal */
/*             to the blocksize.  K >= 0. */

/*     N       (input)  INTEGER */
/*             The number of blocks in T.  N >= 0. */

/*     NRHS    (input)  INTEGER */
/*             The number of right hand sides.  NRHS >= 0. */

/*     T       (input/output)  DOUBLE PRECISION array, dimension */
/*             (LDT,N*K) / (LDT,K) */
/*             On entry, the leading K-by-N*K / N*K-by-K part of this */
/*             array must contain the first block row / column of an */
/*             s.p.d. block Toeplitz matrix. */
/*             On exit, if  INFO = 0  and  NRHS > 0,  then the leading */
/*             K-by-N*K / N*K-by-K part of this array contains the last */
/*             row / column of the Cholesky factor of inv(T). */

/*     LDT     INTEGER */
/*             The leading dimension of the array T. */
/*             LDT >= MAX(1,K),    if TYPET = 'R'; */
/*             LDT >= MAX(1,N*K),  if TYPET = 'C'. */

/*     B       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDB,N*K) / (LDB,NRHS) */
/*             On entry, the leading NRHS-by-N*K / N*K-by-NRHS part of */
/*             this array must contain the right hand side matrix B. */
/*             On exit, the leading NRHS-by-N*K / N*K-by-NRHS part of */
/*             this array contains the solution matrix X. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= MAX(1,NRHS),  if TYPET = 'R'; */
/*             LDB >= MAX(1,N*K),   if TYPET = 'C'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -10,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,N*K*K+(N+2)*K). */
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

/*     Householder transformations, modified hyperbolic rotations and */
/*     block Gaussian eliminations are used in the Schur algorithm [1], */
/*     [2]. */

/*     REFERENCES */

/*     [1] Kailath, T. and Sayed, A. */
/*         Fast Reliable Algorithms for Matrices with Structure. */
/*         SIAM Publications, Philadelphia, 1999. */

/*     [2] Kressner, D. and Van Dooren, P. */
/*         Factorizations and linear system solvers for matrices with */
/*         Toeplitz structure. */
/*         SLICOT Working Note 2000-2, 2000. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is numerically equivalent with forming */
/*     the Cholesky factor R and the inverse Cholesky factor of T, using */
/*     the generalized Schur algorithm, and solving the systems of */
/*     equations  R*X = L*B  or  X*R = B*L by a blocked backward */
/*     substitution algorithm. */
/*                               3 2    2 2 */
/*     The algorithm requires 0(K N  + K N NRHS) floating point */
/*     operations. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Chemnitz, Germany, December 2000. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000, */
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

#line 174 "MB02ED.f"
    /* Parameter adjustments */
#line 174 "MB02ED.f"
    t_dim1 = *ldt;
#line 174 "MB02ED.f"
    t_offset = 1 + t_dim1;
#line 174 "MB02ED.f"
    t -= t_offset;
#line 174 "MB02ED.f"
    b_dim1 = *ldb;
#line 174 "MB02ED.f"
    b_offset = 1 + b_dim1;
#line 174 "MB02ED.f"
    b -= b_offset;
#line 174 "MB02ED.f"
    --dwork;
#line 174 "MB02ED.f"

#line 174 "MB02ED.f"
    /* Function Body */
#line 174 "MB02ED.f"
    *info = 0;
#line 175 "MB02ED.f"
    isrow = lsame_(typet, "R", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

#line 179 "MB02ED.f"
    if (! (isrow || lsame_(typet, "C", (ftnlen)1, (ftnlen)1))) {
#line 180 "MB02ED.f"
	*info = -1;
#line 181 "MB02ED.f"
    } else if (*k < 0) {
#line 182 "MB02ED.f"
	*info = -2;
#line 183 "MB02ED.f"
    } else if (*n < 0) {
#line 184 "MB02ED.f"
	*info = -3;
#line 185 "MB02ED.f"
    } else if (*nrhs < 0) {
#line 186 "MB02ED.f"
	*info = -4;
#line 187 "MB02ED.f"
    } else if (*ldt < 1 || isrow && *ldt < *k || ! isrow && *ldt < *n * *k) {
#line 189 "MB02ED.f"
	*info = -6;
#line 190 "MB02ED.f"
    } else if (*ldb < 1 || isrow && *ldb < *nrhs || ! isrow && *ldb < *n * *k)
	     {
#line 192 "MB02ED.f"
	*info = -8;
#line 193 "MB02ED.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 193 "MB02ED.f"
	i__1 = 1, i__2 = *n * *k * *k + (*n + 2) * *k;
#line 193 "MB02ED.f"
	if (*ldwork < max(i__1,i__2)) {
/* Computing MAX */
#line 194 "MB02ED.f"
	    i__1 = 1, i__2 = *n * *k * *k + (*n + 2) * *k;
#line 194 "MB02ED.f"
	    dwork[1] = (doublereal) max(i__1,i__2);
#line 195 "MB02ED.f"
	    *info = -10;
#line 196 "MB02ED.f"
	}
#line 196 "MB02ED.f"
    }

/*     Return if there were illegal values. */

#line 200 "MB02ED.f"
    if (*info != 0) {
#line 201 "MB02ED.f"
	i__1 = -(*info);
#line 201 "MB02ED.f"
	xerbla_("MB02ED", &i__1, (ftnlen)6);
#line 202 "MB02ED.f"
	return 0;
#line 203 "MB02ED.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 207 "MB02ED.f"
    i__1 = min(*k,*n);
#line 207 "MB02ED.f"
    if (min(i__1,*nrhs) == 0) {
#line 208 "MB02ED.f"
	dwork[1] = 1.;
#line 209 "MB02ED.f"
	return 0;
#line 210 "MB02ED.f"
    }

#line 212 "MB02ED.f"
    maxwrk = 0;
#line 213 "MB02ED.f"
    startn = 1;
#line 214 "MB02ED.f"
    startt = *n * *k * *k + 1;
#line 215 "MB02ED.f"
    starth = startt + *k * 3;

#line 217 "MB02ED.f"
    if (isrow) {

/*        T is the first block row of a block Toeplitz matrix. */
/*        Bring T to proper form by triangularizing its first block. */

#line 222 "MB02ED.f"
	dpotrf_("Upper", k, &t[t_offset], ldt, &ierr, (ftnlen)5);
#line 223 "MB02ED.f"
	if (ierr != 0) {

/*           Error return:  The matrix is not positive definite. */

#line 227 "MB02ED.f"
	    *info = 1;
#line 228 "MB02ED.f"
	    return 0;
#line 229 "MB02ED.f"
	}

#line 231 "MB02ED.f"
	if (*n > 1) {
#line 231 "MB02ED.f"
	    i__1 = (*n - 1) * *k;
#line 231 "MB02ED.f"
	    dtrsm_("Left", "Upper", "Transpose", "NonUnit", k, &i__1, &c_b10, 
		    &t[t_offset], ldt, &t[(*k + 1) * t_dim1 + 1], ldt, (
		    ftnlen)4, (ftnlen)5, (ftnlen)9, (ftnlen)7);
#line 231 "MB02ED.f"
	}

/*        Initialize the generator, do the first Schur step and set */
/*        B = -B. */
/*        T contains the nonzero blocks of the positive parts in the */
/*        generator and the inverse generator. */
/*        DWORK(STARTN) contains the nonzero blocks of the negative parts */
/*        in the generator and the inverse generator. */

#line 242 "MB02ED.f"
	dtrsm_("Right", "Upper", "NonTranspose", "NonUnit", nrhs, k, &c_b10, &
		t[t_offset], ldt, &b[b_offset], ldb, (ftnlen)5, (ftnlen)5, (
		ftnlen)12, (ftnlen)7);
#line 244 "MB02ED.f"
	if (*n > 1) {
#line 244 "MB02ED.f"
	    i__1 = (*n - 1) * *k;
#line 244 "MB02ED.f"
	    dgemm_("NonTranspose", "NonTranspose", nrhs, &i__1, k, &c_b10, &b[
		    b_offset], ldb, &t[(*k + 1) * t_dim1 + 1], ldt, &c_b19, &
		    b[(*k + 1) * b_dim1 + 1], ldb, (ftnlen)12, (ftnlen)12);
#line 244 "MB02ED.f"
	}

#line 249 "MB02ED.f"
	dlaset_("All", k, k, &c_b21, &c_b10, &dwork[startn], k, (ftnlen)3);
#line 250 "MB02ED.f"
	dtrsm_("Left", "Upper", "Transpose", "NonUnit", k, k, &c_b10, &t[
		t_offset], ldt, &dwork[startn], k, (ftnlen)4, (ftnlen)5, (
		ftnlen)9, (ftnlen)7);
#line 252 "MB02ED.f"
	if (*n > 1) {
#line 252 "MB02ED.f"
	    i__1 = (*n - 1) * *k;
#line 252 "MB02ED.f"
	    dlacpy_("All", k, &i__1, &t[(*k + 1) * t_dim1 + 1], ldt, &dwork[
		    startn + *k * *k], k, (ftnlen)3);
#line 252 "MB02ED.f"
	}
#line 255 "MB02ED.f"
	dlacpy_("All", k, k, &dwork[startn], k, &t[((*n - 1) * *k + 1) * 
		t_dim1 + 1], ldt, (ftnlen)3);

#line 258 "MB02ED.f"
	dtrmm_("Right", "Lower", "NonTranspose", "NonUnit", nrhs, k, &c_b10, &
		t[((*n - 1) * *k + 1) * t_dim1 + 1], ldt, &b[b_offset], ldb, (
		ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)7);

/*        Processing the generator. */

#line 263 "MB02ED.f"
	i__1 = *n;
#line 263 "MB02ED.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 264 "MB02ED.f"
	    startr = (i__ - 1) * *k + 1;
#line 265 "MB02ED.f"
	    starti = (*n - i__) * *k + 1;

/*           Transform the generator of T to proper form. */

#line 269 "MB02ED.f"
	    i__2 = *k * 3;
#line 269 "MB02ED.f"
	    i__3 = *ldwork - starth + 1;
#line 269 "MB02ED.f"
	    mb02cx_("Row", k, k, k, &t[t_offset], ldt, &dwork[startn + (i__ - 
		    1) * *k * *k], k, &dwork[startt], &i__2, &dwork[starth], &
		    i__3, &ierr, (ftnlen)3);

#line 273 "MB02ED.f"
	    if (ierr != 0) {

/*              Error return:  The matrix is not positive definite. */

#line 277 "MB02ED.f"
		*info = 1;
#line 278 "MB02ED.f"
		return 0;
#line 279 "MB02ED.f"
	    }

/* Computing MAX */
#line 281 "MB02ED.f"
	    i__2 = maxwrk, i__3 = (integer) dwork[starth];
#line 281 "MB02ED.f"
	    maxwrk = max(i__2,i__3);
#line 282 "MB02ED.f"
	    i__2 = (*n - i__) * *k;
#line 282 "MB02ED.f"
	    i__3 = *k * 3;
#line 282 "MB02ED.f"
	    i__4 = *ldwork - starth + 1;
#line 282 "MB02ED.f"
	    mb02cy_("Row", "NoStructure", k, k, &i__2, k, &t[(*k + 1) * 
		    t_dim1 + 1], ldt, &dwork[startn + i__ * *k * *k], k, &
		    dwork[startn + (i__ - 1) * *k * *k], k, &dwork[startt], &
		    i__3, &dwork[starth], &i__4, &ierr, (ftnlen)3, (ftnlen)11)
		    ;
/* Computing MAX */
#line 286 "MB02ED.f"
	    i__2 = maxwrk, i__3 = (integer) dwork[starth];
#line 286 "MB02ED.f"
	    maxwrk = max(i__2,i__3);

/*           Block Gaussian eliminates the i-th block in B. */

#line 290 "MB02ED.f"
	    dtrsm_("Right", "Upper", "NonTranspose", "NonUnit", nrhs, k, &
		    c_b19, &t[t_offset], ldt, &b[startr * b_dim1 + 1], ldb, (
		    ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)7);
#line 292 "MB02ED.f"
	    if (*n > i__) {
#line 292 "MB02ED.f"
		i__2 = (*n - i__) * *k;
#line 292 "MB02ED.f"
		dgemm_("NonTranspose", "NonTranspose", nrhs, &i__2, k, &c_b10,
			 &b[startr * b_dim1 + 1], ldb, &t[(*k + 1) * t_dim1 + 
			1], ldt, &c_b10, &b[(startr + *k) * b_dim1 + 1], ldb, 
			(ftnlen)12, (ftnlen)12);
#line 292 "MB02ED.f"
	    }

/*           Apply hyperbolic transformations on the negative generator. */

#line 299 "MB02ED.f"
	    dlaset_("All", k, k, &c_b21, &c_b21, &t[starti * t_dim1 + 1], ldt,
		     (ftnlen)3);
#line 300 "MB02ED.f"
	    i__2 = (i__ - 1) * *k;
#line 300 "MB02ED.f"
	    i__3 = *k * 3;
#line 300 "MB02ED.f"
	    i__4 = *ldwork - starth + 1;
#line 300 "MB02ED.f"
	    mb02cy_("Row", "NoStructure", k, k, &i__2, k, &t[starti * t_dim1 
		    + 1], ldt, &dwork[startn], k, &dwork[startn + (i__ - 1) * 
		    *k * *k], k, &dwork[startt], &i__3, &dwork[starth], &i__4,
		     &ierr, (ftnlen)3, (ftnlen)11);
/* Computing MAX */
#line 304 "MB02ED.f"
	    i__2 = maxwrk, i__3 = (integer) dwork[starth];
#line 304 "MB02ED.f"
	    maxwrk = max(i__2,i__3);

/*           Note that  DWORK(STARTN+(I-1)*K*K)  serves simultaneously */
/*           as the transformation container as well as the new block in */
/*           the negative generator. */

#line 310 "MB02ED.f"
	    i__2 = *k * 3;
#line 310 "MB02ED.f"
	    i__3 = *ldwork - starth + 1;
#line 310 "MB02ED.f"
	    mb02cy_("Row", "Triangular", k, k, k, k, &t[((*n - 1) * *k + 1) * 
		    t_dim1 + 1], ldt, &dwork[startn + (i__ - 1) * *k * *k], k,
		     &dwork[startn + (i__ - 1) * *k * *k], k, &dwork[startt], 
		    &i__2, &dwork[starth], &i__3, &ierr, (ftnlen)3, (ftnlen)
		    10);
/* Computing MAX */
#line 314 "MB02ED.f"
	    i__2 = maxwrk, i__3 = (integer) dwork[starth];
#line 314 "MB02ED.f"
	    maxwrk = max(i__2,i__3);

/*           Finally the Gaussian elimination is applied on the inverse */
/*           generator. */

#line 319 "MB02ED.f"
	    i__2 = (i__ - 1) * *k;
#line 319 "MB02ED.f"
	    dgemm_("NonTranspose", "NonTranspose", nrhs, &i__2, k, &c_b10, &b[
		    startr * b_dim1 + 1], ldb, &t[starti * t_dim1 + 1], ldt, &
		    c_b10, &b[b_offset], ldb, (ftnlen)12, (ftnlen)12);
#line 322 "MB02ED.f"
	    dtrmm_("Right", "Lower", "NonTranspose", "NonUnit", nrhs, k, &
		    c_b10, &t[((*n - 1) * *k + 1) * t_dim1 + 1], ldt, &b[
		    startr * b_dim1 + 1], ldb, (ftnlen)5, (ftnlen)5, (ftnlen)
		    12, (ftnlen)7);
#line 325 "MB02ED.f"
/* L10: */
#line 325 "MB02ED.f"
	}

#line 327 "MB02ED.f"
    } else {

/*        T is the first block column of a block Toeplitz matrix. */
/*        Bring T to proper form by triangularizing its first block. */

#line 332 "MB02ED.f"
	dpotrf_("Lower", k, &t[t_offset], ldt, &ierr, (ftnlen)5);
#line 333 "MB02ED.f"
	if (ierr != 0) {

/*           Error return:  The matrix is not positive definite. */

#line 337 "MB02ED.f"
	    *info = 1;
#line 338 "MB02ED.f"
	    return 0;
#line 339 "MB02ED.f"
	}

#line 341 "MB02ED.f"
	if (*n > 1) {
#line 341 "MB02ED.f"
	    i__1 = (*n - 1) * *k;
#line 341 "MB02ED.f"
	    dtrsm_("Right", "Lower", "Transpose", "NonUnit", &i__1, k, &c_b10,
		     &t[t_offset], ldt, &t[*k + 1 + t_dim1], ldt, (ftnlen)5, (
		    ftnlen)5, (ftnlen)9, (ftnlen)7);
#line 341 "MB02ED.f"
	}

/*        Initialize the generator, do the first Schur step and set */
/*        B = -B. */
/*        T contains the nonzero blocks of the positive parts in the */
/*        generator and the inverse generator. */
/*        DWORK(STARTN) contains the nonzero blocks of the negative parts */
/*        in the generator and the inverse generator. */

#line 352 "MB02ED.f"
	dtrsm_("Left", "Lower", "NonTranspose", "NonUnit", k, nrhs, &c_b10, &
		t[t_offset], ldt, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
		ftnlen)12, (ftnlen)7);
#line 354 "MB02ED.f"
	if (*n > 1) {
#line 354 "MB02ED.f"
	    i__1 = (*n - 1) * *k;
#line 354 "MB02ED.f"
	    dgemm_("NonTranspose", "NonTranspose", &i__1, nrhs, k, &c_b10, &t[
		    *k + 1 + t_dim1], ldt, &b[b_offset], ldb, &c_b19, &b[*k + 
		    1 + b_dim1], ldb, (ftnlen)12, (ftnlen)12);
#line 354 "MB02ED.f"
	}

#line 359 "MB02ED.f"
	i__1 = *n * *k;
#line 359 "MB02ED.f"
	dlaset_("All", k, k, &c_b21, &c_b10, &dwork[startn], &i__1, (ftnlen)3)
		;
#line 360 "MB02ED.f"
	i__1 = *n * *k;
#line 360 "MB02ED.f"
	dtrsm_("Right", "Lower", "Transpose", "NonUnit", k, k, &c_b10, &t[
		t_offset], ldt, &dwork[startn], &i__1, (ftnlen)5, (ftnlen)5, (
		ftnlen)9, (ftnlen)7);
#line 362 "MB02ED.f"
	if (*n > 1) {
#line 362 "MB02ED.f"
	    i__1 = (*n - 1) * *k;
#line 362 "MB02ED.f"
	    i__2 = *n * *k;
#line 362 "MB02ED.f"
	    dlacpy_("All", &i__1, k, &t[*k + 1 + t_dim1], ldt, &dwork[startn 
		    + *k], &i__2, (ftnlen)3);
#line 362 "MB02ED.f"
	}
#line 365 "MB02ED.f"
	i__1 = *n * *k;
#line 365 "MB02ED.f"
	dlacpy_("All", k, k, &dwork[startn], &i__1, &t[(*n - 1) * *k + 1 + 
		t_dim1], ldt, (ftnlen)3);

#line 368 "MB02ED.f"
	dtrmm_("Left", "Upper", "NonTranspose", "NonUnit", k, nrhs, &c_b10, &
		t[(*n - 1) * *k + 1 + t_dim1], ldt, &b[b_offset], ldb, (
		ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)7);

/*        Processing the generator. */

#line 373 "MB02ED.f"
	i__1 = *n;
#line 373 "MB02ED.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 374 "MB02ED.f"
	    startr = (i__ - 1) * *k + 1;
#line 375 "MB02ED.f"
	    starti = (*n - i__) * *k + 1;

/*           Transform the generator of T to proper form. */

#line 379 "MB02ED.f"
	    i__2 = *n * *k;
#line 379 "MB02ED.f"
	    i__3 = *k * 3;
#line 379 "MB02ED.f"
	    i__4 = *ldwork - starth + 1;
#line 379 "MB02ED.f"
	    mb02cx_("Column", k, k, k, &t[t_offset], ldt, &dwork[startn + (
		    i__ - 1) * *k], &i__2, &dwork[startt], &i__3, &dwork[
		    starth], &i__4, &ierr, (ftnlen)6);

#line 383 "MB02ED.f"
	    if (ierr != 0) {

/*              Error return:  The matrix is not positive definite. */

#line 387 "MB02ED.f"
		*info = 1;
#line 388 "MB02ED.f"
		return 0;
#line 389 "MB02ED.f"
	    }

/* Computing MAX */
#line 391 "MB02ED.f"
	    i__2 = maxwrk, i__3 = (integer) dwork[starth];
#line 391 "MB02ED.f"
	    maxwrk = max(i__2,i__3);
#line 392 "MB02ED.f"
	    i__2 = (*n - i__) * *k;
#line 392 "MB02ED.f"
	    i__3 = *n * *k;
#line 392 "MB02ED.f"
	    i__4 = *n * *k;
#line 392 "MB02ED.f"
	    i__5 = *k * 3;
#line 392 "MB02ED.f"
	    i__6 = *ldwork - starth + 1;
#line 392 "MB02ED.f"
	    mb02cy_("Column", "NoStructure", k, k, &i__2, k, &t[*k + 1 + 
		    t_dim1], ldt, &dwork[startn + i__ * *k], &i__3, &dwork[
		    startn + (i__ - 1) * *k], &i__4, &dwork[startt], &i__5, &
		    dwork[starth], &i__6, &ierr, (ftnlen)6, (ftnlen)11);
/* Computing MAX */
#line 396 "MB02ED.f"
	    i__2 = maxwrk, i__3 = (integer) dwork[starth];
#line 396 "MB02ED.f"
	    maxwrk = max(i__2,i__3);

/*           Block Gaussian eliminates the i-th block in B. */

#line 400 "MB02ED.f"
	    dtrsm_("Left", "Lower", "NonTranspose", "NonUnit", k, nrhs, &
		    c_b19, &t[t_offset], ldt, &b[startr + b_dim1], ldb, (
		    ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)7);
#line 402 "MB02ED.f"
	    if (*n > i__) {
#line 402 "MB02ED.f"
		i__2 = (*n - i__) * *k;
#line 402 "MB02ED.f"
		dgemm_("NonTranspose", "NonTranspose", &i__2, nrhs, k, &c_b10,
			 &t[*k + 1 + t_dim1], ldt, &b[startr + b_dim1], ldb, &
			c_b10, &b[startr + *k + b_dim1], ldb, (ftnlen)12, (
			ftnlen)12);
#line 402 "MB02ED.f"
	    }

/*           Apply hyperbolic transformations on the negative generator. */

#line 409 "MB02ED.f"
	    dlaset_("All", k, k, &c_b21, &c_b21, &t[starti + t_dim1], ldt, (
		    ftnlen)3);
#line 410 "MB02ED.f"
	    i__2 = (i__ - 1) * *k;
#line 410 "MB02ED.f"
	    i__3 = *n * *k;
#line 410 "MB02ED.f"
	    i__4 = *n * *k;
#line 410 "MB02ED.f"
	    i__5 = *k * 3;
#line 410 "MB02ED.f"
	    i__6 = *ldwork - starth + 1;
#line 410 "MB02ED.f"
	    mb02cy_("Column", "NoStructure", k, k, &i__2, k, &t[starti + 
		    t_dim1], ldt, &dwork[startn], &i__3, &dwork[startn + (i__ 
		    - 1) * *k], &i__4, &dwork[startt], &i__5, &dwork[starth], 
		    &i__6, &ierr, (ftnlen)6, (ftnlen)11);
/* Computing MAX */
#line 414 "MB02ED.f"
	    i__2 = maxwrk, i__3 = (integer) dwork[starth];
#line 414 "MB02ED.f"
	    maxwrk = max(i__2,i__3);

/*           Note that  DWORK(STARTN+(I-1)*K)  serves simultaneously */
/*           as the transformation container as well as the new block in */
/*           the negative generator. */

#line 420 "MB02ED.f"
	    i__2 = *n * *k;
#line 420 "MB02ED.f"
	    i__3 = *n * *k;
#line 420 "MB02ED.f"
	    i__4 = *k * 3;
#line 420 "MB02ED.f"
	    i__5 = *ldwork - starth + 1;
#line 420 "MB02ED.f"
	    mb02cy_("Column", "Triangular", k, k, k, k, &t[(*n - 1) * *k + 1 
		    + t_dim1], ldt, &dwork[startn + (i__ - 1) * *k], &i__2, &
		    dwork[startn + (i__ - 1) * *k], &i__3, &dwork[startt], &
		    i__4, &dwork[starth], &i__5, &ierr, (ftnlen)6, (ftnlen)10)
		    ;
/* Computing MAX */
#line 424 "MB02ED.f"
	    i__2 = maxwrk, i__3 = (integer) dwork[starth];
#line 424 "MB02ED.f"
	    maxwrk = max(i__2,i__3);

/*           Finally the Gaussian elimination is applied on the inverse */
/*           generator. */

#line 429 "MB02ED.f"
	    i__2 = (i__ - 1) * *k;
#line 429 "MB02ED.f"
	    dgemm_("NonTranspose", "NonTranspose", &i__2, nrhs, k, &c_b10, &t[
		    starti + t_dim1], ldt, &b[startr + b_dim1], ldb, &c_b10, &
		    b[b_offset], ldb, (ftnlen)12, (ftnlen)12);
#line 432 "MB02ED.f"
	    dtrmm_("Left", "Upper", "NonTranspose", "NonUnit", k, nrhs, &
		    c_b10, &t[(*n - 1) * *k + 1 + t_dim1], ldt, &b[startr + 
		    b_dim1], ldb, (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)7)
		    ;

#line 436 "MB02ED.f"
/* L20: */
#line 436 "MB02ED.f"
	}

#line 438 "MB02ED.f"
    }

/* Computing MAX */
#line 440 "MB02ED.f"
    i__1 = 1, i__2 = starth - 1 + maxwrk;
#line 440 "MB02ED.f"
    dwork[1] = (doublereal) max(i__1,i__2);

#line 442 "MB02ED.f"
    return 0;

/* *** Last line of MB02ED *** */
} /* mb02ed_ */
