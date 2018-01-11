#line 1 "MB02FD.f"
/* MB02FD.f -- translated by f2c (version 20100827).
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

#line 1 "MB02FD.f"
/* Table of constant values */

static doublereal c_b10 = 1.;

/* Subroutine */ int mb02fd_(char *typet, integer *k, integer *n, integer *p, 
	integer *s, doublereal *t, integer *ldt, doublereal *r__, integer *
	ldr, doublereal *dwork, integer *ldwork, integer *info, ftnlen 
	typet_len)
{
    /* System generated locals */
    integer r_dim1, r_offset, t_dim1, t_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, st, ierr;
    extern /* Subroutine */ int mb02cx_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen), mb02cy_(
	    char *, char *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical isrow;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dpotrf_(char *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);
    static integer maxwrk, countr, startr;


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

/*     To compute the incomplete Cholesky (ICC) factor of a symmetric */
/*     positive definite (s.p.d.) block Toeplitz matrix T, defined by */
/*     either its first block row, or its first block column, depending */
/*     on the routine parameter TYPET. */

/*     By subsequent calls of this routine, further rows / columns of */
/*     the Cholesky factor can be added. */
/*     Furthermore, the generator of the Schur complement of the leading */
/*     (P+S)*K-by-(P+S)*K block in T is available, which can be used, */
/*     e.g., for measuring the quality of the ICC factorization. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYPET   CHARACTER*1 */
/*             Specifies the type of T, as follows: */
/*             = 'R':  T contains the first block row of an s.p.d. block */
/*                     Toeplitz matrix; the ICC factor R is upper */
/*                     trapezoidal; */
/*             = 'C':  T contains the first block column of an s.p.d. */
/*                     block Toeplitz matrix; the ICC factor R is lower */
/*                     trapezoidal; this choice leads to better */
/*                     localized memory references and hence a faster */
/*                     algorithm. */
/*             Note:   in the sequel, the notation x / y means that */
/*                     x corresponds to TYPET = 'R' and y corresponds to */
/*                     TYPET = 'C'. */

/*     Input/Output Parameters */

/*     K       (input)  INTEGER */
/*             The number of rows / columns in T, which should be equal */
/*             to the blocksize.  K >= 0. */

/*     N       (input)  INTEGER */
/*             The number of blocks in T.  N >= 0. */

/*     P       (input)  INTEGER */
/*             The number of previously computed block rows / columns */
/*             of R.  0 <= P <= N. */

/*     S       (input)  INTEGER */
/*             The number of block rows / columns of R to compute. */
/*             0 <= S <= N-P. */

/*     T       (input/output)  DOUBLE PRECISION array, dimension */
/*             (LDT,(N-P)*K) / (LDT,K) */
/*             On entry, if P = 0, then the leading K-by-N*K / N*K-by-K */
/*             part of this array must contain the first block row / */
/*             column of an s.p.d. block Toeplitz matrix. */
/*             If P > 0, the leading K-by-(N-P)*K / (N-P)*K-by-K must */
/*             contain the negative generator of the Schur complement of */
/*             the leading P*K-by-P*K part in T, computed from previous */
/*             calls of this routine. */
/*             On exit, if INFO = 0, then the leading K-by-(N-P)*K / */
/*             (N-P)*K-by-K part of this array contains, in the first */
/*             K-by-K block, the upper / lower Cholesky factor of */
/*             T(1:K,1:K), in the following S-1 K-by-K blocks, the */
/*             Householder transformations applied during the process, */
/*             and in the remaining part, the negative generator of the */
/*             Schur complement of the leading (P+S)*K-by(P+S)*K part */
/*             in T. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T. */
/*             LDT >= MAX(1,K),        if TYPET = 'R'; */
/*             LDT >= MAX(1,(N-P)*K),  if TYPET = 'C'. */

/*     R       (input/output)  DOUBLE PRECISION array, dimension */
/*             (LDR, N*K)       / (LDR, S*K )     if P = 0; */
/*             (LDR, (N-P+1)*K) / (LDR, (S+1)*K ) if P > 0. */
/*             On entry, if P > 0, then the leading K-by-(N-P+1)*K / */
/*             (N-P+1)*K-by-K part of this array must contain the */
/*             nonzero blocks of the last block row / column in the */
/*             ICC factor from a previous call of this routine. Note that */
/*             this part is identical with the positive generator of */
/*             the Schur complement of the leading P*K-by-P*K part in T. */
/*             If P = 0, then R is only an output parameter. */
/*             On exit, if INFO = 0 and P = 0, then the leading */
/*             S*K-by-N*K / N*K-by-S*K part of this array contains the */
/*             upper / lower trapezoidal ICC factor. */
/*             On exit, if INFO = 0 and P > 0, then the leading */
/*             (S+1)*K-by-(N-P+1)*K / (N-P+1)*K-by-(S+1)*K part of this */
/*             array contains the upper / lower trapezoidal part of the */
/*             P-th to (P+S)-th block rows / columns of the ICC factor. */
/*             The elements in the strictly lower / upper trapezoidal */
/*             part are not referenced. */

/*     LDR     INTEGER */
/*             The leading dimension of the array R. */
/*             LDR >= MAX(1, S*K ),        if TYPET = 'R' and P = 0; */
/*             LDR >= MAX(1, (S+1)*K ),    if TYPET = 'R' and P > 0; */
/*             LDR >= MAX(1, N*K ),        if TYPET = 'C' and P = 0; */
/*             LDR >= MAX(1, (N-P+1)*K ),  if TYPET = 'C' and P > 0. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -11,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,(N+1)*K,4*K),   if P = 0; */
/*             LDWORK >= MAX(1,(N-P+2)*K,4*K), if P > 0. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the reduction algorithm failed; the Toeplitz matrix */
/*                   associated with T is not (numerically) positive */
/*                   definite in its leading (P+S)*K-by-(P+S)*K part. */

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
/*                               3 */
/*     The algorithm requires 0(K S (N-P)) floating point operations. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, April 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2001, */
/*     Mar. 2004. */

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

#line 205 "MB02FD.f"
    /* Parameter adjustments */
#line 205 "MB02FD.f"
    t_dim1 = *ldt;
#line 205 "MB02FD.f"
    t_offset = 1 + t_dim1;
#line 205 "MB02FD.f"
    t -= t_offset;
#line 205 "MB02FD.f"
    r_dim1 = *ldr;
#line 205 "MB02FD.f"
    r_offset = 1 + r_dim1;
#line 205 "MB02FD.f"
    r__ -= r_offset;
#line 205 "MB02FD.f"
    --dwork;
#line 205 "MB02FD.f"

#line 205 "MB02FD.f"
    /* Function Body */
#line 205 "MB02FD.f"
    *info = 0;
#line 206 "MB02FD.f"
    isrow = lsame_(typet, "R", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

#line 210 "MB02FD.f"
    if (! (isrow || lsame_(typet, "C", (ftnlen)1, (ftnlen)1))) {
#line 211 "MB02FD.f"
	*info = -1;
#line 212 "MB02FD.f"
    } else if (*k < 0) {
#line 213 "MB02FD.f"
	*info = -2;
#line 214 "MB02FD.f"
    } else if (*n < 0) {
#line 215 "MB02FD.f"
	*info = -3;
#line 216 "MB02FD.f"
    } else if (*p < 0 || *p > *n) {
#line 217 "MB02FD.f"
	*info = -4;
#line 218 "MB02FD.f"
    } else if (*s < 0 || *s > *n - *p) {
#line 219 "MB02FD.f"
	*info = -5;
#line 220 "MB02FD.f"
    } else if (*ldt < 1 || isrow && *ldt < *k || ! isrow && *ldt < (*n - *p) *
	     *k) {
#line 222 "MB02FD.f"
	*info = -7;
#line 223 "MB02FD.f"
    } else if (*ldr < 1 || isrow && *p == 0 && *ldr < *s * *k || isrow && *p 
	    > 0 && *ldr < (*s + 1) * *k || ! isrow && *p == 0 && *ldr < *n * *
	    k || ! isrow && *p > 0 && *ldr < (*n - *p + 1) * *k) {
#line 228 "MB02FD.f"
	*info = -9;
#line 229 "MB02FD.f"
    } else {
#line 230 "MB02FD.f"
	if (*p == 0) {
#line 231 "MB02FD.f"
	    countr = (*n + 1) * *k;
#line 232 "MB02FD.f"
	} else {
#line 233 "MB02FD.f"
	    countr = (*n - *p + 2) * *k;
#line 234 "MB02FD.f"
	}
/* Computing MAX */
#line 235 "MB02FD.f"
	i__1 = countr, i__2 = *k << 2;
#line 235 "MB02FD.f"
	countr = max(i__1,i__2);
#line 236 "MB02FD.f"
	if (*ldwork < max(1,countr)) {
#line 237 "MB02FD.f"
	    dwork[1] = (doublereal) max(1,countr);
#line 238 "MB02FD.f"
	    *info = -11;
#line 239 "MB02FD.f"
	}
#line 240 "MB02FD.f"
    }

/*     Return if there were illegal values. */

#line 244 "MB02FD.f"
    if (*info != 0) {
#line 245 "MB02FD.f"
	i__1 = -(*info);
#line 245 "MB02FD.f"
	xerbla_("MB02FD", &i__1, (ftnlen)6);
#line 246 "MB02FD.f"
	return 0;
#line 247 "MB02FD.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 251 "MB02FD.f"
    i__1 = min(*k,*n);
#line 251 "MB02FD.f"
    if (min(i__1,*s) == 0) {
#line 252 "MB02FD.f"
	dwork[1] = 1.;
#line 253 "MB02FD.f"
	return 0;
#line 254 "MB02FD.f"
    }

#line 256 "MB02FD.f"
    maxwrk = 1;

#line 258 "MB02FD.f"
    if (isrow) {

#line 260 "MB02FD.f"
	if (*p == 0) {

/*           T is the first block row of a block Toeplitz matrix. */
/*           Bring T to proper form by triangularizing its first block. */

#line 265 "MB02FD.f"
	    dpotrf_("Upper", k, &t[t_offset], ldt, &ierr, (ftnlen)5);
#line 266 "MB02FD.f"
	    if (ierr != 0) {

/*              Error return:  The matrix is not positive definite. */

#line 270 "MB02FD.f"
		*info = 1;
#line 271 "MB02FD.f"
		return 0;
#line 272 "MB02FD.f"
	    }

#line 274 "MB02FD.f"
	    if (*n > 1) {
#line 274 "MB02FD.f"
		i__1 = (*n - 1) * *k;
#line 274 "MB02FD.f"
		dtrsm_("Left", "Upper", "Transpose", "NonUnit", k, &i__1, &
			c_b10, &t[t_offset], ldt, &t[(*k + 1) * t_dim1 + 1], 
			ldt, (ftnlen)4, (ftnlen)5, (ftnlen)9, (ftnlen)7);
#line 274 "MB02FD.f"
	    }
#line 277 "MB02FD.f"
	    i__1 = *n * *k;
#line 277 "MB02FD.f"
	    dlacpy_("Upper", k, &i__1, &t[t_offset], ldt, &r__[r_offset], ldr,
		     (ftnlen)5);

#line 279 "MB02FD.f"
	    if (*s == 1) {
#line 280 "MB02FD.f"
		dwork[1] = 1.;
#line 281 "MB02FD.f"
		return 0;
#line 282 "MB02FD.f"
	    }

#line 284 "MB02FD.f"
	    st = 2;
#line 285 "MB02FD.f"
	    countr = (*n - 1) * *k;
#line 286 "MB02FD.f"
	} else {
#line 287 "MB02FD.f"
	    st = 1;
#line 288 "MB02FD.f"
	    countr = (*n - *p) * *k;
#line 289 "MB02FD.f"
	}

#line 291 "MB02FD.f"
	startr = 1;

#line 293 "MB02FD.f"
	i__1 = *s;
#line 293 "MB02FD.f"
	for (i__ = st; i__ <= i__1; ++i__) {
#line 294 "MB02FD.f"
	    dlacpy_("Upper", k, &countr, &r__[startr + startr * r_dim1], ldr, 
		    &r__[startr + *k + (startr + *k) * r_dim1], ldr, (ftnlen)
		    5);
#line 296 "MB02FD.f"
	    startr += *k;
#line 297 "MB02FD.f"
	    countr -= *k;
#line 298 "MB02FD.f"
	    i__2 = *k * 3;
#line 298 "MB02FD.f"
	    i__3 = *ldwork - *k * 3;
#line 298 "MB02FD.f"
	    mb02cx_("Row", k, k, k, &r__[startr + startr * r_dim1], ldr, &t[
		    startr * t_dim1 + 1], ldt, &dwork[1], &i__2, &dwork[*k * 
		    3 + 1], &i__3, &ierr, (ftnlen)3);
#line 301 "MB02FD.f"
	    if (ierr != 0) {

/*              Error return:  The matrix is not positive definite. */

#line 305 "MB02FD.f"
		*info = 1;
#line 306 "MB02FD.f"
		return 0;
#line 307 "MB02FD.f"
	    }

/* Computing MAX */
#line 309 "MB02FD.f"
	    i__2 = maxwrk, i__3 = (integer) dwork[*k * 3 + 1] + *k * 3;
#line 309 "MB02FD.f"
	    maxwrk = max(i__2,i__3);
#line 310 "MB02FD.f"
	    i__2 = *k * 3;
#line 310 "MB02FD.f"
	    i__3 = *ldwork - *k * 3;
#line 310 "MB02FD.f"
	    mb02cy_("Row", "NoStructure", k, k, &countr, k, &r__[startr + (
		    startr + *k) * r_dim1], ldr, &t[(startr + *k) * t_dim1 + 
		    1], ldt, &t[startr * t_dim1 + 1], ldt, &dwork[1], &i__2, &
		    dwork[*k * 3 + 1], &i__3, &ierr, (ftnlen)3, (ftnlen)11);
/* Computing MAX */
#line 314 "MB02FD.f"
	    i__2 = maxwrk, i__3 = (integer) dwork[*k * 3 + 1] + *k * 3;
#line 314 "MB02FD.f"
	    maxwrk = max(i__2,i__3);
#line 315 "MB02FD.f"
/* L10: */
#line 315 "MB02FD.f"
	}

#line 317 "MB02FD.f"
    } else {

#line 319 "MB02FD.f"
	if (*p == 0) {

/*           T is the first block column of a block Toeplitz matrix. */
/*           Bring T to proper form by triangularizing its first block. */

#line 324 "MB02FD.f"
	    dpotrf_("Lower", k, &t[t_offset], ldt, &ierr, (ftnlen)5);
#line 325 "MB02FD.f"
	    if (ierr != 0) {

/*              Error return:  The matrix is not positive definite. */

#line 329 "MB02FD.f"
		*info = 1;
#line 330 "MB02FD.f"
		return 0;
#line 331 "MB02FD.f"
	    }

#line 333 "MB02FD.f"
	    if (*n > 1) {
#line 333 "MB02FD.f"
		i__1 = (*n - 1) * *k;
#line 333 "MB02FD.f"
		dtrsm_("Right", "Lower", "Transpose", "NonUnit", &i__1, k, &
			c_b10, &t[t_offset], ldt, &t[*k + 1 + t_dim1], ldt, (
			ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)7);
#line 333 "MB02FD.f"
	    }
#line 336 "MB02FD.f"
	    i__1 = *n * *k;
#line 336 "MB02FD.f"
	    dlacpy_("Lower", &i__1, k, &t[t_offset], ldt, &r__[r_offset], ldr,
		     (ftnlen)5);

#line 338 "MB02FD.f"
	    if (*s == 1) {
#line 339 "MB02FD.f"
		dwork[1] = 1.;
#line 340 "MB02FD.f"
		return 0;
#line 341 "MB02FD.f"
	    }

#line 343 "MB02FD.f"
	    st = 2;
#line 344 "MB02FD.f"
	    countr = (*n - 1) * *k;
#line 345 "MB02FD.f"
	} else {
#line 346 "MB02FD.f"
	    st = 1;
#line 347 "MB02FD.f"
	    countr = (*n - *p) * *k;
#line 348 "MB02FD.f"
	}

#line 350 "MB02FD.f"
	startr = 1;

#line 352 "MB02FD.f"
	i__1 = *s;
#line 352 "MB02FD.f"
	for (i__ = st; i__ <= i__1; ++i__) {
#line 353 "MB02FD.f"
	    dlacpy_("Lower", &countr, k, &r__[startr + startr * r_dim1], ldr, 
		    &r__[startr + *k + (startr + *k) * r_dim1], ldr, (ftnlen)
		    5);
#line 355 "MB02FD.f"
	    startr += *k;
#line 356 "MB02FD.f"
	    countr -= *k;
#line 357 "MB02FD.f"
	    i__2 = *k * 3;
#line 357 "MB02FD.f"
	    i__3 = *ldwork - *k * 3;
#line 357 "MB02FD.f"
	    mb02cx_("Column", k, k, k, &r__[startr + startr * r_dim1], ldr, &
		    t[startr + t_dim1], ldt, &dwork[1], &i__2, &dwork[*k * 3 
		    + 1], &i__3, &ierr, (ftnlen)6);
#line 360 "MB02FD.f"
	    if (ierr != 0) {

/*              Error return:  The matrix is not positive definite. */

#line 364 "MB02FD.f"
		*info = 1;
#line 365 "MB02FD.f"
		return 0;
#line 366 "MB02FD.f"
	    }

/* Computing MAX */
#line 368 "MB02FD.f"
	    i__2 = maxwrk, i__3 = (integer) dwork[*k * 3 + 1] + *k * 3;
#line 368 "MB02FD.f"
	    maxwrk = max(i__2,i__3);
#line 369 "MB02FD.f"
	    i__2 = *k * 3;
#line 369 "MB02FD.f"
	    i__3 = *ldwork - *k * 3;
#line 369 "MB02FD.f"
	    mb02cy_("Column", "NoStructure", k, k, &countr, k, &r__[startr + *
		    k + startr * r_dim1], ldr, &t[startr + *k + t_dim1], ldt, 
		    &t[startr + t_dim1], ldt, &dwork[1], &i__2, &dwork[*k * 3 
		    + 1], &i__3, &ierr, (ftnlen)6, (ftnlen)11);
/* Computing MAX */
#line 373 "MB02FD.f"
	    i__2 = maxwrk, i__3 = (integer) dwork[*k * 3 + 1] + *k * 3;
#line 373 "MB02FD.f"
	    maxwrk = max(i__2,i__3);
#line 374 "MB02FD.f"
/* L20: */
#line 374 "MB02FD.f"
	}

#line 376 "MB02FD.f"
    }

#line 378 "MB02FD.f"
    dwork[1] = (doublereal) maxwrk;

#line 380 "MB02FD.f"
    return 0;

/* *** Last line of MB02FD *** */
} /* mb02fd_ */

