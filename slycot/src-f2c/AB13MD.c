#line 1 "AB13MD.f"
/* AB13MD.f -- translated by f2c (version 20100827).
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

#line 1 "AB13MD.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static doublereal c_b11 = 1.;
static doublereal c_b15 = 0.;
static doublereal c_b136 = 2.;

/* Subroutine */ int ab13md_(char *fact, integer *n, doublecomplex *z__, 
	integer *ldz, integer *m, integer *nblock, integer *itype, doublereal 
	*x, doublereal *bound, doublereal *d__, doublereal *g, integer *iwork,
	 doublereal *dwork, integer *ldwork, doublecomplex *zwork, integer *
	lzwork, integer *info, ftnlen fact_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double sqrt(doublereal);
    void d_cnjg(doublecomplex *, doublecomplex *);
    double pow_di(doublereal *, integer *), log(doublereal);

    /* Local variables */
    static doublereal c__, e;
    static integer i__, j, k, l;
    static doublereal t1, t2, t3, hn;
    static integer mr;
    static doublereal pp;
    static integer mt, iw2, iw3, iw4, iw5, iw6, iw7, iw8, iw9, iz2, iz3, iz4, 
	    iz5, iz6, iz7, iz8, iz9, iw10, iw11, iw12, iw13, iw14, iw15, iw16,
	     iw17, iw18, iw19, iw20, iw21, iw22, iw23, iw24, iw25, iw26, iw27,
	     iw28, iw29, iw30, iw31, iw32, iw33, iz10, iz11, iz12, iz13, iz14,
	     iz15, iz16, iz17, iz18, iz19, iz20, iz21, iz22, iz23, iz24, lwa, 
	    lza;
    static doublereal eps, phi, rat, tau, tol;
    static logical pos;
    static doublereal tol2, tol3, tol4, tol5;
    static doublecomplex detf;
    static doublereal emin, emax;
    static integer sdim;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer iter;
    static doublereal prod, temp;
    static integer iwrk, isum, nsum, info2;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal scale, delta;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal rcond;
    static logical xfact;
    extern /* Subroutine */ int zgees_(char *, char *, L_fp, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, logical *, integer *, ftnlen, ftnlen), dcopy_(
	    integer *, doublereal *, integer *, doublereal *, integer *), 
	    zgemm_(char *, char *, integer *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen, 
	    ftnlen);
    static doublereal hnorm, svlam;
    static logical bwork[1], gtest;
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static doublereal snorm;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static doublereal znorm;
    static integer izwrk;
    extern /* Subroutine */ int dsysv_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);
    static doublereal hnorm1, dlambd, ynorm1, ynorm2, znorm2;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    extern logical select_();
    static doublereal regpar;
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static doublecomplex tempij;
    static integer lwamax;
    static doublecomplex tempji;
    extern /* Subroutine */ int zlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, integer *, ftnlen);
    static integer lzamax;
    extern /* Subroutine */ int dsycon_(char *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, ftnlen), zgetrf_(integer *, integer *, 
	    doublecomplex *, integer *, integer *, integer *);
    static doublereal colsum;
    extern /* Subroutine */ int zgesvd_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, integer *, ftnlen, ftnlen), zgetri_(integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *), zlacpy_(char *, integer *, integer *, doublecomplex *
	    , integer *, doublecomplex *, integer *, ftnlen);
    static integer minwrk, minzrk;
    extern /* Subroutine */ int dsytrf_(char *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen);
    static doublereal stsize;
    extern /* Subroutine */ int dsytrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static doublereal rowsum;


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

/*     To compute an upper bound on the structured singular value for a */
/*     given square complex matrix and a given block structure of the */
/*     uncertainty. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     FACT    CHARACTER*1 */
/*             Specifies whether or not an information from the */
/*             previous call is supplied in the vector X. */
/*             = 'F':  On entry, X contains information from the */
/*                     previous call. */
/*             = 'N':  On entry, X does not contain an information from */
/*                     the previous call. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix Z.  N >= 0. */

/*     Z       (input) COMPLEX*16 array, dimension (LDZ,N) */
/*             The leading N-by-N part of this array must contain the */
/*             complex matrix Z for which the upper bound on the */
/*             structured singular value is to be computed. */

/*     LDZ     INTEGER */
/*             The leading dimension of the array Z.  LDZ >= max(1,N). */

/*     M       (input) INTEGER */
/*             The number of diagonal blocks in the block structure of */
/*             the uncertainty.  M >= 1. */

/*     NBLOCK  (input) INTEGER array, dimension (M) */
/*             The vector of length M containing the block structure */
/*             of the uncertainty. NBLOCK(I), I = 1:M, is the size of */
/*             each block. */

/*     ITYPE   (input) INTEGER array, dimension (M) */
/*             The vector of length M indicating the type of each block. */
/*             For I = 1:M, */
/*             ITYPE(I) = 1 indicates that the corresponding block is a */
/*                          real block, and */
/*             ITYPE(I) = 2 indicates that the corresponding block is a */
/*                          complex block. */
/*             NBLOCK(I) must be equal to 1 if ITYPE(I) is equal to 1. */

/*     X       (input/output) DOUBLE PRECISION array, dimension */
/*             ( M + MR - 1 ), where MR is the number of the real blocks. */
/*             On entry, if FACT = 'F' and NBLOCK(1) < N, this array */
/*             must contain information from the previous call to AB13MD. */
/*             If NBLOCK(1) = N, this array is not used. */
/*             On exit, if NBLOCK(1) < N, this array contains information */
/*             that can be used in the next call to AB13MD for a matrix */
/*             close to Z. */

/*     BOUND   (output) DOUBLE PRECISION */
/*             The upper bound on the structured singular value. */

/*     D, G    (output) DOUBLE PRECISION arrays, dimension (N) */
/*             The vectors of length N containing the diagonal entries */
/*             of the diagonal N-by-N matrices D and G, respectively, */
/*             such that the matrix */
/*             Z'*D^2*Z + sqrt(-1)*(G*Z-Z'*G) - BOUND^2*D^2 */
/*             is negative semidefinite. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension MAX(4*M-2,N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= 2*N*N*M - N*N + 9*M*M + N*M + 11*N + 33*M - 11. */
/*             For best performance */
/*             LDWORK >= 2*N*N*M - N*N + 9*M*M + N*M + 6*N + 33*M - 11 + */
/*                       MAX( 5*N,2*N*NB ) */
/*             where NB is the optimal blocksize returned by ILAENV. */

/*     ZWORK   COMPLEX*16 array, dimension (LZWORK) */
/*             On exit, if INFO = 0, ZWORK(1) contains the optimal value */
/*             of LZWORK. */

/*     LZWORK  INTEGER */
/*             The dimension of the array ZWORK. */
/*             LZWORK >= 6*N*N*M + 12*N*N + 6*M + 6*N - 3. */
/*             For best performance */
/*             LZWORK >= 6*N*N*M + 12*N*N + 6*M + 3*N - 3 + */
/*                       MAX( 3*N,N*NB ) */
/*             where NB is the optimal blocksize returned by ILAENV. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the block sizes must be positive integers; */
/*             = 2:  the sum of block sizes must be equal to N; */
/*             = 3:  the size of a real block must be equal to 1; */
/*             = 4:  the block type must be either 1 or 2; */
/*             = 5:  errors in solving linear equations or in matrix */
/*                   inversion; */
/*             = 6:  errors in computing eigenvalues or singular values. */

/*     METHOD */

/*     The routine computes the upper bound proposed in [1]. */

/*     REFERENCES */

/*     [1] Fan, M.K.H., Tits, A.L., and Doyle, J.C. */
/*         Robustness in the presence of mixed parametric uncertainty */
/*         and unmodeled dynamics. */
/*         IEEE Trans. Automatic Control, vol. AC-36, 1991, pp. 25-38. */

/*     NUMERICAL ASPECTS */

/*     The accuracy and speed of computation depend on the value of */
/*     the internal threshold TOL. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, F. Delebecque, D.W. Gu, M.M. Konstantinov and */
/*     S. Steer with the assistance of V. Sima, September 2000. */

/*     REVISIONS */

/*     V. Sima, Katholieke Universiteit Leuven, February 2001. */

/*     KEYWORDS */

/*     H-infinity optimal control, Robust control, Structured singular */
/*     value. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Compute workspace. */

#line 233 "AB13MD.f"
    /* Parameter adjustments */
#line 233 "AB13MD.f"
    z_dim1 = *ldz;
#line 233 "AB13MD.f"
    z_offset = 1 + z_dim1;
#line 233 "AB13MD.f"
    z__ -= z_offset;
#line 233 "AB13MD.f"
    --nblock;
#line 233 "AB13MD.f"
    --itype;
#line 233 "AB13MD.f"
    --x;
#line 233 "AB13MD.f"
    --d__;
#line 233 "AB13MD.f"
    --g;
#line 233 "AB13MD.f"
    --iwork;
#line 233 "AB13MD.f"
    --dwork;
#line 233 "AB13MD.f"
    --zwork;
#line 233 "AB13MD.f"

#line 233 "AB13MD.f"
    /* Function Body */
#line 233 "AB13MD.f"
    minwrk = (*n << 1) * *n * *m - *n * *n + *m * 9 * *m + *n * *m + *n * 11 
	    + *m * 33 - 11;
#line 234 "AB13MD.f"
    minzrk = *n * 6 * *n * *m + *n * 12 * *n + *m * 6 + *n * 6 - 3;

/*     Decode and Test input parameters. */

#line 238 "AB13MD.f"
    *info = 0;
#line 239 "AB13MD.f"
    xfact = lsame_(fact, "F", (ftnlen)1, (ftnlen)1);
#line 240 "AB13MD.f"
    if (! xfact && ! lsame_(fact, "N", (ftnlen)1, (ftnlen)1)) {
#line 241 "AB13MD.f"
	*info = -1;
#line 242 "AB13MD.f"
    } else if (*n < 0) {
#line 243 "AB13MD.f"
	*info = -2;
#line 244 "AB13MD.f"
    } else if (*ldz < max(1,*n)) {
#line 245 "AB13MD.f"
	*info = -4;
#line 246 "AB13MD.f"
    } else if (*m < 1) {
#line 247 "AB13MD.f"
	*info = -5;
#line 248 "AB13MD.f"
    } else if (*ldwork < minwrk) {
#line 249 "AB13MD.f"
	*info = -14;
#line 250 "AB13MD.f"
    } else if (*lzwork < minzrk) {
#line 251 "AB13MD.f"
	*info = -16;
#line 252 "AB13MD.f"
    }
#line 253 "AB13MD.f"
    if (*info != 0) {
#line 254 "AB13MD.f"
	i__1 = -(*info);
#line 254 "AB13MD.f"
	xerbla_("AB13MD", &i__1, (ftnlen)6);
#line 255 "AB13MD.f"
	return 0;
#line 256 "AB13MD.f"
    }

#line 258 "AB13MD.f"
    nsum = 0;
#line 259 "AB13MD.f"
    isum = 0;
#line 260 "AB13MD.f"
    mr = 0;
#line 261 "AB13MD.f"
    i__1 = *m;
#line 261 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 262 "AB13MD.f"
	if (nblock[i__] < 1) {
#line 263 "AB13MD.f"
	    *info = 1;
#line 264 "AB13MD.f"
	    return 0;
#line 265 "AB13MD.f"
	}
#line 266 "AB13MD.f"
	if (itype[i__] == 1 && nblock[i__] > 1) {
#line 267 "AB13MD.f"
	    *info = 3;
#line 268 "AB13MD.f"
	    return 0;
#line 269 "AB13MD.f"
	}
#line 270 "AB13MD.f"
	nsum += nblock[i__];
#line 271 "AB13MD.f"
	if (itype[i__] == 1) {
#line 271 "AB13MD.f"
	    ++mr;
#line 271 "AB13MD.f"
	}
#line 272 "AB13MD.f"
	if (itype[i__] == 1 || itype[i__] == 2) {
#line 272 "AB13MD.f"
	    ++isum;
#line 272 "AB13MD.f"
	}
#line 273 "AB13MD.f"
/* L10: */
#line 273 "AB13MD.f"
    }
#line 274 "AB13MD.f"
    if (nsum != *n) {
#line 275 "AB13MD.f"
	*info = 2;
#line 276 "AB13MD.f"
	return 0;
#line 277 "AB13MD.f"
    }
#line 278 "AB13MD.f"
    if (isum != *m) {
#line 279 "AB13MD.f"
	*info = 4;
#line 280 "AB13MD.f"
	return 0;
#line 281 "AB13MD.f"
    }
#line 282 "AB13MD.f"
    mt = *m + mr - 1;

#line 284 "AB13MD.f"
    lwamax = 0;
#line 285 "AB13MD.f"
    lzamax = 0;

/*     Set D = In, G = 0. */

#line 289 "AB13MD.f"
    dlaset_("Full", n, &c__1, &c_b11, &c_b11, &d__[1], n, (ftnlen)4);
#line 290 "AB13MD.f"
    dlaset_("Full", n, &c__1, &c_b15, &c_b15, &g[1], n, (ftnlen)4);

/*     Quick return if possible. */

#line 294 "AB13MD.f"
    znorm = zlange_("F", n, n, &z__[z_offset], ldz, &dwork[1], (ftnlen)1);
#line 295 "AB13MD.f"
    if (znorm == 0.) {
#line 296 "AB13MD.f"
	*bound = 0.;
#line 297 "AB13MD.f"
	dwork[1] = 1.;
#line 298 "AB13MD.f"
	zwork[1].r = 1., zwork[1].i = 0.;
#line 299 "AB13MD.f"
	return 0;
#line 300 "AB13MD.f"
    }

/*     Copy Z into ZWORK. */

#line 304 "AB13MD.f"
    zlacpy_("Full", n, n, &z__[z_offset], ldz, &zwork[1], n, (ftnlen)4);

/*     Exact bound for the case NBLOCK( 1 ) = N. */

#line 308 "AB13MD.f"
    if (nblock[1] == *n) {
#line 309 "AB13MD.f"
	if (itype[1] == 1) {

/*           1-by-1 real block. */

#line 313 "AB13MD.f"
	    *bound = 0.;
#line 314 "AB13MD.f"
	    dwork[1] = 1.;
#line 315 "AB13MD.f"
	    zwork[1].r = 1., zwork[1].i = 0.;
#line 316 "AB13MD.f"
	} else {

/*           N-by-N complex block. */

#line 320 "AB13MD.f"
	    zgesvd_("N", "N", n, n, &zwork[1], n, &dwork[1], &zwork[1], &c__1,
		     &zwork[1], &c__1, &zwork[*n * *n + 1], lzwork, &dwork[*n 
		    + 1], &info2, (ftnlen)1, (ftnlen)1);
#line 323 "AB13MD.f"
	    if (info2 > 0) {
#line 324 "AB13MD.f"
		*info = 6;
#line 325 "AB13MD.f"
		return 0;
#line 326 "AB13MD.f"
	    }
#line 327 "AB13MD.f"
	    *bound = dwork[1];
#line 328 "AB13MD.f"
	    i__1 = *n * *n + 1;
#line 328 "AB13MD.f"
	    lza = *n * *n + (integer) zwork[i__1].r;
#line 329 "AB13MD.f"
	    dwork[1] = (doublereal) (*n * 5);
#line 330 "AB13MD.f"
	    z__1.r = (doublereal) lza, z__1.i = 0.;
#line 330 "AB13MD.f"
	    zwork[1].r = z__1.r, zwork[1].i = z__1.i;
#line 331 "AB13MD.f"
	}
#line 332 "AB13MD.f"
	return 0;
#line 333 "AB13MD.f"
    }

/*     Get machine precision. */

#line 337 "AB13MD.f"
    eps = dlamch_("P", (ftnlen)1);

/*     Set tolerances. */

#line 341 "AB13MD.f"
    tol = sqrt(eps) * 100.;
#line 342 "AB13MD.f"
    tol2 = eps * 1e4;
#line 343 "AB13MD.f"
    tol3 = eps * 10.;
#line 344 "AB13MD.f"
    tol4 = .001;
#line 345 "AB13MD.f"
    tol5 = .001;
#line 346 "AB13MD.f"
    regpar = eps * 1e3;

/*     Real workspace usage. */

#line 350 "AB13MD.f"
    iw2 = *m * *m;
#line 351 "AB13MD.f"
    iw3 = iw2 + *m;
#line 352 "AB13MD.f"
    iw4 = iw3 + *n;
#line 353 "AB13MD.f"
    iw5 = iw4 + *m;
#line 354 "AB13MD.f"
    iw6 = iw5 + *m;
#line 355 "AB13MD.f"
    iw7 = iw6 + *n;
#line 356 "AB13MD.f"
    iw8 = iw7 + *n;
#line 357 "AB13MD.f"
    iw9 = iw8 + *n * (*m - 1);
#line 358 "AB13MD.f"
    iw10 = iw9 + *n * *n * mt;
#line 359 "AB13MD.f"
    iw11 = iw10 + mt;
#line 360 "AB13MD.f"
    iw12 = iw11 + mt * mt;
#line 361 "AB13MD.f"
    iw13 = iw12 + *n;
#line 362 "AB13MD.f"
    iw14 = iw13 + mt + 1;
#line 363 "AB13MD.f"
    iw15 = iw14 + mt + 1;
#line 364 "AB13MD.f"
    iw16 = iw15 + mt + 1;
#line 365 "AB13MD.f"
    iw17 = iw16 + mt + 1;
#line 366 "AB13MD.f"
    iw18 = iw17 + mt + 1;
#line 367 "AB13MD.f"
    iw19 = iw18 + mt;
#line 368 "AB13MD.f"
    iw20 = iw19 + mt;
#line 369 "AB13MD.f"
    iw21 = iw20 + mt;
#line 370 "AB13MD.f"
    iw22 = iw21 + *n;
#line 371 "AB13MD.f"
    iw23 = iw22 + *m - 1;
#line 372 "AB13MD.f"
    iw24 = iw23 + mr;
#line 373 "AB13MD.f"
    iw25 = iw24 + *n;
#line 374 "AB13MD.f"
    iw26 = iw25 + (mt << 1);
#line 375 "AB13MD.f"
    iw27 = iw26 + mt;
#line 376 "AB13MD.f"
    iw28 = iw27 + mt;
#line 377 "AB13MD.f"
    iw29 = iw28 + *m - 1;
#line 378 "AB13MD.f"
    iw30 = iw29 + mr;
#line 379 "AB13MD.f"
    iw31 = iw30 + *n + (mt << 1);
#line 380 "AB13MD.f"
    iw32 = iw31 + mt * mt;
#line 381 "AB13MD.f"
    iw33 = iw32 + mt;
#line 382 "AB13MD.f"
    iwrk = iw33 + mt + 1;

/*     Double complex workspace usage. */

#line 386 "AB13MD.f"
    iz2 = *n * *n;
#line 387 "AB13MD.f"
    iz3 = iz2 + *n * *n;
#line 388 "AB13MD.f"
    iz4 = iz3 + *n * *n;
#line 389 "AB13MD.f"
    iz5 = iz4 + *n * *n;
#line 390 "AB13MD.f"
    iz6 = iz5 + *n * *n;
#line 391 "AB13MD.f"
    iz7 = iz6 + *n * *n * mt;
#line 392 "AB13MD.f"
    iz8 = iz7 + *n * *n;
#line 393 "AB13MD.f"
    iz9 = iz8 + *n * *n;
#line 394 "AB13MD.f"
    iz10 = iz9 + *n * *n;
#line 395 "AB13MD.f"
    iz11 = iz10 + mt;
#line 396 "AB13MD.f"
    iz12 = iz11 + *n * *n;
#line 397 "AB13MD.f"
    iz13 = iz12 + *n;
#line 398 "AB13MD.f"
    iz14 = iz13 + *n * *n;
#line 399 "AB13MD.f"
    iz15 = iz14 + *n;
#line 400 "AB13MD.f"
    iz16 = iz15 + *n * *n;
#line 401 "AB13MD.f"
    iz17 = iz16 + *n;
#line 402 "AB13MD.f"
    iz18 = iz17 + *n * *n;
#line 403 "AB13MD.f"
    iz19 = iz18 + *n * *n * mt;
#line 404 "AB13MD.f"
    iz20 = iz19 + mt;
#line 405 "AB13MD.f"
    iz21 = iz20 + *n * *n * mt;
#line 406 "AB13MD.f"
    iz22 = iz21 + *n * *n;
#line 407 "AB13MD.f"
    iz23 = iz22 + *n * *n;
#line 408 "AB13MD.f"
    iz24 = iz23 + *n * *n;
#line 409 "AB13MD.f"
    izwrk = iz24 + mt;

/*     Compute the cumulative sums of blocks dimensions. */

#line 413 "AB13MD.f"
    iwork[1] = 0;
#line 414 "AB13MD.f"
    i__1 = *m + 1;
#line 414 "AB13MD.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 415 "AB13MD.f"
	iwork[i__] = iwork[i__ - 1] + nblock[i__ - 1];
#line 416 "AB13MD.f"
/* L20: */
#line 416 "AB13MD.f"
    }

/*     Find Osborne scaling if initial scaling is not given. */

#line 420 "AB13MD.f"
    if (! xfact) {
#line 421 "AB13MD.f"
	dlaset_("Full", m, m, &c_b15, &c_b15, &dwork[1], m, (ftnlen)4);
#line 422 "AB13MD.f"
	dlaset_("Full", m, &c__1, &c_b11, &c_b11, &dwork[iw2 + 1], m, (ftnlen)
		4);
#line 423 "AB13MD.f"
	znorm = zlange_("F", n, n, &zwork[1], n, &dwork[1], (ftnlen)1);
#line 424 "AB13MD.f"
	i__1 = *m;
#line 424 "AB13MD.f"
	for (j = 1; j <= i__1; ++j) {
#line 425 "AB13MD.f"
	    i__2 = *m;
#line 425 "AB13MD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 426 "AB13MD.f"
		if (i__ != j) {
#line 427 "AB13MD.f"
		    i__3 = iwork[i__ + 1] - iwork[i__];
#line 427 "AB13MD.f"
		    i__4 = iwork[j + 1] - iwork[j];
#line 427 "AB13MD.f"
		    zlacpy_("Full", &i__3, &i__4, &z__[iwork[i__] + 1 + (
			    iwork[j] + 1) * z_dim1], ldz, &zwork[iz2 + 1], n, 
			    (ftnlen)4);
#line 431 "AB13MD.f"
		    i__3 = iwork[i__ + 1] - iwork[i__];
#line 431 "AB13MD.f"
		    i__4 = iwork[j + 1] - iwork[j];
#line 431 "AB13MD.f"
		    i__5 = *lzwork - izwrk;
#line 431 "AB13MD.f"
		    zgesvd_("N", "N", &i__3, &i__4, &zwork[iz2 + 1], n, &
			    dwork[iw3 + 1], &zwork[1], &c__1, &zwork[1], &
			    c__1, &zwork[izwrk + 1], &i__5, &dwork[iwrk + 1], 
			    &info2, (ftnlen)1, (ftnlen)1);
#line 436 "AB13MD.f"
		    if (info2 > 0) {
#line 437 "AB13MD.f"
			*info = 6;
#line 438 "AB13MD.f"
			return 0;
#line 439 "AB13MD.f"
		    }
#line 440 "AB13MD.f"
		    i__3 = izwrk + 1;
#line 440 "AB13MD.f"
		    lza = (integer) zwork[i__3].r;
#line 441 "AB13MD.f"
		    lzamax = max(lza,lzamax);
#line 442 "AB13MD.f"
		    znorm2 = dwork[iw3 + 1];
#line 443 "AB13MD.f"
		    dwork[i__ + (j - 1) * *m] = znorm2 + znorm * tol2;
#line 444 "AB13MD.f"
		}
#line 445 "AB13MD.f"
/* L30: */
#line 445 "AB13MD.f"
	    }
#line 446 "AB13MD.f"
/* L40: */
#line 446 "AB13MD.f"
	}
#line 447 "AB13MD.f"
	dlaset_("Full", m, &c__1, &c_b15, &c_b15, &dwork[iw4 + 1], m, (ftnlen)
		4);
#line 448 "AB13MD.f"
L50:
#line 448 "AB13MD.f"
	i__1 = *m;
#line 448 "AB13MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 449 "AB13MD.f"
	    dwork[iw5 + i__] = dwork[iw4 + i__] - 1.;
#line 450 "AB13MD.f"
/* L60: */
#line 450 "AB13MD.f"
	}
#line 451 "AB13MD.f"
	hnorm = dlange_("F", m, &c__1, &dwork[iw5 + 1], m, &dwork[1], (ftnlen)
		1);
#line 452 "AB13MD.f"
	if (hnorm <= tol2) {
#line 452 "AB13MD.f"
	    goto L120;
#line 452 "AB13MD.f"
	}
#line 453 "AB13MD.f"
	i__1 = *m;
#line 453 "AB13MD.f"
	for (k = 1; k <= i__1; ++k) {
#line 454 "AB13MD.f"
	    colsum = 0.;
#line 455 "AB13MD.f"
	    i__2 = *m;
#line 455 "AB13MD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 456 "AB13MD.f"
		colsum += dwork[i__ + (k - 1) * *m];
#line 457 "AB13MD.f"
/* L70: */
#line 457 "AB13MD.f"
	    }
#line 458 "AB13MD.f"
	    rowsum = 0.;
#line 459 "AB13MD.f"
	    i__2 = *m;
#line 459 "AB13MD.f"
	    for (j = 1; j <= i__2; ++j) {
#line 460 "AB13MD.f"
		rowsum += dwork[k + (j - 1) * *m];
#line 461 "AB13MD.f"
/* L80: */
#line 461 "AB13MD.f"
	    }
#line 462 "AB13MD.f"
	    rat = sqrt(colsum / rowsum);
#line 463 "AB13MD.f"
	    dwork[iw4 + k] = rat;
#line 464 "AB13MD.f"
	    i__2 = *m;
#line 464 "AB13MD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 465 "AB13MD.f"
		dwork[i__ + (k - 1) * *m] /= rat;
#line 466 "AB13MD.f"
/* L90: */
#line 466 "AB13MD.f"
	    }
#line 467 "AB13MD.f"
	    i__2 = *m;
#line 467 "AB13MD.f"
	    for (j = 1; j <= i__2; ++j) {
#line 468 "AB13MD.f"
		dwork[k + (j - 1) * *m] *= rat;
#line 469 "AB13MD.f"
/* L100: */
#line 469 "AB13MD.f"
	    }
#line 470 "AB13MD.f"
	    dwork[iw2 + k] *= rat;
#line 471 "AB13MD.f"
/* L110: */
#line 471 "AB13MD.f"
	}
#line 472 "AB13MD.f"
	goto L50;
#line 473 "AB13MD.f"
L120:
#line 473 "AB13MD.f"
	scale = 1. / dwork[iw2 + 1];
#line 474 "AB13MD.f"
	dscal_(m, &scale, &dwork[iw2 + 1], &c__1);
#line 475 "AB13MD.f"
    } else {
#line 476 "AB13MD.f"
	dwork[iw2 + 1] = 1.;
#line 477 "AB13MD.f"
	i__1 = *m;
#line 477 "AB13MD.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 478 "AB13MD.f"
	    dwork[iw2 + i__] = sqrt(x[i__ - 1]);
#line 479 "AB13MD.f"
/* L130: */
#line 479 "AB13MD.f"
	}
#line 480 "AB13MD.f"
    }
#line 481 "AB13MD.f"
    i__1 = *m;
#line 481 "AB13MD.f"
    for (j = 1; j <= i__1; ++j) {
#line 482 "AB13MD.f"
	i__2 = *m;
#line 482 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 483 "AB13MD.f"
	    if (i__ != j) {
#line 484 "AB13MD.f"
		i__3 = iwork[i__ + 1] - iwork[i__];
#line 484 "AB13MD.f"
		i__4 = iwork[j + 1] - iwork[j];
#line 484 "AB13MD.f"
		zlascl_("G", m, m, &dwork[iw2 + j], &dwork[iw2 + i__], &i__3, 
			&i__4, &zwork[iwork[i__] + 1 + iwork[j] * *n], n, &
			info2, (ftnlen)1);
#line 489 "AB13MD.f"
	    }
#line 490 "AB13MD.f"
/* L140: */
#line 490 "AB13MD.f"
	}
#line 491 "AB13MD.f"
/* L150: */
#line 491 "AB13MD.f"
    }

/*     Scale Z by its 2-norm. */

#line 495 "AB13MD.f"
    zlacpy_("Full", n, n, &zwork[1], n, &zwork[iz2 + 1], n, (ftnlen)4);
#line 496 "AB13MD.f"
    i__1 = *lzwork - izwrk;
#line 496 "AB13MD.f"
    zgesvd_("N", "N", n, n, &zwork[iz2 + 1], n, &dwork[iw3 + 1], &zwork[1], &
	    c__1, &zwork[1], &c__1, &zwork[izwrk + 1], &i__1, &dwork[iwrk + 1]
	    , &info2, (ftnlen)1, (ftnlen)1);
#line 499 "AB13MD.f"
    if (info2 > 0) {
#line 500 "AB13MD.f"
	*info = 6;
#line 501 "AB13MD.f"
	return 0;
#line 502 "AB13MD.f"
    }
#line 503 "AB13MD.f"
    i__1 = izwrk + 1;
#line 503 "AB13MD.f"
    lza = (integer) zwork[i__1].r;
#line 504 "AB13MD.f"
    lzamax = max(lza,lzamax);
#line 505 "AB13MD.f"
    znorm = dwork[iw3 + 1];
#line 506 "AB13MD.f"
    zlascl_("G", m, m, &znorm, &c_b11, n, n, &zwork[1], n, &info2, (ftnlen)1);

/*     Set BB. */

#line 510 "AB13MD.f"
    i__1 = *n * *n;
#line 510 "AB13MD.f"
    i__2 = *n * *n;
#line 510 "AB13MD.f"
    dlaset_("Full", &i__1, &mt, &c_b15, &c_b15, &dwork[iw9 + 1], &i__2, (
	    ftnlen)4);

/*     Set P. */

#line 514 "AB13MD.f"
    i__1 = nblock[1];
#line 514 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 515 "AB13MD.f"
	dwork[iw6 + i__] = 1.;
#line 516 "AB13MD.f"
/* L160: */
#line 516 "AB13MD.f"
    }
#line 517 "AB13MD.f"
    i__1 = *n;
#line 517 "AB13MD.f"
    for (i__ = nblock[1] + 1; i__ <= i__1; ++i__) {
#line 518 "AB13MD.f"
	dwork[iw6 + i__] = 0.;
#line 519 "AB13MD.f"
/* L170: */
#line 519 "AB13MD.f"
    }

/*     Compute P*Z. */

#line 523 "AB13MD.f"
    i__1 = *n;
#line 523 "AB13MD.f"
    for (j = 1; j <= i__1; ++j) {
#line 524 "AB13MD.f"
	i__2 = *n;
#line 524 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 525 "AB13MD.f"
	    i__3 = iz3 + i__ + (j - 1) * *n;
#line 525 "AB13MD.f"
	    i__4 = iw6 + i__;
#line 525 "AB13MD.f"
	    z__2.r = dwork[i__4], z__2.i = 0.;
#line 525 "AB13MD.f"
	    i__5 = i__ + (j - 1) * *n;
#line 525 "AB13MD.f"
	    z__1.r = z__2.r * zwork[i__5].r - z__2.i * zwork[i__5].i, z__1.i =
		     z__2.r * zwork[i__5].i + z__2.i * zwork[i__5].r;
#line 525 "AB13MD.f"
	    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 527 "AB13MD.f"
/* L180: */
#line 527 "AB13MD.f"
	}
#line 528 "AB13MD.f"
/* L190: */
#line 528 "AB13MD.f"
    }

/*     Compute Z'*P*Z. */

#line 532 "AB13MD.f"
    zgemm_("C", "N", n, n, n, &c_b2, &zwork[1], n, &zwork[iz3 + 1], n, &c_b1, 
	    &zwork[iz4 + 1], n, (ftnlen)1, (ftnlen)1);

/*     Copy Z'*P*Z into A0. */

#line 537 "AB13MD.f"
    zlacpy_("Full", n, n, &zwork[iz4 + 1], n, &zwork[iz5 + 1], n, (ftnlen)4);

/*     Copy diag(P) into B0d. */

#line 541 "AB13MD.f"
    dcopy_(n, &dwork[iw6 + 1], &c__1, &dwork[iw7 + 1], &c__1);

#line 543 "AB13MD.f"
    i__1 = *m;
#line 543 "AB13MD.f"
    for (k = 2; k <= i__1; ++k) {

/*        Set P. */

#line 547 "AB13MD.f"
	i__2 = iwork[k];
#line 547 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 548 "AB13MD.f"
	    dwork[iw6 + i__] = 0.;
#line 549 "AB13MD.f"
/* L200: */
#line 549 "AB13MD.f"
	}
#line 550 "AB13MD.f"
	i__2 = iwork[k] + nblock[k];
#line 550 "AB13MD.f"
	for (i__ = iwork[k] + 1; i__ <= i__2; ++i__) {
#line 551 "AB13MD.f"
	    dwork[iw6 + i__] = 1.;
#line 552 "AB13MD.f"
/* L210: */
#line 552 "AB13MD.f"
	}
#line 553 "AB13MD.f"
	if (k < *m) {
#line 554 "AB13MD.f"
	    i__2 = *n;
#line 554 "AB13MD.f"
	    for (i__ = iwork[k + 1] + 1; i__ <= i__2; ++i__) {
#line 555 "AB13MD.f"
		dwork[iw6 + i__] = 0.;
#line 556 "AB13MD.f"
/* L220: */
#line 556 "AB13MD.f"
	    }
#line 557 "AB13MD.f"
	}

/*        Compute P*Z. */

#line 561 "AB13MD.f"
	i__2 = *n;
#line 561 "AB13MD.f"
	for (j = 1; j <= i__2; ++j) {
#line 562 "AB13MD.f"
	    i__3 = *n;
#line 562 "AB13MD.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 563 "AB13MD.f"
		i__4 = iz3 + i__ + (j - 1) * *n;
#line 563 "AB13MD.f"
		i__5 = iw6 + i__;
#line 563 "AB13MD.f"
		z__2.r = dwork[i__5], z__2.i = 0.;
#line 563 "AB13MD.f"
		i__6 = i__ + (j - 1) * *n;
#line 563 "AB13MD.f"
		z__1.r = z__2.r * zwork[i__6].r - z__2.i * zwork[i__6].i, 
			z__1.i = z__2.r * zwork[i__6].i + z__2.i * zwork[i__6]
			.r;
#line 563 "AB13MD.f"
		zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
#line 565 "AB13MD.f"
/* L230: */
#line 565 "AB13MD.f"
	    }
#line 566 "AB13MD.f"
/* L240: */
#line 566 "AB13MD.f"
	}

/*        Compute t = Z'*P*Z. */

#line 570 "AB13MD.f"
	zgemm_("C", "N", n, n, n, &c_b2, &zwork[1], n, &zwork[iz3 + 1], n, &
		c_b1, &zwork[iz4 + 1], n, (ftnlen)1, (ftnlen)1);

/*        Copy t(:) into the (k-1)-th column of AA. */

#line 575 "AB13MD.f"
	i__2 = *n * *n;
#line 575 "AB13MD.f"
	zcopy_(&i__2, &zwork[iz4 + 1], &c__1, &zwork[iz6 + 1 + (k - 2) * *n * 
		*n], &c__1);

/*        Copy diag(P) into the (k-1)-th column of BBd. */

#line 580 "AB13MD.f"
	dcopy_(n, &dwork[iw6 + 1], &c__1, &dwork[iw8 + 1 + (k - 2) * *n], &
		c__1);

/*        Copy P(:) into the (k-1)-th column of BB. */

#line 584 "AB13MD.f"
	i__2 = *n;
#line 584 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 585 "AB13MD.f"
	    dwork[iw9 + i__ + (i__ - 1) * *n + (k - 2) * *n * *n] = dwork[iw6 
		    + i__];
#line 586 "AB13MD.f"
/* L260: */
#line 586 "AB13MD.f"
	}
#line 587 "AB13MD.f"
/* L270: */
#line 587 "AB13MD.f"
    }

#line 589 "AB13MD.f"
    l = 0;

#line 591 "AB13MD.f"
    i__1 = *m;
#line 591 "AB13MD.f"
    for (k = 1; k <= i__1; ++k) {
#line 592 "AB13MD.f"
	if (itype[k] == 1) {
#line 593 "AB13MD.f"
	    ++l;

/*           Set P. */

#line 597 "AB13MD.f"
	    i__2 = iwork[k];
#line 597 "AB13MD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 598 "AB13MD.f"
		dwork[iw6 + i__] = 0.;
#line 599 "AB13MD.f"
/* L280: */
#line 599 "AB13MD.f"
	    }
#line 600 "AB13MD.f"
	    i__2 = iwork[k] + nblock[k];
#line 600 "AB13MD.f"
	    for (i__ = iwork[k] + 1; i__ <= i__2; ++i__) {
#line 601 "AB13MD.f"
		dwork[iw6 + i__] = 1.;
#line 602 "AB13MD.f"
/* L290: */
#line 602 "AB13MD.f"
	    }
#line 603 "AB13MD.f"
	    if (k < *m) {
#line 604 "AB13MD.f"
		i__2 = *n;
#line 604 "AB13MD.f"
		for (i__ = iwork[k + 1] + 1; i__ <= i__2; ++i__) {
#line 605 "AB13MD.f"
		    dwork[iw6 + i__] = 0.;
#line 606 "AB13MD.f"
/* L300: */
#line 606 "AB13MD.f"
		}
#line 607 "AB13MD.f"
	    }

/*           Compute P*Z. */

#line 611 "AB13MD.f"
	    i__2 = *n;
#line 611 "AB13MD.f"
	    for (j = 1; j <= i__2; ++j) {
#line 612 "AB13MD.f"
		i__3 = *n;
#line 612 "AB13MD.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 613 "AB13MD.f"
		    i__4 = iz3 + i__ + (j - 1) * *n;
#line 613 "AB13MD.f"
		    i__5 = iw6 + i__;
#line 613 "AB13MD.f"
		    z__2.r = dwork[i__5], z__2.i = 0.;
#line 613 "AB13MD.f"
		    i__6 = i__ + (j - 1) * *n;
#line 613 "AB13MD.f"
		    z__1.r = z__2.r * zwork[i__6].r - z__2.i * zwork[i__6].i, 
			    z__1.i = z__2.r * zwork[i__6].i + z__2.i * zwork[
			    i__6].r;
#line 613 "AB13MD.f"
		    zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
#line 615 "AB13MD.f"
/* L310: */
#line 615 "AB13MD.f"
		}
#line 616 "AB13MD.f"
/* L320: */
#line 616 "AB13MD.f"
	    }

/*           Compute t = sqrt(-1)*( P*Z - Z'*P ). */

#line 620 "AB13MD.f"
	    i__2 = *n;
#line 620 "AB13MD.f"
	    for (j = 1; j <= i__2; ++j) {
#line 621 "AB13MD.f"
		i__3 = j;
#line 621 "AB13MD.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 622 "AB13MD.f"
		    i__4 = iz3 + i__ + (j - 1) * *n;
#line 622 "AB13MD.f"
		    tempij.r = zwork[i__4].r, tempij.i = zwork[i__4].i;
#line 623 "AB13MD.f"
		    i__4 = iz3 + j + (i__ - 1) * *n;
#line 623 "AB13MD.f"
		    tempji.r = zwork[i__4].r, tempji.i = zwork[i__4].i;
#line 624 "AB13MD.f"
		    i__4 = iz4 + i__ + (j - 1) * *n;
#line 624 "AB13MD.f"
		    d_cnjg(&z__3, &tempji);
#line 624 "AB13MD.f"
		    z__2.r = tempij.r - z__3.r, z__2.i = tempij.i - z__3.i;
#line 624 "AB13MD.f"
		    z__1.r = z__2.r * 0. - z__2.i * 1., z__1.i = z__2.i * 0. 
			    + z__2.r * 1.;
#line 624 "AB13MD.f"
		    zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
#line 626 "AB13MD.f"
		    i__4 = iz4 + j + (i__ - 1) * *n;
#line 626 "AB13MD.f"
		    d_cnjg(&z__3, &tempij);
#line 626 "AB13MD.f"
		    z__2.r = tempji.r - z__3.r, z__2.i = tempji.i - z__3.i;
#line 626 "AB13MD.f"
		    z__1.r = z__2.r * 0. - z__2.i * 1., z__1.i = z__2.i * 0. 
			    + z__2.r * 1.;
#line 626 "AB13MD.f"
		    zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
#line 628 "AB13MD.f"
/* L330: */
#line 628 "AB13MD.f"
		}
#line 629 "AB13MD.f"
/* L340: */
#line 629 "AB13MD.f"
	    }

/*           Copy t(:) into the (m-1+l)-th column of AA. */

#line 633 "AB13MD.f"
	    i__2 = *n * *n;
#line 633 "AB13MD.f"
	    zcopy_(&i__2, &zwork[iz4 + 1], &c__1, &zwork[iz6 + 1 + (*m - 2 + 
		    l) * *n * *n], &c__1);
#line 635 "AB13MD.f"
	}
#line 636 "AB13MD.f"
/* L350: */
#line 636 "AB13MD.f"
    }

/*     Set initial X. */

#line 640 "AB13MD.f"
    i__1 = *m - 1;
#line 640 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 641 "AB13MD.f"
	x[i__] = 1.;
#line 642 "AB13MD.f"
/* L360: */
#line 642 "AB13MD.f"
    }
#line 643 "AB13MD.f"
    if (mr > 0) {
#line 644 "AB13MD.f"
	if (! xfact) {
#line 645 "AB13MD.f"
	    i__1 = mr;
#line 645 "AB13MD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 646 "AB13MD.f"
		x[*m - 1 + i__] = 0.;
#line 647 "AB13MD.f"
/* L370: */
#line 647 "AB13MD.f"
	    }
#line 648 "AB13MD.f"
	} else {
#line 649 "AB13MD.f"
	    l = 0;
#line 650 "AB13MD.f"
	    i__1 = *m;
#line 650 "AB13MD.f"
	    for (k = 1; k <= i__1; ++k) {
#line 651 "AB13MD.f"
		if (itype[k] == 1) {
#line 652 "AB13MD.f"
		    ++l;
/* Computing 2nd power */
#line 653 "AB13MD.f"
		    d__1 = dwork[iw2 + k];
#line 653 "AB13MD.f"
		    x[*m - 1 + l] /= d__1 * d__1;
#line 654 "AB13MD.f"
		}
#line 655 "AB13MD.f"
/* L380: */
#line 655 "AB13MD.f"
	    }
#line 656 "AB13MD.f"
	}
#line 657 "AB13MD.f"
    }

/*     Set constants. */

#line 661 "AB13MD.f"
    svlam = 1. / eps;
#line 662 "AB13MD.f"
    c__ = 1.;

/*     Set H. */

#line 666 "AB13MD.f"
    dlaset_("Full", &mt, &mt, &c_b15, &c_b11, &dwork[iw11 + 1], &mt, (ftnlen)
	    4);

#line 668 "AB13MD.f"
    iter = -1;

/*     Main iteration loop. */

#line 672 "AB13MD.f"
L390:
#line 672 "AB13MD.f"
    ++iter;

/*        Compute A(:) = A0 + AA*x. */

#line 676 "AB13MD.f"
    i__1 = mt;
#line 676 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 677 "AB13MD.f"
	i__2 = iz10 + i__;
#line 677 "AB13MD.f"
	i__3 = i__;
#line 677 "AB13MD.f"
	z__1.r = x[i__3], z__1.i = 0.;
#line 677 "AB13MD.f"
	zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
#line 678 "AB13MD.f"
/* L400: */
#line 678 "AB13MD.f"
    }
#line 679 "AB13MD.f"
    i__1 = *n * *n;
#line 679 "AB13MD.f"
    zcopy_(&i__1, &zwork[iz5 + 1], &c__1, &zwork[iz7 + 1], &c__1);
#line 680 "AB13MD.f"
    i__1 = *n * *n;
#line 680 "AB13MD.f"
    i__2 = *n * *n;
#line 680 "AB13MD.f"
    zgemv_("N", &i__1, &mt, &c_b2, &zwork[iz6 + 1], &i__2, &zwork[iz10 + 1], &
	    c__1, &c_b2, &zwork[iz7 + 1], &c__1, (ftnlen)1);

/*        Compute diag( Binv ). */

#line 685 "AB13MD.f"
    dcopy_(n, &dwork[iw7 + 1], &c__1, &dwork[iw12 + 1], &c__1);
#line 686 "AB13MD.f"
    i__1 = *m - 1;
#line 686 "AB13MD.f"
    dgemv_("N", n, &i__1, &c_b11, &dwork[iw8 + 1], n, &x[1], &c__1, &c_b11, &
	    dwork[iw12 + 1], &c__1, (ftnlen)1);
#line 688 "AB13MD.f"
    i__1 = *n;
#line 688 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 689 "AB13MD.f"
	dwork[iw12 + i__] = 1. / dwork[iw12 + i__];
#line 690 "AB13MD.f"
/* L410: */
#line 690 "AB13MD.f"
    }

/*        Compute Binv*A. */

#line 694 "AB13MD.f"
    i__1 = *n;
#line 694 "AB13MD.f"
    for (j = 1; j <= i__1; ++j) {
#line 695 "AB13MD.f"
	i__2 = *n;
#line 695 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 696 "AB13MD.f"
	    i__3 = iz11 + i__ + (j - 1) * *n;
#line 696 "AB13MD.f"
	    i__4 = iw12 + i__;
#line 696 "AB13MD.f"
	    z__2.r = dwork[i__4], z__2.i = 0.;
#line 696 "AB13MD.f"
	    i__5 = iz7 + i__ + (j - 1) * *n;
#line 696 "AB13MD.f"
	    z__1.r = z__2.r * zwork[i__5].r - z__2.i * zwork[i__5].i, z__1.i =
		     z__2.r * zwork[i__5].i + z__2.i * zwork[i__5].r;
#line 696 "AB13MD.f"
	    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 698 "AB13MD.f"
/* L420: */
#line 698 "AB13MD.f"
	}
#line 699 "AB13MD.f"
/* L430: */
#line 699 "AB13MD.f"
    }

/*        Compute eig( Binv*A ). */

#line 703 "AB13MD.f"
    i__1 = *lzwork - izwrk;
#line 703 "AB13MD.f"
    zgees_("N", "N", (L_fp)select_, n, &zwork[iz11 + 1], n, &sdim, &zwork[
	    iz12 + 1], &zwork[1], n, &zwork[izwrk + 1], &i__1, &dwork[iwrk + 
	    1], bwork, &info2, (ftnlen)1, (ftnlen)1);
#line 706 "AB13MD.f"
    if (info2 > 0) {
#line 707 "AB13MD.f"
	*info = 6;
#line 708 "AB13MD.f"
	return 0;
#line 709 "AB13MD.f"
    }
#line 710 "AB13MD.f"
    i__1 = izwrk + 1;
#line 710 "AB13MD.f"
    lza = (integer) zwork[i__1].r;
#line 711 "AB13MD.f"
    lzamax = max(lza,lzamax);
#line 712 "AB13MD.f"
    i__1 = iz12 + 1;
#line 712 "AB13MD.f"
    e = zwork[i__1].r;
#line 713 "AB13MD.f"
    if (*n > 1) {
#line 714 "AB13MD.f"
	i__1 = *n;
#line 714 "AB13MD.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 715 "AB13MD.f"
	    i__2 = iz12 + i__;
#line 715 "AB13MD.f"
	    if (zwork[i__2].r > e) {
#line 715 "AB13MD.f"
		i__3 = iz12 + i__;
#line 715 "AB13MD.f"
		e = zwork[i__3].r;
#line 715 "AB13MD.f"
	    }
#line 717 "AB13MD.f"
/* L440: */
#line 717 "AB13MD.f"
	}
#line 718 "AB13MD.f"
    }

/*        Set tau. */

#line 722 "AB13MD.f"
    if (mr > 0) {
#line 723 "AB13MD.f"
	snorm = (d__1 = x[*m], abs(d__1));
#line 724 "AB13MD.f"
	if (mr > 1) {
#line 725 "AB13MD.f"
	    i__1 = mt;
#line 725 "AB13MD.f"
	    for (i__ = *m + 1; i__ <= i__1; ++i__) {
#line 726 "AB13MD.f"
		if ((d__1 = x[i__], abs(d__1)) > snorm) {
#line 726 "AB13MD.f"
		    snorm = (d__2 = x[i__], abs(d__2));
#line 726 "AB13MD.f"
		}
#line 727 "AB13MD.f"
/* L450: */
#line 727 "AB13MD.f"
	    }
#line 728 "AB13MD.f"
	}
#line 729 "AB13MD.f"
	if (snorm > 40.) {
#line 730 "AB13MD.f"
	    tau = 100.;
#line 731 "AB13MD.f"
	} else if (snorm > 8.) {
#line 732 "AB13MD.f"
	    tau = 50.;
#line 733 "AB13MD.f"
	} else if (snorm > 4.) {
#line 734 "AB13MD.f"
	    tau = 10.;
#line 735 "AB13MD.f"
	} else if (snorm > 1.) {
#line 736 "AB13MD.f"
	    tau = 5.;
#line 737 "AB13MD.f"
	} else {
#line 738 "AB13MD.f"
	    tau = 2.;
#line 739 "AB13MD.f"
	}
#line 740 "AB13MD.f"
    }
#line 741 "AB13MD.f"
    if (iter == 0) {
#line 742 "AB13MD.f"
	dlambd = e + .001;
#line 743 "AB13MD.f"
    } else {
#line 744 "AB13MD.f"
	dwork[iw13 + 1] = e;
#line 745 "AB13MD.f"
	dcopy_(&mt, &x[1], &c__1, &dwork[iw13 + 2], &c__1);
#line 746 "AB13MD.f"
	dlambd = dwork[iw13 + 1] * .98999999999999999 + dwork[iw14 + 1] * .01;
#line 748 "AB13MD.f"
	dcopy_(&mt, &dwork[iw13 + 2], &c__1, &dwork[iw18 + 1], &c__1);
#line 749 "AB13MD.f"
	dcopy_(&mt, &dwork[iw14 + 2], &c__1, &dwork[iw19 + 1], &c__1);
#line 750 "AB13MD.f"
	l = 0;
#line 751 "AB13MD.f"
L460:
#line 751 "AB13MD.f"
	i__1 = mt;
#line 751 "AB13MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 752 "AB13MD.f"
	    x[i__] = (1. - .01 / pow_di(&c_b136, &l)) * dwork[iw18 + i__] + 
		    .01 / pow_di(&c_b136, &l) * dwork[iw19 + i__];
#line 754 "AB13MD.f"
/* L470: */
#line 754 "AB13MD.f"
	}

/*           Compute At(:) = A0 + AA*x. */

#line 758 "AB13MD.f"
	i__1 = mt;
#line 758 "AB13MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 759 "AB13MD.f"
	    i__2 = iz10 + i__;
#line 759 "AB13MD.f"
	    i__3 = i__;
#line 759 "AB13MD.f"
	    z__1.r = x[i__3], z__1.i = 0.;
#line 759 "AB13MD.f"
	    zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
#line 760 "AB13MD.f"
/* L480: */
#line 760 "AB13MD.f"
	}
#line 761 "AB13MD.f"
	i__1 = *n * *n;
#line 761 "AB13MD.f"
	zcopy_(&i__1, &zwork[iz5 + 1], &c__1, &zwork[iz9 + 1], &c__1);
#line 762 "AB13MD.f"
	i__1 = *n * *n;
#line 762 "AB13MD.f"
	i__2 = *n * *n;
#line 762 "AB13MD.f"
	zgemv_("N", &i__1, &mt, &c_b2, &zwork[iz6 + 1], &i__2, &zwork[iz10 + 
		1], &c__1, &c_b2, &zwork[iz9 + 1], &c__1, (ftnlen)1);

/*           Compute diag(Bt). */

#line 767 "AB13MD.f"
	dcopy_(n, &dwork[iw7 + 1], &c__1, &dwork[iw21 + 1], &c__1);
#line 768 "AB13MD.f"
	i__1 = *m - 1;
#line 768 "AB13MD.f"
	dgemv_("N", n, &i__1, &c_b11, &dwork[iw8 + 1], n, &x[1], &c__1, &
		c_b11, &dwork[iw21 + 1], &c__1, (ftnlen)1);

/*           Compute W. */

#line 773 "AB13MD.f"
	i__1 = *n;
#line 773 "AB13MD.f"
	for (j = 1; j <= i__1; ++j) {
#line 774 "AB13MD.f"
	    i__2 = *n;
#line 774 "AB13MD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 775 "AB13MD.f"
		if (i__ == j) {
#line 776 "AB13MD.f"
		    i__3 = iz13 + i__ + (i__ - 1) * *n;
#line 776 "AB13MD.f"
		    d__1 = (dwork[iw14 + 1] - dwork[iw13 + 1]) * 1e-4 / 2. - 
			    dlambd * dwork[iw21 + i__];
#line 776 "AB13MD.f"
		    z__2.r = d__1, z__2.i = 0.;
#line 776 "AB13MD.f"
		    i__4 = iz9 + i__ + (i__ - 1) * *n;
#line 776 "AB13MD.f"
		    z__1.r = z__2.r + zwork[i__4].r, z__1.i = z__2.i + zwork[
			    i__4].i;
#line 776 "AB13MD.f"
		    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 780 "AB13MD.f"
		} else {
#line 781 "AB13MD.f"
		    i__3 = iz13 + i__ + (j - 1) * *n;
#line 781 "AB13MD.f"
		    i__4 = iz9 + i__ + (j - 1) * *n;
#line 781 "AB13MD.f"
		    zwork[i__3].r = zwork[i__4].r, zwork[i__3].i = zwork[i__4]
			    .i;
#line 782 "AB13MD.f"
		}
#line 783 "AB13MD.f"
/* L490: */
#line 783 "AB13MD.f"
	    }
#line 784 "AB13MD.f"
/* L500: */
#line 784 "AB13MD.f"
	}

/*           Compute eig( W ). */

#line 788 "AB13MD.f"
	i__1 = *lzwork - izwrk;
#line 788 "AB13MD.f"
	zgees_("N", "N", (L_fp)select_, n, &zwork[iz13 + 1], n, &sdim, &zwork[
		iz14 + 1], &zwork[1], n, &zwork[izwrk + 1], &i__1, &dwork[
		iwrk + 1], bwork, &info2, (ftnlen)1, (ftnlen)1);
#line 791 "AB13MD.f"
	if (info2 > 0) {
#line 792 "AB13MD.f"
	    *info = 6;
#line 793 "AB13MD.f"
	    return 0;
#line 794 "AB13MD.f"
	}
#line 795 "AB13MD.f"
	i__1 = izwrk + 1;
#line 795 "AB13MD.f"
	lza = (integer) zwork[i__1].r;
#line 796 "AB13MD.f"
	lzamax = max(lza,lzamax);
#line 797 "AB13MD.f"
	i__1 = iz14 + 1;
#line 797 "AB13MD.f"
	emax = zwork[i__1].r;
#line 798 "AB13MD.f"
	if (*n > 1) {
#line 799 "AB13MD.f"
	    i__1 = *n;
#line 799 "AB13MD.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 800 "AB13MD.f"
		i__2 = iz14 + i__;
#line 800 "AB13MD.f"
		if (zwork[i__2].r > emax) {
#line 800 "AB13MD.f"
		    i__3 = iz14 + i__;
#line 800 "AB13MD.f"
		    emax = zwork[i__3].r;
#line 800 "AB13MD.f"
		}
#line 802 "AB13MD.f"
/* L510: */
#line 802 "AB13MD.f"
	    }
#line 803 "AB13MD.f"
	}
#line 804 "AB13MD.f"
	if (emax <= 0.) {
#line 805 "AB13MD.f"
	    goto L515;
#line 806 "AB13MD.f"
	} else {
#line 807 "AB13MD.f"
	    ++l;
#line 808 "AB13MD.f"
	    goto L460;
#line 809 "AB13MD.f"
	}
#line 810 "AB13MD.f"
    }

/*        Set y. */

#line 814 "AB13MD.f"
L515:
#line 814 "AB13MD.f"
    dwork[iw13 + 1] = dlambd;
#line 815 "AB13MD.f"
    dcopy_(&mt, &x[1], &c__1, &dwork[iw13 + 2], &c__1);

#line 817 "AB13MD.f"
    if (svlam - dlambd < tol) {
#line 818 "AB13MD.f"
	*bound = sqrt((max(e,0.))) * znorm;
#line 819 "AB13MD.f"
	i__1 = *m - 1;
#line 819 "AB13MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
#line 820 "AB13MD.f"
	    d__1 = dwork[iw2 + i__ + 1];
#line 820 "AB13MD.f"
	    x[i__] *= d__1 * d__1;
#line 821 "AB13MD.f"
/* L520: */
#line 821 "AB13MD.f"
	}

/*           Compute sqrt( x ). */

#line 825 "AB13MD.f"
	i__1 = *m - 1;
#line 825 "AB13MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 826 "AB13MD.f"
	    dwork[iw20 + i__] = sqrt(x[i__]);
#line 827 "AB13MD.f"
/* L530: */
#line 827 "AB13MD.f"
	}

/*           Compute diag( D ). */

#line 831 "AB13MD.f"
	dcopy_(n, &dwork[iw7 + 1], &c__1, &d__[1], &c__1);
#line 832 "AB13MD.f"
	i__1 = *m - 1;
#line 832 "AB13MD.f"
	dgemv_("N", n, &i__1, &c_b11, &dwork[iw8 + 1], n, &dwork[iw20 + 1], &
		c__1, &c_b11, &d__[1], &c__1, (ftnlen)1);

/*           Compute diag( G ). */

#line 837 "AB13MD.f"
	j = 0;
#line 838 "AB13MD.f"
	l = 0;
#line 839 "AB13MD.f"
	i__1 = *m;
#line 839 "AB13MD.f"
	for (k = 1; k <= i__1; ++k) {
#line 840 "AB13MD.f"
	    j += nblock[k];
#line 841 "AB13MD.f"
	    if (itype[k] == 1) {
#line 842 "AB13MD.f"
		++l;
/* Computing 2nd power */
#line 843 "AB13MD.f"
		d__1 = dwork[iw2 + k];
#line 843 "AB13MD.f"
		x[*m - 1 + l] *= d__1 * d__1;
#line 844 "AB13MD.f"
		g[j] = x[*m - 1 + l];
#line 845 "AB13MD.f"
	    }
#line 846 "AB13MD.f"
/* L540: */
#line 846 "AB13MD.f"
	}
#line 847 "AB13MD.f"
	dscal_(n, &znorm, &g[1], &c__1);
#line 848 "AB13MD.f"
	dwork[1] = (doublereal) (minwrk - *n * 5 + lwamax);
#line 849 "AB13MD.f"
	i__1 = minzrk - *n * 3 + lzamax;
#line 849 "AB13MD.f"
	z__1.r = (doublereal) i__1, z__1.i = 0.;
#line 849 "AB13MD.f"
	zwork[1].r = z__1.r, zwork[1].i = z__1.i;
#line 850 "AB13MD.f"
	return 0;
#line 851 "AB13MD.f"
    }
#line 852 "AB13MD.f"
    svlam = dlambd;
#line 853 "AB13MD.f"
    i__1 = *m;
#line 853 "AB13MD.f"
    for (k = 1; k <= i__1; ++k) {

/*           Store xD. */

#line 857 "AB13MD.f"
	i__2 = *m - 1;
#line 857 "AB13MD.f"
	dcopy_(&i__2, &x[1], &c__1, &dwork[iw22 + 1], &c__1);
#line 858 "AB13MD.f"
	if (mr > 0) {

/*              Store xG. */

#line 862 "AB13MD.f"
	    dcopy_(&mr, &x[*m], &c__1, &dwork[iw23 + 1], &c__1);
#line 863 "AB13MD.f"
	}

/*           Compute A(:) = A0 + AA*x. */

#line 867 "AB13MD.f"
	i__2 = mt;
#line 867 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 868 "AB13MD.f"
	    i__3 = iz10 + i__;
#line 868 "AB13MD.f"
	    i__4 = i__;
#line 868 "AB13MD.f"
	    z__1.r = x[i__4], z__1.i = 0.;
#line 868 "AB13MD.f"
	    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 869 "AB13MD.f"
/* L550: */
#line 869 "AB13MD.f"
	}
#line 870 "AB13MD.f"
	i__2 = *n * *n;
#line 870 "AB13MD.f"
	zcopy_(&i__2, &zwork[iz5 + 1], &c__1, &zwork[iz7 + 1], &c__1);
#line 871 "AB13MD.f"
	i__2 = *n * *n;
#line 871 "AB13MD.f"
	i__3 = *n * *n;
#line 871 "AB13MD.f"
	zgemv_("N", &i__2, &mt, &c_b2, &zwork[iz6 + 1], &i__3, &zwork[iz10 + 
		1], &c__1, &c_b2, &zwork[iz7 + 1], &c__1, (ftnlen)1);

/*           Compute B = B0d + BBd*xD. */

#line 876 "AB13MD.f"
	dcopy_(n, &dwork[iw7 + 1], &c__1, &dwork[iw24 + 1], &c__1);
#line 877 "AB13MD.f"
	i__2 = *m - 1;
#line 877 "AB13MD.f"
	dgemv_("N", n, &i__2, &c_b11, &dwork[iw8 + 1], n, &dwork[iw22 + 1], &
		c__1, &c_b11, &dwork[iw24 + 1], &c__1, (ftnlen)1);

/*           Compute F. */

#line 882 "AB13MD.f"
	i__2 = *n;
#line 882 "AB13MD.f"
	for (j = 1; j <= i__2; ++j) {
#line 883 "AB13MD.f"
	    i__3 = *n;
#line 883 "AB13MD.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 884 "AB13MD.f"
		if (i__ == j) {
#line 885 "AB13MD.f"
		    i__4 = iz15 + i__ + (i__ - 1) * *n;
#line 885 "AB13MD.f"
		    d__1 = dlambd * dwork[iw24 + i__];
#line 885 "AB13MD.f"
		    z__2.r = d__1, z__2.i = 0.;
#line 885 "AB13MD.f"
		    i__5 = iz7 + i__ + (i__ - 1) * *n;
#line 885 "AB13MD.f"
		    z__1.r = z__2.r - zwork[i__5].r, z__1.i = z__2.i - zwork[
			    i__5].i;
#line 885 "AB13MD.f"
		    zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
#line 887 "AB13MD.f"
		} else {
#line 888 "AB13MD.f"
		    i__4 = iz15 + i__ + (j - 1) * *n;
#line 888 "AB13MD.f"
		    i__5 = iz7 + i__ + (j - 1) * *n;
#line 888 "AB13MD.f"
		    z__1.r = -zwork[i__5].r, z__1.i = -zwork[i__5].i;
#line 888 "AB13MD.f"
		    zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
#line 889 "AB13MD.f"
		}
#line 890 "AB13MD.f"
/* L555: */
#line 890 "AB13MD.f"
	    }
#line 891 "AB13MD.f"
/* L556: */
#line 891 "AB13MD.f"
	}
#line 892 "AB13MD.f"
	zlacpy_("Full", n, n, &zwork[iz15 + 1], n, &zwork[iz17 + 1], n, (
		ftnlen)4);

/*           Compute det( F ). */

#line 897 "AB13MD.f"
	i__2 = *lzwork - izwrk;
#line 897 "AB13MD.f"
	zgees_("N", "N", (L_fp)select_, n, &zwork[iz15 + 1], n, &sdim, &zwork[
		iz16 + 1], &zwork[1], n, &zwork[izwrk + 1], &i__2, &dwork[
		iwrk + 1], bwork, &info2, (ftnlen)1, (ftnlen)1);
#line 900 "AB13MD.f"
	if (info2 > 0) {
#line 901 "AB13MD.f"
	    *info = 6;
#line 902 "AB13MD.f"
	    return 0;
#line 903 "AB13MD.f"
	}
#line 904 "AB13MD.f"
	i__2 = izwrk + 1;
#line 904 "AB13MD.f"
	lza = (integer) zwork[i__2].r;
#line 905 "AB13MD.f"
	lzamax = max(lza,lzamax);
#line 906 "AB13MD.f"
	detf.r = 1., detf.i = 0.;
#line 907 "AB13MD.f"
	i__2 = *n;
#line 907 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 908 "AB13MD.f"
	    i__3 = iz16 + i__;
#line 908 "AB13MD.f"
	    z__1.r = detf.r * zwork[i__3].r - detf.i * zwork[i__3].i, z__1.i =
		     detf.r * zwork[i__3].i + detf.i * zwork[i__3].r;
#line 908 "AB13MD.f"
	    detf.r = z__1.r, detf.i = z__1.i;
#line 909 "AB13MD.f"
/* L560: */
#line 909 "AB13MD.f"
	}

/*           Compute Finv. */

#line 913 "AB13MD.f"
	zgetrf_(n, n, &zwork[iz17 + 1], n, &iwork[1], &info2);
#line 914 "AB13MD.f"
	if (info2 > 0) {
#line 915 "AB13MD.f"
	    *info = 5;
#line 916 "AB13MD.f"
	    return 0;
#line 917 "AB13MD.f"
	}
#line 918 "AB13MD.f"
	i__2 = *ldwork - iwrk;
#line 918 "AB13MD.f"
	zgetri_(n, &zwork[iz17 + 1], n, &iwork[1], &zwork[izwrk + 1], &i__2, &
		info2);
#line 920 "AB13MD.f"
	i__2 = izwrk + 1;
#line 920 "AB13MD.f"
	lza = (integer) zwork[i__2].r;
#line 921 "AB13MD.f"
	lzamax = max(lza,lzamax);

/*           Compute phi. */

#line 925 "AB13MD.f"
	i__2 = *m - 1;
#line 925 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 926 "AB13MD.f"
	    dwork[iw25 + i__] = dwork[iw22 + i__] - .01;
#line 927 "AB13MD.f"
	    dwork[iw25 + *m - 1 + i__] = 100. - dwork[iw22 + i__];
#line 928 "AB13MD.f"
/* L570: */
#line 928 "AB13MD.f"
	}
#line 929 "AB13MD.f"
	if (mr > 0) {
#line 930 "AB13MD.f"
	    i__2 = mr;
#line 930 "AB13MD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 931 "AB13MD.f"
		dwork[iw25 + (*m - 1 << 1) + i__] = dwork[iw23 + i__] + tau;
#line 932 "AB13MD.f"
		dwork[iw25 + (*m - 1 << 1) + mr + i__] = tau - dwork[iw23 + 
			i__];
#line 933 "AB13MD.f"
/* L580: */
#line 933 "AB13MD.f"
	    }
#line 934 "AB13MD.f"
	}
#line 935 "AB13MD.f"
	prod = 1.;
#line 936 "AB13MD.f"
	i__2 = mt << 1;
#line 936 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 937 "AB13MD.f"
	    prod *= dwork[iw25 + i__];
#line 938 "AB13MD.f"
/* L590: */
#line 938 "AB13MD.f"
	}
#line 939 "AB13MD.f"
	temp = detf.r;
#line 940 "AB13MD.f"
	if (temp < eps) {
#line 940 "AB13MD.f"
	    temp = eps;
#line 940 "AB13MD.f"
	}
#line 941 "AB13MD.f"
	phi = -log(temp) - log(prod);

/*           Compute g. */

#line 945 "AB13MD.f"
	i__2 = mt;
#line 945 "AB13MD.f"
	for (j = 1; j <= i__2; ++j) {
#line 946 "AB13MD.f"
	    i__3 = *n * *n;
#line 946 "AB13MD.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 947 "AB13MD.f"
		i__4 = iz18 + i__ + (j - 1) * *n * *n;
#line 947 "AB13MD.f"
		d__1 = dlambd * dwork[iw9 + i__ + (j - 1) * *n * *n];
#line 947 "AB13MD.f"
		z__2.r = d__1, z__2.i = 0.;
#line 947 "AB13MD.f"
		i__5 = iz6 + i__ + (j - 1) * *n * *n;
#line 947 "AB13MD.f"
		z__1.r = z__2.r - zwork[i__5].r, z__1.i = z__2.i - zwork[i__5]
			.i;
#line 947 "AB13MD.f"
		zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
#line 949 "AB13MD.f"
/* L600: */
#line 949 "AB13MD.f"
	    }
#line 950 "AB13MD.f"
/* L610: */
#line 950 "AB13MD.f"
	}
#line 951 "AB13MD.f"
	i__2 = *n * *n;
#line 951 "AB13MD.f"
	i__3 = *n * *n;
#line 951 "AB13MD.f"
	zgemv_("C", &i__2, &mt, &c_b2, &zwork[iz18 + 1], &i__3, &zwork[iz17 + 
		1], &c__1, &c_b1, &zwork[iz19 + 1], &c__1, (ftnlen)1);
#line 953 "AB13MD.f"
	i__2 = *m - 1;
#line 953 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 954 "AB13MD.f"
	    dwork[iw26 + i__] = 1. / (dwork[iw22 + i__] - .01) - 1. / (100. - 
		    dwork[iw22 + i__]);
#line 956 "AB13MD.f"
/* L620: */
#line 956 "AB13MD.f"
	}
#line 957 "AB13MD.f"
	if (mr > 0) {
#line 958 "AB13MD.f"
	    i__2 = mr;
#line 958 "AB13MD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 959 "AB13MD.f"
		dwork[iw26 + *m - 1 + i__] = 1. / (dwork[iw23 + i__] + tau) - 
			1. / (tau - dwork[iw23 + i__]);
#line 961 "AB13MD.f"
/* L630: */
#line 961 "AB13MD.f"
	    }
#line 962 "AB13MD.f"
	}
#line 963 "AB13MD.f"
	i__2 = mt;
#line 963 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 964 "AB13MD.f"
	    i__3 = iz19 + i__;
#line 964 "AB13MD.f"
	    dwork[iw26 + i__] = -zwork[i__3].r - dwork[iw26 + i__];
#line 966 "AB13MD.f"
/* L640: */
#line 966 "AB13MD.f"
	}

/*           Compute h. */

#line 970 "AB13MD.f"
	dlacpy_("Full", &mt, &mt, &dwork[iw11 + 1], &mt, &dwork[iw31 + 1], &
		mt, (ftnlen)4);
#line 972 "AB13MD.f"
	dcopy_(&mt, &dwork[iw26 + 1], &c__1, &dwork[iw27 + 1], &c__1);
#line 973 "AB13MD.f"
	i__2 = *ldwork - iwrk;
#line 973 "AB13MD.f"
	dsysv_("U", &mt, &c__1, &dwork[iw31 + 1], &mt, &iwork[1], &dwork[iw27 
		+ 1], &mt, &dwork[iwrk + 1], &i__2, &info2, (ftnlen)1);
#line 976 "AB13MD.f"
	if (info2 > 0) {
#line 977 "AB13MD.f"
	    *info = 5;
#line 978 "AB13MD.f"
	    return 0;
#line 979 "AB13MD.f"
	}
#line 980 "AB13MD.f"
	lwa = (integer) dwork[iwrk + 1];
#line 981 "AB13MD.f"
	lwamax = max(lwa,lwamax);
#line 982 "AB13MD.f"
	stsize = 1.;

/*           Store hD. */

#line 986 "AB13MD.f"
	i__2 = *m - 1;
#line 986 "AB13MD.f"
	dcopy_(&i__2, &dwork[iw27 + 1], &c__1, &dwork[iw28 + 1], &c__1);

/*           Determine stepsize. */

#line 990 "AB13MD.f"
	l = 0;
#line 991 "AB13MD.f"
	i__2 = *m - 1;
#line 991 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 992 "AB13MD.f"
	    if (dwork[iw28 + i__] > 0.) {
#line 993 "AB13MD.f"
		++l;
#line 994 "AB13MD.f"
		if (l == 1) {
#line 995 "AB13MD.f"
		    temp = (dwork[iw22 + i__] - .01) / dwork[iw28 + i__];
#line 996 "AB13MD.f"
		} else {
/* Computing MIN */
#line 997 "AB13MD.f"
		    d__1 = temp, d__2 = (dwork[iw22 + i__] - .01) / dwork[
			    iw28 + i__];
#line 997 "AB13MD.f"
		    temp = min(d__1,d__2);
#line 999 "AB13MD.f"
		}
#line 1000 "AB13MD.f"
	    }
#line 1001 "AB13MD.f"
/* L650: */
#line 1001 "AB13MD.f"
	}
#line 1002 "AB13MD.f"
	if (l > 0) {
#line 1002 "AB13MD.f"
	    stsize = min(stsize,temp);
#line 1002 "AB13MD.f"
	}
#line 1003 "AB13MD.f"
	l = 0;
#line 1004 "AB13MD.f"
	i__2 = *m - 1;
#line 1004 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 1005 "AB13MD.f"
	    if (dwork[iw28 + i__] < 0.) {
#line 1006 "AB13MD.f"
		++l;
#line 1007 "AB13MD.f"
		if (l == 1) {
#line 1008 "AB13MD.f"
		    temp = (100. - dwork[iw22 + i__]) / (-dwork[iw28 + i__]);
#line 1010 "AB13MD.f"
		} else {
/* Computing MIN */
#line 1011 "AB13MD.f"
		    d__1 = temp, d__2 = (100. - dwork[iw22 + i__]) / (-dwork[
			    iw28 + i__]);
#line 1011 "AB13MD.f"
		    temp = min(d__1,d__2);
#line 1013 "AB13MD.f"
		}
#line 1014 "AB13MD.f"
	    }
#line 1015 "AB13MD.f"
/* L660: */
#line 1015 "AB13MD.f"
	}
#line 1016 "AB13MD.f"
	if (l > 0) {
#line 1016 "AB13MD.f"
	    stsize = min(stsize,temp);
#line 1016 "AB13MD.f"
	}
#line 1017 "AB13MD.f"
	if (mr > 0) {

/*              Store hG. */

#line 1021 "AB13MD.f"
	    dcopy_(&mr, &dwork[iw27 + *m], &c__1, &dwork[iw29 + 1], &c__1);

/*              Determine stepsize. */

#line 1025 "AB13MD.f"
	    l = 0;
#line 1026 "AB13MD.f"
	    i__2 = mr;
#line 1026 "AB13MD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 1027 "AB13MD.f"
		if (dwork[iw29 + i__] > 0.) {
#line 1028 "AB13MD.f"
		    ++l;
#line 1029 "AB13MD.f"
		    if (l == 1) {
#line 1030 "AB13MD.f"
			temp = (dwork[iw23 + i__] + tau) / dwork[iw29 + i__];
#line 1032 "AB13MD.f"
		    } else {
/* Computing MIN */
#line 1033 "AB13MD.f"
			d__1 = temp, d__2 = (dwork[iw23 + i__] + tau) / dwork[
				iw29 + i__];
#line 1033 "AB13MD.f"
			temp = min(d__1,d__2);
#line 1035 "AB13MD.f"
		    }
#line 1036 "AB13MD.f"
		}
#line 1037 "AB13MD.f"
/* L670: */
#line 1037 "AB13MD.f"
	    }
#line 1038 "AB13MD.f"
	    if (l > 0) {
#line 1038 "AB13MD.f"
		stsize = min(stsize,temp);
#line 1038 "AB13MD.f"
	    }
#line 1039 "AB13MD.f"
	    l = 0;
#line 1040 "AB13MD.f"
	    i__2 = mr;
#line 1040 "AB13MD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 1041 "AB13MD.f"
		if (dwork[iw29 + i__] < 0.) {
#line 1042 "AB13MD.f"
		    ++l;
#line 1043 "AB13MD.f"
		    if (l == 1) {
#line 1044 "AB13MD.f"
			temp = (tau - dwork[iw23 + i__]) / (-dwork[iw29 + i__]
				);
#line 1046 "AB13MD.f"
		    } else {
/* Computing MIN */
#line 1047 "AB13MD.f"
			d__1 = temp, d__2 = (tau - dwork[iw23 + i__]) / (
				-dwork[iw29 + i__]);
#line 1047 "AB13MD.f"
			temp = min(d__1,d__2);
#line 1049 "AB13MD.f"
		    }
#line 1050 "AB13MD.f"
		}
#line 1051 "AB13MD.f"
/* L680: */
#line 1051 "AB13MD.f"
	    }
#line 1052 "AB13MD.f"
	}
#line 1053 "AB13MD.f"
	if (l > 0) {
#line 1053 "AB13MD.f"
	    stsize = min(stsize,temp);
#line 1053 "AB13MD.f"
	}
#line 1054 "AB13MD.f"
	stsize *= .9;
#line 1055 "AB13MD.f"
	if (stsize >= tol4) {

/*              Compute x_new. */

#line 1059 "AB13MD.f"
	    i__2 = mt;
#line 1059 "AB13MD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 1060 "AB13MD.f"
		dwork[iw20 + i__] = x[i__] - stsize * dwork[iw27 + i__];
#line 1061 "AB13MD.f"
/* L700: */
#line 1061 "AB13MD.f"
	    }

/*              Store xD. */

#line 1065 "AB13MD.f"
	    i__2 = *m - 1;
#line 1065 "AB13MD.f"
	    dcopy_(&i__2, &dwork[iw20 + 1], &c__1, &dwork[iw22 + 1], &c__1);
#line 1066 "AB13MD.f"
	    if (mr > 0) {

/*                 Store xG. */

#line 1070 "AB13MD.f"
		dcopy_(&mr, &dwork[iw20 + *m], &c__1, &dwork[iw23 + 1], &c__1)
			;
#line 1072 "AB13MD.f"
	    }

/*              Compute A(:) = A0 + AA*x_new. */

#line 1076 "AB13MD.f"
	    i__2 = mt;
#line 1076 "AB13MD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 1077 "AB13MD.f"
		i__3 = iz10 + i__;
#line 1077 "AB13MD.f"
		i__4 = iw20 + i__;
#line 1077 "AB13MD.f"
		z__1.r = dwork[i__4], z__1.i = 0.;
#line 1077 "AB13MD.f"
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 1078 "AB13MD.f"
/* L710: */
#line 1078 "AB13MD.f"
	    }
#line 1079 "AB13MD.f"
	    i__2 = *n * *n;
#line 1079 "AB13MD.f"
	    zcopy_(&i__2, &zwork[iz5 + 1], &c__1, &zwork[iz7 + 1], &c__1);
#line 1080 "AB13MD.f"
	    i__2 = *n * *n;
#line 1080 "AB13MD.f"
	    i__3 = *n * *n;
#line 1080 "AB13MD.f"
	    zgemv_("N", &i__2, &mt, &c_b2, &zwork[iz6 + 1], &i__3, &zwork[
		    iz10 + 1], &c__1, &c_b2, &zwork[iz7 + 1], &c__1, (ftnlen)
		    1);

/*              Compute B = B0d + BBd*xD. */

#line 1085 "AB13MD.f"
	    dcopy_(n, &dwork[iw7 + 1], &c__1, &dwork[iw24 + 1], &c__1);
#line 1086 "AB13MD.f"
	    i__2 = *m - 1;
#line 1086 "AB13MD.f"
	    dgemv_("N", n, &i__2, &c_b11, &dwork[iw8 + 1], n, &dwork[iw22 + 1]
		    , &c__1, &c_b11, &dwork[iw24 + 1], &c__1, (ftnlen)1);

/*              Compute lambda*diag(B) - A. */

#line 1091 "AB13MD.f"
	    i__2 = *n;
#line 1091 "AB13MD.f"
	    for (j = 1; j <= i__2; ++j) {
#line 1092 "AB13MD.f"
		i__3 = *n;
#line 1092 "AB13MD.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 1093 "AB13MD.f"
		    if (i__ == j) {
#line 1094 "AB13MD.f"
			i__4 = iz15 + i__ + (i__ - 1) * *n;
#line 1094 "AB13MD.f"
			d__1 = dlambd * dwork[iw24 + i__];
#line 1094 "AB13MD.f"
			z__2.r = d__1, z__2.i = 0.;
#line 1094 "AB13MD.f"
			i__5 = iz7 + i__ + (i__ - 1) * *n;
#line 1094 "AB13MD.f"
			z__1.r = z__2.r - zwork[i__5].r, z__1.i = z__2.i - 
				zwork[i__5].i;
#line 1094 "AB13MD.f"
			zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
#line 1096 "AB13MD.f"
		    } else {
#line 1097 "AB13MD.f"
			i__4 = iz15 + i__ + (j - 1) * *n;
#line 1097 "AB13MD.f"
			i__5 = iz7 + i__ + (j - 1) * *n;
#line 1097 "AB13MD.f"
			z__1.r = -zwork[i__5].r, z__1.i = -zwork[i__5].i;
#line 1097 "AB13MD.f"
			zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
#line 1099 "AB13MD.f"
		    }
#line 1100 "AB13MD.f"
/* L720: */
#line 1100 "AB13MD.f"
		}
#line 1101 "AB13MD.f"
/* L730: */
#line 1101 "AB13MD.f"
	    }

/*              Compute eig( lambda*diag(B)-A ). */

#line 1105 "AB13MD.f"
	    i__2 = *lzwork - izwrk;
#line 1105 "AB13MD.f"
	    zgees_("N", "N", (L_fp)select_, n, &zwork[iz15 + 1], n, &sdim, &
		    zwork[iz16 + 1], &zwork[1], n, &zwork[izwrk + 1], &i__2, &
		    dwork[iwrk + 1], bwork, &info2, (ftnlen)1, (ftnlen)1);
#line 1109 "AB13MD.f"
	    if (info2 > 0) {
#line 1110 "AB13MD.f"
		*info = 6;
#line 1111 "AB13MD.f"
		return 0;
#line 1112 "AB13MD.f"
	    }
#line 1113 "AB13MD.f"
	    i__2 = izwrk + 1;
#line 1113 "AB13MD.f"
	    lza = (integer) zwork[i__2].r;
#line 1114 "AB13MD.f"
	    lzamax = max(lza,lzamax);
#line 1115 "AB13MD.f"
	    i__2 = iz16 + 1;
#line 1115 "AB13MD.f"
	    emin = zwork[i__2].r;
#line 1116 "AB13MD.f"
	    if (*n > 1) {
#line 1117 "AB13MD.f"
		i__2 = *n;
#line 1117 "AB13MD.f"
		for (i__ = 2; i__ <= i__2; ++i__) {
#line 1118 "AB13MD.f"
		    i__3 = iz16 + i__;
#line 1118 "AB13MD.f"
		    if (zwork[i__3].r < emin) {
#line 1118 "AB13MD.f"
			i__4 = iz16 + i__;
#line 1118 "AB13MD.f"
			emin = zwork[i__4].r;
#line 1118 "AB13MD.f"
		    }
#line 1120 "AB13MD.f"
/* L740: */
#line 1120 "AB13MD.f"
		}
#line 1121 "AB13MD.f"
	    }
#line 1122 "AB13MD.f"
	    i__2 = *n;
#line 1122 "AB13MD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 1123 "AB13MD.f"
		i__3 = iz16 + i__;
#line 1123 "AB13MD.f"
		dwork[iw30 + i__] = zwork[i__3].r;
#line 1124 "AB13MD.f"
/* L750: */
#line 1124 "AB13MD.f"
	    }
#line 1125 "AB13MD.f"
	    i__2 = *m - 1;
#line 1125 "AB13MD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 1126 "AB13MD.f"
		dwork[iw30 + *n + i__] = dwork[iw22 + i__] - .01;
#line 1127 "AB13MD.f"
		dwork[iw30 + *n + *m - 1 + i__] = 100. - dwork[iw22 + i__];
#line 1128 "AB13MD.f"
/* L760: */
#line 1128 "AB13MD.f"
	    }
#line 1129 "AB13MD.f"
	    if (mr > 0) {
#line 1130 "AB13MD.f"
		i__2 = mr;
#line 1130 "AB13MD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 1131 "AB13MD.f"
		    dwork[iw30 + *n + (*m - 1 << 1) + i__] = dwork[iw23 + i__]
			     + tau;
#line 1132 "AB13MD.f"
		    dwork[iw30 + *n + (*m - 1 << 1) + mr + i__] = tau - dwork[
			    iw23 + i__];
#line 1134 "AB13MD.f"
/* L770: */
#line 1134 "AB13MD.f"
		}
#line 1135 "AB13MD.f"
	    }
#line 1136 "AB13MD.f"
	    prod = 1.;
#line 1137 "AB13MD.f"
	    i__2 = *n + (mt << 1);
#line 1137 "AB13MD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 1138 "AB13MD.f"
		prod *= dwork[iw30 + i__];
#line 1139 "AB13MD.f"
/* L780: */
#line 1139 "AB13MD.f"
	    }
#line 1140 "AB13MD.f"
	    if (emin <= 0. || -log(prod) >= phi) {
#line 1141 "AB13MD.f"
		stsize /= 10.;
#line 1142 "AB13MD.f"
	    } else {
#line 1143 "AB13MD.f"
		dcopy_(&mt, &dwork[iw20 + 1], &c__1, &x[1], &c__1);
#line 1144 "AB13MD.f"
	    }
#line 1145 "AB13MD.f"
	}
#line 1146 "AB13MD.f"
	if (stsize < tol4) {
#line 1146 "AB13MD.f"
	    goto L810;
#line 1146 "AB13MD.f"
	}
#line 1147 "AB13MD.f"
/* L800: */
#line 1147 "AB13MD.f"
    }

#line 1149 "AB13MD.f"
L810:

/*           Store xD. */

#line 1153 "AB13MD.f"
    i__1 = *m - 1;
#line 1153 "AB13MD.f"
    dcopy_(&i__1, &x[1], &c__1, &dwork[iw22 + 1], &c__1);
#line 1154 "AB13MD.f"
    if (mr > 0) {

/*              Store xG. */

#line 1158 "AB13MD.f"
	dcopy_(&mr, &x[*m], &c__1, &dwork[iw23 + 1], &c__1);
#line 1159 "AB13MD.f"
    }

/*           Compute A(:) = A0 + AA*x. */

#line 1163 "AB13MD.f"
    i__1 = mt;
#line 1163 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1164 "AB13MD.f"
	i__2 = iz10 + i__;
#line 1164 "AB13MD.f"
	i__3 = i__;
#line 1164 "AB13MD.f"
	z__1.r = x[i__3], z__1.i = 0.;
#line 1164 "AB13MD.f"
	zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
#line 1165 "AB13MD.f"
/* L820: */
#line 1165 "AB13MD.f"
    }
#line 1166 "AB13MD.f"
    i__1 = *n * *n;
#line 1166 "AB13MD.f"
    zcopy_(&i__1, &zwork[iz5 + 1], &c__1, &zwork[iz7 + 1], &c__1);
#line 1167 "AB13MD.f"
    i__1 = *n * *n;
#line 1167 "AB13MD.f"
    i__2 = *n * *n;
#line 1167 "AB13MD.f"
    zgemv_("N", &i__1, &mt, &c_b2, &zwork[iz6 + 1], &i__2, &zwork[iz10 + 1], &
	    c__1, &c_b2, &zwork[iz7 + 1], &c__1, (ftnlen)1);

/*           Compute diag( B ) = B0d + BBd*xD. */

#line 1172 "AB13MD.f"
    dcopy_(n, &dwork[iw7 + 1], &c__1, &dwork[iw24 + 1], &c__1);
#line 1173 "AB13MD.f"
    i__1 = *m - 1;
#line 1173 "AB13MD.f"
    dgemv_("N", n, &i__1, &c_b11, &dwork[iw8 + 1], n, &dwork[iw22 + 1], &c__1,
	     &c_b11, &dwork[iw24 + 1], &c__1, (ftnlen)1);

/*           Compute F. */

#line 1178 "AB13MD.f"
    i__1 = *n;
#line 1178 "AB13MD.f"
    for (j = 1; j <= i__1; ++j) {
#line 1179 "AB13MD.f"
	i__2 = *n;
#line 1179 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 1180 "AB13MD.f"
	    if (i__ == j) {
#line 1181 "AB13MD.f"
		i__3 = iz15 + i__ + (i__ - 1) * *n;
#line 1181 "AB13MD.f"
		d__1 = dlambd * dwork[iw24 + i__];
#line 1181 "AB13MD.f"
		z__2.r = d__1, z__2.i = 0.;
#line 1181 "AB13MD.f"
		i__4 = iz7 + i__ + (i__ - 1) * *n;
#line 1181 "AB13MD.f"
		z__1.r = z__2.r - zwork[i__4].r, z__1.i = z__2.i - zwork[i__4]
			.i;
#line 1181 "AB13MD.f"
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 1183 "AB13MD.f"
	    } else {
#line 1184 "AB13MD.f"
		i__3 = iz15 + i__ + (j - 1) * *n;
#line 1184 "AB13MD.f"
		i__4 = iz7 + i__ + (j - 1) * *n;
#line 1184 "AB13MD.f"
		z__1.r = -zwork[i__4].r, z__1.i = -zwork[i__4].i;
#line 1184 "AB13MD.f"
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 1185 "AB13MD.f"
	    }
#line 1186 "AB13MD.f"
/* L830: */
#line 1186 "AB13MD.f"
	}
#line 1187 "AB13MD.f"
/* L840: */
#line 1187 "AB13MD.f"
    }
#line 1188 "AB13MD.f"
    zlacpy_("Full", n, n, &zwork[iz15 + 1], n, &zwork[iz17 + 1], n, (ftnlen)4)
	    ;

/*           Compute det( F ). */

#line 1193 "AB13MD.f"
    i__1 = *lzwork - izwrk;
#line 1193 "AB13MD.f"
    zgees_("N", "N", (L_fp)select_, n, &zwork[iz15 + 1], n, &sdim, &zwork[
	    iz16 + 1], &zwork[1], n, &zwork[izwrk + 1], &i__1, &dwork[iwrk + 
	    1], bwork, &info2, (ftnlen)1, (ftnlen)1);
#line 1196 "AB13MD.f"
    if (info2 > 0) {
#line 1197 "AB13MD.f"
	*info = 6;
#line 1198 "AB13MD.f"
	return 0;
#line 1199 "AB13MD.f"
    }
#line 1200 "AB13MD.f"
    i__1 = izwrk + 1;
#line 1200 "AB13MD.f"
    lza = (integer) zwork[i__1].r;
#line 1201 "AB13MD.f"
    lzamax = max(lza,lzamax);
#line 1202 "AB13MD.f"
    detf.r = 1., detf.i = 0.;
#line 1203 "AB13MD.f"
    i__1 = *n;
#line 1203 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1204 "AB13MD.f"
	i__2 = iz16 + i__;
#line 1204 "AB13MD.f"
	z__1.r = detf.r * zwork[i__2].r - detf.i * zwork[i__2].i, z__1.i = 
		detf.r * zwork[i__2].i + detf.i * zwork[i__2].r;
#line 1204 "AB13MD.f"
	detf.r = z__1.r, detf.i = z__1.i;
#line 1205 "AB13MD.f"
/* L850: */
#line 1205 "AB13MD.f"
    }

/*           Compute Finv. */

#line 1209 "AB13MD.f"
    zgetrf_(n, n, &zwork[iz17 + 1], n, &iwork[1], &info2);
#line 1210 "AB13MD.f"
    if (info2 > 0) {
#line 1211 "AB13MD.f"
	*info = 5;
#line 1212 "AB13MD.f"
	return 0;
#line 1213 "AB13MD.f"
    }
#line 1214 "AB13MD.f"
    i__1 = *ldwork - iwrk;
#line 1214 "AB13MD.f"
    zgetri_(n, &zwork[iz17 + 1], n, &iwork[1], &zwork[izwrk + 1], &i__1, &
	    info2);
#line 1216 "AB13MD.f"
    i__1 = izwrk + 1;
#line 1216 "AB13MD.f"
    lza = (integer) zwork[i__1].r;
#line 1217 "AB13MD.f"
    lzamax = max(lza,lzamax);

/*           Compute the barrier function. */

#line 1221 "AB13MD.f"
    i__1 = *m - 1;
#line 1221 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1222 "AB13MD.f"
	dwork[iw25 + i__] = dwork[iw22 + i__] - .01;
#line 1223 "AB13MD.f"
	dwork[iw25 + *m - 1 + i__] = 100. - dwork[iw22 + i__];
#line 1224 "AB13MD.f"
/* L860: */
#line 1224 "AB13MD.f"
    }
#line 1225 "AB13MD.f"
    if (mr > 0) {
#line 1226 "AB13MD.f"
	i__1 = mr;
#line 1226 "AB13MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 1227 "AB13MD.f"
	    dwork[iw25 + (*m - 1 << 1) + i__] = dwork[iw23 + i__] + tau;
#line 1228 "AB13MD.f"
	    dwork[iw25 + (*m - 1 << 1) + mr + i__] = tau - dwork[iw23 + i__];
#line 1229 "AB13MD.f"
/* L870: */
#line 1229 "AB13MD.f"
	}
#line 1230 "AB13MD.f"
    }
#line 1231 "AB13MD.f"
    prod = 1.;
#line 1232 "AB13MD.f"
    i__1 = mt << 1;
#line 1232 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1233 "AB13MD.f"
	prod *= dwork[iw25 + i__];
#line 1234 "AB13MD.f"
/* L880: */
#line 1234 "AB13MD.f"
    }
#line 1235 "AB13MD.f"
    temp = detf.r;
#line 1236 "AB13MD.f"
    if (temp < eps) {
#line 1236 "AB13MD.f"
	temp = eps;
#line 1236 "AB13MD.f"
    }
#line 1237 "AB13MD.f"
    phi = -log(temp) - log(prod);

/*           Compute the gradient of the barrier function. */

#line 1241 "AB13MD.f"
    i__1 = mt;
#line 1241 "AB13MD.f"
    for (j = 1; j <= i__1; ++j) {
#line 1242 "AB13MD.f"
	i__2 = *n * *n;
#line 1242 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 1243 "AB13MD.f"
	    i__3 = iz18 + i__ + (j - 1) * *n * *n;
#line 1243 "AB13MD.f"
	    d__1 = dlambd * dwork[iw9 + i__ + (j - 1) * *n * *n];
#line 1243 "AB13MD.f"
	    z__2.r = d__1, z__2.i = 0.;
#line 1243 "AB13MD.f"
	    i__4 = iz6 + i__ + (j - 1) * *n * *n;
#line 1243 "AB13MD.f"
	    z__1.r = z__2.r - zwork[i__4].r, z__1.i = z__2.i - zwork[i__4].i;
#line 1243 "AB13MD.f"
	    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 1245 "AB13MD.f"
/* L890: */
#line 1245 "AB13MD.f"
	}
#line 1246 "AB13MD.f"
/* L900: */
#line 1246 "AB13MD.f"
    }
#line 1247 "AB13MD.f"
    i__1 = *n * *n;
#line 1247 "AB13MD.f"
    i__2 = *n * *n;
#line 1247 "AB13MD.f"
    zgemv_("C", &i__1, &mt, &c_b2, &zwork[iz18 + 1], &i__2, &zwork[iz17 + 1], 
	    &c__1, &c_b1, &zwork[iz19 + 1], &c__1, (ftnlen)1);
#line 1249 "AB13MD.f"
    i__1 = *m - 1;
#line 1249 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1250 "AB13MD.f"
	dwork[iw26 + i__] = 1. / (dwork[iw22 + i__] - .01) - 1. / (100. - 
		dwork[iw22 + i__]);
#line 1252 "AB13MD.f"
/* L910: */
#line 1252 "AB13MD.f"
    }
#line 1253 "AB13MD.f"
    if (mr > 0) {
#line 1254 "AB13MD.f"
	i__1 = mr;
#line 1254 "AB13MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 1255 "AB13MD.f"
	    dwork[iw26 + *m - 1 + i__] = 1. / (dwork[iw23 + i__] + tau) - 1. /
		     (tau - dwork[iw23 + i__]);
#line 1257 "AB13MD.f"
/* L920: */
#line 1257 "AB13MD.f"
	}
#line 1258 "AB13MD.f"
    }
#line 1259 "AB13MD.f"
    i__1 = mt;
#line 1259 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1260 "AB13MD.f"
	i__2 = iz19 + i__;
#line 1260 "AB13MD.f"
	dwork[iw26 + i__] = -zwork[i__2].r - dwork[iw26 + i__];
#line 1262 "AB13MD.f"
/* L925: */
#line 1262 "AB13MD.f"
    }

/*           Compute the Hessian of the barrier function. */

#line 1266 "AB13MD.f"
    i__1 = *n * mt;
#line 1266 "AB13MD.f"
    zgemm_("N", "N", n, &i__1, n, &c_b2, &zwork[iz17 + 1], n, &zwork[iz18 + 1]
	    , n, &c_b1, &zwork[iz20 + 1], n, (ftnlen)1, (ftnlen)1);
#line 1269 "AB13MD.f"
    dlaset_("Full", &mt, &mt, &c_b15, &c_b15, &dwork[iw11 + 1], &mt, (ftnlen)
	    4);
#line 1271 "AB13MD.f"
    i__1 = mt;
#line 1271 "AB13MD.f"
    for (k = 1; k <= i__1; ++k) {
#line 1272 "AB13MD.f"
	i__2 = *n * *n;
#line 1272 "AB13MD.f"
	zcopy_(&i__2, &zwork[iz20 + 1 + (k - 1) * *n * *n], &c__1, &zwork[
		iz22 + 1], &c__1);
#line 1274 "AB13MD.f"
	i__2 = *n;
#line 1274 "AB13MD.f"
	for (j = 1; j <= i__2; ++j) {
#line 1275 "AB13MD.f"
	    i__3 = *n;
#line 1275 "AB13MD.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1276 "AB13MD.f"
		i__4 = iz23 + i__ + (j - 1) * *n;
#line 1276 "AB13MD.f"
		d_cnjg(&z__1, &zwork[iz22 + j + (i__ - 1) * *n]);
#line 1276 "AB13MD.f"
		zwork[i__4].r = z__1.r, zwork[i__4].i = z__1.i;
#line 1278 "AB13MD.f"
/* L930: */
#line 1278 "AB13MD.f"
	    }
#line 1279 "AB13MD.f"
/* L940: */
#line 1279 "AB13MD.f"
	}
#line 1280 "AB13MD.f"
	i__2 = *n * *n;
#line 1280 "AB13MD.f"
	i__3 = *n * *n;
#line 1280 "AB13MD.f"
	zgemv_("C", &i__2, &k, &c_b2, &zwork[iz20 + 1], &i__3, &zwork[iz23 + 
		1], &c__1, &c_b1, &zwork[iz24 + 1], &c__1, (ftnlen)1);
#line 1283 "AB13MD.f"
	i__2 = k;
#line 1283 "AB13MD.f"
	for (j = 1; j <= i__2; ++j) {
#line 1284 "AB13MD.f"
	    d_cnjg(&z__1, &zwork[iz24 + j]);
#line 1284 "AB13MD.f"
	    dwork[iw11 + k + (j - 1) * mt] = z__1.r;
#line 1286 "AB13MD.f"
/* L950: */
#line 1286 "AB13MD.f"
	}
#line 1287 "AB13MD.f"
/* L960: */
#line 1287 "AB13MD.f"
    }
#line 1288 "AB13MD.f"
    i__1 = *m - 1;
#line 1288 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
#line 1289 "AB13MD.f"
	d__1 = dwork[iw22 + i__] - .01;
/* Computing 2nd power */
#line 1289 "AB13MD.f"
	d__2 = 100. - dwork[iw22 + i__];
#line 1289 "AB13MD.f"
	dwork[iw10 + i__] = 1. / (d__1 * d__1) + 1. / (d__2 * d__2);
#line 1291 "AB13MD.f"
/* L970: */
#line 1291 "AB13MD.f"
    }
#line 1292 "AB13MD.f"
    if (mr > 0) {
#line 1293 "AB13MD.f"
	i__1 = mr;
#line 1293 "AB13MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
#line 1294 "AB13MD.f"
	    d__1 = dwork[iw23 + i__] + tau;
/* Computing 2nd power */
#line 1294 "AB13MD.f"
	    d__2 = tau - dwork[iw23 + i__];
#line 1294 "AB13MD.f"
	    dwork[iw10 + *m - 1 + i__] = 1. / (d__1 * d__1) + 1. / (d__2 * 
		    d__2);
#line 1297 "AB13MD.f"
/* L980: */
#line 1297 "AB13MD.f"
	}
#line 1298 "AB13MD.f"
    }
#line 1299 "AB13MD.f"
    i__1 = mt;
#line 1299 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1300 "AB13MD.f"
	dwork[iw11 + i__ + (i__ - 1) * mt] += dwork[iw10 + i__];
#line 1302 "AB13MD.f"
/* L990: */
#line 1302 "AB13MD.f"
    }
#line 1303 "AB13MD.f"
    i__1 = mt;
#line 1303 "AB13MD.f"
    for (j = 1; j <= i__1; ++j) {
#line 1304 "AB13MD.f"
	i__2 = j;
#line 1304 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 1305 "AB13MD.f"
	    if (i__ != j) {
#line 1306 "AB13MD.f"
		t1 = dwork[iw11 + i__ + (j - 1) * mt];
#line 1307 "AB13MD.f"
		t2 = dwork[iw11 + j + (i__ - 1) * mt];
#line 1308 "AB13MD.f"
		dwork[iw11 + i__ + (j - 1) * mt] = t1 + t2;
#line 1309 "AB13MD.f"
		dwork[iw11 + j + (i__ - 1) * mt] = t1 + t2;
#line 1310 "AB13MD.f"
	    }
#line 1311 "AB13MD.f"
/* L1000: */
#line 1311 "AB13MD.f"
	}
#line 1312 "AB13MD.f"
/* L1100: */
#line 1312 "AB13MD.f"
    }

/*           Compute norm( H ). */

#line 1316 "AB13MD.f"
L1110:
#line 1316 "AB13MD.f"
    hnorm = dlange_("F", &mt, &mt, &dwork[iw11 + 1], &mt, &dwork[1], (ftnlen)
	    1);

/*           Compute rcond( H ). */

#line 1320 "AB13MD.f"
    dlacpy_("Full", &mt, &mt, &dwork[iw11 + 1], &mt, &dwork[iw31 + 1], &mt, (
	    ftnlen)4);
#line 1322 "AB13MD.f"
    hnorm1 = dlange_("1", &mt, &mt, &dwork[iw31 + 1], &mt, &dwork[1], (ftnlen)
	    1);
#line 1323 "AB13MD.f"
    i__1 = *ldwork - iwrk;
#line 1323 "AB13MD.f"
    dsytrf_("U", &mt, &dwork[iw31 + 1], &mt, &iwork[1], &dwork[iwrk + 1], &
	    i__1, &info2, (ftnlen)1);
#line 1325 "AB13MD.f"
    if (info2 > 0) {
#line 1326 "AB13MD.f"
	*info = 5;
#line 1327 "AB13MD.f"
	return 0;
#line 1328 "AB13MD.f"
    }
#line 1329 "AB13MD.f"
    lwa = (integer) dwork[iwrk + 1];
#line 1330 "AB13MD.f"
    lwamax = max(lwa,lwamax);
#line 1331 "AB13MD.f"
    dsycon_("U", &mt, &dwork[iw31 + 1], &mt, &iwork[1], &hnorm1, &rcond, &
	    dwork[iwrk + 1], &iwork[mt + 1], &info2, (ftnlen)1);
#line 1333 "AB13MD.f"
    if (rcond < tol3) {
#line 1334 "AB13MD.f"
	i__1 = mt;
#line 1334 "AB13MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 1335 "AB13MD.f"
	    dwork[iw11 + i__ + (i__ - 1) * mt] += hnorm * regpar;
#line 1337 "AB13MD.f"
/* L1120: */
#line 1337 "AB13MD.f"
	}
#line 1338 "AB13MD.f"
	goto L1110;
#line 1339 "AB13MD.f"
    }

/*           Compute the tangent line to path of center. */

#line 1343 "AB13MD.f"
    dcopy_(&mt, &dwork[iw26 + 1], &c__1, &dwork[iw27 + 1], &c__1);
#line 1344 "AB13MD.f"
    dsytrs_("U", &mt, &c__1, &dwork[iw31 + 1], &mt, &iwork[1], &dwork[iw27 + 
	    1], &mt, &info2, (ftnlen)1);

/*           Check if x-h satisfies the Goldstein test. */

#line 1349 "AB13MD.f"
    gtest = FALSE_;
#line 1350 "AB13MD.f"
    i__1 = mt;
#line 1350 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1351 "AB13MD.f"
	dwork[iw20 + i__] = x[i__] - dwork[iw27 + i__];
#line 1352 "AB13MD.f"
/* L1130: */
#line 1352 "AB13MD.f"
    }

/*           Store xD. */

#line 1356 "AB13MD.f"
    i__1 = *m - 1;
#line 1356 "AB13MD.f"
    dcopy_(&i__1, &dwork[iw20 + 1], &c__1, &dwork[iw22 + 1], &c__1);
#line 1357 "AB13MD.f"
    if (mr > 0) {

/*              Store xG. */

#line 1361 "AB13MD.f"
	dcopy_(&mr, &dwork[iw20 + *m], &c__1, &dwork[iw23 + 1], &c__1);
#line 1362 "AB13MD.f"
    }

/*           Compute A(:) = A0 + AA*x_new. */

#line 1366 "AB13MD.f"
    i__1 = mt;
#line 1366 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1367 "AB13MD.f"
	i__2 = iz10 + i__;
#line 1367 "AB13MD.f"
	i__3 = iw20 + i__;
#line 1367 "AB13MD.f"
	z__1.r = dwork[i__3], z__1.i = 0.;
#line 1367 "AB13MD.f"
	zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
#line 1368 "AB13MD.f"
/* L1140: */
#line 1368 "AB13MD.f"
    }
#line 1369 "AB13MD.f"
    i__1 = *n * *n;
#line 1369 "AB13MD.f"
    zcopy_(&i__1, &zwork[iz5 + 1], &c__1, &zwork[iz7 + 1], &c__1);
#line 1370 "AB13MD.f"
    i__1 = *n * *n;
#line 1370 "AB13MD.f"
    i__2 = *n * *n;
#line 1370 "AB13MD.f"
    zgemv_("N", &i__1, &mt, &c_b2, &zwork[iz6 + 1], &i__2, &zwork[iz10 + 1], &
	    c__1, &c_b2, &zwork[iz7 + 1], &c__1, (ftnlen)1);

/*           Compute diag( B ) = B0d + BBd*xD. */

#line 1375 "AB13MD.f"
    dcopy_(n, &dwork[iw7 + 1], &c__1, &dwork[iw24 + 1], &c__1);
#line 1376 "AB13MD.f"
    i__1 = *m - 1;
#line 1376 "AB13MD.f"
    dgemv_("N", n, &i__1, &c_b11, &dwork[iw8 + 1], n, &dwork[iw22 + 1], &c__1,
	     &c_b11, &dwork[iw24 + 1], &c__1, (ftnlen)1);

/*           Compute lambda*diag(B) - A. */

#line 1381 "AB13MD.f"
    i__1 = *n;
#line 1381 "AB13MD.f"
    for (j = 1; j <= i__1; ++j) {
#line 1382 "AB13MD.f"
	i__2 = *n;
#line 1382 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 1383 "AB13MD.f"
	    if (i__ == j) {
#line 1384 "AB13MD.f"
		i__3 = iz15 + i__ + (i__ - 1) * *n;
#line 1384 "AB13MD.f"
		d__1 = dlambd * dwork[iw24 + i__];
#line 1384 "AB13MD.f"
		z__2.r = d__1, z__2.i = 0.;
#line 1384 "AB13MD.f"
		i__4 = iz7 + i__ + (i__ - 1) * *n;
#line 1384 "AB13MD.f"
		z__1.r = z__2.r - zwork[i__4].r, z__1.i = z__2.i - zwork[i__4]
			.i;
#line 1384 "AB13MD.f"
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 1386 "AB13MD.f"
	    } else {
#line 1387 "AB13MD.f"
		i__3 = iz15 + i__ + (j - 1) * *n;
#line 1387 "AB13MD.f"
		i__4 = iz7 + i__ + (j - 1) * *n;
#line 1387 "AB13MD.f"
		z__1.r = -zwork[i__4].r, z__1.i = -zwork[i__4].i;
#line 1387 "AB13MD.f"
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 1388 "AB13MD.f"
	    }
#line 1389 "AB13MD.f"
/* L1150: */
#line 1389 "AB13MD.f"
	}
#line 1390 "AB13MD.f"
/* L1160: */
#line 1390 "AB13MD.f"
    }

/*           Compute eig( lambda*diag(B)-A ). */

#line 1394 "AB13MD.f"
    i__1 = *lzwork - izwrk;
#line 1394 "AB13MD.f"
    zgees_("N", "N", (L_fp)select_, n, &zwork[iz15 + 1], n, &sdim, &zwork[
	    iz16 + 1], &zwork[1], n, &zwork[izwrk + 1], &i__1, &dwork[iwrk + 
	    1], bwork, &info2, (ftnlen)1, (ftnlen)1);
#line 1397 "AB13MD.f"
    if (info2 > 0) {
#line 1398 "AB13MD.f"
	*info = 6;
#line 1399 "AB13MD.f"
	return 0;
#line 1400 "AB13MD.f"
    }
#line 1401 "AB13MD.f"
    i__1 = izwrk + 1;
#line 1401 "AB13MD.f"
    lza = (integer) zwork[i__1].r;
#line 1402 "AB13MD.f"
    lzamax = max(lza,lzamax);
#line 1403 "AB13MD.f"
    i__1 = *n;
#line 1403 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1404 "AB13MD.f"
	i__2 = iz16 + i__;
#line 1404 "AB13MD.f"
	dwork[iw30 + i__] = zwork[i__2].r;
#line 1405 "AB13MD.f"
/* L1190: */
#line 1405 "AB13MD.f"
    }
#line 1406 "AB13MD.f"
    i__1 = *m - 1;
#line 1406 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1407 "AB13MD.f"
	dwork[iw30 + *n + i__] = dwork[iw22 + i__] - .01;
#line 1408 "AB13MD.f"
	dwork[iw30 + *n + *m - 1 + i__] = 100. - dwork[iw22 + i__];
#line 1409 "AB13MD.f"
/* L1200: */
#line 1409 "AB13MD.f"
    }
#line 1410 "AB13MD.f"
    if (mr > 0) {
#line 1411 "AB13MD.f"
	i__1 = mr;
#line 1411 "AB13MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 1412 "AB13MD.f"
	    dwork[iw30 + *n + (*m - 1 << 1) + i__] = dwork[iw23 + i__] + tau;
#line 1413 "AB13MD.f"
	    dwork[iw30 + *n + (*m - 1 << 1) + mr + i__] = tau - dwork[iw23 + 
		    i__];
#line 1414 "AB13MD.f"
/* L1210: */
#line 1414 "AB13MD.f"
	}
#line 1415 "AB13MD.f"
    }
#line 1416 "AB13MD.f"
    emin = dwork[iw30 + 1];
#line 1417 "AB13MD.f"
    i__1 = *n + (mt << 1);
#line 1417 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1418 "AB13MD.f"
	if (dwork[iw30 + i__] < emin) {
#line 1418 "AB13MD.f"
	    emin = dwork[iw30 + i__];
#line 1418 "AB13MD.f"
	}
#line 1419 "AB13MD.f"
/* L1220: */
#line 1419 "AB13MD.f"
    }
#line 1420 "AB13MD.f"
    if (emin <= 0.) {
#line 1421 "AB13MD.f"
	gtest = FALSE_;
#line 1422 "AB13MD.f"
    } else {
#line 1423 "AB13MD.f"
	pp = ddot_(&mt, &dwork[iw26 + 1], &c__1, &dwork[iw27 + 1], &c__1);
#line 1424 "AB13MD.f"
	prod = 1.;
#line 1425 "AB13MD.f"
	i__1 = *n + (mt << 1);
#line 1425 "AB13MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 1426 "AB13MD.f"
	    prod *= dwork[iw30 + i__];
#line 1427 "AB13MD.f"
/* L1230: */
#line 1427 "AB13MD.f"
	}
#line 1428 "AB13MD.f"
	t1 = -log(prod);
#line 1429 "AB13MD.f"
	t2 = phi - pp * .01;
#line 1430 "AB13MD.f"
	t3 = phi - pp * .9;
#line 1431 "AB13MD.f"
	if (t1 >= t3 && t1 < t2) {
#line 1431 "AB13MD.f"
	    gtest = TRUE_;
#line 1431 "AB13MD.f"
	}
#line 1432 "AB13MD.f"
    }

/*           Use x-h if Goldstein test is satisfied. Otherwise use */
/*           Nesterov-Nemirovsky's stepsize length. */

#line 1437 "AB13MD.f"
    pp = ddot_(&mt, &dwork[iw26 + 1], &c__1, &dwork[iw27 + 1], &c__1);
#line 1438 "AB13MD.f"
    delta = sqrt(pp);
#line 1439 "AB13MD.f"
    if (gtest || delta <= .25) {
#line 1440 "AB13MD.f"
	i__1 = mt;
#line 1440 "AB13MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 1441 "AB13MD.f"
	    x[i__] -= dwork[iw27 + i__];
#line 1442 "AB13MD.f"
/* L1240: */
#line 1442 "AB13MD.f"
	}
#line 1443 "AB13MD.f"
    } else {
#line 1444 "AB13MD.f"
	i__1 = mt;
#line 1444 "AB13MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 1445 "AB13MD.f"
	    x[i__] -= dwork[iw27 + i__] / (delta + 1.);
#line 1446 "AB13MD.f"
/* L1250: */
#line 1446 "AB13MD.f"
	}
#line 1447 "AB13MD.f"
    }

/*           Analytic center is found if delta is sufficiently small. */

#line 1451 "AB13MD.f"
    if (delta < tol5) {
#line 1451 "AB13MD.f"
	goto L1260;
#line 1451 "AB13MD.f"
    }
#line 1452 "AB13MD.f"
    goto L810;

/*        Set yf. */

#line 1456 "AB13MD.f"
L1260:
#line 1456 "AB13MD.f"
    dwork[iw14 + 1] = dlambd;
#line 1457 "AB13MD.f"
    dcopy_(&mt, &x[1], &c__1, &dwork[iw14 + 2], &c__1);

/*        Set yw. */

#line 1461 "AB13MD.f"
    i__1 = mt + 1;
#line 1461 "AB13MD.f"
    dcopy_(&i__1, &dwork[iw14 + 1], &c__1, &dwork[iw15 + 1], &c__1);

/*        Compute Fb. */

#line 1465 "AB13MD.f"
    i__1 = *n;
#line 1465 "AB13MD.f"
    for (j = 1; j <= i__1; ++j) {
#line 1466 "AB13MD.f"
	i__2 = *n;
#line 1466 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 1467 "AB13MD.f"
	    i__3 = iz21 + i__ + (j - 1) * *n;
#line 1467 "AB13MD.f"
	    i__4 = iw24 + i__;
#line 1467 "AB13MD.f"
	    z__2.r = dwork[i__4], z__2.i = 0.;
#line 1467 "AB13MD.f"
	    d_cnjg(&z__3, &zwork[iz17 + j + (i__ - 1) * *n]);
#line 1467 "AB13MD.f"
	    z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = z__2.r * 
		    z__3.i + z__2.i * z__3.r;
#line 1467 "AB13MD.f"
	    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 1469 "AB13MD.f"
/* L1270: */
#line 1469 "AB13MD.f"
	}
#line 1470 "AB13MD.f"
/* L1280: */
#line 1470 "AB13MD.f"
    }
#line 1471 "AB13MD.f"
    i__1 = *n * *n;
#line 1471 "AB13MD.f"
    i__2 = *n * *n;
#line 1471 "AB13MD.f"
    zgemv_("C", &i__1, &mt, &c_b2, &zwork[iz20 + 1], &i__2, &zwork[iz21 + 1], 
	    &c__1, &c_b1, &zwork[iz24 + 1], &c__1, (ftnlen)1);
#line 1473 "AB13MD.f"
    i__1 = mt;
#line 1473 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1474 "AB13MD.f"
	i__2 = iz24 + i__;
#line 1474 "AB13MD.f"
	dwork[iw32 + i__] = zwork[i__2].r;
#line 1475 "AB13MD.f"
/* L1300: */
#line 1475 "AB13MD.f"
    }

/*        Compute h1. */

#line 1479 "AB13MD.f"
    dlacpy_("Full", &mt, &mt, &dwork[iw11 + 1], &mt, &dwork[iw31 + 1], &mt, (
	    ftnlen)4);
#line 1481 "AB13MD.f"
    i__1 = *ldwork - iwrk;
#line 1481 "AB13MD.f"
    dsysv_("U", &mt, &c__1, &dwork[iw31 + 1], &mt, &iwork[1], &dwork[iw32 + 1]
	    , &mt, &dwork[iwrk + 1], &i__1, &info2, (ftnlen)1);
#line 1484 "AB13MD.f"
    if (info2 > 0) {
#line 1485 "AB13MD.f"
	*info = 5;
#line 1486 "AB13MD.f"
	return 0;
#line 1487 "AB13MD.f"
    }
#line 1488 "AB13MD.f"
    lwa = (integer) dwork[iwrk + 1];
#line 1489 "AB13MD.f"
    lwamax = max(lwa,lwamax);

/*        Compute hn. */

#line 1493 "AB13MD.f"
    hn = dlange_("F", &mt, &c__1, &dwork[iw32 + 1], &mt, &dwork[1], (ftnlen)1)
	    ;

/*        Compute y. */

#line 1497 "AB13MD.f"
    dwork[iw13 + 1] = dlambd - c__ / hn;
#line 1498 "AB13MD.f"
    i__1 = mt;
#line 1498 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1499 "AB13MD.f"
	dwork[iw13 + 1 + i__] = x[i__] + c__ * dwork[iw32 + i__] / hn;
#line 1500 "AB13MD.f"
/* L1310: */
#line 1500 "AB13MD.f"
    }

/*        Store xD. */

#line 1504 "AB13MD.f"
    i__1 = *m - 1;
#line 1504 "AB13MD.f"
    dcopy_(&i__1, &dwork[iw13 + 2], &c__1, &dwork[iw22 + 1], &c__1);
#line 1505 "AB13MD.f"
    if (mr > 0) {

/*           Store xG. */

#line 1509 "AB13MD.f"
	dcopy_(&mr, &dwork[iw13 + *m + 1], &c__1, &dwork[iw23 + 1], &c__1);
#line 1510 "AB13MD.f"
    }

/*        Compute A(:) = A0 + AA*y(2:mt+1). */

#line 1514 "AB13MD.f"
    i__1 = mt;
#line 1514 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1515 "AB13MD.f"
	i__2 = iz10 + i__;
#line 1515 "AB13MD.f"
	i__3 = iw13 + 1 + i__;
#line 1515 "AB13MD.f"
	z__1.r = dwork[i__3], z__1.i = 0.;
#line 1515 "AB13MD.f"
	zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
#line 1516 "AB13MD.f"
/* L1320: */
#line 1516 "AB13MD.f"
    }
#line 1517 "AB13MD.f"
    i__1 = *n * *n;
#line 1517 "AB13MD.f"
    zcopy_(&i__1, &zwork[iz5 + 1], &c__1, &zwork[iz7 + 1], &c__1);
#line 1518 "AB13MD.f"
    i__1 = *n * *n;
#line 1518 "AB13MD.f"
    i__2 = *n * *n;
#line 1518 "AB13MD.f"
    zgemv_("N", &i__1, &mt, &c_b2, &zwork[iz6 + 1], &i__2, &zwork[iz10 + 1], &
	    c__1, &c_b2, &zwork[iz7 + 1], &c__1, (ftnlen)1);

/*        Compute B = B0d + BBd*xD. */

#line 1523 "AB13MD.f"
    dcopy_(n, &dwork[iw7 + 1], &c__1, &dwork[iw24 + 1], &c__1);
#line 1524 "AB13MD.f"
    i__1 = *m - 1;
#line 1524 "AB13MD.f"
    dgemv_("N", n, &i__1, &c_b11, &dwork[iw8 + 1], n, &dwork[iw22 + 1], &c__1,
	     &c_b11, &dwork[iw24 + 1], &c__1, (ftnlen)1);

/*        Compute y(1)*diag(B) - A. */

#line 1529 "AB13MD.f"
    i__1 = *n;
#line 1529 "AB13MD.f"
    for (j = 1; j <= i__1; ++j) {
#line 1530 "AB13MD.f"
	i__2 = *n;
#line 1530 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 1531 "AB13MD.f"
	    if (i__ == j) {
#line 1532 "AB13MD.f"
		i__3 = iz15 + i__ + (i__ - 1) * *n;
#line 1532 "AB13MD.f"
		d__1 = dwork[iw13 + 1] * dwork[iw24 + i__];
#line 1532 "AB13MD.f"
		z__2.r = d__1, z__2.i = 0.;
#line 1532 "AB13MD.f"
		i__4 = iz7 + i__ + (i__ - 1) * *n;
#line 1532 "AB13MD.f"
		z__1.r = z__2.r - zwork[i__4].r, z__1.i = z__2.i - zwork[i__4]
			.i;
#line 1532 "AB13MD.f"
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 1534 "AB13MD.f"
	    } else {
#line 1535 "AB13MD.f"
		i__3 = iz15 + i__ + (j - 1) * *n;
#line 1535 "AB13MD.f"
		i__4 = iz7 + i__ + (j - 1) * *n;
#line 1535 "AB13MD.f"
		z__1.r = -zwork[i__4].r, z__1.i = -zwork[i__4].i;
#line 1535 "AB13MD.f"
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 1536 "AB13MD.f"
	    }
#line 1537 "AB13MD.f"
/* L1330: */
#line 1537 "AB13MD.f"
	}
#line 1538 "AB13MD.f"
/* L1340: */
#line 1538 "AB13MD.f"
    }

/*        Compute eig( y(1)*diag(B)-A ). */

#line 1542 "AB13MD.f"
    i__1 = *lzwork - izwrk;
#line 1542 "AB13MD.f"
    zgees_("N", "N", (L_fp)select_, n, &zwork[iz15 + 1], n, &sdim, &zwork[
	    iz16 + 1], &zwork[1], n, &zwork[izwrk + 1], &i__1, &dwork[iwrk + 
	    1], bwork, &info2, (ftnlen)1, (ftnlen)1);
#line 1545 "AB13MD.f"
    if (info2 > 0) {
#line 1546 "AB13MD.f"
	*info = 6;
#line 1547 "AB13MD.f"
	return 0;
#line 1548 "AB13MD.f"
    }
#line 1549 "AB13MD.f"
    i__1 = izwrk + 1;
#line 1549 "AB13MD.f"
    lza = (integer) zwork[i__1].r;
#line 1550 "AB13MD.f"
    lzamax = max(lza,lzamax);
#line 1551 "AB13MD.f"
    i__1 = iz16 + 1;
#line 1551 "AB13MD.f"
    emin = zwork[i__1].r;
#line 1552 "AB13MD.f"
    if (*n > 1) {
#line 1553 "AB13MD.f"
	i__1 = *n;
#line 1553 "AB13MD.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 1554 "AB13MD.f"
	    i__2 = iz16 + i__;
#line 1554 "AB13MD.f"
	    if (zwork[i__2].r < emin) {
#line 1554 "AB13MD.f"
		i__3 = iz16 + i__;
#line 1554 "AB13MD.f"
		emin = zwork[i__3].r;
#line 1554 "AB13MD.f"
	    }
#line 1556 "AB13MD.f"
/* L1350: */
#line 1556 "AB13MD.f"
	}
#line 1557 "AB13MD.f"
    }
#line 1558 "AB13MD.f"
    pos = TRUE_;
#line 1559 "AB13MD.f"
    i__1 = *m - 1;
#line 1559 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1560 "AB13MD.f"
	dwork[iw25 + i__] = dwork[iw22 + i__] - .01;
#line 1561 "AB13MD.f"
	dwork[iw25 + *m - 1 + i__] = 100. - dwork[iw22 + i__];
#line 1562 "AB13MD.f"
/* L1360: */
#line 1562 "AB13MD.f"
    }
#line 1563 "AB13MD.f"
    if (mr > 0) {
#line 1564 "AB13MD.f"
	i__1 = mr;
#line 1564 "AB13MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 1565 "AB13MD.f"
	    dwork[iw25 + (*m - 1 << 1) + i__] = dwork[iw23 + i__] + tau;
#line 1566 "AB13MD.f"
	    dwork[iw25 + (*m - 1 << 1) + mr + i__] = tau - dwork[iw23 + i__];
#line 1567 "AB13MD.f"
/* L1370: */
#line 1567 "AB13MD.f"
	}
#line 1568 "AB13MD.f"
    }
#line 1569 "AB13MD.f"
    temp = dwork[iw25 + 1];
#line 1570 "AB13MD.f"
    i__1 = mt << 1;
#line 1570 "AB13MD.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 1571 "AB13MD.f"
	if (dwork[iw25 + i__] < temp) {
#line 1571 "AB13MD.f"
	    temp = dwork[iw25 + i__];
#line 1571 "AB13MD.f"
	}
#line 1572 "AB13MD.f"
/* L1380: */
#line 1572 "AB13MD.f"
    }
#line 1573 "AB13MD.f"
    if (temp <= 0. || emin <= 0.) {
#line 1573 "AB13MD.f"
	pos = FALSE_;
#line 1573 "AB13MD.f"
    }
#line 1574 "AB13MD.f"
L1390:
#line 1574 "AB13MD.f"
    if (pos) {

/*           Set y2 = y. */

#line 1578 "AB13MD.f"
	i__1 = mt + 1;
#line 1578 "AB13MD.f"
	dcopy_(&i__1, &dwork[iw13 + 1], &c__1, &dwork[iw17 + 1], &c__1);

/*           Compute y = y + 1.5*( y - yw ). */

#line 1582 "AB13MD.f"
	i__1 = mt + 1;
#line 1582 "AB13MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 1583 "AB13MD.f"
	    dwork[iw13 + i__] += (dwork[iw13 + i__] - dwork[iw15 + i__]) * 
		    1.5;
#line 1585 "AB13MD.f"
/* L1400: */
#line 1585 "AB13MD.f"
	}

/*           Store xD. */

#line 1589 "AB13MD.f"
	i__1 = *m - 1;
#line 1589 "AB13MD.f"
	dcopy_(&i__1, &dwork[iw13 + 2], &c__1, &dwork[iw22 + 1], &c__1);
#line 1590 "AB13MD.f"
	if (mr > 0) {

/*              Store xG. */

#line 1594 "AB13MD.f"
	    dcopy_(&mr, &dwork[iw13 + *m + 1], &c__1, &dwork[iw23 + 1], &c__1)
		    ;
#line 1596 "AB13MD.f"
	}

/*           Compute A(:) = A0 + AA*y(2:mt+1). */

#line 1600 "AB13MD.f"
	i__1 = mt;
#line 1600 "AB13MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 1601 "AB13MD.f"
	    i__2 = iz10 + i__;
#line 1601 "AB13MD.f"
	    i__3 = iw13 + 1 + i__;
#line 1601 "AB13MD.f"
	    z__1.r = dwork[i__3], z__1.i = 0.;
#line 1601 "AB13MD.f"
	    zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
#line 1602 "AB13MD.f"
/* L1420: */
#line 1602 "AB13MD.f"
	}
#line 1603 "AB13MD.f"
	i__1 = *n * *n;
#line 1603 "AB13MD.f"
	zcopy_(&i__1, &zwork[iz5 + 1], &c__1, &zwork[iz7 + 1], &c__1);
#line 1604 "AB13MD.f"
	i__1 = *n * *n;
#line 1604 "AB13MD.f"
	i__2 = *n * *n;
#line 1604 "AB13MD.f"
	zgemv_("N", &i__1, &mt, &c_b2, &zwork[iz6 + 1], &i__2, &zwork[iz10 + 
		1], &c__1, &c_b2, &zwork[iz7 + 1], &c__1, (ftnlen)1);

/*           Compute diag( B ) = B0d + BBd*xD. */

#line 1609 "AB13MD.f"
	dcopy_(n, &dwork[iw7 + 1], &c__1, &dwork[iw24 + 1], &c__1);
#line 1610 "AB13MD.f"
	i__1 = *m - 1;
#line 1610 "AB13MD.f"
	dgemv_("N", n, &i__1, &c_b11, &dwork[iw8 + 1], n, &dwork[iw22 + 1], &
		c__1, &c_b11, &dwork[iw24 + 1], &c__1, (ftnlen)1);

/*           Set yw = y2. */

#line 1615 "AB13MD.f"
	i__1 = mt + 1;
#line 1615 "AB13MD.f"
	dcopy_(&i__1, &dwork[iw17 + 1], &c__1, &dwork[iw15 + 1], &c__1);

/*           Compute y(1)*diag(B) - A. */

#line 1619 "AB13MD.f"
	i__1 = *n;
#line 1619 "AB13MD.f"
	for (j = 1; j <= i__1; ++j) {
#line 1620 "AB13MD.f"
	    i__2 = *n;
#line 1620 "AB13MD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 1621 "AB13MD.f"
		if (i__ == j) {
#line 1622 "AB13MD.f"
		    i__3 = iz15 + i__ + (i__ - 1) * *n;
#line 1622 "AB13MD.f"
		    d__1 = dwork[iw13 + 1] * dwork[iw24 + i__];
#line 1622 "AB13MD.f"
		    z__2.r = d__1, z__2.i = 0.;
#line 1622 "AB13MD.f"
		    i__4 = iz7 + i__ + (i__ - 1) * *n;
#line 1622 "AB13MD.f"
		    z__1.r = z__2.r - zwork[i__4].r, z__1.i = z__2.i - zwork[
			    i__4].i;
#line 1622 "AB13MD.f"
		    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 1624 "AB13MD.f"
		} else {
#line 1625 "AB13MD.f"
		    i__3 = iz15 + i__ + (j - 1) * *n;
#line 1625 "AB13MD.f"
		    i__4 = iz7 + i__ + (j - 1) * *n;
#line 1625 "AB13MD.f"
		    z__1.r = -zwork[i__4].r, z__1.i = -zwork[i__4].i;
#line 1625 "AB13MD.f"
		    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 1626 "AB13MD.f"
		}
#line 1627 "AB13MD.f"
/* L1430: */
#line 1627 "AB13MD.f"
	    }
#line 1628 "AB13MD.f"
/* L1440: */
#line 1628 "AB13MD.f"
	}

/*           Compute eig( y(1)*diag(B)-A ). */

#line 1632 "AB13MD.f"
	i__1 = *lzwork - izwrk;
#line 1632 "AB13MD.f"
	zgees_("N", "N", (L_fp)select_, n, &zwork[iz15 + 1], n, &sdim, &zwork[
		iz16 + 1], &zwork[1], n, &zwork[izwrk + 1], &i__1, &dwork[
		iwrk + 1], bwork, &info2, (ftnlen)1, (ftnlen)1);
#line 1635 "AB13MD.f"
	if (info2 > 0) {
#line 1636 "AB13MD.f"
	    *info = 6;
#line 1637 "AB13MD.f"
	    return 0;
#line 1638 "AB13MD.f"
	}
#line 1639 "AB13MD.f"
	i__1 = izwrk + 1;
#line 1639 "AB13MD.f"
	lza = (integer) zwork[i__1].r;
#line 1640 "AB13MD.f"
	lzamax = max(lza,lzamax);
#line 1641 "AB13MD.f"
	i__1 = iz16 + 1;
#line 1641 "AB13MD.f"
	emin = zwork[i__1].r;
#line 1642 "AB13MD.f"
	if (*n > 1) {
#line 1643 "AB13MD.f"
	    i__1 = *n;
#line 1643 "AB13MD.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 1644 "AB13MD.f"
		i__2 = iz16 + i__;
#line 1644 "AB13MD.f"
		if (zwork[i__2].r < emin) {
#line 1644 "AB13MD.f"
		    i__3 = iz16 + i__;
#line 1644 "AB13MD.f"
		    emin = zwork[i__3].r;
#line 1644 "AB13MD.f"
		}
#line 1646 "AB13MD.f"
/* L1450: */
#line 1646 "AB13MD.f"
	    }
#line 1647 "AB13MD.f"
	}
#line 1648 "AB13MD.f"
	pos = TRUE_;
#line 1649 "AB13MD.f"
	i__1 = *m - 1;
#line 1649 "AB13MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 1650 "AB13MD.f"
	    dwork[iw25 + i__] = dwork[iw22 + i__] - .01;
#line 1651 "AB13MD.f"
	    dwork[iw25 + *m - 1 + i__] = 100. - dwork[iw22 + i__];
#line 1652 "AB13MD.f"
/* L1460: */
#line 1652 "AB13MD.f"
	}
#line 1653 "AB13MD.f"
	if (mr > 0) {
#line 1654 "AB13MD.f"
	    i__1 = mr;
#line 1654 "AB13MD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1655 "AB13MD.f"
		dwork[iw25 + (*m - 1 << 1) + i__] = dwork[iw23 + i__] + tau;
#line 1656 "AB13MD.f"
		dwork[iw25 + (*m - 1 << 1) + mr + i__] = tau - dwork[iw23 + 
			i__];
#line 1657 "AB13MD.f"
/* L1470: */
#line 1657 "AB13MD.f"
	    }
#line 1658 "AB13MD.f"
	}
#line 1659 "AB13MD.f"
	temp = dwork[iw25 + 1];
#line 1660 "AB13MD.f"
	i__1 = mt << 1;
#line 1660 "AB13MD.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 1661 "AB13MD.f"
	    if (dwork[iw25 + i__] < temp) {
#line 1661 "AB13MD.f"
		temp = dwork[iw25 + i__];
#line 1661 "AB13MD.f"
	    }
#line 1662 "AB13MD.f"
/* L1480: */
#line 1662 "AB13MD.f"
	}
#line 1663 "AB13MD.f"
	if (temp <= 0. || emin <= 0.) {
#line 1663 "AB13MD.f"
	    pos = FALSE_;
#line 1663 "AB13MD.f"
	}
#line 1664 "AB13MD.f"
	goto L1390;
#line 1665 "AB13MD.f"
    }
#line 1666 "AB13MD.f"
L1490:

/*        Set y1 = ( y + yw ) / 2. */

#line 1670 "AB13MD.f"
    i__1 = mt + 1;
#line 1670 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1671 "AB13MD.f"
	dwork[iw16 + i__] = (dwork[iw13 + i__] + dwork[iw15 + i__]) / 2.;
#line 1673 "AB13MD.f"
/* L1500: */
#line 1673 "AB13MD.f"
    }

/*        Store xD. */

#line 1677 "AB13MD.f"
    i__1 = *m - 1;
#line 1677 "AB13MD.f"
    dcopy_(&i__1, &dwork[iw16 + 2], &c__1, &dwork[iw22 + 1], &c__1);
#line 1678 "AB13MD.f"
    if (mr > 0) {

/*           Store xG. */

#line 1682 "AB13MD.f"
	dcopy_(&mr, &dwork[iw16 + *m + 1], &c__1, &dwork[iw23 + 1], &c__1);
#line 1683 "AB13MD.f"
    }

/*        Compute A(:) = A0 + AA*y1(2:mt+1). */

#line 1687 "AB13MD.f"
    i__1 = mt;
#line 1687 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1688 "AB13MD.f"
	i__2 = iz10 + i__;
#line 1688 "AB13MD.f"
	i__3 = iw16 + 1 + i__;
#line 1688 "AB13MD.f"
	z__1.r = dwork[i__3], z__1.i = 0.;
#line 1688 "AB13MD.f"
	zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
#line 1689 "AB13MD.f"
/* L1510: */
#line 1689 "AB13MD.f"
    }
#line 1690 "AB13MD.f"
    i__1 = *n * *n;
#line 1690 "AB13MD.f"
    zcopy_(&i__1, &zwork[iz5 + 1], &c__1, &zwork[iz7 + 1], &c__1);
#line 1691 "AB13MD.f"
    i__1 = *n * *n;
#line 1691 "AB13MD.f"
    i__2 = *n * *n;
#line 1691 "AB13MD.f"
    zgemv_("N", &i__1, &mt, &c_b2, &zwork[iz6 + 1], &i__2, &zwork[iz10 + 1], &
	    c__1, &c_b2, &zwork[iz7 + 1], &c__1, (ftnlen)1);

/*        Compute diag( B ) = B0d + BBd*xD. */

#line 1696 "AB13MD.f"
    dcopy_(n, &dwork[iw7 + 1], &c__1, &dwork[iw24 + 1], &c__1);
#line 1697 "AB13MD.f"
    i__1 = *m - 1;
#line 1697 "AB13MD.f"
    dgemv_("N", n, &i__1, &c_b11, &dwork[iw8 + 1], n, &dwork[iw22 + 1], &c__1,
	     &c_b11, &dwork[iw24 + 1], &c__1, (ftnlen)1);

/*        Compute y1(1)*diag(B) - A. */

#line 1702 "AB13MD.f"
    i__1 = *n;
#line 1702 "AB13MD.f"
    for (j = 1; j <= i__1; ++j) {
#line 1703 "AB13MD.f"
	i__2 = *n;
#line 1703 "AB13MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 1704 "AB13MD.f"
	    if (i__ == j) {
#line 1705 "AB13MD.f"
		i__3 = iz15 + i__ + (i__ - 1) * *n;
#line 1705 "AB13MD.f"
		d__1 = dwork[iw16 + 1] * dwork[iw24 + i__];
#line 1705 "AB13MD.f"
		z__2.r = d__1, z__2.i = 0.;
#line 1705 "AB13MD.f"
		i__4 = iz7 + i__ + (i__ - 1) * *n;
#line 1705 "AB13MD.f"
		z__1.r = z__2.r - zwork[i__4].r, z__1.i = z__2.i - zwork[i__4]
			.i;
#line 1705 "AB13MD.f"
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 1707 "AB13MD.f"
	    } else {
#line 1708 "AB13MD.f"
		i__3 = iz15 + i__ + (j - 1) * *n;
#line 1708 "AB13MD.f"
		i__4 = iz7 + i__ + (j - 1) * *n;
#line 1708 "AB13MD.f"
		z__1.r = -zwork[i__4].r, z__1.i = -zwork[i__4].i;
#line 1708 "AB13MD.f"
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 1709 "AB13MD.f"
	    }
#line 1710 "AB13MD.f"
/* L1520: */
#line 1710 "AB13MD.f"
	}
#line 1711 "AB13MD.f"
/* L1530: */
#line 1711 "AB13MD.f"
    }

/*        Compute eig( y1(1)*diag(B)-A ). */

#line 1715 "AB13MD.f"
    i__1 = *lzwork - izwrk;
#line 1715 "AB13MD.f"
    zgees_("N", "N", (L_fp)select_, n, &zwork[iz15 + 1], n, &sdim, &zwork[
	    iz16 + 1], &zwork[1], n, &zwork[izwrk + 1], &i__1, &dwork[iwrk + 
	    1], bwork, &info2, (ftnlen)1, (ftnlen)1);
#line 1718 "AB13MD.f"
    if (info2 > 0) {
#line 1719 "AB13MD.f"
	*info = 6;
#line 1720 "AB13MD.f"
	return 0;
#line 1721 "AB13MD.f"
    }
#line 1722 "AB13MD.f"
    i__1 = izwrk + 1;
#line 1722 "AB13MD.f"
    lza = (integer) zwork[i__1].r;
#line 1723 "AB13MD.f"
    lzamax = max(lza,lzamax);
#line 1724 "AB13MD.f"
    i__1 = iz16 + 1;
#line 1724 "AB13MD.f"
    emin = zwork[i__1].r;
#line 1725 "AB13MD.f"
    if (*n > 1) {
#line 1726 "AB13MD.f"
	i__1 = *n;
#line 1726 "AB13MD.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 1727 "AB13MD.f"
	    i__2 = iz16 + i__;
#line 1727 "AB13MD.f"
	    if (zwork[i__2].r < emin) {
#line 1727 "AB13MD.f"
		i__3 = iz16 + i__;
#line 1727 "AB13MD.f"
		emin = zwork[i__3].r;
#line 1727 "AB13MD.f"
	    }
#line 1729 "AB13MD.f"
/* L1540: */
#line 1729 "AB13MD.f"
	}
#line 1730 "AB13MD.f"
    }
#line 1731 "AB13MD.f"
    pos = TRUE_;
#line 1732 "AB13MD.f"
    i__1 = *m - 1;
#line 1732 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1733 "AB13MD.f"
	dwork[iw25 + i__] = dwork[iw22 + i__] - .01;
#line 1734 "AB13MD.f"
	dwork[iw25 + *m - 1 + i__] = 100. - dwork[iw22 + i__];
#line 1735 "AB13MD.f"
/* L1550: */
#line 1735 "AB13MD.f"
    }
#line 1736 "AB13MD.f"
    if (mr > 0) {
#line 1737 "AB13MD.f"
	i__1 = mr;
#line 1737 "AB13MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 1738 "AB13MD.f"
	    dwork[iw25 + (*m - 1 << 1) + i__] = dwork[iw23 + i__] + tau;
#line 1739 "AB13MD.f"
	    dwork[iw25 + (*m - 1 << 1) + mr + i__] = tau - dwork[iw23 + i__];
#line 1740 "AB13MD.f"
/* L1560: */
#line 1740 "AB13MD.f"
	}
#line 1741 "AB13MD.f"
    }
#line 1742 "AB13MD.f"
    temp = dwork[iw25 + 1];
#line 1743 "AB13MD.f"
    i__1 = mt << 1;
#line 1743 "AB13MD.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 1744 "AB13MD.f"
	if (dwork[iw25 + i__] < temp) {
#line 1744 "AB13MD.f"
	    temp = dwork[iw25 + i__];
#line 1744 "AB13MD.f"
	}
#line 1745 "AB13MD.f"
/* L1570: */
#line 1745 "AB13MD.f"
    }
#line 1746 "AB13MD.f"
    if (temp <= 0. || emin <= 0.) {
#line 1746 "AB13MD.f"
	pos = FALSE_;
#line 1746 "AB13MD.f"
    }
#line 1747 "AB13MD.f"
    if (pos) {

/*           Set yw = y1. */

#line 1751 "AB13MD.f"
	i__1 = mt + 1;
#line 1751 "AB13MD.f"
	dcopy_(&i__1, &dwork[iw16 + 1], &c__1, &dwork[iw15 + 1], &c__1);
#line 1752 "AB13MD.f"
    } else {

/*           Set y = y1. */

#line 1756 "AB13MD.f"
	i__1 = mt + 1;
#line 1756 "AB13MD.f"
	dcopy_(&i__1, &dwork[iw16 + 1], &c__1, &dwork[iw13 + 1], &c__1);
#line 1757 "AB13MD.f"
    }
#line 1758 "AB13MD.f"
    i__1 = mt + 1;
#line 1758 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1759 "AB13MD.f"
	dwork[iw33 + i__] = dwork[iw13 + i__] - dwork[iw15 + i__];
#line 1760 "AB13MD.f"
/* L1580: */
#line 1760 "AB13MD.f"
    }
#line 1761 "AB13MD.f"
    i__1 = mt + 1;
#line 1761 "AB13MD.f"
    i__2 = mt + 1;
#line 1761 "AB13MD.f"
    ynorm1 = dlange_("F", &i__1, &c__1, &dwork[iw33 + 1], &i__2, &dwork[1], (
	    ftnlen)1);
#line 1762 "AB13MD.f"
    i__1 = mt + 1;
#line 1762 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1763 "AB13MD.f"
	dwork[iw33 + i__] = dwork[iw13 + i__] - dwork[iw14 + i__];
#line 1764 "AB13MD.f"
/* L1590: */
#line 1764 "AB13MD.f"
    }
#line 1765 "AB13MD.f"
    i__1 = mt + 1;
#line 1765 "AB13MD.f"
    i__2 = mt + 1;
#line 1765 "AB13MD.f"
    ynorm2 = dlange_("F", &i__1, &c__1, &dwork[iw33 + 1], &i__2, &dwork[1], (
	    ftnlen)1);
#line 1766 "AB13MD.f"
    if (ynorm1 < ynorm2 * .01) {
#line 1766 "AB13MD.f"
	goto L1600;
#line 1766 "AB13MD.f"
    }
#line 1767 "AB13MD.f"
    goto L1490;

/*        Compute c. */

#line 1771 "AB13MD.f"
L1600:
#line 1771 "AB13MD.f"
    i__1 = mt + 1;
#line 1771 "AB13MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1772 "AB13MD.f"
	dwork[iw33 + i__] = dwork[iw15 + i__] - dwork[iw14 + i__];
#line 1773 "AB13MD.f"
/* L1610: */
#line 1773 "AB13MD.f"
    }
#line 1774 "AB13MD.f"
    i__1 = mt + 1;
#line 1774 "AB13MD.f"
    i__2 = mt + 1;
#line 1774 "AB13MD.f"
    c__ = dlange_("F", &i__1, &c__1, &dwork[iw33 + 1], &i__2, &dwork[1], (
	    ftnlen)1);

/*        Set x = yw(2:mt+1). */

#line 1778 "AB13MD.f"
    dcopy_(&mt, &dwork[iw15 + 2], &c__1, &x[1], &c__1);
#line 1779 "AB13MD.f"
    goto L390;

/* *** Last line of AB13MD *** */
} /* ab13md_ */

