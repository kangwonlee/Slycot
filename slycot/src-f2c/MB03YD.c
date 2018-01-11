#line 1 "MB03YD.f"
/* MB03YD.f -- translated by f2c (version 20100827).
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

#line 1 "MB03YD.f"
/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;

/* Subroutine */ int mb03yd_(logical *wantt, logical *wantq, logical *wantz, 
	integer *n, integer *ilo, integer *ihi, integer *iloq, integer *ihiq, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	q, integer *ldq, doublereal *z__, integer *ldz, doublereal *alphar, 
	doublereal *alphai, doublereal *beta, doublereal *dwork, integer *
	ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, k, l;
    static doublereal v[3], w[3];
    static integer i1, i2, kk, nh, nq, nr;
    static doublereal cs1, cs2, cs3, sn1, sn2, sn3;
    static integer itn, its;
    static doublereal ulp, tst, temp, ovfl, unfl;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal tauv, tauw, gamma, alpha, delta;
    static integer iseed[4];
    extern /* Subroutine */ int mb03ya_(logical *, logical *, logical *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *);
    static doublereal betax;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), mb03yt_(doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *), dlabad_(
	    doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    extern doublereal dlanhs_(char *, integer *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), xerbla_(char *, 
	    integer *, ftnlen), dlarfx_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     ftnlen), dlarnv_(integer *, integer *, integer *, doublereal *);
    static doublereal smlnum;


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

/*     To deal with small subtasks of the product eigenvalue problem. */

/*     MB03YD is an auxiliary routine called by SLICOT Library routine */
/*     MB03XP. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     WANTT   LOGICAL */
/*             Indicates whether the user wishes to compute the full */
/*             Schur form or the eigenvalues only, as follows: */
/*             = .TRUE. :  Compute the full Schur form; */
/*             = .FALSE.:  compute the eigenvalues only. */

/*     WANTQ   LOGICAL */
/*             Indicates whether or not the user wishes to accumulate */
/*             the matrix Q as follows: */
/*             = .TRUE. :  The matrix Q is updated; */
/*             = .FALSE.:  the matrix Q is not required. */

/*     WANTZ   LOGICAL */
/*             Indicates whether or not the user wishes to accumulate */
/*             the matrix Z as follows: */
/*             = .TRUE. :  The matrix Z is updated; */
/*             = .FALSE.:  the matrix Z is not required. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A and B. N >= 0. */

/*     ILO     (input) INTEGER */
/*     IHI     (input) INTEGER */
/*             It is assumed that the matrices A and B are already */
/*             (quasi) upper triangular in rows and columns 1:ILO-1 and */
/*             IHI+1:N. The routine works primarily with the submatrices */
/*             in rows and columns ILO to IHI, but applies the */
/*             transformations to all the rows and columns of the */
/*             matrices A and B, if WANTT = .TRUE.. */
/*             1 <= ILO <= max(1,N); min(ILO,N) <= IHI <= N. */

/*     ILOQ    (input) INTEGER */
/*     IHIQ    (input) INTEGER */
/*             Specify the rows of Q and Z to which transformations */
/*             must be applied if WANTQ = .TRUE. and WANTZ = .TRUE., */
/*             respectively. */
/*             1 <= ILOQ <= ILO; IHI <= IHIQ <= N. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper Hessenberg matrix A. */
/*             On exit, if WANTT = .TRUE., the leading N-by-N part of */
/*             this array is upper quasi-triangular in rows and columns */
/*             ILO:IHI. */
/*             If WANTT = .FALSE., the diagonal elements and 2-by-2 */
/*             diagonal blocks of A will be correct, but the remaining */
/*             parts of A are unspecified on exit. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper triangular matrix B. */
/*             On exit, if WANTT = .TRUE., the leading N-by-N part of */
/*             this array contains the transformed upper triangular */
/*             matrix. 2-by-2 blocks in B corresponding to 2-by-2 blocks */
/*             in A will be reduced to positive diagonal form. (I.e., if */
/*             A(j+1,j) is non-zero, then B(j+1,j)=B(j,j+1)=0 and B(j,j) */
/*             and B(j+1,j+1) will be positive.) */
/*             If WANTT = .FALSE., the elements corresponding to diagonal */
/*             elements and 2-by-2 diagonal blocks in A will be correct, */
/*             but the remaining parts of B are unspecified on exit. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             On entry, if WANTQ = .TRUE., then the leading N-by-N part */
/*             of this array must contain the current matrix Q of */
/*             transformations accumulated by MB03XP. */
/*             On exit, if WANTQ = .TRUE., then the leading N-by-N part */
/*             of this array contains the matrix Q updated in the */
/*             submatrix Q(ILOQ:IHIQ,ILO:IHI). */
/*             If WANTQ = .FALSE., Q is not referenced. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q.  LDQ >= 1. */
/*             If WANTQ = .TRUE., LDQ >= MAX(1,N). */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             On entry, if WANTZ = .TRUE., then the leading N-by-N part */
/*             of this array must contain the current matrix Z of */
/*             transformations accumulated by MB03XP. */
/*             On exit, if WANTZ = .TRUE., then the leading N-by-N part */
/*             of this array contains the matrix Z updated in the */
/*             submatrix Z(ILOQ:IHIQ,ILO:IHI). */
/*             If WANTZ = .FALSE., Z is not referenced. */

/*     LDZ     INTEGER */
/*             The leading dimension of the array Z.  LDZ >= 1. */
/*             If WANTZ = .TRUE., LDZ >= MAX(1,N). */

/*     ALPHAR  (output) DOUBLE PRECISION array, dimension (N) */
/*     ALPHAI  (output) DOUBLE PRECISION array, dimension (N) */
/*     BETA    (output) DOUBLE PRECISION array, dimension (N) */
/*             The i-th (ILO <= i <= IHI) computed eigenvalue is given */
/*             by BETA(I) * ( ALPHAR(I) + sqrt(-1)*ALPHAI(I) ). If two */
/*             eigenvalues are computed as a complex conjugate pair, */
/*             they are stored in consecutive elements of ALPHAR, ALPHAI */
/*             and BETA. If WANTT = .TRUE., the eigenvalues are stored in */
/*             the same order as on the diagonals of the Schur forms of */
/*             A and B. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = -19,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, then MB03YD failed to compute the Schur */
/*                   form in a total of 30*(IHI-ILO+1) iterations; */
/*                   elements i+1:n of ALPHAR, ALPHAI and BETA contain */
/*                   successfully computed eigenvalues. */

/*     METHOD */

/*     The implemented algorithm is a double-shift version of the */
/*     periodic QR algorithm described in [1,3] with some minor */
/*     modifications [2]. The eigenvalues are computed via an implicit */
/*     complex single shift algorithm. */

/*     REFERENCES */

/*     [1] Bojanczyk, A.W., Golub, G.H., and Van Dooren, P. */
/*         The periodic Schur decomposition: Algorithms and applications. */
/*         Proc. of the SPIE Conference (F.T. Luk, Ed.), 1770, pp. 31-42, */
/*         1992. */

/*     [2] Kressner, D. */
/*         An efficient and reliable implementation of the periodic QZ */
/*         algorithm. Proc. of the IFAC Workshop on Periodic Control */
/*         Systems, pp. 187-192, 2001. */

/*     [3] Van Loan, C. */
/*         Generalized Singular Values with Algorithms and Applications. */
/*         Ph. D. Thesis, University of Michigan, 1973. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires O(N**3) floating point operations and is */
/*     backward stable. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLAPQR). */

/*     KEYWORDS */

/*     Eigenvalue, eigenvalue decomposition, Hessenberg form, orthogonal */
/*     transformation, (periodic) Schur form */

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

/*     Check the scalar input parameters. */

#line 235 "MB03YD.f"
    /* Parameter adjustments */
#line 235 "MB03YD.f"
    a_dim1 = *lda;
#line 235 "MB03YD.f"
    a_offset = 1 + a_dim1;
#line 235 "MB03YD.f"
    a -= a_offset;
#line 235 "MB03YD.f"
    b_dim1 = *ldb;
#line 235 "MB03YD.f"
    b_offset = 1 + b_dim1;
#line 235 "MB03YD.f"
    b -= b_offset;
#line 235 "MB03YD.f"
    q_dim1 = *ldq;
#line 235 "MB03YD.f"
    q_offset = 1 + q_dim1;
#line 235 "MB03YD.f"
    q -= q_offset;
#line 235 "MB03YD.f"
    z_dim1 = *ldz;
#line 235 "MB03YD.f"
    z_offset = 1 + z_dim1;
#line 235 "MB03YD.f"
    z__ -= z_offset;
#line 235 "MB03YD.f"
    --alphar;
#line 235 "MB03YD.f"
    --alphai;
#line 235 "MB03YD.f"
    --beta;
#line 235 "MB03YD.f"
    --dwork;
#line 235 "MB03YD.f"

#line 235 "MB03YD.f"
    /* Function Body */
#line 235 "MB03YD.f"
    *info = 0;
#line 236 "MB03YD.f"
    nh = *ihi - *ilo + 1;
#line 237 "MB03YD.f"
    nq = *ihiq - *iloq + 1;
#line 238 "MB03YD.f"
    if (*n < 0) {
#line 239 "MB03YD.f"
	*info = -4;
#line 240 "MB03YD.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 241 "MB03YD.f"
	*info = -5;
#line 242 "MB03YD.f"
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
#line 243 "MB03YD.f"
	*info = -6;
#line 244 "MB03YD.f"
    } else if (*iloq < 1 || *iloq > *ilo) {
#line 245 "MB03YD.f"
	*info = -7;
#line 246 "MB03YD.f"
    } else if (*ihiq < *ihi || *ihiq > *n) {
#line 247 "MB03YD.f"
	*info = -8;
#line 248 "MB03YD.f"
    } else if (*lda < max(1,*n)) {
#line 249 "MB03YD.f"
	*info = -10;
#line 250 "MB03YD.f"
    } else if (*ldb < max(1,*n)) {
#line 251 "MB03YD.f"
	*info = -12;
#line 252 "MB03YD.f"
    } else if (*ldq < 1 || *wantq && *ldq < *n) {
#line 253 "MB03YD.f"
	*info = -14;
#line 254 "MB03YD.f"
    } else if (*ldz < 1 || *wantz && *ldz < *n) {
#line 255 "MB03YD.f"
	*info = -16;
#line 256 "MB03YD.f"
    } else if (*ldwork < max(1,*n)) {
#line 257 "MB03YD.f"
	dwork[1] = (doublereal) max(1,*n);
#line 258 "MB03YD.f"
	*info = -21;
#line 259 "MB03YD.f"
    }

/*     Return if there were illegal values. */

#line 263 "MB03YD.f"
    if (*info != 0) {
#line 264 "MB03YD.f"
	i__1 = -(*info);
#line 264 "MB03YD.f"
	xerbla_("MB03YD", &i__1, (ftnlen)6);
#line 265 "MB03YD.f"
	return 0;
#line 266 "MB03YD.f"
    }

/*     Quick return if possible. */

#line 270 "MB03YD.f"
    if (*n == 0) {
#line 270 "MB03YD.f"
	return 0;
#line 270 "MB03YD.f"
    }

/*     Set machine-dependent constants for the stopping criterion. */

#line 275 "MB03YD.f"
    unfl = dlamch_("Safe minimum", (ftnlen)12);
#line 276 "MB03YD.f"
    ovfl = 1. / unfl;
#line 277 "MB03YD.f"
    dlabad_(&unfl, &ovfl);
#line 278 "MB03YD.f"
    ulp = dlamch_("Precision", (ftnlen)9);
#line 279 "MB03YD.f"
    smlnum = unfl * (nh / ulp);

/*     I1 and I2 are the indices of the first rows and last columns of */
/*     A and B to which transformations must be applied. */

#line 284 "MB03YD.f"
    i1 = 1;
#line 285 "MB03YD.f"
    i2 = *n;
#line 286 "MB03YD.f"
    iseed[0] = 1;
#line 287 "MB03YD.f"
    iseed[1] = 0;
#line 288 "MB03YD.f"
    iseed[2] = 0;
#line 289 "MB03YD.f"
    iseed[3] = 1;

/*     ITN is the maximal number of QR iterations. */

#line 293 "MB03YD.f"
    itn = nh * 30;

/*     Main loop. Eigenvalues I+1:IHI have converged. Either L = ILO */
/*     or A(L,L-1) is negligible. */

#line 298 "MB03YD.f"
    i__ = *ihi;
#line 299 "MB03YD.f"
L10:
#line 300 "MB03YD.f"
    l = *ilo;
#line 301 "MB03YD.f"
    if (i__ < *ilo) {
#line 301 "MB03YD.f"
	goto L120;
#line 301 "MB03YD.f"
    }

/*     Perform periodic QR iteration on rows and columns ILO to I of A */
/*     and B until a submatrix of order 1 or 2 splits off at the bottom. */

#line 307 "MB03YD.f"
    i__1 = itn;
#line 307 "MB03YD.f"
    for (its = 0; its <= i__1; ++its) {

/*        Look for deflations in A. */

#line 311 "MB03YD.f"
	i__2 = l + 1;
#line 311 "MB03YD.f"
	for (k = i__; k >= i__2; --k) {
#line 312 "MB03YD.f"
	    tst = (d__1 = a[k - 1 + (k - 1) * a_dim1], abs(d__1)) + (d__2 = a[
		    k + k * a_dim1], abs(d__2));
#line 313 "MB03YD.f"
	    if (tst == 0.) {
#line 313 "MB03YD.f"
		i__3 = i__ - l + 1;
#line 313 "MB03YD.f"
		tst = dlanhs_("1", &i__3, &a[l + l * a_dim1], lda, &dwork[1], 
			(ftnlen)1);
#line 313 "MB03YD.f"
	    }
/* Computing MAX */
#line 315 "MB03YD.f"
	    d__2 = ulp * tst;
#line 315 "MB03YD.f"
	    if ((d__1 = a[k + (k - 1) * a_dim1], abs(d__1)) <= max(d__2,
		    smlnum)) {
#line 315 "MB03YD.f"
		goto L30;
#line 315 "MB03YD.f"
	    }
#line 317 "MB03YD.f"
/* L20: */
#line 317 "MB03YD.f"
	}
#line 318 "MB03YD.f"
L30:

/*        Look for deflation in B if problem size is greater than 1. */

#line 322 "MB03YD.f"
	if (i__ - k >= 1) {
#line 323 "MB03YD.f"
	    i__2 = k;
#line 323 "MB03YD.f"
	    for (kk = i__; kk >= i__2; --kk) {
#line 324 "MB03YD.f"
		if (kk == i__) {
#line 325 "MB03YD.f"
		    tst = (d__1 = b[kk - 1 + kk * b_dim1], abs(d__1));
#line 326 "MB03YD.f"
		} else if (kk == k) {
#line 327 "MB03YD.f"
		    tst = (d__1 = b[kk + (kk + 1) * b_dim1], abs(d__1));
#line 328 "MB03YD.f"
		} else {
#line 329 "MB03YD.f"
		    tst = (d__1 = b[kk - 1 + kk * b_dim1], abs(d__1)) + (d__2 
			    = b[kk + (kk + 1) * b_dim1], abs(d__2));
#line 330 "MB03YD.f"
		}
#line 331 "MB03YD.f"
		if (tst == 0.) {
#line 331 "MB03YD.f"
		    i__3 = i__ - k + 1;
#line 331 "MB03YD.f"
		    tst = dlanhs_("1", &i__3, &b[k + k * b_dim1], ldb, &dwork[
			    1], (ftnlen)1);
#line 331 "MB03YD.f"
		}
/* Computing MAX */
#line 333 "MB03YD.f"
		d__2 = ulp * tst;
#line 333 "MB03YD.f"
		if ((d__1 = b[kk + kk * b_dim1], abs(d__1)) <= max(d__2,
			smlnum)) {
#line 333 "MB03YD.f"
		    goto L50;
#line 333 "MB03YD.f"
		}
#line 335 "MB03YD.f"
/* L40: */
#line 335 "MB03YD.f"
	    }
#line 336 "MB03YD.f"
	} else {
#line 337 "MB03YD.f"
	    kk = k - 1;
#line 338 "MB03YD.f"
	}
#line 339 "MB03YD.f"
L50:
#line 340 "MB03YD.f"
	if (kk >= k) {

/*           B has an element close to zero at position (KK,KK). */

#line 344 "MB03YD.f"
	    b[kk + kk * b_dim1] = 0.;
#line 345 "MB03YD.f"
	    mb03ya_(wantt, wantq, wantz, n, &k, &i__, iloq, ihiq, &kk, &a[
		    a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &
		    z__[z_offset], ldz, info);
#line 347 "MB03YD.f"
	    k = kk + 1;
#line 348 "MB03YD.f"
	}
#line 349 "MB03YD.f"
	l = k;
#line 350 "MB03YD.f"
	if (l > *ilo) {

/*           A(L,L-1) is negligible. */

#line 354 "MB03YD.f"
	    a[l + (l - 1) * a_dim1] = 0.;
#line 355 "MB03YD.f"
	}

/*        Exit from loop if a submatrix of order 1 or 2 has split off. */

#line 359 "MB03YD.f"
	if (l >= i__ - 1) {
#line 359 "MB03YD.f"
	    goto L80;
#line 359 "MB03YD.f"
	}

/*        The active submatrices are now in rows and columns L:I. */

#line 364 "MB03YD.f"
	if (! (*wantt)) {
#line 365 "MB03YD.f"
	    i1 = l;
#line 366 "MB03YD.f"
	    i2 = i__;
#line 367 "MB03YD.f"
	}
#line 368 "MB03YD.f"
	if (its == 10 || its == 20) {

/*           Exceptional shift. The first column of the shift polynomial */
/*           is a pseudo-random vector. */

#line 373 "MB03YD.f"
	    dlarnv_(&c__3, iseed, &c__3, v);
#line 374 "MB03YD.f"
	} else {

/*           The implicit double shift is constructed via a partial */
/*           product QR factorization [2]. */

#line 379 "MB03YD.f"
	    dlartg_(&b[l + l * b_dim1], &b[i__ + i__ * b_dim1], &cs2, &sn2, &
		    temp);
#line 380 "MB03YD.f"
	    dlartg_(&temp, &b[i__ - 1 + i__ * b_dim1], &cs1, &sn1, &alpha);

#line 382 "MB03YD.f"
	    alpha = a[l + l * a_dim1] * cs2 - a[i__ + i__ * a_dim1] * sn2;
#line 383 "MB03YD.f"
	    betax = cs1 * (cs2 * a[l + 1 + l * a_dim1]);
#line 384 "MB03YD.f"
	    gamma = cs1 * (sn2 * a[i__ - 1 + i__ * a_dim1]) + sn1 * a[i__ - 1 
		    + (i__ - 1) * a_dim1];
#line 385 "MB03YD.f"
	    alpha = alpha * cs1 - a[i__ + (i__ - 1) * a_dim1] * sn1;
#line 386 "MB03YD.f"
	    dlartg_(&alpha, &betax, &cs1, &sn1, &temp);

#line 388 "MB03YD.f"
	    dlartg_(&temp, &gamma, &cs2, &sn2, &alpha);
#line 389 "MB03YD.f"
	    alpha = cs2;
#line 390 "MB03YD.f"
	    gamma = a[i__ - 1 + (i__ - 1) * a_dim1] * cs1 * cs2 + a[i__ + (
		    i__ - 1) * a_dim1] * sn2;
#line 391 "MB03YD.f"
	    delta = a[i__ - 1 + (i__ - 1) * a_dim1] * sn1 * cs2;
#line 392 "MB03YD.f"
	    dlartg_(&gamma, &delta, &cs3, &sn3, &temp);
#line 393 "MB03YD.f"
	    dlartg_(&alpha, &temp, &cs2, &sn2, &alpha);

#line 395 "MB03YD.f"
	    alpha = (b[l + l * b_dim1] * cs1 + b[l + (l + 1) * b_dim1] * sn1) 
		    * cs2;
#line 396 "MB03YD.f"
	    betax = b[l + 1 + (l + 1) * b_dim1] * sn1 * cs2;
#line 397 "MB03YD.f"
	    gamma = b[i__ - 1 + (i__ - 1) * b_dim1] * sn2;
#line 398 "MB03YD.f"
	    dlartg_(&alpha, &betax, &cs1, &sn1, &temp);
#line 399 "MB03YD.f"
	    dlartg_(&temp, &gamma, &cs2, &sn2, &alpha);

#line 401 "MB03YD.f"
	    alpha = cs1 * a[l + l * a_dim1] + sn1 * a[l + (l + 1) * a_dim1];
#line 402 "MB03YD.f"
	    betax = cs1 * a[l + 1 + l * a_dim1] + sn1 * a[l + 1 + (l + 1) * 
		    a_dim1];
#line 403 "MB03YD.f"
	    gamma = sn1 * a[l + 2 + (l + 1) * a_dim1];

#line 405 "MB03YD.f"
	    v[0] = cs2 * alpha - sn2 * cs3;
#line 406 "MB03YD.f"
	    v[1] = cs2 * betax - sn2 * sn3;
#line 407 "MB03YD.f"
	    v[2] = gamma * cs2;
#line 408 "MB03YD.f"
	}

/*        Double-shift QR step */

#line 412 "MB03YD.f"
	i__2 = i__ - 1;
#line 412 "MB03YD.f"
	for (k = l; k <= i__2; ++k) {

/* Computing MIN */
#line 414 "MB03YD.f"
	    i__3 = 3, i__4 = i__ - k + 1;
#line 414 "MB03YD.f"
	    nr = min(i__3,i__4);
#line 415 "MB03YD.f"
	    if (k > l) {
#line 415 "MB03YD.f"
		dcopy_(&nr, &a[k + (k - 1) * a_dim1], &c__1, v, &c__1);
#line 415 "MB03YD.f"
	    }
#line 417 "MB03YD.f"
	    dlarfg_(&nr, v, &v[1], &c__1, &tauv);
#line 418 "MB03YD.f"
	    if (k > l) {
#line 419 "MB03YD.f"
		a[k + (k - 1) * a_dim1] = v[0];
#line 420 "MB03YD.f"
		a[k + 1 + (k - 1) * a_dim1] = 0.;
#line 421 "MB03YD.f"
		if (k < i__ - 1) {
#line 421 "MB03YD.f"
		    a[k + 2 + (k - 1) * a_dim1] = 0.;
#line 421 "MB03YD.f"
		}
#line 423 "MB03YD.f"
	    }

/*           Apply reflector V from the right to B in rows I1:min(K+2,I). */

#line 427 "MB03YD.f"
	    v[0] = 1.;
/* Computing MIN */
#line 428 "MB03YD.f"
	    i__4 = k + 2;
#line 428 "MB03YD.f"
	    i__3 = min(i__4,i__) - i1 + 1;
#line 428 "MB03YD.f"
	    dlarfx_("Right", &i__3, &nr, v, &tauv, &b[i1 + k * b_dim1], ldb, &
		    dwork[1], (ftnlen)5);

/*           Annihilate the introduced nonzeros in the K-th column. */

#line 433 "MB03YD.f"
	    dcopy_(&nr, &b[k + k * b_dim1], &c__1, w, &c__1);
#line 434 "MB03YD.f"
	    dlarfg_(&nr, w, &w[1], &c__1, &tauw);
#line 435 "MB03YD.f"
	    b[k + k * b_dim1] = w[0];
#line 436 "MB03YD.f"
	    b[k + 1 + k * b_dim1] = 0.;
#line 437 "MB03YD.f"
	    if (k < i__ - 1) {
#line 437 "MB03YD.f"
		b[k + 2 + k * b_dim1] = 0.;
#line 437 "MB03YD.f"
	    }

/*           Apply reflector W from the left to transform the rows of the */
/*           matrix B in columns K+1:I2. */

#line 443 "MB03YD.f"
	    w[0] = 1.;
#line 444 "MB03YD.f"
	    i__3 = i2 - k;
#line 444 "MB03YD.f"
	    dlarfx_("Left", &nr, &i__3, w, &tauw, &b[k + (k + 1) * b_dim1], 
		    ldb, &dwork[1], (ftnlen)4);

/*           Apply reflector V from the left to transform the rows of the */
/*           matrix A in columns K:I2. */

#line 450 "MB03YD.f"
	    i__3 = i2 - k + 1;
#line 450 "MB03YD.f"
	    dlarfx_("Left", &nr, &i__3, v, &tauv, &a[k + k * a_dim1], lda, &
		    dwork[1], (ftnlen)4);

/*           Apply reflector W from the right to transform the columns of */
/*           the matrix A in rows I1:min(K+3,I). */

/* Computing MIN */
#line 456 "MB03YD.f"
	    i__4 = k + 3;
#line 456 "MB03YD.f"
	    i__3 = min(i__4,i__) - i1 + 1;
#line 456 "MB03YD.f"
	    dlarfx_("Right", &i__3, &nr, w, &tauw, &a[i1 + k * a_dim1], lda, &
		    dwork[1], (ftnlen)5);

/*           Accumulate transformations in the matrices Q and Z. */

#line 461 "MB03YD.f"
	    if (*wantq) {
#line 461 "MB03YD.f"
		dlarfx_("Right", &nq, &nr, v, &tauv, &q[*iloq + k * q_dim1], 
			ldq, &dwork[1], (ftnlen)5);
#line 461 "MB03YD.f"
	    }
#line 464 "MB03YD.f"
	    if (*wantz) {
#line 464 "MB03YD.f"
		dlarfx_("Right", &nq, &nr, w, &tauw, &z__[*iloq + k * z_dim1],
			 ldz, &dwork[1], (ftnlen)5);
#line 464 "MB03YD.f"
	    }
#line 467 "MB03YD.f"
/* L60: */
#line 467 "MB03YD.f"
	}
#line 468 "MB03YD.f"
/* L70: */
#line 468 "MB03YD.f"
    }

/*     Failure to converge. */

#line 472 "MB03YD.f"
    *info = i__;
#line 473 "MB03YD.f"
    return 0;

#line 475 "MB03YD.f"
L80:

/*     Compute 1-by-1 or 2-by-2 subproblem. */

#line 479 "MB03YD.f"
    if (l == i__) {

/*        Standardize B, set ALPHAR, ALPHAI and BETA. */

#line 483 "MB03YD.f"
	if (b[i__ + i__ * b_dim1] < 0.) {
#line 484 "MB03YD.f"
	    if (*wantt) {
#line 485 "MB03YD.f"
		i__1 = i__;
#line 485 "MB03YD.f"
		for (k = i1; k <= i__1; ++k) {
#line 486 "MB03YD.f"
		    b[k + i__ * b_dim1] = -b[k + i__ * b_dim1];
#line 487 "MB03YD.f"
/* L90: */
#line 487 "MB03YD.f"
		}
#line 488 "MB03YD.f"
		i__1 = i2;
#line 488 "MB03YD.f"
		for (k = i__; k <= i__1; ++k) {
#line 489 "MB03YD.f"
		    a[i__ + k * a_dim1] = -a[i__ + k * a_dim1];
#line 490 "MB03YD.f"
/* L100: */
#line 490 "MB03YD.f"
		}
#line 491 "MB03YD.f"
	    } else {
#line 492 "MB03YD.f"
		b[i__ + i__ * b_dim1] = -b[i__ + i__ * b_dim1];
#line 493 "MB03YD.f"
		a[i__ + i__ * a_dim1] = -a[i__ + i__ * a_dim1];
#line 494 "MB03YD.f"
	    }
#line 495 "MB03YD.f"
	    if (*wantq) {
#line 496 "MB03YD.f"
		i__1 = *ihiq;
#line 496 "MB03YD.f"
		for (k = *iloq; k <= i__1; ++k) {
#line 497 "MB03YD.f"
		    q[k + i__ * q_dim1] = -q[k + i__ * q_dim1];
#line 498 "MB03YD.f"
/* L110: */
#line 498 "MB03YD.f"
		}
#line 499 "MB03YD.f"
	    }
#line 500 "MB03YD.f"
	}
#line 501 "MB03YD.f"
	alphar[i__] = a[i__ + i__ * a_dim1];
#line 502 "MB03YD.f"
	alphai[i__] = 0.;
#line 503 "MB03YD.f"
	beta[i__] = b[i__ + i__ * b_dim1];
#line 504 "MB03YD.f"
    } else if (l == i__ - 1) {

/*        A double block has converged. */
/*        Compute eigenvalues and standardize double block. */

#line 509 "MB03YD.f"
	mb03yt_(&a[i__ - 1 + (i__ - 1) * a_dim1], lda, &b[i__ - 1 + (i__ - 1) 
		* b_dim1], ldb, &alphar[i__ - 1], &alphai[i__ - 1], &beta[i__ 
		- 1], &cs1, &sn1, &cs2, &sn2);

/*        Apply transformation to rest of A and B. */

#line 514 "MB03YD.f"
	if (i2 > i__) {
#line 514 "MB03YD.f"
	    i__1 = i2 - i__;
#line 514 "MB03YD.f"
	    drot_(&i__1, &a[i__ - 1 + (i__ + 1) * a_dim1], lda, &a[i__ + (i__ 
		    + 1) * a_dim1], lda, &cs1, &sn1);
#line 514 "MB03YD.f"
	}
#line 516 "MB03YD.f"
	i__1 = i__ - i1 - 1;
#line 516 "MB03YD.f"
	drot_(&i__1, &a[i1 + (i__ - 1) * a_dim1], &c__1, &a[i1 + i__ * a_dim1]
		, &c__1, &cs2, &sn2);
#line 517 "MB03YD.f"
	if (i2 > i__) {
#line 517 "MB03YD.f"
	    i__1 = i2 - i__;
#line 517 "MB03YD.f"
	    drot_(&i__1, &b[i__ - 1 + (i__ + 1) * b_dim1], ldb, &b[i__ + (i__ 
		    + 1) * b_dim1], ldb, &cs2, &sn2);
#line 517 "MB03YD.f"
	}
#line 519 "MB03YD.f"
	i__1 = i__ - i1 - 1;
#line 519 "MB03YD.f"
	drot_(&i__1, &b[i1 + (i__ - 1) * b_dim1], &c__1, &b[i1 + i__ * b_dim1]
		, &c__1, &cs1, &sn1);

/*        Apply transformation to rest of Q and Z if desired. */

#line 523 "MB03YD.f"
	if (*wantq) {
#line 523 "MB03YD.f"
	    drot_(&nq, &q[*iloq + (i__ - 1) * q_dim1], &c__1, &q[*iloq + i__ *
		     q_dim1], &c__1, &cs1, &sn1);
#line 523 "MB03YD.f"
	}
#line 525 "MB03YD.f"
	if (*wantz) {
#line 525 "MB03YD.f"
	    drot_(&nq, &z__[*iloq + (i__ - 1) * z_dim1], &c__1, &z__[*iloq + 
		    i__ * z_dim1], &c__1, &cs2, &sn2);
#line 525 "MB03YD.f"
	}
#line 527 "MB03YD.f"
    }

/*     Decrement number of remaining iterations, and return to start of */
/*     the main loop with new value of I. */

#line 532 "MB03YD.f"
    itn -= its;
#line 533 "MB03YD.f"
    i__ = l - 1;
#line 534 "MB03YD.f"
    goto L10;

#line 536 "MB03YD.f"
L120:
#line 537 "MB03YD.f"
    dwork[1] = (doublereal) max(1,*n);
#line 538 "MB03YD.f"
    return 0;
/* *** Last line of MB03YD *** */
} /* mb03yd_ */

