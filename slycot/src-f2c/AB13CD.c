#line 1 "AB13CD.f"
/* AB13CD.f -- translated by f2c (version 20100827).
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

#line 1 "AB13CD.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublereal c_b16 = -1.;
static doublereal c_b17 = 1.;
static doublereal c_b38 = 0.;

doublereal ab13cd_(integer *n, integer *m, integer *np, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, doublereal *c__, integer *
	ldc, doublereal *d__, integer *ldd, doublereal *tol, integer *iwork, 
	doublereal *dwork, integer *ldwork, doublecomplex *cwork, integer *
	lcwork, logical *bwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal ret_val, d__1, d__2;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__, j, k, l, iw2, iw3, iw4, iw5, iw6, iw7, iw8, iw9;
    static doublereal den;
    static integer iw10, iw11, iw12;
    static doublereal rat;
    static integer icw2, icw3, icw4, sdim, iter;
    static doublereal temp;
    static integer iwrk, info2;
    extern /* Subroutine */ int ma02ed_(char *, integer *, doublereal *, 
	    integer *, ftnlen);
    static doublereal gamma, fpeak, omega;
    extern /* Subroutine */ int dgees_(char *, char *, L_fp, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen), dgemm_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    extern logical sb02cx_();
    extern /* Subroutine */ int dgesv_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *);
    extern logical sb02mv_();
    extern /* Subroutine */ int mb01rx_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static integer icwrk;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static doublereal wimax, wrmin;
    extern /* Subroutine */ int dposv_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen), dsyrk_(char *, char *, integer *, integer *, doublereal *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen), zgesv_(integer *, integer *, doublecomplex *, 
	    integer *, integer *, doublecomplex *, integer *, integer *);
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static doublereal gammal, gammau;
    extern /* Subroutine */ int dgesvd_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dlacpy_(char *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    static integer lwamax, lcwamx;
    static doublereal ratmax;
    extern /* Subroutine */ int dpotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer mincwr;
    static logical complx;
    extern /* Subroutine */ int zgesvd_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, integer *, ftnlen, ftnlen);
    static integer minwrk;
    extern /* Subroutine */ int dpotrs_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);


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

/*     To compute the H-infinity norm of the continuous-time stable */
/*     system */

/*                          | A | B | */
/*                   G(s) = |---|---| . */
/*                          | C | D | */

/*     FUNCTION VALUE */

/*     AB13CD  DOUBLE PRECISION */
/*             If INFO = 0, the H-infinity norm of the system, HNORM, */
/*             i.e., the peak gain of the frequency response (as measured */
/*             by the largest singular value in the MIMO case). */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The column size of the matrix B.  M >= 0. */

/*     NP      (input) INTEGER */
/*             The row size of the matrix C.  NP >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             system state matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             system input matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading NP-by-N part of this array must contain the */
/*             system output matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= max(1,NP). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading NP-by-M part of this array must contain the */
/*             system input/output matrix D. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= max(1,NP). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             Tolerance used to set the accuracy in determining the */
/*             norm. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension N */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal value */
/*             of LDWORK, and DWORK(2) contains the frequency where the */
/*             gain of the frequency response achieves its peak value */
/*             HNORM. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= max(2,4*N*N+2*M*M+3*M*N+M*NP+2*(N+NP)*NP+10*N+ */
/*                             6*max(M,NP)). */
/*             For good performance, LDWORK must generally be larger. */

/*     CWORK   COMPLEX*16 array, dimension (LCWORK) */
/*             On exit, if INFO = 0, CWORK(1) contains the optimal value */
/*             of LCWORK. */

/*     LCWORK  INTEGER */
/*             The dimension of the array CWORK. */
/*             LCWORK >= max(1,(N+M)*(N+NP)+3*max(M,NP)). */
/*             For good performance, LCWORK must generally be larger. */

/*     BWORK   LOGICAL array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the system is unstable; */
/*             = 2:  the tolerance is too small (the algorithm for */
/*                   computing the H-infinity norm did not converge); */
/*             = 3:  errors in computing the eigenvalues of A or of the */
/*                   Hamiltonian matrix (the QR algorithm did not */
/*                   converge); */
/*             = 4:  errors in computing singular values. */

/*     METHOD */

/*     The routine implements the method presented in [1]. */

/*     REFERENCES */

/*     [1] Bruinsma, N.A. and Steinbuch, M. */
/*         A fast algorithm to compute the Hinfinity-norm of a transfer */
/*         function matrix. */
/*         Systems & Control Letters, vol. 14, pp. 287-293, 1990. */

/*     NUMERICAL ASPECTS */

/*     If the algorithm does not converge (INFO = 2), the tolerance must */
/*     be increased. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, May 1999. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 1999, */
/*     Oct. 2000. */
/*     P.Hr. Petkov, October 2000. */
/*     A. Varga, October 2000. */
/*     Oct. 2001, V. Sima, Research Institute for Informatics, Bucharest. */

/*     KEYWORDS */

/*     H-infinity optimal control, robust control, system norm. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */

/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input scalar parameters. */

#line 211 "AB13CD.f"
    /* Parameter adjustments */
#line 211 "AB13CD.f"
    a_dim1 = *lda;
#line 211 "AB13CD.f"
    a_offset = 1 + a_dim1;
#line 211 "AB13CD.f"
    a -= a_offset;
#line 211 "AB13CD.f"
    b_dim1 = *ldb;
#line 211 "AB13CD.f"
    b_offset = 1 + b_dim1;
#line 211 "AB13CD.f"
    b -= b_offset;
#line 211 "AB13CD.f"
    c_dim1 = *ldc;
#line 211 "AB13CD.f"
    c_offset = 1 + c_dim1;
#line 211 "AB13CD.f"
    c__ -= c_offset;
#line 211 "AB13CD.f"
    d_dim1 = *ldd;
#line 211 "AB13CD.f"
    d_offset = 1 + d_dim1;
#line 211 "AB13CD.f"
    d__ -= d_offset;
#line 211 "AB13CD.f"
    --iwork;
#line 211 "AB13CD.f"
    --dwork;
#line 211 "AB13CD.f"
    --cwork;
#line 211 "AB13CD.f"
    --bwork;
#line 211 "AB13CD.f"

#line 211 "AB13CD.f"
    /* Function Body */
#line 211 "AB13CD.f"
    *info = 0;
#line 212 "AB13CD.f"
    if (*n < 0) {
#line 213 "AB13CD.f"
	*info = -1;
#line 214 "AB13CD.f"
    } else if (*m < 0) {
#line 215 "AB13CD.f"
	*info = -2;
#line 216 "AB13CD.f"
    } else if (*np < 0) {
#line 217 "AB13CD.f"
	*info = -3;
#line 218 "AB13CD.f"
    } else if (*lda < max(1,*n)) {
#line 219 "AB13CD.f"
	*info = -5;
#line 220 "AB13CD.f"
    } else if (*ldb < max(1,*n)) {
#line 221 "AB13CD.f"
	*info = -7;
#line 222 "AB13CD.f"
    } else if (*ldc < max(1,*np)) {
#line 223 "AB13CD.f"
	*info = -9;
#line 224 "AB13CD.f"
    } else if (*ldd < max(1,*np)) {
#line 225 "AB13CD.f"
	*info = -11;
#line 226 "AB13CD.f"
    }

/*     Compute workspace. */

/* Computing MAX */
#line 230 "AB13CD.f"
    i__1 = 2, i__2 = (*n << 2) * *n + (*m << 1) * *m + *m * 3 * *n + *m * *np 
	    + (*n + *np << 1) * *np + *n * 10 + max(*m,*np) * 6;
#line 230 "AB13CD.f"
    minwrk = max(i__1,i__2);
#line 232 "AB13CD.f"
    if (*ldwork < minwrk) {
#line 233 "AB13CD.f"
	*info = -15;
#line 234 "AB13CD.f"
    }
/* Computing MAX */
#line 235 "AB13CD.f"
    i__1 = 1, i__2 = (*n + *m) * (*n + *np) + max(*m,*np) * 3;
#line 235 "AB13CD.f"
    mincwr = max(i__1,i__2);
#line 236 "AB13CD.f"
    if (*lcwork < mincwr) {
#line 237 "AB13CD.f"
	*info = -17;
#line 238 "AB13CD.f"
    }
#line 239 "AB13CD.f"
    if (*info != 0) {
#line 240 "AB13CD.f"
	i__1 = -(*info);
#line 240 "AB13CD.f"
	xerbla_("AB13CD", &i__1, (ftnlen)6);
#line 241 "AB13CD.f"
	return ret_val;
#line 242 "AB13CD.f"
    }

/*     Quick return if possible. */

#line 246 "AB13CD.f"
    if (*m == 0 || *np == 0) {
#line 246 "AB13CD.f"
	return ret_val;
#line 246 "AB13CD.f"
    }

/*     Workspace usage. */

#line 250 "AB13CD.f"
    iw2 = *n;
#line 251 "AB13CD.f"
    iw3 = iw2 + *n;
#line 252 "AB13CD.f"
    iw4 = iw3 + *n * *n;
#line 253 "AB13CD.f"
    iw5 = iw4 + *n * *m;
#line 254 "AB13CD.f"
    iw6 = iw5 + *np * *m;
#line 255 "AB13CD.f"
    iwrk = iw6 + min(*np,*m);

/*     Determine the maximum singular value of G(infinity) = D . */

#line 259 "AB13CD.f"
    dlacpy_("Full", np, m, &d__[d_offset], ldd, &dwork[iw5 + 1], np, (ftnlen)
	    4);
#line 260 "AB13CD.f"
    i__1 = *ldwork - iwrk;
#line 260 "AB13CD.f"
    dgesvd_("N", "N", np, m, &dwork[iw5 + 1], np, &dwork[iw6 + 1], &dwork[1], 
	    np, &dwork[1], m, &dwork[iwrk + 1], &i__1, &info2, (ftnlen)1, (
	    ftnlen)1);
#line 263 "AB13CD.f"
    if (info2 > 0) {
#line 264 "AB13CD.f"
	*info = 4;
#line 265 "AB13CD.f"
	return ret_val;
#line 266 "AB13CD.f"
    }
#line 267 "AB13CD.f"
    gammal = dwork[iw6 + 1];
#line 268 "AB13CD.f"
    fpeak = 1e30;
#line 269 "AB13CD.f"
    lwamax = (integer) dwork[iwrk + 1] + iwrk;

/*     Quick return if N = 0 . */

#line 273 "AB13CD.f"
    if (*n == 0) {
#line 274 "AB13CD.f"
	ret_val = gammal;
#line 275 "AB13CD.f"
	dwork[1] = 2.;
#line 276 "AB13CD.f"
	dwork[2] = 0.;
#line 277 "AB13CD.f"
	cwork[1].r = 1., cwork[1].i = 0.;
#line 278 "AB13CD.f"
	return ret_val;
#line 279 "AB13CD.f"
    }

/*     Stability check. */

#line 283 "AB13CD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[iw3 + 1], n, (ftnlen)4);
#line 284 "AB13CD.f"
    i__1 = *ldwork - iwrk;
#line 284 "AB13CD.f"
    dgees_("N", "S", (L_fp)sb02mv_, n, &dwork[iw3 + 1], n, &sdim, &dwork[1], &
	    dwork[iw2 + 1], &dwork[1], n, &dwork[iwrk + 1], &i__1, &bwork[1], 
	    &info2, (ftnlen)1, (ftnlen)1);
#line 287 "AB13CD.f"
    if (info2 > 0) {
#line 288 "AB13CD.f"
	*info = 3;
#line 289 "AB13CD.f"
	return ret_val;
#line 290 "AB13CD.f"
    }
#line 291 "AB13CD.f"
    if (sdim < *n) {
#line 292 "AB13CD.f"
	*info = 1;
#line 293 "AB13CD.f"
	return ret_val;
#line 294 "AB13CD.f"
    }
/* Computing MAX */
#line 295 "AB13CD.f"
    i__1 = (integer) dwork[iwrk + 1] + iwrk;
#line 295 "AB13CD.f"
    lwamax = max(i__1,lwamax);

/*     Determine the maximum singular value of G(0) = -C*inv(A)*B + D . */

#line 299 "AB13CD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[iw3 + 1], n, (ftnlen)4);
#line 300 "AB13CD.f"
    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[iw4 + 1], n, (ftnlen)4);
#line 301 "AB13CD.f"
    dlacpy_("Full", np, m, &d__[d_offset], ldd, &dwork[iw5 + 1], np, (ftnlen)
	    4);
#line 302 "AB13CD.f"
    dgesv_(n, m, &dwork[iw3 + 1], n, &iwork[1], &dwork[iw4 + 1], n, &info2);
#line 304 "AB13CD.f"
    if (info2 > 0) {
#line 305 "AB13CD.f"
	*info = 1;
#line 306 "AB13CD.f"
	return ret_val;
#line 307 "AB13CD.f"
    }
#line 308 "AB13CD.f"
    dgemm_("N", "N", np, m, n, &c_b16, &c__[c_offset], ldc, &dwork[iw4 + 1], 
	    n, &c_b17, &dwork[iw5 + 1], np, (ftnlen)1, (ftnlen)1);
#line 310 "AB13CD.f"
    i__1 = *ldwork - iwrk;
#line 310 "AB13CD.f"
    dgesvd_("N", "N", np, m, &dwork[iw5 + 1], np, &dwork[iw6 + 1], &dwork[1], 
	    np, &dwork[1], m, &dwork[iwrk + 1], &i__1, &info2, (ftnlen)1, (
	    ftnlen)1);
#line 313 "AB13CD.f"
    if (info2 > 0) {
#line 314 "AB13CD.f"
	*info = 4;
#line 315 "AB13CD.f"
	return ret_val;
#line 316 "AB13CD.f"
    }
#line 317 "AB13CD.f"
    if (gammal < dwork[iw6 + 1]) {
#line 318 "AB13CD.f"
	gammal = dwork[iw6 + 1];
#line 319 "AB13CD.f"
	fpeak = 0.;
#line 320 "AB13CD.f"
    }
/* Computing MAX */
#line 321 "AB13CD.f"
    i__1 = (integer) dwork[iwrk + 1] + iwrk;
#line 321 "AB13CD.f"
    lwamax = max(i__1,lwamax);

/*     Find a frequency which is close to the peak frequency. */

#line 325 "AB13CD.f"
    complx = FALSE_;
#line 326 "AB13CD.f"
    i__1 = *n;
#line 326 "AB13CD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 327 "AB13CD.f"
	if (dwork[iw2 + i__] != 0.) {
#line 327 "AB13CD.f"
	    complx = TRUE_;
#line 327 "AB13CD.f"
	}
#line 328 "AB13CD.f"
/* L10: */
#line 328 "AB13CD.f"
    }
#line 329 "AB13CD.f"
    if (! complx) {
#line 330 "AB13CD.f"
	wrmin = abs(dwork[1]);
#line 331 "AB13CD.f"
	i__1 = *n;
#line 331 "AB13CD.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 332 "AB13CD.f"
	    if (wrmin > (d__1 = dwork[i__], abs(d__1))) {
#line 332 "AB13CD.f"
		wrmin = (d__2 = dwork[i__], abs(d__2));
#line 332 "AB13CD.f"
	    }
#line 333 "AB13CD.f"
/* L20: */
#line 333 "AB13CD.f"
	}
#line 334 "AB13CD.f"
	omega = wrmin;
#line 335 "AB13CD.f"
    } else {
#line 336 "AB13CD.f"
	ratmax = 0.;
#line 337 "AB13CD.f"
	i__1 = *n;
#line 337 "AB13CD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 338 "AB13CD.f"
	    den = dlapy2_(&dwork[i__], &dwork[iw2 + i__]);
#line 339 "AB13CD.f"
	    rat = (d__1 = dwork[iw2 + i__] / dwork[i__] / den, abs(d__1));
#line 340 "AB13CD.f"
	    if (ratmax < rat) {
#line 341 "AB13CD.f"
		ratmax = rat;
#line 342 "AB13CD.f"
		wimax = den;
#line 343 "AB13CD.f"
	    }
#line 344 "AB13CD.f"
/* L30: */
#line 344 "AB13CD.f"
	}
#line 345 "AB13CD.f"
	omega = wimax;
#line 346 "AB13CD.f"
    }

/*     Workspace usage. */

#line 350 "AB13CD.f"
    icw2 = *n * *n;
#line 351 "AB13CD.f"
    icw3 = icw2 + *n * *m;
#line 352 "AB13CD.f"
    icw4 = icw3 + *np * *n;
#line 353 "AB13CD.f"
    icwrk = icw4 + *np * *m;

/*     Determine the maximum singular value of */
/*     G(omega) = C*inv(j*omega*In - A)*B + D . */

#line 358 "AB13CD.f"
    i__1 = *n;
#line 358 "AB13CD.f"
    for (j = 1; j <= i__1; ++j) {
#line 359 "AB13CD.f"
	i__2 = *n;
#line 359 "AB13CD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 360 "AB13CD.f"
	    i__3 = i__ + (j - 1) * *n;
#line 360 "AB13CD.f"
	    d__1 = -a[i__ + j * a_dim1];
#line 360 "AB13CD.f"
	    cwork[i__3].r = d__1, cwork[i__3].i = 0.;
#line 361 "AB13CD.f"
/* L40: */
#line 361 "AB13CD.f"
	}
#line 362 "AB13CD.f"
	i__2 = j + (j - 1) * *n;
#line 362 "AB13CD.f"
	z__2.r = omega * 0., z__2.i = omega * 1.;
#line 362 "AB13CD.f"
	i__3 = j + j * a_dim1;
#line 362 "AB13CD.f"
	z__1.r = z__2.r - a[i__3], z__1.i = z__2.i;
#line 362 "AB13CD.f"
	cwork[i__2].r = z__1.r, cwork[i__2].i = z__1.i;
#line 363 "AB13CD.f"
/* L50: */
#line 363 "AB13CD.f"
    }
#line 364 "AB13CD.f"
    i__1 = *m;
#line 364 "AB13CD.f"
    for (j = 1; j <= i__1; ++j) {
#line 365 "AB13CD.f"
	i__2 = *n;
#line 365 "AB13CD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 366 "AB13CD.f"
	    i__3 = icw2 + i__ + (j - 1) * *n;
#line 366 "AB13CD.f"
	    i__4 = i__ + j * b_dim1;
#line 366 "AB13CD.f"
	    cwork[i__3].r = b[i__4], cwork[i__3].i = 0.;
#line 367 "AB13CD.f"
/* L60: */
#line 367 "AB13CD.f"
	}
#line 368 "AB13CD.f"
/* L70: */
#line 368 "AB13CD.f"
    }
#line 369 "AB13CD.f"
    i__1 = *n;
#line 369 "AB13CD.f"
    for (j = 1; j <= i__1; ++j) {
#line 370 "AB13CD.f"
	i__2 = *np;
#line 370 "AB13CD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 371 "AB13CD.f"
	    i__3 = icw3 + i__ + (j - 1) * *np;
#line 371 "AB13CD.f"
	    i__4 = i__ + j * c_dim1;
#line 371 "AB13CD.f"
	    cwork[i__3].r = c__[i__4], cwork[i__3].i = 0.;
#line 372 "AB13CD.f"
/* L80: */
#line 372 "AB13CD.f"
	}
#line 373 "AB13CD.f"
/* L90: */
#line 373 "AB13CD.f"
    }
#line 374 "AB13CD.f"
    i__1 = *m;
#line 374 "AB13CD.f"
    for (j = 1; j <= i__1; ++j) {
#line 375 "AB13CD.f"
	i__2 = *np;
#line 375 "AB13CD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 376 "AB13CD.f"
	    i__3 = icw4 + i__ + (j - 1) * *np;
#line 376 "AB13CD.f"
	    i__4 = i__ + j * d_dim1;
#line 376 "AB13CD.f"
	    cwork[i__3].r = d__[i__4], cwork[i__3].i = 0.;
#line 377 "AB13CD.f"
/* L100: */
#line 377 "AB13CD.f"
	}
#line 378 "AB13CD.f"
/* L110: */
#line 378 "AB13CD.f"
    }
#line 379 "AB13CD.f"
    zgesv_(n, m, &cwork[1], n, &iwork[1], &cwork[icw2 + 1], n, &info2);
#line 380 "AB13CD.f"
    if (info2 > 0) {
#line 381 "AB13CD.f"
	*info = 1;
#line 382 "AB13CD.f"
	return ret_val;
#line 383 "AB13CD.f"
    }
#line 384 "AB13CD.f"
    zgemm_("N", "N", np, m, n, &c_b1, &cwork[icw3 + 1], np, &cwork[icw2 + 1], 
	    n, &c_b1, &cwork[icw4 + 1], np, (ftnlen)1, (ftnlen)1);
#line 386 "AB13CD.f"
    i__1 = *lcwork - icwrk;
#line 386 "AB13CD.f"
    zgesvd_("N", "N", np, m, &cwork[icw4 + 1], np, &dwork[iw6 + 1], &cwork[1],
	     np, &cwork[1], m, &cwork[icwrk + 1], &i__1, &dwork[iwrk + 1], &
	    info2, (ftnlen)1, (ftnlen)1);
#line 389 "AB13CD.f"
    if (info2 > 0) {
#line 390 "AB13CD.f"
	*info = 4;
#line 391 "AB13CD.f"
	return ret_val;
#line 392 "AB13CD.f"
    }
#line 393 "AB13CD.f"
    if (gammal < dwork[iw6 + 1]) {
#line 394 "AB13CD.f"
	gammal = dwork[iw6 + 1];
#line 395 "AB13CD.f"
	fpeak = omega;
#line 396 "AB13CD.f"
    }
#line 397 "AB13CD.f"
    i__1 = icwrk + 1;
#line 397 "AB13CD.f"
    lcwamx = (integer) cwork[i__1].r + icwrk;

/*     Workspace usage. */

#line 401 "AB13CD.f"
    iw2 = *m * *n;
#line 402 "AB13CD.f"
    iw3 = iw2 + *m * *m;
#line 403 "AB13CD.f"
    iw4 = iw3 + *np * *np;
#line 404 "AB13CD.f"
    iw5 = iw4 + *m * *m;
#line 405 "AB13CD.f"
    iw6 = iw5 + *m * *n;
#line 406 "AB13CD.f"
    iw7 = iw6 + *m * *n;
#line 407 "AB13CD.f"
    iw8 = iw7 + *np * *np;
#line 408 "AB13CD.f"
    iw9 = iw8 + *np * *n;
#line 409 "AB13CD.f"
    iw10 = iw9 + (*n << 2) * *n;
#line 410 "AB13CD.f"
    iw11 = iw10 + (*n << 1);
#line 411 "AB13CD.f"
    iw12 = iw11 + (*n << 1);
#line 412 "AB13CD.f"
    iwrk = iw12 + min(*np,*m);

/*     Compute D'*C . */

#line 416 "AB13CD.f"
    dgemm_("T", "N", m, n, np, &c_b17, &d__[d_offset], ldd, &c__[c_offset], 
	    ldc, &c_b38, &dwork[1], m, (ftnlen)1, (ftnlen)1);

/*     Compute D'*D . */

#line 421 "AB13CD.f"
    dsyrk_("U", "T", m, np, &c_b17, &d__[d_offset], ldd, &c_b38, &dwork[iw2 + 
	    1], m, (ftnlen)1, (ftnlen)1);

/*     Compute D*D' . */

#line 426 "AB13CD.f"
    dsyrk_("U", "N", np, m, &c_b17, &d__[d_offset], ldd, &c_b38, &dwork[iw3 + 
	    1], np, (ftnlen)1, (ftnlen)1);

/*     Main iteration loop for gamma. */

#line 431 "AB13CD.f"
    iter = 0;
#line 432 "AB13CD.f"
L120:
#line 432 "AB13CD.f"
    ++iter;
#line 433 "AB13CD.f"
    if (iter > 10) {
#line 434 "AB13CD.f"
	*info = 2;
#line 435 "AB13CD.f"
	return ret_val;
#line 436 "AB13CD.f"
    }
#line 437 "AB13CD.f"
    gamma = (*tol * 2. + 1.) * gammal;

/*     Compute R = GAMMA^2*Im - D'*D . */

#line 441 "AB13CD.f"
    i__1 = *m;
#line 441 "AB13CD.f"
    for (j = 1; j <= i__1; ++j) {
#line 442 "AB13CD.f"
	i__2 = j;
#line 442 "AB13CD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 443 "AB13CD.f"
	    dwork[iw4 + i__ + (j - 1) * *m] = -dwork[iw2 + i__ + (j - 1) * *m]
		    ;
#line 444 "AB13CD.f"
/* L130: */
#line 444 "AB13CD.f"
	}
/* Computing 2nd power */
#line 445 "AB13CD.f"
	d__1 = gamma;
#line 445 "AB13CD.f"
	dwork[iw4 + j + (j - 1) * *m] = d__1 * d__1 - dwork[iw2 + j + (j - 1) 
		* *m];
#line 446 "AB13CD.f"
/* L140: */
#line 446 "AB13CD.f"
    }

/*     Compute inv(R)*D'*C . */

#line 450 "AB13CD.f"
    dlacpy_("Full", m, n, &dwork[1], m, &dwork[iw5 + 1], m, (ftnlen)4);
#line 451 "AB13CD.f"
    dpotrf_("U", m, &dwork[iw4 + 1], m, &info2, (ftnlen)1);
#line 452 "AB13CD.f"
    if (info2 > 0) {
#line 453 "AB13CD.f"
	*info = 2;
#line 454 "AB13CD.f"
	return ret_val;
#line 455 "AB13CD.f"
    }
#line 456 "AB13CD.f"
    dpotrs_("U", m, n, &dwork[iw4 + 1], m, &dwork[iw5 + 1], m, &info2, (
	    ftnlen)1);

/*     Compute inv(R)*B' . */

#line 461 "AB13CD.f"
    i__1 = *n;
#line 461 "AB13CD.f"
    for (j = 1; j <= i__1; ++j) {
#line 462 "AB13CD.f"
	i__2 = *m;
#line 462 "AB13CD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 463 "AB13CD.f"
	    dwork[iw6 + i__ + (j - 1) * *m] = b[j + i__ * b_dim1];
#line 464 "AB13CD.f"
/* L150: */
#line 464 "AB13CD.f"
	}
#line 465 "AB13CD.f"
/* L160: */
#line 465 "AB13CD.f"
    }
#line 466 "AB13CD.f"
    dpotrs_("U", m, n, &dwork[iw4 + 1], m, &dwork[iw6 + 1], m, &info2, (
	    ftnlen)1);

/*     Compute S = GAMMA^2*Ip - D*D' . */

#line 471 "AB13CD.f"
    i__1 = *np;
#line 471 "AB13CD.f"
    for (j = 1; j <= i__1; ++j) {
#line 472 "AB13CD.f"
	i__2 = j;
#line 472 "AB13CD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 473 "AB13CD.f"
	    dwork[iw7 + i__ + (j - 1) * *np] = -dwork[iw3 + i__ + (j - 1) * *
		    np];
#line 474 "AB13CD.f"
/* L170: */
#line 474 "AB13CD.f"
	}
/* Computing 2nd power */
#line 475 "AB13CD.f"
	d__1 = gamma;
#line 475 "AB13CD.f"
	dwork[iw7 + j + (j - 1) * *np] = d__1 * d__1 - dwork[iw3 + j + (j - 1)
		 * *np];
#line 476 "AB13CD.f"
/* L180: */
#line 476 "AB13CD.f"
    }

/*     Compute inv(S)*C . */

#line 480 "AB13CD.f"
    dlacpy_("Full", np, n, &c__[c_offset], ldc, &dwork[iw8 + 1], np, (ftnlen)
	    4);
#line 481 "AB13CD.f"
    dposv_("U", np, n, &dwork[iw7 + 1], np, &dwork[iw8 + 1], np, &info2, (
	    ftnlen)1);
#line 483 "AB13CD.f"
    if (info2 > 0) {
#line 484 "AB13CD.f"
	*info = 2;
#line 485 "AB13CD.f"
	return ret_val;
#line 486 "AB13CD.f"
    }

/*     Construct the Hamiltonian matrix . */

#line 490 "AB13CD.f"
    i__1 = *n << 1;
#line 490 "AB13CD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[iw9 + 1], &i__1, (ftnlen)
	    4);
#line 491 "AB13CD.f"
    i__1 = *n << 1;
#line 491 "AB13CD.f"
    dgemm_("N", "N", n, n, m, &c_b17, &b[b_offset], ldb, &dwork[iw5 + 1], m, &
	    c_b17, &dwork[iw9 + 1], &i__1, (ftnlen)1, (ftnlen)1);
#line 493 "AB13CD.f"
    d__1 = -gamma;
#line 493 "AB13CD.f"
    i__1 = *n << 1;
#line 493 "AB13CD.f"
    mb01rx_("Left", "Upper", "Transpose", n, np, &c_b38, &d__1, &dwork[iw9 + *
	    n + 1], &i__1, &c__[c_offset], ldc, &dwork[iw8 + 1], np, &info2, (
	    ftnlen)4, (ftnlen)5, (ftnlen)9);
#line 496 "AB13CD.f"
    i__1 = *n << 1;
#line 496 "AB13CD.f"
    ma02ed_("Upper", n, &dwork[iw9 + *n + 1], &i__1, (ftnlen)5);
#line 497 "AB13CD.f"
    i__1 = *n << 1;
#line 497 "AB13CD.f"
    mb01rx_("Left", "Upper", "NoTranspose", n, m, &c_b38, &gamma, &dwork[iw9 
	    + (*n << 1) * *n + 1], &i__1, &b[b_offset], ldb, &dwork[iw6 + 1], 
	    m, &info2, (ftnlen)4, (ftnlen)5, (ftnlen)11);
#line 500 "AB13CD.f"
    i__1 = *n << 1;
#line 500 "AB13CD.f"
    ma02ed_("Upper", n, &dwork[iw9 + (*n << 1) * *n + 1], &i__1, (ftnlen)5);
#line 501 "AB13CD.f"
    i__1 = *n;
#line 501 "AB13CD.f"
    for (j = 1; j <= i__1; ++j) {
#line 502 "AB13CD.f"
	i__2 = *n;
#line 502 "AB13CD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 503 "AB13CD.f"
	    dwork[iw9 + (*n << 1) * *n + *n + i__ + (j - 1 << 1) * *n] = 
		    -dwork[iw9 + j + (i__ - 1 << 1) * *n];
#line 504 "AB13CD.f"
/* L190: */
#line 504 "AB13CD.f"
	}
#line 505 "AB13CD.f"
/* L200: */
#line 505 "AB13CD.f"
    }

/*     Compute the eigenvalues of the Hamiltonian matrix. */

#line 509 "AB13CD.f"
    i__1 = *n << 1;
#line 509 "AB13CD.f"
    i__2 = *n << 1;
#line 509 "AB13CD.f"
    i__3 = *n << 1;
#line 509 "AB13CD.f"
    i__4 = *ldwork - iwrk;
#line 509 "AB13CD.f"
    dgees_("N", "S", (L_fp)sb02cx_, &i__1, &dwork[iw9 + 1], &i__2, &sdim, &
	    dwork[iw10 + 1], &dwork[iw11 + 1], &dwork[1], &i__3, &dwork[iwrk 
	    + 1], &i__4, &bwork[1], &info2, (ftnlen)1, (ftnlen)1);
#line 512 "AB13CD.f"
    if (info2 > 0) {
#line 513 "AB13CD.f"
	*info = 3;
#line 514 "AB13CD.f"
	return ret_val;
#line 515 "AB13CD.f"
    }
/* Computing MAX */
#line 516 "AB13CD.f"
    i__1 = (integer) dwork[iwrk + 1] + iwrk;
#line 516 "AB13CD.f"
    lwamax = max(i__1,lwamax);

#line 518 "AB13CD.f"
    if (sdim == 0) {
#line 519 "AB13CD.f"
	gammau = gamma;
#line 520 "AB13CD.f"
	goto L330;
#line 521 "AB13CD.f"
    }

/*     Store the positive imaginary parts. */

#line 525 "AB13CD.f"
    j = 0;
#line 526 "AB13CD.f"
    i__1 = sdim - 1;
#line 526 "AB13CD.f"
    for (i__ = 1; i__ <= i__1; i__ += 2) {
#line 527 "AB13CD.f"
	++j;
#line 528 "AB13CD.f"
	dwork[iw10 + j] = dwork[iw11 + i__];
#line 529 "AB13CD.f"
/* L210: */
#line 529 "AB13CD.f"
    }
#line 530 "AB13CD.f"
    k = j;

#line 532 "AB13CD.f"
    if (k >= 2) {

/*        Reorder the imaginary parts. */

#line 536 "AB13CD.f"
	i__1 = k - 1;
#line 536 "AB13CD.f"
	for (j = 1; j <= i__1; ++j) {
#line 537 "AB13CD.f"
	    i__2 = k;
#line 537 "AB13CD.f"
	    for (l = j + 1; l <= i__2; ++l) {
#line 538 "AB13CD.f"
		if (dwork[iw10 + j] <= dwork[iw10 + l]) {
#line 538 "AB13CD.f"
		    goto L220;
#line 538 "AB13CD.f"
		}
#line 539 "AB13CD.f"
		temp = dwork[iw10 + j];
#line 540 "AB13CD.f"
		dwork[iw10 + j] = dwork[iw10 + l];
#line 541 "AB13CD.f"
		dwork[iw10 + l] = temp;
#line 542 "AB13CD.f"
L220:
#line 542 "AB13CD.f"
		;
#line 542 "AB13CD.f"
	    }
#line 543 "AB13CD.f"
/* L230: */
#line 543 "AB13CD.f"
	}

/*        Determine the next frequency. */

#line 547 "AB13CD.f"
	i__1 = k - 1;
#line 547 "AB13CD.f"
	for (l = 1; l <= i__1; ++l) {
#line 548 "AB13CD.f"
	    omega = (dwork[iw10 + l] + dwork[iw10 + l + 1]) / 2.;
#line 549 "AB13CD.f"
	    i__2 = *n;
#line 549 "AB13CD.f"
	    for (j = 1; j <= i__2; ++j) {
#line 550 "AB13CD.f"
		i__3 = *n;
#line 550 "AB13CD.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 551 "AB13CD.f"
		    i__4 = i__ + (j - 1) * *n;
#line 551 "AB13CD.f"
		    d__1 = -a[i__ + j * a_dim1];
#line 551 "AB13CD.f"
		    cwork[i__4].r = d__1, cwork[i__4].i = 0.;
#line 552 "AB13CD.f"
/* L240: */
#line 552 "AB13CD.f"
		}
#line 553 "AB13CD.f"
		i__3 = j + (j - 1) * *n;
#line 553 "AB13CD.f"
		z__2.r = omega * 0., z__2.i = omega * 1.;
#line 553 "AB13CD.f"
		i__4 = j + j * a_dim1;
#line 553 "AB13CD.f"
		z__1.r = z__2.r - a[i__4], z__1.i = z__2.i;
#line 553 "AB13CD.f"
		cwork[i__3].r = z__1.r, cwork[i__3].i = z__1.i;
#line 554 "AB13CD.f"
/* L250: */
#line 554 "AB13CD.f"
	    }
#line 555 "AB13CD.f"
	    i__2 = *m;
#line 555 "AB13CD.f"
	    for (j = 1; j <= i__2; ++j) {
#line 556 "AB13CD.f"
		i__3 = *n;
#line 556 "AB13CD.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 557 "AB13CD.f"
		    i__4 = icw2 + i__ + (j - 1) * *n;
#line 557 "AB13CD.f"
		    i__5 = i__ + j * b_dim1;
#line 557 "AB13CD.f"
		    cwork[i__4].r = b[i__5], cwork[i__4].i = 0.;
#line 558 "AB13CD.f"
/* L260: */
#line 558 "AB13CD.f"
		}
#line 559 "AB13CD.f"
/* L270: */
#line 559 "AB13CD.f"
	    }
#line 560 "AB13CD.f"
	    i__2 = *n;
#line 560 "AB13CD.f"
	    for (j = 1; j <= i__2; ++j) {
#line 561 "AB13CD.f"
		i__3 = *np;
#line 561 "AB13CD.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 562 "AB13CD.f"
		    i__4 = icw3 + i__ + (j - 1) * *np;
#line 562 "AB13CD.f"
		    i__5 = i__ + j * c_dim1;
#line 562 "AB13CD.f"
		    cwork[i__4].r = c__[i__5], cwork[i__4].i = 0.;
#line 563 "AB13CD.f"
/* L280: */
#line 563 "AB13CD.f"
		}
#line 564 "AB13CD.f"
/* L290: */
#line 564 "AB13CD.f"
	    }
#line 565 "AB13CD.f"
	    i__2 = *m;
#line 565 "AB13CD.f"
	    for (j = 1; j <= i__2; ++j) {
#line 566 "AB13CD.f"
		i__3 = *np;
#line 566 "AB13CD.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 567 "AB13CD.f"
		    i__4 = icw4 + i__ + (j - 1) * *np;
#line 567 "AB13CD.f"
		    i__5 = i__ + j * d_dim1;
#line 567 "AB13CD.f"
		    cwork[i__4].r = d__[i__5], cwork[i__4].i = 0.;
#line 568 "AB13CD.f"
/* L300: */
#line 568 "AB13CD.f"
		}
#line 569 "AB13CD.f"
/* L310: */
#line 569 "AB13CD.f"
	    }
#line 570 "AB13CD.f"
	    zgesv_(n, m, &cwork[1], n, &iwork[1], &cwork[icw2 + 1], n, &info2)
		    ;
#line 572 "AB13CD.f"
	    if (info2 > 0) {
#line 573 "AB13CD.f"
		*info = 1;
#line 574 "AB13CD.f"
		return ret_val;
#line 575 "AB13CD.f"
	    }
#line 576 "AB13CD.f"
	    zgemm_("N", "N", np, m, n, &c_b1, &cwork[icw3 + 1], np, &cwork[
		    icw2 + 1], n, &c_b1, &cwork[icw4 + 1], np, (ftnlen)1, (
		    ftnlen)1);
#line 578 "AB13CD.f"
	    i__2 = *lcwork - icwrk;
#line 578 "AB13CD.f"
	    zgesvd_("N", "N", np, m, &cwork[icw4 + 1], np, &dwork[iw6 + 1], &
		    cwork[1], np, &cwork[1], m, &cwork[icwrk + 1], &i__2, &
		    dwork[iwrk + 1], &info2, (ftnlen)1, (ftnlen)1);
#line 582 "AB13CD.f"
	    if (info2 > 0) {
#line 583 "AB13CD.f"
		*info = 4;
#line 584 "AB13CD.f"
		return ret_val;
#line 585 "AB13CD.f"
	    }
#line 586 "AB13CD.f"
	    if (gammal < dwork[iw6 + 1]) {
#line 587 "AB13CD.f"
		gammal = dwork[iw6 + 1];
#line 588 "AB13CD.f"
		fpeak = omega;
#line 589 "AB13CD.f"
	    }
/* Computing MAX */
#line 590 "AB13CD.f"
	    i__3 = icwrk + 1;
#line 590 "AB13CD.f"
	    i__2 = (integer) cwork[i__3].r + icwrk;
#line 590 "AB13CD.f"
	    lcwamx = max(i__2,lcwamx);
#line 591 "AB13CD.f"
/* L320: */
#line 591 "AB13CD.f"
	}
#line 592 "AB13CD.f"
    }
#line 593 "AB13CD.f"
    goto L120;
#line 594 "AB13CD.f"
L330:
#line 594 "AB13CD.f"
    ret_val = (gammal + gammau) / 2.;

#line 596 "AB13CD.f"
    dwork[1] = (doublereal) lwamax;
#line 597 "AB13CD.f"
    dwork[2] = fpeak;
#line 598 "AB13CD.f"
    cwork[1].r = (doublereal) lcwamx, cwork[1].i = 0.;
#line 599 "AB13CD.f"
    return ret_val;
/* *** End of AB13CD *** */
} /* ab13cd_ */

