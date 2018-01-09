#line 1 "AB13FD.f"
/* AB13FD.f -- translated by f2c (version 20100827).
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

#line 1 "AB13FD.f"
/* Table of constant values */

static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b21 = 1.;
static doublereal c_b22 = 0.;

/* Subroutine */ int ab13fd_(integer *n, doublereal *a, integer *lda, 
	doublereal *beta, doublereal *omega, doublereal *tol, doublereal *
	dwork, integer *ldwork, doublecomplex *cwork, integer *lcwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, n2;
    static doublereal om, sv;
    static integer ia2;
    static doublereal om1, om2;
    static integer iaa, igf, ihi, ilo, kom;
    static doublereal eps;
    static integer iwi;
    static doublereal tau;
    static integer iwk, iwr;
    static doublereal low, tol1, sfmn, temp;
    extern /* Subroutine */ int ma02ed_(char *, integer *, doublereal *, 
	    integer *, ftnlen), dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal sigma;
    extern /* Subroutine */ int mb04zd_(char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen);
    extern doublereal mb03ny_(integer *, doublereal *, doublereal *, integer *
	    , doublereal *, doublereal *, integer *, doublecomplex *, integer 
	    *, integer *);
    static integer lbest;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dsymm_(char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal dummy[1];
    static integer itnum, jwork;
    extern /* Subroutine */ int dsymv_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    static doublereal dummy2[1]	/* was [1][1] */;
    extern /* Subroutine */ int dgebal_(char *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dhseqr_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, ftnlen, ftnlen);
    static integer minwrk;
    static logical sufwrk;


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

/*     To compute beta(A), the 2-norm distance from a real matrix A to */
/*     the nearest complex matrix with an eigenvalue on the imaginary */
/*     axis. If A is stable in the sense that all eigenvalues of A lie */
/*     in the open left half complex plane, then beta(A) is the complex */
/*     stability radius, i.e., the distance to the nearest unstable */
/*     complex matrix. The value of beta(A) is the minimum of the */
/*     smallest singular value of (A - jwI), taken over all real w. */
/*     The value of w corresponding to the minimum is also computed. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     BETA    (output) DOUBLE PRECISION */
/*             The computed value of beta(A), which actually is an upper */
/*             bound. */

/*     OMEGA   (output) DOUBLE PRECISION */
/*             The value of w such that the smallest singular value of */
/*             (A - jwI) equals beta(A). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             Specifies the accuracy with which beta(A) is to be */
/*             calculated. (See the Numerical Aspects section below.) */
/*             If the user sets TOL to be less than EPS, where EPS is the */
/*             machine precision (see LAPACK Library Routine DLAMCH), */
/*             then the tolerance is taken to be EPS. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */
/*             If DWORK(1) is not needed, the first 2*N*N entries of */
/*             DWORK may overlay CWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 1, 3*N*(N+2) ). */
/*             For optimum performance LDWORK should be larger. */

/*     CWORK   COMPLEX*16 array, dimension (LCWORK) */
/*             On exit, if INFO = 0, CWORK(1) returns the optimal value */
/*             of LCWORK. */
/*             If CWORK(1) is not needed, the first N*N entries of */
/*             CWORK may overlay DWORK. */

/*     LCWORK  INTEGER */
/*             The length of the array CWORK. */
/*             LCWORK >= MAX( 1, N*(N+3) ). */
/*             For optimum performance LCWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the routine fails to compute beta(A) within the */
/*                   specified tolerance. Nevertheless, the returned */
/*                   value is an upper bound on beta(A); */
/*             = 2:  either the QR or SVD algorithm (LAPACK Library */
/*                   routines DHSEQR, DGESVD or ZGESVD) fails to */
/*                   converge; this error is very rare. */

/*     METHOD */

/*     AB13FD combines the methods of [1] and [2] into a provably */
/*     reliable, quadratically convergent algorithm. It uses the simple */
/*     bisection strategy of [1] to find an interval which contains */
/*     beta(A), and then switches to the modified bisection strategy of */
/*     [2] which converges quadratically to a minimizer. Note that the */
/*     efficiency of the strategy degrades if there are several local */
/*     minima that are near or equal the global minimum. */

/*     REFERENCES */

/*     [1] Byers, R. */
/*         A bisection method for measuring the distance of a stable */
/*         matrix to the unstable matrices. */
/*         SIAM J. Sci. Stat. Comput., Vol. 9, No. 5, pp. 875-880, 1988. */

/*     [2] Boyd, S. and Balakrishnan, K. */
/*         A regularity result for the singular values of a transfer */
/*         matrix and a quadratically convergent algorithm for computing */
/*         its L-infinity norm. */
/*         Systems and Control Letters, Vol. 15, pp. 1-7, 1990. */

/*     NUMERICAL ASPECTS */

/*     In the presence of rounding errors, the computed function value */
/*     BETA  satisfies */

/*           beta(A) <= BETA + epsilon, */

/*           BETA/(1+TOL) - delta <= MAX(beta(A), SQRT(2*N*EPS)*norm(A)), */

/*     where norm(A) is the Frobenius norm of A, */

/*           epsilon = p(N) * EPS * norm(A), */
/*     and */
/*           delta   = p(N) * SQRT(EPS) * norm(A), */

/*     and p(N) is a low degree polynomial. It is recommended to choose */
/*     TOL greater than SQRT(EPS). Although rounding errors can cause */
/*     AB13FD to fail for smaller values of TOL, nevertheless, it usually */
/*     succeeds. Regardless of success or failure, the first inequality */
/*     holds. */

/*     CONTRIBUTORS */

/*     R. Byers, the routines QSEC and QSEC0 (January, 1995). */

/*     REVISIONS */

/*     Release 4.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1999. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2002, */
/*     Jan. 2003. */

/*     KEYWORDS */

/*     complex stability radius, distances, eigenvalue, eigenvalue */
/*     perturbation, norms. */

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

/*     Test the input scalar arguments. */

#line 195 "AB13FD.f"
    /* Parameter adjustments */
#line 195 "AB13FD.f"
    a_dim1 = *lda;
#line 195 "AB13FD.f"
    a_offset = 1 + a_dim1;
#line 195 "AB13FD.f"
    a -= a_offset;
#line 195 "AB13FD.f"
    --dwork;
#line 195 "AB13FD.f"
    --cwork;
#line 195 "AB13FD.f"

#line 195 "AB13FD.f"
    /* Function Body */
#line 195 "AB13FD.f"
    *info = 0;
#line 196 "AB13FD.f"
    minwrk = *n * 3 * (*n + 2);

#line 198 "AB13FD.f"
    if (*n < 0) {
#line 199 "AB13FD.f"
	*info = -1;
#line 200 "AB13FD.f"
    } else if (*lda < max(1,*n)) {
#line 201 "AB13FD.f"
	*info = -3;
#line 202 "AB13FD.f"
    } else if (*ldwork < max(1,minwrk)) {
#line 203 "AB13FD.f"
	*info = -8;
#line 204 "AB13FD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 204 "AB13FD.f"
	i__1 = 1, i__2 = *n * (*n + 3);
#line 204 "AB13FD.f"
	if (*lcwork < max(i__1,i__2)) {
#line 205 "AB13FD.f"
	    *info = -10;
#line 206 "AB13FD.f"
	}
#line 206 "AB13FD.f"
    }

#line 208 "AB13FD.f"
    if (*info != 0) {

/*        Error return. */

#line 212 "AB13FD.f"
	i__1 = -(*info);
#line 212 "AB13FD.f"
	xerbla_("AB13FD", &i__1, (ftnlen)6);
#line 213 "AB13FD.f"
	return 0;
#line 214 "AB13FD.f"
    }

/*     Quick return if possible. */

#line 218 "AB13FD.f"
    *omega = 0.;
#line 219 "AB13FD.f"
    if (*n == 0) {
#line 220 "AB13FD.f"
	*beta = 0.;
#line 221 "AB13FD.f"
	dwork[1] = 1.;
#line 222 "AB13FD.f"
	cwork[1].r = 1., cwork[1].i = 0.;
#line 223 "AB13FD.f"
	return 0;
#line 224 "AB13FD.f"
    }

/*     Indices for splitting the work array. */
/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of workspace needed at that point in the code, */
/*     as well as the preferred amount for good performance.) */

#line 231 "AB13FD.f"
    n2 = *n * *n;
#line 232 "AB13FD.f"
    igf = 1;
#line 233 "AB13FD.f"
    ia2 = igf + n2 + *n;
#line 234 "AB13FD.f"
    iaa = ia2 + n2;
#line 235 "AB13FD.f"
    iwk = iaa + n2;
#line 236 "AB13FD.f"
    iwr = iaa;
#line 237 "AB13FD.f"
    iwi = iwr + *n;

#line 239 "AB13FD.f"
    sufwrk = *ldwork - iwk >= n2;

/*     Computation of the tolerances. EPS is the machine precision. */

#line 243 "AB13FD.f"
    sfmn = dlamch_("Safe minimum", (ftnlen)12);
#line 244 "AB13FD.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 245 "AB13FD.f"
    tol1 = sqrt(eps * (doublereal) (*n << 1)) * dlange_("Frobenius", n, n, &a[
	    a_offset], lda, &dwork[1], (ftnlen)9);
#line 247 "AB13FD.f"
    tau = max(*tol,eps) + 1.;

/*     Initialization, upper bound at known critical point. */
/*     Workspace: need N*(N+1)+5*N; prefer larger. */

#line 252 "AB13FD.f"
    kom = 2;
#line 253 "AB13FD.f"
    low = 0.;
#line 254 "AB13FD.f"
    dlacpy_("All", n, n, &a[a_offset], lda, &dwork[igf], n, (ftnlen)3);
#line 255 "AB13FD.f"
    i__1 = *ldwork - ia2;
#line 255 "AB13FD.f"
    *beta = mb03ny_(n, omega, &dwork[igf], n, &dwork[igf + n2], &dwork[ia2], &
	    i__1, &cwork[1], lcwork, info);
#line 257 "AB13FD.f"
    if (*info != 0) {
#line 257 "AB13FD.f"
	return 0;
#line 257 "AB13FD.f"
    }
/* Computing MAX */
#line 259 "AB13FD.f"
    i__1 = minwrk, i__2 = (integer) dwork[ia2] - ia2 + 1, i__1 = max(i__1,
	    i__2), i__2 = (n2 << 2) + *n;
#line 259 "AB13FD.f"
    lbest = max(i__1,i__2);

#line 261 "AB13FD.f"
    itnum = 1;
/*     WHILE ( ITNUM <= MAXIT and BETA > TAU*MAX( TOL1, LOW ) ) DO */
#line 263 "AB13FD.f"
L10:
#line 263 "AB13FD.f"
    if (itnum <= 50 && *beta > tau * max(tol1,low)) {
#line 265 "AB13FD.f"
	if (kom == 2) {
#line 266 "AB13FD.f"
	    sigma = *beta / tau;
#line 267 "AB13FD.f"
	} else {
#line 268 "AB13FD.f"
	    sigma = sqrt(*beta) * sqrt((max(tol1,low)));
#line 269 "AB13FD.f"
	}

/*        Set up H(sigma). */
/*        Workspace: N*(N+1)+2*N*N. */

#line 274 "AB13FD.f"
	dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[iaa], n, (ftnlen)4);
#line 275 "AB13FD.f"
	dwork[igf] = sigma;
#line 276 "AB13FD.f"
	dwork[igf + *n] = -sigma;
#line 277 "AB13FD.f"
	dummy[0] = 0.;
#line 278 "AB13FD.f"
	i__1 = *n - 1;
#line 278 "AB13FD.f"
	dcopy_(&i__1, dummy, &c__0, &dwork[igf + 1], &c__1);

#line 280 "AB13FD.f"
	i__1 = ia2 - *n - 2;
#line 280 "AB13FD.f"
	i__2 = *n + 1;
#line 280 "AB13FD.f"
	for (i__ = igf; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 281 "AB13FD.f"
	    i__3 = *n + 1;
#line 281 "AB13FD.f"
	    dcopy_(&i__3, &dwork[i__], &c__1, &dwork[i__ + *n + 1], &c__1);
#line 282 "AB13FD.f"
/* L20: */
#line 282 "AB13FD.f"
	}

/*        Computation of the eigenvalues by the square reduced algorithm. */
/*        Workspace: N*(N+1)+2*N*N+2*N. */

#line 287 "AB13FD.f"
	mb04zd_("No vectors", n, &dwork[iaa], n, &dwork[igf], n, dummy2, &
		c__1, &dwork[iwk], info, (ftnlen)10);

/*        Form the matrix A*A + F*G. */
/*        Workspace: need   N*(N+1)+2*N*N+N; */
/*                   prefer N*(N+1)+3*N*N. */

#line 294 "AB13FD.f"
	jwork = ia2;
#line 295 "AB13FD.f"
	if (sufwrk) {
#line 295 "AB13FD.f"
	    jwork = iwk;
#line 295 "AB13FD.f"
	}

#line 298 "AB13FD.f"
	dlacpy_("Lower", n, n, &dwork[igf], n, &dwork[jwork], n, (ftnlen)5);
#line 299 "AB13FD.f"
	ma02ed_("Lower", n, &dwork[jwork], n, (ftnlen)5);

#line 301 "AB13FD.f"
	if (sufwrk) {

/*           Use BLAS 3 calculation. */

#line 305 "AB13FD.f"
	    dsymm_("Left", "Upper", n, n, &c_b21, &dwork[igf + *n], n, &dwork[
		    jwork], n, &c_b22, &dwork[ia2], n, (ftnlen)4, (ftnlen)5);
#line 307 "AB13FD.f"
	} else {

/*           Use BLAS 2 calculation. */

#line 311 "AB13FD.f"
	    i__2 = *n;
#line 311 "AB13FD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 312 "AB13FD.f"
		dsymv_("Upper", n, &c_b21, &dwork[igf + *n], n, &dwork[ia2 + *
			n * (i__ - 1)], &c__1, &c_b22, &dwork[iwk], &c__1, (
			ftnlen)5);
#line 314 "AB13FD.f"
		dcopy_(n, &dwork[iwk], &c__1, &dwork[ia2 + *n * (i__ - 1)], &
			c__1);
#line 315 "AB13FD.f"
/* L30: */
#line 315 "AB13FD.f"
	    }

#line 317 "AB13FD.f"
	}

#line 319 "AB13FD.f"
	dgemm_("NoTranspose", "NoTranspose", n, n, n, &c_b21, &dwork[iaa], n, 
		&dwork[iaa], n, &c_b21, &dwork[ia2], n, (ftnlen)11, (ftnlen)
		11);

/*        Find the eigenvalues of A*A + F*G. */
/*        Workspace: N*(N+1)+N*N+3*N. */

#line 325 "AB13FD.f"
	jwork = iwi + *n;
#line 326 "AB13FD.f"
	dgebal_("Scale", n, &dwork[ia2], n, &ilo, &ihi, &dwork[jwork], &i__, (
		ftnlen)5);
#line 328 "AB13FD.f"
	dhseqr_("Eigenvalues", "NoSchurVectors", n, &ilo, &ihi, &dwork[ia2], 
		n, &dwork[iwr], &dwork[iwi], dummy2, &c__1, &dwork[jwork], n, 
		info, (ftnlen)11, (ftnlen)14);

#line 332 "AB13FD.f"
	if (*info != 0) {
#line 333 "AB13FD.f"
	    *info = 2;
#line 334 "AB13FD.f"
	    return 0;
#line 335 "AB13FD.f"
	}

/*        Count negative real axis squared eigenvalues. If there are two, */
/*        then the valley is isolated, and next approximate minimizer is */
/*        mean of the square roots. */

#line 341 "AB13FD.f"
	kom = 0;
#line 342 "AB13FD.f"
	i__2 = *n - 1;
#line 342 "AB13FD.f"
	for (i__ = 0; i__ <= i__2; ++i__) {
#line 343 "AB13FD.f"
	    temp = (d__1 = dwork[iwi + i__], abs(d__1));
#line 344 "AB13FD.f"
	    if (tol1 > sfmn) {
#line 344 "AB13FD.f"
		temp /= tol1;
#line 344 "AB13FD.f"
	    }
#line 345 "AB13FD.f"
	    if (dwork[iwr + i__] < 0. && temp <= tol1) {
#line 346 "AB13FD.f"
		++kom;
#line 347 "AB13FD.f"
		om = sqrt(-dwork[iwr + i__]);
#line 348 "AB13FD.f"
		if (kom == 1) {
#line 348 "AB13FD.f"
		    om1 = om;
#line 348 "AB13FD.f"
		}
#line 349 "AB13FD.f"
		if (kom == 2) {
#line 349 "AB13FD.f"
		    om2 = om;
#line 349 "AB13FD.f"
		}
#line 350 "AB13FD.f"
	    }
#line 351 "AB13FD.f"
/* L40: */
#line 351 "AB13FD.f"
	}

#line 353 "AB13FD.f"
	if (kom == 0) {
#line 354 "AB13FD.f"
	    low = sigma;
#line 355 "AB13FD.f"
	} else {

/*           In exact arithmetic KOM = 1 is impossible, but if tau is */
/*           close enough to one, MB04ZD may miss the initial near zero */
/*           eigenvalue. */
/*           Workspace, real:    need   3*N*(N+2);  prefer larger; */
/*                      complex: need     N*(N+3);  prefer larger. */

#line 363 "AB13FD.f"
	    if (kom == 2) {
#line 364 "AB13FD.f"
		om = om1 + (om2 - om1) / 2.;
#line 365 "AB13FD.f"
	    } else if (kom == 1 && itnum == 1) {
#line 366 "AB13FD.f"
		om = om1 / 2.;
#line 367 "AB13FD.f"
		kom = 2;
#line 368 "AB13FD.f"
	    }

#line 370 "AB13FD.f"
	    dlacpy_("All", n, n, &a[a_offset], lda, &dwork[igf], n, (ftnlen)3)
		    ;
#line 371 "AB13FD.f"
	    i__2 = *ldwork - ia2;
#line 371 "AB13FD.f"
	    sv = mb03ny_(n, &om, &dwork[igf], n, &dwork[igf + n2], &dwork[ia2]
		    , &i__2, &cwork[1], lcwork, info);
#line 373 "AB13FD.f"
	    if (*info != 0) {
#line 373 "AB13FD.f"
		return 0;
#line 373 "AB13FD.f"
	    }
#line 375 "AB13FD.f"
	    if (*beta > sv) {
#line 376 "AB13FD.f"
		*beta = sv;
#line 377 "AB13FD.f"
		*omega = om;
#line 378 "AB13FD.f"
	    } else {
#line 379 "AB13FD.f"
		*info = 1;
#line 380 "AB13FD.f"
		return 0;
#line 381 "AB13FD.f"
	    }
#line 382 "AB13FD.f"
	}
#line 383 "AB13FD.f"
	++itnum;
#line 384 "AB13FD.f"
	goto L10;
/*        END WHILE 10 */
#line 386 "AB13FD.f"
    }

#line 388 "AB13FD.f"
    if (*beta > tau * max(tol1,low)) {

/*        Failed to meet bounds within MAXIT iterations. */

#line 392 "AB13FD.f"
	*info = 1;
#line 393 "AB13FD.f"
	return 0;
#line 394 "AB13FD.f"
    }

/*     Set optimal real workspace dimension (complex workspace is already */
/*     set by MB03NY). */

#line 399 "AB13FD.f"
    dwork[1] = (doublereal) lbest;

#line 401 "AB13FD.f"
    return 0;
/* *** Last line of AB13FD *** */
} /* ab13fd_ */
