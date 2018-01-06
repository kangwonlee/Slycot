#line 1 "AB13ED.f"
/* AB13ED.f -- translated by f2c (version 20100827).
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

#line 1 "AB13ED.f"
/* Table of constant values */

static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b19 = 1.;
static doublereal c_b20 = 0.;

/* Subroutine */ int ab13ed_(integer *n, doublereal *a, integer *lda, 
	doublereal *low, doublereal *high, doublereal *tol, doublereal *dwork,
	 integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, n2, ia2, iaa, igf, ihi, ilo, iwi;
    static doublereal tau;
    static integer iwk, iwr;
    static doublereal tol1, tol2;
    static logical rneg;
    static doublereal anrm, sfmn, temp, seps;
    extern /* Subroutine */ int ma02ed_(char *, integer *, doublereal *, 
	    integer *, ftnlen), dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal sigma;
    extern /* Subroutine */ int mb04zd_(char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dcopy_(integer *, doublereal *, 
	    integer *, doublereal *, integer *), dsymm_(char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen);
    static doublereal dummy[1];
    static integer jwork;
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

/*     To estimate beta(A), the 2-norm distance from a real matrix A to */
/*     the nearest complex matrix with an eigenvalue on the imaginary */
/*     axis. The estimate is given as */

/*            LOW <= beta(A) <= HIGH, */

/*     where either */

/*            (1 + TOL) * LOW >= HIGH, */

/*     or */

/*            LOW = 0   and   HIGH = delta, */

/*     and delta is a small number approximately equal to the square root */
/*     of machine precision times the Frobenius norm (Euclidean norm) */
/*     of A. If A is stable in the sense that all eigenvalues of A lie */
/*     in the open left half complex plane, then beta(A) is the distance */
/*     to the nearest unstable complex matrix, i.e., the complex */
/*     stability radius. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     LOW     (output) DOUBLE PRECISION */
/*             A lower bound for beta(A). */

/*     HIGH    (output) DOUBLE PRECISION */
/*             An upper bound for beta(A). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             Specifies the accuracy with which LOW and HIGH approximate */
/*             beta(A). If the user sets TOL to be less than SQRT(EPS), */
/*             where EPS is the machine precision (see LAPACK Library */
/*             Routine DLAMCH), then the tolerance is taken to be */
/*             SQRT(EPS). */
/*             The recommended value is TOL = 9, which gives an estimate */
/*             of beta(A) correct to within an order of magnitude. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 1, 3*N*(N+1) ). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the QR algorithm (LAPACK Library routine DHSEQR) */
/*                   fails to converge; this error is very rare. */

/*     METHOD */

/*     Let beta(A) be the 2-norm distance from a real matrix A to the */
/*     nearest complex matrix with an eigenvalue on the imaginary axis. */
/*     It is known that beta(A) = minimum of the smallest singular */
/*     value of (A - jwI), where I is the identity matrix and j**2 = -1, */
/*     and the minimum is taken over all real w. */
/*     The algorithm computes a lower bound LOW and an upper bound HIGH */
/*     for beta(A) by a bisection method in the following way. Given a */
/*     non-negative real number sigma, the Hamiltonian matrix H(sigma) */
/*     is constructed: */

/*                       |   A      -sigma*I |     | A   G  | */
/*           H(sigma) =  |                   | :=  |        | . */
/*                       | sigma*I    -A'    |     | F  -A' | */

/*     It can be shown [1] that H(sigma) has an eigenvalue whose real */
/*     part is zero if and only if sigma >= beta. Any lower and upper */
/*     bounds on beta(A) can be improved by choosing a number between */
/*     them and checking to see if H(sigma) has an eigenvalue with zero */
/*     real part.  This decision is made by computing the eigenvalues of */
/*     H(sigma) using the square reduced algorithm of Van Loan [2]. */

/*     REFERENCES */

/*     [1] Byers, R. */
/*         A bisection method for measuring the distance of a stable */
/*         matrix to the unstable matrices. */
/*         SIAM J. Sci. Stat. Comput., Vol. 9, No. 5, pp. 875-880, 1988. */

/*     [2] Van Loan, C.F. */
/*         A symplectic method for approximating all the eigenvalues of a */
/*         Hamiltonian matrix. */
/*         Linear Algebra and its Applications, Vol 61, 233-251, 1984. */

/*     NUMERICAL ASPECTS */

/*     Due to rounding errors the computed values of LOW and HIGH can be */
/*     proven to satisfy */

/*            LOW - p(n) * sqrt(e) * norm(A) <= beta(A) */
/*     and */
/*            beta(A) <= HIGH + p(n) * sqrt(e) * norm(A), */

/*     where p(n) is a modest polynomial of degree 3, e is the machine */
/*     precision and norm(A) is the Frobenius norm of A, see [1]. */
/*     The recommended value for TOL is 9 which gives an estimate of */
/*     beta(A) correct to within an order of magnitude. */
/*     AB13ED requires approximately 38*N**3 flops for TOL = 9. */

/*     CONTRIBUTOR */

/*     R. Byers, the routines BISEC and BISEC0 (January, 1995). */

/*     REVISIONS */

/*     Release 4.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1999. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Jan. 2003. */

/*     KEYWORDS */

/*     Distances, eigenvalue, eigenvalue perturbation, norms, stability */
/*     radius. */

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

#line 189 "AB13ED.f"
    /* Parameter adjustments */
#line 189 "AB13ED.f"
    a_dim1 = *lda;
#line 189 "AB13ED.f"
    a_offset = 1 + a_dim1;
#line 189 "AB13ED.f"
    a -= a_offset;
#line 189 "AB13ED.f"
    --dwork;
#line 189 "AB13ED.f"

#line 189 "AB13ED.f"
    /* Function Body */
#line 189 "AB13ED.f"
    *info = 0;
#line 190 "AB13ED.f"
    minwrk = *n * 3 * (*n + 1);

#line 192 "AB13ED.f"
    if (*n < 0) {
#line 193 "AB13ED.f"
	*info = -1;
#line 194 "AB13ED.f"
    } else if (*lda < max(1,*n)) {
#line 195 "AB13ED.f"
	*info = -3;
#line 196 "AB13ED.f"
    } else if (*ldwork < max(1,minwrk)) {
#line 197 "AB13ED.f"
	*info = -8;
#line 198 "AB13ED.f"
    }

#line 200 "AB13ED.f"
    if (*info != 0) {

/*        Error return. */

#line 204 "AB13ED.f"
	i__1 = -(*info);
#line 204 "AB13ED.f"
	xerbla_("AB13ED", &i__1, (ftnlen)6);
#line 205 "AB13ED.f"
	return 0;
#line 206 "AB13ED.f"
    }

/*     Quick return if possible. */

#line 210 "AB13ED.f"
    *low = 0.;
#line 211 "AB13ED.f"
    if (*n == 0) {
#line 212 "AB13ED.f"
	*high = 0.;
#line 213 "AB13ED.f"
	dwork[1] = 1.;
#line 214 "AB13ED.f"
	return 0;
#line 215 "AB13ED.f"
    }

/*     Indices for splitting the work array. */
/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance.) */

#line 222 "AB13ED.f"
    n2 = *n * *n;
#line 223 "AB13ED.f"
    igf = 1;
#line 224 "AB13ED.f"
    ia2 = igf + n2 + *n;
#line 225 "AB13ED.f"
    iaa = ia2 + n2;
#line 226 "AB13ED.f"
    iwk = iaa + n2;
#line 227 "AB13ED.f"
    iwr = iaa;
#line 228 "AB13ED.f"
    iwi = iwr + *n;

#line 230 "AB13ED.f"
    sufwrk = *ldwork - iwk >= n2;

/*     Computation of the tolerances and the treshold for termination of */
/*     the bisection method. SEPS is the square root of the machine */
/*     precision. */

#line 236 "AB13ED.f"
    sfmn = dlamch_("Safe minimum", (ftnlen)12);
#line 237 "AB13ED.f"
    seps = sqrt(dlamch_("Epsilon", (ftnlen)7));
#line 238 "AB13ED.f"
    tau = max(*tol,seps) + 1.;
#line 239 "AB13ED.f"
    anrm = dlange_("Frobenius", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)9)
	    ;
#line 240 "AB13ED.f"
    tol1 = seps * anrm;
#line 241 "AB13ED.f"
    tol2 = tol1 * (doublereal) (*n << 1);

/*     Initialization of the bisection method. */

#line 245 "AB13ED.f"
    *high = anrm;

/*     WHILE ( HIGH > TAU*MAX( TOL1, LOW ) ) DO */
#line 248 "AB13ED.f"
L10:
#line 248 "AB13ED.f"
    if (*high > tau * max(tol1,*low)) {
#line 249 "AB13ED.f"
	sigma = sqrt(*high) * sqrt((max(tol1,*low)));

/*        Set up H(sigma). */
/*        Workspace: N*(N+1)+2*N*N. */

#line 254 "AB13ED.f"
	dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[iaa], n, (ftnlen)4);
#line 255 "AB13ED.f"
	dwork[igf] = sigma;
#line 256 "AB13ED.f"
	dwork[igf + *n] = -sigma;
#line 257 "AB13ED.f"
	dummy[0] = 0.;
#line 258 "AB13ED.f"
	i__1 = *n - 1;
#line 258 "AB13ED.f"
	dcopy_(&i__1, dummy, &c__0, &dwork[igf + 1], &c__1);

#line 260 "AB13ED.f"
	i__1 = ia2 - *n - 2;
#line 260 "AB13ED.f"
	i__2 = *n + 1;
#line 260 "AB13ED.f"
	for (i__ = igf; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 261 "AB13ED.f"
	    i__3 = *n + 1;
#line 261 "AB13ED.f"
	    dcopy_(&i__3, &dwork[i__], &c__1, &dwork[i__ + *n + 1], &c__1);
#line 262 "AB13ED.f"
/* L20: */
#line 262 "AB13ED.f"
	}

/*        Computation of the eigenvalues by the square reduced algorithm. */
/*        Workspace: N*(N+1)+2*N*N+2*N. */

#line 267 "AB13ED.f"
	mb04zd_("No vectors", n, &dwork[iaa], n, &dwork[igf], n, dummy2, &
		c__1, &dwork[iwk], info, (ftnlen)10);

/*        Form the matrix A*A + F*G. */
/*        Workspace: need   N*(N+1)+2*N*N+N; */
/*                   prefer N*(N+1)+3*N*N. */

#line 274 "AB13ED.f"
	jwork = ia2;
#line 275 "AB13ED.f"
	if (sufwrk) {
#line 275 "AB13ED.f"
	    jwork = iwk;
#line 275 "AB13ED.f"
	}

#line 278 "AB13ED.f"
	dlacpy_("Lower", n, n, &dwork[igf], n, &dwork[jwork], n, (ftnlen)5);
#line 279 "AB13ED.f"
	ma02ed_("Lower", n, &dwork[jwork], n, (ftnlen)5);

#line 281 "AB13ED.f"
	if (sufwrk) {

/*           Use BLAS 3 calculation. */

#line 285 "AB13ED.f"
	    dsymm_("Left", "Upper", n, n, &c_b19, &dwork[igf + *n], n, &dwork[
		    jwork], n, &c_b20, &dwork[ia2], n, (ftnlen)4, (ftnlen)5);
#line 287 "AB13ED.f"
	} else {

/*           Use BLAS 2 calculation. */

#line 291 "AB13ED.f"
	    i__2 = *n;
#line 291 "AB13ED.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 292 "AB13ED.f"
		dsymv_("Upper", n, &c_b19, &dwork[igf + *n], n, &dwork[ia2 + *
			n * (i__ - 1)], &c__1, &c_b20, &dwork[iwk], &c__1, (
			ftnlen)5);
#line 294 "AB13ED.f"
		dcopy_(n, &dwork[iwk], &c__1, &dwork[ia2 + *n * (i__ - 1)], &
			c__1);
#line 295 "AB13ED.f"
/* L30: */
#line 295 "AB13ED.f"
	    }

#line 297 "AB13ED.f"
	}

#line 299 "AB13ED.f"
	dgemm_("NoTranspose", "NoTranspose", n, n, n, &c_b19, &dwork[iaa], n, 
		&dwork[iaa], n, &c_b19, &dwork[ia2], n, (ftnlen)11, (ftnlen)
		11);

/*        Find the eigenvalues of A*A + F*G. */
/*        Workspace: N*(N+1)+N*N+3*N. */

#line 305 "AB13ED.f"
	jwork = iwi + *n;
#line 306 "AB13ED.f"
	dgebal_("Scale", n, &dwork[ia2], n, &ilo, &ihi, &dwork[jwork], &i__, (
		ftnlen)5);
#line 308 "AB13ED.f"
	dhseqr_("Eigenvalues", "NoSchurVectors", n, &ilo, &ihi, &dwork[ia2], 
		n, &dwork[iwr], &dwork[iwi], dummy2, &c__1, &dwork[jwork], n, 
		info, (ftnlen)11, (ftnlen)14);

#line 312 "AB13ED.f"
	if (*info != 0) {
#line 313 "AB13ED.f"
	    *info = 1;
#line 314 "AB13ED.f"
	    return 0;
#line 315 "AB13ED.f"
	}

/*        (DWORK(IWR+i), DWORK(IWI+i)), i = 0,...,N-1, contain the */
/*        squares of the eigenvalues of H(sigma). */

#line 320 "AB13ED.f"
	i__ = 0;
#line 321 "AB13ED.f"
	rneg = FALSE_;
/*        WHILE ( ( DWORK(IWR+i),DWORK(IWI+i) ) not real positive */
/*                .AND. I < N ) DO */
#line 324 "AB13ED.f"
L40:
#line 324 "AB13ED.f"
	if (! rneg && i__ < *n) {
#line 325 "AB13ED.f"
	    temp = (d__1 = dwork[iwi + i__], abs(d__1));
#line 326 "AB13ED.f"
	    if (tol1 > sfmn) {
#line 326 "AB13ED.f"
		temp /= tol1;
#line 326 "AB13ED.f"
	    }
#line 327 "AB13ED.f"
	    rneg = dwork[iwr + i__] < 0. && temp <= tol2;
#line 328 "AB13ED.f"
	    ++i__;
#line 329 "AB13ED.f"
	    goto L40;
/*           END WHILE 40 */
#line 331 "AB13ED.f"
	}
#line 333 "AB13ED.f"
	if (rneg) {
#line 334 "AB13ED.f"
	    *high = sigma;
#line 335 "AB13ED.f"
	} else {
#line 336 "AB13ED.f"
	    *low = sigma;
#line 337 "AB13ED.f"
	}
#line 338 "AB13ED.f"
	goto L10;
/*        END WHILE 10 */
#line 340 "AB13ED.f"
    }

/*     Set optimal workspace dimension. */

/* Computing MAX */
#line 344 "AB13ED.f"
    i__2 = (n2 << 2) + *n;
#line 344 "AB13ED.f"
    dwork[1] = (doublereal) max(i__2,minwrk);

/* *** Last line of AB13ED *** */
#line 347 "AB13ED.f"
    return 0;
} /* ab13ed_ */

