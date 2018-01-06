#line 1 "SB08DD.f"
/* SB08DD.f -- translated by f2c (version 20100827).
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

#line 1 "SB08DD.f"
/* Table of constant values */

static doublereal c_b6 = 0.;
static doublereal c_b7 = 1.;
static integer c__1 = 1;
static integer c__4 = 4;
static logical c_true = TRUE_;

/* Subroutine */ int sb08dd_(char *dico, integer *n, integer *m, integer *p, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *d__, integer *ldd, integer *nq, 
	integer *nr, doublereal *cr, integer *ldcr, doublereal *dr, integer *
	lddr, doublereal *tol, doublereal *dwork, integer *ldwork, integer *
	iwarn, integer *info, ftnlen dico_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, cr_dim1, 
	    cr_offset, d_dim1, d_offset, dr_dim1, dr_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal x, y, z__[16]	/* was [4][4] */;
    static integer l1, ib, nb;
    static doublereal cs;
    static integer kv, kw;
    static doublereal pr, sm;
    static integer kz;
    static doublereal sn;
    static integer ib1, kfi, nfp, kwi, kwr, ncur;
    static doublereal rmax;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer nlow, nsup;
    static doublereal alpha;
    extern /* Subroutine */ int tb01ld_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, ftnlen, ftnlen, ftnlen), dgemm_(char *, 
	    char *, integer *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    extern /* Subroutine */ int sb01fy_(logical *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *);
    static doublereal bnorm;
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static doublereal toler;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen), dlanv2_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlaexc_(logical *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, doublereal *, integer *), dlacpy_(char *, integer *, integer *,
	     doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static doublereal wrkopt;


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

/*     To construct, for a given system G = (A,B,C,D), a feedback matrix */
/*     F, an orthogonal transformation matrix Z, and a gain matrix V, */
/*     such that the systems */

/*          Q = (Z'*(A+B*F)*Z, Z'*B*V, (C+D*F)*Z, D*V) */
/*     and */
/*          R = (Z'*(A+B*F)*Z, Z'*B*V, F*Z, V) */

/*     provide a stable right coprime factorization of G in the form */
/*                       -1 */
/*              G = Q * R  , */

/*     where G, Q and R are the corresponding transfer-function matrices */
/*     and the denominator R is inner, that is, R'(-s)*R(s) = I in the */
/*     continuous-time case, or R'(1/z)*R(z) = I in the discrete-time */
/*     case. The Z matrix is not explicitly computed. */

/*     Note: G must have no controllable poles on the imaginary axis */
/*     for a continuous-time system, or on the unit circle for a */
/*     discrete-time system. If the given state-space representation */
/*     is not stabilizable, the unstabilizable part of the original */
/*     system is automatically deflated and the order of the systems */
/*     Q and R is accordingly reduced. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the original system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The dimension of the state vector, i.e. the order of the */
/*             matrix A, and also the number of rows of the matrix B and */
/*             the number of columns of the matrices C and CR.  N >= 0. */

/*     M       (input) INTEGER */
/*             The dimension of input vector, i.e. the number of columns */
/*             of the matrices B, D and DR and the number of rows of the */
/*             matrices CR and DR.  M >= 0. */

/*     P       (input) INTEGER */
/*             The dimension of output vector, i.e. the number of rows */
/*             of the matrices C and D.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix A. The matrix A must not */
/*             have controllable eigenvalues on the imaginary axis, if */
/*             DICO = 'C', or on the unit circle, if DICO = 'D'. */
/*             On exit, the leading NQ-by-NQ part of this array contains */
/*             the leading NQ-by-NQ part of the matrix Z'*(A+B*F)*Z, the */
/*             state dynamics matrix of the numerator factor Q, in a */
/*             real Schur form. The trailing NR-by-NR part of this matrix */
/*             represents the state dynamics matrix of a minimal */
/*             realization of the denominator factor R. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input/state matrix. */
/*             On exit, the leading NQ-by-M part of this array contains */
/*             the leading NQ-by-M part of the matrix Z'*B*V, the */
/*             input/state matrix of the numerator factor Q. The last */
/*             NR rows of this matrix form the input/state matrix of */
/*             a minimal realization of the denominator factor R. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C. */
/*             On exit, the leading P-by-NQ part of this array contains */
/*             the leading P-by-NQ part of the matrix (C+D*F)*Z, */
/*             the state/output matrix of the numerator factor Q. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the input/output matrix. */
/*             On exit, the leading P-by-M part of this array contains */
/*             the matrix D*V representing the input/output matrix */
/*             of the numerator factor Q. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P). */

/*     NQ      (output) INTEGER */
/*             The order of the resulting factors Q and R. */
/*             Generally, NQ = N - NS, where NS is the number of */
/*             uncontrollable eigenvalues outside the stability region. */

/*     NR      (output) INTEGER */
/*             The order of the minimal realization of the factor R. */
/*             Generally, NR is the number of controllable eigenvalues */
/*             of A outside the stability region (the number of modified */
/*             eigenvalues). */

/*     CR      (output) DOUBLE PRECISION array, dimension (LDCR,N) */
/*             The leading M-by-NQ part of this array contains the */
/*             leading M-by-NQ part of the feedback matrix F*Z, which */
/*             reflects the eigenvalues of A lying outside the stable */
/*             region to values which are symmetric with respect to the */
/*             imaginary axis (if DICO = 'C') or the unit circle (if */
/*             DICO = 'D').  The last NR columns of this matrix form the */
/*             state/output matrix of a minimal realization of the */
/*             denominator factor R. */

/*     LDCR    INTEGER */
/*             The leading dimension of array CR.  LDCR >= MAX(1,M). */

/*     DR      (output) DOUBLE PRECISION array, dimension (LDDR,M) */
/*             The leading M-by-M part of this array contains the upper */
/*             triangular matrix V of order M representing the */
/*             input/output matrix of the denominator factor R. */

/*     LDDR    INTEGER */
/*             The leading dimension of array DR.  LDDR >= MAX(1,M). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The absolute tolerance level below which the elements of */
/*             B are considered zero (used for controllability tests). */
/*             If the user sets TOL <= 0, then an implicitly computed, */
/*             default tolerance, defined by  TOLDEF = N*EPS*NORM(B), */
/*             is used instead, where EPS is the machine precision */
/*             (see LAPACK Library routine DLAMCH) and NORM(B) denotes */
/*             the 1-norm of B. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of working array DWORK. */
/*             LDWORK >= MAX( 1, N*(N+5), M*(M+2), 4*M, 4*P ). */
/*             For optimum performance LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = K:  K violations of the numerical stability condition */
/*                   NORM(F) <= 10*NORM(A)/NORM(B) occured during the */
/*                   assignment of eigenvalues. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the reduction of A to a real Schur form failed; */
/*             = 2:  a failure was detected during the ordering of the */
/*                   real Schur form of A, or in the iterative process */
/*                   for reordering the eigenvalues of Z'*(A + B*F)*Z */
/*                   along the diagonal; */
/*             = 3:  if DICO = 'C' and the matrix A has a controllable */
/*                   eigenvalue on the imaginary axis, or DICO = 'D' */
/*                   and A has a controllable eigenvalue on the unit */
/*                   circle. */

/*     METHOD */

/*     The subroutine is based on the factorization algorithm of [1]. */

/*     REFERENCES */

/*     [1] Varga A. */
/*         A Schur method for computing coprime factorizations with inner */
/*         denominators and applications in model reduction. */
/*         Proc. ACC'93, San Francisco, CA, pp. 2130-2131, 1993. */

/*     NUMERICAL ASPECTS */
/*                                            3 */
/*     The algorithm requires no more than 14N  floating point */
/*     operations. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, July 1998. */
/*     Based on the RASP routine RCFID. */

/*     REVISIONS */

/*     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven. */
/*     Feb. 1999, May 2003, A. Varga, DLR Oberpfaffenhofen. */

/*     KEYWORDS */

/*     Coprime factorization, eigenvalue, eigenvalue assignment, */
/*     feedback control, pole placement, state-space model. */

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

#line 264 "SB08DD.f"
    /* Parameter adjustments */
#line 264 "SB08DD.f"
    a_dim1 = *lda;
#line 264 "SB08DD.f"
    a_offset = 1 + a_dim1;
#line 264 "SB08DD.f"
    a -= a_offset;
#line 264 "SB08DD.f"
    b_dim1 = *ldb;
#line 264 "SB08DD.f"
    b_offset = 1 + b_dim1;
#line 264 "SB08DD.f"
    b -= b_offset;
#line 264 "SB08DD.f"
    c_dim1 = *ldc;
#line 264 "SB08DD.f"
    c_offset = 1 + c_dim1;
#line 264 "SB08DD.f"
    c__ -= c_offset;
#line 264 "SB08DD.f"
    d_dim1 = *ldd;
#line 264 "SB08DD.f"
    d_offset = 1 + d_dim1;
#line 264 "SB08DD.f"
    d__ -= d_offset;
#line 264 "SB08DD.f"
    cr_dim1 = *ldcr;
#line 264 "SB08DD.f"
    cr_offset = 1 + cr_dim1;
#line 264 "SB08DD.f"
    cr -= cr_offset;
#line 264 "SB08DD.f"
    dr_dim1 = *lddr;
#line 264 "SB08DD.f"
    dr_offset = 1 + dr_dim1;
#line 264 "SB08DD.f"
    dr -= dr_offset;
#line 264 "SB08DD.f"
    --dwork;
#line 264 "SB08DD.f"

#line 264 "SB08DD.f"
    /* Function Body */
#line 264 "SB08DD.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 265 "SB08DD.f"
    *iwarn = 0;
#line 266 "SB08DD.f"
    *info = 0;

/*     Check the scalar input parameters. */

#line 270 "SB08DD.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 271 "SB08DD.f"
	*info = -1;
#line 272 "SB08DD.f"
    } else if (*n < 0) {
#line 273 "SB08DD.f"
	*info = -2;
#line 274 "SB08DD.f"
    } else if (*m < 0) {
#line 275 "SB08DD.f"
	*info = -3;
#line 276 "SB08DD.f"
    } else if (*p < 0) {
#line 277 "SB08DD.f"
	*info = -4;
#line 278 "SB08DD.f"
    } else if (*lda < max(1,*n)) {
#line 279 "SB08DD.f"
	*info = -6;
#line 280 "SB08DD.f"
    } else if (*ldb < max(1,*n)) {
#line 281 "SB08DD.f"
	*info = -8;
#line 282 "SB08DD.f"
    } else if (*ldc < max(1,*p)) {
#line 283 "SB08DD.f"
	*info = -10;
#line 284 "SB08DD.f"
    } else if (*ldd < max(1,*p)) {
#line 285 "SB08DD.f"
	*info = -12;
#line 286 "SB08DD.f"
    } else if (*ldcr < max(1,*m)) {
#line 287 "SB08DD.f"
	*info = -16;
#line 288 "SB08DD.f"
    } else if (*lddr < max(1,*m)) {
#line 289 "SB08DD.f"
	*info = -18;
#line 290 "SB08DD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 290 "SB08DD.f"
	i__1 = 1, i__2 = *n * (*n + 5), i__1 = max(i__1,i__2), i__2 = *m * (*
		m + 2), i__1 = max(i__1,i__2), i__2 = *m << 2, i__1 = max(
		i__1,i__2), i__2 = *p << 2;
#line 290 "SB08DD.f"
	if (*ldwork < max(i__1,i__2)) {
#line 291 "SB08DD.f"
	    *info = -21;
#line 292 "SB08DD.f"
	}
#line 292 "SB08DD.f"
    }
#line 293 "SB08DD.f"
    if (*info != 0) {

/*        Error return. */

#line 297 "SB08DD.f"
	i__1 = -(*info);
#line 297 "SB08DD.f"
	xerbla_("SB08DD", &i__1, (ftnlen)6);
#line 298 "SB08DD.f"
	return 0;
#line 299 "SB08DD.f"
    }

/*     Set DR = I and quick return if possible. */

#line 303 "SB08DD.f"
    *nr = 0;
#line 304 "SB08DD.f"
    if (min(*m,*p) > 0) {
#line 304 "SB08DD.f"
	dlaset_("Full", m, m, &c_b6, &c_b7, &dr[dr_offset], lddr, (ftnlen)4);
#line 304 "SB08DD.f"
    }
#line 306 "SB08DD.f"
    if (min(*n,*m) == 0) {
#line 307 "SB08DD.f"
	*nq = 0;
#line 308 "SB08DD.f"
	dwork[1] = 1.;
#line 309 "SB08DD.f"
	return 0;
#line 310 "SB08DD.f"
    }

/*     Set F = 0 in the array CR. */

#line 314 "SB08DD.f"
    dlaset_("Full", m, n, &c_b6, &c_b6, &cr[cr_offset], ldcr, (ftnlen)4);

/*     Compute the norm of B and set the default tolerance if necessary. */

#line 318 "SB08DD.f"
    bnorm = dlange_("1-norm", n, m, &b[b_offset], ldb, &dwork[1], (ftnlen)6);
#line 319 "SB08DD.f"
    toler = *tol;
#line 320 "SB08DD.f"
    if (toler <= 0.) {
#line 320 "SB08DD.f"
	toler = (doublereal) (*n) * bnorm * dlamch_("Epsilon", (ftnlen)7);
#line 320 "SB08DD.f"
    }
#line 322 "SB08DD.f"
    if (bnorm <= toler) {
#line 323 "SB08DD.f"
	*nq = 0;
#line 324 "SB08DD.f"
	dwork[1] = 1.;
#line 325 "SB08DD.f"
	return 0;
#line 326 "SB08DD.f"
    }

/*     Compute the bound for the numerical stability condition. */

#line 330 "SB08DD.f"
    rmax = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)6) * 
	    10. / bnorm;

/*     Allocate working storage. */

#line 334 "SB08DD.f"
    kz = 1;
#line 335 "SB08DD.f"
    kwr = kz + *n * *n;
#line 336 "SB08DD.f"
    kwi = kwr + *n;
#line 337 "SB08DD.f"
    kw = kwi + *n;

/*     Reduce A to an ordered real Schur form using an orthogonal */
/*     similarity transformation A <- Z'*A*Z and accumulate the */
/*     transformations in Z.  The separation of spectrum of A is */
/*     performed such that the leading NFP-by-NFP submatrix of A */
/*     corresponds to the "stable" eigenvalues which will be not */
/*     modified. The bottom (N-NFP)-by-(N-NFP) diagonal block of A */
/*     corresponds to the "unstable" eigenvalues to be modified. */
/*     Apply the transformation to B and C: B <- Z'*B and C <- C*Z. */

/*     Workspace needed:      N*(N+2); */
/*     Additional workspace:  need   3*N; */
/*                            prefer larger. */

#line 352 "SB08DD.f"
    if (discr) {
#line 353 "SB08DD.f"
	alpha = 1.;
#line 354 "SB08DD.f"
    } else {
#line 355 "SB08DD.f"
	alpha = 0.;
#line 356 "SB08DD.f"
    }
#line 357 "SB08DD.f"
    i__1 = *ldwork - kw + 1;
#line 357 "SB08DD.f"
    tb01ld_(dico, "Stable", "General", n, m, p, &alpha, &a[a_offset], lda, &b[
	    b_offset], ldb, &c__[c_offset], ldc, &nfp, &dwork[kz], n, &dwork[
	    kwr], &dwork[kwi], &dwork[kw], &i__1, info, (ftnlen)1, (ftnlen)6, 
	    (ftnlen)7);
#line 360 "SB08DD.f"
    if (*info != 0) {
#line 360 "SB08DD.f"
	return 0;
#line 360 "SB08DD.f"
    }

#line 363 "SB08DD.f"
    wrkopt = dwork[kw] + (doublereal) (kw - 1);

/*     Perform the pole assignment if there exist "unstable" eigenvalues. */

#line 367 "SB08DD.f"
    *nq = *n;
#line 368 "SB08DD.f"
    if (nfp < *n) {
#line 369 "SB08DD.f"
	kv = 1;
#line 370 "SB08DD.f"
	kfi = kv + *m * *m;
#line 371 "SB08DD.f"
	kw = kfi + (*m << 1);

/*        Set the limits for the bottom diagonal block. */

#line 375 "SB08DD.f"
	nlow = nfp + 1;
#line 376 "SB08DD.f"
	nsup = *n;

/*        WHILE (NLOW <= NSUP) DO */
#line 379 "SB08DD.f"
L10:
#line 379 "SB08DD.f"
	if (nlow <= nsup) {

/*           Main loop for assigning one or two poles. */

/*           Determine the dimension of the last block. */

#line 385 "SB08DD.f"
	    ib = 1;
#line 386 "SB08DD.f"
	    if (nlow < nsup) {
#line 387 "SB08DD.f"
		if (a[nsup + (nsup - 1) * a_dim1] != 0.) {
#line 387 "SB08DD.f"
		    ib = 2;
#line 387 "SB08DD.f"
		}
#line 388 "SB08DD.f"
	    }
#line 389 "SB08DD.f"
	    l = nsup - ib + 1;

/*           Check the controllability of the last block. */

#line 393 "SB08DD.f"
	    if (dlange_("1-norm", &ib, m, &b[l + b_dim1], ldb, &dwork[kw], (
		    ftnlen)6) <= toler) {

/*              Deflate the uncontrollable block and resume the main */
/*              loop. */

#line 399 "SB08DD.f"
		nsup -= ib;
#line 400 "SB08DD.f"
	    } else {

/*              Determine the M-by-IB feedback matrix FI which assigns */
/*              the selected IB poles for the pair */
/*              ( A(L:L+IB-1,L:L+IB-1), B(L:L+IB-1,1:M) ). */

/*              Workspace needed: M*(M+2). */

#line 408 "SB08DD.f"
		sb01fy_(&discr, &ib, m, &a[l + l * a_dim1], lda, &b[l + 
			b_dim1], ldb, &dwork[kfi], m, &dwork[kv], m, info);
#line 410 "SB08DD.f"
		if (*info == 2) {
#line 411 "SB08DD.f"
		    *info = 3;
#line 412 "SB08DD.f"
		    return 0;
#line 413 "SB08DD.f"
		}

/*              Check for possible numerical instability. */

#line 417 "SB08DD.f"
		if (dlange_("1-norm", m, &ib, &dwork[kfi], m, &dwork[kw], (
			ftnlen)6) > rmax) {
#line 417 "SB08DD.f"
		    ++(*iwarn);
#line 417 "SB08DD.f"
		}

/*              Update the state matrix A <-- A + B*[0 FI]. */

#line 422 "SB08DD.f"
		dgemm_("NoTranspose", "NoTranspose", &nsup, &ib, m, &c_b7, &b[
			b_offset], ldb, &dwork[kfi], m, &c_b7, &a[l * a_dim1 
			+ 1], lda, (ftnlen)11, (ftnlen)11);

/*              Update the feedback matrix F <-- F + V*[0 FI] in CR. */

#line 428 "SB08DD.f"
		if (discr) {
#line 428 "SB08DD.f"
		    dtrmm_("Left", "Upper", "NoTranspose", "NonUnit", m, &ib, 
			    &c_b7, &dr[dr_offset], lddr, &dwork[kfi], m, (
			    ftnlen)4, (ftnlen)5, (ftnlen)11, (ftnlen)7);
#line 428 "SB08DD.f"
		}
#line 431 "SB08DD.f"
		k = kfi;
#line 432 "SB08DD.f"
		i__1 = l + ib - 1;
#line 432 "SB08DD.f"
		for (j = l; j <= i__1; ++j) {
#line 433 "SB08DD.f"
		    i__2 = *m;
#line 433 "SB08DD.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 434 "SB08DD.f"
			cr[i__ + j * cr_dim1] += dwork[k];
#line 435 "SB08DD.f"
			++k;
#line 436 "SB08DD.f"
/* L20: */
#line 436 "SB08DD.f"
		    }
#line 437 "SB08DD.f"
/* L30: */
#line 437 "SB08DD.f"
		}

#line 439 "SB08DD.f"
		if (discr) {

/*                 Update the input matrix B <-- B*V. */

#line 443 "SB08DD.f"
		    dtrmm_("Right", "Upper", "NoTranspose", "NonUnit", n, m, &
			    c_b7, &dwork[kv], m, &b[b_offset], ldb, (ftnlen)5,
			     (ftnlen)5, (ftnlen)11, (ftnlen)7);

/*                 Update the feedthrough matrix DR <-- DR*V. */

#line 449 "SB08DD.f"
		    k = kv;
#line 450 "SB08DD.f"
		    i__1 = *m;
#line 450 "SB08DD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 451 "SB08DD.f"
			i__2 = *m - i__ + 1;
#line 451 "SB08DD.f"
			dtrmv_("Upper", "Transpose", "NonUnit", &i__2, &dwork[
				k], m, &dr[i__ + i__ * dr_dim1], lddr, (
				ftnlen)5, (ftnlen)9, (ftnlen)7);
#line 453 "SB08DD.f"
			k = k + *m + 1;
#line 454 "SB08DD.f"
/* L40: */
#line 454 "SB08DD.f"
		    }
#line 455 "SB08DD.f"
		}

#line 457 "SB08DD.f"
		if (ib == 2) {

/*                 Put the 2x2 block in a standard form. */

#line 461 "SB08DD.f"
		    l1 = l + 1;
#line 462 "SB08DD.f"
		    dlanv2_(&a[l + l * a_dim1], &a[l + l1 * a_dim1], &a[l1 + 
			    l * a_dim1], &a[l1 + l1 * a_dim1], &x, &y, &pr, &
			    sm, &cs, &sn);

/*                 Apply the transformation to A, B, C and F. */

#line 467 "SB08DD.f"
		    if (l1 < nsup) {
#line 467 "SB08DD.f"
			i__1 = nsup - l1;
#line 467 "SB08DD.f"
			drot_(&i__1, &a[l + (l1 + 1) * a_dim1], lda, &a[l1 + (
				l1 + 1) * a_dim1], lda, &cs, &sn);
#line 467 "SB08DD.f"
		    }
#line 470 "SB08DD.f"
		    i__1 = l - 1;
#line 470 "SB08DD.f"
		    drot_(&i__1, &a[l * a_dim1 + 1], &c__1, &a[l1 * a_dim1 + 
			    1], &c__1, &cs, &sn);
#line 471 "SB08DD.f"
		    drot_(m, &b[l + b_dim1], ldb, &b[l1 + b_dim1], ldb, &cs, &
			    sn);
#line 472 "SB08DD.f"
		    if (*p > 0) {
#line 472 "SB08DD.f"
			drot_(p, &c__[l * c_dim1 + 1], &c__1, &c__[l1 * 
				c_dim1 + 1], &c__1, &cs, &sn);
#line 472 "SB08DD.f"
		    }
#line 474 "SB08DD.f"
		    drot_(m, &cr[l * cr_dim1 + 1], &c__1, &cr[l1 * cr_dim1 + 
			    1], &c__1, &cs, &sn);
#line 475 "SB08DD.f"
		}
#line 476 "SB08DD.f"
		if (nlow + ib <= nsup) {

/*                 Move the last block(s) to the leading position(s) of */
/*                 the bottom block. */

/*                 Workspace:     need MAX(4*N, 4*M, 4*P). */

#line 483 "SB08DD.f"
		    ncur = nsup - ib;
/*                 WHILE (NCUR >= NLOW) DO */
#line 485 "SB08DD.f"
L50:
#line 485 "SB08DD.f"
		    if (ncur >= nlow) {

/*                    Loop for positioning of the last block. */

/*                    Determine the dimension of the current block. */

#line 491 "SB08DD.f"
			ib1 = 1;
#line 492 "SB08DD.f"
			if (ncur > nlow) {
#line 493 "SB08DD.f"
			    if (a[ncur + (ncur - 1) * a_dim1] != 0.) {
#line 493 "SB08DD.f"
				ib1 = 2;
#line 493 "SB08DD.f"
			    }
#line 494 "SB08DD.f"
			}
#line 495 "SB08DD.f"
			nb = ib1 + ib;

/*                    Initialize the local transformation matrix Z. */

#line 499 "SB08DD.f"
			dlaset_("Full", &nb, &nb, &c_b6, &c_b7, z__, &c__4, (
				ftnlen)4);
#line 500 "SB08DD.f"
			l = ncur - ib1 + 1;

/*                    Exchange two adjacent blocks and accumulate the */
/*                    transformations in Z. */

#line 505 "SB08DD.f"
			dlaexc_(&c_true, &nb, &a[l + l * a_dim1], lda, z__, &
				c__4, &c__1, &ib1, &ib, &dwork[1], info);
#line 507 "SB08DD.f"
			if (*info != 0) {
#line 508 "SB08DD.f"
			    *info = 2;
#line 509 "SB08DD.f"
			    return 0;
#line 510 "SB08DD.f"
			}

/*                    Apply the transformation to the rest of A. */

#line 514 "SB08DD.f"
			l1 = l + nb;
#line 515 "SB08DD.f"
			if (l1 <= nsup) {
#line 516 "SB08DD.f"
			    i__1 = nsup - l1 + 1;
#line 516 "SB08DD.f"
			    dgemm_("Transpose", "NoTranspose", &nb, &i__1, &
				    nb, &c_b7, z__, &c__4, &a[l + l1 * a_dim1]
				    , lda, &c_b6, &dwork[1], &nb, (ftnlen)9, (
				    ftnlen)11);
#line 519 "SB08DD.f"
			    i__1 = nsup - l1 + 1;
#line 519 "SB08DD.f"
			    dlacpy_("Full", &nb, &i__1, &dwork[1], &nb, &a[l 
				    + l1 * a_dim1], lda, (ftnlen)4);
#line 521 "SB08DD.f"
			}
#line 522 "SB08DD.f"
			i__1 = l - 1;
#line 522 "SB08DD.f"
			dgemm_("NoTranspose", "NoTranspose", &i__1, &nb, &nb, 
				&c_b7, &a[l * a_dim1 + 1], lda, z__, &c__4, &
				c_b6, &dwork[1], n, (ftnlen)11, (ftnlen)11);
#line 525 "SB08DD.f"
			i__1 = l - 1;
#line 525 "SB08DD.f"
			dlacpy_("Full", &i__1, &nb, &dwork[1], n, &a[l * 
				a_dim1 + 1], lda, (ftnlen)4);

/*                    Apply the transformation to B, C and F. */

#line 530 "SB08DD.f"
			dgemm_("Transpose", "NoTranspose", &nb, m, &nb, &c_b7,
				 z__, &c__4, &b[l + b_dim1], ldb, &c_b6, &
				dwork[1], &nb, (ftnlen)9, (ftnlen)11);
#line 533 "SB08DD.f"
			dlacpy_("Full", &nb, m, &dwork[1], &nb, &b[l + b_dim1]
				, ldb, (ftnlen)4);

#line 536 "SB08DD.f"
			if (*p > 0) {
#line 537 "SB08DD.f"
			    dgemm_("NoTranspose", "NoTranspose", p, &nb, &nb, 
				    &c_b7, &c__[l * c_dim1 + 1], ldc, z__, &
				    c__4, &c_b6, &dwork[1], p, (ftnlen)11, (
				    ftnlen)11);
#line 540 "SB08DD.f"
			    dlacpy_("Full", p, &nb, &dwork[1], p, &c__[l * 
				    c_dim1 + 1], ldc, (ftnlen)4);
#line 542 "SB08DD.f"
			}

#line 544 "SB08DD.f"
			dgemm_("NoTranspose", "NoTranspose", m, &nb, &nb, &
				c_b7, &cr[l * cr_dim1 + 1], ldcr, z__, &c__4, 
				&c_b6, &dwork[1], m, (ftnlen)11, (ftnlen)11);
#line 547 "SB08DD.f"
			dlacpy_("Full", m, &nb, &dwork[1], m, &cr[l * cr_dim1 
				+ 1], ldcr, (ftnlen)4);

#line 550 "SB08DD.f"
			ncur -= ib1;
#line 551 "SB08DD.f"
			goto L50;
#line 552 "SB08DD.f"
		    }
/*                 END WHILE 50 */

#line 555 "SB08DD.f"
		}
#line 556 "SB08DD.f"
		nlow += ib;
#line 557 "SB08DD.f"
	    }
#line 558 "SB08DD.f"
	    goto L10;
#line 559 "SB08DD.f"
	}
/*        END WHILE 10 */

#line 562 "SB08DD.f"
	*nq = nsup;
#line 563 "SB08DD.f"
	*nr = nsup - nfp;

/*        Annihilate the elements below the first subdiagonal of A. */

#line 567 "SB08DD.f"
	if (*nq > 2) {
#line 567 "SB08DD.f"
	    i__1 = *nq - 2;
#line 567 "SB08DD.f"
	    i__2 = *nq - 2;
#line 567 "SB08DD.f"
	    dlaset_("Lower", &i__1, &i__2, &c_b6, &c_b6, &a[a_dim1 + 3], lda, 
		    (ftnlen)5);
#line 567 "SB08DD.f"
	}
#line 569 "SB08DD.f"
    }

/*     Compute C <-- CQ = C + D*F and D <-- DQ = D*DR. */

#line 573 "SB08DD.f"
    dgemm_("NoTranspose", "NoTranspose", p, nq, m, &c_b7, &d__[d_offset], ldd,
	     &cr[cr_offset], ldcr, &c_b7, &c__[c_offset], ldc, (ftnlen)11, (
	    ftnlen)11);
#line 575 "SB08DD.f"
    if (discr) {
#line 575 "SB08DD.f"
	dtrmm_("Right", "Upper", "NoTranspose", "NonUnit", p, m, &c_b7, &dr[
		dr_offset], lddr, &d__[d_offset], ldd, (ftnlen)5, (ftnlen)5, (
		ftnlen)11, (ftnlen)7);
#line 575 "SB08DD.f"
    }

/* Computing MAX */
/* Computing MAX */
#line 579 "SB08DD.f"
    i__1 = *m * (*m + 2), i__2 = *m << 2, i__1 = max(i__1,i__2), i__2 = *p << 
	    2;
#line 579 "SB08DD.f"
    d__1 = wrkopt, d__2 = (doublereal) max(i__1,i__2);
#line 579 "SB08DD.f"
    dwork[1] = max(d__1,d__2);

#line 581 "SB08DD.f"
    return 0;
/* *** Last line of SB08DD *** */
} /* sb08dd_ */

