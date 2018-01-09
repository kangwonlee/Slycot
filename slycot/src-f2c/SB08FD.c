#line 1 "SB08FD.f"
/* SB08FD.f -- translated by f2c (version 20100827).
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

#line 1 "SB08FD.f"
/* Table of constant values */

static doublereal c_b6 = 0.;
static doublereal c_b7 = 1.;
static integer c__1 = 1;
static integer c__4 = 4;
static logical c_true = TRUE_;

/* Subroutine */ int sb08fd_(char *dico, integer *n, integer *m, integer *p, 
	doublereal *alpha, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *d__, integer 
	*ldd, integer *nq, integer *nr, doublereal *cr, integer *ldcr, 
	doublereal *dr, integer *lddr, doublereal *tol, doublereal *dwork, 
	integer *ldwork, integer *iwarn, integer *info, ftnlen dico_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, cr_dim1, 
	    cr_offset, d_dim1, d_offset, dr_dim1, dr_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal x, y, z__[16]	/* was [4][4] */, a2[4]	/* was [2][2] 
	    */;
    static integer l1, ib, nb, kg;
    static doublereal cs, sm;
    static integer kw;
    static doublereal pr, sn;
    static integer kz, ib1, kfi, nfp, kwi, kwr, ncur;
    static doublereal rmax;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer nlow, nsup, ncur1;
    extern /* Subroutine */ int tb01ld_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, ftnlen, ftnlen, ftnlen), dgemm_(char *, 
	    char *, integer *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen), sb01by_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    static doublereal bnorm, toler;
    extern /* Subroutine */ int dlanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen), dlange_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlaexc_(logical *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, doublereal *, integer *), dlacpy_(char *, integer *, integer *,
	     doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static integer nmoves;
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

/*     To construct, for a given system G = (A,B,C,D), a feedback */
/*     matrix F and an orthogonal transformation matrix Z, such that */
/*     the systems */

/*          Q = (Z'*(A+B*F)*Z, Z'*B, (C+D*F)*Z, D) */
/*     and */
/*          R = (Z'*(A+B*F)*Z, Z'*B, F*Z, I) */

/*     provide a stable right coprime factorization of G in the form */
/*                       -1 */
/*              G = Q * R  , */

/*     where G, Q and R are the corresponding transfer-function matrices. */
/*     The resulting state dynamics matrix of the systems Q and R has */
/*     eigenvalues lying inside a given stability domain. */
/*     The Z matrix is not explicitly computed. */

/*     Note: If the given state-space representation is not stabilizable, */
/*     the unstabilizable part of the original system is automatically */
/*     deflated and the order of the systems Q and R is accordingly */
/*     reduced. */

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

/*     ALPHA   (input) DOUBLE PRECISION array, dimension (2) */
/*             ALPHA(1) contains the desired stability degree to be */
/*             assigned for the eigenvalues of A+B*F, and ALPHA(2) */
/*             the stability margin. The eigenvalues outside the */
/*             ALPHA(2)-stability region will be assigned to have the */
/*             real parts equal to ALPHA(1) < 0 and unmodified */
/*             imaginary parts for a continuous-time system */
/*             (DICO = 'C'), or moduli equal to 0 <= ALPHA(2) < 1 */
/*             for a discrete-time system (DICO = 'D'). */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix A. */
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
/*             the leading NQ-by-M part of the matrix Z'*B, the */
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

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading P-by-M part of this array must contain the */
/*             input/output matrix. D represents also the input/output */
/*             matrix of the numerator factor Q. */

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
/*             moves the eigenvalues of A lying outside the ALPHA-stable */
/*             region to values which are on the ALPHA-stability */
/*             boundary.  The last NR columns of this matrix form the */
/*             state/output matrix of a minimal realization of the */
/*             denominator factor R. */

/*     LDCR    INTEGER */
/*             The leading dimension of array CR.  LDCR >= MAX(1,M). */

/*     DR      (output) DOUBLE PRECISION array, dimension (LDDR,M) */
/*             The leading M-by-M part of this array contains an */
/*             identity matrix representing the input/output matrix */
/*             of the denominator factor R. */

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
/*             LWORK >= MAX( 1, N*(N+5), 5*M, 4*P ). */
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
/*                   along the diagonal. */

/*     METHOD */

/*     The subroutine is based on the factorization algorithm of [1]. */

/*     REFERENCES */

/*     [1] Varga A. */
/*         Coprime factors model reduction method based on */
/*         square-root balancing-free techniques. */
/*         System Analysis, Modelling and Simulation, */
/*         vol. 11, pp. 303-311, 1993. */

/*     NUMERICAL ASPECTS */
/*                                            3 */
/*     The algorithm requires no more than 14N  floating point */
/*     operations. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, July 1998. */
/*     Based on the RASP routine RCFS. */

/*     REVISIONS */

/*     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven. */
/*     Mar. 2003, May 2003, A. Varga, German Aerospace Center. */
/*     May 2003, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Sep. 2005, A. Varga, German Aerospace Center. */

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

#line 266 "SB08FD.f"
    /* Parameter adjustments */
#line 266 "SB08FD.f"
    --alpha;
#line 266 "SB08FD.f"
    a_dim1 = *lda;
#line 266 "SB08FD.f"
    a_offset = 1 + a_dim1;
#line 266 "SB08FD.f"
    a -= a_offset;
#line 266 "SB08FD.f"
    b_dim1 = *ldb;
#line 266 "SB08FD.f"
    b_offset = 1 + b_dim1;
#line 266 "SB08FD.f"
    b -= b_offset;
#line 266 "SB08FD.f"
    c_dim1 = *ldc;
#line 266 "SB08FD.f"
    c_offset = 1 + c_dim1;
#line 266 "SB08FD.f"
    c__ -= c_offset;
#line 266 "SB08FD.f"
    d_dim1 = *ldd;
#line 266 "SB08FD.f"
    d_offset = 1 + d_dim1;
#line 266 "SB08FD.f"
    d__ -= d_offset;
#line 266 "SB08FD.f"
    cr_dim1 = *ldcr;
#line 266 "SB08FD.f"
    cr_offset = 1 + cr_dim1;
#line 266 "SB08FD.f"
    cr -= cr_offset;
#line 266 "SB08FD.f"
    dr_dim1 = *lddr;
#line 266 "SB08FD.f"
    dr_offset = 1 + dr_dim1;
#line 266 "SB08FD.f"
    dr -= dr_offset;
#line 266 "SB08FD.f"
    --dwork;
#line 266 "SB08FD.f"

#line 266 "SB08FD.f"
    /* Function Body */
#line 266 "SB08FD.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 267 "SB08FD.f"
    *iwarn = 0;
#line 268 "SB08FD.f"
    *info = 0;

/*     Check the scalar input parameters. */

#line 272 "SB08FD.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 273 "SB08FD.f"
	*info = -1;
#line 274 "SB08FD.f"
    } else if (*n < 0) {
#line 275 "SB08FD.f"
	*info = -2;
#line 276 "SB08FD.f"
    } else if (*m < 0) {
#line 277 "SB08FD.f"
	*info = -3;
#line 278 "SB08FD.f"
    } else if (*p < 0) {
#line 279 "SB08FD.f"
	*info = -4;
#line 280 "SB08FD.f"
    } else if (discr && (alpha[1] < 0. || alpha[1] >= 1. || alpha[2] < 0. || 
	    alpha[2] >= 1.) || ! discr && (alpha[1] >= 0. || alpha[2] >= 0.)) 
	    {
#line 285 "SB08FD.f"
	*info = -5;
#line 286 "SB08FD.f"
    } else if (*lda < max(1,*n)) {
#line 287 "SB08FD.f"
	*info = -7;
#line 288 "SB08FD.f"
    } else if (*ldb < max(1,*n)) {
#line 289 "SB08FD.f"
	*info = -9;
#line 290 "SB08FD.f"
    } else if (*ldc < max(1,*p)) {
#line 291 "SB08FD.f"
	*info = -11;
#line 292 "SB08FD.f"
    } else if (*ldd < max(1,*p)) {
#line 293 "SB08FD.f"
	*info = -13;
#line 294 "SB08FD.f"
    } else if (*ldcr < max(1,*m)) {
#line 295 "SB08FD.f"
	*info = -17;
#line 296 "SB08FD.f"
    } else if (*lddr < max(1,*m)) {
#line 297 "SB08FD.f"
	*info = -19;
#line 298 "SB08FD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 298 "SB08FD.f"
	i__1 = 1, i__2 = *n * (*n + 5), i__1 = max(i__1,i__2), i__2 = *m * 5, 
		i__1 = max(i__1,i__2), i__2 = *p << 2;
#line 298 "SB08FD.f"
	if (*ldwork < max(i__1,i__2)) {
#line 299 "SB08FD.f"
	    *info = -22;
#line 300 "SB08FD.f"
	}
#line 300 "SB08FD.f"
    }
#line 301 "SB08FD.f"
    if (*info != 0) {

/*        Error return. */

#line 305 "SB08FD.f"
	i__1 = -(*info);
#line 305 "SB08FD.f"
	xerbla_("SB08FD", &i__1, (ftnlen)6);
#line 306 "SB08FD.f"
	return 0;
#line 307 "SB08FD.f"
    }

/*     Set DR = I and quick return if possible. */

#line 311 "SB08FD.f"
    *nr = 0;
#line 312 "SB08FD.f"
    dlaset_("Full", m, m, &c_b6, &c_b7, &dr[dr_offset], lddr, (ftnlen)4);
#line 313 "SB08FD.f"
    if (min(*n,*m) == 0) {
#line 314 "SB08FD.f"
	*nq = 0;
#line 315 "SB08FD.f"
	dwork[1] = 1.;
#line 316 "SB08FD.f"
	return 0;
#line 317 "SB08FD.f"
    }

/*     Set F = 0 in the array CR. */

#line 321 "SB08FD.f"
    dlaset_("Full", m, n, &c_b6, &c_b6, &cr[cr_offset], ldcr, (ftnlen)4);

/*     Compute the norm of B and set the default tolerance if necessary. */

#line 325 "SB08FD.f"
    bnorm = dlange_("1-norm", n, m, &b[b_offset], ldb, &dwork[1], (ftnlen)6);
#line 326 "SB08FD.f"
    toler = *tol;
#line 327 "SB08FD.f"
    if (toler <= 0.) {
#line 327 "SB08FD.f"
	toler = (doublereal) (*n) * bnorm * dlamch_("Epsilon", (ftnlen)7);
#line 327 "SB08FD.f"
    }
#line 329 "SB08FD.f"
    if (bnorm <= toler) {
#line 330 "SB08FD.f"
	*nq = 0;
#line 331 "SB08FD.f"
	dwork[1] = 1.;
#line 332 "SB08FD.f"
	return 0;
#line 333 "SB08FD.f"
    }

/*     Compute the bound for the numerical stability condition. */

#line 337 "SB08FD.f"
    rmax = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)6) * 
	    10. / bnorm;

/*     Allocate working storage. */

#line 341 "SB08FD.f"
    kz = 1;
#line 342 "SB08FD.f"
    kwr = kz + *n * *n;
#line 343 "SB08FD.f"
    kwi = kwr + *n;
#line 344 "SB08FD.f"
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

#line 359 "SB08FD.f"
    i__1 = *ldwork - kw + 1;
#line 359 "SB08FD.f"
    tb01ld_(dico, "Stable", "General", n, m, p, &alpha[2], &a[a_offset], lda, 
	    &b[b_offset], ldb, &c__[c_offset], ldc, &nfp, &dwork[kz], n, &
	    dwork[kwr], &dwork[kwi], &dwork[kw], &i__1, info, (ftnlen)1, (
	    ftnlen)6, (ftnlen)7);
#line 362 "SB08FD.f"
    if (*info != 0) {
#line 362 "SB08FD.f"
	return 0;
#line 362 "SB08FD.f"
    }

#line 365 "SB08FD.f"
    wrkopt = dwork[kw] + (doublereal) (kw - 1);

/*     Perform the pole assignment if there exist "unstable" eigenvalues. */

#line 369 "SB08FD.f"
    *nq = *n;
#line 370 "SB08FD.f"
    if (nfp < *n) {
#line 371 "SB08FD.f"
	kg = 1;
#line 372 "SB08FD.f"
	kfi = kg + (*m << 1);
#line 373 "SB08FD.f"
	kw = kfi + (*m << 1);

/*        Set the limits for the bottom diagonal block. */

#line 377 "SB08FD.f"
	nlow = nfp + 1;
#line 378 "SB08FD.f"
	nsup = *n;

/*        WHILE (NLOW <= NSUP) DO */
#line 381 "SB08FD.f"
L10:
#line 381 "SB08FD.f"
	if (nlow <= nsup) {

/*           Main loop for assigning one or two poles. */

/*           Determine the dimension of the last block. */

#line 387 "SB08FD.f"
	    ib = 1;
#line 388 "SB08FD.f"
	    if (nlow < nsup) {
#line 389 "SB08FD.f"
		if (a[nsup + (nsup - 1) * a_dim1] != 0.) {
#line 389 "SB08FD.f"
		    ib = 2;
#line 389 "SB08FD.f"
		}
#line 390 "SB08FD.f"
	    }
#line 391 "SB08FD.f"
	    l = nsup - ib + 1;

/*           Save the last IB rows of B in G. */

#line 395 "SB08FD.f"
	    dlacpy_("Full", &ib, m, &b[l + b_dim1], ldb, &dwork[kg], &ib, (
		    ftnlen)4);

/*           Check the controllability of the last block. */

#line 399 "SB08FD.f"
	    if (dlange_("1-norm", &ib, m, &dwork[kg], &ib, &dwork[kw], (
		    ftnlen)6) <= toler) {

/*              Deflate the uncontrollable block and resume the */
/*              main loop. */

#line 405 "SB08FD.f"
		nsup -= ib;
#line 406 "SB08FD.f"
	    } else {

/*              Form the IBxIB matrix A2 from the last diagonal block and */
/*              set the pole(s) to be assigned. */

#line 411 "SB08FD.f"
		a2[0] = a[l + l * a_dim1];
#line 412 "SB08FD.f"
		if (ib == 1) {
#line 413 "SB08FD.f"
		    sm = alpha[1];
#line 414 "SB08FD.f"
		    if (discr) {
#line 414 "SB08FD.f"
			sm = d_sign(&alpha[1], a2);
#line 414 "SB08FD.f"
		    }
#line 415 "SB08FD.f"
		    pr = alpha[1];
#line 416 "SB08FD.f"
		} else {
#line 417 "SB08FD.f"
		    a2[2] = a[l + nsup * a_dim1];
#line 418 "SB08FD.f"
		    a2[1] = a[nsup + l * a_dim1];
#line 419 "SB08FD.f"
		    a2[3] = a[nsup + nsup * a_dim1];
#line 420 "SB08FD.f"
		    sm = alpha[1] + alpha[1];
#line 421 "SB08FD.f"
		    pr = alpha[1] * alpha[1];
#line 422 "SB08FD.f"
		    if (discr) {
#line 423 "SB08FD.f"
			x = a2[0];
#line 424 "SB08FD.f"
			y = sqrt((d__1 = a2[2] * a2[1], abs(d__1)));
#line 425 "SB08FD.f"
			sm = sm * x / dlapy2_(&x, &y);
#line 426 "SB08FD.f"
		    } else {
#line 427 "SB08FD.f"
			pr -= a2[2] * a2[1];
#line 428 "SB08FD.f"
		    }
#line 429 "SB08FD.f"
		}

/*              Determine the M-by-IB feedback matrix FI which assigns */
/*              the selected IB poles for the pair (A2,G). */

/*              Workspace needed: 5*M. */

#line 436 "SB08FD.f"
		sb01by_(&ib, m, &sm, &pr, a2, &dwork[kg], &dwork[kfi], &toler,
			 &dwork[kw], info);
#line 438 "SB08FD.f"
		if (*info != 0) {

/*                 Uncontrollable 2x2 block with double real eigenvalues */
/*                 which due to roundoff appear as a pair of complex */
/*                 conjugated eigenvalues. */
/*                 One of them can be elliminated using the information */
/*                 in DWORK(KFI) and DWORK(KFI+M). */

#line 446 "SB08FD.f"
		    cs = dwork[kfi];
#line 447 "SB08FD.f"
		    sn = -dwork[kfi + *m];

/*                 Apply the Givens transformation to A, B, C and F. */

#line 451 "SB08FD.f"
		    l1 = l + 1;
#line 452 "SB08FD.f"
		    i__1 = nsup - l + 1;
#line 452 "SB08FD.f"
		    drot_(&i__1, &a[l1 + l * a_dim1], lda, &a[l + l * a_dim1],
			     lda, &cs, &sn);
#line 454 "SB08FD.f"
		    drot_(&l1, &a[l1 * a_dim1 + 1], &c__1, &a[l * a_dim1 + 1],
			     &c__1, &cs, &sn);
#line 455 "SB08FD.f"
		    drot_(m, &b[l1 + b_dim1], ldb, &b[l + b_dim1], ldb, &cs, &
			    sn);
#line 456 "SB08FD.f"
		    if (*p > 0) {
#line 456 "SB08FD.f"
			drot_(p, &c__[l1 * c_dim1 + 1], &c__1, &c__[l * 
				c_dim1 + 1], &c__1, &cs, &sn);
#line 456 "SB08FD.f"
		    }
#line 458 "SB08FD.f"
		    drot_(m, &cr[l1 * cr_dim1 + 1], &c__1, &cr[l * cr_dim1 + 
			    1], &c__1, &cs, &sn);

/*                 Deflate the uncontrollable block and resume the */
/*                 main loop. */

#line 463 "SB08FD.f"
		    a[l1 + l * a_dim1] = 0.;
#line 464 "SB08FD.f"
		    --nsup;
#line 465 "SB08FD.f"
		    *info = 0;
#line 466 "SB08FD.f"
		    goto L10;
#line 467 "SB08FD.f"
		}

/*              Check for possible numerical instability. */

#line 471 "SB08FD.f"
		if (dlange_("1-norm", m, &ib, &dwork[kfi], m, &dwork[kw], (
			ftnlen)6) > rmax) {
#line 471 "SB08FD.f"
		    ++(*iwarn);
#line 471 "SB08FD.f"
		}

/*              Update the feedback matrix F <-- F + [0 FI] in CR. */

#line 476 "SB08FD.f"
		k = kfi;
#line 477 "SB08FD.f"
		i__1 = l + ib - 1;
#line 477 "SB08FD.f"
		for (j = l; j <= i__1; ++j) {
#line 478 "SB08FD.f"
		    i__2 = *m;
#line 478 "SB08FD.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 479 "SB08FD.f"
			cr[i__ + j * cr_dim1] += dwork[k];
#line 480 "SB08FD.f"
			++k;
#line 481 "SB08FD.f"
/* L20: */
#line 481 "SB08FD.f"
		    }
#line 482 "SB08FD.f"
/* L30: */
#line 482 "SB08FD.f"
		}

/*              Update the state matrix A <-- A + B*[0 FI]. */

#line 486 "SB08FD.f"
		dgemm_("NoTranspose", "NoTranspose", &nsup, &ib, m, &c_b7, &b[
			b_offset], ldb, &dwork[kfi], m, &c_b7, &a[l * a_dim1 
			+ 1], lda, (ftnlen)11, (ftnlen)11);
#line 489 "SB08FD.f"
		if (ib == 2) {

/*                 Try to split the 2x2 block and standardize it. */

#line 493 "SB08FD.f"
		    l1 = l + 1;
#line 494 "SB08FD.f"
		    dlanv2_(&a[l + l * a_dim1], &a[l + l1 * a_dim1], &a[l1 + 
			    l * a_dim1], &a[l1 + l1 * a_dim1], &x, &y, &pr, &
			    sm, &cs, &sn);

/*                 Apply the transformation to A, B, C and F. */

#line 499 "SB08FD.f"
		    if (l1 < nsup) {
#line 499 "SB08FD.f"
			i__1 = nsup - l1;
#line 499 "SB08FD.f"
			drot_(&i__1, &a[l + (l1 + 1) * a_dim1], lda, &a[l1 + (
				l1 + 1) * a_dim1], lda, &cs, &sn);
#line 499 "SB08FD.f"
		    }
#line 502 "SB08FD.f"
		    i__1 = l - 1;
#line 502 "SB08FD.f"
		    drot_(&i__1, &a[l * a_dim1 + 1], &c__1, &a[l1 * a_dim1 + 
			    1], &c__1, &cs, &sn);
#line 503 "SB08FD.f"
		    drot_(m, &b[l + b_dim1], ldb, &b[l1 + b_dim1], ldb, &cs, &
			    sn);
#line 504 "SB08FD.f"
		    if (*p > 0) {
#line 504 "SB08FD.f"
			drot_(p, &c__[l * c_dim1 + 1], &c__1, &c__[l1 * 
				c_dim1 + 1], &c__1, &cs, &sn);
#line 504 "SB08FD.f"
		    }
#line 506 "SB08FD.f"
		    drot_(m, &cr[l * cr_dim1 + 1], &c__1, &cr[l1 * cr_dim1 + 
			    1], &c__1, &cs, &sn);
#line 507 "SB08FD.f"
		}
#line 508 "SB08FD.f"
		if (nlow + ib <= nsup) {

/*                 Move the last block(s) to the leading position(s) of */
/*                 the bottom block. */

/*                 Workspace:     need MAX(4*N, 4*M, 4*P). */

#line 515 "SB08FD.f"
		    ncur1 = nsup - ib;
#line 516 "SB08FD.f"
		    nmoves = 1;
#line 517 "SB08FD.f"
		    if (ib == 2 && a[nsup + (nsup - 1) * a_dim1] == 0.) {
#line 518 "SB08FD.f"
			ib = 1;
#line 519 "SB08FD.f"
			nmoves = 2;
#line 520 "SB08FD.f"
		    }

/*                 WHILE (NMOVES > 0) DO */
#line 523 "SB08FD.f"
L40:
#line 523 "SB08FD.f"
		    if (nmoves > 0) {
#line 524 "SB08FD.f"
			ncur = ncur1;

/*                    WHILE (NCUR >= NLOW) DO */
#line 527 "SB08FD.f"
L50:
#line 527 "SB08FD.f"
			if (ncur >= nlow) {

/*                       Loop for positioning of the last block. */

/*                       Determine the dimension of the current block. */

#line 533 "SB08FD.f"
			    ib1 = 1;
#line 534 "SB08FD.f"
			    if (ncur > nlow) {
#line 535 "SB08FD.f"
				if (a[ncur + (ncur - 1) * a_dim1] != 0.) {
#line 535 "SB08FD.f"
				    ib1 = 2;
#line 535 "SB08FD.f"
				}
#line 536 "SB08FD.f"
			    }
#line 537 "SB08FD.f"
			    nb = ib1 + ib;

/*                       Initialize the local transformation matrix Z. */

#line 541 "SB08FD.f"
			    dlaset_("Full", &nb, &nb, &c_b6, &c_b7, z__, &
				    c__4, (ftnlen)4);
#line 542 "SB08FD.f"
			    l = ncur - ib1 + 1;

/*                       Exchange two adjacent blocks and accumulate the */
/*                       transformations in Z. */

#line 547 "SB08FD.f"
			    dlaexc_(&c_true, &nb, &a[l + l * a_dim1], lda, 
				    z__, &c__4, &c__1, &ib1, &ib, &dwork[1], 
				    info);
#line 549 "SB08FD.f"
			    if (*info != 0) {
#line 550 "SB08FD.f"
				*info = 2;
#line 551 "SB08FD.f"
				return 0;
#line 552 "SB08FD.f"
			    }

/*                       Apply the transformation to the rest of A. */

#line 556 "SB08FD.f"
			    l1 = l + nb;
#line 557 "SB08FD.f"
			    if (l1 <= nsup) {
#line 558 "SB08FD.f"
				i__1 = nsup - l1 + 1;
#line 558 "SB08FD.f"
				dgemm_("Transpose", "NoTranspose", &nb, &i__1,
					 &nb, &c_b7, z__, &c__4, &a[l + l1 * 
					a_dim1], lda, &c_b6, &dwork[1], &nb, (
					ftnlen)9, (ftnlen)11);
#line 561 "SB08FD.f"
				i__1 = nsup - l1 + 1;
#line 561 "SB08FD.f"
				dlacpy_("Full", &nb, &i__1, &dwork[1], &nb, &
					a[l + l1 * a_dim1], lda, (ftnlen)4);
#line 563 "SB08FD.f"
			    }
#line 564 "SB08FD.f"
			    i__1 = l - 1;
#line 564 "SB08FD.f"
			    dgemm_("NoTranspose", "NoTranspose", &i__1, &nb, &
				    nb, &c_b7, &a[l * a_dim1 + 1], lda, z__, &
				    c__4, &c_b6, &dwork[1], n, (ftnlen)11, (
				    ftnlen)11);
#line 567 "SB08FD.f"
			    i__1 = l - 1;
#line 567 "SB08FD.f"
			    dlacpy_("Full", &i__1, &nb, &dwork[1], n, &a[l * 
				    a_dim1 + 1], lda, (ftnlen)4);

/*                       Apply the transformation to B, C and F. */

#line 572 "SB08FD.f"
			    dgemm_("Transpose", "NoTranspose", &nb, m, &nb, &
				    c_b7, z__, &c__4, &b[l + b_dim1], ldb, &
				    c_b6, &dwork[1], &nb, (ftnlen)9, (ftnlen)
				    11);
#line 575 "SB08FD.f"
			    dlacpy_("Full", &nb, m, &dwork[1], &nb, &b[l + 
				    b_dim1], ldb, (ftnlen)4);

#line 578 "SB08FD.f"
			    if (*p > 0) {
#line 579 "SB08FD.f"
				dgemm_("NoTranspose", "NoTranspose", p, &nb, &
					nb, &c_b7, &c__[l * c_dim1 + 1], ldc, 
					z__, &c__4, &c_b6, &dwork[1], p, (
					ftnlen)11, (ftnlen)11);
#line 582 "SB08FD.f"
				dlacpy_("Full", p, &nb, &dwork[1], p, &c__[l *
					 c_dim1 + 1], ldc, (ftnlen)4);
#line 584 "SB08FD.f"
			    }

#line 586 "SB08FD.f"
			    dgemm_("NoTranspose", "NoTranspose", m, &nb, &nb, 
				    &c_b7, &cr[l * cr_dim1 + 1], ldcr, z__, &
				    c__4, &c_b6, &dwork[1], m, (ftnlen)11, (
				    ftnlen)11);
#line 589 "SB08FD.f"
			    dlacpy_("Full", m, &nb, &dwork[1], m, &cr[l * 
				    cr_dim1 + 1], ldcr, (ftnlen)4);

#line 592 "SB08FD.f"
			    ncur -= ib1;
#line 593 "SB08FD.f"
			    goto L50;
#line 594 "SB08FD.f"
			}
/*                    END WHILE 50 */

#line 597 "SB08FD.f"
			--nmoves;
#line 598 "SB08FD.f"
			++ncur1;
#line 599 "SB08FD.f"
			nlow += ib;
#line 600 "SB08FD.f"
			goto L40;
#line 601 "SB08FD.f"
		    }
/*                 END WHILE 40 */

#line 604 "SB08FD.f"
		} else {
#line 605 "SB08FD.f"
		    nlow += ib;
#line 606 "SB08FD.f"
		}
#line 607 "SB08FD.f"
	    }
#line 608 "SB08FD.f"
	    goto L10;
#line 609 "SB08FD.f"
	}
/*        END WHILE 10 */

#line 612 "SB08FD.f"
	*nq = nsup;
#line 613 "SB08FD.f"
	*nr = nsup - nfp;

/*        Annihilate the elements below the first subdiagonal of A. */

#line 617 "SB08FD.f"
	if (*nq > 2) {
#line 617 "SB08FD.f"
	    i__1 = *nq - 2;
#line 617 "SB08FD.f"
	    i__2 = *nq - 2;
#line 617 "SB08FD.f"
	    dlaset_("Lower", &i__1, &i__2, &c_b6, &c_b6, &a[a_dim1 + 3], lda, 
		    (ftnlen)5);
#line 617 "SB08FD.f"
	}
#line 619 "SB08FD.f"
    }

/*     Compute C <-- CQ = C + D*F. */

#line 623 "SB08FD.f"
    dgemm_("NoTranspose", "NoTranspose", p, nq, m, &c_b7, &d__[d_offset], ldd,
	     &cr[cr_offset], ldcr, &c_b7, &c__[c_offset], ldc, (ftnlen)11, (
	    ftnlen)11);

/* Computing MAX */
/* Computing MAX */
#line 626 "SB08FD.f"
    i__1 = *m * 5, i__2 = *p << 2;
#line 626 "SB08FD.f"
    d__1 = wrkopt, d__2 = (doublereal) max(i__1,i__2);
#line 626 "SB08FD.f"
    dwork[1] = max(d__1,d__2);

#line 628 "SB08FD.f"
    return 0;
/* *** Last line of SB08FD *** */
} /* sb08fd_ */

