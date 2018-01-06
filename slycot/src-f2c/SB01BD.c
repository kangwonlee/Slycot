#line 1 "SB01BD.f"
/* SB01BD.f -- translated by f2c (version 20100827).
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

#line 1 "SB01BD.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b14 = 0.;
static doublereal c_b21 = 1.;
static logical c_true = TRUE_;
static integer c__2 = 2;
static logical c_false = FALSE_;

/* Subroutine */ int sb01bd_(char *dico, integer *n, integer *m, integer *np, 
	doublereal *alpha, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *wr, doublereal *wi, integer *nfp, integer *
	nap, integer *nup, doublereal *f, integer *ldf, doublereal *z__, 
	integer *ldz, doublereal *tol, doublereal *dwork, integer *ldwork, 
	integer *iwarn, integer *info, ftnlen dico_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, f_dim1, f_offset, z_dim1, 
	    z_offset, i__1, i__2;

    /* Local variables */
    static doublereal c__;
    static integer i__, j, k;
    static doublereal p, s, x, y, a2[4]	/* was [2][2] */;
    static integer ib, kg, nl, kw, ib1, kfi, ipc, npc, kwi, npr, kwr;
    static logical ceig;
    static integer ierr, ncur;
    static doublereal rmax;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer nlow, nsup, ncur1;
    extern /* Subroutine */ int mb03qd_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen), dgees_(char *, char *, L_fp, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen), dgemm_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen), sb01bx_(logical *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), sb01by_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    extern /* Subroutine */ int mb03qy_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal anorm, bnorm;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical bwork[1];
    static doublereal toler;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlaexc_(logical *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, doublereal *, integer *), dlaset_(char *, integer *, integer *,
	     doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern logical select_();
    static logical simplb;
    static doublereal tolerb;
    static integer nmoves, wrkopt;


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

/*     To determine the state feedback matrix F for a given system (A,B) */
/*     such that the closed-loop state matrix A+B*F has specified */
/*     eigenvalues. */

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
/*             the number of columns of the matrix F.  N >= 0. */

/*     M       (input) INTEGER */
/*             The dimension of input vector, i.e. the number of columns */
/*             of the matrix B and the number of rows of the matrix F. */
/*             M >= 0. */

/*     NP      (input) INTEGER */
/*             The number of given eigenvalues. At most N eigenvalues */
/*             can be assigned.  0 <= NP. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             Specifies the maximum admissible value, either for real */
/*             parts, if DICO = 'C', or for moduli, if DICO = 'D', */
/*             of the eigenvalues of A which will not be modified by */
/*             the eigenvalue assignment algorithm. */
/*             ALPHA >= 0 if DICO = 'D'. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the matrix Z'*(A+B*F)*Z in a real Schur form. */
/*             The leading NFP-by-NFP diagonal block of A corresponds */
/*             to the fixed (unmodified) eigenvalues having real parts */
/*             less than ALPHA, if DICO = 'C', or moduli less than ALPHA, */
/*             if DICO = 'D'. The trailing NUP-by-NUP diagonal block of A */
/*             corresponds to the uncontrollable eigenvalues detected by */
/*             the eigenvalue assignment algorithm. The elements under */
/*             the first subdiagonal are set to zero. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             input/state matrix. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     WR,WI   (input/output) DOUBLE PRECISION array, dimension (NP) */
/*             On entry, these arrays must contain the real and imaginary */
/*             parts, respectively, of the desired eigenvalues of the */
/*             closed-loop system state-matrix A+B*F. The eigenvalues */
/*             can be unordered, except that complex conjugate pairs */
/*             must appear consecutively in these arrays. */
/*             On exit, if INFO = 0, the leading NAP elements of these */
/*             arrays contain the real and imaginary parts, respectively, */
/*             of the assigned eigenvalues. The trailing NP-NAP elements */
/*             contain the unassigned eigenvalues. */

/*     NFP     (output) INTEGER */
/*             The number of eigenvalues of A having real parts less than */
/*             ALPHA, if DICO = 'C', or moduli less than ALPHA, if */
/*             DICO = 'D'. These eigenvalues are not modified by the */
/*             eigenvalue assignment algorithm. */

/*     NAP     (output) INTEGER */
/*             The number of assigned eigenvalues. If INFO = 0 on exit, */
/*             then NAP = N-NFP-NUP. */

/*     NUP     (output) INTEGER */
/*             The number of uncontrollable eigenvalues detected by the */
/*             eigenvalue assignment algorithm (see METHOD). */

/*     F       (output) DOUBLE PRECISION array, dimension (LDF,N) */
/*             The leading M-by-N part of this array contains the state */
/*             feedback F, which assigns NAP closed-loop eigenvalues and */
/*             keeps unaltered N-NAP open-loop eigenvalues. */

/*     LDF     INTEGER */
/*             The leading dimension of array F.  LDF >= MAX(1,M). */

/*     Z       (output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             The leading N-by-N part of this array contains the */
/*             orthogonal matrix Z which reduces the closed-loop */
/*             system state matrix A + B*F to upper real Schur form. */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z.  LDZ >= MAX(1,N). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The absolute tolerance level below which the elements of A */
/*             or B are considered zero (used for controllability tests). */
/*             If the user sets TOL <= 0, then the default tolerance */
/*             TOL = N * EPS * max(NORM(A),NORM(B)) is used, where EPS is */
/*             the machine precision (see LAPACK Library routine DLAMCH) */
/*             and NORM(A) denotes the 1-norm of A. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of working array DWORK. */
/*             LDWORK >= MAX( 1,5*M,5*N,2*N+4*M ). */
/*             For optimum performance LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = K:  K violations of the numerical stability condition */
/*                   NORM(F) <= 100*NORM(A)/NORM(B) occured during the */
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
/*             = 3:  the number of eigenvalues to be assigned is less */
/*                   than the number of possibly assignable eigenvalues; */
/*                   NAP eigenvalues have been properly assigned, */
/*                   but some assignable eigenvalues remain unmodified. */
/*             = 4:  an attempt is made to place a complex conjugate */
/*                   pair on the location of a real eigenvalue. This */
/*                   situation can only appear when N-NFP is odd, */
/*                   NP > N-NFP-NUP is even, and for the last real */
/*                   eigenvalue to be modified there exists no available */
/*                   real eigenvalue to be assigned. However, NAP */
/*                   eigenvalues have been already properly assigned. */

/*     METHOD */

/*     SB01BD is based on the factorization algorithm of [1]. */
/*     Given the matrices A and B of dimensions N-by-N and N-by-M, */
/*     respectively, this subroutine constructs an M-by-N matrix F such */
/*     that A + BF has eigenvalues as follows. */
/*     Let NFP eigenvalues of A have real parts less than ALPHA, if */
/*     DICO = 'C', or moduli less then ALPHA, if DICO = 'D'. Then: */
/*     1) If the pair (A,B) is controllable, then A + B*F has */
/*        NAP = MIN(NP,N-NFP) eigenvalues assigned from those specified */
/*        by WR + j*WI and N-NAP unmodified eigenvalues; */
/*     2) If the pair (A,B) is uncontrollable, then the number of */
/*        assigned eigenvalues NAP satifies generally the condition */
/*        NAP <= MIN(NP,N-NFP). */

/*     At the beginning of the algorithm, F = 0 and the matrix A is */
/*     reduced to an ordered real Schur form by separating its spectrum */
/*     in two parts. The leading NFP-by-NFP part of the Schur form of */
/*     A corresponds to the eigenvalues which will not be modified. */
/*     These eigenvalues have real parts less than ALPHA, if */
/*     DICO = 'C', or moduli less than ALPHA, if DICO = 'D'. */
/*     The performed orthogonal transformations are accumulated in Z. */
/*     After this preliminary reduction, the algorithm proceeds */
/*     recursively. */

/*     Let F be the feedback matrix at the beginning of a typical step i. */
/*     At each step of the algorithm one real eigenvalue or two complex */
/*     conjugate eigenvalues are placed by a feedback Fi of rank 1 or */
/*     rank 2, respectively. Since the feedback Fi affects only the */
/*     last 1 or 2 columns of Z'*(A+B*F)*Z, the matrix Z'*(A+B*F+B*Fi)*Z */
/*     therefore remains in real Schur form. The assigned eigenvalue(s) */
/*     is (are) then moved to another diagonal position of the real */
/*     Schur form using reordering techniques and a new block is */
/*     transfered in the last diagonal position. The feedback matrix F */
/*     is updated as F <-- F + Fi. The eigenvalue(s) to be assigned at */
/*     each step is (are) chosen such that the norm of each Fi is */
/*     minimized. */

/*     If uncontrollable eigenvalues are encountered in the last diagonal */
/*     position of the real Schur matrix Z'*(A+B*F)*Z, the algorithm */
/*     deflates them at the bottom of the real Schur form and redefines */
/*     accordingly the position of the "last" block. */

/*     Note: Not all uncontrollable eigenvalues of the pair (A,B) are */
/*     necessarily detected by the eigenvalue assignment algorithm. */
/*     Undetected uncontrollable eigenvalues may exist if NFP > 0 and/or */
/*     NP < N-NFP. */

/*     REFERENCES */

/*     [1] Varga A. */
/*         A Schur method for pole assignment. */
/*         IEEE Trans. Autom. Control, Vol. AC-26, pp. 517-519, 1981. */

/*     NUMERICAL ASPECTS */
/*                                            3 */
/*     The algorithm requires no more than 14N  floating point */
/*     operations. Although no proof of numerical stability is known, */
/*     the algorithm has always been observed to yield reliable */
/*     numerical results. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     February 1999. Based on the RASP routine SB01BD. */

/*     REVISIONS */

/*     March 30, 1999, V. Sima, Research Institute for Informatics, */
/*     Bucharest. */
/*     April 4, 1999. A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen. */
/*     May 18, 2003. A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen. */
/*     Feb. 15, 2004, V. Sima, Research Institute for Informatics, */
/*     Bucharest. */
/*     May 12, 2005. A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen. */

/*     KEYWORDS */

/*     Eigenvalues, eigenvalue assignment, feedback control, */
/*     pole placement, state-space model. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 296 "SB01BD.f"
    /* Parameter adjustments */
#line 296 "SB01BD.f"
    a_dim1 = *lda;
#line 296 "SB01BD.f"
    a_offset = 1 + a_dim1;
#line 296 "SB01BD.f"
    a -= a_offset;
#line 296 "SB01BD.f"
    b_dim1 = *ldb;
#line 296 "SB01BD.f"
    b_offset = 1 + b_dim1;
#line 296 "SB01BD.f"
    b -= b_offset;
#line 296 "SB01BD.f"
    --wr;
#line 296 "SB01BD.f"
    --wi;
#line 296 "SB01BD.f"
    f_dim1 = *ldf;
#line 296 "SB01BD.f"
    f_offset = 1 + f_dim1;
#line 296 "SB01BD.f"
    f -= f_offset;
#line 296 "SB01BD.f"
    z_dim1 = *ldz;
#line 296 "SB01BD.f"
    z_offset = 1 + z_dim1;
#line 296 "SB01BD.f"
    z__ -= z_offset;
#line 296 "SB01BD.f"
    --dwork;
#line 296 "SB01BD.f"

#line 296 "SB01BD.f"
    /* Function Body */
#line 296 "SB01BD.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 297 "SB01BD.f"
    *iwarn = 0;
#line 298 "SB01BD.f"
    *info = 0;

/*     Check the scalar input parameters. */

#line 302 "SB01BD.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 303 "SB01BD.f"
	*info = -1;
#line 304 "SB01BD.f"
    } else if (*n < 0) {
#line 305 "SB01BD.f"
	*info = -2;
#line 306 "SB01BD.f"
    } else if (*m < 0) {
#line 307 "SB01BD.f"
	*info = -3;
#line 308 "SB01BD.f"
    } else if (*np < 0) {
#line 309 "SB01BD.f"
	*info = -4;
#line 310 "SB01BD.f"
    } else if (discr && *alpha < 0.) {
#line 311 "SB01BD.f"
	*info = -5;
#line 312 "SB01BD.f"
    } else if (*lda < max(1,*n)) {
#line 313 "SB01BD.f"
	*info = -7;
#line 314 "SB01BD.f"
    } else if (*ldb < max(1,*n)) {
#line 315 "SB01BD.f"
	*info = -9;
#line 316 "SB01BD.f"
    } else if (*ldf < max(1,*m)) {
#line 317 "SB01BD.f"
	*info = -16;
#line 318 "SB01BD.f"
    } else if (*ldz < max(1,*n)) {
#line 319 "SB01BD.f"
	*info = -18;
#line 320 "SB01BD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 320 "SB01BD.f"
	i__1 = 1, i__2 = *m * 5, i__1 = max(i__1,i__2), i__2 = *n * 5, i__1 = 
		max(i__1,i__2), i__2 = (*n << 1) + (*m << 2);
#line 320 "SB01BD.f"
	if (*ldwork < max(i__1,i__2)) {
#line 321 "SB01BD.f"
	    *info = -21;
#line 322 "SB01BD.f"
	}
#line 322 "SB01BD.f"
    }
#line 323 "SB01BD.f"
    if (*info != 0) {

/*        Error return. */

#line 327 "SB01BD.f"
	i__1 = -(*info);
#line 327 "SB01BD.f"
	xerbla_("SB01BD", &i__1, (ftnlen)6);
#line 328 "SB01BD.f"
	return 0;
#line 329 "SB01BD.f"
    }

/*     Quick return if possible. */

#line 333 "SB01BD.f"
    if (*n == 0) {
#line 334 "SB01BD.f"
	*nfp = 0;
#line 335 "SB01BD.f"
	*nap = 0;
#line 336 "SB01BD.f"
	*nup = 0;
#line 337 "SB01BD.f"
	dwork[1] = 1.;
#line 338 "SB01BD.f"
	return 0;
#line 339 "SB01BD.f"
    }

/*     Compute the norms of A and B, and set default tolerances */
/*     if necessary. */

#line 344 "SB01BD.f"
    anorm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)6);
#line 345 "SB01BD.f"
    bnorm = dlange_("1-norm", n, m, &b[b_offset], ldb, &dwork[1], (ftnlen)6);
#line 346 "SB01BD.f"
    if (*tol <= 0.) {
#line 347 "SB01BD.f"
	x = dlamch_("Epsilon", (ftnlen)7);
#line 348 "SB01BD.f"
	toler = (doublereal) (*n) * max(anorm,bnorm) * x;
#line 349 "SB01BD.f"
	tolerb = (doublereal) (*n) * bnorm * x;
#line 350 "SB01BD.f"
    } else {
#line 351 "SB01BD.f"
	toler = *tol;
#line 352 "SB01BD.f"
	tolerb = *tol;
#line 353 "SB01BD.f"
    }

/*     Allocate working storage. */

#line 357 "SB01BD.f"
    kwr = 1;
#line 358 "SB01BD.f"
    kwi = kwr + *n;
#line 359 "SB01BD.f"
    kw = kwi + *n;

/*     Reduce A to real Schur form using an orthogonal similarity */
/*     transformation A <- Z'*A*Z and accumulate the transformation in Z. */

/*     Workspace:  need   5*N; */
/*                 prefer larger. */

#line 367 "SB01BD.f"
    i__1 = *ldwork - kw + 1;
#line 367 "SB01BD.f"
    dgees_("Vectors", "No ordering", (L_fp)select_, n, &a[a_offset], lda, &
	    ncur, &dwork[kwr], &dwork[kwi], &z__[z_offset], ldz, &dwork[kw], &
	    i__1, bwork, info, (ftnlen)7, (ftnlen)11);
#line 370 "SB01BD.f"
    wrkopt = kw - 1 + (integer) dwork[kw];
#line 371 "SB01BD.f"
    if (*info != 0) {
#line 372 "SB01BD.f"
	*info = 1;
#line 373 "SB01BD.f"
	return 0;
#line 374 "SB01BD.f"
    }

/*     Reduce A to an ordered real Schur form using an orthogonal */
/*     similarity transformation A <- Z'*A*Z and accumulate the */
/*     transformations in Z. The separation of the spectrum of A is */
/*     performed such that the leading NFP-by-NFP submatrix of A */
/*     corresponds to the "good" eigenvalues which will not be */
/*     modified. The bottom (N-NFP)-by-(N-NFP) diagonal block of A */
/*     corresponds to the "bad" eigenvalues to be modified. */

/*     Workspace needed:  N. */

#line 386 "SB01BD.f"
    mb03qd_(dico, "Stable", "Update", n, &c__1, n, alpha, &a[a_offset], lda, &
	    z__[z_offset], ldz, nfp, &dwork[1], info, (ftnlen)1, (ftnlen)6, (
	    ftnlen)6);
#line 388 "SB01BD.f"
    if (*info != 0) {
#line 388 "SB01BD.f"
	return 0;
#line 388 "SB01BD.f"
    }

/*     Set F = 0. */

#line 393 "SB01BD.f"
    dlaset_("Full", m, n, &c_b14, &c_b14, &f[f_offset], ldf, (ftnlen)4);

/*     Return if B is negligible (uncontrollable system). */

#line 397 "SB01BD.f"
    if (bnorm <= tolerb) {
#line 398 "SB01BD.f"
	*nap = 0;
#line 399 "SB01BD.f"
	*nup = *n;
#line 400 "SB01BD.f"
	dwork[1] = (doublereal) wrkopt;
#line 401 "SB01BD.f"
	return 0;
#line 402 "SB01BD.f"
    }

/*     Compute the bound for the numerical stability condition. */

#line 406 "SB01BD.f"
    rmax = anorm * 100. / bnorm;

/*     Perform eigenvalue assignment if there exist "bad" eigenvalues. */

#line 410 "SB01BD.f"
    *nap = 0;
#line 411 "SB01BD.f"
    *nup = 0;
#line 412 "SB01BD.f"
    if (*nfp < *n) {
#line 413 "SB01BD.f"
	kg = 1;
#line 414 "SB01BD.f"
	kfi = kg + (*m << 1);
#line 415 "SB01BD.f"
	kw = kfi + (*m << 1);

/*        Set the limits for the bottom diagonal block. */

#line 419 "SB01BD.f"
	nlow = *nfp + 1;
#line 420 "SB01BD.f"
	nsup = *n;

/*        Separate and count real and complex eigenvalues to be assigned. */

#line 424 "SB01BD.f"
	npr = 0;
#line 425 "SB01BD.f"
	i__1 = *np;
#line 425 "SB01BD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 426 "SB01BD.f"
	    if (wi[i__] == 0.) {
#line 427 "SB01BD.f"
		++npr;
#line 428 "SB01BD.f"
		k = i__ - npr;
#line 429 "SB01BD.f"
		if (k > 0) {
#line 430 "SB01BD.f"
		    s = wr[i__];
#line 431 "SB01BD.f"
		    i__2 = npr;
#line 431 "SB01BD.f"
		    for (j = npr + k - 1; j >= i__2; --j) {
#line 432 "SB01BD.f"
			wr[j + 1] = wr[j];
#line 433 "SB01BD.f"
			wi[j + 1] = wi[j];
#line 434 "SB01BD.f"
/* L5: */
#line 434 "SB01BD.f"
		    }
#line 435 "SB01BD.f"
		    wr[npr] = s;
#line 436 "SB01BD.f"
		    wi[npr] = 0.;
#line 437 "SB01BD.f"
		}
#line 438 "SB01BD.f"
	    }
#line 439 "SB01BD.f"
/* L10: */
#line 439 "SB01BD.f"
	}
#line 440 "SB01BD.f"
	npc = *np - npr;

/*        The first NPR elements of WR and WI contain the real */
/*        eigenvalues, the last NPC elements contain the complex */
/*        eigenvalues. Set the pointer to complex eigenvalues. */

#line 446 "SB01BD.f"
	ipc = npr + 1;

/*        Main loop for assigning one or two eigenvalues. */

/*        Terminate if all eigenvalues were assigned, or if there */
/*        are no more eigenvalues to be assigned, or if a non-fatal */
/*        error condition was set. */

/*        WHILE (NLOW <= NSUP and INFO = 0) DO */

#line 456 "SB01BD.f"
L20:
#line 456 "SB01BD.f"
	if (nlow <= nsup && *info == 0) {

/*           Determine the dimension of the last block. */

#line 460 "SB01BD.f"
	    ib = 1;
#line 461 "SB01BD.f"
	    if (nlow < nsup) {
#line 462 "SB01BD.f"
		if (a[nsup + (nsup - 1) * a_dim1] != 0.) {
#line 462 "SB01BD.f"
		    ib = 2;
#line 462 "SB01BD.f"
		}
#line 463 "SB01BD.f"
	    }

/*           Compute G, the current last IB rows of Z'*B. */

#line 467 "SB01BD.f"
	    nl = nsup - ib + 1;
#line 468 "SB01BD.f"
	    dgemm_("Transpose", "NoTranspose", &ib, m, n, &c_b21, &z__[nl * 
		    z_dim1 + 1], ldz, &b[b_offset], ldb, &c_b14, &dwork[kg], &
		    ib, (ftnlen)9, (ftnlen)11);

/*           Check the controllability for a simple block. */

#line 473 "SB01BD.f"
	    if (dlange_("1", &ib, m, &dwork[kg], &ib, &dwork[kw], (ftnlen)1) 
		    <= tolerb) {

/*              Deflate the uncontrollable block and resume the */
/*              main loop. */

#line 479 "SB01BD.f"
		nsup -= ib;
#line 480 "SB01BD.f"
		*nup += ib;
#line 481 "SB01BD.f"
		goto L20;
#line 482 "SB01BD.f"
	    }

/*           Test for termination with INFO = 3. */

#line 486 "SB01BD.f"
	    if (*nap == *np) {
#line 487 "SB01BD.f"
		*info = 3;

/*              Test for compatibility. Terminate if an attempt occurs */
/*              to place a complex conjugate pair on a 1x1 block. */

#line 492 "SB01BD.f"
	    } else if (ib == 1 && npr == 0 && nlow == nsup) {
#line 493 "SB01BD.f"
		*info = 4;
#line 494 "SB01BD.f"
	    } else {

/*              Set the simple block flag. */

#line 498 "SB01BD.f"
		simplb = TRUE_;

/*              Form a 2-by-2 block if necessary from two 1-by-1 blocks. */
/*              Consider special case IB = 1, NPR = 1 and */
/*              NPR+NPC > NSUP-NLOW+1 to avoid incompatibility. */

#line 504 "SB01BD.f"
		if (ib == 1 && npr == 0 || ib == 1 && npr == 1 && nsup > nlow 
			&& npr + npc > nsup - nlow + 1) {
#line 507 "SB01BD.f"
		    if (nsup > 2) {
#line 508 "SB01BD.f"
			if (a[nsup - 1 + (nsup - 2) * a_dim1] != 0.) {

/*                       Interchange with the adjacent 2x2 block. */

/*                       Workspace needed: N. */

#line 514 "SB01BD.f"
			    i__1 = nsup - 2;
#line 514 "SB01BD.f"
			    dlaexc_(&c_true, n, &a[a_offset], lda, &z__[
				    z_offset], ldz, &i__1, &c__2, &c__1, &
				    dwork[kw], info);
#line 516 "SB01BD.f"
			    if (*info != 0) {
#line 517 "SB01BD.f"
				*info = 2;
#line 518 "SB01BD.f"
				return 0;
#line 519 "SB01BD.f"
			    }
#line 520 "SB01BD.f"
			} else {

/*                       Form a non-simple block by extending the last */
/*                       block with a 1x1 block. */

#line 525 "SB01BD.f"
			    simplb = FALSE_;
#line 526 "SB01BD.f"
			}
#line 527 "SB01BD.f"
		    } else {
#line 528 "SB01BD.f"
			simplb = FALSE_;
#line 529 "SB01BD.f"
		    }
#line 530 "SB01BD.f"
		    ib = 2;
#line 531 "SB01BD.f"
		}
#line 532 "SB01BD.f"
		nl = nsup - ib + 1;

/*              Compute G, the current last IB rows of Z'*B. */

#line 536 "SB01BD.f"
		dgemm_("Transpose", "NoTranspose", &ib, m, n, &c_b21, &z__[nl 
			* z_dim1 + 1], ldz, &b[b_offset], ldb, &c_b14, &dwork[
			kg], &ib, (ftnlen)9, (ftnlen)11);

/*              Check the controllability for the current block. */

#line 541 "SB01BD.f"
		if (dlange_("1", &ib, m, &dwork[kg], &ib, &dwork[kw], (ftnlen)
			1) <= tolerb) {

/*                 Deflate the uncontrollable block and resume the */
/*                 main loop. */

#line 547 "SB01BD.f"
		    nsup -= ib;
#line 548 "SB01BD.f"
		    *nup += ib;
#line 549 "SB01BD.f"
		    goto L20;
#line 550 "SB01BD.f"
		}

#line 552 "SB01BD.f"
		if (*nap + ib > *np) {

/*                 No sufficient eigenvalues to be assigned. */

#line 556 "SB01BD.f"
		    *info = 3;
#line 557 "SB01BD.f"
		} else {
#line 558 "SB01BD.f"
		    if (ib == 1) {

/*                    A 1-by-1 block. */

/*                    Assign the real eigenvalue nearest to A(NSUP,NSUP). */

#line 564 "SB01BD.f"
			x = a[nsup + nsup * a_dim1];
#line 565 "SB01BD.f"
			sb01bx_(&c_true, &npr, &x, &x, &wr[1], &x, &s, &p);
#line 566 "SB01BD.f"
			--npr;
#line 567 "SB01BD.f"
			ceig = FALSE_;
#line 568 "SB01BD.f"
		    } else {

/*                    A 2-by-2 block. */

#line 572 "SB01BD.f"
			if (simplb) {

/*                       Simple 2-by-2 block with complex eigenvalues. */
/*                       Compute the eigenvalues of the last block. */

#line 577 "SB01BD.f"
			    mb03qy_(n, &nl, &a[a_offset], lda, &z__[z_offset],
				     ldz, &x, &y, info);
#line 578 "SB01BD.f"
			    if (npc > 1) {
#line 579 "SB01BD.f"
				sb01bx_(&c_false, &npc, &x, &y, &wr[ipc], &wi[
					ipc], &s, &p);
#line 581 "SB01BD.f"
				npc += -2;
#line 582 "SB01BD.f"
				ceig = TRUE_;
#line 583 "SB01BD.f"
			    } else {

/*                          Choose the nearest two real eigenvalues. */

#line 587 "SB01BD.f"
				sb01bx_(&c_true, &npr, &x, &x, &wr[1], &x, &s,
					 &p);
#line 588 "SB01BD.f"
				i__1 = npr - 1;
#line 588 "SB01BD.f"
				sb01bx_(&c_true, &i__1, &x, &x, &wr[1], &x, &
					y, &p);
#line 590 "SB01BD.f"
				p = s * y;
#line 591 "SB01BD.f"
				s += y;
#line 592 "SB01BD.f"
				npr += -2;
#line 593 "SB01BD.f"
				ceig = FALSE_;
#line 594 "SB01BD.f"
			    }
#line 595 "SB01BD.f"
			} else {

/*                       Non-simple 2x2 block with real eigenvalues. */
/*                       Choose the nearest pair of complex eigenvalues. */

#line 600 "SB01BD.f"
			    x = (a[nl + nl * a_dim1] + a[nsup + nsup * a_dim1]
				    ) / 2.;
#line 601 "SB01BD.f"
			    sb01bx_(&c_false, &npc, &x, &c_b14, &wr[ipc], &wi[
				    ipc], &s, &p);
#line 603 "SB01BD.f"
			    npc += -2;
#line 604 "SB01BD.f"
			}
#line 605 "SB01BD.f"
		    }

/*                 Form the IBxIB matrix A2 from the current diagonal */
/*                 block. */

#line 610 "SB01BD.f"
		    a2[0] = a[nl + nl * a_dim1];
#line 611 "SB01BD.f"
		    if (ib > 1) {
#line 612 "SB01BD.f"
			a2[2] = a[nl + nsup * a_dim1];
#line 613 "SB01BD.f"
			a2[1] = a[nsup + nl * a_dim1];
#line 614 "SB01BD.f"
			a2[3] = a[nsup + nsup * a_dim1];
#line 615 "SB01BD.f"
		    }

/*                 Determine the M-by-IB feedback matrix FI which */
/*                 assigns the chosen IB eigenvalues for the pair (A2,G). */

/*                 Workspace needed: 5*M. */

#line 622 "SB01BD.f"
		    sb01by_(&ib, m, &s, &p, a2, &dwork[kg], &dwork[kfi], &
			    toler, &dwork[kw], &ierr);
#line 624 "SB01BD.f"
		    if (ierr != 0) {
#line 625 "SB01BD.f"
			if (ib == 1 || simplb) {

/*                       The simple 1x1 block is uncontrollable. */

#line 629 "SB01BD.f"
			    nsup -= ib;
#line 630 "SB01BD.f"
			    if (ceig) {
#line 631 "SB01BD.f"
				npc += ib;
#line 632 "SB01BD.f"
			    } else {
#line 633 "SB01BD.f"
				npr += ib;
#line 634 "SB01BD.f"
			    }
#line 635 "SB01BD.f"
			    *nup += ib;
#line 636 "SB01BD.f"
			} else {

/*                       The non-simple 2x2 block is uncontrollable. */
/*                       Eliminate its uncontrollable part by using */
/*                       the information in elements FI(1,1) and F(1,2). */

#line 642 "SB01BD.f"
			    c__ = dwork[kfi];
#line 643 "SB01BD.f"
			    s = dwork[kfi + ib];

/*                       Apply the transformation to A and accumulate it */
/*                       in Z. */

#line 648 "SB01BD.f"
			    i__1 = *n - nl + 1;
#line 648 "SB01BD.f"
			    drot_(&i__1, &a[nl + nl * a_dim1], lda, &a[nsup + 
				    nl * a_dim1], lda, &c__, &s);
#line 650 "SB01BD.f"
			    drot_(n, &a[nl * a_dim1 + 1], &c__1, &a[nsup * 
				    a_dim1 + 1], &c__1, &c__, &s);
#line 651 "SB01BD.f"
			    drot_(n, &z__[nl * z_dim1 + 1], &c__1, &z__[nsup *
				     z_dim1 + 1], &c__1, &c__, &s);

/*                       Annihilate the subdiagonal element of the last */
/*                       block, redefine the upper limit for the bottom */
/*                       block and resume the main loop. */

#line 657 "SB01BD.f"
			    a[nsup + nl * a_dim1] = 0.;
#line 658 "SB01BD.f"
			    nsup = nl;
#line 659 "SB01BD.f"
			    ++(*nup);
#line 660 "SB01BD.f"
			    npc += 2;
#line 661 "SB01BD.f"
			}
#line 662 "SB01BD.f"
		    } else {

/*                    Successful assignment of IB eigenvalues. */

/*                    Update the feedback matrix F <-- F + [0 FI]*Z'. */

#line 668 "SB01BD.f"
			dgemm_("NoTranspose", "Transpose", m, n, &ib, &c_b21, 
				&dwork[kfi], m, &z__[nl * z_dim1 + 1], ldz, &
				c_b21, &f[f_offset], ldf, (ftnlen)11, (ftnlen)
				9);

/*                    Check for possible numerical instability. */

#line 674 "SB01BD.f"
			if (dlange_("1", m, &ib, &dwork[kfi], m, &dwork[kw], (
				ftnlen)1) > rmax) {
#line 674 "SB01BD.f"
			    ++(*iwarn);
#line 674 "SB01BD.f"
			}

/*                    Update the state matrix A <-- A + Z'*B*[0 FI]. */
/*                    Workspace needed: 2*N+4*M. */

#line 680 "SB01BD.f"
			dgemm_("NoTranspose", "NoTranspose", n, &ib, m, &
				c_b21, &b[b_offset], ldb, &dwork[kfi], m, &
				c_b14, &dwork[kw], n, (ftnlen)11, (ftnlen)11);
#line 683 "SB01BD.f"
			dgemm_("Transpose", "NoTranspose", &nsup, &ib, n, &
				c_b21, &z__[z_offset], ldz, &dwork[kw], n, &
				c_b21, &a[nl * a_dim1 + 1], lda, (ftnlen)9, (
				ftnlen)11);

/*                    Try to split the 2x2 block. */

#line 689 "SB01BD.f"
			if (ib == 2) {
#line 689 "SB01BD.f"
			    mb03qy_(n, &nl, &a[a_offset], lda, &z__[z_offset],
				     ldz, &x, &y, info);
#line 689 "SB01BD.f"
			}
#line 692 "SB01BD.f"
			*nap += ib;
#line 693 "SB01BD.f"
			if (nlow + ib <= nsup) {

/*                       Move the last block(s) to the leading */
/*                       position(s) of the bottom block. */

#line 698 "SB01BD.f"
			    ncur1 = nsup - ib;
#line 699 "SB01BD.f"
			    nmoves = 1;
#line 700 "SB01BD.f"
			    if (ib == 2 && a[nsup + (nsup - 1) * a_dim1] == 
				    0.) {
#line 701 "SB01BD.f"
				ib = 1;
#line 702 "SB01BD.f"
				nmoves = 2;
#line 703 "SB01BD.f"
			    }

/*                       WHILE (NMOVES > 0) DO */
#line 706 "SB01BD.f"
L30:
#line 706 "SB01BD.f"
			    if (nmoves > 0) {
#line 707 "SB01BD.f"
				ncur = ncur1;

/*                          WHILE (NCUR >= NLOW) DO */
#line 710 "SB01BD.f"
L40:
#line 710 "SB01BD.f"
				if (ncur >= nlow) {

/*                             Loop for the last block positioning. */

#line 714 "SB01BD.f"
				    ib1 = 1;
#line 715 "SB01BD.f"
				    if (ncur > nlow) {
#line 716 "SB01BD.f"
					if (a[ncur + (ncur - 1) * a_dim1] != 
						0.) {
#line 716 "SB01BD.f"
					    ib1 = 2;
#line 716 "SB01BD.f"
					}
#line 717 "SB01BD.f"
				    }
#line 718 "SB01BD.f"
				    i__1 = ncur - ib1 + 1;
#line 718 "SB01BD.f"
				    dlaexc_(&c_true, n, &a[a_offset], lda, &
					    z__[z_offset], ldz, &i__1, &ib1, &
					    ib, &dwork[kw], info);
#line 721 "SB01BD.f"
				    if (*info != 0) {
#line 722 "SB01BD.f"
					*info = 2;
#line 723 "SB01BD.f"
					return 0;
#line 724 "SB01BD.f"
				    }
#line 725 "SB01BD.f"
				    ncur -= ib1;
#line 726 "SB01BD.f"
				    goto L40;
#line 727 "SB01BD.f"
				}

/*                          END WHILE 40 */

#line 731 "SB01BD.f"
				--nmoves;
#line 732 "SB01BD.f"
				++ncur1;
#line 733 "SB01BD.f"
				nlow += ib;
#line 734 "SB01BD.f"
				goto L30;
#line 735 "SB01BD.f"
			    }

/*                       END WHILE 30 */

#line 739 "SB01BD.f"
			} else {
#line 740 "SB01BD.f"
			    nlow += ib;
#line 741 "SB01BD.f"
			}
#line 742 "SB01BD.f"
		    }
#line 743 "SB01BD.f"
		}
#line 744 "SB01BD.f"
	    }
#line 745 "SB01BD.f"
	    if (*info == 0) {
#line 745 "SB01BD.f"
		goto L20;
#line 745 "SB01BD.f"
	    }

/*        END WHILE 20 */

#line 749 "SB01BD.f"
	}

/* Computing MAX */
#line 751 "SB01BD.f"
	i__1 = wrkopt, i__2 = *m * 5, i__1 = max(i__1,i__2), i__2 = (*n << 1) 
		+ (*m << 2);
#line 751 "SB01BD.f"
	wrkopt = max(i__1,i__2);
#line 752 "SB01BD.f"
    }

/*     Annihilate the elements below the first subdiagonal of A. */

#line 756 "SB01BD.f"
    if (*n > 2) {
#line 756 "SB01BD.f"
	i__1 = *n - 2;
#line 756 "SB01BD.f"
	i__2 = *n - 2;
#line 756 "SB01BD.f"
	dlaset_("L", &i__1, &i__2, &c_b14, &c_b14, &a[a_dim1 + 3], lda, (
		ftnlen)1);
#line 756 "SB01BD.f"
    }
#line 758 "SB01BD.f"
    if (*nap > 0) {

/*        Move the assigned eigenvalues in the first NAP positions of */
/*        WR and WI. */

#line 763 "SB01BD.f"
	k = ipc - npr - 1;
#line 764 "SB01BD.f"
	if (k > 0) {
#line 764 "SB01BD.f"
	    dswap_(&k, &wr[npr + 1], &c__1, &wr[1], &c__1);
#line 764 "SB01BD.f"
	}
#line 765 "SB01BD.f"
	j = *nap - k;
#line 766 "SB01BD.f"
	if (j > 0) {
#line 767 "SB01BD.f"
	    dswap_(&j, &wr[ipc + npc], &c__1, &wr[k + 1], &c__1);
#line 768 "SB01BD.f"
	    dswap_(&j, &wi[ipc + npc], &c__1, &wi[k + 1], &c__1);
#line 769 "SB01BD.f"
	}
#line 770 "SB01BD.f"
    }

#line 772 "SB01BD.f"
    dwork[1] = (doublereal) wrkopt;

#line 774 "SB01BD.f"
    return 0;
/* *** Last line of SB01BD *** */
} /* sb01bd_ */

