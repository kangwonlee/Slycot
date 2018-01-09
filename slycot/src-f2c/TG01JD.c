#line 1 "TG01JD.f"
/* TG01JD.f -- translated by f2c (version 20100827).
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

#line 1 "TG01JD.f"
/* Table of constant values */

static doublereal c_b12 = 0.;
static integer c__1 = 1;

/* Subroutine */ int tg01jd_(char *job, char *systyp, char *equil, integer *n,
	 integer *m, integer *p, doublereal *a, integer *lda, doublereal *e, 
	integer *lde, doublereal *b, integer *ldb, doublereal *c__, integer *
	ldc, integer *nr, integer *infred, doublereal *tol, integer *iwork, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen job_len, 
	ftnlen systyp_len, ftnlen equil_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, e_dim1, 
	    e_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer m1, n1, p1, nc, lba, lbe, ldm, ldp, ldq, kwa, kwb, kwc;
    static doublereal dum[1];
    static integer kwe, ldz;
    static char jobq[1], jobz[1];
    extern /* Subroutine */ int ma02cd_(integer *, integer *, integer *, 
	    doublereal *, integer *), tg01ad_(char *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    static logical ljobc;
    static integer nblck;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01xd_(char *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static logical ljobo;
    extern /* Subroutine */ int tg01hx_(char *, char *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen);
    static integer maxmp;
    static logical lsysp, lsysr, lsyss, lspace, fincon, infcon;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical finobs, infobs, ljobir, lequil;


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

/*     To find a reduced (controllable, observable, or irreducible) */
/*     descriptor representation (Ar-lambda*Er,Br,Cr) for an original */
/*     descriptor representation (A-lambda*E,B,C). */
/*     The pencil Ar-lambda*Er is in an upper block Hessenberg form, with */
/*     either Ar or Er upper triangular. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Indicates whether the user wishes to remove the */
/*             uncontrollable and/or unobservable parts as follows: */
/*             = 'I':  Remove both the uncontrollable and unobservable */
/*                     parts to get an irreducible descriptor */
/*                     representation; */
/*             = 'C':  Remove the uncontrollable part only to get a */
/*                     controllable descriptor representation; */
/*             = 'O':  Remove the unobservable part only to get an */
/*                     observable descriptor representation. */

/*     SYSTYP  CHARACTER*1 */
/*             Indicates the type of descriptor system algorithm */
/*             to be applied according to the assumed */
/*             transfer-function matrix as follows: */
/*             = 'R':  Rational transfer-function matrix; */
/*             = 'S':  Proper (standard) transfer-function matrix; */
/*             = 'P':  Polynomial transfer-function matrix. */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to preliminarily scale */
/*             the system (A-lambda*E,B,C) as follows: */
/*             = 'S':  Perform scaling; */
/*             = 'N':  Do not perform scaling. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The dimension of the descriptor state vector; also the */
/*             order of square matrices A and E, the number of rows of */
/*             matrix B, and the number of columns of matrix C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The dimension of descriptor system input vector; also the */
/*             number of columns of matrix B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The dimension of descriptor system output vector; also the */
/*             number of rows of matrix C.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original state matrix A. */
/*             On exit, the leading NR-by-NR part of this array contains */
/*             the reduced order state matrix Ar of an irreducible, */
/*             controllable, or observable realization for the original */
/*             system, depending on the value of JOB, JOB = 'I', */
/*             JOB = 'C', or JOB = 'O', respectively. */
/*             The matrix Ar is upper triangular if SYSTYP = 'R' or 'P'. */
/*             If SYSTYP = 'S' and JOB = 'C', the matrix [Br Ar] */
/*             is in a controllable staircase form (see TG01HD). */
/*             If SYSTYP = 'S' and JOB = 'I' or 'O', the matrix ( Ar ) */
/*                                                              ( Cr ) */
/*             is in an observable staircase form (see TG01HD). */
/*             The block structure of staircase forms is contained */
/*             in the leading INFRED(7) elements of IWORK. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original descriptor matrix E. */
/*             On exit, the leading NR-by-NR part of this array contains */
/*             the reduced order descriptor matrix Er of an irreducible, */
/*             controllable, or observable realization for the original */
/*             system, depending on the value of JOB, JOB = 'I', */
/*             JOB = 'C', or JOB = 'O', respectively. */
/*             The resulting Er has INFRED(6) nonzero sub-diagonals. */
/*             If at least for one k = 1,...,4, INFRED(k) >= 0, then the */
/*             resulting Er is structured being either upper triangular */
/*             or block Hessenberg, in accordance to the last */
/*             performed order reduction phase (see METHOD). */
/*             The block structure of staircase forms is contained */
/*             in the leading INFRED(7) elements of IWORK. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M), */
/*             if JOB = 'C', or (LDB,MAX(M,P)), otherwise. */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the original input matrix B; if JOB = 'I', */
/*             or JOB = 'O', the remainder of the leading N-by-MAX(M,P) */
/*             part is used as internal workspace. */
/*             On exit, the leading NR-by-M part of this array contains */
/*             the reduced input matrix Br of an irreducible, */
/*             controllable, or observable realization for the original */
/*             system, depending on the value of JOB, JOB = 'I', */
/*             JOB = 'C', or JOB = 'O', respectively. */
/*             If JOB = 'C', only the first IWORK(1) rows of B are */
/*             nonzero. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the original output matrix C; if JOB = 'I', */
/*             or JOB = 'O', the remainder of the leading MAX(M,P)-by-N */
/*             part is used as internal workspace. */
/*             On exit, the leading P-by-NR part of this array contains */
/*             the transformed state/output matrix Cr of an irreducible, */
/*             controllable, or observable realization for the original */
/*             system, depending on the value of JOB, JOB = 'I', */
/*             JOB = 'C', or JOB = 'O', respectively. */
/*             If JOB = 'I', or JOB = 'O', only the last IWORK(1) columns */
/*             (in the first NR columns) of C are nonzero. */

/*     LDC     INTEGER */
/*             The leading dimension of array C. */
/*             LDC >= MAX(1,M,P) if N > 0. */
/*             LDC >= 1          if N = 0. */

/*     NR      (output) INTEGER */
/*             The order of the reduced descriptor representation */
/*             (Ar-lambda*Er,Br,Cr) of an irreducible, controllable, */
/*             or observable realization for the original system, */
/*             depending on JOB = 'I', JOB = 'C', or JOB = 'O', */
/*             respectively. */

/*     INFRED  (output) INTEGER array, dimension 7 */
/*             This array contains information on performed reduction */
/*             and on structure of resulting system matrices as follows: */
/*             INFRED(k) >= 0 (k = 1, 2, 3, or 4) if Phase k of reduction */
/*                            (see METHOD) has been performed. In this */
/*                            case, INFRED(k) is the achieved order */
/*                            reduction in Phase k. */
/*             INFRED(k) < 0  (k = 1, 2, 3, or 4) if Phase k was not */
/*                            performed. */
/*             INFRED(5)  -   the number of nonzero sub-diagonals of A. */
/*             INFRED(6)  -   the number of nonzero sub-diagonals of E. */
/*             INFRED(7)  -   the number of blocks in the resulting */
/*                            staircase form at last performed reduction */
/*                            phase. The block dimensions are contained */
/*                            in the first INFRED(7) elements of IWORK. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in rank determinations when */
/*             transforming (A-lambda*E,B,C). If the user sets TOL > 0, */
/*             then the given value of TOL is used as a lower bound for */
/*             reciprocal condition numbers in rank determinations; a */
/*             (sub)matrix whose estimated condition number is less than */
/*             1/TOL is considered to be of full rank.  If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance, */
/*             defined by  TOLDEF = N*N*EPS,  is used instead, where */
/*             EPS is the machine precision (see LAPACK Library routine */
/*             DLAMCH).  TOL < 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension N+MAX(M,P) */
/*             On exit, if INFO = 0, the leading INFRED(7) elements of */
/*             IWORK contain the orders of the diagonal blocks of */
/*             Ar-lambda*Er. */

/*     DWORK   DOUBLE PRECISION array, dimension LDWORK */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(8*N,2*M,2*P), if EQUIL = 'S'; */
/*             LDWORK >= MAX(N,2*M,2*P),   if EQUIL = 'N'. */
/*             If LDWORK >= MAX(2*N*N+N*M+N*P)+MAX(N,2*M,2*P) then more */
/*             accurate results are to be expected by performing only */
/*             those reductions phases (see METHOD), where effective */
/*             order reduction occurs. This is achieved by saving the */
/*             system matrices before each phase and restoring them if no */
/*             order reduction took place. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The subroutine is based on the reduction algorithms of [1]. */
/*     The order reduction is performed in 4 phases: */
/*     Phase 1: Eliminate all finite uncontrolable eigenvalues. */
/*              The resulting matrix ( Br Ar ) is in a controllable */
/*              staircase form (see SLICOT Library routine TG01HD), and */
/*              Er is upper triangular. */
/*              This phase is performed if JOB = 'I' or 'C' and */
/*              SYSTYP = 'R' or 'S'. */
/*     Phase 2: Eliminate all infinite and finite nonzero uncontrollable */
/*              eigenvalues. The resulting matrix ( Br Er ) is in a */
/*              controllable staircase form (see TG01HD), and Ar is */
/*              upper triangular. */
/*              This phase is performed if JOB = 'I' or 'C' and */
/*              SYSTYP = 'R' or 'P'. */
/*     Phase 3: Eliminate all finite unobservable eigenvalues. */
/*              The resulting matrix ( Ar ) is in an observable */
/*                                   ( Cr ) */
/*              staircase form (see SLICOT Library routine TG01ID), and */
/*              Er is upper triangular. */
/*              This phase is performed if JOB = 'I' or 'O' and */
/*              SYSTYP = 'R' or 'S'. */
/*     Phase 4: Eliminate all infinite and finite nonzero unobservable */
/*              eigenvalues. The resulting matrix ( Er ) is in an */
/*                                                ( Cr ) */
/*              observable staircase form (see TG01ID), and Ar is */
/*              upper triangular. */
/*              This phase is performed if JOB = 'I' or 'O' and */
/*              SYSTYP = 'R' or 'P'. */

/*     REFERENCES */

/*     [1] A. Varga */
/*         Computation of Irreducible Generalized State-Space */
/*         Realizations. */
/*         Kybernetika, vol. 26, pp. 89-106, 1990. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically backward stable and requires */
/*     0( N**3 )  floating point operations. */

/*     FURTHER COMMENTS */

/*     If the pencil (A-lambda*E) has no zero eigenvalues, then an */
/*     irreducible realization can be computed skipping Phases 1 and 3 */
/*     by using the setting: JOB = 'I' and SYSTYP = 'P'. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     April 1999. Based on the RASP routine RPDSIR. */

/*     REVISIONS */

/*     July 1999, V. Sima, Research Institute for Informatics, Bucharest. */
/*     May 2003, A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     May 2003, March 2004, V. Sima. */

/*     KEYWORDS */

/*     Controllability, irreducible realization, observability, */
/*     orthogonal canonical form, orthogonal transformation. */

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

#line 307 "TG01JD.f"
    /* Parameter adjustments */
#line 307 "TG01JD.f"
    a_dim1 = *lda;
#line 307 "TG01JD.f"
    a_offset = 1 + a_dim1;
#line 307 "TG01JD.f"
    a -= a_offset;
#line 307 "TG01JD.f"
    e_dim1 = *lde;
#line 307 "TG01JD.f"
    e_offset = 1 + e_dim1;
#line 307 "TG01JD.f"
    e -= e_offset;
#line 307 "TG01JD.f"
    b_dim1 = *ldb;
#line 307 "TG01JD.f"
    b_offset = 1 + b_dim1;
#line 307 "TG01JD.f"
    b -= b_offset;
#line 307 "TG01JD.f"
    c_dim1 = *ldc;
#line 307 "TG01JD.f"
    c_offset = 1 + c_dim1;
#line 307 "TG01JD.f"
    c__ -= c_offset;
#line 307 "TG01JD.f"
    --infred;
#line 307 "TG01JD.f"
    --iwork;
#line 307 "TG01JD.f"
    --dwork;
#line 307 "TG01JD.f"

#line 307 "TG01JD.f"
    /* Function Body */
#line 307 "TG01JD.f"
    *info = 0;
#line 308 "TG01JD.f"
    maxmp = max(*m,*p);
#line 309 "TG01JD.f"
    n1 = max(1,*n);

/*     Decode JOB. */

#line 313 "TG01JD.f"
    ljobir = lsame_(job, "I", (ftnlen)1, (ftnlen)1);
#line 314 "TG01JD.f"
    ljobc = ljobir || lsame_(job, "C", (ftnlen)1, (ftnlen)1);
#line 315 "TG01JD.f"
    ljobo = ljobir || lsame_(job, "O", (ftnlen)1, (ftnlen)1);

/*     Decode SYSTYP. */

#line 319 "TG01JD.f"
    lsysr = lsame_(systyp, "R", (ftnlen)1, (ftnlen)1);
#line 320 "TG01JD.f"
    lsyss = lsysr || lsame_(systyp, "S", (ftnlen)1, (ftnlen)1);
#line 321 "TG01JD.f"
    lsysp = lsysr || lsame_(systyp, "P", (ftnlen)1, (ftnlen)1);

#line 323 "TG01JD.f"
    lequil = lsame_(equil, "S", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 327 "TG01JD.f"
    if (! ljobc && ! ljobo) {
#line 328 "TG01JD.f"
	*info = -1;
#line 329 "TG01JD.f"
    } else if (! lsyss && ! lsysp) {
#line 330 "TG01JD.f"
	*info = -2;
#line 331 "TG01JD.f"
    } else if (! lequil && ! lsame_(equil, "N", (ftnlen)1, (ftnlen)1)) {
#line 332 "TG01JD.f"
	*info = -3;
#line 333 "TG01JD.f"
    } else if (*n < 0) {
#line 334 "TG01JD.f"
	*info = -4;
#line 335 "TG01JD.f"
    } else if (*m < 0) {
#line 336 "TG01JD.f"
	*info = -5;
#line 337 "TG01JD.f"
    } else if (*p < 0) {
#line 338 "TG01JD.f"
	*info = -6;
#line 339 "TG01JD.f"
    } else if (*lda < n1) {
#line 340 "TG01JD.f"
	*info = -8;
#line 341 "TG01JD.f"
    } else if (*lde < n1) {
#line 342 "TG01JD.f"
	*info = -10;
#line 343 "TG01JD.f"
    } else if (*ldb < n1) {
#line 344 "TG01JD.f"
	*info = -12;
#line 345 "TG01JD.f"
    } else if (*ldc < 1 || *n > 0 && *ldc < maxmp) {
#line 346 "TG01JD.f"
	*info = -14;
#line 347 "TG01JD.f"
    } else if (*tol >= 1.) {
#line 348 "TG01JD.f"
	*info = -17;
#line 349 "TG01JD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 349 "TG01JD.f"
	i__1 = *n, i__2 = maxmp << 1;
/* Computing MAX */
#line 349 "TG01JD.f"
	i__3 = *n << 3, i__4 = maxmp << 1;
#line 349 "TG01JD.f"
	if (! lequil && *ldwork < max(i__1,i__2) || lequil && *ldwork < max(
		i__3,i__4)) {
#line 351 "TG01JD.f"
	    *info = -20;
#line 352 "TG01JD.f"
	}
#line 352 "TG01JD.f"
    }

#line 354 "TG01JD.f"
    if (*info != 0) {

/*        Error return. */

#line 358 "TG01JD.f"
	i__1 = -(*info);
#line 358 "TG01JD.f"
	xerbla_("TG01JD", &i__1, (ftnlen)6);
#line 359 "TG01JD.f"
	return 0;
#line 360 "TG01JD.f"
    }

/*     Quick return if possible. */

#line 364 "TG01JD.f"
    infred[1] = -1;
#line 365 "TG01JD.f"
    infred[2] = -1;
#line 366 "TG01JD.f"
    infred[3] = -1;
#line 367 "TG01JD.f"
    infred[4] = -1;
#line 368 "TG01JD.f"
    infred[5] = 0;
#line 369 "TG01JD.f"
    infred[6] = 0;
#line 370 "TG01JD.f"
    infred[7] = 0;

#line 372 "TG01JD.f"
    if (max(*n,maxmp) == 0) {
#line 373 "TG01JD.f"
	*nr = 0;
#line 374 "TG01JD.f"
	return 0;
#line 375 "TG01JD.f"
    }

#line 377 "TG01JD.f"
    m1 = max(1,*m);
#line 378 "TG01JD.f"
    p1 = max(1,*p);
#line 379 "TG01JD.f"
    ldm = max(*ldc,*m);
#line 380 "TG01JD.f"
    ldp = max(*ldc,*p);

/*     Set controllability/observability determination options. */

#line 384 "TG01JD.f"
    fincon = ljobc && lsyss;
#line 385 "TG01JD.f"
    infcon = ljobc && lsysp;
#line 386 "TG01JD.f"
    finobs = ljobo && lsyss;
#line 387 "TG01JD.f"
    infobs = ljobo && lsysp;

/*     Set large workspace option and determine offsets. */

/* Computing MAX */
#line 391 "TG01JD.f"
    i__1 = *n, i__2 = maxmp << 1;
#line 391 "TG01JD.f"
    lspace = *ldwork >= *n * ((*n << 1) + *m + *p) + max(i__1,i__2);
/* Computing MAX */
#line 392 "TG01JD.f"
    i__1 = *n, i__2 = maxmp << 1;
#line 392 "TG01JD.f"
    kwa = max(i__1,i__2) + 1;
#line 393 "TG01JD.f"
    kwe = kwa + *n * *n;
#line 394 "TG01JD.f"
    kwb = kwe + *n * *n;
#line 395 "TG01JD.f"
    kwc = kwb + *n * *m;

/*     If required, scale the system (A-lambda*E,B,C). */
/*     Workspace: need 8*N. */

#line 400 "TG01JD.f"
    if (lequil) {
#line 401 "TG01JD.f"
	tg01ad_("All", n, n, m, p, &c_b12, &a[a_offset], lda, &e[e_offset], 
		lde, &b[b_offset], ldb, &c__[c_offset], &ldp, &dwork[1], &
		dwork[*n + 1], &dwork[(*n << 1) + 1], info, (ftnlen)3);
#line 403 "TG01JD.f"
    }

#line 405 "TG01JD.f"
    *(unsigned char *)jobq = 'N';
#line 406 "TG01JD.f"
    *(unsigned char *)jobz = 'N';
#line 407 "TG01JD.f"
    ldq = 1;
#line 408 "TG01JD.f"
    ldz = 1;
/* Computing MAX */
#line 409 "TG01JD.f"
    i__1 = 0, i__2 = *n - 1;
#line 409 "TG01JD.f"
    lba = max(i__1,i__2);
#line 410 "TG01JD.f"
    lbe = lba;
#line 411 "TG01JD.f"
    nc = *n;
#line 412 "TG01JD.f"
    *nr = *n;

#line 414 "TG01JD.f"
    if (fincon) {

/*        Phase 1: Eliminate all finite uncontrolable eigenvalues. */

#line 418 "TG01JD.f"
	if (lspace) {

/*           Save system matrices. */

#line 422 "TG01JD.f"
	    dlacpy_("Full", &nc, &nc, &a[a_offset], lda, &dwork[kwa], &n1, (
		    ftnlen)4);
#line 423 "TG01JD.f"
	    dlacpy_("Full", &nc, &nc, &e[e_offset], lde, &dwork[kwe], &n1, (
		    ftnlen)4);
#line 424 "TG01JD.f"
	    dlacpy_("Full", &nc, m, &b[b_offset], ldb, &dwork[kwb], &n1, (
		    ftnlen)4);
#line 425 "TG01JD.f"
	    dlacpy_("Full", p, &nc, &c__[c_offset], ldc, &dwork[kwc], &p1, (
		    ftnlen)4);
#line 426 "TG01JD.f"
	}

/*        Perform finite controllability form reduction. */
/*        Workspace: need   MAX(N,2*M). */

#line 431 "TG01JD.f"
	tg01hx_(jobq, jobz, &nc, &nc, m, p, &nc, &lbe, &a[a_offset], lda, &e[
		e_offset], lde, &b[b_offset], ldb, &c__[c_offset], &ldp, dum, 
		&ldq, dum, &ldz, nr, &nblck, &iwork[1], tol, &iwork[*n + 1], &
		dwork[1], info, (ftnlen)1, (ftnlen)1);
#line 434 "TG01JD.f"
	if (*nr < nc || ! lspace) {
#line 435 "TG01JD.f"
	    if (nblck > 1) {
#line 436 "TG01JD.f"
		lba = iwork[1] + iwork[2] - 1;
#line 437 "TG01JD.f"
	    } else if (nblck == 1) {
#line 438 "TG01JD.f"
		lba = iwork[1] - 1;
#line 439 "TG01JD.f"
	    } else {
#line 440 "TG01JD.f"
		lba = 0;
#line 441 "TG01JD.f"
	    }
#line 442 "TG01JD.f"
	    lbe = 0;
#line 443 "TG01JD.f"
	    infred[1] = nc - *nr;
#line 444 "TG01JD.f"
	    infred[7] = nblck;
#line 445 "TG01JD.f"
	    nc = *nr;
#line 446 "TG01JD.f"
	} else {

/*           Restore system matrices. */

#line 450 "TG01JD.f"
	    dlacpy_("Full", &nc, &nc, &dwork[kwa], &n1, &a[a_offset], lda, (
		    ftnlen)4);
#line 451 "TG01JD.f"
	    dlacpy_("Full", &nc, &nc, &dwork[kwe], &n1, &e[e_offset], lde, (
		    ftnlen)4);
#line 452 "TG01JD.f"
	    dlacpy_("Full", &nc, m, &dwork[kwb], &n1, &b[b_offset], ldb, (
		    ftnlen)4);
#line 453 "TG01JD.f"
	    dlacpy_("Full", p, &nc, &dwork[kwc], &p1, &c__[c_offset], ldc, (
		    ftnlen)4);
#line 454 "TG01JD.f"
	}
#line 455 "TG01JD.f"
    }

#line 457 "TG01JD.f"
    if (infcon) {

/*        Phase 2: Eliminate all infinite and all finite nonzero */
/*                 uncontrolable eigenvalues. */

#line 462 "TG01JD.f"
	if (lspace) {

/*           Save system matrices. */

#line 466 "TG01JD.f"
	    dlacpy_("Full", &nc, &nc, &a[a_offset], lda, &dwork[kwa], &n1, (
		    ftnlen)4);
#line 467 "TG01JD.f"
	    dlacpy_("Full", &nc, &nc, &e[e_offset], lde, &dwork[kwe], &n1, (
		    ftnlen)4);
#line 468 "TG01JD.f"
	    dlacpy_("Full", &nc, m, &b[b_offset], ldb, &dwork[kwb], &n1, (
		    ftnlen)4);
#line 469 "TG01JD.f"
	    dlacpy_("Full", p, &nc, &c__[c_offset], ldc, &dwork[kwc], &p1, (
		    ftnlen)4);
#line 470 "TG01JD.f"
	}

/*        Perform infinite controllability form reduction. */
/*        Workspace: need   MAX(N,2*M). */

#line 475 "TG01JD.f"
	tg01hx_(jobq, jobz, &nc, &nc, m, p, &nc, &lba, &e[e_offset], lde, &a[
		a_offset], lda, &b[b_offset], ldb, &c__[c_offset], &ldp, dum, 
		&ldq, dum, &ldz, nr, &nblck, &iwork[1], tol, &iwork[*n + 1], &
		dwork[1], info, (ftnlen)1, (ftnlen)1);
#line 478 "TG01JD.f"
	if (*nr < nc || ! lspace) {
#line 479 "TG01JD.f"
	    if (nblck > 1) {
#line 480 "TG01JD.f"
		lbe = iwork[1] + iwork[2] - 1;
#line 481 "TG01JD.f"
	    } else if (nblck == 1) {
#line 482 "TG01JD.f"
		lbe = iwork[1] - 1;
#line 483 "TG01JD.f"
	    } else {
#line 484 "TG01JD.f"
		lbe = 0;
#line 485 "TG01JD.f"
	    }
#line 486 "TG01JD.f"
	    lba = 0;
#line 487 "TG01JD.f"
	    infred[2] = nc - *nr;
#line 488 "TG01JD.f"
	    infred[7] = nblck;
#line 489 "TG01JD.f"
	    nc = *nr;
#line 490 "TG01JD.f"
	} else {

/*           Restore system matrices. */

#line 494 "TG01JD.f"
	    dlacpy_("Full", &nc, &nc, &dwork[kwa], &n1, &a[a_offset], lda, (
		    ftnlen)4);
#line 495 "TG01JD.f"
	    dlacpy_("Full", &nc, &nc, &dwork[kwe], &n1, &e[e_offset], lde, (
		    ftnlen)4);
#line 496 "TG01JD.f"
	    dlacpy_("Full", &nc, m, &dwork[kwb], &n1, &b[b_offset], ldb, (
		    ftnlen)4);
#line 497 "TG01JD.f"
	    dlacpy_("Full", p, &nc, &dwork[kwc], &p1, &c__[c_offset], ldc, (
		    ftnlen)4);
#line 498 "TG01JD.f"
	}
#line 499 "TG01JD.f"
    }

#line 501 "TG01JD.f"
    if (finobs || infobs) {

/*        Compute the pertransposed dual system exploiting matrix shapes. */

/* Computing MAX */
#line 505 "TG01JD.f"
	i__2 = 0, i__3 = nc - 1;
#line 505 "TG01JD.f"
	i__1 = max(i__2,i__3);
#line 505 "TG01JD.f"
	tb01xd_("Z", &nc, m, p, &lba, &i__1, &a[a_offset], lda, &b[b_offset], 
		ldb, &c__[c_offset], ldc, dum, &c__1, info, (ftnlen)1);
/* Computing MAX */
#line 507 "TG01JD.f"
	i__2 = 0, i__3 = nc - 1;
#line 507 "TG01JD.f"
	i__1 = max(i__2,i__3);
#line 507 "TG01JD.f"
	ma02cd_(&nc, &lbe, &i__1, &e[e_offset], lde);
#line 508 "TG01JD.f"
    }

#line 510 "TG01JD.f"
    if (finobs) {

/*        Phase 3: Eliminate all finite unobservable eigenvalues. */

#line 514 "TG01JD.f"
	if (lspace) {

/*           Save system matrices. */

#line 518 "TG01JD.f"
	    dlacpy_("Full", &nc, &nc, &a[a_offset], lda, &dwork[kwa], &n1, (
		    ftnlen)4);
#line 519 "TG01JD.f"
	    dlacpy_("Full", &nc, &nc, &e[e_offset], lde, &dwork[kwe], &n1, (
		    ftnlen)4);
#line 520 "TG01JD.f"
	    dlacpy_("Full", &nc, p, &b[b_offset], ldb, &dwork[kwc], &n1, (
		    ftnlen)4);
#line 521 "TG01JD.f"
	    dlacpy_("Full", m, &nc, &c__[c_offset], ldc, &dwork[kwb], &m1, (
		    ftnlen)4);
#line 522 "TG01JD.f"
	}

/*        Perform finite observability form reduction. */
/*        Workspace: need   MAX(N,2*P). */

#line 527 "TG01JD.f"
	tg01hx_(jobz, jobq, &nc, &nc, p, m, &nc, &lbe, &a[a_offset], lda, &e[
		e_offset], lde, &b[b_offset], ldb, &c__[c_offset], &ldm, dum, 
		&ldz, dum, &ldq, nr, &nblck, &iwork[1], tol, &iwork[*n + 1], &
		dwork[1], info, (ftnlen)1, (ftnlen)1);
#line 530 "TG01JD.f"
	if (*nr < nc || ! lspace) {
#line 531 "TG01JD.f"
	    if (nblck > 1) {
#line 532 "TG01JD.f"
		lba = iwork[1] + iwork[2] - 1;
#line 533 "TG01JD.f"
	    } else if (nblck == 1) {
#line 534 "TG01JD.f"
		lba = iwork[1] - 1;
#line 535 "TG01JD.f"
	    } else {
#line 536 "TG01JD.f"
		lba = 0;
#line 537 "TG01JD.f"
	    }
#line 538 "TG01JD.f"
	    lbe = 0;
#line 539 "TG01JD.f"
	    infred[3] = nc - *nr;
#line 540 "TG01JD.f"
	    infred[7] = nblck;
#line 541 "TG01JD.f"
	    nc = *nr;
#line 542 "TG01JD.f"
	} else {

/*           Restore system matrices. */

#line 546 "TG01JD.f"
	    dlacpy_("Full", &nc, &nc, &dwork[kwa], &n1, &a[a_offset], lda, (
		    ftnlen)4);
#line 547 "TG01JD.f"
	    dlacpy_("Full", &nc, &nc, &dwork[kwe], &n1, &e[e_offset], lde, (
		    ftnlen)4);
#line 548 "TG01JD.f"
	    dlacpy_("Full", &nc, p, &dwork[kwc], &n1, &b[b_offset], ldb, (
		    ftnlen)4);
#line 549 "TG01JD.f"
	    dlacpy_("Full", m, &nc, &dwork[kwb], &m1, &c__[c_offset], ldc, (
		    ftnlen)4);
#line 550 "TG01JD.f"
	}
#line 551 "TG01JD.f"
    }

#line 553 "TG01JD.f"
    if (infobs) {

/*        Phase 4: Eliminate all infinite and all finite nonzero */
/*                 unobservable eigenvalues. */

#line 558 "TG01JD.f"
	if (lspace) {

/*           Save system matrices. */

#line 562 "TG01JD.f"
	    dlacpy_("Full", &nc, &nc, &a[a_offset], lda, &dwork[kwa], &n1, (
		    ftnlen)4);
#line 563 "TG01JD.f"
	    dlacpy_("Full", &nc, &nc, &e[e_offset], lde, &dwork[kwe], &n1, (
		    ftnlen)4);
#line 564 "TG01JD.f"
	    dlacpy_("Full", &nc, p, &b[b_offset], ldb, &dwork[kwc], &n1, (
		    ftnlen)4);
#line 565 "TG01JD.f"
	    dlacpy_("Full", m, &nc, &c__[c_offset], ldc, &dwork[kwb], &m1, (
		    ftnlen)4);
#line 566 "TG01JD.f"
	}

/*        Perform infinite observability form reduction. */
/*        Workspace: need   MAX(N,2*P). */

#line 571 "TG01JD.f"
	tg01hx_(jobz, jobq, &nc, &nc, p, m, &nc, &lba, &e[e_offset], lde, &a[
		a_offset], lda, &b[b_offset], ldb, &c__[c_offset], &ldm, dum, 
		&ldz, dum, &ldq, nr, &nblck, &iwork[1], tol, &iwork[*n + 1], &
		dwork[1], info, (ftnlen)1, (ftnlen)1);
#line 574 "TG01JD.f"
	if (*nr < nc || ! lspace) {
#line 575 "TG01JD.f"
	    if (nblck > 1) {
#line 576 "TG01JD.f"
		lbe = iwork[1] + iwork[2] - 1;
#line 577 "TG01JD.f"
	    } else if (nblck == 1) {
#line 578 "TG01JD.f"
		lbe = iwork[1] - 1;
#line 579 "TG01JD.f"
	    } else {
#line 580 "TG01JD.f"
		lbe = 0;
#line 581 "TG01JD.f"
	    }
#line 582 "TG01JD.f"
	    lba = 0;
#line 583 "TG01JD.f"
	    infred[4] = nc - *nr;
#line 584 "TG01JD.f"
	    infred[7] = nblck;
#line 585 "TG01JD.f"
	    nc = *nr;
#line 586 "TG01JD.f"
	} else {

/*           Restore system matrices. */

#line 590 "TG01JD.f"
	    dlacpy_("Full", &nc, &nc, &dwork[kwa], &n1, &a[a_offset], lda, (
		    ftnlen)4);
#line 591 "TG01JD.f"
	    dlacpy_("Full", &nc, &nc, &dwork[kwe], &n1, &e[e_offset], lde, (
		    ftnlen)4);
#line 592 "TG01JD.f"
	    dlacpy_("Full", &nc, p, &dwork[kwc], &n1, &b[b_offset], ldb, (
		    ftnlen)4);
#line 593 "TG01JD.f"
	    dlacpy_("Full", m, &nc, &dwork[kwb], &m1, &c__[c_offset], ldc, (
		    ftnlen)4);
#line 594 "TG01JD.f"
	}
#line 595 "TG01JD.f"
    }

#line 597 "TG01JD.f"
    if (finobs || infobs) {

/*        Compute the pertransposed dual system exploiting matrix shapes. */

/* Computing MAX */
#line 601 "TG01JD.f"
	i__2 = 0, i__3 = nc - 1;
#line 601 "TG01JD.f"
	i__1 = max(i__2,i__3);
#line 601 "TG01JD.f"
	tb01xd_("Z", &nc, p, m, &lba, &i__1, &a[a_offset], lda, &b[b_offset], 
		ldb, &c__[c_offset], ldc, dum, &c__1, info, (ftnlen)1);
/* Computing MAX */
#line 603 "TG01JD.f"
	i__2 = 0, i__3 = nc - 1;
#line 603 "TG01JD.f"
	i__1 = max(i__2,i__3);
#line 603 "TG01JD.f"
	ma02cd_(&nc, &lbe, &i__1, &e[e_offset], lde);
#line 604 "TG01JD.f"
    }

/*     Set structural information on A and E. */

#line 608 "TG01JD.f"
    infred[5] = lba;
#line 609 "TG01JD.f"
    infred[6] = lbe;

#line 611 "TG01JD.f"
    return 0;
/* *** Last line of TG01JD *** */
} /* tg01jd_ */

