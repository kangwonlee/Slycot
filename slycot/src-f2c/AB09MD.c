#line 1 "AB09MD.f"
/* AB09MD.f -- translated by f2c (version 20100827).
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

#line 1 "AB09MD.f"
/* Subroutine */ int ab09md_(char *dico, char *job, char *equil, char *ordsel,
	 integer *n, integer *m, integer *p, integer *nr, doublereal *alpha, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, integer *ns, doublereal *hsv, doublereal *tol, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *iwarn, 
	integer *info, ftnlen dico_len, ftnlen job_len, ftnlen equil_len, 
	ftnlen ordsel_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer nn, kt, ku, kw, nu, nu1, nra, kti, kwi, kwr, lwr, ierr;
    extern /* Subroutine */ int tb01id_(char *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    ab09ax_(char *, char *, char *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, ftnlen, ftnlen, ftnlen), tb01kd_(
	    char *, char *, char *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal maxred;
    static logical fixord;
    static integer iwarnl;
    static doublereal alpwrk;
    static integer wrkopt;


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

/*     To compute a reduced order model (Ar,Br,Cr) for an original */
/*     state-space representation (A,B,C) by using either the square-root */
/*     or the balancing-free square-root Balance & Truncate (B & T) */
/*     model reduction method for the ALPHA-stable part of the system. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the original system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     JOB     CHARACTER*1 */
/*             Specifies the model reduction approach to be used */
/*             as follows: */
/*             = 'B':  use the square-root Balance & Truncate method; */
/*             = 'N':  use the balancing-free square-root */
/*                     Balance & Truncate method. */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to preliminarily */
/*             equilibrate the triplet (A,B,C) as follows: */
/*             = 'S':  perform equilibration (scaling); */
/*             = 'N':  do not perform equilibration. */

/*     ORDSEL  CHARACTER*1 */
/*             Specifies the order selection method as follows: */
/*             = 'F':  the resulting order NR is fixed; */
/*             = 'A':  the resulting order NR is automatically determined */
/*                     on basis of the given tolerance TOL. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the original state-space representation, i.e. */
/*             the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     NR      (input/output) INTEGER */
/*             On entry with ORDSEL = 'F', NR is the desired order of the */
/*             resulting reduced order system.  0 <= NR <= N. */
/*             On exit, if INFO = 0, NR is the order of the resulting */
/*             reduced order model. For a system with NU ALPHA-unstable */
/*             eigenvalues and NS ALPHA-stable eigenvalues (NU+NS = N), */
/*             NR is set as follows: if ORDSEL = 'F', NR is equal to */
/*             NU+MIN(MAX(0,NR-NU),NMIN), where NR is the desired order */
/*             on entry, and NMIN is the order of a minimal realization */
/*             of the ALPHA-stable part of the given system; NMIN is */
/*             determined as the number of Hankel singular values greater */
/*             than NS*EPS*HNORM(As,Bs,Cs), where EPS is the machine */
/*             precision (see LAPACK Library Routine DLAMCH) and */
/*             HNORM(As,Bs,Cs) is the Hankel norm of the ALPHA-stable */
/*             part of the given system (computed in HSV(1)); */
/*             if ORDSEL = 'A', NR is the sum of NU and the number of */
/*             Hankel singular values greater than */
/*             MAX(TOL,NS*EPS*HNORM(As,Bs,Cs)). */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             Specifies the ALPHA-stability boundary for the eigenvalues */
/*             of the state dynamics matrix A. For a continuous-time */
/*             system (DICO = 'C'), ALPHA <= 0 is the boundary value for */
/*             the real parts of eigenvalues, while for a discrete-time */
/*             system (DICO = 'D'), 0 <= ALPHA <= 1 represents the */
/*             boundary value for the moduli of eigenvalues. */
/*             The ALPHA-stability domain does not include the boundary. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix A. */
/*             On exit, if INFO = 0, the leading NR-by-NR part of this */
/*             array contains the state dynamics matrix Ar of the reduced */
/*             order system. */
/*             The resulting A has a block-diagonal form with two blocks. */
/*             For a system with NU ALPHA-unstable eigenvalues and */
/*             NS ALPHA-stable eigenvalues (NU+NS = N), the leading */
/*             NU-by-NU block contains the unreduced part of A */
/*             corresponding to ALPHA-unstable eigenvalues in an */
/*             upper real Schur form. */
/*             The trailing (NR+NS-N)-by-(NR+NS-N) block contains */
/*             the reduced part of A corresponding to ALPHA-stable */
/*             eigenvalues. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the original input/state matrix B. */
/*             On exit, if INFO = 0, the leading NR-by-M part of this */
/*             array contains the input/state matrix Br of the reduced */
/*             order system. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the original state/output matrix C. */
/*             On exit, if INFO = 0, the leading P-by-NR part of this */
/*             array contains the state/output matrix Cr of the reduced */
/*             order system. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     NS      (output) INTEGER */
/*             The dimension of the ALPHA-stable subsystem. */

/*     HSV     (output) DOUBLE PRECISION array, dimension (N) */
/*             If INFO = 0, the leading NS elements of HSV contain the */
/*             Hankel singular values of the ALPHA-stable part of the */
/*             original system ordered decreasingly. */
/*             HSV(1) is the Hankel norm of the ALPHA-stable subsystem. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If ORDSEL = 'A', TOL contains the tolerance for */
/*             determining the order of reduced system. */
/*             For model reduction, the recommended value is */
/*             TOL = c*HNORM(As,Bs,Cs), where c is a constant in the */
/*             interval [0.00001,0.001], and HNORM(As,Bs,Cs) is the */
/*             Hankel-norm of the ALPHA-stable part of the given system */
/*             (computed in HSV(1)). */
/*             If TOL <= 0 on entry, the used default value is */
/*             TOL = NS*EPS*HNORM(As,Bs,Cs), where NS is the number of */
/*             ALPHA-stable eigenvalues of A and EPS is the machine */
/*             precision (see LAPACK Library Routine DLAMCH). */
/*             This value is appropriate to compute a minimal realization */
/*             of the ALPHA-stable part. */
/*             If ORDSEL = 'F', the value of TOL is ignored. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK = 0, if JOB = 'B'; */
/*             LIWORK = N, if JOB = 'N'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,N*(2*N+MAX(N,M,P)+5) + N*(N+1)/2). */
/*             For optimum performance LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  with ORDSEL = 'F', the selected order NR is greater */
/*                   than NSMIN, the sum of the order of the */
/*                   ALPHA-unstable part and the order of a minimal */
/*                   realization of the ALPHA-stable part of the given */
/*                   system. In this case, the resulting NR is set equal */
/*                   to NSMIN. */
/*             = 2:  with ORDSEL = 'F', the selected order NR is less */
/*                   than the order of the ALPHA-unstable part of the */
/*                   given system. In this case NR is set equal to the */
/*                   order of the ALPHA-unstable part. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the computation of the ordered real Schur form of A */
/*                   failed; */
/*             = 2:  the separation of the ALPHA-stable/unstable diagonal */
/*                   blocks failed because of very close eigenvalues; */
/*             = 3:  the computation of Hankel singular values failed. */

/*     METHOD */

/*     Let be the following linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t)                               (1) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system. The subroutine AB09MD determines for */
/*     the given system (1), the matrices of a reduced order system */

/*          d[z(t)] = Ar*z(t) + Br*u(t) */
/*          yr(t)   = Cr*z(t)                             (2) */

/*     such that */

/*     HSV(NR+NS-N) <= INFNORM(G-Gr) <= 2*[HSV(NR+NS-N+1)+...+HSV(NS)], */

/*     where G and Gr are transfer-function matrices of the systems */
/*     (A,B,C) and (Ar,Br,Cr), respectively, and INFNORM(G) is the */
/*     infinity-norm of G. */

/*     The following procedure is used to reduce a given G: */

/*     1) Decompose additively G as */

/*          G = G1 + G2 */

/*        such that G1 = (As,Bs,Cs) has only ALPHA-stable poles and */
/*        G2 = (Au,Bu,Cu) has only ALPHA-unstable poles. */

/*     2) Determine G1r, a reduced order approximation of the */
/*        ALPHA-stable part G1. */

/*     3) Assemble the reduced model Gr as */

/*           Gr = G1r + G2. */

/*     To reduce the ALPHA-stable part G1, if JOB = 'B', the square-root */
/*     Balance & Truncate method of [1] is used, and for an ALPHA-stable */
/*     continuous-time system (DICO = 'C'), the resulting reduced model */
/*     is balanced. For ALPHA-stable systems, setting TOL < 0, the */
/*     routine can be used to compute balanced minimal state-space */
/*     realizations. */

/*     If JOB = 'N', the balancing-free square-root version of the */
/*     Balance & Truncate method [2] is used to reduce the ALPHA-stable */
/*     part G1. */

/*     REFERENCES */

/*     [1] Tombs M.S. and Postlethwaite I. */
/*         Truncated balanced realization of stable, non-minimal */
/*         state-space systems. */
/*         Int. J. Control, Vol. 46, pp. 1319-1330, 1987. */

/*     [2] Varga A. */
/*         Efficient minimal realization procedure based on balancing. */
/*         Proc. of IMACS/IFAC Symp. MCTS, Lille, France, May 1991, */
/*         A. El Moudui, P. Borne, S. G. Tzafestas (Eds.), */
/*         Vol. 2, pp. 42-46. */

/*     NUMERICAL ASPECTS */

/*     The implemented methods rely on accuracy enhancing square-root or */
/*     balancing-free square-root techniques. */
/*                                         3 */
/*     The algorithms require less than 30N  floating point operations. */

/*     CONTRIBUTOR */

/*     C. Oara, University "Politehnica" Bucharest. */
/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     February 1999. Based on the RASP routines SADSDC, SRBT and SRBFT. */

/*     REVISIONS */

/*     Mar. 1999, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Nov. 2000, A. Varga, DLR Oberpfaffenhofen. */

/*     KEYWORDS */

/*     Balancing, minimal realization, model reduction, multivariable */
/*     system, state-space model, state-space representation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 319 "AB09MD.f"
    /* Parameter adjustments */
#line 319 "AB09MD.f"
    a_dim1 = *lda;
#line 319 "AB09MD.f"
    a_offset = 1 + a_dim1;
#line 319 "AB09MD.f"
    a -= a_offset;
#line 319 "AB09MD.f"
    b_dim1 = *ldb;
#line 319 "AB09MD.f"
    b_offset = 1 + b_dim1;
#line 319 "AB09MD.f"
    b -= b_offset;
#line 319 "AB09MD.f"
    c_dim1 = *ldc;
#line 319 "AB09MD.f"
    c_offset = 1 + c_dim1;
#line 319 "AB09MD.f"
    c__ -= c_offset;
#line 319 "AB09MD.f"
    --hsv;
#line 319 "AB09MD.f"
    --iwork;
#line 319 "AB09MD.f"
    --dwork;
#line 319 "AB09MD.f"

#line 319 "AB09MD.f"
    /* Function Body */
#line 319 "AB09MD.f"
    *info = 0;
#line 320 "AB09MD.f"
    *iwarn = 0;
#line 321 "AB09MD.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 322 "AB09MD.f"
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 326 "AB09MD.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 327 "AB09MD.f"
	*info = -1;
#line 328 "AB09MD.f"
    } else if (! (lsame_(job, "B", (ftnlen)1, (ftnlen)1) || lsame_(job, "N", (
	    ftnlen)1, (ftnlen)1))) {
#line 329 "AB09MD.f"
	*info = -2;
#line 330 "AB09MD.f"
    } else if (! (lsame_(equil, "S", (ftnlen)1, (ftnlen)1) || lsame_(equil, 
	    "N", (ftnlen)1, (ftnlen)1))) {
#line 332 "AB09MD.f"
	*info = -3;
#line 333 "AB09MD.f"
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
#line 334 "AB09MD.f"
	*info = -4;
#line 335 "AB09MD.f"
    } else if (*n < 0) {
#line 336 "AB09MD.f"
	*info = -5;
#line 337 "AB09MD.f"
    } else if (*m < 0) {
#line 338 "AB09MD.f"
	*info = -6;
#line 339 "AB09MD.f"
    } else if (*p < 0) {
#line 340 "AB09MD.f"
	*info = -7;
#line 341 "AB09MD.f"
    } else if (fixord && (*nr < 0 || *nr > *n)) {
#line 342 "AB09MD.f"
	*info = -8;
#line 343 "AB09MD.f"
    } else if (discr && (*alpha < 0. || *alpha > 1.) || ! discr && *alpha > 
	    0.) {
#line 345 "AB09MD.f"
	*info = -9;
#line 346 "AB09MD.f"
    } else if (*lda < max(1,*n)) {
#line 347 "AB09MD.f"
	*info = -11;
#line 348 "AB09MD.f"
    } else if (*ldb < max(1,*n)) {
#line 349 "AB09MD.f"
	*info = -13;
#line 350 "AB09MD.f"
    } else if (*ldc < max(1,*p)) {
#line 351 "AB09MD.f"
	*info = -15;
#line 352 "AB09MD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
#line 352 "AB09MD.f"
	i__3 = max(*n,*m);
#line 352 "AB09MD.f"
	i__1 = 1, i__2 = *n * ((*n << 1) + max(i__3,*p) + 5) + *n * (*n + 1) /
		 2;
#line 352 "AB09MD.f"
	if (*ldwork < max(i__1,i__2)) {
#line 354 "AB09MD.f"
	    *info = -21;
#line 355 "AB09MD.f"
	}
#line 355 "AB09MD.f"
    }

#line 357 "AB09MD.f"
    if (*info != 0) {

/*        Error return. */

#line 361 "AB09MD.f"
	i__1 = -(*info);
#line 361 "AB09MD.f"
	xerbla_("AB09MD", &i__1, (ftnlen)6);
#line 362 "AB09MD.f"
	return 0;
#line 363 "AB09MD.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 367 "AB09MD.f"
    i__1 = min(*n,*m);
#line 367 "AB09MD.f"
    if (min(i__1,*p) == 0) {
#line 368 "AB09MD.f"
	*nr = 0;
#line 369 "AB09MD.f"
	dwork[1] = 1.;
#line 370 "AB09MD.f"
	return 0;
#line 371 "AB09MD.f"
    }

#line 373 "AB09MD.f"
    if (lsame_(equil, "S", (ftnlen)1, (ftnlen)1)) {

/*        Scale simultaneously the matrices A, B and C: */
/*        A <- inv(D)*A*D,  B <- inv(D)*B and C <- C*D, where D is a */
/*        diagonal matrix. */
/*        Workspace: N. */

#line 380 "AB09MD.f"
	maxred = 100.;
#line 381 "AB09MD.f"
	tb01id_("All", n, m, p, &maxred, &a[a_offset], lda, &b[b_offset], ldb,
		 &c__[c_offset], ldc, &dwork[1], info, (ftnlen)3);
#line 383 "AB09MD.f"
    }

/*     Correct the value of ALPHA to ensure stability. */

#line 387 "AB09MD.f"
    alpwrk = *alpha;
#line 388 "AB09MD.f"
    if (discr) {
#line 389 "AB09MD.f"
	if (*alpha == 1.) {
#line 389 "AB09MD.f"
	    alpwrk = 1. - sqrt(dlamch_("E", (ftnlen)1));
#line 389 "AB09MD.f"
	}
#line 390 "AB09MD.f"
    } else {
#line 391 "AB09MD.f"
	if (*alpha == 0.) {
#line 391 "AB09MD.f"
	    alpwrk = -sqrt(dlamch_("E", (ftnlen)1));
#line 391 "AB09MD.f"
	}
#line 392 "AB09MD.f"
    }

/*     Allocate working storage. */

#line 396 "AB09MD.f"
    nn = *n * *n;
#line 397 "AB09MD.f"
    ku = 1;
#line 398 "AB09MD.f"
    kwr = ku + nn;
#line 399 "AB09MD.f"
    kwi = kwr + *n;
#line 400 "AB09MD.f"
    kw = kwi + *n;
#line 401 "AB09MD.f"
    lwr = *ldwork - kw + 1;

/*     Reduce A to a block-diagonal real Schur form, with the */
/*     ALPHA-unstable part in the leading diagonal position, using a */
/*     non-orthogonal similarity transformation A <- inv(T)*A*T and */
/*     apply the transformation to B and C: B <- inv(T)*B and C <- C*T. */

/*     Workspace needed:      N*(N+2); */
/*     Additional workspace:  need   3*N; */
/*                            prefer larger. */

#line 412 "AB09MD.f"
    tb01kd_(dico, "Unstable", "General", n, m, p, &alpwrk, &a[a_offset], lda, 
	    &b[b_offset], ldb, &c__[c_offset], ldc, &nu, &dwork[ku], n, &
	    dwork[kwr], &dwork[kwi], &dwork[kw], &lwr, &ierr, (ftnlen)1, (
	    ftnlen)8, (ftnlen)7);

#line 416 "AB09MD.f"
    if (ierr != 0) {
#line 417 "AB09MD.f"
	if (ierr != 3) {
#line 418 "AB09MD.f"
	    *info = 1;
#line 419 "AB09MD.f"
	} else {
#line 420 "AB09MD.f"
	    *info = 2;
#line 421 "AB09MD.f"
	}
#line 422 "AB09MD.f"
	return 0;
#line 423 "AB09MD.f"
    }

#line 425 "AB09MD.f"
    wrkopt = (integer) (dwork[kw] + (doublereal) (kw - 1));

#line 427 "AB09MD.f"
    iwarnl = 0;
#line 428 "AB09MD.f"
    *ns = *n - nu;
#line 429 "AB09MD.f"
    if (fixord) {
/* Computing MAX */
#line 430 "AB09MD.f"
	i__1 = 0, i__2 = *nr - nu;
#line 430 "AB09MD.f"
	nra = max(i__1,i__2);
#line 431 "AB09MD.f"
	if (*nr < nu) {
#line 431 "AB09MD.f"
	    iwarnl = 2;
#line 431 "AB09MD.f"
	}
#line 433 "AB09MD.f"
    } else {
#line 434 "AB09MD.f"
	nra = 0;
#line 435 "AB09MD.f"
    }

/*     Finish if only unstable part is present. */

#line 439 "AB09MD.f"
    if (*ns == 0) {
#line 440 "AB09MD.f"
	*nr = nu;
#line 441 "AB09MD.f"
	dwork[1] = (doublereal) wrkopt;
#line 442 "AB09MD.f"
	return 0;
#line 443 "AB09MD.f"
    }

#line 445 "AB09MD.f"
    nu1 = nu + 1;

/*     Allocate working storage. */

#line 449 "AB09MD.f"
    kt = 1;
#line 450 "AB09MD.f"
    kti = kt + nn;
#line 451 "AB09MD.f"
    kw = kti + nn;

/*     Compute a B & T approximation of the stable part. */
/*     Workspace: need   N*(2*N+MAX(N,M,P)+5) + N*(N+1)/2; */
/*                prefer larger. */

#line 457 "AB09MD.f"
    i__1 = *ldwork - kw + 1;
#line 457 "AB09MD.f"
    ab09ax_(dico, job, ordsel, ns, m, p, &nra, &a[nu1 + nu1 * a_dim1], lda, &
	    b[nu1 + b_dim1], ldb, &c__[nu1 * c_dim1 + 1], ldc, &hsv[1], &
	    dwork[kt], n, &dwork[kti], n, tol, &iwork[1], &dwork[kw], &i__1, 
	    iwarn, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 461 "AB09MD.f"
    *iwarn = max(*iwarn,iwarnl);

#line 463 "AB09MD.f"
    if (ierr != 0) {
#line 464 "AB09MD.f"
	*info = ierr + 1;
#line 465 "AB09MD.f"
	return 0;
#line 466 "AB09MD.f"
    }

#line 468 "AB09MD.f"
    *nr = nra + nu;

/* Computing MAX */
#line 470 "AB09MD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 470 "AB09MD.f"
    dwork[1] = (doublereal) max(i__1,i__2);

#line 472 "AB09MD.f"
    return 0;
/* *** Last line of AB09MD *** */
} /* ab09md_ */

