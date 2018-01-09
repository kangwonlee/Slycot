#line 1 "AB09ED.f"
/* AB09ED.f -- translated by f2c (version 20100827).
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

#line 1 "AB09ED.f"
/* Subroutine */ int ab09ed_(char *dico, char *equil, char *ordsel, integer *
	n, integer *m, integer *p, integer *nr, doublereal *alpha, doublereal 
	*a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, doublereal *d__, integer *ldd, integer *ns, doublereal *
	hsv, doublereal *tol1, doublereal *tol2, integer *iwork, doublereal *
	dwork, integer *ldwork, integer *iwarn, integer *info, ftnlen 
	dico_len, ftnlen equil_len, ftnlen ordsel_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2, i__3, i__4, i__5;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer ki, kl, ku, kw, nu, nu1, nra, ierr;
    extern /* Subroutine */ int tb01id_(char *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    ab09cx_(char *, char *, integer *, integer *, integer *, integer *
	    , doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, ftnlen, ftnlen), tb01kd_(char *, char *, char *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal maxred;
    static logical fixord;
    static integer iwarnl;
    static doublereal alpwrk, wrkopt;


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

/*     To compute a reduced order model (Ar,Br,Cr,Dr) for an original */
/*     state-space representation (A,B,C,D) by using the optimal */
/*     Hankel-norm approximation method in conjunction with square-root */
/*     balancing for the ALPHA-stable part of the system. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the original system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to preliminarily */
/*             equilibrate the triplet (A,B,C) as follows: */
/*             = 'S':  perform equilibration (scaling); */
/*             = 'N':  do not perform equilibration. */

/*     ORDSEL  CHARACTER*1 */
/*             Specifies the order selection method as follows: */
/*             = 'F':  the resulting order NR is fixed; */
/*             = 'A':  the resulting order NR is automatically determined */
/*                     on basis of the given tolerance TOL1. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the original state-space representation, i.e. */
/*             the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     NR      (input/output) INTEGER */
/*             On entry with ORDSEL = 'F', NR is the desired order of */
/*             the resulting reduced order system.  0 <= NR <= N. */
/*             On exit, if INFO = 0, NR is the order of the resulting */
/*             reduced order model. For a system with NU ALPHA-unstable */
/*             eigenvalues and NS ALPHA-stable eigenvalues (NU+NS = N), */
/*             NR is set as follows: if ORDSEL = 'F', NR is equal to */
/*             NU+MIN(MAX(0,NR-NU-KR+1),NMIN), where KR is the */
/*             multiplicity of the Hankel singular value HSV(NR-NU+1), */
/*             NR is the desired order on entry, and NMIN is the order */
/*             of a minimal realization of the ALPHA-stable part of the */
/*             given system; NMIN is determined as the number of Hankel */
/*             singular values greater than NS*EPS*HNORM(As,Bs,Cs), where */
/*             EPS is the machine precision (see LAPACK Library Routine */
/*             DLAMCH) and HNORM(As,Bs,Cs) is the Hankel norm of the */
/*             ALPHA-stable part of the given system (computed in */
/*             HSV(1)); */
/*             if ORDSEL = 'A', NR is the sum of NU and the number of */
/*             Hankel singular values greater than */
/*             MAX(TOL1,NS*EPS*HNORM(As,Bs,Cs)). */

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
/*             array contains the state dynamics matrix Ar of the */
/*             reduced order system in a real Schur form. */
/*             The resulting A has a block-diagonal form with two blocks. */
/*             For a system with NU ALPHA-unstable eigenvalues and */
/*             NS ALPHA-stable eigenvalues (NU+NS = N), the leading */
/*             NU-by-NU block contains the unreduced part of A */
/*             corresponding to ALPHA-unstable eigenvalues. */
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

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the original input/output matrix D. */
/*             On exit, if INFO = 0, the leading P-by-M part of this */
/*             array contains the input/output matrix Dr of the reduced */
/*             order system. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P). */

/*     NS      (output) INTEGER */
/*             The dimension of the ALPHA-stable subsystem. */

/*     HSV     (output) DOUBLE PRECISION array, dimension (N) */
/*             If INFO = 0, the leading NS elements of HSV contain the */
/*             Hankel singular values of the ALPHA-stable part of the */
/*             original system ordered decreasingly. */
/*             HSV(1) is the Hankel norm of the ALPHA-stable subsystem. */

/*     Tolerances */

/*     TOL1    DOUBLE PRECISION */
/*             If ORDSEL = 'A', TOL1 contains the tolerance for */
/*             determining the order of reduced system. */
/*             For model reduction, the recommended value is */
/*             TOL1 = c*HNORM(As,Bs,Cs), where c is a constant in the */
/*             interval [0.00001,0.001], and HNORM(As,Bs,Cs) is the */
/*             Hankel-norm of the ALPHA-stable part of the given system */
/*             (computed in HSV(1)). */
/*             If TOL1 <= 0 on entry, the used default value is */
/*             TOL1 = NS*EPS*HNORM(As,Bs,Cs), where NS is the number of */
/*             ALPHA-stable eigenvalues of A and EPS is the machine */
/*             precision (see LAPACK Library Routine DLAMCH). */
/*             This value is appropriate to compute a minimal realization */
/*             of the ALPHA-stable part. */
/*             If ORDSEL = 'F', the value of TOL1 is ignored. */

/*     TOL2    DOUBLE PRECISION */
/*             The tolerance for determining the order of a minimal */
/*             realization of the ALPHA-stable part of the given system. */
/*             The recommended value is TOL2 = NS*EPS*HNORM(As,Bs,Cs). */
/*             This value is used by default if TOL2 <= 0 on entry. */
/*             If TOL2 > 0, then TOL2 <= TOL1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK = MAX(1,M),   if DICO = 'C'; */
/*             LIWORK = MAX(1,N,M), if DICO = 'D'. */
/*             On exit, if INFO = 0, IWORK(1) contains NMIN, the order of */
/*             the computed minimal realization. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( LDW1, LDW2 ), where */
/*             LDW1 = N*(2*N + MAX(N,M,P) + 5) + N*(N+1)/2, */
/*             LDW2 = N*(M+P+2) + 2*M*P + MIN(N,M) + */
/*                    MAX( 3*M+1, MIN(N,M)+P ). */
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
/*             = 3:  the computed ALPHA-stable part is just stable, */
/*                   having stable eigenvalues very near to the imaginary */
/*                   axis (if DICO = 'C') or to the unit circle */
/*                   (if DICO = 'D'); */
/*             = 4:  the computation of Hankel singular values failed; */
/*             = 5:  the computation of stable projection in the */
/*                   Hankel-norm approximation algorithm failed; */
/*             = 6:  the order of computed stable projection in the */
/*                   Hankel-norm approximation algorithm differs */
/*                   from the order of Hankel-norm approximation. */

/*     METHOD */

/*     Let be the following linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t) + Du(t)                           (1) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system. The subroutine AB09ED determines for */
/*     the given system (1), the matrices of a reduced order system */

/*          d[z(t)] = Ar*z(t) + Br*u(t) */
/*          yr(t)   = Cr*z(t) + Dr*u(t)                       (2) */

/*     such that */

/*     HSV(NR+NS-N) <= INFNORM(G-Gr) <= 2*[HSV(NR+NS-N+1)+...+HSV(NS)], */

/*     where G and Gr are transfer-function matrices of the systems */
/*     (A,B,C,D) and (Ar,Br,Cr,Dr), respectively, and INFNORM(G) is the */
/*     infinity-norm of G. */

/*     The following procedure is used to reduce a given G: */

/*     1) Decompose additively G as */

/*          G = G1 + G2 */

/*        such that G1 = (As,Bs,Cs,D) has only ALPHA-stable poles and */
/*        G2 = (Au,Bu,Cu,0) has only ALPHA-unstable poles. */

/*     2) Determine G1r, a reduced order approximation of the */
/*        ALPHA-stable part G1. */

/*     3) Assemble the reduced model Gr as */

/*           Gr = G1r + G2. */

/*     To reduce the ALPHA-stable part G1, the optimal Hankel-norm */
/*     approximation method of [1], based on the square-root */
/*     balancing projection formulas of [2], is employed. */

/*     REFERENCES */

/*     [1] Glover, K. */
/*         All optimal Hankel norm approximation of linear */
/*         multivariable systems and their L-infinity error bounds. */
/*         Int. J. Control, Vol. 36, pp. 1145-1193, 1984. */

/*     [2] Tombs M.S. and Postlethwaite I. */
/*         Truncated balanced realization of stable, non-minimal */
/*         state-space systems. */
/*         Int. J. Control, Vol. 46, pp. 1319-1330, 1987. */

/*     NUMERICAL ASPECTS */

/*     The implemented methods rely on an accuracy enhancing square-root */
/*     technique. */
/*                                         3 */
/*     The algorithms require less than 30N  floating point operations. */

/*     CONTRIBUTOR */

/*     C. Oara and A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, July 1998. */
/*     Based on the RASP routines SADSDC and OHNAP. */

/*     REVISIONS */

/*     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven. */
/*     Nov. 2000, A. Varga, DLR Oberpfaffenhofen. */
/*     March 26, 2005, V. Sima, Research Institute for Informatics. */

/*     KEYWORDS */

/*     Balancing, Hankel-norm approximation, model reduction, */
/*     multivariable system, state-space model. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 338 "AB09ED.f"
    /* Parameter adjustments */
#line 338 "AB09ED.f"
    a_dim1 = *lda;
#line 338 "AB09ED.f"
    a_offset = 1 + a_dim1;
#line 338 "AB09ED.f"
    a -= a_offset;
#line 338 "AB09ED.f"
    b_dim1 = *ldb;
#line 338 "AB09ED.f"
    b_offset = 1 + b_dim1;
#line 338 "AB09ED.f"
    b -= b_offset;
#line 338 "AB09ED.f"
    c_dim1 = *ldc;
#line 338 "AB09ED.f"
    c_offset = 1 + c_dim1;
#line 338 "AB09ED.f"
    c__ -= c_offset;
#line 338 "AB09ED.f"
    d_dim1 = *ldd;
#line 338 "AB09ED.f"
    d_offset = 1 + d_dim1;
#line 338 "AB09ED.f"
    d__ -= d_offset;
#line 338 "AB09ED.f"
    --hsv;
#line 338 "AB09ED.f"
    --iwork;
#line 338 "AB09ED.f"
    --dwork;
#line 338 "AB09ED.f"

#line 338 "AB09ED.f"
    /* Function Body */
#line 338 "AB09ED.f"
    *info = 0;
#line 339 "AB09ED.f"
    *iwarn = 0;
#line 340 "AB09ED.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 341 "AB09ED.f"
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);

/*     Check the input scalar arguments. */

#line 345 "AB09ED.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 346 "AB09ED.f"
	*info = -1;
#line 347 "AB09ED.f"
    } else if (! (lsame_(equil, "S", (ftnlen)1, (ftnlen)1) || lsame_(equil, 
	    "N", (ftnlen)1, (ftnlen)1))) {
#line 349 "AB09ED.f"
	*info = -2;
#line 350 "AB09ED.f"
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
#line 351 "AB09ED.f"
	*info = -3;
#line 352 "AB09ED.f"
    } else if (*n < 0) {
#line 353 "AB09ED.f"
	*info = -4;
#line 354 "AB09ED.f"
    } else if (*m < 0) {
#line 355 "AB09ED.f"
	*info = -5;
#line 356 "AB09ED.f"
    } else if (*p < 0) {
#line 357 "AB09ED.f"
	*info = -6;
#line 358 "AB09ED.f"
    } else if (fixord && (*nr < 0 || *nr > *n)) {
#line 359 "AB09ED.f"
	*info = -7;
#line 360 "AB09ED.f"
    } else if (discr && (*alpha < 0. || *alpha > 1.) || ! discr && *alpha > 
	    0.) {
#line 362 "AB09ED.f"
	*info = -8;
#line 363 "AB09ED.f"
    } else if (*lda < max(1,*n)) {
#line 364 "AB09ED.f"
	*info = -10;
#line 365 "AB09ED.f"
    } else if (*ldb < max(1,*n)) {
#line 366 "AB09ED.f"
	*info = -12;
#line 367 "AB09ED.f"
    } else if (*ldc < max(1,*p)) {
#line 368 "AB09ED.f"
	*info = -14;
#line 369 "AB09ED.f"
    } else if (*ldd < max(1,*p)) {
#line 370 "AB09ED.f"
	*info = -16;
#line 371 "AB09ED.f"
    } else if (*tol2 > 0. && *tol2 > *tol1) {
#line 372 "AB09ED.f"
	*info = -20;
#line 373 "AB09ED.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
#line 373 "AB09ED.f"
	i__3 = max(*n,*m);
/* Computing MAX */
#line 373 "AB09ED.f"
	i__4 = *m * 3 + 1, i__5 = min(*n,*m) + *p;
#line 373 "AB09ED.f"
	i__1 = *n * ((*n << 1) + max(i__3,*p) + 5) + *n * (*n + 1) / 2, i__2 =
		 *n * (*m + *p + 2) + (*m << 1) * *p + min(*n,*m) + max(i__4,
		i__5);
#line 373 "AB09ED.f"
	if (*ldwork < max(i__1,i__2)) {
#line 377 "AB09ED.f"
	    *info = -23;
#line 378 "AB09ED.f"
	}
#line 378 "AB09ED.f"
    }

#line 380 "AB09ED.f"
    if (*info != 0) {

/*        Error return. */

#line 384 "AB09ED.f"
	i__1 = -(*info);
#line 384 "AB09ED.f"
	xerbla_("AB09ED", &i__1, (ftnlen)6);
#line 385 "AB09ED.f"
	return 0;
#line 386 "AB09ED.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 390 "AB09ED.f"
    i__1 = min(*n,*m);
#line 390 "AB09ED.f"
    if (min(i__1,*p) == 0) {
#line 391 "AB09ED.f"
	*nr = 0;
#line 392 "AB09ED.f"
	*ns = 0;
#line 393 "AB09ED.f"
	iwork[1] = 0;
#line 394 "AB09ED.f"
	dwork[1] = 1.;
#line 395 "AB09ED.f"
	return 0;
#line 396 "AB09ED.f"
    }

#line 398 "AB09ED.f"
    if (lsame_(equil, "S", (ftnlen)1, (ftnlen)1)) {

/*        Scale simultaneously the matrices A, B and C: */
/*        A <- inv(D)*A*D,  B <- inv(D)*B and C <- C*D, where D is a */
/*        diagonal matrix. */
/*        Workspace: N. */

#line 405 "AB09ED.f"
	maxred = 100.;
#line 406 "AB09ED.f"
	tb01id_("All", n, m, p, &maxred, &a[a_offset], lda, &b[b_offset], ldb,
		 &c__[c_offset], ldc, &dwork[1], info, (ftnlen)3);
#line 408 "AB09ED.f"
    }

/*     Correct the value of ALPHA to ensure stability. */

#line 412 "AB09ED.f"
    alpwrk = *alpha;
#line 413 "AB09ED.f"
    if (discr) {
#line 414 "AB09ED.f"
	if (*alpha == 1.) {
#line 414 "AB09ED.f"
	    alpwrk = 1. - sqrt(dlamch_("E", (ftnlen)1));
#line 414 "AB09ED.f"
	}
#line 415 "AB09ED.f"
    } else {
#line 416 "AB09ED.f"
	if (*alpha == 0.) {
#line 416 "AB09ED.f"
	    alpwrk = -sqrt(dlamch_("E", (ftnlen)1));
#line 416 "AB09ED.f"
	}
#line 417 "AB09ED.f"
    }

/*     Allocate working storage. */

#line 421 "AB09ED.f"
    ku = 1;
#line 422 "AB09ED.f"
    kl = ku + *n * *n;
#line 423 "AB09ED.f"
    ki = kl + *n;
#line 424 "AB09ED.f"
    kw = ki + *n;

/*     Reduce A to a block-diagonal real Schur form, with the */
/*     ALPHA-unstable part in the leading diagonal position, using a */
/*     non-orthogonal similarity transformation A <- inv(T)*A*T and */
/*     apply the transformation to B and C: B <- inv(T)*B and C <- C*T. */

/*     Workspace needed:      N*(N+2); */
/*     Additional workspace:  need   3*N; */
/*                            prefer larger. */

#line 435 "AB09ED.f"
    i__1 = *ldwork - kw + 1;
#line 435 "AB09ED.f"
    tb01kd_(dico, "Unstable", "General", n, m, p, &alpwrk, &a[a_offset], lda, 
	    &b[b_offset], ldb, &c__[c_offset], ldc, &nu, &dwork[ku], n, &
	    dwork[kl], &dwork[ki], &dwork[kw], &i__1, &ierr, (ftnlen)1, (
	    ftnlen)8, (ftnlen)7);

#line 439 "AB09ED.f"
    if (ierr != 0) {
#line 440 "AB09ED.f"
	if (ierr != 3) {
#line 441 "AB09ED.f"
	    *info = 1;
#line 442 "AB09ED.f"
	} else {
#line 443 "AB09ED.f"
	    *info = 2;
#line 444 "AB09ED.f"
	}
#line 445 "AB09ED.f"
	return 0;
#line 446 "AB09ED.f"
    }

#line 448 "AB09ED.f"
    wrkopt = dwork[kw] + (doublereal) (kw - 1);

/*     Determine a reduced order approximation of the ALPHA-stable part. */

/*     Workspace: need   MAX( LDW1, LDW2 ), */
/*                LDW1 = N*(2*N + MAX(N,M,P) + 5) + N*(N+1)/2, */
/*                LDW2 = N*(M+P+2) + 2*M*P + MIN(N,M) + */
/*                       MAX( 3*M+1, MIN(N,M)+P ); */
/*                prefer larger. */

#line 458 "AB09ED.f"
    iwarnl = 0;
#line 459 "AB09ED.f"
    *ns = *n - nu;
#line 460 "AB09ED.f"
    if (fixord) {
/* Computing MAX */
#line 461 "AB09ED.f"
	i__1 = 0, i__2 = *nr - nu;
#line 461 "AB09ED.f"
	nra = max(i__1,i__2);
#line 462 "AB09ED.f"
	if (*nr < nu) {
#line 462 "AB09ED.f"
	    iwarnl = 2;
#line 462 "AB09ED.f"
	}
#line 464 "AB09ED.f"
    } else {
#line 465 "AB09ED.f"
	nra = 0;
#line 466 "AB09ED.f"
    }

/*     Finish if only unstable part is present. */

#line 470 "AB09ED.f"
    if (*ns == 0) {
#line 471 "AB09ED.f"
	*nr = nu;
#line 472 "AB09ED.f"
	dwork[1] = wrkopt;
#line 473 "AB09ED.f"
	return 0;
#line 474 "AB09ED.f"
    }

#line 476 "AB09ED.f"
    nu1 = nu + 1;
#line 477 "AB09ED.f"
    ab09cx_(dico, ordsel, ns, m, p, &nra, &a[nu1 + nu1 * a_dim1], lda, &b[nu1 
	    + b_dim1], ldb, &c__[nu1 * c_dim1 + 1], ldc, &d__[d_offset], ldd, 
	    &hsv[1], tol1, tol2, &iwork[1], &dwork[1], ldwork, iwarn, &ierr, (
	    ftnlen)1, (ftnlen)1);

#line 481 "AB09ED.f"
    *iwarn = max(*iwarn,iwarnl);
#line 482 "AB09ED.f"
    if (ierr != 0) {
#line 483 "AB09ED.f"
	*info = ierr + 2;
#line 484 "AB09ED.f"
	return 0;
#line 485 "AB09ED.f"
    }

#line 487 "AB09ED.f"
    *nr = nra + nu;

#line 489 "AB09ED.f"
    dwork[1] = max(wrkopt,dwork[1]);

#line 491 "AB09ED.f"
    return 0;
/* *** Last line of AB09ED *** */
} /* ab09ed_ */
