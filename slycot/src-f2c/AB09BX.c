#line 1 "AB09BX.f"
/* AB09BX.f -- translated by f2c (version 20100827).
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

#line 1 "AB09BX.f"
/* Table of constant values */

static logical c_true = TRUE_;
static logical c_false = FALSE_;
static integer c__1 = 1;
static doublereal c_b33 = 1.;
static doublereal c_b52 = 0.;

/* Subroutine */ int ab09bx_(char *dico, char *job, char *ordsel, integer *n, 
	integer *m, integer *p, integer *nr, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *hsv, doublereal *t, 
	integer *ldt, doublereal *ti, integer *ldti, doublereal *tol1, 
	doublereal *tol2, integer *iwork, doublereal *dwork, integer *ldwork, 
	integer *iwarn, integer *info, ftnlen dico_len, ftnlen job_len, 
	ftnlen ordsel_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, t_dim1, t_offset, ti_dim1, ti_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, k, ij, ku, kv, kw, ns, nr1;
    static logical bal;
    static integer ldw;
    static doublereal atol;
    static integer ierr, ktau;
    static doublereal temp, rtol;
    extern /* Subroutine */ int ab09dd_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), ma02ad_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen), ma02dd_(char *, char *, integer *, doublereal 
	    *, integer *, doublereal *, ftnlen, ftnlen), dscal_(integer *, 
	    doublereal *, doublereal *, integer *), dgemm_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen), mb03ud_(char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static logical discr;
    static doublereal rcond;
    extern /* Subroutine */ int sb03ou_(logical *, logical *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, integer *);
    static integer nminr;
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dtpmv_(
	    char *, char *, char *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen), dtrmv_(char *, char *, char *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static logical packed;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal scalec, scaleo;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dgetrf_(integer *, integer *, doublereal *, integer *, integer *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen), dgetrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static logical fixord;
    extern /* Subroutine */ int dorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
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

/*     To compute a reduced order model (Ar,Br,Cr,Dr) for a stable */
/*     original state-space representation (A,B,C,D) by using either the */
/*     square-root or the balancing-free square-root */
/*     Singular Perturbation Approximation (SPA) model reduction method. */
/*     The state dynamics matrix A of the original system is an upper */
/*     quasi-triangular matrix in real Schur canonical form. The matrices */
/*     of a minimal realization are computed using the truncation */
/*     formulas: */

/*          Am = TI * A * T ,  Bm = TI * B ,  Cm = C * T .      (1) */

/*     Am, Bm, Cm and D serve further for computing the SPA of the given */
/*     system. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the original system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     JOB     CHARACTER*1 */
/*             Specifies the model reduction approach to be used */
/*             as follows: */
/*             = 'B':  use the square-root SPA method; */
/*             = 'N':  use the balancing-free square-root SPA method. */

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
/*             reduced order model. NR is set as follows: */
/*             if ORDSEL = 'F', NR is equal to MIN(NR,NMIN), where NR */
/*             is the desired order on entry and NMIN is the order of a */
/*             minimal realization of the given system; NMIN is */
/*             determined as the number of Hankel singular values greater */
/*             than N*EPS*HNORM(A,B,C), where EPS is the machine */
/*             precision (see LAPACK Library Routine DLAMCH) and */
/*             HNORM(A,B,C) is the Hankel norm of the system (computed */
/*             in HSV(1)); */
/*             if ORDSEL = 'A', NR is equal to the number of Hankel */
/*             singular values greater than MAX(TOL1,N*EPS*HNORM(A,B,C)). */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix A in a real Schur */
/*             canonical form. */
/*             On exit, if INFO = 0, the leading NR-by-NR part of this */
/*             array contains the state dynamics matrix Ar of the */
/*             reduced order system. */

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

/*     HSV     (output) DOUBLE PRECISION array, dimension (N) */
/*             If INFO = 0, it contains the Hankel singular values of */
/*             the original system ordered decreasingly. HSV(1) is the */
/*             Hankel norm of the system. */

/*     T       (output) DOUBLE PRECISION array, dimension (LDT,N) */
/*             If INFO = 0 and NR > 0, the leading N-by-NR part of this */
/*             array contains the right truncation matrix T in (1). */

/*     LDT     INTEGER */
/*             The leading dimension of array T.  LDT >= MAX(1,N). */

/*     TI      (output) DOUBLE PRECISION array, dimension (LDTI,N) */
/*             If INFO = 0 and NR > 0, the leading NR-by-N part of this */
/*             array contains the left truncation matrix TI in (1). */

/*     LDTI    INTEGER */
/*             The leading dimension of array TI.  LDTI >= MAX(1,N). */

/*     Tolerances */

/*     TOL1    DOUBLE PRECISION */
/*             If ORDSEL = 'A', TOL1 contains the tolerance for */
/*             determining the order of reduced system. */
/*             For model reduction, the recommended value is */
/*             TOL1 = c*HNORM(A,B,C), where c is a constant in the */
/*             interval [0.00001,0.001], and HNORM(A,B,C) is the */
/*             Hankel-norm of the given system (computed in HSV(1)). */
/*             For computing a minimal realization, the recommended */
/*             value is TOL1 = N*EPS*HNORM(A,B,C), where EPS is the */
/*             machine precision (see LAPACK Library Routine DLAMCH). */
/*             This value is used by default if TOL1 <= 0 on entry. */
/*             If ORDSEL = 'F', the value of TOL1 is ignored. */

/*     TOL2    DOUBLE PRECISION */
/*             The tolerance for determining the order of a minimal */
/*             realization of the given system. The recommended value is */
/*             TOL2 = N*EPS*HNORM(A,B,C). This value is used by default */
/*             if TOL2 <= 0 on entry. */
/*             If TOL2 > 0, then TOL2 <= TOL1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension MAX(1,2*N) */
/*             On exit with INFO = 0, IWORK(1) contains the order of the */
/*             minimal realization of the system. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,N*(MAX(N,M,P)+5) + N*(N+1)/2). */
/*             For optimum performance LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  with ORDSEL = 'F', the selected order NR is greater */
/*                   than the order of a minimal realization of the */
/*                   given system. In this case, the resulting NR is */
/*                   set automatically to a value corresponding to the */
/*                   order of a minimal realization of the system. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the state matrix A is not stable (if DICO = 'C') */
/*                   or not convergent (if DICO = 'D'); */
/*             = 2:  the computation of Hankel singular values failed. */

/*     METHOD */

/*     Let be the stable linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t) + Du(t)                              (2) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system. The subroutine AB09BX determines for */
/*     the given system (1), the matrices of a reduced NR order system */

/*          d[z(t)] = Ar*z(t) + Br*u(t) */
/*          yr(t)   = Cr*z(t) + Dr*u(t)                          (3) */

/*     such that */

/*           HSV(NR) <= INFNORM(G-Gr) <= 2*[HSV(NR+1) + ... + HSV(N)], */

/*     where G and Gr are transfer-function matrices of the systems */
/*     (A,B,C,D) and (Ar,Br,Cr,Dr), respectively, and INFNORM(G) is the */
/*     infinity-norm of G. */

/*     If JOB = 'B', the balancing-based square-root SPA method of [1] */
/*     is used and the resulting model is balanced. */

/*     If JOB = 'N', the balancing-free square-root SPA method of [2] */
/*     is used. */
/*     By setting TOL1 = TOL2, the routine can be also used to compute */
/*     Balance & Truncate approximations. */

/*     REFERENCES */

/*     [1] Liu Y. and Anderson B.D.O. */
/*         Singular Perturbation Approximation of Balanced Systems, */
/*         Int. J. Control, Vol. 50, pp. 1379-1405, 1989. */

/*     [2] Varga A. */
/*         Balancing-free square-root algorithm for computing singular */
/*         perturbation approximations. */
/*         Proc. 30-th IEEE CDC,  Brighton, Dec. 11-13, 1991, */
/*         Vol. 2, pp. 1062-1065. */

/*     NUMERICAL ASPECTS */

/*     The implemented methods rely on accuracy enhancing square-root or */
/*     balancing-free square-root techniques. */
/*                                         3 */
/*     The algorithms require less than 30N  floating point operations. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, March 1998. */
/*     Based on the RASP routine SRBFP1. */

/*     REVISIONS */

/*     May 2, 1998. */
/*     November 11, 1998, V. Sima, Research Institute for Informatics, */
/*     Bucharest. */
/*     December 1998, V. Sima, Katholieke Univ. Leuven, Leuven. */
/*     February 14, 1999, A. Varga, German Aerospace Center. */
/*     February 22, 1999, V. Sima, Research Institute for Informatics. */
/*     February 27, 2000, V. Sima, Research Institute for Informatics. */
/*     May 26, 2000, A. Varga, German Aerospace Center. */

/*     KEYWORDS */

/*     Balancing, minimal state-space representation, model reduction, */
/*     multivariable system, singular perturbation approximation, */
/*     state-space model. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 308 "AB09BX.f"
    /* Parameter adjustments */
#line 308 "AB09BX.f"
    a_dim1 = *lda;
#line 308 "AB09BX.f"
    a_offset = 1 + a_dim1;
#line 308 "AB09BX.f"
    a -= a_offset;
#line 308 "AB09BX.f"
    b_dim1 = *ldb;
#line 308 "AB09BX.f"
    b_offset = 1 + b_dim1;
#line 308 "AB09BX.f"
    b -= b_offset;
#line 308 "AB09BX.f"
    c_dim1 = *ldc;
#line 308 "AB09BX.f"
    c_offset = 1 + c_dim1;
#line 308 "AB09BX.f"
    c__ -= c_offset;
#line 308 "AB09BX.f"
    d_dim1 = *ldd;
#line 308 "AB09BX.f"
    d_offset = 1 + d_dim1;
#line 308 "AB09BX.f"
    d__ -= d_offset;
#line 308 "AB09BX.f"
    --hsv;
#line 308 "AB09BX.f"
    t_dim1 = *ldt;
#line 308 "AB09BX.f"
    t_offset = 1 + t_dim1;
#line 308 "AB09BX.f"
    t -= t_offset;
#line 308 "AB09BX.f"
    ti_dim1 = *ldti;
#line 308 "AB09BX.f"
    ti_offset = 1 + ti_dim1;
#line 308 "AB09BX.f"
    ti -= ti_offset;
#line 308 "AB09BX.f"
    --iwork;
#line 308 "AB09BX.f"
    --dwork;
#line 308 "AB09BX.f"

#line 308 "AB09BX.f"
    /* Function Body */
#line 308 "AB09BX.f"
    *info = 0;
#line 309 "AB09BX.f"
    *iwarn = 0;
#line 310 "AB09BX.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 311 "AB09BX.f"
    bal = lsame_(job, "B", (ftnlen)1, (ftnlen)1);
#line 312 "AB09BX.f"
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 316 "AB09BX.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 317 "AB09BX.f"
	*info = -1;
#line 318 "AB09BX.f"
    } else if (! (bal || lsame_(job, "N", (ftnlen)1, (ftnlen)1))) {
#line 319 "AB09BX.f"
	*info = -2;
#line 320 "AB09BX.f"
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
#line 321 "AB09BX.f"
	*info = -3;
#line 322 "AB09BX.f"
    } else if (*n < 0) {
#line 323 "AB09BX.f"
	*info = -4;
#line 324 "AB09BX.f"
    } else if (*m < 0) {
#line 325 "AB09BX.f"
	*info = -5;
#line 326 "AB09BX.f"
    } else if (*p < 0) {
#line 327 "AB09BX.f"
	*info = -6;
#line 328 "AB09BX.f"
    } else if (fixord && (*nr < 0 || *nr > *n)) {
#line 329 "AB09BX.f"
	*info = -7;
#line 330 "AB09BX.f"
    } else if (*lda < max(1,*n)) {
#line 331 "AB09BX.f"
	*info = -9;
#line 332 "AB09BX.f"
    } else if (*ldb < max(1,*n)) {
#line 333 "AB09BX.f"
	*info = -11;
#line 334 "AB09BX.f"
    } else if (*ldc < max(1,*p)) {
#line 335 "AB09BX.f"
	*info = -13;
#line 336 "AB09BX.f"
    } else if (*ldd < max(1,*p)) {
#line 337 "AB09BX.f"
	*info = -15;
#line 338 "AB09BX.f"
    } else if (*ldt < max(1,*n)) {
#line 339 "AB09BX.f"
	*info = -18;
#line 340 "AB09BX.f"
    } else if (*ldti < max(1,*n)) {
#line 341 "AB09BX.f"
	*info = -20;
#line 342 "AB09BX.f"
    } else if (*tol2 > 0. && *tol2 > *tol1) {
#line 343 "AB09BX.f"
	*info = -22;
#line 344 "AB09BX.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
#line 344 "AB09BX.f"
	i__3 = max(*n,*m);
#line 344 "AB09BX.f"
	i__1 = 1, i__2 = *n * (max(i__3,*p) + 5) + *n * (*n + 1) / 2;
#line 344 "AB09BX.f"
	if (*ldwork < max(i__1,i__2)) {
#line 346 "AB09BX.f"
	    *info = -25;
#line 347 "AB09BX.f"
	}
#line 347 "AB09BX.f"
    }

#line 349 "AB09BX.f"
    if (*info != 0) {

/*        Error return. */

#line 353 "AB09BX.f"
	i__1 = -(*info);
#line 353 "AB09BX.f"
	xerbla_("AB09BX", &i__1, (ftnlen)6);
#line 354 "AB09BX.f"
	return 0;
#line 355 "AB09BX.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 359 "AB09BX.f"
    i__1 = min(*n,*m);
#line 359 "AB09BX.f"
    if (min(i__1,*p) == 0) {
#line 360 "AB09BX.f"
	*nr = 0;
#line 361 "AB09BX.f"
	iwork[1] = 0;
#line 362 "AB09BX.f"
	dwork[1] = 1.;
#line 363 "AB09BX.f"
	return 0;
#line 364 "AB09BX.f"
    }

#line 366 "AB09BX.f"
    rtol = (doublereal) (*n) * dlamch_("Epsilon", (ftnlen)7);

/*     Allocate N*MAX(N,M,P) and N working storage for the matrices U */
/*     and TAU, respectively. */

#line 371 "AB09BX.f"
    ku = 1;
/* Computing MAX */
#line 372 "AB09BX.f"
    i__1 = max(*n,*m);
#line 372 "AB09BX.f"
    ktau = ku + *n * max(i__1,*p);
#line 373 "AB09BX.f"
    kw = ktau + *n;
#line 374 "AB09BX.f"
    ldw = *ldwork - kw + 1;

/*     Copy B in U. */

#line 378 "AB09BX.f"
    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[ku], n, (ftnlen)4);

/*     If DISCR = .FALSE., solve for Su the Lyapunov equation */
/*                                      2 */
/*     A*(Su*Su') + (Su*Su')*A' + scalec *B*B' = 0 . */

/*     If DISCR = .TRUE., solve for Su the Lyapunov equation */
/*                           2 */
/*     A*(Su*Su')*A' + scalec *B*B' = Su*Su' . */

/*     Workspace:  need   N*(MAX(N,M,P) + 5); */
/*                 prefer larger. */

#line 391 "AB09BX.f"
    sb03ou_(&discr, &c_true, n, m, &a[a_offset], lda, &dwork[ku], n, &dwork[
	    ktau], &ti[ti_offset], ldti, &scalec, &dwork[kw], &ldw, &ierr);
#line 393 "AB09BX.f"
    if (ierr != 0) {
#line 394 "AB09BX.f"
	*info = 1;
#line 395 "AB09BX.f"
	return 0;
#line 396 "AB09BX.f"
    }
#line 397 "AB09BX.f"
    wrkopt = (integer) dwork[kw] + kw - 1;

/*     Copy C in U. */

#line 401 "AB09BX.f"
    dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[ku], p, (ftnlen)4);

/*     If DISCR = .FALSE., solve for Ru the Lyapunov equation */
/*                                      2 */
/*     A'*(Ru'*Ru) + (Ru'*Ru)*A + scaleo  * C'*C = 0 . */

/*     If DISCR = .TRUE., solve for Ru the Lyapunov equation */
/*                           2 */
/*     A'*(Ru'*Ru)*A + scaleo  * C'*C = Ru'*Ru . */

/*     Workspace:  need   N*(MAX(N,M,P) + 5); */
/*                 prefer larger. */

#line 414 "AB09BX.f"
    sb03ou_(&discr, &c_false, n, p, &a[a_offset], lda, &dwork[ku], p, &dwork[
	    ktau], &t[t_offset], ldt, &scaleo, &dwork[kw], &ldw, &ierr);
/* Computing MAX */
#line 416 "AB09BX.f"
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 416 "AB09BX.f"
    wrkopt = max(i__1,i__2);

/*     Allocate N*(N+1)/2 (or, if possible, N*N) working storage for the */
/*     matrix V, a packed (or unpacked) copy of Su, and save Su in V. */
/*     (The locations for TAU are reused here.) */

#line 422 "AB09BX.f"
    kv = ktau;
#line 423 "AB09BX.f"
    if (*ldwork - kv + 1 < *n * (*n + 5)) {
#line 424 "AB09BX.f"
	packed = TRUE_;
#line 425 "AB09BX.f"
	ma02dd_("Pack", "Upper", n, &ti[ti_offset], ldti, &dwork[kv], (ftnlen)
		4, (ftnlen)5);
#line 426 "AB09BX.f"
	kw = kv + *n * (*n + 1) / 2;
#line 427 "AB09BX.f"
    } else {
#line 428 "AB09BX.f"
	packed = FALSE_;
#line 429 "AB09BX.f"
	dlacpy_("Upper", n, n, &ti[ti_offset], ldti, &dwork[kv], n, (ftnlen)5)
		;
#line 430 "AB09BX.f"
	kw = kv + *n * *n;
#line 431 "AB09BX.f"
    }
/*                               | x x | */
/*     Compute Ru*Su in the form | 0 x | in TI. */

#line 435 "AB09BX.f"
    i__1 = *n;
#line 435 "AB09BX.f"
    for (j = 1; j <= i__1; ++j) {
#line 436 "AB09BX.f"
	dtrmv_("Upper", "NoTranspose", "NonUnit", &j, &t[t_offset], ldt, &ti[
		j * ti_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)11, (ftnlen)7);
#line 438 "AB09BX.f"
/* L10: */
#line 438 "AB09BX.f"
    }

/*     Compute the singular value decomposition Ru*Su = V*S*UT */
/*     of the upper triangular matrix Ru*Su, with UT in TI and V in U. */

/*     Workspace:  need   N*MAX(N,M,P) + N*(N+1)/2 + 5*N; */
/*                 prefer larger. */

#line 446 "AB09BX.f"
    i__1 = *ldwork - kw + 1;
#line 446 "AB09BX.f"
    mb03ud_("Vectors", "Vectors", n, &ti[ti_offset], ldti, &dwork[ku], n, &
	    hsv[1], &dwork[kw], &i__1, &ierr, (ftnlen)7, (ftnlen)7);
#line 448 "AB09BX.f"
    if (ierr != 0) {
#line 449 "AB09BX.f"
	*info = 2;
#line 450 "AB09BX.f"
	return 0;
#line 451 "AB09BX.f"
    }
/* Computing MAX */
#line 452 "AB09BX.f"
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 452 "AB09BX.f"
    wrkopt = max(i__1,i__2);

/*     Scale singular values. */

#line 456 "AB09BX.f"
    d__1 = 1. / scalec / scaleo;
#line 456 "AB09BX.f"
    dscal_(n, &d__1, &hsv[1], &c__1);

/*     Partition S, U and V conformally as: */

/*     S = diag(S1,S2,S3),  U = [U1,U2,U3] (U' in TI) and V = [V1,V2,V3] */
/*     (in U). */

/*     Compute the order NR of reduced system, as the order of S1. */

#line 465 "AB09BX.f"
    atol = rtol * hsv[1];
#line 466 "AB09BX.f"
    if (fixord) {
#line 467 "AB09BX.f"
	if (*nr > 0) {
#line 468 "AB09BX.f"
	    if (hsv[*nr] <= atol) {
#line 469 "AB09BX.f"
		*nr = 0;
#line 470 "AB09BX.f"
		*iwarn = 1;
#line 471 "AB09BX.f"
		fixord = FALSE_;
#line 472 "AB09BX.f"
	    }
#line 473 "AB09BX.f"
	}
#line 474 "AB09BX.f"
    } else {
#line 475 "AB09BX.f"
	atol = max(*tol1,atol);
#line 476 "AB09BX.f"
	*nr = 0;
#line 477 "AB09BX.f"
    }
#line 478 "AB09BX.f"
    if (! fixord) {
#line 479 "AB09BX.f"
	i__1 = *n;
#line 479 "AB09BX.f"
	for (j = 1; j <= i__1; ++j) {
#line 480 "AB09BX.f"
	    if (hsv[j] <= atol) {
#line 480 "AB09BX.f"
		goto L30;
#line 480 "AB09BX.f"
	    }
#line 481 "AB09BX.f"
	    ++(*nr);
#line 482 "AB09BX.f"
/* L20: */
#line 482 "AB09BX.f"
	}
#line 483 "AB09BX.f"
L30:
#line 484 "AB09BX.f"
	;
#line 484 "AB09BX.f"
    }

/*     Finish if the order of the reduced model is zero. */

#line 488 "AB09BX.f"
    if (*nr == 0) {

/*       Compute only Dr using singular perturbation formulas. */
/*       Workspace:  need real    4*N; */
/*                   need integer 2*N. */

#line 494 "AB09BX.f"
	ab09dd_(dico, n, m, p, nr, &a[a_offset], lda, &b[b_offset], ldb, &c__[
		c_offset], ldc, &d__[d_offset], ldd, &rcond, &iwork[1], &
		dwork[1], &ierr, (ftnlen)1);
#line 496 "AB09BX.f"
	iwork[1] = 0;
#line 497 "AB09BX.f"
	dwork[1] = (doublereal) wrkopt;
#line 498 "AB09BX.f"
	return 0;
#line 499 "AB09BX.f"
    }

/*     Compute the order of minimal realization as the order of [S1 S2]. */

#line 503 "AB09BX.f"
    nr1 = *nr + 1;
#line 504 "AB09BX.f"
    nminr = *nr;
#line 505 "AB09BX.f"
    if (*nr < *n) {
/* Computing MAX */
#line 506 "AB09BX.f"
	d__1 = *tol2, d__2 = rtol * hsv[1];
#line 506 "AB09BX.f"
	atol = max(d__1,d__2);
#line 507 "AB09BX.f"
	i__1 = *n;
#line 507 "AB09BX.f"
	for (j = nr1; j <= i__1; ++j) {
#line 508 "AB09BX.f"
	    if (hsv[j] <= atol) {
#line 508 "AB09BX.f"
		goto L50;
#line 508 "AB09BX.f"
	    }
#line 509 "AB09BX.f"
	    ++nminr;
#line 510 "AB09BX.f"
/* L40: */
#line 510 "AB09BX.f"
	}
#line 511 "AB09BX.f"
L50:
#line 512 "AB09BX.f"
	;
#line 512 "AB09BX.f"
    }

/*     Compute the order of S2. */

#line 516 "AB09BX.f"
    ns = nminr - *nr;

/*     Compute the truncation matrices. */

/*     Compute TI' = | TI1' TI2' | = Ru'*| V1 V2 | in U. */

#line 522 "AB09BX.f"
    dtrmm_("Left", "Upper", "Transpose", "NonUnit", n, &nminr, &c_b33, &t[
	    t_offset], ldt, &dwork[ku], n, (ftnlen)4, (ftnlen)5, (ftnlen)9, (
	    ftnlen)7);

/*     Compute  T = | T1 T2 | = Su*| U1 U2 | */
/*     (with Su packed, if not enough workspace). */

#line 528 "AB09BX.f"
    ma02ad_("Full", &nminr, n, &ti[ti_offset], ldti, &t[t_offset], ldt, (
	    ftnlen)4);
#line 529 "AB09BX.f"
    if (packed) {
#line 530 "AB09BX.f"
	i__1 = nminr;
#line 530 "AB09BX.f"
	for (j = 1; j <= i__1; ++j) {
#line 531 "AB09BX.f"
	    dtpmv_("Upper", "NoTranspose", "NonUnit", n, &dwork[kv], &t[j * 
		    t_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)11, (ftnlen)7);
#line 533 "AB09BX.f"
/* L60: */
#line 533 "AB09BX.f"
	}
#line 534 "AB09BX.f"
    } else {
#line 535 "AB09BX.f"
	dtrmm_("Left", "Upper", "NoTranspose", "NonUnit", n, &nminr, &c_b33, &
		dwork[kv], n, &t[t_offset], ldt, (ftnlen)4, (ftnlen)5, (
		ftnlen)11, (ftnlen)7);
#line 537 "AB09BX.f"
    }

#line 539 "AB09BX.f"
    if (bal) {
#line 540 "AB09BX.f"
	ij = ku;

/*        Square-Root SPA method. */

/*        Compute the truncation matrices for balancing */
/*                    -1/2            -1/2 */
/*               T1*S1     and TI1'*S1 */

#line 548 "AB09BX.f"
	i__1 = *nr;
#line 548 "AB09BX.f"
	for (j = 1; j <= i__1; ++j) {
#line 549 "AB09BX.f"
	    temp = 1. / sqrt(hsv[j]);
#line 550 "AB09BX.f"
	    dscal_(n, &temp, &t[j * t_dim1 + 1], &c__1);
#line 551 "AB09BX.f"
	    dscal_(n, &temp, &dwork[ij], &c__1);
#line 552 "AB09BX.f"
	    ij += *n;
#line 553 "AB09BX.f"
/* L70: */
#line 553 "AB09BX.f"
	}
#line 554 "AB09BX.f"
    } else {

/*        Balancing-Free SPA method. */

/*        Compute orthogonal bases for the images of matrices T1 and */
/*        TI1'. */

/*        Workspace:  need   N*MAX(N,M,P) + 2*NR; */
/*                    prefer N*MAX(N,M,P) + NR*(NB+1) */
/*                           (NB determined by ILAENV for DGEQRF). */

#line 565 "AB09BX.f"
	kw = ktau + *nr;
#line 566 "AB09BX.f"
	ldw = *ldwork - kw + 1;
#line 567 "AB09BX.f"
	dgeqrf_(n, nr, &t[t_offset], ldt, &dwork[ktau], &dwork[kw], &ldw, &
		ierr);
#line 568 "AB09BX.f"
	dorgqr_(n, nr, nr, &t[t_offset], ldt, &dwork[ktau], &dwork[kw], &ldw, 
		&ierr);
#line 570 "AB09BX.f"
	dgeqrf_(n, nr, &dwork[ku], n, &dwork[ktau], &dwork[kw], &ldw, &ierr);
/* Computing MAX */
#line 572 "AB09BX.f"
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 572 "AB09BX.f"
	wrkopt = max(i__1,i__2);
#line 573 "AB09BX.f"
	dorgqr_(n, nr, nr, &dwork[ku], n, &dwork[ktau], &dwork[kw], &ldw, &
		ierr);
/* Computing MAX */
#line 575 "AB09BX.f"
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 575 "AB09BX.f"
	wrkopt = max(i__1,i__2);
#line 576 "AB09BX.f"
    }
#line 577 "AB09BX.f"
    if (ns > 0) {

/*        Compute orthogonal bases for the images of matrices T2 and */
/*        TI2'. */

/*        Workspace:  need   N*MAX(N,M,P) + 2*NS; */
/*                    prefer N*MAX(N,M,P) + NS*(NB+1) */
/*                           (NB determined by ILAENV for DGEQRF). */
#line 585 "AB09BX.f"
	kw = ktau + ns;
#line 586 "AB09BX.f"
	ldw = *ldwork - kw + 1;
#line 587 "AB09BX.f"
	dgeqrf_(n, &ns, &t[nr1 * t_dim1 + 1], ldt, &dwork[ktau], &dwork[kw], &
		ldw, &ierr);
#line 589 "AB09BX.f"
	dorgqr_(n, &ns, &ns, &t[nr1 * t_dim1 + 1], ldt, &dwork[ktau], &dwork[
		kw], &ldw, &ierr);
#line 591 "AB09BX.f"
	dgeqrf_(n, &ns, &dwork[ku + *n * *nr], n, &dwork[ktau], &dwork[kw], &
		ldw, &ierr);
/* Computing MAX */
#line 593 "AB09BX.f"
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 593 "AB09BX.f"
	wrkopt = max(i__1,i__2);
#line 594 "AB09BX.f"
	dorgqr_(n, &ns, &ns, &dwork[ku + *n * *nr], n, &dwork[ktau], &dwork[
		kw], &ldw, &ierr);
/* Computing MAX */
#line 596 "AB09BX.f"
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 596 "AB09BX.f"
	wrkopt = max(i__1,i__2);
#line 597 "AB09BX.f"
    }

/*     Transpose TI' in TI. */

#line 601 "AB09BX.f"
    ma02ad_("Full", n, &nminr, &dwork[ku], n, &ti[ti_offset], ldti, (ftnlen)4)
	    ;

#line 603 "AB09BX.f"
    if (! bal) {
/*                        -1 */
/*        Compute (TI1*T1)  *TI1 in TI. */

#line 607 "AB09BX.f"
	dgemm_("NoTranspose", "NoTranspose", nr, nr, n, &c_b33, &ti[ti_offset]
		, ldti, &t[t_offset], ldt, &c_b52, &dwork[ku], n, (ftnlen)11, 
		(ftnlen)11);
#line 609 "AB09BX.f"
	dgetrf_(nr, nr, &dwork[ku], n, &iwork[1], &ierr);
#line 610 "AB09BX.f"
	dgetrs_("NoTranspose", nr, n, &dwork[ku], n, &iwork[1], &ti[ti_offset]
		, ldti, &ierr, (ftnlen)11);

#line 613 "AB09BX.f"
	if (ns > 0) {
/*                           -1 */
/*           Compute (TI2*T2)  *TI2 in TI2. */

#line 617 "AB09BX.f"
	    dgemm_("NoTranspose", "NoTranspose", &ns, &ns, n, &c_b33, &ti[nr1 
		    + ti_dim1], ldti, &t[nr1 * t_dim1 + 1], ldt, &c_b52, &
		    dwork[ku], n, (ftnlen)11, (ftnlen)11);
#line 620 "AB09BX.f"
	    dgetrf_(&ns, &ns, &dwork[ku], n, &iwork[1], &ierr);
#line 621 "AB09BX.f"
	    dgetrs_("NoTranspose", &ns, n, &dwork[ku], n, &iwork[1], &ti[nr1 
		    + ti_dim1], ldti, &ierr, (ftnlen)11);
#line 623 "AB09BX.f"
	}
#line 624 "AB09BX.f"
    }

/*     Compute TI*A*T (A is in RSF). */

#line 628 "AB09BX.f"
    ij = ku;
#line 629 "AB09BX.f"
    i__1 = *n;
#line 629 "AB09BX.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 630 "AB09BX.f"
	i__2 = j + 1;
#line 630 "AB09BX.f"
	k = min(i__2,*n);
#line 631 "AB09BX.f"
	dgemv_("NoTranspose", &nminr, &k, &c_b33, &ti[ti_offset], ldti, &a[j *
		 a_dim1 + 1], &c__1, &c_b52, &dwork[ij], &c__1, (ftnlen)11);
#line 633 "AB09BX.f"
	ij += *n;
#line 634 "AB09BX.f"
/* L80: */
#line 634 "AB09BX.f"
    }
#line 635 "AB09BX.f"
    dgemm_("NoTranspose", "NoTranspose", &nminr, &nminr, n, &c_b33, &dwork[ku]
	    , n, &t[t_offset], ldt, &c_b52, &a[a_offset], lda, (ftnlen)11, (
	    ftnlen)11);

/*     Compute TI*B and C*T. */

#line 640 "AB09BX.f"
    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[ku], n, (ftnlen)4);
#line 641 "AB09BX.f"
    dgemm_("NoTranspose", "NoTranspose", &nminr, m, n, &c_b33, &ti[ti_offset],
	     ldti, &dwork[ku], n, &c_b52, &b[b_offset], ldb, (ftnlen)11, (
	    ftnlen)11);

#line 644 "AB09BX.f"
    dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[ku], p, (ftnlen)4);
#line 645 "AB09BX.f"
    dgemm_("NoTranspose", "NoTranspose", p, &nminr, n, &c_b33, &dwork[ku], p, 
	    &t[t_offset], ldt, &c_b52, &c__[c_offset], ldc, (ftnlen)11, (
	    ftnlen)11);

/*     Compute the singular perturbation approximation if possible. */
/*     Note that IERR = 1 on exit from AB09DD cannot appear here. */

/*     Workspace:  need real    4*(NMINR-NR); */
/*                 need integer 2*(NMINR-NR). */

#line 654 "AB09BX.f"
    ab09dd_(dico, &nminr, m, p, nr, &a[a_offset], lda, &b[b_offset], ldb, &
	    c__[c_offset], ldc, &d__[d_offset], ldd, &rcond, &iwork[1], &
	    dwork[1], &ierr, (ftnlen)1);

#line 657 "AB09BX.f"
    iwork[1] = nminr;
#line 658 "AB09BX.f"
    dwork[1] = (doublereal) wrkopt;

#line 660 "AB09BX.f"
    return 0;
/* *** Last line of AB09BX *** */
} /* ab09bx_ */
