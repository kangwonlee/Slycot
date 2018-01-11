#line 1 "AB09AX.f"
/* AB09AX.f -- translated by f2c (version 20100827).
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

#line 1 "AB09AX.f"
/* Table of constant values */

static logical c_true = TRUE_;
static logical c_false = FALSE_;
static integer c__1 = 1;
static doublereal c_b31 = 1.;
static doublereal c_b50 = 0.;

/* Subroutine */ int ab09ax_(char *dico, char *job, char *ordsel, integer *n, 
	integer *m, integer *p, integer *nr, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *hsv, doublereal *t, integer *ldt, doublereal *ti, integer 
	*ldti, doublereal *tol, integer *iwork, doublereal *dwork, integer *
	ldwork, integer *iwarn, integer *info, ftnlen dico_len, ftnlen 
	job_len, ftnlen ordsel_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, t_dim1, 
	    t_offset, ti_dim1, ti_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, k, ij, ku, kv, kw;
    static logical bal;
    static integer ldw;
    static doublereal atol;
    static integer ierr, ktau;
    static doublereal temp, rtol;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    ma02dd_(char *, char *, integer *, doublereal *, integer *, 
	    doublereal *, ftnlen, ftnlen), dscal_(integer *, doublereal *, 
	    doublereal *, integer *), dgemm_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen), mb03ud_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static logical discr;
    extern /* Subroutine */ int sb03ou_(logical *, logical *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, integer *), dtrmm_(char *, char *, char *, char *, 
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

/*     To compute a reduced order model (Ar,Br,Cr) for a stable original */
/*     state-space representation (A,B,C) by using either the square-root */
/*     or the balancing-free square-root Balance & Truncate model */
/*     reduction method. The state dynamics matrix A of the original */
/*     system is an upper quasi-triangular matrix in real Schur canonical */
/*     form. The matrices of the reduced order system are computed using */
/*     the truncation formulas: */

/*          Ar = TI * A * T ,  Br = TI * B ,  Cr = C * T . */

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
/*             singular values greater than MAX(TOL,N*EPS*HNORM(A,B,C)). */

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

/*     HSV     (output) DOUBLE PRECISION array, dimension (N) */
/*             If INFO = 0, it contains the Hankel singular values of */
/*             the original system ordered decreasingly. HSV(1) is the */
/*             Hankel norm of the system. */

/*     T       (output) DOUBLE PRECISION array, dimension (LDT,N) */
/*             If INFO = 0 and NR > 0, the leading N-by-NR part of this */
/*             array contains the right truncation matrix T. */

/*     LDT     INTEGER */
/*             The leading dimension of array T.  LDT >= MAX(1,N). */

/*     TI      (output) DOUBLE PRECISION array, dimension (LDTI,N) */
/*             If INFO = 0 and NR > 0, the leading NR-by-N part of this */
/*             array contains the left truncation matrix TI. */

/*     LDTI    INTEGER */
/*             The leading dimension of array TI.  LDTI >= MAX(1,N). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If ORDSEL = 'A', TOL contains the tolerance for */
/*             determining the order of reduced system. */
/*             For model reduction, the recommended value is */
/*             TOL = c*HNORM(A,B,C), where c is a constant in the */
/*             interval [0.00001,0.001], and HNORM(A,B,C) is the */
/*             Hankel-norm of the given system (computed in HSV(1)). */
/*             For computing a minimal realization, the recommended */
/*             value is TOL = N*EPS*HNORM(A,B,C), where EPS is the */
/*             machine precision (see LAPACK Library Routine DLAMCH). */
/*             This value is used by default if TOL <= 0 on entry. */
/*             If ORDSEL = 'F', the value of TOL is ignored. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK = 0, if JOB = 'B', or */
/*             LIWORK = N, if JOB = 'N'. */

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
/*          y(t)    = Cx(t)                               (1) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system. The subroutine AB09AX determines for */
/*     the given system (1), the matrices of a reduced NR order system */

/*          d[z(t)] = Ar*z(t) + Br*u(t) */
/*          yr(t)   = Cr*z(t)                             (2) */

/*     such that */

/*           HSV(NR) <= INFNORM(G-Gr) <= 2*[HSV(NR+1) + ... + HSV(N)], */

/*     where G and Gr are transfer-function matrices of the systems */
/*     (A,B,C) and (Ar,Br,Cr), respectively, and INFNORM(G) is the */
/*     infinity-norm of G. */

/*     If JOB = 'B', the square-root Balance & Truncate method of [1] */
/*     is used and, for DICO = 'C', the resulting model is balanced. */
/*     By setting TOL <= 0, the routine can be used to compute balanced */
/*     minimal state-space realizations of stable systems. */

/*     If JOB = 'N', the balancing-free square-root version of the */
/*     Balance & Truncate method [2] is used. */
/*     By setting TOL <= 0, the routine can be used to compute minimal */
/*     state-space realizations of stable systems. */

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

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, March 1998. */
/*     Based on the RASP routines SRBT1 and SRBFT1. */

/*     REVISIONS */

/*     May 2, 1998. */
/*     November 11, 1998, V. Sima, Research Institute for Informatics, */
/*     Bucharest. */
/*     December 1998, V. Sima, Katholieke Univ. Leuven, Leuven. */
/*     February 14, 1999, A. Varga, German Aerospace Center. */
/*     February 22, 1999, V. Sima, Research Institute for Informatics. */
/*     February 27, 2000, V. Sima, Research Institute for Informatics. */

/*     KEYWORDS */

/*     Balancing, minimal state-space representation, model reduction, */
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

#line 288 "AB09AX.f"
    /* Parameter adjustments */
#line 288 "AB09AX.f"
    a_dim1 = *lda;
#line 288 "AB09AX.f"
    a_offset = 1 + a_dim1;
#line 288 "AB09AX.f"
    a -= a_offset;
#line 288 "AB09AX.f"
    b_dim1 = *ldb;
#line 288 "AB09AX.f"
    b_offset = 1 + b_dim1;
#line 288 "AB09AX.f"
    b -= b_offset;
#line 288 "AB09AX.f"
    c_dim1 = *ldc;
#line 288 "AB09AX.f"
    c_offset = 1 + c_dim1;
#line 288 "AB09AX.f"
    c__ -= c_offset;
#line 288 "AB09AX.f"
    --hsv;
#line 288 "AB09AX.f"
    t_dim1 = *ldt;
#line 288 "AB09AX.f"
    t_offset = 1 + t_dim1;
#line 288 "AB09AX.f"
    t -= t_offset;
#line 288 "AB09AX.f"
    ti_dim1 = *ldti;
#line 288 "AB09AX.f"
    ti_offset = 1 + ti_dim1;
#line 288 "AB09AX.f"
    ti -= ti_offset;
#line 288 "AB09AX.f"
    --iwork;
#line 288 "AB09AX.f"
    --dwork;
#line 288 "AB09AX.f"

#line 288 "AB09AX.f"
    /* Function Body */
#line 288 "AB09AX.f"
    *info = 0;
#line 289 "AB09AX.f"
    *iwarn = 0;
#line 290 "AB09AX.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 291 "AB09AX.f"
    bal = lsame_(job, "B", (ftnlen)1, (ftnlen)1);
#line 292 "AB09AX.f"
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 296 "AB09AX.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 297 "AB09AX.f"
	*info = -1;
#line 298 "AB09AX.f"
    } else if (! (bal || lsame_(job, "N", (ftnlen)1, (ftnlen)1))) {
#line 299 "AB09AX.f"
	*info = -2;
#line 300 "AB09AX.f"
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
#line 301 "AB09AX.f"
	*info = -3;
#line 302 "AB09AX.f"
    } else if (*n < 0) {
#line 303 "AB09AX.f"
	*info = -4;
#line 304 "AB09AX.f"
    } else if (*m < 0) {
#line 305 "AB09AX.f"
	*info = -5;
#line 306 "AB09AX.f"
    } else if (*p < 0) {
#line 307 "AB09AX.f"
	*info = -6;
#line 308 "AB09AX.f"
    } else if (fixord && (*nr < 0 || *nr > *n)) {
#line 309 "AB09AX.f"
	*info = -7;
#line 310 "AB09AX.f"
    } else if (*lda < max(1,*n)) {
#line 311 "AB09AX.f"
	*info = -9;
#line 312 "AB09AX.f"
    } else if (*ldb < max(1,*n)) {
#line 313 "AB09AX.f"
	*info = -11;
#line 314 "AB09AX.f"
    } else if (*ldc < max(1,*p)) {
#line 315 "AB09AX.f"
	*info = -13;
#line 316 "AB09AX.f"
    } else if (*ldt < max(1,*n)) {
#line 317 "AB09AX.f"
	*info = -16;
#line 318 "AB09AX.f"
    } else if (*ldti < max(1,*n)) {
#line 319 "AB09AX.f"
	*info = -18;
#line 320 "AB09AX.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
#line 320 "AB09AX.f"
	i__3 = max(*n,*m);
#line 320 "AB09AX.f"
	i__1 = 1, i__2 = *n * (max(i__3,*p) + 5) + *n * (*n + 1) / 2;
#line 320 "AB09AX.f"
	if (*ldwork < max(i__1,i__2)) {
#line 322 "AB09AX.f"
	    *info = -22;
#line 323 "AB09AX.f"
	}
#line 323 "AB09AX.f"
    }

#line 325 "AB09AX.f"
    if (*info != 0) {

/*        Error return. */

#line 329 "AB09AX.f"
	i__1 = -(*info);
#line 329 "AB09AX.f"
	xerbla_("AB09AX", &i__1, (ftnlen)6);
#line 330 "AB09AX.f"
	return 0;
#line 331 "AB09AX.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 335 "AB09AX.f"
    i__1 = min(*n,*m);
#line 335 "AB09AX.f"
    if (min(i__1,*p) == 0 || fixord && *nr == 0) {
#line 336 "AB09AX.f"
	*nr = 0;
#line 337 "AB09AX.f"
	dwork[1] = 1.;
#line 338 "AB09AX.f"
	return 0;
#line 339 "AB09AX.f"
    }

#line 341 "AB09AX.f"
    rtol = (doublereal) (*n) * dlamch_("Epsilon", (ftnlen)7);

/*     Allocate N*MAX(N,M,P) and N working storage for the matrices U */
/*     and TAU, respectively. */

#line 346 "AB09AX.f"
    ku = 1;
/* Computing MAX */
#line 347 "AB09AX.f"
    i__1 = max(*n,*m);
#line 347 "AB09AX.f"
    ktau = ku + *n * max(i__1,*p);
#line 348 "AB09AX.f"
    kw = ktau + *n;
#line 349 "AB09AX.f"
    ldw = *ldwork - kw + 1;

/*     Copy B in U. */

#line 353 "AB09AX.f"
    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[ku], n, (ftnlen)4);

/*     If DISCR = .FALSE., solve for Su the Lyapunov equation */
/*                                      2 */
/*     A*(Su*Su') + (Su*Su')*A' + scalec *B*B' = 0 . */

/*     If DISCR = .TRUE., solve for Su the Lyapunov equation */
/*                           2 */
/*     A*(Su*Su')*A' + scalec *B*B' = Su*Su' . */

/*     Workspace:  need   N*(MAX(N,M,P) + 5); */
/*                 prefer larger. */

#line 366 "AB09AX.f"
    sb03ou_(&discr, &c_true, n, m, &a[a_offset], lda, &dwork[ku], n, &dwork[
	    ktau], &ti[ti_offset], ldti, &scalec, &dwork[kw], &ldw, &ierr);
#line 368 "AB09AX.f"
    if (ierr != 0) {
#line 369 "AB09AX.f"
	*info = 1;
#line 370 "AB09AX.f"
	return 0;
#line 371 "AB09AX.f"
    }
#line 372 "AB09AX.f"
    wrkopt = (integer) dwork[kw] + kw - 1;

/*     Copy C in U. */

#line 376 "AB09AX.f"
    dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[ku], p, (ftnlen)4);

/*     If DISCR = .FALSE., solve for Ru the Lyapunov equation */
/*                                      2 */
/*     A'*(Ru'*Ru) + (Ru'*Ru)*A + scaleo  * C'*C = 0 . */

/*     If DISCR = .TRUE., solve for Ru the Lyapunov equation */
/*                           2 */
/*     A'*(Ru'*Ru)*A + scaleo  * C'*C = Ru'*Ru . */

/*     Workspace:  need   N*(MAX(N,M,P) + 5); */
/*                 prefer larger. */

#line 389 "AB09AX.f"
    sb03ou_(&discr, &c_false, n, p, &a[a_offset], lda, &dwork[ku], p, &dwork[
	    ktau], &t[t_offset], ldt, &scaleo, &dwork[kw], &ldw, &ierr);
/* Computing MAX */
#line 391 "AB09AX.f"
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 391 "AB09AX.f"
    wrkopt = max(i__1,i__2);

/*     Allocate N*(N+1)/2 (or, if possible, N*N) working storage for the */
/*     matrix V, a packed (or unpacked) copy of Su, and save Su in V. */
/*     (The locations for TAU are reused here.) */

#line 397 "AB09AX.f"
    kv = ktau;
#line 398 "AB09AX.f"
    if (*ldwork - kv + 1 < *n * (*n + 5)) {
#line 399 "AB09AX.f"
	packed = TRUE_;
#line 400 "AB09AX.f"
	ma02dd_("Pack", "Upper", n, &ti[ti_offset], ldti, &dwork[kv], (ftnlen)
		4, (ftnlen)5);
#line 401 "AB09AX.f"
	kw = kv + *n * (*n + 1) / 2;
#line 402 "AB09AX.f"
    } else {
#line 403 "AB09AX.f"
	packed = FALSE_;
#line 404 "AB09AX.f"
	dlacpy_("Upper", n, n, &ti[ti_offset], ldti, &dwork[kv], n, (ftnlen)5)
		;
#line 405 "AB09AX.f"
	kw = kv + *n * *n;
#line 406 "AB09AX.f"
    }
/*                               | x x | */
/*     Compute Ru*Su in the form | 0 x | in TI. */

#line 410 "AB09AX.f"
    i__1 = *n;
#line 410 "AB09AX.f"
    for (j = 1; j <= i__1; ++j) {
#line 411 "AB09AX.f"
	dtrmv_("Upper", "NoTranspose", "NonUnit", &j, &t[t_offset], ldt, &ti[
		j * ti_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)11, (ftnlen)7);
#line 413 "AB09AX.f"
/* L10: */
#line 413 "AB09AX.f"
    }

/*     Compute the singular value decomposition Ru*Su = V*S*UT */
/*     of the upper triangular matrix Ru*Su, with UT in TI and V in U. */

/*     Workspace:  need   N*MAX(N,M,P) + N*(N+1)/2 + 5*N; */
/*                 prefer larger. */

#line 421 "AB09AX.f"
    i__1 = *ldwork - kw + 1;
#line 421 "AB09AX.f"
    mb03ud_("Vectors", "Vectors", n, &ti[ti_offset], ldti, &dwork[ku], n, &
	    hsv[1], &dwork[kw], &i__1, &ierr, (ftnlen)7, (ftnlen)7);
#line 423 "AB09AX.f"
    if (ierr != 0) {
#line 424 "AB09AX.f"
	*info = 2;
#line 425 "AB09AX.f"
	return 0;
#line 426 "AB09AX.f"
    }
/* Computing MAX */
#line 427 "AB09AX.f"
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 427 "AB09AX.f"
    wrkopt = max(i__1,i__2);

/*     Scale singular values. */

#line 431 "AB09AX.f"
    d__1 = 1. / scalec / scaleo;
#line 431 "AB09AX.f"
    dscal_(n, &d__1, &hsv[1], &c__1);

/*     Partition S, U and V conformally as: */

/*     S = diag(S1,S2),  U = [U1,U2] (U' in TI) and V = [V1,V2] (in U). */

/*     Compute the order of reduced system, as the order of S1. */

#line 439 "AB09AX.f"
    atol = rtol * hsv[1];
#line 440 "AB09AX.f"
    if (fixord) {
#line 441 "AB09AX.f"
	if (*nr > 0) {
#line 442 "AB09AX.f"
	    if (hsv[*nr] <= atol) {
#line 443 "AB09AX.f"
		*nr = 0;
#line 444 "AB09AX.f"
		*iwarn = 1;
#line 445 "AB09AX.f"
		fixord = FALSE_;
#line 446 "AB09AX.f"
	    }
#line 447 "AB09AX.f"
	}
#line 448 "AB09AX.f"
    } else {
#line 449 "AB09AX.f"
	atol = max(*tol,atol);
#line 450 "AB09AX.f"
	*nr = 0;
#line 451 "AB09AX.f"
    }
#line 452 "AB09AX.f"
    if (! fixord) {
#line 453 "AB09AX.f"
	i__1 = *n;
#line 453 "AB09AX.f"
	for (j = 1; j <= i__1; ++j) {
#line 454 "AB09AX.f"
	    if (hsv[j] <= atol) {
#line 454 "AB09AX.f"
		goto L30;
#line 454 "AB09AX.f"
	    }
#line 455 "AB09AX.f"
	    ++(*nr);
#line 456 "AB09AX.f"
/* L20: */
#line 456 "AB09AX.f"
	}
#line 457 "AB09AX.f"
L30:
#line 458 "AB09AX.f"
	;
#line 458 "AB09AX.f"
    }

#line 460 "AB09AX.f"
    if (*nr == 0) {
#line 461 "AB09AX.f"
	dwork[1] = (doublereal) wrkopt;
#line 462 "AB09AX.f"
	return 0;
#line 463 "AB09AX.f"
    }

/*     Compute the truncation matrices. */

/*     Compute TI' =  Ru'*V1 in U. */

#line 469 "AB09AX.f"
    dtrmm_("Left", "Upper", "Transpose", "NonUnit", n, nr, &c_b31, &t[
	    t_offset], ldt, &dwork[ku], n, (ftnlen)4, (ftnlen)5, (ftnlen)9, (
	    ftnlen)7);

/*     Compute T = Su*U1 (with Su packed, if not enough workspace). */

#line 474 "AB09AX.f"
    ma02ad_("Full", nr, n, &ti[ti_offset], ldti, &t[t_offset], ldt, (ftnlen)4)
	    ;
#line 475 "AB09AX.f"
    if (packed) {
#line 476 "AB09AX.f"
	i__1 = *nr;
#line 476 "AB09AX.f"
	for (j = 1; j <= i__1; ++j) {
#line 477 "AB09AX.f"
	    dtpmv_("Upper", "NoTranspose", "NonUnit", n, &dwork[kv], &t[j * 
		    t_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)11, (ftnlen)7);
#line 479 "AB09AX.f"
/* L40: */
#line 479 "AB09AX.f"
	}
#line 480 "AB09AX.f"
    } else {
#line 481 "AB09AX.f"
	dtrmm_("Left", "Upper", "NoTranspose", "NonUnit", n, nr, &c_b31, &
		dwork[kv], n, &t[t_offset], ldt, (ftnlen)4, (ftnlen)5, (
		ftnlen)11, (ftnlen)7);
#line 483 "AB09AX.f"
    }

#line 485 "AB09AX.f"
    if (bal) {
#line 486 "AB09AX.f"
	ij = ku;

/*        Square-Root B & T method. */

/*        Compute the truncation matrices for balancing */
/*                    -1/2           -1/2 */
/*                T*S1     and TI'*S1 */

#line 494 "AB09AX.f"
	i__1 = *nr;
#line 494 "AB09AX.f"
	for (j = 1; j <= i__1; ++j) {
#line 495 "AB09AX.f"
	    temp = 1. / sqrt(hsv[j]);
#line 496 "AB09AX.f"
	    dscal_(n, &temp, &t[j * t_dim1 + 1], &c__1);
#line 497 "AB09AX.f"
	    dscal_(n, &temp, &dwork[ij], &c__1);
#line 498 "AB09AX.f"
	    ij += *n;
#line 499 "AB09AX.f"
/* L50: */
#line 499 "AB09AX.f"
	}
#line 500 "AB09AX.f"
    } else {

/*        Balancing-Free B & T method. */

/*        Compute orthogonal bases for the images of matrices T and TI'. */

/*        Workspace:  need   N*MAX(N,M,P) + 2*NR; */
/*                    prefer N*MAX(N,M,P) + NR*(NB+1) */
/*                           (NB determined by ILAENV for DGEQRF). */

#line 510 "AB09AX.f"
	kw = ktau + *nr;
#line 511 "AB09AX.f"
	ldw = *ldwork - kw + 1;
#line 512 "AB09AX.f"
	dgeqrf_(n, nr, &t[t_offset], ldt, &dwork[ktau], &dwork[kw], &ldw, &
		ierr);
#line 513 "AB09AX.f"
	dorgqr_(n, nr, nr, &t[t_offset], ldt, &dwork[ktau], &dwork[kw], &ldw, 
		&ierr);
#line 515 "AB09AX.f"
	dgeqrf_(n, nr, &dwork[ku], n, &dwork[ktau], &dwork[kw], &ldw, &ierr);
/* Computing MAX */
#line 517 "AB09AX.f"
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 517 "AB09AX.f"
	wrkopt = max(i__1,i__2);
#line 518 "AB09AX.f"
	dorgqr_(n, nr, nr, &dwork[ku], n, &dwork[ktau], &dwork[kw], &ldw, &
		ierr);
/* Computing MAX */
#line 520 "AB09AX.f"
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 520 "AB09AX.f"
	wrkopt = max(i__1,i__2);
#line 521 "AB09AX.f"
    }

/*     Transpose TI' to obtain TI. */

#line 525 "AB09AX.f"
    ma02ad_("Full", n, nr, &dwork[ku], n, &ti[ti_offset], ldti, (ftnlen)4);

#line 527 "AB09AX.f"
    if (! bal) {
/*                      -1 */
/*        Compute (TI*T)  *TI in TI. */

#line 531 "AB09AX.f"
	dgemm_("NoTranspose", "NoTranspose", nr, nr, n, &c_b31, &ti[ti_offset]
		, ldti, &t[t_offset], ldt, &c_b50, &dwork[ku], n, (ftnlen)11, 
		(ftnlen)11);
#line 533 "AB09AX.f"
	dgetrf_(nr, nr, &dwork[ku], n, &iwork[1], &ierr);
#line 534 "AB09AX.f"
	dgetrs_("NoTranspose", nr, n, &dwork[ku], n, &iwork[1], &ti[ti_offset]
		, ldti, &ierr, (ftnlen)11);
#line 536 "AB09AX.f"
    }

/*     Compute TI*A*T (A is in RSF). */

#line 540 "AB09AX.f"
    ij = ku;
#line 541 "AB09AX.f"
    i__1 = *n;
#line 541 "AB09AX.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 542 "AB09AX.f"
	i__2 = j + 1;
#line 542 "AB09AX.f"
	k = min(i__2,*n);
#line 543 "AB09AX.f"
	dgemv_("NoTranspose", nr, &k, &c_b31, &ti[ti_offset], ldti, &a[j * 
		a_dim1 + 1], &c__1, &c_b50, &dwork[ij], &c__1, (ftnlen)11);
#line 545 "AB09AX.f"
	ij += *n;
#line 546 "AB09AX.f"
/* L60: */
#line 546 "AB09AX.f"
    }
#line 547 "AB09AX.f"
    dgemm_("NoTranspose", "NoTranspose", nr, nr, n, &c_b31, &dwork[ku], n, &t[
	    t_offset], ldt, &c_b50, &a[a_offset], lda, (ftnlen)11, (ftnlen)11)
	    ;

/*     Compute TI*B and C*T. */

#line 552 "AB09AX.f"
    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[ku], n, (ftnlen)4);
#line 553 "AB09AX.f"
    dgemm_("NoTranspose", "NoTranspose", nr, m, n, &c_b31, &ti[ti_offset], 
	    ldti, &dwork[ku], n, &c_b50, &b[b_offset], ldb, (ftnlen)11, (
	    ftnlen)11);

#line 556 "AB09AX.f"
    dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[ku], p, (ftnlen)4);
#line 557 "AB09AX.f"
    dgemm_("NoTranspose", "NoTranspose", p, nr, n, &c_b31, &dwork[ku], p, &t[
	    t_offset], ldt, &c_b50, &c__[c_offset], ldc, (ftnlen)11, (ftnlen)
	    11);

#line 560 "AB09AX.f"
    dwork[1] = (doublereal) wrkopt;

#line 562 "AB09AX.f"
    return 0;
/* *** Last line of AB09AX *** */
} /* ab09ax_ */

