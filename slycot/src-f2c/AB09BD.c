#line 1 "AB09BD.f"
/* AB09BD.f -- translated by f2c (version 20100827).
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

#line 1 "AB09BD.f"
/* Subroutine */ int ab09bd_(char *dico, char *job, char *equil, char *ordsel,
	 integer *n, integer *m, integer *p, integer *nr, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, doublereal *c__, integer *
	ldc, doublereal *d__, integer *ldd, doublereal *hsv, doublereal *tol1,
	 doublereal *tol2, integer *iwork, doublereal *dwork, integer *ldwork,
	 integer *iwarn, integer *info, ftnlen dico_len, ftnlen job_len, 
	ftnlen equil_len, ftnlen ordsel_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer ki, nn, kr, kt, kw, kti, ierr;
    extern /* Subroutine */ int tb01id_(char *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    ab09bx_(char *, char *, char *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01wd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *), xerbla_(char *, integer *, 
	    ftnlen);
    static doublereal maxred;
    static logical fixord;
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

/*     To compute a reduced order model (Ar,Br,Cr,Dr) for a stable */
/*     original state-space representation (A,B,C,D) by using either the */
/*     square-root or the balancing-free square-root Singular */
/*     Perturbation Approximation (SPA) model reduction method. */

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
/*             contain the state dynamics matrix A. */
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
/*             LDWORK >= MAX(1,N*(2*N+MAX(N,M,P)+5)+N*(N+1)/2). */
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
/*             = 1:  the reduction of A to the real Schur form failed; */
/*             = 2:  the state matrix A is not stable (if DICO = 'C') */
/*                   or not convergent (if DICO = 'D'); */
/*             = 3:  the computation of Hankel singular values failed. */

/*     METHOD */

/*     Let be the stable linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t) + Du(t)                           (1) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system. The subroutine AB09BD determines for */
/*     the given system (1), the matrices of a reduced order system */

/*          d[z(t)] = Ar*z(t) + Br*u(t) */
/*          yr(t)   = Cr*z(t) + Dr*u(t)                       (2) */

/*     such that */

/*           HSV(NR) <= INFNORM(G-Gr) <= 2*[HSV(NR+1) + ... + HSV(N)], */

/*     where G and Gr are transfer-function matrices of the systems */
/*     (A,B,C,D) and (Ar,Br,Cr,Dr), respectively, and INFNORM(G) is the */
/*     infinity-norm of G. */

/*     If JOB = 'B', the balancing-based square-root SPA method of [1] */
/*     is used and the resulting model is balanced. */

/*     If JOB = 'N', the balancing-free square-root SPA method of [2] */
/*     is used. */
/*     By setting TOL1 = TOL2, the routine can be used to compute */
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

/*     C. Oara and A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, March 1998. */
/*     Based on the RASP routine SRBFSP. */

/*     REVISIONS */

/*     May 2, 1998. */
/*     November 11, 1998, V. Sima, Research Institute for Informatics, */
/*     Bucharest. */

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

#line 282 "AB09BD.f"
    /* Parameter adjustments */
#line 282 "AB09BD.f"
    a_dim1 = *lda;
#line 282 "AB09BD.f"
    a_offset = 1 + a_dim1;
#line 282 "AB09BD.f"
    a -= a_offset;
#line 282 "AB09BD.f"
    b_dim1 = *ldb;
#line 282 "AB09BD.f"
    b_offset = 1 + b_dim1;
#line 282 "AB09BD.f"
    b -= b_offset;
#line 282 "AB09BD.f"
    c_dim1 = *ldc;
#line 282 "AB09BD.f"
    c_offset = 1 + c_dim1;
#line 282 "AB09BD.f"
    c__ -= c_offset;
#line 282 "AB09BD.f"
    d_dim1 = *ldd;
#line 282 "AB09BD.f"
    d_offset = 1 + d_dim1;
#line 282 "AB09BD.f"
    d__ -= d_offset;
#line 282 "AB09BD.f"
    --hsv;
#line 282 "AB09BD.f"
    --iwork;
#line 282 "AB09BD.f"
    --dwork;
#line 282 "AB09BD.f"

#line 282 "AB09BD.f"
    /* Function Body */
#line 282 "AB09BD.f"
    *info = 0;
#line 283 "AB09BD.f"
    *iwarn = 0;
#line 284 "AB09BD.f"
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 288 "AB09BD.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || lsame_(dico, "D", (
	    ftnlen)1, (ftnlen)1))) {
#line 289 "AB09BD.f"
	*info = -1;
#line 290 "AB09BD.f"
    } else if (! (lsame_(job, "B", (ftnlen)1, (ftnlen)1) || lsame_(job, "N", (
	    ftnlen)1, (ftnlen)1))) {
#line 291 "AB09BD.f"
	*info = -2;
#line 292 "AB09BD.f"
    } else if (! (lsame_(equil, "S", (ftnlen)1, (ftnlen)1) || lsame_(equil, 
	    "N", (ftnlen)1, (ftnlen)1))) {
#line 294 "AB09BD.f"
	*info = -3;
#line 295 "AB09BD.f"
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
#line 296 "AB09BD.f"
	*info = -4;
#line 297 "AB09BD.f"
    } else if (*n < 0) {
#line 298 "AB09BD.f"
	*info = -5;
#line 299 "AB09BD.f"
    } else if (*m < 0) {
#line 300 "AB09BD.f"
	*info = -6;
#line 301 "AB09BD.f"
    } else if (*p < 0) {
#line 302 "AB09BD.f"
	*info = -7;
#line 303 "AB09BD.f"
    } else if (fixord && (*nr < 0 || *nr > *n)) {
#line 304 "AB09BD.f"
	*info = -8;
#line 305 "AB09BD.f"
    } else if (*lda < max(1,*n)) {
#line 306 "AB09BD.f"
	*info = -10;
#line 307 "AB09BD.f"
    } else if (*ldb < max(1,*n)) {
#line 308 "AB09BD.f"
	*info = -12;
#line 309 "AB09BD.f"
    } else if (*ldc < max(1,*p)) {
#line 310 "AB09BD.f"
	*info = -14;
#line 311 "AB09BD.f"
    } else if (*ldd < max(1,*p)) {
#line 312 "AB09BD.f"
	*info = -16;
#line 313 "AB09BD.f"
    } else if (*tol2 > 0. && *tol2 > *tol1) {
#line 314 "AB09BD.f"
	*info = -19;
#line 315 "AB09BD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
#line 315 "AB09BD.f"
	i__3 = max(*n,*m);
#line 315 "AB09BD.f"
	i__1 = 1, i__2 = *n * ((*n << 1) + max(i__3,*p) + 5) + *n * (*n + 1) /
		 2;
#line 315 "AB09BD.f"
	if (*ldwork < max(i__1,i__2)) {
#line 317 "AB09BD.f"
	    *info = -22;
#line 318 "AB09BD.f"
	}
#line 318 "AB09BD.f"
    }

#line 320 "AB09BD.f"
    if (*info != 0) {

/*        Error return. */

#line 324 "AB09BD.f"
	i__1 = -(*info);
#line 324 "AB09BD.f"
	xerbla_("AB09BD", &i__1, (ftnlen)6);
#line 325 "AB09BD.f"
	return 0;
#line 326 "AB09BD.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 330 "AB09BD.f"
    i__1 = min(*n,*m);
#line 330 "AB09BD.f"
    if (min(i__1,*p) == 0) {
#line 331 "AB09BD.f"
	*nr = 0;
#line 332 "AB09BD.f"
	iwork[1] = 0;
#line 333 "AB09BD.f"
	dwork[1] = 1.;
#line 334 "AB09BD.f"
	return 0;
#line 335 "AB09BD.f"
    }

/*     Allocate working storage. */

#line 339 "AB09BD.f"
    nn = *n * *n;
#line 340 "AB09BD.f"
    kt = 1;
#line 341 "AB09BD.f"
    kr = kt + nn;
#line 342 "AB09BD.f"
    ki = kr + *n;
#line 343 "AB09BD.f"
    kw = ki + *n;

#line 345 "AB09BD.f"
    if (lsame_(equil, "S", (ftnlen)1, (ftnlen)1)) {

/*        Scale simultaneously the matrices A, B and C: */
/*        A <- inv(D)*A*D,  B <- inv(D)*B and C <- C*D, where D is a */
/*        diagonal matrix. */

#line 351 "AB09BD.f"
	maxred = 100.;
#line 352 "AB09BD.f"
	tb01id_("All", n, m, p, &maxred, &a[a_offset], lda, &b[b_offset], ldb,
		 &c__[c_offset], ldc, &dwork[1], info, (ftnlen)3);
#line 354 "AB09BD.f"
    }

/*     Reduce A to the real Schur form using an orthogonal similarity */
/*     transformation A <- T'*A*T and apply the transformation to */
/*     B and C: B <- T'*B and C <- C*T. */

#line 360 "AB09BD.f"
    i__1 = *ldwork - kw + 1;
#line 360 "AB09BD.f"
    tb01wd_(n, m, p, &a[a_offset], lda, &b[b_offset], ldb, &c__[c_offset], 
	    ldc, &dwork[kt], n, &dwork[kr], &dwork[ki], &dwork[kw], &i__1, &
	    ierr);
#line 362 "AB09BD.f"
    if (ierr != 0) {
#line 363 "AB09BD.f"
	*info = 1;
#line 364 "AB09BD.f"
	return 0;
#line 365 "AB09BD.f"
    }

#line 367 "AB09BD.f"
    wrkopt = dwork[kw] + (doublereal) (kw - 1);

#line 369 "AB09BD.f"
    kti = kt + nn;
#line 370 "AB09BD.f"
    kw = kti + nn;
#line 371 "AB09BD.f"
    i__1 = *ldwork - kw + 1;
#line 371 "AB09BD.f"
    ab09bx_(dico, job, ordsel, n, m, p, nr, &a[a_offset], lda, &b[b_offset], 
	    ldb, &c__[c_offset], ldc, &d__[d_offset], ldd, &hsv[1], &dwork[kt]
	    , n, &dwork[kti], n, tol1, tol2, &iwork[1], &dwork[kw], &i__1, 
	    iwarn, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 376 "AB09BD.f"
    if (ierr != 0) {
#line 377 "AB09BD.f"
	*info = ierr + 1;
#line 378 "AB09BD.f"
	return 0;
#line 379 "AB09BD.f"
    }

/* Computing MAX */
#line 381 "AB09BD.f"
    d__1 = wrkopt, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 381 "AB09BD.f"
    dwork[1] = max(d__1,d__2);

#line 383 "AB09BD.f"
    return 0;
/* *** Last line of AB09BD *** */
} /* ab09bd_ */

