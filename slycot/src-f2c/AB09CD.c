#line 1 "AB09CD.f"
/* AB09CD.f -- translated by f2c (version 20100827).
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

#line 1 "AB09CD.f"
/* Subroutine */ int ab09cd_(char *dico, char *equil, char *ordsel, integer *
	n, integer *m, integer *p, integer *nr, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *hsv, doublereal *tol1, 
	doublereal *tol2, integer *iwork, doublereal *dwork, integer *ldwork, 
	integer *iwarn, integer *info, ftnlen dico_len, ftnlen equil_len, 
	ftnlen ordsel_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer ki, kl, kt, kw, ierr;
    extern /* Subroutine */ int tb01id_(char *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    ab09cx_(char *, char *, integer *, integer *, integer *, integer *
	    , doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, ftnlen, ftnlen);
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
/*     original state-space representation (A,B,C,D) by using the */
/*     optimal Hankel-norm approximation method in conjunction with */
/*     square-root balancing. */

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
/*             reduced order model. NR is set as follows: */
/*             if ORDSEL = 'F', NR is equal to MIN(MAX(0,NR-KR+1),NMIN), */
/*             where KR is the multiplicity of the Hankel singular value */
/*             HSV(NR+1), NR is the desired order on entry, and NMIN is */
/*             the order of a minimal realization of the given system; */
/*             NMIN is determined as the number of Hankel singular values */
/*             greater than N*EPS*HNORM(A,B,C), where EPS is the machine */
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
/*             reduced order system in a real Schur form. */

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
/*             LDW1 = N*(2*N+MAX(N,M,P)+5) + N*(N+1)/2, */
/*             LDW2 = N*(M+P+2) + 2*M*P + MIN(N,M) + */
/*                    MAX( 3*M+1, MIN(N,M)+P ). */
/*             For optimum performance LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  with ORDSEL = 'F', the selected order NR is greater */
/*                   than the order of a minimal realization of the */
/*                   given system. In this case, the resulting NR is set */
/*                   automatically to a value corresponding to the order */
/*                   of a minimal realization of the system. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the reduction of A to the real Schur form failed; */
/*             = 2:  the state matrix A is not stable (if DICO = 'C') */
/*                   or not convergent (if DICO = 'D'); */
/*             = 3:  the computation of Hankel singular values failed; */
/*             = 4:  the computation of stable projection failed; */
/*             = 5:  the order of computed stable projection differs */
/*                   from the order of Hankel-norm approximation. */

/*     METHOD */

/*     Let be the stable linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t) + Du(t)                           (1) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system. The subroutine AB09CD determines for */
/*     the given system (1), the matrices of a reduced order system */

/*          d[z(t)] = Ar*z(t) + Br*u(t) */
/*          yr(t)   = Cr*z(t) + Dr*u(t)                       (2) */

/*     such that */

/*           HSV(NR) <= INFNORM(G-Gr) <= 2*[HSV(NR+1) + ... + HSV(N)], */

/*     where G and Gr are transfer-function matrices of the systems */
/*     (A,B,C,D) and (Ar,Br,Cr,Dr), respectively, and INFNORM(G) is the */
/*     infinity-norm of G. */

/*     The optimal Hankel-norm approximation method of [1], based on the */
/*     square-root balancing projection formulas of [2], is employed. */

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
/*     DLR Oberpfaffenhofen, April 1998. */
/*     Based on the RASP routine OHNAP. */

/*     REVISIONS */

/*     November 11, 1998, V. Sima, Research Institute for Informatics, */
/*     Bucharest. */
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

#line 279 "AB09CD.f"
    /* Parameter adjustments */
#line 279 "AB09CD.f"
    a_dim1 = *lda;
#line 279 "AB09CD.f"
    a_offset = 1 + a_dim1;
#line 279 "AB09CD.f"
    a -= a_offset;
#line 279 "AB09CD.f"
    b_dim1 = *ldb;
#line 279 "AB09CD.f"
    b_offset = 1 + b_dim1;
#line 279 "AB09CD.f"
    b -= b_offset;
#line 279 "AB09CD.f"
    c_dim1 = *ldc;
#line 279 "AB09CD.f"
    c_offset = 1 + c_dim1;
#line 279 "AB09CD.f"
    c__ -= c_offset;
#line 279 "AB09CD.f"
    d_dim1 = *ldd;
#line 279 "AB09CD.f"
    d_offset = 1 + d_dim1;
#line 279 "AB09CD.f"
    d__ -= d_offset;
#line 279 "AB09CD.f"
    --hsv;
#line 279 "AB09CD.f"
    --iwork;
#line 279 "AB09CD.f"
    --dwork;
#line 279 "AB09CD.f"

#line 279 "AB09CD.f"
    /* Function Body */
#line 279 "AB09CD.f"
    *info = 0;
#line 280 "AB09CD.f"
    *iwarn = 0;
#line 281 "AB09CD.f"
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);

/*     Check the input scalar arguments. */

#line 285 "AB09CD.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || lsame_(dico, "D", (
	    ftnlen)1, (ftnlen)1))) {
#line 286 "AB09CD.f"
	*info = -1;
#line 287 "AB09CD.f"
    } else if (! (lsame_(equil, "S", (ftnlen)1, (ftnlen)1) || lsame_(equil, 
	    "N", (ftnlen)1, (ftnlen)1))) {
#line 289 "AB09CD.f"
	*info = -2;
#line 290 "AB09CD.f"
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
#line 291 "AB09CD.f"
	*info = -3;
#line 292 "AB09CD.f"
    } else if (*n < 0) {
#line 293 "AB09CD.f"
	*info = -4;
#line 294 "AB09CD.f"
    } else if (*m < 0) {
#line 295 "AB09CD.f"
	*info = -5;
#line 296 "AB09CD.f"
    } else if (*p < 0) {
#line 297 "AB09CD.f"
	*info = -6;
#line 298 "AB09CD.f"
    } else if (fixord && (*nr < 0 || *nr > *n)) {
#line 299 "AB09CD.f"
	*info = -7;
#line 300 "AB09CD.f"
    } else if (*lda < max(1,*n)) {
#line 301 "AB09CD.f"
	*info = -9;
#line 302 "AB09CD.f"
    } else if (*ldb < max(1,*n)) {
#line 303 "AB09CD.f"
	*info = -11;
#line 304 "AB09CD.f"
    } else if (*ldc < max(1,*p)) {
#line 305 "AB09CD.f"
	*info = -13;
#line 306 "AB09CD.f"
    } else if (*ldd < max(1,*p)) {
#line 307 "AB09CD.f"
	*info = -15;
#line 308 "AB09CD.f"
    } else if (*tol2 > 0. && *tol2 > *tol1) {
#line 309 "AB09CD.f"
	*info = -18;
#line 310 "AB09CD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
#line 310 "AB09CD.f"
	i__3 = max(*n,*m);
/* Computing MAX */
#line 310 "AB09CD.f"
	i__4 = *m * 3 + 1, i__5 = min(*n,*m) + *p;
#line 310 "AB09CD.f"
	i__1 = *n * ((*n << 1) + max(i__3,*p) + 5) + *n * (*n + 1) / 2, i__2 =
		 *n * (*m + *p + 2) + (*m << 1) * *p + min(*n,*m) + max(i__4,
		i__5);
#line 310 "AB09CD.f"
	if (*ldwork < max(i__1,i__2)) {
#line 314 "AB09CD.f"
	    *info = -21;
#line 315 "AB09CD.f"
	}
#line 315 "AB09CD.f"
    }

#line 317 "AB09CD.f"
    if (*info != 0) {

/*        Error return. */

#line 321 "AB09CD.f"
	i__1 = -(*info);
#line 321 "AB09CD.f"
	xerbla_("AB09CD", &i__1, (ftnlen)6);
#line 322 "AB09CD.f"
	return 0;
#line 323 "AB09CD.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 327 "AB09CD.f"
    i__1 = min(*n,*m);
#line 327 "AB09CD.f"
    if (min(i__1,*p) == 0) {
#line 328 "AB09CD.f"
	*nr = 0;
#line 329 "AB09CD.f"
	iwork[1] = 0;
#line 330 "AB09CD.f"
	dwork[1] = 1.;
#line 331 "AB09CD.f"
	return 0;
#line 332 "AB09CD.f"
    }

#line 334 "AB09CD.f"
    if (lsame_(equil, "S", (ftnlen)1, (ftnlen)1)) {

/*        Scale simultaneously the matrices A, B and C: */
/*        A <- inv(D)*A*D, B <- inv(D)*B and C <- C*D, where D is a */
/*        diagonal matrix. */

#line 340 "AB09CD.f"
	maxred = 100.;
#line 341 "AB09CD.f"
	tb01id_("All", n, m, p, &maxred, &a[a_offset], lda, &b[b_offset], ldb,
		 &c__[c_offset], ldc, &dwork[1], info, (ftnlen)3);
#line 343 "AB09CD.f"
    }

/*     Reduce A to the real Schur form using an orthogonal similarity */
/*     transformation A <- T'*A*T and apply the transformation to B */
/*     and C: B <- T'*B and C <- C*T. */

#line 349 "AB09CD.f"
    kt = 1;
#line 350 "AB09CD.f"
    kl = kt + *n * *n;
#line 351 "AB09CD.f"
    ki = kl + *n;
#line 352 "AB09CD.f"
    kw = ki + *n;
#line 353 "AB09CD.f"
    i__1 = *ldwork - kw + 1;
#line 353 "AB09CD.f"
    tb01wd_(n, m, p, &a[a_offset], lda, &b[b_offset], ldb, &c__[c_offset], 
	    ldc, &dwork[kt], n, &dwork[kl], &dwork[ki], &dwork[kw], &i__1, &
	    ierr);
#line 355 "AB09CD.f"
    if (ierr != 0) {
#line 356 "AB09CD.f"
	*info = 1;
#line 357 "AB09CD.f"
	return 0;
#line 358 "AB09CD.f"
    }

#line 360 "AB09CD.f"
    wrkopt = dwork[kw] + (doublereal) (kw - 1);

#line 362 "AB09CD.f"
    ab09cx_(dico, ordsel, n, m, p, nr, &a[a_offset], lda, &b[b_offset], ldb, &
	    c__[c_offset], ldc, &d__[d_offset], ldd, &hsv[1], tol1, tol2, &
	    iwork[1], &dwork[1], ldwork, iwarn, &ierr, (ftnlen)1, (ftnlen)1);

#line 366 "AB09CD.f"
    if (ierr != 0) {
#line 367 "AB09CD.f"
	*info = ierr + 1;
#line 368 "AB09CD.f"
	return 0;
#line 369 "AB09CD.f"
    }

#line 371 "AB09CD.f"
    dwork[1] = max(wrkopt,dwork[1]);

#line 373 "AB09CD.f"
    return 0;
/* *** Last line of AB09CD *** */
} /* ab09cd_ */

