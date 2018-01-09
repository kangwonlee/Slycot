#line 1 "SB16CD.f"
/* SB16CD.f -- translated by f2c (version 20100827).
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

#line 1 "SB16CD.f"
/* Table of constant values */

static doublereal c_b17 = 1.;
static doublereal c_b18 = 0.;

/* Subroutine */ int sb16cd_(char *dico, char *jobd, char *jobmr, char *jobcf,
	 char *ordsel, integer *n, integer *m, integer *p, integer *ncr, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *d__, integer *ldd, doublereal *f, 
	integer *ldf, doublereal *g, integer *ldg, doublereal *hsv, 
	doublereal *tol, integer *iwork, doublereal *dwork, integer *ldwork, 
	integer *iwarn, integer *info, ftnlen dico_len, ftnlen jobd_len, 
	ftnlen jobmr_len, ftnlen jobcf_len, ftnlen ordsel_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, f_dim1, f_offset, g_dim1, g_offset, i__1, i__2;

    /* Local variables */
    static integer mp, kt, kw, lw;
    static logical bal;
    static integer kti, nmr;
    static logical left;
    static integer ierr;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     ab09ix_(char *, char *, char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    extern /* Subroutine */ int sb16cy_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static logical withd;
    static doublereal scalec, scaleo;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical fixord;
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

/*     To compute, for a given open-loop model (A,B,C,D), and for */
/*     given state feedback gain F and full observer gain G, */
/*     such that A+B*F and A+G*C are stable, a reduced order */
/*     controller model (Ac,Bc,Cc) using a coprime factorization */
/*     based controller reduction approach. For reduction of */
/*     coprime factors, a stability enforcing frequency-weighted */
/*     model reduction is performed using either the square-root or */
/*     the balancing-free square-root versions of the Balance & Truncate */
/*     (B&T) model reduction method. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the open-loop system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     JOBD    CHARACTER*1 */
/*             Specifies whether or not a non-zero matrix D appears */
/*             in the given state space model, as follows: */
/*             = 'D':  D is present; */
/*             = 'Z':  D is assumed a zero matrix. */

/*     JOBMR   CHARACTER*1 */
/*             Specifies the model reduction approach to be used */
/*             as follows: */
/*             = 'B':  use the square-root B&T method; */
/*             = 'F':  use the balancing-free square-root B&T method. */

/*     JOBCF   CHARACTER*1 */
/*             Specifies whether left or right coprime factorization */
/*             of the controller is to be used as follows: */
/*             = 'L':  use left coprime factorization; */
/*             = 'R':  use right coprime factorization. */

/*     ORDSEL  CHARACTER*1 */
/*             Specifies the order selection method as follows: */
/*             = 'F':  the resulting controller order NCR is fixed; */
/*             = 'A':  the resulting controller order NCR is */
/*                     automatically determined on basis of the given */
/*                     tolerance TOL. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the original state-space representation, i.e. */
/*             the order of the matrix A.  N >= 0. */
/*             N also represents the order of the original state-feedback */
/*             controller. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     NCR     (input/output) INTEGER */
/*             On entry with ORDSEL = 'F', NCR is the desired order of */
/*             the resulting reduced order controller.  0 <= NCR <= N. */
/*             On exit, if INFO = 0, NCR is the order of the resulting */
/*             reduced order controller. NCR is set as follows: */
/*             if ORDSEL = 'F', NCR is equal to MIN(NCR,NCRMIN), where */
/*             NCR is the desired order on entry, and NCRMIN is the */
/*             number of Hankel-singular values greater than N*EPS*S1, */
/*             where EPS is the machine precision (see LAPACK Library */
/*             Routine DLAMCH) and S1 is the largest Hankel singular */
/*             value (computed in HSV(1)); NCR can be further reduced */
/*             to ensure HSV(NCR) > HSV(NCR+1); */
/*             if ORDSEL = 'A', NCR is equal to the number of Hankel */
/*             singular values greater than MAX(TOL,N*EPS*S1). */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original state dynamics matrix A. */
/*             On exit, if INFO = 0, the leading NCR-by-NCR part of this */
/*             array contains the state dynamics matrix Ac of the reduced */
/*             controller. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the open-loop system input/state matrix B. */
/*             On exit, this array is overwritten with a NCR-by-M */
/*             B&T approximation of the matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the open-loop system state/output matrix C. */
/*             On exit, this array is overwritten with a P-by-NCR */
/*             B&T approximation of the matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, if JOBD = 'D', the leading P-by-M part of this */
/*             array must contain the system direct input/output */
/*             transmission matrix D. */
/*             The array D is not referenced if JOBD = 'Z'. */

/*     LDD     INTEGER */
/*             The leading dimension of array D. */
/*             LDD >= MAX(1,P), if JOBD = 'D'; */
/*             LDD >= 1,        if JOBD = 'Z'. */

/*     F       (input/output) DOUBLE PRECISION array, dimension (LDF,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain a stabilizing state feedback matrix. */
/*             On exit, if INFO = 0, the leading M-by-NCR part of this */
/*             array contains the output/state matrix Cc of the reduced */
/*             controller. */

/*     LDF     INTEGER */
/*             The leading dimension of array F.  LDF >= MAX(1,M). */

/*     G       (input/output) DOUBLE PRECISION array, dimension (LDG,P) */
/*             On entry, the leading N-by-P part of this array must */
/*             contain a stabilizing observer gain matrix. */
/*             On exit, if INFO = 0, the leading NCR-by-P part of this */
/*             array contains the input/state matrix Bc of the reduced */
/*             controller. */

/*     LDG     INTEGER */
/*             The leading dimension of array G.  LDG >= MAX(1,N). */

/*     HSV     (output) DOUBLE PRECISION array, dimension (N) */
/*             If INFO = 0, HSV contains the N frequency-weighted */
/*             Hankel singular values ordered decreasingly (see METHOD). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If ORDSEL = 'A', TOL contains the tolerance for */
/*             determining the order of reduced controller. */
/*             The recommended value is TOL = c*S1, where c is a constant */
/*             in the interval [0.00001,0.001], and S1 is the largest */
/*             Hankel singular value (computed in HSV(1)). */
/*             The value TOL = N*EPS*S1 is used by default if */
/*             TOL <= 0 on entry, where EPS is the machine precision */
/*             (see LAPACK Library Routine DLAMCH). */
/*             If ORDSEL = 'F', the value of TOL is ignored. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension LIWORK, where */
/*             LIWORK = 0,   if JOBMR = 'B'; */
/*             LIWORK = N,   if JOBMR = 'F'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 2*N*N + MAX( 1, 2*N*N + 5*N, N*MAX(M,P), */
/*                                    N*(N + MAX(N,MP) + MIN(N,MP) + 6)), */
/*             where     MP = M, if JOBCF = 'L'; */
/*                       MP = P, if JOBCF = 'R'. */
/*             For optimum performance LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  with ORDSEL = 'F', the selected order NCR is */
/*                   greater than the order of a minimal realization */
/*                   of the controller; */
/*             = 2:  with ORDSEL = 'F', the selected order NCR */
/*                   corresponds to repeated singular values, which are */
/*                   neither all included nor all excluded from the */
/*                   reduced controller. In this case, the resulting NCR */
/*                   is set automatically to the largest value such that */
/*                   HSV(NCR) > HSV(NCR+1). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  eigenvalue computation failure; */
/*             = 2:  the matrix A+G*C is not stable; */
/*             = 3:  the matrix A+B*F is not stable; */
/*             = 4:  the Lyapunov equation for computing the */
/*                   observability Grammian is (nearly) singular; */
/*             = 5:  the Lyapunov equation for computing the */
/*                   controllability Grammian is (nearly) singular; */
/*             = 6:  the computation of Hankel singular values failed. */

/*     METHOD */

/*     Let be the linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t) + Du(t),                             (1) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system, and let Go(d) be the open-loop */
/*     transfer-function matrix */
/*                          -1 */
/*          Go(d) = C*(d*I-A) *B + D . */

/*     Let F and G be the state feedback and observer gain matrices, */
/*     respectively, chosen such that A+BF and A+GC are stable matrices. */
/*     The controller has a transfer-function matrix K(d) given by */
/*                                       -1 */
/*          K(d) = F*(d*I-A-B*F-G*C-G*D*F) *G . */

/*     The closed-loop transfer function matrix is given by */
/*                                    -1 */
/*          Gcl(d) = Go(d)(I+K(d)Go(d)) . */

/*     K(d) can be expressed as a left coprime factorization (LCF) */
/*                         -1 */
/*          K(d) = M_left(d) *N_left(d), */

/*     or as a right coprime factorization (RCF) */
/*                                     -1 */
/*          K(d) = N_right(d)*M_right(d) , */

/*     where M_left(d), N_left(d), N_right(d), and M_right(d) are */
/*     stable transfer-function matrices. */

/*     The subroutine SB16CD determines the matrices of a reduced */
/*     controller */

/*          d[z(t)] = Ac*z(t) + Bc*y(t) */
/*          u(t)    = Cc*z(t),                                   (2) */

/*     with the transfer-function matrix Kr, using the following */
/*     stability enforcing approach proposed in [1]: */

/*     (1) If JOBCF = 'L', the frequency-weighted approximation problem */
/*         is solved */

/*         min||[M_left(d)-M_leftr(d)  N_left(d)-N_leftr(d)][-Y(d)]|| , */
/*                                                          [ X(d)] */
/*         where */
/*                              -1 */
/*               G(d) = Y(d)*X(d) */

/*         is a RCF of the open-loop system transfer-function matrix. */
/*         The B&T model reduction technique is used in conjunction */
/*         with the method proposed in [1]. */

/*     (2) If JOBCF = 'R', the frequency-weighted approximation problem */
/*         is solved */

/*         min || [ -U(d) V(d) ] [ N_right(d)-N_rightr(d) ] || , */
/*                               [ M_right(d)-M_rightr(d) ] */
/*         where */
/*                         -1 */
/*               G(d) = V(d) *U(d) */

/*         is a LCF of the open-loop system transfer-function matrix. */
/*         The B&T model reduction technique is used in conjunction */
/*         with the method proposed in [1]. */

/*     If ORDSEL = 'A', the order of the controller is determined by */
/*     computing the number of Hankel singular values greater than */
/*     the given tolerance TOL. The Hankel singular values are */
/*     the square roots of the eigenvalues of the product of */
/*     two frequency-weighted Grammians P and Q, defined as follows. */

/*     If JOBCF = 'L', then P is the controllability Grammian of a system */
/*     of the form (A+BF,B,*,*), and Q is the observability Grammian of a */
/*     system of the form (A+GC,*,F,*). This choice corresponds to an */
/*     input frequency-weighted order reduction of left coprime */
/*     factors [1]. */

/*     If JOBCF = 'R', then P is the controllability Grammian of a system */
/*     of the form (A+BF,G,*,*), and Q is the observability Grammian of a */
/*     system of the form (A+GC,*,C,*). This choice corresponds to an */
/*     output frequency-weighted order reduction of right coprime */
/*     factors [1]. */

/*     For the computation of truncation matrices, the B&T approach */
/*     is used in conjunction with accuracy enhancing techniques. */
/*     If JOBMR = 'B', the square-root B&T method of [2,4] is used. */
/*     If JOBMR = 'F', the balancing-free square-root version of the */
/*     B&T method [3,4] is used. */

/*     REFERENCES */

/*     [1] Liu, Y., Anderson, B.D.O. and Ly, O.L. */
/*         Coprime factorization controller reduction with Bezout */
/*         identity induced frequency weighting. */
/*         Automatica, vol. 26, pp. 233-249, 1990. */

/*     [2] Tombs, M.S. and Postlethwaite I. */
/*         Truncated balanced realization of stable, non-minimal */
/*         state-space systems. */
/*         Int. J. Control, Vol. 46, pp. 1319-1330, 1987. */

/*     [3] Varga, A. */
/*         Efficient minimal realization procedure based on balancing. */
/*         Proc. of IMACS/IFAC Symp. MCTS, Lille, France, May 1991, */
/*         A. El Moudui, P. Borne, S. G. Tzafestas (Eds.), Vol. 2, */
/*         pp. 42-46, 1991. */

/*     [4] Varga, A. */
/*         Coprime factors model reduction method based on square-root */
/*         balancing-free techniques. */
/*         System Analysis, Modelling and Simulation, Vol. 11, */
/*         pp. 303-311, 1993. */

/*     NUMERICAL ASPECTS */

/*     The implemented methods rely on accuracy enhancing square-root or */
/*     balancing-free square-root techniques. */
/*                                         3 */
/*     The algorithms require less than 30N  floating point operations. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, October 2000. */
/*     D. Sima, University of Bucharest, October 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2000. */

/*     REVISIONS */

/*     A. Varga, Australian National University, Canberra, November 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2001. */

/*     KEYWORDS */

/*     Controller reduction, coprime factorization, frequency weighting, */
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

#line 388 "SB16CD.f"
    /* Parameter adjustments */
#line 388 "SB16CD.f"
    a_dim1 = *lda;
#line 388 "SB16CD.f"
    a_offset = 1 + a_dim1;
#line 388 "SB16CD.f"
    a -= a_offset;
#line 388 "SB16CD.f"
    b_dim1 = *ldb;
#line 388 "SB16CD.f"
    b_offset = 1 + b_dim1;
#line 388 "SB16CD.f"
    b -= b_offset;
#line 388 "SB16CD.f"
    c_dim1 = *ldc;
#line 388 "SB16CD.f"
    c_offset = 1 + c_dim1;
#line 388 "SB16CD.f"
    c__ -= c_offset;
#line 388 "SB16CD.f"
    d_dim1 = *ldd;
#line 388 "SB16CD.f"
    d_offset = 1 + d_dim1;
#line 388 "SB16CD.f"
    d__ -= d_offset;
#line 388 "SB16CD.f"
    f_dim1 = *ldf;
#line 388 "SB16CD.f"
    f_offset = 1 + f_dim1;
#line 388 "SB16CD.f"
    f -= f_offset;
#line 388 "SB16CD.f"
    g_dim1 = *ldg;
#line 388 "SB16CD.f"
    g_offset = 1 + g_dim1;
#line 388 "SB16CD.f"
    g -= g_offset;
#line 388 "SB16CD.f"
    --hsv;
#line 388 "SB16CD.f"
    --iwork;
#line 388 "SB16CD.f"
    --dwork;
#line 388 "SB16CD.f"

#line 388 "SB16CD.f"
    /* Function Body */
#line 388 "SB16CD.f"
    *info = 0;
#line 389 "SB16CD.f"
    *iwarn = 0;
#line 390 "SB16CD.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 391 "SB16CD.f"
    withd = lsame_(jobd, "D", (ftnlen)1, (ftnlen)1);
#line 392 "SB16CD.f"
    bal = lsame_(jobmr, "B", (ftnlen)1, (ftnlen)1);
#line 393 "SB16CD.f"
    left = lsame_(jobcf, "L", (ftnlen)1, (ftnlen)1);
#line 394 "SB16CD.f"
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);
#line 395 "SB16CD.f"
    if (left) {
#line 396 "SB16CD.f"
	mp = *m;
#line 397 "SB16CD.f"
    } else {
#line 398 "SB16CD.f"
	mp = *p;
#line 399 "SB16CD.f"
    }
/* Computing MAX */
#line 400 "SB16CD.f"
    i__1 = 1, i__2 = (*n << 1) * *n + *n * 5, i__1 = max(i__1,i__2), i__2 = *
	    n * max(*m,*p), i__1 = max(i__1,i__2), i__2 = *n * (*n + max(*n,
	    mp) + min(*n,mp) + 6);
#line 400 "SB16CD.f"
    lw = (*n << 1) * *n + max(i__1,i__2);

/*     Test the input scalar arguments. */

#line 405 "SB16CD.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 406 "SB16CD.f"
	*info = -1;
#line 407 "SB16CD.f"
    } else if (! (withd || lsame_(jobd, "Z", (ftnlen)1, (ftnlen)1))) {
#line 408 "SB16CD.f"
	*info = -2;
#line 409 "SB16CD.f"
    } else if (! (bal || lsame_(jobmr, "F", (ftnlen)1, (ftnlen)1))) {
#line 410 "SB16CD.f"
	*info = -3;
#line 411 "SB16CD.f"
    } else if (! (left || lsame_(jobcf, "R", (ftnlen)1, (ftnlen)1))) {
#line 412 "SB16CD.f"
	*info = -4;
#line 413 "SB16CD.f"
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
#line 414 "SB16CD.f"
	*info = -5;
#line 415 "SB16CD.f"
    } else if (*n < 0) {
#line 416 "SB16CD.f"
	*info = -6;
#line 417 "SB16CD.f"
    } else if (*m < 0) {
#line 418 "SB16CD.f"
	*info = -7;
#line 419 "SB16CD.f"
    } else if (*p < 0) {
#line 420 "SB16CD.f"
	*info = -8;
#line 421 "SB16CD.f"
    } else if (fixord && (*ncr < 0 || *ncr > *n)) {
#line 422 "SB16CD.f"
	*info = -9;
#line 423 "SB16CD.f"
    } else if (*lda < max(1,*n)) {
#line 424 "SB16CD.f"
	*info = -11;
#line 425 "SB16CD.f"
    } else if (*ldb < max(1,*n)) {
#line 426 "SB16CD.f"
	*info = -13;
#line 427 "SB16CD.f"
    } else if (*ldc < max(1,*p)) {
#line 428 "SB16CD.f"
	*info = -15;
#line 429 "SB16CD.f"
    } else if (*ldd < 1 || withd && *ldd < *p) {
#line 430 "SB16CD.f"
	*info = -17;
#line 431 "SB16CD.f"
    } else if (*ldf < max(1,*m)) {
#line 432 "SB16CD.f"
	*info = -19;
#line 433 "SB16CD.f"
    } else if (*ldg < max(1,*n)) {
#line 434 "SB16CD.f"
	*info = -21;
#line 435 "SB16CD.f"
    } else if (*ldwork < lw) {
#line 436 "SB16CD.f"
	*info = -26;
#line 437 "SB16CD.f"
    }

#line 439 "SB16CD.f"
    if (*info != 0) {

/*        Error return. */

#line 443 "SB16CD.f"
	i__1 = -(*info);
#line 443 "SB16CD.f"
	xerbla_("SB16CD", &i__1, (ftnlen)6);
#line 444 "SB16CD.f"
	return 0;
#line 445 "SB16CD.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 449 "SB16CD.f"
    i__1 = min(*n,*m);
#line 449 "SB16CD.f"
    if (min(i__1,*p) == 0 || fixord && *ncr == 0) {
#line 451 "SB16CD.f"
	*ncr = 0;
#line 452 "SB16CD.f"
	dwork[1] = 1.;
#line 453 "SB16CD.f"
	return 0;
#line 454 "SB16CD.f"
    }

/*     Allocate working storage. */

#line 458 "SB16CD.f"
    kt = 1;
#line 459 "SB16CD.f"
    kti = kt + *n * *n;
#line 460 "SB16CD.f"
    kw = kti + *n * *n;

/*     Compute in DWORK(KTI) and DWORK(KT) the Cholesky factors Su and Ru */
/*     of the frequency-weighted controllability and observability */
/*     Grammians, respectively. */

/*     Workspace:   need 2*N*N + MAX(1, N*(N + MAX(N,M) + MIN(N,M) + 6)), */
/*                                                        if JOBCF = 'L'; */
/*                       2*N*N + MAX(1, N*(N + MAX(N,P) + MIN(N,P) + 6)), */
/*                                                        if JOBCF = 'R'. */
/*                  prefer larger. */

#line 472 "SB16CD.f"
    i__1 = *ldwork - kw + 1;
#line 472 "SB16CD.f"
    sb16cy_(dico, jobcf, n, m, p, &a[a_offset], lda, &b[b_offset], ldb, &c__[
	    c_offset], ldc, &f[f_offset], ldf, &g[g_offset], ldg, &scalec, &
	    scaleo, &dwork[kti], n, &dwork[kt], n, &dwork[kw], &i__1, info, (
	    ftnlen)1, (ftnlen)1);

#line 476 "SB16CD.f"
    if (*info != 0) {
#line 476 "SB16CD.f"
	return 0;
#line 476 "SB16CD.f"
    }
#line 478 "SB16CD.f"
    wrkopt = (integer) dwork[kw] + kw - 1;

/*     Compute a B&T approximation (Ar,Br,Cr) of (A,B,C) and */
/*     the corresponding truncation matrices TI and T. */

/*     Real workspace:  need   2*N*N + MAX( 1, 2*N*N+5*N, N*MAX(M,P) ); */
/*                      prefer larger. */
/*     Integer workspace:  0,  if JOBMR = 'B'; */
/*                         N,  if JOBMR = 'F'. */

#line 488 "SB16CD.f"
    i__1 = *ldwork - kw + 1;
#line 488 "SB16CD.f"
    ab09ix_(dico, jobmr, "NotSchur", ordsel, n, m, p, ncr, &scalec, &scaleo, &
	    a[a_offset], lda, &b[b_offset], ldb, &c__[c_offset], ldc, &d__[
	    d_offset], ldd, &dwork[kti], n, &dwork[kt], n, &nmr, &hsv[1], tol,
	     tol, &iwork[1], &dwork[kw], &i__1, iwarn, &ierr, (ftnlen)1, (
	    ftnlen)1, (ftnlen)8, (ftnlen)1);
#line 492 "SB16CD.f"
    if (ierr != 0) {
#line 493 "SB16CD.f"
	*info = 6;
#line 494 "SB16CD.f"
	return 0;
#line 495 "SB16CD.f"
    }
/* Computing MAX */
#line 496 "SB16CD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 496 "SB16CD.f"
    wrkopt = max(i__1,i__2);

/*     Compute reduced gains Bc = Gr = TI*G and Cc = Fr = F*T. */
/*     Workspace:  need   N*(2*N+MAX(M,P)). */

#line 501 "SB16CD.f"
    dlacpy_("Full", n, p, &g[g_offset], ldg, &dwork[kw], n, (ftnlen)4);
#line 502 "SB16CD.f"
    dgemm_("NoTranspose", "NoTranspose", ncr, p, n, &c_b17, &dwork[kti], n, &
	    dwork[kw], n, &c_b18, &g[g_offset], ldg, (ftnlen)11, (ftnlen)11);

#line 505 "SB16CD.f"
    dlacpy_("Full", m, n, &f[f_offset], ldf, &dwork[kw], m, (ftnlen)4);
#line 506 "SB16CD.f"
    dgemm_("NoTranspose", "NoTranspose", m, ncr, n, &c_b17, &dwork[kw], m, &
	    dwork[kt], n, &c_b18, &f[f_offset], ldf, (ftnlen)11, (ftnlen)11);

/*     Form the reduced controller state matrix, */
/*     Ac = Ar + Br*Fr + Gr*Cr + Gr*D*Fr = Ar + Br*Fr + Gr*(Cr+D*Fr) . */

/*     Workspace:    need  P*N. */

#line 514 "SB16CD.f"
    dlacpy_("Full", p, ncr, &c__[c_offset], ldc, &dwork[1], p, (ftnlen)4);
#line 515 "SB16CD.f"
    if (withd) {
#line 515 "SB16CD.f"
	dgemm_("NoTranspose", "NoTranspose", p, ncr, m, &c_b17, &d__[d_offset]
		, ldd, &f[f_offset], ldf, &c_b17, &dwork[1], p, (ftnlen)11, (
		ftnlen)11);
#line 515 "SB16CD.f"
    }
#line 517 "SB16CD.f"
    dgemm_("NoTranspose", "NoTranspose", ncr, ncr, p, &c_b17, &g[g_offset], 
	    ldg, &dwork[1], p, &c_b17, &a[a_offset], lda, (ftnlen)11, (ftnlen)
	    11);
#line 519 "SB16CD.f"
    dgemm_("NoTranspose", "NoTranspose", ncr, ncr, m, &c_b17, &b[b_offset], 
	    ldb, &f[f_offset], ldf, &c_b17, &a[a_offset], lda, (ftnlen)11, (
	    ftnlen)11);

#line 522 "SB16CD.f"
    dwork[1] = (doublereal) wrkopt;

#line 524 "SB16CD.f"
    return 0;
/* *** Last line of SB16CD *** */
} /* sb16cd_ */

