#line 1 "SB16BD.f"
/* SB16BD.f -- translated by f2c (version 20100827).
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

#line 1 "SB16BD.f"
/* Table of constant values */

static doublereal c_b22 = 1.;
static doublereal c_b45 = 0.;

/* Subroutine */ int sb16bd_(char *dico, char *jobd, char *jobmr, char *jobcf,
	 char *equil, char *ordsel, integer *n, integer *m, integer *p, 
	integer *ncr, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *f, integer *ldf, doublereal *g, integer *ldg, doublereal *
	dc, integer *lddc, doublereal *hsv, doublereal *tol1, doublereal *
	tol2, integer *iwork, doublereal *dwork, integer *ldwork, integer *
	iwarn, integer *info, ftnlen dico_len, ftnlen jobd_len, ftnlen 
	jobmr_len, ftnlen jobcf_len, ftnlen equil_len, ftnlen ordsel_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, dc_dim1, dc_offset, f_dim1, f_offset, g_dim1, g_offset, 
	    i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer kw, lw1, lw2;
    static logical bal;
    static integer kbe;
    static logical bta;
    static integer kce, kde;
    static char job[1];
    static logical spa;
    static integer lwr, ldbe, ldce, ldde;
    static logical left;
    extern /* Subroutine */ int ab09ad_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen), ab09bd_(char *, char *
	    , char *, char *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen), sb08gd_(integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *), sb08hd_(integer *, integer *, integer *,
	     doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *), 
	    dgemm_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr, withd;
    static integer maxmp;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static logical fixord, lequil;
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
/*     controller model (Ac,Bc,Cc,Dc) using a coprime factorization */
/*     based controller reduction approach. For reduction, */
/*     either the square-root or the balancing-free square-root */
/*     versions of the Balance & Truncate (B&T) or Singular Perturbation */
/*     Approximation (SPA) model reduction methods are used in */
/*     conjunction with stable coprime factorization techniques. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the open-loop system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     JOBD    CHARACTER*1 */
/*             Specifies whether or not a non-zero matrix D appears */
/*             in the given state space model: */
/*             = 'D':  D is present; */
/*             = 'Z':  D is assumed a zero matrix. */

/*     JOBMR   CHARACTER*1 */
/*             Specifies the model reduction approach to be used */
/*             as follows: */
/*             = 'B':  use the square-root B&T method; */
/*             = 'F':  use the balancing-free square-root B&T method; */
/*             = 'S':  use the square-root SPA method; */
/*             = 'P':  use the balancing-free square-root SPA method. */

/*     JOBCF   CHARACTER*1 */
/*             Specifies whether left or right coprime factorization is */
/*             to be used as follows: */
/*             = 'L':  use left coprime factorization; */
/*             = 'R':  use right coprime factorization. */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to perform a */
/*             preliminary equilibration before performing */
/*             order reduction as follows: */
/*             = 'S':  perform equilibration (scaling); */
/*             = 'N':  do not perform equilibration. */

/*     ORDSEL  CHARACTER*1 */
/*             Specifies the order selection method as follows: */
/*             = 'F':  the resulting controller order NCR is fixed; */
/*             = 'A':  the resulting controller order NCR is */
/*                     automatically determined on basis of the given */
/*                     tolerance TOL1. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the open-loop state-space representation, */
/*             i.e., the order of the matrix A.  N >= 0. */
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
/*             if ORDSEL = 'F', NCR is equal to MIN(NCR,NMIN), where NCR */
/*             is the desired order on entry, and NMIN is the order of a */
/*             minimal realization of an extended system Ge (see METHOD); */
/*             NMIN is determined as the number of */
/*             Hankel singular values greater than N*EPS*HNORM(Ge), */
/*             where EPS is the machine precision (see LAPACK Library */
/*             Routine DLAMCH) and HNORM(Ge) is the Hankel norm of the */
/*             extended system (computed in HSV(1)); */
/*             if ORDSEL = 'A', NCR is equal to the number of Hankel */
/*             singular values greater than MAX(TOL1,N*EPS*HNORM(Ge)). */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original state dynamics matrix A. */
/*             On exit, if INFO = 0, the leading NCR-by-NCR part of this */
/*             array contains the state dynamics matrix Ac of the reduced */
/*             controller. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must */
/*             contain the original input/state matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array must */
/*             contain the original state/output matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             If JOBD = 'D', the leading P-by-M part of this */
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
/*             array contains the state/output matrix Cc of the reduced */
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

/*     DC      (output) DOUBLE PRECISION array, dimension (LDDC,P) */
/*             If INFO = 0, the leading M-by-P part of this array */
/*             contains the input/output matrix Dc of the reduced */
/*             controller. */

/*     LDDC    INTEGER */
/*             The leading dimension of array DC.  LDDC >= MAX(1,M). */

/*     HSV     (output) DOUBLE PRECISION array, dimension (N) */
/*             If INFO = 0, it contains the N Hankel singular values */
/*             of the extended system ordered decreasingly (see METHOD). */

/*     Tolerances */

/*     TOL1    DOUBLE PRECISION */
/*             If ORDSEL = 'A', TOL1 contains the tolerance for */
/*             determining the order of the reduced extended system. */
/*             For model reduction, the recommended value is */
/*             TOL1 = c*HNORM(Ge), where c is a constant in the */
/*             interval [0.00001,0.001], and HNORM(Ge) is the */
/*             Hankel norm of the extended system (computed in HSV(1)). */
/*             The value TOL1 = N*EPS*HNORM(Ge) is used by default if */
/*             TOL1 <= 0 on entry, where EPS is the machine precision */
/*             (see LAPACK Library Routine DLAMCH). */
/*             If ORDSEL = 'F', the value of TOL1 is ignored. */

/*     TOL2    DOUBLE PRECISION */
/*             The tolerance for determining the order of a minimal */
/*             realization of the coprime factorization controller */
/*             (see METHOD). The recommended value is */
/*             TOL2 = N*EPS*HNORM(Ge) (see METHOD). */
/*             This value is used by default if TOL2 <= 0 on entry. */
/*             If TOL2 > 0 and ORDSEL = 'A', then TOL2 <= TOL1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK = 0,         if ORDSEL = 'F' and NCR = N. */
/*                                                 Otherwise, */
/*             LIWORK = MAX(PM,M), if JOBCF = 'L', */
/*             LIWORK = MAX(PM,P), if JOBCF = 'R', where */
/*             PM = 0,             if JOBMR = 'B', */
/*             PM = N,             if JOBMR = 'F', */
/*             PM = MAX(1,2*N),    if JOBMR = 'S' or 'P'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= P*N, if ORDSEL = 'F' and NCR = N. Otherwise, */
/*             LDWORK >= (N+M)*(M+P) + MAX(LWR,4*M), if JOBCF = 'L', */
/*             LDWORK >= (N+P)*(M+P) + MAX(LWR,4*P), if JOBCF = 'R', */
/*             where LWR = MAX(1,N*(2*N+MAX(N,M+P)+5)+N*(N+1)/2). */
/*             For optimum performance LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  with ORDSEL = 'F', the selected order NCR is */
/*                   greater than the order of a minimal */
/*                   realization of the controller. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the reduction of A+G*C to a real Schur form */
/*                   failed; */
/*             = 2:  the matrix A+G*C is not stable (if DICO = 'C'), */
/*                   or not convergent (if DICO = 'D'); */
/*             = 3:  the computation of Hankel singular values failed; */
/*             = 4:  the reduction of A+B*F to a real Schur form */
/*                   failed; */
/*             = 5:  the matrix A+B*F is not stable (if DICO = 'C'), */
/*                   or not convergent (if DICO = 'D'). */

/*     METHOD */

/*     Let be the linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t) + Du(t),                             (1) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system, and let Go(d) be the open-loop */
/*     transfer-function matrix */
/*                           -1 */
/*          Go(d) = C*(d*I-A) *B + D . */

/*     Let F and G be the state feedback and observer gain matrices, */
/*     respectively, chosen so that A+B*F and A+G*C are stable matrices. */
/*     The controller has a transfer-function matrix K(d) given by */
/*                                        -1 */
/*          K(d) = F*(d*I-A-B*F-G*C-G*D*F) *G . */

/*     The closed-loop transfer-function matrix is given by */
/*                                     -1 */
/*          Gcl(d) = Go(d)(I+K(d)Go(d)) . */

/*     K(d) can be expressed as a left coprime factorization (LCF), */
/*                          -1 */
/*          K(d) = M_left(d) *N_left(d) , */

/*     or as a right coprime factorization (RCF), */
/*                                      -1 */
/*          K(d) = N_right(d)*M_right(d) , */

/*     where M_left(d), N_left(d), N_right(d), and M_right(d) are */
/*     stable transfer-function matrices. */

/*     The subroutine SB16BD determines the matrices of a reduced */
/*     controller */

/*          d[z(t)] = Ac*z(t) + Bc*y(t) */
/*          u(t)    = Cc*z(t) + Dc*y(t),                           (2) */

/*     with the transfer-function matrix Kr as follows: */

/*     (1) If JOBCF = 'L', the extended system */
/*         Ge(d)  = [ N_left(d) M_left(d) ] is reduced to */
/*         Ger(d) = [ N_leftr(d) M_leftr(d) ] by using either the */
/*         B&T or SPA methods. The reduced order controller Kr(d) */
/*         is computed as */
/*                           -1 */
/*         Kr(d) = M_leftr(d) *N_leftr(d) ; */

/*     (2) If JOBCF = 'R', the extended system */
/*         Ge(d) = [ N_right(d) ] is reduced to */
/*                 [ M_right(d) ] */
/*         Ger(d) = [ N_rightr(d) ] by using either the */
/*                  [ M_rightr(d) ] */
/*         B&T or SPA methods. The reduced order controller Kr(d) */
/*         is computed as */
/*                                         -1 */
/*         Kr(d) = N_rightr(d)* M_rightr(d) . */

/*     If ORDSEL = 'A', the order of the controller is determined by */
/*     computing the number of Hankel singular values greater than */
/*     the given tolerance TOL1. The Hankel singular values are */
/*     the square roots of the eigenvalues of the product of */
/*     the controllability and observability Grammians of the */
/*     extended system Ge. */

/*     If JOBMR = 'B', the square-root B&T method of [1] is used. */

/*     If JOBMR = 'F', the balancing-free square-root version of the */
/*     B&T method [1] is used. */

/*     If JOBMR = 'S', the square-root version of the SPA method [2,3] */
/*     is used. */

/*     If JOBMR = 'P', the balancing-free square-root version of the */
/*     SPA method [2,3] is used. */

/*     REFERENCES */

/*     [1] Tombs, M.S. and Postlethwaite, I. */
/*         Truncated balanced realization of stable, non-minimal */
/*         state-space systems. */
/*         Int. J. Control, Vol. 46, pp. 1319-1330, 1987. */

/*     [2] Varga, A. */
/*         Efficient minimal realization procedure based on balancing. */
/*         Proc. of IMACS/IFAC Symp. MCTS, Lille, France, May 1991, */
/*         A. El Moudui, P. Borne, S. G. Tzafestas (Eds.), Vol. 2, */
/*         pp. 42-46, 1991. */

/*     [3] Varga, A. */
/*         Coprime factors model reduction method based on square-root */
/*         balancing-free techniques. */
/*         System Analysis, Modelling and Simulation, Vol. 11, */
/*         pp. 303-311, 1993. */

/*     [4] Liu, Y., Anderson, B.D.O. and Ly, O.L. */
/*         Coprime factorization controller reduction with Bezout */
/*         identity induced frequency weighting. */
/*         Automatica, vol. 26, pp. 233-249, 1990. */

/*     NUMERICAL ASPECTS */

/*     The implemented methods rely on accuracy enhancing square-root or */
/*     balancing-free square-root techniques. */
/*                                         3 */
/*     The algorithms require less than 30N  floating point operations. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, August 2000. */
/*     D. Sima, University of Bucharest, August 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2000. */

/*     REVISIONS */

/*     A. Varga, Australian National University, Canberra, November 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000, */
/*              Aug. 2001. */

/*     KEYWORDS */

/*     Balancing, controller reduction, coprime factorization, */
/*     minimal realization, multivariable system, state-space model. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 401 "SB16BD.f"
    /* Parameter adjustments */
#line 401 "SB16BD.f"
    a_dim1 = *lda;
#line 401 "SB16BD.f"
    a_offset = 1 + a_dim1;
#line 401 "SB16BD.f"
    a -= a_offset;
#line 401 "SB16BD.f"
    b_dim1 = *ldb;
#line 401 "SB16BD.f"
    b_offset = 1 + b_dim1;
#line 401 "SB16BD.f"
    b -= b_offset;
#line 401 "SB16BD.f"
    c_dim1 = *ldc;
#line 401 "SB16BD.f"
    c_offset = 1 + c_dim1;
#line 401 "SB16BD.f"
    c__ -= c_offset;
#line 401 "SB16BD.f"
    d_dim1 = *ldd;
#line 401 "SB16BD.f"
    d_offset = 1 + d_dim1;
#line 401 "SB16BD.f"
    d__ -= d_offset;
#line 401 "SB16BD.f"
    f_dim1 = *ldf;
#line 401 "SB16BD.f"
    f_offset = 1 + f_dim1;
#line 401 "SB16BD.f"
    f -= f_offset;
#line 401 "SB16BD.f"
    g_dim1 = *ldg;
#line 401 "SB16BD.f"
    g_offset = 1 + g_dim1;
#line 401 "SB16BD.f"
    g -= g_offset;
#line 401 "SB16BD.f"
    dc_dim1 = *lddc;
#line 401 "SB16BD.f"
    dc_offset = 1 + dc_dim1;
#line 401 "SB16BD.f"
    dc -= dc_offset;
#line 401 "SB16BD.f"
    --hsv;
#line 401 "SB16BD.f"
    --iwork;
#line 401 "SB16BD.f"
    --dwork;
#line 401 "SB16BD.f"

#line 401 "SB16BD.f"
    /* Function Body */
#line 401 "SB16BD.f"
    *info = 0;
#line 402 "SB16BD.f"
    *iwarn = 0;
#line 403 "SB16BD.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 404 "SB16BD.f"
    withd = lsame_(jobd, "D", (ftnlen)1, (ftnlen)1);
#line 405 "SB16BD.f"
    bta = lsame_(jobmr, "B", (ftnlen)1, (ftnlen)1) || lsame_(jobmr, "F", (
	    ftnlen)1, (ftnlen)1);
#line 406 "SB16BD.f"
    spa = lsame_(jobmr, "S", (ftnlen)1, (ftnlen)1) || lsame_(jobmr, "P", (
	    ftnlen)1, (ftnlen)1);
#line 407 "SB16BD.f"
    bal = lsame_(jobmr, "B", (ftnlen)1, (ftnlen)1) || lsame_(jobmr, "S", (
	    ftnlen)1, (ftnlen)1);
#line 408 "SB16BD.f"
    left = lsame_(jobcf, "L", (ftnlen)1, (ftnlen)1);
#line 409 "SB16BD.f"
    lequil = lsame_(equil, "S", (ftnlen)1, (ftnlen)1);
#line 410 "SB16BD.f"
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);
#line 411 "SB16BD.f"
    maxmp = max(*m,*p);

/* Computing MAX */
/* Computing MAX */
#line 413 "SB16BD.f"
    i__3 = *n, i__4 = *m + *p;
#line 413 "SB16BD.f"
    i__1 = 1, i__2 = *n * ((*n << 1) + max(i__3,i__4) + 5) + *n * (*n + 1) / 
	    2;
#line 413 "SB16BD.f"
    lwr = max(i__1,i__2);
/* Computing MAX */
#line 414 "SB16BD.f"
    i__1 = lwr, i__2 = *m << 2;
#line 414 "SB16BD.f"
    lw1 = (*n + *m) * (*m + *p) + max(i__1,i__2);
/* Computing MAX */
#line 415 "SB16BD.f"
    i__1 = lwr, i__2 = *p << 2;
#line 415 "SB16BD.f"
    lw2 = (*n + *p) * (*m + *p) + max(i__1,i__2);

/*     Test the input scalar arguments. */

#line 419 "SB16BD.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 420 "SB16BD.f"
	*info = -1;
#line 421 "SB16BD.f"
    } else if (! (withd || lsame_(jobd, "Z", (ftnlen)1, (ftnlen)1))) {
#line 422 "SB16BD.f"
	*info = -2;
#line 423 "SB16BD.f"
    } else if (! (bta || spa)) {
#line 424 "SB16BD.f"
	*info = -3;
#line 425 "SB16BD.f"
    } else if (! (left || lsame_(jobcf, "R", (ftnlen)1, (ftnlen)1))) {
#line 426 "SB16BD.f"
	*info = -4;
#line 427 "SB16BD.f"
    } else if (! (lequil || lsame_(equil, "N", (ftnlen)1, (ftnlen)1))) {
#line 428 "SB16BD.f"
	*info = -5;
#line 429 "SB16BD.f"
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
#line 430 "SB16BD.f"
	*info = -6;
#line 431 "SB16BD.f"
    } else if (*n < 0) {
#line 432 "SB16BD.f"
	*info = -7;
#line 433 "SB16BD.f"
    } else if (*m < 0) {
#line 434 "SB16BD.f"
	*info = -8;
#line 435 "SB16BD.f"
    } else if (*p < 0) {
#line 436 "SB16BD.f"
	*info = -9;
#line 437 "SB16BD.f"
    } else if (fixord && (*ncr < 0 || *ncr > *n)) {
#line 438 "SB16BD.f"
	*info = -10;
#line 439 "SB16BD.f"
    } else if (*lda < max(1,*n)) {
#line 440 "SB16BD.f"
	*info = -12;
#line 441 "SB16BD.f"
    } else if (*ldb < max(1,*n)) {
#line 442 "SB16BD.f"
	*info = -14;
#line 443 "SB16BD.f"
    } else if (*ldc < max(1,*p)) {
#line 444 "SB16BD.f"
	*info = -16;
#line 445 "SB16BD.f"
    } else if (*ldd < 1 || withd && *ldd < *p) {
#line 446 "SB16BD.f"
	*info = -18;
#line 447 "SB16BD.f"
    } else if (*ldf < max(1,*m)) {
#line 448 "SB16BD.f"
	*info = -20;
#line 449 "SB16BD.f"
    } else if (*ldg < max(1,*n)) {
#line 450 "SB16BD.f"
	*info = -22;
#line 451 "SB16BD.f"
    } else if (*lddc < max(1,*m)) {
#line 452 "SB16BD.f"
	*info = -24;
#line 453 "SB16BD.f"
    } else if (! fixord && *tol2 > 0. && *tol2 > *tol1) {
#line 454 "SB16BD.f"
	*info = -27;
#line 455 "SB16BD.f"
    } else if ((! fixord || *ncr < *n) && (left && *ldwork < lw1) || ! left &&
	     *ldwork < lw2 || fixord && *ncr == *n && *ldwork < *p * *n) {
#line 459 "SB16BD.f"
	*info = -30;
#line 460 "SB16BD.f"
    }

#line 462 "SB16BD.f"
    if (*info != 0) {

/*        Error return. */

#line 466 "SB16BD.f"
	i__1 = -(*info);
#line 466 "SB16BD.f"
	xerbla_("SB16BD", &i__1, (ftnlen)6);
#line 467 "SB16BD.f"
	return 0;
#line 468 "SB16BD.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 472 "SB16BD.f"
    i__1 = min(*n,*m);
#line 472 "SB16BD.f"
    if (min(i__1,*p) == 0 || fixord && bta && *ncr == 0) {
#line 474 "SB16BD.f"
	*ncr = 0;
#line 475 "SB16BD.f"
	dwork[1] = 1.;
#line 476 "SB16BD.f"
	return 0;
#line 477 "SB16BD.f"
    }

#line 479 "SB16BD.f"
    if (*ncr == *n) {

/*        Form the controller state matrix, */
/*        Ac = A + B*F + G*C + G*D*F = A + B*F + G*(C+D*F) . */
/*        Real workspace:    need  P*N. */
/*        Integer workspace: need  0. */

#line 486 "SB16BD.f"
	dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[1], p, (ftnlen)4);
#line 487 "SB16BD.f"
	if (withd) {
#line 487 "SB16BD.f"
	    dgemm_("NoTranspose", "NoTranspose", p, n, m, &c_b22, &d__[
		    d_offset], ldd, &f[f_offset], ldf, &c_b22, &dwork[1], p, (
		    ftnlen)11, (ftnlen)11);
#line 487 "SB16BD.f"
	}
#line 490 "SB16BD.f"
	dgemm_("NoTranspose", "NoTranspose", n, n, p, &c_b22, &g[g_offset], 
		ldg, &dwork[1], p, &c_b22, &a[a_offset], lda, (ftnlen)11, (
		ftnlen)11);
#line 492 "SB16BD.f"
	dgemm_("NoTranspose", "NoTranspose", n, n, m, &c_b22, &b[b_offset], 
		ldb, &f[f_offset], ldf, &c_b22, &a[a_offset], lda, (ftnlen)11,
		 (ftnlen)11);

#line 495 "SB16BD.f"
	dwork[1] = (doublereal) (*p * *n);
#line 496 "SB16BD.f"
	return 0;
#line 497 "SB16BD.f"
    }

#line 499 "SB16BD.f"
    if (bal) {
#line 500 "SB16BD.f"
	*(unsigned char *)job = 'B';
#line 501 "SB16BD.f"
    } else {
#line 502 "SB16BD.f"
	*(unsigned char *)job = 'N';
#line 503 "SB16BD.f"
    }

/*     Reduce the coprime factors. */

#line 507 "SB16BD.f"
    if (left) {

/*        Form Ge(d) = [ N_left(d) M_left(d) ] as */

/*             ( A+G*C |  G  B+GD ) */
/*             (------------------) */
/*             (   F   |  0   I   ) */

/*        Real workspace:    need  (N+M)*(M+P). */
/*        Integer workspace: need  0. */

#line 518 "SB16BD.f"
	dgemm_("NoTranspose", "NoTranspose", n, n, p, &c_b22, &g[g_offset], 
		ldg, &c__[c_offset], ldc, &c_b22, &a[a_offset], lda, (ftnlen)
		11, (ftnlen)11);
#line 520 "SB16BD.f"
	kbe = 1;
#line 521 "SB16BD.f"
	kde = kbe + *n * (*p + *m);
#line 522 "SB16BD.f"
	ldbe = max(1,*n);
#line 523 "SB16BD.f"
	ldde = *m;
#line 524 "SB16BD.f"
	dlacpy_("Full", n, p, &g[g_offset], ldg, &dwork[kbe], &ldbe, (ftnlen)
		4);
#line 525 "SB16BD.f"
	dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[kbe + *n * *p], &ldbe,
		 (ftnlen)4);
#line 526 "SB16BD.f"
	if (withd) {
#line 526 "SB16BD.f"
	    dgemm_("NoTranspose", "NoTranspose", n, m, p, &c_b22, &g[g_offset]
		    , ldg, &d__[d_offset], ldd, &c_b22, &dwork[kbe + *n * *p],
		     &ldbe, (ftnlen)11, (ftnlen)11);
#line 526 "SB16BD.f"
	}
#line 529 "SB16BD.f"
	dlaset_("Full", m, p, &c_b45, &c_b45, &dwork[kde], &ldde, (ftnlen)4);
#line 530 "SB16BD.f"
	dlaset_("Full", m, m, &c_b45, &c_b22, &dwork[kde + *m * *p], &ldde, (
		ftnlen)4);

/*        Compute the reduced coprime factors, */
/*             Ger(d) = [ N_leftr(d) M_leftr(d) ] , */
/*        by using either the B&T or SPA methods. */

/*        Real workspace:    need  (N+M)*(M+P) + */
/*                                 MAX(1,N*(2*N+MAX(N,M+P)+5)+N*(N+1)/2). */
/*        Integer workspace: need  0,         if JOBMR = 'B', */
/*                                 N,         if JOBMR = 'F', and */
/*                                 MAX(1,2*N) if JOBMR = 'S' or 'P'. */

#line 542 "SB16BD.f"
	kw = kde + *m * (*p + *m);
#line 543 "SB16BD.f"
	if (bta) {
#line 544 "SB16BD.f"
	    i__1 = *m + *p;
#line 544 "SB16BD.f"
	    i__2 = *ldwork - kw + 1;
#line 544 "SB16BD.f"
	    ab09ad_(dico, job, equil, ordsel, n, &i__1, m, ncr, &a[a_offset], 
		    lda, &dwork[kbe], &ldbe, &f[f_offset], ldf, &hsv[1], tol1,
		     &iwork[1], &dwork[kw], &i__2, iwarn, info, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 547 "SB16BD.f"
	} else {
#line 548 "SB16BD.f"
	    i__1 = *m + *p;
#line 548 "SB16BD.f"
	    i__2 = *ldwork - kw + 1;
#line 548 "SB16BD.f"
	    ab09bd_(dico, job, equil, ordsel, n, &i__1, m, ncr, &a[a_offset], 
		    lda, &dwork[kbe], &ldbe, &f[f_offset], ldf, &dwork[kde], &
		    ldde, &hsv[1], tol1, tol2, &iwork[1], &dwork[kw], &i__2, 
		    iwarn, info, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 552 "SB16BD.f"
	}
#line 553 "SB16BD.f"
	if (*info != 0) {
#line 553 "SB16BD.f"
	    return 0;
#line 553 "SB16BD.f"
	}

#line 556 "SB16BD.f"
	wrkopt = (integer) dwork[kw] + kw - 1;

/*        Compute the reduced order controller, */
/*                             -1 */
/*           Kr(d) = M_leftr(d)  *N_leftr(d). */

/*        Real workspace:    need  (N+M)*(M+P) + MAX(1,4*M). */
/*        Integer workspace: need  M. */

#line 565 "SB16BD.f"
	sb08gd_(ncr, p, m, &a[a_offset], lda, &dwork[kbe], &ldbe, &f[f_offset]
		, ldf, &dwork[kde], &ldde, &dwork[kbe + *n * *p], &ldbe, &
		dwork[kde + *m * *p], &ldde, &iwork[1], &dwork[kw], info);

/*        Copy the reduced system matrices Bc and Dc. */

#line 571 "SB16BD.f"
	dlacpy_("Full", ncr, p, &dwork[kbe], &ldbe, &g[g_offset], ldg, (
		ftnlen)4);
#line 572 "SB16BD.f"
	dlacpy_("Full", m, p, &dwork[kde], &ldde, &dc[dc_offset], lddc, (
		ftnlen)4);

#line 574 "SB16BD.f"
    } else {

/*        Form Ge(d) = [ N_right(d) ] */
/*                     [ M_right(d) ] as */

/*             ( A+B*F | G ) */
/*             (-----------) */
/*             (   F   | 0 ) */
/*             ( C+D*F | I ) */

/*        Real workspace:    need  (N+P)*(M+P). */
/*        Integer workspace: need  0. */

#line 587 "SB16BD.f"
	dgemm_("NoTranspose", "NoTranspose", n, n, m, &c_b22, &b[b_offset], 
		ldb, &f[f_offset], ldf, &c_b22, &a[a_offset], lda, (ftnlen)11,
		 (ftnlen)11);
#line 589 "SB16BD.f"
	kce = 1;
#line 590 "SB16BD.f"
	kde = kce + *n * (*p + *m);
#line 591 "SB16BD.f"
	ldce = *m + *p;
#line 592 "SB16BD.f"
	ldde = ldce;
#line 593 "SB16BD.f"
	dlacpy_("Full", m, n, &f[f_offset], ldf, &dwork[kce], &ldce, (ftnlen)
		4);
#line 594 "SB16BD.f"
	dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[kce + *m], &ldce, (
		ftnlen)4);
#line 595 "SB16BD.f"
	if (withd) {
#line 595 "SB16BD.f"
	    dgemm_("NoTranspose", "NoTranspose", p, n, m, &c_b22, &d__[
		    d_offset], ldd, &f[f_offset], ldf, &c_b22, &dwork[kce + *
		    m], &ldce, (ftnlen)11, (ftnlen)11);
#line 595 "SB16BD.f"
	}
#line 598 "SB16BD.f"
	dlaset_("Full", m, p, &c_b45, &c_b45, &dwork[kde], &ldde, (ftnlen)4);
#line 599 "SB16BD.f"
	dlaset_("Full", p, p, &c_b45, &c_b22, &dwork[kde + *m], &ldde, (
		ftnlen)4);

/*        Compute the reduced coprime factors, */
/*             Ger(d) = [ N_rightr(d) ] */
/*                      [ M_rightr(d) ], */
/*        by using either the B&T or SPA methods. */

/*        Real workspace:    need  (N+P)*(M+P) + */
/*                                 MAX(1,N*(2*N+MAX(N,M+P)+5)+N*(N+1)/2). */
/*        Integer workspace: need  0,         if JOBMR = 'B', */
/*                                 N,         if JOBMR = 'F', and */
/*                                 MAX(1,2*N) if JOBMR = 'S' or 'P'. */

#line 612 "SB16BD.f"
	kw = kde + *p * (*p + *m);
#line 613 "SB16BD.f"
	if (bta) {
#line 614 "SB16BD.f"
	    i__1 = *m + *p;
#line 614 "SB16BD.f"
	    i__2 = *ldwork - kw + 1;
#line 614 "SB16BD.f"
	    ab09ad_(dico, job, equil, ordsel, n, p, &i__1, ncr, &a[a_offset], 
		    lda, &g[g_offset], ldg, &dwork[kce], &ldce, &hsv[1], tol1,
		     &iwork[1], &dwork[kw], &i__2, iwarn, info, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 617 "SB16BD.f"
	} else {
#line 618 "SB16BD.f"
	    i__1 = *m + *p;
#line 618 "SB16BD.f"
	    i__2 = *ldwork - kw + 1;
#line 618 "SB16BD.f"
	    ab09bd_(dico, job, equil, ordsel, n, p, &i__1, ncr, &a[a_offset], 
		    lda, &g[g_offset], ldg, &dwork[kce], &ldce, &dwork[kde], &
		    ldde, &hsv[1], tol1, tol2, &iwork[1], &dwork[kw], &i__2, 
		    iwarn, info, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 622 "SB16BD.f"
	}
#line 623 "SB16BD.f"
	if (*info != 0) {
#line 624 "SB16BD.f"
	    if (*info != 3) {
#line 624 "SB16BD.f"
		*info += 3;
#line 624 "SB16BD.f"
	    }
#line 625 "SB16BD.f"
	    return 0;
#line 626 "SB16BD.f"
	}

#line 628 "SB16BD.f"
	wrkopt = (integer) dwork[kw] + kw - 1;

/*        Compute the reduced order controller, */
/*                                        -1 */
/*           Kr(d) = N_rightr(d)*M_rightr(d) . */

/*        Real workspace:    need  (N+P)*(M+P) + MAX(1,4*P). */
/*        Integer workspace: need  P. */

#line 637 "SB16BD.f"
	sb08hd_(ncr, p, m, &a[a_offset], lda, &g[g_offset], ldg, &dwork[kce], 
		&ldce, &dwork[kde], &ldde, &dwork[kce + *m], &ldce, &dwork[
		kde + *m], &ldde, &iwork[1], &dwork[kw], info);

/*        Copy the reduced system matrices Cc and Dc. */

#line 643 "SB16BD.f"
	dlacpy_("Full", m, ncr, &dwork[kce], &ldce, &f[f_offset], ldf, (
		ftnlen)4);
#line 644 "SB16BD.f"
	dlacpy_("Full", m, p, &dwork[kde], &ldde, &dc[dc_offset], lddc, (
		ftnlen)4);

#line 646 "SB16BD.f"
    }

#line 648 "SB16BD.f"
    dwork[1] = (doublereal) wrkopt;

#line 650 "SB16BD.f"
    return 0;
/* *** Last line of SB16BD *** */
} /* sb16bd_ */

