#line 1 "SB10AD.f"
/* SB10AD.f -- translated by f2c (version 20100827).
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

#line 1 "SB10AD.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int sb10ad_(integer *job, integer *n, integer *m, integer *
	np, integer *ncon, integer *nmeas, doublereal *gamma, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, doublereal *c__, integer *
	ldc, doublereal *d__, integer *ldd, doublereal *ak, integer *ldak, 
	doublereal *bk, integer *ldbk, doublereal *ck, integer *ldck, 
	doublereal *dk, integer *lddk, doublereal *ac, integer *ldac, 
	doublereal *bc, integer *ldbc, doublereal *cc, integer *ldcc, 
	doublereal *dc, integer *lddc, doublereal *rcond, doublereal *gtol, 
	doublereal *actol, integer *iwork, integer *liwork, doublereal *dwork,
	 integer *ldwork, logical *bwork, integer *lbwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ac_dim1, ac_offset, ak_dim1, ak_offset, b_dim1, 
	    b_offset, bc_dim1, bc_offset, bk_dim1, bk_offset, c_dim1, 
	    c_offset, cc_dim1, cc_offset, ck_dim1, ck_offset, d_dim1, 
	    d_offset, dc_dim1, dc_offset, dk_dim1, dk_offset, i__1, i__2, 
	    i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, m1, m2, m11, np1, np2, lw1, lw2, lw3, lw4, lw5, lw6, 
	    lw7, inf, np11, iwc, iwd, iwf, iwh, iwx, iwy, iwd1;
    static doublereal tol2;
    static integer iws1, iws2, iwac, mode, iwre, iwrk, iwwi, iwtu, iwwr, iwty,
	     info2, info3;
    extern /* Subroutine */ int sb10ld_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *), dgees_(
	    char *, char *, L_fp, integer *, doublereal *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, logical *, integer *, ftnlen, ftnlen), 
	    sb10pd_(integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *), sb10qd_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, logical *, integer *), 
	    sb10rd_(integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *);
    static doublereal gtoll, stepg;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal gamabs, mineac, gamamn, gamamx;
    extern /* Subroutine */ int dgesvd_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dlacpy_(char *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    extern logical select_();
    static integer lwamax, minwrk;


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

/*     To compute the matrices of an H-infinity optimal n-state */
/*     controller */

/*              | AK | BK | */
/*          K = |----|----|, */
/*              | CK | DK | */

/*     using modified Glover's and Doyle's 1988 formulas, for the system */

/*              | A  | B1  B2  |   | A | B | */
/*          P = |----|---------| = |---|---| */
/*              | C1 | D11 D12 |   | C | D | */
/*              | C2 | D21 D22 | */

/*     and for the estimated minimal possible value of gamma with respect */
/*     to GTOL, where B2 has as column size the number of control inputs */
/*     (NCON) and C2 has as row size the number of measurements (NMEAS) */
/*     being provided to the controller, and then to compute the matrices */
/*     of the closed-loop system */

/*              | AC | BC | */
/*          G = |----|----|, */
/*              | CC | DC | */

/*     if the stabilizing controller exists. */

/*     It is assumed that */

/*     (A1) (A,B2) is stabilizable and (C2,A) is detectable, */

/*     (A2) D12 is full column rank and D21 is full row rank, */

/*     (A3) | A-j*omega*I  B2  | has full column rank for all omega, */
/*          |    C1        D12 | */

/*     (A4) | A-j*omega*I  B1  |  has full row rank for all omega. */
/*          |    C2        D21 | */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     JOB     (input) INTEGER */
/*             Indicates the strategy for reducing the GAMMA value, as */
/*             follows: */
/*             = 1: Use bisection method for decreasing GAMMA from GAMMA */
/*                  to GAMMAMIN until the closed-loop system leaves */
/*                  stability. */
/*             = 2: Scan from GAMMA to 0 trying to find the minimal GAMMA */
/*                  for which the closed-loop system retains stability. */
/*             = 3: First bisection, then scanning. */
/*             = 4: Find suboptimal controller only. */

/*     N       (input) INTEGER */
/*             The order of the system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The column size of the matrix B.  M >= 0. */

/*     NP      (input) INTEGER */
/*             The row size of the matrix C.  NP >= 0. */

/*     NCON    (input) INTEGER */
/*             The number of control inputs (M2).  M >= NCON >= 0, */
/*             NP-NMEAS >= NCON. */

/*     NMEAS   (input) INTEGER */
/*             The number of measurements (NP2).  NP >= NMEAS >= 0, */
/*             M-NCON >= NMEAS. */

/*     GAMMA   (input/output) DOUBLE PRECISION */
/*             The initial value of gamma on input. It is assumed that */
/*             gamma is sufficiently large so that the controller is */
/*             admissible. GAMMA >= 0. */
/*             On output it contains the minimal estimated gamma. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             system state matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             system input matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading NP-by-N part of this array must contain the */
/*             system output matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= max(1,NP). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading NP-by-M part of this array must contain the */
/*             system input/output matrix D. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= max(1,NP). */

/*     AK      (output) DOUBLE PRECISION array, dimension (LDAK,N) */
/*             The leading N-by-N part of this array contains the */
/*             controller state matrix AK. */

/*     LDAK    INTEGER */
/*             The leading dimension of the array AK.  LDAK >= max(1,N). */

/*     BK      (output) DOUBLE PRECISION array, dimension (LDBK,NMEAS) */
/*             The leading N-by-NMEAS part of this array contains the */
/*             controller input matrix BK. */

/*     LDBK    INTEGER */
/*             The leading dimension of the array BK.  LDBK >= max(1,N). */

/*     CK      (output) DOUBLE PRECISION array, dimension (LDCK,N) */
/*             The leading NCON-by-N part of this array contains the */
/*             controller output matrix CK. */

/*     LDCK    INTEGER */
/*             The leading dimension of the array CK. */
/*             LDCK >= max(1,NCON). */

/*     DK      (output) DOUBLE PRECISION array, dimension (LDDK,NMEAS) */
/*             The leading NCON-by-NMEAS part of this array contains the */
/*             controller input/output matrix DK. */

/*     LDDK    INTEGER */
/*             The leading dimension of the array DK. */
/*             LDDK >= max(1,NCON). */

/*     AC      (output) DOUBLE PRECISION array, dimension (LDAC,2*N) */
/*             The leading 2*N-by-2*N part of this array contains the */
/*             closed-loop system state matrix AC. */

/*     LDAC    INTEGER */
/*             The leading dimension of the array AC. */
/*             LDAC >= max(1,2*N). */

/*     BC      (output) DOUBLE PRECISION array, dimension (LDBC,M-NCON) */
/*             The leading 2*N-by-(M-NCON) part of this array contains */
/*             the closed-loop system input matrix BC. */

/*     LDBC    INTEGER */
/*             The leading dimension of the array BC. */
/*             LDBC >= max(1,2*N). */

/*     CC      (output) DOUBLE PRECISION array, dimension (LDCC,2*N) */
/*             The leading (NP-NMEAS)-by-2*N part of this array contains */
/*             the closed-loop system output matrix CC. */

/*     LDCC    INTEGER */
/*             The leading dimension of the array CC. */
/*             LDCC >= max(1,NP-NMEAS). */

/*     DC      (output) DOUBLE PRECISION array, dimension (LDDC,M-NCON) */
/*             The leading (NP-NMEAS)-by-(M-NCON) part of this array */
/*             contains the closed-loop system input/output matrix DC. */

/*     LDDC    INTEGER */
/*             The leading dimension of the array DC. */
/*             LDDC >= max(1,NP-NMEAS). */

/*     RCOND   (output) DOUBLE PRECISION array, dimension (4) */
/*                      For the last successful step: */
/*             RCOND(1) contains the reciprocal condition number of the */
/*                      control transformation matrix; */
/*             RCOND(2) contains the reciprocal condition number of the */
/*                      measurement transformation matrix; */
/*             RCOND(3) contains an estimate of the reciprocal condition */
/*                      number of the X-Riccati equation; */
/*             RCOND(4) contains an estimate of the reciprocal condition */
/*                      number of the Y-Riccati equation. */

/*     Tolerances */

/*     GTOL    DOUBLE PRECISION */
/*             Tolerance used for controlling the accuracy of GAMMA */
/*             and its distance to the estimated minimal possible */
/*             value of GAMMA. */
/*             If GTOL <= 0, then a default value equal to sqrt(EPS) */
/*             is used, where EPS is the relative machine precision. */

/*     ACTOL   DOUBLE PRECISION */
/*             Upper bound for the poles of the closed-loop system */
/*             used for determining if it is stable. */
/*             ACTOL <= 0 for stable systems. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */

/*     LIWORK  INTEGER */
/*             The dimension of the array IWORK. */
/*             LIWORK >= max(2*max(N,M-NCON,NP-NMEAS,NCON,NMEAS),N*N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= LW1 + max(1,LW2,LW3,LW4,LW5 + MAX(LW6,LW7)), */
/*             where */
/*             LW1 = N*M + NP*N + NP*M + M2*M2 + NP2*NP2; */
/*             LW2 = max( ( N + NP1 + 1 )*( N + M2 ) + */
/*                          max( 3*( N + M2 ) + N + NP1, 5*( N + M2 ) ), */
/*                        ( N + NP2 )*( N + M1 + 1 ) + */
/*                          max( 3*( N + NP2 ) + N + M1, 5*( N + NP2 ) ), */
/*                        M2 + NP1*NP1 + max( NP1*max( N, M1 ), */
/*                                            3*M2 + NP1, 5*M2 ), */
/*                        NP2 + M1*M1 +  max( max( N, NP1 )*M1, */
/*                                            3*NP2 + M1, 5*NP2 ) ); */
/*             LW3 = max( ND1*M1 + max( 4*min( ND1, M1 ) + max( ND1,M1 ), */
/*                                      6*min( ND1, M1 ) ), */
/*                        NP1*ND2 + max( 4*min( NP1, ND2 ) + */
/*                                                        max( NP1,ND2 ), */
/*                                       6*min( NP1, ND2 ) ) ); */
/*             LW4 = 2*M*M + NP*NP + 2*M*N + M*NP + 2*N*NP; */
/*             LW5 = 2*N*N + M*N + N*NP; */
/*             LW6 = max( M*M   + max( 2*M1, 3*N*N + */
/*                                     max( N*M, 10*N*N + 12*N + 5 ) ), */
/*                        NP*NP + max( 2*NP1, 3*N*N + */
/*                                     max( N*NP, 10*N*N + 12*N + 5 ) )); */
/*             LW7 = M2*NP2 + NP2*NP2 + M2*M2 + */
/*                   max( ND1*ND1 + max( 2*ND1, ( ND1 + ND2 )*NP2 ), */
/*                        ND2*ND2 + max( 2*ND2, ND2*M2 ), 3*N, */
/*                        N*( 2*NP2 + M2 ) + */
/*                        max( 2*N*M2, M2*NP2 + */
/*                                     max( M2*M2 + 3*M2, NP2*( 2*NP2 + */
/*                                          M2 + max( NP2, N ) ) ) ) ); */
/*             M1  = M   - M2, NP1 = NP - NP2, */
/*             ND1 = NP1 - M2, ND2 = M1 - NP2. */
/*             For good performance, LDWORK must generally be larger. */

/*     BWORK   LOGICAL array, dimension (LBWORK) */

/*     LBWORK  INTEGER */
/*             The dimension of the array BWORK.  LBWORK >= 2*N. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the matrix | A-j*omega*I  B2  | had not full */
/*                                 |    C1        D12 | */
/*                   column rank in respect to the tolerance EPS; */
/*             = 2:  if the matrix | A-j*omega*I  B1  |  had not full row */
/*                                 |    C2        D21 | */
/*                   rank in respect to the tolerance EPS; */
/*             = 3:  if the matrix D12 had not full column rank in */
/*                   respect to the tolerance SQRT(EPS); */
/*             = 4:  if the matrix D21 had not full row rank in respect */
/*                   to the tolerance SQRT(EPS); */
/*             = 5:  if the singular value decomposition (SVD) algorithm */
/*                   did not converge (when computing the SVD of one of */
/*                   the matrices |A   B2 |, |A   B1 |, D12 or D21); */
/*                                |C1  D12|  |C2  D21| */
/*             = 6:  if the controller is not admissible (too small value */
/*                   of gamma); */
/*             = 7:  if the X-Riccati equation was not solved */
/*                   successfully (the controller is not admissible or */
/*                   there are numerical difficulties); */
/*             = 8:  if the Y-Riccati equation was not solved */
/*                   successfully (the controller is not admissible or */
/*                   there are numerical difficulties); */
/*             = 9:  if the determinant of Im2 + Tu*D11HAT*Ty*D22 is */
/*                   zero [3]; */
/*             = 10: if there are numerical problems when estimating */
/*                   singular values of D1111, D1112, D1111', D1121'; */
/*             = 11: if the matrices Inp2 - D22*DK or Im2 - DK*D22 */
/*                   are singular to working precision; */
/*             = 12: if a stabilizing controller cannot be found. */

/*     METHOD */

/*     The routine implements the Glover's and Doyle's 1988 formulas [1], */
/*     [2], modified to improve the efficiency as described in [3]. */

/*     JOB = 1: It tries with a decreasing value of GAMMA, starting with */
/*     the given, and with the newly obtained controller estimates of the */
/*     closed-loop system. If it is stable, (i.e., max(eig(AC)) < ACTOL) */
/*     the iterations can be continued until the given tolerance between */
/*     GAMMA and the estimated GAMMAMIN is reached. Otherwise, in the */
/*     next step GAMMA is increased. The step in the all next iterations */
/*     is step = step/2. The closed-loop system is obtained by the */
/*     formulas given in [2]. */

/*     JOB = 2: The same as for JOB = 1, but with non-varying step till */
/*     GAMMA = 0, step = max(0.1, GTOL). */

/*     JOB = 3: Combines the JOB = 1 and JOB = 2 cases for a quicker */
/*     procedure. */

/*     JOB = 4: Suboptimal controller for current GAMMA only. */

/*     REFERENCES */

/*     [1] Glover, K. and Doyle, J.C. */
/*         State-space formulae for all stabilizing controllers that */
/*         satisfy an Hinf norm bound and relations to risk sensitivity. */
/*         Systems and Control Letters, vol. 11, pp. 167-172, 1988. */

/*     [2] Balas, G.J., Doyle, J.C., Glover, K., Packard, A., and */
/*         Smith, R. */
/*         mu-Analysis and Synthesis Toolbox. */
/*         The MathWorks Inc., Natick, MA, 1995. */

/*     [3] Petkov, P.Hr., Gu, D.W., and Konstantinov, M.M. */
/*         Fortran 77 routines for Hinf and H2 design of continuous-time */
/*         linear control systems. */
/*         Rep. 98-14, Department of Engineering, Leicester University, */
/*         Leicester, U.K., 1998. */

/*     NUMERICAL ASPECTS */

/*     The accuracy of the result depends on the condition numbers of the */
/*     input and output transformations and on the condition numbers of */
/*     the two Riccati equations, as given by the values of RCOND(1), */
/*     RCOND(2), RCOND(3) and RCOND(4), respectively. */
/*     This approach by estimating the closed-loop system and checking */
/*     its poles seems to be reliable. */

/*     CONTRIBUTORS */

/*     A. Markovski, P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, */
/*     July 2003. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2003. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, H-infinity optimal control, robust */
/*     control. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and test input parameters. */

#line 416 "SB10AD.f"
    /* Parameter adjustments */
#line 416 "SB10AD.f"
    a_dim1 = *lda;
#line 416 "SB10AD.f"
    a_offset = 1 + a_dim1;
#line 416 "SB10AD.f"
    a -= a_offset;
#line 416 "SB10AD.f"
    b_dim1 = *ldb;
#line 416 "SB10AD.f"
    b_offset = 1 + b_dim1;
#line 416 "SB10AD.f"
    b -= b_offset;
#line 416 "SB10AD.f"
    c_dim1 = *ldc;
#line 416 "SB10AD.f"
    c_offset = 1 + c_dim1;
#line 416 "SB10AD.f"
    c__ -= c_offset;
#line 416 "SB10AD.f"
    d_dim1 = *ldd;
#line 416 "SB10AD.f"
    d_offset = 1 + d_dim1;
#line 416 "SB10AD.f"
    d__ -= d_offset;
#line 416 "SB10AD.f"
    ak_dim1 = *ldak;
#line 416 "SB10AD.f"
    ak_offset = 1 + ak_dim1;
#line 416 "SB10AD.f"
    ak -= ak_offset;
#line 416 "SB10AD.f"
    bk_dim1 = *ldbk;
#line 416 "SB10AD.f"
    bk_offset = 1 + bk_dim1;
#line 416 "SB10AD.f"
    bk -= bk_offset;
#line 416 "SB10AD.f"
    ck_dim1 = *ldck;
#line 416 "SB10AD.f"
    ck_offset = 1 + ck_dim1;
#line 416 "SB10AD.f"
    ck -= ck_offset;
#line 416 "SB10AD.f"
    dk_dim1 = *lddk;
#line 416 "SB10AD.f"
    dk_offset = 1 + dk_dim1;
#line 416 "SB10AD.f"
    dk -= dk_offset;
#line 416 "SB10AD.f"
    ac_dim1 = *ldac;
#line 416 "SB10AD.f"
    ac_offset = 1 + ac_dim1;
#line 416 "SB10AD.f"
    ac -= ac_offset;
#line 416 "SB10AD.f"
    bc_dim1 = *ldbc;
#line 416 "SB10AD.f"
    bc_offset = 1 + bc_dim1;
#line 416 "SB10AD.f"
    bc -= bc_offset;
#line 416 "SB10AD.f"
    cc_dim1 = *ldcc;
#line 416 "SB10AD.f"
    cc_offset = 1 + cc_dim1;
#line 416 "SB10AD.f"
    cc -= cc_offset;
#line 416 "SB10AD.f"
    dc_dim1 = *lddc;
#line 416 "SB10AD.f"
    dc_offset = 1 + dc_dim1;
#line 416 "SB10AD.f"
    dc -= dc_offset;
#line 416 "SB10AD.f"
    --rcond;
#line 416 "SB10AD.f"
    --iwork;
#line 416 "SB10AD.f"
    --dwork;
#line 416 "SB10AD.f"
    --bwork;
#line 416 "SB10AD.f"

#line 416 "SB10AD.f"
    /* Function Body */
#line 416 "SB10AD.f"
    m1 = *m - *ncon;
#line 417 "SB10AD.f"
    m2 = *ncon;
#line 418 "SB10AD.f"
    np1 = *np - *nmeas;
#line 419 "SB10AD.f"
    np2 = *nmeas;
#line 420 "SB10AD.f"
    np11 = np1 - m2;
#line 421 "SB10AD.f"
    m11 = m1 - np2;

#line 423 "SB10AD.f"
    *info = 0;
#line 424 "SB10AD.f"
    if (*job < 1 || *job > 4) {
#line 425 "SB10AD.f"
	*info = -1;
#line 426 "SB10AD.f"
    } else if (*n < 0) {
#line 427 "SB10AD.f"
	*info = -2;
#line 428 "SB10AD.f"
    } else if (*m < 0) {
#line 429 "SB10AD.f"
	*info = -3;
#line 430 "SB10AD.f"
    } else if (*np < 0) {
#line 431 "SB10AD.f"
	*info = -4;
#line 432 "SB10AD.f"
    } else if (*ncon < 0 || m1 < 0 || m2 > np1) {
#line 433 "SB10AD.f"
	*info = -5;
#line 434 "SB10AD.f"
    } else if (*nmeas < 0 || np1 < 0 || np2 > m1) {
#line 435 "SB10AD.f"
	*info = -6;
#line 436 "SB10AD.f"
    } else if (*gamma < 0.) {
#line 437 "SB10AD.f"
	*info = -7;
#line 438 "SB10AD.f"
    } else if (*lda < max(1,*n)) {
#line 439 "SB10AD.f"
	*info = -9;
#line 440 "SB10AD.f"
    } else if (*ldb < max(1,*n)) {
#line 441 "SB10AD.f"
	*info = -11;
#line 442 "SB10AD.f"
    } else if (*ldc < max(1,*np)) {
#line 443 "SB10AD.f"
	*info = -13;
#line 444 "SB10AD.f"
    } else if (*ldd < max(1,*np)) {
#line 445 "SB10AD.f"
	*info = -15;
#line 446 "SB10AD.f"
    } else if (*ldak < max(1,*n)) {
#line 447 "SB10AD.f"
	*info = -17;
#line 448 "SB10AD.f"
    } else if (*ldbk < max(1,*n)) {
#line 449 "SB10AD.f"
	*info = -19;
#line 450 "SB10AD.f"
    } else if (*ldck < max(1,m2)) {
#line 451 "SB10AD.f"
	*info = -21;
#line 452 "SB10AD.f"
    } else if (*lddk < max(1,m2)) {
#line 453 "SB10AD.f"
	*info = -23;
#line 454 "SB10AD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 454 "SB10AD.f"
	i__1 = 1, i__2 = *n << 1;
#line 454 "SB10AD.f"
	if (*ldac < max(i__1,i__2)) {
#line 455 "SB10AD.f"
	    *info = -25;
#line 456 "SB10AD.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 456 "SB10AD.f"
	    i__1 = 1, i__2 = *n << 1;
#line 456 "SB10AD.f"
	    if (*ldbc < max(i__1,i__2)) {
#line 457 "SB10AD.f"
		*info = -27;
#line 458 "SB10AD.f"
	    } else if (*ldcc < max(1,np1)) {
#line 459 "SB10AD.f"
		*info = -29;
#line 460 "SB10AD.f"
	    } else if (*lddc < max(1,np1)) {
#line 461 "SB10AD.f"
		*info = -31;
#line 462 "SB10AD.f"
	    } else {

/*        Compute workspace. */

#line 466 "SB10AD.f"
		lw1 = *n * *m + *np * *n + *np * *m + m2 * m2 + np2 * np2;
/* Computing MAX */
/* Computing MAX */
#line 467 "SB10AD.f"
		i__3 = (*n + m2) * 3 + *n + np1, i__4 = (*n + m2) * 5;
/* Computing MAX */
#line 467 "SB10AD.f"
		i__5 = (*n + np2) * 3 + *n + m1, i__6 = (*n + np2) * 5;
/* Computing MAX */
#line 467 "SB10AD.f"
		i__7 = np1 * max(*n,m1), i__8 = m2 * 3 + np1, i__7 = max(i__7,
			i__8), i__8 = m2 * 5;
/* Computing MAX */
#line 467 "SB10AD.f"
		i__9 = max(*n,np1) * m1, i__10 = np2 * 3 + m1, i__9 = max(
			i__9,i__10), i__10 = np2 * 5;
#line 467 "SB10AD.f"
		i__1 = (*n + np1 + 1) * (*n + m2) + max(i__3,i__4), i__2 = (*
			n + np2) * (*n + m1 + 1) + max(i__5,i__6), i__1 = max(
			i__1,i__2), i__2 = m2 + np1 * np1 + max(i__7,i__8), 
			i__1 = max(i__1,i__2), i__2 = np2 + m1 * m1 + max(
			i__9,i__10);
#line 467 "SB10AD.f"
		lw2 = max(i__1,i__2);
/* Computing MAX */
/* Computing MAX */
#line 475 "SB10AD.f"
		i__3 = (min(np11,m1) << 2) + max(np11,m1), i__4 = min(np11,m1)
			 * 6;
/* Computing MAX */
#line 475 "SB10AD.f"
		i__5 = (min(np1,m11) << 2) + max(np1,m11), i__6 = min(np1,m11)
			 * 6;
#line 475 "SB10AD.f"
		i__1 = np11 * m1 + max(i__3,i__4), i__2 = np1 * m11 + max(
			i__5,i__6);
#line 475 "SB10AD.f"
		lw3 = max(i__1,i__2);
#line 479 "SB10AD.f"
		lw4 = (*m << 1) * *m + *np * *np + (*m << 1) * *n + *m * *np 
			+ (*n << 1) * *np;
#line 480 "SB10AD.f"
		lw5 = (*n << 1) * *n + *m * *n + *n * *np;
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
#line 481 "SB10AD.f"
		i__5 = *n * *m, i__6 = *n * 10 * *n + *n * 12 + 5;
#line 481 "SB10AD.f"
		i__3 = m1 << 1, i__4 = *n * 3 * *n + max(i__5,i__6);
/* Computing MAX */
/* Computing MAX */
#line 481 "SB10AD.f"
		i__9 = *n * *np, i__10 = *n * 10 * *n + *n * 12 + 5;
#line 481 "SB10AD.f"
		i__7 = np1 << 1, i__8 = *n * 3 * *n + max(i__9,i__10);
#line 481 "SB10AD.f"
		i__1 = *m * *m + max(i__3,i__4), i__2 = *np * *np + max(i__7,
			i__8);
#line 481 "SB10AD.f"
		lw6 = max(i__1,i__2);
/* Computing MAX */
/* Computing MAX */
#line 485 "SB10AD.f"
		i__3 = np11 << 1, i__4 = (np11 + m11) * np2;
/* Computing MAX */
#line 485 "SB10AD.f"
		i__5 = m11 << 1, i__6 = m11 * m2;
/* Computing MAX */
/* Computing MAX */
#line 485 "SB10AD.f"
		i__9 = m2 * m2 + m2 * 3, i__10 = np2 * ((np2 << 1) + m2 + max(
			np2,*n));
#line 485 "SB10AD.f"
		i__7 = (*n << 1) * m2, i__8 = m2 * np2 + max(i__9,i__10);
#line 485 "SB10AD.f"
		i__1 = np11 * np11 + max(i__3,i__4), i__2 = m11 * m11 + max(
			i__5,i__6), i__1 = max(i__1,i__2), i__2 = *n * 3, 
			i__1 = max(i__1,i__2), i__2 = *n * ((np2 << 1) + m2) 
			+ max(i__7,i__8);
#line 485 "SB10AD.f"
		lw7 = m2 * np2 + np2 * np2 + m2 * m2 + max(i__1,i__2);
/* Computing MAX */
#line 492 "SB10AD.f"
		i__1 = max(1,lw2), i__1 = max(i__1,lw3), i__1 = max(i__1,lw4),
			 i__2 = lw5 + max(lw6,lw7);
#line 492 "SB10AD.f"
		minwrk = lw1 + max(i__1,i__2);
#line 493 "SB10AD.f"
		if (*ldwork < minwrk) {
#line 494 "SB10AD.f"
		    *info = -38;
#line 495 "SB10AD.f"
		} else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
#line 495 "SB10AD.f"
		    i__3 = max(*n,m1), i__3 = max(i__3,np1), i__3 = max(i__3,
			    m2);
#line 495 "SB10AD.f"
		    i__1 = max(i__3,np2) << 1, i__2 = *n * *n;
#line 495 "SB10AD.f"
		    if (*liwork < max(i__1,i__2)) {
#line 497 "SB10AD.f"
			*info = -36;
#line 498 "SB10AD.f"
		    } else if (*lbwork < *n << 1) {
#line 499 "SB10AD.f"
			*info = -40;
#line 500 "SB10AD.f"
		    }
#line 500 "SB10AD.f"
		}
#line 501 "SB10AD.f"
	    }
#line 501 "SB10AD.f"
	}
#line 501 "SB10AD.f"
    }

#line 503 "SB10AD.f"
    if (*info != 0) {

/*        Error return. */

#line 507 "SB10AD.f"
	i__1 = -(*info);
#line 507 "SB10AD.f"
	xerbla_("SB10AD", &i__1, (ftnlen)6);
#line 508 "SB10AD.f"
	return 0;
#line 509 "SB10AD.f"
    }

/*     Quick return if possible. */

#line 513 "SB10AD.f"
    if (*n == 0 || *m == 0 || *np == 0 || m1 == 0 || m2 == 0 || np1 == 0 || 
	    np2 == 0) {
#line 515 "SB10AD.f"
	rcond[1] = 1.;
#line 516 "SB10AD.f"
	rcond[2] = 1.;
#line 517 "SB10AD.f"
	rcond[3] = 1.;
#line 518 "SB10AD.f"
	rcond[4] = 1.;
#line 519 "SB10AD.f"
	dwork[1] = 1.;
#line 520 "SB10AD.f"
	return 0;
#line 521 "SB10AD.f"
    }

#line 523 "SB10AD.f"
    mode = *job;
#line 524 "SB10AD.f"
    if (mode > 2) {
#line 524 "SB10AD.f"
	mode = 1;
#line 524 "SB10AD.f"
    }
#line 526 "SB10AD.f"
    gtoll = *gtol;
#line 527 "SB10AD.f"
    if (gtoll <= 0.) {

/*        Set the default value of the tolerance for GAMMA. */

#line 531 "SB10AD.f"
	gtoll = sqrt(dlamch_("Epsilon", (ftnlen)7));
#line 532 "SB10AD.f"
    }

/*     Workspace usage 1. */

#line 536 "SB10AD.f"
    iwc = *n * *m + 1;
#line 537 "SB10AD.f"
    iwd = iwc + *np * *n;
#line 538 "SB10AD.f"
    iwtu = iwd + *np * *m;
#line 539 "SB10AD.f"
    iwty = iwtu + m2 * m2;
#line 540 "SB10AD.f"
    iwrk = iwty + np2 * np2;

#line 542 "SB10AD.f"
    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[1], n, (ftnlen)4);

#line 544 "SB10AD.f"
    dlacpy_("Full", np, n, &c__[c_offset], ldc, &dwork[iwc], np, (ftnlen)4);

#line 546 "SB10AD.f"
    dlacpy_("Full", np, m, &d__[d_offset], ldd, &dwork[iwd], np, (ftnlen)4);

/*     Transform the system so that D12 and D21 satisfy the formulas */
/*     in the computation of the Hinf optimal controller. */
/*     Workspace:  need   LW1 + MAX(1,LWP1,LWP2,LWP3,LWP4), */
/*                 prefer larger, */
/*                 where */
/*             LW1  = N*M + NP*N + NP*M + M2*M2 + NP2*NP2 */
/*             LWP1 = (N+NP1+1)*(N+M2) + MAX(3*(N+M2)+N+NP1,5*(N+M2)), */
/*             LWP2 = (N+NP2)*(N+M1+1) + MAX(3*(N+NP2)+N+M1,5*(N+NP2)), */
/*             LWP3 = M2 + NP1*NP1 + MAX(NP1*MAX(N,M1),3*M2+NP1,5*M2), */
/*             LWP4 = NP2 + M1*M1 + MAX(MAX(N,NP1)*M1,3*NP2+M1,5*NP2), */
/*             with M1 = M - M2 and NP1 = NP - NP2. */
/*             Denoting Q = MAX(M1,M2,NP1,NP2), an upper bound is */
/*             LW1 + MAX(1,(N+Q)*(N+Q+6),Q*(Q+MAX(N,Q,5)+1). */

#line 562 "SB10AD.f"
    tol2 = -1.;

#line 564 "SB10AD.f"
    i__1 = *ldwork - iwrk + 1;
#line 564 "SB10AD.f"
    sb10pd_(n, m, np, ncon, nmeas, &a[a_offset], lda, &dwork[1], n, &dwork[
	    iwc], np, &dwork[iwd], np, &dwork[iwtu], &m2, &dwork[iwty], &np2, 
	    &rcond[1], &tol2, &dwork[iwrk], &i__1, &info2);

#line 569 "SB10AD.f"
    lwamax = (integer) dwork[iwrk] + iwrk - 1;

#line 571 "SB10AD.f"
    if (info2 != 0) {
#line 572 "SB10AD.f"
	*info = info2;
#line 573 "SB10AD.f"
	return 0;
#line 574 "SB10AD.f"
    }

/*     Workspace usage 2. */

#line 578 "SB10AD.f"
    iwd1 = iwrk;
#line 579 "SB10AD.f"
    iws1 = iwd1 + np11 * m1;

/*     Check if GAMMA < max(sigma[D1111,D1112],sigma[D1111',D1121']). */
/*     Workspace:  need   LW1 + MAX(1, LWS1, LWS2), */
/*                 prefer larger, */
/*                 where */
/*     LWS1 = NP11*M1 + MAX(4*MIN(NP11,M1)+MAX(NP11,M1),6*MIN(NP11,M1)) */
/*     LWS2 = NP1*M11 + MAX(4*MIN(NP1,M11)+MAX(NP1,M11),6*MIN(NP1,M11)) */

#line 588 "SB10AD.f"
    info2 = 0;
#line 589 "SB10AD.f"
    info3 = 0;

#line 591 "SB10AD.f"
    if (np11 != 0 && m1 != 0) {
#line 592 "SB10AD.f"
	iwrk = iws1 + min(np11,m1);
#line 593 "SB10AD.f"
	dlacpy_("Full", &np11, &m1, &dwork[iwd], ldd, &dwork[iwd1], &np11, (
		ftnlen)4);
#line 595 "SB10AD.f"
	i__1 = *ldwork - iwrk + 1;
#line 595 "SB10AD.f"
	dgesvd_("N", "N", &np11, &m1, &dwork[iwd1], &np11, &dwork[iws1], &
		dwork[iws1], &c__1, &dwork[iws1], &c__1, &dwork[iwrk], &i__1, 
		&info2, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 598 "SB10AD.f"
	i__1 = lwamax, i__2 = (integer) dwork[iwrk] + iwrk - 1;
#line 598 "SB10AD.f"
	lwamax = max(i__1,i__2);
#line 599 "SB10AD.f"
    } else {
#line 600 "SB10AD.f"
	dwork[iws1] = 0.;
#line 601 "SB10AD.f"
    }

#line 603 "SB10AD.f"
    iws2 = iwd1 + np1 * m11;
#line 604 "SB10AD.f"
    if (np1 != 0 && m11 != 0) {
#line 605 "SB10AD.f"
	iwrk = iws2 + min(np1,m11);
#line 606 "SB10AD.f"
	dlacpy_("Full", &np1, &m11, &dwork[iwd], ldd, &dwork[iwd1], &np1, (
		ftnlen)4);
#line 608 "SB10AD.f"
	i__1 = *ldwork - iwrk + 1;
#line 608 "SB10AD.f"
	dgesvd_("N", "N", &np1, &m11, &dwork[iwd1], &np1, &dwork[iws2], &
		dwork[iws2], &c__1, &dwork[iws2], &c__1, &dwork[iwrk], &i__1, 
		&info3, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 611 "SB10AD.f"
	i__1 = lwamax, i__2 = (integer) dwork[iwrk] + iwrk - 1;
#line 611 "SB10AD.f"
	lwamax = max(i__1,i__2);
#line 612 "SB10AD.f"
    } else {
#line 613 "SB10AD.f"
	dwork[iws2] = 0.;
#line 614 "SB10AD.f"
    }

/* Computing MAX */
#line 616 "SB10AD.f"
    d__1 = dwork[iws1], d__2 = dwork[iws2];
#line 616 "SB10AD.f"
    gamamn = max(d__1,d__2);

#line 618 "SB10AD.f"
    if (info2 > 0 || info3 > 0) {
#line 619 "SB10AD.f"
	*info = 10;
#line 620 "SB10AD.f"
	return 0;
#line 621 "SB10AD.f"
    } else if (*gamma <= gamamn) {
#line 622 "SB10AD.f"
	*info = 6;
#line 623 "SB10AD.f"
	return 0;
#line 624 "SB10AD.f"
    }

/*     Workspace usage 3. */

#line 628 "SB10AD.f"
    iwx = iwd1;
#line 629 "SB10AD.f"
    iwy = iwx + *n * *n;
#line 630 "SB10AD.f"
    iwf = iwy + *n * *n;
#line 631 "SB10AD.f"
    iwh = iwf + *m * *n;
#line 632 "SB10AD.f"
    iwrk = iwh + *n * *np;
#line 633 "SB10AD.f"
    iwac = iwd1;
#line 634 "SB10AD.f"
    iwwr = iwac + (*n << 2) * *n;
#line 635 "SB10AD.f"
    iwwi = iwwr + (*n << 1);
#line 636 "SB10AD.f"
    iwre = iwwi + (*n << 1);

/*     Prepare some auxiliary variables for the gamma iteration. */

#line 640 "SB10AD.f"
    stepg = *gamma - gamamn;
#line 641 "SB10AD.f"
    gamabs = *gamma;
#line 642 "SB10AD.f"
    gamamx = *gamma;
#line 643 "SB10AD.f"
    inf = 0;

/*     ############################################################### */

/*     Begin the gamma iteration. */

#line 649 "SB10AD.f"
L10:
#line 650 "SB10AD.f"
    stepg /= 2.;

/*        Try to compute the state feedback and output injection */
/*        matrices for the current GAMMA. */

#line 655 "SB10AD.f"
    i__1 = *ldwork - iwrk + 1;
#line 655 "SB10AD.f"
    sb10qd_(n, m, np, ncon, nmeas, gamma, &a[a_offset], lda, &dwork[1], n, &
	    dwork[iwc], np, &dwork[iwd], np, &dwork[iwf], m, &dwork[iwh], n, &
	    dwork[iwx], n, &dwork[iwy], n, &rcond[3], &iwork[1], &dwork[iwrk],
	     &i__1, &bwork[1], &info2);

#line 661 "SB10AD.f"
    if (info2 != 0) {
#line 661 "SB10AD.f"
	goto L30;
#line 661 "SB10AD.f"
    }

/*        Try to compute the Hinf suboptimal (yet) controller. */

#line 665 "SB10AD.f"
    i__1 = *ldwork - iwrk + 1;
#line 665 "SB10AD.f"
    sb10rd_(n, m, np, ncon, nmeas, gamma, &a[a_offset], lda, &dwork[1], n, &
	    dwork[iwc], np, &dwork[iwd], np, &dwork[iwf], m, &dwork[iwh], n, &
	    dwork[iwtu], &m2, &dwork[iwty], &np2, &dwork[iwx], n, &dwork[iwy],
	     n, &ak[ak_offset], ldak, &bk[bk_offset], ldbk, &ck[ck_offset], 
	    ldck, &dk[dk_offset], lddk, &iwork[1], &dwork[iwrk], &i__1, &
	    info2);

#line 672 "SB10AD.f"
    if (info2 != 0) {
#line 672 "SB10AD.f"
	goto L30;
#line 672 "SB10AD.f"
    }

/*        Compute the closed-loop system. */
/*        Workspace: need   LW1 + 2*M*M + NP*NP + 2*M*N + M*NP + 2*N*NP; */
/*                   prefer larger. */

#line 678 "SB10AD.f"
    i__1 = *ldwork - iwd1 + 1;
#line 678 "SB10AD.f"
    sb10ld_(n, m, np, ncon, nmeas, &a[a_offset], lda, &b[b_offset], ldb, &c__[
	    c_offset], ldc, &d__[d_offset], ldd, &ak[ak_offset], ldak, &bk[
	    bk_offset], ldbk, &ck[ck_offset], ldck, &dk[dk_offset], lddk, &ac[
	    ac_offset], ldac, &bc[bc_offset], ldbc, &cc[cc_offset], ldcc, &dc[
	    dc_offset], lddc, &iwork[1], &dwork[iwd1], &i__1, &info2);

#line 683 "SB10AD.f"
    if (info2 != 0) {
#line 683 "SB10AD.f"
	goto L30;
#line 683 "SB10AD.f"
    }

/* Computing MAX */
#line 685 "SB10AD.f"
    i__1 = lwamax, i__2 = (integer) dwork[iwd1] + iwd1 - 1;
#line 685 "SB10AD.f"
    lwamax = max(i__1,i__2);

/*        Compute the poles of the closed-loop system. */
/*        Workspace:  need   LW1 + 4*N*N + 4*N + max(1,6*N); */
/*                    prefer larger. */

#line 691 "SB10AD.f"
    i__1 = *n << 1;
#line 691 "SB10AD.f"
    i__2 = *n << 1;
#line 691 "SB10AD.f"
    i__3 = *n << 1;
#line 691 "SB10AD.f"
    dlacpy_("Full", &i__1, &i__2, &ac[ac_offset], ldac, &dwork[iwac], &i__3, (
	    ftnlen)4);

#line 693 "SB10AD.f"
    i__1 = *n << 1;
#line 693 "SB10AD.f"
    i__2 = *n << 1;
#line 693 "SB10AD.f"
    i__3 = *ldwork - iwre + 1;
#line 693 "SB10AD.f"
    dgees_("N", "N", (L_fp)select_, &i__1, &dwork[iwac], &i__2, &iwork[1], &
	    dwork[iwwr], &dwork[iwwi], &dwork[iwre], &c__1, &dwork[iwre], &
	    i__3, &bwork[1], &info2, (ftnlen)1, (ftnlen)1);

/* Computing MAX */
#line 697 "SB10AD.f"
    i__1 = lwamax, i__2 = (integer) dwork[iwre] + iwre - 1;
#line 697 "SB10AD.f"
    lwamax = max(i__1,i__2);

/*        Now DWORK(IWWR+I)=Re(Lambda), DWORK(IWWI+I)=Im(Lambda), */
/*        for I=0,2*N-1. */

#line 702 "SB10AD.f"
    mineac = -1e3;

#line 704 "SB10AD.f"
    i__1 = (*n << 1) - 1;
#line 704 "SB10AD.f"
    for (i__ = 0; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 705 "SB10AD.f"
	d__1 = mineac, d__2 = dwork[iwwr + i__];
#line 705 "SB10AD.f"
	mineac = max(d__1,d__2);
#line 706 "SB10AD.f"
/* L20: */
#line 706 "SB10AD.f"
    }

/*        Check if the closed-loop system is stable. */

#line 710 "SB10AD.f"
L30:
#line 710 "SB10AD.f"
    if (mode == 1) {
#line 711 "SB10AD.f"
	if (info2 == 0 && mineac < *actol) {
#line 712 "SB10AD.f"
	    gamabs = *gamma;
#line 713 "SB10AD.f"
	    *gamma -= stepg;
#line 714 "SB10AD.f"
	    inf = 1;
#line 715 "SB10AD.f"
	} else {
/* Computing MIN */
#line 716 "SB10AD.f"
	    d__1 = *gamma + stepg;
#line 716 "SB10AD.f"
	    *gamma = min(d__1,gamamx);
#line 717 "SB10AD.f"
	}
#line 718 "SB10AD.f"
    } else if (mode == 2) {
#line 719 "SB10AD.f"
	if (info2 == 0 && mineac < *actol) {
#line 720 "SB10AD.f"
	    gamabs = *gamma;
#line 721 "SB10AD.f"
	    inf = 1;
#line 722 "SB10AD.f"
	}
#line 723 "SB10AD.f"
	*gamma -= max(.1,gtoll);
#line 724 "SB10AD.f"
    }

/*        More iterations? */

#line 728 "SB10AD.f"
    if (mode == 1 && *job == 3 && stepg * 2. < gtoll) {
#line 729 "SB10AD.f"
	mode = 2;
#line 730 "SB10AD.f"
	*gamma = gamabs;
#line 731 "SB10AD.f"
    }

#line 733 "SB10AD.f"
    if (*job != 4 && (mode == 1 && stepg * 2. >= gtoll || mode == 2 && *gamma 
	    > 0.)) {
#line 736 "SB10AD.f"
	goto L10;
#line 737 "SB10AD.f"
    }

/*     ############################################################### */

/*     End of the gamma iteration - Return if no stabilizing controller */
/*     was found. */

#line 744 "SB10AD.f"
    if (inf == 0) {
#line 745 "SB10AD.f"
	*info = 12;
#line 746 "SB10AD.f"
	return 0;
#line 747 "SB10AD.f"
    }

/*     Now compute the state feedback and output injection matrices */
/*     using GAMABS. */

#line 752 "SB10AD.f"
    *gamma = gamabs;

/*     Integer workspace:  need   max(2*max(N,M-NCON,NP-NMEAS),N*N). */
/*     Workspace: need   LW1P + */
/*                       max(1,M*M + max(2*M1,3*N*N + */
/*                                       max(N*M,10*N*N+12*N+5)), */
/*                           NP*NP + max(2*NP1,3*N*N + */
/*                                       max(N*NP,10*N*N+12*N+5))); */
/*                prefer larger, */
/*             where LW1P = LW1 + 2*N*N + M*N + N*NP. */
/*             An upper bound of the second term after LW1P is */
/*             max(1,4*Q*Q+max(2*Q,3*N*N + max(2*N*Q,10*N*N+12*N+5))). */

#line 765 "SB10AD.f"
    i__1 = *ldwork - iwrk + 1;
#line 765 "SB10AD.f"
    sb10qd_(n, m, np, ncon, nmeas, gamma, &a[a_offset], lda, &dwork[1], n, &
	    dwork[iwc], np, &dwork[iwd], np, &dwork[iwf], m, &dwork[iwh], n, &
	    dwork[iwx], n, &dwork[iwy], n, &rcond[3], &iwork[1], &dwork[iwrk],
	     &i__1, &bwork[1], &info2);

/* Computing MAX */
#line 771 "SB10AD.f"
    i__1 = lwamax, i__2 = (integer) dwork[iwrk] + iwrk - 1;
#line 771 "SB10AD.f"
    lwamax = max(i__1,i__2);

#line 773 "SB10AD.f"
    if (info2 > 0) {
#line 774 "SB10AD.f"
	*info = info2 + 5;
#line 775 "SB10AD.f"
	return 0;
#line 776 "SB10AD.f"
    }

/*     Compute the Hinf optimal controller. */
/*     Integer workspace:  need   max(2*(max(NP,M)-M2-NP2,M2,N),NP2). */
/*     Workspace: need   LW1P + */
/*                       max(1, M2*NP2 + NP2*NP2 + M2*M2 + */
/*                           max(D1*D1 + max(2*D1, (D1+D2)*NP2), */
/*                               D2*D2 + max(2*D2, D2*M2), 3*N, */
/*                               N*(2*NP2 + M2) + */
/*                               max(2*N*M2, M2*NP2 + */
/*                                           max(M2*M2+3*M2, NP2*(2*NP2+ */
/*                                                  M2+max(NP2,N)))))) */
/*                       where D1 = NP1 - M2 = NP11, D2 = M1 - NP2 = M11; */
/*                prefer larger. */
/*             An upper bound of the second term after LW1P is */
/*             max( 1, Q*(3*Q + 3*N + max(2*N, 4*Q + max(Q, N)))). */

#line 793 "SB10AD.f"
    i__1 = *ldwork - iwrk + 1;
#line 793 "SB10AD.f"
    sb10rd_(n, m, np, ncon, nmeas, gamma, &a[a_offset], lda, &dwork[1], n, &
	    dwork[iwc], np, &dwork[iwd], np, &dwork[iwf], m, &dwork[iwh], n, &
	    dwork[iwtu], &m2, &dwork[iwty], &np2, &dwork[iwx], n, &dwork[iwy],
	     n, &ak[ak_offset], ldak, &bk[bk_offset], ldbk, &ck[ck_offset], 
	    ldck, &dk[dk_offset], lddk, &iwork[1], &dwork[iwrk], &i__1, &
	    info2);

/* Computing MAX */
#line 800 "SB10AD.f"
    i__1 = lwamax, i__2 = (integer) dwork[iwrk] + iwrk - 1;
#line 800 "SB10AD.f"
    lwamax = max(i__1,i__2);

#line 802 "SB10AD.f"
    if (info2 == 1) {
#line 803 "SB10AD.f"
	*info = 6;
#line 804 "SB10AD.f"
	return 0;
#line 805 "SB10AD.f"
    } else if (info2 == 2) {
#line 806 "SB10AD.f"
	*info = 9;
#line 807 "SB10AD.f"
	return 0;
#line 808 "SB10AD.f"
    }

/*     Integer workspace:  need   2*max(NCON,NMEAS). */
/*     Workspace: need   2*M*M + NP*NP + 2*M*N + M*NP + 2*N*NP; */
/*                prefer larger. */

#line 814 "SB10AD.f"
    sb10ld_(n, m, np, ncon, nmeas, &a[a_offset], lda, &b[b_offset], ldb, &c__[
	    c_offset], ldc, &d__[d_offset], ldd, &ak[ak_offset], ldak, &bk[
	    bk_offset], ldbk, &ck[ck_offset], ldck, &dk[dk_offset], lddk, &ac[
	    ac_offset], ldac, &bc[bc_offset], ldbc, &cc[cc_offset], ldcc, &dc[
	    dc_offset], lddc, &iwork[1], &dwork[1], ldwork, &info2);

#line 819 "SB10AD.f"
    if (info2 > 0) {
#line 820 "SB10AD.f"
	*info = 11;
#line 821 "SB10AD.f"
	return 0;
#line 822 "SB10AD.f"
    }

#line 824 "SB10AD.f"
    dwork[1] = (doublereal) lwamax;
#line 825 "SB10AD.f"
    return 0;
/* *** Last line of SB10AD *** */
} /* sb10ad_ */

