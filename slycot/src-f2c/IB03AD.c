#line 1 "IB03AD.f"
/* IB03AD.f -- translated by f2c (version 20100827).
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

#line 1 "IB03AD.f"
/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;
static doublereal c_b24 = -1.;
static integer c__0 = 0;
static integer c__5 = 5;

/* Subroutine */ int ib03ad_(char *init, char *alg, char *stor, integer *nobr,
	 integer *m, integer *l, integer *nsmp, integer *n, integer *nn, 
	integer *itmax1, integer *itmax2, integer *nprint, doublereal *u, 
	integer *ldu, doublereal *y, integer *ldy, doublereal *x, integer *lx,
	 doublereal *tol1, doublereal *tol2, integer *iwork, doublereal *
	dwork, integer *ldwork, integer *iwarn, integer *info, ftnlen 
	init_len, ftnlen alg_len, ftnlen stor_len)
{
    /* System generated locals */
    integer u_dim1, u_offset, y_dim1, y_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, z__, n2, ac, bd, ia, ib, ik, ml, iq, ir, is, iv, 
	    ix, ns, nx, iw1, iw2, ix0, ldr, bsn, mno, isv, iry, ldac, isad;
    static doublereal seed[4];
    static logical chol;
    static doublereal rcnd[16];
    static integer ipar[7], nfev, njev, lnol;
    static logical full;
    static integer nsml, lths, nths;
    static doublereal work[5];
    static logical init1, init2;
    extern /* Subroutine */ int ib01ad_(char *, char *, char *, char *, char *
	    , char *, integer *, integer *, integer *, integer *, doublereal *
	    , integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen, ftnlen, ftnlen), ib01bd_(char *, char *, char *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, logical *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen), ib01cd_(char *, char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    extern /* Subroutine */ int nf01ba_();
    extern /* Subroutine */ int md03ad_(char *, char *, char *, char *, U_fp, 
	    U_fp, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen);
    extern /* Subroutine */ int nf01bb_(), nf01bu_(), nf01bv_(), nf01bw_(), 
	    nf01bx_();
    static integer ircnd;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01vd_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen);
    static integer infol, lipar;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), tf01mx_(integer *, integer *, integer *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);
    static logical bwork[1];
    extern /* Subroutine */ int tb01vy_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen);
    static integer jwork, ircndb;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer iwarnl, wrkopt;


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

/*     To compute a set of parameters for approximating a Wiener system */
/*     in a least-squares sense, using a neural network approach and a */
/*     Levenberg-Marquardt algorithm. Conjugate gradients (CG) or */
/*     Cholesky algorithms are used to solve linear systems of equations. */
/*     The Wiener system is represented as */

/*        x(t+1) = A*x(t) + B*u(t) */
/*        z(t)   = C*x(t) + D*u(t), */

/*        y(t)   = f(z(t),wb(1:L)), */

/*     where t = 1, 2, ..., NSMP, and f is a nonlinear function, */
/*     evaluated by the SLICOT Library routine NF01AY. The parameter */
/*     vector X is partitioned as X = ( wb(1), ..., wb(L), theta ), */
/*     where wb(i), i = 1 : L, correspond to the nonlinear part, and */
/*     theta corresponds to the linear part. See SLICOT Library routine */
/*     NF01AD for further details. */

/*     The sum of squares of the error functions, defined by */

/*        e(t) = y(t) - Y(t),  t = 1, 2, ..., NSMP, */

/*     is minimized, where Y(t) is the measured output vector. The */
/*     functions and their Jacobian matrices are evaluated by SLICOT */
/*     Library routine NF01BB (the FCN routine in the call of MD03AD). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     INIT    CHARACTER*1 */
/*             Specifies which parts have to be initialized, as follows: */
/*             = 'L' : initialize the linear part only, X already */
/*                     contains an initial approximation of the */
/*                     nonlinearity; */
/*             = 'S' : initialize the static nonlinearity only, X */
/*                     already contains an initial approximation of the */
/*                     linear part; */
/*             = 'B' : initialize both linear and nonlinear parts; */
/*             = 'N' : do not initialize anything, X already contains */
/*                     an initial approximation. */
/*             If INIT = 'S' or 'B', the error functions for the */
/*             nonlinear part, and their Jacobian matrices, are evaluated */
/*             by SLICOT Library routine NF01BA (used as a second FCN */
/*             routine in the MD03AD call for the initialization step, */
/*             see METHOD). */

/*     ALG     CHARACTER*1 */
/*             Specifies the algorithm used for solving the linear */
/*             systems involving a Jacobian matrix J, as follows: */
/*             = 'D' :  a direct algorithm, which computes the Cholesky */
/*                      factor of the matrix J'*J + par*I is used, where */
/*                      par is the Levenberg factor; */
/*             = 'I' :  an iterative Conjugate Gradients algorithm, which */
/*                      only needs the matrix J, is used. */
/*             In both cases, matrix J is stored in a compressed form. */

/*     STOR    CHARACTER*1 */
/*             If ALG = 'D', specifies the storage scheme for the */
/*             symmetric matrix J'*J, as follows: */
/*             = 'F' :  full storage is used; */
/*             = 'P' :  packed storage is used. */
/*             The option STOR = 'F' usually ensures a faster execution. */
/*             This parameter is not relevant if ALG = 'I'. */

/*     Input/Output Parameters */

/*     NOBR    (input) INTEGER */
/*             If INIT = 'L' or 'B', NOBR is the number of block rows, s, */
/*             in the input and output block Hankel matrices to be */
/*             processed for estimating the linear part.  NOBR > 0. */
/*             (In the MOESP theory,  NOBR  should be larger than  n, */
/*             the estimated dimension of state vector.) */
/*             This parameter is ignored if INIT is 'S' or 'N'. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     L       (input) INTEGER */
/*             The number of system outputs.  L >= 0, and L > 0, if */
/*             INIT = 'L' or 'B'. */

/*     NSMP    (input) INTEGER */
/*             The number of input and output samples, t.  NSMP >= 0, and */
/*             NSMP >= 2*(M+L+1)*NOBR - 1, if INIT = 'L' or 'B'. */

/*     N       (input/output) INTEGER */
/*             The order of the linear part. */
/*             If INIT = 'L' or 'B', and N < 0 on entry, the order is */
/*             assumed unknown and it will be found by the routine. */
/*             Otherwise, the input value will be used. If INIT = 'S' */
/*             or 'N', N must be non-negative. The values N >= NOBR, */
/*             or N = 0, are not acceptable if INIT = 'L' or 'B'. */

/*     NN      (input) INTEGER */
/*             The number of neurons which shall be used to approximate */
/*             the nonlinear part.  NN >= 0. */

/*     ITMAX1  (input) INTEGER */
/*             The maximum number of iterations for the initialization of */
/*             the static nonlinearity. */
/*             This parameter is ignored if INIT is 'N' or 'L'. */
/*             Otherwise, ITMAX1 >= 0. */

/*     ITMAX2  (input) INTEGER */
/*             The maximum number of iterations.  ITMAX2 >= 0. */

/*     NPRINT  (input) INTEGER */
/*             This parameter enables controlled printing of iterates if */
/*             it is positive. In this case, FCN is called with IFLAG = 0 */
/*             at the beginning of the first iteration and every NPRINT */
/*             iterations thereafter and immediately prior to return, */
/*             and the current error norm is printed. Other intermediate */
/*             results could be printed by modifying the corresponding */
/*             FCN routine (NF01BA and/or NF01BB). If NPRINT <= 0, no */
/*             special calls of FCN with IFLAG = 0 are made. */

/*     U       (input) DOUBLE PRECISION array, dimension (LDU, M) */
/*             The leading NSMP-by-M part of this array must contain the */
/*             set of input samples, */
/*             U = ( U(1,1),...,U(1,M); ...; U(NSMP,1),...,U(NSMP,M) ). */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= MAX(1,NSMP). */

/*     Y       (input) DOUBLE PRECISION array, dimension (LDY, L) */
/*             The leading NSMP-by-L part of this array must contain the */
/*             set of output samples, */
/*             Y = ( Y(1,1),...,Y(1,L); ...; Y(NSMP,1),...,Y(NSMP,L) ). */

/*     LDY     INTEGER */
/*             The leading dimension of array Y.  LDY >= MAX(1,NSMP). */

/*     X       (input/output) DOUBLE PRECISION array dimension (LX) */
/*             On entry, if INIT = 'L', the leading (NN*(L+2) + 1)*L part */
/*             of this array must contain the initial parameters for */
/*             the nonlinear part of the system. */
/*             On entry, if INIT = 'S', the elements lin1 : lin2 of this */
/*             array must contain the initial parameters for the linear */
/*             part of the system, corresponding to the output normal */
/*             form, computed by SLICOT Library routine TB01VD, where */
/*                lin1 = (NN*(L+2) + 1)*L + 1; */
/*                lin2 = (NN*(L+2) + 1)*L + N*(L+M+1) + L*M. */
/*             On entry, if INIT = 'N', the elements 1 : lin2 of this */
/*             array must contain the initial parameters for the */
/*             nonlinear part followed by the initial parameters for the */
/*             linear part of the system, as specified above. */
/*             This array need not be set on entry if INIT = 'B'. */
/*             On exit, the elements 1 : lin2 of this array contain the */
/*             optimal parameters for the nonlinear part followed by the */
/*             optimal parameters for the linear part of the system, as */
/*             specified above. */

/*     LX      (input/output) INTEGER */
/*             On entry, this parameter must contain the intended length */
/*             of X. If N >= 0, then LX >= NX := lin2 (see parameter X). */
/*             If N is unknown (N < 0 on entry), a large enough estimate */
/*             of N should be used in the formula of lin2. */
/*             On exit, if N < 0 on entry, but LX is not large enough, */
/*             then this parameter contains the actual length of X, */
/*             corresponding to the computed N. Otherwise, its value */
/*             is unchanged. */

/*     Tolerances */

/*     TOL1    DOUBLE PRECISION */
/*             If INIT = 'S' or 'B' and TOL1 >= 0, TOL1 is the tolerance */
/*             which measures the relative error desired in the sum of */
/*             squares, for the initialization step of nonlinear part. */
/*             Termination occurs when the actual relative reduction in */
/*             the sum of squares is at most TOL1. In addition, if */
/*             ALG = 'I', TOL1 also measures the relative residual of */
/*             the solutions computed by the CG algorithm (for the */
/*             initialization step). Termination of a CG process occurs */
/*             when the relative residual is at most TOL1. */
/*             If the user sets  TOL1 < 0,  then  SQRT(EPS)  is used */
/*             instead TOL1, where EPS is the machine precision */
/*             (see LAPACK Library routine DLAMCH). */
/*             This parameter is ignored if INIT is 'N' or 'L'. */

/*     TOL2    DOUBLE PRECISION */
/*             If TOL2 >= 0, TOL2 is the tolerance which measures the */
/*             relative error desired in the sum of squares, for the */
/*             whole optimization process. Termination occurs when the */
/*             actual relative reduction in the sum of squares is at */
/*             most TOL2. */
/*             If ALG = 'I', TOL2 also measures the relative residual of */
/*             the solutions computed by the CG algorithm (for the whole */
/*             optimization). Termination of a CG process occurs when the */
/*             relative residual is at most TOL2. */
/*             If the user sets  TOL2 < 0,  then  SQRT(EPS)  is used */
/*             instead TOL2. This default value could require many */
/*             iterations, especially if TOL1 is larger. If INIT = 'S' */
/*             or 'B', it is advisable that TOL2 be larger than TOL1, */
/*             and spend more time with cheaper iterations. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (MAX( 3, LIW1, LIW2 )), where */
/*             LIW1 = LIW2 = 0,  if INIT = 'S' or 'N'; otherwise, */
/*             LIW1 = M+L; */
/*             LIW2 = MAX(M*NOBR+N,M*(N+L)). */
/*             On output, if INFO = 0, IWORK(1) and IWORK(2) return the */
/*             (total) number of function and Jacobian evaluations, */
/*             respectively (including the initialization step, if it was */
/*             performed), and if INIT = 'L' or INIT = 'B', IWORK(3) */
/*             specifies how many locations of DWORK contain reciprocal */
/*             condition number estimates (see below); otherwise, */
/*             IWORK(3) = 0. */

/*     DWORK   DOUBLE PRECISION array dimesion (LDWORK) */
/*             On entry, if desired, and if INIT = 'S' or 'B', the */
/*             entries DWORK(1:4) are set to initialize the random */
/*             numbers generator for the nonlinear part parameters (see */
/*             the description of the argument XINIT of SLICOT Library */
/*             routine MD03AD); this enables to obtain reproducible */
/*             results. The same seed is used for all outputs. */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK, DWORK(2) returns the residual error norm (the */
/*             sum of squares), DWORK(3) returns the number of iterations */
/*             performed, DWORK(4) returns the number of conjugate */
/*             gradients iterations performed, and DWORK(5) returns the */
/*             final Levenberg factor, for optimizing the parameters of */
/*             both the linear part and the static nonlinearity part. */
/*             If INIT = 'S' or INIT = 'B' and INFO = 0, then the */
/*             elements DWORK(6) to DWORK(10) contain the corresponding */
/*             five values for the initialization step (see METHOD). */
/*             (If L > 1, DWORK(10) contains the maximum of the Levenberg */
/*             factors for all outputs.) If INIT = 'L' or INIT = 'B', and */
/*             INFO = 0, DWORK(11) to DWORK(10+IWORK(3)) contain */
/*             reciprocal condition number estimates set by SLICOT */
/*             Library routines IB01AD, IB01BD, and IB01CD. */
/*             On exit, if  INFO = -23,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             In the formulas below, N should be taken not larger than */
/*             NOBR - 1, if N < 0 on entry. */
/*             LDWORK = MAX( LW1, LW2, LW3, LW4 ), where */
/*             LW1 = 0, if INIT = 'S' or 'N'; otherwise, */
/*             LW1 = MAX( 2*(M+L)*NOBR*(2*(M+L)*(NOBR+1)+3) + L*NOBR, */
/*                        4*(M+L)*NOBR*(M+L)*NOBR + (N+L)*(N+M) + */
/*                        MAX( LDW1, LDW2 ), */
/*                        (N+L)*(N+M) + N + N*N + 2 + N*(N+M+L) + */
/*                        MAX( 5*N, 2, MIN( LDW3, LDW4 ), LDW5, LDW6 ), */
/*                 where, */
/*                 LDW1 >= MAX( 2*(L*NOBR-L)*N+2*N, (L*NOBR-L)*N+N*N+7*N, */
/*                              L*NOBR*N + */
/*                              MAX( (L*NOBR-L)*N+2*N + (2*M+L)*NOBR+L, */
/*                                   2*(L*NOBR-L)*N+N*N+8*N, */
/*                                   N+4*(M*NOBR+N)+1, M*NOBR+3*N+L ) ) */
/*                 LDW2 >= 0,                                  if M = 0; */
/*                 LDW2 >= L*NOBR*N + M*NOBR*(N+L)*(M*(N+L)+1) + */
/*                         MAX( (N+L)**2, 4*M*(N+L)+1 ),       if M > 0; */
/*                 LDW3 = NSMP*L*(N+1) + 2*N + MAX( 2*N*N, 4*N ), */
/*                 LDW4 = N*(N+1) + 2*N + */
/*                        MAX( N*L*(N+1) + 2*N*N + L*N, 4*N ); */
/*                 LDW5 = NSMP*L + (N+L)*(N+M) + 3*N+M+L; */
/*                 LDW6 = NSMP*L + (N+L)*(N+M) + N + */
/*                        MAX(1, N*N*L + N*L + N, N*N + */
/*                            MAX(N*N + N*MAX(N,L) + 6*N + MIN(N,L), */
/*                                N*M)); */
/*             LW2 = LW3 = 0, if INIT = 'L' or 'N'; otherwise, */
/*             LW2 = NSMP*L + */
/*                   MAX( 5, NSMP + 2*BSN + NSMP*BSN + */
/*                           MAX( 2*NN + BSN, LDW7 ) ); */
/*                 LDW7 = BSN*BSN,       if ALG = 'D' and STOR = 'F'; */
/*                 LDW7 = BSN*(BSN+1)/2, if ALG = 'D' and STOR = 'P'; */
/*                 LDW7 = 3*BSN + NSMP,  if ALG = 'I'; */
/*             LW3 = MAX( LDW8, NSMP*L + (N+L)*(2*N+M) + 2*N ); */
/*                 LDW8 = NSMP*L + (N+L)*(N+M) + 3*N+M+L,  if M > 0; */
/*                 LDW8 = NSMP*L + (N+L)*N + 2*N+L,        if M = 0; */
/*             LW4 = MAX( 5, NSMP*L + 2*NX + NSMP*L*( BSN + LTHS ) + */
/*                           MAX( L1 + NX, NSMP*L + L1, L2 ) ), */
/*                  L0 = MAX( N*(N+L), N+M+L ),    if M > 0; */
/*                  L0 = MAX( N*(N+L), L ),        if M = 0; */
/*                  L1 = NSMP*L + MAX( 2*NN, (N+L)*(N+M) + 2*N + L0); */
/*                  L2 = NX*NX,          if ALG = 'D' and STOR = 'F'; */
/*                  L2 = NX*(NX+1)/2,    if ALG = 'D' and STOR = 'P'; */
/*                  L2 = 3*NX + NSMP*L,  if ALG = 'I', */
/*                  with BSN  = NN*( L + 2 ) + 1, */
/*                       LTHS = N*( L + M + 1 ) + L*M. */
/*             For optimum performance LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             < 0:  the user set IFLAG = IWARN in (one of) the */
/*                   subroutine(s) FCN, i.e., NF01BA, if INIT = 'S' */
/*                   or 'B', and/or NF01BB; this value cannot be returned */
/*                   without changing the FCN routine(s); */
/*                   otherwise, IWARN has the value k*100 + j*10 + i, */
/*                   where k is defined below, i refers to the whole */
/*                   optimization process, and j refers to the */
/*                   initialization step (j = 0, if INIT = 'L' or 'N'), */
/*                   and the possible values for i and j have the */
/*                   following meaning (where TOL* denotes TOL1 or TOL2, */
/*                   and similarly for ITMAX*): */
/*             = 1:  the number of iterations has reached ITMAX* without */
/*                   satisfying the convergence condition; */
/*             = 2:  if alg = 'I' and in an iteration of the Levenberg- */
/*                   Marquardt algorithm, the CG algorithm finished */
/*                   after 3*NX iterations (or 3*(lin1-1) iterations, for */
/*                   the initialization phase), without achieving the */
/*                   precision required in the call; */
/*             = 3:  the cosine of the angle between the vector of error */
/*                   function values and any column of the Jacobian is at */
/*                   most FACTOR*EPS in absolute value (FACTOR = 100); */
/*             = 4:  TOL* is too small: no further reduction in the sum */
/*                   of squares is possible. */
/*             The digit k is normally 0, but if INIT = 'L' or 'B', it */
/*             can have a value in the range 1 to 6 (see IB01AD, IB01BD */
/*             and IB01CD). In all these cases, the entries DWORK(1:5), */
/*             DWORK(6:10) (if INIT = 'S' or 'B'), and */
/*             DWORK(11:10+IWORK(3)) (if INIT = 'L' or 'B'), are set as */
/*             described above. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*                   otherwise, INFO has the value k*100 + j*10 + i, */
/*                   where k is defined below, i refers to the whole */
/*                   optimization process, and j refers to the */
/*                   initialization step (j = 0, if INIT = 'L' or 'N'), */
/*                   and the possible values for i and j have the */
/*                   following meaning: */
/*             = 1:  the routine FCN returned with INFO <> 0 for */
/*                   IFLAG = 1; */
/*             = 2:  the routine FCN returned with INFO <> 0 for */
/*                   IFLAG = 2; */
/*             = 3:  ALG = 'D' and SLICOT Library routines MB02XD or */
/*                   NF01BU (or NF01BV, if INIT = 'S' or 'B') or */
/*                   ALG = 'I' and SLICOT Library routines MB02WD or */
/*                   NF01BW (or NF01BX, if INIT = 'S' or 'B') returned */
/*                   with INFO <> 0. */
/*             In addition, if INIT = 'L' or 'B', i could also be */
/*             = 4:  if a Lyapunov equation could not be solved; */
/*             = 5:  if the identified linear system is unstable; */
/*             = 6:  if the QR algorithm failed on the state matrix */
/*                   of the identified linear system. */
/*             The digit k is normally 0, but if INIT = 'L' or 'B', it */
/*             can have a value in the range 1 to 10 (see IB01AD/IB01BD). */

/*     METHOD */

/*     If INIT = 'L' or 'B', the linear part of the system is */
/*     approximated using the combined MOESP and N4SID algorithm. If */
/*     necessary, this algorithm can also choose the order, but it is */
/*     advantageous if the order is already known. */

/*     If INIT = 'S' or 'B', the output of the approximated linear part */
/*     is computed and used to calculate an approximation of the static */
/*     nonlinearity using the Levenberg-Marquardt algorithm [1]. */
/*     This step is referred to as the (nonlinear) initialization step. */

/*     As last step, the Levenberg-Marquardt algorithm is used again to */
/*     optimize the parameters of the linear part and the static */
/*     nonlinearity as a whole. Therefore, it is necessary to parametrise */
/*     the matrices of the linear part. The output normal form [2] */
/*     parameterisation is used. */

/*     The Jacobian is computed analytically, for the nonlinear part, and */
/*     numerically, for the linear part. */

/*     REFERENCES */

/*     [1] Kelley, C.T. */
/*         Iterative Methods for Optimization. */
/*         Society for Industrial and Applied Mathematics (SIAM), */
/*         Philadelphia (Pa.), 1999. */

/*     [2] Peeters, R.L.M., Hanzon, B., and Olivi, M. */
/*         Balanced realizations of discrete-time stable all-pass */
/*         systems and the tangential Schur algorithm. */
/*         Proceedings of the European Control Conference, */
/*         31 August - 3 September 1999, Karlsruhe, Germany. */
/*         Session CP-6, Discrete-time Systems, 1999. */

/*     CONTRIBUTORS */

/*     A. Riedel, R. Schneider, Chemnitz University of Technology, */
/*     Oct. 2000, during a stay at University of Twente, NL. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001, */
/*     Mar. 2002, Apr. 2002, Feb. 2004, March 2005, Nov. 2005. */

/*     KEYWORDS */

/*     Conjugate gradients, least-squares approximation, */
/*     Levenberg-Marquardt algorithm, matrix operations, optimization. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     The upper triangular part is used in MD03AD; */
/*     For INIT = 'L' or 'B', additional parameters are set: */
/*     The following six parameters are used in the call of IB01AD; */
/*     The following three parameters are used in the call of IB01BD; */
/*     The following two parameters are used in the call of IB01CD; */
/*     TOLN controls the estimated order in IB01AD (default value); */
/*     RCOND controls the rank decisions in IB01AD, IB01BD, and IB01CD */
/*     (default); */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 485 "IB03AD.f"
    /* Parameter adjustments */
#line 485 "IB03AD.f"
    u_dim1 = *ldu;
#line 485 "IB03AD.f"
    u_offset = 1 + u_dim1;
#line 485 "IB03AD.f"
    u -= u_offset;
#line 485 "IB03AD.f"
    y_dim1 = *ldy;
#line 485 "IB03AD.f"
    y_offset = 1 + y_dim1;
#line 485 "IB03AD.f"
    y -= y_offset;
#line 485 "IB03AD.f"
    --x;
#line 485 "IB03AD.f"
    --iwork;
#line 485 "IB03AD.f"
    --dwork;
#line 485 "IB03AD.f"

#line 485 "IB03AD.f"
    /* Function Body */
#line 485 "IB03AD.f"
    chol = lsame_(alg, "D", (ftnlen)1, (ftnlen)1);
#line 486 "IB03AD.f"
    full = lsame_(stor, "F", (ftnlen)1, (ftnlen)1);
#line 487 "IB03AD.f"
    init1 = lsame_(init, "B", (ftnlen)1, (ftnlen)1) || lsame_(init, "L", (
	    ftnlen)1, (ftnlen)1);
#line 488 "IB03AD.f"
    init2 = lsame_(init, "B", (ftnlen)1, (ftnlen)1) || lsame_(init, "S", (
	    ftnlen)1, (ftnlen)1);

#line 490 "IB03AD.f"
    ml = *m + *l;
#line 491 "IB03AD.f"
    *info = 0;
#line 492 "IB03AD.f"
    *iwarn = 0;
#line 493 "IB03AD.f"
    if (! (init1 || init2 || lsame_(init, "N", (ftnlen)1, (ftnlen)1))) {
#line 494 "IB03AD.f"
	*info = -1;
#line 495 "IB03AD.f"
    } else if (! (chol || lsame_(alg, "I", (ftnlen)1, (ftnlen)1))) {
#line 496 "IB03AD.f"
	*info = -2;
#line 497 "IB03AD.f"
    } else if (chol && ! (full || lsame_(stor, "P", (ftnlen)1, (ftnlen)1))) {
#line 498 "IB03AD.f"
	*info = -3;
#line 499 "IB03AD.f"
    } else if (init1 && *nobr <= 0) {
#line 500 "IB03AD.f"
	*info = -4;
#line 501 "IB03AD.f"
    } else if (*m < 0) {
#line 502 "IB03AD.f"
	*info = -5;
#line 503 "IB03AD.f"
    } else if (*l < 0 || init1 && *l == 0) {
#line 504 "IB03AD.f"
	*info = -6;
#line 505 "IB03AD.f"
    } else if (*nsmp < 0 || init1 && *nsmp < (ml + 1 << 1) * *nobr - 1) {
#line 507 "IB03AD.f"
	*info = -7;
#line 508 "IB03AD.f"
    } else if (*n < 0 && ! init1 || (*n == 0 || *n >= *nobr) && init1) {
#line 510 "IB03AD.f"
	*info = -8;
#line 511 "IB03AD.f"
    } else if (*nn < 0) {
#line 512 "IB03AD.f"
	*info = -9;
#line 513 "IB03AD.f"
    } else if (init2 && *itmax1 < 0) {
#line 514 "IB03AD.f"
	*info = -10;
#line 515 "IB03AD.f"
    } else if (*itmax2 < 0) {
#line 516 "IB03AD.f"
	*info = -11;
#line 517 "IB03AD.f"
    } else if (*ldu < max(1,*nsmp)) {
#line 518 "IB03AD.f"
	*info = -14;
#line 519 "IB03AD.f"
    } else if (*ldy < max(1,*nsmp)) {
#line 520 "IB03AD.f"
	*info = -16;
#line 521 "IB03AD.f"
    } else {
#line 522 "IB03AD.f"
	lnol = *l * *nobr - *l;
#line 523 "IB03AD.f"
	mno = *m * *nobr;
#line 524 "IB03AD.f"
	bsn = *nn * (*l + 2) + 1;
#line 525 "IB03AD.f"
	nths = bsn * *l;
#line 526 "IB03AD.f"
	nsml = *nsmp * *l;
#line 527 "IB03AD.f"
	if (*n > 0) {
#line 528 "IB03AD.f"
	    ldac = *n + *l;
#line 529 "IB03AD.f"
	    isad = ldac * (*n + *m);
#line 530 "IB03AD.f"
	    n2 = *n * *n;
#line 531 "IB03AD.f"
	}

/*        Check the workspace size. */

#line 535 "IB03AD.f"
	jwork = 0;
#line 536 "IB03AD.f"
	if (init1) {
/*           Workspace for IB01AD. */
#line 538 "IB03AD.f"
	    jwork = (ml << 1) * *nobr * ((ml << 1) * (*nobr + 1) + 3) + *l * *
		    nobr;
#line 539 "IB03AD.f"
	    if (*n > 0) {
/*              Workspace for IB01BD. */
/* Computing MAX */
/* Computing MAX */
#line 541 "IB03AD.f"
		i__3 = lnol * *n + (*n << 1) + (*m + ml) * *nobr + *l, i__4 = 
			(lnol << 1) * *n + n2 + (*n << 3), i__3 = max(i__3,
			i__4), i__4 = *n + (mno + *n << 2) + 1, i__3 = max(
			i__3,i__4), i__4 = mno + *n * 3 + *l;
#line 541 "IB03AD.f"
		i__1 = (lnol << 1) * *n + (*n << 1), i__2 = lnol * *n + n2 + *
			n * 7, i__1 = max(i__1,i__2), i__2 = *l * *nobr * *n 
			+ max(i__3,i__4);
#line 541 "IB03AD.f"
		iw1 = max(i__1,i__2);
#line 545 "IB03AD.f"
		if (*m > 0) {
/* Computing MAX */
/* Computing 2nd power */
#line 546 "IB03AD.f"
		    i__3 = ldac;
#line 546 "IB03AD.f"
		    i__1 = i__3 * i__3, i__2 = (*m << 2) * ldac + 1;
#line 546 "IB03AD.f"
		    iw2 = *l * *nobr * *n + mno * ldac * (*m * ldac + 1) + 
			    max(i__1,i__2);
#line 548 "IB03AD.f"
		} else {
#line 549 "IB03AD.f"
		    iw2 = 0;
#line 550 "IB03AD.f"
		}
/* Computing MAX */
/* Computing 2nd power */
#line 551 "IB03AD.f"
		i__3 = (ml << 1) * *nobr;
#line 551 "IB03AD.f"
		i__1 = jwork, i__2 = i__3 * i__3 + isad + max(iw1,iw2);
#line 551 "IB03AD.f"
		jwork = max(i__1,i__2);
/*              Workspace for IB01CD. */
/* Computing MAX */
#line 554 "IB03AD.f"
		i__1 = n2 << 1, i__2 = *n << 2;
#line 554 "IB03AD.f"
		iw1 = nsml * (*n + 1) + (*n << 1) + max(i__1,i__2);
/* Computing MAX */
#line 555 "IB03AD.f"
		i__1 = *n * *l * (*n + 1) + (n2 << 1) + *l * *n, i__2 = *n << 
			2;
#line 555 "IB03AD.f"
		iw2 = *n * (*n + 1) + (*n << 1) + max(i__1,i__2);
/* Computing MAX */
/* Computing MAX */
#line 557 "IB03AD.f"
		i__3 = *n * 5, i__3 = max(i__3,2), i__4 = min(iw1,iw2);
#line 557 "IB03AD.f"
		i__1 = jwork, i__2 = isad + 2 + *n * (*n + 1 + ldac + *m) + 
			max(i__3,i__4);
#line 557 "IB03AD.f"
		jwork = max(i__1,i__2);
/*              Workspace for TF01MX. */
/* Computing MAX */
#line 560 "IB03AD.f"
		i__1 = jwork, i__2 = nsml + isad + ldac + (*n << 1) + *m;
#line 560 "IB03AD.f"
		jwork = max(i__1,i__2);
/*              Workspace for TB01VD. */
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
#line 562 "IB03AD.f"
		i__5 = n2 + *n * max(*n,*l) + *n * 6 + min(*n,*l), i__6 = *n *
			 *m;
#line 562 "IB03AD.f"
		i__3 = 1, i__4 = n2 * *l + *n * *l + *n, i__3 = max(i__3,i__4)
			, i__4 = n2 + max(i__5,i__6);
#line 562 "IB03AD.f"
		i__1 = jwork, i__2 = nsml + isad + *n + max(i__3,i__4);
#line 562 "IB03AD.f"
		jwork = max(i__1,i__2);
#line 566 "IB03AD.f"
	    }
#line 567 "IB03AD.f"
	}

#line 569 "IB03AD.f"
	if (init2) {
/*           Workspace for MD03AD (initialization of the nonlinear part). */
#line 571 "IB03AD.f"
	    if (chol) {
#line 572 "IB03AD.f"
		if (full) {
/* Computing 2nd power */
#line 573 "IB03AD.f"
		    i__1 = bsn;
#line 573 "IB03AD.f"
		    iw1 = i__1 * i__1;
#line 574 "IB03AD.f"
		} else {
#line 575 "IB03AD.f"
		    iw1 = bsn * (bsn + 1) / 2;
#line 576 "IB03AD.f"
		}
#line 577 "IB03AD.f"
	    } else {
#line 578 "IB03AD.f"
		iw1 = bsn * 3 + *nsmp;
#line 579 "IB03AD.f"
	    }
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
#line 580 "IB03AD.f"
	    i__5 = (*nn << 1) + bsn;
#line 580 "IB03AD.f"
	    i__3 = 5, i__4 = *nsmp + (bsn << 1) + *nsmp * bsn + max(i__5,iw1);
#line 580 "IB03AD.f"
	    i__1 = jwork, i__2 = nsml + max(i__3,i__4);
#line 580 "IB03AD.f"
	    jwork = max(i__1,i__2);
#line 583 "IB03AD.f"
	    if (*n > 0 && ! init1) {
/*              Workspace for TB01VY. */
/* Computing MAX */
#line 585 "IB03AD.f"
		i__1 = jwork, i__2 = nsml + ldac * ((*n << 1) + *m) + (*n << 
			1);
#line 585 "IB03AD.f"
		jwork = max(i__1,i__2);
/*              Workspace for TF01MX. */
#line 587 "IB03AD.f"
		if (*m > 0) {
#line 588 "IB03AD.f"
		    iw1 = *n + *m;
#line 589 "IB03AD.f"
		} else {
#line 590 "IB03AD.f"
		    iw1 = 0;
#line 591 "IB03AD.f"
		}
/* Computing MAX */
#line 592 "IB03AD.f"
		i__1 = jwork, i__2 = nsml + isad + iw1 + ldac + *n;
#line 592 "IB03AD.f"
		jwork = max(i__1,i__2);
#line 593 "IB03AD.f"
	    }
#line 594 "IB03AD.f"
	}

#line 596 "IB03AD.f"
	if (*n >= 0) {

/*           Find the number of parameters. */

#line 600 "IB03AD.f"
	    lths = *n * (ml + 1) + *l * *m;
#line 601 "IB03AD.f"
	    nx = nths + lths;

#line 603 "IB03AD.f"
	    if (*lx < nx) {
#line 604 "IB03AD.f"
		*info = -18;
#line 605 "IB03AD.f"
		i__1 = -(*info);
#line 605 "IB03AD.f"
		xerbla_("IB03AD", &i__1, (ftnlen)6);
#line 606 "IB03AD.f"
		return 0;
#line 607 "IB03AD.f"
	    }

/*           Workspace for MD03AD (whole optimization). */

#line 611 "IB03AD.f"
	    if (*m > 0) {
#line 612 "IB03AD.f"
		iw1 = ldac + *m;
#line 613 "IB03AD.f"
	    } else {
#line 614 "IB03AD.f"
		iw1 = *l;
#line 615 "IB03AD.f"
	    }
/* Computing MAX */
/* Computing MAX */
#line 616 "IB03AD.f"
	    i__3 = *n * ldac;
#line 616 "IB03AD.f"
	    i__1 = *nn << 1, i__2 = isad + (*n << 1) + max(i__3,iw1);
#line 616 "IB03AD.f"
	    iw1 = nsml + max(i__1,i__2);
#line 617 "IB03AD.f"
	    if (chol) {
#line 618 "IB03AD.f"
		if (full) {
/* Computing 2nd power */
#line 619 "IB03AD.f"
		    i__1 = nx;
#line 619 "IB03AD.f"
		    iw2 = i__1 * i__1;
#line 620 "IB03AD.f"
		} else {
#line 621 "IB03AD.f"
		    iw2 = nx * (nx + 1) / 2;
#line 622 "IB03AD.f"
		}
#line 623 "IB03AD.f"
	    } else {
#line 624 "IB03AD.f"
		iw2 = nx * 3 + nsml;
#line 625 "IB03AD.f"
	    }
/* Computing MAX */
/* Computing MAX */
#line 626 "IB03AD.f"
	    i__3 = iw1 + nx, i__4 = nsml + iw1, i__3 = max(i__3,i__4);
#line 626 "IB03AD.f"
	    i__1 = max(jwork,5), i__2 = nsml + (nx << 1) + nsml * (bsn + lths)
		     + max(i__3,iw2);
#line 626 "IB03AD.f"
	    jwork = max(i__1,i__2);
#line 629 "IB03AD.f"
	}

#line 631 "IB03AD.f"
	if (*ldwork < jwork) {
#line 632 "IB03AD.f"
	    *info = -23;
#line 633 "IB03AD.f"
	    dwork[1] = (doublereal) jwork;
#line 634 "IB03AD.f"
	}
#line 635 "IB03AD.f"
    }

#line 637 "IB03AD.f"
    if (*info != 0) {
#line 638 "IB03AD.f"
	i__1 = -(*info);
#line 638 "IB03AD.f"
	xerbla_("IB03AD", &i__1, (ftnlen)6);
#line 639 "IB03AD.f"
	return 0;
#line 640 "IB03AD.f"
    }

/*     Initialize the pointers to system matrices and save the possible */
/*     seed for random numbers generation. */

#line 645 "IB03AD.f"
    z__ = 1;
#line 646 "IB03AD.f"
    ac = z__ + nsml;
#line 647 "IB03AD.f"
    dcopy_(&c__4, &dwork[1], &c__1, seed, &c__1);

#line 649 "IB03AD.f"
    wrkopt = 1;

#line 651 "IB03AD.f"
    if (init1) {

/*        Initialize the linear part. */
/*        If N < 0, the order of the system is determined by IB01AD; */
/*        otherwise, the given order will be used. */
/*        The workspace needed is defined for the options set above */
/*        in the PARAMETER statements. */
/*        Workspace:  need:   2*(M+L)*NOBR*(2*(M+L)*(NOBR+1)+3) + L*NOBR; */
/*                    prefer: larger. */
/*        Integer workspace:  M+L. (If METH = 'N', (M+L)*NOBR.) */

#line 662 "IB03AD.f"
	ns = *n;
#line 663 "IB03AD.f"
	ir = 1;
#line 664 "IB03AD.f"
	isv = (ml << 1) * *nobr;
#line 665 "IB03AD.f"
	ldr = isv;
#line 666 "IB03AD.f"
	if (lsame_("N", "M", (ftnlen)1, (ftnlen)1)) {
/* Computing MAX */
#line 666 "IB03AD.f"
	    i__1 = ldr, i__2 = mno * 3;
#line 666 "IB03AD.f"
	    ldr = max(i__1,i__2);
#line 666 "IB03AD.f"
	}
#line 668 "IB03AD.f"
	isv = ir + ldr * isv;
#line 669 "IB03AD.f"
	jwork = isv + *l * *nobr;

#line 671 "IB03AD.f"
	i__1 = *ldwork - jwork + 1;
#line 671 "IB03AD.f"
	ib01ad_("M", "F", "N", "O", "N", "N", nobr, m, l, nsmp, &u[u_offset], 
		ldu, &y[y_offset], ldy, n, &dwork[ir], &ldr, &dwork[isv], &
		c_b24, &c_b24, &iwork[1], &dwork[jwork], &i__1, &iwarnl, &
		infol, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
		(ftnlen)1);

#line 676 "IB03AD.f"
	if (infol != 0) {
#line 677 "IB03AD.f"
	    *info = infol * 100;
#line 678 "IB03AD.f"
	    return 0;
#line 679 "IB03AD.f"
	}
#line 680 "IB03AD.f"
	if (iwarnl != 0) {
#line 680 "IB03AD.f"
	    *iwarn = iwarnl * 100;
#line 680 "IB03AD.f"
	}
/* Computing MAX */
#line 682 "IB03AD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 682 "IB03AD.f"
	wrkopt = max(i__1,i__2);
#line 683 "IB03AD.f"
	ircnd = 0;
#line 684 "IB03AD.f"
	if (lsame_("M", "N", (ftnlen)1, (ftnlen)1)) {
#line 685 "IB03AD.f"
	    ircnd = 2;
#line 686 "IB03AD.f"
	    dcopy_(&ircnd, &dwork[jwork + 1], &c__1, rcnd, &c__1);
#line 687 "IB03AD.f"
	}

#line 689 "IB03AD.f"
	if (ns >= 0) {
#line 690 "IB03AD.f"
	    *n = ns;
#line 691 "IB03AD.f"
	} else {

/*           Find the number of parameters. */

#line 695 "IB03AD.f"
	    ldac = *n + *l;
#line 696 "IB03AD.f"
	    isad = ldac * (*n + *m);
#line 697 "IB03AD.f"
	    n2 = *n * *n;
#line 698 "IB03AD.f"
	    lths = *n * (ml + 1) + *l * *m;
#line 699 "IB03AD.f"
	    nx = nths + lths;

#line 701 "IB03AD.f"
	    if (*lx < nx) {
#line 702 "IB03AD.f"
		*lx = nx;
#line 703 "IB03AD.f"
		*info = -18;
#line 704 "IB03AD.f"
		i__1 = -(*info);
#line 704 "IB03AD.f"
		xerbla_("IB03AD", &i__1, (ftnlen)6);
#line 705 "IB03AD.f"
		return 0;
#line 706 "IB03AD.f"
	    }
/*           Workspace for IB01BD. */
/* Computing MAX */
/* Computing MAX */
#line 708 "IB03AD.f"
	    i__3 = lnol * *n + (*n << 1) + (*m + ml) * *nobr + *l, i__4 = (
		    lnol << 1) * *n + n2 + (*n << 3), i__3 = max(i__3,i__4), 
		    i__4 = *n + (mno + *n << 2) + 1, i__3 = max(i__3,i__4), 
		    i__4 = mno + *n * 3 + *l;
#line 708 "IB03AD.f"
	    i__1 = (lnol << 1) * *n + (*n << 1), i__2 = lnol * *n + n2 + *n * 
		    7, i__1 = max(i__1,i__2), i__2 = *l * *nobr * *n + max(
		    i__3,i__4);
#line 708 "IB03AD.f"
	    iw1 = max(i__1,i__2);
#line 712 "IB03AD.f"
	    if (*m > 0) {
/* Computing MAX */
/* Computing 2nd power */
#line 713 "IB03AD.f"
		i__3 = ldac;
#line 713 "IB03AD.f"
		i__1 = i__3 * i__3, i__2 = (*m << 2) * ldac + 1;
#line 713 "IB03AD.f"
		iw2 = *l * *nobr * *n + mno * ldac * (*m * ldac + 1) + max(
			i__1,i__2);
#line 715 "IB03AD.f"
	    } else {
#line 716 "IB03AD.f"
		iw2 = 0;
#line 717 "IB03AD.f"
	    }
#line 718 "IB03AD.f"
	    jwork = isv + isad + max(iw1,iw2);
/*           Workspace for IB01CD. */
/* Computing MAX */
#line 720 "IB03AD.f"
	    i__1 = n2 << 1, i__2 = *n << 2;
#line 720 "IB03AD.f"
	    iw1 = nsml * (*n + 1) + (*n << 1) + max(i__1,i__2);
/* Computing MAX */
#line 721 "IB03AD.f"
	    i__1 = *n * *l * (*n + 1) + (n2 << 1) + *l * *n, i__2 = *n << 2;
#line 721 "IB03AD.f"
	    iw2 = *n * (*n + 1) + (*n << 1) + max(i__1,i__2);
/* Computing MAX */
/* Computing MAX */
#line 723 "IB03AD.f"
	    i__3 = *n * 5, i__3 = max(i__3,2), i__4 = min(iw1,iw2);
#line 723 "IB03AD.f"
	    i__1 = jwork, i__2 = isad + 2 + *n * (*n + 1 + ldac + *m) + max(
		    i__3,i__4);
#line 723 "IB03AD.f"
	    jwork = max(i__1,i__2);
/*           Workspace for TF01MX. */
/* Computing MAX */
#line 726 "IB03AD.f"
	    i__1 = jwork, i__2 = nsml + isad + ldac + (*n << 1) + *m;
#line 726 "IB03AD.f"
	    jwork = max(i__1,i__2);
/*           Workspace for TB01VD. */
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
#line 728 "IB03AD.f"
	    i__5 = n2 + *n * max(*n,*l) + *n * 6 + min(*n,*l), i__6 = *n * *m;
#line 728 "IB03AD.f"
	    i__3 = 1, i__4 = n2 * *l + *n * *l + *n, i__3 = max(i__3,i__4), 
		    i__4 = n2 + max(i__5,i__6);
#line 728 "IB03AD.f"
	    i__1 = jwork, i__2 = nsml + isad + *n + max(i__3,i__4);
#line 728 "IB03AD.f"
	    jwork = max(i__1,i__2);
/*           Workspace for MD03AD (whole optimization). */
#line 733 "IB03AD.f"
	    if (*m > 0) {
#line 734 "IB03AD.f"
		iw1 = ldac + *m;
#line 735 "IB03AD.f"
	    } else {
#line 736 "IB03AD.f"
		iw1 = *l;
#line 737 "IB03AD.f"
	    }
/* Computing MAX */
/* Computing MAX */
#line 738 "IB03AD.f"
	    i__3 = *n * ldac;
#line 738 "IB03AD.f"
	    i__1 = *nn << 1, i__2 = isad + (*n << 1) + max(i__3,iw1);
#line 738 "IB03AD.f"
	    iw1 = nsml + max(i__1,i__2);
#line 739 "IB03AD.f"
	    if (chol) {
#line 740 "IB03AD.f"
		if (full) {
/* Computing 2nd power */
#line 741 "IB03AD.f"
		    i__1 = nx;
#line 741 "IB03AD.f"
		    iw2 = i__1 * i__1;
#line 742 "IB03AD.f"
		} else {
#line 743 "IB03AD.f"
		    iw2 = nx * (nx + 1) / 2;
#line 744 "IB03AD.f"
		}
#line 745 "IB03AD.f"
	    } else {
#line 746 "IB03AD.f"
		iw2 = nx * 3 + nsml;
#line 747 "IB03AD.f"
	    }
/* Computing MAX */
/* Computing MAX */
#line 748 "IB03AD.f"
	    i__3 = iw1 + nx, i__4 = nsml + iw1, i__3 = max(i__3,i__4);
#line 748 "IB03AD.f"
	    i__1 = max(jwork,5), i__2 = nsml + (nx << 1) + nsml * (bsn + lths)
		     + max(i__3,iw2);
#line 748 "IB03AD.f"
	    jwork = max(i__1,i__2);
#line 751 "IB03AD.f"
	    if (*ldwork < jwork) {
#line 752 "IB03AD.f"
		*info = -23;
#line 753 "IB03AD.f"
		dwork[1] = (doublereal) jwork;
#line 754 "IB03AD.f"
		i__1 = -(*info);
#line 754 "IB03AD.f"
		xerbla_("IB03AD", &i__1, (ftnlen)6);
#line 755 "IB03AD.f"
		return 0;
#line 756 "IB03AD.f"
	    }
#line 757 "IB03AD.f"
	}

#line 759 "IB03AD.f"
	bd = ac + ldac * *n;
#line 760 "IB03AD.f"
	ix = bd + ldac * *m;
#line 761 "IB03AD.f"
	ia = isv;
#line 762 "IB03AD.f"
	ib = ia + ldac * *n;
#line 763 "IB03AD.f"
	iq = ib + ldac * *m;
#line 764 "IB03AD.f"
	if (lsame_("N", "N", (ftnlen)1, (ftnlen)1)) {
#line 765 "IB03AD.f"
	    iry = iq;
#line 766 "IB03AD.f"
	    is = iq;
#line 767 "IB03AD.f"
	    ik = iq;
#line 768 "IB03AD.f"
	    jwork = iq;
#line 769 "IB03AD.f"
	} else {
#line 770 "IB03AD.f"
	    iry = iq + n2;
#line 771 "IB03AD.f"
	    is = iry + *l * *l;
#line 772 "IB03AD.f"
	    ik = is + *n * *l;
#line 773 "IB03AD.f"
	    jwork = ik + *n * *l;
#line 774 "IB03AD.f"
	}

/*        The workspace needed is defined for the options set above */
/*        in the PARAMETER statements. */
/*        Workspace: */
/*          need:  4*(M+L)*NOBR*(M+L)*NOBR + (N+L)*(N+M) + */
/*                 max( LDW1,LDW2 ), where, */
/*                 LDW1 >= max( 2*(L*NOBR-L)*N+2*N, (L*NOBR-L)*N+N*N+7*N, */
/*                              L*NOBR*N + */
/*                              max( (L*NOBR-L)*N+2*N + (2*M+L)*NOBR+L, */
/*                                   2*(L*NOBR-L)*N+N*N+8*N, */
/*                                   N+4*(M*NOBR+N)+1, M*NOBR+3*N+L ) ) */
/*                 LDW2 >= 0,                                  if M = 0; */
/*                 LDW2 >= L*NOBR*N+M*NOBR*(N+L)*(M*(N+L)+1)+ */
/*                         max( (N+L)**2, 4*M*(N+L)+1 ),       if M > 0; */
/*          prefer: larger. */
/*        Integer workspace:  MAX(M*NOBR+N,M*(N+L)). */

#line 792 "IB03AD.f"
	i__1 = *ldwork - jwork + 1;
#line 792 "IB03AD.f"
	ib01bd_("C", "A", "N", nobr, n, m, l, nsmp, &dwork[ir], &ldr, &dwork[
		ia], &ldac, &dwork[ia + *n], &ldac, &dwork[ib], &ldac, &dwork[
		ib + *n], &ldac, &dwork[iq], n, &dwork[iry], l, &dwork[is], n,
		 &dwork[ik], n, &c_b24, &iwork[1], &dwork[jwork], &i__1, 
		bwork, &iwarnl, &infol, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 799 "IB03AD.f"
	if (infol == -30) {
#line 800 "IB03AD.f"
	    *info = -23;
#line 801 "IB03AD.f"
	    dwork[1] = dwork[jwork];
#line 802 "IB03AD.f"
	    i__1 = -(*info);
#line 802 "IB03AD.f"
	    xerbla_("IB03AD", &i__1, (ftnlen)6);
#line 803 "IB03AD.f"
	    return 0;
#line 804 "IB03AD.f"
	}
#line 805 "IB03AD.f"
	if (infol != 0) {
#line 806 "IB03AD.f"
	    *info = infol * 100;
#line 807 "IB03AD.f"
	    return 0;
#line 808 "IB03AD.f"
	}
#line 809 "IB03AD.f"
	if (iwarnl != 0) {
#line 809 "IB03AD.f"
	    *iwarn = iwarnl * 100;
#line 809 "IB03AD.f"
	}
/* Computing MAX */
#line 811 "IB03AD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 811 "IB03AD.f"
	wrkopt = max(i__1,i__2);
#line 812 "IB03AD.f"
	ircndb = 4;
#line 813 "IB03AD.f"
	if (lsame_("N", "K", (ftnlen)1, (ftnlen)1)) {
#line 813 "IB03AD.f"
	    ircndb += 8;
#line 813 "IB03AD.f"
	}
#line 815 "IB03AD.f"
	dcopy_(&ircndb, &dwork[jwork + 1], &c__1, &rcnd[ircnd], &c__1);
#line 816 "IB03AD.f"
	ircnd += ircndb;

/*        Copy the system matrices to the beginning of DWORK, to save */
/*        space, and redefine the pointers. */

#line 821 "IB03AD.f"
	dcopy_(&isad, &dwork[ia], &c__1, &dwork[1], &c__1);
#line 822 "IB03AD.f"
	ia = 1;
#line 823 "IB03AD.f"
	ib = ia + ldac * *n;
#line 824 "IB03AD.f"
	ix0 = ib + ldac * *m;
#line 825 "IB03AD.f"
	iv = ix0 + *n;

/*        Compute the initial condition of the system. On normal exit, */
/*           DWORK(i), i = JWORK+2:JWORK+1+N*N, */
/*           DWORK(j), j = JWORK+2+N*N:JWORK+1+N*N+L*N,  and */
/*           DWORK(k), k = JWORK+2+N*N+L*N:JWORK+1+N*N+L*N+N*M, */
/*        contain the transformed system matrices  At, Ct, and Bt, */
/*        respectively, corresponding to the real Schur form of the */
/*        estimated system state matrix  A. The transformation matrix is */
/*        stored in DWORK(IV:IV+N*N-1). */
/*        The workspace needed is defined for the options set above */
/*        in the PARAMETER statements. */
/*        Workspace: */
/*          need:   (N+L)*(N+M) + N + N*N + 2 + N*( N + M + L ) + */
/*                  max( 5*N, 2, min( LDW1, LDW2 ) ), where, */
/*                  LDW1 = NSMP*L*(N + 1) + 2*N + max( 2*N*N, 4*N), */
/*                  LDW2 = N*(N + 1) + 2*N + */
/*                         max( N*L*(N + 1) + 2*N*N + L*N, 4*N); */
/*          prefer: larger. */
/*        Integer workspace:  N. */

#line 846 "IB03AD.f"
	jwork = iv + n2;
#line 847 "IB03AD.f"
	i__1 = *ldwork - jwork + 1;
#line 847 "IB03AD.f"
	ib01cd_("X needed", "U", "D", n, m, l, nsmp, &dwork[ia], &ldac, &
		dwork[ib], &ldac, &dwork[ia + *n], &ldac, &dwork[ib + *n], &
		ldac, &u[u_offset], ldu, &y[y_offset], ldy, &dwork[ix0], &
		dwork[iv], n, &c_b24, &iwork[1], &dwork[jwork], &i__1, &
		iwarnl, &infol, (ftnlen)8, (ftnlen)1, (ftnlen)1);

#line 853 "IB03AD.f"
	if (infol == -26) {
#line 854 "IB03AD.f"
	    *info = -23;
#line 855 "IB03AD.f"
	    dwork[1] = dwork[jwork];
#line 856 "IB03AD.f"
	    i__1 = -(*info);
#line 856 "IB03AD.f"
	    xerbla_("IB03AD", &i__1, (ftnlen)6);
#line 857 "IB03AD.f"
	    return 0;
#line 858 "IB03AD.f"
	}
#line 859 "IB03AD.f"
	if (infol == 1) {
#line 859 "IB03AD.f"
	    infol = 10;
#line 859 "IB03AD.f"
	}
#line 861 "IB03AD.f"
	if (infol != 0) {
#line 862 "IB03AD.f"
	    *info = infol * 100;
#line 863 "IB03AD.f"
	    return 0;
#line 864 "IB03AD.f"
	}
#line 865 "IB03AD.f"
	if (iwarnl != 0) {
#line 865 "IB03AD.f"
	    *iwarn = iwarnl * 100;
#line 865 "IB03AD.f"
	}
/* Computing MAX */
#line 867 "IB03AD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 867 "IB03AD.f"
	wrkopt = max(i__1,i__2);
#line 868 "IB03AD.f"
	++ircnd;
#line 869 "IB03AD.f"
	rcnd[ircnd - 1] = dwork[jwork + 1];

/*        Now, save the system matrices and x0 in the final location. */

#line 873 "IB03AD.f"
	if (iv < ac) {
#line 874 "IB03AD.f"
	    i__1 = isad + *n;
#line 874 "IB03AD.f"
	    dcopy_(&i__1, &dwork[ia], &c__1, &dwork[ac], &c__1);
#line 875 "IB03AD.f"
	} else {
#line 876 "IB03AD.f"
	    i__1 = ac;
#line 876 "IB03AD.f"
	    for (j = ac + isad + *n - 1; j >= i__1; --j) {
#line 877 "IB03AD.f"
		dwork[j] = dwork[ia + j - ac];
#line 878 "IB03AD.f"
/* L5: */
#line 878 "IB03AD.f"
	    }
#line 879 "IB03AD.f"
	}

/*        Compute the output of the linear part. */
/*        Workspace: need   NSMP*L + (N + L)*(N + M) + 3*N + M + L, */
/*                                                              if M > 0; */
/*                          NSMP*L + (N + L)*N + 2*N + L,       if M = 0; */
/*                   prefer larger. */

#line 887 "IB03AD.f"
	jwork = ix + *n;
#line 888 "IB03AD.f"
	dcopy_(n, &dwork[ix], &c__1, &x[nths + 1], &c__1);
#line 889 "IB03AD.f"
	i__1 = *ldwork - jwork + 1;
#line 889 "IB03AD.f"
	tf01mx_(n, m, l, nsmp, &dwork[ac], &ldac, &u[u_offset], ldu, &x[nths 
		+ 1], &dwork[z__], nsmp, &dwork[jwork], &i__1, info);

/*        Convert the state-space representation to output normal form. */
/*        Workspace: */
/*          need:   NSMP*L + (N + L)*(N + M) + N + */
/*                  MAX(1, N*N*L + N*L + N, N*N + */
/*                      MAX(N*N + N*MAX(N,L) + 6*N + MIN(N,L), N*M)); */
/*          prefer: larger. */

#line 900 "IB03AD.f"
	i__1 = *ldwork - jwork + 1;
#line 900 "IB03AD.f"
	tb01vd_("Apply", n, m, l, &dwork[ac], &ldac, &dwork[bd], &ldac, &
		dwork[ac + *n], &ldac, &dwork[bd + *n], &ldac, &dwork[ix], &x[
		nths + 1], &lths, &dwork[jwork], &i__1, &infol, (ftnlen)5);

#line 905 "IB03AD.f"
	if (infol > 0) {
#line 906 "IB03AD.f"
	    *info = infol + 3;
#line 907 "IB03AD.f"
	    return 0;
#line 908 "IB03AD.f"
	}
/* Computing MAX */
#line 909 "IB03AD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 909 "IB03AD.f"
	wrkopt = max(i__1,i__2);

#line 911 "IB03AD.f"
    }

#line 913 "IB03AD.f"
    lipar = 7;
#line 914 "IB03AD.f"
    iw1 = 0;
#line 915 "IB03AD.f"
    iw2 = 0;

#line 917 "IB03AD.f"
    if (init2) {

/*        Initialize the nonlinear part. */

#line 921 "IB03AD.f"
	if (! init1) {
#line 922 "IB03AD.f"
	    bd = ac + ldac * *n;
#line 923 "IB03AD.f"
	    ix = bd + ldac * *m;

/*           Convert the output normal form to state-space model. */
/*           Workspace: need NSMP*L + (N + L)*(2*N + M) + 2*N. */
/*           (NSMP*L locations are reserved for the output of the linear */
/*           part.) */

#line 930 "IB03AD.f"
	    jwork = ix + *n;
#line 931 "IB03AD.f"
	    i__1 = *ldwork - jwork + 1;
#line 931 "IB03AD.f"
	    tb01vy_("Apply", n, m, l, &x[nths + 1], &lths, &dwork[ac], &ldac, 
		    &dwork[bd], &ldac, &dwork[ac + *n], &ldac, &dwork[bd + *n]
		    , &ldac, &dwork[ix], &dwork[jwork], &i__1, info, (ftnlen)
		    5);

/*           Compute the output of the linear part. */
/*           Workspace: need   NSMP*L + (N + L)*(N + M) + 3*N + M + L, */
/*                                                              if M > 0; */
/*                             NSMP*L + (N + L)*N + 2*N + L,    if M = 0; */
/*                      prefer larger. */

#line 942 "IB03AD.f"
	    i__1 = *ldwork - jwork + 1;
#line 942 "IB03AD.f"
	    tf01mx_(n, m, l, nsmp, &dwork[ac], &ldac, &u[u_offset], ldu, &
		    dwork[ix], &dwork[z__], nsmp, &dwork[jwork], &i__1, info);
#line 945 "IB03AD.f"
	}

/*        Optimize the parameters of the nonlinear part. */
/*        Workspace: */
/*          need   NSMP*L + */
/*                 MAX( 5, NSMP + 2*BSN + NSMP*BSN + */
/*                         MAX( 2*NN + BSN, DW( sol ) ) ), */
/*                 where, if ALG = 'D', */
/*                      DW( sol ) = BSN*BSN,        if STOR = 'F'; */
/*                      DW( sol ) = BSN*(BSN+1)/2,  if STOR = 'P'; */
/*                 and  DW( sol ) = 3*BSN + NSMP,   if ALG  = 'I'; */
/*          prefer larger. */

#line 958 "IB03AD.f"
	jwork = ac;
#line 959 "IB03AD.f"
	work[0] = 0.;
#line 960 "IB03AD.f"
	dcopy_(&c__4, work, &c__0, &work[1], &c__1);

/*        Set the integer parameters needed, including the number of */
/*        neurons. */

#line 965 "IB03AD.f"
	ipar[0] = *nsmp;
#line 966 "IB03AD.f"
	ipar[1] = *l;
#line 967 "IB03AD.f"
	ipar[2] = *nn;

#line 969 "IB03AD.f"
	i__1 = *l - 1;
#line 969 "IB03AD.f"
	for (i__ = 0; i__ <= i__1; ++i__) {
#line 970 "IB03AD.f"
	    dcopy_(&c__4, seed, &c__1, &dwork[jwork], &c__1);
#line 971 "IB03AD.f"
	    if (chol) {
#line 972 "IB03AD.f"
		i__2 = *ldwork - jwork + 1;
#line 972 "IB03AD.f"
		md03ad_("Random initialization", alg, stor, "U", (U_fp)
			nf01ba_, (U_fp)nf01bv_, nsmp, &bsn, itmax1, nprint, 
			ipar, &lipar, &dwork[z__], nsmp, &y[(i__ + 1) * 
			y_dim1 + 1], ldy, &x[i__ * bsn + 1], &nfev, &njev, 
			tol1, tol1, &dwork[jwork], &i__2, &iwarnl, &infol, (
			ftnlen)21, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 978 "IB03AD.f"
	    } else {
#line 979 "IB03AD.f"
		i__2 = *ldwork - jwork + 1;
#line 979 "IB03AD.f"
		md03ad_("Random initialization", alg, stor, "U", (U_fp)
			nf01ba_, (U_fp)nf01bx_, nsmp, &bsn, itmax1, nprint, 
			ipar, &lipar, &dwork[z__], nsmp, &y[(i__ + 1) * 
			y_dim1 + 1], ldy, &x[i__ * bsn + 1], &nfev, &njev, 
			tol1, tol1, &dwork[jwork], &i__2, &iwarnl, &infol, (
			ftnlen)21, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 985 "IB03AD.f"
	    }

#line 987 "IB03AD.f"
	    if (infol != 0) {
#line 988 "IB03AD.f"
		*info = infol * 10;
#line 989 "IB03AD.f"
		return 0;
#line 990 "IB03AD.f"
	    }
#line 991 "IB03AD.f"
	    if (iwarnl < 0) {
#line 992 "IB03AD.f"
		*info = infol;
#line 993 "IB03AD.f"
		*iwarn = iwarnl;
#line 994 "IB03AD.f"
		goto L20;
#line 995 "IB03AD.f"
	    } else if (iwarnl > 0) {
#line 996 "IB03AD.f"
		if (*iwarn > 100) {
/* Computing MAX */
#line 997 "IB03AD.f"
		    i__2 = *iwarn, i__3 = *iwarn / 100 * 100 + iwarnl * 10;
#line 997 "IB03AD.f"
		    *iwarn = max(i__2,i__3);
#line 998 "IB03AD.f"
		} else {
/* Computing MAX */
#line 999 "IB03AD.f"
		    i__2 = *iwarn, i__3 = iwarnl * 10;
#line 999 "IB03AD.f"
		    *iwarn = max(i__2,i__3);
#line 1000 "IB03AD.f"
		}
#line 1001 "IB03AD.f"
	    }
/* Computing MAX */
#line 1002 "IB03AD.f"
	    d__1 = work[0], d__2 = dwork[jwork];
#line 1002 "IB03AD.f"
	    work[0] = max(d__1,d__2);
/* Computing MAX */
#line 1003 "IB03AD.f"
	    d__1 = work[1], d__2 = dwork[jwork + 1];
#line 1003 "IB03AD.f"
	    work[1] = max(d__1,d__2);
/* Computing MAX */
#line 1004 "IB03AD.f"
	    d__1 = work[4], d__2 = dwork[jwork + 4];
#line 1004 "IB03AD.f"
	    work[4] = max(d__1,d__2);
#line 1005 "IB03AD.f"
	    work[2] += dwork[jwork + 2];
#line 1006 "IB03AD.f"
	    work[3] += dwork[jwork + 3];
#line 1007 "IB03AD.f"
	    iw1 = nfev + iw1;
#line 1008 "IB03AD.f"
	    iw2 = njev + iw2;
#line 1009 "IB03AD.f"
/* L10: */
#line 1009 "IB03AD.f"
	}

#line 1011 "IB03AD.f"
    }

/*     Main iteration. */
/*     Workspace: need   MAX( 5, NFUN + 2*NX + NFUN*( BSN + LTHS ) + */
/*                            MAX( LDW1 + NX, NFUN + LDW1, DW( sol ) ) ), */
/*                       where NFUN = NSMP*L, and */
/*                       LDW1 = NFUN + MAX( 2*NN, (N + L)*(N + M) + 2*N + */
/*                                          MAX( N*(N + L), N + M + L )), */
/*                                                              if M > 0, */
/*                       LDW1 = NFUN + MAX( 2*NN, (N + L)*N + 2*N + */
/*                                          MAX( N*(N + L), L ) ), */
/*                                                              if M = 0; */
/*                       if ALG = 'D', */
/*                             DW( sol ) = NX*NX,        if STOR = 'F'; */
/*                             DW( sol ) = NX*(NX+1)/2,  if STOR = 'P'; */
/*                       and   DW( sol ) = 3*NX + NFUN,  if ALG  = 'I', */
/*                       and DW( f ) is the workspace needed by the */
/*                       subroutine f; */
/*                prefer larger. */

/*     Set the integer parameters describing the Jacobian structure */
/*     and the number of neurons. */

#line 1034 "IB03AD.f"
    ipar[0] = lths;
#line 1035 "IB03AD.f"
    ipar[1] = *l;
#line 1036 "IB03AD.f"
    ipar[2] = *nsmp;
#line 1037 "IB03AD.f"
    ipar[3] = bsn;
#line 1038 "IB03AD.f"
    ipar[4] = *m;
#line 1039 "IB03AD.f"
    ipar[5] = *n;
#line 1040 "IB03AD.f"
    ipar[6] = *nn;

#line 1042 "IB03AD.f"
    if (chol) {
#line 1043 "IB03AD.f"
	md03ad_("Given initialization", alg, stor, "U", (U_fp)nf01bb_, (U_fp)
		nf01bu_, &nsml, &nx, itmax2, nprint, ipar, &lipar, &u[
		u_offset], ldu, &y[y_offset], ldy, &x[1], &nfev, &njev, tol2, 
		tol2, &dwork[1], ldwork, &iwarnl, info, (ftnlen)20, (ftnlen)1,
		 (ftnlen)1, (ftnlen)1);
#line 1047 "IB03AD.f"
    } else {
#line 1048 "IB03AD.f"
	md03ad_("Given initialization", alg, stor, "U", (U_fp)nf01bb_, (U_fp)
		nf01bw_, &nsml, &nx, itmax2, nprint, ipar, &lipar, &u[
		u_offset], ldu, &y[y_offset], ldy, &x[1], &nfev, &njev, tol2, 
		tol2, &dwork[1], ldwork, &iwarnl, info, (ftnlen)20, (ftnlen)1,
		 (ftnlen)1, (ftnlen)1);
#line 1052 "IB03AD.f"
    }

#line 1054 "IB03AD.f"
    if (*info != 0) {
#line 1054 "IB03AD.f"
	return 0;
#line 1054 "IB03AD.f"
    }

#line 1057 "IB03AD.f"
L20:
#line 1058 "IB03AD.f"
    iwork[1] = iw1 + nfev;
#line 1059 "IB03AD.f"
    iwork[2] = iw2 + njev;
#line 1060 "IB03AD.f"
    if (iwarnl < 0) {
#line 1061 "IB03AD.f"
	*iwarn = iwarnl;
#line 1062 "IB03AD.f"
    } else {
#line 1063 "IB03AD.f"
	*iwarn += iwarnl;
#line 1064 "IB03AD.f"
    }
#line 1065 "IB03AD.f"
    if (init2) {
#line 1065 "IB03AD.f"
	dcopy_(&c__5, work, &c__1, &dwork[6], &c__1);
#line 1065 "IB03AD.f"
    }
#line 1067 "IB03AD.f"
    if (init1) {
#line 1068 "IB03AD.f"
	iwork[3] = ircnd;
#line 1069 "IB03AD.f"
	dcopy_(&ircnd, rcnd, &c__1, &dwork[11], &c__1);
#line 1070 "IB03AD.f"
    } else {
#line 1071 "IB03AD.f"
	iwork[3] = 0;
#line 1072 "IB03AD.f"
    }
#line 1073 "IB03AD.f"
    return 0;

/* *** Last line of IB03AD *** */
} /* ib03ad_ */
