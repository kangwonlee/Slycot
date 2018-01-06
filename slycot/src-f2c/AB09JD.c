#line 1 "AB09JD.f"
/* AB09JD.f -- translated by f2c (version 20100827).
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

#line 1 "AB09JD.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int ab09jd_(char *jobv, char *jobw, char *jobinv, char *dico,
	 char *equil, char *ordsel, integer *n, integer *nv, integer *nw, 
	integer *m, integer *p, integer *nr, doublereal *alpha, doublereal *a,
	 integer *lda, doublereal *b, integer *ldb, doublereal *c__, integer *
	ldc, doublereal *d__, integer *ldd, doublereal *av, integer *ldav, 
	doublereal *bv, integer *ldbv, doublereal *cv, integer *ldcv, 
	doublereal *dv, integer *lddv, doublereal *aw, integer *ldaw, 
	doublereal *bw, integer *ldbw, doublereal *cw, integer *ldcw, 
	doublereal *dw, integer *lddw, integer *ns, doublereal *hsv, 
	doublereal *tol1, doublereal *tol2, integer *iwork, doublereal *dwork,
	 integer *ldwork, integer *iwarn, integer *info, ftnlen jobv_len, 
	ftnlen jobw_len, ftnlen jobinv_len, ftnlen dico_len, ftnlen equil_len,
	 ftnlen ordsel_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, av_dim1, av_offset, aw_dim1, aw_offset, b_dim1, 
	    b_offset, bv_dim1, bv_offset, bw_dim1, bw_offset, c_dim1, 
	    c_offset, cv_dim1, cv_offset, cw_dim1, cw_offset, d_dim1, 
	    d_offset, dv_dim1, dv_offset, dw_dim1, dw_offset, i__1, i__2, 
	    i__3, i__4, i__5, i__6, i__7, i__8;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer ki, kl, ku, kw, lw, nu, nu1, nra, kav, kaw, kbv, kbw, kcv, 
	    kcw, kdv, kdw, kev, kew;
    static doublereal tol;
    static integer nwm, nvp, rank, ierr;
    static doublereal temp[1];
    extern /* Subroutine */ int ag07bd_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen), ab07nd_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *), ab08md_(char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen), tb01id_(char *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
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
    static integer ldabv, ldabw;
    extern /* Subroutine */ int ab09jv_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer ldcdv, ldcdw;
    extern /* Subroutine */ int ab09jw_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lefti, discr;
    static doublereal rcond;
    static char jobvl[1], jobwl[1];
    static logical conjv, conjw, leftw, invfr, autom;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static doublereal maxred;
    static logical righti, fixord;
    static integer iwarnl;
    static doublereal alpwrk;
    static logical frwght, rightw;
    static doublereal sqreps, wrkopt;


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
/*     state-space representation (A,B,C,D) by using the frequency */
/*     weighted optimal Hankel-norm approximation method. */
/*     The Hankel norm of the weighted error */

/*           op(V)*(G-Gr)*op(W) */

/*     is minimized, where G and Gr are the transfer-function matrices */
/*     of the original and reduced systems, respectively, V and W are */
/*     invertible transfer-function matrices representing the left and */
/*     right frequency weights, and op(X) denotes X, inv(X), conj(X) or */
/*     conj(inv(X)). V and W are specified by their state space */
/*     realizations (AV,BV,CV,DV) and (AW,BW,CW,DW), respectively. */
/*     When minimizing ||V*(G-Gr)*W||, V and W must be antistable. */
/*     When minimizing inv(V)*(G-Gr)*inv(W), V and W must have only */
/*     antistable zeros. */
/*     When minimizing conj(V)*(G-Gr)*conj(W), V and W must be stable. */
/*     When minimizing conj(inv(V))*(G-Gr)*conj(inv(W)), V and W must */
/*     be minimum-phase. */
/*     If the original system is unstable, then the frequency weighted */
/*     Hankel-norm approximation is computed only for the */
/*     ALPHA-stable part of the system. */

/*     For a transfer-function matrix G, conj(G) denotes the conjugate */
/*     of G given by G'(-s) for a continuous-time system or G'(1/z) */
/*     for a discrete-time system. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBV    CHARACTER*1 */
/*             Specifies the left frequency-weighting as follows: */
/*             = 'N':  V = I; */
/*             = 'V':  op(V) = V; */
/*             = 'I':  op(V) = inv(V); */
/*             = 'C':  op(V) = conj(V); */
/*             = 'R':  op(V) = conj(inv(V)). */

/*     JOBW    CHARACTER*1 */
/*             Specifies the right frequency-weighting as follows: */
/*             = 'N':  W = I; */
/*             = 'W':  op(W) = W; */
/*             = 'I':  op(W) = inv(W); */
/*             = 'C':  op(W) = conj(W); */
/*             = 'R':  op(W) = conj(inv(W)). */

/*     JOBINV  CHARACTER*1 */
/*             Specifies the computational approach to be used as */
/*             follows: */
/*             = 'N':  use the inverse free descriptor system approach; */
/*             = 'I':  use the inversion based standard approach; */
/*             = 'A':  switch automatically to the inverse free */
/*                     descriptor approach in case of badly conditioned */
/*                     feedthrough matrices in V or W (see METHOD). */

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
/*             The order of the original state-space representation, */
/*             i.e., the order of the matrix A.  N >= 0. */

/*     NV      (input) INTEGER */
/*             The order of the realization of the left frequency */
/*             weighting V, i.e., the order of the matrix AV.  NV >= 0. */

/*     NW      (input) INTEGER */
/*             The order of the realization of the right frequency */
/*             weighting W, i.e., the order of the matrix AW.  NW >= 0. */

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
/*             ALPHA-stable part of the weighted system (computed in */
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

/*     AV      (input/output) DOUBLE PRECISION array, dimension (LDAV,NV) */
/*             On entry, if JOBV <> 'N', the leading NV-by-NV part of */
/*             this array must contain the state matrix AV of a state */
/*             space realization of the left frequency weighting V. */
/*             On exit, if JOBV <> 'N', and INFO = 0, the leading */
/*             NV-by-NV part of this array contains the real Schur form */
/*             of AV. */
/*             AV is not referenced if JOBV = 'N'. */

/*     LDAV    INTEGER */
/*             The leading dimension of the array AV. */
/*             LDAV >= MAX(1,NV), if JOBV <> 'N'; */
/*             LDAV >= 1,         if JOBV =  'N'. */

/*     BV      (input/output) DOUBLE PRECISION array, dimension (LDBV,P) */
/*             On entry, if JOBV <> 'N', the leading NV-by-P part of */
/*             this array must contain the input matrix BV of a state */
/*             space realization of the left frequency weighting V. */
/*             On exit, if JOBV <> 'N', and INFO = 0, the leading */
/*             NV-by-P part of this array contains the transformed */
/*             input matrix BV corresponding to the transformed AV. */
/*             BV is not referenced if JOBV = 'N'. */

/*     LDBV    INTEGER */
/*             The leading dimension of the array BV. */
/*             LDBV >= MAX(1,NV), if JOBV <> 'N'; */
/*             LDBV >= 1,         if JOBV =  'N'. */

/*     CV      (input/output) DOUBLE PRECISION array, dimension (LDCV,NV) */
/*             On entry, if JOBV <> 'N', the leading P-by-NV part of */
/*             this array must contain the output matrix CV of a state */
/*             space realization of the left frequency weighting V. */
/*             On exit, if JOBV <> 'N', and INFO = 0, the leading */
/*             P-by-NV part of this array contains the transformed output */
/*             matrix CV corresponding to the transformed AV. */
/*             CV is not referenced if JOBV = 'N'. */

/*     LDCV    INTEGER */
/*             The leading dimension of the array CV. */
/*             LDCV >= MAX(1,P), if JOBV <> 'N'; */
/*             LDCV >= 1,        if JOBV =  'N'. */

/*     DV      (input) DOUBLE PRECISION array, dimension (LDDV,P) */
/*             If JOBV <> 'N', the leading P-by-P part of this array */
/*             must contain the feedthrough matrix DV of a state space */
/*             realization of the left frequency weighting V. */
/*             DV is not referenced if JOBV = 'N'. */

/*     LDDV    INTEGER */
/*             The leading dimension of the array DV. */
/*             LDDV >= MAX(1,P), if JOBV <> 'N'; */
/*             LDDV >= 1,        if JOBV =  'N'. */

/*     AW      (input/output) DOUBLE PRECISION array, dimension (LDAW,NW) */
/*             On entry, if JOBW <> 'N', the leading NW-by-NW part of */
/*             this array must contain the state matrix AW of a state */
/*             space realization of the right frequency weighting W. */
/*             On exit, if JOBW <> 'N', and INFO = 0, the leading */
/*             NW-by-NW part of this array contains the real Schur form */
/*             of AW. */
/*             AW is not referenced if JOBW = 'N'. */

/*     LDAW    INTEGER */
/*             The leading dimension of the array AW. */
/*             LDAW >= MAX(1,NW), if JOBW <> 'N'; */
/*             LDAW >= 1,         if JOBW =  'N'. */

/*     BW      (input/output) DOUBLE PRECISION array, dimension (LDBW,M) */
/*             On entry, if JOBW <> 'N', the leading NW-by-M part of */
/*             this array must contain the input matrix BW of a state */
/*             space realization of the right frequency weighting W. */
/*             On exit, if JOBW <> 'N', and INFO = 0, the leading */
/*             NW-by-M part of this array contains the transformed */
/*             input matrix BW corresponding to the transformed AW. */
/*             BW is not referenced if JOBW = 'N'. */

/*     LDBW    INTEGER */
/*             The leading dimension of the array BW. */
/*             LDBW >= MAX(1,NW), if JOBW <> 'N'; */
/*             LDBW >= 1,         if JOBW =  'N'. */

/*     CW      (input/output) DOUBLE PRECISION array, dimension (LDCW,NW) */
/*             On entry, if JOBW <> 'N', the leading M-by-NW part of */
/*             this array must contain the output matrix CW of a state */
/*             space realization of the right frequency weighting W. */
/*             On exit, if JOBW <> 'N', and INFO = 0, the leading */
/*             M-by-NW part of this array contains the transformed output */
/*             matrix CW corresponding to the transformed AW. */
/*             CW is not referenced if JOBW = 'N'. */

/*     LDCW    INTEGER */
/*             The leading dimension of the array CW. */
/*             LDCW >= MAX(1,M), if JOBW <> 'N'; */
/*             LDCW >= 1,        if JOBW =  'N'. */

/*     DW      (input) DOUBLE PRECISION array, dimension (LDDW,M) */
/*             If JOBW <> 'N', the leading M-by-M part of this array */
/*             must contain the feedthrough matrix DW of a state space */
/*             realization of the right frequency weighting W. */
/*             DW is not referenced if JOBW = 'N'. */

/*     LDDW    INTEGER */
/*             The leading dimension of the array DW. */
/*             LDDW >= MAX(1,M), if JOBW <> 'N'; */
/*             LDDW >= 1,        if JOBW =  'N'. */

/*     NS      (output) INTEGER */
/*             The dimension of the ALPHA-stable subsystem. */

/*     HSV     (output) DOUBLE PRECISION array, dimension (N) */
/*             If INFO = 0, the leading NS elements of this array contain */
/*             the Hankel singular values, ordered decreasingly, of the */
/*             projection G1s of op(V)*G1*op(W) (see METHOD), where G1 */
/*             is the ALPHA-stable part of the original system. */

/*     Tolerances */

/*     TOL1    DOUBLE PRECISION */
/*             If ORDSEL = 'A', TOL1 contains the tolerance for */
/*             determining the order of reduced system. */
/*             For model reduction, the recommended value is */
/*             TOL1 = c*HNORM(G1s), where c is a constant in the */
/*             interval [0.00001,0.001], and HNORM(G1s) is the */
/*             Hankel-norm of the projection G1s of op(V)*G1*op(W) */
/*             (see METHOD), computed in HSV(1). */
/*             If TOL1 <= 0 on entry, the used default value is */
/*             TOL1 = NS*EPS*HNORM(G1s), where NS is the number of */
/*             ALPHA-stable eigenvalues of A and EPS is the machine */
/*             precision (see LAPACK Library Routine DLAMCH). */
/*             If ORDSEL = 'F', the value of TOL1 is ignored. */
/*             TOL1 < 1. */

/*     TOL2    DOUBLE PRECISION */
/*             The tolerance for determining the order of a minimal */
/*             realization of the ALPHA-stable part of the given system. */
/*             The recommended value is TOL2 = NS*EPS*HNORM(G1s). */
/*             This value is used by default if TOL2 <= 0 on entry. */
/*             If TOL2 > 0 and ORDSEL = 'A', then TOL2 <= TOL1. */
/*             TOL2 < 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK = MAX(1,M,c,d),    if DICO = 'C', */
/*             LIWORK = MAX(1,N,M,c,d),  if DICO = 'D', where */
/*                c = 0,                          if JOBV =  'N', */
/*                c = MAX(2*P,NV+P+N+6,2*NV+P+2), if JOBV <> 'N', */
/*                d = 0,                          if JOBW =  'N', */
/*                d = MAX(2*M,NW+M+N+6,2*NW+M+2), if JOBW <> 'N'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( LDW1, LDW2, LDW3, LDW4 ), where */
/*             for NVP = NV+P and NWM = NW+M we have */
/*             LDW1 = 0 if JOBV =  'N' and */
/*             LDW1 = 2*NVP*(NVP+P) + P*P + */
/*                    MAX( 2*NVP*NVP + MAX( 11*NVP+16, P*NVP ), */
/*                          NVP*N + MAX( NVP*N+N*N, P*N, P*M ) ) */
/*                      if JOBV <> 'N', */
/*             LDW2 = 0 if JOBW =  'N' and */
/*             LDW2 = 2*NWM*(NWM+M) + M*M + */
/*                    MAX( 2*NWM*NWM + MAX( 11*NWM+16, M*NWM ), */
/*                          NWM*N + MAX( NWM*N+N*N, M*N, P*M ) ) */
/*                      if JOBW <> 'N', */
/*             LDW3 = N*(2*N + MAX(N,M,P) + 5) + N*(N+1)/2, */
/*             LDW4 = N*(M+P+2) + 2*M*P + MIN(N,M) + */
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
/*             =  0:  successful exit; */
/*             <  0:  if INFO = -i, the i-th argument had an illegal */
/*                    value; */
/*             =  1:  the computation of the ordered real Schur form of A */
/*                    failed; */
/*             =  2:  the separation of the ALPHA-stable/unstable */
/*                    diagonal blocks failed because of very close */
/*                    eigenvalues; */
/*             =  3:  the reduction of AV to a real Schur form failed; */
/*             =  4:  the reduction of AW to a real Schur form failed; */
/*             =  5:  the reduction to generalized Schur form of the */
/*                    descriptor pair corresponding to the inverse of V */
/*                    failed; */
/*             =  6:  the reduction to generalized Schur form of the */
/*                    descriptor pair corresponding to the inverse of W */
/*                    failed; */
/*             =  7:  the computation of Hankel singular values failed; */
/*             =  8:  the computation of stable projection in the */
/*                    Hankel-norm approximation algorithm failed; */
/*             =  9:  the order of computed stable projection in the */
/*                    Hankel-norm approximation algorithm differs */
/*                    from the order of Hankel-norm approximation; */
/*             = 10:  the reduction of AV-BV*inv(DV)*CV to a */
/*                    real Schur form failed; */
/*             = 11:  the reduction of AW-BW*inv(DW)*CW to a */
/*                    real Schur form failed; */
/*             = 12:  the solution of the Sylvester equation failed */
/*                    because the poles of V (if JOBV = 'V') or of */
/*                    conj(V) (if JOBV = 'C') are not distinct from */
/*                    the poles of G1 (see METHOD); */
/*             = 13:  the solution of the Sylvester equation failed */
/*                    because the poles of W (if JOBW = 'W') or of */
/*                    conj(W) (if JOBW = 'C') are not distinct from */
/*                    the poles of G1 (see METHOD); */
/*             = 14:  the solution of the Sylvester equation failed */
/*                    because the zeros of V (if JOBV = 'I') or of */
/*                    conj(V) (if JOBV = 'R') are not distinct from */
/*                    the poles of G1sr (see METHOD); */
/*             = 15:  the solution of the Sylvester equation failed */
/*                    because the zeros of W (if JOBW = 'I') or of */
/*                    conj(W) (if JOBW = 'R') are not distinct from */
/*                    the poles of G1sr (see METHOD); */
/*             = 16:  the solution of the generalized Sylvester system */
/*                    failed because the zeros of V (if JOBV = 'I') or */
/*                    of conj(V) (if JOBV = 'R') are not distinct from */
/*                    the poles of G1sr (see METHOD); */
/*             = 17:  the solution of the generalized Sylvester system */
/*                    failed because the zeros of W (if JOBW = 'I') or */
/*                    of conj(W) (if JOBW = 'R') are not distinct from */
/*                    the poles of G1sr (see METHOD); */
/*             = 18:  op(V) is not antistable; */
/*             = 19:  op(W) is not antistable; */
/*             = 20:  V is not invertible; */
/*             = 21:  W is not invertible. */

/*     METHOD */

/*     Let G be the transfer-function matrix of the original */
/*     linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t) + Du(t),                          (1) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system. The subroutine AB09JD determines */
/*     the matrices of a reduced order system */

/*          d[z(t)] = Ar*z(t) + Br*u(t) */
/*          yr(t)   = Cr*z(t) + Dr*u(t),                      (2) */

/*     such that the corresponding transfer-function matrix Gr minimizes */
/*     the Hankel-norm of the frequency-weighted error */

/*             op(V)*(G-Gr)*op(W).                            (3) */

/*     For minimizing (3) with op(V) = V and op(W) = W, V and W are */
/*     assumed to have poles distinct from those of G, while with */
/*     op(V) = conj(V) and op(W) = conj(W), conj(V) and conj(W) are */
/*     assumed to have poles distinct from those of G. For minimizing (3) */
/*     with op(V) = inv(V) and op(W) = inv(W), V and W are assumed to */
/*     have zeros distinct from the poles of G, while with */
/*     op(V) = conj(inv(V)) and op(W) = conj(inv(W)), conj(V) and conj(W) */
/*     are assumed to have zeros distinct from the poles of G. */

/*     Note: conj(G) = G'(-s) for a continuous-time system and */
/*           conj(G) = G'(1/z) for a discrete-time system. */

/*     The following procedure is used to reduce G (see [1]): */

/*     1) Decompose additively G as */

/*          G = G1 + G2, */

/*        such that G1 = (A1,B1,C1,D) has only ALPHA-stable poles and */
/*        G2 = (A2,B2,C2,0) has only ALPHA-unstable poles. */

/*     2) Compute G1s, the projection of op(V)*G1*op(W) containing the */
/*        poles of G1, using explicit formulas [4] or the inverse-free */
/*        descriptor system formulas of [5]. */

/*     3) Determine G1sr, the optimal Hankel-norm approximation of G1s, */
/*        of order r. */

/*     4) Compute G1r, the projection of inv(op(V))*G1sr*inv(op(W)) */
/*        containing the poles of G1sr, using explicit formulas [4] */
/*        or the inverse-free descriptor system formulas of [5]. */

/*     5) Assemble the reduced model Gr as */

/*           Gr = G1r + G2. */

/*     To reduce the weighted ALPHA-stable part G1s at step 3, the */
/*     optimal Hankel-norm approximation method of [2], based on the */
/*     square-root balancing projection formulas of [3], is employed. */

/*     The optimal weighted approximation error satisfies */

/*          HNORM[op(V)*(G-Gr)*op(W)] >= S(r+1), */

/*     where S(r+1) is the (r+1)-th Hankel singular value of G1s, the */
/*     transfer-function matrix computed at step 2 of the above */
/*     procedure, and HNORM(.) denotes the Hankel-norm. */

/*     REFERENCES */

/*     [1] Latham, G.A. and Anderson, B.D.O. */
/*         Frequency-weighted optimal Hankel-norm approximation of stable */
/*         transfer functions. */
/*         Systems & Control Letters, Vol. 5, pp. 229-236, 1985. */

/*     [2] Glover, K. */
/*         All optimal Hankel norm approximation of linear */
/*         multivariable systems and their L-infinity error bounds. */
/*         Int. J. Control, Vol. 36, pp. 1145-1193, 1984. */

/*     [3] Tombs, M.S. and Postlethwaite, I. */
/*         Truncated balanced realization of stable, non-minimal */
/*         state-space systems. */
/*         Int. J. Control, Vol. 46, pp. 1319-1330, 1987. */

/*     [4] Varga, A. */
/*         Explicit formulas for an efficient implementation */
/*         of the frequency-weighting model reduction approach. */
/*         Proc. 1993 European Control Conference, Groningen, NL, */
/*         pp. 693-696, 1993. */

/*     [5] Varga, A. */
/*         Efficient and numerically reliable implementation of the */
/*         frequency-weighted Hankel-norm approximation model reduction */
/*         approach. */
/*         Proc. 2001 ECC, Porto, Portugal, 2001. */

/*     NUMERICAL ASPECTS */

/*     The implemented methods rely on an accuracy enhancing square-root */
/*     technique. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, March 2001. */
/*     D. Sima, University of Bucharest, April 2001. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2001. */

/*     REVISIONS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2001. */
/*     V. Sima, Research Institute for Informatics, Bucharest, June 2001, */
/*     March 2005. */

/*     KEYWORDS */

/*     Frequency weighting, model reduction, multivariable system, */
/*     state-space model, state-space representation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 604 "AB09JD.f"
    /* Parameter adjustments */
#line 604 "AB09JD.f"
    a_dim1 = *lda;
#line 604 "AB09JD.f"
    a_offset = 1 + a_dim1;
#line 604 "AB09JD.f"
    a -= a_offset;
#line 604 "AB09JD.f"
    b_dim1 = *ldb;
#line 604 "AB09JD.f"
    b_offset = 1 + b_dim1;
#line 604 "AB09JD.f"
    b -= b_offset;
#line 604 "AB09JD.f"
    c_dim1 = *ldc;
#line 604 "AB09JD.f"
    c_offset = 1 + c_dim1;
#line 604 "AB09JD.f"
    c__ -= c_offset;
#line 604 "AB09JD.f"
    d_dim1 = *ldd;
#line 604 "AB09JD.f"
    d_offset = 1 + d_dim1;
#line 604 "AB09JD.f"
    d__ -= d_offset;
#line 604 "AB09JD.f"
    av_dim1 = *ldav;
#line 604 "AB09JD.f"
    av_offset = 1 + av_dim1;
#line 604 "AB09JD.f"
    av -= av_offset;
#line 604 "AB09JD.f"
    bv_dim1 = *ldbv;
#line 604 "AB09JD.f"
    bv_offset = 1 + bv_dim1;
#line 604 "AB09JD.f"
    bv -= bv_offset;
#line 604 "AB09JD.f"
    cv_dim1 = *ldcv;
#line 604 "AB09JD.f"
    cv_offset = 1 + cv_dim1;
#line 604 "AB09JD.f"
    cv -= cv_offset;
#line 604 "AB09JD.f"
    dv_dim1 = *lddv;
#line 604 "AB09JD.f"
    dv_offset = 1 + dv_dim1;
#line 604 "AB09JD.f"
    dv -= dv_offset;
#line 604 "AB09JD.f"
    aw_dim1 = *ldaw;
#line 604 "AB09JD.f"
    aw_offset = 1 + aw_dim1;
#line 604 "AB09JD.f"
    aw -= aw_offset;
#line 604 "AB09JD.f"
    bw_dim1 = *ldbw;
#line 604 "AB09JD.f"
    bw_offset = 1 + bw_dim1;
#line 604 "AB09JD.f"
    bw -= bw_offset;
#line 604 "AB09JD.f"
    cw_dim1 = *ldcw;
#line 604 "AB09JD.f"
    cw_offset = 1 + cw_dim1;
#line 604 "AB09JD.f"
    cw -= cw_offset;
#line 604 "AB09JD.f"
    dw_dim1 = *lddw;
#line 604 "AB09JD.f"
    dw_offset = 1 + dw_dim1;
#line 604 "AB09JD.f"
    dw -= dw_offset;
#line 604 "AB09JD.f"
    --hsv;
#line 604 "AB09JD.f"
    --iwork;
#line 604 "AB09JD.f"
    --dwork;
#line 604 "AB09JD.f"

#line 604 "AB09JD.f"
    /* Function Body */
#line 604 "AB09JD.f"
    *info = 0;
#line 605 "AB09JD.f"
    *iwarn = 0;
#line 606 "AB09JD.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 607 "AB09JD.f"
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);
#line 608 "AB09JD.f"
    lefti = lsame_(jobv, "I", (ftnlen)1, (ftnlen)1) || lsame_(jobv, "R", (
	    ftnlen)1, (ftnlen)1);
#line 609 "AB09JD.f"
    leftw = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1) || lsame_(jobv, "C", (
	    ftnlen)1, (ftnlen)1) || lefti;
#line 610 "AB09JD.f"
    conjv = lsame_(jobv, "C", (ftnlen)1, (ftnlen)1) || lsame_(jobv, "R", (
	    ftnlen)1, (ftnlen)1);
#line 611 "AB09JD.f"
    righti = lsame_(jobw, "I", (ftnlen)1, (ftnlen)1) || lsame_(jobw, "R", (
	    ftnlen)1, (ftnlen)1);
#line 612 "AB09JD.f"
    rightw = lsame_(jobw, "W", (ftnlen)1, (ftnlen)1) || lsame_(jobw, "C", (
	    ftnlen)1, (ftnlen)1) || righti;
#line 613 "AB09JD.f"
    conjw = lsame_(jobw, "C", (ftnlen)1, (ftnlen)1) || lsame_(jobw, "R", (
	    ftnlen)1, (ftnlen)1);
#line 614 "AB09JD.f"
    frwght = leftw || rightw;
#line 615 "AB09JD.f"
    invfr = lsame_(jobinv, "N", (ftnlen)1, (ftnlen)1);
#line 616 "AB09JD.f"
    autom = lsame_(jobinv, "A", (ftnlen)1, (ftnlen)1);

#line 618 "AB09JD.f"
    lw = 1;
#line 619 "AB09JD.f"
    if (leftw) {
#line 620 "AB09JD.f"
	nvp = *nv + *p;
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
#line 621 "AB09JD.f"
	i__5 = nvp * 11 + 16, i__6 = *p * nvp;
/* Computing MAX */
#line 621 "AB09JD.f"
	i__7 = nvp * *n + *n * *n, i__8 = *p * *n, i__7 = max(i__7,i__8), 
		i__8 = *p * *m;
#line 621 "AB09JD.f"
	i__3 = (nvp << 1) * nvp + max(i__5,i__6), i__4 = nvp * *n + max(i__7,
		i__8);
#line 621 "AB09JD.f"
	i__1 = lw, i__2 = (nvp << 1) * (nvp + *p) + *p * *p + max(i__3,i__4);
#line 621 "AB09JD.f"
	lw = max(i__1,i__2);
#line 624 "AB09JD.f"
    }
#line 625 "AB09JD.f"
    if (rightw) {
#line 626 "AB09JD.f"
	nwm = *nw + *m;
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
#line 627 "AB09JD.f"
	i__5 = nwm * 11 + 16, i__6 = *m * nwm;
/* Computing MAX */
#line 627 "AB09JD.f"
	i__7 = nwm * *n + *n * *n, i__8 = *m * *n, i__7 = max(i__7,i__8), 
		i__8 = *p * *m;
#line 627 "AB09JD.f"
	i__3 = (nwm << 1) * nwm + max(i__5,i__6), i__4 = nwm * *n + max(i__7,
		i__8);
#line 627 "AB09JD.f"
	i__1 = lw, i__2 = (nwm << 1) * (nwm + *m) + *m * *m + max(i__3,i__4);
#line 627 "AB09JD.f"
	lw = max(i__1,i__2);
#line 630 "AB09JD.f"
    }
/* Computing MAX */
/* Computing MAX */
#line 631 "AB09JD.f"
    i__3 = max(*n,*m);
#line 631 "AB09JD.f"
    i__1 = lw, i__2 = *n * ((*n << 1) + max(i__3,*p) + 5) + *n * (*n + 1) / 2;
#line 631 "AB09JD.f"
    lw = max(i__1,i__2);
/* Computing MAX */
/* Computing MAX */
#line 632 "AB09JD.f"
    i__3 = *m * 3 + 1, i__4 = min(*n,*m) + *p;
#line 632 "AB09JD.f"
    i__1 = lw, i__2 = *n * (*m + *p + 2) + (*m << 1) * *p + min(*n,*m) + max(
	    i__3,i__4);
#line 632 "AB09JD.f"
    lw = max(i__1,i__2);

/*     Check the input scalar arguments. */

#line 637 "AB09JD.f"
    if (! (lsame_(jobv, "N", (ftnlen)1, (ftnlen)1) || leftw)) {
#line 638 "AB09JD.f"
	*info = -1;
#line 639 "AB09JD.f"
    } else if (! (lsame_(jobw, "N", (ftnlen)1, (ftnlen)1) || rightw)) {
#line 640 "AB09JD.f"
	*info = -2;
#line 641 "AB09JD.f"
    } else if (! (invfr || autom || lsame_(jobinv, "I", (ftnlen)1, (ftnlen)1))
	    ) {
#line 643 "AB09JD.f"
	*info = -3;
#line 644 "AB09JD.f"
    } else if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 645 "AB09JD.f"
	*info = -4;
#line 646 "AB09JD.f"
    } else if (! (lsame_(equil, "S", (ftnlen)1, (ftnlen)1) || lsame_(equil, 
	    "N", (ftnlen)1, (ftnlen)1))) {
#line 648 "AB09JD.f"
	*info = -5;
#line 649 "AB09JD.f"
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
#line 650 "AB09JD.f"
	*info = -6;
#line 651 "AB09JD.f"
    } else if (*n < 0) {
#line 652 "AB09JD.f"
	*info = -7;
#line 653 "AB09JD.f"
    } else if (*nv < 0) {
#line 654 "AB09JD.f"
	*info = -8;
#line 655 "AB09JD.f"
    } else if (*nw < 0) {
#line 656 "AB09JD.f"
	*info = -9;
#line 657 "AB09JD.f"
    } else if (*m < 0) {
#line 658 "AB09JD.f"
	*info = -10;
#line 659 "AB09JD.f"
    } else if (*p < 0) {
#line 660 "AB09JD.f"
	*info = -11;
#line 661 "AB09JD.f"
    } else if (fixord && (*nr < 0 || *nr > *n)) {
#line 662 "AB09JD.f"
	*info = -12;
#line 663 "AB09JD.f"
    } else if (discr && (*alpha < 0. || *alpha > 1.) || ! discr && *alpha > 
	    0.) {
#line 665 "AB09JD.f"
	*info = -13;
#line 666 "AB09JD.f"
    } else if (*lda < max(1,*n)) {
#line 667 "AB09JD.f"
	*info = -15;
#line 668 "AB09JD.f"
    } else if (*ldb < max(1,*n)) {
#line 669 "AB09JD.f"
	*info = -17;
#line 670 "AB09JD.f"
    } else if (*ldc < max(1,*p)) {
#line 671 "AB09JD.f"
	*info = -19;
#line 672 "AB09JD.f"
    } else if (*ldd < max(1,*p)) {
#line 673 "AB09JD.f"
	*info = -21;
#line 674 "AB09JD.f"
    } else if (*ldav < 1 || leftw && *ldav < *nv) {
#line 675 "AB09JD.f"
	*info = -23;
#line 676 "AB09JD.f"
    } else if (*ldbv < 1 || leftw && *ldbv < *nv) {
#line 677 "AB09JD.f"
	*info = -25;
#line 678 "AB09JD.f"
    } else if (*ldcv < 1 || leftw && *ldcv < *p) {
#line 679 "AB09JD.f"
	*info = -27;
#line 680 "AB09JD.f"
    } else if (*lddv < 1 || leftw && *lddv < *p) {
#line 681 "AB09JD.f"
	*info = -29;
#line 682 "AB09JD.f"
    } else if (*ldaw < 1 || rightw && *ldaw < *nw) {
#line 683 "AB09JD.f"
	*info = -31;
#line 684 "AB09JD.f"
    } else if (*ldbw < 1 || rightw && *ldbw < *nw) {
#line 685 "AB09JD.f"
	*info = -33;
#line 686 "AB09JD.f"
    } else if (*ldcw < 1 || rightw && *ldcw < *m) {
#line 687 "AB09JD.f"
	*info = -35;
#line 688 "AB09JD.f"
    } else if (*lddw < 1 || rightw && *lddw < *m) {
#line 689 "AB09JD.f"
	*info = -37;
#line 690 "AB09JD.f"
    } else if (*tol1 >= 1.) {
#line 691 "AB09JD.f"
	*info = -40;
#line 692 "AB09JD.f"
    } else if (*tol2 > 0. && ! fixord && *tol2 > *tol1 || *tol2 >= 1.) {
#line 694 "AB09JD.f"
	*info = -41;
#line 695 "AB09JD.f"
    } else if (*ldwork < lw) {
#line 696 "AB09JD.f"
	*info = -44;
#line 697 "AB09JD.f"
    }

#line 699 "AB09JD.f"
    if (*info != 0) {

/*        Error return. */

#line 703 "AB09JD.f"
	i__1 = -(*info);
#line 703 "AB09JD.f"
	xerbla_("AB09JD", &i__1, (ftnlen)6);
#line 704 "AB09JD.f"
	return 0;
#line 705 "AB09JD.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 709 "AB09JD.f"
    i__1 = min(*n,*m);
#line 709 "AB09JD.f"
    if (min(i__1,*p) == 0) {
#line 710 "AB09JD.f"
	*nr = 0;
#line 711 "AB09JD.f"
	*ns = 0;
#line 712 "AB09JD.f"
	dwork[1] = 1.;
#line 713 "AB09JD.f"
	return 0;
#line 714 "AB09JD.f"
    }

#line 716 "AB09JD.f"
    if (lsame_(equil, "S", (ftnlen)1, (ftnlen)1)) {

/*        Scale simultaneously the matrices A, B and C: */
/*        A <- inv(D)*A*D,  B <- inv(D)*B  and  C <- C*D,  where D is a */
/*        diagonal matrix. */
/*        Workspace: N. */

#line 723 "AB09JD.f"
	maxred = 100.;
#line 724 "AB09JD.f"
	tb01id_("All", n, m, p, &maxred, &a[a_offset], lda, &b[b_offset], ldb,
		 &c__[c_offset], ldc, &dwork[1], info, (ftnlen)3);
#line 726 "AB09JD.f"
    }

/*     Correct the value of ALPHA to ensure stability. */

#line 730 "AB09JD.f"
    alpwrk = *alpha;
#line 731 "AB09JD.f"
    sqreps = sqrt(dlamch_("E", (ftnlen)1));
#line 732 "AB09JD.f"
    if (discr) {
#line 733 "AB09JD.f"
	if (*alpha == 1.) {
#line 733 "AB09JD.f"
	    alpwrk = 1. - sqreps;
#line 733 "AB09JD.f"
	}
#line 734 "AB09JD.f"
    } else {
#line 735 "AB09JD.f"
	if (*alpha == 0.) {
#line 735 "AB09JD.f"
	    alpwrk = -sqreps;
#line 735 "AB09JD.f"
	}
#line 736 "AB09JD.f"
    }

/*     Allocate working storage. */

#line 740 "AB09JD.f"
    ku = 1;
#line 741 "AB09JD.f"
    kl = ku + *n * *n;
#line 742 "AB09JD.f"
    ki = kl + *n;
#line 743 "AB09JD.f"
    kw = ki + *n;

/*     Compute an additive decomposition G = G1 + G2, where G1 */
/*     is the ALPHA-stable projection of G. */

/*     Reduce A to a block-diagonal real Schur form, with the NU-th order */
/*     ALPHA-unstable part in the leading diagonal position, using a */
/*     non-orthogonal similarity transformation A <- inv(T)*A*T and */
/*     apply the transformation to B and C: B <- inv(T)*B and C <- C*T. */

/*     Workspace needed:      N*(N+2); */
/*     Additional workspace:  need   3*N; */
/*                            prefer larger. */

#line 757 "AB09JD.f"
    i__1 = *ldwork - kw + 1;
#line 757 "AB09JD.f"
    tb01kd_(dico, "Unstable", "General", n, m, p, &alpwrk, &a[a_offset], lda, 
	    &b[b_offset], ldb, &c__[c_offset], ldc, &nu, &dwork[ku], n, &
	    dwork[kl], &dwork[ki], &dwork[kw], &i__1, &ierr, (ftnlen)1, (
	    ftnlen)8, (ftnlen)7);

#line 761 "AB09JD.f"
    if (ierr != 0) {
#line 762 "AB09JD.f"
	if (ierr != 3) {
#line 763 "AB09JD.f"
	    *info = 1;
#line 764 "AB09JD.f"
	} else {
#line 765 "AB09JD.f"
	    *info = 2;
#line 766 "AB09JD.f"
	}
#line 767 "AB09JD.f"
	return 0;
#line 768 "AB09JD.f"
    }

#line 770 "AB09JD.f"
    wrkopt = dwork[kw] + (doublereal) (kw - 1);
#line 771 "AB09JD.f"
    iwarnl = 0;

#line 773 "AB09JD.f"
    *ns = *n - nu;
#line 774 "AB09JD.f"
    if (fixord) {
/* Computing MAX */
#line 775 "AB09JD.f"
	i__1 = 0, i__2 = *nr - nu;
#line 775 "AB09JD.f"
	nra = max(i__1,i__2);
#line 776 "AB09JD.f"
	if (*nr < nu) {
#line 776 "AB09JD.f"
	    iwarnl = 2;
#line 776 "AB09JD.f"
	}
#line 778 "AB09JD.f"
    } else {
#line 779 "AB09JD.f"
	nra = 0;
#line 780 "AB09JD.f"
    }

/*     Finish if only unstable part is present. */

#line 784 "AB09JD.f"
    if (*ns == 0) {
#line 785 "AB09JD.f"
	*nr = nu;
#line 786 "AB09JD.f"
	dwork[1] = wrkopt;
#line 787 "AB09JD.f"
	return 0;
#line 788 "AB09JD.f"
    }

#line 790 "AB09JD.f"
    nu1 = nu + 1;
#line 791 "AB09JD.f"
    if (conjv) {
#line 792 "AB09JD.f"
	*(unsigned char *)jobvl = 'C';
#line 793 "AB09JD.f"
    } else {
#line 794 "AB09JD.f"
	*(unsigned char *)jobvl = 'V';
#line 795 "AB09JD.f"
    }
#line 796 "AB09JD.f"
    if (conjw) {
#line 797 "AB09JD.f"
	*(unsigned char *)jobwl = 'C';
#line 798 "AB09JD.f"
    } else {
#line 799 "AB09JD.f"
	*(unsigned char *)jobwl = 'W';
#line 800 "AB09JD.f"
    }
#line 801 "AB09JD.f"
    if (leftw) {

/*        Check if V is invertible. */
/*        Real workspace:    need   (NV+P)**2 + MAX( P + MAX(3*P,NV), */
/*                                  MIN(P+1,NV) + MAX(3*(P+1),NV+P) ); */
/*                           prefer larger. */
/*        Integer workspace: need   2*NV+P+2. */

#line 809 "AB09JD.f"
	tol = 0.;
#line 810 "AB09JD.f"
	ab08md_("S", nv, p, p, &av[av_offset], ldav, &bv[bv_offset], ldbv, &
		cv[cv_offset], ldcv, &dv[dv_offset], lddv, &rank, &tol, &
		iwork[1], &dwork[1], ldwork, &ierr, (ftnlen)1);
#line 813 "AB09JD.f"
	if (rank != *p) {
#line 814 "AB09JD.f"
	    *info = 20;
#line 815 "AB09JD.f"
	    return 0;
#line 816 "AB09JD.f"
	}
#line 817 "AB09JD.f"
	wrkopt = max(wrkopt,dwork[1]);

#line 819 "AB09JD.f"
	if (lefti) {
#line 820 "AB09JD.f"
	    if (invfr) {
#line 821 "AB09JD.f"
		ierr = 1;
#line 822 "AB09JD.f"
	    } else {

/*              Allocate storage for a standard inverse of V. */
/*              Workspace: need  NV*(NV+2*P) + P*P. */

#line 827 "AB09JD.f"
		kav = 1;
#line 828 "AB09JD.f"
		kbv = kav + *nv * *nv;
#line 829 "AB09JD.f"
		kcv = kbv + *nv * *p;
#line 830 "AB09JD.f"
		kdv = kcv + *p * *nv;
#line 831 "AB09JD.f"
		kw = kdv + *p * *p;

#line 833 "AB09JD.f"
		ldabv = max(*nv,1);
#line 834 "AB09JD.f"
		ldcdv = *p;
#line 835 "AB09JD.f"
		dlacpy_("Full", nv, nv, &av[av_offset], ldav, &dwork[kav], &
			ldabv, (ftnlen)4);
#line 837 "AB09JD.f"
		dlacpy_("Full", nv, p, &bv[bv_offset], ldbv, &dwork[kbv], &
			ldabv, (ftnlen)4);
#line 839 "AB09JD.f"
		dlacpy_("Full", p, nv, &cv[cv_offset], ldcv, &dwork[kcv], &
			ldcdv, (ftnlen)4);
#line 841 "AB09JD.f"
		dlacpy_("Full", p, p, &dv[dv_offset], lddv, &dwork[kdv], &
			ldcdv, (ftnlen)4);

/*              Compute the standard inverse of V. */
/*              Additional real workspace:   need   MAX(1,4*P); */
/*                                           prefer larger. */
/*              Integer workspace:           need   2*P. */

#line 849 "AB09JD.f"
		i__1 = *ldwork - kw + 1;
#line 849 "AB09JD.f"
		ab07nd_(nv, p, &dwork[kav], &ldabv, &dwork[kbv], &ldabv, &
			dwork[kcv], &ldcdv, &dwork[kdv], &ldcdv, &rcond, &
			iwork[1], &dwork[kw], &i__1, &ierr);
/* Computing MAX */
#line 852 "AB09JD.f"
		d__1 = wrkopt, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 852 "AB09JD.f"
		wrkopt = max(d__1,d__2);

/*              Check if inversion is accurate. */

#line 856 "AB09JD.f"
		if (autom) {
#line 857 "AB09JD.f"
		    if (ierr == 0 && rcond <= 1e-4) {
#line 857 "AB09JD.f"
			ierr = 1;
#line 857 "AB09JD.f"
		    }
#line 858 "AB09JD.f"
		} else {
#line 859 "AB09JD.f"
		    if (ierr == 0 && rcond <= sqreps) {
#line 859 "AB09JD.f"
			ierr = 1;
#line 859 "AB09JD.f"
		    }
#line 860 "AB09JD.f"
		}
#line 861 "AB09JD.f"
		if (ierr != 0 && *nv == 0) {
#line 862 "AB09JD.f"
		    *info = 20;
#line 863 "AB09JD.f"
		    return 0;
#line 864 "AB09JD.f"
		}
#line 865 "AB09JD.f"
	    }

#line 867 "AB09JD.f"
	    if (ierr != 0) {

/*              Allocate storage for a descriptor inverse of V. */

#line 871 "AB09JD.f"
		kav = 1;
#line 872 "AB09JD.f"
		kev = kav + nvp * nvp;
#line 873 "AB09JD.f"
		kbv = kev + nvp * nvp;
#line 874 "AB09JD.f"
		kcv = kbv + nvp * *p;
#line 875 "AB09JD.f"
		kdv = kcv + *p * nvp;
#line 876 "AB09JD.f"
		kw = kdv + *p * *p;

#line 878 "AB09JD.f"
		ldabv = max(nvp,1);
#line 879 "AB09JD.f"
		ldcdv = *p;

/*              DV is singular or ill-conditioned. */
/*              Form a descriptor inverse of V. */
/*              Workspace: need  2*(NV+P)*(NV+2*P) + P*P. */

#line 885 "AB09JD.f"
		ag07bd_("I", nv, p, &av[av_offset], ldav, temp, &c__1, &bv[
			bv_offset], ldbv, &cv[cv_offset], ldcv, &dv[dv_offset]
			, lddv, &dwork[kav], &ldabv, &dwork[kev], &ldabv, &
			dwork[kbv], &ldabv, &dwork[kcv], &ldcdv, &dwork[kdv], 
			&ldcdv, &ierr, (ftnlen)1);

/*              Compute the projection containing the poles of weighted */
/*              reduced ALPHA-stable part using descriptor inverse of V */
/*              of order NVP = NV + P. */
/*              Additional real workspace: need */
/*                 MAX( 2*NVP*NVP + MAX( 11*NVP+16, P*NVP ), */
/*                      NVP*N + MAX( NVP*N+N*N, P*N, P*M ) ); */
/*                 prefer larger. */
/*              Integer workspace: need NVP+N+6. */

#line 899 "AB09JD.f"
		i__1 = *ldwork - kw + 1;
#line 899 "AB09JD.f"
		ab09jv_(jobvl, dico, "G", "C", ns, m, p, &nvp, p, &a[nu1 + 
			nu1 * a_dim1], lda, &b[nu1 + b_dim1], ldb, &c__[nu1 * 
			c_dim1 + 1], ldc, &d__[d_offset], ldd, &dwork[kav], &
			ldabv, &dwork[kev], &ldabv, &dwork[kbv], &ldabv, &
			dwork[kcv], &ldcdv, &dwork[kdv], &ldcdv, &iwork[1], &
			dwork[kw], &i__1, &ierr, (ftnlen)1, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 906 "AB09JD.f"
		if (ierr != 0) {
#line 907 "AB09JD.f"
		    if (ierr == 1) {
#line 908 "AB09JD.f"
			*info = 5;
#line 909 "AB09JD.f"
		    } else if (ierr == 2) {
#line 910 "AB09JD.f"
			*info = 16;
#line 911 "AB09JD.f"
		    } else if (ierr == 4) {
#line 912 "AB09JD.f"
			*info = 18;
#line 913 "AB09JD.f"
		    }
#line 914 "AB09JD.f"
		    return 0;
#line 915 "AB09JD.f"
		}
#line 916 "AB09JD.f"
	    } else {

/*              Compute the projection containing the poles of weighted */
/*              reduced ALPHA-stable part using explicit inverse of V. */
/*              Additional real workspace: need */
/*                 MAX( NV*(NV+5), NV*N + MAX( a, P*N, P*M ) ) */
/*                      a = 0,    if DICO = 'C' or  JOBVL = 'V', */
/*                      a = 2*NV, if DICO = 'D' and JOBVL = 'C'; */
/*                 prefer larger. */

#line 926 "AB09JD.f"
		i__1 = *ldwork - kw + 1;
#line 926 "AB09JD.f"
		ab09jv_(jobvl, dico, "I", "C", ns, m, p, nv, p, &a[nu1 + nu1 *
			 a_dim1], lda, &b[nu1 + b_dim1], ldb, &c__[nu1 * 
			c_dim1 + 1], ldc, &d__[d_offset], ldd, &dwork[kav], &
			ldabv, temp, &c__1, &dwork[kbv], &ldabv, &dwork[kcv], 
			&ldcdv, &dwork[kdv], &ldcdv, &iwork[1], &dwork[kw], &
			i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
			1);
#line 932 "AB09JD.f"
		if (ierr != 0) {
#line 933 "AB09JD.f"
		    if (ierr == 1) {
#line 934 "AB09JD.f"
			*info = 10;
#line 935 "AB09JD.f"
		    } else if (ierr == 3) {
#line 936 "AB09JD.f"
			*info = 14;
#line 937 "AB09JD.f"
		    } else if (ierr == 4) {
#line 938 "AB09JD.f"
			*info = 18;
#line 939 "AB09JD.f"
		    }
#line 940 "AB09JD.f"
		    return 0;
#line 941 "AB09JD.f"
		}
#line 942 "AB09JD.f"
	    }

/* Computing MAX */
#line 944 "AB09JD.f"
	    d__1 = wrkopt, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 944 "AB09JD.f"
	    wrkopt = max(d__1,d__2);
#line 945 "AB09JD.f"
	} else {

/*           Compute the projection of V*G1 or conj(V)*G1 containing the */
/*           poles of G. */

/*           Workspace need: */
/*           real   MAX( 1, NV*(NV+5), NV*N + MAX( a, P*N, P*M ) ) */
/*                       a = 0,    if DICO = 'C' or  JOBVL = 'V', */
/*                       a = 2*NV, if DICO = 'D' and JOBVL = 'C'; */
/*           prefer larger. */

#line 956 "AB09JD.f"
	    ab09jv_(jobvl, dico, "I", "C", ns, m, p, nv, p, &a[nu1 + nu1 * 
		    a_dim1], lda, &b[nu1 + b_dim1], ldb, &c__[nu1 * c_dim1 + 
		    1], ldc, &d__[d_offset], ldd, &av[av_offset], ldav, temp, 
		    &c__1, &bv[bv_offset], ldbv, &cv[cv_offset], ldcv, &dv[
		    dv_offset], lddv, &iwork[1], &dwork[1], ldwork, &ierr, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 961 "AB09JD.f"
	    if (ierr != 0) {
#line 962 "AB09JD.f"
		if (ierr == 1) {
#line 963 "AB09JD.f"
		    *info = 3;
#line 964 "AB09JD.f"
		} else if (ierr == 3) {
#line 965 "AB09JD.f"
		    *info = 12;
#line 966 "AB09JD.f"
		} else if (ierr == 4) {
#line 967 "AB09JD.f"
		    *info = 18;
#line 968 "AB09JD.f"
		}
#line 969 "AB09JD.f"
		return 0;
#line 970 "AB09JD.f"
	    }

#line 972 "AB09JD.f"
	    wrkopt = max(wrkopt,dwork[1]);
#line 973 "AB09JD.f"
	}
#line 974 "AB09JD.f"
    }

#line 976 "AB09JD.f"
    if (rightw) {

/*        Check if W is invertible. */
/*        Real workspace:    need   (NW+M)**2 + MAX( M + MAX(3*M,NW), */
/*                                  MIN(M+1,NW) + MAX(3*(M+1),NW+M) ); */
/*                           prefer larger. */
/*        Integer workspace: need   2*NW+M+2. */

#line 984 "AB09JD.f"
	tol = 0.;
#line 985 "AB09JD.f"
	ab08md_("S", nw, m, m, &aw[aw_offset], ldaw, &bw[bw_offset], ldbw, &
		cw[cw_offset], ldcw, &dw[dw_offset], lddw, &rank, &tol, &
		iwork[1], &dwork[1], ldwork, &ierr, (ftnlen)1);
#line 988 "AB09JD.f"
	if (rank != *m) {
#line 989 "AB09JD.f"
	    *info = 21;
#line 990 "AB09JD.f"
	    return 0;
#line 991 "AB09JD.f"
	}
#line 992 "AB09JD.f"
	wrkopt = max(wrkopt,dwork[1]);

#line 994 "AB09JD.f"
	if (righti) {
#line 995 "AB09JD.f"
	    if (invfr) {
#line 996 "AB09JD.f"
		ierr = 1;
#line 997 "AB09JD.f"
	    } else {

/*              Allocate storage for a standard inverse of W. */
/*              Workspace: need  NW*(NW+2*M) + M*M. */

#line 1002 "AB09JD.f"
		kaw = 1;
#line 1003 "AB09JD.f"
		kbw = kaw + *nw * *nw;
#line 1004 "AB09JD.f"
		kcw = kbw + *nw * *m;
#line 1005 "AB09JD.f"
		kdw = kcw + *m * *nw;
#line 1006 "AB09JD.f"
		kw = kdw + *m * *m;

#line 1008 "AB09JD.f"
		ldabw = max(*nw,1);
#line 1009 "AB09JD.f"
		ldcdw = *m;
#line 1010 "AB09JD.f"
		dlacpy_("Full", nw, nw, &aw[aw_offset], ldaw, &dwork[kaw], &
			ldabw, (ftnlen)4);
#line 1012 "AB09JD.f"
		dlacpy_("Full", nw, m, &bw[bw_offset], ldbw, &dwork[kbw], &
			ldabw, (ftnlen)4);
#line 1014 "AB09JD.f"
		dlacpy_("Full", m, nw, &cw[cw_offset], ldcw, &dwork[kcw], &
			ldcdw, (ftnlen)4);
#line 1016 "AB09JD.f"
		dlacpy_("Full", m, m, &dw[dw_offset], lddw, &dwork[kdw], &
			ldcdw, (ftnlen)4);

/*              Compute the standard inverse of W. */
/*              Additional real workspace:   need   MAX(1,4*M); */
/*                                           prefer larger. */
/*              Integer workspace:           need   2*M. */

#line 1024 "AB09JD.f"
		i__1 = *ldwork - kw + 1;
#line 1024 "AB09JD.f"
		ab07nd_(nw, m, &dwork[kaw], &ldabw, &dwork[kbw], &ldabw, &
			dwork[kcw], &ldcdw, &dwork[kdw], &ldcdw, &rcond, &
			iwork[1], &dwork[kw], &i__1, &ierr);
/* Computing MAX */
#line 1027 "AB09JD.f"
		d__1 = wrkopt, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 1027 "AB09JD.f"
		wrkopt = max(d__1,d__2);

/*              Check if inversion is accurate. */

#line 1031 "AB09JD.f"
		if (autom) {
#line 1032 "AB09JD.f"
		    if (ierr == 0 && rcond <= 1e-4) {
#line 1032 "AB09JD.f"
			ierr = 1;
#line 1032 "AB09JD.f"
		    }
#line 1033 "AB09JD.f"
		} else {
#line 1034 "AB09JD.f"
		    if (ierr == 0 && rcond <= sqreps) {
#line 1034 "AB09JD.f"
			ierr = 1;
#line 1034 "AB09JD.f"
		    }
#line 1035 "AB09JD.f"
		}
#line 1036 "AB09JD.f"
		if (ierr != 0 && *nw == 0) {
#line 1037 "AB09JD.f"
		    *info = 21;
#line 1038 "AB09JD.f"
		    return 0;
#line 1039 "AB09JD.f"
		}
#line 1040 "AB09JD.f"
	    }

#line 1042 "AB09JD.f"
	    if (ierr != 0) {

/*              Allocate storage for a descriptor inverse of W. */

#line 1046 "AB09JD.f"
		kaw = 1;
#line 1047 "AB09JD.f"
		kew = kaw + nwm * nwm;
#line 1048 "AB09JD.f"
		kbw = kew + nwm * nwm;
#line 1049 "AB09JD.f"
		kcw = kbw + nwm * *m;
#line 1050 "AB09JD.f"
		kdw = kcw + *m * nwm;
#line 1051 "AB09JD.f"
		kw = kdw + *m * *m;

#line 1053 "AB09JD.f"
		ldabw = max(nwm,1);
#line 1054 "AB09JD.f"
		ldcdw = *m;

/*              DW is singular or ill-conditioned. */
/*              Form the descriptor inverse of W. */
/*              Workspace: need  2*(NW+M)*(NW+2*M) + M*M. */

#line 1060 "AB09JD.f"
		ag07bd_("I", nw, m, &aw[aw_offset], ldaw, temp, &c__1, &bw[
			bw_offset], ldbw, &cw[cw_offset], ldcw, &dw[dw_offset]
			, lddw, &dwork[kaw], &ldabw, &dwork[kew], &ldabw, &
			dwork[kbw], &ldabw, &dwork[kcw], &ldcdw, &dwork[kdw], 
			&ldcdw, &ierr, (ftnlen)1);

/*              Compute the projection containing the poles of weighted */
/*              reduced ALPHA-stable part using descriptor inverse of W */
/*              of order NWM = NW + M. */
/*              Additional real workspace: need */
/*                 MAX( 2*NWM*NWM + MAX( 11*NWM+16, M*NWM ), */
/*                      NWM*N + MAX( NWM*N+N*N, M*N, P*M ) ); */
/*                 prefer larger. */
/*              Integer workspace: need NWM+N+6. */

#line 1074 "AB09JD.f"
		i__1 = *ldwork - kw + 1;
#line 1074 "AB09JD.f"
		ab09jw_(jobwl, dico, "G", "C", ns, m, p, &nwm, m, &a[nu1 + 
			nu1 * a_dim1], lda, &b[nu1 + b_dim1], ldb, &c__[nu1 * 
			c_dim1 + 1], ldc, &d__[d_offset], ldd, &dwork[kaw], &
			ldabw, &dwork[kew], &ldabw, &dwork[kbw], &ldabw, &
			dwork[kcw], &ldcdw, &dwork[kdw], &ldcdw, &iwork[1], &
			dwork[kw], &i__1, &ierr, (ftnlen)1, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 1080 "AB09JD.f"
		if (ierr != 0) {
#line 1081 "AB09JD.f"
		    if (ierr == 1) {
#line 1082 "AB09JD.f"
			*info = 6;
#line 1083 "AB09JD.f"
		    } else if (ierr == 2) {
#line 1084 "AB09JD.f"
			*info = 17;
#line 1085 "AB09JD.f"
		    } else if (ierr == 4) {
#line 1086 "AB09JD.f"
			*info = 19;
#line 1087 "AB09JD.f"
		    }
#line 1088 "AB09JD.f"
		    return 0;
#line 1089 "AB09JD.f"
		}
#line 1090 "AB09JD.f"
	    } else {

/*              Compute the projection containing the poles of weighted */
/*              reduced ALPHA-stable part using explicit inverse of W. */
/*              Additional real workspace: need */
/*                 MAX( NW*(NW+5), NW*N + MAX( a, M*N, P*M ) ) */
/*                      a = 0,    if DICO = 'C' or  JOBWL = 'W', */
/*                      a = 2*NW, if DICO = 'D' and JOBWL = 'C'; */
/*                 prefer larger. */

#line 1100 "AB09JD.f"
		i__1 = *ldwork - kw + 1;
#line 1100 "AB09JD.f"
		ab09jw_(jobwl, dico, "I", "C", ns, m, p, nw, m, &a[nu1 + nu1 *
			 a_dim1], lda, &b[nu1 + b_dim1], ldb, &c__[nu1 * 
			c_dim1 + 1], ldc, &d__[d_offset], ldd, &dwork[kaw], &
			ldabw, temp, &c__1, &dwork[kbw], &ldabw, &dwork[kcw], 
			&ldcdw, &dwork[kdw], &ldcdw, &iwork[1], &dwork[kw], &
			i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
			1);
#line 1106 "AB09JD.f"
		if (ierr != 0) {
#line 1107 "AB09JD.f"
		    if (ierr == 1) {
#line 1108 "AB09JD.f"
			*info = 11;
#line 1109 "AB09JD.f"
		    } else if (ierr == 3) {
#line 1110 "AB09JD.f"
			*info = 15;
#line 1111 "AB09JD.f"
		    } else if (ierr == 4) {
#line 1112 "AB09JD.f"
			*info = 19;
#line 1113 "AB09JD.f"
		    }
#line 1114 "AB09JD.f"
		    return 0;
#line 1115 "AB09JD.f"
		}
#line 1116 "AB09JD.f"
	    }

/* Computing MAX */
#line 1118 "AB09JD.f"
	    d__1 = wrkopt, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 1118 "AB09JD.f"
	    wrkopt = max(d__1,d__2);
#line 1119 "AB09JD.f"
	} else {

/*           Compute the projection G1s of V*G1*W or conj(V)*G1*conj(W) */
/*           containing the poles of G. */

/*           Workspace need: */
/*           real   MAX( 1, NW*(NW+5), NW*N + MAX( b, M*N, P*M ) ) */
/*                    b = 0,    if DICO = 'C' or  JOBWL = 'W', */
/*                    b = 2*NW, if DICO = 'D' and JOBWL = 'C'; */
/*           prefer larger. */

#line 1130 "AB09JD.f"
	    ab09jw_(jobwl, dico, "I", "C", ns, m, p, nw, m, &a[nu1 + nu1 * 
		    a_dim1], lda, &b[nu1 + b_dim1], ldb, &c__[nu1 * c_dim1 + 
		    1], ldc, &d__[d_offset], ldd, &aw[aw_offset], ldaw, temp, 
		    &c__1, &bw[bw_offset], ldbw, &cw[cw_offset], ldcw, &dw[
		    dw_offset], lddw, &iwork[1], &dwork[1], ldwork, &ierr, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1134 "AB09JD.f"
	    if (ierr != 0) {
#line 1135 "AB09JD.f"
		if (ierr == 1) {
#line 1136 "AB09JD.f"
		    *info = 4;
#line 1137 "AB09JD.f"
		} else if (ierr == 3) {
#line 1138 "AB09JD.f"
		    *info = 13;
#line 1139 "AB09JD.f"
		} else if (ierr == 4) {
#line 1140 "AB09JD.f"
		    *info = 19;
#line 1141 "AB09JD.f"
		}
#line 1142 "AB09JD.f"
		return 0;
#line 1143 "AB09JD.f"
	    }

#line 1145 "AB09JD.f"
	    wrkopt = max(wrkopt,dwork[1]);
#line 1146 "AB09JD.f"
	}
#line 1147 "AB09JD.f"
    }

/*     Determine a reduced order approximation G1sr of G1s using the */
/*     Hankel-norm approximation method. The resulting A(NU1:N,NU1:N) */
/*     is further in a real Schur form. */

/*     Workspace: need   MAX( LDW3, LDW4 ), */
/*                LDW3 = N*(2*N + MAX(N,M,P) + 5) + N*(N+1)/2, */
/*                LDW4 = N*(M+P+2) + 2*M*P + MIN(N,M) + */
/*                       MAX( 3*M+1, MIN(N,M)+P ); */
/*                prefer larger. */

#line 1159 "AB09JD.f"
    ab09cx_(dico, ordsel, ns, m, p, &nra, &a[nu1 + nu1 * a_dim1], lda, &b[nu1 
	    + b_dim1], ldb, &c__[nu1 * c_dim1 + 1], ldc, &d__[d_offset], ldd, 
	    &hsv[1], tol1, tol2, &iwork[1], &dwork[1], ldwork, iwarn, &ierr, (
	    ftnlen)1, (ftnlen)1);

#line 1163 "AB09JD.f"
    if (ierr != 0) {

/*        Set INFO = 7, 8 or 9. */

#line 1167 "AB09JD.f"
	*info = ierr + 5;
#line 1168 "AB09JD.f"
	return 0;
#line 1169 "AB09JD.f"
    }

#line 1171 "AB09JD.f"
    *iwarn = max(iwarnl,*iwarn);
#line 1172 "AB09JD.f"
    wrkopt = max(wrkopt,dwork[1]);

#line 1174 "AB09JD.f"
    if (leftw) {
#line 1175 "AB09JD.f"
	if (! lefti) {
#line 1176 "AB09JD.f"
	    if (invfr) {
#line 1177 "AB09JD.f"
		ierr = 1;
#line 1178 "AB09JD.f"
	    } else {

/*              Allocate storage for a standard inverse of V. */
/*              Workspace: need  NV*(NV+2*P) + P*P. */

#line 1183 "AB09JD.f"
		kav = 1;
#line 1184 "AB09JD.f"
		kbv = kav + *nv * *nv;
#line 1185 "AB09JD.f"
		kcv = kbv + *nv * *p;
#line 1186 "AB09JD.f"
		kdv = kcv + *p * *nv;
#line 1187 "AB09JD.f"
		kw = kdv + *p * *p;

#line 1189 "AB09JD.f"
		ldabv = max(*nv,1);
#line 1190 "AB09JD.f"
		ldcdv = *p;
#line 1191 "AB09JD.f"
		dlacpy_("Full", nv, nv, &av[av_offset], ldav, &dwork[kav], &
			ldabv, (ftnlen)4);
#line 1193 "AB09JD.f"
		dlacpy_("Full", nv, p, &bv[bv_offset], ldbv, &dwork[kbv], &
			ldabv, (ftnlen)4);
#line 1195 "AB09JD.f"
		dlacpy_("Full", p, nv, &cv[cv_offset], ldcv, &dwork[kcv], &
			ldcdv, (ftnlen)4);
#line 1197 "AB09JD.f"
		dlacpy_("Full", p, p, &dv[dv_offset], lddv, &dwork[kdv], &
			ldcdv, (ftnlen)4);

/*              Compute the standard inverse of V. */
/*              Additional real workspace:   need   MAX(1,4*P); */
/*                                           prefer larger. */
/*              Integer workspace:           need   2*P. */

#line 1205 "AB09JD.f"
		i__1 = *ldwork - kw + 1;
#line 1205 "AB09JD.f"
		ab07nd_(nv, p, &dwork[kav], &ldabv, &dwork[kbv], &ldabv, &
			dwork[kcv], &ldcdv, &dwork[kdv], &ldcdv, &rcond, &
			iwork[1], &dwork[kw], &i__1, &ierr);
/* Computing MAX */
#line 1208 "AB09JD.f"
		d__1 = wrkopt, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 1208 "AB09JD.f"
		wrkopt = max(d__1,d__2);

/*              Check if inversion is accurate. */

#line 1212 "AB09JD.f"
		if (autom) {
#line 1213 "AB09JD.f"
		    if (ierr == 0 && rcond <= 1e-4) {
#line 1213 "AB09JD.f"
			ierr = 1;
#line 1213 "AB09JD.f"
		    }
#line 1214 "AB09JD.f"
		} else {
#line 1215 "AB09JD.f"
		    if (ierr == 0 && rcond <= sqreps) {
#line 1215 "AB09JD.f"
			ierr = 1;
#line 1215 "AB09JD.f"
		    }
#line 1216 "AB09JD.f"
		}
#line 1217 "AB09JD.f"
		if (ierr != 0 && *nv == 0) {
#line 1218 "AB09JD.f"
		    *info = 20;
#line 1219 "AB09JD.f"
		    return 0;
#line 1220 "AB09JD.f"
		}
#line 1221 "AB09JD.f"
	    }

#line 1223 "AB09JD.f"
	    if (ierr != 0) {

/*              Allocate storage for a descriptor inverse of V. */

#line 1227 "AB09JD.f"
		kav = 1;
#line 1228 "AB09JD.f"
		kev = kav + nvp * nvp;
#line 1229 "AB09JD.f"
		kbv = kev + nvp * nvp;
#line 1230 "AB09JD.f"
		kcv = kbv + nvp * *p;
#line 1231 "AB09JD.f"
		kdv = kcv + *p * nvp;
#line 1232 "AB09JD.f"
		kw = kdv + *p * *p;

#line 1234 "AB09JD.f"
		ldabv = max(nvp,1);
#line 1235 "AB09JD.f"
		ldcdv = *p;

/*              DV is singular or ill-conditioned. */
/*              Form a descriptor inverse of V. */
/*              Workspace: need  2*(NV+P)*(NV+2*P) + P*P. */

#line 1241 "AB09JD.f"
		ag07bd_("I", nv, p, &av[av_offset], ldav, temp, &c__1, &bv[
			bv_offset], ldbv, &cv[cv_offset], ldcv, &dv[dv_offset]
			, lddv, &dwork[kav], &ldabv, &dwork[kev], &ldabv, &
			dwork[kbv], &ldabv, &dwork[kcv], &ldcdv, &dwork[kdv], 
			&ldcdv, &ierr, (ftnlen)1);

/*              Compute the projection containing the poles of weighted */
/*              reduced ALPHA-stable part using descriptor inverse of V */
/*              of order NVP = NV + P. */
/*              Additional real workspace: need */
/*                 MAX( 2*NVP*NVP + MAX( 11*NVP+16, P*NVP ), */
/*                      NVP*N + MAX( NVP*N+N*N, P*N, P*M ) ); */
/*                 prefer larger. */
/*              Integer workspace: need NVP+N+6. */

#line 1255 "AB09JD.f"
		i__1 = *ldwork - kw + 1;
#line 1255 "AB09JD.f"
		ab09jv_(jobvl, dico, "G", "N", &nra, m, p, &nvp, p, &a[nu1 + 
			nu1 * a_dim1], lda, &b[nu1 + b_dim1], ldb, &c__[nu1 * 
			c_dim1 + 1], ldc, &d__[d_offset], ldd, &dwork[kav], &
			ldabv, &dwork[kev], &ldabv, &dwork[kbv], &ldabv, &
			dwork[kcv], &ldcdv, &dwork[kdv], &ldcdv, &iwork[1], &
			dwork[kw], &i__1, &ierr, (ftnlen)1, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 1262 "AB09JD.f"
		if (ierr != 0) {
#line 1263 "AB09JD.f"
		    if (ierr == 1) {
#line 1264 "AB09JD.f"
			*info = 5;
#line 1265 "AB09JD.f"
		    } else if (ierr == 2) {
#line 1266 "AB09JD.f"
			*info = 16;
#line 1267 "AB09JD.f"
		    }
#line 1268 "AB09JD.f"
		    return 0;
#line 1269 "AB09JD.f"
		}
#line 1270 "AB09JD.f"
	    } else {

/*              Compute the projection containing the poles of weighted */
/*              reduced ALPHA-stable part using explicit inverse of V. */
/*              Additional real workspace: need */
/*                 MAX( NV*(NV+5), NV*N + MAX( a, P*N, P*M ) ) */
/*                      a = 0,    if DICO = 'C' or  JOBVL = 'V', */
/*                      a = 2*NV, if DICO = 'D' and JOBVL = 'C'; */
/*                 prefer larger. */

#line 1280 "AB09JD.f"
		i__1 = *ldwork - kw + 1;
#line 1280 "AB09JD.f"
		ab09jv_(jobvl, dico, "I", "N", &nra, m, p, nv, p, &a[nu1 + 
			nu1 * a_dim1], lda, &b[nu1 + b_dim1], ldb, &c__[nu1 * 
			c_dim1 + 1], ldc, &d__[d_offset], ldd, &dwork[kav], &
			ldabv, temp, &c__1, &dwork[kbv], &ldabv, &dwork[kcv], 
			&ldcdv, &dwork[kdv], &ldcdv, &iwork[1], &dwork[kw], &
			i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
			1);
#line 1286 "AB09JD.f"
		if (ierr != 0) {
#line 1287 "AB09JD.f"
		    if (ierr == 1) {
#line 1288 "AB09JD.f"
			*info = 10;
#line 1289 "AB09JD.f"
		    } else if (ierr == 3) {
#line 1290 "AB09JD.f"
			*info = 14;
#line 1291 "AB09JD.f"
		    }
#line 1292 "AB09JD.f"
		    return 0;
#line 1293 "AB09JD.f"
		}
#line 1294 "AB09JD.f"
	    }

/* Computing MAX */
#line 1296 "AB09JD.f"
	    d__1 = wrkopt, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 1296 "AB09JD.f"
	    wrkopt = max(d__1,d__2);
#line 1297 "AB09JD.f"
	} else {

/*           Compute the projection of V*G1sr or conj(V)*G1sr containing */
/*           the poles of G. */

/*           Workspace need: */
/*           real    MAX( 1, NV*(NV+5), NV*N + MAX( a, P*N, P*M ) ) */
/*                        a = 0,    if DICO = 'C' or  JOBVL = 'V', */
/*                        a = 2*NV, if DICO = 'D' and JOBVL = 'C'; */
/*           prefer larger. */

#line 1308 "AB09JD.f"
	    ab09jv_(jobvl, dico, "I", "N", &nra, m, p, nv, p, &a[nu1 + nu1 * 
		    a_dim1], lda, &b[nu1 + b_dim1], ldb, &c__[nu1 * c_dim1 + 
		    1], ldc, &d__[d_offset], ldd, &av[av_offset], ldav, temp, 
		    &c__1, &bv[bv_offset], ldbv, &cv[cv_offset], ldcv, &dv[
		    dv_offset], lddv, &iwork[1], &dwork[1], ldwork, &ierr, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1313 "AB09JD.f"
	    if (ierr != 0) {
#line 1314 "AB09JD.f"
		if (ierr == 1) {
#line 1315 "AB09JD.f"
		    *info = 3;
#line 1316 "AB09JD.f"
		} else if (ierr == 3) {
#line 1317 "AB09JD.f"
		    *info = 12;
#line 1318 "AB09JD.f"
		}
#line 1319 "AB09JD.f"
		return 0;
#line 1320 "AB09JD.f"
	    }

#line 1322 "AB09JD.f"
	    wrkopt = max(wrkopt,dwork[1]);
#line 1323 "AB09JD.f"
	}
#line 1324 "AB09JD.f"
    }

#line 1326 "AB09JD.f"
    if (rightw) {
#line 1327 "AB09JD.f"
	if (! righti) {
#line 1328 "AB09JD.f"
	    if (invfr) {
#line 1329 "AB09JD.f"
		ierr = 1;
#line 1330 "AB09JD.f"
	    } else {

/*              Allocate storage for a standard inverse of W. */
/*              Workspace: need  NW*(NW+2*M) + M*M. */

#line 1335 "AB09JD.f"
		kaw = 1;
#line 1336 "AB09JD.f"
		kbw = kaw + *nw * *nw;
#line 1337 "AB09JD.f"
		kcw = kbw + *nw * *m;
#line 1338 "AB09JD.f"
		kdw = kcw + *m * *nw;
#line 1339 "AB09JD.f"
		kw = kdw + *m * *m;

#line 1341 "AB09JD.f"
		ldabw = max(*nw,1);
#line 1342 "AB09JD.f"
		ldcdw = *m;
#line 1343 "AB09JD.f"
		dlacpy_("Full", nw, nw, &aw[aw_offset], ldaw, &dwork[kaw], &
			ldabw, (ftnlen)4);
#line 1345 "AB09JD.f"
		dlacpy_("Full", nw, m, &bw[bw_offset], ldbw, &dwork[kbw], &
			ldabw, (ftnlen)4);
#line 1347 "AB09JD.f"
		dlacpy_("Full", m, nw, &cw[cw_offset], ldcw, &dwork[kcw], &
			ldcdw, (ftnlen)4);
#line 1349 "AB09JD.f"
		dlacpy_("Full", m, m, &dw[dw_offset], lddw, &dwork[kdw], &
			ldcdw, (ftnlen)4);

/*              Compute the standard inverse of W. */
/*              Additional real workspace:   need   MAX(1,4*M); */
/*                                           prefer larger. */
/*              Integer workspace:           need   2*M. */

#line 1357 "AB09JD.f"
		i__1 = *ldwork - kw + 1;
#line 1357 "AB09JD.f"
		ab07nd_(nw, m, &dwork[kaw], &ldabw, &dwork[kbw], &ldabw, &
			dwork[kcw], &ldcdw, &dwork[kdw], &ldcdw, &rcond, &
			iwork[1], &dwork[kw], &i__1, &ierr);
/* Computing MAX */
#line 1360 "AB09JD.f"
		d__1 = wrkopt, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 1360 "AB09JD.f"
		wrkopt = max(d__1,d__2);

/*              Check if inversion is accurate. */

#line 1364 "AB09JD.f"
		if (autom) {
#line 1365 "AB09JD.f"
		    if (ierr == 0 && rcond <= 1e-4) {
#line 1365 "AB09JD.f"
			ierr = 1;
#line 1365 "AB09JD.f"
		    }
#line 1366 "AB09JD.f"
		} else {
#line 1367 "AB09JD.f"
		    if (ierr == 0 && rcond <= sqreps) {
#line 1367 "AB09JD.f"
			ierr = 1;
#line 1367 "AB09JD.f"
		    }
#line 1368 "AB09JD.f"
		}
#line 1369 "AB09JD.f"
		if (ierr != 0 && *nw == 0) {
#line 1370 "AB09JD.f"
		    *info = 21;
#line 1371 "AB09JD.f"
		    return 0;
#line 1372 "AB09JD.f"
		}
#line 1373 "AB09JD.f"
	    }

#line 1375 "AB09JD.f"
	    if (ierr != 0) {

/*              Allocate storage for a descriptor inverse of W. */

#line 1379 "AB09JD.f"
		kaw = 1;
#line 1380 "AB09JD.f"
		kew = kaw + nwm * nwm;
#line 1381 "AB09JD.f"
		kbw = kew + nwm * nwm;
#line 1382 "AB09JD.f"
		kcw = kbw + nwm * *m;
#line 1383 "AB09JD.f"
		kdw = kcw + *m * nwm;
#line 1384 "AB09JD.f"
		kw = kdw + *m * *m;

#line 1386 "AB09JD.f"
		ldabw = max(nwm,1);
#line 1387 "AB09JD.f"
		ldcdw = *m;

/*              DW is singular or ill-conditioned. */
/*              Form the descriptor inverse of W. */
/*              Workspace: need  2*(NW+M)*(NW+2*M) + M*M. */

#line 1393 "AB09JD.f"
		ag07bd_("I", nw, m, &aw[aw_offset], ldaw, temp, &c__1, &bw[
			bw_offset], ldbw, &cw[cw_offset], ldcw, &dw[dw_offset]
			, lddw, &dwork[kaw], &ldabw, &dwork[kew], &ldabw, &
			dwork[kbw], &ldabw, &dwork[kcw], &ldcdw, &dwork[kdw], 
			&ldcdw, &ierr, (ftnlen)1);

/*              Compute the projection containing the poles of weighted */
/*              reduced ALPHA-stable part using descriptor inverse of W */
/*              of order NWM = NW + M. */
/*              Additional real workspace: need */
/*                 MAX( 2*NWM*NWM + MAX( 11*NWM+16, M*NWM ), */
/*                      NWM*N + MAX( NWM*N+N*N, M*N, P*M ) ); */
/*                 prefer larger. */
/*              Integer workspace: need NWM+N+6. */

#line 1407 "AB09JD.f"
		i__1 = *ldwork - kw + 1;
#line 1407 "AB09JD.f"
		ab09jw_(jobwl, dico, "G", "N", &nra, m, p, &nwm, m, &a[nu1 + 
			nu1 * a_dim1], lda, &b[nu1 + b_dim1], ldb, &c__[nu1 * 
			c_dim1 + 1], ldc, &d__[d_offset], ldd, &dwork[kaw], &
			ldabw, &dwork[kew], &ldabw, &dwork[kbw], &ldabw, &
			dwork[kcw], &ldcdw, &dwork[kdw], &ldcdw, &iwork[1], &
			dwork[kw], &i__1, &ierr, (ftnlen)1, (ftnlen)1, (
			ftnlen)1, (ftnlen)1);
#line 1413 "AB09JD.f"
		if (ierr != 0) {
#line 1414 "AB09JD.f"
		    if (ierr == 1) {
#line 1415 "AB09JD.f"
			*info = 6;
#line 1416 "AB09JD.f"
		    } else if (ierr == 2) {
#line 1417 "AB09JD.f"
			*info = 17;
#line 1418 "AB09JD.f"
		    }
#line 1419 "AB09JD.f"
		    return 0;
#line 1420 "AB09JD.f"
		}
#line 1421 "AB09JD.f"
	    } else {

/*              Compute the projection containing the poles of weighted */
/*              reduced ALPHA-stable part using explicit inverse of W. */
/*              Additional real workspace: need */
/*                 MAX( NW*(NW+5), NW*N + MAX( a, M*N, P*M ) ) */
/*                      a = 0,    if DICO = 'C' or  JOBWL = 'W', */
/*                      a = 2*NW, if DICO = 'D' and JOBWL = 'C'; */
/*                 prefer larger. */

#line 1431 "AB09JD.f"
		i__1 = *ldwork - kw + 1;
#line 1431 "AB09JD.f"
		ab09jw_(jobwl, dico, "I", "N", &nra, m, p, nw, m, &a[nu1 + 
			nu1 * a_dim1], lda, &b[nu1 + b_dim1], ldb, &c__[nu1 * 
			c_dim1 + 1], ldc, &d__[d_offset], ldd, &dwork[kaw], &
			ldabw, temp, &c__1, &dwork[kbw], &ldabw, &dwork[kcw], 
			&ldcdw, &dwork[kdw], &ldcdw, &iwork[1], &dwork[kw], &
			i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
			1);
#line 1437 "AB09JD.f"
		if (ierr != 0) {
#line 1438 "AB09JD.f"
		    if (ierr == 1) {
#line 1439 "AB09JD.f"
			*info = 11;
#line 1440 "AB09JD.f"
		    } else if (ierr == 3) {
#line 1441 "AB09JD.f"
			*info = 15;
#line 1442 "AB09JD.f"
		    }
#line 1443 "AB09JD.f"
		    return 0;
#line 1444 "AB09JD.f"
		}
#line 1445 "AB09JD.f"
	    }

/* Computing MAX */
#line 1447 "AB09JD.f"
	    d__1 = wrkopt, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 1447 "AB09JD.f"
	    wrkopt = max(d__1,d__2);
#line 1448 "AB09JD.f"
	} else {

/*           Compute the projection G1r of V*G1sr*W or */
/*           conj(V)*G1sr*conj(W) containing the poles of G. */

/*           Workspace need: */
/*           real   MAX( 1, NW*(NW+5), NW*N + MAX( b, M*N, P*M ) ) */
/*                    b = 0,    if DICO = 'C' or  JOBWL = 'W', */
/*                    b = 2*NW, if DICO = 'D' and JOBWL = 'C'; */
/*           prefer larger. */

#line 1459 "AB09JD.f"
	    ab09jw_(jobwl, dico, "I", "N", &nra, m, p, nw, m, &a[nu1 + nu1 * 
		    a_dim1], lda, &b[nu1 + b_dim1], ldb, &c__[nu1 * c_dim1 + 
		    1], ldc, &d__[d_offset], ldd, &aw[aw_offset], ldaw, temp, 
		    &c__1, &bw[bw_offset], ldbw, &cw[cw_offset], ldcw, &dw[
		    dw_offset], lddw, &iwork[1], &dwork[1], ldwork, &ierr, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 1464 "AB09JD.f"
	    if (ierr != 0) {
#line 1465 "AB09JD.f"
		if (ierr == 1) {
#line 1466 "AB09JD.f"
		    *info = 4;
#line 1467 "AB09JD.f"
		} else if (ierr == 3) {
#line 1468 "AB09JD.f"
		    *info = 13;
#line 1469 "AB09JD.f"
		}
#line 1470 "AB09JD.f"
		return 0;
#line 1471 "AB09JD.f"
	    }

#line 1473 "AB09JD.f"
	    wrkopt = max(wrkopt,dwork[1]);
#line 1474 "AB09JD.f"
	}
#line 1475 "AB09JD.f"
    }

#line 1477 "AB09JD.f"
    *nr = nra + nu;
#line 1478 "AB09JD.f"
    dwork[1] = wrkopt;

#line 1480 "AB09JD.f"
    return 0;
/* *** Last line of AB09JD *** */
} /* ab09jd_ */

