#line 1 "AB09KD.f"
/* AB09KD.f -- translated by f2c (version 20100827).
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

#line 1 "AB09KD.f"
/* Subroutine */ int ab09kd_(char *job, char *dico, char *weight, char *equil,
	 char *ordsel, integer *n, integer *nv, integer *nw, integer *m, 
	integer *p, integer *nr, doublereal *alpha, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *av, integer *ldav, 
	doublereal *bv, integer *ldbv, doublereal *cv, integer *ldcv, 
	doublereal *dv, integer *lddv, doublereal *aw, integer *ldaw, 
	doublereal *bw, integer *ldbw, doublereal *cw, integer *ldcw, 
	doublereal *dw, integer *lddw, integer *ns, doublereal *hsv, 
	doublereal *tol1, doublereal *tol2, integer *iwork, doublereal *dwork,
	 integer *ldwork, integer *iwarn, integer *info, ftnlen job_len, 
	ftnlen dico_len, ftnlen weight_len, ftnlen equil_len, ftnlen 
	ordsel_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, av_dim1, av_offset, aw_dim1, aw_offset, b_dim1, 
	    b_offset, bv_dim1, bv_offset, bw_dim1, bw_offset, c_dim1, 
	    c_offset, cv_dim1, cv_offset, cw_dim1, cw_offset, d_dim1, 
	    d_offset, dv_dim1, dv_offset, dw_dim1, dw_offset, i__1, i__2, 
	    i__3, i__4, i__5, i__6;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer ia, ib, ki, kl, ku, kw, lw, nu, nu1, nra, ierr, nmin;
    extern /* Subroutine */ int ab07nd_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *), tb01id_(char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), ab09cx_(char *, char *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, ftnlen, ftnlen), tb01kd_(char *, char *, 
	    char *, integer *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen), 
	    ab09kx_(char *, char *, char *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, ftnlen, 
	    ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    static doublereal rcond;
    static logical conjs, leftw;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal maxred;
    static logical fixord;
    static integer iwarnl;
    static doublereal alpwrk;
    static logical frwght, rightw;
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

/*     To compute a reduced order model (Ar,Br,Cr,Dr) for an original */
/*     state-space representation (A,B,C,D) by using the frequency */
/*     weighted optimal Hankel-norm approximation method. */
/*     The Hankel norm of the weighted error */

/*           V*(G-Gr)*W    or    conj(V)*(G-Gr)*conj(W) */

/*     is minimized, where G and Gr are the transfer-function matrices */
/*     of the original and reduced systems, respectively, and V and W */
/*     are the transfer-function matrices of the left and right frequency */
/*     weights, specified by their state space realizations (AV,BV,CV,DV) */
/*     and (AW,BW,CW,DW), respectively. When minimizing the weighted */
/*     error V*(G-Gr)*W, V and W must be antistable transfer-function */
/*     matrices. When minimizing conj(V)*(G-Gr)*conj(W), V and W must be */
/*     stable transfer-function matrices. */
/*     Additionally, V and W must be invertible transfer-function */
/*     matrices, with the feedthrough matrices DV and DW invertible. */
/*     If the original system is unstable, then the frequency weighted */
/*     Hankel-norm approximation is computed only for the */
/*     ALPHA-stable part of the system. */

/*     For a transfer-function matrix G, conj(G) denotes the conjugate */
/*     of G given by G'(-s) for a continuous-time system or G'(1/z) */
/*     for a discrete-time system. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the frequency-weighting problem as follows: */
/*             = 'N':  solve min||V*(G-Gr)*W||_H; */
/*             = 'C':  solve min||conj(V)*(G-Gr)*conj(W)||_H. */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the original system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     WEIGHT  CHARACTER*1 */
/*             Specifies the type of frequency weighting, as follows: */
/*             = 'N':  no weightings are used (V = I, W = I); */
/*             = 'L':  only left weighting V is used (W = I); */
/*             = 'R':  only right weighting W is used (V = I); */
/*             = 'B':  both left and right weightings V and W are used. */

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
/*             On entry, if WEIGHT = 'L' or 'B', the leading NV-by-NV */
/*             part of this array must contain the state matrix AV of a */
/*             state space realization of the left frequency weighting V. */
/*             On exit, if WEIGHT = 'L' or 'B', and INFO = 0, the leading */
/*             NV-by-NV part of this array contains a real Schur form */
/*             of the state matrix of a state space realization of the */
/*             inverse of V. */
/*             AV is not referenced if WEIGHT = 'R' or 'N'. */

/*     LDAV    INTEGER */
/*             The leading dimension of the array AV. */
/*             LDAV >= MAX(1,NV), if WEIGHT = 'L' or 'B'; */
/*             LDAV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     BV      (input/output) DOUBLE PRECISION array, dimension (LDBV,P) */
/*             On entry, if WEIGHT = 'L' or 'B', the leading NV-by-P part */
/*             of this array must contain the input matrix BV of a state */
/*             space realization of the left frequency weighting V. */
/*             On exit, if WEIGHT = 'L' or 'B', and INFO = 0, the leading */
/*             NV-by-P part of this array contains the input matrix of a */
/*             state space realization of the inverse of V. */
/*             BV is not referenced if WEIGHT = 'R' or 'N'. */

/*     LDBV    INTEGER */
/*             The leading dimension of the array BV. */
/*             LDBV >= MAX(1,NV), if WEIGHT = 'L' or 'B'; */
/*             LDBV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     CV      (input/output) DOUBLE PRECISION array, dimension (LDCV,NV) */
/*             On entry, if WEIGHT = 'L' or 'B', the leading P-by-NV part */
/*             of this array must contain the output matrix CV of a state */
/*             space realization of the left frequency weighting V. */
/*             On exit, if WEIGHT = 'L' or 'B', and INFO = 0, the leading */
/*             P-by-NV part of this array contains the output matrix of a */
/*             state space realization of the inverse of V. */
/*             CV is not referenced if WEIGHT = 'R' or 'N'. */

/*     LDCV    INTEGER */
/*             The leading dimension of the array CV. */
/*             LDCV >= MAX(1,P), if WEIGHT = 'L' or 'B'; */
/*             LDCV >= 1,        if WEIGHT = 'R' or 'N'. */

/*     DV      (input/output) DOUBLE PRECISION array, dimension (LDDV,P) */
/*             On entry, if WEIGHT = 'L' or 'B', the leading P-by-P part */
/*             of this array must contain the feedthrough matrix DV of a */
/*             state space realization of the left frequency weighting V. */
/*             On exit, if WEIGHT = 'L' or 'B', and INFO = 0, the leading */
/*             P-by-P part of this array contains the feedthrough matrix */
/*             of a state space realization of the inverse of V. */
/*             DV is not referenced if WEIGHT = 'R' or 'N'. */

/*     LDDV    INTEGER */
/*             The leading dimension of the array DV. */
/*             LDDV >= MAX(1,P), if WEIGHT = 'L' or 'B'; */
/*             LDDV >= 1,        if WEIGHT = 'R' or 'N'. */

/*     AW      (input/output) DOUBLE PRECISION array, dimension (LDAW,NW) */
/*             On entry, if WEIGHT = 'R' or 'B', the leading NW-by-NW */
/*             part of this array must contain the state matrix AW of */
/*             a state space realization of the right frequency */
/*             weighting W. */
/*             On exit, if WEIGHT = 'R' or 'B', and INFO = 0, the leading */
/*             NW-by-NW part of this array contains a real Schur form of */
/*             the state matrix of a state space realization of the */
/*             inverse of W. */
/*             AW is not referenced if WEIGHT = 'L' or 'N'. */

/*     LDAW    INTEGER */
/*             The leading dimension of the array AW. */
/*             LDAW >= MAX(1,NW), if WEIGHT = 'R' or 'B'; */
/*             LDAW >= 1,         if WEIGHT = 'L' or 'N'. */

/*     BW      (input/output) DOUBLE PRECISION array, dimension (LDBW,M) */
/*             On entry, if WEIGHT = 'R' or 'B', the leading NW-by-M part */
/*             of this array must contain the input matrix BW of a state */
/*             space realization of the right frequency weighting W. */
/*             On exit, if WEIGHT = 'R' or 'B', and INFO = 0, the leading */
/*             NW-by-M part of this array contains the input matrix of a */
/*             state space realization of the inverse of W. */
/*             BW is not referenced if WEIGHT = 'L' or 'N'. */

/*     LDBW    INTEGER */
/*             The leading dimension of the array BW. */
/*             LDBW >= MAX(1,NW), if WEIGHT = 'R' or 'B'; */
/*             LDBW >= 1,         if WEIGHT = 'L' or 'N'. */

/*     CW      (input/output) DOUBLE PRECISION array, dimension (LDCW,NW) */
/*             On entry, if WEIGHT = 'R' or 'B', the leading M-by-NW part */
/*             of this array must contain the output matrix CW of a state */
/*             space realization of the right frequency weighting W. */
/*             On exit, if WEIGHT = 'R' or 'B', and INFO = 0, the leading */
/*             M-by-NW part of this array contains the output matrix of a */
/*             state space realization of the inverse of W. */
/*             CW is not referenced if WEIGHT = 'L' or 'N'. */

/*     LDCW    INTEGER */
/*             The leading dimension of the array CW. */
/*             LDCW >= MAX(1,M), if WEIGHT = 'R' or 'B'; */
/*             LDCW >= 1,        if WEIGHT = 'L' or 'N'. */

/*     DW      (input/output) DOUBLE PRECISION array, dimension (LDDW,M) */
/*             On entry, if WEIGHT = 'R' or 'B', the leading M-by-M part */
/*             of this array must contain the feedthrough matrix DW of */
/*             a state space realization of the right frequency */
/*             weighting W. */
/*             On exit, if WEIGHT = 'R' or 'B', and INFO = 0, the leading */
/*             M-by-M part of this array contains the feedthrough matrix */
/*             of a state space realization of the inverse of W. */
/*             DW is not referenced if WEIGHT = 'L' or 'N'. */

/*     LDDW    INTEGER */
/*             The leading dimension of the array DW. */
/*             LDDW >= MAX(1,M), if WEIGHT = 'R' or 'B'; */
/*             LDDW >= 1,        if WEIGHT = 'L' or 'N'. */

/*     NS      (output) INTEGER */
/*             The dimension of the ALPHA-stable subsystem. */

/*     HSV     (output) DOUBLE PRECISION array, dimension (N) */
/*             If INFO = 0, the leading NS elements of this array contain */
/*             the Hankel singular values, ordered decreasingly, of the */
/*             ALPHA-stable part of the weighted original system. */
/*             HSV(1) is the Hankel norm of the ALPHA-stable weighted */
/*             subsystem. */

/*     Tolerances */

/*     TOL1    DOUBLE PRECISION */
/*             If ORDSEL = 'A', TOL1 contains the tolerance for */
/*             determining the order of reduced system. */
/*             For model reduction, the recommended value is */
/*             TOL1 = c*HNORM(As,Bs,Cs), where c is a constant in the */
/*             interval [0.00001,0.001], and HNORM(As,Bs,Cs) is the */
/*             Hankel-norm of the ALPHA-stable part of the weighted */
/*             original system (computed in HSV(1)). */
/*             If TOL1 <= 0 on entry, the used default value is */
/*             TOL1 = NS*EPS*HNORM(As,Bs,Cs), where NS is the number of */
/*             ALPHA-stable eigenvalues of A and EPS is the machine */
/*             precision (see LAPACK Library Routine DLAMCH). */
/*             If ORDSEL = 'F', the value of TOL1 is ignored. */

/*     TOL2    DOUBLE PRECISION */
/*             The tolerance for determining the order of a minimal */
/*             realization of the ALPHA-stable part of the given system. */
/*             The recommended value is TOL2 = NS*EPS*HNORM(As,Bs,Cs). */
/*             This value is used by default if TOL2 <= 0 on entry. */
/*             If TOL2 > 0 and ORDSEL = 'A', then TOL2 <= TOL1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK = MAX(1,M,c),      if DICO = 'C', */
/*             LIWORK = MAX(1,N,M,c),    if DICO = 'D', */
/*             where  c = 0,             if WEIGHT = 'N', */
/*                    c = 2*P,           if WEIGHT = 'L', */
/*                    c = 2*M,           if WEIGHT = 'R', */
/*                    c = MAX(2*M,2*P),  if WEIGHT = 'B'. */
/*             On exit, if INFO = 0, IWORK(1) contains NMIN, the order of */
/*             the computed minimal realization. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( LDW1, LDW2, LDW3, LDW4 ), where */
/*             LDW1 = 0 if WEIGHT = 'R' or 'N' and */
/*             LDW1 = MAX( NV*(NV+5), NV*N + MAX( a, P*N, P*M ) ) */
/*                    if WEIGHT = 'L' or WEIGHT = 'B', */
/*             LDW2 = 0 if WEIGHT = 'L' or 'N' and */
/*             LDW2 = MAX( NW*(NW+5), NW*N + MAX( b, M*N, P*M ) ) */
/*                    if WEIGHT = 'R' or WEIGHT = 'B', with */
/*                a = 0,    b = 0,     if DICO = 'C' or  JOB = 'N', */
/*                a = 2*NV, b = 2*NW,  if DICO = 'D' and JOB = 'C'; */
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
/*                   system; in this case, the resulting NR is set equal */
/*                   to NSMIN; */
/*             = 2:  with ORDSEL = 'F', the selected order NR is less */
/*                   than the order of the ALPHA-unstable part of the */
/*                   given system; in this case NR is set equal to the */
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
/*             =  3:  the reduction of AV or AV-BV*inv(DV)*CV to a */
/*                    real Schur form failed; */
/*             =  4:  the reduction of AW or AW-BW*inv(DW)*CW to a */
/*                    real Schur form failed; */
/*             =  5:  JOB = 'N' and AV is not antistable, or */
/*                    JOB = 'C' and AV is not stable; */
/*             =  6:  JOB = 'N' and AW is not antistable, or */
/*                    JOB = 'C' and AW is not stable; */
/*             =  7:  the computation of Hankel singular values failed; */
/*             =  8:  the computation of stable projection in the */
/*                    Hankel-norm approximation algorithm failed; */
/*             =  9:  the order of computed stable projection in the */
/*                    Hankel-norm approximation algorithm differs */
/*                    from the order of Hankel-norm approximation; */
/*             = 10:  DV is singular; */
/*             = 11:  DW is singular; */
/*             = 12:  the solution of the Sylvester equation failed */
/*                    because the zeros of V (if JOB = 'N') or of conj(V) */
/*                    (if JOB = 'C') are not distinct from the poles */
/*                    of G1sr (see METHOD); */
/*             = 13:  the solution of the Sylvester equation failed */
/*                    because the zeros of W (if JOB = 'N') or of conj(W) */
/*                    (if JOB = 'C') are not distinct from the poles */
/*                    of G1sr (see METHOD). */

/*     METHOD */

/*     Let G be the transfer-function matrix of the original */
/*     linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t) + Du(t),                          (1) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system. The subroutine AB09KD determines */
/*     the matrices of a reduced order system */

/*          d[z(t)] = Ar*z(t) + Br*u(t) */
/*          yr(t)   = Cr*z(t) + Dr*u(t),                      (2) */

/*     such that the corresponding transfer-function matrix Gr minimizes */
/*     the Hankel-norm of the frequency-weighted error */

/*             V*(G-Gr)*W,                                    (3) */
/*     or */
/*             conj(V)*(G-Gr)*conj(W).                        (4) */

/*     For minimizing (3), V and W are assumed to be antistable, while */
/*     for minimizing (4), V and W are assumed to be stable transfer- */
/*     function matrices. */

/*     Note: conj(G) = G'(-s) for a continuous-time system and */
/*           conj(G) = G'(1/z) for a discrete-time system. */

/*     The following procedure is used to reduce G (see [1]): */

/*     1) Decompose additively G as */

/*          G = G1 + G2, */

/*        such that G1 = (A1,B1,C1,D) has only ALPHA-stable poles and */
/*        G2 = (A2,B2,C2,0) has only ALPHA-unstable poles. */

/*     2) Compute G1s, the stable projection of V*G1*W or */
/*        conj(V)*G1*conj(W), using explicit formulas [4]. */

/*     3) Determine G1sr, the optimal Hankel-norm approximation of G1s */
/*        of order r. */

/*     4) Compute G1r, the stable projection of either inv(V)*G1sr*inv(W) */
/*        or conj(inv(V))*G1sr*conj(inv(W)), using explicit formulas [4]. */

/*     5) Assemble the reduced model Gr as */

/*           Gr = G1r + G2. */

/*     To reduce the weighted ALPHA-stable part G1s at step 3, the */
/*     optimal Hankel-norm approximation method of [2], based on the */
/*     square-root balancing projection formulas of [3], is employed. */

/*     The optimal weighted approximation error satisfies */

/*          HNORM[V*(G-Gr)*W] = S(r+1), */
/*     or */
/*          HNORM[conj(V)*(G-Gr)*conj(W)] = S(r+1), */

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

/*     [3] Tombs M.S. and Postlethwaite I. */
/*         Truncated balanced realization of stable, non-minimal */
/*         state-space systems. */
/*         Int. J. Control, Vol. 46, pp. 1319-1330, 1987. */

/*     [4] Varga A. */
/*         Explicit formulas for an efficient implementation */
/*         of the frequency-weighting model reduction approach. */
/*         Proc. 1993 European Control Conference, Groningen, NL, */
/*         pp. 693-696, 1993. */

/*     NUMERICAL ASPECTS */

/*     The implemented methods rely on an accuracy enhancing square-root */
/*     technique. */
/*                                         3 */
/*     The algorithms require less than 30N  floating point operations. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, April 2000. */
/*     D. Sima, University of Bucharest, May 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, May 2000. */
/*     Based on the RASP routines SFRLW, SFRLW1, SFRRW and SFRRW1, */
/*     by A. Varga, 1992. */

/*     REVISIONS */

/*     A. Varga, Australian National University, Canberra, November 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000. */
/*              Oct. 2001, March 2005. */

/*     KEYWORDS */

/*     Frequency weighting, model reduction, multivariable system, */
/*     state-space model, state-space representation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 563 "AB09KD.f"
    /* Parameter adjustments */
#line 563 "AB09KD.f"
    a_dim1 = *lda;
#line 563 "AB09KD.f"
    a_offset = 1 + a_dim1;
#line 563 "AB09KD.f"
    a -= a_offset;
#line 563 "AB09KD.f"
    b_dim1 = *ldb;
#line 563 "AB09KD.f"
    b_offset = 1 + b_dim1;
#line 563 "AB09KD.f"
    b -= b_offset;
#line 563 "AB09KD.f"
    c_dim1 = *ldc;
#line 563 "AB09KD.f"
    c_offset = 1 + c_dim1;
#line 563 "AB09KD.f"
    c__ -= c_offset;
#line 563 "AB09KD.f"
    d_dim1 = *ldd;
#line 563 "AB09KD.f"
    d_offset = 1 + d_dim1;
#line 563 "AB09KD.f"
    d__ -= d_offset;
#line 563 "AB09KD.f"
    av_dim1 = *ldav;
#line 563 "AB09KD.f"
    av_offset = 1 + av_dim1;
#line 563 "AB09KD.f"
    av -= av_offset;
#line 563 "AB09KD.f"
    bv_dim1 = *ldbv;
#line 563 "AB09KD.f"
    bv_offset = 1 + bv_dim1;
#line 563 "AB09KD.f"
    bv -= bv_offset;
#line 563 "AB09KD.f"
    cv_dim1 = *ldcv;
#line 563 "AB09KD.f"
    cv_offset = 1 + cv_dim1;
#line 563 "AB09KD.f"
    cv -= cv_offset;
#line 563 "AB09KD.f"
    dv_dim1 = *lddv;
#line 563 "AB09KD.f"
    dv_offset = 1 + dv_dim1;
#line 563 "AB09KD.f"
    dv -= dv_offset;
#line 563 "AB09KD.f"
    aw_dim1 = *ldaw;
#line 563 "AB09KD.f"
    aw_offset = 1 + aw_dim1;
#line 563 "AB09KD.f"
    aw -= aw_offset;
#line 563 "AB09KD.f"
    bw_dim1 = *ldbw;
#line 563 "AB09KD.f"
    bw_offset = 1 + bw_dim1;
#line 563 "AB09KD.f"
    bw -= bw_offset;
#line 563 "AB09KD.f"
    cw_dim1 = *ldcw;
#line 563 "AB09KD.f"
    cw_offset = 1 + cw_dim1;
#line 563 "AB09KD.f"
    cw -= cw_offset;
#line 563 "AB09KD.f"
    dw_dim1 = *lddw;
#line 563 "AB09KD.f"
    dw_offset = 1 + dw_dim1;
#line 563 "AB09KD.f"
    dw -= dw_offset;
#line 563 "AB09KD.f"
    --hsv;
#line 563 "AB09KD.f"
    --iwork;
#line 563 "AB09KD.f"
    --dwork;
#line 563 "AB09KD.f"

#line 563 "AB09KD.f"
    /* Function Body */
#line 563 "AB09KD.f"
    *info = 0;
#line 564 "AB09KD.f"
    *iwarn = 0;
#line 565 "AB09KD.f"
    conjs = lsame_(job, "C", (ftnlen)1, (ftnlen)1);
#line 566 "AB09KD.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 567 "AB09KD.f"
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);
#line 568 "AB09KD.f"
    leftw = lsame_(weight, "L", (ftnlen)1, (ftnlen)1) || lsame_(weight, "B", (
	    ftnlen)1, (ftnlen)1);
#line 569 "AB09KD.f"
    rightw = lsame_(weight, "R", (ftnlen)1, (ftnlen)1) || lsame_(weight, 
	    "B", (ftnlen)1, (ftnlen)1);
#line 570 "AB09KD.f"
    frwght = leftw || rightw;

#line 572 "AB09KD.f"
    if (discr && conjs) {
#line 573 "AB09KD.f"
	ia = *nv << 1;
#line 574 "AB09KD.f"
	ib = *nw << 1;
#line 575 "AB09KD.f"
    } else {
#line 576 "AB09KD.f"
	ia = 0;
#line 577 "AB09KD.f"
	ib = 0;
#line 578 "AB09KD.f"
    }
#line 579 "AB09KD.f"
    lw = 1;
#line 580 "AB09KD.f"
    if (leftw) {
/* Computing MAX */
/* Computing MAX */
#line 580 "AB09KD.f"
	i__3 = ia, i__4 = *p * *n, i__3 = max(i__3,i__4), i__4 = *p * *m;
#line 580 "AB09KD.f"
	i__1 = lw, i__2 = *nv * (*nv + 5), i__1 = max(i__1,i__2), i__2 = *nv *
		 *n + max(i__3,i__4);
#line 580 "AB09KD.f"
	lw = max(i__1,i__2);
#line 580 "AB09KD.f"
    }
#line 582 "AB09KD.f"
    if (rightw) {
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
#line 582 "AB09KD.f"
	i__5 = ib, i__6 = *m * *n, i__5 = max(i__5,i__6), i__6 = *p * *m;
#line 582 "AB09KD.f"
	i__3 = *nw * (*nw + 5), i__4 = *nw * *n + max(i__5,i__6);
#line 582 "AB09KD.f"
	i__1 = lw, i__2 = max(i__3,i__4);
#line 582 "AB09KD.f"
	lw = max(i__1,i__2);
#line 582 "AB09KD.f"
    }
/* Computing MAX */
/* Computing MAX */
#line 584 "AB09KD.f"
    i__3 = max(*n,*m);
#line 584 "AB09KD.f"
    i__1 = lw, i__2 = *n * ((*n << 1) + max(i__3,*p) + 5) + *n * (*n + 1) / 2;
#line 584 "AB09KD.f"
    lw = max(i__1,i__2);
/* Computing MAX */
/* Computing MAX */
#line 585 "AB09KD.f"
    i__3 = *m * 3 + 1, i__4 = min(*n,*m) + *p;
#line 585 "AB09KD.f"
    i__1 = lw, i__2 = *n * (*m + *p + 2) + (*m << 1) * *p + min(*n,*m) + max(
	    i__3,i__4);
#line 585 "AB09KD.f"
    lw = max(i__1,i__2);

/*     Check the input scalar arguments. */

#line 590 "AB09KD.f"
    if (! (lsame_(job, "N", (ftnlen)1, (ftnlen)1) || conjs)) {
#line 591 "AB09KD.f"
	*info = -1;
#line 592 "AB09KD.f"
    } else if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 593 "AB09KD.f"
	*info = -2;
#line 594 "AB09KD.f"
    } else if (! (frwght || lsame_(weight, "N", (ftnlen)1, (ftnlen)1))) {
#line 595 "AB09KD.f"
	*info = -3;
#line 596 "AB09KD.f"
    } else if (! (lsame_(equil, "S", (ftnlen)1, (ftnlen)1) || lsame_(equil, 
	    "N", (ftnlen)1, (ftnlen)1))) {
#line 598 "AB09KD.f"
	*info = -4;
#line 599 "AB09KD.f"
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
#line 600 "AB09KD.f"
	*info = -5;
#line 601 "AB09KD.f"
    } else if (*n < 0) {
#line 602 "AB09KD.f"
	*info = -6;
#line 603 "AB09KD.f"
    } else if (*nv < 0) {
#line 604 "AB09KD.f"
	*info = -7;
#line 605 "AB09KD.f"
    } else if (*nw < 0) {
#line 606 "AB09KD.f"
	*info = -8;
#line 607 "AB09KD.f"
    } else if (*m < 0) {
#line 608 "AB09KD.f"
	*info = -9;
#line 609 "AB09KD.f"
    } else if (*p < 0) {
#line 610 "AB09KD.f"
	*info = -10;
#line 611 "AB09KD.f"
    } else if (fixord && (*nr < 0 || *nr > *n)) {
#line 612 "AB09KD.f"
	*info = -11;
#line 613 "AB09KD.f"
    } else if (discr && (*alpha < 0. || *alpha > 1.) || ! discr && *alpha > 
	    0.) {
#line 615 "AB09KD.f"
	*info = -12;
#line 616 "AB09KD.f"
    } else if (*lda < max(1,*n)) {
#line 617 "AB09KD.f"
	*info = -14;
#line 618 "AB09KD.f"
    } else if (*ldb < max(1,*n)) {
#line 619 "AB09KD.f"
	*info = -16;
#line 620 "AB09KD.f"
    } else if (*ldc < max(1,*p)) {
#line 621 "AB09KD.f"
	*info = -18;
#line 622 "AB09KD.f"
    } else if (*ldd < max(1,*p)) {
#line 623 "AB09KD.f"
	*info = -20;
#line 624 "AB09KD.f"
    } else if (*ldav < 1 || leftw && *ldav < *nv) {
#line 625 "AB09KD.f"
	*info = -22;
#line 626 "AB09KD.f"
    } else if (*ldbv < 1 || leftw && *ldbv < *nv) {
#line 627 "AB09KD.f"
	*info = -24;
#line 628 "AB09KD.f"
    } else if (*ldcv < 1 || leftw && *ldcv < *p) {
#line 629 "AB09KD.f"
	*info = -26;
#line 630 "AB09KD.f"
    } else if (*lddv < 1 || leftw && *lddv < *p) {
#line 631 "AB09KD.f"
	*info = -28;
#line 632 "AB09KD.f"
    } else if (*ldaw < 1 || rightw && *ldaw < *nw) {
#line 633 "AB09KD.f"
	*info = -30;
#line 634 "AB09KD.f"
    } else if (*ldbw < 1 || rightw && *ldbw < *nw) {
#line 635 "AB09KD.f"
	*info = -32;
#line 636 "AB09KD.f"
    } else if (*ldcw < 1 || rightw && *ldcw < *m) {
#line 637 "AB09KD.f"
	*info = -34;
#line 638 "AB09KD.f"
    } else if (*lddw < 1 || rightw && *lddw < *m) {
#line 639 "AB09KD.f"
	*info = -36;
#line 640 "AB09KD.f"
    } else if (*tol2 > 0. && ! fixord && *tol2 > *tol1) {
#line 641 "AB09KD.f"
	*info = -40;
#line 642 "AB09KD.f"
    } else if (*ldwork < lw) {
#line 643 "AB09KD.f"
	*info = -43;
#line 644 "AB09KD.f"
    }

#line 646 "AB09KD.f"
    if (*info != 0) {

/*        Error return. */

#line 650 "AB09KD.f"
	i__1 = -(*info);
#line 650 "AB09KD.f"
	xerbla_("AB09KD", &i__1, (ftnlen)6);
#line 651 "AB09KD.f"
	return 0;
#line 652 "AB09KD.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 656 "AB09KD.f"
    i__1 = min(*n,*m);
#line 656 "AB09KD.f"
    if (min(i__1,*p) == 0) {
#line 657 "AB09KD.f"
	*nr = 0;
#line 658 "AB09KD.f"
	*ns = 0;
#line 659 "AB09KD.f"
	iwork[1] = 0;
#line 660 "AB09KD.f"
	dwork[1] = 1.;
#line 661 "AB09KD.f"
	return 0;
#line 662 "AB09KD.f"
    }

#line 664 "AB09KD.f"
    if (lsame_(equil, "S", (ftnlen)1, (ftnlen)1)) {

/*        Scale simultaneously the matrices A, B and C: */
/*        A <- inv(D)*A*D, B <- inv(D)*B and C <- C*D, where D is a */
/*        diagonal matrix. */
/*        Workspace: N. */

#line 671 "AB09KD.f"
	maxred = 100.;
#line 672 "AB09KD.f"
	tb01id_("All", n, m, p, &maxred, &a[a_offset], lda, &b[b_offset], ldb,
		 &c__[c_offset], ldc, &dwork[1], info, (ftnlen)3);
#line 674 "AB09KD.f"
    }

/*     Correct the value of ALPHA to ensure stability. */

#line 678 "AB09KD.f"
    alpwrk = *alpha;
#line 679 "AB09KD.f"
    if (discr) {
#line 680 "AB09KD.f"
	if (*alpha == 1.) {
#line 680 "AB09KD.f"
	    alpwrk = 1. - sqrt(dlamch_("E", (ftnlen)1));
#line 680 "AB09KD.f"
	}
#line 681 "AB09KD.f"
    } else {
#line 682 "AB09KD.f"
	if (*alpha == 0.) {
#line 682 "AB09KD.f"
	    alpwrk = -sqrt(dlamch_("E", (ftnlen)1));
#line 682 "AB09KD.f"
	}
#line 683 "AB09KD.f"
    }

/*     Allocate working storage. */

#line 687 "AB09KD.f"
    ku = 1;
#line 688 "AB09KD.f"
    kl = ku + *n * *n;
#line 689 "AB09KD.f"
    ki = kl + *n;
#line 690 "AB09KD.f"
    kw = ki + *n;

/*     Reduce A to a block-diagonal real Schur form, with the */
/*     ALPHA-unstable part in the leading diagonal position, using a */
/*     non-orthogonal similarity transformation, A <- inv(T)*A*T, and */
/*     apply the transformation to B and C: B <- inv(T)*B and C <- C*T. */

/*     Workspace needed:      N*(N+2); */
/*     Additional workspace:  need   3*N; */
/*                            prefer larger. */

#line 701 "AB09KD.f"
    i__1 = *ldwork - kw + 1;
#line 701 "AB09KD.f"
    tb01kd_(dico, "Unstable", "General", n, m, p, &alpwrk, &a[a_offset], lda, 
	    &b[b_offset], ldb, &c__[c_offset], ldc, &nu, &dwork[ku], n, &
	    dwork[kl], &dwork[ki], &dwork[kw], &i__1, &ierr, (ftnlen)1, (
	    ftnlen)8, (ftnlen)7);

#line 705 "AB09KD.f"
    if (ierr != 0) {
#line 706 "AB09KD.f"
	if (ierr != 3) {
#line 707 "AB09KD.f"
	    *info = 1;
#line 708 "AB09KD.f"
	} else {
#line 709 "AB09KD.f"
	    *info = 2;
#line 710 "AB09KD.f"
	}
#line 711 "AB09KD.f"
	return 0;
#line 712 "AB09KD.f"
    }

#line 714 "AB09KD.f"
    wrkopt = dwork[kw] + (doublereal) (kw - 1);

/*     Compute the stable projection of the weighted ALPHA-stable part. */

/*     Workspace: need   MAX( 1, LDW1, LDW2 ), */
/*                LDW1 = 0 if WEIGHT = 'R' or 'N' and */
/*                LDW1 = MAX( NV*(NV+5), NV*N + MAX( a, P*N, P*M ) ) */
/*                       if WEIGHT = 'L' or 'B', */
/*                LDW2 = 0 if WEIGHT = 'L' or 'N' and */
/*                LDW2 = MAX( NW*(NW+5), NW*N + MAX( b, M*N, P*M ) ) */
/*                       if WEIGHT = 'R' or 'B', */
/*                where  a = 0,    b = 0,    if DICO = 'C' or  JOB = 'N', */
/*                       a = 2*NV, b = 2*NW, if DICO = 'D' and JOB = 'C'; */
/*                prefer larger. */

#line 729 "AB09KD.f"
    *ns = *n - nu;

/*     Finish if only unstable part is present. */

#line 733 "AB09KD.f"
    if (*ns == 0) {
#line 734 "AB09KD.f"
	*nr = nu;
#line 735 "AB09KD.f"
	iwork[1] = 0;
#line 736 "AB09KD.f"
	dwork[1] = wrkopt;
#line 737 "AB09KD.f"
	return 0;
#line 738 "AB09KD.f"
    }

#line 740 "AB09KD.f"
    nu1 = nu + 1;
#line 741 "AB09KD.f"
    if (frwght) {
#line 742 "AB09KD.f"
	ab09kx_(job, dico, weight, ns, nv, nw, m, p, &a[nu1 + nu1 * a_dim1], 
		lda, &b[nu1 + b_dim1], ldb, &c__[nu1 * c_dim1 + 1], ldc, &d__[
		d_offset], ldd, &av[av_offset], ldav, &bv[bv_offset], ldbv, &
		cv[cv_offset], ldcv, &dv[dv_offset], lddv, &aw[aw_offset], 
		ldaw, &bw[bw_offset], ldbw, &cw[cw_offset], ldcw, &dw[
		dw_offset], lddw, &dwork[1], ldwork, &iwarnl, &ierr, (ftnlen)
		1, (ftnlen)1, (ftnlen)1);

#line 748 "AB09KD.f"
	if (ierr != 0) {

/*           Note: Only IERR = 1 or IERR = 2 are possible. */
/*           Set INFO to 3 or 4. */

#line 753 "AB09KD.f"
	    *info = ierr + 2;
#line 754 "AB09KD.f"
	    return 0;
#line 755 "AB09KD.f"
	}

#line 757 "AB09KD.f"
	if (iwarnl != 0) {

/*           Stability/antistability of V and W are compulsory. */

#line 761 "AB09KD.f"
	    if (iwarnl == 1 || iwarnl == 3) {
#line 762 "AB09KD.f"
		*info = 5;
#line 763 "AB09KD.f"
	    } else {
#line 764 "AB09KD.f"
		*info = 6;
#line 765 "AB09KD.f"
	    }
#line 766 "AB09KD.f"
	    return 0;
#line 767 "AB09KD.f"
	}

#line 769 "AB09KD.f"
	dwork[1] = max(wrkopt,dwork[1]);
#line 770 "AB09KD.f"
    }

/*     Determine a reduced order approximation of the ALPHA-stable part. */

/*     Workspace: need   MAX( LDW3, LDW4 ), */
/*                LDW3 = N*(2*N + MAX(N,M,P) + 5) + N*(N+1)/2, */
/*                LDW4 = N*(M+P+2) + 2*M*P + MIN(N,M) + */
/*                       MAX( 3*M+1, MIN(N,M)+P ); */
/*                prefer larger. */

#line 780 "AB09KD.f"
    iwarnl = 0;
#line 781 "AB09KD.f"
    if (fixord) {
/* Computing MAX */
#line 782 "AB09KD.f"
	i__1 = 0, i__2 = *nr - nu;
#line 782 "AB09KD.f"
	nra = max(i__1,i__2);
#line 783 "AB09KD.f"
	if (nra == 0) {
#line 783 "AB09KD.f"
	    iwarnl = 2;
#line 783 "AB09KD.f"
	}
#line 785 "AB09KD.f"
    } else {
#line 786 "AB09KD.f"
	nra = 0;
#line 787 "AB09KD.f"
    }
#line 788 "AB09KD.f"
    ab09cx_(dico, ordsel, ns, m, p, &nra, &a[nu1 + nu1 * a_dim1], lda, &b[nu1 
	    + b_dim1], ldb, &c__[nu1 * c_dim1 + 1], ldc, &d__[d_offset], ldd, 
	    &hsv[1], tol1, tol2, &iwork[1], &dwork[1], ldwork, iwarn, &ierr, (
	    ftnlen)1, (ftnlen)1);

#line 792 "AB09KD.f"
    *iwarn = max(*iwarn,iwarnl);
#line 793 "AB09KD.f"
    if (ierr != 0) {

/*        Set INFO = 7, 8 or 9. */

#line 797 "AB09KD.f"
	*info = ierr + 5;
#line 798 "AB09KD.f"
	return 0;
#line 799 "AB09KD.f"
    }

#line 801 "AB09KD.f"
    wrkopt = max(wrkopt,dwork[1]);
#line 802 "AB09KD.f"
    nmin = iwork[1];

/*     Compute the state space realizations of the inverses of V and W. */

/*     Integer workspace: need   c, */
/*     Real workspace:    need   MAX(1,2*c), */
/*                        where  c = 0,             if WEIGHT = 'N', */
/*                               c = 2*P,           if WEIGHT = 'L', */
/*                               c = 2*M,           if WEIGHT = 'R', */
/*                               c = MAX(2*M,2*P),  if WEIGHT = 'B'. */

#line 813 "AB09KD.f"
    if (leftw) {
#line 814 "AB09KD.f"
	ab07nd_(nv, p, &av[av_offset], ldav, &bv[bv_offset], ldbv, &cv[
		cv_offset], ldcv, &dv[dv_offset], lddv, &rcond, &iwork[1], &
		dwork[1], ldwork, &ierr);
#line 816 "AB09KD.f"
	if (ierr != 0) {
#line 817 "AB09KD.f"
	    *info = 10;
#line 818 "AB09KD.f"
	    return 0;
#line 819 "AB09KD.f"
	}
#line 820 "AB09KD.f"
    }
#line 821 "AB09KD.f"
    if (rightw) {
#line 822 "AB09KD.f"
	ab07nd_(nw, m, &aw[aw_offset], ldaw, &bw[bw_offset], ldbw, &cw[
		cw_offset], ldcw, &dw[dw_offset], lddw, &rcond, &iwork[1], &
		dwork[1], ldwork, &ierr);
#line 824 "AB09KD.f"
	if (ierr != 0) {
#line 825 "AB09KD.f"
	    *info = 11;
#line 826 "AB09KD.f"
	    return 0;
#line 827 "AB09KD.f"
	}
#line 828 "AB09KD.f"
    }

#line 830 "AB09KD.f"
    wrkopt = max(wrkopt,dwork[1]);

/*     Compute the stable projection of weighted reduced ALPHA-stable */
/*     part. */

#line 835 "AB09KD.f"
    if (frwght) {
#line 836 "AB09KD.f"
	ab09kx_(job, dico, weight, &nra, nv, nw, m, p, &a[nu1 + nu1 * a_dim1],
		 lda, &b[nu1 + b_dim1], ldb, &c__[nu1 * c_dim1 + 1], ldc, &
		d__[d_offset], ldd, &av[av_offset], ldav, &bv[bv_offset], 
		ldbv, &cv[cv_offset], ldcv, &dv[dv_offset], lddv, &aw[
		aw_offset], ldaw, &bw[bw_offset], ldbw, &cw[cw_offset], ldcw, 
		&dw[dw_offset], lddw, &dwork[1], ldwork, &iwarnl, &ierr, (
		ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 842 "AB09KD.f"
	if (ierr != 0) {
#line 843 "AB09KD.f"
	    if (ierr <= 2) {

/*              Set INFO to 3 or 4. */

#line 847 "AB09KD.f"
		*info = ierr + 2;
#line 848 "AB09KD.f"
	    } else {

/*              Set INFO to 12 or 13. */

#line 852 "AB09KD.f"
		*info = ierr + 9;
#line 853 "AB09KD.f"
	    }
#line 854 "AB09KD.f"
	    return 0;
#line 855 "AB09KD.f"
	}
#line 856 "AB09KD.f"
    }

#line 858 "AB09KD.f"
    *nr = nra + nu;
#line 859 "AB09KD.f"
    iwork[1] = nmin;
#line 860 "AB09KD.f"
    dwork[1] = max(wrkopt,dwork[1]);

#line 862 "AB09KD.f"
    return 0;
/* *** Last line of AB09KD *** */
} /* ab09kd_ */

