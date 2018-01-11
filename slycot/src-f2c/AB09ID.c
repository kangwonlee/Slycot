#line 1 "AB09ID.f"
/* AB09ID.f -- translated by f2c (version 20100827).
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

#line 1 "AB09ID.f"
/* Table of constant values */

static doublereal c_b31 = 0.;

/* Subroutine */ int ab09id_(char *dico, char *jobc, char *jobo, char *job, 
	char *weight, char *equil, char *ordsel, integer *n, integer *m, 
	integer *p, integer *nv, integer *pv, integer *nw, integer *mw, 
	integer *nr, doublereal *alpha, doublereal *alphac, doublereal *
	alphao, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *av, integer *ldav, doublereal *bv, integer *ldbv, 
	doublereal *cv, integer *ldcv, doublereal *dv, integer *lddv, 
	doublereal *aw, integer *ldaw, doublereal *bw, integer *ldbw, 
	doublereal *cw, integer *ldcw, doublereal *dw, integer *lddw, integer 
	*ns, doublereal *hsv, doublereal *tol1, doublereal *tol2, integer *
	iwork, doublereal *dwork, integer *ldwork, integer *iwarn, integer *
	info, ftnlen dico_len, ftnlen jobc_len, ftnlen jobo_len, ftnlen 
	job_len, ftnlen weight_len, ftnlen equil_len, ftnlen ordsel_len)
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
    static integer ki, kl, kt, ku, kw, lw, nn, nu, nu1;
    static logical bal;
    static integer lcf;
    static logical bta;
    static integer kbr, kcr, kdr, kbv;
    static logical spa;
    static integer kbw, kcv, kcw, kdv, kti, ldw, nra, nmr, nnq, nnr, nnv, nnw,
	     nvr, nwr, ppv, ierr;
    extern /* Subroutine */ int sb08cd_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, ftnlen), sb08dd_(
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, ftnlen), tb01id_(char *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);
    static logical scale;
    extern /* Subroutine */ int tb01kd_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, ftnlen, ftnlen, ftnlen), tb01pd_(char *, 
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), ab09ix_(char *, char *, char *, char *, integer *
	    , integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen), ab09iy_(char *, char *
	    , char *, char *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr, leftw;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal scalec, scaleo;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static doublereal maxred;
    static logical fixord;
    static integer iwarnl;
    static doublereal alpwrk;
    static logical frwght, rightw;
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

/*     To compute a reduced order model (Ar,Br,Cr,Dr) for an original */
/*     state-space representation (A,B,C,D) by using the frequency */
/*     weighted square-root or balancing-free square-root */
/*     Balance & Truncate (B&T) or Singular Perturbation Approximation */
/*     (SPA) model reduction methods. The algorithm tries to minimize */
/*     the norm of the frequency-weighted error */

/*           ||V*(G-Gr)*W|| */

/*     where G and Gr are the transfer-function matrices of the original */
/*     and reduced order models, respectively, and V and W are */
/*     frequency-weighting transfer-function matrices. V and W must not */
/*     have poles on the imaginary axis for a continuous-time */
/*     system or on the unit circle for a discrete-time system. */
/*     If G is unstable, only the ALPHA-stable part of G is reduced. */
/*     In case of possible pole-zero cancellations in V*G and/or G*W, */
/*     the absolute values of parameters ALPHAO and/or ALPHAC must be */
/*     different from 1. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the original system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     JOBC    CHARACTER*1 */
/*             Specifies the choice of frequency-weighted controllability */
/*             Grammian as follows: */
/*             = 'S': choice corresponding to a combination method [4] */
/*                    of the approaches of Enns [1] and Lin-Chiu [2,3]; */
/*             = 'E': choice corresponding to the stability enhanced */
/*                    modified combination method of [4]. */

/*     JOBO    CHARACTER*1 */
/*             Specifies the choice of frequency-weighted observability */
/*             Grammian as follows: */
/*             = 'S': choice corresponding to a combination method [4] */
/*                    of the approaches of Enns [1] and Lin-Chiu [2,3]; */
/*             = 'E': choice corresponding to the stability enhanced */
/*                    modified combination method of [4]. */

/*     JOB     CHARACTER*1 */
/*             Specifies the model reduction approach to be used */
/*             as follows: */
/*             = 'B':  use the square-root Balance & Truncate method; */
/*             = 'F':  use the balancing-free square-root */
/*                     Balance & Truncate method; */
/*             = 'S':  use the square-root Singular Perturbation */
/*                     Approximation method; */
/*             = 'P':  use the balancing-free square-root */
/*                     Singular Perturbation Approximation method. */

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

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     NV      (input) INTEGER */
/*             The order of the matrix AV. Also the number of rows of */
/*             the matrix BV and the number of columns of the matrix CV. */
/*             NV represents the dimension of the state vector of the */
/*             system with the transfer-function matrix V.  NV >= 0. */

/*     PV      (input) INTEGER */
/*             The number of rows of the matrices CV and DV.  PV >= 0. */
/*             PV represents the dimension of the output vector of the */
/*             system with the transfer-function matrix V. */

/*     NW      (input) INTEGER */
/*             The order of the matrix AW. Also the number of rows of */
/*             the matrix BW and the number of columns of the matrix CW. */
/*             NW represents the dimension of the state vector of the */
/*             system with the transfer-function matrix W.  NW >= 0. */

/*     MW      (input) INTEGER */
/*             The number of columns of the matrices BW and DW.  MW >= 0. */
/*             MW represents the dimension of the input vector of the */
/*             system with the transfer-function matrix W. */

/*     NR      (input/output) INTEGER */
/*             On entry with ORDSEL = 'F', NR is the desired order of the */
/*             resulting reduced order system.  0 <= NR <= N. */
/*             On exit, if INFO = 0, NR is the order of the resulting */
/*             reduced order model. For a system with NU ALPHA-unstable */
/*             eigenvalues and NS ALPHA-stable eigenvalues (NU+NS = N), */
/*             NR is set as follows: if ORDSEL = 'F', NR is equal to */
/*             NU+MIN(MAX(0,NR-NU),NMIN), where NR is the desired order */
/*             on entry, NMIN is the number of frequency-weighted Hankel */
/*             singular values greater than NS*EPS*S1, EPS is the */
/*             machine precision (see LAPACK Library Routine DLAMCH) */
/*             and S1 is the largest Hankel singular value (computed */
/*             in HSV(1)); NR can be further reduced to ensure */
/*             HSV(NR-NU) > HSV(NR+1-NU); */
/*             if ORDSEL = 'A', NR is the sum of NU and the number of */
/*             Hankel singular values greater than MAX(TOL1,NS*EPS*S1). */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             Specifies the ALPHA-stability boundary for the eigenvalues */
/*             of the state dynamics matrix A. For a continuous-time */
/*             system (DICO = 'C'), ALPHA <= 0 is the boundary value for */
/*             the real parts of eigenvalues, while for a discrete-time */
/*             system (DICO = 'D'), 0 <= ALPHA <= 1 represents the */
/*             boundary value for the moduli of eigenvalues. */
/*             The ALPHA-stability domain does not include the boundary. */

/*     ALPHAC  (input) DOUBLE PRECISION */
/*             Combination method parameter for defining the */
/*             frequency-weighted controllability Grammian (see METHOD); */
/*             ABS(ALPHAC) <= 1. */

/*     ALPHAO  (input) DOUBLE PRECISION */
/*             Combination method parameter for defining the */
/*             frequency-weighted observability Grammian (see METHOD); */
/*             ABS(ALPHAO) <= 1. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix A. */
/*             On exit, if INFO = 0, the leading NR-by-NR part of this */
/*             array contains the state dynamics matrix Ar of the */
/*             reduced order system. */
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
/*             part of this array must contain the state matrix AV of */
/*             the system with the transfer-function matrix V. */
/*             On exit, if WEIGHT = 'L' or 'B', MIN(N,M,P) > 0 and */
/*             INFO = 0, the leading NVR-by-NVR part of this array */
/*             contains the state matrix of a minimal realization of V */
/*             in a real Schur form. NVR is returned in IWORK(2). */
/*             AV is not referenced if WEIGHT = 'R' or 'N', */
/*             or MIN(N,M,P) = 0. */

/*     LDAV    INTEGER */
/*             The leading dimension of array AV. */
/*             LDAV >= MAX(1,NV), if WEIGHT = 'L' or 'B'; */
/*             LDAV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     BV      (input/output) DOUBLE PRECISION array, dimension (LDBV,P) */
/*             On entry, if WEIGHT = 'L' or 'B', the leading NV-by-P part */
/*             of this array must contain the input matrix BV of the */
/*             system with the transfer-function matrix V. */
/*             On exit, if WEIGHT = 'L' or 'B', MIN(N,M,P) > 0 and */
/*             INFO = 0, the leading NVR-by-P part of this array contains */
/*             the input matrix of a minimal realization of V. */
/*             BV is not referenced if WEIGHT = 'R' or 'N', */
/*             or MIN(N,M,P) = 0. */

/*     LDBV    INTEGER */
/*             The leading dimension of array BV. */
/*             LDBV >= MAX(1,NV), if WEIGHT = 'L' or 'B'; */
/*             LDBV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     CV      (input/output) DOUBLE PRECISION array, dimension (LDCV,NV) */
/*             On entry, if WEIGHT = 'L' or 'B', the leading PV-by-NV */
/*             part of this array must contain the output matrix CV of */
/*             the system with the transfer-function matrix V. */
/*             On exit, if WEIGHT = 'L' or 'B', MIN(N,M,P) > 0 and */
/*             INFO = 0, the leading PV-by-NVR part of this array */
/*             contains the output matrix of a minimal realization of V. */
/*             CV is not referenced if WEIGHT = 'R' or 'N', */
/*             or MIN(N,M,P) = 0. */

/*     LDCV    INTEGER */
/*             The leading dimension of array CV. */
/*             LDCV >= MAX(1,PV), if WEIGHT = 'L' or 'B'; */
/*             LDCV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     DV      (input) DOUBLE PRECISION array, dimension (LDDV,P) */
/*             If WEIGHT = 'L' or 'B', the leading PV-by-P part of this */
/*             array must contain the feedthrough matrix DV of the system */
/*             with the transfer-function matrix V. */
/*             DV is not referenced if WEIGHT = 'R' or 'N', */
/*             or MIN(N,M,P) = 0. */

/*     LDDV    INTEGER */
/*             The leading dimension of array DV. */
/*             LDDV >= MAX(1,PV), if WEIGHT = 'L' or 'B'; */
/*             LDDV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     AW      (input/output) DOUBLE PRECISION array, dimension (LDAW,NW) */
/*             On entry, if WEIGHT = 'R' or 'B', the leading NW-by-NW */
/*             part of this array must contain the state matrix AW of */
/*             the system with the transfer-function matrix W. */
/*             On exit, if WEIGHT = 'R' or 'B', MIN(N,M,P) > 0 and */
/*             INFO = 0, the leading NWR-by-NWR part of this array */
/*             contains the state matrix of a minimal realization of W */
/*             in a real Schur form. NWR is returned in IWORK(3). */
/*             AW is not referenced if WEIGHT = 'L' or 'N', */
/*             or MIN(N,M,P) = 0. */

/*     LDAW    INTEGER */
/*             The leading dimension of array AW. */
/*             LDAW >= MAX(1,NW), if WEIGHT = 'R' or 'B'; */
/*             LDAW >= 1,         if WEIGHT = 'L' or 'N'. */

/*     BW      (input/output) DOUBLE PRECISION array, dimension (LDBW,MW) */
/*             On entry, if WEIGHT = 'R' or 'B', the leading NW-by-MW */
/*             part of this array must contain the input matrix BW of the */
/*             system with the transfer-function matrix W. */
/*             On exit, if WEIGHT = 'R' or 'B', MIN(N,M,P) > 0 and */
/*             INFO = 0, the leading NWR-by-MW part of this array */
/*             contains the input matrix of a minimal realization of W. */
/*             BW is not referenced if WEIGHT = 'L' or 'N', */
/*             or MIN(N,M,P) = 0. */

/*     LDBW    INTEGER */
/*             The leading dimension of array BW. */
/*             LDBW >= MAX(1,NW), if WEIGHT = 'R' or 'B'; */
/*             LDBW >= 1,         if WEIGHT = 'L' or 'N'. */

/*     CW      (input/output) DOUBLE PRECISION array, dimension (LDCW,NW) */
/*             On entry, if WEIGHT = 'R' or 'B', the leading M-by-NW part */
/*             of this array must contain the output matrix CW of the */
/*             system with the transfer-function matrix W. */
/*             On exit, if WEIGHT = 'R' or 'B', MIN(N,M,P) > 0 and */
/*             INFO = 0, the leading M-by-NWR part of this array contains */
/*             the output matrix of a minimal realization of W. */
/*             CW is not referenced if WEIGHT = 'L' or 'N', */
/*             or MIN(N,M,P) = 0. */

/*     LDCW    INTEGER */
/*             The leading dimension of array CW. */
/*             LDCW >= MAX(1,M), if WEIGHT = 'R' or 'B'; */
/*             LDCW >= 1,        if WEIGHT = 'L' or 'N'. */

/*     DW      (input) DOUBLE PRECISION array, dimension (LDDW,MW) */
/*             If WEIGHT = 'R' or 'B', the leading M-by-MW part of this */
/*             array must contain the feedthrough matrix DW of the system */
/*             with the transfer-function matrix W. */
/*             DW is not referenced if WEIGHT = 'L' or 'N', */
/*             or MIN(N,M,P) = 0. */

/*     LDDW    INTEGER */
/*             The leading dimension of array DW. */
/*             LDDW >= MAX(1,M), if WEIGHT = 'R' or 'B'; */
/*             LDDW >= 1,        if WEIGHT = 'L' or 'N'. */

/*     NS      (output) INTEGER */
/*             The dimension of the ALPHA-stable subsystem. */

/*     HSV     (output) DOUBLE PRECISION array, dimension (N) */
/*             If INFO = 0, the leading NS elements of this array contain */
/*             the frequency-weighted Hankel singular values, ordered */
/*             decreasingly, of the ALPHA-stable part of the original */
/*             system. */

/*     Tolerances */

/*     TOL1    DOUBLE PRECISION */
/*             If ORDSEL = 'A', TOL1 contains the tolerance for */
/*             determining the order of reduced system. */
/*             For model reduction, the recommended value is */
/*             TOL1 = c*S1, where c is a constant in the */
/*             interval [0.00001,0.001], and S1 is the largest */
/*             frequency-weighted Hankel singular value of the */
/*             ALPHA-stable part of the original system (computed */
/*             in HSV(1)). */
/*             If TOL1 <= 0 on entry, the used default value is */
/*             TOL1 = NS*EPS*S1, where NS is the number of */
/*             ALPHA-stable eigenvalues of A and EPS is the machine */
/*             precision (see LAPACK Library Routine DLAMCH). */
/*             If ORDSEL = 'F', the value of TOL1 is ignored. */

/*     TOL2    DOUBLE PRECISION */
/*             The tolerance for determining the order of a minimal */
/*             realization of the ALPHA-stable part of the given system. */
/*             The recommended value is TOL2 = NS*EPS*S1. */
/*             This value is used by default if TOL2 <= 0 on entry. */
/*             If TOL2 > 0 and ORDSEL = 'A', then TOL2 <= TOL1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension */
/*             ( MAX( 3, LIWRK1, LIWRK2, LIWRK3 ) ), where */
/*             LIWRK1 = 0,             if JOB = 'B'; */
/*             LIWRK1 = N,             if JOB = 'F'; */
/*             LIWRK1 = 2*N,           if JOB = 'S' or 'P'; */
/*             LIWRK2 = 0,             if WEIGHT = 'R' or 'N' or  NV = 0; */
/*             LIWRK2 = NV+MAX(P,PV),  if WEIGHT = 'L' or 'B' and NV > 0; */
/*             LIWRK3 = 0,             if WEIGHT = 'L' or 'N' or  NW = 0; */
/*             LIWRK3 = NW+MAX(M,MW),  if WEIGHT = 'R' or 'B' and NW > 0. */
/*             On exit, if INFO = 0, IWORK(1) contains the order of a */
/*             minimal realization of the stable part of the system, */
/*             IWORK(2) and IWORK(3) contain the actual orders */
/*             of the state space realizations of V and W, respectively. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( LMINL, LMINR, LRCF, */
/*                            2*N*N + MAX( 1, LLEFT, LRIGHT, 2*N*N+5*N, */
/*                                         N*MAX(M,P) ) ), */
/*             where */
/*             LMINL  = 0, if WEIGHT = 'R' or 'N' or NV = 0; otherwise, */
/*             LMINL  = MAX(LLCF,NV+MAX(NV,3*P))           if P =  PV; */
/*             LMINL  = MAX(P,PV)*(2*NV+MAX(P,PV))+ */
/*                      MAX(LLCF,NV+MAX(NV,3*P,3*PV))      if P <> PV; */
/*             LRCF   = 0, and */
/*             LMINR  = 0, if WEIGHT = 'L' or 'N' or NW = 0; otherwise, */
/*             LMINR  = NW+MAX(NW,3*M)                     if M =  MW; */
/*             LMINR  = 2*NW*MAX(M,MW)+NW+MAX(NW,3*M,3*MW) if M <> MW; */
/*             LLCF   = PV*(NV+PV)+PV*NV+MAX(NV*(NV+5), PV*(PV+2), */
/*                                           4*PV, 4*P); */
/*             LRCF   = MW*(NW+MW)+MAX(NW*(NW+5),MW*(MW+2),4*MW,4*M) */
/*             LLEFT  = (N+NV)*(N+NV+MAX(N+NV,PV)+5) */
/*                              if WEIGHT = 'L' or 'B' and PV > 0; */
/*             LLEFT  = N*(P+5) if WEIGHT = 'R' or 'N' or  PV = 0; */
/*             LRIGHT = (N+NW)*(N+NW+MAX(N+NW,MW)+5) */
/*                              if WEIGHT = 'R' or 'B' and MW > 0; */
/*             LRIGHT = N*(M+5) if WEIGHT = 'L' or 'N' or  MW = 0. */
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
/*             = 2:  with ORDSEL = 'F', the selected order NR corresponds */
/*                   to repeated singular values for the ALPHA-stable */
/*                   part, which are neither all included nor all */
/*                   excluded from the reduced model; in this case, the */
/*                   resulting NR is automatically decreased to exclude */
/*                   all repeated singular values; */
/*             = 3:  with ORDSEL = 'F', the selected order NR is less */
/*                   than the order of the ALPHA-unstable part of the */
/*                   given system; in this case NR is set equal to the */
/*                   order of the ALPHA-unstable part. */
/*             = 10+K:  K violations of the numerical stability condition */
/*                   occured during the assignment of eigenvalues in the */
/*                   SLICOT Library routines SB08CD and/or SB08DD. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the computation of the ordered real Schur form of A */
/*                   failed; */
/*             = 2:  the separation of the ALPHA-stable/unstable */
/*                   diagonal blocks failed because of very close */
/*                   eigenvalues; */
/*             = 3:  the reduction to a real Schur form of the state */
/*                   matrix of a minimal realization of V failed; */
/*             = 4:  a failure was detected during the ordering of the */
/*                   real Schur form of the state matrix of a minimal */
/*                   realization of V or in the iterative process to */
/*                   compute a left coprime factorization with inner */
/*                   denominator; */
/*             = 5:  if DICO = 'C' and the matrix AV has an observable */
/*                   eigenvalue on the imaginary axis, or DICO = 'D' and */
/*                   AV has an observable eigenvalue on the unit circle; */
/*             = 6:  the reduction to a real Schur form of the state */
/*                   matrix of a minimal realization of W failed; */
/*             = 7:  a failure was detected during the ordering of the */
/*                   real Schur form of the state matrix of a minimal */
/*                   realization of W or in the iterative process to */
/*                   compute a right coprime factorization with inner */
/*                   denominator; */
/*             = 8:  if DICO = 'C' and the matrix AW has a controllable */
/*                   eigenvalue on the imaginary axis, or DICO = 'D' and */
/*                   AW has a controllable eigenvalue on the unit circle; */
/*             = 9:  the computation of eigenvalues failed; */
/*             = 10: the computation of Hankel singular values failed. */

/*     METHOD */

/*     Let G be the transfer-function matrix of the original */
/*     linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t) + Du(t),                          (1) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system. The subroutine AB09ID determines */
/*     the matrices of a reduced order system */

/*          d[z(t)] = Ar*z(t) + Br*u(t) */
/*          yr(t)   = Cr*z(t) + Dr*u(t),                      (2) */

/*     such that the corresponding transfer-function matrix Gr minimizes */
/*     the norm of the frequency-weighted error */

/*             V*(G-Gr)*W,                                    (3) */

/*     where V and W are transfer-function matrices without poles on the */
/*     imaginary axis in continuous-time case or on the unit circle in */
/*     discrete-time case. */

/*     The following procedure is used to reduce G: */

/*     1) Decompose additively G, of order N, as */

/*          G = G1 + G2, */

/*        such that G1 = (A1,B1,C1,D) has only ALPHA-stable poles and */
/*        G2 = (A2,B2,C2,0), of order NU, has only ALPHA-unstable poles. */

/*     2) Compute for G1 a B&T or SPA frequency-weighted approximation */
/*        G1r of order NR-NU using the combination method or the */
/*        modified combination method of [4]. */

/*     3) Assemble the reduced model Gr as */

/*           Gr = G1r + G2. */

/*     For the frequency-weighted reduction of the ALPHA-stable part, */
/*     several methods described in [4] can be employed in conjunction */
/*     with the combination method and modified combination method */
/*     proposed in [4]. */

/*     If JOB = 'B', the square-root B&T method is used. */
/*     If JOB = 'F', the balancing-free square-root version of the */
/*     B&T method is used. */
/*     If JOB = 'S', the square-root version of the SPA method is used. */
/*     If JOB = 'P', the balancing-free square-root version of the */
/*     SPA method is used. */

/*     For each of these methods, left and right truncation matrices */
/*     are determined using the Cholesky factors of an input */
/*     frequency-weighted controllability Grammian P and an output */
/*     frequency-weighted observability Grammian Q. */
/*     P and Q are computed from the controllability Grammian Pi of G*W */
/*     and the observability Grammian Qo of V*G. Using special */
/*     realizations of G*W and V*G, Pi and Qo are computed in the */
/*     partitioned forms */

/*           Pi = ( P11  P12 )   and    Qo = ( Q11  Q12 ) , */
/*                ( P12' P22 )               ( Q12' Q22 ) */

/*     where P11 and Q11 are the leading N-by-N parts of Pi and Qo, */
/*     respectively. Let P0 and Q0 be non-negative definite matrices */
/*     defined below */
/*                                        -1 */
/*            P0 = P11 - ALPHAC**2*P12*P22 *P21 , */
/*                                        -1 */
/*            Q0 = Q11 - ALPHAO**2*Q12*Q22 *Q21. */

/*     The frequency-weighted controllability and observability */
/*     Grammians, P and Q, respectively, are defined as follows: */
/*     P = P0 if JOBC = 'S' (standard combination method [4]); */
/*     P = P1 >= P0 if JOBC = 'E', where P1 is the controllability */
/*     Grammian defined to enforce stability for a modified combination */
/*     method of [4]; */
/*     Q = Q0 if JOBO = 'S' (standard combination method [4]); */
/*     Q = Q1 >= Q0 if JOBO = 'E', where Q1 is the observability */
/*     Grammian defined to enforce stability for a modified combination */
/*     method of [4]. */

/*     If JOBC = JOBO = 'S' and ALPHAC = ALPHAO = 0, the choice of */
/*     Grammians corresponds to the method of Enns [1], while if */
/*     ALPHAC = ALPHAO = 1, the choice of Grammians corresponds */
/*     to the method of Lin and Chiu [2,3]. */

/*     If JOBC = 'S' and ALPHAC = 1, no pole-zero cancellations must */
/*     occur in G*W. If JOBO = 'S' and ALPHAO = 1, no pole-zero */
/*     cancellations must occur in V*G. The presence of pole-zero */
/*     cancellations leads to meaningless results and must be avoided. */

/*     The frequency-weighted Hankel singular values HSV(1), ...., */
/*     HSV(N) are computed as the square roots of the eigenvalues */
/*     of the product P*Q. */

/*     REFERENCES */

/*     [1] Enns, D. */
/*         Model reduction with balanced realizations: An error bound */
/*         and a frequency weighted generalization. */
/*         Proc. 23-th CDC, Las Vegas, pp. 127-132, 1984. */

/*     [2] Lin, C.-A. and Chiu, T.-Y. */
/*         Model reduction via frequency-weighted balanced realization. */
/*         Control Theory and Advanced Technology, vol. 8, */
/*         pp. 341-351, 1992. */

/*     [3] Sreeram, V., Anderson, B.D.O and Madievski, A.G. */
/*         New results on frequency weighted balanced reduction */
/*         technique. */
/*         Proc. ACC, Seattle, Washington, pp. 4004-4009, 1995. */

/*     [4] Varga, A. and Anderson, B.D.O. */
/*         Square-root balancing-free methods for the frequency-weighted */
/*         balancing related model reduction. */
/*         (report in preparation) */

/*     NUMERICAL ASPECTS */

/*     The implemented methods rely on accuracy enhancing square-root */
/*     techniques. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, August 2000. */
/*     D. Sima, University of Bucharest, August 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2000. */

/*     REVISIONS */

/*     A. Varga, Australian National University, Canberra, November 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000, */
/*              Sep. 2001. */

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

#line 652 "AB09ID.f"
    /* Parameter adjustments */
#line 652 "AB09ID.f"
    a_dim1 = *lda;
#line 652 "AB09ID.f"
    a_offset = 1 + a_dim1;
#line 652 "AB09ID.f"
    a -= a_offset;
#line 652 "AB09ID.f"
    b_dim1 = *ldb;
#line 652 "AB09ID.f"
    b_offset = 1 + b_dim1;
#line 652 "AB09ID.f"
    b -= b_offset;
#line 652 "AB09ID.f"
    c_dim1 = *ldc;
#line 652 "AB09ID.f"
    c_offset = 1 + c_dim1;
#line 652 "AB09ID.f"
    c__ -= c_offset;
#line 652 "AB09ID.f"
    d_dim1 = *ldd;
#line 652 "AB09ID.f"
    d_offset = 1 + d_dim1;
#line 652 "AB09ID.f"
    d__ -= d_offset;
#line 652 "AB09ID.f"
    av_dim1 = *ldav;
#line 652 "AB09ID.f"
    av_offset = 1 + av_dim1;
#line 652 "AB09ID.f"
    av -= av_offset;
#line 652 "AB09ID.f"
    bv_dim1 = *ldbv;
#line 652 "AB09ID.f"
    bv_offset = 1 + bv_dim1;
#line 652 "AB09ID.f"
    bv -= bv_offset;
#line 652 "AB09ID.f"
    cv_dim1 = *ldcv;
#line 652 "AB09ID.f"
    cv_offset = 1 + cv_dim1;
#line 652 "AB09ID.f"
    cv -= cv_offset;
#line 652 "AB09ID.f"
    dv_dim1 = *lddv;
#line 652 "AB09ID.f"
    dv_offset = 1 + dv_dim1;
#line 652 "AB09ID.f"
    dv -= dv_offset;
#line 652 "AB09ID.f"
    aw_dim1 = *ldaw;
#line 652 "AB09ID.f"
    aw_offset = 1 + aw_dim1;
#line 652 "AB09ID.f"
    aw -= aw_offset;
#line 652 "AB09ID.f"
    bw_dim1 = *ldbw;
#line 652 "AB09ID.f"
    bw_offset = 1 + bw_dim1;
#line 652 "AB09ID.f"
    bw -= bw_offset;
#line 652 "AB09ID.f"
    cw_dim1 = *ldcw;
#line 652 "AB09ID.f"
    cw_offset = 1 + cw_dim1;
#line 652 "AB09ID.f"
    cw -= cw_offset;
#line 652 "AB09ID.f"
    dw_dim1 = *lddw;
#line 652 "AB09ID.f"
    dw_offset = 1 + dw_dim1;
#line 652 "AB09ID.f"
    dw -= dw_offset;
#line 652 "AB09ID.f"
    --hsv;
#line 652 "AB09ID.f"
    --iwork;
#line 652 "AB09ID.f"
    --dwork;
#line 652 "AB09ID.f"

#line 652 "AB09ID.f"
    /* Function Body */
#line 652 "AB09ID.f"
    *info = 0;
#line 653 "AB09ID.f"
    *iwarn = 0;
#line 654 "AB09ID.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 655 "AB09ID.f"
    bta = lsame_(job, "B", (ftnlen)1, (ftnlen)1) || lsame_(job, "F", (ftnlen)
	    1, (ftnlen)1);
#line 656 "AB09ID.f"
    spa = lsame_(job, "S", (ftnlen)1, (ftnlen)1) || lsame_(job, "P", (ftnlen)
	    1, (ftnlen)1);
#line 657 "AB09ID.f"
    bal = lsame_(job, "B", (ftnlen)1, (ftnlen)1) || lsame_(job, "S", (ftnlen)
	    1, (ftnlen)1);
#line 658 "AB09ID.f"
    scale = lsame_(equil, "S", (ftnlen)1, (ftnlen)1);
#line 659 "AB09ID.f"
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);
#line 660 "AB09ID.f"
    leftw = lsame_(weight, "L", (ftnlen)1, (ftnlen)1) || lsame_(weight, "B", (
	    ftnlen)1, (ftnlen)1);
#line 661 "AB09ID.f"
    rightw = lsame_(weight, "R", (ftnlen)1, (ftnlen)1) || lsame_(weight, 
	    "B", (ftnlen)1, (ftnlen)1);
#line 662 "AB09ID.f"
    frwght = leftw || rightw;

#line 664 "AB09ID.f"
    lw = 1;
#line 665 "AB09ID.f"
    nn = *n * *n;
#line 666 "AB09ID.f"
    nnv = *n + *nv;
#line 667 "AB09ID.f"
    nnw = *n + *nw;
#line 668 "AB09ID.f"
    ppv = max(*p,*pv);

#line 670 "AB09ID.f"
    if (leftw && *pv > 0) {
/* Computing MAX */
#line 671 "AB09ID.f"
	i__1 = lw, i__2 = nnv * (nnv + max(nnv,*pv) + 5);
#line 671 "AB09ID.f"
	lw = max(i__1,i__2);
#line 672 "AB09ID.f"
    } else {
/* Computing MAX */
#line 673 "AB09ID.f"
	i__1 = lw, i__2 = *n * (*p + 5);
#line 673 "AB09ID.f"
	lw = max(i__1,i__2);
#line 674 "AB09ID.f"
    }

#line 676 "AB09ID.f"
    if (rightw && *mw > 0) {
/* Computing MAX */
#line 677 "AB09ID.f"
	i__1 = lw, i__2 = nnw * (nnw + max(nnw,*mw) + 5);
#line 677 "AB09ID.f"
	lw = max(i__1,i__2);
#line 678 "AB09ID.f"
    } else {
/* Computing MAX */
#line 679 "AB09ID.f"
	i__1 = lw, i__2 = *n * (*m + 5);
#line 679 "AB09ID.f"
	lw = max(i__1,i__2);
#line 680 "AB09ID.f"
    }
/* Computing MAX */
#line 681 "AB09ID.f"
    i__1 = lw, i__2 = (nn << 1) + *n * 5, i__1 = max(i__1,i__2), i__2 = *n * 
	    max(*m,*p);
#line 681 "AB09ID.f"
    lw = (nn << 1) + max(i__1,i__2);

#line 683 "AB09ID.f"
    if (leftw && *nv > 0) {
/* Computing MAX */
#line 684 "AB09ID.f"
	i__1 = *nv * (*nv + 5), i__2 = *pv * (*pv + 2), i__1 = max(i__1,i__2),
		 i__2 = ppv << 2;
#line 684 "AB09ID.f"
	lcf = *pv * (*nv + *pv) + *pv * *nv + max(i__1,i__2);
#line 686 "AB09ID.f"
	if (*pv == *p) {
/* Computing MAX */
/* Computing MAX */
#line 687 "AB09ID.f"
	    i__3 = *nv, i__4 = *p * 3;
#line 687 "AB09ID.f"
	    i__1 = max(lw,lcf), i__2 = *nv + max(i__3,i__4);
#line 687 "AB09ID.f"
	    lw = max(i__1,i__2);
#line 688 "AB09ID.f"
	} else {
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
#line 689 "AB09ID.f"
	    i__5 = *nv, i__6 = ppv * 3;
#line 689 "AB09ID.f"
	    i__3 = lcf, i__4 = *nv + max(i__5,i__6);
#line 689 "AB09ID.f"
	    i__1 = lw, i__2 = ppv * ((*nv << 1) + ppv) + max(i__3,i__4);
#line 689 "AB09ID.f"
	    lw = max(i__1,i__2);
#line 691 "AB09ID.f"
	}
#line 692 "AB09ID.f"
    }

#line 694 "AB09ID.f"
    if (rightw && *nw > 0) {
#line 695 "AB09ID.f"
	if (*mw == *m) {
/* Computing MAX */
/* Computing MAX */
#line 696 "AB09ID.f"
	    i__3 = *nw, i__4 = *m * 3;
#line 696 "AB09ID.f"
	    i__1 = lw, i__2 = *nw + max(i__3,i__4);
#line 696 "AB09ID.f"
	    lw = max(i__1,i__2);
#line 697 "AB09ID.f"
	} else {
/* Computing MAX */
/* Computing MAX */
#line 698 "AB09ID.f"
	    i__3 = *nw, i__4 = *m * 3, i__3 = max(i__3,i__4), i__4 = *mw * 3;
#line 698 "AB09ID.f"
	    i__1 = lw, i__2 = (*nw << 1) * max(*m,*mw) + *nw + max(i__3,i__4);
#line 698 "AB09ID.f"
	    lw = max(i__1,i__2);
#line 700 "AB09ID.f"
	}
/* Computing MAX */
/* Computing MAX */
#line 701 "AB09ID.f"
	i__3 = *nw * (*nw + 5), i__4 = *mw * (*mw + 2), i__3 = max(i__3,i__4),
		 i__4 = *mw << 2, i__3 = max(i__3,i__4), i__4 = *m << 2;
#line 701 "AB09ID.f"
	i__1 = lw, i__2 = *mw * (*nw + *mw) + max(i__3,i__4);
#line 701 "AB09ID.f"
	lw = max(i__1,i__2);
#line 703 "AB09ID.f"
    }

/*     Check the input scalar arguments. */

#line 707 "AB09ID.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 708 "AB09ID.f"
	*info = -1;
#line 709 "AB09ID.f"
    } else if (! (lsame_(jobc, "S", (ftnlen)1, (ftnlen)1) || lsame_(jobc, 
	    "E", (ftnlen)1, (ftnlen)1))) {
#line 711 "AB09ID.f"
	*info = -2;
#line 712 "AB09ID.f"
    } else if (! (lsame_(jobo, "S", (ftnlen)1, (ftnlen)1) || lsame_(jobo, 
	    "E", (ftnlen)1, (ftnlen)1))) {
#line 714 "AB09ID.f"
	*info = -3;
#line 715 "AB09ID.f"
    } else if (! (bta || spa)) {
#line 716 "AB09ID.f"
	*info = -4;
#line 717 "AB09ID.f"
    } else if (! (frwght || lsame_(weight, "N", (ftnlen)1, (ftnlen)1))) {
#line 718 "AB09ID.f"
	*info = -5;
#line 719 "AB09ID.f"
    } else if (! (scale || lsame_(equil, "N", (ftnlen)1, (ftnlen)1))) {
#line 720 "AB09ID.f"
	*info = -6;
#line 721 "AB09ID.f"
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
#line 722 "AB09ID.f"
	*info = -7;
#line 723 "AB09ID.f"
    } else if (*n < 0) {
#line 724 "AB09ID.f"
	*info = -8;
#line 725 "AB09ID.f"
    } else if (*m < 0) {
#line 726 "AB09ID.f"
	*info = -9;
#line 727 "AB09ID.f"
    } else if (*p < 0) {
#line 728 "AB09ID.f"
	*info = -10;
#line 729 "AB09ID.f"
    } else if (*nv < 0) {
#line 730 "AB09ID.f"
	*info = -11;
#line 731 "AB09ID.f"
    } else if (*pv < 0) {
#line 732 "AB09ID.f"
	*info = -12;
#line 733 "AB09ID.f"
    } else if (*nw < 0) {
#line 734 "AB09ID.f"
	*info = -13;
#line 735 "AB09ID.f"
    } else if (*mw < 0) {
#line 736 "AB09ID.f"
	*info = -14;
#line 737 "AB09ID.f"
    } else if (fixord && (*nr < 0 || *nr > *n)) {
#line 738 "AB09ID.f"
	*info = -15;
#line 739 "AB09ID.f"
    } else if (discr && (*alpha < 0. || *alpha > 1.) || ! discr && *alpha > 
	    0.) {
#line 741 "AB09ID.f"
	*info = -16;
#line 742 "AB09ID.f"
    } else if (abs(*alphac) > 1.) {
#line 743 "AB09ID.f"
	*info = -17;
#line 744 "AB09ID.f"
    } else if (abs(*alphao) > 1.) {
#line 745 "AB09ID.f"
	*info = -18;
#line 746 "AB09ID.f"
    } else if (*lda < max(1,*n)) {
#line 747 "AB09ID.f"
	*info = -20;
#line 748 "AB09ID.f"
    } else if (*ldb < max(1,*n)) {
#line 749 "AB09ID.f"
	*info = -22;
#line 750 "AB09ID.f"
    } else if (*ldc < max(1,*p)) {
#line 751 "AB09ID.f"
	*info = -24;
#line 752 "AB09ID.f"
    } else if (*ldd < max(1,*p)) {
#line 753 "AB09ID.f"
	*info = -26;
#line 754 "AB09ID.f"
    } else if (*ldav < 1 || leftw && *ldav < *nv) {
#line 755 "AB09ID.f"
	*info = -28;
#line 756 "AB09ID.f"
    } else if (*ldbv < 1 || leftw && *ldbv < *nv) {
#line 757 "AB09ID.f"
	*info = -30;
#line 758 "AB09ID.f"
    } else if (*ldcv < 1 || leftw && *ldcv < *pv) {
#line 759 "AB09ID.f"
	*info = -32;
#line 760 "AB09ID.f"
    } else if (*lddv < 1 || leftw && *lddv < *pv) {
#line 761 "AB09ID.f"
	*info = -34;
#line 762 "AB09ID.f"
    } else if (*ldaw < 1 || rightw && *ldaw < *nw) {
#line 763 "AB09ID.f"
	*info = -36;
#line 764 "AB09ID.f"
    } else if (*ldbw < 1 || rightw && *ldbw < *nw) {
#line 765 "AB09ID.f"
	*info = -38;
#line 766 "AB09ID.f"
    } else if (*ldcw < 1 || rightw && *ldcw < *m) {
#line 767 "AB09ID.f"
	*info = -40;
#line 768 "AB09ID.f"
    } else if (*lddw < 1 || rightw && *lddw < *m) {
#line 769 "AB09ID.f"
	*info = -42;
#line 770 "AB09ID.f"
    } else if (*tol2 > 0. && ! fixord && *tol2 > *tol1) {
#line 771 "AB09ID.f"
	*info = -46;
#line 772 "AB09ID.f"
    } else if (*ldwork < lw) {
#line 773 "AB09ID.f"
	*info = -49;
#line 774 "AB09ID.f"
    }

#line 776 "AB09ID.f"
    if (*info != 0) {

/*        Error return. */

#line 780 "AB09ID.f"
	i__1 = -(*info);
#line 780 "AB09ID.f"
	xerbla_("AB09ID", &i__1, (ftnlen)6);
#line 781 "AB09ID.f"
	return 0;
#line 782 "AB09ID.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 786 "AB09ID.f"
    i__1 = min(*n,*m);
#line 786 "AB09ID.f"
    if (min(i__1,*p) == 0) {
#line 787 "AB09ID.f"
	*nr = 0;
#line 788 "AB09ID.f"
	*ns = 0;
#line 789 "AB09ID.f"
	iwork[1] = 0;
#line 790 "AB09ID.f"
	iwork[2] = *nv;
#line 791 "AB09ID.f"
	iwork[3] = *nw;
#line 792 "AB09ID.f"
	dwork[1] = 1.;
#line 793 "AB09ID.f"
	return 0;
#line 794 "AB09ID.f"
    }

#line 796 "AB09ID.f"
    if (scale) {

/*        Scale simultaneously the matrices A, B and C: */
/*        A <- inv(D)*A*D, B <- inv(D)*B and C <- C*D, where D is a */
/*        diagonal matrix. */
/*        Workspace: N. */

#line 803 "AB09ID.f"
	maxred = 100.;
#line 804 "AB09ID.f"
	tb01id_("All", n, m, p, &maxred, &a[a_offset], lda, &b[b_offset], ldb,
		 &c__[c_offset], ldc, &dwork[1], info, (ftnlen)3);
#line 806 "AB09ID.f"
    }

/*     Correct the value of ALPHA to ensure stability. */

#line 810 "AB09ID.f"
    alpwrk = *alpha;
#line 811 "AB09ID.f"
    if (discr) {
#line 812 "AB09ID.f"
	if (*alpha == 1.) {
#line 812 "AB09ID.f"
	    alpwrk = 1. - sqrt(dlamch_("E", (ftnlen)1));
#line 812 "AB09ID.f"
	}
#line 813 "AB09ID.f"
    } else {
#line 814 "AB09ID.f"
	if (*alpha == 0.) {
#line 814 "AB09ID.f"
	    alpwrk = -sqrt(dlamch_("E", (ftnlen)1));
#line 814 "AB09ID.f"
	}
#line 815 "AB09ID.f"
    }

/*     Allocate working storage. */

#line 819 "AB09ID.f"
    ku = 1;
#line 820 "AB09ID.f"
    kl = ku + nn;
#line 821 "AB09ID.f"
    ki = kl + *n;
#line 822 "AB09ID.f"
    kw = ki + *n;

/*     Reduce A to a block-diagonal real Schur form, with the */
/*     ALPHA-unstable part in the leading diagonal position, using a */
/*     non-orthogonal similarity transformation, A <- inv(T)*A*T, and */
/*     apply the transformation to B and C: B <- inv(T)*B and C <- C*T. */

/*     Workspace needed:      N*(N+2); */
/*     Additional workspace:  need   3*N; */
/*                            prefer larger. */

#line 833 "AB09ID.f"
    i__1 = *ldwork - kw + 1;
#line 833 "AB09ID.f"
    tb01kd_(dico, "Unstable", "General", n, m, p, &alpwrk, &a[a_offset], lda, 
	    &b[b_offset], ldb, &c__[c_offset], ldc, &nu, &dwork[ku], n, &
	    dwork[kl], &dwork[ki], &dwork[kw], &i__1, &ierr, (ftnlen)1, (
	    ftnlen)8, (ftnlen)7);

#line 837 "AB09ID.f"
    if (ierr != 0) {
#line 838 "AB09ID.f"
	if (ierr != 3) {
#line 839 "AB09ID.f"
	    *info = 1;
#line 840 "AB09ID.f"
	} else {
#line 841 "AB09ID.f"
	    *info = 2;
#line 842 "AB09ID.f"
	}
#line 843 "AB09ID.f"
	return 0;
#line 844 "AB09ID.f"
    }

#line 846 "AB09ID.f"
    wrkopt = (integer) dwork[kw] + kw - 1;

/*     Determine NRA, the desired order for the reduction of stable part. */

#line 850 "AB09ID.f"
    iwarnl = 0;
#line 851 "AB09ID.f"
    *ns = *n - nu;
#line 852 "AB09ID.f"
    if (fixord) {
/* Computing MAX */
#line 853 "AB09ID.f"
	i__1 = 0, i__2 = *nr - nu;
#line 853 "AB09ID.f"
	nra = max(i__1,i__2);
#line 854 "AB09ID.f"
	if (*nr < nu) {
#line 854 "AB09ID.f"
	    iwarnl = 3;
#line 854 "AB09ID.f"
	}
#line 856 "AB09ID.f"
    } else {
#line 857 "AB09ID.f"
	nra = 0;
#line 858 "AB09ID.f"
    }

/*     Finish if only unstable part is present. */

#line 862 "AB09ID.f"
    if (*ns == 0) {
#line 863 "AB09ID.f"
	*nr = nu;
#line 864 "AB09ID.f"
	dwork[1] = (doublereal) wrkopt;
#line 865 "AB09ID.f"
	iwork[1] = 0;
#line 866 "AB09ID.f"
	iwork[2] = *nv;
#line 867 "AB09ID.f"
	iwork[3] = *nw;
#line 868 "AB09ID.f"
	return 0;
#line 869 "AB09ID.f"
    }

#line 871 "AB09ID.f"
    nvr = *nv;
#line 872 "AB09ID.f"
    if (leftw && *nv > 0) {

/*        Compute a left-coprime factorization with inner denominator */
/*        of a minimal realization of V. The resulting AV is in */
/*        real Schur form. */
/*        Workspace needed:   real  LV+MAX( 1, LCF, */
/*                                          NV + MAX( NV, 3*P, 3*PV ) ), */
/*                                  where */
/*                                  LV = 0 if P = PV and */
/*                                  LV = MAX(P,PV)*(2*NV+MAX(P,PV)) */
/*                                         otherwise; */
/*                                  LCF = PV*(NV+PV) + */
/*                                        MAX( 1, PV*NV + MAX( NV*(NV+5), */
/*                                             PV*(PV+2),4*PV,4*P ) ); */
/*                                  prefer larger; */
/*                          integer NV + MAX(P,PV). */

#line 889 "AB09ID.f"
	if (*p == *pv) {
#line 890 "AB09ID.f"
	    kw = 1;
#line 891 "AB09ID.f"
	    tb01pd_("Minimal", "Scale", nv, p, pv, &av[av_offset], ldav, &bv[
		    bv_offset], ldbv, &cv[cv_offset], ldcv, &nvr, &c_b31, &
		    iwork[1], &dwork[1], ldwork, info, (ftnlen)7, (ftnlen)5);
/* Computing MAX */
#line 894 "AB09ID.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 894 "AB09ID.f"
	    wrkopt = max(i__1,i__2);
#line 895 "AB09ID.f"
	    kbr = 1;
#line 896 "AB09ID.f"
	    kdr = kbr + *pv * nvr;
#line 897 "AB09ID.f"
	    kw = kdr + *pv * *pv;
#line 898 "AB09ID.f"
	    i__1 = max(1,nvr);
#line 898 "AB09ID.f"
	    i__2 = *ldwork - kw + 1;
#line 898 "AB09ID.f"
	    sb08cd_(dico, &nvr, p, pv, &av[av_offset], ldav, &bv[bv_offset], 
		    ldbv, &cv[cv_offset], ldcv, &dv[dv_offset], lddv, &nnq, &
		    nnr, &dwork[kbr], &i__1, &dwork[kdr], pv, &c_b31, &dwork[
		    kw], &i__2, iwarn, &ierr, (ftnlen)1);
#line 902 "AB09ID.f"
	} else {
#line 903 "AB09ID.f"
	    ldw = max(*p,*pv);
#line 904 "AB09ID.f"
	    kbv = 1;
#line 905 "AB09ID.f"
	    kcv = kbv + *nv * ldw;
#line 906 "AB09ID.f"
	    kw = kcv + *nv * ldw;
#line 907 "AB09ID.f"
	    dlacpy_("Full", nv, p, &bv[bv_offset], ldbv, &dwork[kbv], nv, (
		    ftnlen)4);
#line 908 "AB09ID.f"
	    dlacpy_("Full", pv, nv, &cv[cv_offset], ldcv, &dwork[kcv], &ldw, (
		    ftnlen)4);
#line 909 "AB09ID.f"
	    i__1 = *ldwork - kw + 1;
#line 909 "AB09ID.f"
	    tb01pd_("Minimal", "Scale", nv, p, pv, &av[av_offset], ldav, &
		    dwork[kbv], nv, &dwork[kcv], &ldw, &nvr, &c_b31, &iwork[1]
		    , &dwork[kw], &i__1, info, (ftnlen)7, (ftnlen)5);
#line 912 "AB09ID.f"
	    kdv = kw;
#line 913 "AB09ID.f"
	    kbr = kdv + ldw * ldw;
#line 914 "AB09ID.f"
	    kdr = kbr + *pv * nvr;
#line 915 "AB09ID.f"
	    kw = kdr + *pv * *pv;
#line 916 "AB09ID.f"
	    dlacpy_("Full", pv, p, &dv[dv_offset], lddv, &dwork[kdv], &ldw, (
		    ftnlen)4);
#line 917 "AB09ID.f"
	    i__1 = max(1,nvr);
#line 917 "AB09ID.f"
	    i__2 = *ldwork - kw + 1;
#line 917 "AB09ID.f"
	    sb08cd_(dico, &nvr, p, pv, &av[av_offset], ldav, &dwork[kbv], nv, 
		    &dwork[kcv], &ldw, &dwork[kdv], &ldw, &nnq, &nnr, &dwork[
		    kbr], &i__1, &dwork[kdr], pv, &c_b31, &dwork[kw], &i__2, 
		    iwarn, &ierr, (ftnlen)1);
#line 921 "AB09ID.f"
	    dlacpy_("Full", &nvr, p, &dwork[kbv], nv, &bv[bv_offset], ldbv, (
		    ftnlen)4);
#line 922 "AB09ID.f"
	    dlacpy_("Full", pv, &nvr, &dwork[kcv], &ldw, &cv[cv_offset], ldcv,
		     (ftnlen)4);
#line 923 "AB09ID.f"
	    dlacpy_("Full", pv, p, &dwork[kdv], &ldw, &dv[dv_offset], lddv, (
		    ftnlen)4);
#line 924 "AB09ID.f"
	}
#line 925 "AB09ID.f"
	if (ierr != 0) {
#line 926 "AB09ID.f"
	    *info = ierr + 2;
#line 927 "AB09ID.f"
	    return 0;
#line 928 "AB09ID.f"
	}
#line 929 "AB09ID.f"
	nvr = nnq;
/* Computing MAX */
#line 930 "AB09ID.f"
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 930 "AB09ID.f"
	wrkopt = max(i__1,i__2);
#line 931 "AB09ID.f"
	if (*iwarn > 0) {
#line 931 "AB09ID.f"
	    *iwarn += 10;
#line 931 "AB09ID.f"
	}
#line 933 "AB09ID.f"
    }

#line 935 "AB09ID.f"
    nwr = *nw;
#line 936 "AB09ID.f"
    if (rightw && *nw > 0) {

/*        Compute a minimal realization of W. */
/*        Workspace needed:   real  LW+MAX(1, NW + MAX(NW, 3*M, 3*MW)); */
/*                                  where */
/*                                  LW = 0,              if M = MW and */
/*                                  LW = 2*NW*MAX(M,MW), otherwise; */
/*                                  prefer larger; */
/*                          integer NW + MAX(M,MW). */

#line 946 "AB09ID.f"
	if (*m == *mw) {
#line 947 "AB09ID.f"
	    kw = 1;
#line 948 "AB09ID.f"
	    tb01pd_("Minimal", "Scale", nw, mw, m, &aw[aw_offset], ldaw, &bw[
		    bw_offset], ldbw, &cw[cw_offset], ldcw, &nwr, &c_b31, &
		    iwork[1], &dwork[1], ldwork, info, (ftnlen)7, (ftnlen)5);
#line 951 "AB09ID.f"
	} else {
#line 952 "AB09ID.f"
	    ldw = max(*m,*mw);
#line 953 "AB09ID.f"
	    kbw = 1;
#line 954 "AB09ID.f"
	    kcw = kbw + *nw * ldw;
#line 955 "AB09ID.f"
	    kw = kcw + *nw * ldw;
#line 956 "AB09ID.f"
	    dlacpy_("Full", nw, mw, &bw[bw_offset], ldbw, &dwork[kbw], nw, (
		    ftnlen)4);
#line 957 "AB09ID.f"
	    dlacpy_("Full", m, nw, &cw[cw_offset], ldcw, &dwork[kcw], &ldw, (
		    ftnlen)4);
#line 958 "AB09ID.f"
	    i__1 = *ldwork - kw + 1;
#line 958 "AB09ID.f"
	    tb01pd_("Minimal", "Scale", nw, mw, m, &aw[aw_offset], ldaw, &
		    dwork[kbw], nw, &dwork[kcw], &ldw, &nwr, &c_b31, &iwork[1]
		    , &dwork[kw], &i__1, info, (ftnlen)7, (ftnlen)5);
#line 961 "AB09ID.f"
	    dlacpy_("Full", &nwr, mw, &dwork[kbw], nw, &bw[bw_offset], ldbw, (
		    ftnlen)4);
#line 962 "AB09ID.f"
	    dlacpy_("Full", m, &nwr, &dwork[kcw], &ldw, &cw[cw_offset], ldcw, 
		    (ftnlen)4);
#line 963 "AB09ID.f"
	}
/* Computing MAX */
#line 964 "AB09ID.f"
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 964 "AB09ID.f"
	wrkopt = max(i__1,i__2);
#line 965 "AB09ID.f"
    }

#line 967 "AB09ID.f"
    if (rightw && nwr > 0) {

/*        Compute a right-coprime factorization with inner denominator */
/*        of the minimal realization of W. The resulting AW is in */
/*        real Schur form. */

/*        Workspace needed:  MW*(NW+MW) + */
/*                           MAX( 1, NW*(NW+5), MW*(MW+2), 4*MW, 4*M ); */
/*                           prefer larger. */

#line 977 "AB09ID.f"
	ldw = max(1,*mw);
#line 978 "AB09ID.f"
	kcr = 1;
#line 979 "AB09ID.f"
	kdr = kcr + nwr * ldw;
#line 980 "AB09ID.f"
	kw = kdr + *mw * ldw;
#line 981 "AB09ID.f"
	i__1 = *ldwork - kw + 1;
#line 981 "AB09ID.f"
	sb08dd_(dico, &nwr, mw, m, &aw[aw_offset], ldaw, &bw[bw_offset], ldbw,
		 &cw[cw_offset], ldcw, &dw[dw_offset], lddw, &nnq, &nnr, &
		dwork[kcr], &ldw, &dwork[kdr], &ldw, &c_b31, &dwork[kw], &
		i__1, iwarn, &ierr, (ftnlen)1);
#line 984 "AB09ID.f"
	if (ierr != 0) {
#line 985 "AB09ID.f"
	    *info = ierr + 5;
#line 986 "AB09ID.f"
	    return 0;
#line 987 "AB09ID.f"
	}
#line 988 "AB09ID.f"
	nwr = nnq;
/* Computing MAX */
#line 989 "AB09ID.f"
	i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 989 "AB09ID.f"
	wrkopt = max(i__1,i__2);
#line 990 "AB09ID.f"
	if (*iwarn > 0) {
#line 990 "AB09ID.f"
	    *iwarn += 10;
#line 990 "AB09ID.f"
	}
#line 992 "AB09ID.f"
    }

#line 994 "AB09ID.f"
    nu1 = nu + 1;

/*     Allocate working storage. */

#line 998 "AB09ID.f"
    kt = 1;
#line 999 "AB09ID.f"
    kti = kt + nn;
#line 1000 "AB09ID.f"
    kw = kti + nn;

/*     Compute in DWORK(KTI) and DWORK(KT) the Cholesky factors S and R */
/*     of the controllability and observability Grammians, respectively. */
/*     Real workspace:    need  2*N*N + MAX( 1, LLEFT, LRIGHT ), */
/*             where */
/*             LLEFT  = (N+NV)*(N+NV+MAX(N+NV,PV)+5) */
/*                              if WEIGHT = 'L' or 'B' and PV > 0; */
/*             LLEFT  = N*(P+5) if WEIGHT = 'R' or 'N' or  PV = 0; */
/*             LRIGHT = (N+NW)*(N+NW+MAX(N+NW,MW)+5) */
/*                              if WEIGHT = 'R' or 'B' and MW > 0; */
/*             LRIGHT = N*(M+5) if WEIGHT = 'L' or 'N' or  MW = 0. */
/*                        prefer larger. */

#line 1014 "AB09ID.f"
    i__1 = *ldwork - kw + 1;
#line 1014 "AB09ID.f"
    ab09iy_(dico, jobc, jobo, weight, ns, m, p, &nvr, pv, &nwr, mw, alphac, 
	    alphao, &a[nu1 + nu1 * a_dim1], lda, &b[nu1 + b_dim1], ldb, &c__[
	    nu1 * c_dim1 + 1], ldc, &av[av_offset], ldav, &bv[bv_offset], 
	    ldbv, &cv[cv_offset], ldcv, &dv[dv_offset], lddv, &aw[aw_offset], 
	    ldaw, &bw[bw_offset], ldbw, &cw[cw_offset], ldcw, &dw[dw_offset], 
	    lddw, &scalec, &scaleo, &dwork[kti], n, &dwork[kt], n, &dwork[kw],
	     &i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1020 "AB09ID.f"
    if (ierr != 0) {
#line 1021 "AB09ID.f"
	*info = 9;
#line 1022 "AB09ID.f"
	return 0;
#line 1023 "AB09ID.f"
    }
/* Computing MAX */
#line 1024 "AB09ID.f"
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 1024 "AB09ID.f"
    wrkopt = max(i__1,i__2);

/*     Compute a BTA or SPA of the stable part. */
/*     Real workspace:  need  2*N*N + MAX( 1, 2*N*N+5*N, N*MAX(M,P) ). */

#line 1029 "AB09ID.f"
    i__1 = *ldwork - kw + 1;
#line 1029 "AB09ID.f"
    ab09ix_(dico, job, "Schur", ordsel, ns, m, p, &nra, &scalec, &scaleo, &a[
	    nu1 + nu1 * a_dim1], lda, &b[nu1 + b_dim1], ldb, &c__[nu1 * 
	    c_dim1 + 1], ldc, &d__[d_offset], ldd, &dwork[kti], n, &dwork[kt],
	     n, &nmr, &hsv[1], tol1, tol2, &iwork[1], &dwork[kw], &i__1, 
	    iwarn, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)5, (ftnlen)1);
#line 1034 "AB09ID.f"
    *iwarn = max(*iwarn,iwarnl);
#line 1035 "AB09ID.f"
    if (ierr != 0) {
#line 1036 "AB09ID.f"
	*info = 10;
#line 1037 "AB09ID.f"
	return 0;
#line 1038 "AB09ID.f"
    }
#line 1039 "AB09ID.f"
    *nr = nra + nu;

/* Computing MAX */
#line 1041 "AB09ID.f"
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 1041 "AB09ID.f"
    dwork[1] = (doublereal) max(i__1,i__2);
#line 1042 "AB09ID.f"
    iwork[1] = nmr;
#line 1043 "AB09ID.f"
    iwork[2] = nvr;
#line 1044 "AB09ID.f"
    iwork[3] = nwr;

#line 1046 "AB09ID.f"
    return 0;
/* *** Last line of AB09ID *** */
} /* ab09id_ */

