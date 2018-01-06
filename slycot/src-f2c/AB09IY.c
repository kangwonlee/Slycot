#line 1 "AB09IY.f"
/* AB09IY.f -- translated by f2c (version 20100827).
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

#line 1 "AB09IY.f"
/* Table of constant values */

static doublereal c_b16 = 0.;
static doublereal c_b20 = 1.;
static logical c_false = FALSE_;
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b43 = -1.;
static logical c_true = TRUE_;

/* Subroutine */ int ab09iy_(char *dico, char *jobc, char *jobo, char *weight,
	 integer *n, integer *m, integer *p, integer *nv, integer *pv, 
	integer *nw, integer *mw, doublereal *alphac, doublereal *alphao, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *av, integer *ldav, doublereal *bv, 
	integer *ldbv, doublereal *cv, integer *ldcv, doublereal *dv, integer 
	*lddv, doublereal *aw, integer *ldaw, doublereal *bw, integer *ldbw, 
	doublereal *cw, integer *ldcw, doublereal *dw, integer *lddw, 
	doublereal *scalec, doublereal *scaleo, doublereal *s, integer *lds, 
	doublereal *r__, integer *ldr, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen dico_len, ftnlen jobc_len, ftnlen jobo_len, 
	ftnlen weight_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, av_dim1, av_offset, aw_dim1, aw_offset, b_dim1, 
	    b_offset, bv_dim1, bv_offset, bw_dim1, bw_offset, c_dim1, 
	    c_offset, cv_dim1, cv_offset, cw_dim1, cw_offset, dv_dim1, 
	    dv_offset, dw_dim1, dw_offset, r_dim1, r_offset, s_dim1, s_offset,
	     i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal t;
    static integer ku, kw, lw, kaw, ldu;
    static doublereal dum[1], tol;
    static integer nnv, nnw, ierr, ktau;
    static doublereal work;
    static integer mbbar;
    extern /* Subroutine */ int mb04nd_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    ftnlen), mb04od_(char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    ftnlen), dscal_(integer *, doublereal *, doublereal *, integer *);
    static integer pcbar;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     mb01wd_(char *, char *, char *, char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    extern /* Subroutine */ int sb03ou_(logical *, logical *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, integer *), dcopy_(integer *, doublereal *, integer *,
	     doublereal *, integer *);
    static logical leftw;
    extern /* Subroutine */ int dsyev_(char *, char *, integer *, doublereal *
	    , integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static logical frwght, rightw;


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

/*     To compute for given state-space representations */
/*     (A,B,C,0), (AV,BV,CV,DV), and (AW,BW,CW,DW) of the */
/*     transfer-function matrices G, V and W, respectively, */
/*     the Cholesky factors of the frequency-weighted */
/*     controllability and observability Grammians corresponding */
/*     to a frequency-weighted model reduction problem. */
/*     G, V and W must be stable transfer-function matrices with */
/*     the state matrices A, AV, and AW in real Schur form. */
/*     It is assumed that the state space realizations (AV,BV,CV,DV) */
/*     and (AW,BW,CW,DW) are minimal. In case of possible pole-zero */
/*     cancellations in forming V*G and/or G*W, the parameters for the */
/*     choice of frequency-weighted Grammians ALPHAO and/or ALPHAC, */
/*     respectively, must be different from 1. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the systems as follows: */
/*             = 'C':  G, V and W are continuous-time systems; */
/*             = 'D':  G, V and W are discrete-time systems. */

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

/*     WEIGHT  CHARACTER*1 */
/*             Specifies the type of frequency weighting, as follows: */
/*             = 'N':  no weightings are used (V = I, W = I); */
/*             = 'L':  only left weighting V is used (W = I); */
/*             = 'R':  only right weighting W is used (V = I); */
/*             = 'B':  both left and right weightings V and W are used. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the state-space representation of G, i.e., */
/*             the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrix B and */
/*             the number of rows of the matrices CW and DW.  M >= 0. */
/*             M represents the dimension of the input vector of the */
/*             system with the transfer-function matrix G and */
/*             also the dimension of the output vector of the system */
/*             with the transfer-function matrix W. */

/*     P       (input) INTEGER */
/*             The number of rows of the matrix C and the */
/*             number of columns of the matrices BV and DV.  P >= 0. */
/*             P represents the dimension of the output vector of the */
/*             system with the transfer-function matrix G and */
/*             also the dimension of the input vector of the system */
/*             with the transfer-function matrix V. */

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

/*     ALPHAC  (input) DOUBLE PRECISION */
/*             Combination method parameter for defining the */
/*             frequency-weighted controllability Grammian (see METHOD); */
/*             ABS(ALPHAC) <= 1. */

/*     ALPHAO  (input) DOUBLE PRECISION */
/*             Combination method parameter for defining the */
/*             frequency-weighted observability Grammian (see METHOD); */
/*             ABS(ALPHAO) <= 1. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must */
/*             contain the state matrix A (of the system with the */
/*             transfer-function matrix G) in a real Schur form. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             input/state matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array must contain the */
/*             state/output matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     AV      (input) DOUBLE PRECISION array, dimension (LDAV,NV) */
/*             If WEIGHT = 'L' or 'B', the leading NV-by-NV part of this */
/*             array must contain the state matrix AV (of the system with */
/*             the transfer-function matrix V) in a real Schur form. */
/*             AV is not referenced if WEIGHT = 'R' or 'N'. */

/*     LDAV    INTEGER */
/*             The leading dimension of array AV. */
/*             LDAV >= MAX(1,NV), if WEIGHT = 'L' or 'B'; */
/*             LDAV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     BV      (input) DOUBLE PRECISION array, dimension (LDBV,P) */
/*             If WEIGHT = 'L' or 'B', the leading NV-by-P part of this */
/*             array must contain the input matrix BV of the system with */
/*             the transfer-function matrix V. */
/*             BV is not referenced if WEIGHT = 'R' or 'N'. */

/*     LDBV    INTEGER */
/*             The leading dimension of array BV. */
/*             LDBV >= MAX(1,NV), if WEIGHT = 'L' or 'B'; */
/*             LDBV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     CV      (input) DOUBLE PRECISION array, dimension (LDCV,NV) */
/*             If WEIGHT = 'L' or 'B', the leading PV-by-NV part of this */
/*             array must contain the output matrix CV of the system with */
/*             the transfer-function matrix V. */
/*             CV is not referenced if WEIGHT = 'R' or 'N'. */

/*     LDCV    INTEGER */
/*             The leading dimension of array CV. */
/*             LDCV >= MAX(1,PV), if WEIGHT = 'L' or 'B'; */
/*             LDCV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     DV      (input) DOUBLE PRECISION array, dimension (LDDV,P) */
/*             If WEIGHT = 'L' or 'B', the leading PV-by-P part of this */
/*             array must contain the feedthrough matrix DV of the system */
/*             with the transfer-function matrix V. */
/*             DV is not referenced if WEIGHT = 'R' or 'N'. */

/*     LDDV    INTEGER */
/*             The leading dimension of array DV. */
/*             LDDV >= MAX(1,PV), if WEIGHT = 'L' or 'B'; */
/*             LDDV >= 1,         if WEIGHT = 'R' or 'N'. */

/*     AW      (input) DOUBLE PRECISION array, dimension (LDAW,NW) */
/*             If WEIGHT = 'R' or 'B', the leading NW-by-NW part of this */
/*             array must contain the state matrix AW (of the system with */
/*             the transfer-function matrix W) in a real Schur form. */
/*             AW is not referenced if WEIGHT = 'L' or 'N'. */

/*     LDAW    INTEGER */
/*             The leading dimension of array AW. */
/*             LDAW >= MAX(1,NW), if WEIGHT = 'R' or 'B'; */
/*             LDAW >= 1,         if WEIGHT = 'L' or 'N'. */

/*     BW      (input) DOUBLE PRECISION array, dimension (LDBW,MW) */
/*             If WEIGHT = 'R' or 'B', the leading NW-by-MW part of this */
/*             array must contain the input matrix BW of the system with */
/*             the transfer-function matrix W. */
/*             BW is not referenced if WEIGHT = 'L' or 'N'. */

/*     LDBW    INTEGER */
/*             The leading dimension of array BW. */
/*             LDBW >= MAX(1,NW), if WEIGHT = 'R' or 'B'; */
/*             LDBW >= 1,         if WEIGHT = 'L' or 'N'. */

/*     CW      (input) DOUBLE PRECISION array, dimension (LDCW,NW) */
/*             If WEIGHT = 'R' or 'B', the leading M-by-NW part of this */
/*             array must contain the output matrix CW of the system with */
/*             the transfer-function matrix W. */
/*             CW is not referenced if WEIGHT = 'L' or 'N'. */

/*     LDCW    INTEGER */
/*             The leading dimension of array CW. */
/*             LDCW >= MAX(1,M), if WEIGHT = 'R' or 'B'; */
/*             LDCW >= 1,        if WEIGHT = 'L' or 'N'. */

/*     DW      (input) DOUBLE PRECISION array, dimension (LDDW,MW) */
/*             If WEIGHT = 'R' or 'B', the leading M-by-MW part of this */
/*             array must contain the feedthrough matrix DW of the system */
/*             with the transfer-function matrix W. */
/*             DW is not referenced if WEIGHT = 'L' or 'N'. */

/*     LDDW    INTEGER */
/*             The leading dimension of array DW. */
/*             LDDW >= MAX(1,M), if WEIGHT = 'R' or 'B'; */
/*             LDDW >= 1,        if WEIGHT = 'L' or 'N'. */

/*     SCALEC  (output) DOUBLE PRECISION */
/*             Scaling factor for the controllability Grammian in (1) */
/*             or (3). See METHOD. */

/*     SCALEO  (output) DOUBLE PRECISION */
/*             Scaling factor for the observability Grammian in (2) */
/*             or (4). See METHOD. */

/*     S       (output) DOUBLE PRECISION array, dimension (LDS,N) */
/*             The leading N-by-N upper triangular part of this array */
/*             contains the Cholesky factor S of the frequency-weighted */
/*             cotrollability Grammian P = S*S'. See METHOD. */

/*     LDS     INTEGER */
/*             The leading dimension of array S.  LDS >= MAX(1,N). */

/*     R       (output) DOUBLE PRECISION array, dimension (LDR,N) */
/*             The leading N-by-N upper triangular part of this array */
/*             contains the Cholesky factor R of the frequency-weighted */
/*             observability Grammian Q = R'*R. See METHOD. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,N). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 1, LLEFT, LRIGHT ), */
/*             where */
/*             LLEFT  = (N+NV)*(N+NV+MAX(N+NV,PV)+5) */
/*                              if WEIGHT = 'L' or 'B' and PV > 0; */
/*             LLEFT  = N*(P+5) if WEIGHT = 'R' or 'N' or  PV = 0; */
/*             LRIGHT = (N+NW)*(N+NW+MAX(N+NW,MW)+5) */
/*                              if WEIGHT = 'R' or 'B' and MW > 0; */
/*             LRIGHT = N*(M+5) if WEIGHT = 'L' or 'N' or  MW = 0. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the state matrices A and/or AV are not stable or */
/*                   not in a real Schur form; */
/*             = 2:  if the state matrices A and/or AW are not stable or */
/*                   not in a real Schur form; */
/*             = 3:  eigenvalues computation failure. */

/*     METHOD */

/*     Let Pi = Si*Si' and Qo = Ro'*Ro be the Cholesky factored */
/*     controllability and observability Grammians satisfying */
/*     in the continuous-time case */

/*            Ai*Pi + Pi*Ai' +  scalec^2*Bi*Bi' = 0,       (1) */

/*            Ao'*Qo + Qo*Ao +  scaleo^2*Co'*Co = 0,       (2) */

/*     and in the discrete-time case */

/*            Ai*Pi*Ai' - Pi +  scalec^2*Bi*Bi' = 0,       (3) */

/*            Ao'*Qo*Ao - Qo +  scaleo^2*Co'*Co = 0,       (4) */

/*     where */

/*           Ai = ( A  B*Cw ) ,   Bi = ( B*Dw ) , */
/*                ( 0   Aw  )          (  Bw  ) */

/*           Ao = (  A   0  ) ,   Co = ( Dv*C  Cv ) . */
/*                ( Bv*C Av ) */

/*     Consider the partitioned Grammians */

/*           Pi = ( P11  P12 )   and    Qo = ( Q11  Q12 ) , */
/*                ( P12' P22 )               ( Q12' Q22 ) */

/*     where P11 and Q11 are the leading N-by-N parts of Pi and Qo, */
/*     respectively, and let P0 and Q0 be non-negative definite matrices */
/*     defined in the combination method [4] */
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
/*     ALPHAC = ALPHAO = 1, the choice of Grammians corresponds to the */
/*     method of Lin and Chiu [2,3]. */

/*     The routine computes directly the Cholesky factors S and R */
/*     such that P = S*S' and Q = R'*R according to formulas */
/*     developed in [4]. No matrix inversions are involved. */

/*     REFERENCES */

/*     [1] Enns, D. */
/*         Model reduction with balanced realizations: An error bound */
/*         and a frequency weighted generalization. */
/*         Proc. CDC, Las Vegas, pp. 127-132, 1984. */

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

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, August 2000. */
/*     D. Sima, University of Bucharest, August 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2000. */

/*     REVISIONS */

/*     A. Varga, Australian National University, Canberra, November 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000. */
/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, August 2001. */

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

#line 423 "AB09IY.f"
    /* Parameter adjustments */
#line 423 "AB09IY.f"
    a_dim1 = *lda;
#line 423 "AB09IY.f"
    a_offset = 1 + a_dim1;
#line 423 "AB09IY.f"
    a -= a_offset;
#line 423 "AB09IY.f"
    b_dim1 = *ldb;
#line 423 "AB09IY.f"
    b_offset = 1 + b_dim1;
#line 423 "AB09IY.f"
    b -= b_offset;
#line 423 "AB09IY.f"
    c_dim1 = *ldc;
#line 423 "AB09IY.f"
    c_offset = 1 + c_dim1;
#line 423 "AB09IY.f"
    c__ -= c_offset;
#line 423 "AB09IY.f"
    av_dim1 = *ldav;
#line 423 "AB09IY.f"
    av_offset = 1 + av_dim1;
#line 423 "AB09IY.f"
    av -= av_offset;
#line 423 "AB09IY.f"
    bv_dim1 = *ldbv;
#line 423 "AB09IY.f"
    bv_offset = 1 + bv_dim1;
#line 423 "AB09IY.f"
    bv -= bv_offset;
#line 423 "AB09IY.f"
    cv_dim1 = *ldcv;
#line 423 "AB09IY.f"
    cv_offset = 1 + cv_dim1;
#line 423 "AB09IY.f"
    cv -= cv_offset;
#line 423 "AB09IY.f"
    dv_dim1 = *lddv;
#line 423 "AB09IY.f"
    dv_offset = 1 + dv_dim1;
#line 423 "AB09IY.f"
    dv -= dv_offset;
#line 423 "AB09IY.f"
    aw_dim1 = *ldaw;
#line 423 "AB09IY.f"
    aw_offset = 1 + aw_dim1;
#line 423 "AB09IY.f"
    aw -= aw_offset;
#line 423 "AB09IY.f"
    bw_dim1 = *ldbw;
#line 423 "AB09IY.f"
    bw_offset = 1 + bw_dim1;
#line 423 "AB09IY.f"
    bw -= bw_offset;
#line 423 "AB09IY.f"
    cw_dim1 = *ldcw;
#line 423 "AB09IY.f"
    cw_offset = 1 + cw_dim1;
#line 423 "AB09IY.f"
    cw -= cw_offset;
#line 423 "AB09IY.f"
    dw_dim1 = *lddw;
#line 423 "AB09IY.f"
    dw_offset = 1 + dw_dim1;
#line 423 "AB09IY.f"
    dw -= dw_offset;
#line 423 "AB09IY.f"
    s_dim1 = *lds;
#line 423 "AB09IY.f"
    s_offset = 1 + s_dim1;
#line 423 "AB09IY.f"
    s -= s_offset;
#line 423 "AB09IY.f"
    r_dim1 = *ldr;
#line 423 "AB09IY.f"
    r_offset = 1 + r_dim1;
#line 423 "AB09IY.f"
    r__ -= r_offset;
#line 423 "AB09IY.f"
    --dwork;
#line 423 "AB09IY.f"

#line 423 "AB09IY.f"
    /* Function Body */
#line 423 "AB09IY.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 424 "AB09IY.f"
    leftw = lsame_(weight, "L", (ftnlen)1, (ftnlen)1) || lsame_(weight, "B", (
	    ftnlen)1, (ftnlen)1);
#line 425 "AB09IY.f"
    rightw = lsame_(weight, "R", (ftnlen)1, (ftnlen)1) || lsame_(weight, 
	    "B", (ftnlen)1, (ftnlen)1);
#line 426 "AB09IY.f"
    frwght = leftw || rightw;

#line 428 "AB09IY.f"
    *info = 0;
#line 429 "AB09IY.f"
    lw = 1;
#line 430 "AB09IY.f"
    nnv = *n + *nv;
#line 431 "AB09IY.f"
    nnw = *n + *nw;
#line 432 "AB09IY.f"
    if (leftw && *pv > 0) {
/* Computing MAX */
#line 433 "AB09IY.f"
	i__1 = lw, i__2 = nnv * (nnv + max(nnv,*pv) + 5);
#line 433 "AB09IY.f"
	lw = max(i__1,i__2);
#line 434 "AB09IY.f"
    } else {
/* Computing MAX */
#line 435 "AB09IY.f"
	i__1 = lw, i__2 = *n * (*p + 5);
#line 435 "AB09IY.f"
	lw = max(i__1,i__2);
#line 436 "AB09IY.f"
    }
#line 437 "AB09IY.f"
    if (rightw && *mw > 0) {
/* Computing MAX */
#line 438 "AB09IY.f"
	i__1 = lw, i__2 = nnw * (nnw + max(nnw,*mw) + 5);
#line 438 "AB09IY.f"
	lw = max(i__1,i__2);
#line 439 "AB09IY.f"
    } else {
/* Computing MAX */
#line 440 "AB09IY.f"
	i__1 = lw, i__2 = *n * (*m + 5);
#line 440 "AB09IY.f"
	lw = max(i__1,i__2);
#line 441 "AB09IY.f"
    }

#line 443 "AB09IY.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 444 "AB09IY.f"
	*info = -1;
#line 445 "AB09IY.f"
    } else if (! (lsame_(jobc, "S", (ftnlen)1, (ftnlen)1) || lsame_(jobc, 
	    "E", (ftnlen)1, (ftnlen)1))) {
#line 447 "AB09IY.f"
	*info = -2;
#line 448 "AB09IY.f"
    } else if (! (lsame_(jobo, "S", (ftnlen)1, (ftnlen)1) || lsame_(jobo, 
	    "E", (ftnlen)1, (ftnlen)1))) {
#line 450 "AB09IY.f"
	*info = -3;
#line 451 "AB09IY.f"
    } else if (! (frwght || lsame_(weight, "N", (ftnlen)1, (ftnlen)1))) {
#line 452 "AB09IY.f"
	*info = -4;
#line 453 "AB09IY.f"
    } else if (*n < 0) {
#line 454 "AB09IY.f"
	*info = -5;
#line 455 "AB09IY.f"
    } else if (*m < 0) {
#line 456 "AB09IY.f"
	*info = -6;
#line 457 "AB09IY.f"
    } else if (*p < 0) {
#line 458 "AB09IY.f"
	*info = -7;
#line 459 "AB09IY.f"
    } else if (*nv < 0) {
#line 460 "AB09IY.f"
	*info = -8;
#line 461 "AB09IY.f"
    } else if (*pv < 0) {
#line 462 "AB09IY.f"
	*info = -9;
#line 463 "AB09IY.f"
    } else if (*nw < 0) {
#line 464 "AB09IY.f"
	*info = -10;
#line 465 "AB09IY.f"
    } else if (*mw < 0) {
#line 466 "AB09IY.f"
	*info = -11;
#line 467 "AB09IY.f"
    } else if (abs(*alphac) > 1.) {
#line 468 "AB09IY.f"
	*info = -12;
#line 469 "AB09IY.f"
    } else if (abs(*alphao) > 1.) {
#line 470 "AB09IY.f"
	*info = -13;
#line 471 "AB09IY.f"
    } else if (*lda < max(1,*n)) {
#line 472 "AB09IY.f"
	*info = -15;
#line 473 "AB09IY.f"
    } else if (*ldb < max(1,*n)) {
#line 474 "AB09IY.f"
	*info = -17;
#line 475 "AB09IY.f"
    } else if (*ldc < max(1,*p)) {
#line 476 "AB09IY.f"
	*info = -19;
#line 477 "AB09IY.f"
    } else if (*ldav < 1 || leftw && *ldav < *nv) {
#line 478 "AB09IY.f"
	*info = -21;
#line 479 "AB09IY.f"
    } else if (*ldbv < 1 || leftw && *ldbv < *nv) {
#line 480 "AB09IY.f"
	*info = -23;
#line 481 "AB09IY.f"
    } else if (*ldcv < 1 || leftw && *ldcv < *pv) {
#line 482 "AB09IY.f"
	*info = -25;
#line 483 "AB09IY.f"
    } else if (*lddv < 1 || leftw && *lddv < *pv) {
#line 484 "AB09IY.f"
	*info = -27;
#line 485 "AB09IY.f"
    } else if (*ldaw < 1 || rightw && *ldaw < *nw) {
#line 486 "AB09IY.f"
	*info = -29;
#line 487 "AB09IY.f"
    } else if (*ldbw < 1 || rightw && *ldbw < *nw) {
#line 488 "AB09IY.f"
	*info = -31;
#line 489 "AB09IY.f"
    } else if (*ldcw < 1 || rightw && *ldcw < *m) {
#line 490 "AB09IY.f"
	*info = -33;
#line 491 "AB09IY.f"
    } else if (*lddw < 1 || rightw && *lddw < *m) {
#line 492 "AB09IY.f"
	*info = -35;
#line 493 "AB09IY.f"
    } else if (*lds < max(1,*n)) {
#line 494 "AB09IY.f"
	*info = -39;
#line 495 "AB09IY.f"
    } else if (*ldr < max(1,*n)) {
#line 496 "AB09IY.f"
	*info = -41;
#line 497 "AB09IY.f"
    } else if (*ldwork < lw) {
#line 498 "AB09IY.f"
	*info = -43;
#line 499 "AB09IY.f"
    }

#line 501 "AB09IY.f"
    if (*info != 0) {

/*        Error return. */

#line 505 "AB09IY.f"
	i__1 = -(*info);
#line 505 "AB09IY.f"
	xerbla_("AB09IY", &i__1, (ftnlen)6);
#line 506 "AB09IY.f"
	return 0;
#line 507 "AB09IY.f"
    }

/*     Quick return if possible. */

#line 511 "AB09IY.f"
    *scalec = 1.;
#line 512 "AB09IY.f"
    *scaleo = 1.;
/* Computing MIN */
#line 513 "AB09IY.f"
    i__1 = min(*n,*m);
#line 513 "AB09IY.f"
    if (min(i__1,*p) == 0) {
#line 514 "AB09IY.f"
	dwork[1] = 1.;
#line 515 "AB09IY.f"
	return 0;
#line 516 "AB09IY.f"
    }

#line 518 "AB09IY.f"
    work = 1.;
#line 519 "AB09IY.f"
    if (leftw && *pv > 0) {

/*        Build the extended permuted matrices */

/*           Ao = ( Av  Bv*C ) ,   Co = ( Cv Dv*C ) . */
/*                ( 0     A  ) */

#line 526 "AB09IY.f"
	kaw = 1;
#line 527 "AB09IY.f"
	ku = kaw + nnv * nnv;
#line 528 "AB09IY.f"
	ldu = max(nnv,*pv);
#line 529 "AB09IY.f"
	dlacpy_("Full", nv, nv, &av[av_offset], ldav, &dwork[kaw], &nnv, (
		ftnlen)4);
#line 530 "AB09IY.f"
	dlaset_("Full", n, nv, &c_b16, &c_b16, &dwork[kaw + *nv], &nnv, (
		ftnlen)4);
#line 531 "AB09IY.f"
	dgemm_("No-transpose", "No-transpose", nv, n, p, &c_b20, &bv[
		bv_offset], ldbv, &c__[c_offset], ldc, &c_b16, &dwork[kaw + 
		nnv * *nv], &nnv, (ftnlen)12, (ftnlen)12);
#line 533 "AB09IY.f"
	dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[kaw + nnv * *nv + *nv]
		, &nnv, (ftnlen)4);

#line 535 "AB09IY.f"
	dlacpy_("Full", pv, nv, &cv[cv_offset], ldcv, &dwork[ku], &ldu, (
		ftnlen)4);
#line 536 "AB09IY.f"
	dgemm_("No-transpose", "No-transpose", pv, n, p, &c_b20, &dv[
		dv_offset], lddv, &c__[c_offset], ldc, &c_b16, &dwork[ku + 
		ldu * *nv], &ldu, (ftnlen)12, (ftnlen)12);

/*        Solve for the Cholesky factor Ro of Qo, Qo = Ro'*Ro, */
/*        the continuous-time Lyapunov equation (if DICO = 'C') */

/*            Ao'*Qo + Qo*Ao  +  scaleo^2*Co'*Co = 0, */

/*        or the discrete-time Lyapunov equation (if DICO = 'D') */

/*            Ao'*Qo*Ao - Qo +  scaleo^2*Co'*Co = 0. */

/*        Workspace:  need   (N+NV)*(N+NV+MAX(N+NV,PV)+5); */
/*                           prefer larger. */

#line 551 "AB09IY.f"
	ktau = ku + ldu * nnv;
#line 552 "AB09IY.f"
	kw = ktau + nnv;

#line 554 "AB09IY.f"
	i__1 = *ldwork - kw + 1;
#line 554 "AB09IY.f"
	sb03ou_(&discr, &c_false, &nnv, pv, &dwork[kaw], &nnv, &dwork[ku], &
		ldu, &dwork[ktau], &dwork[ku], &ldu, scaleo, &dwork[kw], &
		i__1, &ierr);

#line 558 "AB09IY.f"
	if (ierr != 0) {
#line 559 "AB09IY.f"
	    *info = 1;
#line 560 "AB09IY.f"
	    return 0;
#line 561 "AB09IY.f"
	}
/* Computing MAX */
#line 562 "AB09IY.f"
	d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 562 "AB09IY.f"
	work = max(d__1,d__2);

/*        Partition Ro as Ro = ( R11 R12 ) and compute R such that */
/*                             (  0  R22 ) */

/*        R'*R = R22'*R22 + (1-ALPHAO**2)*R12'*R12. */

#line 569 "AB09IY.f"
	kw = ku + ldu * *nv + *nv;
#line 570 "AB09IY.f"
	dlacpy_("Upper", n, n, &dwork[kw], &ldu, &r__[r_offset], ldr, (ftnlen)
		5);
#line 571 "AB09IY.f"
	if (*alphao != 0.) {
#line 572 "AB09IY.f"
	    t = sqrt(1. - *alphao * *alphao);
#line 573 "AB09IY.f"
	    i__1 = ku + ldu * (nnv - 1);
#line 573 "AB09IY.f"
	    i__2 = ldu;
#line 573 "AB09IY.f"
	    for (j = ku + ldu * *nv; i__2 < 0 ? j >= i__1 : j <= i__1; j += 
		    i__2) {
#line 574 "AB09IY.f"
		dscal_(nv, &t, &dwork[j], &c__1);
#line 575 "AB09IY.f"
/* L10: */
#line 575 "AB09IY.f"
	    }
#line 576 "AB09IY.f"
	}
#line 577 "AB09IY.f"
	if (*alphao < 1. && *nv > 0) {
#line 578 "AB09IY.f"
	    ktau = 1;
#line 579 "AB09IY.f"
	    mb04od_("Full", n, &c__0, nv, &r__[r_offset], ldr, &dwork[ku + 
		    ldu * *nv], &ldu, dum, &c__1, dum, &c__1, &dwork[ktau], &
		    dwork[kw], (ftnlen)4);

#line 582 "AB09IY.f"
	    i__2 = *n;
#line 582 "AB09IY.f"
	    for (j = 1; j <= i__2; ++j) {
#line 583 "AB09IY.f"
		dwork[j] = r__[j + j * r_dim1];
#line 584 "AB09IY.f"
		i__1 = j;
#line 584 "AB09IY.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 585 "AB09IY.f"
		    if (dwork[i__] < 0.) {
#line 585 "AB09IY.f"
			r__[i__ + j * r_dim1] = -r__[i__ + j * r_dim1];
#line 585 "AB09IY.f"
		    }
#line 586 "AB09IY.f"
/* L20: */
#line 586 "AB09IY.f"
		}
#line 587 "AB09IY.f"
/* L30: */
#line 587 "AB09IY.f"
	    }

#line 589 "AB09IY.f"
	}

#line 591 "AB09IY.f"
	if (lsame_(jobo, "E", (ftnlen)1, (ftnlen)1) && *alphao < 1.) {

/*           Form Y = -A'*(R'*R)-(R'*R)*A if DICO = 'C', or */
/*                Y = -A'*(R'*R)*A+(R'*R) if DICO = 'D'. */

#line 596 "AB09IY.f"
	    dlacpy_("Upper", n, n, &r__[r_offset], ldr, &dwork[ku], n, (
		    ftnlen)5);
#line 597 "AB09IY.f"
	    mb01wd_(dico, "Upper", "No-transpose", "Hessenberg", n, &c_b43, &
		    c_b16, &r__[r_offset], ldr, &dwork[kaw + nnv * *nv + *nv],
		     &nnv, &dwork[ku], n, &ierr, (ftnlen)1, (ftnlen)5, (
		    ftnlen)12, (ftnlen)10);

/*           Compute the eigendecomposition of Y as Y = Z*Sigma*Z'. */

#line 603 "AB09IY.f"
	    ku = *n + 1;
#line 604 "AB09IY.f"
	    i__2 = *ldwork - *n;
#line 604 "AB09IY.f"
	    dsyev_("Vectors", "Upper", n, &r__[r_offset], ldr, &dwork[1], &
		    dwork[ku], &i__2, &ierr, (ftnlen)7, (ftnlen)5);
#line 606 "AB09IY.f"
	    if (ierr > 0) {
#line 607 "AB09IY.f"
		*info = 3;
#line 608 "AB09IY.f"
		return 0;
#line 609 "AB09IY.f"
	    }
/* Computing MAX */
#line 610 "AB09IY.f"
	    d__1 = work, d__2 = dwork[ku] + (doublereal) (*n);
#line 610 "AB09IY.f"
	    work = max(d__1,d__2);

/*           Partition Sigma = (Sigma1,Sigma2), such that */
/*           Sigma1 <= 0, Sigma2 > 0. */
/*           Partition correspondingly Z = [Z1 Z2]. */

/* Computing MAX */
#line 616 "AB09IY.f"
	    d__2 = abs(dwork[1]), d__3 = (d__1 = dwork[*n], abs(d__1));
#line 616 "AB09IY.f"
	    tol = max(d__2,d__3) * dlamch_("Epsilon", (ftnlen)7);
/*                _ */
/*           Form C = [ sqrt(Sigma2)*Z2' ] */

#line 621 "AB09IY.f"
	    pcbar = 0;
#line 622 "AB09IY.f"
	    i__2 = *n;
#line 622 "AB09IY.f"
	    for (j = 1; j <= i__2; ++j) {
#line 623 "AB09IY.f"
		if (dwork[j] > tol) {
#line 624 "AB09IY.f"
		    d__1 = sqrt(dwork[j]);
#line 624 "AB09IY.f"
		    dscal_(n, &d__1, &r__[j * r_dim1 + 1], &c__1);
#line 625 "AB09IY.f"
		    dcopy_(n, &r__[j * r_dim1 + 1], &c__1, &dwork[ku + pcbar],
			     n);
#line 626 "AB09IY.f"
		    ++pcbar;
#line 627 "AB09IY.f"
		}
#line 628 "AB09IY.f"
/* L40: */
#line 628 "AB09IY.f"
	    }

/*           Solve for the Cholesky factor R of Q, Q = R'*R, */
/*           the continuous-time Lyapunov equation (if DICO = 'C') */
/*                                      _  _ */
/*                   A'*Q + Q*A  +  t^2*C'*C = 0, */

/*           or the discrete-time Lyapunov equation (if DICO = 'D') */
/*                                      _  _ */
/*                   A'*Q*A - Q  +  t^2*C'*C = 0. */

/*           Workspace:  need   N*(N + 6); */
/*                              prefer larger. */

#line 642 "AB09IY.f"
	    ktau = ku + *n * *n;
#line 643 "AB09IY.f"
	    kw = ktau + *n;

#line 645 "AB09IY.f"
	    i__2 = *ldwork - kw + 1;
#line 645 "AB09IY.f"
	    sb03ou_(&discr, &c_false, n, &pcbar, &a[a_offset], lda, &dwork[ku]
		    , n, &dwork[ktau], &r__[r_offset], ldr, &t, &dwork[kw], &
		    i__2, &ierr);
#line 648 "AB09IY.f"
	    if (ierr != 0) {
#line 649 "AB09IY.f"
		*info = 1;
#line 650 "AB09IY.f"
		return 0;
#line 651 "AB09IY.f"
	    }
#line 652 "AB09IY.f"
	    *scaleo *= t;
/* Computing MAX */
#line 653 "AB09IY.f"
	    d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 653 "AB09IY.f"
	    work = max(d__1,d__2);
#line 654 "AB09IY.f"
	}

#line 656 "AB09IY.f"
    } else {

/*        Solve for the Cholesky factor R of Q, Q = R'*R, */
/*        the continuous-time Lyapunov equation (if DICO = 'C') */

/*            A'*Q + Q*A  +  scaleo^2*C'*C = 0, */

/*        or the discrete-time Lyapunov equation (if DICO = 'D') */

/*            A'*Q*A - Q +  scaleo^2*C'*C = 0. */

/*        Workspace:  need   N*(P + 5); */
/*                           prefer larger. */

#line 670 "AB09IY.f"
	ku = 1;
#line 671 "AB09IY.f"
	ktau = ku + *p * *n;
#line 672 "AB09IY.f"
	kw = ktau + *n;

#line 674 "AB09IY.f"
	dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[ku], p, (ftnlen)4);
#line 675 "AB09IY.f"
	i__2 = *ldwork - kw + 1;
#line 675 "AB09IY.f"
	sb03ou_(&discr, &c_false, n, p, &a[a_offset], lda, &dwork[ku], p, &
		dwork[ktau], &r__[r_offset], ldr, scaleo, &dwork[kw], &i__2, &
		ierr);
#line 678 "AB09IY.f"
	if (ierr != 0) {
#line 679 "AB09IY.f"
	    *info = 1;
#line 680 "AB09IY.f"
	    return 0;
#line 681 "AB09IY.f"
	}
/* Computing MAX */
#line 682 "AB09IY.f"
	d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 682 "AB09IY.f"
	work = max(d__1,d__2);
#line 683 "AB09IY.f"
    }

#line 685 "AB09IY.f"
    if (rightw && *mw > 0) {

/*        Build the extended matrices */

/*           Ai = ( A  B*Cw ) ,   Bi = ( B*Dw ) . */
/*                ( 0   Aw  )          (  Bw  ) */

#line 692 "AB09IY.f"
	kaw = 1;
#line 693 "AB09IY.f"
	ku = kaw + nnw * nnw;
#line 694 "AB09IY.f"
	dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[kaw], &nnw, (ftnlen)4)
		;
#line 695 "AB09IY.f"
	dlaset_("Full", nw, n, &c_b16, &c_b16, &dwork[kaw + *n], &nnw, (
		ftnlen)4);
#line 696 "AB09IY.f"
	dgemm_("No-transpose", "No-transpose", n, nw, m, &c_b20, &b[b_offset],
		 ldb, &cw[cw_offset], ldcw, &c_b16, &dwork[kaw + nnw * *n], &
		nnw, (ftnlen)12, (ftnlen)12);
#line 698 "AB09IY.f"
	dlacpy_("Full", nw, nw, &aw[aw_offset], ldaw, &dwork[kaw + nnw * *n + 
		*n], &nnw, (ftnlen)4);

#line 701 "AB09IY.f"
	dgemm_("No-transpose", "No-transpose", n, mw, m, &c_b20, &b[b_offset],
		 ldb, &dw[dw_offset], lddw, &c_b16, &dwork[ku], &nnw, (ftnlen)
		12, (ftnlen)12);
#line 703 "AB09IY.f"
	dlacpy_("Full", nw, mw, &bw[bw_offset], ldbw, &dwork[ku + *n], &nnw, (
		ftnlen)4);

/*        Solve for the Cholesky factor Si of Pi, Pi = Si*Si', */
/*        the continuous-time Lyapunov equation (if DICO = 'C') */

/*            Ai*Pi + Pi*Ai' +  scalec^2*Bi*Bi' = 0, */

/*        or the discrete-time Lyapunov equation (if DICO = 'D') */

/*            Ai*Pi*Ai' - Pi +  scalec^2*Bi*Bi' = 0. */

/*        Workspace:  need   (N+NW)*(N+NW+MAX(N+NW,MW)+5); */
/*                           prefer larger. */

#line 717 "AB09IY.f"
	ktau = ku + nnw * max(nnw,*mw);
#line 718 "AB09IY.f"
	kw = ktau + nnw;

#line 720 "AB09IY.f"
	i__2 = *ldwork - kw + 1;
#line 720 "AB09IY.f"
	sb03ou_(&discr, &c_true, &nnw, mw, &dwork[kaw], &nnw, &dwork[ku], &
		nnw, &dwork[ktau], &dwork[ku], &nnw, scalec, &dwork[kw], &
		i__2, &ierr);

#line 724 "AB09IY.f"
	if (ierr != 0) {
#line 725 "AB09IY.f"
	    *info = 2;
#line 726 "AB09IY.f"
	    return 0;
#line 727 "AB09IY.f"
	}
/* Computing MAX */
#line 728 "AB09IY.f"
	d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 728 "AB09IY.f"
	work = max(d__1,d__2);

/*        Partition Si as Si = ( S11 S12 ) and compute S such that */
/*                             (  0  S22 ) */

/*        S*S' = S11*S11' + (1-ALPHAC**2)*S12*S12'. */

#line 735 "AB09IY.f"
	dlacpy_("Upper", n, n, &dwork[ku], &nnw, &s[s_offset], lds, (ftnlen)5)
		;
#line 736 "AB09IY.f"
	if (*alphac != 0.) {
#line 737 "AB09IY.f"
	    t = sqrt(1. - *alphac * *alphac);
#line 738 "AB09IY.f"
	    i__2 = ku + nnw * (nnw - 1);
#line 738 "AB09IY.f"
	    i__1 = nnw;
#line 738 "AB09IY.f"
	    for (j = ku + nnw * *n; i__1 < 0 ? j >= i__2 : j <= i__2; j += 
		    i__1) {
#line 739 "AB09IY.f"
		dscal_(n, &t, &dwork[j], &c__1);
#line 740 "AB09IY.f"
/* L50: */
#line 740 "AB09IY.f"
	    }
#line 741 "AB09IY.f"
	}
#line 742 "AB09IY.f"
	if (*alphac < 1. && *nw > 0) {
#line 743 "AB09IY.f"
	    ktau = *n * nnw + 1;
#line 744 "AB09IY.f"
	    kw = ktau + *n;
#line 745 "AB09IY.f"
	    mb04nd_("Full", n, &c__0, nw, &s[s_offset], lds, &dwork[ku + nnw *
		     *n], &nnw, dum, &c__1, dum, &c__1, &dwork[ktau], &dwork[
		    kw], (ftnlen)4);

#line 748 "AB09IY.f"
	    i__1 = *n;
#line 748 "AB09IY.f"
	    for (j = 1; j <= i__1; ++j) {
#line 749 "AB09IY.f"
		if (s[j + j * s_dim1] < 0.) {
#line 750 "AB09IY.f"
		    i__2 = j;
#line 750 "AB09IY.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 751 "AB09IY.f"
			s[i__ + j * s_dim1] = -s[i__ + j * s_dim1];
#line 752 "AB09IY.f"
/* L60: */
#line 752 "AB09IY.f"
		    }
#line 753 "AB09IY.f"
		}
#line 754 "AB09IY.f"
/* L70: */
#line 754 "AB09IY.f"
	    }
#line 755 "AB09IY.f"
	}

#line 757 "AB09IY.f"
	if (lsame_(jobc, "E", (ftnlen)1, (ftnlen)1) && *alphac < 1.) {

/*           Form X = -A*(S*S')-(S*S')*A' if DICO = 'C', or */
/*                X = -A*(S*S')*A'+(S*S') if DICO = 'D'. */

#line 762 "AB09IY.f"
	    dlacpy_("Upper", n, n, &s[s_offset], lds, &dwork[ku], n, (ftnlen)
		    5);
#line 763 "AB09IY.f"
	    mb01wd_(dico, "Upper", "Transpose", "Hessenberg", n, &c_b43, &
		    c_b16, &s[s_offset], lds, &dwork[kaw], &nnw, &dwork[ku], 
		    n, &ierr, (ftnlen)1, (ftnlen)5, (ftnlen)9, (ftnlen)10);

/*           Compute the eigendecomposition of X as X = Z*Sigma*Z'. */

#line 769 "AB09IY.f"
	    ku = *n + 1;
#line 770 "AB09IY.f"
	    i__1 = *ldwork - *n;
#line 770 "AB09IY.f"
	    dsyev_("Vectors", "Upper", n, &s[s_offset], lds, &dwork[1], &
		    dwork[ku], &i__1, &ierr, (ftnlen)7, (ftnlen)5);
#line 772 "AB09IY.f"
	    if (ierr > 0) {
#line 773 "AB09IY.f"
		*info = 3;
#line 774 "AB09IY.f"
		return 0;
#line 775 "AB09IY.f"
	    }
/* Computing MAX */
#line 776 "AB09IY.f"
	    d__1 = work, d__2 = dwork[ku] + (doublereal) (*n);
#line 776 "AB09IY.f"
	    work = max(d__1,d__2);

/*           Partition Sigma = (Sigma1,Sigma2), such that */
/*           Sigma1 =< 0, Sigma2 > 0. */
/*           Partition correspondingly Z = [Z1 Z2]. */

/* Computing MAX */
#line 782 "AB09IY.f"
	    d__2 = abs(dwork[1]), d__3 = (d__1 = dwork[*n], abs(d__1));
#line 782 "AB09IY.f"
	    tol = max(d__2,d__3) * dlamch_("Epsilon", (ftnlen)7);
/*                _ */
/*           Form B = [ Z2*sqrt(Sigma2) ] */

#line 787 "AB09IY.f"
	    mbbar = 0;
#line 788 "AB09IY.f"
	    i__ = ku;
#line 789 "AB09IY.f"
	    i__1 = *n;
#line 789 "AB09IY.f"
	    for (j = 1; j <= i__1; ++j) {
#line 790 "AB09IY.f"
		if (dwork[j] > tol) {
#line 791 "AB09IY.f"
		    ++mbbar;
#line 792 "AB09IY.f"
		    d__1 = sqrt(dwork[j]);
#line 792 "AB09IY.f"
		    dscal_(n, &d__1, &s[j * s_dim1 + 1], &c__1);
#line 793 "AB09IY.f"
		    dcopy_(n, &s[j * s_dim1 + 1], &c__1, &dwork[i__], &c__1);
#line 794 "AB09IY.f"
		    i__ += *n;
#line 795 "AB09IY.f"
		}
#line 796 "AB09IY.f"
/* L80: */
#line 796 "AB09IY.f"
	    }

/*           Solve for the Cholesky factor S of P, P = S*S', */
/*           the continuous-time Lyapunov equation (if DICO = 'C') */
/*                                      _ _ */
/*                   A*P + P*A'  +  t^2*B*B' = 0, */

/*           or the discrete-time Lyapunov equation (if DICO = 'D') */
/*                                      _ _ */
/*                   A*P*A' - P  +  t^2*B*B' = 0. */

/*           Workspace:  need   maximum N*(N + 6); */
/*                              prefer larger. */

#line 810 "AB09IY.f"
	    ktau = ku + mbbar * *n;
#line 811 "AB09IY.f"
	    kw = ktau + *n;

#line 813 "AB09IY.f"
	    i__1 = *ldwork - kw + 1;
#line 813 "AB09IY.f"
	    sb03ou_(&discr, &c_true, n, &mbbar, &a[a_offset], lda, &dwork[ku],
		     n, &dwork[ktau], &s[s_offset], lds, &t, &dwork[kw], &
		    i__1, &ierr);
#line 816 "AB09IY.f"
	    if (ierr != 0) {
#line 817 "AB09IY.f"
		*info = 2;
#line 818 "AB09IY.f"
		return 0;
#line 819 "AB09IY.f"
	    }
#line 820 "AB09IY.f"
	    *scalec *= t;
/* Computing MAX */
#line 821 "AB09IY.f"
	    d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 821 "AB09IY.f"
	    work = max(d__1,d__2);
#line 822 "AB09IY.f"
	}

#line 824 "AB09IY.f"
    } else {

/*        Solve for the Cholesky factor S of P, P = S*S', */
/*        the continuous-time Lyapunov equation (if DICO = 'C') */

/*            A*P + P*A' +  scalec^2*B*B' = 0, */

/*        or the discrete-time Lyapunov equation (if DICO = 'D') */

/*            A*P*A' - P +  scalec^2*B*B' = 0. */

/*        Workspace:  need   N*(M+5); */
/*                           prefer larger. */

#line 838 "AB09IY.f"
	ku = 1;
#line 839 "AB09IY.f"
	ktau = ku + *n * *m;
#line 840 "AB09IY.f"
	kw = ktau + *n;

#line 842 "AB09IY.f"
	dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[ku], n, (ftnlen)4);
#line 843 "AB09IY.f"
	i__1 = *ldwork - kw + 1;
#line 843 "AB09IY.f"
	sb03ou_(&discr, &c_true, n, m, &a[a_offset], lda, &dwork[ku], n, &
		dwork[ktau], &s[s_offset], lds, scalec, &dwork[kw], &i__1, &
		ierr);
#line 846 "AB09IY.f"
	if (ierr != 0) {
#line 847 "AB09IY.f"
	    *info = 2;
#line 848 "AB09IY.f"
	    return 0;
#line 849 "AB09IY.f"
	}
/* Computing MAX */
#line 850 "AB09IY.f"
	d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 850 "AB09IY.f"
	work = max(d__1,d__2);
#line 851 "AB09IY.f"
    }

/*     Save optimal workspace. */

#line 855 "AB09IY.f"
    dwork[1] = work;

#line 857 "AB09IY.f"
    return 0;
/* *** Last line of AB09IY *** */
} /* ab09iy_ */

