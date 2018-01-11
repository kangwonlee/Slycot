#line 1 "AB09JV.f"
/* AB09JV.f -- translated by f2c (version 20100827).
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

#line 1 "AB09JV.f"
/* Table of constant values */

static doublereal c_b17 = -1.;
static doublereal c_b18 = 0.;
static integer c_n1 = -1;
static doublereal c_b24 = 1.;
static integer c__1 = 1;
static integer c__0 = 0;

/* Subroutine */ int ab09jv_(char *job, char *dico, char *jobev, char *stbchk,
	 integer *n, integer *m, integer *p, integer *nv, integer *pv, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *d__, integer *ldd, doublereal *av, 
	integer *ldav, doublereal *ev, integer *ldev, doublereal *bv, integer 
	*ldbv, doublereal *cv, integer *ldcv, doublereal *dv, integer *lddv, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen job_len, ftnlen dico_len, ftnlen jobev_len, ftnlen stbchk_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, av_dim1, av_offset, b_dim1, b_offset, bv_dim1, 
	    bv_offset, c_dim1, c_offset, cv_dim1, cv_offset, d_dim1, d_offset,
	     dv_dim1, dv_offset, ev_dim1, ev_offset, i__1, i__2, i__3, i__4, 
	    i__5, i__6;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, ia, kb, kc, ke, kf, kq, kw, lw, kz;
    static doublereal dif;
    static integer kai, kar, ldw, sdim, ierr, ldwn;
    static doublereal work, alpha, scale;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     dgges_(char *, char *, char *, L_fp, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, logical *, integer *, ftnlen,
	     ftnlen, ftnlen), ab09jx_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01wd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);
    static logical discr, conjs;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sb04py_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen);
    static logical bwork[1];
    static char stdom[1];
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern logical delctg_();
    static logical stabck;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static doublereal tolinf;
    extern /* Subroutine */ int dtgsyl_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, ftnlen);
    static logical unitev;
    static char evtype[1];
    extern /* Subroutine */ int dtrsyl_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen);


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

/*     To construct a state-space representation (A,BS,CS,DS) of the */
/*     projection of V*G or conj(V)*G containing the poles of G, from the */
/*     state-space representations (A,B,C,D) and (AV-lambda*EV,BV,CV,DV), */
/*     of the transfer-function matrices G and V, respectively. */
/*     G is assumed to be a stable transfer-function matrix and */
/*     the state matrix A must be in a real Schur form. */
/*     When computing the stable projection of V*G, it is assumed */
/*     that G and V have completely distinct poles. */
/*     When computing the stable projection of conj(V)*G, it is assumed */
/*     that G and conj(V) have completely distinct poles. */

/*     Note: For a transfer-function matrix G, conj(G) denotes the */
/*     conjugate of G given by G'(-s) for a continuous-time system or */
/*     G'(1/z) for a discrete-time system. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the projection to be computed as follows: */
/*             = 'V':  compute the projection of V*G containing */
/*                     the poles of G; */
/*             = 'C':  compute the projection of conj(V)*G containing */
/*                     the poles of G. */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the systems as follows: */
/*             = 'C':  G and V are continuous-time systems; */
/*             = 'D':  G and V are discrete-time systems. */

/*     JOBEV   CHARACTER*1 */
/*             Specifies whether EV is a general square or an identity */
/*             matrix as follows: */
/*             = 'G':  EV is a general square matrix; */
/*             = 'I':  EV is the identity matrix. */

/*     STBCHK  CHARACTER*1 */
/*             Specifies whether stability/antistability of V is to be */
/*             checked as follows: */
/*             = 'C':  check stability if JOB = 'C' or antistability if */
/*                     JOB = 'V'; */
/*             = 'N':  do not check stability or antistability. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The dimension of the state vector of the system with */
/*             the transfer-function matrix G.  N >= 0. */

/*     M       (input) INTEGER */
/*             The dimension of the input vector of the system with */
/*             the transfer-function matrix G.  M >= 0. */

/*     P       (input) INTEGER */
/*             The dimension of the output vector of the system with the */
/*             transfer-function matrix G, and also the dimension of */
/*             the input vector if JOB = 'V', or of the output vector */
/*             if JOB = 'C', of the system with the transfer-function */
/*             matrix V.  P >= 0. */

/*     NV      (input) INTEGER */
/*             The dimension of the state vector of the system with */
/*             the transfer-function matrix V.  NV >= 0. */

/*     PV      (input) INTEGER */
/*             The dimension of the output vector, if JOB = 'V', or */
/*             of the input vector, if JOB = 'C', of the system with */
/*             the transfer-function matrix V.  PV >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             state matrix A of the system with the transfer-function */
/*             matrix G in a real Schur form. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain */
/*             the input/state matrix B of the system with the */
/*             transfer-function matrix G. The matrix BS is equal to B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the output matrix C of the system with the */
/*             transfer-function matrix G. */
/*             On exit, if INFO = 0, the leading PV-by-N part of this */
/*             array contains the output matrix CS of the projection of */
/*             V*G, if JOB = 'V', or of conj(V)*G, if JOB = 'C'. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= MAX(1,P,PV). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the feedthrough matrix D of the system with the */
/*             transfer-function matrix G. */
/*             On exit, if INFO = 0, the leading PV-by-M part of */
/*             this array contains the feedthrough matrix DS of the */
/*             projection of V*G, if JOB = 'V', or of conj(V)*G, */
/*             if JOB = 'C'. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= MAX(1,P,PV). */

/*     AV      (input/output) DOUBLE PRECISION array, dimension (LDAV,NV) */
/*             On entry, the leading NV-by-NV part of this array must */
/*             contain the state matrix AV of the system with the */
/*             transfer-function matrix V. */
/*             On exit, if INFO = 0, the leading NV-by-NV part of this */
/*             array contains a condensed matrix as follows: */
/*             if JOBEV = 'I', it contains the real Schur form of AV; */
/*             if JOBEV = 'G' and JOB = 'V', it contains a quasi-upper */
/*             triangular matrix representing the real Schur matrix */
/*             in the real generalized Schur form of the pair (AV,EV); */
/*             if JOBEV = 'G', JOB = 'C' and DICO = 'C', it contains a */
/*             quasi-upper triangular matrix corresponding to the */
/*             generalized real Schur form of the pair (AV',EV'); */
/*             if JOBEV = 'G', JOB = 'C' and DICO = 'D', it contains an */
/*             upper triangular matrix corresponding to the generalized */
/*             real Schur form of the pair (EV',AV'). */

/*     LDAV    INTEGER */
/*             The leading dimension of the array AV.  LDAV >= MAX(1,NV). */

/*     EV      (input/output) DOUBLE PRECISION array, dimension (LDEV,NV) */
/*             On entry, if JOBEV = 'G', the leading NV-by-NV part of */
/*             this array must contain the descriptor matrix EV of the */
/*             system with the transfer-function matrix V. */
/*             If JOBEV = 'I', EV is assumed to be an identity matrix */
/*             and is not referenced. */
/*             On exit, if INFO = 0 and JOBEV = 'G', the leading NV-by-NV */
/*             part of this array contains a condensed matrix as follows: */
/*             if JOB = 'V', it contains an upper triangular matrix */
/*             corresponding to the real generalized Schur form of the */
/*             pair (AV,EV); */
/*             if JOB = 'C' and DICO = 'C', it contains an upper */
/*             triangular matrix corresponding to the generalized real */
/*             Schur form of the pair (AV',EV'); */
/*             if JOB = 'C' and DICO = 'D', it contains a quasi-upper */
/*             triangular matrix corresponding to the generalized */
/*             real Schur form of the pair (EV',AV'). */

/*     LDEV    INTEGER */
/*             The leading dimension of the array EV. */
/*             LDEV >= MAX(1,NV), if JOBEV = 'G'; */
/*             LDEV >= 1,         if JOBEV = 'I'. */

/*     BV      (input/output) DOUBLE PRECISION array, */
/*             dimension (LDBV,MBV), where MBV = P, if JOB = 'V', and */
/*             MBV = PV, if JOB = 'C'. */
/*             On entry, the leading NV-by-MBV part of this array must */
/*             contain the input matrix BV of the system with the */
/*             transfer-function matrix V. */
/*             On exit, if INFO = 0, the leading NV-by-MBV part of this */
/*             array contains Q'*BV, where Q is the orthogonal matrix */
/*             that reduces AV to the real Schur form or the left */
/*             orthogonal matrix used to reduce the pair (AV,EV), */
/*             (AV',EV') or (EV',AV') to the generalized real Schur form. */

/*     LDBV    INTEGER */
/*             The leading dimension of the array BV.  LDBV >= MAX(1,NV). */

/*     CV      (input/output) DOUBLE PRECISION array, dimension (LDCV,NV) */
/*             On entry, the leading PCV-by-NV part of this array must */
/*             contain the output matrix CV of the system with the */
/*             transfer-function matrix V, where PCV = PV, if JOB = 'V', */
/*             or PCV = P, if JOB = 'C'. */
/*             On exit, if INFO = 0, the leading PCV-by-NV part of this */
/*             array contains CV*Q, where Q is the orthogonal matrix that */
/*             reduces AV to the real Schur form, or CV*Z, where Z is the */
/*             right orthogonal matrix used to reduce the pair (AV,EV), */
/*             (AV',EV') or (EV',AV') to the generalized real Schur form. */

/*     LDCV    INTEGER */
/*             The leading dimension of the array CV. */
/*             LDCV >= MAX(1,PV) if JOB = 'V'; */
/*             LDCV >= MAX(1,P)  if JOB = 'C'. */

/*     DV      (input) DOUBLE PRECISION array, */
/*             dimension (LDDV,MBV), where MBV = P, if JOB = 'V', and */
/*             MBV = PV, if JOB = 'C'. */
/*             The leading PCV-by-MBV part of this array must contain */
/*             the feedthrough matrix DV of the system with the */
/*             transfer-function matrix V, where PCV = PV, if JOB = 'V', */
/*             or PCV = P, if JOB = 'C'. */

/*     LDDV    INTEGER */
/*             The leading dimension of the array DV. */
/*             LDDV >= MAX(1,PV) if JOB = 'V'; */
/*             LDDV >= MAX(1,P)  if JOB = 'C'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK =   0,    if JOBEV = 'I'; */
/*             LIWORK = NV+N+6, if JOBEV = 'G'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= LW1, if JOBEV = 'I', */
/*             LDWORK >= LW2, if JOBEV = 'G', where */
/*               LW1 = MAX( 1, NV*(NV+5), NV*N + MAX( a, PV*N, PV*M ) ) */
/*                     a = 0,    if DICO = 'C' or  JOB = 'V', */
/*                     a = 2*NV, if DICO = 'D' and JOB = 'C'; */
/*               LW2 = MAX( 2*NV*NV + MAX( 11*NV+16, P*NV, PV*NV ), */
/*                          NV*N + MAX( NV*N+N*N, PV*N, PV*M ) ). */
/*             For good performance, LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             =  0:  successful exit; */
/*             <  0:  if INFO = -i, the i-th argument had an illegal */
/*                    value; */
/*             =  1:  the reduction of the pair (AV,EV) to the real */
/*                    generalized Schur form failed (JOBEV = 'G'), */
/*                    or the reduction of the matrix AV to the real */
/*                    Schur form failed (JOBEV = 'I); */
/*             =  2:  the solution of the Sylvester equation failed */
/*                    because the matrix A and the pencil AV-lambda*EV */
/*                    have common eigenvalues (if JOB = 'V'), or the */
/*                    pencil -AV-lambda*EV and A have common eigenvalues */
/*                    (if JOB = 'C' and DICO = 'C'), or the pencil */
/*                    AV-lambda*EV has an eigenvalue which is the */
/*                    reciprocal of one of eigenvalues of A */
/*                    (if JOB = 'C' and DICO = 'D'); */
/*             =  3:  the solution of the Sylvester equation failed */
/*                    because the matrices A and AV have common */
/*                    eigenvalues (if JOB = 'V'), or the matrices A */
/*                    and -AV have common eigenvalues (if JOB = 'C' and */
/*                    DICO = 'C'), or the matrix A has an eigenvalue */
/*                    which is the reciprocal of one of eigenvalues of AV */
/*                    (if JOB = 'C' and DICO = 'D'); */
/*             =  4:  JOB = 'V' and the pair (AV,EV) has not completely */
/*                    unstable generalized eigenvalues, or JOB = 'C' and */
/*                    the pair (AV,EV) has not completely stable */
/*                    generalized eigenvalues. */

/*     METHOD */

/*     If JOB = 'V', the matrices of the stable projection of V*G are */
/*     computed as */

/*       BS = B,  CS = CV*X + DV*C,  DS = DV*D, */

/*     where X satisfies the generalized Sylvester equation */

/*       AV*X - EV*X*A + BV*C = 0. */

/*     If JOB = 'C', the matrices of the stable projection of conj(V)*G */
/*     are computed using the following formulas: */

/*     - for a continuous-time system, the matrices BS, CS and DS of */
/*       the stable projection are computed as */

/*         BS = B,  CS = BV'*X + DV'*C,  DS = DV'*D, */

/*       where X satisfies the generalized Sylvester equation */

/*         AV'*X + EV'*X*A + CV'*C = 0. */

/*     - for a discrete-time system, the matrices BS, CS and DS of */
/*       the stable projection are computed as */

/*         BS = B,  CS = BV'*X*A + DV'*C,  DS = DV'*D + BV'*X*B, */

/*       where X satisfies the generalized Sylvester equation */

/*         EV'*X - AV'*X*A = CV'*C. */

/*     REFERENCES */

/*     [1] Varga, A. */
/*         Efficient and numerically reliable implementation of the */
/*         frequency-weighted Hankel-norm approximation model reduction */
/*         approach. */
/*         Proc. 2001 ECC, Porto, Portugal, 2001. */

/*     [2] Zhou, K. */
/*         Frequency-weighted H-infinity norm and optimal Hankel norm */
/*         model reduction. */
/*         IEEE Trans. Autom. Control, vol. 40, pp. 1687-1699, 1995. */

/*     NUMERICAL ASPECTS */

/*     The implemented methods rely on numerically stable algorithms. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, July 2000. */
/*     D. Sima, University of Bucharest, March 2001. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001. */

/*     REVISIONS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2001. */
/*     V. Sima, Research Institute for Informatics, Bucharest, June 2001. */
/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, Nov. 2003. */

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

#line 373 "AB09JV.f"
    /* Parameter adjustments */
#line 373 "AB09JV.f"
    a_dim1 = *lda;
#line 373 "AB09JV.f"
    a_offset = 1 + a_dim1;
#line 373 "AB09JV.f"
    a -= a_offset;
#line 373 "AB09JV.f"
    b_dim1 = *ldb;
#line 373 "AB09JV.f"
    b_offset = 1 + b_dim1;
#line 373 "AB09JV.f"
    b -= b_offset;
#line 373 "AB09JV.f"
    c_dim1 = *ldc;
#line 373 "AB09JV.f"
    c_offset = 1 + c_dim1;
#line 373 "AB09JV.f"
    c__ -= c_offset;
#line 373 "AB09JV.f"
    d_dim1 = *ldd;
#line 373 "AB09JV.f"
    d_offset = 1 + d_dim1;
#line 373 "AB09JV.f"
    d__ -= d_offset;
#line 373 "AB09JV.f"
    av_dim1 = *ldav;
#line 373 "AB09JV.f"
    av_offset = 1 + av_dim1;
#line 373 "AB09JV.f"
    av -= av_offset;
#line 373 "AB09JV.f"
    ev_dim1 = *ldev;
#line 373 "AB09JV.f"
    ev_offset = 1 + ev_dim1;
#line 373 "AB09JV.f"
    ev -= ev_offset;
#line 373 "AB09JV.f"
    bv_dim1 = *ldbv;
#line 373 "AB09JV.f"
    bv_offset = 1 + bv_dim1;
#line 373 "AB09JV.f"
    bv -= bv_offset;
#line 373 "AB09JV.f"
    cv_dim1 = *ldcv;
#line 373 "AB09JV.f"
    cv_offset = 1 + cv_dim1;
#line 373 "AB09JV.f"
    cv -= cv_offset;
#line 373 "AB09JV.f"
    dv_dim1 = *lddv;
#line 373 "AB09JV.f"
    dv_offset = 1 + dv_dim1;
#line 373 "AB09JV.f"
    dv -= dv_offset;
#line 373 "AB09JV.f"
    --iwork;
#line 373 "AB09JV.f"
    --dwork;
#line 373 "AB09JV.f"

#line 373 "AB09JV.f"
    /* Function Body */
#line 373 "AB09JV.f"
    conjs = lsame_(job, "C", (ftnlen)1, (ftnlen)1);
#line 374 "AB09JV.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 375 "AB09JV.f"
    unitev = lsame_(jobev, "I", (ftnlen)1, (ftnlen)1);
#line 376 "AB09JV.f"
    stabck = lsame_(stbchk, "C", (ftnlen)1, (ftnlen)1);

#line 378 "AB09JV.f"
    *info = 0;
#line 379 "AB09JV.f"
    if (unitev) {
#line 380 "AB09JV.f"
	if (discr && conjs) {
#line 381 "AB09JV.f"
	    ia = *nv << 1;
#line 382 "AB09JV.f"
	} else {
#line 383 "AB09JV.f"
	    ia = 0;
#line 384 "AB09JV.f"
	}
/* Computing MAX */
/* Computing MAX */
#line 385 "AB09JV.f"
	i__3 = ia, i__4 = *pv * *n, i__3 = max(i__3,i__4), i__4 = *pv * *m;
#line 385 "AB09JV.f"
	i__1 = 1, i__2 = *nv * (*nv + 5), i__1 = max(i__1,i__2), i__2 = *nv * 
		*n + max(i__3,i__4);
#line 385 "AB09JV.f"
	lw = max(i__1,i__2);
#line 386 "AB09JV.f"
    } else {
/* Computing MAX */
/* Computing MAX */
#line 387 "AB09JV.f"
	i__3 = *nv * 11 + 16, i__4 = *p * *nv, i__3 = max(i__3,i__4), i__4 = *
		pv * *nv;
/* Computing MAX */
#line 387 "AB09JV.f"
	i__5 = *nv * *n + *n * *n, i__6 = *pv * *n, i__5 = max(i__5,i__6), 
		i__6 = *pv * *m;
#line 387 "AB09JV.f"
	i__1 = (*nv << 1) * *nv + max(i__3,i__4), i__2 = *nv * *n + max(i__5,
		i__6);
#line 387 "AB09JV.f"
	lw = max(i__1,i__2);
#line 389 "AB09JV.f"
    }

/*     Test the input scalar arguments. */

#line 393 "AB09JV.f"
    ldwn = max(1,*n);
#line 394 "AB09JV.f"
    ldw = max(1,*nv);
#line 395 "AB09JV.f"
    if (! (lsame_(job, "V", (ftnlen)1, (ftnlen)1) || conjs)) {
#line 396 "AB09JV.f"
	*info = -1;
#line 397 "AB09JV.f"
    } else if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 398 "AB09JV.f"
	*info = -2;
#line 399 "AB09JV.f"
    } else if (! (lsame_(jobev, "G", (ftnlen)1, (ftnlen)1) || unitev)) {
#line 400 "AB09JV.f"
	*info = -3;
#line 401 "AB09JV.f"
    } else if (! (lsame_(stbchk, "N", (ftnlen)1, (ftnlen)1) || stabck)) {
#line 402 "AB09JV.f"
	*info = -4;
#line 403 "AB09JV.f"
    } else if (*n < 0) {
#line 404 "AB09JV.f"
	*info = -5;
#line 405 "AB09JV.f"
    } else if (*m < 0) {
#line 406 "AB09JV.f"
	*info = -6;
#line 407 "AB09JV.f"
    } else if (*p < 0) {
#line 408 "AB09JV.f"
	*info = -7;
#line 409 "AB09JV.f"
    } else if (*nv < 0) {
#line 410 "AB09JV.f"
	*info = -8;
#line 411 "AB09JV.f"
    } else if (*pv < 0) {
#line 412 "AB09JV.f"
	*info = -9;
#line 413 "AB09JV.f"
    } else if (*lda < ldwn) {
#line 414 "AB09JV.f"
	*info = -11;
#line 415 "AB09JV.f"
    } else if (*ldb < ldwn) {
#line 416 "AB09JV.f"
	*info = -13;
#line 417 "AB09JV.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 417 "AB09JV.f"
	i__1 = max(1,*p);
#line 417 "AB09JV.f"
	if (*ldc < max(i__1,*pv)) {
#line 418 "AB09JV.f"
	    *info = -15;
#line 419 "AB09JV.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 419 "AB09JV.f"
	    i__1 = max(1,*p);
#line 419 "AB09JV.f"
	    if (*ldd < max(i__1,*pv)) {
#line 420 "AB09JV.f"
		*info = -17;
#line 421 "AB09JV.f"
	    } else if (*ldav < ldw) {
#line 422 "AB09JV.f"
		*info = -19;
#line 423 "AB09JV.f"
	    } else if (*ldev < 1 || ! unitev && *ldev < *nv) {
#line 424 "AB09JV.f"
		*info = -21;
#line 425 "AB09JV.f"
	    } else if (*ldbv < ldw) {
#line 426 "AB09JV.f"
		*info = -23;
#line 427 "AB09JV.f"
	    } else if (! conjs && *ldcv < max(1,*pv) || conjs && *ldcv < max(
		    1,*p)) {
#line 429 "AB09JV.f"
		*info = -25;
#line 430 "AB09JV.f"
	    } else if (! conjs && *lddv < max(1,*pv) || conjs && *lddv < max(
		    1,*p)) {
#line 432 "AB09JV.f"
		*info = -27;
#line 433 "AB09JV.f"
	    } else if (*ldwork < lw) {
#line 434 "AB09JV.f"
		*info = -30;
#line 435 "AB09JV.f"
	    }
#line 435 "AB09JV.f"
	}
#line 435 "AB09JV.f"
    }

#line 437 "AB09JV.f"
    if (*info != 0) {

/*        Error return. */

#line 441 "AB09JV.f"
	i__1 = -(*info);
#line 441 "AB09JV.f"
	xerbla_("AB09JV", &i__1, (ftnlen)6);
#line 442 "AB09JV.f"
	return 0;
#line 443 "AB09JV.f"
    }

/*     Quick return if possible. */

#line 447 "AB09JV.f"
    if (*p == 0 || *pv == 0) {
#line 448 "AB09JV.f"
	dwork[1] = 1.;
#line 449 "AB09JV.f"
	return 0;
#line 450 "AB09JV.f"
    }

/*     Set options for stability/antistability checking. */

#line 454 "AB09JV.f"
    if (discr) {
#line 455 "AB09JV.f"
	alpha = 1.;
#line 456 "AB09JV.f"
    } else {
#line 457 "AB09JV.f"
	alpha = 0.;
#line 458 "AB09JV.f"
    }

#line 460 "AB09JV.f"
    work = 1.;
#line 461 "AB09JV.f"
    tolinf = dlamch_("Epsilon", (ftnlen)7);

#line 463 "AB09JV.f"
    if (unitev) {

/*        EV is the identity matrix. */

#line 467 "AB09JV.f"
	if (*nv > 0) {

/*           Reduce AV to the real Schur form using an orthogonal */
/*           similarity transformation AV <- Q'*AV*Q and apply the */
/*           transformation to BV and CV: BV <- Q'*BV and CV <- CV*Q. */

/*           Workspace needed:  NV*(NV+5); */
/*                              prefer larger. */

#line 476 "AB09JV.f"
	    kw = *nv * (*nv + 2) + 1;
#line 477 "AB09JV.f"
	    if (conjs) {
#line 478 "AB09JV.f"
		*(unsigned char *)stdom = 'S';
#line 479 "AB09JV.f"
		alpha += sqrt(tolinf);
#line 480 "AB09JV.f"
		i__1 = *ldwork - kw + 1;
#line 480 "AB09JV.f"
		tb01wd_(nv, pv, p, &av[av_offset], ldav, &bv[bv_offset], ldbv,
			 &cv[cv_offset], ldcv, &dwork[(*nv << 1) + 1], nv, &
			dwork[1], &dwork[*nv + 1], &dwork[kw], &i__1, &ierr);
#line 483 "AB09JV.f"
	    } else {
#line 484 "AB09JV.f"
		*(unsigned char *)stdom = 'U';
#line 485 "AB09JV.f"
		alpha -= sqrt(tolinf);
#line 486 "AB09JV.f"
		i__1 = *ldwork - kw + 1;
#line 486 "AB09JV.f"
		tb01wd_(nv, p, pv, &av[av_offset], ldav, &bv[bv_offset], ldbv,
			 &cv[cv_offset], ldcv, &dwork[(*nv << 1) + 1], nv, &
			dwork[1], &dwork[*nv + 1], &dwork[kw], &i__1, &ierr);
#line 489 "AB09JV.f"
	    }
#line 490 "AB09JV.f"
	    if (ierr != 0) {
#line 491 "AB09JV.f"
		*info = 1;
#line 492 "AB09JV.f"
		return 0;
#line 493 "AB09JV.f"
	    }
#line 494 "AB09JV.f"
	    if (stabck) {

/*              Check stability/antistability of eigenvalues of AV. */

#line 498 "AB09JV.f"
		ab09jx_(dico, stdom, "S", nv, &alpha, &dwork[1], &dwork[*nv + 
			1], &dwork[1], &tolinf, &ierr, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 500 "AB09JV.f"
		if (ierr != 0) {
#line 501 "AB09JV.f"
		    *info = 4;
#line 502 "AB09JV.f"
		    return 0;
#line 503 "AB09JV.f"
		}
#line 504 "AB09JV.f"
	    }

/* Computing MAX */
#line 506 "AB09JV.f"
	    d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 506 "AB09JV.f"
	    work = max(d__1,d__2);

#line 508 "AB09JV.f"
	}

#line 510 "AB09JV.f"
	kw = *nv * *n + 1;
#line 511 "AB09JV.f"
	if (conjs) {

/*           Compute the projection of conj(V)*G. */

/*           Total workspace needed:  NV*N + MAX( a, PV*N, PV*M ), where */
/*                                    a = 0,    if DICO = 'C', */
/*                                    a = 2*NV, if DICO = 'D'. */

/*           Compute -CV'*C. */
/*           Workspace needed: NV*N. */

#line 522 "AB09JV.f"
	    dgemm_("T", "N", nv, n, p, &c_b17, &cv[cv_offset], ldcv, &c__[
		    c_offset], ldc, &c_b18, &dwork[1], &ldw, (ftnlen)1, (
		    ftnlen)1);

#line 525 "AB09JV.f"
	    if (discr) {

/*              Compute X and SCALE satisfying */

/*              AV'*X*A - X = -SCALE*CV'*C. */

/*              Additional workspace needed: 2*NV. */

#line 533 "AB09JV.f"
		sb04py_("T", "N", &c_n1, nv, n, &av[av_offset], ldav, &a[
			a_offset], lda, &dwork[1], &ldw, &scale, &dwork[kw], &
			ierr, (ftnlen)1, (ftnlen)1);
#line 535 "AB09JV.f"
		if (ierr != 0) {
#line 536 "AB09JV.f"
		    *info = 3;
#line 537 "AB09JV.f"
		    return 0;
#line 538 "AB09JV.f"
		}

/*              Construct CS = DV'*C + BV'*X*A/SCALE, */
/*                        DS = DV'*D + BV'*X*B/SCALE. */

/*              Additional workspace needed: MAX( PV*N, PV*M ). */

/*              C <- DV'*C. */

#line 547 "AB09JV.f"
		dgemm_("T", "N", pv, n, p, &c_b24, &dv[dv_offset], lddv, &c__[
			c_offset], ldc, &c_b18, &dwork[kw], pv, (ftnlen)1, (
			ftnlen)1);
#line 549 "AB09JV.f"
		dlacpy_("Full", pv, n, &dwork[kw], pv, &c__[c_offset], ldc, (
			ftnlen)4);

/*              D <- DV'*D. */

#line 553 "AB09JV.f"
		dgemm_("T", "N", pv, m, p, &c_b24, &dv[dv_offset], lddv, &d__[
			d_offset], ldd, &c_b18, &dwork[kw], pv, (ftnlen)1, (
			ftnlen)1);
#line 555 "AB09JV.f"
		dlacpy_("Full", pv, m, &dwork[kw], pv, &d__[d_offset], ldd, (
			ftnlen)4);

/*              C <- C + BV'*X*A/SCALE. */

#line 559 "AB09JV.f"
		d__1 = 1. / scale;
#line 559 "AB09JV.f"
		dgemm_("T", "N", pv, n, nv, &d__1, &bv[bv_offset], ldbv, &
			dwork[1], &ldw, &c_b18, &dwork[kw], pv, (ftnlen)1, (
			ftnlen)1);
#line 561 "AB09JV.f"
		dgemm_("N", "N", pv, n, n, &c_b24, &dwork[kw], pv, &a[
			a_offset], lda, &c_b24, &c__[c_offset], ldc, (ftnlen)
			1, (ftnlen)1);

/*              D <- D + BV'*X*B/SCALE. */

#line 566 "AB09JV.f"
		dgemm_("N", "N", pv, m, n, &c_b24, &dwork[kw], pv, &b[
			b_offset], ldb, &c_b24, &d__[d_offset], ldd, (ftnlen)
			1, (ftnlen)1);
#line 568 "AB09JV.f"
	    } else {

/*              Compute X and SCALE satisfying */

/*              AV'*X + X*A + SCALE*CV'*C = 0. */

#line 574 "AB09JV.f"
		if (*n > 0) {
#line 575 "AB09JV.f"
		    dtrsyl_("T", "N", &c__1, nv, n, &av[av_offset], ldav, &a[
			    a_offset], lda, &dwork[1], &ldw, &scale, &ierr, (
			    ftnlen)1, (ftnlen)1);
#line 577 "AB09JV.f"
		    if (ierr != 0) {
#line 578 "AB09JV.f"
			*info = 3;
#line 579 "AB09JV.f"
			return 0;
#line 580 "AB09JV.f"
		    }
#line 581 "AB09JV.f"
		}

/*              Construct CS = DV'*C + BV'*X/SCALE, */
/*                        DS = DV'*D. */
/*              Additional workspace needed: MAX( PV*N, PV*M ). */

/*              Construct C <- DV'*C + BV'*X/SCALE. */

#line 589 "AB09JV.f"
		dgemm_("T", "N", pv, n, p, &c_b24, &dv[dv_offset], lddv, &c__[
			c_offset], ldc, &c_b18, &dwork[kw], pv, (ftnlen)1, (
			ftnlen)1);
#line 591 "AB09JV.f"
		dlacpy_("Full", pv, n, &dwork[kw], pv, &c__[c_offset], ldc, (
			ftnlen)4);
#line 592 "AB09JV.f"
		d__1 = 1. / scale;
#line 592 "AB09JV.f"
		dgemm_("T", "N", pv, n, nv, &d__1, &bv[bv_offset], ldbv, &
			dwork[1], &ldw, &c_b24, &c__[c_offset], ldc, (ftnlen)
			1, (ftnlen)1);

/*              Construct D <- DV'*D. */

#line 597 "AB09JV.f"
		dgemm_("T", "N", pv, m, p, &c_b24, &dv[dv_offset], lddv, &d__[
			d_offset], ldd, &c_b18, &dwork[kw], pv, (ftnlen)1, (
			ftnlen)1);
#line 599 "AB09JV.f"
		dlacpy_("Full", pv, m, &dwork[kw], pv, &d__[d_offset], ldd, (
			ftnlen)4);
#line 600 "AB09JV.f"
	    }
#line 601 "AB09JV.f"
	} else {

/*           Compute the projection of V*G. */

/*           Total workspace needed:  NV*N + MAX( PV*N, PV*M ). */

/*           Compute -BV*C. */
/*           Workspace needed: NV*N. */

#line 610 "AB09JV.f"
	    dgemm_("N", "N", nv, n, p, &c_b17, &bv[bv_offset], ldbv, &c__[
		    c_offset], ldc, &c_b18, &dwork[1], &ldw, (ftnlen)1, (
		    ftnlen)1);

/*           Compute X and SCALE satisfying */

/*           AV*X - X*A + SCALE*BV*C = 0. */

#line 617 "AB09JV.f"
	    if (*n > 0) {
#line 618 "AB09JV.f"
		dtrsyl_("N", "N", &c_n1, nv, n, &av[av_offset], ldav, &a[
			a_offset], lda, &dwork[1], &ldw, &scale, &ierr, (
			ftnlen)1, (ftnlen)1);
#line 620 "AB09JV.f"
		if (ierr != 0) {
#line 621 "AB09JV.f"
		    *info = 3;
#line 622 "AB09JV.f"
		    return 0;
#line 623 "AB09JV.f"
		}
#line 624 "AB09JV.f"
	    }

/*           Construct CS = DV*C + CV*X/SCALE, */
/*                     DS = DV*D. */
/*           Additional workspace needed: MAX( PV*N, PV*M ). */

/*           Construct C <- DV*C + CV*X/SCALE. */

#line 632 "AB09JV.f"
	    dgemm_("N", "N", pv, n, p, &c_b24, &dv[dv_offset], lddv, &c__[
		    c_offset], ldc, &c_b18, &dwork[kw], pv, (ftnlen)1, (
		    ftnlen)1);
#line 634 "AB09JV.f"
	    dlacpy_("Full", pv, n, &dwork[kw], pv, &c__[c_offset], ldc, (
		    ftnlen)4);
#line 635 "AB09JV.f"
	    d__1 = 1. / scale;
#line 635 "AB09JV.f"
	    dgemm_("N", "N", pv, n, nv, &d__1, &cv[cv_offset], ldcv, &dwork[1]
		    , &ldw, &c_b24, &c__[c_offset], ldc, (ftnlen)1, (ftnlen)1)
		    ;

/*           Construct D <- DV*D. */

#line 640 "AB09JV.f"
	    dgemm_("N", "N", pv, m, p, &c_b24, &dv[dv_offset], lddv, &d__[
		    d_offset], ldd, &c_b18, &dwork[kw], pv, (ftnlen)1, (
		    ftnlen)1);
#line 642 "AB09JV.f"
	    dlacpy_("Full", pv, m, &dwork[kw], pv, &d__[d_offset], ldd, (
		    ftnlen)4);
#line 643 "AB09JV.f"
	}
#line 644 "AB09JV.f"
    } else {

/*        EV is a general matrix. */

#line 648 "AB09JV.f"
	if (*nv > 0) {
#line 649 "AB09JV.f"
	    tolinf *= dlange_("1", nv, nv, &ev[ev_offset], ldev, &dwork[1], (
		    ftnlen)1);

/*           Reduce (AV,EV), or (AV',EV') or (EV',AV') to a generalized */
/*           real Schur form using an orthogonal equivalence */
/*           transformation and apply the orthogonal transformation */
/*           appropriately to BV and CV, or CV' and BV'. */

/*           Workspace needed:  2*NV*NV + MAX( 11*NV+16, NV*P, NV*PV ); */
/*                              prefer larger. */

#line 659 "AB09JV.f"
	    kq = 1;
#line 660 "AB09JV.f"
	    kz = kq + *nv * *nv;
#line 661 "AB09JV.f"
	    kar = kz + *nv * *nv;
#line 662 "AB09JV.f"
	    kai = kar + *nv;
#line 663 "AB09JV.f"
	    kb = kai + *nv;
#line 664 "AB09JV.f"
	    kw = kb + *nv;

#line 666 "AB09JV.f"
	    if (conjs) {
#line 667 "AB09JV.f"
		*(unsigned char *)stdom = 'S';
#line 668 "AB09JV.f"
		alpha += sqrt(tolinf);

/*              Transpose AV and EV, if non-scalar. */

#line 672 "AB09JV.f"
		i__1 = *nv - 1;
#line 672 "AB09JV.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 673 "AB09JV.f"
		    i__2 = *nv - i__;
#line 673 "AB09JV.f"
		    dswap_(&i__2, &av[i__ + 1 + i__ * av_dim1], &c__1, &av[
			    i__ + (i__ + 1) * av_dim1], ldav);
#line 674 "AB09JV.f"
		    i__2 = *nv - i__;
#line 674 "AB09JV.f"
		    dswap_(&i__2, &ev[i__ + 1 + i__ * ev_dim1], &c__1, &ev[
			    i__ + (i__ + 1) * ev_dim1], ldev);
#line 675 "AB09JV.f"
/* L10: */
#line 675 "AB09JV.f"
		}

#line 677 "AB09JV.f"
		if (discr) {

/*                 Reduce (EV',AV') to a generalized real Schur form */
/*                 using orthogonal transformation matrices Q and Z */
/*                 such that Q'*EV'*Z results in a quasi-triangular form */
/*                 and Q'*AV'*Z results upper triangular. */
/*                 Total workspace needed: 2*NV*NV + 11*NV + 16. */

#line 685 "AB09JV.f"
		    *(unsigned char *)evtype = 'R';
#line 686 "AB09JV.f"
		    i__1 = *ldwork - kw + 1;
#line 686 "AB09JV.f"
		    dgges_("Vectors", "Vectors", "Not ordered", (L_fp)delctg_,
			     nv, &ev[ev_offset], ldev, &av[av_offset], ldav, &
			    sdim, &dwork[kar], &dwork[kai], &dwork[kb], &
			    dwork[kq], &ldw, &dwork[kz], &ldw, &dwork[kw], &
			    i__1, bwork, &ierr, (ftnlen)7, (ftnlen)7, (ftnlen)
			    11);
#line 691 "AB09JV.f"
		} else {

/*                 Reduce (AV',EV') to a generalized real Schur form */
/*                 using orthogonal transformation matrices Q and Z */
/*                 such that Q'*AV'*Z results in a quasi-triangular form */
/*                 and Q'*EV'*Z results upper triangular. */
/*                 Total workspace needed: 2*NV*NV + 11*NV + 16. */

#line 699 "AB09JV.f"
		    *(unsigned char *)evtype = 'G';
#line 700 "AB09JV.f"
		    i__1 = *ldwork - kw + 1;
#line 700 "AB09JV.f"
		    dgges_("Vectors", "Vectors", "Not ordered", (L_fp)delctg_,
			     nv, &av[av_offset], ldav, &ev[ev_offset], ldev, &
			    sdim, &dwork[kar], &dwork[kai], &dwork[kb], &
			    dwork[kq], &ldw, &dwork[kz], &ldw, &dwork[kw], &
			    i__1, bwork, &ierr, (ftnlen)7, (ftnlen)7, (ftnlen)
			    11);
#line 705 "AB09JV.f"
		}
#line 706 "AB09JV.f"
		if (ierr != 0) {
#line 707 "AB09JV.f"
		    *info = 1;
#line 708 "AB09JV.f"
		    return 0;
#line 709 "AB09JV.f"
		}
#line 710 "AB09JV.f"
		if (stabck) {

/*                 Check stability/antistability of generalized */
/*                 eigenvalues of the pair (AV,EV). */

#line 715 "AB09JV.f"
		    ab09jx_(dico, stdom, evtype, nv, &alpha, &dwork[kar], &
			    dwork[kai], &dwork[kb], &tolinf, &ierr, (ftnlen)1,
			     (ftnlen)1, (ftnlen)1);
#line 718 "AB09JV.f"
		    if (ierr != 0) {
#line 719 "AB09JV.f"
			*info = 4;
#line 720 "AB09JV.f"
			return 0;
#line 721 "AB09JV.f"
		    }
#line 722 "AB09JV.f"
		}
/* Computing MAX */
#line 723 "AB09JV.f"
		d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 723 "AB09JV.f"
		work = max(d__1,d__2);

/*              Compute Z'*BV and CV*Q. */
/*              Total workspace needed: 2*NV*NV + NV*MAX(P,PV). */

#line 728 "AB09JV.f"
		kw = kar;
#line 729 "AB09JV.f"
		dlacpy_("Full", nv, pv, &bv[bv_offset], ldbv, &dwork[kw], &
			ldw, (ftnlen)4);
#line 730 "AB09JV.f"
		dgemm_("T", "N", nv, pv, nv, &c_b24, &dwork[kz], &ldw, &dwork[
			kw], &ldw, &c_b18, &bv[bv_offset], ldbv, (ftnlen)1, (
			ftnlen)1);
#line 732 "AB09JV.f"
		dlacpy_("Full", p, nv, &cv[cv_offset], ldcv, &dwork[kw], p, (
			ftnlen)4);
#line 733 "AB09JV.f"
		dgemm_("N", "N", p, nv, nv, &c_b24, &dwork[kw], p, &dwork[kq],
			 &ldw, &c_b18, &cv[cv_offset], ldcv, (ftnlen)1, (
			ftnlen)1);
#line 735 "AB09JV.f"
	    } else {

/*              Reduce (AV,EV) to a generalized real Schur form */
/*              using orthogonal transformation matrices Q and Z */
/*              such that Q'*AV*Z results in a quasi-triangular form */
/*              and Q'*EV*Z results upper triangular. */
/*              Total workspace needed: 2*NV*NV + 11*NV + 16. */

#line 743 "AB09JV.f"
		*(unsigned char *)stdom = 'U';
#line 744 "AB09JV.f"
		*(unsigned char *)evtype = 'G';
#line 745 "AB09JV.f"
		alpha -= sqrt(tolinf);
#line 746 "AB09JV.f"
		i__1 = *ldwork - kw + 1;
#line 746 "AB09JV.f"
		dgges_("Vectors", "Vectors", "Not ordered", (L_fp)delctg_, nv,
			 &av[av_offset], ldav, &ev[ev_offset], ldev, &sdim, &
			dwork[kar], &dwork[kai], &dwork[kb], &dwork[kq], &ldw,
			 &dwork[kz], &ldw, &dwork[kw], &i__1, bwork, &ierr, (
			ftnlen)7, (ftnlen)7, (ftnlen)11);
#line 751 "AB09JV.f"
		if (ierr != 0) {
#line 752 "AB09JV.f"
		    *info = 1;
#line 753 "AB09JV.f"
		    return 0;
#line 754 "AB09JV.f"
		}
#line 755 "AB09JV.f"
		if (stabck) {

/*                 Check stability/antistability of generalized */
/*                 eigenvalues of the pair (AV,EV). */

#line 760 "AB09JV.f"
		    ab09jx_(dico, stdom, evtype, nv, &alpha, &dwork[kar], &
			    dwork[kai], &dwork[kb], &tolinf, &ierr, (ftnlen)1,
			     (ftnlen)1, (ftnlen)1);
#line 763 "AB09JV.f"
		    if (ierr != 0) {
#line 764 "AB09JV.f"
			*info = 4;
#line 765 "AB09JV.f"
			return 0;
#line 766 "AB09JV.f"
		    }
#line 767 "AB09JV.f"
		}
/* Computing MAX */
#line 768 "AB09JV.f"
		d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 768 "AB09JV.f"
		work = max(d__1,d__2);

/*              Compute Q'*BV and CV*Z. */
/*              Total workspace needed: 2*NV*NV + NV*MAX(P,PV). */

#line 773 "AB09JV.f"
		kw = kar;
#line 774 "AB09JV.f"
		dlacpy_("Full", nv, p, &bv[bv_offset], ldbv, &dwork[kw], &ldw,
			 (ftnlen)4);
#line 775 "AB09JV.f"
		dgemm_("T", "N", nv, p, nv, &c_b24, &dwork[kq], &ldw, &dwork[
			kw], &ldw, &c_b18, &bv[bv_offset], ldbv, (ftnlen)1, (
			ftnlen)1);
#line 777 "AB09JV.f"
		dlacpy_("Full", pv, nv, &cv[cv_offset], ldcv, &dwork[kw], pv, 
			(ftnlen)4);
#line 778 "AB09JV.f"
		dgemm_("N", "N", pv, nv, nv, &c_b24, &dwork[kw], pv, &dwork[
			kz], &ldw, &c_b18, &cv[cv_offset], ldcv, (ftnlen)1, (
			ftnlen)1);
#line 780 "AB09JV.f"
	    }
/* Computing MAX */
#line 781 "AB09JV.f"
	    d__1 = work, d__2 = (doublereal) ((*nv << 1) * *nv + *nv * max(*p,
		    *pv));
#line 781 "AB09JV.f"
	    work = max(d__1,d__2);

#line 783 "AB09JV.f"
	}

#line 785 "AB09JV.f"
	kc = 1;
#line 786 "AB09JV.f"
	kf = kc + *nv * *n;
#line 787 "AB09JV.f"
	ke = kf + *nv * *n;
#line 788 "AB09JV.f"
	kw = ke + *n * *n;
#line 789 "AB09JV.f"
	dlaset_("Full", nv, n, &c_b18, &c_b18, &dwork[kf], &ldw, (ftnlen)4);

#line 791 "AB09JV.f"
	if (conjs) {

/*           Compute the projection of conj(V)*G. */

/*           Total workspace needed: NV*N + MAX( NV*N+N*N, PV*N, PV*M ) */

/*           Compute CV'*C. */
/*           Workspace needed: NV*N. */

#line 800 "AB09JV.f"
	    dgemm_("T", "N", nv, n, p, &c_b24, &cv[cv_offset], ldcv, &c__[
		    c_offset], ldc, &c_b18, &dwork[kc], &ldw, (ftnlen)1, (
		    ftnlen)1);

#line 803 "AB09JV.f"
	    if (discr) {

/*              Compute X and SCALE satisfying */

/*              EV'*X - AV'*X*A = SCALE*CV'*C by solving equivalently */

/*              EV'*X - Y*A = SCALE*CV'*C, */
/*              AV'*X - Y   = 0. */

/*              Additional workspace needed: */
/*              real    NV*N + N*N; */
/*              integer NV+N+6. */

#line 816 "AB09JV.f"
		if (*n > 0) {
#line 817 "AB09JV.f"
		    dlaset_("Full", n, n, &c_b18, &c_b24, &dwork[ke], &ldwn, (
			    ftnlen)4);
#line 819 "AB09JV.f"
		    i__1 = *ldwork - kw + 1;
#line 819 "AB09JV.f"
		    dtgsyl_("N", &c__0, nv, n, &ev[ev_offset], ldev, &a[
			    a_offset], lda, &dwork[kc], &ldw, &av[av_offset], 
			    ldav, &dwork[ke], &ldwn, &dwork[kf], &ldw, &scale,
			     &dif, &dwork[kw], &i__1, &iwork[1], &ierr, (
			    ftnlen)1);
#line 823 "AB09JV.f"
		    if (ierr != 0) {
#line 824 "AB09JV.f"
			*info = 2;
#line 825 "AB09JV.f"
			return 0;
#line 826 "AB09JV.f"
		    }
#line 827 "AB09JV.f"
		}

/*              Construct C <- DV'*C + BV'*X*A/SCALE, */
/*                        D <- DV'*D + BV'*X*B/SCALE. */

/*              Additional workspace needed: MAX( PV*N, PV*M ). */

/*              C <- DV'*C. */

#line 836 "AB09JV.f"
		kw = kf;
#line 837 "AB09JV.f"
		dgemm_("T", "N", pv, n, p, &c_b24, &dv[dv_offset], lddv, &c__[
			c_offset], ldc, &c_b18, &dwork[kw], pv, (ftnlen)1, (
			ftnlen)1);
#line 839 "AB09JV.f"
		dlacpy_("Full", pv, n, &dwork[kw], pv, &c__[c_offset], ldc, (
			ftnlen)4);

/*              D <- DV'*D. */

#line 843 "AB09JV.f"
		dgemm_("T", "N", pv, m, p, &c_b24, &dv[dv_offset], lddv, &d__[
			d_offset], ldd, &c_b18, &dwork[kw], pv, (ftnlen)1, (
			ftnlen)1);
#line 845 "AB09JV.f"
		dlacpy_("Full", pv, m, &dwork[kw], pv, &d__[d_offset], ldd, (
			ftnlen)4);

/*              C <- C + BV'*X*A/SCALE. */

#line 849 "AB09JV.f"
		d__1 = 1. / scale;
#line 849 "AB09JV.f"
		dgemm_("T", "N", pv, n, nv, &d__1, &bv[bv_offset], ldbv, &
			dwork[kc], &ldw, &c_b18, &dwork[kw], pv, (ftnlen)1, (
			ftnlen)1);
#line 851 "AB09JV.f"
		dgemm_("N", "N", pv, n, n, &c_b24, &dwork[kw], pv, &a[
			a_offset], lda, &c_b24, &c__[c_offset], ldc, (ftnlen)
			1, (ftnlen)1);

/*              D <- D + BV'*X*B/SCALE. */

#line 856 "AB09JV.f"
		dgemm_("N", "N", pv, m, n, &c_b24, &dwork[kw], pv, &b[
			b_offset], ldb, &c_b24, &d__[d_offset], ldd, (ftnlen)
			1, (ftnlen)1);
#line 858 "AB09JV.f"
	    } else {

/*              Compute X and SCALE satisfying */

/*              AV'*X + EV'*X*A + SCALE*CV'*C = 0 by solving equivalently */

/*              AV'*X - Y*A    = -SCALE*CV'*C, */
/*              EV'*X - Y*(-I) = 0. */

/*              Additional workspace needed: */
/*              real    NV*N+N*N; */
/*              integer NV+N+6. */

#line 871 "AB09JV.f"
		if (*n > 0) {
#line 872 "AB09JV.f"
		    dlaset_("Full", n, n, &c_b18, &c_b17, &dwork[ke], &ldwn, (
			    ftnlen)4);
#line 874 "AB09JV.f"
		    i__1 = *ldwork - kw + 1;
#line 874 "AB09JV.f"
		    dtgsyl_("N", &c__0, nv, n, &av[av_offset], ldav, &a[
			    a_offset], lda, &dwork[kc], &ldw, &ev[ev_offset], 
			    ldev, &dwork[ke], &ldwn, &dwork[kf], &ldw, &scale,
			     &dif, &dwork[kw], &i__1, &iwork[1], &ierr, (
			    ftnlen)1);

/*                 Note that the computed solution in DWORK(KC) is -X. */

#line 881 "AB09JV.f"
		    if (ierr != 0) {
#line 882 "AB09JV.f"
			*info = 2;
#line 883 "AB09JV.f"
			return 0;
#line 884 "AB09JV.f"
		    }
#line 885 "AB09JV.f"
		}

/*              Construct C <- DV'*C + BV'*X/SCALE. */

#line 889 "AB09JV.f"
		kw = kf;
#line 890 "AB09JV.f"
		dgemm_("T", "N", pv, n, p, &c_b24, &dv[dv_offset], lddv, &c__[
			c_offset], ldc, &c_b18, &dwork[kw], pv, (ftnlen)1, (
			ftnlen)1);
#line 892 "AB09JV.f"
		dlacpy_("Full", pv, n, &dwork[kw], pv, &c__[c_offset], ldc, (
			ftnlen)4);
#line 893 "AB09JV.f"
		d__1 = -1. / scale;
#line 893 "AB09JV.f"
		dgemm_("T", "N", pv, n, nv, &d__1, &bv[bv_offset], ldbv, &
			dwork[kc], &ldw, &c_b24, &c__[c_offset], ldc, (ftnlen)
			1, (ftnlen)1);

/*              Construct D <- DV'*D. */

#line 898 "AB09JV.f"
		dgemm_("T", "N", pv, m, p, &c_b24, &dv[dv_offset], lddv, &d__[
			d_offset], ldd, &c_b18, &dwork[kw], pv, (ftnlen)1, (
			ftnlen)1);
#line 900 "AB09JV.f"
		dlacpy_("Full", pv, m, &dwork[kw], pv, &d__[d_offset], ldd, (
			ftnlen)4);
#line 901 "AB09JV.f"
	    }
#line 902 "AB09JV.f"
	} else {

/*           Compute the projection of V*G. */

/*           Total workspace needed: NV*N + MAX( NV*N+N*N, PV*N, PV*M ) */

/*           Compute -BV*C. */
/*           Workspace needed: NV*N. */

#line 911 "AB09JV.f"
	    dgemm_("N", "N", nv, n, p, &c_b17, &bv[bv_offset], ldbv, &c__[
		    c_offset], ldc, &c_b18, &dwork[1], &ldw, (ftnlen)1, (
		    ftnlen)1);

/*           Compute X and SCALE satisfying */

/*           AV*X - EV*X*A + SCALE*BV*C = 0 by solving equivalently */

/*           AV*X - Y*A = -SCALE*BV*C, */
/*           EV*X - Y   = 0. */

/*           Additional workspace needed: */
/*           real    NV*N + N*N; */
/*           integer NV+N+6. */

#line 925 "AB09JV.f"
	    if (*n > 0) {
#line 926 "AB09JV.f"
		dlaset_("Full", n, n, &c_b18, &c_b24, &dwork[ke], &ldwn, (
			ftnlen)4);
#line 927 "AB09JV.f"
		i__1 = *ldwork - kw + 1;
#line 927 "AB09JV.f"
		dtgsyl_("N", &c__0, nv, n, &av[av_offset], ldav, &a[a_offset],
			 lda, &dwork[kc], &ldw, &ev[ev_offset], ldev, &dwork[
			ke], &ldwn, &dwork[kf], &ldw, &scale, &dif, &dwork[kw]
			, &i__1, &iwork[1], &ierr, (ftnlen)1);
#line 931 "AB09JV.f"
		if (ierr != 0) {
#line 932 "AB09JV.f"
		    *info = 2;
#line 933 "AB09JV.f"
		    return 0;
#line 934 "AB09JV.f"
		}
#line 935 "AB09JV.f"
	    }

/*           Construct C <- DV*C + CV*X/SCALE. */

#line 939 "AB09JV.f"
	    kw = kf;
#line 940 "AB09JV.f"
	    dgemm_("N", "N", pv, n, p, &c_b24, &dv[dv_offset], lddv, &c__[
		    c_offset], ldc, &c_b18, &dwork[kw], pv, (ftnlen)1, (
		    ftnlen)1);
#line 942 "AB09JV.f"
	    dlacpy_("Full", pv, n, &dwork[kw], pv, &c__[c_offset], ldc, (
		    ftnlen)4);
#line 943 "AB09JV.f"
	    d__1 = 1. / scale;
#line 943 "AB09JV.f"
	    dgemm_("N", "N", pv, n, nv, &d__1, &cv[cv_offset], ldcv, &dwork[1]
		    , &ldw, &c_b24, &c__[c_offset], ldc, (ftnlen)1, (ftnlen)1)
		    ;

/*           Construct D <- DV*D. */

#line 948 "AB09JV.f"
	    dgemm_("N", "N", pv, m, p, &c_b24, &dv[dv_offset], lddv, &d__[
		    d_offset], ldd, &c_b18, &dwork[kw], pv, (ftnlen)1, (
		    ftnlen)1);
#line 950 "AB09JV.f"
	    dlacpy_("Full", pv, m, &dwork[kw], pv, &d__[d_offset], ldd, (
		    ftnlen)4);
#line 951 "AB09JV.f"
	}
#line 952 "AB09JV.f"
    }

/* Computing MAX */
#line 954 "AB09JV.f"
    d__1 = work, d__2 = (doublereal) lw;
#line 954 "AB09JV.f"
    dwork[1] = max(d__1,d__2);

#line 956 "AB09JV.f"
    return 0;
/* *** Last line of AB09JV *** */
} /* ab09jv_ */

