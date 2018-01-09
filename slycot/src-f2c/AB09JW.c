#line 1 "AB09JW.f"
/* AB09JW.f -- translated by f2c (version 20100827).
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

#line 1 "AB09JW.f"
/* Table of constant values */

static doublereal c_b12 = 0.;
static doublereal c_b23 = -1.;
static integer c_n1 = -1;
static doublereal c_b30 = 1.;
static integer c__1 = 1;
static integer c__0 = 0;

/* Subroutine */ int ab09jw_(char *job, char *dico, char *jobew, char *stbchk,
	 integer *n, integer *m, integer *p, integer *nw, integer *mw, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *d__, integer *ldd, doublereal *aw, 
	integer *ldaw, doublereal *ew, integer *ldew, doublereal *bw, integer 
	*ldbw, doublereal *cw, integer *ldcw, doublereal *dw, integer *lddw, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen job_len, ftnlen dico_len, ftnlen jobew_len, ftnlen stbchk_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, aw_dim1, aw_offset, b_dim1, b_offset, bw_dim1, 
	    bw_offset, c_dim1, c_offset, cw_dim1, cw_offset, d_dim1, d_offset,
	     dw_dim1, dw_offset, ew_dim1, ew_offset, i__1, i__2, i__3, i__4, 
	    i__5, i__6;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, ia, kb, kc, ke, kf, kq, kw, lw, kz;
    static doublereal dif;
    static integer kai, kar, ldw, sdim, ierr, ldwm, ldwn, ldwp;
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
    static logical unitew;
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
/*     projection of G*W or G*conj(W) containing the poles of G, from the */
/*     state-space representations (A,B,C,D) and (AW-lambda*EW,BW,CW,DW), */
/*     of the transfer-function matrices G and W, respectively. */
/*     G is assumed to be a stable transfer-function matrix and */
/*     the state matrix A must be in a real Schur form. */
/*     When computing the stable projection of G*W, it is assumed */
/*     that G and W have completely distinct poles. */
/*     When computing the stable projection of G*conj(W), it is assumed */
/*     that G and conj(W) have completely distinct poles. */

/*     Note: For a transfer-function matrix G, conj(G) denotes the */
/*     conjugate of G given by G'(-s) for a continuous-time system or */
/*     G'(1/z) for a discrete-time system. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the projection to be computed as follows: */
/*             = 'W':  compute the projection of G*W containing */
/*                     the poles of G; */
/*             = 'C':  compute the projection of G*conj(W) containing */
/*                     the poles of G. */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the systems as follows: */
/*             = 'C':  G and W are continuous-time systems; */
/*             = 'D':  G and W are discrete-time systems. */

/*     JOBEW   CHARACTER*1 */
/*             Specifies whether EW is a general square or an identity */
/*             matrix as follows: */
/*             = 'G':  EW is a general square matrix; */
/*             = 'I':  EW is the identity matrix. */

/*     STBCHK  CHARACTER*1 */
/*             Specifies whether stability/antistability of W is to be */
/*             checked as follows: */
/*             = 'C':  check stability if JOB = 'C' or antistability if */
/*                     JOB = 'W'; */
/*             = 'N':  do not check stability or antistability. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The dimension of the state vector of the system with */
/*             the transfer-function matrix G.  N >= 0. */

/*     M       (input) INTEGER */
/*             The dimension of the input vector of the system with */
/*             the transfer-function matrix G, and also the dimension */
/*             of the output vector if JOB = 'W', or of the input vector */
/*             if JOB = 'C', of the system with the transfer-function */
/*             matrix W.  M >= 0. */

/*     P       (input) INTEGER */
/*             The dimension of the output vector of the system with the */
/*             transfer-function matrix G.  P >= 0. */

/*     NW      (input) INTEGER */
/*             The dimension of the state vector of the system with the */
/*             transfer-function matrix W.  NW >= 0. */

/*     MW      (input) INTEGER */
/*             The dimension of the input vector, if JOB = 'W', or of */
/*             the output vector, if JOB = 'C', of the system with the */
/*             transfer-function matrix W.  MW >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             state matrix A of the system with the transfer-function */
/*             matrix G in a real Schur form. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, */
/*             dimension (LDB,MAX(M,MW)) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input matrix B of the system with the */
/*             transfer-function matrix G. */
/*             On exit, if INFO = 0, the leading N-by-MW part of this */
/*             array contains the input matrix BS of the projection of */
/*             G*W, if JOB = 'W', or of G*conj(W), if JOB = 'C'. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array must contain */
/*             the output/state matrix C of the system with the */
/*             transfer-function matrix G. The matrix CS is equal to C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= MAX(1,P). */

/*     D       (input/output) DOUBLE PRECISION array, */
/*             dimension (LDB,MAX(M,MW)) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the feedthrough matrix D of the system with */
/*             the transfer-function matrix G. */
/*             On exit, if INFO = 0, the leading P-by-MW part of */
/*             this array contains the feedthrough matrix DS of the */
/*             projection of G*W, if JOB = 'W', or of G*conj(W), */
/*             if JOB = 'C'. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= MAX(1,P). */

/*     AW      (input/output) DOUBLE PRECISION array, dimension (LDAW,NW) */
/*             On entry, the leading NW-by-NW part of this array must */
/*             contain the state matrix AW of the system with the */
/*             transfer-function matrix W. */
/*             On exit, if INFO = 0, the leading NW-by-NW part of this */
/*             array contains a condensed matrix as follows: */
/*             if JOBEW = 'I', it contains the real Schur form of AW; */
/*             if JOBEW = 'G' and JOB = 'W', it contains a quasi-upper */
/*             triangular matrix representing the real Schur matrix */
/*             in the real generalized Schur form of the pair (AW,EW); */
/*             if JOBEW = 'G', JOB = 'C' and DICO = 'C', it contains a */
/*             quasi-upper triangular matrix corresponding to the */
/*             generalized real Schur form of the pair (AW',EW'); */
/*             if JOBEW = 'G', JOB = 'C' and DICO = 'D', it contains an */
/*             upper triangular matrix corresponding to the generalized */
/*             real Schur form of the pair (EW',AW'). */

/*     LDAW    INTEGER */
/*             The leading dimension of the array AW.  LDAW >= MAX(1,NW). */

/*     EW      (input/output) DOUBLE PRECISION array, dimension (LDEW,NW) */
/*             On entry, if JOBEW = 'G', the leading NW-by-NW part of */
/*             this array must contain the descriptor matrix EW of the */
/*             system with the transfer-function matrix W. */
/*             If JOBEW = 'I', EW is assumed to be an identity matrix */
/*             and is not referenced. */
/*             On exit, if INFO = 0 and JOBEW = 'G', the leading NW-by-NW */
/*             part of this array contains a condensed matrix as follows: */
/*             if JOB = 'W', it contains an upper triangular matrix */
/*             corresponding to the real generalized Schur form of the */
/*             pair (AW,EW); */
/*             if JOB = 'C' and DICO = 'C', it contains an upper */
/*             triangular matrix corresponding to the generalized real */
/*             Schur form of the pair (AW',EW'); */
/*             if JOB = 'C' and DICO = 'D', it contains a quasi-upper */
/*             triangular matrix corresponding to the generalized */
/*             real Schur form of the pair (EW',AW'). */

/*     LDEW    INTEGER */
/*             The leading dimension of the array EW. */
/*             LDEW >= MAX(1,NW), if JOBEW = 'G'; */
/*             LDEW >= 1,         if JOBEW = 'I'. */

/*     BW      (input/output) DOUBLE PRECISION array, */
/*             dimension (LDBW,MBW), where MBW = MW, if JOB = 'W', and */
/*             MBW = M, if JOB = 'C'. */
/*             On entry, the leading NW-by-MBW part of this array must */
/*             contain the input matrix BW of the system with the */
/*             transfer-function matrix W. */
/*             On exit, if INFO = 0, the leading NW-by-MBW part of this */
/*             array contains Q'*BW, where Q is the orthogonal matrix */
/*             that reduces AW to the real Schur form or the left */
/*             orthogonal matrix used to reduce the pair (AW,EW), */
/*             (AW',EW') or (EW',AW') to the generalized real Schur form. */

/*     LDBW    INTEGER */
/*             The leading dimension of the array BW.  LDBW >= MAX(1,NW). */

/*     CW      (input/output) DOUBLE PRECISION array, dimension (LDCW,NW) */
/*             On entry, the leading PCW-by-NW part of this array must */
/*             contain the output matrix CW of the system with the */
/*             transfer-function matrix W, where PCW = M if JOB = 'W' or */
/*             PCW = MW if JOB = 'C'. */
/*             On exit, if INFO = 0, the leading PCW-by-NW part of this */
/*             array contains CW*Q, where Q is the orthogonal matrix that */
/*             reduces AW to the real Schur form, or CW*Z, where Z is the */
/*             right orthogonal matrix used to reduce the pair (AW,EW), */
/*             (AW',EW') or (EW',AW') to the generalized real Schur form. */

/*     LDCW    INTEGER */
/*             The leading dimension of the array CW. */
/*             LDCW >= MAX(1,PCW), where PCW = M if JOB = 'W', or */
/*             PCW = MW if JOB = 'C'. */

/*     DW      (input) DOUBLE PRECISION array, */
/*             dimension (LDDW,MBW), where MBW = MW if JOB = 'W', and */
/*             MBW = M if JOB = 'C'. */
/*             The leading PCW-by-MBW part of this array must contain */
/*             the feedthrough matrix DW of the system with the */
/*             transfer-function matrix W, where PCW = M if JOB = 'W', */
/*             or PCW = MW if JOB = 'C'. */

/*     LDDW    INTEGER */
/*             LDDW >= MAX(1,PCW), where PCW = M if JOB = 'W', or */
/*             PCW = MW if JOB = 'C'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK =   0,    if JOBEW = 'I'; */
/*             LIWORK = NW+N+6, if JOBEW = 'G'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= LW1, if JOBEW = 'I', */
/*             LDWORK >= LW2, if JOBEW = 'G', where */
/*               LW1 = MAX( 1, NW*(NW+5), NW*N + MAX( a, N*MW, P*MW ) ) */
/*                     a = 0,    if DICO = 'C' or  JOB = 'W', */
/*                     a = 2*NW, if DICO = 'D' and JOB = 'C'; */
/*               LW2 = MAX( 2*NW*NW + MAX( 11*NW+16, NW*M, MW*NW ), */
/*                          NW*N + MAX( NW*N+N*N, MW*N, P*MW ) ). */
/*             For good performance, LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             =  0:  successful exit; */
/*             <  0:  if INFO = -i, the i-th argument had an illegal */
/*                    value; */
/*             =  1:  the reduction of the pair (AW,EW) to the real */
/*                    generalized Schur form failed (JOBEW = 'G'), */
/*                    or the reduction of the matrix AW to the real */
/*                    Schur form failed (JOBEW = 'I); */
/*             =  2:  the solution of the Sylvester equation failed */
/*                    because the matrix A and the pencil AW-lambda*EW */
/*                    have common eigenvalues (if JOB = 'W'), or the */
/*                    pencil -AW-lambda*EW and A have common eigenvalues */
/*                    (if JOB = 'C' and DICO = 'C'), or the pencil */
/*                    AW-lambda*EW has an eigenvalue which is the */
/*                    reciprocal of one of eigenvalues of A */
/*                    (if JOB = 'C' and DICO = 'D'); */
/*             =  3:  the solution of the Sylvester equation failed */
/*                    because the matrices A and AW have common */
/*                    eigenvalues (if JOB = 'W'), or the matrices A */
/*                    and -AW have common eigenvalues (if JOB = 'C' and */
/*                    DICO = 'C'), or the matrix A has an eigenvalue */
/*                    which is the reciprocal of one of eigenvalues of AW */
/*                    (if JOB = 'C' and DICO = 'D'); */
/*             =  4:  JOB = 'W' and the pair (AW,EW) has not completely */
/*                    unstable generalized eigenvalues, or JOB = 'C' and */
/*                    the pair (AW,EW) has not completely stable */
/*                    generalized eigenvalues. */

/*     METHOD */

/*     If JOB = 'W', the matrices of the stable projection of G*W are */
/*     computed as */

/*       BS = B*DW + Y*BW,  CS = C,  DS = D*DW, */

/*     where Y satisfies the generalized Sylvester equation */

/*       -A*Y*EW + Y*AW + B*CW = 0. */

/*     If JOB = 'C', the matrices of the stable projection of G*conj(W) */
/*     are computed using the following formulas: */

/*     - for a continuous-time system, the matrices BS, CS and DS of */
/*       the stable projection are computed as */

/*         BS = B*DW' + Y*CW',  CS = C,  DS = D*DW', */

/*       where Y satisfies the generalized Sylvester equation */

/*         A*Y*EW' + Y*AW' + B*BW' = 0. */

/*     - for a discrete-time system, the matrices BS, CS and DS of */
/*       the stable projection are computed as */

/*         BS = B*DW' + A*Y*CW',  CS = C,  DS = D*DW' + C*Y*CW', */

/*       where Y satisfies the generalized Sylvester equation */

/*         Y*EW' - A*Y*AW' = B*BW'. */

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

#line 374 "AB09JW.f"
    /* Parameter adjustments */
#line 374 "AB09JW.f"
    a_dim1 = *lda;
#line 374 "AB09JW.f"
    a_offset = 1 + a_dim1;
#line 374 "AB09JW.f"
    a -= a_offset;
#line 374 "AB09JW.f"
    b_dim1 = *ldb;
#line 374 "AB09JW.f"
    b_offset = 1 + b_dim1;
#line 374 "AB09JW.f"
    b -= b_offset;
#line 374 "AB09JW.f"
    c_dim1 = *ldc;
#line 374 "AB09JW.f"
    c_offset = 1 + c_dim1;
#line 374 "AB09JW.f"
    c__ -= c_offset;
#line 374 "AB09JW.f"
    d_dim1 = *ldd;
#line 374 "AB09JW.f"
    d_offset = 1 + d_dim1;
#line 374 "AB09JW.f"
    d__ -= d_offset;
#line 374 "AB09JW.f"
    aw_dim1 = *ldaw;
#line 374 "AB09JW.f"
    aw_offset = 1 + aw_dim1;
#line 374 "AB09JW.f"
    aw -= aw_offset;
#line 374 "AB09JW.f"
    ew_dim1 = *ldew;
#line 374 "AB09JW.f"
    ew_offset = 1 + ew_dim1;
#line 374 "AB09JW.f"
    ew -= ew_offset;
#line 374 "AB09JW.f"
    bw_dim1 = *ldbw;
#line 374 "AB09JW.f"
    bw_offset = 1 + bw_dim1;
#line 374 "AB09JW.f"
    bw -= bw_offset;
#line 374 "AB09JW.f"
    cw_dim1 = *ldcw;
#line 374 "AB09JW.f"
    cw_offset = 1 + cw_dim1;
#line 374 "AB09JW.f"
    cw -= cw_offset;
#line 374 "AB09JW.f"
    dw_dim1 = *lddw;
#line 374 "AB09JW.f"
    dw_offset = 1 + dw_dim1;
#line 374 "AB09JW.f"
    dw -= dw_offset;
#line 374 "AB09JW.f"
    --iwork;
#line 374 "AB09JW.f"
    --dwork;
#line 374 "AB09JW.f"

#line 374 "AB09JW.f"
    /* Function Body */
#line 374 "AB09JW.f"
    conjs = lsame_(job, "C", (ftnlen)1, (ftnlen)1);
#line 375 "AB09JW.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 376 "AB09JW.f"
    unitew = lsame_(jobew, "I", (ftnlen)1, (ftnlen)1);
#line 377 "AB09JW.f"
    stabck = lsame_(stbchk, "C", (ftnlen)1, (ftnlen)1);

#line 379 "AB09JW.f"
    *info = 0;
#line 380 "AB09JW.f"
    if (unitew) {
#line 381 "AB09JW.f"
	if (discr && conjs) {
#line 382 "AB09JW.f"
	    ia = *nw << 1;
#line 383 "AB09JW.f"
	} else {
#line 384 "AB09JW.f"
	    ia = 0;
#line 385 "AB09JW.f"
	}
/* Computing MAX */
/* Computing MAX */
#line 386 "AB09JW.f"
	i__3 = ia, i__4 = *n * *mw, i__3 = max(i__3,i__4), i__4 = *p * *mw;
#line 386 "AB09JW.f"
	i__1 = 1, i__2 = *nw * (*nw + 5), i__1 = max(i__1,i__2), i__2 = *nw * 
		*n + max(i__3,i__4);
#line 386 "AB09JW.f"
	lw = max(i__1,i__2);
#line 387 "AB09JW.f"
    } else {
/* Computing MAX */
/* Computing MAX */
#line 388 "AB09JW.f"
	i__3 = *nw * 11 + 16, i__4 = *nw * *m, i__3 = max(i__3,i__4), i__4 = *
		mw * *nw;
/* Computing MAX */
#line 388 "AB09JW.f"
	i__5 = *nw * *n + *n * *n, i__6 = *mw * *n, i__5 = max(i__5,i__6), 
		i__6 = *p * *mw;
#line 388 "AB09JW.f"
	i__1 = (*nw << 1) * *nw + max(i__3,i__4), i__2 = *nw * *n + max(i__5,
		i__6);
#line 388 "AB09JW.f"
	lw = max(i__1,i__2);
#line 390 "AB09JW.f"
    }

/*     Test the input scalar arguments. */

#line 394 "AB09JW.f"
    ldw = max(1,*nw);
#line 395 "AB09JW.f"
    ldwm = max(1,*mw);
#line 396 "AB09JW.f"
    ldwn = max(1,*n);
#line 397 "AB09JW.f"
    ldwp = max(1,*p);
#line 398 "AB09JW.f"
    if (! (lsame_(job, "W", (ftnlen)1, (ftnlen)1) || conjs)) {
#line 399 "AB09JW.f"
	*info = -1;
#line 400 "AB09JW.f"
    } else if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 401 "AB09JW.f"
	*info = -2;
#line 402 "AB09JW.f"
    } else if (! (lsame_(jobew, "G", (ftnlen)1, (ftnlen)1) || unitew)) {
#line 403 "AB09JW.f"
	*info = -3;
#line 404 "AB09JW.f"
    } else if (! (lsame_(stbchk, "N", (ftnlen)1, (ftnlen)1) || stabck)) {
#line 405 "AB09JW.f"
	*info = -4;
#line 406 "AB09JW.f"
    } else if (*n < 0) {
#line 407 "AB09JW.f"
	*info = -5;
#line 408 "AB09JW.f"
    } else if (*m < 0) {
#line 409 "AB09JW.f"
	*info = -6;
#line 410 "AB09JW.f"
    } else if (*p < 0) {
#line 411 "AB09JW.f"
	*info = -7;
#line 412 "AB09JW.f"
    } else if (*nw < 0) {
#line 413 "AB09JW.f"
	*info = -8;
#line 414 "AB09JW.f"
    } else if (*mw < 0) {
#line 415 "AB09JW.f"
	*info = -9;
#line 416 "AB09JW.f"
    } else if (*lda < ldwn) {
#line 417 "AB09JW.f"
	*info = -11;
#line 418 "AB09JW.f"
    } else if (*ldb < ldwn) {
#line 419 "AB09JW.f"
	*info = -13;
#line 420 "AB09JW.f"
    } else if (*ldc < ldwp) {
#line 421 "AB09JW.f"
	*info = -15;
#line 422 "AB09JW.f"
    } else if (*ldd < ldwp) {
#line 423 "AB09JW.f"
	*info = -17;
#line 424 "AB09JW.f"
    } else if (*ldaw < ldw) {
#line 425 "AB09JW.f"
	*info = -19;
#line 426 "AB09JW.f"
    } else if (*ldew < 1 || ! unitew && *ldew < *nw) {
#line 427 "AB09JW.f"
	*info = -21;
#line 428 "AB09JW.f"
    } else if (*ldbw < ldw) {
#line 429 "AB09JW.f"
	*info = -23;
#line 430 "AB09JW.f"
    } else if (! conjs && *ldcw < max(1,*m) || conjs && *ldcw < ldwm) {
#line 432 "AB09JW.f"
	*info = -25;
#line 433 "AB09JW.f"
    } else if (! conjs && *lddw < max(1,*m) || conjs && *lddw < ldwm) {
#line 435 "AB09JW.f"
	*info = -27;
#line 436 "AB09JW.f"
    } else if (*ldwork < lw) {
#line 437 "AB09JW.f"
	*info = -30;
#line 438 "AB09JW.f"
    }

#line 440 "AB09JW.f"
    if (*info != 0) {

/*        Error return. */

#line 444 "AB09JW.f"
	i__1 = -(*info);
#line 444 "AB09JW.f"
	xerbla_("AB09JW", &i__1, (ftnlen)6);
#line 445 "AB09JW.f"
	return 0;
#line 446 "AB09JW.f"
    }

/*     Quick return if possible. */

#line 450 "AB09JW.f"
    if (*m == 0) {
#line 451 "AB09JW.f"
	dlaset_("Full", n, mw, &c_b12, &c_b12, &b[b_offset], ldb, (ftnlen)4);
#line 452 "AB09JW.f"
	dlaset_("Full", p, mw, &c_b12, &c_b12, &d__[d_offset], ldd, (ftnlen)4)
		;
#line 453 "AB09JW.f"
	dwork[1] = 1.;
#line 454 "AB09JW.f"
	return 0;
#line 455 "AB09JW.f"
    }

/*     Set options for stability/antistability checking. */

#line 459 "AB09JW.f"
    if (discr) {
#line 460 "AB09JW.f"
	alpha = 1.;
#line 461 "AB09JW.f"
    } else {
#line 462 "AB09JW.f"
	alpha = 0.;
#line 463 "AB09JW.f"
    }

#line 465 "AB09JW.f"
    work = 1.;
#line 466 "AB09JW.f"
    tolinf = dlamch_("Epsilon", (ftnlen)7);

#line 468 "AB09JW.f"
    if (unitew) {

/*        EW is the identity matrix. */

#line 472 "AB09JW.f"
	if (*nw > 0) {

/*           Reduce AW to the real Schur form using an orthogonal */
/*           similarity transformation AW <- Q'*AW*Q and apply the */
/*           transformation to BW and CW: BW <- Q'*BW and CW <- CW*Q. */

/*           Workspace needed:  NW*(NW+5); */
/*                              prefer larger. */

#line 481 "AB09JW.f"
	    kw = *nw * (*nw + 2) + 1;
#line 482 "AB09JW.f"
	    if (conjs) {
#line 483 "AB09JW.f"
		*(unsigned char *)stdom = 'S';
#line 484 "AB09JW.f"
		alpha += sqrt(tolinf);
#line 485 "AB09JW.f"
		i__1 = *ldwork - kw + 1;
#line 485 "AB09JW.f"
		tb01wd_(nw, m, mw, &aw[aw_offset], ldaw, &bw[bw_offset], ldbw,
			 &cw[cw_offset], ldcw, &dwork[(*nw << 1) + 1], nw, &
			dwork[1], &dwork[*nw + 1], &dwork[kw], &i__1, &ierr);
#line 488 "AB09JW.f"
	    } else {
#line 489 "AB09JW.f"
		*(unsigned char *)stdom = 'U';
#line 490 "AB09JW.f"
		alpha -= sqrt(tolinf);
#line 491 "AB09JW.f"
		i__1 = *ldwork - kw + 1;
#line 491 "AB09JW.f"
		tb01wd_(nw, mw, m, &aw[aw_offset], ldaw, &bw[bw_offset], ldbw,
			 &cw[cw_offset], ldcw, &dwork[(*nw << 1) + 1], nw, &
			dwork[1], &dwork[*nw + 1], &dwork[kw], &i__1, &ierr);
#line 494 "AB09JW.f"
	    }
#line 495 "AB09JW.f"
	    if (ierr != 0) {
#line 496 "AB09JW.f"
		*info = 1;
#line 497 "AB09JW.f"
		return 0;
#line 498 "AB09JW.f"
	    }
#line 499 "AB09JW.f"
	    if (stabck) {

/*              Check stability/antistability of eigenvalues of AV. */

#line 503 "AB09JW.f"
		ab09jx_(dico, stdom, "S", nw, &alpha, &dwork[1], &dwork[*nw + 
			1], &dwork[1], &tolinf, &ierr, (ftnlen)1, (ftnlen)1, (
			ftnlen)1);
#line 505 "AB09JW.f"
		if (ierr != 0) {
#line 506 "AB09JW.f"
		    *info = 4;
#line 507 "AB09JW.f"
		    return 0;
#line 508 "AB09JW.f"
		}
#line 509 "AB09JW.f"
	    }

/* Computing MAX */
#line 511 "AB09JW.f"
	    d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 511 "AB09JW.f"
	    work = max(d__1,d__2);

#line 513 "AB09JW.f"
	}

#line 515 "AB09JW.f"
	kw = *nw * *n + 1;
#line 516 "AB09JW.f"
	if (conjs) {

/*           Compute the projection of G*conj(W). */

/*           Total workspace needed:  NW*N + MAX( a, N*MW, P*MW ), where */
/*                                    a = 0,    if DICO = 'C', */
/*                                    a = 2*NW, if DICO = 'D'. */

/*           Compute -BW*B'. */
/*           Workspace needed: NW*N. */

#line 527 "AB09JW.f"
	    dgemm_("N", "T", nw, n, m, &c_b23, &bw[bw_offset], ldbw, &b[
		    b_offset], ldb, &c_b12, &dwork[1], &ldw, (ftnlen)1, (
		    ftnlen)1);

#line 530 "AB09JW.f"
	    if (discr) {

/*              Compute Y' and SCALE satisfying */

/*              AW*Y'*A' - Y' = -SCALE*BW*B'. */

/*              Additional workspace needed: 2*NW. */

#line 538 "AB09JW.f"
		sb04py_("N", "T", &c_n1, nw, n, &aw[aw_offset], ldaw, &a[
			a_offset], lda, &dwork[1], &ldw, &scale, &dwork[kw], &
			ierr, (ftnlen)1, (ftnlen)1);
#line 540 "AB09JW.f"
		if (ierr != 0) {
#line 541 "AB09JW.f"
		    *info = 3;
#line 542 "AB09JW.f"
		    return 0;
#line 543 "AB09JW.f"
		}

/*              Construct BS = B*DW' + A*Y*CW'/SCALE, */
/*                        DS = D*DW' + C*Y*CW'/SCALE. */

/*              Additional workspace needed: MAX( N*MW, P*MW ). */

/*              B <- B*DW'. */

#line 552 "AB09JW.f"
		dgemm_("N", "T", n, mw, m, &c_b30, &b[b_offset], ldb, &dw[
			dw_offset], lddw, &c_b12, &dwork[kw], &ldwn, (ftnlen)
			1, (ftnlen)1);
#line 554 "AB09JW.f"
		dlacpy_("Full", n, mw, &dwork[kw], &ldwn, &b[b_offset], ldb, (
			ftnlen)4);

/*              D <- D*DW'. */

#line 558 "AB09JW.f"
		dgemm_("N", "T", p, mw, m, &c_b30, &d__[d_offset], ldd, &dw[
			dw_offset], lddw, &c_b12, &dwork[kw], &ldwp, (ftnlen)
			1, (ftnlen)1);
#line 560 "AB09JW.f"
		dlacpy_("Full", p, mw, &dwork[kw], &ldwp, &d__[d_offset], ldd,
			 (ftnlen)4);

/*              B <- B + A*Y*CW'/SCALE. */

#line 564 "AB09JW.f"
		d__1 = 1. / scale;
#line 564 "AB09JW.f"
		dgemm_("T", "T", n, mw, nw, &d__1, &dwork[1], &ldw, &cw[
			cw_offset], ldcw, &c_b12, &dwork[kw], &ldwn, (ftnlen)
			1, (ftnlen)1);
#line 566 "AB09JW.f"
		dgemm_("N", "N", n, mw, n, &c_b30, &a[a_offset], lda, &dwork[
			kw], &ldwn, &c_b30, &b[b_offset], ldb, (ftnlen)1, (
			ftnlen)1);

/*              D <- D + C*Y*CW'/SCALE. */

#line 571 "AB09JW.f"
		dgemm_("N", "N", p, mw, n, &c_b30, &c__[c_offset], ldc, &
			dwork[kw], &ldwn, &c_b30, &d__[d_offset], ldd, (
			ftnlen)1, (ftnlen)1);
#line 573 "AB09JW.f"
	    } else {

/*              Compute Y' and SCALE satisfying */

/*              AW*Y' + Y'*A' + SCALE*BW*B' = 0. */

#line 579 "AB09JW.f"
		if (*n > 0) {
#line 580 "AB09JW.f"
		    dtrsyl_("N", "T", &c__1, nw, n, &aw[aw_offset], ldaw, &a[
			    a_offset], lda, &dwork[1], &ldw, &scale, &ierr, (
			    ftnlen)1, (ftnlen)1);
#line 582 "AB09JW.f"
		    if (ierr != 0) {
#line 583 "AB09JW.f"
			*info = 3;
#line 584 "AB09JW.f"
			return 0;
#line 585 "AB09JW.f"
		    }
#line 586 "AB09JW.f"
		}

/*              Construct BS = B*DW' + Y*CW'/SCALE, */
/*                        DS = D*DW'. */

/*              Additional workspace needed: MAX( N*MW, P*MW ). */

/*              Construct B <- B*DW' + Y*CW'/SCALE. */

#line 595 "AB09JW.f"
		dgemm_("N", "T", n, mw, m, &c_b30, &b[b_offset], ldb, &dw[
			dw_offset], lddw, &c_b12, &dwork[kw], &ldwn, (ftnlen)
			1, (ftnlen)1);
#line 597 "AB09JW.f"
		dlacpy_("Full", n, mw, &dwork[kw], &ldwn, &b[b_offset], ldb, (
			ftnlen)4);
#line 598 "AB09JW.f"
		d__1 = 1. / scale;
#line 598 "AB09JW.f"
		dgemm_("T", "T", n, mw, nw, &d__1, &dwork[1], &ldw, &cw[
			cw_offset], ldcw, &c_b30, &b[b_offset], ldb, (ftnlen)
			1, (ftnlen)1);

/*              D <- D*DW'. */

#line 603 "AB09JW.f"
		dgemm_("N", "T", p, mw, m, &c_b30, &d__[d_offset], ldd, &dw[
			dw_offset], lddw, &c_b12, &dwork[kw], &ldwp, (ftnlen)
			1, (ftnlen)1);
#line 605 "AB09JW.f"
		dlacpy_("Full", p, mw, &dwork[kw], &ldwp, &d__[d_offset], ldd,
			 (ftnlen)4);
#line 606 "AB09JW.f"
	    }
#line 607 "AB09JW.f"
	} else {

/*           Compute the projection of G*W. */

/*           Total workspace needed:  NW*N + MAX( N*MW, P*MW ). */

/*           Compute B*CW. */
/*           Workspace needed: N*NW. */

#line 616 "AB09JW.f"
	    dgemm_("N", "N", n, nw, m, &c_b30, &b[b_offset], ldb, &cw[
		    cw_offset], ldcw, &c_b12, &dwork[1], &ldwn, (ftnlen)1, (
		    ftnlen)1);

/*           Compute Y and SCALE satisfying */

/*           A*Y - Y*AW - SCALE*B*CW = 0. */

#line 623 "AB09JW.f"
	    if (*n > 0) {
#line 624 "AB09JW.f"
		dtrsyl_("N", "N", &c_n1, n, nw, &a[a_offset], lda, &aw[
			aw_offset], ldaw, &dwork[1], &ldwn, &scale, &ierr, (
			ftnlen)1, (ftnlen)1);
#line 626 "AB09JW.f"
		if (ierr != 0) {
#line 627 "AB09JW.f"
		    *info = 3;
#line 628 "AB09JW.f"
		    return 0;
#line 629 "AB09JW.f"
		}
#line 630 "AB09JW.f"
	    }

/*           Construct BS = B*DW + Y*BW/SCALE, */
/*                     DS = D*DW. */

/*           Additional workspace needed: MAX( N*MW, P*MW ). */
/*           Construct B <- B*DW + Y*BW/SCALE. */

#line 638 "AB09JW.f"
	    dgemm_("N", "N", n, mw, m, &c_b30, &b[b_offset], ldb, &dw[
		    dw_offset], lddw, &c_b12, &dwork[kw], &ldwn, (ftnlen)1, (
		    ftnlen)1);
#line 640 "AB09JW.f"
	    dlacpy_("Full", n, mw, &dwork[kw], &ldwn, &b[b_offset], ldb, (
		    ftnlen)4);
#line 641 "AB09JW.f"
	    d__1 = 1. / scale;
#line 641 "AB09JW.f"
	    dgemm_("N", "N", n, mw, nw, &d__1, &dwork[1], &ldwn, &bw[
		    bw_offset], ldbw, &c_b30, &b[b_offset], ldb, (ftnlen)1, (
		    ftnlen)1);

/*           D <- D*DW. */

#line 646 "AB09JW.f"
	    dgemm_("N", "N", p, mw, m, &c_b30, &d__[d_offset], ldd, &dw[
		    dw_offset], lddw, &c_b12, &dwork[kw], &ldwp, (ftnlen)1, (
		    ftnlen)1);
#line 648 "AB09JW.f"
	    dlacpy_("Full", p, mw, &dwork[kw], &ldwp, &d__[d_offset], ldd, (
		    ftnlen)4);
#line 649 "AB09JW.f"
	}
#line 650 "AB09JW.f"
    } else {

/*        EW is a general matrix. */

#line 654 "AB09JW.f"
	if (*nw > 0) {
#line 655 "AB09JW.f"
	    tolinf *= dlange_("1", nw, nw, &ew[ew_offset], ldew, &dwork[1], (
		    ftnlen)1);

/*           Reduce (AW,EW), or (AW',EW') or (EW',AW') to a generalized */
/*           real Schur form using an orthogonal equivalence */
/*           transformation and apply the orthogonal transformation */
/*           appropriately to BW and CW, or CW' and BW'. */

/*           Workspace needed:  2*NW*NW + MAX( 11*NW+16, NW*M, MW*NW ); */
/*                              prefer larger. */

#line 665 "AB09JW.f"
	    kq = 1;
#line 666 "AB09JW.f"
	    kz = kq + *nw * *nw;
#line 667 "AB09JW.f"
	    kar = kz + *nw * *nw;
#line 668 "AB09JW.f"
	    kai = kar + *nw;
#line 669 "AB09JW.f"
	    kb = kai + *nw;
#line 670 "AB09JW.f"
	    kw = kb + *nw;

#line 672 "AB09JW.f"
	    if (conjs) {
#line 673 "AB09JW.f"
		*(unsigned char *)stdom = 'S';
#line 674 "AB09JW.f"
		alpha += sqrt(tolinf);

/*              Transpose AW and EW, if non-scalar. */

#line 678 "AB09JW.f"
		i__1 = *nw - 1;
#line 678 "AB09JW.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 679 "AB09JW.f"
		    i__2 = *nw - i__;
#line 679 "AB09JW.f"
		    dswap_(&i__2, &aw[i__ + 1 + i__ * aw_dim1], &c__1, &aw[
			    i__ + (i__ + 1) * aw_dim1], ldaw);
#line 680 "AB09JW.f"
		    i__2 = *nw - i__;
#line 680 "AB09JW.f"
		    dswap_(&i__2, &ew[i__ + 1 + i__ * ew_dim1], &c__1, &ew[
			    i__ + (i__ + 1) * ew_dim1], ldew);
#line 681 "AB09JW.f"
/* L10: */
#line 681 "AB09JW.f"
		}

#line 683 "AB09JW.f"
		if (discr) {

/*                 Reduce (EW',AW') to a generalized real Schur form */
/*                 using orthogonal transformation matrices Q and Z */
/*                 such that Q'*EW'*Z results in a quasi-triangular form */
/*                 and Q'*AW'*Z results upper triangular. */
/*                 Total workspace needed: 2*NW*NW + 11*NW + 16. */

#line 691 "AB09JW.f"
		    *(unsigned char *)evtype = 'R';
#line 692 "AB09JW.f"
		    i__1 = *ldwork - kw + 1;
#line 692 "AB09JW.f"
		    dgges_("Vectors", "Vectors", "Not ordered", (L_fp)delctg_,
			     nw, &ew[ew_offset], ldew, &aw[aw_offset], ldaw, &
			    sdim, &dwork[kar], &dwork[kai], &dwork[kb], &
			    dwork[kq], &ldw, &dwork[kz], &ldw, &dwork[kw], &
			    i__1, bwork, &ierr, (ftnlen)7, (ftnlen)7, (ftnlen)
			    11);
#line 697 "AB09JW.f"
		} else {

/*                 Reduce (AW',EW') to a generalized real Schur form */
/*                 using orthogonal transformation matrices Q and Z */
/*                 such that Q'*AW'*Z results in a quasi-triangular form */
/*                 and Q'*EW'*Z results upper triangular. */
/*                 Total workspace needed: 2*NW*NW + 11*NW + 16. */

#line 705 "AB09JW.f"
		    *(unsigned char *)evtype = 'G';
#line 706 "AB09JW.f"
		    i__1 = *ldwork - kw + 1;
#line 706 "AB09JW.f"
		    dgges_("Vectors", "Vectors", "Not ordered", (L_fp)delctg_,
			     nw, &aw[aw_offset], ldaw, &ew[ew_offset], ldew, &
			    sdim, &dwork[kar], &dwork[kai], &dwork[kb], &
			    dwork[kq], &ldw, &dwork[kz], &ldw, &dwork[kw], &
			    i__1, bwork, &ierr, (ftnlen)7, (ftnlen)7, (ftnlen)
			    11);
#line 711 "AB09JW.f"
		}
#line 712 "AB09JW.f"
		if (ierr != 0) {
#line 713 "AB09JW.f"
		    *info = 1;
#line 714 "AB09JW.f"
		    return 0;
#line 715 "AB09JW.f"
		}
#line 716 "AB09JW.f"
		if (stabck) {

/*                 Check stability/antistability of generalized */
/*                 eigenvalues of the pair (AV,EV). */

#line 721 "AB09JW.f"
		    ab09jx_(dico, stdom, evtype, nw, &alpha, &dwork[kar], &
			    dwork[kai], &dwork[kb], &tolinf, &ierr, (ftnlen)1,
			     (ftnlen)1, (ftnlen)1);
#line 724 "AB09JW.f"
		    if (ierr != 0) {
#line 725 "AB09JW.f"
			*info = 4;
#line 726 "AB09JW.f"
			return 0;
#line 727 "AB09JW.f"
		    }
#line 728 "AB09JW.f"
		}
/* Computing MAX */
#line 729 "AB09JW.f"
		d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 729 "AB09JW.f"
		work = max(d__1,d__2);

/*              Compute Z'*BW and CW*Q. */
/*              Total workspace needed: 2*NW*NW + NW*MAX(M,MW). */

#line 734 "AB09JW.f"
		kw = kar;
#line 735 "AB09JW.f"
		dlacpy_("Full", nw, m, &bw[bw_offset], ldbw, &dwork[kw], &ldw,
			 (ftnlen)4);
#line 736 "AB09JW.f"
		dgemm_("T", "N", nw, m, nw, &c_b30, &dwork[kz], &ldw, &dwork[
			kw], &ldw, &c_b12, &bw[bw_offset], ldbw, (ftnlen)1, (
			ftnlen)1);
#line 738 "AB09JW.f"
		dlacpy_("Full", mw, nw, &cw[cw_offset], ldcw, &dwork[kw], &
			ldwm, (ftnlen)4);
#line 739 "AB09JW.f"
		dgemm_("N", "N", mw, nw, nw, &c_b30, &dwork[kw], &ldwm, &
			dwork[kq], &ldw, &c_b12, &cw[cw_offset], ldcw, (
			ftnlen)1, (ftnlen)1);
#line 741 "AB09JW.f"
	    } else {

/*              Reduce (AW,EW) to a generalized real Schur form */
/*              using orthogonal transformation matrices Q and Z */
/*              such that Q'*AW*Z results in a quasi-triangular form */
/*              and Q'*EW*Z results upper triangular. */
/*              Total workspace needed: 2*NW*NW + 11*NW + 16. */

#line 749 "AB09JW.f"
		*(unsigned char *)stdom = 'U';
#line 750 "AB09JW.f"
		*(unsigned char *)evtype = 'G';
#line 751 "AB09JW.f"
		alpha -= sqrt(tolinf);
#line 752 "AB09JW.f"
		i__1 = *ldwork - kw + 1;
#line 752 "AB09JW.f"
		dgges_("Vectors", "Vectors", "Not ordered", (L_fp)delctg_, nw,
			 &aw[aw_offset], ldaw, &ew[ew_offset], ldew, &sdim, &
			dwork[kar], &dwork[kai], &dwork[kb], &dwork[kq], &ldw,
			 &dwork[kz], &ldw, &dwork[kw], &i__1, bwork, &ierr, (
			ftnlen)7, (ftnlen)7, (ftnlen)11);
#line 757 "AB09JW.f"
		if (ierr != 0) {
#line 758 "AB09JW.f"
		    *info = 1;
#line 759 "AB09JW.f"
		    return 0;
#line 760 "AB09JW.f"
		}
#line 761 "AB09JW.f"
		if (stabck) {

/*                 Check stability/antistability of generalized */
/*                 eigenvalues of the pair (AV,EV). */

#line 766 "AB09JW.f"
		    ab09jx_(dico, stdom, evtype, nw, &alpha, &dwork[kar], &
			    dwork[kai], &dwork[kb], &tolinf, &ierr, (ftnlen)1,
			     (ftnlen)1, (ftnlen)1);
#line 769 "AB09JW.f"
		    if (ierr != 0) {
#line 770 "AB09JW.f"
			*info = 4;
#line 771 "AB09JW.f"
			return 0;
#line 772 "AB09JW.f"
		    }
#line 773 "AB09JW.f"
		}
/* Computing MAX */
#line 774 "AB09JW.f"
		d__1 = work, d__2 = dwork[kw] + (doublereal) (kw - 1);
#line 774 "AB09JW.f"
		work = max(d__1,d__2);

/*              Compute Q'*BW and CW*Z. */
/*              Total workspace needed: 2*NW*NW + NW*MAX(M,MW). */

#line 779 "AB09JW.f"
		kw = kar;
#line 780 "AB09JW.f"
		dlacpy_("Full", nw, mw, &bw[bw_offset], ldbw, &dwork[kw], &
			ldw, (ftnlen)4);
#line 781 "AB09JW.f"
		dgemm_("T", "N", nw, mw, nw, &c_b30, &dwork[kq], &ldw, &dwork[
			kw], &ldw, &c_b12, &bw[bw_offset], ldbw, (ftnlen)1, (
			ftnlen)1);
#line 783 "AB09JW.f"
		dlacpy_("Full", m, nw, &cw[cw_offset], ldcw, &dwork[kw], m, (
			ftnlen)4);
#line 784 "AB09JW.f"
		dgemm_("N", "N", m, nw, nw, &c_b30, &dwork[kw], m, &dwork[kz],
			 &ldw, &c_b12, &cw[cw_offset], ldcw, (ftnlen)1, (
			ftnlen)1);
#line 786 "AB09JW.f"
	    }
/* Computing MAX */
#line 787 "AB09JW.f"
	    d__1 = work, d__2 = (doublereal) ((*nw << 1) * *nw + *nw * max(*m,
		    *mw));
#line 787 "AB09JW.f"
	    work = max(d__1,d__2);

#line 789 "AB09JW.f"
	}

#line 791 "AB09JW.f"
	kc = 1;
#line 792 "AB09JW.f"
	kf = kc + *nw * *n;
#line 793 "AB09JW.f"
	ke = kf + *nw * *n;
#line 794 "AB09JW.f"
	kw = ke + *n * *n;
#line 795 "AB09JW.f"
	dlaset_("Full", n, nw, &c_b12, &c_b12, &dwork[kf], &ldwn, (ftnlen)4);

#line 797 "AB09JW.f"
	if (conjs) {

/*           Compute the projection of G*conj(W). */

/*           Total workspace needed: NW*N + MAX( NW*N+N*N, MW*N, P*MW ) */

/*           Compute B*BW'. */
/*           Workspace needed: N*NW. */

#line 806 "AB09JW.f"
	    dgemm_("N", "T", n, nw, m, &c_b30, &b[b_offset], ldb, &bw[
		    bw_offset], ldbw, &c_b12, &dwork[kc], &ldwn, (ftnlen)1, (
		    ftnlen)1);

#line 809 "AB09JW.f"
	    if (discr) {

/*              Compute Y and SCALE satisfying */

/*              Y*EW' - A*Y*AW' = SCALE*B*BW' by solving equivalently */

/*              A*X - Y*EW' = -SCALE*B*BW', */
/*              X   - Y*AW' = 0. */

/*              Additional workspace needed: */
/*              real    N*NW + N*N; */
/*              integer NW+N+6. */


#line 823 "AB09JW.f"
		if (*n > 0) {
#line 824 "AB09JW.f"
		    dlaset_("Full", n, n, &c_b12, &c_b30, &dwork[ke], &ldwn, (
			    ftnlen)4);
#line 826 "AB09JW.f"
		    i__1 = *ldwork - kw + 1;
#line 826 "AB09JW.f"
		    dtgsyl_("N", &c__0, n, nw, &a[a_offset], lda, &ew[
			    ew_offset], ldew, &dwork[kc], &ldwn, &dwork[ke], &
			    ldwn, &aw[aw_offset], ldaw, &dwork[kf], &ldwn, &
			    scale, &dif, &dwork[kw], &i__1, &iwork[1], &ierr, 
			    (ftnlen)1);

/*                 Note that the computed solution in DWORK(KC) is -Y. */

#line 833 "AB09JW.f"
		    if (ierr != 0) {
#line 834 "AB09JW.f"
			*info = 2;
#line 835 "AB09JW.f"
			return 0;
#line 836 "AB09JW.f"
		    }
#line 837 "AB09JW.f"
		}

/*              Construct BS = B*DW' + A*Y*CW'/SCALE, */
/*                        DS = D*DW' + C*Y*CW'/SCALE. */

/*              Additional workspace needed: MAX( N*MW, P*MW ). */

/*              B <- B*DW'. */

#line 846 "AB09JW.f"
		dgemm_("N", "T", n, mw, m, &c_b30, &b[b_offset], ldb, &dw[
			dw_offset], lddw, &c_b12, &dwork[kw], &ldwn, (ftnlen)
			1, (ftnlen)1);
#line 848 "AB09JW.f"
		dlacpy_("Full", n, mw, &dwork[kw], &ldwn, &b[b_offset], ldb, (
			ftnlen)4);

/*              D <- D*DW'. */

#line 852 "AB09JW.f"
		dgemm_("N", "T", p, mw, m, &c_b30, &d__[d_offset], ldd, &dw[
			dw_offset], lddw, &c_b12, &dwork[kw], &ldwp, (ftnlen)
			1, (ftnlen)1);
#line 854 "AB09JW.f"
		dlacpy_("Full", p, mw, &dwork[kw], &ldwp, &d__[d_offset], ldd,
			 (ftnlen)4);

/*              B <- B + A*Y*CW'/SCALE. */

#line 858 "AB09JW.f"
		d__1 = -1. / scale;
#line 858 "AB09JW.f"
		dgemm_("N", "T", n, mw, nw, &d__1, &dwork[kf], &ldwn, &cw[
			cw_offset], ldcw, &c_b12, &dwork[kw], &ldwn, (ftnlen)
			1, (ftnlen)1);
#line 861 "AB09JW.f"
		dgemm_("N", "N", n, mw, n, &c_b30, &a[a_offset], lda, &dwork[
			kw], &ldwn, &c_b30, &b[b_offset], ldb, (ftnlen)1, (
			ftnlen)1);

/*              D <- D + C*Y*CW'/SCALE. */

#line 866 "AB09JW.f"
		dgemm_("N", "N", p, mw, n, &c_b30, &c__[c_offset], ldc, &
			dwork[kw], &ldwn, &c_b30, &d__[d_offset], ldd, (
			ftnlen)1, (ftnlen)1);
#line 868 "AB09JW.f"
	    } else {

/*              Compute Y and SCALE satisfying */

/*              A*Y*EW' + Y*AW' + SCALE*B*BW' = 0 by solving equivalently */

/*              A*X    - Y*AW' = SCALE*B*BW', */
/*              (-I)*X - Y*EW' = 0. */

/*              Additional workspace needed: */
/*              real    N*NW+N*N; */
/*              integer NW+N+6. */

#line 881 "AB09JW.f"
		if (*n > 0) {
#line 882 "AB09JW.f"
		    dlaset_("Full", n, n, &c_b12, &c_b23, &dwork[ke], &ldwn, (
			    ftnlen)4);
#line 884 "AB09JW.f"
		    i__1 = *ldwork - kw + 1;
#line 884 "AB09JW.f"
		    dtgsyl_("N", &c__0, n, nw, &a[a_offset], lda, &aw[
			    aw_offset], ldaw, &dwork[kc], &ldwn, &dwork[ke], &
			    ldwn, &ew[ew_offset], ldew, &dwork[kf], &ldwn, &
			    scale, &dif, &dwork[kw], &i__1, &iwork[1], &ierr, 
			    (ftnlen)1);
#line 888 "AB09JW.f"
		    if (ierr != 0) {
#line 889 "AB09JW.f"
			*info = 2;
#line 890 "AB09JW.f"
			return 0;
#line 891 "AB09JW.f"
		    }
#line 892 "AB09JW.f"
		}

/*              Construct BS = B*DW' + Y*CW'/SCALE, */
/*                        DS = D*DW'. */

/*              Additional workspace needed: MAX( N*MW, P*MW ). */

/*              Construct B <- B*DW' + Y*CW'/SCALE. */

#line 901 "AB09JW.f"
		dgemm_("N", "T", n, mw, m, &c_b30, &b[b_offset], ldb, &dw[
			dw_offset], lddw, &c_b12, &dwork[kw], &ldwn, (ftnlen)
			1, (ftnlen)1);
#line 903 "AB09JW.f"
		dlacpy_("Full", n, mw, &dwork[kw], &ldwn, &b[b_offset], ldb, (
			ftnlen)4);
#line 904 "AB09JW.f"
		d__1 = 1. / scale;
#line 904 "AB09JW.f"
		dgemm_("N", "T", n, mw, nw, &d__1, &dwork[kf], &ldwn, &cw[
			cw_offset], ldcw, &c_b30, &b[b_offset], ldb, (ftnlen)
			1, (ftnlen)1);

/*              D <- D*DW'. */

#line 909 "AB09JW.f"
		dgemm_("N", "T", p, mw, m, &c_b30, &d__[d_offset], ldd, &dw[
			dw_offset], lddw, &c_b12, &dwork[kw], &ldwp, (ftnlen)
			1, (ftnlen)1);
#line 911 "AB09JW.f"
		dlacpy_("Full", p, mw, &dwork[kw], &ldwp, &d__[d_offset], ldd,
			 (ftnlen)4);
#line 912 "AB09JW.f"
	    }
#line 913 "AB09JW.f"
	} else {

/*           Compute the projection of G*W. */

/*           Total workspace needed: NW*N + MAX( NW*N+N*N, MW*N, P*MW ) */

/*           Compute B*CW. */
/*           Workspace needed: N*NW. */

#line 922 "AB09JW.f"
	    dgemm_("N", "N", n, nw, m, &c_b30, &b[b_offset], ldb, &cw[
		    cw_offset], ldcw, &c_b12, &dwork[kc], &ldwn, (ftnlen)1, (
		    ftnlen)1);

/*           Compute Y and SCALE satisfying */

/*           -A*Y*EW + Y*AW + B*CW = 0 by solving equivalently */

/*           A*X - Y*AW = SCALE*B*CW, */
/*           X   - Y*EW = 0. */

/*           Additional workspace needed: */
/*           real    N*NW + N*N; */
/*           integer NW+N+6. */

#line 936 "AB09JW.f"
	    if (*n > 0) {
#line 937 "AB09JW.f"
		dlaset_("Full", n, n, &c_b12, &c_b30, &dwork[ke], &ldwn, (
			ftnlen)4);
#line 938 "AB09JW.f"
		i__1 = *ldwork - kw + 1;
#line 938 "AB09JW.f"
		dtgsyl_("N", &c__0, n, nw, &a[a_offset], lda, &aw[aw_offset], 
			ldaw, &dwork[kc], &ldwn, &dwork[ke], &ldwn, &ew[
			ew_offset], ldew, &dwork[kf], &ldwn, &scale, &dif, &
			dwork[kw], &i__1, &iwork[1], &ierr, (ftnlen)1);
#line 942 "AB09JW.f"
		if (ierr != 0) {
#line 943 "AB09JW.f"
		    *info = 2;
#line 944 "AB09JW.f"
		    return 0;
#line 945 "AB09JW.f"
		}
#line 946 "AB09JW.f"
	    }

/*           Construct BS = B*DW + Y*BW/SCALE, */
/*                     DS = D*DW. */

/*           Additional workspace needed: MAX( N*MW, P*MW ). */
/*           Construct B <- B*DW + Y*BW/SCALE. */

#line 954 "AB09JW.f"
	    dgemm_("N", "N", n, mw, m, &c_b30, &b[b_offset], ldb, &dw[
		    dw_offset], lddw, &c_b12, &dwork[kw], &ldwn, (ftnlen)1, (
		    ftnlen)1);
#line 956 "AB09JW.f"
	    dlacpy_("Full", n, mw, &dwork[kw], &ldwn, &b[b_offset], ldb, (
		    ftnlen)4);
#line 957 "AB09JW.f"
	    d__1 = 1. / scale;
#line 957 "AB09JW.f"
	    dgemm_("N", "N", n, mw, nw, &d__1, &dwork[kf], &ldwn, &bw[
		    bw_offset], ldbw, &c_b30, &b[b_offset], ldb, (ftnlen)1, (
		    ftnlen)1);

/*           D <- D*DW. */

#line 962 "AB09JW.f"
	    dgemm_("N", "N", p, mw, m, &c_b30, &d__[d_offset], ldd, &dw[
		    dw_offset], lddw, &c_b12, &dwork[kw], &ldwp, (ftnlen)1, (
		    ftnlen)1);
#line 964 "AB09JW.f"
	    dlacpy_("Full", p, mw, &dwork[kw], &ldwp, &d__[d_offset], ldd, (
		    ftnlen)4);
#line 965 "AB09JW.f"
	}
#line 966 "AB09JW.f"
    }

/* Computing MAX */
#line 968 "AB09JW.f"
    d__1 = work, d__2 = (doublereal) lw;
#line 968 "AB09JW.f"
    dwork[1] = max(d__1,d__2);

#line 970 "AB09JW.f"
    return 0;
/* *** Last line of AB09JW *** */
} /* ab09jw_ */

