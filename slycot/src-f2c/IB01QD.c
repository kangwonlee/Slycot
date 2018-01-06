#line 1 "IB01QD.f"
/* IB01QD.f -- translated by f2c (version 20100827).
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

#line 1 "IB01QD.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b45 = 1.;
static doublereal c_b47 = 0.;
static integer c__2 = 2;
static doublereal c_b113 = .66666666666666663;
static doublereal c_b130 = -1.;

/* Subroutine */ int ib01qd_(char *jobx0, char *job, integer *n, integer *m, 
	integer *l, integer *nsmp, doublereal *a, integer *lda, doublereal *
	c__, integer *ldc, doublereal *u, integer *ldu, doublereal *y, 
	integer *ldy, doublereal *x0, doublereal *b, integer *ldb, doublereal 
	*d__, integer *ldd, doublereal *tol, integer *iwork, doublereal *
	dwork, integer *ldwork, integer *iwarn, integer *info, ftnlen 
	jobx0_len, ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, u_dim1, u_offset, y_dim1, y_offset, i__1, i__2, i__3, 
	    i__4, i__5;

    /* Builtin functions */
    double log(doublereal);
    integer pow_ii(integer *, integer *);
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, k, i2, ia, ic, ie, ig, lm, iq, ln, nm, nn, ix, iy, 
	    n2m, ias, ini, igs, ldr;
    static doublereal dum[1];
    static integer iny, isv, ncp1, ldw2, ldw3, igam, inih, lddw, irem, lnob, 
	    ncol, rank;
    static logical ncyc;
    static integer ierr, inir, inis, itau, irhs, nobs;
    static doublereal toll;
    static integer nrow;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    mb04od_(char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, ftnlen), 
	    mb01td_(integer *, doublereal *, integer *, doublereal *, integer 
	    *, doublereal *, integer *), mb02ud_(char *, char *, char *, char 
	    *, integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal rcond;
    extern /* Subroutine */ int mb04oy_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *);
    static logical withb;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical withd;
    static integer isize;
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer itauu;
    static logical first;
    static integer nsmpl;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer jwork;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static integer iupnt, iypnt;
    static logical power2, withx0;
    extern doublereal dlamch_(char *, ftnlen);
    static integer icycle;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen);
    static integer ncycle;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dgelss_(integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);
    static integer inygam;
    extern /* Subroutine */ int dtrcon_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static doublereal rcondu;
    static integer ixsave, iexpon, ixinit, minsmp;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer minwrk, maxwrk, minwls;


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

/*     To estimate the initial state and the system matrices  B  and  D */
/*     of a linear time-invariant (LTI) discrete-time system, given the */
/*     matrix pair  (A,C)  and the input and output trajectories of the */
/*     system. The model structure is : */

/*           x(k+1) = Ax(k) + Bu(k),   k >= 0, */
/*           y(k)   = Cx(k) + Du(k), */

/*     where  x(k)  is the  n-dimensional state vector (at time k), */
/*            u(k)  is the  m-dimensional input vector, */
/*            y(k)  is the  l-dimensional output vector, */
/*     and  A, B, C, and D  are real matrices of appropriate dimensions. */
/*     Matrix  A  is assumed to be in a real Schur form. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBX0   CHARACTER*1 */
/*             Specifies whether or not the initial state should be */
/*             computed, as follows: */
/*             = 'X':  compute the initial state x(0); */
/*             = 'N':  do not compute the initial state (x(0) is known */
/*                     to be zero). */

/*     JOB     CHARACTER*1 */
/*             Specifies which matrices should be computed, as follows: */
/*             = 'B':  compute the matrix B only (D is known to be zero); */
/*             = 'D':  compute the matrices B and D. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     L       (input) INTEGER */
/*             The number of system outputs.  L > 0. */

/*     NSMP    (input) INTEGER */
/*             The number of rows of matrices  U  and  Y  (number of */
/*             samples,  t). */
/*             NSMP >= N*M + a + e,  where */
/*             a = 0,  if  JOBX0 = 'N'; */
/*             a = N,  if  JOBX0 = 'X'; */
/*             e = 0,  if  JOBX0 = 'X'  and  JOB = 'B'; */
/*             e = 1,  if  JOBX0 = 'N'  and  JOB = 'B'; */
/*             e = M,  if  JOB   = 'D'. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             system state matrix  A  in a real Schur form. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading L-by-N part of this array must contain the */
/*             system output matrix  C  (corresponding to the real Schur */
/*             form of  A). */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= L. */

/*     U       (input/output) DOUBLE PRECISION array, dimension (LDU,M) */
/*             On entry, the leading NSMP-by-M part of this array must */
/*             contain the t-by-m input-data sequence matrix  U, */
/*             U = [u_1 u_2 ... u_m].  Column  j  of  U  contains the */
/*             NSMP  values of the j-th input component for consecutive */
/*             time increments. */
/*             On exit, if  JOB = 'D',  the leading NSMP-by-M part of */
/*             this array contains details of the QR factorization of */
/*             the t-by-m matrix  U, possibly computed sequentially */
/*             (see METHOD). */
/*             If  JOB = 'B',  this array is unchanged on exit. */
/*             If M = 0, this array is not referenced. */

/*     LDU     INTEGER */
/*             The leading dimension of the array U. */
/*             LDU >= MAX(1,NSMP),  if M > 0; */
/*             LDU >= 1,            if M = 0. */

/*     Y       (input) DOUBLE PRECISION array, dimension (LDY,L) */
/*             The leading NSMP-by-L part of this array must contain the */
/*             t-by-l output-data sequence matrix  Y, */
/*             Y = [y_1 y_2 ... y_l].  Column  j  of  Y  contains the */
/*             NSMP  values of the j-th output component for consecutive */
/*             time increments. */

/*     LDY     INTEGER */
/*             The leading dimension of the array Y.  LDY >= MAX(1,NSMP). */

/*     X0      (output) DOUBLE PRECISION array, dimension (N) */
/*             If  JOBX0 = 'X',  the estimated initial state of the */
/*             system,  x(0). */
/*             If  JOBX0 = 'N',  x(0)  is set to zero without any */
/*             calculations. */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             If  N > 0,  M > 0,  and  INFO = 0,  the leading N-by-M */
/*             part of this array contains the system input matrix  B */
/*             in the coordinates corresponding to the real Schur form */
/*             of  A. */
/*             If  N = 0  or  M = 0,  this array is not referenced. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= N,  if  N > 0  and  M > 0; */
/*             LDB >= 1,  if  N = 0  or   M = 0. */

/*     D       (output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             If  M > 0,  JOB = 'D',  and  INFO = 0,  the leading */
/*             L-by-M part of this array contains the system input-output */
/*             matrix  D. */
/*             If  M = 0  or  JOB = 'B',  this array is not referenced. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D. */
/*             LDD >= L,  if  M > 0  and  JOB = 'D'; */
/*             LDD >= 1,  if  M = 0  or   JOB = 'B'. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used for estimating the rank of */
/*             matrices. If the user sets  TOL > 0,  then the given value */
/*             of  TOL  is used as a lower bound for the reciprocal */
/*             condition number;  a matrix whose estimated condition */
/*             number is less than  1/TOL  is considered to be of full */
/*             rank.  If the user sets  TOL <= 0,  then  EPS  is used */
/*             instead, where  EPS  is the relative machine precision */
/*             (see LAPACK Library routine DLAMCH).  TOL <= 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK), where */
/*             LIWORK >= N*M + a,            if  JOB = 'B', */
/*             LIWORK >= max( N*M + a, M ),  if  JOB = 'D', */
/*             with  a = 0,  if  JOBX0 = 'N'; */
/*                   a = N,  if  JOBX0 = 'X'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = 0,  DWORK(1) returns the optimal value */
/*             of LDWORK;  DWORK(2)  contains the reciprocal condition */
/*             number of the triangular factor of the QR factorization of */
/*             the matrix  W2  (see METHOD); if  M > 0  and  JOB = 'D', */
/*             DWORK(3)  contains the reciprocal condition number of the */
/*             triangular factor of the QR factorization of  U. */
/*             On exit, if  INFO = -23,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= max( LDW1, min( LDW2, LDW3 ) ),  where */
/*             LDW1 = 2,          if  M = 0  or   JOB = 'B', */
/*             LDW1 = 3,          if  M > 0  and  JOB = 'D', */
/*             LDWa = t*L*(r + 1) + max( N + max( d, f ), 6*r ), */
/*             LDW2 = LDWa,       if  M = 0  or  JOB = 'B', */
/*             LDW2 = max( LDWa, t*L*(r + 1) + 2*M*M + 6*M ), */
/*                                if  M > 0  and JOB = 'D', */
/*             LDWb = (b + r)*(r + 1) + */
/*                     max( q*(r + 1) + N*N*M + c + max( d, f ), 6*r ), */
/*             LDW3 = LDWb,       if  M = 0  or  JOB = 'B', */
/*             LDW3 = max( LDWb, (b + r)*(r + 1) + 2*M*M + 6*M ), */
/*                                if  M > 0  and JOB = 'D', */
/*                r = N*M + a, */
/*                a = 0,                  if  JOBX0 = 'N', */
/*                a = N,                  if  JOBX0 = 'X'; */
/*                b = 0,                  if  JOB   = 'B', */
/*                b = L*M,                if  JOB   = 'D'; */
/*                c = 0,                  if  JOBX0 = 'N', */
/*                c = L*N,                if  JOBX0 = 'X'; */
/*                d = 0,                  if  JOBX0 = 'N', */
/*                d = 2*N*N + N,          if  JOBX0 = 'X'; */
/*                f = 2*r,                if  JOB   = 'B'   or  M = 0, */
/*                f = M + max( 2*r, M ),  if  JOB   = 'D'  and  M > 0; */
/*                q = b + r*L. */
/*             For good performance,  LDWORK  should be larger. */
/*             If  LDWORK >= LDW2  or */
/*                 LDWORK >= t*L*(r + 1) + (b + r)*(r + 1) + N*N*M + c + */
/*                           max( d, f ), */
/*             then standard QR factorizations of the matrices  U  and/or */
/*             W2  (see METHOD) are used. */
/*             Otherwise, the QR factorizations are computed sequentially */
/*             by performing  NCYCLE  cycles, each cycle (except possibly */
/*             the last one) processing  s < t  samples, where  s  is */
/*             chosen from the equation */
/*               LDWORK = s*L*(r + 1) + (b + r)*(r + 1) + N*N*M + c + */
/*                        max( d, f ). */
/*             (s  is at least  N*M+a+e,  the minimum value of  NSMP.) */
/*             The computational effort may increase and the accuracy may */
/*             decrease with the decrease of  s.  Recommended value is */
/*             LDWORK = LDW2,  assuming a large enough cache size, to */
/*             also accommodate  A,  C,  U,  and  Y. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 4:  the least squares problem to be solved has a */
/*                   rank-deficient coefficient matrix. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 2:  the singular value decomposition (SVD) algorithm did */
/*                   not converge. */

/*     METHOD */

/*     An extension and refinement of the method in [1,2] is used. */
/*     Specifically, denoting */

/*           X = [ vec(D')' vec(B)' x0' ]', */

/*     where  vec(M)  is the vector obtained by stacking the columns of */
/*     the matrix  M,  then  X  is the least squares solution of the */
/*     system  S*X = vec(Y),  with the matrix  S = [ diag(U)  W ], */
/*     defined by */

/*           ( U         |     | ... |     |     | ... |     |         ) */
/*           (   U       |  11 | ... |  n1 |  12 | ... |  nm |         ) */
/*       S = (     :     | y   | ... | y   | y   | ... | y   | P*Gamma ), */
/*           (       :   |     | ... |     |     | ... |     |         ) */
/*           (         U |     | ... |     |     | ... |     |         ) */
/*                                                                     ij */
/*     diag(U)  having  L  block rows and columns.  In this formula,  y */
/*     are the outputs of the system for zero initial state computed */
/*     using the following model, for j = 1:m, and for i = 1:n, */
/*            ij          ij                    ij */
/*           x  (k+1) = Ax  (k) + e_i u_j(k),  x  (0) = 0, */

/*            ij          ij */
/*           y  (k)   = Cx  (k), */

/*     where  e_i  is the i-th n-dimensional unit vector,  Gamma  is */
/*     given by */

/*                (     C     ) */
/*                (    C*A    ) */
/*        Gamma = (   C*A^2   ), */
/*                (     :     ) */
/*                ( C*A^(t-1) ) */

/*     and  P  is a permutation matrix that groups together the rows of */
/*     Gamma  depending on the same row of  C,  namely */
/*     [ c_j;  c_j*A;  c_j*A^2; ...  c_j*A^(t-1) ],  for j = 1:L. */
/*     The first block column,  diag(U),  is not explicitly constructed, */
/*     but its structure is exploited. The last block column is evaluated */
/*     using powers of A with exponents 2^k. No interchanges are applied. */
/*     A special QR decomposition of the matrix  S  is computed. Let */
/*     U = q*[ r' 0 ]'  be the QR decomposition of  U,  if  M > 0,  where */
/*     r  is  M-by-M.   Then,  diag(q')  is applied to  W  and  vec(Y). */
/*     The block-rows of  S  and  vec(Y)  are implicitly permuted so that */
/*     matrix  S  becomes */

/*        ( diag(r)  W1 ) */
/*        (    0     W2 ), */

/*     where  W1  has L*M rows. Then, the QR decomposition of  W2 is */
/*     computed (sequentially, if  M > 0) and used to obtain  B  and  x0. */
/*     The intermediate results and the QR decomposition of  U  are */
/*     needed to find  D.  If a triangular factor is too ill conditioned, */
/*     then singular value decomposition (SVD) is employed. SVD is not */
/*     generally needed if the input sequence is sufficiently */
/*     persistently exciting and  NSMP  is large enough. */
/*     If the matrix  W  cannot be stored in the workspace (i.e., */
/*     LDWORK < LDW2),  the QR decompositions of  W2  and  U  are */
/*     computed sequentially. */

/*     REFERENCES */

/*     [1] Verhaegen M., and Varga, A. */
/*         Some Experience with the MOESP Class of Subspace Model */
/*         Identification Methods in Identifying the BO105 Helicopter. */
/*         Report TR R165-94, DLR Oberpfaffenhofen, 1994. */

/*     [2] Sima, V., and Varga, A. */
/*         RASP-IDENT : Subspace Model Identification Programs. */
/*         Deutsche Forschungsanstalt fur Luft- und Raumfahrt e. V., */
/*         Report TR R888-94, DLR Oberpfaffenhofen, Oct. 1994. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is numerically stable. */

/*     FURTHER COMMENTS */

/*     The algorithm for computing the system matrices  B  and  D  is */
/*     less efficient than the MOESP or N4SID algorithms implemented in */
/*     SLICOT Library routine IB01PD, because a large least squares */
/*     problem has to be solved, but the accuracy is better, as the */
/*     computed matrices  B  and  D  are fitted to the input and output */
/*     trajectories. However, if matrix  A  is unstable, the computed */
/*     matrices  B  and  D  could be inaccurate. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Identification methods; least squares solutions; multivariable */
/*     systems; QR decomposition; singular value decomposition. */

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

/*     Check the input parameters. */

#line 381 "IB01QD.f"
    /* Parameter adjustments */
#line 381 "IB01QD.f"
    a_dim1 = *lda;
#line 381 "IB01QD.f"
    a_offset = 1 + a_dim1;
#line 381 "IB01QD.f"
    a -= a_offset;
#line 381 "IB01QD.f"
    c_dim1 = *ldc;
#line 381 "IB01QD.f"
    c_offset = 1 + c_dim1;
#line 381 "IB01QD.f"
    c__ -= c_offset;
#line 381 "IB01QD.f"
    u_dim1 = *ldu;
#line 381 "IB01QD.f"
    u_offset = 1 + u_dim1;
#line 381 "IB01QD.f"
    u -= u_offset;
#line 381 "IB01QD.f"
    y_dim1 = *ldy;
#line 381 "IB01QD.f"
    y_offset = 1 + y_dim1;
#line 381 "IB01QD.f"
    y -= y_offset;
#line 381 "IB01QD.f"
    --x0;
#line 381 "IB01QD.f"
    b_dim1 = *ldb;
#line 381 "IB01QD.f"
    b_offset = 1 + b_dim1;
#line 381 "IB01QD.f"
    b -= b_offset;
#line 381 "IB01QD.f"
    d_dim1 = *ldd;
#line 381 "IB01QD.f"
    d_offset = 1 + d_dim1;
#line 381 "IB01QD.f"
    d__ -= d_offset;
#line 381 "IB01QD.f"
    --iwork;
#line 381 "IB01QD.f"
    --dwork;
#line 381 "IB01QD.f"

#line 381 "IB01QD.f"
    /* Function Body */
#line 381 "IB01QD.f"
    withd = lsame_(job, "D", (ftnlen)1, (ftnlen)1);
#line 382 "IB01QD.f"
    withb = lsame_(job, "B", (ftnlen)1, (ftnlen)1) || withd;
#line 383 "IB01QD.f"
    withx0 = lsame_(jobx0, "X", (ftnlen)1, (ftnlen)1);

#line 385 "IB01QD.f"
    *iwarn = 0;
#line 386 "IB01QD.f"
    *info = 0;
#line 387 "IB01QD.f"
    lm = *l * *m;
#line 388 "IB01QD.f"
    ln = *l * *n;
#line 389 "IB01QD.f"
    nn = *n * *n;
#line 390 "IB01QD.f"
    nm = *n * *m;
#line 391 "IB01QD.f"
    n2m = *n * nm;
#line 392 "IB01QD.f"
    ncol = nm;
#line 393 "IB01QD.f"
    if (withx0) {
#line 393 "IB01QD.f"
	ncol += *n;
#line 393 "IB01QD.f"
    }
#line 395 "IB01QD.f"
    minsmp = ncol;
#line 396 "IB01QD.f"
    if (withd) {
#line 397 "IB01QD.f"
	minsmp += *m;
#line 398 "IB01QD.f"
	iq = minsmp;
#line 399 "IB01QD.f"
    } else if (! withx0) {
#line 400 "IB01QD.f"
	iq = minsmp;
#line 401 "IB01QD.f"
	++minsmp;
#line 402 "IB01QD.f"
    } else {
#line 403 "IB01QD.f"
	iq = minsmp;
#line 404 "IB01QD.f"
    }

#line 406 "IB01QD.f"
    if (! (withx0 || lsame_(jobx0, "N", (ftnlen)1, (ftnlen)1))) {
#line 407 "IB01QD.f"
	*info = -1;
#line 408 "IB01QD.f"
    } else if (! withb) {
#line 409 "IB01QD.f"
	*info = -2;
#line 410 "IB01QD.f"
    } else if (*n < 0) {
#line 411 "IB01QD.f"
	*info = -3;
#line 412 "IB01QD.f"
    } else if (*m < 0) {
#line 413 "IB01QD.f"
	*info = -4;
#line 414 "IB01QD.f"
    } else if (*l <= 0) {
#line 415 "IB01QD.f"
	*info = -5;
#line 416 "IB01QD.f"
    } else if (*nsmp < minsmp) {
#line 417 "IB01QD.f"
	*info = -6;
#line 418 "IB01QD.f"
    } else if (*lda < max(1,*n)) {
#line 419 "IB01QD.f"
	*info = -8;
#line 420 "IB01QD.f"
    } else if (*ldc < *l) {
#line 421 "IB01QD.f"
	*info = -10;
#line 422 "IB01QD.f"
    } else if (*ldu < 1 || *m > 0 && *ldu < *nsmp) {
#line 423 "IB01QD.f"
	*info = -12;
#line 424 "IB01QD.f"
    } else if (*ldy < max(1,*nsmp)) {
#line 425 "IB01QD.f"
	*info = -14;
#line 426 "IB01QD.f"
    } else if (*ldb < 1 || *ldb < *n && *m > 0) {
#line 428 "IB01QD.f"
	*info = -17;
#line 429 "IB01QD.f"
    } else if (*ldd < 1 || withd && *ldd < *l && *m > 0) {
#line 431 "IB01QD.f"
	*info = -19;
#line 432 "IB01QD.f"
    } else if (*tol > 1.) {
#line 433 "IB01QD.f"
	*info = -20;
#line 434 "IB01QD.f"
    }

/*     Compute workspace. */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

#line 443 "IB01QD.f"
    nsmpl = *nsmp * *l;
#line 444 "IB01QD.f"
    iq *= *l;
#line 445 "IB01QD.f"
    ncp1 = ncol + 1;
#line 446 "IB01QD.f"
    isize = nsmpl * ncp1;
#line 447 "IB01QD.f"
    if (*n > 0 && withx0) {
#line 448 "IB01QD.f"
	ic = (nn << 1) + *n;
#line 449 "IB01QD.f"
    } else {
#line 450 "IB01QD.f"
	ic = 0;
#line 451 "IB01QD.f"
    }
#line 452 "IB01QD.f"
    minwls = ncol * ncp1;
#line 453 "IB01QD.f"
    if (withd) {
#line 453 "IB01QD.f"
	minwls += lm * ncp1;
#line 453 "IB01QD.f"
    }
#line 455 "IB01QD.f"
    if (*m > 0 && withd) {
/* Computing MAX */
#line 456 "IB01QD.f"
	i__1 = ncol << 1;
#line 456 "IB01QD.f"
	ia = *m + max(i__1,*m);
#line 457 "IB01QD.f"
    } else {
#line 458 "IB01QD.f"
	ia = ncol << 1;
#line 459 "IB01QD.f"
    }
#line 460 "IB01QD.f"
    itau = n2m + max(ic,ia);
#line 461 "IB01QD.f"
    if (withx0) {
#line 461 "IB01QD.f"
	itau += ln;
#line 461 "IB01QD.f"
    }
/* Computing MAX */
#line 463 "IB01QD.f"
    i__1 = *n + max(ic,ia), i__2 = ncol * 6;
#line 463 "IB01QD.f"
    ldw2 = isize + max(i__1,i__2);
/* Computing MAX */
#line 464 "IB01QD.f"
    i__1 = iq * ncp1 + itau, i__2 = ncol * 6;
#line 464 "IB01QD.f"
    ldw3 = minwls + max(i__1,i__2);
#line 465 "IB01QD.f"
    if (*m > 0 && withd) {
/* Computing MAX */
#line 466 "IB01QD.f"
	i__1 = ldw2, i__2 = isize + (*m << 1) * *m + *m * 6;
#line 466 "IB01QD.f"
	ldw2 = max(i__1,i__2);
/* Computing MAX */
#line 467 "IB01QD.f"
	i__1 = ldw3, i__2 = minwls + (*m << 1) * *m + *m * 6;
#line 467 "IB01QD.f"
	ldw3 = max(i__1,i__2);
#line 468 "IB01QD.f"
    }
#line 469 "IB01QD.f"
    minwrk = min(ldw2,ldw3);
#line 470 "IB01QD.f"
    minwrk = max(minwrk,2);
#line 471 "IB01QD.f"
    if (*m > 0 && withd) {
#line 471 "IB01QD.f"
	minwrk = max(minwrk,3);
#line 471 "IB01QD.f"
    }
#line 473 "IB01QD.f"
    if (*info == 0 && *ldwork >= minwrk) {
#line 474 "IB01QD.f"
	if (*m > 0 && withd) {
/* Computing MAX */
#line 475 "IB01QD.f"
	    i__3 = *nsmp - *m;
#line 475 "IB01QD.f"
	    i__1 = *m * ilaenv_(&c__1, "DGEQRF", " ", nsmp, m, &c_n1, &c_n1, (
		    ftnlen)6, (ftnlen)1), i__2 = ncol + ncol * ilaenv_(&c__1, 
		    "DGEQRF", " ", &i__3, &ncol, &c_n1, &c_n1, (ftnlen)6, (
		    ftnlen)1);
#line 475 "IB01QD.f"
	    maxwrk = isize + *n + *m + max(i__1,i__2);
/* Computing MAX */
/* Computing MAX */
#line 479 "IB01QD.f"
	    i__5 = *nsmp - *m;
#line 479 "IB01QD.f"
	    i__3 = ncp1 * ilaenv_(&c__1, "DORMQR", "LT", nsmp, &ncp1, m, &
		    c_n1, (ftnlen)6, (ftnlen)2), i__4 = ncol + ilaenv_(&c__1, 
		    "DORMQR", "LT", &i__5, &c__1, &ncol, &c_n1, (ftnlen)6, (
		    ftnlen)2);
#line 479 "IB01QD.f"
	    i__1 = maxwrk, i__2 = isize + *n + *m + max(i__3,i__4);
#line 479 "IB01QD.f"
	    maxwrk = max(i__1,i__2);
#line 484 "IB01QD.f"
	} else {
/* Computing MAX */
#line 485 "IB01QD.f"
	    i__1 = ncol * ilaenv_(&c__1, "DGEQRF", " ", &nsmpl, &ncol, &c_n1, 
		    &c_n1, (ftnlen)6, (ftnlen)1), i__2 = ilaenv_(&c__1, "DOR"\
		    "MQR", "LT", &nsmpl, &c__1, &ncol, &c_n1, (ftnlen)6, (
		    ftnlen)2);
#line 485 "IB01QD.f"
	    maxwrk = isize + *n + ncol + max(i__1,i__2);
#line 490 "IB01QD.f"
	}
#line 491 "IB01QD.f"
	maxwrk = max(maxwrk,minwrk);
#line 492 "IB01QD.f"
    }

#line 494 "IB01QD.f"
    if (*info == 0 && *ldwork < minwrk) {
#line 495 "IB01QD.f"
	*info = -23;
#line 496 "IB01QD.f"
	dwork[1] = (doublereal) minwrk;
#line 497 "IB01QD.f"
    }

/*     Return if there are illegal arguments. */

#line 501 "IB01QD.f"
    if (*info != 0) {
#line 502 "IB01QD.f"
	i__1 = -(*info);
#line 502 "IB01QD.f"
	xerbla_("IB01QD", &i__1, (ftnlen)6);
#line 503 "IB01QD.f"
	return 0;
#line 504 "IB01QD.f"
    }

/*     Quick return if possible. */

#line 508 "IB01QD.f"
    if (max(*n,*m) == 0) {
#line 509 "IB01QD.f"
	dwork[2] = 1.;
#line 510 "IB01QD.f"
	if (*m > 0 && withd) {
#line 511 "IB01QD.f"
	    dwork[1] = 3.;
#line 512 "IB01QD.f"
	    dwork[3] = 1.;
#line 513 "IB01QD.f"
	} else {
#line 514 "IB01QD.f"
	    dwork[1] = 2.;
#line 515 "IB01QD.f"
	}
#line 516 "IB01QD.f"
	return 0;
#line 517 "IB01QD.f"
    }

/*     Set up the least squares problem, either directly, if enough */
/*     workspace, or sequentially, otherwise. */

#line 522 "IB01QD.f"
    iypnt = 1;
#line 523 "IB01QD.f"
    iupnt = 1;
#line 524 "IB01QD.f"
    lddw = (*ldwork - minwls - itau) / ncp1;
/* Computing MIN */
#line 525 "IB01QD.f"
    i__1 = *nsmp, i__2 = lddw / *l;
#line 525 "IB01QD.f"
    nobs = min(i__1,i__2);

#line 527 "IB01QD.f"
    if (*ldwork >= ldw2 || *nsmp <= nobs) {

/*        Enough workspace for solving the problem directly. */

#line 531 "IB01QD.f"
	ncycle = 1;
#line 532 "IB01QD.f"
	nobs = *nsmp;
#line 533 "IB01QD.f"
	lddw = max(1,nsmpl);
#line 534 "IB01QD.f"
	if (withd) {
#line 535 "IB01QD.f"
	    inir = *m + 1;
#line 536 "IB01QD.f"
	} else {
#line 537 "IB01QD.f"
	    inir = 1;
#line 538 "IB01QD.f"
	}
#line 539 "IB01QD.f"
	iny = 1;
#line 540 "IB01QD.f"
	inis = 1;
#line 541 "IB01QD.f"
    } else {

/*        NCYCLE > 1  cycles are needed for solving the problem */
/*        sequentially, taking  NOBS  samples in each cycle (or the */
/*        remaining samples in the last cycle). */

#line 547 "IB01QD.f"
	lnob = *l * nobs;
#line 548 "IB01QD.f"
	lddw = max(1,lnob);
#line 549 "IB01QD.f"
	ncycle = *nsmp / nobs;
#line 550 "IB01QD.f"
	if (*nsmp % nobs != 0) {
#line 550 "IB01QD.f"
	    ++ncycle;
#line 550 "IB01QD.f"
	}
#line 552 "IB01QD.f"
	inir = 1;
#line 553 "IB01QD.f"
	inih = inir + ncol * ncol;
#line 554 "IB01QD.f"
	inis = inih + ncol;
#line 555 "IB01QD.f"
	if (withd) {
#line 556 "IB01QD.f"
	    iny = inis + lm * ncp1;
#line 557 "IB01QD.f"
	} else {
#line 558 "IB01QD.f"
	    iny = inis;
#line 559 "IB01QD.f"
	}
#line 560 "IB01QD.f"
    }

#line 562 "IB01QD.f"
    ncyc = ncycle > 1;
#line 563 "IB01QD.f"
    inygam = iny + lddw * nm;
#line 564 "IB01QD.f"
    irhs = iny + lddw * ncol;
#line 565 "IB01QD.f"
    ixinit = irhs + lddw;
#line 566 "IB01QD.f"
    if (ncyc) {
#line 567 "IB01QD.f"
	ic = ixinit + n2m;
#line 568 "IB01QD.f"
	if (withx0) {
#line 569 "IB01QD.f"
	    ia = ic + ln;
#line 570 "IB01QD.f"
	} else {
#line 571 "IB01QD.f"
	    ia = ic;
#line 572 "IB01QD.f"
	}
#line 573 "IB01QD.f"
	ldr = max(1,ncol);
#line 574 "IB01QD.f"
	ie = iny;
#line 575 "IB01QD.f"
    } else {
#line 576 "IB01QD.f"
	if (withd) {
#line 577 "IB01QD.f"
	    inih = irhs + *m;
#line 578 "IB01QD.f"
	} else {
#line 579 "IB01QD.f"
	    inih = irhs;
#line 580 "IB01QD.f"
	}
#line 581 "IB01QD.f"
	ia = ixinit + *n;
#line 582 "IB01QD.f"
	ldr = lddw;
#line 583 "IB01QD.f"
	ie = ixinit;
#line 584 "IB01QD.f"
    }
#line 585 "IB01QD.f"
    if (*n > 0 && withx0) {
#line 585 "IB01QD.f"
	ias = ia + nn;
#line 585 "IB01QD.f"
    }

#line 588 "IB01QD.f"
    itauu = ia;
#line 589 "IB01QD.f"
    if (withd) {
#line 590 "IB01QD.f"
	itau = itauu + *m;
#line 591 "IB01QD.f"
    } else {
#line 592 "IB01QD.f"
	itau = itauu;
#line 593 "IB01QD.f"
    }
#line 594 "IB01QD.f"
    dum[0] = 0.;

#line 596 "IB01QD.f"
    i__1 = ncycle;
#line 596 "IB01QD.f"
    for (icycle = 1; icycle <= i__1; ++icycle) {
#line 597 "IB01QD.f"
	first = icycle == 1;
#line 598 "IB01QD.f"
	if (! first) {
#line 599 "IB01QD.f"
	    if (icycle == ncycle) {
#line 600 "IB01QD.f"
		nobs = *nsmp - (ncycle - 1) * nobs;
#line 601 "IB01QD.f"
		lnob = *l * nobs;
#line 602 "IB01QD.f"
	    }
#line 603 "IB01QD.f"
	}

#line 605 "IB01QD.f"
	iy = iny;
#line 606 "IB01QD.f"
	ixsave = ixinit;

/*        Compute the  M*N  output trajectories for zero initial state */
/*        or for the saved final state value of the previous cycle. */
/*        This can be performed in parallel. */
/*        Workspace: need  s*L*(r + 1) + b + w, */
/*                   where r = M*N + a,  s = NOBS, */
/*                         a = 0,             if JOBX0 = 'N'; */
/*                         a = N,             if JOBX0 = 'X'; */
/*                         b = N,             if NCYCLE = 1; */
/*                         b = N*N*M,         if NCYCLE > 1; */
/*                         w = 0,             if NCYCLE = 1; */
/*                         w = r*(r+1),       if NCYCLE > 1,  JOB = 'B'; */
/*                         w = (M*L+r)*(r+1), if NCYCLE > 1,  JOB = 'D'. */

#line 621 "IB01QD.f"
	i__2 = *m;
#line 621 "IB01QD.f"
	for (j = 1; j <= i__2; ++j) {
#line 622 "IB01QD.f"
	    i__3 = *n;
#line 622 "IB01QD.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
/*                            ij */
/*              Compute the  y    trajectory and put the vectorized form */
/*              of it in an appropriate column of  DWORK.  To gain in */
/*              efficiency, a specialization of SLICOT Library routine */
/*              TF01ND is used. */

#line 629 "IB01QD.f"
		if (first) {
#line 629 "IB01QD.f"
		    dcopy_(n, dum, &c__0, &dwork[ixsave], &c__1);
#line 629 "IB01QD.f"
		}
#line 631 "IB01QD.f"
		dcopy_(n, &dwork[ixsave], &c__1, &x0[1], &c__1);
#line 632 "IB01QD.f"
		ini = iy;

#line 634 "IB01QD.f"
		i__4 = nobs;
#line 634 "IB01QD.f"
		for (k = 1; k <= i__4; ++k) {
#line 635 "IB01QD.f"
		    dgemv_("No transpose", l, n, &c_b45, &c__[c_offset], ldc, 
			    &x0[1], &c__1, &c_b47, &dwork[iy], &nobs, (ftnlen)
			    12);
#line 637 "IB01QD.f"
		    ++iy;
#line 638 "IB01QD.f"
		    dtrmv_("Upper", "No transpose", "Non-unit", n, &a[
			    a_offset], lda, &x0[1], &c__1, (ftnlen)5, (ftnlen)
			    12, (ftnlen)8);

#line 641 "IB01QD.f"
		    i__5 = *n;
#line 641 "IB01QD.f"
		    for (ix = 2; ix <= i__5; ++ix) {
#line 642 "IB01QD.f"
			x0[ix] += a[ix + (ix - 1) * a_dim1] * dwork[ixsave + 
				ix - 2];
#line 643 "IB01QD.f"
/* L10: */
#line 643 "IB01QD.f"
		    }

#line 645 "IB01QD.f"
		    x0[i__] += u[iupnt + k - 1 + j * u_dim1];
#line 646 "IB01QD.f"
		    dcopy_(n, &x0[1], &c__1, &dwork[ixsave], &c__1);
#line 647 "IB01QD.f"
/* L20: */
#line 647 "IB01QD.f"
		}

#line 649 "IB01QD.f"
		if (ncyc) {
#line 649 "IB01QD.f"
		    ixsave += *n;
#line 649 "IB01QD.f"
		}
#line 651 "IB01QD.f"
		iy = ini + lddw;
#line 652 "IB01QD.f"
/* L30: */
#line 652 "IB01QD.f"
	    }

#line 654 "IB01QD.f"
/* L40: */
#line 654 "IB01QD.f"
	}

#line 656 "IB01QD.f"
	if (*n > 0 && withx0) {

/*           Compute the permuted extended observability matrix  Gamma */
/*                                                                ij */
/*           in the following  N  columns of  DWORK  (after the  y */
/*           trajectories).  Gamma  is directly constructed in the */
/*           required row structure. */
/*           Workspace: need  s*L*(r + 1) + 2*N*N + N + b + c + w, */
/*                      where c = 0,   if NCYCLE = 1; */
/*                            c = L*N, if NCYCLE > 1. */

#line 667 "IB01QD.f"
	    jwork = ias + nn;
#line 668 "IB01QD.f"
	    ig = inygam;
#line 669 "IB01QD.f"
	    iexpon = (integer) (log((doublereal) nobs) / log(2.));
#line 670 "IB01QD.f"
	    irem = nobs - pow_ii(&c__2, &iexpon);
#line 671 "IB01QD.f"
	    power2 = irem == 0;
#line 672 "IB01QD.f"
	    if (! power2) {
#line 672 "IB01QD.f"
		++iexpon;
#line 672 "IB01QD.f"
	    }

#line 675 "IB01QD.f"
	    if (first) {

#line 677 "IB01QD.f"
		i__2 = *n;
#line 677 "IB01QD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 678 "IB01QD.f"
		    dcopy_(l, &c__[i__ * c_dim1 + 1], &c__1, &dwork[ig], &
			    nobs);
#line 679 "IB01QD.f"
		    ig += lddw;
#line 680 "IB01QD.f"
/* L50: */
#line 680 "IB01QD.f"
		}

#line 682 "IB01QD.f"
	    } else {

#line 684 "IB01QD.f"
		i__2 = ic + ln - 1;
#line 684 "IB01QD.f"
		i__3 = *l;
#line 684 "IB01QD.f"
		for (i__ = ic; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += 
			i__3) {
#line 685 "IB01QD.f"
		    dcopy_(l, &dwork[i__], &c__1, &dwork[ig], &nobs);
#line 686 "IB01QD.f"
		    ig += lddw;
#line 687 "IB01QD.f"
/* L60: */
#line 687 "IB01QD.f"
		}

#line 689 "IB01QD.f"
	    }
/*                                          p */
/*           Use powers of the matrix  A:  A ,  p = 2**(J-1). */

#line 693 "IB01QD.f"
	    dlacpy_("Upper", n, n, &a[a_offset], lda, &dwork[ia], n, (ftnlen)
		    5);
#line 694 "IB01QD.f"
	    if (*n > 1) {
#line 694 "IB01QD.f"
		i__3 = *n - 1;
#line 694 "IB01QD.f"
		i__2 = *lda + 1;
#line 694 "IB01QD.f"
		i__4 = *n + 1;
#line 694 "IB01QD.f"
		dcopy_(&i__3, &a[a_dim1 + 2], &i__2, &dwork[ia + 1], &i__4);
#line 694 "IB01QD.f"
	    }
#line 696 "IB01QD.f"
	    i2 = 1;
#line 697 "IB01QD.f"
	    nrow = 0;

#line 699 "IB01QD.f"
	    i__3 = iexpon;
#line 699 "IB01QD.f"
	    for (j = 1; j <= i__3; ++j) {
#line 700 "IB01QD.f"
		igam = inygam;
#line 701 "IB01QD.f"
		if (j < iexpon || power2) {
#line 702 "IB01QD.f"
		    nrow = i2;
#line 703 "IB01QD.f"
		} else {
#line 704 "IB01QD.f"
		    nrow = irem;
#line 705 "IB01QD.f"
		}

#line 707 "IB01QD.f"
		i__2 = *l;
#line 707 "IB01QD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 708 "IB01QD.f"
		    dlacpy_("Full", &nrow, n, &dwork[igam], &lddw, &dwork[
			    igam + i2], &lddw, (ftnlen)4);
#line 710 "IB01QD.f"
		    dtrmm_("Right", "Upper", "No Transpose", "Non Unit", &
			    nrow, n, &c_b45, &dwork[ia], n, &dwork[igam + i2],
			     &lddw, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)
			    8);
#line 713 "IB01QD.f"
		    ig = igam;
/*                                                                  p */
/*                 Compute the contribution of the subdiagonal of  A */
/*                 to the product. */

#line 718 "IB01QD.f"
		    i__4 = *n - 1;
#line 718 "IB01QD.f"
		    for (ix = 1; ix <= i__4; ++ix) {
#line 719 "IB01QD.f"
			daxpy_(&nrow, &dwork[ia + (ix - 1) * *n + ix], &dwork[
				ig + lddw], &c__1, &dwork[ig + i2], &c__1);
#line 721 "IB01QD.f"
			ig += lddw;
#line 722 "IB01QD.f"
/* L70: */
#line 722 "IB01QD.f"
		    }

#line 724 "IB01QD.f"
		    igam += nobs;
#line 725 "IB01QD.f"
/* L80: */
#line 725 "IB01QD.f"
		}

#line 727 "IB01QD.f"
		if (j < iexpon) {
#line 728 "IB01QD.f"
		    dlacpy_("Upper", n, n, &dwork[ia], n, &dwork[ias], n, (
			    ftnlen)5);
#line 730 "IB01QD.f"
		    if (*n > 1) {
#line 730 "IB01QD.f"
			i__2 = *n - 1;
#line 730 "IB01QD.f"
			i__4 = *n + 1;
#line 730 "IB01QD.f"
			i__5 = *n + 1;
#line 730 "IB01QD.f"
			dcopy_(&i__2, &dwork[ia + 1], &i__4, &dwork[ias + 1], 
				&i__5);
#line 730 "IB01QD.f"
		    }
#line 733 "IB01QD.f"
		    mb01td_(n, &dwork[ias], n, &dwork[ia], n, &dwork[jwork], &
			    ierr);
#line 735 "IB01QD.f"
		    i2 <<= 1;
#line 736 "IB01QD.f"
		}
#line 737 "IB01QD.f"
/* L90: */
#line 737 "IB01QD.f"
	    }

#line 739 "IB01QD.f"
	    if (ncyc && icycle < ncycle) {
#line 740 "IB01QD.f"
		ig = inygam + i2 + nrow - 1;
#line 741 "IB01QD.f"
		igs = ig;

#line 743 "IB01QD.f"
		i__3 = ic + ln - 1;
#line 743 "IB01QD.f"
		i__2 = *l;
#line 743 "IB01QD.f"
		for (i__ = ic; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ += 
			i__2) {
#line 744 "IB01QD.f"
		    dcopy_(l, &dwork[ig], &nobs, &dwork[i__], &c__1);
#line 745 "IB01QD.f"
		    ig += lddw;
#line 746 "IB01QD.f"
/* L100: */
#line 746 "IB01QD.f"
		}

#line 748 "IB01QD.f"
		dtrmm_("Right", "Upper", "No Transpose", "Non Unit", l, n, &
			c_b45, &a[a_offset], lda, &dwork[ic], l, (ftnlen)5, (
			ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 750 "IB01QD.f"
		ig = igs;

/*              Compute the contribution of the subdiagonal of  A  to the */
/*              product. */

#line 755 "IB01QD.f"
		i__2 = *n - 1;
#line 755 "IB01QD.f"
		for (ix = 1; ix <= i__2; ++ix) {
#line 756 "IB01QD.f"
		    daxpy_(l, &a[ix + 1 + ix * a_dim1], &dwork[ig + lddw], &
			    nobs, &dwork[ic + (ix - 1) * *l], &c__1);
#line 758 "IB01QD.f"
		    ig += lddw;
#line 759 "IB01QD.f"
/* L110: */
#line 759 "IB01QD.f"
		}

#line 761 "IB01QD.f"
	    }
#line 762 "IB01QD.f"
	}

/*        Setup (part of) the right hand side of the least squares */
/*        problem. */

#line 767 "IB01QD.f"
	iy = irhs;

#line 769 "IB01QD.f"
	i__2 = *l;
#line 769 "IB01QD.f"
	for (k = 1; k <= i__2; ++k) {
#line 770 "IB01QD.f"
	    dcopy_(&nobs, &y[iypnt + k * y_dim1], &c__1, &dwork[iy], &c__1);
#line 771 "IB01QD.f"
	    iy += nobs;
#line 772 "IB01QD.f"
/* L120: */
#line 772 "IB01QD.f"
	}

/*        Compress the data using a special QR factorization. */
/*        Workspace: need   v + y, */
/*                   where  v = s*L*(r + 1) + b + c + w + x, */
/*                          x = M,  y = max( 2*r, M ), */
/*                                             if  JOB = 'D'  and  M > 0, */
/*                          x = 0,  y = 2*r,   if  JOB = 'B'  or   M = 0. */

#line 781 "IB01QD.f"
	if (*m > 0 && withd) {

/*           Case 1:  D  is requested. */

#line 785 "IB01QD.f"
	    jwork = itau;
#line 786 "IB01QD.f"
	    if (first) {
#line 787 "IB01QD.f"
		ini = iny + *m;

/*              Compress the first or single segment of  U,  U1 = Q1*R1. */
/*              Workspace: need   v + M; */
/*                         prefer v + M*NB. */

#line 793 "IB01QD.f"
		i__2 = *ldwork - jwork + 1;
#line 793 "IB01QD.f"
		dgeqrf_(&nobs, m, &u[u_offset], ldu, &dwork[itauu], &dwork[
			jwork], &i__2, &ierr);
/*                                                  ij */
/*              Apply  diag(Q1')  to the matrix  [ y   Gamma Y ]. */
/*              Workspace: need   v + r + 1, */
/*                         prefer v + (r + 1)*NB. */

#line 800 "IB01QD.f"
		i__2 = *l;
#line 800 "IB01QD.f"
		for (k = 1; k <= i__2; ++k) {
#line 801 "IB01QD.f"
		    i__3 = *ldwork - jwork + 1;
#line 801 "IB01QD.f"
		    dormqr_("Left", "Transpose", &nobs, &ncp1, m, &u[u_offset]
			    , ldu, &dwork[itauu], &dwork[iny + (k - 1) * nobs]
			    , &lddw, &dwork[jwork], &i__3, &ierr, (ftnlen)4, (
			    ftnlen)9);
#line 805 "IB01QD.f"
/* L130: */
#line 805 "IB01QD.f"
		}

#line 807 "IB01QD.f"
		if (ncol > 0) {

/*                 Compress the first part of the first data segment of */
/*                    ij */
/*                 [ y   Gamma ]. */
/*                 Workspace: need   v + 2*r, */
/*                            prefer v + r + r*NB. */

#line 815 "IB01QD.f"
		    jwork = itau + ncol;
#line 816 "IB01QD.f"
		    i__2 = nobs - *m;
#line 816 "IB01QD.f"
		    i__3 = *ldwork - jwork + 1;
#line 816 "IB01QD.f"
		    dgeqrf_(&i__2, &ncol, &dwork[ini], &lddw, &dwork[itau], &
			    dwork[jwork], &i__3, &ierr);

/*                 Apply the transformation to the corresponding right */
/*                 hand side part. */
/*                 Workspace: need   v + r + 1, */
/*                            prefer v + r + NB. */

#line 825 "IB01QD.f"
		    i__2 = nobs - *m;
#line 825 "IB01QD.f"
		    i__3 = *ldwork - jwork + 1;
#line 825 "IB01QD.f"
		    dormqr_("Left", "Transpose", &i__2, &c__1, &ncol, &dwork[
			    ini], &lddw, &dwork[itau], &dwork[irhs + *m], &
			    lddw, &dwork[jwork], &i__3, &ierr, (ftnlen)4, (
			    ftnlen)9);

/*                 Compress the remaining parts of the first data segment */
/*                        ij */
/*                 of  [ y   Gamma ]. */
/*                 Workspace: need   v + r - 1. */

#line 835 "IB01QD.f"
		    i__2 = *l;
#line 835 "IB01QD.f"
		    for (k = 2; k <= i__2; ++k) {
#line 836 "IB01QD.f"
			i__3 = nobs - *m;
#line 836 "IB01QD.f"
			mb04od_("Full", &ncol, &c__1, &i__3, &dwork[ini], &
				lddw, &dwork[ini + (k - 1) * nobs], &lddw, &
				dwork[irhs + *m], &lddw, &dwork[irhs + *m + (
				k - 1) * nobs], &lddw, &dwork[itau], &dwork[
				jwork], (ftnlen)4);
#line 841 "IB01QD.f"
/* L140: */
#line 841 "IB01QD.f"
		    }

#line 843 "IB01QD.f"
		}

#line 845 "IB01QD.f"
		if (ncyc) {
/*                                                   ij */
/*                 Save the triangular factor of  [ y   Gamma ],  the */
/*                 corresponding right hand side, and the first  M  rows */
/*                 in each  NOBS  group of rows. */
/*                 Workspace: need   v. */

#line 852 "IB01QD.f"
		    dlacpy_("Upper", &ncol, &ncp1, &dwork[ini], &lddw, &dwork[
			    inir], &ldr, (ftnlen)5);

#line 855 "IB01QD.f"
		    i__2 = *l;
#line 855 "IB01QD.f"
		    for (k = 1; k <= i__2; ++k) {
#line 856 "IB01QD.f"
			dlacpy_("Full", m, &ncp1, &dwork[iny + (k - 1) * nobs]
				, &lddw, &dwork[inis + (k - 1) * *m], &lm, (
				ftnlen)4);
#line 859 "IB01QD.f"
/* L150: */
#line 859 "IB01QD.f"
		    }

#line 861 "IB01QD.f"
		}
#line 862 "IB01QD.f"
	    } else {

/*              Compress the current data segment of  U,  Ui = Qi*Ri, */
/*              i = ICYCLE. */
/*              Workspace: need   v + r + 1. */

#line 868 "IB01QD.f"
		mb04od_("Full", m, &ncp1, &nobs, &u[u_offset], ldu, &u[iupnt 
			+ u_dim1], ldu, &dwork[inis], &lm, &dwork[iny], &lddw,
			 &dwork[itauu], &dwork[jwork], (ftnlen)4);

/*              Apply  diag(Qi')  to the appropriate part of the matrix */
/*                 ij */
/*              [ y   Gamma Y ]. */
/*              Workspace: need   v + r + 1. */

#line 877 "IB01QD.f"
		i__2 = *l;
#line 877 "IB01QD.f"
		for (k = 2; k <= i__2; ++k) {

#line 879 "IB01QD.f"
		    i__3 = *m;
#line 879 "IB01QD.f"
		    for (ix = 1; ix <= i__3; ++ix) {
#line 880 "IB01QD.f"
			mb04oy_(&nobs, &ncp1, &u[iupnt + ix * u_dim1], &dwork[
				itauu + ix - 1], &dwork[inis + (k - 1) * *m + 
				ix - 1], &lm, &dwork[iny + (k - 1) * nobs], &
				lddw, &dwork[jwork]);
#line 885 "IB01QD.f"
/* L160: */
#line 885 "IB01QD.f"
		    }

#line 887 "IB01QD.f"
/* L170: */
#line 887 "IB01QD.f"
		}

#line 889 "IB01QD.f"
		if (ncol > 0) {

#line 891 "IB01QD.f"
		    jwork = itau + ncol;

/*                 Compress the current (but not the first) data segment */
/*                        ij */
/*                 of  [ y   Gamma ]. */
/*                 Workspace: need   v + r - 1. */

#line 898 "IB01QD.f"
		    i__2 = *l;
#line 898 "IB01QD.f"
		    for (k = 1; k <= i__2; ++k) {
#line 899 "IB01QD.f"
			mb04od_("Full", &ncol, &c__1, &nobs, &dwork[inir], &
				ldr, &dwork[iny + (k - 1) * nobs], &lddw, &
				dwork[inih], &ldr, &dwork[irhs + (k - 1) * 
				nobs], &lddw, &dwork[itau], &dwork[jwork], (
				ftnlen)4);
#line 904 "IB01QD.f"
/* L180: */
#line 904 "IB01QD.f"
		    }

#line 906 "IB01QD.f"
		}
#line 907 "IB01QD.f"
	    }

#line 909 "IB01QD.f"
	} else if (ncol > 0) {

/*           Case 2:  D  is known to be zero. */

#line 913 "IB01QD.f"
	    jwork = itau + ncol;
#line 914 "IB01QD.f"
	    if (first) {

/*              Compress the first or single data segment of */
/*                 ij */
/*              [ y   Gamma ]. */
/*              Workspace: need   v + 2*r, */
/*                         prefer v + r + r*NB. */

#line 922 "IB01QD.f"
		i__2 = *ldwork - jwork + 1;
#line 922 "IB01QD.f"
		dgeqrf_(&lddw, &ncol, &dwork[iny], &lddw, &dwork[itau], &
			dwork[jwork], &i__2, &ierr);

/*              Apply the transformation to the right hand side. */
/*              Workspace: need   v + r + 1, */
/*                         prefer v + r + NB. */

#line 929 "IB01QD.f"
		i__2 = *ldwork - jwork + 1;
#line 929 "IB01QD.f"
		dormqr_("Left", "Transpose", &lddw, &c__1, &ncol, &dwork[iny],
			 &lddw, &dwork[itau], &dwork[irhs], &lddw, &dwork[
			jwork], &i__2, &ierr, (ftnlen)4, (ftnlen)9);
#line 932 "IB01QD.f"
		if (ncyc) {
/*                                                   ij */
/*                 Save the triangular factor of  [ y   Gamma ]  and the */
/*                 corresponding right hand side. */
/*                 Workspace: need   v. */

#line 938 "IB01QD.f"
		    dlacpy_("Upper", &ncol, &ncp1, &dwork[iny], &lddw, &dwork[
			    inir], &ldr, (ftnlen)5);
#line 940 "IB01QD.f"
		}
#line 941 "IB01QD.f"
	    } else {

/*              Compress the current (but not the first) data segment. */
/*              Workspace: need   v + r - 1. */

#line 946 "IB01QD.f"
		mb04od_("Full", &ncol, &c__1, &lnob, &dwork[inir], &ldr, &
			dwork[iny], &lddw, &dwork[inih], &ldr, &dwork[irhs], &
			lddw, &dwork[itau], &dwork[jwork], (ftnlen)4);
#line 950 "IB01QD.f"
	    }
#line 951 "IB01QD.f"
	}

#line 953 "IB01QD.f"
	iupnt += nobs;
#line 954 "IB01QD.f"
	iypnt += nobs;
#line 955 "IB01QD.f"
/* L190: */
#line 955 "IB01QD.f"
    }

/*     Estimate the reciprocal condition number of the triangular factor */
/*     of the QR decomposition. */
/*     Workspace: need  u + 3*r, where */
/*                      u = t*L*(r + 1), if NCYCLE = 1; */
/*                      u = w,           if NCYCLE > 1. */

#line 963 "IB01QD.f"
    dtrcon_("1-norm", "Upper", "No Transpose", &ncol, &dwork[inir], &ldr, &
	    rcond, &dwork[ie], &iwork[1], &ierr, (ftnlen)6, (ftnlen)5, (
	    ftnlen)12);

#line 966 "IB01QD.f"
    toll = *tol;
#line 967 "IB01QD.f"
    if (toll <= 0.) {
#line 967 "IB01QD.f"
	toll = dlamch_("Precision", (ftnlen)9);
#line 967 "IB01QD.f"
    }
#line 969 "IB01QD.f"
    if (rcond <= pow_dd(&toll, &c_b113)) {
#line 970 "IB01QD.f"
	*iwarn = 4;

/*        The least squares problem is ill-conditioned. */
/*        Use SVD to solve it. */
/*        Workspace: need   u + 6*r; */
/*                   prefer larger. */

#line 977 "IB01QD.f"
	if (ncol > 1) {
#line 977 "IB01QD.f"
	    i__1 = ncol - 1;
#line 977 "IB01QD.f"
	    i__2 = ncol - 1;
#line 977 "IB01QD.f"
	    dlaset_("Lower", &i__1, &i__2, &c_b47, &c_b47, &dwork[inir + 1], &
		    ldr, (ftnlen)5);
#line 977 "IB01QD.f"
	}
#line 980 "IB01QD.f"
	isv = ie;
#line 981 "IB01QD.f"
	jwork = isv + ncol;
#line 982 "IB01QD.f"
	i__1 = *ldwork - jwork + 1;
#line 982 "IB01QD.f"
	dgelss_(&ncol, &ncol, &c__1, &dwork[inir], &ldr, &dwork[inih], &ldr, &
		dwork[isv], &toll, &rank, &dwork[jwork], &i__1, &ierr);
#line 985 "IB01QD.f"
	if (ierr > 0) {

/*           Return if SVD algorithm did not converge. */

#line 989 "IB01QD.f"
	    *info = 2;
#line 990 "IB01QD.f"
	    return 0;
#line 991 "IB01QD.f"
	}
/* Computing MAX */
#line 992 "IB01QD.f"
	i__1 = maxwrk, i__2 = (integer) dwork[jwork] - jwork + 1;
#line 992 "IB01QD.f"
	maxwrk = max(i__1,i__2);
#line 993 "IB01QD.f"
    } else {

/*        Find the least squares solution using QR decomposition only. */

#line 997 "IB01QD.f"
	dtrsm_("Left", "Upper", "No Transpose", "Non Unit", &ncol, &c__1, &
		c_b45, &dwork[inir], &ldr, &dwork[inih], &ldr, (ftnlen)4, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 999 "IB01QD.f"
    }

/*     Setup the estimated n-by-m input matrix  B,  and the estimated */
/*     initial state of the system  x0. */

#line 1004 "IB01QD.f"
    dlacpy_("Full", n, m, &dwork[inih], n, &b[b_offset], ldb, (ftnlen)4);

#line 1006 "IB01QD.f"
    if (*n > 0 && withx0) {
#line 1007 "IB01QD.f"
	dcopy_(n, &dwork[inih + nm], &c__1, &x0[1], &c__1);
#line 1008 "IB01QD.f"
    } else {
#line 1009 "IB01QD.f"
	dcopy_(n, dum, &c__0, &x0[1], &c__1);
#line 1010 "IB01QD.f"
    }

#line 1012 "IB01QD.f"
    if (*m > 0 && withd) {

/*        Compute the estimated l-by-m input/output matrix  D. */

#line 1016 "IB01QD.f"
	if (ncyc) {
#line 1017 "IB01QD.f"
	    irhs = inis + lm * ncol;
#line 1018 "IB01QD.f"
	    dgemv_("No Transpose", &lm, &ncol, &c_b130, &dwork[inis], &lm, &
		    dwork[inih], &c__1, &c_b45, &dwork[irhs], &c__1, (ftnlen)
		    12);
#line 1020 "IB01QD.f"
	} else {

#line 1022 "IB01QD.f"
	    i__1 = *l;
#line 1022 "IB01QD.f"
	    for (k = 1; k <= i__1; ++k) {
#line 1023 "IB01QD.f"
		dgemv_("No Transpose", m, &ncol, &c_b130, &dwork[inis + (k - 
			1) * nobs], &lddw, &dwork[inih], &c__1, &c_b45, &
			dwork[irhs + (k - 1) * nobs], &c__1, (ftnlen)12);
#line 1026 "IB01QD.f"
/* L200: */
#line 1026 "IB01QD.f"
	    }

#line 1028 "IB01QD.f"
	    i__1 = *l;
#line 1028 "IB01QD.f"
	    for (k = 2; k <= i__1; ++k) {
#line 1029 "IB01QD.f"
		dcopy_(m, &dwork[irhs + (k - 1) * nobs], &c__1, &dwork[irhs + 
			(k - 1) * *m], &c__1);
#line 1031 "IB01QD.f"
/* L210: */
#line 1031 "IB01QD.f"
	    }

#line 1033 "IB01QD.f"
	}

/*        Estimate the reciprocal condition number of the triangular */
/*        factor of the QR decomposition of the matrix U. */
/*        Workspace: need  u + 3*M. */

#line 1039 "IB01QD.f"
	dtrcon_("1-norm", "Upper", "No Transpose", m, &u[u_offset], ldu, &
		rcondu, &dwork[ie], &iwork[1], &ierr, (ftnlen)6, (ftnlen)5, (
		ftnlen)12);
#line 1041 "IB01QD.f"
	if (rcondu <= pow_dd(&toll, &c_b113)) {
#line 1042 "IB01QD.f"
	    *iwarn = 4;

/*           The least squares problem is ill-conditioned. */
/*           Use SVD to solve it. (QR decomposition of  U  is preserved.) */
/*           Workspace: need   u + 2*M*M + 6*M; */
/*                      prefer larger. */

#line 1049 "IB01QD.f"
	    iq = ie + *m * *m;
#line 1050 "IB01QD.f"
	    isv = iq + *m * *m;
#line 1051 "IB01QD.f"
	    jwork = isv + *m;
#line 1052 "IB01QD.f"
	    dlacpy_("Upper", m, m, &u[u_offset], ldu, &dwork[ie], m, (ftnlen)
		    5);
#line 1053 "IB01QD.f"
	    i__1 = *ldwork - jwork + 1;
#line 1053 "IB01QD.f"
	    mb02ud_("Not Factored", "Left", "No Transpose", "No Pinv", m, l, &
		    c_b45, &toll, &rank, &dwork[ie], m, &dwork[iq], m, &dwork[
		    isv], &dwork[irhs], m, dum, &c__1, &dwork[jwork], &i__1, &
		    ierr, (ftnlen)12, (ftnlen)4, (ftnlen)12, (ftnlen)7);
#line 1057 "IB01QD.f"
	    if (ierr > 0) {

/*              Return if SVD algorithm did not converge. */

#line 1061 "IB01QD.f"
		*info = 2;
#line 1062 "IB01QD.f"
		return 0;
#line 1063 "IB01QD.f"
	    }
/* Computing MAX */
#line 1064 "IB01QD.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] - jwork + 1;
#line 1064 "IB01QD.f"
	    maxwrk = max(i__1,i__2);
#line 1065 "IB01QD.f"
	} else {
#line 1066 "IB01QD.f"
	    dtrsm_("Left", "Upper", "No Transpose", "Non Unit", m, l, &c_b45, 
		    &u[u_offset], ldu, &dwork[irhs], m, (ftnlen)4, (ftnlen)5, 
		    (ftnlen)12, (ftnlen)8);
#line 1068 "IB01QD.f"
	}
#line 1069 "IB01QD.f"
	ma02ad_("Full", m, l, &dwork[irhs], m, &d__[d_offset], ldd, (ftnlen)4)
		;

#line 1071 "IB01QD.f"
    }

#line 1073 "IB01QD.f"
    dwork[1] = (doublereal) maxwrk;
#line 1074 "IB01QD.f"
    dwork[2] = rcond;
#line 1075 "IB01QD.f"
    if (*m > 0 && withd) {
#line 1075 "IB01QD.f"
	dwork[3] = rcondu;
#line 1075 "IB01QD.f"
    }

#line 1078 "IB01QD.f"
    return 0;

/* *** End of IB01QD *** */
} /* ib01qd_ */

