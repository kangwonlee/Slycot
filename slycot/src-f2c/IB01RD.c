#line 1 "IB01RD.f"
/* IB01RD.f -- translated by f2c (version 20100827).
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

#line 1 "IB01RD.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static doublereal c_b26 = 1.;
static integer c__0 = 0;
static doublereal c_b49 = -1.;
static doublereal c_b150 = .66666666666666663;
static doublereal c_b152 = 0.;

/* Subroutine */ int ib01rd_(char *job, integer *n, integer *m, integer *l, 
	integer *nsmp, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *u, integer *ldu, doublereal *y, integer *ldy, doublereal *
	x0, doublereal *tol, integer *iwork, doublereal *dwork, integer *
	ldwork, integer *iwarn, integer *info, ftnlen job_len)
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
    static integer j, k, i2, ia, ic, ie, ig, nc, iq, nn, iu, ix, iy, ias, ldr;
    static doublereal dum[1];
    static integer isv, iut, ncp1, ldw1, ldw2, inih, lddw, irem, nrbl, rank;
    static logical ncyc;
    static integer ierr, inir, init, itau, irhs, nobs;
    static doublereal toll;
    static integer nrow;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    mb04od_(char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, ftnlen), 
	    mb01td_(integer *, doublereal *, integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical block;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal rcond;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical withd;
    static integer isize;
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *);
    static logical first;
    static integer nsmpl, jwork;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static integer iupnt;
    extern /* Subroutine */ int dtrsv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static integer iypnt;
    static logical power2;
    extern doublereal dlamch_(char *, ftnlen);
    static integer inigam, icycle;
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
	    integer *), dtrcon_(char *, char *, char *, integer *, doublereal 
	    *, integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static logical switch__;
    static integer iexpon, iutran, ixinit, minsmp, minwrk;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer maxwrk, minwls;


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

/*     To estimate the initial state of a linear time-invariant (LTI) */
/*     discrete-time system, given the system matrices  (A,B,C,D)  and */
/*     the input and output trajectories of the system. The model */
/*     structure is : */

/*           x(k+1) = Ax(k) + Bu(k),   k >= 0, */
/*           y(k)   = Cx(k) + Du(k), */

/*     where  x(k)  is the  n-dimensional state vector (at time k), */
/*            u(k)  is the  m-dimensional input vector, */
/*            y(k)  is the  l-dimensional output vector, */
/*     and  A, B, C, and D  are real matrices of appropriate dimensions. */
/*     Matrix  A  is assumed to be in a real Schur form. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies whether or not the matrix D is zero, as follows: */
/*             = 'Z':  the matrix  D  is zero; */
/*             = 'N':  the matrix  D  is not zero. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     L       (input) INTEGER */
/*             The number of system outputs.  L > 0. */

/*     NSMP    (input) INTEGER */
/*             The number of rows of matrices  U  and  Y  (number of */
/*             samples used,  t).  NSMP >= N. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             system state matrix  A  in a real Schur form. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             system input matrix  B  (corresponding to the real Schur */
/*             form of  A). */
/*             If  N = 0  or  M = 0,  this array is not referenced. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= N,  if  N > 0  and  M > 0; */
/*             LDB >= 1,  if  N = 0  or   M = 0. */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading L-by-N part of this array must contain the */
/*             system output matrix  C  (corresponding to the real Schur */
/*             form of  A). */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= L. */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading L-by-M part of this array must contain the */
/*             system input-output matrix. */
/*             If  M = 0  or  JOB = 'Z',  this array is not referenced. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D. */
/*             LDD >= L,  if  M > 0  and  JOB = 'N'; */
/*             LDD >= 1,  if  M = 0  or   JOB = 'Z'. */

/*     U       (input) DOUBLE PRECISION array, dimension (LDU,M) */
/*             If  M > 0,  the leading NSMP-by-M part of this array must */
/*             contain the t-by-m input-data sequence matrix  U, */
/*             U = [u_1 u_2 ... u_m].  Column  j  of  U  contains the */
/*             NSMP  values of the j-th input component for consecutive */
/*             time increments. */
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
/*             The estimated initial state of the system,  x(0). */

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

/*     IWORK   INTEGER array, dimension (N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = 0,  DWORK(1) returns the optimal value */
/*             of LDWORK and  DWORK(2)  contains the reciprocal condition */
/*             number of the triangular factor of the QR factorization of */
/*             the matrix  Gamma  (see METHOD). */
/*             On exit, if  INFO = -22,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= max( 2, min( LDW1, LDW2 ) ),  where */
/*             LDW1 = t*L*(N + 1) + 2*N + max( 2*N*N, 4*N ), */
/*             LDW2 =   N*(N + 1) + 2*N + */
/*                      max( q*(N + 1) + 2*N*N + L*N, 4*N ), */
/*                q = N*L. */
/*             For good performance,  LDWORK  should be larger. */
/*             If  LDWORK >= LDW1,  then standard QR factorization of */
/*             the matrix  Gamma  (see METHOD) is used. Otherwise, the */
/*             QR factorization is computed sequentially by performing */
/*             NCYCLE  cycles, each cycle (except possibly the last one) */
/*             processing  s  samples, where  s  is chosen by equating */
/*             LDWORK  to  LDW2,  for  q  replaced by  s*L. */
/*             The computational effort may increase and the accuracy may */
/*             decrease with the decrease of  s.  Recommended value is */
/*             LDRWRK = LDW1,  assuming a large enough cache size, to */
/*             also accommodate  A, B, C, D, U,  and  Y. */

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

/*     An extension and refinement of the method in [1] is used. */
/*     Specifically, the output y0(k) of the system for zero initial */
/*     state is computed for k = 0, 1, ...,  t-1 using the given model. */
/*     Then the following least squares problem is solved for x(0) */

/*                         (     C     )            (   y(0) - y0(0)   ) */
/*                         (    C*A    )            (   y(1) - y0(1)   ) */
/*        Gamma * x(0)  =  (     :     ) * x(0)  =  (        :         ). */
/*                         (     :     )            (        :         ) */
/*                         ( C*A^(t-1) )            ( y(t-1) - y0(t-1) ) */

/*     The coefficient matrix  Gamma  is evaluated using powers of A with */
/*     exponents 2^k. The QR decomposition of this matrix is computed. */
/*     If its triangular factor  R  is too ill conditioned, then singular */
/*     value decomposition of  R  is used. */

/*     If the coefficient matrix cannot be stored in the workspace (i.e., */
/*     LDWORK < LDW1),  the QR decomposition is computed sequentially. */

/*     REFERENCES */

/*     [1] Verhaegen M., and Varga, A. */
/*         Some Experience with the MOESP Class of Subspace Model */
/*         Identification Methods in Identifying the BO105 Helicopter. */
/*         Report TR R165-94, DLR Oberpfaffenhofen, 1994. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is numerically stable. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2004. */

/*     KEYWORDS */

/*     Identification methods; least squares solutions; multivariable */
/*     systems; QR decomposition; singular value decomposition. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     IBLOCK is a threshold value for switching to a block algorithm */
/*     for  U  (to avoid row by row passing through  U). */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Check the input parameters. */

#line 274 "IB01RD.f"
    /* Parameter adjustments */
#line 274 "IB01RD.f"
    a_dim1 = *lda;
#line 274 "IB01RD.f"
    a_offset = 1 + a_dim1;
#line 274 "IB01RD.f"
    a -= a_offset;
#line 274 "IB01RD.f"
    b_dim1 = *ldb;
#line 274 "IB01RD.f"
    b_offset = 1 + b_dim1;
#line 274 "IB01RD.f"
    b -= b_offset;
#line 274 "IB01RD.f"
    c_dim1 = *ldc;
#line 274 "IB01RD.f"
    c_offset = 1 + c_dim1;
#line 274 "IB01RD.f"
    c__ -= c_offset;
#line 274 "IB01RD.f"
    d_dim1 = *ldd;
#line 274 "IB01RD.f"
    d_offset = 1 + d_dim1;
#line 274 "IB01RD.f"
    d__ -= d_offset;
#line 274 "IB01RD.f"
    u_dim1 = *ldu;
#line 274 "IB01RD.f"
    u_offset = 1 + u_dim1;
#line 274 "IB01RD.f"
    u -= u_offset;
#line 274 "IB01RD.f"
    y_dim1 = *ldy;
#line 274 "IB01RD.f"
    y_offset = 1 + y_dim1;
#line 274 "IB01RD.f"
    y -= y_offset;
#line 274 "IB01RD.f"
    --x0;
#line 274 "IB01RD.f"
    --iwork;
#line 274 "IB01RD.f"
    --dwork;
#line 274 "IB01RD.f"

#line 274 "IB01RD.f"
    /* Function Body */
#line 274 "IB01RD.f"
    withd = lsame_(job, "N", (ftnlen)1, (ftnlen)1);
#line 275 "IB01RD.f"
    *iwarn = 0;
#line 276 "IB01RD.f"
    *info = 0;
#line 277 "IB01RD.f"
    nn = *n * *n;
#line 278 "IB01RD.f"
    minsmp = *n;

#line 280 "IB01RD.f"
    if (! (lsame_(job, "Z", (ftnlen)1, (ftnlen)1) || withd)) {
#line 281 "IB01RD.f"
	*info = -1;
#line 282 "IB01RD.f"
    } else if (*n < 0) {
#line 283 "IB01RD.f"
	*info = -2;
#line 284 "IB01RD.f"
    } else if (*m < 0) {
#line 285 "IB01RD.f"
	*info = -3;
#line 286 "IB01RD.f"
    } else if (*l <= 0) {
#line 287 "IB01RD.f"
	*info = -4;
#line 288 "IB01RD.f"
    } else if (*nsmp < minsmp) {
#line 289 "IB01RD.f"
	*info = -5;
#line 290 "IB01RD.f"
    } else if (*lda < max(1,*n)) {
#line 291 "IB01RD.f"
	*info = -7;
#line 292 "IB01RD.f"
    } else if (*ldb < 1 || *ldb < *n && *m > 0) {
#line 293 "IB01RD.f"
	*info = -9;
#line 294 "IB01RD.f"
    } else if (*ldc < *l) {
#line 295 "IB01RD.f"
	*info = -11;
#line 296 "IB01RD.f"
    } else if (*ldd < 1 || withd && *ldd < *l && *m > 0) {
#line 298 "IB01RD.f"
	*info = -13;
#line 299 "IB01RD.f"
    } else if (*ldu < 1 || *m > 0 && *ldu < *nsmp) {
#line 300 "IB01RD.f"
	*info = -15;
#line 301 "IB01RD.f"
    } else if (*ldy < max(1,*nsmp)) {
#line 302 "IB01RD.f"
	*info = -17;
#line 303 "IB01RD.f"
    } else if (*tol > 1.) {
#line 304 "IB01RD.f"
	*info = -19;
#line 305 "IB01RD.f"
    }

/*     Compute workspace. */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

#line 314 "IB01RD.f"
    nsmpl = *nsmp * *l;
#line 315 "IB01RD.f"
    iq = minsmp * *l;
#line 316 "IB01RD.f"
    ncp1 = *n + 1;
#line 317 "IB01RD.f"
    isize = nsmpl * ncp1;
#line 318 "IB01RD.f"
    ic = nn << 1;
#line 319 "IB01RD.f"
    minwls = minsmp * ncp1;
#line 320 "IB01RD.f"
    itau = ic + *l * *n;
/* Computing MAX */
#line 321 "IB01RD.f"
    i__1 = ic, i__2 = *n << 2;
#line 321 "IB01RD.f"
    ldw1 = isize + (*n << 1) + max(i__1,i__2);
/* Computing MAX */
#line 322 "IB01RD.f"
    i__1 = iq * ncp1 + itau, i__2 = *n << 2;
#line 322 "IB01RD.f"
    ldw2 = minwls + (*n << 1) + max(i__1,i__2);
/* Computing MAX */
#line 323 "IB01RD.f"
    i__1 = min(ldw1,ldw2);
#line 323 "IB01RD.f"
    minwrk = max(i__1,2);
#line 324 "IB01RD.f"
    if (*info == 0 && *ldwork >= minwrk) {
/* Computing MAX */
#line 325 "IB01RD.f"
	i__1 = *n * ilaenv_(&c__1, "DGEQRF", " ", &nsmpl, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1), i__2 = ilaenv_(&c__1, "DORMQR", "LT", &
		nsmpl, &c__1, n, &c_n1, (ftnlen)6, (ftnlen)2);
#line 325 "IB01RD.f"
	maxwrk = isize + (*n << 1) + max(i__1,i__2);
#line 329 "IB01RD.f"
	maxwrk = max(maxwrk,minwrk);
#line 330 "IB01RD.f"
    }

#line 332 "IB01RD.f"
    if (*info == 0 && *ldwork < minwrk) {
#line 333 "IB01RD.f"
	*info = -22;
#line 334 "IB01RD.f"
	dwork[1] = (doublereal) minwrk;
#line 335 "IB01RD.f"
    }

/*     Return if there are illegal arguments. */

#line 339 "IB01RD.f"
    if (*info != 0) {
#line 340 "IB01RD.f"
	i__1 = -(*info);
#line 340 "IB01RD.f"
	xerbla_("IB01RD", &i__1, (ftnlen)6);
#line 341 "IB01RD.f"
	return 0;
#line 342 "IB01RD.f"
    }

/*     Quick return if possible. */

#line 346 "IB01RD.f"
    if (*n == 0) {
#line 347 "IB01RD.f"
	dwork[1] = 2.;
#line 348 "IB01RD.f"
	dwork[2] = 1.;
#line 349 "IB01RD.f"
	return 0;
#line 350 "IB01RD.f"
    }

/*     Set up the least squares problem, either directly, if enough */
/*     workspace, or sequentially, otherwise. */

#line 355 "IB01RD.f"
    iypnt = 1;
#line 356 "IB01RD.f"
    iupnt = 1;
#line 357 "IB01RD.f"
    inir = 1;
#line 358 "IB01RD.f"
    if (*ldwork >= ldw1) {

/*        Enough workspace for solving the problem directly. */

#line 362 "IB01RD.f"
	ncycle = 1;
#line 363 "IB01RD.f"
	nobs = *nsmp;
#line 364 "IB01RD.f"
	lddw = nsmpl;
#line 365 "IB01RD.f"
	inigam = 1;
#line 366 "IB01RD.f"
    } else {

/*        NCYCLE > 1  cycles are needed for solving the problem */
/*        sequentially, taking  NOBS  samples in each cycle (or the */
/*        remaining samples in the last cycle). */

#line 372 "IB01RD.f"
	jwork = *ldwork - minwls - (*n << 1) - itau;
#line 373 "IB01RD.f"
	lddw = jwork / ncp1;
#line 374 "IB01RD.f"
	nobs = lddw / *l;
#line 375 "IB01RD.f"
	lddw = *l * nobs;
#line 376 "IB01RD.f"
	ncycle = *nsmp / nobs;
#line 377 "IB01RD.f"
	if (*nsmp % nobs != 0) {
#line 377 "IB01RD.f"
	    ++ncycle;
#line 377 "IB01RD.f"
	}
#line 379 "IB01RD.f"
	inih = inir + nn;
#line 380 "IB01RD.f"
	inigam = inih + *n;
#line 381 "IB01RD.f"
    }

#line 383 "IB01RD.f"
    ncyc = ncycle > 1;
#line 384 "IB01RD.f"
    irhs = inigam + lddw * *n;
#line 385 "IB01RD.f"
    ixinit = irhs + lddw;
#line 386 "IB01RD.f"
    ic = ixinit + *n;
#line 387 "IB01RD.f"
    if (ncyc) {
#line 388 "IB01RD.f"
	ia = ic + *l * *n;
#line 389 "IB01RD.f"
	ldr = *n;
#line 390 "IB01RD.f"
	ie = inigam;
#line 391 "IB01RD.f"
    } else {
#line 392 "IB01RD.f"
	inih = irhs;
#line 393 "IB01RD.f"
	ia = ic;
#line 394 "IB01RD.f"
	ldr = lddw;
#line 395 "IB01RD.f"
	ie = ixinit;
#line 396 "IB01RD.f"
    }
#line 397 "IB01RD.f"
    iutran = ia;
#line 398 "IB01RD.f"
    ias = ia + nn;
#line 399 "IB01RD.f"
    itau = ia;
#line 400 "IB01RD.f"
    dum[0] = 0.;

/*     Set block parameters for passing through the array  U. */

#line 404 "IB01RD.f"
    block = *m > 1 && *nsmp * *m >= 16384;
#line 405 "IB01RD.f"
    if (block) {
#line 406 "IB01RD.f"
	nrbl = (*ldwork - iutran + 1) / *m;
#line 407 "IB01RD.f"
	nc = nobs / nrbl;
#line 408 "IB01RD.f"
	if (nobs % nrbl != 0) {
#line 408 "IB01RD.f"
	    ++nc;
#line 408 "IB01RD.f"
	}
#line 410 "IB01RD.f"
	init = (nc - 1) * nrbl;
#line 411 "IB01RD.f"
	block = block && nrbl > 1;
#line 412 "IB01RD.f"
    }

/*     Perform direct of sequential compression of the matrix  Gamma. */

#line 416 "IB01RD.f"
    i__1 = ncycle;
#line 416 "IB01RD.f"
    for (icycle = 1; icycle <= i__1; ++icycle) {
#line 417 "IB01RD.f"
	first = icycle == 1;
#line 418 "IB01RD.f"
	if (! first) {
#line 419 "IB01RD.f"
	    if (icycle == ncycle) {
#line 420 "IB01RD.f"
		nobs = *nsmp - (ncycle - 1) * nobs;
#line 421 "IB01RD.f"
		lddw = *l * nobs;
#line 422 "IB01RD.f"
		if (block) {
#line 423 "IB01RD.f"
		    nc = nobs / nrbl;
#line 424 "IB01RD.f"
		    if (nobs % nrbl != 0) {
#line 424 "IB01RD.f"
			++nc;
#line 424 "IB01RD.f"
		    }
#line 426 "IB01RD.f"
		    init = (nc - 1) * nrbl;
#line 427 "IB01RD.f"
		}
#line 428 "IB01RD.f"
	    }
#line 429 "IB01RD.f"
	}

/*        Compute the extended observability matrix  Gamma. */
/*        Workspace: need   s*L*(N + 1) + 2*N*N + 2*N + a + w, */
/*                   where  s = NOBS, */
/*                          a = 0,   w = 0,          if NCYCLE = 1, */
/*                          a = L*N, w = N*(N + 1),  if NCYCLE > 1; */
/*                   prefer as above, with  s = t,  a = w = 0. */

#line 438 "IB01RD.f"
	jwork = ias + nn;
#line 439 "IB01RD.f"
	iexpon = (integer) (log((doublereal) nobs) / log(2.));
#line 440 "IB01RD.f"
	irem = *l * (nobs - pow_ii(&c__2, &iexpon));
#line 441 "IB01RD.f"
	power2 = irem == 0;
#line 442 "IB01RD.f"
	if (! power2) {
#line 442 "IB01RD.f"
	    ++iexpon;
#line 442 "IB01RD.f"
	}

#line 445 "IB01RD.f"
	if (first) {
#line 446 "IB01RD.f"
	    dlacpy_("Full", l, n, &c__[c_offset], ldc, &dwork[inigam], &lddw, 
		    (ftnlen)4);
#line 447 "IB01RD.f"
	} else {
#line 448 "IB01RD.f"
	    dlacpy_("Full", l, n, &dwork[ic], l, &dwork[inigam], &lddw, (
		    ftnlen)4);
#line 450 "IB01RD.f"
	}
/*                                       p */
/*        Use powers of the matrix  A:  A ,  p = 2**(J-1). */

#line 454 "IB01RD.f"
	dlacpy_("Upper", n, n, &a[a_offset], lda, &dwork[ia], n, (ftnlen)5);
#line 455 "IB01RD.f"
	if (*n > 1) {
#line 455 "IB01RD.f"
	    i__2 = *n - 1;
#line 455 "IB01RD.f"
	    i__3 = *lda + 1;
#line 455 "IB01RD.f"
	    i__4 = *n + 1;
#line 455 "IB01RD.f"
	    dcopy_(&i__2, &a[a_dim1 + 2], &i__3, &dwork[ia + 1], &i__4);
#line 455 "IB01RD.f"
	}
#line 457 "IB01RD.f"
	i2 = *l;
#line 458 "IB01RD.f"
	nrow = 0;

#line 460 "IB01RD.f"
	i__2 = iexpon;
#line 460 "IB01RD.f"
	for (j = 1; j <= i__2; ++j) {
#line 461 "IB01RD.f"
	    ig = inigam;
#line 462 "IB01RD.f"
	    if (j < iexpon || power2) {
#line 463 "IB01RD.f"
		nrow = i2;
#line 464 "IB01RD.f"
	    } else {
#line 465 "IB01RD.f"
		nrow = irem;
#line 466 "IB01RD.f"
	    }

#line 468 "IB01RD.f"
	    dlacpy_("Full", &nrow, n, &dwork[ig], &lddw, &dwork[ig + i2], &
		    lddw, (ftnlen)4);
#line 470 "IB01RD.f"
	    dtrmm_("Right", "Upper", "No Transpose", "Non Unit", &nrow, n, &
		    c_b26, &dwork[ia], n, &dwork[ig + i2], &lddw, (ftnlen)5, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);
/*                                                            p */
/*           Compute the contribution of the subdiagonal of  A   to the */
/*           product. */

#line 477 "IB01RD.f"
	    i__3 = *n - 1;
#line 477 "IB01RD.f"
	    for (ix = 1; ix <= i__3; ++ix) {
#line 478 "IB01RD.f"
		daxpy_(&nrow, &dwork[ia + (ix - 1) * *n + ix], &dwork[ig + 
			lddw], &c__1, &dwork[ig + i2], &c__1);
#line 480 "IB01RD.f"
		ig += lddw;
#line 481 "IB01RD.f"
/* L10: */
#line 481 "IB01RD.f"
	    }

#line 483 "IB01RD.f"
	    if (j < iexpon) {
#line 484 "IB01RD.f"
		dlacpy_("Upper", n, n, &dwork[ia], n, &dwork[ias], n, (ftnlen)
			5);
#line 485 "IB01RD.f"
		i__3 = *n - 1;
#line 485 "IB01RD.f"
		i__4 = *n + 1;
#line 485 "IB01RD.f"
		i__5 = *n + 1;
#line 485 "IB01RD.f"
		dcopy_(&i__3, &dwork[ia + 1], &i__4, &dwork[ias + 1], &i__5);
#line 486 "IB01RD.f"
		mb01td_(n, &dwork[ias], n, &dwork[ia], n, &dwork[jwork], &
			ierr);
#line 488 "IB01RD.f"
		i2 <<= 1;
#line 489 "IB01RD.f"
	    }
#line 490 "IB01RD.f"
/* L20: */
#line 490 "IB01RD.f"
	}

#line 492 "IB01RD.f"
	if (ncyc) {
#line 493 "IB01RD.f"
	    ig = inigam + i2 + nrow - *l;
#line 494 "IB01RD.f"
	    dlacpy_("Full", l, n, &dwork[ig], &lddw, &dwork[ic], l, (ftnlen)4)
		    ;
#line 495 "IB01RD.f"
	    dtrmm_("Right", "Upper", "No Transpose", "Non Unit", l, n, &c_b26,
		     &a[a_offset], lda, &dwork[ic], l, (ftnlen)5, (ftnlen)5, (
		    ftnlen)12, (ftnlen)8);

/*           Compute the contribution of the subdiagonal of  A  to the */
/*           product. */

#line 501 "IB01RD.f"
	    i__2 = *n - 1;
#line 501 "IB01RD.f"
	    for (ix = 1; ix <= i__2; ++ix) {
#line 502 "IB01RD.f"
		daxpy_(l, &a[ix + 1 + ix * a_dim1], &dwork[ig + lddw], &c__1, 
			&dwork[ic + (ix - 1) * *l], &c__1);
#line 504 "IB01RD.f"
		ig += lddw;
#line 505 "IB01RD.f"
/* L30: */
#line 505 "IB01RD.f"
	    }

#line 507 "IB01RD.f"
	}

/*        Setup (part of) the right hand side of the least squares */
/*        problem starting from  DWORK(IRHS);  use the estimated output */
/*        trajectory for zero initial state, or for the saved final state */
/*        value of the previous cycle. */
/*        A specialization of SLICOT Library routine TF01ND is used. */
/*        For large input sets  (NSMP*M >= IBLOCK),  chunks of  U  are */
/*        transposed, to reduce the number of row-wise passes. */
/*        Workspace: need   s*L*(N + 1) + N + w; */
/*                   prefer as above, with  s = t,  w = 0. */

#line 519 "IB01RD.f"
	if (first) {
#line 519 "IB01RD.f"
	    dcopy_(n, dum, &c__0, &dwork[ixinit], &c__1);
#line 519 "IB01RD.f"
	}
#line 521 "IB01RD.f"
	dcopy_(n, &dwork[ixinit], &c__1, &x0[1], &c__1);
#line 522 "IB01RD.f"
	iy = irhs;

#line 524 "IB01RD.f"
	i__2 = *l;
#line 524 "IB01RD.f"
	for (j = 1; j <= i__2; ++j) {
#line 525 "IB01RD.f"
	    dcopy_(&nobs, &y[iypnt + j * y_dim1], &c__1, &dwork[iy], l);
#line 526 "IB01RD.f"
	    ++iy;
#line 527 "IB01RD.f"
/* L40: */
#line 527 "IB01RD.f"
	}

#line 529 "IB01RD.f"
	iy = irhs;
#line 530 "IB01RD.f"
	iu = iupnt;
#line 531 "IB01RD.f"
	if (*m > 0) {
#line 532 "IB01RD.f"
	    if (withd) {

#line 534 "IB01RD.f"
		if (block) {
#line 535 "IB01RD.f"
		    switch__ = TRUE_;
#line 536 "IB01RD.f"
		    nrow = nrbl;

#line 538 "IB01RD.f"
		    i__2 = nobs;
#line 538 "IB01RD.f"
		    for (k = 1; k <= i__2; ++k) {
#line 539 "IB01RD.f"
			if ((k - 1) % nrow == 0 && switch__) {
#line 540 "IB01RD.f"
			    iut = iutran;
#line 541 "IB01RD.f"
			    if (k > init) {
#line 542 "IB01RD.f"
				nrow = nobs - init;
#line 543 "IB01RD.f"
				switch__ = FALSE_;
#line 544 "IB01RD.f"
			    }
#line 545 "IB01RD.f"
			    ma02ad_("Full", &nrow, m, &u[iu + u_dim1], ldu, &
				    dwork[iut], m, (ftnlen)4);
#line 547 "IB01RD.f"
			    iu += nrow;
#line 548 "IB01RD.f"
			}
#line 549 "IB01RD.f"
			dgemv_("No transpose", l, n, &c_b49, &c__[c_offset], 
				ldc, &x0[1], &c__1, &c_b26, &dwork[iy], &c__1,
				 (ftnlen)12);
#line 551 "IB01RD.f"
			dgemv_("No transpose", l, m, &c_b49, &d__[d_offset], 
				ldd, &dwork[iut], &c__1, &c_b26, &dwork[iy], &
				c__1, (ftnlen)12);
#line 553 "IB01RD.f"
			dtrmv_("Upper", "No transpose", "Non-unit", n, &a[
				a_offset], lda, &x0[1], &c__1, (ftnlen)5, (
				ftnlen)12, (ftnlen)8);

#line 556 "IB01RD.f"
			i__3 = *n;
#line 556 "IB01RD.f"
			for (ix = 2; ix <= i__3; ++ix) {
#line 557 "IB01RD.f"
			    x0[ix] += a[ix + (ix - 1) * a_dim1] * dwork[
				    ixinit + ix - 2];
#line 558 "IB01RD.f"
/* L50: */
#line 558 "IB01RD.f"
			}

#line 560 "IB01RD.f"
			dgemv_("No transpose", n, m, &c_b26, &b[b_offset], 
				ldb, &dwork[iut], &c__1, &c_b26, &x0[1], &
				c__1, (ftnlen)12);
#line 562 "IB01RD.f"
			dcopy_(n, &x0[1], &c__1, &dwork[ixinit], &c__1);
#line 563 "IB01RD.f"
			iy += *l;
#line 564 "IB01RD.f"
			iut += *m;
#line 565 "IB01RD.f"
/* L60: */
#line 565 "IB01RD.f"
		    }

#line 567 "IB01RD.f"
		} else {

#line 569 "IB01RD.f"
		    i__2 = nobs;
#line 569 "IB01RD.f"
		    for (k = 1; k <= i__2; ++k) {
#line 570 "IB01RD.f"
			dgemv_("No transpose", l, n, &c_b49, &c__[c_offset], 
				ldc, &x0[1], &c__1, &c_b26, &dwork[iy], &c__1,
				 (ftnlen)12);
#line 572 "IB01RD.f"
			dgemv_("No transpose", l, m, &c_b49, &d__[d_offset], 
				ldd, &u[iu + u_dim1], ldu, &c_b26, &dwork[iy],
				 &c__1, (ftnlen)12);
#line 574 "IB01RD.f"
			dtrmv_("Upper", "No transpose", "Non-unit", n, &a[
				a_offset], lda, &x0[1], &c__1, (ftnlen)5, (
				ftnlen)12, (ftnlen)8);

#line 577 "IB01RD.f"
			i__3 = *n;
#line 577 "IB01RD.f"
			for (ix = 2; ix <= i__3; ++ix) {
#line 578 "IB01RD.f"
			    x0[ix] += a[ix + (ix - 1) * a_dim1] * dwork[
				    ixinit + ix - 2];
#line 579 "IB01RD.f"
/* L70: */
#line 579 "IB01RD.f"
			}

#line 581 "IB01RD.f"
			dgemv_("No transpose", n, m, &c_b26, &b[b_offset], 
				ldb, &u[iu + u_dim1], ldu, &c_b26, &x0[1], &
				c__1, (ftnlen)12);
#line 583 "IB01RD.f"
			dcopy_(n, &x0[1], &c__1, &dwork[ixinit], &c__1);
#line 584 "IB01RD.f"
			iy += *l;
#line 585 "IB01RD.f"
			++iu;
#line 586 "IB01RD.f"
/* L80: */
#line 586 "IB01RD.f"
		    }

#line 588 "IB01RD.f"
		}

#line 590 "IB01RD.f"
	    } else {

#line 592 "IB01RD.f"
		if (block) {
#line 593 "IB01RD.f"
		    switch__ = TRUE_;
#line 594 "IB01RD.f"
		    nrow = nrbl;

#line 596 "IB01RD.f"
		    i__2 = nobs;
#line 596 "IB01RD.f"
		    for (k = 1; k <= i__2; ++k) {
#line 597 "IB01RD.f"
			if ((k - 1) % nrow == 0 && switch__) {
#line 598 "IB01RD.f"
			    iut = iutran;
#line 599 "IB01RD.f"
			    if (k > init) {
#line 600 "IB01RD.f"
				nrow = nobs - init;
#line 601 "IB01RD.f"
				switch__ = FALSE_;
#line 602 "IB01RD.f"
			    }
#line 603 "IB01RD.f"
			    ma02ad_("Full", &nrow, m, &u[iu + u_dim1], ldu, &
				    dwork[iut], m, (ftnlen)4);
#line 605 "IB01RD.f"
			    iu += nrow;
#line 606 "IB01RD.f"
			}
#line 607 "IB01RD.f"
			dgemv_("No transpose", l, n, &c_b49, &c__[c_offset], 
				ldc, &x0[1], &c__1, &c_b26, &dwork[iy], &c__1,
				 (ftnlen)12);
#line 609 "IB01RD.f"
			dtrmv_("Upper", "No transpose", "Non-unit", n, &a[
				a_offset], lda, &x0[1], &c__1, (ftnlen)5, (
				ftnlen)12, (ftnlen)8);

#line 612 "IB01RD.f"
			i__3 = *n;
#line 612 "IB01RD.f"
			for (ix = 2; ix <= i__3; ++ix) {
#line 613 "IB01RD.f"
			    x0[ix] += a[ix + (ix - 1) * a_dim1] * dwork[
				    ixinit + ix - 2];
#line 614 "IB01RD.f"
/* L90: */
#line 614 "IB01RD.f"
			}

#line 616 "IB01RD.f"
			dgemv_("No transpose", n, m, &c_b26, &b[b_offset], 
				ldb, &dwork[iut], &c__1, &c_b26, &x0[1], &
				c__1, (ftnlen)12);
#line 618 "IB01RD.f"
			dcopy_(n, &x0[1], &c__1, &dwork[ixinit], &c__1);
#line 619 "IB01RD.f"
			iy += *l;
#line 620 "IB01RD.f"
			iut += *m;
#line 621 "IB01RD.f"
/* L100: */
#line 621 "IB01RD.f"
		    }

#line 623 "IB01RD.f"
		} else {

#line 625 "IB01RD.f"
		    i__2 = nobs;
#line 625 "IB01RD.f"
		    for (k = 1; k <= i__2; ++k) {
#line 626 "IB01RD.f"
			dgemv_("No transpose", l, n, &c_b49, &c__[c_offset], 
				ldc, &x0[1], &c__1, &c_b26, &dwork[iy], &c__1,
				 (ftnlen)12);
#line 628 "IB01RD.f"
			dtrmv_("Upper", "No transpose", "Non-unit", n, &a[
				a_offset], lda, &x0[1], &c__1, (ftnlen)5, (
				ftnlen)12, (ftnlen)8);

#line 631 "IB01RD.f"
			i__3 = *n;
#line 631 "IB01RD.f"
			for (ix = 2; ix <= i__3; ++ix) {
#line 632 "IB01RD.f"
			    x0[ix] += a[ix + (ix - 1) * a_dim1] * dwork[
				    ixinit + ix - 2];
#line 633 "IB01RD.f"
/* L110: */
#line 633 "IB01RD.f"
			}

#line 635 "IB01RD.f"
			dgemv_("No transpose", n, m, &c_b26, &b[b_offset], 
				ldb, &u[iu + u_dim1], ldu, &c_b26, &x0[1], &
				c__1, (ftnlen)12);
#line 637 "IB01RD.f"
			dcopy_(n, &x0[1], &c__1, &dwork[ixinit], &c__1);
#line 638 "IB01RD.f"
			iy += *l;
#line 639 "IB01RD.f"
			++iu;
#line 640 "IB01RD.f"
/* L120: */
#line 640 "IB01RD.f"
		    }

#line 642 "IB01RD.f"
		}

#line 644 "IB01RD.f"
	    }

#line 646 "IB01RD.f"
	} else {

#line 648 "IB01RD.f"
	    i__2 = nobs;
#line 648 "IB01RD.f"
	    for (k = 1; k <= i__2; ++k) {
#line 649 "IB01RD.f"
		dgemv_("No transpose", l, n, &c_b49, &c__[c_offset], ldc, &x0[
			1], &c__1, &c_b26, &dwork[iy], &c__1, (ftnlen)12);
#line 651 "IB01RD.f"
		dtrmv_("Upper", "No transpose", "Non-unit", n, &a[a_offset], 
			lda, &x0[1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);

#line 654 "IB01RD.f"
		i__3 = *n;
#line 654 "IB01RD.f"
		for (ix = 2; ix <= i__3; ++ix) {
#line 655 "IB01RD.f"
		    x0[ix] += a[ix + (ix - 1) * a_dim1] * dwork[ixinit + ix - 
			    2];
#line 656 "IB01RD.f"
/* L130: */
#line 656 "IB01RD.f"
		}

#line 658 "IB01RD.f"
		dcopy_(n, &x0[1], &c__1, &dwork[ixinit], &c__1);
#line 659 "IB01RD.f"
		iy += *l;
#line 660 "IB01RD.f"
/* L140: */
#line 660 "IB01RD.f"
	    }

#line 662 "IB01RD.f"
	}

/*        Compress the data using (sequential) QR factorization. */
/*        Workspace: need   v + 2*N; */
/*                   where  v = s*L*(N + 1) + N + a + w. */

#line 668 "IB01RD.f"
	jwork = itau + *n;
#line 669 "IB01RD.f"
	if (first) {

/*           Compress the first data segment of  Gamma. */
/*           Workspace: need   v + 2*N, */
/*                      prefer v + N + N*NB. */

#line 675 "IB01RD.f"
	    i__2 = *ldwork - jwork + 1;
#line 675 "IB01RD.f"
	    dgeqrf_(&lddw, n, &dwork[inigam], &lddw, &dwork[itau], &dwork[
		    jwork], &i__2, &ierr);

/*           Apply the transformation to the right hand side part. */
/*           Workspace: need   v + N + 1, */
/*                      prefer v + N + NB. */

#line 682 "IB01RD.f"
	    i__2 = *ldwork - jwork + 1;
#line 682 "IB01RD.f"
	    dormqr_("Left", "Transpose", &lddw, &c__1, n, &dwork[inigam], &
		    lddw, &dwork[itau], &dwork[irhs], &lddw, &dwork[jwork], &
		    i__2, &ierr, (ftnlen)4, (ftnlen)9);

#line 686 "IB01RD.f"
	    if (ncyc) {

/*              Save the triangular factor of  Gamma  and the */
/*              corresponding right hand side. */

#line 691 "IB01RD.f"
		dlacpy_("Upper", n, &ncp1, &dwork[inigam], &lddw, &dwork[inir]
			, &ldr, (ftnlen)5);
#line 693 "IB01RD.f"
	    }
#line 694 "IB01RD.f"
	} else {

/*           Compress the current (but not the first) data segment of */
/*           Gamma. */
/*           Workspace: need   v + N - 1. */

#line 700 "IB01RD.f"
	    mb04od_("Full", n, &c__1, &lddw, &dwork[inir], &ldr, &dwork[
		    inigam], &lddw, &dwork[inih], &ldr, &dwork[irhs], &lddw, &
		    dwork[itau], &dwork[jwork], (ftnlen)4);
#line 703 "IB01RD.f"
	}

#line 705 "IB01RD.f"
	iupnt += nobs;
#line 706 "IB01RD.f"
	iypnt += nobs;
#line 707 "IB01RD.f"
/* L150: */
#line 707 "IB01RD.f"
    }

/*     Estimate the reciprocal condition number of the triangular factor */
/*     of the QR decomposition. */
/*     Workspace: need  u + 3*N, where */
/*                      u = t*L*(N + 1), if NCYCLE = 1; */
/*                      u = w,           if NCYCLE > 1. */

#line 715 "IB01RD.f"
    dtrcon_("1-norm", "Upper", "No Transpose", n, &dwork[inir], &ldr, &rcond, 
	    &dwork[ie], &iwork[1], &ierr, (ftnlen)6, (ftnlen)5, (ftnlen)12);

#line 718 "IB01RD.f"
    toll = *tol;
#line 719 "IB01RD.f"
    if (toll <= 0.) {
#line 719 "IB01RD.f"
	toll = dlamch_("Precision", (ftnlen)9);
#line 719 "IB01RD.f"
    }
#line 721 "IB01RD.f"
    if (rcond <= pow_dd(&toll, &c_b150)) {
#line 722 "IB01RD.f"
	*iwarn = 4;

/*        The least squares problem is ill-conditioned. */
/*        Use SVD to solve it. */
/*        Workspace: need   u + 6*N; */
/*                   prefer larger. */

#line 729 "IB01RD.f"
	i__1 = *n - 1;
#line 729 "IB01RD.f"
	i__2 = *n - 1;
#line 729 "IB01RD.f"
	dlaset_("Lower", &i__1, &i__2, &c_b152, &c_b152, &dwork[inir + 1], &
		ldr, (ftnlen)5);
#line 731 "IB01RD.f"
	isv = ie;
#line 732 "IB01RD.f"
	jwork = isv + *n;
#line 733 "IB01RD.f"
	i__1 = *ldwork - jwork + 1;
#line 733 "IB01RD.f"
	dgelss_(n, n, &c__1, &dwork[inir], &ldr, &dwork[inih], &ldr, &dwork[
		isv], &toll, &rank, &dwork[jwork], &i__1, &ierr);
#line 736 "IB01RD.f"
	if (ierr > 0) {

/*           Return if SVD algorithm did not converge. */

#line 740 "IB01RD.f"
	    *info = 2;
#line 741 "IB01RD.f"
	    return 0;
#line 742 "IB01RD.f"
	}
/* Computing MAX */
#line 743 "IB01RD.f"
	i__1 = maxwrk, i__2 = (integer) dwork[jwork] - jwork + 1;
#line 743 "IB01RD.f"
	maxwrk = max(i__1,i__2);
#line 744 "IB01RD.f"
    } else {

/*        Find the least squares solution using QR decomposition only. */

#line 748 "IB01RD.f"
	dtrsv_("Upper", "No Transpose", "Non Unit", n, &dwork[inir], &ldr, &
		dwork[inih], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 750 "IB01RD.f"
    }

/*     Return the estimated initial state of the system  x0. */

#line 754 "IB01RD.f"
    dcopy_(n, &dwork[inih], &c__1, &x0[1], &c__1);

#line 756 "IB01RD.f"
    dwork[1] = (doublereal) maxwrk;
#line 757 "IB01RD.f"
    dwork[2] = rcond;

#line 759 "IB01RD.f"
    return 0;

/* *** End of IB01RD *** */
} /* ib01rd_ */

