#line 1 "FB01TD.f"
/* FB01TD.f -- translated by f2c (version 20100827).
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

#line 1 "FB01TD.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b23 = 1.;

/* Subroutine */ int fb01td_(char *jobx, char *multrc, integer *n, integer *m,
	 integer *p, doublereal *sinv, integer *ldsinv, doublereal *ainv, 
	integer *ldainv, doublereal *ainvb, integer *ldainb, doublereal *rinv,
	 integer *ldrinv, doublereal *c__, integer *ldc, doublereal *qinv, 
	integer *ldqinv, doublereal *x, doublereal *rinvy, doublereal *z__, 
	doublereal *e, doublereal *tol, integer *iwork, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen jobx_len, ftnlen multrc_len)
{
    /* System generated locals */
    integer ainv_dim1, ainv_offset, ainvb_dim1, ainvb_offset, c_dim1, 
	    c_offset, qinv_dim1, qinv_offset, rinv_dim1, rinv_offset, 
	    sinv_dim1, sinv_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Local variables */
    static integer i__, m1, n1, i12, i13, i23, i32, i33, ii, ij, nm, np, mp1, 
	    ldw;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer itau;
    extern /* Subroutine */ int mb04id_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *), mb04kd_(char *,
	     integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, ftnlen), mb02od_(char *, 
	    char *, char *, char *, char *, integer *, integer *, doublereal *
	    , doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal rcond;
    static logical ljobx;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer jwork;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    static logical lmultr;
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

/*     To calculate a combined measurement and time update of one */
/*     iteration of the time-invariant Kalman filter. This update is */
/*     given for the square root information filter, using the condensed */
/*     controller Hessenberg form. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBX    CHARACTER*1 */
/*             Indicates whether X    is to be computed as follows: */
/*                                i+1 */
/*             = 'X':  X    is computed and stored in array X; */
/*                      i+1 */
/*             = 'N':  X    is not required. */
/*                      i+1 */

/*     MULTRC  CHARACTER*1             -1/2 */
/*             Indicates how matrices R     and C    are to be passed to */
/*                                     i+1       i+1 */
/*             the routine as follows: */
/*             = 'P':  Array RINV is not used and the array C must */
/*                                          -1/2 */
/*                     contain the product R    C   ; */
/*                                          i+1  i+1 */
/*             = 'N':  Arrays RINV and C must contain the matrices */
/*                     as described below. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The actual state dimension, i.e., the order of the */
/*                       -1      -1 */
/*             matrices S   and A  .  N >= 0. */
/*                       i */

/*     M       (input) INTEGER */
/*             The actual input dimension, i.e., the order of the matrix */
/*              -1/2 */
/*             Q    .  M >= 0. */
/*              i */

/*     P       (input) INTEGER */
/*             The actual output dimension, i.e., the order of the matrix */
/*              -1/2 */
/*             R    .  P >= 0. */
/*              i+1 */

/*     SINV    (input/output) DOUBLE PRECISION array, dimension */
/*             (LDSINV,N) */
/*             On entry, the leading N-by-N upper triangular part of this */
/*                                 -1 */
/*             array must contain S  , the inverse of the square root */
/*                                 i */
/*             (right Cholesky factor) of the state covariance matrix */
/*             P    (hence the information square root) at instant i. */
/*              i|i */
/*             On exit, the leading N-by-N upper triangular part of this */
/*                             -1 */
/*             array contains S   , the inverse of the square root (right */
/*                             i+1 */
/*             Cholesky factor) of the state covariance matrix P */
/*                                                              i+1|i+1 */
/*             (hence the information square root) at instant i+1. */
/*             The strict lower triangular part of this array is not */
/*             referenced. */

/*     LDSINV  INTEGER */
/*             The leading dimension of array SINV.  LDSINV >= MAX(1,N). */

/*     AINV    (input) DOUBLE PRECISION array, dimension (LDAINV,N) */
/*                                                                 -1 */
/*             The leading N-by-N part of this array must contain A  , */
/*             the inverse of the state transition matrix of the discrete */
/*             system in controller Hessenberg form (e.g., as produced by */
/*             SLICOT Library Routine TB01MD). */

/*     LDAINV  INTEGER */
/*             The leading dimension of array AINV.  LDAINV >= MAX(1,N). */

/*     AINVB   (input) DOUBLE PRECISION array, dimension (LDAINB,M) */
/*                                                                  -1 */
/*             The leading N-by-M part of this array must contain  A  B, */
/*                             -1 */
/*             the product of A   and the input weight matrix B of the */
/*             discrete system, in upper controller Hessenberg form */
/*             (e.g., as produced by SLICOT Library Routine TB01MD). */

/*     LDAINB  INTEGER */
/*             The leading dimension of array AINVB.  LDAINB >= MAX(1,N). */

/*     RINV    (input) DOUBLE PRECISION array, dimension (LDRINV,*) */
/*             If MULTRC = 'N', then the leading P-by-P upper triangular */
/*                                              -1/2 */
/*             part of this array must contain R    , the inverse of the */
/*                                              i+1 */
/*             covariance square root (right Cholesky factor) of the */
/*             output (measurement) noise (hence the information square */
/*             root) at instant i+1. */
/*             The strict lower triangular part of this array is not */
/*             referenced. */
/*             Otherwise, RINV is not referenced and can be supplied as a */
/*             dummy array (i.e., set parameter LDRINV = 1 and declare */
/*             this array to be RINV(1,1) in the calling program). */

/*     LDRINV  INTEGER */
/*             The leading dimension of array RINV. */
/*             LDRINV >= MAX(1,P) if MULTRC = 'N'; */
/*             LDRINV >= 1        if MULTRC = 'P'. */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array must contain C   , */
/*                                                       -1/2      i+1 */
/*             the output weight matrix (or the product R    C    if */
/*                                                       i+1  i+1 */
/*             MULTRC = 'P') of the discrete system at instant i+1. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     QINV    (input/output) DOUBLE PRECISION array, dimension */
/*             (LDQINV,M) */
/*             On entry, the leading M-by-M upper triangular part of this */
/*                                 -1/2 */
/*             array must contain Q    , the inverse of the covariance */
/*                                 i */
/*             square root (right Cholesky factor) of the input (process) */
/*             noise (hence the information square root) at instant i. */
/*             On exit, the leading M-by-M upper triangular part of this */
/*                                    -1/2 */
/*             array contains (QINOV )    , the inverse of the covariance */
/*                                  i */
/*             square root (right Cholesky factor) of the process noise */
/*             innovation (hence the information square root) at */
/*             instant i. */
/*             The strict lower triangular part of this array is not */
/*             referenced. */

/*     LDQINV  INTEGER */
/*             The leading dimension of array QINV.  LDQINV >= MAX(1,M). */

/*     X       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain X , the estimated */
/*                                                i */
/*             filtered state at instant i. */
/*             On exit, if JOBX = 'X', and INFO = 0, then this array */
/*             contains X   , the estimated filtered state at */
/*                       i+1 */
/*             instant i+1. */
/*             On exit, if JOBX = 'N', or JOBX = 'X' and INFO = 1, then */
/*                                  -1 */
/*             this array contains S   X   . */
/*                                  i+1 i+1 */

/*     RINVY   (input) DOUBLE PRECISION array, dimension (P) */
/*                                      -1/2 */
/*             This array must contain R    Y   , the product of the */
/*                                      i+1  i+1 */
/*                                      -1/2 */
/*             upper triangular matrix R     and the measured output */
/*                                      i+1 */
/*             vector Y    at instant i+1. */
/*                     i+1 */

/*     Z       (input) DOUBLE PRECISION array, dimension (M) */
/*             This array must contain Z , the mean value of the state */
/*                                      i */
/*             process noise at instant i. */

/*     E       (output) DOUBLE PRECISION array, dimension (P) */
/*             This array contains E   , the estimated error at instant */
/*                                  i+1 */
/*             i+1. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If JOBX = 'X', then TOL is used to test for near */
/*                                        -1 */
/*             singularity of the matrix S   . If the user sets */
/*                                        i+1 */
/*             TOL > 0, then the given value of TOL is used as a */
/*             lower bound for the reciprocal condition number of that */
/*             matrix; a matrix whose estimated condition number is less */
/*             than 1/TOL is considered to be nonsingular. If the user */
/*             sets TOL <= 0, then an implicitly computed, default */
/*             tolerance, defined by TOLDEF = N*N*EPS, is used instead, */
/*             where EPS is the machine precision (see LAPACK Library */
/*             routine DLAMCH). */
/*             Otherwise, TOL is not referenced. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             where LIWORK = N if JOBX = 'X', */
/*             and   LIWORK = 1 otherwise. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK.  If INFO = 0 and JOBX = 'X', DWORK(2) returns */
/*             an estimate of the reciprocal of the condition number */
/*                                 -1 */
/*             (in the 1-norm) of S   . */
/*                                 i+1 */

/*     LDWORK  The length of the array DWORK. */
/*             LDWORK >= MAX(1,N*(N+2*M)+3*M,(N+P)*(N+1)+N+MAX(N-1,M+1)), */
/*                                 if JOBX = 'N'; */
/*             LDWORK >= MAX(2,N*(N+2*M)+3*M,(N+P)*(N+1)+N+MAX(N-1,M+1), */
/*                           3*N), if JOBX = 'X'. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value;                        -1 */
/*             = 1:  if JOBX = 'X' and the matrix S    is singular, */
/*                                                 i+1       -1 */
/*                   i.e., the condition number estimate of S    (in the */
/*                                                           i+1 */
/*                                                         -1    -1/2 */
/*                   1-norm) exceeds 1/TOL.  The matrices S   , Q */
/*                                                         i+1   i */
/*                   and E have been computed. */

/*     METHOD */

/*     The routine performs one recursion of the square root information */
/*     filter algorithm, summarized as follows: */

/*       |    -1/2             -1/2    |     |         -1/2             | */
/*       |   Q         0      Q    Z   |     | (QINOV )     *     *     | */
/*       |    i                i    i  |     |       i                  | */
/*       |                             |     |                          | */
/*       |           -1/2      -1/2    |     |             -1    -1     | */
/*     T |    0     R    C    R    Y   |  =  |    0       S     S   X   | */
/*       |           i+1  i+1  i+1  i+1|     |             i+1   i+1 i+1| */
/*       |                             |     |                          | */
/*       |  -1 -1     -1 -1    -1      |     |                          | */
/*       | S  A  B   S  A     S  X     |     |    0         0     E     | */
/*       |  i         i        i  i    |     |                     i+1  | */

/*                   (Pre-array)                      (Post-array) */

/*     where T is an orthogonal transformation triangularizing the */
/*                        -1/2 */
/*     pre-array, (QINOV )     is the inverse of the covariance square */
/*                      i */
/*     root (right Cholesky factor) of the process noise innovation */
/*                                                            -1  -1 */
/*     (hence the information square root) at instant i and (A  ,A  B) is */
/*     in upper controller Hessenberg form. */

/*     An example of the pre-array is given below (where N = 6, M = 2, */
/*     and P = 3): */

/*         |x x |             | x| */
/*         |  x |             | x| */
/*         _______________________ */
/*         |    | x x x x x x | x| */
/*         |    | x x x x x x | x| */
/*         |    | x x x x x x | x| */
/*         _______________________ */
/*         |x x | x x x x x x | x| */
/*         |  x | x x x x x x | x| */
/*         |    | x x x x x x | x| */
/*         |    |   x x x x x | x| */
/*         |    |     x x x x | x| */
/*         |    |       x x x | x| */

/*     The inverse of the corresponding state covariance matrix P */
/*                                                               i+1|i+1 */
/*     (hence the information matrix I) is then factorized as */

/*                    -1         -1     -1 */
/*         I       = P       = (S   )' S */
/*          i+1|i+1   i+1|i+1    i+1    i+1 */

/*     and one combined time and measurement update for the state is */
/*     given by X   . */
/*               i+1 */

/*     The triangularization is done entirely via Householder */
/*     transformations exploiting the zero pattern of the pre-array. */

/*     REFERENCES */

/*     [1] Anderson, B.D.O. and Moore, J.B. */
/*         Optimal Filtering. */
/*         Prentice Hall, Englewood Cliffs, New Jersey, 1979. */

/*     [2] Van Dooren, P. and Verhaegen, M.H.G. */
/*         Condensed Forms for Efficient Time-Invariant Kalman Filtering. */
/*         SIAM J. Sci. Stat. Comp., 9. pp. 516-530, 1988. */

/*     [3] Verhaegen, M.H.G. and Van Dooren, P. */
/*         Numerical Aspects of Different Kalman Filter Implementations. */
/*         IEEE Trans. Auto. Contr., AC-31, pp. 907-917, Oct. 1986. */

/*     [4] Vanbegin, M., Van Dooren, P., and Verhaegen, M.H.G. */
/*         Algorithm 675: FORTRAN Subroutines for Computing the Square */
/*         Root Covariance Filter and Square Root Information Filter in */
/*         Dense or Hessenberg Forms. */
/*         ACM Trans. Math. Software, 15, pp. 243-256, 1989. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires approximately */

/*           3    2                           2          3 */
/*     (1/6)N  + N x (3/2 x M + P) + 2 x N x M  + 2/3 x M */

/*     operations and is backward stable (see [3]). */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */
/*     Supersedes Release 2.0 routine FB01HD by M. Vanbegin, */
/*     P. Van Dooren, and M.H.G. Verhaegen. */

/*     REVISIONS */

/*     February 20, 1998, November 20, 2003, February 14, 2004. */

/*     KEYWORDS */

/*     Controller Hessenberg form, Kalman filtering, optimal filtering, */
/*     orthogonal transformation, recursive estimation, square-root */
/*     filtering, square-root information filtering. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 389 "FB01TD.f"
    /* Parameter adjustments */
#line 389 "FB01TD.f"
    sinv_dim1 = *ldsinv;
#line 389 "FB01TD.f"
    sinv_offset = 1 + sinv_dim1;
#line 389 "FB01TD.f"
    sinv -= sinv_offset;
#line 389 "FB01TD.f"
    ainv_dim1 = *ldainv;
#line 389 "FB01TD.f"
    ainv_offset = 1 + ainv_dim1;
#line 389 "FB01TD.f"
    ainv -= ainv_offset;
#line 389 "FB01TD.f"
    ainvb_dim1 = *ldainb;
#line 389 "FB01TD.f"
    ainvb_offset = 1 + ainvb_dim1;
#line 389 "FB01TD.f"
    ainvb -= ainvb_offset;
#line 389 "FB01TD.f"
    rinv_dim1 = *ldrinv;
#line 389 "FB01TD.f"
    rinv_offset = 1 + rinv_dim1;
#line 389 "FB01TD.f"
    rinv -= rinv_offset;
#line 389 "FB01TD.f"
    c_dim1 = *ldc;
#line 389 "FB01TD.f"
    c_offset = 1 + c_dim1;
#line 389 "FB01TD.f"
    c__ -= c_offset;
#line 389 "FB01TD.f"
    qinv_dim1 = *ldqinv;
#line 389 "FB01TD.f"
    qinv_offset = 1 + qinv_dim1;
#line 389 "FB01TD.f"
    qinv -= qinv_offset;
#line 389 "FB01TD.f"
    --x;
#line 389 "FB01TD.f"
    --rinvy;
#line 389 "FB01TD.f"
    --z__;
#line 389 "FB01TD.f"
    --e;
#line 389 "FB01TD.f"
    --iwork;
#line 389 "FB01TD.f"
    --dwork;
#line 389 "FB01TD.f"

#line 389 "FB01TD.f"
    /* Function Body */
#line 389 "FB01TD.f"
    np = *n + *p;
#line 390 "FB01TD.f"
    nm = *n + *m;
#line 391 "FB01TD.f"
    n1 = max(1,*n);
#line 392 "FB01TD.f"
    m1 = max(1,*m);
#line 393 "FB01TD.f"
    mp1 = *m + 1;
#line 394 "FB01TD.f"
    *info = 0;
#line 395 "FB01TD.f"
    ljobx = lsame_(jobx, "X", (ftnlen)1, (ftnlen)1);
#line 396 "FB01TD.f"
    lmultr = lsame_(multrc, "P", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 400 "FB01TD.f"
    if (! ljobx && ! lsame_(jobx, "N", (ftnlen)1, (ftnlen)1)) {
#line 401 "FB01TD.f"
	*info = -1;
#line 402 "FB01TD.f"
    } else if (! lmultr && ! lsame_(multrc, "N", (ftnlen)1, (ftnlen)1)) {
#line 403 "FB01TD.f"
	*info = -2;
#line 404 "FB01TD.f"
    } else if (*n < 0) {
#line 405 "FB01TD.f"
	*info = -3;
#line 406 "FB01TD.f"
    } else if (*m < 0) {
#line 407 "FB01TD.f"
	*info = -4;
#line 408 "FB01TD.f"
    } else if (*p < 0) {
#line 409 "FB01TD.f"
	*info = -5;
#line 410 "FB01TD.f"
    } else if (*ldsinv < n1) {
#line 411 "FB01TD.f"
	*info = -7;
#line 412 "FB01TD.f"
    } else if (*ldainv < n1) {
#line 413 "FB01TD.f"
	*info = -9;
#line 414 "FB01TD.f"
    } else if (*ldainb < n1) {
#line 415 "FB01TD.f"
	*info = -11;
#line 416 "FB01TD.f"
    } else if (*ldrinv < 1 || ! lmultr && *ldrinv < *p) {
#line 417 "FB01TD.f"
	*info = -13;
#line 418 "FB01TD.f"
    } else if (*ldc < max(1,*p)) {
#line 419 "FB01TD.f"
	*info = -15;
#line 420 "FB01TD.f"
    } else if (*ldqinv < m1) {
#line 421 "FB01TD.f"
	*info = -17;
#line 422 "FB01TD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
#line 422 "FB01TD.f"
	i__3 = *n - 1;
#line 422 "FB01TD.f"
	i__1 = 2, i__2 = *n * (nm + *m) + *m * 3, i__1 = max(i__1,i__2), i__2 
		= np * (*n + 1) + *n + max(i__3,mp1), i__1 = max(i__1,i__2), 
		i__2 = *n * 3;
/* Computing MAX */
/* Computing MAX */
#line 422 "FB01TD.f"
	i__6 = *n - 1;
#line 422 "FB01TD.f"
	i__4 = 1, i__5 = *n * (nm + *m) + *m * 3, i__4 = max(i__4,i__5), i__5 
		= np * (*n + 1) + *n + max(i__6,mp1);
#line 422 "FB01TD.f"
	if (ljobx && *ldwork < max(i__1,i__2) || ! ljobx && *ldwork < max(
		i__4,i__5)) {
#line 429 "FB01TD.f"
	    *info = -25;
#line 430 "FB01TD.f"
	}
#line 430 "FB01TD.f"
    }

#line 432 "FB01TD.f"
    if (*info != 0) {

/*        Error return. */

#line 436 "FB01TD.f"
	i__1 = -(*info);
#line 436 "FB01TD.f"
	xerbla_("FB01TD", &i__1, (ftnlen)6);
#line 437 "FB01TD.f"
	return 0;
#line 438 "FB01TD.f"
    }

/*     Quick return if possible. */

#line 442 "FB01TD.f"
    if (max(*n,*p) == 0) {
#line 443 "FB01TD.f"
	if (ljobx) {
#line 444 "FB01TD.f"
	    dwork[1] = 2.;
#line 445 "FB01TD.f"
	    dwork[2] = 1.;
#line 446 "FB01TD.f"
	} else {
#line 447 "FB01TD.f"
	    dwork[1] = 1.;
#line 448 "FB01TD.f"
	}
#line 449 "FB01TD.f"
	return 0;
#line 450 "FB01TD.f"
    }

/*     Construction of the needed part of the pre-array in DWORK. */
/*     To save workspace, only the blocks (1,3), (3,1)-(3,3), (2,2), and */
/*     (2,3) will be constructed when needed as shown below. */

/*     Storing SINV x AINVB and SINV x AINV in the (1,1) and (1,2) */
/*     blocks of DWORK, respectively. The upper trapezoidal structure of */
/*     [ AINVB AINV ] is fully exploited. Specifically, if M <= N, the */
/*     following partition is used: */

/*       [ S1  S2 ] [ B1  A1 A3 ] */
/*       [ 0   S3 ] [ 0   A2 A4 ], */

/*     where B1, A3, and S1 are M-by-M matrices, A1 and S2 are */
/*     M-by-(N-M), A2 and S3 are (N-M)-by-(N-M), A4 is (N-M)-by-M, and */
/*     B1, S1, A2, and S3 are upper triangular. The right hand side */
/*     matrix above is stored in the workspace. If M > N, the partition */
/*     is [ SINV ] [ B1 B2  A ], where B1 is N-by-N, B2 is N-by-(M-N), */
/*     and B1 and SINV are upper triangular. */
/*     The variables called Ixy define the starting positions where the */
/*     (x,y) blocks of the pre-array are initially stored in DWORK. */
/*     Workspace: need N*(M+N). */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 480 "FB01TD.f"
    ldw = n1;
#line 481 "FB01TD.f"
    i32 = *n * *m + 1;

#line 483 "FB01TD.f"
    dlacpy_("Upper", n, m, &ainvb[ainvb_offset], ldainb, &dwork[1], &ldw, (
	    ftnlen)5);
#line 484 "FB01TD.f"
    i__1 = min(*m,*n);
#line 484 "FB01TD.f"
    dlacpy_("Full", &i__1, n, &ainv[ainv_offset], ldainv, &dwork[i32], &ldw, (
	    ftnlen)4);
#line 486 "FB01TD.f"
    if (*n > *m) {
#line 486 "FB01TD.f"
	i__1 = *n - *m;
#line 486 "FB01TD.f"
	dlacpy_("Upper", &i__1, n, &ainv[mp1 + ainv_dim1], ldainv, &dwork[i32 
		+ *m], &ldw, (ftnlen)5);
#line 486 "FB01TD.f"
    }

/*                    [ B1  A1 ] */
/*     Compute SINV x [ 0   A2 ] or SINV x B1 as a product of upper */
/*     triangular matrices. */
/*     Workspace: need N*(M+N+1). */

#line 495 "FB01TD.f"
    ii = 1;
#line 496 "FB01TD.f"
    i13 = *n * nm + 1;
/* Computing MAX */
#line 497 "FB01TD.f"
    i__1 = 1, i__2 = *n * nm + *n;
#line 497 "FB01TD.f"
    wrkopt = max(i__1,i__2);

#line 499 "FB01TD.f"
    i__1 = *n;
#line 499 "FB01TD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 500 "FB01TD.f"
	dcopy_(&i__, &dwork[ii], &c__1, &dwork[i13], &c__1);
#line 501 "FB01TD.f"
	dtrmv_("Upper", "No transpose", "Non-unit", &i__, &sinv[sinv_offset], 
		ldsinv, &dwork[i13], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 503 "FB01TD.f"
	dcopy_(&i__, &dwork[i13], &c__1, &dwork[ii], &c__1);
#line 504 "FB01TD.f"
	ii += *n;
#line 505 "FB01TD.f"
/* L10: */
#line 505 "FB01TD.f"
    }

/*                    [ A3 ] */
/*     Compute SINV x [ A4 ] or SINV x [ B2 A ]. */

#line 510 "FB01TD.f"
    dtrmm_("Left", "Upper", "No transpose", "Non-unit", n, m, &c_b23, &sinv[
	    sinv_offset], ldsinv, &dwork[ii], &ldw, (ftnlen)4, (ftnlen)5, (
	    ftnlen)12, (ftnlen)8);

/*     Storing the process noise mean value in (1,3) block of DWORK. */
/*     Workspace: need N*(M+N) + M. */

#line 516 "FB01TD.f"
    dcopy_(m, &z__[1], &c__1, &dwork[i13], &c__1);
#line 517 "FB01TD.f"
    dtrmv_("Upper", "No transpose", "Non-unit", m, &qinv[qinv_offset], ldqinv,
	     &dwork[i13], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*     Computing SINV x X in X. */

#line 522 "FB01TD.f"
    dtrmv_("Upper", "No transpose", "Non-unit", n, &sinv[sinv_offset], ldsinv,
	     &x[1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*     Triangularization (2 steps). */

/*     Step 1: annihilate the matrix SINV x AINVB. */
/*     Workspace: need N*(N+2*M) + 3*M. */

#line 530 "FB01TD.f"
    i12 = i13 + *m;
#line 531 "FB01TD.f"
    itau = i12 + *m * *n;
#line 532 "FB01TD.f"
    jwork = itau + *m;

#line 534 "FB01TD.f"
    mb04kd_("Upper", m, n, n, &qinv[qinv_offset], ldqinv, &dwork[1], &ldw, &
	    dwork[i32], &ldw, &dwork[i12], &m1, &dwork[itau], &dwork[jwork], (
	    ftnlen)5);
/* Computing MAX */
#line 537 "FB01TD.f"
    i__1 = wrkopt, i__2 = *n * (nm + *m) + *m * 3;
#line 537 "FB01TD.f"
    wrkopt = max(i__1,i__2);

#line 539 "FB01TD.f"
    if (*n == 0) {
#line 540 "FB01TD.f"
	dcopy_(p, &rinvy[1], &c__1, &e[1], &c__1);
#line 541 "FB01TD.f"
	if (ljobx) {
#line 541 "FB01TD.f"
	    dwork[2] = 1.;
#line 541 "FB01TD.f"
	}
#line 543 "FB01TD.f"
	dwork[1] = (doublereal) wrkopt;
#line 544 "FB01TD.f"
	return 0;
#line 545 "FB01TD.f"
    }

/*     Apply the transformations to the last column of the pre-array. */
/*     (Only the updated (3,3) block is now needed.) */

#line 550 "FB01TD.f"
    ij = 1;

#line 552 "FB01TD.f"
    i__1 = *m;
#line 552 "FB01TD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 553 "FB01TD.f"
	i__2 = min(i__,*n);
#line 553 "FB01TD.f"
	i__3 = min(i__,*n);
#line 553 "FB01TD.f"
	d__1 = -dwork[itau + i__ - 1] * (dwork[i13 + i__ - 1] + ddot_(&i__3, &
		dwork[ij], &c__1, &x[1], &c__1));
#line 553 "FB01TD.f"
	daxpy_(&i__2, &d__1, &dwork[ij], &c__1, &x[1], &c__1);
#line 556 "FB01TD.f"
	ij += *n;
#line 557 "FB01TD.f"
/* L20: */
#line 557 "FB01TD.f"
    }

/*     Now, the workspace for SINV x AINVB, as well as for the updated */
/*     (1,2) block of the pre-array, are no longer needed. */
/*     Move the computed (3,2) and (3,3) blocks of the pre-array in the */
/*     (1,1) and (1,2) block positions of DWORK, to save space for the */
/*     following computations. */
/*     Then, adjust the implicitly defined leading dimension of DWORK, */
/*     to make space for storing the (2,2) and (2,3) blocks of the */
/*     pre-array. */
/*     Workspace: need (P+N)*(N+1). */

#line 569 "FB01TD.f"
    i__1 = min(*m,*n);
#line 569 "FB01TD.f"
    dlacpy_("Full", &i__1, n, &dwork[i32], &ldw, &dwork[1], &ldw, (ftnlen)4);
#line 570 "FB01TD.f"
    if (*n > *m) {
#line 570 "FB01TD.f"
	i__1 = *n - *m;
#line 570 "FB01TD.f"
	dlacpy_("Upper", &i__1, n, &dwork[i32 + *m], &ldw, &dwork[mp1], &ldw, 
		(ftnlen)5);
#line 570 "FB01TD.f"
    }
#line 573 "FB01TD.f"
    ldw = max(1,np);

#line 575 "FB01TD.f"
    for (i__ = *n; i__ >= 1; --i__) {
/* Computing MIN */
#line 576 "FB01TD.f"
	i__1 = *n, i__2 = i__ + *m;
#line 576 "FB01TD.f"
	for (ij = min(i__1,i__2); ij >= 1; --ij) {
#line 577 "FB01TD.f"
	    dwork[np * (i__ - 1) + *p + ij] = dwork[*n * (i__ - 1) + ij];
#line 578 "FB01TD.f"
/* L30: */
#line 578 "FB01TD.f"
	}
#line 579 "FB01TD.f"
/* L40: */
#line 579 "FB01TD.f"
    }

/*     Copy of RINV x C in the (1,1) block of DWORK. */

#line 583 "FB01TD.f"
    dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[1], &ldw, (ftnlen)4);
#line 584 "FB01TD.f"
    if (! lmultr) {
#line 584 "FB01TD.f"
	dtrmm_("Left", "Upper", "No transpose", "Non-unit", p, n, &c_b23, &
		rinv[rinv_offset], ldrinv, &dwork[1], &ldw, (ftnlen)4, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 584 "FB01TD.f"
    }

/*     Copy the inclusion measurement in the (1,2) block and the updated */
/*     X in the (2,2) block of DWORK. */

#line 591 "FB01TD.f"
    i23 = np * *n + 1;
#line 592 "FB01TD.f"
    i33 = i23 + *p;
#line 593 "FB01TD.f"
    dcopy_(p, &rinvy[1], &c__1, &dwork[i23], &c__1);
#line 594 "FB01TD.f"
    dcopy_(n, &x[1], &c__1, &dwork[i33], &c__1);
/* Computing MAX */
#line 595 "FB01TD.f"
    i__1 = wrkopt, i__2 = np * (*n + 1);
#line 595 "FB01TD.f"
    wrkopt = max(i__1,i__2);

/*     Step 2: QR factorization of the first block column of the matrix */

/*        [ RINV x C     RINV x Y ], */
/*        [ SINV x AINV  SINV x X ] */

/*     where the second block row was modified at Step 1. */
/*     Workspace: need   (P+N)*(N+1) + N + MAX(N-1,M+1); */
/*                prefer (P+N)*(N+1) + N + (M+1)*NB, where NB is the */
/*                       optimal block size for DGEQRF called in MB04ID. */

#line 607 "FB01TD.f"
    itau = i23 + np;
#line 608 "FB01TD.f"
    jwork = itau + *n;

/* Computing MAX */
#line 610 "FB01TD.f"
    i__2 = *n - mp1;
#line 610 "FB01TD.f"
    i__1 = max(i__2,0);
#line 610 "FB01TD.f"
    i__3 = *ldwork - jwork + 1;
#line 610 "FB01TD.f"
    mb04id_(&np, n, &i__1, &c__1, &dwork[1], &ldw, &dwork[i23], &ldw, &dwork[
	    itau], &dwork[jwork], &i__3, info);
/* Computing MAX */
#line 613 "FB01TD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 613 "FB01TD.f"
    wrkopt = max(i__1,i__2);

/*     Output SINV, X, and E and set the optimal workspace dimension */
/*     (and the reciprocal of the condition number estimate). */

#line 618 "FB01TD.f"
    dlacpy_("Upper", n, n, &dwork[1], &ldw, &sinv[sinv_offset], ldsinv, (
	    ftnlen)5);
#line 619 "FB01TD.f"
    dcopy_(n, &dwork[i23], &c__1, &x[1], &c__1);
#line 620 "FB01TD.f"
    if (*p > 0) {
#line 620 "FB01TD.f"
	dcopy_(p, &dwork[i23 + *n], &c__1, &e[1], &c__1);
#line 620 "FB01TD.f"
    }

#line 623 "FB01TD.f"
    if (ljobx) {

/*        Compute X. */
/*        Workspace: need 3*N. */

#line 628 "FB01TD.f"
	mb02od_("Left", "Upper", "No transpose", "Non-unit", "1-norm", n, &
		c__1, &c_b23, &sinv[sinv_offset], ldsinv, &x[1], n, &rcond, 
		tol, &iwork[1], &dwork[1], info, (ftnlen)4, (ftnlen)5, (
		ftnlen)12, (ftnlen)8, (ftnlen)6);
#line 631 "FB01TD.f"
	if (*info == 0) {
/* Computing MAX */
#line 632 "FB01TD.f"
	    i__1 = wrkopt, i__2 = *n * 3;
#line 632 "FB01TD.f"
	    wrkopt = max(i__1,i__2);
#line 633 "FB01TD.f"
	    dwork[2] = rcond;
#line 634 "FB01TD.f"
	}
#line 635 "FB01TD.f"
    }

#line 637 "FB01TD.f"
    dwork[1] = (doublereal) wrkopt;

#line 639 "FB01TD.f"
    return 0;
/* *** Last line of FB01TD*** */
} /* fb01td_ */

