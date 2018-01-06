#line 1 "MB04YD.f"
/* MB04YD.f -- translated by f2c (version 20100827).
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

#line 1 "MB04YD.f"
/* Table of constant values */

static doublereal c_b13 = -.125;
static integer c__1 = 1;
static doublereal c_b18 = 0.;
static doublereal c_b19 = 1.;

/* Subroutine */ int mb04yd_(char *jobu, char *jobv, integer *m, integer *n, 
	integer *rank, doublereal *theta, doublereal *q, doublereal *e, 
	doublereal *u, integer *ldu, doublereal *v, integer *ldv, logical *
	inul, doublereal *tol, doublereal *reltol, doublereal *dwork, integer 
	*ldwork, integer *iwarn, integer *info, ftnlen jobu_len, ftnlen 
	jobv_len)
{
    /* System generated locals */
    integer u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, p, r__;
    static doublereal x;
    static integer i1;
    static doublereal eps;
    static logical noc12;
    static integer oldi, oldk;
    static doublereal cosl;
    static integer iter;
    static doublereal rmin, cosr, rmax, sinl, smax, sinr;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static logical qrit;
    static integer info1;
    extern /* Subroutine */ int mb03md_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *);
    extern integer mb03nd_(integer *, doublereal *, doublereal *, doublereal *
	    , doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer iascl;
    static doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb02ny_(logical *, logical *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *);
    static doublereal shift, sigmn;
    static integer maxit;
    extern /* Subroutine */ int mb04yw_(logical *, logical *, logical *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *);
    static doublereal sigmx;
    extern /* Subroutine */ int dlasv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal thetac, safemn;
    static logical ljobua, ljobva;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical ljobui, ljobvi;
    static integer numeig;
    static doublereal tolabs, thresh, pivmin, tolrel, smlnum;


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

/*     To partially diagonalize the bidiagonal matrix */

/*               |q(1) e(1)  0    ...       0      | */
/*               | 0   q(2) e(2)            .      | */
/*           J = | .                        .      |                  (1) */
/*               | .                  e(MIN(M,N)-1)| */
/*               | 0   ...        ...  q(MIN(M,N)) | */

/*     using QR or QL iterations in such a way that J is split into */
/*     unreduced bidiagonal submatrices whose singular values are either */
/*     all larger than a given bound or are all smaller than (or equal */
/*     to) this bound. The left- and right-hand Givens rotations */
/*     performed on J (corresponding to each QR or QL iteration step) may */
/*     be optionally accumulated in the arrays U and V. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBU    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix U the left-hand Givens rotations, as follows: */
/*             = 'N':  Do not form U; */
/*             = 'I':  U is initialized to the M-by-MIN(M,N) submatrix of */
/*                     the unit matrix and the left-hand Givens rotations */
/*                     are accumulated in U; */
/*             = 'U':  The given matrix U is updated by the left-hand */
/*                     Givens rotations used in the calculation. */

/*     JOBV    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix V the right-hand Givens rotations, as follows: */
/*             = 'N':  Do not form V; */
/*             = 'I':  V is initialized to the N-by-MIN(M,N) submatrix of */
/*                     the unit matrix and the right-hand Givens */
/*                     rotations are accumulated in V; */
/*             = 'U':  The given matrix V is updated by the right-hand */
/*                     Givens rotations used in the calculation. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows in matrix U.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of rows in matrix V.  N >= 0. */

/*     RANK    (input/output) INTEGER */
/*             On entry, if RANK < 0, then the rank of matrix J is */
/*             computed by the routine as the number of singular values */
/*             larger than THETA. */
/*             Otherwise, RANK must specify the rank of matrix J. */
/*             RANK <= MIN(M,N). */
/*             On exit, if RANK < 0 on entry, then RANK contains the */
/*             computed rank of J. That is, the number of singular */
/*             values of J larger than THETA. */
/*             Otherwise, the user-supplied value of RANK may be */
/*             changed by the routine on exit if the RANK-th and the */
/*             (RANK+1)-th singular values of J are considered to be */
/*             equal. See also the parameter TOL. */

/*     THETA   (input/output) DOUBLE PRECISION */
/*             On entry, if RANK < 0, then THETA must specify an upper */
/*             bound on the smallest singular values of J. THETA >= 0.0. */
/*             Otherwise, THETA must specify an initial estimate (t say) */
/*             for computing an upper bound such that precisely RANK */
/*             singular values are greater than this bound. */
/*             If THETA < 0.0, then t is computed by the routine. */
/*             On exit, if RANK >= 0 on entry, then THETA contains the */
/*             computed upper bound such that precisely RANK singular */
/*             values of J are greater than THETA + TOL. */
/*             Otherwise, THETA is unchanged. */

/*     Q       (input/output) DOUBLE PRECISION array, dimension */
/*             (MIN(M,N)) */
/*             On entry, this array must contain the diagonal elements */
/*             q(1),q(2),...,q(MIN(M,N)) of the bidiagonal matrix J. That */
/*             is, Q(i) = J(i,i) for i = 1,2,...,MIN(M,N). */
/*             On exit, this array contains the leading diagonal of the */
/*             transformed bidiagonal matrix J. */

/*     E       (input/output) DOUBLE PRECISION array, dimension */
/*             (MIN(M,N)-1) */
/*             On entry, this array must contain the superdiagonal */
/*             elements e(1),e(2),...,e(MIN(M,N)-1) of the bidiagonal */
/*             matrix J. That is, E(k) = J(k,k+1) for k = 1,2,..., */
/*             MIN(M,N)-1. */
/*             On exit, this array contains the superdiagonal of the */
/*             transformed bidiagonal matrix J. */

/*     U       (input/output) DOUBLE PRECISION array, dimension (LDU,*) */
/*             On entry, if JOBU = 'U', the leading M-by-MIN(M,N) part */
/*             of this array must contain a left transformation matrix */
/*             applied to the original matrix of the problem, and */
/*             on exit, the leading M-by-MIN(M,N) part of this array */
/*             contains the product of the input matrix U and the */
/*             left-hand Givens rotations. */
/*             On exit, if JOBU = 'I', then the leading M-by-MIN(M,N) */
/*             part of this array contains the matrix of accumulated */
/*             left-hand Givens rotations used. */
/*             If JOBU = 'N', the array U is not referenced and can be */
/*             supplied as a dummy array (i.e. set parameter LDU = 1 and */
/*             declare this array to be U(1,1) in the calling program). */

/*     LDU     INTEGER */
/*             The leading dimension of array U. If JOBU = 'U' or */
/*             JOBU = 'I', LDU >= MAX(1,M); if JOBU = 'N', LDU >= 1. */

/*     V       (input/output) DOUBLE PRECISION array, dimension (LDV,*) */
/*             On entry, if JOBV = 'U', the leading N-by-MIN(M,N) part */
/*             of this array must contain a right transformation matrix */
/*             applied to the original matrix of the problem, and */
/*             on exit, the leading N-by-MIN(M,N) part of this array */
/*             contains the product of the input matrix V and the */
/*             right-hand Givens rotations. */
/*             On exit, if JOBV = 'I', then the leading N-by-MIN(M,N) */
/*             part of this array contains the matrix of accumulated */
/*             right-hand Givens rotations used. */
/*             If JOBV = 'N', the array V is not referenced and can be */
/*             supplied as a dummy array (i.e. set parameter LDV = 1 and */
/*             declare this array to be V(1,1) in the calling program). */

/*     LDV     INTEGER */
/*             The leading dimension of array V. If JOBV = 'U' or */
/*             JOBV = 'I', LDV >= MAX(1,N); if JOBV = 'N', LDV >= 1. */

/*     INUL    (input/output) LOGICAL array, dimension (MIN(M,N)) */
/*             On entry, the leading MIN(M,N) elements of this array must */
/*             be set to .FALSE. unless the i-th columns of U (if JOBU = */
/*             'U') and V (if JOBV = 'U') already contain a computed base */
/*             vector of the desired singular subspace of the original */
/*             matrix, in which case INUL(i) must be set to .TRUE. */
/*             for 1 <= i <= MIN(M,N). */
/*             On exit, the indices of the elements of this array with */
/*             value .TRUE. indicate the indices of the diagonal entries */
/*             of J which belong to those bidiagonal submatrices whose */
/*             singular values are all less than or equal to THETA. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             This parameter defines the multiplicity of singular values */
/*             by considering all singular values within an interval of */
/*             length TOL as coinciding. TOL is used in checking how many */
/*             singular values are less than or equal to THETA. Also in */
/*             computing an appropriate upper bound THETA by a bisection */
/*             method, TOL is used as a stopping criterion defining the */
/*             minimum (absolute) subinterval width. TOL is also taken */
/*             as an absolute tolerance for negligible elements in the */
/*             QR/QL iterations. If the user sets TOL to be less than or */
/*             equal to 0, then the tolerance is taken as */
/*             EPS * MAX(ABS(Q(i)), ABS(E(k))), where EPS is the */
/*             machine precision (see LAPACK Library routine DLAMCH), */
/*             i = 1,2,...,MIN(M,N) and k = 1,2,...,MIN(M,N)-1. */

/*     RELTOL  DOUBLE PRECISION */
/*             This parameter specifies the minimum relative width of an */
/*             interval. When an interval is narrower than TOL, or than */
/*             RELTOL times the larger (in magnitude) endpoint, then it */
/*             is considered to be sufficiently small and bisection has */
/*             converged. If the user sets RELTOL to be less than */
/*             BASE * EPS, where BASE is machine radix and EPS is machine */
/*             precision (see LAPACK Library routine DLAMCH), then the */
/*             tolerance is taken as BASE * EPS. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,6*MIN(M,N)-5), if JOBU = 'I' or 'U', or */
/*                                               JOBV = 'I' or 'U'; */
/*             LDWORK >= MAX(1,4*MIN(M,N)-3), if JOBU = 'N' and */
/*                                               JOBV = 'N'. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  if the rank of the bidiagonal matrix J (as specified */
/*                   by the user) has been lowered because a singular */
/*                   value of multiplicity larger than 1 was found. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; this includes values like RANK > MIN(M,N), or */
/*                   THETA < 0.0 and RANK < 0; */
/*             = 1:  if the maximum number of QR/QL iteration steps */
/*                   (30*MIN(M,N)) has been exceeded. */

/*     METHOD */

/*     If the upper bound THETA is not specified by the user, then it is */
/*     computed by the routine (using a bisection method) such that */
/*     precisely (MIN(M,N) - RANK) singular values of J are less than or */
/*     equal to THETA + TOL. */

/*     The method used by the routine (see [1]) then proceeds as follows. */

/*     The unreduced bidiagonal submatrices of J(j), where J(j) is the */
/*     transformed bidiagonal matrix after the j-th iteration step, are */
/*     classified into the following three classes: */

/*     - C1 contains the bidiagonal submatrices with all singular values */
/*       > THETA, */
/*     - C2 contains the bidiagonal submatrices with all singular values */
/*       <= THETA and */
/*     - C3 contains the bidiagonal submatrices with singular values */
/*       > THETA and also singular values <= THETA. */

/*     If C3 is empty, then the partial diagonalization is complete, and */
/*     RANK is the sum of the dimensions of the bidiagonal submatrices of */
/*     C1. */
/*     Otherwise, QR or QL iterations are performed on each bidiagonal */
/*     submatrix of C3, until this bidiagonal submatrix has been split */
/*     into two bidiagonal submatrices. These two submatrices are then */
/*     classified and the iterations are restarted. */
/*     If the upper left diagonal element of the bidiagonal submatrix is */
/*     larger than its lower right diagonal element, then QR iterations */
/*     are performed, else QL iterations are used. The shift is taken as */
/*     the smallest diagonal element of the bidiagonal submatrix (in */
/*     magnitude) unless its value exceeds THETA, in which case it is */
/*     taken as zero. */

/*     REFERENCES */

/*     [1] Van Huffel, S., Vandewalle, J. and Haegemans, A. */
/*         An efficient and reliable algorithm for computing the */
/*         singular subspace of a matrix associated with its smallest */
/*         singular values. */
/*         J. Comput. and Appl. Math., 19, pp. 313-330, 1987. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     To avoid overflow, matrix J is scaled so that its largest element */
/*     is no greater than  overflow**(1/2) * underflow**(1/4) in absolute */
/*     value (and not much smaller than that, for maximal accuracy). */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1997. */
/*     Supersedes Release 2.0 routine MB04QD by S. Van Huffel, Katholieke */
/*     University Leuven, Belgium. */

/*     REVISIONS */

/*     July 10, 1997. V. Sima. */
/*     November 25, 1997. V. Sima: Setting INUL(K) = .TRUE. when handling */
/*                                 2-by-2 submatrix. */

/*     KEYWORDS */

/*     Bidiagonal matrix, orthogonal transformation, singular values. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 321 "MB04YD.f"
    /* Parameter adjustments */
#line 321 "MB04YD.f"
    --q;
#line 321 "MB04YD.f"
    --e;
#line 321 "MB04YD.f"
    u_dim1 = *ldu;
#line 321 "MB04YD.f"
    u_offset = 1 + u_dim1;
#line 321 "MB04YD.f"
    u -= u_offset;
#line 321 "MB04YD.f"
    v_dim1 = *ldv;
#line 321 "MB04YD.f"
    v_offset = 1 + v_dim1;
#line 321 "MB04YD.f"
    v -= v_offset;
#line 321 "MB04YD.f"
    --inul;
#line 321 "MB04YD.f"
    --dwork;
#line 321 "MB04YD.f"

#line 321 "MB04YD.f"
    /* Function Body */
#line 321 "MB04YD.f"
    p = min(*m,*n);
#line 322 "MB04YD.f"
    *info = 0;
#line 323 "MB04YD.f"
    *iwarn = 0;
#line 324 "MB04YD.f"
    ljobui = lsame_(jobu, "I", (ftnlen)1, (ftnlen)1);
#line 325 "MB04YD.f"
    ljobvi = lsame_(jobv, "I", (ftnlen)1, (ftnlen)1);
#line 326 "MB04YD.f"
    ljobua = ljobui || lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);
#line 327 "MB04YD.f"
    ljobva = ljobvi || lsame_(jobv, "U", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 331 "MB04YD.f"
    if (! ljobua && ! lsame_(jobu, "N", (ftnlen)1, (ftnlen)1)) {
#line 332 "MB04YD.f"
	*info = -1;
#line 333 "MB04YD.f"
    } else if (! ljobva && ! lsame_(jobv, "N", (ftnlen)1, (ftnlen)1)) {
#line 334 "MB04YD.f"
	*info = -2;
#line 335 "MB04YD.f"
    } else if (*m < 0) {
#line 336 "MB04YD.f"
	*info = -3;
#line 337 "MB04YD.f"
    } else if (*n < 0) {
#line 338 "MB04YD.f"
	*info = -4;
#line 339 "MB04YD.f"
    } else if (*rank > p) {
#line 340 "MB04YD.f"
	*info = -5;
#line 341 "MB04YD.f"
    } else if (*rank < 0 && *theta < 0.) {
#line 342 "MB04YD.f"
	*info = -6;
#line 343 "MB04YD.f"
    } else if (! ljobua && *ldu < 1 || ljobua && *ldu < max(1,*m)) {
#line 345 "MB04YD.f"
	*info = -10;
#line 346 "MB04YD.f"
    } else if (! ljobva && *ldv < 1 || ljobva && *ldv < max(1,*n)) {
#line 348 "MB04YD.f"
	*info = -12;
#line 349 "MB04YD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 349 "MB04YD.f"
	i__1 = 1, i__2 = p * 6 - 5;
/* Computing MAX */
#line 349 "MB04YD.f"
	i__3 = 1, i__4 = (p << 2) - 3;
#line 349 "MB04YD.f"
	if ((ljobua || ljobva) && *ldwork < max(i__1,i__2) || ! (ljobua || 
		ljobva) && *ldwork < max(i__3,i__4)) {
#line 352 "MB04YD.f"
	    *info = -17;
#line 353 "MB04YD.f"
	}
#line 353 "MB04YD.f"
    }

#line 355 "MB04YD.f"
    if (*info != 0) {

/*        Error return. */

#line 359 "MB04YD.f"
	i__1 = -(*info);
#line 359 "MB04YD.f"
	xerbla_("MB04YD", &i__1, (ftnlen)6);
#line 360 "MB04YD.f"
	return 0;
#line 361 "MB04YD.f"
    }

/*     Quick return if possible. */

#line 365 "MB04YD.f"
    if (p == 0) {
#line 366 "MB04YD.f"
	if (*rank >= 0) {
#line 366 "MB04YD.f"
	    *theta = 0.;
#line 366 "MB04YD.f"
	}
#line 368 "MB04YD.f"
	*rank = 0;
#line 369 "MB04YD.f"
	return 0;
#line 370 "MB04YD.f"
    }

/*     Set tolerances and machine parameters. */

#line 374 "MB04YD.f"
    tolabs = *tol;
#line 375 "MB04YD.f"
    tolrel = *reltol;
#line 376 "MB04YD.f"
    smax = (d__1 = q[p], abs(d__1));

#line 378 "MB04YD.f"
    i__1 = p - 1;
#line 378 "MB04YD.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 379 "MB04YD.f"
	d__3 = smax, d__4 = (d__1 = q[j], abs(d__1)), d__3 = max(d__3,d__4), 
		d__4 = (d__2 = e[j], abs(d__2));
#line 379 "MB04YD.f"
	smax = max(d__3,d__4);
#line 380 "MB04YD.f"
/* L20: */
#line 380 "MB04YD.f"
    }

#line 382 "MB04YD.f"
    safemn = dlamch_("Safe minimum", (ftnlen)12);
#line 383 "MB04YD.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 384 "MB04YD.f"
    if (tolabs <= 0.) {
#line 384 "MB04YD.f"
	tolabs = eps * smax;
#line 384 "MB04YD.f"
    }
#line 385 "MB04YD.f"
    x = dlamch_("Base", (ftnlen)4) * eps;
#line 386 "MB04YD.f"
    if (tolrel <= x) {
#line 386 "MB04YD.f"
	tolrel = x;
#line 386 "MB04YD.f"
    }
/* Computing MAX */
/* Computing MIN */
#line 387 "MB04YD.f"
    d__3 = 100., d__4 = pow_dd(&eps, &c_b13);
#line 387 "MB04YD.f"
    d__1 = 10., d__2 = min(d__3,d__4);
#line 387 "MB04YD.f"
    thresh = max(d__1,d__2) * eps;
#line 388 "MB04YD.f"
    smlnum = safemn / eps;
#line 389 "MB04YD.f"
    rmin = sqrt(smlnum);
/* Computing MIN */
#line 390 "MB04YD.f"
    d__1 = 1. / rmin, d__2 = 1. / sqrt(sqrt(safemn));
#line 390 "MB04YD.f"
    rmax = min(d__1,d__2);
#line 391 "MB04YD.f"
    thetac = *theta;

/*     Scale the matrix to allowable range, if necessary, and set PIVMIN, */
/*     using the squares of Q and E (saved in DWORK). */

#line 396 "MB04YD.f"
    iascl = 0;
#line 397 "MB04YD.f"
    if (smax > 0. && smax < rmin) {
#line 398 "MB04YD.f"
	iascl = 1;
#line 399 "MB04YD.f"
	sigma = rmin / smax;
#line 400 "MB04YD.f"
    } else if (smax > rmax) {
#line 401 "MB04YD.f"
	iascl = 1;
#line 402 "MB04YD.f"
	sigma = rmax / smax;
#line 403 "MB04YD.f"
    }
#line 404 "MB04YD.f"
    if (iascl == 1) {
#line 405 "MB04YD.f"
	dscal_(&p, &sigma, &q[1], &c__1);
#line 406 "MB04YD.f"
	i__1 = p - 1;
#line 406 "MB04YD.f"
	dscal_(&i__1, &sigma, &e[1], &c__1);
#line 407 "MB04YD.f"
	thetac = sigma * *theta;
#line 408 "MB04YD.f"
	tolabs = sigma * tolabs;
#line 409 "MB04YD.f"
    }

/* Computing 2nd power */
#line 411 "MB04YD.f"
    d__1 = q[p];
#line 411 "MB04YD.f"
    pivmin = d__1 * d__1;
#line 412 "MB04YD.f"
    dwork[p] = pivmin;

#line 414 "MB04YD.f"
    i__1 = p - 1;
#line 414 "MB04YD.f"
    for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
#line 415 "MB04YD.f"
	d__1 = q[j];
#line 415 "MB04YD.f"
	dwork[j] = d__1 * d__1;
/* Computing 2nd power */
#line 416 "MB04YD.f"
	d__1 = e[j];
#line 416 "MB04YD.f"
	dwork[p + j] = d__1 * d__1;
/* Computing MAX */
#line 417 "MB04YD.f"
	d__1 = pivmin, d__2 = dwork[j], d__1 = max(d__1,d__2), d__2 = dwork[p 
		+ j];
#line 417 "MB04YD.f"
	pivmin = max(d__1,d__2);
#line 418 "MB04YD.f"
/* L40: */
#line 418 "MB04YD.f"
    }

/* Computing MAX */
#line 420 "MB04YD.f"
    d__1 = pivmin * safemn;
#line 420 "MB04YD.f"
    pivmin = max(d__1,safemn);

/*     Initialize U and/or V to the identity matrix, if needed. */

#line 424 "MB04YD.f"
    if (ljobui) {
#line 424 "MB04YD.f"
	dlaset_("Full", m, &p, &c_b18, &c_b19, &u[u_offset], ldu, (ftnlen)4);
#line 424 "MB04YD.f"
    }
#line 426 "MB04YD.f"
    if (ljobvi) {
#line 426 "MB04YD.f"
	dlaset_("Full", n, &p, &c_b18, &c_b19, &v[v_offset], ldv, (ftnlen)4);
#line 426 "MB04YD.f"
    }

/*     Estimate THETA (if not fixed by the user), and set R. */

#line 431 "MB04YD.f"
    if (*rank >= 0) {
#line 432 "MB04YD.f"
	j = p - *rank;
#line 433 "MB04YD.f"
	mb03md_(&p, &j, &thetac, &q[1], &e[1], &dwork[1], &dwork[p + 1], &
		pivmin, &tolabs, &tolrel, iwarn, &info1);
#line 435 "MB04YD.f"
	*theta = thetac;
#line 436 "MB04YD.f"
	if (iascl == 1) {
#line 436 "MB04YD.f"
	    *theta /= sigma;
#line 436 "MB04YD.f"
	}
#line 437 "MB04YD.f"
	if (j <= 0) {
#line 437 "MB04YD.f"
	    return 0;
#line 437 "MB04YD.f"
	}
#line 439 "MB04YD.f"
	r__ = p - j;
#line 440 "MB04YD.f"
    } else {
#line 441 "MB04YD.f"
	r__ = p - mb03nd_(&p, &thetac, &dwork[1], &dwork[p + 1], &pivmin, &
		info1);
#line 442 "MB04YD.f"
    }

#line 444 "MB04YD.f"
    *rank = p;

#line 446 "MB04YD.f"
    i__1 = p;
#line 446 "MB04YD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 447 "MB04YD.f"
	if (inul[i__]) {
#line 447 "MB04YD.f"
	    --(*rank);
#line 447 "MB04YD.f"
	}
#line 448 "MB04YD.f"
/* L60: */
#line 448 "MB04YD.f"
    }

/*     From now on K is the smallest known index such that the elements */
/*     of the bidiagonal matrix J with indices larger than K belong to C1 */
/*     or C2. */
/*     RANK = P - SUM(dimensions of known bidiagonal matrices of C2). */

#line 455 "MB04YD.f"
    k = p;
#line 456 "MB04YD.f"
    oldi = -1;
#line 457 "MB04YD.f"
    oldk = -1;
#line 458 "MB04YD.f"
    iter = 0;
#line 459 "MB04YD.f"
    maxit = p * 30;
/*     WHILE ( C3 NOT EMPTY ) DO */
#line 461 "MB04YD.f"
L80:
#line 461 "MB04YD.f"
    if (*rank > r__ && k > 0) {
/*        WHILE ( K.GT.0 .AND. INUL(K) ) DO */

/*        Search for the rightmost index of a bidiagonal submatrix, */
/*        not yet classified. */

#line 467 "MB04YD.f"
L100:
#line 467 "MB04YD.f"
	if (k > 0) {
#line 468 "MB04YD.f"
	    if (inul[k]) {
#line 469 "MB04YD.f"
		--k;
#line 470 "MB04YD.f"
		goto L100;
#line 471 "MB04YD.f"
	    }
#line 472 "MB04YD.f"
	}
/*        END WHILE 100 */

#line 475 "MB04YD.f"
	if (k == 0) {
#line 475 "MB04YD.f"
	    return 0;
#line 475 "MB04YD.f"
	}

#line 478 "MB04YD.f"
	noc12 = TRUE_;
/*        WHILE ((ITER < MAXIT).AND.(No bidiagonal matrix of C1 or */
/*                C2 found)) DO */
#line 481 "MB04YD.f"
L120:
#line 481 "MB04YD.f"
	if (iter < maxit && noc12) {

/*           Search for negligible Q(I) or E(I-1) (for I > 1) and find */
/*           the shift. */

#line 486 "MB04YD.f"
	    i__ = k;
#line 487 "MB04YD.f"
	    x = (d__1 = q[i__], abs(d__1));
#line 488 "MB04YD.f"
	    shift = x;
/*           WHILE ABS( Q(I) ) > TOLABS .AND. ABS( E(I-1) ) > TOLABS ) DO */
#line 490 "MB04YD.f"
L140:
#line 490 "MB04YD.f"
	    if (i__ > 1) {
#line 491 "MB04YD.f"
		if (x > tolabs && (d__1 = e[i__ - 1], abs(d__1)) > tolabs) {
#line 493 "MB04YD.f"
		    --i__;
#line 494 "MB04YD.f"
		    x = (d__1 = q[i__], abs(d__1));
#line 495 "MB04YD.f"
		    if (x < shift) {
#line 495 "MB04YD.f"
			shift = x;
#line 495 "MB04YD.f"
		    }
#line 496 "MB04YD.f"
		    goto L140;
#line 497 "MB04YD.f"
		}
#line 498 "MB04YD.f"
	    }
/*           END WHILE 140 */

/*           Classify the bidiagonal submatrix (of order J) found. */

#line 503 "MB04YD.f"
	    j = k - i__ + 1;
#line 504 "MB04YD.f"
	    if (x <= tolabs || k == i__) {
#line 505 "MB04YD.f"
		noc12 = FALSE_;
#line 506 "MB04YD.f"
	    } else {
#line 507 "MB04YD.f"
		numeig = mb03nd_(&j, &thetac, &dwork[i__], &dwork[p + i__], &
			pivmin, &info1);
#line 509 "MB04YD.f"
		if (numeig >= j || numeig <= 0) {
#line 509 "MB04YD.f"
		    noc12 = FALSE_;
#line 509 "MB04YD.f"
		}
#line 510 "MB04YD.f"
	    }
#line 511 "MB04YD.f"
	    if (noc12) {
#line 512 "MB04YD.f"
		if (j == 2) {

/*                 Handle separately the 2-by-2 submatrix. */

#line 516 "MB04YD.f"
		    dlasv2_(&q[i__], &e[i__], &q[k], &sigmn, &sigmx, &sinr, &
			    cosr, &sinl, &cosl);
#line 518 "MB04YD.f"
		    q[i__] = sigmx;
#line 519 "MB04YD.f"
		    q[k] = sigmn;
#line 520 "MB04YD.f"
		    e[i__] = 0.;
#line 521 "MB04YD.f"
		    --(*rank);
#line 522 "MB04YD.f"
		    inul[k] = TRUE_;
#line 523 "MB04YD.f"
		    noc12 = FALSE_;

/*                 Update U and/or V, if needed. */

#line 527 "MB04YD.f"
		    if (ljobua) {
#line 527 "MB04YD.f"
			drot_(m, &u[i__ * u_dim1 + 1], &c__1, &u[k * u_dim1 + 
				1], &c__1, &cosl, &sinl);
#line 527 "MB04YD.f"
		    }
#line 529 "MB04YD.f"
		    if (ljobva) {
#line 529 "MB04YD.f"
			drot_(n, &v[i__ * v_dim1 + 1], &c__1, &v[k * v_dim1 + 
				1], &c__1, &cosr, &sinr);
#line 529 "MB04YD.f"
		    }
#line 531 "MB04YD.f"
		} else {

/*                 If working on new submatrix, choose QR or */
/*                 QL iteration. */

#line 536 "MB04YD.f"
		    if (i__ != oldi || k != oldk) {
#line 536 "MB04YD.f"
			qrit = (d__1 = q[i__], abs(d__1)) >= (d__2 = q[k], 
				abs(d__2));
#line 536 "MB04YD.f"
		    }
#line 538 "MB04YD.f"
		    oldi = i__;
#line 539 "MB04YD.f"
		    if (qrit) {
#line 540 "MB04YD.f"
			if ((d__2 = e[k - 1], abs(d__2)) <= thresh * (d__1 = 
				q[k], abs(d__1))) {
#line 540 "MB04YD.f"
			    e[k - 1] = 0.;
#line 540 "MB04YD.f"
			}
#line 542 "MB04YD.f"
		    } else {
#line 543 "MB04YD.f"
			if ((d__2 = e[i__], abs(d__2)) <= thresh * (d__1 = q[
				i__], abs(d__1))) {
#line 543 "MB04YD.f"
			    e[i__] = 0.;
#line 543 "MB04YD.f"
			}
#line 545 "MB04YD.f"
		    }

#line 547 "MB04YD.f"
		    mb04yw_(&qrit, &ljobua, &ljobva, m, n, &i__, &k, &shift, &
			    q[1], &e[1], &u[u_offset], ldu, &v[v_offset], ldv,
			     &dwork[p * 2]);

#line 550 "MB04YD.f"
		    if (qrit) {
#line 551 "MB04YD.f"
			if ((d__1 = e[k - 1], abs(d__1)) <= tolabs) {
#line 551 "MB04YD.f"
			    e[k - 1] = 0.;
#line 551 "MB04YD.f"
			}
#line 552 "MB04YD.f"
		    } else {
#line 553 "MB04YD.f"
			if ((d__1 = e[i__], abs(d__1)) <= tolabs) {
#line 553 "MB04YD.f"
			    e[i__] = 0.;
#line 553 "MB04YD.f"
			}
#line 554 "MB04YD.f"
		    }
/* Computing 2nd power */
#line 555 "MB04YD.f"
		    d__1 = q[k];
#line 555 "MB04YD.f"
		    dwork[k] = d__1 * d__1;

#line 557 "MB04YD.f"
		    i__1 = k - 1;
#line 557 "MB04YD.f"
		    for (i1 = i__; i1 <= i__1; ++i1) {
/* Computing 2nd power */
#line 558 "MB04YD.f"
			d__1 = q[i1];
#line 558 "MB04YD.f"
			dwork[i1] = d__1 * d__1;
/* Computing 2nd power */
#line 559 "MB04YD.f"
			d__1 = e[i1];
#line 559 "MB04YD.f"
			dwork[p + i1] = d__1 * d__1;
#line 560 "MB04YD.f"
/* L160: */
#line 560 "MB04YD.f"
		    }

#line 562 "MB04YD.f"
		    ++iter;
#line 563 "MB04YD.f"
		}
#line 564 "MB04YD.f"
	    }
#line 565 "MB04YD.f"
	    goto L120;
#line 566 "MB04YD.f"
	}
/*        END WHILE 120 */

#line 569 "MB04YD.f"
	if (iter >= maxit) {
#line 570 "MB04YD.f"
	    *info = 1;
#line 571 "MB04YD.f"
	    goto L200;
#line 572 "MB04YD.f"
	}

#line 574 "MB04YD.f"
	if (x <= tolabs) {

/*           Split at negligible diagonal element ABS( Q(I) ) <= TOLABS. */

#line 578 "MB04YD.f"
	    mb02ny_(&ljobua, &ljobva, m, n, &i__, &k, &q[1], &e[1], &u[
		    u_offset], ldu, &v[v_offset], ldv, &dwork[p * 2]);
#line 580 "MB04YD.f"
	    inul[i__] = TRUE_;
#line 581 "MB04YD.f"
	    --(*rank);
#line 582 "MB04YD.f"
	} else {

/*           A negligible superdiagonal element ABS( E(I-1) ) <= TOL */
/*           has been found, the corresponding bidiagonal submatrix */
/*           belongs to C1 or C2. Treat this bidiagonal submatrix. */

#line 588 "MB04YD.f"
	    if (j >= 2) {
#line 589 "MB04YD.f"
		if (numeig == j) {

#line 591 "MB04YD.f"
		    i__1 = k;
#line 591 "MB04YD.f"
		    for (i1 = i__; i1 <= i__1; ++i1) {
#line 592 "MB04YD.f"
			inul[i1] = TRUE_;
#line 593 "MB04YD.f"
/* L180: */
#line 593 "MB04YD.f"
		    }

#line 595 "MB04YD.f"
		    *rank -= j;
#line 596 "MB04YD.f"
		    k -= j;
#line 597 "MB04YD.f"
		} else {
#line 598 "MB04YD.f"
		    k = i__ - 1;
#line 599 "MB04YD.f"
		}
#line 600 "MB04YD.f"
	    } else {
#line 601 "MB04YD.f"
		if (x <= thetac + tolabs) {
#line 602 "MB04YD.f"
		    inul[i__] = TRUE_;
#line 603 "MB04YD.f"
		    --(*rank);
#line 604 "MB04YD.f"
		}
#line 605 "MB04YD.f"
		--k;
#line 606 "MB04YD.f"
	    }
#line 607 "MB04YD.f"
	    oldk = k;
#line 608 "MB04YD.f"
	}
#line 609 "MB04YD.f"
	goto L80;
#line 610 "MB04YD.f"
    }
/*     END WHILE 80 */

/*     If matrix was scaled, then rescale Q and E appropriately. */

#line 615 "MB04YD.f"
L200:
#line 616 "MB04YD.f"
    if (iascl == 1) {
#line 617 "MB04YD.f"
	d__1 = 1. / sigma;
#line 617 "MB04YD.f"
	dscal_(&p, &d__1, &q[1], &c__1);
#line 618 "MB04YD.f"
	i__1 = p - 1;
#line 618 "MB04YD.f"
	d__1 = 1. / sigma;
#line 618 "MB04YD.f"
	dscal_(&i__1, &d__1, &e[1], &c__1);
#line 619 "MB04YD.f"
    }

#line 621 "MB04YD.f"
    return 0;
/* *** Last line of MB04YD *** */
} /* mb04yd_ */

