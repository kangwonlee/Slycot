#line 1 "MB04XD.f"
/* MB04XD.f -- translated by f2c (version 20100827).
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

#line 1 "MB04XD.f"
/* Table of constant values */

static integer c__6 = 6;
static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b22 = 0.;
static doublereal c_b23 = 1.;

/* Subroutine */ int mb04xd_(char *jobu, char *jobv, integer *m, integer *n, 
	integer *rank, doublereal *theta, doublereal *a, integer *lda, 
	doublereal *u, integer *ldu, doublereal *v, integer *ldv, doublereal *
	q, logical *inul, doublereal *tol, doublereal *reltol, doublereal *
	dwork, integer *ldwork, integer *iwarn, integer *info, ftnlen 
	jobu_len, ftnlen jobv_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, j, k, p, ma, ij;
    static doublereal cs;
    static integer ju, jv;
    static doublereal sn;
    static logical qr;
    static integer pp1;
    static logical all;
    static integer ldw, ldy, itau;
    static doublereal temp;
    extern /* Subroutine */ int mb04yd_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, logical *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlasr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen), dcopy_(integer *, doublereal *, integer *
	    , doublereal *, integer *);
    static integer itaup, itauq;
    extern /* Subroutine */ int mb04xy_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, logical *, integer *, ftnlen,
	     ftnlen);
    static char jobuy[1], jobvy[1];
    static integer jwork;
    static logical wantu, wantv;
    extern /* Subroutine */ int dgebrd_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *), dgeqrf_(integer *, integer *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static logical ljobua, ljobva;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), dlaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen), xerbla_(char *, integer *, ftnlen);
    static logical ljobus, ljobvs;
    static integer ihoush, wrkopt;


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

/*     To compute a basis for the left and/or right singular subspace of */
/*     an M-by-N matrix A corresponding to its smallest singular values. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBU    CHARACTER*1 */
/*             Specifies whether to compute the left singular subspace */
/*             as follows: */
/*             = 'N':  Do not compute the left singular subspace; */
/*             = 'A':  Return the (M - RANK) base vectors of the desired */
/*                     left singular subspace in U; */
/*             = 'S':  Return the first (min(M,N) - RANK) base vectors */
/*                     of the desired left singular subspace in U. */

/*     JOBV    CHARACTER*1 */
/*             Specifies whether to compute the right singular subspace */
/*             as follows: */
/*             = 'N':  Do not compute the right singular subspace; */
/*             = 'A':  Return the (N - RANK) base vectors of the desired */
/*                     right singular subspace in V; */
/*             = 'S':  Return the first (min(M,N) - RANK) base vectors */
/*                     of the desired right singular subspace in V. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows in matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns in matrix A.  N >= 0. */

/*     RANK    (input/output) INTEGER */
/*             On entry, if RANK < 0, then the rank of matrix A is */
/*             computed by the routine as the number of singular values */
/*             greater than THETA. */
/*             Otherwise, RANK must specify the rank of matrix A. */
/*             RANK <= min(M,N). */
/*             On exit, if RANK < 0 on entry, then RANK contains the */
/*             computed rank of matrix A. That is, the number of singular */
/*             values of A greater than THETA. */
/*             Otherwise, the user-supplied value of RANK may be changed */
/*             by the routine on exit if the RANK-th and the (RANK+1)-th */
/*             singular values of A are considered to be equal. */
/*             See also the description of parameter TOL below. */

/*     THETA   (input/output) DOUBLE PRECISION */
/*             On entry, if RANK < 0, then THETA must specify an upper */
/*             bound on the smallest singular values of A corresponding */
/*             to the singular subspace to be computed.  THETA >= 0.0. */
/*             Otherwise, THETA must specify an initial estimate (t say) */
/*             for computing an upper bound on the (min(M,N) - RANK) */
/*             smallest singular values of A. If THETA < 0.0, then t is */
/*             computed by the routine. */
/*             On exit, if RANK >= 0 on entry, then THETA contains the */
/*             computed upper bound such that precisely RANK singular */
/*             values of A are greater than THETA + TOL. */
/*             Otherwise, THETA is unchanged. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading M-by-N part of this array must contain the */
/*             matrix A from which the basis of a desired singular */
/*             subspace is to be computed. */
/*             NOTE that this array is destroyed. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= max(1,M). */

/*     U       (output) DOUBLE PRECISION array, dimension (LDU,*) */
/*             If JOBU = 'A', then the leading M-by-M part of this array */
/*             contains the (M - RANK) M-dimensional base vectors of the */
/*             desired left singular subspace of A corresponding to its */
/*             singular values less than or equal to THETA. These vectors */
/*             are stored in the i-th column(s) of U for which */
/*             INUL(i) = .TRUE., where i = 1,2,...,M. */

/*             If JOBU = 'S', then the leading M-by-min(M,N) part of this */
/*             array contains the first (min(M,N) - RANK) M-dimensional */
/*             base vectors of the desired left singular subspace of A */
/*             corresponding to its singular values less than or equal to */
/*             THETA. These vectors are stored in the i-th column(s) of U */
/*             for which INUL(i) = .TRUE., where i = 1,2,..., min(M,N). */

/*             Otherwise, U is not referenced (since JOBU = 'N') and can */
/*             be supplied as a dummy array (i.e. set parameter LDU = 1 */
/*             and declare this array to be U(1,1) in the calling */
/*             program). */

/*     LDU     INTEGER */
/*             The leading dimension of array U. */
/*             LDU >= max(1,M) if JOBU = 'A' or JOBU = 'S', */
/*             LDU >= 1        if JOBU = 'N'. */

/*     V       (output) DOUBLE PRECISION array, dimension (LDV,*) */
/*             If JOBV = 'A', then the leading N-by-N part of this array */
/*             contains the (N - RANK) N-dimensional base vectors of the */
/*             desired right singular subspace of A corresponding to its */
/*             singular values less than or equal to THETA. These vectors */
/*             are stored in the i-th column(s) of V for which */
/*             INUL(i) = .TRUE., where i = 1,2,...,N. */

/*             If JOBV = 'S', then the leading N-by-min(M,N) part of this */
/*             array contains the first (min(M,N) - RANK) N-dimensional */
/*             base vectors of the desired right singular subspace of A */
/*             corresponding to its singular values less than or equal to */
/*             THETA. These vectors are stored in the i-th column(s) of V */
/*             for which INUL(i) = .TRUE., where i = 1,2,...,MIN( M,N). */

/*             Otherwise, V is not referenced (since JOBV = 'N') and can */
/*             be supplied as a dummy array (i.e. set parameter LDV = 1 */
/*             and declare this array to be V(1,1) in the calling */
/*             program). */

/*     LDV     INTEGER */
/*             The leading dimension of array V. */
/*             LDV >= max(1,N) if JOBV = 'A' or JOBV = 'S', */
/*             LDV >= 1        if JOBV = 'N'. */

/*     Q       (output) DOUBLE PRECISION array, dimension (2*min(M,N)-1) */
/*             This array contains the partially diagonalized bidiagonal */
/*             matrix J computed from A, at the moment that the desired */
/*             singular subspace has been found. Specifically, the */
/*             leading p = min(M,N) entries of Q contain the diagonal */
/*             elements q(1),q(2),...,q(p) and the entries Q(p+1), */
/*             Q(p+2),...,Q(2*p-1) contain the superdiagonal elements */
/*             e(1),e(2),...,e(p-1) of J. */

/*     INUL    (output) LOGICAL array, dimension (max(M,N)) */
/*             If JOBU <> 'N' or JOBV <> 'N', then the indices of the */
/*             elements of this array with value .TRUE. indicate the */
/*             columns in U and/or V containing the base vectors of the */
/*             desired left and/or right singular subspace of A. They */
/*             also equal the indices of the diagonal elements of the */
/*             bidiagonal submatrices in the array Q, which correspond */
/*             to the computed singular subspaces. */

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
/*             equal to 0, then the tolerance is taken as specified in */
/*             SLICOT Library routine MB04YD document. */

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
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK = max(1, LDW + max(2*P + max(M,N), LDY)), where */
/*                  P = min(M,N); */
/*                LDW = max(2*N, N*(N+1)/2), if JOBU <> 'N' and M large */
/*                                                        enough than N; */
/*                LDW = 0,                   otherwise; */
/*                LDY = 8*P - 5, if JOBU <> 'N' or  JOBV <> 'N'; */
/*                LDY = 6*P - 3, if JOBU =  'N' and JOBV =  'N'. */
/*             For optimum performance LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  if the rank of matrix A (as specified by the user) */
/*                   has been lowered because a singular value of */
/*                   multiplicity greater than 1 was found. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the maximum number of QR/QL iteration steps */
/*                   (30*MIN(M,N)) has been exceeded. */

/*     METHOD */

/*     The method used is the Partial Singular Value Decomposition (PSVD) */
/*     approach proposed by Van Huffel, Vandewalle and Haegemans, which */
/*     is an efficient technique (see [1]) for computing the singular */
/*     subspace of a matrix corresponding to its smallest singular */
/*     values. It differs from the classical SVD algorithm [3] at three */
/*     points, which results in high efficiency. Firstly, the Householder */
/*     transformations of the bidiagonalization need only to be applied */
/*     on the base vectors of the desired singular subspaces; secondly, */
/*     the bidiagonal matrix need only be partially diagonalized; and */
/*     thirdly, the convergence rate of the iterative diagonalization can */
/*     be improved by an appropriate choice between QL and QR iterations. */
/*     (Note, however, that LAPACK Library routine DGESVD, for computing */
/*     SVD, also uses either QL and QR iterations.) Depending on the gap, */
/*     the desired numerical accuracy and the dimension of the desired */
/*     singular subspace, the PSVD can be up to three times faster than */
/*     the classical SVD algorithm. */

/*     The PSVD algorithm [1-2] for an M-by-N matrix A proceeds as */
/*     follows: */

/*     Step 1: Bidiagonalization phase */
/*             ----------------------- */
/*      (a) If M is large enough than N, transform A into upper */
/*          triangular form R. */

/*      (b) Transform A (or R) into bidiagonal form: */

/*                |q(1) e(1)  0   ...  0   | */
/*           (0)  | 0   q(2) e(2)      .   | */
/*          J   = | .                  .   | */
/*                | .                e(N-1)| */
/*                | 0            ...  q(N) | */

/*     if M >= N, or */

/*                |q(1)  0    0   ...  0     0   | */
/*           (0)  |e(1) q(2)  0        .     .   | */
/*          J   = | .                  .     .   | */
/*                | .                 q(M-1) .   | */
/*                | 0             ... e(M-1) q(M)| */

/*     if M < N, using Householder transformations. */
/*     In the second case, transform the matrix to the upper bidiagonal */
/*     form by applying Givens rotations. */

/*      (c) If U is requested, initialize U with the identity matrix. */
/*          If V is requested, initialize V with the identity matrix. */

/*     Step 2: Partial diagonalization phase */
/*             ----------------------------- */
/*     If the upper bound THETA is not given, then compute THETA such */
/*     that precisely (min(M,N) - RANK) singular values of the bidiagonal */
/*     matrix are less than or equal to THETA, using a bisection method */
/*     [4]. Diagonalize the given bidiagonal matrix J partially, using */
/*     either QR iterations (if the upper left diagonal element of the */
/*     considered bidiagonal submatrix is larger than the lower right */
/*     diagonal element) or QL iterations, such that J is split into */
/*     unreduced bidiagonal submatrices whose singular values are either */
/*     all larger than THETA or all less than or equal to THETA. */
/*     Accumulate the Givens rotations in U and/or V (if desired). */

/*     Step 3: Back transformation phase */
/*             ------------------------- */
/*      (a) Apply the Householder transformations of Step 1(b) onto the */
/*          columns of U and/or V associated with the bidiagonal */
/*          submatrices with all singular values less than or equal to */
/*          THETA (if U and/or V is desired). */

/*      (b) If M is large enough than N, and U is desired, then apply the */
/*          Householder transformations of Step 1(a) onto each computed */
/*          column of U in Step 3(a). */

/*     REFERENCES */

/*     [1] Van Huffel, S., Vandewalle, J. and Haegemans, A. */
/*         An efficient and reliable algorithm for computing the singular */
/*         subspace of a matrix associated with its smallest singular */
/*         values. */
/*         J. Comput. and Appl. Math., 19, pp. 313-330, 1987. */

/*     [2] Van Huffel, S. */
/*         Analysis of the total least squares problem and its use in */
/*         parameter estimation. */
/*         Doctoral dissertation, Dept. of Electr. Eng., Katholieke */
/*         Universiteit Leuven, Belgium, June 1987. */

/*     [3] Chan, T.F. */
/*         An improved algorithm for computing the singular value */
/*         decomposition. */
/*         ACM TOMS, 8, pp. 72-83, 1982. */

/*     [4] Van Huffel, S. and Vandewalle, J. */
/*         The partial total least squares algorithm. */
/*         J. Comput. and Appl. Math., 21, pp. 333-341, 1988. */

/*     NUMERICAL ASPECTS */

/*     Using the PSVD a large reduction in computation time can be */
/*     gained in total least squares applications (cf [2 - 4]), in the */
/*     computation of the null space of a matrix and in solving */
/*     (non)homogeneous linear equations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1997. */
/*     Supersedes Release 2.0 routine MB04PD by S. Van Huffel, Katholieke */
/*     University Leuven, Belgium. */

/*     REVISIONS */

/*     July 10, 1997. */

/*     KEYWORDS */

/*     Bidiagonalization, singular subspace, singular value */
/*     decomposition, singular values. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 370 "MB04XD.f"
    /* Parameter adjustments */
#line 370 "MB04XD.f"
    a_dim1 = *lda;
#line 370 "MB04XD.f"
    a_offset = 1 + a_dim1;
#line 370 "MB04XD.f"
    a -= a_offset;
#line 370 "MB04XD.f"
    u_dim1 = *ldu;
#line 370 "MB04XD.f"
    u_offset = 1 + u_dim1;
#line 370 "MB04XD.f"
    u -= u_offset;
#line 370 "MB04XD.f"
    v_dim1 = *ldv;
#line 370 "MB04XD.f"
    v_offset = 1 + v_dim1;
#line 370 "MB04XD.f"
    v -= v_offset;
#line 370 "MB04XD.f"
    --q;
#line 370 "MB04XD.f"
    --inul;
#line 370 "MB04XD.f"
    --dwork;
#line 370 "MB04XD.f"

#line 370 "MB04XD.f"
    /* Function Body */
#line 370 "MB04XD.f"
    *iwarn = 0;
#line 371 "MB04XD.f"
    *info = 0;
#line 372 "MB04XD.f"
    p = min(*m,*n);
#line 373 "MB04XD.f"
    k = max(*m,*n);

/*     Determine whether U and/or V are/is to be computed. */

#line 377 "MB04XD.f"
    ljobua = lsame_(jobu, "A", (ftnlen)1, (ftnlen)1);
#line 378 "MB04XD.f"
    ljobus = lsame_(jobu, "S", (ftnlen)1, (ftnlen)1);
#line 379 "MB04XD.f"
    ljobva = lsame_(jobv, "A", (ftnlen)1, (ftnlen)1);
#line 380 "MB04XD.f"
    ljobvs = lsame_(jobv, "S", (ftnlen)1, (ftnlen)1);
#line 381 "MB04XD.f"
    wantu = ljobua || ljobus;
#line 382 "MB04XD.f"
    wantv = ljobva || ljobvs;
#line 383 "MB04XD.f"
    all = ljobua && *m > *n || ljobva && *m < *n;
#line 384 "MB04XD.f"
    qr = *m >= ilaenv_(&c__6, "DGESVD", "NN", m, n, &c__0, &c__0, (ftnlen)6, (
	    ftnlen)2);
#line 385 "MB04XD.f"
    if (qr && wantu) {
/* Computing MAX */
#line 386 "MB04XD.f"
	i__1 = *n << 1, i__2 = *n * (*n + 1) / 2;
#line 386 "MB04XD.f"
	ldw = max(i__1,i__2);
#line 387 "MB04XD.f"
    } else {
#line 388 "MB04XD.f"
	ldw = 0;
#line 389 "MB04XD.f"
    }
#line 390 "MB04XD.f"
    if (wantu || wantv) {
#line 391 "MB04XD.f"
	ldy = (p << 3) - 5;
#line 392 "MB04XD.f"
    } else {
#line 393 "MB04XD.f"
	ldy = p * 6 - 3;
#line 394 "MB04XD.f"
    }

/*     Test the input scalar arguments. */

#line 398 "MB04XD.f"
    if (! wantu && ! lsame_(jobu, "N", (ftnlen)1, (ftnlen)1)) {
#line 399 "MB04XD.f"
	*info = -1;
#line 400 "MB04XD.f"
    } else if (! wantv && ! lsame_(jobv, "N", (ftnlen)1, (ftnlen)1)) {
#line 401 "MB04XD.f"
	*info = -2;
#line 402 "MB04XD.f"
    } else if (*m < 0) {
#line 403 "MB04XD.f"
	*info = -3;
#line 404 "MB04XD.f"
    } else if (*n < 0) {
#line 405 "MB04XD.f"
	*info = -4;
#line 406 "MB04XD.f"
    } else if (*rank > p) {
#line 407 "MB04XD.f"
	*info = -5;
#line 408 "MB04XD.f"
    } else if (*rank < 0 && *theta < 0.) {
#line 409 "MB04XD.f"
	*info = -6;
#line 410 "MB04XD.f"
    } else if (*lda < max(1,*m)) {
#line 411 "MB04XD.f"
	*info = -8;
#line 412 "MB04XD.f"
    } else if (! wantu && *ldu < 1 || wantu && *ldu < max(1,*m)) {
#line 414 "MB04XD.f"
	*info = -10;
#line 415 "MB04XD.f"
    } else if (! wantv && *ldv < 1 || wantv && *ldv < max(1,*n)) {
#line 417 "MB04XD.f"
	*info = -12;
#line 418 "MB04XD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
#line 418 "MB04XD.f"
	i__3 = (p << 1) + k;
#line 418 "MB04XD.f"
	i__1 = 1, i__2 = ldw + max(i__3,ldy);
#line 418 "MB04XD.f"
	if (*ldwork < max(i__1,i__2)) {
#line 419 "MB04XD.f"
	    *info = -18;
#line 420 "MB04XD.f"
	}
#line 420 "MB04XD.f"
    }

#line 422 "MB04XD.f"
    if (*info != 0) {

/*        Error return. */

#line 426 "MB04XD.f"
	i__1 = -(*info);
#line 426 "MB04XD.f"
	xerbla_("MB04XD", &i__1, (ftnlen)6);
#line 427 "MB04XD.f"
	return 0;
#line 428 "MB04XD.f"
    }

/*     Quick return if possible. */

#line 432 "MB04XD.f"
    if (p == 0) {
#line 433 "MB04XD.f"
	if (*rank >= 0) {
#line 433 "MB04XD.f"
	    *theta = 0.;
#line 433 "MB04XD.f"
	}
#line 435 "MB04XD.f"
	*rank = 0;
#line 436 "MB04XD.f"
	return 0;
#line 437 "MB04XD.f"
    }

/*     Initializations. */

#line 441 "MB04XD.f"
    pp1 = p + 1;

#line 443 "MB04XD.f"
    if (all && ! qr) {

#line 445 "MB04XD.f"
	i__1 = p;
#line 445 "MB04XD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 446 "MB04XD.f"
	    inul[i__] = FALSE_;
#line 447 "MB04XD.f"
/* L20: */
#line 447 "MB04XD.f"
	}

#line 449 "MB04XD.f"
	i__1 = k;
#line 449 "MB04XD.f"
	for (i__ = pp1; i__ <= i__1; ++i__) {
#line 450 "MB04XD.f"
	    inul[i__] = TRUE_;
#line 451 "MB04XD.f"
/* L40: */
#line 451 "MB04XD.f"
	}

#line 453 "MB04XD.f"
    } else {

#line 455 "MB04XD.f"
	i__1 = k;
#line 455 "MB04XD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 456 "MB04XD.f"
	    inul[i__] = FALSE_;
#line 457 "MB04XD.f"
/* L60: */
#line 457 "MB04XD.f"
	}

#line 459 "MB04XD.f"
    }

/*     Step 1: Bidiagonalization phase */
/*             ----------------------- */
/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 469 "MB04XD.f"
    if (qr) {

/*        1.a.: M is large enough than N; transform A into upper */
/*              triangular form R by Householder transformations. */

/*        Workspace: need 2*N;  prefer N + N*NB. */

#line 476 "MB04XD.f"
	itau = 1;
#line 477 "MB04XD.f"
	jwork = itau + *n;
#line 478 "MB04XD.f"
	i__1 = *ldwork - jwork + 1;
#line 478 "MB04XD.f"
	dgeqrf_(m, n, &a[a_offset], lda, &dwork[itau], &dwork[jwork], &i__1, 
		info);
#line 480 "MB04XD.f"
	wrkopt = (integer) dwork[jwork] + jwork - 1;

/*        If (WANTU), store information on the Householder */
/*        transformations performed on the columns of A in N*(N+1)/2 */
/*        extra storage locations DWORK(K), for K = 1,2,...,N*(N+1)/2. */
/*        (The first N locations store the scalar factors of Householder */
/*        transformations.) */

/*        Workspace: LDW = max(2*N, N*(N+1)/2). */

#line 490 "MB04XD.f"
	if (wantu) {
#line 491 "MB04XD.f"
	    ihoush = jwork;
#line 492 "MB04XD.f"
	    k = ihoush;
#line 493 "MB04XD.f"
	    i__ = *n;
#line 494 "MB04XD.f"
	} else {
#line 495 "MB04XD.f"
	    k = 1;
#line 496 "MB04XD.f"
	}

#line 498 "MB04XD.f"
	i__1 = *n - 1;
#line 498 "MB04XD.f"
	for (j = 1; j <= i__1; ++j) {
#line 499 "MB04XD.f"
	    if (wantu) {
#line 500 "MB04XD.f"
		--i__;
#line 501 "MB04XD.f"
		dcopy_(&i__, &a[j + 1 + j * a_dim1], &c__1, &dwork[k], &c__1);
#line 502 "MB04XD.f"
		k += i__;
#line 503 "MB04XD.f"
	    }

#line 505 "MB04XD.f"
	    i__2 = *n;
#line 505 "MB04XD.f"
	    for (ij = j + 1; ij <= i__2; ++ij) {
#line 506 "MB04XD.f"
		a[ij + j * a_dim1] = 0.;
#line 507 "MB04XD.f"
/* L80: */
#line 507 "MB04XD.f"
	    }

#line 509 "MB04XD.f"
/* L100: */
#line 509 "MB04XD.f"
	}

#line 511 "MB04XD.f"
	ma = *n;
#line 512 "MB04XD.f"
	wrkopt = max(wrkopt,k);
#line 513 "MB04XD.f"
    } else {

/*        Workspace: LDW = 0. */

#line 517 "MB04XD.f"
	k = 1;
#line 518 "MB04XD.f"
	ma = *m;
#line 519 "MB04XD.f"
	wrkopt = 1;
#line 520 "MB04XD.f"
    }

/*     1.b.: Transform A (or R) into bidiagonal form Q using Householder */
/*           transformations. */

/*     Workspace: need   LDW + 2*min(M,N) + max(M,N); */
/*                prefer LDW + 2*min(M,N) + (M+N)*NB. */

#line 528 "MB04XD.f"
    itauq = k;
#line 529 "MB04XD.f"
    itaup = itauq + p;
#line 530 "MB04XD.f"
    jwork = itaup + p;
#line 531 "MB04XD.f"
    i__1 = *ldwork - jwork + 1;
#line 531 "MB04XD.f"
    dgebrd_(&ma, n, &a[a_offset], lda, &q[1], &q[pp1], &dwork[itauq], &dwork[
	    itaup], &dwork[jwork], &i__1, info);
/* Computing MAX */
#line 533 "MB04XD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 533 "MB04XD.f"
    wrkopt = max(i__1,i__2);

/*     1.c.: Initialize U (if WANTU) and V (if WANTV) with the identity */
/*           matrix. */

#line 538 "MB04XD.f"
    if (wantu) {
#line 539 "MB04XD.f"
	if (all) {
#line 540 "MB04XD.f"
	    ju = *m;
#line 541 "MB04XD.f"
	} else {
#line 542 "MB04XD.f"
	    ju = p;
#line 543 "MB04XD.f"
	}
#line 544 "MB04XD.f"
	dlaset_("Full", m, &ju, &c_b22, &c_b23, &u[u_offset], ldu, (ftnlen)4);
#line 545 "MB04XD.f"
	*(unsigned char *)jobuy = 'U';
#line 546 "MB04XD.f"
    } else {
#line 547 "MB04XD.f"
	*(unsigned char *)jobuy = 'N';
#line 548 "MB04XD.f"
    }
#line 549 "MB04XD.f"
    if (wantv) {
#line 550 "MB04XD.f"
	if (all) {
#line 551 "MB04XD.f"
	    jv = *n;
#line 552 "MB04XD.f"
	} else {
#line 553 "MB04XD.f"
	    jv = p;
#line 554 "MB04XD.f"
	}
#line 555 "MB04XD.f"
	dlaset_("Full", n, &jv, &c_b22, &c_b23, &v[v_offset], ldv, (ftnlen)4);
#line 556 "MB04XD.f"
	*(unsigned char *)jobvy = 'U';
#line 557 "MB04XD.f"
    } else {
#line 558 "MB04XD.f"
	*(unsigned char *)jobvy = 'N';
#line 559 "MB04XD.f"
    }

/*     If the matrix is lower bidiagonal, rotate to be upper bidiagonal */
/*     by applying Givens rotations on the left. */

#line 564 "MB04XD.f"
    if (*m < *n) {

#line 566 "MB04XD.f"
	i__1 = p - 1;
#line 566 "MB04XD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 567 "MB04XD.f"
	    dlartg_(&q[i__], &q[p + i__], &cs, &sn, &temp);
#line 568 "MB04XD.f"
	    q[i__] = temp;
#line 569 "MB04XD.f"
	    q[p + i__] = sn * q[i__ + 1];
#line 570 "MB04XD.f"
	    q[i__ + 1] = cs * q[i__ + 1];
#line 571 "MB04XD.f"
	    if (wantu) {

/*              Workspace: LDW + 4*min(M,N) - 2. */

#line 575 "MB04XD.f"
		dwork[jwork + i__ - 1] = cs;
#line 576 "MB04XD.f"
		dwork[jwork + p + i__ - 2] = sn;
#line 577 "MB04XD.f"
	    }
#line 578 "MB04XD.f"
/* L120: */
#line 578 "MB04XD.f"
	}

/*        Update left singular vectors if desired. */

#line 582 "MB04XD.f"
	if (wantu) {
#line 582 "MB04XD.f"
	    dlasr_("Right", "Variable pivot", "Forward", m, &ju, &dwork[jwork]
		    , &dwork[jwork + p - 1], &u[u_offset], ldu, (ftnlen)5, (
		    ftnlen)14, (ftnlen)7);
#line 582 "MB04XD.f"
	}

#line 586 "MB04XD.f"
    }

/*     Step 2: Partial diagonalization phase. */
/*             ----------------------------- */
/*             Diagonalize the bidiagonal Q partially until convergence */
/*             to  the desired left and/or right singular subspace. */

/*              Workspace: LDW + 8*min(M,N) - 5, if WANTU or WANTV; */
/*              Workspace: LDW + 6*min(M,N) - 3, if JOBU = JOBV = 'N'. */

#line 596 "MB04XD.f"
    i__1 = *ldwork - jwork + 1;
#line 596 "MB04XD.f"
    mb04yd_(jobuy, jobvy, m, n, rank, theta, &q[1], &q[pp1], &u[u_offset], 
	    ldu, &v[v_offset], ldv, &inul[1], tol, reltol, &dwork[jwork], &
	    i__1, iwarn, info, (ftnlen)1, (ftnlen)1);
#line 599 "MB04XD.f"
    if (wantu || wantv) {
/* Computing MAX */
#line 600 "MB04XD.f"
	i__1 = wrkopt, i__2 = jwork - 6 + (p << 3);
#line 600 "MB04XD.f"
	wrkopt = max(i__1,i__2);
#line 601 "MB04XD.f"
    } else {
/* Computing MAX */
#line 602 "MB04XD.f"
	i__1 = wrkopt, i__2 = jwork - 4 + p * 6;
#line 602 "MB04XD.f"
	wrkopt = max(i__1,i__2);
#line 603 "MB04XD.f"
    }
#line 604 "MB04XD.f"
    if (*info > 0) {
#line 604 "MB04XD.f"
	return 0;
#line 604 "MB04XD.f"
    }

/*     Step 3: Back transformation phase. */
/*             ------------------------- */
/*     3.a.: Apply the Householder transformations of the bidiagonaliza- */
/*           tion onto the base vectors associated with the desired */
/*           bidiagonal submatrices. */

/*           Workspace: LDW + 2*min(M,N). */

#line 615 "MB04XD.f"
    mb04xy_(jobu, jobv, &ma, n, &a[a_offset], lda, &dwork[itauq], &dwork[
	    itaup], &u[u_offset], ldu, &v[v_offset], ldv, &inul[1], info, (
	    ftnlen)1, (ftnlen)1);

/*     3.b.: If A was reduced to upper triangular form R and JOBU = 'A' */
/*           or JOBU = 'S' apply the Householder transformations of the */
/*           triangularization of A onto the desired base vectors. */

#line 622 "MB04XD.f"
    if (qr && wantu) {
#line 623 "MB04XD.f"
	if (all) {

#line 625 "MB04XD.f"
	    i__1 = *m;
#line 625 "MB04XD.f"
	    for (i__ = pp1; i__ <= i__1; ++i__) {
#line 626 "MB04XD.f"
		inul[i__] = TRUE_;
#line 627 "MB04XD.f"
/* L140: */
#line 627 "MB04XD.f"
	    }

#line 629 "MB04XD.f"
	}
#line 630 "MB04XD.f"
	k = ihoush;
#line 631 "MB04XD.f"
	i__ = *n;

#line 633 "MB04XD.f"
	i__1 = *n - 1;
#line 633 "MB04XD.f"
	for (j = 1; j <= i__1; ++j) {
#line 634 "MB04XD.f"
	    --i__;
#line 635 "MB04XD.f"
	    dcopy_(&i__, &dwork[k], &c__1, &a[j + 1 + j * a_dim1], &c__1);
#line 636 "MB04XD.f"
	    k += i__;
#line 637 "MB04XD.f"
/* L160: */
#line 637 "MB04XD.f"
	}

/*        Workspace: MIN(M,N) + 1. */

#line 641 "MB04XD.f"
	jwork = pp1;
#line 642 "MB04XD.f"
	mb04xy_(jobu, "No V", m, n, &a[a_offset], lda, &dwork[itau], &dwork[
		itau], &u[u_offset], ldu, &dwork[jwork], &c__1, &inul[1], 
		info, (ftnlen)1, (ftnlen)4);
#line 644 "MB04XD.f"
	wrkopt = max(wrkopt,pp1);
#line 645 "MB04XD.f"
    }

/*     Set the optimal workspace. */

#line 649 "MB04XD.f"
    dwork[1] = (doublereal) wrkopt;
#line 650 "MB04XD.f"
    return 0;
/* *** Last line of MB04XD *** */
} /* mb04xd_ */

