#line 1 "MB02ND.f"
/* MB02ND.f -- translated by f2c (version 20100827).
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

#line 1 "MB02ND.f"
/* Table of constant values */

static doublereal c_b4 = 0.;
static doublereal c_b5 = 1.;
static integer c__6 = 6;
static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b85 = -1.;

/* Subroutine */ int mb02nd_(integer *m, integer *n, integer *l, integer *
	rank, doublereal *theta, doublereal *c__, integer *ldc, doublereal *x,
	 integer *ldx, doublereal *q, logical *inul, doublereal *tol, 
	doublereal *reltol, integer *iwork, doublereal *dwork, integer *
	ldwork, logical *bwork, integer *iwarn, integer *info)
{
    /* System generated locals */
    integer c_dim1, c_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j, k, p, i1, j1, n1, jf, kf, mc, ij;
    static doublereal hh;
    static integer kj;
    static doublereal cs;
    static integer mj, nj, nl, jv;
    static doublereal sn;
    static integer lw, ldf, mnl;
    static doublereal eps;
    static integer ioff;
    static doublereal temp;
    static integer ifail;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen), mb04yd_(char *, char *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, ftnlen, ftnlen);
    static doublereal rcond;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer iwarm;
    static doublereal fnorm;
    static integer itaup, itauq;
    static doublereal first;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer jwork;
    static doublereal dummy[1];
    extern /* Subroutine */ int dgebrd_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), dgeqrf_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dgerqf_(integer *, integer *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *, integer *), dlaset_(char *
	    , integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen), dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dlacpy_(
	    char *, integer *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, ftnlen);
    extern doublereal dlantr_(char *, char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    extern /* Subroutine */ int dormbr_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen), dtrcon_(char *, char *, char *, integer *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static doublereal inprod;
    static integer ihoush;
    static logical lfirst;
    extern /* Subroutine */ int dormrq_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static logical sufwrk;
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

/*     To solve the Total Least Squares (TLS) problem using a Partial */
/*     Singular Value Decomposition (PSVD) approach. */
/*     The TLS problem assumes an overdetermined set of linear equations */
/*     AX = B, where both the data matrix A as well as the observation */
/*     matrix B are inaccurate. The routine also solves determined and */
/*     underdetermined sets of equations by computing the minimum norm */
/*     solution. */
/*     It is assumed that all preprocessing measures (scaling, coordinate */
/*     transformations, whitening, ... ) of the data have been performed */
/*     in advance. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows in the data matrix A and the */
/*             observation matrix B.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns in the data matrix A.  N >= 0. */

/*     L       (input) INTEGER */
/*             The number of columns in the observation matrix B. */
/*             L >= 0. */

/*     RANK    (input/output) INTEGER */
/*             On entry, if RANK < 0, then the rank of the TLS */
/*             approximation [A+DA|B+DB] (r say) is computed by the */
/*             routine. */
/*             Otherwise, RANK must specify the value of r. */
/*             RANK <= min(M,N). */
/*             On exit, if RANK < 0 on entry and INFO = 0, then RANK */
/*             contains the computed rank of the TLS approximation */
/*             [A+DA|B+DB]. */
/*             Otherwise, the user-supplied value of RANK may be */
/*             changed by the routine on exit if the RANK-th and the */
/*             (RANK+1)-th singular values of C = [A|B] are considered */
/*             to be equal, or if the upper triangular matrix F (as */
/*             defined in METHOD) is (numerically) singular. */

/*     THETA   (input/output) DOUBLE PRECISION */
/*             On entry, if RANK < 0, then the rank of the TLS */
/*             approximation [A+DA|B+DB] is computed using THETA as */
/*             (min(M,N+L) - d), where d is the number of singular */
/*             values of [A|B] <= THETA. THETA >= 0.0. */
/*             Otherwise, THETA is an initial estimate (t say) for */
/*             computing a lower bound on the RANK largest singular */
/*             values of [A|B]. If THETA < 0.0 on entry however, then */
/*             t is computed by the routine. */
/*             On exit, if RANK >= 0 on entry, then THETA contains the */
/*             computed bound such that precisely RANK singular values */
/*             of C = [A|B] are greater than THETA + TOL. */
/*             Otherwise, THETA is unchanged. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N+L) */
/*             On entry, the leading M-by-(N+L) part of this array must */
/*             contain the matrices A and B. Specifically, the first N */
/*             columns must contain the data matrix A and the last L */
/*             columns the observation matrix B (right-hand sides). */
/*             On exit, if INFO = 0, the first N+L components of the */
/*             columns of this array whose index i corresponds with */
/*             INUL(i) = .TRUE., are the possibly transformed (N+L-RANK) */
/*             base vectors of the right singular subspace corresponding */
/*             to the singular values of C = [A|B] which are less than or */
/*             equal to THETA. Specifically, if L = 0, or if RANK = 0 and */
/*             IWARN <> 2, these vectors are indeed the base vectors */
/*             above. Otherwise, these vectors form the matrix V2, */
/*             transformed as described in Step 4 of the PTLS algorithm */
/*             (see METHOD). The TLS solution is computed from these */
/*             vectors. The other columns of array C contain no useful */
/*             information. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= max(1,M,N+L). */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,L) */
/*             If INFO = 0, the leading N-by-L part of this array */
/*             contains the solution X to the TLS problem specified by */
/*             A and B. */

/*     LDX     INTEGER */
/*             The leading dimension of array X.  LDX >= max(1,N). */

/*     Q       (output) DOUBLE PRECISION array, dimension */
/*             (max(1,2*min(M,N+L)-1)) */
/*             This array contains the partially diagonalized bidiagonal */
/*             matrix J computed from C, at the moment that the desired */
/*             singular subspace has been found. Specifically, the */
/*             leading p = min(M,N+L) entries of Q contain the diagonal */
/*             elements q(1),q(2),...,q(p) and the entries Q(p+1),Q(p+2), */
/*             ...,Q(2*p-1) contain the superdiagonal elements e(1),e(2), */
/*             ...,e(p-1) of J. */

/*     INUL    (output) LOGICAL array, dimension (N+L) */
/*             The indices of the elements of this array with value */
/*             .TRUE. indicate the columns in C containing the base */
/*             vectors of the right singular subspace of C from which */
/*             the TLS solution has been computed. */

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

/*     IWORK   INTEGER array, dimension (N+2*L) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK, and DWORK(2) returns the reciprocal of the */
/*             condition number of the matrix F. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK = max(2, max(M,N+L) + 2*min(M,N+L), */
/*                          min(M,N+L) + LW + max(6*(N+L)-5, */
/*                                                L*L+max(N+L,3*L)), */
/*             where */
/*             LW = (N+L)*(N+L-1)/2,  if M >= N+L, */
/*             LW = M*(N+L-(M-1)/2),  if M <  N+L. */
/*             For optimum performance LDWORK should be larger. */

/*     BWORK   LOGICAL array, dimension (N+L) */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warnings; */
/*             = 1:  if the rank of matrix C has been lowered because a */
/*                   singular value of multiplicity greater than 1 was */
/*                   found; */
/*             = 2:  if the rank of matrix C has been lowered because the */
/*                   upper triangular matrix F is (numerically) singular. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the maximum number of QR/QL iteration steps */
/*                   (30*MIN(M,N)) has been exceeded; */
/*             = 2:  if the computed rank of the TLS approximation */
/*                   [A+DA|B+DB] exceeds MIN(M,N). Try increasing the */
/*                   value of THETA or set the value of RANK to min(M,N). */

/*     METHOD */

/*     The method used is the Partial Total Least Squares (PTLS) approach */
/*     proposed by Van Huffel and Vandewalle [5]. */

/*     Let C = [A|B] denote the matrix formed by adjoining the columns of */
/*     B to the columns of A on the right. */

/*     Total Least Squares (TLS) definition: */
/*     ------------------------------------- */

/*       Given matrices A and B, find a matrix X satisfying */

/*            (A + DA) X = B + DB, */

/*       where A and DA are M-by-N matrices, B and DB are M-by-L matrices */
/*       and X is an N-by-L matrix. */
/*       The solution X must be such that the Frobenius norm of [DA|DB] */
/*       is a minimum and each column of B + DB is in the range of */
/*       A + DA. Whenever the solution is not unique, the routine singles */
/*       out the minimum norm solution X. */

/*     Let V denote the right singular subspace of C. Since the TLS */
/*     solution can be computed from any orthogonal basis of the subspace */
/*     of V corresponding to the smallest singular values of C, the */
/*     Partial Singular Value Decomposition (PSVD) can be used instead of */
/*     the classical SVD. The dimension of this subspace of V may be */
/*     determined by the rank of C or by an upper bound for those */
/*     smallest singular values. */

/*     The PTLS algorithm proceeds as follows (see [2 - 5]): */

/*     Step 1: Bidiagonalization phase */
/*             ----------------------- */
/*      (a) If M is large enough than N + L, transform C into upper */
/*          triangular form R by Householder transformations. */
/*      (b) Transform C (or R) into upper bidiagonal form */
/*          (p = min(M,N+L)): */

/*                     |q(1) e(1)  0   ...  0   | */
/*                (0)  | 0   q(2) e(2)      .   | */
/*               J   = | .                  .   | */
/*                     | .                e(p-1)| */
/*                     | 0             ... q(p) | */

/*          if M >= N + L, or lower bidiagonal form: */

/*                     |q(1)  0    0   ...  0     0   | */
/*                (0)  |e(1) q(2)  0        .     .   | */
/*               J   = | .                  .     .   | */
/*                     | .                 q(p)   .   | */
/*                     | 0             ... e(p-1) q(p)| */

/*          if M < N + L, using Householder transformations. */
/*          In the second case, transform the matrix to the upper */
/*          bidiagonal form by applying Givens rotations. */
/*      (c) Initialize the right singular base matrix with the identity */
/*          matrix. */

/*     Step 2: Partial diagonalization phase */
/*             ----------------------------- */
/*     If the upper bound THETA is not given, then compute THETA such */
/*     that precisely p - RANK singular values (p=min(M,N+L)) of the */
/*     bidiagonal matrix are less than or equal to THETA, using a */
/*     bisection method [5]. Diagonalize the given bidiagonal matrix J */
/*     partially, using either QL iterations (if the upper left diagonal */
/*     element of the considered bidiagonal submatrix is smaller than the */
/*     lower right diagonal element) or QR iterations, such that J is */
/*     split into unreduced bidiagonal submatrices whose singular values */
/*     are either all larger than THETA or are all less than or equal */
/*     to THETA. Accumulate the Givens rotations in V. */

/*     Step 3: Back transformation phase */
/*             ------------------------- */
/*     Apply the Householder transformations of Step 1(b) onto the base */
/*     vectors of V associated with the bidiagonal submatrices with all */
/*     singular values less than or equal to THETA. */

/*     Step 4: Computation of F and Y */
/*             ---------------------- */
/*     Let V2 be the matrix of the columns of V corresponding to the */
/*     (N + L - RANK) smallest singular values of C. */
/*     Compute with Householder transformations the matrices F and Y */
/*     such that: */

/*                       |VH   Y| */
/*              V2 x Q = |      | */
/*                       |0    F| */

/*     where Q is an orthogonal matrix, VH is an N-by-(N-RANK) matrix, */
/*     Y is an N-by-L matrix and F is an L-by-L upper triangular matrix. */
/*     If F is singular, then reduce the value of RANK by one and repeat */
/*     Steps 2, 3 and 4. */

/*     Step 5: Computation of the TLS solution */
/*             ------------------------------- */
/*     If F is non-singular then the solution X is obtained by solving */
/*     the following equations by forward elimination: */

/*              X F = -Y. */

/*     Notes: */
/*     If RANK is lowered in Step 4, some additional base vectors must */
/*     be computed in Step 2. The additional computations are kept to */
/*     a minimum. */
/*     If RANK is lowered in Step 4 but the multiplicity of the RANK-th */
/*     singular value is larger than 1, then the value of RANK is further */
/*     lowered with its multiplicity defined by the parameter TOL. This */
/*     is done at the beginning of Step 2 by calling SLICOT Library */
/*     routine MB03MD (from MB04YD), which estimates THETA using a */
/*     bisection method. If F in Step 4 is singular, then the computed */
/*     solution is infinite and hence does not satisfy the second TLS */
/*     criterion (see TLS definition). For these cases, Golub and */
/*     Van Loan [1] claim that the TLS problem has no solution. The */
/*     properties of these so-called nongeneric problems are described */
/*     in [6] and the TLS computations are generalized in order to solve */
/*     them. As proven in [6], the proposed generalization satisfies the */
/*     TLS criteria for any number L of observation vectors in B provided */
/*     that, in addition, the solution | X| is constrained to be */
/*                                     |-I| */
/*     orthogonal to all vectors of the form |w| which belong to the */
/*                                           |0| */
/*     space generated by the columns of the submatrix |Y|. */
/*                                                     |F| */

/*     REFERENCES */

/*     [1] Golub, G.H. and Van Loan, C.F. */
/*         An Analysis of the Total Least-Squares Problem. */
/*         SIAM J. Numer. Anal., 17, pp. 883-893, 1980. */

/*     [2] Van Huffel, S., Vandewalle, J. and Haegemans, A. */
/*         An Efficient and Reliable Algorithm for Computing the */
/*         Singular Subspace of a Matrix Associated with its Smallest */
/*         Singular Values. */
/*         J. Comput. and Appl. Math., 19, pp. 313-330, 1987. */

/*     [3] Van Huffel, S. */
/*         Analysis of the Total Least Squares Problem and its Use in */
/*         Parameter Estimation. */
/*         Doctoral dissertation, Dept. of Electr. Eng., Katholieke */
/*         Universiteit Leuven, Belgium, June 1987. */

/*     [4] Chan, T.F. */
/*         An Improved Algorithm for Computing the Singular Value */
/*         Decomposition. */
/*         ACM TOMS, 8, pp. 72-83, 1982. */

/*     [5] Van Huffel, S. and Vandewalle, J. */
/*         The Partial Total Least Squares Algorithm. */
/*         J. Comput. Appl. Math., 21, pp. 333-341, 1988. */

/*     [6] Van Huffel, S. and Vandewalle, J. */
/*         Analysis and Solution of the Nongeneric Total Least Squares */
/*         Problem. */
/*         SIAM J. Matr. Anal. and Appl., 9, pp. 360-372, 1988. */

/*     NUMERICAL ASPECTS */

/*     The computational efficiency of the PTLS algorithm compared with */
/*     the classical TLS algorithm (see [2 - 5]) is obtained by making */
/*     use of PSVD (see [1]) instead of performing the entire SVD. */
/*     Depending on the gap between the RANK-th and the (RANK+1)-th */
/*     singular values of C, the number (N + L - RANK) of base vectors to */
/*     be computed with respect to the column dimension (N + L) of C and */
/*     the desired accuracy RELTOL, the algorithm used by this routine is */
/*     approximately twice as fast as the classical TLS algorithm at the */
/*     expense of extra storage requirements, namely: */
/*       (N + L) x (N + L - 1)/2  if M >= N + L or */
/*       M x (N + L - (M - 1)/2)  if M <  N + L. */
/*     This is because the Householder transformations performed on the */
/*     rows of C in the bidiagonalization phase (see Step 1) must be kept */
/*     until the end (Step 5). */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997. */
/*     Supersedes Release 2.0 routine MB02BD by S. Van Huffel, Katholieke */
/*     University, Leuven, Belgium. */

/*     REVISIONS */

/*     June 30, 1997, Oct. 19, 2003, Feb. 15, 2004. */

/*     KEYWORDS */

/*     Least-squares approximation, singular subspace, singular value */
/*     decomposition, singular values, total least-squares. */

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

#line 416 "MB02ND.f"
    /* Parameter adjustments */
#line 416 "MB02ND.f"
    c_dim1 = *ldc;
#line 416 "MB02ND.f"
    c_offset = 1 + c_dim1;
#line 416 "MB02ND.f"
    c__ -= c_offset;
#line 416 "MB02ND.f"
    x_dim1 = *ldx;
#line 416 "MB02ND.f"
    x_offset = 1 + x_dim1;
#line 416 "MB02ND.f"
    x -= x_offset;
#line 416 "MB02ND.f"
    --q;
#line 416 "MB02ND.f"
    --inul;
#line 416 "MB02ND.f"
    --iwork;
#line 416 "MB02ND.f"
    --dwork;
#line 416 "MB02ND.f"
    --bwork;
#line 416 "MB02ND.f"

#line 416 "MB02ND.f"
    /* Function Body */
#line 416 "MB02ND.f"
    *iwarn = 0;
#line 417 "MB02ND.f"
    *info = 0;
#line 418 "MB02ND.f"
    nl = *n + *l;
#line 419 "MB02ND.f"
    k = max(*m,nl);
#line 420 "MB02ND.f"
    p = min(*m,nl);
#line 421 "MB02ND.f"
    if (*m >= nl) {
#line 422 "MB02ND.f"
	lw = nl * (nl - 1) / 2;
#line 423 "MB02ND.f"
    } else {
#line 424 "MB02ND.f"
	lw = *m * nl - *m * (*m - 1) / 2;
#line 425 "MB02ND.f"
    }
/* Computing MAX */
/* Computing MAX */
#line 426 "MB02ND.f"
    i__3 = nl, i__4 = *l * 3;
#line 426 "MB02ND.f"
    i__1 = nl * 6 - 5, i__2 = *l * *l + max(i__3,i__4);
#line 426 "MB02ND.f"
    jv = p + lw + max(i__1,i__2);

/*     Test the input scalar arguments. */

#line 430 "MB02ND.f"
    if (*m < 0) {
#line 431 "MB02ND.f"
	*info = -1;
#line 432 "MB02ND.f"
    } else if (*n < 0) {
#line 433 "MB02ND.f"
	*info = -2;
#line 434 "MB02ND.f"
    } else if (*l < 0) {
#line 435 "MB02ND.f"
	*info = -3;
#line 436 "MB02ND.f"
    } else if (*rank > min(*m,*n)) {
#line 437 "MB02ND.f"
	*info = -4;
#line 438 "MB02ND.f"
    } else if (*rank < 0 && *theta < 0.) {
#line 439 "MB02ND.f"
	*info = -5;
#line 440 "MB02ND.f"
    } else if (*ldc < max(1,k)) {
#line 441 "MB02ND.f"
	*info = -7;
#line 442 "MB02ND.f"
    } else if (*ldx < max(1,*n)) {
#line 443 "MB02ND.f"
	*info = -9;
#line 444 "MB02ND.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 444 "MB02ND.f"
	i__1 = 2, i__2 = k + (p << 1), i__1 = max(i__1,i__2);
#line 444 "MB02ND.f"
	if (*ldwork < max(i__1,jv)) {
#line 445 "MB02ND.f"
	    *info = -16;
#line 446 "MB02ND.f"
	}
#line 446 "MB02ND.f"
    }

#line 448 "MB02ND.f"
    if (*info != 0) {

/*        Error return. */

#line 452 "MB02ND.f"
	i__1 = -(*info);
#line 452 "MB02ND.f"
	xerbla_("MB02ND", &i__1, (ftnlen)6);
#line 453 "MB02ND.f"
	return 0;
#line 454 "MB02ND.f"
    }

/*     Quick return if possible. */

#line 458 "MB02ND.f"
    if (min(*m,nl) == 0) {
#line 459 "MB02ND.f"
	if (*m == 0) {
#line 460 "MB02ND.f"
	    dlaset_("Full", &nl, &nl, &c_b4, &c_b5, &c__[c_offset], ldc, (
		    ftnlen)4);
#line 461 "MB02ND.f"
	    dlaset_("Full", n, l, &c_b4, &c_b4, &x[x_offset], ldx, (ftnlen)4);

#line 463 "MB02ND.f"
	    i__1 = nl;
#line 463 "MB02ND.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 464 "MB02ND.f"
		inul[i__] = TRUE_;
#line 465 "MB02ND.f"
/* L10: */
#line 465 "MB02ND.f"
	    }

#line 467 "MB02ND.f"
	}
#line 468 "MB02ND.f"
	if (*rank >= 0) {
#line 468 "MB02ND.f"
	    *theta = 0.;
#line 468 "MB02ND.f"
	}
#line 470 "MB02ND.f"
	*rank = 0;
#line 471 "MB02ND.f"
	dwork[1] = 2.;
#line 472 "MB02ND.f"
	dwork[2] = 1.;
#line 473 "MB02ND.f"
	return 0;
#line 474 "MB02ND.f"
    }

#line 476 "MB02ND.f"
    wrkopt = 2;
#line 477 "MB02ND.f"
    n1 = *n + 1;

#line 479 "MB02ND.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 480 "MB02ND.f"
    lfirst = TRUE_;

/*     Initializations. */

#line 484 "MB02ND.f"
    i__1 = p;
#line 484 "MB02ND.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 485 "MB02ND.f"
	inul[i__] = FALSE_;
#line 486 "MB02ND.f"
	bwork[i__] = FALSE_;
#line 487 "MB02ND.f"
/* L20: */
#line 487 "MB02ND.f"
    }

#line 489 "MB02ND.f"
    i__1 = nl;
#line 489 "MB02ND.f"
    for (i__ = p + 1; i__ <= i__1; ++i__) {
#line 490 "MB02ND.f"
	inul[i__] = TRUE_;
#line 491 "MB02ND.f"
	bwork[i__] = FALSE_;
#line 492 "MB02ND.f"
/* L40: */
#line 492 "MB02ND.f"
    }

/*     Subroutine MB02ND solves a set of linear equations by a Total */
/*     Least Squares Approximation, based on the Partial SVD. */

/*     Step 1: Bidiagonalization phase */
/*             ----------------------- */
/*     1.a): If M is large enough than N+L, transform C into upper */
/*           triangular form R by Householder transformations. */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

/* Computing MAX */
#line 508 "MB02ND.f"
    i__1 = nl, i__2 = ilaenv_(&c__6, "DGESVD", "NN", m, &nl, &c__0, &c__0, (
	    ftnlen)6, (ftnlen)2);
#line 508 "MB02ND.f"
    if (*m >= max(i__1,i__2)) {

/*        Workspace: need   2*(N+L), */
/*                   prefer N+L + (N+L)*NB. */

#line 515 "MB02ND.f"
	itauq = 1;
#line 516 "MB02ND.f"
	jwork = itauq + nl;
#line 517 "MB02ND.f"
	i__1 = *ldwork - jwork + 1;
#line 517 "MB02ND.f"
	dgeqrf_(m, &nl, &c__[c_offset], ldc, &dwork[itauq], &dwork[jwork], &
		i__1, &ifail);
/* Computing MAX */
#line 519 "MB02ND.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 519 "MB02ND.f"
	wrkopt = max(i__1,i__2);
#line 520 "MB02ND.f"
	if (nl > 1) {
#line 520 "MB02ND.f"
	    i__1 = nl - 1;
#line 520 "MB02ND.f"
	    i__2 = nl - 1;
#line 520 "MB02ND.f"
	    dlaset_("Lower", &i__1, &i__2, &c_b4, &c_b4, &c__[c_dim1 + 2], 
		    ldc, (ftnlen)5);
#line 520 "MB02ND.f"
	}
#line 522 "MB02ND.f"
	mnl = nl;
#line 523 "MB02ND.f"
    } else {
#line 524 "MB02ND.f"
	mnl = *m;
#line 525 "MB02ND.f"
    }

/*     1.b): Transform C (or R) into bidiagonal form Q using Householder */
/*           transformations. */
/*     Workspace: need   2*min(M,N+L) + max(M,N+L), */
/*                prefer 2*min(M,N+L) + (M+N+L)*NB. */

#line 532 "MB02ND.f"
    itaup = 1;
#line 533 "MB02ND.f"
    itauq = itaup + p;
#line 534 "MB02ND.f"
    jwork = itauq + p;
#line 535 "MB02ND.f"
    i__1 = *ldwork - jwork + 1;
#line 535 "MB02ND.f"
    dgebrd_(&mnl, &nl, &c__[c_offset], ldc, &q[1], &q[p + 1], &dwork[itauq], &
	    dwork[itaup], &dwork[jwork], &i__1, &ifail);
/* Computing MAX */
#line 537 "MB02ND.f"
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 537 "MB02ND.f"
    wrkopt = max(i__1,i__2);

/*     If the matrix is lower bidiagonal, rotate to be upper bidiagonal */
/*     by applying Givens rotations on the left. */

#line 542 "MB02ND.f"
    if (*m < nl) {
#line 543 "MB02ND.f"
	ioff = 0;

#line 545 "MB02ND.f"
	i__1 = p - 1;
#line 545 "MB02ND.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 546 "MB02ND.f"
	    dlartg_(&q[i__], &q[p + i__], &cs, &sn, &temp);
#line 547 "MB02ND.f"
	    q[i__] = temp;
#line 548 "MB02ND.f"
	    q[p + i__] = sn * q[i__ + 1];
#line 549 "MB02ND.f"
	    q[i__ + 1] = cs * q[i__ + 1];
#line 550 "MB02ND.f"
/* L60: */
#line 550 "MB02ND.f"
	}

#line 552 "MB02ND.f"
    } else {
#line 553 "MB02ND.f"
	ioff = 1;
#line 554 "MB02ND.f"
    }

/*     Store the Householder transformations performed onto the rows of C */
/*     in the extra storage locations DWORK(IHOUSH). */
/*     Workspace: need   LDW = min(M,N+L) + (N+L)*(N+L-1)/2, if M >= N+L, */
/*                       LDW = min(M,N+L) + M*(N+L-(M-1)/2), if M <  N+L; */
/*                prefer LDW = min(M,N+L) + (N+L)**2,        if M >= N+L, */
/*                       LDW = min(M,N+L) + M*(N+L),         if M <  N+L. */

#line 563 "MB02ND.f"
    ihoush = itauq;
#line 564 "MB02ND.f"
    mc = nl - ioff;
#line 565 "MB02ND.f"
    kf = ihoush + p * nl;
/* Computing MAX */
/* Computing 2nd power */
#line 566 "MB02ND.f"
    i__3 = nl;
/* Computing MAX */
#line 566 "MB02ND.f"
    i__4 = nl, i__5 = *l * 3;
#line 566 "MB02ND.f"
    i__1 = (*n + *l) * 6 - 5, i__2 = i__3 * i__3 + max(i__4,i__5) - 1;
#line 566 "MB02ND.f"
    sufwrk = *ldwork >= kf + max(i__1,i__2);
#line 568 "MB02ND.f"
    if (sufwrk) {

/*        Enough workspace for a fast algorithm. */

#line 572 "MB02ND.f"
	dlacpy_("Upper", &p, &nl, &c__[c_offset], ldc, &dwork[ihoush], &p, (
		ftnlen)5);
#line 573 "MB02ND.f"
	kj = kf;
/* Computing MAX */
#line 574 "MB02ND.f"
	i__1 = wrkopt, i__2 = kf - 1;
#line 574 "MB02ND.f"
	wrkopt = max(i__1,i__2);
#line 575 "MB02ND.f"
    } else {

/*        Not enough workspace for a fast algorithm. */

#line 579 "MB02ND.f"
	kj = ihoush;

#line 581 "MB02ND.f"
	i__1 = min(p,mc);
#line 581 "MB02ND.f"
	for (nj = 1; nj <= i__1; ++nj) {
#line 582 "MB02ND.f"
	    j = mc - nj + 1;
#line 583 "MB02ND.f"
	    dcopy_(&j, &c__[nj + (nj + ioff) * c_dim1], ldc, &dwork[kj], &
		    c__1);
#line 584 "MB02ND.f"
	    kj += j;
#line 585 "MB02ND.f"
/* L80: */
#line 585 "MB02ND.f"
	}

#line 587 "MB02ND.f"
    }

/*     1.c): Initialize the right singular base matrix V with the */
/*           identity matrix (V overwrites C). */

#line 592 "MB02ND.f"
    dlaset_("Full", &nl, &nl, &c_b4, &c_b5, &c__[c_offset], ldc, (ftnlen)4);
#line 593 "MB02ND.f"
    jv = kj;
#line 594 "MB02ND.f"
    iwarm = 0;

/*     REPEAT */

/*     Compute the Householder matrix Q and matrices F and Y such that */
/*     F is nonsingular. */

/*     Step 2: Partial diagonalization phase. */
/*             ----------------------------- */
/*     Diagonalize the bidiagonal Q partially until convergence to */
/*     the desired right singular subspace. */
/*     Workspace: LDW + 6*(N+L)-5. */

#line 607 "MB02ND.f"
L100:
#line 608 "MB02ND.f"
    jwork = jv;
#line 609 "MB02ND.f"
    i__1 = *ldwork - jwork + 1;
#line 609 "MB02ND.f"
    mb04yd_("No U", "Update V", &p, &nl, rank, theta, &q[1], &q[p + 1], dummy,
	     &c__1, &c__[c_offset], ldc, &inul[1], tol, reltol, &dwork[jwork],
	     &i__1, iwarn, info, (ftnlen)4, (ftnlen)8);
/* Computing MAX */
#line 612 "MB02ND.f"
    i__1 = wrkopt, i__2 = jwork + nl * 6 - 6;
#line 612 "MB02ND.f"
    wrkopt = max(i__1,i__2);

#line 614 "MB02ND.f"
    *iwarn = max(*iwarn,iwarm);
#line 615 "MB02ND.f"
    if (*info > 0) {
#line 615 "MB02ND.f"
	return 0;
#line 615 "MB02ND.f"
    }

/*     Set pointers to the selected base vectors in the right singular */
/*     matrix of C. */

#line 621 "MB02ND.f"
    k = 0;

#line 623 "MB02ND.f"
    i__1 = nl;
#line 623 "MB02ND.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 624 "MB02ND.f"
	if (inul[i__]) {
#line 625 "MB02ND.f"
	    ++k;
#line 626 "MB02ND.f"
	    iwork[k] = i__;
#line 627 "MB02ND.f"
	}
#line 628 "MB02ND.f"
/* L120: */
#line 628 "MB02ND.f"
    }

#line 630 "MB02ND.f"
    if (k < *l) {

/*        Rank of the TLS approximation is larger than min(M,N). */

#line 634 "MB02ND.f"
	*info = 2;
#line 635 "MB02ND.f"
	return 0;
#line 636 "MB02ND.f"
    }

/*     Step 3: Back transformation phase. */
/*             ------------------------- */
/*     Apply in backward order the Householder transformations (stored */
/*     in DWORK(IHOUSH)) performed onto the rows of C during the */
/*     bidiagonalization phase, to the selected base vectors (specified */
/*     by INUL(I) = .TRUE.). Already transformed vectors are those for */
/*     which BWORK(I) = .TRUE.. */

#line 646 "MB02ND.f"
    kf = k;
#line 647 "MB02ND.f"
    if (sufwrk && lfirst) {

/*        Enough workspace for a fast algorithm and first pass. */

#line 651 "MB02ND.f"
	ij = jv;

#line 653 "MB02ND.f"
	i__1 = k;
#line 653 "MB02ND.f"
	for (j = 1; j <= i__1; ++j) {
#line 654 "MB02ND.f"
	    dcopy_(&nl, &c__[iwork[j] * c_dim1 + 1], &c__1, &dwork[ij], &c__1)
		    ;
#line 655 "MB02ND.f"
	    ij += nl;
#line 656 "MB02ND.f"
/* L140: */
#line 656 "MB02ND.f"
	}

/*        Workspace: need   LDW + (N+L)*K + K, */
/*                   prefer LDW + (N+L)*K + K*NB. */

#line 661 "MB02ND.f"
	ij = jv;
#line 662 "MB02ND.f"
	jwork = ij + nl * k;
#line 663 "MB02ND.f"
	i__1 = *ldwork - jwork + 1;
#line 663 "MB02ND.f"
	dormbr_("P vectors", "Left", "No transpose", &nl, &k, &mnl, &dwork[
		ihoush], &p, &dwork[itaup], &dwork[ij], &nl, &dwork[jwork], &
		i__1, &ifail, (ftnlen)9, (ftnlen)4, (ftnlen)12);
/* Computing MAX */
#line 666 "MB02ND.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 666 "MB02ND.f"
	wrkopt = max(i__1,i__2);

#line 668 "MB02ND.f"
	i__1 = nl;
#line 668 "MB02ND.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 669 "MB02ND.f"
	    if (inul[i__] && ! bwork[i__]) {
#line 669 "MB02ND.f"
		bwork[i__] = TRUE_;
#line 669 "MB02ND.f"
	    }
#line 671 "MB02ND.f"
/* L160: */
#line 671 "MB02ND.f"
	}

#line 673 "MB02ND.f"
    } else {

/*        Not enough workspace for a fast algorithm or subsequent passes. */

#line 677 "MB02ND.f"
	i__1 = nl;
#line 677 "MB02ND.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 678 "MB02ND.f"
	    if (inul[i__] && ! bwork[i__]) {
#line 679 "MB02ND.f"
		kj = jv;

#line 681 "MB02ND.f"
		for (nj = min(p,mc); nj >= 1; --nj) {
#line 682 "MB02ND.f"
		    j = mc - nj + 1;
#line 683 "MB02ND.f"
		    kj -= j;
#line 684 "MB02ND.f"
		    first = dwork[kj];
#line 685 "MB02ND.f"
		    dwork[kj] = 1.;
#line 686 "MB02ND.f"
		    dlarf_("Left", &j, &c__1, &dwork[kj], &c__1, &dwork[itaup 
			    + nj - 1], &c__[nj + ioff + i__ * c_dim1], ldc, &
			    dwork[jwork], (ftnlen)4);
#line 689 "MB02ND.f"
		    dwork[kj] = first;
#line 690 "MB02ND.f"
/* L170: */
#line 690 "MB02ND.f"
		}

#line 692 "MB02ND.f"
		bwork[i__] = TRUE_;
#line 693 "MB02ND.f"
	    }
#line 694 "MB02ND.f"
/* L180: */
#line 694 "MB02ND.f"
	}
#line 695 "MB02ND.f"
    }

#line 697 "MB02ND.f"
    if (*rank <= 0) {
#line 697 "MB02ND.f"
	*rank = 0;
#line 697 "MB02ND.f"
    }
#line 699 "MB02ND.f"
    if (min(*rank,*l) == 0) {
#line 700 "MB02ND.f"
	if (sufwrk && lfirst) {
#line 700 "MB02ND.f"
	    dlacpy_("Full", &nl, &k, &dwork[jv], &nl, &c__[c_offset], ldc, (
		    ftnlen)4);
#line 700 "MB02ND.f"
	}
#line 702 "MB02ND.f"
	dwork[1] = (doublereal) wrkopt;
#line 703 "MB02ND.f"
	dwork[2] = 1.;
#line 704 "MB02ND.f"
	return 0;
#line 705 "MB02ND.f"
    }

/*     Step 4: Compute matrices F and Y */
/*             ------------------------ */
/*             using Householder transformation Q. */

/*     Compute the orthogonal matrix Q (in factorized form) and the */
/*     matrices F and Y using RQ factorization. It is assumed that, */
/*     generically, the last L rows of V2 matrix have full rank. */
/*     The code could not be the most efficient when RANK has been */
/*     lowered, because the already created zero pattern of the last */
/*     L rows of V2 matrix is not exploited. */

#line 718 "MB02ND.f"
    if (sufwrk && lfirst) {

/*        Enough workspace for a fast algorithm and first pass. */
/*        Workspace: need   LDW1 + 2*L, */
/*                   prefer LDW1 + L + L*NB, where */
/*                          LDW1 = LDW + (N+L)*K; */

#line 725 "MB02ND.f"
	itauq = jwork;
#line 726 "MB02ND.f"
	jwork = itauq + *l;
#line 727 "MB02ND.f"
	i__1 = *ldwork - jwork + 1;
#line 727 "MB02ND.f"
	dgerqf_(l, &k, &dwork[jv + *n], &nl, &dwork[itauq], &dwork[jwork], &
		i__1, info);
/* Computing MAX */
#line 729 "MB02ND.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 729 "MB02ND.f"
	wrkopt = max(i__1,i__2);

/*        Workspace: need   LDW1 + N+L, */
/*                   prefer LDW1 + L + N*NB. */

#line 734 "MB02ND.f"
	i__1 = *ldwork - jwork + 1;
#line 734 "MB02ND.f"
	dormrq_("Right", "Transpose", n, &k, l, &dwork[jv + *n], &nl, &dwork[
		itauq], &dwork[jv], &nl, &dwork[jwork], &i__1, info, (ftnlen)
		5, (ftnlen)9);
/* Computing MAX */
#line 737 "MB02ND.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 737 "MB02ND.f"
	wrkopt = max(i__1,i__2);

#line 739 "MB02ND.f"
	jf = jv + nl * (k - *l) + *n;
#line 740 "MB02ND.f"
	ldf = nl;
#line 741 "MB02ND.f"
	jwork = jf + ldf * *l - *n;
#line 742 "MB02ND.f"
	i__1 = k - *l;
#line 742 "MB02ND.f"
	dlaset_("Full", l, &i__1, &c_b4, &c_b4, &dwork[jv + *n], &ldf, (
		ftnlen)4);
#line 743 "MB02ND.f"
	if (*l > 1) {
#line 743 "MB02ND.f"
	    i__1 = *l - 1;
#line 743 "MB02ND.f"
	    i__2 = *l - 1;
#line 743 "MB02ND.f"
	    dlaset_("Lower", &i__1, &i__2, &c_b4, &c_b4, &dwork[jf + 1], &ldf,
		     (ftnlen)5);
#line 743 "MB02ND.f"
	}
#line 746 "MB02ND.f"
	ij = jv;

#line 748 "MB02ND.f"
	i__1 = k;
#line 748 "MB02ND.f"
	for (j = 1; j <= i__1; ++j) {
#line 749 "MB02ND.f"
	    dcopy_(&nl, &dwork[ij], &c__1, &c__[iwork[j] * c_dim1 + 1], &c__1)
		    ;
#line 750 "MB02ND.f"
	    ij += nl;
#line 751 "MB02ND.f"
/* L200: */
#line 751 "MB02ND.f"
	}

#line 753 "MB02ND.f"
    } else {

/*        Not enough workspace for a fast algorithm or subsequent passes. */
/*        Workspace: LDW2 + N+L, where LDW2 = LDW + L*L. */

#line 758 "MB02ND.f"
	i__ = nl;
#line 759 "MB02ND.f"
	jf = jv;
#line 760 "MB02ND.f"
	ldf = *l;
#line 761 "MB02ND.f"
	jwork = jf + ldf * *l;
/* Computing MAX */
#line 762 "MB02ND.f"
	i__1 = wrkopt, i__2 = jwork + nl - 1;
#line 762 "MB02ND.f"
	wrkopt = max(i__1,i__2);

/*        WHILE ( ( K >= 1 ) .AND. ( I > N ) ) DO */
#line 765 "MB02ND.f"
L220:
#line 766 "MB02ND.f"
	if (k >= 1 && i__ > *n) {

#line 768 "MB02ND.f"
	    i__1 = k;
#line 768 "MB02ND.f"
	    for (j = 1; j <= i__1; ++j) {
#line 769 "MB02ND.f"
		dwork[jwork + j - 1] = c__[i__ + iwork[j] * c_dim1];
#line 770 "MB02ND.f"
/* L240: */
#line 770 "MB02ND.f"
	    }

/*           Compute Householder transformation. */

#line 774 "MB02ND.f"
	    dlarfg_(&k, &dwork[jwork + k - 1], &dwork[jwork], &c__1, &temp);
#line 775 "MB02ND.f"
	    c__[i__ + iwork[k] * c_dim1] = dwork[jwork + k - 1];
#line 776 "MB02ND.f"
	    if (temp != 0.) {

/*              Apply Householder transformation onto the selected base */
/*              vectors. */

#line 781 "MB02ND.f"
		i__1 = i__ - 1;
#line 781 "MB02ND.f"
		for (i1 = 1; i1 <= i__1; ++i1) {
#line 782 "MB02ND.f"
		    inprod = c__[i1 + iwork[k] * c_dim1];

#line 784 "MB02ND.f"
		    i__2 = k - 1;
#line 784 "MB02ND.f"
		    for (j = 1; j <= i__2; ++j) {
#line 785 "MB02ND.f"
			inprod += dwork[jwork + j - 1] * c__[i1 + iwork[j] * 
				c_dim1];
#line 786 "MB02ND.f"
/* L260: */
#line 786 "MB02ND.f"
		    }

#line 788 "MB02ND.f"
		    hh = inprod * temp;
#line 789 "MB02ND.f"
		    c__[i1 + iwork[k] * c_dim1] -= hh;

#line 791 "MB02ND.f"
		    i__2 = k - 1;
#line 791 "MB02ND.f"
		    for (j = 1; j <= i__2; ++j) {
#line 792 "MB02ND.f"
			j1 = iwork[j];
#line 793 "MB02ND.f"
			c__[i1 + j1 * c_dim1] -= dwork[jwork + j - 1] * hh;
#line 794 "MB02ND.f"
			c__[i__ + j1 * c_dim1] = 0.;
#line 795 "MB02ND.f"
/* L280: */
#line 795 "MB02ND.f"
		    }

#line 797 "MB02ND.f"
/* L300: */
#line 797 "MB02ND.f"
		}

#line 799 "MB02ND.f"
	    }
#line 800 "MB02ND.f"
	    i__1 = i__ - *n;
#line 800 "MB02ND.f"
	    dcopy_(&i__1, &c__[n1 + iwork[k] * c_dim1], &c__1, &dwork[jf + (
		    i__ - *n - 1) * *l], &c__1);
#line 801 "MB02ND.f"
	    --k;
#line 802 "MB02ND.f"
	    --i__;
#line 803 "MB02ND.f"
	    goto L220;
#line 804 "MB02ND.f"
	}
/*        END WHILE 220 */
#line 806 "MB02ND.f"
    }

/*     Estimate the reciprocal condition number of the matrix F. */
/*     If F singular, lower the rank of the TLS approximation. */
/*     Workspace: LDW1 + 3*L or */
/*                LDW2 + 3*L. */

#line 813 "MB02ND.f"
    dtrcon_("1-norm", "Upper", "Non-unit", l, &dwork[jf], &ldf, &rcond, &
	    dwork[jwork], &iwork[kf + 1], info, (ftnlen)6, (ftnlen)5, (ftnlen)
	    8);
/* Computing MAX */
#line 815 "MB02ND.f"
    i__1 = wrkopt, i__2 = jwork + *l * 3 - 1;
#line 815 "MB02ND.f"
    wrkopt = max(i__1,i__2);

#line 817 "MB02ND.f"
    i__1 = *l;
#line 817 "MB02ND.f"
    for (j = 1; j <= i__1; ++j) {
#line 818 "MB02ND.f"
	dcopy_(n, &c__[iwork[kf - *l + j] * c_dim1 + 1], &c__1, &x[j * x_dim1 
		+ 1], &c__1);
#line 819 "MB02ND.f"
/* L320: */
#line 819 "MB02ND.f"
    }

#line 821 "MB02ND.f"
    fnorm = dlantr_("1-norm", "Upper", "Non-unit", l, l, &dwork[jf], &ldf, &
	    dwork[jwork], (ftnlen)6, (ftnlen)5, (ftnlen)8);
#line 823 "MB02ND.f"
    if (rcond <= eps * fnorm) {
#line 824 "MB02ND.f"
	--(*rank);
#line 825 "MB02ND.f"
	goto L340;
#line 826 "MB02ND.f"
    }
#line 827 "MB02ND.f"
    if (fnorm <= eps * dlange_("1-norm", n, l, &x[x_offset], ldx, &dwork[
	    jwork], (ftnlen)6)) {
#line 829 "MB02ND.f"
	*rank -= *l;
#line 830 "MB02ND.f"
	goto L340;
#line 831 "MB02ND.f"
    } else {
#line 832 "MB02ND.f"
	goto L400;
#line 833 "MB02ND.f"
    }

#line 835 "MB02ND.f"
L340:
#line 836 "MB02ND.f"
    iwarm = 2;
#line 837 "MB02ND.f"
    *theta = -1.;
#line 838 "MB02ND.f"
    if (sufwrk && lfirst) {

/*           Rearrange the stored Householder transformations for */
/*           subsequent passes, taking care to avoid overwriting. */

#line 843 "MB02ND.f"
	if (p < nl) {
#line 844 "MB02ND.f"
	    kj = ihoush + nl * (nl - 1);
#line 845 "MB02ND.f"
	    mj = ihoush + p * (nl - 1);

#line 847 "MB02ND.f"
	    i__1 = nl;
#line 847 "MB02ND.f"
	    for (nj = 1; nj <= i__1; ++nj) {
#line 848 "MB02ND.f"
		for (j = p - 1; j >= 0; --j) {
#line 849 "MB02ND.f"
		    dwork[kj + j] = dwork[mj + j];
#line 850 "MB02ND.f"
/* L350: */
#line 850 "MB02ND.f"
		}
#line 851 "MB02ND.f"
		kj -= nl;
#line 852 "MB02ND.f"
		mj -= p;
#line 853 "MB02ND.f"
/* L360: */
#line 853 "MB02ND.f"
	    }

#line 855 "MB02ND.f"
	}
#line 856 "MB02ND.f"
	kj = ihoush;
#line 857 "MB02ND.f"
	mj = ihoush + nl * ioff;

#line 859 "MB02ND.f"
	i__1 = min(p,mc);
#line 859 "MB02ND.f"
	for (nj = 1; nj <= i__1; ++nj) {
#line 860 "MB02ND.f"
	    i__2 = mc - nj;
#line 860 "MB02ND.f"
	    for (j = 0; j <= i__2; ++j) {
#line 861 "MB02ND.f"
		dwork[kj] = dwork[mj + j * p];
#line 862 "MB02ND.f"
		++kj;
#line 863 "MB02ND.f"
/* L370: */
#line 863 "MB02ND.f"
	    }
#line 864 "MB02ND.f"
	    mj = mj + nl + 1;
#line 865 "MB02ND.f"
/* L380: */
#line 865 "MB02ND.f"
	}

#line 867 "MB02ND.f"
	jv = kj;
#line 868 "MB02ND.f"
	lfirst = FALSE_;
#line 869 "MB02ND.f"
    }
#line 870 "MB02ND.f"
    goto L100;
/*     UNTIL ( F nonsingular, i.e., RCOND.GT.EPS*FNORM or */
/*                                  FNORM.GT.EPS*norm(Y) ) */
#line 873 "MB02ND.f"
L400:

/*     Step 5: Compute TLS solution. */
/*             -------------------- */
/*     Solve X F = -Y  by forward elimination  (F is upper triangular). */

#line 879 "MB02ND.f"
    dtrsm_("Right", "Upper", "No transpose", "Non-unit", n, l, &c_b85, &dwork[
	    jf], &ldf, &x[x_offset], ldx, (ftnlen)5, (ftnlen)5, (ftnlen)12, (
	    ftnlen)8);

/*     Set the optimal workspace and reciprocal condition number of F. */

#line 884 "MB02ND.f"
    dwork[1] = (doublereal) wrkopt;
#line 885 "MB02ND.f"
    dwork[2] = rcond;

#line 887 "MB02ND.f"
    return 0;
/* *** Last line of MB02ND *** */
} /* mb02nd_ */

