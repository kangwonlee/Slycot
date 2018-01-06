#line 1 "MB02MD.f"
/* MB02MD.f -- translated by f2c (version 20100827).
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

#line 1 "MB02MD.f"
/* Table of constant values */

static doublereal c_b8 = 0.;
static doublereal c_b9 = 1.;
static integer c__1 = 1;
static doublereal c_b52 = -1.;

/* Subroutine */ int mb02md_(char *job, integer *m, integer *n, integer *l, 
	integer *rank, doublereal *c__, integer *ldc, doublereal *s, 
	doublereal *x, integer *ldx, doublereal *tol, integer *iwork, 
	doublereal *dwork, integer *ldwork, integer *iwarn, integer *info, 
	ftnlen job_len)
{
    /* System generated locals */
    integer c_dim1, c_offset, x_dim1, x_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, k, p, n1, r1, nl, ldw;
    static logical ctol;
    static integer itau;
    static doublereal smax;
    static logical crank;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical ljobn;
    static doublereal rcond;
    static logical ljobr, ljobt;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal fnorm;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer jwork;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgerqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dgesvd_(char *, char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    extern doublereal dlantr_(char *, char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    extern /* Subroutine */ int dtrcon_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer minmnl;
    extern /* Subroutine */ int dormrq_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static doublereal toltmp;
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

/*     To solve the Total Least Squares (TLS) problem using a Singular */
/*     Value Decomposition (SVD) approach. */
/*     The TLS problem assumes an overdetermined set of linear equations */
/*     AX = B, where both the data matrix A as well as the observation */
/*     matrix B are inaccurate. The routine also solves determined and */
/*     underdetermined sets of equations by computing the minimum norm */
/*     solution. */
/*     It is assumed that all preprocessing measures (scaling, coordinate */
/*     transformations, whitening, ... ) of the data have been performed */
/*     in advance. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Determines whether the values of the parameters RANK and */
/*             TOL are to be specified by the user or computed by the */
/*             routine as follows: */
/*             = 'R':  Compute RANK only; */
/*             = 'T':  Compute TOL only; */
/*             = 'B':  Compute both RANK and TOL; */
/*             = 'N':  Compute neither RANK nor TOL. */

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
/*             On entry, if JOB = 'T' or JOB = 'N', then RANK must */
/*             specify r, the rank of the TLS approximation [A+DA|B+DB]. */
/*             RANK <= min(M,N). */
/*             Otherwise, r is computed by the routine. */
/*             On exit, if JOB = 'R' or JOB = 'B', and INFO = 0, then */
/*             RANK contains the computed (effective) rank of the TLS */
/*             approximation [A+DA|B+DB]. */
/*             Otherwise, the user-supplied value of RANK may be */
/*             changed by the routine on exit if the RANK-th and the */
/*             (RANK+1)-th singular values of C = [A|B] are considered */
/*             to be equal, or if the upper triangular matrix F (as */
/*             defined in METHOD) is (numerically) singular. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N+L) */
/*             On entry, the leading M-by-(N+L) part of this array must */
/*             contain the matrices A and B. Specifically, the first N */
/*             columns must contain the data matrix A and the last L */
/*             columns the observation matrix B (right-hand sides). */
/*             On exit, the leading (N+L)-by-(N+L) part of this array */
/*             contains the (transformed) right singular vectors, */
/*             including null space vectors, if any, of C = [A|B]. */
/*             Specifically, the leading (N+L)-by-RANK part of this array */
/*             always contains the first RANK right singular vectors, */
/*             corresponding to the largest singular values of C. If */
/*             L = 0, or if RANK = 0 and IWARN <> 2, the remaining */
/*             (N+L)-by-(N+L-RANK) top-right part of this array contains */
/*             the remaining N+L-RANK right singular vectors. Otherwise, */
/*             this part contains the matrix V2 transformed as described */
/*             in Step 3 of the TLS algorithm (see METHOD). */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= max(1,M,N+L). */

/*     S       (output) DOUBLE PRECISION array, dimension (min(M,N+L)) */
/*             If INFO = 0, the singular values of matrix C, ordered */
/*             such that S(1) >= S(2) >= ... >= S(p-1) >= S(p) >= 0, */
/*             where p = min(M,N+L). */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,L) */
/*             If INFO = 0, the leading N-by-L part of this array */
/*             contains the solution X to the TLS problem specified */
/*             by A and B. */

/*     LDX     INTEGER */
/*             The leading dimension of array X.  LDX >= max(1,N). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             A tolerance used to determine the rank of the TLS */
/*             approximation [A+DA|B+DB] and to check the multiplicity */
/*             of the singular values of matrix C. Specifically, S(i) */
/*             and S(j) (i < j) are considered to be equal if */
/*             SQRT(S(i)**2 - S(j)**2) <= TOL, and the TLS approximation */
/*             [A+DA|B+DB] has rank r if S(i) > TOL*S(1) (or S(i) > TOL, */
/*             if TOL specifies sdev (see below)), for i = 1,2,...,r. */
/*             TOL is also used to check the singularity of the upper */
/*             triangular matrix F (as defined in METHOD). */
/*             If JOB = 'R' or JOB = 'N', then TOL must specify the */
/*             desired tolerance. If the user sets TOL to be less than or */
/*             equal to 0, the tolerance is taken as EPS, where EPS is */
/*             the machine precision (see LAPACK Library routine DLAMCH). */
/*             Otherwise, the tolerance is computed by the routine and */
/*             the user must supply the non-negative value sdev, i.e. the */
/*             estimated standard deviation of the error on each element */
/*             of the matrix C, as input value of TOL. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (L) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK, and DWORK(2) returns the reciprocal of the */
/*             condition number of the matrix F. */
/*             If INFO > 0, DWORK(1:min(M,N+L)-1) contain the unconverged */
/*             non-diagonal elements of the bidiagonal matrix whose */
/*             diagonal is in S (see LAPACK Library routine DGESVD). */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK = max(2, 3*(N+L) + M, 5*(N+L)),       if M >= N+L; */
/*             LDWORK = max(2, M*(N+L) + max( 3M+N+L, 5*M), 3*L), */
/*                                                          if M <  N+L. */
/*             For optimum performance LDWORK should be larger. */

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
/*             > 0:  if the SVD algorithm (in LAPACK Library routine */
/*                   DBDSQR) has failed to converge. In this case, S(1), */
/*                   S(2), ..., S(INFO) may not have been found */
/*                   correctly and the remaining singular values may */
/*                   not be the smallest. This failure is not likely */
/*                   to occur. */

/*     METHOD */

/*     The method used is an extension (see [3,4,5]) of the classical */
/*     TLS algorithm proposed by Golub and Van Loan [1]. */

/*     Let [A|B] denote the matrix formed by adjoining the columns of B */
/*     to the columns of A on the right. */

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

/*     Define matrix C = [A|B] and s(i) as its i-th singular value for */
/*     i = 1,2,...,min(M,NL), where NL = N + L. If M < NL, then s(j) = 0 */
/*     for j = M+1,...,NL. */

/*     The Classical TLS algorithm proceeds as follows (see [3,4,5]): */

/*     Step 1: Compute part of the singular value decomposition (SVD) */
/*             USV' of C = [A|B], namely compute S and V'. (An initial */
/*             QR factorization of C is used when M is larger enough */
/*             than NL.) */

/*     Step 2: If not fixed by the user, compute the rank r0 of the data */
/*             [A|B] based on TOL as follows: if JOB = 'R' or JOB = 'N', */

/*                s(1) >= ... >= s(r0) > TOL*s(1) >= ... >= s(NL). */

/*             Otherwise, using [2], TOL can be computed from the */
/*             standard deviation sdev of the errors on [A|B]: */

/*                TOL = SQRT(2 * max(M,NL)) * sdev, */

/*             and the rank r0 is determined (if JOB = 'R' or 'B') using */

/*                s(1) >= ... >= s(r0) > TOL >= ... >= s(NL). */

/*             The rank r of the approximation [A+DA|B+DB] is then equal */
/*             to the minimum of N and r0. */

/*     Step 3: Let V2 be the matrix of the columns of V corresponding to */
/*             the (NL - r) smallest singular values of C, i.e. the last */
/*             (NL - r) columns of V. */
/*             Compute with Householder transformations the orthogonal */
/*             matrix Q such that: */

/*                       |VH   Y| */
/*              V2 x Q = |      | */
/*                       |0    F| */

/*             where VH is an N-by-(N - r) matrix, Y is an N-by-L matrix */
/*             and F is an L-by-L upper triangular matrix. */
/*             If F is singular, then lower the rank r with the */
/*             multiplicity of s(r) and repeat this step. */

/*     Step 4: If F is nonsingular then the solution X is obtained by */
/*             solving the following equations by forward elimination: */

/*                X F = -Y. */

/*     Notes : */
/*     The TLS solution is unique if r = N, F is nonsingular and */
/*     s(N) > s(N+1). */
/*     If F is singular, however, then the computed solution is infinite */
/*     and hence does not satisfy the second TLS criterion (see TLS */
/*     definition). For these cases, Golub and Van Loan [1] claim that */
/*     the TLS problem has no solution. The properties of these so-called */
/*     nongeneric problems are described in [4] and the TLS computations */
/*     are generalized in order to solve them. As proven in [4], the */
/*     proposed generalization satisfies the TLS criteria for any */
/*     number L of observation vectors in B provided that, in addition, */
/*     the solution | X| is constrained to be orthogonal to all vectors */
/*                  |-I| */
/*     of the form |w| which belong to the space generated by the columns */
/*                 |0| */
/*     of the submatrix |Y|. */
/*                      |F| */

/*     REFERENCES */

/*     [1] Golub, G.H. and Van Loan, C.F. */
/*         An Analysis of the Total Least-Squares Problem. */
/*         SIAM J. Numer. Anal., 17, pp. 883-893, 1980. */

/*     [2] Staar, J., Vandewalle, J. and Wemans, M. */
/*         Realization of Truncated Impulse Response Sequences with */
/*         Prescribed Uncertainty. */
/*         Proc. 8th IFAC World Congress, Kyoto, I, pp. 7-12, 1981. */

/*     [3] Van Huffel, S. */
/*         Analysis of the Total Least Squares Problem and its Use in */
/*         Parameter Estimation. */
/*         Doctoral dissertation, Dept. of Electr. Eng., Katholieke */
/*         Universiteit Leuven, Belgium, June 1987. */

/*     [4] Van Huffel, S. and Vandewalle, J. */
/*         Analysis and Solution of the Nongeneric Total Least Squares */
/*         Problem. */
/*         SIAM J. Matr. Anal. and Appl., 9, pp. 360-372, 1988. */

/*     [5] Van Huffel, S. and Vandewalle, J. */
/*         The Total Least Squares Problem: Computational Aspects and */
/*         Analysis. */
/*         Series "Frontiers in Applied Mathematics", Vol. 9, */
/*         SIAM, Philadelphia, 1991. */

/*     NUMERICAL ASPECTS */

/*     The algorithm consists in (backward) stable steps. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997. */
/*     Supersedes Release 2.0 routine MB02AD by S. Van Huffel, Katholieke */
/*     University, Leuven, Belgium. */

/*     REVISIONS */

/*     June 24, 1997, Feb. 27, 2000, Oct. 19, 2003, Feb. 21, 2004. */

/*     KEYWORDS */

/*     Least-squares approximation, singular subspace, singular value */
/*     decomposition, singular values, total least-squares. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 334 "MB02MD.f"
    /* Parameter adjustments */
#line 334 "MB02MD.f"
    c_dim1 = *ldc;
#line 334 "MB02MD.f"
    c_offset = 1 + c_dim1;
#line 334 "MB02MD.f"
    c__ -= c_offset;
#line 334 "MB02MD.f"
    --s;
#line 334 "MB02MD.f"
    x_dim1 = *ldx;
#line 334 "MB02MD.f"
    x_offset = 1 + x_dim1;
#line 334 "MB02MD.f"
    x -= x_offset;
#line 334 "MB02MD.f"
    --iwork;
#line 334 "MB02MD.f"
    --dwork;
#line 334 "MB02MD.f"

#line 334 "MB02MD.f"
    /* Function Body */
#line 334 "MB02MD.f"
    *iwarn = 0;
#line 335 "MB02MD.f"
    *info = 0;
#line 336 "MB02MD.f"
    nl = *n + *l;
#line 337 "MB02MD.f"
    k = max(*m,nl);
#line 338 "MB02MD.f"
    p = min(*m,*n);
#line 339 "MB02MD.f"
    minmnl = min(*m,nl);
/* Computing MAX */
#line 340 "MB02MD.f"
    i__1 = minmnl * 3 + k, i__2 = minmnl * 5;
#line 340 "MB02MD.f"
    ldw = max(i__1,i__2);
#line 341 "MB02MD.f"
    ljobr = lsame_(job, "R", (ftnlen)1, (ftnlen)1);
#line 342 "MB02MD.f"
    ljobt = lsame_(job, "T", (ftnlen)1, (ftnlen)1);
#line 343 "MB02MD.f"
    ljobn = lsame_(job, "N", (ftnlen)1, (ftnlen)1);

/*     Determine whether RANK or/and TOL is/are to be computed. */

#line 347 "MB02MD.f"
    crank = ! ljobt && ! ljobn;
#line 348 "MB02MD.f"
    ctol = ! ljobr && ! ljobn;

/*     Test the input scalar arguments. */

#line 352 "MB02MD.f"
    if (ctol && crank && ! lsame_(job, "B", (ftnlen)1, (ftnlen)1)) {
#line 353 "MB02MD.f"
	*info = -1;
#line 354 "MB02MD.f"
    } else if (*m < 0) {
#line 355 "MB02MD.f"
	*info = -2;
#line 356 "MB02MD.f"
    } else if (*n < 0) {
#line 357 "MB02MD.f"
	*info = -3;
#line 358 "MB02MD.f"
    } else if (*l < 0) {
#line 359 "MB02MD.f"
	*info = -4;
#line 360 "MB02MD.f"
    } else if (! crank && *rank > p) {
#line 361 "MB02MD.f"
	*info = -5;
#line 362 "MB02MD.f"
    } else if (*ldc < max(1,k)) {
#line 363 "MB02MD.f"
	*info = -7;
#line 364 "MB02MD.f"
    } else if (*ldx < max(1,*n)) {
#line 365 "MB02MD.f"
	*info = -10;
#line 366 "MB02MD.f"
    } else if (ctol && *tol < 0.) {
#line 367 "MB02MD.f"
	*info = -11;
#line 368 "MB02MD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 368 "MB02MD.f"
	i__1 = 2, i__2 = *m * nl + ldw, i__1 = max(i__1,i__2), i__2 = *l * 3;
#line 368 "MB02MD.f"
	if (*m >= nl && *ldwork < max(2,ldw) || *m < nl && *ldwork < max(i__1,
		i__2)) {
#line 371 "MB02MD.f"
	    *info = -14;
#line 372 "MB02MD.f"
	}
#line 372 "MB02MD.f"
    }

#line 374 "MB02MD.f"
    if (*info != 0) {

/*        Error return. */

#line 378 "MB02MD.f"
	i__1 = -(*info);
#line 378 "MB02MD.f"
	xerbla_("MB02MD", &i__1, (ftnlen)6);
#line 379 "MB02MD.f"
	return 0;
#line 380 "MB02MD.f"
    }

/*     Quick return if possible. */

#line 384 "MB02MD.f"
    if (crank) {
#line 384 "MB02MD.f"
	*rank = p;
#line 384 "MB02MD.f"
    }
#line 386 "MB02MD.f"
    if (min(*m,nl) == 0) {
#line 387 "MB02MD.f"
	if (*m == 0) {
#line 388 "MB02MD.f"
	    dlaset_("Full", &nl, &nl, &c_b8, &c_b9, &c__[c_offset], ldc, (
		    ftnlen)4);
#line 389 "MB02MD.f"
	    dlaset_("Full", n, l, &c_b8, &c_b8, &x[x_offset], ldx, (ftnlen)4);
#line 390 "MB02MD.f"
	}
#line 391 "MB02MD.f"
	dwork[1] = 2.;
#line 392 "MB02MD.f"
	dwork[2] = 1.;
#line 393 "MB02MD.f"
	return 0;
#line 394 "MB02MD.f"
    }

/*     Subroutine MB02MD solves a set of linear equations by a Total */
/*     Least Squares Approximation. */

/*     Step 1: Compute part of the singular value decomposition (SVD) */
/*             USV' of C = [A   |B   ], namely compute S and V'. */
/*                           M,N  M,L */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 409 "MB02MD.f"
    if (*m >= nl) {

/*        M >= N + L:  Overwrite V' on C. */
/*        Workspace: need max(3*min(M,N+L) + max(M,N+L), 5*min(M,N+L)). */

#line 414 "MB02MD.f"
	jwork = 1;
#line 415 "MB02MD.f"
	i__1 = *ldwork - jwork + 1;
#line 415 "MB02MD.f"
	dgesvd_("No left vectors", "Overwritten on C", m, &nl, &c__[c_offset],
		 ldc, &s[1], &dwork[1], &c__1, &dwork[1], &c__1, &dwork[jwork]
		, &i__1, info, (ftnlen)15, (ftnlen)16);
#line 418 "MB02MD.f"
    } else {

/*        M < N + L:  Save C in the workspace and compute V' in C. */
/*        Note that the previous DGESVD call cannot be used in this case. */
/*        Workspace: need M*(N+L) + max(3*min(M,N+L) + max(M,N+L), */
/*                                      5*min(M,N+L)). */

#line 425 "MB02MD.f"
	dlacpy_("Full", m, &nl, &c__[c_offset], ldc, &dwork[1], m, (ftnlen)4);
#line 426 "MB02MD.f"
	jwork = *m * nl + 1;
#line 427 "MB02MD.f"
	i__1 = *ldwork - jwork + 1;
#line 427 "MB02MD.f"
	dgesvd_("No left vectors", "All right vectors", m, &nl, &dwork[1], m, 
		&s[1], &dwork[1], &c__1, &c__[c_offset], ldc, &dwork[jwork], &
		i__1, info, (ftnlen)15, (ftnlen)17);
#line 430 "MB02MD.f"
    }

#line 432 "MB02MD.f"
    if (*info > 0) {

/*        Save the unconverged non-diagonal elements of the bidiagonal */
/*        matrix and exit. */

#line 437 "MB02MD.f"
	i__1 = minmnl - 1;
#line 437 "MB02MD.f"
	for (j = 1; j <= i__1; ++j) {
#line 438 "MB02MD.f"
	    dwork[j] = dwork[jwork + j];
#line 439 "MB02MD.f"
/* L10: */
#line 439 "MB02MD.f"
	}

#line 441 "MB02MD.f"
	return 0;
#line 442 "MB02MD.f"
    }
/* Computing MAX */
#line 443 "MB02MD.f"
    i__1 = 2, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 443 "MB02MD.f"
    wrkopt = max(i__1,i__2);

/*     Transpose V' in-situ (in C). */

#line 447 "MB02MD.f"
    i__1 = nl;
#line 447 "MB02MD.f"
    for (j = 2; j <= i__1; ++j) {
#line 448 "MB02MD.f"
	i__2 = j - 1;
#line 448 "MB02MD.f"
	dswap_(&i__2, &c__[j + c_dim1], ldc, &c__[j * c_dim1 + 1], &c__1);
#line 449 "MB02MD.f"
/* L20: */
#line 449 "MB02MD.f"
    }

/*     Step 2: Compute the rank of the approximation [A+DA|B+DB]. */

#line 453 "MB02MD.f"
    if (ctol) {
#line 454 "MB02MD.f"
	toltmp = sqrt((doublereal) k * 2.) * *tol;
#line 455 "MB02MD.f"
	smax = toltmp;
#line 456 "MB02MD.f"
    } else {
#line 457 "MB02MD.f"
	toltmp = *tol;
#line 458 "MB02MD.f"
	if (toltmp <= 0.) {
#line 458 "MB02MD.f"
	    toltmp = dlamch_("Precision", (ftnlen)9);
#line 458 "MB02MD.f"
	}
/* Computing MAX */
#line 459 "MB02MD.f"
	d__1 = toltmp * s[1], d__2 = dlamch_("Safe minimum", (ftnlen)12);
#line 459 "MB02MD.f"
	smax = max(d__1,d__2);
#line 460 "MB02MD.f"
    }

#line 462 "MB02MD.f"
    if (crank) {
/*        WHILE ( RANK .GT. 0 ) .AND. ( S(RANK) .LE. SMAX ) DO */
#line 464 "MB02MD.f"
L40:
#line 464 "MB02MD.f"
	if (*rank > 0) {
#line 465 "MB02MD.f"
	    if (s[*rank] <= smax) {
#line 466 "MB02MD.f"
		--(*rank);
#line 467 "MB02MD.f"
		goto L40;
#line 468 "MB02MD.f"
	    }
#line 469 "MB02MD.f"
	}
/*        END WHILE 40 */
#line 471 "MB02MD.f"
    }

#line 473 "MB02MD.f"
    if (*l == 0) {
#line 474 "MB02MD.f"
	dwork[1] = (doublereal) wrkopt;
#line 475 "MB02MD.f"
	dwork[2] = 1.;
#line 476 "MB02MD.f"
	return 0;
#line 477 "MB02MD.f"
    }

#line 479 "MB02MD.f"
    n1 = *n + 1;
#line 480 "MB02MD.f"
    itau = 1;
#line 481 "MB02MD.f"
    jwork = itau + *l;

/*     Step 3: Compute the orthogonal matrix Q and matrices F and Y */
/*     such that F is nonsingular. */

/*     REPEAT */

/*        Adjust the rank if S(RANK) has multiplicity greater than 1. */

#line 490 "MB02MD.f"
L60:
#line 491 "MB02MD.f"
    r1 = *rank + 1;
#line 492 "MB02MD.f"
    if (*rank < minmnl) {
/*           WHILE RANK.GT.0 .AND. S(RANK)**2 - S(R1)**2.LE.TOL**2 DO */
#line 494 "MB02MD.f"
L80:
#line 494 "MB02MD.f"
	if (*rank > 0) {
/* Computing 2nd power */
#line 495 "MB02MD.f"
	    d__1 = s[r1] / s[*rank];
/* Computing 2nd power */
#line 495 "MB02MD.f"
	    d__2 = toltmp / s[*rank];
#line 495 "MB02MD.f"
	    if (1. - d__1 * d__1 <= d__2 * d__2) {
#line 497 "MB02MD.f"
		--(*rank);
#line 498 "MB02MD.f"
		*iwarn = 1;
#line 499 "MB02MD.f"
		goto L80;
#line 500 "MB02MD.f"
	    }
#line 501 "MB02MD.f"
	}
/*           END WHILE 80 */
#line 503 "MB02MD.f"
    }

#line 505 "MB02MD.f"
    if (*rank == 0) {

/*           Return zero solution. */

#line 509 "MB02MD.f"
	dlaset_("Full", n, l, &c_b8, &c_b8, &x[x_offset], ldx, (ftnlen)4);
#line 510 "MB02MD.f"
	dwork[1] = (doublereal) wrkopt;
#line 511 "MB02MD.f"
	dwork[2] = 1.;
#line 512 "MB02MD.f"
	return 0;
#line 513 "MB02MD.f"
    }

/*        Compute the orthogonal matrix Q (in factorized form) and the */
/*        matrices F and Y using RQ factorization. It is assumed that, */
/*        generically, the last L rows of V2 matrix have full rank. */
/*        The code could not be the most efficient one when RANK has been */
/*        lowered, because the already created zero pattern of the last */
/*        L rows of V2 matrix is not exploited. */
/*        Workspace: need 2*L;  prefer L + L*NB. */

#line 523 "MB02MD.f"
    r1 = *rank + 1;
#line 524 "MB02MD.f"
    i__1 = nl - *rank;
#line 524 "MB02MD.f"
    i__2 = *ldwork - jwork + 1;
#line 524 "MB02MD.f"
    dgerqf_(l, &i__1, &c__[n1 + r1 * c_dim1], ldc, &dwork[itau], &dwork[jwork]
	    , &i__2, info);
/* Computing MAX */
#line 526 "MB02MD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 526 "MB02MD.f"
    wrkopt = max(i__1,i__2);

/*        Workspace: need N+L;  prefer L + N*NB. */

#line 530 "MB02MD.f"
    i__1 = nl - *rank;
#line 530 "MB02MD.f"
    i__2 = *ldwork - jwork + 1;
#line 530 "MB02MD.f"
    dormrq_("Right", "Transpose", n, &i__1, l, &c__[n1 + r1 * c_dim1], ldc, &
	    dwork[itau], &c__[r1 * c_dim1 + 1], ldc, &dwork[jwork], &i__2, 
	    info, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
#line 533 "MB02MD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 533 "MB02MD.f"
    wrkopt = max(i__1,i__2);

#line 535 "MB02MD.f"
    i__1 = *n - *rank;
#line 535 "MB02MD.f"
    dlaset_("Full", l, &i__1, &c_b8, &c_b8, &c__[n1 + r1 * c_dim1], ldc, (
	    ftnlen)4);
#line 536 "MB02MD.f"
    if (*l > 1) {
#line 536 "MB02MD.f"
	i__1 = *l - 1;
#line 536 "MB02MD.f"
	i__2 = *l - 1;
#line 536 "MB02MD.f"
	dlaset_("Lower", &i__1, &i__2, &c_b8, &c_b8, &c__[n1 + 1 + n1 * 
		c_dim1], ldc, (ftnlen)5);
#line 536 "MB02MD.f"
    }

/*        Estimate the reciprocal condition number of the matrix F, */
/*        and lower the rank if F can be considered as singular. */
/*        Workspace: need 3*L. */

#line 544 "MB02MD.f"
    dtrcon_("1-norm", "Upper", "Non-unit", l, &c__[n1 + n1 * c_dim1], ldc, &
	    rcond, &dwork[1], &iwork[1], info, (ftnlen)6, (ftnlen)5, (ftnlen)
	    8);
/* Computing MAX */
#line 546 "MB02MD.f"
    i__1 = wrkopt, i__2 = *l * 3;
#line 546 "MB02MD.f"
    wrkopt = max(i__1,i__2);

#line 548 "MB02MD.f"
    fnorm = dlantr_("1-norm", "Upper", "Non-unit", l, l, &c__[n1 + n1 * 
	    c_dim1], ldc, &dwork[1], (ftnlen)6, (ftnlen)5, (ftnlen)8);
#line 550 "MB02MD.f"
    if (rcond <= toltmp * fnorm) {
#line 551 "MB02MD.f"
	--(*rank);
#line 552 "MB02MD.f"
	*iwarn = 2;
#line 553 "MB02MD.f"
	goto L60;
#line 554 "MB02MD.f"
    } else if (fnorm <= toltmp * dlange_("1-norm", n, l, &c__[n1 * c_dim1 + 1]
	    , ldc, &dwork[1], (ftnlen)6)) {
#line 556 "MB02MD.f"
	*rank -= *l;
#line 557 "MB02MD.f"
	*iwarn = 2;
#line 558 "MB02MD.f"
	goto L60;
#line 559 "MB02MD.f"
    }
/*     UNTIL ( F nonsingular, i.e., RCOND.GT.TOL*FNORM or */
/*                                  FNORM.GT.TOL*norm(Y) ) */

/*     Step 4: Solve X F = -Y by forward elimination, */
/*             (F is upper triangular). */

#line 566 "MB02MD.f"
    dlacpy_("Full", n, l, &c__[n1 * c_dim1 + 1], ldc, &x[x_offset], ldx, (
	    ftnlen)4);
#line 567 "MB02MD.f"
    dtrsm_("Right", "Upper", "No transpose", "Non-unit", n, l, &c_b52, &c__[
	    n1 + n1 * c_dim1], ldc, &x[x_offset], ldx, (ftnlen)5, (ftnlen)5, (
	    ftnlen)12, (ftnlen)8);

/*     Set the optimal workspace and reciprocal condition number of F. */

#line 572 "MB02MD.f"
    dwork[1] = (doublereal) wrkopt;
#line 573 "MB02MD.f"
    dwork[2] = rcond;

#line 575 "MB02MD.f"
    return 0;
/* *** Last line of MB02MD *** */
} /* mb02md_ */

