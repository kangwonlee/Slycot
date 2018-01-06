#line 1 "SB04OD.f"
/* SB04OD.f -- translated by f2c (version 20100827).
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

#line 1 "SB04OD.f"
/* Table of constant values */

static integer c__0 = 0;
static doublereal c_b52 = 1.;
static doublereal c_b53 = 0.;
static integer c__1 = 1;

/* Subroutine */ int sb04od_(char *reduce, char *trans, char *jobd, integer *
	m, integer *n, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *e, integer *lde, doublereal *f, integer *ldf, doublereal *
	scale, doublereal *dif, doublereal *p, integer *ldp, doublereal *q, 
	integer *ldq, doublereal *u, integer *ldu, doublereal *v, integer *
	ldv, integer *iwork, doublereal *dwork, integer *ldwork, integer *
	info, ftnlen reduce_len, ftnlen trans_len, ftnlen jobd_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, e_dim1, e_offset, f_dim1, f_offset, p_dim1, p_offset, 
	    q_dim1, q_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, mn, ijob;
    static doublereal anrm, bnrm, dnrm;
    static integer ierr;
    static doublereal enrm;
    static logical ljob1, ljob2;
    extern /* Subroutine */ int dgegs_(char *, char *, integer *, doublereal *
	    , integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen);
    static logical ljobd, ljobf;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), dlabad_(
	    doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    static logical ljobdf;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static logical ilascl, ilbscl, ildscl, ilescl, lredra, lredrb, lredua, 
	    lredub, lreduc;
    static doublereal bignum, safmax, safmin;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical lredur, ltrann;
    static doublereal anrmto, bnrmto, dnrmto, enrmto;
    extern /* Subroutine */ int dtgsyl_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, ftnlen);
    static integer minwrk;
    static doublereal smlnum;
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

/*     To solve for R and L one of the generalized Sylvester equations */

/*        A * R - L * B = scale * C ) */
/*                                  )                                 (1) */
/*        D * R - L * E = scale * F ) */

/*     or */

/*        A' * R + D' * L = scale * C    ) */
/*                                       )                            (2) */
/*        R * B' + L * E' = scale * (-F) ) */

/*     where A and D are M-by-M matrices, B and E are N-by-N matrices and */
/*     C, F, R and L are M-by-N matrices. */

/*     The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an */
/*     output scaling factor chosen to avoid overflow. */

/*     The routine also optionally computes a Dif estimate, which */
/*     measures the separation of the spectrum of the matrix pair (A,D) */
/*     from the spectrum of the matrix pair (B,E), Dif[(A,D),(B,E)]. */

/*     ARGUMENTS */

/*     MODE PARAMETERS */

/*     REDUCE  CHARACTER*1 */
/*             Indicates whether the matrix pairs (A,D) and/or (B,E) are */
/*             to be reduced to generalized Schur form as follows: */
/*             = 'R':  The matrix pairs (A,D) and (B,E) are to be reduced */
/*                     to generalized (real) Schur canonical form; */
/*             = 'A':  The matrix pair (A,D) only is to be reduced */
/*                     to generalized (real) Schur canonical form, */
/*                     and the matrix pair (B,E) already is in this form; */
/*             = 'B':  The matrix pair (B,E) only is to be reduced */
/*                     to generalized (real) Schur canonical form, */
/*                     and the matrix pair (A,D) already is in this form; */
/*             = 'N':  The matrix pairs (A,D) and (B,E) are already in */
/*                     generalized (real) Schur canonical form, as */
/*                     produced by LAPACK routine DGEES. */

/*     TRANS   CHARACTER*1 */
/*             Indicates which of the equations, (1) or (2), is to be */
/*             solved as follows: */
/*             = 'N':  The generalized Sylvester equation (1) is to be */
/*                     solved; */
/*             = 'T':  The "transposed" generalized Sylvester equation */
/*                     (2) is to be solved. */

/*     JOBD    CHARACTER*1 */
/*             Indicates whether the Dif estimator is to be computed as */
/*             follows: */
/*             = '1':  Only the one-norm-based Dif estimate is computed */
/*                     and stored in DIF; */
/*             = '2':  Only the Frobenius norm-based Dif estimate is */
/*                     computed and stored in DIF; */
/*             = 'D':  The equation (1) is solved and the one-norm-based */
/*                     Dif estimate is computed and stored in DIF; */
/*             = 'F':  The equation (1) is solved and the Frobenius norm- */
/*                     based Dif estimate is computed and stored in DIF; */
/*             = 'N':  The Dif estimator is not required and hence DIF is */
/*                     not referenced. (Solve either (1) or (2) only.) */
/*             JOBD is not referenced if TRANS = 'T'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of the matrices A and D and the number of rows */
/*             of the matrices C, F, R and L.  M >= 0. */

/*     N       (input) INTEGER */
/*             The order of the matrices B and E and the number of */
/*             columns of the matrices C, F, R and L.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,M) */
/*             On entry, the leading M-by-M part of this array must */
/*             contain the coefficient matrix A of the equation; A must */
/*             be in upper quasi-triangular form if REDUCE = 'B' or 'N'. */
/*             On exit, the leading M-by-M part of this array contains */
/*             the upper quasi-triangular form of A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the coefficient matrix B of the equation; B must */
/*             be in upper quasi-triangular form if REDUCE = 'A' or 'N'. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the upper quasi-triangular form of B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the right-hand side matrix C of the first equation */
/*             in (1) or (2). */
/*             On exit, if JOBD = 'N', 'D' or 'F', the leading M-by-N */
/*             part of this array contains the solution matrix R of the */
/*             problem; if JOBD = '1' or '2' and TRANS = 'N', the leading */
/*             M-by-N part of this array contains the solution matrix R */
/*             achieved during the computation of the Dif estimate. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,M). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading M-by-M part of this array must */
/*             contain the coefficient matrix D of the equation; D must */
/*             be in upper triangular form if REDUCE = 'B' or 'N'. */
/*             On exit, the leading M-by-M part of this array contains */
/*             the upper triangular form of D. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,M). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the coefficient matrix E of the equation; E must */
/*             be in upper triangular form if REDUCE = 'A' or 'N'. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the upper triangular form of E. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,N). */

/*     F       (input/output) DOUBLE PRECISION array, dimension (LDF,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the right-hand side matrix F of the second */
/*             equation in (1) or (2). */
/*             On exit, if JOBD = 'N', 'D' or 'F', the leading M-by-N */
/*             part of this array contains the solution matrix L of the */
/*             problem; if JOBD = '1' or '2' and TRANS = 'N', the leading */
/*             M-by-N part of this array contains the solution matrix L */
/*             achieved during the computation of the Dif estimate. */

/*     LDF     INTEGER */
/*             The leading dimension of array F.  LDF >= MAX(1,M). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scaling factor in (1) or (2). If 0 < SCALE < 1, C and */
/*             F hold the solutions R and L, respectively, to a slightly */
/*             perturbed system (but the input or computed generalized */
/*             (real) Schur canonical form matrices A, B, D, and E */
/*             have not been changed). If SCALE = 0, C and F hold the */
/*             solutions R and L, respectively, to the homogeneous system */
/*             with C = F = 0. Normally, SCALE = 1. */

/*     DIF     (output) DOUBLE PRECISION */
/*             If TRANS = 'N' and JOBD <> 'N', then DIF contains the */
/*             value of the Dif estimator, which is an upper bound of */
/*                                                    -1 */
/*             Dif[(A,D),(B,E)] = sigma_min(Z) = 1/||Z  ||, in either the */
/*             one-norm, or Frobenius norm, respectively (see METHOD). */
/*             Otherwise, DIF is not referenced. */

/*     P       (output) DOUBLE PRECISION array, dimension (LDP,*) */
/*             If REDUCE = 'R' or 'A', then the leading M-by-M part of */
/*             this array contains the (left) transformation matrix used */
/*             to reduce (A,D) to generalized Schur form. */
/*             Otherwise, P is not referenced and can be supplied as a */
/*             dummy array (i.e. set parameter LDP = 1 and declare this */
/*             array to be P(1,1) in the calling program). */

/*     LDP     INTEGER */
/*             The leading dimension of array P. */
/*             LDP >= MAX(1,M) if REDUCE = 'R' or 'A', */
/*             LDP >= 1        if REDUCE = 'B' or 'N'. */

/*     Q       (output) DOUBLE PRECISION array, dimension (LDQ,*) */
/*             If REDUCE = 'R' or 'A', then the leading M-by-M part of */
/*             this array contains the (right) transformation matrix used */
/*             to reduce (A,D) to generalized Schur form. */
/*             Otherwise, Q is not referenced and can be supplied as a */
/*             dummy array (i.e. set parameter LDQ = 1 and declare this */
/*             array to be Q(1,1) in the calling program). */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. */
/*             LDQ >= MAX(1,M) if REDUCE = 'R' or 'A', */
/*             LDQ >= 1        if REDUCE = 'B' or 'N'. */

/*     U       (output) DOUBLE PRECISION array, dimension (LDU,*) */
/*             If REDUCE = 'R' or 'B', then the leading N-by-N part of */
/*             this array contains the (left) transformation matrix used */
/*             to reduce (B,E) to generalized Schur form. */
/*             Otherwise, U is not referenced and can be supplied as a */
/*             dummy array (i.e. set parameter LDU = 1 and declare this */
/*             array to be U(1,1) in the calling program). */

/*     LDU     INTEGER */
/*             The leading dimension of array U. */
/*             LDU >= MAX(1,N) if REDUCE = 'R' or 'B', */
/*             LDU >= 1        if REDUCE = 'A' or 'N'. */

/*     V       (output) DOUBLE PRECISION array, dimension (LDV,*) */
/*             If REDUCE = 'R' or 'B', then the leading N-by-N part of */
/*             this array contains the (right) transformation matrix used */
/*             to reduce (B,E) to generalized Schur form. */
/*             Otherwise, V is not referenced and can be supplied as a */
/*             dummy array (i.e. set parameter LDV = 1 and declare this */
/*             array to be V(1,1) in the calling program). */

/*     LDV     INTEGER */
/*             The leading dimension of array V. */
/*             LDV >= MAX(1,N) if REDUCE = 'R' or 'B', */
/*             LDV >= 1        if REDUCE = 'A' or 'N'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M+N+6) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             If TRANS = 'N' and JOBD = 'D' or 'F', then */
/*                LDWORK = MAX(1,7*M,7*N,2*M*N) if REDUCE = 'R'; */
/*                LDWORK = MAX(1,7*M,2*M*N)     if REDUCE = 'A'; */
/*                LDWORK = MAX(1,7*N,2*M*N)     if REDUCE = 'B'; */
/*                LDWORK = MAX(1,2*M*N)         if REDUCE = 'N'. */
/*             Otherwise, the term 2*M*N above should be omitted. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if REDUCE <> 'N' and either (A,D) and/or (B,E) */
/*                   cannot be reduced to generalized Schur form; */
/*             = 2:  if REDUCE = 'N' and either A or B is not in */
/*                   upper quasi-triangular form; */
/*             = 3:  if a singular matrix was encountered during the */
/*                   computation of the solution matrices R and L, that */
/*                   is (A,D) and (B,E) have common or close eigenvalues. */

/*     METHOD */

/*     For the case TRANS = 'N', and REDUCE = 'R' or 'N', the algorithm */
/*     used by the routine consists of four steps (see [1] and [2]) as */
/*     follows: */

/*        (a) if REDUCE = 'R', then the matrix pairs (A,D) and (B,E) are */
/*            transformed to generalized Schur form, i.e. orthogonal */
/*            matrices P, Q, U and V are computed such that P' * A * Q */
/*            and U' * B * V are in upper quasi-triangular form and */
/*            P' * D * Q and U' * E * V are in upper triangular form; */
/*        (b) if REDUCE = 'R', then the matrices C and F are transformed */
/*            to give P' * C * V and P' * F * V respectively; */
/*        (c) if REDUCE = 'R', then the transformed system */

/*            P' * A * Q * R1 - L1 * U' * B * V = scale * P' * C * V */
/*            P' * D * Q * R1 - L1 * U' * E * V = scale * P' * F * V */

/*            is solved to give R1 and L1; otherwise, equation (1) is */
/*            solved to give R and L directly. The Dif estimator */
/*            is also computed if JOBD <> 'N'. */
/*        (d) if REDUCE = 'R', then the solution is transformed back */
/*            to give R = Q * R1 * V' and L = P * L1 * U'. */

/*     By using Kronecker products, equation (1) can also be written as */
/*     the system of linear equations Z * x = scale*y (see [1]), where */

/*            | I*A    I*D  | */
/*        Z = |             |. */
/*            |-B'*I  -E'*I | */

/*                                              -1 */
/*     If JOBD <> 'N', then a lower bound on ||Z  ||, in either the one- */
/*     norm or Frobenius norm, is computed, which in most cases is */
/*     a reliable estimate of the true value. Notice that since Z is a */
/*     matrix of order 2 * M * N, the exact value of Dif (i.e., in the */
/*     Frobenius norm case, the smallest singular value of Z) may be very */
/*     expensive to compute. */

/*     The case TRANS = 'N', and REDUCE = 'A' or 'B', is similar, but */
/*     only one of the matrix pairs should be reduced and the */
/*     calculations simplify. */

/*     For the case TRANS = 'T', and REDUCE = 'R' or 'N', the algorithm */
/*     is similar, but the steps (b), (c), and (d) are as follows: */

/*        (b) if REDUCE = 'R', then the matrices C and F are transformed */
/*            to give Q' * C * V and P' * F * U respectively; */
/*        (c) if REDUCE = 'R', then the transformed system */

/*            Q' * A' * P * R1 + Q' * D' * P * L1 =  scale * Q' * C * V */
/*            R1 * V' * B' * U + L1 * V' * E' * U = -scale * P' * F * U */

/*            is solved to give R1 and L1; otherwise, equation (2) is */
/*            solved to give R and L directly. */
/*        (d) if REDUCE = 'R', then the solution is transformed back */
/*            to give R = P * R1 * V' and L = P * L1 * V'. */

/*     REFERENCES */

/*     [1] Kagstrom, B. and Westin, L. */
/*         Generalized Schur Methods with Condition Estimators for */
/*         Solving the Generalized Sylvester Equation. */
/*         IEEE Trans. Auto. Contr., 34, pp. 745-751, 1989. */
/*     [2] Kagstrom, B. and Westin, L. */
/*         GSYLV - Fortran Routines for the Generalized Schur Method with */
/*         Dif Estimators for Solving the Generalized Sylvester */
/*         Equation. */
/*         Report UMINF-132.86, Institute of Information Processing, */
/*         Univ. of Umea, Sweden, July 1987. */
/*     [3] Golub, G.H., Nash, S. and Van Loan, C.F. */
/*         A Hessenberg-Schur Method for the Problem AX + XB = C. */
/*         IEEE Trans. Auto. Contr., AC-24, pp. 909-913, 1979. */
/*     [4] Kagstrom, B. and Van Dooren, P. */
/*         Additive Decomposition of a Transfer Function with respect to */
/*         a Specified Region. */
/*         In: "Signal Processing, Scattering and Operator Theory, and */
/*         Numerical Methods" (Eds. M.A. Kaashoek et al.). */
/*         Proceedings of MTNS-89, Vol. 3, pp. 469-477, Birkhauser Boston */
/*         Inc., 1990. */
/*     [5] Kagstrom, B. and Van Dooren, P. */
/*         A Generalized State-space Approach for the Additive */
/*         Decomposition of a Transfer Matrix. */
/*         Report UMINF-91.12, Institute of Information Processing, Univ. */
/*         of Umea, Sweden, April 1991. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. A reliable estimate for the */
/*     condition number of Z in the Frobenius norm, is (see [1]) */

/*        K(Z) = SQRT(  ||A||**2 + ||B||**2 + ||C||**2 + ||D||**2 )/DIF. */

/*     If mu is an upper bound on the relative error of the elements of */
/*     the matrices A, B, C, D, E and F, then the relative error in the */
/*     actual solution is approximately mu * K(Z). */

/*     The relative error in the computed solution (due to rounding */
/*     errors) is approximately EPS * K(Z), where EPS is the machine */
/*     precision (see LAPACK Library routine DLAMCH). */

/*     FURTHER COMMENTS */

/*     For applications of the generalized Sylvester equation in control */
/*     theory, see [4] and [5]. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB04CD by Bo Kagstrom and Lars */
/*     Westin. */

/*     REVISIONS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, May 1999, Dec. 1999, */
/*     May 2009. */

/*     KEYWORDS */

/*     Generalized eigenvalue problem, orthogonal transformation, real */
/*     Schur form, Sylvester equation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 422 "SB04OD.f"
    /* Parameter adjustments */
#line 422 "SB04OD.f"
    a_dim1 = *lda;
#line 422 "SB04OD.f"
    a_offset = 1 + a_dim1;
#line 422 "SB04OD.f"
    a -= a_offset;
#line 422 "SB04OD.f"
    b_dim1 = *ldb;
#line 422 "SB04OD.f"
    b_offset = 1 + b_dim1;
#line 422 "SB04OD.f"
    b -= b_offset;
#line 422 "SB04OD.f"
    c_dim1 = *ldc;
#line 422 "SB04OD.f"
    c_offset = 1 + c_dim1;
#line 422 "SB04OD.f"
    c__ -= c_offset;
#line 422 "SB04OD.f"
    d_dim1 = *ldd;
#line 422 "SB04OD.f"
    d_offset = 1 + d_dim1;
#line 422 "SB04OD.f"
    d__ -= d_offset;
#line 422 "SB04OD.f"
    e_dim1 = *lde;
#line 422 "SB04OD.f"
    e_offset = 1 + e_dim1;
#line 422 "SB04OD.f"
    e -= e_offset;
#line 422 "SB04OD.f"
    f_dim1 = *ldf;
#line 422 "SB04OD.f"
    f_offset = 1 + f_dim1;
#line 422 "SB04OD.f"
    f -= f_offset;
#line 422 "SB04OD.f"
    p_dim1 = *ldp;
#line 422 "SB04OD.f"
    p_offset = 1 + p_dim1;
#line 422 "SB04OD.f"
    p -= p_offset;
#line 422 "SB04OD.f"
    q_dim1 = *ldq;
#line 422 "SB04OD.f"
    q_offset = 1 + q_dim1;
#line 422 "SB04OD.f"
    q -= q_offset;
#line 422 "SB04OD.f"
    u_dim1 = *ldu;
#line 422 "SB04OD.f"
    u_offset = 1 + u_dim1;
#line 422 "SB04OD.f"
    u -= u_offset;
#line 422 "SB04OD.f"
    v_dim1 = *ldv;
#line 422 "SB04OD.f"
    v_offset = 1 + v_dim1;
#line 422 "SB04OD.f"
    v -= v_offset;
#line 422 "SB04OD.f"
    --iwork;
#line 422 "SB04OD.f"
    --dwork;
#line 422 "SB04OD.f"

#line 422 "SB04OD.f"
    /* Function Body */
#line 422 "SB04OD.f"
    *info = 0;
#line 423 "SB04OD.f"
    mn = max(*m,*n);
#line 424 "SB04OD.f"
    lredur = lsame_(reduce, "R", (ftnlen)1, (ftnlen)1);
#line 425 "SB04OD.f"
    lredua = lsame_(reduce, "A", (ftnlen)1, (ftnlen)1);
#line 426 "SB04OD.f"
    lredub = lsame_(reduce, "B", (ftnlen)1, (ftnlen)1);
#line 427 "SB04OD.f"
    lredra = lredur || lredua;
#line 428 "SB04OD.f"
    lredrb = lredur || lredub;
#line 429 "SB04OD.f"
    lreduc = lredra || lredub;
#line 430 "SB04OD.f"
    if (lredur) {
/* Computing MAX */
#line 431 "SB04OD.f"
	i__1 = 1, i__2 = mn * 7;
#line 431 "SB04OD.f"
	minwrk = max(i__1,i__2);
#line 432 "SB04OD.f"
    } else if (lredua) {
/* Computing MAX */
#line 433 "SB04OD.f"
	i__1 = 1, i__2 = *m * 7;
#line 433 "SB04OD.f"
	minwrk = max(i__1,i__2);
#line 434 "SB04OD.f"
    } else if (lredub) {
/* Computing MAX */
#line 435 "SB04OD.f"
	i__1 = 1, i__2 = *n * 7;
#line 435 "SB04OD.f"
	minwrk = max(i__1,i__2);
#line 436 "SB04OD.f"
    } else {
#line 437 "SB04OD.f"
	minwrk = 1;
#line 438 "SB04OD.f"
    }
#line 439 "SB04OD.f"
    ltrann = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 440 "SB04OD.f"
    if (ltrann) {
#line 441 "SB04OD.f"
	ljob1 = lsame_(jobd, "1", (ftnlen)1, (ftnlen)1);
#line 442 "SB04OD.f"
	ljob2 = lsame_(jobd, "2", (ftnlen)1, (ftnlen)1);
#line 443 "SB04OD.f"
	ljobd = lsame_(jobd, "D", (ftnlen)1, (ftnlen)1);
#line 444 "SB04OD.f"
	ljobf = lsame_(jobd, "F", (ftnlen)1, (ftnlen)1);
#line 445 "SB04OD.f"
	ljobdf = ljob1 || ljob2 || ljobd || ljobf;
#line 446 "SB04OD.f"
	if (ljobd || ljobf) {
/* Computing MAX */
#line 446 "SB04OD.f"
	    i__1 = minwrk, i__2 = (*m << 1) * *n;
#line 446 "SB04OD.f"
	    minwrk = max(i__1,i__2);
#line 446 "SB04OD.f"
	}
#line 447 "SB04OD.f"
    }

/*     Test the input scalar arguments. */

#line 451 "SB04OD.f"
    if (! lreduc && ! lsame_(reduce, "N", (ftnlen)1, (ftnlen)1)) {
#line 452 "SB04OD.f"
	*info = -1;
#line 453 "SB04OD.f"
    } else if (! ltrann && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
#line 454 "SB04OD.f"
	*info = -2;
#line 455 "SB04OD.f"
    } else if (ltrann) {
#line 456 "SB04OD.f"
	if (! ljobdf && ! lsame_(jobd, "N", (ftnlen)1, (ftnlen)1)) {
#line 456 "SB04OD.f"
	    *info = -3;
#line 456 "SB04OD.f"
	}
#line 458 "SB04OD.f"
    }
#line 459 "SB04OD.f"
    if (*m < 0) {
#line 460 "SB04OD.f"
	*info = -4;
#line 461 "SB04OD.f"
    } else if (*n < 0) {
#line 462 "SB04OD.f"
	*info = -5;
#line 463 "SB04OD.f"
    } else if (*lda < max(1,*m)) {
#line 464 "SB04OD.f"
	*info = -7;
#line 465 "SB04OD.f"
    } else if (*ldb < max(1,*n)) {
#line 466 "SB04OD.f"
	*info = -9;
#line 467 "SB04OD.f"
    } else if (*ldc < max(1,*m)) {
#line 468 "SB04OD.f"
	*info = -11;
#line 469 "SB04OD.f"
    } else if (*ldd < max(1,*m)) {
#line 470 "SB04OD.f"
	*info = -13;
#line 471 "SB04OD.f"
    } else if (*lde < max(1,*n)) {
#line 472 "SB04OD.f"
	*info = -15;
#line 473 "SB04OD.f"
    } else if (*ldf < max(1,*m)) {
#line 474 "SB04OD.f"
	*info = -17;
#line 475 "SB04OD.f"
    } else if (! lredra && *ldp < 1 || lredra && *ldp < max(1,*m)) {
#line 477 "SB04OD.f"
	*info = -21;
#line 478 "SB04OD.f"
    } else if (! lredra && *ldq < 1 || lredra && *ldq < max(1,*m)) {
#line 480 "SB04OD.f"
	*info = -23;
#line 481 "SB04OD.f"
    } else if (! lredrb && *ldu < 1 || lredrb && *ldu < max(1,*n)) {
#line 483 "SB04OD.f"
	*info = -25;
#line 484 "SB04OD.f"
    } else if (! lredrb && *ldv < 1 || lredrb && *ldv < max(1,*n)) {
#line 486 "SB04OD.f"
	*info = -27;
#line 487 "SB04OD.f"
    } else if (*ldwork < minwrk) {
#line 488 "SB04OD.f"
	*info = -30;
#line 489 "SB04OD.f"
    }

#line 491 "SB04OD.f"
    if (*info != 0) {

/*        Error return. */

#line 495 "SB04OD.f"
	i__1 = -(*info);
#line 495 "SB04OD.f"
	xerbla_("SB04OD", &i__1, (ftnlen)6);
#line 496 "SB04OD.f"
	return 0;
#line 497 "SB04OD.f"
    }

/*     Quick return if possible. */

#line 501 "SB04OD.f"
    if (*n == 0 || *m == 0) {
#line 502 "SB04OD.f"
	*scale = 1.;
#line 503 "SB04OD.f"
	dwork[1] = 1.;
#line 504 "SB04OD.f"
	if (ltrann) {
#line 505 "SB04OD.f"
	    if (ljobdf) {
#line 505 "SB04OD.f"
		*dif = 1.;
#line 505 "SB04OD.f"
	    }
#line 506 "SB04OD.f"
	}
#line 507 "SB04OD.f"
	return 0;
#line 508 "SB04OD.f"
    }
#line 509 "SB04OD.f"
    wrkopt = 1;
#line 510 "SB04OD.f"
    sufwrk = *ldwork >= *m * *n;

/*     STEP 1: Reduce (A,D) and/or (B,E) to generalized Schur form. */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 520 "SB04OD.f"
    if (lreduc) {

/*        Get machine constants. */

#line 524 "SB04OD.f"
	safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 525 "SB04OD.f"
	safmax = 1. / safmin;
#line 526 "SB04OD.f"
	dlabad_(&safmin, &safmax);
#line 527 "SB04OD.f"
	smlnum = sqrt(safmin) / dlamch_("Precision", (ftnlen)9);
#line 528 "SB04OD.f"
	bignum = 1. / smlnum;

#line 530 "SB04OD.f"
	if (! lredub) {

/*           Scale A if max element outside range [SMLNUM,BIGNUM]. */

#line 534 "SB04OD.f"
	    anrm = dlange_("M", m, m, &a[a_offset], lda, &dwork[1], (ftnlen)1)
		    ;
#line 535 "SB04OD.f"
	    ilascl = FALSE_;
#line 536 "SB04OD.f"
	    if (anrm > 0. && anrm < smlnum) {
#line 537 "SB04OD.f"
		anrmto = smlnum;
#line 538 "SB04OD.f"
		ilascl = TRUE_;
#line 539 "SB04OD.f"
	    } else if (anrm > bignum) {
#line 540 "SB04OD.f"
		anrmto = bignum;
#line 541 "SB04OD.f"
		ilascl = TRUE_;
#line 542 "SB04OD.f"
	    }
#line 543 "SB04OD.f"
	    if (ilascl) {
#line 543 "SB04OD.f"
		dlascl_("G", &c__0, &c__0, &anrm, &anrmto, m, m, &a[a_offset],
			 lda, &ierr, (ftnlen)1);
#line 543 "SB04OD.f"
	    }

/*           Scale D if max element outside range [SMLNUM,BIGNUM] */

#line 549 "SB04OD.f"
	    dnrm = dlange_("M", m, m, &d__[d_offset], ldd, &dwork[1], (ftnlen)
		    1);
#line 550 "SB04OD.f"
	    ildscl = FALSE_;
#line 551 "SB04OD.f"
	    if (dnrm > 0. && dnrm < smlnum) {
#line 552 "SB04OD.f"
		dnrmto = smlnum;
#line 553 "SB04OD.f"
		ildscl = TRUE_;
#line 554 "SB04OD.f"
	    } else if (dnrm > bignum) {
#line 555 "SB04OD.f"
		dnrmto = bignum;
#line 556 "SB04OD.f"
		ildscl = TRUE_;
#line 557 "SB04OD.f"
	    }
#line 558 "SB04OD.f"
	    if (ildscl) {
#line 558 "SB04OD.f"
		dlascl_("G", &c__0, &c__0, &dnrm, &dnrmto, m, m, &d__[
			d_offset], ldd, &ierr, (ftnlen)1);
#line 558 "SB04OD.f"
	    }

/*           Reduce (A,D) to generalized Schur form. */
/*           Workspace:  need   7*M; */
/*                       prefer 5*M + M*(NB+1). */

#line 566 "SB04OD.f"
	    i__1 = *ldwork - *m * 3;
#line 566 "SB04OD.f"
	    dgegs_("Vectors left", "Vectors right", m, &a[a_offset], lda, &
		    d__[d_offset], ldd, &dwork[1], &dwork[*m + 1], &dwork[(*m 
		    << 1) + 1], &p[p_offset], ldp, &q[q_offset], ldq, &dwork[*
		    m * 3 + 1], &i__1, info, (ftnlen)12, (ftnlen)13);

/*           Undo scaling */

#line 572 "SB04OD.f"
	    if (ilascl) {
#line 572 "SB04OD.f"
		dlascl_("H", &c__0, &c__0, &anrmto, &anrm, m, m, &a[a_offset],
			 lda, &ierr, (ftnlen)1);
#line 572 "SB04OD.f"
	    }

#line 576 "SB04OD.f"
	    if (ildscl) {
#line 576 "SB04OD.f"
		dlascl_("U", &c__0, &c__0, &dnrmto, &dnrm, m, m, &d__[
			d_offset], ldd, &ierr, (ftnlen)1);
#line 576 "SB04OD.f"
	    }

#line 580 "SB04OD.f"
	    if (*info != 0) {
#line 581 "SB04OD.f"
		*info = 1;
#line 582 "SB04OD.f"
		return 0;
#line 583 "SB04OD.f"
	    }
/* Computing MAX */
#line 584 "SB04OD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[*m * 3 + 1] + *m * 3;
#line 584 "SB04OD.f"
	    wrkopt = max(i__1,i__2);
#line 585 "SB04OD.f"
	}
#line 586 "SB04OD.f"
	if (! lredua) {

/*           Scale B if max element outside range [SMLNUM,BIGNUM] */

#line 590 "SB04OD.f"
	    bnrm = dlange_("M", n, n, &b[b_offset], ldb, &dwork[1], (ftnlen)1)
		    ;
#line 591 "SB04OD.f"
	    ilbscl = FALSE_;
#line 592 "SB04OD.f"
	    if (bnrm > 0. && bnrm < smlnum) {
#line 593 "SB04OD.f"
		bnrmto = smlnum;
#line 594 "SB04OD.f"
		ilbscl = TRUE_;
#line 595 "SB04OD.f"
	    } else if (bnrm > bignum) {
#line 596 "SB04OD.f"
		bnrmto = bignum;
#line 597 "SB04OD.f"
		ilbscl = TRUE_;
#line 598 "SB04OD.f"
	    }
#line 599 "SB04OD.f"
	    if (ilbscl) {
#line 599 "SB04OD.f"
		dlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset],
			 ldb, &ierr, (ftnlen)1);
#line 599 "SB04OD.f"
	    }

/*           Scale E if max element outside range [SMLNUM,BIGNUM] */

#line 605 "SB04OD.f"
	    enrm = dlange_("M", n, n, &e[e_offset], lde, &dwork[1], (ftnlen)1)
		    ;
#line 606 "SB04OD.f"
	    ilescl = FALSE_;
#line 607 "SB04OD.f"
	    if (enrm > 0. && enrm < smlnum) {
#line 608 "SB04OD.f"
		enrmto = smlnum;
#line 609 "SB04OD.f"
		ilescl = TRUE_;
#line 610 "SB04OD.f"
	    } else if (enrm > bignum) {
#line 611 "SB04OD.f"
		enrmto = bignum;
#line 612 "SB04OD.f"
		ilescl = TRUE_;
#line 613 "SB04OD.f"
	    }
#line 614 "SB04OD.f"
	    if (ilescl) {
#line 614 "SB04OD.f"
		dlascl_("G", &c__0, &c__0, &enrm, &enrmto, n, n, &e[e_offset],
			 lde, &ierr, (ftnlen)1);
#line 614 "SB04OD.f"
	    }

/*           Reduce (B,E) to generalized Schur form. */
/*           Workspace:  need   7*N; */
/*                       prefer 5*N + N*(NB+1). */

#line 622 "SB04OD.f"
	    i__1 = *ldwork - *n * 3;
#line 622 "SB04OD.f"
	    dgegs_("Vectors left", "Vectors right", n, &b[b_offset], ldb, &e[
		    e_offset], lde, &dwork[1], &dwork[*n + 1], &dwork[(*n << 
		    1) + 1], &u[u_offset], ldu, &v[v_offset], ldv, &dwork[*n *
		     3 + 1], &i__1, info, (ftnlen)12, (ftnlen)13);

/*           Undo scaling */

#line 628 "SB04OD.f"
	    if (ilbscl) {
#line 628 "SB04OD.f"
		dlascl_("H", &c__0, &c__0, &bnrmto, &bnrm, n, n, &b[b_offset],
			 ldb, &ierr, (ftnlen)1);
#line 628 "SB04OD.f"
	    }

#line 632 "SB04OD.f"
	    if (ilescl) {
#line 632 "SB04OD.f"
		dlascl_("U", &c__0, &c__0, &enrmto, &enrm, n, n, &e[e_offset],
			 lde, &ierr, (ftnlen)1);
#line 632 "SB04OD.f"
	    }

#line 636 "SB04OD.f"
	    if (*info != 0) {
#line 637 "SB04OD.f"
		*info = 1;
#line 638 "SB04OD.f"
		return 0;
#line 639 "SB04OD.f"
	    }
/* Computing MAX */
#line 640 "SB04OD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[*n * 3 + 1] + *n * 3;
#line 640 "SB04OD.f"
	    wrkopt = max(i__1,i__2);
#line 641 "SB04OD.f"
	}
#line 642 "SB04OD.f"
    }

#line 644 "SB04OD.f"
    if (! lredur) {

/*        Set INFO = 2 if A and/or B are/is not in quasi-triangular form. */

#line 648 "SB04OD.f"
	if (! lredua) {
#line 649 "SB04OD.f"
	    i__ = 1;

#line 651 "SB04OD.f"
L20:
#line 652 "SB04OD.f"
	    if (i__ <= *m - 2) {
#line 653 "SB04OD.f"
		if (a[i__ + 1 + i__ * a_dim1] != 0.) {
#line 654 "SB04OD.f"
		    if (a[i__ + 2 + (i__ + 1) * a_dim1] != 0.) {
#line 655 "SB04OD.f"
			*info = 2;
#line 656 "SB04OD.f"
			return 0;
#line 657 "SB04OD.f"
		    } else {
#line 658 "SB04OD.f"
			++i__;
#line 659 "SB04OD.f"
		    }
#line 660 "SB04OD.f"
		}
#line 661 "SB04OD.f"
		++i__;
#line 662 "SB04OD.f"
		goto L20;
#line 663 "SB04OD.f"
	    }
#line 664 "SB04OD.f"
	}

#line 666 "SB04OD.f"
	if (! lredub) {
#line 667 "SB04OD.f"
	    i__ = 1;

#line 669 "SB04OD.f"
L40:
#line 670 "SB04OD.f"
	    if (i__ <= *n - 2) {
#line 671 "SB04OD.f"
		if (b[i__ + 1 + i__ * b_dim1] != 0.) {
#line 672 "SB04OD.f"
		    if (b[i__ + 2 + (i__ + 1) * b_dim1] != 0.) {
#line 673 "SB04OD.f"
			*info = 2;
#line 674 "SB04OD.f"
			return 0;
#line 675 "SB04OD.f"
		    } else {
#line 676 "SB04OD.f"
			++i__;
#line 677 "SB04OD.f"
		    }
#line 678 "SB04OD.f"
		}
#line 679 "SB04OD.f"
		++i__;
#line 680 "SB04OD.f"
		goto L40;
#line 681 "SB04OD.f"
	    }
#line 682 "SB04OD.f"
	}
#line 683 "SB04OD.f"
    }

/*     STEP 2: Modify right hand sides (C,F). */

#line 687 "SB04OD.f"
    if (lreduc) {
/* Computing MAX */
#line 688 "SB04OD.f"
	i__1 = wrkopt, i__2 = *m * *n;
#line 688 "SB04OD.f"
	wrkopt = max(i__1,i__2);
#line 689 "SB04OD.f"
	if (sufwrk) {

/*           Enough workspace for a BLAS 3 calculation. */

#line 693 "SB04OD.f"
	    if (ltrann) {

/*              Equation (1). */

#line 697 "SB04OD.f"
		if (! lredub) {
#line 698 "SB04OD.f"
		    dgemm_("Transpose", "No transpose", m, n, m, &c_b52, &p[
			    p_offset], ldp, &c__[c_offset], ldc, &c_b53, &
			    dwork[1], m, (ftnlen)9, (ftnlen)12);
#line 700 "SB04OD.f"
		} else {
#line 701 "SB04OD.f"
		    dlacpy_("Full", m, n, &c__[c_offset], ldc, &dwork[1], m, (
			    ftnlen)4);
#line 702 "SB04OD.f"
		}
#line 703 "SB04OD.f"
		if (! lredua) {
#line 704 "SB04OD.f"
		    dgemm_("No transpose", "No transpose", m, n, n, &c_b52, &
			    dwork[1], m, &v[v_offset], ldv, &c_b53, &c__[
			    c_offset], ldc, (ftnlen)12, (ftnlen)12);
#line 706 "SB04OD.f"
		} else {
#line 707 "SB04OD.f"
		    dlacpy_("Full", m, n, &dwork[1], m, &c__[c_offset], ldc, (
			    ftnlen)4);
#line 708 "SB04OD.f"
		}
#line 709 "SB04OD.f"
		if (! lredub) {
#line 710 "SB04OD.f"
		    dgemm_("Transpose", "No transpose", m, n, m, &c_b52, &p[
			    p_offset], ldp, &f[f_offset], ldf, &c_b53, &dwork[
			    1], m, (ftnlen)9, (ftnlen)12);
#line 712 "SB04OD.f"
		} else {
#line 713 "SB04OD.f"
		    dlacpy_("Full", m, n, &f[f_offset], ldf, &dwork[1], m, (
			    ftnlen)4);
#line 714 "SB04OD.f"
		}
#line 715 "SB04OD.f"
		if (! lredua) {
#line 716 "SB04OD.f"
		    dgemm_("No transpose", "No transpose", m, n, n, &c_b52, &
			    dwork[1], m, &v[v_offset], ldv, &c_b53, &f[
			    f_offset], ldf, (ftnlen)12, (ftnlen)12);
#line 718 "SB04OD.f"
		} else {
#line 719 "SB04OD.f"
		    dlacpy_("Full", m, n, &dwork[1], m, &f[f_offset], ldf, (
			    ftnlen)4);
#line 720 "SB04OD.f"
		}
#line 721 "SB04OD.f"
	    } else {

/*              Equation (2). */

#line 725 "SB04OD.f"
		if (! lredub) {
#line 726 "SB04OD.f"
		    dgemm_("Transpose", "No transpose", m, n, m, &c_b52, &q[
			    q_offset], ldq, &c__[c_offset], ldc, &c_b53, &
			    dwork[1], m, (ftnlen)9, (ftnlen)12);
#line 728 "SB04OD.f"
		} else {
#line 729 "SB04OD.f"
		    dlacpy_("Full", m, n, &c__[c_offset], ldc, &dwork[1], m, (
			    ftnlen)4);
#line 730 "SB04OD.f"
		}
#line 731 "SB04OD.f"
		if (! lredua) {
#line 732 "SB04OD.f"
		    dgemm_("No transpose", "No transpose", m, n, n, &c_b52, &
			    dwork[1], m, &v[v_offset], ldv, &c_b53, &c__[
			    c_offset], ldc, (ftnlen)12, (ftnlen)12);
#line 734 "SB04OD.f"
		} else {
#line 735 "SB04OD.f"
		    dlacpy_("Full", m, n, &dwork[1], m, &c__[c_offset], ldc, (
			    ftnlen)4);
#line 736 "SB04OD.f"
		}
#line 737 "SB04OD.f"
		if (! lredub) {
#line 738 "SB04OD.f"
		    dgemm_("Transpose", "No transpose", m, n, m, &c_b52, &p[
			    p_offset], ldp, &f[f_offset], ldf, &c_b53, &dwork[
			    1], m, (ftnlen)9, (ftnlen)12);
#line 740 "SB04OD.f"
		} else {
#line 741 "SB04OD.f"
		    dlacpy_("Full", m, n, &f[f_offset], ldf, &dwork[1], m, (
			    ftnlen)4);
#line 742 "SB04OD.f"
		}
#line 743 "SB04OD.f"
		if (! lredua) {
#line 744 "SB04OD.f"
		    dgemm_("No transpose", "No transpose", m, n, n, &c_b52, &
			    dwork[1], m, &u[u_offset], ldu, &c_b53, &f[
			    f_offset], ldf, (ftnlen)12, (ftnlen)12);
#line 746 "SB04OD.f"
		} else {
#line 747 "SB04OD.f"
		    dlacpy_("Full", m, n, &dwork[1], m, &f[f_offset], ldf, (
			    ftnlen)4);
#line 748 "SB04OD.f"
		}
#line 749 "SB04OD.f"
	    }
#line 750 "SB04OD.f"
	} else {

/*           Use a BLAS 2 calculation. */

#line 754 "SB04OD.f"
	    if (ltrann) {

/*              Equation (1). */

#line 758 "SB04OD.f"
		if (! lredub) {

#line 760 "SB04OD.f"
		    i__1 = *n;
#line 760 "SB04OD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 761 "SB04OD.f"
			dgemv_("Transpose", m, m, &c_b52, &p[p_offset], ldp, &
				c__[i__ * c_dim1 + 1], &c__1, &c_b53, &dwork[
				1], &c__1, (ftnlen)9);
#line 763 "SB04OD.f"
			dcopy_(m, &dwork[1], &c__1, &c__[i__ * c_dim1 + 1], &
				c__1);
#line 764 "SB04OD.f"
/* L60: */
#line 764 "SB04OD.f"
		    }

#line 766 "SB04OD.f"
		}
#line 767 "SB04OD.f"
		if (! lredua) {

#line 769 "SB04OD.f"
		    i__1 = *m;
#line 769 "SB04OD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 770 "SB04OD.f"
			dgemv_("Transpose", n, n, &c_b52, &v[v_offset], ldv, &
				c__[i__ + c_dim1], ldc, &c_b53, &dwork[1], &
				c__1, (ftnlen)9);
#line 772 "SB04OD.f"
			dcopy_(n, &dwork[1], &c__1, &c__[i__ + c_dim1], ldc);
#line 773 "SB04OD.f"
/* L80: */
#line 773 "SB04OD.f"
		    }

#line 775 "SB04OD.f"
		}
#line 776 "SB04OD.f"
		if (! lredub) {

#line 778 "SB04OD.f"
		    i__1 = *n;
#line 778 "SB04OD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 779 "SB04OD.f"
			dgemv_("Transpose", m, m, &c_b52, &p[p_offset], ldp, &
				f[i__ * f_dim1 + 1], &c__1, &c_b53, &dwork[1],
				 &c__1, (ftnlen)9);
#line 781 "SB04OD.f"
			dcopy_(m, &dwork[1], &c__1, &f[i__ * f_dim1 + 1], &
				c__1);
#line 782 "SB04OD.f"
/* L100: */
#line 782 "SB04OD.f"
		    }

#line 784 "SB04OD.f"
		}
#line 785 "SB04OD.f"
		if (! lredua) {

#line 787 "SB04OD.f"
		    i__1 = *m;
#line 787 "SB04OD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 788 "SB04OD.f"
			dgemv_("Transpose", n, n, &c_b52, &v[v_offset], ldv, &
				f[i__ + f_dim1], ldf, &c_b53, &dwork[1], &
				c__1, (ftnlen)9);
#line 790 "SB04OD.f"
			dcopy_(n, &dwork[1], &c__1, &f[i__ + f_dim1], ldf);
#line 791 "SB04OD.f"
/* L120: */
#line 791 "SB04OD.f"
		    }

#line 793 "SB04OD.f"
		}
#line 794 "SB04OD.f"
	    } else {

/*              Equation (2). */

#line 798 "SB04OD.f"
		if (! lredub) {

#line 800 "SB04OD.f"
		    i__1 = *n;
#line 800 "SB04OD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 801 "SB04OD.f"
			dgemv_("Transpose", m, m, &c_b52, &q[q_offset], ldq, &
				c__[i__ * c_dim1 + 1], &c__1, &c_b53, &dwork[
				1], &c__1, (ftnlen)9);
#line 803 "SB04OD.f"
			dcopy_(m, &dwork[1], &c__1, &c__[i__ * c_dim1 + 1], &
				c__1);
#line 804 "SB04OD.f"
/* L140: */
#line 804 "SB04OD.f"
		    }

#line 806 "SB04OD.f"
		}
#line 807 "SB04OD.f"
		if (! lredua) {

#line 809 "SB04OD.f"
		    i__1 = *m;
#line 809 "SB04OD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 810 "SB04OD.f"
			dgemv_("Transpose", n, n, &c_b52, &v[v_offset], ldv, &
				c__[i__ + c_dim1], ldc, &c_b53, &dwork[1], &
				c__1, (ftnlen)9);
#line 812 "SB04OD.f"
			dcopy_(n, &dwork[1], &c__1, &c__[i__ + c_dim1], ldc);
#line 813 "SB04OD.f"
/* L160: */
#line 813 "SB04OD.f"
		    }

#line 815 "SB04OD.f"
		}
#line 816 "SB04OD.f"
		if (! lredub) {

#line 818 "SB04OD.f"
		    i__1 = *n;
#line 818 "SB04OD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 819 "SB04OD.f"
			dgemv_("Transpose", m, m, &c_b52, &p[p_offset], ldp, &
				f[i__ * f_dim1 + 1], &c__1, &c_b53, &dwork[1],
				 &c__1, (ftnlen)9);
#line 821 "SB04OD.f"
			dcopy_(m, &dwork[1], &c__1, &f[i__ * f_dim1 + 1], &
				c__1);
#line 822 "SB04OD.f"
/* L180: */
#line 822 "SB04OD.f"
		    }

#line 824 "SB04OD.f"
		}
#line 825 "SB04OD.f"
		if (! lredua) {

#line 827 "SB04OD.f"
		    i__1 = *m;
#line 827 "SB04OD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 828 "SB04OD.f"
			dgemv_("Transpose", n, n, &c_b52, &u[u_offset], ldu, &
				f[i__ + f_dim1], ldf, &c_b53, &dwork[1], &
				c__1, (ftnlen)9);
#line 830 "SB04OD.f"
			dcopy_(n, &dwork[1], &c__1, &f[i__ + f_dim1], ldf);
#line 831 "SB04OD.f"
/* L200: */
#line 831 "SB04OD.f"
		    }

#line 833 "SB04OD.f"
		}
#line 834 "SB04OD.f"
	    }
#line 835 "SB04OD.f"
	}
#line 836 "SB04OD.f"
    }

/*     STEP 3: Solve the transformed system and compute the Dif */
/*             estimator. */

#line 841 "SB04OD.f"
    if (ltrann) {
#line 842 "SB04OD.f"
	if (ljobd) {
#line 843 "SB04OD.f"
	    ijob = 1;
#line 844 "SB04OD.f"
	} else if (ljobf) {
#line 845 "SB04OD.f"
	    ijob = 2;
#line 846 "SB04OD.f"
	} else if (ljob1) {
#line 847 "SB04OD.f"
	    ijob = 3;
#line 848 "SB04OD.f"
	} else if (ljob2) {
#line 849 "SB04OD.f"
	    ijob = 4;
#line 850 "SB04OD.f"
	} else {
#line 851 "SB04OD.f"
	    ijob = 0;
#line 852 "SB04OD.f"
	}
#line 853 "SB04OD.f"
    } else {
#line 854 "SB04OD.f"
	ijob = 0;
#line 855 "SB04OD.f"
    }

/*     Workspace:  need 2*M*N if TRANS = 'N' and JOBD = 'D' or 'F'; */
/*                      1, otherwise. */

#line 860 "SB04OD.f"
    dtgsyl_(trans, &ijob, m, n, &a[a_offset], lda, &b[b_offset], ldb, &c__[
	    c_offset], ldc, &d__[d_offset], ldd, &e[e_offset], lde, &f[
	    f_offset], ldf, scale, dif, &dwork[1], ldwork, &iwork[1], info, (
	    ftnlen)1);
#line 863 "SB04OD.f"
    if (*info != 0) {
#line 864 "SB04OD.f"
	*info = 3;
#line 865 "SB04OD.f"
	return 0;
#line 866 "SB04OD.f"
    }
#line 867 "SB04OD.f"
    if (ltrann) {
#line 868 "SB04OD.f"
	if (ljobd || ljobf) {
/* Computing MAX */
#line 868 "SB04OD.f"
	    i__1 = wrkopt, i__2 = (*m << 1) * *n;
#line 868 "SB04OD.f"
	    wrkopt = max(i__1,i__2);
#line 868 "SB04OD.f"
	}
#line 870 "SB04OD.f"
    }

/*     STEP 4: Back transformation of the solution. */

#line 874 "SB04OD.f"
    if (lreduc) {
#line 875 "SB04OD.f"
	if (sufwrk) {

/*           Enough workspace for a BLAS 3 calculation. */

#line 879 "SB04OD.f"
	    if (ltrann) {

/*              Equation (1). */

#line 883 "SB04OD.f"
		if (! lredub) {
#line 884 "SB04OD.f"
		    dgemm_("No transpose", "No transpose", m, n, m, &c_b52, &
			    q[q_offset], ldq, &c__[c_offset], ldc, &c_b53, &
			    dwork[1], m, (ftnlen)12, (ftnlen)12);
#line 886 "SB04OD.f"
		} else {
#line 887 "SB04OD.f"
		    dlacpy_("Full", m, n, &c__[c_offset], ldc, &dwork[1], m, (
			    ftnlen)4);
#line 888 "SB04OD.f"
		}
#line 889 "SB04OD.f"
		if (! lredua) {
#line 890 "SB04OD.f"
		    dgemm_("No transpose", "Transpose", m, n, n, &c_b52, &
			    dwork[1], m, &v[v_offset], ldv, &c_b53, &c__[
			    c_offset], ldc, (ftnlen)12, (ftnlen)9);
#line 892 "SB04OD.f"
		} else {
#line 893 "SB04OD.f"
		    dlacpy_("Full", m, n, &dwork[1], m, &c__[c_offset], ldc, (
			    ftnlen)4);
#line 894 "SB04OD.f"
		}
#line 895 "SB04OD.f"
		if (! lredub) {
#line 896 "SB04OD.f"
		    dgemm_("No transpose", "No transpose", m, n, m, &c_b52, &
			    p[p_offset], ldp, &f[f_offset], ldf, &c_b53, &
			    dwork[1], m, (ftnlen)12, (ftnlen)12);
#line 898 "SB04OD.f"
		} else {
#line 899 "SB04OD.f"
		    dlacpy_("Full", m, n, &f[f_offset], ldf, &dwork[1], m, (
			    ftnlen)4);
#line 900 "SB04OD.f"
		}
#line 901 "SB04OD.f"
		if (! lredua) {
#line 902 "SB04OD.f"
		    dgemm_("No transpose", "Transpose", m, n, n, &c_b52, &
			    dwork[1], m, &u[u_offset], ldu, &c_b53, &f[
			    f_offset], ldf, (ftnlen)12, (ftnlen)9);
#line 904 "SB04OD.f"
		} else {
#line 905 "SB04OD.f"
		    dlacpy_("Full", m, n, &dwork[1], m, &f[f_offset], ldf, (
			    ftnlen)4);
#line 906 "SB04OD.f"
		}
#line 907 "SB04OD.f"
	    } else {

/*              Equation (2). */

#line 911 "SB04OD.f"
		if (! lredub) {
#line 912 "SB04OD.f"
		    dgemm_("No transpose", "No transpose", m, n, m, &c_b52, &
			    p[p_offset], ldp, &c__[c_offset], ldc, &c_b53, &
			    dwork[1], m, (ftnlen)12, (ftnlen)12);
#line 914 "SB04OD.f"
		} else {
#line 915 "SB04OD.f"
		    dlacpy_("Full", m, n, &c__[c_offset], ldc, &dwork[1], m, (
			    ftnlen)4);
#line 916 "SB04OD.f"
		}
#line 917 "SB04OD.f"
		if (! lredua) {
#line 918 "SB04OD.f"
		    dgemm_("No transpose", "Transpose", m, n, n, &c_b52, &
			    dwork[1], m, &v[v_offset], ldv, &c_b53, &c__[
			    c_offset], ldc, (ftnlen)12, (ftnlen)9);
#line 920 "SB04OD.f"
		} else {
#line 921 "SB04OD.f"
		    dlacpy_("Full", m, n, &dwork[1], m, &c__[c_offset], ldc, (
			    ftnlen)4);
#line 922 "SB04OD.f"
		}
#line 923 "SB04OD.f"
		if (! lredub) {
#line 924 "SB04OD.f"
		    dgemm_("No transpose", "No transpose", m, n, m, &c_b52, &
			    p[p_offset], ldp, &f[f_offset], ldf, &c_b53, &
			    dwork[1], m, (ftnlen)12, (ftnlen)12);
#line 926 "SB04OD.f"
		} else {
#line 927 "SB04OD.f"
		    dlacpy_("Full", m, n, &f[f_offset], ldf, &dwork[1], m, (
			    ftnlen)4);
#line 928 "SB04OD.f"
		}
#line 929 "SB04OD.f"
		if (! lredua) {
#line 930 "SB04OD.f"
		    dgemm_("No transpose", "Transpose", m, n, n, &c_b52, &
			    dwork[1], m, &v[v_offset], ldv, &c_b53, &f[
			    f_offset], ldf, (ftnlen)12, (ftnlen)9);
#line 932 "SB04OD.f"
		} else {
#line 933 "SB04OD.f"
		    dlacpy_("Full", m, n, &dwork[1], m, &f[f_offset], ldf, (
			    ftnlen)4);
#line 934 "SB04OD.f"
		}
#line 935 "SB04OD.f"
	    }
#line 936 "SB04OD.f"
	} else {

/*           Use a BLAS 2 calculation. */

#line 940 "SB04OD.f"
	    if (ltrann) {

/*              Equation (1). */

#line 944 "SB04OD.f"
		if (! lredub) {

#line 946 "SB04OD.f"
		    i__1 = *n;
#line 946 "SB04OD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 947 "SB04OD.f"
			dgemv_("No transpose", m, m, &c_b52, &q[q_offset], 
				ldq, &c__[i__ * c_dim1 + 1], &c__1, &c_b53, &
				dwork[1], &c__1, (ftnlen)12);
#line 949 "SB04OD.f"
			dcopy_(m, &dwork[1], &c__1, &c__[i__ * c_dim1 + 1], &
				c__1);
#line 950 "SB04OD.f"
/* L220: */
#line 950 "SB04OD.f"
		    }

#line 952 "SB04OD.f"
		}
#line 953 "SB04OD.f"
		if (! lredua) {

#line 955 "SB04OD.f"
		    i__1 = *m;
#line 955 "SB04OD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 956 "SB04OD.f"
			dgemv_("No transpose", n, n, &c_b52, &v[v_offset], 
				ldv, &c__[i__ + c_dim1], ldc, &c_b53, &dwork[
				1], &c__1, (ftnlen)12);
#line 958 "SB04OD.f"
			dcopy_(n, &dwork[1], &c__1, &c__[i__ + c_dim1], ldc);
#line 959 "SB04OD.f"
/* L240: */
#line 959 "SB04OD.f"
		    }

#line 961 "SB04OD.f"
		}
#line 962 "SB04OD.f"
		if (! lredub) {

#line 964 "SB04OD.f"
		    i__1 = *n;
#line 964 "SB04OD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 965 "SB04OD.f"
			dgemv_("No transpose", m, m, &c_b52, &p[p_offset], 
				ldp, &f[i__ * f_dim1 + 1], &c__1, &c_b53, &
				dwork[1], &c__1, (ftnlen)12);
#line 967 "SB04OD.f"
			dcopy_(m, &dwork[1], &c__1, &f[i__ * f_dim1 + 1], &
				c__1);
#line 968 "SB04OD.f"
/* L260: */
#line 968 "SB04OD.f"
		    }

#line 970 "SB04OD.f"
		}
#line 971 "SB04OD.f"
		if (! lredua) {

#line 973 "SB04OD.f"
		    i__1 = *m;
#line 973 "SB04OD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 974 "SB04OD.f"
			dgemv_("No transpose", n, n, &c_b52, &u[u_offset], 
				ldu, &f[i__ + f_dim1], ldf, &c_b53, &dwork[1],
				 &c__1, (ftnlen)12);
#line 976 "SB04OD.f"
			dcopy_(n, &dwork[1], &c__1, &f[i__ + f_dim1], ldf);
#line 977 "SB04OD.f"
/* L280: */
#line 977 "SB04OD.f"
		    }

#line 979 "SB04OD.f"
		}
#line 980 "SB04OD.f"
	    } else {

/*              Equation (2). */

#line 984 "SB04OD.f"
		if (! lredub) {

#line 986 "SB04OD.f"
		    i__1 = *n;
#line 986 "SB04OD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 987 "SB04OD.f"
			dgemv_("No transpose", m, m, &c_b52, &p[p_offset], 
				ldp, &c__[i__ * c_dim1 + 1], &c__1, &c_b53, &
				dwork[1], &c__1, (ftnlen)12);
#line 989 "SB04OD.f"
			dcopy_(m, &dwork[1], &c__1, &c__[i__ * c_dim1 + 1], &
				c__1);
#line 990 "SB04OD.f"
/* L300: */
#line 990 "SB04OD.f"
		    }

#line 992 "SB04OD.f"
		}
#line 993 "SB04OD.f"
		if (! lredua) {

#line 995 "SB04OD.f"
		    i__1 = *m;
#line 995 "SB04OD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 996 "SB04OD.f"
			dgemv_("No transpose", n, n, &c_b52, &v[v_offset], 
				ldv, &c__[i__ + c_dim1], ldc, &c_b53, &dwork[
				1], &c__1, (ftnlen)12);
#line 998 "SB04OD.f"
			dcopy_(n, &dwork[1], &c__1, &c__[i__ + c_dim1], ldc);
#line 999 "SB04OD.f"
/* L320: */
#line 999 "SB04OD.f"
		    }

#line 1001 "SB04OD.f"
		}
#line 1002 "SB04OD.f"
		if (! lredub) {

#line 1004 "SB04OD.f"
		    i__1 = *n;
#line 1004 "SB04OD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1005 "SB04OD.f"
			dgemv_("No transpose", m, m, &c_b52, &p[p_offset], 
				ldp, &f[i__ * f_dim1 + 1], &c__1, &c_b53, &
				dwork[1], &c__1, (ftnlen)12);
#line 1007 "SB04OD.f"
			dcopy_(m, &dwork[1], &c__1, &f[i__ * f_dim1 + 1], &
				c__1);
#line 1008 "SB04OD.f"
/* L340: */
#line 1008 "SB04OD.f"
		    }

#line 1010 "SB04OD.f"
		}
#line 1011 "SB04OD.f"
		if (! lredua) {

#line 1013 "SB04OD.f"
		    i__1 = *m;
#line 1013 "SB04OD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1014 "SB04OD.f"
			dgemv_("No transpose", n, n, &c_b52, &v[v_offset], 
				ldv, &f[i__ + f_dim1], ldf, &c_b53, &dwork[1],
				 &c__1, (ftnlen)12);
#line 1016 "SB04OD.f"
			dcopy_(n, &dwork[1], &c__1, &f[i__ + f_dim1], ldf);
#line 1017 "SB04OD.f"
/* L360: */
#line 1017 "SB04OD.f"
		    }

#line 1019 "SB04OD.f"
		}
#line 1020 "SB04OD.f"
	    }
#line 1021 "SB04OD.f"
	}
#line 1022 "SB04OD.f"
    }

#line 1024 "SB04OD.f"
    dwork[1] = (doublereal) wrkopt;

#line 1026 "SB04OD.f"
    return 0;
/* *** Last line of SB04OD *** */
} /* sb04od_ */

