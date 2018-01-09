#line 1 "SB02OD.f"
/* SB02OD.f -- translated by f2c (version 20100827).
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

#line 1 "SB02OD.f"
/* Table of constant values */

static integer c__0 = 0;
static doublereal c_b26 = 1.;
static integer c__1 = 1;
static doublereal c_b66 = -1.;
static doublereal c_b104 = 0.;

/* Subroutine */ int sb02od_(char *dico, char *jobb, char *fact, char *uplo, 
	char *jobl, char *sort, integer *n, integer *m, integer *p, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	q, integer *ldq, doublereal *r__, integer *ldr, doublereal *l, 
	integer *ldl, doublereal *rcond, doublereal *x, integer *ldx, 
	doublereal *alfar, doublereal *alfai, doublereal *beta, doublereal *s,
	 integer *lds, doublereal *t, integer *ldt, doublereal *u, integer *
	ldu, doublereal *tol, integer *iwork, doublereal *dwork, integer *
	ldwork, logical *bwork, integer *info, ftnlen dico_len, ftnlen 
	jobb_len, ftnlen fact_len, ftnlen uplo_len, ftnlen jobl_len, ftnlen 
	sort_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, l_dim1, l_offset, q_dim1, 
	    q_offset, r_dim1, r_offset, s_dim1, s_offset, t_dim1, t_offset, 
	    u_dim1, u_offset, x_dim1, x_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, nn, mp, np, np1, ldw, nnm;
    static doublereal dum[1];
    static integer ndim;
    static logical lscl;
    static integer info1;
    static logical lfacb, lfacn, lfacq, lfacr, ljobb;
    static doublereal scale;
    extern /* Subroutine */ int dgees_(char *, char *, L_fp, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen), dgges_(char *, char *, char *, L_fp, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    logical *, integer *, ftnlen, ftnlen, ftnlen), dscal_(integer *, 
	    doublereal *, doublereal *, integer *);
    static logical lscal;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical ljobl;
    static doublereal qscal;
    static logical discr;
    static doublereal rscal;
    extern logical sb02mr_(), sb02mv_(), sb02ou_(), sb02ov_(), sb02ow_();
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), sb02oy_(char *, char *, char *, char 
	    *, char *, char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen)
	    , daxpy_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical luplo;
    static doublereal rnorm, unorm;
    static char qtype[1], rtype[1];
    static logical lsort;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgecon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dgetrf_(integer *, integer *, 
	    doublereal *, integer *, integer *, integer *), dlacpy_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical ljobln;
    static doublereal rcondl;
    extern /* Subroutine */ int dgetrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
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

/*     To solve for X either the continuous-time algebraic Riccati */
/*     equation */
/*                              -1 */
/*        Q + A'X + XA - (L+XB)R  (L+XB)' = 0                       (1) */

/*     or the discrete-time algebraic Riccati equation */
/*                                     -1 */
/*        X = A'XA - (L+A'XB)(R + B'XB)  (L+A'XB)' + Q              (2) */

/*     where A, B, Q, R, and L are N-by-N, N-by-M, N-by-N, M-by-M and */
/*     N-by-M matrices, respectively, such that Q = C'C, R = D'D and */
/*     L = C'D; X is an N-by-N symmetric matrix. */
/*     The routine also returns the computed values of the closed-loop */
/*     spectrum of the system, i.e., the stable eigenvalues lambda(1), */
/*     ..., lambda(N) of the corresponding Hamiltonian or symplectic */
/*     pencil, in the continuous-time case or discrete-time case, */
/*     respectively. */
/*                              -1 */
/*     Optionally, matrix G = BR  B' may be given instead of B and R. */
/*     Other options include the case with Q and/or R given in a */
/*     factored form, Q = C'C, R = D'D, and with L a zero matrix. */

/*     The routine uses the method of deflating subspaces, based on */
/*     reordering the eigenvalues in a generalized Schur matrix pair. */
/*     A standard eigenproblem is solved in the continuous-time case */
/*     if G is given. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of Riccati equation to be solved as */
/*             follows: */
/*             = 'C':  Equation (1), continuous-time case; */
/*             = 'D':  Equation (2), discrete-time case. */

/*     JOBB    CHARACTER*1 */
/*             Specifies whether or not the matrix G is given, instead */
/*             of the matrices B and R, as follows: */
/*             = 'B':  B and R are given; */
/*             = 'G':  G is given. */

/*     FACT    CHARACTER*1 */
/*             Specifies whether or not the matrices Q and/or R (if */
/*             JOBB = 'B') are factored, as follows: */
/*             = 'N':  Not factored, Q and R are given; */
/*             = 'C':  C is given, and Q = C'C; */
/*             = 'D':  D is given, and R = D'D; */
/*             = 'B':  Both factors C and D are given, Q = C'C, R = D'D. */

/*     UPLO    CHARACTER*1 */
/*             If JOBB = 'G', or FACT = 'N', specifies which triangle of */
/*             the matrices G and Q (if FACT = 'N'), or Q and R (if */
/*             JOBB = 'B'), is stored, as follows: */
/*             = 'U':  Upper triangle is stored; */
/*             = 'L':  Lower triangle is stored. */

/*     JOBL    CHARACTER*1 */
/*             Specifies whether or not the matrix L is zero, as follows: */
/*             = 'Z':  L is zero; */
/*             = 'N':  L is nonzero. */
/*             JOBL is not used if JOBB = 'G' and JOBL = 'Z' is assumed. */
/*             SLICOT Library routine SB02MT should be called just before */
/*             SB02OD, for obtaining the results when JOBB = 'G' and */
/*             JOBL = 'N'. */

/*     SORT    CHARACTER*1 */
/*             Specifies which eigenvalues should be obtained in the top */
/*             of the generalized Schur form, as follows: */
/*             = 'S':  Stable   eigenvalues come first; */
/*             = 'U':  Unstable eigenvalues come first. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The actual state dimension, i.e. the order of the matrices */
/*             A, Q, and X, and the number of rows of the matrices B */
/*             and L.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs. If JOBB = 'B', M is the */
/*             order of the matrix R, and the number of columns of the */
/*             matrix B.  M >= 0. */
/*             M is not used if JOBB = 'G'. */

/*     P       (input) INTEGER */
/*             The number of system outputs. If FACT = 'C' or 'D' or 'B', */
/*             P is the number of rows of the matrices C and/or D. */
/*             P >= 0. */
/*             Otherwise, P is not used. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             state matrix A of the system. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,*) */
/*             If JOBB = 'B', the leading N-by-M part of this array must */
/*             contain the input matrix B of the system. */
/*             If JOBB = 'G', the leading N-by-N upper triangular part */
/*             (if UPLO = 'U') or lower triangular part (if UPLO = 'L') */
/*             of this array must contain the upper triangular part or */
/*             lower triangular part, respectively, of the matrix */
/*                   -1 */
/*             G = BR  B'. The stricly lower triangular part (if */
/*             UPLO = 'U') or stricly upper triangular part (if */
/*             UPLO = 'L') is not referenced. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     Q       (input) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             If FACT = 'N' or 'D', the leading N-by-N upper triangular */
/*             part (if UPLO = 'U') or lower triangular part (if UPLO = */
/*             'L') of this array must contain the upper triangular part */
/*             or lower triangular part, respectively, of the symmetric */
/*             state weighting matrix Q. The stricly lower triangular */
/*             part (if UPLO = 'U') or stricly upper triangular part (if */
/*             UPLO = 'L') is not referenced. */
/*             If JOBB = 'B', the triangular part of this array defined */
/*             by UPLO is modified internally, but is restored on exit. */
/*             If FACT = 'C' or 'B', the leading P-by-N part of this */
/*             array must contain the output matrix C of the system. */
/*             If JOBB = 'B', this part is modified internally, but is */
/*             restored on exit. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. */
/*             LDQ >= MAX(1,N) if FACT = 'N' or 'D', */
/*             LDQ >= MAX(1,P) if FACT = 'C' or 'B'. */

/*     R       (input) DOUBLE PRECISION array, dimension (LDR,M) */
/*             If FACT = 'N' or 'C', the leading M-by-M upper triangular */
/*             part (if UPLO = 'U') or lower triangular part (if UPLO = */
/*             'L') of this array must contain the upper triangular part */
/*             or lower triangular part, respectively, of the symmetric */
/*             input weighting matrix R. The stricly lower triangular */
/*             part (if UPLO = 'U') or stricly upper triangular part (if */
/*             UPLO = 'L') is not referenced. */
/*             The triangular part of this array defined by UPLO is */
/*             modified internally, but is restored on exit. */
/*             If FACT = 'D' or 'B', the leading P-by-M part of this */
/*             array must contain the direct transmission matrix D of the */
/*             system. This part is modified internally, but is restored */
/*             on exit. */
/*             If JOBB = 'G', this array is not referenced. */

/*     LDR     INTEGER */
/*             The leading dimension of array R. */
/*             LDR >= MAX(1,M) if JOBB = 'B' and FACT = 'N' or 'C'; */
/*             LDR >= MAX(1,P) if JOBB = 'B' and FACT = 'D' or 'B'; */
/*             LDR >= 1        if JOBB = 'G'. */

/*     L       (input) DOUBLE PRECISION array, dimension (LDL,M) */
/*             If JOBL = 'N' (and JOBB = 'B'), the leading N-by-M part of */
/*             this array must contain the cross weighting matrix L. */
/*             This part is modified internally, but is restored on exit. */
/*             If JOBL = 'Z' or JOBB = 'G', this array is not referenced. */

/*     LDL     INTEGER */
/*             The leading dimension of array L. */
/*             LDL >= MAX(1,N) if JOBL = 'N' and JOBB = 'B'; */
/*             LDL >= 1        if JOBL = 'Z' or  JOBB = 'G'. */

/*     RCOND   (output) DOUBLE PRECISION */
/*             An estimate of the reciprocal of the condition number (in */
/*             the 1-norm) of the N-th order system of algebraic */
/*             equations from which the solution matrix X is obtained. */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             The leading N-by-N part of this array contains the */
/*             solution matrix X of the problem. */

/*     LDX     INTEGER */
/*             The leading dimension of array X.  LDX >= MAX(1,N). */

/*     ALFAR   (output) DOUBLE PRECISION array, dimension (2*N) */
/*     ALFAI   (output) DOUBLE PRECISION array, dimension (2*N) */
/*     BETA    (output) DOUBLE PRECISION array, dimension (2*N) */
/*             The generalized eigenvalues of the 2N-by-2N matrix pair, */
/*             ordered as specified by SORT (if INFO = 0). For instance, */
/*             if SORT = 'S', the leading N elements of these arrays */
/*             contain the closed-loop spectrum of the system matrix */
/*             A - BF, where F is the optimal feedback matrix computed */
/*             based on the solution matrix X. Specifically, */
/*                lambda(k) = [ALFAR(k)+j*ALFAI(k)]/BETA(k) for */
/*             k = 1,2,...,N. */
/*             If DICO = 'C' and JOBB = 'G', the elements of BETA are */
/*             set to 1. */

/*     S       (output) DOUBLE PRECISION array, dimension (LDS,*) */
/*             The leading 2N-by-2N part of this array contains the */
/*             ordered real Schur form S of the first matrix in the */
/*             reduced matrix pencil associated to the optimal problem, */
/*             or of the corresponding Hamiltonian matrix, if DICO = 'C' */
/*             and JOBB = 'G'. That is, */

/*                    (S   S  ) */
/*                    ( 11  12) */
/*                S = (       ), */
/*                    (0   S  ) */
/*                    (     22) */

/*             where S  , S   and S   are N-by-N matrices. */
/*                    11   12      22 */
/*             Array S must have 2*N+M columns if JOBB = 'B', and 2*N */
/*             columns, otherwise. */

/*     LDS     INTEGER */
/*             The leading dimension of array S. */
/*             LDS >= MAX(1,2*N+M) if JOBB = 'B', */
/*             LDS >= MAX(1,2*N)   if JOBB = 'G'. */

/*     T       (output) DOUBLE PRECISION array, dimension (LDT,2*N) */
/*             If DICO = 'D' or JOBB = 'B', the leading 2N-by-2N part of */
/*             this array contains the ordered upper triangular form T of */
/*             the second matrix in the reduced matrix pencil associated */
/*             to the optimal problem. That is, */

/*                    (T   T  ) */
/*                    ( 11  12) */
/*                T = (       ), */
/*                    (0   T  ) */
/*                    (     22) */

/*             where T  , T   and T   are N-by-N matrices. */
/*                    11   12      22 */
/*             If DICO = 'C' and JOBB = 'G' this array is not referenced. */

/*     LDT     INTEGER */
/*             The leading dimension of array T. */
/*             LDT >= MAX(1,2*N+M) if JOBB = 'B', */
/*             LDT >= MAX(1,2*N)   if JOBB = 'G' and DICO = 'D', */
/*             LDT >= 1            if JOBB = 'G' and DICO = 'C'. */

/*     U       (output) DOUBLE PRECISION array, dimension (LDU,2*N) */
/*             The leading 2N-by-2N part of this array contains the right */
/*             transformation matrix U which reduces the 2N-by-2N matrix */
/*             pencil to the ordered generalized real Schur form (S,T), */
/*             or the Hamiltonian matrix to the ordered real Schur */
/*             form S, if DICO = 'C' and JOBB = 'G'. That is, */

/*                    (U   U  ) */
/*                    ( 11  12) */
/*                U = (       ), */
/*                    (U   U  ) */
/*                    ( 21  22) */

/*             where U  , U  , U   and U   are N-by-N matrices. */
/*                    11   12   21      22 */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= MAX(1,2*N). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used to test for near singularity of */
/*             the original matrix pencil, specifically of the triangular */
/*             factor obtained during the reduction process. If the user */
/*             sets TOL > 0, then the given value of TOL is used as a */
/*             lower bound for the reciprocal condition number of that */
/*             matrix; a matrix whose estimated condition number is less */
/*             than 1/TOL is considered to be nonsingular. If the user */
/*             sets TOL <= 0, then a default tolerance, defined by */
/*             TOLDEF = EPS, is used instead, where EPS is the machine */
/*             precision (see LAPACK Library routine DLAMCH). */
/*             This parameter is not referenced if JOBB = 'G'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK >= MAX(1,M,2*N) if JOBB = 'B', */
/*             LIWORK >= MAX(1,2*N)   if JOBB = 'G'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. If JOBB = 'B' and N > 0, DWORK(2) returns the */
/*             reciprocal of the condition number of the M-by-M lower */
/*             triangular matrix obtained after compressing the matrix */
/*             pencil of order 2N+M to obtain a pencil of order 2N. */
/*             If INFO = 0 or INFO = 6, DWORK(3) returns the scaling */
/*             factor used internally, which should multiply the */
/*             submatrix Y2 to recover X from the first N columns of U */
/*             (see METHOD). */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(3,6*N),                       if JOBB = 'G', */
/*                                                            DICO = 'C'; */
/*             LDWORK >= MAX(7*(2*N+1)+16,16*N),           if JOBB = 'G', */
/*                                                            DICO = 'D'; */
/*             LDWORK >= MAX(7*(2*N+1)+16,16*N,2*N+M,3*M), if JOBB = 'B'. */
/*             For optimum performance LDWORK should be larger. */

/*     BWORK   LOGICAL array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the computed extended matrix pencil is singular, */
/*                   possibly due to rounding errors; */
/*             = 2:  if the QZ (or QR) algorithm failed; */
/*             = 3:  if reordering of the (generalized) eigenvalues */
/*                   failed; */
/*             = 4:  if after reordering, roundoff changed values of */
/*                   some complex eigenvalues so that leading eigenvalues */
/*                   in the (generalized) Schur form no longer satisfy */
/*                   the stability condition; this could also be caused */
/*                   due to scaling; */
/*             = 5:  if the computed dimension of the solution does not */
/*                   equal N; */
/*             = 6:  if a singular matrix was encountered during the */
/*                   computation of the solution matrix X. */

/*     METHOD */

/*     The routine uses a variant of the method of deflating subspaces */
/*     proposed by van Dooren [1]. See also [2], [3]. */
/*     It is assumed that (A,B) is stabilizable and (C,A) is detectable. */
/*     Under these assumptions the algebraic Riccati equation is known to */
/*     have a unique non-negative definite solution. */
/*     The first step in the method of deflating subspaces is to form the */
/*     extended Hamiltonian matrices, dimension 2N + M given by */

/*           discrete-time                   continuous-time */

/*     |A   0   B|     |I   0   0|    |A   0   B|     |I   0   0| */
/*     |Q  -I   L| - z |0  -A'  0|,   |Q   A'  L| - s |0  -I   0|. */
/*     |L'  0   R|     |0  -B'  0|    |L'  B'  R|     |0   0   0| */

/*     Next, these pencils are compressed to a form (see [1]) */

/*        lambda x A  - B . */
/*                  f    f */

/*     This generalized eigenvalue problem is then solved using the QZ */
/*     algorithm and the stable deflating subspace Ys is determined. */
/*     If [Y1'|Y2']' is a basis for Ys, then the required solution is */
/*                       -1 */
/*            X = Y2 x Y1  . */
/*     A standard eigenvalue problem is solved using the QR algorithm in */
/*     the continuous-time case when G is given (DICO = 'C', JOBB = 'G'). */

/*     REFERENCES */

/*     [1] Van Dooren, P. */
/*         A Generalized Eigenvalue Approach for Solving Riccati */
/*         Equations. */
/*         SIAM J. Sci. Stat. Comp., 2, pp. 121-135, 1981. */

/*     [2] Mehrmann, V. */
/*         The Autonomous Linear Quadratic Control Problem. Theory and */
/*         Numerical Solution. */
/*         Lect. Notes in Control and Information Sciences, vol. 163, */
/*         Springer-Verlag, Berlin, 1991. */

/*     [3] Sima, V. */
/*         Algorithms for Linear-Quadratic Optimization. */
/*         Pure and Applied Mathematics: A Series of Monographs and */
/*         Textbooks, vol. 200, Marcel Dekker, Inc., New York, 1996. */

/*     NUMERICAL ASPECTS */

/*     This routine is particularly suited for systems where the matrix R */
/*     is ill-conditioned. Internal scaling is used. */

/*     FURTHER COMMENTS */

/*     To obtain a stabilizing solution of the algebraic Riccati */
/*     equations set SORT = 'S'. */

/*     The routine can also compute the anti-stabilizing solutions of */
/*     the algebraic Riccati equations, by specifying SORT = 'U'. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Sep. 1997. */
/*     Supersedes Release 2.0 routine SB02CD by T.G.J. Beelen, Philips, */
/*     Eindhoven, Holland. */

/*     REVISIONS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, May 1999, June 2002, */
/*     December 2002, January 2005. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, closed loop system, continuous-time */
/*     system, discrete-time system, optimal regulator, Schur form. */

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

#line 462 "SB02OD.f"
    /* Parameter adjustments */
#line 462 "SB02OD.f"
    a_dim1 = *lda;
#line 462 "SB02OD.f"
    a_offset = 1 + a_dim1;
#line 462 "SB02OD.f"
    a -= a_offset;
#line 462 "SB02OD.f"
    b_dim1 = *ldb;
#line 462 "SB02OD.f"
    b_offset = 1 + b_dim1;
#line 462 "SB02OD.f"
    b -= b_offset;
#line 462 "SB02OD.f"
    q_dim1 = *ldq;
#line 462 "SB02OD.f"
    q_offset = 1 + q_dim1;
#line 462 "SB02OD.f"
    q -= q_offset;
#line 462 "SB02OD.f"
    r_dim1 = *ldr;
#line 462 "SB02OD.f"
    r_offset = 1 + r_dim1;
#line 462 "SB02OD.f"
    r__ -= r_offset;
#line 462 "SB02OD.f"
    l_dim1 = *ldl;
#line 462 "SB02OD.f"
    l_offset = 1 + l_dim1;
#line 462 "SB02OD.f"
    l -= l_offset;
#line 462 "SB02OD.f"
    x_dim1 = *ldx;
#line 462 "SB02OD.f"
    x_offset = 1 + x_dim1;
#line 462 "SB02OD.f"
    x -= x_offset;
#line 462 "SB02OD.f"
    --alfar;
#line 462 "SB02OD.f"
    --alfai;
#line 462 "SB02OD.f"
    --beta;
#line 462 "SB02OD.f"
    s_dim1 = *lds;
#line 462 "SB02OD.f"
    s_offset = 1 + s_dim1;
#line 462 "SB02OD.f"
    s -= s_offset;
#line 462 "SB02OD.f"
    t_dim1 = *ldt;
#line 462 "SB02OD.f"
    t_offset = 1 + t_dim1;
#line 462 "SB02OD.f"
    t -= t_offset;
#line 462 "SB02OD.f"
    u_dim1 = *ldu;
#line 462 "SB02OD.f"
    u_offset = 1 + u_dim1;
#line 462 "SB02OD.f"
    u -= u_offset;
#line 462 "SB02OD.f"
    --iwork;
#line 462 "SB02OD.f"
    --dwork;
#line 462 "SB02OD.f"
    --bwork;
#line 462 "SB02OD.f"

#line 462 "SB02OD.f"
    /* Function Body */
#line 462 "SB02OD.f"
    *info = 0;
#line 463 "SB02OD.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 464 "SB02OD.f"
    ljobb = lsame_(jobb, "B", (ftnlen)1, (ftnlen)1);
#line 465 "SB02OD.f"
    lfacn = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 466 "SB02OD.f"
    lfacq = lsame_(fact, "C", (ftnlen)1, (ftnlen)1);
#line 467 "SB02OD.f"
    lfacr = lsame_(fact, "D", (ftnlen)1, (ftnlen)1);
#line 468 "SB02OD.f"
    lfacb = lsame_(fact, "B", (ftnlen)1, (ftnlen)1);
#line 469 "SB02OD.f"
    luplo = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 470 "SB02OD.f"
    lsort = lsame_(sort, "S", (ftnlen)1, (ftnlen)1);

#line 472 "SB02OD.f"
    nn = *n << 1;
#line 473 "SB02OD.f"
    if (ljobb) {
#line 474 "SB02OD.f"
	ljobl = lsame_(jobl, "Z", (ftnlen)1, (ftnlen)1);
#line 475 "SB02OD.f"
	ljobln = lsame_(jobl, "N", (ftnlen)1, (ftnlen)1);
#line 476 "SB02OD.f"
	nnm = nn + *m;
/* Computing MAX */
#line 477 "SB02OD.f"
	i__1 = nnm, i__2 = *m * 3;
#line 477 "SB02OD.f"
	ldw = max(i__1,i__2);
#line 478 "SB02OD.f"
    } else {
#line 479 "SB02OD.f"
	nnm = nn;
#line 480 "SB02OD.f"
	ldw = 1;
#line 481 "SB02OD.f"
    }
#line 482 "SB02OD.f"
    np1 = *n + 1;

/*     Test the input scalar arguments. */

#line 486 "SB02OD.f"
    if (! discr && ! lsame_(dico, "C", (ftnlen)1, (ftnlen)1)) {
#line 487 "SB02OD.f"
	*info = -1;
#line 488 "SB02OD.f"
    } else if (! ljobb && ! lsame_(jobb, "G", (ftnlen)1, (ftnlen)1)) {
#line 489 "SB02OD.f"
	*info = -2;
#line 490 "SB02OD.f"
    } else if (! lfacq && ! lfacr && ! lfacb && ! lfacn) {
#line 492 "SB02OD.f"
	*info = -3;
#line 493 "SB02OD.f"
    } else if (! ljobb || lfacn) {
#line 494 "SB02OD.f"
	if (! luplo && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 494 "SB02OD.f"
	    *info = -4;
#line 494 "SB02OD.f"
	}
#line 496 "SB02OD.f"
    }
#line 497 "SB02OD.f"
    if (*info == 0 && ljobb) {
#line 498 "SB02OD.f"
	if (! ljobl && ! ljobln) {
#line 498 "SB02OD.f"
	    *info = -5;
#line 498 "SB02OD.f"
	}
#line 500 "SB02OD.f"
    }
#line 501 "SB02OD.f"
    if (*info == 0) {
#line 502 "SB02OD.f"
	if (! lsort && ! lsame_(sort, "U", (ftnlen)1, (ftnlen)1)) {
#line 503 "SB02OD.f"
	    *info = -6;
#line 504 "SB02OD.f"
	} else if (*n < 0) {
#line 505 "SB02OD.f"
	    *info = -7;
#line 506 "SB02OD.f"
	} else if (ljobb) {
#line 507 "SB02OD.f"
	    if (*m < 0) {
#line 507 "SB02OD.f"
		*info = -8;
#line 507 "SB02OD.f"
	    }
#line 509 "SB02OD.f"
	}
#line 510 "SB02OD.f"
    }
#line 511 "SB02OD.f"
    if (*info == 0 && ! lfacn) {
#line 512 "SB02OD.f"
	if (*p < 0) {
#line 512 "SB02OD.f"
	    *info = -9;
#line 512 "SB02OD.f"
	}
#line 514 "SB02OD.f"
    }
#line 515 "SB02OD.f"
    if (*info == 0) {
#line 516 "SB02OD.f"
	if (*lda < max(1,*n)) {
#line 517 "SB02OD.f"
	    *info = -11;
#line 518 "SB02OD.f"
	} else if (*ldb < max(1,*n)) {
#line 519 "SB02OD.f"
	    *info = -13;
#line 520 "SB02OD.f"
	} else if ((lfacn || lfacr) && *ldq < max(1,*n) || (lfacq || lfacb) &&
		 *ldq < max(1,*p)) {
#line 522 "SB02OD.f"
	    *info = -15;
#line 523 "SB02OD.f"
	} else if (*ldr < 1) {
#line 524 "SB02OD.f"
	    *info = -17;
#line 525 "SB02OD.f"
	} else if (*ldl < 1) {
#line 526 "SB02OD.f"
	    *info = -19;
#line 527 "SB02OD.f"
	} else if (ljobb) {
#line 528 "SB02OD.f"
	    if ((lfacn || lfacq) && *ldr < *m || (lfacr || lfacb) && *ldr < *
		    p) {
#line 530 "SB02OD.f"
		*info = -17;
#line 531 "SB02OD.f"
	    } else if (ljobln && *ldl < *n) {
#line 532 "SB02OD.f"
		*info = -19;
#line 533 "SB02OD.f"
	    }
#line 534 "SB02OD.f"
	}
#line 535 "SB02OD.f"
    }
#line 536 "SB02OD.f"
    if (*info == 0) {
#line 537 "SB02OD.f"
	if (*ldx < max(1,*n)) {
#line 538 "SB02OD.f"
	    *info = -22;
#line 539 "SB02OD.f"
	} else if (*lds < max(1,nnm)) {
#line 540 "SB02OD.f"
	    *info = -27;
#line 541 "SB02OD.f"
	} else if (*ldt < 1) {
#line 542 "SB02OD.f"
	    *info = -29;
#line 543 "SB02OD.f"
	} else if (*ldu < max(1,nn)) {
#line 544 "SB02OD.f"
	    *info = -31;
#line 545 "SB02OD.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 545 "SB02OD.f"
	    i__1 = 3, i__2 = *n * 6;
#line 545 "SB02OD.f"
	    if (*ldwork < max(i__1,i__2)) {
#line 546 "SB02OD.f"
		*info = -35;
#line 547 "SB02OD.f"
	    } else if (discr || ljobb) {
#line 548 "SB02OD.f"
		if (*ldt < nnm) {
#line 549 "SB02OD.f"
		    *info = -29;
#line 550 "SB02OD.f"
		} else /* if(complicated condition) */ {
/* Computing MAX */
#line 550 "SB02OD.f"
		    i__1 = *n * 14 + 23, i__2 = *n << 4, i__1 = max(i__1,i__2)
			    ;
#line 550 "SB02OD.f"
		    if (*ldwork < max(i__1,ldw)) {
#line 551 "SB02OD.f"
			*info = -35;
#line 552 "SB02OD.f"
		    }
#line 552 "SB02OD.f"
		}
#line 553 "SB02OD.f"
	    }
#line 553 "SB02OD.f"
	}
#line 554 "SB02OD.f"
    }

#line 556 "SB02OD.f"
    if (*info != 0) {

/*        Error return. */

#line 560 "SB02OD.f"
	i__1 = -(*info);
#line 560 "SB02OD.f"
	xerbla_("SB02OD", &i__1, (ftnlen)6);
#line 561 "SB02OD.f"
	return 0;
#line 562 "SB02OD.f"
    }

/*     Quick return if possible. */

#line 566 "SB02OD.f"
    if (*n == 0) {
#line 567 "SB02OD.f"
	*rcond = 1.;
#line 568 "SB02OD.f"
	dwork[1] = 3.;
#line 569 "SB02OD.f"
	dwork[3] = 1.;
#line 570 "SB02OD.f"
	return 0;
#line 571 "SB02OD.f"
    }

/*     Always scale the matrix pencil. */

#line 575 "SB02OD.f"
    lscal = TRUE_;

/*     Start computations. */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 585 "SB02OD.f"
    if (lscal && ljobb) {

/*        Scale the matrices Q, R, and L so that */
/*           norm(Q) + norm(R) + norm(L) = 1, */
/*        using the 1-norm. If Q and/or R are factored, the norms of */
/*        the factors are used. */
/*        Workspace: need   max(N,M), if FACT = 'N'; */
/*                          N,        if FACT = 'D'; */
/*                          M,        if FACT = 'C'. */

#line 595 "SB02OD.f"
	if (lfacn || lfacr) {
#line 596 "SB02OD.f"
	    scale = dlansy_("1-norm", uplo, n, &q[q_offset], ldq, &dwork[1], (
		    ftnlen)6, (ftnlen)1);
#line 597 "SB02OD.f"
	    *(unsigned char *)qtype = *(unsigned char *)uplo;
#line 598 "SB02OD.f"
	    np = *n;
#line 599 "SB02OD.f"
	} else {
#line 600 "SB02OD.f"
	    scale = dlange_("1-norm", p, n, &q[q_offset], ldq, &dwork[1], (
		    ftnlen)6);
#line 601 "SB02OD.f"
	    *(unsigned char *)qtype = 'G';
#line 602 "SB02OD.f"
	    np = *p;
#line 603 "SB02OD.f"
	}

#line 605 "SB02OD.f"
	if (lfacn || lfacq) {
#line 606 "SB02OD.f"
	    rnorm = dlansy_("1-norm", uplo, m, &r__[r_offset], ldr, &dwork[1],
		     (ftnlen)6, (ftnlen)1);
#line 607 "SB02OD.f"
	    *(unsigned char *)rtype = *(unsigned char *)uplo;
#line 608 "SB02OD.f"
	    mp = *m;
#line 609 "SB02OD.f"
	} else {
#line 610 "SB02OD.f"
	    rnorm = dlange_("1-norm", p, m, &r__[r_offset], ldr, &dwork[1], (
		    ftnlen)6);
#line 611 "SB02OD.f"
	    *(unsigned char *)rtype = 'G';
#line 612 "SB02OD.f"
	    mp = *p;
#line 613 "SB02OD.f"
	}
#line 614 "SB02OD.f"
	scale += rnorm;

#line 616 "SB02OD.f"
	if (ljobln) {
#line 616 "SB02OD.f"
	    scale += dlange_("1-norm", n, m, &l[l_offset], ldl, &dwork[1], (
		    ftnlen)6);
#line 616 "SB02OD.f"
	}
#line 618 "SB02OD.f"
	if (scale == 0.) {
#line 618 "SB02OD.f"
	    scale = 1.;
#line 618 "SB02OD.f"
	}

#line 621 "SB02OD.f"
	if (lfacn || lfacr) {
#line 622 "SB02OD.f"
	    qscal = scale;
#line 623 "SB02OD.f"
	} else {
#line 624 "SB02OD.f"
	    qscal = sqrt(scale);
#line 625 "SB02OD.f"
	}

#line 627 "SB02OD.f"
	if (lfacn || lfacq) {
#line 628 "SB02OD.f"
	    rscal = scale;
#line 629 "SB02OD.f"
	} else {
#line 630 "SB02OD.f"
	    rscal = sqrt(scale);
#line 631 "SB02OD.f"
	}

#line 633 "SB02OD.f"
	dlascl_(qtype, &c__0, &c__0, &qscal, &c_b26, &np, n, &q[q_offset], 
		ldq, &info1, (ftnlen)1);
#line 634 "SB02OD.f"
	dlascl_(rtype, &c__0, &c__0, &rscal, &c_b26, &mp, m, &r__[r_offset], 
		ldr, &info1, (ftnlen)1);
#line 635 "SB02OD.f"
	if (ljobln) {
#line 635 "SB02OD.f"
	    dlascl_("G", &c__0, &c__0, &scale, &c_b26, n, m, &l[l_offset], 
		    ldl, &info1, (ftnlen)1);
#line 635 "SB02OD.f"
	}
#line 637 "SB02OD.f"
    }

/*     Construct the extended matrix pair. */

/*     Workspace: need   1,                if JOBB = 'G', */
/*                       max(1,2*N+M,3*M), if JOBB = 'B'; */
/*                prefer larger. */

#line 645 "SB02OD.f"
    sb02oy_("Optimal control", dico, jobb, fact, uplo, jobl, "Identity E", n, 
	    m, p, &a[a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &
	    r__[r_offset], ldr, &l[l_offset], ldl, &u[u_offset], &c__1, &s[
	    s_offset], lds, &t[t_offset], ldt, tol, &iwork[1], &dwork[1], 
	    ldwork, info, (ftnlen)15, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
	    ftnlen)1, (ftnlen)1, (ftnlen)10);

#line 650 "SB02OD.f"
    if (lscal && ljobb) {

/*        Undo scaling of the data arrays. */

#line 654 "SB02OD.f"
	dlascl_(qtype, &c__0, &c__0, &c_b26, &qscal, &np, n, &q[q_offset], 
		ldq, &info1, (ftnlen)1);
#line 655 "SB02OD.f"
	dlascl_(rtype, &c__0, &c__0, &c_b26, &rscal, &mp, m, &r__[r_offset], 
		ldr, &info1, (ftnlen)1);
#line 656 "SB02OD.f"
	if (ljobln) {
#line 656 "SB02OD.f"
	    dlascl_("G", &c__0, &c__0, &c_b26, &scale, n, m, &l[l_offset], 
		    ldl, &info1, (ftnlen)1);
#line 656 "SB02OD.f"
	}
#line 658 "SB02OD.f"
    }

#line 660 "SB02OD.f"
    if (*info != 0) {
#line 660 "SB02OD.f"
	return 0;
#line 660 "SB02OD.f"
    }
#line 662 "SB02OD.f"
    wrkopt = (integer) dwork[1];
#line 663 "SB02OD.f"
    if (ljobb) {
#line 663 "SB02OD.f"
	rcondl = dwork[2];
#line 663 "SB02OD.f"
    }

#line 665 "SB02OD.f"
    if (lscal && ! ljobb) {

/*        This part of the code is used when G is given (JOBB = 'G'). */
/*        A standard eigenproblem is solved in the continuous-time case. */
/*        Scale the Hamiltonian matrix S, if DICO = 'C', or the */
/*        symplectic pencil (S,T), if DICO = 'D', using the square roots */
/*        of the norms of the matrices Q and G. */
/*        Workspace: need   N. */

#line 674 "SB02OD.f"
	if (lfacn || lfacr) {
#line 675 "SB02OD.f"
	    scale = sqrt(dlansy_("1-norm", uplo, n, &q[q_offset], ldq, &dwork[
		    1], (ftnlen)6, (ftnlen)1));
#line 676 "SB02OD.f"
	} else {
#line 677 "SB02OD.f"
	    scale = dlange_("1-norm", p, n, &q[q_offset], ldq, &dwork[1], (
		    ftnlen)6);
#line 678 "SB02OD.f"
	}
#line 679 "SB02OD.f"
	rnorm = sqrt(dlansy_("1-norm", uplo, n, &b[b_offset], ldb, &dwork[1], 
		(ftnlen)6, (ftnlen)1));

#line 681 "SB02OD.f"
	lscl = min(scale,rnorm) > 0. && scale != rnorm;

#line 683 "SB02OD.f"
	if (lscl) {
#line 684 "SB02OD.f"
	    if (discr) {
#line 685 "SB02OD.f"
		dlascl_("G", &c__0, &c__0, &scale, &rnorm, n, n, &s[np1 + 
			s_dim1], lds, &info1, (ftnlen)1);
#line 687 "SB02OD.f"
		dlascl_("G", &c__0, &c__0, &rnorm, &scale, n, n, &t[np1 * 
			t_dim1 + 1], ldt, &info1, (ftnlen)1);
#line 689 "SB02OD.f"
	    } else {
#line 690 "SB02OD.f"
		d__1 = -rnorm;
#line 690 "SB02OD.f"
		dlascl_("G", &c__0, &c__0, &scale, &d__1, n, n, &s[np1 + 
			s_dim1], lds, &info1, (ftnlen)1);
#line 692 "SB02OD.f"
		dlascl_("G", &c__0, &c__0, &rnorm, &scale, n, n, &s[np1 * 
			s_dim1 + 1], lds, &info1, (ftnlen)1);
#line 694 "SB02OD.f"
		dlascl_("G", &c__0, &c__0, &c_b26, &c_b66, n, n, &s[np1 + np1 
			* s_dim1], lds, &info1, (ftnlen)1);
#line 696 "SB02OD.f"
	    }
#line 697 "SB02OD.f"
	} else {
#line 698 "SB02OD.f"
	    if (! discr) {
#line 699 "SB02OD.f"
		dlascl_("G", &c__0, &c__0, &c_b26, &c_b66, n, &nn, &s[np1 + 
			s_dim1], lds, &info1, (ftnlen)1);
#line 701 "SB02OD.f"
	    }
#line 702 "SB02OD.f"
	}
#line 703 "SB02OD.f"
    } else {
#line 704 "SB02OD.f"
	lscl = FALSE_;
#line 705 "SB02OD.f"
    }

/*     Workspace: need   max(7*(2*N+1)+16,16*N), */
/*                                          if JOBB = 'B' or  DICO = 'D'; */
/*                       6*N,               if JOBB = 'G' and DICO = 'C'; */
/*                prefer larger. */

#line 712 "SB02OD.f"
    if (discr) {
#line 713 "SB02OD.f"
	if (lsort) {

/*           The natural tendency of the QZ algorithm to get the largest */
/*           eigenvalues in the leading part of the matrix pair is */
/*           exploited, by computing the unstable eigenvalues of the */
/*           permuted matrix pair. */

#line 720 "SB02OD.f"
	    dgges_("No vectors", "Vectors", "Sort", (L_fp)sb02ov_, &nn, &t[
		    t_offset], ldt, &s[s_offset], lds, &ndim, &alfar[1], &
		    alfai[1], &beta[1], &u[u_offset], ldu, &u[u_offset], ldu, 
		    &dwork[1], ldwork, &bwork[1], &info1, (ftnlen)10, (ftnlen)
		    7, (ftnlen)4);
#line 723 "SB02OD.f"
	    dswap_(n, &alfar[np1], &c__1, &alfar[1], &c__1);
#line 724 "SB02OD.f"
	    dswap_(n, &alfai[np1], &c__1, &alfai[1], &c__1);
#line 725 "SB02OD.f"
	    dswap_(n, &beta[np1], &c__1, &beta[1], &c__1);
#line 726 "SB02OD.f"
	} else {
#line 727 "SB02OD.f"
	    dgges_("No vectors", "Vectors", "Sort", (L_fp)sb02ov_, &nn, &s[
		    s_offset], lds, &t[t_offset], ldt, &ndim, &alfar[1], &
		    alfai[1], &beta[1], &u[u_offset], ldu, &u[u_offset], ldu, 
		    &dwork[1], ldwork, &bwork[1], &info1, (ftnlen)10, (ftnlen)
		    7, (ftnlen)4);
#line 730 "SB02OD.f"
	}
#line 731 "SB02OD.f"
    } else {
#line 732 "SB02OD.f"
	if (ljobb) {
#line 733 "SB02OD.f"
	    if (lsort) {
#line 734 "SB02OD.f"
		dgges_("No vectors", "Vectors", "Sort", (L_fp)sb02ow_, &nn, &
			s[s_offset], lds, &t[t_offset], ldt, &ndim, &alfar[1],
			 &alfai[1], &beta[1], &u[u_offset], ldu, &u[u_offset],
			 ldu, &dwork[1], ldwork, &bwork[1], &info1, (ftnlen)
			10, (ftnlen)7, (ftnlen)4);
#line 737 "SB02OD.f"
	    } else {
#line 738 "SB02OD.f"
		dgges_("No vectors", "Vectors", "Sort", (L_fp)sb02ou_, &nn, &
			s[s_offset], lds, &t[t_offset], ldt, &ndim, &alfar[1],
			 &alfai[1], &beta[1], &u[u_offset], ldu, &u[u_offset],
			 ldu, &dwork[1], ldwork, &bwork[1], &info1, (ftnlen)
			10, (ftnlen)7, (ftnlen)4);
#line 741 "SB02OD.f"
	    }
#line 742 "SB02OD.f"
	} else {
#line 743 "SB02OD.f"
	    if (lsort) {
#line 744 "SB02OD.f"
		dgees_("Vectors", "Sort", (L_fp)sb02mv_, &nn, &s[s_offset], 
			lds, &ndim, &alfar[1], &alfai[1], &u[u_offset], ldu, &
			dwork[1], ldwork, &bwork[1], &info1, (ftnlen)7, (
			ftnlen)4);
#line 747 "SB02OD.f"
	    } else {
#line 748 "SB02OD.f"
		dgees_("Vectors", "Sort", (L_fp)sb02mr_, &nn, &s[s_offset], 
			lds, &ndim, &alfar[1], &alfai[1], &u[u_offset], ldu, &
			dwork[1], ldwork, &bwork[1], &info1, (ftnlen)7, (
			ftnlen)4);
#line 751 "SB02OD.f"
	    }
#line 752 "SB02OD.f"
	    dum[0] = 1.;
#line 753 "SB02OD.f"
	    dcopy_(&nn, dum, &c__0, &beta[1], &c__1);
#line 754 "SB02OD.f"
	}
#line 755 "SB02OD.f"
    }
#line 756 "SB02OD.f"
    if (info1 > 0 && info1 <= nn + 1) {
#line 757 "SB02OD.f"
	*info = 2;
#line 758 "SB02OD.f"
    } else if (info1 == nn + 2) {
#line 759 "SB02OD.f"
	*info = 4;
#line 760 "SB02OD.f"
    } else if (info1 == nn + 3) {
#line 761 "SB02OD.f"
	*info = 3;
#line 762 "SB02OD.f"
    } else if (ndim != *n) {
#line 763 "SB02OD.f"
	*info = 5;
#line 764 "SB02OD.f"
    }
#line 765 "SB02OD.f"
    if (*info != 0) {
#line 765 "SB02OD.f"
	return 0;
#line 765 "SB02OD.f"
    }
/* Computing MAX */
#line 767 "SB02OD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 767 "SB02OD.f"
    wrkopt = max(i__1,i__2);

/*     Select submatrices U1 and U2 out of the array U which define the */
/*     solution X = U2 x inv(U1). */
/*     Since X = X' we may obtain X as the solution of the system of */
/*     linear equations U1' x X = U2', where */
/*        U1 = U(1:n, 1:n), */
/*        U2 = U(n+1:2n, 1:n). */
/*     Use the (2,1) block of S as a workspace for factoring U1. */

#line 777 "SB02OD.f"
    i__1 = *n;
#line 777 "SB02OD.f"
    for (j = 1; j <= i__1; ++j) {
#line 778 "SB02OD.f"
	dcopy_(n, &u[np1 + j * u_dim1], &c__1, &x[j + x_dim1], ldx);
#line 779 "SB02OD.f"
/* L20: */
#line 779 "SB02OD.f"
    }

#line 781 "SB02OD.f"
    dlacpy_("Full", n, n, &u[u_offset], ldu, &s[np1 + s_dim1], lds, (ftnlen)4)
	    ;

/*     Check if U1 is singular. */

#line 785 "SB02OD.f"
    unorm = dlange_("1-norm", n, n, &s[np1 + s_dim1], lds, &dwork[1], (ftnlen)
	    6);

/*     Solve the system U1' x X = U2'. */

#line 789 "SB02OD.f"
    dgetrf_(n, n, &s[np1 + s_dim1], lds, &iwork[1], &info1);
#line 790 "SB02OD.f"
    if (info1 != 0) {
#line 791 "SB02OD.f"
	*info = 6;
#line 792 "SB02OD.f"
	dwork[3] = 1.;
#line 793 "SB02OD.f"
	if (lscal) {
#line 794 "SB02OD.f"
	    if (ljobb) {
#line 795 "SB02OD.f"
		dwork[3] = scale;
#line 796 "SB02OD.f"
	    } else if (lscl) {
#line 797 "SB02OD.f"
		dwork[3] = scale / rnorm;
#line 798 "SB02OD.f"
	    }
#line 799 "SB02OD.f"
	}
#line 800 "SB02OD.f"
	return 0;
#line 801 "SB02OD.f"
    } else {

/*        Estimate the reciprocal condition of U1. */
/*        Workspace: need 3*N. */

#line 806 "SB02OD.f"
	dgecon_("1-norm", n, &s[np1 + s_dim1], lds, &unorm, rcond, &dwork[1], 
		&iwork[np1], info, (ftnlen)6);

#line 809 "SB02OD.f"
	if (*rcond < dlamch_("Epsilon", (ftnlen)7)) {

/*           Nearly singular matrix.  Set INFO for error return. */

#line 813 "SB02OD.f"
	    *info = 6;
#line 814 "SB02OD.f"
	    return 0;
#line 815 "SB02OD.f"
	}
/* Computing MAX */
#line 816 "SB02OD.f"
	i__1 = wrkopt, i__2 = *n * 3;
#line 816 "SB02OD.f"
	wrkopt = max(i__1,i__2);
#line 817 "SB02OD.f"
	dgetrs_("Transpose", n, n, &s[np1 + s_dim1], lds, &iwork[1], &x[
		x_offset], ldx, &info1, (ftnlen)9);

/*        Set S(2,1) to zero. */

#line 822 "SB02OD.f"
	dlaset_("Full", n, n, &c_b104, &c_b104, &s[np1 + s_dim1], lds, (
		ftnlen)4);

#line 824 "SB02OD.f"
	if (lscal) {

/*           Prepare to undo scaling for the solution X. */

#line 828 "SB02OD.f"
	    if (! ljobb) {
#line 829 "SB02OD.f"
		if (lscl) {
#line 830 "SB02OD.f"
		    scale /= rnorm;
#line 831 "SB02OD.f"
		} else {
#line 832 "SB02OD.f"
		    scale = 1.;
#line 833 "SB02OD.f"
		}
#line 834 "SB02OD.f"
	    }
#line 835 "SB02OD.f"
	    dwork[3] = scale;
#line 836 "SB02OD.f"
	    scale *= .5;
#line 837 "SB02OD.f"
	} else {
#line 838 "SB02OD.f"
	    dwork[3] = 1.;
#line 839 "SB02OD.f"
	    scale = .5;
#line 840 "SB02OD.f"
	}

/*        Make sure the solution matrix X is symmetric. */

#line 844 "SB02OD.f"
	i__1 = *n;
#line 844 "SB02OD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 845 "SB02OD.f"
	    i__2 = *n - i__ + 1;
#line 845 "SB02OD.f"
	    daxpy_(&i__2, &c_b26, &x[i__ + i__ * x_dim1], ldx, &x[i__ + i__ * 
		    x_dim1], &c__1);
#line 846 "SB02OD.f"
	    i__2 = *n - i__ + 1;
#line 846 "SB02OD.f"
	    dscal_(&i__2, &scale, &x[i__ + i__ * x_dim1], &c__1);
#line 847 "SB02OD.f"
	    i__2 = *n - i__ + 1;
#line 847 "SB02OD.f"
	    dcopy_(&i__2, &x[i__ + i__ * x_dim1], &c__1, &x[i__ + i__ * 
		    x_dim1], ldx);
#line 848 "SB02OD.f"
/* L40: */
#line 848 "SB02OD.f"
	}
#line 849 "SB02OD.f"
    }

#line 851 "SB02OD.f"
    dwork[1] = (doublereal) wrkopt;
#line 852 "SB02OD.f"
    if (ljobb) {
#line 852 "SB02OD.f"
	dwork[2] = rcondl;
#line 852 "SB02OD.f"
    }

#line 854 "SB02OD.f"
    return 0;
/* *** Last line of SB02OD *** */
} /* sb02od_ */
