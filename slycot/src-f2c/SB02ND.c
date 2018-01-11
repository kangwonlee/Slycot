#line 1 "SB02ND.f"
/* SB02ND.f -- translated by f2c (version 20100827).
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

#line 1 "SB02ND.f"
/* Table of constant values */

static doublereal c_b16 = 1.;
static doublereal c_b17 = 0.;
static integer c__1 = 1;
static integer c__0 = 0;

/* Subroutine */ int sb02nd_(char *dico, char *fact, char *uplo, char *jobl, 
	integer *n, integer *m, integer *p, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *r__, integer *ldr, integer *
	ipiv, doublereal *l, integer *ldl, doublereal *x, integer *ldx, 
	doublereal *rnorm, doublereal *f, integer *ldf, integer *oufact, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen dico_len, ftnlen fact_len, ftnlen uplo_len, ftnlen jobl_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, f_dim1, f_offset, l_dim1, 
	    l_offset, r_dim1, r_offset, x_dim1, x_offset, i__1, i__2, i__3, 
	    i__4, i__5, i__6, i__7, i__8;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, jw, jz;
    static doublereal eps;
    static integer itau;
    static doublereal temp;
    extern /* Subroutine */ int mb04kd_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    static integer ifail;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static logical discr;
    static doublereal rcond;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical withl;
    extern /* Subroutine */ int dsyev_(char *, char *, integer *, doublereal *
	    , integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static integer jwork;
    static doublereal dummy[1];
    extern doublereal dlamch_(char *, ftnlen);
    static logical lfacta, lfactc, lfactd;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static logical lfactu;
    extern /* Subroutine */ int dpocon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dtrcon_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int dpotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dsycon_(char *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen), dpotrs_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dsytrf_(char *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static doublereal rnormp;
    static logical luplou;
    static integer wrkopt;
    extern /* Subroutine */ int dsytrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);


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

/*     To compute the optimal feedback matrix F for the problem of */
/*     optimal control given by */

/*                        -1 */
/*          F = (R + B'XB)  (B'XA + L')                           (1) */

/*     in the discrete-time case and */

/*               -1 */
/*          F = R  (B'X + L')                                     (2) */

/*     in the continuous-time case, where A, B and L are N-by-N, N-by-M */
/*     and N-by-M matrices respectively; R and X are M-by-M and N-by-N */
/*     symmetric matrices respectively. */

/*     Optionally, matrix R may be specified in a factored form, and L */
/*     may be zero. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the equation from which F is to be determined, */
/*             as follows: */
/*             = 'D':  Equation (1), discrete-time case; */
/*             = 'C':  Equation (2), continuous-time case. */

/*     FACT    CHARACTER*1 */
/*             Specifies how the matrix R is given (factored or not), as */
/*             follows: */
/*             = 'N':  Array R contains the matrix R; */
/*             = 'D':  Array R contains a P-by-M matrix D, where R = D'D; */
/*             = 'C':  Array R contains the Cholesky factor of R; */
/*             = 'U':  Array R contains the symmetric indefinite UdU' or */
/*                     LdL' factorization of R. This option is not */
/*                     available for DICO = 'D'. */

/*     UPLO    CHARACTER*1 */
/*             Specifies which triangle of the possibly factored matrix R */
/*             (or R + B'XB, on exit) is or should be stored, as follows: */
/*             = 'U':  Upper triangle is stored; */
/*             = 'L':  Lower triangle is stored. */

/*     JOBL    CHARACTER*1 */
/*             Specifies whether or not the matrix L is zero, as follows: */
/*             = 'Z':  L is zero; */
/*             = 'N':  L is nonzero. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A and X.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */
/*             This parameter must be specified only for FACT = 'D'. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             If DICO = 'D', the leading N-by-N part of this array must */
/*             contain the state matrix A of the system. */
/*             If DICO = 'C', this array is not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of array A. */
/*             LDA >= MAX(1,N) if DICO = 'D'; */
/*             LDA >= 1        if DICO = 'C'. */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             input matrix B of the system. */
/*             If DICO = 'D' and FACT = 'D' or 'C', the contents of this */
/*             array is destroyed. */
/*             Otherwise, B is unchanged on exit. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR,M) */
/*             On entry, if FACT = 'N', the leading M-by-M upper */
/*             triangular part (if UPLO = 'U') or lower triangular part */
/*             (if UPLO = 'L') of this array must contain the upper */
/*             triangular part or lower triangular part, respectively, */
/*             of the symmetric input weighting matrix R. */
/*             On entry, if FACT = 'D', the leading P-by-M part of this */
/*             array must contain the direct transmission matrix D of the */
/*             system. */
/*             On entry, if FACT = 'C', the leading M-by-M upper */
/*             triangular part (if UPLO = 'U') or lower triangular part */
/*             (if UPLO = 'L') of this array must contain the Cholesky */
/*             factor of the positive definite input weighting matrix R */
/*             (as produced by LAPACK routine DPOTRF). */
/*             On entry, if DICO = 'C' and FACT = 'U', the leading M-by-M */
/*             upper triangular part (if UPLO = 'U') or lower triangular */
/*             part (if UPLO = 'L') of this array must contain the */
/*             factors of the UdU' or LdL' factorization, respectively, */
/*             of the symmetric indefinite input weighting matrix R (as */
/*             produced by LAPACK routine DSYTRF). */
/*             The stricly lower triangular part (if UPLO = 'U') or */
/*             stricly upper triangular part (if UPLO = 'L') of this */
/*             array is used as workspace. */
/*             On exit, if OUFACT(1) = 1, and INFO = 0 (or INFO = M+1), */
/*             the leading M-by-M upper triangular part (if UPLO = 'U') */
/*             or lower triangular part (if UPLO = 'L') of this array */
/*             contains the Cholesky factor of the given input weighting */
/*             matrix (for DICO = 'C'), or that of the matrix R + B'XB */
/*             (for DICO = 'D'). */
/*             On exit, if OUFACT(1) = 2, and INFO = 0 (or INFO = M+1), */
/*             the leading M-by-M upper triangular part (if UPLO = 'U') */
/*             or lower triangular part (if UPLO = 'L') of this array */
/*             contains the factors of the UdU' or LdL' factorization, */
/*             respectively, of the given input weighting matrix */
/*             (for DICO = 'C'), or that of the matrix R + B'XB */
/*             (for DICO = 'D'). */
/*             On exit R is unchanged if FACT = 'U'. */

/*     LDR     INTEGER. */
/*             The leading dimension of the array R. */
/*             LDR >= MAX(1,M)   if FACT <> 'D'; */
/*             LDR >= MAX(1,M,P) if FACT =  'D'. */

/*     IPIV    (input/output) INTEGER array, dimension (M) */
/*             On entry, if FACT = 'U', this array must contain details */
/*             of the interchanges performed and the block structure of */
/*             the d factor in the UdU' or LdL' factorization of matrix R */
/*             (as produced by LAPACK routine DSYTRF). */
/*             On exit, if OUFACT(1) = 2, this array contains details of */
/*             the interchanges performed and the block structure of the */
/*             d factor in the UdU' or LdL' factorization of matrix R (or */
/*             D'D) or R + B'XB (or D'D + B'XB), as produced by LAPACK */
/*             routine DSYTRF. */
/*             This array is not referenced for DICO = 'D' or FACT = 'D', */
/*             or 'C'. */

/*     L       (input) DOUBLE PRECISION array, dimension (LDL,M) */
/*             If JOBL = 'N', the leading N-by-M part of this array must */
/*             contain the cross weighting matrix L. */
/*             If JOBL = 'Z', this array is not referenced. */

/*     LDL     INTEGER */
/*             The leading dimension of array L. */
/*             LDL >= MAX(1,N) if JOBL = 'N'; */
/*             LDL >= 1        if JOBL = 'Z'. */

/*     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the solution matrix X of the algebraic Riccati */
/*             equation as produced by SLICOT Library routines SB02MD or */
/*             SB02OD. Matrix X is assumed non-negative definite. */
/*             On exit, if DICO = 'D', FACT = 'D' or 'C', OUFACT(2) = 1, */
/*             and INFO = 0, the N-by-N upper triangular part of this */
/*             array contains the Cholesky factor of the given matrix X, */
/*             which is found to be positive definite. */
/*             On exit, if DICO = 'D', FACT = 'D' or 'C', OUFACT(2) = 2, */
/*             and INFO = 0, the leading N-by-N part of this array */
/*             contains the matrix of orthonormal eigenvectors of X. */
/*             On exit X is unchanged if DICO = 'C' or FACT = 'N'. */

/*     LDX     INTEGER */
/*             The leading dimension of array X.  LDX >= MAX(1,N). */

/*     RNORM   (input) DOUBLE PRECISION */
/*             If FACT = 'U', this parameter must contain the 1-norm of */
/*             the original matrix R (before factoring it). */
/*             Otherwise, this parameter is not used. */

/*     F       (output) DOUBLE PRECISION array, dimension (LDF,N) */
/*             The leading M-by-N part of this array contains the */
/*             optimal feedback matrix F. */

/*     LDF     INTEGER */
/*             The leading dimension of array F.  LDF >= MAX(1,M). */

/*     OUFACT  (output) INTEGER array, dimension (2) */
/*             Information about the factorization finally used. */
/*             OUFACT(1) = 1:  Cholesky factorization of R (or R + B'XB) */
/*                             has been used; */
/*             OUFACT(1) = 2:  UdU' (if UPLO = 'U') or LdL' (if UPLO = */
/*                             'L') factorization of R (or R + B'XB) */
/*                             has been used; */
/*             OUFACT(2) = 1:  Cholesky factorization of X has been used; */
/*             OUFACT(2) = 2:  Spectral factorization of X has been used. */
/*             The value of OUFACT(2) is not set for DICO = 'C' or for */
/*             DICO = 'D' and FACT = 'N'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK, and DWORK(2) contains the reciprocal condition */
/*             number of the matrix R (for DICO = 'C') or of R + B'XB */
/*             (for DICO = 'D'). */
/*             If on exit INFO = 0, and OUFACT(2) = 2, then DWORK(3),..., */
/*             DWORK(N+2) contain the eigenvalues of X, in ascending */
/*             order. */

/*     LDWORK  INTEGER */
/*             Dimension of working array DWORK. */
/*             LDWORK >= max(2,3*M)         if FACT = 'N'; */
/*             LDWORK >= max(2,2*M)         if FACT = 'U'; */
/*             LDWORK >= max(2,3*M)         if FACT = 'C', DICO = 'C'; */
/*             LDWORK >= N+3*M+2            if FACT = 'C', DICO = 'D'; */
/*             LDWORK >= max(2,min(P,M)+M)  if FACT = 'D', DICO = 'C'; */
/*             LDWORK >= max(N+3*M+2,4*N+1) if FACT = 'D', DICO = 'D'. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = i:  if the i-th element of the d factor is exactly zero; */
/*                   the UdU' (or LdL') factorization has been completed, */
/*                   but the block diagonal matrix d is exactly singular; */
/*             = M+1:  if the matrix R (if DICO = 'C'), or R + B'XB */
/*                   (if DICO = 'D') is numerically singular (to working */
/*                   precision); */
/*             = M+2:  if one or more of the eigenvalues of X has not */
/*                   converged. */

/*     METHOD */

/*     The optimal feedback matrix F is obtained as the solution to the */
/*     system of linear equations */

/*        (R + B'XB) * F = B'XA + L' */

/*     in the discrete-time case and */

/*        R * F = B'X + L' */

/*     in the continuous-time case, with R replaced by D'D if FACT = 'D'. */
/*     The factored form of R, specified by FACT <> 'N', is taken into */
/*     account. If FACT = 'N', Cholesky factorization is tried first, but */
/*     if the coefficient matrix is not positive definite, then UdU' (or */
/*     LdL') factorization is used. The discrete-time case involves */
/*     updating of a triangular factorization of R (or D'D); Cholesky or */
/*     symmetric spectral factorization of X is employed to avoid */
/*     squaring of the condition number of the matrix. When D is given, */
/*     its QR factorization is determined, and the triangular factor is */
/*     used as described above. */

/*     NUMERICAL ASPECTS */

/*     The algorithm consists of numerically stable steps. */
/*                                    3     2 */
/*     For DICO = 'C', it requires O(m  + mn ) floating point operations */
/*                           2 */
/*     if FACT = 'N' and O(mn ) floating point operations, otherwise. */
/*     For DICO = 'D', the operation counts are similar, but additional */
/*        3 */
/*     O(n ) floating point operations may be needed in the worst case. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Sep. 1997. */
/*     Supersedes Release 2.0 routine SB02BD by M. Vanbegin, and */
/*     P. Van Dooren, Philips Research Laboratory, Brussels, Belgium. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Algebraic Riccati equation, closed loop system, continuous-time */
/*     system, discrete-time system, matrix algebra, optimal control, */
/*     optimal regulator. */

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

#line 333 "SB02ND.f"
    /* Parameter adjustments */
#line 333 "SB02ND.f"
    a_dim1 = *lda;
#line 333 "SB02ND.f"
    a_offset = 1 + a_dim1;
#line 333 "SB02ND.f"
    a -= a_offset;
#line 333 "SB02ND.f"
    b_dim1 = *ldb;
#line 333 "SB02ND.f"
    b_offset = 1 + b_dim1;
#line 333 "SB02ND.f"
    b -= b_offset;
#line 333 "SB02ND.f"
    r_dim1 = *ldr;
#line 333 "SB02ND.f"
    r_offset = 1 + r_dim1;
#line 333 "SB02ND.f"
    r__ -= r_offset;
#line 333 "SB02ND.f"
    --ipiv;
#line 333 "SB02ND.f"
    l_dim1 = *ldl;
#line 333 "SB02ND.f"
    l_offset = 1 + l_dim1;
#line 333 "SB02ND.f"
    l -= l_offset;
#line 333 "SB02ND.f"
    x_dim1 = *ldx;
#line 333 "SB02ND.f"
    x_offset = 1 + x_dim1;
#line 333 "SB02ND.f"
    x -= x_offset;
#line 333 "SB02ND.f"
    f_dim1 = *ldf;
#line 333 "SB02ND.f"
    f_offset = 1 + f_dim1;
#line 333 "SB02ND.f"
    f -= f_offset;
#line 333 "SB02ND.f"
    --oufact;
#line 333 "SB02ND.f"
    --iwork;
#line 333 "SB02ND.f"
    --dwork;
#line 333 "SB02ND.f"

#line 333 "SB02ND.f"
    /* Function Body */
#line 333 "SB02ND.f"
    *info = 0;
#line 334 "SB02ND.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 335 "SB02ND.f"
    lfactc = lsame_(fact, "C", (ftnlen)1, (ftnlen)1);
#line 336 "SB02ND.f"
    lfactd = lsame_(fact, "D", (ftnlen)1, (ftnlen)1);
#line 337 "SB02ND.f"
    lfactu = lsame_(fact, "U", (ftnlen)1, (ftnlen)1);
#line 338 "SB02ND.f"
    luplou = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 339 "SB02ND.f"
    withl = lsame_(jobl, "N", (ftnlen)1, (ftnlen)1);
#line 340 "SB02ND.f"
    lfacta = lfactc || lfactd || lfactu;

/*     Test the input scalar arguments. */

#line 344 "SB02ND.f"
    if (! discr && ! lsame_(dico, "C", (ftnlen)1, (ftnlen)1)) {
#line 345 "SB02ND.f"
	*info = -1;
#line 346 "SB02ND.f"
    } else if (! lfacta && ! lsame_(fact, "N", (ftnlen)1, (ftnlen)1) || discr 
	    && lfactu) {
#line 348 "SB02ND.f"
	*info = -2;
#line 349 "SB02ND.f"
    } else if (! luplou && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 350 "SB02ND.f"
	*info = -3;
#line 351 "SB02ND.f"
    } else if (! withl && ! lsame_(jobl, "Z", (ftnlen)1, (ftnlen)1)) {
#line 352 "SB02ND.f"
	*info = -4;
#line 353 "SB02ND.f"
    } else if (*n < 0) {
#line 354 "SB02ND.f"
	*info = -5;
#line 355 "SB02ND.f"
    } else if (*m < 0) {
#line 356 "SB02ND.f"
	*info = -6;
#line 357 "SB02ND.f"
    } else if (*p < 0) {
#line 358 "SB02ND.f"
	*info = -7;
#line 359 "SB02ND.f"
    } else if (! discr && *lda < 1 || discr && *lda < max(1,*n)) {
#line 361 "SB02ND.f"
	*info = -9;
#line 362 "SB02ND.f"
    } else if (*ldb < max(1,*n)) {
#line 363 "SB02ND.f"
	*info = -11;
#line 364 "SB02ND.f"
    } else if (*ldr < max(1,*m) || lfactd && *ldr < max(1,*p)) {
#line 366 "SB02ND.f"
	*info = -13;
#line 367 "SB02ND.f"
    } else if (! withl && *ldl < 1 || withl && *ldl < max(1,*n)) {
#line 369 "SB02ND.f"
	*info = -16;
#line 370 "SB02ND.f"
    } else if (*ldx < max(1,*n)) {
#line 371 "SB02ND.f"
	*info = -18;
#line 372 "SB02ND.f"
    } else if (lfactu) {
#line 373 "SB02ND.f"
	if (*rnorm < 0.) {
#line 373 "SB02ND.f"
	    *info = -19;
#line 373 "SB02ND.f"
	}
#line 375 "SB02ND.f"
    }
#line 376 "SB02ND.f"
    if (*ldf < max(1,*m)) {
#line 377 "SB02ND.f"
	*info = -21;
#line 378 "SB02ND.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 378 "SB02ND.f"
	i__1 = 2, i__2 = *m * 3;
/* Computing MAX */
#line 378 "SB02ND.f"
	i__3 = 2, i__4 = *m << 1;
/* Computing MAX */
#line 378 "SB02ND.f"
	i__5 = 2, i__6 = min(*p,*m) + *m;
/* Computing MAX */
#line 378 "SB02ND.f"
	i__7 = *n + *m * 3 + 2, i__8 = (*n << 2) + 1;
#line 378 "SB02ND.f"
	if ((! lfacta || lfactc && ! discr) && *ldwork < max(i__1,i__2) || 
		lfactu && *ldwork < max(i__3,i__4) || discr && lfactc && *
		ldwork < *n + *m * 3 + 2 || ! discr && lfactd && *ldwork < 
		max(i__5,i__6) || discr && lfactd && *ldwork < max(i__7,i__8))
		 {
#line 386 "SB02ND.f"
	    *info = -25;
#line 387 "SB02ND.f"
	}
#line 387 "SB02ND.f"
    }

#line 389 "SB02ND.f"
    if (*info != 0) {

/*        Error return. */

#line 393 "SB02ND.f"
	i__1 = -(*info);
#line 393 "SB02ND.f"
	xerbla_("SB02ND", &i__1, (ftnlen)6);
#line 394 "SB02ND.f"
	return 0;
#line 395 "SB02ND.f"
    }

/*     Quick return if possible. */

#line 399 "SB02ND.f"
    if (*n == 0 || *m == 0 || lfactd && *p == 0) {
#line 400 "SB02ND.f"
	dwork[1] = 1.;
#line 401 "SB02ND.f"
	dwork[2] = 1.;
#line 402 "SB02ND.f"
	return 0;
#line 403 "SB02ND.f"
    }

#line 405 "SB02ND.f"
    wrkopt = 1;
#line 406 "SB02ND.f"
    eps = dlamch_("Epsilon", (ftnlen)7);

/*     Determine the right-hand side of the matrix equation. */
/*     Compute  B'X  in F. */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 417 "SB02ND.f"
    dgemm_("Transpose", "No transpose", m, n, n, &c_b16, &b[b_offset], ldb, &
	    x[x_offset], ldx, &c_b17, &f[f_offset], ldf, (ftnlen)9, (ftnlen)
	    12);

#line 420 "SB02ND.f"
    if (! lfacta) {
#line 421 "SB02ND.f"
	if (discr) {

/*           Discrete-time case with R not factored. Compute R + B'XB. */

#line 425 "SB02ND.f"
	    if (luplou) {

#line 427 "SB02ND.f"
		i__1 = *m;
#line 427 "SB02ND.f"
		for (j = 1; j <= i__1; ++j) {
#line 428 "SB02ND.f"
		    dgemv_("No transpose", &j, n, &c_b16, &f[f_offset], ldf, &
			    b[j * b_dim1 + 1], &c__1, &c_b16, &r__[j * r_dim1 
			    + 1], &c__1, (ftnlen)12);
#line 430 "SB02ND.f"
/* L10: */
#line 430 "SB02ND.f"
		}

#line 432 "SB02ND.f"
	    } else {

#line 434 "SB02ND.f"
		i__1 = *m;
#line 434 "SB02ND.f"
		for (j = 1; j <= i__1; ++j) {
#line 435 "SB02ND.f"
		    dgemv_("Transpose", n, &j, &c_b16, &b[b_offset], ldb, &f[
			    j + f_dim1], ldf, &c_b16, &r__[j + r_dim1], ldr, (
			    ftnlen)9);
#line 437 "SB02ND.f"
/* L20: */
#line 437 "SB02ND.f"
		}

#line 439 "SB02ND.f"
	    }
#line 440 "SB02ND.f"
	}

/*        Compute the 1-norm of the matrix  R  or  R + B'XB. */
/*        Workspace: need M. */

#line 445 "SB02ND.f"
	rnormp = dlansy_("1-norm", uplo, m, &r__[r_offset], ldr, &dwork[1], (
		ftnlen)6, (ftnlen)1);
#line 446 "SB02ND.f"
	wrkopt = max(wrkopt,*m);
#line 447 "SB02ND.f"
    }

#line 449 "SB02ND.f"
    if (discr) {

/*        For discrete-time case, postmultiply B'X by A. */
/*        Workspace: need N. */

#line 454 "SB02ND.f"
	i__1 = *m;
#line 454 "SB02ND.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 455 "SB02ND.f"
	    dcopy_(n, &f[i__ + f_dim1], ldf, &dwork[1], &c__1);
#line 456 "SB02ND.f"
	    dgemv_("Transpose", n, n, &c_b16, &a[a_offset], lda, &dwork[1], &
		    c__1, &c_b17, &f[i__ + f_dim1], ldf, (ftnlen)9);
#line 458 "SB02ND.f"
/* L30: */
#line 458 "SB02ND.f"
	}

#line 460 "SB02ND.f"
	wrkopt = max(wrkopt,*n);
#line 461 "SB02ND.f"
    }

#line 463 "SB02ND.f"
    if (withl) {

/*        Add L'. */

#line 467 "SB02ND.f"
	i__1 = *m;
#line 467 "SB02ND.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 469 "SB02ND.f"
	    i__2 = *n;
#line 469 "SB02ND.f"
	    for (j = 1; j <= i__2; ++j) {
#line 470 "SB02ND.f"
		f[i__ + j * f_dim1] += l[j + i__ * l_dim1];
#line 471 "SB02ND.f"
/* L40: */
#line 471 "SB02ND.f"
	    }

#line 473 "SB02ND.f"
/* L50: */
#line 473 "SB02ND.f"
	}

#line 475 "SB02ND.f"
    }

/*     Solve the matrix equation. */

#line 479 "SB02ND.f"
    if (lfacta) {

/*        Case 1: Matrix R is given in a factored form. */

#line 483 "SB02ND.f"
	if (lfactd) {

/*           Use QR factorization of D. */
/*           Workspace: need   min(P,M) + M, */
/*                      prefer min(P,M) + M*NB. */

#line 489 "SB02ND.f"
	    itau = 1;
#line 490 "SB02ND.f"
	    jwork = itau + min(*p,*m);
#line 491 "SB02ND.f"
	    i__1 = *ldwork - jwork + 1;
#line 491 "SB02ND.f"
	    dgeqrf_(p, m, &r__[r_offset], ldr, &dwork[itau], &dwork[jwork], &
		    i__1, &ifail);
/* Computing MAX */
#line 493 "SB02ND.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 493 "SB02ND.f"
	    wrkopt = max(i__1,i__2);

/*           Make positive the diagonal elements of the triangular */
/*           factor. Construct the strictly lower triangle, if requested. */

#line 498 "SB02ND.f"
	    i__1 = *m;
#line 498 "SB02ND.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 499 "SB02ND.f"
		if (r__[i__ + i__ * r_dim1] < 0.) {

#line 501 "SB02ND.f"
		    i__2 = *m;
#line 501 "SB02ND.f"
		    for (j = i__; j <= i__2; ++j) {
#line 502 "SB02ND.f"
			r__[i__ + j * r_dim1] = -r__[i__ + j * r_dim1];
#line 503 "SB02ND.f"
/* L60: */
#line 503 "SB02ND.f"
		    }

#line 505 "SB02ND.f"
		}
#line 506 "SB02ND.f"
		if (! luplou) {
#line 506 "SB02ND.f"
		    i__2 = i__ - 1;
#line 506 "SB02ND.f"
		    dcopy_(&i__2, &r__[i__ * r_dim1 + 1], &c__1, &r__[i__ + 
			    r_dim1], ldr);
#line 506 "SB02ND.f"
		}
#line 508 "SB02ND.f"
/* L70: */
#line 508 "SB02ND.f"
	    }

#line 510 "SB02ND.f"
	    if (*p < *m) {
#line 511 "SB02ND.f"
		i__1 = *m - *p;
#line 511 "SB02ND.f"
		dlaset_("Full", &i__1, m, &c_b17, &c_b17, &r__[*p + 1 + 
			r_dim1], ldr, (ftnlen)4);
#line 512 "SB02ND.f"
		if (! discr) {
#line 513 "SB02ND.f"
		    dwork[2] = 0.;
#line 514 "SB02ND.f"
		    *info = *m + 1;
#line 515 "SB02ND.f"
		    return 0;
#line 516 "SB02ND.f"
		}
#line 517 "SB02ND.f"
	    }
#line 518 "SB02ND.f"
	}

#line 520 "SB02ND.f"
	jw = 1;
#line 521 "SB02ND.f"
	if (discr) {

/*           Discrete-time case. Update the factorization for B'XB. */
/*           Try first the Cholesky factorization of X, saving the */
/*           diagonal of X, in order to recover it, if X is not positive */
/*           definite. In the later case, use spectral factorization. */
/*           Workspace: need N. */
/*           Define     JW = 1   for Cholesky factorization of X, */
/*                      JW = N+3 for spectral factorization of X. */

#line 531 "SB02ND.f"
	    i__1 = *ldx + 1;
#line 531 "SB02ND.f"
	    dcopy_(n, &x[x_offset], &i__1, &dwork[1], &c__1);
#line 532 "SB02ND.f"
	    dpotrf_("Upper", n, &x[x_offset], ldx, &ifail, (ftnlen)5);
#line 533 "SB02ND.f"
	    if (ifail == 0) {

/*              Use Cholesky factorization of X to compute chol(X)*B. */

#line 537 "SB02ND.f"
		oufact[2] = 1;
#line 538 "SB02ND.f"
		dtrmm_("Left", "Upper", "No transpose", "Non unit", n, m, &
			c_b16, &x[x_offset], ldx, &b[b_offset], ldb, (ftnlen)
			4, (ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 540 "SB02ND.f"
	    } else {

/*              Use spectral factorization of X, X = UVU'. */
/*              Workspace: need   4*N+1, */
/*                         prefer N*(NB+2)+N+2. */

#line 546 "SB02ND.f"
		jw = *n + 3;
#line 547 "SB02ND.f"
		oufact[2] = 2;
#line 548 "SB02ND.f"
		i__1 = *ldx + 1;
#line 548 "SB02ND.f"
		dcopy_(n, &dwork[1], &c__1, &x[x_offset], &i__1);
#line 549 "SB02ND.f"
		i__1 = *ldwork - jw + 1;
#line 549 "SB02ND.f"
		dsyev_("Vectors", "Lower", n, &x[x_offset], ldx, &dwork[3], &
			dwork[jw], &i__1, &ifail, (ftnlen)7, (ftnlen)5);
#line 551 "SB02ND.f"
		if (ifail > 0) {
#line 552 "SB02ND.f"
		    *info = *m + 2;
#line 553 "SB02ND.f"
		    return 0;
#line 554 "SB02ND.f"
		}
/* Computing MAX */
#line 555 "SB02ND.f"
		i__1 = wrkopt, i__2 = (integer) dwork[jw] + jw - 1;
#line 555 "SB02ND.f"
		wrkopt = max(i__1,i__2);
#line 556 "SB02ND.f"
		temp = (d__1 = dwork[*n + 2], abs(d__1)) * eps;

/*              Count the negligible eigenvalues and compute sqrt(V)U'B. */
/*              Workspace: need 2*N+2. */

#line 561 "SB02ND.f"
		jz = 0;

#line 563 "SB02ND.f"
L80:
#line 564 "SB02ND.f"
		if ((d__1 = dwork[jz + 3], abs(d__1)) <= temp) {
#line 565 "SB02ND.f"
		    ++jz;
#line 566 "SB02ND.f"
		    if (jz < *n) {
#line 566 "SB02ND.f"
			goto L80;
#line 566 "SB02ND.f"
		    }
#line 567 "SB02ND.f"
		}

#line 569 "SB02ND.f"
		i__1 = *m;
#line 569 "SB02ND.f"
		for (j = 1; j <= i__1; ++j) {
#line 570 "SB02ND.f"
		    dcopy_(n, &b[j * b_dim1 + 1], &c__1, &dwork[jw], &c__1);
#line 571 "SB02ND.f"
		    dgemv_("Transpose", n, n, &c_b16, &x[x_offset], ldx, &
			    dwork[jw], &c__1, &c_b17, &b[j * b_dim1 + 1], &
			    c__1, (ftnlen)9);
#line 573 "SB02ND.f"
/* L90: */
#line 573 "SB02ND.f"
		}

#line 575 "SB02ND.f"
		i__1 = *n;
#line 575 "SB02ND.f"
		for (i__ = jz + 1; i__ <= i__1; ++i__) {
#line 576 "SB02ND.f"
		    d__2 = sqrt((d__1 = dwork[i__ + 2], abs(d__1)));
#line 576 "SB02ND.f"
		    dscal_(m, &d__2, &b[i__ + b_dim1], ldb);
#line 578 "SB02ND.f"
/* L100: */
#line 578 "SB02ND.f"
		}

#line 580 "SB02ND.f"
		if (jz > 0) {
#line 580 "SB02ND.f"
		    dlaset_("Full", &jz, m, &c_b17, &c_b17, &b[b_offset], ldb,
			     (ftnlen)4);
#line 580 "SB02ND.f"
		}
#line 582 "SB02ND.f"
	    }

/*           Update the triangular factorization. */

#line 586 "SB02ND.f"
	    if (! luplou) {

/*              For efficiency, use the transposed of the lower triangle. */

#line 590 "SB02ND.f"
		i__1 = *m;
#line 590 "SB02ND.f"
		for (i__ = 2; i__ <= i__1; ++i__) {
#line 591 "SB02ND.f"
		    i__2 = i__ - 1;
#line 591 "SB02ND.f"
		    dcopy_(&i__2, &r__[i__ + r_dim1], ldr, &r__[i__ * r_dim1 
			    + 1], &c__1);
#line 592 "SB02ND.f"
/* L110: */
#line 592 "SB02ND.f"
		}

#line 594 "SB02ND.f"
	    }

/*           Workspace: need JW+2*M-1. */

#line 598 "SB02ND.f"
	    mb04kd_("Full", m, &c__0, n, &r__[r_offset], ldr, &b[b_offset], 
		    ldb, dummy, n, dummy, m, &dwork[jw], &dwork[jw + *n], (
		    ftnlen)4);
/* Computing MAX */
#line 600 "SB02ND.f"
	    i__1 = wrkopt, i__2 = jw + (*m << 1) - 1;
#line 600 "SB02ND.f"
	    wrkopt = max(i__1,i__2);

/*           Make positive the diagonal elements of the triangular */
/*           factor. */

#line 605 "SB02ND.f"
	    i__1 = *m;
#line 605 "SB02ND.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 606 "SB02ND.f"
		if (r__[i__ + i__ * r_dim1] < 0.) {

#line 608 "SB02ND.f"
		    i__2 = *m;
#line 608 "SB02ND.f"
		    for (j = i__; j <= i__2; ++j) {
#line 609 "SB02ND.f"
			r__[i__ + j * r_dim1] = -r__[i__ + j * r_dim1];
#line 610 "SB02ND.f"
/* L120: */
#line 610 "SB02ND.f"
		    }

#line 612 "SB02ND.f"
		}
#line 613 "SB02ND.f"
/* L130: */
#line 613 "SB02ND.f"
	    }

#line 615 "SB02ND.f"
	    if (! luplou) {

/*              Construct the lower triangle. */

#line 619 "SB02ND.f"
		i__1 = *m;
#line 619 "SB02ND.f"
		for (i__ = 2; i__ <= i__1; ++i__) {
#line 620 "SB02ND.f"
		    i__2 = i__ - 1;
#line 620 "SB02ND.f"
		    dcopy_(&i__2, &r__[i__ * r_dim1 + 1], &c__1, &r__[i__ + 
			    r_dim1], ldr);
#line 621 "SB02ND.f"
/* L140: */
#line 621 "SB02ND.f"
		}

#line 623 "SB02ND.f"
	    }
#line 624 "SB02ND.f"
	}

/*        Compute the condition number of the coefficient matrix. */

#line 628 "SB02ND.f"
	if (! lfactu) {

/*           Workspace: need JW+3*M-1. */

#line 632 "SB02ND.f"
	    dtrcon_("1-norm", uplo, "Non unit", m, &r__[r_offset], ldr, &
		    rcond, &dwork[jw], &iwork[1], &ifail, (ftnlen)6, (ftnlen)
		    1, (ftnlen)8);
#line 634 "SB02ND.f"
	    oufact[1] = 1;
/* Computing MAX */
#line 635 "SB02ND.f"
	    i__1 = wrkopt, i__2 = jw + *m * 3 - 1;
#line 635 "SB02ND.f"
	    wrkopt = max(i__1,i__2);
#line 636 "SB02ND.f"
	} else {

/*           Workspace: need 2*M. */

#line 640 "SB02ND.f"
	    dsycon_(uplo, m, &r__[r_offset], ldr, &ipiv[1], rnorm, &rcond, &
		    dwork[1], &iwork[1], info, (ftnlen)1);
#line 642 "SB02ND.f"
	    oufact[1] = 2;
/* Computing MAX */
#line 643 "SB02ND.f"
	    i__1 = wrkopt, i__2 = *m << 1;
#line 643 "SB02ND.f"
	    wrkopt = max(i__1,i__2);
#line 644 "SB02ND.f"
	}
#line 645 "SB02ND.f"
	dwork[2] = rcond;
#line 646 "SB02ND.f"
	if (rcond < eps) {
#line 647 "SB02ND.f"
	    *info = *m + 1;
#line 648 "SB02ND.f"
	    return 0;
#line 649 "SB02ND.f"
	}

#line 651 "SB02ND.f"
    } else {

/*        Case 2: Matrix R is given in an unfactored form. */

/*        Save the given triangle of  R  or  R + B'XB  in the other */
/*        strict triangle and the diagonal in the workspace, and try */
/*        Cholesky factorization. */
/*        Workspace: need M. */

#line 660 "SB02ND.f"
	i__1 = *ldr + 1;
#line 660 "SB02ND.f"
	dcopy_(m, &r__[r_offset], &i__1, &dwork[1], &c__1);
#line 661 "SB02ND.f"
	if (luplou) {

#line 663 "SB02ND.f"
	    i__1 = *m;
#line 663 "SB02ND.f"
	    for (j = 2; j <= i__1; ++j) {
#line 664 "SB02ND.f"
		i__2 = j - 1;
#line 664 "SB02ND.f"
		dcopy_(&i__2, &r__[j * r_dim1 + 1], &c__1, &r__[j + r_dim1], 
			ldr);
#line 665 "SB02ND.f"
/* L150: */
#line 665 "SB02ND.f"
	    }

#line 667 "SB02ND.f"
	} else {

#line 669 "SB02ND.f"
	    i__1 = *m;
#line 669 "SB02ND.f"
	    for (j = 2; j <= i__1; ++j) {
#line 670 "SB02ND.f"
		i__2 = j - 1;
#line 670 "SB02ND.f"
		dcopy_(&i__2, &r__[j + r_dim1], ldr, &r__[j * r_dim1 + 1], &
			c__1);
#line 671 "SB02ND.f"
/* L160: */
#line 671 "SB02ND.f"
	    }

#line 673 "SB02ND.f"
	}
#line 674 "SB02ND.f"
	dpotrf_(uplo, m, &r__[r_offset], ldr, info, (ftnlen)1);
#line 675 "SB02ND.f"
	oufact[1] = 1;
#line 676 "SB02ND.f"
	if (*info == 0) {

/*           Compute the reciprocal of the condition number of R. */
/*           Workspace: need 3*M. */

#line 681 "SB02ND.f"
	    dpocon_(uplo, m, &r__[r_offset], ldr, &rnormp, &rcond, &dwork[1], 
		    &iwork[1], info, (ftnlen)1);

/*           Return if the matrix is singular to working precision. */

#line 686 "SB02ND.f"
	    dwork[2] = rcond;
#line 687 "SB02ND.f"
	    if (rcond < eps) {
#line 688 "SB02ND.f"
		*info = *m + 1;
#line 689 "SB02ND.f"
		return 0;
#line 690 "SB02ND.f"
	    }
/* Computing MAX */
#line 691 "SB02ND.f"
	    i__1 = wrkopt, i__2 = *m * 3;
#line 691 "SB02ND.f"
	    wrkopt = max(i__1,i__2);
#line 692 "SB02ND.f"
	} else {

/*           Use UdU' or LdL' factorization, first restoring the saved */
/*           triangle. */

#line 697 "SB02ND.f"
	    i__1 = *ldr + 1;
#line 697 "SB02ND.f"
	    dcopy_(m, &dwork[1], &c__1, &r__[r_offset], &i__1);
#line 698 "SB02ND.f"
	    if (luplou) {

#line 700 "SB02ND.f"
		i__1 = *m;
#line 700 "SB02ND.f"
		for (j = 2; j <= i__1; ++j) {
#line 701 "SB02ND.f"
		    i__2 = j - 1;
#line 701 "SB02ND.f"
		    dcopy_(&i__2, &r__[j + r_dim1], ldr, &r__[j * r_dim1 + 1],
			     &c__1);
#line 702 "SB02ND.f"
/* L170: */
#line 702 "SB02ND.f"
		}

#line 704 "SB02ND.f"
	    } else {

#line 706 "SB02ND.f"
		i__1 = *m;
#line 706 "SB02ND.f"
		for (j = 2; j <= i__1; ++j) {
#line 707 "SB02ND.f"
		    i__2 = j - 1;
#line 707 "SB02ND.f"
		    dcopy_(&i__2, &r__[j * r_dim1 + 1], &c__1, &r__[j + 
			    r_dim1], ldr);
#line 708 "SB02ND.f"
/* L180: */
#line 708 "SB02ND.f"
		}

#line 710 "SB02ND.f"
	    }

/*           Workspace: need   1, */
/*                      prefer M*NB. */

#line 715 "SB02ND.f"
	    dsytrf_(uplo, m, &r__[r_offset], ldr, &ipiv[1], &dwork[1], ldwork,
		     info, (ftnlen)1);
#line 716 "SB02ND.f"
	    oufact[1] = 2;
#line 717 "SB02ND.f"
	    if (*info > 0) {
#line 717 "SB02ND.f"
		return 0;
#line 717 "SB02ND.f"
	    }
/* Computing MAX */
#line 719 "SB02ND.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 719 "SB02ND.f"
	    wrkopt = max(i__1,i__2);

/*           Compute the reciprocal of the condition number of R. */
/*           Workspace: need   2*M. */

#line 724 "SB02ND.f"
	    dsycon_(uplo, m, &r__[r_offset], ldr, &ipiv[1], &rnormp, &rcond, &
		    dwork[1], &iwork[1], info, (ftnlen)1);

/*           Return if the matrix is singular to working precision. */

#line 729 "SB02ND.f"
	    dwork[2] = rcond;
#line 730 "SB02ND.f"
	    if (rcond < eps) {
#line 731 "SB02ND.f"
		*info = *m + 1;
#line 732 "SB02ND.f"
		return 0;
#line 733 "SB02ND.f"
	    }
#line 734 "SB02ND.f"
	}
#line 735 "SB02ND.f"
    }

#line 737 "SB02ND.f"
    if (oufact[1] == 1) {

/*        Solve the positive definite linear system. */

#line 741 "SB02ND.f"
	dpotrs_(uplo, m, n, &r__[r_offset], ldr, &f[f_offset], ldf, info, (
		ftnlen)1);
#line 742 "SB02ND.f"
    } else {

/*        Solve the indefinite linear system. */

#line 746 "SB02ND.f"
	dsytrs_(uplo, m, n, &r__[r_offset], ldr, &ipiv[1], &f[f_offset], ldf, 
		info, (ftnlen)1);
#line 747 "SB02ND.f"
    }

/*     Set the optimal workspace. */

#line 751 "SB02ND.f"
    dwork[1] = (doublereal) wrkopt;

#line 753 "SB02ND.f"
    return 0;
/* *** Last line of SB02ND *** */
} /* sb02nd_ */

