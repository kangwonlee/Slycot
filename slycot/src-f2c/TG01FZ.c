#line 1 "TG01FZ.f"
/* TG01FZ.f -- translated by f2c (version 20100827).
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

#line 1 "TG01FZ.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;

/* Subroutine */ int tg01fz_(char *compq, char *compz, char *joba, integer *l,
	 integer *n, integer *m, integer *p, doublecomplex *a, integer *lda, 
	doublecomplex *e, integer *lde, doublecomplex *b, integer *ldb, 
	doublecomplex *c__, integer *ldc, doublecomplex *q, integer *ldq, 
	doublecomplex *z__, integer *ldz, integer *ranke, integer *rnka22, 
	doublereal *tol, integer *iwork, doublereal *dwork, doublecomplex *
	zwork, integer *lzwork, integer *info, ftnlen compq_len, ftnlen 
	compz_len, ftnlen joba_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, e_dim1, 
	    e_offset, q_dim1, q_offset, z_dim1, z_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, nb, lh, ln, kw, ir1, la22, na22;
    static logical ilq, ilz;
    static integer lwr, ire1;
    static logical reda;
    static doublereal sval[3];
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical withb, withc, redtr;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), mb3oyz_(integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, doublecomplex *, doublereal *, 
	    doublecomplex *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal toldef;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static integer icompq, icompz;
    extern /* Subroutine */ int zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen);
    static doublereal svlmax;
    static logical lquery;
    static integer wrkopt;
    extern /* Subroutine */ int zunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), zunmrz_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen), ztzrzf_(integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, integer *)
	    ;


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

/*     To compute for the descriptor system (A-lambda E,B,C) */
/*     the unitary transformation matrices Q and Z such that the */
/*     transformed system (Q'*A*Z-lambda Q'*E*Z, Q'*B, C*Z) is */
/*     in a SVD-like coordinate form with */

/*                  ( A11  A12 )             ( Er  0 ) */
/*         Q'*A*Z = (          ) ,  Q'*E*Z = (       ) , */
/*                  ( A21  A22 )             (  0  0 ) */

/*     where Er is an upper triangular invertible matrix, and ' denotes */
/*     the conjugate transpose. Optionally, the A22 matrix can be further */
/*     reduced to the form */

/*                  ( Ar  X ) */
/*            A22 = (       ) , */
/*                  (  0  0 ) */

/*     with Ar an upper triangular invertible matrix, and X either a full */
/*     or a zero matrix. */
/*     The left and/or right unitary transformations performed */
/*     to reduce E and A22 can be optionally accumulated. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     COMPQ   CHARACTER*1 */
/*             = 'N':  do not compute Q; */
/*             = 'I':  Q is initialized to the unit matrix, and the */
/*                     unitary matrix Q is returned; */
/*             = 'U':  Q must contain a unitary matrix Q1 on entry, */
/*                     and the product Q1*Q is returned. */

/*     COMPZ   CHARACTER*1 */
/*             = 'N':  do not compute Z; */
/*             = 'I':  Z is initialized to the unit matrix, and the */
/*                     unitary matrix Z is returned; */
/*             = 'U':  Z must contain a unitary matrix Z1 on entry, */
/*                     and the product Z1*Z is returned. */

/*     JOBA    CHARACTER*1 */
/*             = 'N':  do not reduce A22; */
/*             = 'R':  reduce A22 to a SVD-like upper triangular form. */
/*             = 'T':  reduce A22 to an upper trapezoidal form. */

/*     Input/Output Parameters */

/*     L       (input) INTEGER */
/*             The number of rows of matrices A, B, and E.  L >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of matrices A, E, and C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of matrix B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of rows of matrix C.  P >= 0. */

/*     A       (input/output) COMPLEX*16 array, dimension (LDA,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the state dynamics matrix A. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the transformed matrix Q'*A*Z. If JOBA = 'T', this matrix */
/*             is in the form */

/*                           ( A11  *   *  ) */
/*                  Q'*A*Z = (  *   Ar  X  ) , */
/*                           (  *   0   0  ) */

/*             where A11 is a RANKE-by-RANKE matrix and Ar is a */
/*             RNKA22-by-RNKA22 invertible upper triangular matrix. */
/*             If JOBA = 'R' then A has the above form with X = 0. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,L). */

/*     E       (input/output) COMPLEX*16 array, dimension (LDE,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the descriptor matrix E. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the transformed matrix Q'*E*Z. */

/*                      ( Er  0 ) */
/*             Q'*E*Z = (       ) , */
/*                      (  0  0 ) */

/*             where Er is a RANKE-by-RANKE upper triangular invertible */
/*             matrix. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,L). */

/*     B       (input/output) COMPLEX*16 array, dimension (LDB,M) */
/*             On entry, the leading L-by-M part of this array must */
/*             contain the input/state matrix B. */
/*             On exit, the leading L-by-M part of this array contains */
/*             the transformed matrix Q'*B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B. */
/*             LDB >= MAX(1,L) if M > 0 or LDB >= 1 if M = 0. */

/*     C       (input/output) COMPLEX*16 array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the transformed matrix C*Z. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     Q       (input/output) COMPLEX*16 array, dimension (LDQ,L) */
/*             If COMPQ = 'N':  Q is not referenced. */
/*             If COMPQ = 'I':  on entry, Q need not be set; */
/*                              on exit, the leading L-by-L part of this */
/*                              array contains the unitary matrix Q, */
/*                              where Q' is the product of Householder */
/*                              transformations which are applied to A, */
/*                              E, and B on the left. */
/*             If COMPQ = 'U':  on entry, the leading L-by-L part of this */
/*                              array must contain a unitary matrix Q1; */
/*                              on exit, the leading L-by-L part of this */
/*                              array contains the unitary matrix Q1*Q. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. */
/*             LDQ >= 1,        if COMPQ = 'N'; */
/*             LDQ >= MAX(1,L), if COMPQ = 'U' or 'I'. */

/*     Z       (input/output) COMPLEX*16 array, dimension (LDZ,N) */
/*             If COMPZ = 'N':  Z is not referenced. */
/*             If COMPZ = 'I':  on entry, Z need not be set; */
/*                              on exit, the leading N-by-N part of this */
/*                              array contains the unitary matrix Z, */
/*                              which is the product of Householder */
/*                              transformations applied to A, E, and C */
/*                              on the right. */
/*             If COMPZ = 'U':  on entry, the leading N-by-N part of this */
/*                              array must contain a unitary matrix Z1; */
/*                              on exit, the leading N-by-N part of this */
/*                              array contains the unitary matrix Z1*Z. */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z. */
/*             LDZ >= 1,        if COMPZ = 'N'; */
/*             LDZ >= MAX(1,N), if COMPZ = 'U' or 'I'. */

/*     RANKE   (output) INTEGER */
/*             The estimated rank of matrix E, and thus also the order */
/*             of the invertible upper triangular submatrix Er. */

/*     RNKA22  (output) INTEGER */
/*             If JOBA = 'R' or 'T', then RNKA22 is the estimated rank of */
/*             matrix A22, and thus also the order of the invertible */
/*             upper triangular submatrix Ar. */
/*             If JOBA = 'N', then RNKA22 is not referenced. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in determining the rank of E */
/*             and of A22. If the user sets TOL > 0, then the given */
/*             value of TOL is used as a lower bound for the */
/*             reciprocal condition numbers of leading submatrices */
/*             of R or R22 in the QR decompositions E * P = Q * R of E */
/*             or A22 * P22 = Q22 * R22 of A22. */
/*             A submatrix whose estimated condition number is less than */
/*             1/TOL is considered to be of full rank.  If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance, */
/*             defined by  TOLDEF = L*N*EPS,  is used instead, where */
/*             EPS is the machine precision (see LAPACK Library routine */
/*             DLAMCH). TOL < 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     DWORK   DOUBLE PRECISION array, dimension (2*N) */

/*     ZWORK   DOUBLE PRECISION array, dimension (LZWORK) */
/*             On exit, if INFO = 0, ZWORK(1) returns the optimal value */
/*             of LZWORK. */

/*     LZWORK  INTEGER */
/*             The length of the array ZWORK. */
/*             LZWORK >= MAX( 1, N+P, MIN(L,N)+MAX(3*N-1,M,L) ). */
/*             For optimal performance, LZWORK should be larger. */

/*             If LZWORK = -1, then a workspace query is assumed; */
/*             the routine only calculates the optimal size of the */
/*             ZWORK array, returns this value as the first entry of */
/*             the ZWORK array, and no error message related to LZWORK */
/*             is issued by XERBLA. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine computes a truncated QR factorization with column */
/*     pivoting of E, in the form */

/*                       ( E11 E12 ) */
/*           E * P = Q * (         ) */
/*                       (  0  E22 ) */

/*     and finds the largest RANKE-by-RANKE leading submatrix E11 whose */
/*     estimated condition number is less than 1/TOL. RANKE defines thus */
/*     the rank of matrix E. Further E22, being negligible, is set to */
/*     zero, and a unitary matrix Y is determined such that */

/*           ( E11 E12 ) = ( Er  0 ) * Y . */

/*     The overal transformation matrix Z results as Z = P * Y' and the */
/*     resulting transformed matrices Q'*A*Z and Q'*E*Z have the form */

/*                          ( Er  0 )                      ( A11  A12 ) */
/*         E <- Q'* E * Z = (       ) ,  A <- Q' * A * Z = (          ) , */
/*                          (  0  0 )                      ( A21  A22 ) */

/*     where Er is an upper triangular invertible matrix. */
/*     If JOBA = 'R' the same reduction is performed on A22 to obtain it */
/*     in the form */

/*                  ( Ar  0 ) */
/*            A22 = (       ) , */
/*                  (  0  0 ) */

/*     with Ar an upper triangular invertible matrix. */
/*     If JOBA = 'T' then A22 is row compressed using the QR */
/*     factorization with column pivoting to the form */

/*                  ( Ar  X ) */
/*            A22 = (       ) */
/*                  (  0  0 ) */

/*     with Ar an upper triangular invertible matrix. */

/*     The transformations are also applied to the rest of system */
/*     matrices */

/*          B <- Q' * B, C <- C * Z. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically backward stable and requires */
/*     0( L*L*N )  floating point operations. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     March 1999. */
/*     Complex version: V. Sima, Research Institute for Informatics, */
/*     Bucharest, Nov. 2008. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Descriptor system, matrix algebra, matrix operations, unitary */
/*     transformation. */

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

/*     Decode COMPQ. */

#line 335 "TG01FZ.f"
    /* Parameter adjustments */
#line 335 "TG01FZ.f"
    a_dim1 = *lda;
#line 335 "TG01FZ.f"
    a_offset = 1 + a_dim1;
#line 335 "TG01FZ.f"
    a -= a_offset;
#line 335 "TG01FZ.f"
    e_dim1 = *lde;
#line 335 "TG01FZ.f"
    e_offset = 1 + e_dim1;
#line 335 "TG01FZ.f"
    e -= e_offset;
#line 335 "TG01FZ.f"
    b_dim1 = *ldb;
#line 335 "TG01FZ.f"
    b_offset = 1 + b_dim1;
#line 335 "TG01FZ.f"
    b -= b_offset;
#line 335 "TG01FZ.f"
    c_dim1 = *ldc;
#line 335 "TG01FZ.f"
    c_offset = 1 + c_dim1;
#line 335 "TG01FZ.f"
    c__ -= c_offset;
#line 335 "TG01FZ.f"
    q_dim1 = *ldq;
#line 335 "TG01FZ.f"
    q_offset = 1 + q_dim1;
#line 335 "TG01FZ.f"
    q -= q_offset;
#line 335 "TG01FZ.f"
    z_dim1 = *ldz;
#line 335 "TG01FZ.f"
    z_offset = 1 + z_dim1;
#line 335 "TG01FZ.f"
    z__ -= z_offset;
#line 335 "TG01FZ.f"
    --iwork;
#line 335 "TG01FZ.f"
    --dwork;
#line 335 "TG01FZ.f"
    --zwork;
#line 335 "TG01FZ.f"

#line 335 "TG01FZ.f"
    /* Function Body */
#line 335 "TG01FZ.f"
    if (lsame_(compq, "N", (ftnlen)1, (ftnlen)1)) {
#line 336 "TG01FZ.f"
	ilq = FALSE_;
#line 337 "TG01FZ.f"
	icompq = 1;
#line 338 "TG01FZ.f"
    } else if (lsame_(compq, "U", (ftnlen)1, (ftnlen)1)) {
#line 339 "TG01FZ.f"
	ilq = TRUE_;
#line 340 "TG01FZ.f"
	icompq = 2;
#line 341 "TG01FZ.f"
    } else if (lsame_(compq, "I", (ftnlen)1, (ftnlen)1)) {
#line 342 "TG01FZ.f"
	ilq = TRUE_;
#line 343 "TG01FZ.f"
	icompq = 3;
#line 344 "TG01FZ.f"
    } else {
#line 345 "TG01FZ.f"
	icompq = 0;
#line 346 "TG01FZ.f"
    }

/*     Decode COMPZ. */

#line 350 "TG01FZ.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 351 "TG01FZ.f"
	ilz = FALSE_;
#line 352 "TG01FZ.f"
	icompz = 1;
#line 353 "TG01FZ.f"
    } else if (lsame_(compz, "U", (ftnlen)1, (ftnlen)1)) {
#line 354 "TG01FZ.f"
	ilz = TRUE_;
#line 355 "TG01FZ.f"
	icompz = 2;
#line 356 "TG01FZ.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 357 "TG01FZ.f"
	ilz = TRUE_;
#line 358 "TG01FZ.f"
	icompz = 3;
#line 359 "TG01FZ.f"
    } else {
#line 360 "TG01FZ.f"
	icompz = 0;
#line 361 "TG01FZ.f"
    }
#line 362 "TG01FZ.f"
    reda = lsame_(joba, "R", (ftnlen)1, (ftnlen)1);
#line 363 "TG01FZ.f"
    redtr = lsame_(joba, "T", (ftnlen)1, (ftnlen)1);
#line 364 "TG01FZ.f"
    withb = *m > 0;
#line 365 "TG01FZ.f"
    withc = *p > 0;
#line 366 "TG01FZ.f"
    lquery = *lzwork == -1;

/*     Test the input parameters. */

#line 370 "TG01FZ.f"
    ln = min(*l,*n);
#line 371 "TG01FZ.f"
    *info = 0;
/* Computing MAX */
/* Computing MAX */
#line 372 "TG01FZ.f"
    i__3 = *n * 3 - 1, i__3 = max(i__3,*m);
#line 372 "TG01FZ.f"
    i__1 = 1, i__2 = *n + *p, i__1 = max(i__1,i__2), i__2 = ln + max(i__3,*l);
#line 372 "TG01FZ.f"
    wrkopt = max(i__1,i__2);
#line 373 "TG01FZ.f"
    if (icompq <= 0) {
#line 374 "TG01FZ.f"
	*info = -1;
#line 375 "TG01FZ.f"
    } else if (icompz <= 0) {
#line 376 "TG01FZ.f"
	*info = -2;
#line 377 "TG01FZ.f"
    } else if (! lsame_(joba, "N", (ftnlen)1, (ftnlen)1) && ! reda && ! redtr)
	     {
#line 379 "TG01FZ.f"
	*info = -3;
#line 380 "TG01FZ.f"
    } else if (*l < 0) {
#line 381 "TG01FZ.f"
	*info = -4;
#line 382 "TG01FZ.f"
    } else if (*n < 0) {
#line 383 "TG01FZ.f"
	*info = -5;
#line 384 "TG01FZ.f"
    } else if (*m < 0) {
#line 385 "TG01FZ.f"
	*info = -6;
#line 386 "TG01FZ.f"
    } else if (*p < 0) {
#line 387 "TG01FZ.f"
	*info = -7;
#line 388 "TG01FZ.f"
    } else if (*lda < max(1,*l)) {
#line 389 "TG01FZ.f"
	*info = -9;
#line 390 "TG01FZ.f"
    } else if (*lde < max(1,*l)) {
#line 391 "TG01FZ.f"
	*info = -11;
#line 392 "TG01FZ.f"
    } else if (*ldb < 1 || withb && *ldb < *l) {
#line 393 "TG01FZ.f"
	*info = -13;
#line 394 "TG01FZ.f"
    } else if (*ldc < max(1,*p)) {
#line 395 "TG01FZ.f"
	*info = -15;
#line 396 "TG01FZ.f"
    } else if (ilq && *ldq < *l || *ldq < 1) {
#line 397 "TG01FZ.f"
	*info = -17;
#line 398 "TG01FZ.f"
    } else if (ilz && *ldz < *n || *ldz < 1) {
#line 399 "TG01FZ.f"
	*info = -19;
#line 400 "TG01FZ.f"
    } else if (*tol >= 1.) {
#line 401 "TG01FZ.f"
	*info = -22;
#line 402 "TG01FZ.f"
    } else {
#line 403 "TG01FZ.f"
	if (lquery) {
/* Computing MIN */
#line 404 "TG01FZ.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMQR", "LC", l, n, &ln, &c_n1,
		     (ftnlen)6, (ftnlen)2);
#line 404 "TG01FZ.f"
	    nb = min(i__1,i__2);
/* Computing MAX */
#line 405 "TG01FZ.f"
	    i__1 = wrkopt, i__2 = ln + *n * nb;
#line 405 "TG01FZ.f"
	    wrkopt = max(i__1,i__2);
#line 406 "TG01FZ.f"
	    if (withb) {
/* Computing MIN */
#line 407 "TG01FZ.f"
		i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMQR", "LC", l, m, &ln, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 407 "TG01FZ.f"
		nb = min(i__1,i__2);
/* Computing MAX */
#line 408 "TG01FZ.f"
		i__1 = wrkopt, i__2 = ln + *m * nb;
#line 408 "TG01FZ.f"
		wrkopt = max(i__1,i__2);
#line 409 "TG01FZ.f"
	    }
#line 410 "TG01FZ.f"
	    if (ilq) {
/* Computing MIN */
#line 411 "TG01FZ.f"
		i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMQR", "RN", l, l, &ln, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 411 "TG01FZ.f"
		nb = min(i__1,i__2);
/* Computing MAX */
#line 412 "TG01FZ.f"
		i__1 = wrkopt, i__2 = ln + *l * nb;
#line 412 "TG01FZ.f"
		wrkopt = max(i__1,i__2);
#line 413 "TG01FZ.f"
	    }
#line 414 "TG01FZ.f"
	    nb = ilaenv_(&c__1, "ZGERQF", " ", l, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
/* Computing MAX */
#line 415 "TG01FZ.f"
	    i__1 = wrkopt, i__2 = ln + *n * nb;
#line 415 "TG01FZ.f"
	    wrkopt = max(i__1,i__2);
/* Computing MIN */
#line 416 "TG01FZ.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMRQ", "RC", l, n, n, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 416 "TG01FZ.f"
	    nb = min(i__1,i__2);
/* Computing MAX */
#line 417 "TG01FZ.f"
	    i__1 = wrkopt, i__2 = *n + max(1,*l) * nb;
#line 417 "TG01FZ.f"
	    wrkopt = max(i__1,i__2);
#line 418 "TG01FZ.f"
	    if (withc) {
/* Computing MIN */
#line 419 "TG01FZ.f"
		i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMRQ", "RC", p, n, n, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 419 "TG01FZ.f"
		nb = min(i__1,i__2);
/* Computing MAX */
#line 420 "TG01FZ.f"
		i__1 = wrkopt, i__2 = *n + max(1,*p) * nb;
#line 420 "TG01FZ.f"
		wrkopt = max(i__1,i__2);
#line 421 "TG01FZ.f"
	    }
#line 422 "TG01FZ.f"
	    if (ilz) {
/* Computing MIN */
#line 423 "TG01FZ.f"
		i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMRQ", "RC", n, n, n, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 423 "TG01FZ.f"
		nb = min(i__1,i__2);
/* Computing MAX */
#line 424 "TG01FZ.f"
		i__1 = wrkopt, i__2 = *n + max(1,*n) * nb;
#line 424 "TG01FZ.f"
		wrkopt = max(i__1,i__2);
#line 425 "TG01FZ.f"
	    }
#line 426 "TG01FZ.f"
	} else if (*lzwork < wrkopt) {
#line 427 "TG01FZ.f"
	    *info = -26;
#line 428 "TG01FZ.f"
	}
#line 429 "TG01FZ.f"
    }
#line 430 "TG01FZ.f"
    if (*info != 0) {
#line 431 "TG01FZ.f"
	i__1 = -(*info);
#line 431 "TG01FZ.f"
	xerbla_("TG01FZ", &i__1, (ftnlen)6);
#line 432 "TG01FZ.f"
	return 0;
#line 433 "TG01FZ.f"
    } else if (lquery) {
#line 434 "TG01FZ.f"
	zwork[1].r = (doublereal) wrkopt, zwork[1].i = 0.;
#line 435 "TG01FZ.f"
	return 0;
#line 436 "TG01FZ.f"
    }

/*     Initialize Q and Z if necessary. */

#line 440 "TG01FZ.f"
    if (icompq == 3) {
#line 440 "TG01FZ.f"
	zlaset_("Full", l, l, &c_b2, &c_b1, &q[q_offset], ldq, (ftnlen)4);
#line 440 "TG01FZ.f"
    }
#line 442 "TG01FZ.f"
    if (icompz == 3) {
#line 442 "TG01FZ.f"
	zlaset_("Full", n, n, &c_b2, &c_b1, &z__[z_offset], ldz, (ftnlen)4);
#line 442 "TG01FZ.f"
    }

/*     Quick return if possible. */

#line 447 "TG01FZ.f"
    if (*l == 0 || *n == 0) {
#line 448 "TG01FZ.f"
	zwork[1].r = 1., zwork[1].i = 0.;
#line 449 "TG01FZ.f"
	*ranke = 0;
#line 450 "TG01FZ.f"
	if (reda || redtr) {
#line 450 "TG01FZ.f"
	    *rnka22 = 0;
#line 450 "TG01FZ.f"
	}
#line 451 "TG01FZ.f"
	return 0;
#line 452 "TG01FZ.f"
    }

#line 454 "TG01FZ.f"
    toldef = *tol;
#line 455 "TG01FZ.f"
    if (toldef <= 0.) {

/*        Use the default tolerance for rank determination. */

#line 459 "TG01FZ.f"
	toldef = (doublereal) (*l * *n) * dlamch_("EPSILON", (ftnlen)7);
#line 460 "TG01FZ.f"
    }

/*     Set the estimate of maximum singular value of E to */
/*     max(||E||,||A||) to detect negligible A or E matrices. */

/* Computing MAX */
#line 465 "TG01FZ.f"
    d__1 = zlange_("F", l, n, &e[e_offset], lde, &dwork[1], (ftnlen)1), d__2 =
	     zlange_("F", l, n, &a[a_offset], lda, &dwork[1], (ftnlen)1);
#line 465 "TG01FZ.f"
    svlmax = max(d__1,d__2);

/*     Compute the rank-revealing QR decomposition of E, */

/*                        ( E11 E12 ) */
/*           E * P = Qr * (         ) , */
/*                        (  0  E22 ) */

/*     and determine the rank of E using incremental condition */
/*     estimation. */
/*     Complex Workspace: MIN(L,N) + 3*N - 1. */
/*     Real Workspace:    2*N. */

#line 479 "TG01FZ.f"
    lwr = *lzwork - ln;
#line 480 "TG01FZ.f"
    kw = ln + 1;

#line 482 "TG01FZ.f"
    mb3oyz_(l, n, &e[e_offset], lde, &toldef, &svlmax, ranke, sval, &iwork[1],
	     &zwork[1], &dwork[1], &zwork[kw], info);

/*     Apply transformation on the rest of matrices. */

#line 487 "TG01FZ.f"
    if (*ranke > 0) {

/*        A <-- Qr' * A. */
/*        Complex Workspace: need   MIN(L,N) + N; */
/*                           prefer MIN(L,N) + N*NB. */

#line 493 "TG01FZ.f"
	zunmqr_("Left", "ConjTranspose", l, n, ranke, &e[e_offset], lde, &
		zwork[1], &a[a_offset], lda, &zwork[kw], &lwr, info, (ftnlen)
		4, (ftnlen)13);
/* Computing MAX */
#line 495 "TG01FZ.f"
	i__3 = kw;
#line 495 "TG01FZ.f"
	i__1 = wrkopt, i__2 = ln + (integer) zwork[i__3].r;
#line 495 "TG01FZ.f"
	wrkopt = max(i__1,i__2);

/*        B <-- Qr' * B. */
/*        Complex Workspace: need   MIN(L,N) + M; */
/*                           prefer MIN(L,N) + M*NB. */

#line 501 "TG01FZ.f"
	if (withb) {
#line 502 "TG01FZ.f"
	    zunmqr_("Left", "ConjTranspose", l, m, ranke, &e[e_offset], lde, &
		    zwork[1], &b[b_offset], ldb, &zwork[kw], &lwr, info, (
		    ftnlen)4, (ftnlen)13);
/* Computing MAX */
#line 504 "TG01FZ.f"
	    i__3 = kw;
#line 504 "TG01FZ.f"
	    i__1 = wrkopt, i__2 = ln + (integer) zwork[i__3].r;
#line 504 "TG01FZ.f"
	    wrkopt = max(i__1,i__2);
#line 505 "TG01FZ.f"
	}

/*        Q <-- Q * Qr. */
/*        Complex Workspace: need   MIN(L,N) + L; */
/*                           prefer MIN(L,N) + L*NB. */

#line 511 "TG01FZ.f"
	if (ilq) {
#line 512 "TG01FZ.f"
	    zunmqr_("Right", "No Transpose", l, l, ranke, &e[e_offset], lde, &
		    zwork[1], &q[q_offset], ldq, &zwork[kw], &lwr, info, (
		    ftnlen)5, (ftnlen)12);
/* Computing MAX */
#line 514 "TG01FZ.f"
	    i__3 = kw;
#line 514 "TG01FZ.f"
	    i__1 = wrkopt, i__2 = ln + (integer) zwork[i__3].r;
#line 514 "TG01FZ.f"
	    wrkopt = max(i__1,i__2);
#line 515 "TG01FZ.f"
	}

/*        Set lower triangle of E to zero. */

#line 519 "TG01FZ.f"
	if (*l >= 2) {
#line 519 "TG01FZ.f"
	    i__1 = *l - 1;
#line 519 "TG01FZ.f"
	    zlaset_("Lower", &i__1, ranke, &c_b2, &c_b2, &e[e_dim1 + 2], lde, 
		    (ftnlen)5);
#line 519 "TG01FZ.f"
	}

/*        Compute A*P, C*P and Z*P by forward permuting the columns of */
/*        A, C and Z based on information in IWORK. */

#line 525 "TG01FZ.f"
	i__1 = *n;
#line 525 "TG01FZ.f"
	for (j = 1; j <= i__1; ++j) {
#line 526 "TG01FZ.f"
	    iwork[j] = -iwork[j];
#line 527 "TG01FZ.f"
/* L10: */
#line 527 "TG01FZ.f"
	}
#line 528 "TG01FZ.f"
	i__1 = *n;
#line 528 "TG01FZ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 529 "TG01FZ.f"
	    if (iwork[i__] < 0) {
#line 530 "TG01FZ.f"
		j = i__;
#line 531 "TG01FZ.f"
		iwork[j] = -iwork[j];
#line 532 "TG01FZ.f"
L20:
#line 533 "TG01FZ.f"
		k = iwork[j];
#line 534 "TG01FZ.f"
		if (iwork[k] < 0) {
#line 535 "TG01FZ.f"
		    zswap_(l, &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
			    c__1);
#line 536 "TG01FZ.f"
		    if (withc) {
#line 536 "TG01FZ.f"
			zswap_(p, &c__[j * c_dim1 + 1], &c__1, &c__[k * 
				c_dim1 + 1], &c__1);
#line 536 "TG01FZ.f"
		    }
#line 538 "TG01FZ.f"
		    if (ilz) {
#line 538 "TG01FZ.f"
			zswap_(n, &z__[j * z_dim1 + 1], &c__1, &z__[k * 
				z_dim1 + 1], &c__1);
#line 538 "TG01FZ.f"
		    }
#line 540 "TG01FZ.f"
		    iwork[k] = -iwork[k];
#line 541 "TG01FZ.f"
		    j = k;
#line 542 "TG01FZ.f"
		    goto L20;
#line 543 "TG01FZ.f"
		}
#line 544 "TG01FZ.f"
	    }
#line 545 "TG01FZ.f"
/* L30: */
#line 545 "TG01FZ.f"
	}

/*        Determine a unitary matrix Y such that */

/*           ( E11 E12 ) = ( Er  0 ) * Y . */

/*        Compute E <-- E*Y', A <-- A*Y', C <-- C*Y', Z <-- Z*Y'. */

#line 553 "TG01FZ.f"
	if (*ranke < *n) {

/*           Complex Workspace: need   2*N; */
/*                              prefer N + N*NB. */

#line 558 "TG01FZ.f"
	    kw = *ranke + 1;
#line 559 "TG01FZ.f"
	    i__1 = *lzwork - kw + 1;
#line 559 "TG01FZ.f"
	    ztzrzf_(ranke, n, &e[e_offset], lde, &zwork[1], &zwork[kw], &i__1,
		     info);
/* Computing MAX */
#line 561 "TG01FZ.f"
	    i__3 = kw;
#line 561 "TG01FZ.f"
	    i__1 = wrkopt, i__2 = (integer) zwork[i__3].r + kw - 1;
#line 561 "TG01FZ.f"
	    wrkopt = max(i__1,i__2);

/*           Complex Workspace: need   N + MAX(L,P,N); */
/*                              prefer N + MAX(L,P,N)*NB. */

#line 566 "TG01FZ.f"
	    lh = *n - *ranke;
#line 567 "TG01FZ.f"
	    i__1 = *lzwork - kw + 1;
#line 567 "TG01FZ.f"
	    zunmrz_("Right", "Conjugate transpose", l, n, ranke, &lh, &e[
		    e_offset], lde, &zwork[1], &a[a_offset], lda, &zwork[kw], 
		    &i__1, info, (ftnlen)5, (ftnlen)19);
/* Computing MAX */
#line 570 "TG01FZ.f"
	    i__3 = kw;
#line 570 "TG01FZ.f"
	    i__1 = wrkopt, i__2 = (integer) zwork[i__3].r + kw - 1;
#line 570 "TG01FZ.f"
	    wrkopt = max(i__1,i__2);
#line 571 "TG01FZ.f"
	    if (withc) {
#line 572 "TG01FZ.f"
		i__1 = *lzwork - kw + 1;
#line 572 "TG01FZ.f"
		zunmrz_("Right", "Conjugate transpose", p, n, ranke, &lh, &e[
			e_offset], lde, &zwork[1], &c__[c_offset], ldc, &
			zwork[kw], &i__1, info, (ftnlen)5, (ftnlen)19);
/* Computing MAX */
#line 575 "TG01FZ.f"
		i__3 = kw;
#line 575 "TG01FZ.f"
		i__1 = wrkopt, i__2 = (integer) zwork[i__3].r + kw - 1;
#line 575 "TG01FZ.f"
		wrkopt = max(i__1,i__2);
#line 576 "TG01FZ.f"
	    }
#line 577 "TG01FZ.f"
	    if (ilz) {
#line 578 "TG01FZ.f"
		i__1 = *lzwork - kw + 1;
#line 578 "TG01FZ.f"
		zunmrz_("Right", "Conjugate transpose", n, n, ranke, &lh, &e[
			e_offset], lde, &zwork[1], &z__[z_offset], ldz, &
			zwork[kw], &i__1, info, (ftnlen)5, (ftnlen)19);
/* Computing MAX */
#line 581 "TG01FZ.f"
		i__3 = kw;
#line 581 "TG01FZ.f"
		i__1 = wrkopt, i__2 = (integer) zwork[i__3].r + kw - 1;
#line 581 "TG01FZ.f"
		wrkopt = max(i__1,i__2);
#line 582 "TG01FZ.f"
	    }

/*           Set E12 and E22 to zero. */

#line 586 "TG01FZ.f"
	    zlaset_("Full", l, &lh, &c_b2, &c_b2, &e[kw * e_dim1 + 1], lde, (
		    ftnlen)4);
#line 587 "TG01FZ.f"
	}
#line 588 "TG01FZ.f"
    } else {
#line 589 "TG01FZ.f"
	zlaset_("Full", l, n, &c_b2, &c_b2, &e[e_offset], lde, (ftnlen)4);
#line 590 "TG01FZ.f"
    }

/*     Reduce A22 if necessary. */

#line 594 "TG01FZ.f"
    if (reda || redtr) {
#line 595 "TG01FZ.f"
	la22 = *l - *ranke;
#line 596 "TG01FZ.f"
	na22 = *n - *ranke;
#line 597 "TG01FZ.f"
	if (min(la22,na22) == 0) {
#line 598 "TG01FZ.f"
	    *rnka22 = 0;
#line 599 "TG01FZ.f"
	} else {

/*           Compute the rank-revealing QR decomposition of A22, */

/*                              ( R11 R12 ) */
/*              A22 * P2 = Q2 * (         ) , */
/*                              (  0  R22 ) */

/*           and determine the rank of A22 using incremental */
/*           condition estimation. */
/*           Complex Workspace: MIN(L,N) + 3*N - 1. */
/*           Real Workspace:    2*N. */

#line 612 "TG01FZ.f"
	    ir1 = *ranke + 1;
#line 613 "TG01FZ.f"
	    mb3oyz_(&la22, &na22, &a[ir1 + ir1 * a_dim1], lda, &toldef, &
		    svlmax, rnka22, sval, &iwork[1], &zwork[1], &dwork[1], &
		    zwork[kw], info);

/*           Apply transformation on the rest of matrices. */

#line 619 "TG01FZ.f"
	    if (*rnka22 > 0) {

/*              A <-- diag(I, Q2') * A */
/*              Complex Workspace: need   MIN(L,N) + N; */
/*                                 prefer MIN(L,N) + N*NB. */

#line 625 "TG01FZ.f"
		zunmqr_("Left", "ConjTranspose", &la22, ranke, rnka22, &a[ir1 
			+ ir1 * a_dim1], lda, &zwork[1], &a[ir1 + a_dim1], 
			lda, &zwork[kw], &lwr, info, (ftnlen)4, (ftnlen)13);

/*              B <-- diag(I, Q2') * B */
/*              Complex Workspace: need   MIN(L,N) + M; */
/*                                 prefer MIN(L,N) + M*NB. */

#line 633 "TG01FZ.f"
		if (withb) {
#line 633 "TG01FZ.f"
		    zunmqr_("Left", "ConjTranspose", &la22, m, rnka22, &a[ir1 
			    + ir1 * a_dim1], lda, &zwork[1], &b[ir1 + b_dim1],
			     ldb, &zwork[kw], &lwr, info, (ftnlen)4, (ftnlen)
			    13);
#line 633 "TG01FZ.f"
		}

/*              Q <-- Q * diag(I, Q2) */
/*              Complex Workspace: need   MIN(L,N) + L; */
/*                                 prefer MIN(L,N) + L*NB. */

#line 642 "TG01FZ.f"
		if (ilq) {
#line 642 "TG01FZ.f"
		    zunmqr_("Right", "No transpose", l, &la22, rnka22, &a[ir1 
			    + ir1 * a_dim1], lda, &zwork[1], &q[ir1 * q_dim1 
			    + 1], ldq, &zwork[kw], &lwr, info, (ftnlen)5, (
			    ftnlen)12);
#line 642 "TG01FZ.f"
		}

/*              Set lower triangle of A22 to zero. */

#line 649 "TG01FZ.f"
		if (la22 >= 2) {
#line 649 "TG01FZ.f"
		    i__1 = la22 - 1;
#line 649 "TG01FZ.f"
		    zlaset_("Lower", &i__1, rnka22, &c_b2, &c_b2, &a[ir1 + 1 
			    + ir1 * a_dim1], lda, (ftnlen)5);
#line 649 "TG01FZ.f"
		}

/*              Compute A*diag(I,P2), C*diag(I,P2) and Z*diag(I,P2) */
/*              by forward permuting the columns of A, C and Z based */
/*              on information in IWORK. */

#line 657 "TG01FZ.f"
		i__1 = na22;
#line 657 "TG01FZ.f"
		for (j = 1; j <= i__1; ++j) {
#line 658 "TG01FZ.f"
		    iwork[j] = -iwork[j];
#line 659 "TG01FZ.f"
/* L40: */
#line 659 "TG01FZ.f"
		}
#line 660 "TG01FZ.f"
		i__1 = na22;
#line 660 "TG01FZ.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 661 "TG01FZ.f"
		    if (iwork[i__] < 0) {
#line 662 "TG01FZ.f"
			j = i__;
#line 663 "TG01FZ.f"
			iwork[j] = -iwork[j];
#line 664 "TG01FZ.f"
L50:
#line 665 "TG01FZ.f"
			k = iwork[j];
#line 666 "TG01FZ.f"
			if (iwork[k] < 0) {
#line 667 "TG01FZ.f"
			    zswap_(ranke, &a[(*ranke + j) * a_dim1 + 1], &
				    c__1, &a[(*ranke + k) * a_dim1 + 1], &
				    c__1);
#line 669 "TG01FZ.f"
			    if (withc) {
#line 669 "TG01FZ.f"
				zswap_(p, &c__[(*ranke + j) * c_dim1 + 1], &
					c__1, &c__[(*ranke + k) * c_dim1 + 1],
					 &c__1);
#line 669 "TG01FZ.f"
			    }
#line 672 "TG01FZ.f"
			    if (ilz) {
#line 672 "TG01FZ.f"
				zswap_(n, &z__[(*ranke + j) * z_dim1 + 1], &
					c__1, &z__[(*ranke + k) * z_dim1 + 1],
					 &c__1);
#line 672 "TG01FZ.f"
			    }
#line 675 "TG01FZ.f"
			    iwork[k] = -iwork[k];
#line 676 "TG01FZ.f"
			    j = k;
#line 677 "TG01FZ.f"
			    goto L50;
#line 678 "TG01FZ.f"
			}
#line 679 "TG01FZ.f"
		    }
#line 680 "TG01FZ.f"
/* L60: */
#line 680 "TG01FZ.f"
		}

#line 682 "TG01FZ.f"
		if (reda && *rnka22 < na22) {

/*                 Determine a unitary matrix Y2 such that */

/*                 ( R11 R12 ) = ( Ar  0 ) * Y2 . */

/*                 Compute A <-- A*diag(I, Y2'), C <-- C*diag(I, Y2'), */
/*                         Z <-- Z*diag(I, Y2'). */

/*                 Complex Workspace: need   2*N; */
/*                                    prefer N + N*NB. */

#line 694 "TG01FZ.f"
		    kw = *ranke + 1;
#line 695 "TG01FZ.f"
		    i__1 = *lzwork - kw + 1;
#line 695 "TG01FZ.f"
		    ztzrzf_(rnka22, &na22, &a[ir1 + ir1 * a_dim1], lda, &
			    zwork[1], &zwork[kw], &i__1, info);
/* Computing MAX */
#line 697 "TG01FZ.f"
		    i__3 = kw;
#line 697 "TG01FZ.f"
		    i__1 = wrkopt, i__2 = (integer) zwork[i__3].r + kw - 1;
#line 697 "TG01FZ.f"
		    wrkopt = max(i__1,i__2);

/*                 Complex Workspace: need   N + MAX(P,N); */
/*                                    prefer N + MAX(P,N)*NB. */

#line 702 "TG01FZ.f"
		    lh = na22 - *rnka22;
#line 703 "TG01FZ.f"
		    if (withc) {
#line 704 "TG01FZ.f"
			i__1 = *lzwork - kw + 1;
#line 704 "TG01FZ.f"
			zunmrz_("Right", "Conjugate transpose", p, n, rnka22, 
				&lh, &a[ir1 + ir1 * a_dim1], lda, &zwork[1], &
				c__[c_offset], ldc, &zwork[kw], &i__1, info, (
				ftnlen)5, (ftnlen)19);
/* Computing MAX */
#line 707 "TG01FZ.f"
			i__3 = kw;
#line 707 "TG01FZ.f"
			i__1 = wrkopt, i__2 = (integer) zwork[i__3].r + kw - 
				1;
#line 707 "TG01FZ.f"
			wrkopt = max(i__1,i__2);
#line 708 "TG01FZ.f"
		    }
#line 709 "TG01FZ.f"
		    if (ilz) {
#line 710 "TG01FZ.f"
			i__1 = *lzwork - kw + 1;
#line 710 "TG01FZ.f"
			zunmrz_("Right", "Conjugate transpose", n, n, rnka22, 
				&lh, &a[ir1 + ir1 * a_dim1], lda, &zwork[1], &
				z__[z_offset], ldz, &zwork[kw], &i__1, info, (
				ftnlen)5, (ftnlen)19);
/* Computing MAX */
#line 713 "TG01FZ.f"
			i__3 = kw;
#line 713 "TG01FZ.f"
			i__1 = wrkopt, i__2 = (integer) zwork[i__3].r + kw - 
				1;
#line 713 "TG01FZ.f"
			wrkopt = max(i__1,i__2);
#line 714 "TG01FZ.f"
		    }
#line 715 "TG01FZ.f"
		    ire1 = *ranke + *rnka22 + 1;

/*                 Set R12 and R22 to zero. */

#line 719 "TG01FZ.f"
		    zlaset_("Full", &la22, &lh, &c_b2, &c_b2, &a[ir1 + ire1 * 
			    a_dim1], lda, (ftnlen)4);
#line 721 "TG01FZ.f"
		}
#line 722 "TG01FZ.f"
	    } else {
#line 723 "TG01FZ.f"
		zlaset_("Full", &la22, &na22, &c_b2, &c_b2, &a[ir1 + ir1 * 
			a_dim1], lda, (ftnlen)4);
#line 725 "TG01FZ.f"
	    }
#line 726 "TG01FZ.f"
	}
#line 727 "TG01FZ.f"
    }

#line 729 "TG01FZ.f"
    zwork[1].r = (doublereal) wrkopt, zwork[1].i = 0.;

#line 731 "TG01FZ.f"
    return 0;
/* *** Last line of TG01FZ *** */
} /* tg01fz_ */

