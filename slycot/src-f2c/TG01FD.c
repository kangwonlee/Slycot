#line 1 "TG01FD.f"
/* TG01FD.f -- translated by f2c (version 20100827).
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

#line 1 "TG01FD.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b42 = 0.;
static doublereal c_b43 = 1.;

/* Subroutine */ int tg01fd_(char *compq, char *compz, char *joba, integer *l,
	 integer *n, integer *m, integer *p, doublereal *a, integer *lda, 
	doublereal *e, integer *lde, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *q, integer *ldq, doublereal *z__, 
	integer *ldz, integer *ranke, integer *rnka22, doublereal *tol, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen compq_len, ftnlen compz_len, ftnlen joba_len)
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
    extern /* Subroutine */ int mb03oy_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *);
    static logical withb, withc, redtr;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static doublereal toldef;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer icompq, icompz;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static doublereal svlmax;
    extern /* Subroutine */ int dormrz_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static logical lquery;
    extern /* Subroutine */ int dtzrzf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
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

/*     To compute for the descriptor system (A-lambda E,B,C) */
/*     the orthogonal transformation matrices Q and Z such that the */
/*     transformed system (Q'*A*Z-lambda Q'*E*Z, Q'*B, C*Z) is */
/*     in a SVD-like coordinate form with */

/*                  ( A11  A12 )             ( Er  0 ) */
/*         Q'*A*Z = (          ) ,  Q'*E*Z = (       ) , */
/*                  ( A21  A22 )             (  0  0 ) */

/*     where Er is an upper triangular invertible matrix. */
/*     Optionally, the A22 matrix can be further reduced to the form */

/*                  ( Ar  X ) */
/*            A22 = (       ) , */
/*                  (  0  0 ) */

/*     with Ar an upper triangular invertible matrix, and X either a full */
/*     or a zero matrix. */
/*     The left and/or right orthogonal transformations performed */
/*     to reduce E and A22 can be optionally accumulated. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     COMPQ   CHARACTER*1 */
/*             = 'N':  do not compute Q; */
/*             = 'I':  Q is initialized to the unit matrix, and the */
/*                     orthogonal matrix Q is returned; */
/*             = 'U':  Q must contain an orthogonal matrix Q1 on entry, */
/*                     and the product Q1*Q is returned. */

/*     COMPZ   CHARACTER*1 */
/*             = 'N':  do not compute Z; */
/*             = 'I':  Z is initialized to the unit matrix, and the */
/*                     orthogonal matrix Z is returned; */
/*             = 'U':  Z must contain an orthogonal matrix Z1 on entry, */
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

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
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

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
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

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading L-by-M part of this array must */
/*             contain the input/state matrix B. */
/*             On exit, the leading L-by-M part of this array contains */
/*             the transformed matrix Q'*B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B. */
/*             LDB >= MAX(1,L) if M > 0 or LDB >= 1 if M = 0. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the transformed matrix C*Z. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,L) */
/*             If COMPQ = 'N':  Q is not referenced. */
/*             If COMPQ = 'I':  on entry, Q need not be set; */
/*                              on exit, the leading L-by-L part of this */
/*                              array contains the orthogonal matrix Q, */
/*                              where Q' is the product of Householder */
/*                              transformations which are applied to A, */
/*                              E, and B on the left. */
/*             If COMPQ = 'U':  on entry, the leading L-by-L part of this */
/*                              array must contain an orthogonal matrix */
/*                              Q1; */
/*                              on exit, the leading L-by-L part of this */
/*                              array contains the orthogonal matrix */
/*                              Q1*Q. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. */
/*             LDQ >= 1,        if COMPQ = 'N'; */
/*             LDQ >= MAX(1,L), if COMPQ = 'U' or 'I'. */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             If COMPZ = 'N':  Z is not referenced. */
/*             If COMPZ = 'I':  on entry, Z need not be set; */
/*                              on exit, the leading N-by-N part of this */
/*                              array contains the orthogonal matrix Z, */
/*                              which is the product of Householder */
/*                              transformations applied to A, E, and C */
/*                              on the right. */
/*             If COMPZ = 'U':  on entry, the leading N-by-N part of this */
/*                              array must contain an orthogonal matrix */
/*                              Z1; */
/*                              on exit, the leading N-by-N part of this */
/*                              array contains the orthogonal matrix */
/*                              Z1*Z. */

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

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 1, N+P, MIN(L,N)+MAX(3*N-1,M,L) ). */
/*             For optimal performance, LDWORK should be larger. */

/*             If LDWORK = -1, then a workspace query is assumed; */
/*             the routine only calculates the optimal size of the */
/*             DWORK array, returns this value as the first entry of */
/*             the DWORK array, and no error message related to LDWORK */
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
/*     zero, and an orthogonal matrix Y is determined such that */

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
/*     March 1999. Based on the RASP routine RPDSSV. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, July 1999, */
/*     May 2003, Jan. 2009. */

/*     KEYWORDS */

/*     Descriptor system, matrix algebra, matrix operations, */
/*     orthogonal transformation. */

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

#line 331 "TG01FD.f"
    /* Parameter adjustments */
#line 331 "TG01FD.f"
    a_dim1 = *lda;
#line 331 "TG01FD.f"
    a_offset = 1 + a_dim1;
#line 331 "TG01FD.f"
    a -= a_offset;
#line 331 "TG01FD.f"
    e_dim1 = *lde;
#line 331 "TG01FD.f"
    e_offset = 1 + e_dim1;
#line 331 "TG01FD.f"
    e -= e_offset;
#line 331 "TG01FD.f"
    b_dim1 = *ldb;
#line 331 "TG01FD.f"
    b_offset = 1 + b_dim1;
#line 331 "TG01FD.f"
    b -= b_offset;
#line 331 "TG01FD.f"
    c_dim1 = *ldc;
#line 331 "TG01FD.f"
    c_offset = 1 + c_dim1;
#line 331 "TG01FD.f"
    c__ -= c_offset;
#line 331 "TG01FD.f"
    q_dim1 = *ldq;
#line 331 "TG01FD.f"
    q_offset = 1 + q_dim1;
#line 331 "TG01FD.f"
    q -= q_offset;
#line 331 "TG01FD.f"
    z_dim1 = *ldz;
#line 331 "TG01FD.f"
    z_offset = 1 + z_dim1;
#line 331 "TG01FD.f"
    z__ -= z_offset;
#line 331 "TG01FD.f"
    --iwork;
#line 331 "TG01FD.f"
    --dwork;
#line 331 "TG01FD.f"

#line 331 "TG01FD.f"
    /* Function Body */
#line 331 "TG01FD.f"
    if (lsame_(compq, "N", (ftnlen)1, (ftnlen)1)) {
#line 332 "TG01FD.f"
	ilq = FALSE_;
#line 333 "TG01FD.f"
	icompq = 1;
#line 334 "TG01FD.f"
    } else if (lsame_(compq, "U", (ftnlen)1, (ftnlen)1)) {
#line 335 "TG01FD.f"
	ilq = TRUE_;
#line 336 "TG01FD.f"
	icompq = 2;
#line 337 "TG01FD.f"
    } else if (lsame_(compq, "I", (ftnlen)1, (ftnlen)1)) {
#line 338 "TG01FD.f"
	ilq = TRUE_;
#line 339 "TG01FD.f"
	icompq = 3;
#line 340 "TG01FD.f"
    } else {
#line 341 "TG01FD.f"
	icompq = 0;
#line 342 "TG01FD.f"
    }

/*     Decode COMPZ. */

#line 346 "TG01FD.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 347 "TG01FD.f"
	ilz = FALSE_;
#line 348 "TG01FD.f"
	icompz = 1;
#line 349 "TG01FD.f"
    } else if (lsame_(compz, "U", (ftnlen)1, (ftnlen)1)) {
#line 350 "TG01FD.f"
	ilz = TRUE_;
#line 351 "TG01FD.f"
	icompz = 2;
#line 352 "TG01FD.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 353 "TG01FD.f"
	ilz = TRUE_;
#line 354 "TG01FD.f"
	icompz = 3;
#line 355 "TG01FD.f"
    } else {
#line 356 "TG01FD.f"
	icompz = 0;
#line 357 "TG01FD.f"
    }
#line 358 "TG01FD.f"
    reda = lsame_(joba, "R", (ftnlen)1, (ftnlen)1);
#line 359 "TG01FD.f"
    redtr = lsame_(joba, "T", (ftnlen)1, (ftnlen)1);
#line 360 "TG01FD.f"
    withb = *m > 0;
#line 361 "TG01FD.f"
    withc = *p > 0;
#line 362 "TG01FD.f"
    lquery = *ldwork == -1;

/*     Test the input parameters. */

#line 366 "TG01FD.f"
    ln = min(*l,*n);
#line 367 "TG01FD.f"
    *info = 0;
/* Computing MAX */
/* Computing MAX */
#line 368 "TG01FD.f"
    i__3 = *n * 3 - 1, i__3 = max(i__3,*m);
#line 368 "TG01FD.f"
    i__1 = 1, i__2 = *n + *p, i__1 = max(i__1,i__2), i__2 = ln + max(i__3,*l);
#line 368 "TG01FD.f"
    wrkopt = max(i__1,i__2);
#line 369 "TG01FD.f"
    if (icompq <= 0) {
#line 370 "TG01FD.f"
	*info = -1;
#line 371 "TG01FD.f"
    } else if (icompz <= 0) {
#line 372 "TG01FD.f"
	*info = -2;
#line 373 "TG01FD.f"
    } else if (! lsame_(joba, "N", (ftnlen)1, (ftnlen)1) && ! reda && ! redtr)
	     {
#line 375 "TG01FD.f"
	*info = -3;
#line 376 "TG01FD.f"
    } else if (*l < 0) {
#line 377 "TG01FD.f"
	*info = -4;
#line 378 "TG01FD.f"
    } else if (*n < 0) {
#line 379 "TG01FD.f"
	*info = -5;
#line 380 "TG01FD.f"
    } else if (*m < 0) {
#line 381 "TG01FD.f"
	*info = -6;
#line 382 "TG01FD.f"
    } else if (*p < 0) {
#line 383 "TG01FD.f"
	*info = -7;
#line 384 "TG01FD.f"
    } else if (*lda < max(1,*l)) {
#line 385 "TG01FD.f"
	*info = -9;
#line 386 "TG01FD.f"
    } else if (*lde < max(1,*l)) {
#line 387 "TG01FD.f"
	*info = -11;
#line 388 "TG01FD.f"
    } else if (*ldb < 1 || withb && *ldb < *l) {
#line 389 "TG01FD.f"
	*info = -13;
#line 390 "TG01FD.f"
    } else if (*ldc < max(1,*p)) {
#line 391 "TG01FD.f"
	*info = -15;
#line 392 "TG01FD.f"
    } else if (ilq && *ldq < *l || *ldq < 1) {
#line 393 "TG01FD.f"
	*info = -17;
#line 394 "TG01FD.f"
    } else if (ilz && *ldz < *n || *ldz < 1) {
#line 395 "TG01FD.f"
	*info = -19;
#line 396 "TG01FD.f"
    } else if (*tol >= 1.) {
#line 397 "TG01FD.f"
	*info = -22;
#line 398 "TG01FD.f"
    } else {
#line 399 "TG01FD.f"
	if (lquery) {
/* Computing MIN */
#line 400 "TG01FD.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "DORMQR", "LC", l, n, &ln, &c_n1,
		     (ftnlen)6, (ftnlen)2);
#line 400 "TG01FD.f"
	    nb = min(i__1,i__2);
/* Computing MAX */
#line 401 "TG01FD.f"
	    i__1 = wrkopt, i__2 = ln + *n * nb;
#line 401 "TG01FD.f"
	    wrkopt = max(i__1,i__2);
#line 402 "TG01FD.f"
	    if (withb) {
/* Computing MIN */
#line 403 "TG01FD.f"
		i__1 = 64, i__2 = ilaenv_(&c__1, "DORMQR", "LC", l, m, &ln, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 403 "TG01FD.f"
		nb = min(i__1,i__2);
/* Computing MAX */
#line 404 "TG01FD.f"
		i__1 = wrkopt, i__2 = ln + *m * nb;
#line 404 "TG01FD.f"
		wrkopt = max(i__1,i__2);
#line 405 "TG01FD.f"
	    }
#line 406 "TG01FD.f"
	    if (ilq) {
/* Computing MIN */
#line 407 "TG01FD.f"
		i__1 = 64, i__2 = ilaenv_(&c__1, "DORMQR", "RN", l, l, &ln, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 407 "TG01FD.f"
		nb = min(i__1,i__2);
/* Computing MAX */
#line 408 "TG01FD.f"
		i__1 = wrkopt, i__2 = ln + *l * nb;
#line 408 "TG01FD.f"
		wrkopt = max(i__1,i__2);
#line 409 "TG01FD.f"
	    }
#line 410 "TG01FD.f"
	    nb = ilaenv_(&c__1, "DGERQF", " ", l, n, &c_n1, &c_n1, (ftnlen)6, 
		    (ftnlen)1);
/* Computing MAX */
#line 411 "TG01FD.f"
	    i__1 = wrkopt, i__2 = ln + *n * nb;
#line 411 "TG01FD.f"
	    wrkopt = max(i__1,i__2);
/* Computing MIN */
#line 412 "TG01FD.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "DORMRQ", "RC", l, n, n, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 412 "TG01FD.f"
	    nb = min(i__1,i__2);
/* Computing MAX */
#line 413 "TG01FD.f"
	    i__1 = wrkopt, i__2 = *n + max(1,*l) * nb;
#line 413 "TG01FD.f"
	    wrkopt = max(i__1,i__2);
#line 414 "TG01FD.f"
	    if (withc) {
/* Computing MIN */
#line 415 "TG01FD.f"
		i__1 = 64, i__2 = ilaenv_(&c__1, "DORMRQ", "RC", p, n, n, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 415 "TG01FD.f"
		nb = min(i__1,i__2);
/* Computing MAX */
#line 416 "TG01FD.f"
		i__1 = wrkopt, i__2 = *n + max(1,*p) * nb;
#line 416 "TG01FD.f"
		wrkopt = max(i__1,i__2);
#line 417 "TG01FD.f"
	    }
#line 418 "TG01FD.f"
	    if (ilz) {
/* Computing MIN */
#line 419 "TG01FD.f"
		i__1 = 64, i__2 = ilaenv_(&c__1, "DORMRQ", "RC", n, n, n, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 419 "TG01FD.f"
		nb = min(i__1,i__2);
/* Computing MAX */
#line 420 "TG01FD.f"
		i__1 = wrkopt, i__2 = *n + max(1,*n) * nb;
#line 420 "TG01FD.f"
		wrkopt = max(i__1,i__2);
#line 421 "TG01FD.f"
	    }
#line 422 "TG01FD.f"
	} else if (*ldwork < wrkopt) {
#line 423 "TG01FD.f"
	    *info = -25;
#line 424 "TG01FD.f"
	}
#line 425 "TG01FD.f"
    }
#line 426 "TG01FD.f"
    if (*info != 0) {
#line 427 "TG01FD.f"
	i__1 = -(*info);
#line 427 "TG01FD.f"
	xerbla_("TG01FD", &i__1, (ftnlen)6);
#line 428 "TG01FD.f"
	return 0;
#line 429 "TG01FD.f"
    } else if (lquery) {
#line 430 "TG01FD.f"
	dwork[1] = (doublereal) wrkopt;
#line 431 "TG01FD.f"
	return 0;
#line 432 "TG01FD.f"
    }

/*     Initialize Q and Z if necessary. */

#line 436 "TG01FD.f"
    if (icompq == 3) {
#line 436 "TG01FD.f"
	dlaset_("Full", l, l, &c_b42, &c_b43, &q[q_offset], ldq, (ftnlen)4);
#line 436 "TG01FD.f"
    }
#line 438 "TG01FD.f"
    if (icompz == 3) {
#line 438 "TG01FD.f"
	dlaset_("Full", n, n, &c_b42, &c_b43, &z__[z_offset], ldz, (ftnlen)4);
#line 438 "TG01FD.f"
    }

/*     Quick return if possible. */

#line 443 "TG01FD.f"
    if (*l == 0 || *n == 0) {
#line 444 "TG01FD.f"
	dwork[1] = 1.;
#line 445 "TG01FD.f"
	*ranke = 0;
#line 446 "TG01FD.f"
	if (reda || redtr) {
#line 446 "TG01FD.f"
	    *rnka22 = 0;
#line 446 "TG01FD.f"
	}
#line 447 "TG01FD.f"
	return 0;
#line 448 "TG01FD.f"
    }

#line 450 "TG01FD.f"
    toldef = *tol;
#line 451 "TG01FD.f"
    if (toldef <= 0.) {

/*        Use the default tolerance for rank determination. */

#line 455 "TG01FD.f"
	toldef = (doublereal) (*l * *n) * dlamch_("EPSILON", (ftnlen)7);
#line 456 "TG01FD.f"
    }

/*     Set the estimate of maximum singular value of E to */
/*     max(||E||,||A||) to detect negligible A or E matrices. */

/* Computing MAX */
#line 461 "TG01FD.f"
    d__1 = dlange_("F", l, n, &e[e_offset], lde, &dwork[1], (ftnlen)1), d__2 =
	     dlange_("F", l, n, &a[a_offset], lda, &dwork[1], (ftnlen)1);
#line 461 "TG01FD.f"
    svlmax = max(d__1,d__2);

/*     Compute the rank-revealing QR decomposition of E, */

/*                        ( E11 E12 ) */
/*           E * P = Qr * (         ) , */
/*                        (  0  E22 ) */

/*     and determine the rank of E using incremental condition */
/*     estimation. */
/*     Workspace: MIN(L,N) + 3*N - 1. */

#line 474 "TG01FD.f"
    lwr = *ldwork - ln;
#line 475 "TG01FD.f"
    kw = ln + 1;

#line 477 "TG01FD.f"
    mb03oy_(l, n, &e[e_offset], lde, &toldef, &svlmax, ranke, sval, &iwork[1],
	     &dwork[1], &dwork[kw], info);

/*     Apply transformation on the rest of matrices. */

#line 482 "TG01FD.f"
    if (*ranke > 0) {

/*        A <-- Qr' * A. */
/*        Workspace: need   MIN(L,N) + N; */
/*                   prefer MIN(L,N) + N*NB. */

#line 488 "TG01FD.f"
	dormqr_("Left", "Transpose", l, n, ranke, &e[e_offset], lde, &dwork[1]
		, &a[a_offset], lda, &dwork[kw], &lwr, info, (ftnlen)4, (
		ftnlen)9);
/* Computing MAX */
#line 490 "TG01FD.f"
	i__1 = wrkopt, i__2 = ln + (integer) dwork[kw];
#line 490 "TG01FD.f"
	wrkopt = max(i__1,i__2);

/*        B <-- Qr' * B. */
/*        Workspace: need   MIN(L,N) + M; */
/*                   prefer MIN(L,N) + M*NB. */

#line 496 "TG01FD.f"
	if (withb) {
#line 497 "TG01FD.f"
	    dormqr_("Left", "Transpose", l, m, ranke, &e[e_offset], lde, &
		    dwork[1], &b[b_offset], ldb, &dwork[kw], &lwr, info, (
		    ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 499 "TG01FD.f"
	    i__1 = wrkopt, i__2 = ln + (integer) dwork[kw];
#line 499 "TG01FD.f"
	    wrkopt = max(i__1,i__2);
#line 500 "TG01FD.f"
	}

/*        Q <-- Q * Qr. */
/*        Workspace: need   MIN(L,N) + L; */
/*                   prefer MIN(L,N) + L*NB. */

#line 506 "TG01FD.f"
	if (ilq) {
#line 507 "TG01FD.f"
	    dormqr_("Right", "No Transpose", l, l, ranke, &e[e_offset], lde, &
		    dwork[1], &q[q_offset], ldq, &dwork[kw], &lwr, info, (
		    ftnlen)5, (ftnlen)12);
/* Computing MAX */
#line 509 "TG01FD.f"
	    i__1 = wrkopt, i__2 = ln + (integer) dwork[kw];
#line 509 "TG01FD.f"
	    wrkopt = max(i__1,i__2);
#line 510 "TG01FD.f"
	}

/*        Set lower triangle of E to zero. */

#line 514 "TG01FD.f"
	if (*l >= 2) {
#line 514 "TG01FD.f"
	    i__1 = *l - 1;
#line 514 "TG01FD.f"
	    dlaset_("Lower", &i__1, ranke, &c_b42, &c_b42, &e[e_dim1 + 2], 
		    lde, (ftnlen)5);
#line 514 "TG01FD.f"
	}

/*        Compute A*P, C*P and Z*P by forward permuting the columns of */
/*        A, C and Z based on information in IWORK. */

#line 520 "TG01FD.f"
	i__1 = *n;
#line 520 "TG01FD.f"
	for (j = 1; j <= i__1; ++j) {
#line 521 "TG01FD.f"
	    iwork[j] = -iwork[j];
#line 522 "TG01FD.f"
/* L10: */
#line 522 "TG01FD.f"
	}
#line 523 "TG01FD.f"
	i__1 = *n;
#line 523 "TG01FD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 524 "TG01FD.f"
	    if (iwork[i__] < 0) {
#line 525 "TG01FD.f"
		j = i__;
#line 526 "TG01FD.f"
		iwork[j] = -iwork[j];
#line 527 "TG01FD.f"
L20:
#line 528 "TG01FD.f"
		k = iwork[j];
#line 529 "TG01FD.f"
		if (iwork[k] < 0) {
#line 530 "TG01FD.f"
		    dswap_(l, &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
			    c__1);
#line 531 "TG01FD.f"
		    if (withc) {
#line 531 "TG01FD.f"
			dswap_(p, &c__[j * c_dim1 + 1], &c__1, &c__[k * 
				c_dim1 + 1], &c__1);
#line 531 "TG01FD.f"
		    }
#line 533 "TG01FD.f"
		    if (ilz) {
#line 533 "TG01FD.f"
			dswap_(n, &z__[j * z_dim1 + 1], &c__1, &z__[k * 
				z_dim1 + 1], &c__1);
#line 533 "TG01FD.f"
		    }
#line 535 "TG01FD.f"
		    iwork[k] = -iwork[k];
#line 536 "TG01FD.f"
		    j = k;
#line 537 "TG01FD.f"
		    goto L20;
#line 538 "TG01FD.f"
		}
#line 539 "TG01FD.f"
	    }
#line 540 "TG01FD.f"
/* L30: */
#line 540 "TG01FD.f"
	}

/*        Determine an orthogonal matrix Y such that */

/*           ( E11 E12 ) = ( Er  0 ) * Y . */

/*        Compute E <-- E*Y', A <-- A*Y', C <-- C*Y', Z <-- Z*Y'. */

#line 548 "TG01FD.f"
	if (*ranke < *n) {

/*           Workspace: need   2*N; */
/*                      prefer N + N*NB. */

#line 553 "TG01FD.f"
	    kw = *ranke + 1;
#line 554 "TG01FD.f"
	    i__1 = *ldwork - kw + 1;
#line 554 "TG01FD.f"
	    dtzrzf_(ranke, n, &e[e_offset], lde, &dwork[1], &dwork[kw], &i__1,
		     info);
/* Computing MAX */
#line 556 "TG01FD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 556 "TG01FD.f"
	    wrkopt = max(i__1,i__2);

/*           Workspace: need   N + MAX(L,P,N); */
/*                      prefer N + MAX(L,P,N)*NB. */

#line 561 "TG01FD.f"
	    lh = *n - *ranke;
#line 562 "TG01FD.f"
	    i__1 = *ldwork - kw + 1;
#line 562 "TG01FD.f"
	    dormrz_("Right", "Transpose", l, n, ranke, &lh, &e[e_offset], lde,
		     &dwork[1], &a[a_offset], lda, &dwork[kw], &i__1, info, (
		    ftnlen)5, (ftnlen)9);
/* Computing MAX */
#line 564 "TG01FD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 564 "TG01FD.f"
	    wrkopt = max(i__1,i__2);
#line 565 "TG01FD.f"
	    if (withc) {
#line 566 "TG01FD.f"
		i__1 = *ldwork - kw + 1;
#line 566 "TG01FD.f"
		dormrz_("Right", "Transpose", p, n, ranke, &lh, &e[e_offset], 
			lde, &dwork[1], &c__[c_offset], ldc, &dwork[kw], &
			i__1, info, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
#line 569 "TG01FD.f"
		i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 569 "TG01FD.f"
		wrkopt = max(i__1,i__2);
#line 570 "TG01FD.f"
	    }
#line 571 "TG01FD.f"
	    if (ilz) {
#line 572 "TG01FD.f"
		i__1 = *ldwork - kw + 1;
#line 572 "TG01FD.f"
		dormrz_("Right", "Transpose", n, n, ranke, &lh, &e[e_offset], 
			lde, &dwork[1], &z__[z_offset], ldz, &dwork[kw], &
			i__1, info, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
#line 575 "TG01FD.f"
		i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 575 "TG01FD.f"
		wrkopt = max(i__1,i__2);
#line 576 "TG01FD.f"
	    }

/*           Set E12 and E22 to zero. */

#line 580 "TG01FD.f"
	    dlaset_("Full", l, &lh, &c_b42, &c_b42, &e[kw * e_dim1 + 1], lde, 
		    (ftnlen)4);
#line 581 "TG01FD.f"
	}
#line 582 "TG01FD.f"
    } else {
#line 583 "TG01FD.f"
	dlaset_("Full", l, n, &c_b42, &c_b42, &e[e_offset], lde, (ftnlen)4);
#line 584 "TG01FD.f"
    }

/*     Reduce A22 if necessary. */

#line 588 "TG01FD.f"
    if (reda || redtr) {
#line 589 "TG01FD.f"
	la22 = *l - *ranke;
#line 590 "TG01FD.f"
	na22 = *n - *ranke;
#line 591 "TG01FD.f"
	if (min(la22,na22) == 0) {
#line 592 "TG01FD.f"
	    *rnka22 = 0;
#line 593 "TG01FD.f"
	} else {

/*           Compute the rank-revealing QR decomposition of A22, */

/*                              ( R11 R12 ) */
/*              A22 * P2 = Q2 * (         ) , */
/*                              (  0  R22 ) */

/*           and determine the rank of A22 using incremental */
/*           condition estimation. */
/*           Workspace: MIN(L,N) + 3*N - 1. */

#line 605 "TG01FD.f"
	    ir1 = *ranke + 1;
#line 606 "TG01FD.f"
	    mb03oy_(&la22, &na22, &a[ir1 + ir1 * a_dim1], lda, &toldef, &
		    svlmax, rnka22, sval, &iwork[1], &dwork[1], &dwork[kw], 
		    info);

/*           Apply transformation on the rest of matrices. */

#line 612 "TG01FD.f"
	    if (*rnka22 > 0) {

/*              A <-- diag(I, Q2') * A */
/*              Workspace: need   MIN(L,N) + N; */
/*                         prefer MIN(L,N) + N*NB. */

#line 618 "TG01FD.f"
		dormqr_("Left", "Transpose", &la22, ranke, rnka22, &a[ir1 + 
			ir1 * a_dim1], lda, &dwork[1], &a[ir1 + a_dim1], lda, 
			&dwork[kw], &lwr, info, (ftnlen)4, (ftnlen)9);

/*              B <-- diag(I, Q2') * B */
/*              Workspace: need   MIN(L,N) + M; */
/*                         prefer MIN(L,N) + M*NB. */

#line 626 "TG01FD.f"
		if (withb) {
#line 626 "TG01FD.f"
		    dormqr_("Left", "Transpose", &la22, m, rnka22, &a[ir1 + 
			    ir1 * a_dim1], lda, &dwork[1], &b[ir1 + b_dim1], 
			    ldb, &dwork[kw], &lwr, info, (ftnlen)4, (ftnlen)9)
			    ;
#line 626 "TG01FD.f"
		}

/*              Q <-- Q * diag(I, Q2) */
/*              Workspace: need   MIN(L,N) + L; */
/*                         prefer MIN(L,N) + L*NB. */

#line 635 "TG01FD.f"
		if (ilq) {
#line 635 "TG01FD.f"
		    dormqr_("Right", "No transpose", l, &la22, rnka22, &a[ir1 
			    + ir1 * a_dim1], lda, &dwork[1], &q[ir1 * q_dim1 
			    + 1], ldq, &dwork[kw], &lwr, info, (ftnlen)5, (
			    ftnlen)12);
#line 635 "TG01FD.f"
		}

/*              Set lower triangle of A22 to zero. */

#line 642 "TG01FD.f"
		if (la22 >= 2) {
#line 642 "TG01FD.f"
		    i__1 = la22 - 1;
#line 642 "TG01FD.f"
		    dlaset_("Lower", &i__1, rnka22, &c_b42, &c_b42, &a[ir1 + 
			    1 + ir1 * a_dim1], lda, (ftnlen)5);
#line 642 "TG01FD.f"
		}

/*              Compute A*diag(I,P2), C*diag(I,P2) and Z*diag(I,P2) */
/*              by forward permuting the columns of A, C and Z based */
/*              on information in IWORK. */

#line 650 "TG01FD.f"
		i__1 = na22;
#line 650 "TG01FD.f"
		for (j = 1; j <= i__1; ++j) {
#line 651 "TG01FD.f"
		    iwork[j] = -iwork[j];
#line 652 "TG01FD.f"
/* L40: */
#line 652 "TG01FD.f"
		}
#line 653 "TG01FD.f"
		i__1 = na22;
#line 653 "TG01FD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 654 "TG01FD.f"
		    if (iwork[i__] < 0) {
#line 655 "TG01FD.f"
			j = i__;
#line 656 "TG01FD.f"
			iwork[j] = -iwork[j];
#line 657 "TG01FD.f"
L50:
#line 658 "TG01FD.f"
			k = iwork[j];
#line 659 "TG01FD.f"
			if (iwork[k] < 0) {
#line 660 "TG01FD.f"
			    dswap_(ranke, &a[(*ranke + j) * a_dim1 + 1], &
				    c__1, &a[(*ranke + k) * a_dim1 + 1], &
				    c__1);
#line 662 "TG01FD.f"
			    if (withc) {
#line 662 "TG01FD.f"
				dswap_(p, &c__[(*ranke + j) * c_dim1 + 1], &
					c__1, &c__[(*ranke + k) * c_dim1 + 1],
					 &c__1);
#line 662 "TG01FD.f"
			    }
#line 665 "TG01FD.f"
			    if (ilz) {
#line 665 "TG01FD.f"
				dswap_(n, &z__[(*ranke + j) * z_dim1 + 1], &
					c__1, &z__[(*ranke + k) * z_dim1 + 1],
					 &c__1);
#line 665 "TG01FD.f"
			    }
#line 668 "TG01FD.f"
			    iwork[k] = -iwork[k];
#line 669 "TG01FD.f"
			    j = k;
#line 670 "TG01FD.f"
			    goto L50;
#line 671 "TG01FD.f"
			}
#line 672 "TG01FD.f"
		    }
#line 673 "TG01FD.f"
/* L60: */
#line 673 "TG01FD.f"
		}

#line 675 "TG01FD.f"
		if (reda && *rnka22 < na22) {

/*                 Determine an orthogonal matrix Y2 such that */

/*                 ( R11 R12 ) = ( Ar  0 ) * Y2 . */

/*                 Compute A <-- A*diag(I, Y2'), C <-- C*diag(I, Y2'), */
/*                         Z <-- Z*diag(I, Y2'). */
/*                 Workspace: need   2*N. */
/*                            prefer N + N*NB. */

#line 686 "TG01FD.f"
		    kw = *ranke + 1;
#line 687 "TG01FD.f"
		    i__1 = *ldwork - kw + 1;
#line 687 "TG01FD.f"
		    dtzrzf_(rnka22, &na22, &a[ir1 + ir1 * a_dim1], lda, &
			    dwork[1], &dwork[kw], &i__1, info);
/* Computing MAX */
#line 689 "TG01FD.f"
		    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 689 "TG01FD.f"
		    wrkopt = max(i__1,i__2);

/*                 Workspace: need   N + MAX(P,N); */
/*                            prefer N + MAX(P,N)*NB. */

#line 694 "TG01FD.f"
		    lh = na22 - *rnka22;
#line 695 "TG01FD.f"
		    if (withc) {
#line 696 "TG01FD.f"
			i__1 = *ldwork - kw + 1;
#line 696 "TG01FD.f"
			dormrz_("Right", "Transpose", p, n, rnka22, &lh, &a[
				ir1 + ir1 * a_dim1], lda, &dwork[1], &c__[
				c_offset], ldc, &dwork[kw], &i__1, info, (
				ftnlen)5, (ftnlen)9);
/* Computing MAX */
#line 699 "TG01FD.f"
			i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 699 "TG01FD.f"
			wrkopt = max(i__1,i__2);
#line 700 "TG01FD.f"
		    }
#line 701 "TG01FD.f"
		    if (ilz) {
#line 702 "TG01FD.f"
			i__1 = *ldwork - kw + 1;
#line 702 "TG01FD.f"
			dormrz_("Right", "Transpose", n, n, rnka22, &lh, &a[
				ir1 + ir1 * a_dim1], lda, &dwork[1], &z__[
				z_offset], ldz, &dwork[kw], &i__1, info, (
				ftnlen)5, (ftnlen)9);
/* Computing MAX */
#line 705 "TG01FD.f"
			i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 705 "TG01FD.f"
			wrkopt = max(i__1,i__2);
#line 706 "TG01FD.f"
		    }
#line 707 "TG01FD.f"
		    ire1 = *ranke + *rnka22 + 1;

/*                 Set R12 and R22 to zero. */

#line 711 "TG01FD.f"
		    dlaset_("Full", &la22, &lh, &c_b42, &c_b42, &a[ir1 + ire1 
			    * a_dim1], lda, (ftnlen)4);
#line 713 "TG01FD.f"
		}
#line 714 "TG01FD.f"
	    } else {
#line 715 "TG01FD.f"
		dlaset_("Full", &la22, &na22, &c_b42, &c_b42, &a[ir1 + ir1 * 
			a_dim1], lda, (ftnlen)4);
#line 717 "TG01FD.f"
	    }
#line 718 "TG01FD.f"
	}
#line 719 "TG01FD.f"
    }

#line 721 "TG01FD.f"
    dwork[1] = (doublereal) wrkopt;

#line 723 "TG01FD.f"
    return 0;
/* *** Last line of TG01FD *** */
} /* tg01fd_ */

