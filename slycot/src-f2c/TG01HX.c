#line 1 "TG01HX.f"
/* TG01HX.f -- translated by f2c (version 20100827).
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

#line 1 "TG01HX.f"
/* Table of constant values */

static doublereal c_b10 = 0.;
static doublereal c_b11 = 1.;
static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int tg01hx_(char *compq, char *compz, integer *l, integer *n,
	 integer *m, integer *p, integer *n1, integer *lbe, doublereal *a, 
	integer *lda, doublereal *e, integer *lde, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *q, integer *ldq, 
	doublereal *z__, integer *ldz, integer *nr, integer *nrblck, integer *
	rtau, doublereal *tol, integer *iwork, doublereal *dwork, integer *
	info, ftnlen compq_len, ftnlen compz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, e_dim1, 
	    e_offset, q_dim1, q_offset, z_dim1, z_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal t, c1, c2, s1, s2;
    static integer ic;
    static doublereal co;
    static integer nf, mn;
    static doublereal si, tt;
    static integer nr1;
    static logical ilq, ilz;
    static integer icol, rank;
    static doublereal smin, smax;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer irow;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal rcond;
    static logical withc;
    static integer ismin, ismax;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dlaic1_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static integer tauim1;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), dlaset_(char *, integer *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    static integer icompq, icompz;
    static doublereal sminpr, smaxpr, svlmax;


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

/*     Given the descriptor system (A-lambda*E,B,C) with the system */
/*     matrices A, E and B of the form */

/*            ( A1 X1 )        ( E1 Y1 )        ( B1 ) */
/*        A = (       ) ,  E = (       ) ,  B = (    ) , */
/*            ( 0  X2 )        ( 0  Y2 )        ( 0  ) */

/*     where */
/*          - B is an L-by-M matrix, with B1 an N1-by-M  submatrix */
/*          - A is an L-by-N matrix, with A1 an N1-by-N1 submatrix */
/*          - E is an L-by-N matrix, with E1 an N1-by-N1 submatrix */
/*              with LBE nonzero sub-diagonals, */
/*     this routine reduces the pair (A1-lambda*E1,B1) to the form */

/*     Qc'*[A1-lambda*E1 B1]*diag(Zc,I) = */

/*                              ( Bc Ac-lambda*Ec      *         ) */
/*                              (                                ) , */
/*                              ( 0     0         Anc-lambda*Enc ) */

/*     where: */
/*     1) the pencil ( Bc Ac-lambda*Ec ) has full row rank NR for */
/*        all finite lambda and is in a staircase form with */
/*                           _      _          _        _ */
/*                         ( A1,0   A1,1  ...  A1,k-1   A1,k  ) */
/*                         (        _          _        _     ) */
/*             ( Bc Ac ) = (  0     A2,1  ...  A2,k-1   A2,k  ) ,  (1) */
/*                         (              ...  _        _     ) */
/*                         (  0       0   ...  Ak,k-1   Ak,k  ) */

/*                           _          _        _ */
/*                         ( E1,1  ...  E1,k-1   E1,k  ) */
/*                         (            _        _     ) */
/*               Ec      = (   0   ...  E2,k-1   E2,k  ) ,         (2) */
/*                         (       ...           _     ) */
/*                         (   0   ...    0      Ek,k  ) */
/*               _ */
/*         where Ai,i-1 is an rtau(i)-by-rtau(i-1) full row rank */
/*                                       _ */
/*         matrix (with rtau(0) = M) and Ei,i is an rtau(i)-by-rtau(i) */
/*         upper triangular matrix. */

/*      2) the pencil Anc-lambda*Enc is regular of order N1-NR with Enc */
/*         upper triangular; this pencil contains the uncontrollable */
/*         finite eigenvalues of the pencil (A1-lambda*E1). */

/*     The transformations are applied to the whole matrices A, E, B */
/*     and C. The left and/or right orthogonal transformations Qc and Zc */
/*     performed to reduce the pencil S(lambda) can be optionally */
/*     accumulated in the matrices Q and Z, respectivelly. */

/*     The reduced order descriptor system (Ac-lambda*Ec,Bc,Cc) has no */
/*     uncontrollable finite eigenvalues and has the same */
/*     transfer-function matrix as the original system (A-lambda*E,B,C). */

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

/*     Input/Output Parameters */

/*     L       (input) INTEGER */
/*             The number of descriptor state equations; also the number */
/*             of rows of matrices A, E and B.  L >= 0. */

/*     N       (input) INTEGER */
/*             The dimension of the descriptor state vector; also the */
/*             number of columns of matrices A, E and C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The dimension of descriptor system input vector; also the */
/*             number of columns of matrix B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The dimension of descriptor system output; also the */
/*             number of rows of matrix C.  P >= 0. */

/*     N1      (input) INTEGER */
/*             The order of subsystem (A1-lambda*E1,B1,C1) to be reduced. */
/*             MIN(L,N) >= N1 >= 0. */

/*     LBE     (input) INTEGER */
/*             The number of nonzero sub-diagonals of submatrix E1. */
/*             MAX(0,N1-1) >= LBE >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the L-by-N state matrix A in the partitioned */
/*             form */
/*                      ( A1 X1 ) */
/*                  A = (       ) , */
/*                      ( 0  X2 ) */

/*             where A1 is N1-by-N1. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the transformed state matrix, */

/*                                  ( Ac  *   * ) */
/*                       Qc'*A*Zc = ( 0  Anc  * ) , */
/*                                  ( 0   0   * ) */

/*             where Ac is NR-by-NR and Anc is (N1-NR)-by-(N1-NR). */
/*             The matrix ( Bc Ac ) is in the controlability */
/*             staircase form (1). */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,L). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the L-by-N descriptor matrix E in the partitioned */
/*             form */
/*                      ( E1 Y1 ) */
/*                  E = (       ) , */
/*                      ( 0  Y2 ) */

/*             where E1 is N1-by-N1 matrix with LBE nonzero */
/*             sub-diagonals. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the transformed descriptor matrix */

/*                                  ( Ec  *   * ) */
/*                       Qc'*E*Zc = ( 0  Enc  * ) , */
/*                                  ( 0   0   * ) */

/*             where Ec is NR-by-NR and Enc is (N1-NR)-by-(N1-NR). */
/*             Both Ec and Enc are upper triangular and Enc is */
/*             nonsingular. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,L). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading L-by-M part of this array must */
/*             contain the L-by-M input matrix B in the partitioned */
/*             form */
/*                      ( B1 ) */
/*                  B = (    ) , */
/*                      ( 0  ) */

/*             where B1 is N1-by-M. */
/*             On exit, the leading L-by-M part of this array contains */
/*             the transformed input matrix */

/*                               ( Bc ) */
/*                       Qc'*B = (    ) , */
/*                               ( 0  ) */

/*             where Bc is NR-by-M. */
/*             The matrix ( Bc Ac ) is in the controlability */
/*             staircase form (1). */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,L). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the transformed matrix C*Zc. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,L) */
/*             If COMPQ = 'N': Q is not referenced. */
/*             If COMPQ = 'I': on entry, Q need not be set; */
/*                             on exit, the leading L-by-L part of this */
/*                             array contains the orthogonal matrix Q, */
/*                             where Q' is the product of transformations */
/*                             which are applied to A, E, and B on */
/*                             the left. */
/*             If COMPQ = 'U': on entry, the leading L-by-L part of this */
/*                             array must contain an orthogonal matrix */
/*                             Qc; */
/*                             on exit, the leading L-by-L part of this */
/*                             array contains the orthogonal matrix */
/*                             Qc*Q. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. */
/*             LDQ >= 1,        if COMPQ = 'N'; */
/*             LDQ >= MAX(1,L), if COMPQ = 'U' or 'I'. */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             If COMPZ = 'N': Z is not referenced. */
/*             If COMPZ = 'I': on entry, Z need not be set; */
/*                             on exit, the leading N-by-N part of this */
/*                             array contains the orthogonal matrix Z, */
/*                             which is the product of transformations */
/*                             applied to A, E, and C on the right. */
/*             If COMPZ = 'U': on entry, the leading N-by-N part of this */
/*                             array must contain an orthogonal matrix */
/*                             Zc; */
/*                             on exit, the leading N-by-N part of this */
/*                             array contains the orthogonal matrix */
/*                             Zc*Z. */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z. */
/*             LDZ >= 1,        if COMPZ = 'N'; */
/*             LDZ >= MAX(1,N), if COMPZ = 'U' or 'I'. */

/*     NR      (output) INTEGER */
/*             The order of the reduced matrices Ac and Ec, and the */
/*             number of rows of the reduced matrix Bc; also the order of */
/*             the controllable part of the pair (B, A-lambda*E). */

/*     NRBLCK  (output) INTEGER                      _ */
/*             The number k, of full row rank blocks Ai,i in the */
/*             staircase form of the pencil (Bc Ac-lambda*Ec) (see (1) */
/*             and (2)). */

/*     RTAU    (output) INTEGER array, dimension (N1) */
/*             RTAU(i), for i = 1, ..., NRBLCK, is the row dimension of */
/*                                     _ */
/*             the full row rank block Ai,i-1 in the staircase form (1). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in rank determinations when */
/*             transforming (A-lambda*E, B). If the user sets TOL > 0, */
/*             then the given value of TOL is used as a lower bound for */
/*             reciprocal condition numbers in rank determinations; a */
/*             (sub)matrix whose estimated condition number is less than */
/*             1/TOL is considered to be of full rank.  If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance, */
/*             defined by  TOLDEF = L*N*EPS,  is used instead, where */
/*             EPS is the machine precision (see LAPACK Library routine */
/*             DLAMCH).  TOL < 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M) */

/*     DWORK   DOUBLE PRECISION array, dimension MAX(N,L,2*M) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The subroutine is based on the reduction algorithm of [1]. */

/*     REFERENCES */

/*     [1] A. Varga */
/*         Computation of Irreducible Generalized State-Space */
/*         Realizations. */
/*         Kybernetika, vol. 26, pp. 89-106, 1990. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically backward stable and requires */
/*     0( N*N1**2 )  floating point operations. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     March 1999. Based on the RASP routine RPDS05. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, July 1999, */
/*     May 2003, Nov. 2003. */
/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, Nov. 2003. */

/*     KEYWORDS */

/*     Controllability, minimal realization, orthogonal canonical form, */
/*     orthogonal transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Decode COMPQ. */

#line 354 "TG01HX.f"
    /* Parameter adjustments */
#line 354 "TG01HX.f"
    a_dim1 = *lda;
#line 354 "TG01HX.f"
    a_offset = 1 + a_dim1;
#line 354 "TG01HX.f"
    a -= a_offset;
#line 354 "TG01HX.f"
    e_dim1 = *lde;
#line 354 "TG01HX.f"
    e_offset = 1 + e_dim1;
#line 354 "TG01HX.f"
    e -= e_offset;
#line 354 "TG01HX.f"
    b_dim1 = *ldb;
#line 354 "TG01HX.f"
    b_offset = 1 + b_dim1;
#line 354 "TG01HX.f"
    b -= b_offset;
#line 354 "TG01HX.f"
    c_dim1 = *ldc;
#line 354 "TG01HX.f"
    c_offset = 1 + c_dim1;
#line 354 "TG01HX.f"
    c__ -= c_offset;
#line 354 "TG01HX.f"
    q_dim1 = *ldq;
#line 354 "TG01HX.f"
    q_offset = 1 + q_dim1;
#line 354 "TG01HX.f"
    q -= q_offset;
#line 354 "TG01HX.f"
    z_dim1 = *ldz;
#line 354 "TG01HX.f"
    z_offset = 1 + z_dim1;
#line 354 "TG01HX.f"
    z__ -= z_offset;
#line 354 "TG01HX.f"
    --rtau;
#line 354 "TG01HX.f"
    --iwork;
#line 354 "TG01HX.f"
    --dwork;
#line 354 "TG01HX.f"

#line 354 "TG01HX.f"
    /* Function Body */
#line 354 "TG01HX.f"
    if (lsame_(compq, "N", (ftnlen)1, (ftnlen)1)) {
#line 355 "TG01HX.f"
	ilq = FALSE_;
#line 356 "TG01HX.f"
	icompq = 1;
#line 357 "TG01HX.f"
    } else if (lsame_(compq, "U", (ftnlen)1, (ftnlen)1)) {
#line 358 "TG01HX.f"
	ilq = TRUE_;
#line 359 "TG01HX.f"
	icompq = 2;
#line 360 "TG01HX.f"
    } else if (lsame_(compq, "I", (ftnlen)1, (ftnlen)1)) {
#line 361 "TG01HX.f"
	ilq = TRUE_;
#line 362 "TG01HX.f"
	icompq = 3;
#line 363 "TG01HX.f"
    } else {
#line 364 "TG01HX.f"
	icompq = 0;
#line 365 "TG01HX.f"
    }

/*     Decode COMPZ. */

#line 369 "TG01HX.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 370 "TG01HX.f"
	ilz = FALSE_;
#line 371 "TG01HX.f"
	icompz = 1;
#line 372 "TG01HX.f"
    } else if (lsame_(compz, "U", (ftnlen)1, (ftnlen)1)) {
#line 373 "TG01HX.f"
	ilz = TRUE_;
#line 374 "TG01HX.f"
	icompz = 2;
#line 375 "TG01HX.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 376 "TG01HX.f"
	ilz = TRUE_;
#line 377 "TG01HX.f"
	icompz = 3;
#line 378 "TG01HX.f"
    } else {
#line 379 "TG01HX.f"
	icompz = 0;
#line 380 "TG01HX.f"
    }

/*     Test the input scalar parameters. */

#line 384 "TG01HX.f"
    *info = 0;
#line 385 "TG01HX.f"
    if (icompq <= 0) {
#line 386 "TG01HX.f"
	*info = -1;
#line 387 "TG01HX.f"
    } else if (icompz <= 0) {
#line 388 "TG01HX.f"
	*info = -2;
#line 389 "TG01HX.f"
    } else if (*l < 0) {
#line 390 "TG01HX.f"
	*info = -3;
#line 391 "TG01HX.f"
    } else if (*n < 0) {
#line 392 "TG01HX.f"
	*info = -4;
#line 393 "TG01HX.f"
    } else if (*m < 0) {
#line 394 "TG01HX.f"
	*info = -5;
#line 395 "TG01HX.f"
    } else if (*p < 0) {
#line 396 "TG01HX.f"
	*info = -6;
#line 397 "TG01HX.f"
    } else if (*n1 < 0 || *n1 > min(*l,*n)) {
#line 398 "TG01HX.f"
	*info = -7;
#line 399 "TG01HX.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 399 "TG01HX.f"
	i__1 = 0, i__2 = *n1 - 1;
#line 399 "TG01HX.f"
	if (*lbe < 0 || *lbe > max(i__1,i__2)) {
#line 400 "TG01HX.f"
	    *info = -8;
#line 401 "TG01HX.f"
	} else if (*lda < max(1,*l)) {
#line 402 "TG01HX.f"
	    *info = -10;
#line 403 "TG01HX.f"
	} else if (*lde < max(1,*l)) {
#line 404 "TG01HX.f"
	    *info = -12;
#line 405 "TG01HX.f"
	} else if (*ldb < max(1,*l)) {
#line 406 "TG01HX.f"
	    *info = -14;
#line 407 "TG01HX.f"
	} else if (*ldc < max(1,*p)) {
#line 408 "TG01HX.f"
	    *info = -16;
#line 409 "TG01HX.f"
	} else if (ilq && *ldq < *l || *ldq < 1) {
#line 410 "TG01HX.f"
	    *info = -18;
#line 411 "TG01HX.f"
	} else if (ilz && *ldz < *n || *ldz < 1) {
#line 412 "TG01HX.f"
	    *info = -20;
#line 413 "TG01HX.f"
	} else if (*tol >= 1.) {
#line 414 "TG01HX.f"
	    *info = -24;
#line 415 "TG01HX.f"
	}
#line 415 "TG01HX.f"
    }
#line 416 "TG01HX.f"
    if (*info != 0) {
#line 417 "TG01HX.f"
	i__1 = -(*info);
#line 417 "TG01HX.f"
	xerbla_("TG01HX", &i__1, (ftnlen)6);
#line 418 "TG01HX.f"
	return 0;
#line 419 "TG01HX.f"
    }

/*     Initialize Q and Z if necessary. */

#line 423 "TG01HX.f"
    if (icompq == 3) {
#line 423 "TG01HX.f"
	dlaset_("Full", l, l, &c_b10, &c_b11, &q[q_offset], ldq, (ftnlen)4);
#line 423 "TG01HX.f"
    }
#line 425 "TG01HX.f"
    if (icompz == 3) {
#line 425 "TG01HX.f"
	dlaset_("Full", n, n, &c_b10, &c_b11, &z__[z_offset], ldz, (ftnlen)4);
#line 425 "TG01HX.f"
    }

/*     Initialize output variables. */

#line 430 "TG01HX.f"
    *nr = 0;
#line 431 "TG01HX.f"
    *nrblck = 0;

/*     Quick return if possible. */

#line 435 "TG01HX.f"
    if (*m == 0 || *n1 == 0) {
#line 436 "TG01HX.f"
	return 0;
#line 437 "TG01HX.f"
    }

#line 439 "TG01HX.f"
    withc = *p > 0;
#line 440 "TG01HX.f"
    d__1 = dlange_("F", l, m, &b[b_offset], ldb, &dwork[1], (ftnlen)1);
#line 440 "TG01HX.f"
    d__2 = dlange_("F", l, n, &a[a_offset], lda, &dwork[1], (ftnlen)1);
#line 440 "TG01HX.f"
    svlmax = dlapy2_(&d__1, &d__2);
#line 442 "TG01HX.f"
    rcond = *tol;
#line 443 "TG01HX.f"
    if (rcond <= 0.) {

/*        Use the default tolerance in controllability determination. */

#line 447 "TG01HX.f"
	rcond = (doublereal) (*l * *n) * dlamch_("EPSILON", (ftnlen)7);
#line 448 "TG01HX.f"
    }

#line 450 "TG01HX.f"
    if (svlmax < rcond) {
#line 450 "TG01HX.f"
	svlmax = 1.;
#line 450 "TG01HX.f"
    }

/*     Reduce E to upper triangular form if necessary. */

#line 455 "TG01HX.f"
    if (*lbe > 0) {
#line 456 "TG01HX.f"
	i__1 = *n1 - 1;
#line 456 "TG01HX.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Generate elementary reflector H(i) to annihilate */
/*           E(i+1:i+lbe,i). */

/* Computing MIN */
#line 461 "TG01HX.f"
	    i__2 = *lbe, i__3 = *n1 - i__;
#line 461 "TG01HX.f"
	    k = min(i__2,i__3) + 1;
#line 462 "TG01HX.f"
	    dlarfg_(&k, &e[i__ + i__ * e_dim1], &e[i__ + 1 + i__ * e_dim1], &
		    c__1, &tt);
#line 463 "TG01HX.f"
	    t = e[i__ + i__ * e_dim1];
#line 464 "TG01HX.f"
	    e[i__ + i__ * e_dim1] = 1.;

/*           Apply H(i) to E(i:n1,i+1:n) from the left. */

#line 468 "TG01HX.f"
	    i__2 = *n - i__;
#line 468 "TG01HX.f"
	    dlarf_("Left", &k, &i__2, &e[i__ + i__ * e_dim1], &c__1, &tt, &e[
		    i__ + (i__ + 1) * e_dim1], lde, &dwork[1], (ftnlen)4);

/*           Apply H(i) to A(i:n1,1:n) from the left. */

#line 473 "TG01HX.f"
	    dlarf_("Left", &k, n, &e[i__ + i__ * e_dim1], &c__1, &tt, &a[i__ 
		    + a_dim1], lda, &dwork[1], (ftnlen)4);

/*           Apply H(i) to B(i:n1,1:m) from the left. */

#line 478 "TG01HX.f"
	    dlarf_("Left", &k, m, &e[i__ + i__ * e_dim1], &c__1, &tt, &b[i__ 
		    + b_dim1], ldb, &dwork[1], (ftnlen)4);
#line 480 "TG01HX.f"
	    if (ilq) {

/*              Apply H(i) to Q(1:l,i:n1) from the right. */

#line 484 "TG01HX.f"
		dlarf_("Right", l, &k, &e[i__ + i__ * e_dim1], &c__1, &tt, &q[
			i__ * q_dim1 + 1], ldq, &dwork[1], (ftnlen)5);
#line 486 "TG01HX.f"
	    }
#line 487 "TG01HX.f"
	    e[i__ + i__ * e_dim1] = t;
#line 488 "TG01HX.f"
/* L10: */
#line 488 "TG01HX.f"
	}
#line 489 "TG01HX.f"
	if (*n1 > 1) {
#line 489 "TG01HX.f"
	    i__1 = *n1 - 1;
#line 489 "TG01HX.f"
	    i__2 = *n1 - 1;
#line 489 "TG01HX.f"
	    dlaset_("Lower", &i__1, &i__2, &c_b10, &c_b10, &e[e_dim1 + 2], 
		    lde, (ftnlen)5);
#line 489 "TG01HX.f"
	}
#line 491 "TG01HX.f"
    }

#line 493 "TG01HX.f"
    ismin = 1;
#line 494 "TG01HX.f"
    ismax = ismin + *m;
#line 495 "TG01HX.f"
    ic = -(*m);
#line 496 "TG01HX.f"
    tauim1 = *m;
#line 497 "TG01HX.f"
    nf = *n1;

#line 499 "TG01HX.f"
L20:
#line 500 "TG01HX.f"
    ++(*nrblck);
#line 501 "TG01HX.f"
    rank = 0;
#line 502 "TG01HX.f"
    if (nf > 0) {

/*        IROW will point to the current pivot line in B, */
/*        ICOL+1 will point to the first active columns of A. */

#line 507 "TG01HX.f"
	icol = ic + tauim1;
#line 508 "TG01HX.f"
	irow = *nr;
#line 509 "TG01HX.f"
	nr1 = *nr + 1;
#line 510 "TG01HX.f"
	if (*nr > 0) {
#line 510 "TG01HX.f"
	    dlacpy_("Full", &nf, &tauim1, &a[nr1 + (ic + 1) * a_dim1], lda, &
		    b[nr1 + b_dim1], ldb, (ftnlen)4);
#line 510 "TG01HX.f"
	}

/*        Perform QR-decomposition with column pivoting on the current B */
/*        while keeping E upper triangular. */
/*        The current B is at first iteration B and for subsequent */
/*        iterations the NF-by-TAUIM1 matrix delimited by rows */
/*        NR + 1 to N1 and columns IC + 1 to IC + TAUIM1 of A. */
/*        The rank of current B is computed in RANK. */

#line 521 "TG01HX.f"
	if (tauim1 > 1) {

/*           Compute column norms. */

#line 525 "TG01HX.f"
	    i__1 = tauim1;
#line 525 "TG01HX.f"
	    for (j = 1; j <= i__1; ++j) {
#line 526 "TG01HX.f"
		dwork[j] = dnrm2_(&nf, &b[nr1 + j * b_dim1], &c__1);
#line 527 "TG01HX.f"
		dwork[*m + j] = dwork[j];
#line 528 "TG01HX.f"
		iwork[j] = j;
#line 529 "TG01HX.f"
/* L30: */
#line 529 "TG01HX.f"
	    }
#line 530 "TG01HX.f"
	}

#line 532 "TG01HX.f"
	mn = min(nf,tauim1);

#line 534 "TG01HX.f"
L40:
#line 535 "TG01HX.f"
	if (rank < mn) {
#line 536 "TG01HX.f"
	    j = rank + 1;
#line 537 "TG01HX.f"
	    ++irow;

/*           Pivot if necessary. */

#line 541 "TG01HX.f"
	    if (j != tauim1) {
#line 542 "TG01HX.f"
		i__1 = tauim1 - j + 1;
#line 542 "TG01HX.f"
		k = j - 1 + idamax_(&i__1, &dwork[j], &c__1);
#line 543 "TG01HX.f"
		if (k != j) {
#line 544 "TG01HX.f"
		    dswap_(&nf, &b[nr1 + j * b_dim1], &c__1, &b[nr1 + k * 
			    b_dim1], &c__1);
#line 545 "TG01HX.f"
		    i__ = iwork[k];
#line 546 "TG01HX.f"
		    iwork[k] = iwork[j];
#line 547 "TG01HX.f"
		    iwork[j] = i__;
#line 548 "TG01HX.f"
		    dwork[k] = dwork[j];
#line 549 "TG01HX.f"
		    dwork[*m + k] = dwork[*m + j];
#line 550 "TG01HX.f"
		}
#line 551 "TG01HX.f"
	    }

/*           Zero elements below the current diagonal element of B. */

#line 555 "TG01HX.f"
	    i__1 = irow;
#line 555 "TG01HX.f"
	    for (i__ = *n1 - 1; i__ >= i__1; --i__) {

/*              Rotate rows I and I+1 to zero B(I+1,J). */

#line 559 "TG01HX.f"
		t = b[i__ + j * b_dim1];
#line 560 "TG01HX.f"
		dlartg_(&t, &b[i__ + 1 + j * b_dim1], &co, &si, &b[i__ + j * 
			b_dim1]);
#line 561 "TG01HX.f"
		b[i__ + 1 + j * b_dim1] = 0.;
#line 562 "TG01HX.f"
		i__2 = *n - i__ + 1;
#line 562 "TG01HX.f"
		drot_(&i__2, &e[i__ + i__ * e_dim1], lde, &e[i__ + 1 + i__ * 
			e_dim1], lde, &co, &si);
#line 563 "TG01HX.f"
		if (j < tauim1) {
#line 563 "TG01HX.f"
		    i__2 = tauim1 - j;
#line 563 "TG01HX.f"
		    drot_(&i__2, &b[i__ + (j + 1) * b_dim1], ldb, &b[i__ + 1 
			    + (j + 1) * b_dim1], ldb, &co, &si);
#line 563 "TG01HX.f"
		}
#line 566 "TG01HX.f"
		i__2 = *n - icol;
#line 566 "TG01HX.f"
		drot_(&i__2, &a[i__ + (icol + 1) * a_dim1], lda, &a[i__ + 1 + 
			(icol + 1) * a_dim1], lda, &co, &si);
#line 568 "TG01HX.f"
		if (ilq) {
#line 568 "TG01HX.f"
		    drot_(l, &q[i__ * q_dim1 + 1], &c__1, &q[(i__ + 1) * 
			    q_dim1 + 1], &c__1, &co, &si);
#line 568 "TG01HX.f"
		}

/*              Rotate columns I, I+1 to zero E(I+1,I). */

#line 572 "TG01HX.f"
		t = e[i__ + 1 + (i__ + 1) * e_dim1];
#line 573 "TG01HX.f"
		dlartg_(&t, &e[i__ + 1 + i__ * e_dim1], &co, &si, &e[i__ + 1 
			+ (i__ + 1) * e_dim1]);
#line 574 "TG01HX.f"
		e[i__ + 1 + i__ * e_dim1] = 0.;
#line 575 "TG01HX.f"
		drot_(&i__, &e[(i__ + 1) * e_dim1 + 1], &c__1, &e[i__ * 
			e_dim1 + 1], &c__1, &co, &si);
#line 576 "TG01HX.f"
		drot_(n1, &a[(i__ + 1) * a_dim1 + 1], &c__1, &a[i__ * a_dim1 
			+ 1], &c__1, &co, &si);
#line 577 "TG01HX.f"
		if (ilz) {
#line 577 "TG01HX.f"
		    drot_(n, &z__[(i__ + 1) * z_dim1 + 1], &c__1, &z__[i__ * 
			    z_dim1 + 1], &c__1, &co, &si);
#line 577 "TG01HX.f"
		}
#line 578 "TG01HX.f"
		if (withc) {
#line 578 "TG01HX.f"
		    drot_(p, &c__[(i__ + 1) * c_dim1 + 1], &c__1, &c__[i__ * 
			    c_dim1 + 1], &c__1, &co, &si);
#line 578 "TG01HX.f"
		}
#line 580 "TG01HX.f"
/* L50: */
#line 580 "TG01HX.f"
	    }

#line 582 "TG01HX.f"
	    if (rank == 0) {

/*              Initialize; exit if matrix is zero (RANK = 0). */

#line 586 "TG01HX.f"
		smax = (d__1 = b[nr1 + b_dim1], abs(d__1));
#line 587 "TG01HX.f"
		if (smax == 0.) {
#line 587 "TG01HX.f"
		    goto L80;
#line 587 "TG01HX.f"
		}
#line 588 "TG01HX.f"
		smin = smax;
#line 589 "TG01HX.f"
		smaxpr = smax;
#line 590 "TG01HX.f"
		sminpr = smin;
#line 591 "TG01HX.f"
		c1 = 1.;
#line 592 "TG01HX.f"
		c2 = 1.;
#line 593 "TG01HX.f"
	    } else {

/*              One step of incremental condition estimation. */

#line 597 "TG01HX.f"
		dlaic1_(&c__2, &rank, &dwork[ismin], &smin, &b[nr1 + j * 
			b_dim1], &b[irow + j * b_dim1], &sminpr, &s1, &c1);
#line 599 "TG01HX.f"
		dlaic1_(&c__1, &rank, &dwork[ismax], &smax, &b[nr1 + j * 
			b_dim1], &b[irow + j * b_dim1], &smaxpr, &s2, &c2);
#line 601 "TG01HX.f"
	    }

/*           Check the rank; finish the loop if rank loss occurs. */

#line 605 "TG01HX.f"
	    if (svlmax * rcond <= smaxpr) {
#line 606 "TG01HX.f"
		if (svlmax * rcond <= sminpr) {
#line 607 "TG01HX.f"
		    if (smaxpr * rcond <= sminpr) {

/*                    Finish the loop if last row. */

#line 611 "TG01HX.f"
			if (irow == *n1) {
#line 612 "TG01HX.f"
			    ++rank;
#line 613 "TG01HX.f"
			    goto L80;
#line 614 "TG01HX.f"
			}

/*                    Update partial column norms. */

#line 618 "TG01HX.f"
			i__1 = tauim1;
#line 618 "TG01HX.f"
			for (i__ = j + 1; i__ <= i__1; ++i__) {
#line 619 "TG01HX.f"
			    if (dwork[i__] != 0.) {
/* Computing 2nd power */
#line 620 "TG01HX.f"
				d__2 = (d__1 = b[irow + i__ * b_dim1], abs(
					d__1)) / dwork[i__];
#line 620 "TG01HX.f"
				t = 1. - d__2 * d__2;
#line 621 "TG01HX.f"
				t = max(t,0.);
/* Computing 2nd power */
#line 622 "TG01HX.f"
				d__1 = dwork[i__] / dwork[*m + i__];
#line 622 "TG01HX.f"
				tt = t * .05 * (d__1 * d__1) + 1.;
#line 623 "TG01HX.f"
				if (tt != 1.) {
#line 624 "TG01HX.f"
				    dwork[i__] *= sqrt(t);
#line 625 "TG01HX.f"
				} else {
#line 626 "TG01HX.f"
				    i__2 = nf - j;
#line 626 "TG01HX.f"
				    dwork[i__] = dnrm2_(&i__2, &b[irow + 1 + 
					    i__ * b_dim1], &c__1);
#line 627 "TG01HX.f"
				    dwork[*m + i__] = dwork[i__];
#line 628 "TG01HX.f"
				}
#line 629 "TG01HX.f"
			    }
#line 630 "TG01HX.f"
/* L60: */
#line 630 "TG01HX.f"
			}

#line 632 "TG01HX.f"
			i__1 = rank;
#line 632 "TG01HX.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 633 "TG01HX.f"
			    dwork[ismin + i__ - 1] = s1 * dwork[ismin + i__ - 
				    1];
#line 634 "TG01HX.f"
			    dwork[ismax + i__ - 1] = s2 * dwork[ismax + i__ - 
				    1];
#line 635 "TG01HX.f"
/* L70: */
#line 635 "TG01HX.f"
			}

#line 637 "TG01HX.f"
			dwork[ismin + rank] = c1;
#line 638 "TG01HX.f"
			dwork[ismax + rank] = c2;
#line 639 "TG01HX.f"
			smin = sminpr;
#line 640 "TG01HX.f"
			smax = smaxpr;
#line 641 "TG01HX.f"
			++rank;
#line 642 "TG01HX.f"
			goto L40;
#line 643 "TG01HX.f"
		    }
#line 644 "TG01HX.f"
		}
#line 645 "TG01HX.f"
	    }
#line 646 "TG01HX.f"
	    if (*nr > 0) {
#line 647 "TG01HX.f"
		i__1 = *n1 - irow + 1;
#line 647 "TG01HX.f"
		i__2 = tauim1 - j + 1;
#line 647 "TG01HX.f"
		dlaset_("Full", &i__1, &i__2, &c_b10, &c_b10, &b[irow + j * 
			b_dim1], ldb, (ftnlen)4);
#line 649 "TG01HX.f"
	    }
#line 650 "TG01HX.f"
	    goto L80;
#line 651 "TG01HX.f"
	}
#line 652 "TG01HX.f"
    }

#line 654 "TG01HX.f"
L80:
#line 654 "TG01HX.f"
    if (rank > 0) {
#line 655 "TG01HX.f"
	rtau[*nrblck] = rank;

/*        Back permute interchanged columns. */

#line 659 "TG01HX.f"
	if (tauim1 > 1) {
#line 660 "TG01HX.f"
	    i__1 = tauim1;
#line 660 "TG01HX.f"
	    for (j = 1; j <= i__1; ++j) {
#line 661 "TG01HX.f"
		if (iwork[j] > 0) {
#line 662 "TG01HX.f"
		    k = iwork[j];
#line 663 "TG01HX.f"
		    iwork[j] = -k;
#line 664 "TG01HX.f"
L90:
#line 665 "TG01HX.f"
		    if (k != j) {
#line 666 "TG01HX.f"
			dswap_(&rank, &b[nr1 + j * b_dim1], &c__1, &b[nr1 + k 
				* b_dim1], &c__1);
#line 667 "TG01HX.f"
			iwork[k] = -iwork[k];
#line 668 "TG01HX.f"
			k = -iwork[k];
#line 669 "TG01HX.f"
			goto L90;
#line 670 "TG01HX.f"
		    }
#line 671 "TG01HX.f"
		}
#line 672 "TG01HX.f"
/* L100: */
#line 672 "TG01HX.f"
	    }
#line 673 "TG01HX.f"
	}
#line 674 "TG01HX.f"
    }
#line 675 "TG01HX.f"
    if (*nr > 0) {
#line 675 "TG01HX.f"
	dlacpy_("Full", &nf, &tauim1, &b[nr1 + b_dim1], ldb, &a[nr1 + (ic + 1)
		 * a_dim1], lda, (ftnlen)4);
#line 675 "TG01HX.f"
    }
#line 678 "TG01HX.f"
    if (rank > 0) {
#line 679 "TG01HX.f"
	*nr += rank;
#line 680 "TG01HX.f"
	nf -= rank;
#line 681 "TG01HX.f"
	ic += tauim1;
#line 682 "TG01HX.f"
	tauim1 = rank;
#line 683 "TG01HX.f"
	goto L20;
#line 684 "TG01HX.f"
    } else {
#line 685 "TG01HX.f"
	--(*nrblck);
#line 686 "TG01HX.f"
    }

#line 688 "TG01HX.f"
    if (*nrblck > 0) {
#line 688 "TG01HX.f"
	rank = rtau[1];
#line 688 "TG01HX.f"
    }
#line 689 "TG01HX.f"
    if (rank < *n1) {
#line 689 "TG01HX.f"
	i__1 = *n1 - rank;
#line 689 "TG01HX.f"
	dlaset_("Full", &i__1, m, &c_b10, &c_b10, &b[rank + 1 + b_dim1], ldb, 
		(ftnlen)4);
#line 689 "TG01HX.f"
    }

#line 692 "TG01HX.f"
    return 0;
/* *** Last line of TG01HX *** */
} /* tg01hx_ */

