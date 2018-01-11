#line 1 "TG01HD.f"
/* TG01HD.f -- translated by f2c (version 20100827).
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

#line 1 "TG01HD.f"
/* Subroutine */ int tg01hd_(char *jobcon, char *compq, char *compz, integer *
	n, integer *m, integer *p, doublereal *a, integer *lda, doublereal *e,
	 integer *lde, doublereal *b, integer *ldb, doublereal *c__, integer *
	ldc, doublereal *q, integer *ldq, doublereal *z__, integer *ldz, 
	integer *ncont, integer *niucon, integer *nrblck, integer *rtau, 
	doublereal *tol, integer *iwork, doublereal *dwork, integer *info, 
	ftnlen jobcon_len, ftnlen compq_len, ftnlen compz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, e_dim1, 
	    e_offset, q_dim1, q_offset, z_dim1, z_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer nr, lba;
    static logical ilq, ilz;
    static char jobq[1], jobz[1];
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tg01hx_(char *, char *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen);
    static logical fincon, infcon;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer icompq, icompz;


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

/*     To compute orthogonal transformation matrices Q and Z which */
/*     reduce the N-th order descriptor system (A-lambda*E,B,C) */
/*     to the form */

/*                ( Ac  *  )             ( Ec  *  )           ( Bc ) */
/*       Q'*A*Z = (        ) ,  Q'*E*Z = (        ) ,  Q'*B = (    ) , */
/*                ( 0  Anc )             ( 0  Enc )           ( 0  ) */

/*          C*Z = ( Cc Cnc ) , */

/*     where the NCONT-th order descriptor system (Ac-lambda*Ec,Bc,Cc) */
/*     is a finite and/or infinite controllable. The pencil */
/*     Anc - lambda*Enc is regular of order N-NCONT and contains the */
/*     uncontrollable finite and/or infinite eigenvalues of the pencil */
/*     A-lambda*E. */

/*     For JOBCON = 'C' or 'I', the pencil ( Bc Ec-lambda*Ac ) has full */
/*     row rank NCONT for all finite lambda and is in a staircase form */
/*     with */
/*                     _      _          _        _ */
/*                   ( E1,0   E1,1  ...  E1,k-1   E1,k  ) */
/*                   (        _          _        _     ) */
/*       ( Bc Ec ) = (  0     E2,1  ...  E2,k-1   E2,k  ) ,  (1) */
/*                   (              ...  _        _     ) */
/*                   (  0       0   ...  Ek,k-1   Ek,k  ) */

/*                     _          _        _ */
/*                   ( A1,1  ...  A1,k-1   A1,k  ) */
/*                   (            _        _     ) */
/*         Ac      = (   0   ...  A2,k-1   A2,k  ) ,         (2) */
/*                   (       ...           _     ) */
/*                   (   0   ...    0      Ak,k  ) */
/*           _ */
/*     where Ei,i-1 is an rtau(i)-by-rtau(i-1) full row rank matrix */
/*                            _ */
/*     (with rtau(0) = M) and Ai,i is an rtau(i)-by-rtau(i) */
/*     upper triangular matrix. */

/*     For JOBCON = 'F', the pencil ( Bc Ac-lambda*Ec ) has full */
/*     row rank NCONT for all finite lambda and is in a staircase form */
/*     with */
/*                     _     _          _        _ */
/*                   ( A1,0  A1,1  ...  A1,k-1   A1,k  ) */
/*                   (       _          _        _     ) */
/*       ( Bc Ac ) = (  0    A2,1  ...  A2,k-1   A2,k  ) ,   (3) */
/*                   (             ...  _        _     ) */
/*                   (  0      0   ...  Ak,k-1   Ak,k  ) */

/*                     _          _        _ */
/*                   ( E1,1  ...  E1,k-1   E1,k  ) */
/*                   (            _        _     ) */
/*         Ec      = (   0   ...  E2,k-1   E2,k  ) ,         (4) */
/*                   (       ...           _     ) */
/*                   (   0   ...    0      Ek,k  ) */
/*           _ */
/*     where Ai,i-1 is an rtau(i)-by-rtau(i-1) full row rank matrix */
/*                            _ */
/*     (with rtau(0) = M) and Ei,i is an rtau(i)-by-rtau(i) */
/*     upper triangular matrix. */

/*     For JOBCON = 'C', the (N-NCONT)-by-(N-NCONT) regular pencil */
/*     Anc - lambda*Enc has the form */

/*                         ( Ainc - lambda*Einc         *          ) */
/*      Anc - lambda*Enc = (                                       ) , */
/*                         (        0           Afnc - lambda*Efnc ) */

/*     where: */
/*       1) the NIUCON-by-NIUCON regular pencil Ainc - lambda*Einc, */
/*          with Ainc upper triangular and nonsingular, contains the */
/*          uncontrollable infinite eigenvalues of A - lambda*E; */
/*       2) the (N-NCONT-NIUCON)-by-(N-NCONT-NIUCON) regular pencil */
/*          Afnc - lambda*Efnc, with Efnc upper triangular and */
/*          nonsingular, contains the uncontrollable finite */
/*          eigenvalues of A - lambda*E. */

/*     Note: The significance of the two diagonal blocks can be */
/*           interchanged by calling the routine with the */
/*           arguments A and E interchanged. In this case, */
/*           Ainc - lambda*Einc contains the uncontrollable zero */
/*           eigenvalues of A - lambda*E, while Afnc - lambda*Efnc */
/*           contains the uncontrollable nonzero finite and infinite */
/*           eigenvalues of A - lambda*E. */

/*     For JOBCON = 'F', the pencil Anc - lambda*Enc has the form */

/*        Anc - lambda*Enc = Afnc - lambda*Efnc , */

/*     where the regular pencil Afnc - lambda*Efnc, with Efnc */
/*     upper triangular and nonsingular, contains the uncontrollable */
/*     finite eigenvalues of A - lambda*E. */

/*     For JOBCON = 'I', the pencil Anc - lambda*Enc has the form */

/*        Anc - lambda*Enc = Ainc - lambda*Einc , */

/*     where the regular pencil Ainc - lambda*Einc, with Ainc */
/*     upper triangular and nonsingular, contains the uncontrollable */
/*     nonzero finite and infinite eigenvalues of A - lambda*E. */

/*     The left and/or right orthogonal transformations Q and Z */
/*     performed to reduce the system matrices can be optionally */
/*     accumulated. */

/*     The reduced order descriptor system (Ac-lambda*Ec,Bc,Cc) has */
/*     the same transfer-function matrix as the original system */
/*     (A-lambda*E,B,C). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBCON  CHARACTER*1 */
/*             = 'C':  separate both finite and infinite uncontrollable */
/*                     eigenvalues; */
/*             = 'F':  separate only finite uncontrollable eigenvalues: */
/*             = 'I':  separate only nonzero finite and infinite */
/*                     uncontrollable eigenvalues. */

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

/*     N       (input) INTEGER */
/*             The dimension of the descriptor state vector; also the */
/*             order of square matrices A and E, the number of rows of */
/*             matrix B, and the number of columns of matrix C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The dimension of descriptor system input vector; also the */
/*             number of columns of matrix B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The dimension of descriptor system output vector; also the */
/*             number of rows of matrix C.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the N-by-N state matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the transformed state matrix Q'*A*Z, */

/*                                ( Ac   *  ) */
/*                       Q'*A*Z = (         ) , */
/*                                ( 0   Anc ) */

/*             where Ac is NCONT-by-NCONT and Anc is */
/*             (N-NCONT)-by-(N-NCONT). */
/*             If JOBCON = 'F', the matrix ( Bc Ac ) is in the */
/*             controllability staircase form (3). */
/*             If JOBCON = 'C' or 'I', the submatrix Ac is upper */
/*             triangular. */
/*             If JOBCON = 'C', the Anc matrix has the form */

/*                             ( Ainc   *  ) */
/*                       Anc = (           ) , */
/*                             (  0   Afnc ) */

/*             where the NIUCON-by-NIUCON matrix Ainc is nonsingular and */
/*             upper triangular. */
/*             If JOBCON = 'I', Anc is nonsingular and upper triangular. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the N-by-N descriptor matrix E. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the transformed descriptor matrix Q'*E*Z, */

/*                                ( Ec   *  ) */
/*                       Q'*E*Z = (         ) , */
/*                                ( 0   Enc ) */

/*             where Ec is NCONT-by-NCONT and Enc is */
/*             (N-NCONT)-by-(N-NCONT). */
/*             If JOBCON = 'C' or 'I', the matrix ( Bc Ec ) is in the */
/*             controllability staircase form (1). */
/*             If JOBCON = 'F', the submatrix Ec is upper triangular. */
/*             If JOBCON = 'C', the Enc matrix has the form */

/*                             ( Einc   *  ) */
/*                       Enc = (           ) , */
/*                             (  0   Efnc ) */

/*             where the NIUCON-by-NIUCON matrix Einc is nilpotent */
/*             and the (N-NCONT-NIUCON)-by-(N-NCONT-NIUCON) matrix Efnc */
/*             is nonsingular and upper triangular. */
/*             If JOBCON = 'F', Enc is nonsingular and upper triangular. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the N-by-M input matrix B. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the transformed input matrix */

/*                              ( Bc ) */
/*                       Q'*B = (    ) , */
/*                              ( 0  ) */

/*              where Bc is NCONT-by-M. */
/*              For JOBCON = 'C' or 'I', the matrix ( Bc Ec ) is in the */
/*              controllability staircase form (1). */
/*              For JOBCON = 'F', the matrix ( Bc Ac ) is in the */
/*              controllability staircase form (3). */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the transformed matrix C*Z. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             If COMPQ = 'N': Q is not referenced. */
/*             If COMPQ = 'I': on entry, Q need not be set; */
/*                             on exit, the leading N-by-N part of this */
/*                             array contains the orthogonal matrix Q, */
/*                             where Q' is the product of transformations */
/*                             which are applied to A, E, and B on */
/*                             the left. */
/*             If COMPQ = 'U': on entry, the leading N-by-N part of this */
/*                             array must contain an orthogonal matrix */
/*                             Qc; */
/*                             on exit, the leading N-by-N part of this */
/*                             array contains the orthogonal matrix */
/*                             Qc*Q. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. */
/*             LDQ >= 1,        if COMPQ = 'N'; */
/*             LDQ >= MAX(1,N), if COMPQ = 'U' or 'I'. */

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

/*     NCONT   (output) INTEGER */
/*             The order of the reduced matrices Ac and Ec, and the */
/*             number of rows of reduced matrix Bc; also the order of */
/*             the controllable part of the pair (A-lambda*E,B). */

/*     NIUCON  (output) INTEGER */
/*             For JOBCON = 'C', the order of the reduced matrices */
/*             Ainc and Einc; also the number of uncontrollable */
/*             infinite eigenvalues of the pencil A - lambda*E. */
/*             For JOBCON = 'F' or 'I', NIUCON has no significance */
/*             and is set to zero. */

/*     NRBLCK  (output) INTEGER */
/*             For JOBCON = 'C' or 'I', the number k, of full row rank */
/*                    _ */
/*             blocks Ei,i in the staircase form of the pencil */
/*             (Bc Ec-lambda*Ac) (see (1) and (2)). */
/*             For JOBCON = 'F', the number k, of full row rank blocks */
/*             _ */
/*             Ai,i in the staircase form of the pencil (Bc Ac-lambda*Ec) */
/*             (see (3) and (4)). */

/*     RTAU    (output) INTEGER array, dimension (N) */
/*             RTAU(i), for i = 1, ..., NRBLCK, is the row dimension of */
/*                                     _         _ */
/*             the full row rank block Ei,i-1 or Ai,i-1 in the staircase */
/*             form (1) or (3) for JOBCON = 'C' or 'I', or */
/*             for JOBCON = 'F', respectively. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in rank determinations when */
/*             transforming (A-lambda*E, B). If the user sets TOL > 0, */
/*             then the given value of TOL is used as a lower bound for */
/*             reciprocal condition numbers in rank determinations; a */
/*             (sub)matrix whose estimated condition number is less than */
/*             1/TOL is considered to be of full rank.  If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance, */
/*             defined by  TOLDEF = N*N*EPS,  is used instead, where EPS */
/*             is the machine precision (see LAPACK Library routine */
/*             DLAMCH).  TOL < 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M) */

/*     DWORK   DOUBLE PRECISION array, dimension MAX(N,2*M) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The subroutine is based on the reduction algorithms of [1]. */

/*     REFERENCES */

/*     [1] A. Varga */
/*         Computation of Irreducible Generalized State-Space */
/*         Realizations. */
/*         Kybernetika, vol. 26, pp. 89-106, 1990. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically backward stable and requires */
/*     0( N**3 )  floating point operations. */

/*     FURTHER COMMENTS */

/*     If the system matrices A, E and B are badly scaled, it is */
/*     generally recommendable to scale them with the SLICOT routine */
/*     TG01AD, before calling TG01HD. */

/*     CONTRIBUTOR */

/*     C. Oara, University "Politehnica" Bucharest. */
/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     March 1999. Based on the RASP routine RPDSCF. */

/*     REVISIONS */

/*     July 1999, V. Sima, Research Institute for Informatics, Bucharest. */

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

/*     Decode JOBCON. */

#line 422 "TG01HD.f"
    /* Parameter adjustments */
#line 422 "TG01HD.f"
    a_dim1 = *lda;
#line 422 "TG01HD.f"
    a_offset = 1 + a_dim1;
#line 422 "TG01HD.f"
    a -= a_offset;
#line 422 "TG01HD.f"
    e_dim1 = *lde;
#line 422 "TG01HD.f"
    e_offset = 1 + e_dim1;
#line 422 "TG01HD.f"
    e -= e_offset;
#line 422 "TG01HD.f"
    b_dim1 = *ldb;
#line 422 "TG01HD.f"
    b_offset = 1 + b_dim1;
#line 422 "TG01HD.f"
    b -= b_offset;
#line 422 "TG01HD.f"
    c_dim1 = *ldc;
#line 422 "TG01HD.f"
    c_offset = 1 + c_dim1;
#line 422 "TG01HD.f"
    c__ -= c_offset;
#line 422 "TG01HD.f"
    q_dim1 = *ldq;
#line 422 "TG01HD.f"
    q_offset = 1 + q_dim1;
#line 422 "TG01HD.f"
    q -= q_offset;
#line 422 "TG01HD.f"
    z_dim1 = *ldz;
#line 422 "TG01HD.f"
    z_offset = 1 + z_dim1;
#line 422 "TG01HD.f"
    z__ -= z_offset;
#line 422 "TG01HD.f"
    --rtau;
#line 422 "TG01HD.f"
    --iwork;
#line 422 "TG01HD.f"
    --dwork;
#line 422 "TG01HD.f"

#line 422 "TG01HD.f"
    /* Function Body */
#line 422 "TG01HD.f"
    if (lsame_(jobcon, "C", (ftnlen)1, (ftnlen)1)) {
#line 423 "TG01HD.f"
	fincon = TRUE_;
#line 424 "TG01HD.f"
	infcon = TRUE_;
#line 425 "TG01HD.f"
    } else if (lsame_(jobcon, "F", (ftnlen)1, (ftnlen)1)) {
#line 426 "TG01HD.f"
	fincon = TRUE_;
#line 427 "TG01HD.f"
	infcon = FALSE_;
#line 428 "TG01HD.f"
    } else if (lsame_(jobcon, "I", (ftnlen)1, (ftnlen)1)) {
#line 429 "TG01HD.f"
	fincon = FALSE_;
#line 430 "TG01HD.f"
	infcon = TRUE_;
#line 431 "TG01HD.f"
    } else {
#line 432 "TG01HD.f"
	fincon = FALSE_;
#line 433 "TG01HD.f"
	infcon = FALSE_;
#line 434 "TG01HD.f"
    }

/*     Decode COMPQ. */

#line 438 "TG01HD.f"
    if (lsame_(compq, "N", (ftnlen)1, (ftnlen)1)) {
#line 439 "TG01HD.f"
	ilq = FALSE_;
#line 440 "TG01HD.f"
	icompq = 1;
#line 441 "TG01HD.f"
    } else if (lsame_(compq, "U", (ftnlen)1, (ftnlen)1)) {
#line 442 "TG01HD.f"
	ilq = TRUE_;
#line 443 "TG01HD.f"
	icompq = 2;
#line 444 "TG01HD.f"
    } else if (lsame_(compq, "I", (ftnlen)1, (ftnlen)1)) {
#line 445 "TG01HD.f"
	ilq = TRUE_;
#line 446 "TG01HD.f"
	icompq = 3;
#line 447 "TG01HD.f"
    } else {
#line 448 "TG01HD.f"
	icompq = 0;
#line 449 "TG01HD.f"
    }

/*     Decode COMPZ. */

#line 453 "TG01HD.f"
    if (lsame_(compz, "N", (ftnlen)1, (ftnlen)1)) {
#line 454 "TG01HD.f"
	ilz = FALSE_;
#line 455 "TG01HD.f"
	icompz = 1;
#line 456 "TG01HD.f"
    } else if (lsame_(compz, "U", (ftnlen)1, (ftnlen)1)) {
#line 457 "TG01HD.f"
	ilz = TRUE_;
#line 458 "TG01HD.f"
	icompz = 2;
#line 459 "TG01HD.f"
    } else if (lsame_(compz, "I", (ftnlen)1, (ftnlen)1)) {
#line 460 "TG01HD.f"
	ilz = TRUE_;
#line 461 "TG01HD.f"
	icompz = 3;
#line 462 "TG01HD.f"
    } else {
#line 463 "TG01HD.f"
	icompz = 0;
#line 464 "TG01HD.f"
    }

/*     Test the input scalar parameters. */

#line 468 "TG01HD.f"
    *info = 0;
#line 469 "TG01HD.f"
    if (! fincon && ! infcon) {
#line 470 "TG01HD.f"
	*info = -1;
#line 471 "TG01HD.f"
    } else if (icompq <= 0) {
#line 472 "TG01HD.f"
	*info = -2;
#line 473 "TG01HD.f"
    } else if (icompz <= 0) {
#line 474 "TG01HD.f"
	*info = -3;
#line 475 "TG01HD.f"
    } else if (*n < 0) {
#line 476 "TG01HD.f"
	*info = -4;
#line 477 "TG01HD.f"
    } else if (*m < 0) {
#line 478 "TG01HD.f"
	*info = -5;
#line 479 "TG01HD.f"
    } else if (*p < 0) {
#line 480 "TG01HD.f"
	*info = -6;
#line 481 "TG01HD.f"
    } else if (*lda < max(1,*n)) {
#line 482 "TG01HD.f"
	*info = -8;
#line 483 "TG01HD.f"
    } else if (*lde < max(1,*n)) {
#line 484 "TG01HD.f"
	*info = -10;
#line 485 "TG01HD.f"
    } else if (*ldb < max(1,*n)) {
#line 486 "TG01HD.f"
	*info = -12;
#line 487 "TG01HD.f"
    } else if (*ldc < max(1,*p)) {
#line 488 "TG01HD.f"
	*info = -14;
#line 489 "TG01HD.f"
    } else if (ilq && *ldq < *n || *ldq < 1) {
#line 490 "TG01HD.f"
	*info = -16;
#line 491 "TG01HD.f"
    } else if (ilz && *ldz < *n || *ldz < 1) {
#line 492 "TG01HD.f"
	*info = -18;
#line 493 "TG01HD.f"
    } else if (*tol >= 1.) {
#line 494 "TG01HD.f"
	*info = -23;
#line 495 "TG01HD.f"
    }
#line 496 "TG01HD.f"
    if (*info != 0) {
#line 497 "TG01HD.f"
	i__1 = -(*info);
#line 497 "TG01HD.f"
	xerbla_("TG01HD", &i__1, (ftnlen)6);
#line 498 "TG01HD.f"
	return 0;
#line 499 "TG01HD.f"
    }

#line 501 "TG01HD.f"
    *(unsigned char *)jobq = *(unsigned char *)compq;
#line 502 "TG01HD.f"
    *(unsigned char *)jobz = *(unsigned char *)compz;

#line 504 "TG01HD.f"
    if (fincon) {

/*        Perform finite controllability form reduction. */

/* Computing MAX */
#line 508 "TG01HD.f"
	i__2 = 0, i__3 = *n - 1;
#line 508 "TG01HD.f"
	i__1 = max(i__2,i__3);
#line 508 "TG01HD.f"
	tg01hx_(jobq, jobz, n, n, m, p, n, &i__1, &a[a_offset], lda, &e[
		e_offset], lde, &b[b_offset], ldb, &c__[c_offset], ldc, &q[
		q_offset], ldq, &z__[z_offset], ldz, &nr, nrblck, &rtau[1], 
		tol, &iwork[1], &dwork[1], info, (ftnlen)1, (ftnlen)1);
#line 511 "TG01HD.f"
	if (*nrblck > 1) {
#line 512 "TG01HD.f"
	    lba = rtau[1] + rtau[2] - 1;
#line 513 "TG01HD.f"
	} else if (*nrblck == 1) {
#line 514 "TG01HD.f"
	    lba = rtau[1] - 1;
#line 515 "TG01HD.f"
	} else {
#line 516 "TG01HD.f"
	    lba = 0;
#line 517 "TG01HD.f"
	}
#line 518 "TG01HD.f"
	if (ilq) {
#line 518 "TG01HD.f"
	    *(unsigned char *)jobq = 'U';
#line 518 "TG01HD.f"
	}
#line 519 "TG01HD.f"
	if (ilz) {
#line 519 "TG01HD.f"
	    *(unsigned char *)jobz = 'U';
#line 519 "TG01HD.f"
	}
#line 520 "TG01HD.f"
    } else {
#line 521 "TG01HD.f"
	nr = *n;
/* Computing MAX */
#line 522 "TG01HD.f"
	i__1 = 0, i__2 = *n - 1;
#line 522 "TG01HD.f"
	lba = max(i__1,i__2);
#line 523 "TG01HD.f"
    }

#line 525 "TG01HD.f"
    if (infcon) {

/*        Perform infinite controllability form reduction. */

#line 529 "TG01HD.f"
	tg01hx_(jobq, jobz, n, n, m, p, &nr, &lba, &e[e_offset], lde, &a[
		a_offset], lda, &b[b_offset], ldb, &c__[c_offset], ldc, &q[
		q_offset], ldq, &z__[z_offset], ldz, ncont, nrblck, &rtau[1], 
		tol, &iwork[1], &dwork[1], info, (ftnlen)1, (ftnlen)1);
#line 532 "TG01HD.f"
	if (fincon) {
#line 533 "TG01HD.f"
	    *niucon = nr - *ncont;
#line 534 "TG01HD.f"
	} else {
#line 535 "TG01HD.f"
	    *niucon = 0;
#line 536 "TG01HD.f"
	}
#line 537 "TG01HD.f"
    } else {
#line 538 "TG01HD.f"
	*ncont = nr;
#line 539 "TG01HD.f"
	*niucon = 0;
#line 540 "TG01HD.f"
    }

#line 542 "TG01HD.f"
    return 0;

/* *** Last line of TG01HD *** */
} /* tg01hd_ */

