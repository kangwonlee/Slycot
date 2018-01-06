#line 1 "MB03ZA.f"
/* MB03ZA.f -- translated by f2c (version 20100827).
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

#line 1 "MB03ZA.f"
/* Table of constant values */

static doublereal c_b20 = 0.;
static doublereal c_b21 = 1.;
static integer c__4 = 4;
static logical c_true = TRUE_;
static integer c__1 = 1;
static integer c__2 = 2;
static doublereal c_b479 = -1.;
static integer c__12 = 12;

/* Subroutine */ int mb03za_(char *compc, char *compu, char *compv, char *
	compw, char *which, logical *select, integer *n, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, doublereal *c__, integer *
	ldc, doublereal *u1, integer *ldu1, doublereal *u2, integer *ldu2, 
	doublereal *v1, integer *ldv1, doublereal *v2, integer *ldv2, 
	doublereal *w, integer *ldw, doublereal *wr, doublereal *wi, integer *
	m, doublereal *dwork, integer *ldwork, integer *info, ftnlen 
	compc_len, ftnlen compu_len, ftnlen compv_len, ftnlen compw_len, 
	ftnlen which_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, u1_dim1, 
	    u1_offset, u2_dim1, u2_offset, v1_dim1, v1_offset, v2_dim1, 
	    v2_offset, w_dim1, w_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, k, l;
    static doublereal q[16]	/* was [4][4] */, t[16]	/* was [4][4] */, z__[
	    16]	/* was [4][4] */;
    static integer nb, mm, ks, pw, nbf, nbl;
    static doublereal dw12[12];
    static integer len, pwc, pwd, pos, here;
    static logical pair;
    static integer idum[1], ierr;
    static logical ldum[1];
    static integer pwck, ifst, pwdl;
    static doublereal temp;
    static logical swap;
    static integer ilst;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgees_(char *, char *, L_fp, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, logical *, integer *, ftnlen, 
	    ftnlen), dgemm_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen), mb03wa_(
	    logical *, logical *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen), lfdum_();
    static logical wantc;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal winew[4];
    static logical initw, wantu, wantv, wantw;
    static doublereal wrnew[4];
    static logical cmpall;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static logical selnew[4];
    static integer nbnext;
    extern /* Subroutine */ int dtrsen_(char *, char *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, integer *, ftnlen, ftnlen);
    static integer wrkmin;


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

/*     1. To compute, for a given matrix pair (A,B) in periodic Schur */
/*        form, orthogonal matrices Ur and Vr so that */

/*            T           [ A11  A12 ]     T           [ B11  B12 ] */
/*          Vr * A * Ur = [          ],  Ur * B * Vr = [          ], (1) */
/*                        [  0   A22 ]                 [  0   B22 ] */

/*        is in periodic Schur form, and the eigenvalues of A11*B11 */
/*        form a selected cluster of eigenvalues. */

/*     2. To compute an orthogonal matrix W so that */

/*                   T  [  0  -A11 ]       [  R11   R12 ] */
/*                  W * [          ] * W = [            ],           (2) */
/*                      [ B11   0  ]       [   0    R22 ] */

/*        where the eigenvalues of R11 and -R22 coincide and have */
/*        positive real part. */

/*     Optionally, the matrix C is overwritten by Ur'*C*Vr. */

/*     All eigenvalues of A11*B11 must either be complex or real and */
/*     negative. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     COMPC   CHARACTER*1 */
/*             = 'U':  update the matrix C; */
/*             = 'N':  do not update C. */

/*     COMPU   CHARACTER*1 */
/*             = 'U':  update the matrices U1 and U2; */
/*             = 'N':  do not update U1 and U2. */
/*             See the description of U1 and U2. */

/*     COMPV   CHARACTER*1 */
/*             = 'U':  update the matrices V1 and V2; */
/*             = 'N':  do not update V1 and V2. */
/*             See the description of V1 and V2. */

/*     COMPW   CHARACTER*1 */
/*             Indicates whether or not the user wishes to accumulate */
/*             the matrix W as follows: */
/*             = 'N':  the matrix W is not required; */
/*             = 'I':  W is initialized to the unit matrix and the */
/*                     orthogonal transformation matrix W is returned; */
/*             = 'V':  W must contain an orthogonal matrix Q on entry, */
/*                     and the product Q*W is returned. */

/*     WHICH   CHARACTER*1 */
/*             = 'A':  select all eigenvalues, this effectively means */
/*                     that Ur and Vr are identity matrices and A11 = A, */
/*                     B11 = B; */
/*             = 'S':  select a cluster of eigenvalues specified by */
/*                     SELECT. */

/*     SELECT  LOGICAL array, dimension (N) */
/*             If WHICH = 'S', then SELECT specifies the eigenvalues of */
/*             A*B in the selected cluster. To select a real eigenvalue */
/*             w(j), SELECT(j) must be set to .TRUE.. To select a complex */
/*             conjugate pair of eigenvalues w(j) and w(j+1), */
/*             corresponding to a 2-by-2 diagonal block in A, both */
/*             SELECT(j) and SELECT(j+1) must be set to .TRUE.; a complex */
/*             conjugate pair of eigenvalues must be either both included */
/*             in the cluster or both excluded. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A. N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper quasi-triangular matrix A of the matrix */
/*             pair (A,B) in periodic Schur form. */
/*             On exit, the leading M-by-M part of this array contains */
/*             the matrix R22 in (2). */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper triangular matrix B of the matrix pair */
/*             (A,B) in periodic Schur form. */
/*             On exit, the leading N-by-N part of this array is */
/*             overwritten. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, if COMPC = 'U', the leading N-by-N part of this */
/*             array must contain a general matrix C. */
/*             On exit, if COMPC = 'U', the leading N-by-N part of this */
/*             array contains the updated matrix Ur'*C*Vr. */
/*             If COMPC = 'N' or WHICH = 'A', this array is not */
/*             referenced. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= 1. */
/*             LDC >= N,  if COMPC = 'U' and WHICH = 'S'. */

/*     U1      (input/output) DOUBLE PRECISION array, dimension (LDU1,N) */
/*             On entry, if COMPU = 'U' and WHICH = 'S', the leading */
/*             N-by-N part of this array must contain U1, the (1,1) */
/*             block of an orthogonal symplectic matrix */
/*             U = [ U1, U2; -U2, U1 ]. */
/*             On exit, if COMPU = 'U' and WHICH = 'S', the leading */
/*             N-by-N part of this array contains U1*Ur. */
/*             If COMPU = 'N' or WHICH = 'A', this array is not */
/*             referenced. */

/*     LDU1    INTEGER */
/*             The leading dimension of the array U1.  LDU1 >= 1. */
/*             LDU1 >= N,  if COMPU = 'U' and WHICH = 'S'. */

/*     U2      (input/output) DOUBLE PRECISION array, dimension (LDU2,N) */
/*             On entry, if COMPU = 'U' and WHICH = 'S', the leading */
/*             N-by-N part of this array must contain U2, the (1,2) */
/*             block of an orthogonal symplectic matrix */
/*             U = [ U1, U2; -U2, U1 ]. */
/*             On exit, if COMPU = 'U' and WHICH = 'S', the leading */
/*             N-by-N part of this array contains U2*Ur. */
/*             If COMPU = 'N' or WHICH = 'A', this array is not */
/*             referenced. */

/*     LDU2    INTEGER */
/*             The leading dimension of the array U2.  LDU2 >= 1. */
/*             LDU2 >= N,  if COMPU = 'U' and WHICH = 'S'. */

/*     V1      (input/output) DOUBLE PRECISION array, dimension (LDV1,N) */
/*             On entry, if COMPV = 'U' and WHICH = 'S', the leading */
/*             N-by-N part of this array must contain V1, the (1,1) */
/*             block of an orthogonal symplectic matrix */
/*             V = [ V1, V2; -V2, V1 ]. */
/*             On exit, if COMPV = 'U' and WHICH = 'S', the leading */
/*             N-by-N part of this array contains V1*Vr. */
/*             If COMPV = 'N' or WHICH = 'A', this array is not */
/*             referenced. */

/*     LDV1    INTEGER */
/*             The leading dimension of the array V1.  LDV1 >= 1. */
/*             LDV1 >= N,  if COMPV = 'U' and WHICH = 'S'. */

/*     V2      (input/output) DOUBLE PRECISION array, dimension (LDV2,N) */
/*             On entry, if COMPV = 'U' and WHICH = 'S', the leading */
/*             N-by-N part of this array must contain V2, the (1,2) */
/*             block of an orthogonal symplectic matrix */
/*             V = [ V1, V2; -V2, V1 ]. */
/*             On exit, if COMPV = 'U' and WHICH = 'S', the leading */
/*             N-by-N part of this array contains V2*Vr. */
/*             If COMPV = 'N' or WHICH = 'A', this array is not */
/*             referenced. */

/*     LDV2    INTEGER */
/*             The leading dimension of the array V2.  LDV2 >= 1. */
/*             LDV2 >= N,  if COMPV = 'U' and WHICH = 'S'. */

/*     W       (input/output) DOUBLE PRECISION array, dimension (LDW,2*M) */
/*             On entry, if COMPW = 'V', then the leading 2*M-by-2*M part */
/*             of this array must contain a matrix W. */
/*             If COMPW = 'I', then W need not be set on entry, W is set */
/*             to the identity matrix. */
/*             On exit, if COMPW = 'I' or 'V' the leading 2*M-by-2*M part */
/*             of this array is post-multiplied by the transformation */
/*             matrix that produced (2). */
/*             If COMPW = 'N', this array is not referenced. */

/*     LDW     INTEGER */
/*             The leading dimension of the array W.  LDW >= 1. */
/*             LDW >= 2*M,  if COMPW = 'I' or COMPW = 'V'. */

/*     WR      (output) DOUBLE PRECISION array, dimension (M) */
/*     WI      (output) DOUBLE PRECISION array, dimension (M) */
/*             The real and imaginary parts, respectively, of the */
/*             eigenvalues of R22. The eigenvalues are stored in the same */
/*             order as on the diagonal of R22, with */
/*             WR(i) = R22(i,i) and, if R22(i:i+1,i:i+1) is a 2-by-2 */
/*             diagonal block, WI(i) > 0 and WI(i+1) = -WI(i). */
/*             In exact arithmetic, these eigenvalue are the positive */
/*             square roots of the selected eigenvalues of the product */
/*             A*B. However, if an eigenvalue is sufficiently */
/*             ill-conditioned, then its value may differ significantly. */

/*     M       (output) INTEGER */
/*             The number of selected eigenvalues. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = -28,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 1, 4*N, 8*M ). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  reordering of the product A*B in Step 1 failed */
/*                   because some eigenvalues are too close to separate; */
/*             = 2:  reordering of some submatrix in Step 2 failed */
/*                   because some eigenvalues are too close to separate; */
/*             = 3:  the QR algorithm failed to compute the Schur form */
/*                   of some submatrix in Step 2; */
/*             = 4:  the condition that all eigenvalues of A11*B11 must */
/*                   either be complex or real and negative is */
/*                   numerically violated. */

/*     METHOD */

/*     Step 1 is performed using a reordering technique analogous to the */
/*     LAPACK routine DTGSEN for reordering matrix pencils [1,2]. Step 2 */
/*     is an implementation of Algorithm 2 in [3]. It requires O(M*N*N) */
/*     floating point operations. */

/*     REFERENCES */

/*     [1] Kagstrom, B. */
/*         A direct method for reordering eigenvalues in the generalized */
/*         real Schur form of a regular matrix pair (A,B), in M.S. Moonen */
/*         et al (eds), Linear Algebra for Large Scale and Real-Time */
/*         Applications, Kluwer Academic Publ., 1993, pp. 195-218. */

/*     [2] Kagstrom, B. and Poromaa P.: */
/*         Computing eigenspaces with specified eigenvalues of a regular */
/*         matrix pair (A, B) and condition estimation: Theory, */
/*         algorithms and software, Numer. Algorithms, 1996, vol. 12, */
/*         pp. 369-407. */

/*     [3] Benner, P., Mehrmann, V., and Xu, H. */
/*         A new method for computing the stable invariant subspace of a */
/*         real Hamiltonian matrix,  J. Comput. Appl. Math., 86, */
/*         pp. 17-43, 1997. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLABMX). */

/*     KEYWORDS */

/*     Hamiltonian matrix, invariant subspace. */

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

/*     Decode and check input parameters */

#line 322 "MB03ZA.f"
    /* Parameter adjustments */
#line 322 "MB03ZA.f"
    --select;
#line 322 "MB03ZA.f"
    a_dim1 = *lda;
#line 322 "MB03ZA.f"
    a_offset = 1 + a_dim1;
#line 322 "MB03ZA.f"
    a -= a_offset;
#line 322 "MB03ZA.f"
    b_dim1 = *ldb;
#line 322 "MB03ZA.f"
    b_offset = 1 + b_dim1;
#line 322 "MB03ZA.f"
    b -= b_offset;
#line 322 "MB03ZA.f"
    c_dim1 = *ldc;
#line 322 "MB03ZA.f"
    c_offset = 1 + c_dim1;
#line 322 "MB03ZA.f"
    c__ -= c_offset;
#line 322 "MB03ZA.f"
    u1_dim1 = *ldu1;
#line 322 "MB03ZA.f"
    u1_offset = 1 + u1_dim1;
#line 322 "MB03ZA.f"
    u1 -= u1_offset;
#line 322 "MB03ZA.f"
    u2_dim1 = *ldu2;
#line 322 "MB03ZA.f"
    u2_offset = 1 + u2_dim1;
#line 322 "MB03ZA.f"
    u2 -= u2_offset;
#line 322 "MB03ZA.f"
    v1_dim1 = *ldv1;
#line 322 "MB03ZA.f"
    v1_offset = 1 + v1_dim1;
#line 322 "MB03ZA.f"
    v1 -= v1_offset;
#line 322 "MB03ZA.f"
    v2_dim1 = *ldv2;
#line 322 "MB03ZA.f"
    v2_offset = 1 + v2_dim1;
#line 322 "MB03ZA.f"
    v2 -= v2_offset;
#line 322 "MB03ZA.f"
    w_dim1 = *ldw;
#line 322 "MB03ZA.f"
    w_offset = 1 + w_dim1;
#line 322 "MB03ZA.f"
    w -= w_offset;
#line 322 "MB03ZA.f"
    --wr;
#line 322 "MB03ZA.f"
    --wi;
#line 322 "MB03ZA.f"
    --dwork;
#line 322 "MB03ZA.f"

#line 322 "MB03ZA.f"
    /* Function Body */
#line 322 "MB03ZA.f"
    wantc = lsame_(compc, "U", (ftnlen)1, (ftnlen)1);
#line 323 "MB03ZA.f"
    wantu = lsame_(compu, "U", (ftnlen)1, (ftnlen)1);
#line 324 "MB03ZA.f"
    wantv = lsame_(compv, "U", (ftnlen)1, (ftnlen)1);
#line 325 "MB03ZA.f"
    initw = lsame_(compw, "I", (ftnlen)1, (ftnlen)1);
#line 326 "MB03ZA.f"
    wantw = initw || lsame_(compw, "V", (ftnlen)1, (ftnlen)1);
#line 327 "MB03ZA.f"
    cmpall = lsame_(which, "A", (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 328 "MB03ZA.f"
    i__1 = 1, i__2 = *n << 2;
#line 328 "MB03ZA.f"
    wrkmin = max(i__1,i__2);

#line 330 "MB03ZA.f"
    *info = 0;
#line 331 "MB03ZA.f"
    if (! wantc && ! lsame_(compc, "N", (ftnlen)1, (ftnlen)1)) {
#line 332 "MB03ZA.f"
	*info = -1;
#line 333 "MB03ZA.f"
    } else if (! wantu && ! lsame_(compu, "N", (ftnlen)1, (ftnlen)1)) {
#line 334 "MB03ZA.f"
	*info = -2;
#line 335 "MB03ZA.f"
    } else if (! wantv && ! lsame_(compv, "N", (ftnlen)1, (ftnlen)1)) {
#line 336 "MB03ZA.f"
	*info = -3;
#line 337 "MB03ZA.f"
    } else if (! wantw && ! lsame_(compw, "N", (ftnlen)1, (ftnlen)1)) {
#line 338 "MB03ZA.f"
	*info = -4;
#line 339 "MB03ZA.f"
    } else if (! cmpall && ! lsame_(which, "S", (ftnlen)1, (ftnlen)1)) {
#line 340 "MB03ZA.f"
	*info = -5;
#line 341 "MB03ZA.f"
    } else {
#line 342 "MB03ZA.f"
	if (cmpall) {
#line 343 "MB03ZA.f"
	    *m = *n;
#line 344 "MB03ZA.f"
	} else {

/*           Set M to the dimension of the specified invariant subspace. */

#line 348 "MB03ZA.f"
	    *m = 0;
#line 349 "MB03ZA.f"
	    pair = FALSE_;
#line 350 "MB03ZA.f"
	    i__1 = *n;
#line 350 "MB03ZA.f"
	    for (k = 1; k <= i__1; ++k) {
#line 351 "MB03ZA.f"
		if (pair) {
#line 352 "MB03ZA.f"
		    pair = FALSE_;
#line 353 "MB03ZA.f"
		} else {
#line 354 "MB03ZA.f"
		    if (k < *n) {
#line 355 "MB03ZA.f"
			if (a[k + 1 + k * a_dim1] == 0.) {
#line 356 "MB03ZA.f"
			    if (select[k]) {
#line 356 "MB03ZA.f"
				++(*m);
#line 356 "MB03ZA.f"
			    }
#line 358 "MB03ZA.f"
			} else {
#line 359 "MB03ZA.f"
			    pair = TRUE_;
#line 360 "MB03ZA.f"
			    if (select[k] || select[k + 1]) {
#line 360 "MB03ZA.f"
				*m += 2;
#line 360 "MB03ZA.f"
			    }
#line 362 "MB03ZA.f"
			}
#line 363 "MB03ZA.f"
		    } else {
#line 364 "MB03ZA.f"
			if (select[*n]) {
#line 364 "MB03ZA.f"
			    ++(*m);
#line 364 "MB03ZA.f"
			}
#line 366 "MB03ZA.f"
		    }
#line 367 "MB03ZA.f"
		}
#line 368 "MB03ZA.f"
/* L10: */
#line 368 "MB03ZA.f"
	    }
#line 369 "MB03ZA.f"
	}

/*        Compute workspace requirements. */

/* Computing MAX */
#line 373 "MB03ZA.f"
	i__1 = wrkmin, i__2 = *m << 3;
#line 373 "MB03ZA.f"
	wrkmin = max(i__1,i__2);

#line 375 "MB03ZA.f"
	if (*n < 0) {
#line 376 "MB03ZA.f"
	    *info = -7;
#line 377 "MB03ZA.f"
	} else if (*lda < max(1,*n)) {
#line 378 "MB03ZA.f"
	    *info = -9;
#line 379 "MB03ZA.f"
	} else if (*ldb < max(1,*n)) {
#line 380 "MB03ZA.f"
	    *info = -11;
#line 381 "MB03ZA.f"
	} else if (*ldc < 1 || wantc && ! cmpall && *ldc < *n) {
#line 383 "MB03ZA.f"
	    *info = -13;
#line 384 "MB03ZA.f"
	} else if (*ldu1 < 1 || wantu && ! cmpall && *ldu1 < *n) {
#line 386 "MB03ZA.f"
	    *info = -15;
#line 387 "MB03ZA.f"
	} else if (*ldu2 < 1 || wantu && ! cmpall && *ldu2 < *n) {
#line 389 "MB03ZA.f"
	    *info = -17;
#line 390 "MB03ZA.f"
	} else if (*ldv1 < 1 || wantv && ! cmpall && *ldv1 < *n) {
#line 392 "MB03ZA.f"
	    *info = -19;
#line 393 "MB03ZA.f"
	} else if (*ldv2 < 1 || wantv && ! cmpall && *ldv2 < *n) {
#line 395 "MB03ZA.f"
	    *info = -21;
#line 396 "MB03ZA.f"
	} else if (*ldw < 1 || wantw && *ldw < *m << 1) {
#line 397 "MB03ZA.f"
	    *info = -23;
#line 398 "MB03ZA.f"
	} else if (*ldwork < wrkmin) {
#line 399 "MB03ZA.f"
	    *info = -28;
#line 400 "MB03ZA.f"
	    dwork[1] = (doublereal) wrkmin;
#line 401 "MB03ZA.f"
	}
#line 402 "MB03ZA.f"
    }

/*     Return if there were illegal values. */

#line 406 "MB03ZA.f"
    if (*info != 0) {
#line 407 "MB03ZA.f"
	i__1 = -(*info);
#line 407 "MB03ZA.f"
	xerbla_("MB03ZA", &i__1, (ftnlen)6);
#line 408 "MB03ZA.f"
	return 0;
#line 409 "MB03ZA.f"
    }

/*     Quick return if possible. */

#line 413 "MB03ZA.f"
    if (*n == 0) {
#line 414 "MB03ZA.f"
	dwork[1] = 1.;
#line 415 "MB03ZA.f"
	return 0;
#line 416 "MB03ZA.f"
    }

/*     Jump immediately to Step 2, if all eigenvalues are requested. */

#line 420 "MB03ZA.f"
    if (cmpall) {
#line 420 "MB03ZA.f"
	goto L50;
#line 420 "MB03ZA.f"
    }

/*     Step 1: Collect the selected blocks at the top-left corner of A*B. */

#line 425 "MB03ZA.f"
    ks = 0;
#line 426 "MB03ZA.f"
    pair = FALSE_;
#line 427 "MB03ZA.f"
    i__1 = *n;
#line 427 "MB03ZA.f"
    for (k = 1; k <= i__1; ++k) {
#line 428 "MB03ZA.f"
	if (pair) {
#line 429 "MB03ZA.f"
	    pair = FALSE_;
#line 430 "MB03ZA.f"
	} else {
#line 431 "MB03ZA.f"
	    swap = select[k];
#line 432 "MB03ZA.f"
	    if (k < *n) {
#line 433 "MB03ZA.f"
		if (a[k + 1 + k * a_dim1] != 0.) {
#line 434 "MB03ZA.f"
		    pair = TRUE_;
#line 435 "MB03ZA.f"
		    swap = swap || select[k + 1];
#line 436 "MB03ZA.f"
		}
#line 437 "MB03ZA.f"
	    }

#line 439 "MB03ZA.f"
	    if (pair) {
#line 440 "MB03ZA.f"
		nbf = 2;
#line 441 "MB03ZA.f"
	    } else {
#line 442 "MB03ZA.f"
		nbf = 1;
#line 443 "MB03ZA.f"
	    }

#line 445 "MB03ZA.f"
	    if (swap) {
#line 446 "MB03ZA.f"
		++ks;
#line 447 "MB03ZA.f"
		ifst = k;

/*              Swap the K-th block to position KS. */

#line 451 "MB03ZA.f"
		ilst = ks;
#line 452 "MB03ZA.f"
		nbl = 1;
#line 453 "MB03ZA.f"
		if (ilst > 1) {
#line 454 "MB03ZA.f"
		    if (a[ilst + (ilst - 1) * a_dim1] != 0.) {
#line 455 "MB03ZA.f"
			--ilst;
#line 456 "MB03ZA.f"
			nbl = 2;
#line 457 "MB03ZA.f"
		    }
#line 458 "MB03ZA.f"
		}

#line 460 "MB03ZA.f"
		if (ilst == ifst) {
#line 460 "MB03ZA.f"
		    goto L30;
#line 460 "MB03ZA.f"
		}

#line 463 "MB03ZA.f"
		here = ifst;
#line 464 "MB03ZA.f"
L20:

/*              Swap block with next one above. */

#line 468 "MB03ZA.f"
		if (nbf == 1 || nbf == 2) {

/*                 Current block either 1-by-1 or 2-by-2. */

#line 472 "MB03ZA.f"
		    nbnext = 1;
#line 473 "MB03ZA.f"
		    if (here >= 3) {
#line 474 "MB03ZA.f"
			if (a[here - 1 + (here - 2) * a_dim1] != 0.) {
#line 474 "MB03ZA.f"
			    nbnext = 2;
#line 474 "MB03ZA.f"
			}
#line 476 "MB03ZA.f"
		    }
#line 477 "MB03ZA.f"
		    pos = here - nbnext;
#line 478 "MB03ZA.f"
		    nb = nbnext + nbf;
#line 479 "MB03ZA.f"
		    dlaset_("All", &nb, &nb, &c_b20, &c_b21, q, &c__4, (
			    ftnlen)3);
#line 480 "MB03ZA.f"
		    dlaset_("All", &nb, &nb, &c_b20, &c_b21, z__, &c__4, (
			    ftnlen)3);

#line 482 "MB03ZA.f"
		    mb03wa_(&c_true, &c_true, &nbnext, &nbf, &a[pos + pos * 
			    a_dim1], lda, &b[pos + pos * b_dim1], ldb, q, &
			    c__4, z__, &c__4, &ierr);

#line 486 "MB03ZA.f"
		    if (ierr != 0) {
#line 487 "MB03ZA.f"
			dwork[1] = (doublereal) wrkmin;
#line 488 "MB03ZA.f"
			*info = 1;
#line 489 "MB03ZA.f"
			return 0;
#line 490 "MB03ZA.f"
		    }

/*                 Update rest of A. */

#line 494 "MB03ZA.f"
		    if (pos > 1) {
#line 495 "MB03ZA.f"
			i__2 = pos - 1;
#line 495 "MB03ZA.f"
			dgemm_("No Transpose", "No Transpose", &i__2, &nb, &
				nb, &c_b21, &a[pos * a_dim1 + 1], lda, z__, &
				c__4, &c_b20, &dwork[1], n, (ftnlen)12, (
				ftnlen)12);
#line 498 "MB03ZA.f"
			i__2 = pos - 1;
#line 498 "MB03ZA.f"
			dlacpy_("All", &i__2, &nb, &dwork[1], n, &a[pos * 
				a_dim1 + 1], lda, (ftnlen)3);
#line 500 "MB03ZA.f"
		    }
#line 501 "MB03ZA.f"
		    if (pos + nb <= *n) {
#line 502 "MB03ZA.f"
			i__2 = *n - pos - nb + 1;
#line 502 "MB03ZA.f"
			dgemm_("Transpose", "No Transpose", &nb, &i__2, &nb, &
				c_b21, q, &c__4, &a[pos + (pos + nb) * a_dim1]
				, lda, &c_b20, &dwork[1], &nb, (ftnlen)9, (
				ftnlen)12);
#line 505 "MB03ZA.f"
			i__2 = *n - pos - nb + 1;
#line 505 "MB03ZA.f"
			dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &a[pos + (
				pos + nb) * a_dim1], lda, (ftnlen)3);
#line 507 "MB03ZA.f"
		    }

/*                 Update rest of B. */

#line 511 "MB03ZA.f"
		    if (pos > 1) {
#line 512 "MB03ZA.f"
			i__2 = pos - 1;
#line 512 "MB03ZA.f"
			dgemm_("No Transpose", "No Transpose", &i__2, &nb, &
				nb, &c_b21, &b[pos * b_dim1 + 1], ldb, q, &
				c__4, &c_b20, &dwork[1], n, (ftnlen)12, (
				ftnlen)12);
#line 515 "MB03ZA.f"
			i__2 = pos - 1;
#line 515 "MB03ZA.f"
			dlacpy_("All", &i__2, &nb, &dwork[1], n, &b[pos * 
				b_dim1 + 1], ldb, (ftnlen)3);
#line 517 "MB03ZA.f"
		    }
#line 518 "MB03ZA.f"
		    if (pos + nb <= *n) {
#line 519 "MB03ZA.f"
			i__2 = *n - pos - nb + 1;
#line 519 "MB03ZA.f"
			dgemm_("Transpose", "No Transpose", &nb, &i__2, &nb, &
				c_b21, z__, &c__4, &b[pos + (pos + nb) * 
				b_dim1], ldb, &c_b20, &dwork[1], &nb, (ftnlen)
				9, (ftnlen)12);
#line 522 "MB03ZA.f"
			i__2 = *n - pos - nb + 1;
#line 522 "MB03ZA.f"
			dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &b[pos + (
				pos + nb) * b_dim1], ldb, (ftnlen)3);
#line 524 "MB03ZA.f"
		    }

/*                 Update C. */

#line 528 "MB03ZA.f"
		    if (wantc) {
#line 529 "MB03ZA.f"
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &c__[pos * c_dim1 + 1], ldc, q, &c__4, 
				&c_b20, &dwork[1], n, (ftnlen)12, (ftnlen)12);
#line 532 "MB03ZA.f"
			dlacpy_("All", n, &nb, &dwork[1], n, &c__[pos * 
				c_dim1 + 1], ldc, (ftnlen)3);
#line 534 "MB03ZA.f"
			dgemm_("Transpose", "No Transpose", &nb, n, &nb, &
				c_b21, z__, &c__4, &c__[pos + c_dim1], ldc, &
				c_b20, &dwork[1], &nb, (ftnlen)9, (ftnlen)12);
#line 537 "MB03ZA.f"
			dlacpy_("All", &nb, n, &dwork[1], &nb, &c__[pos + 
				c_dim1], ldc, (ftnlen)3);
#line 539 "MB03ZA.f"
		    }

/*                 Update U. */

#line 543 "MB03ZA.f"
		    if (wantu) {
#line 544 "MB03ZA.f"
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &u1[pos * u1_dim1 + 1], ldu1, z__, &
				c__4, &c_b20, &dwork[1], n, (ftnlen)12, (
				ftnlen)12);
#line 547 "MB03ZA.f"
			dlacpy_("All", n, &nb, &dwork[1], n, &u1[pos * 
				u1_dim1 + 1], ldu1, (ftnlen)3);
#line 549 "MB03ZA.f"
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &u2[pos * u2_dim1 + 1], ldu2, z__, &
				c__4, &c_b20, &dwork[1], n, (ftnlen)12, (
				ftnlen)12);
#line 552 "MB03ZA.f"
			dlacpy_("All", n, &nb, &dwork[1], n, &u2[pos * 
				u2_dim1 + 1], ldu2, (ftnlen)3);
#line 554 "MB03ZA.f"
		    }

/*                 Update V. */

#line 558 "MB03ZA.f"
		    if (wantv) {
#line 559 "MB03ZA.f"
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &v1[pos * v1_dim1 + 1], ldv1, q, &c__4,
				 &c_b20, &dwork[1], n, (ftnlen)12, (ftnlen)12)
				;
#line 562 "MB03ZA.f"
			dlacpy_("All", n, &nb, &dwork[1], n, &v1[pos * 
				v1_dim1 + 1], ldv1, (ftnlen)3);
#line 564 "MB03ZA.f"
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &v2[pos * v2_dim1 + 1], ldv2, q, &c__4,
				 &c_b20, &dwork[1], n, (ftnlen)12, (ftnlen)12)
				;
#line 567 "MB03ZA.f"
			dlacpy_("All", n, &nb, &dwork[1], n, &v2[pos * 
				v2_dim1 + 1], ldv2, (ftnlen)3);
#line 569 "MB03ZA.f"
		    }

#line 571 "MB03ZA.f"
		    here -= nbnext;

/*                 Test if 2-by-2 block breaks into two 1-by-1 blocks. */

#line 575 "MB03ZA.f"
		    if (nbf == 2) {
#line 576 "MB03ZA.f"
			if (a[here + 1 + here * a_dim1] == 0.) {
#line 576 "MB03ZA.f"
			    nbf = 3;
#line 576 "MB03ZA.f"
			}
#line 578 "MB03ZA.f"
		    }

#line 580 "MB03ZA.f"
		} else {

/*                 Current block consists of two 1 by 1 blocks each of */
/*                 which must be swapped individually. */

#line 585 "MB03ZA.f"
		    nbnext = 1;
#line 586 "MB03ZA.f"
		    if (here >= 3) {
#line 587 "MB03ZA.f"
			if (a[here - 1 + (here - 2) * a_dim1] != 0.) {
#line 587 "MB03ZA.f"
			    nbnext = 2;
#line 587 "MB03ZA.f"
			}
#line 589 "MB03ZA.f"
		    }
#line 590 "MB03ZA.f"
		    pos = here - nbnext;
#line 591 "MB03ZA.f"
		    nb = nbnext + 1;
#line 592 "MB03ZA.f"
		    dlaset_("All", &nb, &nb, &c_b20, &c_b21, q, &c__4, (
			    ftnlen)3);
#line 593 "MB03ZA.f"
		    dlaset_("All", &nb, &nb, &c_b20, &c_b21, z__, &c__4, (
			    ftnlen)3);

#line 595 "MB03ZA.f"
		    mb03wa_(&c_true, &c_true, &nbnext, &c__1, &a[pos + pos * 
			    a_dim1], lda, &b[pos + pos * b_dim1], ldb, q, &
			    c__4, z__, &c__4, &ierr);

#line 599 "MB03ZA.f"
		    if (ierr != 0) {
#line 600 "MB03ZA.f"
			dwork[1] = (doublereal) wrkmin;
#line 601 "MB03ZA.f"
			*info = 1;
#line 602 "MB03ZA.f"
			return 0;
#line 603 "MB03ZA.f"
		    }

/*                 Update rest of A. */

#line 607 "MB03ZA.f"
		    if (pos > 1) {
#line 608 "MB03ZA.f"
			i__2 = pos - 1;
#line 608 "MB03ZA.f"
			dgemm_("No Transpose", "No Transpose", &i__2, &nb, &
				nb, &c_b21, &a[pos * a_dim1 + 1], lda, z__, &
				c__4, &c_b20, &dwork[1], n, (ftnlen)12, (
				ftnlen)12);
#line 611 "MB03ZA.f"
			i__2 = pos - 1;
#line 611 "MB03ZA.f"
			dlacpy_("All", &i__2, &nb, &dwork[1], n, &a[pos * 
				a_dim1 + 1], lda, (ftnlen)3);
#line 613 "MB03ZA.f"
		    }
#line 614 "MB03ZA.f"
		    if (pos + nb <= *n) {
#line 615 "MB03ZA.f"
			i__2 = *n - pos - nb + 1;
#line 615 "MB03ZA.f"
			dgemm_("Transpose", "No Transpose", &nb, &i__2, &nb, &
				c_b21, q, &c__4, &a[pos + (pos + nb) * a_dim1]
				, lda, &c_b20, &dwork[1], &nb, (ftnlen)9, (
				ftnlen)12);
#line 618 "MB03ZA.f"
			i__2 = *n - pos - nb + 1;
#line 618 "MB03ZA.f"
			dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &a[pos + (
				pos + nb) * a_dim1], lda, (ftnlen)3);
#line 620 "MB03ZA.f"
		    }

/*                 Update rest of B. */

#line 624 "MB03ZA.f"
		    if (pos > 1) {
#line 625 "MB03ZA.f"
			i__2 = pos - 1;
#line 625 "MB03ZA.f"
			dgemm_("No Transpose", "No Transpose", &i__2, &nb, &
				nb, &c_b21, &b[pos * b_dim1 + 1], ldb, q, &
				c__4, &c_b20, &dwork[1], n, (ftnlen)12, (
				ftnlen)12);
#line 628 "MB03ZA.f"
			i__2 = pos - 1;
#line 628 "MB03ZA.f"
			dlacpy_("All", &i__2, &nb, &dwork[1], n, &b[pos * 
				b_dim1 + 1], ldb, (ftnlen)3);
#line 630 "MB03ZA.f"
		    }
#line 631 "MB03ZA.f"
		    if (pos + nb <= *n) {
#line 632 "MB03ZA.f"
			i__2 = *n - pos - nb + 1;
#line 632 "MB03ZA.f"
			dgemm_("Transpose", "No Transpose", &nb, &i__2, &nb, &
				c_b21, z__, &c__4, &b[pos + (pos + nb) * 
				b_dim1], ldb, &c_b20, &dwork[1], &nb, (ftnlen)
				9, (ftnlen)12);
#line 635 "MB03ZA.f"
			i__2 = *n - pos - nb + 1;
#line 635 "MB03ZA.f"
			dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &b[pos + (
				pos + nb) * b_dim1], ldb, (ftnlen)3);
#line 637 "MB03ZA.f"
		    }

/*                 Update C. */

#line 641 "MB03ZA.f"
		    if (wantc) {
#line 642 "MB03ZA.f"
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &c__[pos * c_dim1 + 1], ldc, q, &c__4, 
				&c_b20, &dwork[1], n, (ftnlen)12, (ftnlen)12);
#line 645 "MB03ZA.f"
			dlacpy_("All", n, &nb, &dwork[1], n, &c__[pos * 
				c_dim1 + 1], ldc, (ftnlen)3);
#line 647 "MB03ZA.f"
			dgemm_("Transpose", "No Transpose", &nb, n, &nb, &
				c_b21, z__, &c__4, &c__[pos + c_dim1], ldc, &
				c_b20, &dwork[1], &nb, (ftnlen)9, (ftnlen)12);
#line 650 "MB03ZA.f"
			dlacpy_("All", &nb, n, &dwork[1], &nb, &c__[pos + 
				c_dim1], ldc, (ftnlen)3);
#line 652 "MB03ZA.f"
		    }

/*                 Update U. */

#line 656 "MB03ZA.f"
		    if (wantu) {
#line 657 "MB03ZA.f"
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &u1[pos * u1_dim1 + 1], ldu1, z__, &
				c__4, &c_b20, &dwork[1], n, (ftnlen)12, (
				ftnlen)12);
#line 660 "MB03ZA.f"
			dlacpy_("All", n, &nb, &dwork[1], n, &u1[pos * 
				u1_dim1 + 1], ldu1, (ftnlen)3);
#line 662 "MB03ZA.f"
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &u2[pos * u2_dim1 + 1], ldu2, z__, &
				c__4, &c_b20, &dwork[1], n, (ftnlen)12, (
				ftnlen)12);
#line 665 "MB03ZA.f"
			dlacpy_("All", n, &nb, &dwork[1], n, &u2[pos * 
				u2_dim1 + 1], ldu2, (ftnlen)3);
#line 667 "MB03ZA.f"
		    }

/*                 Update V. */

#line 671 "MB03ZA.f"
		    if (wantv) {
#line 672 "MB03ZA.f"
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &v1[pos * v1_dim1 + 1], ldv1, q, &c__4,
				 &c_b20, &dwork[1], n, (ftnlen)12, (ftnlen)12)
				;
#line 675 "MB03ZA.f"
			dlacpy_("All", n, &nb, &dwork[1], n, &v1[pos * 
				v1_dim1 + 1], ldv1, (ftnlen)3);
#line 677 "MB03ZA.f"
			dgemm_("No Transpose", "No Transpose", n, &nb, &nb, &
				c_b21, &v2[pos * v2_dim1 + 1], ldv2, q, &c__4,
				 &c_b20, &dwork[1], n, (ftnlen)12, (ftnlen)12)
				;
#line 680 "MB03ZA.f"
			dlacpy_("All", n, &nb, &dwork[1], n, &v2[pos * 
				v2_dim1 + 1], ldv2, (ftnlen)3);
#line 682 "MB03ZA.f"
		    }

#line 684 "MB03ZA.f"
		    if (nbnext == 1) {

/*                    Swap two 1-by-1 blocks. */

#line 688 "MB03ZA.f"
			pos = here;
#line 689 "MB03ZA.f"
			nb = nbnext + 1;
#line 690 "MB03ZA.f"
			dlaset_("All", &nb, &nb, &c_b20, &c_b21, q, &c__4, (
				ftnlen)3);
#line 691 "MB03ZA.f"
			dlaset_("All", &nb, &nb, &c_b20, &c_b21, z__, &c__4, (
				ftnlen)3);

#line 693 "MB03ZA.f"
			mb03wa_(&c_true, &c_true, &nbnext, &c__1, &a[pos + 
				pos * a_dim1], lda, &b[pos + pos * b_dim1], 
				ldb, q, &c__4, z__, &c__4, &ierr);

#line 697 "MB03ZA.f"
			if (ierr != 0) {
#line 698 "MB03ZA.f"
			    dwork[1] = (doublereal) wrkmin;
#line 699 "MB03ZA.f"
			    *info = 1;
#line 700 "MB03ZA.f"
			    return 0;
#line 701 "MB03ZA.f"
			}

/*                    Update rest of A. */

#line 705 "MB03ZA.f"
			if (pos > 1) {
#line 706 "MB03ZA.f"
			    i__2 = pos - 1;
#line 706 "MB03ZA.f"
			    dgemm_("No Transpose", "No Transpose", &i__2, &nb,
				     &nb, &c_b21, &a[pos * a_dim1 + 1], lda, 
				    z__, &c__4, &c_b20, &dwork[1], n, (ftnlen)
				    12, (ftnlen)12);
#line 709 "MB03ZA.f"
			    i__2 = pos - 1;
#line 709 "MB03ZA.f"
			    dlacpy_("All", &i__2, &nb, &dwork[1], n, &a[pos * 
				    a_dim1 + 1], lda, (ftnlen)3);
#line 711 "MB03ZA.f"
			}
#line 712 "MB03ZA.f"
			if (pos + nb <= *n) {
#line 713 "MB03ZA.f"
			    i__2 = *n - pos - nb + 1;
#line 713 "MB03ZA.f"
			    dgemm_("Transpose", "No Transpose", &nb, &i__2, &
				    nb, &c_b21, q, &c__4, &a[pos + (pos + nb) 
				    * a_dim1], lda, &c_b20, &dwork[1], &nb, (
				    ftnlen)9, (ftnlen)12);
#line 717 "MB03ZA.f"
			    i__2 = *n - pos - nb + 1;
#line 717 "MB03ZA.f"
			    dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &a[pos 
				    + (pos + nb) * a_dim1], lda, (ftnlen)3);
#line 719 "MB03ZA.f"
			}

/*                    Update rest of B. */

#line 723 "MB03ZA.f"
			if (pos > 1) {
#line 724 "MB03ZA.f"
			    i__2 = pos - 1;
#line 724 "MB03ZA.f"
			    dgemm_("No Transpose", "No Transpose", &i__2, &nb,
				     &nb, &c_b21, &b[pos * b_dim1 + 1], ldb, 
				    q, &c__4, &c_b20, &dwork[1], n, (ftnlen)
				    12, (ftnlen)12);
#line 727 "MB03ZA.f"
			    i__2 = pos - 1;
#line 727 "MB03ZA.f"
			    dlacpy_("All", &i__2, &nb, &dwork[1], n, &b[pos * 
				    b_dim1 + 1], ldb, (ftnlen)3);
#line 729 "MB03ZA.f"
			}
#line 730 "MB03ZA.f"
			if (pos + nb <= *n) {
#line 731 "MB03ZA.f"
			    i__2 = *n - pos - nb + 1;
#line 731 "MB03ZA.f"
			    dgemm_("Transpose", "No Transpose", &nb, &i__2, &
				    nb, &c_b21, z__, &c__4, &b[pos + (pos + 
				    nb) * b_dim1], ldb, &c_b20, &dwork[1], &
				    nb, (ftnlen)9, (ftnlen)12);
#line 735 "MB03ZA.f"
			    i__2 = *n - pos - nb + 1;
#line 735 "MB03ZA.f"
			    dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &b[pos 
				    + (pos + nb) * b_dim1], ldb, (ftnlen)3);
#line 737 "MB03ZA.f"
			}

/*                    Update C. */

#line 741 "MB03ZA.f"
			if (wantc) {
#line 742 "MB03ZA.f"
			    dgemm_("No Transpose", "No Transpose", n, &nb, &
				    nb, &c_b21, &c__[pos * c_dim1 + 1], ldc, 
				    q, &c__4, &c_b20, &dwork[1], n, (ftnlen)
				    12, (ftnlen)12);
#line 745 "MB03ZA.f"
			    dlacpy_("All", n, &nb, &dwork[1], n, &c__[pos * 
				    c_dim1 + 1], ldc, (ftnlen)3);
#line 747 "MB03ZA.f"
			    dgemm_("Transpose", "No Transpose", &nb, n, &nb, &
				    c_b21, z__, &c__4, &c__[pos + c_dim1], 
				    ldc, &c_b20, &dwork[1], &nb, (ftnlen)9, (
				    ftnlen)12);
#line 750 "MB03ZA.f"
			    dlacpy_("All", &nb, n, &dwork[1], &nb, &c__[pos + 
				    c_dim1], ldc, (ftnlen)3);
#line 752 "MB03ZA.f"
			}

/*                    Update U. */

#line 756 "MB03ZA.f"
			if (wantu) {
#line 757 "MB03ZA.f"
			    dgemm_("No Transpose", "No Transpose", n, &nb, &
				    nb, &c_b21, &u1[pos * u1_dim1 + 1], ldu1, 
				    z__, &c__4, &c_b20, &dwork[1], n, (ftnlen)
				    12, (ftnlen)12);
#line 760 "MB03ZA.f"
			    dlacpy_("All", n, &nb, &dwork[1], n, &u1[pos * 
				    u1_dim1 + 1], ldu1, (ftnlen)3);
#line 762 "MB03ZA.f"
			    dgemm_("No Transpose", "No Transpose", n, &nb, &
				    nb, &c_b21, &u2[pos * u2_dim1 + 1], ldu2, 
				    z__, &c__4, &c_b20, &dwork[1], n, (ftnlen)
				    12, (ftnlen)12);
#line 765 "MB03ZA.f"
			    dlacpy_("All", n, &nb, &dwork[1], n, &u2[pos * 
				    u2_dim1 + 1], ldu2, (ftnlen)3);
#line 767 "MB03ZA.f"
			}

/*                    Update V. */

#line 771 "MB03ZA.f"
			if (wantv) {
#line 772 "MB03ZA.f"
			    dgemm_("No Transpose", "No Transpose", n, &nb, &
				    nb, &c_b21, &v1[pos * v1_dim1 + 1], ldv1, 
				    q, &c__4, &c_b20, &dwork[1], n, (ftnlen)
				    12, (ftnlen)12);
#line 775 "MB03ZA.f"
			    dlacpy_("All", n, &nb, &dwork[1], n, &v1[pos * 
				    v1_dim1 + 1], ldv1, (ftnlen)3);
#line 777 "MB03ZA.f"
			    dgemm_("No Transpose", "No Transpose", n, &nb, &
				    nb, &c_b21, &v2[pos * v2_dim1 + 1], ldv2, 
				    q, &c__4, &c_b20, &dwork[1], n, (ftnlen)
				    12, (ftnlen)12);
#line 780 "MB03ZA.f"
			    dlacpy_("All", n, &nb, &dwork[1], n, &v2[pos * 
				    v2_dim1 + 1], ldv2, (ftnlen)3);
#line 782 "MB03ZA.f"
			}

#line 784 "MB03ZA.f"
			--here;
#line 785 "MB03ZA.f"
		    } else {

/*                    Recompute NBNEXT in case 2-by-2 split. */

#line 789 "MB03ZA.f"
			if (a[here + (here - 1) * a_dim1] == 0.) {
#line 789 "MB03ZA.f"
			    nbnext = 1;
#line 789 "MB03ZA.f"
			}

#line 792 "MB03ZA.f"
			if (nbnext == 2) {

/*                       2-by-2 block did not split. */

#line 796 "MB03ZA.f"
			    pos = here - 1;
#line 797 "MB03ZA.f"
			    nb = 3;
#line 798 "MB03ZA.f"
			    dlaset_("All", &nb, &nb, &c_b20, &c_b21, q, &c__4,
				     (ftnlen)3);
#line 799 "MB03ZA.f"
			    dlaset_("All", &nb, &nb, &c_b20, &c_b21, z__, &
				    c__4, (ftnlen)3);

#line 801 "MB03ZA.f"
			    mb03wa_(&c_true, &c_true, &c__2, &c__1, &a[pos + 
				    pos * a_dim1], lda, &b[pos + pos * b_dim1]
				    , ldb, q, &c__4, z__, &c__4, &ierr);

#line 805 "MB03ZA.f"
			    if (ierr != 0) {
#line 806 "MB03ZA.f"
				dwork[1] = (doublereal) wrkmin;
#line 807 "MB03ZA.f"
				*info = 1;
#line 808 "MB03ZA.f"
				return 0;
#line 809 "MB03ZA.f"
			    }

/*                       Update rest of A. */

#line 813 "MB03ZA.f"
			    if (pos > 1) {
#line 814 "MB03ZA.f"
				i__2 = pos - 1;
#line 814 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", &i__2, 
					&nb, &nb, &c_b21, &a[pos * a_dim1 + 1]
					, lda, z__, &c__4, &c_b20, &dwork[1], 
					n, (ftnlen)12, (ftnlen)12);
#line 817 "MB03ZA.f"
				i__2 = pos - 1;
#line 817 "MB03ZA.f"
				dlacpy_("All", &i__2, &nb, &dwork[1], n, &a[
					pos * a_dim1 + 1], lda, (ftnlen)3);
#line 819 "MB03ZA.f"
			    }
#line 820 "MB03ZA.f"
			    if (pos + nb <= *n) {
#line 821 "MB03ZA.f"
				i__2 = *n - pos - nb + 1;
#line 821 "MB03ZA.f"
				dgemm_("Transpose", "No Transpose", &nb, &
					i__2, &nb, &c_b21, q, &c__4, &a[pos + 
					(pos + nb) * a_dim1], lda, &c_b20, &
					dwork[1], &nb, (ftnlen)9, (ftnlen)12);
#line 825 "MB03ZA.f"
				i__2 = *n - pos - nb + 1;
#line 825 "MB03ZA.f"
				dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &a[
					pos + (pos + nb) * a_dim1], lda, (
					ftnlen)3);
#line 827 "MB03ZA.f"
			    }

/*                       Update rest of B. */

#line 831 "MB03ZA.f"
			    if (pos > 1) {
#line 832 "MB03ZA.f"
				i__2 = pos - 1;
#line 832 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", &i__2, 
					&nb, &nb, &c_b21, &b[pos * b_dim1 + 1]
					, ldb, q, &c__4, &c_b20, &dwork[1], n,
					 (ftnlen)12, (ftnlen)12);
#line 835 "MB03ZA.f"
				i__2 = pos - 1;
#line 835 "MB03ZA.f"
				dlacpy_("All", &i__2, &nb, &dwork[1], n, &b[
					pos * b_dim1 + 1], ldb, (ftnlen)3);
#line 837 "MB03ZA.f"
			    }
#line 838 "MB03ZA.f"
			    if (pos + nb <= *n) {
#line 839 "MB03ZA.f"
				i__2 = *n - pos - nb + 1;
#line 839 "MB03ZA.f"
				dgemm_("Transpose", "No Transpose", &nb, &
					i__2, &nb, &c_b21, z__, &c__4, &b[pos 
					+ (pos + nb) * b_dim1], ldb, &c_b20, &
					dwork[1], &nb, (ftnlen)9, (ftnlen)12);
#line 843 "MB03ZA.f"
				i__2 = *n - pos - nb + 1;
#line 843 "MB03ZA.f"
				dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &b[
					pos + (pos + nb) * b_dim1], ldb, (
					ftnlen)3);
#line 845 "MB03ZA.f"
			    }

/*                       Update C. */

#line 849 "MB03ZA.f"
			    if (wantc) {
#line 850 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &c__[pos * c_dim1 + 1], 
					ldc, q, &c__4, &c_b20, &dwork[1], n, (
					ftnlen)12, (ftnlen)12);
#line 853 "MB03ZA.f"
				dlacpy_("All", n, &nb, &dwork[1], n, &c__[pos 
					* c_dim1 + 1], ldc, (ftnlen)3);
#line 855 "MB03ZA.f"
				dgemm_("Transpose", "No Transpose", &nb, n, &
					nb, &c_b21, z__, &c__4, &c__[pos + 
					c_dim1], ldc, &c_b20, &dwork[1], &nb, 
					(ftnlen)9, (ftnlen)12);
#line 858 "MB03ZA.f"
				dlacpy_("All", &nb, n, &dwork[1], &nb, &c__[
					pos + c_dim1], ldc, (ftnlen)3);
#line 860 "MB03ZA.f"
			    }

/*                       Update U. */

#line 864 "MB03ZA.f"
			    if (wantu) {
#line 865 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &u1[pos * u1_dim1 + 1], 
					ldu1, z__, &c__4, &c_b20, &dwork[1], 
					n, (ftnlen)12, (ftnlen)12);
#line 868 "MB03ZA.f"
				dlacpy_("All", n, &nb, &dwork[1], n, &u1[pos *
					 u1_dim1 + 1], ldu1, (ftnlen)3);
#line 870 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &u2[pos * u2_dim1 + 1], 
					ldu2, z__, &c__4, &c_b20, &dwork[1], 
					n, (ftnlen)12, (ftnlen)12);
#line 873 "MB03ZA.f"
				dlacpy_("All", n, &nb, &dwork[1], n, &u2[pos *
					 u2_dim1 + 1], ldu2, (ftnlen)3);
#line 875 "MB03ZA.f"
			    }

/*                       Update V. */

#line 879 "MB03ZA.f"
			    if (wantv) {
#line 880 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &v1[pos * v1_dim1 + 1], 
					ldv1, q, &c__4, &c_b20, &dwork[1], n, 
					(ftnlen)12, (ftnlen)12);
#line 883 "MB03ZA.f"
				dlacpy_("All", n, &nb, &dwork[1], n, &v1[pos *
					 v1_dim1 + 1], ldv1, (ftnlen)3);
#line 885 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &v2[pos * v2_dim1 + 1], 
					ldv2, q, &c__4, &c_b20, &dwork[1], n, 
					(ftnlen)12, (ftnlen)12);
#line 888 "MB03ZA.f"
				dlacpy_("All", n, &nb, &dwork[1], n, &v2[pos *
					 v2_dim1 + 1], ldv2, (ftnlen)3);
#line 890 "MB03ZA.f"
			    }

#line 892 "MB03ZA.f"
			    here += -2;
#line 893 "MB03ZA.f"
			} else {

/*                       2-by-2 block did split. */

#line 897 "MB03ZA.f"
			    pos = here;
#line 898 "MB03ZA.f"
			    nb = 2;
#line 899 "MB03ZA.f"
			    dlaset_("All", &nb, &nb, &c_b20, &c_b21, q, &c__4,
				     (ftnlen)3);
#line 900 "MB03ZA.f"
			    dlaset_("All", &nb, &nb, &c_b20, &c_b21, z__, &
				    c__4, (ftnlen)3);

#line 902 "MB03ZA.f"
			    mb03wa_(&c_true, &c_true, &c__2, &c__1, &a[pos + 
				    pos * a_dim1], lda, &b[pos + pos * b_dim1]
				    , ldb, q, &c__4, z__, &c__4, &ierr);

#line 906 "MB03ZA.f"
			    if (ierr != 0) {
#line 907 "MB03ZA.f"
				dwork[1] = (doublereal) wrkmin;
#line 908 "MB03ZA.f"
				*info = 1;
#line 909 "MB03ZA.f"
				return 0;
#line 910 "MB03ZA.f"
			    }

/*                       Update rest of A. */

#line 914 "MB03ZA.f"
			    if (pos > 1) {
#line 915 "MB03ZA.f"
				i__2 = pos - 1;
#line 915 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", &i__2, 
					&nb, &nb, &c_b21, &a[pos * a_dim1 + 1]
					, lda, z__, &c__4, &c_b20, &dwork[1], 
					n, (ftnlen)12, (ftnlen)12);
#line 918 "MB03ZA.f"
				i__2 = pos - 1;
#line 918 "MB03ZA.f"
				dlacpy_("All", &i__2, &nb, &dwork[1], n, &a[
					pos * a_dim1 + 1], lda, (ftnlen)3);
#line 920 "MB03ZA.f"
			    }
#line 921 "MB03ZA.f"
			    if (pos + nb <= *n) {
#line 922 "MB03ZA.f"
				i__2 = *n - pos - nb + 1;
#line 922 "MB03ZA.f"
				dgemm_("Transpose", "No Transpose", &nb, &
					i__2, &nb, &c_b21, q, &c__4, &a[pos + 
					(pos + nb) * a_dim1], lda, &c_b20, &
					dwork[1], &nb, (ftnlen)9, (ftnlen)12);
#line 926 "MB03ZA.f"
				i__2 = *n - pos - nb + 1;
#line 926 "MB03ZA.f"
				dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &a[
					pos + (pos + nb) * a_dim1], lda, (
					ftnlen)3);
#line 928 "MB03ZA.f"
			    }

/*                       Update rest of B. */

#line 932 "MB03ZA.f"
			    if (pos > 1) {
#line 933 "MB03ZA.f"
				i__2 = pos - 1;
#line 933 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", &i__2, 
					&nb, &nb, &c_b21, &b[pos * b_dim1 + 1]
					, ldb, q, &c__4, &c_b20, &dwork[1], n,
					 (ftnlen)12, (ftnlen)12);
#line 936 "MB03ZA.f"
				i__2 = pos - 1;
#line 936 "MB03ZA.f"
				dlacpy_("All", &i__2, &nb, &dwork[1], n, &b[
					pos * b_dim1 + 1], ldb, (ftnlen)3);
#line 938 "MB03ZA.f"
			    }
#line 939 "MB03ZA.f"
			    if (pos + nb <= *n) {
#line 940 "MB03ZA.f"
				i__2 = *n - pos - nb + 1;
#line 940 "MB03ZA.f"
				dgemm_("Transpose", "No Transpose", &nb, &
					i__2, &nb, &c_b21, z__, &c__4, &b[pos 
					+ (pos + nb) * b_dim1], ldb, &c_b20, &
					dwork[1], &nb, (ftnlen)9, (ftnlen)12);
#line 944 "MB03ZA.f"
				i__2 = *n - pos - nb + 1;
#line 944 "MB03ZA.f"
				dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &b[
					pos + (pos + nb) * b_dim1], ldb, (
					ftnlen)3);
#line 946 "MB03ZA.f"
			    }

/*                       Update C. */

#line 950 "MB03ZA.f"
			    if (wantc) {
#line 951 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &c__[pos * c_dim1 + 1], 
					ldc, q, &c__4, &c_b20, &dwork[1], n, (
					ftnlen)12, (ftnlen)12);
#line 954 "MB03ZA.f"
				dlacpy_("All", n, &nb, &dwork[1], n, &c__[pos 
					* c_dim1 + 1], ldc, (ftnlen)3);
#line 956 "MB03ZA.f"
				dgemm_("Transpose", "No Transpose", &nb, n, &
					nb, &c_b21, z__, &c__4, &c__[pos + 
					c_dim1], ldc, &c_b20, &dwork[1], &nb, 
					(ftnlen)9, (ftnlen)12);
#line 959 "MB03ZA.f"
				dlacpy_("All", &nb, n, &dwork[1], &nb, &c__[
					pos + c_dim1], ldc, (ftnlen)3);
#line 961 "MB03ZA.f"
			    }

/*                       Update U. */

#line 965 "MB03ZA.f"
			    if (wantu) {
#line 966 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &u1[pos * u1_dim1 + 1], 
					ldu1, z__, &c__4, &c_b20, &dwork[1], 
					n, (ftnlen)12, (ftnlen)12);
#line 969 "MB03ZA.f"
				dlacpy_("All", n, &nb, &dwork[1], n, &u1[pos *
					 u1_dim1 + 1], ldu1, (ftnlen)3);
#line 971 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &u2[pos * u2_dim1 + 1], 
					ldu2, z__, &c__4, &c_b20, &dwork[1], 
					n, (ftnlen)12, (ftnlen)12);
#line 974 "MB03ZA.f"
				dlacpy_("All", n, &nb, &dwork[1], n, &u2[pos *
					 u2_dim1 + 1], ldu2, (ftnlen)3);
#line 976 "MB03ZA.f"
			    }

/*                       Update V. */

#line 980 "MB03ZA.f"
			    if (wantv) {
#line 981 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &v1[pos * v1_dim1 + 1], 
					ldv1, q, &c__4, &c_b20, &dwork[1], n, 
					(ftnlen)12, (ftnlen)12);
#line 984 "MB03ZA.f"
				dlacpy_("All", n, &nb, &dwork[1], n, &v1[pos *
					 v1_dim1 + 1], ldv1, (ftnlen)3);
#line 986 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &v2[pos * v2_dim1 + 1], 
					ldv2, q, &c__4, &c_b20, &dwork[1], n, 
					(ftnlen)12, (ftnlen)12);
#line 989 "MB03ZA.f"
				dlacpy_("All", n, &nb, &dwork[1], n, &v2[pos *
					 v2_dim1 + 1], ldv2, (ftnlen)3);
#line 991 "MB03ZA.f"
			    }

#line 993 "MB03ZA.f"
			    pos = here - 1;
#line 994 "MB03ZA.f"
			    nb = 2;
#line 995 "MB03ZA.f"
			    dlaset_("All", &nb, &nb, &c_b20, &c_b21, q, &c__4,
				     (ftnlen)3);
#line 996 "MB03ZA.f"
			    dlaset_("All", &nb, &nb, &c_b20, &c_b21, z__, &
				    c__4, (ftnlen)3);

#line 998 "MB03ZA.f"
			    mb03wa_(&c_true, &c_true, &c__2, &c__1, &a[pos + 
				    pos * a_dim1], lda, &b[pos + pos * b_dim1]
				    , ldb, q, &c__4, z__, &c__4, &ierr);

#line 1002 "MB03ZA.f"
			    if (ierr != 0) {
#line 1003 "MB03ZA.f"
				dwork[1] = (doublereal) wrkmin;
#line 1004 "MB03ZA.f"
				*info = 1;
#line 1005 "MB03ZA.f"
				return 0;
#line 1006 "MB03ZA.f"
			    }

/*                       Update rest of A. */

#line 1010 "MB03ZA.f"
			    if (pos > 1) {
#line 1011 "MB03ZA.f"
				i__2 = pos - 1;
#line 1011 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", &i__2, 
					&nb, &nb, &c_b21, &a[pos * a_dim1 + 1]
					, lda, z__, &c__4, &c_b20, &dwork[1], 
					n, (ftnlen)12, (ftnlen)12);
#line 1014 "MB03ZA.f"
				i__2 = pos - 1;
#line 1014 "MB03ZA.f"
				dlacpy_("All", &i__2, &nb, &dwork[1], n, &a[
					pos * a_dim1 + 1], lda, (ftnlen)3);
#line 1016 "MB03ZA.f"
			    }
#line 1017 "MB03ZA.f"
			    if (pos + nb <= *n) {
#line 1018 "MB03ZA.f"
				i__2 = *n - pos - nb + 1;
#line 1018 "MB03ZA.f"
				dgemm_("Transpose", "No Transpose", &nb, &
					i__2, &nb, &c_b21, q, &c__4, &a[pos + 
					(pos + nb) * a_dim1], lda, &c_b20, &
					dwork[1], &nb, (ftnlen)9, (ftnlen)12);
#line 1022 "MB03ZA.f"
				i__2 = *n - pos - nb + 1;
#line 1022 "MB03ZA.f"
				dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &a[
					pos + (pos + nb) * a_dim1], lda, (
					ftnlen)3);
#line 1024 "MB03ZA.f"
			    }

/*                       Update rest of B. */

#line 1028 "MB03ZA.f"
			    if (pos > 1) {
#line 1029 "MB03ZA.f"
				i__2 = pos - 1;
#line 1029 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", &i__2, 
					&nb, &nb, &c_b21, &b[pos * b_dim1 + 1]
					, ldb, q, &c__4, &c_b20, &dwork[1], n,
					 (ftnlen)12, (ftnlen)12);
#line 1032 "MB03ZA.f"
				i__2 = pos - 1;
#line 1032 "MB03ZA.f"
				dlacpy_("All", &i__2, &nb, &dwork[1], n, &b[
					pos * b_dim1 + 1], ldb, (ftnlen)3);
#line 1034 "MB03ZA.f"
			    }
#line 1035 "MB03ZA.f"
			    if (pos + nb <= *n) {
#line 1036 "MB03ZA.f"
				i__2 = *n - pos - nb + 1;
#line 1036 "MB03ZA.f"
				dgemm_("Transpose", "No Transpose", &nb, &
					i__2, &nb, &c_b21, z__, &c__4, &b[pos 
					+ (pos + nb) * b_dim1], ldb, &c_b20, &
					dwork[1], &nb, (ftnlen)9, (ftnlen)12);
#line 1040 "MB03ZA.f"
				i__2 = *n - pos - nb + 1;
#line 1040 "MB03ZA.f"
				dlacpy_("All", &nb, &i__2, &dwork[1], &nb, &b[
					pos + (pos + nb) * b_dim1], ldb, (
					ftnlen)3);
#line 1042 "MB03ZA.f"
			    }

/*                       Update C. */

#line 1046 "MB03ZA.f"
			    if (wantc) {
#line 1047 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &c__[pos * c_dim1 + 1], 
					ldc, q, &c__4, &c_b20, &dwork[1], n, (
					ftnlen)12, (ftnlen)12);
#line 1050 "MB03ZA.f"
				dlacpy_("All", n, &nb, &dwork[1], n, &c__[pos 
					* c_dim1 + 1], ldc, (ftnlen)3);
#line 1052 "MB03ZA.f"
				dgemm_("Transpose", "No Transpose", &nb, n, &
					nb, &c_b21, z__, &c__4, &c__[pos + 
					c_dim1], ldc, &c_b20, &dwork[1], &nb, 
					(ftnlen)9, (ftnlen)12);
#line 1055 "MB03ZA.f"
				dlacpy_("All", &nb, n, &dwork[1], &nb, &c__[
					pos + c_dim1], ldc, (ftnlen)3);
#line 1057 "MB03ZA.f"
			    }

/*                       Update U. */

#line 1061 "MB03ZA.f"
			    if (wantu) {
#line 1062 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &u1[pos * u1_dim1 + 1], 
					ldu1, z__, &c__4, &c_b20, &dwork[1], 
					n, (ftnlen)12, (ftnlen)12);
#line 1065 "MB03ZA.f"
				dlacpy_("All", n, &nb, &dwork[1], n, &u1[pos *
					 u1_dim1 + 1], ldu1, (ftnlen)3);
#line 1067 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &u2[pos * u2_dim1 + 1], 
					ldu2, z__, &c__4, &c_b20, &dwork[1], 
					n, (ftnlen)12, (ftnlen)12);
#line 1070 "MB03ZA.f"
				dlacpy_("All", n, &nb, &dwork[1], n, &u2[pos *
					 u2_dim1 + 1], ldu2, (ftnlen)3);
#line 1072 "MB03ZA.f"
			    }

/*                       Update V. */

#line 1076 "MB03ZA.f"
			    if (wantv) {
#line 1077 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &v1[pos * v1_dim1 + 1], 
					ldv1, q, &c__4, &c_b20, &dwork[1], n, 
					(ftnlen)12, (ftnlen)12);
#line 1080 "MB03ZA.f"
				dlacpy_("All", n, &nb, &dwork[1], n, &v1[pos *
					 v1_dim1 + 1], ldv1, (ftnlen)3);
#line 1082 "MB03ZA.f"
				dgemm_("No Transpose", "No Transpose", n, &nb,
					 &nb, &c_b21, &v2[pos * v2_dim1 + 1], 
					ldv2, q, &c__4, &c_b20, &dwork[1], n, 
					(ftnlen)12, (ftnlen)12);
#line 1085 "MB03ZA.f"
				dlacpy_("All", n, &nb, &dwork[1], n, &v2[pos *
					 v2_dim1 + 1], ldv2, (ftnlen)3);
#line 1087 "MB03ZA.f"
			    }

#line 1089 "MB03ZA.f"
			    here += -2;
#line 1090 "MB03ZA.f"
			}
#line 1091 "MB03ZA.f"
		    }
#line 1092 "MB03ZA.f"
		}

#line 1094 "MB03ZA.f"
		if (here > ilst) {
#line 1094 "MB03ZA.f"
		    goto L20;
#line 1094 "MB03ZA.f"
		}

#line 1097 "MB03ZA.f"
L30:
#line 1098 "MB03ZA.f"
		if (pair) {
#line 1098 "MB03ZA.f"
		    ++ks;
#line 1098 "MB03ZA.f"
		}
#line 1100 "MB03ZA.f"
	    }
#line 1101 "MB03ZA.f"
	}
#line 1102 "MB03ZA.f"
/* L40: */
#line 1102 "MB03ZA.f"
    }

#line 1104 "MB03ZA.f"
L50:

/*     Step 2: Compute an ordered Schur decomposition of */
/*             [ 0, -A11; B11, 0 ]. */

#line 1109 "MB03ZA.f"
    if (initw) {
#line 1109 "MB03ZA.f"
	i__1 = *m << 1;
#line 1109 "MB03ZA.f"
	i__2 = *m << 1;
#line 1109 "MB03ZA.f"
	dlaset_("All", &i__1, &i__2, &c_b20, &c_b21, &w[w_offset], ldw, (
		ftnlen)3);
#line 1109 "MB03ZA.f"
    }
#line 1111 "MB03ZA.f"
    pwc = 1;
#line 1112 "MB03ZA.f"
    pwd = pwc + (*m << 1);
#line 1113 "MB03ZA.f"
    pw = pwd + (*m << 1);
#line 1114 "MB03ZA.f"
    pair = FALSE_;
#line 1115 "MB03ZA.f"
    nb = 1;

#line 1117 "MB03ZA.f"
    i__1 = *m;
#line 1117 "MB03ZA.f"
    for (k = 1; k <= i__1; ++k) {
#line 1118 "MB03ZA.f"
	if (pair) {
#line 1119 "MB03ZA.f"
	    pair = FALSE_;
#line 1120 "MB03ZA.f"
	    nb = 1;
#line 1121 "MB03ZA.f"
	} else {
#line 1122 "MB03ZA.f"
	    if (k < *n) {
#line 1123 "MB03ZA.f"
		if (a[k + 1 + k * a_dim1] != 0.) {
#line 1124 "MB03ZA.f"
		    pair = TRUE_;
#line 1125 "MB03ZA.f"
		    nb = 2;
#line 1126 "MB03ZA.f"
		}
#line 1127 "MB03ZA.f"
	    }
#line 1128 "MB03ZA.f"
	    pwck = pwc + (k - 1 << 1);
#line 1129 "MB03ZA.f"
	    pwdl = pwd + (k - 1 << 1);
#line 1130 "MB03ZA.f"
	    i__2 = *m - k + 1;
#line 1130 "MB03ZA.f"
	    dlaset_("All", &nb, &i__2, &c_b20, &c_b20, &dwork[pwck], &c__2, (
		    ftnlen)3);
#line 1131 "MB03ZA.f"
	    i__2 = *m - k + 1;
#line 1131 "MB03ZA.f"
	    dlacpy_("All", &nb, &i__2, &a[k + k * a_dim1], lda, &dwork[pwdl], 
		    &c__2, (ftnlen)3);
#line 1132 "MB03ZA.f"
	    i__2 = *m - k + 1;
#line 1132 "MB03ZA.f"
	    dlaset_("All", &nb, &i__2, &c_b20, &c_b20, &a[k + k * a_dim1], 
		    lda, (ftnlen)3);

#line 1134 "MB03ZA.f"
	    l = k;

/*           WHILE L >= 1 DO */

#line 1138 "MB03ZA.f"
L60:

#line 1140 "MB03ZA.f"
	    if (k == l) {

/*                 Annihilate B(k,k). */

#line 1144 "MB03ZA.f"
		nbl = nb;
#line 1145 "MB03ZA.f"
		i__2 = nb + nbl;
#line 1145 "MB03ZA.f"
		i__3 = nb + nbl;
#line 1145 "MB03ZA.f"
		dlaset_("All", &i__2, &i__3, &c_b20, &c_b20, t, &c__4, (
			ftnlen)3);
#line 1147 "MB03ZA.f"
		dlacpy_("Upper", &nbl, &nbl, &b[l + l * b_dim1], ldb, &t[nb], 
			&c__4, (ftnlen)5);
#line 1149 "MB03ZA.f"
		if (nb == 1) {
#line 1150 "MB03ZA.f"
		    dwork[pwdl] = -dwork[pwdl];
#line 1151 "MB03ZA.f"
		} else {
#line 1152 "MB03ZA.f"
		    i__2 = nb << 1;
#line 1152 "MB03ZA.f"
		    dscal_(&i__2, &c_b479, &dwork[pwdl], &c__1);
#line 1153 "MB03ZA.f"
		}
#line 1154 "MB03ZA.f"
		dlacpy_("All", &nb, &nb, &dwork[pwdl], &c__2, &t[(nb + 1 << 2)
			 - 4], &c__4, (ftnlen)3);
#line 1156 "MB03ZA.f"
	    } else {

/*                 Annihilate B(l,k). */

#line 1160 "MB03ZA.f"
		i__2 = nbl + nb;
#line 1160 "MB03ZA.f"
		i__3 = nbl + nb;
#line 1160 "MB03ZA.f"
		dlaset_("All", &i__2, &i__3, &c_b20, &c_b20, t, &c__4, (
			ftnlen)3);
#line 1162 "MB03ZA.f"
		dlacpy_("All", &nbl, &nbl, &a[l + l * a_dim1], lda, t, &c__4, 
			(ftnlen)3);
#line 1163 "MB03ZA.f"
		dlacpy_("All", &nbl, &nb, &b[l + k * b_dim1], ldb, &t[(nbl + 
			1 << 2) - 4], &c__4, (ftnlen)3);
#line 1165 "MB03ZA.f"
		dlacpy_("All", &nb, &nb, &dwork[pwck], &c__2, &t[nbl + 1 + (
			nbl + 1 << 2) - 5], &c__4, (ftnlen)3);
#line 1167 "MB03ZA.f"
		pwdl = pwd + (l - 1 << 1);
#line 1168 "MB03ZA.f"
	    }

#line 1170 "MB03ZA.f"
	    i__2 = nb + nbl;
#line 1170 "MB03ZA.f"
	    dgees_("V", "Not Sorted", (L_fp)lfdum_, &i__2, t, &c__4, &mm, 
		    wrnew, winew, q, &c__4, dw12, &c__12, ldum, &ierr, (
		    ftnlen)1, (ftnlen)10);
#line 1173 "MB03ZA.f"
	    if (ierr != 0) {
#line 1174 "MB03ZA.f"
		dwork[1] = (doublereal) wrkmin;
#line 1175 "MB03ZA.f"
		*info = 3;
#line 1176 "MB03ZA.f"
		return 0;
#line 1177 "MB03ZA.f"
	    }

/*              Reorder Schur form. */

#line 1181 "MB03ZA.f"
	    mm = 0;
#line 1182 "MB03ZA.f"
	    i__2 = nb + nbl;
#line 1182 "MB03ZA.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 1183 "MB03ZA.f"
		if (wrnew[i__ - 1] > 0.) {
#line 1184 "MB03ZA.f"
		    ++mm;
#line 1185 "MB03ZA.f"
		    selnew[i__ - 1] = TRUE_;
#line 1186 "MB03ZA.f"
		} else {
#line 1187 "MB03ZA.f"
		    selnew[i__ - 1] = FALSE_;
#line 1188 "MB03ZA.f"
		}
#line 1189 "MB03ZA.f"
/* L70: */
#line 1189 "MB03ZA.f"
	    }
#line 1190 "MB03ZA.f"
	    if (mm < nb) {
#line 1191 "MB03ZA.f"
		dwork[1] = (doublereal) wrkmin;
#line 1192 "MB03ZA.f"
		*info = 4;
#line 1193 "MB03ZA.f"
		return 0;
#line 1194 "MB03ZA.f"
	    }
#line 1195 "MB03ZA.f"
	    i__2 = nb + nbl;
#line 1195 "MB03ZA.f"
	    dtrsen_("None", "V", selnew, &i__2, t, &c__4, q, &c__4, wrnew, 
		    winew, &mm, &temp, &temp, dw12, &c__4, idum, &c__1, &ierr,
		     (ftnlen)4, (ftnlen)1);
#line 1198 "MB03ZA.f"
	    if (ierr != 0) {
#line 1199 "MB03ZA.f"
		dwork[1] = (doublereal) wrkmin;
#line 1200 "MB03ZA.f"
		*info = 2;
#line 1201 "MB03ZA.f"
		return 0;
#line 1202 "MB03ZA.f"
	    }

/*              Permute Q if necessary. */

#line 1206 "MB03ZA.f"
	    if (k != l) {
#line 1207 "MB03ZA.f"
		i__2 = nb + nbl;
#line 1207 "MB03ZA.f"
		dlacpy_("All", &nbl, &i__2, q, &c__4, &z__[nb], &c__4, (
			ftnlen)3);
#line 1209 "MB03ZA.f"
		i__2 = nb + nbl;
#line 1209 "MB03ZA.f"
		dlacpy_("All", &nb, &i__2, &q[nbl], &c__4, z__, &c__4, (
			ftnlen)3);
#line 1211 "MB03ZA.f"
		i__2 = nb + nbl;
#line 1211 "MB03ZA.f"
		i__3 = nb + nbl;
#line 1211 "MB03ZA.f"
		dlacpy_("All", &i__2, &i__3, z__, &c__4, q, &c__4, (ftnlen)3);
#line 1212 "MB03ZA.f"
	    }

/*              Update "diagonal" blocks. */

#line 1216 "MB03ZA.f"
	    dlacpy_("All", &nb, &nb, t, &c__4, &dwork[pwck], &c__2, (ftnlen)3)
		    ;
#line 1217 "MB03ZA.f"
	    dlacpy_("All", &nb, &nbl, &t[(nb + 1 << 2) - 4], &c__4, &dwork[
		    pwdl], &c__2, (ftnlen)3);
#line 1219 "MB03ZA.f"
	    if (nb == 1) {
#line 1220 "MB03ZA.f"
		dscal_(&nbl, &c_b479, &dwork[pwdl], &c__2);
#line 1221 "MB03ZA.f"
	    } else {
#line 1222 "MB03ZA.f"
		i__2 = nbl << 1;
#line 1222 "MB03ZA.f"
		dscal_(&i__2, &c_b479, &dwork[pwdl], &c__1);
#line 1223 "MB03ZA.f"
	    }
#line 1224 "MB03ZA.f"
	    dlacpy_("All", &nbl, &nbl, &t[nb + 1 + (nb + 1 << 2) - 5], &c__4, 
		    &a[l + l * a_dim1], lda, (ftnlen)3);

/*              Update block columns of A and B. */

#line 1229 "MB03ZA.f"
	    len = l - 1;
#line 1230 "MB03ZA.f"
	    if (len > 0) {
#line 1231 "MB03ZA.f"
		dgemm_("No Transpose", "No Transpose", &len, &nb, &nb, &c_b21,
			 &b[k * b_dim1 + 1], ldb, q, &c__4, &c_b20, &dwork[pw]
			, m, (ftnlen)12, (ftnlen)12);
#line 1234 "MB03ZA.f"
		dgemm_("No Transpose", "No Transpose", &len, &nbl, &nb, &
			c_b21, &b[k * b_dim1 + 1], ldb, &q[(nb + 1 << 2) - 4],
			 &c__4, &c_b20, &dwork[pw + (*m << 1)], m, (ftnlen)12,
			 (ftnlen)12);
#line 1237 "MB03ZA.f"
		dgemm_("No Transpose", "No Transpose", &len, &nb, &nbl, &
			c_b21, &a[l * a_dim1 + 1], lda, &q[nb], &c__4, &c_b21,
			 &dwork[pw], m, (ftnlen)12, (ftnlen)12);
#line 1240 "MB03ZA.f"
		dlacpy_("All", &len, &nb, &dwork[pw], m, &b[k * b_dim1 + 1], 
			ldb, (ftnlen)3);
#line 1242 "MB03ZA.f"
		dgemm_("No Transpose", "No Transpose", &len, &nbl, &nbl, &
			c_b21, &a[l * a_dim1 + 1], lda, &q[nb + 1 + (nb + 1 <<
			 2) - 5], &c__4, &c_b21, &dwork[pw + (*m << 1)], m, (
			ftnlen)12, (ftnlen)12);
#line 1245 "MB03ZA.f"
		dlacpy_("All", &len, &nbl, &dwork[pw + (*m << 1)], m, &a[l * 
			a_dim1 + 1], lda, (ftnlen)3);
#line 1247 "MB03ZA.f"
	    }

/*              Update block column of A. */

#line 1251 "MB03ZA.f"
	    len = *m - l - nbl + 1;
#line 1252 "MB03ZA.f"
	    if (len > 0) {
#line 1253 "MB03ZA.f"
		dgemm_("Transpose", "No Transpose", &nb, &len, &nb, &c_b21, q,
			 &c__4, &dwork[pwdl + (nbl << 1)], &c__2, &c_b20, &
			dwork[pw], &c__2, (ftnlen)9, (ftnlen)12);
#line 1256 "MB03ZA.f"
		dgemm_("Transpose", "No Transpose", &nbl, &len, &nb, &c_b479, 
			&q[(nb + 1 << 2) - 4], &c__4, &dwork[pwdl + (nbl << 1)
			], &c__2, &c_b20, &dwork[pw + (*m << 1)], &c__2, (
			ftnlen)9, (ftnlen)12);
#line 1259 "MB03ZA.f"
		dgemm_("Transpose", "No Transpose", &nb, &len, &nbl, &c_b479, 
			&q[nb], &c__4, &a[l + (l + nbl) * a_dim1], lda, &
			c_b21, &dwork[pw], &c__2, (ftnlen)9, (ftnlen)12);
#line 1262 "MB03ZA.f"
		dlacpy_("All", &nb, &len, &dwork[pw], &c__2, &dwork[pwdl + (
			nbl << 1)], &c__2, (ftnlen)3);
#line 1264 "MB03ZA.f"
		dgemm_("Transpose", "No Transpose", &nbl, &len, &nbl, &c_b21, 
			&q[nb + 1 + (nb + 1 << 2) - 5], &c__4, &a[l + (l + 
			nbl) * a_dim1], lda, &c_b21, &dwork[pw + (*m << 1)], &
			c__2, (ftnlen)9, (ftnlen)12);
#line 1267 "MB03ZA.f"
		dlacpy_("All", &nbl, &len, &dwork[pw + (*m << 1)], &c__2, &a[
			l + (l + nbl) * a_dim1], lda, (ftnlen)3);
#line 1269 "MB03ZA.f"
	    }

/*              Update block row of B. */

#line 1273 "MB03ZA.f"
	    len = *m - k - nb + 1;
#line 1274 "MB03ZA.f"
	    if (len > 0) {
#line 1275 "MB03ZA.f"
		dgemm_("Transpose", "No Transpose", &nb, &len, &nb, &c_b21, q,
			 &c__4, &dwork[pwck + (nb << 1)], &c__2, &c_b20, &
			dwork[pw], &c__2, (ftnlen)9, (ftnlen)12);
#line 1278 "MB03ZA.f"
		dgemm_("Transpose", "No Transpose", &nbl, &len, &nb, &c_b21, &
			q[(nb + 1 << 2) - 4], &c__4, &dwork[pwck + (nb << 1)],
			 &c__2, &c_b20, &dwork[pw + (*m << 1)], &c__2, (
			ftnlen)9, (ftnlen)12);
#line 1281 "MB03ZA.f"
		dgemm_("Transpose", "No Transpose", &nb, &len, &nbl, &c_b21, &
			q[nb], &c__4, &b[l + (k + nb) * b_dim1], ldb, &c_b21, 
			&dwork[pw], &c__2, (ftnlen)9, (ftnlen)12);
#line 1284 "MB03ZA.f"
		dlacpy_("All", &nb, &len, &dwork[pw], &c__2, &dwork[pwck + (
			nb << 1)], &c__2, (ftnlen)3);
#line 1286 "MB03ZA.f"
		dgemm_("Transpose", "No Transpose", &nbl, &len, &nbl, &c_b21, 
			&q[nb + 1 + (nb + 1 << 2) - 5], &c__4, &b[l + (k + nb)
			 * b_dim1], ldb, &c_b21, &dwork[pw + (*m << 1)], &
			c__2, (ftnlen)9, (ftnlen)12);
#line 1289 "MB03ZA.f"
		dlacpy_("All", &nbl, &len, &dwork[pw + (*m << 1)], &c__2, &b[
			l + (k + nb) * b_dim1], ldb, (ftnlen)3);
#line 1291 "MB03ZA.f"
	    }

/*              Update W. */

#line 1295 "MB03ZA.f"
	    if (wantw) {
#line 1296 "MB03ZA.f"
		if (initw) {
#line 1297 "MB03ZA.f"
		    pos = l;
#line 1298 "MB03ZA.f"
		    len = k + nb - l;
#line 1299 "MB03ZA.f"
		} else {
#line 1300 "MB03ZA.f"
		    pos = 1;
#line 1301 "MB03ZA.f"
		    len = *m;
#line 1302 "MB03ZA.f"
		}
#line 1303 "MB03ZA.f"
		dgemm_("No Transpose", "No Transpose", &len, &nb, &nb, &c_b21,
			 &w[pos + k * w_dim1], ldw, q, &c__4, &c_b20, &dwork[
			pw], m, (ftnlen)12, (ftnlen)12);
#line 1306 "MB03ZA.f"
		dgemm_("No Transpose", "No Transpose", &len, &nbl, &nb, &
			c_b21, &w[pos + k * w_dim1], ldw, &q[(nb + 1 << 2) - 
			4], &c__4, &c_b20, &dwork[pw + (*m << 1)], m, (ftnlen)
			12, (ftnlen)12);
#line 1309 "MB03ZA.f"
		dgemm_("No Transpose", "No Transpose", &len, &nb, &nbl, &
			c_b21, &w[pos + (*m + l) * w_dim1], ldw, &q[nb], &
			c__4, &c_b21, &dwork[pw], m, (ftnlen)12, (ftnlen)12);
#line 1312 "MB03ZA.f"
		dlacpy_("All", &len, &nb, &dwork[pw], m, &w[pos + k * w_dim1],
			 ldw, (ftnlen)3);
#line 1314 "MB03ZA.f"
		dgemm_("No Transpose", "No Transpose", &len, &nbl, &nbl, &
			c_b21, &w[pos + (*m + l) * w_dim1], ldw, &q[nb + 1 + (
			nb + 1 << 2) - 5], &c__4, &c_b21, &dwork[pw + (*m << 
			1)], m, (ftnlen)12, (ftnlen)12);
#line 1317 "MB03ZA.f"
		dlacpy_("All", &len, &nbl, &dwork[pw + (*m << 1)], m, &w[pos 
			+ (*m + l) * w_dim1], ldw, (ftnlen)3);

#line 1320 "MB03ZA.f"
		dgemm_("No Transpose", "No Transpose", &len, &nb, &nb, &c_b21,
			 &w[*m + pos + k * w_dim1], ldw, q, &c__4, &c_b20, &
			dwork[pw], m, (ftnlen)12, (ftnlen)12);
#line 1323 "MB03ZA.f"
		dgemm_("No Transpose", "No Transpose", &len, &nbl, &nb, &
			c_b21, &w[*m + pos + k * w_dim1], ldw, &q[(nb + 1 << 
			2) - 4], &c__4, &c_b20, &dwork[pw + (*m << 1)], m, (
			ftnlen)12, (ftnlen)12);
#line 1326 "MB03ZA.f"
		dgemm_("No Transpose", "No Transpose", &len, &nb, &nbl, &
			c_b21, &w[*m + pos + (*m + l) * w_dim1], ldw, &q[nb], 
			&c__4, &c_b21, &dwork[pw], m, (ftnlen)12, (ftnlen)12);
#line 1329 "MB03ZA.f"
		dlacpy_("All", &len, &nb, &dwork[pw], m, &w[*m + pos + k * 
			w_dim1], ldw, (ftnlen)3);
#line 1331 "MB03ZA.f"
		dgemm_("No Transpose", "No Transpose", &len, &nbl, &nbl, &
			c_b21, &w[*m + pos + (*m + l) * w_dim1], ldw, &q[nb + 
			1 + (nb + 1 << 2) - 5], &c__4, &c_b21, &dwork[pw + (*
			m << 1)], m, (ftnlen)12, (ftnlen)12);
#line 1334 "MB03ZA.f"
		dlacpy_("All", &len, &nbl, &dwork[pw + (*m << 1)], m, &w[*m + 
			pos + (*m + l) * w_dim1], ldw, (ftnlen)3);
#line 1336 "MB03ZA.f"
	    }

#line 1338 "MB03ZA.f"
	    --l;
#line 1339 "MB03ZA.f"
	    nbl = 1;
#line 1340 "MB03ZA.f"
	    if (l > 1) {
#line 1341 "MB03ZA.f"
		if (a[l + (l - 1) * a_dim1] != 0.) {
#line 1342 "MB03ZA.f"
		    nbl = 2;
#line 1343 "MB03ZA.f"
		    --l;
#line 1344 "MB03ZA.f"
		}
#line 1345 "MB03ZA.f"
	    }

/*           END WHILE L >= 1 DO */

#line 1349 "MB03ZA.f"
	    if (l >= 1) {
#line 1349 "MB03ZA.f"
		goto L60;
#line 1349 "MB03ZA.f"
	    }

/*           Copy recomputed eigenvalues. */

#line 1354 "MB03ZA.f"
	    dcopy_(&nb, wrnew, &c__1, &wr[k], &c__1);
#line 1355 "MB03ZA.f"
	    dcopy_(&nb, winew, &c__1, &wi[k], &c__1);
#line 1356 "MB03ZA.f"
	}
#line 1357 "MB03ZA.f"
/* L80: */
#line 1357 "MB03ZA.f"
    }
#line 1358 "MB03ZA.f"
    dwork[1] = (doublereal) wrkmin;
#line 1359 "MB03ZA.f"
    return 0;
/* *** Last line of MB03ZA *** */
} /* mb03za_ */


logical lfdum_(doublereal *x, doublereal *y)
{
    /* System generated locals */
    logical ret_val;


/*     Void logical function for DGEES. */

#line 1368 "MB03ZA.f"
    ret_val = FALSE_;
#line 1369 "MB03ZA.f"
    return ret_val;
/* *** Last line of LFDUM *** */
} /* lfdum_ */

