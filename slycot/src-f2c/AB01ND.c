#line 1 "AB01ND.f"
/* AB01ND.f -- translated by f2c (version 20100827).
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

#line 1 "AB01ND.f"
/* Table of constant values */

static doublereal c_b7 = 0.;
static doublereal c_b8 = 1.;
static integer c__1 = 1;
static integer c__0 = 0;
static logical c_false = FALSE_;

/* Subroutine */ int ab01nd_(char *jobz, integer *n, integer *m, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, integer *ncont, integer 
	*indcon, integer *nblk, doublereal *z__, integer *ldz, doublereal *
	tau, doublereal *tol, integer *iwork, doublereal *dwork, integer *
	ldwork, integer *info, ftnlen jobz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, z_dim1, z_offset, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j, ni, nj, nbl, iqr, rank, itau;
    static doublereal fnrm;
    static integer mcrt, ncrt;
    static doublereal sval[3];
    extern /* Subroutine */ int mb01pd_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen);
    static logical ljobf, ljobi;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb03oy_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *);
    static doublereal anorm, bnorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical ljobz;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen), dlange_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dlapmt_(logical *, integer *, 
	    integer *, doublereal *, integer *, integer *), dorgqr_(integer *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *), dormqr_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
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

/*     To find a controllable realization for the linear time-invariant */
/*     multi-input system */

/*             dX/dt = A * X + B * U, */

/*     where A and B are N-by-N and N-by-M matrices, respectively, */
/*     which are reduced by this routine to orthogonal canonical form */
/*     using (and optionally accumulating) orthogonal similarity */
/*     transformations.  Specifically, the pair (A, B) is reduced to */
/*     the pair (Ac, Bc),  Ac = Z' * A * Z,  Bc = Z' * B,  given by */

/*             [ Acont     *    ]         [ Bcont ] */
/*        Ac = [                ],   Bc = [       ], */
/*             [   0    Auncont ]         [   0   ] */

/*        and */

/*                [ A11 A12  . . .  A1,p-1 A1p ]         [ B1 ] */
/*                [ A21 A22  . . .  A2,p-1 A2p ]         [ 0  ] */
/*                [  0  A32  . . .  A3,p-1 A3p ]         [ 0  ] */
/*        Acont = [  .   .   . . .    .     .  ],   Bc = [ .  ], */
/*                [  .   .     . .    .     .  ]         [ .  ] */
/*                [  .   .       .    .     .  ]         [ .  ] */
/*                [  0   0   . . .  Ap,p-1 App ]         [ 0  ] */

/*     where the blocks  B1, A21, ..., Ap,p-1  have full row ranks and */
/*     p is the controllability index of the pair.  The size of the */
/*     block  Auncont is equal to the dimension of the uncontrollable */
/*     subspace of the pair (A, B). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBZ    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix Z the orthogonal similarity transformations for */
/*             reducing the system, as follows: */
/*             = 'N':  Do not form Z and do not store the orthogonal */
/*                     transformations; */
/*             = 'F':  Do not form Z, but store the orthogonal */
/*                     transformations in the factored form; */
/*             = 'I':  Z is initialized to the unit matrix and the */
/*                     orthogonal transformation matrix Z is returned. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the original state-space representation, */
/*             i.e. the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs, or of columns of B.  M >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original state dynamics matrix A. */
/*             On exit, the leading NCONT-by-NCONT part contains the */
/*             upper block Hessenberg state dynamics matrix Acont in Ac, */
/*             given by Z' * A * Z, of a controllable realization for */
/*             the original system. The elements below the first block- */
/*             subdiagonal are set to zero. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input matrix B. */
/*             On exit, the leading NCONT-by-M part of this array */
/*             contains the transformed input matrix Bcont in Bc, given */
/*             by Z' * B, with all elements but the first block set to */
/*             zero. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     NCONT   (output) INTEGER */
/*             The order of the controllable state-space representation. */

/*     INDCON  (output) INTEGER */
/*             The controllability index of the controllable part of the */
/*             system representation. */

/*     NBLK    (output) INTEGER array, dimension (N) */
/*             The leading INDCON elements of this array contain the */
/*             the orders of the diagonal blocks of Acont. */

/*     Z       (output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             If JOBZ = 'I', then the leading N-by-N part of this */
/*             array contains the matrix of accumulated orthogonal */
/*             similarity transformations which reduces the given system */
/*             to orthogonal canonical form. */
/*             If JOBZ = 'F', the elements below the diagonal, with the */
/*             array TAU, represent the orthogonal transformation matrix */
/*             as a product of elementary reflectors. The transformation */
/*             matrix can then be obtained by calling the LAPACK Library */
/*             routine DORGQR. */
/*             If JOBZ = 'N', the array Z is not referenced and can be */
/*             supplied as a dummy array (i.e. set parameter LDZ = 1 and */
/*             declare this array to be Z(1,1) in the calling program). */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z. If JOBZ = 'I' or */
/*             JOBZ = 'F', LDZ >= MAX(1,N); if JOBZ = 'N', LDZ >= 1. */

/*     TAU     (output) DOUBLE PRECISION array, dimension (N) */
/*             The elements of TAU contain the scalar factors of the */
/*             elementary reflectors used in the reduction of B and A. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in rank determination when */
/*             transforming (A, B). If the user sets TOL > 0, then */
/*             the given value of TOL is used as a lower bound for the */
/*             reciprocal condition number (see the description of the */
/*             argument RCOND in the SLICOT routine MB03OD);  a */
/*             (sub)matrix whose estimated condition number is less than */
/*             1/TOL is considered to be of full rank.  If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance, */
/*             defined by  TOLDEF = N*N*EPS,  is used instead, where EPS */
/*             is the machine precision (see LAPACK Library routine */
/*             DLAMCH). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1, N, 3*M). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Matrix B is first QR-decomposed and the appropriate orthogonal */
/*     similarity transformation applied to the matrix A. Leaving the */
/*     first rank(B) states unchanged, the remaining lower left block */
/*     of A is then QR-decomposed and the new orthogonal matrix, Q1, */
/*     is also applied to the right of A to complete the similarity */
/*     transformation. By continuing in this manner, a completely */
/*     controllable state-space pair (Acont, Bcont) is found for the */
/*     given (A, B), where Acont is upper block Hessenberg with each */
/*     subdiagonal block of full row rank, and Bcont is zero apart from */
/*     its (independent) first rank(B) rows. */
/*     NOTE that the system controllability indices are easily */
/*     calculated from the dimensions of the blocks of Acont. */

/*     REFERENCES */

/*     [1] Konstantinov, M.M., Petkov, P.Hr. and Christov, N.D. */
/*         Orthogonal Invariants and Canonical Forms for Linear */
/*         Controllable Systems. */
/*         Proc. 8th IFAC World Congress, Kyoto, 1, pp. 49-54, 1981. */

/*     [2] Paige, C.C. */
/*         Properties of numerical algorithms related to computing */
/*         controllablity. */
/*         IEEE Trans. Auto. Contr., AC-26, pp. 130-138, 1981. */

/*     [3] Petkov, P.Hr., Konstantinov, M.M., Gu, D.W. and */
/*         Postlethwaite, I. */
/*         Optimal Pole Assignment Design of Linear Multi-Input Systems. */
/*         Leicester University, Report 99-11, May 1996. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations and is backward stable. */

/*     FURTHER COMMENTS */

/*     If the system matrices A and B are badly scaled, it would be */
/*     useful to scale them with SLICOT routine TB01ID, before calling */
/*     the routine. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996. */
/*     Supersedes Release 2.0 routine AB01BD by P.Hr. Petkov. */

/*     REVISIONS */

/*     January 14, 1997, June 4, 1997, February 13, 1998, */
/*     September 22, 2003, February 29, 2004. */

/*     KEYWORDS */

/*     Controllability, minimal realization, orthogonal canonical form, */
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
/*     .. */
/*     .. Executable Statements .. */

#line 258 "AB01ND.f"
    /* Parameter adjustments */
#line 258 "AB01ND.f"
    a_dim1 = *lda;
#line 258 "AB01ND.f"
    a_offset = 1 + a_dim1;
#line 258 "AB01ND.f"
    a -= a_offset;
#line 258 "AB01ND.f"
    b_dim1 = *ldb;
#line 258 "AB01ND.f"
    b_offset = 1 + b_dim1;
#line 258 "AB01ND.f"
    b -= b_offset;
#line 258 "AB01ND.f"
    --nblk;
#line 258 "AB01ND.f"
    z_dim1 = *ldz;
#line 258 "AB01ND.f"
    z_offset = 1 + z_dim1;
#line 258 "AB01ND.f"
    z__ -= z_offset;
#line 258 "AB01ND.f"
    --tau;
#line 258 "AB01ND.f"
    --iwork;
#line 258 "AB01ND.f"
    --dwork;
#line 258 "AB01ND.f"

#line 258 "AB01ND.f"
    /* Function Body */
#line 258 "AB01ND.f"
    *info = 0;
#line 259 "AB01ND.f"
    ljobf = lsame_(jobz, "F", (ftnlen)1, (ftnlen)1);
#line 260 "AB01ND.f"
    ljobi = lsame_(jobz, "I", (ftnlen)1, (ftnlen)1);
#line 261 "AB01ND.f"
    ljobz = ljobf || ljobi;

/*     Test the input scalar arguments. */

#line 265 "AB01ND.f"
    if (! ljobz && ! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 266 "AB01ND.f"
	*info = -1;
#line 267 "AB01ND.f"
    } else if (*n < 0) {
#line 268 "AB01ND.f"
	*info = -2;
#line 269 "AB01ND.f"
    } else if (*m < 0) {
#line 270 "AB01ND.f"
	*info = -3;
#line 271 "AB01ND.f"
    } else if (*lda < max(1,*n)) {
#line 272 "AB01ND.f"
	*info = -5;
#line 273 "AB01ND.f"
    } else if (*ldb < max(1,*n)) {
#line 274 "AB01ND.f"
	*info = -7;
#line 275 "AB01ND.f"
    } else if (*ldz < 1 || ljobz && *ldz < *n) {
#line 276 "AB01ND.f"
	*info = -12;
#line 277 "AB01ND.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 277 "AB01ND.f"
	i__1 = max(1,*n), i__2 = *m * 3;
#line 277 "AB01ND.f"
	if (*ldwork < max(i__1,i__2)) {
#line 278 "AB01ND.f"
	    *info = -17;
#line 279 "AB01ND.f"
	}
#line 279 "AB01ND.f"
    }

#line 281 "AB01ND.f"
    if (*info != 0) {

/*        Error return. */

#line 285 "AB01ND.f"
	i__1 = -(*info);
#line 285 "AB01ND.f"
	xerbla_("AB01ND", &i__1, (ftnlen)6);
#line 286 "AB01ND.f"
	return 0;
#line 287 "AB01ND.f"
    }

#line 289 "AB01ND.f"
    *ncont = 0;
#line 290 "AB01ND.f"
    *indcon = 0;

/*     Quick return if possible. */

#line 294 "AB01ND.f"
    if (min(*n,*m) == 0) {
#line 295 "AB01ND.f"
	if (*n > 0) {
#line 296 "AB01ND.f"
	    if (ljobi) {
#line 297 "AB01ND.f"
		dlaset_("Full", n, n, &c_b7, &c_b8, &z__[z_offset], ldz, (
			ftnlen)4);
#line 298 "AB01ND.f"
	    } else if (ljobf) {
#line 299 "AB01ND.f"
		dlaset_("Full", n, n, &c_b7, &c_b7, &z__[z_offset], ldz, (
			ftnlen)4);
#line 300 "AB01ND.f"
		dlaset_("Full", n, &c__1, &c_b7, &c_b7, &tau[1], n, (ftnlen)4)
			;
#line 301 "AB01ND.f"
	    }
#line 302 "AB01ND.f"
	}
#line 303 "AB01ND.f"
	dwork[1] = 1.;
#line 304 "AB01ND.f"
	return 0;
#line 305 "AB01ND.f"
    }

/*     Calculate the absolute norms of A and B (used for scaling). */

#line 309 "AB01ND.f"
    anorm = dlange_("M", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)1);
#line 310 "AB01ND.f"
    bnorm = dlange_("M", n, m, &b[b_offset], ldb, &dwork[1], (ftnlen)1);

/*     Return if matrix B is zero. */

#line 314 "AB01ND.f"
    if (bnorm == 0.) {
#line 315 "AB01ND.f"
	if (ljobi) {
#line 316 "AB01ND.f"
	    dlaset_("Full", n, n, &c_b7, &c_b8, &z__[z_offset], ldz, (ftnlen)
		    4);
#line 317 "AB01ND.f"
	} else if (ljobf) {
#line 318 "AB01ND.f"
	    dlaset_("Full", n, n, &c_b7, &c_b7, &z__[z_offset], ldz, (ftnlen)
		    4);
#line 319 "AB01ND.f"
	    dlaset_("Full", n, &c__1, &c_b7, &c_b7, &tau[1], n, (ftnlen)4);
#line 320 "AB01ND.f"
	}
#line 321 "AB01ND.f"
	dwork[1] = 1.;
#line 322 "AB01ND.f"
	return 0;
#line 323 "AB01ND.f"
    }

/*     Scale (if needed) the matrices A and B. */

#line 327 "AB01ND.f"
    mb01pd_("Scale", "G", n, n, &c__0, &c__0, &anorm, &c__0, &nblk[1], &a[
	    a_offset], lda, info, (ftnlen)5, (ftnlen)1);
#line 329 "AB01ND.f"
    mb01pd_("Scale", "G", n, m, &c__0, &c__0, &bnorm, &c__0, &nblk[1], &b[
	    b_offset], ldb, info, (ftnlen)5, (ftnlen)1);

/*     Compute the Frobenius norm of [ B  A ] (used for rank estimation). */

#line 334 "AB01ND.f"
    d__1 = dlange_("F", n, m, &b[b_offset], ldb, &dwork[1], (ftnlen)1);
#line 334 "AB01ND.f"
    d__2 = dlange_("F", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)1);
#line 334 "AB01ND.f"
    fnrm = dlapy2_(&d__1, &d__2);

#line 337 "AB01ND.f"
    toldef = *tol;
#line 338 "AB01ND.f"
    if (toldef <= 0.) {

/*        Use the default tolerance in controllability determination. */

#line 342 "AB01ND.f"
	toldef = (doublereal) (*n * *n) * dlamch_("EPSILON", (ftnlen)7);
#line 343 "AB01ND.f"
    }

#line 345 "AB01ND.f"
    wrkopt = 1;
#line 346 "AB01ND.f"
    ni = 0;
#line 347 "AB01ND.f"
    itau = 1;
#line 348 "AB01ND.f"
    ncrt = *n;
#line 349 "AB01ND.f"
    mcrt = *m;
#line 350 "AB01ND.f"
    iqr = 1;

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 358 "AB01ND.f"
L10:

/*        Rank-revealing QR decomposition with column pivoting. */
/*        The calculation is performed in NCRT rows of B starting from */
/*        the row IQR (initialized to 1 and then set to rank(B)+1). */
/*        Workspace: 3*MCRT. */

#line 365 "AB01ND.f"
    mb03oy_(&ncrt, &mcrt, &b[iqr + b_dim1], ldb, &toldef, &fnrm, &rank, sval, 
	    &iwork[1], &tau[itau], &dwork[1], info);

#line 368 "AB01ND.f"
    if (rank != 0) {
#line 369 "AB01ND.f"
	nj = ni;
#line 370 "AB01ND.f"
	ni = *ncont;
#line 371 "AB01ND.f"
	*ncont += rank;
#line 372 "AB01ND.f"
	++(*indcon);
#line 373 "AB01ND.f"
	nblk[*indcon] = rank;

/*           Premultiply and postmultiply the appropriate block row */
/*           and block column of A by Q' and Q, respectively. */
/*           Workspace: need   NCRT; */
/*                      prefer NCRT*NB. */

#line 380 "AB01ND.f"
	dormqr_("Left", "Transpose", &ncrt, &ncrt, &rank, &b[iqr + b_dim1], 
		ldb, &tau[itau], &a[ni + 1 + (ni + 1) * a_dim1], lda, &dwork[
		1], ldwork, info, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 383 "AB01ND.f"
	i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 383 "AB01ND.f"
	wrkopt = max(i__1,i__2);

/*           Workspace: need   N; */
/*                      prefer N*NB. */

#line 388 "AB01ND.f"
	dormqr_("Right", "No transpose", n, &ncrt, &rank, &b[iqr + b_dim1], 
		ldb, &tau[itau], &a[(ni + 1) * a_dim1 + 1], lda, &dwork[1], 
		ldwork, info, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
#line 391 "AB01ND.f"
	i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 391 "AB01ND.f"
	wrkopt = max(i__1,i__2);

/*           If required, save transformations. */

#line 395 "AB01ND.f"
	if (ljobz && ncrt > 1) {
#line 396 "AB01ND.f"
	    i__1 = ncrt - 1;
/* Computing MIN */
#line 396 "AB01ND.f"
	    i__3 = rank, i__4 = ncrt - 1;
#line 396 "AB01ND.f"
	    i__2 = min(i__3,i__4);
#line 396 "AB01ND.f"
	    dlacpy_("L", &i__1, &i__2, &b[iqr + 1 + b_dim1], ldb, &z__[ni + 2 
		    + itau * z_dim1], ldz, (ftnlen)1);
#line 398 "AB01ND.f"
	}

/*           Zero the subdiagonal elements of the current matrix. */

#line 402 "AB01ND.f"
	if (rank > 1) {
#line 402 "AB01ND.f"
	    i__1 = rank - 1;
#line 402 "AB01ND.f"
	    i__2 = rank - 1;
#line 402 "AB01ND.f"
	    dlaset_("L", &i__1, &i__2, &c_b7, &c_b7, &b[iqr + 1 + b_dim1], 
		    ldb, (ftnlen)1);
#line 402 "AB01ND.f"
	}

/*           Backward permutation of the columns of B or A. */

#line 408 "AB01ND.f"
	if (*indcon == 1) {
#line 409 "AB01ND.f"
	    dlapmt_(&c_false, &rank, m, &b[iqr + b_dim1], ldb, &iwork[1]);
#line 410 "AB01ND.f"
	    iqr = rank + 1;
#line 411 "AB01ND.f"
	} else {
#line 412 "AB01ND.f"
	    i__1 = mcrt;
#line 412 "AB01ND.f"
	    for (j = 1; j <= i__1; ++j) {
#line 413 "AB01ND.f"
		dcopy_(&rank, &b[iqr + j * b_dim1], &c__1, &a[ni + 1 + (nj + 
			iwork[j]) * a_dim1], &c__1);
#line 415 "AB01ND.f"
/* L20: */
#line 415 "AB01ND.f"
	    }
#line 416 "AB01ND.f"
	}

#line 418 "AB01ND.f"
	itau += rank;
#line 419 "AB01ND.f"
	if (rank != ncrt) {
#line 420 "AB01ND.f"
	    mcrt = rank;
#line 421 "AB01ND.f"
	    ncrt -= rank;
#line 422 "AB01ND.f"
	    dlacpy_("G", &ncrt, &mcrt, &a[*ncont + 1 + (ni + 1) * a_dim1], 
		    lda, &b[iqr + b_dim1], ldb, (ftnlen)1);
#line 424 "AB01ND.f"
	    dlaset_("G", &ncrt, &mcrt, &c_b7, &c_b7, &a[*ncont + 1 + (ni + 1) 
		    * a_dim1], lda, (ftnlen)1);
#line 426 "AB01ND.f"
	    goto L10;
#line 427 "AB01ND.f"
	}
#line 428 "AB01ND.f"
    }

/*     If required, accumulate transformations. */
/*     Workspace: need N;  prefer N*NB. */

#line 433 "AB01ND.f"
    if (ljobi) {
/* Computing MAX */
#line 434 "AB01ND.f"
	i__2 = 1, i__3 = itau - 1;
#line 434 "AB01ND.f"
	i__1 = max(i__2,i__3);
#line 434 "AB01ND.f"
	dorgqr_(n, n, &i__1, &z__[z_offset], ldz, &tau[1], &dwork[1], ldwork, 
		info);
/* Computing MAX */
#line 436 "AB01ND.f"
	i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 436 "AB01ND.f"
	wrkopt = max(i__1,i__2);
#line 437 "AB01ND.f"
    }

/*     Annihilate the trailing blocks of B. */

#line 441 "AB01ND.f"
    if (*n >= iqr) {
#line 441 "AB01ND.f"
	i__1 = *n - iqr + 1;
#line 441 "AB01ND.f"
	dlaset_("G", &i__1, m, &c_b7, &c_b7, &b[iqr + b_dim1], ldb, (ftnlen)1)
		;
#line 441 "AB01ND.f"
    }

/*     Annihilate the trailing elements of TAU, if JOBZ = 'F'. */

#line 446 "AB01ND.f"
    if (ljobf) {
#line 447 "AB01ND.f"
	i__1 = *n;
#line 447 "AB01ND.f"
	for (j = itau; j <= i__1; ++j) {
#line 448 "AB01ND.f"
	    tau[j] = 0.;
#line 449 "AB01ND.f"
/* L30: */
#line 449 "AB01ND.f"
	}
#line 450 "AB01ND.f"
    }

/*     Undo scaling of A and B. */

#line 454 "AB01ND.f"
    if (*indcon < *n) {
#line 455 "AB01ND.f"
	nbl = *indcon + 1;
#line 456 "AB01ND.f"
	nblk[nbl] = *n - *ncont;
#line 457 "AB01ND.f"
    } else {
#line 458 "AB01ND.f"
	nbl = 0;
#line 459 "AB01ND.f"
    }
#line 460 "AB01ND.f"
    mb01pd_("Undo", "H", n, n, &c__0, &c__0, &anorm, &nbl, &nblk[1], &a[
	    a_offset], lda, info, (ftnlen)4, (ftnlen)1);
#line 462 "AB01ND.f"
    mb01pd_("Undo", "G", &nblk[1], m, &c__0, &c__0, &bnorm, &c__0, &nblk[1], &
	    b[b_offset], ldb, info, (ftnlen)4, (ftnlen)1);

/*     Set optimal workspace dimension. */

#line 467 "AB01ND.f"
    dwork[1] = (doublereal) wrkopt;
#line 468 "AB01ND.f"
    return 0;
/* *** Last line of AB01ND *** */
} /* ab01nd_ */

