#line 1 "MB03XU.f"
/* MB03XU.f -- translated by f2c (version 20100827).
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

#line 1 "MB03XU.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b7 = 1.;
static doublereal c_b9 = 0.;

/* Subroutine */ int mb03xu_(logical *ltra, logical *ltrb, integer *n, 
	integer *k, integer *nb, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *g, integer *ldg, doublereal *q, integer *
	ldq, doublereal *xa, integer *ldxa, doublereal *xb, integer *ldxb, 
	doublereal *xg, integer *ldxg, doublereal *xq, integer *ldxq, 
	doublereal *ya, integer *ldya, doublereal *yb, integer *ldyb, 
	doublereal *yg, integer *ldyg, doublereal *yq, integer *ldyq, 
	doublereal *csl, doublereal *csr, doublereal *taul, doublereal *taur, 
	doublereal *dwork)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, g_dim1, g_offset, q_dim1, 
	    q_offset, xa_dim1, xa_offset, xb_dim1, xb_offset, xg_dim1, 
	    xg_offset, xq_dim1, xq_offset, ya_dim1, ya_offset, yb_dim1, 
	    yb_offset, yg_dim1, yg_offset, yq_dim1, yq_offset, i__1, i__2, 
	    i__3;
    doublereal d__1;

    /* Local variables */
    static doublereal c__;
    static integer i__, j;
    static doublereal s;
    static integer nb1, nb2, nb3, pdw;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp, tauq;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), dlarfg_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *), 
	    dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);


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

/*     To reduce 2*nb columns and rows of a real (k+2n)-by-(k+2n) */
/*     matrix H: */

/*             [ op(A)   G   ] */
/*         H = [             ], */
/*             [  Q    op(B) ] */

/*     so that elements in the first nb columns below the k-th */
/*     subdiagonal of the (k+n)-by-n matrix op(A), in the first nb */
/*     columns and rows of the n-by-n matrix Q and in the first nb rows */
/*     above the diagonal of the n-by-(k+n) matrix op(B) are zero. */
/*     The reduction is performed by orthogonal symplectic */
/*     transformations UU'*H*VV and matrices U, V, YA, YB, YG, YQ, XA, */
/*     XB, XG, and XQ are returned so that */

/*                    [ op(Aout)+U*YA'+XA*V'     G+U*YG'+XG*V'    ] */
/*         UU' H VV = [                                           ]. */
/*                    [   Qout+U*YQ'+XQ*V'   op(Bout)+U*YB'+XB*V' ] */

/*     This is an auxiliary routine called by MB04TB. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     LTRA    LOGICAL */
/*             Specifies the form of op( A ) as follows: */
/*             = .FALSE.:  op( A ) = A; */
/*             = .TRUE.:   op( A ) = A'. */

/*     LTRB    LOGICAL */
/*             Specifies the form of op( B ) as follows: */
/*             = .FALSE.:  op( B ) = B; */
/*             = .TRUE.:   op( B ) = B'. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix Q. N >= 0. */

/*     K       (input) INTEGER */
/*             The offset of the reduction. Elements below the K-th */
/*             subdiagonal in the first NB columns of op(A) are */
/*             reduced to zero. K >= 0. */

/*     NB      (input) INTEGER */
/*             The number of columns/rows to be reduced. N > NB >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDA,N)     if LTRA = .FALSE. */
/*                     (LDA,K+N)   if LTRA = .TRUE. */
/*             On entry with LTRA = .FALSE., the leading (K+N)-by-N part */
/*             of this array must contain the matrix A. */
/*             On entry with LTRA = .TRUE., the leading N-by-(K+N) part */
/*             of this array must contain the matrix A. */
/*             On exit with LTRA = .FALSE., the leading (K+N)-by-N part */
/*             of this array contains the matrix Aout and, in the zero */
/*             parts, information about the elementary reflectors used to */
/*             compute the reduction. */
/*             On exit with LTRA = .TRUE., the leading N-by-(K+N) part of */
/*             this array contains the matrix Aout and in the zero parts */
/*             information about the elementary reflectors. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. */
/*             LDA >= MAX(1,K+N),  if LTRA = .FALSE.; */
/*             LDA >= MAX(1,N),    if LTRA = .TRUE.. */

/*     B       (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDB,K+N)   if LTRB = .FALSE. */
/*                     (LDB,N)     if LTRB = .TRUE. */
/*             On entry with LTRB = .FALSE., the leading N-by-(K+N) part */
/*             of this array must contain the matrix B. */
/*             On entry with LTRB = .TRUE., the leading (K+N)-by-N part */
/*             of this array must contain the matrix B. */
/*             On exit with LTRB = .FALSE., the leading N-by-(K+N) part */
/*             of this array contains the matrix Bout and, in the zero */
/*             parts, information about the elementary reflectors used to */
/*             compute the reduction. */
/*             On exit with LTRB = .TRUE., the leading (K+N)-by-N part of */
/*             this array contains the matrix Bout and in the zero parts */
/*             information about the elementary reflectors. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= MAX(1,N),    if LTRB = .FALSE.; */
/*             LDB >= MAX(1,K+N),  if LTRB = .TRUE.. */

/*     G       (input/output) DOUBLE PRECISION array, dimension (LDG,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix G. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the matrix Gout. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G.  LDG >= MAX(1,N). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix Q. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the matrix Qout and in the zero parts information about */
/*             the elementary reflectors used to compute the reduction. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q.  LDQ >= MAX(1,N). */

/*     XA      (output) DOUBLE PRECISION array, dimension (LDXA,2*NB) */
/*             On exit, the leading N-by-(2*NB) part of this array */
/*             contains the matrix XA. */

/*     LDXA    INTEGER */
/*             The leading dimension of the array XA.  LDXA >= MAX(1,N). */

/*     XB      (output) DOUBLE PRECISION array, dimension (LDXB,2*NB) */
/*             On exit, the leading (K+N)-by-(2*NB) part of this array */
/*             contains the matrix XB. */

/*     LDXB    INTEGER */
/*             The leading dimension of the array XB. LDXB >= MAX(1,K+N). */

/*     XG      (output) DOUBLE PRECISION array, dimension (LDXG,2*NB) */
/*             On exit, the leading (K+N)-by-(2*NB) part of this array */
/*             contains the matrix XG. */

/*     LDXG    INTEGER */
/*             The leading dimension of the array XG. LDXG >= MAX(1,K+N). */

/*     XQ      (output) DOUBLE PRECISION array, dimension (LDXQ,2*NB) */
/*             On exit, the leading N-by-(2*NB) part of this array */
/*             contains the matrix XQ. */

/*     LDXQ    INTEGER */
/*             The leading dimension of the array XQ.  LDXQ >= MAX(1,N). */

/*     YA      (output) DOUBLE PRECISION array, dimension (LDYA,2*NB) */
/*             On exit, the leading (K+N)-by-(2*NB) part of this array */
/*             contains the matrix YA. */

/*     LDYA    INTEGER */
/*             The leading dimension of the array YA. LDYA >= MAX(1,K+N). */

/*     YB      (output) DOUBLE PRECISION array, dimension (LDYB,2*NB) */
/*             On exit, the leading N-by-(2*NB) part of this array */
/*             contains the matrix YB. */

/*     LDYB    INTEGER */
/*             The leading dimension of the array YB.  LDYB >= MAX(1,N). */

/*     YG      (output) DOUBLE PRECISION array, dimension (LDYG,2*NB) */
/*             On exit, the leading (K+N)-by-(2*NB) part of this array */
/*             contains the matrix YG. */

/*     LDYG    INTEGER */
/*             The leading dimension of the array YG. LDYG >= MAX(1,K+N). */

/*     YQ      (output) DOUBLE PRECISION array, dimension (LDYQ,2*NB) */
/*             On exit, the leading N-by-(2*NB) part of this array */
/*             contains the matrix YQ. */

/*     LDYQ    INTEGER */
/*             The leading dimension of the array YQ.  LDYQ >= MAX(1,N). */

/*     CSL     (output) DOUBLE PRECISION array, dimension (2*NB) */
/*             On exit, the first 2NB elements of this array contain the */
/*             cosines and sines of the symplectic Givens rotations */
/*             applied from the left-hand side used to compute the */
/*             reduction. */

/*     CSR     (output) DOUBLE PRECISION array, dimension (2*NB) */
/*             On exit, the first 2NB-2 elements of this array contain */
/*             the cosines and sines of the symplectic Givens rotations */
/*             applied from the right-hand side used to compute the */
/*             reduction. */

/*     TAUL    (output) DOUBLE PRECISION array, dimension (NB) */
/*             On exit, the first NB elements of this array contain the */
/*             scalar factors of some of the elementary reflectors */
/*             applied form the left-hand side. */

/*     TAUR    (output) DOUBLE PRECISION array, dimension (NB) */
/*             On exit, the first NB-1 elements of this array contain the */
/*             scalar factors of some of the elementary reflectors */
/*             applied form the right-hand side. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (5*NB) */

/*     METHOD */

/*     For details regarding the representation of the orthogonal */
/*     symplectic matrices UU and VV within the arrays A, B, CSL, CSR, Q, */
/*     TAUL and TAUR see the description of MB04TB. */

/*     The contents of A, B, G and Q on exit are illustrated by the */
/*     following example with op(A) = A, op(B) = B, n = 5, k = 2 and */
/*     nb = 2: */

/*          ( a  r  r  a  a  )       ( g  g  g  r  r  g  g  ) */
/*          ( a  r  r  a  a  )       ( g  g  g  r  r  g  g  ) */
/*          ( r  r  r  r  r  )       ( r  r  r  r  r  r  r  ) */
/*      A = ( u2 r  r  r  r  ),  G = ( r  r  r  r  r  r  r  ), */
/*          ( u2 u2 r  a  a  )       ( g  g  g  r  r  g  g  ) */
/*          ( u2 u2 r  a  a  )       ( g  g  g  r  r  g  g  ) */
/*          ( u2 u2 r  a  a  )       ( g  g  g  r  r  g  g  ) */

/*          ( t  t  v1 v1 v1 )       ( r  r  r  r  r  v2 v2 ) */
/*          ( u1 t  t  v1 v1 )       ( r  r  r  r  r  r  v2 ) */
/*      Q = ( u1 u1 r  q  q  ),  B = ( b  b  b  r  r  b  b  ). */
/*          ( u1 u1 r  q  q  )       ( b  b  b  r  r  b  b  ) */
/*          ( u1 u1 r  q  q  )       ( b  b  b  r  r  b  b  ) */

/*     where a, b, g and q denote elements of the original matrices, r */
/*     denotes a modified element, t denotes a scalar factor of an */
/*     applied elementary reflector, ui and vi denote elements of the */
/*     matrices U and V, respectively. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires ( 16*K + 32*N + 42 )*N*NB + */
/*     ( 16*K + 112*N - 208/3*NB - 69 )*NB*NB - 29/3*NB floating point */
/*     operations and is numerically backward stable. */

/*     REFERENCES */

/*     [1] Benner, P., Mehrmann, V., and Xu, H. */
/*         A numerically stable, structure preserving method for */
/*         computing the eigenvalues of real Hamiltonian or symplectic */
/*         pencils. */
/*         Numer. Math., Vol. 78 (3), pp. 329-358, 1998. */

/*     [2] Kressner, D. */
/*         Block algorithms for orthogonal symplectic factorizations. */
/*         BIT Numerical Mathematics, 43 (4), pp. 775-790, 2003. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLASUB). */

/*     KEYWORDS */

/*     Elementary matrix operations, Matrix decompositions, Hamiltonian */
/*     matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */

/*     .. Executable Statements .. */

/*     Quick return if possible. */

#line 303 "MB03XU.f"
    /* Parameter adjustments */
#line 303 "MB03XU.f"
    a_dim1 = *lda;
#line 303 "MB03XU.f"
    a_offset = 1 + a_dim1;
#line 303 "MB03XU.f"
    a -= a_offset;
#line 303 "MB03XU.f"
    b_dim1 = *ldb;
#line 303 "MB03XU.f"
    b_offset = 1 + b_dim1;
#line 303 "MB03XU.f"
    b -= b_offset;
#line 303 "MB03XU.f"
    g_dim1 = *ldg;
#line 303 "MB03XU.f"
    g_offset = 1 + g_dim1;
#line 303 "MB03XU.f"
    g -= g_offset;
#line 303 "MB03XU.f"
    q_dim1 = *ldq;
#line 303 "MB03XU.f"
    q_offset = 1 + q_dim1;
#line 303 "MB03XU.f"
    q -= q_offset;
#line 303 "MB03XU.f"
    xa_dim1 = *ldxa;
#line 303 "MB03XU.f"
    xa_offset = 1 + xa_dim1;
#line 303 "MB03XU.f"
    xa -= xa_offset;
#line 303 "MB03XU.f"
    xb_dim1 = *ldxb;
#line 303 "MB03XU.f"
    xb_offset = 1 + xb_dim1;
#line 303 "MB03XU.f"
    xb -= xb_offset;
#line 303 "MB03XU.f"
    xg_dim1 = *ldxg;
#line 303 "MB03XU.f"
    xg_offset = 1 + xg_dim1;
#line 303 "MB03XU.f"
    xg -= xg_offset;
#line 303 "MB03XU.f"
    xq_dim1 = *ldxq;
#line 303 "MB03XU.f"
    xq_offset = 1 + xq_dim1;
#line 303 "MB03XU.f"
    xq -= xq_offset;
#line 303 "MB03XU.f"
    ya_dim1 = *ldya;
#line 303 "MB03XU.f"
    ya_offset = 1 + ya_dim1;
#line 303 "MB03XU.f"
    ya -= ya_offset;
#line 303 "MB03XU.f"
    yb_dim1 = *ldyb;
#line 303 "MB03XU.f"
    yb_offset = 1 + yb_dim1;
#line 303 "MB03XU.f"
    yb -= yb_offset;
#line 303 "MB03XU.f"
    yg_dim1 = *ldyg;
#line 303 "MB03XU.f"
    yg_offset = 1 + yg_dim1;
#line 303 "MB03XU.f"
    yg -= yg_offset;
#line 303 "MB03XU.f"
    yq_dim1 = *ldyq;
#line 303 "MB03XU.f"
    yq_offset = 1 + yq_dim1;
#line 303 "MB03XU.f"
    yq -= yq_offset;
#line 303 "MB03XU.f"
    --csl;
#line 303 "MB03XU.f"
    --csr;
#line 303 "MB03XU.f"
    --taul;
#line 303 "MB03XU.f"
    --taur;
#line 303 "MB03XU.f"
    --dwork;
#line 303 "MB03XU.f"

#line 303 "MB03XU.f"
    /* Function Body */
#line 303 "MB03XU.f"
    if (*n + *k <= 0) {
#line 304 "MB03XU.f"
	return 0;
#line 305 "MB03XU.f"
    }

#line 307 "MB03XU.f"
    nb1 = *nb + 1;
#line 308 "MB03XU.f"
    nb2 = *nb + *nb;
#line 309 "MB03XU.f"
    nb3 = nb2 + *nb;
#line 310 "MB03XU.f"
    pdw = nb3 + *nb + 1;

#line 312 "MB03XU.f"
    if (*ltra && *ltrb) {
#line 313 "MB03XU.f"
	i__1 = *nb;
#line 313 "MB03XU.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Transform first row/column of A and Q. See routine MB04TS. */

#line 317 "MB03XU.f"
	    alpha = q[i__ + i__ * q_dim1];
#line 318 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 318 "MB03XU.f"
	    dlarfg_(&i__2, &alpha, &q[i__ + 1 + i__ * q_dim1], &c__1, &tauq);
#line 319 "MB03XU.f"
	    q[i__ + i__ * q_dim1] = 1.;
#line 320 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 320 "MB03XU.f"
	    temp = -tauq * ddot_(&i__2, &q[i__ + i__ * q_dim1], &c__1, &a[i__ 
		    + (*k + i__) * a_dim1], lda);
#line 321 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 321 "MB03XU.f"
	    daxpy_(&i__2, &temp, &q[i__ + i__ * q_dim1], &c__1, &a[i__ + (*k 
		    + i__) * a_dim1], lda);
#line 322 "MB03XU.f"
	    temp = a[i__ + (*k + i__) * a_dim1];
#line 323 "MB03XU.f"
	    dlartg_(&temp, &alpha, &c__, &s, &a[i__ + (*k + i__) * a_dim1]);
#line 324 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 324 "MB03XU.f"
	    dlarfg_(&i__2, &a[i__ + (*k + i__) * a_dim1], &a[i__ + (*k + i__ 
		    + 1) * a_dim1], lda, &taul[i__]);
#line 325 "MB03XU.f"
	    temp = a[i__ + (*k + i__) * a_dim1];
#line 326 "MB03XU.f"
	    a[i__ + (*k + i__) * a_dim1] = 1.;

/*           Update XQ with first Householder reflection. */

#line 330 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 330 "MB03XU.f"
	    i__3 = *n - i__;
#line 330 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[
		    i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)9);
#line 332 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 332 "MB03XU.f"
	    i__3 = i__ - 1;
#line 332 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + q_dim1], ldq, &
		    q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[1], &c__1, (
		    ftnlen)9);
#line 334 "MB03XU.f"
	    i__2 = *n - i__;
#line 334 "MB03XU.f"
	    i__3 = i__ - 1;
#line 334 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + xq_dim1]
		    , ldxq, &dwork[1], &c__1, &c_b7, &xq[i__ + 1 + i__ * 
		    xq_dim1], &c__1, (ftnlen)12);
#line 336 "MB03XU.f"
	    i__2 = i__ - 1;
#line 336 "MB03XU.f"
	    i__3 = *n - i__ + 1;
#line 336 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[(*k + i__) * 
		    a_dim1 + 1], lda, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &
		    dwork[*nb + 1], &c__1, (ftnlen)12);
#line 338 "MB03XU.f"
	    i__2 = *n - i__;
#line 338 "MB03XU.f"
	    i__3 = i__ - 1;
#line 338 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[*nb + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + i__ * xq_dim1], &c__1, (ftnlen)12);
#line 340 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 340 "MB03XU.f"
	    i__3 = i__ - 1;
#line 340 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + yq_dim1], ldyq,
		     &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[i__ * xq_dim1 
		    + 1], &c__1, (ftnlen)9);
#line 342 "MB03XU.f"
	    i__2 = i__ - 1;
#line 342 "MB03XU.f"
	    i__3 = *n - i__;
#line 342 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xq[i__ * xq_dim1 + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + i__ * xq_dim1], &c__1, (ftnlen)9);
#line 344 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 344 "MB03XU.f"
	    i__3 = i__ - 1;
#line 344 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + nb1 * yq_dim1],
		     ldyq, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[(i__ + *
		    nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
#line 346 "MB03XU.f"
	    i__2 = *n - i__;
#line 346 "MB03XU.f"
	    i__3 = i__ - 1;
#line 346 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &
		    c_b7, &xq[i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)12);
#line 348 "MB03XU.f"
	    i__2 = *n - i__;
#line 348 "MB03XU.f"
	    d__1 = -tauq;
#line 348 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + i__ * xq_dim1], &c__1);

/*           Update Q(i,i+1:n). */

#line 352 "MB03XU.f"
	    i__2 = *n - i__;
#line 352 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1],
		     ldxq, &q[i__ + q_dim1], ldq, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, (ftnlen)12);
#line 354 "MB03XU.f"
	    i__2 = *n - i__;
#line 354 "MB03XU.f"
	    i__3 = i__ - 1;
#line 354 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &a[(*k + i__) * a_dim1 + 1], &c__1, &c_b7,
		     &q[i__ + (i__ + 1) * q_dim1], ldq, (ftnlen)12);
#line 356 "MB03XU.f"
	    i__2 = i__ - 1;
#line 356 "MB03XU.f"
	    i__3 = *n - i__;
#line 356 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yq[i__ + yq_dim1], ldyq, &c_b7, &q[i__ + (i__ + 
		    1) * q_dim1], ldq, (ftnlen)9);
#line 358 "MB03XU.f"
	    i__2 = *n - i__;
#line 358 "MB03XU.f"
	    i__3 = i__ - 1;
#line 358 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &yq[i__ + nb1 * yq_dim1], ldyq, &c_b7, &q[
		    i__ + (i__ + 1) * q_dim1], ldq, (ftnlen)12);

/*           Update XA with first Householder reflection. */

#line 363 "MB03XU.f"
	    i__2 = *n - i__;
#line 363 "MB03XU.f"
	    i__3 = *n - i__ + 1;
#line 363 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[i__ + 1 + (*k + 
		    i__) * a_dim1], lda, &q[i__ + i__ * q_dim1], &c__1, &c_b9,
		     &xa[i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)12);
#line 365 "MB03XU.f"
	    i__2 = *n - i__;
#line 365 "MB03XU.f"
	    i__3 = i__ - 1;
#line 365 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + xa_dim1]
		    , ldxa, &dwork[1], &c__1, &c_b7, &xa[i__ + 1 + i__ * 
		    xa_dim1], &c__1, (ftnlen)12);
#line 367 "MB03XU.f"
	    i__2 = *n - i__;
#line 367 "MB03XU.f"
	    i__3 = i__ - 1;
#line 367 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[*nb + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + i__ * xa_dim1], &c__1, (ftnlen)12);
#line 369 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 369 "MB03XU.f"
	    i__3 = i__ - 1;
#line 369 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + ya_dim1], 
		    ldya, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[i__ * 
		    xa_dim1 + 1], &c__1, (ftnlen)9);
#line 371 "MB03XU.f"
	    i__2 = i__ - 1;
#line 371 "MB03XU.f"
	    i__3 = *n - i__;
#line 371 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xa[i__ * xa_dim1 + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + i__ * xa_dim1], &c__1, (ftnlen)9);
#line 373 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 373 "MB03XU.f"
	    i__3 = i__ - 1;
#line 373 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + nb1 * 
		    ya_dim1], ldya, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[
		    i__ * xa_dim1 + 1], &c__1, (ftnlen)9);
#line 375 "MB03XU.f"
	    i__2 = *n - i__;
#line 375 "MB03XU.f"
	    i__3 = i__ - 1;
#line 375 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &xa[i__ * xa_dim1 + 1], &c__1, &c_b7, &xa[
		    i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)12);
#line 377 "MB03XU.f"
	    i__2 = *n - i__;
#line 377 "MB03XU.f"
	    d__1 = -tauq;
#line 377 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + i__ * xa_dim1], &c__1);

/*           Update A(i+1:n,k+i). */

#line 381 "MB03XU.f"
	    i__2 = *n - i__;
#line 381 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1],
		     ldxa, &q[i__ + q_dim1], ldq, &c_b7, &a[i__ + 1 + (*k + 
		    i__) * a_dim1], &c__1, (ftnlen)12);
#line 383 "MB03XU.f"
	    i__2 = *n - i__;
#line 383 "MB03XU.f"
	    i__3 = i__ - 1;
#line 383 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &a[(*k + i__) * a_dim1 + 1], &c__1, &c_b7,
		     &a[i__ + 1 + (*k + i__) * a_dim1], &c__1, (ftnlen)12);
#line 385 "MB03XU.f"
	    i__2 = i__ - 1;
#line 385 "MB03XU.f"
	    i__3 = *n - i__;
#line 385 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &ya[*k + i__ + ya_dim1], ldya, &c_b7, &a[i__ + 1 
		    + (*k + i__) * a_dim1], &c__1, (ftnlen)9);
#line 387 "MB03XU.f"
	    i__2 = *n - i__;
#line 387 "MB03XU.f"
	    i__3 = i__ - 1;
#line 387 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &ya[*k + i__ + nb1 * ya_dim1], ldya, &c_b7, 
		    &a[i__ + 1 + (*k + i__) * a_dim1], &c__1, (ftnlen)12);

/*           Apply rotation to [ A(i+1:n,k+i)'; Q(i,i+1:n) ]. */

#line 392 "MB03XU.f"
	    i__2 = *n - i__;
#line 392 "MB03XU.f"
	    drot_(&i__2, &a[i__ + 1 + (*k + i__) * a_dim1], &c__1, &q[i__ + (
		    i__ + 1) * q_dim1], ldq, &c__, &s);

/*           Update XQ with second Householder reflection. */

#line 396 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 396 "MB03XU.f"
	    i__3 = *n - i__;
#line 396 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, &a[i__ + (*k + i__) * a_dim1], lda, &c_b9, &
		    xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)9);
#line 398 "MB03XU.f"
	    i__2 = *n - i__;
#line 398 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], ldq,
		     &a[i__ + (*k + i__ + 1) * a_dim1], lda, &c_b9, &dwork[
		    nb2 + 1], &c__1, (ftnlen)9);
#line 400 "MB03XU.f"
	    i__2 = *n - i__;
#line 400 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1],
		     ldxq, &dwork[nb2 + 1], &c__1, &c_b7, &xq[i__ + 1 + (i__ 
		    + *nb) * xq_dim1], &c__1, (ftnlen)12);
#line 402 "MB03XU.f"
	    i__2 = i__ - 1;
#line 402 "MB03XU.f"
	    i__3 = *n - i__;
#line 402 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &a[i__ + (*k + i__ + 1) * a_dim1], lda, 
		    &c_b9, &dwork[nb3 + 1], &c__1, (ftnlen)12);
#line 404 "MB03XU.f"
	    i__2 = *n - i__;
#line 404 "MB03XU.f"
	    i__3 = i__ - 1;
#line 404 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[nb3 + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)12);
#line 406 "MB03XU.f"
	    i__2 = *n - i__;
#line 406 "MB03XU.f"
	    i__3 = i__ - 1;
#line 406 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + yq_dim1], 
		    ldyq, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &c_b9, &xq[(
		    i__ + *nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
#line 408 "MB03XU.f"
	    i__2 = i__ - 1;
#line 408 "MB03XU.f"
	    i__3 = *n - i__;
#line 408 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &c_b7, &
		    xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)9);
#line 410 "MB03XU.f"
	    i__2 = *n - i__;
#line 410 "MB03XU.f"
	    i__3 = i__ - 1;
#line 410 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
#line 412 "MB03XU.f"
	    i__2 = *n - i__;
#line 412 "MB03XU.f"
	    i__3 = i__ - 1;
#line 412 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &
		    c_b7, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (
		    ftnlen)12);
#line 414 "MB03XU.f"
	    i__2 = *n - i__;
#line 414 "MB03XU.f"
	    d__1 = -taul[i__];
#line 414 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1);

/*           Update Q(i,i+1:n). */

#line 418 "MB03XU.f"
	    i__2 = *n - i__;
#line 418 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, 
		    &q[i__ + (i__ + 1) * q_dim1], ldq);

/*           Update XA with second Householder reflection. */

#line 422 "MB03XU.f"
	    i__2 = *n - i__;
#line 422 "MB03XU.f"
	    i__3 = *n - i__ + 1;
#line 422 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[i__ + 1 + (*k + 
		    i__) * a_dim1], lda, &a[i__ + (*k + i__) * a_dim1], lda, &
		    c_b9, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (
		    ftnlen)12);
#line 424 "MB03XU.f"
	    i__2 = *n - i__;
#line 424 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1],
		     ldxa, &dwork[nb2 + 1], &c__1, &c_b7, &xa[i__ + 1 + (i__ 
		    + *nb) * xa_dim1], &c__1, (ftnlen)12);
#line 426 "MB03XU.f"
	    i__2 = *n - i__;
#line 426 "MB03XU.f"
	    i__3 = i__ - 1;
#line 426 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[nb3 + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)12);
#line 428 "MB03XU.f"
	    i__2 = *n - i__;
#line 428 "MB03XU.f"
	    i__3 = i__ - 1;
#line 428 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + 1 + 
		    ya_dim1], ldya, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, (ftnlen)9);
#line 430 "MB03XU.f"
	    i__2 = i__ - 1;
#line 430 "MB03XU.f"
	    i__3 = *n - i__;
#line 430 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, &c_b7, &
		    xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)9);
#line 432 "MB03XU.f"
	    i__2 = *n - i__;
#line 432 "MB03XU.f"
	    i__3 = i__ - 1;
#line 432 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + 1 + nb1 * 
		    ya_dim1], ldya, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, (ftnlen)9);
#line 434 "MB03XU.f"
	    i__2 = *n - i__;
#line 434 "MB03XU.f"
	    i__3 = i__ - 1;
#line 434 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, &
		    c_b7, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (
		    ftnlen)12);
#line 436 "MB03XU.f"
	    i__2 = *n - i__;
#line 436 "MB03XU.f"
	    d__1 = -taul[i__];
#line 436 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1);

/*           Update A(i+1:n,k+i). */

#line 440 "MB03XU.f"
	    i__2 = *n - i__;
#line 440 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, 
		    &a[i__ + 1 + (*k + i__) * a_dim1], &c__1);

/*           Update XG with first Householder reflection. */

#line 444 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 444 "MB03XU.f"
	    i__3 = *k + *n;
#line 444 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &g[*k + i__ + g_dim1], 
		    ldg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xg[i__ * 
		    xg_dim1 + 1], &c__1, (ftnlen)9);
#line 446 "MB03XU.f"
	    i__2 = *k + *n;
#line 446 "MB03XU.f"
	    i__3 = i__ - 1;
#line 446 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[xg_offset], ldxg, 
		    &dwork[1], &c__1, &c_b7, &xg[i__ * xg_dim1 + 1], &c__1, (
		    ftnlen)12);
#line 448 "MB03XU.f"
	    i__2 = *k + *n;
#line 448 "MB03XU.f"
	    i__3 = i__ - 1;
#line 448 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &dwork[*nb + 1], &c__1, &c_b7, &xg[i__ * xg_dim1 
		    + 1], &c__1, (ftnlen)12);
#line 450 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 450 "MB03XU.f"
	    i__3 = i__ - 1;
#line 450 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + yg_dim1], 
		    ldyg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], &
		    c__1, (ftnlen)9);
#line 452 "MB03XU.f"
	    i__2 = i__ - 1;
#line 452 "MB03XU.f"
	    i__3 = *n - i__;
#line 452 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 + 
		    i__ * xg_dim1], &c__1, (ftnlen)9);
#line 454 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 454 "MB03XU.f"
	    i__3 = i__ - 1;
#line 454 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + nb1 * 
		    yg_dim1], ldyg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
#line 456 "MB03XU.f"
	    i__2 = *n - i__;
#line 456 "MB03XU.f"
	    i__3 = i__ - 1;
#line 456 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 
		    + i__ * xg_dim1], &c__1, (ftnlen)12);
#line 458 "MB03XU.f"
	    i__2 = *k + *n;
#line 458 "MB03XU.f"
	    d__1 = -tauq;
#line 458 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xg[i__ * xg_dim1 + 1], &c__1);

/*           Update G(k+i,:). */

#line 462 "MB03XU.f"
	    i__2 = *k + *n;
#line 462 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xg[xg_offset], ldxg, &
		    q[i__ + q_dim1], ldq, &c_b7, &g[*k + i__ + g_dim1], ldg, (
		    ftnlen)12);
#line 464 "MB03XU.f"
	    i__2 = *k + *n;
#line 464 "MB03XU.f"
	    i__3 = i__ - 1;
#line 464 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &a[(*k + i__) * a_dim1 + 1], &c__1, &c_b7, &g[*k 
		    + i__ + g_dim1], ldg, (ftnlen)12);
#line 466 "MB03XU.f"
	    i__2 = i__ - 1;
#line 466 "MB03XU.f"
	    i__3 = *n - i__;
#line 466 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yg[*k + i__ + yg_dim1], ldyg, &c_b7, &g[*k + 
		    i__ + (*k + i__ + 1) * g_dim1], ldg, (ftnlen)9);
#line 468 "MB03XU.f"
	    i__2 = *n - i__;
#line 468 "MB03XU.f"
	    i__3 = i__ - 1;
#line 468 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &yg[*k + i__ + nb1 * yg_dim1], ldyg, &c_b7, 
		    &g[*k + i__ + (*k + i__ + 1) * g_dim1], ldg, (ftnlen)12);

/*           Update XB with first Householder reflection. */

#line 473 "MB03XU.f"
	    i__2 = *k + *n;
#line 473 "MB03XU.f"
	    i__3 = *n - i__ + 1;
#line 473 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[i__ * b_dim1 + 1], 
		    ldb, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xb[i__ * 
		    xb_dim1 + 1], &c__1, (ftnlen)12);
#line 475 "MB03XU.f"
	    i__2 = *k + *n;
#line 475 "MB03XU.f"
	    i__3 = i__ - 1;
#line 475 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[xb_offset], ldxb, 
		    &dwork[1], &c__1, &c_b7, &xb[i__ * xb_dim1 + 1], &c__1, (
		    ftnlen)12);
#line 477 "MB03XU.f"
	    i__2 = *k + *n;
#line 477 "MB03XU.f"
	    i__3 = i__ - 1;
#line 477 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &dwork[*nb + 1], &c__1, &c_b7, &xb[i__ * xb_dim1 
		    + 1], &c__1, (ftnlen)12);
#line 479 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 479 "MB03XU.f"
	    i__3 = i__ - 1;
#line 479 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + yb_dim1], ldyb,
		     &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], &c__1,
		     (ftnlen)9);
#line 481 "MB03XU.f"
	    i__2 = i__ - 1;
#line 481 "MB03XU.f"
	    i__3 = *n - i__;
#line 481 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 + 
		    i__ * xb_dim1], &c__1, (ftnlen)9);
#line 483 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 483 "MB03XU.f"
	    i__3 = i__ - 1;
#line 483 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + nb1 * yb_dim1],
		     ldyb, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], 
		    &c__1, (ftnlen)9);
#line 485 "MB03XU.f"
	    i__2 = *n - i__;
#line 485 "MB03XU.f"
	    i__3 = i__ - 1;
#line 485 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 
		    + i__ * xb_dim1], &c__1, (ftnlen)12);
#line 487 "MB03XU.f"
	    i__2 = *k + *n;
#line 487 "MB03XU.f"
	    d__1 = -tauq;
#line 487 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xb[i__ * xb_dim1 + 1], &c__1);

/*           Update B(:,i). */

#line 491 "MB03XU.f"
	    i__2 = *k + *n;
#line 491 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xb[xb_offset], ldxb, &
		    q[i__ + q_dim1], ldq, &c_b7, &b[i__ * b_dim1 + 1], &c__1, 
		    (ftnlen)12);
#line 493 "MB03XU.f"
	    i__2 = *k + *n;
#line 493 "MB03XU.f"
	    i__3 = i__ - 1;
#line 493 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &a[(*k + i__) * a_dim1 + 1], &c__1, &c_b7, &b[i__ 
		    * b_dim1 + 1], &c__1, (ftnlen)12);
#line 495 "MB03XU.f"
	    i__2 = i__ - 1;
#line 495 "MB03XU.f"
	    i__3 = *n - i__;
#line 495 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yb[i__ + yb_dim1], ldyb, &c_b7, &b[*k + i__ + 1 
		    + i__ * b_dim1], &c__1, (ftnlen)9);
#line 497 "MB03XU.f"
	    i__2 = *n - i__;
#line 497 "MB03XU.f"
	    i__3 = i__ - 1;
#line 497 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &yb[i__ + nb1 * yb_dim1], ldyb, &c_b7, &b[*
		    k + i__ + 1 + i__ * b_dim1], &c__1, (ftnlen)12);

/*           Apply rotation to [ G(k+i,:); B(:,i)' ]. */

#line 502 "MB03XU.f"
	    i__2 = *k + *n;
#line 502 "MB03XU.f"
	    drot_(&i__2, &g[*k + i__ + g_dim1], ldg, &b[i__ * b_dim1 + 1], &
		    c__1, &c__, &s);

#line 504 "MB03XU.f"
	    i__2 = i__ - 1;
#line 504 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 505 "MB03XU.f"
		yg[*k + i__ + j * yg_dim1] = 0.;
#line 506 "MB03XU.f"
/* L10: */
#line 506 "MB03XU.f"
	    }
#line 507 "MB03XU.f"
	    i__2 = i__ - 1;
#line 507 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 508 "MB03XU.f"
		yg[*k + i__ + (*nb + j) * yg_dim1] = 0.;
#line 509 "MB03XU.f"
/* L20: */
#line 509 "MB03XU.f"
	    }
#line 510 "MB03XU.f"
	    i__2 = i__ - 1;
#line 510 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 511 "MB03XU.f"
		ya[*k + i__ + j * ya_dim1] = 0.;
#line 512 "MB03XU.f"
/* L30: */
#line 512 "MB03XU.f"
	    }
#line 513 "MB03XU.f"
	    i__2 = i__ - 1;
#line 513 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 514 "MB03XU.f"
		ya[*k + i__ + (*nb + j) * ya_dim1] = 0.;
#line 515 "MB03XU.f"
/* L40: */
#line 515 "MB03XU.f"
	    }

/*           Update XG with second Householder reflection. */

#line 519 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 519 "MB03XU.f"
	    i__3 = *k + *n;
#line 519 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &g[*k + i__ + g_dim1], 
		    ldg, &a[i__ + (*k + i__) * a_dim1], lda, &c_b9, &xg[(i__ 
		    + *nb) * xg_dim1 + 1], &c__1, (ftnlen)9);
#line 521 "MB03XU.f"
	    i__2 = *k + *n;
#line 521 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xg[xg_offset], ldxg, &
		    dwork[nb2 + 1], &c__1, &c_b7, &xg[(i__ + *nb) * xg_dim1 + 
		    1], &c__1, (ftnlen)12);
#line 523 "MB03XU.f"
	    i__2 = *k + *n;
#line 523 "MB03XU.f"
	    i__3 = i__ - 1;
#line 523 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &dwork[nb3 + 1], &c__1, &c_b7, &xg[(i__ + *nb) * 
		    xg_dim1 + 1], &c__1, (ftnlen)12);
#line 525 "MB03XU.f"
	    i__2 = *n - i__;
#line 525 "MB03XU.f"
	    i__3 = i__ - 1;
#line 525 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + 1 + 
		    yg_dim1], ldyg, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 527 "MB03XU.f"
	    i__2 = i__ - 1;
#line 527 "MB03XU.f"
	    i__3 = *n - i__;
#line 527 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 + (
		    i__ + *nb) * xg_dim1], &c__1, (ftnlen)9);
#line 529 "MB03XU.f"
	    i__2 = *n - i__;
#line 529 "MB03XU.f"
	    i__3 = i__ - 1;
#line 529 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + 1 + nb1 * 
		    yg_dim1], ldyg, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 531 "MB03XU.f"
	    i__2 = *n - i__;
#line 531 "MB03XU.f"
	    i__3 = i__ - 1;
#line 531 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 
		    + (i__ + *nb) * xg_dim1], &c__1, (ftnlen)12);
#line 533 "MB03XU.f"
	    i__2 = *k + *n;
#line 533 "MB03XU.f"
	    d__1 = -taul[i__];
#line 533 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xg[(i__ + *nb) * xg_dim1 + 1], &c__1);

/*           Update G(k+i,:). */

#line 537 "MB03XU.f"
	    i__2 = *k + *n;
#line 537 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &xg[(i__ + *nb) * xg_dim1 + 1], &c__1, &g[*k 
		    + i__ + g_dim1], ldg);

/*           Update XB with second Householder reflection. */

#line 541 "MB03XU.f"
	    i__2 = *k + *n;
#line 541 "MB03XU.f"
	    i__3 = *n - i__ + 1;
#line 541 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[i__ * b_dim1 + 1], 
		    ldb, &a[i__ + (*k + i__) * a_dim1], lda, &c_b9, &xb[(i__ 
		    + *nb) * xb_dim1 + 1], &c__1, (ftnlen)12);
#line 543 "MB03XU.f"
	    i__2 = *k + *n;
#line 543 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xb[xb_offset], ldxb, &
		    dwork[nb2 + 1], &c__1, &c_b7, &xb[(i__ + *nb) * xb_dim1 + 
		    1], &c__1, (ftnlen)12);
#line 545 "MB03XU.f"
	    i__2 = *k + *n;
#line 545 "MB03XU.f"
	    i__3 = i__ - 1;
#line 545 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &dwork[nb3 + 1], &c__1, &c_b7, &xb[(i__ + *nb) * 
		    xb_dim1 + 1], &c__1, (ftnlen)12);
#line 547 "MB03XU.f"
	    i__2 = *n - i__;
#line 547 "MB03XU.f"
	    i__3 = i__ - 1;
#line 547 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + yb_dim1], 
		    ldyb, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
#line 549 "MB03XU.f"
	    i__2 = i__ - 1;
#line 549 "MB03XU.f"
	    i__3 = *n - i__;
#line 549 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 + (
		    i__ + *nb) * xb_dim1], &c__1, (ftnlen)9);
#line 551 "MB03XU.f"
	    i__2 = *n - i__;
#line 551 "MB03XU.f"
	    i__3 = i__ - 1;
#line 551 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 553 "MB03XU.f"
	    i__2 = *n - i__;
#line 553 "MB03XU.f"
	    i__3 = i__ - 1;
#line 553 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 
		    + (i__ + *nb) * xb_dim1], &c__1, (ftnlen)12);
#line 555 "MB03XU.f"
	    i__2 = *k + *n;
#line 555 "MB03XU.f"
	    d__1 = -taul[i__];
#line 555 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xb[(i__ + *nb) * xb_dim1 + 1], &c__1);

/*           Update B(:,i). */

#line 559 "MB03XU.f"
	    i__2 = *k + *n;
#line 559 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &xb[(i__ + *nb) * xb_dim1 + 1], &c__1, &b[
		    i__ * b_dim1 + 1], &c__1);

#line 561 "MB03XU.f"
	    a[i__ + (*k + i__) * a_dim1] = temp;
#line 562 "MB03XU.f"
	    q[i__ + i__ * q_dim1] = tauq;
#line 563 "MB03XU.f"
	    csl[(i__ << 1) - 1] = c__;
#line 564 "MB03XU.f"
	    csl[i__ * 2] = s;

/*           Transform first row/column of Q and B. */

#line 568 "MB03XU.f"
	    alpha = q[i__ + (i__ + 1) * q_dim1];
#line 569 "MB03XU.f"
	    i__2 = *n - i__;
#line 569 "MB03XU.f"
	    dlarfg_(&i__2, &alpha, &q[i__ + (i__ + 2) * q_dim1], ldq, &tauq);
#line 570 "MB03XU.f"
	    q[i__ + (i__ + 1) * q_dim1] = 1.;
#line 571 "MB03XU.f"
	    i__2 = *n - i__;
#line 571 "MB03XU.f"
	    temp = -tauq * ddot_(&i__2, &q[i__ + (i__ + 1) * q_dim1], ldq, &b[
		    *k + i__ + 1 + i__ * b_dim1], &c__1);
#line 572 "MB03XU.f"
	    i__2 = *n - i__;
#line 572 "MB03XU.f"
	    daxpy_(&i__2, &temp, &q[i__ + (i__ + 1) * q_dim1], ldq, &b[*k + 
		    i__ + 1 + i__ * b_dim1], &c__1);
#line 573 "MB03XU.f"
	    temp = b[*k + i__ + 1 + i__ * b_dim1];
#line 574 "MB03XU.f"
	    dlartg_(&temp, &alpha, &c__, &s, &b[*k + i__ + 1 + i__ * b_dim1]);
#line 575 "MB03XU.f"
	    s = -s;
#line 576 "MB03XU.f"
	    i__2 = *n - i__;
#line 576 "MB03XU.f"
	    dlarfg_(&i__2, &b[*k + i__ + 1 + i__ * b_dim1], &b[*k + i__ + 2 + 
		    i__ * b_dim1], &c__1, &taur[i__]);
#line 577 "MB03XU.f"
	    temp = b[*k + i__ + 1 + i__ * b_dim1];
#line 578 "MB03XU.f"
	    b[*k + i__ + 1 + i__ * b_dim1] = 1.;

/*           Update YB with first Householder reflection. */

#line 582 "MB03XU.f"
	    i__2 = *n - i__;
#line 582 "MB03XU.f"
	    i__3 = *n - i__;
#line 582 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + (i__ + 
		    1) * b_dim1], ldb, &q[i__ + (i__ + 1) * q_dim1], ldq, &
		    c_b9, &yb[i__ + 1 + i__ * yb_dim1], &c__1, (ftnlen)9);
#line 584 "MB03XU.f"
	    i__2 = *n - i__;
#line 584 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 1 + 
		    xb_dim1], ldxb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yb[i__ * yb_dim1 + 1], &c__1, (ftnlen)9);
#line 586 "MB03XU.f"
	    i__2 = *n - i__;
#line 586 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yb[i__ * yb_dim1 + 1], &c__1, &c_b7, &yb[i__ + 1 + 
		    i__ * yb_dim1], &c__1, (ftnlen)12);
#line 588 "MB03XU.f"
	    i__2 = *n - i__;
#line 588 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 1 + nb1 * 
		    xb_dim1], ldxb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yb[i__ * yb_dim1 + 1], &c__1, (ftnlen)9);
#line 590 "MB03XU.f"
	    i__2 = *n - i__;
#line 590 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &yb[i__ * yb_dim1 + 1], &c__1, &c_b7, &
		    yb[i__ + 1 + i__ * yb_dim1], &c__1, (ftnlen)9);
#line 592 "MB03XU.f"
	    i__2 = i__ - 1;
#line 592 "MB03XU.f"
	    i__3 = *n - i__;
#line 592 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 
		    + 1], ldq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &
		    dwork[1], &c__1, (ftnlen)12);
#line 594 "MB03XU.f"
	    i__2 = *n - i__;
#line 594 "MB03XU.f"
	    i__3 = i__ - 1;
#line 594 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + yb_dim1]
		    , ldyb, &dwork[1], &c__1, &c_b7, &yb[i__ + 1 + i__ * 
		    yb_dim1], &c__1, (ftnlen)12);
#line 596 "MB03XU.f"
	    i__2 = *n - i__;
#line 596 "MB03XU.f"
	    i__3 = i__ - 1;
#line 596 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + b_dim1]
		    , ldb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &dwork[*
		    nb + 1], &c__1, (ftnlen)9);
#line 598 "MB03XU.f"
	    i__2 = *n - i__;
#line 598 "MB03XU.f"
	    i__3 = i__ - 1;
#line 598 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &dwork[*nb + 1], &c__1, &c_b7, &yb[i__ + 
		    1 + i__ * yb_dim1], &c__1, (ftnlen)12);
#line 600 "MB03XU.f"
	    i__2 = *n - i__;
#line 600 "MB03XU.f"
	    d__1 = -tauq;
#line 600 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yb[i__ + 1 + i__ * yb_dim1], &c__1);

/*           Update B(k+i+1,i+1:n). */

#line 604 "MB03XU.f"
	    i__2 = *n - i__;
#line 604 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xb[*k + i__ + 1 + xb_dim1], ldxb, &c_b7, &b[*k + 
		    i__ + 1 + (i__ + 1) * b_dim1], ldb, (ftnlen)12);
#line 606 "MB03XU.f"
	    i__2 = *n - i__;
#line 606 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &xb[*k + i__ + 1 + nb1 * xb_dim1], ldxb,
		     &c_b7, &b[*k + i__ + 1 + (i__ + 1) * b_dim1], ldb, (
		    ftnlen)9);
#line 608 "MB03XU.f"
	    i__2 = *n - i__;
#line 608 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yb[i__ + 1 + yb_dim1],
		     ldyb, &q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &b[*k + 
		    i__ + 1 + (i__ + 1) * b_dim1], ldb, (ftnlen)12);
#line 610 "MB03XU.f"
	    i__2 = *n - i__;
#line 610 "MB03XU.f"
	    i__3 = i__ - 1;
#line 610 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &b[*k + i__ + 1 + b_dim1], ldb, &c_b7, &b[
		    *k + i__ + 1 + (i__ + 1) * b_dim1], ldb, (ftnlen)12);

/*           Update YQ with first Householder reflection. */

#line 615 "MB03XU.f"
	    i__2 = *n - i__;
#line 615 "MB03XU.f"
	    i__3 = *n - i__;
#line 615 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[i__ + 1 + (i__ + 1)
		     * q_dim1], ldq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9,
		     &yq[i__ + 1 + i__ * yq_dim1], &c__1, (ftnlen)12);
#line 617 "MB03XU.f"
	    i__2 = *n - i__;
#line 617 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1], 
		    ldxq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &yq[i__ * 
		    yq_dim1 + 1], &c__1, (ftnlen)9);
#line 619 "MB03XU.f"
	    i__2 = *n - i__;
#line 619 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yq[i__ * yq_dim1 + 1], &c__1, &c_b7, &yq[i__ + 1 + 
		    i__ * yq_dim1], &c__1, (ftnlen)12);
#line 621 "MB03XU.f"
	    i__2 = *n - i__;
#line 621 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yq[i__ * yq_dim1 + 1], &c__1, (ftnlen)9);
#line 623 "MB03XU.f"
	    i__2 = *n - i__;
#line 623 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &yq[i__ * yq_dim1 + 1], &c__1, &c_b7, &
		    yq[i__ + 1 + i__ * yq_dim1], &c__1, (ftnlen)9);
#line 625 "MB03XU.f"
	    i__2 = *n - i__;
#line 625 "MB03XU.f"
	    i__3 = i__ - 1;
#line 625 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + yq_dim1]
		    , ldyq, &dwork[1], &c__1, &c_b7, &yq[i__ + 1 + i__ * 
		    yq_dim1], &c__1, (ftnlen)12);
#line 627 "MB03XU.f"
	    i__2 = *n - i__;
#line 627 "MB03XU.f"
	    i__3 = i__ - 1;
#line 627 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &dwork[*nb + 1], &c__1, &c_b7, &yq[i__ + 
		    1 + i__ * yq_dim1], &c__1, (ftnlen)12);
#line 629 "MB03XU.f"
	    i__2 = *n - i__;
#line 629 "MB03XU.f"
	    d__1 = -tauq;
#line 629 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yq[i__ + 1 + i__ * yq_dim1], &c__1);

/*           Update Q(i+1:n,i+1). */

#line 633 "MB03XU.f"
	    i__2 = *n - i__;
#line 633 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xq[i__ + 1 + xq_dim1], ldxq, &c_b7, &q[i__ + 1 + (
		    i__ + 1) * q_dim1], &c__1, (ftnlen)12);
#line 635 "MB03XU.f"
	    i__2 = *n - i__;
#line 635 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &xq[i__ + 1 + nb1 * xq_dim1], ldxq, &
		    c_b7, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)9);
#line 637 "MB03XU.f"
	    i__2 = *n - i__;
#line 637 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yq[i__ + 1 + yq_dim1],
		     ldyq, &q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &q[i__ + 
		    1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12);
#line 639 "MB03XU.f"
	    i__2 = *n - i__;
#line 639 "MB03XU.f"
	    i__3 = i__ - 1;
#line 639 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &b[*k + i__ + 1 + b_dim1], ldb, &c_b7, &q[
		    i__ + 1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12);

/*           Apply rotation to [ Q(i+1:n,i+1), B(k+i+1,i+1:n)' ]. */

#line 644 "MB03XU.f"
	    i__2 = *n - i__;
#line 644 "MB03XU.f"
	    drot_(&i__2, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, &b[*k + i__ 
		    + 1 + (i__ + 1) * b_dim1], ldb, &c__, &s);
#line 645 "MB03XU.f"
	    i__2 = i__;
#line 645 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 646 "MB03XU.f"
		xb[*k + i__ + 1 + j * xb_dim1] = 0.;
#line 647 "MB03XU.f"
/* L50: */
#line 647 "MB03XU.f"
	    }
#line 648 "MB03XU.f"
	    i__2 = i__;
#line 648 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 649 "MB03XU.f"
		xb[*k + i__ + 1 + (*nb + j) * xb_dim1] = 0.;
#line 650 "MB03XU.f"
/* L60: */
#line 650 "MB03XU.f"
	    }

/*           Update YB with second Householder reflection. */

#line 654 "MB03XU.f"
	    i__2 = *n - i__;
#line 654 "MB03XU.f"
	    i__3 = *n - i__;
#line 654 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + (i__ + 
		    1) * b_dim1], ldb, &b[*k + i__ + 1 + i__ * b_dim1], &c__1,
		     &c_b9, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (
		    ftnlen)9);
#line 656 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 656 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 2 + 
		    xb_dim1], ldxb, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, (ftnlen)9);
#line 658 "MB03XU.f"
	    i__2 = *n - i__;
#line 658 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, &c_b7, &yb[
		    i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (ftnlen)12);
#line 660 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 660 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 2 + nb1 * 
		    xb_dim1], ldxb, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, (ftnlen)9);
#line 662 "MB03XU.f"
	    i__2 = *n - i__;
#line 662 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, &
		    c_b7, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (
		    ftnlen)9);
#line 664 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 664 "MB03XU.f"
	    dgemv_("No transpose", &i__, &i__2, &c_b7, &q[(i__ + 2) * q_dim1 
		    + 1], ldq, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &c_b9, 
		    &dwork[nb2 + 1], &c__1, (ftnlen)12);
#line 666 "MB03XU.f"
	    i__2 = *n - i__;
#line 666 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yb[i__ + 1 + yb_dim1],
		     ldyb, &dwork[nb2 + 1], &c__1, &c_b7, &yb[i__ + 1 + (i__ 
		    + *nb) * yb_dim1], &c__1, (ftnlen)12);
#line 668 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 668 "MB03XU.f"
	    i__3 = i__ - 1;
#line 668 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 2 + b_dim1]
		    , ldq, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &c_b9, &
		    dwork[nb3 + 1], &c__1, (ftnlen)9);
#line 670 "MB03XU.f"
	    i__2 = *n - i__;
#line 670 "MB03XU.f"
	    i__3 = i__ - 1;
#line 670 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &dwork[nb3 + 1], &c__1, &c_b7, &yb[i__ + 
		    1 + (i__ + *nb) * yb_dim1], &c__1, (ftnlen)12);
#line 672 "MB03XU.f"
	    i__2 = *n - i__;
#line 672 "MB03XU.f"
	    d__1 = -taur[i__];
#line 672 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1);

/*           Update B(k+i+1,i+1:n). */

#line 676 "MB03XU.f"
	    i__2 = *n - i__;
#line 676 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, 
		    &b[*k + i__ + 1 + (i__ + 1) * b_dim1], ldb);

/*           Update YQ with second Householder reflection. */

#line 680 "MB03XU.f"
	    i__2 = *n - i__;
#line 680 "MB03XU.f"
	    i__3 = *n - i__;
#line 680 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[i__ + 1 + (i__ + 1)
		     * q_dim1], ldq, &b[*k + i__ + 1 + i__ * b_dim1], &c__1, &
		    c_b9, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (
		    ftnlen)12);
#line 682 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 682 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 2 + xq_dim1], 
		    ldxq, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &c_b9, &yq[(
		    i__ + *nb) * yq_dim1 + 1], &c__1, (ftnlen)9);
#line 684 "MB03XU.f"
	    i__2 = *n - i__;
#line 684 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, &c_b7, &yq[
		    i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (ftnlen)12);
#line 686 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 686 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 2 + nb1 * 
		    xq_dim1], ldxq, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, (ftnlen)9);
#line 688 "MB03XU.f"
	    i__2 = *n - i__;
#line 688 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, &
		    c_b7, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (
		    ftnlen)9);
#line 690 "MB03XU.f"
	    i__2 = *n - i__;
#line 690 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yq[i__ + 1 + yq_dim1],
		     ldyq, &dwork[nb2 + 1], &c__1, &c_b7, &yq[i__ + 1 + (i__ 
		    + *nb) * yq_dim1], &c__1, (ftnlen)12);
#line 692 "MB03XU.f"
	    i__2 = *n - i__;
#line 692 "MB03XU.f"
	    i__3 = i__ - 1;
#line 692 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &dwork[nb3 + 1], &c__1, &c_b7, &yq[i__ + 
		    1 + (i__ + *nb) * yq_dim1], &c__1, (ftnlen)12);
#line 694 "MB03XU.f"
	    i__2 = *n - i__;
#line 694 "MB03XU.f"
	    d__1 = -taur[i__];
#line 694 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1);

/*           Update Q(i+1:n,i+1). */

#line 698 "MB03XU.f"
	    i__2 = *n - i__;
#line 698 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, 
		    &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1);

/*           Update YA with first Householder reflection. */

#line 702 "MB03XU.f"
	    i__2 = *n - i__;
#line 702 "MB03XU.f"
	    i__3 = *k + *n;
#line 702 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[i__ + 1 + a_dim1], 
		    lda, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &ya[i__ * 
		    ya_dim1 + 1], &c__1, (ftnlen)9);
#line 704 "MB03XU.f"
	    i__2 = *n - i__;
#line 704 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1], 
		    ldxa, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &dwork[
		    pdw], &c__1, (ftnlen)9);
#line 706 "MB03XU.f"
	    i__2 = *n - i__;
#line 706 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 + i__ * 
		    ya_dim1], &c__1, (ftnlen)12);
#line 708 "MB03XU.f"
	    i__2 = *n - i__;
#line 708 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
#line 710 "MB03XU.f"
	    i__2 = *n - i__;
#line 710 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ 
		    + 1 + i__ * ya_dim1], &c__1, (ftnlen)9);
#line 712 "MB03XU.f"
	    i__2 = *k + *n;
#line 712 "MB03XU.f"
	    i__3 = i__ - 1;
#line 712 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[ya_offset], ldya, 
		    &dwork[1], &c__1, &c_b7, &ya[i__ * ya_dim1 + 1], &c__1, (
		    ftnlen)12);
#line 714 "MB03XU.f"
	    i__2 = *k + *n;
#line 714 "MB03XU.f"
	    i__3 = i__ - 1;
#line 714 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &dwork[*nb + 1], &c__1, &c_b7, &ya[i__ * ya_dim1 
		    + 1], &c__1, (ftnlen)12);
#line 716 "MB03XU.f"
	    i__2 = *k + *n;
#line 716 "MB03XU.f"
	    d__1 = -tauq;
#line 716 "MB03XU.f"
	    dscal_(&i__2, &d__1, &ya[i__ * ya_dim1 + 1], &c__1);

/*           Update A(i+1,1:k+n). */

#line 720 "MB03XU.f"
	    i__2 = *n - i__;
#line 720 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xa[i__ + 1 + xa_dim1], ldxa, &c_b7, &a[i__ + 1 + (*
		    k + i__ + 1) * a_dim1], lda, (ftnlen)12);
#line 722 "MB03XU.f"
	    i__2 = *n - i__;
#line 722 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &xa[i__ + 1 + nb1 * xa_dim1], ldxa, &
		    c_b7, &a[i__ + 1 + (*k + i__ + 1) * a_dim1], lda, (ftnlen)
		    9);
#line 724 "MB03XU.f"
	    i__2 = *k + *n;
#line 724 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &ya[ya_offset], ldya, &
		    q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &a[i__ + 1 + 
		    a_dim1], lda, (ftnlen)12);
#line 726 "MB03XU.f"
	    i__2 = *k + *n;
#line 726 "MB03XU.f"
	    i__3 = i__ - 1;
#line 726 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &b[*k + i__ + 1 + b_dim1], ldb, &c_b7, &a[i__ + 1 
		    + a_dim1], lda, (ftnlen)12);

/*           Update YG with first Householder reflection. */

#line 731 "MB03XU.f"
	    i__2 = *k + *n;
#line 731 "MB03XU.f"
	    i__3 = *n - i__;
#line 731 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], ldg, &q[i__ + (i__ + 1) * q_dim1], ldq, &
		    c_b9, &yg[i__ * yg_dim1 + 1], &c__1, (ftnlen)12);
#line 733 "MB03XU.f"
	    i__2 = *n - i__;
#line 733 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 1 + 
		    xg_dim1], ldxg, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
#line 735 "MB03XU.f"
	    i__2 = *n - i__;
#line 735 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 + i__ * 
		    yg_dim1], &c__1, (ftnlen)12);
#line 737 "MB03XU.f"
	    i__2 = *n - i__;
#line 737 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 1 + nb1 * 
		    xg_dim1], ldxg, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
#line 739 "MB03XU.f"
	    i__2 = *n - i__;
#line 739 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ 
		    + 1 + i__ * yg_dim1], &c__1, (ftnlen)9);
#line 741 "MB03XU.f"
	    i__2 = *k + *n;
#line 741 "MB03XU.f"
	    i__3 = i__ - 1;
#line 741 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[yg_offset], ldyg, 
		    &dwork[1], &c__1, &c_b7, &yg[i__ * yg_dim1 + 1], &c__1, (
		    ftnlen)12);
#line 743 "MB03XU.f"
	    i__2 = *k + *n;
#line 743 "MB03XU.f"
	    i__3 = i__ - 1;
#line 743 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &dwork[*nb + 1], &c__1, &c_b7, &yg[i__ * yg_dim1 
		    + 1], &c__1, (ftnlen)12);
#line 745 "MB03XU.f"
	    i__2 = *k + *n;
#line 745 "MB03XU.f"
	    d__1 = -tauq;
#line 745 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yg[i__ * yg_dim1 + 1], &c__1);

/*           Update G(1:k+n,k+i+1). */

#line 749 "MB03XU.f"
	    i__2 = *n - i__;
#line 749 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xg[*k + i__ + 1 + xg_dim1], ldxg, &c_b7, &g[*k + 
		    i__ + 1 + (*k + i__ + 1) * g_dim1], &c__1, (ftnlen)12);
#line 751 "MB03XU.f"
	    i__2 = *n - i__;
#line 751 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &xg[*k + i__ + 1 + nb1 * xg_dim1], ldxg,
		     &c_b7, &g[*k + i__ + 1 + (*k + i__ + 1) * g_dim1], &c__1,
		     (ftnlen)9);
#line 753 "MB03XU.f"
	    i__2 = *k + *n;
#line 753 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yg[yg_offset], ldyg, &
		    q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &g[(*k + i__ + 1)
		     * g_dim1 + 1], &c__1, (ftnlen)12);
#line 755 "MB03XU.f"
	    i__2 = *k + *n;
#line 755 "MB03XU.f"
	    i__3 = i__ - 1;
#line 755 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &b[*k + i__ + 1 + b_dim1], ldb, &c_b7, &g[(*k + 
		    i__ + 1) * g_dim1 + 1], &c__1, (ftnlen)12);
#line 757 "MB03XU.f"
	    i__2 = i__;
#line 757 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 758 "MB03XU.f"
		xg[*k + i__ + 1 + j * xg_dim1] = 0.;
#line 759 "MB03XU.f"
/* L70: */
#line 759 "MB03XU.f"
	    }
#line 760 "MB03XU.f"
	    i__2 = i__;
#line 760 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 761 "MB03XU.f"
		xg[*k + i__ + 1 + (*nb + j) * xg_dim1] = 0.;
#line 762 "MB03XU.f"
/* L80: */
#line 762 "MB03XU.f"
	    }

/*           Apply rotation to [ A(i+1,1:k+n)', G(1:k+n,k+i+1) ]. */

#line 766 "MB03XU.f"
	    i__2 = *k + *n;
#line 766 "MB03XU.f"
	    drot_(&i__2, &a[i__ + 1 + a_dim1], lda, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], &c__1, &c__, &s);

/*           Update YA with second Householder reflection. */

#line 770 "MB03XU.f"
	    i__2 = *n - i__;
#line 770 "MB03XU.f"
	    i__3 = *k + *n;
#line 770 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[i__ + 1 + a_dim1], 
		    lda, &b[*k + i__ + 1 + i__ * b_dim1], &c__1, &c_b9, &ya[(
		    i__ + *nb) * ya_dim1 + 1], &c__1, (ftnlen)9);
#line 772 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 772 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 2 + xa_dim1], 
		    ldxa, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
#line 774 "MB03XU.f"
	    i__2 = *n - i__;
#line 774 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 + (i__ + 
		    *nb) * ya_dim1], &c__1, (ftnlen)12);
#line 776 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 776 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 2 + nb1 * 
		    xa_dim1], ldxa, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 778 "MB03XU.f"
	    i__2 = *n - i__;
#line 778 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ 
		    + 1 + (i__ + *nb) * ya_dim1], &c__1, (ftnlen)9);
#line 780 "MB03XU.f"
	    i__2 = *k + *n;
#line 780 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &ya[ya_offset], ldya, &
		    dwork[nb2 + 1], &c__1, &c_b7, &ya[(i__ + *nb) * ya_dim1 + 
		    1], &c__1, (ftnlen)12);
#line 782 "MB03XU.f"
	    i__2 = *k + *n;
#line 782 "MB03XU.f"
	    i__3 = i__ - 1;
#line 782 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &dwork[nb3 + 1], &c__1, &c_b7, &ya[(i__ + *nb) * 
		    ya_dim1 + 1], &c__1, (ftnlen)12);
#line 784 "MB03XU.f"
	    i__2 = *k + *n;
#line 784 "MB03XU.f"
	    d__1 = -taur[i__];
#line 784 "MB03XU.f"
	    dscal_(&i__2, &d__1, &ya[(i__ + *nb) * ya_dim1 + 1], &c__1);

/*           Update A(i+1,1:k+n). */

#line 788 "MB03XU.f"
	    i__2 = *k + *n;
#line 788 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &ya[(i__ + *nb) * ya_dim1 + 1], &c__1, &a[
		    i__ + 1 + a_dim1], lda);

/*           Update YG with second Householder reflection. */

#line 792 "MB03XU.f"
	    i__2 = *k + *n;
#line 792 "MB03XU.f"
	    i__3 = *n - i__;
#line 792 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], ldg, &b[*k + i__ + 1 + i__ * b_dim1], &c__1, 
		    &c_b9, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1, (ftnlen)12);
#line 794 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 794 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 2 + 
		    xg_dim1], ldxg, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 796 "MB03XU.f"
	    i__2 = *n - i__;
#line 796 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 + (i__ + 
		    *nb) * yg_dim1], &c__1, (ftnlen)12);
#line 798 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 798 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 2 + nb1 * 
		    xg_dim1], ldxg, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 800 "MB03XU.f"
	    i__2 = *n - i__;
#line 800 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ 
		    + 1 + (i__ + *nb) * yg_dim1], &c__1, (ftnlen)9);
#line 802 "MB03XU.f"
	    i__2 = *k + *n;
#line 802 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yg[yg_offset], ldyg, &
		    dwork[nb2 + 1], &c__1, &c_b7, &yg[(i__ + *nb) * yg_dim1 + 
		    1], &c__1, (ftnlen)12);
#line 804 "MB03XU.f"
	    i__2 = *k + *n;
#line 804 "MB03XU.f"
	    i__3 = i__ - 1;
#line 804 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &dwork[nb3 + 1], &c__1, &c_b7, &yg[(i__ + *nb) * 
		    yg_dim1 + 1], &c__1, (ftnlen)12);
#line 806 "MB03XU.f"
	    i__2 = *k + *n;
#line 806 "MB03XU.f"
	    d__1 = -taur[i__];
#line 806 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1);

/*           Update G(1:k+n,k+i+1). */

#line 810 "MB03XU.f"
	    i__2 = *k + *n;
#line 810 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1, &g[(*
		    k + i__ + 1) * g_dim1 + 1], &c__1);

#line 812 "MB03XU.f"
	    b[*k + i__ + 1 + i__ * b_dim1] = temp;
#line 813 "MB03XU.f"
	    q[i__ + (i__ + 1) * q_dim1] = tauq;
#line 814 "MB03XU.f"
	    csr[(i__ << 1) - 1] = c__;
#line 815 "MB03XU.f"
	    csr[i__ * 2] = s;
#line 816 "MB03XU.f"
/* L90: */
#line 816 "MB03XU.f"
	}
#line 817 "MB03XU.f"
    } else if (*ltra) {
#line 818 "MB03XU.f"
	i__1 = *nb;
#line 818 "MB03XU.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Transform first row/column of A and Q. See routine MB04TS. */

#line 822 "MB03XU.f"
	    alpha = q[i__ + i__ * q_dim1];
#line 823 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 823 "MB03XU.f"
	    dlarfg_(&i__2, &alpha, &q[i__ + 1 + i__ * q_dim1], &c__1, &tauq);
#line 824 "MB03XU.f"
	    q[i__ + i__ * q_dim1] = 1.;
#line 825 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 825 "MB03XU.f"
	    temp = -tauq * ddot_(&i__2, &q[i__ + i__ * q_dim1], &c__1, &a[i__ 
		    + (*k + i__) * a_dim1], lda);
#line 826 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 826 "MB03XU.f"
	    daxpy_(&i__2, &temp, &q[i__ + i__ * q_dim1], &c__1, &a[i__ + (*k 
		    + i__) * a_dim1], lda);
#line 827 "MB03XU.f"
	    temp = a[i__ + (*k + i__) * a_dim1];
#line 828 "MB03XU.f"
	    dlartg_(&temp, &alpha, &c__, &s, &a[i__ + (*k + i__) * a_dim1]);
#line 829 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 829 "MB03XU.f"
	    dlarfg_(&i__2, &a[i__ + (*k + i__) * a_dim1], &a[i__ + (*k + i__ 
		    + 1) * a_dim1], lda, &taul[i__]);
#line 830 "MB03XU.f"
	    temp = a[i__ + (*k + i__) * a_dim1];
#line 831 "MB03XU.f"
	    a[i__ + (*k + i__) * a_dim1] = 1.;

/*           Update XQ with first Householder reflection. */

#line 835 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 835 "MB03XU.f"
	    i__3 = *n - i__;
#line 835 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[
		    i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)9);
#line 837 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 837 "MB03XU.f"
	    i__3 = i__ - 1;
#line 837 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + q_dim1], ldq, &
		    q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[1], &c__1, (
		    ftnlen)9);
#line 839 "MB03XU.f"
	    i__2 = *n - i__;
#line 839 "MB03XU.f"
	    i__3 = i__ - 1;
#line 839 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + xq_dim1]
		    , ldxq, &dwork[1], &c__1, &c_b7, &xq[i__ + 1 + i__ * 
		    xq_dim1], &c__1, (ftnlen)12);
#line 841 "MB03XU.f"
	    i__2 = i__ - 1;
#line 841 "MB03XU.f"
	    i__3 = *n - i__ + 1;
#line 841 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[(*k + i__) * 
		    a_dim1 + 1], lda, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &
		    dwork[*nb + 1], &c__1, (ftnlen)12);
#line 843 "MB03XU.f"
	    i__2 = *n - i__;
#line 843 "MB03XU.f"
	    i__3 = i__ - 1;
#line 843 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[*nb + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + i__ * xq_dim1], &c__1, (ftnlen)12);
#line 845 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 845 "MB03XU.f"
	    i__3 = i__ - 1;
#line 845 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + yq_dim1], ldyq,
		     &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[i__ * xq_dim1 
		    + 1], &c__1, (ftnlen)9);
#line 847 "MB03XU.f"
	    i__2 = i__ - 1;
#line 847 "MB03XU.f"
	    i__3 = *n - i__;
#line 847 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xq[i__ * xq_dim1 + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + i__ * xq_dim1], &c__1, (ftnlen)9);
#line 849 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 849 "MB03XU.f"
	    i__3 = i__ - 1;
#line 849 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + nb1 * yq_dim1],
		     ldyq, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[(i__ + *
		    nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
#line 851 "MB03XU.f"
	    i__2 = i__ - 1;
#line 851 "MB03XU.f"
	    i__3 = *n - i__;
#line 851 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &
		    c_b7, &xq[i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)9);
#line 853 "MB03XU.f"
	    i__2 = *n - i__;
#line 853 "MB03XU.f"
	    d__1 = -tauq;
#line 853 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + i__ * xq_dim1], &c__1);

/*           Update Q(i,i+1:n). */

#line 857 "MB03XU.f"
	    i__2 = *n - i__;
#line 857 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1],
		     ldxq, &q[i__ + q_dim1], ldq, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, (ftnlen)12);
#line 859 "MB03XU.f"
	    i__2 = *n - i__;
#line 859 "MB03XU.f"
	    i__3 = i__ - 1;
#line 859 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &a[(*k + i__) * a_dim1 + 1], &c__1, &c_b7,
		     &q[i__ + (i__ + 1) * q_dim1], ldq, (ftnlen)12);
#line 861 "MB03XU.f"
	    i__2 = i__ - 1;
#line 861 "MB03XU.f"
	    i__3 = *n - i__;
#line 861 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yq[i__ + yq_dim1], ldyq, &c_b7, &q[i__ + (i__ + 
		    1) * q_dim1], ldq, (ftnlen)9);
#line 863 "MB03XU.f"
	    i__2 = i__ - 1;
#line 863 "MB03XU.f"
	    i__3 = *n - i__;
#line 863 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &yq[i__ + nb1 * yq_dim1], ldyq, &c_b7, &
		    q[i__ + (i__ + 1) * q_dim1], ldq, (ftnlen)9);

/*           Update XA with first Householder reflection. */

#line 868 "MB03XU.f"
	    i__2 = *n - i__;
#line 868 "MB03XU.f"
	    i__3 = *n - i__ + 1;
#line 868 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[i__ + 1 + (*k + 
		    i__) * a_dim1], lda, &q[i__ + i__ * q_dim1], &c__1, &c_b9,
		     &xa[i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)12);
#line 870 "MB03XU.f"
	    i__2 = *n - i__;
#line 870 "MB03XU.f"
	    i__3 = i__ - 1;
#line 870 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + xa_dim1]
		    , ldxa, &dwork[1], &c__1, &c_b7, &xa[i__ + 1 + i__ * 
		    xa_dim1], &c__1, (ftnlen)12);
#line 872 "MB03XU.f"
	    i__2 = *n - i__;
#line 872 "MB03XU.f"
	    i__3 = i__ - 1;
#line 872 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[*nb + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + i__ * xa_dim1], &c__1, (ftnlen)12);
#line 874 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 874 "MB03XU.f"
	    i__3 = i__ - 1;
#line 874 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + ya_dim1], 
		    ldya, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[i__ * 
		    xa_dim1 + 1], &c__1, (ftnlen)9);
#line 876 "MB03XU.f"
	    i__2 = i__ - 1;
#line 876 "MB03XU.f"
	    i__3 = *n - i__;
#line 876 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xa[i__ * xa_dim1 + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + i__ * xa_dim1], &c__1, (ftnlen)9);
#line 878 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 878 "MB03XU.f"
	    i__3 = i__ - 1;
#line 878 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + nb1 * 
		    ya_dim1], ldya, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[
		    i__ * xa_dim1 + 1], &c__1, (ftnlen)9);
#line 880 "MB03XU.f"
	    i__2 = i__ - 1;
#line 880 "MB03XU.f"
	    i__3 = *n - i__;
#line 880 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &xa[i__ * xa_dim1 + 1], &c__1, &c_b7, &
		    xa[i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)9);
#line 882 "MB03XU.f"
	    i__2 = *n - i__;
#line 882 "MB03XU.f"
	    d__1 = -tauq;
#line 882 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + i__ * xa_dim1], &c__1);

/*           Update A(i+1:n,k+i). */

#line 886 "MB03XU.f"
	    i__2 = *n - i__;
#line 886 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1],
		     ldxa, &q[i__ + q_dim1], ldq, &c_b7, &a[i__ + 1 + (*k + 
		    i__) * a_dim1], &c__1, (ftnlen)12);
#line 888 "MB03XU.f"
	    i__2 = *n - i__;
#line 888 "MB03XU.f"
	    i__3 = i__ - 1;
#line 888 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &a[(*k + i__) * a_dim1 + 1], &c__1, &c_b7,
		     &a[i__ + 1 + (*k + i__) * a_dim1], &c__1, (ftnlen)12);
#line 890 "MB03XU.f"
	    i__2 = i__ - 1;
#line 890 "MB03XU.f"
	    i__3 = *n - i__;
#line 890 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &ya[*k + i__ + ya_dim1], ldya, &c_b7, &a[i__ + 1 
		    + (*k + i__) * a_dim1], &c__1, (ftnlen)9);
#line 892 "MB03XU.f"
	    i__2 = i__ - 1;
#line 892 "MB03XU.f"
	    i__3 = *n - i__;
#line 892 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &ya[*k + i__ + nb1 * ya_dim1], ldya, &
		    c_b7, &a[i__ + 1 + (*k + i__) * a_dim1], &c__1, (ftnlen)9)
		    ;

/*           Apply rotation to [ A(i+1:n,k+i)'; Q(i,i+1:n) ]. */

#line 897 "MB03XU.f"
	    i__2 = *n - i__;
#line 897 "MB03XU.f"
	    drot_(&i__2, &a[i__ + 1 + (*k + i__) * a_dim1], &c__1, &q[i__ + (
		    i__ + 1) * q_dim1], ldq, &c__, &s);

/*           Update XQ with second Householder reflection. */

#line 901 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 901 "MB03XU.f"
	    i__3 = *n - i__;
#line 901 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, &a[i__ + (*k + i__) * a_dim1], lda, &c_b9, &
		    xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)9);
#line 903 "MB03XU.f"
	    i__2 = *n - i__;
#line 903 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], ldq,
		     &a[i__ + (*k + i__ + 1) * a_dim1], lda, &c_b9, &dwork[
		    nb2 + 1], &c__1, (ftnlen)9);
#line 905 "MB03XU.f"
	    i__2 = *n - i__;
#line 905 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1],
		     ldxq, &dwork[nb2 + 1], &c__1, &c_b7, &xq[i__ + 1 + (i__ 
		    + *nb) * xq_dim1], &c__1, (ftnlen)12);
#line 907 "MB03XU.f"
	    i__2 = i__ - 1;
#line 907 "MB03XU.f"
	    i__3 = *n - i__;
#line 907 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &a[i__ + (*k + i__ + 1) * a_dim1], lda, 
		    &c_b9, &dwork[nb3 + 1], &c__1, (ftnlen)12);
#line 909 "MB03XU.f"
	    i__2 = *n - i__;
#line 909 "MB03XU.f"
	    i__3 = i__ - 1;
#line 909 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[nb3 + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)12);
#line 911 "MB03XU.f"
	    i__2 = *n - i__;
#line 911 "MB03XU.f"
	    i__3 = i__ - 1;
#line 911 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + yq_dim1], 
		    ldyq, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &c_b9, &xq[(
		    i__ + *nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
#line 913 "MB03XU.f"
	    i__2 = i__ - 1;
#line 913 "MB03XU.f"
	    i__3 = *n - i__;
#line 913 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &c_b7, &
		    xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)9);
#line 915 "MB03XU.f"
	    i__2 = *n - i__;
#line 915 "MB03XU.f"
	    i__3 = i__ - 1;
#line 915 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
#line 917 "MB03XU.f"
	    i__2 = i__ - 1;
#line 917 "MB03XU.f"
	    i__3 = *n - i__;
#line 917 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &
		    c_b7, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (
		    ftnlen)9);
#line 919 "MB03XU.f"
	    i__2 = *n - i__;
#line 919 "MB03XU.f"
	    d__1 = -taul[i__];
#line 919 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1);

/*           Update Q(i,i+1:n). */

#line 923 "MB03XU.f"
	    i__2 = *n - i__;
#line 923 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, 
		    &q[i__ + (i__ + 1) * q_dim1], ldq);

/*           Update XA with second Householder reflection. */

#line 927 "MB03XU.f"
	    i__2 = *n - i__;
#line 927 "MB03XU.f"
	    i__3 = *n - i__ + 1;
#line 927 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[i__ + 1 + (*k + 
		    i__) * a_dim1], lda, &a[i__ + (*k + i__) * a_dim1], lda, &
		    c_b9, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (
		    ftnlen)12);
#line 929 "MB03XU.f"
	    i__2 = *n - i__;
#line 929 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1],
		     ldxa, &dwork[nb2 + 1], &c__1, &c_b7, &xa[i__ + 1 + (i__ 
		    + *nb) * xa_dim1], &c__1, (ftnlen)12);
#line 931 "MB03XU.f"
	    i__2 = *n - i__;
#line 931 "MB03XU.f"
	    i__3 = i__ - 1;
#line 931 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[nb3 + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)12);
#line 933 "MB03XU.f"
	    i__2 = *n - i__;
#line 933 "MB03XU.f"
	    i__3 = i__ - 1;
#line 933 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + 1 + 
		    ya_dim1], ldya, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, (ftnlen)9);
#line 935 "MB03XU.f"
	    i__2 = i__ - 1;
#line 935 "MB03XU.f"
	    i__3 = *n - i__;
#line 935 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, &c_b7, &
		    xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)9);
#line 937 "MB03XU.f"
	    i__2 = *n - i__;
#line 937 "MB03XU.f"
	    i__3 = i__ - 1;
#line 937 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + 1 + nb1 * 
		    ya_dim1], ldya, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, (ftnlen)9);
#line 939 "MB03XU.f"
	    i__2 = i__ - 1;
#line 939 "MB03XU.f"
	    i__3 = *n - i__;
#line 939 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, &
		    c_b7, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (
		    ftnlen)9);
#line 941 "MB03XU.f"
	    i__2 = *n - i__;
#line 941 "MB03XU.f"
	    d__1 = -taul[i__];
#line 941 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1);

/*           Update A(i+1:n,k+i). */

#line 945 "MB03XU.f"
	    i__2 = *n - i__;
#line 945 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, 
		    &a[i__ + 1 + (*k + i__) * a_dim1], &c__1);

/*           Update XG with first Householder reflection. */

#line 949 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 949 "MB03XU.f"
	    i__3 = *k + *n;
#line 949 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &g[*k + i__ + g_dim1], 
		    ldg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xg[i__ * 
		    xg_dim1 + 1], &c__1, (ftnlen)9);
#line 951 "MB03XU.f"
	    i__2 = *k + *n;
#line 951 "MB03XU.f"
	    i__3 = i__ - 1;
#line 951 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[xg_offset], ldxg, 
		    &dwork[1], &c__1, &c_b7, &xg[i__ * xg_dim1 + 1], &c__1, (
		    ftnlen)12);
#line 953 "MB03XU.f"
	    i__2 = *k + *n;
#line 953 "MB03XU.f"
	    i__3 = i__ - 1;
#line 953 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &dwork[*nb + 1], &c__1, &c_b7, &xg[i__ * xg_dim1 
		    + 1], &c__1, (ftnlen)12);
#line 955 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 955 "MB03XU.f"
	    i__3 = i__ - 1;
#line 955 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + yg_dim1], 
		    ldyg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], &
		    c__1, (ftnlen)9);
#line 957 "MB03XU.f"
	    i__2 = i__ - 1;
#line 957 "MB03XU.f"
	    i__3 = *n - i__;
#line 957 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 + 
		    i__ * xg_dim1], &c__1, (ftnlen)9);
#line 959 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 959 "MB03XU.f"
	    i__3 = i__ - 1;
#line 959 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + nb1 * 
		    yg_dim1], ldyg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
#line 961 "MB03XU.f"
	    i__2 = i__ - 1;
#line 961 "MB03XU.f"
	    i__3 = *n - i__;
#line 961 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ 
		    + 1 + i__ * xg_dim1], &c__1, (ftnlen)9);
#line 963 "MB03XU.f"
	    i__2 = *k + *n;
#line 963 "MB03XU.f"
	    d__1 = -tauq;
#line 963 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xg[i__ * xg_dim1 + 1], &c__1);

/*           Update G(k+i,:). */

#line 967 "MB03XU.f"
	    i__2 = *k + *n;
#line 967 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xg[xg_offset], ldxg, &
		    q[i__ + q_dim1], ldq, &c_b7, &g[*k + i__ + g_dim1], ldg, (
		    ftnlen)12);
#line 969 "MB03XU.f"
	    i__2 = *k + *n;
#line 969 "MB03XU.f"
	    i__3 = i__ - 1;
#line 969 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &a[(*k + i__) * a_dim1 + 1], &c__1, &c_b7, &g[*k 
		    + i__ + g_dim1], ldg, (ftnlen)12);
#line 971 "MB03XU.f"
	    i__2 = i__ - 1;
#line 971 "MB03XU.f"
	    i__3 = *n - i__;
#line 971 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yg[*k + i__ + yg_dim1], ldyg, &c_b7, &g[*k + 
		    i__ + (*k + i__ + 1) * g_dim1], ldg, (ftnlen)9);
#line 973 "MB03XU.f"
	    i__2 = i__ - 1;
#line 973 "MB03XU.f"
	    i__3 = *n - i__;
#line 973 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &yg[*k + i__ + nb1 * yg_dim1], ldyg, &
		    c_b7, &g[*k + i__ + (*k + i__ + 1) * g_dim1], ldg, (
		    ftnlen)9);

/*           Update XB with first Householder reflection. */

#line 978 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 978 "MB03XU.f"
	    i__3 = *k + *n;
#line 978 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[i__ + b_dim1], ldb, &
		    q[i__ + i__ * q_dim1], &c__1, &c_b9, &xb[i__ * xb_dim1 + 
		    1], &c__1, (ftnlen)9);
#line 980 "MB03XU.f"
	    i__2 = *k + *n;
#line 980 "MB03XU.f"
	    i__3 = i__ - 1;
#line 980 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[xb_offset], ldxb, 
		    &dwork[1], &c__1, &c_b7, &xb[i__ * xb_dim1 + 1], &c__1, (
		    ftnlen)12);
#line 982 "MB03XU.f"
	    i__2 = *k + *n;
#line 982 "MB03XU.f"
	    i__3 = i__ - 1;
#line 982 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &dwork[*nb + 1], &c__1, &c_b7, &xb[i__ * xb_dim1 
		    + 1], &c__1, (ftnlen)12);
#line 984 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 984 "MB03XU.f"
	    i__3 = i__ - 1;
#line 984 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + yb_dim1], ldyb,
		     &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], &c__1,
		     (ftnlen)9);
#line 986 "MB03XU.f"
	    i__2 = i__ - 1;
#line 986 "MB03XU.f"
	    i__3 = *n - i__;
#line 986 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 + 
		    i__ * xb_dim1], &c__1, (ftnlen)9);
#line 988 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 988 "MB03XU.f"
	    i__3 = i__ - 1;
#line 988 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + nb1 * yb_dim1],
		     ldyb, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], 
		    &c__1, (ftnlen)9);
#line 990 "MB03XU.f"
	    i__2 = i__ - 1;
#line 990 "MB03XU.f"
	    i__3 = *n - i__;
#line 990 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ 
		    + 1 + i__ * xb_dim1], &c__1, (ftnlen)9);
#line 992 "MB03XU.f"
	    i__2 = *k + *n;
#line 992 "MB03XU.f"
	    d__1 = -tauq;
#line 992 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xb[i__ * xb_dim1 + 1], &c__1);

/*           Update B(i,:). */

#line 996 "MB03XU.f"
	    i__2 = *k + *n;
#line 996 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xb[xb_offset], ldxb, &
		    q[i__ + q_dim1], ldq, &c_b7, &b[i__ + b_dim1], ldb, (
		    ftnlen)12);
#line 998 "MB03XU.f"
	    i__2 = *k + *n;
#line 998 "MB03XU.f"
	    i__3 = i__ - 1;
#line 998 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &a[(*k + i__) * a_dim1 + 1], &c__1, &c_b7, &b[i__ 
		    + b_dim1], ldb, (ftnlen)12);
#line 1000 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1000 "MB03XU.f"
	    i__3 = *n - i__;
#line 1000 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yb[i__ + yb_dim1], ldyb, &c_b7, &b[i__ + (*k + 
		    i__ + 1) * b_dim1], ldb, (ftnlen)9);
#line 1002 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1002 "MB03XU.f"
	    i__3 = *n - i__;
#line 1002 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &yb[i__ + nb1 * yb_dim1], ldyb, &c_b7, &
		    b[i__ + (*k + i__ + 1) * b_dim1], ldb, (ftnlen)9);

/*           Apply rotation to [ G(k+i,:); B(i,:) ]. */

#line 1007 "MB03XU.f"
	    i__2 = *k + *n;
#line 1007 "MB03XU.f"
	    drot_(&i__2, &g[*k + i__ + g_dim1], ldg, &b[i__ + b_dim1], ldb, &
		    c__, &s);

#line 1009 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1009 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 1010 "MB03XU.f"
		yg[*k + i__ + j * yg_dim1] = 0.;
#line 1011 "MB03XU.f"
/* L100: */
#line 1011 "MB03XU.f"
	    }
#line 1012 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1012 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 1013 "MB03XU.f"
		yg[*k + i__ + (*nb + j) * yg_dim1] = 0.;
#line 1014 "MB03XU.f"
/* L110: */
#line 1014 "MB03XU.f"
	    }
#line 1015 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1015 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 1016 "MB03XU.f"
		ya[*k + i__ + j * ya_dim1] = 0.;
#line 1017 "MB03XU.f"
/* L120: */
#line 1017 "MB03XU.f"
	    }
#line 1018 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1018 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 1019 "MB03XU.f"
		ya[*k + i__ + (*nb + j) * ya_dim1] = 0.;
#line 1020 "MB03XU.f"
/* L130: */
#line 1020 "MB03XU.f"
	    }

/*           Update XG with second Householder reflection. */

#line 1024 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1024 "MB03XU.f"
	    i__3 = *k + *n;
#line 1024 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &g[*k + i__ + g_dim1], 
		    ldg, &a[i__ + (*k + i__) * a_dim1], lda, &c_b9, &xg[(i__ 
		    + *nb) * xg_dim1 + 1], &c__1, (ftnlen)9);
#line 1026 "MB03XU.f"
	    i__2 = *k + *n;
#line 1026 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xg[xg_offset], ldxg, &
		    dwork[nb2 + 1], &c__1, &c_b7, &xg[(i__ + *nb) * xg_dim1 + 
		    1], &c__1, (ftnlen)12);
#line 1028 "MB03XU.f"
	    i__2 = *k + *n;
#line 1028 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1028 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &dwork[nb3 + 1], &c__1, &c_b7, &xg[(i__ + *nb) * 
		    xg_dim1 + 1], &c__1, (ftnlen)12);
#line 1030 "MB03XU.f"
	    i__2 = *n - i__;
#line 1030 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1030 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + 1 + 
		    yg_dim1], ldyg, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 1032 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1032 "MB03XU.f"
	    i__3 = *n - i__;
#line 1032 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 + (
		    i__ + *nb) * xg_dim1], &c__1, (ftnlen)9);
#line 1034 "MB03XU.f"
	    i__2 = *n - i__;
#line 1034 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1034 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + 1 + nb1 * 
		    yg_dim1], ldyg, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 1036 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1036 "MB03XU.f"
	    i__3 = *n - i__;
#line 1036 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ 
		    + 1 + (i__ + *nb) * xg_dim1], &c__1, (ftnlen)9);
#line 1038 "MB03XU.f"
	    i__2 = *k + *n;
#line 1038 "MB03XU.f"
	    d__1 = -taul[i__];
#line 1038 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xg[(i__ + *nb) * xg_dim1 + 1], &c__1);

/*           Update G(k+i,:). */

#line 1042 "MB03XU.f"
	    i__2 = *k + *n;
#line 1042 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &xg[(i__ + *nb) * xg_dim1 + 1], &c__1, &g[*k 
		    + i__ + g_dim1], ldg);

/*           Update XB with second Householder reflection. */

#line 1046 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1046 "MB03XU.f"
	    i__3 = *k + *n;
#line 1046 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[i__ + b_dim1], ldb, &
		    a[i__ + (*k + i__) * a_dim1], lda, &c_b9, &xb[(i__ + *nb) 
		    * xb_dim1 + 1], &c__1, (ftnlen)9);
#line 1048 "MB03XU.f"
	    i__2 = *k + *n;
#line 1048 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xb[xb_offset], ldxb, &
		    dwork[nb2 + 1], &c__1, &c_b7, &xb[(i__ + *nb) * xb_dim1 + 
		    1], &c__1, (ftnlen)12);
#line 1050 "MB03XU.f"
	    i__2 = *k + *n;
#line 1050 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1050 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &dwork[nb3 + 1], &c__1, &c_b7, &xb[(i__ + *nb) * 
		    xb_dim1 + 1], &c__1, (ftnlen)12);
#line 1052 "MB03XU.f"
	    i__2 = *n - i__;
#line 1052 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1052 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + yb_dim1], 
		    ldyb, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
#line 1054 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1054 "MB03XU.f"
	    i__3 = *n - i__;
#line 1054 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 + (
		    i__ + *nb) * xb_dim1], &c__1, (ftnlen)9);
#line 1056 "MB03XU.f"
	    i__2 = *n - i__;
#line 1056 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1056 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &a[i__ + (*k + i__ + 1) * a_dim1], lda, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 1058 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1058 "MB03XU.f"
	    i__3 = *n - i__;
#line 1058 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ 
		    + 1 + (i__ + *nb) * xb_dim1], &c__1, (ftnlen)9);
#line 1060 "MB03XU.f"
	    i__2 = *k + *n;
#line 1060 "MB03XU.f"
	    d__1 = -taul[i__];
#line 1060 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xb[(i__ + *nb) * xb_dim1 + 1], &c__1);

/*           Update B(i,:). */

#line 1064 "MB03XU.f"
	    i__2 = *k + *n;
#line 1064 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &xb[(i__ + *nb) * xb_dim1 + 1], &c__1, &b[
		    i__ + b_dim1], ldb);

#line 1066 "MB03XU.f"
	    a[i__ + (*k + i__) * a_dim1] = temp;
#line 1067 "MB03XU.f"
	    q[i__ + i__ * q_dim1] = tauq;
#line 1068 "MB03XU.f"
	    csl[(i__ << 1) - 1] = c__;
#line 1069 "MB03XU.f"
	    csl[i__ * 2] = s;

/*           Transform first rows of Q and B. */

#line 1073 "MB03XU.f"
	    alpha = q[i__ + (i__ + 1) * q_dim1];
#line 1074 "MB03XU.f"
	    i__2 = *n - i__;
#line 1074 "MB03XU.f"
	    dlarfg_(&i__2, &alpha, &q[i__ + (i__ + 2) * q_dim1], ldq, &tauq);
#line 1075 "MB03XU.f"
	    q[i__ + (i__ + 1) * q_dim1] = 1.;
#line 1076 "MB03XU.f"
	    i__2 = *n - i__;
#line 1076 "MB03XU.f"
	    temp = -tauq * ddot_(&i__2, &q[i__ + (i__ + 1) * q_dim1], ldq, &b[
		    i__ + (*k + i__ + 1) * b_dim1], ldb);
#line 1077 "MB03XU.f"
	    i__2 = *n - i__;
#line 1077 "MB03XU.f"
	    daxpy_(&i__2, &temp, &q[i__ + (i__ + 1) * q_dim1], ldq, &b[i__ + (
		    *k + i__ + 1) * b_dim1], ldb);
#line 1078 "MB03XU.f"
	    temp = b[i__ + (*k + i__ + 1) * b_dim1];
#line 1079 "MB03XU.f"
	    dlartg_(&temp, &alpha, &c__, &s, &b[i__ + (*k + i__ + 1) * b_dim1]
		    );
#line 1080 "MB03XU.f"
	    s = -s;
#line 1081 "MB03XU.f"
	    i__2 = *n - i__;
#line 1081 "MB03XU.f"
	    dlarfg_(&i__2, &b[i__ + (*k + i__ + 1) * b_dim1], &b[i__ + (*k + 
		    i__ + 2) * b_dim1], ldb, &taur[i__]);
#line 1082 "MB03XU.f"
	    temp = b[i__ + (*k + i__ + 1) * b_dim1];
#line 1083 "MB03XU.f"
	    b[i__ + (*k + i__ + 1) * b_dim1] = 1.;

/*           Update YB with first Householder reflection. */

#line 1087 "MB03XU.f"
	    i__2 = *n - i__;
#line 1087 "MB03XU.f"
	    i__3 = *n - i__;
#line 1087 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[i__ + 1 + (*k + 
		    i__ + 1) * b_dim1], ldb, &q[i__ + (i__ + 1) * q_dim1], 
		    ldq, &c_b9, &yb[i__ + 1 + i__ * yb_dim1], &c__1, (ftnlen)
		    12);
#line 1089 "MB03XU.f"
	    i__2 = *n - i__;
#line 1089 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 1 + 
		    xb_dim1], ldxb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yb[i__ * yb_dim1 + 1], &c__1, (ftnlen)9);
#line 1091 "MB03XU.f"
	    i__2 = *n - i__;
#line 1091 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yb[i__ * yb_dim1 + 1], &c__1, &c_b7, &yb[i__ + 1 + 
		    i__ * yb_dim1], &c__1, (ftnlen)12);
#line 1093 "MB03XU.f"
	    i__2 = *n - i__;
#line 1093 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 1 + nb1 * 
		    xb_dim1], ldxb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yb[i__ * yb_dim1 + 1], &c__1, (ftnlen)9);
#line 1095 "MB03XU.f"
	    i__2 = *n - i__;
#line 1095 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &yb[i__ * yb_dim1 + 1], &c__1, &c_b7, &
		    yb[i__ + 1 + i__ * yb_dim1], &c__1, (ftnlen)9);
#line 1097 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1097 "MB03XU.f"
	    i__3 = *n - i__;
#line 1097 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 
		    + 1], ldq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &
		    dwork[1], &c__1, (ftnlen)12);
#line 1099 "MB03XU.f"
	    i__2 = *n - i__;
#line 1099 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1099 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + yb_dim1]
		    , ldyb, &dwork[1], &c__1, &c_b7, &yb[i__ + 1 + i__ * 
		    yb_dim1], &c__1, (ftnlen)12);
#line 1101 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1101 "MB03XU.f"
	    i__3 = *n - i__;
#line 1101 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &q[i__ + (i__ + 1) * q_dim1], ldq, &
		    c_b9, &dwork[*nb + 1], &c__1, (ftnlen)12);
#line 1103 "MB03XU.f"
	    i__2 = *n - i__;
#line 1103 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1103 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &dwork[*nb + 1], &c__1, &c_b7, &yb[i__ + 
		    1 + i__ * yb_dim1], &c__1, (ftnlen)12);
#line 1105 "MB03XU.f"
	    i__2 = *n - i__;
#line 1105 "MB03XU.f"
	    d__1 = -tauq;
#line 1105 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yb[i__ + 1 + i__ * yb_dim1], &c__1);

/*           Update B(i+1:n,k+i+1). */

#line 1109 "MB03XU.f"
	    i__2 = *n - i__;
#line 1109 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xb[*k + i__ + 1 + xb_dim1], ldxb, &c_b7, &b[i__ + 1 
		    + (*k + i__ + 1) * b_dim1], &c__1, (ftnlen)12);
#line 1111 "MB03XU.f"
	    i__2 = *n - i__;
#line 1111 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &xb[*k + i__ + 1 + nb1 * xb_dim1], ldxb,
		     &c_b7, &b[i__ + 1 + (*k + i__ + 1) * b_dim1], &c__1, (
		    ftnlen)9);
#line 1113 "MB03XU.f"
	    i__2 = *n - i__;
#line 1113 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yb[i__ + 1 + yb_dim1],
		     ldyb, &q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &b[i__ + 
		    1 + (*k + i__ + 1) * b_dim1], &c__1, (ftnlen)12);
#line 1115 "MB03XU.f"
	    i__2 = *n - i__;
#line 1115 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1115 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &b[(*k + i__ + 1) * b_dim1 + 1], &c__1, &
		    c_b7, &b[i__ + 1 + (*k + i__ + 1) * b_dim1], &c__1, (
		    ftnlen)12);

/*           Update YQ with first Householder reflection. */

#line 1120 "MB03XU.f"
	    i__2 = *n - i__;
#line 1120 "MB03XU.f"
	    i__3 = *n - i__;
#line 1120 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[i__ + 1 + (i__ + 1)
		     * q_dim1], ldq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9,
		     &yq[i__ + 1 + i__ * yq_dim1], &c__1, (ftnlen)12);
#line 1122 "MB03XU.f"
	    i__2 = *n - i__;
#line 1122 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1], 
		    ldxq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &yq[i__ * 
		    yq_dim1 + 1], &c__1, (ftnlen)9);
#line 1124 "MB03XU.f"
	    i__2 = *n - i__;
#line 1124 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yq[i__ * yq_dim1 + 1], &c__1, &c_b7, &yq[i__ + 1 + 
		    i__ * yq_dim1], &c__1, (ftnlen)12);
#line 1126 "MB03XU.f"
	    i__2 = *n - i__;
#line 1126 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yq[i__ * yq_dim1 + 1], &c__1, (ftnlen)9);
#line 1128 "MB03XU.f"
	    i__2 = *n - i__;
#line 1128 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &yq[i__ * yq_dim1 + 1], &c__1, &c_b7, &
		    yq[i__ + 1 + i__ * yq_dim1], &c__1, (ftnlen)9);
#line 1130 "MB03XU.f"
	    i__2 = *n - i__;
#line 1130 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1130 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + yq_dim1]
		    , ldyq, &dwork[1], &c__1, &c_b7, &yq[i__ + 1 + i__ * 
		    yq_dim1], &c__1, (ftnlen)12);
#line 1132 "MB03XU.f"
	    i__2 = *n - i__;
#line 1132 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1132 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &dwork[*nb + 1], &c__1, &c_b7, &yq[i__ + 
		    1 + i__ * yq_dim1], &c__1, (ftnlen)12);
#line 1134 "MB03XU.f"
	    i__2 = *n - i__;
#line 1134 "MB03XU.f"
	    d__1 = -tauq;
#line 1134 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yq[i__ + 1 + i__ * yq_dim1], &c__1);

/*           Update Q(i+1:n,i+1). */

#line 1138 "MB03XU.f"
	    i__2 = *n - i__;
#line 1138 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xq[i__ + 1 + xq_dim1], ldxq, &c_b7, &q[i__ + 1 + (
		    i__ + 1) * q_dim1], &c__1, (ftnlen)12);
#line 1140 "MB03XU.f"
	    i__2 = *n - i__;
#line 1140 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &xq[i__ + 1 + nb1 * xq_dim1], ldxq, &
		    c_b7, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)9);
#line 1142 "MB03XU.f"
	    i__2 = *n - i__;
#line 1142 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yq[i__ + 1 + yq_dim1],
		     ldyq, &q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &q[i__ + 
		    1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12);
#line 1144 "MB03XU.f"
	    i__2 = *n - i__;
#line 1144 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1144 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &b[(*k + i__ + 1) * b_dim1 + 1], &c__1, &
		    c_b7, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12)
		    ;

/*           Apply rotation to [ Q(i+1:n,i+1), B(i+1:n,k+i+1) ]. */

#line 1149 "MB03XU.f"
	    i__2 = *n - i__;
#line 1149 "MB03XU.f"
	    drot_(&i__2, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, &b[i__ + 1 
		    + (*k + i__ + 1) * b_dim1], &c__1, &c__, &s);
#line 1150 "MB03XU.f"
	    i__2 = i__;
#line 1150 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 1151 "MB03XU.f"
		xb[*k + i__ + 1 + j * xb_dim1] = 0.;
#line 1152 "MB03XU.f"
/* L140: */
#line 1152 "MB03XU.f"
	    }
#line 1153 "MB03XU.f"
	    i__2 = i__;
#line 1153 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 1154 "MB03XU.f"
		xb[*k + i__ + 1 + (*nb + j) * xb_dim1] = 0.;
#line 1155 "MB03XU.f"
/* L150: */
#line 1155 "MB03XU.f"
	    }

/*           Update YB with second Householder reflection. */

#line 1159 "MB03XU.f"
	    i__2 = *n - i__;
#line 1159 "MB03XU.f"
	    i__3 = *n - i__;
#line 1159 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[i__ + 1 + (*k + 
		    i__ + 1) * b_dim1], ldb, &b[i__ + (*k + i__ + 1) * b_dim1]
		    , ldb, &c_b9, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1,
		     (ftnlen)12);
#line 1161 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1161 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 2 + 
		    xb_dim1], ldxb, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, (ftnlen)9);
#line 1163 "MB03XU.f"
	    i__2 = *n - i__;
#line 1163 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, &c_b7, &yb[
		    i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (ftnlen)12);
#line 1165 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1165 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 2 + nb1 * 
		    xb_dim1], ldxb, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, (ftnlen)9);
#line 1167 "MB03XU.f"
	    i__2 = *n - i__;
#line 1167 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, &
		    c_b7, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (
		    ftnlen)9);
#line 1169 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1169 "MB03XU.f"
	    dgemv_("No transpose", &i__, &i__2, &c_b7, &q[(i__ + 2) * q_dim1 
		    + 1], ldq, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &c_b9, 
		    &dwork[nb2 + 1], &c__1, (ftnlen)12);
#line 1171 "MB03XU.f"
	    i__2 = *n - i__;
#line 1171 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yb[i__ + 1 + yb_dim1],
		     ldyb, &dwork[nb2 + 1], &c__1, &c_b7, &yb[i__ + 1 + (i__ 
		    + *nb) * yb_dim1], &c__1, (ftnlen)12);
#line 1173 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1173 "MB03XU.f"
	    i__3 = *n - i__ - 1;
#line 1173 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 2) * 
		    b_dim1 + 1], ldq, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, 
		    &c_b9, &dwork[nb3 + 1], &c__1, (ftnlen)12);
#line 1175 "MB03XU.f"
	    i__2 = *n - i__;
#line 1175 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1175 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &dwork[nb3 + 1], &c__1, &c_b7, &yb[i__ + 
		    1 + (i__ + *nb) * yb_dim1], &c__1, (ftnlen)12);
#line 1177 "MB03XU.f"
	    i__2 = *n - i__;
#line 1177 "MB03XU.f"
	    d__1 = -taur[i__];
#line 1177 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1);

/*           Update B(i+1:n,k+i+1). */

#line 1181 "MB03XU.f"
	    i__2 = *n - i__;
#line 1181 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, 
		    &b[i__ + 1 + (*k + i__ + 1) * b_dim1], &c__1);

/*           Update YQ with second Householder reflection. */

#line 1185 "MB03XU.f"
	    i__2 = *n - i__;
#line 1185 "MB03XU.f"
	    i__3 = *n - i__;
#line 1185 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[i__ + 1 + (i__ + 1)
		     * q_dim1], ldq, &b[i__ + (*k + i__ + 1) * b_dim1], ldb, &
		    c_b9, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (
		    ftnlen)12);
#line 1187 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1187 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 2 + xq_dim1], 
		    ldxq, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &c_b9, &yq[(
		    i__ + *nb) * yq_dim1 + 1], &c__1, (ftnlen)9);
#line 1189 "MB03XU.f"
	    i__2 = *n - i__;
#line 1189 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, &c_b7, &yq[
		    i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (ftnlen)12);
#line 1191 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1191 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 2 + nb1 * 
		    xq_dim1], ldxq, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, (ftnlen)9);
#line 1193 "MB03XU.f"
	    i__2 = *n - i__;
#line 1193 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, &
		    c_b7, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (
		    ftnlen)9);
#line 1195 "MB03XU.f"
	    i__2 = *n - i__;
#line 1195 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yq[i__ + 1 + yq_dim1],
		     ldyq, &dwork[nb2 + 1], &c__1, &c_b7, &yq[i__ + 1 + (i__ 
		    + *nb) * yq_dim1], &c__1, (ftnlen)12);
#line 1197 "MB03XU.f"
	    i__2 = *n - i__;
#line 1197 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1197 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &dwork[nb3 + 1], &c__1, &c_b7, &yq[i__ + 
		    1 + (i__ + *nb) * yq_dim1], &c__1, (ftnlen)12);
#line 1199 "MB03XU.f"
	    i__2 = *n - i__;
#line 1199 "MB03XU.f"
	    d__1 = -taur[i__];
#line 1199 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1);

/*           Update Q(i+1:n,i+1). */

#line 1203 "MB03XU.f"
	    i__2 = *n - i__;
#line 1203 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, 
		    &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1);

/*           Update YA with first Householder reflection. */

#line 1207 "MB03XU.f"
	    i__2 = *n - i__;
#line 1207 "MB03XU.f"
	    i__3 = *k + *n;
#line 1207 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[i__ + 1 + a_dim1], 
		    lda, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &ya[i__ * 
		    ya_dim1 + 1], &c__1, (ftnlen)9);
#line 1209 "MB03XU.f"
	    i__2 = *n - i__;
#line 1209 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1], 
		    ldxa, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &dwork[
		    pdw], &c__1, (ftnlen)9);
#line 1211 "MB03XU.f"
	    i__2 = *n - i__;
#line 1211 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 + i__ * 
		    ya_dim1], &c__1, (ftnlen)12);
#line 1213 "MB03XU.f"
	    i__2 = *n - i__;
#line 1213 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
#line 1215 "MB03XU.f"
	    i__2 = *n - i__;
#line 1215 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ 
		    + 1 + i__ * ya_dim1], &c__1, (ftnlen)9);
#line 1217 "MB03XU.f"
	    i__2 = *k + *n;
#line 1217 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1217 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[ya_offset], ldya, 
		    &dwork[1], &c__1, &c_b7, &ya[i__ * ya_dim1 + 1], &c__1, (
		    ftnlen)12);
#line 1219 "MB03XU.f"
	    i__2 = *k + *n;
#line 1219 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1219 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &dwork[*nb + 1], &c__1, &c_b7, &ya[i__ * ya_dim1 
		    + 1], &c__1, (ftnlen)12);
#line 1221 "MB03XU.f"
	    i__2 = *k + *n;
#line 1221 "MB03XU.f"
	    d__1 = -tauq;
#line 1221 "MB03XU.f"
	    dscal_(&i__2, &d__1, &ya[i__ * ya_dim1 + 1], &c__1);

/*           Update A(i+1,1:k+n). */

#line 1225 "MB03XU.f"
	    i__2 = *n - i__;
#line 1225 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xa[i__ + 1 + xa_dim1], ldxa, &c_b7, &a[i__ + 1 + (*
		    k + i__ + 1) * a_dim1], lda, (ftnlen)12);
#line 1227 "MB03XU.f"
	    i__2 = *n - i__;
#line 1227 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &xa[i__ + 1 + nb1 * xa_dim1], ldxa, &
		    c_b7, &a[i__ + 1 + (*k + i__ + 1) * a_dim1], lda, (ftnlen)
		    9);
#line 1229 "MB03XU.f"
	    i__2 = *k + *n;
#line 1229 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &ya[ya_offset], ldya, &
		    q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &a[i__ + 1 + 
		    a_dim1], lda, (ftnlen)12);
#line 1231 "MB03XU.f"
	    i__2 = *k + *n;
#line 1231 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1231 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &b[(*k + i__ + 1) * b_dim1 + 1], &c__1, &c_b7, &a[
		    i__ + 1 + a_dim1], lda, (ftnlen)12);

/*           Update YG with first Householder reflection. */

#line 1236 "MB03XU.f"
	    i__2 = *k + *n;
#line 1236 "MB03XU.f"
	    i__3 = *n - i__;
#line 1236 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], ldg, &q[i__ + (i__ + 1) * q_dim1], ldq, &
		    c_b9, &yg[i__ * yg_dim1 + 1], &c__1, (ftnlen)12);
#line 1238 "MB03XU.f"
	    i__2 = *n - i__;
#line 1238 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 1 + 
		    xg_dim1], ldxg, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
#line 1240 "MB03XU.f"
	    i__2 = *n - i__;
#line 1240 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 + i__ * 
		    yg_dim1], &c__1, (ftnlen)12);
#line 1242 "MB03XU.f"
	    i__2 = *n - i__;
#line 1242 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 1 + nb1 * 
		    xg_dim1], ldxg, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
#line 1244 "MB03XU.f"
	    i__2 = *n - i__;
#line 1244 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ 
		    + 1 + i__ * yg_dim1], &c__1, (ftnlen)9);
#line 1246 "MB03XU.f"
	    i__2 = *k + *n;
#line 1246 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1246 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[yg_offset], ldyg, 
		    &dwork[1], &c__1, &c_b7, &yg[i__ * yg_dim1 + 1], &c__1, (
		    ftnlen)12);
#line 1248 "MB03XU.f"
	    i__2 = *k + *n;
#line 1248 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1248 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &dwork[*nb + 1], &c__1, &c_b7, &yg[i__ * yg_dim1 
		    + 1], &c__1, (ftnlen)12);
#line 1250 "MB03XU.f"
	    i__2 = *k + *n;
#line 1250 "MB03XU.f"
	    d__1 = -tauq;
#line 1250 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yg[i__ * yg_dim1 + 1], &c__1);

/*           Update G(1:k+n,k+i+1). */

#line 1254 "MB03XU.f"
	    i__2 = *n - i__;
#line 1254 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xg[*k + i__ + 1 + xg_dim1], ldxg, &c_b7, &g[*k + 
		    i__ + 1 + (*k + i__ + 1) * g_dim1], &c__1, (ftnlen)12);
#line 1256 "MB03XU.f"
	    i__2 = *n - i__;
#line 1256 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &xg[*k + i__ + 1 + nb1 * xg_dim1], ldxg,
		     &c_b7, &g[*k + i__ + 1 + (*k + i__ + 1) * g_dim1], &c__1,
		     (ftnlen)9);
#line 1258 "MB03XU.f"
	    i__2 = *k + *n;
#line 1258 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yg[yg_offset], ldyg, &
		    q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &g[(*k + i__ + 1)
		     * g_dim1 + 1], &c__1, (ftnlen)12);
#line 1260 "MB03XU.f"
	    i__2 = *k + *n;
#line 1260 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1260 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &b[(*k + i__ + 1) * b_dim1 + 1], &c__1, &c_b7, &g[
		    (*k + i__ + 1) * g_dim1 + 1], &c__1, (ftnlen)12);
#line 1262 "MB03XU.f"
	    i__2 = i__;
#line 1262 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 1263 "MB03XU.f"
		xg[*k + i__ + 1 + j * xg_dim1] = 0.;
#line 1264 "MB03XU.f"
/* L160: */
#line 1264 "MB03XU.f"
	    }
#line 1265 "MB03XU.f"
	    i__2 = i__;
#line 1265 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 1266 "MB03XU.f"
		xg[*k + i__ + 1 + (*nb + j) * xg_dim1] = 0.;
#line 1267 "MB03XU.f"
/* L170: */
#line 1267 "MB03XU.f"
	    }

/*           Apply rotation to [ A(i+1,1:k+n)', G(1:k+n,k+i+1) ]. */

#line 1271 "MB03XU.f"
	    i__2 = *k + *n;
#line 1271 "MB03XU.f"
	    drot_(&i__2, &a[i__ + 1 + a_dim1], lda, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], &c__1, &c__, &s);

/*           Update YA with second Householder reflection. */

#line 1275 "MB03XU.f"
	    i__2 = *n - i__;
#line 1275 "MB03XU.f"
	    i__3 = *k + *n;
#line 1275 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[i__ + 1 + a_dim1], 
		    lda, &b[i__ + (*k + i__ + 1) * b_dim1], ldb, &c_b9, &ya[(
		    i__ + *nb) * ya_dim1 + 1], &c__1, (ftnlen)9);
#line 1277 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1277 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 2 + xa_dim1], 
		    ldxa, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
#line 1279 "MB03XU.f"
	    i__2 = *n - i__;
#line 1279 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 + (i__ + 
		    *nb) * ya_dim1], &c__1, (ftnlen)12);
#line 1281 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1281 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 2 + nb1 * 
		    xa_dim1], ldxa, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 1283 "MB03XU.f"
	    i__2 = *n - i__;
#line 1283 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ 
		    + 1 + (i__ + *nb) * ya_dim1], &c__1, (ftnlen)9);
#line 1285 "MB03XU.f"
	    i__2 = *k + *n;
#line 1285 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &ya[ya_offset], ldya, &
		    dwork[nb2 + 1], &c__1, &c_b7, &ya[(i__ + *nb) * ya_dim1 + 
		    1], &c__1, (ftnlen)12);
#line 1287 "MB03XU.f"
	    i__2 = *k + *n;
#line 1287 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1287 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &dwork[nb3 + 1], &c__1, &c_b7, &ya[(i__ + *nb) * 
		    ya_dim1 + 1], &c__1, (ftnlen)12);
#line 1289 "MB03XU.f"
	    i__2 = *k + *n;
#line 1289 "MB03XU.f"
	    d__1 = -taur[i__];
#line 1289 "MB03XU.f"
	    dscal_(&i__2, &d__1, &ya[(i__ + *nb) * ya_dim1 + 1], &c__1);

/*           Update A(i+1,1:k+n). */

#line 1293 "MB03XU.f"
	    i__2 = *k + *n;
#line 1293 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &ya[(i__ + *nb) * ya_dim1 + 1], &c__1, &a[
		    i__ + 1 + a_dim1], lda);

/*           Update YG with second Householder reflection. */

#line 1297 "MB03XU.f"
	    i__2 = *k + *n;
#line 1297 "MB03XU.f"
	    i__3 = *n - i__;
#line 1297 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], ldg, &b[i__ + (*k + i__ + 1) * b_dim1], ldb, 
		    &c_b9, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1, (ftnlen)12);
#line 1299 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1299 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 2 + 
		    xg_dim1], ldxg, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 1301 "MB03XU.f"
	    i__2 = *n - i__;
#line 1301 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 + (i__ + 
		    *nb) * yg_dim1], &c__1, (ftnlen)12);
#line 1303 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1303 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 2 + nb1 * 
		    xg_dim1], ldxg, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 1305 "MB03XU.f"
	    i__2 = *n - i__;
#line 1305 "MB03XU.f"
	    dgemv_("Transpose", &i__, &i__2, &c_b7, &a[(*k + i__ + 1) * 
		    a_dim1 + 1], lda, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ 
		    + 1 + (i__ + *nb) * yg_dim1], &c__1, (ftnlen)9);
#line 1307 "MB03XU.f"
	    i__2 = *k + *n;
#line 1307 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yg[yg_offset], ldyg, &
		    dwork[nb2 + 1], &c__1, &c_b7, &yg[(i__ + *nb) * yg_dim1 + 
		    1], &c__1, (ftnlen)12);
#line 1309 "MB03XU.f"
	    i__2 = *k + *n;
#line 1309 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1309 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &dwork[nb3 + 1], &c__1, &c_b7, &yg[(i__ + *nb) * 
		    yg_dim1 + 1], &c__1, (ftnlen)12);
#line 1311 "MB03XU.f"
	    i__2 = *k + *n;
#line 1311 "MB03XU.f"
	    d__1 = -taur[i__];
#line 1311 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1);

/*           Update G(1:k+n,k+i+1). */

#line 1315 "MB03XU.f"
	    i__2 = *k + *n;
#line 1315 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1, &g[(*
		    k + i__ + 1) * g_dim1 + 1], &c__1);

#line 1317 "MB03XU.f"
	    b[i__ + (*k + i__ + 1) * b_dim1] = temp;
#line 1318 "MB03XU.f"
	    q[i__ + (i__ + 1) * q_dim1] = tauq;
#line 1319 "MB03XU.f"
	    csr[(i__ << 1) - 1] = c__;
#line 1320 "MB03XU.f"
	    csr[i__ * 2] = s;
#line 1321 "MB03XU.f"
/* L180: */
#line 1321 "MB03XU.f"
	}

#line 1323 "MB03XU.f"
    } else if (*ltrb) {
#line 1324 "MB03XU.f"
	i__1 = *nb;
#line 1324 "MB03XU.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Transform first columns of A and Q. See routine MB04TS. */

#line 1328 "MB03XU.f"
	    alpha = q[i__ + i__ * q_dim1];
#line 1329 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1329 "MB03XU.f"
	    dlarfg_(&i__2, &alpha, &q[i__ + 1 + i__ * q_dim1], &c__1, &tauq);
#line 1330 "MB03XU.f"
	    q[i__ + i__ * q_dim1] = 1.;
#line 1331 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1331 "MB03XU.f"
	    temp = -tauq * ddot_(&i__2, &q[i__ + i__ * q_dim1], &c__1, &a[*k 
		    + i__ + i__ * a_dim1], &c__1);
#line 1332 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1332 "MB03XU.f"
	    daxpy_(&i__2, &temp, &q[i__ + i__ * q_dim1], &c__1, &a[*k + i__ + 
		    i__ * a_dim1], &c__1);
#line 1333 "MB03XU.f"
	    temp = a[*k + i__ + i__ * a_dim1];
#line 1334 "MB03XU.f"
	    dlartg_(&temp, &alpha, &c__, &s, &a[*k + i__ + i__ * a_dim1]);
#line 1335 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1335 "MB03XU.f"
	    dlarfg_(&i__2, &a[*k + i__ + i__ * a_dim1], &a[*k + i__ + 1 + i__ 
		    * a_dim1], &c__1, &taul[i__]);
#line 1336 "MB03XU.f"
	    temp = a[*k + i__ + i__ * a_dim1];
#line 1337 "MB03XU.f"
	    a[*k + i__ + i__ * a_dim1] = 1.;

/*           Update XQ with first Householder reflection. */

#line 1341 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1341 "MB03XU.f"
	    i__3 = *n - i__;
#line 1341 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[
		    i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)9);
#line 1343 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1343 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1343 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + q_dim1], ldq, &
		    q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[1], &c__1, (
		    ftnlen)9);
#line 1345 "MB03XU.f"
	    i__2 = *n - i__;
#line 1345 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1345 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + xq_dim1]
		    , ldxq, &dwork[1], &c__1, &c_b7, &xq[i__ + 1 + i__ * 
		    xq_dim1], &c__1, (ftnlen)12);
#line 1347 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1347 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1347 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[*k + i__ + a_dim1], 
		    lda, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[*nb + 1]
		    , &c__1, (ftnlen)9);
#line 1349 "MB03XU.f"
	    i__2 = *n - i__;
#line 1349 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1349 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[*nb + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + i__ * xq_dim1], &c__1, (ftnlen)12);
#line 1351 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1351 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1351 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + yq_dim1], ldyq,
		     &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[i__ * xq_dim1 
		    + 1], &c__1, (ftnlen)9);
#line 1353 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1353 "MB03XU.f"
	    i__3 = *n - i__;
#line 1353 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xq[i__ * xq_dim1 + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + i__ * xq_dim1], &c__1, (ftnlen)9);
#line 1355 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1355 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1355 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + nb1 * yq_dim1],
		     ldyq, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[(i__ + *
		    nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
#line 1357 "MB03XU.f"
	    i__2 = *n - i__;
#line 1357 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1357 "MB03XU.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &
		    c_b7, &xq[i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)12);
#line 1359 "MB03XU.f"
	    i__2 = *n - i__;
#line 1359 "MB03XU.f"
	    d__1 = -tauq;
#line 1359 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + i__ * xq_dim1], &c__1);

/*           Update Q(i,i+1:n). */

#line 1363 "MB03XU.f"
	    i__2 = *n - i__;
#line 1363 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1],
		     ldxq, &q[i__ + q_dim1], ldq, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, (ftnlen)12);
#line 1365 "MB03XU.f"
	    i__2 = *n - i__;
#line 1365 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1365 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &a[*k + i__ + a_dim1], lda, &c_b7, &q[i__ 
		    + (i__ + 1) * q_dim1], ldq, (ftnlen)12);
#line 1367 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1367 "MB03XU.f"
	    i__3 = *n - i__;
#line 1367 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yq[i__ + yq_dim1], ldyq, &c_b7, &q[i__ + (i__ + 
		    1) * q_dim1], ldq, (ftnlen)9);
#line 1369 "MB03XU.f"
	    i__2 = *n - i__;
#line 1369 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1369 "MB03XU.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &yq[i__ + nb1 * yq_dim1], ldyq, &c_b7, &q[
		    i__ + (i__ + 1) * q_dim1], ldq, (ftnlen)12);

/*           Update XA with first Householder reflection. */

#line 1374 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1374 "MB03XU.f"
	    i__3 = *n - i__;
#line 1374 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[*k + i__ + (i__ + 1) *
		     a_dim1], lda, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[
		    i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)9);
#line 1376 "MB03XU.f"
	    i__2 = *n - i__;
#line 1376 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1376 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + xa_dim1]
		    , ldxa, &dwork[1], &c__1, &c_b7, &xa[i__ + 1 + i__ * 
		    xa_dim1], &c__1, (ftnlen)12);
#line 1378 "MB03XU.f"
	    i__2 = *n - i__;
#line 1378 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1378 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[*nb + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + i__ * xa_dim1], &c__1, (ftnlen)12);
#line 1380 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1380 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1380 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + ya_dim1], 
		    ldya, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[i__ * 
		    xa_dim1 + 1], &c__1, (ftnlen)9);
#line 1382 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1382 "MB03XU.f"
	    i__3 = *n - i__;
#line 1382 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xa[i__ * xa_dim1 + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + i__ * xa_dim1], &c__1, (ftnlen)9);
#line 1384 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1384 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1384 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + nb1 * 
		    ya_dim1], ldya, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[
		    i__ * xa_dim1 + 1], &c__1, (ftnlen)9);
#line 1386 "MB03XU.f"
	    i__2 = *n - i__;
#line 1386 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1386 "MB03XU.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &xa[i__ * xa_dim1 + 1], &c__1, &c_b7, &xa[
		    i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)12);
#line 1388 "MB03XU.f"
	    i__2 = *n - i__;
#line 1388 "MB03XU.f"
	    d__1 = -tauq;
#line 1388 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + i__ * xa_dim1], &c__1);

/*           Update A(k+i,i+1:n). */

#line 1392 "MB03XU.f"
	    i__2 = *n - i__;
#line 1392 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1],
		     ldxa, &q[i__ + q_dim1], ldq, &c_b7, &a[*k + i__ + (i__ + 
		    1) * a_dim1], lda, (ftnlen)12);
#line 1394 "MB03XU.f"
	    i__2 = *n - i__;
#line 1394 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1394 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &a[*k + i__ + a_dim1], lda, &c_b7, &a[*k 
		    + i__ + (i__ + 1) * a_dim1], lda, (ftnlen)12);
#line 1396 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1396 "MB03XU.f"
	    i__3 = *n - i__;
#line 1396 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &ya[*k + i__ + ya_dim1], ldya, &c_b7, &a[*k + 
		    i__ + (i__ + 1) * a_dim1], lda, (ftnlen)9);
#line 1398 "MB03XU.f"
	    i__2 = *n - i__;
#line 1398 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1398 "MB03XU.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &ya[*k + i__ + nb1 * ya_dim1], ldya, &c_b7, 
		    &a[*k + i__ + (i__ + 1) * a_dim1], lda, (ftnlen)12);

/*           Apply rotation to [ A(k+i,i+1:n); Q(i,i+1:n) ]. */

#line 1403 "MB03XU.f"
	    i__2 = *n - i__;
#line 1403 "MB03XU.f"
	    drot_(&i__2, &a[*k + i__ + (i__ + 1) * a_dim1], lda, &q[i__ + (
		    i__ + 1) * q_dim1], ldq, &c__, &s);

/*           Update XQ with second Householder reflection. */

#line 1407 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1407 "MB03XU.f"
	    i__3 = *n - i__;
#line 1407 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b9, &
		    xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)9);
#line 1409 "MB03XU.f"
	    i__2 = *n - i__;
#line 1409 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], ldq,
		     &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &c_b9, &dwork[
		    nb2 + 1], &c__1, (ftnlen)9);
#line 1411 "MB03XU.f"
	    i__2 = *n - i__;
#line 1411 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1],
		     ldxq, &dwork[nb2 + 1], &c__1, &c_b7, &xq[i__ + 1 + (i__ 
		    + *nb) * xq_dim1], &c__1, (ftnlen)12);
#line 1413 "MB03XU.f"
	    i__2 = *n - i__;
#line 1413 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1413 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[*k + i__ + 1 + a_dim1]
		    , lda, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &c_b9, &
		    dwork[nb3 + 1], &c__1, (ftnlen)9);
#line 1415 "MB03XU.f"
	    i__2 = *n - i__;
#line 1415 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1415 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[nb3 + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)12);
#line 1417 "MB03XU.f"
	    i__2 = *n - i__;
#line 1417 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1417 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + yq_dim1], 
		    ldyq, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &c_b9, &xq[(
		    i__ + *nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
#line 1419 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1419 "MB03XU.f"
	    i__3 = *n - i__;
#line 1419 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &c_b7, &
		    xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)9);
#line 1421 "MB03XU.f"
	    i__2 = *n - i__;
#line 1421 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1421 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
#line 1423 "MB03XU.f"
	    i__2 = *n - i__;
#line 1423 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1423 "MB03XU.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &
		    c_b7, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (
		    ftnlen)12);
#line 1425 "MB03XU.f"
	    i__2 = *n - i__;
#line 1425 "MB03XU.f"
	    d__1 = -taul[i__];
#line 1425 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1);

/*           Update Q(i,i+1:n). */

#line 1429 "MB03XU.f"
	    i__2 = *n - i__;
#line 1429 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, 
		    &q[i__ + (i__ + 1) * q_dim1], ldq);

/*           Update XA with second Householder reflection. */

#line 1433 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1433 "MB03XU.f"
	    i__3 = *n - i__;
#line 1433 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[*k + i__ + (i__ + 1) *
		     a_dim1], lda, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b9, 
		    &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)9);
#line 1435 "MB03XU.f"
	    i__2 = *n - i__;
#line 1435 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1],
		     ldxa, &dwork[nb2 + 1], &c__1, &c_b7, &xa[i__ + 1 + (i__ 
		    + *nb) * xa_dim1], &c__1, (ftnlen)12);
#line 1437 "MB03XU.f"
	    i__2 = *n - i__;
#line 1437 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1437 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[nb3 + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)12);
#line 1439 "MB03XU.f"
	    i__2 = *n - i__;
#line 1439 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1439 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + 1 + 
		    ya_dim1], ldya, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, (ftnlen)9);
#line 1441 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1441 "MB03XU.f"
	    i__3 = *n - i__;
#line 1441 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, &c_b7, &
		    xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)9);
#line 1443 "MB03XU.f"
	    i__2 = *n - i__;
#line 1443 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1443 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + 1 + nb1 * 
		    ya_dim1], ldya, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, (ftnlen)9);
#line 1445 "MB03XU.f"
	    i__2 = *n - i__;
#line 1445 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1445 "MB03XU.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, &
		    c_b7, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (
		    ftnlen)12);
#line 1447 "MB03XU.f"
	    i__2 = *n - i__;
#line 1447 "MB03XU.f"
	    d__1 = -taul[i__];
#line 1447 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1);

/*           Update A(k+i,i+1:n). */

#line 1451 "MB03XU.f"
	    i__2 = *n - i__;
#line 1451 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, 
		    &a[*k + i__ + (i__ + 1) * a_dim1], lda);

/*           Update XG with first Householder reflection. */

#line 1455 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1455 "MB03XU.f"
	    i__3 = *k + *n;
#line 1455 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &g[*k + i__ + g_dim1], 
		    ldg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xg[i__ * 
		    xg_dim1 + 1], &c__1, (ftnlen)9);
#line 1457 "MB03XU.f"
	    i__2 = *k + *n;
#line 1457 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1457 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[xg_offset], ldxg, 
		    &dwork[1], &c__1, &c_b7, &xg[i__ * xg_dim1 + 1], &c__1, (
		    ftnlen)12);
#line 1459 "MB03XU.f"
	    i__2 = *k + *n;
#line 1459 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1459 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &dwork[*nb + 1], &c__1, &c_b7, &xg[i__ * xg_dim1 
		    + 1], &c__1, (ftnlen)12);
#line 1461 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1461 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1461 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + yg_dim1], 
		    ldyg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], &
		    c__1, (ftnlen)9);
#line 1463 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1463 "MB03XU.f"
	    i__3 = *n - i__;
#line 1463 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 + 
		    i__ * xg_dim1], &c__1, (ftnlen)9);
#line 1465 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1465 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1465 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + nb1 * 
		    yg_dim1], ldyg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
#line 1467 "MB03XU.f"
	    i__2 = *n - i__;
#line 1467 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1467 "MB03XU.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 
		    + i__ * xg_dim1], &c__1, (ftnlen)12);
#line 1469 "MB03XU.f"
	    i__2 = *k + *n;
#line 1469 "MB03XU.f"
	    d__1 = -tauq;
#line 1469 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xg[i__ * xg_dim1 + 1], &c__1);

/*           Update G(k+i,:). */

#line 1473 "MB03XU.f"
	    i__2 = *k + *n;
#line 1473 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xg[xg_offset], ldxg, &
		    q[i__ + q_dim1], ldq, &c_b7, &g[*k + i__ + g_dim1], ldg, (
		    ftnlen)12);
#line 1475 "MB03XU.f"
	    i__2 = *k + *n;
#line 1475 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1475 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &a[*k + i__ + a_dim1], lda, &c_b7, &g[*k + i__ + 
		    g_dim1], ldg, (ftnlen)12);
#line 1477 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1477 "MB03XU.f"
	    i__3 = *n - i__;
#line 1477 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yg[*k + i__ + yg_dim1], ldyg, &c_b7, &g[*k + 
		    i__ + (*k + i__ + 1) * g_dim1], ldg, (ftnlen)9);
#line 1479 "MB03XU.f"
	    i__2 = *n - i__;
#line 1479 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1479 "MB03XU.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &yg[*k + i__ + nb1 * yg_dim1], ldyg, &c_b7, 
		    &g[*k + i__ + (*k + i__ + 1) * g_dim1], ldg, (ftnlen)12);

/*           Update XB with first Householder reflection. */

#line 1484 "MB03XU.f"
	    i__2 = *k + *n;
#line 1484 "MB03XU.f"
	    i__3 = *n - i__ + 1;
#line 1484 "MB03XU.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[i__ * b_dim1 + 1], 
		    ldb, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xb[i__ * 
		    xb_dim1 + 1], &c__1, (ftnlen)12);
#line 1486 "MB03XU.f"
	    i__2 = *k + *n;
#line 1486 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1486 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[xb_offset], ldxb, 
		    &dwork[1], &c__1, &c_b7, &xb[i__ * xb_dim1 + 1], &c__1, (
		    ftnlen)12);
#line 1488 "MB03XU.f"
	    i__2 = *k + *n;
#line 1488 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1488 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &dwork[*nb + 1], &c__1, &c_b7, &xb[i__ * xb_dim1 
		    + 1], &c__1, (ftnlen)12);
#line 1490 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1490 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1490 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + yb_dim1], ldyb,
		     &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], &c__1,
		     (ftnlen)9);
#line 1492 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1492 "MB03XU.f"
	    i__3 = *n - i__;
#line 1492 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 + 
		    i__ * xb_dim1], &c__1, (ftnlen)9);
#line 1494 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1494 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1494 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + nb1 * yb_dim1],
		     ldyb, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], 
		    &c__1, (ftnlen)9);
#line 1496 "MB03XU.f"
	    i__2 = *n - i__;
#line 1496 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1496 "MB03XU.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 
		    + i__ * xb_dim1], &c__1, (ftnlen)12);
#line 1498 "MB03XU.f"
	    i__2 = *k + *n;
#line 1498 "MB03XU.f"
	    d__1 = -tauq;
#line 1498 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xb[i__ * xb_dim1 + 1], &c__1);

/*           Update B(:,i). */

#line 1502 "MB03XU.f"
	    i__2 = *k + *n;
#line 1502 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xb[xb_offset], ldxb, &
		    q[i__ + q_dim1], ldq, &c_b7, &b[i__ * b_dim1 + 1], &c__1, 
		    (ftnlen)12);
#line 1504 "MB03XU.f"
	    i__2 = *k + *n;
#line 1504 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1504 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &a[*k + i__ + a_dim1], lda, &c_b7, &b[i__ * 
		    b_dim1 + 1], &c__1, (ftnlen)12);
#line 1506 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1506 "MB03XU.f"
	    i__3 = *n - i__;
#line 1506 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yb[i__ + yb_dim1], ldyb, &c_b7, &b[*k + i__ + 1 
		    + i__ * b_dim1], &c__1, (ftnlen)9);
#line 1508 "MB03XU.f"
	    i__2 = *n - i__;
#line 1508 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1508 "MB03XU.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &yb[i__ + nb1 * yb_dim1], ldyb, &c_b7, &b[*
		    k + i__ + 1 + i__ * b_dim1], &c__1, (ftnlen)12);

/*           Apply rotation to [ G(k+i,:); B(:,i)' ]. */

#line 1513 "MB03XU.f"
	    i__2 = *k + *n;
#line 1513 "MB03XU.f"
	    drot_(&i__2, &g[*k + i__ + g_dim1], ldg, &b[i__ * b_dim1 + 1], &
		    c__1, &c__, &s);

#line 1515 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1515 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 1516 "MB03XU.f"
		yg[*k + i__ + j * yg_dim1] = 0.;
#line 1517 "MB03XU.f"
/* L190: */
#line 1517 "MB03XU.f"
	    }
#line 1518 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1518 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 1519 "MB03XU.f"
		yg[*k + i__ + (*nb + j) * yg_dim1] = 0.;
#line 1520 "MB03XU.f"
/* L200: */
#line 1520 "MB03XU.f"
	    }
#line 1521 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1521 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 1522 "MB03XU.f"
		ya[*k + i__ + j * ya_dim1] = 0.;
#line 1523 "MB03XU.f"
/* L210: */
#line 1523 "MB03XU.f"
	    }
#line 1524 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1524 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 1525 "MB03XU.f"
		ya[*k + i__ + (*nb + j) * ya_dim1] = 0.;
#line 1526 "MB03XU.f"
/* L220: */
#line 1526 "MB03XU.f"
	    }

/*           Update XG with second Householder reflection. */

#line 1530 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1530 "MB03XU.f"
	    i__3 = *k + *n;
#line 1530 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &g[*k + i__ + g_dim1], 
		    ldg, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b9, &xg[(i__ 
		    + *nb) * xg_dim1 + 1], &c__1, (ftnlen)9);
#line 1532 "MB03XU.f"
	    i__2 = *k + *n;
#line 1532 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xg[xg_offset], ldxg, &
		    dwork[nb2 + 1], &c__1, &c_b7, &xg[(i__ + *nb) * xg_dim1 + 
		    1], &c__1, (ftnlen)12);
#line 1534 "MB03XU.f"
	    i__2 = *k + *n;
#line 1534 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1534 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &dwork[nb3 + 1], &c__1, &c_b7, &xg[(i__ + *nb) * 
		    xg_dim1 + 1], &c__1, (ftnlen)12);
#line 1536 "MB03XU.f"
	    i__2 = *n - i__;
#line 1536 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1536 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + 1 + 
		    yg_dim1], ldyg, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 1538 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1538 "MB03XU.f"
	    i__3 = *n - i__;
#line 1538 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 + (
		    i__ + *nb) * xg_dim1], &c__1, (ftnlen)9);
#line 1540 "MB03XU.f"
	    i__2 = *n - i__;
#line 1540 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1540 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + 1 + nb1 * 
		    yg_dim1], ldyg, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 1542 "MB03XU.f"
	    i__2 = *n - i__;
#line 1542 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1542 "MB03XU.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 
		    + (i__ + *nb) * xg_dim1], &c__1, (ftnlen)12);
#line 1544 "MB03XU.f"
	    i__2 = *k + *n;
#line 1544 "MB03XU.f"
	    d__1 = -taul[i__];
#line 1544 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xg[(i__ + *nb) * xg_dim1 + 1], &c__1);

/*           Update G(k+i,:). */

#line 1548 "MB03XU.f"
	    i__2 = *k + *n;
#line 1548 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &xg[(i__ + *nb) * xg_dim1 + 1], &c__1, &g[*k 
		    + i__ + g_dim1], ldg);

/*           Update XB with second Householder reflection. */

#line 1552 "MB03XU.f"
	    i__2 = *k + *n;
#line 1552 "MB03XU.f"
	    i__3 = *n - i__ + 1;
#line 1552 "MB03XU.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[i__ * b_dim1 + 1], 
		    ldb, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b9, &xb[(i__ 
		    + *nb) * xb_dim1 + 1], &c__1, (ftnlen)12);
#line 1554 "MB03XU.f"
	    i__2 = *k + *n;
#line 1554 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xb[xb_offset], ldxb, &
		    dwork[nb2 + 1], &c__1, &c_b7, &xb[(i__ + *nb) * xb_dim1 + 
		    1], &c__1, (ftnlen)12);
#line 1556 "MB03XU.f"
	    i__2 = *k + *n;
#line 1556 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1556 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &dwork[nb3 + 1], &c__1, &c_b7, &xb[(i__ + *nb) * 
		    xb_dim1 + 1], &c__1, (ftnlen)12);
#line 1558 "MB03XU.f"
	    i__2 = *n - i__;
#line 1558 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1558 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + yb_dim1], 
		    ldyb, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
#line 1560 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1560 "MB03XU.f"
	    i__3 = *n - i__;
#line 1560 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 + (
		    i__ + *nb) * xb_dim1], &c__1, (ftnlen)9);
#line 1562 "MB03XU.f"
	    i__2 = *n - i__;
#line 1562 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1562 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 1564 "MB03XU.f"
	    i__2 = *n - i__;
#line 1564 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1564 "MB03XU.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + 
		    b_dim1], ldb, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 
		    + (i__ + *nb) * xb_dim1], &c__1, (ftnlen)12);
#line 1566 "MB03XU.f"
	    i__2 = *k + *n;
#line 1566 "MB03XU.f"
	    d__1 = -taul[i__];
#line 1566 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xb[(i__ + *nb) * xb_dim1 + 1], &c__1);

/*           Update B(:,i). */

#line 1570 "MB03XU.f"
	    i__2 = *k + *n;
#line 1570 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &xb[(i__ + *nb) * xb_dim1 + 1], &c__1, &b[
		    i__ * b_dim1 + 1], &c__1);

#line 1572 "MB03XU.f"
	    a[*k + i__ + i__ * a_dim1] = temp;
#line 1573 "MB03XU.f"
	    q[i__ + i__ * q_dim1] = tauq;
#line 1574 "MB03XU.f"
	    csl[(i__ << 1) - 1] = c__;
#line 1575 "MB03XU.f"
	    csl[i__ * 2] = s;

/*           Transform first rows of Q and B. */

#line 1579 "MB03XU.f"
	    alpha = q[i__ + (i__ + 1) * q_dim1];
#line 1580 "MB03XU.f"
	    i__2 = *n - i__;
#line 1580 "MB03XU.f"
	    dlarfg_(&i__2, &alpha, &q[i__ + (i__ + 2) * q_dim1], ldq, &tauq);
#line 1581 "MB03XU.f"
	    q[i__ + (i__ + 1) * q_dim1] = 1.;
#line 1582 "MB03XU.f"
	    i__2 = *n - i__;
#line 1582 "MB03XU.f"
	    temp = -tauq * ddot_(&i__2, &q[i__ + (i__ + 1) * q_dim1], ldq, &b[
		    *k + i__ + 1 + i__ * b_dim1], &c__1);
#line 1583 "MB03XU.f"
	    i__2 = *n - i__;
#line 1583 "MB03XU.f"
	    daxpy_(&i__2, &temp, &q[i__ + (i__ + 1) * q_dim1], ldq, &b[*k + 
		    i__ + 1 + i__ * b_dim1], &c__1);
#line 1584 "MB03XU.f"
	    temp = b[*k + i__ + 1 + i__ * b_dim1];
#line 1585 "MB03XU.f"
	    dlartg_(&temp, &alpha, &c__, &s, &b[*k + i__ + 1 + i__ * b_dim1]);
#line 1586 "MB03XU.f"
	    s = -s;
#line 1587 "MB03XU.f"
	    i__2 = *n - i__;
#line 1587 "MB03XU.f"
	    dlarfg_(&i__2, &b[*k + i__ + 1 + i__ * b_dim1], &b[*k + i__ + 2 + 
		    i__ * b_dim1], &c__1, &taur[i__]);
#line 1588 "MB03XU.f"
	    temp = b[*k + i__ + 1 + i__ * b_dim1];
#line 1589 "MB03XU.f"
	    b[*k + i__ + 1 + i__ * b_dim1] = 1.;

/*           Update YB with first Householder reflection. */

#line 1593 "MB03XU.f"
	    i__2 = *n - i__;
#line 1593 "MB03XU.f"
	    i__3 = *n - i__;
#line 1593 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + (i__ + 
		    1) * b_dim1], ldb, &q[i__ + (i__ + 1) * q_dim1], ldq, &
		    c_b9, &yb[i__ + 1 + i__ * yb_dim1], &c__1, (ftnlen)9);
#line 1595 "MB03XU.f"
	    i__2 = *n - i__;
#line 1595 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 1 + 
		    xb_dim1], ldxb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yb[i__ * yb_dim1 + 1], &c__1, (ftnlen)9);
#line 1597 "MB03XU.f"
	    i__2 = *n - i__;
#line 1597 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yb[i__ * yb_dim1 + 1], &c__1, &c_b7, &yb[i__ + 1 + 
		    i__ * yb_dim1], &c__1, (ftnlen)12);
#line 1599 "MB03XU.f"
	    i__2 = *n - i__;
#line 1599 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 1 + nb1 * 
		    xb_dim1], ldxb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yb[i__ * yb_dim1 + 1], &c__1, (ftnlen)9);
#line 1601 "MB03XU.f"
	    i__2 = *n - i__;
#line 1601 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &yb[i__ * yb_dim1 + 1], &c__1, &c_b7, &yb[
		    i__ + 1 + i__ * yb_dim1], &c__1, (ftnlen)12);
#line 1603 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1603 "MB03XU.f"
	    i__3 = *n - i__;
#line 1603 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 
		    + 1], ldq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &
		    dwork[1], &c__1, (ftnlen)12);
#line 1605 "MB03XU.f"
	    i__2 = *n - i__;
#line 1605 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1605 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + yb_dim1]
		    , ldyb, &dwork[1], &c__1, &c_b7, &yb[i__ + 1 + i__ * 
		    yb_dim1], &c__1, (ftnlen)12);
#line 1607 "MB03XU.f"
	    i__2 = *n - i__;
#line 1607 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1607 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + b_dim1]
		    , ldb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &dwork[*
		    nb + 1], &c__1, (ftnlen)9);
#line 1609 "MB03XU.f"
	    i__2 = *n - i__;
#line 1609 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1609 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &dwork[*nb + 1], &c__1, &c_b7, &yb[i__ + 
		    1 + i__ * yb_dim1], &c__1, (ftnlen)12);
#line 1611 "MB03XU.f"
	    i__2 = *n - i__;
#line 1611 "MB03XU.f"
	    d__1 = -tauq;
#line 1611 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yb[i__ + 1 + i__ * yb_dim1], &c__1);

/*           Update B(k+i+1,i+1:n). */

#line 1615 "MB03XU.f"
	    i__2 = *n - i__;
#line 1615 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xb[*k + i__ + 1 + xb_dim1], ldxb, &c_b7, &b[*k + 
		    i__ + 1 + (i__ + 1) * b_dim1], ldb, (ftnlen)12);
#line 1617 "MB03XU.f"
	    i__2 = *n - i__;
#line 1617 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xb[*k + i__ + 1 + nb1 * xb_dim1], ldxb, &
		    c_b7, &b[*k + i__ + 1 + (i__ + 1) * b_dim1], ldb, (ftnlen)
		    12);
#line 1619 "MB03XU.f"
	    i__2 = *n - i__;
#line 1619 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yb[i__ + 1 + yb_dim1],
		     ldyb, &q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &b[*k + 
		    i__ + 1 + (i__ + 1) * b_dim1], ldb, (ftnlen)12);
#line 1621 "MB03XU.f"
	    i__2 = *n - i__;
#line 1621 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1621 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &b[*k + i__ + 1 + b_dim1], ldb, &c_b7, &b[
		    *k + i__ + 1 + (i__ + 1) * b_dim1], ldb, (ftnlen)12);

/*           Update YQ with first Householder reflection. */

#line 1626 "MB03XU.f"
	    i__2 = *n - i__;
#line 1626 "MB03XU.f"
	    i__3 = *n - i__;
#line 1626 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[i__ + 1 + (i__ + 1)
		     * q_dim1], ldq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9,
		     &yq[i__ + 1 + i__ * yq_dim1], &c__1, (ftnlen)12);
#line 1628 "MB03XU.f"
	    i__2 = *n - i__;
#line 1628 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1], 
		    ldxq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &yq[i__ * 
		    yq_dim1 + 1], &c__1, (ftnlen)9);
#line 1630 "MB03XU.f"
	    i__2 = *n - i__;
#line 1630 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yq[i__ * yq_dim1 + 1], &c__1, &c_b7, &yq[i__ + 1 + 
		    i__ * yq_dim1], &c__1, (ftnlen)12);
#line 1632 "MB03XU.f"
	    i__2 = *n - i__;
#line 1632 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yq[i__ * yq_dim1 + 1], &c__1, (ftnlen)9);
#line 1634 "MB03XU.f"
	    i__2 = *n - i__;
#line 1634 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &yq[i__ * yq_dim1 + 1], &c__1, &c_b7, &yq[
		    i__ + 1 + i__ * yq_dim1], &c__1, (ftnlen)12);
#line 1636 "MB03XU.f"
	    i__2 = *n - i__;
#line 1636 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1636 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + yq_dim1]
		    , ldyq, &dwork[1], &c__1, &c_b7, &yq[i__ + 1 + i__ * 
		    yq_dim1], &c__1, (ftnlen)12);
#line 1638 "MB03XU.f"
	    i__2 = *n - i__;
#line 1638 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1638 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &dwork[*nb + 1], &c__1, &c_b7, &yq[i__ + 
		    1 + i__ * yq_dim1], &c__1, (ftnlen)12);
#line 1640 "MB03XU.f"
	    i__2 = *n - i__;
#line 1640 "MB03XU.f"
	    d__1 = -tauq;
#line 1640 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yq[i__ + 1 + i__ * yq_dim1], &c__1);

/*           Update Q(i+1:n,i+1). */

#line 1644 "MB03XU.f"
	    i__2 = *n - i__;
#line 1644 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xq[i__ + 1 + xq_dim1], ldxq, &c_b7, &q[i__ + 1 + (
		    i__ + 1) * q_dim1], &c__1, (ftnlen)12);
#line 1646 "MB03XU.f"
	    i__2 = *n - i__;
#line 1646 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xq[i__ + 1 + nb1 * xq_dim1], ldxq, &c_b7, &
		    q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12);
#line 1648 "MB03XU.f"
	    i__2 = *n - i__;
#line 1648 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yq[i__ + 1 + yq_dim1],
		     ldyq, &q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &q[i__ + 
		    1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12);
#line 1650 "MB03XU.f"
	    i__2 = *n - i__;
#line 1650 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1650 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &b[*k + i__ + 1 + b_dim1], ldb, &c_b7, &q[
		    i__ + 1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12);

/*           Apply rotation to [ Q(i+1:n,i+1), B(k+i+1,i+1:n)' ]. */

#line 1655 "MB03XU.f"
	    i__2 = *n - i__;
#line 1655 "MB03XU.f"
	    drot_(&i__2, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, &b[*k + i__ 
		    + 1 + (i__ + 1) * b_dim1], ldb, &c__, &s);
#line 1656 "MB03XU.f"
	    i__2 = i__;
#line 1656 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 1657 "MB03XU.f"
		xb[*k + i__ + 1 + j * xb_dim1] = 0.;
#line 1658 "MB03XU.f"
/* L230: */
#line 1658 "MB03XU.f"
	    }
#line 1659 "MB03XU.f"
	    i__2 = i__;
#line 1659 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 1660 "MB03XU.f"
		xb[*k + i__ + 1 + (*nb + j) * xb_dim1] = 0.;
#line 1661 "MB03XU.f"
/* L240: */
#line 1661 "MB03XU.f"
	    }

/*           Update YB with second Householder reflection. */

#line 1665 "MB03XU.f"
	    i__2 = *n - i__;
#line 1665 "MB03XU.f"
	    i__3 = *n - i__;
#line 1665 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 1 + (i__ + 
		    1) * b_dim1], ldb, &b[*k + i__ + 1 + i__ * b_dim1], &c__1,
		     &c_b9, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (
		    ftnlen)9);
#line 1667 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1667 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 2 + 
		    xb_dim1], ldxb, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, (ftnlen)9);
#line 1669 "MB03XU.f"
	    i__2 = *n - i__;
#line 1669 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, &c_b7, &yb[
		    i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (ftnlen)12);
#line 1671 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1671 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 2 + nb1 * 
		    xb_dim1], ldxb, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, (ftnlen)9);
#line 1673 "MB03XU.f"
	    i__2 = *n - i__;
#line 1673 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, &
		    c_b7, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (
		    ftnlen)12);
#line 1675 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1675 "MB03XU.f"
	    dgemv_("No transpose", &i__, &i__2, &c_b7, &q[(i__ + 2) * q_dim1 
		    + 1], ldq, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &c_b9, 
		    &dwork[nb2 + 1], &c__1, (ftnlen)12);
#line 1677 "MB03XU.f"
	    i__2 = *n - i__;
#line 1677 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yb[i__ + 1 + yb_dim1],
		     ldyb, &dwork[nb2 + 1], &c__1, &c_b7, &yb[i__ + 1 + (i__ 
		    + *nb) * yb_dim1], &c__1, (ftnlen)12);
#line 1679 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1679 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1679 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[*k + i__ + 2 + b_dim1]
		    , ldq, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &c_b9, &
		    dwork[nb3 + 1], &c__1, (ftnlen)9);
#line 1681 "MB03XU.f"
	    i__2 = *n - i__;
#line 1681 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1681 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &dwork[nb3 + 1], &c__1, &c_b7, &yb[i__ + 
		    1 + (i__ + *nb) * yb_dim1], &c__1, (ftnlen)12);
#line 1683 "MB03XU.f"
	    i__2 = *n - i__;
#line 1683 "MB03XU.f"
	    d__1 = -taur[i__];
#line 1683 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1);

/*           Update B(k+i+1,i+1:n). */

#line 1687 "MB03XU.f"
	    i__2 = *n - i__;
#line 1687 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, 
		    &b[*k + i__ + 1 + (i__ + 1) * b_dim1], ldb);

/*           Update YQ with second Householder reflection. */

#line 1691 "MB03XU.f"
	    i__2 = *n - i__;
#line 1691 "MB03XU.f"
	    i__3 = *n - i__;
#line 1691 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[i__ + 1 + (i__ + 1)
		     * q_dim1], ldq, &b[*k + i__ + 1 + i__ * b_dim1], &c__1, &
		    c_b9, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (
		    ftnlen)12);
#line 1693 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1693 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 2 + xq_dim1], 
		    ldxq, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &c_b9, &yq[(
		    i__ + *nb) * yq_dim1 + 1], &c__1, (ftnlen)9);
#line 1695 "MB03XU.f"
	    i__2 = *n - i__;
#line 1695 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, &c_b7, &yq[
		    i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (ftnlen)12);
#line 1697 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1697 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 2 + nb1 * 
		    xq_dim1], ldxq, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, (ftnlen)9);
#line 1699 "MB03XU.f"
	    i__2 = *n - i__;
#line 1699 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, &
		    c_b7, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (
		    ftnlen)12);
#line 1701 "MB03XU.f"
	    i__2 = *n - i__;
#line 1701 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yq[i__ + 1 + yq_dim1],
		     ldyq, &dwork[nb2 + 1], &c__1, &c_b7, &yq[i__ + 1 + (i__ 
		    + *nb) * yq_dim1], &c__1, (ftnlen)12);
#line 1703 "MB03XU.f"
	    i__2 = *n - i__;
#line 1703 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1703 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &dwork[nb3 + 1], &c__1, &c_b7, &yq[i__ + 
		    1 + (i__ + *nb) * yq_dim1], &c__1, (ftnlen)12);
#line 1705 "MB03XU.f"
	    i__2 = *n - i__;
#line 1705 "MB03XU.f"
	    d__1 = -taur[i__];
#line 1705 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1);

/*           Update Q(i+1:n,i+1). */

#line 1709 "MB03XU.f"
	    i__2 = *n - i__;
#line 1709 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, 
		    &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1);

/*           Update YA with first Householder reflection. */

#line 1713 "MB03XU.f"
	    i__2 = *k + *n;
#line 1713 "MB03XU.f"
	    i__3 = *n - i__;
#line 1713 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[(i__ + 1) * a_dim1 
		    + 1], lda, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &ya[
		    i__ * ya_dim1 + 1], &c__1, (ftnlen)12);
#line 1715 "MB03XU.f"
	    i__2 = *n - i__;
#line 1715 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1], 
		    ldxa, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &dwork[
		    pdw], &c__1, (ftnlen)9);
#line 1717 "MB03XU.f"
	    i__2 = *n - i__;
#line 1717 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 + i__ * 
		    ya_dim1], &c__1, (ftnlen)12);
#line 1719 "MB03XU.f"
	    i__2 = *n - i__;
#line 1719 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
#line 1721 "MB03XU.f"
	    i__2 = *n - i__;
#line 1721 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 
		    + i__ * ya_dim1], &c__1, (ftnlen)12);
#line 1723 "MB03XU.f"
	    i__2 = *k + *n;
#line 1723 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1723 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[ya_offset], ldya, 
		    &dwork[1], &c__1, &c_b7, &ya[i__ * ya_dim1 + 1], &c__1, (
		    ftnlen)12);
#line 1725 "MB03XU.f"
	    i__2 = *k + *n;
#line 1725 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1725 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &dwork[*nb + 1], &c__1, &c_b7, &ya[i__ * ya_dim1 
		    + 1], &c__1, (ftnlen)12);
#line 1727 "MB03XU.f"
	    i__2 = *k + *n;
#line 1727 "MB03XU.f"
	    d__1 = -tauq;
#line 1727 "MB03XU.f"
	    dscal_(&i__2, &d__1, &ya[i__ * ya_dim1 + 1], &c__1);

/*           Update A(1:k+n,i+1). */

#line 1731 "MB03XU.f"
	    i__2 = *n - i__;
#line 1731 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xa[i__ + 1 + xa_dim1], ldxa, &c_b7, &a[*k + i__ + 1 
		    + (i__ + 1) * a_dim1], &c__1, (ftnlen)12);
#line 1733 "MB03XU.f"
	    i__2 = *n - i__;
#line 1733 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xa[i__ + 1 + nb1 * xa_dim1], ldxa, &c_b7, &
		    a[*k + i__ + 1 + (i__ + 1) * a_dim1], &c__1, (ftnlen)12);
#line 1735 "MB03XU.f"
	    i__2 = *k + *n;
#line 1735 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &ya[ya_offset], ldya, &
		    q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &a[(i__ + 1) * 
		    a_dim1 + 1], &c__1, (ftnlen)12);
#line 1737 "MB03XU.f"
	    i__2 = *k + *n;
#line 1737 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1737 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &b[*k + i__ + 1 + b_dim1], ldb, &c_b7, &a[(i__ + 
		    1) * a_dim1 + 1], &c__1, (ftnlen)12);

/*           Update YG with first Householder reflection. */

#line 1742 "MB03XU.f"
	    i__2 = *k + *n;
#line 1742 "MB03XU.f"
	    i__3 = *n - i__;
#line 1742 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], ldg, &q[i__ + (i__ + 1) * q_dim1], ldq, &
		    c_b9, &yg[i__ * yg_dim1 + 1], &c__1, (ftnlen)12);
#line 1744 "MB03XU.f"
	    i__2 = *n - i__;
#line 1744 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 1 + 
		    xg_dim1], ldxg, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
#line 1746 "MB03XU.f"
	    i__2 = *n - i__;
#line 1746 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 + i__ * 
		    yg_dim1], &c__1, (ftnlen)12);
#line 1748 "MB03XU.f"
	    i__2 = *n - i__;
#line 1748 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 1 + nb1 * 
		    xg_dim1], ldxg, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
#line 1750 "MB03XU.f"
	    i__2 = *n - i__;
#line 1750 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 
		    + i__ * yg_dim1], &c__1, (ftnlen)12);
#line 1752 "MB03XU.f"
	    i__2 = *k + *n;
#line 1752 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1752 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[yg_offset], ldyg, 
		    &dwork[1], &c__1, &c_b7, &yg[i__ * yg_dim1 + 1], &c__1, (
		    ftnlen)12);
#line 1754 "MB03XU.f"
	    i__2 = *k + *n;
#line 1754 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1754 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &dwork[*nb + 1], &c__1, &c_b7, &yg[i__ * yg_dim1 
		    + 1], &c__1, (ftnlen)12);
#line 1756 "MB03XU.f"
	    i__2 = *k + *n;
#line 1756 "MB03XU.f"
	    d__1 = -tauq;
#line 1756 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yg[i__ * yg_dim1 + 1], &c__1);

/*           Update G(1:k+n,k+i+1). */

#line 1760 "MB03XU.f"
	    i__2 = *n - i__;
#line 1760 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xg[*k + i__ + 1 + xg_dim1], ldxg, &c_b7, &g[*k + 
		    i__ + 1 + (*k + i__ + 1) * g_dim1], &c__1, (ftnlen)12);
#line 1762 "MB03XU.f"
	    i__2 = *n - i__;
#line 1762 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xg[*k + i__ + 1 + nb1 * xg_dim1], ldxg, &
		    c_b7, &g[*k + i__ + 1 + (*k + i__ + 1) * g_dim1], &c__1, (
		    ftnlen)12);
#line 1764 "MB03XU.f"
	    i__2 = *k + *n;
#line 1764 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yg[yg_offset], ldyg, &
		    q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &g[(*k + i__ + 1)
		     * g_dim1 + 1], &c__1, (ftnlen)12);
#line 1766 "MB03XU.f"
	    i__2 = *k + *n;
#line 1766 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1766 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &b[*k + i__ + 1 + b_dim1], ldb, &c_b7, &g[(*k + 
		    i__ + 1) * g_dim1 + 1], &c__1, (ftnlen)12);
#line 1768 "MB03XU.f"
	    i__2 = i__;
#line 1768 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 1769 "MB03XU.f"
		xg[*k + i__ + 1 + j * xg_dim1] = 0.;
#line 1770 "MB03XU.f"
/* L250: */
#line 1770 "MB03XU.f"
	    }
#line 1771 "MB03XU.f"
	    i__2 = i__;
#line 1771 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 1772 "MB03XU.f"
		xg[*k + i__ + 1 + (*nb + j) * xg_dim1] = 0.;
#line 1773 "MB03XU.f"
/* L260: */
#line 1773 "MB03XU.f"
	    }

/*           Apply rotation to [ A(1:k+n,i+1), G(1:k+n,k+i+1) ]. */

#line 1777 "MB03XU.f"
	    i__2 = *k + *n;
#line 1777 "MB03XU.f"
	    drot_(&i__2, &a[(i__ + 1) * a_dim1 + 1], &c__1, &g[(*k + i__ + 1) 
		    * g_dim1 + 1], &c__1, &c__, &s);

/*           Update YA with second Householder reflection. */

#line 1781 "MB03XU.f"
	    i__2 = *k + *n;
#line 1781 "MB03XU.f"
	    i__3 = *n - i__;
#line 1781 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[(i__ + 1) * a_dim1 
		    + 1], lda, &b[*k + i__ + 1 + i__ * b_dim1], &c__1, &c_b9, 
		    &ya[(i__ + *nb) * ya_dim1 + 1], &c__1, (ftnlen)12);
#line 1783 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1783 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 2 + xa_dim1], 
		    ldxa, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
#line 1785 "MB03XU.f"
	    i__2 = *n - i__;
#line 1785 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 + (i__ + 
		    *nb) * ya_dim1], &c__1, (ftnlen)12);
#line 1787 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1787 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 2 + nb1 * 
		    xa_dim1], ldxa, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 1789 "MB03XU.f"
	    i__2 = *n - i__;
#line 1789 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 
		    + (i__ + *nb) * ya_dim1], &c__1, (ftnlen)12);
#line 1791 "MB03XU.f"
	    i__2 = *k + *n;
#line 1791 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &ya[ya_offset], ldya, &
		    dwork[nb2 + 1], &c__1, &c_b7, &ya[(i__ + *nb) * ya_dim1 + 
		    1], &c__1, (ftnlen)12);
#line 1793 "MB03XU.f"
	    i__2 = *k + *n;
#line 1793 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1793 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &dwork[nb3 + 1], &c__1, &c_b7, &ya[(i__ + *nb) * 
		    ya_dim1 + 1], &c__1, (ftnlen)12);
#line 1795 "MB03XU.f"
	    i__2 = *k + *n;
#line 1795 "MB03XU.f"
	    d__1 = -taur[i__];
#line 1795 "MB03XU.f"
	    dscal_(&i__2, &d__1, &ya[(i__ + *nb) * ya_dim1 + 1], &c__1);

/*           Update A(1:k+n,i+1). */

#line 1799 "MB03XU.f"
	    i__2 = *k + *n;
#line 1799 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &ya[(i__ + *nb) * ya_dim1 + 1], &c__1, &a[(
		    i__ + 1) * a_dim1 + 1], &c__1);

/*           Update YG with second Householder reflection. */

#line 1803 "MB03XU.f"
	    i__2 = *k + *n;
#line 1803 "MB03XU.f"
	    i__3 = *n - i__;
#line 1803 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], ldg, &b[*k + i__ + 1 + i__ * b_dim1], &c__1, 
		    &c_b9, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1, (ftnlen)12);
#line 1805 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1805 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 2 + 
		    xg_dim1], ldxg, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 1807 "MB03XU.f"
	    i__2 = *n - i__;
#line 1807 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 + (i__ + 
		    *nb) * yg_dim1], &c__1, (ftnlen)12);
#line 1809 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 1809 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 2 + nb1 * 
		    xg_dim1], ldxg, &b[*k + i__ + 2 + i__ * b_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 1811 "MB03XU.f"
	    i__2 = *n - i__;
#line 1811 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 
		    + (i__ + *nb) * yg_dim1], &c__1, (ftnlen)12);
#line 1813 "MB03XU.f"
	    i__2 = *k + *n;
#line 1813 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yg[yg_offset], ldyg, &
		    dwork[nb2 + 1], &c__1, &c_b7, &yg[(i__ + *nb) * yg_dim1 + 
		    1], &c__1, (ftnlen)12);
#line 1815 "MB03XU.f"
	    i__2 = *k + *n;
#line 1815 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1815 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &dwork[nb3 + 1], &c__1, &c_b7, &yg[(i__ + *nb) * 
		    yg_dim1 + 1], &c__1, (ftnlen)12);
#line 1817 "MB03XU.f"
	    i__2 = *k + *n;
#line 1817 "MB03XU.f"
	    d__1 = -taur[i__];
#line 1817 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1);

/*           Update G(1:k+n,k+i+1). */

#line 1821 "MB03XU.f"
	    i__2 = *k + *n;
#line 1821 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1, &g[(*
		    k + i__ + 1) * g_dim1 + 1], &c__1);

#line 1823 "MB03XU.f"
	    b[*k + i__ + 1 + i__ * b_dim1] = temp;
#line 1824 "MB03XU.f"
	    q[i__ + (i__ + 1) * q_dim1] = tauq;
#line 1825 "MB03XU.f"
	    csr[(i__ << 1) - 1] = c__;
#line 1826 "MB03XU.f"
	    csr[i__ * 2] = s;
#line 1827 "MB03XU.f"
/* L270: */
#line 1827 "MB03XU.f"
	}

#line 1829 "MB03XU.f"
    } else {
#line 1830 "MB03XU.f"
	i__1 = *nb;
#line 1830 "MB03XU.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Transform first columns of A and Q. See routine MB04TS. */

#line 1834 "MB03XU.f"
	    alpha = q[i__ + i__ * q_dim1];
#line 1835 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1835 "MB03XU.f"
	    dlarfg_(&i__2, &alpha, &q[i__ + 1 + i__ * q_dim1], &c__1, &tauq);
#line 1836 "MB03XU.f"
	    q[i__ + i__ * q_dim1] = 1.;
#line 1837 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1837 "MB03XU.f"
	    temp = -tauq * ddot_(&i__2, &q[i__ + i__ * q_dim1], &c__1, &a[*k 
		    + i__ + i__ * a_dim1], &c__1);
#line 1838 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1838 "MB03XU.f"
	    daxpy_(&i__2, &temp, &q[i__ + i__ * q_dim1], &c__1, &a[*k + i__ + 
		    i__ * a_dim1], &c__1);
#line 1839 "MB03XU.f"
	    temp = a[*k + i__ + i__ * a_dim1];
#line 1840 "MB03XU.f"
	    dlartg_(&temp, &alpha, &c__, &s, &a[*k + i__ + i__ * a_dim1]);
#line 1841 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1841 "MB03XU.f"
	    dlarfg_(&i__2, &a[*k + i__ + i__ * a_dim1], &a[*k + i__ + 1 + i__ 
		    * a_dim1], &c__1, &taul[i__]);
#line 1842 "MB03XU.f"
	    temp = a[*k + i__ + i__ * a_dim1];
#line 1843 "MB03XU.f"
	    a[*k + i__ + i__ * a_dim1] = 1.;

/*           Update XQ with first Householder reflection. */

#line 1847 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1847 "MB03XU.f"
	    i__3 = *n - i__;
#line 1847 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[
		    i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)9);
#line 1849 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1849 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1849 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + q_dim1], ldq, &
		    q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[1], &c__1, (
		    ftnlen)9);
#line 1851 "MB03XU.f"
	    i__2 = *n - i__;
#line 1851 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1851 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + xq_dim1]
		    , ldxq, &dwork[1], &c__1, &c_b7, &xq[i__ + 1 + i__ * 
		    xq_dim1], &c__1, (ftnlen)12);
#line 1853 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1853 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1853 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[*k + i__ + a_dim1], 
		    lda, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[*nb + 1]
		    , &c__1, (ftnlen)9);
#line 1855 "MB03XU.f"
	    i__2 = *n - i__;
#line 1855 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1855 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[*nb + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + i__ * xq_dim1], &c__1, (ftnlen)12);
#line 1857 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1857 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1857 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + yq_dim1], ldyq,
		     &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[i__ * xq_dim1 
		    + 1], &c__1, (ftnlen)9);
#line 1859 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1859 "MB03XU.f"
	    i__3 = *n - i__;
#line 1859 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xq[i__ * xq_dim1 + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + i__ * xq_dim1], &c__1, (ftnlen)9);
#line 1861 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1861 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1861 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + nb1 * yq_dim1],
		     ldyq, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xq[(i__ + *
		    nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
#line 1863 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1863 "MB03XU.f"
	    i__3 = *n - i__;
#line 1863 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &
		    c_b7, &xq[i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)9);
#line 1865 "MB03XU.f"
	    i__2 = *n - i__;
#line 1865 "MB03XU.f"
	    d__1 = -tauq;
#line 1865 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + i__ * xq_dim1], &c__1);

/*           Update Q(i,i+1:n). */

#line 1869 "MB03XU.f"
	    i__2 = *n - i__;
#line 1869 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1],
		     ldxq, &q[i__ + q_dim1], ldq, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, (ftnlen)12);
#line 1871 "MB03XU.f"
	    i__2 = *n - i__;
#line 1871 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1871 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &a[*k + i__ + a_dim1], lda, &c_b7, &q[i__ 
		    + (i__ + 1) * q_dim1], ldq, (ftnlen)12);
#line 1873 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1873 "MB03XU.f"
	    i__3 = *n - i__;
#line 1873 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yq[i__ + yq_dim1], ldyq, &c_b7, &q[i__ + (i__ + 
		    1) * q_dim1], ldq, (ftnlen)9);
#line 1875 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1875 "MB03XU.f"
	    i__3 = *n - i__;
#line 1875 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &yq[i__ + nb1 * yq_dim1], ldyq, &c_b7, &
		    q[i__ + (i__ + 1) * q_dim1], ldq, (ftnlen)9);

/*           Update XA with first Householder reflection. */

#line 1880 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1880 "MB03XU.f"
	    i__3 = *n - i__;
#line 1880 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[*k + i__ + (i__ + 1) *
		     a_dim1], lda, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[
		    i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)9);
#line 1882 "MB03XU.f"
	    i__2 = *n - i__;
#line 1882 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1882 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + xa_dim1]
		    , ldxa, &dwork[1], &c__1, &c_b7, &xa[i__ + 1 + i__ * 
		    xa_dim1], &c__1, (ftnlen)12);
#line 1884 "MB03XU.f"
	    i__2 = *n - i__;
#line 1884 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1884 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[*nb + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + i__ * xa_dim1], &c__1, (ftnlen)12);
#line 1886 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1886 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1886 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + ya_dim1], 
		    ldya, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[i__ * 
		    xa_dim1 + 1], &c__1, (ftnlen)9);
#line 1888 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1888 "MB03XU.f"
	    i__3 = *n - i__;
#line 1888 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xa[i__ * xa_dim1 + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + i__ * xa_dim1], &c__1, (ftnlen)9);
#line 1890 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1890 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1890 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + nb1 * 
		    ya_dim1], ldya, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xa[
		    i__ * xa_dim1 + 1], &c__1, (ftnlen)9);
#line 1892 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1892 "MB03XU.f"
	    i__3 = *n - i__;
#line 1892 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &xa[i__ * xa_dim1 + 1], &c__1, &c_b7, &
		    xa[i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)9);
#line 1894 "MB03XU.f"
	    i__2 = *n - i__;
#line 1894 "MB03XU.f"
	    d__1 = -tauq;
#line 1894 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + i__ * xa_dim1], &c__1);

/*           Update A(k+i,i+1:n). */

#line 1898 "MB03XU.f"
	    i__2 = *n - i__;
#line 1898 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1],
		     ldxa, &q[i__ + q_dim1], ldq, &c_b7, &a[*k + i__ + (i__ + 
		    1) * a_dim1], lda, (ftnlen)12);
#line 1900 "MB03XU.f"
	    i__2 = *n - i__;
#line 1900 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1900 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &a[*k + i__ + a_dim1], lda, &c_b7, &a[*k 
		    + i__ + (i__ + 1) * a_dim1], lda, (ftnlen)12);
#line 1902 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1902 "MB03XU.f"
	    i__3 = *n - i__;
#line 1902 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &ya[*k + i__ + ya_dim1], ldya, &c_b7, &a[*k + 
		    i__ + (i__ + 1) * a_dim1], lda, (ftnlen)9);
#line 1904 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1904 "MB03XU.f"
	    i__3 = *n - i__;
#line 1904 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &ya[*k + i__ + nb1 * ya_dim1], ldya, &
		    c_b7, &a[*k + i__ + (i__ + 1) * a_dim1], lda, (ftnlen)9);

/*           Apply rotation to [ A(k+i,i+1:n); Q(i,i+1:n) ]. */

#line 1909 "MB03XU.f"
	    i__2 = *n - i__;
#line 1909 "MB03XU.f"
	    drot_(&i__2, &a[*k + i__ + (i__ + 1) * a_dim1], lda, &q[i__ + (
		    i__ + 1) * q_dim1], ldq, &c__, &s);

/*           Update XQ with second Householder reflection. */

#line 1913 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1913 "MB03XU.f"
	    i__3 = *n - i__;
#line 1913 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[i__ + (i__ + 1) * 
		    q_dim1], ldq, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b9, &
		    xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)9);
#line 1915 "MB03XU.f"
	    i__2 = *n - i__;
#line 1915 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], ldq,
		     &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &c_b9, &dwork[
		    nb2 + 1], &c__1, (ftnlen)9);
#line 1917 "MB03XU.f"
	    i__2 = *n - i__;
#line 1917 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1],
		     ldxq, &dwork[nb2 + 1], &c__1, &c_b7, &xq[i__ + 1 + (i__ 
		    + *nb) * xq_dim1], &c__1, (ftnlen)12);
#line 1919 "MB03XU.f"
	    i__2 = *n - i__;
#line 1919 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1919 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[*k + i__ + 1 + a_dim1]
		    , lda, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &c_b9, &
		    dwork[nb3 + 1], &c__1, (ftnlen)9);
#line 1921 "MB03XU.f"
	    i__2 = *n - i__;
#line 1921 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1921 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[nb3 + 1], &c__1, &c_b7, &xq[i__ + 
		    1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)12);
#line 1923 "MB03XU.f"
	    i__2 = *n - i__;
#line 1923 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1923 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + yq_dim1], 
		    ldyq, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &c_b9, &xq[(
		    i__ + *nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
#line 1925 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1925 "MB03XU.f"
	    i__3 = *n - i__;
#line 1925 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &c_b7, &
		    xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (ftnlen)9);
#line 1927 "MB03XU.f"
	    i__2 = *n - i__;
#line 1927 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1927 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, (ftnlen)9);
#line 1929 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1929 "MB03XU.f"
	    i__3 = *n - i__;
#line 1929 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &xq[(i__ + *nb) * xq_dim1 + 1], &c__1, &
		    c_b7, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, (
		    ftnlen)9);
#line 1931 "MB03XU.f"
	    i__2 = *n - i__;
#line 1931 "MB03XU.f"
	    d__1 = -taul[i__];
#line 1931 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1);

/*           Update Q(i,i+1:n). */

#line 1935 "MB03XU.f"
	    i__2 = *n - i__;
#line 1935 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &xq[i__ + 1 + (i__ + *nb) * xq_dim1], &c__1, 
		    &q[i__ + (i__ + 1) * q_dim1], ldq);

/*           Update XA with second Householder reflection. */

#line 1939 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1939 "MB03XU.f"
	    i__3 = *n - i__;
#line 1939 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &a[*k + i__ + (i__ + 1) *
		     a_dim1], lda, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b9, 
		    &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)9);
#line 1941 "MB03XU.f"
	    i__2 = *n - i__;
#line 1941 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1],
		     ldxa, &dwork[nb2 + 1], &c__1, &c_b7, &xa[i__ + 1 + (i__ 
		    + *nb) * xa_dim1], &c__1, (ftnlen)12);
#line 1943 "MB03XU.f"
	    i__2 = *n - i__;
#line 1943 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1943 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[nb3 + 1], &c__1, &c_b7, &xa[i__ + 
		    1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)12);
#line 1945 "MB03XU.f"
	    i__2 = *n - i__;
#line 1945 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1945 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + 1 + 
		    ya_dim1], ldya, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, (ftnlen)9);
#line 1947 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1947 "MB03XU.f"
	    i__3 = *n - i__;
#line 1947 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, &c_b7, &
		    xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (ftnlen)9);
#line 1949 "MB03XU.f"
	    i__2 = *n - i__;
#line 1949 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1949 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &ya[*k + i__ + 1 + nb1 * 
		    ya_dim1], ldya, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, (ftnlen)9);
#line 1951 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1951 "MB03XU.f"
	    i__3 = *n - i__;
#line 1951 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &xa[(i__ + *nb) * xa_dim1 + 1], &c__1, &
		    c_b7, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, (
		    ftnlen)9);
#line 1953 "MB03XU.f"
	    i__2 = *n - i__;
#line 1953 "MB03XU.f"
	    d__1 = -taul[i__];
#line 1953 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1);

/*           Update A(k+i,i+1:n). */

#line 1957 "MB03XU.f"
	    i__2 = *n - i__;
#line 1957 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &xa[i__ + 1 + (i__ + *nb) * xa_dim1], &c__1, 
		    &a[*k + i__ + (i__ + 1) * a_dim1], lda);

/*           Update XG with first Householder reflection. */

#line 1961 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1961 "MB03XU.f"
	    i__3 = *k + *n;
#line 1961 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &g[*k + i__ + g_dim1], 
		    ldg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &xg[i__ * 
		    xg_dim1 + 1], &c__1, (ftnlen)9);
#line 1963 "MB03XU.f"
	    i__2 = *k + *n;
#line 1963 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1963 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[xg_offset], ldxg, 
		    &dwork[1], &c__1, &c_b7, &xg[i__ * xg_dim1 + 1], &c__1, (
		    ftnlen)12);
#line 1965 "MB03XU.f"
	    i__2 = *k + *n;
#line 1965 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1965 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &dwork[*nb + 1], &c__1, &c_b7, &xg[i__ * xg_dim1 
		    + 1], &c__1, (ftnlen)12);
#line 1967 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1967 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1967 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + yg_dim1], 
		    ldyg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], &
		    c__1, (ftnlen)9);
#line 1969 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1969 "MB03XU.f"
	    i__3 = *n - i__;
#line 1969 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 + 
		    i__ * xg_dim1], &c__1, (ftnlen)9);
#line 1971 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1971 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1971 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + nb1 * 
		    yg_dim1], ldyg, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
#line 1973 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1973 "MB03XU.f"
	    i__3 = *n - i__;
#line 1973 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ 
		    + 1 + i__ * xg_dim1], &c__1, (ftnlen)9);
#line 1975 "MB03XU.f"
	    i__2 = *k + *n;
#line 1975 "MB03XU.f"
	    d__1 = -tauq;
#line 1975 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xg[i__ * xg_dim1 + 1], &c__1);

/*           Update G(k+i,:). */

#line 1979 "MB03XU.f"
	    i__2 = *k + *n;
#line 1979 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xg[xg_offset], ldxg, &
		    q[i__ + q_dim1], ldq, &c_b7, &g[*k + i__ + g_dim1], ldg, (
		    ftnlen)12);
#line 1981 "MB03XU.f"
	    i__2 = *k + *n;
#line 1981 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1981 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &a[*k + i__ + a_dim1], lda, &c_b7, &g[*k + i__ + 
		    g_dim1], ldg, (ftnlen)12);
#line 1983 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1983 "MB03XU.f"
	    i__3 = *n - i__;
#line 1983 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yg[*k + i__ + yg_dim1], ldyg, &c_b7, &g[*k + 
		    i__ + (*k + i__ + 1) * g_dim1], ldg, (ftnlen)9);
#line 1985 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1985 "MB03XU.f"
	    i__3 = *n - i__;
#line 1985 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &yg[*k + i__ + nb1 * yg_dim1], ldyg, &
		    c_b7, &g[*k + i__ + (*k + i__ + 1) * g_dim1], ldg, (
		    ftnlen)9);

/*           Update XB with first Householder reflection. */

#line 1990 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1990 "MB03XU.f"
	    i__3 = *k + *n;
#line 1990 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[i__ + b_dim1], ldb, &
		    q[i__ + i__ * q_dim1], &c__1, &c_b9, &xb[i__ * xb_dim1 + 
		    1], &c__1, (ftnlen)9);
#line 1992 "MB03XU.f"
	    i__2 = *k + *n;
#line 1992 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1992 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[xb_offset], ldxb, 
		    &dwork[1], &c__1, &c_b7, &xb[i__ * xb_dim1 + 1], &c__1, (
		    ftnlen)12);
#line 1994 "MB03XU.f"
	    i__2 = *k + *n;
#line 1994 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1994 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &dwork[*nb + 1], &c__1, &c_b7, &xb[i__ * xb_dim1 
		    + 1], &c__1, (ftnlen)12);
#line 1996 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 1996 "MB03XU.f"
	    i__3 = i__ - 1;
#line 1996 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + yb_dim1], ldyb,
		     &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], &c__1,
		     (ftnlen)9);
#line 1998 "MB03XU.f"
	    i__2 = i__ - 1;
#line 1998 "MB03XU.f"
	    i__3 = *n - i__;
#line 1998 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 + 
		    i__ * xb_dim1], &c__1, (ftnlen)9);
#line 2000 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 2000 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2000 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + nb1 * yb_dim1],
		     ldyb, &q[i__ + i__ * q_dim1], &c__1, &c_b9, &dwork[pdw], 
		    &c__1, (ftnlen)9);
#line 2002 "MB03XU.f"
	    i__2 = i__ - 1;
#line 2002 "MB03XU.f"
	    i__3 = *n - i__;
#line 2002 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ 
		    + 1 + i__ * xb_dim1], &c__1, (ftnlen)9);
#line 2004 "MB03XU.f"
	    i__2 = *k + *n;
#line 2004 "MB03XU.f"
	    d__1 = -tauq;
#line 2004 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xb[i__ * xb_dim1 + 1], &c__1);

/*           Update B(i,:). */

#line 2008 "MB03XU.f"
	    i__2 = *k + *n;
#line 2008 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xb[xb_offset], ldxb, &
		    q[i__ + q_dim1], ldq, &c_b7, &b[i__ + b_dim1], ldb, (
		    ftnlen)12);
#line 2010 "MB03XU.f"
	    i__2 = *k + *n;
#line 2010 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2010 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &a[*k + i__ + a_dim1], lda, &c_b7, &b[i__ + 
		    b_dim1], ldb, (ftnlen)12);
#line 2012 "MB03XU.f"
	    i__2 = i__ - 1;
#line 2012 "MB03XU.f"
	    i__3 = *n - i__;
#line 2012 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &yb[i__ + yb_dim1], ldyb, &c_b7, &b[i__ + (*k + 
		    i__ + 1) * b_dim1], ldb, (ftnlen)9);
#line 2014 "MB03XU.f"
	    i__2 = i__ - 1;
#line 2014 "MB03XU.f"
	    i__3 = *n - i__;
#line 2014 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &yb[i__ + nb1 * yb_dim1], ldyb, &c_b7, &
		    b[i__ + (*k + i__ + 1) * b_dim1], ldb, (ftnlen)9);

/*           Apply rotation to [ G(k+i,:); B(i,:) ]. */

#line 2019 "MB03XU.f"
	    i__2 = *k + *n;
#line 2019 "MB03XU.f"
	    drot_(&i__2, &g[*k + i__ + g_dim1], ldg, &b[i__ + b_dim1], ldb, &
		    c__, &s);

#line 2021 "MB03XU.f"
	    i__2 = i__ - 1;
#line 2021 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 2022 "MB03XU.f"
		yg[*k + i__ + j * yg_dim1] = 0.;
#line 2023 "MB03XU.f"
/* L280: */
#line 2023 "MB03XU.f"
	    }
#line 2024 "MB03XU.f"
	    i__2 = i__ - 1;
#line 2024 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 2025 "MB03XU.f"
		yg[*k + i__ + (*nb + j) * yg_dim1] = 0.;
#line 2026 "MB03XU.f"
/* L290: */
#line 2026 "MB03XU.f"
	    }
#line 2027 "MB03XU.f"
	    i__2 = i__ - 1;
#line 2027 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 2028 "MB03XU.f"
		ya[*k + i__ + j * ya_dim1] = 0.;
#line 2029 "MB03XU.f"
/* L300: */
#line 2029 "MB03XU.f"
	    }
#line 2030 "MB03XU.f"
	    i__2 = i__ - 1;
#line 2030 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 2031 "MB03XU.f"
		ya[*k + i__ + (*nb + j) * ya_dim1] = 0.;
#line 2032 "MB03XU.f"
/* L310: */
#line 2032 "MB03XU.f"
	    }

/*           Update XG with second Householder reflection. */

#line 2036 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 2036 "MB03XU.f"
	    i__3 = *k + *n;
#line 2036 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &g[*k + i__ + g_dim1], 
		    ldg, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b9, &xg[(i__ 
		    + *nb) * xg_dim1 + 1], &c__1, (ftnlen)9);
#line 2038 "MB03XU.f"
	    i__2 = *k + *n;
#line 2038 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xg[xg_offset], ldxg, &
		    dwork[nb2 + 1], &c__1, &c_b7, &xg[(i__ + *nb) * xg_dim1 + 
		    1], &c__1, (ftnlen)12);
#line 2040 "MB03XU.f"
	    i__2 = *k + *n;
#line 2040 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2040 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xg[nb1 * xg_dim1 + 1]
		    , ldxg, &dwork[nb3 + 1], &c__1, &c_b7, &xg[(i__ + *nb) * 
		    xg_dim1 + 1], &c__1, (ftnlen)12);
#line 2042 "MB03XU.f"
	    i__2 = *n - i__;
#line 2042 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2042 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + 1 + 
		    yg_dim1], ldyg, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 2044 "MB03XU.f"
	    i__2 = i__ - 1;
#line 2044 "MB03XU.f"
	    i__3 = *n - i__;
#line 2044 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ + 1 + (
		    i__ + *nb) * xg_dim1], &c__1, (ftnlen)9);
#line 2046 "MB03XU.f"
	    i__2 = *n - i__;
#line 2046 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2046 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yg[*k + i__ + 1 + nb1 * 
		    yg_dim1], ldyg, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 2048 "MB03XU.f"
	    i__2 = i__ - 1;
#line 2048 "MB03XU.f"
	    i__3 = *n - i__;
#line 2048 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &dwork[pdw], &c__1, &c_b7, &xg[*k + i__ 
		    + 1 + (i__ + *nb) * xg_dim1], &c__1, (ftnlen)9);
#line 2050 "MB03XU.f"
	    i__2 = *k + *n;
#line 2050 "MB03XU.f"
	    d__1 = -taul[i__];
#line 2050 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xg[(i__ + *nb) * xg_dim1 + 1], &c__1);

/*           Update G(k+i,:). */

#line 2054 "MB03XU.f"
	    i__2 = *k + *n;
#line 2054 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &xg[(i__ + *nb) * xg_dim1 + 1], &c__1, &g[*k 
		    + i__ + g_dim1], ldg);

/*           Update XB with second Householder reflection. */

#line 2058 "MB03XU.f"
	    i__2 = *n - i__ + 1;
#line 2058 "MB03XU.f"
	    i__3 = *k + *n;
#line 2058 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[i__ + b_dim1], ldb, &
		    a[*k + i__ + i__ * a_dim1], &c__1, &c_b9, &xb[(i__ + *nb) 
		    * xb_dim1 + 1], &c__1, (ftnlen)9);
#line 2060 "MB03XU.f"
	    i__2 = *k + *n;
#line 2060 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &xb[xb_offset], ldxb, &
		    dwork[nb2 + 1], &c__1, &c_b7, &xb[(i__ + *nb) * xb_dim1 + 
		    1], &c__1, (ftnlen)12);
#line 2062 "MB03XU.f"
	    i__2 = *k + *n;
#line 2062 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2062 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &xb[nb1 * xb_dim1 + 1]
		    , ldxb, &dwork[nb3 + 1], &c__1, &c_b7, &xb[(i__ + *nb) * 
		    xb_dim1 + 1], &c__1, (ftnlen)12);
#line 2064 "MB03XU.f"
	    i__2 = *n - i__;
#line 2064 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2064 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + yb_dim1], 
		    ldyb, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
#line 2066 "MB03XU.f"
	    i__2 = i__ - 1;
#line 2066 "MB03XU.f"
	    i__3 = *n - i__;
#line 2066 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 + 
		    1], ldq, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ + 1 + (
		    i__ + *nb) * xb_dim1], &c__1, (ftnlen)9);
#line 2068 "MB03XU.f"
	    i__2 = *n - i__;
#line 2068 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2068 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 2070 "MB03XU.f"
	    i__2 = i__ - 1;
#line 2070 "MB03XU.f"
	    i__3 = *n - i__;
#line 2070 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &dwork[pdw], &c__1, &c_b7, &xb[*k + i__ 
		    + 1 + (i__ + *nb) * xb_dim1], &c__1, (ftnlen)9);
#line 2072 "MB03XU.f"
	    i__2 = *k + *n;
#line 2072 "MB03XU.f"
	    d__1 = -taul[i__];
#line 2072 "MB03XU.f"
	    dscal_(&i__2, &d__1, &xb[(i__ + *nb) * xb_dim1 + 1], &c__1);

/*           Update B(i,:). */

#line 2076 "MB03XU.f"
	    i__2 = *k + *n;
#line 2076 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &xb[(i__ + *nb) * xb_dim1 + 1], &c__1, &b[
		    i__ + b_dim1], ldb);

#line 2078 "MB03XU.f"
	    a[*k + i__ + i__ * a_dim1] = temp;
#line 2079 "MB03XU.f"
	    q[i__ + i__ * q_dim1] = tauq;
#line 2080 "MB03XU.f"
	    csl[(i__ << 1) - 1] = c__;
#line 2081 "MB03XU.f"
	    csl[i__ * 2] = s;

/*           Transform first rows of Q and B. */

#line 2085 "MB03XU.f"
	    alpha = q[i__ + (i__ + 1) * q_dim1];
#line 2086 "MB03XU.f"
	    i__2 = *n - i__;
#line 2086 "MB03XU.f"
	    dlarfg_(&i__2, &alpha, &q[i__ + (i__ + 2) * q_dim1], ldq, &tauq);
#line 2087 "MB03XU.f"
	    q[i__ + (i__ + 1) * q_dim1] = 1.;
#line 2088 "MB03XU.f"
	    i__2 = *n - i__;
#line 2088 "MB03XU.f"
	    temp = -tauq * ddot_(&i__2, &q[i__ + (i__ + 1) * q_dim1], ldq, &b[
		    i__ + (*k + i__ + 1) * b_dim1], ldb);
#line 2089 "MB03XU.f"
	    i__2 = *n - i__;
#line 2089 "MB03XU.f"
	    daxpy_(&i__2, &temp, &q[i__ + (i__ + 1) * q_dim1], ldq, &b[i__ + (
		    *k + i__ + 1) * b_dim1], ldb);
#line 2090 "MB03XU.f"
	    temp = b[i__ + (*k + i__ + 1) * b_dim1];
#line 2091 "MB03XU.f"
	    dlartg_(&temp, &alpha, &c__, &s, &b[i__ + (*k + i__ + 1) * b_dim1]
		    );
#line 2092 "MB03XU.f"
	    s = -s;
#line 2093 "MB03XU.f"
	    i__2 = *n - i__;
#line 2093 "MB03XU.f"
	    dlarfg_(&i__2, &b[i__ + (*k + i__ + 1) * b_dim1], &b[i__ + (*k + 
		    i__ + 2) * b_dim1], ldb, &taur[i__]);
#line 2094 "MB03XU.f"
	    temp = b[i__ + (*k + i__ + 1) * b_dim1];
#line 2095 "MB03XU.f"
	    b[i__ + (*k + i__ + 1) * b_dim1] = 1.;

/*           Update YB with first Householder reflection. */

#line 2099 "MB03XU.f"
	    i__2 = *n - i__;
#line 2099 "MB03XU.f"
	    i__3 = *n - i__;
#line 2099 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[i__ + 1 + (*k + 
		    i__ + 1) * b_dim1], ldb, &q[i__ + (i__ + 1) * q_dim1], 
		    ldq, &c_b9, &yb[i__ + 1 + i__ * yb_dim1], &c__1, (ftnlen)
		    12);
#line 2101 "MB03XU.f"
	    i__2 = *n - i__;
#line 2101 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 1 + 
		    xb_dim1], ldxb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yb[i__ * yb_dim1 + 1], &c__1, (ftnlen)9);
#line 2103 "MB03XU.f"
	    i__2 = *n - i__;
#line 2103 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yb[i__ * yb_dim1 + 1], &c__1, &c_b7, &yb[i__ + 1 + 
		    i__ * yb_dim1], &c__1, (ftnlen)12);
#line 2105 "MB03XU.f"
	    i__2 = *n - i__;
#line 2105 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 1 + nb1 * 
		    xb_dim1], ldxb, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yb[i__ * yb_dim1 + 1], &c__1, (ftnlen)9);
#line 2107 "MB03XU.f"
	    i__2 = *n - i__;
#line 2107 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &yb[i__ * yb_dim1 + 1], &c__1, &c_b7, &yb[
		    i__ + 1 + i__ * yb_dim1], &c__1, (ftnlen)12);
#line 2109 "MB03XU.f"
	    i__2 = i__ - 1;
#line 2109 "MB03XU.f"
	    i__3 = *n - i__;
#line 2109 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[(i__ + 1) * q_dim1 
		    + 1], ldq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &
		    dwork[1], &c__1, (ftnlen)12);
#line 2111 "MB03XU.f"
	    i__2 = *n - i__;
#line 2111 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2111 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + yb_dim1]
		    , ldyb, &dwork[1], &c__1, &c_b7, &yb[i__ + 1 + i__ * 
		    yb_dim1], &c__1, (ftnlen)12);
#line 2113 "MB03XU.f"
	    i__2 = i__ - 1;
#line 2113 "MB03XU.f"
	    i__3 = *n - i__;
#line 2113 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 1) * 
		    b_dim1 + 1], ldb, &q[i__ + (i__ + 1) * q_dim1], ldq, &
		    c_b9, &dwork[*nb + 1], &c__1, (ftnlen)12);
#line 2115 "MB03XU.f"
	    i__2 = *n - i__;
#line 2115 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2115 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &dwork[*nb + 1], &c__1, &c_b7, &yb[i__ + 
		    1 + i__ * yb_dim1], &c__1, (ftnlen)12);
#line 2117 "MB03XU.f"
	    i__2 = *n - i__;
#line 2117 "MB03XU.f"
	    d__1 = -tauq;
#line 2117 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yb[i__ + 1 + i__ * yb_dim1], &c__1);

/*           Update B(i+1:n,k+i+1). */

#line 2121 "MB03XU.f"
	    i__2 = *n - i__;
#line 2121 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xb[*k + i__ + 1 + xb_dim1], ldxb, &c_b7, &b[i__ + 1 
		    + (*k + i__ + 1) * b_dim1], &c__1, (ftnlen)12);
#line 2123 "MB03XU.f"
	    i__2 = *n - i__;
#line 2123 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xb[*k + i__ + 1 + nb1 * xb_dim1], ldxb, &
		    c_b7, &b[i__ + 1 + (*k + i__ + 1) * b_dim1], &c__1, (
		    ftnlen)12);
#line 2125 "MB03XU.f"
	    i__2 = *n - i__;
#line 2125 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yb[i__ + 1 + yb_dim1],
		     ldyb, &q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &b[i__ + 
		    1 + (*k + i__ + 1) * b_dim1], &c__1, (ftnlen)12);
#line 2127 "MB03XU.f"
	    i__2 = *n - i__;
#line 2127 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2127 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &b[(*k + i__ + 1) * b_dim1 + 1], &c__1, &
		    c_b7, &b[i__ + 1 + (*k + i__ + 1) * b_dim1], &c__1, (
		    ftnlen)12);

/*           Update YQ with first Householder reflection. */

#line 2132 "MB03XU.f"
	    i__2 = *n - i__;
#line 2132 "MB03XU.f"
	    i__3 = *n - i__;
#line 2132 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[i__ + 1 + (i__ + 1)
		     * q_dim1], ldq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9,
		     &yq[i__ + 1 + i__ * yq_dim1], &c__1, (ftnlen)12);
#line 2134 "MB03XU.f"
	    i__2 = *n - i__;
#line 2134 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + xq_dim1], 
		    ldxq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &yq[i__ * 
		    yq_dim1 + 1], &c__1, (ftnlen)9);
#line 2136 "MB03XU.f"
	    i__2 = *n - i__;
#line 2136 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yq[i__ * yq_dim1 + 1], &c__1, &c_b7, &yq[i__ + 1 + 
		    i__ * yq_dim1], &c__1, (ftnlen)12);
#line 2138 "MB03XU.f"
	    i__2 = *n - i__;
#line 2138 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &yq[i__ * yq_dim1 + 1], &c__1, (ftnlen)9);
#line 2140 "MB03XU.f"
	    i__2 = *n - i__;
#line 2140 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &yq[i__ * yq_dim1 + 1], &c__1, &c_b7, &yq[
		    i__ + 1 + i__ * yq_dim1], &c__1, (ftnlen)12);
#line 2142 "MB03XU.f"
	    i__2 = *n - i__;
#line 2142 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2142 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + yq_dim1]
		    , ldyq, &dwork[1], &c__1, &c_b7, &yq[i__ + 1 + i__ * 
		    yq_dim1], &c__1, (ftnlen)12);
#line 2144 "MB03XU.f"
	    i__2 = *n - i__;
#line 2144 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2144 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &dwork[*nb + 1], &c__1, &c_b7, &yq[i__ + 
		    1 + i__ * yq_dim1], &c__1, (ftnlen)12);
#line 2146 "MB03XU.f"
	    i__2 = *n - i__;
#line 2146 "MB03XU.f"
	    d__1 = -tauq;
#line 2146 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yq[i__ + 1 + i__ * yq_dim1], &c__1);

/*           Update Q(i+1:n,i+1). */

#line 2150 "MB03XU.f"
	    i__2 = *n - i__;
#line 2150 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xq[i__ + 1 + xq_dim1], ldxq, &c_b7, &q[i__ + 1 + (
		    i__ + 1) * q_dim1], &c__1, (ftnlen)12);
#line 2152 "MB03XU.f"
	    i__2 = *n - i__;
#line 2152 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xq[i__ + 1 + nb1 * xq_dim1], ldxq, &c_b7, &
		    q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12);
#line 2154 "MB03XU.f"
	    i__2 = *n - i__;
#line 2154 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yq[i__ + 1 + yq_dim1],
		     ldyq, &q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &q[i__ + 
		    1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12);
#line 2156 "MB03XU.f"
	    i__2 = *n - i__;
#line 2156 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2156 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &b[(*k + i__ + 1) * b_dim1 + 1], &c__1, &
		    c_b7, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, (ftnlen)12)
		    ;

/*           Apply rotation to [ Q(i+1:n,i+1), B(i+1:n,k+i+1) ]. */

#line 2161 "MB03XU.f"
	    i__2 = *n - i__;
#line 2161 "MB03XU.f"
	    drot_(&i__2, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, &b[i__ + 1 
		    + (*k + i__ + 1) * b_dim1], &c__1, &c__, &s);
#line 2162 "MB03XU.f"
	    i__2 = i__;
#line 2162 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 2163 "MB03XU.f"
		xb[*k + i__ + 1 + j * xb_dim1] = 0.;
#line 2164 "MB03XU.f"
/* L320: */
#line 2164 "MB03XU.f"
	    }
#line 2165 "MB03XU.f"
	    i__2 = i__;
#line 2165 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 2166 "MB03XU.f"
		xb[*k + i__ + 1 + (*nb + j) * xb_dim1] = 0.;
#line 2167 "MB03XU.f"
/* L330: */
#line 2167 "MB03XU.f"
	    }

/*           Update YB with second Householder reflection. */

#line 2171 "MB03XU.f"
	    i__2 = *n - i__;
#line 2171 "MB03XU.f"
	    i__3 = *n - i__;
#line 2171 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[i__ + 1 + (*k + 
		    i__ + 1) * b_dim1], ldb, &b[i__ + (*k + i__ + 1) * b_dim1]
		    , ldb, &c_b9, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1,
		     (ftnlen)12);
#line 2173 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 2173 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 2 + 
		    xb_dim1], ldxb, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, (ftnlen)9);
#line 2175 "MB03XU.f"
	    i__2 = *n - i__;
#line 2175 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, &c_b7, &yb[
		    i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (ftnlen)12);
#line 2177 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 2177 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xb[*k + i__ + 2 + nb1 * 
		    xb_dim1], ldxb, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, (ftnlen)9);
#line 2179 "MB03XU.f"
	    i__2 = *n - i__;
#line 2179 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &yb[(i__ + *nb) * yb_dim1 + 1], &c__1, &
		    c_b7, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, (
		    ftnlen)12);
#line 2181 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 2181 "MB03XU.f"
	    dgemv_("No transpose", &i__, &i__2, &c_b7, &q[(i__ + 2) * q_dim1 
		    + 1], ldq, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &c_b9, 
		    &dwork[nb2 + 1], &c__1, (ftnlen)12);
#line 2183 "MB03XU.f"
	    i__2 = *n - i__;
#line 2183 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yb[i__ + 1 + yb_dim1],
		     ldyb, &dwork[nb2 + 1], &c__1, &c_b7, &yb[i__ + 1 + (i__ 
		    + *nb) * yb_dim1], &c__1, (ftnlen)12);
#line 2185 "MB03XU.f"
	    i__2 = i__ - 1;
#line 2185 "MB03XU.f"
	    i__3 = *n - i__ - 1;
#line 2185 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &b[(*k + i__ + 2) * 
		    b_dim1 + 1], ldq, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, 
		    &c_b9, &dwork[nb3 + 1], &c__1, (ftnlen)12);
#line 2187 "MB03XU.f"
	    i__2 = *n - i__;
#line 2187 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2187 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yb[i__ + 1 + nb1 * 
		    yb_dim1], ldyb, &dwork[nb3 + 1], &c__1, &c_b7, &yb[i__ + 
		    1 + (i__ + *nb) * yb_dim1], &c__1, (ftnlen)12);
#line 2189 "MB03XU.f"
	    i__2 = *n - i__;
#line 2189 "MB03XU.f"
	    d__1 = -taur[i__];
#line 2189 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1);

/*           Update B(i+1:n,k+i+1). */

#line 2193 "MB03XU.f"
	    i__2 = *n - i__;
#line 2193 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &yb[i__ + 1 + (i__ + *nb) * yb_dim1], &c__1, 
		    &b[i__ + 1 + (*k + i__ + 1) * b_dim1], &c__1);

/*           Update YQ with second Householder reflection. */

#line 2197 "MB03XU.f"
	    i__2 = *n - i__;
#line 2197 "MB03XU.f"
	    i__3 = *n - i__;
#line 2197 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &q[i__ + 1 + (i__ + 1)
		     * q_dim1], ldq, &b[i__ + (*k + i__ + 1) * b_dim1], ldb, &
		    c_b9, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (
		    ftnlen)12);
#line 2199 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 2199 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 2 + xq_dim1], 
		    ldxq, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &c_b9, &yq[(
		    i__ + *nb) * yq_dim1 + 1], &c__1, (ftnlen)9);
#line 2201 "MB03XU.f"
	    i__2 = *n - i__;
#line 2201 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, &c_b7, &yq[
		    i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (ftnlen)12);
#line 2203 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 2203 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xq[i__ + 2 + nb1 * 
		    xq_dim1], ldxq, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, (ftnlen)9);
#line 2205 "MB03XU.f"
	    i__2 = *n - i__;
#line 2205 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &yq[(i__ + *nb) * yq_dim1 + 1], &c__1, &
		    c_b7, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, (
		    ftnlen)12);
#line 2207 "MB03XU.f"
	    i__2 = *n - i__;
#line 2207 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yq[i__ + 1 + yq_dim1],
		     ldyq, &dwork[nb2 + 1], &c__1, &c_b7, &yq[i__ + 1 + (i__ 
		    + *nb) * yq_dim1], &c__1, (ftnlen)12);
#line 2209 "MB03XU.f"
	    i__2 = *n - i__;
#line 2209 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2209 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yq[i__ + 1 + nb1 * 
		    yq_dim1], ldyq, &dwork[nb3 + 1], &c__1, &c_b7, &yq[i__ + 
		    1 + (i__ + *nb) * yq_dim1], &c__1, (ftnlen)12);
#line 2211 "MB03XU.f"
	    i__2 = *n - i__;
#line 2211 "MB03XU.f"
	    d__1 = -taur[i__];
#line 2211 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1);

/*           Update Q(i+1:n,i+1). */

#line 2215 "MB03XU.f"
	    i__2 = *n - i__;
#line 2215 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &yq[i__ + 1 + (i__ + *nb) * yq_dim1], &c__1, 
		    &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1);

/*           Update YA with first Householder reflection. */

#line 2219 "MB03XU.f"
	    i__2 = *k + *n;
#line 2219 "MB03XU.f"
	    i__3 = *n - i__;
#line 2219 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[(i__ + 1) * a_dim1 
		    + 1], lda, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &ya[
		    i__ * ya_dim1 + 1], &c__1, (ftnlen)12);
#line 2221 "MB03XU.f"
	    i__2 = *n - i__;
#line 2221 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + xa_dim1], 
		    ldxa, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, &dwork[
		    pdw], &c__1, (ftnlen)9);
#line 2223 "MB03XU.f"
	    i__2 = *n - i__;
#line 2223 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 + i__ * 
		    ya_dim1], &c__1, (ftnlen)12);
#line 2225 "MB03XU.f"
	    i__2 = *n - i__;
#line 2225 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
#line 2227 "MB03XU.f"
	    i__2 = *n - i__;
#line 2227 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 
		    + i__ * ya_dim1], &c__1, (ftnlen)12);
#line 2229 "MB03XU.f"
	    i__2 = *k + *n;
#line 2229 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2229 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[ya_offset], ldya, 
		    &dwork[1], &c__1, &c_b7, &ya[i__ * ya_dim1 + 1], &c__1, (
		    ftnlen)12);
#line 2231 "MB03XU.f"
	    i__2 = *k + *n;
#line 2231 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2231 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &dwork[*nb + 1], &c__1, &c_b7, &ya[i__ * ya_dim1 
		    + 1], &c__1, (ftnlen)12);
#line 2233 "MB03XU.f"
	    i__2 = *k + *n;
#line 2233 "MB03XU.f"
	    d__1 = -tauq;
#line 2233 "MB03XU.f"
	    dscal_(&i__2, &d__1, &ya[i__ * ya_dim1 + 1], &c__1);

/*           Update A(1:k+n,i+1). */

#line 2237 "MB03XU.f"
	    i__2 = *n - i__;
#line 2237 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xa[i__ + 1 + xa_dim1], ldxa, &c_b7, &a[*k + i__ + 1 
		    + (i__ + 1) * a_dim1], &c__1, (ftnlen)12);
#line 2239 "MB03XU.f"
	    i__2 = *n - i__;
#line 2239 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xa[i__ + 1 + nb1 * xa_dim1], ldxa, &c_b7, &
		    a[*k + i__ + 1 + (i__ + 1) * a_dim1], &c__1, (ftnlen)12);
#line 2241 "MB03XU.f"
	    i__2 = *k + *n;
#line 2241 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &ya[ya_offset], ldya, &
		    q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &a[(i__ + 1) * 
		    a_dim1 + 1], &c__1, (ftnlen)12);
#line 2243 "MB03XU.f"
	    i__2 = *k + *n;
#line 2243 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2243 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &b[(*k + i__ + 1) * b_dim1 + 1], &c__1, &c_b7, &a[
		    (i__ + 1) * a_dim1 + 1], &c__1, (ftnlen)12);

/*           Update YG with first Householder reflection. */

#line 2248 "MB03XU.f"
	    i__2 = *k + *n;
#line 2248 "MB03XU.f"
	    i__3 = *n - i__;
#line 2248 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], ldg, &q[i__ + (i__ + 1) * q_dim1], ldq, &
		    c_b9, &yg[i__ * yg_dim1 + 1], &c__1, (ftnlen)12);
#line 2250 "MB03XU.f"
	    i__2 = *n - i__;
#line 2250 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 1 + 
		    xg_dim1], ldxg, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
#line 2252 "MB03XU.f"
	    i__2 = *n - i__;
#line 2252 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 + i__ * 
		    yg_dim1], &c__1, (ftnlen)12);
#line 2254 "MB03XU.f"
	    i__2 = *n - i__;
#line 2254 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 1 + nb1 * 
		    xg_dim1], ldxg, &q[i__ + (i__ + 1) * q_dim1], ldq, &c_b9, 
		    &dwork[pdw], &c__1, (ftnlen)9);
#line 2256 "MB03XU.f"
	    i__2 = *n - i__;
#line 2256 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 
		    + i__ * yg_dim1], &c__1, (ftnlen)12);
#line 2258 "MB03XU.f"
	    i__2 = *k + *n;
#line 2258 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2258 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[yg_offset], ldyg, 
		    &dwork[1], &c__1, &c_b7, &yg[i__ * yg_dim1 + 1], &c__1, (
		    ftnlen)12);
#line 2260 "MB03XU.f"
	    i__2 = *k + *n;
#line 2260 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2260 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &dwork[*nb + 1], &c__1, &c_b7, &yg[i__ * yg_dim1 
		    + 1], &c__1, (ftnlen)12);
#line 2262 "MB03XU.f"
	    i__2 = *k + *n;
#line 2262 "MB03XU.f"
	    d__1 = -tauq;
#line 2262 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yg[i__ * yg_dim1 + 1], &c__1);

/*           Update G(1:k+n,k+i+1). */

#line 2266 "MB03XU.f"
	    i__2 = *n - i__;
#line 2266 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &xg[*k + i__ + 1 + xg_dim1], ldxg, &c_b7, &g[*k + 
		    i__ + 1 + (*k + i__ + 1) * g_dim1], &c__1, (ftnlen)12);
#line 2268 "MB03XU.f"
	    i__2 = *n - i__;
#line 2268 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xg[*k + i__ + 1 + nb1 * xg_dim1], ldxg, &
		    c_b7, &g[*k + i__ + 1 + (*k + i__ + 1) * g_dim1], &c__1, (
		    ftnlen)12);
#line 2270 "MB03XU.f"
	    i__2 = *k + *n;
#line 2270 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yg[yg_offset], ldyg, &
		    q[(i__ + 1) * q_dim1 + 1], &c__1, &c_b7, &g[(*k + i__ + 1)
		     * g_dim1 + 1], &c__1, (ftnlen)12);
#line 2272 "MB03XU.f"
	    i__2 = *k + *n;
#line 2272 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2272 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &b[(*k + i__ + 1) * b_dim1 + 1], &c__1, &c_b7, &g[
		    (*k + i__ + 1) * g_dim1 + 1], &c__1, (ftnlen)12);
#line 2274 "MB03XU.f"
	    i__2 = i__;
#line 2274 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 2275 "MB03XU.f"
		xg[*k + i__ + 1 + j * xg_dim1] = 0.;
#line 2276 "MB03XU.f"
/* L340: */
#line 2276 "MB03XU.f"
	    }
#line 2277 "MB03XU.f"
	    i__2 = i__;
#line 2277 "MB03XU.f"
	    for (j = 1; j <= i__2; ++j) {
#line 2278 "MB03XU.f"
		xg[*k + i__ + 1 + (*nb + j) * xg_dim1] = 0.;
#line 2279 "MB03XU.f"
/* L350: */
#line 2279 "MB03XU.f"
	    }

/*           Apply rotation to [ A(1:k+n,i+1), G(1:k+n,k+i+1) ]. */

#line 2283 "MB03XU.f"
	    i__2 = *k + *n;
#line 2283 "MB03XU.f"
	    drot_(&i__2, &a[(i__ + 1) * a_dim1 + 1], &c__1, &g[(*k + i__ + 1) 
		    * g_dim1 + 1], &c__1, &c__, &s);

/*           Update YA with second Householder reflection. */

#line 2287 "MB03XU.f"
	    i__2 = *k + *n;
#line 2287 "MB03XU.f"
	    i__3 = *n - i__;
#line 2287 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &a[(i__ + 1) * a_dim1 
		    + 1], lda, &b[i__ + (*k + i__ + 1) * b_dim1], ldb, &c_b9, 
		    &ya[(i__ + *nb) * ya_dim1 + 1], &c__1, (ftnlen)12);
#line 2289 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 2289 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 2 + xa_dim1], 
		    ldxa, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &c_b9, &
		    dwork[pdw], &c__1, (ftnlen)9);
#line 2291 "MB03XU.f"
	    i__2 = *n - i__;
#line 2291 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 + (i__ + 
		    *nb) * ya_dim1], &c__1, (ftnlen)12);
#line 2293 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 2293 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xa[i__ + 2 + nb1 * 
		    xa_dim1], ldxa, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 2295 "MB03XU.f"
	    i__2 = *n - i__;
#line 2295 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[pdw], &c__1, &c_b7, &ya[*k + i__ + 1 
		    + (i__ + *nb) * ya_dim1], &c__1, (ftnlen)12);
#line 2297 "MB03XU.f"
	    i__2 = *k + *n;
#line 2297 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &ya[ya_offset], ldya, &
		    dwork[nb2 + 1], &c__1, &c_b7, &ya[(i__ + *nb) * ya_dim1 + 
		    1], &c__1, (ftnlen)12);
#line 2299 "MB03XU.f"
	    i__2 = *k + *n;
#line 2299 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2299 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &ya[nb1 * ya_dim1 + 1]
		    , ldya, &dwork[nb3 + 1], &c__1, &c_b7, &ya[(i__ + *nb) * 
		    ya_dim1 + 1], &c__1, (ftnlen)12);
#line 2301 "MB03XU.f"
	    i__2 = *k + *n;
#line 2301 "MB03XU.f"
	    d__1 = -taur[i__];
#line 2301 "MB03XU.f"
	    dscal_(&i__2, &d__1, &ya[(i__ + *nb) * ya_dim1 + 1], &c__1);

/*           Update A(1:k+n,i+1). */

#line 2305 "MB03XU.f"
	    i__2 = *k + *n;
#line 2305 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &ya[(i__ + *nb) * ya_dim1 + 1], &c__1, &a[(
		    i__ + 1) * a_dim1 + 1], &c__1);

/*           Update YG with second Householder reflection. */

#line 2309 "MB03XU.f"
	    i__2 = *k + *n;
#line 2309 "MB03XU.f"
	    i__3 = *n - i__;
#line 2309 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &g[(*k + i__ + 1) * 
		    g_dim1 + 1], ldg, &b[i__ + (*k + i__ + 1) * b_dim1], ldb, 
		    &c_b9, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1, (ftnlen)12);
#line 2311 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 2311 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 2 + 
		    xg_dim1], ldxg, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 2313 "MB03XU.f"
	    i__2 = *n - i__;
#line 2313 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &q[i__ + 1 + q_dim1], 
		    ldq, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 + (i__ + 
		    *nb) * yg_dim1], &c__1, (ftnlen)12);
#line 2315 "MB03XU.f"
	    i__2 = *n - i__ - 1;
#line 2315 "MB03XU.f"
	    dgemv_("Transpose", &i__2, &i__, &c_b7, &xg[*k + i__ + 2 + nb1 * 
		    xg_dim1], ldxg, &b[i__ + (*k + i__ + 2) * b_dim1], ldb, &
		    c_b9, &dwork[pdw], &c__1, (ftnlen)9);
#line 2317 "MB03XU.f"
	    i__2 = *n - i__;
#line 2317 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[pdw], &c__1, &c_b7, &yg[*k + i__ + 1 
		    + (i__ + *nb) * yg_dim1], &c__1, (ftnlen)12);
#line 2319 "MB03XU.f"
	    i__2 = *k + *n;
#line 2319 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b7, &yg[yg_offset], ldyg, &
		    dwork[nb2 + 1], &c__1, &c_b7, &yg[(i__ + *nb) * yg_dim1 + 
		    1], &c__1, (ftnlen)12);
#line 2321 "MB03XU.f"
	    i__2 = *k + *n;
#line 2321 "MB03XU.f"
	    i__3 = i__ - 1;
#line 2321 "MB03XU.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b7, &yg[nb1 * yg_dim1 + 1]
		    , ldyg, &dwork[nb3 + 1], &c__1, &c_b7, &yg[(i__ + *nb) * 
		    yg_dim1 + 1], &c__1, (ftnlen)12);
#line 2323 "MB03XU.f"
	    i__2 = *k + *n;
#line 2323 "MB03XU.f"
	    d__1 = -taur[i__];
#line 2323 "MB03XU.f"
	    dscal_(&i__2, &d__1, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1);

/*           Update G(1:k+n,k+i+1). */

#line 2327 "MB03XU.f"
	    i__2 = *k + *n;
#line 2327 "MB03XU.f"
	    daxpy_(&i__2, &c_b7, &yg[(i__ + *nb) * yg_dim1 + 1], &c__1, &g[(*
		    k + i__ + 1) * g_dim1 + 1], &c__1);

#line 2329 "MB03XU.f"
	    b[i__ + (*k + i__ + 1) * b_dim1] = temp;
#line 2330 "MB03XU.f"
	    q[i__ + (i__ + 1) * q_dim1] = tauq;
#line 2331 "MB03XU.f"
	    csr[(i__ << 1) - 1] = c__;
#line 2332 "MB03XU.f"
	    csr[i__ * 2] = s;
#line 2333 "MB03XU.f"
/* L360: */
#line 2333 "MB03XU.f"
	}
#line 2334 "MB03XU.f"
    }

#line 2336 "MB03XU.f"
    return 0;
/* *** Last line of MB03XU *** */
} /* mb03xu_ */

