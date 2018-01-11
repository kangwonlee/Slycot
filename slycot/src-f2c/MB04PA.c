#line 1 "MB04PA.f"
/* MB04PA.f -- translated by f2c (version 20100827).
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

#line 1 "MB04PA.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b10 = 1.;
static doublereal c_b12 = 0.;
static doublereal c_b586 = -1.;

/* Subroutine */ int mb04pa_(logical *lham, integer *n, integer *k, integer *
	nb, doublereal *a, integer *lda, doublereal *qg, integer *ldqg, 
	doublereal *xa, integer *ldxa, doublereal *xg, integer *ldxg, 
	doublereal *xq, integer *ldxq, doublereal *ya, integer *ldya, 
	doublereal *cs, doublereal *tau, doublereal *dwork)
{
    /* System generated locals */
    integer a_dim1, a_offset, qg_dim1, qg_offset, xa_dim1, xa_offset, xg_dim1,
	     xg_offset, xq_dim1, xq_offset, ya_dim1, ya_offset, i__1, i__2, 
	    i__3;
    doublereal d__1;

    /* Local variables */
    static doublereal c__;
    static integer i__, j;
    static doublereal s;
    static integer nb1, nb2;
    static doublereal aki;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal tauq;
    extern /* Subroutine */ int mb01md_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    static doublereal alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal ttemp;
    extern /* Subroutine */ int dsymv_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), dlarfg_(integer *, doublereal *,
	     doublereal *, integer *, doublereal *), dlartg_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);


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

/*     To reduce a Hamiltonian like matrix */

/*                   [  A   G  ]           T          T */
/*              H =  [       T ] ,    G = G ,    Q = Q, */
/*                   [  Q  -A  ] */

/*     or a skew-Hamiltonian like matrix */

/*                   [  A   G  ]            T          T */
/*              W =  [       T ] ,    G = -G ,   Q = -Q, */
/*                   [  Q   A  ] */

/*     so that elements below the (k+1)-th subdiagonal in the first nb */
/*     columns of the (k+n)-by-n matrix A, and offdiagonal elements */
/*     in the first nb columns and rows of the n-by-n matrix Q are zero. */

/*     The reduction is performed by an orthogonal symplectic */
/*     transformation UU'*H*UU and matrices U, XA, XG, XQ, and YA are */
/*     returned so that */

/*                    [ Aout + U*XA'+ YA*U'   Gout + U*XG'+ XG*U' ] */
/*         UU'*H*UU = [                                           ]. */
/*                    [ Qout + U*XQ'+ XQ*U'  -Aout'- XA*U'- U*YA' ] */

/*     Similarly, */

/*                    [ Aout + U*XA'+ YA*U'   Gout + U*XG'- XG*U' ] */
/*         UU'*W*UU = [                                           ]. */
/*                    [ Qout + U*XQ'- XQ*U'   Aout'+ XA*U'+ U*YA' ] */

/*     This is an auxiliary routine called by MB04PB. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     LHAM    LOGICAL */
/*             Specifies the type of matrix to be reduced: */
/*             = .FALSE. :  skew-Hamiltonian like W; */
/*             = .TRUE.  :  Hamiltonian like H. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     K       (input) INTEGER */
/*             The offset of the reduction. Elements below the (K+1)-th */
/*             subdiagonal in the first NB columns of A are reduced */
/*             to zero.  K >= 0. */

/*     NB      (input) INTEGER */
/*             The number of columns/rows to be reduced.  N > NB >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading (K+N)-by-N part of this array must */
/*             contain the matrix A. */
/*             On exit, the leading (K+N)-by-N part of this array */
/*             contains the matrix Aout and in the zero part */
/*             information about the elementary reflectors used to */
/*             compute the reduction. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,K+N). */

/*     QG      (input/output) DOUBLE PRECISION array, dimension */
/*                            (LDQG,N+1) */
/*             On entry, the leading N+K-by-N+1 part of this array must */
/*             contain in the bottom left part the lower triangular part */
/*             of the N-by-N matrix Q and in the remainder the upper */
/*             trapezoidal part of the last N columns of the N+K-by-N+K */
/*             matrix G. */
/*             On exit, the leading N+K-by-N+1 part of this array */
/*             contains parts of the matrices Q and G in the same fashion */
/*             as on entry only that the zero parts of Q contain */
/*             information about the elementary reflectors used to */
/*             compute the reduction. Note that if LHAM = .FALSE. then */
/*             the (K-1)-th and K-th subdiagonals are not referenced. */

/*     LDQG    INTEGER */
/*             The leading dimension of the array QG. LDQG >= MAX(1,N+K). */

/*     XA      (output) DOUBLE PRECISION array, dimension (LDXA,2*NB) */
/*             On exit, the leading N-by-(2*NB) part of this array */
/*             contains the matrix XA. */

/*     LDXA    INTEGER */
/*             The leading dimension of the array XA.  LDXA >= MAX(1,N). */

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

/*     CS      (output) DOUBLE PRECISION array, dimension (2*NB) */
/*             On exit, the first 2*NB elements of this array contain the */
/*             cosines and sines of the symplectic Givens rotations used */
/*             to compute the reduction. */

/*     TAU     (output) DOUBLE PRECISION array, dimension (NB) */
/*             On exit, the first NB elements of this array contain the */
/*             scalar factors of some of the elementary reflectors. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (3*NB) */

/*     METHOD */

/*     For details regarding the representation of the orthogonal */
/*     symplectic matrix UU within the arrays A, QG, CS, TAU see the */
/*     description of MB04PU. */

/*     The contents of A and QG on exit are illustrated by the following */
/*     example with n = 5, k = 2 and nb = 2: */

/*           ( a  r  r  a  a  )         ( g  g  r  r  g  g  ) */
/*           ( a  r  r  a  a  )         ( g  g  r  r  g  g  ) */
/*           ( a  r  r  a  a  )         ( q  g  r  r  g  g  ) */
/*       A = ( r  r  r  r  r  ),   QG = ( t  r  r  r  r  r  ), */
/*           ( u2 r  r  r  r  )         ( u1 t  r  r  r  r  ) */
/*           ( u2 u2 r  a  a  )         ( u1 u1 r  q  g  g  ) */
/*           ( u2 u2 r  a  a  )         ( u1 u1 r  q  q  g  ) */

/*     where a, g and q denote elements of the original matrices, r */
/*     denotes a modified element, t denotes a scalar factor of an */
/*     applied elementary reflector and ui denote elements of the */
/*     matrix U. */

/*     REFERENCES */

/*     [1] C. F. VAN LOAN: */
/*         A symplectic method for approximating all the eigenvalues of */
/*         a Hamiltonian matrix. */
/*         Linear Algebra and its Applications, 61, pp. 233-251, 1984. */

/*     [2] D. KRESSNER: */
/*         Block algorithms for orthogonal symplectic factorizations. */
/*         BIT, 43 (4), pp. 775-790, 2003. */

/*     CONTRIBUTORS */

/*     D. Kressner (Technical Univ. Berlin, Germany) and */
/*     P. Benner (Technical Univ. Chemnitz, Germany), December 2003. */

/*     REVISIONS */

/*     V. Sima, Nov. 2008 (SLICOT version of the HAPACK routine DLAPVL). */

/*     KEYWORDS */

/*     Elementary matrix operations, Hamiltonian matrix, */
/*     skew-Hamiltonian matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Quick return if possible. */

#line 221 "MB04PA.f"
    /* Parameter adjustments */
#line 221 "MB04PA.f"
    a_dim1 = *lda;
#line 221 "MB04PA.f"
    a_offset = 1 + a_dim1;
#line 221 "MB04PA.f"
    a -= a_offset;
#line 221 "MB04PA.f"
    qg_dim1 = *ldqg;
#line 221 "MB04PA.f"
    qg_offset = 1 + qg_dim1;
#line 221 "MB04PA.f"
    qg -= qg_offset;
#line 221 "MB04PA.f"
    xa_dim1 = *ldxa;
#line 221 "MB04PA.f"
    xa_offset = 1 + xa_dim1;
#line 221 "MB04PA.f"
    xa -= xa_offset;
#line 221 "MB04PA.f"
    xg_dim1 = *ldxg;
#line 221 "MB04PA.f"
    xg_offset = 1 + xg_dim1;
#line 221 "MB04PA.f"
    xg -= xg_offset;
#line 221 "MB04PA.f"
    xq_dim1 = *ldxq;
#line 221 "MB04PA.f"
    xq_offset = 1 + xq_dim1;
#line 221 "MB04PA.f"
    xq -= xq_offset;
#line 221 "MB04PA.f"
    ya_dim1 = *ldya;
#line 221 "MB04PA.f"
    ya_offset = 1 + ya_dim1;
#line 221 "MB04PA.f"
    ya -= ya_offset;
#line 221 "MB04PA.f"
    --cs;
#line 221 "MB04PA.f"
    --tau;
#line 221 "MB04PA.f"
    --dwork;
#line 221 "MB04PA.f"

#line 221 "MB04PA.f"
    /* Function Body */
#line 221 "MB04PA.f"
    if (*n + *k <= 0) {
#line 222 "MB04PA.f"
	dwork[1] = 1.;
#line 223 "MB04PA.f"
	return 0;
#line 224 "MB04PA.f"
    }

#line 226 "MB04PA.f"
    nb1 = *nb + 1;
#line 227 "MB04PA.f"
    nb2 = *nb + nb1;

#line 229 "MB04PA.f"
    if (*lham) {
#line 230 "MB04PA.f"
	i__1 = *nb;
#line 230 "MB04PA.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Transform i-th columns of A and Q. See routine MB04PU. */

#line 234 "MB04PA.f"
	    alpha = qg[*k + i__ + 1 + i__ * qg_dim1];
#line 235 "MB04PA.f"
	    i__2 = *n - i__;
/* Computing MIN */
#line 235 "MB04PA.f"
	    i__3 = i__ + 2;
#line 235 "MB04PA.f"
	    dlarfg_(&i__2, &alpha, &qg[*k + min(i__3,*n) + i__ * qg_dim1], &
		    c__1, &tauq);
#line 236 "MB04PA.f"
	    qg[*k + i__ + 1 + i__ * qg_dim1] = 1.;
#line 237 "MB04PA.f"
	    i__2 = *n - i__;
#line 237 "MB04PA.f"
	    temp = -tauq * ddot_(&i__2, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &a[*k + i__ + 1 + i__ * a_dim1], &c__1);
#line 238 "MB04PA.f"
	    i__2 = *n - i__;
#line 238 "MB04PA.f"
	    daxpy_(&i__2, &temp, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &a[
		    *k + i__ + 1 + i__ * a_dim1], &c__1);
#line 239 "MB04PA.f"
	    aki = a[*k + i__ + 1 + i__ * a_dim1];
#line 240 "MB04PA.f"
	    dlartg_(&aki, &alpha, &c__, &s, &a[*k + i__ + 1 + i__ * a_dim1]);
#line 241 "MB04PA.f"
	    aki = a[*k + i__ + 1 + i__ * a_dim1];
#line 242 "MB04PA.f"
	    i__2 = *n - i__;
/* Computing MIN */
#line 242 "MB04PA.f"
	    i__3 = i__ + 2;
#line 242 "MB04PA.f"
	    dlarfg_(&i__2, &aki, &a[*k + min(i__3,*n) + i__ * a_dim1], &c__1, 
		    &tau[i__]);
#line 243 "MB04PA.f"
	    a[*k + i__ + 1 + i__ * a_dim1] = 1.;

/*           Update XA with first Householder reflection. */

/*           xa = H(1:n,1:n)'*u1 */
#line 248 "MB04PA.f"
	    i__2 = *n - i__;
#line 248 "MB04PA.f"
	    i__3 = *n - i__;
#line 248 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + (i__ 
		    + 1) * a_dim1], lda, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &c_b12, &xa[i__ + 1 + i__ * xa_dim1], &c__1, (
		    ftnlen)9);
/*           w1 = U1'*u1 */
#line 251 "MB04PA.f"
	    i__2 = *n - i__;
#line 251 "MB04PA.f"
	    i__3 = i__ - 1;
#line 251 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &dwork[1], &c__1, (ftnlen)9);
/*           xa = xa + XA1*w1 */
#line 254 "MB04PA.f"
	    i__2 = *n - i__;
#line 254 "MB04PA.f"
	    i__3 = i__ - 1;
#line 254 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xa[i__ + 1 + 
		    xa_dim1], ldxa, &dwork[1], &c__1, &c_b10, &xa[i__ + 1 + 
		    i__ * xa_dim1], &c__1, (ftnlen)12);
/*           w2 = U2'*u1 */
#line 257 "MB04PA.f"
	    i__2 = *n - i__;
#line 257 "MB04PA.f"
	    i__3 = i__ - 1;
#line 257 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &
		    c_b12, &dwork[nb1], &c__1, (ftnlen)9);
/*           xa = xa + XA2*w2 */
#line 260 "MB04PA.f"
	    i__2 = *n - i__;
#line 260 "MB04PA.f"
	    i__3 = i__ - 1;
#line 260 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[nb1], &c__1, &c_b10, &xa[i__ + 1 + 
		    i__ * xa_dim1], &c__1, (ftnlen)12);
/*           temp = YA1'*u1 */
#line 263 "MB04PA.f"
	    i__2 = *n - i__;
#line 263 "MB04PA.f"
	    i__3 = i__ - 1;
#line 263 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &ya[*k + i__ + 1 + 
		    ya_dim1], ldya, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &xa[i__ * xa_dim1 + 1], &c__1, (ftnlen)9);
/*           xa = xa + U1*temp */
#line 266 "MB04PA.f"
	    i__2 = *n - i__;
#line 266 "MB04PA.f"
	    i__3 = i__ - 1;
#line 266 "MB04PA.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &xa[i__ * xa_dim1 + 1], &c__1, &c_b10, &
		    xa[i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)12);
/*           temp = YA2'*u1 */
#line 269 "MB04PA.f"
	    i__2 = *n - i__;
#line 269 "MB04PA.f"
	    i__3 = i__ - 1;
#line 269 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &ya[*k + i__ + 1 + nb1 *
		     ya_dim1], ldya, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1,
		     &c_b12, &xa[i__ * xa_dim1 + 1], &c__1, (ftnlen)9);
/*           xa = xa + U2*temp */
#line 272 "MB04PA.f"
	    i__2 = *n - i__;
#line 272 "MB04PA.f"
	    i__3 = i__ - 1;
#line 272 "MB04PA.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xa[i__ * xa_dim1 + 1], &c__1, &c_b10, &xa[
		    i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)12);
/*           xa = -tauq*xa */
#line 275 "MB04PA.f"
	    i__2 = *n - i__;
#line 275 "MB04PA.f"
	    d__1 = -tauq;
#line 275 "MB04PA.f"
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + i__ * xa_dim1], &c__1);

/*           Update YA with first Householder reflection. */

/*           ya = H(1:n,1:n)*u1 */
#line 280 "MB04PA.f"
	    i__2 = *k + *n;
#line 280 "MB04PA.f"
	    i__3 = *n - i__;
#line 280 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[(i__ + 1) * 
		    a_dim1 + 1], lda, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &c_b12, &ya[i__ * ya_dim1 + 1], &c__1, (ftnlen)12);
/*           temp = XA1'*u1 */
#line 283 "MB04PA.f"
	    i__2 = *n - i__;
#line 283 "MB04PA.f"
	    i__3 = i__ - 1;
#line 283 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xa[i__ + 1 + xa_dim1], 
		    ldxa, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &c_b12, &
		    dwork[nb2], &c__1, (ftnlen)9);
/*           ya = ya + U1*temp */
#line 286 "MB04PA.f"
	    i__2 = *n - i__;
#line 286 "MB04PA.f"
	    i__3 = i__ - 1;
#line 286 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &dwork[nb2], &c__1, &c_b10, &ya[*k + i__ 
		    + 1 + i__ * ya_dim1], &c__1, (ftnlen)12);
/*           temp = XA2'*u1 */
#line 289 "MB04PA.f"
	    i__2 = *n - i__;
#line 289 "MB04PA.f"
	    i__3 = i__ - 1;
#line 289 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*           ya = ya + U2*temp */
#line 292 "MB04PA.f"
	    i__2 = *n - i__;
#line 292 "MB04PA.f"
	    i__3 = i__ - 1;
#line 292 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[nb2], &c__1, &c_b10, &ya[*k + i__ + 
		    1 + i__ * ya_dim1], &c__1, (ftnlen)12);
/*           ya = ya + YA1*w1 */
#line 295 "MB04PA.f"
	    i__2 = *k + *n;
#line 295 "MB04PA.f"
	    i__3 = i__ - 1;
#line 295 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &ya[ya_offset], ldya,
		     &dwork[1], &c__1, &c_b10, &ya[i__ * ya_dim1 + 1], &c__1, 
		    (ftnlen)12);
/*           ya = ya + YA2*w2 */
#line 298 "MB04PA.f"
	    i__2 = *k + *n;
#line 298 "MB04PA.f"
	    i__3 = i__ - 1;
#line 298 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &ya[nb1 * ya_dim1 + 
		    1], ldya, &dwork[nb1], &c__1, &c_b10, &ya[i__ * ya_dim1 + 
		    1], &c__1, (ftnlen)12);
/*           ya = -tauq*ya */
#line 301 "MB04PA.f"
	    i__2 = *k + *n;
#line 301 "MB04PA.f"
	    d__1 = -tauq;
#line 301 "MB04PA.f"
	    dscal_(&i__2, &d__1, &ya[i__ * ya_dim1 + 1], &c__1);
/*           temp = -tauq*ya'*u1 */
#line 303 "MB04PA.f"
	    i__2 = *n - i__;
#line 303 "MB04PA.f"
	    temp = -tauq * ddot_(&i__2, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &ya[*k + i__ + 1 + i__ * ya_dim1], &c__1);
/*           ya = ya + temp*u1 */
#line 305 "MB04PA.f"
	    i__2 = *n - i__;
#line 305 "MB04PA.f"
	    daxpy_(&i__2, &temp, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &
		    ya[*k + i__ + 1 + i__ * ya_dim1], &c__1);

/*           Update (i+1)-th column of A. */

/*           A(:,i+1) = A(:,i+1) + U1 * XA1(i+1,:)'; */
#line 310 "MB04PA.f"
	    i__2 = *n - i__;
#line 310 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &xa[i__ + 1 + xa_dim1], ldxa, &c_b10, &a[*
		    k + i__ + 1 + (i__ + 1) * a_dim1], &c__1, (ftnlen)12);
/*           A(:,i+1) = A(:,i+1) + U2 * XA2(i+1,:)'; */
#line 313 "MB04PA.f"
	    i__2 = *n - i__;
#line 313 "MB04PA.f"
	    i__3 = i__ - 1;
#line 313 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xa[i__ + 1 + nb1 * xa_dim1], ldxa, &c_b10, 
		    &a[*k + i__ + 1 + (i__ + 1) * a_dim1], &c__1, (ftnlen)12);
/*           A(:,i+1) = A(:,i+1) + YA1 * U1(i+1,:)'; */
#line 316 "MB04PA.f"
	    i__2 = *n + *k;
#line 316 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &ya[ya_offset], ldya, 
		    &qg[*k + i__ + 1 + qg_dim1], ldqg, &c_b10, &a[(i__ + 1) * 
		    a_dim1 + 1], &c__1, (ftnlen)12);
/*           A(:,i+1) = A(:,i+1) + YA2 * U2(i+1,:)'; */
#line 319 "MB04PA.f"
	    i__2 = *n + *k;
#line 319 "MB04PA.f"
	    i__3 = i__ - 1;
#line 319 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &ya[nb1 * ya_dim1 + 
		    1], ldya, &a[*k + i__ + 1 + a_dim1], lda, &c_b10, &a[(i__ 
		    + 1) * a_dim1 + 1], &c__1, (ftnlen)12);

/*           Update (i+1)-th row of A. */

#line 324 "MB04PA.f"
	    if (*n > i__ + 1) {
/*              A(i+1,i+2:n) = A(i+1,i+2:n) + U1(i+1,:)*XA1(i+2:n,:)' */
#line 326 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 326 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__, &c_b10, &xa[i__ + 2 + 
			xa_dim1], ldxa, &qg[*k + i__ + 1 + qg_dim1], ldqg, &
			c_b10, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda, (
			ftnlen)12);
/*              A(i+1,i+2:n) = A(i+1,i+2:n) + U2(i+1,:)*XA2(i+2:n,:)' */
#line 330 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 330 "MB04PA.f"
		i__3 = i__ - 1;
#line 330 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &xa[i__ + 2 + 
			nb1 * xa_dim1], ldxa, &a[*k + i__ + 1 + a_dim1], lda, 
			&c_b10, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda, (
			ftnlen)12);
/*              A(i+1,i+2:n) = A(i+1,i+2:n) + YA1(i+1,:) * U1(i+2:n,:)' */
#line 334 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 334 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &ya[*k + i__ + 1 + ya_dim1], ldya, &
			c_b10, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda, (
			ftnlen)12);
/*              A(i+1,i+2:n) = A(i+1,i+2:n) + YA2(i+1,:) * U2(i+2:n,:)' */
#line 338 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 338 "MB04PA.f"
		i__3 = i__ - 1;
#line 338 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 
			+ a_dim1], lda, &ya[*k + i__ + 1 + nb1 * ya_dim1], 
			ldya, &c_b10, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], 
			lda, (ftnlen)12);
#line 341 "MB04PA.f"
	    }

/*           Annihilate updated parts in YA. */

#line 345 "MB04PA.f"
	    i__2 = i__;
#line 345 "MB04PA.f"
	    for (j = 1; j <= i__2; ++j) {
#line 346 "MB04PA.f"
		ya[*k + i__ + 1 + j * ya_dim1] = 0.;
#line 347 "MB04PA.f"
/* L10: */
#line 347 "MB04PA.f"
	    }
#line 348 "MB04PA.f"
	    i__2 = i__ - 1;
#line 348 "MB04PA.f"
	    for (j = 1; j <= i__2; ++j) {
#line 349 "MB04PA.f"
		ya[*k + i__ + 1 + (*nb + j) * ya_dim1] = 0.;
#line 350 "MB04PA.f"
/* L20: */
#line 350 "MB04PA.f"
	    }

/*           Update XQ with first Householder reflection. */

/*           xq = Q*u1 */
#line 355 "MB04PA.f"
	    i__2 = *n - i__;
#line 355 "MB04PA.f"
	    dsymv_("Lower", &i__2, &c_b10, &qg[*k + i__ + 1 + (i__ + 1) * 
		    qg_dim1], ldqg, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &xq[i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)5);
/*           xq = xq + XQ1*w1 */
#line 358 "MB04PA.f"
	    i__2 = *n - i__;
#line 358 "MB04PA.f"
	    i__3 = i__ - 1;
#line 358 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xq[i__ + 1 + 
		    xq_dim1], ldxq, &dwork[1], &c__1, &c_b10, &xq[i__ + 1 + 
		    i__ * xq_dim1], &c__1, (ftnlen)12);
/*           xq = xq + XQ2*w2 */
#line 361 "MB04PA.f"
	    i__2 = *n - i__;
#line 361 "MB04PA.f"
	    i__3 = i__ - 1;
#line 361 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[nb1], &c__1, &c_b10, &xq[i__ + 1 + 
		    i__ * xq_dim1], &c__1, (ftnlen)12);
/*           temp = XQ1'*u1 */
#line 364 "MB04PA.f"
	    i__2 = *n - i__;
#line 364 "MB04PA.f"
	    i__3 = i__ - 1;
#line 364 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xq[i__ + 1 + xq_dim1], 
		    ldxq, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &c_b12, &
		    xq[i__ * xq_dim1 + 1], &c__1, (ftnlen)9);
/*           xq = xq + U1*temp */
#line 367 "MB04PA.f"
	    i__2 = *n - i__;
#line 367 "MB04PA.f"
	    i__3 = i__ - 1;
#line 367 "MB04PA.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &xq[i__ * xq_dim1 + 1], &c__1, &c_b10, &
		    xq[i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)12);
/*           temp = XQ2'*u1 */
#line 370 "MB04PA.f"
	    i__2 = *n - i__;
#line 370 "MB04PA.f"
	    i__3 = i__ - 1;
#line 370 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &xq[i__ * xq_dim1 + 1], &c__1, (ftnlen)9);
/*           xq = xq + U2*temp */
#line 373 "MB04PA.f"
	    i__2 = *n - i__;
#line 373 "MB04PA.f"
	    i__3 = i__ - 1;
#line 373 "MB04PA.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xq[i__ * xq_dim1 + 1], &c__1, &c_b10, &xq[
		    i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)12);
/*           xq = -tauq*xq */
#line 376 "MB04PA.f"
	    i__2 = *n - i__;
#line 376 "MB04PA.f"
	    d__1 = -tauq;
#line 376 "MB04PA.f"
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + i__ * xq_dim1], &c__1);
/*           temp = -tauq/2*xq'*u1 */
#line 378 "MB04PA.f"
	    i__2 = *n - i__;
#line 378 "MB04PA.f"
	    temp = tauq * -.5 * ddot_(&i__2, &qg[*k + i__ + 1 + i__ * qg_dim1]
		    , &c__1, &xq[i__ + 1 + i__ * xq_dim1], &c__1);
/*           xq = xq + temp*u1 */
#line 380 "MB04PA.f"
	    i__2 = *n - i__;
#line 380 "MB04PA.f"
	    daxpy_(&i__2, &temp, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &
		    xq[i__ + 1 + i__ * xq_dim1], &c__1);

/*           Update (i+1)-th column and row of Q. */

/*           Q(:,i+1) = Q(:,i+1) + U1 * XQ1(i+1,:)'; */
#line 385 "MB04PA.f"
	    i__2 = *n - i__;
#line 385 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &xq[i__ + 1 + xq_dim1], ldxq, &c_b10, &qg[
		    *k + i__ + 1 + (i__ + 1) * qg_dim1], &c__1, (ftnlen)12);
/*           Q(:,i+1) = Q(:,i+1) + U2 * XQ2(i+1,:)'; */
#line 388 "MB04PA.f"
	    i__2 = *n - i__;
#line 388 "MB04PA.f"
	    i__3 = i__ - 1;
#line 388 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xq[i__ + 1 + nb1 * xq_dim1], ldxq, &c_b10, 
		    &qg[*k + i__ + 1 + (i__ + 1) * qg_dim1], &c__1, (ftnlen)
		    12);
/*           Q(:,i+1) = Q(:,i+1) + XQ1 * U1(i+1,:)'; */
#line 391 "MB04PA.f"
	    i__2 = *n - i__;
#line 391 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &xq[i__ + 1 + xq_dim1]
		    , ldxq, &qg[*k + i__ + 1 + qg_dim1], ldqg, &c_b10, &qg[*k 
		    + i__ + 1 + (i__ + 1) * qg_dim1], &c__1, (ftnlen)12);
/*           Q(:,i+1) = Q(:,i+1) + XQ2 * U2(i+1,:)'; */
#line 394 "MB04PA.f"
	    i__2 = *n - i__;
#line 394 "MB04PA.f"
	    i__3 = i__ - 1;
#line 394 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &a[*k + i__ + 1 + a_dim1], lda, &c_b10, &
		    qg[*k + i__ + 1 + (i__ + 1) * qg_dim1], &c__1, (ftnlen)12)
		    ;

/*           Update XG with first Householder reflection. */

/*           xg = G*u1 */
#line 400 "MB04PA.f"
	    i__2 = *k + i__;
#line 400 "MB04PA.f"
	    i__3 = *n - i__;
#line 400 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &qg[(i__ + 2) * 
		    qg_dim1 + 1], ldqg, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &c_b12, &xg[i__ * xg_dim1 + 1], &c__1, (ftnlen)12);
#line 402 "MB04PA.f"
	    i__2 = *n - i__;
#line 402 "MB04PA.f"
	    dsymv_("Upper", &i__2, &c_b10, &qg[*k + i__ + 1 + (i__ + 2) * 
		    qg_dim1], ldqg, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &xg[*k + i__ + 1 + i__ * xg_dim1], &c__1, (ftnlen)
		    5);
/*           xg = xg + XG1*w1 */
#line 405 "MB04PA.f"
	    i__2 = *k + *n;
#line 405 "MB04PA.f"
	    i__3 = i__ - 1;
#line 405 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xg[xg_offset], ldxg,
		     &dwork[1], &c__1, &c_b10, &xg[i__ * xg_dim1 + 1], &c__1, 
		    (ftnlen)12);
/*           xg = xg + XG2*w2 */
#line 408 "MB04PA.f"
	    i__2 = *k + *n;
#line 408 "MB04PA.f"
	    i__3 = i__ - 1;
#line 408 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xg[nb1 * xg_dim1 + 
		    1], ldxg, &dwork[nb1], &c__1, &c_b10, &xg[i__ * xg_dim1 + 
		    1], &c__1, (ftnlen)12);
/*           temp = XG1'*u1 */
#line 411 "MB04PA.f"
	    i__2 = *n - i__;
#line 411 "MB04PA.f"
	    i__3 = i__ - 1;
#line 411 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xg[*k + i__ + 1 + 
		    xg_dim1], ldxq, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*           xg = xg + U1*temp */
#line 414 "MB04PA.f"
	    i__2 = *n - i__;
#line 414 "MB04PA.f"
	    i__3 = i__ - 1;
#line 414 "MB04PA.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &dwork[nb2], &c__1, &c_b10, &xg[*k + i__ 
		    + 1 + i__ * xg_dim1], &c__1, (ftnlen)12);
/*           temp = XG2'*u1 */
#line 417 "MB04PA.f"
	    i__2 = *n - i__;
#line 417 "MB04PA.f"
	    i__3 = i__ - 1;
#line 417 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xg[*k + i__ + 1 + nb1 *
		     xg_dim1], ldxq, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1,
		     &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*           xg = xg + U2*temp */
#line 420 "MB04PA.f"
	    i__2 = *n - i__;
#line 420 "MB04PA.f"
	    i__3 = i__ - 1;
#line 420 "MB04PA.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[nb2], &c__1, &c_b10, &xg[*k + i__ + 
		    1 + i__ * xg_dim1], &c__1, (ftnlen)12);
/*           xg = -tauq*xg */
#line 423 "MB04PA.f"
	    i__2 = *n + *k;
#line 423 "MB04PA.f"
	    d__1 = -tauq;
#line 423 "MB04PA.f"
	    dscal_(&i__2, &d__1, &xg[i__ * xg_dim1 + 1], &c__1);
/*           temp = -tauq/2*xq'*u1 */
#line 425 "MB04PA.f"
	    i__2 = *n - i__;
#line 425 "MB04PA.f"
	    temp = tauq * -.5 * ddot_(&i__2, &qg[*k + i__ + 1 + i__ * qg_dim1]
		    , &c__1, &xg[*k + i__ + 1 + i__ * xg_dim1], &c__1);
/*           xg = xg + temp*u1 */
#line 428 "MB04PA.f"
	    i__2 = *n - i__;
#line 428 "MB04PA.f"
	    daxpy_(&i__2, &temp, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &
		    xg[*k + i__ + 1 + i__ * xg_dim1], &c__1);

/*           Update (i+1)-th column and row of G. */

/*           G(:,i+1) = G(:,i+1) + XG1 * U1(i+1,:)'; */
#line 433 "MB04PA.f"
	    i__2 = *k + i__;
#line 433 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &xg[xg_offset], ldxg, 
		    &qg[*k + i__ + 1 + qg_dim1], ldqg, &c_b10, &qg[(i__ + 2) *
		     qg_dim1 + 1], &c__1, (ftnlen)12);
/*           G(:,i+1) = G(:,i+1) + XG2 * U2(i+1,:)'; */
#line 436 "MB04PA.f"
	    i__2 = *k + i__;
#line 436 "MB04PA.f"
	    i__3 = i__ - 1;
#line 436 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xg[nb1 * xg_dim1 + 
		    1], ldxg, &a[*k + i__ + 1 + a_dim1], lda, &c_b10, &qg[(
		    i__ + 2) * qg_dim1 + 1], &c__1, (ftnlen)12);
/*           G(:,i+1) = G(:,i+1) + XG1 * U1(i+1,:)'; */
#line 439 "MB04PA.f"
	    i__2 = *n - i__;
#line 439 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &xg[*k + i__ + 1 + 
		    xg_dim1], ldxg, &qg[*k + i__ + 1 + qg_dim1], ldqg, &c_b10,
		     &qg[*k + i__ + 1 + (i__ + 2) * qg_dim1], ldqg, (ftnlen)
		    12);
/*           G(:,i+1) = G(:,i+1) + XG2 * U2(i+1,:)'; */
#line 442 "MB04PA.f"
	    i__2 = *n - i__;
#line 442 "MB04PA.f"
	    i__3 = i__ - 1;
#line 442 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xg[*k + i__ + 1 + 
		    nb1 * xg_dim1], ldxg, &a[*k + i__ + 1 + a_dim1], lda, &
		    c_b10, &qg[*k + i__ + 1 + (i__ + 2) * qg_dim1], ldqg, (
		    ftnlen)12);
/*           G(:,i+1) = G(:,i+1) + U1 * XG1(i+1,:)'; */
#line 446 "MB04PA.f"
	    i__2 = *n - i__;
#line 446 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &xg[*k + i__ + 1 + xg_dim1], ldxg, &c_b10,
		     &qg[*k + i__ + 1 + (i__ + 2) * qg_dim1], ldqg, (ftnlen)
		    12);
/*           G(:,i+1) = G(:,i+1) + U2 * XG2(i+1,:)'; */
#line 449 "MB04PA.f"
	    i__2 = *n - i__;
#line 449 "MB04PA.f"
	    i__3 = i__ - 1;
#line 449 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xg[*k + i__ + 1 + nb1 * xg_dim1], ldxg, &
		    c_b10, &qg[*k + i__ + 1 + (i__ + 2) * qg_dim1], ldqg, (
		    ftnlen)12);

/*           Annihilate updated parts in XG. */

#line 454 "MB04PA.f"
	    i__2 = i__;
#line 454 "MB04PA.f"
	    for (j = 1; j <= i__2; ++j) {
#line 455 "MB04PA.f"
		xg[*k + i__ + 1 + j * xg_dim1] = 0.;
#line 456 "MB04PA.f"
/* L30: */
#line 456 "MB04PA.f"
	    }
#line 457 "MB04PA.f"
	    i__2 = i__ - 1;
#line 457 "MB04PA.f"
	    for (j = 1; j <= i__2; ++j) {
#line 458 "MB04PA.f"
		xg[*k + i__ + 1 + (*nb + j) * xg_dim1] = 0.;
#line 459 "MB04PA.f"
/* L40: */
#line 459 "MB04PA.f"
	    }

/*           Apply orthogonal symplectic Givens rotation. */

#line 463 "MB04PA.f"
	    i__2 = *k + i__;
#line 463 "MB04PA.f"
	    drot_(&i__2, &a[(i__ + 1) * a_dim1 + 1], &c__1, &qg[(i__ + 2) * 
		    qg_dim1 + 1], &c__1, &c__, &s);
#line 464 "MB04PA.f"
	    if (*n > i__ + 1) {
#line 465 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 465 "MB04PA.f"
		drot_(&i__2, &a[*k + i__ + 2 + (i__ + 1) * a_dim1], &c__1, &
			qg[*k + i__ + 1 + (i__ + 3) * qg_dim1], ldqg, &c__, &
			s);
#line 467 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 467 "MB04PA.f"
		drot_(&i__2, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda, &qg[*
			k + i__ + 2 + (i__ + 1) * qg_dim1], &c__1, &c__, &s);
#line 469 "MB04PA.f"
	    }
#line 470 "MB04PA.f"
	    temp = a[*k + i__ + 1 + (i__ + 1) * a_dim1];
#line 471 "MB04PA.f"
	    ttemp = qg[*k + i__ + 1 + (i__ + 2) * qg_dim1];
#line 472 "MB04PA.f"
	    a[*k + i__ + 1 + (i__ + 1) * a_dim1] = c__ * temp + s * qg[*k + 
		    i__ + 1 + (i__ + 1) * qg_dim1];
#line 473 "MB04PA.f"
	    qg[*k + i__ + 1 + (i__ + 2) * qg_dim1] = c__ * ttemp - s * temp;
#line 474 "MB04PA.f"
	    qg[*k + i__ + 1 + (i__ + 1) * qg_dim1] = -s * temp + c__ * qg[*k 
		    + i__ + 1 + (i__ + 1) * qg_dim1];
#line 475 "MB04PA.f"
	    ttemp = -s * ttemp - c__ * temp;
#line 476 "MB04PA.f"
	    temp = a[*k + i__ + 1 + (i__ + 1) * a_dim1];
#line 477 "MB04PA.f"
	    qg[*k + i__ + 1 + (i__ + 1) * qg_dim1] = c__ * qg[*k + i__ + 1 + (
		    i__ + 1) * qg_dim1] + s * ttemp;
#line 478 "MB04PA.f"
	    a[*k + i__ + 1 + (i__ + 1) * a_dim1] = c__ * temp + s * qg[*k + 
		    i__ + 1 + (i__ + 2) * qg_dim1];
#line 479 "MB04PA.f"
	    qg[*k + i__ + 1 + (i__ + 2) * qg_dim1] = -s * temp + c__ * qg[*k 
		    + i__ + 1 + (i__ + 2) * qg_dim1];
#line 480 "MB04PA.f"
	    cs[(i__ << 1) - 1] = c__;
#line 481 "MB04PA.f"
	    cs[i__ * 2] = s;
#line 482 "MB04PA.f"
	    qg[*k + i__ + 1 + i__ * qg_dim1] = tauq;

/*           Update XA with second Householder reflection. */

/*           xa = H(1:n,1:n)'*u2 */
#line 487 "MB04PA.f"
	    i__2 = *n - i__;
#line 487 "MB04PA.f"
	    i__3 = *n - i__;
#line 487 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + (i__ 
		    + 1) * a_dim1], lda, &a[*k + i__ + 1 + i__ * a_dim1], &
		    c__1, &c_b12, &xa[i__ + 1 + (*nb + i__) * xa_dim1], &c__1,
		     (ftnlen)9);
#line 489 "MB04PA.f"
	    if (*n > i__ + 1) {
/*              w1 = U1'*u2 */
#line 491 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 491 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 + 
			qg_dim1], ldqg, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &dwork[1], &c__1, (ftnlen)9);
/*              xa = xa + XA1*w1 */
#line 494 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 494 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__, &c_b10, &xa[i__ + 2 + 
			xa_dim1], ldxa, &dwork[1], &c__1, &c_b10, &xa[i__ + 2 
			+ (*nb + i__) * xa_dim1], &c__1, (ftnlen)12);
/*              w2 = U2'*u2 */
#line 497 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 497 "MB04PA.f"
		i__3 = i__ - 1;
#line 497 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 + 
			a_dim1], lda, &a[*k + i__ + 2 + i__ * a_dim1], &c__1, 
			&c_b12, &dwork[nb1], &c__1, (ftnlen)9);
/*              xa = xa + XA2*w2 */
#line 500 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 500 "MB04PA.f"
		i__3 = i__ - 1;
#line 500 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &xa[i__ + 2 + 
			nb1 * xa_dim1], ldxa, &dwork[nb1], &c__1, &c_b10, &xa[
			i__ + 2 + (*nb + i__) * xa_dim1], &c__1, (ftnlen)12);
/*              temp = YA1'*u2 */
#line 503 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 503 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__, &c_b10, &ya[*k + i__ + 2 + 
			ya_dim1], ldya, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &xa[(*nb + i__) * xa_dim1 + 1], &c__1, (
			ftnlen)9);
/*              xa = xa + U1*temp */
#line 506 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 506 "MB04PA.f"
		dgemv_("No Transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &xa[(*nb + i__) * xa_dim1 + 1], &
			c__1, &c_b10, &xa[i__ + 2 + (*nb + i__) * xa_dim1], &
			c__1, (ftnlen)12);
/*              temp = YA2'*u1 */
#line 509 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 509 "MB04PA.f"
		i__3 = i__ - 1;
#line 509 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &ya[*k + i__ + 2 + 
			nb1 * ya_dim1], ldya, &a[*k + i__ + 2 + i__ * a_dim1],
			 &c__1, &c_b12, &xa[(*nb + i__) * xa_dim1 + 1], &c__1,
			 (ftnlen)9);
/*              xa = xa + U2*temp */
#line 512 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 512 "MB04PA.f"
		i__3 = i__ - 1;
#line 512 "MB04PA.f"
		dgemv_("No Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 
			+ a_dim1], lda, &xa[(*nb + i__) * xa_dim1 + 1], &c__1,
			 &c_b10, &xa[i__ + 2 + (*nb + i__) * xa_dim1], &c__1, 
			(ftnlen)12);
#line 514 "MB04PA.f"
	    }
/*           xa = -tau*xa */
#line 516 "MB04PA.f"
	    i__2 = *n - i__;
#line 516 "MB04PA.f"
	    d__1 = -tau[i__];
#line 516 "MB04PA.f"
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + (*nb + i__) * xa_dim1], &c__1);

/*           Update YA with second Householder reflection. */

/*           ya = H(1:n,1:n)*u2 */
#line 521 "MB04PA.f"
	    i__2 = *k + *n;
#line 521 "MB04PA.f"
	    i__3 = *n - i__;
#line 521 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[(i__ + 1) * 
		    a_dim1 + 1], lda, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, 
		    &c_b12, &ya[(*nb + i__) * ya_dim1 + 1], &c__1, (ftnlen)12)
		    ;
#line 523 "MB04PA.f"
	    if (*n > i__ + 1) {
/*              temp = XA1'*u2 */
#line 525 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 525 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__, &c_b10, &xa[i__ + 2 + 
			xa_dim1], ldxa, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*              ya = ya + U1*temp */
#line 528 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 528 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &dwork[nb2], &c__1, &c_b10, &ya[*k 
			+ i__ + 2 + (*nb + i__) * ya_dim1], &c__1, (ftnlen)12)
			;
/*              temp = XA2'*u1 */
#line 531 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 531 "MB04PA.f"
		i__3 = i__ - 1;
#line 531 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &xa[i__ + 2 + nb1 * 
			xa_dim1], ldxa, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*              ya = ya + U2*temp */
#line 534 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 534 "MB04PA.f"
		i__3 = i__ - 1;
#line 534 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 
			+ a_dim1], lda, &dwork[nb2], &c__1, &c_b10, &ya[*k + 
			i__ + 2 + (*nb + i__) * ya_dim1], &c__1, (ftnlen)12);
#line 536 "MB04PA.f"
	    }
/*           ya = ya + YA1*w1 */
#line 538 "MB04PA.f"
	    i__2 = *k + *n;
#line 538 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &ya[ya_offset], ldya, 
		    &dwork[1], &c__1, &c_b10, &ya[(*nb + i__) * ya_dim1 + 1], 
		    &c__1, (ftnlen)12);
/*           ya = ya + YA2*w2 */
#line 541 "MB04PA.f"
	    i__2 = *k + *n;
#line 541 "MB04PA.f"
	    i__3 = i__ - 1;
#line 541 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &ya[nb1 * ya_dim1 + 
		    1], ldya, &dwork[nb1], &c__1, &c_b10, &ya[(*nb + i__) * 
		    ya_dim1 + 1], &c__1, (ftnlen)12);
/*           ya = -tau*ya */
#line 544 "MB04PA.f"
	    i__2 = *k + *n;
#line 544 "MB04PA.f"
	    d__1 = -tau[i__];
#line 544 "MB04PA.f"
	    dscal_(&i__2, &d__1, &ya[(*nb + i__) * ya_dim1 + 1], &c__1);
/*           temp = -tau*ya'*u2 */
#line 546 "MB04PA.f"
	    i__2 = *n - i__;
#line 546 "MB04PA.f"
	    temp = -tau[i__] * ddot_(&i__2, &a[*k + i__ + 1 + i__ * a_dim1], &
		    c__1, &ya[*k + i__ + 1 + (*nb + i__) * ya_dim1], &c__1);
/*           ya = ya + temp*u2 */
#line 548 "MB04PA.f"
	    i__2 = *n - i__;
#line 548 "MB04PA.f"
	    daxpy_(&i__2, &temp, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &ya[*
		    k + i__ + 1 + (*nb + i__) * ya_dim1], &c__1);

/*           Update (i+1)-th column of A. */

/*           H(1:n,i+1) = H(1:n,i+1) + ya */
#line 553 "MB04PA.f"
	    i__2 = *k + *n;
#line 553 "MB04PA.f"
	    daxpy_(&i__2, &c_b10, &ya[(*nb + i__) * ya_dim1 + 1], &c__1, &a[(
		    i__ + 1) * a_dim1 + 1], &c__1);
/*           H(1:n,i+1) = H(1:n,i+1) + xa(i+1)*u2 */
#line 555 "MB04PA.f"
	    i__2 = *n - i__;
#line 555 "MB04PA.f"
	    daxpy_(&i__2, &xa[i__ + 1 + (*nb + i__) * xa_dim1], &a[*k + i__ + 
		    1 + i__ * a_dim1], &c__1, &a[*k + i__ + 1 + (i__ + 1) * 
		    a_dim1], &c__1);

/*           Update (i+1)-th row of A. */

#line 560 "MB04PA.f"
	    if (*n > i__ + 1) {
/*              H(i+1,i+2:n) = H(i+1,i+2:n) + xa(i+2:n)'; */
#line 562 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 562 "MB04PA.f"
		daxpy_(&i__2, &c_b10, &xa[i__ + 2 + (*nb + i__) * xa_dim1], &
			c__1, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda);
/*              H(i+1,i+2:n) = H(i+1,i+2:n) + YA(i+1,:) * U(i+2:n,:)' */
#line 565 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 565 "MB04PA.f"
		daxpy_(&i__2, &ya[*k + i__ + 1 + (*nb + i__) * ya_dim1], &a[*
			k + i__ + 2 + i__ * a_dim1], &c__1, &a[*k + i__ + 1 + 
			(i__ + 2) * a_dim1], lda);
#line 567 "MB04PA.f"
	    }

/*           Annihilate updated parts in YA. */

#line 571 "MB04PA.f"
	    ya[*k + i__ + 1 + (*nb + i__) * ya_dim1] = 0.;

/*           Update XQ with second Householder reflection. */

/*           xq = Q*u2 */
#line 576 "MB04PA.f"
	    i__2 = *n - i__;
#line 576 "MB04PA.f"
	    dsymv_("Lower", &i__2, &c_b10, &qg[*k + i__ + 1 + (i__ + 1) * 
		    qg_dim1], ldqg, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b12, &xq[i__ + 1 + (*nb + i__) * xq_dim1], &c__1, (
		    ftnlen)5);
#line 578 "MB04PA.f"
	    if (*n > i__ + 1) {
/*              xq = xq + XQ1*w1 */
#line 580 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 580 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__, &c_b10, &xq[i__ + 2 + 
			xq_dim1], ldxq, &dwork[1], &c__1, &c_b10, &xq[i__ + 2 
			+ (*nb + i__) * xq_dim1], &c__1, (ftnlen)12);
/*              xq = xq + XQ2*w2 */
#line 583 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 583 "MB04PA.f"
		i__3 = i__ - 1;
#line 583 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &xq[i__ + 2 + 
			nb1 * xq_dim1], ldxq, &dwork[nb1], &c__1, &c_b10, &xq[
			i__ + 2 + (*nb + i__) * xq_dim1], &c__1, (ftnlen)12);
/*              temp = XQ1'*u2 */
#line 586 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 586 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__, &c_b10, &xq[i__ + 2 + 
			xq_dim1], ldxq, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &xq[(*nb + i__) * xq_dim1 + 1], &c__1, (
			ftnlen)9);
/*              xq = xq + U1*temp */
#line 589 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 589 "MB04PA.f"
		dgemv_("No Transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &xq[(*nb + i__) * xq_dim1 + 1], &
			c__1, &c_b10, &xq[i__ + 2 + (*nb + i__) * xq_dim1], &
			c__1, (ftnlen)12);
/*              temp = XQ2'*u2 */
#line 592 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 592 "MB04PA.f"
		i__3 = i__ - 1;
#line 592 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &xq[i__ + 2 + nb1 * 
			xq_dim1], ldxq, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &xq[(*nb + i__) * xq_dim1 + 1], &c__1, (
			ftnlen)9);
/*              xq = xq + U2*temp */
#line 595 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 595 "MB04PA.f"
		i__3 = i__ - 1;
#line 595 "MB04PA.f"
		dgemv_("No Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 
			+ a_dim1], lda, &xq[(*nb + i__) * xq_dim1 + 1], &c__1,
			 &c_b10, &xq[i__ + 2 + (*nb + i__) * xq_dim1], &c__1, 
			(ftnlen)12);
#line 597 "MB04PA.f"
	    }
/*           xq = -tauq*xq */
#line 599 "MB04PA.f"
	    i__2 = *n - i__;
#line 599 "MB04PA.f"
	    d__1 = -tau[i__];
#line 599 "MB04PA.f"
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + (*nb + i__) * xq_dim1], &c__1);
/*           temp = -tauq/2*xq'*u2 */
#line 601 "MB04PA.f"
	    i__2 = *n - i__;
#line 601 "MB04PA.f"
	    temp = tau[i__] * -.5 * ddot_(&i__2, &a[*k + i__ + 1 + i__ * 
		    a_dim1], &c__1, &xq[i__ + 1 + (*nb + i__) * xq_dim1], &
		    c__1);
/*           xq = xq + temp*u2 */
#line 604 "MB04PA.f"
	    i__2 = *n - i__;
#line 604 "MB04PA.f"
	    daxpy_(&i__2, &temp, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &xq[
		    i__ + 1 + (*nb + i__) * xq_dim1], &c__1);

/*           Update (i+1)-th column and row of Q. */

#line 608 "MB04PA.f"
	    i__2 = *n - i__;
#line 608 "MB04PA.f"
	    daxpy_(&i__2, &c_b10, &xq[i__ + 1 + (*nb + i__) * xq_dim1], &c__1,
		     &qg[*k + i__ + 1 + (i__ + 1) * qg_dim1], &c__1);
/*           H(1:n,n+i+1) = H(1:n,n+i+1) + U * XQ(i+1,:)'; */
#line 610 "MB04PA.f"
	    i__2 = *n - i__;
#line 610 "MB04PA.f"
	    daxpy_(&i__2, &xq[i__ + 1 + (*nb + i__) * xq_dim1], &a[*k + i__ + 
		    1 + i__ * a_dim1], &c__1, &qg[*k + i__ + 1 + (i__ + 1) * 
		    qg_dim1], &c__1);

/*           Update XG with second Householder reflection. */

/*           xg = G*u2 */
#line 616 "MB04PA.f"
	    i__2 = *k + i__;
#line 616 "MB04PA.f"
	    i__3 = *n - i__;
#line 616 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &qg[(i__ + 2) * 
		    qg_dim1 + 1], ldqg, &a[*k + i__ + 1 + i__ * a_dim1], &
		    c__1, &c_b12, &xg[(*nb + i__) * xg_dim1 + 1], &c__1, (
		    ftnlen)12);
#line 618 "MB04PA.f"
	    i__2 = *n - i__;
#line 618 "MB04PA.f"
	    dsymv_("Upper", &i__2, &c_b10, &qg[*k + i__ + 1 + (i__ + 2) * 
		    qg_dim1], ldqg, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b12, &xg[*k + i__ + 1 + (*nb + i__) * xg_dim1], &c__1, (
		    ftnlen)5);
/*           xg = xg + XG1*w1 */
#line 621 "MB04PA.f"
	    i__2 = *k + *n;
#line 621 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &xg[xg_offset], ldxg, 
		    &dwork[1], &c__1, &c_b10, &xg[(*nb + i__) * xg_dim1 + 1], 
		    &c__1, (ftnlen)12);
/*           xg = xg + XG2*w2 */
#line 624 "MB04PA.f"
	    i__2 = *k + *n;
#line 624 "MB04PA.f"
	    i__3 = i__ - 1;
#line 624 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xg[nb1 * xg_dim1 + 
		    1], ldxg, &dwork[nb1], &c__1, &c_b10, &xg[(*nb + i__) * 
		    xg_dim1 + 1], &c__1, (ftnlen)12);
#line 626 "MB04PA.f"
	    if (*n > i__ + 1) {
/*              temp = XG1'*u2 */
#line 628 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 628 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__, &c_b10, &xg[*k + i__ + 2 + 
			xg_dim1], ldxq, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*              xg = xg + U1*temp */
#line 631 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 631 "MB04PA.f"
		dgemv_("No Transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &dwork[nb2], &c__1, &c_b10, &xg[*k 
			+ i__ + 2 + (*nb + i__) * xg_dim1], &c__1, (ftnlen)12)
			;
/*              temp = XG2'*u2 */
#line 634 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 634 "MB04PA.f"
		i__3 = i__ - 1;
#line 634 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &xg[*k + i__ + 2 + 
			nb1 * xg_dim1], ldxq, &a[*k + i__ + 2 + i__ * a_dim1],
			 &c__1, &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*              xg = xg + U2*temp */
#line 637 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 637 "MB04PA.f"
		i__3 = i__ - 1;
#line 637 "MB04PA.f"
		dgemv_("No Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 
			+ a_dim1], lda, &dwork[nb2], &c__1, &c_b10, &xg[*k + 
			i__ + 2 + (*nb + i__) * xg_dim1], &c__1, (ftnlen)12);
#line 639 "MB04PA.f"
	    }
/*           xg = -tauq*xg */
#line 641 "MB04PA.f"
	    i__2 = *n + *k;
#line 641 "MB04PA.f"
	    d__1 = -tau[i__];
#line 641 "MB04PA.f"
	    dscal_(&i__2, &d__1, &xg[(*nb + i__) * xg_dim1 + 1], &c__1);
/*           temp = -tauq/2*xg'*u1 */
#line 643 "MB04PA.f"
	    i__2 = *n - i__;
#line 643 "MB04PA.f"
	    temp = tau[i__] * -.5 * ddot_(&i__2, &a[*k + i__ + 1 + i__ * 
		    a_dim1], &c__1, &xg[*k + i__ + 1 + (*nb + i__) * xg_dim1],
		     &c__1);
/*           xg = xg + temp*u1 */
#line 646 "MB04PA.f"
	    i__2 = *n - i__;
#line 646 "MB04PA.f"
	    daxpy_(&i__2, &temp, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &xg[*
		    k + i__ + 1 + (*nb + i__) * xg_dim1], &c__1);

/*           Update (i+1)-th column and row of G. */

#line 650 "MB04PA.f"
	    i__2 = *k + i__;
#line 650 "MB04PA.f"
	    daxpy_(&i__2, &c_b10, &xg[(*nb + i__) * xg_dim1 + 1], &c__1, &qg[(
		    i__ + 2) * qg_dim1 + 1], &c__1);
#line 651 "MB04PA.f"
	    i__2 = *n - i__;
#line 651 "MB04PA.f"
	    daxpy_(&i__2, &c_b10, &xg[*k + i__ + 1 + (*nb + i__) * xg_dim1], &
		    c__1, &qg[*k + i__ + 1 + (i__ + 2) * qg_dim1], ldqg);
#line 653 "MB04PA.f"
	    i__2 = *n - i__;
#line 653 "MB04PA.f"
	    daxpy_(&i__2, &xg[*k + i__ + 1 + (*nb + i__) * xg_dim1], &a[*k + 
		    i__ + 1 + i__ * a_dim1], &c__1, &qg[*k + i__ + 1 + (i__ + 
		    2) * qg_dim1], ldqg);

/*           Annihilate updated parts in XG. */

#line 658 "MB04PA.f"
	    xg[*k + i__ + 1 + (*nb + i__) * xg_dim1] = 0.;

#line 660 "MB04PA.f"
	    a[*k + i__ + 1 + i__ * a_dim1] = aki;
#line 661 "MB04PA.f"
/* L50: */
#line 661 "MB04PA.f"
	}
#line 662 "MB04PA.f"
    } else {
#line 663 "MB04PA.f"
	i__1 = *nb;
#line 663 "MB04PA.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Transform i-th columns of A and Q. */

#line 667 "MB04PA.f"
	    alpha = qg[*k + i__ + 1 + i__ * qg_dim1];
#line 668 "MB04PA.f"
	    i__2 = *n - i__;
/* Computing MIN */
#line 668 "MB04PA.f"
	    i__3 = i__ + 2;
#line 668 "MB04PA.f"
	    dlarfg_(&i__2, &alpha, &qg[*k + min(i__3,*n) + i__ * qg_dim1], &
		    c__1, &tauq);
#line 669 "MB04PA.f"
	    qg[*k + i__ + 1 + i__ * qg_dim1] = 1.;
#line 670 "MB04PA.f"
	    i__2 = *n - i__;
#line 670 "MB04PA.f"
	    temp = -tauq * ddot_(&i__2, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &a[*k + i__ + 1 + i__ * a_dim1], &c__1);
#line 671 "MB04PA.f"
	    i__2 = *n - i__;
#line 671 "MB04PA.f"
	    daxpy_(&i__2, &temp, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &a[
		    *k + i__ + 1 + i__ * a_dim1], &c__1);
#line 672 "MB04PA.f"
	    aki = a[*k + i__ + 1 + i__ * a_dim1];
#line 673 "MB04PA.f"
	    dlartg_(&aki, &alpha, &c__, &s, &a[*k + i__ + 1 + i__ * a_dim1]);
#line 674 "MB04PA.f"
	    aki = a[*k + i__ + 1 + i__ * a_dim1];
#line 675 "MB04PA.f"
	    i__2 = *n - i__;
/* Computing MIN */
#line 675 "MB04PA.f"
	    i__3 = i__ + 2;
#line 675 "MB04PA.f"
	    dlarfg_(&i__2, &aki, &a[*k + min(i__3,*n) + i__ * a_dim1], &c__1, 
		    &tau[i__]);
#line 676 "MB04PA.f"
	    a[*k + i__ + 1 + i__ * a_dim1] = 1.;

/*           Update XA with first Householder reflection. */

/*           xa = H(1:n,1:n)'*u1 */
#line 681 "MB04PA.f"
	    i__2 = *n - i__;
#line 681 "MB04PA.f"
	    i__3 = *n - i__;
#line 681 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + (i__ 
		    + 1) * a_dim1], lda, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &c_b12, &xa[i__ + 1 + i__ * xa_dim1], &c__1, (
		    ftnlen)9);
/*           w1 = U1'*u1 */
#line 684 "MB04PA.f"
	    i__2 = *n - i__;
#line 684 "MB04PA.f"
	    i__3 = i__ - 1;
#line 684 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &dwork[1], &c__1, (ftnlen)9);
/*           xa = xa + XA1*w1 */
#line 687 "MB04PA.f"
	    i__2 = *n - i__;
#line 687 "MB04PA.f"
	    i__3 = i__ - 1;
#line 687 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xa[i__ + 1 + 
		    xa_dim1], ldxa, &dwork[1], &c__1, &c_b10, &xa[i__ + 1 + 
		    i__ * xa_dim1], &c__1, (ftnlen)12);
/*           w2 = U2'*u1 */
#line 690 "MB04PA.f"
	    i__2 = *n - i__;
#line 690 "MB04PA.f"
	    i__3 = i__ - 1;
#line 690 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &
		    c_b12, &dwork[nb1], &c__1, (ftnlen)9);
/*           xa = xa + XA2*w2 */
#line 693 "MB04PA.f"
	    i__2 = *n - i__;
#line 693 "MB04PA.f"
	    i__3 = i__ - 1;
#line 693 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &dwork[nb1], &c__1, &c_b10, &xa[i__ + 1 + 
		    i__ * xa_dim1], &c__1, (ftnlen)12);
/*           temp = YA1'*u1 */
#line 696 "MB04PA.f"
	    i__2 = *n - i__;
#line 696 "MB04PA.f"
	    i__3 = i__ - 1;
#line 696 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &ya[*k + i__ + 1 + 
		    ya_dim1], ldya, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &xa[i__ * xa_dim1 + 1], &c__1, (ftnlen)9);
/*           xa = xa + U1*temp */
#line 699 "MB04PA.f"
	    i__2 = *n - i__;
#line 699 "MB04PA.f"
	    i__3 = i__ - 1;
#line 699 "MB04PA.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &xa[i__ * xa_dim1 + 1], &c__1, &c_b10, &
		    xa[i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)12);
/*           temp = YA2'*u1 */
#line 702 "MB04PA.f"
	    i__2 = *n - i__;
#line 702 "MB04PA.f"
	    i__3 = i__ - 1;
#line 702 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &ya[*k + i__ + 1 + nb1 *
		     ya_dim1], ldya, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1,
		     &c_b12, &xa[i__ * xa_dim1 + 1], &c__1, (ftnlen)9);
/*           xa = xa + U2*temp */
#line 705 "MB04PA.f"
	    i__2 = *n - i__;
#line 705 "MB04PA.f"
	    i__3 = i__ - 1;
#line 705 "MB04PA.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xa[i__ * xa_dim1 + 1], &c__1, &c_b10, &xa[
		    i__ + 1 + i__ * xa_dim1], &c__1, (ftnlen)12);
/*           xa = -tauq*xa */
#line 708 "MB04PA.f"
	    i__2 = *n - i__;
#line 708 "MB04PA.f"
	    d__1 = -tauq;
#line 708 "MB04PA.f"
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + i__ * xa_dim1], &c__1);

/*           Update YA with first Householder reflection. */

/*           ya = H(1:n,1:n)*u1 */
#line 713 "MB04PA.f"
	    i__2 = *k + *n;
#line 713 "MB04PA.f"
	    i__3 = *n - i__;
#line 713 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[(i__ + 1) * 
		    a_dim1 + 1], lda, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &c_b12, &ya[i__ * ya_dim1 + 1], &c__1, (ftnlen)12);
/*           temp = XA1'*u1 */
#line 716 "MB04PA.f"
	    i__2 = *n - i__;
#line 716 "MB04PA.f"
	    i__3 = i__ - 1;
#line 716 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xa[i__ + 1 + xa_dim1], 
		    ldxa, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &c_b12, &
		    dwork[nb2], &c__1, (ftnlen)9);
/*           ya = ya + U1*temp */
#line 719 "MB04PA.f"
	    i__2 = *n - i__;
#line 719 "MB04PA.f"
	    i__3 = i__ - 1;
#line 719 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &dwork[nb2], &c__1, &c_b10, &ya[*k + i__ 
		    + 1 + i__ * ya_dim1], &c__1, (ftnlen)12);
/*           temp = XA2'*u1 */
#line 722 "MB04PA.f"
	    i__2 = *n - i__;
#line 722 "MB04PA.f"
	    i__3 = i__ - 1;
#line 722 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xa[i__ + 1 + nb1 * 
		    xa_dim1], ldxa, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*           ya = ya + U2*temp */
#line 725 "MB04PA.f"
	    i__2 = *n - i__;
#line 725 "MB04PA.f"
	    i__3 = i__ - 1;
#line 725 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[nb2], &c__1, &c_b10, &ya[*k + i__ + 
		    1 + i__ * ya_dim1], &c__1, (ftnlen)12);
/*           ya = ya + YA1*w1 */
#line 728 "MB04PA.f"
	    i__2 = *k + *n;
#line 728 "MB04PA.f"
	    i__3 = i__ - 1;
#line 728 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &ya[ya_offset], ldya,
		     &dwork[1], &c__1, &c_b10, &ya[i__ * ya_dim1 + 1], &c__1, 
		    (ftnlen)12);
/*           ya = ya + YA2*w2 */
#line 731 "MB04PA.f"
	    i__2 = *k + *n;
#line 731 "MB04PA.f"
	    i__3 = i__ - 1;
#line 731 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &ya[nb1 * ya_dim1 + 
		    1], ldya, &dwork[nb1], &c__1, &c_b10, &ya[i__ * ya_dim1 + 
		    1], &c__1, (ftnlen)12);
/*           ya = -tauq*ya */
#line 734 "MB04PA.f"
	    i__2 = *k + *n;
#line 734 "MB04PA.f"
	    d__1 = -tauq;
#line 734 "MB04PA.f"
	    dscal_(&i__2, &d__1, &ya[i__ * ya_dim1 + 1], &c__1);
/*           temp = -tauq*ya'*u1 */
#line 736 "MB04PA.f"
	    i__2 = *n - i__;
#line 736 "MB04PA.f"
	    temp = -tauq * ddot_(&i__2, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &ya[*k + i__ + 1 + i__ * ya_dim1], &c__1);
/*           ya = ya + temp*u1 */
#line 738 "MB04PA.f"
	    i__2 = *n - i__;
#line 738 "MB04PA.f"
	    daxpy_(&i__2, &temp, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &
		    ya[*k + i__ + 1 + i__ * ya_dim1], &c__1);

/*           Update (i+1)-th column of A. */

/*           A(:,i+1) = A(:,i+1) + U1 * XA1(i+1,:)'; */
#line 743 "MB04PA.f"
	    i__2 = *n - i__;
#line 743 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &xa[i__ + 1 + xa_dim1], ldxa, &c_b10, &a[*
		    k + i__ + 1 + (i__ + 1) * a_dim1], &c__1, (ftnlen)12);
/*           A(:,i+1) = A(:,i+1) + U2 * XA2(i+1,:)'; */
#line 746 "MB04PA.f"
	    i__2 = *n - i__;
#line 746 "MB04PA.f"
	    i__3 = i__ - 1;
#line 746 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xa[i__ + 1 + nb1 * xa_dim1], ldxa, &c_b10, 
		    &a[*k + i__ + 1 + (i__ + 1) * a_dim1], &c__1, (ftnlen)12);
/*           A(:,i+1) = A(:,i+1) + YA1 * U1(i+1,:)'; */
#line 749 "MB04PA.f"
	    i__2 = *n + *k;
#line 749 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &ya[ya_offset], ldya, 
		    &qg[*k + i__ + 1 + qg_dim1], ldqg, &c_b10, &a[(i__ + 1) * 
		    a_dim1 + 1], &c__1, (ftnlen)12);
/*           A(:,i+1) = A(:,i+1) + YA2 * U2(i+1,:)'; */
#line 752 "MB04PA.f"
	    i__2 = *n + *k;
#line 752 "MB04PA.f"
	    i__3 = i__ - 1;
#line 752 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &ya[nb1 * ya_dim1 + 
		    1], ldya, &a[*k + i__ + 1 + a_dim1], lda, &c_b10, &a[(i__ 
		    + 1) * a_dim1 + 1], &c__1, (ftnlen)12);

/*           Update (i+1)-th row of A. */

#line 757 "MB04PA.f"
	    if (*n > i__ + 1) {
/*              A(i+1,i+2:n) = A(i+1,i+2:n) + U1(i+1,:)*XA1(i+2:n,:)' */
#line 759 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 759 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__, &c_b10, &xa[i__ + 2 + 
			xa_dim1], ldxa, &qg[*k + i__ + 1 + qg_dim1], ldqg, &
			c_b10, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda, (
			ftnlen)12);
/*              A(i+1,i+2:n) = A(i+1,i+2:n) + U2(i+1,:)*XA2(i+2:n,:)' */
#line 763 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 763 "MB04PA.f"
		i__3 = i__ - 1;
#line 763 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &xa[i__ + 2 + 
			nb1 * xa_dim1], ldxa, &a[*k + i__ + 1 + a_dim1], lda, 
			&c_b10, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda, (
			ftnlen)12);
/*              A(i+1,i+2:n) = A(i+1,i+2:n) + YA1(i+1,:) * U1(i+2:n,:)' */
#line 767 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 767 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &ya[*k + i__ + 1 + ya_dim1], ldya, &
			c_b10, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda, (
			ftnlen)12);
/*              A(i+1,i+2:n) = A(i+1,i+2:n) + YA2(i+1,:) * U2(i+2:n,:)' */
#line 771 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 771 "MB04PA.f"
		i__3 = i__ - 1;
#line 771 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 
			+ a_dim1], lda, &ya[*k + i__ + 1 + nb1 * ya_dim1], 
			ldya, &c_b10, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], 
			lda, (ftnlen)12);
#line 774 "MB04PA.f"
	    }

/*           Annihilate updated parts in YA. */

#line 778 "MB04PA.f"
	    i__2 = i__;
#line 778 "MB04PA.f"
	    for (j = 1; j <= i__2; ++j) {
#line 779 "MB04PA.f"
		ya[*k + i__ + 1 + j * ya_dim1] = 0.;
#line 780 "MB04PA.f"
/* L60: */
#line 780 "MB04PA.f"
	    }
#line 781 "MB04PA.f"
	    i__2 = i__ - 1;
#line 781 "MB04PA.f"
	    for (j = 1; j <= i__2; ++j) {
#line 782 "MB04PA.f"
		ya[*k + i__ + 1 + (*nb + j) * ya_dim1] = 0.;
#line 783 "MB04PA.f"
/* L70: */
#line 783 "MB04PA.f"
	    }

/*           Update XQ with first Householder reflection. */

/*           xq = Q*u1 */
#line 788 "MB04PA.f"
	    i__2 = *n - i__;
#line 788 "MB04PA.f"
	    mb01md_("Lower", &i__2, &c_b10, &qg[*k + i__ + 1 + (i__ + 1) * 
		    qg_dim1], ldqg, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &xq[i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)5);
/*           xq = xq + XQ1*w1 */
#line 791 "MB04PA.f"
	    i__2 = *n - i__;
#line 791 "MB04PA.f"
	    i__3 = i__ - 1;
#line 791 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xq[i__ + 1 + 
		    xq_dim1], ldxq, &dwork[1], &c__1, &c_b10, &xq[i__ + 1 + 
		    i__ * xq_dim1], &c__1, (ftnlen)12);
/*           xq = xq + XQ2*w2 */
#line 794 "MB04PA.f"
	    i__2 = *n - i__;
#line 794 "MB04PA.f"
	    i__3 = i__ - 1;
#line 794 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &dwork[nb1], &c__1, &c_b10, &xq[i__ + 1 + 
		    i__ * xq_dim1], &c__1, (ftnlen)12);
/*           temp = XQ1'*u1 */
#line 797 "MB04PA.f"
	    i__2 = *n - i__;
#line 797 "MB04PA.f"
	    i__3 = i__ - 1;
#line 797 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xq[i__ + 1 + xq_dim1], 
		    ldxq, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &c_b12, &
		    xq[i__ * xq_dim1 + 1], &c__1, (ftnlen)9);
/*           xq = xq - U1*temp */
#line 800 "MB04PA.f"
	    i__2 = *n - i__;
#line 800 "MB04PA.f"
	    i__3 = i__ - 1;
#line 800 "MB04PA.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b586, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &xq[i__ * xq_dim1 + 1], &c__1, &c_b10, &
		    xq[i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)12);
/*           temp = XQ2'*u1 */
#line 803 "MB04PA.f"
	    i__2 = *n - i__;
#line 803 "MB04PA.f"
	    i__3 = i__ - 1;
#line 803 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xq[i__ + 1 + nb1 * 
		    xq_dim1], ldxq, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &xq[i__ * xq_dim1 + 1], &c__1, (ftnlen)9);
/*           xq = xq - U2*temp */
#line 806 "MB04PA.f"
	    i__2 = *n - i__;
#line 806 "MB04PA.f"
	    i__3 = i__ - 1;
#line 806 "MB04PA.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b586, &a[*k + i__ + 1 + 
		    a_dim1], lda, &xq[i__ * xq_dim1 + 1], &c__1, &c_b10, &xq[
		    i__ + 1 + i__ * xq_dim1], &c__1, (ftnlen)12);
/*           xq = -tauq*xq */
#line 809 "MB04PA.f"
	    i__2 = *n - i__;
#line 809 "MB04PA.f"
	    d__1 = -tauq;
#line 809 "MB04PA.f"
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + i__ * xq_dim1], &c__1);
/*           temp = -tauq/2*xq'*u1 */
#line 811 "MB04PA.f"
	    i__2 = *n - i__;
#line 811 "MB04PA.f"
	    temp = tauq * -.5 * ddot_(&i__2, &qg[*k + i__ + 1 + i__ * qg_dim1]
		    , &c__1, &xq[i__ + 1 + i__ * xq_dim1], &c__1);
/*           xq = xq + temp*u1 */
#line 813 "MB04PA.f"
	    i__2 = *n - i__;
#line 813 "MB04PA.f"
	    daxpy_(&i__2, &temp, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &
		    xq[i__ + 1 + i__ * xq_dim1], &c__1);

/*           Update (i+1)-th column and row of Q. */

#line 817 "MB04PA.f"
	    if (*n > i__ + 1) {
/*              Q(:,i+1) = Q(:,i+1) - U1 * XQ1(i+1,:)'; */
#line 819 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 819 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__, &c_b586, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &xq[i__ + 1 + xq_dim1], ldxq, &
			c_b10, &qg[*k + i__ + 2 + (i__ + 1) * qg_dim1], &c__1,
			 (ftnlen)12);
/*              Q(:,i+1) = Q(:,i+1) - U2 * XQ2(i+1,:)'; */
#line 823 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 823 "MB04PA.f"
		i__3 = i__ - 1;
#line 823 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b586, &a[*k + i__ + 2 
			+ a_dim1], lda, &xq[i__ + 1 + nb1 * xq_dim1], ldxq, &
			c_b10, &qg[*k + i__ + 2 + (i__ + 1) * qg_dim1], &c__1,
			 (ftnlen)12);
/*              Q(:,i+1) = Q(:,i+1) + XQ1 * U1(i+1,:)'; */
#line 827 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 827 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__, &c_b10, &xq[i__ + 2 + 
			xq_dim1], ldxq, &qg[*k + i__ + 1 + qg_dim1], ldqg, &
			c_b10, &qg[*k + i__ + 2 + (i__ + 1) * qg_dim1], &c__1,
			 (ftnlen)12);
/*              Q(:,i+1) = Q(:,i+1) + XQ2 * U2(i+1,:)'; */
#line 831 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 831 "MB04PA.f"
		i__3 = i__ - 1;
#line 831 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &xq[i__ + 2 + 
			nb1 * xq_dim1], ldxq, &a[*k + i__ + 1 + a_dim1], lda, 
			&c_b10, &qg[*k + i__ + 2 + (i__ + 1) * qg_dim1], &
			c__1, (ftnlen)12);
#line 834 "MB04PA.f"
	    }

/*           Update XG with first Householder reflection. */

/*           xg = G*u1 */
#line 839 "MB04PA.f"
	    i__2 = *k + i__;
#line 839 "MB04PA.f"
	    i__3 = *n - i__;
#line 839 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &qg[(i__ + 2) * 
		    qg_dim1 + 1], ldqg, &qg[*k + i__ + 1 + i__ * qg_dim1], &
		    c__1, &c_b12, &xg[i__ * xg_dim1 + 1], &c__1, (ftnlen)12);
#line 841 "MB04PA.f"
	    i__2 = *n - i__;
#line 841 "MB04PA.f"
	    mb01md_("Upper", &i__2, &c_b10, &qg[*k + i__ + 1 + (i__ + 2) * 
		    qg_dim1], ldqg, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &xg[*k + i__ + 1 + i__ * xg_dim1], &c__1, (ftnlen)
		    5);
/*           xg = xg + XG1*w1 */
#line 844 "MB04PA.f"
	    i__2 = *k + *n;
#line 844 "MB04PA.f"
	    i__3 = i__ - 1;
#line 844 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xg[xg_offset], ldxg,
		     &dwork[1], &c__1, &c_b10, &xg[i__ * xg_dim1 + 1], &c__1, 
		    (ftnlen)12);
/*           xg = xg + XG2*w2 */
#line 847 "MB04PA.f"
	    i__2 = *k + *n;
#line 847 "MB04PA.f"
	    i__3 = i__ - 1;
#line 847 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xg[nb1 * xg_dim1 + 
		    1], ldxg, &dwork[nb1], &c__1, &c_b10, &xg[i__ * xg_dim1 + 
		    1], &c__1, (ftnlen)12);
/*           temp = XG1'*u1 */
#line 850 "MB04PA.f"
	    i__2 = *n - i__;
#line 850 "MB04PA.f"
	    i__3 = i__ - 1;
#line 850 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xg[*k + i__ + 1 + 
		    xg_dim1], ldxq, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, 
		    &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*           xg = xg - U1*temp */
#line 853 "MB04PA.f"
	    i__2 = *n - i__;
#line 853 "MB04PA.f"
	    i__3 = i__ - 1;
#line 853 "MB04PA.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b586, &qg[*k + i__ + 1 + 
		    qg_dim1], ldqg, &dwork[nb2], &c__1, &c_b10, &xg[*k + i__ 
		    + 1 + i__ * xg_dim1], &c__1, (ftnlen)12);
/*           temp = XG2'*u1 */
#line 856 "MB04PA.f"
	    i__2 = *n - i__;
#line 856 "MB04PA.f"
	    i__3 = i__ - 1;
#line 856 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &xg[*k + i__ + 1 + nb1 *
		     xg_dim1], ldxq, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1,
		     &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*           xg = xg - U2*temp */
#line 859 "MB04PA.f"
	    i__2 = *n - i__;
#line 859 "MB04PA.f"
	    i__3 = i__ - 1;
#line 859 "MB04PA.f"
	    dgemv_("No Transpose", &i__2, &i__3, &c_b586, &a[*k + i__ + 1 + 
		    a_dim1], lda, &dwork[nb2], &c__1, &c_b10, &xg[*k + i__ + 
		    1 + i__ * xg_dim1], &c__1, (ftnlen)12);
/*           xg = -tauq*xg */
#line 862 "MB04PA.f"
	    i__2 = *n + *k;
#line 862 "MB04PA.f"
	    d__1 = -tauq;
#line 862 "MB04PA.f"
	    dscal_(&i__2, &d__1, &xg[i__ * xg_dim1 + 1], &c__1);
/*           temp = -tauq/2*xq'*u1 */
#line 864 "MB04PA.f"
	    i__2 = *n - i__;
#line 864 "MB04PA.f"
	    temp = tauq * -.5 * ddot_(&i__2, &qg[*k + i__ + 1 + i__ * qg_dim1]
		    , &c__1, &xg[*k + i__ + 1 + i__ * xg_dim1], &c__1);
/*           xg = xg + temp*u1 */
#line 867 "MB04PA.f"
	    i__2 = *n - i__;
#line 867 "MB04PA.f"
	    daxpy_(&i__2, &temp, &qg[*k + i__ + 1 + i__ * qg_dim1], &c__1, &
		    xg[*k + i__ + 1 + i__ * xg_dim1], &c__1);

/*           Update (i+1)-th column and row of G. */

/*           G(:,i+1) = G(:,i+1) + XG1 * U1(i+1,:)'; */
#line 872 "MB04PA.f"
	    i__2 = *k + i__;
#line 872 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &xg[xg_offset], ldxg, 
		    &qg[*k + i__ + 1 + qg_dim1], ldqg, &c_b10, &qg[(i__ + 2) *
		     qg_dim1 + 1], &c__1, (ftnlen)12);
/*           G(:,i+1) = G(:,i+1) + XG2 * U2(i+1,:)'; */
#line 875 "MB04PA.f"
	    i__2 = *k + i__;
#line 875 "MB04PA.f"
	    i__3 = i__ - 1;
#line 875 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xg[nb1 * xg_dim1 + 
		    1], ldxg, &a[*k + i__ + 1 + a_dim1], lda, &c_b10, &qg[(
		    i__ + 2) * qg_dim1 + 1], &c__1, (ftnlen)12);
#line 877 "MB04PA.f"
	    if (*n > i__ + 1) {
/*              G(:,i+1) = G(:,i+1) + XG1 * U1(i+1,:)'; */
#line 879 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 879 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__, &c_b586, &xg[*k + i__ + 2 
			+ xg_dim1], ldxg, &qg[*k + i__ + 1 + qg_dim1], ldqg, &
			c_b10, &qg[*k + i__ + 1 + (i__ + 3) * qg_dim1], ldqg, 
			(ftnlen)12);
/*              G(:,i+1) = G(:,i+1) + XG2 * U2(i+1,:)'; */
#line 883 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 883 "MB04PA.f"
		i__3 = i__ - 1;
#line 883 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b586, &xg[*k + i__ + 
			2 + nb1 * xg_dim1], ldxg, &a[*k + i__ + 1 + a_dim1], 
			lda, &c_b10, &qg[*k + i__ + 1 + (i__ + 3) * qg_dim1], 
			ldqg, (ftnlen)12);
/*              G(:,i+1) = G(:,i+1) + U1 * XG1(i+1,:)'; */
#line 887 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 887 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &xg[*k + i__ + 1 + xg_dim1], ldxg, &
			c_b10, &qg[*k + i__ + 1 + (i__ + 3) * qg_dim1], ldqg, 
			(ftnlen)12);
/*              G(:,i+1) = G(:,i+1) + U2 * XG2(i+1,:)'; */
#line 891 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 891 "MB04PA.f"
		i__3 = i__ - 1;
#line 891 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 
			+ a_dim1], lda, &xg[*k + i__ + 1 + nb1 * xg_dim1], 
			ldxg, &c_b10, &qg[*k + i__ + 1 + (i__ + 3) * qg_dim1],
			 ldqg, (ftnlen)12);
#line 894 "MB04PA.f"
	    }

/*           Annihilate updated parts in XG. */

#line 898 "MB04PA.f"
	    i__2 = i__;
#line 898 "MB04PA.f"
	    for (j = 1; j <= i__2; ++j) {
#line 899 "MB04PA.f"
		xg[*k + i__ + 1 + j * xg_dim1] = 0.;
#line 900 "MB04PA.f"
/* L80: */
#line 900 "MB04PA.f"
	    }
#line 901 "MB04PA.f"
	    i__2 = i__ - 1;
#line 901 "MB04PA.f"
	    for (j = 1; j <= i__2; ++j) {
#line 902 "MB04PA.f"
		xg[*k + i__ + 1 + (*nb + j) * xg_dim1] = 0.;
#line 903 "MB04PA.f"
/* L90: */
#line 903 "MB04PA.f"
	    }

/*           Apply orthogonal symplectic Givens rotation. */

#line 907 "MB04PA.f"
	    i__2 = *k + i__;
#line 907 "MB04PA.f"
	    drot_(&i__2, &a[(i__ + 1) * a_dim1 + 1], &c__1, &qg[(i__ + 2) * 
		    qg_dim1 + 1], &c__1, &c__, &s);
#line 908 "MB04PA.f"
	    if (*n > i__ + 1) {
#line 909 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 909 "MB04PA.f"
		d__1 = -s;
#line 909 "MB04PA.f"
		drot_(&i__2, &a[*k + i__ + 2 + (i__ + 1) * a_dim1], &c__1, &
			qg[*k + i__ + 1 + (i__ + 3) * qg_dim1], ldqg, &c__, &
			d__1);
#line 911 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 911 "MB04PA.f"
		d__1 = -s;
#line 911 "MB04PA.f"
		drot_(&i__2, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda, &qg[*
			k + i__ + 2 + (i__ + 1) * qg_dim1], &c__1, &c__, &
			d__1);
#line 913 "MB04PA.f"
	    }
#line 914 "MB04PA.f"
	    cs[(i__ << 1) - 1] = c__;
#line 915 "MB04PA.f"
	    cs[i__ * 2] = s;
#line 916 "MB04PA.f"
	    qg[*k + i__ + 1 + i__ * qg_dim1] = tauq;

/*           Update XA with second Householder reflection. */

/*           xa = H(1:n,1:n)'*u2 */
#line 921 "MB04PA.f"
	    i__2 = *n - i__;
#line 921 "MB04PA.f"
	    i__3 = *n - i__;
#line 921 "MB04PA.f"
	    dgemv_("Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 1 + (i__ 
		    + 1) * a_dim1], lda, &a[*k + i__ + 1 + i__ * a_dim1], &
		    c__1, &c_b12, &xa[i__ + 1 + (*nb + i__) * xa_dim1], &c__1,
		     (ftnlen)9);
#line 923 "MB04PA.f"
	    if (*n > i__ + 1) {
/*              w1 = U1'*u2 */
#line 925 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 925 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 + 
			qg_dim1], ldqg, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &dwork[1], &c__1, (ftnlen)9);
/*              xa = xa + XA1*w1 */
#line 928 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 928 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__, &c_b10, &xa[i__ + 2 + 
			xa_dim1], ldxa, &dwork[1], &c__1, &c_b10, &xa[i__ + 2 
			+ (*nb + i__) * xa_dim1], &c__1, (ftnlen)12);
/*              w2 = U2'*u2 */
#line 931 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 931 "MB04PA.f"
		i__3 = i__ - 1;
#line 931 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 + 
			a_dim1], lda, &a[*k + i__ + 2 + i__ * a_dim1], &c__1, 
			&c_b12, &dwork[nb1], &c__1, (ftnlen)9);
/*              xa = xa + XA2*w2 */
#line 934 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 934 "MB04PA.f"
		i__3 = i__ - 1;
#line 934 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &xa[i__ + 2 + 
			nb1 * xa_dim1], ldxa, &dwork[nb1], &c__1, &c_b10, &xa[
			i__ + 2 + (*nb + i__) * xa_dim1], &c__1, (ftnlen)12);
/*              temp = YA1'*u2 */
#line 937 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 937 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__, &c_b10, &ya[*k + i__ + 2 + 
			ya_dim1], ldya, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &xa[(*nb + i__) * xa_dim1 + 1], &c__1, (
			ftnlen)9);
/*              xa = xa + U1*temp */
#line 940 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 940 "MB04PA.f"
		dgemv_("No Transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &xa[(*nb + i__) * xa_dim1 + 1], &
			c__1, &c_b10, &xa[i__ + 2 + (*nb + i__) * xa_dim1], &
			c__1, (ftnlen)12);
/*              temp = YA2'*u1 */
#line 943 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 943 "MB04PA.f"
		i__3 = i__ - 1;
#line 943 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &ya[*k + i__ + 2 + 
			nb1 * ya_dim1], ldya, &a[*k + i__ + 2 + i__ * a_dim1],
			 &c__1, &c_b12, &xa[(*nb + i__) * xa_dim1 + 1], &c__1,
			 (ftnlen)9);
/*              xa = xa + U2*temp */
#line 946 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 946 "MB04PA.f"
		i__3 = i__ - 1;
#line 946 "MB04PA.f"
		dgemv_("No Transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 
			+ a_dim1], lda, &xa[(*nb + i__) * xa_dim1 + 1], &c__1,
			 &c_b10, &xa[i__ + 2 + (*nb + i__) * xa_dim1], &c__1, 
			(ftnlen)12);
#line 948 "MB04PA.f"
	    }
/*           xa = -tau*xa */
#line 950 "MB04PA.f"
	    i__2 = *n - i__;
#line 950 "MB04PA.f"
	    d__1 = -tau[i__];
#line 950 "MB04PA.f"
	    dscal_(&i__2, &d__1, &xa[i__ + 1 + (*nb + i__) * xa_dim1], &c__1);

/*           Update YA with second Householder reflection. */

/*           ya = H(1:n,1:n)*u2 */
#line 955 "MB04PA.f"
	    i__2 = *k + *n;
#line 955 "MB04PA.f"
	    i__3 = *n - i__;
#line 955 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[(i__ + 1) * 
		    a_dim1 + 1], lda, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, 
		    &c_b12, &ya[(*nb + i__) * ya_dim1 + 1], &c__1, (ftnlen)12)
		    ;
#line 957 "MB04PA.f"
	    if (*n > i__ + 1) {
/*              temp = XA1'*u2 */
#line 959 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 959 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__, &c_b10, &xa[i__ + 2 + 
			xa_dim1], ldxa, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*              ya = ya + U1*temp */
#line 962 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 962 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__, &c_b10, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &dwork[nb2], &c__1, &c_b10, &ya[*k 
			+ i__ + 2 + (*nb + i__) * ya_dim1], &c__1, (ftnlen)12)
			;
/*              temp = XA2'*u1 */
#line 965 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 965 "MB04PA.f"
		i__3 = i__ - 1;
#line 965 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &xa[i__ + 2 + nb1 * 
			xa_dim1], ldxa, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*              ya = ya + U2*temp */
#line 968 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 968 "MB04PA.f"
		i__3 = i__ - 1;
#line 968 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &a[*k + i__ + 2 
			+ a_dim1], lda, &dwork[nb2], &c__1, &c_b10, &ya[*k + 
			i__ + 2 + (*nb + i__) * ya_dim1], &c__1, (ftnlen)12);
#line 970 "MB04PA.f"
	    }
/*           ya = ya + YA1*w1 */
#line 972 "MB04PA.f"
	    i__2 = *k + *n;
#line 972 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &ya[ya_offset], ldya, 
		    &dwork[1], &c__1, &c_b10, &ya[(*nb + i__) * ya_dim1 + 1], 
		    &c__1, (ftnlen)12);
/*           ya = ya + YA2*w2 */
#line 975 "MB04PA.f"
	    i__2 = *k + *n;
#line 975 "MB04PA.f"
	    i__3 = i__ - 1;
#line 975 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &ya[nb1 * ya_dim1 + 
		    1], ldya, &dwork[nb1], &c__1, &c_b10, &ya[(*nb + i__) * 
		    ya_dim1 + 1], &c__1, (ftnlen)12);
/*           ya = -tau*ya */
#line 978 "MB04PA.f"
	    i__2 = *k + *n;
#line 978 "MB04PA.f"
	    d__1 = -tau[i__];
#line 978 "MB04PA.f"
	    dscal_(&i__2, &d__1, &ya[(*nb + i__) * ya_dim1 + 1], &c__1);
/*           temp = -tau*ya'*u2 */
#line 980 "MB04PA.f"
	    i__2 = *n - i__;
#line 980 "MB04PA.f"
	    temp = -tau[i__] * ddot_(&i__2, &a[*k + i__ + 1 + i__ * a_dim1], &
		    c__1, &ya[*k + i__ + 1 + (*nb + i__) * ya_dim1], &c__1);
/*           ya = ya + temp*u2 */
#line 982 "MB04PA.f"
	    i__2 = *n - i__;
#line 982 "MB04PA.f"
	    daxpy_(&i__2, &temp, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &ya[*
		    k + i__ + 1 + (*nb + i__) * ya_dim1], &c__1);

/*           Update (i+1)-th column of A. */

/*           H(1:n,i+1) = H(1:n,i+1) + ya */
#line 987 "MB04PA.f"
	    i__2 = *k + *n;
#line 987 "MB04PA.f"
	    daxpy_(&i__2, &c_b10, &ya[(*nb + i__) * ya_dim1 + 1], &c__1, &a[(
		    i__ + 1) * a_dim1 + 1], &c__1);
/*           H(1:n,i+1) = H(1:n,i+1) + xa(i+1)*u2 */
#line 989 "MB04PA.f"
	    i__2 = *n - i__;
#line 989 "MB04PA.f"
	    daxpy_(&i__2, &xa[i__ + 1 + (*nb + i__) * xa_dim1], &a[*k + i__ + 
		    1 + i__ * a_dim1], &c__1, &a[*k + i__ + 1 + (i__ + 1) * 
		    a_dim1], &c__1);

/*           Update (i+1)-th row of A. */

#line 994 "MB04PA.f"
	    if (*n > i__ + 1) {
/*              H(i+1,i+2:n) = H(i+1,i+2:n) + xa(i+2:n)'; */
#line 996 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 996 "MB04PA.f"
		daxpy_(&i__2, &c_b10, &xa[i__ + 2 + (*nb + i__) * xa_dim1], &
			c__1, &a[*k + i__ + 1 + (i__ + 2) * a_dim1], lda);
/*              H(i+1,i+2:n) = H(i+1,i+2:n) + YA(i+1,:) * U(i+2:n,:)' */
#line 999 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 999 "MB04PA.f"
		daxpy_(&i__2, &ya[*k + i__ + 1 + (*nb + i__) * ya_dim1], &a[*
			k + i__ + 2 + i__ * a_dim1], &c__1, &a[*k + i__ + 1 + 
			(i__ + 2) * a_dim1], lda);
#line 1001 "MB04PA.f"
	    }

/*           Annihilate updated parts in YA. */

#line 1005 "MB04PA.f"
	    ya[*k + i__ + 1 + (*nb + i__) * ya_dim1] = 0.;

/*           Update XQ with second Householder reflection. */

/*           xq = Q*u2 */
#line 1010 "MB04PA.f"
	    i__2 = *n - i__;
#line 1010 "MB04PA.f"
	    mb01md_("Lower", &i__2, &c_b10, &qg[*k + i__ + 1 + (i__ + 1) * 
		    qg_dim1], ldqg, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b12, &xq[i__ + 1 + (*nb + i__) * xq_dim1], &c__1, (
		    ftnlen)5);
#line 1012 "MB04PA.f"
	    if (*n > i__ + 1) {
/*              xq = xq + XQ1*w1 */
#line 1014 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 1014 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__, &c_b10, &xq[i__ + 2 + 
			xq_dim1], ldxq, &dwork[1], &c__1, &c_b10, &xq[i__ + 2 
			+ (*nb + i__) * xq_dim1], &c__1, (ftnlen)12);
/*              xq = xq + XQ2*w2 */
#line 1017 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 1017 "MB04PA.f"
		i__3 = i__ - 1;
#line 1017 "MB04PA.f"
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &xq[i__ + 2 + 
			nb1 * xq_dim1], ldxq, &dwork[nb1], &c__1, &c_b10, &xq[
			i__ + 2 + (*nb + i__) * xq_dim1], &c__1, (ftnlen)12);
/*              temp = XQ1'*u2 */
#line 1020 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 1020 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__, &c_b10, &xq[i__ + 2 + 
			xq_dim1], ldxq, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &xq[(*nb + i__) * xq_dim1 + 1], &c__1, (
			ftnlen)9);
/*              xq = xq - U1*temp */
#line 1023 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 1023 "MB04PA.f"
		dgemv_("No Transpose", &i__2, &i__, &c_b586, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &xq[(*nb + i__) * xq_dim1 + 1], &
			c__1, &c_b10, &xq[i__ + 2 + (*nb + i__) * xq_dim1], &
			c__1, (ftnlen)12);
/*              temp = XQ2'*u2 */
#line 1026 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 1026 "MB04PA.f"
		i__3 = i__ - 1;
#line 1026 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &xq[i__ + 2 + nb1 * 
			xq_dim1], ldxq, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &xq[(*nb + i__) * xq_dim1 + 1], &c__1, (
			ftnlen)9);
/*              xq = xq - U2*temp */
#line 1029 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 1029 "MB04PA.f"
		i__3 = i__ - 1;
#line 1029 "MB04PA.f"
		dgemv_("No Transpose", &i__2, &i__3, &c_b586, &a[*k + i__ + 2 
			+ a_dim1], lda, &xq[(*nb + i__) * xq_dim1 + 1], &c__1,
			 &c_b10, &xq[i__ + 2 + (*nb + i__) * xq_dim1], &c__1, 
			(ftnlen)12);
#line 1031 "MB04PA.f"
	    }
/*           xq = -tauq*xq */
#line 1033 "MB04PA.f"
	    i__2 = *n - i__;
#line 1033 "MB04PA.f"
	    d__1 = -tau[i__];
#line 1033 "MB04PA.f"
	    dscal_(&i__2, &d__1, &xq[i__ + 1 + (*nb + i__) * xq_dim1], &c__1);
/*           temp = -tauq/2*xq'*u2 */
#line 1035 "MB04PA.f"
	    i__2 = *n - i__;
#line 1035 "MB04PA.f"
	    temp = tau[i__] * -.5 * ddot_(&i__2, &a[*k + i__ + 1 + i__ * 
		    a_dim1], &c__1, &xq[i__ + 1 + (*nb + i__) * xq_dim1], &
		    c__1);
/*           xq = xq + temp*u2 */
#line 1038 "MB04PA.f"
	    i__2 = *n - i__;
#line 1038 "MB04PA.f"
	    daxpy_(&i__2, &temp, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &xq[
		    i__ + 1 + (*nb + i__) * xq_dim1], &c__1);

/*           Update (i+1)-th column and row of Q. */

#line 1042 "MB04PA.f"
	    if (*n > i__ + 1) {
#line 1043 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 1043 "MB04PA.f"
		daxpy_(&i__2, &c_b10, &xq[i__ + 2 + (*nb + i__) * xq_dim1], &
			c__1, &qg[*k + i__ + 2 + (i__ + 1) * qg_dim1], &c__1);
/*              H(1:n,n+i+1) = H(1:n,n+i+1) - U * XQ(i+1,:)'; */
#line 1046 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 1046 "MB04PA.f"
		d__1 = -xq[i__ + 1 + (*nb + i__) * xq_dim1];
#line 1046 "MB04PA.f"
		daxpy_(&i__2, &d__1, &a[*k + i__ + 2 + i__ * a_dim1], &c__1, &
			qg[*k + i__ + 2 + (i__ + 1) * qg_dim1], &c__1);
#line 1048 "MB04PA.f"
	    }

/*           Update XG with second Householder reflection. */

/*           xg = G*u2 */
#line 1053 "MB04PA.f"
	    i__2 = *k + i__;
#line 1053 "MB04PA.f"
	    i__3 = *n - i__;
#line 1053 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &qg[(i__ + 2) * 
		    qg_dim1 + 1], ldqg, &a[*k + i__ + 1 + i__ * a_dim1], &
		    c__1, &c_b12, &xg[(*nb + i__) * xg_dim1 + 1], &c__1, (
		    ftnlen)12);
#line 1055 "MB04PA.f"
	    i__2 = *n - i__;
#line 1055 "MB04PA.f"
	    mb01md_("Upper", &i__2, &c_b10, &qg[*k + i__ + 1 + (i__ + 2) * 
		    qg_dim1], ldqg, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &
		    c_b12, &xg[*k + i__ + 1 + (*nb + i__) * xg_dim1], &c__1, (
		    ftnlen)5);
/*           xg = xg + XG1*w1 */
#line 1058 "MB04PA.f"
	    i__2 = *k + *n;
#line 1058 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__, &c_b10, &xg[xg_offset], ldxg, 
		    &dwork[1], &c__1, &c_b10, &xg[(*nb + i__) * xg_dim1 + 1], 
		    &c__1, (ftnlen)12);
/*           xg = xg + XG2*w2 */
#line 1061 "MB04PA.f"
	    i__2 = *k + *n;
#line 1061 "MB04PA.f"
	    i__3 = i__ - 1;
#line 1061 "MB04PA.f"
	    dgemv_("No transpose", &i__2, &i__3, &c_b10, &xg[nb1 * xg_dim1 + 
		    1], ldxg, &dwork[nb1], &c__1, &c_b10, &xg[(*nb + i__) * 
		    xg_dim1 + 1], &c__1, (ftnlen)12);
#line 1063 "MB04PA.f"
	    if (*n > i__ + 1) {
/*              temp = XG1'*u2 */
#line 1065 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 1065 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__, &c_b10, &xg[*k + i__ + 2 + 
			xg_dim1], ldxq, &a[*k + i__ + 2 + i__ * a_dim1], &
			c__1, &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*              xg = xg - U1*temp */
#line 1068 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 1068 "MB04PA.f"
		dgemv_("No Transpose", &i__2, &i__, &c_b586, &qg[*k + i__ + 2 
			+ qg_dim1], ldqg, &dwork[nb2], &c__1, &c_b10, &xg[*k 
			+ i__ + 2 + (*nb + i__) * xg_dim1], &c__1, (ftnlen)12)
			;
/*              temp = XG2'*u2 */
#line 1071 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 1071 "MB04PA.f"
		i__3 = i__ - 1;
#line 1071 "MB04PA.f"
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &xg[*k + i__ + 2 + 
			nb1 * xg_dim1], ldxq, &a[*k + i__ + 2 + i__ * a_dim1],
			 &c__1, &c_b12, &dwork[nb2], &c__1, (ftnlen)9);
/*              xg = xg - U2*temp */
#line 1074 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 1074 "MB04PA.f"
		i__3 = i__ - 1;
#line 1074 "MB04PA.f"
		dgemv_("No Transpose", &i__2, &i__3, &c_b586, &a[*k + i__ + 2 
			+ a_dim1], lda, &dwork[nb2], &c__1, &c_b10, &xg[*k + 
			i__ + 2 + (*nb + i__) * xg_dim1], &c__1, (ftnlen)12);
#line 1076 "MB04PA.f"
	    }
/*           xg = -tauq*xg */
#line 1078 "MB04PA.f"
	    i__2 = *n + *k;
#line 1078 "MB04PA.f"
	    d__1 = -tau[i__];
#line 1078 "MB04PA.f"
	    dscal_(&i__2, &d__1, &xg[(*nb + i__) * xg_dim1 + 1], &c__1);
/*           temp = -tauq/2*xg'*u1 */
#line 1080 "MB04PA.f"
	    i__2 = *n - i__;
#line 1080 "MB04PA.f"
	    temp = tau[i__] * -.5 * ddot_(&i__2, &a[*k + i__ + 1 + i__ * 
		    a_dim1], &c__1, &xg[*k + i__ + 1 + (*nb + i__) * xg_dim1],
		     &c__1);
/*           xg = xg + temp*u1 */
#line 1083 "MB04PA.f"
	    i__2 = *n - i__;
#line 1083 "MB04PA.f"
	    daxpy_(&i__2, &temp, &a[*k + i__ + 1 + i__ * a_dim1], &c__1, &xg[*
		    k + i__ + 1 + (*nb + i__) * xg_dim1], &c__1);

/*           Update (i+1)-th column and row of G. */

#line 1087 "MB04PA.f"
	    i__2 = *k + i__;
#line 1087 "MB04PA.f"
	    daxpy_(&i__2, &c_b10, &xg[(*nb + i__) * xg_dim1 + 1], &c__1, &qg[(
		    i__ + 2) * qg_dim1 + 1], &c__1);
#line 1088 "MB04PA.f"
	    if (*n > i__ + 1) {
#line 1089 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 1089 "MB04PA.f"
		daxpy_(&i__2, &c_b586, &xg[*k + i__ + 2 + (*nb + i__) * 
			xg_dim1], &c__1, &qg[*k + i__ + 1 + (i__ + 3) * 
			qg_dim1], ldqg);
#line 1091 "MB04PA.f"
		i__2 = *n - i__ - 1;
#line 1091 "MB04PA.f"
		daxpy_(&i__2, &xg[*k + i__ + 1 + (*nb + i__) * xg_dim1], &a[*
			k + i__ + 2 + i__ * a_dim1], &c__1, &qg[*k + i__ + 1 
			+ (i__ + 3) * qg_dim1], ldqg);
#line 1093 "MB04PA.f"
	    }

/*           Annihilate updated parts in XG. */

#line 1097 "MB04PA.f"
	    xg[*k + i__ + 1 + (*nb + i__) * xg_dim1] = 0.;

#line 1099 "MB04PA.f"
	    a[*k + i__ + 1 + i__ * a_dim1] = aki;
#line 1100 "MB04PA.f"
/* L100: */
#line 1100 "MB04PA.f"
	}
#line 1101 "MB04PA.f"
    }

#line 1103 "MB04PA.f"
    return 0;
/* *** Last line of MB04PA *** */
} /* mb04pa_ */

