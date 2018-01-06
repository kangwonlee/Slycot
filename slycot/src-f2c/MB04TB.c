#line 1 "MB04TB.f"
/* MB04TB.f -- translated by f2c (version 20100827).
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

#line 1 "MB04TB.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__3 = 3;
static doublereal c_b25 = 1.;

/* Subroutine */ int mb04tb_(char *trana, char *tranb, integer *n, integer *
	ilo, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *g, integer *ldg, doublereal *q, integer *ldq, doublereal *
	csl, doublereal *csr, doublereal *taul, doublereal *taur, doublereal *
	dwork, integer *ldwork, integer *info, ftnlen trana_len, ftnlen 
	tranb_len)
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, b_dim1, b_offset, g_dim1, g_offset, q_dim1, 
	    q_offset, i__1, i__2[2], i__3, i__4, i__5;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer i__, ib, nb, nh, nx, nib, nnb, pxa, pxb, pdw, pya, pyb, 
	    pxg, pyg, pxq, pyq, ierr;
    static logical ltra, ltrb;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern integer ue01md_(integer *, char *, char *, integer *, integer *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin;
    extern /* Subroutine */ int mb04ts_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), mb03xu_(logical *, logical *, integer *, integer 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), xerbla_(char *, integer *, ftnlen);
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

/*     To compute a symplectic URV (SURV) decomposition of a real */
/*     2N-by-2N matrix H, */

/*            [ op(A)   G   ]                 [ op(R11)   R12   ] */
/*        H = [             ] = U R V'  = U * [                 ] * V' , */
/*            [  Q    op(B) ]                 [   0     op(R22) ] */

/*     where A, B, G, Q, R12 are real N-by-N matrices, op(R11) is a real */
/*     N-by-N upper triangular matrix, op(R22) is a real N-by-N lower */
/*     Hessenberg matrix and U, V are 2N-by-2N orthogonal symplectic */
/*     matrices. Blocked version. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op( A ) as follows: */
/*             = 'N':  op( A ) = A; */
/*             = 'T':  op( A ) = A'; */
/*             = 'C':  op( A ) = A'. */

/*     TRANB   CHARACTER*1 */
/*             Specifies the form of op( B ) as follows: */
/*             = 'N':  op( B ) = B; */
/*             = 'T':  op( B ) = B'; */
/*             = 'C':  op( B ) = B'. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A. N >= 0. */

/*     ILO     (input) INTEGER */
/*             It is assumed that op(A) is already upper triangular, */
/*             op(B) is lower triangular and Q is zero in rows and */
/*             columns 1:ILO-1. ILO is normally set by a previous call */
/*             to MB04DD; otherwise it should be set to 1. */
/*             1 <= ILO <= N, if N > 0; ILO=1, if N=0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the triangular matrix R11, and in the zero part */
/*             information about the elementary reflectors used to */
/*             compute the SURV decomposition. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix B. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the Hessenberg matrix R22, and in the zero part */
/*             information about the elementary reflectors used to */
/*             compute the SURV decomposition. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     G       (input/output) DOUBLE PRECISION array, dimension (LDG,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix G. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the matrix R12. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G.  LDG >= MAX(1,N). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix Q. */
/*             On exit, the leading N-by-N part of this array contains */
/*             information about the elementary reflectors used to */
/*             compute the SURV decomposition. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q.  LDQ >= MAX(1,N). */

/*     CSL     (output) DOUBLE PRECISION array, dimension (2N) */
/*             On exit, the first 2N elements of this array contain the */
/*             cosines and sines of the symplectic Givens rotations */
/*             applied from the left-hand side used to compute the SURV */
/*             decomposition. */

/*     CSR     (output) DOUBLE PRECISION array, dimension (2N-2) */
/*             On exit, the first 2N-2 elements of this array contain the */
/*             cosines and sines of the symplectic Givens rotations */
/*             applied from the right-hand side used to compute the SURV */
/*             decomposition. */

/*     TAUL    (output) DOUBLE PRECISION array, dimension (N) */
/*             On exit, the first N elements of this array contain the */
/*             scalar factors of some of the elementary reflectors */
/*             applied form the left-hand side. */

/*     TAUR    (output) DOUBLE PRECISION array, dimension (N-1) */
/*             On exit, the first N-1 elements of this array contain the */
/*             scalar factors of some of the elementary reflectors */
/*             applied form the right-hand side. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK, (16*N + 5)*NB, where NB is the optimal */
/*             block size determined by the function UE01MD. */
/*             On exit, if  INFO = -16,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The matrices U and V are represented as products of symplectic */
/*     reflectors and Givens rotators */

/*     U = diag( HU(1),HU(1) )  GU(1)  diag( FU(1),FU(1) ) */
/*         diag( HU(2),HU(2) )  GU(2)  diag( FU(2),FU(2) ) */
/*                              .... */
/*         diag( HU(n),HU(n) )  GU(n)  diag( FU(n),FU(n) ), */

/*     V = diag( HV(1),HV(1) )       GV(1)   diag( FV(1),FV(1) ) */
/*         diag( HV(2),HV(2) )       GV(2)   diag( FV(2),FV(2) ) */
/*                                   .... */
/*         diag( HV(n-1),HV(n-1) )  GV(n-1)  diag( FV(n-1),FV(n-1) ). */

/*     Each HU(i) has the form */

/*           HU(i) = I - tau * v * v' */

/*     where tau is a real scalar, and v is a real vector with */
/*     v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in */
/*     Q(i+1:n,i), and tau in Q(i,i). */

/*     Each FU(i) has the form */

/*           FU(i) = I - nu * w * w' */

/*     where nu is a real scalar, and w is a real vector with */
/*     w(1:i-1) = 0 and w(i) = 1; w(i+1:n) is stored on exit in */
/*     A(i+1:n,i), if op(A) = 'N', and in A(i,i+1:n), otherwise. The */
/*     scalar nu is stored in TAUL(i). */

/*     Each GU(i) is a Givens rotator acting on rows i and n+i, */
/*     where the cosine is stored in CSL(2*i-1) and the sine in */
/*     CSL(2*i). */

/*     Each HV(i) has the form */

/*           HV(i) = I - tau * v * v' */

/*     where tau is a real scalar, and v is a real vector with */
/*     v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in */
/*     Q(i,i+2:n), and tau in Q(i,i+1). */

/*     Each FV(i) has the form */

/*           FV(i) = I - nu * w * w' */

/*     where nu is a real scalar, and w is a real vector with */
/*     w(1:i) = 0 and w(i+1) = 1; w(i+2:n) is stored on exit in */
/*     B(i,i+2:n), if op(B) = 'N', and in B(i+2:n,i), otherwise. */
/*     The scalar nu is stored in TAUR(i). */

/*     Each GV(i) is a Givens rotator acting on columns i+1 and n+i+1, */
/*     where the cosine is stored in CSR(2*i-1) and the sine in */
/*     CSR(2*i). */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 80/3*N**3 + ( 64*NB + 77 )*N**2 + */
/*     ( -16*NB + 48 )*NB*N + O(N) floating point operations, where */
/*     NB is the used block size, and is numerically backward stable. */

/*     REFERENCES */

/*     [1] Benner, P., Mehrmann, V., and Xu, H. */
/*         A numerically stable, structure preserving method for */
/*         computing the eigenvalues of real Hamiltonian or symplectic */
/*         pencils. Numer. Math., Vol 78 (3), pp. 329-358, 1998. */

/*     [2] Kressner, D. */
/*         Block algorithms for orthogonal symplectic factorizations. */
/*         BIT, 43 (4), pp. 775-790, 2003. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DGESUB). */

/*     KEYWORDS */

/*     Elementary matrix operations, Matrix decompositions, Hamiltonian */
/*     matrix */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

#line 263 "MB04TB.f"
    /* Parameter adjustments */
#line 263 "MB04TB.f"
    a_dim1 = *lda;
#line 263 "MB04TB.f"
    a_offset = 1 + a_dim1;
#line 263 "MB04TB.f"
    a -= a_offset;
#line 263 "MB04TB.f"
    b_dim1 = *ldb;
#line 263 "MB04TB.f"
    b_offset = 1 + b_dim1;
#line 263 "MB04TB.f"
    b -= b_offset;
#line 263 "MB04TB.f"
    g_dim1 = *ldg;
#line 263 "MB04TB.f"
    g_offset = 1 + g_dim1;
#line 263 "MB04TB.f"
    g -= g_offset;
#line 263 "MB04TB.f"
    q_dim1 = *ldq;
#line 263 "MB04TB.f"
    q_offset = 1 + q_dim1;
#line 263 "MB04TB.f"
    q -= q_offset;
#line 263 "MB04TB.f"
    --csl;
#line 263 "MB04TB.f"
    --csr;
#line 263 "MB04TB.f"
    --taul;
#line 263 "MB04TB.f"
    --taur;
#line 263 "MB04TB.f"
    --dwork;
#line 263 "MB04TB.f"

#line 263 "MB04TB.f"
    /* Function Body */
#line 263 "MB04TB.f"
    *info = 0;
#line 264 "MB04TB.f"
    ltra = lsame_(trana, "T", (ftnlen)1, (ftnlen)1) || lsame_(trana, "C", (
	    ftnlen)1, (ftnlen)1);
#line 265 "MB04TB.f"
    ltrb = lsame_(tranb, "T", (ftnlen)1, (ftnlen)1) || lsame_(tranb, "C", (
	    ftnlen)1, (ftnlen)1);
#line 266 "MB04TB.f"
    if (! ltra && ! lsame_(trana, "N", (ftnlen)1, (ftnlen)1)) {
#line 267 "MB04TB.f"
	*info = -1;
#line 268 "MB04TB.f"
    } else if (! ltrb && ! lsame_(tranb, "N", (ftnlen)1, (ftnlen)1)) {
#line 269 "MB04TB.f"
	*info = -2;
#line 270 "MB04TB.f"
    } else if (*n < 0) {
#line 271 "MB04TB.f"
	*info = -3;
#line 272 "MB04TB.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 273 "MB04TB.f"
	*info = -4;
#line 274 "MB04TB.f"
    } else if (*lda < max(1,*n)) {
#line 275 "MB04TB.f"
	*info = -6;
#line 276 "MB04TB.f"
    } else if (*ldb < max(1,*n)) {
#line 277 "MB04TB.f"
	*info = -8;
#line 278 "MB04TB.f"
    } else if (*ldg < max(1,*n)) {
#line 279 "MB04TB.f"
	*info = -10;
#line 280 "MB04TB.f"
    } else if (*ldq < max(1,*n)) {
#line 281 "MB04TB.f"
	*info = -12;
#line 282 "MB04TB.f"
    } else if (*ldwork < max(1,*n)) {
#line 283 "MB04TB.f"
	dwork[1] = (doublereal) max(1,*n);
#line 284 "MB04TB.f"
	*info = -18;
#line 285 "MB04TB.f"
    }

/*     Return if there were illegal values. */

#line 289 "MB04TB.f"
    if (*info != 0) {
#line 290 "MB04TB.f"
	i__1 = -(*info);
#line 290 "MB04TB.f"
	xerbla_("MB04TB", &i__1, (ftnlen)6);
#line 291 "MB04TB.f"
	return 0;
#line 292 "MB04TB.f"
    }

/*     Set elements 1:ILO-1 of CSL, CSR, TAUL and TAUR to their default */
/*     values. */

#line 297 "MB04TB.f"
    i__1 = *ilo - 1;
#line 297 "MB04TB.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 298 "MB04TB.f"
	csl[(i__ << 1) - 1] = 1.;
#line 299 "MB04TB.f"
	csl[i__ * 2] = 0.;
#line 300 "MB04TB.f"
	csr[(i__ << 1) - 1] = 1.;
#line 301 "MB04TB.f"
	csr[i__ * 2] = 0.;
#line 302 "MB04TB.f"
	taul[i__] = 0.;
#line 303 "MB04TB.f"
	taur[i__] = 0.;
#line 304 "MB04TB.f"
/* L10: */
#line 304 "MB04TB.f"
    }

/*     Quick return if possible. */

#line 308 "MB04TB.f"
    nh = *n - *ilo + 1;
#line 309 "MB04TB.f"
    if (nh == 0) {
#line 310 "MB04TB.f"
	dwork[1] = 1.;
#line 311 "MB04TB.f"
	return 0;
#line 312 "MB04TB.f"
    }

/*     Determine the block size. */

/* Writing concatenation */
#line 316 "MB04TB.f"
    i__2[0] = 1, a__1[0] = trana;
#line 316 "MB04TB.f"
    i__2[1] = 1, a__1[1] = tranb;
#line 316 "MB04TB.f"
    s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)2);
#line 316 "MB04TB.f"
    nb = ue01md_(&c__1, "MB04TB", ch__1, n, ilo, &c_n1, (ftnlen)6, (ftnlen)2);
#line 317 "MB04TB.f"
    nbmin = 2;
#line 318 "MB04TB.f"
    wrkopt = *n;
#line 319 "MB04TB.f"
    if (nb > 1 && nb < nh) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
/* Writing concatenation */
#line 323 "MB04TB.f"
	i__2[0] = 1, a__1[0] = trana;
#line 323 "MB04TB.f"
	i__2[1] = 1, a__1[1] = tranb;
#line 323 "MB04TB.f"
	s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)2);
#line 323 "MB04TB.f"
	i__1 = nb, i__3 = ue01md_(&c__3, "MB04TB", ch__1, n, ilo, &c_n1, (
		ftnlen)6, (ftnlen)2);
#line 323 "MB04TB.f"
	nx = max(i__1,i__3);
#line 325 "MB04TB.f"
	if (nx < nh) {

/*           Check whether workspace is large enough for blocked code. */

#line 329 "MB04TB.f"
	    wrkopt = (*n << 4) * nb + nb * 5;
#line 330 "MB04TB.f"
	    if (*ldwork < wrkopt) {

/*              Not enough workspace available. Determine minimum value */
/*              of NB, and reduce NB. */

/* Computing MAX */
/* Writing concatenation */
#line 335 "MB04TB.f"
		i__2[0] = 1, a__1[0] = trana;
#line 335 "MB04TB.f"
		i__2[1] = 1, a__1[1] = tranb;
#line 335 "MB04TB.f"
		s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)2);
#line 335 "MB04TB.f"
		i__1 = 2, i__3 = ue01md_(&c__2, "MB04TB", ch__1, n, ilo, &
			c_n1, (ftnlen)6, (ftnlen)2);
#line 335 "MB04TB.f"
		nbmin = max(i__1,i__3);
#line 337 "MB04TB.f"
		nb = *ldwork / ((*n << 4) + 5);
#line 338 "MB04TB.f"
	    }
#line 339 "MB04TB.f"
	}
#line 340 "MB04TB.f"
    }

#line 342 "MB04TB.f"
    nnb = *n * nb;
#line 343 "MB04TB.f"
    pyb = 1;
#line 344 "MB04TB.f"
    pyq = pyb + (nnb << 1);
#line 345 "MB04TB.f"
    pya = pyq + (nnb << 1);
#line 346 "MB04TB.f"
    pyg = pya + (nnb << 1);
#line 347 "MB04TB.f"
    pxq = pyg + (nnb << 1);
#line 348 "MB04TB.f"
    pxa = pxq + (nnb << 1);
#line 349 "MB04TB.f"
    pxg = pxa + (nnb << 1);
#line 350 "MB04TB.f"
    pxb = pxg + (nnb << 1);
#line 351 "MB04TB.f"
    pdw = pxb + (nnb << 1);

#line 353 "MB04TB.f"
    if (nb < nbmin || nb >= nh) {

/*        Use unblocked code. */

#line 357 "MB04TB.f"
	i__ = *ilo;

#line 359 "MB04TB.f"
    } else if (ltra && ltrb) {
#line 360 "MB04TB.f"
	i__1 = *n - nx - 1;
#line 360 "MB04TB.f"
	i__3 = nb;
#line 360 "MB04TB.f"
	for (i__ = *ilo; i__3 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__3) {
/* Computing MIN */
#line 361 "MB04TB.f"
	    i__4 = nb, i__5 = *n - i__;
#line 361 "MB04TB.f"
	    ib = min(i__4,i__5);
#line 362 "MB04TB.f"
	    nib = *n * ib;

/*           Reduce rows and columns i:i+nb-1 to symplectic URV form and */
/*           return the matrices XA, XB, XG, XQ, YA, YB, YG and YQ which */
/*           are needed to update the unreduced parts of the matrices. */

#line 368 "MB04TB.f"
	    i__4 = *n - i__ + 1;
#line 368 "MB04TB.f"
	    i__5 = i__ - 1;
#line 368 "MB04TB.f"
	    mb03xu_(&ltra, &ltrb, &i__4, &i__5, &ib, &a[i__ + a_dim1], lda, &
		    b[i__ * b_dim1 + 1], ldb, &g[g_offset], ldg, &q[i__ + i__ 
		    * q_dim1], ldq, &dwork[pxa], n, &dwork[pxb], n, &dwork[
		    pxg], n, &dwork[pxq], n, &dwork[pya], n, &dwork[pyb], n, &
		    dwork[pyg], n, &dwork[pyq], n, &csl[(i__ << 1) - 1], &csr[
		    (i__ << 1) - 1], &taul[i__], &taur[i__], &dwork[pdw]);

/*           Update the submatrix A(i+1+ib:n,1:n). */

#line 377 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 377 "MB04TB.f"
	    i__5 = *n - i__ - ib + 1;
#line 377 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &
		    dwork[pxa + nb + 1], n, &q[i__ + ib + i__ * q_dim1], ldq, 
		    &c_b25, &a[i__ + ib + 1 + (i__ + ib) * a_dim1], lda, (
		    ftnlen)12, (ftnlen)9);
#line 380 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 380 "MB04TB.f"
	    i__5 = *n - i__ - ib + 1;
#line 380 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", &i__4, &i__5, &ib, &c_b25, 
		    &dwork[pxa + nib + nb + 1], n, &a[i__ + (i__ + ib) * 
		    a_dim1], lda, &c_b25, &a[i__ + ib + 1 + (i__ + ib) * 
		    a_dim1], lda, (ftnlen)12, (ftnlen)12);
#line 383 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 383 "MB04TB.f"
	    dgemm_("Transpose", "Transpose", &i__4, n, &ib, &c_b25, &q[i__ + (
		    i__ + ib + 1) * q_dim1], ldq, &dwork[pya], n, &c_b25, &a[
		    i__ + ib + 1 + a_dim1], lda, (ftnlen)9, (ftnlen)9);
#line 386 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 386 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, n, &ib, &c_b25, &b[i__ 
		    + ib + 1 + i__ * b_dim1], ldb, &dwork[pya + nib], n, &
		    c_b25, &a[i__ + ib + 1 + a_dim1], lda, (ftnlen)12, (
		    ftnlen)9);

/*           Update the submatrix Q(i+ib:n,i+1+ib:n). */

#line 392 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 392 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 392 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &q[
		    i__ + ib + i__ * q_dim1], ldq, &dwork[pxq + nb + 1], n, &
		    c_b25, &q[i__ + ib + (i__ + ib + 1) * q_dim1], ldq, (
		    ftnlen)12, (ftnlen)9);
#line 395 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 395 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 395 "MB04TB.f"
	    dgemm_("Transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &a[
		    i__ + (i__ + ib) * a_dim1], lda, &dwork[pxq + nib + nb + 
		    1], n, &c_b25, &q[i__ + ib + (i__ + ib + 1) * q_dim1], 
		    ldq, (ftnlen)9, (ftnlen)9);
#line 398 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 398 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 398 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", &i__4, &i__5, &ib, &c_b25, 
		    &dwork[pyq + nb], n, &q[i__ + (i__ + ib + 1) * q_dim1], 
		    ldq, &c_b25, &q[i__ + ib + (i__ + ib + 1) * q_dim1], ldq, 
		    (ftnlen)12, (ftnlen)12);
#line 401 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 401 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 401 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &
		    dwork[pyq + nib + nb], n, &b[i__ + ib + 1 + i__ * b_dim1],
		     ldb, &c_b25, &q[i__ + ib + (i__ + ib + 1) * q_dim1], ldq,
		     (ftnlen)12, (ftnlen)9);

/*           Update the matrix G. */

#line 407 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 407 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, n, &ib, &c_b25, &q[i__ 
		    + ib + i__ * q_dim1], ldq, &dwork[pxg], n, &c_b25, &g[i__ 
		    + ib + g_dim1], ldg, (ftnlen)12, (ftnlen)9);
#line 410 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 410 "MB04TB.f"
	    dgemm_("Transpose", "Transpose", &i__4, n, &ib, &c_b25, &a[i__ + (
		    i__ + ib) * a_dim1], lda, &dwork[pxg + nib], n, &c_b25, &
		    g[i__ + ib + g_dim1], ldg, (ftnlen)9, (ftnlen)9);
#line 413 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 413 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", n, &i__4, &ib, &c_b25, &
		    dwork[pyg], n, &q[i__ + (i__ + ib + 1) * q_dim1], ldq, &
		    c_b25, &g[(i__ + ib + 1) * g_dim1 + 1], ldg, (ftnlen)12, (
		    ftnlen)12);
#line 416 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 416 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", n, &i__4, &ib, &c_b25, &dwork[
		    pyg + nib], n, &b[i__ + ib + 1 + i__ * b_dim1], ldb, &
		    c_b25, &g[(i__ + ib + 1) * g_dim1 + 1], ldg, (ftnlen)12, (
		    ftnlen)9);

/*           Update the submatrix B(1:n,i+ib:n). */

#line 422 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 422 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", n, &i__4, &ib, &c_b25, &dwork[
		    pxb], n, &q[i__ + ib + i__ * q_dim1], ldq, &c_b25, &b[(
		    i__ + ib) * b_dim1 + 1], ldb, (ftnlen)12, (ftnlen)9);
#line 425 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 425 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", n, &i__4, &ib, &c_b25, &
		    dwork[pxb + nib], n, &a[i__ + (i__ + ib) * a_dim1], lda, &
		    c_b25, &b[(i__ + ib) * b_dim1 + 1], ldb, (ftnlen)12, (
		    ftnlen)12);
#line 428 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 428 "MB04TB.f"
	    i__5 = *n - i__ - ib + 1;
#line 428 "MB04TB.f"
	    dgemm_("Transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &q[
		    i__ + (i__ + ib + 1) * q_dim1], ldq, &dwork[pyb + nb], n, 
		    &c_b25, &b[i__ + ib + 1 + (i__ + ib) * b_dim1], ldb, (
		    ftnlen)9, (ftnlen)9);
#line 431 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 431 "MB04TB.f"
	    i__5 = *n - i__ - ib + 1;
#line 431 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &b[
		    i__ + ib + 1 + i__ * b_dim1], ldb, &dwork[pyb + nib + nb],
		     n, &c_b25, &b[i__ + ib + 1 + (i__ + ib) * b_dim1], ldb, (
		    ftnlen)12, (ftnlen)9);
#line 434 "MB04TB.f"
/* L20: */
#line 434 "MB04TB.f"
	}

#line 436 "MB04TB.f"
    } else if (ltra) {
#line 437 "MB04TB.f"
	i__3 = *n - nx - 1;
#line 437 "MB04TB.f"
	i__1 = nb;
#line 437 "MB04TB.f"
	for (i__ = *ilo; i__1 < 0 ? i__ >= i__3 : i__ <= i__3; i__ += i__1) {
/* Computing MIN */
#line 438 "MB04TB.f"
	    i__4 = nb, i__5 = *n - i__;
#line 438 "MB04TB.f"
	    ib = min(i__4,i__5);
#line 439 "MB04TB.f"
	    nib = *n * ib;

/*           Reduce rows and columns i:i+nb-1 to symplectic URV form and */
/*           return the matrices XA, XB, XG, XQ, YA, YB, YG and YQ which */
/*           are needed to update the unreduced parts of the matrices. */

#line 445 "MB04TB.f"
	    i__4 = *n - i__ + 1;
#line 445 "MB04TB.f"
	    i__5 = i__ - 1;
#line 445 "MB04TB.f"
	    mb03xu_(&ltra, &ltrb, &i__4, &i__5, &ib, &a[i__ + a_dim1], lda, &
		    b[i__ + b_dim1], ldb, &g[g_offset], ldg, &q[i__ + i__ * 
		    q_dim1], ldq, &dwork[pxa], n, &dwork[pxb], n, &dwork[pxg],
		     n, &dwork[pxq], n, &dwork[pya], n, &dwork[pyb], n, &
		    dwork[pyg], n, &dwork[pyq], n, &csl[(i__ << 1) - 1], &csr[
		    (i__ << 1) - 1], &taul[i__], &taur[i__], &dwork[pdw]);

/*           Update the submatrix A(i+1+ib:n,1:n). */

#line 454 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 454 "MB04TB.f"
	    i__5 = *n - i__ - ib + 1;
#line 454 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &
		    dwork[pxa + nb + 1], n, &q[i__ + ib + i__ * q_dim1], ldq, 
		    &c_b25, &a[i__ + ib + 1 + (i__ + ib) * a_dim1], lda, (
		    ftnlen)12, (ftnlen)9);
#line 457 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 457 "MB04TB.f"
	    i__5 = *n - i__ - ib + 1;
#line 457 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", &i__4, &i__5, &ib, &c_b25, 
		    &dwork[pxa + nib + nb + 1], n, &a[i__ + (i__ + ib) * 
		    a_dim1], lda, &c_b25, &a[i__ + ib + 1 + (i__ + ib) * 
		    a_dim1], lda, (ftnlen)12, (ftnlen)12);
#line 460 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 460 "MB04TB.f"
	    dgemm_("Transpose", "Transpose", &i__4, n, &ib, &c_b25, &q[i__ + (
		    i__ + ib + 1) * q_dim1], ldq, &dwork[pya], n, &c_b25, &a[
		    i__ + ib + 1 + a_dim1], lda, (ftnlen)9, (ftnlen)9);
#line 463 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 463 "MB04TB.f"
	    dgemm_("Transpose", "Transpose", &i__4, n, &ib, &c_b25, &b[i__ + (
		    i__ + ib + 1) * b_dim1], ldb, &dwork[pya + nib], n, &
		    c_b25, &a[i__ + ib + 1 + a_dim1], lda, (ftnlen)9, (ftnlen)
		    9);

/*           Update the submatrix Q(i+ib:n,i+1+ib:n). */

#line 469 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 469 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 469 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &q[
		    i__ + ib + i__ * q_dim1], ldq, &dwork[pxq + nb + 1], n, &
		    c_b25, &q[i__ + ib + (i__ + ib + 1) * q_dim1], ldq, (
		    ftnlen)12, (ftnlen)9);
#line 472 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 472 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 472 "MB04TB.f"
	    dgemm_("Transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &a[
		    i__ + (i__ + ib) * a_dim1], lda, &dwork[pxq + nib + nb + 
		    1], n, &c_b25, &q[i__ + ib + (i__ + ib + 1) * q_dim1], 
		    ldq, (ftnlen)9, (ftnlen)9);
#line 475 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 475 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 475 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", &i__4, &i__5, &ib, &c_b25, 
		    &dwork[pyq + nb], n, &q[i__ + (i__ + ib + 1) * q_dim1], 
		    ldq, &c_b25, &q[i__ + ib + (i__ + ib + 1) * q_dim1], ldq, 
		    (ftnlen)12, (ftnlen)12);
#line 478 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 478 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 478 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", &i__4, &i__5, &ib, &c_b25, 
		    &dwork[pyq + nib + nb], n, &b[i__ + (i__ + ib + 1) * 
		    b_dim1], ldb, &c_b25, &q[i__ + ib + (i__ + ib + 1) * 
		    q_dim1], ldq, (ftnlen)12, (ftnlen)12);

/*           Update the matrix G. */

#line 484 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 484 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, n, &ib, &c_b25, &q[i__ 
		    + ib + i__ * q_dim1], ldq, &dwork[pxg], n, &c_b25, &g[i__ 
		    + ib + g_dim1], ldg, (ftnlen)12, (ftnlen)9);
#line 487 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 487 "MB04TB.f"
	    dgemm_("Transpose", "Transpose", &i__4, n, &ib, &c_b25, &a[i__ + (
		    i__ + ib) * a_dim1], lda, &dwork[pxg + nib], n, &c_b25, &
		    g[i__ + ib + g_dim1], ldg, (ftnlen)9, (ftnlen)9);
#line 490 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 490 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", n, &i__4, &ib, &c_b25, &
		    dwork[pyg], n, &q[i__ + (i__ + ib + 1) * q_dim1], ldq, &
		    c_b25, &g[(i__ + ib + 1) * g_dim1 + 1], ldg, (ftnlen)12, (
		    ftnlen)12);
#line 493 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 493 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", n, &i__4, &ib, &c_b25, &
		    dwork[pyg + nib], n, &b[i__ + (i__ + ib + 1) * b_dim1], 
		    ldb, &c_b25, &g[(i__ + ib + 1) * g_dim1 + 1], ldg, (
		    ftnlen)12, (ftnlen)12);

/*           Update the submatrix B(i+ib:n,1:n). */

#line 499 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 499 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, n, &ib, &c_b25, &q[i__ 
		    + ib + i__ * q_dim1], ldq, &dwork[pxb], n, &c_b25, &b[i__ 
		    + ib + b_dim1], ldb, (ftnlen)12, (ftnlen)9);
#line 502 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 502 "MB04TB.f"
	    dgemm_("Transpose", "Transpose", &i__4, n, &ib, &c_b25, &a[i__ + (
		    i__ + ib) * a_dim1], lda, &dwork[pxb + nib], n, &c_b25, &
		    b[i__ + ib + b_dim1], ldb, (ftnlen)9, (ftnlen)9);
#line 505 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 505 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 505 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", &i__4, &i__5, &ib, &c_b25, 
		    &dwork[pyb + nb], n, &q[i__ + (i__ + ib + 1) * q_dim1], 
		    ldq, &c_b25, &b[i__ + ib + (i__ + ib + 1) * b_dim1], ldb, 
		    (ftnlen)12, (ftnlen)12);
#line 508 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 508 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 508 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", &i__4, &i__5, &ib, &c_b25, 
		    &dwork[pyb + nib + nb], n, &b[i__ + (i__ + ib + 1) * 
		    b_dim1], ldb, &c_b25, &b[i__ + ib + (i__ + ib + 1) * 
		    b_dim1], ldb, (ftnlen)12, (ftnlen)12);
#line 511 "MB04TB.f"
/* L30: */
#line 511 "MB04TB.f"
	}

#line 513 "MB04TB.f"
    } else if (ltrb) {
#line 514 "MB04TB.f"
	i__1 = *n - nx - 1;
#line 514 "MB04TB.f"
	i__3 = nb;
#line 514 "MB04TB.f"
	for (i__ = *ilo; i__3 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__3) {
/* Computing MIN */
#line 515 "MB04TB.f"
	    i__4 = nb, i__5 = *n - i__;
#line 515 "MB04TB.f"
	    ib = min(i__4,i__5);
#line 516 "MB04TB.f"
	    nib = *n * ib;

/*           Reduce rows and columns i:i+nb-1 to symplectic URV form and */
/*           return the matrices XA, XB, XG, XQ, YA, YB, YG and YQ which */
/*           are needed to update the unreduced parts of the matrices. */

#line 522 "MB04TB.f"
	    i__4 = *n - i__ + 1;
#line 522 "MB04TB.f"
	    i__5 = i__ - 1;
#line 522 "MB04TB.f"
	    mb03xu_(&ltra, &ltrb, &i__4, &i__5, &ib, &a[i__ * a_dim1 + 1], 
		    lda, &b[i__ * b_dim1 + 1], ldb, &g[g_offset], ldg, &q[i__ 
		    + i__ * q_dim1], ldq, &dwork[pxa], n, &dwork[pxb], n, &
		    dwork[pxg], n, &dwork[pxq], n, &dwork[pya], n, &dwork[pyb]
		    , n, &dwork[pyg], n, &dwork[pyq], n, &csl[(i__ << 1) - 1],
		     &csr[(i__ << 1) - 1], &taul[i__], &taur[i__], &dwork[pdw]
		    );

/*           Update the submatrix A(1:n,i+1+ib:n). */

#line 531 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 531 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 531 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &q[
		    i__ + ib + i__ * q_dim1], ldq, &dwork[pxa + nb + 1], n, &
		    c_b25, &a[i__ + ib + (i__ + ib + 1) * a_dim1], lda, (
		    ftnlen)12, (ftnlen)9);
#line 534 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 534 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 534 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &a[
		    i__ + ib + i__ * a_dim1], lda, &dwork[pxa + nib + nb + 1],
		     n, &c_b25, &a[i__ + ib + (i__ + ib + 1) * a_dim1], lda, (
		    ftnlen)12, (ftnlen)9);
#line 537 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 537 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", n, &i__4, &ib, &c_b25, &
		    dwork[pya], n, &q[i__ + (i__ + ib + 1) * q_dim1], ldq, &
		    c_b25, &a[(i__ + ib + 1) * a_dim1 + 1], lda, (ftnlen)12, (
		    ftnlen)12);
#line 540 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 540 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", n, &i__4, &ib, &c_b25, &dwork[
		    pya + nib], n, &b[i__ + ib + 1 + i__ * b_dim1], ldb, &
		    c_b25, &a[(i__ + ib + 1) * a_dim1 + 1], lda, (ftnlen)12, (
		    ftnlen)9);

/*           Update the submatrix Q(i+ib:n,i+1+ib:n). */

#line 546 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 546 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 546 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &q[
		    i__ + ib + i__ * q_dim1], ldq, &dwork[pxq + nb + 1], n, &
		    c_b25, &q[i__ + ib + (i__ + ib + 1) * q_dim1], ldq, (
		    ftnlen)12, (ftnlen)9);
#line 549 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 549 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 549 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &a[
		    i__ + ib + i__ * a_dim1], lda, &dwork[pxq + nib + nb + 1],
		     n, &c_b25, &q[i__ + ib + (i__ + ib + 1) * q_dim1], ldq, (
		    ftnlen)12, (ftnlen)9);
#line 552 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 552 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 552 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", &i__4, &i__5, &ib, &c_b25, 
		    &dwork[pyq + nb], n, &q[i__ + (i__ + ib + 1) * q_dim1], 
		    ldq, &c_b25, &q[i__ + ib + (i__ + ib + 1) * q_dim1], ldq, 
		    (ftnlen)12, (ftnlen)12);
#line 555 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 555 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 555 "MB04TB.f"
	    dgemm_("No Transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &
		    dwork[pyq + nib + nb], n, &b[i__ + ib + 1 + i__ * b_dim1],
		     ldb, &c_b25, &q[i__ + ib + (i__ + ib + 1) * q_dim1], ldq,
		     (ftnlen)12, (ftnlen)9);

/*           Update the matrix G. */

#line 561 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 561 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, n, &ib, &c_b25, &q[i__ 
		    + ib + i__ * q_dim1], ldq, &dwork[pxg], n, &c_b25, &g[i__ 
		    + ib + g_dim1], ldg, (ftnlen)12, (ftnlen)9);
#line 564 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 564 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, n, &ib, &c_b25, &a[i__ 
		    + ib + i__ * a_dim1], lda, &dwork[pxg + nib], n, &c_b25, &
		    g[i__ + ib + g_dim1], ldg, (ftnlen)12, (ftnlen)9);
#line 567 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 567 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", n, &i__4, &ib, &c_b25, &
		    dwork[pyg], n, &q[i__ + (i__ + ib + 1) * q_dim1], ldq, &
		    c_b25, &g[(i__ + ib + 1) * g_dim1 + 1], ldg, (ftnlen)12, (
		    ftnlen)12);
#line 570 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 570 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", n, &i__4, &ib, &c_b25, &dwork[
		    pyg + nib], n, &b[i__ + ib + 1 + i__ * b_dim1], ldb, &
		    c_b25, &g[(i__ + ib + 1) * g_dim1 + 1], ldg, (ftnlen)12, (
		    ftnlen)9);

/*           Update the submatrix B(1:n,i+ib:n). */

#line 576 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 576 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", n, &i__4, &ib, &c_b25, &dwork[
		    pxb], n, &q[i__ + ib + i__ * q_dim1], ldq, &c_b25, &b[(
		    i__ + ib) * b_dim1 + 1], ldb, (ftnlen)12, (ftnlen)9);
#line 579 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 579 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", n, &i__4, &ib, &c_b25, &dwork[
		    pxb + nib], n, &a[i__ + ib + i__ * a_dim1], lda, &c_b25, &
		    b[(i__ + ib) * b_dim1 + 1], ldb, (ftnlen)12, (ftnlen)9);
#line 582 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 582 "MB04TB.f"
	    i__5 = *n - i__ - ib + 1;
#line 582 "MB04TB.f"
	    dgemm_("Transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &q[
		    i__ + (i__ + ib + 1) * q_dim1], ldq, &dwork[pyb + nb], n, 
		    &c_b25, &b[i__ + ib + 1 + (i__ + ib) * b_dim1], ldb, (
		    ftnlen)9, (ftnlen)9);
#line 585 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 585 "MB04TB.f"
	    i__5 = *n - i__ - ib + 1;
#line 585 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &b[
		    i__ + ib + 1 + i__ * b_dim1], ldb, &dwork[pyb + nib + nb],
		     n, &c_b25, &b[i__ + ib + 1 + (i__ + ib) * b_dim1], ldb, (
		    ftnlen)12, (ftnlen)9);
#line 588 "MB04TB.f"
/* L40: */
#line 588 "MB04TB.f"
	}

#line 590 "MB04TB.f"
    } else {
#line 591 "MB04TB.f"
	i__3 = *n - nx - 1;
#line 591 "MB04TB.f"
	i__1 = nb;
#line 591 "MB04TB.f"
	for (i__ = *ilo; i__1 < 0 ? i__ >= i__3 : i__ <= i__3; i__ += i__1) {
/* Computing MIN */
#line 592 "MB04TB.f"
	    i__4 = nb, i__5 = *n - i__;
#line 592 "MB04TB.f"
	    ib = min(i__4,i__5);
#line 593 "MB04TB.f"
	    nib = *n * ib;

/*           Reduce rows and columns i:i+nb-1 to symplectic URV form and */
/*           return the matrices XA, XB, XG, XQ, YA, YB, YG and YQ which */
/*           are needed to update the unreduced parts of the matrices. */

#line 599 "MB04TB.f"
	    i__4 = *n - i__ + 1;
#line 599 "MB04TB.f"
	    i__5 = i__ - 1;
#line 599 "MB04TB.f"
	    mb03xu_(&ltra, &ltrb, &i__4, &i__5, &ib, &a[i__ * a_dim1 + 1], 
		    lda, &b[i__ + b_dim1], ldb, &g[g_offset], ldg, &q[i__ + 
		    i__ * q_dim1], ldq, &dwork[pxa], n, &dwork[pxb], n, &
		    dwork[pxg], n, &dwork[pxq], n, &dwork[pya], n, &dwork[pyb]
		    , n, &dwork[pyg], n, &dwork[pyq], n, &csl[(i__ << 1) - 1],
		     &csr[(i__ << 1) - 1], &taul[i__], &taur[i__], &dwork[pdw]
		    );

/*           Update the submatrix A(1:n,i+1+ib:n). */

#line 608 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 608 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 608 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &q[
		    i__ + ib + i__ * q_dim1], ldq, &dwork[pxa + nb + 1], n, &
		    c_b25, &a[i__ + ib + (i__ + ib + 1) * a_dim1], lda, (
		    ftnlen)12, (ftnlen)9);
#line 611 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 611 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 611 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &a[
		    i__ + ib + i__ * a_dim1], lda, &dwork[pxa + nib + nb + 1],
		     n, &c_b25, &a[i__ + ib + (i__ + ib + 1) * a_dim1], lda, (
		    ftnlen)12, (ftnlen)9);
#line 614 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 614 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", n, &i__4, &ib, &c_b25, &
		    dwork[pya], n, &q[i__ + (i__ + ib + 1) * q_dim1], ldq, &
		    c_b25, &a[(i__ + ib + 1) * a_dim1 + 1], lda, (ftnlen)12, (
		    ftnlen)12);
#line 617 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 617 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", n, &i__4, &ib, &c_b25, &
		    dwork[pya + nib], n, &b[i__ + (i__ + ib + 1) * b_dim1], 
		    ldb, &c_b25, &a[(i__ + ib + 1) * a_dim1 + 1], lda, (
		    ftnlen)12, (ftnlen)12);

/*           Update the submatrix Q(i+ib:n,i+1+ib:n). */

#line 623 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 623 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 623 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &q[
		    i__ + ib + i__ * q_dim1], ldq, &dwork[pxq + nb + 1], n, &
		    c_b25, &q[i__ + ib + (i__ + ib + 1) * q_dim1], ldq, (
		    ftnlen)12, (ftnlen)9);
#line 626 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 626 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 626 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, &i__5, &ib, &c_b25, &a[
		    i__ + ib + i__ * a_dim1], lda, &dwork[pxq + nib + nb + 1],
		     n, &c_b25, &q[i__ + ib + (i__ + ib + 1) * q_dim1], ldq, (
		    ftnlen)12, (ftnlen)9);
#line 629 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 629 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 629 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", &i__4, &i__5, &ib, &c_b25, 
		    &dwork[pyq + nb], n, &q[i__ + (i__ + ib + 1) * q_dim1], 
		    ldq, &c_b25, &q[i__ + ib + (i__ + ib + 1) * q_dim1], ldq, 
		    (ftnlen)12, (ftnlen)12);
#line 632 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 632 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 632 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", &i__4, &i__5, &ib, &c_b25, 
		    &dwork[pyq + nib + nb], n, &b[i__ + (i__ + ib + 1) * 
		    b_dim1], ldb, &c_b25, &q[i__ + ib + (i__ + ib + 1) * 
		    q_dim1], ldq, (ftnlen)12, (ftnlen)12);

/*           Update the matrix G. */

#line 638 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 638 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, n, &ib, &c_b25, &q[i__ 
		    + ib + i__ * q_dim1], ldq, &dwork[pxg], n, &c_b25, &g[i__ 
		    + ib + g_dim1], ldg, (ftnlen)12, (ftnlen)9);
#line 641 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 641 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, n, &ib, &c_b25, &a[i__ 
		    + ib + i__ * a_dim1], lda, &dwork[pxg + nib], n, &c_b25, &
		    g[i__ + ib + g_dim1], ldg, (ftnlen)12, (ftnlen)9);
#line 644 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 644 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", n, &i__4, &ib, &c_b25, &
		    dwork[pyg], n, &q[i__ + (i__ + ib + 1) * q_dim1], ldq, &
		    c_b25, &g[(i__ + ib + 1) * g_dim1 + 1], ldg, (ftnlen)12, (
		    ftnlen)12);
#line 647 "MB04TB.f"
	    i__4 = *n - i__ - ib;
#line 647 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", n, &i__4, &ib, &c_b25, &
		    dwork[pyg + nib], n, &b[i__ + (i__ + ib + 1) * b_dim1], 
		    ldb, &c_b25, &g[(i__ + ib + 1) * g_dim1 + 1], ldg, (
		    ftnlen)12, (ftnlen)12);

/*           Update the submatrix B(i+ib:n,1:n). */

#line 653 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 653 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, n, &ib, &c_b25, &q[i__ 
		    + ib + i__ * q_dim1], ldq, &dwork[pxb], n, &c_b25, &b[i__ 
		    + ib + b_dim1], ldb, (ftnlen)12, (ftnlen)9);
#line 656 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 656 "MB04TB.f"
	    dgemm_("No transpose", "Transpose", &i__4, n, &ib, &c_b25, &a[i__ 
		    + ib + i__ * a_dim1], lda, &dwork[pxb + nib], n, &c_b25, &
		    b[i__ + ib + b_dim1], ldb, (ftnlen)12, (ftnlen)9);
#line 659 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 659 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 659 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", &i__4, &i__5, &ib, &c_b25, 
		    &dwork[pyb + nb], n, &q[i__ + (i__ + ib + 1) * q_dim1], 
		    ldq, &c_b25, &b[i__ + ib + (i__ + ib + 1) * b_dim1], ldb, 
		    (ftnlen)12, (ftnlen)12);
#line 662 "MB04TB.f"
	    i__4 = *n - i__ - ib + 1;
#line 662 "MB04TB.f"
	    i__5 = *n - i__ - ib;
#line 662 "MB04TB.f"
	    dgemm_("No transpose", "No transpose", &i__4, &i__5, &ib, &c_b25, 
		    &dwork[pyb + nib + nb], n, &b[i__ + (i__ + ib + 1) * 
		    b_dim1], ldb, &c_b25, &b[i__ + ib + (i__ + ib + 1) * 
		    b_dim1], ldb, (ftnlen)12, (ftnlen)12);
#line 665 "MB04TB.f"
/* L50: */
#line 665 "MB04TB.f"
	}
#line 666 "MB04TB.f"
    }

/*     Unblocked code to reduce the rest of the matrices. */

#line 670 "MB04TB.f"
    mb04ts_(trana, tranb, n, &i__, &a[a_offset], lda, &b[b_offset], ldb, &g[
	    g_offset], ldg, &q[q_offset], ldq, &csl[1], &csr[1], &taul[1], &
	    taur[1], &dwork[1], ldwork, &ierr, (ftnlen)1, (ftnlen)1);

#line 673 "MB04TB.f"
    dwork[1] = (doublereal) wrkopt;

#line 675 "MB04TB.f"
    return 0;
/* *** Last line of MB04TB *** */
} /* mb04tb_ */

