#line 1 "MB04TS.f"
/* MB04TS.f -- translated by f2c (version 20100827).
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

#line 1 "MB04TS.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mb04ts_(char *trana, char *tranb, integer *n, integer *
	ilo, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *g, integer *ldg, doublereal *q, integer *ldq, doublereal *
	csl, doublereal *csr, doublereal *taul, doublereal *taur, doublereal *
	dwork, integer *ldwork, integer *info, ftnlen trana_len, ftnlen 
	tranb_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, g_dim1, g_offset, q_dim1, 
	    q_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal c__;
    static integer i__;
    static doublereal s, nu;
    static logical ltra, ltrb;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal alpha;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), xerbla_(char *, 
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

/*     To compute a symplectic URV (SURV) decomposition of a real */
/*     2N-by-2N matrix H: */

/*             [ op(A)   G   ]        T       [ op(R11)   R12   ]    T */
/*         H = [             ] = U R V  = U * [                 ] * V , */
/*             [  Q    op(B) ]                [   0     op(R22) ] */

/*     where A, B, G, Q, R12 are real N-by-N matrices, op(R11) is a real */
/*     N-by-N upper triangular matrix, op(R22) is a real N-by-N lower */
/*     Hessenberg matrix and U, V are 2N-by-2N orthogonal symplectic */
/*     matrices. Unblocked version. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op( A ) as follows: */
/*             = 'N': op( A ) = A; */
/*             = 'T': op( A ) = A'; */
/*             = 'C': op( A ) = A'. */

/*     TRANB   CHARACTER*1 */
/*             Specifies the form of op( B ) as follows: */
/*             = 'N': op( B ) = B; */
/*             = 'T': op( B ) = B'; */
/*             = 'C': op( B ) = B'. */

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
/*             The leading dimension of the array Q.  LDG >= MAX(1,N). */

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
/*             applied from the left-hand side. */

/*     TAUR    (output) DOUBLE PRECISION array, dimension (N-1) */
/*             On exit, the first N-1 elements of this array contain the */
/*             scalar factors of some of the elementary reflectors */
/*             applied from the right-hand side. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
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

/*     The algorithm requires 80/3 N**3 + 20 N**2 + O(N) floating point */
/*     operations and is numerically backward stable. */

/*     REFERENCES */

/*     [1] Benner, P., Mehrmann, V., and Xu, H. */
/*         A numerically stable, structure preserving method for */
/*         computing the eigenvalues of real Hamiltonian or symplectic */
/*         pencils. Numer. Math., Vol 78 (3), pp. 329-358, 1998. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DGESUV). */

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

#line 256 "MB04TS.f"
    /* Parameter adjustments */
#line 256 "MB04TS.f"
    a_dim1 = *lda;
#line 256 "MB04TS.f"
    a_offset = 1 + a_dim1;
#line 256 "MB04TS.f"
    a -= a_offset;
#line 256 "MB04TS.f"
    b_dim1 = *ldb;
#line 256 "MB04TS.f"
    b_offset = 1 + b_dim1;
#line 256 "MB04TS.f"
    b -= b_offset;
#line 256 "MB04TS.f"
    g_dim1 = *ldg;
#line 256 "MB04TS.f"
    g_offset = 1 + g_dim1;
#line 256 "MB04TS.f"
    g -= g_offset;
#line 256 "MB04TS.f"
    q_dim1 = *ldq;
#line 256 "MB04TS.f"
    q_offset = 1 + q_dim1;
#line 256 "MB04TS.f"
    q -= q_offset;
#line 256 "MB04TS.f"
    --csl;
#line 256 "MB04TS.f"
    --csr;
#line 256 "MB04TS.f"
    --taul;
#line 256 "MB04TS.f"
    --taur;
#line 256 "MB04TS.f"
    --dwork;
#line 256 "MB04TS.f"

#line 256 "MB04TS.f"
    /* Function Body */
#line 256 "MB04TS.f"
    *info = 0;
#line 257 "MB04TS.f"
    ltra = lsame_(trana, "T", (ftnlen)1, (ftnlen)1) || lsame_(trana, "C", (
	    ftnlen)1, (ftnlen)1);
#line 258 "MB04TS.f"
    ltrb = lsame_(tranb, "T", (ftnlen)1, (ftnlen)1) || lsame_(tranb, "C", (
	    ftnlen)1, (ftnlen)1);
#line 259 "MB04TS.f"
    if (! ltra && ! lsame_(trana, "N", (ftnlen)1, (ftnlen)1)) {
#line 260 "MB04TS.f"
	*info = -1;
#line 261 "MB04TS.f"
    } else if (! ltrb && ! lsame_(tranb, "N", (ftnlen)1, (ftnlen)1)) {
#line 262 "MB04TS.f"
	*info = -2;
#line 263 "MB04TS.f"
    } else if (*n < 0) {
#line 264 "MB04TS.f"
	*info = -3;
#line 265 "MB04TS.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 266 "MB04TS.f"
	*info = -4;
#line 267 "MB04TS.f"
    } else if (*lda < max(1,*n)) {
#line 268 "MB04TS.f"
	*info = -6;
#line 269 "MB04TS.f"
    } else if (*ldb < max(1,*n)) {
#line 270 "MB04TS.f"
	*info = -8;
#line 271 "MB04TS.f"
    } else if (*ldg < max(1,*n)) {
#line 272 "MB04TS.f"
	*info = -10;
#line 273 "MB04TS.f"
    } else if (*ldq < max(1,*n)) {
#line 274 "MB04TS.f"
	*info = -12;
#line 275 "MB04TS.f"
    } else if (*ldwork < max(1,*n)) {
#line 276 "MB04TS.f"
	dwork[1] = (doublereal) max(1,*n);
#line 277 "MB04TS.f"
	*info = -18;
#line 278 "MB04TS.f"
    }

/*     Return if there were illegal values. */

#line 282 "MB04TS.f"
    if (*info != 0) {
#line 283 "MB04TS.f"
	i__1 = -(*info);
#line 283 "MB04TS.f"
	xerbla_("MB04TS", &i__1, (ftnlen)6);
#line 284 "MB04TS.f"
	return 0;
#line 285 "MB04TS.f"
    }

/*     Quick return if possible. */

#line 289 "MB04TS.f"
    if (*n == 0) {
#line 290 "MB04TS.f"
	dwork[1] = 1.;
#line 291 "MB04TS.f"
	return 0;
#line 292 "MB04TS.f"
    }

#line 294 "MB04TS.f"
    i__1 = *n;
#line 294 "MB04TS.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 295 "MB04TS.f"
	alpha = q[i__ + i__ * q_dim1];
#line 296 "MB04TS.f"
	if (i__ < *n) {

/*           Generate elementary reflector HU(i) to annihilate Q(i+1:n,i) */

#line 300 "MB04TS.f"
	    i__2 = *n - i__ + 1;
#line 300 "MB04TS.f"
	    dlarfg_(&i__2, &alpha, &q[i__ + 1 + i__ * q_dim1], &c__1, &nu);

/*           Apply HU(i) from the left. */

#line 304 "MB04TS.f"
	    q[i__ + i__ * q_dim1] = 1.;
#line 305 "MB04TS.f"
	    i__2 = *n - i__ + 1;
#line 305 "MB04TS.f"
	    i__3 = *n - i__;
#line 305 "MB04TS.f"
	    dlarf_("Left", &i__2, &i__3, &q[i__ + i__ * q_dim1], &c__1, &nu, &
		    q[i__ + (i__ + 1) * q_dim1], ldq, &dwork[1], (ftnlen)4);
#line 307 "MB04TS.f"
	    if (ltra) {
#line 308 "MB04TS.f"
		i__2 = *n - i__ + 1;
#line 308 "MB04TS.f"
		i__3 = *n - i__ + 1;
#line 308 "MB04TS.f"
		dlarf_("Right", &i__2, &i__3, &q[i__ + i__ * q_dim1], &c__1, &
			nu, &a[i__ + i__ * a_dim1], lda, &dwork[1], (ftnlen)5)
			;
#line 310 "MB04TS.f"
	    } else {
#line 311 "MB04TS.f"
		i__2 = *n - i__ + 1;
#line 311 "MB04TS.f"
		i__3 = *n - i__ + 1;
#line 311 "MB04TS.f"
		dlarf_("Left", &i__2, &i__3, &q[i__ + i__ * q_dim1], &c__1, &
			nu, &a[i__ + i__ * a_dim1], lda, &dwork[1], (ftnlen)4)
			;
#line 313 "MB04TS.f"
	    }
#line 314 "MB04TS.f"
	    if (ltrb) {
#line 315 "MB04TS.f"
		i__2 = *n - i__ + 1;
#line 315 "MB04TS.f"
		dlarf_("Right", n, &i__2, &q[i__ + i__ * q_dim1], &c__1, &nu, 
			&b[i__ * b_dim1 + 1], ldb, &dwork[1], (ftnlen)5);
#line 317 "MB04TS.f"
	    } else {
#line 318 "MB04TS.f"
		i__2 = *n - i__ + 1;
#line 318 "MB04TS.f"
		dlarf_("Left", &i__2, n, &q[i__ + i__ * q_dim1], &c__1, &nu, &
			b[i__ + b_dim1], ldb, &dwork[1], (ftnlen)4);
#line 320 "MB04TS.f"
	    }
#line 321 "MB04TS.f"
	    i__2 = *n - i__ + 1;
#line 321 "MB04TS.f"
	    dlarf_("Left", &i__2, n, &q[i__ + i__ * q_dim1], &c__1, &nu, &g[
		    i__ + g_dim1], ldg, &dwork[1], (ftnlen)4);
#line 323 "MB04TS.f"
	    q[i__ + i__ * q_dim1] = nu;
#line 324 "MB04TS.f"
	} else {
#line 325 "MB04TS.f"
	    q[i__ + i__ * q_dim1] = 0.;
#line 326 "MB04TS.f"
	}

/*        Generate symplectic Givens rotator GU(i) to annihilate Q(i,i). */

#line 330 "MB04TS.f"
	temp = a[i__ + i__ * a_dim1];
#line 331 "MB04TS.f"
	dlartg_(&temp, &alpha, &c__, &s, &a[i__ + i__ * a_dim1]);

/*        Apply G(i) from the left. */

#line 335 "MB04TS.f"
	if (ltra) {
#line 336 "MB04TS.f"
	    i__2 = *n - i__;
#line 336 "MB04TS.f"
	    drot_(&i__2, &a[i__ + 1 + i__ * a_dim1], &c__1, &q[i__ + (i__ + 1)
		     * q_dim1], ldq, &c__, &s);
#line 337 "MB04TS.f"
	} else {
#line 338 "MB04TS.f"
	    i__2 = *n - i__;
#line 338 "MB04TS.f"
	    drot_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda, &q[i__ + (i__ + 1)
		     * q_dim1], ldq, &c__, &s);
#line 339 "MB04TS.f"
	}
#line 340 "MB04TS.f"
	if (ltrb) {
#line 341 "MB04TS.f"
	    drot_(n, &g[i__ + g_dim1], ldg, &b[i__ * b_dim1 + 1], &c__1, &c__,
		     &s);
#line 342 "MB04TS.f"
	} else {
#line 343 "MB04TS.f"
	    drot_(n, &g[i__ + g_dim1], ldg, &b[i__ + b_dim1], ldb, &c__, &s);
#line 344 "MB04TS.f"
	}
#line 345 "MB04TS.f"
	csl[(i__ << 1) - 1] = c__;
#line 346 "MB04TS.f"
	csl[i__ * 2] = s;

#line 348 "MB04TS.f"
	if (i__ < *n) {
#line 349 "MB04TS.f"
	    if (ltra) {

/*              Generate elementary reflector FU(i) to annihilate */
/*              A(i,i+1:n). */

#line 354 "MB04TS.f"
		i__2 = *n - i__ + 1;
#line 354 "MB04TS.f"
		dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + (i__ + 1) * 
			a_dim1], lda, &taul[i__]);

/*              Apply FU(i) from the left. */

#line 358 "MB04TS.f"
		temp = a[i__ + i__ * a_dim1];
#line 359 "MB04TS.f"
		a[i__ + i__ * a_dim1] = 1.;
#line 360 "MB04TS.f"
		i__2 = *n - i__;
#line 360 "MB04TS.f"
		i__3 = *n - i__ + 1;
#line 360 "MB04TS.f"
		dlarf_("Right", &i__2, &i__3, &a[i__ + i__ * a_dim1], lda, &
			taul[i__], &a[i__ + 1 + i__ * a_dim1], lda, &dwork[1],
			 (ftnlen)5);
#line 362 "MB04TS.f"
		i__2 = *n - i__ + 1;
#line 362 "MB04TS.f"
		i__3 = *n - i__;
#line 362 "MB04TS.f"
		dlarf_("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], lda, &
			taul[i__], &q[i__ + (i__ + 1) * q_dim1], ldq, &dwork[
			1], (ftnlen)4);
#line 364 "MB04TS.f"
		if (ltrb) {
#line 365 "MB04TS.f"
		    i__2 = *n - i__ + 1;
#line 365 "MB04TS.f"
		    dlarf_("Right", n, &i__2, &a[i__ + i__ * a_dim1], lda, &
			    taul[i__], &b[i__ * b_dim1 + 1], ldb, &dwork[1], (
			    ftnlen)5);
#line 367 "MB04TS.f"
		} else {
#line 368 "MB04TS.f"
		    i__2 = *n - i__ + 1;
#line 368 "MB04TS.f"
		    dlarf_("Left", &i__2, n, &a[i__ + i__ * a_dim1], lda, &
			    taul[i__], &b[i__ + b_dim1], ldb, &dwork[1], (
			    ftnlen)4);
#line 370 "MB04TS.f"
		}
#line 371 "MB04TS.f"
		i__2 = *n - i__ + 1;
#line 371 "MB04TS.f"
		dlarf_("Left", &i__2, n, &a[i__ + i__ * a_dim1], lda, &taul[
			i__], &g[i__ + g_dim1], ldg, &dwork[1], (ftnlen)4);
#line 373 "MB04TS.f"
		a[i__ + i__ * a_dim1] = temp;
#line 374 "MB04TS.f"
	    } else {

/*              Generate elementary reflector FU(i) to annihilate */
/*              A(i+1:n,i). */

#line 379 "MB04TS.f"
		i__2 = *n - i__ + 1;
#line 379 "MB04TS.f"
		dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + 1 + i__ * 
			a_dim1], &c__1, &taul[i__]);

/*              Apply FU(i) from the left. */

#line 383 "MB04TS.f"
		temp = a[i__ + i__ * a_dim1];
#line 384 "MB04TS.f"
		a[i__ + i__ * a_dim1] = 1.;
#line 385 "MB04TS.f"
		i__2 = *n - i__ + 1;
#line 385 "MB04TS.f"
		i__3 = *n - i__;
#line 385 "MB04TS.f"
		dlarf_("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &
			taul[i__], &a[i__ + (i__ + 1) * a_dim1], lda, &dwork[
			1], (ftnlen)4);
#line 387 "MB04TS.f"
		i__2 = *n - i__ + 1;
#line 387 "MB04TS.f"
		i__3 = *n - i__;
#line 387 "MB04TS.f"
		dlarf_("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &
			taul[i__], &q[i__ + (i__ + 1) * q_dim1], ldq, &dwork[
			1], (ftnlen)4);
#line 389 "MB04TS.f"
		if (ltrb) {
#line 390 "MB04TS.f"
		    i__2 = *n - i__ + 1;
#line 390 "MB04TS.f"
		    dlarf_("Right", n, &i__2, &a[i__ + i__ * a_dim1], &c__1, &
			    taul[i__], &b[i__ * b_dim1 + 1], ldb, &dwork[1], (
			    ftnlen)5);
#line 392 "MB04TS.f"
		} else {
#line 393 "MB04TS.f"
		    i__2 = *n - i__ + 1;
#line 393 "MB04TS.f"
		    dlarf_("Left", &i__2, n, &a[i__ + i__ * a_dim1], &c__1, &
			    taul[i__], &b[i__ + b_dim1], ldb, &dwork[1], (
			    ftnlen)4);
#line 395 "MB04TS.f"
		}
#line 396 "MB04TS.f"
		i__2 = *n - i__ + 1;
#line 396 "MB04TS.f"
		dlarf_("Left", &i__2, n, &a[i__ + i__ * a_dim1], &c__1, &taul[
			i__], &g[i__ + g_dim1], ldg, &dwork[1], (ftnlen)4);
#line 398 "MB04TS.f"
		a[i__ + i__ * a_dim1] = temp;
#line 399 "MB04TS.f"
	    }
#line 400 "MB04TS.f"
	} else {
#line 401 "MB04TS.f"
	    taul[i__] = 0.;
#line 402 "MB04TS.f"
	}
#line 403 "MB04TS.f"
	if (i__ < *n) {
#line 403 "MB04TS.f"
	    alpha = q[i__ + (i__ + 1) * q_dim1];
#line 403 "MB04TS.f"
	}
#line 405 "MB04TS.f"
	if (i__ < *n - 1) {

/*           Generate elementary reflector HV(i) to annihilate Q(i,i+2:n) */

#line 409 "MB04TS.f"
	    i__2 = *n - i__;
#line 409 "MB04TS.f"
	    dlarfg_(&i__2, &alpha, &q[i__ + (i__ + 2) * q_dim1], ldq, &nu);

/*           Apply HV(i) from the right. */

#line 413 "MB04TS.f"
	    q[i__ + (i__ + 1) * q_dim1] = 1.;
#line 414 "MB04TS.f"
	    i__2 = *n - i__;
#line 414 "MB04TS.f"
	    i__3 = *n - i__;
#line 414 "MB04TS.f"
	    dlarf_("Right", &i__2, &i__3, &q[i__ + (i__ + 1) * q_dim1], ldq, &
		    nu, &q[i__ + 1 + (i__ + 1) * q_dim1], ldq, &dwork[1], (
		    ftnlen)5);
#line 416 "MB04TS.f"
	    if (ltra) {
#line 417 "MB04TS.f"
		i__2 = *n - i__;
#line 417 "MB04TS.f"
		dlarf_("Left", &i__2, n, &q[i__ + (i__ + 1) * q_dim1], ldq, &
			nu, &a[i__ + 1 + a_dim1], lda, &dwork[1], (ftnlen)4);
#line 419 "MB04TS.f"
	    } else {
#line 420 "MB04TS.f"
		i__2 = *n - i__;
#line 420 "MB04TS.f"
		dlarf_("Right", n, &i__2, &q[i__ + (i__ + 1) * q_dim1], ldq, &
			nu, &a[(i__ + 1) * a_dim1 + 1], lda, &dwork[1], (
			ftnlen)5);
#line 422 "MB04TS.f"
	    }
#line 423 "MB04TS.f"
	    if (ltrb) {
#line 424 "MB04TS.f"
		i__2 = *n - i__;
#line 424 "MB04TS.f"
		i__3 = *n - i__ + 1;
#line 424 "MB04TS.f"
		dlarf_("Left", &i__2, &i__3, &q[i__ + (i__ + 1) * q_dim1], 
			ldq, &nu, &b[i__ + 1 + i__ * b_dim1], ldb, &dwork[1], 
			(ftnlen)4);
#line 426 "MB04TS.f"
	    } else {
#line 427 "MB04TS.f"
		i__2 = *n - i__ + 1;
#line 427 "MB04TS.f"
		i__3 = *n - i__;
#line 427 "MB04TS.f"
		dlarf_("Right", &i__2, &i__3, &q[i__ + (i__ + 1) * q_dim1], 
			ldq, &nu, &b[i__ + (i__ + 1) * b_dim1], ldb, &dwork[1]
			, (ftnlen)5);
#line 429 "MB04TS.f"
	    }
#line 430 "MB04TS.f"
	    i__2 = *n - i__;
#line 430 "MB04TS.f"
	    dlarf_("Right", n, &i__2, &q[i__ + (i__ + 1) * q_dim1], ldq, &nu, 
		    &g[(i__ + 1) * g_dim1 + 1], ldg, &dwork[1], (ftnlen)5);
#line 432 "MB04TS.f"
	    q[i__ + (i__ + 1) * q_dim1] = nu;
#line 433 "MB04TS.f"
	} else if (i__ < *n) {
#line 434 "MB04TS.f"
	    q[i__ + (i__ + 1) * q_dim1] = 0.;
#line 435 "MB04TS.f"
	}
#line 436 "MB04TS.f"
	if (i__ < *n) {

/*           Generate symplectic Givens rotator GV(i) to annihilate */
/*           Q(i,i+1). */

#line 441 "MB04TS.f"
	    if (ltrb) {
#line 442 "MB04TS.f"
		temp = b[i__ + 1 + i__ * b_dim1];
#line 443 "MB04TS.f"
		dlartg_(&temp, &alpha, &c__, &s, &b[i__ + 1 + i__ * b_dim1]);
#line 444 "MB04TS.f"
		s = -s;
#line 445 "MB04TS.f"
		i__2 = *n - i__;
#line 445 "MB04TS.f"
		drot_(&i__2, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, &b[i__ 
			+ 1 + (i__ + 1) * b_dim1], ldb, &c__, &s);
#line 446 "MB04TS.f"
	    } else {
#line 447 "MB04TS.f"
		temp = b[i__ + (i__ + 1) * b_dim1];
#line 448 "MB04TS.f"
		dlartg_(&temp, &alpha, &c__, &s, &b[i__ + (i__ + 1) * b_dim1])
			;
#line 449 "MB04TS.f"
		s = -s;
#line 450 "MB04TS.f"
		i__2 = *n - i__;
#line 450 "MB04TS.f"
		drot_(&i__2, &q[i__ + 1 + (i__ + 1) * q_dim1], &c__1, &b[i__ 
			+ 1 + (i__ + 1) * b_dim1], &c__1, &c__, &s);
#line 451 "MB04TS.f"
	    }
#line 452 "MB04TS.f"
	    if (ltra) {
#line 453 "MB04TS.f"
		drot_(n, &a[i__ + 1 + a_dim1], lda, &g[(i__ + 1) * g_dim1 + 1]
			, &c__1, &c__, &s);
#line 454 "MB04TS.f"
	    } else {
#line 455 "MB04TS.f"
		drot_(n, &a[(i__ + 1) * a_dim1 + 1], &c__1, &g[(i__ + 1) * 
			g_dim1 + 1], &c__1, &c__, &s);
#line 456 "MB04TS.f"
	    }
#line 457 "MB04TS.f"
	    csr[(i__ << 1) - 1] = c__;
#line 458 "MB04TS.f"
	    csr[i__ * 2] = s;
#line 459 "MB04TS.f"
	}
#line 460 "MB04TS.f"
	if (i__ < *n - 1) {
#line 461 "MB04TS.f"
	    if (ltrb) {

/*              Generate elementary reflector FV(i) to annihilate */
/*              B(i+2:n,i). */

#line 466 "MB04TS.f"
		i__2 = *n - i__;
#line 466 "MB04TS.f"
		dlarfg_(&i__2, &b[i__ + 1 + i__ * b_dim1], &b[i__ + 2 + i__ * 
			b_dim1], &c__1, &taur[i__]);

/*              Apply FV(i) from the right. */

#line 470 "MB04TS.f"
		temp = b[i__ + 1 + i__ * b_dim1];
#line 471 "MB04TS.f"
		b[i__ + 1 + i__ * b_dim1] = 1.;
#line 472 "MB04TS.f"
		i__2 = *n - i__;
#line 472 "MB04TS.f"
		i__3 = *n - i__;
#line 472 "MB04TS.f"
		dlarf_("Left", &i__2, &i__3, &b[i__ + 1 + i__ * b_dim1], &
			c__1, &taur[i__], &b[i__ + 1 + (i__ + 1) * b_dim1], 
			ldb, &dwork[1], (ftnlen)4);
#line 474 "MB04TS.f"
		i__2 = *n - i__;
#line 474 "MB04TS.f"
		i__3 = *n - i__;
#line 474 "MB04TS.f"
		dlarf_("Right", &i__2, &i__3, &b[i__ + 1 + i__ * b_dim1], &
			c__1, &taur[i__], &q[i__ + 1 + (i__ + 1) * q_dim1], 
			ldq, &dwork[1], (ftnlen)5);
#line 476 "MB04TS.f"
		if (ltra) {
#line 477 "MB04TS.f"
		    i__2 = *n - i__;
#line 477 "MB04TS.f"
		    dlarf_("Left", &i__2, n, &b[i__ + 1 + i__ * b_dim1], &
			    c__1, &taur[i__], &a[i__ + 1 + a_dim1], lda, &
			    dwork[1], (ftnlen)4);
#line 479 "MB04TS.f"
		} else {
#line 480 "MB04TS.f"
		    i__2 = *n - i__;
#line 480 "MB04TS.f"
		    dlarf_("Right", n, &i__2, &b[i__ + 1 + i__ * b_dim1], &
			    c__1, &taur[i__], &a[(i__ + 1) * a_dim1 + 1], lda,
			     &dwork[1], (ftnlen)5);
#line 482 "MB04TS.f"
		}
#line 483 "MB04TS.f"
		i__2 = *n - i__;
#line 483 "MB04TS.f"
		dlarf_("Right", n, &i__2, &b[i__ + 1 + i__ * b_dim1], &c__1, &
			taur[i__], &g[(i__ + 1) * g_dim1 + 1], ldg, &dwork[1],
			 (ftnlen)5);
#line 485 "MB04TS.f"
		b[i__ + 1 + i__ * b_dim1] = temp;
#line 486 "MB04TS.f"
	    } else {

/*              Generate elementary reflector FV(i) to annihilate */
/*              B(i,i+2:n). */

#line 491 "MB04TS.f"
		i__2 = *n - i__;
#line 491 "MB04TS.f"
		dlarfg_(&i__2, &b[i__ + (i__ + 1) * b_dim1], &b[i__ + (i__ + 
			2) * b_dim1], ldb, &taur[i__]);

/*              Apply FV(i) from the right. */

#line 495 "MB04TS.f"
		temp = b[i__ + (i__ + 1) * b_dim1];
#line 496 "MB04TS.f"
		b[i__ + (i__ + 1) * b_dim1] = 1.;
#line 497 "MB04TS.f"
		i__2 = *n - i__;
#line 497 "MB04TS.f"
		i__3 = *n - i__;
#line 497 "MB04TS.f"
		dlarf_("Right", &i__2, &i__3, &b[i__ + (i__ + 1) * b_dim1], 
			ldb, &taur[i__], &b[i__ + 1 + (i__ + 1) * b_dim1], 
			ldb, &dwork[1], (ftnlen)5);
#line 499 "MB04TS.f"
		i__2 = *n - i__;
#line 499 "MB04TS.f"
		i__3 = *n - i__;
#line 499 "MB04TS.f"
		dlarf_("Right", &i__2, &i__3, &b[i__ + (i__ + 1) * b_dim1], 
			ldb, &taur[i__], &q[i__ + 1 + (i__ + 1) * q_dim1], 
			ldq, &dwork[1], (ftnlen)5);
#line 501 "MB04TS.f"
		if (ltra) {
#line 502 "MB04TS.f"
		    i__2 = *n - i__;
#line 502 "MB04TS.f"
		    dlarf_("Left", &i__2, n, &b[i__ + (i__ + 1) * b_dim1], 
			    ldb, &taur[i__], &a[i__ + 1 + a_dim1], lda, &
			    dwork[1], (ftnlen)4);
#line 504 "MB04TS.f"
		} else {
#line 505 "MB04TS.f"
		    i__2 = *n - i__;
#line 505 "MB04TS.f"
		    dlarf_("Right", n, &i__2, &b[i__ + (i__ + 1) * b_dim1], 
			    ldb, &taur[i__], &a[(i__ + 1) * a_dim1 + 1], lda, 
			    &dwork[1], (ftnlen)5);
#line 507 "MB04TS.f"
		}
#line 508 "MB04TS.f"
		i__2 = *n - i__;
#line 508 "MB04TS.f"
		dlarf_("Right", n, &i__2, &b[i__ + (i__ + 1) * b_dim1], ldb, &
			taur[i__], &g[(i__ + 1) * g_dim1 + 1], ldg, &dwork[1],
			 (ftnlen)5);
#line 510 "MB04TS.f"
		b[i__ + (i__ + 1) * b_dim1] = temp;
#line 511 "MB04TS.f"
	    }
#line 512 "MB04TS.f"
	} else if (i__ < *n) {
#line 513 "MB04TS.f"
	    taur[i__] = 0.;
#line 514 "MB04TS.f"
	}
#line 515 "MB04TS.f"
/* L10: */
#line 515 "MB04TS.f"
    }
#line 516 "MB04TS.f"
    dwork[1] = (doublereal) max(1,*n);
#line 517 "MB04TS.f"
    return 0;
/* *** Last line of MB04TS *** */
} /* mb04ts_ */

