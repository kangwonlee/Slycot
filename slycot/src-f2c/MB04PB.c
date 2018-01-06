#line 1 "MB04PB.f"
/* MB04PB.f -- translated by f2c (version 20100827).
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

#line 1 "MB04PB.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
static logical c_true = TRUE_;
static doublereal c_b20 = 1.;

/* Subroutine */ int mb04pb_(integer *n, integer *ilo, doublereal *a, integer 
	*lda, doublereal *qg, integer *ldqg, doublereal *cs, doublereal *tau, 
	doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, qg_dim1, qg_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, ib, nb, nh, nx, nib, nnb, pxa, pdw, pya, pxg, pxq, 
	    ierr;
    extern /* Subroutine */ int mb04pa_(logical *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *), dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern integer ue01md_(integer *, char *, char *, integer *, integer *, 
	    integer *, ftnlen, ftnlen);
    static integer nbmin;
    extern /* Subroutine */ int mb04pu_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *), dsyr2k_(char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
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

/*     To reduce a Hamiltonian matrix, */

/*                   [  A   G  ] */
/*              H =  [       T ] , */
/*                   [  Q  -A  ] */

/*     where A is an N-by-N matrix and G,Q are N-by-N symmetric matrices, */
/*     to Paige/Van Loan (PVL) form. That is, an orthogonal symplectic U */
/*     is computed so that */

/*               T       [  Aout   Gout  ] */
/*              U H U =  [             T ] , */
/*                       [  Qout  -Aout  ] */

/*     where Aout is upper Hessenberg and Qout is diagonal. */
/*     Blocked version. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     ILO     (input) INTEGER */
/*             It is assumed that A is already upper triangular and Q is */
/*             zero in rows and columns 1:ILO-1. ILO is normally set by a */
/*             previous call to MB04DD; otherwise it should be set to 1. */
/*             1 <= ILO <= N, if N > 0; ILO = 1, if N = 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the matrix Aout and, in the zero part of Aout, */
/*             information about the elementary reflectors used to */
/*             compute the PVL factorization. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     QG      (input/output) DOUBLE PRECISION array, dimension */
/*                            (LDQG,N+1) */
/*             On entry, the leading N-by-N+1 part of this array must */
/*             contain the lower triangular part of the matrix Q and */
/*             the upper triangular part of the matrix G. */
/*             On exit, the leading N-by-N+1 part of this array contains */
/*             the diagonal of the matrix Qout, the upper triangular part */
/*             of the matrix Gout and, in the zero parts of Qout, */
/*             information about the elementary reflectors used to */
/*             compute the PVL factorization. */

/*     LDQG    INTEGER */
/*             The leading dimension of the array QG.  LDQG >= MAX(1,N). */

/*     CS      (output) DOUBLE PRECISION array, dimension (2N-2) */
/*             On exit, the first 2N-2 elements of this array contain the */
/*             cosines and sines of the symplectic Givens rotations used */
/*             to compute the PVL factorization. */

/*     TAU     (output) DOUBLE PRECISION array, dimension (N-1) */
/*             On exit, the first N-1 elements of this array contain the */
/*             scalar factors of some of the elementary reflectors. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK, 8*N*NB + 3*NB, where NB is the optimal */
/*             block size determined by the function UE01MD. */
/*             On exit, if  INFO = -10,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,N-1). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The matrix U is represented as a product of symplectic reflectors */
/*     and Givens rotators */

/*     U = diag( H(1),H(1) )     G(1)   diag( F(1),F(1) ) */
/*         diag( H(2),H(2) )     G(2)   diag( F(2),F(2) ) */
/*                                .... */
/*         diag( H(n-1),H(n-1) ) G(n-1) diag( F(n-1),F(n-1) ). */

/*     Each H(i) has the form */

/*           H(i) = I - tau * v * v' */

/*     where tau is a real scalar, and v is a real vector with */
/*     v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in */
/*     QG(i+2:n,i), and tau in QG(i+1,i). */

/*     Each F(i) has the form */

/*           F(i) = I - nu * w * w' */

/*     where nu is a real scalar, and w is a real vector with */
/*     w(1:i) = 0 and w(i+1) = 1; w(i+2:n) is stored on exit in */
/*     A(i+2:n,i), and nu in TAU(i). */

/*     Each G(i) is a Givens rotator acting on rows i+1 and n+i+1, */
/*     where the cosine is stored in CS(2*i-1) and the sine in */
/*     CS(2*i). */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires O(N**3) floating point operations and is */
/*     strongly backward stable. */

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

/*     V. Sima, Nov. 2008 (SLICOT version of the HAPACK routine DHAPVB). */

/*     KEYWORDS */

/*     Elementary matrix operations, Hamiltonian matrix. */

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

#line 190 "MB04PB.f"
    /* Parameter adjustments */
#line 190 "MB04PB.f"
    a_dim1 = *lda;
#line 190 "MB04PB.f"
    a_offset = 1 + a_dim1;
#line 190 "MB04PB.f"
    a -= a_offset;
#line 190 "MB04PB.f"
    qg_dim1 = *ldqg;
#line 190 "MB04PB.f"
    qg_offset = 1 + qg_dim1;
#line 190 "MB04PB.f"
    qg -= qg_offset;
#line 190 "MB04PB.f"
    --cs;
#line 190 "MB04PB.f"
    --tau;
#line 190 "MB04PB.f"
    --dwork;
#line 190 "MB04PB.f"

#line 190 "MB04PB.f"
    /* Function Body */
#line 190 "MB04PB.f"
    *info = 0;
#line 191 "MB04PB.f"
    if (*n < 0) {
#line 192 "MB04PB.f"
	*info = -1;
#line 193 "MB04PB.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 194 "MB04PB.f"
	*info = -2;
#line 195 "MB04PB.f"
    } else if (*lda < max(1,*n)) {
#line 196 "MB04PB.f"
	*info = -4;
#line 197 "MB04PB.f"
    } else if (*ldqg < max(1,*n)) {
#line 198 "MB04PB.f"
	*info = -6;
#line 199 "MB04PB.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 199 "MB04PB.f"
	i__1 = 1, i__2 = *n - 1;
#line 199 "MB04PB.f"
	if (*ldwork < max(i__1,i__2)) {
/* Computing MAX */
#line 200 "MB04PB.f"
	    i__1 = 1, i__2 = *n - 1;
#line 200 "MB04PB.f"
	    dwork[1] = (doublereal) max(i__1,i__2);
#line 201 "MB04PB.f"
	    *info = -10;
#line 202 "MB04PB.f"
	}
#line 202 "MB04PB.f"
    }

/*     Return if there were illegal values. */

#line 206 "MB04PB.f"
    if (*info != 0) {
#line 207 "MB04PB.f"
	i__1 = -(*info);
#line 207 "MB04PB.f"
	xerbla_("MB04PB", &i__1, (ftnlen)6);
#line 208 "MB04PB.f"
	return 0;
#line 209 "MB04PB.f"
    }

/*     Set elements 1:ILO-1 of TAU and CS. */

#line 213 "MB04PB.f"
    i__1 = *ilo - 1;
#line 213 "MB04PB.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 214 "MB04PB.f"
	tau[i__] = 0.;
#line 215 "MB04PB.f"
	cs[(i__ << 1) - 1] = 1.;
#line 216 "MB04PB.f"
	cs[i__ * 2] = 0.;
#line 217 "MB04PB.f"
/* L10: */
#line 217 "MB04PB.f"
    }

/*     Quick return if possible. */

#line 221 "MB04PB.f"
    if (*n <= *ilo) {
#line 222 "MB04PB.f"
	dwork[1] = 1.;
#line 223 "MB04PB.f"
	return 0;
#line 224 "MB04PB.f"
    }

/*     Determine the block size. */

#line 228 "MB04PB.f"
    nh = *n - *ilo + 1;
#line 229 "MB04PB.f"
    nb = ue01md_(&c__1, "MB04PB", " ", n, ilo, &c_n1, (ftnlen)6, (ftnlen)1);
#line 230 "MB04PB.f"
    nbmin = 2;
#line 231 "MB04PB.f"
    wrkopt = *n - 1;
#line 232 "MB04PB.f"
    if (nb > 1 && nb < nh) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
#line 236 "MB04PB.f"
	i__1 = nb, i__2 = ue01md_(&c__3, "MB04PB", " ", n, ilo, &c_n1, (
		ftnlen)6, (ftnlen)1);
#line 236 "MB04PB.f"
	nx = max(i__1,i__2);
#line 237 "MB04PB.f"
	if (nx < nh) {

/*           Check whether workspace is large enough for blocked code. */

#line 241 "MB04PB.f"
	    wrkopt = (*n << 3) * nb + nb * 3;
#line 242 "MB04PB.f"
	    if (*ldwork < wrkopt) {

/*              Not enough workspace available. Determine minimum value */
/*              of NB, and reduce NB. */

/* Computing MAX */
#line 247 "MB04PB.f"
		i__1 = 2, i__2 = ue01md_(&c__2, "MB04PB", " ", n, ilo, &c_n1, 
			(ftnlen)6, (ftnlen)1);
#line 247 "MB04PB.f"
		nbmin = max(i__1,i__2);
#line 248 "MB04PB.f"
		nb = *ldwork / ((*n << 3) + 3);
#line 249 "MB04PB.f"
	    }
#line 250 "MB04PB.f"
	}
#line 251 "MB04PB.f"
    }

#line 253 "MB04PB.f"
    nnb = *n * nb;
#line 254 "MB04PB.f"
    pxa = 1;
#line 255 "MB04PB.f"
    pya = pxa + (nnb << 1);
#line 256 "MB04PB.f"
    pxq = pya + (nnb << 1);
#line 257 "MB04PB.f"
    pxg = pxq + (nnb << 1);
#line 258 "MB04PB.f"
    pdw = pxg + (nnb << 1);

#line 260 "MB04PB.f"
    if (nb < nbmin || nb >= nh) {

/*        Use unblocked code. */

#line 264 "MB04PB.f"
	i__ = *ilo;

#line 266 "MB04PB.f"
    } else {
#line 267 "MB04PB.f"
	i__1 = *n - nx - 1;
#line 267 "MB04PB.f"
	i__2 = nb;
#line 267 "MB04PB.f"
	for (i__ = *ilo; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 268 "MB04PB.f"
	    i__3 = nb, i__4 = *n - i__;
#line 268 "MB04PB.f"
	    ib = min(i__3,i__4);
#line 269 "MB04PB.f"
	    nib = *n * ib;

/*           Reduce rows and columns i:i+nb-1 to PVL form and return the */
/*           matrices XA, XG, XQ, and YA which are needed to update the */
/*           unreduced parts of the matrices. */

#line 275 "MB04PB.f"
	    i__3 = *n - i__ + 1;
#line 275 "MB04PB.f"
	    i__4 = i__ - 1;
#line 275 "MB04PB.f"
	    mb04pa_(&c_true, &i__3, &i__4, &ib, &a[i__ * a_dim1 + 1], lda, &
		    qg[i__ * qg_dim1 + 1], ldqg, &dwork[pxa], n, &dwork[pxg], 
		    n, &dwork[pxq], n, &dwork[pya], n, &cs[(i__ << 1) - 1], &
		    tau[i__], &dwork[pdw]);
#line 279 "MB04PB.f"
	    if (*n > i__ + ib) {

/*              Update the submatrix A(1:n,i+ib+1:n). */

#line 283 "MB04PB.f"
		i__3 = *n - i__ - ib;
#line 283 "MB04PB.f"
		i__4 = *n - i__ - ib;
#line 283 "MB04PB.f"
		dgemm_("No transpose", "Transpose", &i__3, &i__4, &ib, &c_b20,
			 &qg[i__ + ib + 1 + i__ * qg_dim1], ldqg, &dwork[pxa 
			+ ib + 1], n, &c_b20, &a[i__ + ib + 1 + (i__ + ib + 1)
			 * a_dim1], lda, (ftnlen)12, (ftnlen)9);
#line 286 "MB04PB.f"
		i__3 = *n - i__ - ib;
#line 286 "MB04PB.f"
		i__4 = *n - i__ - ib;
#line 286 "MB04PB.f"
		dgemm_("No transpose", "Transpose", &i__3, &i__4, &ib, &c_b20,
			 &a[i__ + ib + 1 + i__ * a_dim1], lda, &dwork[pxa + 
			nib + ib + 1], n, &c_b20, &a[i__ + ib + 1 + (i__ + ib 
			+ 1) * a_dim1], lda, (ftnlen)12, (ftnlen)9);
#line 290 "MB04PB.f"
		i__3 = *n - i__ - ib;
#line 290 "MB04PB.f"
		dgemm_("No transpose", "Transpose", n, &i__3, &ib, &c_b20, &
			dwork[pya], n, &qg[i__ + ib + 1 + i__ * qg_dim1], 
			ldqg, &c_b20, &a[(i__ + ib + 1) * a_dim1 + 1], lda, (
			ftnlen)12, (ftnlen)9);
#line 293 "MB04PB.f"
		i__3 = *n - i__ - ib;
#line 293 "MB04PB.f"
		dgemm_("No transpose", "Transpose", n, &i__3, &ib, &c_b20, &
			dwork[pya + nib], n, &a[i__ + ib + 1 + i__ * a_dim1], 
			lda, &c_b20, &a[(i__ + ib + 1) * a_dim1 + 1], lda, (
			ftnlen)12, (ftnlen)9);

/*              Update the submatrix Q(i+ib+1:n,i+ib+1:n). */

#line 299 "MB04PB.f"
		i__3 = *n - i__ - ib;
#line 299 "MB04PB.f"
		dsyr2k_("Lower", "No Transpose", &i__3, &ib, &c_b20, &dwork[
			pxq + ib + 1], n, &qg[i__ + ib + 1 + i__ * qg_dim1], 
			ldqg, &c_b20, &qg[i__ + ib + 1 + (i__ + ib + 1) * 
			qg_dim1], ldqg, (ftnlen)5, (ftnlen)12);
#line 302 "MB04PB.f"
		i__3 = *n - i__ - ib;
#line 302 "MB04PB.f"
		dsyr2k_("Lower", "No Transpose", &i__3, &ib, &c_b20, &dwork[
			pxq + nib + ib + 1], n, &a[i__ + ib + 1 + i__ * 
			a_dim1], lda, &c_b20, &qg[i__ + ib + 1 + (i__ + ib + 
			1) * qg_dim1], ldqg, (ftnlen)5, (ftnlen)12);

/*              Update the submatrix G(1:n,1:n). */

#line 308 "MB04PB.f"
		i__3 = i__ + ib;
#line 308 "MB04PB.f"
		i__4 = *n - i__ - ib;
#line 308 "MB04PB.f"
		dgemm_("No transpose", "Transpose", &i__3, &i__4, &ib, &c_b20,
			 &dwork[pxg], n, &qg[i__ + ib + 1 + i__ * qg_dim1], 
			ldqg, &c_b20, &qg[(i__ + ib + 2) * qg_dim1 + 1], ldqg,
			 (ftnlen)12, (ftnlen)9);
#line 311 "MB04PB.f"
		i__3 = i__ + ib;
#line 311 "MB04PB.f"
		i__4 = *n - i__ - ib;
#line 311 "MB04PB.f"
		dgemm_("No transpose", "Transpose", &i__3, &i__4, &ib, &c_b20,
			 &dwork[pxg + nib], n, &a[i__ + ib + 1 + i__ * a_dim1]
			, lda, &c_b20, &qg[(i__ + ib + 2) * qg_dim1 + 1], 
			ldqg, (ftnlen)12, (ftnlen)9);
#line 314 "MB04PB.f"
		i__3 = *n - i__ - ib;
#line 314 "MB04PB.f"
		dsyr2k_("Upper", "No Transpose", &i__3, &ib, &c_b20, &dwork[
			pxg + ib + i__], n, &qg[i__ + ib + 1 + i__ * qg_dim1],
			 ldqg, &c_b20, &qg[i__ + ib + 1 + (i__ + ib + 2) * 
			qg_dim1], ldqg, (ftnlen)5, (ftnlen)12);
#line 317 "MB04PB.f"
		i__3 = *n - i__ - ib;
#line 317 "MB04PB.f"
		dsyr2k_("Upper", "No Transpose", &i__3, &ib, &c_b20, &dwork[
			pxg + nib + ib + i__], n, &a[i__ + ib + 1 + i__ * 
			a_dim1], lda, &c_b20, &qg[i__ + ib + 1 + (i__ + ib + 
			2) * qg_dim1], ldqg, (ftnlen)5, (ftnlen)12);
#line 320 "MB04PB.f"
	    }
#line 321 "MB04PB.f"
/* L20: */
#line 321 "MB04PB.f"
	}
#line 322 "MB04PB.f"
    }

/*     Unblocked code to reduce the rest of the matrices. */

#line 326 "MB04PB.f"
    mb04pu_(n, &i__, &a[a_offset], lda, &qg[qg_offset], ldqg, &cs[1], &tau[1],
	     &dwork[1], ldwork, &ierr);

#line 329 "MB04PB.f"
    dwork[1] = (doublereal) wrkopt;

#line 331 "MB04PB.f"
    return 0;
/* *** Last line of MB04PB *** */
} /* mb04pb_ */

