#line 1 "MB04PU.f"
/* MB04PU.f -- translated by f2c (version 20100827).
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

#line 1 "MB04PU.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b7 = 0.;
static doublereal c_b14 = -1.;

/* Subroutine */ int mb04pu_(integer *n, integer *ilo, doublereal *a, integer 
	*lda, doublereal *qg, integer *ldqg, doublereal *cs, doublereal *tau, 
	doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, qg_dim1, qg_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal c__;
    static integer i__;
    static doublereal s, mu, nu;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dsyr2_(char 
	    *, integer *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, ftnlen);
    static doublereal alpha;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen), daxpy_(integer *, doublereal *, doublereal 
	    *, integer *, doublereal *, integer *);
    static doublereal ttemp;
    extern /* Subroutine */ int dsymv_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), dlarfg_(integer *, doublereal *,
	     doublereal *, integer *, doublereal *), dlartg_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), xerbla_(
	    char *, integer *, ftnlen);


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
/*     Unblocked version. */

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
/*             value of LDWORK. */
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

/*     The algorithm requires 40/3 N**3 + O(N) floating point operations */
/*     and is strongly backward stable. */

/*     REFERENCES */

/*     [1] C. F. VAN LOAN: */
/*         A symplectic method for approximating all the eigenvalues of */
/*         a Hamiltonian matrix. */
/*         Linear Algebra and its Applications, 61, pp. 233-251, 1984. */

/*     CONTRIBUTORS */

/*     D. Kressner (Technical Univ. Berlin, Germany) and */
/*     P. Benner (Technical Univ. Chemnitz, Germany), December 2003. */

/*     REVISIONS */

/*     V. Sima, Nov. 2008 (SLICOT version of the HAPACK routine DHAPVL). */

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

#line 186 "MB04PU.f"
    /* Parameter adjustments */
#line 186 "MB04PU.f"
    a_dim1 = *lda;
#line 186 "MB04PU.f"
    a_offset = 1 + a_dim1;
#line 186 "MB04PU.f"
    a -= a_offset;
#line 186 "MB04PU.f"
    qg_dim1 = *ldqg;
#line 186 "MB04PU.f"
    qg_offset = 1 + qg_dim1;
#line 186 "MB04PU.f"
    qg -= qg_offset;
#line 186 "MB04PU.f"
    --cs;
#line 186 "MB04PU.f"
    --tau;
#line 186 "MB04PU.f"
    --dwork;
#line 186 "MB04PU.f"

#line 186 "MB04PU.f"
    /* Function Body */
#line 186 "MB04PU.f"
    *info = 0;
#line 187 "MB04PU.f"
    if (*n < 0) {
#line 188 "MB04PU.f"
	*info = -1;
#line 189 "MB04PU.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 190 "MB04PU.f"
	*info = -2;
#line 191 "MB04PU.f"
    } else if (*lda < max(1,*n)) {
#line 192 "MB04PU.f"
	*info = -4;
#line 193 "MB04PU.f"
    } else if (*ldqg < max(1,*n)) {
#line 194 "MB04PU.f"
	*info = -6;
#line 195 "MB04PU.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 195 "MB04PU.f"
	i__1 = 1, i__2 = *n - 1;
#line 195 "MB04PU.f"
	if (*ldwork < max(i__1,i__2)) {
/* Computing MAX */
#line 196 "MB04PU.f"
	    i__1 = 1, i__2 = *n - 1;
#line 196 "MB04PU.f"
	    dwork[1] = (doublereal) max(i__1,i__2);
#line 197 "MB04PU.f"
	    *info = -10;
#line 198 "MB04PU.f"
	}
#line 198 "MB04PU.f"
    }

/*     Return if there were illegal values. */

#line 202 "MB04PU.f"
    if (*info != 0) {
#line 203 "MB04PU.f"
	i__1 = -(*info);
#line 203 "MB04PU.f"
	xerbla_("MB04PU", &i__1, (ftnlen)6);
#line 204 "MB04PU.f"
	return 0;
#line 205 "MB04PU.f"
    }

/*     Quick return if possible. */

#line 209 "MB04PU.f"
    if (*n <= *ilo) {
#line 210 "MB04PU.f"
	dwork[1] = 1.;
#line 211 "MB04PU.f"
	return 0;
#line 212 "MB04PU.f"
    }

#line 214 "MB04PU.f"
    i__1 = *n - 1;
#line 214 "MB04PU.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {

/*        Generate elementary reflector H(i) to annihilate QG(i+2:n,i). */

#line 218 "MB04PU.f"
	alpha = qg[i__ + 1 + i__ * qg_dim1];
#line 219 "MB04PU.f"
	i__2 = *n - i__;
/* Computing MIN */
#line 219 "MB04PU.f"
	i__3 = i__ + 2;
#line 219 "MB04PU.f"
	dlarfg_(&i__2, &alpha, &qg[min(i__3,*n) + i__ * qg_dim1], &c__1, &nu);
#line 220 "MB04PU.f"
	if (nu != 0.) {
#line 221 "MB04PU.f"
	    qg[i__ + 1 + i__ * qg_dim1] = 1.;

/*           Apply H(i) from both sides to QG(i+1:n,i+1:n). */
/*           Compute  x := nu * QG(i+1:n,i+1:n) * v. */

#line 226 "MB04PU.f"
	    i__2 = *n - i__;
#line 226 "MB04PU.f"
	    dsymv_("Lower", &i__2, &nu, &qg[i__ + 1 + (i__ + 1) * qg_dim1], 
		    ldqg, &qg[i__ + 1 + i__ * qg_dim1], &c__1, &c_b7, &dwork[
		    1], &c__1, (ftnlen)5);

/*           Compute  w := x - 1/2 * nu * (x'*v) * v. */

#line 231 "MB04PU.f"
	    i__2 = *n - i__;
#line 231 "MB04PU.f"
	    mu = nu * -.5 * ddot_(&i__2, &dwork[1], &c__1, &qg[i__ + 1 + i__ *
		     qg_dim1], &c__1);
#line 232 "MB04PU.f"
	    i__2 = *n - i__;
#line 232 "MB04PU.f"
	    daxpy_(&i__2, &mu, &qg[i__ + 1 + i__ * qg_dim1], &c__1, &dwork[1],
		     &c__1);

/*           Apply the transformation as a rank-2 update: */
/*                QG := QG - v * w' - w * v'. */

#line 237 "MB04PU.f"
	    i__2 = *n - i__;
#line 237 "MB04PU.f"
	    dsyr2_("Lower", &i__2, &c_b14, &qg[i__ + 1 + i__ * qg_dim1], &
		    c__1, &dwork[1], &c__1, &qg[i__ + 1 + (i__ + 1) * qg_dim1]
		    , ldqg, (ftnlen)5);

/*           Apply H(i) from the right hand side to QG(1:i,i+2:n+1). */

#line 242 "MB04PU.f"
	    i__2 = *n - i__;
#line 242 "MB04PU.f"
	    dlarf_("Right", &i__, &i__2, &qg[i__ + 1 + i__ * qg_dim1], &c__1, 
		    &nu, &qg[(i__ + 2) * qg_dim1 + 1], ldqg, &dwork[1], (
		    ftnlen)5);

/*           Apply H(i) from both sides to QG(i+1:n,i+2:n+1). */
/*           Compute  x := nu * QG(i+1:n,i+2:n+1) * v. */

#line 248 "MB04PU.f"
	    i__2 = *n - i__;
#line 248 "MB04PU.f"
	    dsymv_("Upper", &i__2, &nu, &qg[i__ + 1 + (i__ + 2) * qg_dim1], 
		    ldqg, &qg[i__ + 1 + i__ * qg_dim1], &c__1, &c_b7, &dwork[
		    1], &c__1, (ftnlen)5);

/*           Compute  w := x - 1/2 * nu * (x'*v) * v. */

#line 253 "MB04PU.f"
	    i__2 = *n - i__;
#line 253 "MB04PU.f"
	    mu = nu * -.5 * ddot_(&i__2, &dwork[1], &c__1, &qg[i__ + 1 + i__ *
		     qg_dim1], &c__1);
#line 254 "MB04PU.f"
	    i__2 = *n - i__;
#line 254 "MB04PU.f"
	    daxpy_(&i__2, &mu, &qg[i__ + 1 + i__ * qg_dim1], &c__1, &dwork[1],
		     &c__1);

/*           Apply the transformation as a rank-2 update: */
/*              QG(i+1:n,i+2:n+1) := QG(i+1:n,i+2:n+1) - v * w' - w * v'. */

#line 259 "MB04PU.f"
	    i__2 = *n - i__;
#line 259 "MB04PU.f"
	    dsyr2_("Upper", &i__2, &c_b14, &qg[i__ + 1 + i__ * qg_dim1], &
		    c__1, &dwork[1], &c__1, &qg[i__ + 1 + (i__ + 2) * qg_dim1]
		    , ldqg, (ftnlen)5);

/*           Apply H(i) from the left hand side to A(i+1:n,i:n). */

#line 264 "MB04PU.f"
	    i__2 = *n - i__;
#line 264 "MB04PU.f"
	    i__3 = *n - i__ + 1;
#line 264 "MB04PU.f"
	    dlarf_("Left", &i__2, &i__3, &qg[i__ + 1 + i__ * qg_dim1], &c__1, 
		    &nu, &a[i__ + 1 + i__ * a_dim1], lda, &dwork[1], (ftnlen)
		    4);

/*           Apply H(i) from the right hand side to A(1:n,i+1:n). */

#line 269 "MB04PU.f"
	    i__2 = *n - i__;
#line 269 "MB04PU.f"
	    dlarf_("Right", n, &i__2, &qg[i__ + 1 + i__ * qg_dim1], &c__1, &
		    nu, &a[(i__ + 1) * a_dim1 + 1], lda, &dwork[1], (ftnlen)5)
		    ;
#line 271 "MB04PU.f"
	}
#line 272 "MB04PU.f"
	qg[i__ + 1 + i__ * qg_dim1] = nu;

/*        Generate symplectic Givens rotation G(i) to annihilate */
/*        QG(i+1,i). */

#line 277 "MB04PU.f"
	temp = a[i__ + 1 + i__ * a_dim1];
#line 278 "MB04PU.f"
	dlartg_(&temp, &alpha, &c__, &s, &a[i__ + 1 + i__ * a_dim1]);

/*        Apply G(i) to [A(I+1,I+2:N); QG(I+2:N,I+1)']. */

#line 282 "MB04PU.f"
	i__2 = *n - i__ - 1;
#line 282 "MB04PU.f"
	drot_(&i__2, &a[i__ + 1 + (i__ + 2) * a_dim1], lda, &qg[i__ + 2 + (
		i__ + 1) * qg_dim1], &c__1, &c__, &s);

/*        Apply G(i) to [A(1:I,I+1) QG(1:I,I+2)]. */

#line 286 "MB04PU.f"
	drot_(&i__, &a[(i__ + 1) * a_dim1 + 1], &c__1, &qg[(i__ + 2) * 
		qg_dim1 + 1], &c__1, &c__, &s);

/*        Apply G(i) to [A(I+2:N,I+1) QG(I+1, I+3:N+1)'] from the right. */

#line 290 "MB04PU.f"
	i__2 = *n - i__ - 1;
#line 290 "MB04PU.f"
	drot_(&i__2, &a[i__ + 2 + (i__ + 1) * a_dim1], &c__1, &qg[i__ + 1 + (
		i__ + 3) * qg_dim1], ldqg, &c__, &s);

/*        Fix the diagonal part. */

#line 294 "MB04PU.f"
	temp = a[i__ + 1 + (i__ + 1) * a_dim1];
#line 295 "MB04PU.f"
	ttemp = qg[i__ + 1 + (i__ + 2) * qg_dim1];
#line 296 "MB04PU.f"
	a[i__ + 1 + (i__ + 1) * a_dim1] = c__ * temp + s * qg[i__ + 1 + (i__ 
		+ 1) * qg_dim1];
#line 297 "MB04PU.f"
	qg[i__ + 1 + (i__ + 2) * qg_dim1] = c__ * ttemp - s * temp;
#line 298 "MB04PU.f"
	qg[i__ + 1 + (i__ + 1) * qg_dim1] = -s * temp + c__ * qg[i__ + 1 + (
		i__ + 1) * qg_dim1];
#line 299 "MB04PU.f"
	ttemp = -s * ttemp - c__ * temp;
#line 300 "MB04PU.f"
	temp = a[i__ + 1 + (i__ + 1) * a_dim1];
#line 301 "MB04PU.f"
	qg[i__ + 1 + (i__ + 1) * qg_dim1] = c__ * qg[i__ + 1 + (i__ + 1) * 
		qg_dim1] + s * ttemp;
#line 302 "MB04PU.f"
	a[i__ + 1 + (i__ + 1) * a_dim1] = c__ * temp + s * qg[i__ + 1 + (i__ 
		+ 2) * qg_dim1];
#line 303 "MB04PU.f"
	qg[i__ + 1 + (i__ + 2) * qg_dim1] = -s * temp + c__ * qg[i__ + 1 + (
		i__ + 2) * qg_dim1];
#line 304 "MB04PU.f"
	cs[(i__ << 1) - 1] = c__;
#line 305 "MB04PU.f"
	cs[i__ * 2] = s;

/*        Generate elementary reflector F(i) to annihilate A(i+2:n,i). */

#line 309 "MB04PU.f"
	i__2 = *n - i__;
/* Computing MIN */
#line 309 "MB04PU.f"
	i__3 = i__ + 2;
#line 309 "MB04PU.f"
	dlarfg_(&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min(i__3,*n) + i__ * 
		a_dim1], &c__1, &nu);
#line 310 "MB04PU.f"
	if (nu != 0.) {
#line 311 "MB04PU.f"
	    temp = a[i__ + 1 + i__ * a_dim1];
#line 312 "MB04PU.f"
	    a[i__ + 1 + i__ * a_dim1] = 1.;

/*           Apply F(i) from the left hand side to A(i+1:n,i+1:n). */

#line 316 "MB04PU.f"
	    i__2 = *n - i__;
#line 316 "MB04PU.f"
	    i__3 = *n - i__;
#line 316 "MB04PU.f"
	    dlarf_("Left", &i__2, &i__3, &a[i__ + 1 + i__ * a_dim1], &c__1, &
		    nu, &a[i__ + 1 + (i__ + 1) * a_dim1], lda, &dwork[1], (
		    ftnlen)4);

/*           Apply G(i) from the right hand side to A(1:n,i+1:n). */

#line 321 "MB04PU.f"
	    i__2 = *n - i__;
#line 321 "MB04PU.f"
	    dlarf_("Right", n, &i__2, &a[i__ + 1 + i__ * a_dim1], &c__1, &nu, 
		    &a[(i__ + 1) * a_dim1 + 1], lda, &dwork[1], (ftnlen)5);

/*           Apply G(i) from both sides to QG(i+1:n,i+1:n). */
/*           Compute  x := nu * QG(i+1:n,i+1:n) * v. */

#line 327 "MB04PU.f"
	    i__2 = *n - i__;
#line 327 "MB04PU.f"
	    dsymv_("Lower", &i__2, &nu, &qg[i__ + 1 + (i__ + 1) * qg_dim1], 
		    ldqg, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b7, &dwork[1],
		     &c__1, (ftnlen)5);

/*           Compute  w := x - 1/2 * tau * (x'*v) * v. */

#line 332 "MB04PU.f"
	    i__2 = *n - i__;
#line 332 "MB04PU.f"
	    mu = nu * -.5 * ddot_(&i__2, &dwork[1], &c__1, &a[i__ + 1 + i__ * 
		    a_dim1], &c__1);
#line 333 "MB04PU.f"
	    i__2 = *n - i__;
#line 333 "MB04PU.f"
	    daxpy_(&i__2, &mu, &a[i__ + 1 + i__ * a_dim1], &c__1, &dwork[1], &
		    c__1);

/*           Apply the transformation as a rank-2 update: */
/*                QG := QG - v * w' - w * v'. */

#line 338 "MB04PU.f"
	    i__2 = *n - i__;
#line 338 "MB04PU.f"
	    dsyr2_("Lower", &i__2, &c_b14, &a[i__ + 1 + i__ * a_dim1], &c__1, 
		    &dwork[1], &c__1, &qg[i__ + 1 + (i__ + 1) * qg_dim1], 
		    ldqg, (ftnlen)5);

/*           Apply G(i) from the right hand side to QG(1:i,i+2:n+1). */

#line 343 "MB04PU.f"
	    i__2 = *n - i__;
#line 343 "MB04PU.f"
	    dlarf_("Right", &i__, &i__2, &a[i__ + 1 + i__ * a_dim1], &c__1, &
		    nu, &qg[(i__ + 2) * qg_dim1 + 1], ldqg, &dwork[1], (
		    ftnlen)5);

/*           Apply G(i) from both sides to QG(i+1:n,i+2:n+1). */
/*           Compute  x := nu * QG(i+1:n,i+2:n+1) * v. */

#line 349 "MB04PU.f"
	    i__2 = *n - i__;
#line 349 "MB04PU.f"
	    dsymv_("Upper", &i__2, &nu, &qg[i__ + 1 + (i__ + 2) * qg_dim1], 
		    ldqg, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b7, &dwork[1],
		     &c__1, (ftnlen)5);

/*           Compute  w := x - 1/2 * tau * (x'*v) * v. */

#line 354 "MB04PU.f"
	    i__2 = *n - i__;
#line 354 "MB04PU.f"
	    mu = nu * -.5 * ddot_(&i__2, &dwork[1], &c__1, &a[i__ + 1 + i__ * 
		    a_dim1], &c__1);
#line 355 "MB04PU.f"
	    i__2 = *n - i__;
#line 355 "MB04PU.f"
	    daxpy_(&i__2, &mu, &a[i__ + 1 + i__ * a_dim1], &c__1, &dwork[1], &
		    c__1);

/*           Apply the transformation as a rank-2 update: */
/*              QG(i+1:n,i+2:n+1) := QG(i+1:n,i+2:n+1) - v * w' - w * v'. */

#line 360 "MB04PU.f"
	    i__2 = *n - i__;
#line 360 "MB04PU.f"
	    dsyr2_("Upper", &i__2, &c_b14, &a[i__ + 1 + i__ * a_dim1], &c__1, 
		    &dwork[1], &c__1, &qg[i__ + 1 + (i__ + 2) * qg_dim1], 
		    ldqg, (ftnlen)5);
#line 362 "MB04PU.f"
	    a[i__ + 1 + i__ * a_dim1] = temp;
#line 363 "MB04PU.f"
	}
#line 364 "MB04PU.f"
	tau[i__] = nu;
#line 365 "MB04PU.f"
/* L10: */
#line 365 "MB04PU.f"
    }
/* Computing MAX */
#line 366 "MB04PU.f"
    i__1 = 1, i__2 = *n - 1;
#line 366 "MB04PU.f"
    dwork[1] = (doublereal) max(i__1,i__2);
#line 367 "MB04PU.f"
    return 0;
/* *** Last line of MB04PU *** */
} /* mb04pu_ */

