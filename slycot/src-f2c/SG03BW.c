#line 1 "SG03BW.f"
/* SG03BW.f -- translated by f2c (version 20100827).
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

#line 1 "SG03BW.f"
/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;
static doublereal c_b18 = 1.;
static doublereal c_b19 = 0.;
static integer c__2 = 2;
static doublereal c_b23 = -1.;

/* Subroutine */ int sg03bw_(char *trans, integer *m, integer *n, doublereal *
	a, integer *lda, doublereal *c__, integer *ldc, doublereal *e, 
	integer *lde, doublereal *d__, integer *ldd, doublereal *x, integer *
	ldx, doublereal *scale, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, d_dim1, d_offset, e_dim1, 
	    e_offset, x_dim1, x_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ma, mb, me;
    static doublereal tm[4]	/* was [2][2] */;
    static integer mai, maj;
    static doublereal mat[16]	/* was [4][4] */, rhs[4];
    static integer piv1[4], piv2[4], info1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb02uu_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *), mb02uv_(
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *);
    static doublereal scale1;
    static integer dimmat;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical notrns;


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

/*     To solve for X the generalized Sylvester equation */

/*         T            T */
/*        A  * X * C + E  * X * D  =  SCALE * Y,                      (1) */

/*     or the transposed equation */

/*                 T            T */
/*        A * X * C  + E * X * D   =  SCALE * Y,                      (2) */

/*     where A and E are real M-by-M matrices, C and D are real N-by-N */
/*     matrices, X and Y are real M-by-N matrices. N is either 1 or 2. */
/*     The pencil A - lambda * E must be in generalized real Schur form */
/*     (A upper quasitriangular, E upper triangular). SCALE is an output */
/*     scale factor, set to avoid overflow in X. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANS   CHARACTER*1 */
/*             Specifies whether the transposed equation is to be solved */
/*             or not: */
/*             = 'N':  Solve equation (1); */
/*             = 'T':  Solve equation (2). */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of the matrices A and E.  M >= 0. */

/*     N       (input) INTEGER */
/*             The order of the matrices C and D.  N = 1 or N = 2. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,M) */
/*             The leading M-by-M part of this array must contain the */
/*             upper quasitriangular matrix A. The elements below the */
/*             upper Hessenberg part are not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,M). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading N-by-N part of this array must contain the */
/*             matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= MAX(1,N). */

/*     E       (input) DOUBLE PRECISION array, dimension (LDE,M) */
/*             The leading M-by-M part of this array must contain the */
/*             upper triangular matrix E. The elements below the main */
/*             diagonal are not referenced. */

/*     LDE     INTEGER */
/*             The leading dimension of the array E.  LDE >= MAX(1,M). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,N) */
/*             The leading N-by-N part of this array must contain the */
/*             matrix D. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= MAX(1,N). */

/*     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the right hand side matrix Y. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the solution matrix X. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X.  LDX >= MAX(1,M). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor set to avoid overflow in X. */
/*             0 < SCALE <= 1. */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the generalized Sylvester equation is (nearly) */
/*                   singular to working precision;  perturbed values */
/*                   were used to solve the equation (but the matrices */
/*                   A, C, D, and E are unchanged). */

/*     METHOD */

/*     The method used by the routine is based on a generalization of the */
/*     algorithm due to Bartels and Stewart [1]. See also [2] and [3] for */
/*     details. */

/*     REFERENCES */

/*     [1] Bartels, R.H., Stewart, G.W. */
/*         Solution of the equation A X + X B = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     [2] Gardiner, J.D., Laub, A.J., Amato, J.J., Moler, C.B. */
/*         Solution of the Sylvester Matrix Equation */
/*         A X B**T + C X D**T = E. */
/*         A.C.M. Trans. Math. Soft., vol. 18, no. 2, pp. 223-231, 1992. */

/*     [3] Penzl, T. */
/*         Numerical solution of generalized Lyapunov equations. */
/*         Advances in Comp. Math., vol. 8, pp. 33-48, 1998. */

/*     NUMERICAL ASPECTS */

/*     The routine requires about 2 * N * M**2 flops. Note that we count */
/*     a single floating point arithmetic operation as one flop. */

/*     The algorithm is backward stable if the eigenvalues of the pencil */
/*     A - lambda * E are real. Otherwise, linear systems of order at */
/*     most 4 are involved into the computation. These systems are solved */
/*     by Gauss elimination with complete pivoting. The loss of stability */
/*     of the Gauss elimination with complete pivoting is rarely */
/*     encountered in practice. */

/*     FURTHER COMMENTS */

/*     When near singularity is detected, perturbed values are used */
/*     to solve the equation (but the given matrices are unchanged). */

/*     CONTRIBUTOR */

/*     T. Penzl, Technical University Chemnitz, Germany, Aug. 1998. */

/*     REVISIONS */

/*     Sep. 1998 (V. Sima). */

/*     KEYWORDS */

/*     Lyapunov equation */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     Decode input parameters. */

#line 190 "SG03BW.f"
    /* Parameter adjustments */
#line 190 "SG03BW.f"
    a_dim1 = *lda;
#line 190 "SG03BW.f"
    a_offset = 1 + a_dim1;
#line 190 "SG03BW.f"
    a -= a_offset;
#line 190 "SG03BW.f"
    c_dim1 = *ldc;
#line 190 "SG03BW.f"
    c_offset = 1 + c_dim1;
#line 190 "SG03BW.f"
    c__ -= c_offset;
#line 190 "SG03BW.f"
    e_dim1 = *lde;
#line 190 "SG03BW.f"
    e_offset = 1 + e_dim1;
#line 190 "SG03BW.f"
    e -= e_offset;
#line 190 "SG03BW.f"
    d_dim1 = *ldd;
#line 190 "SG03BW.f"
    d_offset = 1 + d_dim1;
#line 190 "SG03BW.f"
    d__ -= d_offset;
#line 190 "SG03BW.f"
    x_dim1 = *ldx;
#line 190 "SG03BW.f"
    x_offset = 1 + x_dim1;
#line 190 "SG03BW.f"
    x -= x_offset;
#line 190 "SG03BW.f"

#line 190 "SG03BW.f"
    /* Function Body */
#line 190 "SG03BW.f"
    notrns = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

#line 194 "SG03BW.f"
    if (! (notrns || lsame_(trans, "T", (ftnlen)1, (ftnlen)1))) {
#line 195 "SG03BW.f"
	*info = -1;
#line 196 "SG03BW.f"
    } else if (*m < 0) {
#line 197 "SG03BW.f"
	*info = -2;
#line 198 "SG03BW.f"
    } else if (*n != 1 && *n != 2) {
#line 199 "SG03BW.f"
	*info = -3;
#line 200 "SG03BW.f"
    } else if (*lda < max(1,*m)) {
#line 201 "SG03BW.f"
	*info = -5;
#line 202 "SG03BW.f"
    } else if (*ldc < max(1,*n)) {
#line 203 "SG03BW.f"
	*info = -7;
#line 204 "SG03BW.f"
    } else if (*lde < max(1,*m)) {
#line 205 "SG03BW.f"
	*info = -9;
#line 206 "SG03BW.f"
    } else if (*ldd < max(1,*n)) {
#line 207 "SG03BW.f"
	*info = -11;
#line 208 "SG03BW.f"
    } else if (*ldx < max(1,*m)) {
#line 209 "SG03BW.f"
	*info = -13;
#line 210 "SG03BW.f"
    } else {
#line 211 "SG03BW.f"
	*info = 0;
#line 212 "SG03BW.f"
    }
#line 213 "SG03BW.f"
    if (*info != 0) {
#line 214 "SG03BW.f"
	i__1 = -(*info);
#line 214 "SG03BW.f"
	xerbla_("SG03BW", &i__1, (ftnlen)6);
#line 215 "SG03BW.f"
	return 0;
#line 216 "SG03BW.f"
    }

#line 218 "SG03BW.f"
    *scale = 1.;

/*     Quick return if possible. */

#line 222 "SG03BW.f"
    if (*m == 0) {
#line 222 "SG03BW.f"
	return 0;
#line 222 "SG03BW.f"
    }

#line 225 "SG03BW.f"
    if (notrns) {

/*        Solve equation (1). */

/*        Compute block row X(MA:ME,:). MB denotes the number of rows in */
/*        this block row. */

#line 232 "SG03BW.f"
	me = 0;
/*        WHILE ( ME .NE. M ) DO */
#line 234 "SG03BW.f"
L20:
#line 234 "SG03BW.f"
	if (me != *m) {
#line 235 "SG03BW.f"
	    ma = me + 1;
#line 236 "SG03BW.f"
	    if (ma == *m) {
#line 237 "SG03BW.f"
		me = *m;
#line 238 "SG03BW.f"
		mb = 1;
#line 239 "SG03BW.f"
	    } else {
#line 240 "SG03BW.f"
		if (a[ma + 1 + ma * a_dim1] == 0.) {
#line 241 "SG03BW.f"
		    me = ma;
#line 242 "SG03BW.f"
		    mb = 1;
#line 243 "SG03BW.f"
		} else {
#line 244 "SG03BW.f"
		    me = ma + 1;
#line 245 "SG03BW.f"
		    mb = 2;
#line 246 "SG03BW.f"
		}
#line 247 "SG03BW.f"
	    }

/*           Assemble Kronecker product system of linear equations with */
/*           matrix */

/*              MAT = kron(C',A(MA:ME,MA:ME)') + kron(D',E(MA:ME,MA:ME)') */

/*           and right hand side */

/*              RHS = vec(X(MA:ME,:)). */

#line 258 "SG03BW.f"
	    if (*n == 1) {
#line 259 "SG03BW.f"
		dimmat = mb;
#line 260 "SG03BW.f"
		i__1 = mb;
#line 260 "SG03BW.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 261 "SG03BW.f"
		    mai = ma + i__ - 1;
#line 262 "SG03BW.f"
		    i__2 = mb;
#line 262 "SG03BW.f"
		    for (j = 1; j <= i__2; ++j) {
#line 263 "SG03BW.f"
			maj = ma + j - 1;
#line 264 "SG03BW.f"
			mat[i__ + (j << 2) - 5] = c__[c_dim1 + 1] * a[maj + 
				mai * a_dim1];
#line 265 "SG03BW.f"
			if (maj <= mai) {
#line 265 "SG03BW.f"
			    mat[i__ + (j << 2) - 5] += d__[d_dim1 + 1] * e[
				    maj + mai * e_dim1];
#line 265 "SG03BW.f"
			}
#line 267 "SG03BW.f"
/* L40: */
#line 267 "SG03BW.f"
		    }
#line 268 "SG03BW.f"
		    rhs[i__ - 1] = x[mai + x_dim1];
#line 269 "SG03BW.f"
/* L60: */
#line 269 "SG03BW.f"
		}
#line 270 "SG03BW.f"
	    } else {
#line 271 "SG03BW.f"
		dimmat = mb << 1;
#line 272 "SG03BW.f"
		i__1 = mb;
#line 272 "SG03BW.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 273 "SG03BW.f"
		    mai = ma + i__ - 1;
#line 274 "SG03BW.f"
		    i__2 = mb;
#line 274 "SG03BW.f"
		    for (j = 1; j <= i__2; ++j) {
#line 275 "SG03BW.f"
			maj = ma + j - 1;
#line 276 "SG03BW.f"
			mat[i__ + (j << 2) - 5] = c__[c_dim1 + 1] * a[maj + 
				mai * a_dim1];
#line 277 "SG03BW.f"
			mat[mb + i__ + (j << 2) - 5] = c__[(c_dim1 << 1) + 1] 
				* a[maj + mai * a_dim1];
#line 278 "SG03BW.f"
			mat[i__ + (mb + j << 2) - 5] = c__[c_dim1 + 2] * a[
				maj + mai * a_dim1];
#line 279 "SG03BW.f"
			mat[mb + i__ + (mb + j << 2) - 5] = c__[(c_dim1 << 1) 
				+ 2] * a[maj + mai * a_dim1];
#line 280 "SG03BW.f"
			if (maj <= mai) {
#line 281 "SG03BW.f"
			    mat[i__ + (j << 2) - 5] += d__[d_dim1 + 1] * e[
				    maj + mai * e_dim1];
#line 282 "SG03BW.f"
			    mat[mb + i__ + (j << 2) - 5] += d__[(d_dim1 << 1) 
				    + 1] * e[maj + mai * e_dim1];
#line 283 "SG03BW.f"
			    mat[i__ + (mb + j << 2) - 5] += d__[d_dim1 + 2] * 
				    e[maj + mai * e_dim1];
#line 284 "SG03BW.f"
			    mat[mb + i__ + (mb + j << 2) - 5] += d__[(d_dim1 
				    << 1) + 2] * e[maj + mai * e_dim1];
#line 286 "SG03BW.f"
			}
#line 287 "SG03BW.f"
/* L80: */
#line 287 "SG03BW.f"
		    }
#line 288 "SG03BW.f"
		    rhs[i__ - 1] = x[mai + x_dim1];
#line 289 "SG03BW.f"
		    rhs[mb + i__ - 1] = x[mai + (x_dim1 << 1)];
#line 290 "SG03BW.f"
/* L100: */
#line 290 "SG03BW.f"
		}
#line 291 "SG03BW.f"
	    }

/*           Solve the system of linear equations. */

#line 295 "SG03BW.f"
	    mb02uv_(&dimmat, mat, &c__4, piv1, piv2, &info1);
#line 296 "SG03BW.f"
	    if (info1 != 0) {
#line 296 "SG03BW.f"
		*info = 1;
#line 296 "SG03BW.f"
	    }
#line 298 "SG03BW.f"
	    mb02uu_(&dimmat, mat, &c__4, rhs, piv1, piv2, &scale1);
#line 299 "SG03BW.f"
	    if (scale1 != 1.) {
#line 300 "SG03BW.f"
		*scale = scale1 * *scale;
#line 301 "SG03BW.f"
		i__1 = *n;
#line 301 "SG03BW.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 302 "SG03BW.f"
		    dscal_(m, &scale1, &x[i__ * x_dim1 + 1], &c__1);
#line 303 "SG03BW.f"
/* L120: */
#line 303 "SG03BW.f"
		}
#line 304 "SG03BW.f"
	    }

#line 306 "SG03BW.f"
	    if (*n == 1) {
#line 307 "SG03BW.f"
		i__1 = mb;
#line 307 "SG03BW.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 308 "SG03BW.f"
		    mai = ma + i__ - 1;
#line 309 "SG03BW.f"
		    x[mai + x_dim1] = rhs[i__ - 1];
#line 310 "SG03BW.f"
/* L140: */
#line 310 "SG03BW.f"
		}
#line 311 "SG03BW.f"
	    } else {
#line 312 "SG03BW.f"
		i__1 = mb;
#line 312 "SG03BW.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 313 "SG03BW.f"
		    mai = ma + i__ - 1;
#line 314 "SG03BW.f"
		    x[mai + x_dim1] = rhs[i__ - 1];
#line 315 "SG03BW.f"
		    x[mai + (x_dim1 << 1)] = rhs[mb + i__ - 1];
#line 316 "SG03BW.f"
/* L160: */
#line 316 "SG03BW.f"
		}
#line 317 "SG03BW.f"
	    }

/*           Update right hand sides. */

/*           X(ME+1:M,:) = X(ME+1:M,:) - A(MA:ME,ME+1:M)'*X(MA:ME,:)*C */

/*           X(ME+1:M,:) = X(ME+1:M,:) - E(MA:ME,ME+1:M)'*X(MA:ME,:)*D */

#line 325 "SG03BW.f"
	    if (me < *m) {
#line 326 "SG03BW.f"
		dgemm_("N", "N", &mb, n, n, &c_b18, &x[ma + x_dim1], ldx, &
			c__[c_offset], ldc, &c_b19, tm, &c__2, (ftnlen)1, (
			ftnlen)1);
#line 328 "SG03BW.f"
		i__1 = *m - me;
#line 328 "SG03BW.f"
		dgemm_("T", "N", &i__1, n, &mb, &c_b23, &a[ma + (me + 1) * 
			a_dim1], lda, tm, &c__2, &c_b18, &x[me + 1 + x_dim1], 
			ldx, (ftnlen)1, (ftnlen)1);
#line 330 "SG03BW.f"
		dgemm_("N", "N", &mb, n, n, &c_b18, &x[ma + x_dim1], ldx, &
			d__[d_offset], ldd, &c_b19, tm, &c__2, (ftnlen)1, (
			ftnlen)1);
#line 332 "SG03BW.f"
		i__1 = *m - me;
#line 332 "SG03BW.f"
		dgemm_("T", "N", &i__1, n, &mb, &c_b23, &e[ma + (me + 1) * 
			e_dim1], lde, tm, &c__2, &c_b18, &x[me + 1 + x_dim1], 
			ldx, (ftnlen)1, (ftnlen)1);
#line 334 "SG03BW.f"
	    }

#line 336 "SG03BW.f"
	    goto L20;
#line 337 "SG03BW.f"
	}
/*        END WHILE 20 */

#line 340 "SG03BW.f"
    } else {

/*        Solve equation (2). */

/*        Compute block row X(MA:ME,:). MB denotes the number of rows in */
/*        this block row. */

#line 347 "SG03BW.f"
	ma = *m + 1;
/*        WHILE ( MA .NE. 1 ) DO */
#line 349 "SG03BW.f"
L180:
#line 349 "SG03BW.f"
	if (ma != 1) {
#line 350 "SG03BW.f"
	    me = ma - 1;
#line 351 "SG03BW.f"
	    if (me == 1) {
#line 352 "SG03BW.f"
		ma = 1;
#line 353 "SG03BW.f"
		mb = 1;
#line 354 "SG03BW.f"
	    } else {
#line 355 "SG03BW.f"
		if (a[me + (me - 1) * a_dim1] == 0.) {
#line 356 "SG03BW.f"
		    ma = me;
#line 357 "SG03BW.f"
		    mb = 1;
#line 358 "SG03BW.f"
		} else {
#line 359 "SG03BW.f"
		    ma = me - 1;
#line 360 "SG03BW.f"
		    mb = 2;
#line 361 "SG03BW.f"
		}
#line 362 "SG03BW.f"
	    }

/*           Assemble Kronecker product system of linear equations with */
/*           matrix */

/*              MAT = kron(C,A(MA:ME,MA:ME)) + kron(D,E(MA:ME,MA:ME)) */

/*           and right hand side */

/*              RHS = vec(X(MA:ME,:)). */

#line 373 "SG03BW.f"
	    if (*n == 1) {
#line 374 "SG03BW.f"
		dimmat = mb;
#line 375 "SG03BW.f"
		i__1 = mb;
#line 375 "SG03BW.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 376 "SG03BW.f"
		    mai = ma + i__ - 1;
#line 377 "SG03BW.f"
		    i__2 = mb;
#line 377 "SG03BW.f"
		    for (j = 1; j <= i__2; ++j) {
#line 378 "SG03BW.f"
			maj = ma + j - 1;
#line 379 "SG03BW.f"
			mat[i__ + (j << 2) - 5] = c__[c_dim1 + 1] * a[mai + 
				maj * a_dim1];
#line 380 "SG03BW.f"
			if (maj >= mai) {
#line 380 "SG03BW.f"
			    mat[i__ + (j << 2) - 5] += d__[d_dim1 + 1] * e[
				    mai + maj * e_dim1];
#line 380 "SG03BW.f"
			}
#line 382 "SG03BW.f"
/* L200: */
#line 382 "SG03BW.f"
		    }
#line 383 "SG03BW.f"
		    rhs[i__ - 1] = x[mai + x_dim1];
#line 384 "SG03BW.f"
/* L220: */
#line 384 "SG03BW.f"
		}
#line 385 "SG03BW.f"
	    } else {
#line 386 "SG03BW.f"
		dimmat = mb << 1;
#line 387 "SG03BW.f"
		i__1 = mb;
#line 387 "SG03BW.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 388 "SG03BW.f"
		    mai = ma + i__ - 1;
#line 389 "SG03BW.f"
		    i__2 = mb;
#line 389 "SG03BW.f"
		    for (j = 1; j <= i__2; ++j) {
#line 390 "SG03BW.f"
			maj = ma + j - 1;
#line 391 "SG03BW.f"
			mat[i__ + (j << 2) - 5] = c__[c_dim1 + 1] * a[mai + 
				maj * a_dim1];
#line 392 "SG03BW.f"
			mat[mb + i__ + (j << 2) - 5] = c__[c_dim1 + 2] * a[
				mai + maj * a_dim1];
#line 393 "SG03BW.f"
			mat[i__ + (mb + j << 2) - 5] = c__[(c_dim1 << 1) + 1] 
				* a[mai + maj * a_dim1];
#line 394 "SG03BW.f"
			mat[mb + i__ + (mb + j << 2) - 5] = c__[(c_dim1 << 1) 
				+ 2] * a[mai + maj * a_dim1];
#line 395 "SG03BW.f"
			if (maj >= mai) {
#line 396 "SG03BW.f"
			    mat[i__ + (j << 2) - 5] += d__[d_dim1 + 1] * e[
				    mai + maj * e_dim1];
#line 397 "SG03BW.f"
			    mat[mb + i__ + (j << 2) - 5] += d__[d_dim1 + 2] * 
				    e[mai + maj * e_dim1];
#line 398 "SG03BW.f"
			    mat[i__ + (mb + j << 2) - 5] += d__[(d_dim1 << 1) 
				    + 1] * e[mai + maj * e_dim1];
#line 399 "SG03BW.f"
			    mat[mb + i__ + (mb + j << 2) - 5] += d__[(d_dim1 
				    << 1) + 2] * e[mai + maj * e_dim1];
#line 401 "SG03BW.f"
			}
#line 402 "SG03BW.f"
/* L240: */
#line 402 "SG03BW.f"
		    }
#line 403 "SG03BW.f"
		    rhs[i__ - 1] = x[mai + x_dim1];
#line 404 "SG03BW.f"
		    rhs[mb + i__ - 1] = x[mai + (x_dim1 << 1)];
#line 405 "SG03BW.f"
/* L260: */
#line 405 "SG03BW.f"
		}
#line 406 "SG03BW.f"
	    }

/*           Solve the system of linear equations. */

#line 410 "SG03BW.f"
	    mb02uv_(&dimmat, mat, &c__4, piv1, piv2, &info1);
#line 411 "SG03BW.f"
	    if (info1 != 0) {
#line 411 "SG03BW.f"
		*info = 1;
#line 411 "SG03BW.f"
	    }
#line 413 "SG03BW.f"
	    mb02uu_(&dimmat, mat, &c__4, rhs, piv1, piv2, &scale1);
#line 414 "SG03BW.f"
	    if (scale1 != 1.) {
#line 415 "SG03BW.f"
		*scale = scale1 * *scale;
#line 416 "SG03BW.f"
		i__1 = *n;
#line 416 "SG03BW.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 417 "SG03BW.f"
		    dscal_(m, &scale1, &x[i__ * x_dim1 + 1], &c__1);
#line 418 "SG03BW.f"
/* L280: */
#line 418 "SG03BW.f"
		}
#line 419 "SG03BW.f"
	    }

#line 421 "SG03BW.f"
	    if (*n == 1) {
#line 422 "SG03BW.f"
		i__1 = mb;
#line 422 "SG03BW.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 423 "SG03BW.f"
		    mai = ma + i__ - 1;
#line 424 "SG03BW.f"
		    x[mai + x_dim1] = rhs[i__ - 1];
#line 425 "SG03BW.f"
/* L300: */
#line 425 "SG03BW.f"
		}
#line 426 "SG03BW.f"
	    } else {
#line 427 "SG03BW.f"
		i__1 = mb;
#line 427 "SG03BW.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 428 "SG03BW.f"
		    mai = ma + i__ - 1;
#line 429 "SG03BW.f"
		    x[mai + x_dim1] = rhs[i__ - 1];
#line 430 "SG03BW.f"
		    x[mai + (x_dim1 << 1)] = rhs[mb + i__ - 1];
#line 431 "SG03BW.f"
/* L320: */
#line 431 "SG03BW.f"
		}
#line 432 "SG03BW.f"
	    }

/*           Update right hand sides. */

/*              X(1:MA-1,:) = X(1:MA-1,:) - A(1:MA-1,MA:ME)*X(MA:ME,:)*C' */

/*              X(1:MA-1,:) = X(1:MA-1,:) - E(1:MA-1,MA:ME)*X(MA:ME,:)*D' */

#line 440 "SG03BW.f"
	    if (ma > 1) {
#line 441 "SG03BW.f"
		dgemm_("N", "T", &mb, n, n, &c_b18, &x[ma + x_dim1], ldx, &
			c__[c_offset], ldc, &c_b19, tm, &c__2, (ftnlen)1, (
			ftnlen)1);
#line 443 "SG03BW.f"
		i__1 = ma - 1;
#line 443 "SG03BW.f"
		dgemm_("N", "N", &i__1, n, &mb, &c_b23, &a[ma * a_dim1 + 1], 
			lda, tm, &c__2, &c_b18, &x[x_offset], ldx, (ftnlen)1, 
			(ftnlen)1);
#line 445 "SG03BW.f"
		dgemm_("N", "T", &mb, n, n, &c_b18, &x[ma + x_dim1], ldx, &
			d__[d_offset], ldd, &c_b19, tm, &c__2, (ftnlen)1, (
			ftnlen)1);
#line 447 "SG03BW.f"
		i__1 = ma - 1;
#line 447 "SG03BW.f"
		dgemm_("N", "N", &i__1, n, &mb, &c_b23, &e[ma * e_dim1 + 1], 
			lde, tm, &c__2, &c_b18, &x[x_offset], ldx, (ftnlen)1, 
			(ftnlen)1);
#line 449 "SG03BW.f"
	    }

#line 451 "SG03BW.f"
	    goto L180;
#line 452 "SG03BW.f"
	}
/*        END WHILE 180 */

#line 455 "SG03BW.f"
    }

#line 457 "SG03BW.f"
    return 0;
/* *** Last line of SG03BW *** */
} /* sg03bw_ */

