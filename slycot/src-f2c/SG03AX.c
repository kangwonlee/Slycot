#line 1 "SG03AX.f"
/* SG03AX.f -- translated by f2c (version 20100827).
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

#line 1 "SG03AX.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b11 = 1.;
static doublereal c_b12 = 0.;
static integer c__2 = 2;
static doublereal c_b16 = -1.;
static integer c__4 = 4;

/* Subroutine */ int sg03ax_(char *trans, integer *n, doublereal *a, integer *
	lda, doublereal *e, integer *lde, doublereal *x, integer *ldx, 
	doublereal *scale, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, e_dim1, e_offset, x_dim1, x_offset, i__1, i__2;

    /* Local variables */
    static integer i__, kb, lb, kh, lh, kl, ll;
    static doublereal tm[4]	/* was [2][2] */, ak11, ak12, ak21, ak22, 
	    al11, al12, al21, al22, ek11, ek12, ek22, el11, el12, el22, mat[
	    16]	/* was [4][4] */, rhs[4];
    static integer piv1[4], piv2[4], info1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), mb02uu_(integer *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *), mb02uv_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *), dcopy_(integer *, doublereal *, 
	    integer *, doublereal *, integer *), daxpy_(integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, integer *);
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

/*     To solve for X either the reduced generalized discrete-time */
/*     Lyapunov equation */

/*         T            T */
/*        A  * X * A - E  * X * E  =  SCALE * Y                       (1) */

/*     or */

/*                 T            T */
/*        A * X * A  - E * X * E   =  SCALE * Y                       (2) */

/*     where the right hand side Y is symmetric. A, E, Y, and the */
/*     solution X are N-by-N matrices. The pencil A - lambda * E must be */
/*     in generalized Schur form (A upper quasitriangular, E upper */
/*     triangular). SCALE is an output scale factor, set to avoid */
/*     overflow in X. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANS   CHARACTER*1 */
/*             Specifies whether the transposed equation is to be solved */
/*             or not: */
/*             = 'N':  Solve equation (1); */
/*             = 'T':  Solve equation (2). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N upper Hessenberg part of this array */
/*             must contain the quasitriangular matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     E       (input) DOUBLE PRECISION array, dimension (LDE,N) */
/*             The leading N-by-N upper triangular part of this array */
/*             must contain the matrix E. */

/*     LDE     INTEGER */
/*             The leading dimension of the array E.  LDE >= MAX(1,N). */

/*     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the right hand side matrix Y of the equation. Only */
/*             the upper triangular part of this matrix need be given. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the solution matrix X of the equation. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X.  LDX >= MAX(1,N). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor set to avoid overflow in X. */
/*             (0 < SCALE <= 1) */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  equation is (almost) singular to working precision; */
/*                   perturbed values were used to solve the equation */
/*                   (but the matrices A and E are unchanged). */

/*     METHOD */

/*     The solution X of (1) or (2) is computed via block back */
/*     substitution or block forward substitution, respectively. (See */
/*     [1] and [2] for details.) */

/*     REFERENCES */

/*     [1] Bartels, R.H., Stewart, G.W. */
/*         Solution of the equation A X + X B = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     [2] Penzl, T. */
/*         Numerical solution of generalized Lyapunov equations. */
/*         Advances in Comp. Math., vol. 8, pp. 33-48, 1998. */

/*     NUMERICAL ASPECTS */

/*     8/3 * N**3 flops are required by the routine. Note that we count a */
/*     single floating point arithmetic operation as one flop. */

/*     The algorithm is backward stable if the eigenvalues of the pencil */
/*     A - lambda * E are real. Otherwise, linear systems of order at */
/*     most 4 are involved into the computation. These systems are solved */
/*     by Gauss elimination with complete pivoting. The loss of stability */
/*     of the Gauss elimination with complete pivoting is rarely */
/*     encountered in practice. */

/*     CONTRIBUTOR */

/*     T. Penzl, Technical University Chemnitz, Germany, Aug. 1998. */

/*     REVISIONS */

/*     Sep. 1998 (V. Sima). */
/*     Dec. 1998 (V. Sima). */

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
/*     .. Executable Statements .. */

/*     Decode input parameter. */

#line 165 "SG03AX.f"
    /* Parameter adjustments */
#line 165 "SG03AX.f"
    a_dim1 = *lda;
#line 165 "SG03AX.f"
    a_offset = 1 + a_dim1;
#line 165 "SG03AX.f"
    a -= a_offset;
#line 165 "SG03AX.f"
    e_dim1 = *lde;
#line 165 "SG03AX.f"
    e_offset = 1 + e_dim1;
#line 165 "SG03AX.f"
    e -= e_offset;
#line 165 "SG03AX.f"
    x_dim1 = *ldx;
#line 165 "SG03AX.f"
    x_offset = 1 + x_dim1;
#line 165 "SG03AX.f"
    x -= x_offset;
#line 165 "SG03AX.f"

#line 165 "SG03AX.f"
    /* Function Body */
#line 165 "SG03AX.f"
    notrns = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

#line 169 "SG03AX.f"
    if (! (notrns || lsame_(trans, "T", (ftnlen)1, (ftnlen)1))) {
#line 170 "SG03AX.f"
	*info = -1;
#line 171 "SG03AX.f"
    } else if (*n < 0) {
#line 172 "SG03AX.f"
	*info = -2;
#line 173 "SG03AX.f"
    } else if (*lda < max(1,*n)) {
#line 174 "SG03AX.f"
	*info = -4;
#line 175 "SG03AX.f"
    } else if (*lde < max(1,*n)) {
#line 176 "SG03AX.f"
	*info = -6;
#line 177 "SG03AX.f"
    } else if (*ldx < max(1,*n)) {
#line 178 "SG03AX.f"
	*info = -8;
#line 179 "SG03AX.f"
    } else {
#line 180 "SG03AX.f"
	*info = 0;
#line 181 "SG03AX.f"
    }
#line 182 "SG03AX.f"
    if (*info != 0) {
#line 183 "SG03AX.f"
	i__1 = -(*info);
#line 183 "SG03AX.f"
	xerbla_("SG03AX", &i__1, (ftnlen)6);
#line 184 "SG03AX.f"
	return 0;
#line 185 "SG03AX.f"
    }

#line 187 "SG03AX.f"
    *scale = 1.;

/*     Quick return if possible. */

#line 191 "SG03AX.f"
    if (*n == 0) {
#line 191 "SG03AX.f"
	return 0;
#line 191 "SG03AX.f"
    }

#line 193 "SG03AX.f"
    if (notrns) {

/*        Solve equation (1). */

/*        Outer Loop. Compute block row X(KL:KH,:). KB denotes the number */
/*        of rows in this block row. */

#line 200 "SG03AX.f"
	kl = 0;
#line 201 "SG03AX.f"
	kb = 1;
/*        WHILE ( KL+KB .LE. N ) DO */
#line 203 "SG03AX.f"
L20:
#line 203 "SG03AX.f"
	if (kl + kb <= *n) {
#line 204 "SG03AX.f"
	    kl += kb;
#line 205 "SG03AX.f"
	    if (kl == *n) {
#line 206 "SG03AX.f"
		kb = 1;
#line 207 "SG03AX.f"
	    } else {
#line 208 "SG03AX.f"
		if (a[kl + 1 + kl * a_dim1] != 0.) {
#line 209 "SG03AX.f"
		    kb = 2;
#line 210 "SG03AX.f"
		} else {
#line 211 "SG03AX.f"
		    kb = 1;
#line 212 "SG03AX.f"
		}
#line 213 "SG03AX.f"
	    }
#line 214 "SG03AX.f"
	    kh = kl + kb - 1;

/*           Copy elements of solution already known by symmetry. */

/*              X(KL:KH,1:KL-1) = X(1:KL-1,KL:KH)' */

#line 220 "SG03AX.f"
	    if (kl > 1) {
#line 221 "SG03AX.f"
		i__1 = kh;
#line 221 "SG03AX.f"
		for (i__ = kl; i__ <= i__1; ++i__) {
#line 222 "SG03AX.f"
		    i__2 = kl - 1;
#line 222 "SG03AX.f"
		    dcopy_(&i__2, &x[i__ * x_dim1 + 1], &c__1, &x[i__ + 
			    x_dim1], ldx);
#line 223 "SG03AX.f"
/* L40: */
#line 223 "SG03AX.f"
		}
#line 224 "SG03AX.f"
	    }

/*           Inner Loop. Compute block X(KL:KH,LL:LH). LB denotes the */
/*           number of columns in this block. */

#line 229 "SG03AX.f"
	    ll = kl - 1;
#line 230 "SG03AX.f"
	    lb = 1;
/*           WHILE ( LL+LB .LE. N ) DO */
#line 232 "SG03AX.f"
L60:
#line 232 "SG03AX.f"
	    if (ll + lb <= *n) {
#line 233 "SG03AX.f"
		ll += lb;
#line 234 "SG03AX.f"
		if (ll == *n) {
#line 235 "SG03AX.f"
		    lb = 1;
#line 236 "SG03AX.f"
		} else {
#line 237 "SG03AX.f"
		    if (a[ll + 1 + ll * a_dim1] != 0.) {
#line 238 "SG03AX.f"
			lb = 2;
#line 239 "SG03AX.f"
		    } else {
#line 240 "SG03AX.f"
			lb = 1;
#line 241 "SG03AX.f"
		    }
#line 242 "SG03AX.f"
		}
#line 243 "SG03AX.f"
		lh = ll + lb - 1;

/*              Update right hand sides (I). */

/*                 X(KL:LH,LL:LH) = X(KL:LH,LL:LH) - */
/*                    A(KL:KH,KL:LH)'*(X(KL:KH,1:LL-1)*A(1:LL-1,LL:LH)) */

/*                 X(KL:LH,LL:LH) = X(KL:LH,LL:LH) + */
/*                    E(KL:KH,KL:LH)'*(X(KL:KH,1:LL-1)*E(1:LL-1,LL:LH)) */

#line 253 "SG03AX.f"
		if (ll > 1) {
#line 254 "SG03AX.f"
		    i__1 = ll - 1;
#line 254 "SG03AX.f"
		    dgemm_("N", "N", &kb, &lb, &i__1, &c_b11, &x[kl + x_dim1],
			     ldx, &a[ll * a_dim1 + 1], lda, &c_b12, tm, &c__2,
			     (ftnlen)1, (ftnlen)1);
#line 256 "SG03AX.f"
		    i__1 = lh - kl + 1;
#line 256 "SG03AX.f"
		    dgemm_("T", "N", &i__1, &lb, &kb, &c_b16, &a[kl + kl * 
			    a_dim1], lda, tm, &c__2, &c_b11, &x[kl + ll * 
			    x_dim1], ldx, (ftnlen)1, (ftnlen)1);
#line 258 "SG03AX.f"
		    i__1 = ll - 1;
#line 258 "SG03AX.f"
		    dgemm_("N", "N", &kb, &lb, &i__1, &c_b11, &x[kl + x_dim1],
			     ldx, &e[ll * e_dim1 + 1], lde, &c_b12, tm, &c__2,
			     (ftnlen)1, (ftnlen)1);
#line 260 "SG03AX.f"
		    i__1 = lh - kh + 1;
#line 260 "SG03AX.f"
		    dgemm_("T", "N", &i__1, &lb, &kb, &c_b11, &e[kl + kh * 
			    e_dim1], lde, tm, &c__2, &c_b11, &x[kh + ll * 
			    x_dim1], ldx, (ftnlen)1, (ftnlen)1);
#line 262 "SG03AX.f"
		    if (kb == 2) {
#line 262 "SG03AX.f"
			daxpy_(&lb, &e[kl + kl * e_dim1], tm, &c__2, &x[kl + 
				ll * x_dim1], ldx);
#line 262 "SG03AX.f"
		    }
#line 264 "SG03AX.f"
		}

/*              Solve small Sylvester equations of order at most (2,2). */

#line 268 "SG03AX.f"
		if (kb == 1 && lb == 1) {

#line 270 "SG03AX.f"
		    dimmat = 1;

#line 272 "SG03AX.f"
		    mat[0] = a[ll + ll * a_dim1] * a[kl + kl * a_dim1] - e[ll 
			    + ll * e_dim1] * e[kl + kl * e_dim1];

#line 274 "SG03AX.f"
		    rhs[0] = x[kl + ll * x_dim1];

#line 276 "SG03AX.f"
		} else if (kb == 2 && lb == 1) {

#line 278 "SG03AX.f"
		    dimmat = 2;

#line 280 "SG03AX.f"
		    ak11 = a[kl + kl * a_dim1];
#line 281 "SG03AX.f"
		    ak12 = a[kl + kh * a_dim1];
#line 282 "SG03AX.f"
		    ak21 = a[kh + kl * a_dim1];
#line 283 "SG03AX.f"
		    ak22 = a[kh + kh * a_dim1];

#line 285 "SG03AX.f"
		    al11 = a[ll + ll * a_dim1];

#line 287 "SG03AX.f"
		    ek11 = e[kl + kl * e_dim1];
#line 288 "SG03AX.f"
		    ek12 = e[kl + kh * e_dim1];
#line 289 "SG03AX.f"
		    ek22 = e[kh + kh * e_dim1];

#line 291 "SG03AX.f"
		    el11 = e[ll + ll * e_dim1];

#line 293 "SG03AX.f"
		    mat[0] = al11 * ak11 - el11 * ek11;
#line 294 "SG03AX.f"
		    mat[4] = al11 * ak21;
#line 295 "SG03AX.f"
		    mat[1] = al11 * ak12 - el11 * ek12;
#line 296 "SG03AX.f"
		    mat[5] = al11 * ak22 - el11 * ek22;

#line 298 "SG03AX.f"
		    rhs[0] = x[kl + ll * x_dim1];
#line 299 "SG03AX.f"
		    rhs[1] = x[kh + ll * x_dim1];

#line 301 "SG03AX.f"
		} else if (kb == 1 && lb == 2) {

#line 303 "SG03AX.f"
		    dimmat = 2;

#line 305 "SG03AX.f"
		    ak11 = a[kl + kl * a_dim1];

#line 307 "SG03AX.f"
		    al11 = a[ll + ll * a_dim1];
#line 308 "SG03AX.f"
		    al12 = a[ll + lh * a_dim1];
#line 309 "SG03AX.f"
		    al21 = a[lh + ll * a_dim1];
#line 310 "SG03AX.f"
		    al22 = a[lh + lh * a_dim1];

#line 312 "SG03AX.f"
		    ek11 = e[kl + kl * e_dim1];

#line 314 "SG03AX.f"
		    el11 = e[ll + ll * e_dim1];
#line 315 "SG03AX.f"
		    el12 = e[ll + lh * e_dim1];
#line 316 "SG03AX.f"
		    el22 = e[lh + lh * e_dim1];

#line 318 "SG03AX.f"
		    mat[0] = al11 * ak11 - el11 * ek11;
#line 319 "SG03AX.f"
		    mat[4] = al21 * ak11;
#line 320 "SG03AX.f"
		    mat[1] = al12 * ak11 - el12 * ek11;
#line 321 "SG03AX.f"
		    mat[5] = al22 * ak11 - el22 * ek11;

#line 323 "SG03AX.f"
		    rhs[0] = x[kl + ll * x_dim1];
#line 324 "SG03AX.f"
		    rhs[1] = x[kl + lh * x_dim1];

#line 326 "SG03AX.f"
		} else {

#line 328 "SG03AX.f"
		    dimmat = 4;

#line 330 "SG03AX.f"
		    ak11 = a[kl + kl * a_dim1];
#line 331 "SG03AX.f"
		    ak12 = a[kl + kh * a_dim1];
#line 332 "SG03AX.f"
		    ak21 = a[kh + kl * a_dim1];
#line 333 "SG03AX.f"
		    ak22 = a[kh + kh * a_dim1];

#line 335 "SG03AX.f"
		    al11 = a[ll + ll * a_dim1];
#line 336 "SG03AX.f"
		    al12 = a[ll + lh * a_dim1];
#line 337 "SG03AX.f"
		    al21 = a[lh + ll * a_dim1];
#line 338 "SG03AX.f"
		    al22 = a[lh + lh * a_dim1];

#line 340 "SG03AX.f"
		    ek11 = e[kl + kl * e_dim1];
#line 341 "SG03AX.f"
		    ek12 = e[kl + kh * e_dim1];
#line 342 "SG03AX.f"
		    ek22 = e[kh + kh * e_dim1];

#line 344 "SG03AX.f"
		    el11 = e[ll + ll * e_dim1];
#line 345 "SG03AX.f"
		    el12 = e[ll + lh * e_dim1];
#line 346 "SG03AX.f"
		    el22 = e[lh + lh * e_dim1];

#line 348 "SG03AX.f"
		    mat[0] = al11 * ak11 - el11 * ek11;
#line 349 "SG03AX.f"
		    mat[4] = al11 * ak21;
#line 350 "SG03AX.f"
		    mat[8] = al21 * ak11;
#line 351 "SG03AX.f"
		    mat[12] = al21 * ak21;

#line 353 "SG03AX.f"
		    mat[1] = al11 * ak12 - el11 * ek12;
#line 354 "SG03AX.f"
		    mat[5] = al11 * ak22 - el11 * ek22;
#line 355 "SG03AX.f"
		    mat[9] = al21 * ak12;
#line 356 "SG03AX.f"
		    mat[13] = al21 * ak22;

#line 358 "SG03AX.f"
		    mat[2] = al12 * ak11 - el12 * ek11;
#line 359 "SG03AX.f"
		    mat[6] = al12 * ak21;
#line 360 "SG03AX.f"
		    mat[10] = al22 * ak11 - el22 * ek11;
#line 361 "SG03AX.f"
		    mat[14] = al22 * ak21;

#line 363 "SG03AX.f"
		    mat[3] = al12 * ak12 - el12 * ek12;
#line 364 "SG03AX.f"
		    mat[7] = al12 * ak22 - el12 * ek22;
#line 365 "SG03AX.f"
		    mat[11] = al22 * ak12 - el22 * ek12;
#line 366 "SG03AX.f"
		    mat[15] = al22 * ak22 - el22 * ek22;

#line 368 "SG03AX.f"
		    rhs[0] = x[kl + ll * x_dim1];
#line 369 "SG03AX.f"
		    if (kl == ll) {
#line 370 "SG03AX.f"
			rhs[1] = x[kl + kh * x_dim1];
#line 371 "SG03AX.f"
		    } else {
#line 372 "SG03AX.f"
			rhs[1] = x[kh + ll * x_dim1];
#line 373 "SG03AX.f"
		    }
#line 374 "SG03AX.f"
		    rhs[2] = x[kl + lh * x_dim1];
#line 375 "SG03AX.f"
		    rhs[3] = x[kh + lh * x_dim1];

#line 377 "SG03AX.f"
		}

#line 379 "SG03AX.f"
		mb02uv_(&dimmat, mat, &c__4, piv1, piv2, &info1);
#line 380 "SG03AX.f"
		if (info1 != 0) {
#line 380 "SG03AX.f"
		    *info = 1;
#line 380 "SG03AX.f"
		}
#line 382 "SG03AX.f"
		mb02uu_(&dimmat, mat, &c__4, rhs, piv1, piv2, &scale1);

/*              Scaling. */

#line 386 "SG03AX.f"
		if (scale1 != 1.) {
#line 387 "SG03AX.f"
		    i__1 = *n;
#line 387 "SG03AX.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 388 "SG03AX.f"
			dscal_(n, &scale1, &x[i__ * x_dim1 + 1], &c__1);
#line 389 "SG03AX.f"
/* L80: */
#line 389 "SG03AX.f"
		    }
#line 390 "SG03AX.f"
		    *scale *= scale1;
#line 391 "SG03AX.f"
		}

#line 393 "SG03AX.f"
		if (lb == 1 && kb == 1) {
#line 394 "SG03AX.f"
		    x[kl + ll * x_dim1] = rhs[0];
#line 395 "SG03AX.f"
		} else if (lb == 1 && kb == 2) {
#line 396 "SG03AX.f"
		    x[kl + ll * x_dim1] = rhs[0];
#line 397 "SG03AX.f"
		    x[kh + ll * x_dim1] = rhs[1];
#line 398 "SG03AX.f"
		} else if (lb == 2 && kb == 1) {
#line 399 "SG03AX.f"
		    x[kl + ll * x_dim1] = rhs[0];
#line 400 "SG03AX.f"
		    x[kl + lh * x_dim1] = rhs[1];
#line 401 "SG03AX.f"
		} else {
#line 402 "SG03AX.f"
		    x[kl + ll * x_dim1] = rhs[0];
#line 403 "SG03AX.f"
		    x[kh + ll * x_dim1] = rhs[1];
#line 404 "SG03AX.f"
		    x[kl + lh * x_dim1] = rhs[2];
#line 405 "SG03AX.f"
		    x[kh + lh * x_dim1] = rhs[3];
#line 406 "SG03AX.f"
		}

/*              Update right hand sides (II). */

/*              X(KH+1:LH,LL:LH) = X(KH+1:LH,LL:LH) - */
/*                 A(KL:KH,KH+1:LH)'*(X(KL:KH,LL:LH)*A(LL:LH,LL:LH)) */

/*              X(KH+1:LH,LL:LH) = X(KH+1:LH,LL:LH) + */
/*                 E(KL:KH,KH+1:LH)'*(X(KL:KH,LL:LH)*E(LL:LH,LL:LH)) */

#line 416 "SG03AX.f"
		if (kl < ll) {
#line 417 "SG03AX.f"
		    dgemm_("N", "N", &kb, &lb, &lb, &c_b11, &x[kl + ll * 
			    x_dim1], ldx, &a[ll + ll * a_dim1], lda, &c_b12, 
			    tm, &c__2, (ftnlen)1, (ftnlen)1);
#line 419 "SG03AX.f"
		    i__1 = lh - kh;
#line 419 "SG03AX.f"
		    dgemm_("T", "N", &i__1, &lb, &kb, &c_b16, &a[kl + (kh + 1)
			     * a_dim1], lda, tm, &c__2, &c_b11, &x[kh + 1 + 
			    ll * x_dim1], ldx, (ftnlen)1, (ftnlen)1);
#line 421 "SG03AX.f"
		    if (lb == 2) {
#line 422 "SG03AX.f"
			dcopy_(&kb, &x[kl + ll * x_dim1], &c__1, tm, &c__1);
#line 423 "SG03AX.f"
			dscal_(&kb, &e[ll + ll * e_dim1], tm, &c__1);
#line 424 "SG03AX.f"
		    }
#line 425 "SG03AX.f"
		    dgemv_("N", &kb, &lb, &c_b11, &x[kl + ll * x_dim1], ldx, &
			    e[ll + lh * e_dim1], &c__1, &c_b12, &tm[(lb << 1) 
			    - 2], &c__1, (ftnlen)1);
#line 427 "SG03AX.f"
		    i__1 = lh - kh;
#line 427 "SG03AX.f"
		    dgemm_("T", "N", &i__1, &lb, &kb, &c_b11, &e[kl + (kh + 1)
			     * e_dim1], lde, tm, &c__2, &c_b11, &x[kh + 1 + 
			    ll * x_dim1], ldx, (ftnlen)1, (ftnlen)1);
#line 429 "SG03AX.f"
		}

#line 431 "SG03AX.f"
		goto L60;
#line 432 "SG03AX.f"
	    }
/*           END WHILE 60 */

#line 435 "SG03AX.f"
	    goto L20;
#line 436 "SG03AX.f"
	}
/*        END WHILE 20 */

#line 439 "SG03AX.f"
    } else {

/*        Solve equation (2). */

/*        Outer Loop. Compute block column X(:,LL:LH). LB denotes the */
/*        number of columns in this block column. */

#line 446 "SG03AX.f"
	ll = *n + 1;
/*        WHILE ( LL .GT. 1 ) DO */
#line 448 "SG03AX.f"
L100:
#line 448 "SG03AX.f"
	if (ll > 1) {
#line 449 "SG03AX.f"
	    lh = ll - 1;
#line 450 "SG03AX.f"
	    if (lh == 1) {
#line 451 "SG03AX.f"
		lb = 1;
#line 452 "SG03AX.f"
	    } else {
#line 453 "SG03AX.f"
		if (a[ll - 1 + (ll - 2) * a_dim1] != 0.) {
#line 454 "SG03AX.f"
		    lb = 2;
#line 455 "SG03AX.f"
		} else {
#line 456 "SG03AX.f"
		    lb = 1;
#line 457 "SG03AX.f"
		}
#line 458 "SG03AX.f"
	    }
#line 459 "SG03AX.f"
	    ll -= lb;

/*           Copy elements of solution already known by symmetry. */

/*              X(LH+1:N,LL:LH) = X(LL:LH,LH+1:N)' */

#line 465 "SG03AX.f"
	    if (lh < *n) {
#line 466 "SG03AX.f"
		i__1 = lh;
#line 466 "SG03AX.f"
		for (i__ = ll; i__ <= i__1; ++i__) {
#line 467 "SG03AX.f"
		    i__2 = *n - lh;
#line 467 "SG03AX.f"
		    dcopy_(&i__2, &x[i__ + (lh + 1) * x_dim1], ldx, &x[lh + 1 
			    + i__ * x_dim1], &c__1);
#line 468 "SG03AX.f"
/* L120: */
#line 468 "SG03AX.f"
		}
#line 469 "SG03AX.f"
	    }

/*           Inner Loop. Compute block X(KL:KH,LL:LH). KB denotes the */
/*           number of rows in this block. */

#line 474 "SG03AX.f"
	    kl = lh + 1;
/*           WHILE ( KL .GT. 1 ) DO */
#line 476 "SG03AX.f"
L140:
#line 476 "SG03AX.f"
	    if (kl > 1) {
#line 477 "SG03AX.f"
		kh = kl - 1;
#line 478 "SG03AX.f"
		if (kh == 1) {
#line 479 "SG03AX.f"
		    kb = 1;
#line 480 "SG03AX.f"
		} else {
#line 481 "SG03AX.f"
		    if (a[kl - 1 + (kl - 2) * a_dim1] != 0.) {
#line 482 "SG03AX.f"
			kb = 2;
#line 483 "SG03AX.f"
		    } else {
#line 484 "SG03AX.f"
			kb = 1;
#line 485 "SG03AX.f"
		    }
#line 486 "SG03AX.f"
		}
#line 487 "SG03AX.f"
		kl -= kb;

/*              Update right hand sides (I). */

/*                 X(KL:KH,KL:LH) = X(KL:KH,KL:LH) - */
/*                    (A(KL:KH,KH+1:N)*X(KH+1:N,LL:LH))*A(KL:LH,LL:LH)' */

/*                 X(KL:KH,KL:LH) = X(KL:KH,KL:LH) + */
/*                    (E(KL:KH,KH+1:N)*X(KH+1:N,LL:LH))*E(KL:LH,LL:LH)' */

#line 497 "SG03AX.f"
		if (kh < *n) {
#line 498 "SG03AX.f"
		    i__1 = *n - kh;
#line 498 "SG03AX.f"
		    dgemm_("N", "N", &kb, &lb, &i__1, &c_b11, &a[kl + (kh + 1)
			     * a_dim1], lda, &x[kh + 1 + ll * x_dim1], ldx, &
			    c_b12, tm, &c__2, (ftnlen)1, (ftnlen)1);
#line 500 "SG03AX.f"
		    i__1 = lh - kl + 1;
#line 500 "SG03AX.f"
		    dgemm_("N", "T", &kb, &i__1, &lb, &c_b16, tm, &c__2, &a[
			    kl + ll * a_dim1], lda, &c_b11, &x[kl + kl * 
			    x_dim1], ldx, (ftnlen)1, (ftnlen)1);
#line 502 "SG03AX.f"
		    i__1 = *n - kh;
#line 502 "SG03AX.f"
		    dgemm_("N", "N", &kb, &lb, &i__1, &c_b11, &e[kl + (kh + 1)
			     * e_dim1], lde, &x[kh + 1 + ll * x_dim1], ldx, &
			    c_b12, tm, &c__2, (ftnlen)1, (ftnlen)1);
#line 504 "SG03AX.f"
		    i__1 = ll - kl + 1;
#line 504 "SG03AX.f"
		    dgemm_("N", "T", &kb, &i__1, &lb, &c_b11, tm, &c__2, &e[
			    kl + ll * e_dim1], lde, &c_b11, &x[kl + kl * 
			    x_dim1], ldx, (ftnlen)1, (ftnlen)1);
#line 506 "SG03AX.f"
		    if (lb == 2) {
#line 506 "SG03AX.f"
			daxpy_(&kb, &e[lh + lh * e_dim1], &tm[2], &c__1, &x[
				kl + lh * x_dim1], &c__1);
#line 506 "SG03AX.f"
		    }
#line 508 "SG03AX.f"
		}

/*              Solve small Sylvester equations of order at most (2,2). */

#line 512 "SG03AX.f"
		if (kb == 1 && lb == 1) {

#line 514 "SG03AX.f"
		    dimmat = 1;

#line 516 "SG03AX.f"
		    mat[0] = a[ll + ll * a_dim1] * a[kl + kl * a_dim1] - e[ll 
			    + ll * e_dim1] * e[kl + kl * e_dim1];

#line 518 "SG03AX.f"
		    rhs[0] = x[kl + ll * x_dim1];

#line 520 "SG03AX.f"
		} else if (kb == 2 && lb == 1) {

#line 522 "SG03AX.f"
		    dimmat = 2;

#line 524 "SG03AX.f"
		    ak11 = a[kl + kl * a_dim1];
#line 525 "SG03AX.f"
		    ak12 = a[kl + kh * a_dim1];
#line 526 "SG03AX.f"
		    ak21 = a[kh + kl * a_dim1];
#line 527 "SG03AX.f"
		    ak22 = a[kh + kh * a_dim1];

#line 529 "SG03AX.f"
		    al11 = a[ll + ll * a_dim1];

#line 531 "SG03AX.f"
		    ek11 = e[kl + kl * e_dim1];
#line 532 "SG03AX.f"
		    ek12 = e[kl + kh * e_dim1];
#line 533 "SG03AX.f"
		    ek22 = e[kh + kh * e_dim1];

#line 535 "SG03AX.f"
		    el11 = e[ll + ll * e_dim1];

#line 537 "SG03AX.f"
		    mat[0] = al11 * ak11 - el11 * ek11;
#line 538 "SG03AX.f"
		    mat[4] = al11 * ak12 - el11 * ek12;
#line 539 "SG03AX.f"
		    mat[1] = al11 * ak21;
#line 540 "SG03AX.f"
		    mat[5] = al11 * ak22 - el11 * ek22;

#line 542 "SG03AX.f"
		    rhs[0] = x[kl + ll * x_dim1];
#line 543 "SG03AX.f"
		    rhs[1] = x[kh + ll * x_dim1];

#line 545 "SG03AX.f"
		} else if (kb == 1 && lb == 2) {

#line 547 "SG03AX.f"
		    dimmat = 2;

#line 549 "SG03AX.f"
		    ak11 = a[kl + kl * a_dim1];

#line 551 "SG03AX.f"
		    al11 = a[ll + ll * a_dim1];
#line 552 "SG03AX.f"
		    al12 = a[ll + lh * a_dim1];
#line 553 "SG03AX.f"
		    al21 = a[lh + ll * a_dim1];
#line 554 "SG03AX.f"
		    al22 = a[lh + lh * a_dim1];

#line 556 "SG03AX.f"
		    ek11 = e[kl + kl * e_dim1];

#line 558 "SG03AX.f"
		    el11 = e[ll + ll * e_dim1];
#line 559 "SG03AX.f"
		    el12 = e[ll + lh * e_dim1];
#line 560 "SG03AX.f"
		    el22 = e[lh + lh * e_dim1];

#line 562 "SG03AX.f"
		    mat[0] = al11 * ak11 - el11 * ek11;
#line 563 "SG03AX.f"
		    mat[4] = al12 * ak11 - el12 * ek11;
#line 564 "SG03AX.f"
		    mat[1] = al21 * ak11;
#line 565 "SG03AX.f"
		    mat[5] = al22 * ak11 - el22 * ek11;

#line 567 "SG03AX.f"
		    rhs[0] = x[kl + ll * x_dim1];
#line 568 "SG03AX.f"
		    rhs[1] = x[kl + lh * x_dim1];

#line 570 "SG03AX.f"
		} else {

#line 572 "SG03AX.f"
		    dimmat = 4;

#line 574 "SG03AX.f"
		    ak11 = a[kl + kl * a_dim1];
#line 575 "SG03AX.f"
		    ak12 = a[kl + kh * a_dim1];
#line 576 "SG03AX.f"
		    ak21 = a[kh + kl * a_dim1];
#line 577 "SG03AX.f"
		    ak22 = a[kh + kh * a_dim1];

#line 579 "SG03AX.f"
		    al11 = a[ll + ll * a_dim1];
#line 580 "SG03AX.f"
		    al12 = a[ll + lh * a_dim1];
#line 581 "SG03AX.f"
		    al21 = a[lh + ll * a_dim1];
#line 582 "SG03AX.f"
		    al22 = a[lh + lh * a_dim1];

#line 584 "SG03AX.f"
		    ek11 = e[kl + kl * e_dim1];
#line 585 "SG03AX.f"
		    ek12 = e[kl + kh * e_dim1];
#line 586 "SG03AX.f"
		    ek22 = e[kh + kh * e_dim1];

#line 588 "SG03AX.f"
		    el11 = e[ll + ll * e_dim1];
#line 589 "SG03AX.f"
		    el12 = e[ll + lh * e_dim1];
#line 590 "SG03AX.f"
		    el22 = e[lh + lh * e_dim1];

#line 592 "SG03AX.f"
		    mat[0] = al11 * ak11 - el11 * ek11;
#line 593 "SG03AX.f"
		    mat[4] = al11 * ak12 - el11 * ek12;
#line 594 "SG03AX.f"
		    mat[8] = al12 * ak11 - el12 * ek11;
#line 595 "SG03AX.f"
		    mat[12] = al12 * ak12 - el12 * ek12;

#line 597 "SG03AX.f"
		    mat[1] = al11 * ak21;
#line 598 "SG03AX.f"
		    mat[5] = al11 * ak22 - el11 * ek22;
#line 599 "SG03AX.f"
		    mat[9] = al12 * ak21;
#line 600 "SG03AX.f"
		    mat[13] = al12 * ak22 - el12 * ek22;

#line 602 "SG03AX.f"
		    mat[2] = al21 * ak11;
#line 603 "SG03AX.f"
		    mat[6] = al21 * ak12;
#line 604 "SG03AX.f"
		    mat[10] = al22 * ak11 - el22 * ek11;
#line 605 "SG03AX.f"
		    mat[14] = al22 * ak12 - el22 * ek12;

#line 607 "SG03AX.f"
		    mat[3] = al21 * ak21;
#line 608 "SG03AX.f"
		    mat[7] = al21 * ak22;
#line 609 "SG03AX.f"
		    mat[11] = al22 * ak21;
#line 610 "SG03AX.f"
		    mat[15] = al22 * ak22 - el22 * ek22;

#line 612 "SG03AX.f"
		    rhs[0] = x[kl + ll * x_dim1];
#line 613 "SG03AX.f"
		    if (kl == ll) {
#line 614 "SG03AX.f"
			rhs[1] = x[kl + kh * x_dim1];
#line 615 "SG03AX.f"
		    } else {
#line 616 "SG03AX.f"
			rhs[1] = x[kh + ll * x_dim1];
#line 617 "SG03AX.f"
		    }
#line 618 "SG03AX.f"
		    rhs[2] = x[kl + lh * x_dim1];
#line 619 "SG03AX.f"
		    rhs[3] = x[kh + lh * x_dim1];

#line 621 "SG03AX.f"
		}

#line 623 "SG03AX.f"
		mb02uv_(&dimmat, mat, &c__4, piv1, piv2, &info1);
#line 624 "SG03AX.f"
		if (info1 != 0) {
#line 624 "SG03AX.f"
		    *info = 1;
#line 624 "SG03AX.f"
		}
#line 626 "SG03AX.f"
		mb02uu_(&dimmat, mat, &c__4, rhs, piv1, piv2, &scale1);

/*              Scaling. */

#line 630 "SG03AX.f"
		if (scale1 != 1.) {
#line 631 "SG03AX.f"
		    i__1 = *n;
#line 631 "SG03AX.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 632 "SG03AX.f"
			dscal_(n, &scale1, &x[i__ * x_dim1 + 1], &c__1);
#line 633 "SG03AX.f"
/* L160: */
#line 633 "SG03AX.f"
		    }
#line 634 "SG03AX.f"
		    *scale *= scale1;
#line 635 "SG03AX.f"
		}

#line 637 "SG03AX.f"
		if (lb == 1 && kb == 1) {
#line 638 "SG03AX.f"
		    x[kl + ll * x_dim1] = rhs[0];
#line 639 "SG03AX.f"
		} else if (lb == 1 && kb == 2) {
#line 640 "SG03AX.f"
		    x[kl + ll * x_dim1] = rhs[0];
#line 641 "SG03AX.f"
		    x[kh + ll * x_dim1] = rhs[1];
#line 642 "SG03AX.f"
		} else if (lb == 2 && kb == 1) {
#line 643 "SG03AX.f"
		    x[kl + ll * x_dim1] = rhs[0];
#line 644 "SG03AX.f"
		    x[kl + lh * x_dim1] = rhs[1];
#line 645 "SG03AX.f"
		} else {
#line 646 "SG03AX.f"
		    x[kl + ll * x_dim1] = rhs[0];
#line 647 "SG03AX.f"
		    x[kh + ll * x_dim1] = rhs[1];
#line 648 "SG03AX.f"
		    x[kl + lh * x_dim1] = rhs[2];
#line 649 "SG03AX.f"
		    x[kh + lh * x_dim1] = rhs[3];
#line 650 "SG03AX.f"
		}

/*              Update right hand sides (II). */

/*                 X(KL:KH,KL:LL-1) = X(KL:KH,KL:LL-1) - */
/*                    (A(KL:KH,KL:KH)*X(KL:KH,LL:LH))*A(KL:LL-1,LL:LH)' */

/*                 X(KL:KH,KL:LL-1) = X(KL:KH,KL:LL-1) + */
/*                    (E(KL:KH,KL:KH)*X(KL:KH,LL:LH))*E(KL:LL-1,LL:LH)' */

#line 660 "SG03AX.f"
		if (kl < ll) {
#line 661 "SG03AX.f"
		    dgemm_("N", "N", &kb, &lb, &kb, &c_b11, &a[kl + kl * 
			    a_dim1], lda, &x[kl + ll * x_dim1], ldx, &c_b12, 
			    tm, &c__2, (ftnlen)1, (ftnlen)1);
#line 663 "SG03AX.f"
		    i__1 = ll - kl;
#line 663 "SG03AX.f"
		    dgemm_("N", "T", &kb, &i__1, &lb, &c_b16, tm, &c__2, &a[
			    kl + ll * a_dim1], lda, &c_b11, &x[kl + kl * 
			    x_dim1], ldx, (ftnlen)1, (ftnlen)1);
#line 665 "SG03AX.f"
		    dgemv_("T", &kb, &lb, &c_b11, &x[kl + ll * x_dim1], ldx, &
			    e[kl + kl * e_dim1], lde, &c_b12, tm, &c__2, (
			    ftnlen)1);
#line 667 "SG03AX.f"
		    if (kb == 2) {
#line 668 "SG03AX.f"
			dcopy_(&lb, &x[kh + ll * x_dim1], ldx, &tm[1], &c__2);
#line 669 "SG03AX.f"
			dscal_(&lb, &e[kh + kh * e_dim1], &tm[1], &c__2);
#line 670 "SG03AX.f"
		    }
#line 671 "SG03AX.f"
		    i__1 = ll - kl;
#line 671 "SG03AX.f"
		    dgemm_("N", "T", &kb, &i__1, &lb, &c_b11, tm, &c__2, &e[
			    kl + ll * e_dim1], lde, &c_b11, &x[kl + kl * 
			    x_dim1], ldx, (ftnlen)1, (ftnlen)1);
#line 673 "SG03AX.f"
		}

#line 675 "SG03AX.f"
		goto L140;
#line 676 "SG03AX.f"
	    }
/*           END WHILE 140 */

#line 679 "SG03AX.f"
	    goto L100;
#line 680 "SG03AX.f"
	}
/*        END WHILE 100 */

#line 683 "SG03AX.f"
    }

#line 685 "SG03AX.f"
    return 0;
/* *** Last line of SG03AX *** */
} /* sg03ax_ */

