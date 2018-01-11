#line 1 "SB04OW.f"
/* SB04OW.f -- translated by f2c (version 20100827).
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

#line 1 "SB04OW.f"
/* Table of constant values */

static integer c__8 = 8;
static integer c__1 = 1;
static doublereal c_b23 = -1.;
static doublereal c_b37 = 1.;
static doublereal c_b51 = 0.;

/* Subroutine */ int sb04ow_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *e, integer *lde, 
	doublereal *f, integer *ldf, doublereal *scale, integer *iwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, e_dim1, e_offset, f_dim1, f_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, p, q;
    static doublereal z__[64]	/* was [8][8] */;
    static integer ie, je, mb, nb, ii, jj, is, js;
    static doublereal rhs[8];
    static integer isp1, jsp1;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ierr, zdim, ipiv[8], jpiv[8];
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen), dgemv_(
	    char *, integer *, integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, integer *,
	     ftnlen), dcopy_(integer *, doublereal *, integer *, doublereal *,
	     integer *), daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dgesc2_(integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *), dgetc2_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *);
    static doublereal scaloc;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);


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

/*     To solve a periodic Sylvester equation */

/*              A * R - L * B = scale * C                           (1) */
/*              D * L - R * E = scale * F, */

/*     using Level 1 and 2 BLAS, where R and L are unknown M-by-N */
/*     matrices, (A, D), (B, E) and (C, F) are given matrix pairs of */
/*     size M-by-M, N-by-N and M-by-N, respectively, with real entries. */
/*     (A, D) and (B, E) must be in periodic Schur form, i.e. A, B are */
/*     upper quasi triangular and D, E are upper triangular. The solution */
/*     (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output scaling */
/*     factor chosen to avoid overflow. */

/*     This routine is largely based on the LAPACK routine DTGSY2 */
/*     developed by Bo Kagstrom and Peter Poromaa. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of A and D, and the row dimension of C, F, R */
/*             and L.  M >= 0. */

/*     N       (input) INTEGER */
/*             The order of B and E, and the column dimension of C, F, R */
/*             and L.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,M) */
/*             On entry, the leading M-by-M part of this array must */
/*             contain the upper quasi triangular matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,M). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper quasi triangular matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the right-hand-side of the first matrix equation */
/*             in (1). */
/*             On exit, the leading M-by-N part of this array contains */
/*             the solution R. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= MAX(1,M). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading M-by-M part of this array must */
/*             contain the upper triangular matrix D. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= MAX(1,M). */

/*     E       (input) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper triangular matrix E. */

/*     LDE     INTEGER */
/*             The leading dimension of the array E.  LDE >= MAX(1,N). */

/*     F       (input/output) DOUBLE PRECISION array, dimension (LDF,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the right-hand-side of the second matrix equation */
/*             in (1). */
/*             On exit, the leading M-by-N part of this array contains */
/*             the solution L. */

/*     LDF     INTEGER */
/*             The leading dimension of the array F.  LDF >= MAX(1,M). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             On exit, 0 <= SCALE <= 1. If 0 < SCALE < 1, the arrays */
/*             C and F will hold the solutions R and L, respectively, to */
/*             a slightly perturbed system but the input matrices A, B, D */
/*             and E have not been changed. If SCALE = 0, C and F will */
/*             hold solutions to the homogeneous system with C = F = 0. */
/*             Normally, SCALE = 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M+N+2) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  the matrix products A*D and B*E have common or very */
/*                   close eigenvalues. */

/*     METHOD */

/*     In matrix notation solving equation (1) corresponds to solving */
/*     Z*x = scale*b, where Z is defined as */

/*         Z = [  kron(In, A)  -kron(B', Im) ]            (2) */
/*             [ -kron(E', Im)  kron(In, D)  ], */

/*     Ik is the identity matrix of size k and X' is the transpose of X. */
/*     kron(X, Y) is the Kronecker product between the matrices X and Y. */
/*     In the process of solving (1), we solve a number of such systems */
/*     where Dim(Im), Dim(In) = 1 or 2. */

/*     REFERENCES */

/*     [1] Kagstrom, B. */
/*         A Direct Method for Reordering Eigenvalues in the Generalized */
/*         Real Schur Form of a Regular Matrix Pair (A,B). M.S. Moonen */
/*         et al (eds.), Linear Algebra for Large Scale and Real-Time */
/*         Applications, Kluwer Academic Publ., pp. 195-218, 1993. */

/*     [2] Sreedhar, J. and Van Dooren, P. */
/*         A Schur approach for solving some periodic matrix equations. */
/*         U. Helmke et al (eds.), Systems and Networks: Mathematical */
/*         Theory and Applications, Akademie Verlag, Berlin, vol. 77, */
/*         pp. 339-362, 1994. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DTGPY2). */

/*     KEYWORDS */

/*     Matrix equation, periodic Sylvester equation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

#line 192 "SB04OW.f"
    /* Parameter adjustments */
#line 192 "SB04OW.f"
    a_dim1 = *lda;
#line 192 "SB04OW.f"
    a_offset = 1 + a_dim1;
#line 192 "SB04OW.f"
    a -= a_offset;
#line 192 "SB04OW.f"
    b_dim1 = *ldb;
#line 192 "SB04OW.f"
    b_offset = 1 + b_dim1;
#line 192 "SB04OW.f"
    b -= b_offset;
#line 192 "SB04OW.f"
    c_dim1 = *ldc;
#line 192 "SB04OW.f"
    c_offset = 1 + c_dim1;
#line 192 "SB04OW.f"
    c__ -= c_offset;
#line 192 "SB04OW.f"
    d_dim1 = *ldd;
#line 192 "SB04OW.f"
    d_offset = 1 + d_dim1;
#line 192 "SB04OW.f"
    d__ -= d_offset;
#line 192 "SB04OW.f"
    e_dim1 = *lde;
#line 192 "SB04OW.f"
    e_offset = 1 + e_dim1;
#line 192 "SB04OW.f"
    e -= e_offset;
#line 192 "SB04OW.f"
    f_dim1 = *ldf;
#line 192 "SB04OW.f"
    f_offset = 1 + f_dim1;
#line 192 "SB04OW.f"
    f -= f_offset;
#line 192 "SB04OW.f"
    --iwork;
#line 192 "SB04OW.f"

#line 192 "SB04OW.f"
    /* Function Body */
#line 192 "SB04OW.f"
    *info = 0;
#line 193 "SB04OW.f"
    ierr = 0;
#line 194 "SB04OW.f"
    if (*m <= 0) {
#line 195 "SB04OW.f"
	*info = -1;
#line 196 "SB04OW.f"
    } else if (*n <= 0) {
#line 197 "SB04OW.f"
	*info = -2;
#line 198 "SB04OW.f"
    } else if (*lda < max(1,*m)) {
#line 199 "SB04OW.f"
	*info = -4;
#line 200 "SB04OW.f"
    } else if (*ldb < max(1,*n)) {
#line 201 "SB04OW.f"
	*info = -6;
#line 202 "SB04OW.f"
    } else if (*ldc < max(1,*m)) {
#line 203 "SB04OW.f"
	*info = -8;
#line 204 "SB04OW.f"
    } else if (*ldd < max(1,*m)) {
#line 205 "SB04OW.f"
	*info = -10;
#line 206 "SB04OW.f"
    } else if (*lde < max(1,*n)) {
#line 207 "SB04OW.f"
	*info = -12;
#line 208 "SB04OW.f"
    } else if (*ldf < max(1,*m)) {
#line 209 "SB04OW.f"
	*info = -14;
#line 210 "SB04OW.f"
    }

/*     Return if there were illegal values. */

#line 214 "SB04OW.f"
    if (*info != 0) {
#line 215 "SB04OW.f"
	i__1 = -(*info);
#line 215 "SB04OW.f"
	xerbla_("SB04OW", &i__1, (ftnlen)6);
#line 216 "SB04OW.f"
	return 0;
#line 217 "SB04OW.f"
    }

/*     Determine block structure of A. */

#line 221 "SB04OW.f"
    p = 0;
#line 222 "SB04OW.f"
    i__ = 1;
#line 223 "SB04OW.f"
L10:
#line 224 "SB04OW.f"
    if (i__ > *m) {
#line 224 "SB04OW.f"
	goto L20;
#line 224 "SB04OW.f"
    }
#line 226 "SB04OW.f"
    ++p;
#line 227 "SB04OW.f"
    iwork[p] = i__;
#line 228 "SB04OW.f"
    if (i__ == *m) {
#line 228 "SB04OW.f"
	goto L20;
#line 228 "SB04OW.f"
    }
#line 230 "SB04OW.f"
    if (a[i__ + 1 + i__ * a_dim1] != 0.) {
#line 231 "SB04OW.f"
	i__ += 2;
#line 232 "SB04OW.f"
    } else {
#line 233 "SB04OW.f"
	++i__;
#line 234 "SB04OW.f"
    }
#line 235 "SB04OW.f"
    goto L10;
#line 236 "SB04OW.f"
L20:
#line 237 "SB04OW.f"
    iwork[p + 1] = *m + 1;

/*     Determine block structure of B. */

#line 241 "SB04OW.f"
    q = p + 1;
#line 242 "SB04OW.f"
    j = 1;
#line 243 "SB04OW.f"
L30:
#line 244 "SB04OW.f"
    if (j > *n) {
#line 244 "SB04OW.f"
	goto L40;
#line 244 "SB04OW.f"
    }
#line 246 "SB04OW.f"
    ++q;
#line 247 "SB04OW.f"
    iwork[q] = j;
#line 248 "SB04OW.f"
    if (j == *n) {
#line 248 "SB04OW.f"
	goto L40;
#line 248 "SB04OW.f"
    }
#line 250 "SB04OW.f"
    if (b[j + 1 + j * b_dim1] != 0.) {
#line 251 "SB04OW.f"
	j += 2;
#line 252 "SB04OW.f"
    } else {
#line 253 "SB04OW.f"
	++j;
#line 254 "SB04OW.f"
    }
#line 255 "SB04OW.f"
    goto L30;
#line 256 "SB04OW.f"
L40:
#line 257 "SB04OW.f"
    iwork[q + 1] = *n + 1;

/*     Solve (I, J) - subsystem */
/*       A(I,I) * R(I,J) - L(I,J) * B(J,J) = C(I,J) */
/*       D(I,I) * L(I,J) - R(I,J) * E(J,J) = F(I,J) */
/*     for I = P, P - 1, ..., 1; J = 1, 2, ..., Q. */

#line 264 "SB04OW.f"
    *scale = 1.;
#line 265 "SB04OW.f"
    scaloc = 1.;
#line 266 "SB04OW.f"
    i__1 = q;
#line 266 "SB04OW.f"
    for (j = p + 2; j <= i__1; ++j) {
#line 267 "SB04OW.f"
	js = iwork[j];
#line 268 "SB04OW.f"
	jsp1 = js + 1;
#line 269 "SB04OW.f"
	je = iwork[j + 1] - 1;
#line 270 "SB04OW.f"
	nb = je - js + 1;
#line 271 "SB04OW.f"
	for (i__ = p; i__ >= 1; --i__) {

#line 273 "SB04OW.f"
	    is = iwork[i__];
#line 274 "SB04OW.f"
	    isp1 = is + 1;
#line 275 "SB04OW.f"
	    ie = iwork[i__ + 1] - 1;
#line 276 "SB04OW.f"
	    mb = ie - is + 1;
#line 277 "SB04OW.f"
	    zdim = mb * nb << 1;

#line 279 "SB04OW.f"
	    if (mb == 1 && nb == 1) {

/*              Build a 2-by-2 system Z * x = RHS. */

#line 283 "SB04OW.f"
		z__[0] = a[is + is * a_dim1];
#line 284 "SB04OW.f"
		z__[1] = -e[js + js * e_dim1];
#line 285 "SB04OW.f"
		z__[8] = -b[js + js * b_dim1];
#line 286 "SB04OW.f"
		z__[9] = d__[is + is * d_dim1];

/*              Set up right hand side(s). */

#line 290 "SB04OW.f"
		rhs[0] = c__[is + js * c_dim1];
#line 291 "SB04OW.f"
		rhs[1] = f[is + js * f_dim1];

/*              Solve Z * x = RHS. */

#line 295 "SB04OW.f"
		dgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 296 "SB04OW.f"
		if (ierr > 0) {
#line 296 "SB04OW.f"
		    *info = ierr;
#line 296 "SB04OW.f"
		}

#line 299 "SB04OW.f"
		dgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 300 "SB04OW.f"
		if (scaloc != 1.) {
#line 301 "SB04OW.f"
		    i__2 = *n;
#line 301 "SB04OW.f"
		    for (k = 1; k <= i__2; ++k) {
#line 302 "SB04OW.f"
			dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 303 "SB04OW.f"
			dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 304 "SB04OW.f"
/* L50: */
#line 304 "SB04OW.f"
		    }
#line 305 "SB04OW.f"
		    *scale *= scaloc;
#line 306 "SB04OW.f"
		}

/*              Unpack solution vector(s). */

#line 310 "SB04OW.f"
		c__[is + js * c_dim1] = rhs[0];
#line 311 "SB04OW.f"
		f[is + js * f_dim1] = rhs[1];

/*              Substitute R(I,J) and L(I,J) into remaining equation. */

#line 315 "SB04OW.f"
		if (i__ > 1) {
#line 316 "SB04OW.f"
		    i__2 = is - 1;
#line 316 "SB04OW.f"
		    d__1 = -rhs[0];
#line 316 "SB04OW.f"
		    daxpy_(&i__2, &d__1, &a[is * a_dim1 + 1], &c__1, &c__[js *
			     c_dim1 + 1], &c__1);
#line 317 "SB04OW.f"
		    i__2 = is - 1;
#line 317 "SB04OW.f"
		    d__1 = -rhs[1];
#line 317 "SB04OW.f"
		    daxpy_(&i__2, &d__1, &d__[is * d_dim1 + 1], &c__1, &f[js *
			     f_dim1 + 1], &c__1);
#line 318 "SB04OW.f"
		}
#line 319 "SB04OW.f"
		if (j < q) {
#line 320 "SB04OW.f"
		    i__2 = *n - je;
#line 320 "SB04OW.f"
		    daxpy_(&i__2, &rhs[1], &b[js + (je + 1) * b_dim1], ldb, &
			    c__[is + (je + 1) * c_dim1], ldc);
#line 322 "SB04OW.f"
		    i__2 = *n - je;
#line 322 "SB04OW.f"
		    daxpy_(&i__2, rhs, &e[js + (je + 1) * e_dim1], lde, &f[is 
			    + (je + 1) * f_dim1], ldf);
#line 324 "SB04OW.f"
		}

#line 326 "SB04OW.f"
	    } else if (mb == 1 && nb == 2) {

/*              Build a 4-by-4 system Z * x = RHS. */

#line 330 "SB04OW.f"
		z__[0] = a[is + is * a_dim1];
#line 331 "SB04OW.f"
		z__[1] = 0.;
#line 332 "SB04OW.f"
		z__[2] = -e[js + js * e_dim1];
#line 333 "SB04OW.f"
		z__[3] = -e[js + jsp1 * e_dim1];

#line 335 "SB04OW.f"
		z__[8] = 0.;
#line 336 "SB04OW.f"
		z__[9] = a[is + is * a_dim1];
#line 337 "SB04OW.f"
		z__[10] = 0.;
#line 338 "SB04OW.f"
		z__[11] = -e[jsp1 + jsp1 * e_dim1];

#line 340 "SB04OW.f"
		z__[16] = -b[js + js * b_dim1];
#line 341 "SB04OW.f"
		z__[17] = -b[js + jsp1 * b_dim1];
#line 342 "SB04OW.f"
		z__[18] = d__[is + is * d_dim1];
#line 343 "SB04OW.f"
		z__[19] = 0.;

#line 345 "SB04OW.f"
		z__[24] = -b[jsp1 + js * b_dim1];
#line 346 "SB04OW.f"
		z__[25] = -b[jsp1 + jsp1 * b_dim1];
#line 347 "SB04OW.f"
		z__[26] = 0.;
#line 348 "SB04OW.f"
		z__[27] = d__[is + is * d_dim1];

/*              Set up right hand side(s). */

#line 352 "SB04OW.f"
		rhs[0] = c__[is + js * c_dim1];
#line 353 "SB04OW.f"
		rhs[1] = c__[is + jsp1 * c_dim1];
#line 354 "SB04OW.f"
		rhs[2] = f[is + js * f_dim1];
#line 355 "SB04OW.f"
		rhs[3] = f[is + jsp1 * f_dim1];

/*              Solve Z * x = RHS. */

#line 359 "SB04OW.f"
		dgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 360 "SB04OW.f"
		if (ierr > 0) {
#line 360 "SB04OW.f"
		    *info = ierr;
#line 360 "SB04OW.f"
		}

#line 363 "SB04OW.f"
		dgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 364 "SB04OW.f"
		if (scaloc != 1.) {
#line 365 "SB04OW.f"
		    i__2 = *n;
#line 365 "SB04OW.f"
		    for (k = 1; k <= i__2; ++k) {
#line 366 "SB04OW.f"
			dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 367 "SB04OW.f"
			dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 368 "SB04OW.f"
/* L60: */
#line 368 "SB04OW.f"
		    }
#line 369 "SB04OW.f"
		    *scale *= scaloc;
#line 370 "SB04OW.f"
		}

/*              Unpack solution vector(s). */

#line 374 "SB04OW.f"
		c__[is + js * c_dim1] = rhs[0];
#line 375 "SB04OW.f"
		c__[is + jsp1 * c_dim1] = rhs[1];
#line 376 "SB04OW.f"
		f[is + js * f_dim1] = rhs[2];
#line 377 "SB04OW.f"
		f[is + jsp1 * f_dim1] = rhs[3];

/*              Substitute R(I,J) and L(I,J) into remaining equation. */

#line 381 "SB04OW.f"
		if (i__ > 1) {
#line 382 "SB04OW.f"
		    i__2 = is - 1;
#line 382 "SB04OW.f"
		    dger_(&i__2, &nb, &c_b23, &a[is * a_dim1 + 1], &c__1, rhs,
			     &c__1, &c__[js * c_dim1 + 1], ldc);
#line 384 "SB04OW.f"
		    i__2 = is - 1;
#line 384 "SB04OW.f"
		    dger_(&i__2, &nb, &c_b23, &d__[is * d_dim1 + 1], &c__1, &
			    rhs[2], &c__1, &f[js * f_dim1 + 1], ldf);
#line 386 "SB04OW.f"
		}
#line 387 "SB04OW.f"
		if (j < q) {
#line 388 "SB04OW.f"
		    i__2 = *n - je;
#line 388 "SB04OW.f"
		    daxpy_(&i__2, &rhs[2], &b[js + (je + 1) * b_dim1], ldb, &
			    c__[is + (je + 1) * c_dim1], ldc);
#line 390 "SB04OW.f"
		    i__2 = *n - je;
#line 390 "SB04OW.f"
		    daxpy_(&i__2, rhs, &e[js + (je + 1) * e_dim1], lde, &f[is 
			    + (je + 1) * f_dim1], ldf);
#line 392 "SB04OW.f"
		    i__2 = *n - je;
#line 392 "SB04OW.f"
		    daxpy_(&i__2, &rhs[3], &b[jsp1 + (je + 1) * b_dim1], ldb, 
			    &c__[is + (je + 1) * c_dim1], ldc);
#line 394 "SB04OW.f"
		    i__2 = *n - je;
#line 394 "SB04OW.f"
		    daxpy_(&i__2, &rhs[1], &e[jsp1 + (je + 1) * e_dim1], lde, 
			    &f[is + (je + 1) * f_dim1], ldf);
#line 396 "SB04OW.f"
		}

#line 398 "SB04OW.f"
	    } else if (mb == 2 && nb == 1) {

/*              Build a 4-by-4 system Z * x = RHS. */

#line 402 "SB04OW.f"
		z__[0] = a[is + is * a_dim1];
#line 403 "SB04OW.f"
		z__[1] = a[isp1 + is * a_dim1];
#line 404 "SB04OW.f"
		z__[2] = -e[js + js * e_dim1];
#line 405 "SB04OW.f"
		z__[3] = 0.;

#line 407 "SB04OW.f"
		z__[8] = a[is + isp1 * a_dim1];
#line 408 "SB04OW.f"
		z__[9] = a[isp1 + isp1 * a_dim1];
#line 409 "SB04OW.f"
		z__[10] = 0.;
#line 410 "SB04OW.f"
		z__[11] = -e[js + js * e_dim1];

#line 412 "SB04OW.f"
		z__[16] = -b[js + js * b_dim1];
#line 413 "SB04OW.f"
		z__[17] = 0.;
#line 414 "SB04OW.f"
		z__[18] = d__[is + is * d_dim1];
#line 415 "SB04OW.f"
		z__[19] = 0.;

#line 417 "SB04OW.f"
		z__[24] = 0.;
#line 418 "SB04OW.f"
		z__[25] = -b[js + js * b_dim1];
#line 419 "SB04OW.f"
		z__[26] = d__[is + isp1 * d_dim1];
#line 420 "SB04OW.f"
		z__[27] = d__[isp1 + isp1 * d_dim1];

/*              Set up right hand side(s). */

#line 424 "SB04OW.f"
		rhs[0] = c__[is + js * c_dim1];
#line 425 "SB04OW.f"
		rhs[1] = c__[isp1 + js * c_dim1];
#line 426 "SB04OW.f"
		rhs[2] = f[is + js * f_dim1];
#line 427 "SB04OW.f"
		rhs[3] = f[isp1 + js * f_dim1];

/*              Solve Z * x = RHS. */

#line 431 "SB04OW.f"
		dgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 432 "SB04OW.f"
		if (ierr > 0) {
#line 432 "SB04OW.f"
		    *info = ierr;
#line 432 "SB04OW.f"
		}

#line 435 "SB04OW.f"
		dgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 436 "SB04OW.f"
		if (scaloc != 1.) {
#line 437 "SB04OW.f"
		    i__2 = *n;
#line 437 "SB04OW.f"
		    for (k = 1; k <= i__2; ++k) {
#line 438 "SB04OW.f"
			dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 439 "SB04OW.f"
			dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 440 "SB04OW.f"
/* L70: */
#line 440 "SB04OW.f"
		    }
#line 441 "SB04OW.f"
		    *scale *= scaloc;
#line 442 "SB04OW.f"
		}

/*              Unpack solution vector(s). */

#line 446 "SB04OW.f"
		c__[is + js * c_dim1] = rhs[0];
#line 447 "SB04OW.f"
		c__[isp1 + js * c_dim1] = rhs[1];
#line 448 "SB04OW.f"
		f[is + js * f_dim1] = rhs[2];
#line 449 "SB04OW.f"
		f[isp1 + js * f_dim1] = rhs[3];

/*              Substitute R(I,J) and L(I,J) into remaining equation. */

#line 453 "SB04OW.f"
		if (i__ > 1) {
#line 454 "SB04OW.f"
		    i__2 = is - 1;
#line 454 "SB04OW.f"
		    dgemv_("N", &i__2, &mb, &c_b23, &a[is * a_dim1 + 1], lda, 
			    rhs, &c__1, &c_b37, &c__[js * c_dim1 + 1], &c__1, 
			    (ftnlen)1);
#line 456 "SB04OW.f"
		    i__2 = is - 1;
#line 456 "SB04OW.f"
		    dgemv_("N", &i__2, &mb, &c_b23, &d__[is * d_dim1 + 1], 
			    ldd, &rhs[2], &c__1, &c_b37, &f[js * f_dim1 + 1], 
			    &c__1, (ftnlen)1);
#line 458 "SB04OW.f"
		}
#line 459 "SB04OW.f"
		if (j < q) {
#line 460 "SB04OW.f"
		    i__2 = *n - je;
#line 460 "SB04OW.f"
		    dger_(&mb, &i__2, &c_b37, &rhs[2], &c__1, &b[js + (je + 1)
			     * b_dim1], ldb, &c__[is + (je + 1) * c_dim1], 
			    ldc);
#line 462 "SB04OW.f"
		    i__2 = *n - je;
#line 462 "SB04OW.f"
		    dger_(&mb, &i__2, &c_b37, rhs, &c__1, &e[js + (je + 1) * 
			    e_dim1], lde, &f[is + (je + 1) * f_dim1], ldf);
#line 464 "SB04OW.f"
		}

#line 466 "SB04OW.f"
	    } else if (mb == 2 && nb == 2) {

/*              Build an 8-by-8 system Z * x = RHS. */

#line 470 "SB04OW.f"
		dlaset_("All", &c__8, &c__8, &c_b51, &c_b51, z__, &c__8, (
			ftnlen)3);

#line 472 "SB04OW.f"
		z__[0] = a[is + is * a_dim1];
#line 473 "SB04OW.f"
		z__[1] = a[isp1 + is * a_dim1];
#line 474 "SB04OW.f"
		z__[4] = -e[js + js * e_dim1];
#line 475 "SB04OW.f"
		z__[6] = -e[js + jsp1 * e_dim1];

#line 477 "SB04OW.f"
		z__[8] = a[is + isp1 * a_dim1];
#line 478 "SB04OW.f"
		z__[9] = a[isp1 + isp1 * a_dim1];
#line 479 "SB04OW.f"
		z__[13] = -e[js + js * e_dim1];
#line 480 "SB04OW.f"
		z__[15] = -e[js + jsp1 * e_dim1];

#line 482 "SB04OW.f"
		z__[18] = a[is + is * a_dim1];
#line 483 "SB04OW.f"
		z__[19] = a[isp1 + is * a_dim1];
#line 484 "SB04OW.f"
		z__[22] = -e[jsp1 + jsp1 * e_dim1];

#line 486 "SB04OW.f"
		z__[26] = a[is + isp1 * a_dim1];
#line 487 "SB04OW.f"
		z__[27] = a[isp1 + isp1 * a_dim1];
#line 488 "SB04OW.f"
		z__[31] = -e[jsp1 + jsp1 * e_dim1];

#line 490 "SB04OW.f"
		z__[32] = -b[js + js * b_dim1];
#line 491 "SB04OW.f"
		z__[34] = -b[js + jsp1 * b_dim1];
#line 492 "SB04OW.f"
		z__[36] = d__[is + is * d_dim1];

#line 494 "SB04OW.f"
		z__[41] = -b[js + js * b_dim1];
#line 495 "SB04OW.f"
		z__[43] = -b[js + jsp1 * b_dim1];
#line 496 "SB04OW.f"
		z__[44] = d__[is + isp1 * d_dim1];
#line 497 "SB04OW.f"
		z__[45] = d__[isp1 + isp1 * d_dim1];

#line 499 "SB04OW.f"
		z__[48] = -b[jsp1 + js * b_dim1];
#line 500 "SB04OW.f"
		z__[50] = -b[jsp1 + jsp1 * b_dim1];
#line 501 "SB04OW.f"
		z__[54] = d__[is + is * d_dim1];

#line 503 "SB04OW.f"
		z__[57] = -b[jsp1 + js * b_dim1];
#line 504 "SB04OW.f"
		z__[59] = -b[jsp1 + jsp1 * b_dim1];

#line 506 "SB04OW.f"
		z__[62] = d__[is + isp1 * d_dim1];
#line 507 "SB04OW.f"
		z__[63] = d__[isp1 + isp1 * d_dim1];

/*              Set up right hand side(s). */

#line 511 "SB04OW.f"
		k = 1;
#line 512 "SB04OW.f"
		ii = mb * nb + 1;
#line 513 "SB04OW.f"
		i__2 = nb - 1;
#line 513 "SB04OW.f"
		for (jj = 0; jj <= i__2; ++jj) {
#line 514 "SB04OW.f"
		    dcopy_(&mb, &c__[is + (js + jj) * c_dim1], &c__1, &rhs[k 
			    - 1], &c__1);
#line 515 "SB04OW.f"
		    dcopy_(&mb, &f[is + (js + jj) * f_dim1], &c__1, &rhs[ii - 
			    1], &c__1);
#line 516 "SB04OW.f"
		    k += mb;
#line 517 "SB04OW.f"
		    ii += mb;
#line 518 "SB04OW.f"
/* L80: */
#line 518 "SB04OW.f"
		}

/*              Solve Z * x = RHS. */

#line 522 "SB04OW.f"
		dgetc2_(&zdim, z__, &c__8, ipiv, jpiv, &ierr);
#line 523 "SB04OW.f"
		if (ierr > 0) {
#line 523 "SB04OW.f"
		    *info = ierr;
#line 523 "SB04OW.f"
		}

#line 526 "SB04OW.f"
		dgesc2_(&zdim, z__, &c__8, rhs, ipiv, jpiv, &scaloc);
#line 527 "SB04OW.f"
		if (scaloc != 1.) {
#line 528 "SB04OW.f"
		    i__2 = *n;
#line 528 "SB04OW.f"
		    for (k = 1; k <= i__2; ++k) {
#line 529 "SB04OW.f"
			dscal_(m, &scaloc, &c__[k * c_dim1 + 1], &c__1);
#line 530 "SB04OW.f"
			dscal_(m, &scaloc, &f[k * f_dim1 + 1], &c__1);
#line 531 "SB04OW.f"
/* L90: */
#line 531 "SB04OW.f"
		    }
#line 532 "SB04OW.f"
		    *scale *= scaloc;
#line 533 "SB04OW.f"
		}

/*              Unpack solution vector(s). */

#line 537 "SB04OW.f"
		k = 1;
#line 538 "SB04OW.f"
		ii = mb * nb + 1;
#line 539 "SB04OW.f"
		i__2 = nb - 1;
#line 539 "SB04OW.f"
		for (jj = 0; jj <= i__2; ++jj) {
#line 540 "SB04OW.f"
		    dcopy_(&mb, &rhs[k - 1], &c__1, &c__[is + (js + jj) * 
			    c_dim1], &c__1);
#line 541 "SB04OW.f"
		    dcopy_(&mb, &rhs[ii - 1], &c__1, &f[is + (js + jj) * 
			    f_dim1], &c__1);
#line 542 "SB04OW.f"
		    k += mb;
#line 543 "SB04OW.f"
		    ii += mb;
#line 544 "SB04OW.f"
/* L100: */
#line 544 "SB04OW.f"
		}

/*              Substitute R(I,J) and L(I,J) into remaining equation. */

#line 548 "SB04OW.f"
		k = mb * nb + 1;
#line 549 "SB04OW.f"
		if (i__ > 1) {
#line 550 "SB04OW.f"
		    i__2 = is - 1;
#line 550 "SB04OW.f"
		    dgemm_("N", "N", &i__2, &nb, &mb, &c_b23, &a[is * a_dim1 
			    + 1], lda, rhs, &mb, &c_b37, &c__[js * c_dim1 + 1]
			    , ldc, (ftnlen)1, (ftnlen)1);
#line 552 "SB04OW.f"
		    i__2 = is - 1;
#line 552 "SB04OW.f"
		    dgemm_("N", "N", &i__2, &nb, &mb, &c_b23, &d__[is * 
			    d_dim1 + 1], ldd, &rhs[k - 1], &mb, &c_b37, &f[js 
			    * f_dim1 + 1], ldf, (ftnlen)1, (ftnlen)1);
#line 554 "SB04OW.f"
		}
#line 555 "SB04OW.f"
		if (j < q) {
#line 556 "SB04OW.f"
		    i__2 = *n - je;
#line 556 "SB04OW.f"
		    dgemm_("N", "N", &mb, &i__2, &nb, &c_b37, &rhs[k - 1], &
			    mb, &b[js + (je + 1) * b_dim1], ldb, &c_b37, &c__[
			    is + (je + 1) * c_dim1], ldc, (ftnlen)1, (ftnlen)
			    1);
#line 558 "SB04OW.f"
		    i__2 = *n - je;
#line 558 "SB04OW.f"
		    dgemm_("N", "N", &mb, &i__2, &nb, &c_b37, rhs, &mb, &e[js 
			    + (je + 1) * e_dim1], lde, &c_b37, &f[is + (je + 
			    1) * f_dim1], ldf, (ftnlen)1, (ftnlen)1);
#line 560 "SB04OW.f"
		}

#line 562 "SB04OW.f"
	    }

#line 564 "SB04OW.f"
/* L110: */
#line 564 "SB04OW.f"
	}
#line 565 "SB04OW.f"
/* L120: */
#line 565 "SB04OW.f"
    }
#line 566 "SB04OW.f"
    return 0;
/* *** Last line of SB04OW *** */
} /* sb04ow_ */

