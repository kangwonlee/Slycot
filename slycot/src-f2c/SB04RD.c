#line 1 "SB04RD.f"
/* SB04RD.f -- translated by f2c (version 20100827).
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

#line 1 "SB04RD.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int sb04rd_(char *abschu, char *ula, char *ulb, integer *n, 
	integer *m, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *c__, integer *ldc, doublereal *tol, integer *iwork, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen abschu_len, 
	ftnlen ula_len, ftnlen ulb_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1;

    /* Local variables */
    static integer i__, fwd, ldw;
    static doublereal tol1;
    static integer ibeg, iend, incr;
    static logical lula, lulb;
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer maxmn;
    extern /* Subroutine */ int sb04py_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen), sb04rv_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, ftnlen, 
	    ftnlen), sb04rw_(char *, char *, integer *, integer *, doublereal 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, ftnlen, ftnlen);
    static integer istep;
    extern /* Subroutine */ int sb04rx_(char *, char *, integer *, doublereal 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, ftnlen, ftnlen), sb04ry_(char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static integer jwork;
    static logical labscb;
    extern doublereal dlamch_(char *, ftnlen);
    static char abschr[1];
    static logical labscs;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer ipincr;


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

/*     To solve for X the discrete-time Sylvester equation */

/*        X + AXB = C, */

/*     with at least one of the matrices A or B in Schur form and the */
/*     other in Hessenberg or Schur form (both either upper or lower); */
/*     A, B, C and X are N-by-N, M-by-M, N-by-M, and N-by-M matrices, */
/*     respectively. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     ABSCHU  CHARACTER*1 */
/*             Indicates whether A and/or B is/are in Schur or */
/*             Hessenberg form as follows: */
/*             = 'A':  A is in Schur form, B is in Hessenberg form; */
/*             = 'B':  B is in Schur form, A is in Hessenberg form; */
/*             = 'S':  Both A and B are in Schur form. */

/*     ULA     CHARACTER*1 */
/*             Indicates whether A is in upper or lower Schur form or */
/*             upper or lower Hessenberg form as follows: */
/*             = 'U':  A is in upper Hessenberg form if ABSCHU = 'B' and */
/*                     upper Schur form otherwise; */
/*             = 'L':  A is in lower Hessenberg form if ABSCHU = 'B' and */
/*                     lower Schur form otherwise. */

/*     ULB     CHARACTER*1 */
/*             Indicates whether B is in upper or lower Schur form or */
/*             upper or lower Hessenberg form as follows: */
/*             = 'U':  B is in upper Hessenberg form if ABSCHU = 'A' and */
/*                     upper Schur form otherwise; */
/*             = 'L':  B is in lower Hessenberg form if ABSCHU = 'A' and */
/*                     lower Schur form otherwise. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The order of the matrix B.  M >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             coefficient matrix A of the equation. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading M-by-M part of this array must contain the */
/*             coefficient matrix B of the equation. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,M). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the coefficient matrix C of the equation. */
/*             On exit, if INFO = 0, the leading N-by-M part of this */
/*             array contains the solution matrix X of the problem. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,N). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used to test for near singularity in */
/*             the Sylvester equation. If the user sets TOL > 0, then the */
/*             given value of TOL is used as a lower bound for the */
/*             reciprocal condition number; a matrix whose estimated */
/*             condition number is less than 1/TOL is considered to be */
/*             nonsingular. If the user sets TOL <= 0, then a default */
/*             tolerance, defined by TOLDEF = EPS, is used instead, where */
/*             EPS is the machine precision (see LAPACK Library routine */
/*             DLAMCH). */
/*             This parameter is not referenced if ABSCHU = 'S', */
/*             ULA = 'U', and ULB = 'U'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (2*MAX(M,N)) */
/*             This parameter is not referenced if ABSCHU = 'S', */
/*             ULA = 'U', and ULB = 'U'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK = 2*N, if ABSCHU = 'S', ULA = 'U', and ULB = 'U'; */
/*             LDWORK = 2*MAX(M,N)*(4 + 2*MAX(M,N)), otherwise. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if a (numerically) singular matrix T was encountered */
/*                   during the computation of the solution matrix X. */
/*                   That is, the estimated reciprocal condition number */
/*                   of T is less than or equal to TOL. */

/*     METHOD */

/*     Matrices A and B are assumed to be in (upper or lower) Hessenberg */
/*     or Schur form (with at least one of them in Schur form). The */
/*     solution matrix X is then computed by rows or columns via the back */
/*     substitution scheme proposed by Golub, Nash and Van Loan (see */
/*     [1]), which involves the solution of triangular systems of */
/*     equations that are constructed recursively and which may be nearly */
/*     singular if A and -B have almost reciprocal eigenvalues. If near */
/*     singularity is detected, then the routine returns with the Error */
/*     Indicator (INFO) set to 1. */

/*     REFERENCES */

/*     [1] Golub, G.H., Nash, S. and Van Loan, C.F. */
/*         A Hessenberg-Schur method for the problem AX + XB = C. */
/*         IEEE Trans. Auto. Contr., AC-24, pp. 909-913, 1979. */

/*     [2] Sima, V. */
/*         Algorithms for Linear-quadratic Optimization. */
/*         Marcel Dekker, Inc., New York, 1996. */

/*     NUMERICAL ASPECTS */
/*                                            2         2 */
/*     The algorithm requires approximately 5M N + 0.5MN  operations in */
/*                            2         2 */
/*     the worst case and 2.5M N + 0.5MN  operations in the best case */
/*     (where M is the order of the matrix in Hessenberg form and N is */
/*     the order of the matrix in Schur form) and is mixed stable (see */
/*     [1]). */

/*     CONTRIBUTORS */

/*     D. Sima, University of Bucharest, May 2000. */

/*     REVISIONS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, June 2000. */

/*     KEYWORDS */

/*     Hessenberg form, orthogonal transformation, real Schur form, */
/*     Sylvester equation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 203 "SB04RD.f"
    /* Parameter adjustments */
#line 203 "SB04RD.f"
    a_dim1 = *lda;
#line 203 "SB04RD.f"
    a_offset = 1 + a_dim1;
#line 203 "SB04RD.f"
    a -= a_offset;
#line 203 "SB04RD.f"
    b_dim1 = *ldb;
#line 203 "SB04RD.f"
    b_offset = 1 + b_dim1;
#line 203 "SB04RD.f"
    b -= b_offset;
#line 203 "SB04RD.f"
    c_dim1 = *ldc;
#line 203 "SB04RD.f"
    c_offset = 1 + c_dim1;
#line 203 "SB04RD.f"
    c__ -= c_offset;
#line 203 "SB04RD.f"
    --iwork;
#line 203 "SB04RD.f"
    --dwork;
#line 203 "SB04RD.f"

#line 203 "SB04RD.f"
    /* Function Body */
#line 203 "SB04RD.f"
    *info = 0;
#line 204 "SB04RD.f"
    maxmn = max(*m,*n);
#line 205 "SB04RD.f"
    labscb = lsame_(abschu, "B", (ftnlen)1, (ftnlen)1);
#line 206 "SB04RD.f"
    labscs = lsame_(abschu, "S", (ftnlen)1, (ftnlen)1);
#line 207 "SB04RD.f"
    lula = lsame_(ula, "U", (ftnlen)1, (ftnlen)1);
#line 208 "SB04RD.f"
    lulb = lsame_(ulb, "U", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 212 "SB04RD.f"
    if (! labscb && ! labscs && ! lsame_(abschu, "A", (ftnlen)1, (ftnlen)1)) {
#line 214 "SB04RD.f"
	*info = -1;
#line 215 "SB04RD.f"
    } else if (! lula && ! lsame_(ula, "L", (ftnlen)1, (ftnlen)1)) {
#line 216 "SB04RD.f"
	*info = -2;
#line 217 "SB04RD.f"
    } else if (! lulb && ! lsame_(ulb, "L", (ftnlen)1, (ftnlen)1)) {
#line 218 "SB04RD.f"
	*info = -3;
#line 219 "SB04RD.f"
    } else if (*n < 0) {
#line 220 "SB04RD.f"
	*info = -4;
#line 221 "SB04RD.f"
    } else if (*m < 0) {
#line 222 "SB04RD.f"
	*info = -5;
#line 223 "SB04RD.f"
    } else if (*lda < max(1,*n)) {
#line 224 "SB04RD.f"
	*info = -7;
#line 225 "SB04RD.f"
    } else if (*ldb < max(1,*m)) {
#line 226 "SB04RD.f"
	*info = -9;
#line 227 "SB04RD.f"
    } else if (*ldc < max(1,*n)) {
#line 228 "SB04RD.f"
	*info = -11;
#line 229 "SB04RD.f"
    } else if (*ldwork < *n << 1 || *ldwork < (maxmn << 1) * ((maxmn << 1) + 
	    4) && ! (labscs && lula && lulb)) {
#line 232 "SB04RD.f"
	*info = -15;
#line 233 "SB04RD.f"
    }

#line 235 "SB04RD.f"
    if (*info != 0) {

/*        Error return. */

#line 239 "SB04RD.f"
	i__1 = -(*info);
#line 239 "SB04RD.f"
	xerbla_("SB04RD", &i__1, (ftnlen)6);
#line 240 "SB04RD.f"
	return 0;
#line 241 "SB04RD.f"
    }

/*     Quick return if possible. */

#line 245 "SB04RD.f"
    if (maxmn == 0) {
#line 245 "SB04RD.f"
	return 0;
#line 245 "SB04RD.f"
    }

#line 248 "SB04RD.f"
    if (labscs && lula && lulb) {

/*        If both matrices are in a real Schur form, use SB04PY. */

#line 252 "SB04RD.f"
	sb04py_("NoTranspose", "NoTranspose", &c__1, n, m, &a[a_offset], lda, 
		&b[b_offset], ldb, &c__[c_offset], ldc, &scale, &dwork[1], 
		info, (ftnlen)11, (ftnlen)11);
#line 254 "SB04RD.f"
	if (scale != 1.) {
#line 254 "SB04RD.f"
	    *info = 1;
#line 254 "SB04RD.f"
	}
#line 256 "SB04RD.f"
	return 0;
#line 257 "SB04RD.f"
    }

#line 259 "SB04RD.f"
    ldw = maxmn << 1;
#line 260 "SB04RD.f"
    jwork = ldw * ldw + ldw * 3 + 1;
#line 261 "SB04RD.f"
    tol1 = *tol;
#line 262 "SB04RD.f"
    if (tol1 <= 0.) {
#line 262 "SB04RD.f"
	tol1 = dlamch_("Epsilon", (ftnlen)7);
#line 262 "SB04RD.f"
    }

/*     Choose the smallest of both matrices as the one in Hessenberg */
/*     form when possible. */

#line 268 "SB04RD.f"
    *(unsigned char *)abschr = *(unsigned char *)abschu;
#line 269 "SB04RD.f"
    if (labscs) {
#line 270 "SB04RD.f"
	if (*n > *m) {
#line 271 "SB04RD.f"
	    *(unsigned char *)abschr = 'A';
#line 272 "SB04RD.f"
	} else {
#line 273 "SB04RD.f"
	    *(unsigned char *)abschr = 'B';
#line 274 "SB04RD.f"
	}
#line 275 "SB04RD.f"
    }
#line 276 "SB04RD.f"
    if (lsame_(abschr, "B", (ftnlen)1, (ftnlen)1)) {

/*        B is in Schur form: recursion on the columns of B. */

#line 280 "SB04RD.f"
	if (lulb) {

/*           B is upper: forward recursion. */

#line 284 "SB04RD.f"
	    ibeg = 1;
#line 285 "SB04RD.f"
	    iend = *m;
#line 286 "SB04RD.f"
	    fwd = 1;
#line 287 "SB04RD.f"
	    incr = 0;
#line 288 "SB04RD.f"
	} else {

/*           B is lower: backward recursion. */

#line 292 "SB04RD.f"
	    ibeg = *m;
#line 293 "SB04RD.f"
	    iend = 1;
#line 294 "SB04RD.f"
	    fwd = -1;
#line 295 "SB04RD.f"
	    incr = -1;
#line 296 "SB04RD.f"
	}
#line 297 "SB04RD.f"
	i__ = ibeg;
/*        WHILE ( ( IEND - I ) * FWD .GE. 0 ) DO */
#line 299 "SB04RD.f"
L20:
#line 299 "SB04RD.f"
	if ((iend - i__) * fwd >= 0) {

/*           Test for 1-by-1 or 2-by-2 diagonal block in the Schur */
/*           form. */

#line 304 "SB04RD.f"
	    if (i__ == iend) {
#line 305 "SB04RD.f"
		istep = 1;
#line 306 "SB04RD.f"
	    } else {
#line 307 "SB04RD.f"
		if (b[i__ + fwd + i__ * b_dim1] == 0.) {
#line 308 "SB04RD.f"
		    istep = 1;
#line 309 "SB04RD.f"
		} else {
#line 310 "SB04RD.f"
		    istep = 2;
#line 311 "SB04RD.f"
		}
#line 312 "SB04RD.f"
	    }

#line 314 "SB04RD.f"
	    if (istep == 1) {
#line 315 "SB04RD.f"
		sb04rw_(abschr, ulb, n, m, &c__[c_offset], ldc, &i__, &b[
			b_offset], ldb, &a[a_offset], lda, &dwork[jwork], &
			dwork[1], (ftnlen)1, (ftnlen)1);
#line 317 "SB04RD.f"
		sb04ry_("R", ula, n, &a[a_offset], lda, &b[i__ + i__ * b_dim1]
			, &dwork[jwork], &tol1, &iwork[1], &dwork[1], &ldw, 
			info, (ftnlen)1, (ftnlen)1);
#line 319 "SB04RD.f"
		if (*info == 1) {
#line 319 "SB04RD.f"
		    return 0;
#line 319 "SB04RD.f"
		}
#line 321 "SB04RD.f"
		dcopy_(n, &dwork[jwork], &c__1, &c__[i__ * c_dim1 + 1], &c__1)
			;
#line 322 "SB04RD.f"
	    } else {
#line 323 "SB04RD.f"
		ipincr = i__ + incr;
#line 324 "SB04RD.f"
		sb04rv_(abschr, ulb, n, m, &c__[c_offset], ldc, &ipincr, &b[
			b_offset], ldb, &a[a_offset], lda, &dwork[jwork], &
			dwork[1], (ftnlen)1, (ftnlen)1);
#line 326 "SB04RD.f"
		sb04rx_("R", ula, n, &a[a_offset], lda, &b[ipincr + ipincr * 
			b_dim1], &b[ipincr + 1 + ipincr * b_dim1], &b[ipincr 
			+ (ipincr + 1) * b_dim1], &b[ipincr + 1 + (ipincr + 1)
			 * b_dim1], &dwork[jwork], &tol1, &iwork[1], &dwork[1]
			, &ldw, info, (ftnlen)1, (ftnlen)1);
#line 330 "SB04RD.f"
		if (*info == 1) {
#line 330 "SB04RD.f"
		    return 0;
#line 330 "SB04RD.f"
		}
#line 332 "SB04RD.f"
		dcopy_(n, &dwork[jwork], &c__2, &c__[ipincr * c_dim1 + 1], &
			c__1);
#line 333 "SB04RD.f"
		dcopy_(n, &dwork[jwork + 1], &c__2, &c__[(ipincr + 1) * 
			c_dim1 + 1], &c__1);
#line 334 "SB04RD.f"
	    }
#line 335 "SB04RD.f"
	    i__ += fwd * istep;
#line 336 "SB04RD.f"
	    goto L20;
#line 337 "SB04RD.f"
	}
/*        END WHILE 20 */
#line 339 "SB04RD.f"
    } else {

/*        A is in Schur form: recursion on the rows of A. */

#line 343 "SB04RD.f"
	if (lula) {

/*           A is upper: backward recursion. */

#line 347 "SB04RD.f"
	    ibeg = *n;
#line 348 "SB04RD.f"
	    iend = 1;
#line 349 "SB04RD.f"
	    fwd = -1;
#line 350 "SB04RD.f"
	    incr = -1;
#line 351 "SB04RD.f"
	} else {

/*           A is lower: forward recursion. */

#line 355 "SB04RD.f"
	    ibeg = 1;
#line 356 "SB04RD.f"
	    iend = *n;
#line 357 "SB04RD.f"
	    fwd = 1;
#line 358 "SB04RD.f"
	    incr = 0;
#line 359 "SB04RD.f"
	}
#line 360 "SB04RD.f"
	i__ = ibeg;
/*        WHILE ( ( IEND - I ) * FWD .GE. 0 ) DO */
#line 362 "SB04RD.f"
L40:
#line 362 "SB04RD.f"
	if ((iend - i__) * fwd >= 0) {

/*           Test for 1-by-1 or 2-by-2 diagonal block in the Schur */
/*           form. */

#line 367 "SB04RD.f"
	    if (i__ == iend) {
#line 368 "SB04RD.f"
		istep = 1;
#line 369 "SB04RD.f"
	    } else {
#line 370 "SB04RD.f"
		if (a[i__ + (i__ + fwd) * a_dim1] == 0.) {
#line 371 "SB04RD.f"
		    istep = 1;
#line 372 "SB04RD.f"
		} else {
#line 373 "SB04RD.f"
		    istep = 2;
#line 374 "SB04RD.f"
		}
#line 375 "SB04RD.f"
	    }

#line 377 "SB04RD.f"
	    if (istep == 1) {
#line 378 "SB04RD.f"
		sb04rw_(abschr, ula, n, m, &c__[c_offset], ldc, &i__, &a[
			a_offset], lda, &b[b_offset], ldb, &dwork[jwork], &
			dwork[1], (ftnlen)1, (ftnlen)1);
#line 380 "SB04RD.f"
		sb04ry_("C", ulb, m, &b[b_offset], ldb, &a[i__ + i__ * a_dim1]
			, &dwork[jwork], &tol1, &iwork[1], &dwork[1], &ldw, 
			info, (ftnlen)1, (ftnlen)1);
#line 382 "SB04RD.f"
		if (*info == 1) {
#line 382 "SB04RD.f"
		    return 0;
#line 382 "SB04RD.f"
		}
#line 384 "SB04RD.f"
		dcopy_(m, &dwork[jwork], &c__1, &c__[i__ + c_dim1], ldc);
#line 385 "SB04RD.f"
	    } else {
#line 386 "SB04RD.f"
		ipincr = i__ + incr;
#line 387 "SB04RD.f"
		sb04rv_(abschr, ula, n, m, &c__[c_offset], ldc, &ipincr, &a[
			a_offset], lda, &b[b_offset], ldb, &dwork[jwork], &
			dwork[1], (ftnlen)1, (ftnlen)1);
#line 389 "SB04RD.f"
		sb04rx_("C", ulb, m, &b[b_offset], ldb, &a[ipincr + ipincr * 
			a_dim1], &a[ipincr + 1 + ipincr * a_dim1], &a[ipincr 
			+ (ipincr + 1) * a_dim1], &a[ipincr + 1 + (ipincr + 1)
			 * a_dim1], &dwork[jwork], &tol1, &iwork[1], &dwork[1]
			, &ldw, info, (ftnlen)1, (ftnlen)1);
#line 393 "SB04RD.f"
		if (*info == 1) {
#line 393 "SB04RD.f"
		    return 0;
#line 393 "SB04RD.f"
		}
#line 395 "SB04RD.f"
		dcopy_(m, &dwork[jwork], &c__2, &c__[ipincr + c_dim1], ldc);
#line 396 "SB04RD.f"
		dcopy_(m, &dwork[jwork + 1], &c__2, &c__[ipincr + 1 + c_dim1],
			 ldc);
#line 397 "SB04RD.f"
	    }
#line 398 "SB04RD.f"
	    i__ += fwd * istep;
#line 399 "SB04RD.f"
	    goto L40;
#line 400 "SB04RD.f"
	}
/*        END WHILE 40 */
#line 402 "SB04RD.f"
    }

#line 404 "SB04RD.f"
    return 0;
/* *** Last line of SB04RD *** */
} /* sb04rd_ */
