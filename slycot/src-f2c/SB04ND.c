#line 1 "SB04ND.f"
/* SB04ND.f -- translated by f2c (version 20100827).
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

#line 1 "SB04ND.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int sb04nd_(char *abschu, char *ula, char *ulb, integer *n, 
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
    extern /* Subroutine */ int sb04nv_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, ftnlen, ftnlen), sb04nw_(char *, char *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), sb04nx_(char *,
	     char *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static integer maxmn;
    extern /* Subroutine */ int sb04ny_(char *, char *, integer *, doublereal 
	    *, integer *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer istep, jwork;
    static logical labscb;
    extern doublereal dlamch_(char *, ftnlen);
    static char abschr[1];
    static logical labscs;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer ipincr;
    extern /* Subroutine */ int dtrsyl_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen);


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

/*     To solve for X the continuous-time Sylvester equation */

/*        AX + XB = C, */

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
/*             This parameter is not referenced if ABSCHU = 'S', */
/*             ULA = 'U', and ULB = 'U'. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK = 0, if ABSCHU = 'S', ULA = 'U', and ULB = 'U'; */
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
/*     singular if A and -B have close eigenvalues. If near singularity */
/*     is detected, then the routine returns with the Error Indicator */
/*     (INFO) set to 1. */

/*     REFERENCES */

/*     [1] Golub, G.H., Nash, S. and Van Loan, C.F. */
/*         A Hessenberg-Schur method for the problem AX + XB = C. */
/*         IEEE Trans. Auto. Contr., AC-24, pp. 909-913, 1979. */

/*     NUMERICAL ASPECTS */
/*                                            2         2 */
/*     The algorithm requires approximately 5M N + 0.5MN  operations in */
/*                            2         2 */
/*     the worst case and 2.5M N + 0.5MN  operations in the best case */
/*     (where M is the order of the matrix in Hessenberg form and N is */
/*     the order of the matrix in Schur form) and is mixed stable (see */
/*     [1]). */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB04BD by M. Vanbegin, and */
/*     P. Van Dooren, Philips Research Laboratory, Brussels, Belgium. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000. */

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

#line 203 "SB04ND.f"
    /* Parameter adjustments */
#line 203 "SB04ND.f"
    a_dim1 = *lda;
#line 203 "SB04ND.f"
    a_offset = 1 + a_dim1;
#line 203 "SB04ND.f"
    a -= a_offset;
#line 203 "SB04ND.f"
    b_dim1 = *ldb;
#line 203 "SB04ND.f"
    b_offset = 1 + b_dim1;
#line 203 "SB04ND.f"
    b -= b_offset;
#line 203 "SB04ND.f"
    c_dim1 = *ldc;
#line 203 "SB04ND.f"
    c_offset = 1 + c_dim1;
#line 203 "SB04ND.f"
    c__ -= c_offset;
#line 203 "SB04ND.f"
    --iwork;
#line 203 "SB04ND.f"
    --dwork;
#line 203 "SB04ND.f"

#line 203 "SB04ND.f"
    /* Function Body */
#line 203 "SB04ND.f"
    *info = 0;
#line 204 "SB04ND.f"
    maxmn = max(*m,*n);
#line 205 "SB04ND.f"
    labscb = lsame_(abschu, "B", (ftnlen)1, (ftnlen)1);
#line 206 "SB04ND.f"
    labscs = lsame_(abschu, "S", (ftnlen)1, (ftnlen)1);
#line 207 "SB04ND.f"
    lula = lsame_(ula, "U", (ftnlen)1, (ftnlen)1);
#line 208 "SB04ND.f"
    lulb = lsame_(ulb, "U", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 212 "SB04ND.f"
    if (! labscb && ! labscs && ! lsame_(abschu, "A", (ftnlen)1, (ftnlen)1)) {
#line 214 "SB04ND.f"
	*info = -1;
#line 215 "SB04ND.f"
    } else if (! lula && ! lsame_(ula, "L", (ftnlen)1, (ftnlen)1)) {
#line 216 "SB04ND.f"
	*info = -2;
#line 217 "SB04ND.f"
    } else if (! lulb && ! lsame_(ulb, "L", (ftnlen)1, (ftnlen)1)) {
#line 218 "SB04ND.f"
	*info = -3;
#line 219 "SB04ND.f"
    } else if (*n < 0) {
#line 220 "SB04ND.f"
	*info = -4;
#line 221 "SB04ND.f"
    } else if (*m < 0) {
#line 222 "SB04ND.f"
	*info = -5;
#line 223 "SB04ND.f"
    } else if (*lda < max(1,*n)) {
#line 224 "SB04ND.f"
	*info = -7;
#line 225 "SB04ND.f"
    } else if (*ldb < max(1,*m)) {
#line 226 "SB04ND.f"
	*info = -9;
#line 227 "SB04ND.f"
    } else if (*ldc < max(1,*n)) {
#line 228 "SB04ND.f"
	*info = -11;
#line 229 "SB04ND.f"
    } else if (*ldwork < 0 || ! (labscs && lula && lulb) && *ldwork < (maxmn 
	    << 1) * ((maxmn << 1) + 4)) {
#line 231 "SB04ND.f"
	*info = -15;
#line 232 "SB04ND.f"
    }

#line 234 "SB04ND.f"
    if (*info != 0) {

/*        Error return. */

#line 238 "SB04ND.f"
	i__1 = -(*info);
#line 238 "SB04ND.f"
	xerbla_("SB04ND", &i__1, (ftnlen)6);
#line 239 "SB04ND.f"
	return 0;
#line 240 "SB04ND.f"
    }

/*     Quick return if possible. */

#line 244 "SB04ND.f"
    if (maxmn == 0) {
#line 244 "SB04ND.f"
	return 0;
#line 244 "SB04ND.f"
    }

#line 247 "SB04ND.f"
    if (labscs && lula && lulb) {

/*        If both matrices are in a real Schur form, use DTRSYL. */

#line 251 "SB04ND.f"
	dtrsyl_("NoTranspose", "NoTranspose", &c__1, n, m, &a[a_offset], lda, 
		&b[b_offset], ldb, &c__[c_offset], ldc, &scale, info, (ftnlen)
		11, (ftnlen)11);
#line 253 "SB04ND.f"
	if (scale != 1.) {
#line 253 "SB04ND.f"
	    *info = 1;
#line 253 "SB04ND.f"
	}
#line 255 "SB04ND.f"
	return 0;
#line 256 "SB04ND.f"
    }

#line 258 "SB04ND.f"
    ldw = maxmn << 1;
#line 259 "SB04ND.f"
    jwork = ldw * ldw + ldw * 3 + 1;
#line 260 "SB04ND.f"
    tol1 = *tol;
#line 261 "SB04ND.f"
    if (tol1 <= 0.) {
#line 261 "SB04ND.f"
	tol1 = dlamch_("Epsilon", (ftnlen)7);
#line 261 "SB04ND.f"
    }

/*     Choose the smallest of both matrices as the one in Hessenberg */
/*     form when possible. */

#line 267 "SB04ND.f"
    *(unsigned char *)abschr = *(unsigned char *)abschu;
#line 268 "SB04ND.f"
    if (labscs) {
#line 269 "SB04ND.f"
	if (*n > *m) {
#line 270 "SB04ND.f"
	    *(unsigned char *)abschr = 'A';
#line 271 "SB04ND.f"
	} else {
#line 272 "SB04ND.f"
	    *(unsigned char *)abschr = 'B';
#line 273 "SB04ND.f"
	}
#line 274 "SB04ND.f"
    }
#line 275 "SB04ND.f"
    if (lsame_(abschr, "B", (ftnlen)1, (ftnlen)1)) {

/*        B is in Schur form: recursion on the columns of B. */

#line 279 "SB04ND.f"
	if (lulb) {

/*           B is upper: forward recursion. */

#line 283 "SB04ND.f"
	    ibeg = 1;
#line 284 "SB04ND.f"
	    iend = *m;
#line 285 "SB04ND.f"
	    fwd = 1;
#line 286 "SB04ND.f"
	    incr = 0;
#line 287 "SB04ND.f"
	} else {

/*           B is lower: backward recursion. */

#line 291 "SB04ND.f"
	    ibeg = *m;
#line 292 "SB04ND.f"
	    iend = 1;
#line 293 "SB04ND.f"
	    fwd = -1;
#line 294 "SB04ND.f"
	    incr = -1;
#line 295 "SB04ND.f"
	}
#line 296 "SB04ND.f"
	i__ = ibeg;
/*        WHILE ( ( IEND - I ) * FWD .GE. 0 ) DO */
#line 298 "SB04ND.f"
L20:
#line 298 "SB04ND.f"
	if ((iend - i__) * fwd >= 0) {

/*           Test for 1-by-1 or 2-by-2 diagonal block in the Schur */
/*           form. */

#line 303 "SB04ND.f"
	    if (i__ == iend) {
#line 304 "SB04ND.f"
		istep = 1;
#line 305 "SB04ND.f"
	    } else {
#line 306 "SB04ND.f"
		if (b[i__ + fwd + i__ * b_dim1] == 0.) {
#line 307 "SB04ND.f"
		    istep = 1;
#line 308 "SB04ND.f"
		} else {
#line 309 "SB04ND.f"
		    istep = 2;
#line 310 "SB04ND.f"
		}
#line 311 "SB04ND.f"
	    }

#line 313 "SB04ND.f"
	    if (istep == 1) {
#line 314 "SB04ND.f"
		sb04nw_(abschr, ulb, n, m, &c__[c_offset], ldc, &i__, &b[
			b_offset], ldb, &dwork[jwork], (ftnlen)1, (ftnlen)1);
#line 316 "SB04ND.f"
		sb04ny_("R", ula, n, &a[a_offset], lda, &b[i__ + i__ * b_dim1]
			, &dwork[jwork], &tol1, &iwork[1], &dwork[1], &ldw, 
			info, (ftnlen)1, (ftnlen)1);
#line 318 "SB04ND.f"
		if (*info == 1) {
#line 318 "SB04ND.f"
		    return 0;
#line 318 "SB04ND.f"
		}
#line 320 "SB04ND.f"
		dcopy_(n, &dwork[jwork], &c__1, &c__[i__ * c_dim1 + 1], &c__1)
			;
#line 321 "SB04ND.f"
	    } else {
#line 322 "SB04ND.f"
		ipincr = i__ + incr;
#line 323 "SB04ND.f"
		sb04nv_(abschr, ulb, n, m, &c__[c_offset], ldc, &ipincr, &b[
			b_offset], ldb, &dwork[jwork], (ftnlen)1, (ftnlen)1);
#line 325 "SB04ND.f"
		sb04nx_("R", ula, n, &a[a_offset], lda, &b[ipincr + ipincr * 
			b_dim1], &b[ipincr + 1 + ipincr * b_dim1], &b[ipincr 
			+ (ipincr + 1) * b_dim1], &b[ipincr + 1 + (ipincr + 1)
			 * b_dim1], &dwork[jwork], &tol1, &iwork[1], &dwork[1]
			, &ldw, info, (ftnlen)1, (ftnlen)1);
#line 329 "SB04ND.f"
		if (*info == 1) {
#line 329 "SB04ND.f"
		    return 0;
#line 329 "SB04ND.f"
		}
#line 331 "SB04ND.f"
		dcopy_(n, &dwork[jwork], &c__2, &c__[ipincr * c_dim1 + 1], &
			c__1);
#line 332 "SB04ND.f"
		dcopy_(n, &dwork[jwork + 1], &c__2, &c__[(ipincr + 1) * 
			c_dim1 + 1], &c__1);
#line 333 "SB04ND.f"
	    }
#line 334 "SB04ND.f"
	    i__ += fwd * istep;
#line 335 "SB04ND.f"
	    goto L20;
#line 336 "SB04ND.f"
	}
/*        END WHILE 20 */
#line 338 "SB04ND.f"
    } else {

/*        A is in Schur form: recursion on the rows of A. */

#line 342 "SB04ND.f"
	if (lula) {

/*           A is upper: backward recursion. */

#line 346 "SB04ND.f"
	    ibeg = *n;
#line 347 "SB04ND.f"
	    iend = 1;
#line 348 "SB04ND.f"
	    fwd = -1;
#line 349 "SB04ND.f"
	    incr = -1;
#line 350 "SB04ND.f"
	} else {

/*           A is lower: forward recursion. */

#line 354 "SB04ND.f"
	    ibeg = 1;
#line 355 "SB04ND.f"
	    iend = *n;
#line 356 "SB04ND.f"
	    fwd = 1;
#line 357 "SB04ND.f"
	    incr = 0;
#line 358 "SB04ND.f"
	}
#line 359 "SB04ND.f"
	i__ = ibeg;
/*        WHILE ( ( IEND - I ) * FWD .GE. 0 ) DO */
#line 361 "SB04ND.f"
L40:
#line 361 "SB04ND.f"
	if ((iend - i__) * fwd >= 0) {

/*           Test for 1-by-1 or 2-by-2 diagonal block in the Schur */
/*           form. */

#line 366 "SB04ND.f"
	    if (i__ == iend) {
#line 367 "SB04ND.f"
		istep = 1;
#line 368 "SB04ND.f"
	    } else {
#line 369 "SB04ND.f"
		if (a[i__ + (i__ + fwd) * a_dim1] == 0.) {
#line 370 "SB04ND.f"
		    istep = 1;
#line 371 "SB04ND.f"
		} else {
#line 372 "SB04ND.f"
		    istep = 2;
#line 373 "SB04ND.f"
		}
#line 374 "SB04ND.f"
	    }

#line 376 "SB04ND.f"
	    if (istep == 1) {
#line 377 "SB04ND.f"
		sb04nw_(abschr, ula, n, m, &c__[c_offset], ldc, &i__, &a[
			a_offset], lda, &dwork[jwork], (ftnlen)1, (ftnlen)1);
#line 379 "SB04ND.f"
		sb04ny_("C", ulb, m, &b[b_offset], ldb, &a[i__ + i__ * a_dim1]
			, &dwork[jwork], &tol1, &iwork[1], &dwork[1], &ldw, 
			info, (ftnlen)1, (ftnlen)1);
#line 381 "SB04ND.f"
		if (*info == 1) {
#line 381 "SB04ND.f"
		    return 0;
#line 381 "SB04ND.f"
		}
#line 383 "SB04ND.f"
		dcopy_(m, &dwork[jwork], &c__1, &c__[i__ + c_dim1], ldc);
#line 384 "SB04ND.f"
	    } else {
#line 385 "SB04ND.f"
		ipincr = i__ + incr;
#line 386 "SB04ND.f"
		sb04nv_(abschr, ula, n, m, &c__[c_offset], ldc, &ipincr, &a[
			a_offset], lda, &dwork[jwork], (ftnlen)1, (ftnlen)1);
#line 388 "SB04ND.f"
		sb04nx_("C", ulb, m, &b[b_offset], ldb, &a[ipincr + ipincr * 
			a_dim1], &a[ipincr + 1 + ipincr * a_dim1], &a[ipincr 
			+ (ipincr + 1) * a_dim1], &a[ipincr + 1 + (ipincr + 1)
			 * a_dim1], &dwork[jwork], &tol1, &iwork[1], &dwork[1]
			, &ldw, info, (ftnlen)1, (ftnlen)1);
#line 392 "SB04ND.f"
		if (*info == 1) {
#line 392 "SB04ND.f"
		    return 0;
#line 392 "SB04ND.f"
		}
#line 394 "SB04ND.f"
		dcopy_(m, &dwork[jwork], &c__2, &c__[ipincr + c_dim1], ldc);
#line 395 "SB04ND.f"
		dcopy_(m, &dwork[jwork + 1], &c__2, &c__[ipincr + 1 + c_dim1],
			 ldc);
#line 396 "SB04ND.f"
	    }
#line 397 "SB04ND.f"
	    i__ += fwd * istep;
#line 398 "SB04ND.f"
	    goto L40;
#line 399 "SB04ND.f"
	}
/*        END WHILE 40 */
#line 401 "SB04ND.f"
    }

#line 403 "SB04ND.f"
    return 0;
/* *** Last line of SB04ND *** */
} /* sb04nd_ */

