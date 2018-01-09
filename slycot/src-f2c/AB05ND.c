#line 1 "AB05ND.f"
/* AB05ND.f -- translated by f2c (version 20100827).
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

#line 1 "AB05ND.f"
/* Table of constant values */

static doublereal c_b6 = 0.;
static doublereal c_b7 = 1.;
static integer c__1 = 1;

/* Subroutine */ int ab05nd_(char *over, integer *n1, integer *m1, integer *
	p1, integer *n2, doublereal *alpha, doublereal *a1, integer *lda1, 
	doublereal *b1, integer *ldb1, doublereal *c1, integer *ldc1, 
	doublereal *d1, integer *ldd1, doublereal *a2, integer *lda2, 
	doublereal *b2, integer *ldb2, doublereal *c2, integer *ldc2, 
	doublereal *d2, integer *ldd2, integer *n, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, integer *iwork, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen over_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, a1_dim1, a1_offset, a2_dim1, a2_offset, b_dim1, 
	    b_offset, b1_dim1, b1_offset, b2_dim1, b2_offset, c_dim1, 
	    c_offset, c1_dim1, c1_offset, c2_dim1, c2_offset, d_dim1, 
	    d_offset, d1_dim1, d1_offset, d2_dim1, d2_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, ldw, ldwm1;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static logical lover;
    extern /* Subroutine */ int dgetrf_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), dlacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), dlaset_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen), dgetrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
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

/*     To obtain the state-space model (A,B,C,D) for the feedback */
/*     inter-connection of two systems, each given in state-space form. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     OVER    CHARACTER*1 */
/*             Indicates whether the user wishes to overlap pairs of */
/*             arrays, as follows: */
/*             = 'N':  Do not overlap; */
/*             = 'O':  Overlap pairs of arrays: A1 and A, B1 and B, */
/*                     C1 and C, and D1 and D, i.e. the same name is */
/*                     effectively used for each pair (for all pairs) */
/*                     in the routine call.  In this case, setting */
/*                     LDA1 = LDA, LDB1 = LDB, LDC1 = LDC, and LDD1 = LDD */
/*                     will give maximum efficiency. */

/*     Input/Output Parameters */

/*     N1      (input) INTEGER */
/*             The number of state variables in the first system, i.e. */
/*             the order of the matrix A1.  N1 >= 0. */

/*     M1      (input) INTEGER */
/*             The number of input variables for the first system and the */
/*             number of output variables from the second system. */
/*             M1 >= 0. */

/*     P1      (input) INTEGER */
/*             The number of output variables from the first system and */
/*             the number of input variables for the second system. */
/*             P1 >= 0. */

/*     N2      (input) INTEGER */
/*             The number of state variables in the second system, i.e. */
/*             the order of the matrix A2.  N2 >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             A coefficient multiplying the transfer-function matrix */
/*             (or the output equation) of the second system. */
/*             ALPHA = +1 corresponds to positive feedback, and */
/*             ALPHA = -1 corresponds to negative feedback. */

/*     A1      (input) DOUBLE PRECISION array, dimension (LDA1,N1) */
/*             The leading N1-by-N1 part of this array must contain the */
/*             state transition matrix A1 for the first system. */

/*     LDA1    INTEGER */
/*             The leading dimension of array A1.  LDA1 >= MAX(1,N1). */

/*     B1      (input) DOUBLE PRECISION array, dimension (LDB1,M1) */
/*             The leading N1-by-M1 part of this array must contain the */
/*             input/state matrix B1 for the first system. */

/*     LDB1    INTEGER */
/*             The leading dimension of array B1.  LDB1 >= MAX(1,N1). */

/*     C1      (input) DOUBLE PRECISION array, dimension (LDC1,N1) */
/*             The leading P1-by-N1 part of this array must contain the */
/*             state/output matrix C1 for the first system. */

/*     LDC1    INTEGER */
/*             The leading dimension of array C1. */
/*             LDC1 >= MAX(1,P1) if N1 > 0. */
/*             LDC1 >= 1 if N1 = 0. */

/*     D1      (input) DOUBLE PRECISION array, dimension (LDD1,M1) */
/*             The leading P1-by-M1 part of this array must contain the */
/*             input/output matrix D1 for the first system. */

/*     LDD1    INTEGER */
/*             The leading dimension of array D1.  LDD1 >= MAX(1,P1). */

/*     A2      (input) DOUBLE PRECISION array, dimension (LDA2,N2) */
/*             The leading N2-by-N2 part of this array must contain the */
/*             state transition matrix A2 for the second system. */

/*     LDA2    INTEGER */
/*             The leading dimension of array A2.  LDA2 >= MAX(1,N2). */

/*     B2      (input) DOUBLE PRECISION array, dimension (LDB2,P1) */
/*             The leading N2-by-P1 part of this array must contain the */
/*             input/state matrix B2 for the second system. */

/*     LDB2    INTEGER */
/*             The leading dimension of array B2.  LDB2 >= MAX(1,N2). */

/*     C2      (input) DOUBLE PRECISION array, dimension (LDC2,N2) */
/*             The leading M1-by-N2 part of this array must contain the */
/*             state/output matrix C2 for the second system. */

/*     LDC2    INTEGER */
/*             The leading dimension of array C2. */
/*             LDC2 >= MAX(1,M1) if N2 > 0. */
/*             LDC2 >= 1 if N2 = 0. */

/*     D2      (input) DOUBLE PRECISION array, dimension (LDD2,P1) */
/*             The leading M1-by-P1 part of this array must contain the */
/*             input/output matrix D2 for the second system. */

/*     LDD2    INTEGER */
/*             The leading dimension of array D2.  LDD2 >= MAX(1,M1). */

/*     N       (output) INTEGER */
/*             The number of state variables (N1 + N2) in the connected */
/*             system, i.e. the order of the matrix A, the number of rows */
/*             of B and the number of columns of C. */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N1+N2) */
/*             The leading N-by-N part of this array contains the state */
/*             transition matrix A for the connected system. */
/*             The array A can overlap A1 if OVER = 'O'. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N1+N2). */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,M1) */
/*             The leading N-by-M1 part of this array contains the */
/*             input/state matrix B for the connected system. */
/*             The array B can overlap B1 if OVER = 'O'. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N1+N2). */

/*     C       (output) DOUBLE PRECISION array, dimension (LDC,N1+N2) */
/*             The leading P1-by-N part of this array contains the */
/*             state/output matrix C for the connected system. */
/*             The array C can overlap C1 if OVER = 'O'. */

/*     LDC     INTEGER */
/*             The leading dimension of array C. */
/*             LDC >= MAX(1,P1) if N1+N2 > 0. */
/*             LDC >= 1 if N1+N2 = 0. */

/*     D       (output) DOUBLE PRECISION array, dimension (LDD,M1) */
/*             The leading P1-by-M1 part of this array contains the */
/*             input/output matrix D for the connected system. */
/*             The array D can overlap D1 if OVER = 'O'. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P1). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (P1) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.             If OVER = 'N', */
/*             LDWORK >= MAX(1, P1*P1, M1*M1, N1*P1), and if OVER = 'O', */
/*             LDWORK >= MAX(1, N1*P1 + MAX( P1*P1, M1*M1, N1*P1) ), */
/*                                                        if M1 <= N*N2; */
/*             LDWORK >= MAX(1, N1*P1 + MAX( P1*P1, M1*(M1+1), N1*P1) ), */
/*                                                        if M1 >  N*N2. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */
/*             > 0:  if INFO = i, 1 <= i <= P1, the system is not */
/*                   completely controllable. That is, the matrix */
/*                   (I + ALPHA*D1*D2) is exactly singular (the element */
/*                   U(i,i) of the upper triangular factor of LU */
/*                   factorization is exactly zero), possibly due to */
/*                   rounding errors. */

/*     METHOD */

/*     After feedback inter-connection of the two systems, */

/*     X1'     = A1*X1 + B1*U1 */
/*     Y1      = C1*X1 + D1*U1 */

/*     X2'     = A2*X2 + B2*U2 */
/*     Y2      = C2*X2 + D2*U2 */

/*     (where  '  denotes differentiation with respect to time) */

/*     the following state-space model will be obtained: */

/*     X'      = A*X  +  B*U */
/*     Y       = C*X  +  D*U */

/*     where       U = U1 + alpha*Y2,    X  =  ( X1 ), */
/*                 Y = Y1 = U2,                ( X2 ) */

/*     matrix  A  has the form */

/*     ( A1  -  alpha*B1*E12*D2*C1       -  alpha*B1*E12*C2    ), */
/*     (        B2*E21*C1            A2  -  alpha*B2*E21*D1*C2 ) */

/*     matrix  B  has the form */

/*     (  B1*E12    ), */
/*     (  B2*E21*D1 ) */

/*     matrix  C  has the form */

/*     (  E21*C1     -  alpha*E21*D1*C2 ), */

/*     matrix D  has the form */

/*     (  E21*D1 ), */

/*     E21  =  ( I + alpha*D1*D2 )-INVERSE and */
/*     E12  =  ( I + alpha*D2*D1 )-INVERSE = I - alpha*D2*E21*D1. */

/*     Taking N1 = 0 and/or N2 = 0 on the routine call will solve the */
/*     constant plant and/or constant feedback cases. */

/*     REFERENCES */

/*     None */

/*     NUMERICAL ASPECTS */

/*     None */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996. */
/*     Supersedes Release 2.0 routine AB05BD by C.J.Benson, Kingston */
/*     Polytechnic, United Kingdom, January 1982. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, July 2003, */
/*     Feb. 2004. */

/*     KEYWORDS */

/*     Continuous-time system, multivariable system, state-space model, */
/*     state-space representation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 294 "AB05ND.f"
    /* Parameter adjustments */
#line 294 "AB05ND.f"
    a1_dim1 = *lda1;
#line 294 "AB05ND.f"
    a1_offset = 1 + a1_dim1;
#line 294 "AB05ND.f"
    a1 -= a1_offset;
#line 294 "AB05ND.f"
    b1_dim1 = *ldb1;
#line 294 "AB05ND.f"
    b1_offset = 1 + b1_dim1;
#line 294 "AB05ND.f"
    b1 -= b1_offset;
#line 294 "AB05ND.f"
    c1_dim1 = *ldc1;
#line 294 "AB05ND.f"
    c1_offset = 1 + c1_dim1;
#line 294 "AB05ND.f"
    c1 -= c1_offset;
#line 294 "AB05ND.f"
    d1_dim1 = *ldd1;
#line 294 "AB05ND.f"
    d1_offset = 1 + d1_dim1;
#line 294 "AB05ND.f"
    d1 -= d1_offset;
#line 294 "AB05ND.f"
    a2_dim1 = *lda2;
#line 294 "AB05ND.f"
    a2_offset = 1 + a2_dim1;
#line 294 "AB05ND.f"
    a2 -= a2_offset;
#line 294 "AB05ND.f"
    b2_dim1 = *ldb2;
#line 294 "AB05ND.f"
    b2_offset = 1 + b2_dim1;
#line 294 "AB05ND.f"
    b2 -= b2_offset;
#line 294 "AB05ND.f"
    c2_dim1 = *ldc2;
#line 294 "AB05ND.f"
    c2_offset = 1 + c2_dim1;
#line 294 "AB05ND.f"
    c2 -= c2_offset;
#line 294 "AB05ND.f"
    d2_dim1 = *ldd2;
#line 294 "AB05ND.f"
    d2_offset = 1 + d2_dim1;
#line 294 "AB05ND.f"
    d2 -= d2_offset;
#line 294 "AB05ND.f"
    a_dim1 = *lda;
#line 294 "AB05ND.f"
    a_offset = 1 + a_dim1;
#line 294 "AB05ND.f"
    a -= a_offset;
#line 294 "AB05ND.f"
    b_dim1 = *ldb;
#line 294 "AB05ND.f"
    b_offset = 1 + b_dim1;
#line 294 "AB05ND.f"
    b -= b_offset;
#line 294 "AB05ND.f"
    c_dim1 = *ldc;
#line 294 "AB05ND.f"
    c_offset = 1 + c_dim1;
#line 294 "AB05ND.f"
    c__ -= c_offset;
#line 294 "AB05ND.f"
    d_dim1 = *ldd;
#line 294 "AB05ND.f"
    d_offset = 1 + d_dim1;
#line 294 "AB05ND.f"
    d__ -= d_offset;
#line 294 "AB05ND.f"
    --iwork;
#line 294 "AB05ND.f"
    --dwork;
#line 294 "AB05ND.f"

#line 294 "AB05ND.f"
    /* Function Body */
#line 294 "AB05ND.f"
    lover = lsame_(over, "O", (ftnlen)1, (ftnlen)1);
#line 295 "AB05ND.f"
    ldwm1 = max(1,*m1);
#line 296 "AB05ND.f"
    *n = *n1 + *n2;
#line 297 "AB05ND.f"
    *info = 0;

/*     Test the input scalar arguments. */

#line 301 "AB05ND.f"
    if (! lover && ! lsame_(over, "N", (ftnlen)1, (ftnlen)1)) {
#line 302 "AB05ND.f"
	*info = -1;
#line 303 "AB05ND.f"
    } else if (*n1 < 0) {
#line 304 "AB05ND.f"
	*info = -2;
#line 305 "AB05ND.f"
    } else if (*m1 < 0) {
#line 306 "AB05ND.f"
	*info = -3;
#line 307 "AB05ND.f"
    } else if (*p1 < 0) {
#line 308 "AB05ND.f"
	*info = -4;
#line 309 "AB05ND.f"
    } else if (*n2 < 0) {
#line 310 "AB05ND.f"
	*info = -5;
#line 311 "AB05ND.f"
    } else if (*lda1 < max(1,*n1)) {
#line 312 "AB05ND.f"
	*info = -8;
#line 313 "AB05ND.f"
    } else if (*ldb1 < max(1,*n1)) {
#line 314 "AB05ND.f"
	*info = -10;
#line 315 "AB05ND.f"
    } else if (*n1 > 0 && *ldc1 < max(1,*p1) || *n1 == 0 && *ldc1 < 1) {
#line 317 "AB05ND.f"
	*info = -12;
#line 318 "AB05ND.f"
    } else if (*ldd1 < max(1,*p1)) {
#line 319 "AB05ND.f"
	*info = -14;
#line 320 "AB05ND.f"
    } else if (*lda2 < max(1,*n2)) {
#line 321 "AB05ND.f"
	*info = -16;
#line 322 "AB05ND.f"
    } else if (*ldb2 < max(1,*n2)) {
#line 323 "AB05ND.f"
	*info = -18;
#line 324 "AB05ND.f"
    } else if (*n2 > 0 && *ldc2 < ldwm1 || *n2 == 0 && *ldc2 < 1) {
#line 326 "AB05ND.f"
	*info = -20;
#line 327 "AB05ND.f"
    } else if (*ldd2 < ldwm1) {
#line 328 "AB05ND.f"
	*info = -22;
#line 329 "AB05ND.f"
    } else if (*lda < max(1,*n)) {
#line 330 "AB05ND.f"
	*info = -25;
#line 331 "AB05ND.f"
    } else if (*ldb < max(1,*n)) {
#line 332 "AB05ND.f"
	*info = -27;
#line 333 "AB05ND.f"
    } else if (*n > 0 && *ldc < max(1,*p1) || *n == 0 && *ldc < 1) {
#line 335 "AB05ND.f"
	*info = -29;
#line 336 "AB05ND.f"
    } else if (*ldd < max(1,*p1)) {
#line 337 "AB05ND.f"
	*info = -31;
#line 338 "AB05ND.f"
    } else {
/* Computing MAX */
#line 339 "AB05ND.f"
	i__1 = *p1 * *p1, i__2 = *m1 * *m1, i__1 = max(i__1,i__2), i__2 = *n1 
		* *p1;
#line 339 "AB05ND.f"
	ldw = max(i__1,i__2);
#line 340 "AB05ND.f"
	if (lover) {
#line 341 "AB05ND.f"
	    if (*m1 > *n * *n2) {
/* Computing MAX */
#line 341 "AB05ND.f"
		i__1 = ldw, i__2 = *m1 * (*m1 + 1);
#line 341 "AB05ND.f"
		ldw = max(i__1,i__2);
#line 341 "AB05ND.f"
	    }
#line 343 "AB05ND.f"
	    ldw = *n1 * *p1 + ldw;
#line 344 "AB05ND.f"
	}
#line 345 "AB05ND.f"
	if (*ldwork < max(1,ldw)) {
#line 345 "AB05ND.f"
	    *info = -34;
#line 345 "AB05ND.f"
	}
#line 347 "AB05ND.f"
    }

#line 349 "AB05ND.f"
    if (*info != 0) {

/*        Error return. */

#line 353 "AB05ND.f"
	i__1 = -(*info);
#line 353 "AB05ND.f"
	xerbla_("AB05ND", &i__1, (ftnlen)6);
#line 354 "AB05ND.f"
	return 0;
#line 355 "AB05ND.f"
    }

/*     Quick return if possible. */

/* Computing MAX */
#line 359 "AB05ND.f"
    i__1 = *n, i__2 = min(*m1,*p1);
#line 359 "AB05ND.f"
    if (max(i__1,i__2) == 0) {
#line 359 "AB05ND.f"
	return 0;
#line 359 "AB05ND.f"
    }

#line 362 "AB05ND.f"
    if (*p1 > 0) {

/*        Form  ( I  +  alpha * D1 * D2 ). */

#line 366 "AB05ND.f"
	dlaset_("F", p1, p1, &c_b6, &c_b7, &dwork[1], p1, (ftnlen)1);
#line 367 "AB05ND.f"
	dgemm_("No transpose", "No transpose", p1, p1, m1, alpha, &d1[
		d1_offset], ldd1, &d2[d2_offset], ldd2, &c_b7, &dwork[1], p1, 
		(ftnlen)12, (ftnlen)12);

/*        Factorize this matrix. */

#line 372 "AB05ND.f"
	dgetrf_(p1, p1, &dwork[1], p1, &iwork[1], info);

#line 374 "AB05ND.f"
	if (*info != 0) {
#line 374 "AB05ND.f"
	    return 0;
#line 374 "AB05ND.f"
	}

/*        Form  E21 * D1. */

#line 379 "AB05ND.f"
	if (lover && *ldd1 <= *ldd) {
#line 380 "AB05ND.f"
	    if (*ldd1 < *ldd) {

#line 382 "AB05ND.f"
		for (j = *m1; j >= 1; --j) {
#line 383 "AB05ND.f"
		    for (i__ = *p1; i__ >= 1; --i__) {
#line 384 "AB05ND.f"
			d__[i__ + j * d_dim1] = d1[i__ + j * d1_dim1];
#line 385 "AB05ND.f"
/* L10: */
#line 385 "AB05ND.f"
		    }
#line 386 "AB05ND.f"
/* L20: */
#line 386 "AB05ND.f"
		}

#line 388 "AB05ND.f"
	    }
#line 389 "AB05ND.f"
	} else {
#line 390 "AB05ND.f"
	    dlacpy_("F", p1, m1, &d1[d1_offset], ldd1, &d__[d_offset], ldd, (
		    ftnlen)1);
#line 391 "AB05ND.f"
	}

#line 393 "AB05ND.f"
	dgetrs_("No transpose", p1, m1, &dwork[1], p1, &iwork[1], &d__[
		d_offset], ldd, info, (ftnlen)12);
#line 395 "AB05ND.f"
	if (*n1 > 0) {

/*           Form  E21 * C1. */

#line 399 "AB05ND.f"
	    if (lover) {

/*              First save  C1. */

#line 403 "AB05ND.f"
		ldw = ldw - *p1 * *n1 + 1;
#line 404 "AB05ND.f"
		dlacpy_("F", p1, n1, &c1[c1_offset], ldc1, &dwork[ldw], p1, (
			ftnlen)1);

#line 406 "AB05ND.f"
		if (*ldc1 != *ldc) {
#line 406 "AB05ND.f"
		    dlacpy_("F", p1, n1, &dwork[ldw], p1, &c__[c_offset], ldc,
			     (ftnlen)1);
#line 406 "AB05ND.f"
		}
#line 408 "AB05ND.f"
	    } else {
#line 409 "AB05ND.f"
		dlacpy_("F", p1, n1, &c1[c1_offset], ldc1, &c__[c_offset], 
			ldc, (ftnlen)1);
#line 410 "AB05ND.f"
	    }

#line 412 "AB05ND.f"
	    dgetrs_("No transpose", p1, n1, &dwork[1], p1, &iwork[1], &c__[
		    c_offset], ldc, info, (ftnlen)12);
#line 414 "AB05ND.f"
	}

/*        Form  E12  =  I  -  alpha * D2 * ( E21 * D1 ). */

#line 418 "AB05ND.f"
	dlaset_("F", m1, m1, &c_b6, &c_b7, &dwork[1], &ldwm1, (ftnlen)1);
#line 419 "AB05ND.f"
	d__1 = -(*alpha);
#line 419 "AB05ND.f"
	dgemm_("No transpose", "No transpose", m1, m1, p1, &d__1, &d2[
		d2_offset], ldd2, &d__[d_offset], ldd, &c_b7, &dwork[1], &
		ldwm1, (ftnlen)12, (ftnlen)12);

#line 422 "AB05ND.f"
    } else {
#line 423 "AB05ND.f"
	dlaset_("F", m1, m1, &c_b6, &c_b7, &dwork[1], &ldwm1, (ftnlen)1);
#line 424 "AB05ND.f"
    }

#line 426 "AB05ND.f"
    if (lover && *lda1 <= *lda) {
#line 427 "AB05ND.f"
	if (*lda1 < *lda) {

#line 429 "AB05ND.f"
	    for (j = *n1; j >= 1; --j) {
#line 430 "AB05ND.f"
		for (i__ = *n1; i__ >= 1; --i__) {
#line 431 "AB05ND.f"
		    a[i__ + j * a_dim1] = a1[i__ + j * a1_dim1];
#line 432 "AB05ND.f"
/* L30: */
#line 432 "AB05ND.f"
		}
#line 433 "AB05ND.f"
/* L40: */
#line 433 "AB05ND.f"
	    }

#line 435 "AB05ND.f"
	}
#line 436 "AB05ND.f"
    } else {
#line 437 "AB05ND.f"
	dlacpy_("F", n1, n1, &a1[a1_offset], lda1, &a[a_offset], lda, (ftnlen)
		1);
#line 438 "AB05ND.f"
    }

#line 440 "AB05ND.f"
    if (*n1 > 0 && *m1 > 0) {

/*        Form  B1 * E12. */

#line 444 "AB05ND.f"
	if (lover) {

/*           Use the blocks (1,2) and (2,2) of A as workspace. */

#line 448 "AB05ND.f"
	    if (*n1 * *m1 <= *n * *n2) {

/*              Use BLAS 3 code. */

#line 452 "AB05ND.f"
		dlacpy_("F", n1, m1, &b1[b1_offset], ldb1, &a[(*n1 + 1) * 
			a_dim1 + 1], n1, (ftnlen)1);
#line 453 "AB05ND.f"
		dgemm_("No transpose", "No transpose", n1, m1, m1, &c_b7, &a[(
			*n1 + 1) * a_dim1 + 1], n1, &dwork[1], &ldwm1, &c_b6, 
			&b[b_offset], ldb, (ftnlen)12, (ftnlen)12);
#line 456 "AB05ND.f"
	    } else if (*ldb1 < *ldb) {

#line 458 "AB05ND.f"
		for (j = *m1; j >= 1; --j) {
#line 459 "AB05ND.f"
		    for (i__ = *n1; i__ >= 1; --i__) {
#line 460 "AB05ND.f"
			b[i__ + j * b_dim1] = b1[i__ + j * b1_dim1];
#line 461 "AB05ND.f"
/* L50: */
#line 461 "AB05ND.f"
		    }
#line 462 "AB05ND.f"
/* L60: */
#line 462 "AB05ND.f"
		}

#line 464 "AB05ND.f"
		if (*m1 <= *n * *n2) {

/*                 Use BLAS 2 code. */

#line 468 "AB05ND.f"
		    i__1 = *n1;
#line 468 "AB05ND.f"
		    for (j = 1; j <= i__1; ++j) {
#line 469 "AB05ND.f"
			dcopy_(m1, &b[j + b_dim1], ldb, &a[(*n1 + 1) * a_dim1 
				+ 1], &c__1);
#line 470 "AB05ND.f"
			dgemv_("Transpose", m1, m1, &c_b7, &dwork[1], &ldwm1, 
				&a[(*n1 + 1) * a_dim1 + 1], &c__1, &c_b6, &b[
				j + b_dim1], ldb, (ftnlen)9);
#line 472 "AB05ND.f"
/* L70: */
#line 472 "AB05ND.f"
		    }

#line 474 "AB05ND.f"
		} else {

/*                 Use additional workspace. */

#line 478 "AB05ND.f"
		    i__1 = *n1;
#line 478 "AB05ND.f"
		    for (j = 1; j <= i__1; ++j) {
#line 479 "AB05ND.f"
			dcopy_(m1, &b[j + b_dim1], ldb, &dwork[*m1 * *m1 + 1],
				 &c__1);
#line 480 "AB05ND.f"
			dgemv_("Transpose", m1, m1, &c_b7, &dwork[1], &ldwm1, 
				&dwork[*m1 * *m1 + 1], &c__1, &c_b6, &b[j + 
				b_dim1], ldb, (ftnlen)9);
#line 482 "AB05ND.f"
/* L80: */
#line 482 "AB05ND.f"
		    }

#line 484 "AB05ND.f"
		}

#line 486 "AB05ND.f"
	    } else if (*m1 <= *n * *n2) {

/*              Use BLAS 2 code. */

#line 490 "AB05ND.f"
		i__1 = *n1;
#line 490 "AB05ND.f"
		for (j = 1; j <= i__1; ++j) {
#line 491 "AB05ND.f"
		    dcopy_(m1, &b1[j + b1_dim1], ldb1, &a[(*n1 + 1) * a_dim1 
			    + 1], &c__1);
#line 492 "AB05ND.f"
		    dgemv_("Transpose", m1, m1, &c_b7, &dwork[1], &ldwm1, &a[(
			    *n1 + 1) * a_dim1 + 1], &c__1, &c_b6, &b[j + 
			    b_dim1], ldb, (ftnlen)9);
#line 494 "AB05ND.f"
/* L90: */
#line 494 "AB05ND.f"
		}

#line 496 "AB05ND.f"
	    } else {

/*              Use additional workspace. */

#line 500 "AB05ND.f"
		i__1 = *n1;
#line 500 "AB05ND.f"
		for (j = 1; j <= i__1; ++j) {
#line 501 "AB05ND.f"
		    dcopy_(m1, &b1[j + b1_dim1], ldb1, &dwork[*m1 * *m1 + 1], 
			    &c__1);
#line 502 "AB05ND.f"
		    dgemv_("Transpose", m1, m1, &c_b7, &dwork[1], &ldwm1, &
			    dwork[*m1 * *m1 + 1], &c__1, &c_b6, &b[j + b_dim1]
			    , ldb, (ftnlen)9);
#line 504 "AB05ND.f"
/* L100: */
#line 504 "AB05ND.f"
		}

#line 506 "AB05ND.f"
	    }
#line 507 "AB05ND.f"
	} else {
#line 508 "AB05ND.f"
	    dgemm_("No transpose", "No transpose", n1, m1, m1, &c_b7, &b1[
		    b1_offset], ldb1, &dwork[1], &ldwm1, &c_b6, &b[b_offset], 
		    ldb, (ftnlen)12, (ftnlen)12);
#line 510 "AB05ND.f"
	}
#line 511 "AB05ND.f"
    }

#line 513 "AB05ND.f"
    if (*n2 > 0) {

/*        Complete matrices  B  and  C. */

#line 517 "AB05ND.f"
	if (*p1 > 0) {
#line 518 "AB05ND.f"
	    dgemm_("No transpose", "No transpose", n2, m1, p1, &c_b7, &b2[
		    b2_offset], ldb2, &d__[d_offset], ldd, &c_b6, &b[*n1 + 1 
		    + b_dim1], ldb, (ftnlen)12, (ftnlen)12);
#line 520 "AB05ND.f"
	    d__1 = -(*alpha);
#line 520 "AB05ND.f"
	    dgemm_("No transpose", "No transpose", p1, n2, m1, &d__1, &d__[
		    d_offset], ldd, &c2[c2_offset], ldc2, &c_b6, &c__[(*n1 + 
		    1) * c_dim1 + 1], ldc, (ftnlen)12, (ftnlen)12);
#line 523 "AB05ND.f"
	} else if (*m1 > 0) {
#line 524 "AB05ND.f"
	    dlaset_("F", n2, m1, &c_b6, &c_b6, &b[*n1 + 1 + b_dim1], ldb, (
		    ftnlen)1);
#line 525 "AB05ND.f"
	}
#line 526 "AB05ND.f"
    }

#line 528 "AB05ND.f"
    if (*n1 > 0 && *p1 > 0) {

/*        Form upper left quadrant of  A. */

#line 532 "AB05ND.f"
	d__1 = -(*alpha);
#line 532 "AB05ND.f"
	dgemm_("No transpose", "No transpose", n1, p1, m1, &d__1, &b[b_offset]
		, ldb, &d2[d2_offset], ldd2, &c_b6, &dwork[1], n1, (ftnlen)12,
		 (ftnlen)12);

#line 535 "AB05ND.f"
	if (lover) {
#line 536 "AB05ND.f"
	    dgemm_("No transpose", "No transpose", n1, n1, p1, &c_b7, &dwork[
		    1], n1, &dwork[ldw], p1, &c_b7, &a[a_offset], lda, (
		    ftnlen)12, (ftnlen)12);
#line 538 "AB05ND.f"
	} else {
#line 539 "AB05ND.f"
	    dgemm_("No transpose", "No transpose", n1, n1, p1, &c_b7, &dwork[
		    1], n1, &c1[c1_offset], ldc1, &c_b7, &a[a_offset], lda, (
		    ftnlen)12, (ftnlen)12);
#line 541 "AB05ND.f"
	}
#line 542 "AB05ND.f"
    }

#line 544 "AB05ND.f"
    if (*n2 > 0) {

/*        Form lower right quadrant of  A. */

#line 548 "AB05ND.f"
	dlacpy_("F", n2, n2, &a2[a2_offset], lda2, &a[*n1 + 1 + (*n1 + 1) * 
		a_dim1], lda, (ftnlen)1);
#line 549 "AB05ND.f"
	if (*m1 > 0) {
#line 549 "AB05ND.f"
	    d__1 = -(*alpha);
#line 549 "AB05ND.f"
	    dgemm_("No transpose", "No transpose", n2, n2, m1, &d__1, &b[*n1 
		    + 1 + b_dim1], ldb, &c2[c2_offset], ldc2, &c_b7, &a[*n1 + 
		    1 + (*n1 + 1) * a_dim1], lda, (ftnlen)12, (ftnlen)12);
#line 549 "AB05ND.f"
	}

/*        Complete the matrix  A. */

#line 556 "AB05ND.f"
	dgemm_("No transpose", "No transpose", n2, n1, p1, &c_b7, &b2[
		b2_offset], ldb2, &c__[c_offset], ldc, &c_b6, &a[*n1 + 1 + 
		a_dim1], lda, (ftnlen)12, (ftnlen)12);
#line 558 "AB05ND.f"
	d__1 = -(*alpha);
#line 558 "AB05ND.f"
	dgemm_("No transpose", "No transpose", n1, n2, m1, &d__1, &b[b_offset]
		, ldb, &c2[c2_offset], ldc2, &c_b6, &a[(*n1 + 1) * a_dim1 + 1]
		, lda, (ftnlen)12, (ftnlen)12);
#line 560 "AB05ND.f"
    }

#line 562 "AB05ND.f"
    return 0;
/* *** Last line of AB05ND *** */
} /* ab05nd_ */
