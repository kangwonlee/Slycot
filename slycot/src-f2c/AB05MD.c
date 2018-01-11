#line 1 "AB05MD.f"
/* AB05MD.f -- translated by f2c (version 20100827).
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

#line 1 "AB05MD.f"
/* Table of constant values */

static doublereal c_b16 = 0.;
static doublereal c_b20 = 1.;

/* Subroutine */ int ab05md_(char *uplo, char *over, integer *n1, integer *m1,
	 integer *p1, integer *n2, integer *p2, doublereal *a1, integer *lda1,
	 doublereal *b1, integer *ldb1, doublereal *c1, integer *ldc1, 
	doublereal *d1, integer *ldd1, doublereal *a2, integer *lda2, 
	doublereal *b2, integer *ldb2, doublereal *c2, integer *ldc2, 
	doublereal *d2, integer *ldd2, integer *n, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen uplo_len, ftnlen over_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, a1_dim1, a1_offset, a2_dim1, a2_offset, b_dim1, 
	    b_offset, b1_dim1, b1_offset, b2_dim1, b2_offset, c_dim1, 
	    c_offset, c1_dim1, c1_offset, c2_dim1, c2_offset, d_dim1, 
	    d_offset, d1_dim1, d1_offset, d2_dim1, d2_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, j, i1, i2, ldwn2, ldwp1, ldwp2;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lover, luplo;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);


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

/*     To obtain the state-space model (A,B,C,D) for the cascaded */
/*     inter-connection of two systems, each given in state-space form. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPLO    CHARACTER*1 */
/*             Indicates whether the user wishes to obtain the matrix A */
/*             in the upper or lower block diagonal form, as follows: */
/*             = 'U':  Obtain A in the upper block diagonal form; */
/*             = 'L':  Obtain A in the lower block diagonal form. */

/*     OVER    CHARACTER*1 */
/*             Indicates whether the user wishes to overlap pairs of */
/*             arrays, as follows: */
/*             = 'N':  Do not overlap; */
/*             = 'O':  Overlap pairs of arrays: A1 and A, B1 and B, */
/*                     C1 and C, and D1 and D (for UPLO = 'L'), or A2 */
/*                     and A, B2 and B, C2 and C, and D2 and D (for */
/*                     UPLO = 'U'), i.e. the same name is effectively */
/*                     used for each pair (for all pairs) in the routine */
/*                     call.  In this case, setting LDA1 = LDA, */
/*                     LDB1 = LDB, LDC1 = LDC, and LDD1 = LDD, or */
/*                     LDA2 = LDA, LDB2 = LDB, LDC2 = LDC, and LDD2 = LDD */
/*                     will give maximum efficiency. */

/*     Input/Output Parameters */

/*     N1      (input) INTEGER */
/*             The number of state variables in the first system, i.e. */
/*             the order of the matrix A1.  N1 >= 0. */

/*     M1      (input) INTEGER */
/*             The number of input variables for the first system. */
/*             M1 >= 0. */

/*     P1      (input) INTEGER */
/*             The number of output variables from the first system and */
/*             the number of input variables for the second system. */
/*             P1 >= 0. */

/*     N2      (input) INTEGER */
/*             The number of state variables in the second system, i.e. */
/*             the order of the matrix A2.  N2 >= 0. */

/*     P2      (input) INTEGER */
/*             The number of output variables from the second system. */
/*             P2 >= 0. */

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
/*             The leading P2-by-N2 part of this array must contain the */
/*             state/output matrix C2 for the second system. */

/*     LDC2    INTEGER */
/*             The leading dimension of array C2. */
/*             LDC2 >= MAX(1,P2) if N2 > 0. */
/*             LDC2 >= 1 if N2 = 0. */

/*     D2      (input) DOUBLE PRECISION array, dimension (LDD2,P1) */
/*             The leading P2-by-P1 part of this array must contain the */
/*             input/output matrix D2 for the second system. */

/*     LDD2    INTEGER */
/*             The leading dimension of array D2.  LDD2 >= MAX(1,P2). */

/*     N       (output) INTEGER */
/*             The number of state variables (N1 + N2) in the resulting */
/*             system, i.e. the order of the matrix A, the number of rows */
/*             of B and the number of columns of C. */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N1+N2) */
/*             The leading N-by-N part of this array contains the state */
/*             transition matrix A for the cascaded system. */
/*             If OVER = 'O', the array A can overlap A1, if UPLO = 'L', */
/*             or A2, if UPLO = 'U'. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N1+N2). */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,M1) */
/*             The leading N-by-M1 part of this array contains the */
/*             input/state matrix B for the cascaded system. */
/*             If OVER = 'O', the array B can overlap B1, if UPLO = 'L', */
/*             or B2, if UPLO = 'U'. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N1+N2). */

/*     C       (output) DOUBLE PRECISION array, dimension (LDC,N1+N2) */
/*             The leading P2-by-N part of this array contains the */
/*             state/output matrix C for the cascaded system. */
/*             If OVER = 'O', the array C can overlap C1, if UPLO = 'L', */
/*             or C2, if UPLO = 'U'. */

/*     LDC     INTEGER */
/*             The leading dimension of array C. */
/*             LDC >= MAX(1,P2) if N1+N2 > 0. */
/*             LDC >= 1 if N1+N2 = 0. */

/*     D       (output) DOUBLE PRECISION array, dimension (LDD,M1) */
/*             The leading P2-by-M1 part of this array contains the */
/*             input/output matrix D for the cascaded system. */
/*             If OVER = 'O', the array D can overlap D1, if UPLO = 'L', */
/*             or D2, if UPLO = 'U'. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P2). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             The array DWORK is not referenced if OVER = 'N'. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 1, P1*MAX(N1, M1, N2, P2) ) if OVER = 'O'. */
/*             LDWORK >= 1 if OVER = 'N'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     After cascaded inter-connection of the two systems */

/*     X1'     = A1*X1 + B1*U */
/*     V       = C1*X1 + D1*U */

/*     X2'     = A2*X2 + B2*V */
/*     Y       = C2*X2 + D2*V */

/*     (where  '  denotes differentiation with respect to time) */

/*     the following state-space model will be obtained: */

/*     X'      = A*X + B*U */
/*     Y       = C*X + D*U */

/*     where matrix  A  has the form   ( A1     0 ), */
/*                                     ( B2*C1  A2) */

/*           matrix  B  has the form  (  B1   ), */
/*                                    ( B2*D1 ) */

/*           matrix  C  has the form  ( D2*C1  C2 ) and */

/*           matrix  D  has the form  ( D2*D1 ). */

/*     This form is returned by the routine when UPLO = 'L'.  Note that */
/*     when A1 and A2 are block lower triangular, the resulting state */
/*     matrix is also block lower triangular. */

/*     By applying a similarity transformation to the system above, */
/*     using the matrix  ( 0  I ),  where  I  is the identity matrix of */
/*                       ( J  0 ) */
/*     order  N2,  and  J  is the identity matrix of order  N1,  the */
/*     system matrices become */

/*           A = ( A2  B2*C1 ), */
/*               ( 0     A1  ) */

/*           B = ( B2*D1 ), */
/*               (  B1   ) */

/*           C = ( C2  D2*C1 ) and */

/*           D = ( D2*D1 ). */

/*     This form is returned by the routine when UPLO = 'U'.  Note that */
/*     when A1 and A2 are block upper triangular (for instance, in the */
/*     real Schur form), the resulting state matrix is also block upper */
/*     triangular. */

/*     REFERENCES */

/*     None */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires P1*(N1+M1)*(N2+P2) operations. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, and */
/*                  A. Varga, German Aerospace Research Establishment, */
/*                  Oberpfaffenhofen, Germany, Nov. 1996. */
/*     Supersedes Release 2.0 routine AB05AD by C.J.Benson, Kingston */
/*     Polytechnic, United Kingdom, January 1982. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, July 2003, */
/*     Feb. 2004. */

/*     KEYWORDS */

/*     Cascade control, continuous-time system, multivariable */
/*     system, state-space model, state-space representation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 301 "AB05MD.f"
    /* Parameter adjustments */
#line 301 "AB05MD.f"
    a1_dim1 = *lda1;
#line 301 "AB05MD.f"
    a1_offset = 1 + a1_dim1;
#line 301 "AB05MD.f"
    a1 -= a1_offset;
#line 301 "AB05MD.f"
    b1_dim1 = *ldb1;
#line 301 "AB05MD.f"
    b1_offset = 1 + b1_dim1;
#line 301 "AB05MD.f"
    b1 -= b1_offset;
#line 301 "AB05MD.f"
    c1_dim1 = *ldc1;
#line 301 "AB05MD.f"
    c1_offset = 1 + c1_dim1;
#line 301 "AB05MD.f"
    c1 -= c1_offset;
#line 301 "AB05MD.f"
    d1_dim1 = *ldd1;
#line 301 "AB05MD.f"
    d1_offset = 1 + d1_dim1;
#line 301 "AB05MD.f"
    d1 -= d1_offset;
#line 301 "AB05MD.f"
    a2_dim1 = *lda2;
#line 301 "AB05MD.f"
    a2_offset = 1 + a2_dim1;
#line 301 "AB05MD.f"
    a2 -= a2_offset;
#line 301 "AB05MD.f"
    b2_dim1 = *ldb2;
#line 301 "AB05MD.f"
    b2_offset = 1 + b2_dim1;
#line 301 "AB05MD.f"
    b2 -= b2_offset;
#line 301 "AB05MD.f"
    c2_dim1 = *ldc2;
#line 301 "AB05MD.f"
    c2_offset = 1 + c2_dim1;
#line 301 "AB05MD.f"
    c2 -= c2_offset;
#line 301 "AB05MD.f"
    d2_dim1 = *ldd2;
#line 301 "AB05MD.f"
    d2_offset = 1 + d2_dim1;
#line 301 "AB05MD.f"
    d2 -= d2_offset;
#line 301 "AB05MD.f"
    a_dim1 = *lda;
#line 301 "AB05MD.f"
    a_offset = 1 + a_dim1;
#line 301 "AB05MD.f"
    a -= a_offset;
#line 301 "AB05MD.f"
    b_dim1 = *ldb;
#line 301 "AB05MD.f"
    b_offset = 1 + b_dim1;
#line 301 "AB05MD.f"
    b -= b_offset;
#line 301 "AB05MD.f"
    c_dim1 = *ldc;
#line 301 "AB05MD.f"
    c_offset = 1 + c_dim1;
#line 301 "AB05MD.f"
    c__ -= c_offset;
#line 301 "AB05MD.f"
    d_dim1 = *ldd;
#line 301 "AB05MD.f"
    d_offset = 1 + d_dim1;
#line 301 "AB05MD.f"
    d__ -= d_offset;
#line 301 "AB05MD.f"
    --dwork;
#line 301 "AB05MD.f"

#line 301 "AB05MD.f"
    /* Function Body */
#line 301 "AB05MD.f"
    lover = lsame_(over, "O", (ftnlen)1, (ftnlen)1);
#line 302 "AB05MD.f"
    luplo = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 303 "AB05MD.f"
    *n = *n1 + *n2;
#line 304 "AB05MD.f"
    *info = 0;

/*     Test the input scalar arguments. */

#line 308 "AB05MD.f"
    if (! luplo && ! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 309 "AB05MD.f"
	*info = -1;
#line 310 "AB05MD.f"
    } else if (! lover && ! lsame_(over, "N", (ftnlen)1, (ftnlen)1)) {
#line 311 "AB05MD.f"
	*info = -2;
#line 312 "AB05MD.f"
    } else if (*n1 < 0) {
#line 313 "AB05MD.f"
	*info = -3;
#line 314 "AB05MD.f"
    } else if (*m1 < 0) {
#line 315 "AB05MD.f"
	*info = -4;
#line 316 "AB05MD.f"
    } else if (*p1 < 0) {
#line 317 "AB05MD.f"
	*info = -5;
#line 318 "AB05MD.f"
    } else if (*n2 < 0) {
#line 319 "AB05MD.f"
	*info = -6;
#line 320 "AB05MD.f"
    } else if (*p2 < 0) {
#line 321 "AB05MD.f"
	*info = -7;
#line 322 "AB05MD.f"
    } else if (*lda1 < max(1,*n1)) {
#line 323 "AB05MD.f"
	*info = -9;
#line 324 "AB05MD.f"
    } else if (*ldb1 < max(1,*n1)) {
#line 325 "AB05MD.f"
	*info = -11;
#line 326 "AB05MD.f"
    } else if (*n1 > 0 && *ldc1 < max(1,*p1) || *n1 == 0 && *ldc1 < 1) {
#line 328 "AB05MD.f"
	*info = -13;
#line 329 "AB05MD.f"
    } else if (*ldd1 < max(1,*p1)) {
#line 330 "AB05MD.f"
	*info = -15;
#line 331 "AB05MD.f"
    } else if (*lda2 < max(1,*n2)) {
#line 332 "AB05MD.f"
	*info = -17;
#line 333 "AB05MD.f"
    } else if (*ldb2 < max(1,*n2)) {
#line 334 "AB05MD.f"
	*info = -19;
#line 335 "AB05MD.f"
    } else if (*n2 > 0 && *ldc2 < max(1,*p2) || *n2 == 0 && *ldc2 < 1) {
#line 337 "AB05MD.f"
	*info = -21;
#line 338 "AB05MD.f"
    } else if (*ldd2 < max(1,*p2)) {
#line 339 "AB05MD.f"
	*info = -23;
#line 340 "AB05MD.f"
    } else if (*lda < max(1,*n)) {
#line 341 "AB05MD.f"
	*info = -26;
#line 342 "AB05MD.f"
    } else if (*ldb < max(1,*n)) {
#line 343 "AB05MD.f"
	*info = -28;
#line 344 "AB05MD.f"
    } else if (*n > 0 && *ldc < max(1,*p2) || *n == 0 && *ldc < 1) {
#line 346 "AB05MD.f"
	*info = -30;
#line 347 "AB05MD.f"
    } else if (*ldd < max(1,*p2)) {
#line 348 "AB05MD.f"
	*info = -32;
#line 349 "AB05MD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
#line 349 "AB05MD.f"
	i__3 = max(*n1,*m1), i__3 = max(i__3,*n2);
#line 349 "AB05MD.f"
	i__1 = 1, i__2 = *p1 * max(i__3,*p2);
#line 349 "AB05MD.f"
	if (lover && *ldwork < max(i__1,i__2) || ! lover && *ldwork < 1) {
#line 351 "AB05MD.f"
	    *info = -34;
#line 352 "AB05MD.f"
	}
#line 352 "AB05MD.f"
    }

#line 354 "AB05MD.f"
    if (*info != 0) {

/*        Error return. */

#line 358 "AB05MD.f"
	i__1 = -(*info);
#line 358 "AB05MD.f"
	xerbla_("AB05MD", &i__1, (ftnlen)6);
#line 359 "AB05MD.f"
	return 0;
#line 360 "AB05MD.f"
    }

/*     Quick return if possible. */

/* Computing MAX */
#line 364 "AB05MD.f"
    i__1 = *n, i__2 = min(*m1,*p2);
#line 364 "AB05MD.f"
    if (max(i__1,i__2) == 0) {
#line 364 "AB05MD.f"
	return 0;
#line 364 "AB05MD.f"
    }

/*     Set row/column indices for storing the results. */

#line 369 "AB05MD.f"
    if (luplo) {
#line 370 "AB05MD.f"
	i1 = 1;
/* Computing MIN */
#line 371 "AB05MD.f"
	i__1 = *n1 + 1;
#line 371 "AB05MD.f"
	i2 = min(i__1,*n);
#line 372 "AB05MD.f"
    } else {
/* Computing MIN */
#line 373 "AB05MD.f"
	i__1 = *n2 + 1;
#line 373 "AB05MD.f"
	i1 = min(i__1,*n);
#line 374 "AB05MD.f"
	i2 = 1;
#line 375 "AB05MD.f"
    }

#line 377 "AB05MD.f"
    ldwn2 = max(1,*n2);
#line 378 "AB05MD.f"
    ldwp1 = max(1,*p1);
#line 379 "AB05MD.f"
    ldwp2 = max(1,*p2);

/*     Construct the cascaded system matrices, taking the desired block */
/*     structure and possible overwriting into account. */

/*     Form the diagonal blocks of matrix  A. */

#line 386 "AB05MD.f"
    if (luplo) {

/*        Lower block diagonal structure. */

#line 390 "AB05MD.f"
	if (lover && *lda1 <= *lda) {
#line 391 "AB05MD.f"
	    if (*lda1 < *lda) {

#line 393 "AB05MD.f"
		for (j = *n1; j >= 1; --j) {
#line 394 "AB05MD.f"
		    for (i__ = *n1; i__ >= 1; --i__) {
#line 395 "AB05MD.f"
			a[i__ + j * a_dim1] = a1[i__ + j * a1_dim1];
#line 396 "AB05MD.f"
/* L10: */
#line 396 "AB05MD.f"
		    }
#line 397 "AB05MD.f"
/* L20: */
#line 397 "AB05MD.f"
		}

#line 399 "AB05MD.f"
	    }
#line 400 "AB05MD.f"
	} else {
#line 401 "AB05MD.f"
	    dlacpy_("F", n1, n1, &a1[a1_offset], lda1, &a[a_offset], lda, (
		    ftnlen)1);
#line 402 "AB05MD.f"
	}
#line 403 "AB05MD.f"
	if (*n2 > 0) {
#line 403 "AB05MD.f"
	    dlacpy_("F", n2, n2, &a2[a2_offset], lda2, &a[i2 + i2 * a_dim1], 
		    lda, (ftnlen)1);
#line 403 "AB05MD.f"
	}
#line 405 "AB05MD.f"
    } else {

/*        Upper block diagonal structure. */

#line 409 "AB05MD.f"
	if (lover && *lda2 <= *lda) {
#line 410 "AB05MD.f"
	    if (*lda2 < *lda) {

#line 412 "AB05MD.f"
		for (j = *n2; j >= 1; --j) {
#line 413 "AB05MD.f"
		    for (i__ = *n2; i__ >= 1; --i__) {
#line 414 "AB05MD.f"
			a[i__ + j * a_dim1] = a2[i__ + j * a2_dim1];
#line 415 "AB05MD.f"
/* L30: */
#line 415 "AB05MD.f"
		    }
#line 416 "AB05MD.f"
/* L40: */
#line 416 "AB05MD.f"
		}

#line 418 "AB05MD.f"
	    }
#line 419 "AB05MD.f"
	} else {
#line 420 "AB05MD.f"
	    dlacpy_("F", n2, n2, &a2[a2_offset], lda2, &a[a_offset], lda, (
		    ftnlen)1);
#line 421 "AB05MD.f"
	}
#line 422 "AB05MD.f"
	if (*n1 > 0) {
#line 422 "AB05MD.f"
	    dlacpy_("F", n1, n1, &a1[a1_offset], lda1, &a[i1 + i1 * a_dim1], 
		    lda, (ftnlen)1);
#line 422 "AB05MD.f"
	}
#line 424 "AB05MD.f"
    }

/*     Form the off-diagonal blocks of matrix  A. */

#line 428 "AB05MD.f"
    if (min(*n1,*n2) > 0) {
#line 429 "AB05MD.f"
	dlaset_("F", n1, n2, &c_b16, &c_b16, &a[i1 + i2 * a_dim1], lda, (
		ftnlen)1);
#line 430 "AB05MD.f"
	dgemm_("No transpose", "No transpose", n2, n1, p1, &c_b20, &b2[
		b2_offset], ldb2, &c1[c1_offset], ldc1, &c_b16, &a[i2 + i1 * 
		a_dim1], lda, (ftnlen)12, (ftnlen)12);
#line 432 "AB05MD.f"
    }

#line 434 "AB05MD.f"
    if (luplo) {

/*        Form the matrix  B. */

#line 438 "AB05MD.f"
	if (lover && *ldb1 <= *ldb) {
#line 439 "AB05MD.f"
	    if (*ldb1 < *ldb) {

#line 441 "AB05MD.f"
		for (j = *m1; j >= 1; --j) {
#line 442 "AB05MD.f"
		    for (i__ = *n1; i__ >= 1; --i__) {
#line 443 "AB05MD.f"
			b[i__ + j * b_dim1] = b1[i__ + j * b1_dim1];
#line 444 "AB05MD.f"
/* L50: */
#line 444 "AB05MD.f"
		    }
#line 445 "AB05MD.f"
/* L60: */
#line 445 "AB05MD.f"
		}

#line 447 "AB05MD.f"
	    }
#line 448 "AB05MD.f"
	} else {
#line 449 "AB05MD.f"
	    dlacpy_("F", n1, m1, &b1[b1_offset], ldb1, &b[b_offset], ldb, (
		    ftnlen)1);
#line 450 "AB05MD.f"
	}

#line 452 "AB05MD.f"
	if (min(*n2,*m1) > 0) {
#line 452 "AB05MD.f"
	    dgemm_("No transpose", "No transpose", n2, m1, p1, &c_b20, &b2[
		    b2_offset], ldb2, &d1[d1_offset], ldd1, &c_b16, &b[i2 + 
		    b_dim1], ldb, (ftnlen)12, (ftnlen)12);
#line 452 "AB05MD.f"
	}

/*        Form the matrix  C. */

#line 458 "AB05MD.f"
	if (*n1 > 0) {
#line 459 "AB05MD.f"
	    if (lover) {

/*              Workspace:  P1*N1. */

#line 463 "AB05MD.f"
		dlacpy_("F", p1, n1, &c1[c1_offset], ldc1, &dwork[1], &ldwp1, 
			(ftnlen)1);
#line 464 "AB05MD.f"
		dgemm_("No transpose", "No transpose", p2, n1, p1, &c_b20, &
			d2[d2_offset], ldd2, &dwork[1], &ldwp1, &c_b16, &c__[
			c_offset], ldc, (ftnlen)12, (ftnlen)12);
#line 466 "AB05MD.f"
	    } else {
#line 467 "AB05MD.f"
		dgemm_("No transpose", "No transpose", p2, n1, p1, &c_b20, &
			d2[d2_offset], ldd2, &c1[c1_offset], ldc1, &c_b16, &
			c__[c_offset], ldc, (ftnlen)12, (ftnlen)12);
#line 469 "AB05MD.f"
	    }
#line 470 "AB05MD.f"
	}

#line 472 "AB05MD.f"
	if (min(*p2,*n2) > 0) {
#line 472 "AB05MD.f"
	    dlacpy_("F", p2, n2, &c2[c2_offset], ldc2, &c__[i2 * c_dim1 + 1], 
		    ldc, (ftnlen)1);
#line 472 "AB05MD.f"
	}

/*        Now form the matrix  D. */

#line 477 "AB05MD.f"
	if (lover) {

/*           Workspace:  P1*M1. */

#line 481 "AB05MD.f"
	    dlacpy_("F", p1, m1, &d1[d1_offset], ldd1, &dwork[1], &ldwp1, (
		    ftnlen)1);
#line 482 "AB05MD.f"
	    dgemm_("No transpose", "No transpose", p2, m1, p1, &c_b20, &d2[
		    d2_offset], ldd2, &dwork[1], &ldwp1, &c_b16, &d__[
		    d_offset], ldd, (ftnlen)12, (ftnlen)12);
#line 484 "AB05MD.f"
	} else {
#line 485 "AB05MD.f"
	    dgemm_("No transpose", "No transpose", p2, m1, p1, &c_b20, &d2[
		    d2_offset], ldd2, &d1[d1_offset], ldd1, &c_b16, &d__[
		    d_offset], ldd, (ftnlen)12, (ftnlen)12);
#line 487 "AB05MD.f"
	}

#line 489 "AB05MD.f"
    } else {

/*        Form the matrix  B. */

#line 493 "AB05MD.f"
	if (lover) {

/*           Workspace:  N2*P1. */

#line 497 "AB05MD.f"
	    dlacpy_("F", n2, p1, &b2[b2_offset], ldb2, &dwork[1], &ldwn2, (
		    ftnlen)1);
#line 498 "AB05MD.f"
	    if (min(*n2,*m1) > 0) {
#line 498 "AB05MD.f"
		dgemm_("No transpose", "No transpose", n2, m1, p1, &c_b20, &
			dwork[1], &ldwn2, &d1[d1_offset], ldd1, &c_b16, &b[i2 
			+ b_dim1], ldb, (ftnlen)12, (ftnlen)12);
#line 498 "AB05MD.f"
	    }
#line 502 "AB05MD.f"
	} else {
#line 503 "AB05MD.f"
	    dgemm_("No transpose", "No transpose", n2, m1, p1, &c_b20, &b2[
		    b2_offset], ldb2, &d1[d1_offset], ldd1, &c_b16, &b[
		    b_offset], ldb, (ftnlen)12, (ftnlen)12);
#line 505 "AB05MD.f"
	}

#line 507 "AB05MD.f"
	if (min(*n1,*m1) > 0) {
#line 507 "AB05MD.f"
	    dlacpy_("F", n1, m1, &b1[b1_offset], ldb1, &b[i1 + b_dim1], ldb, (
		    ftnlen)1);
#line 507 "AB05MD.f"
	}

/*        Form the matrix  C. */

#line 512 "AB05MD.f"
	if (lover && *ldc2 <= *ldc) {
#line 513 "AB05MD.f"
	    if (*ldc2 < *ldc) {

#line 515 "AB05MD.f"
		for (j = *n2; j >= 1; --j) {
#line 516 "AB05MD.f"
		    for (i__ = *p2; i__ >= 1; --i__) {
#line 517 "AB05MD.f"
			c__[i__ + j * c_dim1] = c2[i__ + j * c2_dim1];
#line 518 "AB05MD.f"
/* L70: */
#line 518 "AB05MD.f"
		    }
#line 519 "AB05MD.f"
/* L80: */
#line 519 "AB05MD.f"
		}

#line 521 "AB05MD.f"
	    }
#line 522 "AB05MD.f"
	} else {
#line 523 "AB05MD.f"
	    dlacpy_("F", p2, n2, &c2[c2_offset], ldc2, &c__[c_offset], ldc, (
		    ftnlen)1);
#line 524 "AB05MD.f"
	}

#line 526 "AB05MD.f"
	if (min(*p2,*n1) > 0) {
#line 526 "AB05MD.f"
	    dgemm_("No transpose", "No transpose", p2, n1, p1, &c_b20, &d2[
		    d2_offset], ldd2, &c1[c1_offset], ldc1, &c_b16, &c__[i1 * 
		    c_dim1 + 1], ldc, (ftnlen)12, (ftnlen)12);
#line 526 "AB05MD.f"
	}

/*        Now form the matrix  D. */

#line 532 "AB05MD.f"
	if (lover) {

/*           Workspace:  P2*P1. */

#line 536 "AB05MD.f"
	    dlacpy_("F", p2, p1, &d2[d2_offset], ldd2, &dwork[1], &ldwp2, (
		    ftnlen)1);
#line 537 "AB05MD.f"
	    dgemm_("No transpose", "No transpose", p2, m1, p1, &c_b20, &dwork[
		    1], &ldwp2, &d1[d1_offset], ldd1, &c_b16, &d__[d_offset], 
		    ldd, (ftnlen)12, (ftnlen)12);
#line 539 "AB05MD.f"
	} else {
#line 540 "AB05MD.f"
	    dgemm_("No transpose", "No transpose", p2, m1, p1, &c_b20, &d2[
		    d2_offset], ldd2, &d1[d1_offset], ldd1, &c_b16, &d__[
		    d_offset], ldd, (ftnlen)12, (ftnlen)12);
#line 542 "AB05MD.f"
	}
#line 543 "AB05MD.f"
    }

#line 545 "AB05MD.f"
    return 0;
/* *** Last line of AB05MD *** */
} /* ab05md_ */

