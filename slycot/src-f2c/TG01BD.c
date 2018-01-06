#line 1 "TG01BD.f"
/* TG01BD.f -- translated by f2c (version 20100827).
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

#line 1 "TG01BD.f"
/* Table of constant values */

static doublereal c_b12 = 0.;
static doublereal c_b13 = 1.;
static integer c__1 = 1;

/* Subroutine */ int tg01bd_(char *jobe, char *compq, char *compz, integer *n,
	 integer *m, integer *p, integer *ilo, integer *ihi, doublereal *a, 
	integer *lda, doublereal *e, integer *lde, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *q, integer *ldq, 
	doublereal *z__, integer *ldz, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen jobe_len, ftnlen compq_len, ftnlen compz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, e_dim1, 
	    e_offset, q_dim1, q_offset, z_dim1, z_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal s, cs;
    static logical ilq, inq, ilz, inz;
    static integer jcol, ierr, itau;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer iwrk, jrow;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical withb, withc, upper;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), dlartg_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), xerbla_(
	    char *, integer *, ftnlen), dormqr_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static integer minwrk, maxwrk;


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

/*     To reduce the matrices A and E of the system pencil */

/*             S =  ( A  B ) - lambda ( E  0 ) , */
/*                  ( C  0 )          ( 0  0 ) */

/*     corresponding to the descriptor triple (A-lambda E,B,C), */
/*     to generalized upper Hessenberg form using orthogonal */
/*     transformations, */

/*          Q' * A * Z = H,   Q' * E * Z = T, */

/*     where H is upper Hessenberg, T is upper triangular, Q and Z */
/*     are orthogonal, and ' means transpose. The corresponding */
/*     transformations, written compactly as diag(Q',I) * S * diag(Z,I), */
/*     are also applied to B and C, getting Q' * B and C * Z. */

/*     The orthogonal matrices Q and Z are determined as products of */
/*     Givens rotations. They may either be formed explicitly, or they */
/*     may be postmultiplied into input matrices Q1 and Z1, so that */

/*          Q1 * A * Z1' = (Q1*Q) * H * (Z1*Z)' */
/*          Q1 * E * Z1' = (Q1*Q) * T * (Z1*Z)'. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBE    CHARACTER*1 */
/*             Specifies whether E is a general square or an upper */
/*             triangular matrix, as follows: */
/*             = 'G':  E is a general square matrix; */
/*             = 'U':  E is an upper triangular matrix. */

/*     COMPQ   CHARACTER*1 */
/*             Indicates what should be done with matrix Q, as follows: */
/*             = 'N':  do not compute Q; */
/*             = 'I':  Q is initialized to the unit matrix, and the */
/*                     orthogonal matrix Q is returned; */
/*             = 'V':  Q must contain an orthogonal matrix Q1 on entry, */
/*                     and the product Q1*Q is returned. */

/*     COMPZ   CHARACTER*1 */
/*             Indicates what should be done with matrix Z, as follows: */
/*             = 'N':  do not compute Z; */
/*             = 'I':  Z is initialized to the unit matrix, and the */
/*                     orthogonal matrix Z is returned; */
/*             = 'V':  Z must contain an orthogonal matrix Z1 on entry, */
/*                     and the product Z1*Z is returned. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, E, and the number of rows of */
/*             the matrix B.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrix B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of rows of the matrix C.  P >= 0. */

/*     ILO     (input) INTEGER */
/*     IHI     (input) INTEGER */
/*             It is assumed that A and E are already upper triangular in */
/*             rows and columns 1:ILO-1 and IHI+1:N.  ILO and IHI could */
/*             normally be set by a previous call to LAPACK Library */
/*             routine DGGBAL; otherwise they should be set to 1 and N, */
/*             respectively. */
/*             1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */
/*             If JOBE = 'U', the matrix E is assumed upper triangular. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the upper Hessenberg matrix H = Q' * A * Z. The elements */
/*             below the first subdiagonal are set to zero. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the descriptor matrix E. If JOBE = 'U', this */
/*             matrix is assumed upper triangular. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the upper triangular matrix T = Q' * E * Z. The elements */
/*             below the diagonal are set to zero. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input/state matrix B. */
/*             On exit, if M > 0, the leading N-by-M part of this array */
/*             contains the transformed matrix Q' * B. */
/*             The array B is not referenced if M = 0. */

/*     LDB     INTEGER */
/*             The leading dimension of array B. */
/*             LDB >= MAX(1,N) if M > 0;  LDB >= 1 if M = 0. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C. */
/*             On exit, if P > 0, the leading P-by-N part of this array */
/*             contains the transformed matrix C * Z. */
/*             The array C is not referenced if P = 0. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             If COMPQ = 'N':  Q is not referenced; */
/*             If COMPQ = 'I':  on entry, Q need not be set, and on exit */
/*                              it contains the orthogonal matrix Q, */
/*                              where Q' is the product of the Givens */
/*                              transformations which are applied to A, */
/*                              E, and B on the left; */
/*             If COMPQ = 'V':  on entry, Q must contain an orthogonal */
/*                              matrix Q1, and on exit this is */
/*                              overwritten by Q1*Q. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. */
/*             LDQ >= 1,        if COMPQ = 'N'; */
/*             LDQ >= MAX(1,N), if COMPQ = 'I' or 'V'. */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             If COMPZ = 'N':  Z is not referenced; */
/*             If COMPZ = 'I':  on entry, Z need not be set, and on exit */
/*                              it contains the orthogonal matrix Z, */
/*                              which is the product of the Givens */
/*                              transformations applied to A, E, and C */
/*                              on the right; */
/*             If COMPZ = 'V':  on entry, Z must contain an orthogonal */
/*                              matrix Z1, and on exit this is */
/*                              overwritten by Z1*Z. */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z. */
/*             LDZ >= 1,        if COMPZ = 'N'; */
/*             LDZ >= MAX(1,N), if COMPZ = 'I' or 'V'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= 1,                          if JOBE = 'U'; */
/*             LDWORK >= MAX(1,IHI+1-ILO+MAX(NI,M)), if JOBE = 'G', where */
/*             NI = N+1-ILO, if COMPQ = 'N', and NI = N, otherwise. */
/*             For good performance, if JOBE = 'G', LDWORK must generally */
/*             be larger, LDWORK >= MAX(1,IHI+1-ILO+MAX(NI,M)*NB), where */
/*             NB is the optimal block size. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit. */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     First, this routine computes the QR factorization of E and applies */
/*     the transformations to A, B, and possibly Q. Then, the routine */
/*     reduces A to upper Hessenberg form, preserving E triangular, by */
/*     an unblocked reduction [1], using two sequences of plane rotations */
/*     applied alternately from the left and from the right. The */
/*     corresponding transformations may be accumulated and/or applied */
/*     to the matrices B and C. If JOBE = 'U', the initial reduction of E */
/*     to upper triangular form is skipped. */

/*     This routine is a modification and extension of the LAPACK Library */
/*     routine DGGHRD [2]. */

/*     REFERENCES */

/*     [1] Golub, G.H. and van Loan, C.F. */
/*         Matrix Computations. Third Edition. */
/*         M. D. Johns Hopkins University Press, Baltimore, 1996. */

/*     [2] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J., */
/*         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A., */
/*         Ostrouchov, S., and Sorensen, D. */
/*         LAPACK Users' Guide: Second Edition. */
/*         SIAM, Philadelphia, 1995. */

/*     CONTRIBUTOR */

/*     D. Sima, University of Bucharest, May 2001. */
/*     V. Sima, Research Institute for Informatics, Bucharest, May 2001. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Eigenvalue, matrix algebra, matrix operations, similarity */
/*     transformation. */

/*  ********************************************************************* */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Test the input scalar parameters. */

#line 261 "TG01BD.f"
    /* Parameter adjustments */
#line 261 "TG01BD.f"
    a_dim1 = *lda;
#line 261 "TG01BD.f"
    a_offset = 1 + a_dim1;
#line 261 "TG01BD.f"
    a -= a_offset;
#line 261 "TG01BD.f"
    e_dim1 = *lde;
#line 261 "TG01BD.f"
    e_offset = 1 + e_dim1;
#line 261 "TG01BD.f"
    e -= e_offset;
#line 261 "TG01BD.f"
    b_dim1 = *ldb;
#line 261 "TG01BD.f"
    b_offset = 1 + b_dim1;
#line 261 "TG01BD.f"
    b -= b_offset;
#line 261 "TG01BD.f"
    c_dim1 = *ldc;
#line 261 "TG01BD.f"
    c_offset = 1 + c_dim1;
#line 261 "TG01BD.f"
    c__ -= c_offset;
#line 261 "TG01BD.f"
    q_dim1 = *ldq;
#line 261 "TG01BD.f"
    q_offset = 1 + q_dim1;
#line 261 "TG01BD.f"
    q -= q_offset;
#line 261 "TG01BD.f"
    z_dim1 = *ldz;
#line 261 "TG01BD.f"
    z_offset = 1 + z_dim1;
#line 261 "TG01BD.f"
    z__ -= z_offset;
#line 261 "TG01BD.f"
    --dwork;
#line 261 "TG01BD.f"

#line 261 "TG01BD.f"
    /* Function Body */
#line 261 "TG01BD.f"
    upper = lsame_(jobe, "U", (ftnlen)1, (ftnlen)1);
#line 262 "TG01BD.f"
    inq = lsame_(compq, "I", (ftnlen)1, (ftnlen)1);
#line 263 "TG01BD.f"
    ilq = lsame_(compq, "V", (ftnlen)1, (ftnlen)1) || inq;
#line 264 "TG01BD.f"
    inz = lsame_(compz, "I", (ftnlen)1, (ftnlen)1);
#line 265 "TG01BD.f"
    ilz = lsame_(compz, "V", (ftnlen)1, (ftnlen)1) || inz;
#line 266 "TG01BD.f"
    withb = *m > 0;
#line 267 "TG01BD.f"
    withc = *p > 0;

#line 269 "TG01BD.f"
    *info = 0;
#line 270 "TG01BD.f"
    if (! (upper || lsame_(jobe, "G", (ftnlen)1, (ftnlen)1))) {
#line 271 "TG01BD.f"
	*info = -1;
#line 272 "TG01BD.f"
    } else if (! (ilq || lsame_(compq, "N", (ftnlen)1, (ftnlen)1))) {
#line 273 "TG01BD.f"
	*info = -2;
#line 274 "TG01BD.f"
    } else if (! (ilz || lsame_(compz, "N", (ftnlen)1, (ftnlen)1))) {
#line 275 "TG01BD.f"
	*info = -3;
#line 276 "TG01BD.f"
    } else if (*n < 0) {
#line 277 "TG01BD.f"
	*info = -4;
#line 278 "TG01BD.f"
    } else if (*m < 0) {
#line 279 "TG01BD.f"
	*info = -5;
#line 280 "TG01BD.f"
    } else if (*p < 0) {
#line 281 "TG01BD.f"
	*info = -6;
#line 282 "TG01BD.f"
    } else if (*ilo < 1) {
#line 283 "TG01BD.f"
	*info = -7;
#line 284 "TG01BD.f"
    } else if (*ihi > *n || *ihi < *ilo - 1) {
#line 285 "TG01BD.f"
	*info = -8;
#line 286 "TG01BD.f"
    } else if (*lda < max(1,*n)) {
#line 287 "TG01BD.f"
	*info = -10;
#line 288 "TG01BD.f"
    } else if (*lde < max(1,*n)) {
#line 289 "TG01BD.f"
	*info = -12;
#line 290 "TG01BD.f"
    } else if (withb && *ldb < *n || *ldb < 1) {
#line 291 "TG01BD.f"
	*info = -14;
#line 292 "TG01BD.f"
    } else if (*ldc < max(1,*p)) {
#line 293 "TG01BD.f"
	*info = -16;
#line 294 "TG01BD.f"
    } else if (ilq && *ldq < *n || *ldq < 1) {
#line 295 "TG01BD.f"
	*info = -18;
#line 296 "TG01BD.f"
    } else if (ilz && *ldz < *n || *ldz < 1) {
#line 297 "TG01BD.f"
	*info = -20;
#line 298 "TG01BD.f"
    } else {
#line 299 "TG01BD.f"
	jrow = *ihi + 1 - *ilo;
#line 300 "TG01BD.f"
	jcol = *n + 1 - *ilo;
#line 301 "TG01BD.f"
	if (upper) {
#line 302 "TG01BD.f"
	    minwrk = 1;
#line 303 "TG01BD.f"
	    maxwrk = 1;
#line 304 "TG01BD.f"
	} else {
#line 305 "TG01BD.f"
	    if (ilq) {
#line 306 "TG01BD.f"
		minwrk = *n;
#line 307 "TG01BD.f"
	    } else {
#line 308 "TG01BD.f"
		minwrk = jcol;
#line 309 "TG01BD.f"
	    }
/* Computing MAX */
#line 310 "TG01BD.f"
	    i__1 = 1, i__2 = jrow + max(minwrk,*m);
#line 310 "TG01BD.f"
	    minwrk = max(i__1,i__2);
#line 311 "TG01BD.f"
	}
#line 312 "TG01BD.f"
	if (*ldwork < minwrk) {
#line 312 "TG01BD.f"
	    *info = -22;
#line 312 "TG01BD.f"
	}
#line 314 "TG01BD.f"
    }

#line 316 "TG01BD.f"
    if (*info != 0) {
#line 317 "TG01BD.f"
	i__1 = -(*info);
#line 317 "TG01BD.f"
	xerbla_("TG01BD", &i__1, (ftnlen)6);
#line 318 "TG01BD.f"
	return 0;
#line 319 "TG01BD.f"
    }

/*     Initialize Q and Z if desired. */

#line 323 "TG01BD.f"
    if (inq) {
#line 323 "TG01BD.f"
	dlaset_("Full", n, n, &c_b12, &c_b13, &q[q_offset], ldq, (ftnlen)4);
#line 323 "TG01BD.f"
    }
#line 325 "TG01BD.f"
    if (inz) {
#line 325 "TG01BD.f"
	dlaset_("Full", n, n, &c_b12, &c_b13, &z__[z_offset], ldz, (ftnlen)4);
#line 325 "TG01BD.f"
    }

/*     Quick return if possible. */

#line 330 "TG01BD.f"
    if (*n <= 1) {
#line 331 "TG01BD.f"
	dwork[1] = 1.;
#line 332 "TG01BD.f"
	return 0;
#line 333 "TG01BD.f"
    }

#line 335 "TG01BD.f"
    if (! upper) {

/*        Reduce E to triangular form (QR decomposition of E). */

/*        (Note: Comments in the code beginning "Workspace:" describe the */
/*        minimal amount of real workspace needed at that point in the */
/*        code, as well as the preferred amount for good performance. */
/*        NB refers to the optimal block size for the immediately */
/*        following subroutine, as returned by ILAENV.) */

/*        Workspace: need   IHI+1-ILO+N+1-ILO; */
/*                   prefer IHI+1-ILO+(N+1-ILO)*NB. */

#line 348 "TG01BD.f"
	itau = 1;
#line 349 "TG01BD.f"
	iwrk = itau + jrow;
#line 350 "TG01BD.f"
	i__1 = *ldwork - iwrk + 1;
#line 350 "TG01BD.f"
	dgeqrf_(&jrow, &jcol, &e[*ilo + *ilo * e_dim1], lde, &dwork[itau], &
		dwork[iwrk], &i__1, &ierr);
/* Computing MAX */
#line 352 "TG01BD.f"
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 352 "TG01BD.f"
	maxwrk = max(i__1,minwrk);

/*        Apply the orthogonal transformation to matrices A, B, and Q. */
/*        Workspace: need   IHI+1-ILO+N+1-ILO; */
/*                   prefer IHI+1-ILO+(N+1-ILO)*NB. */

#line 358 "TG01BD.f"
	i__1 = *ldwork - iwrk + 1;
#line 358 "TG01BD.f"
	dormqr_("Left", "Transpose", &jrow, &jcol, &jrow, &e[*ilo + *ilo * 
		e_dim1], lde, &dwork[itau], &a[*ilo + *ilo * a_dim1], lda, &
		dwork[iwrk], &i__1, &ierr, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 361 "TG01BD.f"
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 361 "TG01BD.f"
	maxwrk = max(i__1,maxwrk);

#line 363 "TG01BD.f"
	if (withb) {

/*           Workspace: need   IHI+1-ILO+M; */
/*                      prefer IHI+1-ILO+M*NB. */

#line 368 "TG01BD.f"
	    i__1 = *ldwork - iwrk + 1;
#line 368 "TG01BD.f"
	    dormqr_("Left", "Transpose", &jrow, m, &jrow, &e[*ilo + *ilo * 
		    e_dim1], lde, &dwork[itau], &b[*ilo + b_dim1], ldb, &
		    dwork[iwrk], &i__1, &ierr, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 371 "TG01BD.f"
	    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 371 "TG01BD.f"
	    maxwrk = max(i__1,maxwrk);
#line 372 "TG01BD.f"
	}

#line 374 "TG01BD.f"
	if (ilq) {

/*           Workspace: need   IHI+1-ILO+N; */
/*                      prefer IHI+1-ILO+N*NB. */

#line 379 "TG01BD.f"
	    i__1 = *ldwork - iwrk + 1;
#line 379 "TG01BD.f"
	    dormqr_("Right", "No Transpose", n, &jrow, &jrow, &e[*ilo + *ilo *
		     e_dim1], lde, &dwork[itau], &q[*ilo * q_dim1 + 1], ldq, &
		    dwork[iwrk], &i__1, &ierr, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
#line 382 "TG01BD.f"
	    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 382 "TG01BD.f"
	    maxwrk = max(i__1,maxwrk);
#line 383 "TG01BD.f"
	}
#line 384 "TG01BD.f"
    }

/*     Zero out lower triangle of E. */

#line 388 "TG01BD.f"
    if (jrow > 1) {
#line 388 "TG01BD.f"
	i__1 = jrow - 1;
#line 388 "TG01BD.f"
	i__2 = jrow - 1;
#line 388 "TG01BD.f"
	dlaset_("Lower", &i__1, &i__2, &c_b12, &c_b12, &e[*ilo + 1 + *ilo * 
		e_dim1], lde, (ftnlen)5);
#line 388 "TG01BD.f"
    }

/*     Reduce A and E and apply the transformations to B, C, Q and Z. */

#line 394 "TG01BD.f"
    i__1 = *ihi - 2;
#line 394 "TG01BD.f"
    for (jcol = *ilo; jcol <= i__1; ++jcol) {

#line 396 "TG01BD.f"
	i__2 = jcol + 2;
#line 396 "TG01BD.f"
	for (jrow = *ihi; jrow >= i__2; --jrow) {

/*           Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL). */

#line 400 "TG01BD.f"
	    temp = a[jrow - 1 + jcol * a_dim1];
#line 401 "TG01BD.f"
	    dlartg_(&temp, &a[jrow + jcol * a_dim1], &cs, &s, &a[jrow - 1 + 
		    jcol * a_dim1]);
#line 403 "TG01BD.f"
	    a[jrow + jcol * a_dim1] = 0.;
#line 404 "TG01BD.f"
	    i__3 = *n - jcol;
#line 404 "TG01BD.f"
	    drot_(&i__3, &a[jrow - 1 + (jcol + 1) * a_dim1], lda, &a[jrow + (
		    jcol + 1) * a_dim1], lda, &cs, &s);
#line 406 "TG01BD.f"
	    i__3 = *n + 2 - jrow;
#line 406 "TG01BD.f"
	    drot_(&i__3, &e[jrow - 1 + (jrow - 1) * e_dim1], lde, &e[jrow + (
		    jrow - 1) * e_dim1], lde, &cs, &s);
#line 408 "TG01BD.f"
	    if (withb) {
#line 408 "TG01BD.f"
		drot_(m, &b[jrow - 1 + b_dim1], ldb, &b[jrow + b_dim1], ldb, &
			cs, &s);
#line 408 "TG01BD.f"
	    }
#line 411 "TG01BD.f"
	    if (ilq) {
#line 411 "TG01BD.f"
		drot_(n, &q[(jrow - 1) * q_dim1 + 1], &c__1, &q[jrow * q_dim1 
			+ 1], &c__1, &cs, &s);
#line 411 "TG01BD.f"
	    }

/*           Step 2: rotate columns JROW, JROW-1 to kill E(JROW,JROW-1). */

#line 416 "TG01BD.f"
	    temp = e[jrow + jrow * e_dim1];
#line 417 "TG01BD.f"
	    dlartg_(&temp, &e[jrow + (jrow - 1) * e_dim1], &cs, &s, &e[jrow + 
		    jrow * e_dim1]);
#line 419 "TG01BD.f"
	    e[jrow + (jrow - 1) * e_dim1] = 0.;
#line 420 "TG01BD.f"
	    drot_(ihi, &a[jrow * a_dim1 + 1], &c__1, &a[(jrow - 1) * a_dim1 + 
		    1], &c__1, &cs, &s);
#line 421 "TG01BD.f"
	    i__3 = jrow - 1;
#line 421 "TG01BD.f"
	    drot_(&i__3, &e[jrow * e_dim1 + 1], &c__1, &e[(jrow - 1) * e_dim1 
		    + 1], &c__1, &cs, &s);
#line 423 "TG01BD.f"
	    if (withc) {
#line 423 "TG01BD.f"
		drot_(p, &c__[jrow * c_dim1 + 1], &c__1, &c__[(jrow - 1) * 
			c_dim1 + 1], &c__1, &cs, &s);
#line 423 "TG01BD.f"
	    }
#line 425 "TG01BD.f"
	    if (ilz) {
#line 425 "TG01BD.f"
		drot_(n, &z__[jrow * z_dim1 + 1], &c__1, &z__[(jrow - 1) * 
			z_dim1 + 1], &c__1, &cs, &s);
#line 425 "TG01BD.f"
	    }
#line 427 "TG01BD.f"
/* L10: */
#line 427 "TG01BD.f"
	}

#line 429 "TG01BD.f"
/* L20: */
#line 429 "TG01BD.f"
    }

#line 431 "TG01BD.f"
    dwork[1] = (doublereal) maxwrk;
#line 432 "TG01BD.f"
    return 0;
/* *** Last line of TG01BD *** */
} /* tg01bd_ */

