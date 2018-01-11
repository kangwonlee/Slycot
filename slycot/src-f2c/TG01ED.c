#line 1 "TG01ED.f"
/* TG01ED.f -- translated by f2c (version 20100827).
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

#line 1 "TG01ED.f"
/* Table of constant values */

static doublereal c_b6 = 0.;
static doublereal c_b7 = 1.;
static integer c__1 = 1;

/* Subroutine */ int tg01ed_(char *joba, integer *l, integer *n, integer *m, 
	integer *p, doublereal *a, integer *lda, doublereal *e, integer *lde, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *q, integer *ldq, doublereal *z__, integer *ldz, integer *
	ranke, integer *rnka22, doublereal *tol, doublereal *dwork, integer *
	ldwork, integer *info, ftnlen joba_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, e_dim1, 
	    e_offset, q_dim1, q_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, ln, kw, ir1, ln2, la22, na22, lwr;
    static logical reda;
    static doublereal epsm;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), dgemm_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen), mb03ud_(char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), dswap_(integer 
	    *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgelqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dgeqrf_(integer *, integer *, doublereal *, integer *, doublereal 
	    *, doublereal *, integer *, integer *), dgesvd_(char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), dormlq_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen);
    static doublereal svemax;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static doublereal svlmax;
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

/*     To compute for the descriptor system (A-lambda E,B,C) */
/*     the orthogonal transformation matrices Q and Z such that the */
/*     transformed system (Q'*A*Z-lambda Q'*E*Z, Q'*B, C*Z) is in an */
/*     SVD (singular value decomposition) coordinate form with */
/*     the system matrices  Q'*A*Z and Q'*E*Z in the form */

/*                  ( A11  A12 )             ( Er  0 ) */
/*         Q'*A*Z = (          ) ,  Q'*E*Z = (       ) , */
/*                  ( A21  A22 )             (  0  0 ) */

/*     where Er is an invertible diagonal matrix having on the diagonal */
/*     the decreasingly ordered nonzero singular values of E. */
/*     Optionally, the A22 matrix can be further reduced to the */
/*     SVD form */

/*                  ( Ar  0 ) */
/*            A22 = (       ) , */
/*                  (  0  0 ) */

/*     where Ar is an invertible diagonal matrix having on the diagonal */
/*     the decreasingly ordered nonzero singular values of A22. */
/*     The left and/or right orthogonal transformations performed */
/*     to reduce E and A22 are accumulated. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBA    CHARACTER*1 */
/*             = 'N':  do not reduce A22; */
/*             = 'R':  reduce A22 to an SVD form. */

/*     Input/Output Parameters */

/*     L       (input) INTEGER */
/*             The number of rows of matrices A, B, and E.  L >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of matrices A, E, and C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of matrix B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of rows of matrix C.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the state dynamics matrix A. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the transformed matrix Q'*A*Z. If JOBA = 'R', this matrix */
/*             is in the form */

/*                           ( A11  *   *  ) */
/*                  Q'*A*Z = (  *   Ar  0  ) , */
/*                           (  *   0   0  ) */

/*             where A11 is a RANKE-by-RANKE matrix and Ar is a */
/*             RNKA22-by-RNKA22 invertible diagonal matrix, with */
/*             decresingly ordered positive diagonal elements. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,L). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the descriptor matrix E. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the transformed matrix Q'*E*Z. */

/*                      ( Er  0 ) */
/*             Q'*E*Z = (       ) , */
/*                      (  0  0 ) */

/*             where Er is a RANKE-by-RANKE invertible diagonal matrix */
/*             having on the diagonal the decreasingly ordered positive */
/*             singular values of E. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,L). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading L-by-M part of this array must */
/*             contain the input/state matrix B. */
/*             On exit, the leading L-by-M part of this array contains */
/*             the transformed matrix Q'*B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B. */
/*             LDB >= MAX(1,L) if M > 0 or LDB >= 1 if M = 0. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the transformed matrix C*Z. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     Q       (output) DOUBLE PRECISION array, dimension (LDQ,L) */
/*             The leading L-by-L part of this array contains the */
/*             orthogonal matrix Q, which is the accumulated product of */
/*             transformations applied to A, E, and B on the left. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q.  LDQ >= MAX(1,L). */

/*     Z       (output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             The leading N-by-N part of this array contains the */
/*             orthogonal matrix Z, which is the accumulated product of */
/*             transformations applied to A, E, and C on the right. */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z.  LDZ >= MAX(1,N). */

/*     RANKE   (output) INTEGER */
/*             The effective rank of matrix E, and thus also the order */
/*             of the invertible diagonal submatrix Er. */
/*             RANKE is computed as the number of singular values of E */
/*             greater than TOL*SVEMAX, where SVEMAX is the maximum */
/*             singular value of E. */

/*     RNKA22  (output) INTEGER */
/*             If JOBA = 'R', then RNKA22 is the effective rank of */
/*             matrix A22, and thus also the order of the invertible */
/*             diagonal submatrix Ar. RNKA22 is computed as the number */
/*             of singular values of A22 greater than TOL*SVAMAX, */
/*             where SVAMAX is an estimate of the maximum singular value */
/*             of A. */
/*             If JOBA = 'N', then RNKA22 is not referenced. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in determining the rank of E */
/*             and of A22. If TOL > 0, then singular values less than */
/*             TOL*SVMAX are treated as zero, where SVMAX is the maximum */
/*             singular value of E or an estimate of it for A and E. */
/*             If TOL <= 0, the default tolerance TOLDEF = EPS*L*N is */
/*             used instead, where EPS is the machine precision */
/*             (see LAPACK Library routine DLAMCH). TOL < 1. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,MIN(L,N) + */
/*                           MAX(3*MIN(L,N)+MAX(L,N), 5*MIN(L,N), M, P)). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  the QR algorithm has failed to converge when computing */
/*                   singular value decomposition. In this case INFO */
/*                   specifies how many superdiagonals did not converge. */
/*                   This failure is not likely to occur. */

/*     METHOD */

/*     The routine computes the singular value decomposition (SVD) of E, */
/*     in the form */

/*                    ( Er  0 ) */
/*           E  = Q * (       ) * Z' */
/*                    (  0  0 ) */

/*     and finds the largest RANKE-by-RANKE leading diagonal submatrix */
/*     Er whose condition number is less than 1/TOL. RANKE defines thus */
/*     the effective rank of matrix E. */
/*     If JOBA = 'R' the same reduction is performed on A22 in the */
/*     partitioned matrix */

/*                  ( A11  A12 ) */
/*         Q'*A*Z = (          ) , */
/*                  ( A21  A22 ) */

/*     to obtain it in the form */

/*                  ( Ar  0 ) */
/*            A22 = (       ) , */
/*                  (  0  0 ) */

/*     with Ar an invertible diagonal matrix. */

/*     The accumulated transformations are also applied to the rest of */
/*     matrices */

/*          B <- Q' * B,  C <- C * Z. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically backward stable and requires */
/*     0( L*L*N )  floating point operations. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     March 1999. Based on the RASP routine RPDSSV. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, July 1999, */
/*     Feb. 2000, Oct. 2001, May 2003. */

/*     KEYWORDS */

/*     Descriptor system, matrix algebra, matrix operations, */
/*     orthogonal transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

#line 273 "TG01ED.f"
    /* Parameter adjustments */
#line 273 "TG01ED.f"
    a_dim1 = *lda;
#line 273 "TG01ED.f"
    a_offset = 1 + a_dim1;
#line 273 "TG01ED.f"
    a -= a_offset;
#line 273 "TG01ED.f"
    e_dim1 = *lde;
#line 273 "TG01ED.f"
    e_offset = 1 + e_dim1;
#line 273 "TG01ED.f"
    e -= e_offset;
#line 273 "TG01ED.f"
    b_dim1 = *ldb;
#line 273 "TG01ED.f"
    b_offset = 1 + b_dim1;
#line 273 "TG01ED.f"
    b -= b_offset;
#line 273 "TG01ED.f"
    c_dim1 = *ldc;
#line 273 "TG01ED.f"
    c_offset = 1 + c_dim1;
#line 273 "TG01ED.f"
    c__ -= c_offset;
#line 273 "TG01ED.f"
    q_dim1 = *ldq;
#line 273 "TG01ED.f"
    q_offset = 1 + q_dim1;
#line 273 "TG01ED.f"
    q -= q_offset;
#line 273 "TG01ED.f"
    z_dim1 = *ldz;
#line 273 "TG01ED.f"
    z_offset = 1 + z_dim1;
#line 273 "TG01ED.f"
    z__ -= z_offset;
#line 273 "TG01ED.f"
    --dwork;
#line 273 "TG01ED.f"

#line 273 "TG01ED.f"
    /* Function Body */
#line 273 "TG01ED.f"
    reda = lsame_(joba, "R", (ftnlen)1, (ftnlen)1);

/*     Test the input parameters. */

#line 277 "TG01ED.f"
    *info = 0;
/* Computing MAX */
#line 278 "TG01ED.f"
    i__1 = max(*m,*p), i__2 = min(*l,*n) * 3 + max(*l,*n), i__1 = max(i__1,
	    i__2), i__2 = min(*l,*n) * 5;
#line 278 "TG01ED.f"
    wrkopt = min(*l,*n) + max(i__1,i__2);
#line 280 "TG01ED.f"
    if (! lsame_(joba, "N", (ftnlen)1, (ftnlen)1) && ! reda) {
#line 281 "TG01ED.f"
	*info = -1;
#line 282 "TG01ED.f"
    } else if (*l < 0) {
#line 283 "TG01ED.f"
	*info = -2;
#line 284 "TG01ED.f"
    } else if (*n < 0) {
#line 285 "TG01ED.f"
	*info = -3;
#line 286 "TG01ED.f"
    } else if (*m < 0) {
#line 287 "TG01ED.f"
	*info = -4;
#line 288 "TG01ED.f"
    } else if (*p < 0) {
#line 289 "TG01ED.f"
	*info = -5;
#line 290 "TG01ED.f"
    } else if (*lda < max(1,*l)) {
#line 291 "TG01ED.f"
	*info = -7;
#line 292 "TG01ED.f"
    } else if (*lde < max(1,*l)) {
#line 293 "TG01ED.f"
	*info = -9;
#line 294 "TG01ED.f"
    } else if (*ldb < 1 || *m > 0 && *ldb < *l) {
#line 295 "TG01ED.f"
	*info = -11;
#line 296 "TG01ED.f"
    } else if (*ldc < max(1,*p)) {
#line 297 "TG01ED.f"
	*info = -13;
#line 298 "TG01ED.f"
    } else if (*ldq < max(1,*l)) {
#line 299 "TG01ED.f"
	*info = -15;
#line 300 "TG01ED.f"
    } else if (*ldz < max(1,*n)) {
#line 301 "TG01ED.f"
	*info = -17;
#line 302 "TG01ED.f"
    } else if (*tol >= 1.) {
#line 303 "TG01ED.f"
	*info = -20;
#line 304 "TG01ED.f"
    } else if (*ldwork < max(1,wrkopt)) {
#line 305 "TG01ED.f"
	*info = -22;
#line 306 "TG01ED.f"
    }
#line 307 "TG01ED.f"
    if (*info != 0) {
#line 308 "TG01ED.f"
	i__1 = -(*info);
#line 308 "TG01ED.f"
	xerbla_("TG01ED", &i__1, (ftnlen)6);
#line 309 "TG01ED.f"
	return 0;
#line 310 "TG01ED.f"
    }

/*     Quick return if possible. */

#line 314 "TG01ED.f"
    if (*l == 0 || *n == 0) {
#line 315 "TG01ED.f"
	if (*l > 0) {
#line 315 "TG01ED.f"
	    dlaset_("Full", l, l, &c_b6, &c_b7, &q[q_offset], ldq, (ftnlen)4);
#line 315 "TG01ED.f"
	}
#line 317 "TG01ED.f"
	if (*n > 0) {
#line 317 "TG01ED.f"
	    dlaset_("Full", n, n, &c_b6, &c_b7, &z__[z_offset], ldz, (ftnlen)
		    4);
#line 317 "TG01ED.f"
	}
#line 319 "TG01ED.f"
	dwork[1] = 1.;
#line 320 "TG01ED.f"
	*ranke = 0;
#line 321 "TG01ED.f"
	if (reda) {
#line 321 "TG01ED.f"
	    *rnka22 = 0;
#line 321 "TG01ED.f"
	}
#line 322 "TG01ED.f"
	return 0;
#line 323 "TG01ED.f"
    }

#line 325 "TG01ED.f"
    ln = min(*l,*n);
#line 326 "TG01ED.f"
    epsm = dlamch_("EPSILON", (ftnlen)7);

#line 328 "TG01ED.f"
    toldef = *tol;
#line 329 "TG01ED.f"
    if (toldef <= 0.) {

/*        Use the default tolerance for rank determination. */

#line 333 "TG01ED.f"
	toldef = epsm * (doublereal) (*l * *n);
#line 334 "TG01ED.f"
    }

/*     Set the estimate of the maximum singular value of E to */
/*     max(||E||,||A||) to detect negligible A or E matrices. */

/* Computing MAX */
#line 339 "TG01ED.f"
    d__1 = dlange_("F", l, n, &e[e_offset], lde, &dwork[1], (ftnlen)1), d__2 =
	     dlange_("F", l, n, &a[a_offset], lda, &dwork[1], (ftnlen)1);
#line 339 "TG01ED.f"
    svlmax = max(d__1,d__2);

/*     Compute the SVD of E */

/*                    ( Er  0 ) */
/*           E = Qr * (       ) * Zr' */
/*                    (  0  0 ) */

/*     Workspace: needed  MIN(L,N) + MAX(3*MIN(L,N)+MAX(L,N),5*MIN(L,N)); */
/*                prefer larger. */

#line 351 "TG01ED.f"
    lwr = *ldwork - ln;
#line 352 "TG01ED.f"
    kw = ln + 1;

#line 354 "TG01ED.f"
    dgesvd_("A", "A", l, n, &e[e_offset], lde, &dwork[1], &q[q_offset], ldq, &
	    z__[z_offset], ldz, &dwork[kw], &lwr, info, (ftnlen)1, (ftnlen)1);
#line 356 "TG01ED.f"
    if (*info > 0) {
#line 356 "TG01ED.f"
	return 0;
#line 356 "TG01ED.f"
    }
/* Computing MAX */
#line 358 "TG01ED.f"
    i__1 = wrkopt, i__2 = ln + (integer) dwork[kw];
#line 358 "TG01ED.f"
    wrkopt = max(i__1,i__2);

/*     Determine the rank of E. */

#line 362 "TG01ED.f"
    *ranke = 0;
#line 363 "TG01ED.f"
    if (dwork[1] > svlmax * epsm) {
#line 364 "TG01ED.f"
	*ranke = 1;
#line 365 "TG01ED.f"
	svemax = dwork[1];
#line 366 "TG01ED.f"
	i__1 = ln;
#line 366 "TG01ED.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 367 "TG01ED.f"
	    if (dwork[i__] < svemax * toldef) {
#line 367 "TG01ED.f"
		goto L20;
#line 367 "TG01ED.f"
	    }
#line 368 "TG01ED.f"
	    ++(*ranke);
#line 369 "TG01ED.f"
/* L10: */
#line 369 "TG01ED.f"
	}

#line 371 "TG01ED.f"
L20:
#line 372 "TG01ED.f"
	;
#line 372 "TG01ED.f"
    }

/*     Apply transformation on the rest of matrices. */

#line 376 "TG01ED.f"
    if (*ranke > 0) {

/*        A <-- Qr' * A * Zr. */

#line 380 "TG01ED.f"
	dgemm_("Transpose", "No transpose", l, n, l, &c_b7, &q[q_offset], ldq,
		 &a[a_offset], lda, &c_b6, &e[e_offset], lde, (ftnlen)9, (
		ftnlen)12);
#line 382 "TG01ED.f"
	dgemm_("No transpose", "Transpose", l, n, n, &c_b7, &e[e_offset], lde,
		 &z__[z_offset], ldz, &c_b6, &a[a_offset], lda, (ftnlen)12, (
		ftnlen)9);

/*        B <-- Qr' * B. */
/*        Workspace: need   L; */
/*                   prefer L*M. */

#line 389 "TG01ED.f"
	if (lwr > *l * *m && *m > 0) {

#line 391 "TG01ED.f"
	    dgemm_("Transpose", "No transpose", l, m, l, &c_b7, &q[q_offset], 
		    ldq, &b[b_offset], ldb, &c_b6, &dwork[kw], l, (ftnlen)9, (
		    ftnlen)12);
#line 393 "TG01ED.f"
	    dlacpy_("Full", l, m, &dwork[kw], l, &b[b_offset], ldb, (ftnlen)4)
		    ;
#line 394 "TG01ED.f"
	} else {
#line 395 "TG01ED.f"
	    i__1 = *m;
#line 395 "TG01ED.f"
	    for (j = 1; j <= i__1; ++j) {
#line 396 "TG01ED.f"
		dgemv_("Transpose", l, l, &c_b7, &q[q_offset], ldq, &b[j * 
			b_dim1 + 1], &c__1, &c_b6, &dwork[kw], &c__1, (ftnlen)
			9);
#line 398 "TG01ED.f"
		dcopy_(l, &dwork[kw], &c__1, &b[j * b_dim1 + 1], &c__1);
#line 399 "TG01ED.f"
/* L30: */
#line 399 "TG01ED.f"
	    }
#line 400 "TG01ED.f"
	}

/*        C <-- C * Zr. */
/*        Workspace: need   N; */
/*                   prefer P*N. */

#line 406 "TG01ED.f"
	if (lwr > *p * *n) {

#line 408 "TG01ED.f"
	    i__1 = max(1,*p);
#line 408 "TG01ED.f"
	    dgemm_("No transpose", "Transpose", p, n, n, &c_b7, &c__[c_offset]
		    , ldc, &z__[z_offset], ldz, &c_b6, &dwork[kw], &i__1, (
		    ftnlen)12, (ftnlen)9);
#line 410 "TG01ED.f"
	    i__1 = max(1,*p);
#line 410 "TG01ED.f"
	    dlacpy_("Full", p, n, &dwork[kw], &i__1, &c__[c_offset], ldc, (
		    ftnlen)4);
#line 411 "TG01ED.f"
	} else {
#line 412 "TG01ED.f"
	    i__1 = *p;
#line 412 "TG01ED.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 413 "TG01ED.f"
		dgemv_("No transpose", n, n, &c_b7, &z__[z_offset], ldz, &c__[
			i__ + c_dim1], ldc, &c_b6, &dwork[kw], &c__1, (ftnlen)
			12);
#line 415 "TG01ED.f"
		dcopy_(n, &dwork[kw], &c__1, &c__[i__ + c_dim1], ldc);
#line 416 "TG01ED.f"
/* L40: */
#line 416 "TG01ED.f"
	    }
#line 417 "TG01ED.f"
	}
/* Computing MAX */
#line 418 "TG01ED.f"
	i__1 = wrkopt, i__2 = *l * *m, i__1 = max(i__1,i__2), i__2 = *p * *n;
#line 418 "TG01ED.f"
	wrkopt = max(i__1,i__2);
#line 419 "TG01ED.f"
    }

/*     Reduce A22 if necessary. */

#line 423 "TG01ED.f"
    if (reda) {
#line 424 "TG01ED.f"
	la22 = *l - *ranke;
#line 425 "TG01ED.f"
	na22 = *n - *ranke;
#line 426 "TG01ED.f"
	ln2 = min(la22,na22);
#line 427 "TG01ED.f"
	if (ln2 == 0) {
#line 428 "TG01ED.f"
	    ir1 = 1;
#line 429 "TG01ED.f"
	    *rnka22 = 0;
#line 430 "TG01ED.f"
	} else {

/*           Compute the SVD of A22 using a storage saving approach. */

#line 434 "TG01ED.f"
	    ir1 = *ranke + 1;
#line 435 "TG01ED.f"
	    if (la22 >= na22) {

/*              Compute the QR decomposition of A22 in the form */

/*              A22 = Q2 * ( R2 ) , */
/*                         ( 0  ) */

/*              where R2 is upper triangular. */
/*              Workspace: need   MIN(L,N) + N; */
/*                         prefer MIN(L,N) + N*NB. */

#line 446 "TG01ED.f"
		dgeqrf_(&la22, &na22, &a[ir1 + ir1 * a_dim1], lda, &dwork[ir1]
			, &dwork[kw], &lwr, info);
/* Computing MAX */
#line 448 "TG01ED.f"
		i__1 = wrkopt, i__2 = ln + (integer) dwork[kw];
#line 448 "TG01ED.f"
		wrkopt = max(i__1,i__2);

/*              Apply transformation Q2 to A, B, and Q. */

/*              A <--diag(I, Q2') * A */
/*              Workspace: need   MIN(L,N) + N; */
/*                         prefer MIN(L,N) + N*NB. */

#line 456 "TG01ED.f"
		dormqr_("Left", "Transpose", &la22, ranke, &ln2, &a[ir1 + ir1 
			* a_dim1], lda, &dwork[ir1], &a[ir1 + a_dim1], lda, &
			dwork[kw], &lwr, info, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 459 "TG01ED.f"
		i__1 = wrkopt, i__2 = ln + (integer) dwork[kw];
#line 459 "TG01ED.f"
		wrkopt = max(i__1,i__2);

/*              B <-- diag(I, Q2') * B */
/*              Workspace: need   MIN(L,N) + M; */
/*                         prefer MIN(L,N) + M*NB. */

#line 465 "TG01ED.f"
		if (*m > 0) {
#line 466 "TG01ED.f"
		    dormqr_("Left", "Transpose", &la22, m, &ln2, &a[ir1 + ir1 
			    * a_dim1], lda, &dwork[ir1], &b[ir1 + b_dim1], 
			    ldb, &dwork[kw], &lwr, info, (ftnlen)4, (ftnlen)9)
			    ;
/* Computing MAX */
#line 469 "TG01ED.f"
		    i__1 = wrkopt, i__2 = ln + (integer) dwork[kw];
#line 469 "TG01ED.f"
		    wrkopt = max(i__1,i__2);
#line 470 "TG01ED.f"
		}

/*              Q <-- Q * diag(I, Q2) */
/*              Workspace: need   MIN(L,N) + L; */
/*                         prefer MIN(L,N) + L*NB. */

#line 476 "TG01ED.f"
		dormqr_("Right", "No transpose", l, &la22, &ln2, &a[ir1 + ir1 
			* a_dim1], lda, &dwork[ir1], &q[ir1 * q_dim1 + 1], 
			ldq, &dwork[kw], &lwr, info, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
#line 479 "TG01ED.f"
		i__1 = wrkopt, i__2 = ln + (integer) dwork[kw];
#line 479 "TG01ED.f"
		wrkopt = max(i__1,i__2);

/*              Compute the SVD of the upper triangular submatrix R2 as */

/*                               ( Ar  0 ) */
/*                    R2 = Q2r * (       ) * Z2r' , */
/*                               (  0  0 ) */

/*              where Q2r is stored in E and Z2r' is stored in A22. */
/*              Workspace: need   MAX(1,5*MIN(L,N)); */
/*                         prefer larger. */

#line 491 "TG01ED.f"
		mb03ud_("Vectors", "Vectors", &ln2, &a[ir1 + ir1 * a_dim1], 
			lda, &e[ir1 + ir1 * e_dim1], lde, &dwork[ir1], &dwork[
			kw], &lwr, info, (ftnlen)7, (ftnlen)7);
#line 494 "TG01ED.f"
		if (*info > 0) {
#line 494 "TG01ED.f"
		    return 0;
#line 494 "TG01ED.f"
		}
/* Computing MAX */
#line 496 "TG01ED.f"
		i__1 = wrkopt, i__2 = ln + (integer) dwork[kw];
#line 496 "TG01ED.f"
		wrkopt = max(i__1,i__2);

/*              Determine the rank of A22. */

#line 500 "TG01ED.f"
		*rnka22 = 0;
#line 501 "TG01ED.f"
		if (dwork[ir1] > svlmax * epsm) {
#line 502 "TG01ED.f"
		    *rnka22 = 1;
#line 503 "TG01ED.f"
		    i__1 = ln;
#line 503 "TG01ED.f"
		    for (i__ = ir1 + 1; i__ <= i__1; ++i__) {
#line 504 "TG01ED.f"
			if (dwork[i__] <= svlmax * toldef) {
#line 504 "TG01ED.f"
			    goto L60;
#line 504 "TG01ED.f"
			}
#line 505 "TG01ED.f"
			++(*rnka22);
#line 506 "TG01ED.f"
/* L50: */
#line 506 "TG01ED.f"
		    }

#line 508 "TG01ED.f"
L60:
#line 509 "TG01ED.f"
		    ;
#line 509 "TG01ED.f"
		}

/*              Apply transformation on the rest of matrices. */

#line 513 "TG01ED.f"
		if (*rnka22 > 0) {

/*                 A <-- diag(I,Q2r') * A * diag(I,Zr2) */

#line 517 "TG01ED.f"
		    dgemm_("Transpose", "No transpose", &ln2, ranke, &ln2, &
			    c_b7, &e[ir1 + ir1 * e_dim1], lde, &a[ir1 + 
			    a_dim1], lda, &c_b6, &e[ir1 + e_dim1], lde, (
			    ftnlen)9, (ftnlen)12);
#line 520 "TG01ED.f"
		    dlacpy_("Full", &ln2, ranke, &e[ir1 + e_dim1], lde, &a[
			    ir1 + a_dim1], lda, (ftnlen)4);
#line 522 "TG01ED.f"
		    dgemm_("No transpose", "Transpose", ranke, &ln2, &ln2, &
			    c_b7, &a[ir1 * a_dim1 + 1], lda, &a[ir1 + ir1 * 
			    a_dim1], lda, &c_b6, &e[ir1 * e_dim1 + 1], lde, (
			    ftnlen)12, (ftnlen)9);
#line 525 "TG01ED.f"
		    dlacpy_("Full", ranke, &ln2, &e[ir1 * e_dim1 + 1], lde, &
			    a[ir1 * a_dim1 + 1], lda, (ftnlen)4);

/*                 B <-- diag(I,Q2r') * B */

#line 530 "TG01ED.f"
		    if (lwr > ln2 * *m && *m > 0) {

#line 532 "TG01ED.f"
			dgemm_("Transpose", "No transpose", &ln2, m, &ln2, &
				c_b7, &e[ir1 + ir1 * e_dim1], lde, &b[ir1 + 
				b_dim1], ldb, &c_b6, &dwork[kw], &ln2, (
				ftnlen)9, (ftnlen)12);
#line 535 "TG01ED.f"
			dlacpy_("Full", &ln2, m, &dwork[kw], &ln2, &b[ir1 + 
				b_dim1], ldb, (ftnlen)4);
#line 537 "TG01ED.f"
		    } else {
#line 538 "TG01ED.f"
			i__1 = *m;
#line 538 "TG01ED.f"
			for (j = 1; j <= i__1; ++j) {
#line 539 "TG01ED.f"
			    dgemv_("Transpose", &ln2, &ln2, &c_b7, &e[ir1 + 
				    ir1 * e_dim1], lde, &b[ir1 + j * b_dim1], 
				    &c__1, &c_b6, &dwork[kw], &c__1, (ftnlen)
				    9);
#line 542 "TG01ED.f"
			    dcopy_(&ln2, &dwork[kw], &c__1, &b[ir1 + j * 
				    b_dim1], &c__1);
#line 543 "TG01ED.f"
/* L70: */
#line 543 "TG01ED.f"
			}
#line 544 "TG01ED.f"
		    }

/*                 C <-- C * diag(I,Zr2) */

#line 548 "TG01ED.f"
		    if (lwr > *p * ln2 && *p > 0) {

#line 550 "TG01ED.f"
			dgemm_("No transpose", "Transpose", p, &ln2, &ln2, &
				c_b7, &c__[ir1 * c_dim1 + 1], ldc, &a[ir1 + 
				ir1 * a_dim1], lda, &c_b6, &dwork[kw], p, (
				ftnlen)12, (ftnlen)9);
#line 553 "TG01ED.f"
			dlacpy_("Full", p, &ln2, &dwork[kw], p, &c__[ir1 * 
				c_dim1 + 1], ldc, (ftnlen)4);
#line 555 "TG01ED.f"
		    } else {
#line 556 "TG01ED.f"
			i__1 = *p;
#line 556 "TG01ED.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 557 "TG01ED.f"
			    dgemv_("No transpose", &ln2, &ln2, &c_b7, &a[ir1 
				    + ir1 * a_dim1], lda, &c__[i__ + ir1 * 
				    c_dim1], ldc, &c_b6, &dwork[kw], &c__1, (
				    ftnlen)12);
#line 560 "TG01ED.f"
			    dcopy_(&ln2, &dwork[kw], &c__1, &c__[i__ + ir1 * 
				    c_dim1], ldc);
#line 561 "TG01ED.f"
/* L80: */
#line 561 "TG01ED.f"
			}
#line 562 "TG01ED.f"
		    }

/*                 Q <-- Q * diag(I, Qr2) */

#line 566 "TG01ED.f"
		    if (lwr > *l * ln2) {

#line 568 "TG01ED.f"
			dgemm_("No transpose", "No transpose", l, &ln2, &ln2, 
				&c_b7, &q[ir1 * q_dim1 + 1], ldq, &e[ir1 + 
				ir1 * e_dim1], lde, &c_b6, &dwork[kw], l, (
				ftnlen)12, (ftnlen)12);
#line 571 "TG01ED.f"
			dlacpy_("Full", l, &ln2, &dwork[kw], l, &q[ir1 * 
				q_dim1 + 1], ldq, (ftnlen)4);
#line 573 "TG01ED.f"
		    } else {
#line 574 "TG01ED.f"
			i__1 = *l;
#line 574 "TG01ED.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 575 "TG01ED.f"
			    dgemv_("Transpose", &ln2, &ln2, &c_b7, &e[ir1 + 
				    ir1 * e_dim1], lde, &q[i__ + ir1 * q_dim1]
				    , ldq, &c_b6, &dwork[kw], &c__1, (ftnlen)
				    9);
#line 578 "TG01ED.f"
			    dcopy_(&ln2, &dwork[kw], &c__1, &q[i__ + ir1 * 
				    q_dim1], ldq);
#line 579 "TG01ED.f"
/* L90: */
#line 579 "TG01ED.f"
			}
#line 580 "TG01ED.f"
		    }

/*                 Z' <-- diag(I, Zr2') * Z' */

#line 584 "TG01ED.f"
		    if (lwr > *n * ln2) {

#line 586 "TG01ED.f"
			dgemm_("No transpose", "No transpose", &ln2, n, &ln2, 
				&c_b7, &a[ir1 + ir1 * a_dim1], lda, &z__[ir1 
				+ z_dim1], ldz, &c_b6, &dwork[kw], &ln2, (
				ftnlen)12, (ftnlen)12);
#line 589 "TG01ED.f"
			dlacpy_("Full", &ln2, n, &dwork[kw], &ln2, &z__[ir1 + 
				z_dim1], ldz, (ftnlen)4);
#line 591 "TG01ED.f"
		    } else {
#line 592 "TG01ED.f"
			i__1 = *n;
#line 592 "TG01ED.f"
			for (j = 1; j <= i__1; ++j) {
#line 593 "TG01ED.f"
			    dgemv_("No transpose", &ln2, &ln2, &c_b7, &a[ir1 
				    + ir1 * a_dim1], lda, &z__[ir1 + j * 
				    z_dim1], &c__1, &c_b6, &dwork[kw], &c__1, 
				    (ftnlen)12);
#line 596 "TG01ED.f"
			    dcopy_(&ln2, &dwork[kw], &c__1, &z__[ir1 + j * 
				    z_dim1], &c__1);
#line 597 "TG01ED.f"
/* L100: */
#line 597 "TG01ED.f"
			}
#line 598 "TG01ED.f"
		    }
#line 599 "TG01ED.f"
		}
#line 600 "TG01ED.f"
	    } else {

/*              Compute the LQ decomposition of A22 in the form */

/*                  A22 = ( L2 0 )* Z2 */

/*              where L2 is lower triangular. */
/*              Workspace: need   MIN(L,N) + L; */
/*                         prefer MIN(L,N) + L*NB. */

#line 610 "TG01ED.f"
		dgelqf_(&la22, &na22, &a[ir1 + ir1 * a_dim1], lda, &dwork[ir1]
			, &dwork[kw], &lwr, info);
/* Computing MAX */
#line 612 "TG01ED.f"
		i__1 = wrkopt, i__2 = ln + (integer) dwork[kw];
#line 612 "TG01ED.f"
		wrkopt = max(i__1,i__2);

/*              Apply transformation Z2 to A, C, and Z. */

/*              A <-- A * diag(I, Z2') */
/*              Workspace: need   2*MIN(L,N); */
/*                         prefer MIN(L,N) + MIN(L,N)*NB. */

#line 620 "TG01ED.f"
		dormlq_("Right", "Transpose", ranke, &na22, &ln2, &a[ir1 + 
			ir1 * a_dim1], lda, &dwork[ir1], &a[ir1 * a_dim1 + 1],
			 lda, &dwork[kw], &lwr, info, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
#line 623 "TG01ED.f"
		i__1 = wrkopt, i__2 = ln + (integer) dwork[kw];
#line 623 "TG01ED.f"
		wrkopt = max(i__1,i__2);

/*              C <-- C * diag(I, Z2') */
/*              Workspace: need   MIN(L,N) + P; */
/*                         prefer MIN(L,N) + P*NB. */

#line 629 "TG01ED.f"
		if (*p > 0) {
#line 630 "TG01ED.f"
		    dormlq_("Right", "Transpose", p, &na22, &ln2, &a[ir1 + 
			    ir1 * a_dim1], lda, &dwork[ir1], &c__[ir1 * 
			    c_dim1 + 1], ldc, &dwork[kw], &lwr, info, (ftnlen)
			    5, (ftnlen)9);
/* Computing MAX */
#line 633 "TG01ED.f"
		    i__1 = wrkopt, i__2 = ln + (integer) dwork[kw];
#line 633 "TG01ED.f"
		    wrkopt = max(i__1,i__2);
#line 634 "TG01ED.f"
		}

/*              Z' <-  diag(I, Z2) * Z' */
/*              Workspace: need   MIN(L,N) + N; */
/*                         prefer MIN(L,N) + N*NB. */

#line 640 "TG01ED.f"
		dormlq_("Left", "No transpose", &na22, n, &ln2, &a[ir1 + ir1 *
			 a_dim1], lda, &dwork[ir1], &z__[ir1 + z_dim1], ldz, &
			dwork[kw], &lwr, info, (ftnlen)4, (ftnlen)12);
/* Computing MAX */
#line 643 "TG01ED.f"
		i__1 = wrkopt, i__2 = ln + (integer) dwork[kw];
#line 643 "TG01ED.f"
		wrkopt = max(i__1,i__2);

/*              Compute the SVD of the lower triangular submatrix L2 as */

/*                              ( Ar  0 ) */
/*                  L2' = Z2r * (       ) * Q2r' */
/*                              (  0  0 ) */

/*              where Q2r' is stored in E and Z2r is stored in A22. */
/*              Workspace: need   MAX(1,5*MIN(L,N)); */
/*                         prefer larger. */

#line 655 "TG01ED.f"
		ma02ad_("Lower", &ln2, &ln2, &a[ir1 + ir1 * a_dim1], lda, &e[
			ir1 + ir1 * e_dim1], lde, (ftnlen)5);
#line 657 "TG01ED.f"
		mb03ud_("Vectors", "Vectors", &ln2, &e[ir1 + ir1 * e_dim1], 
			lde, &a[ir1 + ir1 * a_dim1], lda, &dwork[ir1], &dwork[
			kw], &lwr, info, (ftnlen)7, (ftnlen)7);
#line 660 "TG01ED.f"
		if (*info > 0) {
#line 660 "TG01ED.f"
		    return 0;
#line 660 "TG01ED.f"
		}
/* Computing MAX */
#line 662 "TG01ED.f"
		i__1 = wrkopt, i__2 = ln + (integer) dwork[kw];
#line 662 "TG01ED.f"
		wrkopt = max(i__1,i__2);

/*              Determine the rank of A22. */

#line 666 "TG01ED.f"
		*rnka22 = 0;
#line 667 "TG01ED.f"
		if (dwork[ir1] > svlmax * epsm) {
#line 668 "TG01ED.f"
		    *rnka22 = 1;
#line 669 "TG01ED.f"
		    i__1 = ln;
#line 669 "TG01ED.f"
		    for (i__ = ir1 + 1; i__ <= i__1; ++i__) {
#line 670 "TG01ED.f"
			if (dwork[i__] <= svlmax * toldef) {
#line 670 "TG01ED.f"
			    goto L120;
#line 670 "TG01ED.f"
			}
#line 671 "TG01ED.f"
			++(*rnka22);
#line 672 "TG01ED.f"
/* L110: */
#line 672 "TG01ED.f"
		    }

#line 674 "TG01ED.f"
L120:
#line 675 "TG01ED.f"
		    ;
#line 675 "TG01ED.f"
		}

/*              Apply transformation on the rest of matrices. */

#line 679 "TG01ED.f"
		if (*rnka22 > 0) {

/*                 A <-- diag(I,Q2r') * A * diag(I,Zr2) */

#line 683 "TG01ED.f"
		    dgemm_("No transpose", "No transpose", &ln2, ranke, &ln2, 
			    &c_b7, &e[ir1 + ir1 * e_dim1], lde, &a[ir1 + 
			    a_dim1], lda, &c_b6, &e[ir1 + e_dim1], lde, (
			    ftnlen)12, (ftnlen)12);
#line 686 "TG01ED.f"
		    dlacpy_("Full", &ln2, ranke, &e[ir1 + e_dim1], lde, &a[
			    ir1 + a_dim1], lda, (ftnlen)4);
#line 688 "TG01ED.f"
		    dgemm_("No transpose", "No transpose", ranke, &ln2, &ln2, 
			    &c_b7, &a[ir1 * a_dim1 + 1], lda, &a[ir1 + ir1 * 
			    a_dim1], lda, &c_b6, &e[ir1 * e_dim1 + 1], lde, (
			    ftnlen)12, (ftnlen)12);
#line 691 "TG01ED.f"
		    dlacpy_("Full", ranke, &ln2, &e[ir1 * e_dim1 + 1], lde, &
			    a[ir1 * a_dim1 + 1], lda, (ftnlen)4);

/*                 B <-- diag(I,Q2r') * B */

#line 696 "TG01ED.f"
		    if (lwr > ln2 * *m && *m > 0) {

#line 698 "TG01ED.f"
			dgemm_("No transpose", "No transpose", &ln2, m, &ln2, 
				&c_b7, &e[ir1 + ir1 * e_dim1], lde, &b[ir1 + 
				b_dim1], ldb, &c_b6, &dwork[kw], &ln2, (
				ftnlen)12, (ftnlen)12);
#line 701 "TG01ED.f"
			dlacpy_("Full", &ln2, m, &dwork[kw], &ln2, &b[ir1 + 
				b_dim1], ldb, (ftnlen)4);
#line 703 "TG01ED.f"
		    } else {
#line 704 "TG01ED.f"
			i__1 = *m;
#line 704 "TG01ED.f"
			for (j = 1; j <= i__1; ++j) {
#line 705 "TG01ED.f"
			    dgemv_("No transpose", &ln2, &ln2, &c_b7, &e[ir1 
				    + ir1 * e_dim1], lde, &b[ir1 + j * b_dim1]
				    , &c__1, &c_b6, &dwork[kw], &c__1, (
				    ftnlen)12);
#line 708 "TG01ED.f"
			    dcopy_(&ln2, &dwork[kw], &c__1, &b[ir1 + j * 
				    b_dim1], &c__1);
#line 709 "TG01ED.f"
/* L130: */
#line 709 "TG01ED.f"
			}
#line 710 "TG01ED.f"
		    }

/*                 C <-- C * diag(I,Zr2) */

#line 714 "TG01ED.f"
		    if (lwr > *p * ln2 && *p > 0) {

#line 716 "TG01ED.f"
			dgemm_("No transpose", "No transpose", p, &ln2, &ln2, 
				&c_b7, &c__[ir1 * c_dim1 + 1], ldc, &a[ir1 + 
				ir1 * a_dim1], lda, &c_b6, &dwork[kw], p, (
				ftnlen)12, (ftnlen)12);
#line 719 "TG01ED.f"
			dlacpy_("Full", p, &ln2, &dwork[kw], p, &c__[ir1 * 
				c_dim1 + 1], ldc, (ftnlen)4);
#line 721 "TG01ED.f"
		    } else {
#line 722 "TG01ED.f"
			i__1 = *p;
#line 722 "TG01ED.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 723 "TG01ED.f"
			    dgemv_("Transpose", &ln2, &ln2, &c_b7, &a[ir1 + 
				    ir1 * a_dim1], lda, &c__[i__ + ir1 * 
				    c_dim1], ldc, &c_b6, &dwork[kw], &c__1, (
				    ftnlen)9);
#line 726 "TG01ED.f"
			    dcopy_(&ln2, &dwork[kw], &c__1, &c__[i__ + ir1 * 
				    c_dim1], ldc);
#line 727 "TG01ED.f"
/* L140: */
#line 727 "TG01ED.f"
			}
#line 728 "TG01ED.f"
		    }

/*                 Q <-- Q * diag(I, Qr2) */

#line 732 "TG01ED.f"
		    if (lwr > *l * ln2) {

#line 734 "TG01ED.f"
			dgemm_("No transpose", "Transpose", l, &ln2, &ln2, &
				c_b7, &q[ir1 * q_dim1 + 1], ldq, &e[ir1 + ir1 
				* e_dim1], lde, &c_b6, &dwork[kw], l, (ftnlen)
				12, (ftnlen)9);
#line 737 "TG01ED.f"
			dlacpy_("Full", l, &ln2, &dwork[kw], l, &q[ir1 * 
				q_dim1 + 1], ldq, (ftnlen)4);
#line 739 "TG01ED.f"
		    } else {
#line 740 "TG01ED.f"
			i__1 = *l;
#line 740 "TG01ED.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 741 "TG01ED.f"
			    dgemv_("No transpose", &ln2, &ln2, &c_b7, &e[ir1 
				    + ir1 * e_dim1], lde, &q[i__ + ir1 * 
				    q_dim1], ldq, &c_b6, &dwork[kw], &c__1, (
				    ftnlen)12);
#line 744 "TG01ED.f"
			    dcopy_(&ln2, &dwork[kw], &c__1, &q[i__ + ir1 * 
				    q_dim1], ldq);
#line 745 "TG01ED.f"
/* L150: */
#line 745 "TG01ED.f"
			}
#line 746 "TG01ED.f"
		    }

/*                 Z' <-- diag(I, Zr2') * Z' */

#line 750 "TG01ED.f"
		    if (lwr > *n * ln2) {

#line 752 "TG01ED.f"
			dgemm_("Transpose", "No transpose", &ln2, n, &ln2, &
				c_b7, &a[ir1 + ir1 * a_dim1], lda, &z__[ir1 + 
				z_dim1], ldz, &c_b6, &dwork[kw], &ln2, (
				ftnlen)9, (ftnlen)12);
#line 755 "TG01ED.f"
			dlacpy_("Full", &ln2, n, &dwork[kw], &ln2, &z__[ir1 + 
				z_dim1], ldz, (ftnlen)4);
#line 757 "TG01ED.f"
		    } else {
#line 758 "TG01ED.f"
			i__1 = *n;
#line 758 "TG01ED.f"
			for (j = 1; j <= i__1; ++j) {
#line 759 "TG01ED.f"
			    dgemv_("Transpose", &ln2, &ln2, &c_b7, &a[ir1 + 
				    ir1 * a_dim1], lda, &z__[ir1 + j * z_dim1]
				    , &c__1, &c_b6, &dwork[kw], &c__1, (
				    ftnlen)9);
#line 762 "TG01ED.f"
			    dcopy_(&ln2, &dwork[kw], &c__1, &z__[ir1 + j * 
				    z_dim1], &c__1);
#line 763 "TG01ED.f"
/* L160: */
#line 763 "TG01ED.f"
			}
#line 764 "TG01ED.f"
		    }
#line 765 "TG01ED.f"
		}
#line 766 "TG01ED.f"
	    }
#line 767 "TG01ED.f"
	}
#line 768 "TG01ED.f"
    }

/*     Set E. */

#line 772 "TG01ED.f"
    dlaset_("Full", l, n, &c_b6, &c_b6, &e[e_offset], lde, (ftnlen)4);
#line 773 "TG01ED.f"
    i__1 = *lde + 1;
#line 773 "TG01ED.f"
    dcopy_(ranke, &dwork[1], &c__1, &e[e_offset], &i__1);

#line 775 "TG01ED.f"
    if (reda) {

/*        Set A22. */

#line 779 "TG01ED.f"
	dlaset_("Full", &la22, &na22, &c_b6, &c_b6, &a[ir1 + ir1 * a_dim1], 
		lda, (ftnlen)4);
#line 780 "TG01ED.f"
	i__1 = *lda + 1;
#line 780 "TG01ED.f"
	dcopy_(rnka22, &dwork[ir1], &c__1, &a[ir1 + ir1 * a_dim1], &i__1);
#line 781 "TG01ED.f"
    }

/*     Transpose Z. */

#line 785 "TG01ED.f"
    i__1 = *n;
#line 785 "TG01ED.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 786 "TG01ED.f"
	i__2 = i__ - 1;
#line 786 "TG01ED.f"
	dswap_(&i__2, &z__[i__ * z_dim1 + 1], &c__1, &z__[i__ + z_dim1], ldz);
#line 787 "TG01ED.f"
/* L170: */
#line 787 "TG01ED.f"
    }

#line 789 "TG01ED.f"
    dwork[1] = (doublereal) wrkopt;

#line 791 "TG01ED.f"
    return 0;
/* *** Last line of TG01ED *** */
} /* tg01ed_ */

