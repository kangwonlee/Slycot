#line 1 "TB05AD.f"
/* TB05AD.f -- translated by f2c (version 20100827).
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

#line 1 "TB05AD.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static integer c__1 = 1;

/* Subroutine */ int tb05ad_(char *baleig, char *inita, integer *n, integer *
	m, integer *p, doublecomplex *freq, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *rcond, doublecomplex *g, integer *ldg, doublereal *evre, 
	doublereal *evim, doublecomplex *hinvb, integer *ldhinv, integer *
	iwork, doublereal *dwork, integer *ldwork, doublecomplex *zwork, 
	integer *lzwork, integer *info, ftnlen baleig_len, ftnlen inita_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, g_dim1, 
	    g_offset, hinvb_dim1, hinvb_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6, i__7;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal t;
    static integer ij, jj, jp, igh, low, itau;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern doublereal dasum_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int mb02rz_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen), mb02sz_(integer *, doublecomplex *, integer *
	    , integer *, integer *), dswap_(integer *, doublereal *, integer *
	    , doublereal *, integer *), mb02tz_(char *, integer *, doublereal 
	    *, doublecomplex *, integer *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, ftnlen);
    static doublereal hnorm;
    static integer jwork;
    static logical lbalba;
    extern /* Subroutine */ int dgebal_(char *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, ftnlen);
    static char balanc[1];
    static logical lbalea, lbaleb, lbalec;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dgehrd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), xerbla_(char *, integer *, ftnlen);
    static logical linita;
    extern /* Subroutine */ int dhseqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dormhr_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen);
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

/*     To find the complex frequency response matrix (transfer matrix) */
/*     G(freq) of the state-space representation (A,B,C) given by */
/*                                   -1 */
/*        G(freq) = C * ((freq*I - A)  ) * B */

/*     where A, B and C are real N-by-N, N-by-M and P-by-N matrices */
/*     respectively and freq is a complex scalar. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     BALEIG  CHARACTER*1 */
/*             Determines whether the user wishes to balance matrix A */
/*             and/or compute its eigenvalues and/or estimate the */
/*             condition number of the problem as follows: */
/*             = 'N':  The matrix A should not be balanced and neither */
/*                     the eigenvalues of A nor the condition number */
/*                     estimate of the problem are to be calculated; */
/*             = 'C':  The matrix A should not be balanced and only an */
/*                     estimate of the condition number of the problem */
/*                     is to be calculated; */
/*             = 'B' or 'E' and INITA = 'G':  The matrix A is to be */
/*                     balanced and its eigenvalues calculated; */
/*             = 'A' and INITA = 'G':  The matrix A is to be balanced, */
/*                     and its eigenvalues and an estimate of the */
/*                     condition number of the problem are to be */
/*                     calculated. */

/*     INITA   CHARACTER*1 */
/*             Specifies whether or not the matrix A is already in upper */
/*             Hessenberg form as follows: */
/*             = 'G':  The matrix A is a general matrix; */
/*             = 'H':  The matrix A is in upper Hessenberg form and */
/*                     neither balancing nor the eigenvalues of A are */
/*                     required. */
/*             INITA must be set to 'G' for the first call to the */
/*             routine, unless the matrix A is already in upper */
/*             Hessenberg form and neither balancing nor the eigenvalues */
/*             of A are required. Thereafter, it must be set to 'H' for */
/*             all subsequent calls. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of states, i.e. the order of the state */
/*             transition matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of inputs, i.e. the number of columns in the */
/*             matrix B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of outputs, i.e. the number of rows in the */
/*             matrix C.  P >= 0. */

/*     FREQ    (input) COMPLEX*16 */
/*             The frequency freq at which the frequency response matrix */
/*             (transfer matrix) is to be evaluated. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state transition matrix A. */
/*             If INITA = 'G', then, on exit, the leading N-by-N part of */
/*             this array contains an upper Hessenberg matrix similar to */
/*             (via an orthogonal matrix consisting of a sequence of */
/*             Householder transformations) the original state transition */
/*             matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input/state matrix B. */
/*             If INITA = 'G', then, on exit, the leading N-by-M part of */
/*             this array contains the product of the transpose of the */
/*             orthogonal transformation matrix used to reduce A to upper */
/*             Hessenberg form and the original input/state matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C. */
/*             If INITA = 'G', then, on exit, the leading P-by-N part of */
/*             this array contains the product of the original output/ */
/*             state matrix C and the orthogonal transformation matrix */
/*             used to reduce A to upper Hessenberg form. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     RCOND   (output) DOUBLE PRECISION */
/*             If BALEIG = 'C' or BALEIG = 'A', then RCOND contains an */
/*             estimate of the reciprocal of the condition number of */
/*             matrix H with respect to inversion (see METHOD). */

/*     G       (output) COMPLEX*16 array, dimension (LDG,M) */
/*             The leading P-by-M part of this array contains the */
/*             frequency response matrix G(freq). */

/*     LDG     INTEGER */
/*             The leading dimension of array G.  LDG >= MAX(1,P). */

/*     EVRE,   (output) DOUBLE PRECISION arrays, dimension (N) */
/*     EVIM    If INITA = 'G' and BALEIG = 'B' or 'E' or BALEIG = 'A', */
/*             then these arrays contain the real and imaginary parts, */
/*             respectively, of the eigenvalues of the matrix A. */
/*             Otherwise, these arrays are not referenced. */

/*     HINVB   (output) COMPLEX*16 array, dimension (LDHINV,M) */
/*             The leading N-by-M part of this array contains the */
/*                      -1 */
/*             product H  B. */

/*     LDHINV  INTEGER */
/*             The leading dimension of array HINVB.  LDHINV >= MAX(1,N). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1, N - 1 + MAX(N,M,P)), */
/*                       if INITA = 'G' and BALEIG = 'N', or 'B', or 'E'; */
/*             LDWORK >= MAX(1, N + MAX(N,M-1,P-1)), */
/*                       if INITA = 'G' and BALEIG = 'C', or 'A'; */
/*             LDWORK >= MAX(1, 2*N), */
/*                       if INITA = 'H' and BALEIG = 'C', or 'A'; */
/*             LDWORK >= 1, otherwise. */
/*             For optimum performance when INITA = 'G' LDWORK should be */
/*             larger. */

/*     ZWORK   COMPLEX*16 array, dimension (LZWORK) */

/*     LZWORK  INTEGER */
/*             The length of the array ZWORK. */
/*             LZWORK >= MAX(1,N*N+2*N), if BALEIG = 'C', or 'A'; */
/*             LZWORK >= MAX(1,N*N),     otherwise. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if more than 30*N iterations are required to */
/*                   isolate all the eigenvalues of the matrix A; the */
/*                   computations are continued; */
/*             = 2:  if either FREQ is too near to an eigenvalue of the */
/*                   matrix A, or RCOND is less than EPS, where EPS is */
/*                   the machine  precision (see LAPACK Library routine */
/*                   DLAMCH). */

/*     METHOD */

/*     The matrix A is first balanced (if BALEIG = 'B' or 'E', or */
/*     BALEIG = 'A') and then reduced to upper Hessenberg form; the same */
/*     transformations are applied to the matrix B and the matrix C. */
/*     The complex Hessenberg matrix  H = (freq*I - A) is then used */
/*                       -1 */
/*     to solve for C * H  * B. */

/*     Depending on the input values of parameters BALEIG and INITA, */
/*     the eigenvalues of matrix A and the condition number of */
/*     matrix H with respect to inversion are also calculated. */

/*     REFERENCES */

/*     [1] Laub, A.J. */
/*         Efficient Calculation of Frequency Response Matrices from */
/*         State-Space Models. */
/*         ACM TOMS, 12, pp. 26-33, 1986. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TB01FD by A.J.Laub, University of */
/*     Southern California, Los Angeles, CA 90089, United States of */
/*     America, June 1982. */

/*     REVISIONS */

/*     V. Sima, February 22, 1998 (changed the name of TB01RD). */
/*     V. Sima, February 12, 1999, August 7, 2003. */
/*     A. Markovski, Technical University of Sofia, September 30, 2003. */
/*     V. Sima, October 1, 2003. */

/*     KEYWORDS */

/*     Frequency response, Hessenberg form, matrix algebra, input output */
/*     description, multivariable system, orthogonal transformation, */
/*     similarity transformation, state-space representation, transfer */
/*     matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 267 "TB05AD.f"
    /* Parameter adjustments */
#line 267 "TB05AD.f"
    a_dim1 = *lda;
#line 267 "TB05AD.f"
    a_offset = 1 + a_dim1;
#line 267 "TB05AD.f"
    a -= a_offset;
#line 267 "TB05AD.f"
    b_dim1 = *ldb;
#line 267 "TB05AD.f"
    b_offset = 1 + b_dim1;
#line 267 "TB05AD.f"
    b -= b_offset;
#line 267 "TB05AD.f"
    c_dim1 = *ldc;
#line 267 "TB05AD.f"
    c_offset = 1 + c_dim1;
#line 267 "TB05AD.f"
    c__ -= c_offset;
#line 267 "TB05AD.f"
    g_dim1 = *ldg;
#line 267 "TB05AD.f"
    g_offset = 1 + g_dim1;
#line 267 "TB05AD.f"
    g -= g_offset;
#line 267 "TB05AD.f"
    --evre;
#line 267 "TB05AD.f"
    --evim;
#line 267 "TB05AD.f"
    hinvb_dim1 = *ldhinv;
#line 267 "TB05AD.f"
    hinvb_offset = 1 + hinvb_dim1;
#line 267 "TB05AD.f"
    hinvb -= hinvb_offset;
#line 267 "TB05AD.f"
    --iwork;
#line 267 "TB05AD.f"
    --dwork;
#line 267 "TB05AD.f"
    --zwork;
#line 267 "TB05AD.f"

#line 267 "TB05AD.f"
    /* Function Body */
#line 267 "TB05AD.f"
    *info = 0;
#line 268 "TB05AD.f"
    lbalec = lsame_(baleig, "C", (ftnlen)1, (ftnlen)1);
#line 269 "TB05AD.f"
    lbaleb = lsame_(baleig, "B", (ftnlen)1, (ftnlen)1) || lsame_(baleig, 
	    "E", (ftnlen)1, (ftnlen)1);
#line 270 "TB05AD.f"
    lbalea = lsame_(baleig, "A", (ftnlen)1, (ftnlen)1);
#line 271 "TB05AD.f"
    lbalba = lbaleb || lbalea;
#line 272 "TB05AD.f"
    linita = lsame_(inita, "G", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 276 "TB05AD.f"
    if (! lbalec && ! lbalba && ! lsame_(baleig, "N", (ftnlen)1, (ftnlen)1)) {
#line 278 "TB05AD.f"
	*info = -1;
#line 279 "TB05AD.f"
    } else if (! linita && ! lsame_(inita, "H", (ftnlen)1, (ftnlen)1)) {
#line 280 "TB05AD.f"
	*info = -2;
#line 281 "TB05AD.f"
    } else if (*n < 0) {
#line 282 "TB05AD.f"
	*info = -3;
#line 283 "TB05AD.f"
    } else if (*m < 0) {
#line 284 "TB05AD.f"
	*info = -4;
#line 285 "TB05AD.f"
    } else if (*p < 0) {
#line 286 "TB05AD.f"
	*info = -5;
#line 287 "TB05AD.f"
    } else if (*lda < max(1,*n)) {
#line 288 "TB05AD.f"
	*info = -8;
#line 289 "TB05AD.f"
    } else if (*ldb < max(1,*n)) {
#line 290 "TB05AD.f"
	*info = -10;
#line 291 "TB05AD.f"
    } else if (*ldc < max(1,*p)) {
#line 292 "TB05AD.f"
	*info = -12;
#line 293 "TB05AD.f"
    } else if (*ldg < max(1,*p)) {
#line 294 "TB05AD.f"
	*info = -15;
#line 295 "TB05AD.f"
    } else if (*ldhinv < max(1,*n)) {
#line 296 "TB05AD.f"
	*info = -19;
#line 297 "TB05AD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 297 "TB05AD.f"
	i__1 = max(*n,*m);
/* Computing MAX */
#line 297 "TB05AD.f"
	i__2 = *n, i__3 = *m - 1, i__2 = max(i__2,i__3), i__3 = *p - 1;
#line 297 "TB05AD.f"
	if (linita && ! lbalec && ! lbalea && *ldwork < *n - 1 + max(i__1,*p) 
		|| linita && (lbalec || lbalea) && *ldwork < *n + max(i__2,
		i__3) || ! linita && (lbalec || lbalea) && *ldwork < *n << 1 
		|| *ldwork < 1) {
#line 303 "TB05AD.f"
	    *info = -22;
#line 304 "TB05AD.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 304 "TB05AD.f"
	    i__1 = 1, i__2 = *n * *n;
#line 304 "TB05AD.f"
	    if ((lbalec || lbalea) && *lzwork < *n * (*n + 2) || *lzwork < 
		    max(i__1,i__2)) {
#line 306 "TB05AD.f"
		*info = -24;
#line 307 "TB05AD.f"
	    }
#line 307 "TB05AD.f"
	}
#line 307 "TB05AD.f"
    }

#line 309 "TB05AD.f"
    if (*info != 0) {

/*        Error return */

#line 313 "TB05AD.f"
	i__1 = -(*info);
#line 313 "TB05AD.f"
	xerbla_("TB05AD", &i__1, (ftnlen)6);
#line 314 "TB05AD.f"
	return 0;
#line 315 "TB05AD.f"
    }

/*     Quick return if possible. */

#line 319 "TB05AD.f"
    if (*n == 0) {
#line 320 "TB05AD.f"
	if (min(*m,*p) > 0) {
#line 320 "TB05AD.f"
	    zlaset_("Full", p, m, &c_b1, &c_b1, &g[g_offset], ldg, (ftnlen)4);
#line 320 "TB05AD.f"
	}
#line 322 "TB05AD.f"
	*rcond = 1.;
#line 323 "TB05AD.f"
	dwork[1] = 1.;
#line 324 "TB05AD.f"
	return 0;
#line 325 "TB05AD.f"
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 333 "TB05AD.f"
    wrkopt = 1;

#line 335 "TB05AD.f"
    if (linita) {
#line 336 "TB05AD.f"
	*(unsigned char *)balanc = 'N';
#line 337 "TB05AD.f"
	if (lbalba) {
#line 337 "TB05AD.f"
	    *(unsigned char *)balanc = 'B';
#line 337 "TB05AD.f"
	}

/*        Workspace: need N. */

#line 341 "TB05AD.f"
	dgebal_(balanc, n, &a[a_offset], lda, &low, &igh, &dwork[1], info, (
		ftnlen)1);
#line 342 "TB05AD.f"
	if (lbalba) {

/*           Adjust B and C matrices based on information in the */
/*           vector DWORK which describes the balancing of A and is */
/*           defined in the subroutine DGEBAL. */

#line 348 "TB05AD.f"
	    i__1 = *n;
#line 348 "TB05AD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 349 "TB05AD.f"
		jj = j;
#line 350 "TB05AD.f"
		if (jj < low || jj > igh) {
#line 351 "TB05AD.f"
		    if (jj < low) {
#line 351 "TB05AD.f"
			jj = low - jj;
#line 351 "TB05AD.f"
		    }
#line 352 "TB05AD.f"
		    jp = (integer) dwork[jj];
#line 353 "TB05AD.f"
		    if (jp != jj) {

/*                    Permute rows of B. */

#line 357 "TB05AD.f"
			if (*m > 0) {
#line 357 "TB05AD.f"
			    dswap_(m, &b[jj + b_dim1], ldb, &b[jp + b_dim1], 
				    ldb);
#line 357 "TB05AD.f"
			}

/*                    Permute columns of C. */

#line 362 "TB05AD.f"
			if (*p > 0) {
#line 362 "TB05AD.f"
			    dswap_(p, &c__[jj * c_dim1 + 1], &c__1, &c__[jp * 
				    c_dim1 + 1], &c__1);
#line 362 "TB05AD.f"
			}
#line 364 "TB05AD.f"
		    }
#line 365 "TB05AD.f"
		}
#line 366 "TB05AD.f"
/* L10: */
#line 366 "TB05AD.f"
	    }

#line 368 "TB05AD.f"
	    if (igh != low) {

#line 370 "TB05AD.f"
		i__1 = igh;
#line 370 "TB05AD.f"
		for (j = low; j <= i__1; ++j) {
#line 371 "TB05AD.f"
		    t = dwork[j];

/*                 Scale rows of permuted B. */

#line 375 "TB05AD.f"
		    if (*m > 0) {
#line 375 "TB05AD.f"
			d__1 = 1. / t;
#line 375 "TB05AD.f"
			dscal_(m, &d__1, &b[j + b_dim1], ldb);
#line 375 "TB05AD.f"
		    }

/*                 Scale columns of permuted C. */

#line 380 "TB05AD.f"
		    if (*p > 0) {
#line 380 "TB05AD.f"
			dscal_(p, &t, &c__[j * c_dim1 + 1], &c__1);
#line 380 "TB05AD.f"
		    }
#line 382 "TB05AD.f"
/* L20: */
#line 382 "TB05AD.f"
		}

#line 384 "TB05AD.f"
	    }
#line 385 "TB05AD.f"
	}

/*        Reduce A to Hessenberg form by orthogonal similarities and */
/*        accumulate the orthogonal transformations into B and C. */
/*        Workspace: need 2*N - 1;  prefer N - 1 + N*NB. */

#line 391 "TB05AD.f"
	itau = 1;
#line 392 "TB05AD.f"
	jwork = itau + *n - 1;
#line 393 "TB05AD.f"
	i__1 = *ldwork - jwork + 1;
#line 393 "TB05AD.f"
	dgehrd_(n, &low, &igh, &a[a_offset], lda, &dwork[itau], &dwork[jwork],
		 &i__1, info);
/* Computing MAX */
#line 395 "TB05AD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 395 "TB05AD.f"
	wrkopt = max(i__1,i__2);

/*        Workspace: need N - 1 + M;  prefer N - 1 + M*NB. */

#line 399 "TB05AD.f"
	i__1 = *ldwork - jwork + 1;
#line 399 "TB05AD.f"
	dormhr_("Left", "Transpose", n, m, &low, &igh, &a[a_offset], lda, &
		dwork[itau], &b[b_offset], ldb, &dwork[jwork], &i__1, info, (
		ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 402 "TB05AD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 402 "TB05AD.f"
	wrkopt = max(i__1,i__2);

/*        Workspace: need N - 1 + P;  prefer N - 1 + P*NB. */

#line 406 "TB05AD.f"
	i__1 = *ldwork - jwork + 1;
#line 406 "TB05AD.f"
	dormhr_("Right", "No transpose", p, n, &low, &igh, &a[a_offset], lda, 
		&dwork[itau], &c__[c_offset], ldc, &dwork[jwork], &i__1, info,
		 (ftnlen)5, (ftnlen)12);
/* Computing MAX */
#line 409 "TB05AD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 409 "TB05AD.f"
	wrkopt = max(i__1,i__2);
#line 410 "TB05AD.f"
	if (lbalba) {

/*           Temporarily store Hessenberg form of A in array ZWORK. */

#line 414 "TB05AD.f"
	    ij = 0;
#line 415 "TB05AD.f"
	    i__1 = *n;
#line 415 "TB05AD.f"
	    for (j = 1; j <= i__1; ++j) {

#line 417 "TB05AD.f"
		i__2 = *n;
#line 417 "TB05AD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 418 "TB05AD.f"
		    ++ij;
#line 419 "TB05AD.f"
		    i__3 = ij;
#line 419 "TB05AD.f"
		    i__4 = i__ + j * a_dim1;
#line 419 "TB05AD.f"
		    z__1.r = a[i__4], z__1.i = 0.;
#line 419 "TB05AD.f"
		    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 420 "TB05AD.f"
/* L30: */
#line 420 "TB05AD.f"
		}

#line 422 "TB05AD.f"
/* L40: */
#line 422 "TB05AD.f"
	    }

/*           Compute the eigenvalues of A if that option is requested. */
/*           Workspace: need N. */

#line 427 "TB05AD.f"
	    dhseqr_("Eigenvalues", "No Schur", n, &low, &igh, &a[a_offset], 
		    lda, &evre[1], &evim[1], &dwork[1], &c__1, &dwork[1], 
		    ldwork, info, (ftnlen)11, (ftnlen)8);

/*           Restore upper Hessenberg form of A. */

#line 432 "TB05AD.f"
	    ij = 0;
#line 433 "TB05AD.f"
	    i__1 = *n;
#line 433 "TB05AD.f"
	    for (j = 1; j <= i__1; ++j) {

#line 435 "TB05AD.f"
		i__2 = *n;
#line 435 "TB05AD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 436 "TB05AD.f"
		    ++ij;
#line 437 "TB05AD.f"
		    i__3 = ij;
#line 437 "TB05AD.f"
		    a[i__ + j * a_dim1] = zwork[i__3].r;
#line 438 "TB05AD.f"
/* L50: */
#line 438 "TB05AD.f"
		}

#line 440 "TB05AD.f"
/* L60: */
#line 440 "TB05AD.f"
	    }

#line 442 "TB05AD.f"
	    if (*info > 0) {

/*              DHSEQR could not evaluate the eigenvalues of A. */

#line 446 "TB05AD.f"
		*info = 1;
#line 447 "TB05AD.f"
	    }
#line 448 "TB05AD.f"
	}
#line 449 "TB05AD.f"
    }

/*     Update  H := (FREQ * I) - A   with appropriate value of FREQ. */

#line 453 "TB05AD.f"
    ij = 0;
#line 454 "TB05AD.f"
    jj = 1;
#line 455 "TB05AD.f"
    i__1 = *n;
#line 455 "TB05AD.f"
    for (j = 1; j <= i__1; ++j) {

#line 457 "TB05AD.f"
	i__2 = *n;
#line 457 "TB05AD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 458 "TB05AD.f"
	    ++ij;
#line 459 "TB05AD.f"
	    i__3 = ij;
#line 459 "TB05AD.f"
	    i__4 = i__ + j * a_dim1;
#line 459 "TB05AD.f"
	    z__2.r = a[i__4], z__2.i = 0.;
#line 459 "TB05AD.f"
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
#line 459 "TB05AD.f"
	    zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 460 "TB05AD.f"
/* L70: */
#line 460 "TB05AD.f"
	}

#line 462 "TB05AD.f"
	i__2 = jj;
#line 462 "TB05AD.f"
	i__3 = jj;
#line 462 "TB05AD.f"
	z__1.r = freq->r + zwork[i__3].r, z__1.i = freq->i + zwork[i__3].i;
#line 462 "TB05AD.f"
	zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
#line 463 "TB05AD.f"
	jj = jj + *n + 1;
#line 464 "TB05AD.f"
/* L80: */
#line 464 "TB05AD.f"
    }

#line 466 "TB05AD.f"
    if (lbalec || lbalea) {

/*        Efficiently compute the 1-norm of the matrix for condition */
/*        estimation. */

#line 471 "TB05AD.f"
	hnorm = 0.;
#line 472 "TB05AD.f"
	jj = 1;

#line 474 "TB05AD.f"
	i__1 = *n;
#line 474 "TB05AD.f"
	for (j = 1; j <= i__1; ++j) {
#line 475 "TB05AD.f"
	    i__2 = j - 1;
#line 475 "TB05AD.f"
	    t = z_abs(&zwork[jj]) + dasum_(&i__2, &a[j * a_dim1 + 1], &c__1);
#line 476 "TB05AD.f"
	    if (j < *n) {
#line 476 "TB05AD.f"
		t += (d__1 = a[j + 1 + j * a_dim1], abs(d__1));
#line 476 "TB05AD.f"
	    }
#line 477 "TB05AD.f"
	    hnorm = max(hnorm,t);
#line 478 "TB05AD.f"
	    jj = jj + *n + 1;
#line 479 "TB05AD.f"
/* L90: */
#line 479 "TB05AD.f"
	}

#line 481 "TB05AD.f"
    }

/*     Factor the complex Hessenberg matrix. */

#line 485 "TB05AD.f"
    mb02sz_(n, &zwork[1], n, &iwork[1], info);
#line 486 "TB05AD.f"
    if (*info != 0) {
#line 486 "TB05AD.f"
	*info = 2;
#line 486 "TB05AD.f"
    }

#line 488 "TB05AD.f"
    if (lbalec || lbalea) {

/*        Estimate the condition of the matrix. */

/*        Workspace: need 2*N. */

#line 494 "TB05AD.f"
	mb02tz_("1-norm", n, &hnorm, &zwork[1], n, &iwork[1], rcond, &dwork[1]
		, &zwork[*n * *n + 1], info, (ftnlen)6);
/* Computing MAX */
#line 496 "TB05AD.f"
	i__1 = wrkopt, i__2 = *n << 1;
#line 496 "TB05AD.f"
	wrkopt = max(i__1,i__2);
#line 497 "TB05AD.f"
	if (*rcond < dlamch_("Epsilon", (ftnlen)7)) {
#line 497 "TB05AD.f"
	    *info = 2;
#line 497 "TB05AD.f"
	}
#line 498 "TB05AD.f"
    }

#line 500 "TB05AD.f"
    if (*info != 0) {

/*        Error return: Linear system is numerically or exactly singular. */

#line 504 "TB05AD.f"
	return 0;
#line 505 "TB05AD.f"
    }

/*     Compute  (H-INVERSE)*B. */

#line 509 "TB05AD.f"
    i__1 = *m;
#line 509 "TB05AD.f"
    for (j = 1; j <= i__1; ++j) {

#line 511 "TB05AD.f"
	i__2 = *n;
#line 511 "TB05AD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 512 "TB05AD.f"
	    i__3 = i__ + j * hinvb_dim1;
#line 512 "TB05AD.f"
	    i__4 = i__ + j * b_dim1;
#line 512 "TB05AD.f"
	    z__1.r = b[i__4], z__1.i = 0.;
#line 512 "TB05AD.f"
	    hinvb[i__3].r = z__1.r, hinvb[i__3].i = z__1.i;
#line 513 "TB05AD.f"
/* L100: */
#line 513 "TB05AD.f"
	}

#line 515 "TB05AD.f"
/* L110: */
#line 515 "TB05AD.f"
    }

#line 517 "TB05AD.f"
    mb02rz_("No transpose", n, m, &zwork[1], n, &iwork[1], &hinvb[
	    hinvb_offset], ldhinv, info, (ftnlen)12);

/*     Compute  C*(H-INVERSE)*B. */

#line 522 "TB05AD.f"
    i__1 = *m;
#line 522 "TB05AD.f"
    for (j = 1; j <= i__1; ++j) {

#line 524 "TB05AD.f"
	i__2 = *p;
#line 524 "TB05AD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 525 "TB05AD.f"
	    i__3 = i__ + j * g_dim1;
#line 525 "TB05AD.f"
	    g[i__3].r = 0., g[i__3].i = 0.;
#line 526 "TB05AD.f"
/* L120: */
#line 526 "TB05AD.f"
	}

#line 528 "TB05AD.f"
	i__2 = *n;
#line 528 "TB05AD.f"
	for (k = 1; k <= i__2; ++k) {

#line 530 "TB05AD.f"
	    i__3 = *p;
#line 530 "TB05AD.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 531 "TB05AD.f"
		i__4 = i__ + j * g_dim1;
#line 531 "TB05AD.f"
		i__5 = i__ + j * g_dim1;
#line 531 "TB05AD.f"
		i__6 = i__ + k * c_dim1;
#line 531 "TB05AD.f"
		z__3.r = c__[i__6], z__3.i = 0.;
#line 531 "TB05AD.f"
		i__7 = k + j * hinvb_dim1;
#line 531 "TB05AD.f"
		z__2.r = z__3.r * hinvb[i__7].r - z__3.i * hinvb[i__7].i, 
			z__2.i = z__3.r * hinvb[i__7].i + z__3.i * hinvb[i__7]
			.r;
#line 531 "TB05AD.f"
		z__1.r = g[i__5].r + z__2.r, z__1.i = g[i__5].i + z__2.i;
#line 531 "TB05AD.f"
		g[i__4].r = z__1.r, g[i__4].i = z__1.i;
#line 532 "TB05AD.f"
/* L130: */
#line 532 "TB05AD.f"
	    }

#line 534 "TB05AD.f"
/* L140: */
#line 534 "TB05AD.f"
	}

#line 536 "TB05AD.f"
/* L150: */
#line 536 "TB05AD.f"
    }

/*     G now contains the desired frequency response matrix. */
/*     Set the optimal workspace. */

#line 541 "TB05AD.f"
    dwork[1] = (doublereal) wrkopt;

#line 543 "TB05AD.f"
    return 0;
/* *** Last line of TB05AD *** */
} /* tb05ad_ */

