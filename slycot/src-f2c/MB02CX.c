#line 1 "MB02CX.f"
/* MB02CX.f -- translated by f2c (version 20100827).
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

#line 1 "MB02CX.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mb02cx_(char *typet, integer *p, integer *q, integer *k, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	cs, integer *lcs, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen typet_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static doublereal c__;
    static integer i__;
    static doublereal s, tau, beta;
    static integer ierr;
    extern /* Subroutine */ int ma02fd_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    static doublereal alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dlarf_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical isrow;
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), dgelqf_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dgeqrf_(integer *, integer *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *, integer *), xerbla_(char *
	    , integer *, ftnlen);
    static doublereal maxwrk;


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

/*     To bring the first blocks of a generator in proper form. */
/*     The columns / rows of the positive and negative generators */
/*     are contained in the arrays A and B, respectively. */
/*     Transformation information will be stored and can be applied */
/*     via SLICOT Library routine MB02CY. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYPET   CHARACTER*1 */
/*             Specifies the type of the generator, as follows: */
/*             = 'R':  A and B are the first blocks of the rows of the */
/*                     positive and negative generators; */
/*             = 'C':  A and B are the first blocks of the columns of the */
/*                     positive and negative generators. */
/*             Note:   in the sequel, the notation x / y means that */
/*                     x corresponds to TYPET = 'R' and y corresponds to */
/*                     TYPET = 'C'. */

/*     Input/Output Parameters */

/*     P       (input)  INTEGER */
/*             The number of rows / columns in A containing the positive */
/*             generators.  P >= 0. */

/*     Q       (input)  INTEGER */
/*             The number of rows / columns in B containing the negative */
/*             generators.  Q >= 0. */

/*     K       (input)  INTEGER */
/*             The number of columns / rows in A and B to be processed. */
/*             Normally, the size of the first block.  P >= K >= 0. */

/*     A       (input/output)  DOUBLE PRECISION array, dimension */
/*             (LDA, K) / (LDA, P) */
/*             On entry, the leading P-by-K upper / K-by-P lower */
/*             triangular part of this array must contain the rows / */
/*             columns of the positive part in the first block of the */
/*             generator. */
/*             On exit, the leading P-by-K upper / K-by-P lower */
/*             triangular part of this array contains the rows / columns */
/*             of the positive part in the first block of the proper */
/*             generator. */
/*             The lower / upper trapezoidal part is not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. */
/*             LDA >= MAX(1,P),    if TYPET = 'R'; */
/*             LDA >= MAX(1,K),    if TYPET = 'C'. */

/*     B       (input/output)  DOUBLE PRECISION array, dimension */
/*             (LDB, K) / (LDB, Q) */
/*             On entry, the leading Q-by-K / K-by-Q part of this array */
/*             must contain the rows / columns of the negative part in */
/*             the first block of the generator. */
/*             On exit, the leading Q-by-K / K-by-Q part of this array */
/*             contains part of the necessary information for the */
/*             Householder transformations. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= MAX(1,Q),    if TYPET = 'R'; */
/*             LDB >= MAX(1,K),    if TYPET = 'C'. */

/*     CS      (output)  DOUBLE PRECISION array, dimension (LCS) */
/*             On exit, the leading 2*K + MIN(K,Q) part of this array */
/*             contains necessary information for the SLICOT Library */
/*             routine MB02CY (modified hyperbolic rotation parameters */
/*             and scalar factors of the Householder transformations). */

/*     LCS     INTEGER */
/*             The length of the array CS.  LCS >= 2*K + MIN(K,Q). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -12,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,K). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  succesful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the reduction algorithm failed. The matrix */
/*                   associated with the generator is not (numerically) */
/*                   positive definite. */

/*     METHOD */

/*     If  TYPET = 'R',  a QR decomposition of B is first computed. */
/*     Then, the elements below the first row of each column i of B */
/*     are annihilated by a Householder transformation modifying the */
/*     first element in that column. This first element, in turn, is */
/*     then annihilated by a modified hyperbolic rotation, acting also */
/*     on the i-th row of A. */

/*     If  TYPET = 'C',  an LQ decomposition of B is first computed. */
/*     Then, the elements on the right of the first column of each row i */
/*     of B are annihilated by a Householder transformation modifying the */
/*     first element in that row. This first element, in turn, is */
/*     then annihilated by a modified hyperbolic rotation, acting also */
/*     on the i-th column of A. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Chemnitz, Germany, June 2000. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, July 2000, */
/*     February 2004. */

/*     KEYWORDS */

/*     Elementary matrix operations, Householder transformation, matrix */
/*     operations, Toeplitz matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

#line 177 "MB02CX.f"
    /* Parameter adjustments */
#line 177 "MB02CX.f"
    a_dim1 = *lda;
#line 177 "MB02CX.f"
    a_offset = 1 + a_dim1;
#line 177 "MB02CX.f"
    a -= a_offset;
#line 177 "MB02CX.f"
    b_dim1 = *ldb;
#line 177 "MB02CX.f"
    b_offset = 1 + b_dim1;
#line 177 "MB02CX.f"
    b -= b_offset;
#line 177 "MB02CX.f"
    --cs;
#line 177 "MB02CX.f"
    --dwork;
#line 177 "MB02CX.f"

#line 177 "MB02CX.f"
    /* Function Body */
#line 177 "MB02CX.f"
    *info = 0;
#line 178 "MB02CX.f"
    isrow = lsame_(typet, "R", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

#line 182 "MB02CX.f"
    if (! (isrow || lsame_(typet, "C", (ftnlen)1, (ftnlen)1))) {
#line 183 "MB02CX.f"
	*info = -1;
#line 184 "MB02CX.f"
    } else if (*p < 0) {
#line 185 "MB02CX.f"
	*info = -2;
#line 186 "MB02CX.f"
    } else if (*q < 0) {
#line 187 "MB02CX.f"
	*info = -3;
#line 188 "MB02CX.f"
    } else if (*k < 0 || *k > *p) {
#line 189 "MB02CX.f"
	*info = -4;
#line 190 "MB02CX.f"
    } else if (*lda < 1 || isrow && *lda < *p || ! isrow && *lda < *k) {
#line 192 "MB02CX.f"
	*info = -6;
#line 193 "MB02CX.f"
    } else if (*ldb < 1 || isrow && *ldb < *q || ! isrow && *ldb < *k) {
#line 195 "MB02CX.f"
	*info = -8;
#line 196 "MB02CX.f"
    } else if (*lcs < (*k << 1) + min(*k,*q)) {
#line 197 "MB02CX.f"
	*info = -10;
#line 198 "MB02CX.f"
    } else if (*ldwork < max(1,*k)) {
#line 199 "MB02CX.f"
	dwork[1] = (doublereal) max(1,*k);
#line 200 "MB02CX.f"
	*info = -12;
#line 201 "MB02CX.f"
    }

/*     Return if there were illegal values. */

#line 205 "MB02CX.f"
    if (*info != 0) {
#line 206 "MB02CX.f"
	i__1 = -(*info);
#line 206 "MB02CX.f"
	xerbla_("MB02CX", &i__1, (ftnlen)6);
#line 207 "MB02CX.f"
	return 0;
#line 208 "MB02CX.f"
    }

/*     Quick return if possible. */

#line 212 "MB02CX.f"
    if (min(*q,*k) == 0) {
#line 213 "MB02CX.f"
	dwork[1] = 1.;
#line 214 "MB02CX.f"
	return 0;
#line 215 "MB02CX.f"
    }

#line 217 "MB02CX.f"
    if (isrow) {

/*        The generator is row wise stored. */

/*        Step 0: Do QR decomposition of B. */

#line 223 "MB02CX.f"
	dgeqrf_(q, k, &b[b_offset], ldb, &cs[(*k << 1) + 1], &dwork[1], 
		ldwork, &ierr);
#line 224 "MB02CX.f"
	maxwrk = dwork[1];

#line 226 "MB02CX.f"
	i__1 = *k;
#line 226 "MB02CX.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Step 1: annihilate the i-th column of B. */

#line 230 "MB02CX.f"
	    if (*q > 1) {
#line 231 "MB02CX.f"
		i__2 = min(i__,*q);
#line 231 "MB02CX.f"
		dlarfg_(&i__2, &b[i__ * b_dim1 + 1], &b[i__ * b_dim1 + 2], &
			c__1, &tau);
#line 232 "MB02CX.f"
		alpha = b[i__ * b_dim1 + 1];
#line 233 "MB02CX.f"
		b[i__ * b_dim1 + 1] = 1.;
#line 234 "MB02CX.f"
		if (*k > i__) {
#line 234 "MB02CX.f"
		    i__2 = min(i__,*q);
#line 234 "MB02CX.f"
		    i__3 = *k - i__;
#line 234 "MB02CX.f"
		    dlarf_("Left", &i__2, &i__3, &b[i__ * b_dim1 + 1], &c__1, 
			    &tau, &b[(i__ + 1) * b_dim1 + 1], ldb, &dwork[1], 
			    (ftnlen)4);
#line 234 "MB02CX.f"
		}
#line 237 "MB02CX.f"
		b[i__ * b_dim1 + 1] = alpha;
#line 238 "MB02CX.f"
	    } else {
#line 239 "MB02CX.f"
		alpha = b[i__ * b_dim1 + 1];
#line 240 "MB02CX.f"
		tau = 0.;
#line 241 "MB02CX.f"
	    }

/*           Step 2: annihilate the top entry of the column. */

#line 245 "MB02CX.f"
	    beta = a[i__ + i__ * a_dim1];
#line 246 "MB02CX.f"
	    ma02fd_(&beta, &alpha, &c__, &s, &ierr);
#line 247 "MB02CX.f"
	    if (ierr != 0) {

/*              Error return:  The matrix is not positive definite. */

#line 251 "MB02CX.f"
		*info = 1;
#line 252 "MB02CX.f"
		return 0;
#line 253 "MB02CX.f"
	    }

#line 255 "MB02CX.f"
	    cs[(i__ << 1) - 1] = c__;
#line 256 "MB02CX.f"
	    cs[i__ * 2] = s;
#line 257 "MB02CX.f"
	    i__2 = *k - i__ + 1;
#line 257 "MB02CX.f"
	    d__1 = 1. / c__;
#line 257 "MB02CX.f"
	    dscal_(&i__2, &d__1, &a[i__ + i__ * a_dim1], lda);
#line 258 "MB02CX.f"
	    i__2 = *k - i__ + 1;
#line 258 "MB02CX.f"
	    d__1 = -s / c__;
#line 258 "MB02CX.f"
	    daxpy_(&i__2, &d__1, &b[i__ * b_dim1 + 1], ldb, &a[i__ + i__ * 
		    a_dim1], lda);
#line 259 "MB02CX.f"
	    i__2 = *k - i__ + 1;
#line 259 "MB02CX.f"
	    dscal_(&i__2, &c__, &b[i__ * b_dim1 + 1], ldb);
#line 260 "MB02CX.f"
	    i__2 = *k - i__ + 1;
#line 260 "MB02CX.f"
	    d__1 = -s;
#line 260 "MB02CX.f"
	    daxpy_(&i__2, &d__1, &a[i__ + i__ * a_dim1], lda, &b[i__ * b_dim1 
		    + 1], ldb);
#line 261 "MB02CX.f"
	    b[i__ * b_dim1 + 1] = tau;
#line 262 "MB02CX.f"
/* L10: */
#line 262 "MB02CX.f"
	}

#line 264 "MB02CX.f"
    } else {

/*        The generator is column wise stored. */

/*        Step 0: Do LQ decomposition of B. */

#line 270 "MB02CX.f"
	dgelqf_(k, q, &b[b_offset], ldb, &cs[(*k << 1) + 1], &dwork[1], 
		ldwork, &ierr);
#line 271 "MB02CX.f"
	maxwrk = dwork[1];

#line 273 "MB02CX.f"
	i__1 = *k;
#line 273 "MB02CX.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Step 1: annihilate the i-th row of B. */

#line 277 "MB02CX.f"
	    if (*q > 1) {
#line 278 "MB02CX.f"
		i__2 = min(i__,*q);
#line 278 "MB02CX.f"
		dlarfg_(&i__2, &b[i__ + b_dim1], &b[i__ + (b_dim1 << 1)], ldb,
			 &tau);
#line 279 "MB02CX.f"
		alpha = b[i__ + b_dim1];
#line 280 "MB02CX.f"
		b[i__ + b_dim1] = 1.;
#line 281 "MB02CX.f"
		if (*k > i__) {
#line 281 "MB02CX.f"
		    i__2 = *k - i__;
#line 281 "MB02CX.f"
		    i__3 = min(i__,*q);
#line 281 "MB02CX.f"
		    dlarf_("Right", &i__2, &i__3, &b[i__ + b_dim1], ldb, &tau,
			     &b[i__ + 1 + b_dim1], ldb, &dwork[1], (ftnlen)5);
#line 281 "MB02CX.f"
		}
#line 284 "MB02CX.f"
		b[i__ + b_dim1] = alpha;
#line 285 "MB02CX.f"
	    } else {
#line 286 "MB02CX.f"
		alpha = b[i__ + b_dim1];
#line 287 "MB02CX.f"
		tau = 0.;
#line 288 "MB02CX.f"
	    }

/*           Step 2: annihilate the left entry of the row. */

#line 292 "MB02CX.f"
	    beta = a[i__ + i__ * a_dim1];
#line 293 "MB02CX.f"
	    ma02fd_(&beta, &alpha, &c__, &s, &ierr);
#line 294 "MB02CX.f"
	    if (ierr != 0) {

/*              Error return:  The matrix is not positive definite. */

#line 298 "MB02CX.f"
		*info = 1;
#line 299 "MB02CX.f"
		return 0;
#line 300 "MB02CX.f"
	    }

#line 302 "MB02CX.f"
	    cs[(i__ << 1) - 1] = c__;
#line 303 "MB02CX.f"
	    cs[i__ * 2] = s;
#line 304 "MB02CX.f"
	    i__2 = *k - i__ + 1;
#line 304 "MB02CX.f"
	    d__1 = 1. / c__;
#line 304 "MB02CX.f"
	    dscal_(&i__2, &d__1, &a[i__ + i__ * a_dim1], &c__1);
#line 305 "MB02CX.f"
	    i__2 = *k - i__ + 1;
#line 305 "MB02CX.f"
	    d__1 = -s / c__;
#line 305 "MB02CX.f"
	    daxpy_(&i__2, &d__1, &b[i__ + b_dim1], &c__1, &a[i__ + i__ * 
		    a_dim1], &c__1);
#line 306 "MB02CX.f"
	    i__2 = *k - i__ + 1;
#line 306 "MB02CX.f"
	    dscal_(&i__2, &c__, &b[i__ + b_dim1], &c__1);
#line 307 "MB02CX.f"
	    i__2 = *k - i__ + 1;
#line 307 "MB02CX.f"
	    d__1 = -s;
#line 307 "MB02CX.f"
	    daxpy_(&i__2, &d__1, &a[i__ + i__ * a_dim1], &c__1, &b[i__ + 
		    b_dim1], &c__1);
#line 308 "MB02CX.f"
	    b[i__ + b_dim1] = tau;
#line 309 "MB02CX.f"
/* L20: */
#line 309 "MB02CX.f"
	}

#line 311 "MB02CX.f"
    }

#line 313 "MB02CX.f"
    dwork[1] = maxwrk;

#line 315 "MB02CX.f"
    return 0;

/* *** Last line of MB02CX *** */
} /* mb02cx_ */

