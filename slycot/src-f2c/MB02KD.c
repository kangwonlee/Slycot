#line 1 "MB02KD.f"
/* MB02KD.f -- translated by f2c (version 20100827).
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

#line 1 "MB02KD.f"
/* Table of constant values */

static doublereal c_b9 = 0.;
static integer c__1 = 1;
static doublereal c_b23 = 1.;

/* Subroutine */ int mb02kd_(char *ldblk, char *trans, integer *k, integer *l,
	 integer *m, integer *n, integer *r__, doublereal *alpha, doublereal *
	beta, doublereal *tc, integer *ldtc, doublereal *tr, integer *ldtr, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen ldblk_len, 
	ftnlen trans_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, c_dim1, c_offset, tc_dim1, tc_offset, tr_dim1, 
	    tr_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Builtin functions */
    double atan(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j, p, p1, p2, q1, q2, r1, r2, s1, s2;
    static doublereal t1, t2, cf;
    static integer pb, pc, jj, kk, ll, mk, ln, ir, nl;
    static doublereal sf, th;
    static integer pp, pt, icp, icq, len, pdw, dimb, dimc;
    static doublereal coef, scal;
    static integer meth, ierr, shft;
    static char wght[1];
    static integer wpos;
    extern /* Subroutine */ int dg01od_(char *, char *, integer *, doublereal 
	    *, doublereal *, integer *, ftnlen, ftnlen), dscal_(integer *, 
	    doublereal *, doublereal *, integer *), dgemm_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen);
    static doublereal param;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical fullc;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical ltran;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical lmult;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
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

/*     To compute the matrix product */

/*               C = alpha*op( T )*B + beta*C, */

/*     where alpha and beta are scalars and T is a block Toeplitz matrix */
/*     specified by its first block column TC and first block row TR; */
/*     B and C are general matrices of appropriate dimensions. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     LDBLK   CHARACTER*1 */
/*             Specifies where the (1,1)-block of T is stored, as */
/*             follows: */
/*             = 'C':  in the first block of TC; */
/*             = 'R':  in the first block of TR. */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of op( T ) to be used in the matrix */
/*             multiplication as follows: */
/*             = 'N':  op( T ) = T; */
/*             = 'T':  op( T ) = T'; */
/*             = 'C':  op( T ) = T'. */

/*     Input/Output Parameters */

/*     K       (input) INTEGER */
/*             The number of rows in the blocks of T.  K >= 0. */

/*     L       (input) INTEGER */
/*             The number of columns in the blocks of T.  L >= 0. */

/*     M       (input) INTEGER */
/*             The number of blocks in the first block column of T. */
/*             M >= 0. */

/*     N       (input) INTEGER */
/*             The number of blocks in the first block row of T.  N >= 0. */

/*     R       (input) INTEGER */
/*             The number of columns in B and C.  R >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar alpha. When alpha is zero then TC, TR and B */
/*             are not referenced. */

/*     BETA    (input) DOUBLE PRECISION */
/*             The scalar beta. When beta is zero then C need not be set */
/*             before entry. */

/*     TC      (input)  DOUBLE PRECISION array, dimension (LDTC,L) */
/*             On entry with LDBLK = 'C', the leading M*K-by-L part of */
/*             this array must contain the first block column of T. */
/*             On entry with LDBLK = 'R', the leading (M-1)*K-by-L part */
/*             of this array must contain the 2nd to the M-th blocks of */
/*             the first block column of T. */

/*     LDTC    INTEGER */
/*             The leading dimension of the array TC. */
/*             LDTC >= MAX(1,M*K),      if LDBLK = 'C'; */
/*             LDTC >= MAX(1,(M-1)*K),  if LDBLK = 'R'. */

/*     TR      (input)  DOUBLE PRECISION array, dimension (LDTR,k) */
/*             where k is (N-1)*L when LDBLK = 'C' and is N*L when */
/*             LDBLK = 'R'. */
/*             On entry with LDBLK = 'C', the leading K-by-(N-1)*L part */
/*             of this array must contain the 2nd to the N-th blocks of */
/*             the first block row of T. */
/*             On entry with LDBLK = 'R', the leading K-by-N*L part of */
/*             this array must contain the first block row of T. */

/*     LDTR    INTEGER */
/*             The leading dimension of the array TR.  LDTR >= MAX(1,K). */

/*     B       (input)  DOUBLE PRECISION array, dimension (LDB,R) */
/*             On entry with TRANS = 'N', the leading N*L-by-R part of */
/*             this array must contain the matrix B. */
/*             On entry with TRANS = 'T' or TRANS = 'C', the leading */
/*             M*K-by-R part of this array must contain the matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= MAX(1,N*L),  if TRANS = 'N'; */
/*             LDB >= MAX(1,M*K),  if TRANS = 'T' or TRANS = 'C'. */

/*     C       (input/output)  DOUBLE PRECISION array, dimension (LDC,R) */
/*             On entry with TRANS = 'N', the leading M*K-by-R part of */
/*             this array must contain the matrix C. */
/*             On entry with TRANS = 'T' or TRANS = 'C', the leading */
/*             N*L-by-R part of this array must contain the matrix C. */
/*             On exit with TRANS = 'N', the leading M*K-by-R part of */
/*             this array contains the updated matrix C. */
/*             On exit with TRANS = 'T' or TRANS = 'C', the leading */
/*             N*L-by-R part of this array contains the updated matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C. */
/*             LDC >= MAX(1,M*K),  if TRANS = 'N'; */
/*             LDC >= MAX(1,N*L),  if TRANS = 'T' or TRANS = 'C'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -19,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= 1. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     For point Toeplitz matrices or sufficiently large block Toeplitz */
/*     matrices, this algorithm uses convolution algorithms based on */
/*     the fast Hartley transforms [1]. Otherwise, TC is copied in */
/*     reversed order into the workspace such that C can be computed from */
/*     barely M matrix-by-matrix multiplications. */

/*     REFERENCES */

/*     [1] Van Loan, Charles. */
/*         Computational frameworks for the fast Fourier transform. */
/*         SIAM, 1992. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires O( (K*L+R*L+K*R)*(N+M)*log(N+M) + K*L*R ) */
/*     floating point operations. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, May 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, June 2001, */
/*     March 2004. */

/*     KEYWORDS */

/*     Convolution, elementary matrix operations, */
/*     fast Hartley transform, Toeplitz matrix. */

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

#line 213 "MB02KD.f"
    /* Parameter adjustments */
#line 213 "MB02KD.f"
    tc_dim1 = *ldtc;
#line 213 "MB02KD.f"
    tc_offset = 1 + tc_dim1;
#line 213 "MB02KD.f"
    tc -= tc_offset;
#line 213 "MB02KD.f"
    tr_dim1 = *ldtr;
#line 213 "MB02KD.f"
    tr_offset = 1 + tr_dim1;
#line 213 "MB02KD.f"
    tr -= tr_offset;
#line 213 "MB02KD.f"
    b_dim1 = *ldb;
#line 213 "MB02KD.f"
    b_offset = 1 + b_dim1;
#line 213 "MB02KD.f"
    b -= b_offset;
#line 213 "MB02KD.f"
    c_dim1 = *ldc;
#line 213 "MB02KD.f"
    c_offset = 1 + c_dim1;
#line 213 "MB02KD.f"
    c__ -= c_offset;
#line 213 "MB02KD.f"
    --dwork;
#line 213 "MB02KD.f"

#line 213 "MB02KD.f"
    /* Function Body */
#line 213 "MB02KD.f"
    *info = 0;
#line 214 "MB02KD.f"
    fullc = lsame_(ldblk, "C", (ftnlen)1, (ftnlen)1);
#line 215 "MB02KD.f"
    ltran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1);
#line 216 "MB02KD.f"
    lmult = *alpha != 0.;
#line 217 "MB02KD.f"
    mk = *m * *k;
#line 218 "MB02KD.f"
    nl = *n * *l;

/*     Check the scalar input parameters. */

#line 222 "MB02KD.f"
    if (! (fullc || lsame_(ldblk, "R", (ftnlen)1, (ftnlen)1))) {
#line 223 "MB02KD.f"
	*info = -1;
#line 224 "MB02KD.f"
    } else if (! (ltran || lsame_(trans, "N", (ftnlen)1, (ftnlen)1))) {
#line 225 "MB02KD.f"
	*info = -2;
#line 226 "MB02KD.f"
    } else if (*k < 0) {
#line 227 "MB02KD.f"
	*info = -3;
#line 228 "MB02KD.f"
    } else if (*l < 0) {
#line 229 "MB02KD.f"
	*info = -4;
#line 230 "MB02KD.f"
    } else if (*m < 0) {
#line 231 "MB02KD.f"
	*info = -5;
#line 232 "MB02KD.f"
    } else if (*n < 0) {
#line 233 "MB02KD.f"
	*info = -6;
#line 234 "MB02KD.f"
    } else if (*r__ < 0) {
#line 235 "MB02KD.f"
	*info = -7;
#line 236 "MB02KD.f"
    } else if (lmult && fullc && *ldtc < max(1,mk)) {
#line 237 "MB02KD.f"
	*info = -11;
#line 238 "MB02KD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 238 "MB02KD.f"
	i__1 = 1, i__2 = (*m - 1) * *k;
#line 238 "MB02KD.f"
	if (lmult && ! fullc && *ldtc < max(i__1,i__2)) {
#line 240 "MB02KD.f"
	    *info = -11;
#line 241 "MB02KD.f"
	} else if (lmult && *ldtr < max(1,*k)) {
#line 242 "MB02KD.f"
	    *info = -13;
#line 243 "MB02KD.f"
	} else if (lmult && ! ltran && *ldb < max(1,nl)) {
#line 244 "MB02KD.f"
	    *info = -15;
#line 245 "MB02KD.f"
	} else if (lmult && ltran && *ldb < max(1,mk)) {
#line 246 "MB02KD.f"
	    *info = -15;
#line 247 "MB02KD.f"
	} else if (! ltran && *ldc < max(1,mk)) {
#line 248 "MB02KD.f"
	    *info = -17;
#line 249 "MB02KD.f"
	} else if (ltran && *ldc < max(1,nl)) {
#line 250 "MB02KD.f"
	    *info = -17;
#line 251 "MB02KD.f"
	} else if (*ldwork < 1) {
#line 252 "MB02KD.f"
	    dwork[1] = 1.;
#line 253 "MB02KD.f"
	    *info = -19;
#line 254 "MB02KD.f"
	}
#line 254 "MB02KD.f"
    }

/*     Return if there were illegal values. */

#line 258 "MB02KD.f"
    if (*info != 0) {
#line 259 "MB02KD.f"
	i__1 = -(*info);
#line 259 "MB02KD.f"
	xerbla_("MB02KD", &i__1, (ftnlen)6);
#line 260 "MB02KD.f"
	return 0;
#line 261 "MB02KD.f"
    }

/*     Scale C beforehand. */

#line 265 "MB02KD.f"
    if (*beta == 0.) {
#line 266 "MB02KD.f"
	if (ltran) {
#line 267 "MB02KD.f"
	    dlaset_("All", &nl, r__, &c_b9, &c_b9, &c__[c_offset], ldc, (
		    ftnlen)3);
#line 268 "MB02KD.f"
	} else {
#line 269 "MB02KD.f"
	    dlaset_("All", &mk, r__, &c_b9, &c_b9, &c__[c_offset], ldc, (
		    ftnlen)3);
#line 270 "MB02KD.f"
	}
#line 271 "MB02KD.f"
    } else if (*beta != 1.) {
#line 272 "MB02KD.f"
	if (ltran) {

#line 274 "MB02KD.f"
	    i__1 = *r__;
#line 274 "MB02KD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 275 "MB02KD.f"
		dscal_(&nl, beta, &c__[i__ * c_dim1 + 1], &c__1);
#line 276 "MB02KD.f"
/* L10: */
#line 276 "MB02KD.f"
	    }

#line 278 "MB02KD.f"
	} else {

#line 280 "MB02KD.f"
	    i__1 = *r__;
#line 280 "MB02KD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 281 "MB02KD.f"
		dscal_(&mk, beta, &c__[i__ * c_dim1 + 1], &c__1);
#line 282 "MB02KD.f"
/* L20: */
#line 282 "MB02KD.f"
	    }

#line 284 "MB02KD.f"
	}
#line 285 "MB02KD.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 289 "MB02KD.f"
    i__1 = min(mk,nl);
#line 289 "MB02KD.f"
    if (! lmult || min(i__1,*r__) == 0) {
#line 290 "MB02KD.f"
	dwork[1] = 1.;
#line 291 "MB02KD.f"
	return 0;
#line 292 "MB02KD.f"
    }

/*     The parameter PARAM is the watershed between conventional */
/*     multiplication and convolution. This is of course depending */
/*     on the used computer architecture. The lower this value is set */
/*     the more likely the routine will use convolution to compute */
/*     op( T )*B. Note that if there is enough workspace available, */
/*     convolution is always used for point Toeplitz matrices. */

#line 301 "MB02KD.f"
    param = 950.;

/*     Decide which method to choose, based on the block sizes and */
/*     the available workspace. */

#line 306 "MB02KD.f"
    len = 1;
#line 307 "MB02KD.f"
    p = 0;

#line 309 "MB02KD.f"
L30:
#line 310 "MB02KD.f"
    if (len < *m + *n - 1) {
#line 311 "MB02KD.f"
	len <<= 1;
#line 312 "MB02KD.f"
	++p;
#line 313 "MB02KD.f"
	goto L30;
#line 314 "MB02KD.f"
    }

#line 316 "MB02KD.f"
    coef = (doublereal) (*m * *n) * 3. * (doublereal) (*k * *l) * (doublereal)
	     (*r__) / (doublereal) (len * (*k * *l + *l * *r__ + *k * *r__));

#line 319 "MB02KD.f"
    if (fullc) {
#line 320 "MB02KD.f"
	p1 = mk * *l;
#line 321 "MB02KD.f"
	shft = 0;
#line 322 "MB02KD.f"
    } else {
#line 323 "MB02KD.f"
	p1 = (*m - 1) * *k * *l;
#line 324 "MB02KD.f"
	shft = 1;
#line 325 "MB02KD.f"
    }
#line 326 "MB02KD.f"
    if (*k * *l == 1 && min(*m,*n) > 1) {
#line 327 "MB02KD.f"
	wrkopt = len * (*r__ + 2) - p;
#line 328 "MB02KD.f"
	meth = 3;
#line 329 "MB02KD.f"
    } else if (len < *m * *n && coef >= param) {
#line 330 "MB02KD.f"
	wrkopt = len * (*k * *l + *k * *r__ + *l * *r__ + 1) - p;
#line 331 "MB02KD.f"
	meth = 3;
#line 332 "MB02KD.f"
    } else {
#line 333 "MB02KD.f"
	meth = 2;
#line 334 "MB02KD.f"
	wrkopt = p1;
#line 335 "MB02KD.f"
    }

#line 337 "MB02KD.f"
    if (*ldwork < wrkopt) {
#line 337 "MB02KD.f"
	--meth;
#line 337 "MB02KD.f"
    }
#line 338 "MB02KD.f"
    if (*ldwork < p1) {
#line 338 "MB02KD.f"
	meth = 1;
#line 338 "MB02KD.f"
    }

/*     Start computations. */

#line 342 "MB02KD.f"
    if (meth == 1 && ! ltran) {

/*        Method 1 is the most unlucky way to multiply Toeplitz matrices */
/*        with vectors. Due to the memory restrictions it is not */
/*        possible to flip TC. */

#line 348 "MB02KD.f"
	pc = 1;

#line 350 "MB02KD.f"
	i__1 = *m;
#line 350 "MB02KD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 351 "MB02KD.f"
	    pt = (i__ - 1 - shft) * *k + 1;
#line 352 "MB02KD.f"
	    pb = 1;

#line 354 "MB02KD.f"
	    i__2 = i__;
#line 354 "MB02KD.f"
	    for (j = shft + 1; j <= i__2; ++j) {
#line 355 "MB02KD.f"
		dgemm_("No Transpose", "No Transpose", k, r__, l, alpha, &tc[
			pt + tc_dim1], ldtc, &b[pb + b_dim1], ldb, &c_b23, &
			c__[pc + c_dim1], ldc, (ftnlen)12, (ftnlen)12);
#line 358 "MB02KD.f"
		pt -= *k;
#line 359 "MB02KD.f"
		pb += *l;
#line 360 "MB02KD.f"
/* L40: */
#line 360 "MB02KD.f"
	    }

#line 362 "MB02KD.f"
	    if (*n > i__ - shft) {
#line 363 "MB02KD.f"
		i__2 = (*n - i__ + shft) * *l;
#line 363 "MB02KD.f"
		dgemm_("No Transpose", "No Transpose", k, r__, &i__2, alpha, &
			tr[tr_offset], ldtr, &b[pb + b_dim1], ldb, &c_b23, &
			c__[pc + c_dim1], ldc, (ftnlen)12, (ftnlen)12);
#line 366 "MB02KD.f"
	    }
#line 367 "MB02KD.f"
	    pc += *k;
#line 368 "MB02KD.f"
/* L50: */
#line 368 "MB02KD.f"
	}

#line 370 "MB02KD.f"
    } else if (meth == 1 && ltran) {

#line 372 "MB02KD.f"
	pb = 1;

#line 374 "MB02KD.f"
	i__1 = *m;
#line 374 "MB02KD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 375 "MB02KD.f"
	    pt = (i__ - 1 - shft) * *k + 1;
#line 376 "MB02KD.f"
	    pc = 1;

#line 378 "MB02KD.f"
	    i__2 = i__;
#line 378 "MB02KD.f"
	    for (j = shft + 1; j <= i__2; ++j) {
#line 379 "MB02KD.f"
		dgemm_("Transpose", "No Transpose", l, r__, k, alpha, &tc[pt 
			+ tc_dim1], ldtc, &b[pb + b_dim1], ldb, &c_b23, &c__[
			pc + c_dim1], ldc, (ftnlen)9, (ftnlen)12);
#line 382 "MB02KD.f"
		pt -= *k;
#line 383 "MB02KD.f"
		pc += *l;
#line 384 "MB02KD.f"
/* L60: */
#line 384 "MB02KD.f"
	    }

#line 386 "MB02KD.f"
	    if (*n > i__ - shft) {
#line 387 "MB02KD.f"
		i__2 = (*n - i__ + shft) * *l;
#line 387 "MB02KD.f"
		dgemm_("Transpose", "No Transpose", &i__2, r__, k, alpha, &tr[
			tr_offset], ldtr, &b[pb + b_dim1], ldb, &c_b23, &c__[
			pc + c_dim1], ldc, (ftnlen)9, (ftnlen)12);
#line 390 "MB02KD.f"
	    }
#line 391 "MB02KD.f"
	    pb += *k;
#line 392 "MB02KD.f"
/* L70: */
#line 392 "MB02KD.f"
	}

#line 394 "MB02KD.f"
    } else if (meth == 2 && ! ltran) {

/*        In method 2 TC is flipped resulting in less calls to the BLAS */
/*        routine DGEMM. Actually this seems often to be the best way to */
/*        multiply with Toeplitz matrices except the point Toeplitz */
/*        case. */

#line 401 "MB02KD.f"
	pt = (*m - 1 - shft) * *k + 1;

#line 403 "MB02KD.f"
	i__1 = (*m - shft) * *k * *l;
#line 403 "MB02KD.f"
	i__2 = *k * *l;
#line 403 "MB02KD.f"
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 404 "MB02KD.f"
	    dlacpy_("All", k, l, &tc[pt + tc_dim1], ldtc, &dwork[i__], k, (
		    ftnlen)3);
#line 405 "MB02KD.f"
	    pt -= *k;
#line 406 "MB02KD.f"
/* L80: */
#line 406 "MB02KD.f"
	}

#line 408 "MB02KD.f"
	pt = (*m - 1) * *k * *l + 1;
#line 409 "MB02KD.f"
	pc = 1;

#line 411 "MB02KD.f"
	i__2 = *m;
#line 411 "MB02KD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MIN */
#line 412 "MB02KD.f"
	    i__3 = i__ - shft;
#line 412 "MB02KD.f"
	    i__1 = min(i__3,*n) * *l;
#line 412 "MB02KD.f"
	    dgemm_("No Transpose", "No Transpose", k, r__, &i__1, alpha, &
		    dwork[pt], k, &b[b_offset], ldb, &c_b23, &c__[pc + c_dim1]
		    , ldc, (ftnlen)12, (ftnlen)12);
#line 415 "MB02KD.f"
	    if (*n > i__ - shft) {
#line 416 "MB02KD.f"
		i__1 = (*n - i__ + shft) * *l;
#line 416 "MB02KD.f"
		dgemm_("No Transpose", "No Transpose", k, r__, &i__1, alpha, &
			tr[tr_offset], ldtr, &b[(i__ - shft) * *l + 1 + 
			b_dim1], ldb, &c_b23, &c__[pc + c_dim1], ldc, (ftnlen)
			12, (ftnlen)12);
#line 419 "MB02KD.f"
	    }
#line 420 "MB02KD.f"
	    pc += *k;
#line 421 "MB02KD.f"
	    pt -= *k * *l;
#line 422 "MB02KD.f"
/* L90: */
#line 422 "MB02KD.f"
	}

#line 424 "MB02KD.f"
    } else if (meth == 2 && ltran) {

#line 426 "MB02KD.f"
	pt = (*m - 1 - shft) * *k + 1;

#line 428 "MB02KD.f"
	i__2 = (*m - shft) * *k * *l;
#line 428 "MB02KD.f"
	i__1 = *k * *l;
#line 428 "MB02KD.f"
	for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
#line 429 "MB02KD.f"
	    dlacpy_("All", k, l, &tc[pt + tc_dim1], ldtc, &dwork[i__], k, (
		    ftnlen)3);
#line 430 "MB02KD.f"
	    pt -= *k;
#line 431 "MB02KD.f"
/* L100: */
#line 431 "MB02KD.f"
	}

#line 433 "MB02KD.f"
	pt = (*m - 1) * *k * *l + 1;
#line 434 "MB02KD.f"
	pb = 1;

#line 436 "MB02KD.f"
	i__1 = *m;
#line 436 "MB02KD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
#line 437 "MB02KD.f"
	    i__3 = i__ - shft;
#line 437 "MB02KD.f"
	    i__2 = min(i__3,*n) * *l;
#line 437 "MB02KD.f"
	    dgemm_("Tranpose", "No Transpose", &i__2, r__, k, alpha, &dwork[
		    pt], k, &b[pb + b_dim1], ldb, &c_b23, &c__[c_offset], ldc,
		     (ftnlen)8, (ftnlen)12);
#line 440 "MB02KD.f"
	    if (*n > i__ - shft) {
#line 441 "MB02KD.f"
		i__2 = (*n - i__ + shft) * *l;
#line 441 "MB02KD.f"
		dgemm_("Transpose", "No Transpose", &i__2, r__, k, alpha, &tr[
			tr_offset], ldtr, &b[pb + b_dim1], ldb, &c_b23, &c__[(
			i__ - shft) * *l + 1 + c_dim1], ldc, (ftnlen)9, (
			ftnlen)12);
#line 444 "MB02KD.f"
	    }
#line 445 "MB02KD.f"
	    pb += *k;
#line 446 "MB02KD.f"
	    pt -= *k * *l;
#line 447 "MB02KD.f"
/* L110: */
#line 447 "MB02KD.f"
	}

#line 449 "MB02KD.f"
    } else if (meth == 3) {

/*        In method 3 the matrix-vector product is computed by a suitable */
/*        block convolution via fast Hartley transforms similar to the */
/*        SLICOT routine DE01PD. */

/*        Step 1: Copy input data into the workspace arrays. */

#line 457 "MB02KD.f"
	pdw = 1;
#line 458 "MB02KD.f"
	if (ltran) {
#line 459 "MB02KD.f"
	    dimb = *k;
#line 460 "MB02KD.f"
	    dimc = *l;
#line 461 "MB02KD.f"
	} else {
#line 462 "MB02KD.f"
	    dimb = *l;
#line 463 "MB02KD.f"
	    dimc = *k;
#line 464 "MB02KD.f"
	}
#line 465 "MB02KD.f"
	pb = len * *k * *l;
#line 466 "MB02KD.f"
	pc = len * (*k * *l + dimb * *r__);
#line 467 "MB02KD.f"
	if (ltran) {
#line 468 "MB02KD.f"
	    if (fullc) {
#line 469 "MB02KD.f"
		i__1 = len * *k;
#line 469 "MB02KD.f"
		dlacpy_("All", k, l, &tc[tc_offset], ldtc, &dwork[1], &i__1, (
			ftnlen)3);
#line 470 "MB02KD.f"
	    }

#line 472 "MB02KD.f"
	    i__1 = *n - 1 + shft;
#line 472 "MB02KD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 473 "MB02KD.f"
		i__2 = len * *k;
#line 473 "MB02KD.f"
		dlacpy_("All", k, l, &tr[((i__ - 1) * *l + 1) * tr_dim1 + 1], 
			ldtr, &dwork[(i__ - shft) * *k + 1], &i__2, (ftnlen)3)
			;
#line 475 "MB02KD.f"
/* L120: */
#line 475 "MB02KD.f"
	    }

#line 477 "MB02KD.f"
	    pdw = *n * *k + 1;
#line 478 "MB02KD.f"
	    r1 = (len - *m - *n + 1) * *k;
#line 479 "MB02KD.f"
	    i__1 = len * *k;
#line 479 "MB02KD.f"
	    dlaset_("All", &r1, l, &c_b9, &c_b9, &dwork[pdw], &i__1, (ftnlen)
		    3);
#line 480 "MB02KD.f"
	    pdw += r1;

#line 482 "MB02KD.f"
	    i__1 = *k - shft * *k + 1;
#line 482 "MB02KD.f"
	    i__2 = -(*k);
#line 482 "MB02KD.f"
	    for (i__ = (*m - 1 - shft) * *k + 1; i__2 < 0 ? i__ >= i__1 : i__ 
		    <= i__1; i__ += i__2) {
#line 483 "MB02KD.f"
		i__3 = len * *k;
#line 483 "MB02KD.f"
		dlacpy_("All", k, l, &tc[i__ + tc_dim1], ldtc, &dwork[pdw], &
			i__3, (ftnlen)3);
#line 485 "MB02KD.f"
		pdw += *k;
#line 486 "MB02KD.f"
/* L130: */
#line 486 "MB02KD.f"
	    }

#line 488 "MB02KD.f"
	    pdw = pb + 1;
#line 489 "MB02KD.f"
	    i__2 = len * *k;
#line 489 "MB02KD.f"
	    dlacpy_("All", &mk, r__, &b[b_offset], ldb, &dwork[pdw], &i__2, (
		    ftnlen)3);
#line 490 "MB02KD.f"
	    pdw += mk;
#line 491 "MB02KD.f"
	    i__2 = (len - *m) * *k;
#line 491 "MB02KD.f"
	    i__1 = len * *k;
#line 491 "MB02KD.f"
	    dlaset_("All", &i__2, r__, &c_b9, &c_b9, &dwork[pdw], &i__1, (
		    ftnlen)3);
#line 493 "MB02KD.f"
	} else {
#line 494 "MB02KD.f"
	    if (! fullc) {
#line 495 "MB02KD.f"
		i__2 = len * *k;
#line 495 "MB02KD.f"
		dlacpy_("All", k, l, &tr[tr_offset], ldtr, &dwork[1], &i__2, (
			ftnlen)3);
#line 496 "MB02KD.f"
	    }
#line 497 "MB02KD.f"
	    i__2 = (*m - shft) * *k;
#line 497 "MB02KD.f"
	    i__1 = len * *k;
#line 497 "MB02KD.f"
	    dlacpy_("All", &i__2, l, &tc[tc_offset], ldtc, &dwork[shft * *k + 
		    1], &i__1, (ftnlen)3);
#line 499 "MB02KD.f"
	    pdw = mk + 1;
#line 500 "MB02KD.f"
	    r1 = (len - *m - *n + 1) * *k;
#line 501 "MB02KD.f"
	    i__2 = len * *k;
#line 501 "MB02KD.f"
	    dlaset_("All", &r1, l, &c_b9, &c_b9, &dwork[pdw], &i__2, (ftnlen)
		    3);
#line 502 "MB02KD.f"
	    pdw += r1;

#line 504 "MB02KD.f"
	    i__2 = shft * *l + 1;
#line 504 "MB02KD.f"
	    i__1 = -(*l);
#line 504 "MB02KD.f"
	    for (i__ = (*n - 2 + shft) * *l + 1; i__1 < 0 ? i__ >= i__2 : i__ 
		    <= i__2; i__ += i__1) {
#line 505 "MB02KD.f"
		i__3 = len * *k;
#line 505 "MB02KD.f"
		dlacpy_("All", k, l, &tr[i__ * tr_dim1 + 1], ldtr, &dwork[pdw]
			, &i__3, (ftnlen)3);
#line 507 "MB02KD.f"
		pdw += *k;
#line 508 "MB02KD.f"
/* L140: */
#line 508 "MB02KD.f"
	    }

#line 510 "MB02KD.f"
	    pdw = pb + 1;
#line 511 "MB02KD.f"
	    i__1 = len * *l;
#line 511 "MB02KD.f"
	    dlacpy_("All", &nl, r__, &b[b_offset], ldb, &dwork[pdw], &i__1, (
		    ftnlen)3);
#line 512 "MB02KD.f"
	    pdw += nl;
#line 513 "MB02KD.f"
	    i__1 = (len - *n) * *l;
#line 513 "MB02KD.f"
	    i__2 = len * *l;
#line 513 "MB02KD.f"
	    dlaset_("All", &i__1, r__, &c_b9, &c_b9, &dwork[pdw], &i__2, (
		    ftnlen)3);
#line 515 "MB02KD.f"
	}

/*        Take point Toeplitz matrices into extra consideration. */

#line 519 "MB02KD.f"
	if (*k * *l == 1) {
#line 520 "MB02KD.f"
	    *(unsigned char *)wght = 'N';
#line 521 "MB02KD.f"
	    dg01od_("OutputScrambled", wght, &len, &dwork[1], &dwork[pc + 1], 
		    &ierr, (ftnlen)15, (ftnlen)1);

#line 524 "MB02KD.f"
	    i__1 = pb + len * *r__ - 1;
#line 524 "MB02KD.f"
	    i__2 = len;
#line 524 "MB02KD.f"
	    for (i__ = pb; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) 
		    {
#line 525 "MB02KD.f"
		dg01od_("OutputScrambled", wght, &len, &dwork[i__ + 1], &
			dwork[pc + 1], &ierr, (ftnlen)15, (ftnlen)1);
#line 527 "MB02KD.f"
		scal = *alpha / (doublereal) len;
#line 528 "MB02KD.f"
		dwork[i__ + 1] = scal * dwork[i__ + 1] * dwork[1];
#line 529 "MB02KD.f"
		dwork[i__ + 2] = scal * dwork[i__ + 2] * dwork[2];
#line 530 "MB02KD.f"
		scal /= 2.;

#line 532 "MB02KD.f"
		ln = 1;

#line 534 "MB02KD.f"
		i__3 = p - 1;
#line 534 "MB02KD.f"
		for (ll = 1; ll <= i__3; ++ll) {
#line 535 "MB02KD.f"
		    ln <<= 1;
#line 536 "MB02KD.f"
		    r1 = ln << 1;

#line 538 "MB02KD.f"
		    i__4 = ln + ln / 2;
#line 538 "MB02KD.f"
		    for (p1 = ln + 1; p1 <= i__4; ++p1) {
#line 539 "MB02KD.f"
			t1 = dwork[p1] + dwork[r1];
#line 540 "MB02KD.f"
			t2 = dwork[p1] - dwork[r1];
#line 541 "MB02KD.f"
			th = t2 * dwork[i__ + p1];
#line 542 "MB02KD.f"
			dwork[i__ + p1] = scal * (t1 * dwork[i__ + p1] + t2 * 
				dwork[i__ + r1]);
#line 544 "MB02KD.f"
			dwork[i__ + r1] = scal * (t1 * dwork[i__ + r1] - th);
#line 545 "MB02KD.f"
			--r1;
#line 546 "MB02KD.f"
/* L150: */
#line 546 "MB02KD.f"
		    }

#line 548 "MB02KD.f"
/* L160: */
#line 548 "MB02KD.f"
		}

#line 550 "MB02KD.f"
		dg01od_("InputScrambled", wght, &len, &dwork[i__ + 1], &dwork[
			pc + 1], &ierr, (ftnlen)14, (ftnlen)1);
#line 552 "MB02KD.f"
/* L170: */
#line 552 "MB02KD.f"
	    }

#line 554 "MB02KD.f"
	    pc = pb;
#line 555 "MB02KD.f"
	    goto L420;
#line 556 "MB02KD.f"
	}

/*        Step 2: Compute the weights for the Hartley transforms. */

#line 560 "MB02KD.f"
	pdw = pc;
#line 561 "MB02KD.f"
	r1 = 1;
#line 562 "MB02KD.f"
	ln = 1;
#line 563 "MB02KD.f"
	th = atan(1.) * 4. / (doublereal) len;

#line 565 "MB02KD.f"
	i__2 = p - 2;
#line 565 "MB02KD.f"
	for (ll = 1; ll <= i__2; ++ll) {
#line 566 "MB02KD.f"
	    ln <<= 1;
#line 567 "MB02KD.f"
	    th *= 2.;
#line 568 "MB02KD.f"
	    cf = cos(th);
#line 569 "MB02KD.f"
	    sf = sin(th);
#line 570 "MB02KD.f"
	    dwork[pdw + r1] = cf;
#line 571 "MB02KD.f"
	    dwork[pdw + r1 + 1] = sf;
#line 572 "MB02KD.f"
	    r1 += 2;

#line 574 "MB02KD.f"
	    i__1 = ln - 2;
#line 574 "MB02KD.f"
	    for (i__ = 1; i__ <= i__1; i__ += 2) {
#line 575 "MB02KD.f"
		dwork[pdw + r1] = cf * dwork[pdw + i__] - sf * dwork[pdw + 
			i__ + 1];
#line 576 "MB02KD.f"
		dwork[pdw + r1 + 1] = sf * dwork[pdw + i__] + cf * dwork[pdw 
			+ i__ + 1];
#line 577 "MB02KD.f"
		r1 += 2;
#line 578 "MB02KD.f"
/* L180: */
#line 578 "MB02KD.f"
	    }

#line 580 "MB02KD.f"
/* L190: */
#line 580 "MB02KD.f"
	}

#line 582 "MB02KD.f"
	p1 = 3;
#line 583 "MB02KD.f"
	q1 = r1 - 2;

#line 585 "MB02KD.f"
	for (ll = p - 2; ll >= 1; --ll) {

#line 587 "MB02KD.f"
	    i__2 = q1;
#line 587 "MB02KD.f"
	    for (i__ = p1; i__ <= i__2; i__ += 4) {
#line 588 "MB02KD.f"
		dwork[pdw + r1] = dwork[pdw + i__];
#line 589 "MB02KD.f"
		dwork[pdw + r1 + 1] = dwork[pdw + i__ + 1];
#line 590 "MB02KD.f"
		r1 += 2;
#line 591 "MB02KD.f"
/* L200: */
#line 591 "MB02KD.f"
	    }

#line 593 "MB02KD.f"
	    p1 = q1 + 4;
#line 594 "MB02KD.f"
	    q1 = r1 - 2;
#line 595 "MB02KD.f"
/* L210: */
#line 595 "MB02KD.f"
	}

/*        Step 3: Compute the Hartley transforms with scrambled output. */

#line 599 "MB02KD.f"
	j = 0;
#line 600 "MB02KD.f"
	kk = *k;

/*        WHILE   J < (L*LEN*K + R*LEN*DIMB), */

#line 604 "MB02KD.f"
L220:

#line 606 "MB02KD.f"
	ln = len;
#line 607 "MB02KD.f"
	wpos = pdw + 1;

#line 609 "MB02KD.f"
	for (pp = p - 1; pp >= 1; --pp) {
#line 610 "MB02KD.f"
	    ln /= 2;
#line 611 "MB02KD.f"
	    p2 = 1;
#line 612 "MB02KD.f"
	    q2 = ln * kk + 1;
#line 613 "MB02KD.f"
	    r2 = ln / 2 * kk + 1;
#line 614 "MB02KD.f"
	    s2 = r2 + q2 - 1;

#line 616 "MB02KD.f"
	    i__2 = len / (ln << 1) - 1;
#line 616 "MB02KD.f"
	    for (i__ = 0; i__ <= i__2; ++i__) {

#line 618 "MB02KD.f"
		i__1 = kk - 1;
#line 618 "MB02KD.f"
		for (ir = 0; ir <= i__1; ++ir) {
#line 619 "MB02KD.f"
		    t1 = dwork[q2 + ir + j];
#line 620 "MB02KD.f"
		    dwork[q2 + ir + j] = dwork[p2 + ir + j] - t1;
#line 621 "MB02KD.f"
		    dwork[p2 + ir + j] += t1;
#line 622 "MB02KD.f"
		    t1 = dwork[s2 + ir + j];
#line 623 "MB02KD.f"
		    dwork[s2 + ir + j] = dwork[r2 + ir + j] - t1;
#line 624 "MB02KD.f"
		    dwork[r2 + ir + j] += t1;
#line 625 "MB02KD.f"
/* L230: */
#line 625 "MB02KD.f"
		}

#line 627 "MB02KD.f"
		p1 = p2 + kk;
#line 628 "MB02KD.f"
		q1 = p1 + ln * kk;
#line 629 "MB02KD.f"
		r1 = q1 - (kk << 1);
#line 630 "MB02KD.f"
		s1 = r1 + ln * kk;

#line 632 "MB02KD.f"
		i__1 = wpos + ln - 3;
#line 632 "MB02KD.f"
		for (jj = wpos; jj <= i__1; jj += 2) {
#line 633 "MB02KD.f"
		    cf = dwork[jj];
#line 634 "MB02KD.f"
		    sf = dwork[jj + 1];

#line 636 "MB02KD.f"
		    i__3 = kk - 1;
#line 636 "MB02KD.f"
		    for (ir = 0; ir <= i__3; ++ir) {
#line 637 "MB02KD.f"
			t1 = dwork[p1 + ir + j] - dwork[q1 + ir + j];
#line 638 "MB02KD.f"
			t2 = dwork[r1 + ir + j] - dwork[s1 + ir + j];
#line 639 "MB02KD.f"
			dwork[p1 + ir + j] += dwork[q1 + ir + j];
#line 641 "MB02KD.f"
			dwork[r1 + ir + j] += dwork[s1 + ir + j];
#line 643 "MB02KD.f"
			dwork[q1 + ir + j] = cf * t1 + sf * t2;
#line 644 "MB02KD.f"
			dwork[s1 + ir + j] = -cf * t2 + sf * t1;
#line 645 "MB02KD.f"
/* L240: */
#line 645 "MB02KD.f"
		    }

#line 647 "MB02KD.f"
		    p1 += kk;
#line 648 "MB02KD.f"
		    q1 += kk;
#line 649 "MB02KD.f"
		    r1 -= kk;
#line 650 "MB02KD.f"
		    s1 -= kk;
#line 651 "MB02KD.f"
/* L250: */
#line 651 "MB02KD.f"
		}

#line 653 "MB02KD.f"
		p2 += (kk << 1) * ln;
#line 654 "MB02KD.f"
		q2 += (kk << 1) * ln;
#line 655 "MB02KD.f"
		r2 += (kk << 1) * ln;
#line 656 "MB02KD.f"
		s2 += (kk << 1) * ln;
#line 657 "MB02KD.f"
/* L260: */
#line 657 "MB02KD.f"
	    }

#line 659 "MB02KD.f"
	    wpos = wpos + ln - 2;
#line 660 "MB02KD.f"
/* L270: */
#line 660 "MB02KD.f"
	}

#line 662 "MB02KD.f"
	i__2 = len * kk;
#line 662 "MB02KD.f"
	i__1 = kk << 1;
#line 662 "MB02KD.f"
	for (icp = kk + 1; i__1 < 0 ? icp >= i__2 : icp <= i__2; icp += i__1) 
		{
#line 663 "MB02KD.f"
	    icq = icp - kk;

#line 665 "MB02KD.f"
	    i__3 = kk - 1;
#line 665 "MB02KD.f"
	    for (ir = 0; ir <= i__3; ++ir) {
#line 666 "MB02KD.f"
		t1 = dwork[icp + ir + j];
#line 667 "MB02KD.f"
		dwork[icp + ir + j] = dwork[icq + ir + j] - t1;
#line 668 "MB02KD.f"
		dwork[icq + ir + j] += t1;
#line 669 "MB02KD.f"
/* L280: */
#line 669 "MB02KD.f"
	    }

#line 671 "MB02KD.f"
/* L290: */
#line 671 "MB02KD.f"
	}

#line 673 "MB02KD.f"
	j += len * kk;
#line 674 "MB02KD.f"
	if (j == *l * len * *k) {
#line 675 "MB02KD.f"
	    kk = dimb;
#line 676 "MB02KD.f"
	}
#line 677 "MB02KD.f"
	if (j < pc) {
#line 677 "MB02KD.f"
	    goto L220;
#line 677 "MB02KD.f"
	}
/*        END WHILE 220 */

/*        Step 4: Compute a Hadamard like product. */

#line 682 "MB02KD.f"
	i__1 = len - p;
#line 682 "MB02KD.f"
	dcopy_(&i__1, &dwork[pdw + 1], &c__1, &dwork[pdw + 1 + *r__ * len * 
		dimc], &c__1);
#line 683 "MB02KD.f"
	pdw += *r__ * len * dimc;
#line 684 "MB02KD.f"
	scal = *alpha / (doublereal) len;
#line 685 "MB02KD.f"
	p1 = 1;
#line 686 "MB02KD.f"
	r1 = len * *k * *l + 1;
#line 687 "MB02KD.f"
	s1 = r1 + len * dimb * *r__;
#line 688 "MB02KD.f"
	if (ltran) {
#line 689 "MB02KD.f"
	    kk = *l;
#line 690 "MB02KD.f"
	    ll = *k;
#line 691 "MB02KD.f"
	} else {
#line 692 "MB02KD.f"
	    kk = *k;
#line 693 "MB02KD.f"
	    ll = *l;
#line 694 "MB02KD.f"
	}
#line 695 "MB02KD.f"
	i__1 = len * *k;
#line 695 "MB02KD.f"
	i__2 = len * dimb;
#line 695 "MB02KD.f"
	i__3 = len * dimc;
#line 695 "MB02KD.f"
	dgemm_(trans, "No Transpose", &kk, r__, &ll, &scal, &dwork[p1], &i__1,
		 &dwork[r1], &i__2, &c_b9, &dwork[s1], &i__3, (ftnlen)1, (
		ftnlen)12);
#line 698 "MB02KD.f"
	p1 += *k;
#line 699 "MB02KD.f"
	r1 += dimb;
#line 700 "MB02KD.f"
	s1 += dimc;
#line 701 "MB02KD.f"
	i__1 = len * *k;
#line 701 "MB02KD.f"
	i__2 = len * dimb;
#line 701 "MB02KD.f"
	i__3 = len * dimc;
#line 701 "MB02KD.f"
	dgemm_(trans, "No Transpose", &kk, r__, &ll, &scal, &dwork[p1], &i__1,
		 &dwork[r1], &i__2, &c_b9, &dwork[s1], &i__3, (ftnlen)1, (
		ftnlen)12);
#line 704 "MB02KD.f"
	scal /= 2.;
#line 705 "MB02KD.f"
	ln = 1;

#line 707 "MB02KD.f"
	i__1 = p - 1;
#line 707 "MB02KD.f"
	for (pp = 1; pp <= i__1; ++pp) {
#line 708 "MB02KD.f"
	    ln <<= 1;
#line 709 "MB02KD.f"
	    p2 = ((ln << 1) - 1) * *k + 1;
#line 710 "MB02KD.f"
	    r1 = pb + ln * dimb + 1;
#line 711 "MB02KD.f"
	    r2 = pb + ((ln << 1) - 1) * dimb + 1;
#line 712 "MB02KD.f"
	    s1 = pc + ln * dimc + 1;
#line 713 "MB02KD.f"
	    s2 = pc + ((ln << 1) - 1) * dimc + 1;

#line 715 "MB02KD.f"
	    i__2 = (ln + ln / 2) * *k;
#line 715 "MB02KD.f"
	    i__3 = *k;
#line 715 "MB02KD.f"
	    for (p1 = ln * *k + 1; i__3 < 0 ? p1 >= i__2 : p1 <= i__2; p1 += 
		    i__3) {

#line 717 "MB02KD.f"
		i__4 = len * *k * (*l - 1);
#line 717 "MB02KD.f"
		i__5 = len * *k;
#line 717 "MB02KD.f"
		for (j = 0; i__5 < 0 ? j >= i__4 : j <= i__4; j += i__5) {

#line 719 "MB02KD.f"
		    i__6 = p1 + *k - 1;
#line 719 "MB02KD.f"
		    for (i__ = p1; i__ <= i__6; ++i__) {
#line 720 "MB02KD.f"
			t1 = dwork[p2];
#line 721 "MB02KD.f"
			dwork[p2] = dwork[j + i__] - t1;
#line 722 "MB02KD.f"
			dwork[j + i__] += t1;
#line 723 "MB02KD.f"
			++p2;
#line 724 "MB02KD.f"
/* L300: */
#line 724 "MB02KD.f"
		    }

#line 726 "MB02KD.f"
		    p2 += (len - 1) * *k;
#line 727 "MB02KD.f"
/* L310: */
#line 727 "MB02KD.f"
		}

#line 729 "MB02KD.f"
		p2 -= len * *k * *l;
#line 730 "MB02KD.f"
		i__5 = len * *k;
#line 730 "MB02KD.f"
		i__4 = len * dimb;
#line 730 "MB02KD.f"
		i__6 = len * dimc;
#line 730 "MB02KD.f"
		dgemm_(trans, "No Transpose", &kk, r__, &ll, &scal, &dwork[p1]
			, &i__5, &dwork[r1], &i__4, &c_b9, &dwork[s1], &i__6, 
			(ftnlen)1, (ftnlen)12);
#line 733 "MB02KD.f"
		i__5 = len * *k;
#line 733 "MB02KD.f"
		i__4 = len * dimb;
#line 733 "MB02KD.f"
		i__6 = len * dimc;
#line 733 "MB02KD.f"
		dgemm_(trans, "No Transpose", &kk, r__, &ll, &scal, &dwork[p2]
			, &i__5, &dwork[r2], &i__4, &c_b23, &dwork[s1], &i__6,
			 (ftnlen)1, (ftnlen)12);
#line 736 "MB02KD.f"
		i__5 = len * *k;
#line 736 "MB02KD.f"
		i__4 = len * dimb;
#line 736 "MB02KD.f"
		i__6 = len * dimc;
#line 736 "MB02KD.f"
		dgemm_(trans, "No Transpose", &kk, r__, &ll, &scal, &dwork[p1]
			, &i__5, &dwork[r2], &i__4, &c_b9, &dwork[s2], &i__6, 
			(ftnlen)1, (ftnlen)12);
#line 739 "MB02KD.f"
		d__1 = -scal;
#line 739 "MB02KD.f"
		i__5 = len * *k;
#line 739 "MB02KD.f"
		i__4 = len * dimb;
#line 739 "MB02KD.f"
		i__6 = len * dimc;
#line 739 "MB02KD.f"
		dgemm_(trans, "No Transpose", &kk, r__, &ll, &d__1, &dwork[p2]
			, &i__5, &dwork[r1], &i__4, &c_b23, &dwork[s2], &i__6,
			 (ftnlen)1, (ftnlen)12);
#line 742 "MB02KD.f"
		p2 -= *k;
#line 743 "MB02KD.f"
		r1 += dimb;
#line 744 "MB02KD.f"
		r2 -= dimb;
#line 745 "MB02KD.f"
		s1 += dimc;
#line 746 "MB02KD.f"
		s2 -= dimc;
#line 747 "MB02KD.f"
/* L320: */
#line 747 "MB02KD.f"
	    }

#line 749 "MB02KD.f"
/* L330: */
#line 749 "MB02KD.f"
	}

/*        Step 5: Hartley transform with scrambled input. */

#line 753 "MB02KD.f"
	i__1 = pc + len * dimc * *r__;
#line 753 "MB02KD.f"
	i__3 = len * dimc;
#line 753 "MB02KD.f"
	for (j = pc; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {

#line 755 "MB02KD.f"
	    i__2 = len * dimc;
#line 755 "MB02KD.f"
	    i__5 = dimc << 1;
#line 755 "MB02KD.f"
	    for (icp = dimc + 1; i__5 < 0 ? icp >= i__2 : icp <= i__2; icp += 
		    i__5) {
#line 756 "MB02KD.f"
		icq = icp - dimc;

#line 758 "MB02KD.f"
		i__4 = dimc - 1;
#line 758 "MB02KD.f"
		for (ir = 0; ir <= i__4; ++ir) {
#line 759 "MB02KD.f"
		    t1 = dwork[icp + ir + j];
#line 760 "MB02KD.f"
		    dwork[icp + ir + j] = dwork[icq + ir + j] - t1;
#line 761 "MB02KD.f"
		    dwork[icq + ir + j] += t1;
#line 762 "MB02KD.f"
/* L340: */
#line 762 "MB02KD.f"
		}

#line 764 "MB02KD.f"
/* L350: */
#line 764 "MB02KD.f"
	    }

#line 766 "MB02KD.f"
	    ln = 1;
#line 767 "MB02KD.f"
	    wpos = pdw + len - (p << 1) + 1;

#line 769 "MB02KD.f"
	    i__5 = p - 1;
#line 769 "MB02KD.f"
	    for (pp = 1; pp <= i__5; ++pp) {
#line 770 "MB02KD.f"
		ln <<= 1;
#line 771 "MB02KD.f"
		p2 = 1;
#line 772 "MB02KD.f"
		q2 = ln * dimc + 1;
#line 773 "MB02KD.f"
		r2 = ln / 2 * dimc + 1;
#line 774 "MB02KD.f"
		s2 = r2 + q2 - 1;

#line 776 "MB02KD.f"
		i__2 = len / (ln << 1) - 1;
#line 776 "MB02KD.f"
		for (i__ = 0; i__ <= i__2; ++i__) {

#line 778 "MB02KD.f"
		    i__4 = dimc - 1;
#line 778 "MB02KD.f"
		    for (ir = 0; ir <= i__4; ++ir) {
#line 779 "MB02KD.f"
			t1 = dwork[q2 + ir + j];
#line 780 "MB02KD.f"
			dwork[q2 + ir + j] = dwork[p2 + ir + j] - t1;
#line 781 "MB02KD.f"
			dwork[p2 + ir + j] += t1;
#line 782 "MB02KD.f"
			t1 = dwork[s2 + ir + j];
#line 783 "MB02KD.f"
			dwork[s2 + ir + j] = dwork[r2 + ir + j] - t1;
#line 784 "MB02KD.f"
			dwork[r2 + ir + j] += t1;
#line 785 "MB02KD.f"
/* L360: */
#line 785 "MB02KD.f"
		    }

#line 787 "MB02KD.f"
		    p1 = p2 + dimc;
#line 788 "MB02KD.f"
		    q1 = p1 + ln * dimc;
#line 789 "MB02KD.f"
		    r1 = q1 - (dimc << 1);
#line 790 "MB02KD.f"
		    s1 = r1 + ln * dimc;

#line 792 "MB02KD.f"
		    i__4 = wpos + ln - 3;
#line 792 "MB02KD.f"
		    for (jj = wpos; jj <= i__4; jj += 2) {
#line 793 "MB02KD.f"
			cf = dwork[jj];
#line 794 "MB02KD.f"
			sf = dwork[jj + 1];

#line 796 "MB02KD.f"
			i__6 = dimc - 1;
#line 796 "MB02KD.f"
			for (ir = 0; ir <= i__6; ++ir) {
#line 797 "MB02KD.f"
			    t1 = cf * dwork[q1 + ir + j] + sf * dwork[s1 + ir 
				    + j];
#line 798 "MB02KD.f"
			    t2 = -cf * dwork[s1 + ir + j] + sf * dwork[q1 + 
				    ir + j];
#line 799 "MB02KD.f"
			    dwork[q1 + ir + j] = dwork[p1 + ir + j] - t1;
#line 800 "MB02KD.f"
			    dwork[p1 + ir + j] += t1;
#line 801 "MB02KD.f"
			    dwork[s1 + ir + j] = dwork[r1 + ir + j] - t2;
#line 802 "MB02KD.f"
			    dwork[r1 + ir + j] += t2;
#line 803 "MB02KD.f"
/* L370: */
#line 803 "MB02KD.f"
			}

#line 805 "MB02KD.f"
			p1 += dimc;
#line 806 "MB02KD.f"
			q1 += dimc;
#line 807 "MB02KD.f"
			r1 -= dimc;
#line 808 "MB02KD.f"
			s1 -= dimc;
#line 809 "MB02KD.f"
/* L380: */
#line 809 "MB02KD.f"
		    }

#line 811 "MB02KD.f"
		    p2 += (dimc << 1) * ln;
#line 812 "MB02KD.f"
		    q2 += (dimc << 1) * ln;
#line 813 "MB02KD.f"
		    r2 += (dimc << 1) * ln;
#line 814 "MB02KD.f"
		    s2 += (dimc << 1) * ln;
#line 815 "MB02KD.f"
/* L390: */
#line 815 "MB02KD.f"
		}

#line 817 "MB02KD.f"
		wpos = wpos - (ln << 1) + 2;
#line 818 "MB02KD.f"
/* L400: */
#line 818 "MB02KD.f"
	    }

#line 820 "MB02KD.f"
/* L410: */
#line 820 "MB02KD.f"
	}

/*        Step 6: Copy data from workspace to output. */

#line 824 "MB02KD.f"
L420:

#line 826 "MB02KD.f"
	if (ltran) {
#line 827 "MB02KD.f"
	    i__ = nl;
#line 828 "MB02KD.f"
	} else {
#line 829 "MB02KD.f"
	    i__ = mk;
#line 830 "MB02KD.f"
	}

#line 832 "MB02KD.f"
	i__3 = *r__ - 1;
#line 832 "MB02KD.f"
	for (j = 0; j <= i__3; ++j) {
#line 833 "MB02KD.f"
	    daxpy_(&i__, &c_b23, &dwork[pc + j * len * dimc + 1], &c__1, &c__[
		    (j + 1) * c_dim1 + 1], &c__1);
#line 835 "MB02KD.f"
/* L430: */
#line 835 "MB02KD.f"
	}

#line 837 "MB02KD.f"
    }
#line 838 "MB02KD.f"
    dwork[1] = (doublereal) max(1,wrkopt);
#line 839 "MB02KD.f"
    return 0;

/* *** Last line of MB02KD *** */
} /* mb02kd_ */

