#line 1 "TG01AZ.f"
/* TG01AZ.f -- translated by f2c (version 20100827).
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

#line 1 "TG01AZ.f"
/* Table of constant values */

static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b18 = 10.;
static doublereal c_b56 = .5;

/* Subroutine */ int tg01az_(char *job, integer *l, integer *n, integer *m, 
	integer *p, doublereal *thresh, doublecomplex *a, integer *lda, 
	doublecomplex *e, integer *lde, doublecomplex *b, integer *ldb, 
	doublecomplex *c__, integer *ldc, doublereal *lscale, doublereal *
	rscale, doublereal *dwork, integer *info, ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, e_dim1, 
	    e_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double d_lg10(doublereal *), d_imag(doublecomplex *), z_abs(doublecomplex 
	    *), d_sign(doublereal *, doublereal *), pow_di(doublereal *, 
	    integer *);

    /* Local variables */
    static integer i__, j;
    static doublereal t;
    static integer jc;
    static doublereal ta, tb, tc, te;
    static integer ir;
    static doublereal ew;
    static integer it, kw1, kw2, kw3, kw4, kw5;
    static doublereal cab, rab, ewc, cor, dum[1], sum;
    static integer nrp2, icab, lcab;
    static doublereal beta, coef;
    static integer irab, lrab;
    static doublereal basl, cmax;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal coef2, coef5, gamma, alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal sfmin;
    static logical withb, withc;
    static doublereal sfmax;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer kount;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal pgamma;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    static integer lsfmin;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    static integer lsfmax;


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

/*     To balance the matrices of the system pencil */

/*             S =  ( A  B ) - lambda ( E  0 ) :=  Q - lambda Z, */
/*                  ( C  0 )          ( 0  0 ) */

/*     corresponding to the descriptor triple (A-lambda E,B,C), */
/*     by balancing. This involves diagonal similarity transformations */
/*     (Dl*A*Dr - lambda Dl*E*Dr, Dl*B, C*Dr) applied to the system */
/*     (A-lambda E,B,C) to make the rows and columns of system pencil */
/*     matrices */

/*                  diag(Dl,I) * S * diag(Dr,I) */

/*     as close in norm as possible. Balancing may reduce the 1-norms */
/*     of the matrices of the system pencil S. */

/*     The balancing can be performed optionally on the following */
/*     particular system pencils */

/*              S = A-lambda E, */

/*              S = ( A-lambda E  B ),    or */

/*              S = ( A-lambda E ). */
/*                  (     C      ) */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Indicates which matrices are involved in balancing, as */
/*             follows: */
/*             = 'A':  All matrices are involved in balancing; */
/*             = 'B':  B, A and E matrices are involved in balancing; */
/*             = 'C':  C, A and E matrices are involved in balancing; */
/*             = 'N':  B and C matrices are not involved in balancing. */

/*     Input/Output Parameters */

/*     L       (input) INTEGER */
/*             The number of rows of matrices A, B, and E.  L >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of matrices A, E, and C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of matrix B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of rows of matrix C.  P >= 0. */

/*     THRESH  (input) DOUBLE PRECISION */
/*             Threshold value for magnitude of elements: */
/*             elements with magnitude less than or equal to */
/*             THRESH are ignored for balancing. THRESH >= 0. */
/*             The magnitude is computed as the sum of the absolute */
/*             values of the real and imaginary parts. */

/*     A       (input/output) COMPLEX*16 array, dimension (LDA,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the state dynamics matrix A. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the balanced matrix Dl*A*Dr. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,L). */

/*     E       (input/output) COMPLEX*16 array, dimension (LDE,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the descriptor matrix E. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the balanced matrix Dl*E*Dr. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,L). */

/*     B       (input/output) COMPLEX*16 array, dimension (LDB,M) */
/*             On entry, the leading L-by-M part of this array must */
/*             contain the input/state matrix B. */
/*             On exit, if M > 0, the leading L-by-M part of this array */
/*             contains the balanced matrix Dl*B. */
/*             The array B is not referenced if M = 0. */

/*     LDB     INTEGER */
/*             The leading dimension of array B. */
/*             LDB >= MAX(1,L) if M > 0 or LDB >= 1 if M = 0. */

/*     C       (input/output) COMPLEX*16 array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C. */
/*             On exit, if P > 0, the leading P-by-N part of this array */
/*             contains the balanced matrix C*Dr. */
/*             The array C is not referenced if P = 0. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     LSCALE  (output) DOUBLE PRECISION array, dimension (L) */
/*             The scaling factors applied to S from left.  If Dl(j) is */
/*             the scaling factor applied to row j, then */
/*             SCALE(j) = Dl(j), for j = 1,...,L. */

/*     RSCALE  (output) DOUBLE PRECISION array, dimension (N) */
/*             The scaling factors applied to S from right.  If Dr(j) is */
/*             the scaling factor applied to column j, then */
/*             SCALE(j) = Dr(j), for j = 1,...,N. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (3*(L+N)) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit. */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Balancing consists of applying a diagonal similarity */
/*     transformation */
/*                            -1 */
/*                  diag(Dl,I)  * S * diag(Dr,I) */

/*     to make the 1-norms of each row of the first L rows of S and its */
/*     corresponding N columns nearly equal. */

/*     Information about the diagonal matrices Dl and Dr are returned in */
/*     the vectors LSCALE and RSCALE, respectively. */

/*     REFERENCES */

/*     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J., */
/*         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A., */
/*         Ostrouchov, S., and Sorensen, D. */
/*         LAPACK Users' Guide: Second Edition. */
/*         SIAM, Philadelphia, 1995. */

/*     [2] R.C. Ward, R. C. */
/*         Balancing the generalized eigenvalue problem. */
/*         SIAM J. Sci. Stat. Comp. 2 (1981), 141-152. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen. */
/*     March 1999. */
/*     Complex version: V. Sima, Research Institute for Informatics, */
/*     Bucharest, Nov. 2008. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Balancing, eigenvalue, matrix algebra, matrix operations, */
/*     similarity transformation. */

/*  ********************************************************************* */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */

/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 234 "TG01AZ.f"
    /* Parameter adjustments */
#line 234 "TG01AZ.f"
    a_dim1 = *lda;
#line 234 "TG01AZ.f"
    a_offset = 1 + a_dim1;
#line 234 "TG01AZ.f"
    a -= a_offset;
#line 234 "TG01AZ.f"
    e_dim1 = *lde;
#line 234 "TG01AZ.f"
    e_offset = 1 + e_dim1;
#line 234 "TG01AZ.f"
    e -= e_offset;
#line 234 "TG01AZ.f"
    b_dim1 = *ldb;
#line 234 "TG01AZ.f"
    b_offset = 1 + b_dim1;
#line 234 "TG01AZ.f"
    b -= b_offset;
#line 234 "TG01AZ.f"
    c_dim1 = *ldc;
#line 234 "TG01AZ.f"
    c_offset = 1 + c_dim1;
#line 234 "TG01AZ.f"
    c__ -= c_offset;
#line 234 "TG01AZ.f"
    --lscale;
#line 234 "TG01AZ.f"
    --rscale;
#line 234 "TG01AZ.f"
    --dwork;
#line 234 "TG01AZ.f"

#line 234 "TG01AZ.f"
    /* Function Body */
#line 234 "TG01AZ.f"
    *info = 0;
#line 235 "TG01AZ.f"
    withb = lsame_(job, "A", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (
	    ftnlen)1, (ftnlen)1);
#line 236 "TG01AZ.f"
    withc = lsame_(job, "A", (ftnlen)1, (ftnlen)1) || lsame_(job, "C", (
	    ftnlen)1, (ftnlen)1);

#line 238 "TG01AZ.f"
    if (! withb && ! withc && ! lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 240 "TG01AZ.f"
	*info = -1;
#line 241 "TG01AZ.f"
    } else if (*l < 0) {
#line 242 "TG01AZ.f"
	*info = -2;
#line 243 "TG01AZ.f"
    } else if (*n < 0) {
#line 244 "TG01AZ.f"
	*info = -3;
#line 245 "TG01AZ.f"
    } else if (*m < 0) {
#line 246 "TG01AZ.f"
	*info = -4;
#line 247 "TG01AZ.f"
    } else if (*p < 0) {
#line 248 "TG01AZ.f"
	*info = -5;
#line 249 "TG01AZ.f"
    } else if (*thresh < 0.) {
#line 250 "TG01AZ.f"
	*info = -6;
#line 251 "TG01AZ.f"
    } else if (*lda < max(1,*l)) {
#line 252 "TG01AZ.f"
	*info = -8;
#line 253 "TG01AZ.f"
    } else if (*lde < max(1,*l)) {
#line 254 "TG01AZ.f"
	*info = -10;
#line 255 "TG01AZ.f"
    } else if (*ldb < 1 || *m > 0 && *ldb < *l) {
#line 256 "TG01AZ.f"
	*info = -12;
#line 257 "TG01AZ.f"
    } else if (*ldc < max(1,*p)) {
#line 258 "TG01AZ.f"
	*info = -14;
#line 259 "TG01AZ.f"
    }
#line 260 "TG01AZ.f"
    if (*info != 0) {
#line 261 "TG01AZ.f"
	i__1 = -(*info);
#line 261 "TG01AZ.f"
	xerbla_("TG01AZ", &i__1, (ftnlen)6);
#line 262 "TG01AZ.f"
	return 0;
#line 263 "TG01AZ.f"
    }

/*     Quick return if possible. */

#line 267 "TG01AZ.f"
    if (*l == 0 || *n == 0) {
#line 268 "TG01AZ.f"
	dum[0] = 1.;
#line 269 "TG01AZ.f"
	if (*l > 0) {
#line 270 "TG01AZ.f"
	    dcopy_(l, dum, &c__0, &lscale[1], &c__1);
#line 271 "TG01AZ.f"
	} else if (*n > 0) {
#line 272 "TG01AZ.f"
	    dcopy_(n, dum, &c__0, &rscale[1], &c__1);
#line 273 "TG01AZ.f"
	}
#line 274 "TG01AZ.f"
	return 0;
#line 275 "TG01AZ.f"
    }

/*     Initialize balancing and allocate work storage. */

#line 279 "TG01AZ.f"
    kw1 = *n;
#line 280 "TG01AZ.f"
    kw2 = kw1 + *l;
#line 281 "TG01AZ.f"
    kw3 = kw2 + *l;
#line 282 "TG01AZ.f"
    kw4 = kw3 + *n;
#line 283 "TG01AZ.f"
    kw5 = kw4 + *l;
#line 284 "TG01AZ.f"
    dum[0] = 0.;
#line 285 "TG01AZ.f"
    dcopy_(l, dum, &c__0, &lscale[1], &c__1);
#line 286 "TG01AZ.f"
    dcopy_(n, dum, &c__0, &rscale[1], &c__1);
#line 287 "TG01AZ.f"
    i__1 = (*l + *n) * 3;
#line 287 "TG01AZ.f"
    dcopy_(&i__1, dum, &c__0, &dwork[1], &c__1);

/*     Compute right side vector in resulting linear equations. */

#line 291 "TG01AZ.f"
    basl = d_lg10(&c_b18);
#line 292 "TG01AZ.f"
    i__1 = *l;
#line 292 "TG01AZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 293 "TG01AZ.f"
	i__2 = *n;
#line 293 "TG01AZ.f"
	for (j = 1; j <= i__2; ++j) {
#line 294 "TG01AZ.f"
	    i__3 = i__ + j * e_dim1;
#line 294 "TG01AZ.f"
	    te = (d__1 = e[i__3].r, abs(d__1)) + (d__2 = d_imag(&e[i__ + j * 
		    e_dim1]), abs(d__2));
#line 295 "TG01AZ.f"
	    i__3 = i__ + j * a_dim1;
#line 295 "TG01AZ.f"
	    ta = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ + j * 
		    a_dim1]), abs(d__2));
#line 296 "TG01AZ.f"
	    if (ta > *thresh) {
#line 297 "TG01AZ.f"
		ta = d_lg10(&ta) / basl;
#line 298 "TG01AZ.f"
	    } else {
#line 299 "TG01AZ.f"
		ta = 0.;
#line 300 "TG01AZ.f"
	    }
#line 301 "TG01AZ.f"
	    if (te > *thresh) {
#line 302 "TG01AZ.f"
		te = d_lg10(&te) / basl;
#line 303 "TG01AZ.f"
	    } else {
#line 304 "TG01AZ.f"
		te = 0.;
#line 305 "TG01AZ.f"
	    }
#line 306 "TG01AZ.f"
	    dwork[i__ + kw4] = dwork[i__ + kw4] - ta - te;
#line 307 "TG01AZ.f"
	    dwork[j + kw5] = dwork[j + kw5] - ta - te;
#line 308 "TG01AZ.f"
/* L10: */
#line 308 "TG01AZ.f"
	}
#line 309 "TG01AZ.f"
/* L20: */
#line 309 "TG01AZ.f"
    }

#line 311 "TG01AZ.f"
    if (*m == 0) {
#line 312 "TG01AZ.f"
	withb = FALSE_;
#line 313 "TG01AZ.f"
	tb = 0.;
#line 314 "TG01AZ.f"
    }
#line 315 "TG01AZ.f"
    if (*p == 0) {
#line 316 "TG01AZ.f"
	withc = FALSE_;
#line 317 "TG01AZ.f"
	tc = 0.;
#line 318 "TG01AZ.f"
    }

#line 320 "TG01AZ.f"
    if (withb) {
#line 321 "TG01AZ.f"
	i__1 = *l;
#line 321 "TG01AZ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 322 "TG01AZ.f"
	    j = izamax_(m, &b[i__ + b_dim1], ldb);
#line 323 "TG01AZ.f"
	    i__2 = i__ + j * b_dim1;
#line 323 "TG01AZ.f"
	    tb = (d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[i__ + j * 
		    b_dim1]), abs(d__2));
#line 324 "TG01AZ.f"
	    if (tb > *thresh) {
#line 325 "TG01AZ.f"
		tb = d_lg10(&tb) / basl;
#line 326 "TG01AZ.f"
		dwork[i__ + kw4] -= tb;
#line 327 "TG01AZ.f"
	    }
#line 328 "TG01AZ.f"
/* L30: */
#line 328 "TG01AZ.f"
	}
#line 329 "TG01AZ.f"
    }

#line 331 "TG01AZ.f"
    if (withc) {
#line 332 "TG01AZ.f"
	i__1 = *n;
#line 332 "TG01AZ.f"
	for (j = 1; j <= i__1; ++j) {
#line 333 "TG01AZ.f"
	    i__ = izamax_(p, &c__[j * c_dim1 + 1], &c__1);
#line 334 "TG01AZ.f"
	    i__2 = i__ + j * c_dim1;
#line 334 "TG01AZ.f"
	    tc = (d__1 = c__[i__2].r, abs(d__1)) + (d__2 = d_imag(&c__[i__ + 
		    j * c_dim1]), abs(d__2));
#line 335 "TG01AZ.f"
	    if (tc > *thresh) {
#line 336 "TG01AZ.f"
		tc = d_lg10(&tc) / basl;
#line 337 "TG01AZ.f"
		dwork[j + kw5] -= tc;
#line 338 "TG01AZ.f"
	    }
#line 339 "TG01AZ.f"
/* L40: */
#line 339 "TG01AZ.f"
	}
#line 340 "TG01AZ.f"
    }

#line 342 "TG01AZ.f"
    coef = 1. / (doublereal) (*l + *n);
#line 343 "TG01AZ.f"
    coef2 = coef * coef;
#line 344 "TG01AZ.f"
    coef5 = coef2 * .5;
#line 345 "TG01AZ.f"
    nrp2 = max(*l,*n) + 2;
#line 346 "TG01AZ.f"
    beta = 0.;
#line 347 "TG01AZ.f"
    it = 1;

/*     Start generalized conjugate gradient iteration. */

#line 351 "TG01AZ.f"
L50:

#line 353 "TG01AZ.f"
    gamma = ddot_(l, &dwork[kw4 + 1], &c__1, &dwork[kw4 + 1], &c__1) + ddot_(
	    n, &dwork[kw5 + 1], &c__1, &dwork[kw5 + 1], &c__1);

#line 356 "TG01AZ.f"
    ew = 0.;
#line 357 "TG01AZ.f"
    i__1 = *l;
#line 357 "TG01AZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 358 "TG01AZ.f"
	ew += dwork[i__ + kw4];
#line 359 "TG01AZ.f"
/* L60: */
#line 359 "TG01AZ.f"
    }

#line 361 "TG01AZ.f"
    ewc = 0.;
#line 362 "TG01AZ.f"
    i__1 = *n;
#line 362 "TG01AZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 363 "TG01AZ.f"
	ewc += dwork[i__ + kw5];
#line 364 "TG01AZ.f"
/* L70: */
#line 364 "TG01AZ.f"
    }

/* Computing 2nd power */
#line 366 "TG01AZ.f"
    d__1 = ew;
/* Computing 2nd power */
#line 366 "TG01AZ.f"
    d__2 = ewc;
/* Computing 2nd power */
#line 366 "TG01AZ.f"
    d__3 = ew - ewc;
#line 366 "TG01AZ.f"
    gamma = coef * gamma - coef2 * (d__1 * d__1 + d__2 * d__2) - coef5 * (
	    d__3 * d__3);
#line 368 "TG01AZ.f"
    if (gamma == 0.) {
#line 368 "TG01AZ.f"
	goto L160;
#line 368 "TG01AZ.f"
    }
#line 370 "TG01AZ.f"
    if (it != 1) {
#line 370 "TG01AZ.f"
	beta = gamma / pgamma;
#line 370 "TG01AZ.f"
    }
#line 372 "TG01AZ.f"
    t = coef5 * (ewc - ew * 3.);
#line 373 "TG01AZ.f"
    tc = coef5 * (ew - ewc * 3.);

#line 375 "TG01AZ.f"
    i__1 = *n + *l;
#line 375 "TG01AZ.f"
    dscal_(&i__1, &beta, &dwork[1], &c__1);

#line 377 "TG01AZ.f"
    daxpy_(l, &coef, &dwork[kw4 + 1], &c__1, &dwork[kw1 + 1], &c__1);
#line 378 "TG01AZ.f"
    daxpy_(n, &coef, &dwork[kw5 + 1], &c__1, &dwork[1], &c__1);

#line 380 "TG01AZ.f"
    i__1 = *n;
#line 380 "TG01AZ.f"
    for (j = 1; j <= i__1; ++j) {
#line 381 "TG01AZ.f"
	dwork[j] += tc;
#line 382 "TG01AZ.f"
/* L80: */
#line 382 "TG01AZ.f"
    }

#line 384 "TG01AZ.f"
    i__1 = *l;
#line 384 "TG01AZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 385 "TG01AZ.f"
	dwork[i__ + kw1] += t;
#line 386 "TG01AZ.f"
/* L90: */
#line 386 "TG01AZ.f"
    }

/*     Apply matrix to vector. */

#line 390 "TG01AZ.f"
    i__1 = *l;
#line 390 "TG01AZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 391 "TG01AZ.f"
	kount = 0;
#line 392 "TG01AZ.f"
	sum = 0.;
#line 393 "TG01AZ.f"
	i__2 = *n;
#line 393 "TG01AZ.f"
	for (j = 1; j <= i__2; ++j) {
#line 394 "TG01AZ.f"
	    i__3 = i__ + j * a_dim1;
#line 394 "TG01AZ.f"
	    if ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ + j * 
		    a_dim1]), abs(d__2)) > *thresh) {
#line 395 "TG01AZ.f"
		++kount;
#line 396 "TG01AZ.f"
		sum += dwork[j];
#line 397 "TG01AZ.f"
	    }
#line 398 "TG01AZ.f"
	    i__3 = i__ + j * e_dim1;
#line 398 "TG01AZ.f"
	    if ((d__1 = e[i__3].r, abs(d__1)) + (d__2 = d_imag(&e[i__ + j * 
		    e_dim1]), abs(d__2)) > *thresh) {
#line 399 "TG01AZ.f"
		++kount;
#line 400 "TG01AZ.f"
		sum += dwork[j];
#line 401 "TG01AZ.f"
	    }
#line 402 "TG01AZ.f"
/* L100: */
#line 402 "TG01AZ.f"
	}
#line 403 "TG01AZ.f"
	if (withb) {
#line 404 "TG01AZ.f"
	    j = izamax_(m, &b[i__ + b_dim1], ldb);
#line 405 "TG01AZ.f"
	    i__2 = i__ + j * b_dim1;
#line 405 "TG01AZ.f"
	    if ((d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[i__ + j * 
		    b_dim1]), abs(d__2)) > *thresh) {
#line 405 "TG01AZ.f"
		++kount;
#line 405 "TG01AZ.f"
	    }
#line 406 "TG01AZ.f"
	}
#line 407 "TG01AZ.f"
	dwork[i__ + kw2] = (doublereal) kount * dwork[i__ + kw1] + sum;
#line 408 "TG01AZ.f"
/* L110: */
#line 408 "TG01AZ.f"
    }

#line 410 "TG01AZ.f"
    i__1 = *n;
#line 410 "TG01AZ.f"
    for (j = 1; j <= i__1; ++j) {
#line 411 "TG01AZ.f"
	kount = 0;
#line 412 "TG01AZ.f"
	sum = 0.;
#line 413 "TG01AZ.f"
	i__2 = *l;
#line 413 "TG01AZ.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 414 "TG01AZ.f"
	    i__3 = i__ + j * a_dim1;
#line 414 "TG01AZ.f"
	    if ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ + j * 
		    a_dim1]), abs(d__2)) > *thresh) {
#line 415 "TG01AZ.f"
		++kount;
#line 416 "TG01AZ.f"
		sum += dwork[i__ + kw1];
#line 417 "TG01AZ.f"
	    }
#line 418 "TG01AZ.f"
	    i__3 = i__ + j * e_dim1;
#line 418 "TG01AZ.f"
	    if ((d__1 = e[i__3].r, abs(d__1)) + (d__2 = d_imag(&e[i__ + j * 
		    e_dim1]), abs(d__2)) > *thresh) {
#line 419 "TG01AZ.f"
		++kount;
#line 420 "TG01AZ.f"
		sum += dwork[i__ + kw1];
#line 421 "TG01AZ.f"
	    }
#line 422 "TG01AZ.f"
/* L120: */
#line 422 "TG01AZ.f"
	}
#line 423 "TG01AZ.f"
	if (withc) {
#line 424 "TG01AZ.f"
	    i__ = izamax_(p, &c__[j * c_dim1 + 1], &c__1);
#line 425 "TG01AZ.f"
	    i__2 = i__ + j * c_dim1;
#line 425 "TG01AZ.f"
	    if ((d__1 = c__[i__2].r, abs(d__1)) + (d__2 = d_imag(&c__[i__ + j 
		    * c_dim1]), abs(d__2)) > *thresh) {
#line 425 "TG01AZ.f"
		++kount;
#line 425 "TG01AZ.f"
	    }
#line 426 "TG01AZ.f"
	}
#line 427 "TG01AZ.f"
	dwork[j + kw3] = (doublereal) kount * dwork[j] + sum;
#line 428 "TG01AZ.f"
/* L130: */
#line 428 "TG01AZ.f"
    }

#line 430 "TG01AZ.f"
    sum = ddot_(l, &dwork[kw1 + 1], &c__1, &dwork[kw2 + 1], &c__1) + ddot_(n, 
	    &dwork[1], &c__1, &dwork[kw3 + 1], &c__1);
#line 432 "TG01AZ.f"
    alpha = gamma / sum;

/*     Determine correction to current iteration. */

#line 436 "TG01AZ.f"
    cmax = 0.;
#line 437 "TG01AZ.f"
    i__1 = *l;
#line 437 "TG01AZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 438 "TG01AZ.f"
	cor = alpha * dwork[i__ + kw1];
#line 439 "TG01AZ.f"
	if (abs(cor) > cmax) {
#line 439 "TG01AZ.f"
	    cmax = abs(cor);
#line 439 "TG01AZ.f"
	}
#line 441 "TG01AZ.f"
	lscale[i__] += cor;
#line 442 "TG01AZ.f"
/* L140: */
#line 442 "TG01AZ.f"
    }

#line 444 "TG01AZ.f"
    i__1 = *n;
#line 444 "TG01AZ.f"
    for (j = 1; j <= i__1; ++j) {
#line 445 "TG01AZ.f"
	cor = alpha * dwork[j];
#line 446 "TG01AZ.f"
	if (abs(cor) > cmax) {
#line 446 "TG01AZ.f"
	    cmax = abs(cor);
#line 446 "TG01AZ.f"
	}
#line 448 "TG01AZ.f"
	rscale[j] += cor;
#line 449 "TG01AZ.f"
/* L150: */
#line 449 "TG01AZ.f"
    }
#line 450 "TG01AZ.f"
    if (cmax < .5) {
#line 450 "TG01AZ.f"
	goto L160;
#line 450 "TG01AZ.f"
    }

#line 453 "TG01AZ.f"
    d__1 = -alpha;
#line 453 "TG01AZ.f"
    daxpy_(l, &d__1, &dwork[kw2 + 1], &c__1, &dwork[kw4 + 1], &c__1);
#line 454 "TG01AZ.f"
    d__1 = -alpha;
#line 454 "TG01AZ.f"
    daxpy_(n, &d__1, &dwork[kw3 + 1], &c__1, &dwork[kw5 + 1], &c__1);

#line 456 "TG01AZ.f"
    pgamma = gamma;
#line 457 "TG01AZ.f"
    ++it;
#line 458 "TG01AZ.f"
    if (it <= nrp2) {
#line 458 "TG01AZ.f"
	goto L50;
#line 458 "TG01AZ.f"
    }

/*     End generalized conjugate gradient iteration. */

#line 463 "TG01AZ.f"
L160:
#line 464 "TG01AZ.f"
    sfmin = dlamch_("Safe minimum", (ftnlen)12);
#line 465 "TG01AZ.f"
    sfmax = 1. / sfmin;
#line 466 "TG01AZ.f"
    lsfmin = (integer) (d_lg10(&sfmin) / basl + 1.);
#line 467 "TG01AZ.f"
    lsfmax = (integer) (d_lg10(&sfmax) / basl);

/*     Compute left diagonal scaling matrix. */

#line 471 "TG01AZ.f"
    i__1 = *l;
#line 471 "TG01AZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 472 "TG01AZ.f"
	irab = izamax_(n, &a[i__ + a_dim1], lda);
#line 473 "TG01AZ.f"
	rab = z_abs(&a[i__ + irab * a_dim1]);
#line 474 "TG01AZ.f"
	irab = izamax_(n, &e[i__ + e_dim1], lde);
/* Computing MAX */
#line 475 "TG01AZ.f"
	d__1 = rab, d__2 = z_abs(&e[i__ + irab * e_dim1]);
#line 475 "TG01AZ.f"
	rab = max(d__1,d__2);
#line 476 "TG01AZ.f"
	if (withb) {
#line 477 "TG01AZ.f"
	    irab = izamax_(m, &b[i__ + b_dim1], ldb);
/* Computing MAX */
#line 478 "TG01AZ.f"
	    d__1 = rab, d__2 = z_abs(&b[i__ + irab * b_dim1]);
#line 478 "TG01AZ.f"
	    rab = max(d__1,d__2);
#line 479 "TG01AZ.f"
	}
#line 480 "TG01AZ.f"
	d__1 = rab + sfmin;
#line 480 "TG01AZ.f"
	lrab = (integer) (d_lg10(&d__1) / basl + 1.);
#line 481 "TG01AZ.f"
	ir = (integer) (lscale[i__] + d_sign(&c_b56, &lscale[i__]));
/* Computing MIN */
#line 482 "TG01AZ.f"
	i__2 = max(ir,lsfmin), i__2 = min(i__2,lsfmax), i__3 = lsfmax - lrab;
#line 482 "TG01AZ.f"
	ir = min(i__2,i__3);
#line 483 "TG01AZ.f"
	lscale[i__] = pow_di(&c_b18, &ir);
#line 484 "TG01AZ.f"
/* L170: */
#line 484 "TG01AZ.f"
    }

/*     Compute right diagonal scaling matrix. */

#line 488 "TG01AZ.f"
    i__1 = *n;
#line 488 "TG01AZ.f"
    for (j = 1; j <= i__1; ++j) {
#line 489 "TG01AZ.f"
	icab = izamax_(l, &a[j * a_dim1 + 1], &c__1);
#line 490 "TG01AZ.f"
	cab = z_abs(&a[icab + j * a_dim1]);
#line 491 "TG01AZ.f"
	icab = izamax_(l, &e[j * e_dim1 + 1], &c__1);
/* Computing MAX */
#line 492 "TG01AZ.f"
	d__1 = cab, d__2 = z_abs(&e[icab + j * e_dim1]);
#line 492 "TG01AZ.f"
	cab = max(d__1,d__2);
#line 493 "TG01AZ.f"
	if (withc) {
#line 494 "TG01AZ.f"
	    icab = izamax_(p, &c__[j * c_dim1 + 1], &c__1);
/* Computing MAX */
#line 495 "TG01AZ.f"
	    d__1 = cab, d__2 = z_abs(&c__[icab + j * c_dim1]);
#line 495 "TG01AZ.f"
	    cab = max(d__1,d__2);
#line 496 "TG01AZ.f"
	}
#line 497 "TG01AZ.f"
	d__1 = cab + sfmin;
#line 497 "TG01AZ.f"
	lcab = (integer) (d_lg10(&d__1) / basl + 1.);
#line 498 "TG01AZ.f"
	jc = (integer) (rscale[j] + d_sign(&c_b56, &rscale[j]));
/* Computing MIN */
#line 499 "TG01AZ.f"
	i__2 = max(jc,lsfmin), i__2 = min(i__2,lsfmax), i__3 = lsfmax - lcab;
#line 499 "TG01AZ.f"
	jc = min(i__2,i__3);
#line 500 "TG01AZ.f"
	rscale[j] = pow_di(&c_b18, &jc);
#line 501 "TG01AZ.f"
/* L180: */
#line 501 "TG01AZ.f"
    }

/*     Row scaling of matrices A, E and B. */

#line 505 "TG01AZ.f"
    i__1 = *l;
#line 505 "TG01AZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 506 "TG01AZ.f"
	zdscal_(n, &lscale[i__], &a[i__ + a_dim1], lda);
#line 507 "TG01AZ.f"
	zdscal_(n, &lscale[i__], &e[i__ + e_dim1], lde);
#line 508 "TG01AZ.f"
	if (withb) {
#line 508 "TG01AZ.f"
	    zdscal_(m, &lscale[i__], &b[i__ + b_dim1], ldb);
#line 508 "TG01AZ.f"
	}
#line 510 "TG01AZ.f"
/* L190: */
#line 510 "TG01AZ.f"
    }

/*     Column scaling of matrices A, E and C. */

#line 514 "TG01AZ.f"
    i__1 = *n;
#line 514 "TG01AZ.f"
    for (j = 1; j <= i__1; ++j) {
#line 515 "TG01AZ.f"
	zdscal_(l, &rscale[j], &a[j * a_dim1 + 1], &c__1);
#line 516 "TG01AZ.f"
	zdscal_(l, &rscale[j], &e[j * e_dim1 + 1], &c__1);
#line 517 "TG01AZ.f"
	if (withc) {
#line 517 "TG01AZ.f"
	    zdscal_(p, &rscale[j], &c__[j * c_dim1 + 1], &c__1);
#line 517 "TG01AZ.f"
	}
#line 519 "TG01AZ.f"
/* L200: */
#line 519 "TG01AZ.f"
    }

#line 521 "TG01AZ.f"
    return 0;
/* *** Last line of TG01AZ *** */
} /* tg01az_ */

