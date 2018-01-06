#line 1 "TG01AD.f"
/* TG01AD.f -- translated by f2c (version 20100827).
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

#line 1 "TG01AD.f"
/* Table of constant values */

static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b18 = 10.;
static doublereal c_b56 = .5;

/* Subroutine */ int tg01ad_(char *job, integer *l, integer *n, integer *m, 
	integer *p, doublereal *thresh, doublereal *a, integer *lda, 
	doublereal *e, integer *lde, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *lscale, doublereal *rscale, doublereal 
	*dwork, integer *info, ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, e_dim1, 
	    e_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double d_lg10(doublereal *), d_sign(doublereal *, doublereal *), pow_di(
	    doublereal *, integer *);

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
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lsfmin, lsfmax;


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

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the state dynamics matrix A. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the balanced matrix Dl*A*Dr. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,L). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the descriptor matrix E. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the balanced matrix Dl*E*Dr. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,L). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading L-by-M part of this array must */
/*             contain the input/state matrix B. */
/*             On exit, if M > 0, the leading L-by-M part of this array */
/*             contains the balanced matrix Dl*B. */
/*             The array B is not referenced if M = 0. */

/*     LDB     INTEGER */
/*             The leading dimension of array B. */
/*             LDB >= MAX(1,L) if M > 0 or LDB >= 1 if M = 0. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
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
/*     March 1999. Based on the LAPACK routine DGGBAL. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, July 1999, */
/*     May 2003, March 2004, Jan. 2009. */

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

/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 224 "TG01AD.f"
    /* Parameter adjustments */
#line 224 "TG01AD.f"
    a_dim1 = *lda;
#line 224 "TG01AD.f"
    a_offset = 1 + a_dim1;
#line 224 "TG01AD.f"
    a -= a_offset;
#line 224 "TG01AD.f"
    e_dim1 = *lde;
#line 224 "TG01AD.f"
    e_offset = 1 + e_dim1;
#line 224 "TG01AD.f"
    e -= e_offset;
#line 224 "TG01AD.f"
    b_dim1 = *ldb;
#line 224 "TG01AD.f"
    b_offset = 1 + b_dim1;
#line 224 "TG01AD.f"
    b -= b_offset;
#line 224 "TG01AD.f"
    c_dim1 = *ldc;
#line 224 "TG01AD.f"
    c_offset = 1 + c_dim1;
#line 224 "TG01AD.f"
    c__ -= c_offset;
#line 224 "TG01AD.f"
    --lscale;
#line 224 "TG01AD.f"
    --rscale;
#line 224 "TG01AD.f"
    --dwork;
#line 224 "TG01AD.f"

#line 224 "TG01AD.f"
    /* Function Body */
#line 224 "TG01AD.f"
    *info = 0;
#line 225 "TG01AD.f"
    withb = lsame_(job, "A", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (
	    ftnlen)1, (ftnlen)1);
#line 226 "TG01AD.f"
    withc = lsame_(job, "A", (ftnlen)1, (ftnlen)1) || lsame_(job, "C", (
	    ftnlen)1, (ftnlen)1);

#line 228 "TG01AD.f"
    if (! withb && ! withc && ! lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 230 "TG01AD.f"
	*info = -1;
#line 231 "TG01AD.f"
    } else if (*l < 0) {
#line 232 "TG01AD.f"
	*info = -2;
#line 233 "TG01AD.f"
    } else if (*n < 0) {
#line 234 "TG01AD.f"
	*info = -3;
#line 235 "TG01AD.f"
    } else if (*m < 0) {
#line 236 "TG01AD.f"
	*info = -4;
#line 237 "TG01AD.f"
    } else if (*p < 0) {
#line 238 "TG01AD.f"
	*info = -5;
#line 239 "TG01AD.f"
    } else if (*thresh < 0.) {
#line 240 "TG01AD.f"
	*info = -6;
#line 241 "TG01AD.f"
    } else if (*lda < max(1,*l)) {
#line 242 "TG01AD.f"
	*info = -8;
#line 243 "TG01AD.f"
    } else if (*lde < max(1,*l)) {
#line 244 "TG01AD.f"
	*info = -10;
#line 245 "TG01AD.f"
    } else if (*ldb < 1 || *m > 0 && *ldb < *l) {
#line 246 "TG01AD.f"
	*info = -12;
#line 247 "TG01AD.f"
    } else if (*ldc < max(1,*p)) {
#line 248 "TG01AD.f"
	*info = -14;
#line 249 "TG01AD.f"
    }
#line 250 "TG01AD.f"
    if (*info != 0) {
#line 251 "TG01AD.f"
	i__1 = -(*info);
#line 251 "TG01AD.f"
	xerbla_("TG01AD", &i__1, (ftnlen)6);
#line 252 "TG01AD.f"
	return 0;
#line 253 "TG01AD.f"
    }

/*     Quick return if possible. */

#line 257 "TG01AD.f"
    if (*l == 0 || *n == 0) {
#line 258 "TG01AD.f"
	dum[0] = 1.;
#line 259 "TG01AD.f"
	if (*l > 0) {
#line 260 "TG01AD.f"
	    dcopy_(l, dum, &c__0, &lscale[1], &c__1);
#line 261 "TG01AD.f"
	} else if (*n > 0) {
#line 262 "TG01AD.f"
	    dcopy_(n, dum, &c__0, &rscale[1], &c__1);
#line 263 "TG01AD.f"
	}
#line 264 "TG01AD.f"
	return 0;
#line 265 "TG01AD.f"
    }

/*     Initialize balancing and allocate work storage. */

#line 269 "TG01AD.f"
    kw1 = *n;
#line 270 "TG01AD.f"
    kw2 = kw1 + *l;
#line 271 "TG01AD.f"
    kw3 = kw2 + *l;
#line 272 "TG01AD.f"
    kw4 = kw3 + *n;
#line 273 "TG01AD.f"
    kw5 = kw4 + *l;
#line 274 "TG01AD.f"
    dum[0] = 0.;
#line 275 "TG01AD.f"
    dcopy_(l, dum, &c__0, &lscale[1], &c__1);
#line 276 "TG01AD.f"
    dcopy_(n, dum, &c__0, &rscale[1], &c__1);
#line 277 "TG01AD.f"
    i__1 = (*l + *n) * 3;
#line 277 "TG01AD.f"
    dcopy_(&i__1, dum, &c__0, &dwork[1], &c__1);

/*     Compute right side vector in resulting linear equations. */

#line 281 "TG01AD.f"
    basl = d_lg10(&c_b18);
#line 282 "TG01AD.f"
    i__1 = *l;
#line 282 "TG01AD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 283 "TG01AD.f"
	i__2 = *n;
#line 283 "TG01AD.f"
	for (j = 1; j <= i__2; ++j) {
#line 284 "TG01AD.f"
	    te = (d__1 = e[i__ + j * e_dim1], abs(d__1));
#line 285 "TG01AD.f"
	    ta = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 286 "TG01AD.f"
	    if (ta > *thresh) {
#line 287 "TG01AD.f"
		ta = d_lg10(&ta) / basl;
#line 288 "TG01AD.f"
	    } else {
#line 289 "TG01AD.f"
		ta = 0.;
#line 290 "TG01AD.f"
	    }
#line 291 "TG01AD.f"
	    if (te > *thresh) {
#line 292 "TG01AD.f"
		te = d_lg10(&te) / basl;
#line 293 "TG01AD.f"
	    } else {
#line 294 "TG01AD.f"
		te = 0.;
#line 295 "TG01AD.f"
	    }
#line 296 "TG01AD.f"
	    dwork[i__ + kw4] = dwork[i__ + kw4] - ta - te;
#line 297 "TG01AD.f"
	    dwork[j + kw5] = dwork[j + kw5] - ta - te;
#line 298 "TG01AD.f"
/* L10: */
#line 298 "TG01AD.f"
	}
#line 299 "TG01AD.f"
/* L20: */
#line 299 "TG01AD.f"
    }

#line 301 "TG01AD.f"
    if (*m == 0) {
#line 302 "TG01AD.f"
	withb = FALSE_;
#line 303 "TG01AD.f"
	tb = 0.;
#line 304 "TG01AD.f"
    }
#line 305 "TG01AD.f"
    if (*p == 0) {
#line 306 "TG01AD.f"
	withc = FALSE_;
#line 307 "TG01AD.f"
	tc = 0.;
#line 308 "TG01AD.f"
    }

#line 310 "TG01AD.f"
    if (withb) {
#line 311 "TG01AD.f"
	i__1 = *l;
#line 311 "TG01AD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 312 "TG01AD.f"
	    j = idamax_(m, &b[i__ + b_dim1], ldb);
#line 313 "TG01AD.f"
	    tb = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 314 "TG01AD.f"
	    if (tb > *thresh) {
#line 315 "TG01AD.f"
		tb = d_lg10(&tb) / basl;
#line 316 "TG01AD.f"
		dwork[i__ + kw4] -= tb;
#line 317 "TG01AD.f"
	    }
#line 318 "TG01AD.f"
/* L30: */
#line 318 "TG01AD.f"
	}
#line 319 "TG01AD.f"
    }

#line 321 "TG01AD.f"
    if (withc) {
#line 322 "TG01AD.f"
	i__1 = *n;
#line 322 "TG01AD.f"
	for (j = 1; j <= i__1; ++j) {
#line 323 "TG01AD.f"
	    i__ = idamax_(p, &c__[j * c_dim1 + 1], &c__1);
#line 324 "TG01AD.f"
	    tc = (d__1 = c__[i__ + j * c_dim1], abs(d__1));
#line 325 "TG01AD.f"
	    if (tc > *thresh) {
#line 326 "TG01AD.f"
		tc = d_lg10(&tc) / basl;
#line 327 "TG01AD.f"
		dwork[j + kw5] -= tc;
#line 328 "TG01AD.f"
	    }
#line 329 "TG01AD.f"
/* L40: */
#line 329 "TG01AD.f"
	}
#line 330 "TG01AD.f"
    }

#line 332 "TG01AD.f"
    coef = 1. / (doublereal) (*l + *n);
#line 333 "TG01AD.f"
    coef2 = coef * coef;
#line 334 "TG01AD.f"
    coef5 = coef2 * .5;
#line 335 "TG01AD.f"
    nrp2 = max(*l,*n) + 2;
#line 336 "TG01AD.f"
    beta = 0.;
#line 337 "TG01AD.f"
    it = 1;

/*     Start generalized conjugate gradient iteration. */

#line 341 "TG01AD.f"
L50:

#line 343 "TG01AD.f"
    gamma = ddot_(l, &dwork[kw4 + 1], &c__1, &dwork[kw4 + 1], &c__1) + ddot_(
	    n, &dwork[kw5 + 1], &c__1, &dwork[kw5 + 1], &c__1);

#line 346 "TG01AD.f"
    ew = 0.;
#line 347 "TG01AD.f"
    i__1 = *l;
#line 347 "TG01AD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 348 "TG01AD.f"
	ew += dwork[i__ + kw4];
#line 349 "TG01AD.f"
/* L60: */
#line 349 "TG01AD.f"
    }

#line 351 "TG01AD.f"
    ewc = 0.;
#line 352 "TG01AD.f"
    i__1 = *n;
#line 352 "TG01AD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 353 "TG01AD.f"
	ewc += dwork[i__ + kw5];
#line 354 "TG01AD.f"
/* L70: */
#line 354 "TG01AD.f"
    }

/* Computing 2nd power */
#line 356 "TG01AD.f"
    d__1 = ew;
/* Computing 2nd power */
#line 356 "TG01AD.f"
    d__2 = ewc;
/* Computing 2nd power */
#line 356 "TG01AD.f"
    d__3 = ew - ewc;
#line 356 "TG01AD.f"
    gamma = coef * gamma - coef2 * (d__1 * d__1 + d__2 * d__2) - coef5 * (
	    d__3 * d__3);
#line 358 "TG01AD.f"
    if (gamma == 0.) {
#line 358 "TG01AD.f"
	goto L160;
#line 358 "TG01AD.f"
    }
#line 360 "TG01AD.f"
    if (it != 1) {
#line 360 "TG01AD.f"
	beta = gamma / pgamma;
#line 360 "TG01AD.f"
    }
#line 362 "TG01AD.f"
    t = coef5 * (ewc - ew * 3.);
#line 363 "TG01AD.f"
    tc = coef5 * (ew - ewc * 3.);

#line 365 "TG01AD.f"
    i__1 = *n + *l;
#line 365 "TG01AD.f"
    dscal_(&i__1, &beta, &dwork[1], &c__1);

#line 367 "TG01AD.f"
    daxpy_(l, &coef, &dwork[kw4 + 1], &c__1, &dwork[kw1 + 1], &c__1);
#line 368 "TG01AD.f"
    daxpy_(n, &coef, &dwork[kw5 + 1], &c__1, &dwork[1], &c__1);

#line 370 "TG01AD.f"
    i__1 = *n;
#line 370 "TG01AD.f"
    for (j = 1; j <= i__1; ++j) {
#line 371 "TG01AD.f"
	dwork[j] += tc;
#line 372 "TG01AD.f"
/* L80: */
#line 372 "TG01AD.f"
    }

#line 374 "TG01AD.f"
    i__1 = *l;
#line 374 "TG01AD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 375 "TG01AD.f"
	dwork[i__ + kw1] += t;
#line 376 "TG01AD.f"
/* L90: */
#line 376 "TG01AD.f"
    }

/*     Apply matrix to vector. */

#line 380 "TG01AD.f"
    i__1 = *l;
#line 380 "TG01AD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 381 "TG01AD.f"
	kount = 0;
#line 382 "TG01AD.f"
	sum = 0.;
#line 383 "TG01AD.f"
	i__2 = *n;
#line 383 "TG01AD.f"
	for (j = 1; j <= i__2; ++j) {
#line 384 "TG01AD.f"
	    if ((d__1 = a[i__ + j * a_dim1], abs(d__1)) > *thresh) {
#line 385 "TG01AD.f"
		++kount;
#line 386 "TG01AD.f"
		sum += dwork[j];
#line 387 "TG01AD.f"
	    }
#line 388 "TG01AD.f"
	    if ((d__1 = e[i__ + j * e_dim1], abs(d__1)) > *thresh) {
#line 389 "TG01AD.f"
		++kount;
#line 390 "TG01AD.f"
		sum += dwork[j];
#line 391 "TG01AD.f"
	    }
#line 392 "TG01AD.f"
/* L100: */
#line 392 "TG01AD.f"
	}
#line 393 "TG01AD.f"
	if (withb) {
#line 394 "TG01AD.f"
	    j = idamax_(m, &b[i__ + b_dim1], ldb);
#line 395 "TG01AD.f"
	    if ((d__1 = b[i__ + j * b_dim1], abs(d__1)) > *thresh) {
#line 395 "TG01AD.f"
		++kount;
#line 395 "TG01AD.f"
	    }
#line 396 "TG01AD.f"
	}
#line 397 "TG01AD.f"
	dwork[i__ + kw2] = (doublereal) kount * dwork[i__ + kw1] + sum;
#line 398 "TG01AD.f"
/* L110: */
#line 398 "TG01AD.f"
    }

#line 400 "TG01AD.f"
    i__1 = *n;
#line 400 "TG01AD.f"
    for (j = 1; j <= i__1; ++j) {
#line 401 "TG01AD.f"
	kount = 0;
#line 402 "TG01AD.f"
	sum = 0.;
#line 403 "TG01AD.f"
	i__2 = *l;
#line 403 "TG01AD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 404 "TG01AD.f"
	    if ((d__1 = a[i__ + j * a_dim1], abs(d__1)) > *thresh) {
#line 405 "TG01AD.f"
		++kount;
#line 406 "TG01AD.f"
		sum += dwork[i__ + kw1];
#line 407 "TG01AD.f"
	    }
#line 408 "TG01AD.f"
	    if ((d__1 = e[i__ + j * e_dim1], abs(d__1)) > *thresh) {
#line 409 "TG01AD.f"
		++kount;
#line 410 "TG01AD.f"
		sum += dwork[i__ + kw1];
#line 411 "TG01AD.f"
	    }
#line 412 "TG01AD.f"
/* L120: */
#line 412 "TG01AD.f"
	}
#line 413 "TG01AD.f"
	if (withc) {
#line 414 "TG01AD.f"
	    i__ = idamax_(p, &c__[j * c_dim1 + 1], &c__1);
#line 415 "TG01AD.f"
	    if ((d__1 = c__[i__ + j * c_dim1], abs(d__1)) > *thresh) {
#line 415 "TG01AD.f"
		++kount;
#line 415 "TG01AD.f"
	    }
#line 416 "TG01AD.f"
	}
#line 417 "TG01AD.f"
	dwork[j + kw3] = (doublereal) kount * dwork[j] + sum;
#line 418 "TG01AD.f"
/* L130: */
#line 418 "TG01AD.f"
    }

#line 420 "TG01AD.f"
    sum = ddot_(l, &dwork[kw1 + 1], &c__1, &dwork[kw2 + 1], &c__1) + ddot_(n, 
	    &dwork[1], &c__1, &dwork[kw3 + 1], &c__1);
#line 422 "TG01AD.f"
    alpha = gamma / sum;

/*     Determine correction to current iteration. */

#line 426 "TG01AD.f"
    cmax = 0.;
#line 427 "TG01AD.f"
    i__1 = *l;
#line 427 "TG01AD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 428 "TG01AD.f"
	cor = alpha * dwork[i__ + kw1];
#line 429 "TG01AD.f"
	if (abs(cor) > cmax) {
#line 429 "TG01AD.f"
	    cmax = abs(cor);
#line 429 "TG01AD.f"
	}
#line 431 "TG01AD.f"
	lscale[i__] += cor;
#line 432 "TG01AD.f"
/* L140: */
#line 432 "TG01AD.f"
    }

#line 434 "TG01AD.f"
    i__1 = *n;
#line 434 "TG01AD.f"
    for (j = 1; j <= i__1; ++j) {
#line 435 "TG01AD.f"
	cor = alpha * dwork[j];
#line 436 "TG01AD.f"
	if (abs(cor) > cmax) {
#line 436 "TG01AD.f"
	    cmax = abs(cor);
#line 436 "TG01AD.f"
	}
#line 438 "TG01AD.f"
	rscale[j] += cor;
#line 439 "TG01AD.f"
/* L150: */
#line 439 "TG01AD.f"
    }
#line 440 "TG01AD.f"
    if (cmax < .5) {
#line 440 "TG01AD.f"
	goto L160;
#line 440 "TG01AD.f"
    }

#line 443 "TG01AD.f"
    d__1 = -alpha;
#line 443 "TG01AD.f"
    daxpy_(l, &d__1, &dwork[kw2 + 1], &c__1, &dwork[kw4 + 1], &c__1);
#line 444 "TG01AD.f"
    d__1 = -alpha;
#line 444 "TG01AD.f"
    daxpy_(n, &d__1, &dwork[kw3 + 1], &c__1, &dwork[kw5 + 1], &c__1);

#line 446 "TG01AD.f"
    pgamma = gamma;
#line 447 "TG01AD.f"
    ++it;
#line 448 "TG01AD.f"
    if (it <= nrp2) {
#line 448 "TG01AD.f"
	goto L50;
#line 448 "TG01AD.f"
    }

/*     End generalized conjugate gradient iteration. */

#line 453 "TG01AD.f"
L160:
#line 454 "TG01AD.f"
    sfmin = dlamch_("Safe minimum", (ftnlen)12);
#line 455 "TG01AD.f"
    sfmax = 1. / sfmin;
#line 456 "TG01AD.f"
    lsfmin = (integer) (d_lg10(&sfmin) / basl + 1.);
#line 457 "TG01AD.f"
    lsfmax = (integer) (d_lg10(&sfmax) / basl);

/*     Compute left diagonal scaling matrix. */

#line 461 "TG01AD.f"
    i__1 = *l;
#line 461 "TG01AD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 462 "TG01AD.f"
	irab = idamax_(n, &a[i__ + a_dim1], lda);
#line 463 "TG01AD.f"
	rab = (d__1 = a[i__ + irab * a_dim1], abs(d__1));
#line 464 "TG01AD.f"
	irab = idamax_(n, &e[i__ + e_dim1], lde);
/* Computing MAX */
#line 465 "TG01AD.f"
	d__2 = rab, d__3 = (d__1 = e[i__ + irab * e_dim1], abs(d__1));
#line 465 "TG01AD.f"
	rab = max(d__2,d__3);
#line 466 "TG01AD.f"
	if (withb) {
#line 467 "TG01AD.f"
	    irab = idamax_(m, &b[i__ + b_dim1], ldb);
/* Computing MAX */
#line 468 "TG01AD.f"
	    d__2 = rab, d__3 = (d__1 = b[i__ + irab * b_dim1], abs(d__1));
#line 468 "TG01AD.f"
	    rab = max(d__2,d__3);
#line 469 "TG01AD.f"
	}
#line 470 "TG01AD.f"
	d__1 = rab + sfmin;
#line 470 "TG01AD.f"
	lrab = (integer) (d_lg10(&d__1) / basl + 1.);
#line 471 "TG01AD.f"
	ir = (integer) (lscale[i__] + d_sign(&c_b56, &lscale[i__]));
/* Computing MIN */
#line 472 "TG01AD.f"
	i__2 = max(ir,lsfmin), i__2 = min(i__2,lsfmax), i__3 = lsfmax - lrab;
#line 472 "TG01AD.f"
	ir = min(i__2,i__3);
#line 473 "TG01AD.f"
	lscale[i__] = pow_di(&c_b18, &ir);
#line 474 "TG01AD.f"
/* L170: */
#line 474 "TG01AD.f"
    }

/*     Compute right diagonal scaling matrix. */

#line 478 "TG01AD.f"
    i__1 = *n;
#line 478 "TG01AD.f"
    for (j = 1; j <= i__1; ++j) {
#line 479 "TG01AD.f"
	icab = idamax_(l, &a[j * a_dim1 + 1], &c__1);
#line 480 "TG01AD.f"
	cab = (d__1 = a[icab + j * a_dim1], abs(d__1));
#line 481 "TG01AD.f"
	icab = idamax_(l, &e[j * e_dim1 + 1], &c__1);
/* Computing MAX */
#line 482 "TG01AD.f"
	d__2 = cab, d__3 = (d__1 = e[icab + j * e_dim1], abs(d__1));
#line 482 "TG01AD.f"
	cab = max(d__2,d__3);
#line 483 "TG01AD.f"
	if (withc) {
#line 484 "TG01AD.f"
	    icab = idamax_(p, &c__[j * c_dim1 + 1], &c__1);
/* Computing MAX */
#line 485 "TG01AD.f"
	    d__2 = cab, d__3 = (d__1 = c__[icab + j * c_dim1], abs(d__1));
#line 485 "TG01AD.f"
	    cab = max(d__2,d__3);
#line 486 "TG01AD.f"
	}
#line 487 "TG01AD.f"
	d__1 = cab + sfmin;
#line 487 "TG01AD.f"
	lcab = (integer) (d_lg10(&d__1) / basl + 1.);
#line 488 "TG01AD.f"
	jc = (integer) (rscale[j] + d_sign(&c_b56, &rscale[j]));
/* Computing MIN */
#line 489 "TG01AD.f"
	i__2 = max(jc,lsfmin), i__2 = min(i__2,lsfmax), i__3 = lsfmax - lcab;
#line 489 "TG01AD.f"
	jc = min(i__2,i__3);
#line 490 "TG01AD.f"
	rscale[j] = pow_di(&c_b18, &jc);
#line 491 "TG01AD.f"
/* L180: */
#line 491 "TG01AD.f"
    }

/*     Row scaling of matrices A, E and B. */

#line 495 "TG01AD.f"
    i__1 = *l;
#line 495 "TG01AD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 496 "TG01AD.f"
	dscal_(n, &lscale[i__], &a[i__ + a_dim1], lda);
#line 497 "TG01AD.f"
	dscal_(n, &lscale[i__], &e[i__ + e_dim1], lde);
#line 498 "TG01AD.f"
	if (withb) {
#line 498 "TG01AD.f"
	    dscal_(m, &lscale[i__], &b[i__ + b_dim1], ldb);
#line 498 "TG01AD.f"
	}
#line 500 "TG01AD.f"
/* L190: */
#line 500 "TG01AD.f"
    }

/*     Column scaling of matrices A, E and C. */

#line 504 "TG01AD.f"
    i__1 = *n;
#line 504 "TG01AD.f"
    for (j = 1; j <= i__1; ++j) {
#line 505 "TG01AD.f"
	dscal_(l, &rscale[j], &a[j * a_dim1 + 1], &c__1);
#line 506 "TG01AD.f"
	dscal_(l, &rscale[j], &e[j * e_dim1 + 1], &c__1);
#line 507 "TG01AD.f"
	if (withc) {
#line 507 "TG01AD.f"
	    dscal_(p, &rscale[j], &c__[j * c_dim1 + 1], &c__1);
#line 507 "TG01AD.f"
	}
#line 509 "TG01AD.f"
/* L200: */
#line 509 "TG01AD.f"
    }

#line 511 "TG01AD.f"
    return 0;
/* *** Last line of TG01AD *** */
} /* tg01ad_ */

