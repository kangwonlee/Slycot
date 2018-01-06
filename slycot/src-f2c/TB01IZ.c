#line 1 "TB01IZ.f"
/* TB01IZ.f -- translated by f2c (version 20100827).
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

#line 1 "TB01IZ.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int tb01iz_(char *job, integer *n, integer *m, integer *p, 
	doublereal *maxred, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *c__, integer *ldc, doublereal *scale, 
	integer *info, ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_imag(doublecomplex *), z_abs(doublecomplex *);

    /* Local variables */
    static doublereal f, g;
    static integer i__, j;
    static doublereal s, ca, co, ra, ro;
    static integer ica, ira;
    static doublereal sred;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical withb, withc;
    static doublereal snorm, sfmin1, sfmin2, sfmax1, sfmax2;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    extern integer izamax_(integer *, doublecomplex *, integer *);
    static logical noconv;
    static doublereal maxnrm;
    extern doublereal dzasum_(integer *, doublecomplex *, integer *);


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

/*     To reduce the 1-norm of a system matrix */

/*             S =  ( A  B ) */
/*                  ( C  0 ) */

/*     corresponding to the triple (A,B,C), by balancing. This involves */
/*     a diagonal similarity transformation inv(D)*A*D applied */
/*     iteratively to A to make the rows and columns of */
/*                           -1 */
/*                  diag(D,I)  * S * diag(D,I) */

/*     as close in norm as possible. */

/*     The balancing can be performed optionally on the following */
/*     particular system matrices */

/*              S = A,    S = ( A  B )    or    S = ( A ) */
/*                                                  ( C ) */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Indicates which matrices are involved in balancing, as */
/*             follows: */
/*             = 'A':  All matrices are involved in balancing; */
/*             = 'B':  B and A matrices are involved in balancing; */
/*             = 'C':  C and A matrices are involved in balancing; */
/*             = 'N':  B and C matrices are not involved in balancing. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A, the number of rows of matrix B */
/*             and the number of columns of matrix C. */
/*             N represents the dimension of the state vector.  N >= 0. */

/*     M       (input) INTEGER. */
/*             The number of columns of matrix B. */
/*             M represents the dimension of input vector.  M >= 0. */

/*     P       (input) INTEGER. */
/*             The number of rows of matrix C. */
/*             P represents the dimension of output vector.  P >= 0. */

/*     MAXRED  (input/output) DOUBLE PRECISION */
/*             On entry, the maximum allowed reduction in the 1-norm of */
/*             S (in an iteration) if zero rows or columns are */
/*             encountered. */
/*             If MAXRED > 0.0, MAXRED must be larger than one (to enable */
/*             the norm reduction). */
/*             If MAXRED <= 0.0, then the value 10.0 for MAXRED is */
/*             used. */
/*             On exit, if the 1-norm of the given matrix S is non-zero, */
/*             the ratio between the 1-norm of the given matrix and the */
/*             1-norm of the balanced matrix. */

/*     A       (input/output) COMPLEX*16 array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the system state matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the balanced matrix inv(D)*A*D. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     B       (input/output) COMPLEX*16 array, dimension (LDB,M) */
/*             On entry, if M > 0, the leading N-by-M part of this array */
/*             must contain the system input matrix B. */
/*             On exit, if M > 0, the leading N-by-M part of this array */
/*             contains the balanced matrix inv(D)*B. */
/*             The array B is not referenced if M = 0. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= MAX(1,N) if M > 0. */
/*             LDB >= 1        if M = 0. */

/*     C       (input/output) COMPLEX*16 array, dimension (LDC,N) */
/*             On entry, if P > 0, the leading P-by-N part of this array */
/*             must contain the system output matrix C. */
/*             On exit, if P > 0, the leading P-by-N part of this array */
/*             contains the balanced matrix C*D. */
/*             The array C is not referenced if P = 0. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= MAX(1,P). */

/*     SCALE   (output) DOUBLE PRECISION array, dimension (N) */
/*             The scaling factors applied to S.  If D(j) is the scaling */
/*             factor applied to row and column j, then SCALE(j) = D(j), */
/*             for j = 1,...,N. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit. */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Balancing consists of applying a diagonal similarity */
/*     transformation */
/*                           -1 */
/*                  diag(D,I)  * S * diag(D,I) */

/*     to make the 1-norms of each row of the first N rows of S and its */
/*     corresponding column nearly equal. */

/*     Information about the diagonal matrix D is returned in the vector */
/*     SCALE. */

/*     REFERENCES */

/*     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J., */
/*         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A., */
/*         Ostrouchov, S., and Sorensen, D. */
/*         LAPACK Users' Guide: Second Edition. */
/*         SIAM, Philadelphia, 1995. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Jan. 1998. */
/*     Complex version: V. Sima, Research Institute for Informatics, */
/*     Bucharest, Nov. 2008. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Balancing, eigenvalue, matrix algebra, matrix operations, */
/*     similarity transformation. */

/*  ********************************************************************* */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the scalar input arguments. */

#line 213 "TB01IZ.f"
    /* Parameter adjustments */
#line 213 "TB01IZ.f"
    a_dim1 = *lda;
#line 213 "TB01IZ.f"
    a_offset = 1 + a_dim1;
#line 213 "TB01IZ.f"
    a -= a_offset;
#line 213 "TB01IZ.f"
    b_dim1 = *ldb;
#line 213 "TB01IZ.f"
    b_offset = 1 + b_dim1;
#line 213 "TB01IZ.f"
    b -= b_offset;
#line 213 "TB01IZ.f"
    c_dim1 = *ldc;
#line 213 "TB01IZ.f"
    c_offset = 1 + c_dim1;
#line 213 "TB01IZ.f"
    c__ -= c_offset;
#line 213 "TB01IZ.f"
    --scale;
#line 213 "TB01IZ.f"

#line 213 "TB01IZ.f"
    /* Function Body */
#line 213 "TB01IZ.f"
    *info = 0;
#line 214 "TB01IZ.f"
    withb = lsame_(job, "A", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (
	    ftnlen)1, (ftnlen)1);
#line 215 "TB01IZ.f"
    withc = lsame_(job, "A", (ftnlen)1, (ftnlen)1) || lsame_(job, "C", (
	    ftnlen)1, (ftnlen)1);

#line 217 "TB01IZ.f"
    if (! withb && ! withc && ! lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 219 "TB01IZ.f"
	*info = -1;
#line 220 "TB01IZ.f"
    } else if (*n < 0) {
#line 221 "TB01IZ.f"
	*info = -2;
#line 222 "TB01IZ.f"
    } else if (*m < 0) {
#line 223 "TB01IZ.f"
	*info = -3;
#line 224 "TB01IZ.f"
    } else if (*p < 0) {
#line 225 "TB01IZ.f"
	*info = -4;
#line 226 "TB01IZ.f"
    } else if (*maxred > 0. && *maxred < 1.) {
#line 227 "TB01IZ.f"
	*info = -5;
#line 228 "TB01IZ.f"
    } else if (*lda < max(1,*n)) {
#line 229 "TB01IZ.f"
	*info = -7;
#line 230 "TB01IZ.f"
    } else if (*m > 0 && *ldb < max(1,*n) || *m == 0 && *ldb < 1) {
#line 232 "TB01IZ.f"
	*info = -9;
#line 233 "TB01IZ.f"
    } else if (*ldc < max(1,*p)) {
#line 234 "TB01IZ.f"
	*info = -11;
#line 235 "TB01IZ.f"
    }
#line 236 "TB01IZ.f"
    if (*info != 0) {
#line 237 "TB01IZ.f"
	i__1 = -(*info);
#line 237 "TB01IZ.f"
	xerbla_("TB01IZ", &i__1, (ftnlen)6);
#line 238 "TB01IZ.f"
	return 0;
#line 239 "TB01IZ.f"
    }

#line 241 "TB01IZ.f"
    if (*n == 0) {
#line 241 "TB01IZ.f"
	return 0;
#line 241 "TB01IZ.f"
    }

/*     Compute the 1-norm of the required part of matrix S and exit if */
/*     it is zero. */

#line 247 "TB01IZ.f"
    snorm = 0.;

#line 249 "TB01IZ.f"
    i__1 = *n;
#line 249 "TB01IZ.f"
    for (j = 1; j <= i__1; ++j) {
#line 250 "TB01IZ.f"
	scale[j] = 1.;
#line 251 "TB01IZ.f"
	co = dzasum_(n, &a[j * a_dim1 + 1], &c__1);
#line 252 "TB01IZ.f"
	if (withc && *p > 0) {
#line 252 "TB01IZ.f"
	    co += dzasum_(p, &c__[j * c_dim1 + 1], &c__1);
#line 252 "TB01IZ.f"
	}
#line 254 "TB01IZ.f"
	snorm = max(snorm,co);
#line 255 "TB01IZ.f"
/* L10: */
#line 255 "TB01IZ.f"
    }

#line 257 "TB01IZ.f"
    if (withb) {

#line 259 "TB01IZ.f"
	i__1 = *m;
#line 259 "TB01IZ.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 260 "TB01IZ.f"
	    d__1 = snorm, d__2 = dzasum_(n, &b[j * b_dim1 + 1], &c__1);
#line 260 "TB01IZ.f"
	    snorm = max(d__1,d__2);
#line 261 "TB01IZ.f"
/* L20: */
#line 261 "TB01IZ.f"
	}

#line 263 "TB01IZ.f"
    }

#line 265 "TB01IZ.f"
    if (snorm == 0.) {
#line 265 "TB01IZ.f"
	return 0;
#line 265 "TB01IZ.f"
    }

/*     Set some machine parameters and the maximum reduction in the */
/*     1-norm of S if zero rows or columns are encountered. */

#line 271 "TB01IZ.f"
    sfmin1 = dlamch_("S", (ftnlen)1) / dlamch_("P", (ftnlen)1);
#line 272 "TB01IZ.f"
    sfmax1 = 1. / sfmin1;
#line 273 "TB01IZ.f"
    sfmin2 = sfmin1 * 10.;
#line 274 "TB01IZ.f"
    sfmax2 = 1. / sfmin2;

#line 276 "TB01IZ.f"
    sred = *maxred;
#line 277 "TB01IZ.f"
    if (sred <= 0.) {
#line 277 "TB01IZ.f"
	sred = 10.;
#line 277 "TB01IZ.f"
    }

/* Computing MAX */
#line 279 "TB01IZ.f"
    d__1 = snorm / sred;
#line 279 "TB01IZ.f"
    maxnrm = max(d__1,sfmin1);

/*     Balance the matrix. */

/*     Iterative loop for norm reduction. */

#line 285 "TB01IZ.f"
L30:
#line 286 "TB01IZ.f"
    noconv = FALSE_;

#line 288 "TB01IZ.f"
    i__1 = *n;
#line 288 "TB01IZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 289 "TB01IZ.f"
	co = 0.;
#line 290 "TB01IZ.f"
	ro = 0.;

#line 292 "TB01IZ.f"
	i__2 = *n;
#line 292 "TB01IZ.f"
	for (j = 1; j <= i__2; ++j) {
#line 293 "TB01IZ.f"
	    if (j == i__) {
#line 293 "TB01IZ.f"
		goto L40;
#line 293 "TB01IZ.f"
	    }
#line 295 "TB01IZ.f"
	    i__3 = j + i__ * a_dim1;
#line 295 "TB01IZ.f"
	    co += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[j + i__ * 
		    a_dim1]), abs(d__2));
#line 296 "TB01IZ.f"
	    i__3 = i__ + j * a_dim1;
#line 296 "TB01IZ.f"
	    ro += (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[i__ + j * 
		    a_dim1]), abs(d__2));
#line 297 "TB01IZ.f"
L40:
#line 297 "TB01IZ.f"
	    ;
#line 297 "TB01IZ.f"
	}

#line 299 "TB01IZ.f"
	ica = izamax_(n, &a[i__ * a_dim1 + 1], &c__1);
#line 300 "TB01IZ.f"
	ca = z_abs(&a[ica + i__ * a_dim1]);
#line 301 "TB01IZ.f"
	ira = izamax_(n, &a[i__ + a_dim1], lda);
#line 302 "TB01IZ.f"
	ra = z_abs(&a[i__ + ira * a_dim1]);

#line 304 "TB01IZ.f"
	if (withc && *p > 0) {
#line 305 "TB01IZ.f"
	    co += dzasum_(p, &c__[i__ * c_dim1 + 1], &c__1);
#line 306 "TB01IZ.f"
	    ica = izamax_(p, &c__[i__ * c_dim1 + 1], &c__1);
/* Computing MAX */
#line 307 "TB01IZ.f"
	    d__1 = ca, d__2 = z_abs(&c__[ica + i__ * c_dim1]);
#line 307 "TB01IZ.f"
	    ca = max(d__1,d__2);
#line 308 "TB01IZ.f"
	}

#line 310 "TB01IZ.f"
	if (withb && *m > 0) {
#line 311 "TB01IZ.f"
	    ro += dzasum_(m, &b[i__ + b_dim1], ldb);
#line 312 "TB01IZ.f"
	    ira = izamax_(m, &b[i__ + b_dim1], ldb);
/* Computing MAX */
#line 313 "TB01IZ.f"
	    d__1 = ra, d__2 = z_abs(&b[i__ + ira * b_dim1]);
#line 313 "TB01IZ.f"
	    ra = max(d__1,d__2);
#line 314 "TB01IZ.f"
	}

/*        Special case of zero CO and/or RO. */

#line 318 "TB01IZ.f"
	if (co == 0. && ro == 0.) {
#line 318 "TB01IZ.f"
	    goto L90;
#line 318 "TB01IZ.f"
	}
#line 320 "TB01IZ.f"
	if (co == 0.) {
#line 321 "TB01IZ.f"
	    if (ro <= maxnrm) {
#line 321 "TB01IZ.f"
		goto L90;
#line 321 "TB01IZ.f"
	    }
#line 323 "TB01IZ.f"
	    co = maxnrm;
#line 324 "TB01IZ.f"
	}
#line 325 "TB01IZ.f"
	if (ro == 0.) {
#line 326 "TB01IZ.f"
	    if (co <= maxnrm) {
#line 326 "TB01IZ.f"
		goto L90;
#line 326 "TB01IZ.f"
	    }
#line 328 "TB01IZ.f"
	    ro = maxnrm;
#line 329 "TB01IZ.f"
	}

/*        Guard against zero CO or RO due to underflow. */

#line 333 "TB01IZ.f"
	g = ro / 10.;
#line 334 "TB01IZ.f"
	f = 1.;
#line 335 "TB01IZ.f"
	s = co + ro;
#line 336 "TB01IZ.f"
L50:
/* Computing MAX */
#line 337 "TB01IZ.f"
	d__1 = max(f,co);
/* Computing MIN */
#line 337 "TB01IZ.f"
	d__2 = min(ro,g);
#line 337 "TB01IZ.f"
	if (co >= g || max(d__1,ca) >= sfmax2 || min(d__2,ra) <= sfmin2) {
#line 337 "TB01IZ.f"
	    goto L60;
#line 337 "TB01IZ.f"
	}
#line 339 "TB01IZ.f"
	f *= 10.;
#line 340 "TB01IZ.f"
	co *= 10.;
#line 341 "TB01IZ.f"
	ca *= 10.;
#line 342 "TB01IZ.f"
	g /= 10.;
#line 343 "TB01IZ.f"
	ro /= 10.;
#line 344 "TB01IZ.f"
	ra /= 10.;
#line 345 "TB01IZ.f"
	goto L50;

#line 347 "TB01IZ.f"
L60:
#line 348 "TB01IZ.f"
	g = co / 10.;
#line 349 "TB01IZ.f"
L70:
/* Computing MIN */
#line 350 "TB01IZ.f"
	d__1 = min(f,co), d__1 = min(d__1,g);
#line 350 "TB01IZ.f"
	if (g < ro || max(ro,ra) >= sfmax2 || min(d__1,ca) <= sfmin2) {
#line 350 "TB01IZ.f"
	    goto L80;
#line 350 "TB01IZ.f"
	}
#line 352 "TB01IZ.f"
	f /= 10.;
#line 353 "TB01IZ.f"
	co /= 10.;
#line 354 "TB01IZ.f"
	ca /= 10.;
#line 355 "TB01IZ.f"
	g /= 10.;
#line 356 "TB01IZ.f"
	ro *= 10.;
#line 357 "TB01IZ.f"
	ra *= 10.;
#line 358 "TB01IZ.f"
	goto L70;

/*        Now balance. */

#line 362 "TB01IZ.f"
L80:
#line 363 "TB01IZ.f"
	if (co + ro >= s * .95) {
#line 363 "TB01IZ.f"
	    goto L90;
#line 363 "TB01IZ.f"
	}
#line 365 "TB01IZ.f"
	if (f < 1. && scale[i__] < 1.) {
#line 366 "TB01IZ.f"
	    if (f * scale[i__] <= sfmin1) {
#line 366 "TB01IZ.f"
		goto L90;
#line 366 "TB01IZ.f"
	    }
#line 368 "TB01IZ.f"
	}
#line 369 "TB01IZ.f"
	if (f > 1. && scale[i__] > 1.) {
#line 370 "TB01IZ.f"
	    if (scale[i__] >= sfmax1 / f) {
#line 370 "TB01IZ.f"
		goto L90;
#line 370 "TB01IZ.f"
	    }
#line 372 "TB01IZ.f"
	}
#line 373 "TB01IZ.f"
	g = 1. / f;
#line 374 "TB01IZ.f"
	scale[i__] *= f;
#line 375 "TB01IZ.f"
	noconv = TRUE_;

#line 377 "TB01IZ.f"
	zdscal_(n, &g, &a[i__ + a_dim1], lda);
#line 378 "TB01IZ.f"
	zdscal_(n, &f, &a[i__ * a_dim1 + 1], &c__1);
#line 379 "TB01IZ.f"
	if (*m > 0) {
#line 379 "TB01IZ.f"
	    zdscal_(m, &g, &b[i__ + b_dim1], ldb);
#line 379 "TB01IZ.f"
	}
#line 380 "TB01IZ.f"
	if (*p > 0) {
#line 380 "TB01IZ.f"
	    zdscal_(p, &f, &c__[i__ * c_dim1 + 1], &c__1);
#line 380 "TB01IZ.f"
	}

#line 382 "TB01IZ.f"
L90:
#line 382 "TB01IZ.f"
	;
#line 382 "TB01IZ.f"
    }

#line 384 "TB01IZ.f"
    if (noconv) {
#line 384 "TB01IZ.f"
	goto L30;
#line 384 "TB01IZ.f"
    }

/*     Set the norm reduction parameter. */

#line 389 "TB01IZ.f"
    *maxred = snorm;
#line 390 "TB01IZ.f"
    snorm = 0.;

#line 392 "TB01IZ.f"
    i__1 = *n;
#line 392 "TB01IZ.f"
    for (j = 1; j <= i__1; ++j) {
#line 393 "TB01IZ.f"
	co = dzasum_(n, &a[j * a_dim1 + 1], &c__1);
#line 394 "TB01IZ.f"
	if (withc && *p > 0) {
#line 394 "TB01IZ.f"
	    co += dzasum_(p, &c__[j * c_dim1 + 1], &c__1);
#line 394 "TB01IZ.f"
	}
#line 396 "TB01IZ.f"
	snorm = max(snorm,co);
#line 397 "TB01IZ.f"
/* L100: */
#line 397 "TB01IZ.f"
    }

#line 399 "TB01IZ.f"
    if (withb) {

#line 401 "TB01IZ.f"
	i__1 = *m;
#line 401 "TB01IZ.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 402 "TB01IZ.f"
	    d__1 = snorm, d__2 = dzasum_(n, &b[j * b_dim1 + 1], &c__1);
#line 402 "TB01IZ.f"
	    snorm = max(d__1,d__2);
#line 403 "TB01IZ.f"
/* L110: */
#line 403 "TB01IZ.f"
	}

#line 405 "TB01IZ.f"
    }
#line 406 "TB01IZ.f"
    *maxred /= snorm;
#line 407 "TB01IZ.f"
    return 0;
/* *** Last line of TB01IZ *** */
} /* tb01iz_ */

