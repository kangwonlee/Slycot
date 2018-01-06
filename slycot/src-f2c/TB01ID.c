#line 1 "TB01ID.f"
/* TB01ID.f -- translated by f2c (version 20100827).
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

#line 1 "TB01ID.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int tb01id_(char *job, integer *n, integer *m, integer *p, 
	doublereal *maxred, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *scale, 
	integer *info, ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static doublereal f, g;
    static integer i__, j;
    static doublereal s, ca, co, ra, ro;
    static integer ica, ira;
    static doublereal sred;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern doublereal dasum_(integer *, doublereal *, integer *);
    static logical withb, withc;
    static doublereal snorm, sfmin1, sfmin2, sfmax1, sfmax2;
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical noconv;
    static doublereal maxnrm;


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

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the system state matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the balanced matrix inv(D)*A*D. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, if M > 0, the leading N-by-M part of this array */
/*             must contain the system input matrix B. */
/*             On exit, if M > 0, the leading N-by-M part of this array */
/*             contains the balanced matrix inv(D)*B. */
/*             The array B is not referenced if M = 0. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= MAX(1,N) if M > 0. */
/*             LDB >= 1        if M = 0. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
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

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Jan. 1998. */
/*     This subroutine is based on LAPACK routine DGEBAL, and routine */
/*     BALABC (A. Varga, German Aerospace Research Establishment, DLR). */

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
/*     .. Executable Statements .. */

/*     Test the scalar input arguments. */

#line 206 "TB01ID.f"
    /* Parameter adjustments */
#line 206 "TB01ID.f"
    a_dim1 = *lda;
#line 206 "TB01ID.f"
    a_offset = 1 + a_dim1;
#line 206 "TB01ID.f"
    a -= a_offset;
#line 206 "TB01ID.f"
    b_dim1 = *ldb;
#line 206 "TB01ID.f"
    b_offset = 1 + b_dim1;
#line 206 "TB01ID.f"
    b -= b_offset;
#line 206 "TB01ID.f"
    c_dim1 = *ldc;
#line 206 "TB01ID.f"
    c_offset = 1 + c_dim1;
#line 206 "TB01ID.f"
    c__ -= c_offset;
#line 206 "TB01ID.f"
    --scale;
#line 206 "TB01ID.f"

#line 206 "TB01ID.f"
    /* Function Body */
#line 206 "TB01ID.f"
    *info = 0;
#line 207 "TB01ID.f"
    withb = lsame_(job, "A", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (
	    ftnlen)1, (ftnlen)1);
#line 208 "TB01ID.f"
    withc = lsame_(job, "A", (ftnlen)1, (ftnlen)1) || lsame_(job, "C", (
	    ftnlen)1, (ftnlen)1);

#line 210 "TB01ID.f"
    if (! withb && ! withc && ! lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 212 "TB01ID.f"
	*info = -1;
#line 213 "TB01ID.f"
    } else if (*n < 0) {
#line 214 "TB01ID.f"
	*info = -2;
#line 215 "TB01ID.f"
    } else if (*m < 0) {
#line 216 "TB01ID.f"
	*info = -3;
#line 217 "TB01ID.f"
    } else if (*p < 0) {
#line 218 "TB01ID.f"
	*info = -4;
#line 219 "TB01ID.f"
    } else if (*maxred > 0. && *maxred < 1.) {
#line 220 "TB01ID.f"
	*info = -5;
#line 221 "TB01ID.f"
    } else if (*lda < max(1,*n)) {
#line 222 "TB01ID.f"
	*info = -7;
#line 223 "TB01ID.f"
    } else if (*m > 0 && *ldb < max(1,*n) || *m == 0 && *ldb < 1) {
#line 225 "TB01ID.f"
	*info = -9;
#line 226 "TB01ID.f"
    } else if (*ldc < max(1,*p)) {
#line 227 "TB01ID.f"
	*info = -11;
#line 228 "TB01ID.f"
    }
#line 229 "TB01ID.f"
    if (*info != 0) {
#line 230 "TB01ID.f"
	i__1 = -(*info);
#line 230 "TB01ID.f"
	xerbla_("TB01ID", &i__1, (ftnlen)6);
#line 231 "TB01ID.f"
	return 0;
#line 232 "TB01ID.f"
    }

#line 234 "TB01ID.f"
    if (*n == 0) {
#line 234 "TB01ID.f"
	return 0;
#line 234 "TB01ID.f"
    }

/*     Compute the 1-norm of the required part of matrix S and exit if */
/*     it is zero. */

#line 240 "TB01ID.f"
    snorm = 0.;

#line 242 "TB01ID.f"
    i__1 = *n;
#line 242 "TB01ID.f"
    for (j = 1; j <= i__1; ++j) {
#line 243 "TB01ID.f"
	scale[j] = 1.;
#line 244 "TB01ID.f"
	co = dasum_(n, &a[j * a_dim1 + 1], &c__1);
#line 245 "TB01ID.f"
	if (withc && *p > 0) {
#line 245 "TB01ID.f"
	    co += dasum_(p, &c__[j * c_dim1 + 1], &c__1);
#line 245 "TB01ID.f"
	}
#line 247 "TB01ID.f"
	snorm = max(snorm,co);
#line 248 "TB01ID.f"
/* L10: */
#line 248 "TB01ID.f"
    }

#line 250 "TB01ID.f"
    if (withb) {

#line 252 "TB01ID.f"
	i__1 = *m;
#line 252 "TB01ID.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 253 "TB01ID.f"
	    d__1 = snorm, d__2 = dasum_(n, &b[j * b_dim1 + 1], &c__1);
#line 253 "TB01ID.f"
	    snorm = max(d__1,d__2);
#line 254 "TB01ID.f"
/* L20: */
#line 254 "TB01ID.f"
	}

#line 256 "TB01ID.f"
    }

#line 258 "TB01ID.f"
    if (snorm == 0.) {
#line 258 "TB01ID.f"
	return 0;
#line 258 "TB01ID.f"
    }

/*     Set some machine parameters and the maximum reduction in the */
/*     1-norm of S if zero rows or columns are encountered. */

#line 264 "TB01ID.f"
    sfmin1 = dlamch_("S", (ftnlen)1) / dlamch_("P", (ftnlen)1);
#line 265 "TB01ID.f"
    sfmax1 = 1. / sfmin1;
#line 266 "TB01ID.f"
    sfmin2 = sfmin1 * 10.;
#line 267 "TB01ID.f"
    sfmax2 = 1. / sfmin2;

#line 269 "TB01ID.f"
    sred = *maxred;
#line 270 "TB01ID.f"
    if (sred <= 0.) {
#line 270 "TB01ID.f"
	sred = 10.;
#line 270 "TB01ID.f"
    }

/* Computing MAX */
#line 272 "TB01ID.f"
    d__1 = snorm / sred;
#line 272 "TB01ID.f"
    maxnrm = max(d__1,sfmin1);

/*     Balance the matrix. */

/*     Iterative loop for norm reduction. */

#line 278 "TB01ID.f"
L30:
#line 279 "TB01ID.f"
    noconv = FALSE_;

#line 281 "TB01ID.f"
    i__1 = *n;
#line 281 "TB01ID.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 282 "TB01ID.f"
	co = 0.;
#line 283 "TB01ID.f"
	ro = 0.;

#line 285 "TB01ID.f"
	i__2 = *n;
#line 285 "TB01ID.f"
	for (j = 1; j <= i__2; ++j) {
#line 286 "TB01ID.f"
	    if (j == i__) {
#line 286 "TB01ID.f"
		goto L40;
#line 286 "TB01ID.f"
	    }
#line 288 "TB01ID.f"
	    co += (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 289 "TB01ID.f"
	    ro += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 290 "TB01ID.f"
L40:
#line 290 "TB01ID.f"
	    ;
#line 290 "TB01ID.f"
	}

#line 292 "TB01ID.f"
	ica = idamax_(n, &a[i__ * a_dim1 + 1], &c__1);
#line 293 "TB01ID.f"
	ca = (d__1 = a[ica + i__ * a_dim1], abs(d__1));
#line 294 "TB01ID.f"
	ira = idamax_(n, &a[i__ + a_dim1], lda);
#line 295 "TB01ID.f"
	ra = (d__1 = a[i__ + ira * a_dim1], abs(d__1));

#line 297 "TB01ID.f"
	if (withc && *p > 0) {
#line 298 "TB01ID.f"
	    co += dasum_(p, &c__[i__ * c_dim1 + 1], &c__1);
#line 299 "TB01ID.f"
	    ica = idamax_(p, &c__[i__ * c_dim1 + 1], &c__1);
/* Computing MAX */
#line 300 "TB01ID.f"
	    d__2 = ca, d__3 = (d__1 = c__[ica + i__ * c_dim1], abs(d__1));
#line 300 "TB01ID.f"
	    ca = max(d__2,d__3);
#line 301 "TB01ID.f"
	}

#line 303 "TB01ID.f"
	if (withb && *m > 0) {
#line 304 "TB01ID.f"
	    ro += dasum_(m, &b[i__ + b_dim1], ldb);
#line 305 "TB01ID.f"
	    ira = idamax_(m, &b[i__ + b_dim1], ldb);
/* Computing MAX */
#line 306 "TB01ID.f"
	    d__2 = ra, d__3 = (d__1 = b[i__ + ira * b_dim1], abs(d__1));
#line 306 "TB01ID.f"
	    ra = max(d__2,d__3);
#line 307 "TB01ID.f"
	}

/*        Special case of zero CO and/or RO. */

#line 311 "TB01ID.f"
	if (co == 0. && ro == 0.) {
#line 311 "TB01ID.f"
	    goto L90;
#line 311 "TB01ID.f"
	}
#line 313 "TB01ID.f"
	if (co == 0.) {
#line 314 "TB01ID.f"
	    if (ro <= maxnrm) {
#line 314 "TB01ID.f"
		goto L90;
#line 314 "TB01ID.f"
	    }
#line 316 "TB01ID.f"
	    co = maxnrm;
#line 317 "TB01ID.f"
	}
#line 318 "TB01ID.f"
	if (ro == 0.) {
#line 319 "TB01ID.f"
	    if (co <= maxnrm) {
#line 319 "TB01ID.f"
		goto L90;
#line 319 "TB01ID.f"
	    }
#line 321 "TB01ID.f"
	    ro = maxnrm;
#line 322 "TB01ID.f"
	}

/*        Guard against zero CO or RO due to underflow. */

#line 326 "TB01ID.f"
	g = ro / 10.;
#line 327 "TB01ID.f"
	f = 1.;
#line 328 "TB01ID.f"
	s = co + ro;
#line 329 "TB01ID.f"
L50:
/* Computing MAX */
#line 330 "TB01ID.f"
	d__1 = max(f,co);
/* Computing MIN */
#line 330 "TB01ID.f"
	d__2 = min(ro,g);
#line 330 "TB01ID.f"
	if (co >= g || max(d__1,ca) >= sfmax2 || min(d__2,ra) <= sfmin2) {
#line 330 "TB01ID.f"
	    goto L60;
#line 330 "TB01ID.f"
	}
#line 332 "TB01ID.f"
	f *= 10.;
#line 333 "TB01ID.f"
	co *= 10.;
#line 334 "TB01ID.f"
	ca *= 10.;
#line 335 "TB01ID.f"
	g /= 10.;
#line 336 "TB01ID.f"
	ro /= 10.;
#line 337 "TB01ID.f"
	ra /= 10.;
#line 338 "TB01ID.f"
	goto L50;

#line 340 "TB01ID.f"
L60:
#line 341 "TB01ID.f"
	g = co / 10.;
#line 342 "TB01ID.f"
L70:
/* Computing MIN */
#line 343 "TB01ID.f"
	d__1 = min(f,co), d__1 = min(d__1,g);
#line 343 "TB01ID.f"
	if (g < ro || max(ro,ra) >= sfmax2 || min(d__1,ca) <= sfmin2) {
#line 343 "TB01ID.f"
	    goto L80;
#line 343 "TB01ID.f"
	}
#line 345 "TB01ID.f"
	f /= 10.;
#line 346 "TB01ID.f"
	co /= 10.;
#line 347 "TB01ID.f"
	ca /= 10.;
#line 348 "TB01ID.f"
	g /= 10.;
#line 349 "TB01ID.f"
	ro *= 10.;
#line 350 "TB01ID.f"
	ra *= 10.;
#line 351 "TB01ID.f"
	goto L70;

/*        Now balance. */

#line 355 "TB01ID.f"
L80:
#line 356 "TB01ID.f"
	if (co + ro >= s * .95) {
#line 356 "TB01ID.f"
	    goto L90;
#line 356 "TB01ID.f"
	}
#line 358 "TB01ID.f"
	if (f < 1. && scale[i__] < 1.) {
#line 359 "TB01ID.f"
	    if (f * scale[i__] <= sfmin1) {
#line 359 "TB01ID.f"
		goto L90;
#line 359 "TB01ID.f"
	    }
#line 361 "TB01ID.f"
	}
#line 362 "TB01ID.f"
	if (f > 1. && scale[i__] > 1.) {
#line 363 "TB01ID.f"
	    if (scale[i__] >= sfmax1 / f) {
#line 363 "TB01ID.f"
		goto L90;
#line 363 "TB01ID.f"
	    }
#line 365 "TB01ID.f"
	}
#line 366 "TB01ID.f"
	g = 1. / f;
#line 367 "TB01ID.f"
	scale[i__] *= f;
#line 368 "TB01ID.f"
	noconv = TRUE_;

#line 370 "TB01ID.f"
	dscal_(n, &g, &a[i__ + a_dim1], lda);
#line 371 "TB01ID.f"
	dscal_(n, &f, &a[i__ * a_dim1 + 1], &c__1);
#line 372 "TB01ID.f"
	if (*m > 0) {
#line 372 "TB01ID.f"
	    dscal_(m, &g, &b[i__ + b_dim1], ldb);
#line 372 "TB01ID.f"
	}
#line 373 "TB01ID.f"
	if (*p > 0) {
#line 373 "TB01ID.f"
	    dscal_(p, &f, &c__[i__ * c_dim1 + 1], &c__1);
#line 373 "TB01ID.f"
	}

#line 375 "TB01ID.f"
L90:
#line 375 "TB01ID.f"
	;
#line 375 "TB01ID.f"
    }

#line 377 "TB01ID.f"
    if (noconv) {
#line 377 "TB01ID.f"
	goto L30;
#line 377 "TB01ID.f"
    }

/*     Set the norm reduction parameter. */

#line 382 "TB01ID.f"
    *maxred = snorm;
#line 383 "TB01ID.f"
    snorm = 0.;

#line 385 "TB01ID.f"
    i__1 = *n;
#line 385 "TB01ID.f"
    for (j = 1; j <= i__1; ++j) {
#line 386 "TB01ID.f"
	co = dasum_(n, &a[j * a_dim1 + 1], &c__1);
#line 387 "TB01ID.f"
	if (withc && *p > 0) {
#line 387 "TB01ID.f"
	    co += dasum_(p, &c__[j * c_dim1 + 1], &c__1);
#line 387 "TB01ID.f"
	}
#line 389 "TB01ID.f"
	snorm = max(snorm,co);
#line 390 "TB01ID.f"
/* L100: */
#line 390 "TB01ID.f"
    }

#line 392 "TB01ID.f"
    if (withb) {

#line 394 "TB01ID.f"
	i__1 = *m;
#line 394 "TB01ID.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 395 "TB01ID.f"
	    d__1 = snorm, d__2 = dasum_(n, &b[j * b_dim1 + 1], &c__1);
#line 395 "TB01ID.f"
	    snorm = max(d__1,d__2);
#line 396 "TB01ID.f"
/* L110: */
#line 396 "TB01ID.f"
	}

#line 398 "TB01ID.f"
    }
#line 399 "TB01ID.f"
    *maxred /= snorm;
#line 400 "TB01ID.f"
    return 0;
/* *** Last line of TB01ID *** */
} /* tb01id_ */

