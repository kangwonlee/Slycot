#line 1 "BB04AD.f"
/* BB04AD.f -- translated by f2c (version 20100827).
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

#line 1 "BB04AD.f"
/* Table of constant values */

static doublereal c_b8 = 0.;
static doublereal c_b9 = 1.;
static integer c__1 = 1;
static doublereal c_b53 = -1.;
static doublereal c_b105 = 2.;

/* Subroutine */ int bb04ad_(char *def, integer *nr, doublereal *dpar, 
	integer *ipar, logical *vec, integer *n, integer *m, doublereal *e, 
	integer *lde, doublereal *a, integer *lda, doublereal *y, integer *
	ldy, doublereal *b, integer *ldb, doublereal *x, integer *ldx, 
	doublereal *u, integer *ldu, char *note, doublereal *dwork, integer *
	ldwork, integer *info, ftnlen def_len, ftnlen note_len)
{
    /* Initialized data */

    static logical vecdef[8] = { TRUE_,TRUE_,FALSE_,TRUE_,TRUE_,FALSE_,FALSE_,
	    FALSE_ };

    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, e_dim1, e_offset, u_dim1, 
	    u_offset, x_dim1, x_offset, y_dim1, y_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double pow_di(doublereal *, integer *), pow_dd(doublereal *, doublereal *)
	    , sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), daxpy_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);
    static doublereal ttemp;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static doublereal twobyn;


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

/*     To generate benchmark examples of (generalized) discrete-time */
/*     Lyapunov equations */

/*        T           T */
/*       A  X  A  -  E  X E  =  Y . */

/*     In some examples, the right hand side has the form */

/*                T */
/*       Y  =  - B  B */

/*     and the solution can be represented as a product of Cholesky */
/*     factors */

/*              T */
/*       X  =  U  U . */

/*     E, A, Y, X, and U are real N-by-N matrices, and B is M-by-N. Note */
/*     that E can be the identity matrix. For some examples, B, X, or U */
/*     are not provided. */

/*     This routine is an implementation of the benchmark library */
/*     DTLEX (Version 1.0) described in [1]. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DEF     CHARACTER*1 */
/*             Specifies the kind of values used as parameters when */
/*             generating parameter-dependent and scalable examples */
/*             (i.e., examples with NR(1) = 2, 3, or 4): */
/*             DEF = 'D' or 'd': Default values are used. */
/*             DEF = 'N' or 'n': Values set in DPAR and IPAR are used. */
/*             This parameter is not referenced if NR(1) = 1. */
/*             Note that the scaling parameter of examples with */
/*             NR(1) = 3 or 4 is considered as a regular parameter in */
/*             this context. */

/*     Input/Output Parameters */

/*     NR      (input) INTEGER array, dimension 2 */
/*             Specifies the index of the desired example according */
/*             to [1]. */
/*             NR(1) defines the group: */
/*                   1 : parameter-free problems of fixed size */
/*                   2 : parameter-dependent problems of fixed size */
/*                   3 : parameter-free problems of scalable size */
/*                   4 : parameter-dependent problems of scalable size */
/*             NR(2) defines the number of the benchmark example */
/*             within a certain group according to [1]. */

/*     DPAR    (input/output) DOUBLE PRECISION array, dimension 2 */
/*             On entry, if DEF = 'N' or 'n' and the desired example */
/*             depends on real parameters, then the array DPAR must */
/*             contain the values for these parameters. */
/*             For an explanation of the parameters see [1]. */
/*             For Example 4.1, DPAR(1) and DPAR(2) define 'r' and 's', */
/*             respectively. */
/*             For Example 4.2, DPAR(1) and DPAR(2) define 'lambda' and */
/*             's', respectively. */
/*             For Examples 4.3 and 4.4, DPAR(1) defines the parameter */
/*             't'. */
/*             On exit, if DEF = 'D' or 'd' and the desired example */
/*             depends on real parameters, then the array DPAR is */
/*             overwritten by the default values given in [1]. */

/*     IPAR    (input/output) INTEGER array of DIMENSION at least 1 */
/*             On entry, if DEF = 'N' or 'n' and the desired example */
/*             depends on integer parameters, then the array IPAR must */
/*             contain the values for these parameters. */
/*             For an explanation of the parameters see [1]. */
/*             For Examples 4.1, 4.2, and 4.3, IPAR(1) defines 'n'. */
/*             For Example 4.4, IPAR(1) defines 'q'. */
/*             On exit, if DEF = 'D' or 'd' and the desired example */
/*             depends on integer parameters, then the array IPAR is */
/*             overwritten by the default values given in [1]. */

/*     VEC     (output) LOGICAL array, dimension 8 */
/*             Flag vector which displays the availability of the output */
/*             data: */
/*             VEC(1) and VEC(2) refer to N and M, respectively, and are */
/*             always .TRUE. */
/*             VEC(3) is .TRUE. iff E is NOT the identity matrix. */
/*             VEC(4) and VEC(5) refer to A and Y, respectively, and are */
/*             always .TRUE. */
/*             VEC(6) is .TRUE. iff B is provided. */
/*             VEC(7) is .TRUE. iff the solution matrix X is provided. */
/*             VEC(8) is .TRUE. iff the Cholesky factor U is provided. */

/*     N       (output) INTEGER */
/*             The actual state dimension, i.e., the order of the */
/*             matrices E and A. */

/*     M       (output) INTEGER */
/*             The number of rows in the matrix B. If B is not provided */
/*             for the desired example, M = 0 is returned. */

/*     E       (output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             The leading N-by-N part of this array contains the */
/*             matrix E. */
/*             NOTE that this array is overwritten (by the identity */
/*             matrix), if VEC(3) = .FALSE. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= N. */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array contains the */
/*             matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= N. */

/*     Y       (output) DOUBLE PRECISION array, dimension (LDY,N) */
/*             The leading N-by-N part of this array contains the */
/*             matrix Y. */

/*     LDY     INTEGER */
/*             The leading dimension of array Y.  LDY >= N. */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             The leading M-by-N part of this array contains the */
/*             matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= M. */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             The leading N-by-N part of this array contains the */
/*             matrix X. */

/*     LDX     INTEGER */
/*             The leading dimension of array X.  LDX >= N. */

/*     U       (output) DOUBLE PRECISION array, dimension (LDU,N) */
/*             The leading N-by-N part of this array contains the */
/*             matrix U. */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= N. */

/*     NOTE    (output) CHARACTER*70 */
/*             String containing short information about the chosen */
/*             example. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             For Examples 4.1 and 4.2., LDWORK >= 2*IPAR(1) is */
/*             required. */
/*             For the other examples, no workspace is needed, i.e., */
/*             LDWORK >= 1. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; in particular, INFO = -3 or -4 indicates */
/*                   that at least one of the parameters in DPAR or */
/*                   IPAR, respectively, has an illegal value. */

/*     REFERENCES */

/*     [1]  D. Kressner, V. Mehrmann, and T. Penzl. */
/*          DTLEX - a Collection of Benchmark Examples for Discrete- */
/*          Time Lyapunov Equations. */
/*          SLICOT Working Note 1999-7, 1999. */

/*     NUMERICAL ASPECTS */

/*     None */

/*     CONTRIBUTOR */

/*     D. Kressner, V. Mehrmann, and T. Penzl (TU Chemnitz) */

/*     For questions concerning the collection or for the submission of */
/*     test examples, please contact Volker Mehrmann */
/*     (Email: volker.mehrmann@mathematik.tu-chemnitz.de). */

/*     REVISIONS */

/*     June 1999, V. Sima. */

/*     KEYWORDS */

/*     discrete-time Lyapunov equations */

/*     ******************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     . BLAS . */
/*     . LAPACK . */
/*     .. External Subroutines .. */
/*     . BLAS . */
/*     . LAPACK . */
/*     .. Intrinsic Functions .. */
/*     .. Data Statements .. */
/*     . default values for availabilities . */
#line 254 "BB04AD.f"
    /* Parameter adjustments */
#line 254 "BB04AD.f"
    --nr;
#line 254 "BB04AD.f"
    --dpar;
#line 254 "BB04AD.f"
    --ipar;
#line 254 "BB04AD.f"
    --vec;
#line 254 "BB04AD.f"
    e_dim1 = *lde;
#line 254 "BB04AD.f"
    e_offset = 1 + e_dim1;
#line 254 "BB04AD.f"
    e -= e_offset;
#line 254 "BB04AD.f"
    a_dim1 = *lda;
#line 254 "BB04AD.f"
    a_offset = 1 + a_dim1;
#line 254 "BB04AD.f"
    a -= a_offset;
#line 254 "BB04AD.f"
    y_dim1 = *ldy;
#line 254 "BB04AD.f"
    y_offset = 1 + y_dim1;
#line 254 "BB04AD.f"
    y -= y_offset;
#line 254 "BB04AD.f"
    b_dim1 = *ldb;
#line 254 "BB04AD.f"
    b_offset = 1 + b_dim1;
#line 254 "BB04AD.f"
    b -= b_offset;
#line 254 "BB04AD.f"
    x_dim1 = *ldx;
#line 254 "BB04AD.f"
    x_offset = 1 + x_dim1;
#line 254 "BB04AD.f"
    x -= x_offset;
#line 254 "BB04AD.f"
    u_dim1 = *ldu;
#line 254 "BB04AD.f"
    u_offset = 1 + u_dim1;
#line 254 "BB04AD.f"
    u -= u_offset;
#line 254 "BB04AD.f"
    --dwork;
#line 254 "BB04AD.f"

#line 254 "BB04AD.f"
    /* Function Body */

/*     .. Executable Statements .. */

#line 259 "BB04AD.f"
    *info = 0;
#line 260 "BB04AD.f"
    for (i__ = 1; i__ <= 8; ++i__) {
#line 261 "BB04AD.f"
	vec[i__] = vecdef[i__ - 1];
#line 262 "BB04AD.f"
/* L10: */
#line 262 "BB04AD.f"
    }

#line 264 "BB04AD.f"
    if (nr[1] == 4) {
#line 265 "BB04AD.f"
	if (! (lsame_(def, "D", (ftnlen)1, (ftnlen)1) || lsame_(def, "N", (
		ftnlen)1, (ftnlen)1))) {
#line 266 "BB04AD.f"
	    *info = -1;
#line 267 "BB04AD.f"
	    return 0;
#line 268 "BB04AD.f"
	}

#line 270 "BB04AD.f"
	if (nr[2] == 1) {
#line 271 "BB04AD.f"
	    s_copy(note, "DTLEX: Example 4.1", (ftnlen)70, (ftnlen)18);
#line 272 "BB04AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 273 "BB04AD.f"
		ipar[1] = 10;
#line 274 "BB04AD.f"
		dpar[1] = 1.5;
#line 275 "BB04AD.f"
		dpar[2] = 1.5;
#line 276 "BB04AD.f"
	    }
#line 277 "BB04AD.f"
	    if (dpar[1] <= 1. || dpar[2] <= 1.) {
#line 277 "BB04AD.f"
		*info = -3;
#line 277 "BB04AD.f"
	    }
#line 278 "BB04AD.f"
	    if (ipar[1] < 2) {
#line 278 "BB04AD.f"
		*info = -4;
#line 278 "BB04AD.f"
	    }
#line 279 "BB04AD.f"
	    *n = ipar[1];
#line 280 "BB04AD.f"
	    *m = 1;
#line 281 "BB04AD.f"
	    if (*lde < *n) {
#line 281 "BB04AD.f"
		*info = -9;
#line 281 "BB04AD.f"
	    }
#line 282 "BB04AD.f"
	    if (*lda < *n) {
#line 282 "BB04AD.f"
		*info = -11;
#line 282 "BB04AD.f"
	    }
#line 283 "BB04AD.f"
	    if (*ldy < *n) {
#line 283 "BB04AD.f"
		*info = -13;
#line 283 "BB04AD.f"
	    }
#line 284 "BB04AD.f"
	    if (*ldb < *m) {
#line 284 "BB04AD.f"
		*info = -15;
#line 284 "BB04AD.f"
	    }
#line 285 "BB04AD.f"
	    if (*ldx < *n) {
#line 285 "BB04AD.f"
		*info = -17;
#line 285 "BB04AD.f"
	    }
#line 286 "BB04AD.f"
	    if (*ldwork < *n << 1) {
#line 286 "BB04AD.f"
		*info = -22;
#line 286 "BB04AD.f"
	    }
#line 287 "BB04AD.f"
	    if (*info != 0) {
#line 287 "BB04AD.f"
		return 0;
#line 287 "BB04AD.f"
	    }

#line 289 "BB04AD.f"
	    vec[6] = TRUE_;
#line 290 "BB04AD.f"
	    vec[7] = TRUE_;
#line 291 "BB04AD.f"
	    twobyn = 2. / (doublereal) (*n);
#line 292 "BB04AD.f"
	    dlaset_("A", n, n, &c_b8, &c_b9, &e[e_offset], lde, (ftnlen)1);
#line 293 "BB04AD.f"
	    dlaset_("A", n, n, &c_b8, &c_b8, &a[a_offset], lda, (ftnlen)1);
#line 294 "BB04AD.f"
	    dlaset_("A", n, n, &c_b8, &c_b8, &y[y_offset], ldy, (ftnlen)1);
#line 295 "BB04AD.f"
	    d__1 = -twobyn;
#line 295 "BB04AD.f"
	    d__2 = 1. - twobyn;
#line 295 "BB04AD.f"
	    dlaset_("A", m, n, &d__1, &d__2, &b[b_offset], ldb, (ftnlen)1);
#line 296 "BB04AD.f"
	    dlaset_("A", n, n, &c_b8, &c_b8, &x[x_offset], ldx, (ftnlen)1);
#line 297 "BB04AD.f"
	    i__1 = *n;
#line 297 "BB04AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 298 "BB04AD.f"
		i__2 = i__ - 1;
#line 298 "BB04AD.f"
		temp = pow_di(&dpar[1], &i__2);
#line 299 "BB04AD.f"
		a[i__ + i__ * a_dim1] = (temp - 1.) / (temp + 1.);
#line 300 "BB04AD.f"
		dwork[i__] = 1.;
#line 301 "BB04AD.f"
/* L20: */
#line 301 "BB04AD.f"
	    }
/*         H1 * A */
#line 303 "BB04AD.f"
	    dgemv_("T", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 304 "BB04AD.f"
	    d__1 = -twobyn;
#line 304 "BB04AD.f"
	    dger_(n, n, &d__1, &dwork[1], &c__1, &dwork[*n + 1], &c__1, &a[
		    a_offset], lda);
/*         A * H1 */
#line 306 "BB04AD.f"
	    dgemv_("N", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 307 "BB04AD.f"
	    d__1 = -twobyn;
#line 307 "BB04AD.f"
	    dger_(n, n, &d__1, &dwork[*n + 1], &c__1, &dwork[1], &c__1, &a[
		    a_offset], lda);
/*         S A INV(S), B INV(S) */
#line 309 "BB04AD.f"
	    i__1 = *n;
#line 309 "BB04AD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 310 "BB04AD.f"
		i__2 = j - 1;
#line 310 "BB04AD.f"
		b[j * b_dim1 + 1] /= pow_di(&dpar[2], &i__2);
#line 311 "BB04AD.f"
		i__2 = *n;
#line 311 "BB04AD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 312 "BB04AD.f"
		    i__3 = i__ - j;
#line 312 "BB04AD.f"
		    a[i__ + j * a_dim1] *= pow_di(&dpar[2], &i__3);
#line 313 "BB04AD.f"
/* L30: */
#line 313 "BB04AD.f"
		}
#line 314 "BB04AD.f"
		dwork[j] = 1. - j % 2 * 2.;
#line 315 "BB04AD.f"
/* L40: */
#line 315 "BB04AD.f"
	    }
/*         H2 * A */
#line 317 "BB04AD.f"
	    dgemv_("T", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 318 "BB04AD.f"
	    d__1 = -twobyn;
#line 318 "BB04AD.f"
	    dger_(n, n, &d__1, &dwork[1], &c__1, &dwork[*n + 1], &c__1, &a[
		    a_offset], lda);
/*         A * H2 */
#line 320 "BB04AD.f"
	    dgemv_("N", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 321 "BB04AD.f"
	    d__1 = -twobyn;
#line 321 "BB04AD.f"
	    dger_(n, n, &d__1, &dwork[*n + 1], &c__1, &dwork[1], &c__1, &a[
		    a_offset], lda);
/*         B * H2 */
#line 323 "BB04AD.f"
	    d__1 = -twobyn * ddot_(n, &b[b_offset], ldb, &dwork[1], &c__1);
#line 323 "BB04AD.f"
	    daxpy_(n, &d__1, &dwork[1], &c__1, &b[b_offset], ldb);
/*         Y = -B' * B */
#line 326 "BB04AD.f"
	    dger_(n, n, &c_b53, &b[b_offset], ldb, &b[b_offset], ldb, &y[
		    y_offset], ldy);
/*         X = -Y */
#line 328 "BB04AD.f"
	    i__1 = *n;
#line 328 "BB04AD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 329 "BB04AD.f"
		daxpy_(n, &c_b53, &y[j * y_dim1 + 1], &c__1, &x[j * x_dim1 + 
			1], &c__1);
#line 330 "BB04AD.f"
/* L50: */
#line 330 "BB04AD.f"
	    }

#line 332 "BB04AD.f"
	} else if (nr[2] == 2) {
#line 333 "BB04AD.f"
	    s_copy(note, "DTLEX: Example 4.2", (ftnlen)70, (ftnlen)18);
#line 334 "BB04AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 335 "BB04AD.f"
		ipar[1] = 10;
#line 336 "BB04AD.f"
		dpar[1] = -.5;
#line 337 "BB04AD.f"
		dpar[2] = 1.5;
#line 338 "BB04AD.f"
	    }
#line 339 "BB04AD.f"
	    if (dpar[1] <= -1. || dpar[1] >= 1. || dpar[2] <= 1.) {
#line 339 "BB04AD.f"
		*info = -3;
#line 339 "BB04AD.f"
	    }
#line 341 "BB04AD.f"
	    if (ipar[1] < 2) {
#line 341 "BB04AD.f"
		*info = -4;
#line 341 "BB04AD.f"
	    }
#line 342 "BB04AD.f"
	    *n = ipar[1];
#line 343 "BB04AD.f"
	    *m = 1;
#line 344 "BB04AD.f"
	    if (*lde < *n) {
#line 344 "BB04AD.f"
		*info = -9;
#line 344 "BB04AD.f"
	    }
#line 345 "BB04AD.f"
	    if (*lda < *n) {
#line 345 "BB04AD.f"
		*info = -11;
#line 345 "BB04AD.f"
	    }
#line 346 "BB04AD.f"
	    if (*ldy < *n) {
#line 346 "BB04AD.f"
		*info = -13;
#line 346 "BB04AD.f"
	    }
#line 347 "BB04AD.f"
	    if (*ldb < *m) {
#line 347 "BB04AD.f"
		*info = -15;
#line 347 "BB04AD.f"
	    }
#line 348 "BB04AD.f"
	    if (*ldwork < *n << 1) {
#line 348 "BB04AD.f"
		*info = -22;
#line 348 "BB04AD.f"
	    }
#line 349 "BB04AD.f"
	    if (*info != 0) {
#line 349 "BB04AD.f"
		return 0;
#line 349 "BB04AD.f"
	    }

#line 351 "BB04AD.f"
	    vec[6] = TRUE_;
#line 352 "BB04AD.f"
	    twobyn = 2. / (doublereal) (*n);
#line 353 "BB04AD.f"
	    dlaset_("A", n, n, &c_b8, &c_b9, &e[e_offset], lde, (ftnlen)1);
#line 354 "BB04AD.f"
	    dlaset_("A", n, n, &c_b8, &dpar[1], &a[a_offset], lda, (ftnlen)1);
#line 355 "BB04AD.f"
	    dlaset_("A", n, n, &c_b8, &c_b8, &y[y_offset], ldy, (ftnlen)1);
#line 356 "BB04AD.f"
	    d__1 = -twobyn;
#line 356 "BB04AD.f"
	    d__2 = 1. - twobyn;
#line 356 "BB04AD.f"
	    dlaset_("A", m, n, &d__1, &d__2, &b[b_offset], ldb, (ftnlen)1);
#line 357 "BB04AD.f"
	    i__1 = *n - 1;
#line 357 "BB04AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 358 "BB04AD.f"
		dwork[i__] = 1.;
#line 359 "BB04AD.f"
		a[i__ + (i__ + 1) * a_dim1] = 1.;
#line 360 "BB04AD.f"
/* L60: */
#line 360 "BB04AD.f"
	    }
#line 361 "BB04AD.f"
	    dwork[*n] = 1.;
/*         H1 * A */
#line 363 "BB04AD.f"
	    dgemv_("T", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 364 "BB04AD.f"
	    d__1 = -twobyn;
#line 364 "BB04AD.f"
	    dger_(n, n, &d__1, &dwork[1], &c__1, &dwork[*n + 1], &c__1, &a[
		    a_offset], lda);
/*         A * H1 */
#line 366 "BB04AD.f"
	    dgemv_("N", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 367 "BB04AD.f"
	    d__1 = -twobyn;
#line 367 "BB04AD.f"
	    dger_(n, n, &d__1, &dwork[*n + 1], &c__1, &dwork[1], &c__1, &a[
		    a_offset], lda);
/*         S A INV(S), B INV(S) */
#line 369 "BB04AD.f"
	    i__1 = *n;
#line 369 "BB04AD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 370 "BB04AD.f"
		i__2 = j - 1;
#line 370 "BB04AD.f"
		b[j * b_dim1 + 1] /= pow_di(&dpar[2], &i__2);
#line 371 "BB04AD.f"
		i__2 = *n;
#line 371 "BB04AD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 372 "BB04AD.f"
		    i__3 = i__ - j;
#line 372 "BB04AD.f"
		    a[i__ + j * a_dim1] *= pow_di(&dpar[2], &i__3);
#line 373 "BB04AD.f"
/* L70: */
#line 373 "BB04AD.f"
		}
#line 374 "BB04AD.f"
		dwork[j] = 1. - j % 2 * 2.;
#line 375 "BB04AD.f"
/* L80: */
#line 375 "BB04AD.f"
	    }
/*         H2 * A */
#line 377 "BB04AD.f"
	    dgemv_("T", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 378 "BB04AD.f"
	    d__1 = -twobyn;
#line 378 "BB04AD.f"
	    dger_(n, n, &d__1, &dwork[1], &c__1, &dwork[*n + 1], &c__1, &a[
		    a_offset], lda);
/*         A * H2 */
#line 380 "BB04AD.f"
	    dgemv_("N", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 381 "BB04AD.f"
	    d__1 = -twobyn;
#line 381 "BB04AD.f"
	    dger_(n, n, &d__1, &dwork[*n + 1], &c__1, &dwork[1], &c__1, &a[
		    a_offset], lda);
/*         B * H2 */
#line 383 "BB04AD.f"
	    d__1 = -twobyn * ddot_(n, &b[b_offset], ldb, &dwork[1], &c__1);
#line 383 "BB04AD.f"
	    daxpy_(n, &d__1, &dwork[1], &c__1, &b[b_offset], ldb);
/*         Y = -B' * B */
#line 386 "BB04AD.f"
	    dger_(n, n, &c_b53, &b[b_offset], ldb, &b[b_offset], ldb, &y[
		    y_offset], ldy);

#line 388 "BB04AD.f"
	} else if (nr[2] == 3) {
#line 389 "BB04AD.f"
	    s_copy(note, "DTLEX: Example 4.3", (ftnlen)70, (ftnlen)18);
#line 390 "BB04AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 391 "BB04AD.f"
		ipar[1] = 10;
#line 392 "BB04AD.f"
		dpar[1] = 10.;
#line 393 "BB04AD.f"
	    }
#line 394 "BB04AD.f"
	    if (dpar[1] < 0.) {
#line 394 "BB04AD.f"
		*info = -3;
#line 394 "BB04AD.f"
	    }
#line 395 "BB04AD.f"
	    if (ipar[1] < 2) {
#line 395 "BB04AD.f"
		*info = -4;
#line 395 "BB04AD.f"
	    }
#line 396 "BB04AD.f"
	    *n = ipar[1];
#line 397 "BB04AD.f"
	    *m = 0;
#line 398 "BB04AD.f"
	    if (*lde < *n) {
#line 398 "BB04AD.f"
		*info = -9;
#line 398 "BB04AD.f"
	    }
#line 399 "BB04AD.f"
	    if (*lda < *n) {
#line 399 "BB04AD.f"
		*info = -11;
#line 399 "BB04AD.f"
	    }
#line 400 "BB04AD.f"
	    if (*ldy < *n) {
#line 400 "BB04AD.f"
		*info = -13;
#line 400 "BB04AD.f"
	    }
#line 401 "BB04AD.f"
	    if (*ldx < *n) {
#line 401 "BB04AD.f"
		*info = -17;
#line 401 "BB04AD.f"
	    }
#line 402 "BB04AD.f"
	    if (*info != 0) {
#line 402 "BB04AD.f"
		return 0;
#line 402 "BB04AD.f"
	    }

#line 404 "BB04AD.f"
	    vec[3] = TRUE_;
#line 405 "BB04AD.f"
	    vec[7] = TRUE_;
#line 406 "BB04AD.f"
	    d__1 = -dpar[1];
#line 406 "BB04AD.f"
	    temp = pow_dd(&c_b105, &d__1);
#line 407 "BB04AD.f"
	    dlaset_("U", n, n, &c_b8, &c_b8, &e[e_offset], lde, (ftnlen)1);
#line 408 "BB04AD.f"
	    dlaset_("L", n, n, &temp, &c_b9, &e[e_offset], lde, (ftnlen)1);
#line 409 "BB04AD.f"
	    dlaset_("L", n, n, &c_b8, &c_b8, &a[a_offset], lda, (ftnlen)1);
#line 410 "BB04AD.f"
	    dlaset_("U", n, n, &c_b9, &c_b8, &a[a_offset], lda, (ftnlen)1);
#line 411 "BB04AD.f"
	    dlaset_("A", n, n, &c_b9, &c_b9, &x[x_offset], ldx, (ftnlen)1);
#line 412 "BB04AD.f"
	    i__1 = *n;
#line 412 "BB04AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 413 "BB04AD.f"
		a[i__ + i__ * a_dim1] = (doublereal) i__ + temp;
#line 414 "BB04AD.f"
/* L90: */
#line 414 "BB04AD.f"
	    }
#line 415 "BB04AD.f"
	    i__1 = *n;
#line 415 "BB04AD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 416 "BB04AD.f"
		i__2 = *n;
#line 416 "BB04AD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 417 "BB04AD.f"
		    y[i__ + j * y_dim1] = temp * temp * (doublereal) (1 - (*n 
			    - i__) * (*n - j)) + temp * (doublereal) ((i__ + 
			    j) * 3 - (*n + 1 << 1)) + (doublereal) (i__ * j) *
			     4. - (doublereal) (i__ + j) * 2.;
#line 420 "BB04AD.f"
/* L100: */
#line 420 "BB04AD.f"
		}
#line 421 "BB04AD.f"
/* L110: */
#line 421 "BB04AD.f"
	    }

#line 423 "BB04AD.f"
	} else if (nr[2] == 4) {
#line 424 "BB04AD.f"
	    s_copy(note, "DTLEX: Example 4.4", (ftnlen)70, (ftnlen)18);
#line 425 "BB04AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 426 "BB04AD.f"
		ipar[1] = 10;
#line 427 "BB04AD.f"
		dpar[1] = 1.5;
#line 428 "BB04AD.f"
	    }
#line 429 "BB04AD.f"
	    if (dpar[1] < 1.) {
#line 429 "BB04AD.f"
		*info = -3;
#line 429 "BB04AD.f"
	    }
#line 430 "BB04AD.f"
	    if (ipar[1] < 1) {
#line 430 "BB04AD.f"
		*info = -4;
#line 430 "BB04AD.f"
	    }
#line 431 "BB04AD.f"
	    *n = ipar[1] * 3;
#line 432 "BB04AD.f"
	    *m = 1;
#line 433 "BB04AD.f"
	    if (*lde < *n) {
#line 433 "BB04AD.f"
		*info = -9;
#line 433 "BB04AD.f"
	    }
#line 434 "BB04AD.f"
	    if (*lda < *n) {
#line 434 "BB04AD.f"
		*info = -11;
#line 434 "BB04AD.f"
	    }
#line 435 "BB04AD.f"
	    if (*ldy < *n) {
#line 435 "BB04AD.f"
		*info = -13;
#line 435 "BB04AD.f"
	    }
#line 436 "BB04AD.f"
	    if (*ldb < *m) {
#line 436 "BB04AD.f"
		*info = -15;
#line 436 "BB04AD.f"
	    }
#line 437 "BB04AD.f"
	    if (*info != 0) {
#line 437 "BB04AD.f"
		return 0;
#line 437 "BB04AD.f"
	    }

#line 439 "BB04AD.f"
	    vec[3] = TRUE_;
#line 440 "BB04AD.f"
	    vec[6] = TRUE_;
#line 441 "BB04AD.f"
	    dlaset_("A", n, n, &c_b8, &c_b8, &e[e_offset], lde, (ftnlen)1);
#line 442 "BB04AD.f"
	    dlaset_("A", n, n, &c_b8, &c_b8, &a[a_offset], lda, (ftnlen)1);
#line 443 "BB04AD.f"
	    i__1 = ipar[1];
#line 443 "BB04AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 444 "BB04AD.f"
		ttemp = 1. - 1. / pow_di(&dpar[1], &i__);
#line 445 "BB04AD.f"
		temp = -ttemp / sqrt(2.);
#line 446 "BB04AD.f"
		i__2 = i__ - 1;
#line 446 "BB04AD.f"
		for (j = 1; j <= i__2; ++j) {
#line 447 "BB04AD.f"
		    for (k = 0; k <= 2; ++k) {
#line 448 "BB04AD.f"
			a[*n - i__ * 3 + 3 + (j * 3 - k) * a_dim1] = ttemp;
#line 449 "BB04AD.f"
			a[*n - i__ * 3 + 2 + (j * 3 - k) * a_dim1] = temp * 
				2.;
#line 450 "BB04AD.f"
/* L120: */
#line 450 "BB04AD.f"
		    }
#line 451 "BB04AD.f"
/* L130: */
#line 451 "BB04AD.f"
		}
#line 452 "BB04AD.f"
		a[*n - i__ * 3 + 3 + (i__ * 3 - 2) * a_dim1] = ttemp;
#line 453 "BB04AD.f"
		a[*n - i__ * 3 + 2 + (i__ * 3 - 2) * a_dim1] = temp * 2.;
#line 454 "BB04AD.f"
		a[*n - i__ * 3 + 2 + (i__ * 3 - 1) * a_dim1] = temp * 2.;
#line 455 "BB04AD.f"
		a[*n - i__ * 3 + 2 + i__ * 3 * a_dim1] = temp;
#line 456 "BB04AD.f"
		a[*n - i__ * 3 + 1 + i__ * 3 * a_dim1] = temp;
#line 457 "BB04AD.f"
/* L140: */
#line 457 "BB04AD.f"
	    }
#line 458 "BB04AD.f"
	    i__1 = *n;
#line 458 "BB04AD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 459 "BB04AD.f"
		if (j > 1) {
#line 459 "BB04AD.f"
		    daxpy_(n, &c_b9, &a[j - 1 + a_dim1], lda, &a[j + a_dim1], 
			    lda);
#line 459 "BB04AD.f"
		}
#line 460 "BB04AD.f"
		b[j * b_dim1 + 1] = (doublereal) j;
#line 461 "BB04AD.f"
		i__2 = *n;
#line 461 "BB04AD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 462 "BB04AD.f"
		    e[i__ + (*n - j + 1) * e_dim1] = (doublereal) min(i__,j);
#line 463 "BB04AD.f"
		    y[i__ + j * y_dim1] = -((doublereal) (i__ * j));
#line 464 "BB04AD.f"
/* L150: */
#line 464 "BB04AD.f"
		}
#line 465 "BB04AD.f"
/* L160: */
#line 465 "BB04AD.f"
	    }

#line 467 "BB04AD.f"
	} else {
#line 468 "BB04AD.f"
	    *info = -2;
#line 469 "BB04AD.f"
	}
#line 470 "BB04AD.f"
    } else {
#line 471 "BB04AD.f"
	*info = -2;
#line 472 "BB04AD.f"
    }

#line 474 "BB04AD.f"
    return 0;
/* *** Last Line of BB04AD *** */
} /* bb04ad_ */

