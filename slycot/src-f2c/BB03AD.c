#line 1 "BB03AD.f"
/* BB03AD.f -- translated by f2c (version 20100827).
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

#line 1 "BB03AD.f"
/* Table of constant values */

static doublereal c_b8 = 0.;
static doublereal c_b9 = 1.;
static integer c__1 = 1;
static doublereal c_b84 = -1.;
static doublereal c_b132 = 2.;

/* Subroutine */ int bb03ad_(char *def, integer *nr, doublereal *dpar, 
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
	    ;

    /* Local variables */
    static integer i__, j, k;
    static doublereal ttm1, ttp1;
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
	    doublereal *, doublereal *, integer *, doublereal *, integer *), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
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

/*     To generate benchmark examples of (generalized) continuous-time */
/*     Lyapunov equations */

/*        T           T */
/*       A  X  E  +  E  X A  =  Y . */

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
/*     CTLEX (Version 1.0) described in [1]. */

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
/*          CTLEX - a Collection of Benchmark Examples for Continuous- */
/*          Time Lyapunov Equations. */
/*          SLICOT Working Note 1999-6, 1999. */

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

/*     continuous-time Lyapunov equations */

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
#line 254 "BB03AD.f"
    /* Parameter adjustments */
#line 254 "BB03AD.f"
    --nr;
#line 254 "BB03AD.f"
    --dpar;
#line 254 "BB03AD.f"
    --ipar;
#line 254 "BB03AD.f"
    --vec;
#line 254 "BB03AD.f"
    e_dim1 = *lde;
#line 254 "BB03AD.f"
    e_offset = 1 + e_dim1;
#line 254 "BB03AD.f"
    e -= e_offset;
#line 254 "BB03AD.f"
    a_dim1 = *lda;
#line 254 "BB03AD.f"
    a_offset = 1 + a_dim1;
#line 254 "BB03AD.f"
    a -= a_offset;
#line 254 "BB03AD.f"
    y_dim1 = *ldy;
#line 254 "BB03AD.f"
    y_offset = 1 + y_dim1;
#line 254 "BB03AD.f"
    y -= y_offset;
#line 254 "BB03AD.f"
    b_dim1 = *ldb;
#line 254 "BB03AD.f"
    b_offset = 1 + b_dim1;
#line 254 "BB03AD.f"
    b -= b_offset;
#line 254 "BB03AD.f"
    x_dim1 = *ldx;
#line 254 "BB03AD.f"
    x_offset = 1 + x_dim1;
#line 254 "BB03AD.f"
    x -= x_offset;
#line 254 "BB03AD.f"
    u_dim1 = *ldu;
#line 254 "BB03AD.f"
    u_offset = 1 + u_dim1;
#line 254 "BB03AD.f"
    u -= u_offset;
#line 254 "BB03AD.f"
    --dwork;
#line 254 "BB03AD.f"

#line 254 "BB03AD.f"
    /* Function Body */

/*     .. Executable Statements .. */

#line 259 "BB03AD.f"
    *info = 0;
#line 260 "BB03AD.f"
    for (i__ = 1; i__ <= 8; ++i__) {
#line 261 "BB03AD.f"
	vec[i__] = vecdef[i__ - 1];
#line 262 "BB03AD.f"
/* L10: */
#line 262 "BB03AD.f"
    }

#line 264 "BB03AD.f"
    if (nr[1] == 4) {
#line 265 "BB03AD.f"
	if (! (lsame_(def, "D", (ftnlen)1, (ftnlen)1) || lsame_(def, "N", (
		ftnlen)1, (ftnlen)1))) {
#line 266 "BB03AD.f"
	    *info = -1;
#line 267 "BB03AD.f"
	    return 0;
#line 268 "BB03AD.f"
	}

#line 270 "BB03AD.f"
	if (nr[2] == 1) {
#line 271 "BB03AD.f"
	    s_copy(note, "CTLEX: Example 4.1", (ftnlen)70, (ftnlen)18);
#line 272 "BB03AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 273 "BB03AD.f"
		ipar[1] = 10;
#line 274 "BB03AD.f"
		dpar[1] = 1.5;
#line 275 "BB03AD.f"
		dpar[2] = 1.5;
#line 276 "BB03AD.f"
	    }
#line 277 "BB03AD.f"
	    if (dpar[1] <= 1. || dpar[2] <= 1.) {
#line 277 "BB03AD.f"
		*info = -3;
#line 277 "BB03AD.f"
	    }
#line 278 "BB03AD.f"
	    if (ipar[1] < 2) {
#line 278 "BB03AD.f"
		*info = -4;
#line 278 "BB03AD.f"
	    }
#line 279 "BB03AD.f"
	    *n = ipar[1];
#line 280 "BB03AD.f"
	    *m = 1;
#line 281 "BB03AD.f"
	    if (*lde < *n) {
#line 281 "BB03AD.f"
		*info = -9;
#line 281 "BB03AD.f"
	    }
#line 282 "BB03AD.f"
	    if (*lda < *n) {
#line 282 "BB03AD.f"
		*info = -11;
#line 282 "BB03AD.f"
	    }
#line 283 "BB03AD.f"
	    if (*ldy < *n) {
#line 283 "BB03AD.f"
		*info = -13;
#line 283 "BB03AD.f"
	    }
#line 284 "BB03AD.f"
	    if (*ldb < *m) {
#line 284 "BB03AD.f"
		*info = -15;
#line 284 "BB03AD.f"
	    }
#line 285 "BB03AD.f"
	    if (*ldx < *n) {
#line 285 "BB03AD.f"
		*info = -17;
#line 285 "BB03AD.f"
	    }
#line 286 "BB03AD.f"
	    if (*ldwork < *n << 1) {
#line 286 "BB03AD.f"
		*info = -22;
#line 286 "BB03AD.f"
	    }
#line 287 "BB03AD.f"
	    if (*info != 0) {
#line 287 "BB03AD.f"
		return 0;
#line 287 "BB03AD.f"
	    }

#line 289 "BB03AD.f"
	    vec[6] = TRUE_;
#line 290 "BB03AD.f"
	    vec[7] = TRUE_;
#line 291 "BB03AD.f"
	    twobyn = 2. / (doublereal) (*n);
#line 292 "BB03AD.f"
	    dlaset_("A", n, n, &c_b8, &c_b9, &e[e_offset], lde, (ftnlen)1);
#line 293 "BB03AD.f"
	    dlaset_("A", n, n, &c_b8, &c_b8, &a[a_offset], lda, (ftnlen)1);
#line 294 "BB03AD.f"
	    dlaset_("A", n, n, &c_b8, &c_b8, &y[y_offset], ldy, (ftnlen)1);
#line 295 "BB03AD.f"
	    dlaset_("A", m, n, &c_b8, &c_b8, &b[b_offset], ldb, (ftnlen)1);
#line 296 "BB03AD.f"
	    dlaset_("A", n, n, &c_b8, &c_b8, &x[x_offset], ldx, (ftnlen)1);
#line 297 "BB03AD.f"
	    i__1 = *n;
#line 297 "BB03AD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 298 "BB03AD.f"
		i__2 = j - 1;
#line 298 "BB03AD.f"
		temp = pow_di(&dpar[1], &i__2);
#line 299 "BB03AD.f"
		a[j + j * a_dim1] = -temp;
#line 300 "BB03AD.f"
		dwork[j] = 1.;
#line 301 "BB03AD.f"
		i__2 = *n;
#line 301 "BB03AD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 302 "BB03AD.f"
		    i__3 = i__ - 1;
#line 302 "BB03AD.f"
		    x[i__ + j * x_dim1] = (doublereal) (i__ * j) / (temp + 
			    pow_di(&dpar[1], &i__3));
#line 303 "BB03AD.f"
/* L20: */
#line 303 "BB03AD.f"
		}
#line 304 "BB03AD.f"
/* L30: */
#line 304 "BB03AD.f"
	    }
/*         H1 * A */
#line 306 "BB03AD.f"
	    dgemv_("T", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 307 "BB03AD.f"
	    d__1 = -twobyn;
#line 307 "BB03AD.f"
	    dger_(n, n, &d__1, &dwork[1], &c__1, &dwork[*n + 1], &c__1, &a[
		    a_offset], lda);
/*         A * H1 */
#line 309 "BB03AD.f"
	    dgemv_("N", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 310 "BB03AD.f"
	    d__1 = -twobyn;
#line 310 "BB03AD.f"
	    dger_(n, n, &d__1, &dwork[*n + 1], &c__1, &dwork[1], &c__1, &a[
		    a_offset], lda);
/*         H1 * X */
#line 312 "BB03AD.f"
	    dgemv_("T", n, n, &c_b9, &x[x_offset], ldx, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 313 "BB03AD.f"
	    d__1 = -twobyn;
#line 313 "BB03AD.f"
	    dger_(n, n, &d__1, &dwork[1], &c__1, &dwork[*n + 1], &c__1, &x[
		    x_offset], ldx);
/*         X * H1 */
#line 315 "BB03AD.f"
	    dgemv_("N", n, n, &c_b9, &x[x_offset], ldx, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 316 "BB03AD.f"
	    d__1 = -twobyn;
#line 316 "BB03AD.f"
	    dger_(n, n, &d__1, &dwork[*n + 1], &c__1, &dwork[1], &c__1, &x[
		    x_offset], ldx);
/*         S A INV(S), INV(S) X INV(S), B INV(S) */
#line 318 "BB03AD.f"
	    i__1 = *n;
#line 318 "BB03AD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 319 "BB03AD.f"
		i__2 = j - 1;
#line 319 "BB03AD.f"
		b[j * b_dim1 + 1] = (doublereal) (j - *n - 1) / pow_di(&dpar[
			2], &i__2);
#line 320 "BB03AD.f"
		i__2 = *n;
#line 320 "BB03AD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 321 "BB03AD.f"
		    i__3 = i__ + j - 2;
#line 321 "BB03AD.f"
		    x[i__ + j * x_dim1] /= pow_di(&dpar[2], &i__3);
#line 322 "BB03AD.f"
		    i__3 = i__ - j;
#line 322 "BB03AD.f"
		    a[i__ + j * a_dim1] *= pow_di(&dpar[2], &i__3);
#line 323 "BB03AD.f"
/* L40: */
#line 323 "BB03AD.f"
		}
#line 324 "BB03AD.f"
		dwork[j] = 1. - j % 2 * 2.;
#line 325 "BB03AD.f"
/* L50: */
#line 325 "BB03AD.f"
	    }
/*         H2 * A */
#line 327 "BB03AD.f"
	    dgemv_("T", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 328 "BB03AD.f"
	    d__1 = -twobyn;
#line 328 "BB03AD.f"
	    dger_(n, n, &d__1, &dwork[1], &c__1, &dwork[*n + 1], &c__1, &a[
		    a_offset], lda);
/*         A * H2 */
#line 330 "BB03AD.f"
	    dgemv_("N", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 331 "BB03AD.f"
	    d__1 = -twobyn;
#line 331 "BB03AD.f"
	    dger_(n, n, &d__1, &dwork[*n + 1], &c__1, &dwork[1], &c__1, &a[
		    a_offset], lda);
/*         H2 * X */
#line 333 "BB03AD.f"
	    dgemv_("T", n, n, &c_b9, &x[x_offset], ldx, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 334 "BB03AD.f"
	    d__1 = -twobyn;
#line 334 "BB03AD.f"
	    dger_(n, n, &d__1, &dwork[1], &c__1, &dwork[*n + 1], &c__1, &x[
		    x_offset], ldx);
/*         X * H2 */
#line 336 "BB03AD.f"
	    dgemv_("N", n, n, &c_b9, &x[x_offset], ldx, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 337 "BB03AD.f"
	    d__1 = -twobyn;
#line 337 "BB03AD.f"
	    dger_(n, n, &d__1, &dwork[*n + 1], &c__1, &dwork[1], &c__1, &x[
		    x_offset], ldx);
/*         B * H2 */
#line 339 "BB03AD.f"
	    d__1 = -twobyn * ddot_(n, &b[b_offset], ldb, &dwork[1], &c__1);
#line 339 "BB03AD.f"
	    daxpy_(n, &d__1, &dwork[1], &c__1, &b[b_offset], ldb);
/*         Y = -B' * B */
#line 342 "BB03AD.f"
	    dger_(n, n, &c_b84, &b[b_offset], ldb, &b[b_offset], ldb, &y[
		    y_offset], ldy);

#line 344 "BB03AD.f"
	} else if (nr[2] == 2) {
#line 345 "BB03AD.f"
	    s_copy(note, "CTLEX: Example 4.2", (ftnlen)70, (ftnlen)18);
#line 346 "BB03AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 347 "BB03AD.f"
		ipar[1] = 10;
#line 348 "BB03AD.f"
		dpar[1] = -.5;
#line 349 "BB03AD.f"
		dpar[2] = 1.5;
#line 350 "BB03AD.f"
	    }
#line 351 "BB03AD.f"
	    if (dpar[1] >= 0. || dpar[2] <= 1.) {
#line 351 "BB03AD.f"
		*info = -3;
#line 351 "BB03AD.f"
	    }
#line 352 "BB03AD.f"
	    if (ipar[1] < 2) {
#line 352 "BB03AD.f"
		*info = -4;
#line 352 "BB03AD.f"
	    }
#line 353 "BB03AD.f"
	    *n = ipar[1];
#line 354 "BB03AD.f"
	    *m = 1;
#line 355 "BB03AD.f"
	    if (*lde < *n) {
#line 355 "BB03AD.f"
		*info = -9;
#line 355 "BB03AD.f"
	    }
#line 356 "BB03AD.f"
	    if (*lda < *n) {
#line 356 "BB03AD.f"
		*info = -11;
#line 356 "BB03AD.f"
	    }
#line 357 "BB03AD.f"
	    if (*ldy < *n) {
#line 357 "BB03AD.f"
		*info = -13;
#line 357 "BB03AD.f"
	    }
#line 358 "BB03AD.f"
	    if (*ldb < *m) {
#line 358 "BB03AD.f"
		*info = -15;
#line 358 "BB03AD.f"
	    }
#line 359 "BB03AD.f"
	    if (*ldwork < *n << 1) {
#line 359 "BB03AD.f"
		*info = -22;
#line 359 "BB03AD.f"
	    }
#line 360 "BB03AD.f"
	    if (*info != 0) {
#line 360 "BB03AD.f"
		return 0;
#line 360 "BB03AD.f"
	    }

#line 362 "BB03AD.f"
	    vec[6] = TRUE_;
#line 363 "BB03AD.f"
	    twobyn = 2. / (doublereal) (*n);
#line 364 "BB03AD.f"
	    dlaset_("A", n, n, &c_b8, &c_b9, &e[e_offset], lde, (ftnlen)1);
#line 365 "BB03AD.f"
	    dlaset_("A", n, n, &c_b8, &dpar[1], &a[a_offset], lda, (ftnlen)1);
#line 366 "BB03AD.f"
	    dlaset_("A", n, n, &c_b8, &c_b8, &y[y_offset], ldy, (ftnlen)1);
#line 367 "BB03AD.f"
	    d__1 = -twobyn;
#line 367 "BB03AD.f"
	    d__2 = 1. - twobyn;
#line 367 "BB03AD.f"
	    dlaset_("A", m, n, &d__1, &d__2, &b[b_offset], ldb, (ftnlen)1);
#line 368 "BB03AD.f"
	    i__1 = *n - 1;
#line 368 "BB03AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 369 "BB03AD.f"
		dwork[i__] = 1.;
#line 370 "BB03AD.f"
		a[i__ + (i__ + 1) * a_dim1] = 1.;
#line 371 "BB03AD.f"
/* L60: */
#line 371 "BB03AD.f"
	    }
#line 372 "BB03AD.f"
	    dwork[*n] = 1.;
/*         H1 * A */
#line 374 "BB03AD.f"
	    dgemv_("T", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 375 "BB03AD.f"
	    d__1 = -twobyn;
#line 375 "BB03AD.f"
	    dger_(n, n, &d__1, &dwork[1], &c__1, &dwork[*n + 1], &c__1, &a[
		    a_offset], lda);
/*         A * H1 */
#line 377 "BB03AD.f"
	    dgemv_("N", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 378 "BB03AD.f"
	    d__1 = -twobyn;
#line 378 "BB03AD.f"
	    dger_(n, n, &d__1, &dwork[*n + 1], &c__1, &dwork[1], &c__1, &a[
		    a_offset], lda);
/*         S A INV(S), B INV(S) */
#line 380 "BB03AD.f"
	    i__1 = *n;
#line 380 "BB03AD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 381 "BB03AD.f"
		i__2 = j - 1;
#line 381 "BB03AD.f"
		b[j * b_dim1 + 1] /= pow_di(&dpar[2], &i__2);
#line 382 "BB03AD.f"
		i__2 = *n;
#line 382 "BB03AD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 383 "BB03AD.f"
		    i__3 = i__ - j;
#line 383 "BB03AD.f"
		    a[i__ + j * a_dim1] *= pow_di(&dpar[2], &i__3);
#line 384 "BB03AD.f"
/* L70: */
#line 384 "BB03AD.f"
		}
#line 385 "BB03AD.f"
		dwork[j] = 1. - j % 2 * 2.;
#line 386 "BB03AD.f"
/* L80: */
#line 386 "BB03AD.f"
	    }
/*         H2 * A */
#line 388 "BB03AD.f"
	    dgemv_("T", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 389 "BB03AD.f"
	    d__1 = -twobyn;
#line 389 "BB03AD.f"
	    dger_(n, n, &d__1, &dwork[1], &c__1, &dwork[*n + 1], &c__1, &a[
		    a_offset], lda);
/*         A * H2 */
#line 391 "BB03AD.f"
	    dgemv_("N", n, n, &c_b9, &a[a_offset], lda, &dwork[1], &c__1, &
		    c_b8, &dwork[*n + 1], &c__1, (ftnlen)1);
#line 392 "BB03AD.f"
	    d__1 = -twobyn;
#line 392 "BB03AD.f"
	    dger_(n, n, &d__1, &dwork[*n + 1], &c__1, &dwork[1], &c__1, &a[
		    a_offset], lda);
/*         B * H2 */
#line 394 "BB03AD.f"
	    d__1 = -twobyn * ddot_(n, &b[b_offset], ldb, &dwork[1], &c__1);
#line 394 "BB03AD.f"
	    daxpy_(n, &d__1, &dwork[1], &c__1, &b[b_offset], ldb);
/*         Y = -B' * B */
#line 397 "BB03AD.f"
	    dger_(n, n, &c_b84, &b[b_offset], ldb, &b[b_offset], ldb, &y[
		    y_offset], ldy);

#line 399 "BB03AD.f"
	} else if (nr[2] == 3) {
#line 400 "BB03AD.f"
	    s_copy(note, "CTLEX: Example 4.3", (ftnlen)70, (ftnlen)18);
#line 401 "BB03AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 402 "BB03AD.f"
		ipar[1] = 10;
#line 403 "BB03AD.f"
		dpar[1] = 10.;
#line 404 "BB03AD.f"
	    }
#line 405 "BB03AD.f"
	    if (dpar[1] < 0.) {
#line 405 "BB03AD.f"
		*info = -3;
#line 405 "BB03AD.f"
	    }
#line 406 "BB03AD.f"
	    if (ipar[1] < 2) {
#line 406 "BB03AD.f"
		*info = -4;
#line 406 "BB03AD.f"
	    }
#line 407 "BB03AD.f"
	    *n = ipar[1];
#line 408 "BB03AD.f"
	    *m = 0;
#line 409 "BB03AD.f"
	    if (*lde < *n) {
#line 409 "BB03AD.f"
		*info = -9;
#line 409 "BB03AD.f"
	    }
#line 410 "BB03AD.f"
	    if (*lda < *n) {
#line 410 "BB03AD.f"
		*info = -11;
#line 410 "BB03AD.f"
	    }
#line 411 "BB03AD.f"
	    if (*ldy < *n) {
#line 411 "BB03AD.f"
		*info = -13;
#line 411 "BB03AD.f"
	    }
#line 412 "BB03AD.f"
	    if (*ldx < *n) {
#line 412 "BB03AD.f"
		*info = -17;
#line 412 "BB03AD.f"
	    }
#line 413 "BB03AD.f"
	    if (*info != 0) {
#line 413 "BB03AD.f"
		return 0;
#line 413 "BB03AD.f"
	    }

#line 415 "BB03AD.f"
	    vec[3] = TRUE_;
#line 416 "BB03AD.f"
	    vec[7] = TRUE_;
#line 417 "BB03AD.f"
	    d__1 = -dpar[1];
#line 417 "BB03AD.f"
	    temp = pow_dd(&c_b132, &d__1);
#line 418 "BB03AD.f"
	    dlaset_("U", n, n, &c_b8, &c_b8, &e[e_offset], lde, (ftnlen)1);
#line 419 "BB03AD.f"
	    dlaset_("L", n, n, &temp, &c_b9, &e[e_offset], lde, (ftnlen)1);
#line 420 "BB03AD.f"
	    dlaset_("L", n, n, &c_b8, &c_b8, &a[a_offset], lda, (ftnlen)1);
#line 421 "BB03AD.f"
	    dlaset_("U", n, n, &c_b9, &c_b8, &a[a_offset], lda, (ftnlen)1);
#line 422 "BB03AD.f"
	    dlaset_("A", n, n, &c_b9, &c_b9, &x[x_offset], ldx, (ftnlen)1);
#line 423 "BB03AD.f"
	    i__1 = *n;
#line 423 "BB03AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 424 "BB03AD.f"
		a[i__ + i__ * a_dim1] = (doublereal) (i__ - 1) + temp;
#line 425 "BB03AD.f"
/* L90: */
#line 425 "BB03AD.f"
	    }
/* Computing 2nd power */
#line 426 "BB03AD.f"
	    d__1 = temp;
#line 426 "BB03AD.f"
	    y[y_dim1 + 1] = temp * 2. + (doublereal) (*n - 1) * 2. * (d__1 * 
		    d__1);
/* Computing 2nd power */
#line 427 "BB03AD.f"
	    d__1 = temp;
#line 427 "BB03AD.f"
	    ttp1 = (doublereal) (*n + 1) * 2. * temp + 2. - d__1 * d__1;
/* Computing 2nd power */
#line 428 "BB03AD.f"
	    d__1 = temp;
#line 428 "BB03AD.f"
	    ttm1 = (doublereal) (*n - 1) * 2. * temp + 2. - d__1 * d__1;
#line 429 "BB03AD.f"
	    i__1 = *n;
#line 429 "BB03AD.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 430 "BB03AD.f"
		y[i__ + y_dim1] = y[y_dim1 + 1] + (doublereal) (i__ - 1) * 
			ttm1;
#line 431 "BB03AD.f"
/* L100: */
#line 431 "BB03AD.f"
	    }
#line 432 "BB03AD.f"
	    i__1 = *n;
#line 432 "BB03AD.f"
	    for (j = 2; j <= i__1; ++j) {
#line 433 "BB03AD.f"
		i__2 = *n;
#line 433 "BB03AD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 434 "BB03AD.f"
		    y[i__ + j * y_dim1] = y[i__ + y_dim1] + (doublereal) (j - 
			    1) * (ttp1 - i__ * 4. * temp);
#line 435 "BB03AD.f"
/* L110: */
#line 435 "BB03AD.f"
		}
#line 436 "BB03AD.f"
/* L120: */
#line 436 "BB03AD.f"
	    }

#line 438 "BB03AD.f"
	} else if (nr[2] == 4) {
#line 439 "BB03AD.f"
	    s_copy(note, "CTLEX: Example 4.4", (ftnlen)70, (ftnlen)18);
#line 440 "BB03AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 441 "BB03AD.f"
		ipar[1] = 10;
#line 442 "BB03AD.f"
		dpar[1] = 1.5;
#line 443 "BB03AD.f"
	    }
#line 444 "BB03AD.f"
	    if (dpar[1] < 1.) {
#line 444 "BB03AD.f"
		*info = -3;
#line 444 "BB03AD.f"
	    }
#line 445 "BB03AD.f"
	    if (ipar[1] < 1) {
#line 445 "BB03AD.f"
		*info = -4;
#line 445 "BB03AD.f"
	    }
#line 446 "BB03AD.f"
	    *n = ipar[1] * 3;
#line 447 "BB03AD.f"
	    *m = 1;
#line 448 "BB03AD.f"
	    if (*lde < *n) {
#line 448 "BB03AD.f"
		*info = -9;
#line 448 "BB03AD.f"
	    }
#line 449 "BB03AD.f"
	    if (*lda < *n) {
#line 449 "BB03AD.f"
		*info = -11;
#line 449 "BB03AD.f"
	    }
#line 450 "BB03AD.f"
	    if (*ldy < *n) {
#line 450 "BB03AD.f"
		*info = -13;
#line 450 "BB03AD.f"
	    }
#line 451 "BB03AD.f"
	    if (*ldb < *m) {
#line 451 "BB03AD.f"
		*info = -15;
#line 451 "BB03AD.f"
	    }
#line 452 "BB03AD.f"
	    if (*info != 0) {
#line 452 "BB03AD.f"
		return 0;
#line 452 "BB03AD.f"
	    }

#line 454 "BB03AD.f"
	    vec[3] = TRUE_;
#line 455 "BB03AD.f"
	    vec[6] = TRUE_;
#line 456 "BB03AD.f"
	    dlaset_("A", n, n, &c_b8, &c_b8, &e[e_offset], lde, (ftnlen)1);
#line 457 "BB03AD.f"
	    dlaset_("A", n, n, &c_b8, &c_b8, &a[a_offset], lda, (ftnlen)1);
#line 458 "BB03AD.f"
	    i__1 = ipar[1];
#line 458 "BB03AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 459 "BB03AD.f"
		temp = -pow_di(&dpar[1], &i__);
#line 460 "BB03AD.f"
		i__2 = i__ - 1;
#line 460 "BB03AD.f"
		for (j = 1; j <= i__2; ++j) {
#line 461 "BB03AD.f"
		    for (k = 0; k <= 2; ++k) {
#line 462 "BB03AD.f"
			a[*n - i__ * 3 + 3 + (j * 3 - k) * a_dim1] = temp;
#line 463 "BB03AD.f"
			a[*n - i__ * 3 + 2 + (j * 3 - k) * a_dim1] = temp * 
				2.;
#line 464 "BB03AD.f"
/* L130: */
#line 464 "BB03AD.f"
		    }
#line 465 "BB03AD.f"
/* L140: */
#line 465 "BB03AD.f"
		}
#line 466 "BB03AD.f"
		a[*n - i__ * 3 + 3 + (i__ * 3 - 2) * a_dim1] = temp;
#line 467 "BB03AD.f"
		a[*n - i__ * 3 + 2 + (i__ * 3 - 2) * a_dim1] = temp * 2.;
#line 468 "BB03AD.f"
		a[*n - i__ * 3 + 2 + (i__ * 3 - 1) * a_dim1] = temp * 2.;
#line 469 "BB03AD.f"
		a[*n - i__ * 3 + 2 + i__ * 3 * a_dim1] = temp;
#line 470 "BB03AD.f"
		a[*n - i__ * 3 + 1 + i__ * 3 * a_dim1] = temp;
#line 471 "BB03AD.f"
/* L150: */
#line 471 "BB03AD.f"
	    }
#line 472 "BB03AD.f"
	    i__1 = *n;
#line 472 "BB03AD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 473 "BB03AD.f"
		if (j > 1) {
#line 473 "BB03AD.f"
		    daxpy_(n, &c_b9, &a[j - 1 + a_dim1], lda, &a[j + a_dim1], 
			    lda);
#line 473 "BB03AD.f"
		}
#line 474 "BB03AD.f"
		b[j * b_dim1 + 1] = (doublereal) j;
#line 475 "BB03AD.f"
		i__2 = *n;
#line 475 "BB03AD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 476 "BB03AD.f"
		    e[i__ + (*n - j + 1) * e_dim1] = (doublereal) min(i__,j);
#line 477 "BB03AD.f"
		    y[i__ + j * y_dim1] = -((doublereal) (i__ * j));
#line 478 "BB03AD.f"
/* L160: */
#line 478 "BB03AD.f"
		}
#line 479 "BB03AD.f"
/* L170: */
#line 479 "BB03AD.f"
	    }

#line 481 "BB03AD.f"
	} else {
#line 482 "BB03AD.f"
	    *info = -2;
#line 483 "BB03AD.f"
	}
#line 484 "BB03AD.f"
    } else {
#line 485 "BB03AD.f"
	*info = -2;
#line 486 "BB03AD.f"
    }

#line 488 "BB03AD.f"
    return 0;
/* *** Last Line of BB03AD *** */
} /* bb03ad_ */

