#line 1 "BD01AD.f"
/* BD01AD.f -- translated by f2c (version 20100827).
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

#line 1 "BD01AD.f"
/* Table of constant values */

static doublereal c_b5 = 0.;
static doublereal c_b6 = 1.;
static integer c__1 = 1;
static integer c__5 = 5;

/* Subroutine */ int bd01ad_(char *def, integer *nr, doublereal *dpar, 
	integer *ipar, logical *vec, integer *n, integer *m, integer *p, 
	doublereal *e, integer *lde, doublereal *a, integer *lda, doublereal *
	b, integer *ldb, doublereal *c__, integer *ldc, doublereal *d__, 
	integer *ldd, char *note, doublereal *dwork, integer *ldwork, integer 
	*info, ftnlen def_len, ftnlen note_len)
{
    /* Initialized data */

    static logical vecdef[8] = { TRUE_,TRUE_,TRUE_,FALSE_,TRUE_,TRUE_,TRUE_,
	    FALSE_ };

    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, e_dim1, e_offset, i__1, i__2;
    doublereal d__1, d__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , f_open(olist *), s_rsle(cilist *), do_lio(integer *, integer *, 
	    char *, ftnlen), e_rsle(void), f_clos(cllist *);

    /* Local variables */
    static integer i__, j, l;
    static doublereal b1, b2, c1, c2, temp;
    static char dataf[12];
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal ttemp, appind;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static integer status;

    /* Fortran I/O blocks */
    static icilist io___4 = { 0, dataf, 0, "(A,I2.2,A)", 11, 1 };
    static cilist io___6 = { 1, 1, 1, 0, 0 };
    static cilist io___8 = { 1, 1, 1, 0, 0 };
    static cilist io___9 = { 1, 1, 1, 0, 0 };
    static cilist io___10 = { 1, 1, 1, 0, 0 };
    static cilist io___11 = { 1, 1, 1, 0, 0 };
    static cilist io___12 = { 1, 1, 1, 0, 0 };
    static icilist io___13 = { 0, dataf, 0, "(A,I1,A)", 12, 1 };
    static cilist io___14 = { 1, 1, 1, 0, 0 };
    static cilist io___15 = { 1, 1, 1, 0, 0 };
    static cilist io___16 = { 1, 1, 1, 0, 0 };
    static cilist io___17 = { 1, 1, 1, 0, 0 };
    static cilist io___18 = { 1, 1, 1, 0, 0 };
    static cilist io___19 = { 1, 1, 1, 0, 0 };
    static cilist io___21 = { 1, 1, 1, 0, 0 };



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

/*     To generate benchmark examples for time-invariant, */
/*     continuous-time dynamical systems */

/*         . */
/*       E x(t) = A x(t) + B u(t) */

/*         y(t) = C x(t) + D u(t) */

/*     E, A are real N-by-N matrices, B is N-by-M, C is P-by-N, and */
/*     D is P-by-M. In many examples, E is the identity matrix and D is */
/*     the zero matrix. */

/*     This routine is an implementation of the benchmark library */
/*     CTDSX (Version 1.0) described in [1]. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DEF     CHARACTER*1 */
/*             Specifies the kind of values used as parameters when */
/*             generating parameter-dependent and scalable examples */
/*             (i.e., examples with NR(1) = 2, 3, or 4): */
/*             = 'D':  Default values defined in [1] are used; */
/*             = 'N':  Values set in DPAR and IPAR are used. */
/*             This parameter is not referenced if NR(1) = 1. */
/*             Note that the scaling parameter of examples with */
/*             NR(1) = 3 or 4 is considered as a regular parameter in */
/*             this context. */

/*     Input/Output Parameters */

/*     NR      (input) INTEGER array, dimension (2) */
/*             Specifies the index of the desired example according */
/*             to [1]. */
/*             NR(1) defines the group: */
/*                   1 : parameter-free problems of fixed size */
/*                   2 : parameter-dependent problems of fixed size */
/*                   3 : parameter-free problems of scalable size */
/*                   4 : parameter-dependent problems of scalable size */
/*             NR(2) defines the number of the benchmark example */
/*             within a certain group according to [1]. */

/*     DPAR    (input/output) DOUBLE PRECISION array, dimension (7) */
/*             On entry, if DEF = 'N' and the desired example depends on */
/*             real parameters, then the array DPAR must contain the */
/*             values for these parameters. */
/*             For an explanation of the parameters see [1]. */
/*             For Examples 2.1 and 2.2, DPAR(1) defines the parameter */
/*             'epsilon'. */
/*             For Example 2.4, DPAR(1), ..., DPAR(7) define 'b', 'mu', */
/*             'r', 'r_c', 'k_l', 'sigma', 'a', respectively. */
/*             For Example 2.7, DPAR(1) and DPAR(2) define 'mu' and 'nu', */
/*             respectively. */
/*             For Example 4.1, DPAR(1), ..., DPAR(7) define 'a', 'b', */
/*             'c', 'beta_1', 'beta_2', 'gamma_1', 'gamma_2', */
/*             respectively. */
/*             For Example 4.2, DPAR(1), ..., DPAR(3) define 'mu', */
/*             'delta', 'kappa', respectively. */
/*             On exit, if DEF = 'D' and the desired example depends on */
/*             real parameters, then the array DPAR is overwritten by the */
/*             default values given in [1]. */

/*     IPAR    (input/output) INTEGER array, dimension (1) */
/*             On entry, if DEF = 'N' and the desired example depends on */
/*             integer parameters, then the array IPAR must contain the */
/*             values for these parameters. */
/*             For an explanation of the parameters see [1]. */
/*             For Examples 2.3, 2.5, and 2.6, IPAR(1) defines the */
/*             parameter 's'. */
/*             For Example 3.1, IPAR(1) defines 'q'. */
/*             For Examples 3.2 and 3.3, IPAR(1) defines 'n'. */
/*             For Example 3.4, IPAR(1) defines 'l'. */
/*             For Example 4.1, IPAR(1) defines 'n'. */
/*             For Example 4.2, IPAR(1) defines 'l'. */
/*             On exit, if DEF = 'D' and the desired example depends on */
/*             integer parameters, then the array IPAR is overwritten by */
/*             the default values given in [1]. */

/*     VEC     (output) LOGICAL array, dimension (8) */
/*             Flag vector which displays the availabilty of the output */
/*             data: */
/*             VEC(1), ..., VEC(3) refer to N, M, and P, respectively, */
/*             and are always .TRUE.. */
/*             VEC(4) is .TRUE. iff E is NOT the identity matrix. */
/*             VEC(5), ..., VEC(7) refer to A, B, and C, respectively, */
/*             and are always .TRUE.. */
/*             VEC(8) is .TRUE. iff D is NOT the zero matrix. */

/*     N       (output) INTEGER */
/*             The actual state dimension, i.e., the order of the */
/*             matrices E and A. */

/*     M       (output) INTEGER */
/*             The number of columns in the matrices B and D. */

/*     P       (output) INTEGER */
/*             The number of rows in the matrices C and D. */

/*     E       (output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             The leading N-by-N part of this array contains the */
/*             matrix E. */
/*             NOTE that this array is overwritten (by the identity */
/*             matrix), if VEC(4) = .FALSE.. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= N. */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array contains the */
/*             matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= N. */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array contains the */
/*             matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= N. */

/*     C       (output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array contains the */
/*             matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= P. */

/*     D       (output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading P-by-M part of this array contains the */
/*             matrix D. */
/*             NOTE that this array is overwritten (by the zero */
/*             matrix), if VEC(8) = .FALSE.. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= P. */

/*     NOTE    (output) CHARACTER*70 */
/*             String containing short information about the chosen */
/*             example. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             For Example 3.4, LDWORK >= 4*IPAR(1) is required. */
/*             For the other examples, no workspace is needed, i.e., */
/*             LDWORK >= 1. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; in particular, INFO = -3 or -4 indicates */
/*                   that at least one of the parameters in DPAR or */
/*                   IPAR, respectively, has an illegal value; */
/*             = 1:  data file can not be opened or has wrong format. */


/*     REFERENCES */

/*     [1]  Kressner, D., Mehrmann, V. and Penzl, T. */
/*          CTDSX - a Collection of Benchmark Examples for State-Space */
/*          Realizations of Continuous-Time Dynamical Systems. */
/*          SLICOT Working Note 1998-9. 1998. */

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

/*     continuous-time dynamical systems */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     . LAPACK . */
/*     .. External Subroutines .. */
/*     . BLAS . */
/*     . LAPACK . */
/*     .. Intrinsic Functions .. */
/*     .. Data Statements .. */
/*     . default values for availabities . */
#line 250 "BD01AD.f"
    /* Parameter adjustments */
#line 250 "BD01AD.f"
    --nr;
#line 250 "BD01AD.f"
    --dpar;
#line 250 "BD01AD.f"
    --ipar;
#line 250 "BD01AD.f"
    --vec;
#line 250 "BD01AD.f"
    e_dim1 = *lde;
#line 250 "BD01AD.f"
    e_offset = 1 + e_dim1;
#line 250 "BD01AD.f"
    e -= e_offset;
#line 250 "BD01AD.f"
    a_dim1 = *lda;
#line 250 "BD01AD.f"
    a_offset = 1 + a_dim1;
#line 250 "BD01AD.f"
    a -= a_offset;
#line 250 "BD01AD.f"
    b_dim1 = *ldb;
#line 250 "BD01AD.f"
    b_offset = 1 + b_dim1;
#line 250 "BD01AD.f"
    b -= b_offset;
#line 250 "BD01AD.f"
    c_dim1 = *ldc;
#line 250 "BD01AD.f"
    c_offset = 1 + c_dim1;
#line 250 "BD01AD.f"
    c__ -= c_offset;
#line 250 "BD01AD.f"
    d_dim1 = *ldd;
#line 250 "BD01AD.f"
    d_offset = 1 + d_dim1;
#line 250 "BD01AD.f"
    d__ -= d_offset;
#line 250 "BD01AD.f"
    --dwork;
#line 250 "BD01AD.f"

#line 250 "BD01AD.f"
    /* Function Body */

/*     .. Executable Statements .. */

#line 255 "BD01AD.f"
    *info = 0;
#line 256 "BD01AD.f"
    for (i__ = 1; i__ <= 8; ++i__) {
#line 257 "BD01AD.f"
	vec[i__] = vecdef[i__ - 1];
#line 258 "BD01AD.f"
/* L10: */
#line 258 "BD01AD.f"
    }

#line 260 "BD01AD.f"
    if (nr[1] == 1) {

#line 262 "BD01AD.f"
	if (nr[2] == 1) {
#line 263 "BD01AD.f"
	    s_copy(note, "Laub 1979, Ex.1", (ftnlen)70, (ftnlen)15);
#line 264 "BD01AD.f"
	    *n = 2;
#line 265 "BD01AD.f"
	    *m = 1;
#line 266 "BD01AD.f"
	    *p = 2;
#line 267 "BD01AD.f"
	    if (*lde < *n) {
#line 267 "BD01AD.f"
		*info = -10;
#line 267 "BD01AD.f"
	    }
#line 268 "BD01AD.f"
	    if (*lda < *n) {
#line 268 "BD01AD.f"
		*info = -12;
#line 268 "BD01AD.f"
	    }
#line 269 "BD01AD.f"
	    if (*ldb < *n) {
#line 269 "BD01AD.f"
		*info = -14;
#line 269 "BD01AD.f"
	    }
#line 270 "BD01AD.f"
	    if (*ldc < *p) {
#line 270 "BD01AD.f"
		*info = -16;
#line 270 "BD01AD.f"
	    }
#line 271 "BD01AD.f"
	    if (*ldd < *p) {
#line 271 "BD01AD.f"
		*info = -18;
#line 271 "BD01AD.f"
	    }
#line 272 "BD01AD.f"
	    if (*info != 0) {
#line 272 "BD01AD.f"
		return 0;
#line 272 "BD01AD.f"
	    }

#line 274 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 275 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
#line 276 "BD01AD.f"
	    a[(a_dim1 << 1) + 1] = 1.;
#line 277 "BD01AD.f"
	    b[b_dim1 + 1] = 0.;
#line 278 "BD01AD.f"
	    b[b_dim1 + 2] = 1.;
#line 279 "BD01AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
#line 280 "BD01AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 282 "BD01AD.f"
	} else if (nr[2] == 2) {
#line 283 "BD01AD.f"
	    s_copy(note, "Laub 1979, Ex.2: uncontrollable-unobservable data", 
		    (ftnlen)70, (ftnlen)49);
#line 284 "BD01AD.f"
	    *n = 2;
#line 285 "BD01AD.f"
	    *m = 1;
#line 286 "BD01AD.f"
	    *p = 1;
#line 287 "BD01AD.f"
	    if (*lde < *n) {
#line 287 "BD01AD.f"
		*info = -10;
#line 287 "BD01AD.f"
	    }
#line 288 "BD01AD.f"
	    if (*lda < *n) {
#line 288 "BD01AD.f"
		*info = -12;
#line 288 "BD01AD.f"
	    }
#line 289 "BD01AD.f"
	    if (*ldb < *n) {
#line 289 "BD01AD.f"
		*info = -14;
#line 289 "BD01AD.f"
	    }
#line 290 "BD01AD.f"
	    if (*ldc < *p) {
#line 290 "BD01AD.f"
		*info = -16;
#line 290 "BD01AD.f"
	    }
#line 291 "BD01AD.f"
	    if (*ldd < *p) {
#line 291 "BD01AD.f"
		*info = -18;
#line 291 "BD01AD.f"
	    }
#line 292 "BD01AD.f"
	    if (*info != 0) {
#line 292 "BD01AD.f"
		return 0;
#line 292 "BD01AD.f"
	    }

#line 294 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 295 "BD01AD.f"
	    a[a_dim1 + 1] = 4.;
#line 296 "BD01AD.f"
	    a[a_dim1 + 2] = -4.5;
#line 297 "BD01AD.f"
	    a[(a_dim1 << 1) + 1] = 3.;
#line 298 "BD01AD.f"
	    a[(a_dim1 << 1) + 2] = -3.5;
#line 299 "BD01AD.f"
	    b[b_dim1 + 1] = 1.;
#line 300 "BD01AD.f"
	    b[b_dim1 + 2] = -1.;
#line 301 "BD01AD.f"
	    c__[c_dim1 + 1] = 3.;
#line 302 "BD01AD.f"
	    c__[(c_dim1 << 1) + 1] = 2.;
#line 303 "BD01AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 305 "BD01AD.f"
	} else if (nr[2] == 3) {
#line 306 "BD01AD.f"
	    s_copy(note, "Beale/Shafai 1989: model of L-1011 aircraft", (
		    ftnlen)70, (ftnlen)43);
#line 307 "BD01AD.f"
	    *n = 4;
#line 308 "BD01AD.f"
	    *m = 2;
#line 309 "BD01AD.f"
	    *p = 4;
#line 310 "BD01AD.f"
	    if (*lde < *n) {
#line 310 "BD01AD.f"
		*info = -10;
#line 310 "BD01AD.f"
	    }
#line 311 "BD01AD.f"
	    if (*lda < *n) {
#line 311 "BD01AD.f"
		*info = -12;
#line 311 "BD01AD.f"
	    }
#line 312 "BD01AD.f"
	    if (*ldb < *n) {
#line 312 "BD01AD.f"
		*info = -14;
#line 312 "BD01AD.f"
	    }
#line 313 "BD01AD.f"
	    if (*ldc < *p) {
#line 313 "BD01AD.f"
		*info = -16;
#line 313 "BD01AD.f"
	    }
#line 314 "BD01AD.f"
	    if (*ldd < *p) {
#line 314 "BD01AD.f"
		*info = -18;
#line 314 "BD01AD.f"
	    }
#line 315 "BD01AD.f"
	    if (*info != 0) {
#line 315 "BD01AD.f"
		return 0;
#line 315 "BD01AD.f"
	    }

#line 317 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 318 "BD01AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
#line 319 "BD01AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 321 "BD01AD.f"
	} else if (nr[2] == 4) {
#line 322 "BD01AD.f"
	    s_copy(note, "Bhattacharyya et al. 1983: binary distillation col"\
		    "umn", (ftnlen)70, (ftnlen)53);
#line 323 "BD01AD.f"
	    *n = 8;
#line 324 "BD01AD.f"
	    *m = 2;
#line 325 "BD01AD.f"
	    *p = 8;
#line 326 "BD01AD.f"
	    if (*lde < *n) {
#line 326 "BD01AD.f"
		*info = -10;
#line 326 "BD01AD.f"
	    }
#line 327 "BD01AD.f"
	    if (*lda < *n) {
#line 327 "BD01AD.f"
		*info = -12;
#line 327 "BD01AD.f"
	    }
#line 328 "BD01AD.f"
	    if (*ldb < *n) {
#line 328 "BD01AD.f"
		*info = -14;
#line 328 "BD01AD.f"
	    }
#line 329 "BD01AD.f"
	    if (*ldc < *p) {
#line 329 "BD01AD.f"
		*info = -16;
#line 329 "BD01AD.f"
	    }
#line 330 "BD01AD.f"
	    if (*ldd < *p) {
#line 330 "BD01AD.f"
		*info = -18;
#line 330 "BD01AD.f"
	    }
#line 331 "BD01AD.f"
	    if (*info != 0) {
#line 331 "BD01AD.f"
		return 0;
#line 331 "BD01AD.f"
	    }

#line 333 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 334 "BD01AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
#line 335 "BD01AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 337 "BD01AD.f"
	} else if (nr[2] == 5) {
#line 338 "BD01AD.f"
	    s_copy(note, "Patnaik et al. 1980: tubular ammonia reactor", (
		    ftnlen)70, (ftnlen)44);
#line 339 "BD01AD.f"
	    *n = 9;
#line 340 "BD01AD.f"
	    *m = 3;
#line 341 "BD01AD.f"
	    *p = 9;
#line 342 "BD01AD.f"
	    if (*lde < *n) {
#line 342 "BD01AD.f"
		*info = -10;
#line 342 "BD01AD.f"
	    }
#line 343 "BD01AD.f"
	    if (*lda < *n) {
#line 343 "BD01AD.f"
		*info = -12;
#line 343 "BD01AD.f"
	    }
#line 344 "BD01AD.f"
	    if (*ldb < *n) {
#line 344 "BD01AD.f"
		*info = -14;
#line 344 "BD01AD.f"
	    }
#line 345 "BD01AD.f"
	    if (*ldc < *p) {
#line 345 "BD01AD.f"
		*info = -16;
#line 345 "BD01AD.f"
	    }
#line 346 "BD01AD.f"
	    if (*ldd < *p) {
#line 346 "BD01AD.f"
		*info = -18;
#line 346 "BD01AD.f"
	    }
#line 347 "BD01AD.f"
	    if (*info != 0) {
#line 347 "BD01AD.f"
		return 0;
#line 347 "BD01AD.f"
	    }

#line 349 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 350 "BD01AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
#line 351 "BD01AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 353 "BD01AD.f"
	} else if (nr[2] == 6) {
#line 354 "BD01AD.f"
	    s_copy(note, "Davison/Gesing 1978: J-100 jet engine", (ftnlen)70, 
		    (ftnlen)37);
#line 355 "BD01AD.f"
	    *n = 30;
#line 356 "BD01AD.f"
	    *m = 3;
#line 357 "BD01AD.f"
	    *p = 5;
#line 358 "BD01AD.f"
	    if (*lde < *n) {
#line 358 "BD01AD.f"
		*info = -10;
#line 358 "BD01AD.f"
	    }
#line 359 "BD01AD.f"
	    if (*lda < *n) {
#line 359 "BD01AD.f"
		*info = -12;
#line 359 "BD01AD.f"
	    }
#line 360 "BD01AD.f"
	    if (*ldb < *n) {
#line 360 "BD01AD.f"
		*info = -14;
#line 360 "BD01AD.f"
	    }
#line 361 "BD01AD.f"
	    if (*ldc < *p) {
#line 361 "BD01AD.f"
		*info = -16;
#line 361 "BD01AD.f"
	    }
#line 362 "BD01AD.f"
	    if (*ldd < *p) {
#line 362 "BD01AD.f"
		*info = -18;
#line 362 "BD01AD.f"
	    }
#line 363 "BD01AD.f"
	    if (*info != 0) {
#line 363 "BD01AD.f"
		return 0;
#line 363 "BD01AD.f"
	    }

#line 365 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 366 "BD01AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 368 "BD01AD.f"
	} else if (nr[2] == 7) {
#line 369 "BD01AD.f"
	    s_copy(note, "Davison 1967: binary distillation column", (ftnlen)
		    70, (ftnlen)40);
#line 370 "BD01AD.f"
	    *n = 11;
#line 371 "BD01AD.f"
	    *m = 3;
#line 372 "BD01AD.f"
	    *p = 3;
#line 373 "BD01AD.f"
	    if (*lde < *n) {
#line 373 "BD01AD.f"
		*info = -10;
#line 373 "BD01AD.f"
	    }
#line 374 "BD01AD.f"
	    if (*lda < *n) {
#line 374 "BD01AD.f"
		*info = -12;
#line 374 "BD01AD.f"
	    }
#line 375 "BD01AD.f"
	    if (*ldb < *n) {
#line 375 "BD01AD.f"
		*info = -14;
#line 375 "BD01AD.f"
	    }
#line 376 "BD01AD.f"
	    if (*ldc < *p) {
#line 376 "BD01AD.f"
		*info = -16;
#line 376 "BD01AD.f"
	    }
#line 377 "BD01AD.f"
	    if (*ldd < *p) {
#line 377 "BD01AD.f"
		*info = -18;
#line 377 "BD01AD.f"
	    }
#line 378 "BD01AD.f"
	    if (*info != 0) {
#line 378 "BD01AD.f"
		return 0;
#line 378 "BD01AD.f"
	    }

#line 380 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 381 "BD01AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
#line 382 "BD01AD.f"
	    c__[c_dim1 + 2] = 1.;
#line 383 "BD01AD.f"
	    c__[c_dim1 * 10 + 1] = 1.;
#line 384 "BD01AD.f"
	    c__[c_dim1 * 11 + 3] = 1.;
#line 385 "BD01AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);
#line 387 "BD01AD.f"
	} else if (nr[2] == 8) {
#line 388 "BD01AD.f"
	    s_copy(note, "Chien/Ergin/Ling/Lee 1958: drum boiler", (ftnlen)70,
		     (ftnlen)38);
#line 389 "BD01AD.f"
	    *n = 9;
#line 390 "BD01AD.f"
	    *m = 3;
#line 391 "BD01AD.f"
	    *p = 2;
#line 392 "BD01AD.f"
	    if (*lde < *n) {
#line 392 "BD01AD.f"
		*info = -10;
#line 392 "BD01AD.f"
	    }
#line 393 "BD01AD.f"
	    if (*lda < *n) {
#line 393 "BD01AD.f"
		*info = -12;
#line 393 "BD01AD.f"
	    }
#line 394 "BD01AD.f"
	    if (*ldb < *n) {
#line 394 "BD01AD.f"
		*info = -14;
#line 394 "BD01AD.f"
	    }
#line 395 "BD01AD.f"
	    if (*ldc < *p) {
#line 395 "BD01AD.f"
		*info = -16;
#line 395 "BD01AD.f"
	    }
#line 396 "BD01AD.f"
	    if (*ldd < *p) {
#line 396 "BD01AD.f"
		*info = -18;
#line 396 "BD01AD.f"
	    }
#line 397 "BD01AD.f"
	    if (*info != 0) {
#line 397 "BD01AD.f"
		return 0;
#line 397 "BD01AD.f"
	    }

#line 399 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 400 "BD01AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
#line 401 "BD01AD.f"
	    c__[c_dim1 * 6 + 1] = 1.;
#line 402 "BD01AD.f"
	    c__[c_dim1 * 9 + 2] = 1.;
#line 403 "BD01AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 405 "BD01AD.f"
	} else if (nr[2] == 9) {
#line 406 "BD01AD.f"
	    s_copy(note, "Ly, Gangsaas 1981: B-767 airplane", (ftnlen)70, (
		    ftnlen)33);
#line 407 "BD01AD.f"
	    *n = 55;
#line 408 "BD01AD.f"
	    *m = 2;
#line 409 "BD01AD.f"
	    *p = 2;
#line 410 "BD01AD.f"
	    if (*lde < *n) {
#line 410 "BD01AD.f"
		*info = -10;
#line 410 "BD01AD.f"
	    }
#line 411 "BD01AD.f"
	    if (*lda < *n) {
#line 411 "BD01AD.f"
		*info = -12;
#line 411 "BD01AD.f"
	    }
#line 412 "BD01AD.f"
	    if (*ldb < *n) {
#line 412 "BD01AD.f"
		*info = -14;
#line 412 "BD01AD.f"
	    }
#line 413 "BD01AD.f"
	    if (*ldc < *p) {
#line 413 "BD01AD.f"
		*info = -16;
#line 413 "BD01AD.f"
	    }
#line 414 "BD01AD.f"
	    if (*ldd < *p) {
#line 414 "BD01AD.f"
		*info = -18;
#line 414 "BD01AD.f"
	    }
#line 415 "BD01AD.f"
	    if (*info != 0) {
#line 415 "BD01AD.f"
		return 0;
#line 415 "BD01AD.f"
	    }

#line 417 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 418 "BD01AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 420 "BD01AD.f"
	} else if (nr[2] == 10) {
#line 421 "BD01AD.f"
	    s_copy(note, "control surface servo for an underwater vehicle", (
		    ftnlen)70, (ftnlen)47);
#line 422 "BD01AD.f"
	    *n = 8;
#line 423 "BD01AD.f"
	    *m = 2;
#line 424 "BD01AD.f"
	    *p = 1;
#line 425 "BD01AD.f"
	    if (*lde < *n) {
#line 425 "BD01AD.f"
		*info = -10;
#line 425 "BD01AD.f"
	    }
#line 426 "BD01AD.f"
	    if (*lda < *n) {
#line 426 "BD01AD.f"
		*info = -12;
#line 426 "BD01AD.f"
	    }
#line 427 "BD01AD.f"
	    if (*ldb < *n) {
#line 427 "BD01AD.f"
		*info = -14;
#line 427 "BD01AD.f"
	    }
#line 428 "BD01AD.f"
	    if (*ldc < *p) {
#line 428 "BD01AD.f"
		*info = -16;
#line 428 "BD01AD.f"
	    }
#line 429 "BD01AD.f"
	    if (*ldd < *p) {
#line 429 "BD01AD.f"
		*info = -18;
#line 429 "BD01AD.f"
	    }
#line 430 "BD01AD.f"
	    if (*info != 0) {
#line 430 "BD01AD.f"
		return 0;
#line 430 "BD01AD.f"
	    }

#line 432 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 433 "BD01AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
#line 434 "BD01AD.f"
	    c__[c_dim1 * 7 + 1] = 1.;
#line 435 "BD01AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);
#line 436 "BD01AD.f"
	} else {
#line 437 "BD01AD.f"
	    *info = -2;
#line 438 "BD01AD.f"
	}

#line 440 "BD01AD.f"
	if (nr[2] >= 3 && nr[2] <= 10) {
/*         .. loading data files */
#line 442 "BD01AD.f"
	    s_wsfi(&io___4);
#line 442 "BD01AD.f"
	    do_fio(&c__1, "BD011", (ftnlen)5);
#line 442 "BD01AD.f"
	    do_fio(&c__1, (char *)&nr[2], (ftnlen)sizeof(integer));
#line 442 "BD01AD.f"
	    do_fio(&c__1, ".dat", (ftnlen)4);
#line 442 "BD01AD.f"
	    e_wsfi();
#line 443 "BD01AD.f"
	    o__1.oerr = 1;
#line 443 "BD01AD.f"
	    o__1.ounit = 1;
#line 443 "BD01AD.f"
	    o__1.ofnmlen = 11;
#line 443 "BD01AD.f"
	    o__1.ofnm = dataf;
#line 443 "BD01AD.f"
	    o__1.orl = 0;
#line 443 "BD01AD.f"
	    o__1.osta = "OLD";
#line 443 "BD01AD.f"
	    o__1.oacc = 0;
#line 443 "BD01AD.f"
	    o__1.ofm = 0;
#line 443 "BD01AD.f"
	    o__1.oblnk = 0;
#line 443 "BD01AD.f"
	    status = f_open(&o__1);
#line 444 "BD01AD.f"
	    if (status != 0) {
#line 445 "BD01AD.f"
		*info = 1;
#line 446 "BD01AD.f"
	    } else {
#line 447 "BD01AD.f"
		i__1 = *n;
#line 447 "BD01AD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 448 "BD01AD.f"
		    status = s_rsle(&io___6);
#line 448 "BD01AD.f"
		    if (status != 0) {
#line 448 "BD01AD.f"
			goto L100001;
#line 448 "BD01AD.f"
		    }
#line 448 "BD01AD.f"
		    i__2 = *n;
#line 448 "BD01AD.f"
		    for (j = 1; j <= i__2; ++j) {
#line 448 "BD01AD.f"
			status = do_lio(&c__5, &c__1, (char *)&a[i__ + j * 
				a_dim1], (ftnlen)sizeof(doublereal));
#line 448 "BD01AD.f"
			if (status != 0) {
#line 448 "BD01AD.f"
			    goto L100001;
#line 448 "BD01AD.f"
			}
#line 448 "BD01AD.f"
		    }
#line 448 "BD01AD.f"
		    status = e_rsle();
#line 448 "BD01AD.f"
L100001:
#line 449 "BD01AD.f"
		    if (status != 0) {
#line 449 "BD01AD.f"
			*info = 1;
#line 449 "BD01AD.f"
		    }
#line 450 "BD01AD.f"
/* L110: */
#line 450 "BD01AD.f"
		}
#line 451 "BD01AD.f"
		i__1 = *n;
#line 451 "BD01AD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 452 "BD01AD.f"
		    status = s_rsle(&io___8);
#line 452 "BD01AD.f"
		    if (status != 0) {
#line 452 "BD01AD.f"
			goto L100002;
#line 452 "BD01AD.f"
		    }
#line 452 "BD01AD.f"
		    i__2 = *m;
#line 452 "BD01AD.f"
		    for (j = 1; j <= i__2; ++j) {
#line 452 "BD01AD.f"
			status = do_lio(&c__5, &c__1, (char *)&b[i__ + j * 
				b_dim1], (ftnlen)sizeof(doublereal));
#line 452 "BD01AD.f"
			if (status != 0) {
#line 452 "BD01AD.f"
			    goto L100002;
#line 452 "BD01AD.f"
			}
#line 452 "BD01AD.f"
		    }
#line 452 "BD01AD.f"
		    status = e_rsle();
#line 452 "BD01AD.f"
L100002:
#line 453 "BD01AD.f"
		    if (status != 0) {
#line 453 "BD01AD.f"
			*info = 1;
#line 453 "BD01AD.f"
		    }
#line 454 "BD01AD.f"
/* L120: */
#line 454 "BD01AD.f"
		}
#line 455 "BD01AD.f"
		if (nr[2] == 6 || nr[2] == 9) {
#line 456 "BD01AD.f"
		    i__1 = *p;
#line 456 "BD01AD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 457 "BD01AD.f"
			status = s_rsle(&io___9);
#line 457 "BD01AD.f"
			if (status != 0) {
#line 457 "BD01AD.f"
			    goto L100003;
#line 457 "BD01AD.f"
			}
#line 457 "BD01AD.f"
			i__2 = *n;
#line 457 "BD01AD.f"
			for (j = 1; j <= i__2; ++j) {
#line 457 "BD01AD.f"
			    status = do_lio(&c__5, &c__1, (char *)&c__[i__ + 
				    j * c_dim1], (ftnlen)sizeof(doublereal));
#line 457 "BD01AD.f"
			    if (status != 0) {
#line 457 "BD01AD.f"
				goto L100003;
#line 457 "BD01AD.f"
			    }
#line 457 "BD01AD.f"
			}
#line 457 "BD01AD.f"
			status = e_rsle();
#line 457 "BD01AD.f"
L100003:
#line 458 "BD01AD.f"
			if (status != 0) {
#line 458 "BD01AD.f"
			    *info = 1;
#line 458 "BD01AD.f"
			}
#line 459 "BD01AD.f"
/* L130: */
#line 459 "BD01AD.f"
		    }
#line 460 "BD01AD.f"
		}
#line 461 "BD01AD.f"
	    }
#line 462 "BD01AD.f"
	    cl__1.cerr = 0;
#line 462 "BD01AD.f"
	    cl__1.cunit = 1;
#line 462 "BD01AD.f"
	    cl__1.csta = 0;
#line 462 "BD01AD.f"
	    f_clos(&cl__1);
#line 463 "BD01AD.f"
	}

#line 465 "BD01AD.f"
    } else if (nr[1] == 2) {
#line 466 "BD01AD.f"
	if (! (lsame_(def, "D", (ftnlen)1, (ftnlen)1) || lsame_(def, "N", (
		ftnlen)1, (ftnlen)1))) {
#line 467 "BD01AD.f"
	    *info = -1;
#line 468 "BD01AD.f"
	    return 0;
#line 469 "BD01AD.f"
	}

#line 471 "BD01AD.f"
	if (nr[2] == 1) {
#line 472 "BD01AD.f"
	    s_copy(note, "Chow/Kokotovic 1976: magnetic tape control system", 
		    (ftnlen)70, (ftnlen)49);
#line 473 "BD01AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 473 "BD01AD.f"
		dpar[1] = 1e-6;
#line 473 "BD01AD.f"
	    }
#line 474 "BD01AD.f"
	    if (dpar[1] == 0.) {
#line 474 "BD01AD.f"
		*info = -3;
#line 474 "BD01AD.f"
	    }
#line 475 "BD01AD.f"
	    *n = 4;
#line 476 "BD01AD.f"
	    *m = 1;
#line 477 "BD01AD.f"
	    *p = 2;
#line 478 "BD01AD.f"
	    if (*lde < *n) {
#line 478 "BD01AD.f"
		*info = -10;
#line 478 "BD01AD.f"
	    }
#line 479 "BD01AD.f"
	    if (*lda < *n) {
#line 479 "BD01AD.f"
		*info = -12;
#line 479 "BD01AD.f"
	    }
#line 480 "BD01AD.f"
	    if (*ldb < *n) {
#line 480 "BD01AD.f"
		*info = -14;
#line 480 "BD01AD.f"
	    }
#line 481 "BD01AD.f"
	    if (*ldc < *p) {
#line 481 "BD01AD.f"
		*info = -16;
#line 481 "BD01AD.f"
	    }
#line 482 "BD01AD.f"
	    if (*ldd < *p) {
#line 482 "BD01AD.f"
		*info = -18;
#line 482 "BD01AD.f"
	    }
#line 483 "BD01AD.f"
	    if (*info != 0) {
#line 483 "BD01AD.f"
		return 0;
#line 483 "BD01AD.f"
	    }

#line 485 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 486 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
#line 487 "BD01AD.f"
	    a[(a_dim1 << 1) + 1] = .4;
#line 488 "BD01AD.f"
	    a[a_dim1 * 3 + 2] = .345;
#line 489 "BD01AD.f"
	    a[(a_dim1 << 1) + 3] = -.524 / dpar[1];
#line 490 "BD01AD.f"
	    a[a_dim1 * 3 + 3] = -.465 / dpar[1];
#line 491 "BD01AD.f"
	    a[(a_dim1 << 2) + 3] = .262 / dpar[1];
#line 492 "BD01AD.f"
	    a[(a_dim1 << 2) + 4] = -1. / dpar[1];
#line 493 "BD01AD.f"
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
#line 494 "BD01AD.f"
	    b[b_dim1 + 4] = 1. / dpar[1];
#line 495 "BD01AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
#line 496 "BD01AD.f"
	    c__[c_dim1 + 1] = 1.;
#line 497 "BD01AD.f"
	    c__[c_dim1 * 3 + 2] = 1.;
#line 498 "BD01AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 500 "BD01AD.f"
	} else if (nr[2] == 2) {
#line 501 "BD01AD.f"
	    s_copy(note, "Arnold/Laub 1984", (ftnlen)70, (ftnlen)16);
#line 502 "BD01AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 502 "BD01AD.f"
		dpar[1] = 1e-6;
#line 502 "BD01AD.f"
	    }
#line 503 "BD01AD.f"
	    *n = 4;
#line 504 "BD01AD.f"
	    *m = 1;
#line 505 "BD01AD.f"
	    *p = 1;
#line 506 "BD01AD.f"
	    if (*lde < *n) {
#line 506 "BD01AD.f"
		*info = -10;
#line 506 "BD01AD.f"
	    }
#line 507 "BD01AD.f"
	    if (*lda < *n) {
#line 507 "BD01AD.f"
		*info = -12;
#line 507 "BD01AD.f"
	    }
#line 508 "BD01AD.f"
	    if (*ldb < *n) {
#line 508 "BD01AD.f"
		*info = -14;
#line 508 "BD01AD.f"
	    }
#line 509 "BD01AD.f"
	    if (*ldc < *p) {
#line 509 "BD01AD.f"
		*info = -16;
#line 509 "BD01AD.f"
	    }
#line 510 "BD01AD.f"
	    if (*ldd < *p) {
#line 510 "BD01AD.f"
		*info = -18;
#line 510 "BD01AD.f"
	    }
#line 511 "BD01AD.f"
	    if (*info != 0) {
#line 511 "BD01AD.f"
		return 0;
#line 511 "BD01AD.f"
	    }

#line 513 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 514 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &dpar[1], &a[a_offset], lda, (ftnlen)1);
#line 515 "BD01AD.f"
	    a[a_dim1 + 1] = -dpar[1];
#line 516 "BD01AD.f"
	    a[a_dim1 + 2] = -1.;
#line 517 "BD01AD.f"
	    a[(a_dim1 << 1) + 1] = 1.;
#line 518 "BD01AD.f"
	    a[(a_dim1 << 1) + 2] = -dpar[1];
#line 519 "BD01AD.f"
	    a[a_dim1 * 3 + 4] = -1.;
#line 520 "BD01AD.f"
	    a[(a_dim1 << 2) + 3] = 1.;
#line 521 "BD01AD.f"
	    dlaset_("A", n, m, &c_b6, &c_b6, &b[b_offset], ldb, (ftnlen)1);
#line 522 "BD01AD.f"
	    dlaset_("A", p, n, &c_b6, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
#line 523 "BD01AD.f"
	    d__[d_dim1 + 1] = 0.;

#line 525 "BD01AD.f"
	} else if (nr[2] == 3) {
#line 526 "BD01AD.f"
	    s_copy(note, "Vertical acceleration of a rigid guided missile", (
		    ftnlen)70, (ftnlen)47);
#line 527 "BD01AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 527 "BD01AD.f"
		ipar[1] = 1;
#line 527 "BD01AD.f"
	    }
#line 528 "BD01AD.f"
	    if (ipar[1] < 1 || ipar[1] > 10) {
#line 528 "BD01AD.f"
		*info = -4;
#line 528 "BD01AD.f"
	    }
#line 529 "BD01AD.f"
	    *n = 3;
#line 530 "BD01AD.f"
	    *m = 1;
#line 531 "BD01AD.f"
	    *p = 1;
#line 532 "BD01AD.f"
	    if (*lde < *n) {
#line 532 "BD01AD.f"
		*info = -10;
#line 532 "BD01AD.f"
	    }
#line 533 "BD01AD.f"
	    if (*lda < *n) {
#line 533 "BD01AD.f"
		*info = -12;
#line 533 "BD01AD.f"
	    }
#line 534 "BD01AD.f"
	    if (*ldb < *n) {
#line 534 "BD01AD.f"
		*info = -14;
#line 534 "BD01AD.f"
	    }
#line 535 "BD01AD.f"
	    if (*ldc < *p) {
#line 535 "BD01AD.f"
		*info = -16;
#line 535 "BD01AD.f"
	    }
#line 536 "BD01AD.f"
	    if (*ldd < *p) {
#line 536 "BD01AD.f"
		*info = -18;
#line 536 "BD01AD.f"
	    }
#line 537 "BD01AD.f"
	    if (*info != 0) {
#line 537 "BD01AD.f"
		return 0;
#line 537 "BD01AD.f"
	    }

#line 539 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 540 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
#line 541 "BD01AD.f"
	    a[a_dim1 + 2] = 1.;
#line 542 "BD01AD.f"
	    a[a_dim1 * 3 + 3] = -190.;
#line 543 "BD01AD.f"
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
#line 544 "BD01AD.f"
	    b[b_dim1 + 3] = 190.;
#line 545 "BD01AD.f"
	    d__[d_dim1 + 1] = 0.;
#line 546 "BD01AD.f"
	    o__1.oerr = 1;
#line 546 "BD01AD.f"
	    o__1.ounit = 1;
#line 546 "BD01AD.f"
	    o__1.ofnmlen = 11;
#line 546 "BD01AD.f"
	    o__1.ofnm = "BD01203.dat";
#line 546 "BD01AD.f"
	    o__1.orl = 0;
#line 546 "BD01AD.f"
	    o__1.osta = "OLD";
#line 546 "BD01AD.f"
	    o__1.oacc = 0;
#line 546 "BD01AD.f"
	    o__1.ofm = 0;
#line 546 "BD01AD.f"
	    o__1.oblnk = 0;
#line 546 "BD01AD.f"
	    status = f_open(&o__1);
#line 547 "BD01AD.f"
	    if (status != 0) {
#line 548 "BD01AD.f"
		*info = 1;
#line 549 "BD01AD.f"
	    } else {
#line 550 "BD01AD.f"
		i__1 = ipar[1];
#line 550 "BD01AD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 551 "BD01AD.f"
		    status = s_rsle(&io___10);
#line 551 "BD01AD.f"
		    if (status != 0) {
#line 551 "BD01AD.f"
			goto L100004;
#line 551 "BD01AD.f"
		    }
#line 551 "BD01AD.f"
		    i__2 = *n;
#line 551 "BD01AD.f"
		    for (j = 1; j <= i__2; ++j) {
#line 551 "BD01AD.f"
			status = do_lio(&c__5, &c__1, (char *)&a[j * a_dim1 + 
				1], (ftnlen)sizeof(doublereal));
#line 551 "BD01AD.f"
			if (status != 0) {
#line 551 "BD01AD.f"
			    goto L100004;
#line 551 "BD01AD.f"
			}
#line 551 "BD01AD.f"
		    }
#line 551 "BD01AD.f"
		    status = e_rsle();
#line 551 "BD01AD.f"
L100004:
#line 552 "BD01AD.f"
		    if (status != 0) {
#line 552 "BD01AD.f"
			*info = 1;
#line 552 "BD01AD.f"
		    }
#line 553 "BD01AD.f"
		    status = s_rsle(&io___11);
#line 553 "BD01AD.f"
		    if (status != 0) {
#line 553 "BD01AD.f"
			goto L100005;
#line 553 "BD01AD.f"
		    }
#line 553 "BD01AD.f"
		    i__2 = *n;
#line 553 "BD01AD.f"
		    for (j = 2; j <= i__2; ++j) {
#line 553 "BD01AD.f"
			status = do_lio(&c__5, &c__1, (char *)&a[j * a_dim1 + 
				2], (ftnlen)sizeof(doublereal));
#line 553 "BD01AD.f"
			if (status != 0) {
#line 553 "BD01AD.f"
			    goto L100005;
#line 553 "BD01AD.f"
			}
#line 553 "BD01AD.f"
		    }
#line 553 "BD01AD.f"
		    status = e_rsle();
#line 553 "BD01AD.f"
L100005:
#line 554 "BD01AD.f"
		    if (status != 0) {
#line 554 "BD01AD.f"
			*info = 1;
#line 554 "BD01AD.f"
		    }
#line 555 "BD01AD.f"
		    status = s_rsle(&io___12);
#line 555 "BD01AD.f"
		    if (status != 0) {
#line 555 "BD01AD.f"
			goto L100006;
#line 555 "BD01AD.f"
		    }
#line 555 "BD01AD.f"
		    i__2 = *n;
#line 555 "BD01AD.f"
		    for (j = 1; j <= i__2; ++j) {
#line 555 "BD01AD.f"
			status = do_lio(&c__5, &c__1, (char *)&c__[j * c_dim1 
				+ 1], (ftnlen)sizeof(doublereal));
#line 555 "BD01AD.f"
			if (status != 0) {
#line 555 "BD01AD.f"
			    goto L100006;
#line 555 "BD01AD.f"
			}
#line 555 "BD01AD.f"
		    }
#line 555 "BD01AD.f"
		    status = e_rsle();
#line 555 "BD01AD.f"
L100006:
#line 556 "BD01AD.f"
		    if (status != 0) {
#line 556 "BD01AD.f"
			*info = 1;
#line 556 "BD01AD.f"
		    }
#line 557 "BD01AD.f"
/* L210: */
#line 557 "BD01AD.f"
		}
#line 558 "BD01AD.f"
	    }
#line 559 "BD01AD.f"
	    cl__1.cerr = 0;
#line 559 "BD01AD.f"
	    cl__1.cunit = 1;
#line 559 "BD01AD.f"
	    cl__1.csta = 0;
#line 559 "BD01AD.f"
	    f_clos(&cl__1);

#line 561 "BD01AD.f"
	} else if (nr[2] == 4) {
#line 562 "BD01AD.f"
	    s_copy(note, "Senning 1980: hydraulic positioning system", (
		    ftnlen)70, (ftnlen)42);
#line 563 "BD01AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 564 "BD01AD.f"
		dpar[1] = 1.4e4;
#line 565 "BD01AD.f"
		dpar[2] = .1287;
#line 566 "BD01AD.f"
		dpar[3] = .15;
#line 567 "BD01AD.f"
		dpar[4] = .01;
#line 568 "BD01AD.f"
		dpar[5] = .002;
#line 569 "BD01AD.f"
		dpar[6] = .24;
#line 570 "BD01AD.f"
		dpar[7] = 10.75;
#line 571 "BD01AD.f"
	    }
#line 572 "BD01AD.f"
	    if (dpar[1] <= 9e3 || dpar[1] >= 1.6e4 || (dpar[2] <= .05 || dpar[
		    2] >= .3) || (dpar[3] <= .05 || dpar[3] >= 5.) || (dpar[4]
		     <= 0. || dpar[4] >= .05) || (dpar[5] <= 1.03e-4 || dpar[
		    5] >= .0035) || (dpar[6] <= .001 || dpar[6] >= 15.) || (
		    dpar[7] <= 10.5 || dpar[7] >= 11.1)) {
#line 579 "BD01AD.f"
		*info = -3;
#line 580 "BD01AD.f"
	    }
#line 581 "BD01AD.f"
	    *n = 3;
#line 582 "BD01AD.f"
	    *m = 1;
#line 583 "BD01AD.f"
	    *p = 1;
#line 584 "BD01AD.f"
	    if (*lde < *n) {
#line 584 "BD01AD.f"
		*info = -10;
#line 584 "BD01AD.f"
	    }
#line 585 "BD01AD.f"
	    if (*lda < *n) {
#line 585 "BD01AD.f"
		*info = -12;
#line 585 "BD01AD.f"
	    }
#line 586 "BD01AD.f"
	    if (*ldb < *n) {
#line 586 "BD01AD.f"
		*info = -14;
#line 586 "BD01AD.f"
	    }
#line 587 "BD01AD.f"
	    if (*ldc < *p) {
#line 587 "BD01AD.f"
		*info = -16;
#line 587 "BD01AD.f"
	    }
#line 588 "BD01AD.f"
	    if (*ldd < *p) {
#line 588 "BD01AD.f"
		*info = -18;
#line 588 "BD01AD.f"
	    }
#line 589 "BD01AD.f"
	    if (*info != 0) {
#line 589 "BD01AD.f"
		return 0;
#line 589 "BD01AD.f"
	    }

#line 591 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 592 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
#line 593 "BD01AD.f"
	    a[(a_dim1 << 1) + 1] = 1.;
#line 594 "BD01AD.f"
	    a[(a_dim1 << 1) + 2] = -(dpar[3] + dpar[4] * 4. / 
		    3.141592653589793) / dpar[2];
#line 595 "BD01AD.f"
	    a[a_dim1 * 3 + 2] = dpar[7] / dpar[2];
#line 596 "BD01AD.f"
	    a[(a_dim1 << 1) + 3] = dpar[7] * -4. * dpar[1] / 874.;
#line 597 "BD01AD.f"
	    a[a_dim1 * 3 + 3] = dpar[1] * -4. * (dpar[6] + dpar[5]) / 874.;
#line 598 "BD01AD.f"
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
#line 599 "BD01AD.f"
	    b[b_dim1 + 3] = dpar[1] * -4. / 874.;
#line 600 "BD01AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
#line 601 "BD01AD.f"
	    d__[d_dim1 + 1] = 0.;

#line 603 "BD01AD.f"
	} else if (nr[2] == 5) {
#line 604 "BD01AD.f"
	    s_copy(note, "Kwakernaak/Westdyk 1985: cascade of inverted pendu"\
		    "la", (ftnlen)70, (ftnlen)52);
#line 605 "BD01AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 605 "BD01AD.f"
		ipar[1] = 1;
#line 605 "BD01AD.f"
	    }
#line 606 "BD01AD.f"
	    if (ipar[1] < 1 || ipar[1] > 7) {
#line 606 "BD01AD.f"
		*info = -4;
#line 606 "BD01AD.f"
	    }
#line 607 "BD01AD.f"
	    if (ipar[1] <= 6) {
#line 608 "BD01AD.f"
		*m = ipar[1];
#line 609 "BD01AD.f"
	    } else {
#line 610 "BD01AD.f"
		*m = 10;
#line 611 "BD01AD.f"
	    }
#line 612 "BD01AD.f"
	    *n = *m << 1;
#line 613 "BD01AD.f"
	    *p = *m;
#line 614 "BD01AD.f"
	    if (*lde < *n) {
#line 614 "BD01AD.f"
		*info = -10;
#line 614 "BD01AD.f"
	    }
#line 615 "BD01AD.f"
	    if (*lda < *n) {
#line 615 "BD01AD.f"
		*info = -12;
#line 615 "BD01AD.f"
	    }
#line 616 "BD01AD.f"
	    if (*ldb < *n) {
#line 616 "BD01AD.f"
		*info = -14;
#line 616 "BD01AD.f"
	    }
#line 617 "BD01AD.f"
	    if (*ldc < *p) {
#line 617 "BD01AD.f"
		*info = -16;
#line 617 "BD01AD.f"
	    }
#line 618 "BD01AD.f"
	    if (*ldd < *p) {
#line 618 "BD01AD.f"
		*info = -18;
#line 618 "BD01AD.f"
	    }
#line 619 "BD01AD.f"
	    if (*info != 0) {
#line 619 "BD01AD.f"
		return 0;
#line 619 "BD01AD.f"
	    }

#line 621 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 622 "BD01AD.f"
	    s_wsfi(&io___13);
#line 622 "BD01AD.f"
	    do_fio(&c__1, "BD01205", (ftnlen)7);
#line 622 "BD01AD.f"
	    do_fio(&c__1, (char *)&ipar[1], (ftnlen)sizeof(integer));
#line 622 "BD01AD.f"
	    do_fio(&c__1, ".dat", (ftnlen)4);
#line 622 "BD01AD.f"
	    e_wsfi();
#line 623 "BD01AD.f"
	    o__1.oerr = 1;
#line 623 "BD01AD.f"
	    o__1.ounit = 1;
#line 623 "BD01AD.f"
	    o__1.ofnmlen = 12;
#line 623 "BD01AD.f"
	    o__1.ofnm = dataf;
#line 623 "BD01AD.f"
	    o__1.orl = 0;
#line 623 "BD01AD.f"
	    o__1.osta = "OLD";
#line 623 "BD01AD.f"
	    o__1.oacc = 0;
#line 623 "BD01AD.f"
	    o__1.ofm = 0;
#line 623 "BD01AD.f"
	    o__1.oblnk = 0;
#line 623 "BD01AD.f"
	    status = f_open(&o__1);
#line 624 "BD01AD.f"
	    if (status != 0) {
#line 625 "BD01AD.f"
		*info = 1;
#line 626 "BD01AD.f"
	    } else {
#line 627 "BD01AD.f"
		i__1 = *n;
#line 627 "BD01AD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 628 "BD01AD.f"
		    status = s_rsle(&io___14);
#line 628 "BD01AD.f"
		    if (status != 0) {
#line 628 "BD01AD.f"
			goto L100007;
#line 628 "BD01AD.f"
		    }
#line 628 "BD01AD.f"
		    i__2 = *n;
#line 628 "BD01AD.f"
		    for (j = 1; j <= i__2; ++j) {
#line 628 "BD01AD.f"
			status = do_lio(&c__5, &c__1, (char *)&a[i__ + j * 
				a_dim1], (ftnlen)sizeof(doublereal));
#line 628 "BD01AD.f"
			if (status != 0) {
#line 628 "BD01AD.f"
			    goto L100007;
#line 628 "BD01AD.f"
			}
#line 628 "BD01AD.f"
		    }
#line 628 "BD01AD.f"
		    status = e_rsle();
#line 628 "BD01AD.f"
L100007:
#line 629 "BD01AD.f"
		    if (status != 0) {
#line 629 "BD01AD.f"
			*info = 1;
#line 629 "BD01AD.f"
		    }
#line 630 "BD01AD.f"
/* L220: */
#line 630 "BD01AD.f"
		}
#line 631 "BD01AD.f"
		i__1 = *n;
#line 631 "BD01AD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 632 "BD01AD.f"
		    status = s_rsle(&io___15);
#line 632 "BD01AD.f"
		    if (status != 0) {
#line 632 "BD01AD.f"
			goto L100008;
#line 632 "BD01AD.f"
		    }
#line 632 "BD01AD.f"
		    i__2 = *m;
#line 632 "BD01AD.f"
		    for (j = 1; j <= i__2; ++j) {
#line 632 "BD01AD.f"
			status = do_lio(&c__5, &c__1, (char *)&b[i__ + j * 
				b_dim1], (ftnlen)sizeof(doublereal));
#line 632 "BD01AD.f"
			if (status != 0) {
#line 632 "BD01AD.f"
			    goto L100008;
#line 632 "BD01AD.f"
			}
#line 632 "BD01AD.f"
		    }
#line 632 "BD01AD.f"
		    status = e_rsle();
#line 632 "BD01AD.f"
L100008:
#line 633 "BD01AD.f"
		    if (status != 0) {
#line 633 "BD01AD.f"
			*info = 1;
#line 633 "BD01AD.f"
		    }
#line 634 "BD01AD.f"
/* L230: */
#line 634 "BD01AD.f"
		}
#line 635 "BD01AD.f"
		i__1 = *p;
#line 635 "BD01AD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 636 "BD01AD.f"
		    status = s_rsle(&io___16);
#line 636 "BD01AD.f"
		    if (status != 0) {
#line 636 "BD01AD.f"
			goto L100009;
#line 636 "BD01AD.f"
		    }
#line 636 "BD01AD.f"
		    i__2 = *n;
#line 636 "BD01AD.f"
		    for (j = 1; j <= i__2; ++j) {
#line 636 "BD01AD.f"
			status = do_lio(&c__5, &c__1, (char *)&c__[i__ + j * 
				c_dim1], (ftnlen)sizeof(doublereal));
#line 636 "BD01AD.f"
			if (status != 0) {
#line 636 "BD01AD.f"
			    goto L100009;
#line 636 "BD01AD.f"
			}
#line 636 "BD01AD.f"
		    }
#line 636 "BD01AD.f"
		    status = e_rsle();
#line 636 "BD01AD.f"
L100009:
#line 637 "BD01AD.f"
		    if (status != 0) {
#line 637 "BD01AD.f"
			*info = 1;
#line 637 "BD01AD.f"
		    }
#line 638 "BD01AD.f"
/* L240: */
#line 638 "BD01AD.f"
		}
#line 639 "BD01AD.f"
	    }
#line 640 "BD01AD.f"
	    cl__1.cerr = 0;
#line 640 "BD01AD.f"
	    cl__1.cunit = 1;
#line 640 "BD01AD.f"
	    cl__1.csta = 0;
#line 640 "BD01AD.f"
	    f_clos(&cl__1);
#line 641 "BD01AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 643 "BD01AD.f"
	} else if (nr[2] == 6) {
#line 644 "BD01AD.f"
	    s_copy(note, "Kallstrom/Astrom 1981: regulation of a ship heading"
		    , (ftnlen)70, (ftnlen)51);
#line 645 "BD01AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 645 "BD01AD.f"
		ipar[1] = 1;
#line 645 "BD01AD.f"
	    }
#line 646 "BD01AD.f"
	    if (ipar[1] < 1 || ipar[1] > 5) {
#line 646 "BD01AD.f"
		*info = -4;
#line 646 "BD01AD.f"
	    }
#line 647 "BD01AD.f"
	    *n = 3;
#line 648 "BD01AD.f"
	    *m = 1;
#line 649 "BD01AD.f"
	    *p = 1;
#line 650 "BD01AD.f"
	    if (*lde < *n) {
#line 650 "BD01AD.f"
		*info = -10;
#line 650 "BD01AD.f"
	    }
#line 651 "BD01AD.f"
	    if (*lda < *n) {
#line 651 "BD01AD.f"
		*info = -12;
#line 651 "BD01AD.f"
	    }
#line 652 "BD01AD.f"
	    if (*ldb < *n) {
#line 652 "BD01AD.f"
		*info = -14;
#line 652 "BD01AD.f"
	    }
#line 653 "BD01AD.f"
	    if (*ldc < *p) {
#line 653 "BD01AD.f"
		*info = -16;
#line 653 "BD01AD.f"
	    }
#line 654 "BD01AD.f"
	    if (*ldd < *p) {
#line 654 "BD01AD.f"
		*info = -18;
#line 654 "BD01AD.f"
	    }
#line 655 "BD01AD.f"
	    if (*info != 0) {
#line 655 "BD01AD.f"
		return 0;
#line 655 "BD01AD.f"
	    }

#line 657 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 658 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
#line 659 "BD01AD.f"
	    a[(a_dim1 << 1) + 3] = 1.;
#line 660 "BD01AD.f"
	    b[b_dim1 + 3] = 0.;
#line 661 "BD01AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
#line 662 "BD01AD.f"
	    c__[c_dim1 * 3 + 1] = 1.;
#line 663 "BD01AD.f"
	    d__[d_dim1 + 1] = 0.;
#line 664 "BD01AD.f"
	    o__1.oerr = 1;
#line 664 "BD01AD.f"
	    o__1.ounit = 1;
#line 664 "BD01AD.f"
	    o__1.ofnmlen = 11;
#line 664 "BD01AD.f"
	    o__1.ofnm = "BD01206.dat";
#line 664 "BD01AD.f"
	    o__1.orl = 0;
#line 664 "BD01AD.f"
	    o__1.osta = "OLD";
#line 664 "BD01AD.f"
	    o__1.oacc = 0;
#line 664 "BD01AD.f"
	    o__1.ofm = 0;
#line 664 "BD01AD.f"
	    o__1.oblnk = 0;
#line 664 "BD01AD.f"
	    status = f_open(&o__1);
#line 665 "BD01AD.f"
	    if (status != 0) {
#line 666 "BD01AD.f"
		*info = 1;
#line 667 "BD01AD.f"
	    } else {
#line 668 "BD01AD.f"
		i__1 = ipar[1];
#line 668 "BD01AD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 669 "BD01AD.f"
		    status = s_rsle(&io___17);
#line 669 "BD01AD.f"
		    if (status != 0) {
#line 669 "BD01AD.f"
			goto L100010;
#line 669 "BD01AD.f"
		    }
#line 669 "BD01AD.f"
		    for (j = 1; j <= 2; ++j) {
#line 669 "BD01AD.f"
			status = do_lio(&c__5, &c__1, (char *)&a[j * a_dim1 + 
				1], (ftnlen)sizeof(doublereal));
#line 669 "BD01AD.f"
			if (status != 0) {
#line 669 "BD01AD.f"
			    goto L100010;
#line 669 "BD01AD.f"
			}
#line 669 "BD01AD.f"
		    }
#line 669 "BD01AD.f"
		    status = e_rsle();
#line 669 "BD01AD.f"
L100010:
#line 670 "BD01AD.f"
		    if (status != 0) {
#line 670 "BD01AD.f"
			*info = 1;
#line 670 "BD01AD.f"
		    }
#line 671 "BD01AD.f"
		    status = s_rsle(&io___18);
#line 671 "BD01AD.f"
		    if (status != 0) {
#line 671 "BD01AD.f"
			goto L100011;
#line 671 "BD01AD.f"
		    }
#line 671 "BD01AD.f"
		    for (j = 1; j <= 2; ++j) {
#line 671 "BD01AD.f"
			status = do_lio(&c__5, &c__1, (char *)&a[j * a_dim1 + 
				2], (ftnlen)sizeof(doublereal));
#line 671 "BD01AD.f"
			if (status != 0) {
#line 671 "BD01AD.f"
			    goto L100011;
#line 671 "BD01AD.f"
			}
#line 671 "BD01AD.f"
		    }
#line 671 "BD01AD.f"
		    status = e_rsle();
#line 671 "BD01AD.f"
L100011:
#line 672 "BD01AD.f"
		    if (status != 0) {
#line 672 "BD01AD.f"
			*info = 1;
#line 672 "BD01AD.f"
		    }
#line 673 "BD01AD.f"
		    status = s_rsle(&io___19);
#line 673 "BD01AD.f"
		    if (status != 0) {
#line 673 "BD01AD.f"
			goto L100012;
#line 673 "BD01AD.f"
		    }
#line 673 "BD01AD.f"
		    for (j = 1; j <= 2; ++j) {
#line 673 "BD01AD.f"
			status = do_lio(&c__5, &c__1, (char *)&b[j + b_dim1], 
				(ftnlen)sizeof(doublereal));
#line 673 "BD01AD.f"
			if (status != 0) {
#line 673 "BD01AD.f"
			    goto L100012;
#line 673 "BD01AD.f"
			}
#line 673 "BD01AD.f"
		    }
#line 673 "BD01AD.f"
		    status = e_rsle();
#line 673 "BD01AD.f"
L100012:
#line 674 "BD01AD.f"
		    if (status != 0) {
#line 674 "BD01AD.f"
			*info = 1;
#line 674 "BD01AD.f"
		    }
#line 675 "BD01AD.f"
/* L250: */
#line 675 "BD01AD.f"
		}
#line 676 "BD01AD.f"
	    }
#line 677 "BD01AD.f"
	    cl__1.cerr = 0;
#line 677 "BD01AD.f"
	    cl__1.cunit = 1;
#line 677 "BD01AD.f"
	    cl__1.csta = 0;
#line 677 "BD01AD.f"
	    f_clos(&cl__1);

#line 679 "BD01AD.f"
	} else if (nr[2] == 7) {
#line 680 "BD01AD.f"
	    s_copy(note, "Ackermann 1989: track-guided bus", (ftnlen)70, (
		    ftnlen)32);
#line 681 "BD01AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 682 "BD01AD.f"
		dpar[1] = 15.;
#line 683 "BD01AD.f"
		dpar[2] = 10.;
#line 684 "BD01AD.f"
	    }
#line 685 "BD01AD.f"
	    if (dpar[1] < 9.95 || dpar[1] > 16.) {
#line 685 "BD01AD.f"
		*info = -3;
#line 685 "BD01AD.f"
	    }
#line 686 "BD01AD.f"
	    if (dpar[1] < 1. || dpar[1] > 20.) {
#line 686 "BD01AD.f"
		*info = -3;
#line 686 "BD01AD.f"
	    }
#line 687 "BD01AD.f"
	    *n = 5;
#line 688 "BD01AD.f"
	    *m = 1;
#line 689 "BD01AD.f"
	    *p = 1;
#line 690 "BD01AD.f"
	    if (*lde < *n) {
#line 690 "BD01AD.f"
		*info = -10;
#line 690 "BD01AD.f"
	    }
#line 691 "BD01AD.f"
	    if (*lda < *n) {
#line 691 "BD01AD.f"
		*info = -12;
#line 691 "BD01AD.f"
	    }
#line 692 "BD01AD.f"
	    if (*ldb < *n) {
#line 692 "BD01AD.f"
		*info = -14;
#line 692 "BD01AD.f"
	    }
#line 693 "BD01AD.f"
	    if (*ldc < *p) {
#line 693 "BD01AD.f"
		*info = -16;
#line 693 "BD01AD.f"
	    }
#line 694 "BD01AD.f"
	    if (*ldd < *p) {
#line 694 "BD01AD.f"
		*info = -18;
#line 694 "BD01AD.f"
	    }
#line 695 "BD01AD.f"
	    if (*info != 0) {
#line 695 "BD01AD.f"
		return 0;
#line 695 "BD01AD.f"
	    }

#line 697 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 698 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
#line 699 "BD01AD.f"
	    a[a_dim1 + 1] = -668. / (dpar[1] * dpar[2]);
/* Computing 2nd power */
#line 700 "BD01AD.f"
	    d__1 = dpar[2];
#line 700 "BD01AD.f"
	    a[(a_dim1 << 1) + 1] = 180.4 / (dpar[1] * (d__1 * d__1)) - 1.;
#line 701 "BD01AD.f"
	    a[a_dim1 + 2] = 180.4 / (dpar[1] * 10.86);
#line 702 "BD01AD.f"
	    a[(a_dim1 << 1) + 2] = -4417.5452 / (dpar[1] * 10.86 * dpar[2]);
#line 703 "BD01AD.f"
	    a[a_dim1 * 5 + 1] = 198 / (dpar[1] * dpar[2]);
#line 704 "BD01AD.f"
	    a[a_dim1 * 5 + 2] = 726.66 / (dpar[1] * 10.86);
#line 705 "BD01AD.f"
	    a[a_dim1 + 3] = dpar[2];
#line 706 "BD01AD.f"
	    a[(a_dim1 << 2) + 3] = dpar[2];
#line 707 "BD01AD.f"
	    a[(a_dim1 << 1) + 4] = 1.;
#line 708 "BD01AD.f"
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
#line 709 "BD01AD.f"
	    b[b_dim1 + 5] = 1.;
#line 710 "BD01AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
#line 711 "BD01AD.f"
	    c__[c_dim1 * 3 + 1] = 1.;
#line 712 "BD01AD.f"
	    c__[(c_dim1 << 2) + 1] = 6.12;
#line 713 "BD01AD.f"
	    d__[d_dim1 + 1] = 0.;

#line 715 "BD01AD.f"
	} else {
#line 716 "BD01AD.f"
	    *info = -2;
#line 717 "BD01AD.f"
	}

#line 719 "BD01AD.f"
    } else if (nr[1] == 3) {
#line 720 "BD01AD.f"
	if (! (lsame_(def, "D", (ftnlen)1, (ftnlen)1) || lsame_(def, "N", (
		ftnlen)1, (ftnlen)1))) {
#line 721 "BD01AD.f"
	    *info = -1;
#line 722 "BD01AD.f"
	    return 0;
#line 723 "BD01AD.f"
	}

#line 725 "BD01AD.f"
	if (nr[2] == 1) {
#line 726 "BD01AD.f"
	    s_copy(note, "Laub 1979, Ex.4: string of high speed vehicles", (
		    ftnlen)70, (ftnlen)46);
#line 727 "BD01AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 727 "BD01AD.f"
		ipar[1] = 20;
#line 727 "BD01AD.f"
	    }
#line 728 "BD01AD.f"
	    if (ipar[1] < 2) {
#line 728 "BD01AD.f"
		*info = -4;
#line 728 "BD01AD.f"
	    }
#line 729 "BD01AD.f"
	    *n = (ipar[1] << 1) - 1;
#line 730 "BD01AD.f"
	    *m = ipar[1];
#line 731 "BD01AD.f"
	    *p = ipar[1] - 1;
#line 732 "BD01AD.f"
	    if (*lde < *n) {
#line 732 "BD01AD.f"
		*info = -10;
#line 732 "BD01AD.f"
	    }
#line 733 "BD01AD.f"
	    if (*lda < *n) {
#line 733 "BD01AD.f"
		*info = -12;
#line 733 "BD01AD.f"
	    }
#line 734 "BD01AD.f"
	    if (*ldb < *n) {
#line 734 "BD01AD.f"
		*info = -14;
#line 734 "BD01AD.f"
	    }
#line 735 "BD01AD.f"
	    if (*ldc < *p) {
#line 735 "BD01AD.f"
		*info = -16;
#line 735 "BD01AD.f"
	    }
#line 736 "BD01AD.f"
	    if (*ldd < *p) {
#line 736 "BD01AD.f"
		*info = -18;
#line 736 "BD01AD.f"
	    }
#line 737 "BD01AD.f"
	    if (*info != 0) {
#line 737 "BD01AD.f"
		return 0;
#line 737 "BD01AD.f"
	    }

#line 739 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 740 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
#line 741 "BD01AD.f"
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
#line 742 "BD01AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
#line 743 "BD01AD.f"
	    i__1 = *n;
#line 743 "BD01AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 744 "BD01AD.f"
		if (i__ % 2 == 1) {
#line 745 "BD01AD.f"
		    a[i__ + i__ * a_dim1] = -1.;
#line 746 "BD01AD.f"
		    b[i__ + (i__ + 1) / 2 * b_dim1] = 1.;
#line 747 "BD01AD.f"
		} else {
#line 748 "BD01AD.f"
		    a[i__ + (i__ - 1) * a_dim1] = 1.;
#line 749 "BD01AD.f"
		    a[i__ + (i__ + 1) * a_dim1] = -1.;
#line 750 "BD01AD.f"
		    c__[i__ / 2 + i__ * c_dim1] = 1.;
#line 751 "BD01AD.f"
		}
#line 752 "BD01AD.f"
/* L310: */
#line 752 "BD01AD.f"
	    }
#line 753 "BD01AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 755 "BD01AD.f"
	} else if (nr[2] == 2) {
#line 756 "BD01AD.f"
	    s_copy(note, "Hodel et al. 1996: heat flow in a thin rod", (
		    ftnlen)70, (ftnlen)42);
#line 757 "BD01AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 757 "BD01AD.f"
		ipar[1] = 100;
#line 757 "BD01AD.f"
	    }
#line 758 "BD01AD.f"
	    if (ipar[1] < 1) {
#line 758 "BD01AD.f"
		*info = -4;
#line 758 "BD01AD.f"
	    }
#line 759 "BD01AD.f"
	    *n = ipar[1];
#line 760 "BD01AD.f"
	    *m = 1;
#line 761 "BD01AD.f"
	    *p = *n;
#line 762 "BD01AD.f"
	    if (*lde < *n) {
#line 762 "BD01AD.f"
		*info = -10;
#line 762 "BD01AD.f"
	    }
#line 763 "BD01AD.f"
	    if (*lda < *n) {
#line 763 "BD01AD.f"
		*info = -12;
#line 763 "BD01AD.f"
	    }
#line 764 "BD01AD.f"
	    if (*ldb < *n) {
#line 764 "BD01AD.f"
		*info = -14;
#line 764 "BD01AD.f"
	    }
#line 765 "BD01AD.f"
	    if (*ldc < *p) {
#line 765 "BD01AD.f"
		*info = -16;
#line 765 "BD01AD.f"
	    }
#line 766 "BD01AD.f"
	    if (*ldd < *p) {
#line 766 "BD01AD.f"
		*info = -18;
#line 766 "BD01AD.f"
	    }
#line 767 "BD01AD.f"
	    if (*info != 0) {
#line 767 "BD01AD.f"
		return 0;
#line 767 "BD01AD.f"
	    }

#line 769 "BD01AD.f"
	    temp = (doublereal) (*n + 1);
#line 770 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 771 "BD01AD.f"
	    d__1 = temp * -2.;
#line 771 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &d__1, &a[a_offset], lda, (ftnlen)1);
#line 772 "BD01AD.f"
	    a[a_dim1 + 1] = -temp;
#line 773 "BD01AD.f"
	    i__1 = *n - 1;
#line 773 "BD01AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 774 "BD01AD.f"
		a[i__ + (i__ + 1) * a_dim1] = temp;
#line 775 "BD01AD.f"
		a[i__ + 1 + i__ * a_dim1] = temp;
#line 776 "BD01AD.f"
/* L320: */
#line 776 "BD01AD.f"
	    }
#line 777 "BD01AD.f"
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
#line 778 "BD01AD.f"
	    b[*n + b_dim1] = temp;
#line 779 "BD01AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
#line 780 "BD01AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 782 "BD01AD.f"
	} else if (nr[2] == 3) {
#line 783 "BD01AD.f"
	    s_copy(note, "Laub 1979, Ex.6", (ftnlen)70, (ftnlen)15);
#line 784 "BD01AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 784 "BD01AD.f"
		ipar[1] = 21;
#line 784 "BD01AD.f"
	    }
#line 785 "BD01AD.f"
	    if (ipar[1] < 1) {
#line 785 "BD01AD.f"
		*info = -4;
#line 785 "BD01AD.f"
	    }
#line 786 "BD01AD.f"
	    *n = ipar[1];
#line 787 "BD01AD.f"
	    *m = 1;
#line 788 "BD01AD.f"
	    *p = 1;
#line 789 "BD01AD.f"
	    if (*lde < *n) {
#line 789 "BD01AD.f"
		*info = -10;
#line 789 "BD01AD.f"
	    }
#line 790 "BD01AD.f"
	    if (*lda < *n) {
#line 790 "BD01AD.f"
		*info = -12;
#line 790 "BD01AD.f"
	    }
#line 791 "BD01AD.f"
	    if (*ldb < *n) {
#line 791 "BD01AD.f"
		*info = -14;
#line 791 "BD01AD.f"
	    }
#line 792 "BD01AD.f"
	    if (*ldc < *p) {
#line 792 "BD01AD.f"
		*info = -16;
#line 792 "BD01AD.f"
	    }
#line 793 "BD01AD.f"
	    if (*ldd < *p) {
#line 793 "BD01AD.f"
		*info = -18;
#line 793 "BD01AD.f"
	    }
#line 794 "BD01AD.f"
	    if (*info != 0) {
#line 794 "BD01AD.f"
		return 0;
#line 794 "BD01AD.f"
	    }

#line 796 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 797 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
#line 798 "BD01AD.f"
	    i__1 = *n - 1;
#line 798 "BD01AD.f"
	    i__2 = *n - 1;
#line 798 "BD01AD.f"
	    dlaset_("A", &i__1, &i__2, &c_b5, &c_b6, &a[(a_dim1 << 1) + 1], 
		    lda, (ftnlen)1);
#line 799 "BD01AD.f"
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
#line 800 "BD01AD.f"
	    b[*n + b_dim1] = 1.;
#line 801 "BD01AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
#line 802 "BD01AD.f"
	    c__[c_dim1 + 1] = 1.;
#line 803 "BD01AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 805 "BD01AD.f"
	} else if (nr[2] == 4) {
#line 806 "BD01AD.f"
	    s_copy(note, "Lang/Penzl 1994: rotating axle", (ftnlen)70, (
		    ftnlen)30);
#line 807 "BD01AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 807 "BD01AD.f"
		ipar[1] = 211;
#line 807 "BD01AD.f"
	    }
#line 808 "BD01AD.f"
	    if (ipar[1] < 1 || ipar[1] > 211) {
#line 808 "BD01AD.f"
		*info = -4;
#line 808 "BD01AD.f"
	    }
#line 809 "BD01AD.f"
	    *n = (ipar[1] << 1) - 1;
#line 810 "BD01AD.f"
	    *m = ipar[1];
#line 811 "BD01AD.f"
	    *p = ipar[1];
#line 812 "BD01AD.f"
	    if (*lde < *n) {
#line 812 "BD01AD.f"
		*info = -10;
#line 812 "BD01AD.f"
	    }
#line 813 "BD01AD.f"
	    if (*lda < *n) {
#line 813 "BD01AD.f"
		*info = -12;
#line 813 "BD01AD.f"
	    }
#line 814 "BD01AD.f"
	    if (*ldb < *n) {
#line 814 "BD01AD.f"
		*info = -14;
#line 814 "BD01AD.f"
	    }
#line 815 "BD01AD.f"
	    if (*ldc < *p) {
#line 815 "BD01AD.f"
		*info = -16;
#line 815 "BD01AD.f"
	    }
#line 816 "BD01AD.f"
	    if (*ldd < *p) {
#line 816 "BD01AD.f"
		*info = -18;
#line 816 "BD01AD.f"
	    }
#line 817 "BD01AD.f"
	    if (*ldwork < *m << 2) {
#line 817 "BD01AD.f"
		*info = -21;
#line 817 "BD01AD.f"
	    }
#line 818 "BD01AD.f"
	    if (*info != 0) {
#line 818 "BD01AD.f"
		return 0;
#line 818 "BD01AD.f"
	    }

#line 820 "BD01AD.f"
	    o__1.oerr = 1;
#line 820 "BD01AD.f"
	    o__1.ounit = 1;
#line 820 "BD01AD.f"
	    o__1.ofnmlen = 11;
#line 820 "BD01AD.f"
	    o__1.ofnm = "BD01304.dat";
#line 820 "BD01AD.f"
	    o__1.orl = 0;
#line 820 "BD01AD.f"
	    o__1.osta = "OLD";
#line 820 "BD01AD.f"
	    o__1.oacc = 0;
#line 820 "BD01AD.f"
	    o__1.ofm = 0;
#line 820 "BD01AD.f"
	    o__1.oblnk = 0;
#line 820 "BD01AD.f"
	    status = f_open(&o__1);
#line 821 "BD01AD.f"
	    if (status != 0) {
#line 822 "BD01AD.f"
		*info = 1;
#line 823 "BD01AD.f"
	    } else {
#line 824 "BD01AD.f"
		i__1 = *m << 2;
#line 824 "BD01AD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 825 "BD01AD.f"
		    status = s_rsle(&io___21);
#line 825 "BD01AD.f"
		    if (status != 0) {
#line 825 "BD01AD.f"
			goto L100013;
#line 825 "BD01AD.f"
		    }
#line 825 "BD01AD.f"
		    status = do_lio(&c__5, &c__1, (char *)&dwork[i__], (
			    ftnlen)sizeof(doublereal));
#line 825 "BD01AD.f"
		    if (status != 0) {
#line 825 "BD01AD.f"
			goto L100013;
#line 825 "BD01AD.f"
		    }
#line 825 "BD01AD.f"
		    status = e_rsle();
#line 825 "BD01AD.f"
L100013:
#line 826 "BD01AD.f"
		    if (status != 0) {
#line 826 "BD01AD.f"
			*info = 1;
#line 826 "BD01AD.f"
		    }
#line 827 "BD01AD.f"
/* L330: */
#line 827 "BD01AD.f"
		}
#line 828 "BD01AD.f"
	    }
#line 829 "BD01AD.f"
	    cl__1.cerr = 0;
#line 829 "BD01AD.f"
	    cl__1.cunit = 1;
#line 829 "BD01AD.f"
	    cl__1.csta = 0;
#line 829 "BD01AD.f"
	    f_clos(&cl__1);
#line 830 "BD01AD.f"
	    if (*info != 0) {
#line 830 "BD01AD.f"
		return 0;
#line 830 "BD01AD.f"
	    }
#line 831 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 832 "BD01AD.f"
	    e[e_dim1 + 1] = dwork[1];
#line 833 "BD01AD.f"
	    i__1 = *m;
#line 833 "BD01AD.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 834 "BD01AD.f"
		e[i__ + (i__ - 1) * e_dim1] = dwork[(i__ - 2 << 2) + 1];
#line 835 "BD01AD.f"
		e[i__ + i__ * e_dim1] = -dwork[(i__ - 1 << 2) + 1];
#line 836 "BD01AD.f"
/* L340: */
#line 836 "BD01AD.f"
	    }
#line 837 "BD01AD.f"
	    e[*m + *m * e_dim1] = -e[*m + *m * e_dim1];
#line 838 "BD01AD.f"
	    for (i__ = *m - 1; i__ >= 1; --i__) {
#line 839 "BD01AD.f"
		i__1 = *m;
#line 839 "BD01AD.f"
		for (j = i__; j <= i__1; ++j) {
#line 840 "BD01AD.f"
		    if (i__ == 1) {
#line 841 "BD01AD.f"
			e[j + i__ * e_dim1] -= e[j + (i__ + 1) * e_dim1];
#line 842 "BD01AD.f"
		    } else {
#line 843 "BD01AD.f"
			e[j + i__ * e_dim1] = e[j + (i__ + 1) * e_dim1] - e[j 
				+ i__ * e_dim1];
#line 844 "BD01AD.f"
		    }
#line 845 "BD01AD.f"
/* L345: */
#line 845 "BD01AD.f"
		}
#line 846 "BD01AD.f"
/* L350: */
#line 846 "BD01AD.f"
	    }
#line 847 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
#line 848 "BD01AD.f"
	    i__1 = *m;
#line 848 "BD01AD.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 849 "BD01AD.f"
		a[i__ - 1 + i__ * a_dim1] = dwork[(i__ - 2 << 2) + 3];
#line 850 "BD01AD.f"
		a[i__ + i__ * a_dim1] = dwork[(i__ - 2 << 2) + 3] * -2. - 
			dwork[(i__ - 1 << 2) + 2];
#line 851 "BD01AD.f"
		a[i__ + a_dim1] = dwork[(i__ - 1 << 2) + 2] - dwork[(i__ - 2 
			<< 2) + 2];
#line 852 "BD01AD.f"
		a[i__ - 1 + (*m + i__ - 1) * a_dim1] = dwork[(i__ - 1) * 4];
#line 853 "BD01AD.f"
		a[i__ + (*m + i__ - 1) * a_dim1] = dwork[(i__ - 1) * 4] * -2.;
#line 854 "BD01AD.f"
		if (i__ < *m) {
#line 855 "BD01AD.f"
		    a[i__ + 1 + i__ * a_dim1] = dwork[(i__ - 2 << 2) + 3];
#line 856 "BD01AD.f"
		    i__2 = *m;
#line 856 "BD01AD.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 857 "BD01AD.f"
			a[j + i__ * a_dim1] = a[j + i__ * a_dim1] + dwork[(j 
				- 2 << 2) + 2] - dwork[(j - 1 << 2) + 2];
#line 859 "BD01AD.f"
/* L355: */
#line 859 "BD01AD.f"
		    }
#line 860 "BD01AD.f"
		    a[i__ + 1 + (*m + i__ - 1) * a_dim1] = dwork[(i__ - 1) * 
			    4];
#line 861 "BD01AD.f"
		}
#line 862 "BD01AD.f"
/* L360: */
#line 862 "BD01AD.f"
	    }
#line 863 "BD01AD.f"
	    a[a_dim1 + 1] = -dwork[2];
#line 864 "BD01AD.f"
	    a[(a_dim1 << 1) + 1] = -dwork[3];
#line 865 "BD01AD.f"
	    a[(*m + 1) * a_dim1 + 1] = -a[(*m + 1) * a_dim1 + 1];
#line 866 "BD01AD.f"
	    i__1 = *m - 1;
#line 866 "BD01AD.f"
	    i__2 = *m - 1;
#line 866 "BD01AD.f"
	    dlaset_("A", &i__1, &i__2, &c_b5, &c_b6, &a[*m + 1 + (a_dim1 << 1)
		    ], lda, (ftnlen)1);
#line 867 "BD01AD.f"
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
#line 868 "BD01AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
#line 869 "BD01AD.f"
	    i__1 = *m;
#line 869 "BD01AD.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 870 "BD01AD.f"
		b[i__ + i__ * b_dim1] = -1.;
#line 871 "BD01AD.f"
		b[i__ + (i__ - 1) * b_dim1] = 1.;
#line 872 "BD01AD.f"
		c__[i__ + i__ * c_dim1] = dwork[(i__ - 2 << 2) + 3];
#line 873 "BD01AD.f"
		c__[i__ + (*m + i__ - 1) * c_dim1] = dwork[(i__ - 1) * 4];
#line 874 "BD01AD.f"
/* L370: */
#line 874 "BD01AD.f"
	    }
#line 875 "BD01AD.f"
	    b[b_dim1 + 1] = 1.;
#line 876 "BD01AD.f"
	    c__[c_dim1 + 1] = 1.;
#line 877 "BD01AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 879 "BD01AD.f"
	} else {
#line 880 "BD01AD.f"
	    *info = -2;
#line 881 "BD01AD.f"
	}

#line 883 "BD01AD.f"
    } else if (nr[1] == 4) {
#line 884 "BD01AD.f"
	if (! (lsame_(def, "D", (ftnlen)1, (ftnlen)1) || lsame_(def, "N", (
		ftnlen)1, (ftnlen)1))) {
#line 885 "BD01AD.f"
	    *info = -1;
#line 886 "BD01AD.f"
	    return 0;
#line 887 "BD01AD.f"
	}

#line 889 "BD01AD.f"
	if (nr[2] == 1) {
#line 890 "BD01AD.f"
	    s_copy(note, "Rosen/Wang 1995: control of 1-dim. heat flow", (
		    ftnlen)70, (ftnlen)44);
#line 891 "BD01AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 892 "BD01AD.f"
		ipar[1] = 100;
#line 893 "BD01AD.f"
		dpar[1] = .01;
#line 894 "BD01AD.f"
		dpar[2] = 1.;
#line 895 "BD01AD.f"
		dpar[3] = 1.;
#line 896 "BD01AD.f"
		dpar[4] = .2;
#line 897 "BD01AD.f"
		dpar[5] = .3;
#line 898 "BD01AD.f"
		dpar[6] = .2;
#line 899 "BD01AD.f"
		dpar[7] = .3;
#line 900 "BD01AD.f"
	    }
#line 901 "BD01AD.f"
	    if (ipar[1] < 2) {
#line 901 "BD01AD.f"
		*info = -4;
#line 901 "BD01AD.f"
	    }
#line 902 "BD01AD.f"
	    *n = ipar[1];
#line 903 "BD01AD.f"
	    *m = 1;
#line 904 "BD01AD.f"
	    *p = 1;
#line 905 "BD01AD.f"
	    if (*lde < *n) {
#line 905 "BD01AD.f"
		*info = -10;
#line 905 "BD01AD.f"
	    }
#line 906 "BD01AD.f"
	    if (*lda < *n) {
#line 906 "BD01AD.f"
		*info = -12;
#line 906 "BD01AD.f"
	    }
#line 907 "BD01AD.f"
	    if (*ldb < *n) {
#line 907 "BD01AD.f"
		*info = -14;
#line 907 "BD01AD.f"
	    }
#line 908 "BD01AD.f"
	    if (*ldc < *p) {
#line 908 "BD01AD.f"
		*info = -16;
#line 908 "BD01AD.f"
	    }
#line 909 "BD01AD.f"
	    if (*ldd < *p) {
#line 909 "BD01AD.f"
		*info = -18;
#line 909 "BD01AD.f"
	    }
#line 910 "BD01AD.f"
	    if (*info != 0) {
#line 910 "BD01AD.f"
		return 0;
#line 910 "BD01AD.f"
	    }

#line 912 "BD01AD.f"
	    vec[4] = TRUE_;
#line 913 "BD01AD.f"
	    appind = (doublereal) (*n + 1);
#line 914 "BD01AD.f"
	    ttemp = -dpar[1] * appind;
#line 915 "BD01AD.f"
	    temp = 1 / (appind * 6.);
#line 916 "BD01AD.f"
	    d__1 = temp * 4.;
#line 916 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &d__1, &e[e_offset], lde, (ftnlen)1);
#line 917 "BD01AD.f"
	    d__1 = ttemp * 2.;
#line 917 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &d__1, &a[a_offset], lda, (ftnlen)1);
#line 918 "BD01AD.f"
	    i__1 = *n - 1;
#line 918 "BD01AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 919 "BD01AD.f"
		a[i__ + 1 + i__ * a_dim1] = -ttemp;
#line 920 "BD01AD.f"
		a[i__ + (i__ + 1) * a_dim1] = -ttemp;
#line 921 "BD01AD.f"
		e[i__ + 1 + i__ * e_dim1] = temp;
#line 922 "BD01AD.f"
		e[i__ + (i__ + 1) * e_dim1] = temp;
#line 923 "BD01AD.f"
/* L410: */
#line 923 "BD01AD.f"
	    }
#line 924 "BD01AD.f"
	    i__1 = *n;
#line 924 "BD01AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 925 "BD01AD.f"
		d__1 = (doublereal) (i__ - 1) / appind;
#line 925 "BD01AD.f"
		b1 = max(d__1,dpar[4]);
/* Computing MIN */
#line 926 "BD01AD.f"
		d__1 = (doublereal) (i__ + 1) / appind;
#line 926 "BD01AD.f"
		b2 = min(d__1,dpar[5]);
/* Computing MAX */
#line 927 "BD01AD.f"
		d__1 = (doublereal) (i__ - 1) / appind;
#line 927 "BD01AD.f"
		c1 = max(d__1,dpar[6]);
/* Computing MIN */
#line 928 "BD01AD.f"
		d__1 = (doublereal) (i__ + 1) / appind;
#line 928 "BD01AD.f"
		c2 = min(d__1,dpar[7]);
#line 929 "BD01AD.f"
		if (b1 >= b2) {
#line 930 "BD01AD.f"
		    b[i__ + b_dim1] = 0.;
#line 931 "BD01AD.f"
		} else {
#line 932 "BD01AD.f"
		    b[i__ + b_dim1] = b2 - b1;
/* Computing MIN */
#line 933 "BD01AD.f"
		    d__1 = b2, d__2 = (doublereal) i__ / appind;
#line 933 "BD01AD.f"
		    temp = min(d__1,d__2);
#line 934 "BD01AD.f"
		    if (b1 < temp) {
/* Computing 2nd power */
#line 935 "BD01AD.f"
			d__1 = temp;
/* Computing 2nd power */
#line 935 "BD01AD.f"
			d__2 = b1;
#line 935 "BD01AD.f"
			b[i__ + b_dim1] += appind * (d__1 * d__1 - d__2 * 
				d__2) / 2.;
#line 936 "BD01AD.f"
			b[i__ + b_dim1] += (doublereal) i__ * (b1 - temp);
#line 937 "BD01AD.f"
		    }
/* Computing MAX */
#line 938 "BD01AD.f"
		    d__1 = b1, d__2 = (doublereal) i__ / appind;
#line 938 "BD01AD.f"
		    temp = max(d__1,d__2);
#line 939 "BD01AD.f"
		    if (temp < b2) {
/* Computing 2nd power */
#line 940 "BD01AD.f"
			d__1 = b2;
/* Computing 2nd power */
#line 940 "BD01AD.f"
			d__2 = temp;
#line 940 "BD01AD.f"
			b[i__ + b_dim1] -= appind * (d__1 * d__1 - d__2 * 
				d__2) / 2.;
#line 941 "BD01AD.f"
			b[i__ + b_dim1] -= (doublereal) i__ * (temp - b2);
#line 942 "BD01AD.f"
		    }
#line 943 "BD01AD.f"
		}
#line 944 "BD01AD.f"
		if (c1 >= c2) {
#line 945 "BD01AD.f"
		    c__[i__ * c_dim1 + 1] = 0.;
#line 946 "BD01AD.f"
		} else {
#line 947 "BD01AD.f"
		    c__[i__ * c_dim1 + 1] = c2 - c1;
/* Computing MIN */
#line 948 "BD01AD.f"
		    d__1 = c2, d__2 = (doublereal) i__ / appind;
#line 948 "BD01AD.f"
		    temp = min(d__1,d__2);
#line 949 "BD01AD.f"
		    if (c1 < temp) {
/* Computing 2nd power */
#line 950 "BD01AD.f"
			d__1 = temp;
/* Computing 2nd power */
#line 950 "BD01AD.f"
			d__2 = c1;
#line 950 "BD01AD.f"
			c__[i__ * c_dim1 + 1] += appind * (d__1 * d__1 - d__2 
				* d__2) / 2.;
#line 951 "BD01AD.f"
			c__[i__ * c_dim1 + 1] += (doublereal) i__ * (c1 - 
				temp);
#line 952 "BD01AD.f"
		    }
/* Computing MAX */
#line 953 "BD01AD.f"
		    d__1 = c1, d__2 = (doublereal) i__ / appind;
#line 953 "BD01AD.f"
		    temp = max(d__1,d__2);
#line 954 "BD01AD.f"
		    if (temp < c2) {
/* Computing 2nd power */
#line 955 "BD01AD.f"
			d__1 = c2;
/* Computing 2nd power */
#line 955 "BD01AD.f"
			d__2 = temp;
#line 955 "BD01AD.f"
			c__[i__ * c_dim1 + 1] -= appind * (d__1 * d__1 - d__2 
				* d__2) / 2.;
#line 956 "BD01AD.f"
			c__[i__ * c_dim1 + 1] -= (doublereal) i__ * (temp - 
				c2);
#line 957 "BD01AD.f"
		    }
#line 958 "BD01AD.f"
		}
#line 959 "BD01AD.f"
/* L420: */
#line 959 "BD01AD.f"
	    }
#line 960 "BD01AD.f"
	    dscal_(n, &dpar[2], &b[b_dim1 + 1], &c__1);
#line 961 "BD01AD.f"
	    dscal_(n, &dpar[3], &c__[c_dim1 + 1], ldc);
#line 962 "BD01AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 964 "BD01AD.f"
	} else if (nr[2] == 2) {
#line 965 "BD01AD.f"
	    s_copy(note, "Hench et al. 1995: coupled springs, dashpots, mass"\
		    "es", (ftnlen)70, (ftnlen)52);
#line 966 "BD01AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 967 "BD01AD.f"
		ipar[1] = 30;
#line 968 "BD01AD.f"
		dpar[1] = 4.;
#line 969 "BD01AD.f"
		dpar[2] = 4.;
#line 970 "BD01AD.f"
		dpar[3] = 1.;
#line 971 "BD01AD.f"
	    }
#line 972 "BD01AD.f"
	    if (ipar[1] < 2) {
#line 972 "BD01AD.f"
		*info = -4;
#line 972 "BD01AD.f"
	    }
#line 973 "BD01AD.f"
	    l = ipar[1];
#line 974 "BD01AD.f"
	    *n = l << 1;
#line 975 "BD01AD.f"
	    *m = 2;
#line 976 "BD01AD.f"
	    *p = l << 1;
#line 977 "BD01AD.f"
	    if (*lde < *n) {
#line 977 "BD01AD.f"
		*info = -10;
#line 977 "BD01AD.f"
	    }
#line 978 "BD01AD.f"
	    if (*lda < *n) {
#line 978 "BD01AD.f"
		*info = -12;
#line 978 "BD01AD.f"
	    }
#line 979 "BD01AD.f"
	    if (*ldb < *n) {
#line 979 "BD01AD.f"
		*info = -14;
#line 979 "BD01AD.f"
	    }
#line 980 "BD01AD.f"
	    if (*ldc < *p) {
#line 980 "BD01AD.f"
		*info = -16;
#line 980 "BD01AD.f"
	    }
#line 981 "BD01AD.f"
	    if (*ldd < *p) {
#line 981 "BD01AD.f"
		*info = -18;
#line 981 "BD01AD.f"
	    }
#line 982 "BD01AD.f"
	    if (*info != 0) {
#line 982 "BD01AD.f"
		return 0;
#line 982 "BD01AD.f"
	    }

#line 984 "BD01AD.f"
	    vec[4] = TRUE_;
#line 985 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &dpar[1], &e[e_offset], lde, (ftnlen)1);
#line 986 "BD01AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
#line 987 "BD01AD.f"
	    temp = dpar[3] * -2.;
#line 988 "BD01AD.f"
	    i__1 = l;
#line 988 "BD01AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 989 "BD01AD.f"
		e[i__ + i__ * e_dim1] = 1.;
#line 990 "BD01AD.f"
		a[i__ + (i__ + l) * a_dim1] = 1.;
#line 991 "BD01AD.f"
		a[i__ + l + (i__ + l) * a_dim1] = -dpar[2];
#line 992 "BD01AD.f"
		if (i__ < l) {
#line 993 "BD01AD.f"
		    a[i__ + l + (i__ + 1) * a_dim1] = dpar[3];
#line 994 "BD01AD.f"
		    a[i__ + l + 1 + i__ * a_dim1] = dpar[3];
#line 995 "BD01AD.f"
		    if (i__ > 1) {
#line 996 "BD01AD.f"
			a[i__ + l + i__ * a_dim1] = temp;
#line 997 "BD01AD.f"
		    }
#line 998 "BD01AD.f"
		}
#line 999 "BD01AD.f"
/* L430: */
#line 999 "BD01AD.f"
	    }
#line 1000 "BD01AD.f"
	    a[l + 1 + a_dim1] = -dpar[3];
#line 1001 "BD01AD.f"
	    a[*n + l * a_dim1] = -dpar[3];
#line 1002 "BD01AD.f"
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
#line 1003 "BD01AD.f"
	    b[l + 1 + b_dim1] = 1.;
#line 1004 "BD01AD.f"
	    b[*n + (b_dim1 << 1)] = -1.;
#line 1005 "BD01AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
#line 1006 "BD01AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 1008 "BD01AD.f"
	} else {
#line 1009 "BD01AD.f"
	    *info = -2;
#line 1010 "BD01AD.f"
	}
#line 1011 "BD01AD.f"
    } else {
#line 1012 "BD01AD.f"
	*info = -2;
#line 1013 "BD01AD.f"
    }

#line 1015 "BD01AD.f"
    return 0;
/* *** Last Line of BD01AD *** */
} /* bd01ad_ */

