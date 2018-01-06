#line 1 "BD02AD.f"
/* BD02AD.f -- translated by f2c (version 20100827).
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

#line 1 "BD02AD.f"
/* Table of constant values */

static doublereal c_b5 = 0.;
static doublereal c_b6 = 1.;
static doublereal c_b8 = -1.;
static integer c__5 = 5;
static integer c__1 = 1;

/* Subroutine */ int bd02ad_(char *def, integer *nr, doublereal *dpar, 
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
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_open(olist *), s_rsle(cilist *), do_lio(integer *, integer *, 
	    char *, ftnlen), e_rsle(void), f_clos(cllist *), s_wsfi(icilist *)
	    , do_fio(integer *, char *, ftnlen), e_wsfi(void);

    /* Local variables */
    static integer i__, j;
    static doublereal temp;
    static char dataf[12];
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    static integer status;

    /* Fortran I/O blocks */
    static cilist io___4 = { 1, 1, 1, 0, 0 };
    static icilist io___7 = { 0, dataf, 0, "(A,I2.2,A)", 11, 1 };
    static cilist io___8 = { 1, 1, 1, 0, 0 };
    static cilist io___9 = { 1, 1, 1, 0, 0 };



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
/*     discrete-time dynamical systems */

/*       E x_k+1 = A x_k + B u_k */

/*           y_k = C x_k + D u_k */

/*     E, A are real N-by-N matrices, B is N-by-M, C is P-by-N, and */
/*     D is P-by-M. In many examples, E is the identity matrix and D is */
/*     the zero matrix. */

/*     This routine is an implementation of the benchmark library */
/*     DTDSX (Version 1.0) described in [1]. */

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
/*             For Example 2.1, DPAR(1), ..., DPAR(3) define the */
/*             parameters 'tau', 'delta', 'K', respectively. */
/*             On exit, if DEF = 'D' and the desired example depends on */
/*             real parameters, then the array DPAR is overwritten by the */
/*             default values given in [1]. */

/*     IPAR    (input/output) INTEGER array, dimension (1) */
/*             On entry, if DEF = 'N' and the desired example depends on */
/*             integer parameters, then the array IPAR must contain the */
/*             values for these parameters. */
/*             For an explanation of the parameters see [1]. */
/*             For Example 3.1, IPAR(1) defines the parameter 'n'. */
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
/*             NOTE that DWORK is not used in the current version */
/*             of BD02AD. */

/*     LDWORK  INTEGER */
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
/*          DTDSX - a Collection of Benchmark Examples for State-Space */
/*          Realizations of Discrete-Time Dynamical Systems. */
/*          SLICOT Working Note 1998-10. 1998. */

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

/*     discrete-time dynamical systems */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     . LAPACK . */
/*     .. External Subroutines .. */
/*     . LAPACK . */
/*     .. Data Statements .. */
/*     . default values for availabities . */
#line 228 "BD02AD.f"
    /* Parameter adjustments */
#line 228 "BD02AD.f"
    --nr;
#line 228 "BD02AD.f"
    --dpar;
#line 228 "BD02AD.f"
    --ipar;
#line 228 "BD02AD.f"
    --vec;
#line 228 "BD02AD.f"
    e_dim1 = *lde;
#line 228 "BD02AD.f"
    e_offset = 1 + e_dim1;
#line 228 "BD02AD.f"
    e -= e_offset;
#line 228 "BD02AD.f"
    a_dim1 = *lda;
#line 228 "BD02AD.f"
    a_offset = 1 + a_dim1;
#line 228 "BD02AD.f"
    a -= a_offset;
#line 228 "BD02AD.f"
    b_dim1 = *ldb;
#line 228 "BD02AD.f"
    b_offset = 1 + b_dim1;
#line 228 "BD02AD.f"
    b -= b_offset;
#line 228 "BD02AD.f"
    c_dim1 = *ldc;
#line 228 "BD02AD.f"
    c_offset = 1 + c_dim1;
#line 228 "BD02AD.f"
    c__ -= c_offset;
#line 228 "BD02AD.f"
    d_dim1 = *ldd;
#line 228 "BD02AD.f"
    d_offset = 1 + d_dim1;
#line 228 "BD02AD.f"
    d__ -= d_offset;
#line 228 "BD02AD.f"
    --dwork;
#line 228 "BD02AD.f"

#line 228 "BD02AD.f"
    /* Function Body */

/*     .. Executable Statements .. */

#line 233 "BD02AD.f"
    *info = 0;
#line 234 "BD02AD.f"
    for (i__ = 1; i__ <= 8; ++i__) {
#line 235 "BD02AD.f"
	vec[i__] = vecdef[i__ - 1];
#line 236 "BD02AD.f"
/* L10: */
#line 236 "BD02AD.f"
    }

#line 238 "BD02AD.f"
    if (nr[1] == 1) {

#line 240 "BD02AD.f"
	if (nr[2] == 1) {
#line 241 "BD02AD.f"
	    s_copy(note, "Laub 1979, Ex. 2: uncontrollable-unobservable data",
		     (ftnlen)70, (ftnlen)50);
#line 242 "BD02AD.f"
	    *n = 2;
#line 243 "BD02AD.f"
	    *m = 1;
#line 244 "BD02AD.f"
	    *p = 1;
#line 245 "BD02AD.f"
	    if (*lde < *n) {
#line 245 "BD02AD.f"
		*info = -10;
#line 245 "BD02AD.f"
	    }
#line 246 "BD02AD.f"
	    if (*lda < *n) {
#line 246 "BD02AD.f"
		*info = -12;
#line 246 "BD02AD.f"
	    }
#line 247 "BD02AD.f"
	    if (*ldb < *n) {
#line 247 "BD02AD.f"
		*info = -14;
#line 247 "BD02AD.f"
	    }
#line 248 "BD02AD.f"
	    if (*ldc < *p) {
#line 248 "BD02AD.f"
		*info = -16;
#line 248 "BD02AD.f"
	    }
#line 249 "BD02AD.f"
	    if (*ldd < *p) {
#line 249 "BD02AD.f"
		*info = -18;
#line 249 "BD02AD.f"
	    }
#line 250 "BD02AD.f"
	    if (*info != 0) {
#line 250 "BD02AD.f"
		return 0;
#line 250 "BD02AD.f"
	    }

#line 252 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 253 "BD02AD.f"
	    a[a_dim1 + 1] = 4.;
#line 254 "BD02AD.f"
	    a[a_dim1 + 2] = -4.5;
#line 255 "BD02AD.f"
	    a[(a_dim1 << 1) + 1] = 3.;
#line 256 "BD02AD.f"
	    a[(a_dim1 << 1) + 2] = -3.5;
#line 257 "BD02AD.f"
	    dlaset_("A", n, m, &c_b8, &c_b6, &b[b_offset], ldb, (ftnlen)1);
#line 258 "BD02AD.f"
	    c__[c_dim1 + 1] = 3.;
#line 259 "BD02AD.f"
	    c__[(c_dim1 << 1) + 1] = 2.;
#line 260 "BD02AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 262 "BD02AD.f"
	} else if (nr[2] == 2) {
#line 263 "BD02AD.f"
	    s_copy(note, "Laub 1979, Ex. 3", (ftnlen)70, (ftnlen)16);
#line 264 "BD02AD.f"
	    *n = 2;
#line 265 "BD02AD.f"
	    *m = 2;
#line 266 "BD02AD.f"
	    *p = 2;
#line 267 "BD02AD.f"
	    if (*lde < *n) {
#line 267 "BD02AD.f"
		*info = -10;
#line 267 "BD02AD.f"
	    }
#line 268 "BD02AD.f"
	    if (*lda < *n) {
#line 268 "BD02AD.f"
		*info = -12;
#line 268 "BD02AD.f"
	    }
#line 269 "BD02AD.f"
	    if (*ldb < *n) {
#line 269 "BD02AD.f"
		*info = -14;
#line 269 "BD02AD.f"
	    }
#line 270 "BD02AD.f"
	    if (*ldc < *p) {
#line 270 "BD02AD.f"
		*info = -16;
#line 270 "BD02AD.f"
	    }
#line 271 "BD02AD.f"
	    if (*ldd < *p) {
#line 271 "BD02AD.f"
		*info = -18;
#line 271 "BD02AD.f"
	    }
#line 272 "BD02AD.f"
	    if (*info != 0) {
#line 272 "BD02AD.f"
		return 0;
#line 272 "BD02AD.f"
	    }

#line 274 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 275 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
#line 276 "BD02AD.f"
	    a[a_dim1 + 1] = .9512;
#line 277 "BD02AD.f"
	    a[(a_dim1 << 1) + 2] = .9048;
#line 278 "BD02AD.f"
	    b[b_dim1 + 1] = 4.877;
#line 279 "BD02AD.f"
	    b[(b_dim1 << 1) + 1] = 4.877;
#line 280 "BD02AD.f"
	    b[b_dim1 + 2] = -1.1895;
#line 281 "BD02AD.f"
	    b[(b_dim1 << 1) + 2] = 3.569;
#line 282 "BD02AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
#line 283 "BD02AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 285 "BD02AD.f"
	} else if (nr[2] == 3) {
#line 286 "BD02AD.f"
	    s_copy(note, "Van Dooren 1981, Ex. II", (ftnlen)70, (ftnlen)23);
#line 287 "BD02AD.f"
	    *n = 2;
#line 288 "BD02AD.f"
	    *m = 1;
#line 289 "BD02AD.f"
	    *p = 1;
#line 290 "BD02AD.f"
	    if (*lde < *n) {
#line 290 "BD02AD.f"
		*info = -10;
#line 290 "BD02AD.f"
	    }
#line 291 "BD02AD.f"
	    if (*lda < *n) {
#line 291 "BD02AD.f"
		*info = -12;
#line 291 "BD02AD.f"
	    }
#line 292 "BD02AD.f"
	    if (*ldb < *n) {
#line 292 "BD02AD.f"
		*info = -14;
#line 292 "BD02AD.f"
	    }
#line 293 "BD02AD.f"
	    if (*ldc < *p) {
#line 293 "BD02AD.f"
		*info = -16;
#line 293 "BD02AD.f"
	    }
#line 294 "BD02AD.f"
	    if (*ldd < *p) {
#line 294 "BD02AD.f"
		*info = -18;
#line 294 "BD02AD.f"
	    }
#line 295 "BD02AD.f"
	    if (*info != 0) {
#line 295 "BD02AD.f"
		return 0;
#line 295 "BD02AD.f"
	    }

#line 297 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 298 "BD02AD.f"
	    a[a_dim1 + 1] = 2.;
#line 299 "BD02AD.f"
	    a[a_dim1 + 2] = 1.;
#line 300 "BD02AD.f"
	    a[(a_dim1 << 1) + 1] = -1.;
#line 301 "BD02AD.f"
	    a[(a_dim1 << 1) + 2] = 0.;
#line 302 "BD02AD.f"
	    dlaset_("A", n, m, &c_b5, &c_b6, &b[b_offset], ldb, (ftnlen)1);
#line 303 "BD02AD.f"
	    dlaset_("A", p, n, &c_b6, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
#line 304 "BD02AD.f"
	    d__[d_dim1 + 1] = 0.;

#line 306 "BD02AD.f"
	} else if (nr[2] == 4) {
#line 307 "BD02AD.f"
	    s_copy(note, "Ionescu/Weiss 1992", (ftnlen)70, (ftnlen)18);
#line 308 "BD02AD.f"
	    *n = 2;
#line 309 "BD02AD.f"
	    *m = 2;
#line 310 "BD02AD.f"
	    *p = 2;
#line 311 "BD02AD.f"
	    if (*lde < *n) {
#line 311 "BD02AD.f"
		*info = -10;
#line 311 "BD02AD.f"
	    }
#line 312 "BD02AD.f"
	    if (*lda < *n) {
#line 312 "BD02AD.f"
		*info = -12;
#line 312 "BD02AD.f"
	    }
#line 313 "BD02AD.f"
	    if (*ldb < *n) {
#line 313 "BD02AD.f"
		*info = -14;
#line 313 "BD02AD.f"
	    }
#line 314 "BD02AD.f"
	    if (*ldc < *p) {
#line 314 "BD02AD.f"
		*info = -16;
#line 314 "BD02AD.f"
	    }
#line 315 "BD02AD.f"
	    if (*ldd < *p) {
#line 315 "BD02AD.f"
		*info = -18;
#line 315 "BD02AD.f"
	    }
#line 316 "BD02AD.f"
	    if (*info != 0) {
#line 316 "BD02AD.f"
		return 0;
#line 316 "BD02AD.f"
	    }

#line 318 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 319 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
#line 320 "BD02AD.f"
	    a[(a_dim1 << 1) + 1] = 1.;
#line 321 "BD02AD.f"
	    a[(a_dim1 << 1) + 2] = -1.;
#line 322 "BD02AD.f"
	    dlaset_("A", n, m, &c_b5, &c_b6, &b[b_offset], ldb, (ftnlen)1);
#line 323 "BD02AD.f"
	    b[b_dim1 + 2] = 2.;
#line 324 "BD02AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
#line 325 "BD02AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 327 "BD02AD.f"
	} else if (nr[2] == 5) {
#line 328 "BD02AD.f"
	    s_copy(note, "Jonckheere 1981", (ftnlen)70, (ftnlen)15);
#line 329 "BD02AD.f"
	    *n = 2;
#line 330 "BD02AD.f"
	    *m = 1;
#line 331 "BD02AD.f"
	    *p = 2;
#line 332 "BD02AD.f"
	    if (*lde < *n) {
#line 332 "BD02AD.f"
		*info = -10;
#line 332 "BD02AD.f"
	    }
#line 333 "BD02AD.f"
	    if (*lda < *n) {
#line 333 "BD02AD.f"
		*info = -12;
#line 333 "BD02AD.f"
	    }
#line 334 "BD02AD.f"
	    if (*ldb < *n) {
#line 334 "BD02AD.f"
		*info = -14;
#line 334 "BD02AD.f"
	    }
#line 335 "BD02AD.f"
	    if (*ldc < *p) {
#line 335 "BD02AD.f"
		*info = -16;
#line 335 "BD02AD.f"
	    }
#line 336 "BD02AD.f"
	    if (*ldd < *p) {
#line 336 "BD02AD.f"
		*info = -18;
#line 336 "BD02AD.f"
	    }
#line 337 "BD02AD.f"
	    if (*info != 0) {
#line 337 "BD02AD.f"
		return 0;
#line 337 "BD02AD.f"
	    }

#line 339 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 340 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
#line 341 "BD02AD.f"
	    a[(a_dim1 << 1) + 1] = 1.;
#line 342 "BD02AD.f"
	    dlaset_("A", n, m, &c_b6, &c_b5, &b[b_offset], ldb, (ftnlen)1);
#line 343 "BD02AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
#line 344 "BD02AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 346 "BD02AD.f"
	} else if (nr[2] == 6) {
#line 347 "BD02AD.f"
	    s_copy(note, "Ackerson/Fu 1970: satellite control problem", (
		    ftnlen)70, (ftnlen)43);
#line 348 "BD02AD.f"
	    *n = 4;
#line 349 "BD02AD.f"
	    *m = 2;
#line 350 "BD02AD.f"
	    *p = 4;
#line 351 "BD02AD.f"
	    if (*lde < *n) {
#line 351 "BD02AD.f"
		*info = -10;
#line 351 "BD02AD.f"
	    }
#line 352 "BD02AD.f"
	    if (*lda < *n) {
#line 352 "BD02AD.f"
		*info = -12;
#line 352 "BD02AD.f"
	    }
#line 353 "BD02AD.f"
	    if (*ldb < *n) {
#line 353 "BD02AD.f"
		*info = -14;
#line 353 "BD02AD.f"
	    }
#line 354 "BD02AD.f"
	    if (*ldc < *p) {
#line 354 "BD02AD.f"
		*info = -16;
#line 354 "BD02AD.f"
	    }
#line 355 "BD02AD.f"
	    if (*ldd < *p) {
#line 355 "BD02AD.f"
		*info = -18;
#line 355 "BD02AD.f"
	    }
#line 356 "BD02AD.f"
	    if (*info != 0) {
#line 356 "BD02AD.f"
		return 0;
#line 356 "BD02AD.f"
	    }

#line 358 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 359 "BD02AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
#line 360 "BD02AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 362 "BD02AD.f"
	} else if (nr[2] == 7) {
#line 363 "BD02AD.f"
	    s_copy(note, "Litkouhi 1983: system with slow and fast modes", (
		    ftnlen)70, (ftnlen)46);
#line 364 "BD02AD.f"
	    *n = 4;
#line 365 "BD02AD.f"
	    *m = 2;
#line 366 "BD02AD.f"
	    *p = 4;
#line 367 "BD02AD.f"
	    if (*lde < *n) {
#line 367 "BD02AD.f"
		*info = -10;
#line 367 "BD02AD.f"
	    }
#line 368 "BD02AD.f"
	    if (*lda < *n) {
#line 368 "BD02AD.f"
		*info = -12;
#line 368 "BD02AD.f"
	    }
#line 369 "BD02AD.f"
	    if (*ldb < *n) {
#line 369 "BD02AD.f"
		*info = -14;
#line 369 "BD02AD.f"
	    }
#line 370 "BD02AD.f"
	    if (*ldc < *p) {
#line 370 "BD02AD.f"
		*info = -16;
#line 370 "BD02AD.f"
	    }
#line 371 "BD02AD.f"
	    if (*ldd < *p) {
#line 371 "BD02AD.f"
		*info = -18;
#line 371 "BD02AD.f"
	    }
#line 372 "BD02AD.f"
	    if (*info != 0) {
#line 372 "BD02AD.f"
		return 0;
#line 372 "BD02AD.f"
	    }

#line 374 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 375 "BD02AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
#line 376 "BD02AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 378 "BD02AD.f"
	} else if (nr[2] == 8) {
#line 379 "BD02AD.f"
	    s_copy(note, "Lu/Lin 1993, Ex. 4.3", (ftnlen)70, (ftnlen)20);
#line 380 "BD02AD.f"
	    *n = 4;
#line 381 "BD02AD.f"
	    *m = 4;
#line 382 "BD02AD.f"
	    *p = 4;
#line 383 "BD02AD.f"
	    if (*lde < *n) {
#line 383 "BD02AD.f"
		*info = -10;
#line 383 "BD02AD.f"
	    }
#line 384 "BD02AD.f"
	    if (*lda < *n) {
#line 384 "BD02AD.f"
		*info = -12;
#line 384 "BD02AD.f"
	    }
#line 385 "BD02AD.f"
	    if (*ldb < *n) {
#line 385 "BD02AD.f"
		*info = -14;
#line 385 "BD02AD.f"
	    }
#line 386 "BD02AD.f"
	    if (*ldc < *p) {
#line 386 "BD02AD.f"
		*info = -16;
#line 386 "BD02AD.f"
	    }
#line 387 "BD02AD.f"
	    if (*ldd < *p) {
#line 387 "BD02AD.f"
		*info = -18;
#line 387 "BD02AD.f"
	    }
#line 388 "BD02AD.f"
	    if (*info != 0) {
#line 388 "BD02AD.f"
		return 0;
#line 388 "BD02AD.f"
	    }

#line 390 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 391 "BD02AD.f"
	    dlaset_("U", p, n, &c_b6, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
#line 392 "BD02AD.f"
	    c__[c_dim1 * 3 + 1] = 2.;
#line 393 "BD02AD.f"
	    c__[(c_dim1 << 2) + 1] = 4.;
#line 394 "BD02AD.f"
	    c__[(c_dim1 << 2) + 2] = 2.;
#line 395 "BD02AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 397 "BD02AD.f"
	} else if (nr[2] == 9) {
#line 398 "BD02AD.f"
	    s_copy(note, "Gajic/Shen 1993, Section 2.7.4: chemical plant", (
		    ftnlen)70, (ftnlen)46);
#line 399 "BD02AD.f"
	    *n = 5;
#line 400 "BD02AD.f"
	    *m = 2;
#line 401 "BD02AD.f"
	    *p = 5;
#line 402 "BD02AD.f"
	    if (*lde < *n) {
#line 402 "BD02AD.f"
		*info = -10;
#line 402 "BD02AD.f"
	    }
#line 403 "BD02AD.f"
	    if (*lda < *n) {
#line 403 "BD02AD.f"
		*info = -12;
#line 403 "BD02AD.f"
	    }
#line 404 "BD02AD.f"
	    if (*ldb < *n) {
#line 404 "BD02AD.f"
		*info = -14;
#line 404 "BD02AD.f"
	    }
#line 405 "BD02AD.f"
	    if (*ldc < *p) {
#line 405 "BD02AD.f"
		*info = -16;
#line 405 "BD02AD.f"
	    }
#line 406 "BD02AD.f"
	    if (*ldd < *p) {
#line 406 "BD02AD.f"
		*info = -18;
#line 406 "BD02AD.f"
	    }
#line 407 "BD02AD.f"
	    if (*info != 0) {
#line 407 "BD02AD.f"
		return 0;
#line 407 "BD02AD.f"
	    }

#line 409 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 410 "BD02AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
#line 411 "BD02AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 413 "BD02AD.f"
	} else if (nr[2] == 10) {
#line 414 "BD02AD.f"
	    s_copy(note, "Davison/Wang 1974", (ftnlen)70, (ftnlen)17);
#line 415 "BD02AD.f"
	    *n = 6;
#line 416 "BD02AD.f"
	    *m = 2;
#line 417 "BD02AD.f"
	    *p = 2;
#line 418 "BD02AD.f"
	    if (*lde < *n) {
#line 418 "BD02AD.f"
		*info = -10;
#line 418 "BD02AD.f"
	    }
#line 419 "BD02AD.f"
	    if (*lda < *n) {
#line 419 "BD02AD.f"
		*info = -12;
#line 419 "BD02AD.f"
	    }
#line 420 "BD02AD.f"
	    if (*ldb < *n) {
#line 420 "BD02AD.f"
		*info = -14;
#line 420 "BD02AD.f"
	    }
#line 421 "BD02AD.f"
	    if (*ldc < *p) {
#line 421 "BD02AD.f"
		*info = -16;
#line 421 "BD02AD.f"
	    }
#line 422 "BD02AD.f"
	    if (*ldd < *p) {
#line 422 "BD02AD.f"
		*info = -18;
#line 422 "BD02AD.f"
	    }
#line 423 "BD02AD.f"
	    if (*info != 0) {
#line 423 "BD02AD.f"
		return 0;
#line 423 "BD02AD.f"
	    }
#line 424 "BD02AD.f"
	    vec[8] = TRUE_;

#line 426 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 427 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
#line 428 "BD02AD.f"
	    a[(a_dim1 << 1) + 1] = 1.;
#line 429 "BD02AD.f"
	    a[a_dim1 * 3 + 2] = 1.;
#line 430 "BD02AD.f"
	    a[a_dim1 * 5 + 4] = 1.;
#line 431 "BD02AD.f"
	    a[a_dim1 * 6 + 5] = 1.;
#line 432 "BD02AD.f"
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
#line 433 "BD02AD.f"
	    b[b_dim1 + 3] = 1.;
#line 434 "BD02AD.f"
	    b[(b_dim1 << 1) + 6] = 1.;
#line 435 "BD02AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
#line 436 "BD02AD.f"
	    c__[c_dim1 + 1] = 1.;
#line 437 "BD02AD.f"
	    c__[(c_dim1 << 1) + 1] = 1.;
#line 438 "BD02AD.f"
	    c__[(c_dim1 << 2) + 2] = 1.;
#line 439 "BD02AD.f"
	    c__[c_dim1 * 5 + 2] = -1.;
#line 440 "BD02AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);
#line 441 "BD02AD.f"
	    d__[d_dim1 + 1] = 1.;
#line 442 "BD02AD.f"
	    d__[d_dim1 + 2] = 1.;

#line 444 "BD02AD.f"
	} else if (nr[2] == 11) {
#line 445 "BD02AD.f"
	    s_copy(note, "Patnaik et al. 1980: tubular ammonia reactor", (
		    ftnlen)70, (ftnlen)44);
#line 446 "BD02AD.f"
	    *n = 9;
#line 447 "BD02AD.f"
	    *m = 3;
#line 448 "BD02AD.f"
	    *p = 2;
#line 449 "BD02AD.f"
	    if (*lde < *n) {
#line 449 "BD02AD.f"
		*info = -10;
#line 449 "BD02AD.f"
	    }
#line 450 "BD02AD.f"
	    if (*lda < *n) {
#line 450 "BD02AD.f"
		*info = -12;
#line 450 "BD02AD.f"
	    }
#line 451 "BD02AD.f"
	    if (*ldb < *n) {
#line 451 "BD02AD.f"
		*info = -14;
#line 451 "BD02AD.f"
	    }
#line 452 "BD02AD.f"
	    if (*ldc < *p) {
#line 452 "BD02AD.f"
		*info = -16;
#line 452 "BD02AD.f"
	    }
#line 453 "BD02AD.f"
	    if (*ldd < *p) {
#line 453 "BD02AD.f"
		*info = -18;
#line 453 "BD02AD.f"
	    }
#line 454 "BD02AD.f"
	    if (*info != 0) {
#line 454 "BD02AD.f"
		return 0;
#line 454 "BD02AD.f"
	    }

#line 456 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 457 "BD02AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
#line 458 "BD02AD.f"
	    c__[c_dim1 + 1] = 1.;
#line 459 "BD02AD.f"
	    c__[c_dim1 * 5 + 2] = 1.;
#line 460 "BD02AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 462 "BD02AD.f"
	} else if (nr[2] == 12) {
#line 463 "BD02AD.f"
	    s_copy(note, "Smith 1969: two-stand cold rolling mill", (ftnlen)
		    70, (ftnlen)39);
#line 464 "BD02AD.f"
	    *n = 10;
#line 465 "BD02AD.f"
	    *m = 3;
#line 466 "BD02AD.f"
	    *p = 5;
#line 467 "BD02AD.f"
	    if (*lde < *n) {
#line 467 "BD02AD.f"
		*info = -10;
#line 467 "BD02AD.f"
	    }
#line 468 "BD02AD.f"
	    if (*lda < *n) {
#line 468 "BD02AD.f"
		*info = -12;
#line 468 "BD02AD.f"
	    }
#line 469 "BD02AD.f"
	    if (*ldb < *n) {
#line 469 "BD02AD.f"
		*info = -14;
#line 469 "BD02AD.f"
	    }
#line 470 "BD02AD.f"
	    if (*ldc < *p) {
#line 470 "BD02AD.f"
		*info = -16;
#line 470 "BD02AD.f"
	    }
#line 471 "BD02AD.f"
	    if (*ldd < *p) {
#line 471 "BD02AD.f"
		*info = -18;
#line 471 "BD02AD.f"
	    }
#line 472 "BD02AD.f"
	    if (*info != 0) {
#line 472 "BD02AD.f"
		return 0;
#line 472 "BD02AD.f"
	    }
#line 473 "BD02AD.f"
	    vec[8] = TRUE_;

#line 475 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 476 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
#line 477 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &a[a_dim1 + 2], lda, (ftnlen)1);
#line 478 "BD02AD.f"
	    a[a_dim1 * 10 + 1] = .112;
#line 479 "BD02AD.f"
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
#line 480 "BD02AD.f"
	    b[b_dim1 + 1] = 2.76;
#line 481 "BD02AD.f"
	    b[(b_dim1 << 1) + 1] = -1.35;
#line 482 "BD02AD.f"
	    b[b_dim1 * 3 + 1] = -.46;
#line 483 "BD02AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
#line 484 "BD02AD.f"
	    c__[c_dim1 + 1] = 1.;
#line 485 "BD02AD.f"
	    c__[c_dim1 * 10 + 2] = .894;
#line 486 "BD02AD.f"
	    c__[c_dim1 * 10 + 3] = -16.93;
#line 487 "BD02AD.f"
	    c__[c_dim1 * 10 + 4] = .07;
#line 488 "BD02AD.f"
	    c__[c_dim1 * 10 + 5] = .398;
#line 489 "BD02AD.f"
	    o__1.oerr = 1;
#line 489 "BD02AD.f"
	    o__1.ounit = 1;
#line 489 "BD02AD.f"
	    o__1.ofnmlen = 11;
#line 489 "BD02AD.f"
	    o__1.ofnm = "BD02112.dat";
#line 489 "BD02AD.f"
	    o__1.orl = 0;
#line 489 "BD02AD.f"
	    o__1.osta = "OLD";
#line 489 "BD02AD.f"
	    o__1.oacc = 0;
#line 489 "BD02AD.f"
	    o__1.ofm = 0;
#line 489 "BD02AD.f"
	    o__1.oblnk = 0;
#line 489 "BD02AD.f"
	    status = f_open(&o__1);
#line 490 "BD02AD.f"
	    if (status != 0) {
#line 491 "BD02AD.f"
		*info = 1;
#line 492 "BD02AD.f"
	    } else {
#line 493 "BD02AD.f"
		i__1 = *p;
#line 493 "BD02AD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 494 "BD02AD.f"
		    status = s_rsle(&io___4);
#line 494 "BD02AD.f"
		    if (status != 0) {
#line 494 "BD02AD.f"
			goto L100001;
#line 494 "BD02AD.f"
		    }
#line 494 "BD02AD.f"
		    i__2 = *m;
#line 494 "BD02AD.f"
		    for (j = 1; j <= i__2; ++j) {
#line 494 "BD02AD.f"
			status = do_lio(&c__5, &c__1, (char *)&d__[i__ + j * 
				d_dim1], (ftnlen)sizeof(doublereal));
#line 494 "BD02AD.f"
			if (status != 0) {
#line 494 "BD02AD.f"
			    goto L100001;
#line 494 "BD02AD.f"
			}
#line 494 "BD02AD.f"
		    }
#line 494 "BD02AD.f"
		    status = e_rsle();
#line 494 "BD02AD.f"
L100001:
#line 495 "BD02AD.f"
		    if (status != 0) {
#line 495 "BD02AD.f"
			*info = 1;
#line 495 "BD02AD.f"
		    }
#line 496 "BD02AD.f"
/* L110: */
#line 496 "BD02AD.f"
		}
#line 497 "BD02AD.f"
	    }
#line 498 "BD02AD.f"
	    cl__1.cerr = 0;
#line 498 "BD02AD.f"
	    cl__1.cunit = 1;
#line 498 "BD02AD.f"
	    cl__1.csta = 0;
#line 498 "BD02AD.f"
	    f_clos(&cl__1);

#line 500 "BD02AD.f"
	} else {
#line 501 "BD02AD.f"
	    *info = -2;
#line 502 "BD02AD.f"
	}

#line 504 "BD02AD.f"
	if (nr[2] >= 6 && nr[2] <= 9 || nr[2] == 11) {
/*         .. loading data files */
#line 507 "BD02AD.f"
	    s_wsfi(&io___7);
#line 507 "BD02AD.f"
	    do_fio(&c__1, "BD021", (ftnlen)5);
#line 507 "BD02AD.f"
	    do_fio(&c__1, (char *)&nr[2], (ftnlen)sizeof(integer));
#line 507 "BD02AD.f"
	    do_fio(&c__1, ".dat", (ftnlen)4);
#line 507 "BD02AD.f"
	    e_wsfi();
#line 508 "BD02AD.f"
	    o__1.oerr = 1;
#line 508 "BD02AD.f"
	    o__1.ounit = 1;
#line 508 "BD02AD.f"
	    o__1.ofnmlen = 11;
#line 508 "BD02AD.f"
	    o__1.ofnm = dataf;
#line 508 "BD02AD.f"
	    o__1.orl = 0;
#line 508 "BD02AD.f"
	    o__1.osta = "OLD";
#line 508 "BD02AD.f"
	    o__1.oacc = 0;
#line 508 "BD02AD.f"
	    o__1.ofm = 0;
#line 508 "BD02AD.f"
	    o__1.oblnk = 0;
#line 508 "BD02AD.f"
	    status = f_open(&o__1);
#line 509 "BD02AD.f"
	    if (status != 0) {
#line 510 "BD02AD.f"
		*info = 1;
#line 511 "BD02AD.f"
	    } else {
#line 512 "BD02AD.f"
		i__1 = *n;
#line 512 "BD02AD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 513 "BD02AD.f"
		    status = s_rsle(&io___8);
#line 513 "BD02AD.f"
		    if (status != 0) {
#line 513 "BD02AD.f"
			goto L100002;
#line 513 "BD02AD.f"
		    }
#line 513 "BD02AD.f"
		    i__2 = *n;
#line 513 "BD02AD.f"
		    for (j = 1; j <= i__2; ++j) {
#line 513 "BD02AD.f"
			status = do_lio(&c__5, &c__1, (char *)&a[i__ + j * 
				a_dim1], (ftnlen)sizeof(doublereal));
#line 513 "BD02AD.f"
			if (status != 0) {
#line 513 "BD02AD.f"
			    goto L100002;
#line 513 "BD02AD.f"
			}
#line 513 "BD02AD.f"
		    }
#line 513 "BD02AD.f"
		    status = e_rsle();
#line 513 "BD02AD.f"
L100002:
#line 514 "BD02AD.f"
		    if (status != 0) {
#line 514 "BD02AD.f"
			*info = 1;
#line 514 "BD02AD.f"
		    }
#line 515 "BD02AD.f"
/* L120: */
#line 515 "BD02AD.f"
		}
#line 516 "BD02AD.f"
		i__1 = *n;
#line 516 "BD02AD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 517 "BD02AD.f"
		    status = s_rsle(&io___9);
#line 517 "BD02AD.f"
		    if (status != 0) {
#line 517 "BD02AD.f"
			goto L100003;
#line 517 "BD02AD.f"
		    }
#line 517 "BD02AD.f"
		    i__2 = *m;
#line 517 "BD02AD.f"
		    for (j = 1; j <= i__2; ++j) {
#line 517 "BD02AD.f"
			status = do_lio(&c__5, &c__1, (char *)&b[i__ + j * 
				b_dim1], (ftnlen)sizeof(doublereal));
#line 517 "BD02AD.f"
			if (status != 0) {
#line 517 "BD02AD.f"
			    goto L100003;
#line 517 "BD02AD.f"
			}
#line 517 "BD02AD.f"
		    }
#line 517 "BD02AD.f"
		    status = e_rsle();
#line 517 "BD02AD.f"
L100003:
#line 518 "BD02AD.f"
		    if (status != 0) {
#line 518 "BD02AD.f"
			*info = 1;
#line 518 "BD02AD.f"
		    }
#line 519 "BD02AD.f"
/* L130: */
#line 519 "BD02AD.f"
		}
#line 520 "BD02AD.f"
	    }
#line 521 "BD02AD.f"
	    cl__1.cerr = 0;
#line 521 "BD02AD.f"
	    cl__1.cunit = 1;
#line 521 "BD02AD.f"
	    cl__1.csta = 0;
#line 521 "BD02AD.f"
	    f_clos(&cl__1);
#line 522 "BD02AD.f"
	}

#line 524 "BD02AD.f"
    } else if (nr[1] == 2) {
#line 525 "BD02AD.f"
	if (! (lsame_(def, "D", (ftnlen)1, (ftnlen)1) || lsame_(def, "N", (
		ftnlen)1, (ftnlen)1))) {
#line 526 "BD02AD.f"
	    *info = -1;
#line 527 "BD02AD.f"
	    return 0;
#line 528 "BD02AD.f"
	}

#line 530 "BD02AD.f"
	if (nr[2] == 1) {
#line 531 "BD02AD.f"
	    s_copy(note, "Pappas et al. 1980: process control of paper machi"\
		    "ne", (ftnlen)70, (ftnlen)52);
#line 532 "BD02AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 533 "BD02AD.f"
		dpar[1] = 1e8;
#line 534 "BD02AD.f"
		dpar[2] = 1.;
#line 535 "BD02AD.f"
		dpar[3] = 1.;
#line 536 "BD02AD.f"
	    }
#line 537 "BD02AD.f"
	    if (dpar[1] == 0.) {
#line 537 "BD02AD.f"
		*info = -3;
#line 537 "BD02AD.f"
	    }
#line 538 "BD02AD.f"
	    *n = 4;
#line 539 "BD02AD.f"
	    *m = 1;
#line 540 "BD02AD.f"
	    *p = 1;
#line 541 "BD02AD.f"
	    if (*lde < *n) {
#line 541 "BD02AD.f"
		*info = -10;
#line 541 "BD02AD.f"
	    }
#line 542 "BD02AD.f"
	    if (*lda < *n) {
#line 542 "BD02AD.f"
		*info = -12;
#line 542 "BD02AD.f"
	    }
#line 543 "BD02AD.f"
	    if (*ldb < *n) {
#line 543 "BD02AD.f"
		*info = -14;
#line 543 "BD02AD.f"
	    }
#line 544 "BD02AD.f"
	    if (*ldc < *p) {
#line 544 "BD02AD.f"
		*info = -16;
#line 544 "BD02AD.f"
	    }
#line 545 "BD02AD.f"
	    if (*ldd < *p) {
#line 545 "BD02AD.f"
		*info = -18;
#line 545 "BD02AD.f"
	    }
#line 546 "BD02AD.f"
	    if (*info != 0) {
#line 546 "BD02AD.f"
		return 0;
#line 546 "BD02AD.f"
	    }

#line 548 "BD02AD.f"
	    temp = dpar[2] / dpar[1];
#line 549 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 550 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
#line 551 "BD02AD.f"
	    i__1 = *n - 1;
#line 551 "BD02AD.f"
	    i__2 = *n - 1;
#line 551 "BD02AD.f"
	    dlaset_("A", &i__1, &i__2, &c_b5, &c_b6, &a[a_dim1 + 2], lda, (
		    ftnlen)1);
#line 552 "BD02AD.f"
	    a[a_dim1 + 1] = 1. - temp;
#line 553 "BD02AD.f"
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
#line 554 "BD02AD.f"
	    b[b_dim1 + 1] = dpar[3] * temp;
#line 555 "BD02AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b5, &c__[c_offset], ldc, (ftnlen)1);
#line 556 "BD02AD.f"
	    c__[(c_dim1 << 2) + 1] = 1.;
#line 557 "BD02AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 559 "BD02AD.f"
	} else {
#line 560 "BD02AD.f"
	    *info = -2;
#line 561 "BD02AD.f"
	}

#line 563 "BD02AD.f"
    } else if (nr[1] == 3) {
#line 564 "BD02AD.f"
	if (! (lsame_(def, "D", (ftnlen)1, (ftnlen)1) || lsame_(def, "N", (
		ftnlen)1, (ftnlen)1))) {
#line 565 "BD02AD.f"
	    *info = -1;
#line 566 "BD02AD.f"
	    return 0;
#line 567 "BD02AD.f"
	}

#line 569 "BD02AD.f"
	if (nr[2] == 1) {
#line 570 "BD02AD.f"
	    s_copy(note, "Pappas et al. 1980, Ex. 3", (ftnlen)70, (ftnlen)25);
#line 571 "BD02AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 571 "BD02AD.f"
		ipar[1] = 100;
#line 571 "BD02AD.f"
	    }
#line 572 "BD02AD.f"
	    if (ipar[1] < 2) {
#line 572 "BD02AD.f"
		*info = -4;
#line 572 "BD02AD.f"
	    }
#line 573 "BD02AD.f"
	    *n = ipar[1];
#line 574 "BD02AD.f"
	    *m = 1;
#line 575 "BD02AD.f"
	    *p = *n;
#line 576 "BD02AD.f"
	    if (*lde < *n) {
#line 576 "BD02AD.f"
		*info = -10;
#line 576 "BD02AD.f"
	    }
#line 577 "BD02AD.f"
	    if (*lda < *n) {
#line 577 "BD02AD.f"
		*info = -12;
#line 577 "BD02AD.f"
	    }
#line 578 "BD02AD.f"
	    if (*ldb < *n) {
#line 578 "BD02AD.f"
		*info = -14;
#line 578 "BD02AD.f"
	    }
#line 579 "BD02AD.f"
	    if (*ldc < *p) {
#line 579 "BD02AD.f"
		*info = -16;
#line 579 "BD02AD.f"
	    }
#line 580 "BD02AD.f"
	    if (*ldd < *p) {
#line 580 "BD02AD.f"
		*info = -18;
#line 580 "BD02AD.f"
	    }
#line 581 "BD02AD.f"
	    if (*info != 0) {
#line 581 "BD02AD.f"
		return 0;
#line 581 "BD02AD.f"
	    }

#line 583 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b6, &e[e_offset], lde, (ftnlen)1);
#line 584 "BD02AD.f"
	    dlaset_("A", n, n, &c_b5, &c_b5, &a[a_offset], lda, (ftnlen)1);
#line 585 "BD02AD.f"
	    i__1 = *n - 1;
#line 585 "BD02AD.f"
	    i__2 = *n - 1;
#line 585 "BD02AD.f"
	    dlaset_("A", &i__1, &i__2, &c_b5, &c_b6, &a[(a_dim1 << 1) + 1], 
		    lda, (ftnlen)1);
#line 586 "BD02AD.f"
	    dlaset_("A", n, m, &c_b5, &c_b5, &b[b_offset], ldb, (ftnlen)1);
#line 587 "BD02AD.f"
	    b[*n + b_dim1] = 1.;
#line 588 "BD02AD.f"
	    dlaset_("A", p, n, &c_b5, &c_b6, &c__[c_offset], ldc, (ftnlen)1);
#line 589 "BD02AD.f"
	    dlaset_("A", p, m, &c_b5, &c_b5, &d__[d_offset], ldd, (ftnlen)1);

#line 591 "BD02AD.f"
	} else {
#line 592 "BD02AD.f"
	    *info = -2;
#line 593 "BD02AD.f"
	}

#line 595 "BD02AD.f"
    } else {
#line 596 "BD02AD.f"
	*info = -2;
#line 597 "BD02AD.f"
    }

#line 599 "BD02AD.f"
    return 0;
/* *** Last Line of BD02AD *** */
} /* bd02ad_ */

