#line 1 "MB02ID.f"
/* MB02ID.f -- translated by f2c (version 20100827).
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

#line 1 "MB02ID.f"
/* Table of constant values */

static doublereal c_b8 = 0.;
static doublereal c_b18 = 1.;
static doublereal c_b42 = -1.;
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;

/* Subroutine */ int mb02id_(char *job, integer *k, integer *l, integer *m, 
	integer *n, integer *rb, integer *rc, doublereal *tc, integer *ldtc, 
	doublereal *tr, integer *ldtr, doublereal *b, integer *ldb, 
	doublereal *c__, integer *ldc, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen job_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, c_dim1, c_offset, tc_dim1, tc_offset, tr_dim1, 
	    tr_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, 
	    i__10, i__11, i__12;

    /* Local variables */
    static integer i__, x, y, nb, kk, pt, pdi, len, pni, ppi, pdw, rnk, pnr, 
	    ppr, ierr, ipvt[1];
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    mb02kd_(char *, char *, integer *, integer *, integer *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    mb02cu_(char *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen), dgemm_(char *, char *
	    , integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen), mb02cv_(char *, char *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dgels_(char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin;
    static logical compo, compu;
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dtrsm_(
	    char *, char *, char *, char *, integer *, integer *, doublereal *
	    , doublereal *, integer *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen), dgeqrf_(integer *, integer *, doublereal 
	    *, integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dorgqr_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);
    static integer wrkmin;
    extern /* Subroutine */ int dtrtri_(char *, char *, integer *, doublereal 
	    *, integer *, integer *, ftnlen, ftnlen);
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

/*     To solve the overdetermined or underdetermined real linear systems */
/*     involving an M*K-by-N*L block Toeplitz matrix T that is specified */
/*     by its first block column and row. It is assumed that T has full */
/*     rank. */
/*     The following options are provided: */

/*     1. If JOB = 'O' or JOB = 'A' :  find the least squares solution of */
/*        an overdetermined system, i.e., solve the least squares problem */

/*                  minimize || B - T*X ||.                           (1) */

/*     2. If JOB = 'U' or JOB = 'A' :  find the minimum norm solution of */
/*        the undetermined system */
/*                   T */
/*                  T * X = C.                                        (2) */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the problem to be solved as follows */
/*             = 'O':  solve the overdetermined system (1); */
/*             = 'U':  solve the underdetermined system (2); */
/*             = 'A':  solve (1) and (2). */

/*     Input/Output Parameters */

/*     K       (input) INTEGER */
/*             The number of rows in the blocks of T.  K >= 0. */

/*     L       (input) INTEGER */
/*             The number of columns in the blocks of T.  L >= 0. */

/*     M       (input) INTEGER */
/*             The number of blocks in the first block column of T. */
/*             M >= 0. */

/*     N       (input) INTEGER */
/*             The number of blocks in the first block row of T. */
/*             0 <= N <= M*K / L. */

/*     RB      (input) INTEGER */
/*             If JOB = 'O' or 'A', the number of columns in B.  RB >= 0. */

/*     RC      (input) INTEGER */
/*             If JOB = 'U' or 'A', the number of columns in C.  RC >= 0. */

/*     TC      (input)  DOUBLE PRECISION array, dimension (LDTC,L) */
/*             On entry, the leading M*K-by-L part of this array must */
/*             contain the first block column of T. */

/*     LDTC    INTEGER */
/*             The leading dimension of the array TC.  LDTC >= MAX(1,M*K) */

/*     TR      (input)  DOUBLE PRECISION array, dimension (LDTR,(N-1)*L) */
/*             On entry, the leading K-by-(N-1)*L part of this array must */
/*             contain the 2nd to the N-th blocks of the first block row */
/*             of T. */

/*     LDTR    INTEGER */
/*             The leading dimension of the array TR.  LDTR >= MAX(1,K). */

/*     B       (input/output)  DOUBLE PRECISION array, dimension (LDB,RB) */
/*             On entry, if JOB = 'O' or JOB = 'A', the leading M*K-by-RB */
/*             part of this array must contain the right hand side */
/*             matrix B of the overdetermined system (1). */
/*             On exit, if JOB = 'O' or JOB = 'A', the leading N*L-by-RB */
/*             part of this array contains the solution of the */
/*             overdetermined system (1). */
/*             This array is not referenced if JOB = 'U'. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= MAX(1,M*K),  if JOB = 'O'  or  JOB = 'A'; */
/*             LDB >= 1,           if JOB = 'U'. */

/*     C       (input)  DOUBLE PRECISION array, dimension (LDC,RC) */
/*             On entry, if JOB = 'U' or JOB = 'A', the leading N*L-by-RC */
/*             part of this array must contain the right hand side */
/*             matrix C of the underdetermined system (2). */
/*             On exit, if JOB = 'U' or JOB = 'A', the leading M*K-by-RC */
/*             part of this array contains the solution of the */
/*             underdetermined system (2). */
/*             This array is not referenced if JOB = 'O'. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C. */
/*             LDB >= 1,           if JOB = 'O'; */
/*             LDB >= MAX(1,M*K),  if JOB = 'U'  or  JOB = 'A'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -17,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             Let x = MAX( 2*N*L*(L+K) + (6+N)*L,(N*L+M*K+1)*L + M*K ) */
/*             and y = N*M*K*L + N*L, then */
/*             if MIN( M,N ) = 1 and JOB = 'O', */
/*                         LDWORK >= MAX( y + MAX( M*K,RB ),1 ); */
/*             if MIN( M,N ) = 1 and JOB = 'U', */
/*                         LDWORK >= MAX( y + MAX( M*K,RC ),1 ); */
/*             if MIN( M,N ) = 1 and JOB = 'A', */
/*                         LDWORK >= MAX( y +MAX( M*K,MAX( RB,RC ),1 ); */
/*             if MIN( M,N ) > 1 and JOB = 'O', */
/*                         LDWORK >= MAX( x,N*L*RB + 1 ); */
/*             if MIN( M,N ) > 1 and JOB = 'U', */
/*                         LDWORK >= MAX( x,N*L*RC + 1 ); */
/*             if MIN( M,N ) > 1 and JOB = 'A', */
/*                         LDWORK >= MAX( x,N*L*MAX( RB,RC ) + 1 ). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the reduction algorithm failed. The Toeplitz matrix */
/*                   associated with T is (numerically) not of full rank. */

/*     METHOD */

/*     Householder transformations and modified hyperbolic rotations */
/*     are used in the Schur algorithm [1], [2]. */

/*     REFERENCES */

/*     [1] Kailath, T. and Sayed, A. */
/*         Fast Reliable Algorithms for Matrices with Structure. */
/*         SIAM Publications, Philadelphia, 1999. */

/*     [2] Kressner, D. and Van Dooren, P. */
/*         Factorizations and linear system solvers for matrices with */
/*         Toeplitz structure. */
/*         SLICOT Working Note 2000-2, 2000. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires O( L*L*K*(N+M)*log(N+M) + N*N*L*L*(L+K) ) */
/*     and additionally */

/*     if JOB = 'O' or JOB = 'A', */
/*                  O( (K*L+RB*L+K*RB)*(N+M)*log(N+M) + N*N*L*L*RB ); */
/*     if JOB = 'U' or JOB = 'A', */
/*                  O( (K*L+RC*L+K*RC)*(N+M)*log(N+M) + N*N*L*L*RC ); */

/*     floating point operations. */

/*     CONTRIBUTOR */

/*     D. Kressner, Technical Univ. Berlin, Germany, May 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, June 2001. */
/*     D. Kressner, Technical Univ. Berlin, Germany, July 2002. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2004. */

/*     KEYWORDS */

/*     Elementary matrix operations, Householder transformation, matrix */
/*     operations, Toeplitz matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

#line 226 "MB02ID.f"
    /* Parameter adjustments */
#line 226 "MB02ID.f"
    tc_dim1 = *ldtc;
#line 226 "MB02ID.f"
    tc_offset = 1 + tc_dim1;
#line 226 "MB02ID.f"
    tc -= tc_offset;
#line 226 "MB02ID.f"
    tr_dim1 = *ldtr;
#line 226 "MB02ID.f"
    tr_offset = 1 + tr_dim1;
#line 226 "MB02ID.f"
    tr -= tr_offset;
#line 226 "MB02ID.f"
    b_dim1 = *ldb;
#line 226 "MB02ID.f"
    b_offset = 1 + b_dim1;
#line 226 "MB02ID.f"
    b -= b_offset;
#line 226 "MB02ID.f"
    c_dim1 = *ldc;
#line 226 "MB02ID.f"
    c_offset = 1 + c_dim1;
#line 226 "MB02ID.f"
    c__ -= c_offset;
#line 226 "MB02ID.f"
    --dwork;
#line 226 "MB02ID.f"

#line 226 "MB02ID.f"
    /* Function Body */
#line 226 "MB02ID.f"
    *info = 0;
#line 227 "MB02ID.f"
    compo = lsame_(job, "O", (ftnlen)1, (ftnlen)1) || lsame_(job, "A", (
	    ftnlen)1, (ftnlen)1);
#line 228 "MB02ID.f"
    compu = lsame_(job, "U", (ftnlen)1, (ftnlen)1) || lsame_(job, "A", (
	    ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 229 "MB02ID.f"
    i__1 = (*n << 1) * *l * (*l + *k) + (*n + 6) * *l, i__2 = (*n * *l + *m * 
	    *k + 1) * *l + *m * *k;
#line 229 "MB02ID.f"
    x = max(i__1,i__2);
#line 231 "MB02ID.f"
    y = *n * *m * *k * *l + *n * *l;
#line 232 "MB02ID.f"
    if (min(*m,*n) == 1) {
/* Computing MAX */
#line 233 "MB02ID.f"
	i__1 = *m * *k;
#line 233 "MB02ID.f"
	wrkmin = max(i__1,1);
#line 234 "MB02ID.f"
	if (compo) {
#line 234 "MB02ID.f"
	    wrkmin = max(wrkmin,*rb);
#line 234 "MB02ID.f"
	}
#line 235 "MB02ID.f"
	if (compu) {
#line 235 "MB02ID.f"
	    wrkmin = max(wrkmin,*rc);
#line 235 "MB02ID.f"
	}
/* Computing MAX */
#line 236 "MB02ID.f"
	i__1 = y + wrkmin;
#line 236 "MB02ID.f"
	wrkmin = max(i__1,1);
#line 237 "MB02ID.f"
    } else {
#line 238 "MB02ID.f"
	wrkmin = x;
#line 239 "MB02ID.f"
	if (compo) {
/* Computing MAX */
#line 239 "MB02ID.f"
	    i__1 = wrkmin, i__2 = *n * *l * *rb + 1;
#line 239 "MB02ID.f"
	    wrkmin = max(i__1,i__2);
#line 239 "MB02ID.f"
	}
#line 240 "MB02ID.f"
	if (compu) {
/* Computing MAX */
#line 240 "MB02ID.f"
	    i__1 = wrkmin, i__2 = *n * *l * *rc + 1;
#line 240 "MB02ID.f"
	    wrkmin = max(i__1,i__2);
#line 240 "MB02ID.f"
	}
#line 241 "MB02ID.f"
    }
#line 242 "MB02ID.f"
    wrkopt = 1;

/*     Check the scalar input parameters. */

#line 246 "MB02ID.f"
    if (! (compo || compu)) {
#line 247 "MB02ID.f"
	*info = -1;
#line 248 "MB02ID.f"
    } else if (*k < 0) {
#line 249 "MB02ID.f"
	*info = -2;
#line 250 "MB02ID.f"
    } else if (*l < 0) {
#line 251 "MB02ID.f"
	*info = -3;
#line 252 "MB02ID.f"
    } else if (*m < 0) {
#line 253 "MB02ID.f"
	*info = -4;
#line 254 "MB02ID.f"
    } else if (*n < 0 || *n * *l > *m * *k) {
#line 255 "MB02ID.f"
	*info = -5;
#line 256 "MB02ID.f"
    } else if (compo && *rb < 0) {
#line 257 "MB02ID.f"
	*info = -6;
#line 258 "MB02ID.f"
    } else if (compu && *rc < 0) {
#line 259 "MB02ID.f"
	*info = -7;
#line 260 "MB02ID.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 260 "MB02ID.f"
	i__1 = 1, i__2 = *m * *k;
#line 260 "MB02ID.f"
	if (*ldtc < max(i__1,i__2)) {
#line 261 "MB02ID.f"
	    *info = -9;
#line 262 "MB02ID.f"
	} else if (*ldtr < max(1,*k)) {
#line 263 "MB02ID.f"
	    *info = -11;
#line 264 "MB02ID.f"
	} else if (*ldb < 1 || compo && *ldb < *m * *k) {
#line 265 "MB02ID.f"
	    *info = -13;
#line 266 "MB02ID.f"
	} else if (*ldc < 1 || compu && *ldc < *m * *k) {
#line 267 "MB02ID.f"
	    *info = -15;
#line 268 "MB02ID.f"
	} else if (*ldwork < wrkmin) {
#line 269 "MB02ID.f"
	    dwork[1] = (doublereal) wrkmin;
#line 270 "MB02ID.f"
	    *info = -17;
#line 271 "MB02ID.f"
	}
#line 271 "MB02ID.f"
    }

/*     Return if there were illegal values. */

#line 275 "MB02ID.f"
    if (*info != 0) {
#line 276 "MB02ID.f"
	i__1 = -(*info);
#line 276 "MB02ID.f"
	xerbla_("MB02ID", &i__1, (ftnlen)6);
#line 277 "MB02ID.f"
	return 0;
#line 278 "MB02ID.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 282 "MB02ID.f"
    i__1 = *n * *l;
#line 282 "MB02ID.f"
    if (compo && min(i__1,*rb) == 0) {
#line 283 "MB02ID.f"
	compo = FALSE_;
#line 284 "MB02ID.f"
    }
/* Computing MIN */
#line 285 "MB02ID.f"
    i__1 = *n * *l;
#line 285 "MB02ID.f"
    if (compu && min(i__1,*rc) == 0) {
#line 286 "MB02ID.f"
	i__1 = *m * *k;
#line 286 "MB02ID.f"
	dlaset_("Full", &i__1, rc, &c_b8, &c_b8, &c__[c_offset], ldc, (ftnlen)
		4);
#line 287 "MB02ID.f"
	compu = FALSE_;
#line 288 "MB02ID.f"
    }
#line 289 "MB02ID.f"
    if (! (compo || compu)) {
#line 290 "MB02ID.f"
	dwork[1] = 1.;
#line 291 "MB02ID.f"
	return 0;
#line 292 "MB02ID.f"
    }

/*     Check cases M = 1 or N = 1. */

#line 296 "MB02ID.f"
    if (min(*m,*n) == 1) {
#line 297 "MB02ID.f"
	pdw = *k * *l * *m * *n;
#line 298 "MB02ID.f"
	if (compo) {
#line 299 "MB02ID.f"
	    i__1 = *m * *k;
#line 299 "MB02ID.f"
	    i__2 = *m * *k;
#line 299 "MB02ID.f"
	    dlacpy_("All", &i__1, l, &tc[tc_offset], ldtc, &dwork[1], &i__2, (
		    ftnlen)3);
#line 300 "MB02ID.f"
	    i__1 = (*n - 1) * *l;
#line 300 "MB02ID.f"
	    i__2 = *m * *k;
#line 300 "MB02ID.f"
	    dlacpy_("All", k, &i__1, &tr[tr_offset], ldtr, &dwork[*k * *l + 1]
		    , &i__2, (ftnlen)3);
#line 302 "MB02ID.f"
	    i__1 = *m * *k;
#line 302 "MB02ID.f"
	    i__2 = *n * *l;
#line 302 "MB02ID.f"
	    i__3 = *m * *k;
#line 302 "MB02ID.f"
	    i__4 = *ldwork - pdw;
#line 302 "MB02ID.f"
	    dgels_("NonTranspose", &i__1, &i__2, rb, &dwork[1], &i__3, &b[
		    b_offset], ldb, &dwork[pdw + 1], &i__4, &ierr, (ftnlen)12)
		    ;
/* Computing MAX */
#line 304 "MB02ID.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[pdw + 1] + pdw;
#line 304 "MB02ID.f"
	    wrkopt = max(i__1,i__2);
#line 305 "MB02ID.f"
	}
#line 306 "MB02ID.f"
	if (compu) {
#line 307 "MB02ID.f"
	    i__1 = *m * *k;
#line 307 "MB02ID.f"
	    i__2 = *m * *k;
#line 307 "MB02ID.f"
	    dlacpy_("All", &i__1, l, &tc[tc_offset], ldtc, &dwork[1], &i__2, (
		    ftnlen)3);
#line 308 "MB02ID.f"
	    i__1 = (*n - 1) * *l;
#line 308 "MB02ID.f"
	    i__2 = *m * *k;
#line 308 "MB02ID.f"
	    dlacpy_("All", k, &i__1, &tr[tr_offset], ldtr, &dwork[*k * *l + 1]
		    , &i__2, (ftnlen)3);
#line 310 "MB02ID.f"
	    i__1 = *m * *k;
#line 310 "MB02ID.f"
	    i__2 = *n * *l;
#line 310 "MB02ID.f"
	    i__3 = *m * *k;
#line 310 "MB02ID.f"
	    i__4 = *ldwork - pdw;
#line 310 "MB02ID.f"
	    dgels_("Transpose", &i__1, &i__2, rc, &dwork[1], &i__3, &c__[
		    c_offset], ldc, &dwork[pdw + 1], &i__4, &ierr, (ftnlen)9);
/* Computing MAX */
#line 312 "MB02ID.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[pdw + 1] + pdw;
#line 312 "MB02ID.f"
	    wrkopt = max(i__1,i__2);
#line 313 "MB02ID.f"
	}
#line 314 "MB02ID.f"
	dwork[1] = (doublereal) wrkopt;
#line 315 "MB02ID.f"
	return 0;
#line 316 "MB02ID.f"
    }

/*     Step 1:  Compute the generator. */

#line 320 "MB02ID.f"
    if (compo) {
#line 321 "MB02ID.f"
	i__1 = *n * *l;
#line 321 "MB02ID.f"
	i__2 = *ldwork - *n * *l * *rb;
#line 321 "MB02ID.f"
	mb02kd_("Column", "Transpose", k, l, m, n, rb, &c_b18, &c_b8, &tc[
		tc_offset], ldtc, &tr[tr_offset], ldtr, &b[b_offset], ldb, &
		dwork[1], &i__1, &dwork[*n * *l * *rb + 1], &i__2, &ierr, (
		ftnlen)6, (ftnlen)9);
/* Computing MAX */
#line 324 "MB02ID.f"
	i__1 = wrkopt, i__2 = (integer) dwork[*n * *l * *rb + 1] + *n * *l * *
		rb;
#line 324 "MB02ID.f"
	wrkopt = max(i__1,i__2);
#line 325 "MB02ID.f"
	i__1 = *n * *l;
#line 325 "MB02ID.f"
	i__2 = *n * *l;
#line 325 "MB02ID.f"
	dlacpy_("All", &i__1, rb, &dwork[1], &i__2, &b[b_offset], ldb, (
		ftnlen)3);
#line 326 "MB02ID.f"
    }

#line 328 "MB02ID.f"
    pdw = *n * *l * *l + 1;
#line 329 "MB02ID.f"
    i__1 = *m * *k;
#line 329 "MB02ID.f"
    i__2 = *m * *k;
#line 329 "MB02ID.f"
    dlacpy_("All", &i__1, l, &tc[tc_offset], ldtc, &dwork[pdw], &i__2, (
	    ftnlen)3);
#line 330 "MB02ID.f"
    i__1 = *m * *k;
#line 330 "MB02ID.f"
    i__2 = *m * *k;
#line 330 "MB02ID.f"
    i__3 = *ldwork - pdw - (*m * *k + 1) * *l - 1;
#line 330 "MB02ID.f"
    dgeqrf_(&i__1, l, &dwork[pdw], &i__2, &dwork[pdw + *m * *k * *l], &dwork[
	    pdw + (*m * *k + 1) * *l], &i__3, &ierr);
/* Computing MAX */
#line 332 "MB02ID.f"
    i__1 = wrkopt, i__2 = (integer) dwork[pdw + (*m * *k + 1) * *l] + pdw + (*
	    m * *k + 1) * *l - 1;
#line 332 "MB02ID.f"
    wrkopt = max(i__1,i__2);

#line 335 "MB02ID.f"
    i__1 = pdw + *m * *k * *l - 1;
#line 335 "MB02ID.f"
    i__2 = *m * *k + 1;
#line 335 "MB02ID.f"
    for (i__ = pdw; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 336 "MB02ID.f"
	if (dwork[i__] == 0.) {
#line 337 "MB02ID.f"
	    *info = 1;
#line 338 "MB02ID.f"
	    return 0;
#line 339 "MB02ID.f"
	}
#line 340 "MB02ID.f"
/* L10: */
#line 340 "MB02ID.f"
    }

#line 342 "MB02ID.f"
    i__2 = *m * *k;
#line 342 "MB02ID.f"
    i__1 = *n * *l;
#line 342 "MB02ID.f"
    ma02ad_("Upper", l, l, &dwork[pdw], &i__2, &dwork[1], &i__1, (ftnlen)5);
#line 343 "MB02ID.f"
    i__2 = *m * *k;
#line 343 "MB02ID.f"
    i__1 = *m * *k;
#line 343 "MB02ID.f"
    i__3 = *ldwork - pdw - (*m * *k + 1) * *l - 1;
#line 343 "MB02ID.f"
    dorgqr_(&i__2, l, l, &dwork[pdw], &i__1, &dwork[pdw + *m * *k * *l], &
	    dwork[pdw + (*m * *k + 1) * *l], &i__3, &ierr);
/* Computing MAX */
#line 345 "MB02ID.f"
    i__2 = wrkopt, i__1 = (integer) dwork[pdw + (*m * *k + 1) * *l] + pdw + (*
	    m * *k + 1) * *l - 1;
#line 345 "MB02ID.f"
    wrkopt = max(i__2,i__1);
#line 347 "MB02ID.f"
    i__2 = *n - 1;
#line 347 "MB02ID.f"
    i__1 = *m * *k;
#line 347 "MB02ID.f"
    i__3 = *n * *l;
#line 347 "MB02ID.f"
    i__4 = *ldwork - pdw - *m * *k * *l + 1;
#line 347 "MB02ID.f"
    mb02kd_("Row", "Transpose", k, l, m, &i__2, l, &c_b18, &c_b8, &tc[
	    tc_offset], ldtc, &tr[tr_offset], ldtr, &dwork[pdw], &i__1, &
	    dwork[*l + 1], &i__3, &dwork[pdw + *m * *k * *l], &i__4, &ierr, (
	    ftnlen)3, (ftnlen)9);
/* Computing MAX */
#line 350 "MB02ID.f"
    i__2 = wrkopt, i__1 = (integer) dwork[pdw + *m * *k * *l] + pdw + *m * *k 
	    * *l - 1;
#line 350 "MB02ID.f"
    wrkopt = max(i__2,i__1);
#line 351 "MB02ID.f"
    ppr = *n * *l * *l + 1;
#line 352 "MB02ID.f"
    pnr = *n * *l * (*l + *k) + 1;
#line 353 "MB02ID.f"
    i__2 = (*n - 1) * *l;
#line 353 "MB02ID.f"
    i__1 = *n * *l;
#line 353 "MB02ID.f"
    ma02ad_("All", k, &i__2, &tr[tr_offset], ldtr, &dwork[ppr + *l], &i__1, (
	    ftnlen)3);
#line 354 "MB02ID.f"
    i__2 = (*n - 1) * *l;
#line 354 "MB02ID.f"
    i__1 = *n * *l;
#line 354 "MB02ID.f"
    i__3 = *n * *l;
#line 354 "MB02ID.f"
    dlacpy_("All", &i__2, l, &dwork[*l + 1], &i__1, &dwork[pnr + *l], &i__3, (
	    ftnlen)3);
#line 356 "MB02ID.f"
    pt = (*m - 1) * *k + 1;
#line 357 "MB02ID.f"
    pdw = pnr + *n * *l * *l + *l;

/* Computing MIN */
#line 359 "MB02ID.f"
    i__1 = *m, i__3 = *n - 1;
#line 359 "MB02ID.f"
    i__2 = min(i__1,i__3);
#line 359 "MB02ID.f"
    for (i__ = 1; i__ <= i__2; ++i__) {
#line 360 "MB02ID.f"
	i__1 = *n * *l;
#line 360 "MB02ID.f"
	ma02ad_("All", k, l, &tc[pt + tc_dim1], ldtc, &dwork[pdw], &i__1, (
		ftnlen)3);
#line 361 "MB02ID.f"
	pt -= *k;
#line 362 "MB02ID.f"
	pdw += *l;
#line 363 "MB02ID.f"
/* L30: */
#line 363 "MB02ID.f"
    }

#line 365 "MB02ID.f"
    pt = 1;

#line 367 "MB02ID.f"
    i__2 = *n - 1;
#line 367 "MB02ID.f"
    for (i__ = *m + 1; i__ <= i__2; ++i__) {
#line 368 "MB02ID.f"
	i__1 = *n * *l;
#line 368 "MB02ID.f"
	ma02ad_("All", k, l, &tr[pt * tr_dim1 + 1], ldtr, &dwork[pdw], &i__1, 
		(ftnlen)3);
#line 369 "MB02ID.f"
	pt += *l;
#line 370 "MB02ID.f"
	pdw += *l;
#line 371 "MB02ID.f"
/* L40: */
#line 371 "MB02ID.f"
    }

#line 373 "MB02ID.f"
    if (compo) {

/*        Apply the first reduction step to T'*B. */

#line 377 "MB02ID.f"
	i__2 = *n * *l;
#line 377 "MB02ID.f"
	dtrsm_("Left", "Lower", "NonTranspose", "NonUnit", l, rb, &c_b18, &
		dwork[1], &i__2, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
		ftnlen)12, (ftnlen)7);
#line 379 "MB02ID.f"
	i__2 = (*n - 1) * *l;
#line 379 "MB02ID.f"
	i__1 = *n * *l;
#line 379 "MB02ID.f"
	dgemm_("NoTranspose", "NoTranspose", &i__2, rb, l, &c_b18, &dwork[*l 
		+ 1], &i__1, &b[b_offset], ldb, &c_b42, &b[*l + 1 + b_dim1], 
		ldb, (ftnlen)11, (ftnlen)11);
#line 381 "MB02ID.f"
	i__2 = *n * *l;
#line 381 "MB02ID.f"
	dtrsm_("Left", "Lower", "Transpose", "NonUnit", l, rb, &c_b18, &dwork[
		1], &i__2, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (ftnlen)9,
		 (ftnlen)7);
#line 383 "MB02ID.f"
    }

#line 385 "MB02ID.f"
    if (compu) {

/*        Apply the first reduction step to C. */

#line 389 "MB02ID.f"
	i__2 = *n * *l;
#line 389 "MB02ID.f"
	dtrsm_("Left", "Lower", "NonTranspose", "NonUnit", l, rc, &c_b18, &
		dwork[1], &i__2, &c__[c_offset], ldc, (ftnlen)4, (ftnlen)5, (
		ftnlen)12, (ftnlen)7);
#line 391 "MB02ID.f"
	i__2 = (*n - 1) * *l;
#line 391 "MB02ID.f"
	i__1 = *n * *l;
#line 391 "MB02ID.f"
	dgemm_("NoTranspose", "NoTranspose", &i__2, rc, l, &c_b18, &dwork[*l 
		+ 1], &i__1, &c__[c_offset], ldc, &c_b42, &c__[*l + 1 + 
		c_dim1], ldc, (ftnlen)11, (ftnlen)11);
#line 393 "MB02ID.f"
	i__2 = *n * *l;
#line 393 "MB02ID.f"
	dtrsm_("Left", "Lower", "Transpose", "NonUnit", l, rc, &c_b18, &dwork[
		1], &i__2, &c__[c_offset], ldc, (ftnlen)4, (ftnlen)5, (ftnlen)
		9, (ftnlen)7);
#line 395 "MB02ID.f"
    }

#line 397 "MB02ID.f"
    pdi = (*n - 1) * *l + 1;
#line 398 "MB02ID.f"
    i__2 = *n * *l;
#line 398 "MB02ID.f"
    i__1 = *n * *l;
#line 398 "MB02ID.f"
    dlacpy_("Lower", l, l, &dwork[1], &i__2, &dwork[pdi], &i__1, (ftnlen)5);
#line 399 "MB02ID.f"
    i__2 = *n * *l;
#line 399 "MB02ID.f"
    dtrtri_("Lower", "NonUnit", l, &dwork[pdi], &i__2, &ierr, (ftnlen)5, (
	    ftnlen)7);
#line 400 "MB02ID.f"
    i__2 = *l - 1;
#line 400 "MB02ID.f"
    i__1 = *n * *l;
#line 400 "MB02ID.f"
    i__3 = *n * *l;
#line 400 "MB02ID.f"
    ma02ad_("Lower", &i__2, l, &dwork[pdi + 1], &i__1, &dwork[((*n << 1) - 1) 
	    * *l + 1], &i__3, (ftnlen)5);
#line 402 "MB02ID.f"
    i__2 = *l - 1;
#line 402 "MB02ID.f"
    i__1 = *n * *l;
#line 402 "MB02ID.f"
    dlaset_("Lower", &i__2, l, &c_b8, &c_b8, &dwork[pdi + 1], &i__1, (ftnlen)
	    5);
#line 403 "MB02ID.f"
    i__2 = *n * *l;
#line 403 "MB02ID.f"
    i__1 = *n * *l;
#line 403 "MB02ID.f"
    dlacpy_("Upper", l, l, &dwork[pdi], &i__2, &dwork[pnr], &i__1, (ftnlen)5);
#line 404 "MB02ID.f"
    i__2 = *l - 1;
#line 404 "MB02ID.f"
    i__1 = *n * *l;
#line 404 "MB02ID.f"
    dlaset_("Lower", &i__2, l, &c_b8, &c_b8, &dwork[pnr + 1], &i__1, (ftnlen)
	    5);
#line 405 "MB02ID.f"
    i__2 = *n * *l;
#line 405 "MB02ID.f"
    dlaset_("All", l, k, &c_b8, &c_b8, &dwork[ppr], &i__2, (ftnlen)3);
#line 406 "MB02ID.f"
    i__2 = *n * *l;
#line 406 "MB02ID.f"
    dlaset_("All", l, k, &c_b8, &c_b8, &dwork[pnr + *n * *l * *l], &i__2, (
	    ftnlen)3);

#line 408 "MB02ID.f"
    ppi = ppr;
#line 409 "MB02ID.f"
    ppr += *l;
#line 410 "MB02ID.f"
    pni = pnr;
#line 411 "MB02ID.f"
    pnr += *l;
#line 412 "MB02ID.f"
    pdw = (*n << 1) * *l * (*l + *k) + 1;
#line 413 "MB02ID.f"
    len = (*n - 1) * *l;

/*     Determine block size for the involved block Householder */
/*     transformations. */

/* Computing MIN */
#line 418 "MB02ID.f"
    i__1 = *n * *l;
#line 418 "MB02ID.f"
    i__2 = ilaenv_(&c__1, "DGELQF", " ", &i__1, l, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);
#line 418 "MB02ID.f"
    nb = min(i__2,*l);
#line 419 "MB02ID.f"
    kk = pdw + *l * 6 - 1;
/* Computing MAX */
#line 420 "MB02ID.f"
    i__2 = wrkopt, i__1 = kk + *n * *l * nb;
#line 420 "MB02ID.f"
    wrkopt = max(i__2,i__1);
#line 421 "MB02ID.f"
    kk = *ldwork - kk;
#line 422 "MB02ID.f"
    if (kk < *n * *l * nb) {
#line 422 "MB02ID.f"
	nb = kk / (*n * *l);
#line 422 "MB02ID.f"
    }
/* Computing MAX */
#line 423 "MB02ID.f"
    i__3 = *n * *l;
#line 423 "MB02ID.f"
    i__2 = 2, i__1 = ilaenv_(&c__2, "DGELQF", " ", &i__3, l, &c_n1, &c_n1, (
	    ftnlen)6, (ftnlen)1);
#line 423 "MB02ID.f"
    nbmin = max(i__2,i__1);
#line 424 "MB02ID.f"
    if (nb < nbmin) {
#line 424 "MB02ID.f"
	nb = 0;
#line 424 "MB02ID.f"
    }

#line 426 "MB02ID.f"
    i__2 = *n * *l;
#line 426 "MB02ID.f"
    i__1 = *l;
#line 426 "MB02ID.f"
    for (i__ = *l + 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
#line 427 "MB02ID.f"
	i__3 = *l + *k;
#line 427 "MB02ID.f"
	i__4 = *l + *k;
#line 427 "MB02ID.f"
	i__5 = *n * *l;
#line 427 "MB02ID.f"
	i__6 = *n * *l;
#line 427 "MB02ID.f"
	i__7 = *n * *l;
#line 427 "MB02ID.f"
	i__8 = *ldwork - pdw - *l * 6 + 1;
#line 427 "MB02ID.f"
	mb02cu_("Column", l, &i__3, &i__4, &nb, &dwork[1], &i__5, &dwork[ppr],
		 &i__6, &dwork[pnr], &i__7, &rnk, ipvt, &dwork[pdw], &c_b8, &
		dwork[pdw + *l * 6], &i__8, &ierr, (ftnlen)6);
#line 430 "MB02ID.f"
	if (ierr != 0) {

/*           Error return:  The rank condition is (numerically) not */
/*                          satisfied. */

#line 435 "MB02ID.f"
	    *info = 1;
#line 436 "MB02ID.f"
	    return 0;
#line 437 "MB02ID.f"
	}
#line 438 "MB02ID.f"
	i__3 = len - *l;
#line 438 "MB02ID.f"
	i__4 = *l + *k;
#line 438 "MB02ID.f"
	i__5 = *l + *k;
#line 438 "MB02ID.f"
	i__6 = *n * *l;
#line 438 "MB02ID.f"
	i__7 = *n * *l;
#line 438 "MB02ID.f"
	i__8 = *n * *l;
#line 438 "MB02ID.f"
	i__9 = *n * *l;
#line 438 "MB02ID.f"
	i__10 = *n * *l;
#line 438 "MB02ID.f"
	i__11 = *n * *l;
#line 438 "MB02ID.f"
	i__12 = *ldwork - pdw - *l * 6 + 1;
#line 438 "MB02ID.f"
	mb02cv_("Column", "NoStructure", l, &i__3, &i__4, &i__5, &nb, &c_n1, &
		dwork[1], &i__6, &dwork[ppr], &i__7, &dwork[pnr], &i__8, &
		dwork[*l + 1], &i__9, &dwork[ppr + *l], &i__10, &dwork[pnr + *
		l], &i__11, &dwork[pdw], &dwork[pdw + *l * 6], &i__12, &ierr, 
		(ftnlen)6, (ftnlen)11);
#line 443 "MB02ID.f"
	pdi -= *l;
#line 444 "MB02ID.f"
	if (compo) {

/*           Block Gaussian elimination to B. */

#line 448 "MB02ID.f"
	    i__3 = *n * *l;
#line 448 "MB02ID.f"
	    dtrsm_("Left", "Lower", "NonTranspose", "NonUnit", l, rb, &c_b42, 
		    &dwork[1], &i__3, &b[i__ + b_dim1], ldb, (ftnlen)4, (
		    ftnlen)5, (ftnlen)12, (ftnlen)7);
#line 450 "MB02ID.f"
	    if (len > *l) {
#line 451 "MB02ID.f"
		i__3 = len - *l;
#line 451 "MB02ID.f"
		i__4 = *n * *l;
#line 451 "MB02ID.f"
		dgemm_("NonTranspose", "NonTranspose", &i__3, rb, l, &c_b18, &
			dwork[*l + 1], &i__4, &b[i__ + b_dim1], ldb, &c_b18, &
			b[i__ + *l + b_dim1], ldb, (ftnlen)12, (ftnlen)12);
#line 454 "MB02ID.f"
	    }
#line 455 "MB02ID.f"
	}
#line 456 "MB02ID.f"
	if (compu) {

/*           Block Gaussian elimination to C. */

#line 460 "MB02ID.f"
	    i__3 = *n * *l;
#line 460 "MB02ID.f"
	    dtrsm_("Left", "Lower", "NonTranspose", "NonUnit", l, rc, &c_b42, 
		    &dwork[1], &i__3, &c__[i__ + c_dim1], ldc, (ftnlen)4, (
		    ftnlen)5, (ftnlen)12, (ftnlen)7);
#line 462 "MB02ID.f"
	    if (len > *l) {
#line 463 "MB02ID.f"
		i__3 = len - *l;
#line 463 "MB02ID.f"
		i__4 = *n * *l;
#line 463 "MB02ID.f"
		dgemm_("NonTranspose", "NonTranspose", &i__3, rc, l, &c_b18, &
			dwork[*l + 1], &i__4, &c__[i__ + c_dim1], ldc, &c_b18,
			 &c__[i__ + *l + c_dim1], ldc, (ftnlen)12, (ftnlen)12)
			;
#line 466 "MB02ID.f"
	    }
#line 467 "MB02ID.f"
	}
#line 468 "MB02ID.f"
	i__3 = *n * *l;
#line 468 "MB02ID.f"
	dlaset_("All", l, l, &c_b8, &c_b8, &dwork[pdi], &i__3, (ftnlen)3);
#line 469 "MB02ID.f"
	i__3 = i__ + *l - 1;
#line 469 "MB02ID.f"
	i__4 = *l + *k;
#line 469 "MB02ID.f"
	i__5 = *l + *k;
#line 469 "MB02ID.f"
	i__6 = *n * *l;
#line 469 "MB02ID.f"
	i__7 = *n * *l;
#line 469 "MB02ID.f"
	i__8 = *n * *l;
#line 469 "MB02ID.f"
	i__9 = *n * *l;
#line 469 "MB02ID.f"
	i__10 = *n * *l;
#line 469 "MB02ID.f"
	i__11 = *n * *l;
#line 469 "MB02ID.f"
	i__12 = *ldwork - pdw - *l * 6 + 1;
#line 469 "MB02ID.f"
	mb02cv_("Column", "Triangular", l, &i__3, &i__4, &i__5, &nb, &c_n1, &
		dwork[1], &i__6, &dwork[ppr], &i__7, &dwork[pnr], &i__8, &
		dwork[pdi], &i__9, &dwork[ppi], &i__10, &dwork[pni], &i__11, &
		dwork[pdw], &dwork[pdw + *l * 6], &i__12, &ierr, (ftnlen)6, (
		ftnlen)10);
#line 474 "MB02ID.f"
	if (compo) {

/*           Apply block Gaussian elimination to B. */

#line 478 "MB02ID.f"
	    i__3 = i__ - 1;
#line 478 "MB02ID.f"
	    i__4 = *n * *l;
#line 478 "MB02ID.f"
	    dgemm_("NoTranspose", "NoTranspose", &i__3, rb, l, &c_b18, &dwork[
		    pdi], &i__4, &b[i__ + b_dim1], ldb, &c_b18, &b[b_offset], 
		    ldb, (ftnlen)11, (ftnlen)11);
#line 480 "MB02ID.f"
	    i__3 = *n * *l;
#line 480 "MB02ID.f"
	    dtrmm_("Left", "Upper", "NonTranspose", "NonUnit", l, rb, &c_b18, 
		    &dwork[(*n - 1) * *l + 1], &i__3, &b[i__ + b_dim1], ldb, (
		    ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)7);
#line 482 "MB02ID.f"
	}
#line 483 "MB02ID.f"
	if (compu) {

/*           Apply block Gaussian elimination to C. */

#line 487 "MB02ID.f"
	    i__3 = i__ - 1;
#line 487 "MB02ID.f"
	    i__4 = *n * *l;
#line 487 "MB02ID.f"
	    dgemm_("NonTranspose", "NonTranspose", &i__3, rc, l, &c_b18, &
		    dwork[pdi], &i__4, &c__[i__ + c_dim1], ldc, &c_b18, &c__[
		    c_offset], ldc, (ftnlen)12, (ftnlen)12);
#line 489 "MB02ID.f"
	    i__3 = *n * *l;
#line 489 "MB02ID.f"
	    dtrmm_("Left", "Upper", "NonTranspose", "NonUnit", l, rc, &c_b18, 
		    &dwork[(*n - 1) * *l + 1], &i__3, &c__[i__ + c_dim1], ldc,
		     (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)7);
#line 491 "MB02ID.f"
	}
#line 492 "MB02ID.f"
	len -= *l;
#line 493 "MB02ID.f"
	pnr += *l;
#line 494 "MB02ID.f"
	ppr += *l;
#line 495 "MB02ID.f"
/* L50: */
#line 495 "MB02ID.f"
    }

#line 497 "MB02ID.f"
    if (compu) {
#line 498 "MB02ID.f"
	i__1 = *m * *k;
#line 498 "MB02ID.f"
	i__2 = *ldwork - *m * *k * *rc;
#line 498 "MB02ID.f"
	mb02kd_("Column", "NonTranspose", k, l, m, n, rc, &c_b18, &c_b8, &tc[
		tc_offset], ldtc, &tr[tr_offset], ldtr, &c__[c_offset], ldc, &
		dwork[1], &i__1, &dwork[*m * *k * *rc + 1], &i__2, &ierr, (
		ftnlen)6, (ftnlen)12);
/* Computing MAX */
#line 501 "MB02ID.f"
	i__1 = wrkopt, i__2 = (integer) dwork[*m * *k * *rc + 1] + *m * *k * *
		rc;
#line 501 "MB02ID.f"
	wrkopt = max(i__1,i__2);
#line 502 "MB02ID.f"
	i__1 = *m * *k;
#line 502 "MB02ID.f"
	i__2 = *m * *k;
#line 502 "MB02ID.f"
	dlacpy_("All", &i__1, rc, &dwork[1], &i__2, &c__[c_offset], ldc, (
		ftnlen)3);
#line 503 "MB02ID.f"
    }
#line 504 "MB02ID.f"
    dwork[1] = (doublereal) wrkopt;
#line 505 "MB02ID.f"
    return 0;

/* *** Last line of MB02ID *** */
} /* mb02id_ */

