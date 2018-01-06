#line 1 "MB03ZD.f"
/* MB03ZD.f -- translated by f2c (version 20100827).
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

#line 1 "MB03ZD.f"
/* Table of constant values */

static doublereal c_b26 = 1.;
static doublereal c_b34 = -1.;
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b272 = 0.;

/* Subroutine */ int mb03zd_(char *which, char *meth, char *stab, char *
	balanc, char *ortbal, logical *select, integer *n, integer *mm, 
	integer *ilo, doublereal *scale, doublereal *s, integer *lds, 
	doublereal *t, integer *ldt, doublereal *g, integer *ldg, doublereal *
	u1, integer *ldu1, doublereal *u2, integer *ldu2, doublereal *v1, 
	integer *ldv1, doublereal *v2, integer *ldv2, integer *m, doublereal *
	wr, doublereal *wi, doublereal *us, integer *ldus, doublereal *uu, 
	integer *lduu, logical *lwork, integer *iwork, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen which_len, ftnlen meth_len, 
	ftnlen stab_len, ftnlen balanc_len, ftnlen ortbal_len)
{
    /* System generated locals */
    integer g_dim1, g_offset, s_dim1, s_offset, t_dim1, t_offset, u1_dim1, 
	    u1_offset, u2_dim1, u2_offset, us_dim1, us_offset, uu_dim1, 
	    uu_offset, v1_dim1, v1_offset, v2_dim1, v2_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, j, k, pw, pdw;
    static logical lus, luu, lbef, lbal, lall, pair;
    static integer ierr;
    static doublereal temp;
    static logical lext;
    extern /* Subroutine */ int mb04di_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen), dscal_(integer *, 
	    doublereal *, doublereal *, integer *), dgemm_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen), mb03td_(char *, char *, logical *, 
	    logical *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), mb03za_(char *, char *, char *, char *
	    , char *, logical *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb01ux_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen), daxpy_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *), dgeqp3_(integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *), dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dgeqrf_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), dlaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen), xerbla_(char *, integer *, ftnlen), dorgqr_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);
    static integer wrkmin, wrkopt;


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

/*     To compute the stable and unstable invariant subspaces for a */
/*     Hamiltonian matrix with no eigenvalues on the imaginary axis, */
/*     using the output of the SLICOT Library routine MB03XD. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     WHICH   CHARACTER*1 */
/*             Specifies the cluster of eigenvalues for which the */
/*             invariant subspaces are computed: */
/*             = 'A':  select all n eigenvalues; */
/*             = 'S':  select a cluster of eigenvalues specified by */
/*                     SELECT. */

/*     METH    CHARACTER*1 */
/*             If WHICH = 'A' this parameter specifies the method to be */
/*             used for computing bases of the invariant subspaces: */
/*             = 'S':  compute the n-dimensional basis from a set of */
/*                     n vectors; */
/*             = 'L':  compute the n-dimensional basis from a set of */
/*                     2*n vectors. */
/*             When in doubt, use METH = 'S'. In some cases, METH = 'L' */
/*             may result in more accurately computed invariant */
/*             subspaces, see [1]. */

/*     STAB    CHARACTER*1 */
/*             Specifies the type of invariant subspaces to be computed: */
/*             = 'S':  compute the stable invariant subspace, i.e., the */
/*                     invariant subspace belonging to those selected */
/*                     eigenvalues that have negative real part; */
/*             = 'U':  compute the unstable invariant subspace, i.e., */
/*                     the invariant subspace belonging to those */
/*                     selected eigenvalues that have positive real */
/*                     part; */
/*             = 'B':  compute both the stable and unstable invariant */
/*                     subspaces. */

/*     BALANC  CHARACTER*1 */
/*             Specifies the type of inverse balancing transformation */
/*             required: */
/*             = 'N':  do nothing; */
/*             = 'P':  do inverse transformation for permutation only; */
/*             = 'S':  do inverse transformation for scaling only; */
/*             = 'B':  do inverse transformations for both permutation */
/*                     and scaling. */
/*             BALANC must be the same as the argument BALANC supplied to */
/*             MB03XD. Note that if the data is further post-processed, */
/*             e.g., for solving an algebraic Riccati equation, it is */
/*             recommended to delay inverse balancing (in particular the */
/*             scaling part) and apply it to the final result only, */
/*             see [2]. */

/*     ORTBAL  CHARACTER*1 */
/*             If BALANC <> 'N', this option specifies how inverse */
/*             balancing is applied to the computed invariant subspaces: */
/*             = 'B':  apply inverse balancing before orthogonal bases */
/*                     for the invariant subspaces are computed; */
/*             = 'A':  apply inverse balancing after orthogonal bases */
/*                     for the invariant subspaces have been computed; */
/*                     this may yield non-orthogonal bases if */
/*                     BALANC = 'S' or BALANC = 'B'. */

/*     SELECT  (input) LOGICAL array, dimension (N) */
/*             If WHICH = 'S', SELECT specifies the eigenvalues */
/*             corresponding to the positive and negative square */
/*             roots of the eigenvalues of S*T in the selected cluster. */
/*             To select a real eigenvalue w(j), SELECT(j) must be set */
/*             to .TRUE.. To select a complex conjugate pair of */
/*             eigenvalues w(j) and w(j+1), corresponding to a 2-by-2 */
/*             diagonal block, both SELECT(j) and SELECT(j+1) must be set */
/*             to .TRUE.; a complex conjugate pair of eigenvalues must be */
/*             either both included in the cluster or both excluded. */
/*             This array is not referenced if WHICH = 'A'. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices S, T and G. N >= 0. */

/*     MM      (input) INTEGER */
/*             The number of columns in the arrays US and/or UU. */
/*             If WHICH = 'A' and METH = 'S',  MM >=   N; */
/*             if WHICH = 'A' and METH = 'L',  MM >= 2*N; */
/*             if WHICH = 'S',                 MM >=   M. */
/*             The minimal values above for MM give the numbers of */
/*             vectors to be used for computing a basis for the */
/*             invariant subspace(s). */

/*     ILO     (input) INTEGER */
/*             If BALANC <> 'N', then ILO is the integer returned by */
/*             MB03XD.  1 <= ILO <= N+1. */

/*     SCALE   (input) DOUBLE PRECISION array, dimension (N) */
/*             If BALANC <> 'N', the leading N elements of this array */
/*             must contain details of the permutation and scaling */
/*             factors, as returned by MB03XD. */
/*             This array is not referenced if BALANC = 'N'. */

/*     S       (input/output) DOUBLE PRECISION array, dimension (LDS,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix S in real Schur form. */
/*             On exit, the leading N-by-N part of this array is */
/*             overwritten. */

/*     LDS     INTEGER */
/*             The leading dimension of the array S.  LDS >= max(1,N). */

/*     T       (input/output) DOUBLE PRECISION array, dimension (LDT,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper triangular matrix T. */
/*             On exit, the leading N-by-N part of this array is */
/*             overwritten. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T.  LDT >= max(1,N). */

/*     G       (input/output) DOUBLE PRECISION array, dimension (LDG,N) */
/*             On entry, if METH = 'L', the leading N-by-N part of this */
/*             array must contain a general matrix G. */
/*             On exit, if METH = 'L', the leading N-by-N part of this */
/*             array is overwritten. */
/*             This array is not referenced if METH = 'S'. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G.  LDG >= 1. */
/*             LDG >= max(1,N) if METH = 'L'. */

/*     U1      (input/output) DOUBLE PRECISION array, dimension (LDU1,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the (1,1) block of an orthogonal symplectic */
/*             matrix U. */
/*             On exit, this array is overwritten. */

/*     LDU1    INTEGER */
/*             The leading dimension of the array U1.  LDU1 >= MAX(1,N). */

/*     U2      (input/output) DOUBLE PRECISION array, dimension (LDU2,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the (2,1) block of an orthogonal symplectic */
/*             matrix U. */
/*             On exit, this array is overwritten. */

/*     LDU2    INTEGER */
/*             The leading dimension of the array U2.  LDU2 >= MAX(1,N). */

/*     V1      (input/output) DOUBLE PRECISION array, dimension (LDV1,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the (1,1) block of an orthogonal symplectic */
/*             matrix V. */
/*             On exit, this array is overwritten. */

/*     LDV1    INTEGER */
/*             The leading dimension of the array V1.  LDV1 >= MAX(1,N). */

/*     V2      (input/output) DOUBLE PRECISION array, dimension (LDV1,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the (2,1) block of an orthogonal symplectic */
/*             matrix V. */
/*             On exit, this array is overwritten. */

/*     LDV2    INTEGER */
/*             The leading dimension of the array V2.  LDV2 >= MAX(1,N). */

/*     M       (output) INTEGER */
/*             The number of selected eigenvalues. */

/*     WR      (output) DOUBLE PRECISION array, dimension (M) */
/*     WI      (output) DOUBLE PRECISION array, dimension (M) */
/*             On exit, the leading M elements of WR and WI contain the */
/*             real and imaginary parts, respectively, of the selected */
/*             eigenvalues that have nonpositive real part. Complex */
/*             conjugate pairs of eigenvalues with real part not equal */
/*             to zero will appear consecutively with the eigenvalue */
/*             having the positive imaginary part first. Note that, due */
/*             to roundoff errors, these numbers may differ from the */
/*             eigenvalues computed by MB03XD. */

/*     US      (output) DOUBLE PRECISION array, dimension (LDUS,MM) */
/*             On exit, if STAB = 'S' or STAB = 'B', the leading 2*N-by-M */
/*             part of this array contains a basis for the stable */
/*             invariant subspace belonging to the selected eigenvalues. */
/*             This basis is orthogonal unless ORTBAL = 'A'. */

/*     LDUS    INTEGER */
/*             The leading dimension of the array US.  LDUS >= 1. */
/*             If STAB = 'S' or STAB = 'B',  LDUS >= 2*N. */

/*     UU      (output) DOUBLE PRECISION array, dimension (LDUU,MM) */
/*             On exit, if STAB = 'U' or STAB = 'B', the leading 2*N-by-M */
/*             part of this array contains a basis for the unstable */
/*             invariant subspace belonging to the selected eigenvalues. */
/*             This basis is orthogonal unless ORTBAL = 'A'. */

/*     LDUU    INTEGER */
/*             The leading dimension of the array UU.  LDUU >= 1. */
/*             If STAB = 'U' or STAB = 'B',  LDUU >= 2*N. */

/*     Workspace */

/*     LWORK   LOGICAL array, dimension (2*N) */
/*             This array is only referenced if WHICH = 'A' and */
/*             METH = 'L'. */

/*     IWORK   INTEGER array, dimension (2*N), */
/*             This array is only referenced if WHICH = 'A' and */
/*             METH = 'L'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -35,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             If WHICH = 'S' or METH = 'S': */
/*                LDWORK >= MAX( 1, 4*M*M + MAX( 8*M, 4*N ) ). */
/*             If WHICH = 'A' and METH = 'L' and */
/*                ( STAB = 'U' or STAB = 'S' ): */
/*                LDWORK >= MAX( 1, 2*N*N + 2*N, 8*N ). */
/*             If WHICH = 'A' and METH = 'L' and STAB = 'B': */
/*                LDWORK >= 8*N + 1. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  some of the selected eigenvalues are on or too close */
/*                   to the imaginary axis; */
/*             = 2:  reordering of the product S*T in routine MB03ZA */
/*                   failed because some eigenvalues are too close to */
/*                   separate; */
/*             = 3:  the QR algorithm failed to compute some Schur form */
/*                   in MB03ZA; */
/*             = 4:  reordering of the Hamiltonian Schur form in routine */
/*                   MB03TD failed because some eigenvalues are too close */
/*                   to separate. */

/*     METHOD */

/*     This is an implementation of Algorithm 1 in [1]. */

/*     NUMERICAL ASPECTS */

/*     The method is strongly backward stable for an embedded */
/*     (skew-)Hamiltonian matrix, see [1]. Although good results have */
/*     been reported if the eigenvalues are not too close to the */
/*     imaginary axis, the method is not backward stable for the original */
/*     Hamiltonian matrix itself. */

/*     REFERENCES */

/*     [1] Benner, P., Mehrmann, V., and Xu, H. */
/*         A new method for computing the stable invariant subspace of a */
/*         real Hamiltonian matrix, J. Comput. Appl. Math., 86, */
/*         pp. 17-43, 1997. */

/*     [2] Benner, P. */
/*         Symplectic balancing of Hamiltonian matrices. */
/*         SIAM J. Sci. Comput., 22 (5), pp. 1885-1904, 2000. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DHASUB). */

/*     KEYWORDS */

/*     Hamiltonian matrix, invariant subspace. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Decode and check input parameters. */

#line 338 "MB03ZD.f"
    /* Parameter adjustments */
#line 338 "MB03ZD.f"
    --select;
#line 338 "MB03ZD.f"
    --scale;
#line 338 "MB03ZD.f"
    s_dim1 = *lds;
#line 338 "MB03ZD.f"
    s_offset = 1 + s_dim1;
#line 338 "MB03ZD.f"
    s -= s_offset;
#line 338 "MB03ZD.f"
    t_dim1 = *ldt;
#line 338 "MB03ZD.f"
    t_offset = 1 + t_dim1;
#line 338 "MB03ZD.f"
    t -= t_offset;
#line 338 "MB03ZD.f"
    g_dim1 = *ldg;
#line 338 "MB03ZD.f"
    g_offset = 1 + g_dim1;
#line 338 "MB03ZD.f"
    g -= g_offset;
#line 338 "MB03ZD.f"
    u1_dim1 = *ldu1;
#line 338 "MB03ZD.f"
    u1_offset = 1 + u1_dim1;
#line 338 "MB03ZD.f"
    u1 -= u1_offset;
#line 338 "MB03ZD.f"
    u2_dim1 = *ldu2;
#line 338 "MB03ZD.f"
    u2_offset = 1 + u2_dim1;
#line 338 "MB03ZD.f"
    u2 -= u2_offset;
#line 338 "MB03ZD.f"
    v1_dim1 = *ldv1;
#line 338 "MB03ZD.f"
    v1_offset = 1 + v1_dim1;
#line 338 "MB03ZD.f"
    v1 -= v1_offset;
#line 338 "MB03ZD.f"
    v2_dim1 = *ldv2;
#line 338 "MB03ZD.f"
    v2_offset = 1 + v2_dim1;
#line 338 "MB03ZD.f"
    v2 -= v2_offset;
#line 338 "MB03ZD.f"
    --wr;
#line 338 "MB03ZD.f"
    --wi;
#line 338 "MB03ZD.f"
    us_dim1 = *ldus;
#line 338 "MB03ZD.f"
    us_offset = 1 + us_dim1;
#line 338 "MB03ZD.f"
    us -= us_offset;
#line 338 "MB03ZD.f"
    uu_dim1 = *lduu;
#line 338 "MB03ZD.f"
    uu_offset = 1 + uu_dim1;
#line 338 "MB03ZD.f"
    uu -= uu_offset;
#line 338 "MB03ZD.f"
    --lwork;
#line 338 "MB03ZD.f"
    --iwork;
#line 338 "MB03ZD.f"
    --dwork;
#line 338 "MB03ZD.f"

#line 338 "MB03ZD.f"
    /* Function Body */
#line 338 "MB03ZD.f"
    lall = lsame_(which, "A", (ftnlen)1, (ftnlen)1);
#line 339 "MB03ZD.f"
    if (lall) {
#line 340 "MB03ZD.f"
	lext = lsame_(meth, "L", (ftnlen)1, (ftnlen)1);
#line 341 "MB03ZD.f"
    } else {
#line 342 "MB03ZD.f"
	lext = FALSE_;
#line 343 "MB03ZD.f"
    }
#line 344 "MB03ZD.f"
    lus = lsame_(stab, "S", (ftnlen)1, (ftnlen)1) || lsame_(stab, "B", (
	    ftnlen)1, (ftnlen)1);
#line 345 "MB03ZD.f"
    luu = lsame_(stab, "U", (ftnlen)1, (ftnlen)1) || lsame_(stab, "B", (
	    ftnlen)1, (ftnlen)1);
#line 346 "MB03ZD.f"
    lbal = lsame_(balanc, "P", (ftnlen)1, (ftnlen)1) || lsame_(balanc, "S", (
	    ftnlen)1, (ftnlen)1) || lsame_(balanc, "B", (ftnlen)1, (ftnlen)1);
#line 348 "MB03ZD.f"
    lbef = FALSE_;
#line 349 "MB03ZD.f"
    if (lbal) {
#line 349 "MB03ZD.f"
	lbef = lsame_(ortbal, "B", (ftnlen)1, (ftnlen)1);
#line 349 "MB03ZD.f"
    }

#line 352 "MB03ZD.f"
    wrkmin = 1;
#line 353 "MB03ZD.f"
    wrkopt = wrkmin;

#line 355 "MB03ZD.f"
    *info = 0;

#line 357 "MB03ZD.f"
    if (! lall && ! lsame_(which, "S", (ftnlen)1, (ftnlen)1)) {
#line 358 "MB03ZD.f"
	*info = -1;
#line 359 "MB03ZD.f"
    } else if (lall && (! lext && ! lsame_(meth, "S", (ftnlen)1, (ftnlen)1))) 
	    {
#line 361 "MB03ZD.f"
	*info = -2;
#line 362 "MB03ZD.f"
    } else if (! lus && ! luu) {
#line 363 "MB03ZD.f"
	*info = -3;
#line 364 "MB03ZD.f"
    } else if (! lbal && ! lsame_(balanc, "N", (ftnlen)1, (ftnlen)1)) {
#line 365 "MB03ZD.f"
	*info = -4;
#line 366 "MB03ZD.f"
    } else if (lbal && (! lbef && ! lsame_(ortbal, "A", (ftnlen)1, (ftnlen)1))
	    ) {
#line 368 "MB03ZD.f"
	*info = -5;
#line 369 "MB03ZD.f"
    } else {
#line 370 "MB03ZD.f"
	if (lall) {
#line 371 "MB03ZD.f"
	    *m = *n;
#line 372 "MB03ZD.f"
	} else {

/*           Set M to the dimension of the specified invariant subspace. */

#line 376 "MB03ZD.f"
	    *m = 0;
#line 377 "MB03ZD.f"
	    pair = FALSE_;
#line 378 "MB03ZD.f"
	    i__1 = *n;
#line 378 "MB03ZD.f"
	    for (k = 1; k <= i__1; ++k) {
#line 379 "MB03ZD.f"
		if (pair) {
#line 380 "MB03ZD.f"
		    pair = FALSE_;
#line 381 "MB03ZD.f"
		} else {
#line 382 "MB03ZD.f"
		    if (k < *n) {
#line 383 "MB03ZD.f"
			if (s[k + 1 + k * s_dim1] == 0.) {
#line 384 "MB03ZD.f"
			    if (select[k]) {
#line 384 "MB03ZD.f"
				++(*m);
#line 384 "MB03ZD.f"
			    }
#line 386 "MB03ZD.f"
			} else {
#line 387 "MB03ZD.f"
			    pair = TRUE_;
#line 388 "MB03ZD.f"
			    if (select[k] || select[k + 1]) {
#line 388 "MB03ZD.f"
				*m += 2;
#line 388 "MB03ZD.f"
			    }
#line 390 "MB03ZD.f"
			}
#line 391 "MB03ZD.f"
		    } else {
#line 392 "MB03ZD.f"
			if (select[*n]) {
#line 392 "MB03ZD.f"
			    ++(*m);
#line 392 "MB03ZD.f"
			}
#line 394 "MB03ZD.f"
		    }
#line 395 "MB03ZD.f"
		}
#line 396 "MB03ZD.f"
/* L10: */
#line 396 "MB03ZD.f"
	    }
#line 397 "MB03ZD.f"
	}

/*        Compute workspace requirements. */

#line 401 "MB03ZD.f"
	if (! lext) {
/* Computing MAX */
/* Computing MAX */
#line 402 "MB03ZD.f"
	    i__3 = *m << 3, i__4 = *n << 2;
#line 402 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (*m << 2) * *m + max(i__3,i__4);
#line 402 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);
#line 403 "MB03ZD.f"
	} else {
#line 404 "MB03ZD.f"
	    if (lus && luu) {
/* Computing MAX */
#line 405 "MB03ZD.f"
		i__1 = wrkopt, i__2 = (*n << 3) + 1;
#line 405 "MB03ZD.f"
		wrkopt = max(i__1,i__2);
#line 406 "MB03ZD.f"
	    } else {
/* Computing MAX */
#line 407 "MB03ZD.f"
		i__1 = wrkopt, i__2 = (*n << 1) * *n + (*n << 1), i__1 = max(
			i__1,i__2), i__2 = *n << 3;
#line 407 "MB03ZD.f"
		wrkopt = max(i__1,i__2);
#line 408 "MB03ZD.f"
	    }
#line 409 "MB03ZD.f"
	}

#line 411 "MB03ZD.f"
	if (*n < 0) {
#line 412 "MB03ZD.f"
	    *info = -7;
#line 413 "MB03ZD.f"
	} else if (*mm < *m || lext && *mm < *n << 1) {
#line 414 "MB03ZD.f"
	    *info = -8;
#line 415 "MB03ZD.f"
	} else if (lbal && (*ilo < 1 || *ilo > *n + 1)) {
#line 416 "MB03ZD.f"
	    *info = -9;
#line 417 "MB03ZD.f"
	} else if (*lds < max(1,*n)) {
#line 418 "MB03ZD.f"
	    *info = -12;
#line 419 "MB03ZD.f"
	} else if (*ldt < max(1,*n)) {
#line 420 "MB03ZD.f"
	    *info = -14;
#line 421 "MB03ZD.f"
	} else if (*ldg < 1 || lext && *ldg < *n) {
#line 422 "MB03ZD.f"
	    *info = -16;
#line 423 "MB03ZD.f"
	} else if (*ldu1 < max(1,*n)) {
#line 424 "MB03ZD.f"
	    *info = -18;
#line 425 "MB03ZD.f"
	} else if (*ldu2 < max(1,*n)) {
#line 426 "MB03ZD.f"
	    *info = -20;
#line 427 "MB03ZD.f"
	} else if (*ldv1 < max(1,*n)) {
#line 428 "MB03ZD.f"
	    *info = -22;
#line 429 "MB03ZD.f"
	} else if (*ldv2 < max(1,*n)) {
#line 430 "MB03ZD.f"
	    *info = -24;
#line 431 "MB03ZD.f"
	} else if (*ldus < 1 || lus && *ldus < *n << 1) {
#line 432 "MB03ZD.f"
	    *info = -29;
#line 433 "MB03ZD.f"
	} else if (*lduu < 1 || luu && *lduu < *n << 1) {
#line 434 "MB03ZD.f"
	    *info = -31;
#line 435 "MB03ZD.f"
	} else if (*ldwork < wrkmin) {
#line 436 "MB03ZD.f"
	    *info = -35;
#line 437 "MB03ZD.f"
	    dwork[1] = (doublereal) wrkmin;
#line 438 "MB03ZD.f"
	}
#line 439 "MB03ZD.f"
    }

/*     Return if there were illegal values. */

#line 443 "MB03ZD.f"
    if (*info != 0) {
#line 444 "MB03ZD.f"
	i__1 = -(*info);
#line 444 "MB03ZD.f"
	xerbla_("MB03ZD", &i__1, (ftnlen)6);
#line 445 "MB03ZD.f"
	return 0;
#line 446 "MB03ZD.f"
    }

/*     Quick return if possible. */

#line 450 "MB03ZD.f"
    if (min(*m,*n) == 0) {
#line 451 "MB03ZD.f"
	dwork[1] = 1.;
#line 452 "MB03ZD.f"
	return 0;
#line 453 "MB03ZD.f"
    }
#line 454 "MB03ZD.f"
    wrkopt = wrkmin;

#line 456 "MB03ZD.f"
    if (! lext) {

/*        Workspace requirements: 4*M*M + MAX( 8*M, 4*N ). */

#line 460 "MB03ZD.f"
	pw = 1;
#line 461 "MB03ZD.f"
	pdw = pw + (*m << 2) * *m;
#line 462 "MB03ZD.f"
	i__1 = *m << 1;
#line 462 "MB03ZD.f"
	i__2 = *ldwork - pdw + 1;
#line 462 "MB03ZD.f"
	mb03za_("No Update", "Update", "Update", "Init", which, &select[1], n,
		 &s[s_offset], lds, &t[t_offset], ldt, &g[g_offset], ldg, &u1[
		u1_offset], ldu1, &u2[u2_offset], ldu2, &v1[v1_offset], ldv1, 
		&v2[v2_offset], ldv2, &dwork[pw], &i__1, &wr[1], &wi[1], m, &
		dwork[pdw], &i__2, &ierr, (ftnlen)9, (ftnlen)6, (ftnlen)6, (
		ftnlen)4, (ftnlen)1);
#line 466 "MB03ZD.f"
	if (ierr != 0) {
#line 466 "MB03ZD.f"
	    goto L250;
#line 466 "MB03ZD.f"
	}

#line 469 "MB03ZD.f"
	pdw = pw + (*m << 1) * *m;
#line 470 "MB03ZD.f"
	i__1 = *m << 1;
#line 470 "MB03ZD.f"
	i__2 = *ldwork - pdw + 1;
#line 470 "MB03ZD.f"
	mb01ux_("Right", "Upper", "No Transpose", n, m, &c_b26, &dwork[pw], &
		i__1, &v1[v1_offset], ldv1, &dwork[pdw], &i__2, &ierr, (
		ftnlen)5, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
#line 473 "MB03ZD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
#line 473 "MB03ZD.f"
	wrkopt = max(i__1,i__2);

#line 475 "MB03ZD.f"
	if (lus) {
#line 475 "MB03ZD.f"
	    dlacpy_("All", n, m, &v1[v1_offset], ldv1, &us[us_offset], ldus, (
		    ftnlen)3);
#line 475 "MB03ZD.f"
	}
#line 477 "MB03ZD.f"
	if (luu) {
#line 477 "MB03ZD.f"
	    dlacpy_("All", n, m, &v1[v1_offset], ldv1, &uu[uu_offset], lduu, (
		    ftnlen)3);
#line 477 "MB03ZD.f"
	}

#line 480 "MB03ZD.f"
	i__1 = *m << 1;
#line 480 "MB03ZD.f"
	i__2 = *ldwork - pdw + 1;
#line 480 "MB03ZD.f"
	mb01ux_("Right", "Upper", "No Transpose", n, m, &c_b26, &dwork[pw + *
		m], &i__1, &u1[u1_offset], ldu1, &dwork[pdw], &i__2, &ierr, (
		ftnlen)5, (ftnlen)5, (ftnlen)12);

#line 484 "MB03ZD.f"
	if (lus) {
#line 485 "MB03ZD.f"
	    i__1 = *m;
#line 485 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 486 "MB03ZD.f"
		daxpy_(n, &c_b34, &u1[j * u1_dim1 + 1], &c__1, &us[j * 
			us_dim1 + 1], &c__1);
#line 487 "MB03ZD.f"
/* L20: */
#line 487 "MB03ZD.f"
	    }
#line 488 "MB03ZD.f"
	}
#line 489 "MB03ZD.f"
	if (luu) {
#line 490 "MB03ZD.f"
	    i__1 = *m;
#line 490 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 491 "MB03ZD.f"
		daxpy_(n, &c_b26, &u1[j * u1_dim1 + 1], &c__1, &uu[j * 
			uu_dim1 + 1], &c__1);
#line 492 "MB03ZD.f"
/* L30: */
#line 492 "MB03ZD.f"
	    }
#line 493 "MB03ZD.f"
	}

#line 495 "MB03ZD.f"
	i__1 = *m << 1;
#line 495 "MB03ZD.f"
	i__2 = *ldwork - pdw + 1;
#line 495 "MB03ZD.f"
	mb01ux_("Right", "Upper", "No Transpose", n, m, &c_b34, &dwork[pw], &
		i__1, &v2[v2_offset], ldv2, &dwork[pdw], &i__2, &ierr, (
		ftnlen)5, (ftnlen)5, (ftnlen)12);

#line 499 "MB03ZD.f"
	if (lus) {
#line 499 "MB03ZD.f"
	    dlacpy_("All", n, m, &v2[v2_offset], ldv2, &us[*n + 1 + us_dim1], 
		    ldus, (ftnlen)3);
#line 499 "MB03ZD.f"
	}
#line 501 "MB03ZD.f"
	if (luu) {
#line 501 "MB03ZD.f"
	    dlacpy_("All", n, m, &v2[v2_offset], ldv2, &uu[*n + 1 + uu_dim1], 
		    lduu, (ftnlen)3);
#line 501 "MB03ZD.f"
	}

#line 504 "MB03ZD.f"
	i__1 = *m << 1;
#line 504 "MB03ZD.f"
	i__2 = *ldwork - pdw + 1;
#line 504 "MB03ZD.f"
	mb01ux_("Right", "Upper", "No Transpose", n, m, &c_b26, &dwork[pw + *
		m], &i__1, &u2[u2_offset], ldu2, &dwork[pdw], &i__2, &ierr, (
		ftnlen)5, (ftnlen)5, (ftnlen)12);

#line 508 "MB03ZD.f"
	if (lus) {
#line 509 "MB03ZD.f"
	    i__1 = *m;
#line 509 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 510 "MB03ZD.f"
		daxpy_(n, &c_b26, &u2[j * u2_dim1 + 1], &c__1, &us[*n + 1 + j 
			* us_dim1], &c__1);
#line 511 "MB03ZD.f"
/* L40: */
#line 511 "MB03ZD.f"
	    }
#line 512 "MB03ZD.f"
	}
#line 513 "MB03ZD.f"
	if (luu) {
#line 514 "MB03ZD.f"
	    i__1 = *m;
#line 514 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 515 "MB03ZD.f"
		daxpy_(n, &c_b34, &u2[j * u2_dim1 + 1], &c__1, &uu[*n + 1 + j 
			* uu_dim1], &c__1);
#line 516 "MB03ZD.f"
/* L50: */
#line 516 "MB03ZD.f"
	    }
#line 517 "MB03ZD.f"
	}

/*        Orthonormalize obtained bases and apply inverse balancing */
/*        transformation. */

#line 522 "MB03ZD.f"
	if (lbal && lbef) {
#line 523 "MB03ZD.f"
	    if (lus) {
#line 523 "MB03ZD.f"
		mb04di_(balanc, "Positive", n, ilo, &scale[1], m, &us[
			us_offset], ldus, &us[*n + 1 + us_dim1], ldus, &ierr, 
			(ftnlen)1, (ftnlen)8);
#line 523 "MB03ZD.f"
	    }
#line 526 "MB03ZD.f"
	    if (luu) {
#line 526 "MB03ZD.f"
		mb04di_(balanc, "Positive", n, ilo, &scale[1], m, &uu[
			uu_offset], lduu, &uu[*n + 1 + uu_dim1], lduu, &ierr, 
			(ftnlen)1, (ftnlen)8);
#line 526 "MB03ZD.f"
	    }
#line 529 "MB03ZD.f"
	}

#line 531 "MB03ZD.f"
	if (lus) {
#line 532 "MB03ZD.f"
	    i__1 = *n << 1;
#line 532 "MB03ZD.f"
	    i__2 = *ldwork - *m;
#line 532 "MB03ZD.f"
	    dgeqrf_(&i__1, m, &us[us_offset], ldus, &dwork[1], &dwork[*m + 1],
		     &i__2, &ierr);
/* Computing MAX */
#line 534 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[*m + 1] + *m;
#line 534 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);
#line 535 "MB03ZD.f"
	    i__1 = *n << 1;
#line 535 "MB03ZD.f"
	    i__2 = *ldwork - *m;
#line 535 "MB03ZD.f"
	    dorgqr_(&i__1, m, m, &us[us_offset], ldus, &dwork[1], &dwork[*m + 
		    1], &i__2, &ierr);
/* Computing MAX */
#line 537 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[*m + 1] + *m;
#line 537 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);
#line 538 "MB03ZD.f"
	}
#line 539 "MB03ZD.f"
	if (luu) {
#line 540 "MB03ZD.f"
	    i__1 = *n << 1;
#line 540 "MB03ZD.f"
	    i__2 = *ldwork - *m;
#line 540 "MB03ZD.f"
	    dgeqrf_(&i__1, m, &uu[uu_offset], lduu, &dwork[1], &dwork[*m + 1],
		     &i__2, &ierr);
/* Computing MAX */
#line 542 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[*m + 1] + *m;
#line 542 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);
#line 543 "MB03ZD.f"
	    i__1 = *n << 1;
#line 543 "MB03ZD.f"
	    i__2 = *ldwork - *m;
#line 543 "MB03ZD.f"
	    dorgqr_(&i__1, m, m, &uu[uu_offset], lduu, &dwork[1], &dwork[*m + 
		    1], &i__2, &ierr);
/* Computing MAX */
#line 545 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[*m + 1] + *m;
#line 545 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);
#line 546 "MB03ZD.f"
	}

#line 548 "MB03ZD.f"
	if (lbal && ! lbef) {
#line 549 "MB03ZD.f"
	    if (lus) {
#line 549 "MB03ZD.f"
		mb04di_(balanc, "Positive", n, ilo, &scale[1], m, &us[
			us_offset], ldus, &us[*n + 1 + us_dim1], ldus, &ierr, 
			(ftnlen)1, (ftnlen)8);
#line 549 "MB03ZD.f"
	    }
#line 552 "MB03ZD.f"
	    if (luu) {
#line 552 "MB03ZD.f"
		mb04di_(balanc, "Positive", n, ilo, &scale[1], m, &uu[
			uu_offset], lduu, &uu[*n + 1 + uu_dim1], lduu, &ierr, 
			(ftnlen)1, (ftnlen)8);
#line 552 "MB03ZD.f"
	    }
#line 555 "MB03ZD.f"
	}

#line 557 "MB03ZD.f"
    } else {

#line 559 "MB03ZD.f"
	i__1 = *n << 1;
#line 559 "MB03ZD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 560 "MB03ZD.f"
	    lwork[i__] = TRUE_;
#line 561 "MB03ZD.f"
/* L60: */
#line 561 "MB03ZD.f"
	}

#line 563 "MB03ZD.f"
	if (lus && ! luu) {

/*           Workspace requirements: MAX( 2*N*N + 2*N, 8*N ) */

#line 567 "MB03ZD.f"
	    mb03za_("Update", "Update", "Update", "Init", which, &select[1], 
		    n, &s[s_offset], lds, &t[t_offset], ldt, &g[g_offset], 
		    ldg, &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1[
		    v1_offset], ldv1, &v2[v2_offset], ldv2, &us[us_offset], 
		    ldus, &wr[1], &wi[1], m, &dwork[1], ldwork, &ierr, (
		    ftnlen)6, (ftnlen)6, (ftnlen)6, (ftnlen)4, (ftnlen)1);
#line 571 "MB03ZD.f"
	    if (ierr != 0) {
#line 571 "MB03ZD.f"
		goto L250;
#line 571 "MB03ZD.f"
	    }

#line 574 "MB03ZD.f"
	    mb01ux_("Left", "Lower", "Transpose", n, n, &c_b26, &us[*n + 1 + (
		    *n + 1) * us_dim1], ldus, &g[g_offset], ldg, &dwork[1], 
		    ldwork, &ierr, (ftnlen)4, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
#line 577 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 577 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);

#line 579 "MB03ZD.f"
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[(*n + 
		    1) * us_dim1 + 1], ldus, &g[g_offset], ldg, &dwork[1], 
		    ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
#line 582 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 582 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);

#line 584 "MB03ZD.f"
	    i__1 = *n;
#line 584 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 585 "MB03ZD.f"
		daxpy_(&j, &c_b26, &g[j + g_dim1], ldg, &g[j * g_dim1 + 1], &
			c__1);
#line 586 "MB03ZD.f"
/* L70: */
#line 586 "MB03ZD.f"
	    }
#line 587 "MB03ZD.f"
	    pdw = (*n << 1) * *n + 1;

/*           DW <- -[V1;V2]*W11 */

#line 591 "MB03ZD.f"
	    i__1 = *n << 1;
#line 591 "MB03ZD.f"
	    dlacpy_("All", n, n, &v1[v1_offset], ldv1, &dwork[1], &i__1, (
		    ftnlen)3);
#line 592 "MB03ZD.f"
	    i__1 = *n << 1;
#line 592 "MB03ZD.f"
	    dlacpy_("All", n, n, &v2[v2_offset], ldv2, &dwork[*n + 1], &i__1, 
		    (ftnlen)3);
#line 593 "MB03ZD.f"
	    i__1 = *n << 1;
#line 593 "MB03ZD.f"
	    i__2 = *n << 1;
#line 593 "MB03ZD.f"
	    i__3 = *ldwork - pdw + 1;
#line 593 "MB03ZD.f"
	    mb01ux_("Right", "Upper", "No Transpose", &i__1, n, &c_b34, &us[
		    us_offset], ldus, &dwork[1], &i__2, &dwork[pdw], &i__3, &
		    ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
#line 596 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
#line 596 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);

/*           DW2 <- DW2 - U2*W21 */

#line 600 "MB03ZD.f"
	    dlacpy_("All", n, n, &u2[u2_offset], ldu2, &us[us_offset], ldus, (
		    ftnlen)3);
#line 601 "MB03ZD.f"
	    i__1 = *ldwork - pdw + 1;
#line 601 "MB03ZD.f"
	    mb01ux_("Right", "Upper", "No Transpose", n, n, &c_b26, &us[*n + 
		    1 + us_dim1], ldus, &us[us_offset], ldus, &dwork[pdw], &
		    i__1, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
#line 604 "MB03ZD.f"
	    i__1 = *n;
#line 604 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 605 "MB03ZD.f"
		daxpy_(n, &c_b26, &us[j * us_dim1 + 1], &c__1, &dwork[*n + (j 
			- 1 << 1) * *n + 1], &c__1);
#line 606 "MB03ZD.f"
/* L80: */
#line 606 "MB03ZD.f"
	    }

/*           US11 <- -U1*W21 - DW1 */

#line 610 "MB03ZD.f"
	    dlacpy_("All", n, n, &u1[u1_offset], ldu1, &us[us_offset], ldus, (
		    ftnlen)3);
#line 611 "MB03ZD.f"
	    i__1 = *ldwork - pdw + 1;
#line 611 "MB03ZD.f"
	    mb01ux_("Right", "Upper", "No Transpose", n, n, &c_b34, &us[*n + 
		    1 + us_dim1], ldus, &us[us_offset], ldus, &dwork[pdw], &
		    i__1, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
#line 614 "MB03ZD.f"
	    i__1 = *n;
#line 614 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 615 "MB03ZD.f"
		daxpy_(n, &c_b34, &dwork[(j - 1 << 1) * *n + 1], &c__1, &us[j 
			* us_dim1 + 1], &c__1);
#line 616 "MB03ZD.f"
/* L90: */
#line 616 "MB03ZD.f"
	    }

/*           US21 <- DW2 */

#line 620 "MB03ZD.f"
	    i__1 = *n << 1;
#line 620 "MB03ZD.f"
	    dlacpy_("All", n, n, &dwork[*n + 1], &i__1, &us[*n + 1 + us_dim1],
		     ldus, (ftnlen)3);

#line 622 "MB03ZD.f"
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[(*n + 
		    1) * us_dim1 + 1], ldus, &v1[v1_offset], ldv1, &dwork[1], 
		    ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
#line 625 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 625 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);
#line 626 "MB03ZD.f"
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[(*n + 
		    1) * us_dim1 + 1], ldus, &v2[v2_offset], ldv2, &dwork[1], 
		    ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
#line 629 "MB03ZD.f"
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[*n + 
		    1 + (*n + 1) * us_dim1], ldus, &u1[u1_offset], ldu1, &
		    dwork[1], ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12)
		    ;
#line 632 "MB03ZD.f"
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[*n + 
		    1 + (*n + 1) * us_dim1], ldus, &u2[u2_offset], ldu2, &
		    dwork[1], ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12)
		    ;
#line 635 "MB03ZD.f"
	    dlacpy_("All", n, n, &v1[v1_offset], ldv1, &us[(*n + 1) * us_dim1 
		    + 1], ldus, (ftnlen)3);
#line 636 "MB03ZD.f"
	    dlacpy_("All", n, n, &v2[v2_offset], ldv2, &us[*n + 1 + (*n + 1) *
		     us_dim1], ldus, (ftnlen)3);
#line 637 "MB03ZD.f"
	    i__1 = *n;
#line 637 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 638 "MB03ZD.f"
		daxpy_(n, &c_b34, &u1[j * u1_dim1 + 1], &c__1, &us[(*n + j) * 
			us_dim1 + 1], &c__1);
#line 639 "MB03ZD.f"
/* L100: */
#line 639 "MB03ZD.f"
	    }
#line 640 "MB03ZD.f"
	    i__1 = *n;
#line 640 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 641 "MB03ZD.f"
		daxpy_(n, &c_b34, &u2[j * u2_dim1 + 1], &c__1, &us[*n + 1 + (*
			n + j) * us_dim1], &c__1);
#line 642 "MB03ZD.f"
/* L110: */
#line 642 "MB03ZD.f"
	    }

#line 644 "MB03ZD.f"
	    mb03td_("Hamiltonian", "Update", &lwork[1], &lwork[*n + 1], n, &s[
		    s_offset], lds, &g[g_offset], ldg, &us[(*n + 1) * us_dim1 
		    + 1], ldus, &us[*n + 1 + (*n + 1) * us_dim1], ldus, &wr[1]
		    , &wi[1], m, &dwork[1], ldwork, &ierr, (ftnlen)11, (
		    ftnlen)6);
#line 647 "MB03ZD.f"
	    if (ierr != 0) {
#line 648 "MB03ZD.f"
		*info = 4;
#line 649 "MB03ZD.f"
		return 0;
#line 650 "MB03ZD.f"
	    }
#line 651 "MB03ZD.f"
	    dlascl_("General", &c__0, &c__0, &c_b26, &c_b34, n, n, &us[*n + 1 
		    + (*n + 1) * us_dim1], ldus, &ierr, (ftnlen)7);

#line 654 "MB03ZD.f"
	} else if (! lus && luu) {

/*           Workspace requirements: MAX( 2*N*N + 2*N, 8*N ) */

#line 658 "MB03ZD.f"
	    mb03za_("Update", "Update", "Update", "Init", which, &select[1], 
		    n, &s[s_offset], lds, &t[t_offset], ldt, &g[g_offset], 
		    ldg, &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1[
		    v1_offset], ldv1, &v2[v2_offset], ldv2, &uu[uu_offset], 
		    lduu, &wr[1], &wi[1], m, &dwork[1], ldwork, &ierr, (
		    ftnlen)6, (ftnlen)6, (ftnlen)6, (ftnlen)4, (ftnlen)1);
#line 662 "MB03ZD.f"
	    if (ierr != 0) {
#line 662 "MB03ZD.f"
		goto L250;
#line 662 "MB03ZD.f"
	    }
#line 664 "MB03ZD.f"
	    mb01ux_("Left", "Lower", "Transpose", n, n, &c_b26, &uu[*n + 1 + (
		    *n + 1) * uu_dim1], lduu, &g[g_offset], ldg, &dwork[1], 
		    ldwork, &ierr, (ftnlen)4, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
#line 667 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 667 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);
#line 668 "MB03ZD.f"
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &uu[(*n + 
		    1) * uu_dim1 + 1], lduu, &g[g_offset], ldg, &dwork[1], 
		    ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
#line 671 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 671 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);
#line 672 "MB03ZD.f"
	    i__1 = *n;
#line 672 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 673 "MB03ZD.f"
		daxpy_(&j, &c_b26, &g[j + g_dim1], ldg, &g[j * g_dim1 + 1], &
			c__1);
#line 674 "MB03ZD.f"
/* L120: */
#line 674 "MB03ZD.f"
	    }
#line 675 "MB03ZD.f"
	    pdw = (*n << 1) * *n + 1;

/*           DW <- -[V1;V2]*W11 */

#line 679 "MB03ZD.f"
	    i__1 = *n << 1;
#line 679 "MB03ZD.f"
	    dlacpy_("All", n, n, &v1[v1_offset], ldv1, &dwork[1], &i__1, (
		    ftnlen)3);
#line 680 "MB03ZD.f"
	    i__1 = *n << 1;
#line 680 "MB03ZD.f"
	    dlacpy_("All", n, n, &v2[v2_offset], ldv2, &dwork[*n + 1], &i__1, 
		    (ftnlen)3);
#line 681 "MB03ZD.f"
	    i__1 = *n << 1;
#line 681 "MB03ZD.f"
	    i__2 = *n << 1;
#line 681 "MB03ZD.f"
	    i__3 = *ldwork - pdw + 1;
#line 681 "MB03ZD.f"
	    mb01ux_("Right", "Upper", "No Transpose", &i__1, n, &c_b34, &uu[
		    uu_offset], lduu, &dwork[1], &i__2, &dwork[pdw], &i__3, &
		    ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
#line 684 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
#line 684 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);

/*           DW2 <- DW2 - U2*W21 */

#line 688 "MB03ZD.f"
	    dlacpy_("All", n, n, &u2[u2_offset], ldu2, &uu[uu_offset], lduu, (
		    ftnlen)3);
#line 689 "MB03ZD.f"
	    i__1 = *ldwork - pdw + 1;
#line 689 "MB03ZD.f"
	    mb01ux_("Right", "Upper", "No Transpose", n, n, &c_b34, &uu[*n + 
		    1 + uu_dim1], lduu, &uu[uu_offset], lduu, &dwork[pdw], &
		    i__1, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
#line 692 "MB03ZD.f"
	    i__1 = *n;
#line 692 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 693 "MB03ZD.f"
		daxpy_(n, &c_b26, &uu[j * uu_dim1 + 1], &c__1, &dwork[*n + (j 
			- 1 << 1) * *n + 1], &c__1);
#line 694 "MB03ZD.f"
/* L130: */
#line 694 "MB03ZD.f"
	    }

/*           UU11 <- U1*W21 - DW1 */

#line 698 "MB03ZD.f"
	    dlacpy_("All", n, n, &u1[u1_offset], ldu1, &uu[uu_offset], lduu, (
		    ftnlen)3);
#line 699 "MB03ZD.f"
	    i__1 = *ldwork - pdw + 1;
#line 699 "MB03ZD.f"
	    mb01ux_("Right", "Upper", "No Transpose", n, n, &c_b26, &uu[*n + 
		    1 + uu_dim1], lduu, &uu[uu_offset], lduu, &dwork[pdw], &
		    i__1, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
#line 702 "MB03ZD.f"
	    i__1 = *n;
#line 702 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 703 "MB03ZD.f"
		daxpy_(n, &c_b34, &dwork[(j - 1 << 1) * *n + 1], &c__1, &uu[j 
			* uu_dim1 + 1], &c__1);
#line 704 "MB03ZD.f"
/* L140: */
#line 704 "MB03ZD.f"
	    }

/*           UU21 <- DW2 */

#line 708 "MB03ZD.f"
	    i__1 = *n << 1;
#line 708 "MB03ZD.f"
	    dlacpy_("All", n, n, &dwork[*n + 1], &i__1, &uu[*n + 1 + uu_dim1],
		     lduu, (ftnlen)3);

#line 710 "MB03ZD.f"
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &uu[(*n + 
		    1) * uu_dim1 + 1], lduu, &v1[v1_offset], ldv1, &dwork[1], 
		    ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
#line 713 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 713 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);
#line 714 "MB03ZD.f"
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &uu[(*n + 
		    1) * uu_dim1 + 1], lduu, &v2[v2_offset], ldv2, &dwork[1], 
		    ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
#line 717 "MB03ZD.f"
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &uu[*n + 
		    1 + (*n + 1) * uu_dim1], lduu, &u1[u1_offset], ldu1, &
		    dwork[1], ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12)
		    ;
#line 720 "MB03ZD.f"
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &uu[*n + 
		    1 + (*n + 1) * uu_dim1], lduu, &u2[u2_offset], ldu2, &
		    dwork[1], ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12)
		    ;
#line 723 "MB03ZD.f"
	    dlacpy_("All", n, n, &v1[v1_offset], ldv1, &uu[(*n + 1) * uu_dim1 
		    + 1], lduu, (ftnlen)3);
#line 724 "MB03ZD.f"
	    dlacpy_("All", n, n, &v2[v2_offset], ldv2, &uu[*n + 1 + (*n + 1) *
		     uu_dim1], lduu, (ftnlen)3);
#line 725 "MB03ZD.f"
	    i__1 = *n;
#line 725 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 726 "MB03ZD.f"
		daxpy_(n, &c_b26, &u1[j * u1_dim1 + 1], &c__1, &uu[(*n + j) * 
			uu_dim1 + 1], &c__1);
#line 727 "MB03ZD.f"
/* L150: */
#line 727 "MB03ZD.f"
	    }
#line 728 "MB03ZD.f"
	    i__1 = *n;
#line 728 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 729 "MB03ZD.f"
		daxpy_(n, &c_b26, &u2[j * u2_dim1 + 1], &c__1, &uu[*n + 1 + (*
			n + j) * uu_dim1], &c__1);
#line 730 "MB03ZD.f"
/* L160: */
#line 730 "MB03ZD.f"
	    }

#line 732 "MB03ZD.f"
	    mb03td_("Hamiltonian", "Update", &lwork[1], &lwork[*n + 1], n, &s[
		    s_offset], lds, &g[g_offset], ldg, &uu[(*n + 1) * uu_dim1 
		    + 1], lduu, &uu[*n + 1 + (*n + 1) * uu_dim1], lduu, &wr[1]
		    , &wi[1], m, &dwork[1], ldwork, &ierr, (ftnlen)11, (
		    ftnlen)6);
#line 735 "MB03ZD.f"
	    if (ierr != 0) {
#line 736 "MB03ZD.f"
		*info = 4;
#line 737 "MB03ZD.f"
		return 0;
#line 738 "MB03ZD.f"
	    }
#line 739 "MB03ZD.f"
	    dlascl_("General", &c__0, &c__0, &c_b26, &c_b34, n, n, &uu[*n + 1 
		    + (*n + 1) * uu_dim1], lduu, &ierr, (ftnlen)7);
#line 741 "MB03ZD.f"
	} else {

/*           Workspace requirements: 8*N */

#line 745 "MB03ZD.f"
	    mb03za_("Update", "Update", "Update", "Init", which, &select[1], 
		    n, &s[s_offset], lds, &t[t_offset], ldt, &g[g_offset], 
		    ldg, &u1[u1_offset], ldu1, &u2[u2_offset], ldu2, &v1[
		    v1_offset], ldv1, &v2[v2_offset], ldv2, &us[us_offset], 
		    ldus, &wr[1], &wi[1], m, &dwork[1], ldwork, &ierr, (
		    ftnlen)6, (ftnlen)6, (ftnlen)6, (ftnlen)4, (ftnlen)1);
#line 749 "MB03ZD.f"
	    if (ierr != 0) {
#line 749 "MB03ZD.f"
		goto L250;
#line 749 "MB03ZD.f"
	    }
#line 751 "MB03ZD.f"
	    mb01ux_("Left", "Lower", "Transpose", n, n, &c_b26, &us[*n + 1 + (
		    *n + 1) * us_dim1], ldus, &g[g_offset], ldg, &dwork[1], 
		    ldwork, &ierr, (ftnlen)4, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
#line 754 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 754 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);
#line 755 "MB03ZD.f"
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[(*n + 
		    1) * us_dim1 + 1], ldus, &g[g_offset], ldg, &dwork[1], 
		    ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
#line 758 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 758 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);
#line 759 "MB03ZD.f"
	    i__1 = *n;
#line 759 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 760 "MB03ZD.f"
		daxpy_(&j, &c_b26, &g[j + g_dim1], ldg, &g[j * g_dim1 + 1], &
			c__1);
#line 761 "MB03ZD.f"
/* L170: */
#line 761 "MB03ZD.f"
	    }

/*           UU = [ V1 -V2; U1 -U2 ]*diag(W11,W21) */

#line 765 "MB03ZD.f"
	    dlacpy_("All", n, n, &v1[v1_offset], ldv1, &uu[uu_offset], lduu, (
		    ftnlen)3);
#line 766 "MB03ZD.f"
	    dlacpy_("All", n, n, &v2[v2_offset], ldv2, &uu[*n + 1 + uu_dim1], 
		    lduu, (ftnlen)3);
#line 767 "MB03ZD.f"
	    i__1 = *n << 1;
#line 767 "MB03ZD.f"
	    mb01ux_("Right", "Upper", "No Transpose", &i__1, n, &c_b26, &us[
		    us_offset], ldus, &uu[uu_offset], lduu, &dwork[1], ldwork,
		     &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
/* Computing MAX */
#line 769 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 769 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);
#line 770 "MB03ZD.f"
	    dlacpy_("All", n, n, &u1[u1_offset], ldu1, &uu[(*n + 1) * uu_dim1 
		    + 1], lduu, (ftnlen)3);
#line 771 "MB03ZD.f"
	    dlacpy_("All", n, n, &u2[u2_offset], ldu2, &uu[*n + 1 + (*n + 1) *
		     uu_dim1], lduu, (ftnlen)3);
#line 772 "MB03ZD.f"
	    i__1 = *n << 1;
#line 772 "MB03ZD.f"
	    mb01ux_("Right", "Upper", "No Transpose", &i__1, n, &c_b26, &us[*
		    n + 1 + us_dim1], ldus, &uu[(*n + 1) * uu_dim1 + 1], lduu,
		     &dwork[1], ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)
		    12);
#line 775 "MB03ZD.f"
	    i__1 = *n << 1;
#line 775 "MB03ZD.f"
	    dlascl_("General", &c__0, &c__0, &c_b26, &c_b34, n, &i__1, &uu[*n 
		    + 1 + uu_dim1], lduu, &ierr, (ftnlen)7);

#line 778 "MB03ZD.f"
	    i__1 = *n << 1;
#line 778 "MB03ZD.f"
	    dlacpy_("All", &i__1, n, &uu[uu_offset], lduu, &us[us_offset], 
		    ldus, (ftnlen)3);
#line 779 "MB03ZD.f"
	    i__1 = *n;
#line 779 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 780 "MB03ZD.f"
		i__2 = *n << 1;
#line 780 "MB03ZD.f"
		daxpy_(&i__2, &c_b34, &uu[(*n + j) * uu_dim1 + 1], &c__1, &us[
			j * us_dim1 + 1], &c__1);
#line 781 "MB03ZD.f"
/* L180: */
#line 781 "MB03ZD.f"
	    }
#line 782 "MB03ZD.f"
	    i__1 = *n;
#line 782 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 783 "MB03ZD.f"
		i__2 = *n << 1;
#line 783 "MB03ZD.f"
		daxpy_(&i__2, &c_b26, &uu[(*n + j) * uu_dim1 + 1], &c__1, &uu[
			j * uu_dim1 + 1], &c__1);
#line 784 "MB03ZD.f"
/* L190: */
#line 784 "MB03ZD.f"
	    }

/*           V1 <- V1*W12-U1*W22 */
/*           U1 <- V1*W12+U1*W22 */
/*           V2 <- V2*W12-U2*W22 */
/*           U2 <- V2*W12+U2*W22 */

#line 791 "MB03ZD.f"
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[(*n + 
		    1) * us_dim1 + 1], ldus, &v1[v1_offset], ldv1, &dwork[1], 
		    ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
#line 794 "MB03ZD.f"
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[(*n + 
		    1) * us_dim1 + 1], ldus, &v2[v2_offset], ldv2, &dwork[1], 
		    ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12);
#line 797 "MB03ZD.f"
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[*n + 
		    1 + (*n + 1) * us_dim1], ldus, &u1[u1_offset], ldu1, &
		    dwork[1], ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12)
		    ;
#line 800 "MB03ZD.f"
	    mb01ux_("Right", "Lower", "No Transpose", n, n, &c_b26, &us[*n + 
		    1 + (*n + 1) * us_dim1], ldus, &u2[u2_offset], ldu2, &
		    dwork[1], ldwork, &ierr, (ftnlen)5, (ftnlen)5, (ftnlen)12)
		    ;
#line 803 "MB03ZD.f"
	    i__1 = *n;
#line 803 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 804 "MB03ZD.f"
		i__2 = *n;
#line 804 "MB03ZD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 805 "MB03ZD.f"
		    temp = v1[i__ + j * v1_dim1];
#line 806 "MB03ZD.f"
		    v1[i__ + j * v1_dim1] = temp - u1[i__ + j * u1_dim1];
#line 807 "MB03ZD.f"
		    u1[i__ + j * u1_dim1] = temp + u1[i__ + j * u1_dim1];
#line 808 "MB03ZD.f"
/* L200: */
#line 808 "MB03ZD.f"
		}
#line 809 "MB03ZD.f"
/* L210: */
#line 809 "MB03ZD.f"
	    }
#line 810 "MB03ZD.f"
	    i__1 = *n;
#line 810 "MB03ZD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 811 "MB03ZD.f"
		i__2 = *n;
#line 811 "MB03ZD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 812 "MB03ZD.f"
		    temp = v2[i__ + j * v2_dim1];
#line 813 "MB03ZD.f"
		    v2[i__ + j * v2_dim1] = temp - u2[i__ + j * u2_dim1];
#line 814 "MB03ZD.f"
		    u2[i__ + j * u2_dim1] = temp + u2[i__ + j * u2_dim1];
#line 815 "MB03ZD.f"
/* L220: */
#line 815 "MB03ZD.f"
		}
#line 816 "MB03ZD.f"
/* L230: */
#line 816 "MB03ZD.f"
	    }

#line 818 "MB03ZD.f"
	    i__1 = *n << 1;
#line 818 "MB03ZD.f"
	    dlaset_("All", &i__1, n, &c_b272, &c_b26, &us[(*n + 1) * us_dim1 
		    + 1], ldus, (ftnlen)3);
#line 819 "MB03ZD.f"
	    mb03td_("Hamiltonian", "Update", &lwork[1], &lwork[*n + 1], n, &s[
		    s_offset], lds, &g[g_offset], ldg, &us[(*n + 1) * us_dim1 
		    + 1], ldus, &us[*n + 1 + (*n + 1) * us_dim1], ldus, &wr[1]
		    , &wi[1], m, &dwork[1], ldwork, &ierr, (ftnlen)11, (
		    ftnlen)6);
#line 822 "MB03ZD.f"
	    if (ierr != 0) {
#line 823 "MB03ZD.f"
		*info = 4;
#line 824 "MB03ZD.f"
		return 0;
#line 825 "MB03ZD.f"
	    }
#line 826 "MB03ZD.f"
	    dgemm_("No Transpose", "No Transpose", n, n, n, &c_b26, &u1[
		    u1_offset], ldu1, &us[(*n + 1) * us_dim1 + 1], ldus, &
		    c_b272, &uu[(*n + 1) * uu_dim1 + 1], lduu, (ftnlen)12, (
		    ftnlen)12);
#line 829 "MB03ZD.f"
	    dgemm_("No Transpose", "No Transpose", n, n, n, &c_b34, &u2[
		    u2_offset], ldu2, &us[*n + 1 + (*n + 1) * us_dim1], ldus, 
		    &c_b26, &uu[(*n + 1) * uu_dim1 + 1], lduu, (ftnlen)12, (
		    ftnlen)12);
#line 832 "MB03ZD.f"
	    dgemm_("No Transpose", "No Transpose", n, n, n, &c_b34, &u1[
		    u1_offset], ldu1, &us[*n + 1 + (*n + 1) * us_dim1], ldus, 
		    &c_b272, &uu[*n + 1 + (*n + 1) * uu_dim1], lduu, (ftnlen)
		    12, (ftnlen)12);
#line 835 "MB03ZD.f"
	    dgemm_("No Transpose", "No Transpose", n, n, n, &c_b34, &u2[
		    u2_offset], ldu2, &us[(*n + 1) * us_dim1 + 1], ldus, &
		    c_b26, &uu[*n + 1 + (*n + 1) * uu_dim1], lduu, (ftnlen)12,
		     (ftnlen)12);
#line 838 "MB03ZD.f"
	    dlacpy_("All", n, n, &us[(*n + 1) * us_dim1 + 1], ldus, &u1[
		    u1_offset], ldu1, (ftnlen)3);
#line 839 "MB03ZD.f"
	    dlacpy_("All", n, n, &us[*n + 1 + (*n + 1) * us_dim1], ldus, &u2[
		    u2_offset], ldu2, (ftnlen)3);
#line 840 "MB03ZD.f"
	    dgemm_("No Transpose", "No Transpose", n, n, n, &c_b26, &v1[
		    v1_offset], ldv1, &u1[u1_offset], ldu1, &c_b272, &us[(*n 
		    + 1) * us_dim1 + 1], ldus, (ftnlen)12, (ftnlen)12);
#line 842 "MB03ZD.f"
	    dgemm_("No Transpose", "No Transpose", n, n, n, &c_b34, &v2[
		    v2_offset], ldv2, &u2[u2_offset], ldu2, &c_b26, &us[(*n + 
		    1) * us_dim1 + 1], ldus, (ftnlen)12, (ftnlen)12);
#line 844 "MB03ZD.f"
	    dgemm_("No Transpose", "No Transpose", n, n, n, &c_b34, &v1[
		    v1_offset], ldv1, &u2[u2_offset], ldu2, &c_b272, &us[*n + 
		    1 + (*n + 1) * us_dim1], ldus, (ftnlen)12, (ftnlen)12);
#line 846 "MB03ZD.f"
	    dgemm_("No Transpose", "No Transpose", n, n, n, &c_b34, &v2[
		    v2_offset], ldv2, &u1[u1_offset], ldu1, &c_b26, &us[*n + 
		    1 + (*n + 1) * us_dim1], ldus, (ftnlen)12, (ftnlen)12);
#line 848 "MB03ZD.f"
	}

/*        Orthonormalize obtained bases and apply inverse balancing */
/*        transformation. */

#line 853 "MB03ZD.f"
	if (lbal && lbef) {
#line 854 "MB03ZD.f"
	    if (lus) {
#line 854 "MB03ZD.f"
		mb04di_(balanc, "Positive", n, ilo, &scale[1], n, &us[
			us_offset], ldus, &us[*n + 1 + us_dim1], ldus, &ierr, 
			(ftnlen)1, (ftnlen)8);
#line 854 "MB03ZD.f"
	    }
#line 857 "MB03ZD.f"
	    if (luu) {
#line 857 "MB03ZD.f"
		mb04di_(balanc, "Positive", n, ilo, &scale[1], n, &uu[
			uu_offset], lduu, &uu[*n + 1 + uu_dim1], lduu, &ierr, 
			(ftnlen)1, (ftnlen)8);
#line 857 "MB03ZD.f"
	    }
#line 860 "MB03ZD.f"
	}

/*        Workspace requirements: 8*N+1 */

#line 864 "MB03ZD.f"
	i__1 = *n << 1;
#line 864 "MB03ZD.f"
	for (j = 1; j <= i__1; ++j) {
#line 865 "MB03ZD.f"
	    iwork[j] = 0;
#line 866 "MB03ZD.f"
/* L240: */
#line 866 "MB03ZD.f"
	}
#line 867 "MB03ZD.f"
	if (lus) {
#line 868 "MB03ZD.f"
	    i__1 = *n << 1;
#line 868 "MB03ZD.f"
	    i__2 = *n << 1;
#line 868 "MB03ZD.f"
	    i__3 = *ldwork - (*n << 1);
#line 868 "MB03ZD.f"
	    dgeqp3_(&i__1, &i__2, &us[us_offset], ldus, &iwork[1], &dwork[1], 
		    &dwork[(*n << 1) + 1], &i__3, &ierr);
/* Computing MAX */
#line 870 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[(*n << 1) + 1] + (*n << 1);
#line 870 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);
#line 871 "MB03ZD.f"
	    i__1 = *n << 1;
#line 871 "MB03ZD.f"
	    i__2 = *n << 1;
#line 871 "MB03ZD.f"
	    i__3 = *ldwork - (*n << 1);
#line 871 "MB03ZD.f"
	    dorgqr_(&i__1, &i__2, n, &us[us_offset], ldus, &dwork[1], &dwork[(
		    *n << 1) + 1], &i__3, &ierr);
/* Computing MAX */
#line 873 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[(*n << 1) + 1] + (*n << 1);
#line 873 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);
#line 874 "MB03ZD.f"
	}
#line 875 "MB03ZD.f"
	if (luu) {
#line 876 "MB03ZD.f"
	    i__1 = *n << 1;
#line 876 "MB03ZD.f"
	    i__2 = *n << 1;
#line 876 "MB03ZD.f"
	    i__3 = *ldwork - (*n << 1);
#line 876 "MB03ZD.f"
	    dgeqp3_(&i__1, &i__2, &uu[uu_offset], lduu, &iwork[1], &dwork[1], 
		    &dwork[(*n << 1) + 1], &i__3, &ierr);
/* Computing MAX */
#line 878 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[(*n << 1) + 1] + (*n << 1);
#line 878 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);
#line 879 "MB03ZD.f"
	    i__1 = *n << 1;
#line 879 "MB03ZD.f"
	    i__2 = *n << 1;
#line 879 "MB03ZD.f"
	    i__3 = *ldwork - (*n << 1);
#line 879 "MB03ZD.f"
	    dorgqr_(&i__1, &i__2, n, &uu[uu_offset], lduu, &dwork[1], &dwork[(
		    *n << 1) + 1], &i__3, &ierr);
/* Computing MAX */
#line 881 "MB03ZD.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[(*n << 1) + 1] + (*n << 1);
#line 881 "MB03ZD.f"
	    wrkopt = max(i__1,i__2);
#line 882 "MB03ZD.f"
	}

#line 884 "MB03ZD.f"
	if (lbal && ! lbef) {
#line 885 "MB03ZD.f"
	    if (lus) {
#line 885 "MB03ZD.f"
		mb04di_(balanc, "Positive", n, ilo, &scale[1], n, &us[
			us_offset], ldus, &us[*n + 1 + us_dim1], ldus, &ierr, 
			(ftnlen)1, (ftnlen)8);
#line 885 "MB03ZD.f"
	    }
#line 888 "MB03ZD.f"
	    if (luu) {
#line 888 "MB03ZD.f"
		mb04di_(balanc, "Positive", n, ilo, &scale[1], n, &uu[
			uu_offset], lduu, &uu[*n + 1 + uu_dim1], lduu, &ierr, 
			(ftnlen)1, (ftnlen)8);
#line 888 "MB03ZD.f"
	    }
#line 891 "MB03ZD.f"
	}
#line 892 "MB03ZD.f"
    }

#line 894 "MB03ZD.f"
    dscal_(m, &c_b34, &wr[1], &c__1);
#line 895 "MB03ZD.f"
    dwork[1] = (doublereal) wrkopt;

#line 897 "MB03ZD.f"
    return 0;
#line 898 "MB03ZD.f"
L250:
#line 899 "MB03ZD.f"
    if (ierr == 1) {
#line 900 "MB03ZD.f"
	*info = 2;
#line 901 "MB03ZD.f"
    } else if (ierr == 2 || ierr == 4) {
#line 902 "MB03ZD.f"
	*info = 1;
#line 903 "MB03ZD.f"
    } else if (ierr == 3) {
#line 904 "MB03ZD.f"
	*info = 3;
#line 905 "MB03ZD.f"
    }
#line 906 "MB03ZD.f"
    return 0;
/* *** Last line of MB03ZD *** */
} /* mb03zd_ */

