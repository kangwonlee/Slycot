#line 1 "TB01ZD.f"
/* TB01ZD.f -- translated by f2c (version 20100827).
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

#line 1 "TB01ZD.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b10 = 0.;
static doublereal c_b18 = 1.;
static integer c__0 = 0;

/* Subroutine */ int tb01zd_(char *jobz, integer *n, integer *p, doublereal *
	a, integer *lda, doublereal *b, doublereal *c__, integer *ldc, 
	integer *ncont, doublereal *z__, integer *ldz, doublereal *tau, 
	doublereal *tol, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen jobz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal h__;
    static integer j;
    static doublereal b1, nblk[1];
    static integer itau;
    extern /* Subroutine */ int mb01pd_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen), dlarf_(char *
	    , integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, ftnlen);
    static logical ljobf, ljobi;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal anorm, bnorm;
    static logical ljobz;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgehrd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dlarfg_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *), dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static doublereal fanorm, fbnorm;
    extern /* Subroutine */ int dormhr_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static doublereal thresh;
    extern /* Subroutine */ int dorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal wrkopt;


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

/*     To find a controllable realization for the linear time-invariant */
/*     single-input system */

/*             dX/dt = A * X + B * U, */
/*                Y  = C * X, */

/*     where A is an N-by-N matrix, B is an N element vector, C is an */
/*     P-by-N matrix, and A and B are reduced by this routine to */
/*     orthogonal canonical form using (and optionally accumulating) */
/*     orthogonal similarity transformations, which are also applied */
/*     to C. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBZ    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix Z the orthogonal similarity transformations for */
/*             reducing the system, as follows: */
/*             = 'N':  Do not form Z and do not store the orthogonal */
/*                     transformations; */
/*             = 'F':  Do not form Z, but store the orthogonal */
/*                     transformations in the factored form; */
/*             = 'I':  Z is initialized to the unit matrix and the */
/*                     orthogonal transformation matrix Z is returned. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the original state-space representation, */
/*             i.e. the order of the matrix A.  N >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs, or of rows of C.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original state dynamics matrix A. */
/*             On exit, the leading NCONT-by-NCONT upper Hessenberg */
/*             part of this array contains the canonical form of the */
/*             state dynamics matrix, given by Z' * A * Z, of a */
/*             controllable realization for the original system. The */
/*             elements below the first subdiagonal are set to zero. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, the original input/state vector B. */
/*             On exit, the leading NCONT elements of this array contain */
/*             canonical form of the input/state vector, given by Z' * B, */
/*             with all elements but B(1) set to zero. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the output/state matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the transformed output/state matrix, given by C * Z, and */
/*             the leading P-by-NCONT part contains the output/state */
/*             matrix of the controllable realization. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     NCONT   (output) INTEGER */
/*             The order of the controllable state-space representation. */

/*     Z       (output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             If JOBZ = 'I', then the leading N-by-N part of this array */
/*             contains the matrix of accumulated orthogonal similarity */
/*             transformations which reduces the given system to */
/*             orthogonal canonical form. */
/*             If JOBZ = 'F', the elements below the diagonal, with the */
/*             array TAU, represent the orthogonal transformation matrix */
/*             as a product of elementary reflectors. The transformation */
/*             matrix can then be obtained by calling the LAPACK Library */
/*             routine DORGQR. */
/*             If JOBZ = 'N', the array Z is not referenced and can be */
/*             supplied as a dummy array (i.e. set parameter LDZ = 1 and */
/*             declare this array to be Z(1,1) in the calling program). */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z. If JOBZ = 'I' or */
/*             JOBZ = 'F', LDZ >= MAX(1,N); if JOBZ = 'N', LDZ >= 1. */

/*     TAU     (output) DOUBLE PRECISION array, dimension (N) */
/*             The elements of TAU contain the scalar factors of the */
/*             elementary reflectors used in the reduction of B and A. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in determining the */
/*             controllability of (A,B). If the user sets TOL > 0, then */
/*             the given value of TOL is used as an absolute tolerance; */
/*             elements with absolute value less than TOL are considered */
/*             neglijible. If the user sets TOL <= 0, then an implicitly */
/*             computed, default tolerance, defined by */
/*             TOLDEF = N*EPS*MAX( NORM(A), NORM(B) ) is used instead, */
/*             where EPS is the machine precision (see LAPACK Library */
/*             routine DLAMCH). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. LDWORK >= MAX(1,N,P). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The Householder matrix which reduces all but the first element */
/*     of vector B to zero is found and this orthogonal similarity */
/*     transformation is applied to the matrix A. The resulting A is then */
/*     reduced to upper Hessenberg form by a sequence of Householder */
/*     transformations. Finally, the order of the controllable state- */
/*     space representation (NCONT) is determined by finding the position */
/*     of the first sub-diagonal element of A which is below an */
/*     appropriate zero threshold, either TOL or TOLDEF (see parameter */
/*     TOL); if NORM(B) is smaller than this threshold, NCONT is set to */
/*     zero, and no computations for reducing the system to orthogonal */
/*     canonical form are performed. */
/*     All orthogonal transformations determined in this process are also */
/*     applied to the matrix C, from the right. */

/*     REFERENCES */

/*     [1] Konstantinov, M.M., Petkov, P.Hr. and Christov, N.D. */
/*         Orthogonal Invariants and Canonical Forms for Linear */
/*         Controllable Systems. */
/*         Proc. 8th IFAC World Congress, Kyoto, 1, pp. 49-54, 1981. */

/*     [2] Hammarling, S.J. */
/*         Notes on the use of orthogonal similarity transformations in */
/*         control. */
/*         NPL Report DITC 8/82, August 1982. */

/*     [3] Paige, C.C */
/*         Properties of numerical algorithms related to computing */
/*         controllability. */
/*         IEEE Trans. Auto. Contr., AC-26, pp. 130-138, 1981. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations and is backward stable. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1998. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2001, */
/*     Sept. 2003. */

/*     KEYWORDS */

/*     Controllability, minimal realization, orthogonal canonical form, */
/*     orthogonal transformation. */

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

#line 225 "TB01ZD.f"
    /* Parameter adjustments */
#line 225 "TB01ZD.f"
    a_dim1 = *lda;
#line 225 "TB01ZD.f"
    a_offset = 1 + a_dim1;
#line 225 "TB01ZD.f"
    a -= a_offset;
#line 225 "TB01ZD.f"
    --b;
#line 225 "TB01ZD.f"
    c_dim1 = *ldc;
#line 225 "TB01ZD.f"
    c_offset = 1 + c_dim1;
#line 225 "TB01ZD.f"
    c__ -= c_offset;
#line 225 "TB01ZD.f"
    z_dim1 = *ldz;
#line 225 "TB01ZD.f"
    z_offset = 1 + z_dim1;
#line 225 "TB01ZD.f"
    z__ -= z_offset;
#line 225 "TB01ZD.f"
    --tau;
#line 225 "TB01ZD.f"
    --dwork;
#line 225 "TB01ZD.f"

#line 225 "TB01ZD.f"
    /* Function Body */
#line 225 "TB01ZD.f"
    *info = 0;
#line 226 "TB01ZD.f"
    ljobf = lsame_(jobz, "F", (ftnlen)1, (ftnlen)1);
#line 227 "TB01ZD.f"
    ljobi = lsame_(jobz, "I", (ftnlen)1, (ftnlen)1);
#line 228 "TB01ZD.f"
    ljobz = ljobf || ljobi;

/*     Test the input scalar arguments. */

#line 232 "TB01ZD.f"
    if (! ljobz && ! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 233 "TB01ZD.f"
	*info = -1;
#line 234 "TB01ZD.f"
    } else if (*n < 0) {
#line 235 "TB01ZD.f"
	*info = -2;
#line 236 "TB01ZD.f"
    } else if (*p < 0) {
#line 237 "TB01ZD.f"
	*info = -3;
#line 238 "TB01ZD.f"
    } else if (*lda < max(1,*n)) {
#line 239 "TB01ZD.f"
	*info = -5;
#line 240 "TB01ZD.f"
    } else if (*ldc < max(1,*p)) {
#line 241 "TB01ZD.f"
	*info = -8;
#line 242 "TB01ZD.f"
    } else if (*ldz < 1 || ljobz && *ldz < *n) {
#line 243 "TB01ZD.f"
	*info = -11;
#line 244 "TB01ZD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 244 "TB01ZD.f"
	i__1 = max(1,*n);
#line 244 "TB01ZD.f"
	if (*ldwork < max(i__1,*p)) {
#line 245 "TB01ZD.f"
	    *info = -15;
#line 246 "TB01ZD.f"
	}
#line 246 "TB01ZD.f"
    }

#line 248 "TB01ZD.f"
    if (*info != 0) {

/*        Error return. */

#line 252 "TB01ZD.f"
	i__1 = -(*info);
#line 252 "TB01ZD.f"
	xerbla_("TB01ZD", &i__1, (ftnlen)6);
#line 253 "TB01ZD.f"
	return 0;
#line 254 "TB01ZD.f"
    }

/*     Quick return if possible. */

#line 258 "TB01ZD.f"
    *ncont = 0;
#line 259 "TB01ZD.f"
    dwork[1] = 1.;
#line 260 "TB01ZD.f"
    if (*n == 0) {
#line 260 "TB01ZD.f"
	return 0;
#line 260 "TB01ZD.f"
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 269 "TB01ZD.f"
    wrkopt = 1.;

/*     Calculate the absolute norms of A and B (used for scaling). */

#line 273 "TB01ZD.f"
    anorm = dlange_("Max", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)3);
#line 274 "TB01ZD.f"
    bnorm = dlange_("Max", n, &c__1, &b[1], n, &dwork[1], (ftnlen)3);

/*     Return if matrix B is zero. */

#line 278 "TB01ZD.f"
    if (bnorm == 0.) {
#line 279 "TB01ZD.f"
	if (ljobf) {
#line 280 "TB01ZD.f"
	    dlaset_("Full", n, n, &c_b10, &c_b10, &z__[z_offset], ldz, (
		    ftnlen)4);
#line 281 "TB01ZD.f"
	    dlaset_("Full", n, &c__1, &c_b10, &c_b10, &tau[1], n, (ftnlen)4);
#line 282 "TB01ZD.f"
	} else if (ljobi) {
#line 283 "TB01ZD.f"
	    dlaset_("Full", n, n, &c_b10, &c_b18, &z__[z_offset], ldz, (
		    ftnlen)4);
#line 284 "TB01ZD.f"
	}
#line 285 "TB01ZD.f"
	return 0;
#line 286 "TB01ZD.f"
    }

/*     Scale (if needed) the matrices A and B. */

#line 290 "TB01ZD.f"
    mb01pd_("S", "G", n, n, &c__0, &c__0, &anorm, &c__0, nblk, &a[a_offset], 
	    lda, info, (ftnlen)1, (ftnlen)1);
#line 291 "TB01ZD.f"
    mb01pd_("S", "G", n, &c__1, &c__0, &c__0, &bnorm, &c__0, nblk, &b[1], n, 
	    info, (ftnlen)1, (ftnlen)1);

/*     Calculate the Frobenius norm of A and the 1-norm of B (used for */
/*     controlability test). */

#line 296 "TB01ZD.f"
    fanorm = dlange_("Frobenius", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)
	    9);
#line 297 "TB01ZD.f"
    fbnorm = dlange_("1-norm", n, &c__1, &b[1], n, &dwork[1], (ftnlen)6);

#line 299 "TB01ZD.f"
    toldef = *tol;
#line 300 "TB01ZD.f"
    if (toldef <= 0.) {

/*        Use the default tolerance in controllability determination. */

#line 304 "TB01ZD.f"
	thresh = (doublereal) (*n) * dlamch_("EPSILON", (ftnlen)7);
#line 305 "TB01ZD.f"
	toldef = thresh * max(fanorm,fbnorm);
#line 306 "TB01ZD.f"
    }

#line 308 "TB01ZD.f"
    itau = 1;
#line 309 "TB01ZD.f"
    if (fbnorm > toldef) {

/*        B is not negligible compared with A. */

#line 313 "TB01ZD.f"
	if (*n > 1) {

/*           Transform B by a Householder matrix Z1: store vector */
/*           describing this temporarily in B and in the local scalar H. */

#line 318 "TB01ZD.f"
	    dlarfg_(n, &b[1], &b[2], &c__1, &h__);

#line 320 "TB01ZD.f"
	    b1 = b[1];
#line 321 "TB01ZD.f"
	    b[1] = 1.;

/*           Form Z1 * A * Z1. */
/*           Workspace: need N. */

#line 326 "TB01ZD.f"
	    dlarf_("Right", n, n, &b[1], &c__1, &h__, &a[a_offset], lda, &
		    dwork[1], (ftnlen)5);
#line 327 "TB01ZD.f"
	    dlarf_("Left", n, n, &b[1], &c__1, &h__, &a[a_offset], lda, &
		    dwork[1], (ftnlen)4);

/*           Form C * Z1. */
/*           Workspace: need P. */

#line 332 "TB01ZD.f"
	    dlarf_("Right", p, n, &b[1], &c__1, &h__, &c__[c_offset], ldc, &
		    dwork[1], (ftnlen)5);

#line 334 "TB01ZD.f"
	    b[1] = b1;
#line 335 "TB01ZD.f"
	    tau[1] = h__;
#line 336 "TB01ZD.f"
	    ++itau;
#line 337 "TB01ZD.f"
	} else {
#line 338 "TB01ZD.f"
	    b1 = b[1];
#line 339 "TB01ZD.f"
	    tau[1] = 0.;
#line 340 "TB01ZD.f"
	}

/*        Reduce modified A to upper Hessenberg form by an orthogonal */
/*        similarity transformation with matrix Z2. */
/*        Workspace: need N;  prefer N*NB. */

#line 346 "TB01ZD.f"
	dgehrd_(n, &c__1, n, &a[a_offset], lda, &tau[itau], &dwork[1], ldwork,
		 info);
#line 347 "TB01ZD.f"
	wrkopt = dwork[1];

/*        Form C * Z2. */
/*        Workspace: need P;  prefer P*NB. */

#line 352 "TB01ZD.f"
	dormhr_("Right", "No transpose", p, n, &c__1, n, &a[a_offset], lda, &
		tau[itau], &c__[c_offset], ldc, &dwork[1], ldwork, info, (
		ftnlen)5, (ftnlen)12);
#line 354 "TB01ZD.f"
	wrkopt = max(wrkopt,dwork[1]);

#line 356 "TB01ZD.f"
	if (ljobz) {

/*           Save the orthogonal transformations used, so that they could */
/*           be accumulated by calling DORGQR routine. */

#line 361 "TB01ZD.f"
	    if (*n > 1) {
#line 361 "TB01ZD.f"
		i__1 = *n - 1;
#line 361 "TB01ZD.f"
		i__2 = *n - 1;
#line 361 "TB01ZD.f"
		dlacpy_("Full", &i__1, &c__1, &b[2], &i__2, &z__[z_dim1 + 2], 
			ldz, (ftnlen)4);
#line 361 "TB01ZD.f"
	    }
#line 363 "TB01ZD.f"
	    if (*n > 2) {
#line 363 "TB01ZD.f"
		i__1 = *n - 2;
#line 363 "TB01ZD.f"
		i__2 = *n - 2;
#line 363 "TB01ZD.f"
		dlacpy_("Lower", &i__1, &i__2, &a[a_dim1 + 3], lda, &z__[(
			z_dim1 << 1) + 3], ldz, (ftnlen)5);
#line 363 "TB01ZD.f"
	    }
#line 366 "TB01ZD.f"
	    if (ljobi) {

/*              Form the orthogonal transformation matrix Z = Z1 * Z2. */
/*              Workspace: need N;  prefer N*NB. */

#line 371 "TB01ZD.f"
		dorgqr_(n, n, n, &z__[z_offset], ldz, &tau[1], &dwork[1], 
			ldwork, info);
#line 372 "TB01ZD.f"
		wrkopt = max(wrkopt,dwork[1]);
#line 373 "TB01ZD.f"
	    }
#line 374 "TB01ZD.f"
	}

/*        Annihilate the lower part of A and B. */

#line 378 "TB01ZD.f"
	if (*n > 2) {
#line 378 "TB01ZD.f"
	    i__1 = *n - 2;
#line 378 "TB01ZD.f"
	    i__2 = *n - 2;
#line 378 "TB01ZD.f"
	    dlaset_("Lower", &i__1, &i__2, &c_b10, &c_b10, &a[a_dim1 + 3], 
		    lda, (ftnlen)5);
#line 378 "TB01ZD.f"
	}
#line 380 "TB01ZD.f"
	if (*n > 1) {
#line 380 "TB01ZD.f"
	    i__1 = *n - 1;
#line 380 "TB01ZD.f"
	    i__2 = *n - 1;
#line 380 "TB01ZD.f"
	    dlaset_("Full", &i__1, &c__1, &c_b10, &c_b10, &b[2], &i__2, (
		    ftnlen)4);
#line 380 "TB01ZD.f"
	}

/*        Find NCONT by checking sizes of the sub-diagonal elements of */
/*        transformed A. */

#line 386 "TB01ZD.f"
	if (*tol <= 0.) {
/* Computing MAX */
#line 386 "TB01ZD.f"
	    d__1 = fanorm, d__2 = abs(b1);
#line 386 "TB01ZD.f"
	    toldef = thresh * max(d__1,d__2);
#line 386 "TB01ZD.f"
	}

#line 389 "TB01ZD.f"
	j = 1;

/*        WHILE ( J < N and ABS( A(J+1,J) ) > TOLDEF ) DO */

#line 393 "TB01ZD.f"
L10:
#line 394 "TB01ZD.f"
	if (j < *n) {
#line 395 "TB01ZD.f"
	    if ((d__1 = a[j + 1 + j * a_dim1], abs(d__1)) > toldef) {
#line 396 "TB01ZD.f"
		++j;
#line 397 "TB01ZD.f"
		goto L10;
#line 398 "TB01ZD.f"
	    }
#line 399 "TB01ZD.f"
	}

/*        END WHILE 10 */

/*        First negligible sub-diagonal element found, if any: set NCONT. */

#line 405 "TB01ZD.f"
	*ncont = j;
#line 406 "TB01ZD.f"
	if (j < *n) {
#line 406 "TB01ZD.f"
	    a[j + 1 + j * a_dim1] = 0.;
#line 406 "TB01ZD.f"
	}

/*        Undo scaling of A and B. */

#line 411 "TB01ZD.f"
	mb01pd_("U", "H", ncont, ncont, &c__0, &c__0, &anorm, &c__0, nblk, &a[
		a_offset], lda, info, (ftnlen)1, (ftnlen)1);
#line 413 "TB01ZD.f"
	mb01pd_("U", "G", &c__1, &c__1, &c__0, &c__0, &bnorm, &c__0, nblk, &b[
		1], n, info, (ftnlen)1, (ftnlen)1);
#line 414 "TB01ZD.f"
	if (*ncont < *n) {
#line 414 "TB01ZD.f"
	    i__1 = *n - *ncont;
#line 414 "TB01ZD.f"
	    mb01pd_("U", "G", n, &i__1, &c__0, &c__0, &anorm, &c__0, nblk, &a[
		    (*ncont + 1) * a_dim1 + 1], lda, info, (ftnlen)1, (ftnlen)
		    1);
#line 414 "TB01ZD.f"
	}
#line 417 "TB01ZD.f"
    } else {

/*        B is negligible compared with A. No computations for reducing */
/*        the system to orthogonal canonical form have been performed, */
/*        except scaling (which is undoed). */

#line 423 "TB01ZD.f"
	mb01pd_("U", "G", n, n, &c__0, &c__0, &anorm, &c__0, nblk, &a[
		a_offset], lda, info, (ftnlen)1, (ftnlen)1);
#line 425 "TB01ZD.f"
	mb01pd_("U", "G", n, &c__1, &c__0, &c__0, &bnorm, &c__0, nblk, &b[1], 
		n, info, (ftnlen)1, (ftnlen)1);
#line 426 "TB01ZD.f"
	if (ljobf) {
#line 427 "TB01ZD.f"
	    dlaset_("Full", n, n, &c_b10, &c_b10, &z__[z_offset], ldz, (
		    ftnlen)4);
#line 428 "TB01ZD.f"
	    dlaset_("Full", n, &c__1, &c_b10, &c_b10, &tau[1], n, (ftnlen)4);
#line 429 "TB01ZD.f"
	} else if (ljobi) {
#line 430 "TB01ZD.f"
	    dlaset_("Full", n, n, &c_b10, &c_b18, &z__[z_offset], ldz, (
		    ftnlen)4);
#line 431 "TB01ZD.f"
	}
#line 432 "TB01ZD.f"
    }

/*     Set optimal workspace dimension. */

#line 436 "TB01ZD.f"
    dwork[1] = wrkopt;

#line 438 "TB01ZD.f"
    return 0;
/* *** Last line of TB01ZD *** */
} /* tb01zd_ */

