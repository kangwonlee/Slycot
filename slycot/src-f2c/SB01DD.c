#line 1 "SB01DD.f"
/* SB01DD.f -- translated by f2c (version 20100827).
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

#line 1 "SB01DD.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b24 = -1.;
static doublereal c_b26 = 1.;
static integer c__2 = 2;
static doublereal c_b99 = 0.;

/* Subroutine */ int sb01dd_(integer *n, integer *m, integer *indcon, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	nblk, doublereal *wr, doublereal *wi, doublereal *z__, integer *ldz, 
	doublereal *y, integer *count, doublereal *g, integer *ldg, 
	doublereal *tol, integer *iwork, doublereal *dwork, integer *ldwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, g_dim1, g_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, k, l;
    static doublereal p, q, r__, s;
    static integer m1, ia, nc, kk, mi, ni, ip, nj, mr, nr, lp1, mp1, np1, mr1,
	     nr1, kmr, rank;
    static doublereal sval[3];
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer iwrk, irmx;
    extern /* Subroutine */ int mb02qd_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen, ftnlen),
	     dscal_(integer *, doublereal *, doublereal *, integer *), dlarf_(
	    char *, integer *, integer *, doublereal *, integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, ftnlen), dgemm_(char *, 
	    char *, integer *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen), dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    extern doublereal dasum_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer indcn1, indcn2;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen), dlange_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    static integer nblkcr;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), dlaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen), xerbla_(char *, integer *, ftnlen);
    static integer indcrt;
    static logical complx;
    static integer maxwrk;
    static doublereal svlmax;


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

/*     To compute for a controllable matrix pair ( A, B ) a matrix G */
/*     such that the matrix A - B*G has the desired eigenstructure, */
/*     specified by desired eigenvalues and free eigenvector elements. */

/*     The pair ( A, B ) should be given in orthogonal canonical form */
/*     as returned by the SLICOT Library routine AB01ND. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A and the number of rows of the */
/*             matrix B.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrix B.  M >= 0. */

/*     INDCON  (input) INTEGER */
/*             The controllability index of the pair ( A, B ). */
/*             0 <= INDCON <= N. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the N-by-N matrix A in orthogonal canonical form, */
/*             as returned by SLICOT Library routine AB01ND. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the real Schur form of the matrix A - B*G. */
/*             The elements below the real Schur form of A are set to */
/*             zero. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the N-by-M matrix B in orthogonal canonical form, */
/*             as returned by SLICOT Library routine AB01ND. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the transformed matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     NBLK    (input) INTEGER array, dimension (N) */
/*             The leading INDCON elements of this array must contain the */
/*             orders of the diagonal blocks in the orthogonal canonical */
/*             form of A, as returned by SLICOT Library routine AB01ND. */
/*             The values of these elements must satisfy the following */
/*             conditions: */
/*             NBLK(1) >= NBLK(2) >= ... >= NBLK(INDCON), */
/*             NBLK(1) + NBLK(2) + ... + NBLK(INDCON) = N. */

/*     WR      (input) DOUBLE PRECISION array, dimension (N) */
/*     WI      (input) DOUBLE PRECISION array, dimension (N) */
/*             These arrays must contain the real and imaginary parts, */
/*             respectively, of the desired poles of the closed-loop */
/*             system, i.e., the eigenvalues of A - B*G. The poles can be */
/*             unordered, except that complex conjugate pairs of poles */
/*             must appear consecutively. */
/*             The elements of WI for complex eigenvalues are modified */
/*             internally, but restored on exit. */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the orthogonal matrix Z generated by SLICOT */
/*             Library routine AB01ND in the reduction of ( A, B ) to */
/*             orthogonal canonical form. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the orthogonal transformation matrix which reduces A - B*G */
/*             to real Schur form. */

/*     LDZ     INTEGER */
/*             The leading dimension of the array Z.  LDZ >= max(1,N). */

/*     Y       (input) DOUBLE PRECISION array, dimension (M*N) */
/*             Y contains elements which are used as free parameters */
/*             in the eigenstructure design. The values of these */
/*             parameters are often set by an external optimization */
/*             procedure. */

/*     COUNT   (output) INTEGER */
/*             The actual number of elements in Y used as free */
/*             eigenvector and feedback matrix elements in the */
/*             eigenstructure design. */

/*     G       (output) DOUBLE PRECISION array, dimension (LDG,N) */
/*             The leading M-by-N part of this array contains the */
/*             feedback matrix which assigns the desired eigenstructure */
/*             of A - B*G. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G.  LDG >= max(1,M). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in rank determination when */
/*             transforming (A, B). If the user sets TOL > 0, then */
/*             the given value of TOL is used as a lower bound for the */
/*             reciprocal condition number (see the description of the */
/*             argument RCOND in the SLICOT routine MB03OD);  a */
/*             (sub)matrix whose estimated condition number is less than */
/*             1/TOL is considered to be of full rank.  If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance, */
/*             defined by  TOLDEF = N*N*EPS,  is used instead, where */
/*             EPS  is the machine precision (see LAPACK Library routine */
/*             DLAMCH). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(M*N,M*M+2*N+4*M+1). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the pair ( A, B ) is not controllable or the free */
/*                   parameters are not set appropriately. */

/*     METHOD */

/*     The routine implements the method proposed in [1], [2]. */

/*     REFERENCES */

/*     [1] Petkov, P.Hr., Konstantinov, M.M., Gu, D.W. and */
/*         Postlethwaite, I. */
/*         Optimal pole assignment design of linear multi-input systems. */
/*         Report 96-11, Department of Engineering, Leicester University, */
/*         1996. */

/*     [2] Petkov, P.Hr., Christov, N.D. and Konstantinov, M.M. */
/*         A computational algorithm for pole assignment of linear multi */
/*         input systems. IEEE Trans. Automatic Control, vol. AC-31, */
/*         pp. 1044-1047, 1986. */

/*     NUMERICAL ASPECTS */

/*     The method implemented is backward stable. */

/*     FURTHER COMMENTS */

/*     The eigenvalues of the real Schur form matrix As, returned in the */
/*     array A, are very close to the desired eigenvalues WR+WI*i. */
/*     However, the eigenvalues of the closed-loop matrix A - B*G, */
/*     computed by the QR algorithm using the matrices A and B, given on */
/*     entry, may be far from WR+WI*i, although the relative error */
/*        norm( Z'*(A - B*G)*Z - As )/norm( As ) */
/*     is close to machine accuracy. This may happen when the eigenvalue */
/*     problem for the matrix A - B*G is ill-conditioned. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, Technical University of Sofia, Oct. 1998. */
/*     V. Sima, Katholieke Universiteit Leuven, Jan. 1999, SLICOT Library */
/*     version. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2005. */

/*     KEYWORDS */

/*     Closed loop spectrum, closed loop systems, eigenvalue assignment, */
/*     orthogonal canonical form, orthogonal transformation, pole */
/*     placement, Schur form. */

/*     ****************************************************************** */

/*     .. Parameters .. */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments. */

#line 247 "SB01DD.f"
    /* Parameter adjustments */
#line 247 "SB01DD.f"
    a_dim1 = *lda;
#line 247 "SB01DD.f"
    a_offset = 1 + a_dim1;
#line 247 "SB01DD.f"
    a -= a_offset;
#line 247 "SB01DD.f"
    b_dim1 = *ldb;
#line 247 "SB01DD.f"
    b_offset = 1 + b_dim1;
#line 247 "SB01DD.f"
    b -= b_offset;
#line 247 "SB01DD.f"
    --nblk;
#line 247 "SB01DD.f"
    --wr;
#line 247 "SB01DD.f"
    --wi;
#line 247 "SB01DD.f"
    z_dim1 = *ldz;
#line 247 "SB01DD.f"
    z_offset = 1 + z_dim1;
#line 247 "SB01DD.f"
    z__ -= z_offset;
#line 247 "SB01DD.f"
    --y;
#line 247 "SB01DD.f"
    g_dim1 = *ldg;
#line 247 "SB01DD.f"
    g_offset = 1 + g_dim1;
#line 247 "SB01DD.f"
    g -= g_offset;
#line 247 "SB01DD.f"
    --iwork;
#line 247 "SB01DD.f"
    --dwork;
#line 247 "SB01DD.f"

#line 247 "SB01DD.f"
    /* Function Body */
#line 247 "SB01DD.f"
    *info = 0;
#line 248 "SB01DD.f"
    nr = 0;
/* Computing MAX */
#line 249 "SB01DD.f"
    i__1 = *m * *n, i__2 = *m * *m + (*n << 1) + (*m << 2) + 1;
#line 249 "SB01DD.f"
    iwrk = max(i__1,i__2);
#line 250 "SB01DD.f"
    i__1 = min(*indcon,*n);
#line 250 "SB01DD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 251 "SB01DD.f"
	nr += nblk[i__];
#line 252 "SB01DD.f"
	if (i__ > 1) {
#line 253 "SB01DD.f"
	    if (nblk[i__ - 1] < nblk[i__]) {
#line 253 "SB01DD.f"
		*info = -8;
#line 253 "SB01DD.f"
	    }
#line 255 "SB01DD.f"
	}
#line 256 "SB01DD.f"
/* L10: */
#line 256 "SB01DD.f"
    }
#line 257 "SB01DD.f"
    if (*n < 0) {
#line 258 "SB01DD.f"
	*info = -1;
#line 259 "SB01DD.f"
    } else if (*m < 0) {
#line 260 "SB01DD.f"
	*info = -2;
#line 261 "SB01DD.f"
    } else if (*indcon < 0 || *indcon > *n) {
#line 262 "SB01DD.f"
	*info = -3;
#line 263 "SB01DD.f"
    } else if (*lda < max(1,*n)) {
#line 264 "SB01DD.f"
	*info = -5;
#line 265 "SB01DD.f"
    } else if (*ldb < max(1,*n)) {
#line 266 "SB01DD.f"
	*info = -7;
#line 267 "SB01DD.f"
    } else if (nr != *n) {
#line 268 "SB01DD.f"
	*info = -8;
#line 269 "SB01DD.f"
    } else if (*ldz < max(1,*n)) {
#line 270 "SB01DD.f"
	*info = -12;
#line 271 "SB01DD.f"
    } else if (*ldg < max(1,*m)) {
#line 272 "SB01DD.f"
	*info = -16;
#line 273 "SB01DD.f"
    } else if (*ldwork < iwrk) {
#line 274 "SB01DD.f"
	*info = -20;
#line 275 "SB01DD.f"
    }

#line 277 "SB01DD.f"
    if (*info != 0) {
#line 278 "SB01DD.f"
	i__1 = -(*info);
#line 278 "SB01DD.f"
	xerbla_("SB01DD", &i__1, (ftnlen)6);
#line 279 "SB01DD.f"
	return 0;
#line 280 "SB01DD.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 284 "SB01DD.f"
    i__1 = min(*m,*n);
#line 284 "SB01DD.f"
    if (min(i__1,*indcon) == 0) {
#line 285 "SB01DD.f"
	*count = 0;
#line 286 "SB01DD.f"
	dwork[1] = 1.;
#line 287 "SB01DD.f"
	return 0;
#line 288 "SB01DD.f"
    }

#line 290 "SB01DD.f"
    maxwrk = iwrk;
#line 291 "SB01DD.f"
    toldef = *tol;
#line 292 "SB01DD.f"
    if (toldef <= 0.) {

/*        Use the default tolerance, based on machine precision. */

#line 296 "SB01DD.f"
	toldef = (doublereal) (*n * *n) * dlamch_("EPSILON", (ftnlen)7);
#line 297 "SB01DD.f"
    }

#line 299 "SB01DD.f"
    irmx = (*n << 1) + 1;
#line 300 "SB01DD.f"
    iwrk = irmx + *m * *m;
#line 301 "SB01DD.f"
    m1 = nblk[1];
#line 302 "SB01DD.f"
    *count = 1;
#line 303 "SB01DD.f"
    indcrt = *indcon;
#line 304 "SB01DD.f"
    nblkcr = nblk[indcrt];

/*     Compute the Frobenius norm of [ B  A ] (used for rank estimation), */
/*     taking into account the structure. */

#line 309 "SB01DD.f"
    nr = m1;
#line 310 "SB01DD.f"
    nc = 1;
#line 311 "SB01DD.f"
    svlmax = dlange_("Frobenius", &m1, m, &b[b_offset], ldb, &dwork[1], (
	    ftnlen)9);

#line 313 "SB01DD.f"
    i__1 = indcrt - 1;
#line 313 "SB01DD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 314 "SB01DD.f"
	nr += nblk[i__ + 1];
#line 315 "SB01DD.f"
	d__1 = dlange_("Frobenius", &nr, &nblk[i__], &a[nc * a_dim1 + 1], lda,
		 &dwork[1], (ftnlen)9);
#line 315 "SB01DD.f"
	svlmax = dlapy2_(&svlmax, &d__1);
#line 318 "SB01DD.f"
	nc += nblk[i__];
#line 319 "SB01DD.f"
/* L20: */
#line 319 "SB01DD.f"
    }

#line 321 "SB01DD.f"
    d__1 = dlange_("Frobenius", n, &nblkcr, &a[nc * a_dim1 + 1], lda, &dwork[
	    1], (ftnlen)9);
#line 321 "SB01DD.f"
    svlmax = dlapy2_(&svlmax, &d__1);
#line 324 "SB01DD.f"
    l = 1;
#line 325 "SB01DD.f"
    mr = nblkcr;
#line 326 "SB01DD.f"
    nr = *n - mr + 1;
#line 327 "SB01DD.f"
L30:
/*     WHILE( INDCRT.GT.1 )LOOP */
#line 329 "SB01DD.f"
    if (indcrt > 1) {

/*        Assign next eigenvalue/eigenvector. */

#line 333 "SB01DD.f"
	lp1 = l + m1;
#line 334 "SB01DD.f"
	indcn1 = indcrt - 1;
#line 335 "SB01DD.f"
	mr1 = nblk[indcn1];
#line 336 "SB01DD.f"
	nr1 = nr - mr1;
#line 337 "SB01DD.f"
	complx = wi[l] != 0.;
#line 338 "SB01DD.f"
	dcopy_(&mr, &y[*count], &c__1, &dwork[nr], &c__1);
#line 339 "SB01DD.f"
	*count += mr;
#line 340 "SB01DD.f"
	nc = 1;
#line 341 "SB01DD.f"
	if (complx) {
#line 342 "SB01DD.f"
	    dcopy_(&mr, &y[*count], &c__1, &dwork[*n + nr], &c__1);
#line 343 "SB01DD.f"
	    *count += mr;
#line 344 "SB01DD.f"
	    wi[l + 1] = wi[l] * wi[l + 1];
#line 345 "SB01DD.f"
	    nc = 2;
#line 346 "SB01DD.f"
	}

/*        Compute and transform eiegenvector. */

#line 350 "SB01DD.f"
	i__1 = indcrt;
#line 350 "SB01DD.f"
	for (ip = 1; ip <= i__1; ++ip) {
#line 351 "SB01DD.f"
	    if (ip != indcrt) {
#line 352 "SB01DD.f"
		dlacpy_("Full", &mr, &mr1, &a[nr + nr1 * a_dim1], lda, &dwork[
			irmx], m, (ftnlen)4);
#line 354 "SB01DD.f"
		if (ip == 1) {
#line 355 "SB01DD.f"
		    mp1 = mr;
#line 356 "SB01DD.f"
		    np1 = nr + mp1;
#line 357 "SB01DD.f"
		} else {
#line 358 "SB01DD.f"
		    mp1 = mr + 1;
#line 359 "SB01DD.f"
		    np1 = nr + mp1;
#line 360 "SB01DD.f"
		    s = dasum_(&mp1, &dwork[nr], &c__1);
#line 361 "SB01DD.f"
		    if (complx) {
#line 361 "SB01DD.f"
			s += dasum_(&mp1, &dwork[*n + nr], &c__1);
#line 361 "SB01DD.f"
		    }
#line 362 "SB01DD.f"
		    if (s != 0.) {

/*                    Scale eigenvector elements. */

#line 366 "SB01DD.f"
			d__1 = 1. / s;
#line 366 "SB01DD.f"
			dscal_(&mp1, &d__1, &dwork[nr], &c__1);
#line 367 "SB01DD.f"
			if (complx) {
#line 368 "SB01DD.f"
			    d__1 = 1. / s;
#line 368 "SB01DD.f"
			    dscal_(&mp1, &d__1, &dwork[*n + nr], &c__1);
#line 369 "SB01DD.f"
			    if (np1 <= *n) {
#line 369 "SB01DD.f"
				dwork[*n + np1] /= s;
#line 369 "SB01DD.f"
			    }
#line 371 "SB01DD.f"
			}
#line 372 "SB01DD.f"
		    }
#line 373 "SB01DD.f"
		}

/*              Compute the right-hand side of the eigenvector equations. */

#line 377 "SB01DD.f"
		dcopy_(&mr, &dwork[nr], &c__1, &dwork[nr1], &c__1);
#line 378 "SB01DD.f"
		dscal_(&mr, &wr[l], &dwork[nr1], &c__1);
#line 379 "SB01DD.f"
		dgemv_("No transpose", &mr, &mp1, &c_b24, &a[nr + nr * a_dim1]
			, lda, &dwork[nr], &c__1, &c_b26, &dwork[nr1], &c__1, 
			(ftnlen)12);
#line 381 "SB01DD.f"
		if (complx) {
#line 382 "SB01DD.f"
		    daxpy_(&mr, &wi[l + 1], &dwork[*n + nr], &c__1, &dwork[
			    nr1], &c__1);
#line 384 "SB01DD.f"
		    dcopy_(&mr, &dwork[nr], &c__1, &dwork[*n + nr1], &c__1);
#line 385 "SB01DD.f"
		    daxpy_(&mr, &wr[l + 1], &dwork[*n + nr], &c__1, &dwork[*n 
			    + nr1], &c__1);
#line 387 "SB01DD.f"
		    dgemv_("No transpose", &mr, &mp1, &c_b24, &a[nr + nr * 
			    a_dim1], lda, &dwork[*n + nr], &c__1, &c_b26, &
			    dwork[*n + nr1], &c__1, (ftnlen)12);
#line 390 "SB01DD.f"
		    if (np1 <= *n) {
#line 390 "SB01DD.f"
			d__1 = -dwork[*n + np1];
#line 390 "SB01DD.f"
			daxpy_(&mr, &d__1, &a[nr + np1 * a_dim1], &c__1, &
				dwork[*n + nr1], &c__1);
#line 390 "SB01DD.f"
		    }
#line 393 "SB01DD.f"
		}

/*              Solve linear equations for eigenvector elements. */

#line 397 "SB01DD.f"
		i__2 = *ldwork - iwrk + 1;
#line 397 "SB01DD.f"
		mb02qd_("FreeElements", "NoPermuting", &mr, &mr1, &nc, &
			toldef, &svlmax, &dwork[irmx], m, &dwork[nr1], n, &y[*
			count], &iwork[1], &rank, sval, &dwork[iwrk], &i__2, 
			info, (ftnlen)12, (ftnlen)11);
/* Computing MAX */
#line 401 "SB01DD.f"
		i__2 = maxwrk, i__3 = (integer) dwork[iwrk] + iwrk - 1;
#line 401 "SB01DD.f"
		maxwrk = max(i__2,i__3);
#line 402 "SB01DD.f"
		if (rank < mr) {
#line 402 "SB01DD.f"
		    goto L80;
#line 402 "SB01DD.f"
		}

#line 404 "SB01DD.f"
		*count += (mr1 - mr) * nc;
#line 405 "SB01DD.f"
		nj = nr1;
#line 406 "SB01DD.f"
	    } else {
#line 407 "SB01DD.f"
		nj = nr;
#line 408 "SB01DD.f"
	    }
#line 409 "SB01DD.f"
	    ni = nr + mr - 1;
#line 410 "SB01DD.f"
	    if (ip == 1) {
#line 411 "SB01DD.f"
		kmr = mr - 1;
#line 412 "SB01DD.f"
	    } else {
#line 413 "SB01DD.f"
		kmr = mr;
#line 414 "SB01DD.f"
		if (ip == 2) {
#line 415 "SB01DD.f"
		    ni += nblkcr;
#line 416 "SB01DD.f"
		} else {
#line 417 "SB01DD.f"
		    ni = ni + nblk[indcrt - ip + 2] + 1;
#line 418 "SB01DD.f"
		    if (complx) {
/* Computing MIN */
#line 418 "SB01DD.f"
			i__2 = ni + 1;
#line 418 "SB01DD.f"
			ni = min(i__2,*n);
#line 418 "SB01DD.f"
		    }
#line 419 "SB01DD.f"
		}
#line 420 "SB01DD.f"
	    }

#line 422 "SB01DD.f"
	    i__2 = kmr;
#line 422 "SB01DD.f"
	    for (kk = 1; kk <= i__2; ++kk) {
#line 423 "SB01DD.f"
		k = nr + mr - kk;
#line 424 "SB01DD.f"
		if (ip == 1) {
#line 424 "SB01DD.f"
		    k = *n - kk;
#line 424 "SB01DD.f"
		}
#line 425 "SB01DD.f"
		dlartg_(&dwork[k], &dwork[k + 1], &p, &q, &r__);
#line 426 "SB01DD.f"
		dwork[k] = r__;
#line 427 "SB01DD.f"
		dwork[k + 1] = 0.;

/*              Transform  A. */

#line 431 "SB01DD.f"
		i__3 = *n - nj + 1;
#line 431 "SB01DD.f"
		drot_(&i__3, &a[k + nj * a_dim1], lda, &a[k + 1 + nj * a_dim1]
			, lda, &p, &q);
#line 433 "SB01DD.f"
		drot_(&ni, &a[k * a_dim1 + 1], &c__1, &a[(k + 1) * a_dim1 + 1]
			, &c__1, &p, &q);

#line 435 "SB01DD.f"
		if (k < lp1) {

/*                 Transform B. */

#line 439 "SB01DD.f"
		    drot_(m, &b[k + b_dim1], ldb, &b[k + 1 + b_dim1], ldb, &p,
			     &q);
#line 440 "SB01DD.f"
		}

/*              Accumulate transformations. */

#line 444 "SB01DD.f"
		drot_(n, &z__[k * z_dim1 + 1], &c__1, &z__[(k + 1) * z_dim1 + 
			1], &c__1, &p, &q);

#line 446 "SB01DD.f"
		if (complx) {
#line 447 "SB01DD.f"
		    drot_(&c__1, &dwork[*n + k], &c__1, &dwork[*n + k + 1], &
			    c__1, &p, &q);
#line 449 "SB01DD.f"
		    ++k;
#line 450 "SB01DD.f"
		    if (k < *n) {
#line 451 "SB01DD.f"
			dlartg_(&dwork[*n + k], &dwork[*n + k + 1], &p, &q, &
				r__);
#line 453 "SB01DD.f"
			dwork[*n + k] = r__;
#line 454 "SB01DD.f"
			dwork[*n + k + 1] = 0.;

/*                    Transform  A. */

#line 458 "SB01DD.f"
			i__3 = *n - nj + 1;
#line 458 "SB01DD.f"
			drot_(&i__3, &a[k + nj * a_dim1], lda, &a[k + 1 + nj *
				 a_dim1], lda, &p, &q);
#line 460 "SB01DD.f"
			drot_(&ni, &a[k * a_dim1 + 1], &c__1, &a[(k + 1) * 
				a_dim1 + 1], &c__1, &p, &q);

#line 462 "SB01DD.f"
			if (k <= lp1) {

/*                       Transform B. */

#line 466 "SB01DD.f"
			    drot_(m, &b[k + b_dim1], ldb, &b[k + 1 + b_dim1], 
				    ldb, &p, &q);
#line 468 "SB01DD.f"
			}

/*                    Accumulate transformations. */

#line 472 "SB01DD.f"
			drot_(n, &z__[k * z_dim1 + 1], &c__1, &z__[(k + 1) * 
				z_dim1 + 1], &c__1, &p, &q);

#line 474 "SB01DD.f"
		    }
#line 475 "SB01DD.f"
		}
#line 476 "SB01DD.f"
/* L40: */
#line 476 "SB01DD.f"
	    }

#line 478 "SB01DD.f"
	    if (ip != indcrt) {
#line 479 "SB01DD.f"
		mr = mr1;
#line 480 "SB01DD.f"
		nr = nr1;
#line 481 "SB01DD.f"
		if (ip != indcn1) {
#line 482 "SB01DD.f"
		    indcn2 = indcrt - ip - 1;
#line 483 "SB01DD.f"
		    mr1 = nblk[indcn2];
#line 484 "SB01DD.f"
		    nr1 -= mr1;
#line 485 "SB01DD.f"
		}
#line 486 "SB01DD.f"
	    }
#line 487 "SB01DD.f"
/* L50: */
#line 487 "SB01DD.f"
	}

#line 489 "SB01DD.f"
	if (! complx) {

/*           Find one column of G. */

#line 493 "SB01DD.f"
	    dlacpy_("Full", &m1, m, &b[l + 1 + b_dim1], ldb, &dwork[irmx], m, 
		    (ftnlen)4);
#line 495 "SB01DD.f"
	    dcopy_(&m1, &a[l + 1 + l * a_dim1], &c__1, &g[l * g_dim1 + 1], &
		    c__1);
#line 496 "SB01DD.f"
	} else {

/*           Find two columns of G. */

#line 500 "SB01DD.f"
	    if (lp1 < *n) {
#line 501 "SB01DD.f"
		++lp1;
#line 502 "SB01DD.f"
		k = l + 2;
#line 503 "SB01DD.f"
	    } else {
#line 504 "SB01DD.f"
		k = l + 1;
#line 505 "SB01DD.f"
	    }
#line 506 "SB01DD.f"
	    dlacpy_("Full", &m1, m, &b[k + b_dim1], ldb, &dwork[irmx], m, (
		    ftnlen)4);
#line 508 "SB01DD.f"
	    dlacpy_("Full", &m1, &c__2, &a[k + l * a_dim1], lda, &g[l * 
		    g_dim1 + 1], ldg, (ftnlen)4);
#line 509 "SB01DD.f"
	    if (k == l + 1) {
#line 510 "SB01DD.f"
		g[l * g_dim1 + 1] -= dwork[*n + l + 1] / dwork[l] * wi[l + 1];
#line 512 "SB01DD.f"
		g[(l + 1) * g_dim1 + 1] = g[(l + 1) * g_dim1 + 1] - wr[l + 1] 
			+ dwork[*n + l] / dwork[l] * wi[l + 1];
#line 514 "SB01DD.f"
	    }
#line 515 "SB01DD.f"
	}

#line 517 "SB01DD.f"
	i__1 = *ldwork - iwrk + 1;
#line 517 "SB01DD.f"
	mb02qd_("FreeElements", "NoPermuting", &m1, m, &nc, &toldef, &svlmax, 
		&dwork[irmx], m, &g[l * g_dim1 + 1], ldg, &y[*count], &iwork[
		1], &rank, sval, &dwork[iwrk], &i__1, info, (ftnlen)12, (
		ftnlen)11);
/* Computing MAX */
#line 521 "SB01DD.f"
	i__1 = maxwrk, i__2 = (integer) dwork[iwrk] + iwrk - 1;
#line 521 "SB01DD.f"
	maxwrk = max(i__1,i__2);
#line 522 "SB01DD.f"
	if (rank < m1) {
#line 522 "SB01DD.f"
	    goto L80;
#line 522 "SB01DD.f"
	}

#line 524 "SB01DD.f"
	*count += (*m - m1) * nc;
#line 525 "SB01DD.f"
	dgemm_("No transpose", "No transpose", &lp1, &nc, m, &c_b24, &b[
		b_offset], ldb, &g[l * g_dim1 + 1], ldg, &c_b26, &a[l * 
		a_dim1 + 1], lda, (ftnlen)12, (ftnlen)12);
#line 527 "SB01DD.f"
	++l;
#line 528 "SB01DD.f"
	--nblkcr;
#line 529 "SB01DD.f"
	if (nblkcr == 0) {
#line 530 "SB01DD.f"
	    --indcrt;
#line 531 "SB01DD.f"
	    nblkcr = nblk[indcrt];
#line 532 "SB01DD.f"
	}
#line 533 "SB01DD.f"
	if (complx) {
#line 534 "SB01DD.f"
	    wi[l] = -wi[l - 1];
#line 535 "SB01DD.f"
	    ++l;
#line 536 "SB01DD.f"
	    --nblkcr;
#line 537 "SB01DD.f"
	    if (nblkcr == 0) {
#line 538 "SB01DD.f"
		--indcrt;
#line 539 "SB01DD.f"
		if (indcrt > 0) {
#line 539 "SB01DD.f"
		    nblkcr = nblk[indcrt];
#line 539 "SB01DD.f"
		}
#line 540 "SB01DD.f"
	    }
#line 541 "SB01DD.f"
	}
#line 542 "SB01DD.f"
	mr = nblkcr;
#line 543 "SB01DD.f"
	nr = *n - mr + 1;
#line 544 "SB01DD.f"
	goto L30;
#line 545 "SB01DD.f"
    }
/*     END WHILE 30 */

#line 548 "SB01DD.f"
    if (l <= *n) {

/*        Find the remaining columns of G. */

/*        QR decomposition of the free eigenvectors. */

#line 554 "SB01DD.f"
	i__1 = mr - 1;
#line 554 "SB01DD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 555 "SB01DD.f"
	    ia = l + i__ - 1;
#line 556 "SB01DD.f"
	    mi = mr - i__ + 1;
#line 557 "SB01DD.f"
	    dcopy_(&mi, &y[*count], &c__1, &dwork[1], &c__1);
#line 558 "SB01DD.f"
	    *count += mi;
#line 559 "SB01DD.f"
	    dlarfg_(&mi, &dwork[1], &dwork[2], &c__1, &r__);
#line 560 "SB01DD.f"
	    dwork[1] = 1.;

/*           Transform A. */

#line 564 "SB01DD.f"
	    dlarf_("Left", &mi, &mr, &dwork[1], &c__1, &r__, &a[ia + l * 
		    a_dim1], lda, &dwork[*n + 1], (ftnlen)4);
#line 566 "SB01DD.f"
	    dlarf_("Right", n, &mi, &dwork[1], &c__1, &r__, &a[ia * a_dim1 + 
		    1], lda, &dwork[*n + 1], (ftnlen)5);

/*           Transform B. */

#line 571 "SB01DD.f"
	    dlarf_("Left", &mi, m, &dwork[1], &c__1, &r__, &b[ia + b_dim1], 
		    ldb, &dwork[*n + 1], (ftnlen)4);

/*           Accumulate transformations. */

#line 576 "SB01DD.f"
	    dlarf_("Right", n, &mi, &dwork[1], &c__1, &r__, &z__[ia * z_dim1 
		    + 1], ldz, &dwork[*n + 1], (ftnlen)5);
#line 578 "SB01DD.f"
/* L60: */
#line 578 "SB01DD.f"
	}

#line 580 "SB01DD.f"
	i__ = 0;
/*        REPEAT */
#line 582 "SB01DD.f"
L70:
#line 583 "SB01DD.f"
	++i__;
#line 584 "SB01DD.f"
	ia = l + i__ - 1;
#line 585 "SB01DD.f"
	if (wi[ia] == 0.) {
#line 586 "SB01DD.f"
	    dcopy_(&mr, &a[ia + l * a_dim1], lda, &g[i__ + l * g_dim1], ldg);
#line 587 "SB01DD.f"
	    i__1 = mr - i__;
#line 587 "SB01DD.f"
	    daxpy_(&i__1, &c_b24, &y[*count], &c__1, &g[i__ + (l + i__) * 
		    g_dim1], ldg);
#line 588 "SB01DD.f"
	    *count = *count + mr - i__;
#line 589 "SB01DD.f"
	    g[i__ + ia * g_dim1] -= wr[ia];
#line 590 "SB01DD.f"
	} else {
#line 591 "SB01DD.f"
	    dlacpy_("Full", &c__2, &mr, &a[ia + l * a_dim1], lda, &g[i__ + l *
		     g_dim1], ldg, (ftnlen)4);
#line 593 "SB01DD.f"
	    i__1 = mr - i__ - 1;
#line 593 "SB01DD.f"
	    daxpy_(&i__1, &c_b24, &y[*count], &c__2, &g[i__ + (l + i__ + 1) * 
		    g_dim1], ldg);
#line 595 "SB01DD.f"
	    i__1 = mr - i__ - 1;
#line 595 "SB01DD.f"
	    daxpy_(&i__1, &c_b24, &y[*count + 1], &c__2, &g[i__ + 1 + (l + 
		    i__ + 1) * g_dim1], ldg);
#line 597 "SB01DD.f"
	    *count += mr - i__ - 1 << 1;
#line 598 "SB01DD.f"
	    g[i__ + ia * g_dim1] -= wr[ia];
#line 599 "SB01DD.f"
	    g[i__ + (ia + 1) * g_dim1] -= wi[ia];
#line 600 "SB01DD.f"
	    g[i__ + 1 + ia * g_dim1] -= wi[ia + 1];
#line 601 "SB01DD.f"
	    g[i__ + 1 + (ia + 1) * g_dim1] -= wr[ia + 1];
#line 602 "SB01DD.f"
	    ++i__;
#line 603 "SB01DD.f"
	}
#line 604 "SB01DD.f"
	if (i__ < mr) {
#line 604 "SB01DD.f"
	    goto L70;
#line 604 "SB01DD.f"
	}
/*        UNTIL I.GE.MR */

#line 607 "SB01DD.f"
	dlacpy_("Full", &mr, m, &b[l + b_dim1], ldb, &dwork[irmx], m, (ftnlen)
		4);
#line 608 "SB01DD.f"
	i__1 = *ldwork - iwrk + 1;
#line 608 "SB01DD.f"
	mb02qd_("FreeElements", "NoPermuting", &mr, m, &mr, &toldef, &svlmax, 
		&dwork[irmx], m, &g[l * g_dim1 + 1], ldg, &y[*count], &iwork[
		1], &rank, sval, &dwork[iwrk], &i__1, info, (ftnlen)12, (
		ftnlen)11);
/* Computing MAX */
#line 612 "SB01DD.f"
	i__1 = maxwrk, i__2 = (integer) dwork[iwrk] + iwrk - 1;
#line 612 "SB01DD.f"
	maxwrk = max(i__1,i__2);
#line 613 "SB01DD.f"
	if (rank < mr) {
#line 613 "SB01DD.f"
	    goto L80;
#line 613 "SB01DD.f"
	}

#line 615 "SB01DD.f"
	*count += (*m - mr) * mr;
#line 616 "SB01DD.f"
	dgemm_("No transpose", "No transpose", n, &mr, m, &c_b24, &b[b_offset]
		, ldb, &g[l * g_dim1 + 1], ldg, &c_b26, &a[l * a_dim1 + 1], 
		lda, (ftnlen)12, (ftnlen)12);
#line 618 "SB01DD.f"
    }

/*     Transform G: */
/*     G := G * Z'. */

#line 623 "SB01DD.f"
    dgemm_("No transpose", "Transpose", m, n, n, &c_b26, &g[g_offset], ldg, &
	    z__[z_offset], ldz, &c_b99, &dwork[1], m, (ftnlen)12, (ftnlen)9);
#line 625 "SB01DD.f"
    dlacpy_("Full", m, n, &dwork[1], m, &g[g_offset], ldg, (ftnlen)4);
#line 626 "SB01DD.f"
    --(*count);

#line 628 "SB01DD.f"
    if (*n > 2) {

/*        Set the elements of A below the Hessenberg part to zero. */

#line 632 "SB01DD.f"
	i__1 = *n - 2;
#line 632 "SB01DD.f"
	i__2 = *n - 2;
#line 632 "SB01DD.f"
	dlaset_("Lower", &i__1, &i__2, &c_b99, &c_b99, &a[a_dim1 + 3], lda, (
		ftnlen)5);
#line 633 "SB01DD.f"
    }
#line 634 "SB01DD.f"
    dwork[1] = (doublereal) maxwrk;
#line 635 "SB01DD.f"
    return 0;

/*     Exit with INFO = 1 if the pair ( A, B ) is not controllable or */
/*     the free parameters are not set appropriately. */

#line 640 "SB01DD.f"
L80:
#line 640 "SB01DD.f"
    *info = 1;
#line 641 "SB01DD.f"
    return 0;
/* *** Last line of SB01DD *** */
} /* sb01dd_ */

