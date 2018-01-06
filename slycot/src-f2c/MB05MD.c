#line 1 "MB05MD.f"
/* MB05MD.f -- translated by f2c (version 20100827).
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

#line 1 "MB05MD.f"
/* Table of constant values */

static doublereal c_b16 = 1.;
static integer c__1 = 1;
static integer c__2 = 2;
static doublereal c_b40 = 0.;

/* Subroutine */ int mb05md_(char *balanc, integer *n, doublereal *delta, 
	doublereal *a, integer *lda, doublereal *v, integer *ldv, doublereal *
	y, integer *ldy, doublereal *valr, doublereal *vali, integer *iwork, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen balanc_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, v_dim1, v_offset, y_dim1, y_offset, i__1, i__2;

    /* Builtin functions */
    double exp(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal tmp[4]	/* was [2][2] */;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical scale;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal rcond;
    extern /* Subroutine */ int mb05my_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);
    static doublereal tempi;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static doublereal tempr;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dgebak_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dtrcon_(char *, char *, char *
	    , integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
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

/*     To compute exp(A*delta) where A is a real N-by-N non-defective */
/*     matrix with real or complex eigenvalues and delta is a scalar */
/*     value. The routine also returns the eigenvalues and eigenvectors */
/*     of A as well as (if all eigenvalues are real) the matrix product */
/*     exp(Lambda*delta) times the inverse of the eigenvector matrix */
/*     of A, where Lambda is the diagonal matrix of eigenvalues. */
/*     Optionally, the routine computes a balancing transformation to */
/*     improve the conditioning of the eigenvalues and eigenvectors. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     BALANC  CHARACTER*1 */
/*             Indicates how the input matrix should be diagonally scaled */
/*             to improve the conditioning of its eigenvalues as follows: */
/*             = 'N':  Do not diagonally scale; */
/*             = 'S':  Diagonally scale the matrix, i.e. replace A by */
/*                     D*A*D**(-1), where D is a diagonal matrix chosen */
/*                     to make the rows and columns of A more equal in */
/*                     norm. Do not permute. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     DELTA   (input) DOUBLE PRECISION */
/*             The scalar value delta of the problem. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A of the problem. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the solution matrix exp(A*delta). */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= max(1,N). */

/*     V       (output) DOUBLE PRECISION array, dimension (LDV,N) */
/*             The leading N-by-N part of this array contains the */
/*             eigenvector matrix for A. */
/*             If the k-th eigenvalue is real the k-th column of the */
/*             eigenvector matrix holds the eigenvector corresponding */
/*             to the k-th eigenvalue. */
/*             Otherwise, the k-th and (k+1)-th eigenvalues form a */
/*             complex conjugate pair and the k-th and (k+1)-th columns */
/*             of the eigenvector matrix hold the real and imaginary */
/*             parts of the eigenvectors corresponding to these */
/*             eigenvalues as follows. */
/*             If p and q denote the k-th and (k+1)-th columns of the */
/*             eigenvector matrix, respectively, then the eigenvector */
/*             corresponding to the complex eigenvalue with positive */
/*             (negative) imaginary value is given by */
/*                                       2 */
/*             p + q*j (p - q*j), where j  = -1. */

/*     LDV     INTEGER */
/*             The leading dimension of array V.  LDV >= max(1,N). */

/*     Y       (output) DOUBLE PRECISION array, dimension (LDY,N) */
/*             The leading N-by-N part of this array contains an */
/*             intermediate result for computing the matrix exponential. */
/*             Specifically, exp(A*delta) is obtained as the product V*Y, */
/*             where V is the matrix stored in the leading N-by-N part of */
/*             the array V. If all eigenvalues of A are real, then the */
/*             leading N-by-N part of this array contains the matrix */
/*             product exp(Lambda*delta) times the inverse of the (right) */
/*             eigenvector matrix of A, where Lambda is the diagonal */
/*             matrix of eigenvalues. */

/*     LDY     INTEGER */
/*             The leading dimension of array Y.  LDY >= max(1,N). */

/*     VALR    (output) DOUBLE PRECISION array, dimension (N) */
/*     VALI    (output) DOUBLE PRECISION array, dimension (N) */
/*             These arrays contain the real and imaginary parts, */
/*             respectively, of the eigenvalues of the matrix A. The */
/*             eigenvalues are unordered except that complex conjugate */
/*             pairs of values appear consecutively with the eigenvalue */
/*             having positive imaginary part first. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK, and if N > 0, DWORK(2) returns the reciprocal */
/*             condition number of the triangular matrix used to obtain */
/*             the inverse of the eigenvector matrix. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= max(1,4*N). */
/*             For good performance, LDWORK must generally be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = i:  if INFO = i, the QR algorithm failed to compute all */
/*                   the eigenvalues; no eigenvectors have been computed; */
/*                   elements i+1:N of VALR and VALI contain eigenvalues */
/*                   which have converged; */
/*             = N+1:  if the inverse of the eigenvector matrix could not */
/*                   be formed due to an attempt to divide by zero, i.e., */
/*                   the eigenvector matrix is singular; */
/*             = N+2:  if the matrix A is defective, possibly due to */
/*                   rounding errors. */

/*     METHOD */

/*     This routine is an implementation of "Method 15" of the set of */
/*     methods described in reference [1], which uses an eigenvalue/ */
/*     eigenvector decomposition technique. A modification of LAPACK */
/*     Library routine DGEEV is used for obtaining the right eigenvector */
/*     matrix. A condition estimate is then employed to determine if the */
/*     matrix A is near defective and hence the exponential solution is */
/*     inaccurate. In this case the routine returns with the Error */
/*     Indicator (INFO) set to N+2, and SLICOT Library routines MB05ND or */
/*     MB05OD are the preferred alternative routines to be used. */

/*     REFERENCES */

/*     [1] Moler, C.B. and Van Loan, C.F. */
/*         Nineteen dubious ways to compute the exponential of a matrix. */
/*         SIAM Review, 20, pp. 801-836, 1978. */

/*     [2] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J., */
/*         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A., */
/*         Ostrouchov, S., and Sorensen, D. */
/*         LAPACK Users' Guide: Second Edition. */
/*         SIAM, Philadelphia, 1995. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997. */
/*     Supersedes Release 2.0 routine MB05AD by M.J. Denham, Kingston */
/*     Polytechnic, March 1981. */

/*     REVISIONS */

/*     V. Sima, June 13, 1997, April 25, 2003, Feb. 15, 2004. */

/*     KEYWORDS */

/*     Eigenvalue, eigenvector decomposition, matrix exponential. */

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

/*     Test the input scalar arguments. */

#line 210 "MB05MD.f"
    /* Parameter adjustments */
#line 210 "MB05MD.f"
    a_dim1 = *lda;
#line 210 "MB05MD.f"
    a_offset = 1 + a_dim1;
#line 210 "MB05MD.f"
    a -= a_offset;
#line 210 "MB05MD.f"
    v_dim1 = *ldv;
#line 210 "MB05MD.f"
    v_offset = 1 + v_dim1;
#line 210 "MB05MD.f"
    v -= v_offset;
#line 210 "MB05MD.f"
    y_dim1 = *ldy;
#line 210 "MB05MD.f"
    y_offset = 1 + y_dim1;
#line 210 "MB05MD.f"
    y -= y_offset;
#line 210 "MB05MD.f"
    --valr;
#line 210 "MB05MD.f"
    --vali;
#line 210 "MB05MD.f"
    --iwork;
#line 210 "MB05MD.f"
    --dwork;
#line 210 "MB05MD.f"

#line 210 "MB05MD.f"
    /* Function Body */
#line 210 "MB05MD.f"
    *info = 0;
#line 211 "MB05MD.f"
    scale = lsame_(balanc, "S", (ftnlen)1, (ftnlen)1);
#line 212 "MB05MD.f"
    if (! (lsame_(balanc, "N", (ftnlen)1, (ftnlen)1) || scale)) {
#line 213 "MB05MD.f"
	*info = -1;
#line 214 "MB05MD.f"
    } else if (*n < 0) {
#line 215 "MB05MD.f"
	*info = -2;
#line 216 "MB05MD.f"
    } else if (*lda < max(1,*n)) {
#line 217 "MB05MD.f"
	*info = -5;
#line 218 "MB05MD.f"
    } else if (*ldv < max(1,*n)) {
#line 219 "MB05MD.f"
	*info = -7;
#line 220 "MB05MD.f"
    } else if (*ldy < max(1,*n)) {
#line 221 "MB05MD.f"
	*info = -9;
#line 222 "MB05MD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 222 "MB05MD.f"
	i__1 = 1, i__2 = *n << 2;
#line 222 "MB05MD.f"
	if (*ldwork < max(i__1,i__2)) {
#line 223 "MB05MD.f"
	    *info = -14;
#line 224 "MB05MD.f"
	}
#line 224 "MB05MD.f"
    }

#line 226 "MB05MD.f"
    if (*info != 0) {

/*        Error return. */

#line 230 "MB05MD.f"
	i__1 = -(*info);
#line 230 "MB05MD.f"
	xerbla_("MB05MD", &i__1, (ftnlen)6);
#line 231 "MB05MD.f"
	return 0;
#line 232 "MB05MD.f"
    }

/*     Quick return if possible. */

#line 236 "MB05MD.f"
    if (*n == 0) {
#line 237 "MB05MD.f"
	dwork[1] = 1.;
#line 238 "MB05MD.f"
	return 0;
#line 239 "MB05MD.f"
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of workspace needed at that point in the code, */
/*     as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

/*     Compute the eigenvalues and right eigenvectors of the real */
/*     nonsymmetric matrix A; optionally, compute a balancing */
/*     transformation. */
/*     Workspace:  need: 4*N. */

#line 252 "MB05MD.f"
    mb05my_(balanc, n, &a[a_offset], lda, &valr[1], &vali[1], &v[v_offset], 
	    ldv, &y[y_offset], ldy, &dwork[1], ldwork, info, (ftnlen)1);

#line 255 "MB05MD.f"
    if (*info > 0) {
#line 255 "MB05MD.f"
	return 0;
#line 255 "MB05MD.f"
    }
#line 257 "MB05MD.f"
    wrkopt = dwork[1];
#line 258 "MB05MD.f"
    if (scale) {
#line 259 "MB05MD.f"
	i__1 = *n;
#line 259 "MB05MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 260 "MB05MD.f"
	    dwork[i__] = dwork[i__ + 1];
#line 261 "MB05MD.f"
/* L10: */
#line 261 "MB05MD.f"
	}
#line 262 "MB05MD.f"
    }

/*     Exit with INFO = N + 1 if V is exactly singular. */

#line 266 "MB05MD.f"
    i__1 = *n;
#line 266 "MB05MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 267 "MB05MD.f"
	if (v[i__ + i__ * v_dim1] == 0.) {
#line 268 "MB05MD.f"
	    *info = *n + 1;
#line 269 "MB05MD.f"
	    return 0;
#line 270 "MB05MD.f"
	}
#line 271 "MB05MD.f"
/* L20: */
#line 271 "MB05MD.f"
    }

/*     Compute the reciprocal condition number of the triangular matrix. */

#line 275 "MB05MD.f"
    dtrcon_("1-norm", "Upper", "Non unit", n, &v[v_offset], ldv, &rcond, &
	    dwork[*n + 1], &iwork[1], info, (ftnlen)6, (ftnlen)5, (ftnlen)8);

/*     Return if the matrix is singular to working precision. */

#line 280 "MB05MD.f"
    if (rcond < dlamch_("Epsilon", (ftnlen)7)) {
#line 281 "MB05MD.f"
	dwork[2] = rcond;
#line 282 "MB05MD.f"
	*info = *n + 2;
#line 283 "MB05MD.f"
	return 0;
#line 284 "MB05MD.f"
    }

/*     Compute the right eigenvector matrix (temporarily) in A. */

#line 288 "MB05MD.f"
    dlacpy_("Full", n, n, &y[y_offset], ldy, &a[a_offset], lda, (ftnlen)4);
#line 289 "MB05MD.f"
    dtrmm_("Right", "Upper", "No transpose", "Non unit", n, n, &c_b16, &v[
	    v_offset], ldv, &a[a_offset], lda, (ftnlen)5, (ftnlen)5, (ftnlen)
	    12, (ftnlen)8);
#line 291 "MB05MD.f"
    if (scale) {
#line 291 "MB05MD.f"
	dgebak_(balanc, "Right", n, &c__1, n, &dwork[1], n, &a[a_offset], lda,
		 info, (ftnlen)1, (ftnlen)5);
#line 291 "MB05MD.f"
    }

/*     Compute the inverse of the right eigenvector matrix, by solving */
/*     a set of linear systems, V * X = Y' (if BALANC = 'N'). */

#line 297 "MB05MD.f"
    i__1 = *n;
#line 297 "MB05MD.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 298 "MB05MD.f"
	i__2 = i__ - 1;
#line 298 "MB05MD.f"
	dswap_(&i__2, &y[i__ + y_dim1], ldy, &y[i__ * y_dim1 + 1], &c__1);
#line 299 "MB05MD.f"
/* L40: */
#line 299 "MB05MD.f"
    }

#line 301 "MB05MD.f"
    dtrsm_("Left", "Upper", "No transpose", "Non unit", n, n, &c_b16, &v[
	    v_offset], ldv, &y[y_offset], ldy, (ftnlen)4, (ftnlen)5, (ftnlen)
	    12, (ftnlen)8);
#line 303 "MB05MD.f"
    if (scale) {

#line 305 "MB05MD.f"
	i__1 = *n;
#line 305 "MB05MD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 306 "MB05MD.f"
	    tempr = 1. / dwork[i__];
#line 307 "MB05MD.f"
	    dscal_(n, &tempr, &y[i__ * y_dim1 + 1], &c__1);
#line 308 "MB05MD.f"
/* L60: */
#line 308 "MB05MD.f"
	}

#line 310 "MB05MD.f"
    }

/*     Save the right eigenvector matrix in V. */

#line 314 "MB05MD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &v[v_offset], ldv, (ftnlen)4);

/*     Premultiply the inverse eigenvector matrix by the exponential of */
/*     quasi-diagonal matrix Lambda * DELTA, where Lambda is the matrix */
/*     of eigenvalues. */
/*     Note that only real arithmetic is used, taking the special storing */
/*     of eigenvalues/eigenvectors into account. */

#line 322 "MB05MD.f"
    i__ = 0;
/*     REPEAT */
#line 324 "MB05MD.f"
L80:
#line 325 "MB05MD.f"
    ++i__;
#line 326 "MB05MD.f"
    if (vali[i__] == 0.) {
#line 327 "MB05MD.f"
	tempr = exp(valr[i__] * *delta);
#line 328 "MB05MD.f"
	dscal_(n, &tempr, &y[i__ + y_dim1], ldy);
#line 329 "MB05MD.f"
    } else {
#line 330 "MB05MD.f"
	tempr = valr[i__] * *delta;
#line 331 "MB05MD.f"
	tempi = vali[i__] * *delta;
#line 332 "MB05MD.f"
	tmp[0] = cos(tempi) * exp(tempr);
#line 333 "MB05MD.f"
	tmp[2] = sin(tempi) * exp(tempr);
#line 334 "MB05MD.f"
	tmp[1] = -tmp[2];
#line 335 "MB05MD.f"
	tmp[3] = tmp[0];
#line 336 "MB05MD.f"
	dlacpy_("Full", &c__2, n, &y[i__ + y_dim1], ldy, &dwork[1], &c__2, (
		ftnlen)4);
#line 337 "MB05MD.f"
	dgemm_("No transpose", "No transpose", &c__2, n, &c__2, &c_b16, tmp, &
		c__2, &dwork[1], &c__2, &c_b40, &y[i__ + y_dim1], ldy, (
		ftnlen)12, (ftnlen)12);
#line 339 "MB05MD.f"
	++i__;
#line 340 "MB05MD.f"
    }
#line 341 "MB05MD.f"
    if (i__ < *n) {
#line 341 "MB05MD.f"
	goto L80;
#line 341 "MB05MD.f"
    }
/*     UNTIL I = N. */

/*     Compute the matrix exponential as the product V * Y. */

#line 346 "MB05MD.f"
    dgemm_("No transpose", "No transpose", n, n, n, &c_b16, &v[v_offset], ldv,
	     &y[y_offset], ldy, &c_b40, &a[a_offset], lda, (ftnlen)12, (
	    ftnlen)12);

/*     Set optimal workspace dimension and reciprocal condition number. */

#line 351 "MB05MD.f"
    dwork[1] = wrkopt;
#line 352 "MB05MD.f"
    dwork[2] = rcond;

#line 354 "MB05MD.f"
    return 0;
/* *** Last line of MB05MD *** */
} /* mb05md_ */

