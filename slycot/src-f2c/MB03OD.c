#line 1 "MB03OD.f"
/* MB03OD.f -- translated by f2c (version 20100827).
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

#line 1 "MB03OD.f"
/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;

/* Subroutine */ int mb03od_(char *jobqr, integer *m, integer *n, doublereal *
	a, integer *lda, integer *jpvt, doublereal *rcond, doublereal *svlmax,
	 doublereal *tau, integer *rank, doublereal *sval, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen jobqr_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal c1, c2, s1, s2;
    static integer mn;
    static doublereal smin, smax;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer ismin, ismax;
    extern /* Subroutine */ int dlaic1_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), dgeqp3_(integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *), xerbla_(char *, integer *, ftnlen);
    static logical ljobqr;
    static integer minwrk;
    static doublereal sminpr;
    static integer maxwrk;
    static doublereal smaxpr;


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

/*     To compute (optionally) a rank-revealing QR factorization of a */
/*     real general M-by-N matrix  A,  which may be rank-deficient, */
/*     and estimate its effective rank using incremental condition */
/*     estimation. */

/*     The routine uses a QR factorization with column pivoting: */
/*        A * P = Q * R,  where  R = [ R11 R12 ], */
/*                                   [  0  R22 ] */
/*     with R11 defined as the largest leading submatrix whose estimated */
/*     condition number is less than 1/RCOND.  The order of R11, RANK, */
/*     is the effective rank of A. */

/*     MB03OD  does not perform any scaling of the matrix A. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBQR   CHARACTER*1 */
/*             = 'Q':  Perform a QR factorization with column pivoting; */
/*             = 'N':  Do not perform the QR factorization (but assume */
/*                     that it has been done outside). */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension */
/*             ( LDA, N ) */
/*             On entry with JOBQR = 'Q', the leading M by N part of this */
/*             array must contain the given matrix A. */
/*             On exit with JOBQR = 'Q', the leading min(M,N) by N upper */
/*             triangular part of A contains the triangular factor R, */
/*             and the elements below the diagonal, with the array TAU, */
/*             represent the orthogonal matrix Q as a product of */
/*             min(M,N) elementary reflectors. */
/*             On entry and on exit with JOBQR = 'N', the leading */
/*             min(M,N) by N upper triangular part of A contains the */
/*             triangular factor R, as determined by the QR factorization */
/*             with pivoting.  The elements below the diagonal of A are */
/*             not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     JPVT    (input/output) INTEGER array, dimension ( N ) */
/*             On entry with JOBQR = 'Q', if JPVT(i) <> 0, the i-th */
/*             column of A is an initial column, otherwise it is a free */
/*             column. Before the QR factorization of A, all initial */
/*             columns are permuted to the leading positions; only the */
/*             remaining free columns are moved as a result of column */
/*             pivoting during the factorization.  For rank determination */
/*             it is preferable that all columns be free. */
/*             On exit with JOBQR = 'Q', if JPVT(i) = k, then the i-th */
/*             column of A*P was the k-th column of A. */
/*             Array JPVT is not referenced when JOBQR = 'N'. */

/*     RCOND   (input) DOUBLE PRECISION */
/*             RCOND is used to determine the effective rank of A, which */
/*             is defined as the order of the largest leading triangular */
/*             submatrix R11 in the QR factorization with pivoting of A, */
/*             whose estimated condition number is less than 1/RCOND. */
/*             RCOND >= 0. */
/*             NOTE that when SVLMAX > 0, the estimated rank could be */
/*             less than that defined above (see SVLMAX). */

/*     SVLMAX  (input) DOUBLE PRECISION */
/*             If A is a submatrix of another matrix B, and the rank */
/*             decision should be related to that matrix, then SVLMAX */
/*             should be an estimate of the largest singular value of B */
/*             (for instance, the Frobenius norm of B).  If this is not */
/*             the case, the input value SVLMAX = 0 should work. */
/*             SVLMAX >= 0. */

/*     TAU     (output) DOUBLE PRECISION array, dimension ( MIN( M, N ) ) */
/*             On exit with JOBQR = 'Q', the leading min(M,N) elements of */
/*             TAU contain the scalar factors of the elementary */
/*             reflectors. */
/*             Array TAU is not referenced when JOBQR = 'N'. */

/*     RANK    (output) INTEGER */
/*             The effective (estimated) rank of A, i.e. the order of */
/*             the submatrix R11. */

/*     SVAL    (output) DOUBLE PRECISION array, dimension ( 3 ) */
/*             The estimates of some of the singular values of the */
/*             triangular factor R: */
/*             SVAL(1): largest singular value of R(1:RANK,1:RANK); */
/*             SVAL(2): smallest singular value of R(1:RANK,1:RANK); */
/*             SVAL(3): smallest singular value of R(1:RANK+1,1:RANK+1), */
/*                      if RANK < MIN( M, N ), or of R(1:RANK,1:RANK), */
/*                      otherwise. */
/*             If the triangular factorization is a rank-revealing one */
/*             (which will be the case if the leading columns were well- */
/*             conditioned), then SVAL(1) will also be an estimate for */
/*             the largest singular value of A, and SVAL(2) and SVAL(3) */
/*             will be estimates for the RANK-th and (RANK+1)-st singular */
/*             values of A, respectively. */
/*             By examining these values, one can confirm that the rank */
/*             is well defined with respect to the chosen value of RCOND. */
/*             The ratio SVAL(1)/SVAL(2) is an estimate of the condition */
/*             number of R(1:RANK,1:RANK). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension ( LDWORK ) */
/*             On exit, if  INFO = 0,  DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 3*N + 1,                 if JOBQR = 'Q'; */
/*             LDWORK >= max( 1, 2*min( M, N ) ), if JOBQR = 'N'. */
/*             For good performance when JOBQR = 'Q', LDWORK should be */
/*             larger. Specifically, LDWORK >= 2*N + ( N + 1 )*NB, where */
/*             NB is the optimal block size for the LAPACK Library */
/*             routine DGEQP3. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine computes or uses a QR factorization with column */
/*     pivoting of A,  A * P = Q * R,  with  R  defined above, and then */
/*     finds the largest leading submatrix whose estimated condition */
/*     number is less than 1/RCOND, taking the possible positive value of */
/*     SVLMAX into account.  This is performed using the LAPACK */
/*     incremental condition estimation scheme and a slightly modified */
/*     rank decision test. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2005. */

/*    ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 200 "MB03OD.f"
    /* Parameter adjustments */
#line 200 "MB03OD.f"
    a_dim1 = *lda;
#line 200 "MB03OD.f"
    a_offset = 1 + a_dim1;
#line 200 "MB03OD.f"
    a -= a_offset;
#line 200 "MB03OD.f"
    --jpvt;
#line 200 "MB03OD.f"
    --tau;
#line 200 "MB03OD.f"
    --sval;
#line 200 "MB03OD.f"
    --dwork;
#line 200 "MB03OD.f"

#line 200 "MB03OD.f"
    /* Function Body */
#line 200 "MB03OD.f"
    ljobqr = lsame_(jobqr, "Q", (ftnlen)1, (ftnlen)1);
#line 201 "MB03OD.f"
    mn = min(*m,*n);
#line 202 "MB03OD.f"
    ismin = 1;
#line 203 "MB03OD.f"
    ismax = mn + 1;
#line 204 "MB03OD.f"
    if (ljobqr) {
#line 205 "MB03OD.f"
	minwrk = *n * 3 + 1;
#line 206 "MB03OD.f"
    } else {
/* Computing MAX */
#line 207 "MB03OD.f"
	i__1 = 1, i__2 = mn << 1;
#line 207 "MB03OD.f"
	minwrk = max(i__1,i__2);
#line 208 "MB03OD.f"
    }
#line 209 "MB03OD.f"
    maxwrk = minwrk;

/*     Test the input scalar arguments. */

#line 213 "MB03OD.f"
    *info = 0;
#line 214 "MB03OD.f"
    if (! ljobqr && ! lsame_(jobqr, "N", (ftnlen)1, (ftnlen)1)) {
#line 215 "MB03OD.f"
	*info = -1;
#line 216 "MB03OD.f"
    } else if (*m < 0) {
#line 217 "MB03OD.f"
	*info = -2;
#line 218 "MB03OD.f"
    } else if (*n < 0) {
#line 219 "MB03OD.f"
	*info = -3;
#line 220 "MB03OD.f"
    } else if (*lda < max(1,*m)) {
#line 221 "MB03OD.f"
	*info = -5;
#line 222 "MB03OD.f"
    } else if (*rcond < 0.) {
#line 223 "MB03OD.f"
	*info = -7;
#line 224 "MB03OD.f"
    } else if (*svlmax < 0.) {
#line 225 "MB03OD.f"
	*info = -8;
#line 226 "MB03OD.f"
    } else if (*ldwork < minwrk) {
#line 227 "MB03OD.f"
	*info = -13;
#line 228 "MB03OD.f"
    }

#line 230 "MB03OD.f"
    if (*info != 0) {
#line 231 "MB03OD.f"
	i__1 = -(*info);
#line 231 "MB03OD.f"
	xerbla_("MB03OD", &i__1, (ftnlen)6);
#line 232 "MB03OD.f"
	return 0;
#line 233 "MB03OD.f"
    }

/*     Quick return if possible */

#line 237 "MB03OD.f"
    if (mn == 0) {
#line 238 "MB03OD.f"
	*rank = 0;
#line 239 "MB03OD.f"
	sval[1] = 0.;
#line 240 "MB03OD.f"
	sval[2] = 0.;
#line 241 "MB03OD.f"
	sval[3] = 0.;
#line 242 "MB03OD.f"
	dwork[1] = 1.;
#line 243 "MB03OD.f"
	return 0;
#line 244 "MB03OD.f"
    }

#line 246 "MB03OD.f"
    if (ljobqr) {

/*        Compute QR factorization with column pivoting of A: */
/*           A * P = Q * R */
/*        Workspace need   3*N + 1; */
/*                  prefer 2*N + (N+1)*NB. */
/*        Details of Householder rotations stored in TAU. */

#line 254 "MB03OD.f"
	dgeqp3_(m, n, &a[a_offset], lda, &jpvt[1], &tau[1], &dwork[1], ldwork,
		 info);
/* Computing MAX */
#line 255 "MB03OD.f"
	i__1 = maxwrk, i__2 = (integer) dwork[1];
#line 255 "MB03OD.f"
	maxwrk = max(i__1,i__2);
#line 256 "MB03OD.f"
    }

/*     Determine RANK using incremental condition estimation */

#line 260 "MB03OD.f"
    dwork[ismin] = 1.;
#line 261 "MB03OD.f"
    dwork[ismax] = 1.;
#line 262 "MB03OD.f"
    smax = (d__1 = a[a_dim1 + 1], abs(d__1));
#line 263 "MB03OD.f"
    smin = smax;
#line 264 "MB03OD.f"
    if (smax == 0. || *svlmax * *rcond > smax) {
#line 265 "MB03OD.f"
	*rank = 0;
#line 266 "MB03OD.f"
	sval[1] = smax;
#line 267 "MB03OD.f"
	sval[2] = 0.;
#line 268 "MB03OD.f"
	sval[3] = 0.;
#line 269 "MB03OD.f"
    } else {
#line 270 "MB03OD.f"
	*rank = 1;
#line 271 "MB03OD.f"
	sminpr = smin;

#line 273 "MB03OD.f"
L10:
#line 274 "MB03OD.f"
	if (*rank < mn) {
#line 275 "MB03OD.f"
	    i__ = *rank + 1;
#line 276 "MB03OD.f"
	    dlaic1_(&c__2, rank, &dwork[ismin], &smin, &a[i__ * a_dim1 + 1], &
		    a[i__ + i__ * a_dim1], &sminpr, &s1, &c1);
#line 278 "MB03OD.f"
	    dlaic1_(&c__1, rank, &dwork[ismax], &smax, &a[i__ * a_dim1 + 1], &
		    a[i__ + i__ * a_dim1], &smaxpr, &s2, &c2);

#line 281 "MB03OD.f"
	    if (*svlmax * *rcond <= smaxpr) {
#line 282 "MB03OD.f"
		if (*svlmax * *rcond <= sminpr) {
#line 283 "MB03OD.f"
		    if (smaxpr * *rcond <= sminpr) {
#line 284 "MB03OD.f"
			i__1 = *rank;
#line 284 "MB03OD.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 285 "MB03OD.f"
			    dwork[ismin + i__ - 1] = s1 * dwork[ismin + i__ - 
				    1];
#line 286 "MB03OD.f"
			    dwork[ismax + i__ - 1] = s2 * dwork[ismax + i__ - 
				    1];
#line 287 "MB03OD.f"
/* L20: */
#line 287 "MB03OD.f"
			}
#line 288 "MB03OD.f"
			dwork[ismin + *rank] = c1;
#line 289 "MB03OD.f"
			dwork[ismax + *rank] = c2;
#line 290 "MB03OD.f"
			smin = sminpr;
#line 291 "MB03OD.f"
			smax = smaxpr;
#line 292 "MB03OD.f"
			++(*rank);
#line 293 "MB03OD.f"
			goto L10;
#line 294 "MB03OD.f"
		    }
#line 295 "MB03OD.f"
		}
#line 296 "MB03OD.f"
	    }
#line 297 "MB03OD.f"
	}
#line 298 "MB03OD.f"
	sval[1] = smax;
#line 299 "MB03OD.f"
	sval[2] = smin;
#line 300 "MB03OD.f"
	sval[3] = sminpr;
#line 301 "MB03OD.f"
    }

#line 303 "MB03OD.f"
    dwork[1] = (doublereal) maxwrk;
#line 304 "MB03OD.f"
    return 0;
/* *** Last line of MB03OD *** */
} /* mb03od_ */

