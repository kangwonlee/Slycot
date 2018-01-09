#line 1 "TB01VY.f"
/* TB01VY.f -- translated by f2c (version 20100827).
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

#line 1 "TB01VY.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b18 = -1.;
static doublereal c_b20 = 0.;
static doublereal c_b24 = .5;
static doublereal c_b27 = 1.;

/* Subroutine */ int tb01vy_(char *apply, integer *n, integer *m, integer *l, 
	doublereal *theta, integer *ltheta, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *x0, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen apply_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double atan(doublereal), sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, ca, in;
    static doublereal ri, ti;
    static integer ldca;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), daxpy_(integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *)
	    ;
    static integer jwork;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static doublereal factor;
    static logical lapply;
    static doublereal tobypi;


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

/*     To convert the linear discrete-time system given as its output */
/*     normal form [1], with parameter vector THETA, into the state-space */
/*     representation (A, B, C, D), with the initial state x0. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     APPLY   CHARACTER*1 */
/*             Specifies whether or not the parameter vector should be */
/*             transformed using a bijective mapping, as follows: */
/*             = 'A' : apply the bijective mapping to the N vectors in */
/*                     THETA corresponding to the matrices A and C; */
/*             = 'N' : do not apply the bijective mapping. */
/*             The transformation performed when APPLY = 'A' allows */
/*             to get rid of the constraints norm(THETAi) < 1, i = 1:N. */
/*             A call of the SLICOT Library routine TB01VD associated to */
/*             a call of TB01VY must use the same value of APPLY. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     L       (input) INTEGER */
/*             The number of system outputs.  L >= 0. */

/*     THETA   (input) DOUBLE PRECISION array, dimension (LTHETA) */
/*             The leading N*(L+M+1)+L*M part of this array must contain */
/*             the parameter vector that defines a system (A, B, C, D), */
/*             with the initial state x0. The parameters are: */

/*             THETA(1:N*L)                      : parameters for A, C; */
/*             THETA(N*L+1:N*(L+M))              : parameters for B; */
/*             THETA(N*(L+M)+1:N*(L+M)+L*M)      : parameters for D; */
/*             THETA(N*(L+M)+L*M+1:N*(L+M+1)+L*M): parameters for x0. */

/*     LTHETA  INTEGER */
/*             The length of array THETA.  LTHETA >= N*(L+M+1)+L*M. */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array contains the system */
/*             state matrix corresponding to the output normal form with */
/*             parameter vector THETA. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array contains the system */
/*             input matrix corresponding to the output normal form with */
/*             parameter vector THETA. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading L-by-N part of this array contains the system */
/*             output matrix corresponding to the output normal form with */
/*             parameter vector THETA. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,L). */

/*     D       (output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading L-by-M part of this array contains the system */
/*             input/output matrix corresponding to the output normal */
/*             form with parameter vector THETA. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,L). */

/*     X0      (output) DOUBLE PRECISION array, dimension (N) */
/*             This array contains the initial state of the system, x0, */
/*             corresponding to the output normal form with parameter */
/*             vector THETA. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= N*(N+L+1). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The parameters characterizing A and C are used to build N */
/*     orthogonal transformations, which are then applied to recover */
/*     these matrices. */

/*     CONTRIBUTORS */

/*     A. Riedel, R. Schneider, Chemnitz University of Technology, */
/*     Oct. 2000, during a stay at University of Twente, NL. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001, */
/*     Feb. 2002, Feb. 2004. */

/*     KEYWORDS */

/*     Asymptotically stable, output normal form, parameter estimation, */
/*     similarity transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

#line 170 "TB01VY.f"
    /* Parameter adjustments */
#line 170 "TB01VY.f"
    --theta;
#line 170 "TB01VY.f"
    a_dim1 = *lda;
#line 170 "TB01VY.f"
    a_offset = 1 + a_dim1;
#line 170 "TB01VY.f"
    a -= a_offset;
#line 170 "TB01VY.f"
    b_dim1 = *ldb;
#line 170 "TB01VY.f"
    b_offset = 1 + b_dim1;
#line 170 "TB01VY.f"
    b -= b_offset;
#line 170 "TB01VY.f"
    c_dim1 = *ldc;
#line 170 "TB01VY.f"
    c_offset = 1 + c_dim1;
#line 170 "TB01VY.f"
    c__ -= c_offset;
#line 170 "TB01VY.f"
    d_dim1 = *ldd;
#line 170 "TB01VY.f"
    d_offset = 1 + d_dim1;
#line 170 "TB01VY.f"
    d__ -= d_offset;
#line 170 "TB01VY.f"
    --x0;
#line 170 "TB01VY.f"
    --dwork;
#line 170 "TB01VY.f"

#line 170 "TB01VY.f"
    /* Function Body */
#line 170 "TB01VY.f"
    lapply = lsame_(apply, "A", (ftnlen)1, (ftnlen)1);

#line 172 "TB01VY.f"
    *info = 0;
#line 173 "TB01VY.f"
    if (! (lapply || lsame_(apply, "N", (ftnlen)1, (ftnlen)1))) {
#line 174 "TB01VY.f"
	*info = -1;
#line 175 "TB01VY.f"
    } else if (*n < 0) {
#line 176 "TB01VY.f"
	*info = -2;
#line 177 "TB01VY.f"
    } else if (*m < 0) {
#line 178 "TB01VY.f"
	*info = -3;
#line 179 "TB01VY.f"
    } else if (*l < 0) {
#line 180 "TB01VY.f"
	*info = -4;
#line 181 "TB01VY.f"
    } else if (*ltheta < *n * (*l + *m + 1) + *l * *m) {
#line 182 "TB01VY.f"
	*info = -6;
#line 183 "TB01VY.f"
    } else if (*lda < max(1,*n)) {
#line 184 "TB01VY.f"
	*info = -8;
#line 185 "TB01VY.f"
    } else if (*ldb < max(1,*n)) {
#line 186 "TB01VY.f"
	*info = -10;
#line 187 "TB01VY.f"
    } else if (*ldc < max(1,*l)) {
#line 188 "TB01VY.f"
	*info = -12;
#line 189 "TB01VY.f"
    } else if (*ldd < max(1,*l)) {
#line 190 "TB01VY.f"
	*info = -14;
#line 191 "TB01VY.f"
    } else if (*ldwork < *n * (*n + *l + 1)) {
#line 192 "TB01VY.f"
	*info = -17;
#line 193 "TB01VY.f"
    }

/*     Return if there are illegal arguments. */

#line 197 "TB01VY.f"
    if (*info != 0) {
#line 198 "TB01VY.f"
	i__1 = -(*info);
#line 198 "TB01VY.f"
	xerbla_("TB01VY", &i__1, (ftnlen)6);
#line 199 "TB01VY.f"
	return 0;
#line 200 "TB01VY.f"
    }

/*     Quick return if possible. */

/* Computing MAX */
#line 204 "TB01VY.f"
    i__1 = max(*n,*m);
#line 204 "TB01VY.f"
    if (max(i__1,*l) == 0) {
#line 204 "TB01VY.f"
	return 0;
#line 204 "TB01VY.f"
    }

#line 207 "TB01VY.f"
    if (*m > 0) {

/*        Copy the matrix B from THETA. */

#line 211 "TB01VY.f"
	dlacpy_("Full", n, m, &theta[*n * *l + 1], n, &b[b_offset], ldb, (
		ftnlen)4);

/*        Copy the matrix D. */

#line 215 "TB01VY.f"
	dlacpy_("Full", l, m, &theta[*n * (*l + *m) + 1], l, &d__[d_offset], 
		ldd, (ftnlen)4);
#line 216 "TB01VY.f"
    }

#line 218 "TB01VY.f"
    if (*n == 0) {
#line 219 "TB01VY.f"
	return 0;
#line 220 "TB01VY.f"
    } else if (*l == 0) {
#line 221 "TB01VY.f"
	dcopy_(n, &theta[*n * *m + 1], &c__1, &x0[1], &c__1);
#line 222 "TB01VY.f"
	return 0;
#line 223 "TB01VY.f"
    }

/*     Initialize the indices in the workspace. */

#line 227 "TB01VY.f"
    ldca = *n + *l;

#line 229 "TB01VY.f"
    ca = 1;

#line 231 "TB01VY.f"
    jwork = ca + *n * ldca;
#line 232 "TB01VY.f"
    tobypi = .5 / atan(1.);

/*     Generate the matrices C and A from their parameters. */
/*     Start with the block matrix [0; I], where 0 is a block of zeros */
/*     of size L-by-N, and I is the identity matrix of order N. */

#line 238 "TB01VY.f"
    dwork[ca] = 0.;
#line 239 "TB01VY.f"
    i__1 = *n * (*l + *n);
#line 239 "TB01VY.f"
    dcopy_(&i__1, &dwork[ca], &c__0, &dwork[ca], &c__1);
#line 240 "TB01VY.f"
    dwork[ca + *l] = 1.;
#line 241 "TB01VY.f"
    i__1 = ldca + 1;
#line 241 "TB01VY.f"
    dcopy_(n, &dwork[ca + *l], &c__0, &dwork[ca + *l], &i__1);

/*     Now, read out THETA(1 : N*L) and perform the transformations */
/*     defined by the parameters in THETA. */

#line 246 "TB01VY.f"
    for (i__ = *n; i__ >= 1; --i__) {

/*        Save THETAi in the first column of C and use the copy for */
/*        further processing. */

#line 251 "TB01VY.f"
	dcopy_(l, &theta[(i__ - 1) * *l + 1], &c__1, &c__[c_offset], &c__1);
#line 252 "TB01VY.f"
	ti = dnrm2_(l, &c__[c_offset], &c__1);
#line 253 "TB01VY.f"
	if (lapply && ti != 0.) {

/*           Apply the bijective mapping which guarantees that TI < 1. */

#line 257 "TB01VY.f"
	    factor = tobypi * atan(ti) / ti;

/*           Scale THETAi and apply the same scaling on TI. */

#line 261 "TB01VY.f"
	    dscal_(l, &factor, &c__[c_offset], &c__1);
#line 262 "TB01VY.f"
	    ti *= factor;
#line 263 "TB01VY.f"
	}

/*        RI = sqrt( 1 - TI**2 ). */

#line 267 "TB01VY.f"
	ri = sqrt((1. - ti) * (ti + 1.));

/*        Multiply a certain part of DWORK(CA) with Ui' from the left, */
/*        where Ui = [ -THETAi, Si; RI, THETAi' ] is (L+1)-by-(L+1), but */
/*        Ui is not stored. */

#line 273 "TB01VY.f"
	dgemv_("Transpose", l, n, &c_b18, &dwork[ca + *n - i__], &ldca, &c__[
		c_offset], &c__1, &c_b20, &dwork[jwork], &c__1, (ftnlen)9);

#line 276 "TB01VY.f"
	if (ti > 0.) {
#line 277 "TB01VY.f"
	    d__1 = (1. - ri) / ti / ti;
#line 277 "TB01VY.f"
	    dger_(l, n, &d__1, &c__[c_offset], &c__1, &dwork[jwork], &c__1, &
		    dwork[ca + *n - i__], &ldca);
#line 279 "TB01VY.f"
	} else {

/*           The call below is for the limiting case. */

#line 283 "TB01VY.f"
	    dger_(l, n, &c_b24, &c__[c_offset], &c__1, &dwork[jwork], &c__1, &
		    dwork[ca + *n - i__], &ldca);
#line 285 "TB01VY.f"
	}

#line 287 "TB01VY.f"
	dger_(l, n, &c_b27, &c__[c_offset], &c__1, &dwork[ca + *n - i__ + *l],
		 &ldca, &dwork[ca + *n - i__], &ldca);
#line 289 "TB01VY.f"
	daxpy_(n, &ri, &dwork[ca + *n - i__ + *l], &ldca, &dwork[jwork], &
		c__1);

/*        Move these results to their appropriate locations. */

#line 293 "TB01VY.f"
	i__1 = *n;
#line 293 "TB01VY.f"
	for (j = 1; j <= i__1; ++j) {
#line 294 "TB01VY.f"
	    in = ca + *n - i__ + (j - 1) * ldca;
#line 295 "TB01VY.f"
	    i__2 = in + 1;
#line 295 "TB01VY.f"
	    for (k = in + *l; k >= i__2; --k) {
#line 296 "TB01VY.f"
		dwork[k] = dwork[k - 1];
#line 297 "TB01VY.f"
/* L10: */
#line 297 "TB01VY.f"
	    }
#line 298 "TB01VY.f"
	    dwork[in] = dwork[jwork + j - 1];
#line 299 "TB01VY.f"
/* L20: */
#line 299 "TB01VY.f"
	}

#line 301 "TB01VY.f"
/* L30: */
#line 301 "TB01VY.f"
    }

/*     Now, DWORK(CA) = [C; A]. Copy to C and A. */

#line 305 "TB01VY.f"
    i__1 = *n;
#line 305 "TB01VY.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 306 "TB01VY.f"
	dcopy_(l, &dwork[ca + (i__ - 1) * ldca], &c__1, &c__[i__ * c_dim1 + 1]
		, &c__1);
#line 307 "TB01VY.f"
	dcopy_(n, &dwork[ca + *l + (i__ - 1) * ldca], &c__1, &a[i__ * a_dim1 
		+ 1], &c__1);
#line 308 "TB01VY.f"
/* L40: */
#line 308 "TB01VY.f"
    }

/*     Copy the initial state x0. */

#line 312 "TB01VY.f"
    dcopy_(n, &theta[*n * (*l + *m) + *l * *m + 1], &c__1, &x0[1], &c__1);

#line 314 "TB01VY.f"
    return 0;

/* *** Last line of TB01VY *** */
} /* tb01vy_ */

