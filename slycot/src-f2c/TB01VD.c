#line 1 "TB01VD.f"
/* TB01VD.f -- translated by f2c (version 20100827).
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

#line 1 "TB01VD.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b18 = 1.;
static doublereal c_b42 = 0.;
static doublereal c_b93 = -.5;
static doublereal c_b96 = -1.;

/* Subroutine */ int tb01vd_(char *apply, integer *n, integer *m, integer *l, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *d__, integer *ldd, doublereal *x0, 
	doublereal *theta, integer *ltheta, doublereal *dwork, integer *
	ldwork, integer *info, ftnlen apply_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Builtin functions */
    double atan(doublereal), tan(doublereal);

    /* Local variables */
    static integer i__, j, k, ca, ia, in, iq;
    static doublereal ri;
    static integer ir;
    static doublereal ti;
    static integer it, iu, ldt, iwi, iwr, ldca;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer itau;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal piby2;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), dscal_(
	    integer *, doublereal *, doublereal *, integer *);
    static doublereal scale;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     sb03od_(char *, char *, char *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), dtrmm_(char *, 
	    char *, char *, char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen), daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dtrsm_(char *, char *, char *
	    , char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    static integer jwork;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen), dgeqrf_(integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *), dlacpy_(char *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical lapply;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
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

/*     To convert the linear discrete-time system given as (A, B, C, D), */
/*     with initial state x0, into the output normal form [1], with */
/*     parameter vector THETA. The matrix A is assumed to be stable. */
/*     The matrices A, B, C, D and the vector x0 are converted, so that */
/*     on exit they correspond to the system defined by THETA. */

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
/*             A call of the SLICOT Library routine TB01VY associated to */
/*             a call of TB01VD must use the same value of APPLY. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     L       (input) INTEGER */
/*             The number of system outputs.  L >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the system state matrix A, assumed to be stable. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the transformed system state matrix corresponding to the */
/*             output normal form with parameter vector THETA. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the system input matrix B. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the transformed system input matrix corresponding to the */
/*             output normal form with parameter vector THETA. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the system output matrix C. */
/*             On exit, the leading L-by-N part of this array contains */
/*             the transformed system output matrix corresponding to the */
/*             output normal form with parameter vector THETA. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,L). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading L-by-M part of this array must contain the */
/*             system input/output matrix D. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,L). */

/*     X0      (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the initial state of the */
/*             system, x0. */
/*             On exit, this array contains the transformed initial state */
/*             of the system, corresponding to the output normal form */
/*             with parameter vector THETA. */

/*     THETA   (output) DOUBLE PRECISION array, dimension (LTHETA) */
/*             The leading N*(L+M+1)+L*M part of this array contains the */
/*             parameter vector that defines a system (A, B, C, D, x0) */
/*             which is equivalent up to a similarity transformation to */
/*             the system given on entry. The parameters are: */

/*             THETA(1:N*L)                      : parameters for A, C; */
/*             THETA(N*L+1:N*(L+M))              : parameters for B; */
/*             THETA(N*(L+M)+1:N*(L+M)+L*M)      : parameters for D; */
/*             THETA(N*(L+M)+L*M+1:N*(L+M+1)+L*M): parameters for x0. */

/*     LTHETA  INTEGER */
/*             The length of array THETA.  LTHETA >= N*(L+M+1)+L*M. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1, N*N*L + N*L + N, */
/*                           N*N + MAX(N*N + N*MAX(N,L) + 6*N + MIN(N,L), */
/*                                     N*M)). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the Lyapunov equation A'*Q*A - Q = -scale^2*C'*C */
/*                   could only be solved with scale = 0; */
/*             = 2:  if matrix A is not discrete-time stable; */
/*             = 3:  if the QR algorithm failed to converge for */
/*                   matrix A. */

/*     METHOD */

/*     The matrices A and C are converted to output normal form. */
/*     First, the Lyapunov equation */

/*        A'*Q*A - Q = -scale^2*C'*C, */

/*     is solved in the Cholesky factor T, T'*T = Q, and then T is used */
/*     to get the transformation matrix. */

/*     The matrix B and the initial state x0 are transformed accordingly. */

/*     Then, the QR factorization of the transposed observability matrix */
/*     is computed, and the matrix Q is used to further transform the */
/*     system matrices. The parameters characterizing A and C are finally */
/*     obtained by applying a set of N orthogonal transformations. */

/*     REFERENCES */

/*     [1] Peeters, R.L.M., Hanzon, B., and Olivi, M. */
/*         Balanced realizations of discrete-time stable all-pass */
/*         systems and the tangential Schur algorithm. */
/*         Proceedings of the European Control Conference, */
/*         31 August - 3 September 1999, Karlsruhe, Germany. */
/*         Session CP-6, Discrete-time Systems, 1999. */

/*     CONTRIBUTORS */

/*     A. Riedel, R. Schneider, Chemnitz University of Technology, */
/*     Oct. 2000, during a stay at University of Twente, NL. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001, */
/*     Feb. 2002, Feb. 2004. */

/*     KEYWORDS */

/*     Asymptotically stable, Lyapunov equation, output normal form, */
/*     parameter estimation, similarity transformation. */

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

#line 213 "TB01VD.f"
    /* Parameter adjustments */
#line 213 "TB01VD.f"
    a_dim1 = *lda;
#line 213 "TB01VD.f"
    a_offset = 1 + a_dim1;
#line 213 "TB01VD.f"
    a -= a_offset;
#line 213 "TB01VD.f"
    b_dim1 = *ldb;
#line 213 "TB01VD.f"
    b_offset = 1 + b_dim1;
#line 213 "TB01VD.f"
    b -= b_offset;
#line 213 "TB01VD.f"
    c_dim1 = *ldc;
#line 213 "TB01VD.f"
    c_offset = 1 + c_dim1;
#line 213 "TB01VD.f"
    c__ -= c_offset;
#line 213 "TB01VD.f"
    d_dim1 = *ldd;
#line 213 "TB01VD.f"
    d_offset = 1 + d_dim1;
#line 213 "TB01VD.f"
    d__ -= d_offset;
#line 213 "TB01VD.f"
    --x0;
#line 213 "TB01VD.f"
    --theta;
#line 213 "TB01VD.f"
    --dwork;
#line 213 "TB01VD.f"

#line 213 "TB01VD.f"
    /* Function Body */
#line 213 "TB01VD.f"
    lapply = lsame_(apply, "A", (ftnlen)1, (ftnlen)1);

#line 215 "TB01VD.f"
    *info = 0;
#line 216 "TB01VD.f"
    if (! (lapply || lsame_(apply, "N", (ftnlen)1, (ftnlen)1))) {
#line 217 "TB01VD.f"
	*info = -1;
#line 218 "TB01VD.f"
    } else if (*n < 0) {
#line 219 "TB01VD.f"
	*info = -2;
#line 220 "TB01VD.f"
    } else if (*m < 0) {
#line 221 "TB01VD.f"
	*info = -3;
#line 222 "TB01VD.f"
    } else if (*l < 0) {
#line 223 "TB01VD.f"
	*info = -4;
#line 224 "TB01VD.f"
    } else if (*lda < max(1,*n)) {
#line 225 "TB01VD.f"
	*info = -6;
#line 226 "TB01VD.f"
    } else if (*ldb < max(1,*n)) {
#line 227 "TB01VD.f"
	*info = -8;
#line 228 "TB01VD.f"
    } else if (*ldc < max(1,*l)) {
#line 229 "TB01VD.f"
	*info = -10;
#line 230 "TB01VD.f"
    } else if (*ldd < max(1,*l)) {
#line 231 "TB01VD.f"
	*info = -12;
#line 232 "TB01VD.f"
    } else if (*ltheta < *n * (*m + *l + 1) + *l * *m) {
#line 233 "TB01VD.f"
	*info = -15;
#line 234 "TB01VD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
#line 234 "TB01VD.f"
	i__3 = *n * (*n + max(*n,*l) + 6) + min(*n,*l), i__4 = *n * *m;
#line 234 "TB01VD.f"
	i__1 = 1, i__2 = *n * *n * *l + *n * *l + *n, i__1 = max(i__1,i__2), 
		i__2 = *n * *n + max(i__3,i__4);
#line 234 "TB01VD.f"
	if (*ldwork < max(i__1,i__2)) {
#line 237 "TB01VD.f"
	    *info = -17;
#line 238 "TB01VD.f"
	}
#line 238 "TB01VD.f"
    }

/*     Return if there are illegal arguments. */

#line 242 "TB01VD.f"
    if (*info != 0) {
#line 243 "TB01VD.f"
	i__1 = -(*info);
#line 243 "TB01VD.f"
	xerbla_("TB01VD", &i__1, (ftnlen)6);
#line 244 "TB01VD.f"
	return 0;
#line 245 "TB01VD.f"
    }

/*     Quick return if possible. */

/* Computing MAX */
#line 249 "TB01VD.f"
    i__1 = max(*n,*m);
#line 249 "TB01VD.f"
    if (max(i__1,*l) == 0) {
#line 250 "TB01VD.f"
	dwork[1] = 1.;
#line 251 "TB01VD.f"
	return 0;
#line 252 "TB01VD.f"
    } else if (*n == 0) {
#line 253 "TB01VD.f"
	i__1 = max(1,*l);
#line 253 "TB01VD.f"
	dlacpy_("Full", l, m, &d__[d_offset], ldd, &theta[1], &i__1, (ftnlen)
		4);
#line 254 "TB01VD.f"
	dwork[1] = 1.;
#line 255 "TB01VD.f"
	return 0;
#line 256 "TB01VD.f"
    } else if (*l == 0) {
#line 257 "TB01VD.f"
	dlacpy_("Full", n, m, &b[b_offset], ldb, &theta[1], n, (ftnlen)4);
#line 258 "TB01VD.f"
	dcopy_(n, &x0[1], &c__1, &theta[*n * *m + 1], &c__1);
#line 259 "TB01VD.f"
	dwork[1] = 1.;
#line 260 "TB01VD.f"
	return 0;
#line 261 "TB01VD.f"
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 269 "TB01VD.f"
    wrkopt = 1;
#line 270 "TB01VD.f"
    piby2 = atan(1.) * 2.;

/*     Convert A and C to output normal form. */
/*     First, solve the Lyapunov equation */
/*        A'*Q*A - Q = -scale^2*C'*C, */
/*     in the Cholesky factor T, T'*T = Q, and use T to get the */
/*     transformation matrix. Copy A and C, to preserve them. */

/*     Workspace: need   N*(2*N + MAX(N,L) + 6) + MIN(N,L). */
/*                prefer larger. */

/*     Initialize the indices in the workspace. */

#line 283 "TB01VD.f"
    ldt = max(*n,*l);
#line 284 "TB01VD.f"
    ca = 1;
#line 285 "TB01VD.f"
    ia = 1;
#line 286 "TB01VD.f"
    it = ia + *n * *n;
#line 287 "TB01VD.f"
    iu = it + ldt * *n;
#line 288 "TB01VD.f"
    iwr = iu + *n * *n;
#line 289 "TB01VD.f"
    iwi = iwr + *n;

#line 291 "TB01VD.f"
    jwork = iwi + *n;

#line 293 "TB01VD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[ia], n, (ftnlen)4);
#line 294 "TB01VD.f"
    dlacpy_("Full", l, n, &c__[c_offset], ldc, &dwork[it], &ldt, (ftnlen)4);

#line 296 "TB01VD.f"
    i__1 = *ldwork - jwork + 1;
#line 296 "TB01VD.f"
    sb03od_("Discrete", "NotFactored", "NoTranspose", n, l, &dwork[ia], n, &
	    dwork[iu], n, &dwork[it], &ldt, &scale, &dwork[iwr], &dwork[iwi], 
	    &dwork[jwork], &i__1, info, (ftnlen)8, (ftnlen)11, (ftnlen)11);
#line 300 "TB01VD.f"
    if (*info != 0) {
#line 301 "TB01VD.f"
	if (*info == 6) {
#line 302 "TB01VD.f"
	    *info = 3;
#line 303 "TB01VD.f"
	} else {
#line 304 "TB01VD.f"
	    *info = 2;
#line 305 "TB01VD.f"
	}
#line 306 "TB01VD.f"
	return 0;
#line 307 "TB01VD.f"
    }
#line 308 "TB01VD.f"
    wrkopt = (integer) dwork[jwork] + jwork - 1;

#line 310 "TB01VD.f"
    if (scale == 0.) {
#line 311 "TB01VD.f"
	*info = 1;
#line 312 "TB01VD.f"
	return 0;
#line 313 "TB01VD.f"
    }

/*     Compute A = T*A*T^(-1). */

#line 317 "TB01VD.f"
    dtrmm_("Left", "Upper", "NoTranspose", "NonUnit", n, n, &c_b18, &dwork[it]
	    , &ldt, &a[a_offset], lda, (ftnlen)4, (ftnlen)5, (ftnlen)11, (
	    ftnlen)7);

#line 320 "TB01VD.f"
    dtrsm_("Right", "Upper", "NoTranspose", "NonUnit", n, n, &c_b18, &dwork[
	    it], &ldt, &a[a_offset], lda, (ftnlen)5, (ftnlen)5, (ftnlen)11, (
	    ftnlen)7);
#line 322 "TB01VD.f"
    if (*m > 0) {

/*        Compute B = (1/scale)*T*B. */

#line 326 "TB01VD.f"
	d__1 = 1. / scale;
#line 326 "TB01VD.f"
	dtrmm_("Left", "Upper", "NoTranspose", "NonUnit", n, m, &d__1, &dwork[
		it], &ldt, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (ftnlen)
		11, (ftnlen)7);
#line 328 "TB01VD.f"
    }

/*     Compute x0 = (1/scale)*T*x0. */

#line 332 "TB01VD.f"
    dtrmv_("Upper", "NoTranspose", "NonUnit", n, &dwork[it], &ldt, &x0[1], &
	    c__1, (ftnlen)5, (ftnlen)11, (ftnlen)7);
#line 334 "TB01VD.f"
    d__1 = 1. / scale;
#line 334 "TB01VD.f"
    dscal_(n, &d__1, &x0[1], &c__1);

/*     Compute C = scale*C*T^(-1). */

#line 338 "TB01VD.f"
    dtrsm_("Right", "Upper", "NoTranspose", "NonUnit", l, n, &scale, &dwork[
	    it], &ldt, &c__[c_offset], ldc, (ftnlen)5, (ftnlen)5, (ftnlen)11, 
	    (ftnlen)7);

/*     Now, the system has been transformed to the output normal form. */
/*     Build the transposed observability matrix in DWORK(CA) and compute */
/*     its QR factorization. */

#line 345 "TB01VD.f"
    ma02ad_("Full", l, n, &c__[c_offset], ldc, &dwork[ca], n, (ftnlen)4);

#line 347 "TB01VD.f"
    i__1 = *n - 1;
#line 347 "TB01VD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 348 "TB01VD.f"
	dgemm_("Transpose", "NoTranspose", n, l, n, &c_b18, &a[a_offset], lda,
		 &dwork[ca + (i__ - 1) * *n * *l], n, &c_b42, &dwork[ca + i__ 
		* *n * *l], n, (ftnlen)9, (ftnlen)11);
#line 350 "TB01VD.f"
/* L10: */
#line 350 "TB01VD.f"
    }

/*     Compute the QR factorization. */

/*     Workspace: need   N*N*L + N + L*N. */
/*                prefer N*N*L + N + NB*L*N. */

#line 357 "TB01VD.f"
    itau = ca + *n * *n * *l;
#line 358 "TB01VD.f"
    jwork = itau + *n;
#line 359 "TB01VD.f"
    i__1 = *l * *n;
#line 359 "TB01VD.f"
    i__2 = *ldwork - jwork + 1;
#line 359 "TB01VD.f"
    dgeqrf_(n, &i__1, &dwork[ca], n, &dwork[itau], &dwork[jwork], &i__2, info)
	    ;
/* Computing MAX */
#line 361 "TB01VD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 361 "TB01VD.f"
    wrkopt = max(i__1,i__2);

/*     Compute Q such that R has all diagonal elements nonnegative. */
/*     Only the first N*N part of R is needed. Move the details */
/*     of the QR factorization process, to gain memory and efficiency. */

/*     Workspace: need   2*N*N + 2*N. */
/*                prefer 2*N*N + N + NB*N. */

#line 370 "TB01VD.f"
    ir = *n * *n + 1;
#line 371 "TB01VD.f"
    if (*l != 2) {
#line 371 "TB01VD.f"
	dcopy_(n, &dwork[itau], &c__1, &dwork[ir + *n * *n], &c__1);
#line 371 "TB01VD.f"
    }
#line 373 "TB01VD.f"
    dlacpy_("Lower", n, n, &dwork[ca], n, &dwork[ir], n, (ftnlen)5);
#line 374 "TB01VD.f"
    itau = ir + *n * *n;
#line 375 "TB01VD.f"
    jwork = itau + *n;

#line 377 "TB01VD.f"
    iq = 1;
#line 378 "TB01VD.f"
    dlaset_("Full", n, n, &c_b42, &c_b18, &dwork[iq], n, (ftnlen)4);

#line 380 "TB01VD.f"
    i__1 = *n;
#line 380 "TB01VD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 381 "TB01VD.f"
	if (dwork[ir + (i__ - 1) * (*n + 1)] < 0.) {
#line 381 "TB01VD.f"
	    dwork[iq + (i__ - 1) * (*n + 1)] = -1.;
#line 381 "TB01VD.f"
	}
#line 383 "TB01VD.f"
/* L20: */
#line 383 "TB01VD.f"
    }

#line 385 "TB01VD.f"
    i__1 = *ldwork - jwork + 1;
#line 385 "TB01VD.f"
    dormqr_("Left", "NoTranspose", n, n, n, &dwork[ir], n, &dwork[itau], &
	    dwork[iq], n, &dwork[jwork], &i__1, info, (ftnlen)4, (ftnlen)11);
/* Computing MAX */
#line 388 "TB01VD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 388 "TB01VD.f"
    wrkopt = max(i__1,i__2);
#line 389 "TB01VD.f"
    jwork = ir;

/*     Now, the transformation matrix Q is in DWORK(IQ). */

/*     Compute A = Q'*A*Q. */

#line 395 "TB01VD.f"
    dgemm_("Transpose", "NoTranspose", n, n, n, &c_b18, &dwork[iq], n, &a[
	    a_offset], lda, &c_b42, &dwork[jwork], n, (ftnlen)9, (ftnlen)11);
#line 397 "TB01VD.f"
    dgemm_("NoTranspose", "NoTranspose", n, n, n, &c_b18, &dwork[jwork], n, &
	    dwork[iq], n, &c_b42, &a[a_offset], lda, (ftnlen)11, (ftnlen)11);

#line 400 "TB01VD.f"
    if (*m > 0) {

/*        Compute B = Q'*B. */
/*        Workspace: need   N*N + N*M. */

#line 405 "TB01VD.f"
	dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[jwork], n, (ftnlen)4);
#line 406 "TB01VD.f"
	dgemm_("Transpose", "NoTranspose", n, m, n, &c_b18, &dwork[iq], n, &
		dwork[jwork], n, &c_b42, &b[b_offset], ldb, (ftnlen)9, (
		ftnlen)11);
#line 408 "TB01VD.f"
    }

/*     Compute C = C*Q. */
/*     Workspace: need   N*N + N*L. */

#line 413 "TB01VD.f"
    dlacpy_("Full", l, n, &c__[c_offset], ldc, &dwork[jwork], l, (ftnlen)4);
#line 414 "TB01VD.f"
    dgemm_("NoTranspose", "NoTranspose", l, n, n, &c_b18, &dwork[jwork], l, &
	    dwork[iq], n, &c_b42, &c__[c_offset], ldc, (ftnlen)11, (ftnlen)11)
	    ;

/*     Compute x0 = Q'*x0. */

#line 419 "TB01VD.f"
    dcopy_(n, &x0[1], &c__1, &dwork[jwork], &c__1);
#line 420 "TB01VD.f"
    dgemv_("Transpose", n, n, &c_b18, &dwork[iq], n, &dwork[jwork], &c__1, &
	    c_b42, &x0[1], &c__1, (ftnlen)9);

/*     Now, copy C and A into the workspace to make it easier to read out */
/*     the corresponding part of THETA, and to apply the transformations. */

#line 426 "TB01VD.f"
    ldca = *n + *l;

#line 428 "TB01VD.f"
    i__1 = *n;
#line 428 "TB01VD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 429 "TB01VD.f"
	dcopy_(l, &c__[i__ * c_dim1 + 1], &c__1, &dwork[ca + (i__ - 1) * ldca]
		, &c__1);
#line 430 "TB01VD.f"
	dcopy_(n, &a[i__ * a_dim1 + 1], &c__1, &dwork[ca + *l + (i__ - 1) * 
		ldca], &c__1);
#line 431 "TB01VD.f"
/* L30: */
#line 431 "TB01VD.f"
    }

#line 433 "TB01VD.f"
    jwork = ca + ldca * *n;

/*     The parameters characterizing A and C are extracted in this loop. */
/*     Workspace: need   N*(N + L + 1). */

#line 438 "TB01VD.f"
    i__1 = *n;
#line 438 "TB01VD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 439 "TB01VD.f"
	dcopy_(l, &dwork[ca + 1 + (*n - i__) * (ldca + 1)], &c__1, &theta[(
		i__ - 1) * *l + 1], &c__1);
#line 441 "TB01VD.f"
	ri = dwork[ca + (*n - i__) * (ldca + 1)];
#line 442 "TB01VD.f"
	ti = dnrm2_(l, &theta[(i__ - 1) * *l + 1], &c__1);

/*        Multiply the part of [C; A] which will be currently transformed */
/*        with Ui = [ -THETAi, Si; RI, THETAi' ] from the left, without */
/*        storing Ui. Ui has the size (L+1)-by-(L+1). */

#line 448 "TB01VD.f"
	dgemv_("Transpose", l, n, &c_b18, &dwork[ca + *n - i__ + 1], &ldca, &
		theta[(i__ - 1) * *l + 1], &c__1, &c_b42, &dwork[jwork], &
		c__1, (ftnlen)9);

#line 451 "TB01VD.f"
	if (ti > 0.) {
#line 452 "TB01VD.f"
	    d__1 = (ri - 1.) / ti / ti;
#line 452 "TB01VD.f"
	    dger_(l, n, &d__1, &theta[(i__ - 1) * *l + 1], &c__1, &dwork[
		    jwork], &c__1, &dwork[ca + *n - i__ + 1], &ldca);
#line 454 "TB01VD.f"
	} else {

/*           The call below is for the limiting case. */

#line 458 "TB01VD.f"
	    dger_(l, n, &c_b93, &theta[(i__ - 1) * *l + 1], &c__1, &dwork[
		    jwork], &c__1, &dwork[ca + *n - i__ + 1], &ldca);
#line 460 "TB01VD.f"
	}

#line 462 "TB01VD.f"
	dger_(l, n, &c_b96, &theta[(i__ - 1) * *l + 1], &c__1, &dwork[ca + *n 
		- i__], &ldca, &dwork[ca + *n - i__ + 1], &ldca);
#line 464 "TB01VD.f"
	daxpy_(n, &ri, &dwork[ca + *n - i__], &ldca, &dwork[jwork], &c__1);

/*        Move these results to their appropriate locations. */

#line 468 "TB01VD.f"
	i__2 = *n;
#line 468 "TB01VD.f"
	for (j = 1; j <= i__2; ++j) {
#line 469 "TB01VD.f"
	    in = ca + *n - i__ + (j - 1) * ldca;
#line 470 "TB01VD.f"
	    i__3 = in + *l;
#line 470 "TB01VD.f"
	    for (k = in + 1; k <= i__3; ++k) {
#line 471 "TB01VD.f"
		dwork[k - 1] = dwork[k];
#line 472 "TB01VD.f"
/* L40: */
#line 472 "TB01VD.f"
	    }
#line 473 "TB01VD.f"
	    dwork[in + *l] = dwork[jwork + j - 1];
#line 474 "TB01VD.f"
/* L50: */
#line 474 "TB01VD.f"
	}

/*        Now, apply the bijective mapping, which allows to get rid */
/*        of the constraint norm(THETAi) < 1. */

#line 479 "TB01VD.f"
	if (lapply && ti != 0.) {
#line 479 "TB01VD.f"
	    d__1 = tan(ti * piby2) / ti;
#line 479 "TB01VD.f"
	    dscal_(l, &d__1, &theta[(i__ - 1) * *l + 1], &c__1);
#line 479 "TB01VD.f"
	}

#line 482 "TB01VD.f"
/* L60: */
#line 482 "TB01VD.f"
    }

#line 484 "TB01VD.f"
    if (*m > 0) {

/*        The next part of THETA is B. */

#line 488 "TB01VD.f"
	dlacpy_("Full", n, m, &b[b_offset], ldb, &theta[*n * *l + 1], n, (
		ftnlen)4);

/*        Copy the matrix D. */

#line 492 "TB01VD.f"
	dlacpy_("Full", l, m, &d__[d_offset], ldd, &theta[*n * (*l + *m) + 1],
		 l, (ftnlen)4);
#line 493 "TB01VD.f"
    }

/*     Copy the initial state x0. */

#line 497 "TB01VD.f"
    dcopy_(n, &x0[1], &c__1, &theta[*n * (*l + *m) + *l * *m + 1], &c__1);

#line 499 "TB01VD.f"
    dwork[1] = (doublereal) wrkopt;
#line 500 "TB01VD.f"
    return 0;

/* *** Last line of TB01VD *** */
} /* tb01vd_ */

